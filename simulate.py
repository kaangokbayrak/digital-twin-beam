"""
Simulation and visualization functions for beam vibration control.
"""

import numpy as np
from scipy.integrate import solve_ivp
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter


def simulate_free(A, x0, t_span, t_eval, C_out):
    """
    Simulate free vibration (no control).
    
    Parameters
    ----------
    A : ndarray
        State matrix
    x0 : ndarray
        Initial state
    t_span : tuple
        Time span (t_start, t_end)
    t_eval : ndarray
        Time points to evaluate
    C_out : ndarray
        Output matrix
    
    Returns
    -------
    t : ndarray
        Time array
    x : ndarray
        State history
    y : ndarray
        Output history
    """
    def dynamics(t, x):
        return A @ x
    
    sol = solve_ivp(dynamics, t_span, x0, t_eval=t_eval, method='RK45', rtol=1e-8)
    
    t = sol.t
    x = sol.y.T
    y = (C_out @ x.T).flatten()
    
    return t, x, y


def simulate_lqr(A, B, K_lqr, x0, t_span, t_eval, C_out):
    """
    Simulate LQR full-state feedback control.
    
    Parameters
    ----------
    A : ndarray
        State matrix
    B : ndarray
        Input matrix
    K_lqr : ndarray
        LQR gain matrix
    x0 : ndarray
        Initial state
    t_span : tuple
        Time span
    t_eval : ndarray
        Time points
    C_out : ndarray
        Output matrix
    
    Returns
    -------
    t : ndarray
        Time array
    x : ndarray
        State history
    y : ndarray
        Output history
    u : ndarray
        Control history
    """
    def dynamics(t, x):
        u = -K_lqr @ x
        return A @ x + B @ u
    
    sol = solve_ivp(dynamics, t_span, x0, t_eval=t_eval, method='RK45', rtol=1e-8)
    
    t = sol.t
    x = sol.y.T
    y = (C_out @ x.T).flatten()
    
    # Compute control for each time point
    u = np.array([(-K_lqr @ x[i]).item() for i in range(len(t))])
    
    return t, x, y, u


def simulate_lqg(A, B, C_out, lqg, x0, t_span, dt, noise_std=0.0):
    """
    Simulate LQG control with noisy measurements using fixed-step RK4 method.
    
    Parameters
    ----------
    A : ndarray
        State matrix
    B : ndarray
        Input matrix
    C_out : ndarray
        Output matrix
    lqg : LQGController
        LQG controller
    x0 : ndarray
        Initial state
    t_span : tuple
        Time span
    dt : float
        Fixed time step
    noise_std : float
        Measurement noise standard deviation
    
    Returns
    -------
    t : ndarray
        Time array
    x : ndarray
        State history
    x_hat : ndarray
        Estimated state history
    y : ndarray
        Output history
    u : ndarray
        Control history
    """
    # Reset LQG controller
    lqg.reset()
    
    # Time array
    t = np.arange(t_span[0], t_span[1], dt)
    n_steps = len(t)
    
    # Initialize histories
    n_states = x0.shape[0]
    x_hist = np.zeros((n_steps, n_states))
    x_hat_hist = np.zeros((n_steps, n_states))
    y_hist = np.zeros(n_steps)
    u_hist = np.zeros(n_steps)
    
    # Initial conditions
    x_hist[0] = x0
    x_hat_hist[0] = lqg.kalman.x_hat
    
    # Integration loop with RK4
    for i in range(n_steps - 1):
        x = x_hist[i]
        
        # Measure output with noise
        y_true = (C_out @ x).item()
        y_meas = y_true + noise_std * np.random.randn()
        y_hist[i] = y_true
        
        # LQG step: compute control and update observer
        u, x_hat = lqg.step(np.array([y_meas]), dt)
        u_val = u.item()
        u_hist[i] = u_val
        x_hat_hist[i] = x_hat
        
        # RK4 integration of plant (consistent with Kalman's continuous-time model)
        def f(x_):
            return A @ x_ + (B * u_val).flatten()
        
        k1 = f(x)
        k2 = f(x + 0.5 * dt * k1)
        k3 = f(x + 0.5 * dt * k2)
        k4 = f(x + dt * k3)
        x_hist[i+1] = x + (dt / 6.0) * (k1 + 2*k2 + 2*k3 + k4)
    
    # Final measurement
    y_hist[-1] = (C_out @ x_hist[-1]).item()
    u_hist[-1] = u_hist[-2]  # Hold last control
    x_hat_hist[-1] = x_hat_hist[-2]
    
    return t, x_hist, x_hat_hist, y_hist, u_hist


def plot_mode_shapes(fem, n_modes=3, filename='mode_shapes.png'):
    """
    Plot mode shapes with Hermite interpolation and analytical comparison.
    
    Parameters
    ----------
    fem : BeamFEM
        FEM model
    n_modes : int
        Number of modes to plot
    filename : str
        Output filename
    """
    freqs_fem, modes = fem.modal_analysis(n_modes)
    freqs_analytical = fem.analytical_frequencies(n_modes)
    
    fig, axes = plt.subplots(n_modes, 1, figsize=(10, 2.5*n_modes))
    if n_modes == 1:
        axes = [axes]
    
    for i in range(n_modes):
        ax = axes[i]
        
        # Get full mode shape
        mode_full = fem.get_full_mode_shape(modes[:, i])
        
        # Normalize mode shape
        mode_full = mode_full / np.max(np.abs(mode_full))
        
        # Interpolate for smooth curve
        x_interp, w_interp = fem.interpolate_mode(mode_full, n_interp=200)
        
        # Plot smooth mode shape
        ax.plot(x_interp, w_interp, 'b-', linewidth=2, label=f'FEM: {freqs_fem[i]:.2f} Hz')
        ax.fill_between(x_interp, 0, w_interp, alpha=0.2, color='blue')
        
        # Plot node markers
        x_nodes = np.linspace(0, fem.L, fem.n_nodes)
        w_nodes = mode_full[::2]  # Extract translational DOFs
        ax.plot(x_nodes, w_nodes, 'ko', markersize=4, label='Nodes')
        
        # Error calculation
        error_pct = 100 * (freqs_fem[i] - freqs_analytical[i]) / freqs_analytical[i]
        
        # Labels
        ax.set_xlabel('Position (m)')
        ax.set_ylabel('Normalized Amplitude')
        ax.set_title(f'Mode {i+1} — FEM: {freqs_fem[i]:.3f} Hz, '
                    f'Analytical: {freqs_analytical[i]:.3f} Hz, '
                    f'Error: {error_pct:.2f}%')
        ax.grid(True, alpha=0.3)
        ax.legend(loc='upper right')
        ax.axhline(0, color='k', linewidth=0.5)
        
        # Mark fixed end
        ax.axvline(0, color='red', linestyle='--', linewidth=1.5, alpha=0.7, label='Fixed')
    
    plt.tight_layout()
    plt.savefig(filename, dpi=200, bbox_inches='tight')
    plt.close()
    print(f"Saved mode shapes to {filename}")


def plot_transient_comparison(t_free, y_free, t_lqr, y_lqr, u_lqr, 
                               t_lqg, y_lqg, u_lqg,
                               filename='transient_response.png'):
    """
    Plot comparison of free, LQR, and LQG responses.
    
    Parameters
    ----------
    t_free, y_free : ndarray
        Free vibration time and output
    t_lqr, y_lqr, u_lqr : ndarray
        LQR time, output, and control
    t_lqg, y_lqg, u_lqg : ndarray
        LQG time, output, and control
    filename : str
        Output filename
    """
    fig, axes = plt.subplots(3, 1, figsize=(12, 10))
    
    # Subplot 1: Tip displacement comparison
    ax = axes[0]
    ax.plot(t_free, y_free*1000, 'k-', linewidth=2, label='Free vibration', alpha=0.7)
    ax.plot(t_lqr, y_lqr*1000, 'b-', linewidth=2, label='LQR control')
    ax.plot(t_lqg, y_lqg*1000, 'r--', linewidth=2, label='LQG control')
    ax.set_xlabel('Time (s)')
    ax.set_ylabel('Tip Displacement (mm)')
    ax.set_title('Transient Response Comparison')
    ax.grid(True, alpha=0.3)
    ax.legend()
    ax.axhline(0, color='gray', linewidth=0.5)
    
    # Subplot 2: Control effort
    ax = axes[1]
    ax.plot(t_lqr, u_lqr, 'b-', linewidth=2, label='LQR control')
    ax.plot(t_lqg, u_lqg, 'r--', linewidth=2, label='LQG control')
    ax.set_xlabel('Time (s)')
    ax.set_ylabel('Control Force (N)')
    ax.set_title('Control Effort')
    ax.grid(True, alpha=0.3)
    ax.legend()
    ax.axhline(0, color='gray', linewidth=0.5)
    
    # Subplot 3: Energy decay (log scale)
    ax = axes[2]
    # Energy proportional to displacement squared
    # Add small floor to prevent log(0) or log(negative) errors
    eps = 1e-30
    energy_free = np.maximum(y_free**2, eps)
    energy_lqr = np.maximum(y_lqr**2, eps)
    energy_lqg = np.maximum(y_lqg**2, eps)
    
    ax.semilogy(t_free, energy_free / energy_free[0], 'k-', linewidth=2, 
                label='Free vibration', alpha=0.7)
    ax.semilogy(t_lqr, energy_lqr / energy_lqr[0], 'b-', linewidth=2, 
                label='LQR control')
    ax.semilogy(t_lqg, energy_lqg / energy_lqg[0], 'r--', linewidth=2, 
                label='LQG control')
    ax.set_xlabel('Time (s)')
    ax.set_ylabel('Normalized Energy (log scale)')
    ax.set_title('Energy Decay')
    ax.grid(True, alpha=0.3, which='both')
    ax.legend()
    
    plt.tight_layout()
    plt.savefig(filename, dpi=200, bbox_inches='tight')
    plt.close()
    print(f"Saved transient comparison to {filename}")


def animate_beam(fem, t_arr, x_hist, filename='beam_animation.gif', 
                 title='Beam Animation', fps=30, duration=None):
    """
    Create animated GIF of beam vibration with Hermite interpolation.
    
    Parameters
    ----------
    fem : BeamFEM
        FEM model
    t_arr : ndarray
        Time array
    x_hist : ndarray
        State history (n_steps x n_states)
    filename : str
        Output filename
    title : str
        Animation title
    fps : int
        Frames per second
    duration : float, optional
        Duration in seconds (if None, uses full t_arr)
    """
    # Downsample if needed for animation
    if duration is not None:
        # Determine frame indices for desired duration
        n_frames = int(duration * fps)
        frame_indices = np.linspace(0, len(t_arr)-1, n_frames, dtype=int)
    else:
        frame_indices = np.arange(0, len(t_arr), max(1, len(t_arr) // 200))
    
    # Setup figure
    fig, ax = plt.subplots(figsize=(12, 6))
    
    # Determine y-axis limits
    y_max = 0
    for idx in frame_indices[::5]:  # Sample every 5th frame for speed
        q = x_hist[idx, :fem.n_dof_free]
        mode_full = fem.get_full_mode_shape(q)
        _, w_interp = fem.interpolate_mode(mode_full, n_interp=100)
        y_max = max(y_max, np.max(np.abs(w_interp)))
    
    y_lim = y_max * 1.2 * 1000  # Convert to mm with margin
    
    # Animation function
    def update(frame_idx):
        ax.clear()
        
        t = t_arr[frame_idx]
        q = x_hist[frame_idx, :fem.n_dof_free]
        mode_full = fem.get_full_mode_shape(q)
        x_interp, w_interp = fem.interpolate_mode(mode_full, n_interp=100)
        
        # Plot beam
        ax.plot(x_interp, w_interp*1000, 'b-', linewidth=3)
        ax.fill_between(x_interp, 0, w_interp*1000, alpha=0.3, color='blue')
        
        # Plot nodes
        x_nodes = np.linspace(0, fem.L, fem.n_nodes)
        w_nodes = mode_full[::2] * 1000
        ax.plot(x_nodes, w_nodes, 'ko', markersize=5)
        
        # Fixed end marker
        ax.axvline(0, color='red', linewidth=4, alpha=0.8)
        ax.fill_betweenx([-y_lim, y_lim], -0.02, 0, color='gray', alpha=0.3, 
                         hatch='///', label='Fixed end')
        
        # Formatting
        ax.set_xlim(-0.05, fem.L + 0.05)
        ax.set_ylim(-y_lim, y_lim)
        ax.set_xlabel('Position (m)', fontsize=12)
        ax.set_ylabel('Displacement (mm)', fontsize=12)
        ax.set_title(f'{title} — Time: {t:.3f} s', fontsize=14)
        ax.grid(True, alpha=0.3)
        ax.axhline(0, color='k', linewidth=0.5)
        
    # Create animation
    anim = FuncAnimation(fig, update, frames=frame_indices, interval=1000/fps)
    
    # Save as GIF
    writer = PillowWriter(fps=fps)
    anim.save(filename, writer=writer)
    plt.close()
    
    print(f"Saved animation to {filename}")
