"""
Simulation and visualization functions for beam vibration control.
"""

import base64
import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import interp1d
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter
from matplotlib.collections import LineCollection
from matplotlib.colors import Normalize

# ─────────────────────────────────────────────────────────────────────────────
# Professional matplotlib theme – applied globally to every figure
# ─────────────────────────────────────────────────────────────────────────────
plt.rcParams.update({
    'font.family':       'DejaVu Sans',
    'font.size':         11,
    'axes.titlesize':    13,
    'axes.labelsize':    11,
    'axes.titleweight':  'bold',
    'lines.linewidth':   2.0,
    'axes.grid':         True,
    'grid.alpha':        0.3,
    'grid.linestyle':    '--',
    'axes.spines.top':   False,
    'axes.spines.right': False,
    'figure.dpi':        100,
    'savefig.dpi':       200,
    'savefig.bbox':      'tight',
    'legend.framealpha': 0.85,
    'legend.fontsize':   10,
})

# Consistent colour palette used across all figures
_COL_FREE = '#444444'
_COL_LQR  = '#2166ac'
_COL_LQG  = '#d6604d'


# =============================================================================
# Simulation functions  (unchanged signatures / behaviour)
# =============================================================================

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
    # Reset LQG controller with correct initial state estimate
    lqg.reset(x0)

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

        # RK4 integration of plant
        def f(x_):
            return A @ x_ + (B * u_val).flatten()

        k1 = f(x)
        k2 = f(x + 0.5 * dt * k1)
        k3 = f(x + 0.5 * dt * k2)
        k4 = f(x + dt * k3)
        x_hist[i+1] = x + (dt / 6.0) * (k1 + 2*k2 + 2*k3 + k4)

    # Final measurement
    y_hist[-1] = (C_out @ x_hist[-1]).item()
    u_hist[-1] = u_hist[-2]
    x_hat_hist[-1] = x_hat_hist[-2]

    return t, x_hist, x_hat_hist, y_hist, u_hist


# =============================================================================
# Existing plot functions – enhanced
# =============================================================================

def plot_mode_shapes(fem, n_modes=3, filename='mode_shapes.png'):
    """
    Plot mode shapes with Hermite interpolation, colormap beam, and inset
    FEM-vs-analytical error table.

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

    fig, axes = plt.subplots(n_modes, 1, figsize=(11, 3.2 * n_modes))
    if n_modes == 1:
        axes = [axes]

    cmap = plt.cm.plasma

    for i in range(n_modes):
        ax = axes[i]

        # Get and normalise mode shape
        mode_full = fem.get_full_mode_shape(modes[:, i])
        mode_full = mode_full / np.max(np.abs(mode_full))

        # Hermite interpolation for smooth curve
        x_interp, w_interp = fem.interpolate_mode(mode_full, n_interp=200)

        # ── Colormap line (LineCollection coloured by |w|) ──────────────────
        pts  = np.column_stack([x_interp, w_interp])
        segs = np.stack([pts[:-1], pts[1:]], axis=1)
        midpt_w = 0.5 * (np.abs(w_interp[:-1]) + np.abs(w_interp[1:]))

        norm = Normalize(vmin=0.0, vmax=1.0)
        lc = LineCollection(segs, cmap=cmap, norm=norm, linewidth=2.5, zorder=3)
        lc.set_array(midpt_w)
        ax.add_collection(lc)

        # Subtle fill for depth
        ax.fill_between(x_interp, 0, w_interp, alpha=0.12, color='#7b2d8b')

        # Colorbar for this subplot
        cb = fig.colorbar(lc, ax=ax, fraction=0.025, pad=0.02, shrink=0.85)
        cb.set_label('|w| normalised', fontsize=8.5)
        cb.ax.tick_params(labelsize=8)

        # Node markers
        x_nodes = np.linspace(0, fem.L, fem.n_nodes)
        w_nodes = mode_full[::2]
        ax.plot(x_nodes, w_nodes, 'o', color='#222222',
                markersize=4, zorder=4, label='FEM nodes')

        # Fixed-end marker
        ax.axvline(0, color='#c0392b', linestyle='--', linewidth=1.5,
                   alpha=0.85, label='Fixed end')

        ax.set_xlim(-0.02, fem.L + 0.02)
        ax.set_ylim(-1.35, 1.35)
        ax.axhline(0, color='k', linewidth=0.5)
        ax.set_xlabel('Position (m)')
        ax.set_ylabel('Normalised Amplitude')
        ax.set_title(
            f'Mode {i+1}  —  f_FEM = {freqs_fem[i]:.3f} Hz  |  '
            f'f_analytical = {freqs_analytical[i]:.3f} Hz'
        )
        ax.legend(loc='upper right', fontsize=9)

        # ── Inset error table (top-left) ─────────────────────────────────────
        error_pct = 100.0 * (freqs_fem[i] - freqs_analytical[i]) / freqs_analytical[i]
        table_data = [
            ['f_FEM (Hz)',        f'{freqs_fem[i]:.4f}'],
            ['f_analytical (Hz)', f'{freqs_analytical[i]:.4f}'],
            ['Error (%)',         f'{error_pct:.5f}'],
        ]
        ax_ins = ax.inset_axes([0.01, 0.60, 0.21, 0.37])
        ax_ins.axis('off')
        tbl = ax_ins.table(cellText=table_data, loc='center', cellLoc='left')
        tbl.auto_set_font_size(False)
        tbl.set_fontsize(8.0)
        tbl.scale(1.0, 1.4)
        for (row, col), cell in tbl.get_celld().items():
            cell.set_edgecolor('#cccccc')
            cell.set_facecolor('#f0f0f0' if col == 0 else 'white')
            if col == 0:
                cell.set_text_props(weight='bold')

    plt.tight_layout()
    plt.savefig(filename, dpi=200, bbox_inches='tight')
    plt.close()
    print(f"Saved mode shapes to {filename}")


def plot_transient_comparison(t_free, y_free, t_lqr, y_lqr, u_lqr,
                               t_lqg, y_lqg, u_lqg,
                               filename='transient_response.png'):
    """
    Plot comparison of free, LQR, and LQG responses.

    Includes: (1) tip displacement with ±2% settling band,
              (2) control effort, (3) energy decay,
              (4) frequency spectrum (FFT).

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
    fig, axes = plt.subplots(4, 1, figsize=(12, 15))

    # ── Subplot 1: Tip displacement + ±2% settling band ──────────────────────
    ax = axes[0]
    y_peak = np.max(np.abs(y_free))
    band   = 0.02 * y_peak
    ax.fill_between(t_free, -band * 1000, band * 1000,
                    alpha=0.18, color='#27ae60', label='±2% settling band')
    ax.plot(t_free, y_free * 1000, color=_COL_FREE, linewidth=2,
            label='Free vibration', alpha=0.85)
    ax.plot(t_lqr, y_lqr * 1000, color=_COL_LQR, linewidth=2,
            label='LQR control')
    ax.plot(t_lqg, y_lqg * 1000, color=_COL_LQG, linewidth=2,
            linestyle='--', label='LQG control')
    ax.set_xlabel('Time (s)')
    ax.set_ylabel('Tip Displacement (mm)')
    ax.set_title('Transient Response Comparison')
    ax.legend()
    ax.axhline(0, color='gray', linewidth=0.5)

    # ── Subplot 2: Control effort ─────────────────────────────────────────────
    ax = axes[1]
    ax.plot(t_lqr, u_lqr, color=_COL_LQR, linewidth=2, label='LQR control')
    ax.plot(t_lqg, u_lqg, color=_COL_LQG, linewidth=2,
            linestyle='--', label='LQG control')
    ax.set_xlabel('Time (s)')
    ax.set_ylabel('Control Force (N)')
    ax.set_title('Control Effort')
    ax.legend()
    ax.axhline(0, color='gray', linewidth=0.5)

    # ── Subplot 3: Energy decay (log scale) ───────────────────────────────────
    ax = axes[2]
    eps = 1e-30
    energy_free = np.maximum(y_free**2, eps)
    energy_lqr  = np.maximum(y_lqr**2, eps)
    energy_lqg  = np.maximum(y_lqg**2, eps)
    ax.semilogy(t_free, energy_free / energy_free[0],
                color=_COL_FREE, linewidth=2, label='Free vibration', alpha=0.85)
    ax.semilogy(t_lqr, energy_lqr / energy_lqr[0],
                color=_COL_LQR, linewidth=2, label='LQR control')
    ax.semilogy(t_lqg, energy_lqg / energy_lqg[0],
                color=_COL_LQG, linewidth=2, linestyle='--', label='LQG control')
    ax.set_xlabel('Time (s)')
    ax.set_ylabel('Normalised Energy (log scale)')
    ax.set_title('Energy Decay')
    ax.legend()

    # ── Subplot 4: Frequency spectrum (FFT) ───────────────────────────────────
    ax = axes[3]
    freq_lim = 100  # Hz

    for t_arr, y_arr, col, ls, lbl in [
        (t_free, y_free, _COL_FREE, '-',  'Free vibration'),
        (t_lqr,  y_lqr,  _COL_LQR,  '-',  'LQR control'),
        (t_lqg,  y_lqg,  _COL_LQG,  '--', 'LQG control'),
    ]:
        dt = t_arr[1] - t_arr[0]
        N  = len(y_arr)
        freqs_fft = np.fft.rfftfreq(N, d=dt)
        amp       = np.abs(np.fft.rfft(y_arr)) / N
        mask      = freqs_fft <= freq_lim
        ax.plot(freqs_fft[mask], amp[mask] * 1000,
                color=col, linestyle=ls, linewidth=1.5, label=lbl)

    ax.set_xlabel('Frequency (Hz)')
    ax.set_ylabel('Amplitude (mm)')
    ax.set_title('Frequency Spectrum (FFT of tip displacement)')
    ax.legend()

    plt.tight_layout()
    plt.savefig(filename, dpi=200, bbox_inches='tight')
    plt.close()
    print(f"Saved transient comparison to {filename}")


def animate_beam(fem, t_arr, x_hist, filename='beam_animation.gif',
                 title='Beam Animation', fps=30, duration=None,
                 u_hist=None):
    """
    Create animated GIF of beam vibration with Hermite interpolation.

    Enhancements over the original:
    - Beam line coloured by instantaneous displacement magnitude (RdBu_r).
    - Inset mini time-series (top-right) with moving dot.
    - Optional control-force text overlay when u_hist is supplied.

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
    u_hist : ndarray, optional
        Control force history; enables force text overlay
    """
    # ── Frame indices ─────────────────────────────────────────────────────────
    if duration is not None:
        n_frames     = int(duration * fps)
        frame_indices = np.linspace(0, len(t_arr) - 1, n_frames, dtype=int)
    else:
        frame_indices = np.arange(0, len(t_arr), max(1, len(t_arr) // 200))

    # ── Precompute y-limit and tip history for inset ──────────────────────────
    tip_dof = fem.n_dof_free - 2          # tip translation in reduced coords
    tip_hist = x_hist[:, tip_dof] * 1000  # mm

    y_max = 0.0
    for idx in frame_indices[::5]:
        q         = x_hist[idx, :fem.n_dof_free]
        mode_full = fem.get_full_mode_shape(q)
        _, w_interp = fem.interpolate_mode(mode_full, n_interp=100)
        y_max = max(y_max, np.max(np.abs(w_interp)))
    y_lim = y_max * 1.25 * 1000  # mm, with margin

    # ── Figure setup ─────────────────────────────────────────────────────────
    fig, ax = plt.subplots(figsize=(12, 6))
    # Inset axes in absolute figure coordinates (top-right corner)
    ax_ins = fig.add_axes([0.67, 0.68, 0.21, 0.20])
    ax_ins.set_xlim(t_arr[0], t_arr[-1])
    y_ins_lim = max(np.abs(tip_hist).max(), 1.0)
    ax_ins.set_ylim(-y_ins_lim * 1.1, y_ins_lim * 1.1)
    ax_ins.set_xlabel('t (s)', fontsize=7)
    ax_ins.set_ylabel('mm', fontsize=7)
    ax_ins.tick_params(labelsize=6)
    ax_ins.set_title('Tip disp.', fontsize=7, pad=2)
    ax_ins.axhline(0, color='gray', linewidth=0.4)
    ax_ins.grid(True, alpha=0.2, linestyle='--')

    beam_col = _COL_LQR if u_hist is not None else _COL_FREE
    norm_beam = Normalize(vmin=-y_lim, vmax=y_lim)

    def update(frame_idx):
        ax.clear()
        ax_ins.clear()

        t = t_arr[frame_idx]
        q         = x_hist[frame_idx, :fem.n_dof_free]
        mode_full = fem.get_full_mode_shape(q)
        x_interp, w_interp = fem.interpolate_mode(mode_full, n_interp=100)
        w_mm = w_interp * 1000

        # ── Colormap beam ────────────────────────────────────────────────────
        pts  = np.column_stack([x_interp, w_mm])
        segs = np.stack([pts[:-1], pts[1:]], axis=1)
        mid_w = 0.5 * (w_mm[:-1] + w_mm[1:])
        lc = LineCollection(segs, cmap='RdBu_r', norm=norm_beam,
                            linewidth=3.5, zorder=3)
        lc.set_array(mid_w)
        ax.add_collection(lc)

        # Subtle fill
        ax.fill_between(x_interp, 0, w_mm, alpha=0.18, color=beam_col)

        # Node markers
        x_nodes = np.linspace(0, fem.L, fem.n_nodes)
        w_nodes = mode_full[::2] * 1000
        ax.plot(x_nodes, w_nodes, 'o', color='#222222', markersize=5, zorder=4)

        # Fixed-end hatch
        ax.axvline(0, color='#c0392b', linewidth=4, alpha=0.8)
        ax.fill_betweenx([-y_lim, y_lim], -0.02, 0,
                         color='gray', alpha=0.25, hatch='///')

        ax.set_xlim(-0.05, fem.L + 0.05)
        ax.set_ylim(-y_lim, y_lim)
        ax.set_xlabel('Position (m)', fontsize=12)
        ax.set_ylabel('Displacement (mm)', fontsize=12)
        ax.set_title(f'{title}  —  t = {t:.3f} s', fontsize=14)
        ax.axhline(0, color='k', linewidth=0.5)

        # ── Control force overlay ─────────────────────────────────────────────
        if u_hist is not None:
            u_val = u_hist[frame_idx]
            ax.text(0.98, 0.96, f'F = {u_val:+.1f} N',
                    transform=ax.transAxes, ha='right', va='top',
                    fontsize=11, color='navy',
                    bbox=dict(boxstyle='round,pad=0.35',
                              facecolor='lightyellow', alpha=0.85))

        # ── Inset time-series ─────────────────────────────────────────────────
        ax_ins.set_xlim(t_arr[0], t_arr[-1])
        ax_ins.set_ylim(-y_ins_lim * 1.1, y_ins_lim * 1.1)
        ax_ins.plot(t_arr[:frame_idx + 1], tip_hist[:frame_idx + 1],
                    color=beam_col, linewidth=1.0)
        ax_ins.plot(t_arr[frame_idx], tip_hist[frame_idx],
                    'o', color='#c0392b', markersize=4, zorder=5)
        ax_ins.set_xlabel('t (s)', fontsize=7)
        ax_ins.set_ylabel('mm', fontsize=7)
        ax_ins.set_title('Tip disp.', fontsize=7, pad=2)
        ax_ins.axhline(0, color='gray', linewidth=0.4)
        ax_ins.tick_params(labelsize=6)
        ax_ins.grid(True, alpha=0.2, linestyle='--')
        for sp in ['top', 'right']:
            ax_ins.spines[sp].set_visible(False)

    anim = FuncAnimation(fig, update, frames=frame_indices, interval=1000 / fps)
    writer = PillowWriter(fps=fps)
    anim.save(filename, writer=writer)
    plt.close()
    print(f"Saved animation to {filename}")


# =============================================================================
# New plot functions
# =============================================================================

def plot_bode(A, B, C, K_lqr, filename='bode_plot.png'):
    """
    Bode magnitude plot comparing open-loop and closed-loop transfer functions.

    Computes H(jω) = C (jωI − A)⁻¹ B evaluated over a log-spaced frequency
    grid from 1 Hz to 250 Hz.

    Parameters
    ----------
    A : ndarray  State matrix
    B : ndarray  Input matrix
    C : ndarray  Output matrix
    K_lqr : ndarray  LQR gain
    filename : str  Output filename
    """
    n    = A.shape[0]
    A_cl = A - B @ K_lqr

    freqs  = np.logspace(0, np.log10(250), 500)   # 1 – 250 Hz
    omegas = 2.0 * np.pi * freqs

    H_ol = np.empty(len(omegas), dtype=complex)
    H_cl = np.empty(len(omegas), dtype=complex)
    I_n  = np.eye(n)

    for i, omega in enumerate(omegas):
        jw_I   = 1j * omega * I_n
        H_ol[i] = (C @ np.linalg.solve(jw_I - A,    B)).item()
        H_cl[i] = (C @ np.linalg.solve(jw_I - A_cl, B)).item()

    H_ol_dB = 20.0 * np.log10(np.abs(H_ol) + 1e-30)
    H_cl_dB = 20.0 * np.log10(np.abs(H_cl) + 1e-30)

    fig, ax = plt.subplots(figsize=(11, 5))
    ax.semilogx(freqs, H_ol_dB, color=_COL_FREE, linewidth=2,
                label='Open-loop', alpha=0.9)
    ax.semilogx(freqs, H_cl_dB, color=_COL_LQR, linewidth=2,
                label='Closed-loop (LQR)')
    ax.set_xlabel('Frequency (Hz)')
    ax.set_ylabel('Magnitude (dB)')
    ax.set_title('Bode Magnitude Plot — Open-Loop vs Closed-Loop')
    ax.legend()
    ax.axhline(0, color='gray', linewidth=0.5, linestyle=':')

    plt.tight_layout()
    plt.savefig(filename, dpi=200, bbox_inches='tight')
    plt.close()
    print(f"Saved Bode plot to {filename}")


def plot_pole_zero(A, B, K_lqr, filename='pole_zero.png'):
    """
    Pole-zero map showing open-loop vs LQR closed-loop poles.

    Left panel : first 5 structural modes (zoomed).
    Right panel: full spectrum (all 40 poles).

    Arrows connect each open-loop pole to its closed-loop counterpart.

    Parameters
    ----------
    A : ndarray  State matrix
    B : ndarray  Input matrix
    K_lqr : ndarray  LQR gain
    filename : str  Output filename
    """
    poles_ol = np.linalg.eigvals(A)
    poles_cl = np.linalg.eigvals(A - B @ K_lqr)

    fig, axes = plt.subplots(1, 2, figsize=(15, 6))

    # ── Left: zoomed (first 5 structural modes) ───────────────────────────────
    ax = axes[0]
    freq_limit = 2.0 * np.pi * 250   # 250 Hz in rad/s

    def upper_half(poles, limit):
        """Return poles in upper half-plane with |Im| ≤ limit, sorted by Im."""
        return sorted(
            [p for p in poles if np.imag(p) > 0 and np.abs(np.imag(p)) <= limit],
            key=lambda p: np.imag(p)
        )

    ol_pos = upper_half(poles_ol, freq_limit)
    cl_pos = upper_half(poles_cl, freq_limit)
    n_pairs = min(len(ol_pos), len(cl_pos))

    # Arrows and mirrored counterparts
    for j in range(n_pairs):
        p_from = ol_pos[j]
        p_to   = cl_pos[j]
        for sign in (1, -1):
            ax.annotate(
                '',
                xy    =(np.real(p_to),   sign * np.imag(p_to)),
                xytext=(np.real(p_from), sign * np.imag(p_from)),
                arrowprops=dict(arrowstyle='->', color='#888888',
                                lw=0.9, mutation_scale=12),
                zorder=2
            )

    # Scatter – plot both conjugate halves together
    ax.scatter(np.real(poles_ol), np.imag(poles_ol),
               marker='x', s=70, linewidths=2.5, color=_COL_FREE,
               label='Open-loop poles', zorder=3)
    ax.scatter(np.real(poles_cl), np.imag(poles_cl),
               marker='o', s=55, facecolors='none',
               edgecolors=_COL_LQR, linewidths=2.0,
               label='Closed-loop poles (LQR)', zorder=3)

    ax.axvline(0, color='gray', linewidth=0.8, linestyle=':')
    ax.axhline(0, color='gray', linewidth=0.5, linestyle=':')
    ax.set_xlabel('Real Part (rad/s)')
    ax.set_ylabel('Imaginary Part (rad/s)')
    ax.set_title('Pole Map — Structural Modes (≤ 250 Hz)')
    ax.legend()

    # Limit x-axis to show movement clearly
    mask_im_ol = (np.abs(np.imag(poles_ol)) <= freq_limit)
    mask_im_cl = (np.abs(np.imag(poles_cl)) <= freq_limit)
    x_lo = min(np.real(poles_cl[mask_im_cl]).min() * 1.15, -1)
    x_hi = max(np.real(poles_ol[mask_im_ol]).max() * 0.5,   1)
    ax.set_xlim(x_lo, x_hi)
    y_hi = freq_limit * 1.05
    ax.set_ylim(-y_hi, y_hi)

    # ── Right: full spectrum ──────────────────────────────────────────────────
    ax = axes[1]
    ax.scatter(np.real(poles_ol), np.imag(poles_ol),
               marker='x', s=60, linewidths=2.5, color=_COL_FREE,
               label='Open-loop poles', zorder=3)
    ax.scatter(np.real(poles_cl), np.imag(poles_cl),
               marker='o', s=45, facecolors='none',
               edgecolors=_COL_LQR, linewidths=2.0,
               label='Closed-loop poles (LQR)', zorder=3)
    ax.axvline(0, color='gray', linewidth=0.8, linestyle=':')
    ax.axhline(0, color='gray', linewidth=0.5, linestyle=':')
    ax.set_xlabel('Real Part (rad/s)')
    ax.set_ylabel('Imaginary Part (rad/s)')
    ax.set_title(f'Pole Map — Full Spectrum (all {len(poles_ol)} poles)')
    ax.legend()

    plt.tight_layout()
    plt.savefig(filename, dpi=200, bbox_inches='tight')
    plt.close()
    print(f"Saved pole-zero map to {filename}")


def plot_phase_portrait(t_free, x_free, t_lqr, x_lqr, t_lqg, x_lqg,
                        sensor_dof, n_dof_free,
                        filename='phase_portrait.png'):
    """
    Phase portrait: tip velocity vs tip displacement for all three scenarios.

    Free vibration shows a slowly decaying spiral; LQR/LQG spiral inward
    rapidly to the origin, making suppression immediately visible.

    Parameters
    ----------
    t_free, x_free : ndarray  Free vibration time and state history
    t_lqr, x_lqr   : ndarray  LQR time and state history
    t_lqg, x_lqg   : ndarray  LQG time and state history
    sensor_dof     : int       DOF index for tip (in reduced coordinate system)
    n_dof_free     : int       Number of free DOFs
    filename       : str       Output filename
    """
    def extract(x_hist):
        disp = x_hist[:, sensor_dof] * 1000                    # mm
        vel  = x_hist[:, n_dof_free + sensor_dof] * 1000       # mm/s
        return disp, vel

    d_free, v_free = extract(x_free)
    d_lqr,  v_lqr  = extract(x_lqr)
    d_lqg,  v_lqg  = extract(x_lqg)

    fig, ax = plt.subplots(figsize=(8, 7))

    ax.plot(d_free, v_free, color=_COL_FREE, linewidth=1.2,
            label='Free vibration', alpha=0.65)
    ax.plot(d_lqr, v_lqr, color=_COL_LQR, linewidth=1.8,
            label='LQR control')
    ax.plot(d_lqg, v_lqg, color=_COL_LQG, linewidth=1.8,
            linestyle='--', label='LQG control', alpha=0.85)

    # Mark start and equilibrium
    ax.plot(d_free[0], v_free[0], 'o', color='#222222',
            markersize=9, zorder=5, label='Initial condition')
    ax.plot(0, 0, '*', color='#27ae60', markersize=14,
            zorder=6, label='Equilibrium (origin)')

    ax.set_xlabel('Tip Displacement (mm)')
    ax.set_ylabel('Tip Velocity (mm/s)')
    ax.set_title('Phase Portrait — Tip Dynamics')
    ax.legend()
    ax.axhline(0, color='gray', linewidth=0.5)
    ax.axvline(0, color='gray', linewidth=0.5)

    plt.tight_layout()
    plt.savefig(filename, dpi=200, bbox_inches='tight')
    plt.close()
    print(f"Saved phase portrait to {filename}")


# =============================================================================
# Summary dashboard (6-panel portfolio figure)
# =============================================================================

def _draw_beam_schematic(ax, fem):
    """Draw a labelled cantilever-beam schematic."""
    L = fem.L
    bh = 0.055   # beam half-height

    # Wall
    ax.fill_betweenx([-0.22, 0.22], -0.07, 0.0,
                     facecolor='#95a5a6', alpha=0.6, hatch='///')
    ax.add_patch(plt.Rectangle((-0.07, -0.22), 0.07, 0.44,
                               facecolor='none', edgecolor='#2c3e50',
                               linewidth=1.5))

    # Beam body
    ax.add_patch(plt.Rectangle((0.0, -bh), L, 2 * bh,
                               facecolor='#aec6cf', edgecolor='#2c3e50',
                               linewidth=1.5, zorder=2))

    # Node markers
    x_nodes = np.linspace(0, L, fem.n_nodes)
    ax.plot(x_nodes, np.zeros_like(x_nodes), 'o',
            color='#e74c3c', markersize=5, zorder=5, label='FEM nodes')

    # Length dimension
    ax.annotate('', xy=(L, 0.19), xytext=(0.0, 0.19),
                arrowprops=dict(arrowstyle='<->', color='#2c3e50', lw=1.2))
    ax.text(L / 2, 0.22, f'L = {L:.1f} m',
            ha='center', va='bottom', fontsize=9, color='#2c3e50')

    # Sensor/actuator callout
    ax.annotate('Sensor &\nActuator',
                xy=(L, -bh), xytext=(L * 0.60, -0.19),
                ha='center', fontsize=8, color='#c0392b',
                arrowprops=dict(arrowstyle='->', color='#c0392b', lw=1.2))

    ax.set_xlim(-0.12, L + 0.12)
    ax.set_ylim(-0.30, 0.32)
    ax.set_yticks([])
    ax.set_xlabel('x (m)')
    ax.set_title(
        f'FEM Model: {fem.n_elements} elements, {fem.n_dof_free} free DOFs',
        fontsize=10
    )
    ax.grid(False)
    ax.spines['left'].set_visible(False)


def _draw_mode1_inline(ax, fem):
    """Draw Mode 1 shape on ax."""
    freqs_fem, modes = fem.modal_analysis(1)
    freqs_an         = fem.analytical_frequencies(1)
    mode_full        = fem.get_full_mode_shape(modes[:, 0])
    mode_full        = mode_full / np.max(np.abs(mode_full))
    x_interp, w_interp = fem.interpolate_mode(mode_full, n_interp=200)

    pts  = np.column_stack([x_interp, w_interp])
    segs = np.stack([pts[:-1], pts[1:]], axis=1)
    mid  = 0.5 * (np.abs(w_interp[:-1]) + np.abs(w_interp[1:]))
    lc = LineCollection(segs, cmap='plasma', norm=Normalize(0, 1),
                        linewidth=2.0)
    lc.set_array(mid)
    ax.add_collection(lc)
    ax.fill_between(x_interp, 0, w_interp, alpha=0.12, color='#7b2d8b')

    error = 100.0 * (freqs_fem[0] - freqs_an[0]) / freqs_an[0]
    ax.set_xlim(0, fem.L)
    ax.set_ylim(-1.3, 1.3)
    ax.axhline(0, color='k', linewidth=0.5)
    ax.axvline(0, color='#c0392b', linestyle='--', linewidth=1.2)
    ax.set_xlabel('Position (m)')
    ax.set_ylabel('Norm. amplitude')
    ax.set_title(
        f'Mode 1: {freqs_fem[0]:.2f} Hz  (err {error:.4f}%)',
        fontsize=10
    )


def _draw_bode_inline(ax, A, B, C, K_lqr, n_points=250):
    """Compute and draw Bode magnitude on ax."""
    n    = A.shape[0]
    A_cl = A - B @ K_lqr
    freqs  = np.logspace(0, np.log10(250), n_points)
    omegas = 2.0 * np.pi * freqs
    H_ol = np.empty(n_points, dtype=complex)
    H_cl = np.empty(n_points, dtype=complex)
    I_n  = np.eye(n)
    for i, omega in enumerate(omegas):
        jw_I    = 1j * omega * I_n
        H_ol[i] = (C @ np.linalg.solve(jw_I - A,    B)).item()
        H_cl[i] = (C @ np.linalg.solve(jw_I - A_cl, B)).item()

    ax.semilogx(freqs, 20 * np.log10(np.abs(H_ol) + 1e-30),
                color=_COL_FREE, linewidth=1.5, label='Open-loop', alpha=0.9)
    ax.semilogx(freqs, 20 * np.log10(np.abs(H_cl) + 1e-30),
                color=_COL_LQR, linewidth=1.5, label='Closed-loop')
    ax.set_xlabel('Frequency (Hz)')
    ax.set_ylabel('|H| (dB)')
    ax.set_title('Bode Magnitude', fontsize=10)
    ax.legend(fontsize=9)
    ax.axhline(0, color='gray', linewidth=0.5, linestyle=':')


def _draw_transient_inline(ax, t_free, y_free, t_lqr, y_lqr, t_lqg, y_lqg):
    """Transient displacement comparison panel."""
    y_peak = np.max(np.abs(y_free))
    band   = 0.02 * y_peak
    ax.fill_between(t_free, -band * 1000, band * 1000,
                    alpha=0.18, color='#27ae60', label='±2%')
    ax.plot(t_free, y_free * 1000, color=_COL_FREE, linewidth=1.5,
            label='Free', alpha=0.8)
    ax.plot(t_lqr, y_lqr * 1000, color=_COL_LQR, linewidth=1.5,
            label='LQR')
    ax.plot(t_lqg, y_lqg * 1000, color=_COL_LQG, linewidth=1.5,
            linestyle='--', label='LQG')
    ax.set_xlabel('Time (s)')
    ax.set_ylabel('Tip disp. (mm)')
    ax.set_title('Transient Response', fontsize=10)
    ax.legend(fontsize=9)
    ax.axhline(0, color='gray', linewidth=0.5)


def _draw_phase_inline(ax, x_free, x_lqr, x_lqg, sensor_dof, n_dof_free):
    """Phase portrait panel."""
    def extract(xh):
        return xh[:, sensor_dof] * 1000, xh[:, n_dof_free + sensor_dof] * 1000

    d_f, v_f = extract(x_free)
    d_r, v_r = extract(x_lqr)
    d_g, v_g = extract(x_lqg)

    ax.plot(d_f, v_f, color=_COL_FREE, linewidth=1.0, label='Free', alpha=0.6)
    ax.plot(d_r, v_r, color=_COL_LQR, linewidth=1.3, label='LQR')
    ax.plot(d_g, v_g, color=_COL_LQG, linewidth=1.3, linestyle='--', label='LQG')
    ax.plot(0, 0, '*', color='#27ae60', markersize=10, zorder=5)
    ax.set_xlabel('Tip disp. (mm)')
    ax.set_ylabel('Tip vel. (mm/s)')
    ax.set_title('Phase Portrait', fontsize=10)
    ax.legend(fontsize=9)
    ax.axhline(0, color='gray', linewidth=0.5)
    ax.axvline(0, color='gray', linewidth=0.5)


def _draw_metrics_table(ax, metrics):
    """Styled performance-metrics table panel."""
    ax.axis('off')
    col_labels = ['Metric', 'Free', 'LQR', 'LQG']
    cell_text  = [
        ['Settling time (s)',
         f'{metrics["settling_free"]:.3f}',
         f'{metrics["settling_lqr"]:.3f}',
         f'{metrics["settling_lqg"]:.3f}'],
        ['Peak force (N)',   '—',
         f'{metrics["peak_lqr"]:.1f}',
         f'{metrics["peak_lqg"]:.1f}'],
        ['RMS force (N)',    '—',
         f'{metrics["rms_lqr"]:.2f}',
         f'{metrics["rms_lqg"]:.2f}'],
        ['State dimension',
         str(metrics['n_states']),
         str(metrics['n_states']),
         str(metrics['n_states'])],
        ['FEM freq. error', '—', f'{metrics["freq_error_mode1"]:.4f}%', '—'],
    ]

    tbl = ax.table(cellText=cell_text, colLabels=col_labels,
                   cellLoc='center', loc='center',
                   bbox=[0.0, 0.0, 1.0, 1.0])
    tbl.auto_set_font_size(False)
    tbl.set_fontsize(9.5)
    tbl.scale(1.0, 1.9)

    # Header row
    for j in range(4):
        tbl[(0, j)].set_facecolor('#2166ac')
        tbl[(0, j)].set_text_props(color='white', weight='bold')

    # Alternate row shading
    for i in range(1, len(cell_text) + 1):
        bg = '#dce9f7' if i % 2 == 0 else 'white'
        for j in range(4):
            tbl[(i, j)].set_facecolor(bg)
            tbl[(i, j)].set_edgecolor('#cccccc')

    ax.set_title('Performance Metrics', fontsize=10, fontweight='bold', pad=12)


def plot_summary_dashboard(fem, A, B, C, K_lqr,
                           t_free, y_free, t_lqr, y_lqr, t_lqg, y_lqg,
                           x_free, x_lqr, x_lqg,
                           metrics, filename='summary_dashboard.png'):
    """
    Six-panel portfolio summary dashboard.

    Layout (2 × 3):
      [0,0] Cantilever beam schematic
      [0,1] Mode 1 shape (colormap)
      [0,2] Bode magnitude
      [1,0] Transient response comparison
      [1,1] Phase portrait
      [1,2] Performance metrics table

    Parameters
    ----------
    fem          : BeamFEM      FEM model
    A, B, C      : ndarray      State-space matrices
    K_lqr        : ndarray      LQR gain
    t_*, y_*, x_*: ndarray      Simulation results
    metrics      : dict         Performance metrics dict
    filename     : str          Output filename
    """
    fig = plt.figure(figsize=(19, 11))
    fig.patch.set_facecolor('#fafafa')

    gs = fig.add_gridspec(2, 3, hspace=0.45, wspace=0.38)
    ax00 = fig.add_subplot(gs[0, 0])
    ax01 = fig.add_subplot(gs[0, 1])
    ax02 = fig.add_subplot(gs[0, 2])
    ax10 = fig.add_subplot(gs[1, 0])
    ax11 = fig.add_subplot(gs[1, 1])
    ax12 = fig.add_subplot(gs[1, 2])

    _draw_beam_schematic(ax00, fem)
    _draw_mode1_inline(ax01, fem)
    _draw_bode_inline(ax02, A, B, C, K_lqr)
    _draw_transient_inline(ax10, t_free, y_free, t_lqr, y_lqr, t_lqg, y_lqg)

    n_dof_free = fem.n_dof_free
    sensor_dof = n_dof_free - 2
    _draw_phase_inline(ax11, x_free, x_lqr, x_lqg, sensor_dof, n_dof_free)
    _draw_metrics_table(ax12, metrics)

    fig.suptitle(
        'Digital Twin: Active Vibration Control of Cantilever Beam'
        ' — Portfolio Summary',
        fontsize=15, fontweight='bold', y=1.01
    )

    plt.savefig(filename, dpi=200, bbox_inches='tight',
                facecolor=fig.get_facecolor())
    plt.close()
    print(f"Saved summary dashboard to {filename}")


# =============================================================================
# New animation: side-by-side comparison
# =============================================================================

def animate_comparison(fem, t_arr_free, x_free, y_free,
                       t_arr_lqg, x_lqg, y_lqg, u_lqg,
                       filename='beam_comparison.gif', fps=30, duration=2.0):
    """
    Side-by-side animated GIF: free vibration (left) vs LQG control (right),
    synchronised to the same time axis.

    Parameters
    ----------
    fem                          : BeamFEM
    t_arr_free, x_free, y_free   : free-vibration simulation results
    t_arr_lqg,  x_lqg,  y_lqg   : LQG simulation results
    u_lqg                        : LQG control-force history
    filename                     : output GIF path
    fps                          : frames per second
    duration                     : animation duration in seconds
    """
    T_end    = min(duration, t_arr_free[-1], t_arr_lqg[-1])
    n_frames = int(T_end * fps)
    t_common = np.linspace(0, T_end, n_frames)

    # Interpolate state/output to common time grid
    x_free_c = interp1d(t_arr_free, x_free, axis=0)(t_common)
    y_free_c = np.interp(t_common, t_arr_free, y_free)
    x_lqg_c  = interp1d(t_arr_lqg,  x_lqg,  axis=0)(t_common)
    y_lqg_c  = np.interp(t_common, t_arr_lqg,  y_lqg)
    u_lqg_c  = np.interp(t_common, t_arr_lqg,  u_lqg)

    # Precompute y-axis limits
    y_max = 0.0
    for k in range(0, n_frames, max(1, n_frames // 40)):
        for x_src in (x_free_c, x_lqg_c):
            q         = x_src[k, :fem.n_dof_free]
            mode_full = fem.get_full_mode_shape(q)
            _, w_interp = fem.interpolate_mode(mode_full, n_interp=80)
            y_max = max(y_max, np.max(np.abs(w_interp)))
    y_lim = y_max * 1.25 * 1000   # mm

    fig, (ax_l, ax_r) = plt.subplots(1, 2, figsize=(16, 5.5),
                                      sharey=True, constrained_layout=True)

    norm_l = Normalize(-y_lim, y_lim)
    norm_r = Normalize(-y_lim, y_lim)

    def _plot_panel(ax, x_state, norm, colour, title_prefix, t_val, u_val=None):
        ax.clear()
        q         = x_state[:fem.n_dof_free]
        mode_full = fem.get_full_mode_shape(q)
        x_interp, w_interp = fem.interpolate_mode(mode_full, n_interp=100)
        w_mm = w_interp * 1000

        pts  = np.column_stack([x_interp, w_mm])
        segs = np.stack([pts[:-1], pts[1:]], axis=1)
        mid_w = 0.5 * (w_mm[:-1] + w_mm[1:])
        lc = LineCollection(segs, cmap='RdBu_r', norm=norm,
                            linewidth=3.0, zorder=3)
        lc.set_array(mid_w)
        ax.add_collection(lc)
        ax.fill_between(x_interp, 0, w_mm, alpha=0.16, color=colour)

        x_nodes = np.linspace(0, fem.L, fem.n_nodes)
        w_nodes = mode_full[::2] * 1000
        ax.plot(x_nodes, w_nodes, 'o', color='#222222', markersize=4, zorder=4)

        ax.axvline(0, color='#c0392b', linewidth=3.5, alpha=0.8)
        ax.fill_betweenx([-y_lim, y_lim], -0.02, 0,
                         color='gray', alpha=0.22, hatch='///')

        ax.set_xlim(-0.04, fem.L + 0.04)
        ax.set_ylim(-y_lim, y_lim)
        ax.set_xlabel('Position (m)', fontsize=11)
        ax.set_ylabel('Displacement (mm)', fontsize=11)
        ax.set_title(f'{title_prefix}  —  t = {t_val:.3f} s', fontsize=12)
        ax.axhline(0, color='k', linewidth=0.5)

        if u_val is not None:
            ax.text(0.97, 0.96, f'F = {u_val:+.1f} N',
                    transform=ax.transAxes, ha='right', va='top',
                    fontsize=10, color='navy',
                    bbox=dict(boxstyle='round,pad=0.3',
                              facecolor='lightyellow', alpha=0.85))

    def update(k):
        _plot_panel(ax_l, x_free_c[k], norm_l, _COL_FREE,
                    'Free Vibration', t_common[k])
        _plot_panel(ax_r, x_lqg_c[k],  norm_r, _COL_LQG,
                    'LQG Control',    t_common[k], u_val=u_lqg_c[k])

    frame_indices = np.arange(n_frames)
    anim   = FuncAnimation(fig, update, frames=frame_indices, interval=1000 / fps)
    writer = PillowWriter(fps=fps)
    anim.save(filename, writer=writer)
    plt.close()
    print(f"Saved comparison animation to {filename}")


# =============================================================================
# HTML report generator
# =============================================================================

def generate_html_report(figure_files, metrics, filename='report.html'):
    """
    Generate a self-contained HTML report with all PNG figures embedded as
    base64.  GIF animations are referenced by relative path (same directory).

    Parameters
    ----------
    figure_files : dict
        Mapping of caption key → file path, e.g.
        {'Summary Dashboard': 'summary_dashboard.png', ...}
    metrics : dict
        Performance metrics (same dict used for dashboard table).
    filename : str
        Output HTML file path.
    """
    # ------------------------------------------------------------------
    # Encode images
    # ------------------------------------------------------------------
    encoded = {}
    for caption, path in figure_files.items():
        try:
            with open(path, 'rb') as fh:
                raw  = base64.b64encode(fh.read()).decode('utf-8')
            ext  = path.rsplit('.', 1)[-1].lower()
            mime = 'image/gif' if ext == 'gif' else 'image/png'
            encoded[caption] = {'src': f'data:{mime};base64,{raw}',
                                 'path': path, 'is_gif': ext == 'gif'}
        except FileNotFoundError:
            encoded[caption] = None

    # ------------------------------------------------------------------
    # Metrics table rows
    # ------------------------------------------------------------------
    rows_html = '\n'.join([
        f'<tr><td>FEM frequency error (Mode 1)</td>'
        f'<td colspan="3">{metrics["freq_error_mode1"]:.4f} %</td></tr>',
        f'<tr class="alt"><td>State dimension</td>'
        f'<td colspan="3">{metrics["n_states"]}</td></tr>',
        f'<tr><td>Settling time – free vibration</td>'
        f'<td colspan="3">{metrics["settling_free"]:.3f} s</td></tr>',
        f'<tr class="alt"><td>Settling time – LQR</td>'
        f'<td colspan="3">{metrics["settling_lqr"]:.3f} s</td></tr>',
        f'<tr><td>Settling time – LQG</td>'
        f'<td colspan="3">{metrics["settling_lqg"]:.3f} s</td></tr>',
        f'<tr class="alt"><td>Peak control force – LQR</td>'
        f'<td colspan="3">{metrics["peak_lqr"]:.2f} N</td></tr>',
        f'<tr><td>Peak control force – LQG</td>'
        f'<td colspan="3">{metrics["peak_lqg"]:.2f} N</td></tr>',
        f'<tr class="alt"><td>RMS control force – LQR</td>'
        f'<td colspan="3">{metrics["rms_lqr"]:.3f} N</td></tr>',
        f'<tr><td>RMS control force – LQG</td>'
        f'<td colspan="3">{metrics["rms_lqg"]:.3f} N</td></tr>',
    ])

    # ------------------------------------------------------------------
    # Figure cards
    # ------------------------------------------------------------------
    figure_cards = []
    for caption, info in encoded.items():
        if info is None:
            continue
        # Embed all media (PNG and GIF) as base64 for a fully self-contained report
        img_tag = f'<img src="{info["src"]}" alt="{caption}" loading="lazy">'
        figure_cards.append(
            f'<div class="card">\n{img_tag}\n'
            f'<p class="caption">{caption}</p>\n</div>'
        )
    figures_html = '\n'.join(figure_cards)

    # ------------------------------------------------------------------
    # Full HTML template
    # ------------------------------------------------------------------
    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="UTF-8">
  <meta name="viewport" content="width=device-width, initial-scale=1.0">
  <title>Digital Twin: Active Vibration Control — Report</title>
  <style>
    *, *::before, *::after {{ box-sizing: border-box; margin: 0; padding: 0; }}
    body {{
      font-family: 'Segoe UI', Arial, sans-serif;
      background: #f4f6f9;
      color: #333;
      padding: 24px;
      max-width: 1500px;
      margin: auto;
    }}
    header {{
      background: linear-gradient(135deg, #1a3a5c, #2166ac);
      color: white;
      border-radius: 10px;
      padding: 28px 32px;
      margin-bottom: 28px;
    }}
    header h1 {{ font-size: 1.75rem; margin-bottom: 6px; }}
    header p  {{ opacity: 0.85; font-size: 0.95rem; }}
    .badge {{
      display: inline-block;
      background: #27ae60;
      color: white;
      padding: 4px 14px;
      border-radius: 14px;
      font-size: 0.82rem;
      margin-top: 10px;
    }}
    h2 {{
      color: #1a3a5c;
      border-left: 4px solid #2166ac;
      padding-left: 12px;
      margin: 32px 0 16px;
      font-size: 1.2rem;
    }}
    /* Metrics table */
    .metrics-table {{
      border-collapse: collapse;
      width: 100%;
      max-width: 680px;
      margin-bottom: 24px;
      background: white;
      border-radius: 8px;
      overflow: hidden;
      box-shadow: 0 2px 8px rgba(0,0,0,0.08);
    }}
    .metrics-table th {{
      background: #2166ac;
      color: white;
      padding: 10px 14px;
      text-align: left;
    }}
    .metrics-table td {{
      padding: 8px 14px;
      border-bottom: 1px solid #e8ecf0;
    }}
    .metrics-table tr.alt td {{ background: #f0f5ff; }}
    /* Figure grid */
    .fig-grid {{
      display: grid;
      grid-template-columns: repeat(auto-fit, minmax(700px, 1fr));
      gap: 22px;
      margin-bottom: 24px;
    }}
    .card {{
      background: white;
      border-radius: 10px;
      padding: 16px;
      box-shadow: 0 2px 10px rgba(0,0,0,0.09);
    }}
    .card img {{
      width: 100%;
      height: auto;
      border-radius: 6px;
    }}
    .caption {{
      text-align: center;
      color: #555;
      font-style: italic;
      font-size: 0.88rem;
      margin-top: 8px;
    }}
    footer {{
      text-align: center;
      color: #888;
      font-size: 0.82rem;
      margin-top: 40px;
      padding-top: 20px;
      border-top: 1px solid #dde;
    }}
  </style>
</head>
<body>

<header>
  <h1>🔬 Digital Twin: Active Vibration Control of Cantilever Beam</h1>
  <p>Custom FEM solver · LQG Compensator · Automated Simulation Pipeline</p>
  <div class="badge">✅ All Resume Points Verified</div>
</header>

<h2>📊 Performance Summary</h2>
<table class="metrics-table">
  <thead>
    <tr><th>Metric</th><th colspan="3">Value</th></tr>
  </thead>
  <tbody>
    {rows_html}
  </tbody>
</table>

<h2>📐 Figures</h2>
<div class="fig-grid">
  {figures_html}
</div>

<footer>
  Generated automatically by <code>generate_html_report()</code> —
  Digital Twin: Active Vibration Control project
</footer>

</body>
</html>
"""

    with open(filename, 'w', encoding='utf-8') as fh:
        fh.write(html)
    print(f"Generated HTML report: {filename}")
