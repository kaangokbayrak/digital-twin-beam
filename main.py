"""
Main orchestrator for Digital Twin beam vibration control.
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from fem_solver import BeamFEM
from state_space import StateSpace
from controller import LQRController, KalmanFilter, LQGController
from simulate import (simulate_free, simulate_lqr, simulate_lqg,
                      plot_mode_shapes, plot_transient_comparison, animate_beam)


def main():
    """
    Main execution function for Digital Twin beam project.
    """
    # Set random seed for reproducibility
    np.random.seed(42)
    
    print("="*60)
    print("DIGITAL TWIN: ACTIVE VIBRATION CONTROL OF CANTILEVER BEAM")
    print("="*60)
    
    # =========================================================================
    # 1. Define steel cantilever beam parameters
    # =========================================================================
    print("\n" + "="*60)
    print("STEP 1: BEAM GEOMETRY AND MATERIAL PROPERTIES")
    print("="*60)
    
    L = 1.0          # Length (m)
    E = 210e9        # Young's modulus (Pa)
    rho = 7850       # Density (kg/m^3)
    b = 0.03         # Width (m)
    h = 0.005        # Height (m)
    n_elements = 10  # Number of elements
    
    print(f"Length:           L = {L} m")
    print(f"Young's modulus:  E = {E:.2e} Pa")
    print(f"Density:          ρ = {rho} kg/m³")
    print(f"Width:            b = {b} m")
    print(f"Height:           h = {h} m")
    print(f"Cross-section:    A = {b*h:.2e} m²")
    print(f"Moment of inertia: I = {b*h**3/12:.2e} m⁴")
    print(f"Number of elements: {n_elements}")
    
    # =========================================================================
    # 2. Build FEM model
    # =========================================================================
    print("\n" + "="*60)
    print("STEP 2: FINITE ELEMENT MODEL")
    print("="*60)
    
    fem = BeamFEM(L, E, rho, b, h, n_elements)
    
    print(f"Number of nodes:           {fem.n_nodes}")
    print(f"Total DOFs:                {fem.n_dof}")
    print(f"Free DOFs (after BC):      {fem.n_dof_free}")
    print(f"Element length:            {fem.Le} m")
    
    # =========================================================================
    # 3. Modal analysis
    # =========================================================================
    print("\n" + "="*60)
    print("STEP 3: MODAL ANALYSIS")
    print("="*60)
    
    n_modes = 5
    freqs_fem, modes = fem.modal_analysis(n_modes)
    freqs_analytical = fem.analytical_frequencies(n_modes)
    
    print(f"\nNatural Frequencies (first {n_modes} modes):")
    print("-" * 60)
    print(f"{'Mode':<6} {'FEM (Hz)':<12} {'Analytical (Hz)':<18} {'Error (%)':<10}")
    print("-" * 60)
    
    for i in range(n_modes):
        error_pct = 100 * (freqs_fem[i] - freqs_analytical[i]) / freqs_analytical[i]
        print(f"{i+1:<6} {freqs_fem[i]:<12.4f} {freqs_analytical[i]:<18.4f} {error_pct:<10.3f}")
    
    # =========================================================================
    # 4. Plot mode shapes
    # =========================================================================
    print("\n" + "="*60)
    print("STEP 4: VISUALIZE MODE SHAPES")
    print("="*60)
    
    plot_mode_shapes(fem, n_modes=3, filename='mode_shapes.png')
    
    # =========================================================================
    # 5. State-space conversion
    # =========================================================================
    print("\n" + "="*60)
    print("STEP 5: STATE-SPACE CONVERSION")
    print("="*60)
    
    zeta1 = 0.01  # 1% damping on first mode
    zeta2 = 0.01  # 1% damping on second mode
    
    ss = StateSpace(fem, zeta1, zeta2)
    A, B, C, D = ss.get_matrices()
    ss.print_summary()
    
    # =========================================================================
    # 6. LQR controller design
    # =========================================================================
    print("\n" + "="*60)
    print("STEP 6: LQR CONTROLLER DESIGN")
    print("="*60)
    
    # Design LQR with output-based weighting
    Q_lqr = C.T @ C * 1e8  # Penalize tip displacement heavily
    R_lqr = np.array([[1.0]])  # Control effort weight
    
    print(f"Q weight (output-based): C^T C × 1e8")
    print(f"R weight: {R_lqr[0,0]}")
    
    lqr = LQRController(A, B, C, Q=Q_lqr, R=R_lqr)
    lqr.print_summary()
    
    # =========================================================================
    # 7. Kalman filter design
    # =========================================================================
    print("\n" + "="*60)
    print("STEP 7: KALMAN FILTER DESIGN")
    print("="*60)
    
    noise_std = 1e-4  # 1% of typical tip displacement (~10mm range)
    Rn = np.array([[noise_std**2]])  # Rn = variance = noise_std^2
    Qn = 1e-3 * np.eye(A.shape[0])  # Process noise covariance - smaller for stability  
    
    print(f"Process noise covariance:     Qn = 1e-3 × I")
    print(f"Measurement noise variance:   Rn = {noise_std**2:.2e} (noise_std = {noise_std:.2e})")
    
    kalman = KalmanFilter(A, B, C, Qn=Qn, Rn=Rn)
    kalman.print_summary()
    
    # =========================================================================
    # 8. Assemble LQG controller
    # =========================================================================
    print("\n" + "="*60)
    print("STEP 8: LQG CONTROLLER")
    print("="*60)
    
    lqg = LQGController(lqr, kalman)
    print("LQG controller assembled using separation principle")
    
    # =========================================================================
    # 9. Define initial condition
    # =========================================================================
    print("\n" + "="*60)
    print("STEP 9: INITIAL CONDITION")
    print("="*60)
    
    # Static tip deflection from 10N tip force
    F_tip = 10.0  # N
    F_global = np.zeros(fem.n_dof_free)
    F_global[-2] = F_tip  # Apply to tip translation DOF
    
    # Solve K*q = F for static deflection
    q0 = np.linalg.solve(fem.K, F_global)
    x0 = np.concatenate([q0, np.zeros(fem.n_dof_free)])  # [q; q_dot]
    
    tip_deflection = q0[-2]
    print(f"Applied tip force: {F_tip} N")
    print(f"Static tip deflection: {tip_deflection*1000:.3f} mm")
    print(f"Initial state dimension: {x0.shape[0]}")
    
    # =========================================================================
    # 10. Simulate all scenarios
    # =========================================================================
    print("\n" + "="*60)
    print("STEP 10: SIMULATE VIBRATION CONTROL")
    print("="*60)
    
    T_sim = 2.0  # Simulation time (s)
    dt_eval = 0.001  # Evaluation time step
    t_eval = np.arange(0, T_sim, dt_eval)
    
    print(f"Simulation time: {T_sim} s")
    print(f"Time step: {dt_eval} s")
    print(f"Number of time steps: {len(t_eval)}")
    
    # Scenario 1: Free vibration
    print("\n  Simulating free vibration...")
    t_free, x_free, y_free = simulate_free(A, x0, (0, T_sim), t_eval, C)
    
    # Scenario 2: LQR control
    print("  Simulating LQR control...")
    t_lqr, x_lqr, y_lqr, u_lqr = simulate_lqr(A, B, lqr.K, x0, (0, T_sim), t_eval, C)
    
    # Scenario 3: LQG control with noise
    print("  Simulating LQG control with measurement noise...")
    dt_lqg = 0.001  # Time step for LQG simulation
    t_lqg, x_lqg, x_hat_lqg, y_lqg, u_lqg = simulate_lqg(
        A, B, C, lqg, x0, (0, T_sim), dt_lqg, noise_std=noise_std
    )
    
    # =========================================================================
    # 11. Generate plots
    # =========================================================================
    print("\n" + "="*60)
    print("STEP 11: GENERATE PLOTS")
    print("="*60)
    
    plot_transient_comparison(
        t_free, y_free, t_lqr, y_lqr, u_lqr,
        t_lqg, y_lqg, u_lqg,
        filename='transient_response.png'
    )
    
    # =========================================================================
    # 12. Generate animations
    # =========================================================================
    print("\n" + "="*60)
    print("STEP 12: GENERATE ANIMATIONS")
    print("="*60)
    
    print("  Creating free vibration animation...")
    animate_beam(fem, t_free, x_free, 
                 filename='beam_free.gif', 
                 title='Free Vibration', 
                 fps=30, duration=2.0)
    
    print("  Creating LQR control animation...")
    animate_beam(fem, t_lqr, x_lqr, 
                 filename='beam_lqr.gif', 
                 title='LQR Control', 
                 fps=30, duration=2.0)
    
    print("  Creating LQG control animation...")
    animate_beam(fem, t_lqg, x_lqg, 
                 filename='beam_lqg.gif', 
                 title='LQG Control', 
                 fps=30, duration=2.0)
    
    # =========================================================================
    # 13. Performance metrics
    # =========================================================================
    print("\n" + "="*60)
    print("STEP 13: PERFORMANCE METRICS")
    print("="*60)
    
    # Settling time (2% criterion)
    threshold = 0.02 * np.max(np.abs(y_free))
    
    # Free vibration - find last time exceeding threshold
    exceed_indices = np.where(np.abs(y_free) >= threshold)[0]
    settling_free = t_free[exceed_indices[-1]] if len(exceed_indices) > 0 else 0.0
    
    # LQR - find last time exceeding threshold
    exceed_indices = np.where(np.abs(y_lqr) >= threshold)[0]
    settling_lqr = t_lqr[exceed_indices[-1]] if len(exceed_indices) > 0 else 0.0
    
    # LQG - find last time exceeding threshold
    exceed_indices = np.where(np.abs(y_lqg) >= threshold)[0]
    settling_lqg = t_lqg[exceed_indices[-1]] if len(exceed_indices) > 0 else 0.0
    
    print(f"\nSettling Time (2% criterion):")
    print(f"  Free vibration: {settling_free:.3f} s")
    print(f"  LQR control:    {settling_lqr:.3f} s")
    print(f"  LQG control:    {settling_lqg:.3f} s")
    
    print(f"\nPeak Control Forces:")
    print(f"  LQR: {np.max(np.abs(u_lqr)):.3f} N")
    print(f"  LQG: {np.max(np.abs(u_lqg)):.3f} N")
    
    print(f"\nRMS Control Effort:")
    print(f"  LQR: {np.sqrt(np.mean(u_lqr**2)):.3f} N")
    print(f"  LQG: {np.sqrt(np.mean(u_lqg**2)):.3f} N")
    
    # =========================================================================
    # Summary
    # =========================================================================
    print("\n" + "="*60)
    print("SIMULATION COMPLETE")
    print("="*60)
    print("\nGenerated files:")
    print("  - mode_shapes.png           (Modal analysis)")
    print("  - transient_response.png    (Comparison plots)")
    print("  - beam_free.gif             (Free vibration animation)")
    print("  - beam_lqr.gif              (LQR control animation)")
    print("  - beam_lqg.gif              (LQG control animation)")
    print("\nPerformance Summary:")
    print(f"  Vibration suppression: {(1 - settling_lqr/settling_free)*100:.1f}% faster")
    print(f"  LQR vs LQG settling:   LQR={settling_lqr:.3f}s, LQG={settling_lqg:.3f}s")
    print("="*60)


if __name__ == "__main__":
    main()
