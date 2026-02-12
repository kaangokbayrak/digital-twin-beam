"""
State-space representation for beam vibration with Rayleigh damping.
"""

import numpy as np


class StateSpace:
    """
    Convert second-order beam dynamics to first-order state-space form.
    
    Includes Rayleigh proportional damping.
    """
    
    def __init__(self, fem, zeta1, zeta2, actuator_dof=None, sensor_dof=None):
        """
        Initialize state-space model.
        
        Parameters
        ----------
        fem : BeamFEM
            Finite element model
        zeta1 : float
            Target damping ratio for first mode
        zeta2 : float
            Target damping ratio for second mode
        actuator_dof : int, optional
            DOF where force is applied (default: tip translation)
        sensor_dof : int, optional
            DOF to measure (default: tip translation)
        """
        self.fem = fem
        self.zeta1 = zeta1
        self.zeta2 = zeta2
        
        # Get first two natural frequencies for Rayleigh damping
        freqs, _ = fem.modal_analysis(2)
        omega1 = 2 * np.pi * freqs[0]
        omega2 = 2 * np.pi * freqs[1]
        
        # Compute Rayleigh damping coefficients
        # C = alpha*M + beta*K
        # Solving system for alpha and beta given zeta1, zeta2
        denom = omega2**2 - omega1**2
        self.alpha = 2 * omega1 * omega2 * (zeta1*omega2 - zeta2*omega1) / denom
        self.beta = 2 * (zeta2*omega2 - zeta1*omega1) / denom
        
        # Construct damping matrix
        C_damp = self.alpha * fem.M + self.beta * fem.K
        
        # State-space dimensions
        n = fem.n_dof_free
        self.n_states = 2 * n
        
        # Default actuator and sensor at tip translation DOF
        if actuator_dof is None:
            actuator_dof = n - 2  # Tip translation (last node, translation DOF)
        if sensor_dof is None:
            sensor_dof = n - 2  # Tip translation
        
        self.actuator_dof = actuator_dof
        self.sensor_dof = sensor_dof
        
        # Compute M^{-1}
        M_inv = np.linalg.inv(fem.M)
        
        # Build state-space matrices
        # State: x = [q; q_dot]
        # A = [[0, I], [-M^{-1}K, -M^{-1}C]]
        A_top = np.hstack([np.zeros((n, n)), np.eye(n)])
        A_bottom = np.hstack([-M_inv @ fem.K, -M_inv @ C_damp])
        self.A = np.vstack([A_top, A_bottom])
        
        # B: force at actuator_dof
        # B = [0; M^{-1} * B_u]
        B_u = np.zeros(n)
        B_u[actuator_dof] = 1.0
        self.B = np.vstack([np.zeros((n, 1)), (M_inv @ B_u).reshape(-1, 1)])
        
        # C: measure displacement at sensor_dof
        C_y = np.zeros(n)
        C_y[sensor_dof] = 1.0
        self.C = np.hstack([C_y, np.zeros(n)]).reshape(1, -1)
        
        # D: no direct feedthrough
        self.D = np.zeros((1, 1))
        
        # Store for reference
        self.C_damp = C_damp
        self.M_inv = M_inv
    
    def get_matrices(self):
        """
        Get state-space matrices.
        
        Returns
        -------
        A : ndarray
            State matrix
        B : ndarray
            Input matrix
        C : ndarray
            Output matrix
        D : ndarray
            Feedthrough matrix
        """
        return self.A, self.B, self.C, self.D
    
    def print_summary(self):
        """
        Print state-space model summary.
        """
        print("\n" + "="*60)
        print("STATE-SPACE MODEL SUMMARY")
        print("="*60)
        print(f"Number of states:     {self.n_states}")
        print(f"Number of inputs:     {self.B.shape[1]}")
        print(f"Number of outputs:    {self.C.shape[0]}")
        print(f"\nRayleigh damping:")
        print(f"  alpha = {self.alpha:.6e} (mass proportional)")
        print(f"  beta  = {self.beta:.6e} (stiffness proportional)")
        print(f"  Target zeta1 = {self.zeta1:.4f}")
        print(f"  Target zeta2 = {self.zeta2:.4f}")
        print(f"\nActuator at DOF: {self.actuator_dof}")
        print(f"Sensor at DOF:   {self.sensor_dof}")
        
        # Check controllability
        # Controllability matrix (simplified rank check for large systems)
        # Full controllability matrix is expensive for large n
        # Check rank of [B, AB, A^2B, ...]
        controllability_rank = 0
        ctrl_mat = self.B.copy()
        for i in range(min(self.n_states, 5)):  # Check first few
            controllability_rank = np.linalg.matrix_rank(ctrl_mat)
            if i < 4:
                ctrl_mat = np.hstack([ctrl_mat, np.linalg.matrix_power(self.A, i+1) @ self.B])
        
        print(f"\nControllability rank (approx): {controllability_rank} / {self.n_states}")
        if controllability_rank == self.n_states:
            print("  System is controllable")
        else:
            print("  System may not be fully controllable")
        
        print("="*60)
