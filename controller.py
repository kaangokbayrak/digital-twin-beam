"""
Control system design: LQR, Kalman Filter, and LQG compensator.
"""

import numpy as np
from scipy.linalg import solve_continuous_are


class LQRController:
    """
    Linear Quadratic Regulator for full-state feedback control.
    """
    
    def __init__(self, A, B, C=None, Q=None, R=None):
        """
        Design LQR controller.
        
        Parameters
        ----------
        A : ndarray
            State matrix
        B : ndarray
            Input matrix
        C : ndarray, optional
            Output matrix for Q = C^T C weighting
        Q : ndarray, optional
            State weighting matrix (default: C^T C * 1e6)
        R : ndarray, optional
            Control weighting matrix (default: identity)
        """
        self.A = A
        self.B = B
        self.C = C
        
        n_states = A.shape[0]
        n_inputs = B.shape[1]
        
        # Default Q: output-based weighting
        if Q is None:
            if C is not None:
                Q = C.T @ C * 1e6
            else:
                Q = np.eye(n_states)
        
        # Default R: unit control cost
        if R is None:
            R = np.eye(n_inputs)
        
        self.Q = Q
        self.R = R
        
        # Solve Continuous Algebraic Riccati Equation (CARE)
        self.P = solve_continuous_are(A, B, Q, R)
        
        # Compute optimal gain: K = R^{-1} B^T P
        self.K = np.linalg.solve(R, B.T @ self.P)
        
        # Closed-loop eigenvalues
        A_cl = A - B @ self.K
        self.eigenvalues_cl = np.linalg.eigvals(A_cl)
    
    def control_law(self, x):
        """
        Compute control input: u = -K*x
        
        Parameters
        ----------
        x : ndarray
            State vector
        
        Returns
        -------
        u : ndarray
            Control input
        """
        return -self.K @ x
    
    def print_summary(self):
        """
        Print LQR controller summary.
        """
        print("\n" + "="*60)
        print("LQR CONTROLLER SUMMARY")
        print("="*60)
        print(f"Gain matrix shape: {self.K.shape}")
        print(f"Gain norm: {np.linalg.norm(self.K):.6e}")
        print(f"\nClosed-loop eigenvalues (real part):")
        real_parts = np.real(self.eigenvalues_cl)
        print(f"  Max: {np.max(real_parts):.6e}")
        print(f"  Min: {np.min(real_parts):.6e}")
        
        if np.all(real_parts < 0):
            print("  ✓ All eigenvalues stable (negative real parts)")
        else:
            print("  ✗ Warning: Some eigenvalues unstable!")
        
        print("="*60)


class KalmanFilter:
    """
    Continuous-time Kalman filter for state estimation.
    """
    
    def __init__(self, A, B, C, Qn=None, Rn=None):
        """
        Design Kalman filter.
        
        Parameters
        ----------
        A : ndarray
            State matrix
        B : ndarray
            Input matrix
        C : ndarray
            Output matrix
        Qn : ndarray, optional
            Process noise covariance (default: 1e-3 * I)
        Rn : ndarray, optional
            Measurement noise covariance (default: 1e-6 * I)
        """
        self.A = A
        self.B = B
        self.C = C
        
        n_states = A.shape[0]
        n_outputs = C.shape[0]
        
        # Default noise covariances
        if Qn is None:
            Qn = 1e-3 * np.eye(n_states)
        if Rn is None:
            Rn = 1e-6 * np.eye(n_outputs)
        
        self.Qn = Qn
        self.Rn = Rn
        
        # Solve dual CARE for Kalman filter
        # Filter CARE: A*Pf + Pf*A^T - Pf*C^T*Rn^{-1}*C*Pf + Qn = 0
        self.Pf = solve_continuous_are(A.T, C.T, Qn, Rn)
        
        # Compute Kalman gain: L = Pf * C^T * Rn^{-1}
        self.L = self.Pf @ C.T @ np.linalg.inv(Rn)
        
        # Observer eigenvalues
        A_obs = A - self.L @ C
        self.eigenvalues_obs = np.linalg.eigvals(A_obs)
        
        # Initialize state estimate
        self.x_hat = None
        self.reset()
    
    def reset(self, x0=None):
        """
        Reset state estimate.
        
        Parameters
        ----------
        x0 : ndarray, optional
            Initial state estimate (default: zeros)
        """
        if x0 is None:
            n_states = self.A.shape[0]
            self.x_hat = np.zeros(n_states)
        else:
            self.x_hat = x0.copy()
    
    def update(self, y, u, dt):
        """
        Update state estimate using measurement.
        
        Uses RK4 integration for numerical stability.
        
        Parameters
        ----------
        y : ndarray
            Measurement
        u : ndarray
            Control input
        dt : float
            Time step
        
        Returns
        -------
        x_hat : ndarray
            Updated state estimate
        """
        # Observer dynamics function: x_hat_dot = A*x_hat + B*u + L*(y - C*x_hat)
        def observer_dynamics(x_hat_):
            y_pred = self.C @ x_hat_
            innovation = y - y_pred
            return self.A @ x_hat_ + self.B @ u + self.L @ innovation
        
        # RK4 integration for numerical stability with large gains
        k1 = observer_dynamics(self.x_hat)
        k2 = observer_dynamics(self.x_hat + 0.5 * dt * k1)
        k3 = observer_dynamics(self.x_hat + 0.5 * dt * k2)
        k4 = observer_dynamics(self.x_hat + dt * k3)
        self.x_hat = self.x_hat + (dt / 6.0) * (k1 + 2*k2 + 2*k3 + k4)
        
        return self.x_hat.copy()
    
    def print_summary(self):
        """
        Print Kalman filter summary.
        """
        print("\n" + "="*60)
        print("KALMAN FILTER SUMMARY")
        print("="*60)
        print(f"Gain matrix shape: {self.L.shape}")
        print(f"Gain norm: {np.linalg.norm(self.L):.6e}")
        print(f"\nObserver eigenvalues (real part):")
        real_parts = np.real(self.eigenvalues_obs)
        print(f"  Max: {np.max(real_parts):.6e}")
        print(f"  Min: {np.min(real_parts):.6e}")
        
        if np.all(real_parts < 0):
            print("  ✓ All eigenvalues stable (negative real parts)")
        else:
            print("  ✗ Warning: Some eigenvalues unstable!")
        
        print("="*60)


class LQGController:
    """
    LQG compensator combining LQR regulator and Kalman filter observer.
    """
    
    def __init__(self, lqr, kalman):
        """
        Assemble LQG controller.
        
        Parameters
        ----------
        lqr : LQRController
            LQR regulator
        kalman : KalmanFilter
            Kalman filter observer
        """
        self.lqr = lqr
        self.kalman = kalman
        
        # Reset observer
        self.kalman.reset()
    
    def step(self, y, dt):
        """
        Compute control from measurement and update observer.
        
        Parameters
        ----------
        y : ndarray
            Measurement
        dt : float
            Time step
        
        Returns
        -------
        u : ndarray
            Control input
        x_hat : ndarray
            State estimate
        """
        # Compute control from estimated state
        u = self.lqr.control_law(self.kalman.x_hat)
        
        # Update observer with measurement and control
        x_hat = self.kalman.update(y, u, dt)
        
        return u, x_hat
    
    def reset(self, x0=None):
        """
        Reset LQG controller.
        
        Parameters
        ----------
        x0 : ndarray, optional
            Initial state estimate
        """
        self.kalman.reset(x0)
    
    def print_summary(self):
        """
        Print LQG controller summary.
        """
        print("\n" + "="*60)
        print("LQG CONTROLLER SUMMARY")
        print("="*60)
        print("Combines LQR regulator with Kalman filter observer")
        print("Using separation principle")
        print("="*60)
        self.lqr.print_summary()
        self.kalman.print_summary()
