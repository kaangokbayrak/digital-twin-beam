"""
Finite Element Method solver for Euler-Bernoulli beam vibration.

Implements 2-node beam elements with cubic Hermite shape functions.
"""

import numpy as np
from scipy.linalg import eigh


class BeamFEM:
    """
    Finite element model for a cantilever beam.
    
    Uses 2-node beam elements with 2 DOFs per node (translation and rotation).
    """
    
    def __init__(self, L, E, rho, b, h, n_elements=10):
        """
        Initialize beam FEM model.
        
        Parameters
        ----------
        L : float
            Beam length (m)
        E : float
            Young's modulus (Pa)
        rho : float
            Material density (kg/m^3)
        b : float
            Beam width (m)
        h : float
            Beam height (m)
        n_elements : int
            Number of finite elements
        """
        self.L = L
        self.E = E
        self.rho = rho
        self.b = b
        self.h = h
        self.n_elements = n_elements
        
        # Geometric properties
        self.A = b * h  # Cross-sectional area
        self.I = b * h**3 / 12  # Second moment of area
        self.Le = L / n_elements  # Element length
        
        # Number of nodes and DOFs
        self.n_nodes = n_elements + 1
        self.n_dof = 2 * self.n_nodes  # 2 DOFs per node
        
        # Global matrices (before BC application)
        self.K_global = None
        self.M_global = None
        
        # Reduced matrices (after BC application)
        self.K = None
        self.M = None
        self.n_dof_free = None
        
        # Build global matrices
        self._assemble_global()
        self._apply_bc()
        
    def _element_stiffness(self):
        """
        Compute element stiffness matrix for a beam element.
        
        Returns
        -------
        Ke : ndarray (4, 4)
            Element stiffness matrix
        """
        EI = self.E * self.I
        Le = self.Le
        
        Ke = (EI / Le**3) * np.array([
            [12,      6*Le,   -12,     6*Le],
            [6*Le,    4*Le**2, -6*Le,  2*Le**2],
            [-12,    -6*Le,    12,    -6*Le],
            [6*Le,    2*Le**2, -6*Le,  4*Le**2]
        ])
        
        return Ke
    
    def _element_mass(self):
        """
        Compute consistent mass matrix for a beam element.
        
        Returns
        -------
        Me : ndarray (4, 4)
            Element consistent mass matrix
        """
        rho = self.rho
        A = self.A
        Le = self.Le
        
        Me = (rho * A * Le / 420) * np.array([
            [156,     22*Le,    54,      -13*Le],
            [22*Le,   4*Le**2,  13*Le,   -3*Le**2],
            [54,      13*Le,    156,     -22*Le],
            [-13*Le, -3*Le**2, -22*Le,   4*Le**2]
        ])
        
        return Me
    
    def _assemble_global(self):
        """
        Assemble global stiffness and mass matrices using scatter-add.
        """
        # Initialize global matrices
        self.K_global = np.zeros((self.n_dof, self.n_dof))
        self.M_global = np.zeros((self.n_dof, self.n_dof))
        
        # Get element matrices
        Ke = self._element_stiffness()
        Me = self._element_mass()
        
        # Loop over elements and assemble
        for e in range(self.n_elements):
            # Node indices for this element
            node1 = e
            node2 = e + 1
            
            # Global DOF indices
            dofs = [2*node1, 2*node1+1, 2*node2, 2*node2+1]
            
            # Scatter-add element matrices into global
            for i in range(4):
                for j in range(4):
                    self.K_global[dofs[i], dofs[j]] += Ke[i, j]
                    self.M_global[dofs[i], dofs[j]] += Me[i, j]
    
    def _apply_bc(self):
        """
        Apply cantilever boundary conditions (fixed at x=0).
        
        Removes DOFs 0 and 1 (translation and rotation at node 0).
        """
        # Fixed DOFs (node 0: translation and rotation)
        fixed_dofs = [0, 1]
        
        # Free DOFs (all others)
        free_dofs = [i for i in range(self.n_dof) if i not in fixed_dofs]
        
        # Extract reduced matrices
        self.K = self.K_global[np.ix_(free_dofs, free_dofs)]
        self.M = self.M_global[np.ix_(free_dofs, free_dofs)]
        self.n_dof_free = len(free_dofs)
    
    def modal_analysis(self, n_modes):
        """
        Perform modal analysis to find natural frequencies and mode shapes.
        
        Solves the generalized eigenvalue problem: K*phi = omega^2 * M * phi
        
        Parameters
        ----------
        n_modes : int
            Number of modes to extract
        
        Returns
        -------
        freqs : ndarray
            Natural frequencies in Hz
        modes : ndarray
            Mode shape matrix (each column is a mode)
        """
        # Solve generalized eigenvalue problem
        eigenvalues, eigenvectors = eigh(self.K, self.M)
        
        # Extract first n_modes
        omega_squared = eigenvalues[:n_modes]
        omega = np.sqrt(omega_squared)
        freqs = omega / (2 * np.pi)
        
        modes = eigenvectors[:, :n_modes]
        
        return freqs, modes
    
    def get_full_mode_shape(self, mode_vec):
        """
        Convert reduced mode shape to full DOF vector (prepend zeros for fixed DOFs).
        
        Parameters
        ----------
        mode_vec : ndarray
            Mode shape in reduced coordinates
        
        Returns
        -------
        mode_full : ndarray
            Mode shape in full coordinates
        """
        mode_full = np.zeros(self.n_dof)
        mode_full[2:] = mode_vec  # Fixed DOFs are 0 and 1
        return mode_full
    
    def interpolate_mode(self, mode_full, n_interp=100):
        """
        Interpolate mode shape using Hermite shape functions for smooth visualization.
        
        Parameters
        ----------
        mode_full : ndarray
            Mode shape in full coordinates
        n_interp : int
            Number of interpolation points
        
        Returns
        -------
        x_interp : ndarray
            Interpolated x-coordinates
        w_interp : ndarray
            Interpolated displacement values
        """
        x_interp = np.linspace(0, self.L, n_interp)
        w_interp = np.zeros(n_interp)
        
        for i, x in enumerate(x_interp):
            # Find which element this point belongs to
            elem_idx = int(x / self.Le)
            if elem_idx >= self.n_elements:
                elem_idx = self.n_elements - 1
            
            # Local coordinate within element [0, 1]
            x_local = x - elem_idx * self.Le
            xi = x_local / self.Le
            
            # Node indices
            node1 = elem_idx
            node2 = elem_idx + 1
            
            # DOF values
            w1 = mode_full[2*node1]
            theta1 = mode_full[2*node1 + 1]
            w2 = mode_full[2*node2]
            theta2 = mode_full[2*node2 + 1]
            
            # Hermite shape functions
            N1 = 1 - 3*xi**2 + 2*xi**3
            N2 = self.Le * (xi - 2*xi**2 + xi**3)
            N3 = 3*xi**2 - 2*xi**3
            N4 = self.Le * (-xi**2 + xi**3)
            
            # Interpolate displacement
            w_interp[i] = N1*w1 + N2*theta1 + N3*w2 + N4*theta2
        
        return x_interp, w_interp
    
    def analytical_frequencies(self, n_modes):
        """
        Compute analytical natural frequencies for cantilever beam.
        
        Uses closed-form solution with beta_n*L values.
        
        Parameters
        ----------
        n_modes : int
            Number of modes
        
        Returns
        -------
        freqs_analytical : ndarray
            Analytical natural frequencies in Hz
        """
        # Beta_n * L values for cantilever beam
        beta_L = np.array([1.8751, 4.6941, 7.8548, 10.9955, 14.1372])
        
        # Ensure we have enough values
        if n_modes > len(beta_L):
            # Extend using asymptotic formula for higher modes
            for n in range(len(beta_L) + 1, n_modes + 1):
                beta_L = np.append(beta_L, (2*n - 1) * np.pi / 2)
        
        # Compute frequencies
        freqs_analytical = np.zeros(n_modes)
        for i in range(n_modes):
            freqs_analytical[i] = (beta_L[i]**2 / (2 * np.pi * self.L**2)) * \
                                  np.sqrt(self.E * self.I / (self.rho * self.A))
        
        return freqs_analytical
