# Theory Manual: Active Vibration Control of Cantilever Beams

## 1. Euler-Bernoulli Beam Theory

### 1.1 Governing Equation

The transverse vibration of a slender beam is governed by the Euler-Bernoulli partial differential equation:

$$
\frac{\partial^2}{\partial x^2}\left(EI\frac{\partial^2 w}{\partial x^2}\right) + \rho A \frac{\partial^2 w}{\partial t^2} = f(x,t)
$$

where:
- $w(x,t)$ = transverse displacement
- $E$ = Young's modulus
- $I$ = second moment of area
- $\rho$ = material density
- $A$ = cross-sectional area
- $f(x,t)$ = distributed load per unit length

For uniform beams with constant $EI$:

$$
EI\frac{\partial^4 w}{\partial x^4} + \rho A \frac{\partial^2 w}{\partial t^2} = f(x,t)
$$

### 1.2 Hermite Shape Functions

The beam element uses cubic Hermite polynomials to interpolate displacement and rotation. For an element of length $L_e$, the normalized coordinate $\xi = x/L_e \in [0,1]$ gives four shape functions:

**Translational DOF at node 1:**
$$
N_1(\xi) = 1 - 3\xi^2 + 2\xi^3
$$

**Rotational DOF at node 1:**
$$
N_2(\xi) = L_e(\xi - 2\xi^2 + \xi^3)
$$

**Translational DOF at node 2:**
$$
N_3(\xi) = 3\xi^2 - 2\xi^3
$$

**Rotational DOF at node 2:**
$$
N_4(\xi) = L_e(-\xi^2 + \xi^3)
$$

The element displacement field is:
$$
w(x) = N_1(\xi)w_1 + N_2(\xi)\theta_1 + N_3(\xi)w_2 + N_4(\xi)\theta_2
$$

### 1.3 Element Stiffness Matrix

The element stiffness matrix is derived from the strain energy:

$$
U_e = \frac{1}{2}\int_0^{L_e} EI \left(\frac{d^2w}{dx^2}\right)^2 dx
$$

Using the shape functions and computing second derivatives:

$$
\mathbf{K}_e = \frac{EI}{L_e^3} \begin{bmatrix}
12 & 6L_e & -12 & 6L_e \\
6L_e & 4L_e^2 & -6L_e & 2L_e^2 \\
-12 & -6L_e & 12 & -6L_e \\
6L_e & 2L_e^2 & -6L_e & 4L_e^2
\end{bmatrix}
$$

### 1.4 Element Consistent Mass Matrix

The kinetic energy of the element is:

$$
T_e = \frac{1}{2}\int_0^{L_e} \rho A \left(\frac{\partial w}{\partial t}\right)^2 dx
$$

Using consistent mass formulation (integrating shape function products):

$$
\mathbf{M}_e = \frac{\rho A L_e}{420} \begin{bmatrix}
156 & 22L_e & 54 & -13L_e \\
22L_e & 4L_e^2 & 13L_e & -3L_e^2 \\
54 & 13L_e & 156 & -22L_e \\
-13L_e & -3L_e^2 & -22L_e & 4L_e^2
\end{bmatrix}
$$

## 2. Global Assembly

### 2.1 Scatter-Add Algorithm

For $n$ elements with $n+1$ nodes, the global system has $2(n+1)$ degrees of freedom (2 per node: translation $w$ and rotation $\theta$).

The **scatter-add** algorithm maps local element DOFs to global DOFs:

```
For element e connecting nodes i and j:
  Local DOFs:  [w_1, θ_1, w_2, θ_2]
  Global DOFs: [w_i, θ_i, w_j, θ_j]
  
  For each local DOF pair (m, n):
    K_global[global_m, global_n] += K_e[local_m, local_n]
    M_global[global_m, global_n] += M_e[local_m, local_n]
```

This assembly process ensures:
- Continuity of displacement and slope across elements
- Equilibrium at interior nodes
- Proper accumulation of contributions from adjacent elements

### 2.2 DOF Mapping

For element $e$ spanning nodes $i$ and $i+1$:
- Local DOF 0 → Global DOF $2i$     (translation at node $i$)
- Local DOF 1 → Global DOF $2i+1$   (rotation at node $i$)
- Local DOF 2 → Global DOF $2(i+1)$ (translation at node $i+1$)
- Local DOF 3 → Global DOF $2(i+1)+1$ (rotation at node $i+1$)

## 3. Boundary Conditions

### 3.1 Cantilever Beam (Fixed-Free)

A cantilever beam is fixed at $x=0$ (node 0) and free at $x=L$ (node $n$).

**At the fixed end:**
- $w(0) = 0$ (zero displacement)
- $\theta(0) = 0$ (zero rotation)

**At the free end:**
- Shear force $V = 0$ (natural BC)
- Bending moment $M = 0$ (natural BC)

### 3.2 DOF Elimination Method

To apply the fixed boundary conditions:

1. **Identify constrained DOFs**: For a cantilever, DOFs 0 and 1 (node 0: translation and rotation)
2. **Remove rows and columns**: Delete rows and columns corresponding to constrained DOFs from $\mathbf{K}$ and $\mathbf{M}$
3. **Result**: Reduced system with $(2n+2-2) = 2n$ free DOFs

After applying BCs:
$$
\mathbf{K}_{reduced} \in \mathbb{R}^{2n \times 2n}, \quad \mathbf{M}_{reduced} \in \mathbb{R}^{2n \times 2n}
$$

## 4. Modal Analysis

### 4.1 Generalized Eigenvalue Problem

Free vibration (no damping, no forcing) leads to:

$$
\mathbf{M}\ddot{\mathbf{q}} + \mathbf{K}\mathbf{q} = \mathbf{0}
$$

Assuming harmonic motion $\mathbf{q}(t) = \boldsymbol{\phi}e^{i\omega t}$:

$$
\mathbf{K}\boldsymbol{\phi} = \omega^2 \mathbf{M}\boldsymbol{\phi}
$$

This is a **generalized eigenvalue problem** where:
- $\omega_i$ = natural frequencies (rad/s)
- $\boldsymbol{\phi}_i$ = mode shapes (eigenvectors)
- $f_i = \omega_i/(2\pi)$ = frequencies in Hz

### 4.2 Orthogonality Properties

Mode shapes satisfy:

$$
\boldsymbol{\phi}_i^T \mathbf{M} \boldsymbol{\phi}_j = \delta_{ij} \quad \text{(mass normalization)}
$$

$$
\boldsymbol{\phi}_i^T \mathbf{K} \boldsymbol{\phi}_j = \omega_i^2 \delta_{ij}
$$

where $\delta_{ij}$ is the Kronecker delta.

### 4.3 Analytical Solution for Cantilever

The exact natural frequencies of a uniform cantilever beam are:

$$
f_n = \frac{(\beta_n L)^2}{2\pi L^2}\sqrt{\frac{EI}{\rho A}}
$$

where $\beta_n L$ are roots of $\cos(\beta_n L)\cosh(\beta_n L) = -1$:

| Mode | $\beta_n L$ | Description |
|------|-------------|-------------|
| 1    | 1.8751      | First bending |
| 2    | 4.6941      | Second bending |
| 3    | 7.8548      | Third bending |
| 4    | 10.9955     | Fourth bending |
| 5    | 14.1372     | Fifth bending |

For higher modes $(n \geq 4)$: $\beta_n L \approx (2n-1)\pi/2$

## 5. State-Space Conversion

### 5.1 Second-Order to First-Order Form

The damped equation of motion is:

$$
\mathbf{M}\ddot{\mathbf{q}} + \mathbf{C}\dot{\mathbf{q}} + \mathbf{K}\mathbf{q} = \mathbf{B}_u\mathbf{u}
$$

Define state vector $\mathbf{x} = [\mathbf{q}^T, \dot{\mathbf{q}}^T]^T$ with $2n$ states:

$$
\dot{\mathbf{x}} = \mathbf{A}\mathbf{x} + \mathbf{B}\mathbf{u}
$$

$$
\mathbf{y} = \mathbf{C}\mathbf{x} + \mathbf{D}\mathbf{u}
$$

where:

$$
\mathbf{A} = \begin{bmatrix}
\mathbf{0} & \mathbf{I} \\
-\mathbf{M}^{-1}\mathbf{K} & -\mathbf{M}^{-1}\mathbf{C}
\end{bmatrix} \in \mathbb{R}^{2n \times 2n}
$$

$$
\mathbf{B} = \begin{bmatrix}
\mathbf{0} \\
\mathbf{M}^{-1}\mathbf{B}_u
\end{bmatrix} \in \mathbb{R}^{2n \times m}
$$

$$
\mathbf{C} = \begin{bmatrix}
\mathbf{C}_y & \mathbf{0}
\end{bmatrix} \in \mathbb{R}^{p \times 2n}
$$

### 5.2 Rayleigh Damping

Rayleigh (proportional) damping provides a simple model:

$$
\mathbf{C} = \alpha \mathbf{M} + \beta \mathbf{K}
$$

The modal damping ratio for mode $i$ is:

$$
\zeta_i = \frac{1}{2}\left(\frac{\alpha}{\omega_i} + \beta\omega_i\right)
$$

Given two target damping ratios $\zeta_1$ and $\zeta_2$ at frequencies $\omega_1$ and $\omega_2$:

$$
\begin{bmatrix}
\zeta_1 \\
\zeta_2
\end{bmatrix} = \frac{1}{2} \begin{bmatrix}
1/\omega_1 & \omega_1 \\
1/\omega_2 & \omega_2
\end{bmatrix} \begin{bmatrix}
\alpha \\
\beta
\end{bmatrix}
$$

Solving for $\alpha$ and $\beta$:

$$
\alpha = \frac{2\omega_1\omega_2(\zeta_1\omega_2 - \zeta_2\omega_1)}{\omega_2^2 - \omega_1^2}
$$

$$
\beta = \frac{2(\zeta_2\omega_2 - \zeta_1\omega_1)}{\omega_2^2 - \omega_1^2}
$$

## 6. Linear Quadratic Regulator (LQR)

### 6.1 Optimal Control Problem

Minimize the quadratic cost functional:

$$
J = \int_0^\infty \left(\mathbf{x}^T\mathbf{Q}\mathbf{x} + \mathbf{u}^T\mathbf{R}\mathbf{u}\right) dt
$$

subject to: $\dot{\mathbf{x}} = \mathbf{A}\mathbf{x} + \mathbf{B}\mathbf{u}$

where:
- $\mathbf{Q} \geq 0$ = state weighting matrix (positive semi-definite)
- $\mathbf{R} > 0$ = control weighting matrix (positive definite)

### 6.2 Continuous-Time Algebraic Riccati Equation (CARE)

The optimal control law is:

$$
\mathbf{u}^* = -\mathbf{K}\mathbf{x}, \quad \mathbf{K} = \mathbf{R}^{-1}\mathbf{B}^T\mathbf{P}
$$

where $\mathbf{P}$ satisfies the CARE:

$$
\mathbf{A}^T\mathbf{P} + \mathbf{P}\mathbf{A} - \mathbf{P}\mathbf{B}\mathbf{R}^{-1}\mathbf{B}^T\mathbf{P} + \mathbf{Q} = \mathbf{0}
$$

### 6.3 Weight Selection

**Output-based weighting:**
$$
\mathbf{Q} = \mathbf{C}^T\mathbf{C} \cdot q_{scale}
$$

This penalizes the measured output (e.g., tip displacement) rather than all states equally.

**Control effort:**
$$
\mathbf{R} = r_{scale} \cdot \mathbf{I}
$$

Increasing $r_{scale}$ reduces control authority, decreasing it increases aggressiveness.

### 6.4 Closed-Loop Stability

The closed-loop system is:

$$
\dot{\mathbf{x}} = (\mathbf{A} - \mathbf{B}\mathbf{K})\mathbf{x}
$$

The LQR solution guarantees:
- $\mathbf{A} - \mathbf{B}\mathbf{K}$ is stable (all eigenvalues have negative real parts)
- Gain margin: $\infty$ to 0.5 ($\pm 6$ dB)
- Phase margin: at least $60°$

## 7. Kalman Filter (Observer)

### 7.1 State Estimation Problem

Given noisy measurements:

$$
\dot{\mathbf{x}} = \mathbf{A}\mathbf{x} + \mathbf{B}\mathbf{u} + \mathbf{w}
$$

$$
\mathbf{y} = \mathbf{C}\mathbf{x} + \mathbf{v}
$$

where:
- $\mathbf{w} \sim \mathcal{N}(0, \mathbf{Q}_n)$ = process noise
- $\mathbf{v} \sim \mathcal{N}(0, \mathbf{R}_n)$ = measurement noise

### 7.2 Kalman Filter Equations

The continuous-time Kalman filter is:

$$
\dot{\hat{\mathbf{x}}} = \mathbf{A}\hat{\mathbf{x}} + \mathbf{B}\mathbf{u} + \mathbf{L}(\mathbf{y} - \mathbf{C}\hat{\mathbf{x}})
$$

where $\mathbf{L}$ is the Kalman gain:

$$
\mathbf{L} = \mathbf{P}_f\mathbf{C}^T\mathbf{R}_n^{-1}
$$

and $\mathbf{P}_f$ satisfies the filter CARE:

$$
\mathbf{A}\mathbf{P}_f + \mathbf{P}_f\mathbf{A}^T - \mathbf{P}_f\mathbf{C}^T\mathbf{R}_n^{-1}\mathbf{C}\mathbf{P}_f + \mathbf{Q}_n = \mathbf{0}
$$

### 7.3 Duality with LQR

The filter problem is the **dual** of the LQR problem:
- $\mathbf{A} \leftrightarrow \mathbf{A}^T$
- $\mathbf{B} \leftrightarrow \mathbf{C}^T$
- $\mathbf{Q} \leftrightarrow \mathbf{Q}_n$
- $\mathbf{R} \leftrightarrow \mathbf{R}_n$

This allows solving the filter CARE using standard CARE solvers.

## 8. Linear Quadratic Gaussian (LQG) Compensator

### 8.1 Separation Principle

The LQG controller combines:
1. **Regulator**: $\mathbf{u} = -\mathbf{K}\hat{\mathbf{x}}$ (LQR gain using estimated state)
2. **Observer**: Kalman filter for state estimation

The combined system is:

$$
\dot{\hat{\mathbf{x}}} = (\mathbf{A} - \mathbf{B}\mathbf{K} - \mathbf{L}\mathbf{C})\hat{\mathbf{x}} + \mathbf{L}\mathbf{y}
$$

$$
\mathbf{u} = -\mathbf{K}\hat{\mathbf{x}}
$$

### 8.2 Stability Properties

The **separation principle** states:
- The regulator and observer can be designed independently
- The closed-loop system is stable if both:
  - $\mathbf{A} - \mathbf{B}\mathbf{K}$ is stable (guaranteed by LQR)
  - $\mathbf{A} - \mathbf{L}\mathbf{C}$ is stable (guaranteed by Kalman design)

### 8.3 Augmented State Dynamics

Define estimation error: $\tilde{\mathbf{x}} = \mathbf{x} - \hat{\mathbf{x}}$

The augmented system is:

$$
\begin{bmatrix}
\dot{\mathbf{x}} \\
\dot{\tilde{\mathbf{x}}}
\end{bmatrix} = \begin{bmatrix}
\mathbf{A} - \mathbf{B}\mathbf{K} & \mathbf{B}\mathbf{K} \\
\mathbf{0} & \mathbf{A} - \mathbf{L}\mathbf{C}
\end{bmatrix} \begin{bmatrix}
\mathbf{x} \\
\tilde{\mathbf{x}}
\end{bmatrix}
$$

The eigenvalues are the union of:
- LQR closed-loop poles (from $\mathbf{A} - \mathbf{B}\mathbf{K}$)
- Observer poles (from $\mathbf{A} - \mathbf{L}\mathbf{C}$)

### 8.4 Performance Considerations

Unlike LQR, LQG does **not** guarantee:
- Gain/phase margins
- Robustness to model uncertainty

However, it is optimal for:
- Minimizing the cost functional under Gaussian noise assumptions
- Separable design of regulator and observer

## 9. Implementation Notes

### 9.1 Numerical Considerations

**Matrix inversion**: Compute $\mathbf{M}^{-1}$ once and reuse:
```python
M_inv = np.linalg.inv(M)
A_bottom = np.hstack([-M_inv @ K, -M_inv @ C])
```

**Eigenvalue solver**: Use `scipy.linalg.eigh` for symmetric matrices (faster and more accurate than `eig`).

**CARE solver**: Use `scipy.linalg.solve_continuous_are` which implements the Schur method.

### 9.2 Units and Scaling

Ensure consistent units throughout:
- Length: meters (m)
- Mass: kilograms (kg)
- Time: seconds (s)
- Force: Newtons (N)
- Stress: Pascals (Pa = N/m²)
- Frequency: Hertz (Hz) or rad/s

### 9.3 Verification Methods

1. **Modal analysis**: Compare FEM frequencies with analytical solution
2. **Orthogonality**: Check $\boldsymbol{\Phi}^T\mathbf{M}\boldsymbol{\Phi} \approx \mathbf{I}$
3. **Stability**: Verify all eigenvalues of $\mathbf{A}-\mathbf{B}\mathbf{K}$ have negative real parts
4. **Energy decay**: Confirm exponential decay in free vibration
5. **Control effort**: Monitor actuator saturation limits

## References

1. Meirovitch, L. (1997). *Principles and Techniques of Vibrations*. Prentice Hall.
2. Bathe, K.J. (1996). *Finite Element Procedures*. Prentice Hall.
3. Lewis, F.L., Vrabie, D., Syrmos, V.L. (2012). *Optimal Control*, 3rd ed. Wiley.
4. Astrom, K.J., Murray, R.M. (2008). *Feedback Systems*. Princeton University Press.

---

*This manual provides the theoretical foundation for the Digital Twin implementation of active vibration control using custom FEM with LQR/LQG control strategies.*
