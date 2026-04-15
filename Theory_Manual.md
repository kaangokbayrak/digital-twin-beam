# Theory Manual: Active Vibration Control of Cantilever Beams

> **Audience**: Undergraduate engineering students (dynamics, vibrations, or controls courses).
> Each section starts with an intuitive explanation before introducing the math.

---

## 1. Euler-Bernoulli Beam Theory

### 1.1 Governing Equation

**Intuition**: When you push down on the free end of a cantilever beam and release it, why does it vibrate? The beam stores energy in two ways: bending stiffness resists deformation (like a spring), and mass resists acceleration (like inertia). The balance between these two produces oscillation.

The **Euler-Bernoulli beam model** assumes that cross-sections remain planar and perpendicular to the neutral axis during bending (true for slender beams where length >> thickness). This gives a 4th-order PDE:

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

**Intuition**: We need to approximate the beam's curved shape using polynomials. A cubic polynomial has 4 coefficients, which we match to 2 DOFs at each end of the element (displacement $w$ and rotation $\theta$). This gives a smooth curve that automatically satisfies continuity at element boundaries.

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

**Intuition**: The stiffness matrix tells you how much force/moment is needed at each node to produce unit displacement/rotation at any other node. It is derived from the beam's strain energy — the energy stored when the beam bends.

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

**Intuition**: The mass matrix tells you how forces at each node relate to accelerations elsewhere. Using the same shape functions for mass as for stiffness gives the "consistent" formulation, which is more accurate than lumping mass at nodes. The factor 1/420 comes from integrating products of cubic shape functions.

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

**Intuition**: Each element "knows" about 4 local DOFs (displacement and rotation at its two nodes). The assembly step maps these into the correct positions in the global matrix — like a jigsaw puzzle where adjacent elements share nodes.

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

**Intuition**: Every structure has natural frequencies at which it vibrates freely. At these frequencies, all parts of the structure move in a fixed pattern (the mode shape) and at the same frequency. Finding these frequencies is critical because a controller must damp them all.

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

**Intuition**: Control theory tools (LQR, Kalman) work on first-order ODE systems $\dot{\mathbf{x}} = \mathbf{A}\mathbf{x} + \mathbf{B}u$. But our FEM gives a second-order system ($\mathbf{M}\ddot{\mathbf{q}} + \ldots$). The trick is simple: define a new "state vector" that stacks both positions and velocities. Then Newton's law becomes a first-order equation in this bigger state.

For 20 physical DOFs: we get a **40-state** system (20 positions + 20 velocities).

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

**Intuition**: Real beams lose energy due to internal friction. Rayleigh damping is the simplest model: damping is a linear combination of mass and stiffness matrices. The mass term damps low frequencies; the stiffness term damps high frequencies. We tune $\alpha$ and $\beta$ to set exactly 1% damping on the first two modes.

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

**Intuition**: We want to push the beam back to zero quickly (small $\mathbf{x}$) without burning out the actuator (small $u$). LQR frames this as a trade-off: minimise a weighted sum of vibration energy and control effort over all future time. The answer turns out to be a simple linear feedback $u = -\mathbf{K}\mathbf{x}$, and the optimal gain matrix $\mathbf{K}$ can be computed offline.

The matrices $\mathbf{Q}$ and $R$ are design choices:
- Larger $\mathbf{Q}$ → heavier penalty on vibration → more aggressive control, faster settling
- Larger $R$ → heavier penalty on actuator force → gentler control, slower settling

Minimize the quadratic cost functional:

$$
J = \int_0^\infty \left(\mathbf{x}^T\mathbf{Q}\mathbf{x} + \mathbf{u}^T\mathbf{R}\mathbf{u}\right) dt
$$

subject to: $\dot{\mathbf{x}} = \mathbf{A}\mathbf{x} + \mathbf{B}\mathbf{u}$

where:
- $\mathbf{Q} \geq 0$ = state weighting matrix (positive semi-definite)
- $\mathbf{R} > 0$ = control weighting matrix (positive definite)

### 6.2 Continuous-Time Algebraic Riccati Equation (CARE)

**Intuition**: The CARE is the key equation that LQR solves. It finds the matrix $\mathbf{P}$ that represents the "cost-to-go" from any state — the minimum future cost if you always act optimally. The optimal gain is then $\mathbf{K} = R^{-1}\mathbf{B}^T\mathbf{P}$. You solve this once (offline), then just multiply $\mathbf{K}\mathbf{x}$ in real-time.

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

**Intuition**: LQR requires the full 40-state vector $\mathbf{x}$, but in practice you only have one noisy sensor (tip displacement). The Kalman filter runs a parallel simulation of the beam in software, continuously correcting its estimate using the sensor measurement. The correction gain $\mathbf{L}$ is chosen to minimise estimation error.

Think of it as: "My model says the beam should be here, but my sensor says it's slightly different — I'll update my estimate by blending the two, weighted by how much I trust each."

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

**Intuition**: The Kalman filter problem is mathematically identical to the LQR problem, just with matrices transposed. This means we can use the same CARE solver — just swap $\mathbf{A} \to \mathbf{A}^T$ and $\mathbf{B} \to \mathbf{C}^T$.

The filter problem is the **dual** of the LQR problem:
- $\mathbf{A} \leftrightarrow \mathbf{A}^T$
- $\mathbf{B} \leftrightarrow \mathbf{C}^T$
- $\mathbf{Q} \leftrightarrow \mathbf{Q}_n$
- $\mathbf{R} \leftrightarrow \mathbf{R}_n$

This allows solving the filter CARE using standard CARE solvers.

## 8. Linear Quadratic Gaussian (LQG) Compensator

### 8.1 Separation Principle

**Intuition**: The separation principle is one of the most powerful results in control theory. It says you can design the regulator (LQR) and the observer (Kalman) completely independently, then plug them together and the combined system is still stable and optimal. This makes the design problem tractable — a 40-state system would be extremely hard to design as a whole, but splitting it into two 40-state problems solved independently is feasible.

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

**Intuition**: LQR has nice robustness guarantees (it can tolerate gain errors of up to 50% or phase errors of up to 60°). LQG loses these guarantees because the observer introduces additional dynamics. However, in practice LQG works very well, especially when the noise model is accurate.

Unlike LQR, LQG does **not** guarantee:
- Gain/phase margins
- Robustness to model uncertainty

However, it is optimal for:
- Minimizing the cost functional under Gaussian noise assumptions
- Separable design of regulator and observer

## 9. Implementation Notes

### 9.1 Numerical Considerations

**Stiff dynamics**: The Kalman filter observer has very fast poles (eigenvalues with large negative real parts, up to ~$-2\times10^5$). This means the observer dynamics are much faster than the beam dynamics. Such "stiff" systems require small time steps — we use $\Delta t = 10^{-5}$ s for the LQG simulation (vs. $10^{-3}$ s for the free/LQR simulations). RK4 integration is used throughout for numerical accuracy.

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
