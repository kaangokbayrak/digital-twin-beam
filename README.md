# Digital Twin: Active Vibration Control of a Cantilever Beam

[![Python](https://img.shields.io/badge/python-3.10%2B-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

A complete **digital twin** implementation for active vibration control of a cantilever beam. This project uses custom finite element methods (FEM) to build a physics-based model of a steel beam, then designs LQR and LQG controllers to suppress unwanted vibrations.

> **New to control theory?** See the [Theory Manual](Theory_Manual.md) for a step-by-step walkthrough of all the math — written to be accessible to undergraduate engineering students.

---

## What Does This Project Do?

Imagine a thin steel beam clamped at one end (like a diving board). If you push the free end and let go, it oscillates for a long time. This project:

1. **Models the beam** using the Finite Element Method (FEM) — divides the beam into 10 small elements and writes equations of motion for each one.
2. **Finds natural frequencies** — the speeds at which the beam naturally wants to vibrate (like a guitar string has specific pitches).
3. **Designs a controller** that applies a force at the tip to kill the vibration quickly, using:
   - **LQR** (Linear Quadratic Regulator): knows the full state of the beam — works perfectly.
   - **LQG** (Linear Quadratic Gaussian): only gets a noisy sensor reading — uses a Kalman filter to estimate the full state first.
4. **Simulates all three scenarios** (free vibration, LQR, LQG) and generates plots and animations.

**Result**: The controller reduces settling time from ~2 seconds down to ~0.2 seconds — roughly **10× faster**.

---

## Quick Start

### Installation

```bash
git clone https://github.com/kaangokbayrak/digital-twin-beam.git
cd digital-twin-beam
pip install -r requirements.txt
```

### Run the Simulation

```bash
python main.py
```

This runs the full 19-step pipeline (takes 2–5 minutes). When done, open `docs/index.html` in your browser to see the interactive report with all plots and animations fully embedded.

---

## Project Structure

```
digital-twin-beam/
│
├── docs/
│   └── index.html             # Self-contained HTML report (open in browser)
│
├── fem_solver.py              # Finite element implementation
│   └── BeamFEM                #   Element matrices, assembly, modal analysis
│
├── state_space.py             # Converts second-order ODE to first-order state-space
│   └── StateSpace             #   Rayleigh damping, A/B/C/D matrices
│
├── controller.py              # Control system design
│   ├── LQRController          #   Optimal full-state feedback (solves CARE)
│   ├── KalmanFilter           #   State estimator (solves dual CARE)
│   └── LQGController          #   Combined LQR + Kalman filter
│
├── simulate.py                # Simulation engine and visualization
│   ├── simulate_free()        #   Free vibration (RK45)
│   ├── simulate_lqr()         #   LQR-controlled simulation (RK45)
│   ├── simulate_lqg()         #   LQG simulation with noise (fixed-step RK4)
│   ├── plot_mode_shapes()     #   Mode shape plots with frequency comparison
│   ├── plot_transient_comparison()  # Tip displacement, control effort, FFT
│   ├── plot_bode()            #   Bode magnitude — open-loop vs closed-loop
│   ├── plot_pole_zero()       #   Pole-zero map showing controller effect
│   ├── plot_phase_portrait()  #   Tip velocity vs displacement spirals
│   ├── plot_summary_dashboard()  # 6-panel portfolio figure
│   ├── animate_beam()         #   Animated GIF of beam deflecting
│   ├── animate_comparison()   #   Side-by-side free vs LQG animation
│   └── generate_html_report() #   Self-contained HTML with all results embedded
│
├── main.py                    # Master orchestrator (runs all 19 steps)
├── Theory_Manual.md           # Detailed math walkthrough for undergrads
├── requirements.txt           # Python dependencies
└── LICENSE
```

---

## The Physics in Plain English

### Step 1 – Beam Model (FEM)

The beam is divided into 10 **elements**, each with 2 nodes. Each node has 2 **degrees of freedom (DOFs)**: vertical displacement $w$ and rotation $\theta$.

- 10 elements → 11 nodes → 22 total DOFs
- Cantilever boundary condition: clamp the left end → remove 2 DOFs → **20 free DOFs**

Each element contributes a **stiffness matrix** $\mathbf{K}_e$ (how hard the element is to bend) and a **mass matrix** $\mathbf{M}_e$ (how much it resists acceleration). These get assembled into global $20\times20$ matrices $\mathbf{K}$ and $\mathbf{M}$.

### Step 2 – Modal Analysis

Solving the eigenvalue problem $\mathbf{K}\boldsymbol{\phi} = \omega^2\mathbf{M}\boldsymbol{\phi}$ gives:
- **Natural frequencies** $f_1, f_2, \ldots$ (in Hz) — how fast each mode oscillates
- **Mode shapes** $\boldsymbol{\phi}_i$ — what the beam looks like when vibrating in that mode

For this beam: $f_1 \approx 4.18$ Hz, $f_2 \approx 26.2$ Hz, $f_3 \approx 73.3$ Hz.

The FEM result matches the analytical formula to within **0.001%** — extremely accurate.

### Step 3 – State-Space Conversion

The 20-DOF second-order equation:
$$\mathbf{M}\ddot{\mathbf{q}} + \mathbf{C}\dot{\mathbf{q}} + \mathbf{K}\mathbf{q} = \mathbf{B}_u u$$

gets rewritten as a 40-dimensional first-order system:
$$\dot{\mathbf{x}} = \mathbf{A}\mathbf{x} + \mathbf{B}u, \qquad y = \mathbf{C}\mathbf{x}$$

where $\mathbf{x} = [\mathbf{q}^T, \dot{\mathbf{q}}^T]^T$ stacks positions and velocities.

Damping uses the **Rayleigh model**: $\mathbf{C} = \alpha\mathbf{M} + \beta\mathbf{K}$, tuned to give 1% damping on the first two modes.

### Step 4 – LQR Control

**Idea**: design a feedback gain $\mathbf{K}$ so that $u = -\mathbf{K}\mathbf{x}$ minimises:
$$J = \int_0^\infty \bigl(\mathbf{x}^T\mathbf{Q}\mathbf{x} + u^T R u\bigr)\,dt$$

Think of $\mathbf{Q}$ as "how much do I care about vibration?" and $R$ as "how much do I care about actuator effort?".

The optimal $\mathbf{K}$ is found by solving the **Continuous Algebraic Riccati Equation (CARE)** — a special matrix equation that `scipy.linalg.solve_continuous_are` handles.

### Step 5 – Kalman Filter

In reality, you only measure tip displacement with noise. The Kalman filter continuously estimates the full 40-state vector $\hat{\mathbf{x}}$ from the noisy measurement:
$$\dot{\hat{\mathbf{x}}} = \mathbf{A}\hat{\mathbf{x}} + \mathbf{B}u + \mathbf{L}(y - \mathbf{C}\hat{\mathbf{x}})$$

The gain $\mathbf{L}$ is chosen to minimise estimation error, and is found by solving another CARE (the "dual" problem).

### Step 6 – LQG = LQR + Kalman

The **separation principle** says: design the regulator and observer independently, then combine them. The resulting LQG controller uses:
$$u = -\mathbf{K}\hat{\mathbf{x}}$$

This works even with 1 mm sensor noise, and still achieves the same ~0.19 s settling time as the ideal LQR.

---

## Simulation Pipeline (19 Steps)

| Step | Description |
|------|-------------|
| 1 | Define beam geometry (L = 1 m, E = 210 GPa, 30 × 5 mm cross-section) |
| 2 | Build FEM model (10 elements → 20 free DOFs) |
| 3 | Modal analysis — extract natural frequencies and mode shapes |
| 4 | Plot mode shapes (`mode_shapes.png`) |
| 5 | Convert to state-space (40 states, Rayleigh damping) |
| 6 | Design LQR controller |
| 7 | Design Kalman filter |
| 8 | Assemble LQG compensator |
| 9 | Set initial condition (static 10 N tip deflection) |
| 10 | Simulate all three scenarios (free, LQR, LQG) |
| 11 | Plot transient response comparison (`transient_response.png`) |
| 12 | Generate beam animations (`beam_free.gif`, `beam_lqr.gif`, `beam_lqg.gif`) |
| 13 | Compute performance metrics (settling time, peak/RMS force) |
| 14 | Bode magnitude plot (`bode_plot.png`) |
| 15 | Pole-zero map (`pole_zero.png`) |
| 16 | Phase portrait (`phase_portrait.png`) |
| 17 | Summary dashboard (`summary_dashboard.png`) |
| 18 | Side-by-side comparison animation (`beam_comparison.gif`) |
| 19 | Generate self-contained HTML report (`docs/index.html`) |

---

## Typical Results

| Metric | Free Vibration | LQR Control | LQG Control |
|--------|---------------|-------------|-------------|
| Settling time (2% criterion) | ~2.0 s | ~0.19 s | ~0.19 s |
| Peak control force | — | ~498 N | ~498 N |
| RMS control force | — | ~17.4 N | ~15.3 N |
| First natural frequency | 4.178 Hz (FEM), 4.178 Hz (analytical) | — | — |
| FEM frequency error (mode 1) | < 0.001% | — | — |
| State dimension | 40 | 40 | 40 |

**~90% faster settling** with active control.

---

## Theory Overview

### Euler-Bernoulli Beam Equation

Governing PDE for transverse vibration of a uniform beam:

$$EI\frac{\partial^4 w}{\partial x^4} + \rho A \frac{\partial^2 w}{\partial t^2} = f(x,t)$$

### Element Stiffness Matrix

$$\mathbf{K}_e = \frac{EI}{L_e^3} \begin{bmatrix} 12 & 6L_e & -12 & 6L_e \\ 6L_e & 4L_e^2 & -6L_e & 2L_e^2 \\ -12 & -6L_e & 12 & -6L_e \\ 6L_e & 2L_e^2 & -6L_e & 4L_e^2 \end{bmatrix}$$

### Consistent Mass Matrix

$$\mathbf{M}_e = \frac{\rho A L_e}{420} \begin{bmatrix} 156 & 22L_e & 54 & -13L_e \\ 22L_e & 4L_e^2 & 13L_e & -3L_e^2 \\ 54 & 13L_e & 156 & -22L_e \\ -13L_e & -3L_e^2 & -22L_e & 4L_e^2 \end{bmatrix}$$

### State-Space Form

$$\mathbf{A} = \begin{bmatrix} \mathbf{0} & \mathbf{I} \\ -\mathbf{M}^{-1}\mathbf{K} & -\mathbf{M}^{-1}\mathbf{C} \end{bmatrix}, \quad \mathbf{B} = \begin{bmatrix} \mathbf{0} \\ \mathbf{M}^{-1}\mathbf{b}_u \end{bmatrix}$$

### LQR Cost Function

$$J = \int_0^\infty \left(\mathbf{x}^T\mathbf{Q}\mathbf{x} + u^T R u\right) dt$$

Optimal gain: $\mathbf{K} = R^{-1}\mathbf{B}^T\mathbf{P}$ where $\mathbf{P}$ solves the CARE.

### Kalman Filter Observer

$$\dot{\hat{\mathbf{x}}} = \mathbf{A}\hat{\mathbf{x}} + \mathbf{B}u + \mathbf{L}(y - \mathbf{C}\hat{\mathbf{x}})$$

Kalman gain: $\mathbf{L} = \mathbf{P}_f\mathbf{C}^T R_n^{-1}$ from the dual CARE.

See [Theory_Manual.md](Theory_Manual.md) for full derivations.

---

## Dependencies

| Package | Purpose |
|---------|---------|
| **NumPy** ≥ 1.24 | Matrix operations, linear algebra |
| **SciPy** ≥ 1.10 | Eigensolvers, CARE solver (`solve_continuous_are`), ODE integrators |
| **Matplotlib** ≥ 3.7 | Plotting, animations, colormaps |
| **Pillow** ≥ 9.0 | GIF animation encoding |

---

## References

1. Meirovitch, L. (1997). *Principles and Techniques of Vibrations*. Prentice Hall.
2. Bathe, K.J. (1996). *Finite Element Procedures*. Prentice Hall.
3. Lewis, F.L., Vrabie, D., Syrmos, V.L. (2012). *Optimal Control*, 3rd ed. Wiley.
4. Åström, K.J., Murray, R.M. (2008). *Feedback Systems: An Introduction for Scientists and Engineers*. Princeton University Press.

---

## License

MIT License — see [LICENSE](LICENSE).

---

*This is an educational/research project demonstrating control theory and FEM applied to structural vibration suppression.*
