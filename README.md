# MOS-Capacitor-Simulator

This is a C++ based simulation tool to compute the **C–V characteristics of a MOS Capacitor**, using physically accurate numerical methods. Python is used for post-processing and plotting.

---
## Table of Contents
- [Features](#features)
- [Physics and Numerical Methods Used](#physics-and-numerical-methods-used)
- [Getting Started](#getting-started)
- [Run Simulation](#run-simulation)
- [Example Plots](#example-plots)
- [References](#references)

---

## Features

- Non-uniform meshing for oxide and semiconductor. More data points near the oxide-semiconductor interface and less points in bulk.
- Poisson equation solver for electrostatics.
- Drift-diffusion solver for carrier concentrations using Scharfetter–Gummel discretization.
- Time-dependent AC simulation using Crank–Nicolson.
- Capacitance extraction using $`\frac{dQ}{dV_G}`$.
- Python scripts for plotting C–V characteristics and profiles.

---

## Physics and Numerical Methods Used

- **Poisson’s Equation**:
  ```math
  \frac{d}{dx} \left( \varepsilon \frac{dV}{dx} \right) = -\rho(x)
  ```
- **Carrier Continuity**:
  
```math
\frac{\partial n}{\partial t} = \frac{1}{q} \frac{\partial J_n}{\partial x} - \frac{n - n_{eq}}{\tau_n}
```

```math
\frac{\partial p}{\partial t} = -\frac{1}{q} \frac{\partial J_p}{\partial x} - \frac{p - p_{eq}}{\tau_p}
```
  
- **Discretization Methods**:
  - Finite difference for Poisson's equation
  - Scharfetter–Gummel for carrier transport
  - Crank–Nicolson for time-domain simulation (AC)
- **Capacitance Calculation**: $`C = \frac{dQ}{dV_G}`$

---

## Getting Started

### Prerequisites

- C++17 or later
- Python 3.x (for post-processing and plotting)
- [matplotlib](https://matplotlib.org/) and [numpy](https://numpy.org/) for plotting

### Build Instructions

1. **Clone the repository**:
   ```bash
   git clone https://github.com/FortuneDevil/MOS-Capacitor-Simulator.git
   cd MOS-Capacitor-Simulator
   ```

2. **Compile the C++ code**:
   ```bash
   make
   ```

---

## Run Simulation

1. **Edit the input parameters**:

   - Modify the `main.c` or relevant input file to set device parameters, mesh, voltages, and simulation settings.

2. **Run the simulator**:

   ```bash
   ./mos_sim
   ```

   - This produces output files such as `CV_data.csv`, `profile_V.csv`, etc.

3. **Post-process and plot (Python)**:

   - Go to the `python` directory (if provided) and run:

     ```bash
     python3 plot/plot_cv.py
     ```

   - You may plot potential, carrier profiles, or C–V curves using the provided scripts.

---

## Example Plots

*(Insert sample plots of C–V characteristic and potential profile here)*

---

## References

---