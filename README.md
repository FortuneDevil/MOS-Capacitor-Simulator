# MOS-Capacitor-Simulator

This is a C++ based simulation tool to compute the **C–V characteristics of a MOS Capacitor**, using physically accurate numerical methods. Python is used for post-processing and plotting.

---
## Table of Contents
- [Features](#features)
- [Physics and Numerical Methods Used](#physics-and-numerical-methods-used)
- [Getting Started](#getting-started)
- [Run Simulation](#run-simulation)

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

More Details are in reports folder.

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

   - Modify the `main.cpp` and `main.hpp' to set device parameters, mesh, voltages, and simulation settings.

2. **Run the simulator**:

   - If you make any changes in any of the hpp files, or if you want to remove all and run again from starting,

     ```bash
     make clean
     make
     ```
   - If you just make changes in .cpp files, there is no need to run make clean.
     ```bash
     make
     ```
   
   ```bash
   ./moscap_sim
   ```

   - This produces output files such as `VG-C.csv`, `eq.csv`, `cc.csv` and `ac.csv` in csv_files folder.

4. **Post-process and plot (Python)**:

   - Go to the directory `plot' and run:

     ```bash
     ./plot.sh
     ```
     This generates the plots of all modes and CV plot in folder figs.

   - You may plot potential, carrier profiles, or C–V curves using the provided scripts individually.

---
