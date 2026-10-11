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
- Small-signal AC analysis (linearised around each DC bias point) at any frequency.
- Low-frequency and high-frequency C–V curves from $`C = \mathrm{Re}(\tilde{Q}/\tilde{v})`$.
- Python scripts for plotting C–V characteristics and profiles.

---

## Physics and Numerical Methods Used

- **Poisson’s Equation**:
  ```math
  \frac{d}{dx} \left( \varepsilon \frac{dV}{dx} \right) = -\rho(x)
  ```
- **Carrier Continuity**:
  
```math
\frac{\partial n}{\partial t} = \frac{1}{q} \frac{\partial J_n}{\partial x} - R
```

```math
\frac{\partial p}{\partial t} = -\frac{1}{q} \frac{\partial J_p}{\partial x} - R
```

- **Recombination** (Shockley–Read–Hall):

```math
R = \frac{np - n_i^2}{\tau_p (n + n_i) + \tau_n (p + n_i)}
```
  
- **Discretization Methods**:
  - Finite volume for Poisson's equation
  - Scharfetter–Gummel for carrier transport
  - Small-signal (phasor) linearisation for the AC analysis, solved together with the linearised Poisson equation
- **Capacitance Calculation**: $`C = \mathrm{Re}\left(\tilde{Q}/\tilde{v}\right)`$, with $`\tilde{Q}`$ the AC gate charge for a gate signal $`\tilde{v}`$

More details are in the equations folder: `eq_corrected.pdf` (source `eq_corrected.tex`) has the equations exactly as implemented.

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

### Windows (MSYS2)

The steps above are for Linux. On Windows, open the **MSYS2 UCRT64** shell and use the makefile in the `windows` folder instead:

```bash
cd windows
make clean   # only needed after changing .hpp files
make run     # builds windows/moscap_sim.exe and runs it from the repo root
```

All Windows build output stays inside `windows/`, so the Linux build files in the repo root are not touched. Plotting works the same way as on Linux (`cd plot && ./plot.sh`).

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

   - This produces these output files in the csv_files folder:
     - `eq.csv`: V, n, p at equilibrium (VG = 0).
     - `dc.csv`: V, n, p at the last DC bias.
     - `ac.csv`: real and imaginary parts of the AC phasors of V, n, p at the last bias, per volt of gate signal.
     - `VG_C.csv`: capacitance per unit area versus gate voltage, one column per frequency (set by `freqs` in `main.cpp`).

4. **Post-process and plot (Python)**:

   - Go to the directory `plot' and run:

     ```bash
     ./plot.sh
     ```
     This generates the plots of all modes and CV plot in folder figs.

   - You may plot potential, carrier profiles, or C–V curves using the provided scripts individually.

---
