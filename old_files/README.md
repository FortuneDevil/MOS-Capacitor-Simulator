# Old files

Nothing in this folder is used by the simulator. It is kept for reference only.

| Folder | What it is |
|---|---|
| `c_version/` | The first version of the project, written in C, with its own makefile, plot scripts and results. |
| `equations/equations_original.pdf` | The original derivation (time-stepping AC analysis). Superseded by `equations/equations.pdf`. |
| `equations/eq_small_signal_draft.pdf` | First draft of the small-signal AC equations. It has sign errors in the DC recombination term and in `c_n`; the corrected version is `equations/equations.pdf`. |
| `tests/` | Small test programs for the old dense matrix solver. That solver has been removed, so they no longer compile against the current `math/` folder. |
| `vscode/` | Editor settings that were scattered over the repo. They contain paths of another machine. |

## What changed from the original C++ version

- **AC analysis:** one Crank–Nicolson time step replaced by a small-signal (phasor) analysis, solved together with the linearised Poisson equation.
- **Recombination:** `(n - n_eq)/tau` replaced by Shockley–Read–Hall.
- **Mesh:** the first semiconductor node was 10 nm from the interface instead of 0.1 nm (missing `- 1` in the mesh formula).
- **Poisson:** charge terms made consistent (factor `q`), interface node given half a control volume.
- **Solvers:** the dense Gaussian elimination replaced by tridiagonal solvers that do not lose the small recombination terms to round-off.
- **Added:** n-type substrates, gate work function, fixed oxide charge, conductance, frequency sweeps, deep depletion and transient simulations.
