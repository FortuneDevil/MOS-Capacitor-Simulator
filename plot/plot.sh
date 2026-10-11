#!/bin/bash
# Makes every figure in figs/ from the data in csv_files/. Can be run from any folder.

set -e
cd "$(dirname "$0")"

echo "Plotting EQUILIBRIUM"
python3 plot_eq.py

echo "Plotting DC"
python3 plot_dc.py

echo "Plotting AC"
python3 plot_ac.py

echo "Plotting VG - C"
python3 cv.py

echo "Plotting VG - G"
python3 gv.py

echo "Plotting f - C, G"
python3 cf.py

echo "Plotting VG - C at different sweep rates"
python3 cv_sweep.py

echo "Plotting t - C"
python3 ct.py
