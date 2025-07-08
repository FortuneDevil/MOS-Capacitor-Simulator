#!/bin/bash

echo "Plotting EQUILIBRIUM"
python3 plot_eq.py

echo "Plotting DC"
python3 plot_dc.py

echo "Plotting AC"
python3 plot_ac.py

echo "Plotting VG - C"
python3 cv.py
