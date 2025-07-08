import numpy as np
import matplotlib.pyplot as plt
import csv

# Function to read data from CSV file
def read_csv(filename):
    VG = []
    C = []

    # Read the CSV file and extract the data
    with open(filename, newline='') as csvfile:
        csvreader = csv.reader(csvfile)
        next(csvreader)  # Skip header row
        for row in csvreader:
            VG.append(float(row[0]))
            C.append(float(row[1]))
    
    return np.array(VG), np.array(C)

vg, c = read_csv("../csv_files/VG_C.csv")

# Plot the data
plt.figure(figsize=(8, 6))
plt.plot(vg, c, marker='o', linestyle='-', color='b', label='C vs VG')

# Customize the plot
plt.title('Capacitance vs Gate Voltage', fontsize=14)
plt.xlabel('Gate Voltage (V)', fontsize=12)
plt.ylabel('Capacitance (F)', fontsize=12)
plt.grid(True, linestyle='--', alpha=0.7)
plt.legend(fontsize=10)
plt.tight_layout()

# Show the plot
# plt.show()
plt.savefig("../figs/cv.png")
