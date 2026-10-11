import numpy as np
import matplotlib.pyplot as plt
import csv

# Function to read data from CSV file
# First column is VG, then one capacitance column per frequency (the header names them)
def read_csv(filename):
    VG = []
    C = []

    # Read the CSV file and extract the data
    with open(filename, newline='') as csvfile:
        csvreader = csv.reader(csvfile)
        labels = next(csvreader)[1:]  # Header row
        for row in csvreader:
            VG.append(float(row[0]))
            C.append([float(value) for value in row[1:]])

    return np.array(VG), np.array(C), labels

vg, c, labels = read_csv("../csv_files/VG_C.csv")

# Plot the data
plt.figure(figsize=(8, 6))
for k, label in enumerate(labels):
    if "deep depletion" in label:
        plt.plot(vg, c[:, k], linestyle='--', color='k', label=label)
    else:
        plt.plot(vg, c[:, k], marker='o', markersize=3, linestyle='-', label=label)

# Customize the plot
plt.title('Capacitance vs Gate Voltage', fontsize=14)
plt.xlabel('Gate Voltage (V)', fontsize=12)
plt.ylabel('Capacitance (F/m^2)', fontsize=12)
plt.ylim(bottom=0)
plt.grid(True, linestyle='--', alpha=0.7)
plt.legend(fontsize=10)
plt.tight_layout()

# Show the plot
# plt.show()
plt.savefig("../figs/cv.png")
