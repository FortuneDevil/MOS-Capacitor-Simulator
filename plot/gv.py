import numpy as np
import matplotlib.pyplot as plt
import csv

# Function to read data from CSV file
# First column is VG, then one column of G/w (conductance over angular frequency) per frequency
def read_csv(filename):
    VG = []
    G = []

    # Read the CSV file and extract the data
    with open(filename, newline='') as csvfile:
        csvreader = csv.reader(csvfile)
        labels = next(csvreader)[1:]  # Header row
        for row in csvreader:
            VG.append(float(row[0]))
            G.append([float(value) for value in row[1:]])

    return np.array(VG), np.array(G), labels

vg, g, labels = read_csv("../csv_files/VG_G.csv")

# Plot the data
plt.figure(figsize=(8, 6))
for k, label in enumerate(labels):
    plt.plot(vg, g[:, k], marker='o', markersize=3, linestyle='-', label=label)

# Customize the plot
plt.title('Conductance / Angular Frequency vs Gate Voltage', fontsize=14)
plt.xlabel('Gate Voltage (V)', fontsize=12)
plt.ylabel('G/w (F/m^2)', fontsize=12)
plt.yscale("log")
plt.ylim(bottom=g.max() * 1e-6, top=g.max() * 3)
plt.grid(True, linestyle='--', alpha=0.7)
plt.legend(fontsize=9)
plt.tight_layout()

# Show the plot
# plt.show()
plt.savefig("../figs/gv.png")
