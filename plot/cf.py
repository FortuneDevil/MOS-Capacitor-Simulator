import numpy as np
import matplotlib.pyplot as plt
import csv

# Function to read data from CSV file
# First column is the frequency, then two columns per gate voltage: C and G/w
def read_csv(filename):
    f = []
    data = []

    # Read the CSV file and extract the data
    with open(filename, newline='') as csvfile:
        csvreader = csv.reader(csvfile)
        labels = next(csvreader)[1:]  # Header row
        for row in csvreader:
            f.append(float(row[0]))
            data.append([float(value) for value in row[1:]])

    return np.array(f), np.array(data), labels

f, data, labels = read_csv("../csv_files/f_C.csv")

# Plot the data
plt.figure(figsize=(8, 8))

# Plot Capacitance vs frequency (even columns)
plt.subplot(2, 1, 1)  # 2 rows, 1 column, 1st plot
for k in range(0, len(labels), 2):
    plt.plot(f, data[:, k], label=labels[k])
plt.title('Capacitance vs Frequency', fontsize=14)
plt.xlabel('Frequency (Hz)', fontsize=12)
plt.ylabel('Capacitance (F/m^2)', fontsize=12)
plt.xscale("log")
plt.ylim(bottom=0)
plt.grid(True, linestyle='--', alpha=0.7)
plt.legend(fontsize=9)

# Plot Conductance / angular frequency vs frequency (odd columns)
plt.subplot(2, 1, 2)  # 2 rows, 1 column, 2nd plot
for k in range(1, len(labels), 2):
    plt.plot(f, data[:, k], label=labels[k])
plt.title('Conductance / Angular Frequency vs Frequency', fontsize=14)
plt.xlabel('Frequency (Hz)', fontsize=12)
plt.ylabel('G/w (F/m^2)', fontsize=12)
plt.xscale("log")
plt.yscale("log")
plt.ylim(bottom=data[:, 1::2].max() * 1e-8)
plt.grid(True, linestyle='--', alpha=0.7)
plt.legend(fontsize=9)

# Adjust layout and show the plot
plt.tight_layout()
# plt.show()
plt.savefig("../figs/cf.png")
