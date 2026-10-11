import numpy as np
import matplotlib.pyplot as plt
import csv

# Function to read data from CSV file
# First column is the time after the gate step, then one capacitance column per gate voltage
def read_csv(filename):
    t = []
    C = []

    # Read the CSV file and extract the data
    with open(filename, newline='') as csvfile:
        csvreader = csv.reader(csvfile)
        labels = next(csvreader)[1:]  # Header row
        for row in csvreader:
            t.append(float(row[0]))
            C.append([float(value) for value in row[1:]])

    return np.array(t), np.array(C), labels

t, c, labels = read_csv("../csv_files/t_C.csv")

# Plot the data
plt.figure(figsize=(8, 6))
for k, label in enumerate(labels):
    plt.plot(t, c[:, k], linestyle='-', label=label)

# Customize the plot
plt.title('Capacitance vs Time after a Gate Step', fontsize=14)
plt.xlabel('Time (s)', fontsize=12)
plt.ylabel('Capacitance (F/m^2)', fontsize=12)
plt.xscale("log")
plt.ylim(bottom=0)
plt.grid(True, linestyle='--', alpha=0.7)
plt.legend(fontsize=10)
plt.tight_layout()

# Show the plot
# plt.show()
plt.savefig("../figs/ct.png")
