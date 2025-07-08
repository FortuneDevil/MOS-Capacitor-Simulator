import numpy as np
import matplotlib.pyplot as plt
import csv

# Function to read data from CSV file
def read_csv(filename):
    x = []
    V = []
    n = []
    p = []
    
    # Read the CSV file and extract the data
    with open(filename, newline='') as csvfile:
        csvreader = csv.reader(csvfile)
        next(csvreader)  # Skip header row
        for row in csvreader:
            x.append(float(row[0]))
            V.append(float(row[1]))
            n.append(float(row[2]))
            p.append(float(row[3]))
    
    return np.array(x), np.array(V), np.array(n), np.array(p)

# Read data from the CSV file
x, V, n, p = read_csv("../csv_files/eq.csv")

# Plot the data
plt.figure(figsize=(10, 6))

# Plot Voltage vs x
plt.subplot(3, 1, 1)  # 3 rows, 1 column, 1st plot
plt.plot(x, V, label="Voltage (V)", color="b")
plt.xlabel("x (cm)")
plt.ylabel("Voltage (V)")
plt.title("Voltage vs x")
plt.grid(True)
plt.legend()

# Plot Electron Concentration vs x
plt.subplot(3, 1, 2)  # 3 rows, 1 column, 2nd plot
plt.plot(x, n, label="Electron Concentration (n)", color="g")
plt.xlabel("x (cm)")
plt.ylabel("Electron Concentration (n)")
plt.title("Electron Concentration vs x")
plt.grid(True)
plt.legend()

# Plot Hole Concentration vs x
plt.subplot(3, 1, 3)  # 3 rows, 1 column, 3rd plot
plt.plot(x, p, label="Hole Concentration (p)", color="r")
plt.xlabel("x (cm)")
plt.ylabel("Hole Concentration (p)")
plt.title("Hole Concentration vs x")
plt.grid(True)
plt.legend()

# Adjust layout and show the plot
plt.tight_layout()
# plt.show()
plt.savefig("../figs/eq.png")
