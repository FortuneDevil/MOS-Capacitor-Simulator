import numpy as np
import matplotlib.pyplot as plt
import csv

# Function to read data from CSV file
# Columns: x, then real and imaginary parts of the AC phasors of V, n, p (per volt of gate signal)
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
            V.append(complex(float(row[1]), float(row[2])))
            n.append(complex(float(row[3]), float(row[4])))
            p.append(complex(float(row[5]), float(row[6])))

    return np.array(x), np.array(V), np.array(n), np.array(p)

# Read data from the CSV file
x, V, n, p = read_csv("../csv_files/ac.csv")

# Plot the data
plt.figure(figsize=(10, 6))

# Plot AC Voltage vs x
plt.subplot(3, 1, 1)  # 3 rows, 1 column, 1st plot
plt.plot(x, V.real, label="Re", color="b")
plt.plot(x, V.imag, label="Im", color="b", linestyle="--")
plt.xlabel("x (m)")
plt.ylabel("AC Voltage (V/V)")
plt.title("AC Voltage vs x")
plt.xscale("log")
plt.grid(True)
plt.legend()

# Plot AC Electron Concentration vs x
plt.subplot(3, 1, 2)  # 3 rows, 1 column, 2nd plot
plt.plot(x, n.real, label="Re", color="g")
plt.plot(x, n.imag, label="Im", color="g", linestyle="--")
plt.xlabel("x (m)")
plt.ylabel("AC n (m^-3/V)")
plt.title("AC Electron Concentration vs x")
plt.xscale("log")
plt.yscale("symlog", linthresh=1e15)
plt.yticks([-1e24, -1e21, -1e18, 0, 1e18, 1e21, 1e24])
plt.grid(True)
plt.legend()

# Plot AC Hole Concentration vs x
plt.subplot(3, 1, 3)  # 3 rows, 1 column, 3rd plot
plt.plot(x, p.real, label="Re", color="r")
plt.plot(x, p.imag, label="Im", color="r", linestyle="--")
plt.xlabel("x (m)")
plt.ylabel("AC p (m^-3/V)")
plt.title("AC Hole Concentration vs x")
plt.xscale("log")
plt.yscale("symlog", linthresh=1e15)
plt.yticks([-1e24, -1e21, -1e18, 0, 1e18, 1e21, 1e24])
plt.grid(True)
plt.legend()

# Adjust layout and show the plot
plt.tight_layout()
#plt.show()
plt.savefig("../figs/ac.png")
