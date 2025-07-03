import matplotlib.pyplot as plt

# Read data from file
vg = []
c = []

with open("C_vs_VG.txt", "r") as file:
    # Skip the header
    next(file)
    
    # Read the data line by line
    for line in file:
        values = line.split()
        vg.append(float(values[0]))
        c.append(float(values[1]))

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
plt.show()

