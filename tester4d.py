import numpy as np
import matplotlib.pyplot as plt

# Define the parameters of the map
K1 = 0.5
K2 = 0.5
C = 0.05 #B
n_steps = 500

# Define the initial conditions of the two particles
"""a=0.55
b=0.10
c=0.005 #0.54 #0.005
d=0.01"""

x1_0 = 0.55
p1_0 = 0.1
x2_0 = 0.54
p2_0 = 0.01

# Define the arrays to store the position and momentum of each particle
x1 = np.zeros(n_steps)
p1 = np.zeros(n_steps)
x2 = np.zeros(n_steps)
p2 = np.zeros(n_steps)

# Initialize the arrays with the initial conditions
x1[0] = x1_0
p1[0] = p1_0
x2[0] = x2_0
p2[0] = p2_0

# Calculate the position and momentum of each particle at each time step
for i in range(1, n_steps):
    p1[i] = p1[i-1] + (K1/2*np.pi)*np.sin(x1[i-1]) + (C/2*np.pi)*(np.sin(p2[i-1] - p1[i-1])) 
    p1[i] = p1[i] %1 
    x1[i] = x1[i-1] + p1[i] %1
    x1[i] = x1[i] %1
    print(p1[i])

    p2[i] = p2[i-1] + (K2/2*np.pi)*np.sin(x2[i-1]) + (C/2*np.pi)*(np.sin(p1[i-1] - p2[i-1])) %1
    p2[i] = p2[i] %1 
    x2[i] = x2[i-1] + p2[i] %1
    x2[i] = x2[i] %1
# Plot the results
fig, ax = plt.subplots()
ax.plot(x1, p1, ".",'g', label='Particle 1')
ax.plot(x2, p2,'.', 'b', label='Particle 2')
ax.set_xlabel('x')
ax.set_ylabel('p')
ax.set_title('Coupled 2D Standard Map')
ax.legend()
plt.show()