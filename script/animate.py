import math
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.animation import FuncAnimation

# Read the output.csv file into a pandas dataframe
df = pd.read_csv("output.csv")

# Extract the time and state data
time = df["t"]
x = df["q1"]
theta = df["q2"]

# Create a new figure and set the aspect ratio
fig, ax = plt.subplots(figsize=(5, 5))
ax.set_aspect('equal', 'box')
ax.set_xlim(-2, 2)
ax.set_ylim(-2, 2)

# Plot the track
ax.plot([-100, 100], [0, 0], '--', lw = 1, color='black')

# Plot the cart as a rectangle
box_width = 0.5
box_height = 0.3
rect = plt.Rectangle((-box_width/2, -box_height/2), box_width, box_height, color='red')
ax.add_patch(rect)

# Plot the pendulum as a line
line, = ax.plot([0, 0], [0, 0], '-o', lw=4, color='black')

# Function to update the animation
def update(i):
    rect.set_x(x[i]-box_width/2)
    line.set_xdata([x[i], x[i] + 1 * math.sin(theta[i])])
    line.set_ydata([0, -1 * math.cos(theta[i])])

# Create animation
ani = FuncAnimation(fig, update, frames=len(time), interval=20, blit=False)

# Set the animation's interval based on the delta time for each frame
for i in range(1, len(time)):
    ani._step(time[i] - time[i-1])

# Show the animation
plt.show()
