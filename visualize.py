import numpy as np
import matplotlib.pyplot as plt

points = np.genfromtxt("points.txt", delimiter=",")
print(points)

ptcs = np.linspace(0, len(points) - 1, len(points))
print(ptcs)
xptcs = points[:,0]; yptcs = points[:,1]
### Visualise the particle and grid

fig, ax = plt.subplots(figsize=(8,8))
plt.scatter(xptcs, yptcs)
plt.xlim = (0,8); plt.ylim = (0,12)
plt.grid()

for i, txt in enumerate(ptcs):
    ax.annotate(txt, (xptcs[i], yptcs[i]), fontsize=12)

#plt.show()
plt.savefig("points.png")
