"""
TestCase Dambreak

"""

boxL = 1.61
boxH = 0.8

fluidL = 0.6
fluidH = 0.3

dp = 0.01
#dp = 0.1
rho0 = 1000.
p0 = 0.
n_layers = 3

import numpy as np

boxL_n = round(boxL/dp)
boxH_n = round(boxH/dp)

fluidL_n = round(fluidL/dp)
fluidH_n = round(fluidH/dp)


### Generate fluid particles
fluid_rx = np.zeros(fluidL_n*fluidH_n)
fluid_ry = np.zeros(fluidL_n*fluidH_n)
fluid_rz = np.zeros(fluidL_n*fluidH_n)


for x in range(fluidL_n):
    for z in range(fluidH_n):
        fluid_rx[x*fluidH_n + z] = dp*(x + 1)
        fluid_ry[x*fluidH_n + z] = 0. #we use only 2D case
        fluid_rz[x*fluidH_n + z] = dp*(z + 1)


### Generate boundary particles
box_rx = []
box_ry = []
box_rz = []

# left wall
for l in range(n_layers):
    for z in range(boxH_n - 1):
        box_rx.append(0. - l*dp)
        box_ry.append(0.) #we use only 2D case
        box_rz.append((z+1)*dp)


# bottom wall
for l in range(n_layers):
    for x in range(boxL_n + (n_layers - 1)*2):
        box_rx.append((x-(n_layers - 1))*dp)
        box_ry.append(0.) #we use only 2D case
        box_rz.append(0. - l*dp)

x_last = box_rx[-1 -(n_layers - 1)] #due to discretisation, we need to save last value of bottom wall

# right wall - layer 1
for l in range(n_layers):
    for z in range(boxH_n -1):
        box_rx.append(x_last + dp*l)
        box_ry.append(0.) #we use only 2D case
        box_rz.append((z+1)*dp)

### Write fluid particles
with open("dambreak.ptcs", "w") as f:
    f.write(str(len(fluid_rx) + len(box_rx)) + "\n")
    for i in range(len(box_rx)):
        f.write(str(round(box_rx[i], 5)) + " " + str(round(box_rz[i], 5)) + " " + str(round(box_ry[i], 5)) + " " + \
                str(0.) + " " + str(0.) + " " + str(0.) + " " + \
                str(round(rho0, 5)) + " " + str(round(p0, 5)) + "\n")
    for i in range(len(fluid_rx)):
        f.write(str(round(fluid_rx[i], 5)) + " " + str(round(fluid_rz[i], 5)) + " " + str(round(fluid_ry[i], 5)) + " " + \
                str(0.) + " " + str(0.) + " " + str(0.) + " " + \
                str(round(rho0, 5)) + " " + str(round(p0, 5)) + "\n")


