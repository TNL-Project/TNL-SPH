import numpy as np
import matplotlib.pyplot as plt

points = np.genfromtxt("points.txt", delimiter=",")
#print(points)

def plotParticles():


    ptcs = np.linspace(0, len(points) - 1, len(points))
    #print(ptcs)
    xptcs = points[:,0]; yptcs = points[:,1]
    ### Visualise the particle and grid

    fig, ax = plt.subplots(figsize=(8,8))
    plt.scatter(xptcs, yptcs)
    plt.xlim = (0,8); plt.ylim = (0,12)
    plt.grid()

    for i, txt in enumerate(ptcs):
        ax.annotate(int(txt), (xptcs[i], yptcs[i]), fontsize=14)

    for i in range(len(ptcs)):
        circle = plt.Circle((xptcs[i], yptcs[i]), 1, color='b', fill=False)
        ax.add_patch(circle)

    ax.set_xlim([0,8])
    ax.set_ylim([0,8])

    plt.savefig("points.png")

    #: delete :txt = np.loadtxt("neighbors.txt", delimiter=",")
    #: delete :neighbors = np.loadtxt("neighbors.txt")
    #: delete :print(neighbors)
    #: delete :neighbors_clean = neighbors[np.isfinite(neighbors)] #for some
    #: delete :print(neighbors_clean)

def plotNeighbors():
    fig, ax = plt.subplots(figsize=(8,8))

    f = open("neighbors.txt", "r")
    txt = f.read()
    nblist = []
    for line in txt.splitlines():
        nblist.append([int(i) for i in line.split()])
    #print(nblist)

    cidx = np.zeros((len(points), len(points)))

    for i in range(len(nblist)):
        #print("i: ", i)
        for j in range(len(nblist[i])):
            if j > 0:
                #print(" nbs: ", nblist[i][j])
                cidx[i, nblist[i][j]] = 1

    plt.imshow(cidx)
    plt.savefig("neighbors.png")

plotParticles()
plotNeighbors()
