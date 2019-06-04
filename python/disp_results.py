import numpy as np
import csv
import os, glob
from numpy.random import normal as normal
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.animation as animation
import matplotlib


def update(ifrm, xa, ya, za):
    sct.set_data(xa[ifrm], ya[ifrm])
    sct.set_3d_properties(za[ifrm])


def read_file(pathfile):
    files = glob.glob(pathfile + "*.csv")
    files.sort(key=os.path.getmtime)

    xs = [[] for i in range(len(files))]
    ys = [[] for i in range(len(files))]
    zs = [[] for i in range(len(files))]

    for idx, f in enumerate(files):
        with open(pathfile + f, newline='') as data:
            reader = csv.reader(data)

            for row in reader:
                if row[0] != 'x':
                    xs[idx].append(float(row[0]))
                    ys[idx].append(float(row[1]))
                    zs[idx].append(float(row[2]))

    return xs, ys, zs


if __name__ == "__main__":
    pathfile = "../output/"
    dirs = os.listdir(pathfile)

    nfr = len(dirs)  # Number of frames
    fps = 90  # Frame per sec
    side = 200

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    sct, = ax.plot([], [], [], "o", markersize=2)

    # Read file
    xs, ys, zs = read_file(pathfile)

    ax.set_xlim(-side, side)
    ax.set_ylim(-side, side)
    ax.set_zlim(-side, side)
    ani = animation.FuncAnimation(fig, update, nfr, fargs=(xs, ys, zs), interval=1000 / fps)

    plt.rcParams['animation.html'] = 'html5'
    plt.show()