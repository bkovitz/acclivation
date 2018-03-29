#!/usr/bin/env python2
import csv
import matplotlib.pyplot as plt
import numpy as np
import sys
from math import sqrt

from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D

def read(filename, x_col, y_col, z_col):
  X = []
  Y = []
  Z = []
  with open(filename, 'r') as csvfile:
    virt_func = csv.reader(csvfile, delimiter=' ')
    for row in virt_func:
      X.append(float(row[x_col]))
      Y.append(float(row[y_col]))
      Z.append(float(row[z_col]))
  return np.array(X), np.array(Y), np.array(Z)

def plot(X, Y, Z, scatter):
  fig = plt.figure()
  ax = Axes3D(fig)
  if scatter:
    plt.ylim((-1,1))
    plt.xlim((-1,1))
    ax.scatter3D(X, Y, Z)
  else:
    dim = int(sqrt(len(X)))
    X = X.reshape((dim, dim))
    Y = Y.reshape((dim, dim))
    Z = Z.reshape((dim, dim))
    ax.plot_surface(X, Y, Z, rcount=dim, ccount=dim, cmap=cm.coolwarm, linewidth=3)
    #ax.plot_trisurf(X, Y, Z, cmap=cm.jet, edgecolor='none', shade=True)
  plt.show()
  #plt.savefig('%s.png', dpi=100)

if __name__=='__main__':
  filename = sys.argv[1]
  x_col = int(sys.argv[2])
  y_col = int(sys.argv[3])
  z_col = int(sys.argv[4])
  scatter = (sys.argv[5] == 'scatter') if len(sys.argv) > 5 else False
  X, Y, Z = read(filename, x_col, y_col, z_col)
  plot(X, Y, Z, scatter)
