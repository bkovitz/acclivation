#!/usr/bin/env python
import csv
import matplotlib.pyplot as plt
import numpy as np
import sys

from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D

def read(filename):
  X = []
  Y = []
  Z = []
  with open(filename, 'r') as csvfile:
    virt_func = csv.reader(csvfile, delimiter=' ')
    for row in virt_func:
      g1, g2, p1, p2, f = row[:]
      X.append(float(g1))
      Y.append(float(g2))
      Z.append(float(f))
  return np.array(X), np.array(Y), np.array(Z)

def plot(X, Y, Z):
  fig = plt.figure()
  ax = Axes3D(fig)
  ax.plot_trisurf(X, Y, Z, cmap=cm.jet, edgecolor='none', shade=True)
  plt.show()
  #plt.savefig('%s.png', dpi=100)

if __name__=='__main__':
  X, Y, Z = read(sys.argv[1])
  plot(X, Y, Z)
