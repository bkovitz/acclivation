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
  with open(filename, 'rb') as csvfile:
    virt_func = csv.reader(csvfile, delimiter=' ')
    for row in virt_func:
      g1, g2, p1, p2, f = row[:]
      X.append(g1)
      Y.append(g2)
      Z.append(f)
  return X, Y, Z

def plot(X, Y, Z):
  fig = plt.figure()
  ax = Axes3D(fig)
  #X = np.arange(-5, 5, 0.25)
  #Y = np.arange(-5, 5, 0.25)
  X, Y = np.meshgrid(X, Y)
  #R = np.sqrt(X**2 + Y**2)
  #Z = np.sin(R)
  ax.plot_surface(X, Y, Z, alpha=.8, rstride=1, cstride=1, linewidth=1, cmap=cm.jet, shade=True)
  plt.show()
  #plt.savefig('h-dist%s.png'%suffix, dpi=100)

if __name__=='__main__':
  X, Y, Z = read(sys.argv[1])
  plot(X, Y, Z)
