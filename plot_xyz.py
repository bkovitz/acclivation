#!/usr/bin/env python2
import csv
import matplotlib.pyplot as plt
import numpy as np
import sys
from math import sqrt

from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D

def parse(lines, x_col, y_col, z_col):
  X = []
  Y = []
  Z = []
  virt_func = csv.reader(lines.splitlines(), delimiter=' ')
  for row in virt_func:
    X.append(float(row[x_col]))
    Y.append(float(row[y_col]))
    Z.append(float(row[z_col]))
  return np.array(X), np.array(Y), np.array(Z)

def read(filename, x_col, y_col, z_col):
  with open(filename, 'r') as csvfile:
    return parse(csvfile.read(), x_col, y_col, z_col)

def plot(X, Y, Z, scatter, show=True, filename='', ylim1=-1, ylim2=1, xlim1=-1,
        xlim2=1, zlim1=0, zlim2=11, cmapname='coolwarm', linewidth=0, elev=0.0,
        azim=180.0):
  fig = plt.figure()
  ax = Axes3D(fig)
  if scatter:
    ax.scatter3D(X, Y, Z)
    ax.set_ylim((-1,1))
    ax.set_xlim((-1,1))
    ax.set_zlim((0,11))
  else:
    dim = int(sqrt(len(X)))
    X = X.reshape((dim, dim))
    Y = Y.reshape((dim, dim))
    Z = Z.reshape((dim, dim))
    #ax.plot_surface(X, Y, Z, rcount=dim, ccount=dim, cmap=cm.coolwarm, linewidth=3)
    #ax.plot_surface(X, Y, Z, rstride=2, cstride=2, cmap=cm.coolwarm, linewidth=0)
    stride = 2 if show else 1
    #print(cm)
    #print(cm.__dict__)
    cmap = getattr(cm, cmapname)
    ax.view_init(elev, azim)
    ax.plot_surface(X, Y, Z, rstride=stride, cstride=stride, cmap=cmap, linewidth=linewidth)
    #ax.plot_surface(X, Y, Z, rstride=stride, cstride=stride, cmap=cm.coolwarm, linewidth=0)
    #ax.set_ylim((-1,1))
    #ax.set_xlim((-1,1))
    #ax.set_zlim((0,11))
    ax.set_ylim((ylim1,ylim2))
    ax.set_xlim((xlim1,xlim2))
    ax.set_zlim((zlim1,zlim2))
    #ax.plot_trisurf(X, Y, Z, cmap=cm.jet, edgecolor='none', shade=True)
  if show:
    plt.show()
  else:
    if filename:
        plt.savefig(filename, dpi=300)

if __name__=='__main__':
  filename = sys.argv[1]
  x_col = int(sys.argv[2])
  y_col = int(sys.argv[3])
  z_col = int(sys.argv[4])
  scatter = (sys.argv[5] == 'scatter') if len(sys.argv) > 5 else False
  X, Y, Z = read(filename, x_col, y_col, z_col)
  plot(X, Y, Z, scatter)
