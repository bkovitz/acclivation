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

#see https://matplotlib.org/tutorials/colors/colormaps.html
#
#binary - little too white
#ocean - nice green/blue
#gist_stern - shows nice detail on the many bumps
#rainbow - slightly better contrast within bumps than default
#gnuplot2_r - white/yellow bumps stand out
#gnuplot_r - everything clear
#
def plot(X, Y, Z, scatter, show=True, filename='', ylim1=-1, ylim2=1, xlim1=-1,
        xlim2=1, zlim1=0, zlim2=11, linewidth=0, elev=33.0, cmapname='coolwarm', 
        azim=33.0, dpi=300, stride=1, pointsize=20, wtitle=''):
  cmap = plt.get_cmap(cmapname)
  fig = plt.figure(num=wtitle)
  ax = Axes3D(fig)
  if scatter:
    #ax.scatter3D(X, Y, Z)
    #ax.scatter3D(X, Y, Z, s=20, c=np.maximum(X+1.0, Y+1.0)/2.0, cmap='plasma', #perceptually uniform fade from origin to 1,1
    ax.scatter3D(X, Y, Z, s=pointsize, c=np.maximum(X+1.0, Y+1.0)/2.0, cmap=cmap,
            edgecolors='none')
  else:
    dim = int(sqrt(len(X)))
    print('dim : %dx%d' % (dim, dim))
    X = X.reshape((dim, dim))
    Y = Y.reshape((dim, dim))
    Z = Z.reshape((dim, dim))
    ax.view_init(elev, azim)
    ax.plot_surface(X, Y, Z, rstride=stride, cstride=stride, cmap=cmap, linewidth=linewidth)
    #ax.plot_trisurf(X, Y, Z, cmap=cm.jet, edgecolor='none', shade=True)
  ax.set_ylim((ylim1,ylim2))
  ax.set_xlim((xlim1,xlim2))
  ax.set_zlim((zlim1,zlim2))
  if show:
    plt.show()
  else:
    if filename:
        print('saving to %s' % filename)
        plt.savefig(filename, dpi=dpi)

if __name__=='__main__':
  filename = sys.argv[1]
  x_col = int(sys.argv[2])
  y_col = int(sys.argv[3])
  z_col = int(sys.argv[4])
  scatter = (sys.argv[5] == 'scatter') if len(sys.argv) > 5 else False
  X, Y, Z = read(filename, x_col, y_col, z_col)
  plot(X, Y, Z, scatter)
