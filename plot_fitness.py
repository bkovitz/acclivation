#!/usr/bin/env python
import csv
import matplotlib.pyplot as plt
import numpy as np
import sys

from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D

def parse_fitness_from_results(filename):
  X = []
  Y = []
  expected_generation = 0
  fitness = None
  last_epoch_num = None
  run_X = []
  run_Y = []
  for line in open(filename):
    line = line.strip()
    if line.startswith('epoch'):
      epoch_num = int(line.split(' ')[1])
      if epoch_num == 0:
        # epoch label preceeds generation data
        if last_epoch_num is not None:
          # start of new run
          X.append(np.array(run_X))
          Y.append(np.array(run_Y))
          run_X = []
          run_Y = []
          expected_generation = 0
        print('start of run', len(X))
        last_epoch_num = 0
        continue
      if fitness == None:
        raise Exception('no fitness for epoch %d' % epoch_num)
      if fitness == -10.0:
        run_Y.append(0.0)
      else:
        run_Y.append(fitness)
      fitness = None
      run_X.append(epoch_num - 1)
      expected_generation = 0
      last_epoch_num = epoch_num
    if line.startswith('generation'):
      generation = int(line.split(' ')[1])
      if generation != expected_generation:
        raise Exception('got generation %d at expected generation %d'
            % (generation, expected_generation))
        return [], []
    if line.startswith('best fitness'):
      fitness = float(line.split(' ')[1].split('=')[-1])
      expected_generation += 1
  if run_X:
    X.append(np.array(run_X))
    Y.append(np.array(run_Y))
  return X, Y

def plot(X, Y):
  for i,x in enumerate(X):
    y = Y[i]
    plt.figure(i, figsize=(10, 2))
    plt.plot(x, y)
    #fig.set_size_inches(10.0, 2)
  plt.show()
  #plt.savefig('%s.png', dpi=100)

if __name__=='__main__':
  filename = sys.argv[1]
  X, Y = parse_fitness_from_results(filename)
  plot(X, Y)
