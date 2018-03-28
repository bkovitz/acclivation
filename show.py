#!/usr/bin/env python2
import re
import sys
from tempfile import NamedTemporaryFile
import subprocess

class Organism:
  def __init__(self, fitness, nodes, edges, g, p):
    self.fitness = fitness
    self.nodes = nodes
    self.edges = edges
    self.g = g
    self.p = p
    self.dot = []

def parse(filename):
  organisms = {}
  orgRe = re.compile('^organism \[(\d+),(\d+),(\d+)\]\s+fitness\s*=\s*([0-9\.]+)\s*nodes\s*=\s*(\d+).*')
  curKey = None
  curOrg = None
  curDot = []
  for line in open(filename):
    if line.startswith('organism'):
      m = orgRe.match(line)
      if m:
        epoch, generation, index, fitness, nodes = m.groups()
        curKey = (int(epoch), int(generation), int(index))
        curOrg = Organism(fitness, nodes, 0, 0, 0)
        curDot = []
    elif line.startswith(' ') or line.startswith('digraph'):
      curDot.append(line.rstrip())
    elif line.startswith('}'):
      curDot.append(line.rstrip())
      curOrg.dot = curDot
      organisms[curKey] = curOrg
      curKey = None
      curOrg = None
  return organisms

def show(organisms, epoch, generation, index):
  o = organisms[(epoch, generation, index)]
  outf = 'e%dg%do%d' % (epoch, generation, index)
  with open(outf + '.dot', 'w') as f:
    for line in o.dot:
      f.write("%s\n" % line)
  subprocess.call(['make', outf + '.pdf'])
  subprocess.call(['evince', outf + '.pdf'])

if __name__=='__main__':
  if len(sys.argv) != 4:
    print 'usage: show <epoch> <generation> <index>'
  else:
    filename = 'ancestors'
    epoch = int(sys.argv[1])
    generation = int(sys.argv[2])
    index = int(sys.argv[3])
    organisms = parse(filename)
    show(organisms, epoch, generation, index)
