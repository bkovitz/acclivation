#!/usr/bin/env python2
import re
import sys
import subprocess
# ----------------------------------------------------------------------------

class Birth:
  traverse = {}
  def get_parents(self):
    return []

class Mutation(Birth):
  def __init__(self, parent, child, changes):
    self.parent = parent
    self.child = child
    self.changes = changes
  def get_parents(self):
    return [self.parent]

class Crossover(Birth):
  def __init__(self, mommy, daddy, child):
    self.mommy = mommy
    self.daddy = daddy
    self.child = child
  def get_parents(self):
    return [self.mommy, self.daddy]

class Organism:
  def __init__(self, me, fitness, nodes, edges, g, p, birth):
    self.me = me
    self.fitness = fitness
    self.nodes = nodes
    self.edges = edges
    self.g = g
    self.p = p
    self.birth = birth
    self.dot = []
  def get_lineage(self, organisms, sofar={}, links=[]):
    if sofar is None:
      sofar = {}
    if sofar.has_key(self.me):
      return
    sofar[self.me] = True
    for p in self.birth.get_parents():
      links.append((p,self.me))
    for p in self.birth.get_parents():
      organisms[p].get_lineage(organisms, sofar, links)
    return links

# ----------------------------------------------------------------------------
def parse(filename):
  organisms = {}
  births = {}
  organismRe = re.compile('^organism\s*\[(\d+),(\d+),(\d+)\]\s+' + \
    'fitness\s*=\s*([-0-9\.]+)\s*' + \
    'nodes\s*=\s*(\d+)\s*' + \
    'edges\s*=\s*(\d+)\s*' + \
    'g-vector\s*=\s*\[\s*([-0-9\.]+)\s+([-0-9\.]+)\s*]\s*' + \
    'phenotype\s*=\s*\[\s*([-0-9\.]+)\s+([-0-9\.]+)\s*\]$')
  crossoverRe = re.compile('^crossover\s*' + \
    '(\d+),(\d+),(\d+)\s+' + \
    '(\d+),(\d+),(\d+)\s+' + \
    '(\d+),(\d+),(\d+)$')
  mutationRe = re.compile('^mutation\s*' + \
    '(\d+),(\d+),(\d+)\s+' + \
    '(\d+),(\d+),(\d+)\s+' + \
    '\[\s*([a-z_ ]+)+\s*\]$')
  curKey = None
  curOrg = None
  curDot = []
  for line in open(filename):
    if line.startswith('organism'):
      m = organismRe.match(line)
      if m:
        epoch, generation, index, fitness, nodes, edges, g1, g2, p1, p2 = m.groups()
        curKey = (int(epoch), int(generation), int(index))
        if epoch == '1' and generation == '1':
          curBirth = Birth()
        else:
          curBirth = births[curKey]
        curOrg = Organism(curKey, fitness, nodes, edges, [g1, g2], [p1, p2], curBirth)
        curDot = []
    elif line.startswith(' ') or line.startswith('digraph'):
      curDot.append(line.rstrip())
    elif line.startswith('}'):
      curDot.append(line.rstrip())
      curOrg.dot = curDot
      organisms[curKey] = curOrg
      curKey = None
      curOrg = None
    elif line.startswith('mutation'):
      m = mutationRe.match(line)
      if m:
        pe, pg, pi, ce, cg, ci, changes = m.groups()
        parentKey = (int(pe), int(pg), int(pi))
        childKey = (int(ce), int(cg), int(ci))
        births[childKey] = Mutation(parentKey, childKey, changes.split())
      else:
        raise Exception('bad mutation line: ' + line)
    elif line.startswith('crossover'):
      m = crossoverRe.match(line)
      if m:
        me, mg, mi, de, dg, di, ce, cg, ci = m.groups()
        momKey = (int(me), int(mg), int(mi))
        dadKey = (int(de), int(dg), int(di))
        childKey = (int(ce), int(cg), int(ci))
        births[childKey] = Crossover(momKey, dadKey, childKey)
      else:
        raise Exception('bad crossover line: ' + line)
  return births, organisms

def lineage(births, organisms, epoch, generation, index):
  o = organisms[(epoch, generation, index)]
  outf = 'e%dg%do%d-lineage' % (epoch, generation, index)
  with open(outf + '.dot', 'w') as f:
    f.write('digraph g {\n')
    for parent,child in o.get_lineage(organisms):
      f.write('  %s -> %s\n' %
          ('e%dg%do%d' % parent,
          'e%dg%do%d' % child))
    f.write('}\n')
  subprocess.call(['make', outf + '.pdf'])
  subprocess.call(['evince', outf + '.pdf'])

def show(births, organisms, epoch, generation, index):
  o = organisms[(epoch, generation, index)]
  outf = 'e%dg%do%d' % (epoch, generation, index)
  with open(outf + '.dot', 'w') as f:
    for line in o.dot:
      f.write("%s\n" % line)
  subprocess.call(['make', outf + '.pdf'])
  subprocess.call(['evince', outf + '.pdf'])

# ----------------------------------------------------------------------------
if __name__=='__main__':
  if len(sys.argv) not in (4, 5):
    print 'usage: show <epoch> <generation> <index> [lineage]'
  else:
    filename = 'ancestors'
    epoch = int(sys.argv[1])
    generation = int(sys.argv[2])
    index = int(sys.argv[3])
    births, organisms = parse(filename)
    if 'lineage' in sys.argv:
      lineage(births, organisms, epoch, generation, index)
    else:
      show(births, organisms, epoch, generation, index)
