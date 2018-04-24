#!/usr/bin/env python2
import atexit
import os
import re
import readline
import sa
import sys

from cStringIO import StringIO
from graphviz import Source
from Queue import Queue
from threading import Thread
from Tkinter import Tk, Frame, Scrollbar, Canvas, mainloop, LAST, ALL, \
    HORIZONTAL, VERTICAL, BOTTOM, TOP, X, Y, BOTH, LEFT, RIGHT, YES, NW

# ----------------------------------------------------------------------------

def parseSimResults(simResults, numSteps, numNodes):
    steps = []
    for line in simResults.split('\n'):
        if line.startswith('initial') or line.startswith('final'):
            activations = [float(a) for a in re.findall('(?P<activation>[0-9.-]+)', line)]
            assert(len(activations) == numNodes)
            steps.append(activations)
    assert(len(steps) == 1 + numSteps) # +1 for initial activation
    return steps


class SASim(Thread):
    # - maybe add some text tag to edges?
    # - maybe click to see history of a node's activations? outputs?
    # - simple way of editing graph
    #   - add/rem node
    #   - add/rem/move edge
    #   - turn knob/control
    #   - alter activation type, input acc, output type
    # - any (easy) way to make arrows reach 'to' node?
    def __init__(self, graphData, numSteps, simResults):
        Thread.__init__(self)
        self.q = Queue()
        self.scale = 90.
        self.margin = 40
        self.canvas = None
        self.halfMargin = self.margin / 2.
        self.graphData = graphData
        self.numSteps = numSteps
        self.simResults = simResults
        self.start()

    def onUi(self, func):
        self.q.put(func)

    def quit(self):
        self.root.destroy()

    def pollLoop(self):
        while True:
            try:
                func = self.q.get(block=False)
            except:
                break
            else:
                self.root.after_idle(func)
        self.root.after(100, self.pollLoop)

    def run(self):
        self.root = Tk()
        self.root.protocol('WM_DELETE_WINDOW', self.quit)
        self.setup(self.graphData, self.numSteps, self.simResults)
        self.root.after(100, self.pollLoop)
        self.root.mainloop()

    def update(self, graphData, numSteps, simResults):
        def do_update():
            self.setup(graphData, numSteps, simResults)
        self.onUi(do_update)

    def setup(self, graphData, numSteps, simResults):
        self.parseGraphData(graphData)
        if self.canvas is None:
            self.buildDisplay()
        else:
            self.canvas.configure(scrollregion=(0, 0, self.canvasWidth, self.canvasHeight))
        self.populateCanvas()
        self.steps = parseSimResults(simResults, numSteps, len(self.nodes))
        self.updateActivations()

    def invertY(self, y):
        return self.graphHeight - y

    def parseGraphData(self, graphData):
        nodes = []
        edges = []
        for line in graphData.split('\n'):
            if line.startswith('graph'):
                m = re.match('^graph (\d+) ([0-9.]+) ([0-9.]+)$', line)
                _id, self.graphWidth, self.graphHeight = [float(x) for x in m.groups()]
            elif line.startswith('node'):
                m = re.match('^node ([a-zA-Z0-9]+) ([0-9.]+) ([0-9.]+) ([0-9.]+) ([0-9.]+) ("[^"]+")', line)
                name, x, y, width, height, label = m.groups()
                x, y, width, height = [float(x) for x in (x, y, width, height)]
                y = self.invertY(y)
                nodes.append([name, x, y, width, height, label.strip('"')])
            elif line.startswith('edge'):
                if 'invis' in line:
                    continue
                m = re.match('^edge ([a-zA-Z0-9]+) ([a-zA-Z0-9]+) ([0-9]+) ([0-9. ]+)', line)
                src, dst, n, coordsStr = m.groups()
                n = int(n)
                coords = [float(c) for c in coordsStr.strip().split()[:2*n]]
                coords = [self.invertY(c) if i&1 else c for i,c in enumerate(coords)]
                edges.append([src, dst, n, coords])
        self.nodes = nodes
        self.edges = edges
        self.canvasWidth = int(self.graphWidth * self.scale + self.margin)
        self.canvasHeight = int(self.graphHeight * self.scale + self.margin)

    def updateActivations(self):
        for n,_id in enumerate(self.nodeIds):
            fullText = self.nodes[n][-1]
            act = self.steps[self.curStep][n]
            actStr = ('%5.4f' % act) if act != -1000.0 else '-'
            m = re.search('^([a-zA-Z0-9]+) (\([^\)]+\)) (i=[0-9.-]+)?\s?([0-9.-]+)?\s?(c=[0-9.-]+)?', fullText)
            if m:
                nm, typ, i, initialAct, c = m.groups()
                s = [nm]
                if typ: s.append(typ)
                if i: s.append(i)
                s.append(actStr)
                if c: s.append(c)
                self.canvas.itemconfigure(_id, text=' '.join(s))
            else:
                # unrecognized format: just use activation as entire label
                self.canvas.itemconfigure(_id, text=actStr)

    def updateStepText(self):
        self.canvas.itemconfigure(self.stepId, text='[%d]' % self.curStep)

    def left(self, event):
        if self.curStep > 0:
            self.curStep -= 1
            self.updateStepText()
        self.updateActivations()

    def right(self, event):
        if self.curStep < len(self.steps)-1:
            self.curStep += 1
            self.updateStepText()
        self.updateActivations()

    def buildDisplay(self):
        winWidth = 1000
        winHeight = 400
        frame = Frame(self.root, width=winWidth, height=winHeight)
        frame.pack(expand=YES, fill=BOTH)

        c = Canvas(frame, width=winWidth, height=winHeight, scrollregion=(0, 0, self.canvasWidth, self.canvasHeight))

        hbar = Scrollbar(frame, orient=HORIZONTAL)
        hbar.pack(side=BOTTOM, fill=X)
        hbar.config(command=c.xview)
        vbar = Scrollbar(frame, orient=VERTICAL)
        vbar.pack(side=RIGHT, fill=Y)
        vbar.config(command=c.yview)

        c.config(xscrollcommand=hbar.set, yscrollcommand=vbar.set)
        c.pack(side=LEFT, expand=YES, fill=BOTH)
        self.canvas = c

    def populateCanvas(self):
        self.canvas.delete('all')
        nodeIds = []
        m = self.halfMargin
        for name, x, y, w, h, label in self.nodes:
            self.canvas.create_oval(m + (x-w/2.)*self.scale, m + (y-h/2.)*self.scale, m + (x+w/2.)*self.scale, m + (y+h/2.)*self.scale, width=3)
            nodeIds.append(self.canvas.create_text(m + x*self.scale, m + y*self.scale, text=label))
        self.nodeIds = nodeIds

        for src, dst, n, coords in self.edges:
            scaledCoords = [m + x*self.scale for x in coords]
            self.canvas.create_line(*scaledCoords, smooth=1, width=2, arrow=LAST, arrowshape=(15,20,10))

        self.canvas.bind_all('<Left>', self.left)
        self.canvas.bind_all('<Right>', self.right)
        self.curStep = 0
        self.stepId = self.canvas.create_text(3, 3, text='[%d]' % self.curStep, anchor=NW)

# ----------------------------------------------------------------------------

class Command(object):
    def __init__(self, help, choices):
        self.help = help
        self.choices = choices


class StringArg(object):
    def __init__(self, s):
        self.s = s

    def tabMatch(self, s):
        if self.s.startswith(s):
            return [self.s + ' ']
        else:
            return []

    def match(self, s):
        return self.s.startswith(s)

    def get(self, s):
        return s


class RegexArg(object):
    def __init__(self, regex, fmt):
        self.regex = regex
        self.fmt = fmt

    def tabMatch(self, s):
        return [(s if re.match(self.regex, s) else self.fmt) + ' ']

    def match(self, s):
        return re.match(self.regex, s) is not None

    def get(self, s):
        return re.match(self.regex, s).groups()


class SelectByEpochGenerationOrganismArg(RegexArg):
    def __init__(self):
        super(SelectByEpochGenerationOrganismArg, self).__init__(
            '^([0-9]+)\.([0-9]+)\.([0-9]+)$', 'e.g.o')

    def get(self, s):
        return tuple(int(x) for x in re.match(self.regex, s).groups())


class SelectByEpochArg(RegexArg):
    def __init__(self):
        super(SelectByEpochArg, self).__init__('^([0-9]+)\.\.$', 'e..')

    def get(self, s):
        epoch = int(re.match(self.regex, s).group(1))
        return (epoch, 1, 0)


class SelectByGenerationArg(RegexArg):
    def __init__(self, runner):
        self.runner = runner
        super(SelectByGenerationArg, self).__init__('^\.([0-9]+)\.$', '.g.')

    def get(self, s):
        generation = int(re.match(self.regex, s).group(1))
        epoch, g, o = self.runner.selected
        return (epoch, generation, 0)


class SelectByOrganismArg(RegexArg):
    def __init__(self, runner):
        self.runner = runner
        super(SelectByOrganismArg, self).__init__('^\.\.([0-9]+)$', '..o')

    def get(self, s):
        organism = int(re.match(self.regex, s).group(1))
        epoch, generation, o = self.runner.selected
        return (epoch, generation, organism)


class LineageArg(object):
    def __init__(self, runner, relMap):
        self.runner = runner
        self.relMap = relMap

    def toStr(self, rel):
        return '%d.%d.%d' % rel

    def tabMatch(self, s):
        cur = self.runner.selected
        if cur is None:
            return []
        if not cur in self.relMap:
            return []
        relStrs = [self.toStr(rel) for rel in self.relMap[cur]]
        return [relStr + ' ' for relStr in relStrs if relStr.startswith(s)]

    def match(self, s):
        return self.tabMatch(s) != []

    def get(self, s):
        s = self.tabMatch(s)[0].strip()
        return tuple(int(x) for x in re.match('^([0-9]+)\.([0-9]+)\.([0-9]+)$', s).groups())


# ----------------------------------------------------------------------------


class Runner(object):
    def __init__(self, filename='ancestors'):
        sys.stderr.write('loading..')
        self.runData = sa.load_ancestor_file(filename)
        sys.stderr.write('done\n')
        self.world = self.runData.w
        sys.stderr.write('mapping..')
        self.buildOrgMap()
        sys.stderr.write('done\n')
        self.sasim = None
        self.selected = (1, 1, 0)
        self.selectedOrg = self.orgMap[self.selected]
        self.commands = {
            'quit': Command('exit', {}),
            'select':
                Command('select an organism by <epoch>.<generation>.<organism-index>', {
                    SelectByEpochGenerationOrganismArg(),
                    SelectByEpochArg(),
                    SelectByGenerationArg(self),
                    SelectByOrganismArg(self),
                }),
            'sa': Command('run sa on current organism with current world params', {}),
            'saplot': Command('plot node activations by timestep', {}),
            'print': Command('print current organism info', {}),
            'printc': Command('print c code for current organism', {}),
            'parent': Command('select parent of current organism', {
                LineageArg(self, self.parentMap) }),
            'child': Command('select child of current organism', {
                LineageArg(self, self.childMap) }),
            'dot': Command('show graph of current organism dot', {}),
            'sim': Command('simulate current organism', {}),
            'plot':
                Command('plot current phenotype fitness or virtual fitness', {
                    StringArg('pfitness'), StringArg('vfitness')
                }),
            'neighborhood': Command('show current organism\'s fitness neighborhood', {}),
            'list': Command('list parent(s) or child(ren) of current organism', {
                StringArg('parents'), StringArg('children')
            }),
            'world': Command('show current world parameters', {}),
            'help': Command('show help', {}),
        }
        readline.parse_and_bind('tab: complete')
        self.loadHistory()
        readline.set_completer(self.makeCompleter(self.commands,
            self.getWorldParams()))

    @staticmethod
    def loadHistory():
        histfile = os.path.join(os.path.expanduser("~"), ".sa_history")
        try:
            readline.read_history_file(histfile)
            readline.set_history_length(500)
        except IOError:
            pass
        atexit.register(readline.write_history_file, histfile)

    @staticmethod
    def makeCompleter(commands, worldParams):
        def customCompleter(text, state):
            linebuf = readline.get_line_buffer()
            parts = linebuf.split()
            if len(parts) >= 1 and linebuf.endswith(' '):
                # start of new subcommand
                parts.append('')
            if len(parts) < 2:
                matches = [w + ' ' for w in commands.keys()
                           if w.startswith(text)] \
                        + [p for p in worldParams
                           if p.startswith(text)] \
                        + [None]
            else:
                command = parts[0]
                matches = []
                for choice in commands[command].choices:
                    try:
                        matches += choice.tabMatch(parts[1])
                    except Exception, e:
                        print e
                matches += [None]
            return matches[state]
        return customCompleter

    def buildOrgMap(self):
        orgMap = {}
        parentMap = {}
        childMap = {}
        for e, epoch in enumerate(self.runData):
            for g, generation in enumerate(epoch):
                for o, organism in enumerate(generation):
                    # organism map
                    orgId = (e+1, g+1, o)
                    orgMap[orgId] = organism
                    # parent/child maps
                    parents = parentMap.setdefault(orgId, [])
                    if organism.birth_info.type == sa.CROSSOVER:
                        m = organism.birth_info.crossover_info.mom
                        d = organism.birth_info.crossover_info.dad
                        momId = (m.epoch, m.generation, m.org_index)
                        dadId = (d.epoch, d.generation, d.org_index)
                        parents.append(momId)
                        parents.append(dadId)
                        children = childMap.setdefault(momId, [])
                        children.append(orgId)
                        childMap[momId] = children
                        children = childMap.setdefault(dadId, [])
                        children.append(orgId)
                        childMap[dadId] = children
                    elif organism.birth_info.type == sa.MUTATION:
                        p = organism.birth_info.mutation_info.parent
                        parentId = (p.epoch, p.generation, p.org_index)
                        parents.append(parentId)
                        children = childMap.setdefault(parentId, [])
                        children.append(orgId)
                        childMap[parentId] = children
                    parentMap[orgId] = parents
        self.orgMap = orgMap
        self.parentMap = parentMap
        self.childMap = childMap

    def getWorldParams(self):
        buf = StringIO()
        sa.print_world_params(self.world, buf)
        worldParamsStr = buf.getvalue()
        worldParams = set()
        parseNameVal = r'(?P<name>\w+) *= *(?P<value>[-+]?[a-zA-Z0-9_.]+)'
        for (name, val) in re.findall(parseNameVal, worldParamsStr):
            worldParams.add('w.' + name + '=')
        return worldParams

    def cmdQuit(self):
        if self.sasim is not None:
            self.sasim.onUi(lambda: self.sasim.quit())
            self.sasim.join()
        sys.exit(0)

    def cmdSelect(self, selection):
        if selection in self.orgMap:
            self.selectedOrg = self.orgMap[selection]
            self.selected = selection
            self.cmdPrint()
        else:
            w = self.world
            print 'org? <1-%d>.<1-%d>.<0-%d>' % \
                    (w.num_epochs, w.generations_per_epoch, w.num_organisms-1)

    def cmdSa(self):
        w = self.world
        o = self.selectedOrg
        buf = StringIO()
        sa.sa(w, o.genotype, buf)
        print sa.genotype_fitness(w, o.genotype) # can't use phenotype_fitness_func yet

    def cmdSim(self):
        # graph data
        buf = StringIO()
        org = self.selectedOrg
        sa.print_dot(self.world, org, buf)
        dot = buf.getvalue()
        g = Source(dot)
        graphData = g.pipe('plain')
        # sa results
        buf = StringIO()
        sa.sa(self.world, org.genotype, buf)
        simResults = buf.getvalue()
        if self.sasim is None or not self.sasim.is_alive():
            self.sasim = SASim(graphData, self.world.sa_timesteps, simResults)
        else:
            self.sasim.update(graphData, self.world.sa_timesteps, simResults)

    def cmdSaplot(self):
        import numpy as np
        import matplotlib.pyplot as plt
        buf = StringIO()
        w = self.world
        org = self.selectedOrg
        sa.sa(w, org.genotype, buf)
        simResults = buf.getvalue()
        steps = parseSimResults(simResults, w.sa_timesteps, org.genotype.num_nodes_in_use)
        activationsByStep = np.array([np.array(step) for step in steps])
        plt.figure(1, figsize=(10,2))
        l = np.array(['n%d' % n for n in range(org.genotype.num_nodes_in_use)])
        plt.plot(activationsByStep)
        plt.ylim((-1.2,1.2))
        plt.legend(l)
        plt.show()

    def cmdPrint(self):
        print self.selectedStr()
        sa.print_organism(self.selectedOrg)

    def cmdPrintc(self):
        buf = StringIO()
        sa.print_genotype_c(self.selectedOrg.genotype, buf)
        print buf.getvalue()

    def cmdParent(self, parent):
        if parent in self.orgMap:
            self.selectedOrg = self.orgMap[parent]
            self.selected = parent
            self.cmdPrint()
        else:
            print parent, 'not in organism map'

    def cmdChild(self, child):
        if child in self.orgMap:
            self.selectedOrg = self.orgMap[child]
            self.selected = child
            self.cmdPrint()
        else:
            print child, 'not in organism map'

    def cmdPlot(self, typ):
        # delay import because it takes several seconds
        from plot_xyz import plot, parse
        if typ == 'pfitness':
            buf = StringIO()
            sa.dump_phenotype_fitness_func(self.world, False, buf)
            csv = buf.getvalue()
            X, Y, Z = parse(csv, 0, 1, 2)
        elif typ == 'vfitness':
            buf = StringIO()
            sa.dump_organism_virtual_fitness_func(self.world, self.selectedOrg, False, buf)
            csv = buf.getvalue()
            X, Y, Z = parse(csv, 0, 1, 4)
        plot(X, Y, Z, False)

    def cmdDot(self):
        buf = StringIO()
        sa.print_dot(self.world, self.selectedOrg, buf)
        dot = buf.getvalue()
        g = Source(dot)
        g.render(view=True)

    def cmdNeighborhood(self):
        sa.dump_organism_fitness_nbhd(self.world, self.selectedOrg)

    def cmdList(self, relative):
        pass

    def cmdWorld(self):
        sa.print_world_params(self.world, sys.stdout)

    def cmdSetWorld(self, cmd):
        try:
            w = self.world
            exec(cmd)
        except Exception, e:
            print e

    def cmdHelp(self):
        for name, command in self.commands.items():
            self.printHelp(name, command.help)
        self.printHelp('w.<param name>=<value>', 'change world parameter')

    def printHelp(self, command, help):
        print '%s\t- %s' % (command, help)

    def selectedStr(self):
        return '.'.join(['%d' % x for x in self.selected])

    def run(self):
        try:
            while True:
                s = raw_input('>> ').strip()
                tokens = s.split()
                if not len(tokens):
                    print self.selectedStr()
                    continue
                command = tokens[0]
                if command in self.commands.keys():
                    argChoices = self.commands[command].choices
                    if not len(argChoices):
                        # command with no args
                        if len(tokens) != 1:
                            self.printHelp(
                                command, self.commands[command].help)
                        else:
                            getattr(self, 'cmd' + command.title())()
                    else:
                        if len(tokens) != 2:
                            self.printHelp(
                                command, self.commands[command].help)
                        else:
                            arg = tokens[1]
                            gotit = False
                            for argMatch in argChoices:
                                if argMatch.match(arg):
                                    gotit = True
                                    getattr(self, 'cmd' +
                                            command.title())(argMatch.get(arg))
                                    break
                            if not gotit:
                                self.printHelp(
                                    command, self.commands[command].help)
                elif command.startswith('w.'):
                    self.cmdSetWorld(command)
                else:
                    print '?'
        except (EOFError, KeyboardInterrupt) as e:
            print '\n'


# ----------------------------------------------------------------------------
if __name__ == '__main__':
    Runner().run()
    sys.exit(0)