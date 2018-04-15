#!/usr/bin/env python2
import atexit
import re
import readline
import sa
import sys
import os

from cStringIO import StringIO
from graphviz import Source
from Tkinter import Tk, Frame, Scrollbar, Canvas, mainloop, LAST, ALL, \
    HORIZONTAL, VERTICAL, BOTTOM, TOP, X, Y, BOTH, LEFT, RIGHT, YES

# ----------------------------------------------------------------------------

class SASim(object):
    # - maybe add some text tag to edges?
    # - maybe click to see history of a node's activations? outputs?
    # - any (easy) way to make arrows reach 'to' node?
    root = None
    def __init__(self, graphData, numSteps, simResults):
        self.scale = 90.
        self.margin = 40
        self.halfMargin = self.margin / 2.
        self.parseGraphData(graphData)
        self.buildDisplay()
        self.parseSimResults(numSteps, simResults)
        self.update()
        self.root.mainloop()

    def parseSimResults(self, numSteps, simResults):
        steps = []
        for line in simResults.split('\n'):
            if line.startswith('initial') or line.startswith('final'):
                activations = [float(a) for a in re.findall('(?P<activation>[0-9.-]+)', line)]
                assert(len(activations) == len(self.nodes))
                steps.append(activations)
        assert(len(steps) == 1 + numSteps) # +1 for initial activation
        self.steps = steps

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
                m = re.match('^edge ([a-zA-Z0-9]+) ([a-zA-Z0-9]+) ([0-9]+) ([0-9. ]+)', line)
                src, dst, n, coordsStr = m.groups()
                n = int(n)
                coords = [float(c) for c in coordsStr.strip().split()[:2*n]]
                coords = [self.invertY(c) if i&1 else c for i,c in enumerate(coords)]
                edges.append([src, dst, n, coords])
        self.nodes = nodes
        self.edges = edges

    def update(self):
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

    def left(self, event):
        if self.curStep > 0:
            self.curStep -= 1
            print self.curStep
        self.update()

    def right(self, event):
        if self.curStep < len(self.steps)-1:
            self.curStep += 1
            print self.curStep
        self.update()

    def buildDisplay(self):
        if self.root is None:
            self.root = Tk()
        winWidth = 1000
        winHeight = 400
        frame = Frame(self.root, width=winWidth, height=winHeight)
        frame.pack(expand=YES, fill=BOTH)

        canvasWidth = int(self.graphWidth * self.scale + self.margin)
        canvasHeight = int(self.graphHeight * self.scale + self.margin)
        c = Canvas(frame, width=winWidth, height=winHeight, scrollregion=(0, 0, canvasWidth, canvasHeight))

        hbar = Scrollbar(frame, orient=HORIZONTAL)
        hbar.pack(side=BOTTOM, fill=X)
        hbar.config(command=c.xview)
        vbar = Scrollbar(frame, orient=VERTICAL)
        vbar.pack(side=RIGHT, fill=Y)
        vbar.config(command=c.yview)

        c.config(xscrollcommand=hbar.set, yscrollcommand=vbar.set)
        c.pack(side=LEFT, expand=YES, fill=BOTH)

        nodeIds = []
        m = self.halfMargin
        for name, x, y, w, h, label in self.nodes:
            c.create_oval(m + (x-w/2.)*self.scale, m + (y-h/2.)*self.scale, m + (x+w/2.)*self.scale, m + (y+h/2.)*self.scale, width=3)
            nodeIds.append(c.create_text(m + x*self.scale, m + y*self.scale, text=label))
        self.nodeIds = nodeIds
        for src, dst, n, coords in self.edges:
            scaledCoords = [m + x*self.scale for x in coords]
            c.create_line(*scaledCoords, smooth=1, width=2, arrow=LAST, arrowshape=(15,20,10))
        c.bind_all('<Left>', self.left)
        c.bind_all('<Right>', self.right)
        #c.config(scrollregion=c.bbox(ALL))
        self.canvas = c
        self.curStep = 0

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


# ----------------------------------------------------------------------------


class Runner(object):
    def __init__(self, filename='ancestors'):
        self.runData = sa.load_ancestor_file(filename)
        self.world = self.runData.w
        self.buildOrgMap()
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
            'print': Command('print current organism info', {}),
            'parent': Command('select parent of current organism', {}),
            'child': Command('select child of current organism', {}),
            'dot': Command('print current organism dot', {}),
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
        for e, epoch in enumerate(self.runData):
            for g, generation in enumerate(epoch):
                for o, organism in enumerate(generation):
                    orgMap[(e+1, g+1, o)] = organism
        self.orgMap = orgMap

    def getWorldParams(self):
        buf = StringIO()
        sa.print_world_params(self.world, buf)
        worldParamsStr = buf.getvalue()
        worldParams = set()
        for line in worldParamsStr.split('\n'):
            if line.startswith('w->'):
                param = line[len('w->'):].split('=')[0]
                worldParams.add('w.' + param + '=')
        return worldParams

    def cmdQuit(self):
        sys.exit(0)

    def cmdSelect(self, selection):
        if selection in self.orgMap:
            self.selectedOrg = self.orgMap[selection]
            self.selected = selection
        else:
            print 'org?'

    def cmdSa(self):
        w = self.world
        o = self.selectedOrg
        sa.sa(o, w.sa_timesteps, w.decay, w.spreading_rate)
        print sa.phenotype_fitness(w, o.genotype) # can't use phenotype_fitness_func yet

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
        SASim(graphData, self.world.sa_timesteps, simResults)

    def cmdPrint(self):
        o = self.selectedOrg
        print self.selectedStr()
        print 'fitness =', o.fitness
        print 'from_turned_knob =', o.from_turned_knob
        print 'num_nodes =', o.genotype.num_nodes
        print 'num_edges =', o.genotype.num_edges

    def cmdParent(self):
        pass

    def cmdChild(self):
        pass

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
