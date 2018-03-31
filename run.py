#!/usr/bin/env python2
import itertools
import subprocess
import re
import sys

num_iterations_per_param_set = 20

parse_name_val = r'(?P<name>\w+) *= *(?P<value>[-+]?[a-zA-Z0-9_.]+)'

def main():
    paramss = [
      ['--bumps=0', '--bumps=1'],
      ['--c2=1.0 --c3=0.0', '--c2=2.0 --c3=-0.45'],
      ['--activation_types=0', '--activation_types=1'],
      ['--crossover_freq=0.1', '--crossover_freq=0.2', '--crossover_freq=0.3'],
      ['--ridge_radius=0.05', '--ridge_radius=0.2'],
      ['--num_epochs=20', '--num_epochs=100', '--num_epochs=200'],
    ]
    for row_num, params in enumerate(itertools.product(*paramss)):
        #print params
        flat_params = []
        for param in params:
            for p in re.split('\s+', param):
                flat_params.append(p)
        for i in range(num_iterations_per_param_set):
            output = subprocess.check_output(['./sa', '--quiet'] + flat_params)

            name_vals = {}
            for line in output.splitlines():
                for (name, val) in re.findall(parse_name_val, line):
                    name_vals[name] = val
                    continue
            sorted_names = sorted(name_vals.keys())
            if row_num == 0:
                print ','.join(sorted_names)
            print ','.join([name_vals[n] for n in sorted_names])
            sys.stdout.flush()

if __name__ == '__main__':
    main()
