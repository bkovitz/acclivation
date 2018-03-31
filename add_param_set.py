#!/usr/bin/env python2
import csv
import sys

def main(csvfile, num_iterations_per_param_set):
    n = 1
    param_set = 1
    first = True
    with open(csvfile, 'rb') as csvfile:
        reader = csv.reader(csvfile, delimiter=',')
        for row in reader:
            if first:
                print ','.join(['param_set'] + row)
                first = False
            else:
                print ','.join(['set%03d' % param_set] + row)
                if (n % num_iterations_per_param_set) == 0:
                    param_set += 1
                n += 1

if __name__ == '__main__':
    csvfile = sys.argv[1]
    if len(sys.argv) > 2:
        num_iterations_per_param_set = int(sys.argv[2])
    else:
        num_iterations_per_param_set = 20
    main(csvfile, num_iterations_per_param_set)
