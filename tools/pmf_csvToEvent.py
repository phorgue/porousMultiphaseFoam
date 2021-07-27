#!/usr/bin/python

import numpy as np
import argparse

# parse argument
parser = argparse.ArgumentParser()
# positional arguments
parser.add_argument("distribution", help='specificy the distribution (x,y,z position and weight)')
parser.add_argument("input", help='specify input file .csv')
parser.add_argument("output", help='specify output file .evt')
# optional arguments
parser.add_argument("-m", help='optional multiplier', type=float, default=1)
parser.add_argument("-c", help='specify csv column (start at 0)', type=int, default=1)
args = parser.parse_args()

surface = args.m
column = args.c

# reading files
distrib = np.genfromtxt(args.distribution)
npoints = np.shape(distrib)[0]

input_csv = np.genfromtxt(args.input)
ntimes = np.shape(input_csv)[0]
time = input_csv[:,0]

data_in = -input_csv[:,args.c]*surface
file_out = open(args.output, 'w')

data_in = data_in.tolist()
distrib_list = distrib.tolist()
time = time.tolist()

iter_time = 0
while iter_time < ntimes: 
    line = "date "+str(time[iter_time])+'\n'
    file_out.write(line)
    iter_points=0
    while iter_points<npoints:
        line = str(distrib_list[1][0])+" "
        line += str(distrib_list[1][1])+" "
        line += str(distrib_list[1][2])+" "
        line += str(data_in[iter_time]*distrib[iter_points,3])+'\n'
        file_out.write(line)
        iter_points += 1
    iter_time += 1 
