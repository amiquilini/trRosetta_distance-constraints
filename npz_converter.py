#!/software/python/anaconda37/bin/python3
"""
Dalton Lab - UNICAMP
Created on Thu Feb 20 2020
@author: amanda/fabio

Read trRosetta distance distribution matrix and print/save files or graphs
"""

import os
import math
import argparse
import numpy as np
import matplotlib.pyplot as plt


# Calculate potential from Probability
def potential(prob):
    x_values = np.arange(2.,20.,0.5).tolist()
    y_values = []
    for dist, prob_dist in enumerate(prob):
        y_values.append(-math.log(prob_dist)+(math.log((math.pow((x_values[dist]/x_values[-1]), 1.57))*prob[-1])))
        y_values_norm=[(((y - min(y_values))/(max(y_values) - min(y_values) + 1e-20))-1.) for y in y_values] # normalize to values between 0 and -1
    return y_values_norm

def read_fasta(fasta): # Parse fasta file
    name, seq = None, []
    for line in fasta:
        line = line.rstrip()
        if line.startswith(">"):
            if name: yield (name, ''.join(seq))
            name, seq = line, []
        else:
            seq.append(line)
            if name: yield (name, ''.join(seq))


# Command line argument parser
parser = argparse.ArgumentParser(description='Convert trRosetta distance matrix to Rosetta constraints')
parser.add_argument('-n', '--npz', nargs='?', type=str,required=True, help='NpzFile')
parser.add_argument('-f', '--fasta', nargs='?', type=str,required=True, help='Fasta file')
parser.add_argument('-t', '--threshold', type=float, nargs='?', help='Minimum threshold (0-1)', default='0.9')
parser.add_argument('-g', '--graph', type=bool, nargs='?', help='Plot functions(True/False)', default=False)
parser.add_argument('-c', '--constraints', type=bool, nargs='?', help='Save constraint files (True/False)', default=True)
parser.add_argument('-a', '--target', type=str,required=True, help='Target name')
parser.add_argument('-p', '--path', type=str,required=True, help='The full path to the spline files must be specified')
args = parser.parse_args()

plt.ioff()
threshold = args.threshold # minimum probability to select constraints (< 20 A)
plot = args.graph
save = args.constraints
prefix = args.target
path = args.path

# Read distance distribution from trRosetta npz file
with open(args.npz, 'rb') as matrix:
    data = np.load(matrix)
    dis = data['dist']

# Find Gly res in sequence
with open (args.fasta, 'r') as fasta:
    for name, seq in read_fasta(fasta):
        G_pos=[pos for pos, char in enumerate(seq) if char == 'G']


os.chdir(path)

xs = np.arange(2.,20.,0.5)
for idx1, res1 in enumerate(dis):
    for idx2, res2 in enumerate(list(res1)):
        if idx2 < idx1:
            if res2[1:37].sum() > threshold:	# position [0] is prob of being > 20 A
                pot = potential(list(res2[1:37]))
                if not ((idx1 in G_pos) or (idx2 in G_pos)): 
                    if plot:
                        fig, ax1 = plt.subplots()
                        color = 'tab:red'
                        ax1.set_xlabel('Dist (A)')
                        ax1.set_ylabel('Prob', color=color)
                        ax1.plot(xs, list(res2[1:37]), color=color, label='trRosetta')
                        ax1.tick_params(axis='y', labelcolor=color)
                        ax2 = ax1.twinx()            
                        color = 'tab:blue'
                        ax2.set_ylabel('Score', color=color)  
                        ax2.plot(xs, pot, color=color, label='Score')
                        ax2.tick_params(axis='y', labelcolor=color)
                        filename = prefix + '_' + str(idx1) + '_'+ str(idx2)
                        print(filename)
                        plt.legend()
                        plt.show()
                        plt.savefig(filename)
                    else:
                        spline = prefix + '_' + str(idx1+1) + '_'+ str(idx2+1)+ '.spline'
                        constraint = prefix + '.cst'
                        print(spline)
                        x_line = '\t'.join(['x_axis'] + [f"{x:7.2f}" for x in xs])
                        y_line = '\t'.join(['y_axis'] + [f"{x:7.2f}" for x in pot])
                        if save:
                            with open(spline, 'w') as sp, open(constraint, 'a+') as cst:
                                sp.write(x_line + '\n')
                                sp.write(y_line + '\n')
                                cst.write('AtomPair CB ' + str(idx1+1) + ' CB ' + str(idx2+1) + \
                                          ' SPLINE TAG '+ path +'/'+str(prefix)+ \
                                          + spline + ' 1.0 1.0 .1\n')
                        else:
                            print(x_line)
                            print(y_line)

print('Done!')
