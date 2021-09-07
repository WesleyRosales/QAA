#!/usr/bin/env python

import gzip
import argparse
import matplotlib as mp
import matplotlib.pyplot as plt

#collect read 1 and read 2 trimmed files as output
def get_args():
    my_parser = argparse.ArgumentParser(description='returns histogram of trimmed read lengths')

    my_parser.add_argument('-pr1',
                        '--pairedread1',
                        action='store',
                        help='enter paired read 1 fastq.gz file',
                        required=True)

    my_parser.add_argument('-ur1',
                        '--unpairedread1',
                        action='store',
                        help='enter unpaired read 1 fastq.gz file',
                        required=True)

    my_parser.add_argument('-pr2',
                        '--pairedread2',
                        action='store',
                        help='enter paired read 2 fastq.gz file',
                        required=True)

    my_parser.add_argument('-ur2',
                        '--unpairedread2',
                        action='store',
                        help='enter unpaired read 2 fastq.gz file',
                        required=True) 

    my_parser.add_argument('-o',
                        '--output',
                        action='store',
                        help='enter output name for histogram',
                        required=True)
    return my_parser.parse_args()

args = get_args()

read1_lengths_list = []
read2_lengths_list = []

r1p = gzip.open(args.pairedread1, "rt")
while True:
    line1 = r1p.readline()
    if line1 == "":
        break
    line2 = r1p.readline().strip()
    line3 = r1p.readline()
    line4 = r1p.readline()
    read1_lengths_list.append(len(line2))
r1p.close()

r1up = gzip.open(args.unpairedread1, "rt")
while True:
    line1 = r1up.readline()
    if line1 == "":
        break
    line2 = r1up.readline().strip()
    line3 = r1up.readline()
    line4 = r1up.readline()
    read1_lengths_list.append(len(line2))
r1up.close()

r2p = gzip.open(args.pairedread2, "rt")
while True:
    line1 = r2p.readline()
    if line1 == "":
        break
    line2 = r2p.readline().strip()
    line3 = r2p.readline()
    line4 = r2p.readline()
    read2_lengths_list.append(len(line2))
r2p.close()

r2up = gzip.open(args.unpairedread2, "rt")
while True:
    line1 = r2up.readline()
    if line1 == "":
        break
    line2 = r2up.readline().strip()
    line3 = r2up.readline()
    line4 = r2up.readline()
    read2_lengths_list.append(len(line2))
r2up.close()

#create histogram
plt.rcParams['figure.figsize'] = [16, 6]
plt.hist(x=read1_lengths_list, bins=50, alpha=0.5, label="Read 1")
plt.hist(x=read2_lengths_list, bins=50, alpha=0.5, label="Read 2")
plt.title("Read length distribution")
plt.xlabel("Read Length (bp)")
plt.ylabel("Frequency")
plt.legend(loc="upper left")
plt.savefig(args.output)
