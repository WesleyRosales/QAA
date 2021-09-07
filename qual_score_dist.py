#!/usr/bin/env python

import Bioinfo
import gzip
import argparse

#collect read 1 and read 2 files as input
def get_args():
    my_parser = argparse.ArgumentParser(description='return per base quality score figures')

    my_parser.add_argument('-r1',
                        '--read1',
                        action='store',
                        help='enter read 1 gz file',
                        required=True)

    my_parser.add_argument('-r2',
                        '--read2',
                        action='store',
                        help='enter read 2 gz file',
                        required=True)
    
    my_parser.add_argument('-o',
                        '--output',
                        action='store',
                        help='enter output directory',
                        required=True)
    return my_parser.parse_args()

args = get_args()

read1_score_sums = [0 for f in range(101)]
read2_score_sums = [0 for f in range(101)]


fl1 = gzip.open(args.read1, "rt")
fl2 = gzip.open(args.read2, "rt")

#just need to have the sums of the quality scores, not all the scores themselves
#this block continuously sums of quality scores for each base position until the end of the file
record_count = 0
while True:
    
    #run through lines in files until quality score is reached (every 4th line)
    read1inp = fl1.readline()
    read1inp = fl1.readline()
    read1inp = fl1.readline()
    #this one will be the quality score
    read1inp = fl1.readline().strip()
    if read1inp == "":
        break

    read2inp = fl2.readline()
    read2inp = fl2.readline()
    read2inp = fl2.readline()
    #this one will be the quality score
    read2inp = fl2.readline().strip()
    
    for h in range(101):
        read1_score_sums[h] += Bioinfo.convert_phred(read1inp[h])
        read2_score_sums[h] += Bioinfo.convert_phred(read2inp[h])

    record_count += 1

#Closing the opened files since we are done with them now
fl1.close()
fl2.close()

#mean the sums
read1_mean = []
for number in read1_score_sums:
    read1_mean.append(number/record_count)
read2_mean = []
for number in read2_score_sums:
    read2_mean.append(number/record_count)

#create lists of read and index base numbers
read_base_num = [num for num in range(101)]

#create a histogram plot for read 1, read 2, index 1, and index 2 base position quality means
import matplotlib as mp
import matplotlib.pyplot as plt
plt.rcParams['figure.figsize'] = [16, 6]
plt.bar(x=read_base_num,height=read1_mean, width=0.5)
plt.title("Mean Quality Score of Base Pairs At Each Position of Read 1")
plt.xlabel("# Base Pair")
plt.ylabel("Mean Quality Score")
plt.ylim(0,45)
plt.savefig(args.output + "/Read1_hist.png")
plt.clf()

plt.rcParams['figure.figsize'] = [16, 6]
plt.bar(x=read_base_num,height=read2_mean, width=0.5)
plt.title("Mean Quality Score of Base Pairs At Each Position of Read 2")
plt.xlabel("# Base Pair")
plt.ylabel("Mean Quality Score")
plt.ylim(0,45)
plt.savefig(args.output + "/Read2_hist.png")
plt.clf()

