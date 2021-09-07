#!/usr/bin/env python

import argparse

def get_args():
    my_parser = argparse.ArgumentParser(description='counts mapped and unmapped reads in a SAM file')

    my_parser.add_argument('-f',
                        '--file',
                        action='store',
                        help='alignment SAM file',
                        required=True)
    return my_parser.parse_args()

args = get_args()

#initialize counts for mapped and unmapped reads
mapped_count = 0
unmapped_count = 0
with open(args.file, "r") as fl:
    while True:
        line = fl.readline().strip()
        if line == "":
            break
        #lots of headers that start with "@", not wanted in this instance
        if not line.startswith("@"):
            inp = int(line.split("\t")[1])
            #do not double count any reads so check for presence of secondary alignment flag and skip those that have it. Only primary alignments counted.
            if ((inp & 256) != 256):
                if ((inp & 4) != 4):
                    mapped_count += 1
                else:
                    unmapped_count += 1

print(mapped_count)
print(unmapped_count)
