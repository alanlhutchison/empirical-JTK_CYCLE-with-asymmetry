#!/usr/bin/env python
"""
This script makes N permutations based on the header of the sample file fed in to it.
Author: Alan L. Hutchison, alanlhutchison@uchicago.edu, Aaron Dinner Group, University of Chicago
"""
import numpy as np
import math
import sys
import random

def main():
    # header is of length 24 with period 12, by 2
    N_perms = 100000
    fn = sys.argv[1]
    fn_out = sys.argv[2]
    header = read_in(fn)
    length = len(header.split())-1
    with open(fn_out,'w') as g:
        g.write(header)
        for i in xrange(0,N_perms):
            alist=range(length)
            random.shuffle(alist)
            g.write("length"+str(length)+"_"+str(i)+"\t"+"\t".join([str(x) for x in alist])+"\n")

#def perm_given_index(alist, apermindex):
#    alist = alist[:]
#    for i in range(len(alist)-1):
#        apermindex, j = divmod(apermindex, len(alist)-i)
#        alist[i], alist[i+j] = alist[i+j], alist[i]
#    return alist

def read_in(fn):
    with open(fn,'r') as f:
        for line in f:
            if line[0]=="#":
                header = line
    return header

if __name__=="__main__":
    main()
