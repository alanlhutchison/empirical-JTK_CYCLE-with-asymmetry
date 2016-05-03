#!/usr/bin/env python
"""
This script takes a jtkout file finds the best P, outputs it with the ending ".best_ps.txt"
Author: Alan L. Hutchison, alanlhutchison@uchicago.edu, Aaron Dinner Group, University of Chicago
"""

#import statsmodels.stats.multitest as ssm
from operator import itemgetter
import sys


def main():
    fns = sys.argv[1:]
    for fn in fns:
        print fn
        fn_out = ".".join(fn.split(".")[0:-1])+".best_ps.txt"
        read_in_write_out(fn,fn_out)

def read_in_write_out(fn,fn_out):
    data = []
    seen = []
    old_pair = []
    previous = ""
    with open(fn,'r') as f:
        with open(fn_out,'w') as g:
            for line in f:
                words = line.strip().split()
                if words[0]=="ID" or words[0]=="#ID":
                    header=words
                else:
                    if words[0] not in seen:
                        seen.append(words[0])
                        pair = [words[0],float(words[-1])]
                        if old_pair!=[]:
                            g.write("{0}\t{1}\n".format(old_pair[0],str(old_pair[1])))
                        old_pair = pair
                    else:
                        if float(words[-1]) < old_pair[1]:
                            old_pair = [words[0],float(words[-1])]
            g.write("{0}\t{1}\n".format(old_pair[0],str(old_pair[1])))

if __name__=="__main__":
    main()

#def get_name_set(data):
#    names = []
#    for item in data:
#        names.append(item[0])
#    names = set(names)
#    return names

#def sort_name(data):
#    good = []
#    bad = []
#    print "in beginning"
#    for line in data:
#        if line[-1]==-10000:
#            bad.append(line)
#        else:
#            line = [float(num) if is_number(num) else num for num in line]
#            good.append(line)
#    print "in sort name"
#    data = None
#    good = sorted(good,key= itemgetter(-1))
#    data= good + bad
#    return data

#def is_number(s):
#    try:
#        float(s)
#        return True
#    except ValueError:
#        return False

#def best_P(data):
#    # data should already be ordered at this point
#    seen=[]
#    data_uniq =[]
#    for datum in data:
#        if datum[0] not in seen:
#            #print "Good", datum
#            seen.append(datum[0])
#            data_uniq.append(datum)
#    return data_uniq,seen

#def write_out_best_P(data,header,fn_out):
#    new_header = "\t".join(header)
#    with open(fn_out,'w') as g:
#        g.write(new_header+"\n")
#        for datum in data:
#            datum = [str(d) for d in datum]
#            g.write("\t".join(datum)+"\n")



