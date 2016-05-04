#!/usr/bin/env python
"""
Created on May 2 2016
@author: Alan L. Hutchison, alanlhutchison@uchicago.edu, Aaron R. Dinner Group, University of Chicago

"""
VERSION="1.0"

import argparse
import sys
import os
#sys.path.insert(1,'bin/')
binpath=os.path.join(os.path.dirname(sys.argv[0]),'bin/')
sys.path.insert(1,binpath)



import eJTK
import CalcP

def main(args):
    fn          = args.filename
    prefix      = args.prefix
    fn_waveform = args.waveform
    fn_period   = args.period
    fn_phase    = args.phase
    fn_width    = args.width
    fn_null     = args.null
    fn_out      = args.output
    fit         = args.fit
    fn_out,fn_null_out = eJTK.main(args)

    args.filename = fn_out
    args.null = fn_null_out

    CalcP.main(args)
    
    
    

def __create_parser__():
    p = argparse.ArgumentParser(
        description="Python script for running empirical JTK_CYCLE with asymmetry search as described in Hutchison, Maienschein-Cline, and Chiang et al. Improved statistical methods enable greater sensitivity in rhythm detection for genome-wide data, PLoS Computational Biology 2015 11(3): e1004094. This script was written by Alan L. Hutchison, alanlhutchison@uchicago.edu, Aaron R. Dinner Group, University of Chicago.",
        epilog="Please contact the correpsonding author if you have any questions.",
        version=VERSION
        )

    analysis = p.add_argument_group(title="JTK_CYCLE analysis options")

    analysis.add_argument("-f", "--filename",
                   dest="filename",
                   action='store',
                   metavar="filename string",
                   type=str,
                   help='This is the filename of the data series you wish to analyze.\
                   The data should be tab-spaced. The first row should contain a # sign followed by the time points with either CT or ZT preceding the time point (such as ZT0 or ZT4). Longer or shorter prefixes will not work. The following rows should contain the gene/series ID followed by the values for every time point. Where values are not available NA should be put in it\'s place.')

    analysis.add_argument("-o", "--output",
                   dest="output",
                   action='store',
                   metavar="filename string",
                   type=str,
                   default = "DEFAULT",
                   help="You want to output something. If you leave this blank, _jtkout.txt will be appended to your filename")

    analysis.add_argument("-x","--prefix",
                          dest="prefix",
                          type=str,
                          metavar="string",
                          action='store',
                          default="",
                          help="string to be inserted in the output filename")


    analysis.add_argument("-w","--waveform",
                          dest="waveform",
                          type=str,
                          metavar="filename string",
                          action='store',
                          default="cosine",
                          help='Should be a file with waveforms  you wish to search for listed in a single column separated by newlines.\
                          Options include cosine (dflt), trough')

    analysis.add_argument("--width", "-a", "--asymmetry",
                          dest="width",
                          type=str,
                          metavar="filename string",
                          action='store',
                          default="widths_02-22.txt",
                              help='Should be a file with asymmetries (widths) you wish to search for listed in a single column separated by newlines.\
                          Provided files include files like "widths_02-22.txt","widths_04-20_by4.txt","widths_04-12-20.txt","widths_08-16.txt","width_12.txt"\nasymmetries=widths')
    
    analysis.add_argument("-s","-ph", "--phase",
                          dest="phase",
                          metavar="filename string",
                          type=str,
                          default="phases_00-22_by2.txt",
                          help='Should be a file with phases you wish to search for listed in a single column separated by newlines.\
                          Example files include "phases_00-22_by2.txt" or "phases_00-22_by4.txt" or "phases_00-20_by4.txt"')


    
    analysis.add_argument("-p","--period",
                          dest="period",
                          metavar="filename string",
                          type=str,
                          action='store',
                          default="24",
                          help='Should be a file with periods you wish to search for listed in a single column separated by newlines.\
                          Provided file is "period_24.txt. Will default to 24"')
    
    analysis.add_argument("-n","--null",
                          dest="null",
                          metavar="filename string",
                          type=str,
                          action='store',
                          default="",
                          help='A jtkout filename with a Tau column from which the null distribution will be calculated in CalcP.py')


    analysis.add_argument("-t","--fit",
                          dest="fit",
                          action='store_true',
                          default=False,
                          help='Boolean option, will use conservative fitting in CalcP if used.')
    
    return p

if __name__=="__main__":
    parser = __create_parser__()
    args = parser.parse_args()
    main(args)
