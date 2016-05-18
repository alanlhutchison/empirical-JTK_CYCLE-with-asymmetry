#!/usr/bin/env python
"""
Created on April 20 2014
Updated on April 13 2016
@author: Alan L. Hutchison, alanlhutchison@uchicago.edu, Aaron R. Dinner Group, University of Chicago

This script is one in a series of scripts for running empirical JTK_CYCLE analysis as described in 

Hutchison, Maienschein-Cline, and Chiang et al. Improved statistical methods enable greater sensitivity in rhythm detection for genome-wide data, PLoS Computational Biology 2015 11(3): e 1004094. doi:10.1371/journal.pcbi.1004094

It is a sped-up version of the original provided script using reference waveform pre-calculation and symmetry arguments to reduce computational time.

The resulting *jtkout.txt file contains only the best-matching waveform for each time series (unlike the original script).
Please use ./jtk7.py -h to see the help screen for further instructions on running this script.

"""
VERSION="1.1"

from scipy.stats import kendalltau
import numpy as np
import argparse
import os.path

#from operator import itemgetter
#import sys
#import itertools as it
#import time
#from scipy.stats import norm

def main(args):
    fn          = args.filename
    prefix      = args.prefix
    fn_waveform = args.waveform
    fn_period   = args.period
    fn_phase    = args.phase
    fn_width    = args.width
    fn_null     = args.null
    fn_out      = args.output

    fn_out        = set_fn_out(fn,prefix,fn_out)
    fnw = fn_waveform.split('/')[-1] if '/' in fn_waveform else fn_waveform
    if fnw== 'cosine':
        waveforms = ['cosine']
    else:
        waveforms     = read_in_list(fn_waveform)

        
        
        
    periods       = np.array(read_in_list(fn_period),dtype=float)
    phases        = np.array(read_in_list(fn_phase),dtype=float)
    widths        = np.array(read_in_list(fn_width),dtype=float)
    header,series   = read_in(fn)
    #header,series = organize_data(header,data)
    out_lines = [[]]*len(series)

    dref = make_references(header,waveforms,periods,phases,widths)

    #print header
    
    for i,serie in enumerate(series):
        if [s for s in serie[1:] if s!="NA"]==[]:
            name = [serie[0]]+["All_NA"]+[-10000]*10+[np.nan,np.nan]
        else:
            mmax,mmaxloc,mmin,mminloc,MAX_AMP=series_char(serie,header)
            sIQR_FC = IQR_FC(serie)
            smean   = series_mean(serie)
            sstd    = series_std(serie)
            sFC     = FC(serie)
        
        best = get_best_match(serie,waveforms,periods,phases,widths,dref)
        geneID,waveform,period,phase,nadir,tau,p = best
        out_line =     [geneID,waveform,period,phase,nadir,smean,sstd,mmax,mmaxloc,mmin,mminloc,MAX_AMP,sFC,sIQR_FC,tau,p]
        out_line =     [str(l) for l in out_line]
        out_lines[i] = out_line

    with open(fn_out,'w') as g:
        g.write("ID\tWaveform\tPeriod\tPhase\tNadir\tMean\tStd_Dev\tMax\tMaxLoc\tMin\tMinLoc\tMax_Amp\tFC\tIQR_FC\tTau\tP\n")
        for out_line in out_lines:
            g.write("\t".join(out_line)+"\n")

    fnn = fn_null.split('/')[-1] if '/' in fn_null else fn_null
    
    if fnn=='':
        fn_out_null = fn_out.replace('jtkout','jtknull1000')
        #print len(header)
        null_size = 1000
        out_lines = [[]]*null_size
        for j in xrange(null_size):
            
            ser = np.random.normal(0,1,size=len(header)+1)
            best = get_best_match(ser,waveforms,periods,phases,widths,dref)
            geneID,waveform,period,phase,nadir,tau,p = best
            out_line =     [geneID,waveform,period,phase,nadir,tau,p]
            out_line =     [str(l) for l in out_line]
            out_lines[j] = out_line

        with open(fn_out_null,'w') as g:
            g.write("ID\tWaveform\tPeriod\tPhase\tNadir\tTau\tP\n")
            for out_line in out_lines:
                g.write("\t".join(out_line)+"\n")
            
        return fn_out,fn_out_null
    else:
        return fn_out,fn_null
            
def get_best_match(serie,waveforms,periods,phases,widths,dref):
    best = [serie[0],'cosine',0.,0.,0.,0.,1.]    
    for waveform in waveforms:
        for period in periods:
            pairs = []
            for phase in phases:
                for width in widths:
                    nadir = (phase+width)%period
                    pair = [nadir,phase]
                    if pair not in pairs:
                        pairs.append([phase,nadir])
                        reference = dref[waveform][period][phase][nadir]
                        geneID,tau,p = generate_mod_series(reference,serie)
                        if min(p,best[-1])==p:
                            best = [geneID,waveform,period,phase,nadir,tau,p]
    if best[-2] <0:
        geneID,waveform,period,phase,nadir,tau,p = best
        best = [geneID,waveform,period,nadir,phase,np.abs(tau),p]
    return best


def make_references(header,waveforms,periods,phases,widths):
    dref ={}
    pairs = []
    for waveform in waveforms:
        dref.setdefault(waveform,{})
        for period in periods:
            dref[waveform].setdefault(period,{})                
            for phase in phases:
                dref[waveform][period].setdefault(phase,{})                    
                for width in widths:
                    nadir = (phase+width)%period                        
                    dref[waveform][period][phase].setdefault(nadir,[])                        
                    pair = [nadir,phase]
                    if pair not in pairs:
                        pairs.append([phase,nadir])                            
                        dref[waveform][period][phase][nadir]= generate_base_reference(header,waveform,period,phase,width)
    return dref

            
def set_fn_out(fn,prefix,fn_out):
    if fn_out == "DEFAULT":
        if ".txt" in fn:
            fn_out = fn.replace(".txt","_"+prefix+"_jtkout.txt")
        else:
            fn_out = fn+"_" +prefix + "_jtkout.txt"

    def f_add_on(fn_out):
        add_on = 1
        while os.path.isfile(fn_out):
            print fn_out, "already exists, take evasive action!!!"
            endstr = '.'+fn_out.split('.')[-1]
            mid = '_'+str(add_on)+endstr
            if add_on ==1:
                fn_out = fn_out.replace(endstr,mid)
            else:
                midendstr = '_'+fn_out.split('_')[-1]
                fn_out = fn_out.replace(midendstr,mid)
            add_on = add_on + 1
        return fn_out
    fn_out = f_add_on(fn_out)                
    return fn_out

            
def append_out(fn_out,line):
    line = [str(l) for l in line]
    with open(fn_out,'a') as g:
        g.write("\t".join(line)+"\n")

def write_out(fn_out,output):
    with open(fn_out,'w') as g:
        for line in output:
            g.write(str(line)+"\n")

def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False


def read_in_list(fn):
    with open(fn,'r') as f:
        lines = f.read().splitlines()
    return lines
        
def read_in(fn):
    """Read in data to header and data"""
    with open(fn,'r') as f:
        data=[]
        start_right=0
        for line in f:
            words = line.strip().split()
            words = [word.strip() for word in words]
            if words[0] == "#":
                start_right = 1
                header = words[1:]
            else:
                if start_right == 0:
                    print "Please enter file with header starting with #"
                elif start_right == 1:
                    data.append(words)
    return header, data

#def organize_data(header,data):
#    """
#    Organize list of lists from such that genes with similar time-series holes match (for null distribution calc)
#    Return a header ['#','ZTX','ZTY'...] and a list of lists [ lists with similar holes (identical null distribution) , [],[],[]] 
#    """
#    L = data
#
#    for i in xrange(1,len(header)):
#        L=sorted(L, key=itemgetter(i))
#    return header,L

def generate_base_reference(header,waveform="cosine",period=24,phase=0,width=12):
    """
    This will generate a waveform with a given phase and period based on the header, 
    """
    ZTs = header
    tpoints = np.zeros(len(ZTs))    
    coef = 2.0 * np.pi / float(period)
    w = float(width) * coef
    for i,ZT in enumerate(ZTs):
        z = float(ZT[2:].split("_")[0])
        tpoints[i] =  (z-phase) * coef

    if waveform == "cosine":
        def cosine(x,w):
            x = x % (2*np.pi)
            w = w % (2*np.pi)
            if x <= w:
                y = np.cos(x/(w/np.pi))
            elif x > w:
                y = np.cos( (x+2.*(np.pi-w))*np.pi/ (2*np.pi - w) )
            return y
        #fcos = lambda tp : cosine(tp,w)        
        reference= [cosine(tpoint,w) for tpoint in tpoints]
    
    elif waveform == "trough":
        def trough(x,w):
            x = x % (2*np.pi)
            w = w % (2*np.pi)
            if x <= w:
                y = 1 + -x/w
            elif x > w:
                y = (x-w)/(2*np.pi - w)
            return y
        
        #ftro = lambda tp : trough(tp,w)        
        reference= [trough(tpoint,w) for tpoint in tpoints]
        
    return reference



def IQR_FC(series):
    qlo = __score_at_percentile__(series, 25)
    qhi = __score_at_percentile__(series, 75)
    if (qlo=="NA" or qhi=="NA"):
        return "NA"
    elif (qhi==0):
        return 0
    elif ( qlo==0):
        return "NA"
    else:
        iqr = qhi/qlo
        return iqr

def FC(series):
    series=[float(s) if s!="NA" else 0 for s in series[1:] if s!="NA"  ]
    if series!=[]:
        mmax = max(series)
        mmin = min(series)
        if mmin==0:
            sFC = -10000
        else:
            sFC = mmax / mmin
    else:
        sFC = "NA"
    return sFC


def series_char(series,header):
    """Uses interquartile range to estimate amplitude of a time series."""
    series=[float(s) for s in series[1:] if s!="NA"]
    head = [header[i] for i,s in enumerate(series) if s!="NA"]
    if series!=[]:
        mmax = max(series)
        mmaxloc = head[series.index(mmax)]
        mmin = min(series)
        mminloc = head[series.index(mmin)]
        diff=mmax-mmin
    else:
        mmax = "NA"
        mmaxloc = "NA"
        mmin = "NA"
        mminloc = "NA"
        diff = "NA"
    return mmax,mmaxloc,mmin,mminloc,diff


def series_mean(series):
    """Finds the mean of a timeseries"""
    series = [float(s) for s in series[1:] if s!="NA"]
    return np.mean(series)

def series_std(series):
    """Finds the std dev of a timeseries"""
    series = [float(s) for s in series[1:] if s!="NA"]
    return np.std(series)

def __score_at_percentile__(ser, per):
    ser = [float(se) for se in ser[1:] if se!="NA"]
    if len(ser)<5:
        score ="NA"
        return score
    else: 
        ser = np.sort(ser)
        i = int(per/100. * len(ser))
        if (i % 1 == 0):
            score = ser[i]
        else:
            interpolate = lambda a,b,frac: a + (b - a)*frac
            score = interpolate(ser[int(i)], ser[int(i) + 1], i % 1)
        return float(score)

def generate_mod_series(reference,series):
    """
    Takes the series from generate_base_null, takes the list from data, and makes a null
    for each gene in data or uses the one previously calculated.
    Then it runs Kendall's Tau on the exp. series against the null
    """
    
    geneID = series[0]
    values = series[1:]

    if len(reference)!=len(values):
        geneID = 'blank'
        return geneID,0,1

    tau,p=kendalltau(values,reference)
    return geneID,tau,p


def __create_parser__():
    p = argparse.ArgumentParser(
        description="Python script for running empirical JTK_CYCLE with asymmetry search as described in Hutchison, Maienschein-Cline, and Chiang et al. Improved statistical methods enable greater sensitivity in rhythm detection for genome-wide data, PLoS Computational Biology 2015 11(3): e1004094. This script was written by Alan L. Hutchison, alanlhutchison@uchicago.edu, Aaron R. Dinner Group, University of Chicago.",
        epilog="Please contact the correpsonding author if you have any questions.",
        version=VERSION
        )


    analysis = p.add_argument_group(title="JTK_CYCLE analysis options")
    analysis.add_argument("-o", "--output",
                   dest="output",
                   action='store',
                   metavar="filename string",
                   type=str,
                   default = "DEFAULT",
                   help="You want to output something. If you leave this blank, _jtkout.txt will be appended to your filename")

    analysis.add_argument("-f", "--filename",
                   dest="filename",
                   action='store',
                   metavar="filename string",
                   type=str,
                   help='This is the filename of the data series you wish to analyze.\
                   The data should be tab-spaced. The first row should contain a # sign followed by the time points with either CT or ZT preceding the time point (such as ZT0 or ZT4). Longer or shorter prefixes will not work. The following rows should contain the gene/series ID followed by the values for every time point. Where values are not available NA should be put in it\'s place.')


    analysis.add_argument("-x","--prefix",
                          dest="prefix",
                          type=str,
                          metavar="string",
                          action='store',
                          default="",
                          help="string to be inserted in the output filename for this run")


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
                          help='Should be a file with phases you wish to search for listed in a single column separated by newlines.\
                          Provided file is "period_24.txt. Will default to 24"')

    analysis.add_argument("-n","--null",
                          dest="null",
                          metavar="filename string",
                          type=str,
                          action='store',
                          default="",
                          help='A name for the file generated by the calculation of the null distrubition of eJTK, which can be used with CalcP.py')

    
    return p



if __name__=="__main__":
    parser = __create_parser__()
    args = parser.parse_args()
    main(args)
