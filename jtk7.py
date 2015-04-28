#!/usr/bin/env python
"""
Created on April 20 2014
@author: Alan L. Hutchison, alanlhutchison@uchicago.edu, Aaron Dinner Group, University of Chicago

This script is one in a series of scripts for running empirical JTK_CYCLE analysis as described in 

Hutchison, Maienschein-Cline, and Chiang et al. Improved statistical methods enable greater sensitivity in rhythm detection for genome-wide data, PLoS Computational Biology 2014 (in review).

Please use ./jtk7.py -h to see the help screen for further instructions on running this script.

"""
VERSION="1.00"

from scipy.stats import kendalltau
from operator import itemgetter
import numpy as np
import sys
import argparse
import itertools as it
import time

#import matplotlib.pyplot as plt
#import matplotlib.cm as cm
from scipy.stats import norm
import os.path

def main(args):
    fn = args.filename
    prefix = args.prefix
    fn_waveform = args.waveform
    fn_period = args.period
    fn_phase = args.phase
    fn_width = args.width
    fn_out = args.output
    if fn_out == "DEFAULT":
        if ".txt" in fn:
            fn_out=fn.replace(".txt","_"+prefix+"_jtkout.txt")
        else:
            fn_out = fn+"_" +prefix + "_jtkout.txt"
        print fn
    add_on = 1
    while os.path.isfile(fn_out):
        print fn_out, "already exists, take evasive action!!!"
        fn_out = ".".join(fn_out.split(".")[0:-1])+"_"+str(add_on)+".txt"
        add_on = add_on + 1

    waveforms = read_in_list(fn_waveform)
    periods = read_in_list(fn_period)
    phases = read_in_list(fn_phase)
    widths = read_in_list(fn_width)
    #fn_out = "\t".join(fn.replace("jtkprepared","").split(".")[0:-1])+"_jtkout_emprical.txt"
    header,data = read_in(fn)
    header,series = organize_data(header,data)
    RealKen = KendallTauP()        
    #output = ["ID\tWaveform\tPeriod\tPhase\tAsymmetry\tMean\tStd_Dev\tMax\tMin\tMax_Amp\tFC\tIQR_FC\tTau\tempP"]
    Ps = []
    with open(fn_out,'w') as g:
        g.write("ID\tWaveform\tPeriod\tPhase\tAsymmetry\tMean\tStd_Dev\tMax\tMin\tMax_Amp\tFC\tIQR_FC\tTau\tP\n")
        for serie in series:
            if [s for s in serie[1:] if s!="NA"]==[]:
                name = [serie[0]]+["All_NA"]+[-10000]*10+[np.nan,np.nan]
            else:
                mmax,mmin,MAX_AMP=series_char(serie)
                sIQR_FC=IQR_FC(serie)
                smean = series_mean(serie)
                sstd = series_std(serie)
                sFC = FC(serie)

            local_ps = [] 
            for waveform in waveforms:
                for period in periods:
                    for phase in phases:
                        for width in widths:
                            reference = generate_base_reference(header,waveform,period,phase,width)
                            geneID,tau,p = generate_mod_series(reference,serie,RealKen)
                            out_line = [geneID,waveform,period,phase,width,smean,sstd,mmax,mmin,MAX_AMP,sFC,sIQR_FC,tau,p]
                            out_line = [str(l) for l in out_line]
                            g.write("\t".join(out_line)+"\n")

         #local_ps = sorted(local_ps)
         #best = min(local_ps)
         #Ps.append(best)
         #append_out(fn_out,best)
         #name = [geneID,waveform,period,phase,width,smean,sstd,mmax,mmin,MAX_AMP,sFC,sIQR_FC,tau,empirical_p]
         #name = [str(n) for n in name]
         #print "\t".join(name)
         #print time.asctime( time.localtime(time.time()) )
         #output.append("\t".join(name))
     #write_out(fn_out,Ps)

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

def organize_data(header,data):
    """
    Organize list of lists from such that genes with similar time-series holes match (for null distribution calc)
    Return a header ['#','ZTX','ZTY'...] and a list of lists [ lists with similar holes (identical null distribution) , [],[],[]] 
    """
    L = data

    for i in xrange(1,len(header)):
        L=sorted(L, key=itemgetter(i))
    return header,L

def generate_base_reference(header,waveform="trough",period=24,phase=0,width=12):
    """
    This will generate a waveform with a given phase and period based on the header, 
    """
    tpoints = []
    ZTs = header
    coef = 2.0 * np.pi / float(period)
    w = float(width) * coef
    for ZT in ZTs:
        z = ZT[2:].split("_")[0]
        tpoints.append( (float(z)-float(phase) ) * coef)


    def ramp_down(x,w=3*np.pi/2):
        x = x % (2*np.pi)
        y=max(-1*x/w + 1.0, 0.0)
        return y
    def ramp_up(x,w=3*np.pi/2):
        x = x % (2*np.pi)
        y =min(x/w, 1.0)
        return y
    def impulse(x,w=3*np.pi/4):
        x = x % (2*np.pi)
        d = min(x,np.abs(np.pi*2-x))
        y = max(-2.*d/w + 1.0,0.0)
        return y
    def step(x,w=np.pi):
        x = x % (2*np.pi)
        y = 1.0 if x < w else 0.0
        return y
    def trough(x,w):
        x = x % (2*np.pi)
        w = w % (2*np.pi)
        if x <= w:
            y = 1 + -x/w
        elif x > w:
            y = (x-w)/(2*np.pi - w)
        return y

    def cosine(x,w):
        x = x % (2*np.pi)
        w = w % (2*np.pi)
        if x <= w:
            y = np.cos(x/(w/np.pi))
        elif x > w:
            y = np.cos( (x+2.*(np.pi-w))*np.pi/ (2*np.pi - w) )
        return y


    #print tpoints
    #print [tpoint/np.pi/2.0 for tpoint in tpoints]
    #print phase
    #print waveform
    if waveform == "cosine":
        reference=[cosine(tpoint,w) for tpoint in tpoints]
    elif waveform == "impulse":
        reference=[impulse(tpoint,w) for tpoint in tpoints]
    elif waveform == "rampup":
        reference=[ramp_up(tpoint,w) for tpoint in tpoints]
    elif waveform == "rampdown":
        reference=[ramp_down(tpoint,w) for tpoint in tpoints]
    elif waveform == "step":
        reference=[step(tpoint,w)     for tpoint in tpoints]
    elif waveform == "trough":
        reference=[trough(tpoint,w) for tpoint in tpoints]
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


#def find_peak(series):
#    """Finds maximum of series"""
#    series=[float(s) if s!="NA" else 0 for s in series[1:] if s!="NA"  ]
#    mmax = max(series)
#    max_index = series.index(mmax)
#    return max_index

#def find_trough(series):
#    """Finds minimum of series"""
#    series=[float(s) if s!="NA" else 0 for s in series[1:] if s!="NA"  ]
#    mmin = min(series)
#    min_index = series.index(mmin)
#    return min_index
    

#def maximum(series):
#    """Finds maximum of series"""
#    series=[float(s) for s in series[1:] if s!="NA"]
#    mmax = max(series)
#    return mmax


def series_char(series):
    """Uses interquartile range to estimate amplitude of a time series."""
    series=[float(s) for s in series[1:] if s!="NA"]
    if series!=[]:
        mmax = max(series)
        mmin = min(series)
        diff=mmax-mmin
    else:
        mmax = "NA"
        mmin = "NA"
        diff = "NA"
    return mmax,mmin,diff


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
        i = (per/100. * len(ser))
        if (i % 1 == 0):
            score = ser[i]
        else:
            interpolate = lambda a,b,frac: a + (b - a)*frac
            score = interpolate(ser[int(i)], ser[int(i) + 1], i % 1)
        return float(score)

def generate_mod_series(reference,series,RealKen):
    """
    Takes the series from generate_base_null, takes the list from data, and makes a null
    for each gene in data or uses the one previously calculated.
    Then it runs Kendall's Tau on the exp. series against the null
    """

    geneID = series[0]
    values = series[1:]
    binary = np.array([1.0 if value!="NA" else np.nan for value in values])
    reference = np.array(reference)
    temp = reference*binary
    mod_reference = [value for value in temp if not np.isnan(value)]
    mod_values = [float(value) for value in values if value!='NA']

    if len(mod_values) < 3:
        tau,p = np.nan,np.nan
    elif mod_values.count(np.nan) == len(mod_values):
        tau,p = np.nan,np.nan
    elif mod_values.count(0) == len(mod_values):
        tau,p = np.nan,np.nan
    #elif sum(mod_values)<0.00001:
    #    tau,p = np.nan,np.nan        
    else:
        tau,p=kendalltau(mod_values,mod_reference)
        if not np.isnan(tau):
            if len(mod_values) < 150:
                pk = RealKen.pval(tau,len(mod_values))
                if pk!=None:
                    p=pk
            else:
                p = p / 2.0
                if tau < 0:
                    p = 1-p
                    

    #print tau,p
    return geneID,tau,p




def __create_parser__():
    p = argparse.ArgumentParser(
        description="Python script for running empirical JTK_CYCLE with asymmetry search as described in Hutchison, Maienschein-Cline, and Chiang et al. Improved statistical methods enable greater sensitivity in rhythm detection for genome-wide data, PLoS Computational Biology 2014 (in review). This script was written by Alan L. Hutchison, alanlhutchison@uchicago.edu, Aaron Dinner Group, University of Chicago.",
        epilog="Please contact the correpsonding author if you have any questions.",
        version=VERSION
        )

                   
    #p.add_argument("-t", "--test",
    #               action='store_true',
    #               default=False,
    #               help="run the Python unittest testing suite")
    p.add_argument("-o", "--output",
                   dest="output",
                   action='store',
                   metavar="filename string",
                   type=str,
                   default = "DEFAULT",
                   help="You want to output something. If you leave this blank, _jtkout.txt will be appended to your filename")



    analysis = p.add_argument_group(title="JTK_CYCLE analysis options")

    analysis.add_argument("-f", "--filename",
                   dest="filename",
                   action='store',
                   metavar="filename string",
                   type=str,
                   help='This is the filename of the data series you wish to analyze.\
                   The data should be tab-spaced. The first row should contain a # sign followed by the time points with either CT or ZT preceding the time point (such as ZT0 or ZT4). Longer or shorter prefixes will not work. The following rows should contain the gene/series ID followed by the values for every time point. Where values are not available NA should be put in it\'s place.')


    analysis.add_argument("--prefix",
                          dest="prefix",
                          type=str,
                          metavar="string",
                          action='store',
                          default="",
                          help="string to be inserted in the output filename for this run")


    analysis.add_argument("--waveform",
                          dest="waveform",
                          type=str,
                          metavar="filename string",
                          action='store',
                          default="cosine",
                          #choices=["wavform_cosine.txt","waveform_rampup.txt","waveform_rampdown.txt","waveform_step.txt","waveform_impulse.txt","waveform_trough.txt"],
                          help='Should be a file with asymmetries (widths) you wish to search for listed in a single column separated by newlines.\
                          Options include cosine (dflt), rampup, rampdown, impulse, step, trough')

    analysis.add_argument("-w", "--width", "-a", "--asymmetry",
                          dest="width",
                          type=str,
                          metavar="filename string",
                          action='store',
                          default="widths_02-22.txt",
                          #choices=["widths_02-22.txt","widths_04-20_by4.txt","widths_04-12-20.txt","widths_08-16.txt","width_12.txt"]
                          help='Should be a file with asymmetries (widths) you wish to search for listed in a single column separated by newlines.\
                          Provided files include "widths_02-22.txt","widths_04-20_by4.txt","widths_04-12-20.txt","widths_08-16.txt","width_12.txt"')
    analysis.add_argument("-ph", "--phase",
                          dest="phase",
                          metavar="filename string",
                          type=str,
                          default="phases_00-22.txt",
                          #choices=["phases_00-22.txt","phases_00-22_by2.txt","phases_00-20.txt","phases_00-20_by4.txt"]
                          help='Should be a file with phases you wish to search for listed in a single column separated by newlines.\
                          Provided files include "phases_00-22.txt" or "phases_00-22_by2.txt" or "phases_00-20.txt" or "phases_00-20_by4.txt"')

    analysis.add_argument("-p","--period",
                          dest="period",
                          metavar="filename string",
                          type=str,
                          action='store',
                          default="period_24.txt",
                          help='Should be a file with phases you wish to search for listed in a single column separated by newlines.\
                          Provide file is "period_24.txt"')
    
    distribution = analysis.add_mutually_exclusive_group(required=False)
    distribution.add_argument("-e", "--exact",
                              dest="harding",
                              action='store_true',
                              default=False,
                              help="use Harding's exact null distribution (dflt)")
    distribution.add_argument("-n", "--normal",
                              dest="normal",
                              action='store_true',
                              default=False,
                              help="use normal approximation to null distribution")
    
    
    return p


# instantiate class to precalculate distribution
# usage: 
#   K = KendallTauP()
#   pval = K.pval(tau,n,two_tailed=True)
class KendallTauP:
    def __init__(self,N=150):        
        # largest number of samples to precompute
        self.N = N
        Nint = self.N*(self.N-1)/2

        # first allocate freq slots for largest sample array
        # as we fill this in we'll save the results for smaller samples

        # total possible number of inversions is Nint + 1
        freqN = np.zeros(Nint + 1)
        freqN[0] = 1.0

        # save results at each step in freqs array
        self.freqs = [np.array([1.0])]
        for i in xrange(1,self.N):
            last = np.copy(freqN)
            for j in xrange(Nint+1):
                # update each entry by summing over i entries to the left
                freqN[j] += sum(last[max(0,j-i):j])
            # copy current state into freqs array
            # the kth entry of freqs should have 1+k*(k-1)/2 entries
            self.freqs.append(np.copy(freqN[0:(1+(i+1)*i/2)]))
            
        # turn freqs into cdfs
        # distributions still with respect to number of inversions
        self.cdfs = []
        for i in xrange(self.N):
            self.cdfs.append(np.copy(self.freqs[i]))
            # turn into cumulative frequencies
            for j in xrange(1,len(self.freqs[i])):
                self.cdfs[i][j] += self.cdfs[i][j-1]
            # convert freqs to probs
            self.cdfs[i] = self.cdfs[i]/sum(self.freqs[i])
            
    # plot exact distribution compared to normal approx
    def plot(self,nlist):
        colors = cm.Set1(np.linspace(0,1,len(nlist)))

        # for plotting gaussian
        x = np.linspace(-1.2,1.2,300)
        # plot pdfs
        plt.figure()
        for i in xrange(len(nlist)):
            ntot = len(self.freqs[nlist[i]-1])-1
            tauvals = (ntot - 2.0*np.arange(len(self.freqs[nlist[i]-1])))/ntot
            probs = ((ntot+1.0)/2.0)*self.freqs[nlist[i]-1]/sum(self.freqs[nlist[i]-1])
            plt.scatter(tauvals,probs,color=colors[i])
            # now plot gaussian comparison
            var = 2.0*(2.0*nlist[i]+5.0)/(nlist[i]*(nlist[i]-1)*9.0)
            plt.plot(x,norm.pdf(x,0.0,np.sqrt(var)),color=colors[i])
        plt.legend(nlist,loc='best')
        # plt.savefig('pdfs.png')
        plt.show()

        # now plot cdfs
        plt.figure()
        for i in xrange(len(nlist)):
            ntot = len(self.freqs[nlist[i]-1])-1
            tauvals = -1.0*(ntot - 2.0*np.arange(len(self.freqs[nlist[i]-1])))/ntot
            probs = self.cdfs[nlist[i]-1]
            plt.scatter(tauvals,probs,color=colors[i])
            # now plot gaussian comparison
            var = 2.0*(2.0*nlist[i]+5.0)/(nlist[i]*(nlist[i]-1)*9.0)
            plt.plot(x,norm.cdf(x,0.0,np.sqrt(var)),color=colors[i])
        plt.legend(nlist,loc='best')
        # plt.savefig('cdfs.png')
        plt.show()

    # use cdfs to return pval
    # default to return two tailed pval
    def pval(self,tau,n,two_tailed=False):
        # enforce tau is between -1 and 1
        if tau <= -1.000001 or tau >= 1.000001:
            sys.stderr.write(str(type(tau))+"\n")
            sys.stderr.write(str(tau)+"\n")
            sys.stderr.write("invalid tau\n")
            #print 'invalid tau'
            return None
        # enforce n is less than our precomputed quantities
        if n > self.N:
            #print 'n is too large'
            sys.stderr.write("n is too large/n")
            return None

        # convert tau to value in terms of number of inversions
        ntot = n*(n-1)/2
        inv_score = int(round((ntot - tau * ntot)/2.0))
        # I'm a little worried about the precision of this,
        # but probably not enough to be really worried for reasonable n
        # since we really only need precision to resolve ntot points

        # if two tailed, we're getting a tail from a symmetric dist
        min_inv_score = min(inv_score,ntot-inv_score)

        if two_tailed:
            pval = self.cdfs[n-1][min_inv_score]*2.0
        else:
            # if one tailed return prob of getting that or fewer inversions
            pval = self.cdfs[n-1][inv_score]

        # if inv_score is 0, might have larger than 0.5 prob
        return min(pval,1.0)



if __name__=="__main__":
    parser = __create_parser__()
    args = parser.parse_args()
    main(args)
