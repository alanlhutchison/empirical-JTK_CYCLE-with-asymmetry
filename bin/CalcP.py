#!/usr/bin/env python
"""
Created on March 1, 2016
@author: Alan L. Hutchison, alanlhutchison@uchicago.edu, Aaron R. Dinner Group, University of Chicago

Use ./CalcP.py -h to see usage

Credit for the arbfit code goes to Nablaquabla
"""
import numpy as np
import matplotlib.pylab as plt
from mpfit import mpfit

VERSION="0.9"

from scipy.stats import kendalltau
from operator import itemgetter
import numpy as np
import scipy.stats as ss
import sys
import argparse
import itertools as it
import time
import arbfit
import pickle
arbFit = arbfit.arbFit
#import matplotlib.pyplot as plt
#import matplotlib.cm as cm
from scipy.stats import norm
import os.path
import scipy.stats as ss
import pandas as pd
import statsmodels.stats.multitest as ssm

def main(args):
    #def calcP(fn_jtk,pkl):

    #fn_ser = sys.argv[1]
    fn_jtk = args.filename
    fn_pkl = args.null
    fit = args.fit
    
    #ser = pd.read_table(fn_ser,index_col=0)
    #NUM = ser.shape[1]
    jtk = pd.read_table(fn_jtk,index_col='ID')

    
    if '.pkl' in fn_pkl:
        params,taus = pickle.load(open(fn_pkl,'rb'))
    else:
        if 'boot' not in fn_pkl:
            taus = pd.read_table(fn_pkl)['Tau']

            keys,intvalues,yerr,p0,limit = prepare(taus)            
        else:
            taus = pd.read_table(fn_pkl)['TauMean']
            keys,intvalues,yerr,p0,limit = prepare(taus)
        if fit:
            for _ in xrange(2):
                params = GammaFit(keys,intvalues,yerr,p0,limit)
                p0 = params[0]
    params = p0
    gd = ss.gamma(params[0],params[1],params[2])

    if 'boot' not in fn_jtk:
        keys = jtk['Tau']
    else:
        keys = jtk['TauMean']
    print p0

    
    empPs = empP(keys,taus)
    jtk['empP']=empPs

    ps = gd.sf(keys)
    jtk['GammaP'] = ps
    jtk['GammaBH'] = list(ssm.multipletests(ps,method='fdr_bh')[1])
    
    fn_out = fn_jtk.replace('.txt','_GammaP.txt')
    jtk.to_csv(fn_out,sep='\t')

    
def empP(taus,emps):
    taus = np.array(taus)
    emps = np.array(emps)
    ps = [(np.sum(emps>=t)+1)/float(len(emps)+1) for t in taus]
    return np.array(ps)


def prepare(taus):
    i = 0
    #print NUM
    #TOT = NUM*(NUM-1)/2
    d_hugh = {}
    for tau in taus:
        if tau not in d_hugh:
            d_hugh[tau] = 0
        d_hugh[tau]+=1
    keys = sorted(d_hugh.keys())
    values = [d_hugh[key] for key in keys]

    #intkeys = [int(np.round((o+1.)/2.*NUM*(NUM-1)/2,0)) for o in keys]
    intvalues = [v/float(np.sum(values)) for v in values]
    gd = lambda x,p: ss.gamma(p[0],p[1],p[2]).cdf(x)
    yerr = [1e-5/(1*i+1)]*(len(intvalues)-sum(np.cumsum(intvalues)>0.9))+[1e-6/(1*i+1)]*sum(np.cumsum(intvalues)>0.9)

    a,b,c = ss.gamma.fit(taus)
    #a = np.mean(taus)**2/np.var(taus)
    #b = 1e-8
    #c = np.var(taus)/(np.mean(taus))
    p0 = [a,b,c]
    ind = list(np.cumsum(intvalues)>0.9).index(1)
    limit = keys[ind]    
    
                    
    return keys,intvalues,yerr,p0,limit

def GammaFit(x,ydata,yerr,p0,limit):

    gd =  lambda x,p : ss.gamma(p[0],p[1],p[2]).cdf(x)

    #=========================================================================#
    #     --- 'Magic' happens here ---
    #=========================================================================#
    # Create a vectorized function that checks whether a model value 
    # is larger than a limit set by you. So if the model is smaller than the limit provided
    # the fit function will return the simple chi^2
    def checkLimit(x,lim,y,model):
        if x > lim and y<model:
            return 1e5
        else:
            return 1.0
    checkLimitVectorized = np.vectorize(checkLimit)

    # Add a limit=None argument to the fitfunc. The model is your function that you
    # want to fit to the data. In the return statement it now checks for each data
    def fitfunc(p, fjac=None, x=None, y=None, err=None, limit=None):
        model = gd(x,p)
        status = 0
        return [status, checkLimitVectorized(x,limit,y,model)*(y-model)/err]

    #=========================================================================#
    #  Initialize fit info dictionary and try to fit function to data
    #=========================================================================# 

    # Create an info dictionary for each parameter in p0
    # If you want to add some bounds to the parameters you can use the limited and
    # limits keyword to tell mpfit that a parameter is actually bound and what
    # bounds it uses so parinfo[2]['limited'] = [1,0] would mean that p0[2] has a
    # lower bound. parinfo[2]['limits'] = [-10,0] would mean that this lower bound
    # is -10.
    #print p0
    parbase1 = {'value':0,'fixed':0,'limited':[1,1],'limits':[0.,1000.]}
    parbase2 = {'value':0,'fixed':0,'limited':[1,1],'limits':[0.,1000.]}
    parbase3 = {'value':0,'fixed':0,'limited':[1,1],'limits':[0.,1000.]}

    parinfo = [parbase1,parbase2,parbase3]
    for i in range(len(p0)):
        parinfo[i]['value'] = p0[i]
        parinfo[i]['limits'] = [0.8*p0[i],1.2*p0[i]]

    # Define the function arguments that you want to pass to your fit. Here you have
    # to adjust the limit keyword properly.
    #print limit
    fa = {'x': x, 'y': ydata, 'err': yerr, 'limit': limit}

    # Perform the fit using mpfit
    #print 'p0 is',p0
    #print 'fitfunc is',fitfunc
    m = mpfit(fitfunc, p0, parinfo=parinfo, functkw=fa,quiet=1)
    #print 'm is', m
    # Get degrees of freedom. This assumes that you don't use any bounds on your fit
    # Otherwise you have to substract them from your dof. The -1 takes care of the
    # additional overall limit that you impose on your fit
    dof = len(x) - len(m.params) - 1

    # Calculate the fit errors
    #print 'm.perror is', m.perror
    #print 'm.fnorm is', m.fnorm
    #print 'dof is', dof
    if m.perror==None:
        m.perror=0
    pcerror = m.perror * np.sqrt(m.fnorm / dof)

    # Convert the parameter output to the pars output format from easyfit
    par = [m.params,m.perror,pcerror]
    if(m.status <=0):
        print 'status = ', m.status
    return par


def __create_parser__():
    p = argparse.ArgumentParser(
        description="Python script for calculating the p-values generated by empirical JTK_CYCLE with asymmetry search, which is described in Hutchison, Maienschein-Cline, and Chiang et al. Improved statistical methods enable greater sensitivity in rhythm detection for genome-wide data, PLoS Computational Biology 2015 11(3): e1004094. This script was written by Alan L. Hutchison, alanlhutchison@uchicago.edu, Aaron R. Dinner Group, University of Chicago.",
        epilog="Please contact the correpsonding author if you have any questions.",
        version=VERSION
        )

    analysis = p.add_argument_group(title="CalcP analysis options")

    analysis.add_argument("-f", "--filename",
                   dest="filename",
                   action='store',
                   metavar="filename string",
                   type=str,
                    help='An output file from eJTK.py containing the time series analyzed')
    
    analysis.add_argument("-n", "--null",
                   dest="null",
                   action='store',
                   metavar="null filename string",
                   type=str,
                    help='An output file from eJTK.py which is generated from Gaussian noise with the same header (time points) as the time series analyzed and input with the -f flag')

    analysis.add_argument("-t","--fit",
                          dest="fit",
                          action='store_true',
                          default=False,
                          help="A flag without arguments indicating that the p-value calculation should use the fitting method gauranteed to produce conservative results. THIS FITTING HAS BEEN SHOWN TO BE UNSTABLE IN CERTAIN SITUATIONS AND WILL BE FIXED IN AN UPCOMING VERSION..")
    return p


if __name__=="__main__":
    parser = __create_parser__()
    args = parser.parse_args()
    main(args)

