# -*- coding: utf-8 -*-
'''====
Easyfit
====

Provides an easy to use wrapper to fit common functions to a data set using the
Levenberg–Marquardt algorithm provided by mpfit. A full description of the
suppoerted functions and how to use the wrapper is given in easyfit.fit
'''
import numpy as np
from mpfit import mpfit
import warnings


def line(x,p):
    '''Parameter: Slope, Intercept\n
    Return
    ------
    >>> p[0]*x + p[1]
    '''
    return p[0]*x + p[1]

def line0(x,p):
    '''Parameter: Slope\n
    Return
    ------
    >>> p[0]*x
    '''
    return p[0]*x

def sine(x,p):
    '''Parameter: Scale, Wavelength, Phase, Offset\n
    Return
    ------
    >>> p[0]*np.sin(2*np.pi*x/p[1]+p[2])+p[3]
    '''
    return p[0]*np.sin(2*np.pi*x/p[1]+p[2])+p[3]
    
def fermi(x,p):
    '''Parameter: Scale, Edge Position, Width, Offset\n
    Return
    ------
    >>> p[0]/(np.exp((x-p[1])/p[2])+1)+p[3]
    '''
    return p[0]/(np.exp((x-p[1])/p[2])+1)+p[3]

def gauss(x,p):
    '''Parameter: Scale, Mean, Std, Offset\n
    Return
    ------
    >>> p[0]*np.exp(-0.5*(x-p[1])**2/p[2]**2)+p[3]
    '''
    return p[0]*np.exp(-0.5*(x-p[1])**2/p[2]**2)+p[3]
    
def exp(x,p):
    '''Parameter: Scale, Decay Time, Offset\n
    Return
    ------
    >>> p[0]*np.exp(-x*p[1])+p[2]
    '''
    return p[0]*np.exp(-x*p[1])+p[2]

def poly(x,p,n):
    '''Parameter: Scale of each power from [0..n]\n
    Return
    ------
    >>> Sum[n=0,n=N] p[n]*x**n
    '''
    y = 0
    for i in np.arange(n+1, dtype='float64'):
        try:
            y+=np.power(x,i)*p[i]
        except FloatingPointError:
            print y
    return y

def ipoly(x,p,n):
    '''Parameter: Scale of each power from [0..n]\n
    Return
    ------
    >>> Sum[n=0,n=N] p[n]*x**-n
    '''
    y = 0
    for i in range(n+1):
        y+=np.power(x,-i)*p[i]
    return y
    
def plaw(x,p):
    '''Parameter: Scale, Exponent
    Return
    ------
    >>> p[0]*x**p[1]
    '''
    return p[0]*x**p[1]
    
def fit(typ='line',x=None, y=None, yerr=None,p0=None):
    '''
    Takes the data and performs a least square fit of the specified type.
    
    Parameters
    ----------
    typ : string
        Predefined function that will be fitted to the data. You can find a
        list of all supported functions below.
    x : array_like or None
        X data. If None is given a fit will be performed, yet it is based on
        an internally created x data set that runs from [0,N] where N is the
        number of y data points provided. Thus all parameters that are not
        independent of your choice of x, e.g. slope, are not to be trusted!
        If you are only interested in parameters that are independent of x such
        as the heigth of a gaussian you'll probably get away without providing
        an adequate set of x data.
    y : array_like
        y data. You have to provide a y array. Otherwise there is nothing to
        fit.
    yerr : array_like or None
        Error in y direction. If None is given the fit will assume a uniform
        weight of 1.
    p0 : array_like or None
        Initial guess of fit parameters. If p0 is None all parameters are
        initalized to one or zero depending on the meaning of the individual
        parameter.
    
    Returns
    -------
    x2 : float
        Reducd chi-square.
    pars : array_like
        Fit parameters returned by mpfit. The meaning of the subarrays are:\n
        pars[0]\tBest fit parameters\n
        pars[1]\tFit errors\n
        pars[2]\tProperly scaled errors\n
        Note that it is assumed that the chi-squared returned is sufficiently
        good to justify the scaling of the fit erros. It is pars[2] = pars[1]*
        sqrt(x2)
    xfit,yfit : array_like
        x and y data that can directly be used within matplotlib or another 
        comparable plotting library to display the fit.
        
    Available functions/fits
    ------------------------
    line
        Straight line, parameters: Slope and intercept\n
        >>>  p[0]*x + p[1]
    line0
        Straight line with designated zero crossing, parameters: Slope\n
        >>> p[0]*x
    sine
        Sine, parameters: Scaling, Wavelength, Phase, Offset\n
        >>> p[0]*sin(2*Pi*x/p[1]+p[2])+p[3]
    fermi
        Fermifunction, parameters: Scaling, Edge Position, Width, Offset\n 
        >>> p[0]/(exp((x-p[1])/p[2])+1)+p[3]
    gauss
        Gaussian, parameters: Scaling, Mean, Std, Offset\n
        >>> p[0]*exp(-0.5*(x-p[1])**2/p[2]**2)+p[3]
    exp
        Exponential, parameters: Scaling, Inverse Decaytime, Offset\n
        >>> p[0]*exp(-x*p[1])+p[2]
    plaw
        Power law, parameters: Scaling, Power\n
        >>> p[0]*x**p[1]
    polyN
        Polynomial of order N. Usage: poly3, poly5, poly10, etc. Parameters: 
        Scalings\n
        >>> Sum[n=0,n=N] p[n]*x**n
    ipolyN
        Inverse polynomial of order N. Usage: ipoly3, poly5, poly10 etc. 
        Parameters: Scalings\n
        >>> Sum[n=0,n=N] p[n]*x**-n
        
    Example
    -------
    
    The following code snippet explains the use of the easyfit wrapper
    
    >>> import matplotlib.pylab as plt
    >>> import numpy as np
    >>> import easyfit as ef
    >>>
    >>> x = np.linspace(0,100,30)
    >>> y = 0.05*x + 2*(np.random.rand(30)-0.5)
    >>>
    >>> p0 = [1]
    >>> x2, pars, xfit, yfit = ef.fit('line0',x,y,None,p0)
    >>>
    >>> plt.scatter(x,y)
    >>> plt.plot(xfit,yfit)
    >>> plt.show()
    
    '''
    #=========================================================================#
    #                   Filter Future Warning From Numpy
    #=========================================================================#
    warnings.filterwarnings("ignore",category=FutureWarning)
    
    #=========================================================================#
    #                           Set default arrays
    #=========================================================================#
    n=0
    if 'ipoly' in typ:
        n = int(typ[5:])
        typ = 'ipoly'
    elif 'poly' in typ:
        n = int(typ[4:])
        typ = 'poly'
    if x==None: x = np.arange(len(y))
    if yerr==None: yerr = np.ones(len(y))
    if p0 == None: 
        if typ == 'line':    p0 = [1,0]
        elif typ == 'line0': p0 = [1]
        elif typ == 'sine':  p0 = [1,1,0,0]
        elif typ == 'fermi': p0 = [1,1,1,0]
        elif typ == 'gauss': p0 = [1,0,1,0]
        elif typ == 'exp':   p0 = [1,1,0]
        elif typ == 'plaw':  p0 = [1,1,0]
        elif typ == 'poly' or typ == 'ipoly':
            p0 = [1]*(n+1)
            
    #=========================================================================#
    #               Ensure that all given arrays are numpy arrays
    #=========================================================================#            
    x = np.array(x)
    y = np.array(y)
    yerr = np.array(yerr)
    p0 = np.array(p0)

    #=========================================================================#
    #                       Setup proper fit function
    #=========================================================================#
    models = {'line': line,
              'line0': line0,
              'sine': sine,
              'fermi': fermi,
              'gauss': gauss,
              'exp': exp,
              'plaw': plaw,
              'poly': lambda x,p: poly(x,p,n),
              'ipoly': lambda x,p: ipoly(x,p,n)}

              
    def fitfunc(p, fjac=None, x=None, y=None, err=None):
        model = models[typ](x,p)
        status = 0
        return [status, (y-model)/err]

    #=========================================================================#
    #  Initialize fit info dictionary and try to fit function to data
    #=========================================================================#    
    parbase = {'value':0,'fixed':0,'limited':[0,0],'limits':[0.,0.]}
    parinfo = [parbase]*len(p0)
    for i in range(len(p0)):
        parinfo[i]['value'] = p0[i]
    fa = {'x': x, 'y': y, 'err': yerr}
    m = mpfit(fitfunc, p0, parinfo=parinfo, functkw=fa,quiet=1)
    dof = len(x) - len(m.params)
    pcerror = m.perror * np.sqrt(m.fnorm / dof)
    par = [m.params,m.perror,pcerror]
    if(m.status <=0):
        print 'status = ', m.status

    #=========================================================================#
    #    Calculate goodness of fit and an easy to plot fitted data set
    #=========================================================================#
    x2 = m.fnorm/dof
    xfit = np.linspace(np.min(x),np.max(x),1000)
    yfit = models[typ](xfit,par[0])
    
    return x2,par,xfit,yfit
    
    
def arbFit(fct=line,x=None, y=None, yerr=None,p0=None):
    '''
    Takes the data and performs a least square fit of the specified type.
    
    Parameters
    ----------
    typ : function
        User defined function that will be fitted to the data. Has to obey the
        following convention for its arguments: F(x,p)
    x : array_like or None
        X data. If None is given a fit will be performed, yet it is based on
        an internally created x data set that runs from [0,N] where N is the
        number of y data points provided. Thus all parameters that are not
        independent of your choice of x, e.g. slope, are not to be trusted!
        If you are only interested in parameters that are independent of x such
        as the heigth of a gaussian you'll probably get away without providing
        an adequate set of x data.
    y : array_like
        y data. You have to provide a y array. Otherwise there is nothing to
        fit.
    yerr : array_like or None
        Error in y direction. If None is given the fit will assume a uniform
        weight of 1.
    p0 : array_like or None
        Initial guess of fit parameters. If p0 is None all parameters are
        initalized to one or zero depending on the meaning of the individual
        parameter.
    
    Returns
    -------
    x2 : float
        Reducd chi-square.
    pars : array_like
        Fit parameters returned by mpfit. The meaning of the subarrays are:\n
        pars[0]\tBest fit parameters\n
        pars[1]\tFit errors\n
        pars[2]\tProperly scaled errors\n
        Note that it is assumed that the chi-squared returned is sufficiently
        good to justify the scaling of the fit erros. It is pars[2] = pars[1]*
        sqrt(x2)
    xfit,yfit : array_like
        x and y data that can directly be used within matplotlib or another 
        comparable plotting library to display the fit.
    Example
    -------
    
    The following code snippet explains the use of the easyfit wrapper
    
    >>> import matplotlib.pylab as plt
    >>> import numpy as np
    >>> import easyfit as ef
    >>>
    >>> def userFct(x,p):
    >>>     return p[0]*x**2 + np.exp(-p[1]*x)
    >>>
    >>> x = np.linspace(0,100,30)
    >>> y = userFct(x,[-0.5,0.25]) + 100*(2*np.random.rand(30)-1)
    >>>
    >>> p0 = [1,0]
    >>> x2, pars, xfit, yfit = ef.arbFit(userFct,x,y,None,p0)
    >>>
    >>> plt.scatter(x,y)
    >>> plt.plot(xfit,yfit)
    >>> plt.show()
    
    '''
    #=========================================================================#
    #                   Filter Future Warning From Numpy
    #=========================================================================#
    warnings.filterwarnings("ignore",category=FutureWarning)
    
    #=========================================================================#
    #                           Set default arrays
    #=========================================================================#
    if x==None: x = np.arange(len(y))
    if yerr==None: yerr = np.ones(len(y))
    if p0 == None: 
        p0 = np.ones(100)
            
    #=========================================================================#
    #               Ensure that all given arrays are numpy arrays
    #=========================================================================#            
    x = np.array(x)
    y = np.array(y)
    yerr = np.array(yerr)
    p0 = np.array(p0)

    #=========================================================================#
    #                       Setup proper fit function
    #=========================================================================#
    def fitfunc(p, fjac=None, x=None, y=None, err=None):
        model = fct(x,p)
        status = 0
        return [status, (y-model)/err]

    #=========================================================================#
    #  Initialize fit info dictionary and try to fit function to data
    #=========================================================================#    
    parbase = {'value':0,'fixed':0,'limited':[0,0],'limits':[0.,0.]}
    parinfo = [parbase]*len(p0)
    for i in range(len(p0)):
        parinfo[i]['value'] = p0[i]
    fa = {'x': x, 'y': y, 'err': yerr}
    m = mpfit(fitfunc, p0, parinfo=parinfo, functkw=fa,quiet=1)
    dof = len(x) - len(m.params)
    if m.perror==None:
        pcerror=None
    else:
        pcerror = m.perror * np.sqrt(m.fnorm / dof)
    par = [m.params,m.perror,pcerror]
    if(m.status <=0):
        print 'status = ', m.status

    #=========================================================================#
    #    Calculate goodness of fit and an easy to plot fitted data set
    #=========================================================================#
    x2 = m.fnorm/dof
    xfit = np.linspace(np.min(x),np.max(x),1000)
    yfit = fct(xfit,par[0])
    
    return x2,par,xfit,yfit