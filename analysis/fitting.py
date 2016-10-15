"""
    This file is part of pi3diamond, a toolkit for
    confocal scanning, anti-bunching, FLIM, pulsed ODMR / NMR,
    and more sophisticated quantum physics experiments,
    typically performed with NV centers in diamond,
    written in python using the enthought traits packages.

    pi3diamond is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    pi3diamond is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with diamond. If not, see <http://www.gnu.org/licenses/>.

    Copyright (C) 2009-2016 Helmut Fedder <helmut@fedder.net>
"""

import numpy as np
from scipy.optimize import leastsq
from scipy.special import gammaincc
from scipy.stats import mode

########################################################
# utility functions 
########################################################

def baseline(y,n=None):
    """
    Returns the baseline of 'y'. 'n' controls the discretization.
    The difference between the maximum and the minimum of y is discretized into 'n' steps.
    """
    if not n: # estimate a useful number for the histogram bins from shot noise
        y_min = y.min()
        y_max = y.max()
        y_typ = np.max(np.abs((y_min,y_max)))
        n = (y_max-y_min)/y_typ**0.5
    hist, bin_edges = np.histogram(y,int(n))
    return bin_edges[hist.argmax()]

def find_edge(y,bins=20):
    """Returns edge of a step function"""
    h,b=np.histogram(y,bins=bins)
    i0 = bins/2
    i  = h[i0:].argmax()+i0
    threshold = 0.5*(b[0]+b[i])
    return np.where(y>threshold)[0][0]

def run_sum(y, n=10):
    """Calculates the running sum over 'y' (1D array) in a window with 'n' samples."""

    N = len(y)
    
    yp = np.empty(N)

    for i in range(N):
        if i+n > N:
            yp[i]=yp[N-n] # pad the last array entries with the last real entry
        else:
            yp[i]=np.sum(y[i:i+n])
            
    return yp

########################################################
# non-linear least square fitting 
########################################################

def fit(x, y, model, estimator):
    """Perform least-squares fit of two dimensional data (x,y) to model 'Model' using Levenberg-Marquardt algorithm.\n
    'Model' is a callable that takes as an argument the model parameters and returns a function representing the model.\n
    'Estimator' can either be an N-tuple containing a starting guess of the fit parameters, or a callable that returns a respective N-tuple for given x and y."""
    if callable(estimator):
        #return leastsq(lambda pp: model(*pp)(x) - y, estimator(x,y), warning=False)[0]
        return leastsq(lambda pp: model(*pp)(x) - y, estimator(x,y))[0]
    else:
        #return leastsq(lambda pp: model(*pp)(x) - y, estimator, warning=False)[0]
        return leastsq(lambda pp: model(*pp)(x) - y, estimator)[0]

def nonlinear_model(x, y, s, model, estimator, message=False):
    """Performs a non-linear least-squares fit of two dimensional data and a primitive error analysis. 
    
    parameters:
    
    x         = x-data
    y         = y-data
    s         = standard deviation of y
    model     = the model to use for the fit. must be a factory function
                that takes as parameters the parameters to fit and returns
                a function y(x)
    estimator = either an n-tuple (or array) containing the starting guess
                of the fit parameters or a callable that takes x and y
                as arguments and returns a starting guess

    return values:
    
    p        = set of parameters that minimizes the chisqr
    
    cov      = covariance matrix
    
    q        = probability of obtaining a chisqr larger than the observed one
    
               if 0.9 > q > 0.1 the fit is credible
               
               if q > 0.001, the fit may be credible if we expect that the
               reason for the small q are non-normal distributed errors
               
               if q < 0.001, the fit must be questioned. Possible causes are
                   (i) the model is not suitable
                   (ii) the standard deviations s are underestimated
                   (iii) the standard deviations s are not normal distributed
                   
               if q > 0.9, the fit must be questioned. Possible causes are
                   (i) the standard deviations are overestimated
                   (ii) the data has been manipulated to fit the model 
    
    chisqr0  = sum over chisqr evaluated at the minimum
    """
    chisqr = lambda p: ( model(*p)(x) - y ) / s
    if callable(estimator):
        p = estimator(x,y)
    else:
        p = estimator
    result = leastsq(chisqr, p, full_output=True)
    
    if message:
        print result[4], result[3]
    p = result[0]
    cov = result[1]
    
    # there are some cases where leastsq doesn't raise an exception, however returns None for
    # the covariance matrix. To prevent 'invalid index' errors in functions that call nonlinear_model,
    # we replace the 'None' by a matrix with right dimension filled with np.NaN.
    if cov is None:
        cov = np.NaN * np.empty( (len(p),len(p)) )
    
    chi0 = result[2]['fvec']
    
    chisqr0 = np.sum(chi0**2)
    nu = len(x) - len(p)
    
    q = gammaincc(0.5*nu,0.5*chisqr0)
    
    return p, cov, q, chisqr0

########################################################
# standard factory function for non-linear fitting 
########################################################

def Cosinus(a, T, c):
    """Returns a Cosinus function.
    
        f = a\cos(2\pi(x-x0)/T)+c
    
    Parameter:
    
    a    = amplitude
    T    = period
    x0   = position
    c    = offset in y-direction
    """
    return lambda x: a*np.cos( 2*np.pi*x/float(T) ) + c

setattr(Cosinus, 'formula', r'$cos(c,a,T;x)=a\cos(2\pi x/T)+c$')

def CosinusEstimator(x, y):
    c = y.mean()
    a = 2**0.5 * np.sqrt( ((y-c)**2).sum() )
    # better to do estimation of period from
    Y = np.fft.fft(y)
    N = len(Y)
    D = float(x[1] - x[0])
    i = abs(Y[1:N/2+1]).argmax()+1
    T = (N * D) / i
    return a, T, c

def CosinusNoOffset(a, T):
    """Returns a Cosinus function without constant offset.
    
        f = a\cos(2\pi(x-x0)/T)
    
    Parameter:
    
    a    = amplitude
    T    = period
    x0   = position
    """
    return lambda x: a*np.cos( 2*np.pi*x/float(T) )

setattr(CosinusNoOffset, 'formula', r'$cos(a,T;x)=a\cos(2\pi x/T)$')

def CosinusNoOffsetEstimator(x, y):
    a = 2**0.5 * np.sqrt( (y**2).sum() )
    # better to do estimation of period from
    Y = np.fft.fft(y)
    N = len(Y)
    D = float(x[1] - x[0])
    i = abs(Y[1:N/2+1]).argmax()+1
    T = (N * D) / i
    return a, T

def ExponentialZero(a, w, c):
    """Exponential centered at zero.
    
        f = a*exp(-x/w) + c
    
    Parameter:
    
    a    = amplitude
    w    = width
    c    = offset in y-direction
    """
    return lambda x: a*np.exp(-x/w)+c

def ExponentialZeroEstimator(x, y): 
    """Exponential Estimator without offset. a*exp(-x/w) + c"""
    c=y[-1]
    a=y[0]-c
    w=x[-1]*0.5
    return a, w, c

def GaussianZero(a, w, c):
    """Gaussian function centered at zero.
    
        f = a*exp(-(x/w)**2) + c
    
    Parameter:
    
    a    = amplitude
    w    = width
    c    = offset in y-direction
    """
    return lambda x: a*np.exp( -(x/w)**2 ) + c

setattr(GaussianZero, 'formula', r'$f(a,w,c;x)=a\exp(-(x/w)^2)+c$')

def GaussianZeroEstimator(x, y): 
    """Estimator for GaussianZero: a*exp(-0.5*(x/w)**2) + c"""
    c=y[-1]
    a=y[0]-c
    w=x[-1]*0.5
    return a, w, c

def Gaussian(c, a, x0, w):
    """Gaussian function.
    
        f = a*exp( -0.5(x-x0)**2 / w**2 ) + c
    
    Parameter:
    
    a    = amplitude
    w    = width
    c    = offset in y-direction
    """
    return lambda x: c + a*np.exp( -0.5*((x-x0)/w)**2   )

setattr(Gaussian, 'formula', r'$f(c,a,x0,w;x)=c+a\exp(-0.5(x-x0)^2/w^2)$')

def ExponentialPowerZero(a, w, p, c):
    """Exponential decay with variable power centered at zero.
    
        f = a*exp(-(x/w)**p) + c
    
    Parameter:
    
    a    = amplitude
    w    = width
    p    = power
    c    = offset in y-direction
    """
    return lambda x: a*np.exp( -(x/w)**p ) + c

setattr(ExponentialPowerZero, 'formula', r'$f(a,w,p,c;x)=a\exp(-(x/w)^p)+c$')

def ExponentialPowerZeroEstimator(x, y): 
    """Estimator for exponential decay with variable offset."""
    c=y[-1]
    a=y[0]-c
    w=x[-1]*0.5
    return a, w, 2, c


def GaussianZeroEstimator(x, y): 
    """Gaussian Estimator without x offset. c+ a*exp( -0.5*(x/w)**2)"""
    a=y.argmax()
    #x0=x[y.argmax()]
    w=x[(len(x)/2)]
    c=(min(y)+max(y))/2
    return a, w, c


def DoubleGaussian(a1, a2, x01, x02, w1, w2):
    """Gaussian function with offset."""
    return lambda x: a1*np.exp( -0.5*((x-x01)/w1)**2   ) + a2*np.exp( -0.5*((x-x02)/w2)**2   )

setattr(DoubleGaussian, 'formula', r'$f(c,a1, a2,x01, x02,w1,w2;x)=a_1\exp(-0.5((x-x_{01})/w_1)^2)+a_2\exp(-0.5((x-x_{02})/w_2)^2)$')

def DoubleGaussianEstimator(x, y):
    center = (x*y).sum() / y.sum()
    ylow = y[x < center]
    yhigh = y[x > center]
    x01 = x[ylow.argmax()]
    x02 = x[len(ylow)+yhigh.argmax()]
    a1 = ylow.max()
    a2 = yhigh.max()
    w1 = w2 = center**0.5
    return a1, a2, x01, x02, w1, w2

# important note: lorentzian can also be parametrized with an a' instead of a,
# such that a' is directly related to the amplitude (a'=f(x=x0)). In this case a'=a/(pi*g)
# and f = a * g**2 / ( (x-x0)**2 + g**2 ) + c.
# However, this results in much poorer fitting success. Probably the g**2 in the numerator
# causes problems in Levenberg-Marquardt algorithm when derivatives
# w.r.t the parameters are evaluated. Therefore it is strongly recommended
# to stick to the parametrization given below.
# The amplitude is a/(pi*g), the area under the curve is 'a'
def Lorentzian(c, x0, g, a):
    """Lorentzian centered at x0, with area a, offset y0 and HWHM g."""
    return lambda x: a / np.pi * (  g / ( (x-x0)**2 + g**2 )  ) + c

setattr(Lorentzian, 'formula', r'$f(x0,g,a,c;x)=a/\pi (g/((x-x_0)^2+g^2)) + c$')



def LorentzianNoOffset(x0, g, a):
    """Lorentzian centered at x0, with amplitude a, and HWHM g."""
    return lambda x: a / np.pi * (  g / ( (x-x0)**2 + g**2 )  )

def Nlorentzians(*p):
    N = (len(p)-1)/3
    def f(x):
        y = p[0]*np.ones(x.shape)
        i = 0
        for i in range(N):
            y += LorentzianNoOffset(*p[i*3+1:i*3+4])(x)
        return y   
    return f

def LorentzianEstimator(x, y):
    c = mode(y)[0][0]
    yp = y - c
    Y = np.sum(yp) * (x[-1] - x[0]) / len(x)
    ymin = yp.min()
    ymax = yp.max()
    if ymax > abs(ymin):
        y0 = ymax
    else:
        y0 = ymin
    x0 = x[y.argmin()]
    g = Y / (np.pi * y0)
    a = y0 * np.pi * g
    return x0, g, a, c

def Antibunching(alpha, c, tau, t0):
    """Antibunching. g(2) accounting for Poissonian background."""
    return lambda t: c*(1-alpha*np.exp(-(t-t0)/tau))

setattr(Antibunching, 'formula', r'$g(\alpha,c,\tau,t_0;t)=c(1 - \alpha \exp(-(t-t_0)/\tau))$')

def FCSTranslationRotation(alpha, tau_r, tau_t, N):
    """Fluorescence Correlation Spectroscopy. g(2) accounting for translational and rotational diffusion."""
    return lambda t: (1 + alpha*np.exp(-t/tau_r) ) / (N * (1 + t/tau_t) )

setattr(FCSTranslationRotation, 'formula', r'$g(\alpha,\tau_R,\tau_T,N;t)=\frac{1 + \alpha \exp(-t/\tau_R)}{N (1 + t/\tau_T)}$')

def FCSTranslation(tau, N):
    """Fluorescence Correlation Spectroscopy. g(2) accounting for translational diffusion."""
    return lambda t: 1. / (N * (1 + t/tau) )

setattr(FCSTranslation, 'formula', r'$g(\tau,N;t)=\frac{1}{N (1 + t/\tau)}$')

def SumOverFunctions( functions ):
    """Creates a factory that returns a function representing the sum over 'functions'.
    'functions' is a list of functions. 
    The resulting factory takes as arguments the parameters to all functions,
    flattened and in the same order as in 'functions'."""
    def function_factory(*args):
        def f(x):
            y = np.zeros(x.shape)
            i = 0
            for func in functions:
                n = func.func_code.co_argcount
                y += func(*args[i,i+n])(x)
                i += n
        return f
    return function_factory

def brot_transitions_upper(B, D, E, phase):
    return lambda theta: 3./2. * B**2/D * np.sin(theta + phase)**2 + ( B**2 * np.cos(theta + phase)**2 + (E + B**2/(2*D) * np.sin(theta+phase)**2)**2)**0.5 + D
    
def brot_transitions_lower(B, D, E, phase):
    return lambda theta: 3./2. * B**2/D * np.sin(theta + phase)**2 - ( B**2 * np.cos(theta + phase)**2 + (E + B**2/(2*D) * np.sin(theta+phase)**2)**2)**0.5 + D


#################################################################
# convenience functions for performing some frequently used fits
#################################################################

from scipy.signal import find_peaks_cwt

def find_peaks(x, y, width, n_peaks=-1, baseline_bins=None, estimator='wavelet', peak_shape='Lorentzian'):
    """
    Find peaks in a noisy 1D data set by applying continuous wavelet transform
    with an expected peak width and fit the determined peaks with lorentzians.
    
    The function determines automatically whether the peaks are positive or negative,
    however all peaks are expected to have the same sign.
    
    The number of peaks can be limited by specifying n_peaks > 0. In this case, the 'n_peaks'
    peaks with largest amplitudes are taken.
    
    By default, the peaks are subsequently least-square-fitted with a multi-lorentzian function.
    This step can be omitted by specifying 'peak_shape'=None. 
    """
    
    # estimate the baseline
    y0 = baseline(y, baseline_bins)
    
    # determine whether extrema are positive or negative by checking the distance of the absolute maximum and absolute minimum w.r.t the baseline
    if np.abs(y.max()-y0) > np.abs(y0-y.min()): # highest maximum larger than smallest minimum
        yp = y - y0
        sign = 1
    else:
        yp = y0 - y
        sign = -1
    
    if estimator == 'wavelet':
        dx = x[1]-x[0]
        peak_indices = np.array(find_peaks_cwt(yp, np.array((width/dx,))))
        peak_amps = yp[peak_indices]
        if n_peaks > 0: # keep only the n_peaks largest
            sort_map = peak_amps.argsort()
            peak_amps = peak_amps[sort_map][-n_peaks:]
            peak_indices = peak_indices[sort_map][-n_peaks:]
    else:
        raise ValueError('Estimator not implemented')
    
    res = {'x0':x[peak_indices], 'y0':y[peak_indices]}
    
    n = len(peak_indices)
    if peak_shape == 'Lorentzian':
        hwhm = 0.5*width
        p = [0.0]
        for i, peak_index in enumerate(peak_indices):
            p.append(x[peak_index])
            p.append(hwhm)
            p.append(peak_amps[i]*np.pi*hwhm)
        r = fit_n_peaks(x,yp,p,LorentzianNoOffset)
        if not (r[-1] == 0):
            p = np.array(r[0])
            #delta = np.diag(r[1])**0.5
            if sign == -1:
                p[0] *= -1
                p[3::3] *= -1
            p[0] += y0
            res['p']=p
    return res

def fit_n_peaks(x,y,p,peak_func):

    N = (len(p)-1)/3

    # chi for N peaks with a common baseline
    def chi(p):
        yp = p[0]-y
        for i in range(N):
            yp += peak_func(*p[i*3+1:i*3+4])(x)
        return yp

    r = leastsq(chi, p, full_output=True)

    return r

    

def find_local_maxima(y,n):
    "Returns the indices of the n largest local maxima of y."

    half = 0.5*y.max()
    mask = y>half
    
    # get left and right edges of connected regions
    
    right_shifted = np.append(False, mask[:-1])
    left_shifted = np.append(mask[1:], False)
    
    left_edges =  np.where( np.logical_and(mask,np.logical_not(right_shifted) ))[0]
    right_edges = np.where( np.logical_and(mask,np.logical_not(left_shifted)) )[0] + 1

    if len(left_edges) < n:
        raise RuntimeError('did not find enough edges')
    
    indices = []
    for k in range(len(left_edges)):
        left = left_edges[k]
        right = right_edges[k]
        indices.append( y[left:right].argmax()+left )
    indices = np.array(indices)
    maxima = y[indices]
    indices = indices[maxima.argsort()][::-1]
    return indices[:n]


"""
def fit_rabi(x, y, s):
    y_offset=y.mean()
    yp = y - y_offset

    p = fit(x, yp, CosinusNoOffset, CosinusNoOffsetEstimator)
    if p[0] < 0:
        p[0] = -p[0]
        p[2] =  ( ( p[2]/p[1] + 0.5 ) % 1 ) * p[1]
        p = fit(x, yp, CosinusNoOffset, p)
    p = (p[0], p[1], p[2], y_offset)
    result = nonlinear_model(x, y, s, Cosinus, p)
    p = result[0]
    if p[2]>0.5*p[1]:
        while(p[2]>0.5*p[1]):
            p[2] -= p[1]
        result = nonlinear_model(x, y, s, Cosinus, p)
    return result
"""

def fit_rabi(x, y, s):
    y_offset=y.mean()
    yp = y - y_offset

    p = fit(x, yp, CosinusNoOffset, CosinusNoOffsetEstimator)
    if p[0] < 0:
        p[0] = -p[0]
        p[2] =  ( ( p[2]/p[1] + 0.5 ) % 1 ) * p[1]
        #p = fit(x, yp, CosinusNoOffset, p)
    p = (p[0], p[1], y_offset)
    return nonlinear_model(x, y, s, Cosinus, p)

def extract_pulses(y):
    """
    Extracts pi, pi/2 and 3pi/2 pulses from a Rabi measurement.
    
    Parameters:
        y = the arry containing y data
        
    Returns:
        f, r, pi, 2pi = arrays containing the indices of the respective pulses and their multiples 
    """
    # The goal is to find local the rising and falling edges and local minima and maxima.
    # First we estimate the 'middle line' by the absolute minimum and maximum.
    # Then we cut the data into sections below and above the middle line.
    # For every section we compute the minimum, respectively maximum.
    # The falling and rising edges mark multiples of pi/2, respectively 3pi/2 pulses.

    # center line
    m=0.5*(y.max()+y.min())
    
    # boolean array containing positive and negative sections
    b = y < m
    
    # indices of rising and falling edges
    # rising edges: last point below center line
    # falling edges: last point above center line
    rising = np.where(b[:-1]&~b[1:])[0]
    falling = np.where(b[1:]&~b[:-1])[0]

    # local minima and maxima
    
    pi = [ y[:rising[0]].argmin() ]
    two_pi = [ y[:falling[0]].argmax() ]
    
    for i in range(1,len(rising)):
        pi.append( rising[i-1] + y[rising[i-1]:rising[i]].argmin() )
        
    for i in range(1,len(falling)):
        two_pi.append(falling[i-1] + y[falling[i-1]:falling[i]].argmax() )

    # For rising edged, we always use the last point below the center line,
    # however due to finite sampling and shot noise, sometimes
    # the first point above the line may be closer to the actual zero crossing
    for i, edge in enumerate(rising):
        if y[edge+1]-m < m-y[edge]:
            rising[i] += 1
    # similarly for the falling edges
    for i, edge in enumerate(falling):
        if m-y[edge+1] < y[edge]-m:
            falling[i] += 1
    
    return np.array(falling), np.array(rising), np.array(pi), np.array(two_pi)


if __name__ == '__main__':
    
    from tools.data_toolbox import load
    
    filename = '2014-11-14_ODMR_08.pys'
    
    d = load(filename)
    
    x = d['frequency']
    y = d['counts']

    #y0 = baseline(y)
    
    #y = y0 - y

    width = 5e5

    r = find_peaks(x,y,width,n_peaks=3)
    
    x0 = r['x0']
    y0 = r['y0']
    
    p = r['p']

    #hwhm = 0.5*width
    
    #p = [0.0]
    #for i, xi in enumerate(x0):
    #    p.append(xi)
    #    p.append(hwhm)
    #    p.append(y0[i]*np.pi*hwhm)
    #r = fit_n_peaks(x,y,p,LorentzianNoOffset)
    
    #pp = r[0]
    
    import pylab
    pylab.figure()
    pylab.plot(x,y,'b-')
    pylab.plot(x0,y0,'r.')
    for i in range(3):
        pi = np.append(p[0],p[i*3+1:i*3+4])
        pylab.plot(x,Lorentzian(*pi)(x),'r-')
    #pylab.plot(x,n_lorentzians(*p)(x),'r.')
    #pylab.plot(x,n_lorentzians(*pp)(x),'g.')
    pylab.show()
    
