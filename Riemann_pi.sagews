# Author: Daniel Hutama
# This is a SageMath file. It will probably work as an .ipynb with some modifications like changing ^ to **, etc. 

# This block creates 3 functions of interest: 
# The first is the main term in Riemann's explicit formula.
# The second is the correction term involving the zeta zeros (using an asymptotic approximation for li(x).
# The third is the explicit formula to count primes up to x. The inputs are described in the readme. 
# Note there is a tradeoff between accuracy and computation time.

# To plot pi(x), the prime step function, call something like: 
# primes = plot(prime_pi(x), (x, 0, plotdistance), rgbcolor='red')
# To plot the explicit formula simply call something like:
# Rx = plot(Riemann_explicit(x, numsums, gamma, numasymterms), (x, 2, plotdistance), rgbcolor = 'blue',xmin=0, xmax=plotdistance, ymax=prime_pi(plotdistance))
# To plot both together use show(primes + Rx)


# not all of these imports are used. 
import numpy as np
import cmath
import matplotlib.pyplot as plt
import math
import cmath
import mpmath as mp
from numpy import polyfit
plt.rcParams['text.usetex'] = True
plt.rcParams['text.latex.unicode'] = True
import scipy

from sage.functions.spike_function import SpikeFunction
from scipy.interpolate import interp1d
from scipy import stats
import scipy.integrate as integrate
import scipy.special as special
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
from mpl_toolkits.axes_grid1.inset_locator import mark_inset

# pi(x) using the asymptotic expansion of li(x)
%time
def Riemann_main(x, numbersums, numberzeros, numberasymterms):
    return sum([sum([moebius(n)/n*li(x^(1/n)) for n in xrange(1, numbersums+1)]), -1/log(x), (1/pi)*arctan(pi/log(x))])

# includes first k+1 summation terms of the asymptotic expansion of li(x)~x/log(x)
def Riemann_correction(x, numbersums, numberzeros, numberasymterms):
    zz = zeta_zeros()[0:numberzeros]
    n_dict = {n:[] for n in xrange(1, numbersums+1)}
    for n in xrange(1, numbersums+1):
        asymterms=[1]
        for k in xrange(1, numberasymterms+1):
            for gamma in zz:
                asymterms.append((factorial(k)*(2*n)^k)/(log(x)^k*(1+4*gamma^2)^k))
        sumasymterms=sum(asymterms)
        if moebius(n)==1:
            n_dict[n] = -sum([((4*x^(1/(2*n))/log(x)*(cos((gamma/n)*log(x))+2*gamma*sin((gamma/n)*log(x)))/(1+4*gamma^2))*sumasymterms) for gamma in zz])
        elif moebius(n)==0:
            n_dict[n] = 0
        elif moebius(n)==-1:
            n_dict[n] = sum([((4*x^(1/(2*n))/log(x)*(cos((gamma/n)*log(x))+2*gamma*sin((gamma/n)*log(x)))/(1+4*gamma^2))*sumasymterms) for gamma in zz])
    return sum(n_dict.values())

def Riemann_explicit(x, numbersums, numberzeros, numberasymterms):
    return sum([Riemann_main(x, numbersums, numberzeros, numberasymterms), Riemann_correction(x, numbersums, numberzeros, numberasymterms)])
