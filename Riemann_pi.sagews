# Author: Daniel Hutama
# This is a SageMath file. It will probably work as an .ipynb with some modifications like changing ^ to **, etc. 

# This block creates 3 functions of interest: 
# The first is the main term in Riemann's explicit formula.
# The second is the correction term involving the zeta zeros (using an asymptotic approximation for li(x)).
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



# pi(x) using R(x)
%time
def Riemann_main(x, numbersums, numberzeros):
    return sum([sum([moebius(n)/n*li(x^(1/n)) for n in xrange(1, numbersums+1)]), -1/log(x), (1/pi)*arctan(pi/log(x))])

def Riemann_correction(x, numbersums, numberzeros):
    zz = zeta_zeros()[0:numberzeros]
    n_dict = {n:[] for n in xrange(1, numbersums+1)}
    for n in xrange(1, numbersums+1):
        if moebius(n)==1:
            n_dict[n] = -sum([(1/n)*2*real(Ei(((1/2)+I*gamma)/n*log(x))) for gamma in zz])
        elif moebius(n)==0:
            n_dict[n] = 0
        elif moebius(n)==-1:
            n_dict[n] = sum([(1/n)*2*real(Ei(((1/2)+I*gamma)/n*log(x))) for gamma in zz])
    return sum(n_dict.values())

def Riemann_explicit(x, numbersums, numberzeros):
    return sum([Riemann_main(x, numbersums, numberzeros), Riemann_correction(x, numbersums, numberzeros)])


# ------ This function saves plots to your folder
def savetofile(x, plotdistance, numsums, minzeros, maxzeros):
    primes = plot(prime_pi(x), (x, 0, plotdistance+1), rgbcolor = 'red')
    for gamma in xrange(minzeros, maxzeros+1):
        myRx = plot(Riemann_explicit(x, numsums, gamma), (x, 2, plotdistance), rgbcolor = 'blue',xmin=0, xmax=plotdistance, ymax=prime_pi(plotdistance))
        if gamma == 1:
            p = plot(primes + myRx + text('{} nontrivial zero'.format(gamma),[10, 22.5], rgbcolor='black'))
            p.save('{}.png'.format(gamma))
        if gamma != 1:
            p = plot(primes + myRx + text('{} nontrivial zeros'.format(gamma),[10, 22.5], rgbcolor='black'))
            p.save('{}.png'.format(gamma))
