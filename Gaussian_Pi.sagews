︠bc2a740d-6aa7-468d-8d2f-1ce544597cd1︠
# Not all imports are used.

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

from sage.libs.lcalc.lcalc_Lfunction import *
︡b880eadd-7e9e-4cf0-b223-ea5f34f79dcc︡{"done":true}︡
︠afeda72e-046b-4d3e-be3d-af524b49c9b9s︠
# These functions count the exact value of pi_G(x)

def gi_of_norm(max_norm):
    Gaussian_primes = {}
    Gaussian_integers = {}
    Gaussian_integers[0] = [(0,0)]
    for x in range(1, ceil(sqrt(max_norm))):
        for y in range(0, ceil(sqrt(max_norm - x^2))):
            N = x^2 + y^2
            if Gaussian_integers.has_key(N):
                Gaussian_integers[N].append((x,y))
            else:
                Gaussian_integers[N] = [(x,y)]
            if(y == 0 and is_prime(x) and x%4==3):
                have_prime = True
            elif is_prime(N) and N%4==1 or N==2:
                have_prime = True
            else:
                have_prime =False
            if have_prime:
                if Gaussian_primes.has_key(N):
                    Gaussian_primes[N].append((x,y))
                else:
                    Gaussian_primes[N] = [(x,y)]
    return Gaussian_primes,Gaussian_integers

def all_gaussian_primes_up_to_norm(N):
    gips = gi_of_norm(N+1)[0]
    return flatten([uniq([(x,y), (-y,x), (y,-x), (-x,-y)]) for x,y in flatten(gips.values(), max_level=1)], max_level=1)

def all_gaussian_integers_up_to_norm(N):
    gis = gi_of_norm(N+1)[1]
    return flatten([uniq([(x,y), (-y,x), (y,-x), (-x,-y)]) for x,y in flatten(gis.values(), max_level=1)], max_level=1)

# returns a list of the cumulative value of Gaussian primes up to a given norm (unique up to multiplcation by units in Z[i])
def prime_counter(maxnorm):
    myrange=range(0,maxnorm+1)
    emptyset = []
    for mynorm in myrange:
        emptyset.append(len(all_gaussian_primes_up_to_norm(mynorm))/4)
    return emptyset
︡9dcdb143-319c-496e-83a3-f2c895f15b40︡{"done":true}︡
︠b05ad548-809d-4d51-a8a1-0f4e2d347be7s︠
# These lines create a list of L function's non-trivial zeros
chi_mod4=DirichletGroup(4).gen() 
L_mod4=Lfunction_from_character(chi_mod4, type="int")
L_mod4_zeros_list = L_mod4.find_zeros(0,1000,.1)

#sorted list of zeros of zeta and L with imaginary parts less than 1000
all_zeros = sorted(L_mod4_zeros_list + zeta_zeros()[0:648])

x = var('x')
def Riemann_main(x, numbersums, numberzeros, numbertrivials):
    return sum([moebius(n)/n*Ei((1/n)*log(x)) for n in xrange(1, numbersums+1)])

def nontrivial_correction(x, numbersums, numberzeros, numbertrivials):
    zz = all_zeros[0:numberzeros]
    n_dict = {n:[] for n in xrange(1, numbersums+1)}
    for n in xrange(1, numbersums+1):
        if moebius(n)==1:
            n_dict[n] = -sum([(1/n)*2*real(Ei((0.5+I*theta)/n*log(x))) for theta in zz])
        elif moebius(n)==0:
            n_dict[n] = 0
        elif moebius(n)==-1:
            n_dict[n] = sum([(1/n)*2*real(Ei((0.5+I*theta)/n*log(x))) for theta in zz])
    return sum(n_dict.values())

def trivial_zeros_correction(x, numbersums, numberzeros, numbertrivials):
    n_dict = {n:[] for n in xrange(1, numbersums + 1)}
    for n in xrange(1, numbersums+1):
        if moebius(n)==1:
            n_dict[n] = -sum([(m/n)*Ei((-m/n)*log(x)) for m in xrange(1, numbertrivials+1)])
        elif moebius(n)==0:
            n_dict[n] = 0
        elif moebius(n)==-1:
            n_dict[n] = sum([(m/n)*Ei((-m/n)*log(x)) for m in xrange(1, numbertrivials+1)])
    return sum(n_dict.values())

def pi_G(x, numbersums, numberzeros,  numbertrivials):
    return sum([Riemann_main(x, numbersums, numberzeros, numbertrivials), nontrivial_correction(x, numbersums, numberzeros,  numbertrivials), trivial_zeros_correction(x, numbersums, numberzeros,  numbertrivials)])
︡ab8fc1bc-a922-47c3-890f-779255ce0b5c︡{"done":true}︡
︠2a28bd5e-5d64-431e-b6e7-0fd218365765s︠
%time
mindistance = 2
maxdistance = 100
nontrivialzeros = 0
trivialzeros = 100

pi_G_plot = plot(pi_G(x, ceil(log(maxdistance)-log(2))+1, nontrivialzeros, trivialzeros), (x, mindistance, maxdistance), rgbcolor = 'blue', legend_label=r'$\pi_0^G(x)$, 0 non-trivial zeros')
G_primes = list_plot(prime_counter(maxdistance), linestyle='steps-post', plotjoined=True, rgbcolor = 'red', legend_label = r'$\pi_G(x)$')

p = plot(G_primes + pi_G_plot)
p.set_legend_options(loc=(0.05, 0.8))
p.show(svg=False, legend_font_size=18, fontsize = 12, ymin=-1)
︡97ae297b-fdf7-444d-bd67-24dd250615b9︡{"file":{"filename":"/projects/67bef305-e752-49d7-a0b5-7dca07ae7d94/.sage/temp/compute5-us/6198/tmp_mZzW5v.png","show":true,"text":null,"uuid":"e71a8862-4cf8-4b57-bae3-612200f47b70"},"once":false}︡{"stdout":"\nCPU time: 6.70 s, Wall time: 7.16 s"}︡{"stdout":"\n"}︡{"done":true}︡
︠768748c2-15e1-40fc-ba31-7c1fedb2a4a4s︠
%time
mindistance = 2
maxdistance = 100
nontrivialzeros = 1
trivialzeros = 100

pi_G_plot = plot(pi_G(x, ceil(log(maxdistance)-log(2))+1, nontrivialzeros, trivialzeros), (x, mindistance, maxdistance), rgbcolor = 'blue', legend_label=r'$\pi_0^G(x)$, 1 non-trivial zero')
G_primes = list_plot(prime_counter(maxdistance), linestyle='steps-post', plotjoined=True, rgbcolor = 'red', legend_label = r'$\pi_G(x)$')

p = plot(G_primes + pi_G_plot)
p.set_legend_options(loc=(0.05, 0.8))
p.show(svg=False, legend_font_size=18, fontsize = 12)
︡545b9415-a90e-49ef-af96-319e9c8d9bf1︡{"file":{"filename":"/projects/67bef305-e752-49d7-a0b5-7dca07ae7d94/.sage/temp/compute5-us/6198/tmp_g4LuIV.png","show":true,"text":null,"uuid":"c8a1a16a-3edb-42fe-a1cd-ba6db2e81435"},"once":false}︡{"stdout":"\nCPU time: 8.50 s, Wall time: 8.84 s"}︡{"stdout":"\n"}︡{"done":true}︡
︠7598bd7f-a630-45d2-84d1-95bff0b345b4s︠
%time
mindistance = 2
maxdistance = 100
nontrivialzeros = 2
trivialzeros = 100

pi_G_plot = plot(pi_G(x, ceil(log(maxdistance)-log(2))+1, nontrivialzeros, trivialzeros), (x, mindistance, maxdistance), rgbcolor = 'blue', legend_label=r'$\pi_0^G(x)$, 2 non-trivial zeros')
G_primes = list_plot(prime_counter(maxdistance), linestyle='steps-post', plotjoined=True, rgbcolor = 'red', legend_label = r'$\pi_G(x)$')

p = plot(G_primes + pi_G_plot)
p.set_legend_options(loc=(0.05, 0.8))
p.show(svg=False, legend_font_size=18, fontsize = 12)
︡ea6650cd-046a-4a87-9ffa-46f6dbd4b73b︡{"file":{"filename":"/projects/67bef305-e752-49d7-a0b5-7dca07ae7d94/.sage/temp/compute5-us/6198/tmp_9AMHxt.png","show":true,"text":null,"uuid":"9c7ec5f9-c096-4259-a406-d3b63147c975"},"once":false}︡{"stdout":"\nCPU time: 9.38 s, Wall time: 9.51 s"}︡{"stdout":"\n"}︡{"done":true}︡
︠52d636b9-6852-4248-9e32-219a1cb2690fs︠
%time
mindistance = 2
maxdistance = 100
nontrivialzeros = 3
trivialzeros = 100

pi_G_plot = plot(pi_G(x, ceil(log(maxdistance)-log(2))+1, nontrivialzeros, trivialzeros), (x, mindistance, maxdistance), rgbcolor = 'blue', legend_label=r'$\pi_0^G(x)$, 3 non-trivial zeros')
G_primes = list_plot(prime_counter(maxdistance), linestyle='steps-post', plotjoined=True, rgbcolor = 'red', legend_label = r'$\pi_G(x)$')

p = plot(G_primes + pi_G_plot)
p.set_legend_options(loc=(0.05, 0.8))
p.show(svg=False, legend_font_size=18, fontsize = 12)
︡aa7e2469-208e-4885-8f7b-bfc5472715f3︡{"file":{"filename":"/projects/67bef305-e752-49d7-a0b5-7dca07ae7d94/.sage/temp/compute5-us/6198/tmp_j22OwT.png","show":true,"text":null,"uuid":"65da6184-5350-4e22-a34b-b8543aae9c07"},"once":false}︡{"stdout":"\nCPU time: 11.00 s, Wall time: 12.25 s"}︡{"stdout":"\n"}︡{"done":true}︡
︠455bb948-bb95-4d57-96b8-a63d88ea071cs︠
%time
mindistance = 2
maxdistance = 100
nontrivialzeros = 10
trivialzeros = 100

pi_G_plot = plot(pi_G(x, ceil(log(maxdistance)-log(2))+1, nontrivialzeros, trivialzeros), (x, mindistance, maxdistance), rgbcolor = 'blue', legend_label=r'$\pi_0^G(x)$, 10 non-trivial zeros')
G_primes = list_plot(prime_counter(maxdistance), linestyle='steps-post', plotjoined=True, rgbcolor = 'red', legend_label = r'$\pi_G(x)$')

p = plot(G_primes + pi_G_plot)
p.set_legend_options(loc=(0.05, 0.8))
p.show(svg=False, legend_font_size=18, fontsize = 12)
︡6ce4a5ec-3852-4f74-a694-f38d4d416eaa︡{"file":{"filename":"/projects/67bef305-e752-49d7-a0b5-7dca07ae7d94/.sage/temp/compute5-us/6198/tmp_AsWUQI.png","show":true,"text":null,"uuid":"03cc5460-a183-49fe-ba88-1e53ce84ec8d"},"once":false}︡{"stdout":"\nCPU time: 26.21 s, Wall time: 27.23 s"}︡{"stdout":"\n"}︡{"done":true}︡
︠d46a862f-4b77-4aab-b00d-68cb1d7c1f7ds︠
%time
mindistance = 2
maxdistance = 100
nontrivialzeros = 50
trivialzeros = 100

pi_G_plot = plot(pi_G(x, ceil(log(maxdistance)-log(2))+1, nontrivialzeros, trivialzeros), (x, mindistance, maxdistance), rgbcolor = 'blue', legend_label=r'$\pi_0^G(x)$, 50 non-trivial zeros')
G_primes = list_plot(prime_counter(maxdistance), linestyle='steps-post', plotjoined=True, rgbcolor = 'red', legend_label = r'$\pi_G(x)$')

p = plot(G_primes + pi_G_plot)
p.set_legend_options(loc=(0.05, 0.8))
p.show(svg=False, legend_font_size=18, fontsize = 12)
︡d6ff2639-7afc-467c-9c2b-09043b45e89d︡{"file":{"filename":"/projects/67bef305-e752-49d7-a0b5-7dca07ae7d94/.sage/temp/compute5-us/6198/tmp_2VDrRX.png","show":true,"text":null,"uuid":"304033e1-a4cd-4faf-80d5-5fb0f41495c1"},"once":false}︡{"stdout":"\nCPU time: 149.98 s, Wall time: 152.68 s"}︡{"stdout":"\n"}︡{"done":true}︡
︠14dbbf58-da9d-4047-9000-0f400d4ed780s︠
%time
mindistance = 2
maxdistance = 100
nontrivialzeros = 100
trivialzeros = 100

pi_G_plot = plot(pi_G(x, ceil(log(maxdistance)-log(2))+1, nontrivialzeros, trivialzeros), (x, mindistance, maxdistance), rgbcolor = 'blue', legend_label=r'$\pi_0^G(x)$, 100 non-trivial zeros')
G_primes = list_plot(prime_counter(maxdistance), linestyle='steps-post', plotjoined=True, rgbcolor = 'red', legend_label = r'$\pi_G(x)$')

p = plot(G_primes + pi_G_plot)
p.set_legend_options(loc=(0.05, 0.8))
p.show(svg=False, legend_font_size=18, fontsize = 12)
︡9497d03b-378e-4074-8ea0-a813f0636d6b︡{"file":{"filename":"/projects/67bef305-e752-49d7-a0b5-7dca07ae7d94/.sage/temp/compute5-us/6198/tmp_k4fYL_.png","show":true,"text":null,"uuid":"cdeb67c3-463b-44a5-a20e-7c88cd9ac4d0"},"once":false}︡{"stdout":"\nCPU time: 326.21 s, Wall time: 330.89 s"}︡{"stdout":"\n"}︡{"done":true}︡
︠4dd495e5-33cb-4120-a940-9353078e6ee4s︠
%time
mindistance = 2
maxdistance = 100
nontrivialzeros = 200
trivialzeros = 100

pi_G_plot = plot(pi_G(x, ceil(log(maxdistance)-log(2))+1, nontrivialzeros, trivialzeros), (x, mindistance, maxdistance), rgbcolor = 'blue', legend_label=r'$\pi_0^G(x)$, 200 non-trivial zeros')
G_primes = list_plot(prime_counter(maxdistance), linestyle='steps-post', plotjoined=True, rgbcolor = 'red', legend_label = r'$\pi_G(x)$')

p = plot(G_primes + pi_G_plot)
p.set_legend_options(loc=(0.05, 0.8))
p.show(svg=False, legend_font_size=18, fontsize = 12)
︡bf654fcd-4162-42ce-8894-5c4a5f285a05︡{"file":{"filename":"/projects/67bef305-e752-49d7-a0b5-7dca07ae7d94/.sage/temp/compute5-us/6198/tmp_zMBFwa.png","show":true,"text":null,"uuid":"62122107-af22-491d-8d94-e34a9d5090a4"},"once":false}︡{"stdout":"\nCPU time: 620.07 s, Wall time: 635.05 s"}︡{"stdout":"\n"}︡{"done":true}︡
︠7a59b044-3273-4094-b0a9-426e536795d7︠









