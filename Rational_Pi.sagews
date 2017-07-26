︠3c027c34-1b79-439a-9e86-2c3f7cd7b5bb︠
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
︡fc3a0a1f-8cfd-4c34-8546-18c45a600f6b︡{"done":true}︡
︠c21cbb44-2aa8-401e-9a43-df35f158f959︠
# These formulas are the components of the explicit formula.
x = var('x')
def Riemann_main(x, numbersums, numberzeros, numbertrivials):
    return sum([moebius(n)/n*Ei((1/n)*log(x)) for n in xrange(1, numbersums+1)])

def zeta_nontrivial_correction(x, numbersums, numberzeros, numbertrivials):
    zz = zeta_zeros()[0:numberzeros]
    n_dict = {n:[] for n in xrange(1, numbersums+1)}
    for n in xrange(1, numbersums+1):
        if moebius(n)==1:
            n_dict[n] = -sum([(1/n)*2*real(Ei((0.5+I*gamma)/n*log(x))) for gamma in zz])
        elif moebius(n)==0:
            n_dict[n] = 0
        elif moebius(n)==-1:
            n_dict[n] = sum([(1/n)*2*real(Ei((0.5+I*gamma)/n*log(x))) for gamma in zz])
    return sum(n_dict.values())

def trivial_zeros_correction(x, numbersums, numberzeros, numbertrivials):
    n_dict = {n:[] for n in xrange(1, numbersums + 1)}
    for n in xrange(1, numbersums+1):
        if moebius(n)==1:
            n_dict[n] = -sum([(1/n)*Ei((-2*m/n)*log(x)) for m in xrange(1, numbertrivials+1)])
        elif moebius(n)==0:
            n_dict[n] = 0
        elif moebius(n)==-1:
            n_dict[n] = sum([(1/n)*Ei((-2*m/n)*log(x)) for m in xrange(1, numbertrivials+1)])
    return sum(n_dict.values())

def pi_0(x, numbersums, numberzeros, numbertrivials):
    return sum([Riemann_main(x, numbersums, numberzeros, numbertrivials), zeta_nontrivial_correction(x, numbersums, numberzeros, numbertrivials), trivial_zeros_correction(x, numbersums, numberzeros, numbertrivials)])
︡edd09c51-e61e-4ce8-aace-4699ad2f7b32︡{"done":true}︡
︠34858163-9c84-4032-b30f-49e49faaadcbs︠
%time
mindistance = 2
maxdistance = 100
nontrivialzeros = 0
trivialzeros = 100

pi_0_plot = plot(pi_0(x, ceil(log(maxdistance)-log(2))+1, nontrivialzeros, trivialzeros), (x, mindistance, maxdistance), rgbcolor = 'blue', legend_label=r'$\pi_0(x)$, 0 non-trivial zeros')
primes = plot(prime_pi(x), (x, 0, maxdistance), rgbcolor = 'red', legend_label = r'$\pi(x)$')

p = plot(primes + pi_0_plot)
p.set_legend_options(loc=(0.05, 0.8))
p.show(svg=False, legend_font_size=18, fontsize = 12)
︡ebf6939e-3c47-4da7-be49-952c23e052e7︡{"file":{"filename":"/projects/67bef305-e752-49d7-a0b5-7dca07ae7d94/.sage/temp/compute5-us/16166/tmp_4OegZp.png","show":true,"text":null,"uuid":"777486e3-e60b-40a9-9c3b-a92ca414b27d"},"once":false}︡{"stdout":"\nCPU time: 16.64 s, Wall time: 17.22 s"}︡{"stdout":"\n"}︡{"done":true}︡
︠87871174-3fa5-4405-b3c0-14690999eee0s︠
%time
mindistance = 2
maxdistance = 100
nontrivialzeros = 1
trivialzeros = 100

pi_0_plot = plot(pi_0(x, ceil(log(maxdistance)-log(2))+1, nontrivialzeros, trivialzeros), (x, mindistance, maxdistance), rgbcolor = 'blue', legend_label=r'$\pi_0(x)$, 1 non-trivial zero')
primes = plot(prime_pi(x), (x, 0, maxdistance), rgbcolor = 'red', legend_label = r'$\pi(x)$')

p = plot(primes + pi_0_plot)
p.set_legend_options(loc=(0.05, 0.8))
p.show(svg=False, legend_font_size=18, fontsize = 12)
︡de0c5f77-e1d9-483c-95dc-80f5f5b05526︡{"file":{"filename":"/projects/67bef305-e752-49d7-a0b5-7dca07ae7d94/.sage/temp/compute5-us/16166/tmp_LZRIPy.png","show":true,"text":null,"uuid":"e994927c-66d5-4d7e-8ad8-aad7c3675015"},"once":false}︡{"stdout":"\nCPU time: 22.50 s, Wall time: 22.72 s"}︡{"stdout":"\n"}︡{"done":true}︡
︠f0a4df7d-3cf3-44b7-8798-56886765c0c9s︠
%time
mindistance = 2
maxdistance = 100
nontrivialzeros = 5
trivialzeros = 100

pi_0_plot = plot(pi_0(x, ceil(log(maxdistance)-log(2))+1, nontrivialzeros, trivialzeros), (x, mindistance, maxdistance), rgbcolor = 'blue', legend_label=r'$\pi_0(x)$, 5 non-trivial zeros')
primes = plot(prime_pi(x), (x, 0, maxdistance), rgbcolor = 'red', legend_label = r'$\pi(x)$')

p = plot(primes + pi_0_plot)
p.set_legend_options(loc=(0.05, 0.8))
p.show(svg=False, legend_font_size=18, fontsize = 12)
︡9e1589ad-3021-4723-a544-87cb0a6cf6be︡{"file":{"filename":"/projects/67bef305-e752-49d7-a0b5-7dca07ae7d94/.sage/temp/compute5-us/16719/tmp_yCt96c.png","show":true,"text":null,"uuid":"e825df33-b196-4fc9-ada3-e144ab31f874"},"once":false}︡{"stdout":"\nCPU time: 37.92 s, Wall time: 39.49 s"}︡{"stdout":"\n"}︡{"done":true}︡
︠11d60643-28d2-46bc-8c71-837ed3715c91s︠
%time
mindistance = 2
maxdistance = 100
nontrivialzeros = 10
trivialzeros = 100

pi_0_plot = plot(pi_0(x, ceil(log(maxdistance)-log(2))+1, nontrivialzeros, trivialzeros), (x, mindistance, maxdistance), rgbcolor = 'blue', legend_label=r'$\pi_0(x)$, 10 non-trivial zeros')
primes = plot(prime_pi(x), (x, 0, maxdistance), rgbcolor = 'red', legend_label = r'$\pi(x)$')

p = plot(primes + pi_0_plot)
p.set_legend_options(loc=(0.05, 0.8))
p.show(svg=False, legend_font_size=18, fontsize = 12)
︡fa65f321-e205-43c6-baf9-3818d1a91251︡{"file":{"filename":"/projects/67bef305-e752-49d7-a0b5-7dca07ae7d94/.sage/temp/compute5-us/16166/tmp_vS2Tg4.png","show":true,"text":null,"uuid":"d41b0a33-8412-462b-8fbc-f6454691d74d"},"once":false}︡{"stdout":"\nCPU time: 57.02 s, Wall time: 58.33 s"}︡{"stdout":"\n"}︡{"done":true}︡
︠77968b68-9624-4a91-b740-46d5e4eb761ds︠
%time
mindistance = 2
maxdistance = 100
nontrivialzeros = 2
trivialzeros = 100

pi_0_plot = plot(pi_0(x, ceil(log(maxdistance)-log(2))+1, nontrivialzeros, trivialzeros), (x, mindistance, maxdistance), rgbcolor = 'blue', legend_label=r'$\pi_0(x)$, 2 non-trivial zeros')
primes = plot(prime_pi(x), (x, 0, maxdistance), rgbcolor = 'red', legend_label = r'$\pi(x)$')

p = plot(primes + pi_0_plot)
p.set_legend_options(loc=(0.05, 0.8))
p.show(svg=False, legend_font_size=18, fontsize = 12)
︡eda2e435-b6b4-46dd-be26-e165da06413b︡{"file":{"filename":"/projects/67bef305-e752-49d7-a0b5-7dca07ae7d94/.sage/temp/compute5-us/16719/tmp_oDldo_.png","show":true,"text":null,"uuid":"aa9041bb-3786-4c4d-82df-172169d0aade"},"once":false}︡{"stdout":"\nCPU time: 27.69 s, Wall time: 28.68 s"}︡{"stdout":"\n"}︡{"done":true}︡
︠6f0ee64d-1271-4a53-af23-ed89f7a5ea9fs︠
%time
mindistance = 2
maxdistance = 100
nontrivialzeros = 3
trivialzeros = 100

pi_0_plot = plot(pi_0(x, ceil(log(maxdistance)-log(2))+1, nontrivialzeros, trivialzeros), (x, mindistance, maxdistance), rgbcolor = 'blue', legend_label=r'$\pi_0(x)$, 3 non-trivial zeros')
primes = plot(prime_pi(x), (x, 0, maxdistance), rgbcolor = 'red', legend_label = r'$\pi(x)$')

p = plot(primes + pi_0_plot)
p.set_legend_options(loc=(0.05, 0.8))
p.show(svg=False, legend_font_size=18, fontsize = 12)
︡1c03862b-d7ad-494a-9573-f7f161f970b1︡{"file":{"filename":"/projects/67bef305-e752-49d7-a0b5-7dca07ae7d94/.sage/temp/compute5-us/16719/tmp_wLmBBK.png","show":true,"text":null,"uuid":"cad33092-def1-4735-93c9-6c1b566658e5"},"once":false}︡{"stdout":"\nCPU time: 30.69 s, Wall time: 31.99 s"}︡{"stdout":"\n"}︡{"done":true}︡
︠fcbaa655-a9af-4237-abbd-003c246f71f8s︠
%time
mindistance = 2
maxdistance = 100
nontrivialzeros = 20
trivialzeros = 100

pi_0_plot = plot(pi_0(x, ceil(log(maxdistance)-log(2))+1, nontrivialzeros, trivialzeros), (x, mindistance, maxdistance), rgbcolor = 'blue', legend_label=r'$\pi_0(x)$, 20 non-trivial zeros')
primes = plot(prime_pi(x), (x, 0, maxdistance), rgbcolor = 'red', legend_label = r'$\pi(x)$')

p = plot(primes + pi_0_plot)
p.set_legend_options(loc=(0.05, 0.8))
p.show(svg=False, legend_font_size=18, fontsize = 12)
︡6a98cad8-95e8-4fb0-8510-c64ad8d7580b︡{"file":{"filename":"/projects/67bef305-e752-49d7-a0b5-7dca07ae7d94/.sage/temp/compute5-us/16719/tmp_eJxm_b.png","show":true,"text":null,"uuid":"3b996f42-cc4d-499a-b481-fbed844acfe5"},"once":false}︡{"stdout":"\nCPU time: 97.55 s, Wall time: 101.08 s"}︡{"stdout":"\n"}︡{"done":true}︡
︠5225a8f0-ad07-4ce0-ab32-ee879f0ae92fs︠
%time
mindistance = 2
maxdistance = 100
nontrivialzeros = 50
trivialzeros = 100

pi_0_plot = plot(pi_0(x, ceil(log(maxdistance)-log(2))+1, nontrivialzeros, trivialzeros), (x, mindistance, maxdistance), rgbcolor = 'blue', legend_label=r'$\pi_0(x)$, 50 non-trivial zeros')
primes = plot(prime_pi(x), (x, 0, maxdistance), rgbcolor = 'red', legend_label = r'$\pi(x)$')

p = plot(primes + pi_0_plot)
p.set_legend_options(loc=(0.05, 0.8))
p.show(svg=False, legend_font_size=18, fontsize = 12)
︡f5cfa1d7-dfea-462b-b5ff-b5c328d3bf53︡{"file":{"filename":"/projects/67bef305-e752-49d7-a0b5-7dca07ae7d94/.sage/temp/compute5-us/16166/tmp__qvg_y.png","show":true,"text":null,"uuid":"00afde4a-c4a1-4dc9-a281-3d04d894e7a1"},"once":false}︡{"stdout":"\nCPU time: 205.14 s, Wall time: 211.44 s"}︡{"stdout":"\n"}︡{"done":true}︡
︠febb244b-bb7e-4cb2-a53a-686ef548b719s︠
%time
mindistance = 2
maxdistance = 100
nontrivialzeros = 100
trivialzeros = 100

pi_0_plot = plot(pi_0(x, ceil(log(maxdistance)-log(2))+1, nontrivialzeros, trivialzeros), (x, mindistance, maxdistance), rgbcolor = 'blue', legend_label=r'$\pi_0(x)$, 100 non-trivial zeros')
primes = plot(prime_pi(x), (x, 0, maxdistance), rgbcolor = 'red', legend_label = r'$\pi(x)$')

p = plot(primes + pi_0_plot)
p.set_legend_options(loc=(0.05, 0.8))
p.show(svg=False, legend_font_size=18, fontsize = 12)
︡c99c8ece-d5ea-4f55-9836-09b0df1c38b7︡{"file":{"filename":"/projects/67bef305-e752-49d7-a0b5-7dca07ae7d94/.sage/temp/compute5-us/16719/tmp_biU6o1.png","show":true,"text":null,"uuid":"f46e7720-daff-4f97-babc-c89ebd943ccb"},"once":false}︡{"stdout":"\nCPU time: 329.88 s, Wall time: 341.17 s"}︡{"stdout":"\n"}︡{"done":true}︡
︠c47c5235-0c29-4cbc-bb64-310996353de1s︠
%time
mindistance = 2
maxdistance = 100
nontrivialzeros = 100
trivialzeros = 200

pi_0_plot = plot(pi_0(x, ceil(log(maxdistance)-log(2))+1, nontrivialzeros, trivialzeros), (x, mindistance, maxdistance), rgbcolor = 'blue', legend_label=r'$\pi_0(x)$, 200 non-trivial zeros')
primes = plot(prime_pi(x), (x, 0, maxdistance), rgbcolor = 'red', legend_label = r'$\pi(x)$')

p = plot(primes + pi_0_plot)
p.set_legend_options(loc=(0.05, 0.8))
p.show(svg=False, legend_font_size=18, fontsize = 12)
︡de35d21b-9499-4d41-848a-2df1add591ef︡{"file":{"filename":"/projects/67bef305-e752-49d7-a0b5-7dca07ae7d94/.sage/temp/compute5-us/16719/tmp_zQumwc.png","show":true,"text":null,"uuid":"80d7307f-9840-49e5-a853-f0e239598713"},"once":false}︡{"stdout":"\nCPU time: 393.60 s, Wall time: 407.66 s"}︡{"stdout":"\n"}︡{"done":true}︡
︠7267f6b6-176b-4abe-ac92-7348ad34c0a1︠

︠b7697370-b831-4518-8c8f-dbcf6d1f30bf︠









