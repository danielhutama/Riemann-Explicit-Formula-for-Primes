︠96934fdb-a680-4508-bb8a-d20ece3bc87c︠
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
︡7a74b981-b280-4a83-b630-aa2034dc1131︡{"done":true}︡
︠472e38ed-f808-49ff-b581-49e6cdb77e6a︠
# This cell defines psi(x) (step function), and psi_0(x) (explicit formula).
def prime_powers_list(maxdistance):
    emptyset=[]
    emptyset.append(1)
    for p in prime_range(0, maxdistance+1):
        for n in range(1, maxdistance+1):
            if p^n<=maxdistance:
                emptyset.append(p^n)
    return sorted(emptyset)

def prime_powers_list_with_zeros(maxdistance):
    emptyset = []
    for x in xrange(1, maxdistance+1):
        if x in prime_powers_list(maxdistance):
            emptyset.append(x)
        else:
            emptyset.append(0)
    return emptyset

def psi_step_maker(maxdistance):
    emptyset = []
    for x in range(0, len(prime_powers_list_with_zeros(maxdistance))):
        emptyset.append(x*0)
    for i,n in enumerate(prime_powers_list_with_zeros(maxdistance)):
        if n==1:
            emptyset[i]=RR(log(2*pi))
        for p in prime_range(0, maxdistance+1):
            for m in xrange(1, maxdistance+1):
                if n==p^m:
                    emptyset[i]=RR(log(p))
    emptyset.insert(0,0)
    return np.cumsum(emptyset)

def pi_step_maker(maxdistance):
    emptyset = []
    for x in xrange(0, maxdistance+1):
        emptyset.append(x*0)
    for i,n in enumerate(range(0, maxdistance+1)):
        if n in prime_range(0, maxdistance+1):
            emptyset[i]=1
    return np.cumsum(emptyset)

x = var('x')
def psi_not(x, numberzeros):
    zz = zeta_zeros()[0:numberzeros]
    correction_sum = []
    for gamma in zz:
        correction_sum.append((cos(gamma*log(x))+2*gamma*sin(gamma*log(x)))/(1+4*gamma^2))
    correction_term = sum(correction_sum)
    return x - 4*x^(1/2)*correction_term #- (1/2)*log(1-x^(-2)) #- log(2*pi)
︡1885006a-db07-49d7-bab6-31b1ad8da26f︡{"done":true}︡
︠53f45faf-07a4-46d5-bf98-cf92066b728es︠

%time
# psi(x)
NumberOfZerosUsed = 0
distance = 20

psi_not_plot = plot(psi_not(x, NumberOfZerosUsed), (x, 1, distance), legend_label = '$\psi_0(x) + \log(2\pi)$, 0 non-trivial zeros')
psi = list_plot(psi_step_maker(distance), plotjoined=True, linestyle = 'steps-post', rgbcolor = "red", legend_label = '$\psi(x) + \log(2\pi)$')
fundamental = plot(x, (x, 0, distance), rgbcolor= 'green')


p = plot(psi + psi_not_plot)
p.set_legend_options(loc=(0.05, 0.8))
p.show(svg=False, legend_font_size=18, fontsize = 12)
︡57b15497-e666-4dbf-bf6f-78f0b45a6f48︡{"file":{"filename":"/projects/67bef305-e752-49d7-a0b5-7dca07ae7d94/.sage/temp/compute5-us/9794/tmp_KNpiAu.png","show":true,"text":null,"uuid":"9778a82b-4af3-4ff8-8c0a-d0ce0a195ca3"},"once":false}︡{"stdout":"\nCPU time: 0.58 s, Wall time: 0.58 s"}︡{"stdout":"\n"}︡{"done":true}︡
︠cff84514-3d76-472b-bafc-6df635c5097ds︠
%time
# psi(x)
NumberOfZerosUsed = 1
distance = 20

psi_not_plot = plot(psi_not(x, NumberOfZerosUsed), (x, 1, distance), legend_label = '$\psi_0(x) + \log(2\pi)$, 1 non-trivial zero')
psi = list_plot(psi_step_maker(distance), plotjoined=True, linestyle = 'steps-post', rgbcolor = "red", legend_label = '$\psi(x) + \log(2\pi)$')
fundamental = plot(x, (x, 0, distance), rgbcolor= 'green')


p = plot(psi + psi_not_plot)
p.set_legend_options(loc=(0.05, 0.8))
p.show(svg=False, legend_font_size=18, fontsize = 12)

︡265aa463-30ae-448d-9247-a3c8f5ae83ff︡{"file":{"filename":"/projects/67bef305-e752-49d7-a0b5-7dca07ae7d94/.sage/temp/compute5-us/9794/tmp_GYNrcG.png","show":true,"text":null,"uuid":"67f1d778-18f2-4efa-962c-2b86b1a1f541"},"once":false}︡{"stdout":"\nCPU time: 0.78 s, Wall time: 0.99 s"}︡{"stdout":"\n"}︡{"done":true}︡
︠14109497-d78f-4f64-9fa2-f89cd0859ee7s︠
%time
# psi(x)
NumberOfZerosUsed = 10
distance = 20

psi_not_plot = plot(psi_not(x, NumberOfZerosUsed), (x, 1, distance), legend_label = '$\psi_0(x) + \log(2\pi)$, 10 non-trivial zeros')
psi = list_plot(psi_step_maker(distance), plotjoined=True, linestyle = 'steps-post', rgbcolor = "red", legend_label = '$\psi(x) + \log(2\pi)$')
fundamental = plot(x, (x, 0, distance), rgbcolor= 'green')

p = plot(psi + psi_not_plot)
p.set_legend_options(loc=(0.05, 0.8))
p.show(svg=False, legend_font_size=18, fontsize = 12)
︡e1aa73e1-8f9c-4fb7-8efa-7da006397460︡{"file":{"filename":"/projects/67bef305-e752-49d7-a0b5-7dca07ae7d94/.sage/temp/compute5-us/9794/tmp_kBwQKb.png","show":true,"text":null,"uuid":"69af3098-97ac-4fcf-b9d4-69e5cdda1a43"},"once":false}︡{"stdout":"\nCPU time: 1.00 s, Wall time: 1.27 s"}︡{"stdout":"\n"}︡{"done":true}︡
︠6685fbaf-20cb-4e2a-b3c5-ebd425f28489s︠
%time
# psi(x)
NumberOfZerosUsed = 100
distance = 20

psi_not_plot = plot(psi_not(x, NumberOfZerosUsed), (x, 1, distance), legend_label = '$\psi_0(x) + \log(2\pi)$, 100 non-trivial zeros')
psi = list_plot(psi_step_maker(distance), plotjoined=True, linestyle = 'steps-post', rgbcolor = "red", legend_label = '$\psi(x) + \log(2\pi)$')
fundamental = plot(x, (x, 0, distance), rgbcolor= 'green')

p = plot(psi + psi_not_plot)
p.set_legend_options(loc=(0.05, 0.8))
p.show(svg=False, legend_font_size=18, fontsize = 12)
︡f92fb03d-50fb-4287-a281-b7a7365237eb︡{"file":{"filename":"/projects/67bef305-e752-49d7-a0b5-7dca07ae7d94/.sage/temp/compute5-us/9794/tmp_3IWO2K.png","show":true,"text":null,"uuid":"b86a0724-25bb-43cd-bc29-47de8ee4230c"},"once":false}︡{"stdout":"\nCPU time: 3.33 s, Wall time: 3.61 s"}︡{"stdout":"\n"}︡{"done":true}︡
︠77e5033d-cb29-43da-8a05-d35bf77dd5f2s︠
%time
# psi(x)
NumberOfZerosUsed = 100
distance = 20

psi_not_plot = plot(psi_not(x, NumberOfZerosUsed), (x, 1, distance), legend_label = '$\psi_0(x) + \log(2\pi)$, 100 non-trivials, 0 trivials')
psi = list_plot(psi_step_maker(distance), plotjoined=True, linestyle = 'steps-post', rgbcolor = "red", legend_label = '$\psi(x) + \log(2\pi)$')
fundamental = plot(x, (x, 0, distance), rgbcolor= 'green')

p = plot(psi + psi_not_plot)
p.set_legend_options(loc=(0.05, 0.8))
p.show(svg=False, legend_font_size=18, fontsize = 12)
︡320a8001-44b9-43f7-a1d7-e1d738c46c73︡{"file":{"filename":"/projects/67bef305-e752-49d7-a0b5-7dca07ae7d94/.sage/temp/compute5-us/9794/tmp_FpmwTJ.png","show":true,"text":null,"uuid":"ffdebd41-67f3-45fc-a75b-5a159e19edad"},"once":false}︡{"stdout":"\nCPU time: 3.29 s, Wall time: 3.49 s"}︡{"stdout":"\n"}︡{"done":true}︡
︠a987cc93-bf5d-4302-9892-64a9e7c17c3e︠
︠1af70f3b-bd76-4e2a-879a-42a1467ca624︠









