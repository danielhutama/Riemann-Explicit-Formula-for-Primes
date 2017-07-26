︠ceb728b0-6af9-4257-9814-6253682f85b0︠
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
︡93b65312-cccc-453f-85e1-b21cfa84895d︡{"done":true}︡
︠86672902-a810-4fcd-8cee-07d4c7fbcfc4s︠
# psi(x)
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
    for rho in zz:
        correction_sum.append(x^((0.5)+I*rho)/(0.5 + I*rho) + x^((0.5)-I*rho)/(0.5 - I*rho))
    correction_term = sum(correction_sum)
    return x - correction_term - (1/2)*log(1-x^(-2)) #- log(2*pi)
︡0ac15da9-df42-4c60-b55f-6b958e5b02b1︡{"done":true}︡
︠6685fbaf-20cb-4e2a-b3c5-ebd425f28489︠
# This cell defines Chi and Lambda. 
# This cell defines functions to plot psi(x, Chi) for a given x

def chi(q, n):
    mygroup = DirichletGroup(q).gen()
    G = mygroup.values()
    return G[n%q]

def von_man_fct(n):
    emptyset = []
    if is_prime(n):
        emptyset = [log(n)]
    else:
        for p in prime_range(1, sqrt(n)+1):
            for k in xrange(1, log(n)/log(2) + 1):
                if p^k == n:
                    emptyset.append(log(p))
                else:
                    emptyset.append(0)
    return sum(emptyset)

def psi_list(maxdistance, q):
    emptyset = []
    for n in xrange(0, maxdistance+1):
        emptyset.append(chi(q,n)*von_man_fct(n))
    return emptyset

def psi_plt(maxdistance, q):
    emptyset = []
    for x in xrange(0, len(psi_list(maxdistance, q))):
        emptyset.append(x*0)
    for i, n in enumerate(psi_list(maxdistance, q)):
        if n!=0:
            emptyset[i]=RR(n)
    return np.cumsum(emptyset)
︡9221623e-dc9a-4a8a-a1b8-55421952a2a4︡{"done":true}︡
︠071cc29f-8262-44aa-b64a-ac73ebe4db21s︠
chi_mod4=DirichletGroup(4).gen() 
L_mod4=Lfunction_from_character(chi_mod4, type="int")
L_mod4_zeros_list = L_mod4.find_zeros(0,1000,.1)
shift = real(L_mod4.value(0, 1)/L_mod4.value(0))

L_mod4_shift = real(L_mod4.value(0, 1)/L_mod4.value(0))

def psi_mod4_cong1(maxdistance):
    my_psi_mod4 = psi_list(maxdistance,4)
    emptyset = []
    for k in xrange(0, len(my_psi_mod4)):
        emptyset.append(0)
    for i,a in enumerate(my_psi_mod4):
        if a>0:
            emptyset[i] = a
    return emptyset

def psi_mod4_cong3(maxdistance):
    my_psi_mod4 = psi_list(maxdistance,4)
    emptyset = []
    for k in xrange(0, len(my_psi_mod4)):
        emptyset.append(0)
    for i,a in enumerate(my_psi_mod4):
        if a<0:
            emptyset[i] = -a
    return emptyset

def list_cumul_sum(mylist):
    return np.cumsum(mylist)

x = var('x')
def psi_mod4(x, numberzeros):
    zz = np.array(L_mod4_zeros_list[0:numberzeros])
    correction_sum = []
    for gamma in zz:
        correction_sum.append(sqrt(x)*(x^(I*gamma)/(0.5+I*gamma) + (x^(-I*gamma)/(0.5-I*gamma))))
    correction_term = sum(correction_sum)
    return -correction_term + arctanh(1/x) - shift #- log(1-x^(-1))
︡a65d613b-2be6-461f-a163-132f31745a6d︡{"done":true}︡
︠aa9691c7-dfdf-47da-a236-775f70072dbfs︠
%time
# psi(x)

NumberOfZerosUsed = 0
distance = 100

psi_not_plot = plot(psi_mod4(x, NumberOfZerosUsed), (x, 3, distance), legend_label = '$\psi_0(x, \chi)$, 0 non-trivial zeros')

psi_x_chi = list_plot(psi_plt(distance, 4), plotjoined=True, linestyle = 'steps-post', rgbcolor = 'red', legend_label = '$\psi(x, \chi)$')

p = plot(psi_x_chi + psi_not_plot)
p.set_legend_options(loc=(0.05, 0.1))
p.show(xmax = 100, ymax=5, ymin= -10, svg=False, legend_font_size=18, fontsize = 12)
︡f1fc7256-0f9f-4399-9449-b2b983210795︡{"file":{"filename":"/projects/67bef305-e752-49d7-a0b5-7dca07ae7d94/.sage/temp/compute5-us/15884/tmp_er5EAW.png","show":true,"text":null,"uuid":"866d1809-8bed-451d-9b47-ac5178edd160"},"once":false}︡{"stdout":"\nCPU time: 1.91 s, Wall time: 1.92 s"}︡{"stdout":"\n"}︡{"done":true}︡
︠3bcd901c-62b9-4a5b-b540-f2b46bbaee99︠
%time
# psi(x)

NumberOfZerosUsed = 1
distance = 100

psi_not_plot = plot(psi_mod4(x, NumberOfZerosUsed), (x, 3, distance), legend_label = '$\psi_0(x, \chi)$, 1 non-trivial zero')

psi_x_chi = list_plot(psi_plt(distance, 4), plotjoined=True, linestyle = 'steps-post', rgbcolor = 'red', legend_label = '$\psi(x, \chi)$')

p = plot(psi_x_chi + psi_not_plot)
p.set_legend_options(loc=(0.05, 0.1))
p.show(xmax = 100, ymax=5, ymin= -10, svg=False, legend_font_size=18, fontsize = 12)
︡88a1f5de-dabf-4c60-b4b2-38aecad40adf︡{"file":{"filename":"/projects/67bef305-e752-49d7-a0b5-7dca07ae7d94/.sage/temp/compute5-us/15884/tmp_byqyuD.png","show":true,"text":null,"uuid":"a7313803-3889-4209-b887-d9b0059ca231"},"once":false}︡{"stdout":"\nCPU time: 2.77 s, Wall time: 2.80 s"}︡{"stdout":"\n"}︡{"done":true}︡
︠5dfc1607-783b-46a1-856a-e88bb25083dds︠
%time
# psi(x)

NumberOfZerosUsed = 2
distance = 100

psi_not_plot = plot(psi_mod4(x, NumberOfZerosUsed), (x, 3, distance), legend_label = '$\psi_0(x, \chi)$, 2 non-trivial zeros')

psi_x_chi = list_plot(psi_plt(distance, 4), plotjoined=True, linestyle = 'steps-post', rgbcolor = 'red', legend_label = '$\psi(x, \chi)$')

p = plot(psi_x_chi + psi_not_plot)
p.set_legend_options(loc=(0.05, 0.1))
p.show(xmax = 100, ymax=5, ymin= -10, svg=False, legend_font_size=18, fontsize = 12)
︡30aba18b-c886-4b09-b2e1-8d2079aba088︡{"file":{"filename":"/projects/67bef305-e752-49d7-a0b5-7dca07ae7d94/.sage/temp/compute5-us/15884/tmp_NOJRAO.png","show":true,"text":null,"uuid":"f522ed5d-b1af-4a99-bc8e-cc9c33e0548c"},"once":false}︡{"stdout":"\nCPU time: 2.00 s, Wall time: 2.06 s"}︡{"stdout":"\n"}︡{"done":true}︡{"stdout":"\nCPU time: 2.00 s, Wall time: 2.06 s"}︡{"stdout":"\n"}︡{"done":true}︡
︠15b33894-8fdb-4932-90a3-7f67fa5ee9c0︠
%time
# psi(x)

NumberOfZerosUsed = 3
distance = 100

psi_not_plot = plot(psi_mod4(x, NumberOfZerosUsed), (x, 3, distance), legend_label = '$\psi_0(x, \chi)$, 3 non-trivial zeros')

psi_x_chi = list_plot(psi_plt(distance, 4), plotjoined=True, linestyle = 'steps-post', rgbcolor = 'red', legend_label = '$\psi(x, \chi)$')

p = plot(psi_x_chi + psi_not_plot)
p.set_legend_options(loc=(0.05, 0.1))
p.show(xmax = 100, ymax=5, ymin= -10, svg=False, legend_font_size=18, fontsize = 12)
︡8aebb4c1-3f13-4a74-a648-8f6499b272c7︡{"file":{"filename":"/projects/67bef305-e752-49d7-a0b5-7dca07ae7d94/.sage/temp/compute5-us/15884/tmp_81eBP5.png","show":true,"text":null,"uuid":"bb60f455-5b9d-4e7c-82ad-720651b1039b"},"once":false}︡{"stdout":"\nCPU time: 1.96 s, Wall time: 1.97 s"}︡{"stdout":"\n"}︡{"done":true}︡{"stdout":"\nCPU time: 1.96 s, Wall time: 1.97 s"}︡{"stdout":"\n"}︡{"done":true}︡
︠557caa17-df5e-46c8-96a3-877c5eabd328s︠
%time
# psi(x)

NumberOfZerosUsed = 10
distance = 100

psi_not_plot = plot(psi_mod4(x, NumberOfZerosUsed), (x, 3, distance), legend_label = '$\psi_0(x, \chi)$, 10 non-trivial zeros')

psi_x_chi = list_plot(psi_plt(distance, 4), plotjoined=True, linestyle = 'steps-post', rgbcolor = 'red', legend_label = '$\psi(x, \chi)$')

p = plot(psi_x_chi + psi_not_plot)
p.set_legend_options(loc=(0.05, 0.1))
p.show(xmax = 100, ymax=5, ymin= -10, svg=False, legend_font_size=18, fontsize = 12)
︡82e8307c-53d1-4109-bfc5-bf4c5dca26da︡{"file":{"filename":"/projects/67bef305-e752-49d7-a0b5-7dca07ae7d94/.sage/temp/compute5-us/15884/tmp_q4yjj5.png","show":true,"text":null,"uuid":"8e86d7ae-e2a9-4657-ae20-04c068a7ef5b"},"once":false}︡{"stdout":"\nCPU time: 2.25 s, Wall time: 2.36 s"}︡{"stdout":"\n"}︡{"done":true}︡
︠57f2d370-bd38-4b6b-b5bc-3ac3851507b6s︠
%time
# psi(x)

NumberOfZerosUsed = 50
distance = 100

psi_not_plot = plot(psi_mod4(x, NumberOfZerosUsed), (x, 3, distance), legend_label = '$\psi_0(x, \chi)$, 50 non-trivial zeros')

psi_x_chi = list_plot(psi_plt(distance, 4), plotjoined=True, linestyle = 'steps-post', rgbcolor = 'red', legend_label = '$\psi(x, \chi)$')

p = plot(psi_x_chi + psi_not_plot)
p.set_legend_options(loc=(0.05, 0.1))
p.show(xmax = 100, ymax=5, ymin= -10, svg=False, legend_font_size=18, fontsize = 12)
︡fd45b2d7-79e3-4e2c-9484-8cf217b62a93︡{"file":{"filename":"/projects/67bef305-e752-49d7-a0b5-7dca07ae7d94/.sage/temp/compute5-us/15884/tmp_0f2wVr.png","show":true,"text":null,"uuid":"8118b7ec-f4c8-416f-94fa-a2b0093823e2"},"once":false}︡{"stdout":"\nCPU time: 2.82 s, Wall time: 2.83 s"}︡{"stdout":"\n"}︡{"done":true}︡
︠6e56c8ef-13d2-468c-a29f-e5805364b607︠
%time
# psi(x)

NumberOfZerosUsed = 100
distance = 100

psi_not_plot = plot(psi_mod4(x, NumberOfZerosUsed), (x, 3, distance), legend_label = '$\psi_0(x, \chi)$, 100 non-trivial zeros')

psi_x_chi = list_plot(psi_plt(distance, 4), plotjoined=True, linestyle = 'steps-post', rgbcolor = 'red', legend_label = '$\psi(x, \chi)$')

p = plot(psi_x_chi + psi_not_plot)
p.set_legend_options(loc=(0.05, 0.1))
p.show(xmax = 100, ymax=5, ymin= -10, svg=False, legend_font_size=18, fontsize = 12)
︡4b224ce6-a9d4-48b2-b367-bac90cb32d5e︡{"file":{"filename":"/projects/67bef305-e752-49d7-a0b5-7dca07ae7d94/.sage/temp/compute5-us/15884/tmp_w2eleS.png","show":true,"text":null,"uuid":"72b0d12f-4811-430f-9665-e14e53d8fa42"},"once":false}︡{"stdout":"\nCPU time: 3.96 s, Wall time: 4.09 s"}︡{"stdout":"\n"}︡{"done":true}︡
︠9e8f8101-6cc9-41b6-bc93-7a7e37046992s︠
%time
# psi(x)

NumberOfZerosUsed = 500
distance = 100

psi_not_plot = plot(psi_mod4(x, NumberOfZerosUsed), (x, 3, distance), legend_label = '$\psi_0(x, \chi)$, 500 non-trivial zeros')

psi_x_chi = list_plot(psi_plt(distance, 4), plotjoined=True, linestyle = 'steps-post', rgbcolor = 'red', legend_label = '$\psi(x, \chi)$')

p = plot(psi_x_chi + psi_not_plot)
p.set_legend_options(loc=(0.05, 0.1))
p.show(xmax = 100, ymax=5, ymin= -10, svg=False, legend_font_size=18, fontsize = 12)
︡ab66df28-b9b9-4b0e-b3c8-e04619d69f45︡{"file":{"filename":"/projects/67bef305-e752-49d7-a0b5-7dca07ae7d94/.sage/temp/compute5-us/15884/tmp_DrbtrT.png","show":true,"text":null,"uuid":"a0f4f153-7633-4bd9-96f3-52a1954edc1b"},"once":false}︡{"stdout":"\nCPU time: 14.44 s, Wall time: 14.97 s"}︡{"stdout":"\n"}︡{"done":true}︡
︠fcdb0906-bd21-49de-ae42-277fbe270d8cs︠

︡abf2e58e-491b-451d-b4ef-0ed497a601f7︡{"done":true}︡
︠314a642b-dce0-4198-9b41-db23b705e2bf︠









