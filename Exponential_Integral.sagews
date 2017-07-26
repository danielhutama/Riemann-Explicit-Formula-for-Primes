︠a47dade0-6010-4770-ad76-ce4d24677ace︠
#Not all imports are used.

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
︡a88f1f57-6f76-4485-9967-d1ecc578cd9f︡{"done":true}︡
︠4d618172-e128-4c00-a35e-cc7f528e9dba︠
# compares x^(1/2+ i*gamma1)/((1/2 + i*gamma1)log(x)) vs Ei((1/2 + i*gamma1)*log(x))

%time
zz = zeta_zeros()[0:10]
gamma = zz[0]

x = var('x')

mindist = 100
maxdist = 150
Ei_arg = (1/2+I*gamma)*log(x)

Ei_rho_logx = plot(real(2*Ei(Ei_arg)), (x, mindist, maxdist), rgbcolor = 'blue', legend_label=r'$2\Re(Ei(\rho\log(x))$')


def li_approx(x):
    return 4*sqrt(x)/log(x)*(cos(gamma*log(x))+2*gamma*sin(gamma*log(x)))/(1+4*gamma^2)

approx_plot = plot(li_approx(x), (x, mindist, maxdist), rgbcolor = 'red', legend_label = r'$2\Re \left(\frac{x^\rho}{\rho \log(x)}\right)$')

plot(Ei_rho_logx + approx_plot).show(svg=False, legend_font_size=18, fontsize = 12)

︡de79f463-fa11-4a8b-bc0d-f2da47d68744︡{"file":{"filename":"/projects/67bef305-e752-49d7-a0b5-7dca07ae7d94/.sage/temp/compute5-us/8435/tmp_MbBuhd.png","show":true,"text":null,"uuid":"de150541-4b2b-4574-8ba6-9f6db1829d15"},"once":false}︡{"stdout":"\nCPU time: 2.00 s, Wall time: 14.93 s"}︡{"stdout":"\n"}︡{"done":true}︡
︠dd798c4a-6b06-4f3a-8940-f57977c5fac9s︠
%time
zz = zeta_zeros()[0:10]
gamma = zz[0]

x = var('x')

mindist = 115
maxdist = 120
Ei_arg = (1/2+I*gamma)*log(x)

Ei_rho_logx = plot(real(2*Ei(Ei_arg)), (x, mindist, maxdist), rgbcolor = 'blue', legend_label=r'$2\Re(Ei(\rho\log(x))$')


def li_approx(x):
    return 4*sqrt(x)/log(x)*(cos(gamma*log(x))+2*gamma*sin(gamma*log(x)))/(1+4*gamma^2)

approx_plot = plot(li_approx(x), (x, mindist, maxdist), rgbcolor = 'red', legend_label = r'$2\Re \left(\frac{x^\rho}{\rho \log(x)}\right)$')

plot(Ei_rho_logx + approx_plot).show(svg=False, legend_font_size=18, fontsize = 12)
︡5e85c93b-0bc4-47c7-ac52-d6ff8960002b︡{"file":{"filename":"/projects/67bef305-e752-49d7-a0b5-7dca07ae7d94/.sage/temp/compute5-us/8435/tmp_qNXLi4.png","show":true,"text":null,"uuid":"5decffa6-b785-485e-a700-3c10540d7e11"},"once":false}︡{"stdout":"\nCPU time: 1.19 s, Wall time: 1.52 s"}︡{"stdout":"\n"}︡{"done":true}︡
︠c9d8da7a-8b5f-40ad-a93d-c179831c6a85s︠
%time
zz = zeta_zeros()[0:11]
gamma = zz[1]

x = var('x')

mindist = 100
maxdist = 150
Ei_arg = (1/2+I*gamma)*log(x)

Ei_rho_logx = plot(real(2*Ei(Ei_arg)), (x, mindist, maxdist), rgbcolor = 'blue', legend_label=r'$2\Re(Ei(\rho\log(x))$')


def li_approx(x):
    return 4*sqrt(x)/log(x)*(cos(gamma*log(x))+2*gamma*sin(gamma*log(x)))/(1+4*gamma^2)

approx_plot = plot(li_approx(x), (x, mindist, maxdist), rgbcolor = 'red', legend_label = r'$2\Re \left(\frac{x^\rho}{\rho \log(x)}\right)$')

p = plot(Ei_rho_logx + approx_plot)
p.set_legend_options(loc=(0.42, 0.08))
p.show(svg=False, legend_font_size=18, fontsize = 12)
︡a5a8d4fa-3b8e-4f2c-8a0b-81414bd2fd25︡{"file":{"filename":"/projects/67bef305-e752-49d7-a0b5-7dca07ae7d94/.sage/temp/compute5-us/8435/tmp_tIEHtL.png","show":true,"text":null,"uuid":"d50222a7-dba3-46d3-bb52-ad1b90b72356"},"once":false}︡{"stdout":"\nCPU time: 0.74 s, Wall time: 0.80 s"}︡{"stdout":"\n"}︡{"done":true}︡
︠bd2147b8-c123-42e1-a1d1-ecc2af319c0bs︠
%time
zz = zeta_zeros()[0:10]
gamma = zz[1]

x = var('x')

mindist = 115
maxdist = 120
Ei_arg = (1/2+I*gamma)*log(x)

Ei_rho_logx = plot(real(2*Ei(Ei_arg)), (x, mindist, maxdist), rgbcolor = 'blue', legend_label=r'$2\Re(Ei(\rho\log(x))$')


def li_approx(x):
    return 4*sqrt(x)/log(x)*(cos(gamma*log(x))+2*gamma*sin(gamma*log(x)))/(1+4*gamma^2)

approx_plot = plot(li_approx(x), (x, mindist, maxdist), rgbcolor = 'red', legend_label = r'$2\Re \left(\frac{x^\rho}{\rho \log(x)}\right)$')

p = plot(Ei_rho_logx + approx_plot)
p.set_legend_options(loc=(0.65, 0.1))
p.show(svg=False, legend_font_size=18, fontsize = 12)
︡5ee56cde-2d44-4a25-bd4a-e4efaf504d25︡{"file":{"filename":"/projects/67bef305-e752-49d7-a0b5-7dca07ae7d94/.sage/temp/compute5-us/8435/tmp_agpZzm.png","show":true,"text":null,"uuid":"c3832032-8d70-4885-98ac-303167ce1bd8"},"once":false}︡{"stdout":"\nCPU time: 0.72 s, Wall time: 0.73 s"}︡{"stdout":"\n"}︡{"done":true}︡
︠a0df8888-35e1-416d-9cab-4cd56662a149s︠
%time
zz = zeta_zeros()[0:11]
gamma = zz[1]

x = var('x')

mindist = 10^300
maxdist = 10^301
Ei_arg = (1/2+I*gamma)*log(x)

Ei_rho_logx = plot(real(2*Ei(Ei_arg)), (x, mindist, maxdist), rgbcolor = 'blue', legend_label=r'$2\Re(Ei(\rho\log(x))$')


def li_approx(x):
    return 4*sqrt(x)/log(x)*(cos(gamma*log(x))+2*gamma*sin(gamma*log(x)))/(1+4*gamma^2)

approx_plot = plot(li_approx(x), (x, mindist, maxdist), rgbcolor = 'red', legend_label = r'$2\Re \left(\frac{x^\rho}{\rho \log(x)}\right)$')

plot(Ei_rho_logx + approx_plot).show(svg=False)
︡55674516-3fbc-492e-83f3-851f2c7aa286︡{"file":{"filename":"/projects/67bef305-e752-49d7-a0b5-7dca07ae7d94/.sage/temp/compute5-us/8435/tmp_XXs48M.png","show":true,"text":null,"uuid":"a5a4e24a-92d0-4bd8-b419-4cad0e139328"},"once":false}︡{"stdout":"\nCPU time: 0.80 s, Wall time: 0.82 s"}︡{"stdout":"\n"}︡{"done":true}︡
︠13ba5dec-4835-460d-be78-9c342b205577s︠
︡c5229894-9196-443a-a8ff-013040e35319︡{"done":true}︡
︠7d0b2017-0bc4-4068-bb84-6bc6bf4c91e6︠










