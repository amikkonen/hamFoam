#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 21 19:20:11 2019

Steady state 1D validation. 

@author: Antti Mikkonen, a.mikkonen@iki.fi
"""

import scipy as sp
from scipy import optimize
import os
from matplotlib import pyplot as plt
from scipy import interpolate


T_l = 293
T_r = 263

#rho_C1 = 574
#cp_C1 =  1100
k_C1 =  0.19
#rho_C11 =  980
#cp_C11 = 2300
k_C11 = 0.4
#rho_D2 = 37
#cp_D2 = 850
k_D2 = 0.0336
#rho_A3 = 73
#cp_A3 = 850
k_A3 = 0.0305

dx = 1e-3
L_C1 = 0.013
L_C11 = 0.001
L_D2 = 0.173
L_A3 = 0.03
xs = [L_C1,
      L_C1 + L_C11,
      L_C1 + L_C11 + L_D2 ,
      ]

L = L_C1 + L_C11 + L_D2 + L_A3




with open(os.path.join("postProcessing_harmonic","sampleLine","1", "out_T_rho_cp_k.xy"), "r") as ifile:
    lines = ifile.readlines()

x_cfd_harmonic = []
T_cfd_harmonic = []
for line in lines:
    parts = [float(part) for part in line.split()]
    x_cfd_harmonic.append(parts[0])
    T_cfd_harmonic.append(parts[1])
x_cfd_harmonic = sp.array(x_cfd_harmonic)    
T_cfd_harmonic = sp.array(T_cfd_harmonic)    

with open(os.path.join("postProcessing_linear","sampleLine","1", "out_T_rho_cp_k.xy"), "r") as ifile:
    lines = ifile.readlines()

x_cfd_linear = []
T_cfd_linear = []
for line in lines:
    parts = [float(part) for part in line.split()]
    x_cfd_linear.append(parts[0])
    T_cfd_linear.append(parts[1])
x_cfd_linear = sp.array(x_cfd_linear)    
T_cfd_linear = sp.array(T_cfd_linear)    



#
def eq(T):
    q0 = -k_C1/L_C1*(T[0]-T_l)
    q1 = -k_C11/L_C11*(T[1]-T[0])
    q2 = -k_D2/L_D2*(T[2]-T[1])
    q3 = -k_A3/L_A3*(T_r-T[2])
    
    return [q0-q1, q0-q2, q0-q3]


T = [290, 288, 270]
sol = optimize.root(eq, T, tol=1e-20)
T = sol.x

x_ana = sp.array([0]+list(xs)+[L])
T_ana = [T_l]+list(T)+[T_r]
T_func = interpolate.interp1d(x_ana, T_ana)

Ta_harmonic = T_func(x_cfd_harmonic)
er_harmonic = sp.absolute(Ta_harmonic - T_cfd_harmonic)
print("Max difference harmonic(C)", er_harmonic.max())

Ta_linear = T_func(x_cfd_linear)
er_linear = sp.absolute(Ta_linear - T_cfd_linear)
print("Max difference linear(C)", er_linear.max())



fig, axes = plt.subplots(1,2, figsize=(7,3))

ax = axes[0]
ax.plot(x_cfd_harmonic*1e3,T_cfd_harmonic, "r-", label="CFD harmonic")

ax.text(0,267, "Harmonic Max diff %f C " % er_harmonic.max())

ax.plot(x_cfd_linear*1e3,T_cfd_linear, "b--", label="CFD linear")
ax.text(0,263, "Linear Max diff %f C " % er_linear.max())
ax.plot(x_ana*1e3,T_ana, "kd", label="Analytical")


    
ax = axes[1]
ax.plot(x_cfd_harmonic*1e3,T_cfd_harmonic, "rd", label="CFD harmonic")
ax.plot(x_cfd_linear*1e3,T_cfd_linear, "bx", label="CFD linear")
ax.plot(x_ana*1e3,T_ana, "kd-", label="Analytical")



ax.set_xlim(xs[0]*1e3-2, xs[1]*1e3+2)
ax.set_ylim(292, 293)

for ax in axes.ravel():
    ax.grid(True)
    ax.set_xlabel("x (mm)")
    ax.set_ylabel("T (C)")
    ax.legend(frameon=False)



fig.tight_layout()
fig.savefig("harmonicSchemeValidation.pdf")

