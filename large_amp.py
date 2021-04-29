#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 26 16:40:44 2021

@author: sanjanavinta
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.optimize import curve_fit

## The following code iswritten to find g using equation 24 for large amplitudes.

##KNOWN VALUES IN THE PENDULUM
D = 1.670 #Total length of rod
p = 0.337 #Position of pivot from top
radius = 0.047 #Radius of discs (It's te same for both)
M = 3.67 #Total mass of whole pendulum
m1 = 1 #Mass of top disc
m2 = 1.4 #Mass od bottom disc
mr = M - m1 -m2 #Mass of rod

# Known distances used for finding COM Of pendulum
d1 = 0.107 #position of Mass 1 from the top
dr = 0.8375 #Position of COM of rod from the top (Half the total distance of rod (167.5/2))

#Known distances used for finding moment of inertia
dr0 = 0.311 #Distance between the pivot and COM of rod
d10 = 0.23 #Distance between the pivot and COM of mass1

#Times (for large amplitude)
times = pd.read_excel('Reversible pendulum all data.xlsx', sheet_name="times")
x2 = times["x2"].dropna(); #x in meters
T_1_2 = times["t1_2"].dropna(); #t1 in seconds
T_2_2 = times["t2_2"].dropna(); #t2 in seconds

l_line = np.linspace(0.15, 0.85, 100)

# This is the deifinition of the fitted period of T_1 and T_2... Since there are closer points, it's better to curve fit them and find the intersect point
def periodT(l, a, b, c):
    return (a*(l**2)) + (b*l) + c

popt_T1, pcov_T1 = curve_fit(periodT, x2, T_1_2)    #Curve fit function.. Fits T_1 nad T_2 to a function
popt_T2, pcov_T2 = curve_fit(periodT, x2, T_2_2)


#intercepts
l1 = 0.7306449629153

print("when l = ", l1, "then T1 = ", periodT(l1, popt_T1[0], popt_T1[1], popt_T1[2])) # Finds the value of T_1 at the intersecption point. .I.e it puts the intercept value into the curve fit equation
print("when l = ", l1, "then T2 = ", periodT(l1, popt_T2[0], popt_T2[1], popt_T2[2]))

T = 2.01
angle = 16*np.pi/180

g = 0.995/((T/(2*np.pi*(1+((1/16)*angle**2)+(11/3072)*angle**4)))**2)   # Modified eqaution of g.. it takes into account large amplitudes.
print(g)