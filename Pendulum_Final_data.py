# -*- coding: utf-8 -*-
"""
Created on Sat Mar  6 16:24:16 2021

@author: Lowri
"""
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.optimize import curve_fit

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

#Times
times = pd.read_excel('Reversible pendulum all data.xlsx', sheet_name="times")
x1 = times["x1"].dropna(); #x in meters
T_1 = times["t1"].dropna(); #t1 in seconds
T_2 = times["t2"].dropna(); #t2 in seconds

# Create a "figure" and set its aspect ratio. Play with the numbers to make it square, or long or short. 
fig = plt.figure() #figsize=(7,7)
# Here we only want one plot in the figure to a 1 by 1 set of plots and this is number 1. Leave alone for now. 
ax = fig.add_subplot(111)
#ax.set_xlim(0.25, 0.35)
#ax.set_ylim(1.98, 2.02)
# This nest bit does a lot of work and plots a graph with error bars.
ax.errorbar(x1,           # x coordinates
             T_1,              # y coordinates
             yerr = (0.2/50),     # y errors
             xerr = 0.002,
             marker='o',             # marker used is a cicle 'o'. Could be crosses 'x', or squares 's', or 'none'
             markersize = 7,        # marker size
             color='r',          # overall colour I think ( red)
             ecolor='r',         # edge colour for you marker (red)
             markerfacecolor='r', # red
             linestyle='-',       # no line joining markers, could be a line '-', or a dashed line '--'
             capsize=6,              # width of the end bit of the error bars, not too small nor too big please.
             )

ax.errorbar(x1,           # x coordinates
             T_2,              # y coordinates
             yerr = (0.2/50),     # y errors
             xerr = 0.002,
             marker='x',             # marker used is a cicle 'o'. Could be crosses 'x', or squares 's', or 'none'
             markersize = 7,        # marker size
             color='black',          # overall colour I think
             ecolor='black',         # edge colour for you marker
             markerfacecolor='black',
             linestyle='-',       # no line joining markers, could be a line '-', or a dashed line '--'
             capsize=6,              # width of the end bit of the error bars, not too small nor too big please.
             )
plt.xlabel('x1 of movable mass (m) ', fontsize = 20, fontname='Times New Roman' )
plt.ylabel('time (seconds)', fontsize = 20, fontname='Times New Roman')

def x1_line(c, g):
    return g*c

results_c=[]

#CHANGE
i = 0
while i < 35:
    x1 = times["x1"].dropna(); #x in meters
    d20 = x1[i] #Distance between the pivot and COM od Mass2
    d2 = 0.337 + d20 #position of mass 2 from top
    
    # Find the Center of Mass of pendulum
    x1 = m1*d1
    x2 = m2*d2
    xr = mr*dr
    X = (x1+x2+xr) / (m1+m2+mr)
    #print ("center of mass of the pendulum", d20, "cm is {}".format(X))
    
    # Find s (distance between pivot and COM of pendulum)
    s =  X - p
    #print ("distance between pivot and COM of pendulum is {}".format(s))
    
    h1 = s #meter
    h2 = 0.995 - s #meter
    
    c = (((T_1[i]**2+T_2[i]**2)/(h1+h2))+((T_1[i]**2-T_2[i]**2)/(h1-h2)))
    
    results_c.append(c)
        
    i += 1

    
popt_g, pcov_g = curve_fit(x1_line, results_c, (8*np.pi**2))
a1 = popt_g[0]
print("g1 =", a1)
print("uncertainty g1 = ", pcov_g[0]**0.5)


def period(l, g):
    2*np.pi*np.sqrt(l/g)

#x_line = np.linspace(7.9,8.1,100)
#popt_1, pcov_1 = curve_fit(period, x1, T_1)
#popt_2, pcov_2 = curve_fit(period, x1, T_2)
#plt.plot(x_line, period(x1, a1), "-k")
#plt.plot(x_line, period(x_line, a1), "-k")

