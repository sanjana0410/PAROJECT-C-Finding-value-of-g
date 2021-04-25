# -*- coding: utf-8 -*-
"""
Created on Sat Mar  6 16:24:16 2021

@author: Lowri
"""
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.optimize import curve_fit

#                 KNOWN VALUES OF THE PENDULUM
L = 1.670 #Total length of rod
p = 0.337 #Position of pivot from top in m
p_bottom = 0.338 #Position of pivot on the other end in m
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

#                                 RAW DATA UPLOAD
# code below written to allocate the time period for each orientation with the respective position of m2
times = pd.read_excel('Reversible pendulum all data.xlsx', sheet_name="times")
x1 = times["x1"].dropna(); #x in meters
T_1 = times["t1"].dropna(); #t1 in seconds
T_2 = times["t2"].dropna(); #t2 in seconds

#                                     PLOT OF RAW DATA AND INTERSECT POINTS
#                                PLOTTING DATA POINTS
# Create a "figure" and set its aspect ratio. Play with the numbers to make it square, or long or short. 
fig = plt.figure() #figsize=(7,7)
# Here we only want one plot in the figure to a 1 by 1 set of plots and this is number 1. Leave alone for now. 
ax = fig.add_subplot(111)
#ax.set_xlim(0.25, 0.35)
#ax.set_ylim(1.98, 2.02)
# This next bit does a lot of work and plots a graph with error bars.
# This plots the T1 for varying positions x1, whilst including the uncertainty of the time period (yerr) and the position
ax.errorbar(x1,           # x coordinates
             T_1,              # y coordinates
             yerr = (0.2/50),     # y errors
             xerr = 0.002,
             marker='.',             # marker used is a cicle 'o'. Could be crosses 'x', or squares 's', or 'none'
             markersize = 1,        # marker size
             color='g',          # overall colour I think ( red)
             ecolor='g',         # edge colour for you marker (red)
             markerfacecolor='g', # red
             linestyle='-',       # no line joining markers, could be a line '-', or a dashed line '--'
             capsize=1,              # width of the end bit of the error bars, not too small nor too big please.
             )

ax.errorbar(x1,           # x coordinates
             T_2,              # y coordinates
             yerr = (0.2/50),     # y errors
             xerr = 0.002,
             marker='.',             # marker used is a cicle 'o'. Could be crosses 'x', or squares 's', or 'none'
             markersize = 1,        # marker size
             color='black',          # overall colour I think
             ecolor='black',         # edge colour for you marker
             markerfacecolor='black',
             linestyle='-',       # no line joining markers, could be a line '-', or a dashed line '--'
             capsize=1,              # width of the end bit of the error bars, not too small nor too big please.
             )

#                     FINDING INTERSECT VALUES
#Finds the intersect values
idx = np.argwhere(np.diff(np.sign(T_1 - T_2))).flatten() #Finds the diff between the T_2 and T_2 values and appends it to idx if they're (approx) equal
plt.plot(x1[idx], T_1[idx], 'ro') #Plots the intersect points on graph
plt.plot(x1[idx], T_2[idx], 'ro')

plt.xlabel('Position of movable mass (m) ', fontsize = 10 )
plt.ylabel('Time period (s)', fontsize = 10)
plt.show()

T_period=[]   #List to add in intersect points
a = 0
while a<3:
    #print("Time period is 1", T_1[idx[a]])
    #print("Time period is 2", T_2[idx[a]])
    T_period.append(T_1[idx[a]])
    T_period.append(T_2[idx[a]])
    a += 1


#                                 FINDING THE VALUE OF g WHERE T=T_1=T_2
#Value of T.. this is the average of T_1 and T_2 at the intersect
avg= sum(T_period)/len(T_period)
print("The time period (average) when T_1=T_2 at the Intersection points is", avg)

# When T_1=T_2, L_eff is the distance between the knife edge
L_eff = L - (p + p_bottom)

def g(T, l):
    return ((4*l*(np.pi**2))/(T**2))
print("The value of g found using the intersect values of T_1 nad T_2 is:", g(avg, L_eff))





#                          FINDING g FOR KATER'S PENDULUM USING ALL T_1 AND T_2
#We have T1 and T2 for each x. Once we find h1 and h2 for each x, we will compute it into a function which will then be curve fit to find an idela value for g
#Definition of 8pi^2 where c represents the RHS of Kater's equation
def x1_line(c, g):
    return g*c  

#List in which the RHS value for each x is stored
results_c=[]

#                            FINDING h1 AND h2 FOR EACH x
# As the position x for mass 2 is varied, the COM of the pendulum changes. Thus the values of h1 and h2 also change. 
#Here the value of h1 and h2 are found for each position of x by finding the new COM
i = 0
while i < 35:
    x1 = times["x1"].dropna(); #x in meters
    d20 = x1[i] #Distance between the pivot and COM of Mass2 #i.e Distance between the pivot and the Mass 2 plus the radii of mass2
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

# Curve fits the values of the RHS of Kater's eq with the LHS 
popt_g, pcov_g = curve_fit(x1_line, results_c, (8*np.pi**2))
a1 = popt_g[0]
print("The value of g found using LKater's equation for all T_1 nad T_2 is: =", a1, "+/-", pcov_g[0]**0.5)


