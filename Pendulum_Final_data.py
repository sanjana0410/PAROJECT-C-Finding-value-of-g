# -*- coding: utf-8 -*-
"""
Created on Sat Mar  6 16:24:16 2021

@author: Sanjana
"""
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.optimize import curve_fit

## The following is the codewritten to find g using equation 24 (for small angles) and equation 22 (for all T_1 and T_2) .
#                 KNOWN VALUES OF THE PENDULUM

radius = 0.047 #Radius of discs (It's the same for both)
M = 3.67 #Total mass of whole pendulum
m1 = 1 #Mass of top disc
m2 = 1.4 #Mass od bottom disc
mr = M - m1 -m2 #Mass of rod


# Known distances used for finding COM Of pendulum
L = 1.670 #Total length of rod
p_k1 = 0.337 #Position of k1 from the top of pendulum in m
p_k2 = 0.338 #Position of k2 from the bottom of the pendulum in m
L_deltak = 0.995 # Distance between the fixed knife edges
p_m1 = 0.107 #position of m1 from the top
p_mr = 0.8375 #Position of COM of rod from the top (Half the total distance of rod (167.5/2))

#Known distances used for finding moment of inertia (about pivot point k1)
L_r0 = 0.5005 #Distance between k1 and COM of rod
L_10 = 0.23 #Distance between k1 and COM of mass1 i.e(L_k1-d1)


#                                 RAW DATA UPLOAD (FOR SMALL AMPLITUDE)
#FOR SMALL AMPLITUDE
# code below written to allocate the time period for each orientation with the respective position of m2
times_small = pd.read_excel('Reversible pendulum all data.xlsx', sheet_name="times")
xs = times_small["x1"].dropna(); #x in meters
T_1_s = times_small["t1"].dropna(); #t1 in seconds
T_2_s = times_small["t2"].dropna(); #t2 in seconds


#                                     PLOT OF RAW DATA AND INTERSECT POINTS (FOR SMALL AMPLITUDES)
#                                PLOTTING DATA POINTS ANDFINIDNIG THE INTERCEPT VALUE OF PLOT
#for small amplitudes
# Create a "figure" and set its aspect ratio. Play with the numbers to make it square, or long or short. 
fig_s = plt.figure() #figsize=(7,7)
ax = fig_s.add_subplot(111)
# This plots the T1 for varying positions x1, whilst including the uncertainty of the time period (yerr) and the position
ax.errorbar(xs,           # x coordinates
             T_1_s,              # y coordinates
             yerr = (0.2/50),     # y errors
             xerr = 0.002,
             marker='.',             # marker used is a cicle 'o'. Could be crosses 'x', or squares 's', or 'none'
             markersize = 1,        
             color='g',          
             ecolor='g',         
             markerfacecolor='g', 
             linestyle='-',       
             capsize=1,              # width of the end bit of the error bars
             )

ax.errorbar(xs,           # x coordinates
             T_2_s,              # y coordinates
             yerr = (0.2/50),     # y errors
             xerr = 0.002,
             marker='.',             # marker used is a cicle 'o'. Could be crosses 'x', or squares 's', or 'none'
             markersize = 1,        
             color='black',          
             ecolor='black',         
             markerfacecolor='black',
             linestyle='-',       
             capsize=1,              # width of the end bit of the error bars, 
             )

#Intersect of data for small amplitude
idx1 = np.argwhere(np.diff(np.sign(T_1_s - T_2_s))).flatten() #Finds the diff between the T_2 and T_2 values and appends it to idx if they're (approx) equal
plt.plot(xs[idx1], T_1_s[idx1], 'ro') #Plots the intersect points on graph
plt.plot(xs[idx1], T_2_s[idx1], 'ro')

plt.xlabel('Position of movable mass (m) ', fontsize = 12 )
plt.ylabel('Time period (s)', fontsize = 12)
plt.show()

T_period_s=[]   #List to add in intersect points
a = 0
while a<3:
    #print("Time period is 1", T_1[idx[a]])
    #print("Time period is 2", T_2[idx[a]])
    T_period_s.append(T_1_s[idx1[a]])
    T_period_s.append(T_2_s[idx1[a]])
    a += 1

#Value of T.. this is the average of T_1 and T_2 at the intersect
avg_s= sum(T_period_s)/len(T_period_s)


#                                 FINDING THE VALUE OF g WHERE T=T_1=T_2

#Value of T
print("The time period (average) when T_1=T_2 at the Intersection points for small amplitude  is", avg_s, "+/- 0.004 s")

# When T_1=T_2, L_eff is the distance between the knife edge
L_eff = L - (p_k1 + p_k2)  #Essentially it's equal to L_deltak

def g(T, l):
    return ((4*l*(np.pi**2))/(T**2))


#                                 FINDING THE UNCERTAINTY OF THE VALUE OF g AND T WHERE T=T_1=T_2
#The uncertainty of g is the square of the partial derivative of g with repsect to T summed with the partial derivative of g with respect to L_eff all rooted

u_T_avg = 0.2/50
u_L_eff = 0.002
dg_dL_eff_s = ((4*(np.pi**2))/(avg_s**2))*u_L_eff
dg_dT_s= ((-8*(np.pi**2)*L_eff)/(avg_s**3))*u_T_avg
u_g_s = ((dg_dL_eff_s**2)+(dg_dT_s**2)**0.5)

print("The value of g found using the intersect values of T_1 nad T_2  for small amplitude is:", g(avg_s, L_eff), "+/-", u_g_s)


#.........................................................................................................................
#.........................................................................................................................


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
    x1 = times_small["x1"].dropna(); #x in meters
    L_20 = x1[i] #Distance between the k1 and the COM of Mass2 
    p_m2 = 0.337 + L_20 #position of mass 2 from top. This is the distance L_20 plus the position of k1 from the top
    
    # Find the Center of Mass of pendulum
    x1 = m1*p_m1 #mass*position of COM of mass
    x2 = m2*p_m2
    xr = mr*p_mr
    X = (x1+x2+xr) / (m1+m2+mr)  
    #print ("center of mass of the pendulum", d20, "cm is {}".format(X))
    
    # Find s (distance between k1 and COM of pendulum)
    s =  X - p_k1
    #print ("distance between k1 and COM of pendulum is {}".format(s))
    
    h1 = s # Distance between k1 and COM of pendulum in meter
    h2 = 0.995 - s #Distance between k2 and the COM of pendulum in meter
    
    c = (((T_1_s[i]**2+T_2_s[i]**2)/(h1+h2))+((T_1_s[i]**2-T_2_s[i]**2)/(h1-h2)))
    results_c.append(c)   
    i += 1

# Curve fits the values of the RHS of Kater's eq with the LHS 
popt_g, pcov_g = curve_fit(x1_line, results_c, (8*np.pi**2))
a1 = popt_g[0]
print("The value of g found using Kater's equation for all T_1 nad T_2 is: =", a1, "+/-", pcov_g[0]**0.5)

