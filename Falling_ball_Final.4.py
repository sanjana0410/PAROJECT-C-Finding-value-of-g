#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 23 12:04:17 2021

@author: sanjanavinta
"""

#import libaries
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import pandas as pd

#Constant values used in defining F_C  (centrifugal force) and F_D (drag force)
θ = (51.3811*2*np.pi)/(360) #latitude of Bath university in radians
re = np.cos(θ) * 6365316 #r for Bath. This is equivalent to 'r' in the lab report
Ve = re/(24*60*60) #veloity of Earth spinning in m/s
η = 1.8347*10**(-5) #viscosity of air 
p = 1.225 #density of air
CD = 0.5 #drag coefficeint
Centri = ((Ve**2) / re)#Centrifugal constant



#            MASSES AND RADII OF BALLS FOR EACH DATASET
#None of the masses used as the same in any dataset
# For dataset_3 (3 diff balls were dropped)
m1_3 = 0.0635; m2_3 = 0.0434; m3_3 =0.0326; # in kg
r1_3 = 0.02496/2; r2_3 = 0.02197/2; r3_3 = 0.01997; # in meters
# For dataset_2 (3 diff balls were dropped)
m1_2 = 0.00703; m2_2 = 0.00406; m3_2 = 0.00209#mass of the ball from dataset 2 in kg
r1_2 = 0.01199/2; r2_2 = 0.00998/2; r3_2 = 0.00798/2 #radius of the ball from dataset 2 in meters
#For dataset_1 (2 diff balls were dropped)
m1_1 = 0.007031; m2_1 = 0.002986#mass of the ball from dataset1 in kg
r1_1 = 0.01198/2; r2_1 = 0.00898/2 #radius of the ball from dataset 1 in meters


#                HEIGHT AT WHICH BALLS WERE DROPPED FROM
#Heights for dataset_2
h1_2 = 0.451; h2_2 = 0.718; h3_2 = 0.970; h4_2 = 1.096; h5_2 = 1.209; h6_2 = 1.473; h7_2 = 1.763 # height at which balls were dropped in meters
#Heights for dataset_1
h3_1 = 1.1175; h2_1 = 0.798; h1_1 = 0.544 # height in meters
# Heights for dataset_3
h1_3 = 0.4; h2_3 = 0.6; h3_3 = 0.8; h4_3 = 1.0; h5_3 = 1.2; h6_3 = 1.4; h7_3 = 1.6; h8_3 = 1.8; h9_3 = 2.0;


# k (constant coefficient in rag formula) for different spheres
# K for each mass in dataset 3
k1_3 = 0.5*p*CD*2*np.pi*r1_3**2; k2_3 = 0.5*p*CD*2*np.pi*r2_3**2; k3_3 = 0.5*p*CD*2*np.pi*r3_3**2
# K for dataset 1
k1_1 = 0.5*p*CD*2*np.pi*r1_1**2; k2_1 = 0.5*p*CD*2*np.pi*r2_1**2;
# K for dataset 2
k1_2 = 0.5*p*CD*2*np.pi*r1_2**2; k2_2 = 0.5*p*CD*2*np.pi*r2_2**2; k3_2 = 0.5*p*CD*2*np.pi*r3_2**2;


#             TIME TAKEN FOR EACH BALL TO FALL AT EACH HEIGHT FOR ALL DATASETS
#The code below reads the excel file and allocates the data within said column into that list.
#Therefore, each list contains the time (all trials) for each mass and height. These will be the x-axis values
#Time for dataset 2
times2 = pd.read_excel('Free_fall_data_MASTER.xlsx', sheet_name="data_2")
m1_h1_2 = times2["m1_h1"].dropna(); m2_h1_2 = times2["m2_h1"].dropna(); 
m1_h2_2 = times2["m1_h2"].dropna(); m2_h2_2 = times2["m2_h2"].dropna(); m3_h2_2 = times2["m3_h2"].dropna();
m1_h3_2 = times2["m1_h3"].dropna(); m2_h3_2 = times2["m2_h3"].dropna(); m3_h3_2 = times2["m3_h3"].dropna();
m1_h4_2 = times2["m1_h4"].dropna(); m2_h4_2 = times2["m2_h4"].dropna();
m1_h5_2 = times2["m1_h5"].dropna(); m2_h5_2 = times2["m2_h5"].dropna(); m3_h5_2 = times2["m3_h5"].dropna();
m1_h6_2 = times2["m1_h6"].dropna(); m2_h6_2 = times2["m2_h6"].dropna(); m3_h6_2 = times2["m3_h6"].dropna();
m1_h7_2 = times2["m1_h7"].dropna(); m2_h7_2 = times2["m2_h7"].dropna(); 
m1_yerr_2 = times2["m1_yerr"].dropna();
m2_yerr_2 = times2["m2_yerr"].dropna();
m3_yerr_2 = times2["m3_yerr"].dropna();

#Time from dataset 1
times1 = pd.read_excel('Free_fall_data_MASTER.xlsx', sheet_name="data_1")
m1_h3_1 = times1["m1_h1"].dropna(); m2_h3_1 = times1["m2_h1"].dropna(); 
m1_h2_1 = times1["m1_h2"].dropna();
m1_h1_1 = times1["m1_h3"].dropna();
m1_yerr_1 = times1["m1_yerr"].dropna();
m2_yerr_1 = times1["m2_yerr"].dropna();

#Time for dataset 3
times3 = pd.read_excel('Free_fall_data_MASTER.xlsx', sheet_name="New")
m1_h1_3 = times3["LB_40"].dropna(); m1_h2_3 = times3["LB_60"].dropna(); m1_h3_3 = times3["LB_80"].dropna(); m1_h4_3 = times3["LB_100"].dropna(); m1_h5_3 = times3["LB_120"].dropna(); m1_h6_3 = times3["LB_140"].dropna(); m1_h7_3 = times3["LB_160"].dropna(); m1_h8_3 = times3["LB_180"].dropna();m1_h9_3 = times3["LB_200"].dropna();
m2_h1_3 = times3["MB_40"].dropna(); m2_h2_3 = times3["MB_60"].dropna(); m2_h3_3 = times3["MB_80"].dropna(); m2_h4_3 = times3["MB_100"].dropna(); m2_h5_3 = times3["MB_120"].dropna(); m2_h6_3 = times3["MB_140"].dropna(); m2_h7_3 = times3["MB_160"].dropna(); m2_h8_3 = times3["MB_180"].dropna();m2_h9_3 = times3["MB_200"].dropna();
m3_h1_3 = times3["SB_40"].dropna(); m3_h2_3 = times3["SB_60"].dropna(); m3_h3_3 = times3["SB_80"].dropna(); m3_h4_3 = times3["SB_100"].dropna(); m3_h5_3 = times3["SB_120"].dropna(); m3_h6_3 = times3["SB_140"].dropna(); m3_h7_3 = times3["SB_160"].dropna(); m3_h8_3 = times3["SB_180"].dropna();m3_h9_3 = times3["SB_200"].dropna();


# This is our expression or formula for x(t) after consiering drag and centrifugal force
#By substituting in the raw data, the curve will be fit to it and the parameter g for this fit curve will be given
def get_g(m, r,  t, g):
    A = g + ((Ve**2)/re)
    B = 0.5*np.pi*(r**2)*CD*p/m
    return (np.log(np.cosh(t*((A*B)**0.5)))/B)





#                                   DATASET 2
#The code below plots the position for each mass in dataset 2 against time. 
#Becuase there are many trials conducted for each initial height, using Mutable sequence and the lens() func.,
#we are able to overplot the points for each trial. i.e. it plots the points for each time in the list for the same height which sohuld make data more accurate

#                             for mass 1 in dataset 2
x1_2 = [ 0,
         *([h1_2-r1_2] * len(m1_h1_2)),
         *([h2_2-r1_2] * len(m1_h2_2)),
         *([h3_2-r1_2] * len(m1_h3_2)),
         *([h4_2-r1_2] * len(m1_h4_2)),
         *([h5_2-r1_2] * len(m1_h5_2)),
         *([h6_2-r1_2] * len(m1_h6_2)),
         *([h7_2-r1_2] * len(m1_h7_2))
]
#Error of the height values
error_y1_2=[0.0005] * len(x1_2)

#t for m1
time_m1_2 = [0, *m1_h1_2, *m1_h2_2, *m1_h3_2, *m1_h4_2, *m1_h5_2, *m1_h6_2, *m1_h7_2]


# Our formula for x(t) for mass 1 in dataset 2
def get_g1_2(t, g):
    return get_g(m1_2, r1_2, t, g)

#Curve fit code for m1
popt1_2, pcov1_2 = curve_fit(get_g1_2, time_m1_2, x1_2, sigma = error_y1_2)
print("g for m1 in dataset 2 =", popt1_2[0], "+/-", pcov1_2[0]**0.5 )


#                                for m2 in data2
x2_2 = [ 0,
         *([h1_2-r2_2] * len(m2_h1_2)),
         *([h2_2-r2_2] * len(m2_h2_2)),
         *([h3_2-r2_2] * len(m2_h3_2)),
         *([h4_2-r2_2] * len(m2_h4_2)),
         *([h5_2-r2_2] * len(m2_h5_2)),
         *([h6_2-r2_2] * len(m2_h6_2)),
         *([h7_2-r2_2] * len(m2_h7_2))
]

error_y2_2=[0.0005] * len(x2_2)

time_m2_2 = [0, *m2_h1_2, *m2_h2_2, *m2_h3_2, *m2_h4_2, *m2_h5_2, *m2_h6_2, *m2_h7_2]

def get_g2_2(t, g):
    return get_g(m2_2, r2_2, t, g)

popt2_2, pcov2_2 = curve_fit(get_g2_2, time_m2_2, x2_2, sigma=error_y2_2)
print("g for m2 in dataset 2 =", popt2_2[0], "+/-", pcov2_2[0]**0.5)


#                               for m3 in dataset2
x3_2 = [ 0,
         *([h2_2-r3_2] * len(m3_h2_2)),
         *([h3_2-r3_2] * len(m3_h3_2)),
         *([h5_2-r3_2] * len(m3_h5_2)),
         *([h6_2-r3_2] * len(m3_h6_2)),
]

error_y3_2=[0.0005] * len(x3_2)

time_m3_2 = [0, *m3_h2_2, *m3_h3_2, *m3_h5_2, *m3_h6_2]

def get_g3_2(t, g):
    return get_g(m3_2, r3_2, t, g)

popt3_2, pcov3_2 = curve_fit(get_g3_2, time_m3_2, x3_2, sigma = error_y3_2)
print("g for m3 in dataset 2 =", popt3_2[0], "+/-", pcov3_2[0]**0.5 )






#                                  DATASET 1
#                             for m1 in dataset 1
x1_1 = [ 0,
         *([h1_1-r1_1] * len(m1_h1_1)),
         *([h2_1-r1_1] * len(m1_h2_1)),
         *([h3_1-r1_1] * len(m1_h3_1))
]

error_y1_1=[0.005] * len(x1_1)

time_m1_1 = [0, *m1_h1_1, *m1_h2_1, *m1_h3_1]


def get_g1_1(t, g):
    return get_g(m1_1, r1_1, t, g)

popt1_1, pcov1_1 = curve_fit(get_g1_1, time_m1_1, x1_1, sigma = error_y1_1)
print("g for m1 in dataset 1 =", popt1_1[0], "+/-", pcov1_1[0]**0.5)


#                              for m2 for data 1
x2_1 = [ 0,
         *([h3_1-r2_1] * len(m2_h3_1))
]

error_y2_1=[0.005] * len(x2_1)

time_m2_1 = [0, *m2_h3_1]


def get_g2_1(t, g):
    return get_g(m2_1, r2_1, t, g)

popt2_1, pcov2_1 = curve_fit(get_g2_1, time_m2_1, x2_1, sigma=error_y2_1)
print("g for m2 for data 1 =", popt2_1[0], "+/-", pcov2_1[0]**0.5)



#                                  DATASET 3
#                           for m1 in dataset 3
x1_3 = [ 0,
         *([h1_3-r1_3] * len(m1_h1_3)),
         *([h2_3-r1_3] * len(m1_h2_3)),
         *([h3_3-r1_3] * len(m1_h3_3)),
         *([h4_3-r1_3] * len(m1_h4_3)),
         *([h5_3-r1_3] * len(m1_h5_3)),
         *([h6_3-r1_3] * len(m1_h6_3)),
         *([h7_3-r1_3] * len(m1_h7_3)),
         *([h8_3-r1_3] * len(m1_h8_3)),
         *([h9_3-r1_3] * len(m1_h9_3))
]

error_y1_3=[0.004] * len(x1_3)

time_m1_3 = [0, *m1_h1_3, *m1_h2_3, *m1_h3_3, *m1_h4_3, *m1_h5_3, *m1_h6_3, *m1_h7_3, *m1_h8_3, *m1_h9_3]


def get_g1_3(t, g):
    return get_g(m1_3, r1_3, t, g)

popt1_3, pcov1_3 = curve_fit(get_g1_3, time_m1_3, x1_3, sigma=error_y1_3)
print("g for m1 in dataset 3 =", popt1_3[0], "+/-", pcov1_3[0]**0.5)


#                           for m2 in datset 3

x2_3 = [ 0,
         *([h1_3-r2_3] * len(m2_h1_3)),
         *([h2_3-r2_3] * len(m2_h2_3)),
         *([h3_3-r2_3] * len(m2_h3_3)),
         *([h4_3-r2_3] * len(m2_h4_3)),
         *([h5_3-r2_3] * len(m2_h5_3)),
         *([h6_3-r2_3] * len(m2_h6_3)),
         *([h7_3-r2_3] * len(m2_h7_3)),
         *([h8_3-r2_3] * len(m2_h8_3)),
         *([h9_3-r2_3] * len(m2_h9_3))
]

error_y2_3=[0.004] * len(x2_3)

time_m2_3 = [0, *m2_h1_3, *m2_h2_3, *m2_h3_3, *m2_h4_3, *m2_h5_3, *m2_h6_3, *m2_h7_3, *m2_h8_3, *m2_h9_3]


def get_g2_3(t, g):
    return get_g(m2_3, r2_3, t, g)

popt2_3, pcov2_3 = curve_fit(get_g2_3, time_m2_3, x2_3, sigma=error_y2_3)
print("g for m2 in dataset 3 =", popt2_3[0], "+/-", pcov2_3[0]**0.5)


#x for m3 for data set 3

x3_3 = [ 0,
         *([h1_3-r3_3] * len(m3_h1_3)),
         *([h2_3-r3_3] * len(m3_h2_3)),
         *([h3_3-r3_3] * len(m3_h3_3)),
         *([h4_3-r3_3] * len(m3_h4_3)),
         *([h5_3-r3_3] * len(m3_h5_3)),
         *([h6_3-r3_3] * len(m3_h6_3)),
         *([h7_3-r3_3] * len(m3_h7_3)),
         *([h8_3-r3_3] * len(m3_h8_3)),
         *([h9_3-r3_3] * len(m3_h9_3))
]

error_y3_3=[0.004] * len(x3_3)

time_m3_3 = [0, *m3_h1_3, *m3_h2_3, *m3_h3_3, *m3_h4_3, *m3_h5_3, *m3_h6_3, *m3_h7_3, *m3_h8_3, *m3_h9_3]


def get_g3_3(t, g):
    return get_g(m3_3, r3_3, t, g)

popt3_3, pcov3_3= curve_fit(get_g3_3, time_m3_3, x3_3, sigma=error_y3_3)
print("g for m3 in dataset 3 =", popt3_3[0], "+/-", pcov3_3[0]**0.5)



#           FINDING THE AVERAGE VALUES OF THE VALUES OF g OBTAINED
print("AVERAGE g:")
print (((popt3_3[0]+popt2_3[0]+popt1_3[0]+popt1_1[0]+popt2_1[0]+popt1_2[0]+popt2_2[0]+popt3_2[0])/8), "+/-", ((pcov1_1[0]**0.5 +pcov2_1[0]**0.5+pcov1_2[0]**0.5 +pcov2_2[0]**0.5+pcov3_2[0]**0.5+pcov1_3[0]**0.5 +pcov2_3[0]**0.5+pcov3_3[0]**0.5)/8))

print("AVERAGE g for dataset 1:")
print(((popt1_1[0]+popt2_1[0])/2), "+/-", ((pcov1_1[0]**0.5 +pcov2_1[0]**0.5)/2))

print("AVERAGE g for dataset 2:")
print(((popt1_2[0]+popt2_2[0]+popt3_2[0])/3), "+/-", ((pcov1_2[0]**0.5 +pcov2_2[0]**0.5+pcov3_2[0]**0.5)/3))

print("AVERAGE g for dataset 3:")
print(((popt3_3[0]+popt2_3[0]+popt1_3[0])/3), "+/-", ((pcov1_3[0]**0.5 +pcov2_3[0]**0.5+pcov3_3[0]**0.5)/3))

#                                  PLOTS OF ALL DATA
#Plots of datapoints
plt.plot(time_m1_1, x1_1, 'x', color='red')
plt.plot(time_m2_1, x2_1, 'x', color='orange')
plt.plot(time_m1_2, x1_2, 'x', color='gold')
plt.plot(time_m2_2, x2_2, 'x', color='blue')
plt.plot(time_m3_2, x3_2, 'x', color='green')
plt.plot(time_m1_3, x1_3, 'x', color='violet')
plt.plot(time_m2_3, x2_3, 'x', color='pink')
plt.plot(time_m3_3, x3_3, 'x', color='brown')

#Plots of fitted curves

x_line = np.linspace(0,0.65,100)
plt.plot(x_line, get_g1_1(x_line, *popt1_1), "-g")
plt.plot(x_line, get_g2_1(x_line, *popt2_1), "-g")
plt.plot(x_line, get_g1_2(x_line, *popt1_2), "-g")
plt.plot(x_line, get_g2_2(x_line, *popt2_2), "-g")
plt.plot(x_line, get_g3_2(x_line, *popt3_2), "-g")
plt.plot(x_line, get_g1_3(x_line, *popt1_3), "-g")
plt.plot(x_line, get_g2_3(x_line, *popt2_3), "-g")
plt.plot(x_line, get_g3_3(x_line, *popt3_3), "-g")


plt.xlabel("Time taken for ball to fall [s]", fontsize = 12)
plt.ylabel("Drop height of the ball [m]", fontsize = 12)

plt.legend(("m1 for data 1","m2 for data 1","m1 for data 2", "m2 for data 2", "m3 for data 2", "m1 for data 3" , "m2 for data 3", "m3 for data 3"))
plt.show()

"""
#                           PLOT STUDYING HOW g VARIES WITH m
m1 = [m1_1, m2_1]
m2 = [m1_2, m2_2, m3_2]
m3 = [m1_3, m2_3, m3_3]
g1 = [popt1_1[0], popt2_1[0]]
g2 = [popt1_2[0], popt2_2[0], popt3_2[0]]
g3 = [popt1_3[0], popt1_3[0], popt3_3[0]]
u_g1 = [0.041, 0.091]
u_g2 = [0.023, 0.044, 0.082]
u_g3 = [0.050, 0.029, 0.099]

fig, ax = plt.subplots()
ax.errorbar(m1, g1, yerr=u_g1, marker = '.', linestyle='', color='black')
ax.errorbar(m2, g2, yerr=u_g2, marker = '.', linestyle='', color='black')
ax.errorbar(m3, g3, yerr=u_g3, marker = '.', linestyle='', color='black')
ax.set_ylabel('gravitational acceleration [ms^-2]', fontsize = 12)
ax.set_xlabel('mass [kg]', fontsize = 12)
plt.show()
"""

