import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from drawnow import drawnow, figure

"""
Function Definations
"""

# Model Defination
def two_dof_cubic(x, t,F0,omg,m,c,k1,k2,a3,a3_12):
    x1 = x[1]
    x1d = F0 * np.cos(omg * t)/m - c * x[1] / m - k1 * x[0] / m - a3 * x[0]**3/m - a3_12*(x[0]-x[2])**3
    x2 = x[3]
    x2d =  - c * x[3] / m - k2 * x[2] / m - a3_12*(x[2]-x[0])**3
    return np.asarray([x1, x1d, x2, x2d])

# RK 4th order function
def my_RK4(two_dof_cubic, init, t,F0,omg,m,c,k1,k2,a3,a3_12):
    xi = np.zeros((len(t),4))
    x1 = np.asarray(init)
    for i in range(len(t)-1):
        ti = t[i]
        tip1 = t[i+1]
        ht = tip1 - ti
        x2 = two_dof_cubic(x1             ,ti      ,F0,omg,m,c,k1,k2,a3,a3_12)
        x3 = two_dof_cubic(x1 + (ht/2)*x2 ,ti+ht/2 ,F0,omg,m,c,k1,k2,a3,a3_12)
        x4 = two_dof_cubic(x1 + (ht/2)*x3 ,ti+ht/2 ,F0,omg,m,c,k1,k2,a3,a3_12)
        x5 = two_dof_cubic(x1 + (ht)*x4   ,ti+ht   ,F0,omg,m,c,k1,k2,a3,a3_12)
        x6 = x1 + (ht/6)*(x2 + 2*x3 + 2*x4 + x5)
        x1 = x6
        xi[i,0] = x6[0]
        xi[i,1] = x6[1]
        xi[i,2] = x6[2]
        xi[i,3] = x6[3]
    return xi

"""
Main
"""
# Model Parameters
m = 1e0
c = 0.01e0
k1 = 9.0e0
k2 = 81
a3=0.005e2
a3_12 = 0.005e2
init = [0.0, 0.0, 0.0, 0.0] # initial conditions

# Time Span and time resolution
t = np.linspace(0.0e0, 500.0e0, 3001)

# Sweep Parameters
Forces = np.linspace(0.05e0, 0.1e0, 1)
omgf = np.linspace(2.5e0, 3.6e0, 500)
omgb = np.flip(omgf)

# Directory to save data
directory = "Force_sweep"

# Force Sweep loop
for F0 in Forces:
    fileName = "Force_"+'{:05.2f}'.format(F0) # name of data file
    if not os.path.exists(directory):
        os.makedirs(directory)
    tPath=os.path.join(directory,fileName)
    outF = open(tPath+"_FWD.csv","w")
    outF.write("Omega   Amp1    Amp2\n")
    outB = open(tPath+"_BKW.csv","w")
    outB.write("Omega   Amp1    Amp2\n")
    # Forward Sweep
    for omgs in omgf:
        x = my_RK4(two_dof_cubic, init, t,F0,omgs,m,c,k1,k2,a3,a3_12)
        tp = 2*np.pi/omgs
        tmin = max(t)-10*tp
        amp1f = max(x[t>tmin,0])
        amp2f = max(x[t>tmin,2])
        outF.write(("{0:5.5f}{1:10.10f}{2:10.10f}\n").format(omgs,amp1f,amp2f))
        init = [amp1f,0,0,0]
    # Backward Sweep
    for omgs in omgb:
        x = my_RK4(two_dof_cubic, init, t,F0,omgs,m,c,k1,k2,a3,a3_12)
        tp = 2*np.pi/omgs
        tmin = max(t)-10*tp
        amp1b = max(x[t>tmin,0])
        amp2b = max(x[t>tmin,2])
        outB.write(("{0:5.5f}{1:10.10f}{2:10.10f}\n").format(omgs,amp1b,amp2b))
        init = [amp1b,0,0,0]
    outF.close()
    outB.close()