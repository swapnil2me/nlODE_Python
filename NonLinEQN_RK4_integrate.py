import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from drawnow import drawnow, figure

def dampfunc(x, t,F0,omg,m,c,k,a3):
    x0 = x[0]
    x1 = x[1]
    x2 = F0 * np.cos(omg * t)/m - c * x1 / m - k * x0 / m - a3 * x0**3/m
    return np.asarray([x1, x2])
    
def my_RK4(dampfunc, init, t,F0,omg,m,c,k,a3):
    xi = np.zeros((len(t),2))
    x1 = np.asarray(init)
    for i in range(len(t)-1):
        ti = t[i]
        tip1 = t[i+1]
        ht = tip1 - ti

        x2 = dampfunc(x1             ,ti      ,F0,omg,m,c,k,a3)
        x3 = dampfunc(x1 + (ht/2)*x2 ,ti+ht/2 ,F0,omg,m,c,k,a3)
        x4 = dampfunc(x1 + (ht/2)*x3 ,ti+ht/2 ,F0,omg,m,c,k,a3)
        x5 = dampfunc(x1 + (ht)*x4   ,ti+ht   ,F0,omg,m,c,k,a3)

        x6 = x1 + (ht/6)*(x2 + 2*x3 + 2*x4 + x5)
        x1 = x6
        xi[i,0] = x6[0]
        xi[i,1] = x6[1]

    return xi

def save2Fold_timeseries_csv(data,directory,fileName):
    if not os.path.exists(directory):
        os.makedirs(directory)
    tPath=os.path.join(directory,fileName)
    np.savetxt(tPath+'.csv',data,delimiter=',')
def draw_fig(omgs, amp):
    #f2 = plt.figure(1)
    #plt.plot(omgs, amp, 'ro')
    plt.show()

F0 = 0.1e0
omgf = np.linspace(2.5e0, 3.6e0, 500e0)
omgb = np.flip(omgf)
m = 1.0e0
c = 0.01e0
k = 9.0e0
a3=0.005e2
dataFold = 'my_data_fold'

init = [0.0, 0.0]
t = np.linspace(0.0, 500.0e0, 3001e0)
plt.ion()
plt.figure()
for omgs in omgf:
    x = my_RK4(dampfunc, init, t,F0,omgs,m,c,k,a3)
    tp = 2*np.pi/omgs
    tmin = max(t)-10*tp
    amp = max(x[t>tmin,0])
    init = [amp,0]
    f2 = plt.figure(1)
    plt.plot(omgs, amp, 'ro')
    #plt.axis([1.5,10.5,0,25])
    #plt.plot(t,x[:,0])
    plt.pause(0.000000001)
    #plt.show()
    #drawnow(draw_fig(omgs, amp))
    #plt.cla()
#    plt.show()
#save2Fold_timeseries_csv(x,'my_data','timeseries')
for omgs in omgb:
    x = my_RK4(dampfunc, init, t,F0,omgs,m,c,k,a3)
    tp = 2*np.pi/omgs
    tmin = max(t)-10*tp
    amp = max(x[t>tmin,0])
    init = [amp,0]
    f2 = plt.figure(1)
    plt.plot(omgs, amp, 'bo')
    plt.pause(0.000000001)
plt.show()
#f2.show()
