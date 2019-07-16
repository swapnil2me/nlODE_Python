import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from drawnow import drawnow, figure


# [x1,dx1,x2,dx2]
def dampfunc(x, t,F0,omg,m,c,k1,k2,a3,a3_12):
    x1 = x[1]
    x1d = F0 * np.cos(omg * t)/m - c * x[1] / m - k1 * x[0] / m - a3 * x[0]**3/m - a3_12*(x[0]-x[2])**3
    x2 = x[3]
    x2d =  - c * x[3] / m - k2 * x[2] / m - a3_12*(x[2]-x[0])**3
    return np.asarray([x1, x1d, x2, x2d])
    
def my_RK4(dampfunc, init, t,F0,omg,m,c,k1,k2,a3,a3_12):
    xi = np.zeros((len(t),4))
    x1 = np.asarray(init)
    for i in range(len(t)-1):
        ti = t[i]
        tip1 = t[i+1]
        ht = tip1 - ti

        x2 = dampfunc(x1             ,ti      ,F0,omg,m,c,k1,k2,a3,a3_12)
        x3 = dampfunc(x1 + (ht/2)*x2 ,ti+ht/2 ,F0,omg,m,c,k1,k2,a3,a3_12)
        x4 = dampfunc(x1 + (ht/2)*x3 ,ti+ht/2 ,F0,omg,m,c,k1,k2,a3,a3_12)
        x5 = dampfunc(x1 + (ht)*x4   ,ti+ht   ,F0,omg,m,c,k1,k2,a3,a3_12)

        x6 = x1 + (ht/6)*(x2 + 2*x3 + 2*x4 + x5)
        x1 = x6
        xi[i,0] = x6[0]
        xi[i,1] = x6[1]
        xi[i,2] = x6[2]
        xi[i,3] = x6[3]

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

Forces = np.linspace(0.05e0, 0.1e0, 5)
omgf = np.linspace(2.5e0, 3.6e0, 500e0)
omgb = np.flip(omgf)
m = 1.0e0
c = 0.01e0
k1 = 9.0e0
k2 = 81
a3=0.005e2
a3_12 = 0.005e2
dataFold = 'my_data_fold'

init = [0.0, 0.0, 0.0, 0.0]
t = np.linspace(0.0, 500.0e0, 3001e0)
plt.ion()
plt.figure()

directory = "Force_freq_sweep"

for F0 in Forces:
    fileName = "Force_"+'{:05.2f}'.format(F0)
    if not os.path.exists(directory):
        os.makedirs(directory)
    tPath=os.path.join(directory,fileName)


    outF = open(tPath+"_FWD.csv","w")
    outF.write("Omega   Amp1    Amp2\n")
    outB = open(tPath+"_BKW.csv","w")
    outB.write("Omega   Amp1    Amp2\n")
    for omgs in omgf:
        x = my_RK4(dampfunc, init, t,F0,omgs,m,c,k1,k2,a3,a3_12)
        tp = 2*np.pi/omgs
        tmin = max(t)-10*tp
        amp1f = max(x[t>tmin,0])
        amp2f = max(x[t>tmin,2])

        outF.write(("{0:10.5f}{1:10.5f}{2:10.5f}\n").format(omgs,amp1f,amp2f))

        init = [amp1f,0,0,0]
        f2 = plt.figure(1)
        plt.plot(omgs, amp1f, 'ro')
        plt.plot(omgs, amp2f, 'ko')
        #plt.axis([1.5,10.5,0,25])
        #plt.plot(t,x[:,0])
        plt.pause(0.000000001)
        #plt.show()
        #drawnow(draw_fig(omgs, amp))
        #plt.cla()
        #plt.show()
    #save2Fold_timeseries_csv(x,'my_data','timeseries')
    for omgs in omgb:
        x = my_RK4(dampfunc, init, t,F0,omgs,m,c,k1,k2,a3,a3_12)
        tp = 2*np.pi/omgs
        tmin = max(t)-10*tp
        amp1b = max(x[t>tmin,0])
        amp2b = max(x[t>tmin,2])

        outB.write(("{0:10.5f}{1:10.5f}{2:10.5f}\n").format(omgs,amp1b,amp2b))

        init = [amp1b,0,0,0]
        f2 = plt.figure(1)
        plt.plot(omgs, amp1b, 'bo')
        plt.plot(omgs, amp2b, 'yo')
        plt.pause(0.000000001)
    plt.show()
    #f2.show()
    outF.close()
    outB.close()