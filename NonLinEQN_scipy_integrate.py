import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from drawnow import drawnow, figure

def dampfunc(x, t,F0,omg,m,c,k,a3):
    x0 = x[0]
    x1 = x[1]
    x2 = F0 * np.cos(omg * t)/m - c * x1 / m - k * x0 / m - a3 * x0**3/m
    return [x1, x2]
def save2Fold_timeseries_csv(data,directory,fileName):
    if not os.path.exists(directory):
        os.makedirs(directory)
    tPath=os.path.join(directory,fileName)
    np.savetxt(tPath+'.csv',data,delimiter=',')
def draw_fig(omgs, amp):
    #f2 = plt.figure(1)
    #plt.plot(omgs, amp, 'ro')
    plt.show()

F0 = 5
omg = np.linspace(1.5, 4.5, 100)
m = 1.0
c = 0.0001
k = 9.0
a3=0.005
dataFold = 'my_data_fold'

init = [0.0, 0.0]
t = np.linspace(0.0, 10.0, 501)
plt.ion()
plt.figure()
for omgs in omg:
    x = odeint(dampfunc, init, t,args=(F0, omgs,m,c,k,a3),mxordn=100)
    tp = 2*np.pi/omgs
    tmin = max(t)-10*tp
    amp = max(x[t>tmin,0])
    f2 = plt.figure(1)
    plt.plot(omgs, amp, 'ro')
    plt.axis([1.5,4.5,0,2])
    plt.pause(0.001)
    #plt.show()
    #drawnow(draw_fig(omgs, amp))
#    plt.hold(True)
#    plt.show()
#save2Fold_timeseries_csv(x,'my_data','timeseries')

plt.show()
#f2.show()
