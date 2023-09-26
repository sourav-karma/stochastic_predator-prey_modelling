import numpy as np
import matplotlib
matplotlib.use('tkagg')
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from sympy import sin, cos, pi


y0 = [5,550]#initial value [cm swad grass, 1 cow = 550 kg]
t=np.linspace(0,1200,num=1000)

r = 0.02;     #grass growth rate (day−1)  #when changed to 2, high change in graph, cows are increasing  
c = 7.2e-4;   #Maximal grass intake (cm·kg metabolic live weight−1·day−1) 
#ha=5.2e-4;
H = 2.8;     #Half saturation grass intake [cm]
b =81;      #Grass conversion efficiency [kg metabolic liveweight·cm−1]
l = 5.3e-4; #Grass requirements for animal maintenance [cm·kg metabolic liveweight−1·day−1]
s = 1;      #stocking rate
kmax = 121;    #Maximal grass carrying capacity 
kmin = 5.45;   #MINImal grass carrying capacity
CC = 1;        #climate factor, remains const.


params = [r,c,H,b,l,s,kmax,kmin,CC]


def sim(variables, t, params):
    x = variables[0]
    y = variables[1]
    
    r = params[0]
    c = params[1]
    H = params[2]
    b = params[3]
    l = params[4]
    s = params[5]
    kmax = params[6]
    kmin = params[7]
    CC = params[8]
    
    k = ((kmax-kmin)/2 *cos((2*pi*t)/365) +  (kmax + kmin)/2 )*CC
   #if t>=150 and t<300: c =c/2 + ha
    dxdt = r*x*(1-x/k) - ((c*s*x**2*y**0.75)/(H**2+x**2))
    dydt = b*y**0.75*((c*x**2/(H**2+x**2))- l) 
    dzdt = [dxdt,dydt]
    return dzdt
    
    
y=odeint(sim, y0, t, args=(params,))

f,(ax1,ax2) = plt.subplots(2)
line1, = ax1.plot(t,y[:,0],color="b")
line2, = ax2.plot(t,y[:,1],color="r")
ax1.set_ylabel("Grass cm")
ax2.set_ylabel("cows")
ax2.set_xlabel("time")
plt.savefig('./predator_prey_nigeria_model.png')
plt.show()


##state diagram
z1 = odeint(sim,[5,550],t, args=(params,))
z2 = odeint(sim,[7,550],t, args=(params,))
z3 = odeint(sim,[9,550],t, args=(params,))
plt.figure
plt.plot(z1[:,0],z1[:,1],'b',label='prey/predator: initial=5/550')
plt.plot(z2[:,0],z2[:,1],'r--',label='prey/predator: initial=7/550')
plt.plot(z3[:,0],z3[:,1],'k-.',label='prey/predator: initial=9/550')
ax = plt.gca()
plt.xlabel('grass (cm)')
plt.ylabel('Cow (kg/animal)')
plt.legend(loc='best')
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
plt.savefig('./predator_prey_nigeria.png')
plt.show()