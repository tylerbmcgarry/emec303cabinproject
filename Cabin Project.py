# -*- coding: utf-8 -*-
"""
Created on Tue Feb 25 13:29:12 2025

@author: sam and tyler
"""
import numpy as np
import matplotlib.pyplot as plt

def euler_multi(f1, f2, x0,y1_0, y2_0, h, L):
    
    
    N = int(L/h) # number of steps
    
    x = np.zeros(N+1)
    y1 = np.zeros(N+1)
    y2 = np.zeros(N+1)
    
    x[0] = x0
    y1[0] = y1_0
    y2[0] = y2_0
    #build x and y arrays
    
    for i in range(N):
        x[i+1] = x[i] + h
        y1[i+1] = y1[i] + h * f1(x[i], y1[i], y2[i]) #euler method
        y2[i+1] = y2[i] + h * f2(x[i], y1[i], y2[i])
    return x,y1,y2

#%% Constants
C_air = 1.005e3  # J/(kg*K)
rho_air = 1.2  # kg/m3
k1 = 0.2  # W/(m2*K) 
k3 = 0.5  # W/(m2*K)
Q0 = 500  # W/m2

# Pine Constants and Oak Constants
kpine = 0.5 #h^-1
npine = 2.1 

koak = 0.2 #h^-1
noak = 1.4 

V_bot = 50 #volume bottom room
V_top = 30  # Top Room volumes in m3

t_vals = np.linspace(0, 48*3600, 100)


#%% Room Dimensions

L1 = 5
W = 6
H1 = 3 
H2 = 3

A_fire = 0.5  # Fireplace area in m2
A_walls =  L1*H1*2 + W*H1*2 #all side walls of cabin in m2
A_roof =  0.5*H2*W*2 + L1*H2*np.sqrt(2)*2 #roof area in m2
A_ceil = L1*W  #ceiling area in m2

#%% Outside Temperature Function
def T_out(t):
    return -10 * np.sin((2 * np.pi * t) / 86400)

# plt.figure(1)
# plt.plot(t_vals/3600, T_out(t_vals))
# plt.show()

#%% Fireplace : Qin


def Q_fire(add_log, add_time, t):
    t = np.asarray(t)
    log_time = t - add_time  # Time since log was added
    return np.where(log_time >= 0, add_log * 500 * (1 + koak * t/3600) ** -noak, 0)


# Define Qin outside of Q_fire
def Qin(t):
    return (
        Q_fire(6, 0 * 3600, t) +
        Q_fire(4, 4 * 3600, t) +
        Q_fire(20, 12 * 3600, t) +
        Q_fire(50, 18 * 3600, t) +
        Q_fire(1, 24 * 3600, t)
    )


def Q_in_fire(t):
    return Qin(t)
    
    

#plt.figure(2)
#plt.plot(t_vals/3600, Q_in_fire(t_vals))
#plt.show()


#%%System of ODE's

k2 = 10 

def dT1_dt(t,T1,T2):  
   
    Q_wall = k1 * (T1 - T_out(t))
    Q_ceil = k2 * (T1 - T2)
    
    return ( (-A_walls * Q_wall) - (A_ceil * Q_ceil) + (A_fire * Q_in_fire(t)) ) / (C_air * rho_air * V_bot)

def dT2_dt(t,T1,T2):
    
    Q_ceil = k2 * (T1 - T2)
    Q_roof = k3 * (T2 - T_out(t))
    
    return ( (A_ceil * Q_ceil) - (A_roof * Q_roof) ) / (C_air * rho_air * V_top)


#%% Euler
# Initial conditions
T1_0 = 7
T2_0 = 5

t_start = 0 
t_end = 2 * 86400  # Simulate for 2 days (in seconds)
t_eval = np.linspace(t_start, t_end)  # Time points for evaluation
t0 = t_start
h = 60
L = t_end
time, T1, T2 = euler_multi(dT1_dt, dT2_dt, t0, T1_0, T2_0, h, L)

#%% Plot
plt.figure(figsize=(10, 8))  # Create the figure


low_limit = 12.5
high_limit = 27



# First subplot - Cabin Temperature Over Time
#plt.subplot(2, 1, 1)
plt.figure(3)
plt.plot(time / 3600, T1, label="T1 (Downstairs)")
plt.plot(time / 3600, T2, label="T2 (Upstairs)")
plt.plot(t_vals / 3600, T_out(t_vals), 'g-', label='Temp Outside')

plt.axhline(y=low_limit, color='black', linestyle='--', label='Temp Boundary')
plt.axhline(y=high_limit, color='black', linestyle='--',)

plt.ylabel("Temperature (°C)")
plt.legend()
plt.title("Cabin Temperature Over Time")
plt.grid()

# Second subplot - Fire Energy Over Time
#plt.subplot(2, 1, 2)  # 2 rows, 1 column, second subplot
# plt.plot(time / 3600, Q_in_fire(time) / 50, 'r-', label='Energy of Fire')

# plt.xlabel("Time (hours)")
# plt.ylabel("Fire Energy (kW)")
# plt.legend()
# plt.title("Fire Energy")
# plt.grid()

























