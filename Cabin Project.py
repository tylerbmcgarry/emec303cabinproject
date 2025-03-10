# -*- coding: utf-8 -*-
"""
Created on Tue Feb 25 13:29:12 2025

@author: sam and tyler
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp


#%% Define Euler Function
def euler(f, x0, y0, h, L):

    
    N = int(L/h)
    
    x = np.zeros(N)
    y = np.zeros(N)
    
    x[0] = x0
    y[0] = y0
    
    #build arrays
    
    for i in range(N-1):
        x[i+1] = x[i] + h 
        y[i+1] = y[i] + h*f(x[i], y[i]) #euler method
        
    return x, y
        

#%% Constants
C_air = 1.005e3  # J/(kg*K)
rho_air = 1.2  # kg/m3
k1 = 0.2  # W/(m2*K)
k2 = 10  
k3 = 0.5  # W/(m2*K)
Q0 = 500  # W/m2

V_bot = 50 #volume bottom room
V_top = 30  # Top Room volumes in m3



#%% Room Dimensions

L = 5
W = 6
H1 = 3 
H2 = 3

A_fire = 0.5  # Fireplace area in m2
A_walls =  L*H1*2 + W*H1*2 #all side walls of cabin in m2
A_roof =  0.5*H2*W*2 + L*H2*np.sqrt(2)*2 #roof area in m2
A_ceil = L*W  #ceiling area in m2
#%%

# Outside temperature function (Eq. 4)
def T_out(t):
    return -10 * np.sin((2 * np.pi * t) / 86400   )





low_limit = 12.5
high_limit = 27







# Fireplace



# Pine Constants
kpine = 0.5 #h^-1
npine = 2.1 

#Oak Constants
koak = 0.2 #h^-1
noak = 1.4 

# Fireplace

logtime = [0, 6, 12, 18, 24, 30, 36, 38]

def Q_fire(t):#, t1, k, n):
    return 500 * np.sin((2 * np.pi * t) / 86400)+1700 #Q0*(1 + kpine*t)**(-npine)  #fire control
    
    #return np.heaviside()




# ODE system 
def cabin_ode(t, T):
    T1, T2 = T  # Unpack temperatures
    Tout = T_out(t)
    
    # Heat fluxes
    Q_wall = k1 * (T1 - Tout)
    Q_ceil = k2 * (T1 - T2)
    Q_roof = k3 * (T2 - Tout)
    
    # Temperature derivatives
    dT1_dt = ( (-A_walls * Q_wall) - (A_ceil * Q_ceil) + (A_fire * Q_fire(t)) ) / (C_air * rho_air * V_bot) ## rate of bot temp change wrt time
    
    dT2_dt = ( (A_ceil * Q_ceil) - (A_roof * Q_roof) ) / (C_air * rho_air * V_top) #rate of top temp change wrt time
    
    
    return [dT1_dt, dT2_dt]

 






# Initial conditions
T1_0 = 7
T2_0 = 5  # Initial temperatures in Celsius

t_start = 0 
t_end = 2 * 86400  # Simulate for 2 days (in seconds)
t_eval = np.linspace(t_start, t_end)  # Time points for evaluation

# Solve the ODEs
sol = solve_ivp(cabin_ode, [t_start, t_end], [T1_0, T2_0], t_eval=t_eval, method='RK45')



plt.figure(figsize=(10, 8))  # Create the figure

# First subplot - Cabin Temperature Over Time
plt.subplot(2, 1, 1)  # 2 rows, 1 column, first subplot
plt.plot(sol.t / 3600, sol.y[0], label="T1 (Downstairs)")
plt.plot(sol.t / 3600, sol.y[1], label="T2 (Upstairs)")
plt.plot(sol.t / 3600, T_out(t_eval), 'g-', label='Temp Outside')

plt.axhline(y=low_limit, color='black', linestyle='--', label='Low Temp Boundary')
plt.axhline(y=high_limit, color='black', linestyle='--', label='High Temp Boundary')

plt.ylabel("Temperature (°C)")
plt.legend()
plt.title("Cabin Temperature Over Time")
plt.grid()

# Second subplot - Fire Energy Over Time
plt.subplot(2, 1, 2)  # 2 rows, 1 column, second subplot
plt.plot(sol.t / 3600, Q_fire(t_eval) / 50, 'r-', label='Energy of Fire')

plt.xlabel("Time (hours)")
plt.ylabel("Fire Energy (kW)")
plt.legend()
plt.title("Fire Energy")
plt.grid()



































