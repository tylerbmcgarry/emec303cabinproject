# -*- coding: utf-8 -*-
"""
Created on Tue Feb 25 13:29:12 2025

@author: tyler
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

# Constants
C_air = 1.005e3  # J/(kg*K)
rho_air = 1.2  # kg/m³
k1 = 0.2  # W/(m²*K)
k3 = 0.5  # W/(m²*K)
A1, A2, A3 =   # Areas in m² 
V1, V2 = 50, 30  # Room volumes in m³ 
A_in = 0.5  # Fireplace area in m²
k2 = 5  # Initial value, adjustable

# Outside temperature function (Eq. 4)
def T_out(t):
    return -10 * np.sin(2 * np.pi * t / 86400)

# Fireplace
def Q_in(t):
    return 0  #fire control

# ODE system 
def cabin_ode(t, T):
    T1, T2 = T  # Unpack temperatures
    Tout = T_out(t)
    
    # Heat fluxes
    Q1 = k1 * (T1 - Tout)
    Q2 = k2 * (T1 - T2)
    Q3 = k3 * (T2 - Tout)
    
    # Temperature derivatives
    dT1_dt = (-A1 * Q1 - A2 * Q2 + A_in * Q_in(t)) / (C_air * rho_air * V1)
    dT2_dt = (A2 * Q2 - A3 * Q3) / (C_air * rho_air * V2)
    
    return [dT1_dt, dT2_dt]

# Initial conditions
T1_0, T2_0 = 7, 5  # Initial temperatures in Celsius
t_start, t_end = 0, 2 * 86400  # Simulate for 2 days (in seconds)
t_eval = np.linspace(t_start, t_end, 1000)  # Time points for evaluation

# Solve the ODEs
sol = solve_ivp(cabin_ode, [t_start, t_end], [T1_0, T2_0], t_eval=t_eval, method='RK45')

# Plot results
plt.figure(figsize=(10, 5))
plt.plot(sol.t / 3600, sol.y[0], label="T1 (Downstairs)")
plt.plot(sol.t / 3600, sol.y[1], label="T2 (Upstairs)")
plt.xlabel("Time (hours)")
plt.ylabel("Temperature (°C)")
plt.legend()
plt.title("Cabin Temperature Over Time (No Fire)")
plt.grid()
plt.show()
