# Pouya Taghipour    -----> 9933014
# Amirhossein Mohebi -----> 9933057

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

# The time-dependent behavior of the potassium conductance can be simulated using a program that integrates the equation for

f = 5  # Hz
t = np.arange(0, 1/f * 1e+3 + 1)  # ms
E_field = 1  # v/m
v = 10  # mv

# Function for the first ODE system
def dXdT_n(t, n, v):
    a_n = 0.01 * (10 - v) / (np.exp((10 - v) / 10) - 1)
    b_n = 0.125 * np.exp(-v / 80)
    dn_dt = a_n * (1 - n) - b_n * n
    return [dn_dt]

# Voltage-clamp simulation for potassium current
v = 10
xo = 0
g_K = 10
sol = solve_ivp(dXdT_n, [0, 12], [xo], args=(v,), method='BDF')
t = sol.t
x = sol.y.T
g = g_K * (x[:, 0] ** 4)

plt.figure()
plt.plot(t, g)
plt.xlabel('t (ms)')
plt.ylabel('g_{K}')
plt.gca().tick_params(axis='both', labelsize=16)

# Data for v=10 mV
data2 = np.array([
    [0.1640, 0.3330, 0.5860, 0.7540, 1.1100, 1.4900, 2.0000, 2.8400, 4.1900, 6.3800, 8.8000, 11.2000],
    [0.2500, 0.2600, 0.2600, 0.3100, 0.3100, 0.3200, 0.4000, 0.5000, 0.7100, 1.0000, 1.3000, 1.6000]
])

v = 10.001
xo = 0.28
sol = solve_ivp(dXdT_n, [0, 12], [xo], args=(v,), method='BDF')
t = sol.t
x = sol.y.T
g = g_K * (x[:, 0] ** 4)

plt.figure()
plt.plot(t, g, 'k-', linewidth=1.5)
plt.plot(data2[0, :], data2[1, :], 'ko', linewidth=1.5, markerfacecolor=[1, 1, 1], markersize=5)
plt.gca().tick_params(axis='both', which='both', labelsize=9, bottom=False, left=False)
plt.box(False)
plt.axis([0, 12, 0, 2.4])
plt.text(12.2, 0.45*2.4, '$v = 10$~mV', fontsize=10)

# Function for the second ODE system
def dXdT_mh(t, x, v):
    m, h = x
    a_m = 0.1 * (25 - v) / (np.exp((25 - v) / 10) - 1)
    b_m = 4 * np.exp(-v / 18)
    a_h = 0.07 * np.exp(-v / 20)
    b_h = 1 / (np.exp((30 - v) / 10) + 1)
    dm_dt = a_m * (1 - m) - b_m * m
    dh_dt = a_h * (1 - h) - b_h * h
    return [dm_dt, dh_dt]

# Voltage-clamp simulation for sodium current
v = 10
xo = [0, 1]
g_Na = 10
sol = solve_ivp(dXdT_mh, [0, 12], xo, args=(v,), method='BDF')
t = sol.t
x = sol.y.T
g = g_Na * (x[:, 0] ** 3) * x[:, 1]

plt.figure()
plt.plot(t, g)
plt.xlabel('t (ms)')
plt.ylabel('g_{Na}')
plt.gca().tick_params(axis='both', labelsize=16)

# Data for v=10 mV
data2 = np.array([
    [0.1490, 0.3400, 0.5100, 0.7010, 1.1000, 1.4900, 1.9700, 2.8000, 4.1600, 6.4100, 8.7700, 11.3000],
    [0.0400, 0.0900, 0.1100, 0.1200, 0.1200, 0.1200, 0.1100, 0.1100, 0.1000, 0.0900, 0.0800, 0.0800]
])

v = 10.0
xo = [0.05, 0.7]
sol = solve_ivp(dXdT_mh, [0, 12], xo, args=(v,), method='BDF')
t = sol.t
x = sol.y.T
g = g_Na * (x[:, 0] ** 3) * x[:, 1]

plt.figure()
plt.plot(t, g, 'k-', linewidth=1.5)
plt.plot(data2[0, :], data2[1, :], 'ko', linewidth=1.5, markerfacecolor=[1, 1, 1], markersize=5)
plt.gca().tick_params(axis='both', which='both', labelsize=9, bottom=False, left=False)
plt.box(False)
plt.axis([0, 12, 0, 0.15])
plt.text(12.2, 0.30*0.15, '$v = 10$~mV', fontsize=10)

# Function for the Hodgkin-Huxley model
def dXdT_HH(t, x, I_app):
    v, m, n, h = x
    V_Na = 115
    V_K = -12
    V_L = 10.6
    g_Na = 120
    g_K = 36
    g_L = 0.3
    C_m = 1e-6
    
    a_m = 0.1 * (25 - v) / (np.exp((25 - v) / 10) - 1)
    b_m = 4 * np.exp(-v / 18)
    a_h = 0.07 * np.exp(-v / 20)
    b_h = 1 / (np.exp((30 - v) / 10) + 1)
    a_n = 0.01 * (10 - v) / (np.exp((10 - v) / 10) - 1)
    b_n = 0.125 * np.exp(-v / 80)
    
    I_Na = (m ** 3) * h * g_Na * (v - V_Na)
    I_K = (n ** 4) * g_K * (v - V_K)
    I_L = g_L * (v - V_L)
    
    dv_dt = (-I_Na - I_K - I_L + I_app) / C_m
    dm_dt = a_m * (1 - m) - b_m * m
    dn_dt = a_n * (1 - n) - b_n * n
    dh_dt = a_h * (1 - h) - b_h * h
    
    g_Na_out = (m ** 3) * h * g_Na
    g_K_out = (n ** 4) * g_K
    g_L_out = g_L
    
    return [dv_dt, dm_dt, dn_dt, dh_dt], [g_Na_out, g_K_out, g_L_out]

# Simulation of Hodgkin-Huxley model with applied current
I_app = 0
sol = solve_ivp(lambda t, x: dXdT_HH(t, x, I_app)[0], [0, 30], [0, 0, 0, 0], method='BDF')
xo = sol.y[:, -1]

# Now Adding nonzero applied current:
I_app = 6.2
sol = solve_ivp(lambda t, x: dXdT_HH(t, x, I_app)[0], [0, 30], xo, method='BDF')
t = sol.t
x = sol.y.T

plt.figure()
plt.plot(t, x[:, 0], 'k-', linewidth=1.5)
plt.gca().tick_params(axis='both', labelsize=18)
plt.box(True)
plt.axis([0, 30, -20, 120])
plt.xlabel('$t$ (ms)', fontsize=20)
plt.ylabel('$v$ (mV)', fontsize=20)

# Obtaining and plotting conductivity values
G3 = np.array([dXdT_HH(0, xi, I_app)[1] for xi in x])

plt.figure()
plt.plot(t, G3[:, 0], 'k', linewidth=1.5)
plt.plot(t, G3[:, 1], 'k--', linewidth=1.5)
plt.gca().tick_params(axis='both', labelsize=18)
plt.box(True)
plt.axis([0, 30, 0, 45])
plt.xlabel('$t$ (ms)', fontsize=20)
plt.ylabel('Conductance (mS$\cdot$cm$^{-2}$)', fontsize=20)
plt.text(16.4, 21, '$g_{Na}$', fontsize=20)
plt.text(19, 8, '$g_{K}$', fontsize=20)

plt.show()
