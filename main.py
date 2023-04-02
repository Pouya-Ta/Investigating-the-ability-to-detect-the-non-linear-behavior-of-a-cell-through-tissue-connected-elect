import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from sympy import symbols, solve

# declaring variables
global gK, gNa, gL, vK, vNa, vL, phi, C, Q

# assigning values
gK = 36  # unit: mS/cm^2
vK = -72  # unit: mV
gNa = 120  # unit: mS/cm^2
vNa = 55  # unit: mV
gL = 0.3  # unit: mS/cm^2
vL = -49.4  # unit: mV
C = 1  # unit: muF/cm^2
phi = 1
Q = 3  # Temperature coefficient-unitless ; phi = Q^((T-6.3)/10) ; here T is taken as 6.3.

tspan = np.array([0, 100])

# Defining opening rate(alpha) and closing rate(beta) of K channel(n), Na channel(m,h) :

an = lambda V: (-0.01*(V+50))/(np.exp(-(V+50)/10)-1)
bn = lambda V: 0.125*np.exp(-(V+60)/80)

am = lambda V: (-0.1*(V+35))/(np.exp(-(V+35)/10)-1)
bm = lambda V: 4*np.exp(-(V+60)/18)

ah = lambda V: 0.07*np.exp(-(V+60)/20)
bh = lambda V: 1/(np.exp(-(V+30)/10)+1)

# Plots of gating Probabilities and membrane potential with various current injection value :

for I in range(0, 11, 1):
    # defining hh model equation :

    def func(p, t):
        return [(1/C)*(I-(gK*p[1]**4*(p[0] - vK)) - (gNa*p[2]**3*p[3]*(p[0]-vNa))-(gL*(p[0]-vL))),
                (an(p[0])*(1-p[1])) - (bn(p[0])*p[1]),
                (am(p[0])*(1-p[2])) - (bm(p[0])*p[2]),
                (ah(p[0])*(1-p[3])) - (bh(p[0])*p[3])]

    Parameters = odeint(func, y0=[-60, 0.317, 0.0529, 0.596], t=tspan)
    
    plt.figure()
    plt.plot(tspan, Parameters[:,0], color='red')
    plt.xlabel('time')
    plt.ylabel('membrane Voltage (mV)')
    plt.title(f'Membrane potential/time, I = {I} μA/cm^2')
    plt.show()

    plt.figure()
    plt.plot(tspan, Parameters[:,1], label='n')
    plt.plot(tspan, Parameters[:,2], label='m')
    plt.plot(tspan, Parameters[:,3], label='h')
    plt.legend()
    plt.xlabel('t')
    plt.title(f'n, m & h/t, I = {I} μA/cm^2')
    plt.show()

#Effect of varying Leakage voltage on membrane potential :

vL_values = [-49.2, -49.3, -49.4, -49.45, -49.5]
for vL in vL_values:
    I = 0
    
    #defining hodgkin-huxley model equation :
    def func(p, t):
        return [(1/C)*(I-(gK*p[1]**4*(p[0] - vK)) - (gNa*p[2]**3*p[3]*(p[0]-vNa))-(gL*(p[0]-vL))),
                (an(p[0])*(1-p[1])) - (bn(p[0])*p[1]),
                (am(p[0])*(1-p[2])) - (bm(p[0])*p[2]),
                (ah(p[0])*(1-p[3])) - (bh(p[0])*p[3])]

    Parameters = odeint(func, y0=[-60, 0.317, 0.0529, 0.596], t=tspan)

    plt.plot(tspan, Parameters[:, 0])
    plt.xlabel('t')
    plt.ylabel('membrane Voltage (mV)')
    plt.title('Effect of varying Leakage voltage on membrane potential')
plt.legend(['El =-49.2 V', 'El=-49.3 V', 'El=-49.4 V', 'El=-49.45 V', 'El=-49.5 V'])
plt.show()

# Evaluating equilibrium points of HH equation :

for I in range(8, 13):
    # assigning sym variables :
    V1, n1, m1, h1 = symbols('V1 n1 m1 h1')

    # equilibrium points :
    eq1 = ((I-gK*n1**4*(V1 - vK) - gNa*m1**3*h1*(V1-vNa)-gL*(V1-vL))/C)
    eq2 = an(V1)*(1-n1) - bn(V1)*n1
    eq3 = am(V1)*(1-m1) - bm(V1)*m1
    eq4 = ah(V1)*(1-h1) - bh(V1)*h1

    equi = solve((eq1, eq2, eq3, eq4), (V1, n1, m1, h1))
    equi_v = float(equi[0][0])
    equi_n = float(equi[0][1])
    equi_m = float(equi[0][2])
    equi_h = float(equi[0][3])

    def func(p, t):
        V, n, m, h = p
        return [(1/C)*(I - (gK*n**4*(V - vK)) - (gNa*m**3*h*(V-vNa)) - (gL*(V-vL))),
                (an(V)*(1-n)) - bn(V)*n,
                (am(V)*(1-m)) - bm(V)*m,
                (ah(V)*(1-h)) - bh(V)*h]

    Parameters = odeint(func, y0=[equi_v, equi_n, equi_m, equi_h], t=tspan)

    plt.plot(tspan, Parameters[:, 0])
    plt.title(f'I={I} μA/cm^2')
    plt.xlabel('t')
    plt.ylabel('membrane Voltage (mV)')
    plt.show()
