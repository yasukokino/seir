# use this cell to call your function again to produce additional visualizations
# add more cells as needed
# write your function to implement the numerical simulation here!

import numpy as np
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 15})
plt.rcParams["figure.figsize"] = [8,5]

def SIR_Euler(b,k,d,e,initial_conds):
    t0 = 0 # t-start
    t_end = 300 # t-end
    days = range(0, 301)

    h = 1 # stepsize
    steps = int((t_end - t0)/h + 1) # number of steps

    # variables:
    t = np.linspace(t0, t_end, steps) # storing t values
    S = np.zeros(steps) # for storing S values
    E = np.zeros(steps) # for storing E values
    I = np.zeros(steps) # for storing I values
    R = np.zeros(steps) # for storing R values
    D = np.zeros(steps) # for storing D values

    # initial conditions:
    S[0] = initial_conds[0]
    E[0] = initial_conds[1]
    I[0] = initial_conds[2] 
    R[0] = initial_conds[3]
    D[0] = initial_conds[4]
    N = S[0] + E[0]+ I[0] + R[0] + D[0]

    for n in range(steps-1): # range(start, stop, step)
        S[n+1] = S[n] - h* (b*S[n]*I[n]/N) 
        E[n+1] = I[n] + h* ((b*S[n]*I[n]/N) - (d*E[n]))
        I[n+1] = I[n] + h* (d*E[n] - (k+e)*I[n])
        R[n+1] = R[n] + h* (k*I[n])
        D[n+1] = D[n] + h* (e*I[n])
    
    plt.style.use('ggplot')
    plt.plot(t,S,linewidth=2,label='S(t)')
    plt.plot(t,E,linewidth=2,label='E(t)')
    plt.plot(t,I,linewidth=2,label='I(t)')
    plt.plot(t,R,linewidth=2,label='R(t)')
    plt.plot(t,D,linewidth=2,label='D(t)')
    plt.xlabel('t [days]')
    plt.ylabel('S, E, I, R, D')
    plt.legend(loc='best')
    plt.show()

    
# parameters:
infection_rate = 0.085
recovery_rate = 0.062
incubation_period = 1/11
fatality_rate = 0.014 

# initial conditions: 
S0 = 1000
E0 = 10
I0 = 10
R0 = 0
D0 = 0
initial_vals = [S0,E0,I0,R0,D0]

# call the function to run the simulation
SIR_Euler(b=infection_rate, k=recovery_rate, d=incubation_period,
          e=fatality_rate, initial_conds=initial_vals)
