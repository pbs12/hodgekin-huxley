#PRANAV SHARMA (pbs12)

import numpy as np
import matplotlib.pyplot as plt

gNa= 120 #mS/cm^2
# maximum conductance of sodium per unit area 
# used in time-dependent equation for sodium conductance                            

gK=  36 #mS/cm^2        
# maximum conductance of potassium
# used in time-dependent equation for sodium conductance                            
                      
GO= 0.3 #mS/cm^2

# conductance of all channels and ion movements other than sodium and potassium
                            
Cm= 1 #µF/cm^2 
# base membrane capacitance in microFarads per unit area
                             
Vr= -1*60 #mV 
# resting membrane potential of squid neuron 
                            
ENa= 115 + Vr #mV
# Nernst potential for Na+ for squid neuron
                      
EK= Vr - 12 #mV                        
# Nernst potential for K+ for squid neuron

EO= 10.613 + Vr  #mV 
# Nernst potential for all other ions for squid neuron
                     
                             
M0= .0529324853
# quantity between 0 & 1 that represents sodium gate activation at t = 0
# changes over time as modeled by differential equation dm/dt, affecting sodium conductance
                         
N0= .3176769141 
# quantity between 0 & 1 that represents potassium gate activation at t = 0
# changes over time as modeled by differential equation dm/dt, affecting potassium conductance
                      
H0= .5961207535                         

# quantity between 0 & 1 that represents sodium gate inactivation at t = 0
# changes over time as modeled by differential equation dh/dt, affecting sodium conductance

Istim= 32
# stimulus current provided to neuron, independent input of experiment                  
Istim2= 51
# second stimulus current provided to neuron, independent input of experiment                  

t0=0
# start time of simulation, t = 0
tEnd=40
# end time of simulation, t = 40
dt=0.01                         
#step of Forward Euler method function   

t=np.arange(t0,tEnd,dt)            
# code creating list object containing all times multiples of dt between t0 and end Time

n=int((tEnd-t0)/dt)
# determine number of time intervals                 

V=np.zeros(n)                      
N=np.zeros(n)                       
M=np.zeros(n)                       
H=np.zeros(n)                       
GNa=np.zeros(n)                     
GK=np.zeros(n)                      
INa=np.zeros(n)                     
IK=np.zeros(n)                      
IO=np.zeros(n)                      
# make arrays for each time-dependent variable


V[0]=Vr                            
M[0]=M0                            
N[0]=N0                            
H[0]=H0   
#set intial conditions for membrane potential, m, n, and h                         

for i in range(0,len(t)-1):
    
    GNa[i]=gNa*H[i]*M[i]**3 
    #calculate sodium conductance for time interval i
    # GNa(max) * h * m^3
    GK[i]=gK*N[i]**4                
    #calculate potassium conductance for time interval i
    #GK*max)*n^4

    INa[i]=GNa[i]*(V[i]-ENa)  
    #calculate sodium current for time interval i
    IK[i]=GK[i]*(V[i]-EK)    
    #calculate potassium current for time interval i     
    IO[i]=GO*(V[i]-EO)   
    #calculate other current for time interval i          
    
    vm=V[i]-Vr                     
    #calculate relative membrane potential for i 
    
    am=(0.1*(25-vm))/(np.exp((25-vm)/10)-1)
    bm=4*np.exp(-vm/18)
    an=(0.01*(10-vm))/(np.exp((10-vm)/10)-1)
    bn=0.125*np.exp(-vm/80)
    ah=0.07*np.exp(-vm/20)
    bh=1/(np.exp((30-vm)/10)+1)
    #calculate alpha and beta for gating variables m,n,h using vm time interval i
    
    N[i+1]=N[i]+dt*(an*(1-N[i])-bn*N[i]);            
    M[i+1]=M[i]+dt*(am*(1-M[i])-bm*M[i]);            
    H[i+1]=H[i]+dt*(ah*(1-H[i])-bh*H[i]); 
    #calculate m,n,h for next time interval so forward Euler process proceeds           
    
    
    if t[i]>=5 and t[i]<=5.2:                              
        V[i+1]=V[i]+dt*(-INa[i]-IK[i]-IO[i]+Istim)/Cm    
    elif t[i]>=18 and t[i]<=18.4:                              
        V[i+1]=V[i]+dt*(-INa[i]-IK[i]-IO[i]+Istim2)/Cm    
    else:
        V[i+1]=V[i]+dt*(-INa[i]-IK[i]-IO[i])/Cm 

    # code that calculates voltage for next time interval based off currents
    # of each ion in present interval
    # adds stimulus current in time interval [7,7.2]         



plt.figure(1)
plt.subplot(231)    
plt.plot(t,V)
plt.xlabel('Time')
plt.ylabel('Membrane Potential (mV)')
plt.title('Membrane Potential vs. Time')

plt.subplot(232)    
plt.plot(t,N)
plt.xlabel('Time')
plt.ylabel('Potassium Activation Gate Probability')
plt.title('Potassium Activation Gate Probability vs. Time')

plt.subplot(233)    
plt.plot(t,M)
plt.xlabel('Time')
plt.ylabel('Sodium Activation Gate Probability')
plt.title('Sodium Activation Gate Probability vs. Time')

plt.subplot(234)    
plt.plot(t,H)
plt.xlabel('Time')
plt.ylabel('Sodium Inactivation Gate Probability')
plt.title('Sodium Inactivation Gate Probability vs. Time')

plt.subplot(235)    
plt.plot(t,GNa)
plt.xlabel('Time')
plt.ylabel('Sodium Conductance (mS/cm^2)')
plt.title('Sodium Conductance vs. Time')

plt.subplot(236)    
plt.plot(t,GK)
plt.xlabel('Time')
plt.ylabel('Potassium Conductance (mS/cm^2)')
plt.title('Potassium Conductance vs. Time')

"""

fig, ax = plt.subplots()
ax.plot(t, INa, '--', label='INa')
ax.plot(t, IK, '-.', label='IK')
ax.plot(t, IO, '-', label='IO')
plt.xlabel('Time')
plt.ylabel('Ion Current (µA)')
plt.title('Ion Current vs. Time')
leg = ax.legend();
    
"""
