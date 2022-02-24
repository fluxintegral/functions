# -*- coding: utf-8 -*-
"""

"""
import numpy as np
import matplotlib.pylab as plt
import fluids.atmosphere as atm

def newtonRaphson(V,a,rho,W,S,rho_star,prop_efficiency,CDo,Power_max):
    dV = 10;
    delta_v = 1e-4
    i = 0
    while abs(dV) > delta_v and i < 100:
        Cdo = cdzero(V,a,CDo)
        K = AK(V,a)
        f = f_metric(Power_max,V,a,rho,Cdo,K,W,S,rho_star,prop_efficiency)
        f_prime = (f_metric(Power_max,(V+delta_v),a, rho, cdzero(V+delta_v,a,CDo),AK(V+delta_v,a),W,S,rho_star,prop_efficiency)-f)/delta_v
        dV = f/f_prime
        V = V-dV
        root = V
        i += 1
    if i > 99:
        print("Newton Raphson Terminated on iteration limit")        
    return(root)
    
def cdzero(V,a,CDo):
    Mach = V/a
    if Mach < 0.4:
        Cdo = CDo
    elif Mach >= 0.4:
        Cdo = CDo + CDo*(Mach-0.4)**2
    return(Cdo)

def AK(V,a):
    Mach = V/a
    if Mach < 0.7:
        K = 0.045
    elif Mach > 0.7 and Mach < 0.95:
        K = 0.325 - np.sqrt(0.0625 - (Mach - 0.7)**2)
    elif Mach > 0.95:
        K = 0.325
    return(K)

def f_metric(Power_max,V,a,rho,Cdo,K,W,S,rho_star,prop_efficiency):
    # Thrust = Power_max*(rho/rho_star)*prop_efficiency/V # Turbo prop
    Thrust = Power_max*prop_efficiency/V # Electric motor
    Drag = (0.5*(Cdo)*rho*S*V**2) +(2*K*(W**2))/(rho*S*V**2)
    f = Thrust - Drag
    return(f)
    

# Sref = 10.4
# Qmax = .5*1.25*(280*0.514)**2 #% N/m2 Base on maximum kCAS of 280
# Mach_max = 0.6;
# CL_max = 1.39;
# air_star = atm.ATMOSPHERE_1976(0);
# rho_star = air_star.rho
# m = 3084  # kg min
# # m = 4536    # kg max
# prop_efficiency = .83
# CDo = 0.023*1.15 # Clean
# Power_max = 1196e3
# flight_ceiling = 30000/3.27

# Sref = 10
# Qmax = .5*1.25*(225/3.6)**2 #% N/m2 Based on maximum 225 km/h max
# Mach_max = 0.3;
# CL_max = 1.2;
# m = 300  # kg min
# prop_efficiency = .7
# CDo = 0.023*1.15 # Clean
# Power_max = 63e3-28*400
# flight_ceiling = 18e3/3.27

air_star = atm.ATMOSPHERE_1976(0);


g = 9.81
W = m*g 
altitude_m = np.arange(0, flight_ceiling, 100)

rho_star = air_star.rho

Vstall = []
VQ_max = []
VMach_max = []
V_low_root = []
V_high_root = []
for alt in altitude_m:
    
    air = atm.ATMOSPHERE_1976(alt)
    Vstall = np.append(Vstall,np.sqrt(2*m*g/(air.rho*Sref*CL_max)))
    VQ_max = np.append(VQ_max,np.sqrt((2*Qmax)/air.rho))
    VMach_max = np.append(VMach_max,Mach_max*air.v_sonic)
    
    V_high = VQ_max[-1]+150 # guess high V
    V_low = Vstall[-1]- 50   # guess low V    
    
    V_low_check = newtonRaphson(V_low, air.v_sonic, air.rho, W, Sref, rho_star, prop_efficiency, CDo, Power_max)
    V_high_check = newtonRaphson(V_high, air.v_sonic, air.rho, W, Sref, rho_star, prop_efficiency, CDo, Power_max)
    
    # Tripper Drop some apples on Newton to prevent upper and lower root being the same
    its = 0
    while np.abs(V_low_check - V_high_check) < 0.1 and its < 1000:
        V_low = V_low*.95
        V_low_check = newtonRaphson(V_low, air.v_sonic, air.rho, W, Sref, rho_star, prop_efficiency, CDo, Power_max)
        V_high_check = newtonRaphson(V_high, air.v_sonic, air.rho, W, Sref, rho_star, prop_efficiency, CDo, Power_max)
        its += 1
    print(its)    
        
    V_low_root = np.append(V_low_root,newtonRaphson(V_low, air.v_sonic, air.rho, W, Sref, rho_star, prop_efficiency, CDo, Power_max))
    V_high_root = np.append(V_high_root,newtonRaphson(V_high, air.v_sonic, air.rho, W, Sref, rho_star, prop_efficiency, CDo, Power_max))

    
plt.figure()
plt.plot(Vstall,altitude_m,'r')    
plt.plot(VMach_max,altitude_m)    
plt.plot(VQ_max,altitude_m)    
plt.plot([Vstall[-1],V_high_root[-1]],[altitude_m[-1],altitude_m[-1]])
plt.plot(V_high_root,altitude_m)    
plt.plot(V_low_root,altitude_m)    
plt.xlabel('V TAS')
plt.ylabel('Altitude TAS')
