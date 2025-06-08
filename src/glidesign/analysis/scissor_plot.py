# Python native imports

# Python third party imports
import numpy as np
import matplotlib.pyplot as plt

# ParaPy imports

# Custom imports
#from glidesign.geometry import Glider

#TODO: make it automatically take inputs from glidesign and outputs from aerodynamic analysis

#Designer input
SM = 0.05                     #Stability margin [-]

#Geometry parameters
x_ac = 0.25                   #Location aerodynamic centre w.r.t. MAC [-]
s_h = 1.4                     #Horizontal tail surface [m^2]
s = 10.5                      #Wing surface area  [m^2]
l_h = 5                       #Tail length (distance AC wing - AC horizontal tail) [m]
c = 1.2                       #Mean aerodynamic wing chord length [m]

#Aerodynamic parameters
cm_ac = -0.05                 #Pitching moment coefficient around AC [-]

cl_a_min_h = 0.7              #Lift coefficient of airplane MINUS horizontal tail [-]
cl_h = -0.2                   #Lift coefficient of horizontal tail [-]

cl_alpha_h = 3.5              #Lift curve slope coefficient of horizontal tail [-]
cl_alpha_a_min_h = 5.5        #Lift curve slope coefficient of airplane MINUS horizontal tail [-]

de_da = 0.4                   #Rate of change of downwash angle at horizontal tail [-]

velocity = 55                 #Velocity in the freestream [m/s]
velocity_h = 54               #Velocity at the horizontal tail [m/s]

#Generate limit curves
def x_cg_stability_limit(sh_s):
    return x_ac + ((cl_alpha_h / cl_alpha_a_min_h) * (1 - de_da) * sh_s * (l_h / c) * (velocity_h / velocity) ** 2) - SM

def x_cg_controllability_limit(sh_s):
    return x_ac - (cm_ac / cl_a_min_h) + ((cl_h / cl_a_min_h) * sh_s * (l_h /c) * (velocity_h / velocity)**2)

s_h_s_range = np.linspace(0.005, 1, 100)
x_cg_stability = x_cg_stability_limit(s_h_s_range)
x_cg_controllability = x_cg_controllability_limit(s_h_s_range)

current_sh_s = s_h / s

#Plot
plt.figure(figsize=(10,6))
plt.plot(x_cg_stability, s_h_s_range, label='Stability Limit', color='blue')
plt.plot(x_cg_controllability, s_h_s_range, label='Controllability Limit', color='red')
plt.axhline(current_sh_s, color='green', linestyle='--', label=f'Current Sh/S = {current_sh_s:.2f}')
plt.xlabel('Xcg / MAC [-]')
plt.ylabel('Sh/S [-]')
plt.title('Scissor Plot')
plt.grid(True)
plt.legend()
plt.show()

