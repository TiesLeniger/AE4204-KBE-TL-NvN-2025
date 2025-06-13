# Python native imports

# Python third party imports
import numpy as np
import matplotlib.pyplot as plt

# ParaPy imports
from parapy.core import Base, Input, Attribute

# Custom imports
#from ..geometry import Glider
#from ..geometry import LiftingSurface
#from glidesign.geometry import Glider

#TODO: make it automatically take inputs from glidesign and outputs from aerodynamic analysis

#Initiate (testing purposes)
#glider = Glider(name="glider_test")

class ScissorPlot(Base):
    #Designer input
    SM = 0.05                     #Stability margin [-]

    #Geometry parameters
    x_ac = 0.25                   #Location aerodynamic centre w.r.t. MAC [-]
    s_h = 1.4                     #Horizontal tail surface [m^2]
    s = 10.5                    #Wing surface area  [m^2]

    l_h = 0                     #Tail length (distance AC wing - AC horizontal tail) [m]
    chord = 1.2                   #Mean aerodynamic wing chord length [m]

    #Aerodynamic parameters (from aerodynamic analysis)
    cm_ac = -0.05                 #Pitching moment coefficient around AC [-]

    cl_a_min_h = 0.7              #Lift coefficient of airplane MINUS horizontal tail [-]
    cl_h = -0.2                   #Lift coefficient of horizontal tail [-]

    cl_alpha_h = 5                #Lift curve slope coefficient of horizontal tail [/rad]
    cl_alpha_a_min_h = 6          #Lift curve slope coefficient of airplane MINUS horizontal tail (assumed equal to cl_alpha_w) [/rad]

    velocity = 55                 #Velocity in the freestream [m/s]
    velocity_h = 55               #Velocity at the horizontal tail [m/s]

    #Datcom de/da estimation
    wingspan = 15                 #Glider wingspan [m]
    m_tv = 1                      #Distance between horizontal tail and vortex shed plane
    r = 2*l_h / wingspan          #Derived from tail length
    sweep_4c = np.deg2rad(10)     #Quarter chord sweep angle [rad]
    AR = 25                       #Aspect ratio [-]

    def k_e_sweep(self):
        a = 0.1124 + 0.1265*self.sweep_4c + 0.1766*self.sweep_4c**2
        return (a/self.r**2) + (0.1024/self.r) + 2

    def k_e_sweep_0(self):
        return (0.1124/self.r**2) + (0.1024/self.r) + 2

    def de_da_est(self):
        a = 1 + (self.r**2 / (self.r**2 + 0.7915 + 5.0734*self.m_tv**2))**0.3113
        b = 1 - np.sqrt(self.m_tv**2 / (1 + self.m_tv**2))
        c = 0.4876 / (np.sqrt(self.r**2 + 0.6319 + self.m_tv**2))
        d = self.r / (self.r**2 + self.m_tv**2)
        return (self.k_e_sweep()/self.k_e_sweep_0()) * (d*c + a*b) * (self.cl_alpha_a_min_h/(np.pi * self.AR))

    #Stability and Controllability limit curves
    def x_cg_stability_limit(self, sh_s):
        return self.x_ac + ((self.cl_alpha_h / self.cl_alpha_a_min_h) * (1 - self.de_da) * sh_s * (self.l_h / self.chord) * (self.velocity_h / self.velocity) ** 2) - self.SM

    def x_cg_controllability_limit(self, sh_s):
        return self.x_ac - (self.cm_ac / self.cl_a_min_h) + ((self.cl_h / self.cl_a_min_h) * sh_s * (self.l_h /self.chord) * (self.velocity_h / self.velocity)**2)

    #Define Sh/S range
    s_h_s_range = np.linspace(0.005, 1, 100)
    #x_cg_stability = x_cg_stability_limit(self, s_h_s_range)
    #x_cg_controllability = x_cg_controllability_limit(s_h_s_range)

    current_sh_s = s_h / s

    #Plot the scissor plot
    def plot_scissor_plot(self):
        plt.figure(figsize=(10,6))
        plt.plot(self.x_cg_stability_limit(self), self.s_h_s_range, label='Stability Limit', color='blue')
        plt.plot(self.x_cg_controllability_limit(self), self.s_h_s_range, label='Controllability Limit', color='red')
        plt.axhline(self.current_sh_s, color='green', linestyle='--', label=f'Current Sh/S = {self.current_sh_s:.2f}')
        plt.xlabel('Xcg / MAC [-]')
        plt.ylabel('Sh/S [-]')
        plt.title('Scissor Plot')
        plt.grid(True)
        plt.legend()
        plt.show()