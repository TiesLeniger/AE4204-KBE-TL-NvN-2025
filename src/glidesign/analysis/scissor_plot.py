# Python native imports

# Python third party imports
import numpy as np
import matplotlib.pyplot as plt

# ParaPy imports
from parapy.core import Base, Input, Attribute

# Custom imports

#TODO: make it automatically take inputs from aerodynamic analysis

class ScissorPlot(Base):
    #Designer input
    SM = Input()                            #Stability margin [-]

    #Center of gravity inputs
    actual_x_cog = Input()
    fwd_limit_x_cog = Input()               #maximum allowed pilot mass
    bwd_limit_x_cog = Input()               #minimum allowed pilot mass

    #Wing location
    wing_x_location = Input()

    #Geometry parameters
    hor_tail_span = Input()
    hor_tail_taper = Input()
    x_ac = Input(0.25)                      # Location aerodynamic centre w.r.t. MAC [-]
    s = Input()                             # Wing surface area  [m^2]

    l_h = Input()                           # Tail length (distance AC wing - AC horizontal tail) [m]
    chord = Input()                         # Mean aerodynamic wing chord length [m]

    #Aerodynamic parameters (from aerodynamic analysis)
    cm_ac = Input()                         # Pitching moment coefficient around AC [-]
    cl_a_min_h = Input()                    # Lift coefficient of airplane MINUS horizontal tail [-]
    cl_h = Input()                          # Lift coefficient of horizontal tail [-]
    cl_alpha_h = Input()                    # Lift curve slope coefficient of horizontal tail [/rad]
    cl_alpha_a_min_h = Input()              # Lift curve slope coefficient of airplane MINUS horizontal tail (assumed equal to cl_alpha_w) [/rad]
    Vh_V = Input(1)                         # Ratio of velocity at the horizontal tail to velocity at the main wing

    #Datcom de/da estimation
    wingspan = Input()                      # Glider wingspan [m]
    m_tv = Input()                          # Distance between horizontal tail and vortex shed plane
    sweep_4c = Input()                      # Quarter chord sweep angle [rad]
    AR = Input()                            # Aspect ratio [-]

    @Attribute
    def r(self):
        return 2 * self.l_h / self.wingspan          # Derived from tail length

    @Attribute
    def k_e_sweep(self):
        a = 0.1124 + 0.1265*self.sweep_4c + 0.1766*self.sweep_4c**2
        return (a/self.r**2) + (0.1024/self.r) + 2
    
    @Attribute
    def k_e_sweep_0(self):
        return (0.1124/self.r**2) + (0.1024/self.r) + 2

    @Attribute
    def de_da_est(self):
        a = 1 + (self.r**2 / (self.r**2 + 0.7915 + 5.0734*self.m_tv**2))**0.3113
        b = 1 - np.sqrt(self.m_tv**2 / (1 + self.m_tv**2))
        c = 0.4876 / (np.sqrt(self.r**2 + 0.6319 + self.m_tv**2))
        d = self.r / (self.r**2 + self.m_tv**2)
        return (self.k_e_sweep/self.k_e_sweep_0()) * (d*c + a*b) * (self.cl_alpha_a_min_h/(np.pi * self.AR))

    #Stability and Controllability limit curves
    def x_cg_stability_limit(self, sh_s):
        return self.x_ac + ((self.cl_alpha_h / self.cl_alpha_a_min_h) * (1 - self.de_da_est) * sh_s * (self.l_h / self.chord) * (self.velocity_h / self.velocity) ** 2) - self.SM

    def x_cg_stability_limit_no_sm(self, sh_s):
        return self.x_ac + ((self.cl_alpha_h / self.cl_alpha_a_min_h) * (1 - self.de_da_est) * sh_s * (self.l_h / self.chord) * (self.velocity_h / self.velocity) ** 2)

    def x_cg_controllability_limit(self, sh_s):
        return self.x_ac - (self.cm_ac / self.cl_a_min_h) + ((self.cl_h / self.cl_a_min_h) * sh_s * (self.l_h /self.chord) * (self.velocity_h / self.velocity)**2)
    
    def convert_cog_abs_to_rel_mac(self, x_cog_abs):
        #make leading edge (should technically be LE of MAC(adjust in glider.py)) the datum
        return (x_cog_abs - self.wing_x_location)/self.chord
    
    def convert_cog_abs_to_rel_mac(self, x_cog_abs):
        #make leading edge (should technically be LE of MAC(adjust in glider.py)) the datum
        return (x_cog_abs - self.wing_x_location)/self.chord
    
    def size_tail(self):
        Sh_S = 0.2
        dSh_S = -0.1
        max_it = 200
        i = 0
        fwd_limit_normalised = self.convert_cog_abs_to_rel_mac(self.fwd_limit_x_cog)
        bwd_limit_normalised = self.convert_cog_abs_to_rel_mac(self.bwd_limit_x_cog)
        while i < max_it:
            stability_limit = self.x_cg_stability_limit(Sh_S)
            controllability_limit = self.x_cg_controllability_limit(Sh_S)
            if fwd_limit_normalised < controllability_limit or bwd_limit_normalised > stability_limit:
                Sh_S += dSh_S
                break
            else:
                Sh_S -= dSh_S
                i += 1
        return Sh_S

    #Define Sh/S range
    s_h_s_range = np.linspace(0.005, 0.5, 100)

    #Plot the scissor plot
    def plot_scissor_plot(self):
        plt.figure(figsize=(10,6))
        plt.plot(self.x_cg_stability_limit(self.s_h_s_range), self.s_h_s_range, label='Stability Limit', color='blue')
        plt.plot(self.x_cg_stability_limit_no_sm(self.s_h_s_range), self.s_h_s_range, label='Stability Limit (no margin)', linestyle = '--', color='blue')
        plt.plot(self.x_cg_controllability_limit(self.s_h_s_range), self.s_h_s_range, label='Controllability Limit', color='red')
        plt.axvline(self.convert_cog_abs_to_rel_mac(self.current_x_cog), color='black', linestyle='-', label=f'Current Xcg')
        plt.axvline(self.convert_cog_abs_to_rel_mac(self.fwd_limit_x_cog), color='red', linestyle='dashdot', label=f'Forward Limit Xcg (max. pilot weight)')
        plt.axvline(self.convert_cog_abs_to_rel_mac(self.bwd_limit_x_cog), color='red', linestyle='dashdot', label=f'Backward Limit Xcg (min. pilot weight)')
        plt.axhline(self.current_sh_s, color='green', linestyle='--', label=f'Current Sh/S = {self.current_sh_s:.2f}')
        plt.xlabel('Xcg / MAC [-]')
        plt.ylabel('Sh/S [-]')
        plt.title('Scissor Plot')
        plt.grid(True)
        plt.legend()
        plt.show()