def size_tail(
        x_cg_fwd_lim: float, 
        x_cg: float, 
        x_cg_bwd_lim: float, 
        mac: float,
        x_LEMAC: float, 
        x_np: float, 
        Cm_ac: float, 
        CL_max: float,
        CL_alpha_h: float,
        CL_alpha_a_min_h: float,
        CL_h: float,
        Vh_V: float,
        l_h: float
        ) -> float:
    """
    Calculates the required horizontal tail surface area ratio (Sh/S) to ensure 
    stability and controllability of the aircraft based on its center of gravity (CG) 
    position, aerodynamic parameters, and design constraints.

    Args:
        x_cg_fwd_lim (float): The forward limit of the center of gravity (CG) location 
                               relative to the leading edge of the mean aerodynamic chord (MAC).
        x_cg (float): The current center of gravity (CG) location relative to the leading 
                      edge of the mean aerodynamic chord (MAC).
        x_cg_bwd_lim (float): The backward limit of the CG location relative to the 
                               leading edge of the MAC.
        mac (float): The mean aerodynamic chord length of the wing (in meters).
        x_LEMAC (float): The location of the leading edge of the mean aerodynamic chord (MAC).
        x_np (float): The neutral point (the point where the pitching moment coefficient is zero).
        Cm_ac (float): The pitching moment coefficient around the aerodynamic center.
        CL_max (float): The maximum lift coefficient of the aircraft.
        CL_alpha_h (float): The lift curve slope coefficient of the horizontal tail.
        CL_alpha_a_min_h (float): The lift curve slope coefficient of the airplane minus the horizontal tail.
        CL_h (float): The lift coefficient of the horizontal tail.
        Vh_V (float): The ratio of velocity at the horizontal tail to velocity at the main wing.
        l_h (float): The distance from the aerodynamic center of the main wing to the aerodynamic center of the horizontal tail.

    Returns:
        float: The required horizontal tail surface area ratio (Sh/S) for stability and controllability.

    Example:
        Sh_S = size_tail(x_cg_fwd_lim=0.1, x_cg=0.2, x_cg_bwd_lim=0.3, mac=2.0, x_LEMAC=0.0, 
                          x_np=0.4, Cm_ac=0.05, CL_max=1.5, CL_alpha_h=0.1, CL_alpha_a_min_h=0.2, 
                          CL_h=0.5, Vh_V=1.0, l_h=5.0)
        print(Sh_S)  # Output: the calculated horizontal tail surface area ratio (Sh/S)
    """
    
    def _convert_cog_abs_to_rel_mac(x_LEMAC: float, x_cog_abs: float, mac: float) -> float:
        """
        Converts the absolute center of gravity (CG) location to a relative position 
        with respect to the mean aerodynamic chord (MAC).

        Args:
            x_LEMAC (float): The location of the leading edge of the mean aerodynamic chord (MAC).
            x_cog_abs (float): The absolute position of the center of gravity (CG).
            mac (float): The mean aerodynamic chord length.

        Returns:
            float: The relative center of gravity (CG) position with respect to MAC.
        """
        # Make the leading edge of the MAC the datum point and return relative CG position
        return (x_cog_abs - x_LEMAC) / mac

    # Convert absolute CG locations to relative positions with respect to the MAC
    x_cg_fwd_lim = _convert_cog_abs_to_rel_mac(x_LEMAC, x_cg_fwd_lim, mac)
    x_cg = _convert_cog_abs_to_rel_mac(x_LEMAC, x_cg, mac)
    x_cg_bwd_lim = _convert_cog_abs_to_rel_mac(x_LEMAC, x_cg_bwd_lim, mac)
    x_np = _convert_cog_abs_to_rel_mac(x_LEMAC, x_np, mac)

    # Calculate the difference between the forward and backward CG limits
    delta_x_cg = x_cg_bwd_lim - x_cg_fwd_lim

    # Calculate the required horizontal tail surface area ratio (Sh/S) using the provided formula
    Sh_S = ((delta_x_cg + (x_np - x_cg)) - (Cm_ac / CL_max)) / (
        ((CL_alpha_h / CL_alpha_a_min_h) - (CL_h / CL_max)) * Vh_V ** 2 * (l_h / mac)
    )

    return Sh_S