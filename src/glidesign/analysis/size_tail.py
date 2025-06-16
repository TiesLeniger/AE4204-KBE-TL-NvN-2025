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
    
    def _convert_cog_abs_to_rel_mac(x_LEMAC: float, x_cog_abs: float, mac: float):
        #make leading edge (should technically be LE of MAC(adjust in glider.py)) the datum
        return (x_cog_abs - x_LEMAC)/mac
    x_cg_fwd_lim = _convert_cog_abs_to_rel_mac(x_LEMAC, x_cg_fwd_lim, mac)
    x_cg = _convert_cog_abs_to_rel_mac(x_LEMAC, x_cg, mac)
    x_cg_bwd_lim = _convert_cog_abs_to_rel_mac(x_LEMAC, x_cg_bwd_lim, mac)
    x_np = _convert_cog_abs_to_rel_mac(x_LEMAC, x_np, mac)

    delta_x_cg = x_cg_bwd_lim - x_cg_fwd_lim
    Sh_S = ((delta_x_cg+(x_np-x_cg))-(Cm_ac/CL_max))/(((CL_alpha_h/CL_alpha_a_min_h)-(CL_h/CL_max))* Vh_V**2 * (l_h/mac))

    return Sh_S