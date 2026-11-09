import numpy as np

def calculate_pressure_traverse(
    depth: float,
    thp: float,
    pwf: float,
    num_steps: int = 50
):
    """
    Generates depth steps and corresponding flowing tubing pressure profile.
    """
    depth_array = np.linspace(0, depth, num_steps)
    
    if depth <= 0:
        return depth_array, np.full(num_steps, thp)
        
    # Pressure profile along depth (surface THP -> bottomhole FBHP)
    # Using smooth curve profile to reflect hydrostatic gas/liquid mixture behavior
    p_profile = thp + (pwf - thp) * (depth_array / depth) ** 0.92
    
    return depth_array, p_profile


def calculate_choke_whp(p_sep: float, choke_64ths: float, q_liq: float, total_gor: float):
    """
    Simplified surface choke model computing WHP and FLP based on downstream separator pressure.
    """
    if q_liq <= 0 or choke_64ths <= 0:
        return p_sep + 50.0, p_sep
        
    # Choke pressure drop approximation based on flow rate and choke size
    choke_ratio = (choke_64ths / 64.0) ** 2
    delta_p_choke = max(15.0, (q_liq / 100.0) / max(choke_ratio, 0.01))
    
    flp = p_sep + 15.0  # Flowline pressure slightly above separator pressure
    whp = flp + delta_p_choke
    
    return round(whp, 1), round(flp, 1)