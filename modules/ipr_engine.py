# modules/ipr_engine.py
from typing import Tuple
import numpy as np


def calculate_ipr(
    fluid_type: str,
    ipr_model: str,
    pr: float,
    tres_f: float = 150.0,
    pi_input: float = 5.8,
    pb: float = 1911.0,
    test_q1: float = 1000.0,
    test_pwf1: float = 1800.0,
    fet_q1: float = 800.0,
    fet_pwf1: float = 2000.0,
    fet_q2: float = 1500.0,
    fet_pwf2: float = 1400.0,
    k_md: float = 100.0,
    h_ft: float = 20.0,
    mu_cp: float = 1.2,
    bo: float = 1.15,
    re_ft: float = 1500.0,
    rw_ft: float = 0.33,
    skin: float = 0.0,
    gas_sg: float = 0.65,
    z_factor_input: float = 0.90,
    c_backpress: float = 0.005,
    n_backpress: float = 0.85,
    num_points: int = 50,
) -> Tuple[np.ndarray, np.ndarray, float, float, float]:
  """Rigorous IPR calculation matching NExT Gas Reservoir Engineering standards."""
  pwf_array = np.linspace(pr, 0.0, num_points)
  q_array = []
  model = ipr_model.upper().strip()

  calculated_a = 0.0
  calculated_b = 0.0

  if "GAS" in fluid_type.upper():
    t_rankine = tres_f + 460.0
    denom = np.log(re_ft / rw_ft) - 0.75 + skin
    if denom <= 0.01:
      denom = 0.01

    # Dynamic gas viscosity calculation (Lee et al. / NExT average)
    mu_g = (
        1e-4
        * (7.77 + 0.0083 * 28.96 * gas_sg)
        * (t_rankine**1.5)
        / (122.4 + 12.9 * 28.96 * gas_sg + t_rankine)
    )
    mu_g = max(0.01, min(0.05, mu_g))

    k_safe = max(0.1, k_md)
    h_safe = max(1.0, h_ft)
    hp_safe = h_safe  # Perforated interval default to pay zone
    z_safe = max(0.5, z_factor_input)

    if "JONES" in model:
      # NExT Slide 16: Turbulence parameter beta (ft^-1) via Firoozabadi & Katz / Cooke
      beta = 2.33e10 / (k_safe**1.201)

      # NExT Slide 11: Jones Laminar 'a' [psi^2 / (Mscf/D)]
      calculated_a = (
          1.424e3 * mu_g * z_safe * t_rankine * denom / (k_safe * h_safe)
      )

      # NExT Slide 11: Jones Turbulent 'b' [psi^2 / (Mscf/D)^2]
      calculated_b = (
          3.16e-12 * beta * gas_sg * t_rankine * z_safe / ((hp_safe**2) * rw_ft)
      )

      # NExT Slide 18: Robust Quadratic formulation: q = 2*(Pr^2 - Pwf^2) / (a + sqrt(a^2 + 4*b*(Pr^2 - Pwf^2)))
      for pwf in pwf_array:
        dp_sq = max(0.0, pr**2 - pwf**2)
        q = (2.0 * dp_sq) / (
            calculated_a + np.sqrt(max(0.0, calculated_a**2 + 4.0 * calculated_b * dp_sq))
        )
        q_array.append(max(0.0, q))

    elif "BACK PRESSURE" in model:
      for pwf in pwf_array:
        dp_sq = max(0.0, pr**2 - pwf**2)
        q = c_backpress * (dp_sq**n_backpress)
        q_array.append(max(0.0, q))

    else:
      # NExT Slide 9: Darcy Gas Law: q = 7.03e-4 * k * h * (Pr^2 - Pwf^2) / (mu_g * T * z * denom)
      for pwf in pwf_array:
        dp_sq = max(0.0, pr**2 - pwf**2)
        q = (7.03e-4 * k_safe * h_safe * dp_sq) / (mu_g * t_rankine * z_safe * denom)
        q_array.append(max(0.0, q))

  else:
    # --- OIL & WATER IPR ---
    j_effective = pi_input

    if "PI ENTRY" in model:
      j_effective = max(0.01, pi_input)

    elif "VOGEL" in model:
      q1 = max(0.1, test_q1)
      pwf1 = max(0.0, min(pr - 1.0, test_pwf1))

      if pwf1 >= pb:
        j_effective = q1 / max(0.1, pr - pwf1)
      else:
        vogel_term = 1.0 - 0.2 * (pwf1 / pb) - 0.8 * ((pwf1 / pb) ** 2)
        denom_v = (pr - pb) + (pb / 1.8) * vogel_term
        j_effective = q1 / max(0.1, denom_v)

    elif "DARCY" in model:
      denom_d = np.log(re_ft / rw_ft) - 0.75 + skin
      if denom_d <= 0.05:
        denom_d = 0.05
      j_effective = (k_md * h_ft) / (
          141.2 * max(bo, 1.0) * max(mu_cp, 0.01) * denom_d
      )

    elif "FETKOVICH" in model:
      q1, pwf1 = max(0.1, fet_q1), max(0.0, min(pr - 1.0, fet_pwf1))
      q2, pwf2 = max(0.1, fet_q2), max(0.0, min(pr - 1.0, fet_pwf2))

      p_diff1 = max(1.0, pr**2 - pwf1**2)
      p_diff2 = max(1.0, pr**2 - pwf2**2)

      if abs(p_diff1 - p_diff2) > 1e-3 and abs(q1 - q2) > 1e-3:
        n_fet = np.log(q1 / q2) / np.log(p_diff1 / p_diff2)
        n_fet = max(0.5, min(1.0, n_fet))
        c_fet = q1 / (p_diff1**n_fet)
      else:
        n_fet = 1.0
        c_fet = q1 / p_diff1

    for pwf in pwf_array:
      if "FETKOVICH" in model:
        dp_sq = max(0.0, pr**2 - pwf**2)
        q = c_fet * (dp_sq**n_fet)
      else:
        if pr <= pb:
          q_max = (j_effective * pr) / 1.8
          q = q_max * (1.0 - 0.2 * (pwf / pr) - 0.8 * ((pwf / pr) ** 2))
        else:
          if pwf >= pb:
            q = j_effective * (pr - pwf)
          else:
            q_b = j_effective * (pr - pb)
            q_v = (j_effective * pb) / 1.8
            q_delta = q_v * (1.0 - 0.2 * (pwf / pb) - 0.8 * ((pwf / pb) ** 2))
            q = q_b + q_delta

      q_array.append(max(0.0, q))

  q_liq_array = np.array(q_array)
  aof = float(q_liq_array[-1])
  return (
      q_liq_array,
      pwf_array,
      round(aof, 1),
      calculated_a,
      calculated_b,
  )