# modules/vlp_engine.py
import math
from typing import Tuple
import numpy as np


def calculate_hb_zfactor(p_psia: float, temp_f: float, gas_sg: float) -> float:
  """Calculates gas Z-factor using DAK EOS."""
  t_rankine = temp_f + 460.0
  tpc = 169.0 + 314.0 * gas_sg
  ppc = 708.75 - 57.5 * gas_sg

  tpr = t_rankine / max(tpc, 1.0)
  ppr = p_psia / max(ppc, 1.0)

  a = [
      -2.462,
      2.97,
      -0.2862,
      0.008054,
      2.808,
      -3.498,
      0.3603,
      -0.01044,
      -0.7933,
      1.396,
      -0.1491,
      0.00441,
      0.08393,
      -0.1864,
      0.02033,
      -0.0006095,
  ]

  t = 1.0 / tpr
  a_term = 0.0707 * ppr * t * math.exp(-0.32 * ppr)
  b_term = (a[0] + a[1] * t + a[2] * t**2) * ppr
  c_term = (a[3] + a[4] * t) * ppr**2
  d_term = a[5] * t * ppr**5
  y_term = a_term + b_term + c_term + d_term

  z = 1.0 + y_term
  return max(0.60, min(1.50, z))


def calculate_vlp(
    fluid_type: str,
    vlp_corr: str,
    q_rate_array: np.ndarray,
    depth: float = 7000.0,
    tubing_id: float = 2.992,
    thp: float = 175.0,
    wc: float = 28.0,
    gor: float = 1800.0,
    q_gl: float = 0.0,
    oil_api: float = 41.7,
    gas_sg: float = 0.65,
    water_sg: float = 1.025,
    wgr: float = 10.0,
    cgr: float = 15.0,
    temp_surf: float = 80.0,
    temp_res: float = 150.0,
    operating_q: float = None,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
  """Multiphase VLP Engine with distinct, active correlation dispatching."""
  pwf_vlp = []
  is_gas_well = "GAS" in fluid_type.upper()

  d_ft = tubing_id / 12.0
  area_ft2 = (math.pi / 4.0) * (d_ft**2)
  roughness_ft = 0.00065
  g_c = 32.174

  num_steps = 25
  dz = depth / num_steps
  corr = vlp_corr.upper().strip()

  sg_oil = 141.5 / (131.5 + max(oil_api, 10.0))
  rho_oil_std = sg_oil * 62.4
  rho_water_std = max(1.0, water_sg) * 62.4

  def compute_dp_dz(q_val, current_p, current_z):
    t_avg = temp_surf + (temp_res - temp_surf) * (current_z + dz / 2.0) / depth
    p_calc = max(current_p, 14.7)

    if is_gas_well:
      q_gas_scfd = q_val * 1000.0
      q_liq_stbd = (q_val / 1000.0) * (wgr + cgr)
      total_liq = wgr + cgr
      wc_frac = (wgr / total_liq) if total_liq > 0 else 0.0
      sigma_L = (30.0 * cgr + 60.0 * wgr) / max(total_liq, 1e-4) if total_liq > 0 else 30.0
      rho_liq_std = wc_frac * rho_water_std + (1.0 - wc_frac) * rho_oil_std
      mu_L = 0.5
    else:
      q_liq_stbd = q_val
      total_gor = gor + (q_gl * 1000.0 / max(q_liq_stbd, 1.0))
      q_gas_scfd = q_liq_stbd * total_gor
      wc_frac = wc / 100.0
      rho_liq_std = wc_frac * rho_water_std + (1.0 - wc_frac) * rho_oil_std
      sigma_L = 30.0
      mu_L = 5.0 * (1.0 - wc_frac) + 0.68 * wc_frac

    v_sl = (q_liq_stbd * 5.615 / 86400.0) / area_ft2

    z_fact = calculate_hb_zfactor(p_calc, t_avg, gas_sg)
    bg = 0.02827 * z_fact * (t_avg + 460.0) / p_calc
    v_sg = (q_gas_scfd * bg / 86400.0) / area_ft2
    v_m = v_sl + v_sg

    C_L = v_sl / max(v_m, 1e-6)
    rho_gas_insitu = (gas_sg * 0.0764) / bg if bg > 0 else 0.01
    rho_m_noslip = rho_liq_std * C_L + rho_gas_insitu * (1.0 - C_L)

    w_mass_lbs = (
        (q_liq_stbd * 5.615 * rho_liq_std) + (q_gas_scfd * 0.0764 * gas_sg)
    ) / 86400.0

    # --------------------------------------------------------------------------
    # 1. GRAY CORRELATION
    # --------------------------------------------------------------------------
    if "GRAY" in corr:
      sigma_english = sigma_L * 0.00220462
      if (wgr + cgr <= 0.0 and is_gas_well) or q_liq_stbd <= 1e-4:
        rho_mix = rho_gas_insitu
        n_re = (rho_mix * v_sg * d_ft) / max(0.018 * 6.7197e-4, 1e-7)
        f_tp = (
            (-1.8 * math.log10((roughness_ft / (3.7 * d_ft)) ** 1.11 + 6.9 / max(n_re, 1.0))) ** -2
            if n_re > 2000
            else 64.0 / max(n_re, 1.0)
        )
      else:
        N_V = (rho_m_noslip**2 * (v_m**4)) / (
            g_c * max(sigma_english, 1e-5) * max(rho_liq_std - rho_gas_insitu, 0.1)
        )
        N_D = (g_c * (rho_liq_std - rho_gas_insitu) * (d_ft**2)) / max(sigma_english, 1e-5)
        R = v_sl / max(v_sg, 1e-6)
        B = 0.0814 * (1.0 - 0.554 * math.log(1.0 + 730.0 / max(N_D, 1e-3)))
        param_A_inner = N_V * (1.0 + 205.0 / max(N_D, 1e-3))
        A = -2.2314 * (max(param_A_inner, 1e-6) ** B)

        H_g = (1.0 - math.exp(A)) / (R + 1.0)
        H_g = max(0.001, min(0.999, H_g))
        H_L = 1.0 - H_g
        rho_mix = rho_liq_std * H_L + rho_gas_insitu * H_g

        eps_prime = max(2.77e-5, (28.5 * sigma_L) / max(rho_m_noslip * (v_m**2), 1e-4))
        rel_rough = min(0.05, eps_prime / d_ft)
        n_re = (2.2e-2 * (w_mass_lbs * 86400.0)) / max(d_ft * 0.018, 1e-6)
        f_tp = (
            (-1.8 * math.log10((rel_rough / 3.7) ** 1.11 + 6.9 / max(n_re, 1.0))) ** -2
            if n_re > 2000
            else 64.0 / max(n_re, 1.0)
        )

      g_h = rho_mix / 144.0
      g_f = (f_tp * rho_mix * (v_m**2)) / (2.0 * g_c * d_ft * 144.0)

    # --------------------------------------------------------------------------
    # 2. DUNS & ROS CORRELATION (With 3 Flow Regimes & Friction Scaling)
    # --------------------------------------------------------------------------
    elif "DUNS" in corr:
      # Dimensionless velocities
      v_sg_star = v_sg * ((rho_liq_std / (g_c * sigma_L)) ** 0.25)
      v_sl_star = v_sl * ((rho_liq_std / (g_c * sigma_L)) ** 0.25)
      N_d_duns = d_ft * ((rho_liq_std * g_c / sigma_L) ** 0.5)

      # Flow regime boundaries L1, L2
      L1 = 0.01 + 0.2 * v_sl_star
      L2 = 1.0 + 0.8 * v_sl_star

      if v_sg_star < L1:
        # Regime I: Bubble Flow
        s_val = 0.98 * (1.0 + v_sl_star) / max(v_sg_star**0.2, 0.1)
        y_L = (v_sl + s_val) / max(v_m, 1e-4)
        f_factor = 1.0
      elif v_sg_star < L2:
        # Regime II: Slug / Froth Flow
        s_val = 1.57 * (1.0 + v_sl_star) / max(v_sg_star**0.1, 0.1)
        y_L = (v_sl + s_val) / max(v_m, 1e-4)
        f_factor = 1.4 + 0.3 * v_sg_star
      else:
        # Regime III: Mist / Annular Flow (Liquid entrainment)
        y_L = C_L  # Droplets move at gas velocity
        f_factor = 2.5 + 0.8 * (v_sg_star - L2)

      y_L = max(C_L, min(0.98, y_L))
      rho_mix = y_L * rho_liq_std + (1.0 - y_L) * rho_gas_insitu

      n_re = (rho_liq_std * v_sl * d_ft) / max(mu_L * 6.7197e-4, 1e-6)
      f_base = (
          (-1.8 * math.log10((roughness_ft / (3.7 * d_ft)) ** 1.11 + 6.9 / max(n_re, 1.0))) ** -2
          if n_re > 2000
          else 64.0 / max(n_re, 1.0)
      )
      f_tp = f_base * f_factor

      g_h = rho_mix / 144.0
      g_f = (f_tp * rho_mix * (v_m**2)) / (2.0 * g_c * d_ft * 144.0)

    # --------------------------------------------------------------------------
    # 3. BEGGS & BRILL CORRELATION
    # --------------------------------------------------------------------------
    elif "BEGGS" in corr:
      n_fr = (v_m**2) / max(g_c * d_ft, 1e-5)
      l1 = 316.0 * (C_L**0.302)
      l2 = 0.0009252 * (C_L**-2.4684)
      l3 = 0.10 * (C_L**-1.4516)

      if (C_L < 0.01 and n_fr < l1) or (C_L >= 0.01 and n_fr < l2):
        y_L0 = 0.98 * (C_L**0.4846) / max(n_fr**0.0868, 1e-4)
      elif (C_L >= 0.01 and l2 <= n_fr <= l3) or (C_L >= 0.001 and n_fr < l1):
        y_L0 = 0.845 * (C_L**0.5351) / max(n_fr**0.0173, 1e-4)
      else:
        y_L0 = 1.065 * (C_L**0.5824) / max(n_fr**0.0609, 1e-4)

      y_L = max(C_L, min(0.98, y_L0))
      rho_mix = y_L * rho_liq_std + (1.0 - y_L) * rho_gas_insitu

      mu_mix = C_L * mu_L + (1.0 - C_L) * 0.018
      n_re = (2.2e-2 * (w_mass_lbs * 86400.0)) / max(d_ft * mu_mix, 1e-6)
      s_val = math.log(max(C_L / max(y_L**2, 1e-4), 1e-3))
      f_ns = (
          (-1.8 * math.log10((roughness_ft / (3.7 * d_ft)) ** 1.11 + 6.9 / max(n_re, 1.0))) ** -2
          if n_re > 2000
          else 16.0 / max(n_re, 1.0)
      )
      f_tp = f_ns * math.exp(s_val)

      g_h = rho_mix / 144.0
      g_f = (f_tp * ((w_mass_lbs * 86400.0) ** 2)) / (1.067472e13 * (d_ft**5) * max(rho_mix, 1.0))

    # --------------------------------------------------------------------------
    # 4. FANCHER & BROWN CORRELATION (No-Slip Baseline Model)
    # --------------------------------------------------------------------------
    elif "FANCHER" in corr:
      y_L = C_L  # Pure No-Slip Holdup
      rho_mix = rho_m_noslip
      mu_mix = C_L * mu_L + (1.0 - C_L) * 0.018
      n_re = (2.2e-2 * (w_mass_lbs * 86400.0)) / max(d_ft * mu_mix, 1e-6)
      f_tp = 0.005 + 0.5 / max(n_re**0.32, 1.0)

      g_h = rho_mix / 144.0
      g_f = (f_tp * ((w_mass_lbs * 86400.0) ** 2)) / (1.067472e13 * (d_ft**5) * max(rho_mix, 1.0))

    # --------------------------------------------------------------------------
    # 5. HAGEDORN & BROWN CORRELATION (Default)
    # --------------------------------------------------------------------------
    else:
      N_vL = 1.938 * v_sl * ((rho_liq_std / sigma_L) ** 0.25)
      N_vg = 1.938 * v_sg * ((rho_liq_std / sigma_L) ** 0.25)
      N_D = 120.872 * d_ft * ((rho_liq_std / sigma_L) ** 0.5)
      N_L = 0.15726 * mu_L * ((1.0 / (rho_liq_std * (sigma_L**3))) ** 0.25)

      x_nl = math.log10(max(N_L, 1e-5)) + 3.0
      x_nl_clamped = max(-2.5, min(2.5, x_nl))
      cnl_coeffs = [-2.69851, 0.15840954, -0.55099756, 0.54784917, -0.12194578]
      log_cnl = sum(c * (x_nl_clamped**i) for i, c in enumerate(cnl_coeffs))
      cnl = 10 ** max(-4.0, min(2.0, log_cnl))

      param_hl = (
          (N_vL / (max(N_vg, 1e-4) ** 0.575))
          * ((p_calc / 14.7) ** 0.1)
          * (cnl / max(N_D, 1e-3))
      )
      hl_coeffs = [-0.10306578, 0.617774, -0.632946, 0.29598, -0.0401]
      x_hl = math.log10(max(param_hl, 1e-4))
      x_hl_clamped = max(-2.0, min(1.5, x_hl))

      if param_hl < 0.01:
        ratio_hl = max(0.0, min(1.0, param_hl))
      else:
        log_ratio_hl = sum(c * (x_hl_clamped**i) for i, c in enumerate(hl_coeffs))
        ratio_hl = 10 ** max(-2.0, min(0.0, log_ratio_hl))

      param_psi = (N_vg * (N_L**0.38)) / (max(N_D, 1e-3) ** 2.14)
      psi_coeffs = [
          0.91162574,
          -4.82175636,
          1232.25036621,
          -22253.57617,
          116174.28125,
      ]
      if param_psi > 0.01:
        x_psi = max(-2.5, min(1.5, math.log10(param_psi)))
        log_psi = sum(c * (x_psi**i) for i, c in enumerate(psi_coeffs))
        psi = 10 ** max(-1.0, min(0.0, log_psi))
      else:
        psi = 1.0

      y_L = max(C_L, min(1.0, ratio_hl * psi))
      rho_mix = y_L * rho_liq_std + (1.0 - y_L) * rho_gas_insitu

      mu_gas = 0.01 + 0.0002 * (p_calc / 100.0)
      mu_two_phase = (mu_L**y_L) * (mu_gas ** (1.0 - y_L))
      n_re = (2.2e-2 * (w_mass_lbs * 86400.0)) / max(d_ft * mu_two_phase, 1e-6)

      f_tp = (
          (-1.8 * math.log10((roughness_ft / (3.7 * d_ft)) ** 1.11 + 6.9 / max(n_re, 1.0))) ** -2
          if n_re > 2000
          else 16.0 / max(n_re, 1.0)
      )

      g_h = rho_mix / 144.0
      g_f = (f_tp * ((w_mass_lbs * 86400.0) ** 2)) / (1.067472e13 * (d_ft**5) * max(rho_mix, 1.0))

    return g_h + g_f

  for q_val in q_rate_array:
    if q_val <= 0.1:
      grad_static = 0.05 if is_gas_well else (rho_oil_std / 144.0)
      pwf_vlp.append(thp + grad_static * depth)
      continue

    current_p = thp
    current_z = 0.0
    for _ in range(num_steps):
      dp_dz = compute_dp_dz(q_val, current_p, current_z)
      current_p += dp_dz * dz
      current_z += dz
    pwf_vlp.append(max(0.0, current_p))

  q_eval = (
      operating_q
      if (operating_q is not None and operating_q > 0)
      else float(q_rate_array[len(q_rate_array) // 2])
  )
  depth_profile = np.linspace(0, depth, num_steps + 1)
  temp_profile = temp_surf + (temp_res - temp_surf) * (depth_profile / depth)

  pressure_profile = [thp]
  curr_p = thp
  curr_z = 0.0
  for _ in range(num_steps):
    dp_dz = compute_dp_dz(q_eval, curr_p, curr_z)
    curr_p += dp_dz * dz
    curr_z += dz
    pressure_profile.append(curr_p)

  return (
      np.array(pwf_vlp),
      depth_profile,
      np.array(pressure_profile),
      np.array(temp_profile),
  )


def find_operating_point(
    q_liq: np.ndarray, pwf_ipr: np.ndarray, pwf_vlp: np.ndarray
):
  """Finds intersection point between IPR and VLP curves."""
  diff = pwf_ipr - pwf_vlp
  idx = np.where(np.diff(np.sign(diff)))[0]

  if len(idx) == 0:
    return None, None

  i = idx[-1]
  q1, q2 = q_liq[i], q_liq[i + 1]
  d1, d2 = diff[i], diff[i + 1]

  if (d1 - d2) == 0:
    return round(float(q1), 1), round(float(pwf_ipr[i]), 1)

  q_op = q1 + (0 - d1) * (q2 - q1) / (d2 - d1)
  pwf_op = np.interp(q_op, q_liq, pwf_ipr)
  return round(float(q_op), 1), round(float(pwf_op), 1)