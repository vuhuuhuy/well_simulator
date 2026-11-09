# modules/glv_engine.py
import math


def calculate_tcf(temp_f: float) -> float:
  """Calculates Temperature Correction Factor for Nitrogen-charged bellows GLVs."""
  dt = max(-50.0, temp_f - 60.0)
  return 1.0 + 0.00214 * dt + 1.25e-6 * (dt**2)


def calculate_valve_pressures(
    valve_type: str,
    ptro: float,
    r_factor: float,
    pt: float,
    pc: float,
    temp_f: float = 140.0,
):
  """Calculates Valve Opening Pressure (Pvo) and Closing Pressure (Pvc) for IPO or PPO valves."""
  tcf = calculate_tcf(temp_f)
  r = max(0.0, min(0.95, r_factor))

  if valve_type.upper() == "IPO":
    # Injection-Pressure-Operated (casing sensitive)
    p_d = (ptro / max(1.0 - r, 0.05)) * tcf
    p_vo = (p_d - pt * r) / max(1.0 - r, 0.05)
    p_vc = p_d
    return round(p_vo, 1), round(p_vc, 1)
  else:
    # Production-Pressure-Operated (tubing sensitive / PPO)
    p_spring = ptro
    p_to = (p_spring - pc * r) / max(1.0 - r, 0.05)
    p_tc = p_spring
    return round(p_to, 1), round(p_tc, 1)


def thornhill_craver_rate(
    p_up_psig: float,
    p_down_psig: float,
    port_size_in: float,
    gas_sg: float = 0.65,
    temp_f: float = 140.0,
) -> float:
  """Calculates Thornhill-Craver Orifice Gas Throughput Rate (MSCF/d)."""
  p1 = p_up_psig + 14.7
  p2 = p_down_psig + 14.7

  if p1 <= 0 or port_size_in <= 0:
    return 0.0

  pr_ratio = max(0.0, min(1.0, p2 / p1))

  # Critical ratio cutoff for natural gas (~0.55)
  if pr_ratio <= 0.55:
    flow_func = 0.388
  else:
    flow_func = math.sqrt(
        max(0.0, pr_ratio ** (2 / 1.3) - pr_ratio ** ((1.3 + 1) / 1.3))
    )

  port_area = math.pi * ((port_size_in / 64.0 if port_size_in > 1.0 else port_size_in) / 2.0) ** 2
  cd = 0.85
  temp_r = max(460.0, temp_f + 460.0)

  q_mscfd = (
      155.5 * cd * port_area * p1 * flow_func / math.sqrt(max(gas_sg * temp_r, 1e-4))
  )
  return round(max(0.0, q_mscfd), 1)