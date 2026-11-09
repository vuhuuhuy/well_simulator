# modules/unloading_engine.py
from typing import List, Tuple


def calculate_initial_fluid_level(
    p_res: float,
    thp: float,
    grad_u_tubing: float,
    grad_gas: float,
    packer_depth: float,
) -> float:
  """Calculates initial static tubing liquid level depth from surface prior to gas-lift unloading."""
  denominator = grad_gas - grad_u_tubing
  if abs(denominator) < 1e-5:
    return 0.0

  fluid_level = (p_res - thp - (grad_u_tubing * packer_depth)) / denominator
  return max(0.0, min(fluid_level, packer_depth))


def calculate_unloading_traverses(
    depth: float,
    pr: float,
    thp: float,
    chp: float,
    grad_tubing_kill: float,
    gas_sg: float,
    step_ft: float = 100.0,
) -> Tuple[List[float], List[float], List[float], float]:
  """Calculates casing pressure and tubing static pressure traverses down to packer depth."""
  grad_gas = 0.1 * gas_sg
  grad_casing_kill = grad_tubing_kill
  fluid_lvl = calculate_initial_fluid_level(
      pr, thp, grad_tubing_kill, grad_gas, depth
  )

  depths = []
  p_tubing = []
  p_casing = []
  current_depth = 0.0

  while current_depth <= depth:
    depths.append(current_depth)
    p_casing.append(chp + (grad_casing_kill * current_depth))

    if current_depth <= fluid_lvl:
      p_tubing.append(thp + (grad_gas * current_depth))
    else:
      p_t_at_fl = thp + (grad_gas * fluid_lvl)
      p_tubing.append(
          p_t_at_fl + (grad_tubing_kill * (current_depth - fluid_lvl))
      )
    current_depth += step_ft

  if len(depths) == 0 or depths[-1] < depth:
    depths.append(depth)
    p_casing.append(chp + (grad_casing_kill * depth))
    p_t_at_fl = thp + (grad_gas * fluid_lvl)
    p_tubing.append(p_t_at_fl + (grad_tubing_kill * (depth - fluid_lvl)))

  return depths, p_casing, p_tubing, fluid_lvl