# modules/schematic_engine.py
import streamlit.components.v1 as components


def render_well_schematic(
    depth: float = 7000,
    tubing_id: float = 2.992,
    thp: float = 200,
    q_liq: float = 1122,
    q_gl: float = 0,
    gl_depth: float = 6000,
    pwf: float = 1957,
):
  """Renders an interactive petroleum wellbore schematic with completely open tubing end (no horizontal dashed line)."""
  is_flowing = q_liq > 0
  fluid_color = "#38bdf8" if is_flowing else "#475569"
  gas_lift_active = q_gl > 0
  gl_status_text = "GLV Active" if gas_lift_active else "SPM (Standby)"
  gl_status_color = "#22c55e" if gas_lift_active else "#eab308"
  p_res = pwf + 200.0

  # Dynamic scaling for Gas Lift Valve position relative to tubing depth (y=80 surface to y=430 packer)
  spm_y = 80 + int(330 * min(gl_depth, depth - 100) / max(depth, 1.0))

  svg_code = f"""
    <svg viewBox="0 0 490 550" width="100%" height="100%" xmlns="http://www.w3.org/2000/svg" style="background-color: #030712; border-radius: 8px; font-family: sans-serif;">
      <defs>
        <!-- High Contrast Surface / Sky Header Gradient -->
        <linearGradient id="surfaceGrad" x1="0%" y1="0%" x2="0%" y2="100%">
          <stop offset="0%" stop-color="#0f172a"/>
          <stop offset="100%" stop-color="#1e293b"/>
        </linearGradient>

        <!-- Subterranean Earth Layer Gradient -->
        <linearGradient id="earthGrad" x1="0%" y1="0%" x2="0%" y2="100%">
          <stop offset="0%" stop-color="#2a1810"/>
          <stop offset="40%" stop-color="#1c100a"/>
          <stop offset="100%" stop-color="#0f0905"/>
        </linearGradient>

        <!-- Full Unbroken Reservoir Layer Gradient -->
        <linearGradient id="reservoirGrad" x1="0%" y1="0%" x2="100%" y2="0%">
          <stop offset="0%" stop-color="#d97706" stop-opacity="0.9"/>
          <stop offset="50%" stop-color="#b45309" stop-opacity="0.95"/>
          <stop offset="100%" stop-color="#d97706" stop-opacity="0.9"/>
        </linearGradient>
      </defs>

      <!-- Top Section (Surface / Xmas Tree Header) -->
      <rect x="0" y="0" width="490" height="80" fill="url(#surfaceGrad)" />
      <line x1="0" y1="80" x2="490" y2="80" stroke="#0284c7" stroke-width="2.5" opacity="0.8"/>

      <!-- Subterranean Earth Layer -->
      <rect x="0" y="80" width="490" height="350" fill="url(#earthGrad)" />

      <!-- Full Continuous Reservoir Layer (y=430 to y=550) -->
      <rect x="0" y="430" width="490" height="120" fill="url(#reservoirGrad)" />

      <!-- Wellhead Christmas Tree -->
      <g id="xmas-tree" stroke="#e2e8f0" stroke-width="3" fill="none">
        <path d="M 120 80 L 120 20 M 95 38 L 145 38 M 95 58 L 145 58 L 190 58" />
        <circle cx="120" cy="20" r="6" fill="#ef4444"/>
        <circle cx="95" cy="38" r="5" fill="#f59e0b"/>
        <circle cx="145" cy="38" r="5" fill="#f59e0b"/>
      </g>
      <text x="195" y="63" fill="#38bdf8" font-size="12" font-weight="bold">THP = {round(thp, 0):.0f} psi</text>

      <!-- Cement Strings Extended to Bottom of Reservoir (y=535) -->
      <rect x="65" y="80" width="15" height="455" fill="#475569" opacity="0.5"/>
      <rect x="160" y="80" width="15" height="455" fill="#475569" opacity="0.5"/>

      <!-- Upper Casing String down to Packer depth (y=430) -->
      <rect x="80" y="80" width="80" height="350" fill="#111827" stroke="#64748b" stroke-width="2.5"/>

      <!-- Wellbore Cavity Below Tubing (Rat Hole x=80 to x=160, y=430 to y=535): Pure Background Fill -->
      <rect x="80" y="430" width="80" height="105" fill="#030712"/>
      <rect x="80" y="430" width="80" height="105" fill="{fluid_color}" opacity="0.3"/>

      <!-- Dashed Perforated Liner Side Walls ONLY (Left x=80 & Right x=160) - NO top stroke across y=430 -->
      <line x1="80" y1="430" x2="80" y2="535" stroke="#64748b" stroke-width="2.5" stroke-dasharray="6,4"/>
      <line x1="160" y1="430" x2="160" y2="535" stroke="#64748b" stroke-width="2.5" stroke-dasharray="6,4"/>

      <!-- Inner Tubing Fluid Column (Above Packers y=80 to y=430) -->
      <rect x="102" y="80" width="36" height="350" fill="#030712"/>
      <rect x="102" y="80" width="36" height="350" fill="{fluid_color}" opacity="0.3"/>

      <!-- Tubing Wall Lines (Terminating cleanly at packer depth y=430 with NO horizontal line) -->
      <line x1="102" y1="80" x2="102" y2="430" stroke="#38bdf8" stroke-width="2.5"/>
      <line x1="138" y1="80" x2="138" y2="430" stroke="#38bdf8" stroke-width="2.5"/>

      <!-- Bottom Rat Hole Base Plate (y=535) -->
      <line x1="80" y1="535" x2="160" y2="535" stroke="#64748b" stroke-width="3"/>

      <!-- Production Packers (Sealing Annulus at top of Reservoir y=410 to y=430) -->
      <g id="packers" fill="#ef4444" stroke="#991b1b" stroke-width="1.5">
        <rect x="80" y="410" width="22" height="20" rx="2"/>
        <line x1="80" y1="410" x2="102" y2="430" stroke="#ffffff" stroke-width="1.5"/>
        <line x1="80" y1="430" x2="102" y2="410" stroke="#ffffff" stroke-width="1.5"/>
        <rect x="138" y="410" width="22" height="20" rx="2"/>
        <line x1="138" y1="410" x2="160" y2="430" stroke="#ffffff" stroke-width="1.5"/>
        <line x1="138" y1="430" x2="160" y2="410" stroke="#ffffff" stroke-width="1.5"/>
      </g>

      <!-- Side Pocket Mandrel & GLV -->
      <g id="spm" transform="translate(130, {spm_y})">
        <rect x="0" y="0" width="10" height="20" fill="{gl_status_color}" rx="2"/>
        <text x="16" y="14" fill="{gl_status_color}" font-size="10" font-weight="bold">{gl_status_text}</text>
      </g>

      <!-- Perforations across Full Reservoir Height (y=430 to y=530) -->
      <g id="perfs" stroke="#f59e0b" stroke-width="2.5">
        <line x1="50" y1="450" x2="80" y2="450"/>
        <line x1="50" y1="475" x2="80" y2="475"/>
        <line x1="50" y1="500" x2="80" y2="500"/>
        <line x1="50" y1="525" x2="80" y2="525"/>
        <line x1="160" y1="450" x2="190" y2="450"/>
        <line x1="160" y1="475" x2="190" y2="475"/>
        <line x1="160" y1="500" x2="190" y2="500"/>
        <line x1="160" y1="525" x2="190" y2="525"/>
      </g>

      <!-- Reservoir & Well Data Overlay Box -->
      <g id="res-box" transform="translate(200, 80)">
        <rect x="0" y="0" width="275" height="165" rx="8" fill="#0f172a" stroke="#0284c7" stroke-width="1.5" opacity="0.95"/>
        <text x="15" y="24" fill="#38bdf8" font-size="12" font-weight="bold">RESERVOIR &amp; WELL DATA</text>
        <line x1="15" y1="32" x2="260" y2="32" stroke="#334155" stroke-width="1"/>
        <text x="15" y="52" fill="#94a3b8" font-size="11">Reservoir Pressure (Pr)</text>
        <text x="260" y="52" fill="#fde047" font-size="11" text-anchor="end" font-weight="bold">{p_res:.0f} psi</text>
        <text x="15" y="74" fill="#94a3b8" font-size="11">Tubing Depth</text>
        <text x="260" y="74" fill="#ffffff" font-size="11" text-anchor="end" font-weight="bold">{depth:.0f} ft</text>
        <text x="15" y="96" fill="#94a3b8" font-size="11">Gas Lift Depth</text>
        <text x="260" y="96" fill="#ffffff" font-size="11" text-anchor="end" font-weight="bold">{gl_depth:.0f} ft</text>
        <text x="15" y="118" fill="#94a3b8" font-size="11">Tubing Inner Diameter</text>
        <text x="260" y="118" fill="#ffffff" font-size="11" text-anchor="end" font-weight="bold">{tubing_id:.3f} in</text>
        <text x="15" y="140" fill="#94a3b8" font-size="11">FBHP (Flowing Pressure)</text>
        <text x="260" y="140" fill="#38bdf8" font-size="12" text-anchor="end" font-weight="bold">{pwf:.0f} psi</text>
      </g>
    </svg>
    """
  components.html(svg_code, height=560)