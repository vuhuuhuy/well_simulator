import numpy as np
import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import streamlit as st

# Local modular imports
from modules.glv_engine import calculate_valve_pressures
from modules.ipr_engine import calculate_ipr
from modules.schematic_engine import render_well_schematic
from modules.unloading_engine import calculate_unloading_traverses
from modules.vlp_engine import calculate_vlp, find_operating_point

# --- STREAMLIT PAGE CONFIGURATION ---
st.set_page_config(
    page_title="Well Simulator — Petroleum Engineering",
    layout="wide",
    initial_sidebar_state="expanded",
)

st.markdown(
    """
    <style>
        .block-container {
            padding-top: 1rem !important;
            padding-bottom: 1rem !important;
        }
        h1, h2, h3, .stHeader {
            font-size: 1.5rem !important;
            font-weight: 700 !important;
            margin-top: 0.2rem !important;
            margin-bottom: 0.5rem !important;
        }
    </style>
""",
    unsafe_allow_html=True,
)

st.title("⚡ Interactive Well Simulator & Gas-Lift Tool")

tab_live, tab_unloading = st.tabs(
    ["📊 Live Well Simulator", "⛽ Gas-Lift Unloading Process"]
)

# ==============================================================================
# SIDEBAR: MODEL CONFIGURATION & INPUT PARAMETERS
# ==============================================================================
with st.sidebar:
  st.header("⚙️ Model Configuration")

  fluid_type = st.radio(
      "Fluid Description", options=["Oil & Water", "Dry & Wet Gas"], index=0
  )

  if fluid_type == "Oil & Water":
    lift_method = st.radio(
        "Artificial Lift Method",
        options=["None (Free Flowing)", "Gas Lift", "ESP (Placeholder)"],
        index=0,
    )
  else:
    lift_method = "None (Free Flowing)"

  st.subheader("1. Fluid PVT & Reservoir")
  pr = st.slider("Reservoir Pressure (Pr) [psi]", 500, 5000, 2031, 50)
  tres = st.slider("Reservoir Temp (°F)", 100, 300, 150, 5)

  if fluid_type == "Oil & Water":
    pb = st.slider("Bubble Point Pressure (Pb) [psi]", 100, 5000, 1911, 50)
    gor = st.slider("Total GOR (scf/STB)", 0, 5000, 1800, 50)
    oil_api = st.slider("Oil Gravity (°API)", 15.0, 50.0, 41.7, 0.5)
    gas_sg = st.slider("Gas Specific Gravity", 0.55, 0.95, 0.65, 0.01)
    water_sg = st.slider(
        "Water Specific Gravity (SG)", 1.00, 1.25, 1.025, 0.005
    )
    wc = st.slider("Water Cut (%)", 0.0, 100.0, 28.0, 1.0)
    wgr, cgr = 0.0, 0.0
    z_factor_input = 0.85

    ipr_options = ["PI Entry", "Vogel", "Darcy", "Fetkovich"]
    vlp_options = [
        "Hagedorn & Brown",
        "Beggs & Brill",
        "Duns and Ros",
        "Fancher Brown (No-Slip)",
    ]
  else:
    pb, gor, oil_api, wc, water_sg = 0.0, 0.0, 0.0, 0.0, 1.0
    gas_sg = st.slider("Gas Specific Gravity", 0.55, 0.95, 0.65, 0.01)
    z_factor_input = st.slider(
        "Gas Compressibility (Z-factor)", 0.70, 1.20, 0.90, 0.01
    )
    wgr = st.number_input("WGR (bbl/MMSCF)", 0.0, 500.0, 10.0, 5.0)
    cgr = st.number_input("CGR (bbl/MMSCF)", 0.0, 200.0, 15.0, 2.5)

    ipr_options = ["Darcy", "Jones", "Back Pressure (C & n)"]
    vlp_options = ["Gray", "Hagedorn & Brown"]

  # --- 2. IPR MODEL DYNAMIC INPUTS ---
  st.subheader("2. IPR Model")
  ipr_model = st.selectbox("IPR Model Selection", options=ipr_options, index=0)

  pi_val = 5.8
  test_q1, test_pwf1 = 1000.0, 1800.0
  fet_q1, fet_pwf1, fet_q2, fet_pwf2 = 800.0, 2000.0, 1500.0, 1400.0
  k_md, h_ft, mu_cp, bo, re_ft, rw_ft, skin = (
      100.0,
      20.0,
      1.2,
      1.15,
      1500.0,
      0.33,
      0.0,
  )
  c_back, n_back = 0.005, 0.85

  if "PI ENTRY" in ipr_model.upper():
    pi_val = st.slider("Productivity Index (PI) [STB/d/psi]", 0.1, 20.0, 5.8, 0.1)

  elif "VOGEL" in ipr_model.upper():
    st.markdown("**Well Test Input Data:**")
    test_q1 = st.number_input(
        "Test Liquid Rate (q1) [STB/d]", 10.0, 10000.0, 1000.0, 50.0
    )
    test_pwf1 = st.number_input(
        "Test FBHP (Pwf1) [psi]",
        0.0,
        float(pr),
        float(min(1800, pr - 50)),
        50.0,
    )

  elif "DARCY" in ipr_model.upper() or "JONES" in ipr_model.upper():
    st.markdown("**Reservoir & Fluid Parameters:**")
    k_md = st.number_input("Permeability (k) [mD]", 1.0, 5000.0, 100.0, 10.0)
    h_ft = st.number_input("Pay Thickness (h) [ft]", 5.0, 500.0, 20.0, 5.0)
    if fluid_type == "Oil & Water":
      bo = st.number_input(
          "Oil Formation Volume Factor (Bo) [rb/STB]", 1.0, 2.0, 1.15, 0.05
      )
      mu_cp = st.number_input("Oil Viscosity (μo) [cp]", 0.1, 50.0, 1.2, 0.1)
    re_ft = st.number_input(
        "Drainage Radius (re) [ft]", 100.0, 5000.0, 1500.0, 100.0
    )
    rw_ft = st.number_input("Wellbore Radius (rw) [ft]", 0.1, 2.0, 0.33, 0.01)
    skin = st.slider("Total Skin Factor (S)", -5.0, 50.0, 0.0, 0.5)

  elif "FETKOVICH" in ipr_model.upper():
    st.markdown("**Two-Point Oil Well Test Data:**")
    c_t1, c_t2 = st.columns(2)
    fet_q1 = c_t1.number_input("Test 1 Rate (q1)", 10.0, 10000.0, 800.0, 50.0)
    fet_pwf1 = c_t1.number_input(
        "Test 1 Pwf1", 0.0, float(pr), float(min(2000, pr - 50)), 50.0
    )
    fet_q2 = c_t2.number_input("Test 2 Rate (q2)", 10.0, 10000.0, 1500.0, 50.0)
    fet_pwf2 = c_t2.number_input(
        "Test 2 Pwf2", 0.0, float(pr), float(min(1400, pr - 100)), 50.0
    )

  elif "BACK PRESSURE" in ipr_model.upper():
    c_back = st.number_input("Back Pressure C", 0.0001, 1.0, 0.005)
    n_back = st.slider("Back Pressure n", 0.5, 1.0, 0.85, 0.05)

  st.subheader("3. VLP Correlation & Tubing")
  vlp_corr = st.selectbox("VLP Correlation", options=vlp_options, index=0)
  depth = st.slider("Tubing Depth (ft)", 1000, 15000, 7000, 100)
  tubing_id = st.selectbox(
      "Tubing ID (in)",
      options=[1.995, 2.441, 2.992, 3.958, 4.892, 6.181],
      index=3,
  )
  thp = st.slider("Tubing Head Pressure (THP) [psi]", 15, 1000, 175, 5)

  q_gl = 0.0
  gl_depth = max(500, depth - 1000)
  if lift_method == "Gas Lift":
    st.subheader("4. Gas Lift Control")
    gl_depth = st.slider(
        "Gas Lift Injection Depth (ft)",
        500,
        depth,
        max(500, depth - 1000),
        100,
    )
    q_gl = st.slider("Gas Lift Injection Rate (MSCF/d)", 0, 5000, 1500, 50)

# ==============================================================================
# COMPUTATIONAL ENGINE DISPATCH
# ==============================================================================
q_liq, pwf_ipr, aof, jones_a, jones_b = calculate_ipr(
    fluid_type=fluid_type,
    ipr_model=ipr_model,
    pr=pr,
    tres_f=tres,
    pi_input=pi_val,
    pb=pb,
    test_q1=test_q1,
    test_pwf1=test_pwf1,
    fet_q1=fet_q1,
    fet_pwf1=fet_pwf1,
    fet_q2=fet_q2,
    fet_pwf2=fet_pwf2,
    k_md=k_md,
    h_ft=h_ft,
    mu_cp=mu_cp,
    bo=bo,
    re_ft=re_ft,
    rw_ft=rw_ft,
    skin=skin,
    gas_sg=gas_sg,
    z_factor_input=z_factor_input,
    c_backpress=c_back,
    n_backpress=n_back,
)

# Preliminary temporary intersection check for P-T traverse evaluation rate
pwf_vlp_temp, _, _, _ = calculate_vlp(
    fluid_type=fluid_type,
    vlp_corr=vlp_corr,
    q_rate_array=q_liq,
    depth=depth,
    tubing_id=tubing_id,
    thp=thp,
    wc=wc,
    gor=gor,
    q_gl=q_gl,
    oil_api=oil_api,
    gas_sg=gas_sg,
    water_sg=water_sg,
    wgr=wgr,
    cgr=cgr,
    temp_surf=80.0,
    temp_res=tres,
)
q_op_temp, _ = find_operating_point(q_liq, pwf_ipr, pwf_vlp_temp)

# Full evaluation returning VLP FBHP curve AND P-T profiles at Operating Point
pwf_vlp, depth_prof, press_prof, temp_prof = calculate_vlp(
    fluid_type=fluid_type,
    vlp_corr=vlp_corr,
    q_rate_array=q_liq,
    depth=depth,
    tubing_id=tubing_id,
    thp=thp,
    wc=wc,
    gor=gor,
    q_gl=q_gl,
    oil_api=oil_api,
    gas_sg=gas_sg,
    water_sg=water_sg,
    wgr=wgr,
    cgr=cgr,
    temp_surf=80.0,
    temp_res=tres,
    operating_q=q_op_temp,
)

q_op, pwf_op = find_operating_point(q_liq, pwf_ipr, pwf_vlp)

# ==============================================================================
# TAB 1: LIVE WELL SIMULATOR & NODAL ANALYSIS
# ==============================================================================
with tab_live:
  c1, c2, c3, c4 = st.columns(4)
  if q_op is not None and q_op > 0:
    if fluid_type == "Oil & Water":
      q_oil = q_op * (1.0 - wc / 100.0)
      q_water = q_op * (wc / 100.0)
      q_gas_total = (q_op * gor / 1000.0) + q_gl
      c1.metric("Liquid Rate", f"{q_op:,.0f} STB/d")
      c2.metric("Oil Rate", f"{q_oil:,.1f} STB/d")
      c3.metric("Gas Rate", f"{q_gas_total:,.1f} MSCF/d")
      c4.metric("Water Rate", f"{q_water:,.1f} STB/d")
    else:
      q_gas_total = q_op
      q_cond = (q_op / 1000.0) * cgr
      q_water = (q_op / 1000.0) * wgr
      c1.metric("Gas Rate", f"{q_gas_total:,.0f} MSCF/d")
      c2.metric("Condensate Rate", f"{q_cond:,.1f} STB/d")
      c3.metric("Water Rate", f"{q_water:,.1f} STB/d")
      c4.metric("AOF Potential", f"{aof:,.0f} MSCF/d")
  else:
    c1.metric("Production Rate", "NO FLOW")
    c2.metric("Liquid / Condensate", "0 STB/d")
    c3.metric("Total Gas Rate", "0 MSCF/d")
    c4.metric("Water Rate", "0 STB/d")

  st.divider()

  # View options toggle
  show_pt_traverse = st.checkbox(
      "🌡️ Show Pressure & Temperature Traverse Profile Plot", value=False
  )

  col_sch, col_nodal = st.columns([0.8, 1.4])

  with col_sch:
    render_well_schematic(
        depth=depth,
        tubing_id=tubing_id,
        thp=thp,
        q_liq=q_op if q_op else 0,
        q_gl=q_gl,
        gl_depth=gl_depth,
        pwf=pwf_op if q_op else 0,
    )

  with col_nodal:
    if "JONES" in ipr_model.upper() and fluid_type == "Dry & Wet Gas":
      st.info(
          f"📐 **NExT Jones Coefficients:**  \n"
          f"**a (Laminar):** `{jones_a:.6f}` psi²/(Mscf/D) | "
          f"**b (Turbulent):** `{jones_b:.10f}` psi²/(Mscf/D)²"
      )

    # Dynamic render key includes vlp_corr to force immediate Plotly canvas redrawing
    dynamic_chart_key = f"nodal_{fluid_type}_{ipr_model}_{vlp_corr}_{tubing_id}_{pr}_{thp}_{k_md}_{skin}_{z_factor_input}_{wc}_{gor}_{q_gl}"

    if not show_pt_traverse:
      # STANDARD NODAL ANALYSIS PLOT
      fig_nodal = go.Figure()
      rate_label = (
          "Liquid Rate (STB/d)"
          if fluid_type == "Oil & Water"
          else "Gas Rate (MSCF/d)"
      )

      fig_nodal.add_trace(
          go.Scatter(
              x=q_liq,
              y=pwf_ipr,
              mode="lines",
              name=f"IPR ({ipr_model})",
              line=dict(color="#0066CC", width=3),
          )
      )
      fig_nodal.add_trace(
          go.Scatter(
              x=q_liq,
              y=pwf_vlp,
              mode="lines",
              name=f"VLP ({vlp_corr})",
              line=dict(color="#FF4B4B", width=3, dash="dash"),
          )
      )

      if q_op is not None and q_op > 0:
        fig_nodal.add_trace(
            go.Scatter(
                x=[q_op],
                y=[pwf_op],
                mode="markers",
                name="Operating Point",
                marker=dict(color="#00FFCC", size=12),
            )
        )

      fig_nodal.add_trace(
          go.Scatter(
              x=[aof],
              y=[0],
              mode="markers",
              name="AOF Potential",
              marker=dict(color="#f59e0b", size=12, symbol="diamond"),
          )
      )
      fig_nodal.add_annotation(
          x=aof,
          y=0,
          text=f"<b>AOF = {aof:,.0f}</b>",
          showarrow=True,
          arrowhead=2,
          ax=0,
          ay=-35,
          bgcolor="#1e293b",
          bordercolor="#f59e0b",
          font=dict(color="#f59e0b", size=11),
      )

      fig_nodal.update_layout(
          title=dict(
              text=f"<b>Nodal Analysis ({fluid_type}) — {vlp_corr}</b>",
              y=0.98,
              x=0.0,
              xanchor="left",
              yanchor="top",
          ),
          legend=dict(
              orientation="h",
              yanchor="bottom",
              y=1.02,
              xanchor="center",
              x=0.5,
              bgcolor="rgba(0,0,0,0)",
              font=dict(size=12),
          ),
          xaxis_title=rate_label,
          yaxis_title="Flowing Bottomhole Pressure FBHP (psi)",
          xaxis=dict(range=[0, max(aof * 1.08, max(q_liq) * 1.05)]),
          yaxis=dict(range=[0, max(pr, max(pwf_vlp)) * 1.1]),
          template="plotly_dark",
          height=560,
          margin=dict(t=60, b=40, l=60, r=20),
      )
      st.plotly_chart(
          fig_nodal, use_container_width=True, key=dynamic_chart_key
      )

    else:
      # PRESSURE & TEMPERATURE TRAVERSE SUBPLOT (2 AXES)
      fig_pt = make_subplots(specs=[[{"secondary_y": True}]])

      fig_pt.add_trace(
          go.Scatter(
              x=press_prof,
              y=depth_prof,
              mode="lines+markers",
              name="Pressure Profile (psi)",
              line=dict(color="#00FFCC", width=3),
          ),
          secondary_y=False,
      )

      fig_pt.add_trace(
          go.Scatter(
              x=temp_prof,
              y=depth_prof,
              mode="lines+markers",
              name="Temperature Profile (°F)",
              line=dict(color="#FF9900", width=3, dash="dash"),
          ),
          secondary_y=True,
      )

      fig_pt.update_layout(
          title=dict(
              text=(
                  f"<b>P-T Traverse along Depth ({vlp_corr}, q = {q_op if q_op else 0:,.0f} units)</b>"
              ),
              y=0.98,
              x=0.0,
              xanchor="left",
              yanchor="top",
          ),
          legend=dict(
              orientation="h",
              yanchor="bottom",
              y=1.02,
              xanchor="center",
              x=0.5,
              bgcolor="rgba(0,0,0,0)",
              font=dict(size=12),
          ),
          xaxis_title="Pressure (psi)",
          yaxis_title="Depth (ft)",
          yaxis=dict(range=[depth * 1.05, -50], autorange=False),
          template="plotly_dark",
          height=560,
          margin=dict(t=60, b=40, l=60, r=20),
      )
      fig_pt.update_yaxes(
          title_text="Temperature (°F)", secondary_y=True, autorange=False
      )
      st.plotly_chart(
          fig_pt, use_container_width=True, key=f"pt_{dynamic_chart_key}"
      )

# ==============================================================================
# TAB 2: GAS-LIFT UNLOADING PROCESS
# ==============================================================================
with tab_unloading:
  st.subheader(
      "⛽ Gas-Lift Unloading Diagnostics (Existing Well with Installed SPMs)"
  )

  u1, u2, u3, u4 = st.columns(4)
  chp_max = u1.number_input("Max Operating CHP (psi)", 500, 3000, 1000, 50)
  kill_grad = u2.number_input(
      "Kill Fluid Gradient (psi/ft)", 0.30, 0.60, 0.40, 0.01
  )
  gas_sg_unl = u3.number_input("Lift Gas Density (SG)", 0.55, 0.90, 0.65, 0.01)
  unloading_thp = u4.number_input(
      "Unloading Surface THP (psi)", 0, 500, 150, 10
  )

  d_steps, casing_p, tubing_p, fluid_lvl = calculate_unloading_traverses(
      depth=depth,
      pr=pr,
      thp=unloading_thp,
      chp=chp_max,
      grad_tubing_kill=kill_grad,
      gas_sg=gas_sg_unl,
  )

  st.divider()

  col_tbl, col_chart = st.columns([1.1, 1.3])

  with col_tbl:
    st.markdown("### 📋 Installed SPM Schedule")
    default_spm = pd.DataFrame([
        {
            "Valve": "V1",
            "Depth (ft)": 2386,
            "Type": "IPO",
            "Port (64th)": 16,
            "R-Factor": 0.08,
            "Ptro (psi)": 1149,
        },
        {
            "Valve": "V2",
            "Depth (ft)": 3889,
            "Type": "IPO",
            "Port (64th)": 16,
            "R-Factor": 0.08,
            "Ptro (psi)": 1105,
        },
        {
            "Valve": "V3",
            "Depth (ft)": 4983,
            "Type": "IPO",
            "Port (64th)": 16,
            "R-Factor": 0.08,
            "Ptro (psi)": 1068,
        },
        {
            "Valve": "V4",
            "Depth (ft)": 5759,
            "Type": "IPO",
            "Port (64th)": 16,
            "R-Factor": 0.08,
            "Ptro (psi)": 1035,
        },
        {
            "Valve": "V5 (Orifice)",
            "Depth (ft)": 6291,
            "Type": "IPO",
            "Port (64th)": 20,
            "R-Factor": 0.00,
            "Ptro (psi)": 980,
        },
    ])

    edited_spm = st.data_editor(default_spm, use_container_width=True)

    valve_results = []
    can_kickoff = False
    for idx, row in edited_spm.iterrows():
      z = float(row["Depth (ft)"])
      pt_at_z = float(np.interp(z, d_steps, tubing_p))
      pc_at_z = float(np.interp(z, d_steps, casing_p))

      pvo, pvc = calculate_valve_pressures(
          row["Type"],
          float(row["Ptro (psi)"]),
          float(row["R-Factor"]),
          pt_at_z,
          pc_at_z,
      )
      if idx == 0:
        can_kickoff = pc_at_z >= pvo
        status = "OPEN (Kick-off)" if can_kickoff else "CLOSED"
      else:
        status = "Standby"

      valve_results.append({
          "Valve": row["Valve"],
          "Depth (ft)": z,
          "Pt (psi)": round(pt_at_z, 1),
          "Pc (psi)": round(pc_at_z, 1),
          "Pvo (psi)": pvo,
          "Pvc (psi)": pvc,
          "Status": status,
      })

    res_df = pd.DataFrame(valve_results)
    if can_kickoff:
      st.success(
          f"✅ **WELL KICK-OFF FEASIBLE:** Surface CHP of {chp_max} psi is"
          f" sufficient to open V1 ({res_df.loc[0, 'Pvo (psi)']} psi required)."
      )
    else:
      st.error(
          f"❌ **KICK-OFF BLOCKED:** Surface CHP ({chp_max} psi) at V1"
          f" ({res_df.loc[0, 'Pc (psi)']:.0f} psi) is less than required Pvo"
          f" ({res_df.loc[0, 'Pvo (psi)']} psi)."
      )

    st.markdown("#### Calculated Opening Pressures")
    st.dataframe(res_df, use_container_width=True)

  with col_chart:
    fig_unl = go.Figure()
    fig_unl.add_trace(
        go.Scatter(
            x=casing_p,
            y=d_steps,
            mode="lines",
            name="Casing Pressure (Pc)",
            line=dict(color="#FF4B4B", width=3),
        )
    )
    fig_unl.add_trace(
        go.Scatter(
            x=tubing_p,
            y=d_steps,
            mode="lines",
            name="Tubing Static (Pt)",
            line=dict(color="#00FFCC", width=3),
        )
    )
    fig_unl.add_hline(
        y=fluid_lvl,
        line_dash="dash",
        line_color="#f59e0b",
        annotation_text=f"Fluid Level = {fluid_lvl:.0f} ft",
    )

    for _, r in res_df.iterrows():
      v_color = (
          "#00FF00"
          if "OPEN" in r["Status"]
          else ("#FFCC00" if "Standby" in r["Status"] else "#FF0000")
      )
      fig_unl.add_trace(
          go.Scatter(
              x=[r["Pvo (psi)"]],
              y=[r["Depth (ft)"]],
              mode="markers+text",
              name=r["Valve"],
              marker=dict(color=v_color, size=11, symbol="diamond"),
              text=[f" {r['Valve']}: Pvo {r['Pvo (psi)']:.0f} psi"],
              textposition="top right",
              showlegend=False,
          )
      )

    fig_unl.update_layout(
        title=dict(
            text=(
                "<b>Gas-Lift Unloading Pressure Traverse & Kick-off Status</b>"
            ),
            y=0.98,
            x=0.0,
            xanchor="left",
            yanchor="top",
        ),
        legend=dict(
            orientation="h",
            yanchor="bottom",
            y=1.02,
            xanchor="center",
            x=0.5,
            bgcolor="rgba(0,0,0,0)",
            font=dict(size=12),
        ),
        xaxis_title="Pressure (psi)",
        yaxis_title="Depth (ft)",
        yaxis=dict(range=[depth * 1.08, -100], autorange=False),
        template="plotly_dark",
        height=560,
        margin=dict(t=60, b=40, l=60, r=20),
    )
    st.plotly_chart(fig_unl, use_container_width=True, key="fig_unloading_chart")