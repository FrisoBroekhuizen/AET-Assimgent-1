import numpy as np
import math 
## ------------------------------------------ Flight condition--------------------------------
Mach = 0.78
Altitude = 11000        # m
T_ambient = 216.5       # K
P_ambient = 22632       # Pa
R_air = 287             # J/kg·K
# Mass flows
mdot_air = 220          # kg/s
Bypass_ratio = 12.5
# Pressure ratios
IPR_inlet = 0.99
FPR_fan = 1.35
PR_LPC = 4.5
PR_HPC = 5.5
PR_combustor = 0.96
# Efficiencies
eta_fan = 0.92
eta_LPC = 0.92
eta_HPC = 0.92
eta_LPT = 0.91
eta_HPT = 0.91
eta_mech = 0.995
eta_combustor = 0.995
eta_nozzle = 0.99
eta_gearbox = 0.995
# Turbine inlet temperature
Tt4 = 1600              # K
# Fuel properties
LHV_fuel = 43e6         # J/kg  (43 MJ/kg)
# Specific heats and specific heat ratios
cp_air = 1000           # J/kg·K
gamma_air = 1.4
cp_gas = 1150           # J/kg·K
gamma_gas = 1.33


##  ----------------------- Equations Engine Cycle -----------------------------











## ------------------- Equations Gas Generator ------------------------------ 










## -------------------------------- Efficiency formulas --------------------------
def ETA_comb(mdot, cp, DeltaT_comb, mdot_f, LHV_f):
    """
    Combustor efficiency:
    """
    return mdot * cp * DeltaT_comb / (mdot_f * LHV_f)
def ETA_thdy(P_gg, mdot, cp, DeltaT_comb):
    """
    Thermodynamic efficiency of gas generator:
    """
    return P_gg / (mdot * cp * DeltaT_comb)
def ETA_jet_gener(mdot, vj, v0, P_gg):
    """
    Jet-generation efficiency:
    """
    mdot = np.asarray(mdot)
    vj   = np.asarray(vj)
    num = 0.5 * np.sum(mdot * (vj**2 - v0**2))
    return num / P_gg
def ETA_prop(mdot, vj, v0):
    """
    Propulsive efficiency:
    """
    mdot = np.asarray(mdot)
    vj   = np.asarray(vj)
    num = np.sum(mdot * (vj - v0) * v0)
    den = 0.5 * np.sum(mdot * (vj**2 - v0**2))
    return num / den
def ETA_thermal(mdot, vj, v0, mdot_f, LHV_f):
    """
    Overall thermal efficiency:
    UNCHOKED
    """
    mdot = np.asarray(mdot)
    vj   = np.asarray(vj)
    num = 0.5 * np.sum(mdot * (vj**2 - v0**2))
    den = mdot_f * LHV_f
    return num / den
def ETA_total(mdot, vj, v0, mdot_f, LHV_f):
    """
    Total efficiency:
    UNCHOKED
    """
    mdot = np.asarray(mdot)
    vj   = np.asarray(vj)
    num = np.sum(mdot * (vj - v0) * v0)
    den = mdot_f * LHV_f
    return num / den
## -- CALCULATIONS 1) Walking trough engine cycle -----------------------------------------------------------------------------
# 1.1) Inlet 0->2

# 1.2) Fan 2->13


## -- CALCULATIONS 2) Gas generator & Efficiencies -----------------------------------------------------------------------------------