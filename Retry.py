import numpy as np


############################### PARAMETERS ########################################

Mach = 0.78
Altitude = 10668        # m
T_ambient = 220       # K
P_ambient = 23842       # Pa
R_air = 287             # J/kg·K

# Mass flows
mdot_air = (8.6 + 1) * 50          # kg/s
Bypass_ratio = 8.6

# Pressure ratios
IPR_inlet = 0.99
FPR_fan = 1.5
PR_LPC = (40/(1.5*10))
PR_HPC = 10
PR_combustor = 0.96

# Efficiencies
eta_fan = 0.92
eta_LPC = 0.92
eta_HPC = 0.92
eta_LPT = 0.92
eta_HPT = 0.92
eta_mech = 0.99
eta_combustor = 0.995
eta_nozzle = 0.99
eta_gearbox = 1

# Turbine inlet temperature
Tt4 = 1450              # K

# Fuel properties
LHV_fuel = 43e6         # J/kg  (43 MJ/kg)

# Specific heats and specific heat ratios
cp_air = 1000           # J/kg·K
gamma_air = 1.4
cp_gas = 1150           # J/kg·K
gamma_gas = 1.33


########################################### Equations Engine Cycle ###########################################
def PR_crit(eta, kappa):
    """Determines the critical pressure ratio for choked flow."""
    PR = (1 - (1)/(eta) * ((kappa - 1)/(kappa + 1)))**((-1 * kappa)/(kappa - 1))
    return PR

def M_to_V(M, gamma, R, T):
    """Calculates velocity from Mach number."""
    V = M * np.sqrt(gamma * R * T)
    return V

def total_temp2(T_s, M, gamma_air):
    """Calculates total temp from static temp and M."""
    T_tot = T_s *(1 + (M**2)*(gamma_air-1)/(2))
    return T_tot

def total_p2(p_s, M, gamma_air):
    """Calculates total p from static p and M."""
    T_tot = p_s *((1 + (M**2)*(gamma_air-1)/(2)))**(gamma_air/(gamma_air - 1))
    return T_tot

def p_tot_after_comp(p_tot_before, PR_comp):
    """Calculates total p after compression."""
    p_tot_after = p_tot_before * PR_comp
    return p_tot_after

def T_tot_after_comp(T_tot_before, eta_is_comp, PR_comp, kappa):
    """Calculates total T after compression."""
    T_tot_after = T_tot_before + (T_tot_before)/(eta_is_comp) * (PR_comp**((kappa - 1)/kappa) - 1)
    return T_tot_after

def work_comp(mfr, c_p_a, T_tot_before, T_tot_after):
    """Calculates work done by compressor."""
    work = mfr * c_p_a * (T_tot_after - T_tot_before)
    return work

def mfr_fuel(mfr, c_p_a, c_p_g, T_tot_before, T_tot_after, LHV_f, eta_cc):
    """Calculates fuel mass flow."""
    mfr_f = ((mfr * c_p_g * T_tot_after) - (mfr * c_p_a * T_tot_before))/(LHV_f * eta_cc - c_p_g * T_tot_after)
    return mfr_f

def work_turb(work_comp, eta_mech):
    """Calculates work that turbine has to deliver."""
    work = work_comp/eta_mech
    return work

def T_tot_turb(T_tot_before, W_turb, mfr, c_p_g):
    """Calculates total T after turbine"""
    T_tot_after = T_tot_before - (W_turb)/(mfr * c_p_g)
    return T_tot_after

def p_tot_turb(p_tot_before, eta_turb, T_tot_before, T_tot_after, kappa_g):
    """Calculates total p after turbine."""
    p_tot_after = p_tot_before * (1 - 1/eta_turb * (1 - (T_tot_after)/(T_tot_before)))**(kappa_g/(kappa_g - 1))
    return p_tot_after

def work_gas_gen(mfr, c_p_g, T_g, T_is_exp_amb):
    """Calculates work output of gas generator."""
    work = mfr * c_p_g * (T_g - T_is_exp_amb)
    return work

def T_is_exp_amb(T_g, p_g, p_amb, kappa_g):
    """Calculates temperature of isentropic expansion."""
    T_tot = T_g*(((p_amb)/(p_g))**((kappa_g-1)/(kappa_g)))
    return T_tot


################################## Efficiency formulas #################################
def ETA_comb(mdot_air, mdot_f, mdot_gas, T3, T4, c_p_g, c_p_a, LHV_f):
    """
    Combustor efficiency:

    """
    return (mdot_gas * c_p_g * T4 - mdot_air * c_p_a * T3) / (mdot_f * LHV_f)

def ETA_thdy(P_gg, mdot_air, mdot_gas, T3, T4, c_p_g, c_p_a):
    """
    Alternative Thermodynamic efficiency of gas generator:
    """
    num = mdot_gas * c_p_g * T4 - mdot_air * c_p_a * T3
    return P_gg / num

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
    """
    mdot = np.asarray(mdot)
    vj   = np.asarray(vj)
    num = 0.5 * np.sum(mdot * (vj**2 - v0**2))
    den = mdot_f * LHV_f
    return num / den

def ETA_total(mdot, vj, v0, mdot_f, LHV_f):
    """
    Total efficiency:
    """
    mdot = np.asarray(mdot)
    vj   = np.asarray(vj)
    num = np.sum(mdot * (vj - v0) * v0)
    den = mdot_f * LHV_f
    return num / den