import numpy as np
import math 
## ------------------------------------------ Flight condition--------------------------------
#constants from TUTORIAL:
# LEAP-1A @ cruise (Mach 0.78, Altitude 10668 m)
'''
Mach = 0.78
Altitude = 10668          # m

# Ambient conditions
T_ambient = 218.8         # K
P_ambient = 23842         # Pa
R_air = 287               # J/kg·K

# Mass flow & bypass
mdot_air = 173            # kg/s
Bypass_ratio = 12

# Pressure ratios
IPR_inlet = 0.98
FPR_fan = 1.4
PR_LPC = 1.7
PR_HPC = 12.5
PR_combustor = 0.96

# Temperatures
Tt4 = 1400                # K (combustor exit)

# Efficiencies
eta_fan = 0.90
eta_LPC = 0.92
eta_HPC = 0.92
eta_LPT = 0.90
eta_HPT = 0.90
eta_mech = 0.99
eta_combustor = 0.995
eta_nozzle = 0.98
eta_gearbox = 0.995

# Fuel properties
LHV_fuel = 43e6           # J/kg  (43 MJ/kg)

# Specific heats & gammas
cp_air = 1000             # J/kg·K
gamma_air = 1.4
cp_gas = 1150             # J/kg·K
gamma_gas = 1.33
'''


#Constant from assignment:

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
def PR_crit(eta, kappa):
    """Determines the critical pressure ratio for choked flow."""
    PR = (1 - (1)/(eta) * ((kappa - 1)/(kappa + 1)))**((-1 * kappa)/(kappa - 1))
    return PR

def M_to_V(M, gamma, R, T):
    """Calculates velocity from Mach number."""
    V = M * np.sqrt(gamma * R * T)
    return V

def total_temp1(T_s, V, c_p_a):
    """Calculates total Temp from static temp and velocity."""
    T_tot = T_s + (V**2)/(2*c_p_a)
    return T_tot

def total_temp2(T_s, M, gamma_air):
    """Calculates total temp from static temp and M."""
    T_tot = T_s *(1 + (M**2)*(gamma_air-1)/(2))
    return T_tot

def total_p1(p_s, V, c_p_a, T_tot, kappa_a):
    """Calculates total p from static p, total T and V"""
    p_tot = p_s*(1 + (V**2/(2 * c_p_a * T_tot)))**(kappa_a/(kappa_a - 1))
    return p_tot

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
    T_tot_after = T_tot_before + (T_tot_before)/(eta_is_comp) * (PR_comp**((kappa - 1)/kappa)-1)
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


## -------------------------------- Efficiency formulas --------------------------
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

## -- CALCULATIONS 1) Walking trough engine cycle -----------------------------------------------------------------------------
# 1.1) Inlet 0->2
# Ambient:
T_tota = total_temp2(T_ambient, Mach, gamma_air)
p_tota = total_p2(P_ambient, Mach, gamma_air)
T_tot2 = T_tota
p_tot2 = p_tot_after_comp(p_tota, IPR_inlet)

mdot_core = mdot_air / (1 + Bypass_ratio)
mdot_bypass = mdot_air - mdot_core

# 1.2) Fan 2->13
p_tot13 = p_tot_after_comp(p_tot2, FPR_fan)
p_tot21 = p_tot13
T_tot13 = T_tot_after_comp(T_tot2, eta_fan, FPR_fan, gamma_air)
T_tot21 = T_tot13

#work done on fan
W_fan = work_comp(mdot_air, cp_air, T_tot2, T_tot13)
print("Work done on fan (W): ", W_fan)
p_tot25 = p_tot_after_comp(p_tot21, PR_LPC)
T_tot25 = T_tot_after_comp(T_tot21, eta_LPC, PR_LPC, gamma_air) 
# work done on LPC
W_LPC = work_comp(mdot_core, cp_air, T_tot21, T_tot25)
print("Work done on LPC (W): ", W_LPC)

#1.3) HPC 25->3
p_tot3 = p_tot_after_comp(p_tot25, PR_HPC)
T_tot3 = T_tot_after_comp(T_tot25, eta_HPC, PR_HPC, gamma_air)
# work done on HPC
W_HPC = work_comp(mdot_core, cp_air, T_tot25, T_tot3)
print("Work done on HPC (W): ", W_HPC)
# 1.4) Combustor 3->4
p_tot4 = p_tot3 * PR_combustor
T_tot4 = Tt4
# fuel mass flow
mdot_fuel = mfr_fuel(mdot_core, cp_air,  cp_gas, T_tot3, T_tot4, LHV_fuel, eta_combustor)
print("Fuel mass flow (kg/s): ", mdot_fuel)
mdot_core_fuel = mdot_core + mdot_fuel
print("Total mass flow at combustor exit (kg/s): ", mdot_core_fuel)
# 1.5) HPT 4->45
W_HPT = W_HPC/eta_mech
T_tot45 = T_tot4-(W_HPT)/(mdot_core_fuel * cp_gas)
print("T_tot after HPT (K): ", T_tot45)
p_tot45 = p_tot_turb(p_tot4, eta_HPT, T_tot4, T_tot45, gamma_gas)
print("p_tot after HPT (Pa): ", p_tot45)
# 1.6) LPT 45->5
W_fangeared = W_fan / eta_gearbox
W_LPT = (W_LPC + W_fangeared)/eta_mech
mdot_45 = mdot_core_fuel
T_tot5 = T_tot45 - (W_LPT)/(mdot_45 * cp_gas)
print("T_tot after LPT (K): ", T_tot5)
p_tot5 = p_tot_turb(p_tot45, eta_LPT, T_tot45, T_tot5, gamma_gas)
print("p_tot after LPT (Pa): ", p_tot5)
# 1.7) Core nozzle 
p_tot7 = p_tot5
T_tot7 = T_tot5
nozlleratio = p_tot7 / P_ambient
print("Nozzle pressure ratio: ", nozlleratio) 
criticalcore = PR_crit(eta_nozzle, gamma_gas) 
if nozlleratio > criticalcore:
    print("CHOKED!")
    T_8 = T_tot7 * (2/(gamma_gas + 1))
    p_8 = p_tot7 / (criticalcore)
    print("p_8 (Pa): ", p_8)
    print("T_8 (K): ", T_8)
    V_8 = np.sqrt(gamma_gas * R_air * T_8)
    print("V_8 (m/s): ", V_8)
    rho_8 = p_8 / (R_air * T_8)
    print("rho_8 (kg/m3): ", rho_8)
    A_8 = mdot_45 / (rho_8 * V_8)
    print("A_8 (m2): ", A_8)
    v_flight = M_to_V(Mach, gamma_air, R_air, T_ambient)
    print("Flight velocity (m/s): ", v_flight)
    F_core = mdot_45 * (V_8 -  v_flight) + A_8 * (p_8 - P_ambient)
    print("Core thrust (N): ", F_core)
else:
    print("UNCHOKED!")
    T_8 = T_tot7 * (2/(gamma_gas + 1))
    p_8 = p_tot7 / (criticalcore)
    print("p_8 (Pa): ", p_8)
    print("T_8 (K): ", T_8)
    v_flight = M_to_V(Mach, gamma_air, R_air, T_ambient)
    print("Flight velocity (m/s): ", v_flight)
    V_8 = np.sqrt(gamma_gas * R_air * T_8)
    print("V_8 (m/s): ", V_8)
    F_core = mdot_45 * (V_8 -  v_flight)
    print("Core thrust (N): ", F_core)
#1.8) Bypass nozzle
p_tot16 = p_tot13
T_tot16 = T_tot13
bypasspressureratio = p_tot16 / P_ambient
print("Bypass nozzle pressure ratio: ", bypasspressureratio)
criticalBypass = PR_crit(eta_nozzle, gamma_air)
if bypasspressureratio > criticalBypass:
    #choked
    print("CHOKED!")
    T_18 = T_tot16 * (2/(gamma_air + 1))
    p_18 = p_tot16 / (criticalBypass)
    print("p_18 (Pa): ", p_18)
    print("T_18 (K): ", T_18)
    V_18 = np.sqrt(gamma_air * R_air * T_18)
    print("V_18 (m/s): ", V_18)
    rho_18 = p_18 / (R_air * T_18)
    print("rho_18 (kg/m3): ", rho_18)
    A_18 = mdot_bypass / (rho_18 * V_18)
    print("A_18 (m2): ", A_18)
F_bypass = mdot_bypass * (V_18 - v_flight) + A_18 * (p_18 - P_ambient)
print("Bypass thrust (N): ", F_bypass)
# 1.9) Total thrust
F_total = F_core + F_bypass
print("Total thrust (N): ", F_total)
TSFC = mdot_fuel / (F_total*0.000001) 
print("TSFC (kg/Ns): ", TSFC)

## -- CALCULATIONS 2) Gas generator & Efficiencies -----------------------------------------------------------------------------------

# Gas generator
W_g = work_comp(mdot_core, cp_air, T_tot21, T_tot25)/eta_mech + work_comp(mdot_core, cp_air, T_tot2, T_tot21)/(eta_mech*eta_gearbox)
print("Power output of gas generator (W): ", W_g)
print("T_45", T_tot45)
T_g = T_tot_turb(T_tot45, W_g, mdot_core_fuel, cp_gas)
print("T_g (K): ", T_g)
print("T_5", T_tot5)

p_g = p_tot_turb(p_tot45, eta_LPT, T_tot45, T_g, gamma_gas)
print("P_g (Pa): ", p_g)

T_8_is_exp_amb = T_is_exp_amb(T_g, p_g, P_ambient, gamma_gas)
print("T_8_is_exp_amb (K): ", T_8_is_exp_amb)

W_gg = work_gas_gen(mdot_core_fuel, cp_gas, T_g, T_8_is_exp_amb)
print("Work output of gas generator (W): ", W_gg)
P_gg = W_gg -0.5*mdot_core*v_flight**2
print("Power output of gas generator (W): ", P_gg)

# Effective jet velocities
v_jet_eff_core = V_8 + (A_8/mdot_core_fuel)*(p_8 - P_ambient) 
V_jet_eff_bypass = V_18 + (A_18/mdot_bypass)*(p_18 - P_ambient)
print("Effective jet velocity for efficiency calculations (m/s): ", v_jet_eff_core, "<- core | bypass ->", V_jet_eff_bypass)

# Combustion efficiency
eta_comb = ETA_comb(mdot_core, mdot_fuel, mdot_core_fuel, T_tot3, T_tot4, cp_gas, cp_air, LHV_fuel)
print("Combustor efficiency: ", eta_comb)

# Thermodynamic efficiency
eta_thdy = ETA_thdy(P_gg, mdot_core, mdot_core_fuel, T_tot3, T_tot4, cp_gas, cp_air)
print("Thermodynamic efficiency: ", eta_thdy)

# Jet-generation efficiency
eta_jet_gen = 0.5*(mdot_core_fuel*(v_jet_eff_core**2-v_flight**2)+mdot_bypass*(V_jet_eff_bypass**2-v_flight**2))/(P_gg)
print("Jet-generation efficiency: ", eta_jet_gen)

# Propulsive efficiency
eta_prop = ETA_prop([mdot_core_fuel, mdot_bypass], [v_jet_eff_core, V_jet_eff_bypass], v_flight)
print("Propulsive efficiency: ", eta_prop)

# Thermal efficiency
# eta_thermal = ((0.5*mdot_core_fuel*(v_jet_eff_core**2-v_flight**2))+(0.5*mdot_bypass*(V_jet_eff_bypass**2-v_flight**2)))/(mdot_fuel * LHV_fuel)
eta_thermal = ETA_thermal([mdot_core_fuel, mdot_bypass], [v_jet_eff_core, V_jet_eff_bypass], v_flight, mdot_fuel, LHV_fuel)
print("Overall thermal efficiency: ", eta_thermal)

# Total efficiency
eta_total = ETA_total([mdot_core_fuel, mdot_bypass], [v_jet_eff_core, V_jet_eff_bypass], v_flight, mdot_fuel, LHV_fuel)
print("Total efficiency: ", eta_total)
print("Thermal efficiency consistency check:", eta_thermal, eta_thdy*eta_jet_gen*eta_comb)
print("Total efficiency consistency check:",eta_total, eta_prop*eta_thermal)



