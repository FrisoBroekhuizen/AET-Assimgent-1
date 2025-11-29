import numpy as np


def M_to_V(M, gamma, R, T):
    V = M * np.sqrt(gamma * R * T)
    return V

def total_temp(T_s, V, c_p_a):
    T_tot = T_s + (V**2)/(2*c_p_a)
    return T_tot

def total_p(p_s, V, c_p_a, T_tot, kappa_a):
    p_tot = p_s*(1 + (V**2/(2 * c_p_a * T_tot)))**(kappa_a/(kappa_a - 1))
    return p_tot

def p_tot_after_comp(p_tot_before, PR_comp):
    p_tot_after = p_tot_before * PR_comp
    return p_tot_after

def T_tot_after_comp(T_tot_before, eta_is_comp, PR_comp, kappa):
    T_tot_after = T_tot_before + (T_tot_before)/(eta_is_comp) * (PR_comp**((kappa - 1)/kappa)-1)
    return T_tot_after

def work_comp(mfr, c_p_a, T_tot_before, T_tot_after):
    work = mfr * c_p_a * (T_tot_after - T_tot_before)
    return work

def mfr_fuel(mfr, c_p_g, T_tot_before, T_tot_after, LHV_f, eta_cc):
    mfr_f = (mfr * c_p_g * (T_tot_after - T_tot_before))/(LHV_f * eta_cc)
    return mfr_f

def work_turb(work_comp, eta_mech):
    work = work_comp/eta_mech
    return work

def T_tot_turb(T_tot_before, W_turb, mfr, c_p_g):
    T_tot_after = T_tot_before - (W_turb)/(mfr * c_p_g)
    return T_tot_after

def p_tot_turb(p_tot_before, eta_turb, T_tot_before, T_tot_after, kappa_g):
    p_tot_after = p_tot_before * (1 - 1/eta_turb * (1 - (T_tot_after)/(T_tot_before)))^(kappa_g/(kappa_g - 1))
    return p_tot_after

def work_gas_gen(mfr, c_p_g, T_g, T_is_exp_amb):
    work = mfr * c_p_g * (T_g - T_is_exp_amb)
    return work

