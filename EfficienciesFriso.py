import numpy as np
import math 
import numpy as np

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
