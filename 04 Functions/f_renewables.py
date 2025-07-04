#!/usr/bin/env python
# coding: utf-8

# # Functions

# In[1]:
import numpy as np

def f_renewables(gamma, T, L, F_e, w_ei, SP_e, k_e, Delta_e, p_s_t, CF_e_t, ptc, parameters, t):
    CF_e = np.mean(CF_e_t) # CF renewable source(s)
    years = np.arange(0, len(T), 1) # calculated here, could also be imported in the function from f_parameters
    
    # Calculation of LCOE
    disc = gamma ** years # re-introduce discount factor, with years (list) instead of T (float) --- import years or calculate new 
    f_e = np.sum(F_e * disc) / (CF_e * L) if CF_e * L != 0 else 0 # replaced gamma ** T with disc # Brauchen wir hier nicht auch ne degradation rate? Also fixed division by 0
    c_e = SP_e / (CF_e * L) if CF_e * L != 0 else 0
    LCOE = w_ei + f_e + Delta_e * c_e
    LCOE_ptc = LCOE - ptc

    # Calculation of co-variance
    if CF_e == 0:  # REMOVE ERROR of dividing by zero
       epsilon_e_t = np.ones(len(CF_e_t)) * 0
    else:
       epsilon_e_t = CF_e_t / CF_e
    p_s = np.mean(p_s_t)
    mu_s_t = p_s_t / p_s
    Gamma_e = np.mean(epsilon_e_t * mu_s_t)
    
    # Contribution Margin Renewable - added w_ei to function (unequal to 0!), also k_e, which should usually just be 1
    CM_e_t = (p_s_t - w_ei + ptc) * CF_e_t * k_e
    CM_e = np.mean(CM_e_t)

    # Return Summary
    renewables = {
        'p_s': p_s,
        'Gamma_e': Gamma_e,
        'CF_e': CF_e,
        'F_e': F_e,
        'SP_e': SP_e,
        'f_e': f_e,
        'c_e': c_e,
        'Delta_e': Delta_e,
        'LCOE': LCOE,
        'LCOE_ptc': LCOE_ptc,
        # 'CM_e_t': CM_e_t,
        'CM_e': CM_e,
        'w_ei': w_ei
    }

    return renewables



# Notes: condition of discarding energy not yet introduced (should never be the case really, unless negative prices and no production hydrogen (will not be the case)

