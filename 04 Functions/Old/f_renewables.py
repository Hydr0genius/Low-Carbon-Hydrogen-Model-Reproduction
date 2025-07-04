#!/usr/bin/env python
# coding: utf-8

# # Functions

# In[1]:
import numpy as np

def f_renewables(gamma, T, L, F_e, SP_e, Delta_e, p_s_t, CF_e_t, ptc, parameters, t):
    CF_e = np.mean(CF_e_t) # CF renewable source(s)

    # Calculation of LCOE
    f_e = np.sum(F_e * parameters["disc"]) / (CF_e * L) # replaced gamma ** T with disc # Brauchen wir hier nicht auch ne degradation rate?
    c_e = SP_e / (CF_e * L)
    LCOE = f_e + Delta_e * c_e
    LCOE_ptc = LCOE - ptc

    # Calculation of co-variance
    if CF_e == 0:  # REMOVE ERROR of dividing by zero
       epsilon_e_t = np.ones(len(CF_e_t)) * 0
    else:
       epsilon_e_t = CF_e_t / CF_e
    p_s = np.mean(p_s_t)
    mu_s_t = p_s_t / p_s
    Gamma_e = np.mean(epsilon_e_t * mu_s_t)

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
        'LCOE_ptc': LCOE_ptc
    }

    return renewables






