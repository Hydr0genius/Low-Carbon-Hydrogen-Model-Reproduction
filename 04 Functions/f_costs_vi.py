#!/usr/bin/env python
# coding: utf-8

# In[10]:

import numpy as np

def f_costs_vi(k_e, k_h, r, T, eta_c, F_h, SP_h, revenues, parameters, renewables):
    Delta_h = parameters['Delta_h']
    L = parameters['L']
    disc_degr = parameters['disc_degr']
    CF_h = revenues['CF_h']
#    CF_c = revenues['CF_c']
    w_o = revenues['w_o']

    p_s = revenues['p_s']
    p_b = revenues['p_b']
    Gamma_e = revenues['Gamma_e']
    CF_e = revenues['CF_e']
    ptc = revenues['ptc']
    LCOE = renewables['LCOE']
    
    Gamma_cs = revenues['Gamma_cs']
    Gamma_cb = revenues['Gamma_cb']
    CF_cs = revenues['CF_cs']
    CF_cb = revenues['CF_cb']
 
    # Parameters needed to import: eta_c, CF_cb, CF_cs, p_b, Gamma_cs, Gamme_cb, 
    
#     # Profit margin RES
#     D_e = p_s * Gamma_e + ptc - LCOE
# #    j_r = 0
#     j_c = 0
#     if D_e < 0:
# #        j_r = (CF_e * k_e)/(CF_r * k_h)
#         j_c = (CF_e * k_e)/(CF_c * k_h)

#     # Calculate LFCH
#     gamma = 1/(1+r)
#     if CF_h == 0:
#         f_h = 0
#         c_h = 0
#     else:

    f_h = (np.sum(F_h * parameters["disc"])) / (CF_h * L) if CF_h * L != 0 else 0 # replaced (gamma**T) with parameters["disc"] as returns TypeError, also fixed division by 0 problem for k_h = 0 
    c_h = SP_h * Delta_h / (CF_h * L) if CF_h * L != 0 else 0 # added Delta already in this term
    LFCH = f_h + (Delta_h * c_h)

#     # Calculate LCOE_h
#     if CF_r == 0:
#         f_r = 0
#         c_r = 0
#         LCOE_h = revenues['w_r']
#     else:
#         f_r = (np.sum(F_h * parameters["disc"]))/(CF_r * L) # Same here
#         c_r = (SP_h)/(CF_r * L)
#         LCOE_h = revenues['w_r'] + lambda_r * (f_r + Delta_h * c_r - j_r * D_e + 1/CF_r * LP)

    # Calculate LCOH - so far PTC not included in cost function, only NPV and revenues
#     if CF_h == 0:
#         f_c = 0
#         c_c = 0
#         LCOH_old = 1/eta_c * (revenues['w_c'])
#     else:
#     f_c = (np.sum(F_h * parameters["disc"]))/(CF_h * L)
#     c_c = (SP_h * Delta_h)/(CF_h * L)  # Added delta_h
#     LCOH_old = 1/eta_c * (revenues['w_c'] + (f_c + Delta_h * c_c - j_c * D_e + 1/CF_c)) # * LP)) 
    # Excluded Lambda_c - cost allocation should not be necessary anymore
        
        # LCOH new
    LCOH = w_o + 1/eta_c * (disc_degr * ((p_s * Gamma_cs * CF_cs) + (p_b * Gamma_cb * CF_cb)) / (disc_degr * CF_h)
           + f_h + c_h) if disc_degr * CF_h != 0 else 0 # Question: Is there a case f_h and c_h are not 0 if CF_h is 0? I don't think so            
    # Idea here: possibly replace p_s with variable OPEX of renewables + Transfer fee! 
    
    # Return summary
    costs = {
        'CF_e': CF_e,
        'F_e': renewables['F_e'],
        'SP_e': renewables['SP_e'],
        'p_s': p_s,
        'Gamma_e': Gamma_e,
        'f_e': renewables['f_e'],
        'c_e': renewables['c_e'],
        'Delta_e': renewables['Delta_e'],
        'LCOE': renewables['LCOE'],
        'LCOE_ptc': renewables['LCOE_ptc'],
#        'D_e': D_e,
        'eta_c': eta_c,
        'F_h': F_h,
        'SP_h': SP_h,
        'f_h': f_h,
        'c_h': c_h,
        'Delta_h': Delta_h,
        'LFCH': LFCH,
        'recon': 1/eta_c,
#        'w_c': revenues['w_c'],
#         'f_c': f_c,
#         'c_c': c_c,
#         'j_c': j_c,
  #      'LPadj_c': 1/CF_c * LP if CF_c != 0 else 0,
        'LCOH': LCOH
    }

    return costs

