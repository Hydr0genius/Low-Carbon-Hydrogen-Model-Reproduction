#!/usr/bin/env python
# coding: utf-8

import numpy as np


def f_npv_vi(k_e, k_h, alpha, eta_c, parameters, revenues, costs, renewables): # excluded p_h for now - maybe replace revenue p_h later
    L = parameters['L']
    disc = parameters['disc']
    x_vect = parameters['x_vect']
    disc_degr = parameters['disc_degr']
    Delta_e = parameters['Delta_e']
    Delta_h = parameters['Delta_h']
    alpha = parameters['alpha']
    ptc_sum = parameters['ptc_sum']
    ptc_h_sum = parameters['ptc_h_sum']
    PTC_vector = parameters['PTC_vector']
    PTC_h_vector = parameters['PTC_h_vector']
    PTC_grid_vector = parameters['PTC_grid_vector']
    m = parameters['m']

    CF_e = revenues['CF_e']
    Gamma_e = revenues['Gamma_e']
    p_s = revenues['p_s']
    ptc = revenues['ptc']
    w_o = revenues['w_o']
    p_b = revenues["p_b"]
    eta_c = revenues['eta_c']
    
    LCOE = costs['LCOE']
    LCOH = costs['LCOH']
#     D_e = costs['D_e']
    SP_h = costs['SP_h']
    F_h = costs['F_h']

    p_h = revenues['p_h']
    PTC_h = revenues["PTC_h"]
    PTC_grid = revenues["PTC_grid"] 
    PTC_grid_vector_actual = revenues["PTC_grid_vector_actual"]
    w_cs = revenues['w_cs'] # Currently still in 
    w_cb = revenues['w_cb']
    Gamma_cs = revenues['Gamma_cs']
    Gamma_cb = revenues['Gamma_cb']
    CF_cs = revenues['CF_cs']
    CF_cb = revenues['CF_cb']
    CF_h = revenues['CF_h']
    CF_e_remaining = revenues['CF_e_remaining']
    
    w_ei = renewables['w_ei']
    F_e = renewables['F_e']
    SP_e = renewables['SP_e']
    CF_e = renewables['CF_e']

#    p_s_PtG = revenues['p_s_PtG']
#    Gamma_r = revenues['Gamma_r']
#    w_r = revenues['w_r']
#    CF_r = revenues['CF_r']

    LFCH = costs['LFCH']
    CF_h = revenues['CF_h']
#     LP = revenues['LP']

    # --------------- New Approach ---------------------------------------------------------#
    # Even necessary? We can calculate the same interim KPIs as before!
#     NPV_h4 = (1 - alpha) * 
    
    # --------------- Old Approach ---------------------------------------------------------#
    
    # Help vectors (as long as power price does not change over the years)
    p_s_vector = np.array([p_s] * len(PTC_vector))
    Gamma_e_vector = np.array([Gamma_e] * len(PTC_vector))
    p_h_vector = np.array([p_h] * len(PTC_h_vector))
    
    # LROE and LROH
#     LROE = m * disc_degr * CF_e * (p_s * Gamma_e + ptc) / (m * disc_degr * CF_e) # v1
    LROE = (
        m * np.sum(disc * x_vect * CF_e_remaining * (p_s_vector * Gamma_e_vector + PTC_vector)) # CF_e_remaining not yet functional!!
        / (m * disc_degr * CF_e_remaining) if m * disc_degr * CF_e != 0 else 0
    ) # v2 Now added the CF for selling the power to the grid, which is obviously much lower than the previous full CF_e 
    
    # Split term into grid conversion and renewable conversion # Add PTC vector term instead (as in analytical model) 
    # -> so calculate for every year over the lifetime and sum up values
#     LROH = m * eta_c * disc_degr * (CF_cs * (p_h + PTC_h) 
#            + CF_cb * (p_h + PTC_grid)) / (m * eta_c * disc_degr * CF_h) # v1
    
    LROH = (
        m * eta_c * np.sum(disc * x_vect * (CF_cs * (p_h_vector + PTC_h_vector) + CF_cb * (p_h_vector + PTC_grid_vector_actual))) 
        / (m * eta_c * disc_degr * CF_h) if m * eta_c * disc_degr * CF_h != 0 else 0
    ) #v2 # changed to actual grid vector
    
#     LROH = (
#         m * eta_c * np.sum(disc * x_vect * (CF_h * (p_h_vector + PTC_grid_vector_actual))) # PTC_grid_vector_actual not yet functional
#         / (m * eta_c * disc_degr * CF_h) if m * eta_c * disc_degr * CF_h != 0 else 0
#     ) #v3 # Now with one common CF as PTC is calculated directly in f_revenues
    
    
    # NPV_v4 @GG how could we adapt that with CM inside? Don't currently see it in terms of PTC vectors etc...
    # NPV_h_new = (1 - alpha) * np.sum(disc * x_vect * CM_total) ...
                                                      
    # NPV_v3
    NPV_h = (1 - alpha) * ((L * CF_e * k_e * (LROE - LCOE)) 
            + L * eta_c * CF_h * k_h * (LROH - LCOH))
        
#     # NPV v2
    
#     NPV_h2 = ((1 - alpha) * (L * (CF_e * k_e * (p_s * Gamma_e + ptc - w_ei) - F_e * k_e)
#             + CF_cs * k_h * eta_c * (p_h + PTC_h - w_o - (p_s / eta_c) * Gamma_cs)
#             + CF_cb * k_h * eta_c * (p_h + PTC_grid - w_o - (p_b / eta_c) * Gamma_cb)
#             - F_h * k_h) 
#             - k_e * SP_e * Delta_e - k_h * SP_h * Delta_h)
        
    # Can we just not the joint disc_degr for small gamma and x as degradation is the same? Even better L (also includes m and disc_degr!)
    # Problem for now: F_e and k_e are also multiplied not only by aggr. disc factor but also degradation  (as one term) - fix later
    # We have to distinguish the h2 conversion from grid and renewable source
     
      # NPV original
#     NPV_h = ((Gamma_e * p_s + ptc - LCOE) * CF_e * k_e
#             + ((eta_c * (p_h + PTC_h) - w_cs * Gamma_cs) * CF_cs # Added PTC_h Component here
#             + (eta_c * p_h + PTC_grid - w_cb * Gamma_cb) * CF_cb
#             # + (p_s_PtG * Gamma_r - w_r) * CF_r
#             - LFCH * CF_h - LP)  * k_h)  # * (1 - alpha) * L 

    # Check
#     PM_e = D_e * CF_e * k_e
#     PM_c = eta_c * (p_h + PTC_h) - eta_c * costs['LCOH'] # and here
# #    PM_r = p_s_PtG * Gamma_r - costs['LCOE_h']
#     if PM_e > 0:
#         NPV_h2 = PM_e + PM_c * CF_h * k_h # + PM_r * CF_r * k_h
#     else:
#         NPV_h2 = 2*PM_e + PM_c * CF_h * k_h # + PM_r * CF_r * k_h

    NPV = {
        'k_e': k_e,
        'k_h': k_h,
        'NPV_h': NPV_h,
#        'NPV_h2': NPV_h2,
#         'PM_e': PM_e,
#         'PM_c': PM_c,
        'con': eta_c * p_h,
        'con_LCOH': eta_c * costs['LCOH'],
#        'PM_r': PM_r,
#        'recon': p_s_PtG * Gamma_r,
#        'LCOE_h': costs['LCOE_h'],
#        'NPV_h2': NPV_h2,
        'LROE': LROE,
        'LROH': LROH,
        "LCOH": LCOH,
        'LCOE': LCOE
    }

    return NPV

