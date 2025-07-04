#!/usr/bin/env python
# coding: utf-8

# In[10]:


def f_npv_vi(k_e, k_h, alpha, eta_c, parameters, revenues, costs):
    L = parameters['L']

    CF_e = revenues['CF_e']
    Gamma_e = revenues['Gamma_e']
    p_s = revenues['p_s']
    ptc = revenues['ptc']
    LCOE = costs['LCOE']
    D_e = costs['D_e']

    p_h = revenues['p_h']
    w_cs = revenues['w_cs']
    w_cb = revenues['w_cb']
    Gamma_cs = revenues['Gamma_cs']
    Gamma_cb = revenues['Gamma_cb']
    CF_cs = revenues['CF_cs']
    CF_cb = revenues['CF_cb']
    CF_c = revenues['CF_c']

    p_s_PtG = revenues['p_s_PtG']
    Gamma_r = revenues['Gamma_r']
    w_r = revenues['w_r']
    CF_r = revenues['CF_r']

    LFCH = costs['LFCH']
    CF_h = revenues['CF_h']
    LP = revenues['LP']

    NPV_h = ((Gamma_e * p_s + ptc - LCOE) * CF_e * k_e
            + ((eta_c * p_h - w_cs * Gamma_cs) * CF_cs
            + (eta_c * p_h - w_cb * Gamma_cb) * CF_cb
            + (p_s_PtG * Gamma_r - w_r) * CF_r
            - LFCH * CF_h - LP)  * k_h)  # * (1 - alpha) * L 

    # Check
    PM_e = D_e * CF_e * k_e
    PM_c = eta_c * p_h - eta_c * costs['LCOH']
    PM_r = p_s_PtG * Gamma_r - costs['LCOE_h']
    if PM_e > 0:
        NPV_h2 = PM_e + PM_c * CF_c * k_h + PM_r * CF_r * k_h
    else:
        NPV_h2 = 2*PM_e + PM_c * CF_c * k_h + PM_r * CF_r * k_h

    NPV = {
        'k_e': k_e,
        'k_h': k_h,
        'NPV_h': NPV_h,
        'PM_e': PM_e,
        'PM_c': PM_c,
        'con': eta_c * p_h,
        'con_LCOH': eta_c * costs['LCOH'],
        'PM_r': PM_r,
        'recon': p_s_PtG * Gamma_r,
        'LCOE_h': costs['LCOE_h'],
        'NPV_h2': NPV_h2,
    }

    return NPV

