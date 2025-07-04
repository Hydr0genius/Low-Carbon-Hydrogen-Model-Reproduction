#!/usr/bin/env python
# coding: utf-8

# In[5]:
import numpy as np

def f_revenues_vi(k_e, k_h, CF_e_t, p_h, eta_c, eta_r, p_s_t, ptc, p_s_t_PtG, p_b_t, w_o, delta_h, LP_t, phi_t, PTC_h):
    # Initialize a dictionary to store the results
    revenues = {}

    # Hourly capacity factors
    p_h_c = eta_c * p_h
    CV_h = eta_c * (p_h - w_o + PTC_h) # included PTC_h here as it affects the conversion value of renewable power directly

#     # Renewable conversion old
#     w_cs_t = p_s_t_PtG + eta_c * w_o
#     if k_h == 0:  # Remove errors from dividing by zero
#         z_t = np.zeros(len(CF_e_t))
#     else:
#         z_t = np.minimum(CF_e_t * k_e / k_h, (1 - phi_t))
#     CF_cs_t = np.zeros(len(w_cs_t))
#     if min(CV_h, p_b_t) > p_s_t_PtG: # removed the loop over i
#         CF_cs_t = z_t

#     # Grid conversion
#     w_cb_t = p_b_t + eta_c * w_o
#     CF_cb_t = np.zeros(len(w_cb_t))
#     if CV_h > p_b_t: # Removed loop over i - returns error
#         CF_cb_t = 1 - CF_cs_t
        
    # Renewable conversion alternative
    w_cs_t = p_s_t_PtG + eta_c * w_o
    if k_h == 0:  # Remove errors from dividing by zero
        z_t = np.zeros(len(CF_e_t))
    else:
        z_t = np.minimum(CF_e_t * k_e / k_h) # ,(1 - phi_t)) # Excluded the phi condition so far
 
    CF_cs_t = np.zeros(len(w_cs_t))
    mask = np.minimum(CV_h, p_b_t) > p_s_t_PtG # new try to make it fit, otherwise the condition only holds for 6632 output values (and 0 without PTC!)
    CF_cs_t[mask] = z_t[mask]
    # CF_cs_t[np.minimum(CV_h, p_b_t) > p_s_t_PtG] = z_t # What is the point of that? 

    # Grid conversion - not needed for now - make unattractive
    w_cb_t = p_b_t + eta_c * w_o
    CF_cb_t = np.zeros(len(w_cb_t))
    CF_cb_t[CV_h > p_b_t] = 1 - CF_cs_t


    # Reconversion - also not needed for now
    w_r = (1 / eta_r) * (p_h + delta_h)
    CF_r_t = p_s_t_PtG > w_r

    # CHECK: complementary slackness
    CF_check_t = CF_cs_t + CF_cb_t + CF_r_t
    CF_check = np.sum(CF_check_t > 1)  # if > 0, then CF_c_t and CF_r_t = 1 at the same time!

    # Average CFs
    CF_cs = np.mean(CF_cs_t)
    CF_cb = np.mean(CF_cb_t)
    CF_c = np.mean(CF_cs_t + CF_cb_t)
    CF_r = np.mean(CF_r_t)
    CF_h = np.mean(CF_check_t)  # CF_cb + CF_cs + CF_r

    beta_cs = CF_cs / CF_h
    beta_cb = CF_cb / CF_h
    beta_c = CF_c / CF_h
    beta_r = CF_r / CF_h

    # Co-variation
    if CF_cs == 0:  # Remove error of dividing by zero
        epsilon_cs_t = np.ones(len(CF_cs_t)) * 0 
    else:
        epsilon_cs_t = CF_cs_t / CF_cs
    w_cs = np.mean(w_cs_t)
    mu_cs_t = w_cs_t / w_cs
    Gamma_cs = np.mean(epsilon_cs_t * mu_cs_t)

    if CF_cb == 0:  # Remove error of dividing by zero
        epsilon_cb_t = np.ones(len(CF_cb_t)) * 0
    else:
        epsilon_cb_t = CF_cb_t / CF_cb
    w_cb = np.mean(w_cb_t)
    mu_cb_t = w_cb_t / w_cb
    Gamma_cb = np.mean(epsilon_cb_t * mu_cb_t)

    if CF_r == 0:  # Remove error of dividing by zero
        epsilon_r_t = np.ones(len(CF_r_t)) * 0
    else:
        epsilon_r_t = CF_r_t / CF_r
    p_s_PtG = np.mean(p_s_t_PtG)
    mu_s_t_PtG = p_s_t_PtG / p_s_PtG
    Gamma_r = np.mean(epsilon_r_t * mu_s_t_PtG)

    CF_e = np.mean(CF_e_t)
    if CF_e == 0:  # Remove error of dividing by zero
        epsilon_e_t = np.ones(len(CF_e_t)) * 0
    else:
        epsilon_e_t = CF_e_t / CF_e
    p_s = np.mean(p_s_t)
    mu_s_t = p_s_t / p_s
    Gamma_e = np.mean(epsilon_e_t * mu_s_t)

    # Contribution margins
    LP = np.mean(LP_t)
    CM_cs = (eta_c * p_h - w_cs * Gamma_cs) * CF_cs * k_h
    CM_cb = (eta_c * p_h - w_cb * Gamma_cb) * CF_cb * k_h
    CM_c = CM_cs + CM_cb
    CM_r = (p_s_PtG * Gamma_r - w_r) * CF_r * k_h
    CM_h = CM_cs + CM_cb + CM_r

    w_c = CF_cs/CF_c * w_cs * Gamma_cs + CF_cb/CF_c * w_cb * Gamma_cb

    # Cost allocation
    if CM_h == 0:
        lambda_c = 0
        lambda_r = 0
    else:
        lambda_c = CM_c / CM_h
        lambda_r = CM_r / CM_h

    # Return summary
    revenues['eta_c'] = eta_c
    revenues['eta_r'] = eta_r
    revenues['p_s'] = p_s
    revenues['Gamma_e'] = Gamma_e
    revenues['CF_e'] = CF_e
    revenues['ptc'] = ptc
    revenues['k_h'] = k_h
    revenues['p_h'] = p_h
    revenues['p_h_c'] = p_h_c
    revenues['delta_h'] = delta_h
    revenues['w_c'] = w_c
    revenues['w_c_c'] = w_c * 1/eta_c
    revenues['w_cs'] = w_cs
    revenues['w_cb'] = w_cb
    revenues['p_s_PtG'] = p_s_PtG
    revenues['w_r'] = w_r
    revenues['CF_h'] = CF_h
    revenues['CF_c'] = CF_c
    revenues['CF_cs'] = CF_cs
    revenues['CF_cb'] = CF_cb
    revenues['CF_r'] = CF_r
    revenues['CF_check'] = CF_check
    revenues['beta_c'] = beta_c
    revenues['beta_cs'] = beta_cs
    revenues['beta_cb'] = beta_cb
    revenues['beta_r'] = beta_r
    revenues['Gamma_cs'] = Gamma_cs
    revenues['Gamma_cb'] = Gamma_cb
    revenues['Gamma_r'] = Gamma_r
    revenues['CM_h'] = CM_h
    revenues['CM_c'] = CM_c
    revenues['CM_cs'] = CM_cs
    revenues['CM_cb'] = CM_cb
    revenues['CM_r'] = CM_r
    revenues['lambda_c'] = lambda_c
    revenues['lambda_r'] = lambda_r
    revenues['LP'] = LP

    return revenues


