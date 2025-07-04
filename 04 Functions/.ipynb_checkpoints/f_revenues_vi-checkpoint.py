#!/usr/bin/env python
# coding: utf-8

# In[5]:
import numpy as np

def f_revenues_vi(k_e, k_h, CF_e_t, p_h, eta_c, p_s_t, ptc, p_b_t, w_o, PTC_h, PTC_grid, c_i_t, CT_h): #
    # Initialize a dictionary to store the results 
    revenues = {}
    
    # Introduce new variables for new approach -------------------------------------------------------# 
#     g_w_t = k_e * rho * CF_e_t_w # hourly generation potential wind [in kwh]
#     g_s_t = k_e * (1 - rho) * CF_e_t_s # hourly generation potential solar [in kwh]
    
    # introduce v1 and v2 as float variables that are dependent on each other - define as a fraction of CF_e_t? 
    # And then     
    
#     # Max contribution margin - first try in main model
#     CM_i_t = (p_s_t + ptc) * v_1_t - v_3_t * p_b_t + eta_c * (v_2_t + v_3_t) * [p_h - w_o + PTC_h]
#     # CM_i_t = (p_s_t + ptc) * v_w_t + (p_s_t + ptc) * v_s_t + eta_c * (c_w_t + c_s_t + c_b_t) * [p_h - w_o] + PTC_h
    
#     CI_t = (v_3_t * c_i_t) / [eta_c * (v_3_t + v_2_t)] # New equivalent) form of scaling the hourly carbon intensity
    
#     PTC_h_values = np.select(
#         [
#             CI_t <= 0.45,
#             (CI_t > 0.45) & (CI_t <= 1.5),
#             (CI_t > 1.5) & (CI_t <= 2.5),
#             (CI_t > 2.5) & (CI_t <= 4),
#             CI_t > 4
#         ],
#         [3, 1, 0.75, 0.6, 0]
#     ) 
    
#     # Constraints: 
#     v_1_t + v_2_t <= CF_e_t * k_e
#     v_2_t + v_3_t <= k_h
    
    # CM_i = np.mean(CM_i_t)
    
    # --------------- Old Approach ---------------------------------------------------------#
    
    # Hourly capacity factors
#     p_h_c = eta_c * p_h # Not needed in this version
#    CV_h = eta_c * (p_h + PTC_h - w_o) # included PTC_h here as it affects the conversion value of renewable power directly. 
    # Attention: CV_p later not used. PTC_h needs to be introduced later as well for LCOH

    # Renewable Conversion NEW TEST
    CV_h = eta_c * (p_h + PTC_h - w_o)
    w_cs_t = [ps + eta_c * w_o for ps in p_s_t]
    z_t = [min(cf * k_e / k_h, 1) for cf in CF_e_t] if k_h != 0 else np.zeros(len(CF_e_t))
    
    CM_condition_renewable = CV_h >= np.array(p_s_t)
#    renewable_conversion_condition = np.array(z_t) * k_h <= np.array(CF_e_t) * k_e
    CF_cs_t = np.where(CM_condition_renewable, z_t, 0) # Take renewable condition out
#     CF_cs_t = np.minimum(CF_cs_t, 1)  # Ensure CF_cs_t is between 0 and 1, should not be necessary after z_t introduction

    # Sell rest of renewables to grid (so can for instance adapt CF_e_t!) @GG
    CF_e_t_remaining = np.array(CF_e_t) * k_e - np.array(CF_cs_t) * k_h
    CF_e_remaining = np.mean(CF_e_t_remaining)

    # Grid Conversion with carbon intensity Measure
    
    CI_full_grid = c_i_t / eta_c # Full Carbon intensity
    
    CI_actual = CI_full_grid * ((1 - CF_cs_t) / k_h) if k_h != 0 else CI_full_grid # For the conversion condition # Newly scaled by k_h @GG # Careful division by 0 
    
    # CT_h = 0.449 # Define Carbon Threshold - rather import from main script
    
    MSE_t = np.where(CT_h >= CI_full_grid, 1, 0) # Max sourceable electricity rate help variable
    
#     MSE_t = np.where(CF_cs_t > 0, CT_h / CI_actual, 0) # Renewable power in mix increases potential grid part 
#     # (alternative: run with CI_full_grid instead to ensure respecting threshold?
#     MSE_t = np.minimum(1 - CF_cs_t, MSE_t) # Dependent on how much we can fill with CF_cs_t - take the lower end
    
    MSE_t = np.where(CF_cs_t > 0, np.minimum(CT_h / CI_actual, 1 - CF_cs_t), 0) # scale by k_h?  
    MSE = np.mean(MSE_t)
    
    # Alternative: Only draw as much as to reach the 0.45kg --- could be calculated as an alternative e.g. "Minimum Carbon"
    
    # Try with updated CI_actual based on our max sourcable grid parameter
    CI_actual_PTC = CI_actual * MSE_t # CT_h Alternative: calculate hourly: CI_full_grid * MSE_t
       
#     # Calculate the maximum allowable CF_cb_t based on the carbon intensity threshold - Maybe add condition before?
#     max_CF_cb_t = 0.45 / carbon_intensity_per_kg_H2
#     max_CF_cb_t = np.where(max_CF_cb_t > 1, 1, max_CF_cb_t)  # Ensure max_CF_cb_t does not exceed 1


    # Determine PTC_grid based on carbon intensity # In current update should be uniform with carbon thresholds
    PTC_grid_values = np.select(
        [
            CI_actual_PTC <= 0.45,
            (CI_actual_PTC > 0.45) & (CI_actual_PTC <= 1.5),
            (CI_actual_PTC > 1.5) & (CI_actual_PTC <= 2.5),
            (CI_actual_PTC > 2.5) & (CI_actual_PTC <= 4),
            CI_actual_PTC > 4
        ],
        [3, 1, 0.75, 0.6, 0]
    )

    # Grid Conversion with updated PTC_grid and carbon intensity condition
    w_cb_t = [pb + eta_c * w_o for pb in p_b_t]
    
    # Introduce new clause for non hourly matching (so add REC penalty + define new CI that triggers that if clause)
    if CT_h == 100:
        CM_condition_grid = eta_c * (p_h + PTC_grid_values - np.array(p_b_t) / eta_c - w_o) * k_h > 0
        CF_cb_t = np.where(CM_condition_grid, MSE_t, 0) # Replaced 1 - CF_cs_t with MSE_t now as well here - should fix jump
        PTC_grid = np.mean(PTC_grid_values)
        CM_cb_t = eta_c * (p_h + PTC_grid_values - p_b_t / eta_c - w_o) * CF_cb_t * k_h
    elif CT_h == 101:
        CM_condition_grid = eta_c * (p_h + PTC_h - np.array(p_b_t) / eta_c - w_o) * k_h > 0 #  + 0.02 e.g. Here think about reasonable surcharge
        CF_cb_t = np.where(CM_condition_grid, 1 - CF_cs_t, 0)
        PTC_grid = PTC_h
        CM_cb_t = eta_c * (p_h + PTC_h - p_b_t / eta_c - w_o) * CF_cb_t * k_h
    else:
        CM_condition_grid = eta_c * (p_h + PTC_grid_values - np.array(p_b_t) / eta_c - w_o) * k_h > 0
        CF_cb_t = np.where(CM_condition_grid, MSE_t, 0)
        PTC_grid = np.mean(PTC_grid_values)
        CM_cb_t = eta_c * (p_h + PTC_grid_values - p_b_t / eta_c - w_o) * CF_cb_t * k_h
    
    # RECENT VERSIONS
    # CM_condition_grid = eta_c * (p_h + PTC_grid_values - np.array(p_b_t) / eta_c - w_o) * k_h > 0
    # CM_condition_grid_2 = eta_c * (p_h - np.array(p_b_t) / eta_c - w_o) * k_h > 0 # Without PTC for hours with low power prices
    
    # CF_cb_t = np.where(CM_condition_grid, MSE_t, 0) # economic check + bound by carbon threshold!
    # CF_cb_t = np.where(CM_condition_grid_2, 1 - CF_cs_t, 0) # Add instances with low power prices even if that means no PTC! Add in analytical model + test numerically 
    
    # CF_cb_t = np.where(CI_actual < CT_h, np.minimum(1 - CF_cs_t, np.minimum(CT_h / CI_full_grid, 1)), 0)
    
    # CF_cb_t = np.minimum(CF_cb_t, max_CF_cb_t)  # Limit CF_cb_t based on carbon intensity condition
    
    # Make dependent on scenario!
#    PTC_grid = np.mean(PTC_grid_values)
    
#     # Grid Conversion 2.0
#     w_cb_t = [pb + eta_c * w_o for pb in p_b_t]
#     CM_condition_grid = eta_c * (p_h + PTC_grid - np.array(p_b_t) / eta_c - w_o) * k_h > 0
#     CF_cb_t = np.where(CM_condition_grid, 1 - CF_cs_t, 0)
#     CF_cb_t = np.minimum(CF_cb_t, 1)  # Ensure CF_cb_t is between 0 and 1
      
        ## Renewable Conversion 2.0
#     w_cs_t = p_s_t + eta_c * w_o 
#     z_t = CF_e_t * k_e / k_h if k_h != 0 else np.zeros(len(CF_e_t))
    
    
#     # Renewable conversion condition 1.0
#     CM_condition_renewable = eta_c * (p_h + PTC_h - (p_s_t / eta_c) - w_o) > 0 # * z_t
#     renewable_conversion_condition = z_t * k_h <= CF_e_t * k_e
#     CF_cs_t = np.where((CV_h > p_s_t) & CM_condition_renewable & renewable_conversion_condition, z_t, 0)
   
#     # Grid conversion 1.0
#     w_cb_t = p_b_t + eta_c * w_o
#     CM_condition_grid = eta_c * (p_h + PTC_grid - p_b_t / eta_c - w_o) * k_h > 0
#     CF_cb_t = np.where(CM_condition_grid, 1 - CF_cs_t, 0)
    
        
    # CHECK: complementary slackness - unnecessary now as no reconversion but leave in in case
    CF_check_t = CF_cs_t + CF_cb_t 
    CF_check = np.sum(CF_check_t > 1) # if > 0, then CF_c_t and CF_r_t = 1 at the same time!
    
#    print(CF_check_t)
    
    # Average CFs
    CF_cs = np.mean(CF_cs_t)
    CF_cb = np.mean(CF_cb_t)
    
    CF_h_t = CF_cs_t + CF_cb_t
    CF_h_t = np.where(CF_h_t > 0.2, CF_h_t, 0) # New condition for shutting off electrolyzer
    CF_h = np.mean(CF_h_t)
    
    # Calculate hourly carbon intensity of H2
    # CI_h_t = np.where(CF_cs_t + CF_cb_t != 0, CI_full_grid * CF_cb_t / (CF_cs_t + CF_cb_t), 0) # * (1 - CF_cs_t) 
    # changed to actual pre-carbon intensity as official CI_h seems too large. Possibly then no scaling by (CF_cs_t + CF_cb_t) necessary?!
    
    # New simpler approach 
    CI_h_t = np.where(CF_cs_t + CF_cb_t != 0, CI_actual * CF_cb_t, 0) # or take full CI and adjust for renewable load?* (1 - CF_cs_t)
    
    CI_h = np.mean(CI_h_t)

    # Co-variation # Probably not necessary as we already factor in co-variation with electricity prices for our renewables
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

    # Capacity Factors + covariation - leave in for now, important for renewables
    CF_e = np.mean(CF_e_t)
    if CF_e == 0:  # Remove error of dividing by zero
        epsilon_e_t = np.ones(len(CF_e_t)) * 0
    else:
        epsilon_e_t = CF_e_t / CF_e
    p_s = np.mean(p_s_t)
    mu_s_t = p_s_t / p_s
    Gamma_e = np.mean(epsilon_e_t * mu_s_t)
    
    p_b = np.mean(p_b_t)

    # Contribution margins - had to move CM_cb_t up in the if clause for no carbon intensity restriction!
#    LP = np.mean(LP_t)
    CM_cs_t = eta_c * (p_h + PTC_h - p_s_t / eta_c - w_o) * CF_cs_t * k_h
#     CM_cb_t = eta_c * (p_h + PTC_grid_values - p_b_t / eta_c - w_o) * CF_cb_t * k_h
    
    CM_cs = np.mean(CM_cs_t)
    CM_cb = np.mean(CM_cb_t)
    
    # Old
#    CM_cs = (eta_c * (p_h + PTC_h) - w_cs * Gamma_cs) * CF_cs * k_h # added PTC_h also here - excluded Gamma_cs in "w_cs * Gamma_cs)"
#    CM_cb = (eta_c * (p_h + PTC_grid) - w_cb * Gamma_cb) * CF_cb * k_h # Here PTC_grid --- excluded Gamma_cb, as above
    CM_h = CM_cs + CM_cb 
    

    # Added PTC grid vector. Link lifetime vectors from "Variables"
    PTC_h_life = 10
    T_max = 20
    PTC_grid_vector_actual = [0 if year == 0 else (PTC_grid if year <= PTC_h_life else 0) for year in range(T_max + 1)]
    

    # Return summary
    revenues['eta_c'] = eta_c
    revenues['p_s'] = p_s
    revenues['p_b'] = p_b
    revenues['Gamma_e'] = Gamma_e
    revenues['CF_e'] = CF_e
    revenues['ptc'] = ptc
    revenues['k_h'] = k_h
    revenues['p_h'] = p_h
#    revenues['p_h_c'] = p_h_c
#    revenues['delta_h'] = delta_h
#     revenues['w_c'] = w_c
#     revenues['w_c_c'] = w_c * 1/eta_c
    revenues['w_cs'] = w_cs
    revenues['w_cb'] = w_cb
    revenues['CF_h'] = CF_h
#     revenues['CF_c'] = CF_c
    revenues['CF_cs'] = CF_cs
    revenues['CF_cb'] = CF_cb
    revenues['CF_check'] = CF_check
#     revenues['beta_c'] = beta_c
#     revenues['beta_cs'] = beta_cs
#     revenues['beta_cb'] = beta_cb
    revenues['Gamma_cs'] = Gamma_cs
    revenues['Gamma_cb'] = Gamma_cb
    revenues['CM_h'] = CM_h
    revenues['CM_cs'] = CM_cs
    revenues['CM_cb'] = CM_cb
#     revenues['LP'] = LP
    revenues["PTC_h"] = PTC_h
    revenues["PTC_grid"] = PTC_grid
    revenues["PTC_grid_vector_actual"] = PTC_grid_vector_actual
    revenues["CM_cs_t"] = CM_cs_t
    revenues["CM_cb_t"] = CM_cb_t
    revenues["w_o"] = w_o
    revenues["p_s_t"] = p_s_t
    revenues["p_b_t"] = p_b_t
    revenues["CI_h"]= CI_h
    revenues["MSE"]= MSE    
    revenues["CF_e_remaining"]= CF_e_remaining
    
    return revenues

