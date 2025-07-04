#!/usr/bin/env python
# coding: utf-8

# In[1]:
import numpy as np

def f_parameters(T, d_i, r, x, alpha_f, alpha_s, m, ITC_SP_e, ITC_SP_h, PTC_vector, PTC_h_vector, PTC_grid_vector, t):
    # Initialize a dictionary to store the results
    parameters = {}
    
    # Discounting
    years = np.arange(0, len(T), 1)  # The end in np.arange is exclusive - add 1
    lifetime_vect = np.zeros(len(T)) # deleted +1
    lifetime_vect[1:] = 1  # Note: regular discounting from year 1 to T

    gamma = 1/(1+r)
    disc = gamma ** years  # Note: from year 0 to end
   
    sum_disc = lifetime_vect @ disc  # '@' operator performs matrix multiplication - result in scalar 

    # Degradation
    x_vect = (1-x) ** (years[0:]-1) # Note: 1st entry is year 1 // Had to manually adjust year 0 value to 0, then start with 1 in year 1
    x_vect[0] = 1 # Was 0 before but then it messes up the sumation of PTCs - no simply move back PTC payments by 1 year!
    PV_x = (disc[0:] * x_vect)
    disc_degr = lifetime_vect[0:] @ PV_x
 

    # Federal Depreciation
    lifetime_vect[0] = 1  # Note: discounting of depreciation from year 0 to T, necessary?
    PV_depr_f = disc * d_i  # [1xn]

    disc_depr = PV_depr_f @ lifetime_vect # Should result in scalar - yes


    # Tax Factor:  previously
#     Delta_e = (1 - alpha * disc_depr) / (1 - alpha)
#     Delta_h = (1 - ITC_SP_h - alpha * disc_depr) / (1 - alpha)
    
    small_delta_e = 0.848 # ITC Capitalization factor - Currently from excel file (maybe calculate in here)
    small_delta_h = 0.848 
    
    # Federal Tax Factor
    Delta_e = (1 - ITC_SP_e - alpha_f * disc_depr * (1 - small_delta_e * ITC_SP_e)) / (1 - alpha_f)
    Delta_h = (1 - ITC_SP_h - alpha_f * disc_depr * (1 - small_delta_h * ITC_SP_h)) / (1 - alpha_f)
    
#     # Blended Tax Factor: New - currently not implemented with different depreciation schedules cause effortful...
#     Delta_e_s = (1 - ITC_SP_e - alpha_s * disc_depr - alpha_f * disc_depr * (1 - ITC_SP_e)) / (1 - alpha_s - alpha_f * (1 - alpha_s))
#     Delta_h_s = (1 - ITC_SP_h - alpha_s * disc_depr - alpha_f * disc_depr * (1 - ITC_SP_h)) / (1 - alpha_s - alpha_f * (1 - alpha_s))

    # Levelization factor
    L = m * disc_degr

    # Production Tax Credit - now different: simply calc the vector term and add in the LROE and LROH part
    
    # Old ptc:
    # ptc = np.sum(PTC * disc * x_vect) / ((1 - alpha_PTC) * disc_degr) # gamma ** T, x_i replaced with x_vect  // returns 0.035 $/kwh
    
    # PTC sumation
    ptc_sum = np.sum(PTC_vector * disc * x_vect)
    ptc_h_sum = np.sum(PTC_h_vector * disc * x_vect) # Still TEST
    PTC_grid_sum = np.sum(PTC_grid_vector * disc * x_vect)
    
    # Also export the current PTC in the given year t (as imported in the other functions)
    # STILL CHANGE PTC levelization !!!  
    
    # Summary
    parameters['gamma'] = gamma
    parameters['m'] = m
    parameters['alpha'] = alpha_f
    parameters['disc_degr'] = disc_degr
    parameters['disc_depr'] = disc_depr
    parameters['disc'] = disc
    parameters['x_vect'] = x_vect
    parameters['Delta_e'] = Delta_e
    parameters['Delta_h'] = Delta_h
    parameters['L'] = L
    parameters['PTC_vector'] = PTC_vector
    parameters['PTC_h_vector'] = PTC_h_vector
    parameters['PTC_grid_vector'] = PTC_grid_vector
    parameters['ptc'] = PTC_vector[t]
    parameters['ptc_h'] = PTC_h_vector[t]
    parameters['ptc_grid'] = PTC_grid_vector[t] # Still introduce this yearly factor for revenue function
    parameters['ptc_sum'] = ptc_sum # Still import to NPV
    parameters['ptc_h_sum'] = ptc_h_sum # Still import to NPV
    parameters['PTC_grid_sum'] = PTC_grid_sum# Still import to NPV
#    parameters['disc'] = disc # added

    return parameters


