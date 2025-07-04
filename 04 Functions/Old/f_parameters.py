#!/usr/bin/env python
# coding: utf-8

# In[1]:
import numpy as np

def f_parameters(T, d_i, r, x, alpha, m, ITC_SP_h, alpha_PTC, PTC, t):
    # Initialize a dictionary to store the results
    parameters = {}

    # Discounting
    years = np.arange(0, len(T), 1)  # The end in np.arange is exclusive, so we add 1 - no we have year 0 as new year. delete +1
    lifetime_vect = np.zeros(len(T)) # deleted +1
    lifetime_vect[1:] = 1  # Note: regular discounting from year 1 to T

    gamma = 1/(1+r)
    disc = gamma ** years  # Note: from year 0 to end
    sum_disc = lifetime_vect @ disc  # '@' operator performs matrix multiplication

    # Degradation
    x_vect = (1-x) ** (years[0:])  # Note: 1st entry is year 1
    PV_x = (disc[0:] * x_vect)
    disc_degr = lifetime_vect[0:] @ PV_x

    # Federal Depreciation
    lifetime_vect[0] = 1  # Note: discounting of depreciation from year 0 to T
    PV_depr_f = disc * d_i  # [1xn]
    disc_depr = PV_depr_f @ lifetime_vect

    # Blended Tax Factor
    Delta_e = (1 - alpha * disc_depr) / (1 - alpha)
    Delta_h = (1 - ITC_SP_h - alpha * disc_depr) / (1 - alpha)

    # Levelization factor
    L = m * disc_degr

    # Production Tax Credit
    ptc = np.sum(PTC * disc * x_vect) / ((1 - alpha_PTC) * disc_degr) # gamma ** T, x_i replaced with x_vect  // returns 0.035 $/kwh
    
    # STILL CHANGE PTC levelization !!! 
    
    # Summary
    parameters['gamma'] = gamma
    parameters['disc_degr'] = disc_degr
    parameters['disc_depr'] = disc_depr
    parameters['Delta_e'] = Delta_e
    parameters['Delta_h'] = Delta_h
    parameters['L'] = L
    parameters['ptc'] = ptc
    parameters['disc'] = disc[t] # added

    return parameters

