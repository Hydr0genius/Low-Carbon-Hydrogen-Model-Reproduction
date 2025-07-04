#!/usr/bin/env python
# coding: utf-8

# # Functions

# In[1]:


def lcoe(m, CF, disc_degr, gamma, x_i, T, W_ei, F_ei, SP_e, Delta):
    """This function calculates the LCOE"""
    w_e = np.sum(W_ei * (gamma**T) * x_i) / disc_degr  # [EUR/kWh] levelized variable operating cost
    c_e = SP_e / (m * CF * disc_degr)  # [EUR/kWh] levelized capacity cost
    f_e = np.sum(F_ei * (gamma**T)) / (m * CF * disc_degr)  # [EUR/kWh] levelized fixed operating cost
    LCOE = w_e + f_e + (Delta * c_e)  # [EUR/kWh] levelized cost of electricity incl. discount PTC
    return LCOE, f_e, c_e

