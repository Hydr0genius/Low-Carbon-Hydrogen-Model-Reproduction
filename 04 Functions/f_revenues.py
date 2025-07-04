#!/usr/bin/env python
# coding: utf-8

# In[2]:


def f_revenues(j, p_h, eta_c, eta_r, p_s_t, p_b_t, w_o, delta_h, LP_t, phi_t):

    # Hourly capacity factors
    w_c_t = p_b_t + eta_c * w_o
    p_b = np.mean(p_b_t)
    p_h_c = eta_c * p_h
    CF_c_t = np.zeros(len(w_c_t)) if j > 1 else w_c_t < p_h_c

    w_r = (1/eta_r) * (p_h + delta_h)
    CF_r_t = p_s_t > w_r

    # CHECK: complementary slackness
    CF_check_t = CF_c_t + CF_r_t
    CF_check = np.sum(CF_check_t > 1)

    # Average CFs
    CF_c = np.mean(CF_c_t)
    CF_r = np.mean(CF_r_t)
    CF_h = np.mean(CF_check_t)

    beta_c = CF_c / CF_h
    beta_r = CF_r / CF_h

    # Co-variation
    epsilon_c_t = np.zeros(len(CF_c_t)) if CF_c == 0 else CF_c_t / CF_c
    w_c = np.mean(w_c_t)
    mu_c_t = w_c_t / w_c
    Gamma_c = np.mean(epsilon_c_t * mu_c_t)

    epsilon_r_t = np.zeros(len(CF_r_t)) if CF_r == 0 else CF_r_t / CF_r
    p_s = np.mean(p_s_t)
    mu_s_t = p_s_t / p_s
    Gamma_r = np.mean(epsilon_r_t * mu_s_t)

    # Contribution margins
    LP = np.mean(LP_t)
    CM_c = (eta_c * p_h - w_c * Gamma_c) * CF_c
    CM_r = (p_s * Gamma_r - w_r) * CF_r
    CM_h = CM_c + CM_r

    # Cost allocation
    lambda_c = CM_c / CM_h
    lambda_r = CM_r / CM_h

    # Return summary
    revenues = {
        'eta_c': eta_c,
        'eta_r': eta_r,
        'p_h': p_h,
        'p_h_c': p_h_c,
        'w_c': w_c,
        'p_s': p_s,
        'p_b': p_b,
        'w_r': w_r,
        'CF_h': CF_h,
        'CF_c': CF_c,
        'CF_r': CF_r,
        'CF_check': CF_check,
        'beta_c': beta_c,
        'beta_r': beta_r,
        'Gamma_c': Gamma_c,
        'Gamma_r': Gamma_r,
        'CM_h': CM_h,
        'CM_c': CM_c,
        'CM_r': CM_r,
        'lambda_c': lambda_c,
        'lambda_r': lambda_r,
        'LP': LP,
    }
    
    return revenues

