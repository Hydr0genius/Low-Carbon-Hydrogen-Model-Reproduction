def f_overview_vi(parameters, revenues, costs, NPV):
    overview = {}

    overview['k_e'] = NPV['k_e']
    overview['k_h'] = NPV['k_h']

    overview['p_s'] = costs['p_s']
    overview['Gamma_e'] = costs['Gamma_e']

    overview['CF_e'] = costs['CF_e']
    overview['F_e'] = costs['F_e']
    overview['SP_e'] = costs['SP_e']

    overview['f_e'] = costs['f_e']
    overview['c_e'] = costs['c_e']
    overview['Delta_e'] = costs['Delta_e']
    overview['LCOE'] = costs['LCOE']
    overview['LCOE_ptc'] = costs['LCOE_ptc']
    overview['D_e'] = costs['D_e']
    overview['j_r'] = costs['j_r']
    overview['j_c'] = costs['j_c']

    overview['p_h'] = revenues['p_h']
    overview['delta_h'] = revenues['delta_h']
    overview['eta_c'] = costs['eta_c']
    overview['eta_r'] = costs['eta_r']
    overview['eta'] = costs['eta_c'] * costs['eta_r']

    overview['CM_h'] = revenues['CM_h']
    overview['beta_cs'] = revenues['beta_cs']
    overview['beta_cb'] = revenues['beta_cb']
    overview['beta_r'] = revenues['beta_r']
    overview['CF_h'] = revenues['CF_h']
    overview['CF_c'] = revenues['CF_c']
    overview['CF_r'] = revenues['CF_r']

    overview['F_h'] = costs['F_h']
    overview['SP_h'] = costs['SP_h']
    overview['f_h'] = costs['f_h']
    overview['c_h'] = costs['c_h']
    overview['Delta_h'] = parameters['Delta_h']
    overview['LFCH'] = costs['LFCH']

    overview['p_h_c'] = revenues['p_h_c']
    overview['w_c'] = revenues['w_c']
    overview['w_c_c'] = revenues['w_c_c']
    overview['Gamma_cs'] = revenues['Gamma_cs']
    overview['Gamma_cb'] = revenues['Gamma_cb']
    overview['lambda_c'] = revenues['lambda_c']
    overview['f_c'] = costs['f_c']
    overview['c_c'] = costs['c_c']
    overview['Delta_h'] = parameters['Delta_h']
    overview['LPadj_c'] = costs['LPadj_c']
    overview['cap_cost_c'] = ((costs['LCOH'] * costs['eta_c']) - revenues['w_c']) / revenues['lambda_c'] # Careful, lambda_c can be 0!
    overview['LCOH'] = costs['LCOH']

    overview['p_s_PtG'] = revenues['p_s_PtG']
    overview['Gamma_r'] = revenues['Gamma_r']
    overview['w_r'] = revenues['w_r']
    overview['lambda_r'] = revenues['lambda_r']
    overview['f_r'] = costs['f_r']
    overview['c_r'] = costs['c_r']
    overview['Delta_h'] = parameters['Delta_h']
    overview['LPadj_r'] = costs['LPadj_r']
    overview['cap_cost_r'] = (costs['LCOE_h'] - revenues['w_r']) / revenues['lambda_r']
    overview['LCOE_h'] = costs['LCOE_h']

    overview['NPV_h'] = NPV['NPV_h']
    overview['NPV_h2'] = NPV['NPV_h2']
    overview['PM_e'] = NPV['PM_e']
    overview['PM_c'] = NPV['PM_c']
    overview['con'] = NPV['con']
    overview['con_LCOH'] = NPV['con_LCOH']
    overview['PM_r'] = NPV['PM_r']
    overview['recon'] = NPV['recon']
    overview['LCOE_h'] = NPV['LCOE_h']

    return overview
