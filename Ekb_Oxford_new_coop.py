#-*- coding: utf-8 -*-
from numpy import exp, log, fabs as abs, fabs, floor, power as pow

import numpy as np
# def Ekb_Oxford(Y, time):
def Ekb_Oxford(Y, time, b_on, b_off, B_tot, a_on, a_off, k_mu, mu, a_eqmin):
    # if (time - int(time)) <= 0.0001:
    # print time
    dY = np.zeros(48)
    alpha_P_lengthening = 16.0   # per_micrometre (in CE_velocity)
    alpha_P_shortening = 16.0   # per_micrometre (in CE_velocity)
    beta_P_lengthening = 0.0015   # millinewton_second_per_micrometre (in CE_velocity)
    beta_P_shortening = 0.0015   # millinewton_second_per_micrometre (in CE_velocity)
    speed_d = 3.0   # dimensionless (in L_type_Ca_channel_d_gate)
    delta_f = 0.0001   # millivolt (in L_type_Ca_channel_f_gate)
    speed_f = 0.3   # dimensionless (in L_type_Ca_channel_f_gate)
    FrICa = 1.0   # dimensionless (in L_type_Ca_channel)
    Km_f2 = 100000.0   # millimolar (in L_type_Ca_channel)
    Km_f2ds = 0.001   # millimolar (in L_type_Ca_channel)
    P_CaK = 0.002   # dimensionless (in L_type_Ca_channel)
    P_CaNa = 0.01   # dimensionless (in L_type_Ca_channel)
    P_Ca_L = 0.1   # nanoA_per_millimolar (in L_type_Ca_channel)
    R_decay = 20.0   # per_second (in L_type_Ca_channel)
    alpha_S_lengthening = 46.0   # per_micrometre (in PE_velocity)
    alpha_S_shortening = 39.0   # per_micrometre (in PE_velocity)
    beta_S_lengthening = 0.0   # millinewton_second_per_micrometre (in PE_velocity)
    beta_S_shortening = 0.0   # millinewton_second_per_micrometre (in PE_velocity)
    g_bca = 0.00025   # microS (in calcium_background_current)
    K_m_Ca_cyt = 0.0005   # millimolar (in calcium_release)
    K_m_Ca_ds = 0.01   # millimolar (in calcium_release)
    K_m_rel = 10000.0   # per_second (in calcium_release)

    # K_m_rel = 5500.0    # from SVYATOSLAV

    SRLeak = 0.05   # per_second (in calcium_release)
    CaS_tot = 40.0   # millimolar (in calcium_translocation)
    a_tr = 15.0   # per_second (in calcium_translocation)
    alpha_CaS = 50000.0   # per_millimolar_second (in calcium_translocation)
    beta_CaS = 32500.0   # per_second (in calcium_translocation)

    # beta_CaS = 65000.0     #from SVYATOSLAV

    g_1 = 0.6   # per_micrometre (in crossbridge_kinetics)
    g_2 = 0.52   # dimensionless (in crossbridge_kinetics)
    Ca_o = 2.0   # millimolar (in extracellular_calcium_concentration)
    K_b = 4.0   # millimolar (in extracellular_potassium_concentration)
    Na_o = 140.0   # millimolar (in extracellular_sodium_concentration)
    delta_m = 1.0e-5   # millivolt (in fast_sodium_current_m_gate)
    g_Na = 2.5   # microS (in fast_sodium_current)
    E_fibro_stretch = 0.0   # millivolt (in fibroblast)
    c_fibro = 1.0e-5   # microF (in fibroblast)
    g_fibro = 2.0e-4   # microS (in fibroblast)
    g_fibro_junct = 2.9e-4   # microS (in fibroblast)
    g_fibro_stretch = 0.0   # microS (in fibroblast)
    A_tot = 0.07   # millimolar (in intracellular_calcium_concentration)
    B_1_tot = 0.08   # millimolar (in intracellular_calcium_concentration)
    B_2_tot = 0.1   # millimolar (in intracellular_calcium_concentration)
    Kdecay = 10.0   # per_second (in intracellular_calcium_concentration)
    V_ds_ratio = 0.1   # dimensionless (in intracellular_calcium_concentration)
    V_e_ratio = 0.4   # dimensionless (in intracellular_calcium_concentration)
    V_i_ratio = 0.49   # dimensionless (in intracellular_calcium_concentration)
    V_rel_ratio = 0.003   # dimensionless (in intracellular_calcium_concentration)
    V_up_ratio = 0.03   # dimensionless (in intracellular_calcium_concentration)
    # a_off = 200.0   # per_second (in intracellular_calcium_concentration)

    # a_on = 70000.0   # per_millimolar_second (in intracellular_calcium_concentration)
    # a_on = 55000.0  # from SVYATOSLAV

    # a_eqmin = 0.0029
    # a_eqmin = 0.001851
    # tau_inf = 100.
    tau_inf = 1000.

    b_1_off = 182.0   # per_second (in intracellular_calcium_concentration)
    b_1_on = 100000.0   # per_millimolar_second (in intracellular_calcium_concentration)
    b_2_off = 3.0   # per_second (in intracellular_calcium_concentration)
    b_2_on = 1000.0   # per_millimolar_second (in intracellular_calcium_concentration)
    length = 74.0   # micrometre (in intracellular_calcium_concentration)
    pi_min = 0.03   # dimensionless (in intracellular_calcium_concentration)
    radius = 12.0   # micrometre (in intracellular_calcium_concentration)
    n_NaK = 1.5   # dimensionless (in intracellular_sodium_concentration)

    F_afterload = 1.2  # millinewton (in isotonic)
    # F_afterload = 8.  # millinewton (in isotonic)

    l_0 = 0.525139356105856   # micrometre (in length)
    Cm = 9.5e-5   # microF (in membrane)
    F = 96485.3415   # coulomb_per_mole (in membrane)
    R = 8314.472   # joule_per_kilomole_kelvin (in membrane)
    T = 310.0   # kelvin (in membrane)
    stim_amplitude = -3.0   # nanoA (in membrane)
    stim_duration = 0.0025   # second (in membrane)
    stim_end = 10000.0   # second (in membrane)
    stim_period = 1.0   # second (in membrane)
    stim_start = 0.06   # second (in membrane)

    # stim_start = 0.006   # from SVYATOSLAV

    S_0 = 1.14   # micrometre (in parameters_izakov_et_al_1991)
    alpha_G = 1.0   # dimensionless (in parameters_izakov_et_al_1991)
    alpha_Q = 10.0   # dimensionless (in parameters_izakov_et_al_1991)
    beta_Q = 5.0   # dimensionless (in parameters_izakov_et_al_1991)
    k_A = 40.0   # per_millimolar (in parameters_izakov_et_al_1991)
    q_1 = 17.3   # per_second (in parameters_izakov_et_al_1991)
    q_2 = 259.0   # per_second (in parameters_izakov_et_al_1991)
    q_3 = 17.3   # per_second (in parameters_izakov_et_al_1991)
    q_4 = 15.0   # per_second (in parameters_izakov_et_al_1991)
    x_st = 0.964285   # dimensionless (in parameters_izakov_et_al_1991)
    a = 0.25   # dimensionless (in parameters)
    alpha_1 = 14.6   # per_micrometre (in parameters)
    alpha_2 = 14.6   # per_micrometre (in parameters)
    alpha_3 = 48.0   # per_micrometre (in parameters)
    alpha_P = 4.0   # dimensionless (in parameters)
    beta_1 = 0.84   # millinewton (in parameters)
    beta_2 = 0.0018   # millinewton (in parameters)
    beta_3 = 0.015   # millinewton (in parameters)
    chi = 0.705   # dimensionless (in parameters)
    chi_0 = 3.0   # dimensionless (in parameters)
    d_h = 0.5   # dimensionless (in parameters)

    isotonic = 0.0   # dimensionless (in parameters)

    # k_mu = 0.6   # dimensionless (in parameters)
    llambda = 30.0   # millinewton (in parameters)
    m_0 = 0.9   # dimensionless (in parameters)
    # mu = 3.0   # dimensionless (in parameters)
    v_max = 5.5   # micrometre_per_second (in parameters)
    g_pna = 0.004   # microS (in persistent_sodium_current)
    g_Kr1 = 0.0021   # microS (in rapid_delayed_rectifier_potassium_current)
    g_Kr2 = 0.0013   # microS (in rapid_delayed_rectifier_potassium_current)
    P_kna = 0.03   # dimensionless (in reversal_potentials)
    K_cyca = 0.00015   # millimolar (in sarcoplasmic_reticulum_calcium_pump)
    K_inh = 4.0   # millimolar (in sarcoplasmic_reticulum_calcium_pump)
    K_srca = 0.5   # millimolar (in sarcoplasmic_reticulum_calcium_pump)
    K_xcs = 0.4   # dimensionless (in sarcoplasmic_reticulum_calcium_pump)
    alpha_up = 1.0   # millimolar_per_second (in sarcoplasmic_reticulum_calcium_pump)
    beta_up = 0.03   # millimolar_per_second (in sarcoplasmic_reticulum_calcium_pump)
    flag_ingib = 0.0   # dimensionless (in sarcoplasmic_reticulum_calcium_pump)
    g_Ks = 0.0026   # microS (in slow_delayed_rectifier_potassium_current)
    K_kna = 20.0   # millimolar (in sodium_activated_potassium_current)
    g_K_Na = 0.0   # microS (in sodium_activated_potassium_current)
    g_bna = 0.0006   # microS (in sodium_background_current)
    FRiNaCa = 0.001   # dimensionless (in sodium_calcium_exchanger)
    d_NaCa = 0.0   # dimensionless (in sodium_calcium_exchanger)
    gamma1 = 0.5   # dimensionless (gamma in sodium_calcium_exchanger)
    k_NaCa = 0.0005   # nanoA (in sodium_calcium_exchanger)
    n_NaCa = 3.0   # dimensionless (in sodium_calcium_exchanger)
    K_mK = 1.0   # millimolar (in sodium_potassium_pump)
    K_mNa = 24.2   # millimolar (in sodium_potassium_pump)
    i_NaK_max = 0.7   # nanoA (in sodium_potassium_pump)
    K_mk1 = 10.0   # millimolar (in time_indepent_potassium_current)
    g_K1 = 0.5   # microS (in time_indepent_potassium_current)
    g_to = 0.006   # microS (in transient_outward_current)
    g_tos = 0.0   # dimensionless (in transient_outward_current)
    # time = 0.1
    v_st = x_st*v_max
    v_1 = v_max/10.0
    gamma2 = a*d_h*(v_1/v_max)**(2.0)/(3.0*a*d_h-(a+1.0)*v_1/v_max)
    case_1 = a*(0.4+0.4*a)/(v_max*((a+1.0)*0.4)**(2.0))
    case_3 = (0.4*a+1.0)/(a*v_max)
    beta = beta_CaS/alpha_CaS
    V_Cell = np.pi*(radius/1000.0)**2.0*length/1000.0
    V_e = V_Cell*V_e_ratio
    V_i = V_Cell*V_i_ratio
    K_1 = K_cyca*K_xcs/K_srca
    if (Y[0] <= 0.0):
      alp_p = alpha_P_lengthening
    else:
      alp_p = alpha_P_shortening
    if (Y[0] <= 0.0):
      k_P_vis = beta_P_lengthening*exp(alpha_P_lengthening*Y[22])
    else:
      k_P_vis = beta_P_shortening*exp(alpha_P_shortening*Y[22])

    if (Y[0] <= 0.0):
      q_v = q_1-q_2*Y[0]/v_max
    elif ((Y[0] <= v_st) and (0.0 < Y[0])):
      q_v = (q_4-q_3)*Y[0]/v_st+q_3
    else:
      q_v = q_4/(1.0+beta_Q*(Y[0]-v_st)/v_max)**alpha_Q

    if (Y[0] <= 0.0):
      P_star = a*(1.0+Y[0]/v_max)/(a-Y[0]/v_max)
    else:
      P_star = 1.0+d_h-(d_h)**2.0*a/(a*d_h/gamma2*(Y[0]/v_max)**2.0+(a+1.0)*Y[0]/v_max+a*d_h)

    if ((-v_max <= Y[0]) and (Y[0] <= 0.0)):
      G_star = 1.0+0.6*Y[0]/v_max
    elif ((0.0 < Y[0]) and (Y[0] <= v_1)):
      G_star = P_star/((0.4*a+1.0)*Y[0]/(a*v_max)+1.0)
    else:
      G_star = P_star*exp(-alpha_G*((Y[0]-v_1)/v_max)**alpha_P)/((0.4*a+1.0)*Y[0]/(a*v_max)+1.0)

    k_p_v = chi*chi_0*q_v*m_0*G_star
    M_A = (Y[13]/A_tot)**mu*(1.0+(k_mu)**mu)/((Y[13]/A_tot)**mu+(k_mu)**mu)

    if (g_1*Y[22]+g_2 < 0.0):
      n_1 = 0.0
    elif (g_1*Y[22]+g_2 < 1.0):
      n_1 = g_1*Y[22]+g_2
    else:
      n_1 = 1.0

    if (Y[22] > 0.55):
      L_oz = (Y[22]+S_0)/(0.46+S_0)
    else:
      L_oz = (S_0+0.55)*1.0
    k_m_v = chi_0*q_v*(1.0-chi*m_0*G_star)
    K_chi = k_p_v*M_A*n_1*L_oz*(1.0-Y[8])-k_m_v*Y[8]
    p_v = P_star/G_star
    case_2 = a*1.0*(1.0+0.4*a+1.2*Y[0]/v_max+0.6*(Y[0]/v_max)**2.0)/(v_max*((a-Y[0]/v_max)*(1.0+0.6*Y[0]/v_max))**2.0)
    case_4 = 1.0/v_max*exp(-alpha_G*(Y[0]/v_max-v_1/v_max)**alpha_P)*((0.4*a+1.0)/a+alpha_G*alpha_P*(1.0+(0.4*a+1.0)*Y[0]/(a*v_max))*(Y[0]/v_max-v_1/v_max)**(alpha_P-1.0))

    if (Y[0] <= -v_max):
      p_prime_v = case_1
    elif ((-v_max < Y[0]) and (Y[0] <= 0.0)):
      p_prime_v = case_2
    elif ((0.0 < Y[0]) and (Y[0] <= v_1)):
      p_prime_v = case_3
    else:
      p_prime_v = case_4

    F_XSE = beta_3*(exp(alpha_3*Y[24])-1.0)
    F_muscle = F_XSE
    l = Y[23]+Y[24]

    if ((isotonic == 1.0) and (F_muscle > F_afterload) and (l <= l_0*(1.0+1.0e-4))):
      isotonic_mode = 1.0
    else:
      isotonic_mode = 0.0

    if (isotonic_mode == 1.0):
      phi_chi = -(llambda*K_chi*p_v+alp_p*k_P_vis*(Y[0])**(2.0)+alpha_2*beta_2*exp(alpha_2*Y[23])*Y[5])/(llambda*Y[8]*p_prime_v+k_P_vis)
    else:
      phi_chi = -(llambda*K_chi*p_v+alp_p*k_P_vis*(Y[0])**(2.0)+(alpha_2*beta_2*exp(alpha_2*Y[23])+alpha_3*beta_3*exp(alpha_3*Y[24]))*Y[5])/(llambda*Y[8]*p_prime_v+k_P_vis)

    if (isotonic_mode == 1.0):
      phi_chi_2 = alpha_1*beta_1*exp(alpha_1*(Y[23]-Y[22]))*Y[0]/(alpha_1*beta_1*exp(alpha_1*(Y[23]-Y[22]))+alpha_2*beta_2*exp(alpha_2*Y[23]))
    else:
      phi_chi_2 = alpha_1*beta_1*exp(alpha_1*(Y[23]-Y[22]))*Y[0]/(alpha_1*beta_1*exp(alpha_1*(Y[23]-Y[22]))+alpha_2*beta_2*exp(alpha_2*Y[23])+alpha_3*beta_3*exp(alpha_3*Y[24]))

    if (Y[5] <= Y[0]):
      k_S_vis = beta_S_lengthening*exp(alpha_S_lengthening*(Y[23]-Y[22]))
    else:
      k_S_vis = beta_S_shortening*exp(alpha_S_shortening*(Y[23]-Y[22]))

    if (k_S_vis == 0.0):
      dY[0] = (alpha_1*beta_1*exp(alpha_1*(Y[23]-Y[22]))*(phi_chi_2-Y[0])-(llambda*K_chi*p_v+alp_p*k_P_vis*(Y[0])**(2.0)))/(llambda*Y[8]*p_prime_v+k_P_vis)
    else:
      dY[0] = phi_chi

    i_Ca_L_Ca_cyt = (1.0-FrICa)*4.0*P_Ca_L*Y[1]*Y[4]*Y[2]*(Y[25]-50.0)*F/(R*T)/(1.0-exp(-(Y[25]-50.0)*F*2.0/(R*T)))*(Y[17]*exp(100.0*F/(R*T))-Ca_o*exp(-(Y[25]-50.0)*F*2.0/(R*T)))
    i_Ca_L_K_cyt = (1.0-FrICa)*P_CaK*P_Ca_L*Y[1]*Y[4]*Y[2]*(Y[25]-50.0)*F/(R*T)/(1.0-exp(-(Y[25]-50.0)*F/(R*T)))*(Y[20]*exp(50.0*F/(R*T))-Y[9]*exp(-(Y[25]-50.0)*F/(R*T)))
    i_Ca_L_Na_cyt = (1.0-FrICa)*P_CaNa*P_Ca_L*Y[1]*Y[4]*Y[2]*(Y[25]-50.0)*F/(R*T)/(1.0-exp(-(Y[25]-50.0)*F/(R*T)))*(Y[21]*exp(50.0*F/(R*T))-Na_o*exp(-(Y[25]-50.0)*F/(R*T)))
    i_Ca_L_Ca_ds = FrICa*4.0*P_Ca_L*Y[1]*Y[4]*Y[3]*(Y[25]-50.0)*F/(R*T)/(1.0-exp(-(Y[25]-50.0)*F*2.0/(R*T)))*(Y[17]*exp(100.0*F/(R*T))-Ca_o*exp(-(Y[25]-50.0)*F*2.0/(R*T)))
    i_Ca_L_K_ds = FrICa*P_CaK*P_Ca_L*Y[1]*Y[4]*Y[3]*(Y[25]-50.0)*F/(R*T)/(1.0-exp(-(Y[25]-50.0)*F/(R*T)))*(Y[20]*exp(50.0*F/(R*T))-Y[9]*exp(-(Y[25]-50.0)*F/(R*T)))
    i_Ca_L_Na_ds = FrICa*P_CaNa*P_Ca_L*Y[1]*Y[4]*Y[3]*(Y[25]-50.0)*F/(R*T)/(1.0-exp(-(Y[25]-50.0)*F/(R*T)))*(Y[21]*exp(50.0*F/(R*T))-Na_o*exp(-(Y[25]-50.0)*F/(R*T)))
    i_Ca_L = i_Ca_L_Ca_cyt+i_Ca_L_K_cyt+i_Ca_L_Na_cyt+i_Ca_L_Ca_ds+i_Ca_L_K_ds+i_Ca_L_Na_ds
    E0_d = Y[25]+24.0-5.0

    if (abs(E0_d) < 0.0001):
      alpha_d = 120.0
    else:
      alpha_d = 30.0*E0_d/(1.0-exp(-E0_d/4.0))

    if (abs(E0_d) < 0.0001):
      beta_d = 120.0
    else:
      beta_d = 12.0*E0_d/(exp(E0_d/10.0)-1.0)

    dY[1] = speed_d*(alpha_d*(1.0-Y[1])-beta_d*Y[1])
    dY[2] = 1.0-1.0*(Y[17]/(Km_f2+Y[17])+Y[2])
    dY[3] = R_decay*(1.0-(Y[16]/(Km_f2ds+Y[16])+Y[3]))
    E0_f = Y[25]+34.0

    if (abs(E0_f) < delta_f):
      alpha_f = 25.0
    else:
      alpha_f = 6.25*E0_f/(exp(E0_f/4.0)-1.0)

    beta_f = 12.0/(1.0+exp(-1.0*(Y[25]+34.0)/4.0))
    dY[4] = speed_f*(alpha_f*(1.0-Y[4])-beta_f*Y[4])

    if (Y[5] <= Y[0]):
      alp_s = alpha_S_lengthening
    else:
      alp_s = alpha_S_shortening

    # if ((isotonic_mode == 1.0) and (k_S_vis != 0.0)):
    #   dY[5] = (k_S_vis*(phi_chi-alp_s*(Y[5]-Y[0])**(2.0))-alpha_1*beta_1*exp(alpha_1*(Y[23]-Y[22]))*(Y[5]-Y[0])-alpha_2*beta_2*exp(alpha_2*Y[23])*Y[5])/k_S_vis
    # elif ((isotonic_mode == 0.0) and (k_S_vis != 0.0)):
    #   dY[5] = phi_chi-alp_s*(Y[5]-Y[0])**(2.0)-(alpha_1*beta_1*exp(alpha_1*(Y[23]-Y[22]))*(Y[5]-Y[0])+(alpha_2*beta_2*exp(alpha_2*Y[23])+alpha_3*beta_3*exp(alpha_3*Y[24]))*Y[5])/k_S_vis
    # elif (k_S_vis == 0.0):
    #   dY[5] = 0.0

    # Добавил из кода SVYATOSLAV

    if ((isotonic_mode == 1.0) and (abs(k_S_vis) > 0.0001)):
        dY[5] = (k_S_vis*(phi_chi-alp_s*(Y[5]-Y[0])**(2.0))-alpha_1*beta_1*exp(alpha_1*(Y[23]-Y[22]))*(Y[5]-Y[0])-alpha_2*beta_2*exp(alpha_2*Y[23])*Y[5])/k_S_vis
    elif ((isotonic_mode == 0.0) and (abs(k_S_vis) > 0.0001)):
        dY[5] = phi_chi-alp_s*(Y[5]-Y[0])**(2.0)-(alpha_1*beta_1*exp(alpha_1*(Y[23]-Y[22]))*(Y[5]-Y[0])+(alpha_2*beta_2*exp(alpha_2*Y[23])+alpha_3*beta_3*exp(alpha_3*Y[24]))*Y[5])/k_S_vis
    elif (k_S_vis == 0.0):
        dY[5] = 0.0

    # print Y[17]
    E_Ca = 0.5*R*T/F*log(Ca_o/Y[17])
    # i_b_Ca = g_bca*ina.08*(Y[25]-40.0))
    i_b_Ca = g_bca*(Y[25]-E_Ca)
    CaiReg = Y[17]/(Y[17]+K_m_Ca_cyt)
    CadsReg = Y[16]/(Y[16]+K_m_Ca_ds)
    RegBindSite = CaiReg+(1.0-CaiReg)*CadsReg
    ActRate = 500.0*(RegBindSite)**(2.0)
    InactRate = 60.0+500.0*(RegBindSite)**(2.0)

    if (Y[25] < -50.0):
      SpeedRel = 5.0
    else:
      SpeedRel = 1.0

    PrecFrac = 1.0-Y[6]-Y[7]
    dY[6] = PrecFrac*SpeedRel*ActRate-Y[6]*SpeedRel*InactRate
    dY[7] = Y[6]*SpeedRel*InactRate-SpeedRel*1.0*Y[7]
    i_rel = ((Y[6]/(Y[6]+0.25))**(2.0)*K_m_rel+SRLeak)*Y[18]
    i_trans = a_tr*(Y[19]-Y[18])
    dY[8] = K_chi
    E_K = R*T/F*log(Y[9]/Y[20])
    i_Kr = (g_Kr1*Y[26]+g_Kr2*Y[27])*1.0/(1.0+exp((Y[25]+9.0)/22.4))*(Y[25]-E_K)
    E_Ks = R*T/F*log((Y[9]+P_kna*Na_o)/(Y[20]+P_kna*Y[21]))
    i_Ks = g_Ks*(Y[28])**(2.0)*(Y[25]-E_Ks)
    i_K1 = g_K1*Y[9]/(Y[9]+K_mk1)*(Y[25]-E_K)/(1.0+exp((Y[25]-E_K-10.0)*F*1.25/(R*T)))
    i_to = g_to*(g_tos+Y[30]*(1.0-g_tos))*Y[29]*(Y[25]-E_K)
    i_NaK = i_NaK_max*Y[9]/(K_mK+Y[9])*Y[21]/(K_mNa+Y[21])
    dY[9] = 1.0*(i_Kr+i_Ks+i_K1+i_to-1.0/(n_NaK-1.0)*i_NaK+i_Ca_L_K_cyt+i_Ca_L_K_ds)/(1.0*V_e*F)-0.7*(Y[9]-K_b)
    E_mh = R*T/F*log((Na_o+0.12*Y[9])/(Y[21]+0.12*Y[20]))
    i_Na = g_Na*(Y[11])**(3.0)*Y[10]*(Y[25]-E_mh)
    alpha_h = 20.0*exp(-0.125*(Y[25]+75.0))
    beta_h = 2000.0/(1.0+320.0*exp(-0.1*(Y[25]+75.0)))
    dY[10] = alpha_h*(1.0-Y[10])-beta_h*Y[10]
    E0_m = Y[25]+41.0

    if (abs(E0_m) < delta_m):
        alpha_m = 2000.0
    else:
        alpha_m = 200.0*E0_m/(1.0-exp(-0.1*E0_m))

    beta_m = 8000.0*exp(-0.056*(Y[25]+66.0))
    dY[11] = alpha_m*(1.0-Y[11])-beta_m*Y[11]
    i_fibro = g_fibro*(Y[12]+20.0)+g_fibro_stretch*(Y[12]-E_fibro_stretch)
    i_fibro_junct = -g_fibro_junct*(Y[25]-Y[12])
    dY[12] = -(i_fibro+i_fibro_junct)/c_fibro
    F_CE = llambda*p_v*Y[8]
    F_SE = beta_1*(exp(alpha_1*(Y[23]-Y[22]))-1.0)
    F_PE = beta_2*(exp(alpha_2*Y[23])-1.0)

    # N_A = A_tot*Y[8]/(L_oz*Y[13])
    N_A = A_tot*Y[8]/(Y[13])    # L_oz убрали, спросить Сашу

    # if (N_A <= 0.0):
    #   pi_N_A = pi_min
    # elif (N_A <= 1.0):
    #   pi_N_A = (pi_min)**(N_A)
    # else:
    #   pi_N_A = 1.0

    if (N_A < 0.0):
      pi_N_A = 1.
    elif (N_A <= 1.0 and N_A >= 0):
      pi_N_A = (pi_min)**(N_A)
    else:
      pi_N_A = pi_min

    # if Y[8] < 0.0:
    #     pi_N_A = 1.
    # elif A_tot*Y[8] >=0. and A_tot*Y[8]<=1.:
    #     pi_N_A = pi_min**(A_tot*Y[8])
    # elif A_tot*Y[8] >= 1:
    #     pi_N_A = pi_min


    A_off_1 = a_off*pi_N_A*exp(-k_A*Y[13])
    A_off_2 = a_on * a_eqmin
    # A_off_2 = a_on
    A_off = Y[47]*A_off_1 + (1-Y[47])*A_off_2
    # dY[13] = a_on*(A_tot-Y[13])*Y[17]-a_off*exp(-k_A*Y[13])*pi_N_A*Y[13]
    dY[13] = a_on*(A_tot-Y[13])*Y[17] - A_off*Y[13]    # Модифицированная кооперативность


    # dY[14] = b_1_on*(B_1_tot-Y[14])*Y[17]-b_1_off*Y[14]
    # dY[15] = b_2_on*(B_2_tot-Y[15])*Y[17]-b_2_off*Y[15]
    dY[14] = b_on[0] * (B_tot[0] - Y[14]) * Y[17] - b_off[0] * Y[14]
    dY[15] = b_on[1] * (B_tot[1] - Y[15]) * Y[17] - b_off[1] * Y[15]
    for i in xrange(2, 16, 1):
        dY[29 + i] = b_on[i] * (B_tot[i] - Y[29 + i]) * Y[17] - b_off[i] * Y[29 + i]

    i_NaCa_cyt = (1.0-FRiNaCa)*k_NaCa*(exp(gamma1*(n_NaCa-2.0)*Y[25]*F/(R*T))*pow(Y[21],n_NaCa)*Ca_o-exp((gamma1-1.0)*(n_NaCa-2.0)*Y[25]*F/(R*T))*pow(Na_o,n_NaCa)*Y[17])/((1.0+d_NaCa*(Y[17]*pow(Na_o,n_NaCa)+Ca_o*pow(Y[21],n_NaCa)))*(1.0+Y[17]/0.0069))
    K_2 = Y[17]+Y[19]*K_1+K_cyca*K_xcs+K_cyca

    if (flag_ingib == 0.0):
      i_up = Y[17]/K_2*alpha_up-Y[19]*K_1/K_2*beta_up
    else:
      i_up = Y[17]/K_2*alpha_up/(1.0+Y[19]/K_inh)-Y[19]*K_1/K_2*beta_up

    # dY[17] = -1.0/(2.0*1.0*V_i*F)*(i_Ca_L_Ca_cyt+i_b_Ca-2.0/(n_NaCa-2.0)*i_NaCa_cyt)+Y[16]*V_ds_ratio*Kdecay+i_rel*V_rel_ratio/V_i_ratio-dY[13]-dY[14]-dY[15]-i_up
    dY[17] = -1.0 / (2.0 * 1.0 * V_i * F) * (i_Ca_L_Ca_cyt + i_b_Ca - 2.0 / (n_NaCa - 2.0) * i_NaCa_cyt) + \
             Y[16] * V_ds_ratio * Kdecay + i_rel * V_rel_ratio / V_i_ratio - dY[13] - dY[14] - dY[15] - dY[31] - dY[32] \
             - dY[33] - dY[34] - dY[35] - dY[36] - dY[37] - dY[38] - dY[39] - dY[40] - dY[41] - dY[43] - dY[43] - \
             dY[44] - dY[45] - dY[46] - i_up
    i_NaCa_ds = FRiNaCa*k_NaCa*(exp(gamma1*(n_NaCa-2.0)*Y[25]*F/(R*T))*(Y[21])**(n_NaCa)*Ca_o-exp((gamma1-1.0)*(n_NaCa-2.0)*Y[25]*F/(R*T))*(Na_o)**(n_NaCa)*Y[16])/((1.0+d_NaCa*(Y[16]*(Na_o)**(n_NaCa)+Ca_o*(Y[21])**(n_NaCa)))*(1.0+Y[16]/0.0069))
    dY[16] = (-1.0*i_Ca_L_Ca_ds+2.0*i_NaCa_ds/(n_NaCa-2.0))/(2.0*1.0*V_ds_ratio*V_i*F)-Y[16]*Kdecay
    dY[19] = V_i_ratio/V_up_ratio*i_up-i_trans
    dY[18] = (V_up_ratio/V_rel_ratio*i_trans-i_rel)/(1.0+beta*CaS_tot/(Y[18]+beta)**(2.0))
    dY[20] = -1.0/(1.0*V_i*F)*(i_K1+i_Kr+i_Ks+i_Ca_L_K_cyt+i_Ca_L_K_ds+i_to-1.0/(n_NaK-1.0)*i_NaK)
    E_Na = R*T/F*log(Na_o/Y[21])
    i_p_Na = g_pna*1.0/(1.0+exp(-(Y[25]+52.0)/8.0))*(Y[25]-E_Na)
    i_b_Na = g_bna*(Y[25]-E_Na)
    i_NaCa = i_NaCa_cyt+i_NaCa_ds
    dY[21] = -1.0/(1.0*V_i*F)*(i_Na+i_p_Na+i_b_Na+i_Ca_L_Na_cyt+i_Ca_L_Na_ds+n_NaK/(n_NaK-1.0)*i_NaK+n_NaCa/(n_NaCa-2.0)*i_NaCa)
    dl_1_dt = Y[0]
    dY[22] = dl_1_dt

    if (k_S_vis == 0.0):
      dl_2_dt = phi_chi_2
    else:
      dl_2_dt = Y[5]

    dY[23] = dl_2_dt

    if (isotonic_mode == 1.0):
      dl_3_dt = 0.0
    elif ((isotonic_mode == 0.0) and (k_S_vis == 0.0)):
      dl_3_dt = -dl_2_dt
    elif ((isotonic_mode == 0.0) and (k_S_vis != 0.0)):
      dl_3_dt = -Y[5]

    dY[24] = dl_3_dt

    if ((time >= stim_start) and (time <= stim_end) and (time-stim_start-floor((time-stim_start)/stim_period)*stim_period*1.0 <= stim_duration)):
      i_Stim = stim_amplitude
    else:
      i_Stim = 0.0

    dY[25] = -1.0/Cm*(i_Stim+i_K1+i_to+i_Kr+i_Ks+i_NaK+i_Na+i_b_Na+i_p_Na+i_Ca_L_Na_cyt+i_Ca_L_Na_ds+i_NaCa_cyt+i_NaCa_ds+i_Ca_L_Ca_cyt+i_Ca_L_Ca_ds+i_Ca_L_K_cyt+i_Ca_L_K_ds+i_b_Ca)
    alpha_xr1 = 50.0/(1.0+exp(-(Y[25]-5.0)/9.0))
    beta_xr1 = 0.05*exp(-(Y[25]-20.0)/15.0)
    dY[26] = alpha_xr1*(1.0-Y[26])-beta_xr1*Y[26]
    alpha_xr2 = 50.0/(1.0+exp(-(Y[25]-5.0)/9.0))
    beta_xr2 = 0.4*exp(-((Y[25]+30.0)/30.0)**(3.0))
    dY[27] = alpha_xr2*(1.0-Y[27])-beta_xr2*Y[27]
    alpha_xs = 14.0/(1.0+exp(-(Y[25]-40.0)/9.0))
    beta_xs = 1.0*exp(-Y[25]/45.0)
    dY[28] = alpha_xs*(1.0-Y[28])-beta_xs*Y[28]
    i_KNa = g_K_Na*Y[21]/(Y[21]+K_kna)*(Y[25]-E_K)
    dY[29] = 333.0*(1.0/(1.0+exp(-(Y[25]+4.0)/5.0))-Y[29])
    alpha_s = 0.033*exp(-Y[25]/17.0)
    beta_s = 33.0/(1.0+exp(-0.125*(Y[25]+10.0)))
    dY[30] = alpha_s*(1.0-Y[30])-beta_s*Y[30]


    if A_off_1<=A_off_2:
        etha = 1./tau_inf
    else:
        etha = 0.
    dY[47] = - etha * Y[47]
    return dY
