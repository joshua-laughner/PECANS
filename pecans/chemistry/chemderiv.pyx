from libc.math cimport exp, sqrt, log, log10

def chem_solver(double dt, double TEMP, double CAIR, double conc_NO, double conc_NO2, double conc_O3, double conc_PAN, double conc_HO, double conc_ACO3, double conc_HNO3):
    cdef double dconc_NO = dNO_dt( TEMP,  CAIR,  conc_NO2,  conc_NO,  conc_O3) * dt
    cdef double dconc_NO2 = dNO2_dt( TEMP,  CAIR,  conc_ACO3,  conc_NO2,  conc_PAN,  conc_O3,  conc_HO,  conc_NO) * dt
    cdef double dconc_O3 = dO3_dt( TEMP,  CAIR,  conc_NO2,  conc_NO,  conc_O3) * dt
    cdef double dconc_PAN = dPAN_dt( TEMP,  CAIR,  conc_ACO3,  conc_NO2,  conc_PAN) * dt
    cdef double dconc_HO = dHO_dt( TEMP,  CAIR,  conc_HO,  conc_NO2) * dt
    cdef double dconc_ACO3 = dACO3_dt( TEMP,  CAIR,  conc_ACO3,  conc_NO2,  conc_PAN) * dt
    cdef double dconc_HNO3 = dHNO3_dt( TEMP,  CAIR,  conc_HO,  conc_NO2) * dt
    return { "conc_NO": conc_NO + dconc_NO, "conc_NO2": conc_NO2 + dconc_NO2, "conc_O3": conc_O3 + dconc_O3, "conc_PAN": conc_PAN + dconc_PAN, "conc_HO": conc_HO + dconc_HO, "conc_ACO3": conc_ACO3 + dconc_ACO3, "conc_HNO3": conc_HNO3 + dconc_HNO3 }

cdef double dNO_dt(double TEMP, double CAIR, double conc_NO2, double conc_NO, double conc_O3):
    cdef double dNO = -1.0*ARR2( 1.40e-12, 1310.0, TEMP )*conc_O3*conc_NO + 1.0*0.01*conc_NO2
    return dNO


cdef double dNO2_dt(double TEMP, double CAIR, double conc_ACO3, double conc_NO2, double conc_PAN, double conc_O3, double conc_HO, double conc_NO):
    cdef double dNO2 = 1.0*ARR2( 1.40e-12, 1310.0, TEMP )*conc_O3*conc_NO + -1.0*0.01*conc_NO2 + -1.0*TROE( 1.49e-30 , 1.8, 2.58e-11 , 0.0, TEMP, CAIR)*conc_HO*conc_NO2 + -1.0*TROE( 9.70e-29 , 5.6 , 9.30e-12 , 1.5 , TEMP, CAIR)*conc_ACO3*conc_NO2 + 1.0*TROEE(1.11e28,14000.0, 9.70e-29 , 5.6 , 9.30e-12 , 1.5 , TEMP, CAIR)*conc_PAN
    return dNO2


cdef double dO3_dt(double TEMP, double CAIR, double conc_NO2, double conc_NO, double conc_O3):
    cdef double dO3 = -1.0*ARR2( 1.40e-12, 1310.0, TEMP )*conc_O3*conc_NO + 1.0*0.01*conc_NO2
    return dO3


cdef double dPAN_dt(double TEMP, double CAIR, double conc_ACO3, double conc_NO2, double conc_PAN):
    cdef double dPAN = 1.0*TROE( 9.70e-29 , 5.6 , 9.30e-12 , 1.5 , TEMP, CAIR)*conc_ACO3*conc_NO2 + -1.0*TROEE(1.11e28,14000.0, 9.70e-29 , 5.6 , 9.30e-12 , 1.5 , TEMP, CAIR)*conc_PAN
    return dPAN


cdef double dHO_dt(double TEMP, double CAIR, double conc_HO, double conc_NO2):
    cdef double dHO = -1.0*TROE( 1.49e-30 , 1.8, 2.58e-11 , 0.0, TEMP, CAIR)*conc_HO*conc_NO2
    return dHO


cdef double dACO3_dt(double TEMP, double CAIR, double conc_ACO3, double conc_NO2, double conc_PAN):
    cdef double dACO3 = -1.0*TROE( 9.70e-29 , 5.6 , 9.30e-12 , 1.5 , TEMP, CAIR)*conc_ACO3*conc_NO2 + 1.0*TROEE(1.11e28,14000.0, 9.70e-29 , 5.6 , 9.30e-12 , 1.5 , TEMP, CAIR)*conc_PAN
    return dACO3


cdef double dHNO3_dt(double TEMP, double CAIR, double conc_HO, double conc_NO2):
    cdef double dHNO3 = 1.0*TROE( 1.49e-30 , 1.8, 2.58e-11 , 0.0, TEMP, CAIR)*conc_HO*conc_NO2
    return dHNO3


####################
# RATE EXPRESSIONS #
####################

cdef ARR2(double A0, double B0, double TEMP):
    return A0 * exp(-B0 / TEMP)


cdef TROE(double k0_300K, double n, double kinf_300K, double m, double TEMP, double CAIR):
    cdef double zt_help
    cdef double k0_T
    cdef double kinf_T
    cdef double k_ratio
    zt_help = 300.0 / TEMP;
    k0_T    = k0_300K   * zt_help ** n * CAIR   # k_0   at current T
    kinf_T  = kinf_300K * zt_help ** m          # k_inf at current T
    k_ratio = k0_T/kinf_T
    return k0_T/(1.0 + k_ratio)*0.6 ** (1.0 / (1.0+log10(k_ratio)**2))


cdef TROEE(double A, double B, double k0_300K, double n, double kinf_300K, double m, double TEMP, double CAIR):
    cdef double zt_help
    cdef double k0_T
    cdef double kinf_T
    cdef double k_ratio
    cdef double troe
    zt_help = 300.0 / TEMP;
    k0_T    = k0_300K   * zt_help ** n * CAIR   # k_0   at current T
    kinf_T  = kinf_300K * zt_help ** m          # k_inf at current T
    k_ratio = k0_T/kinf_T
    troe = k0_T/(1.0 + k_ratio)*0.6 ** (1.0 / (1.0+log10(k_ratio)**2))
    return A * exp(- B / TEMP) * troe


