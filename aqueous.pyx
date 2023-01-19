import numpy as np
cimport numpy as cnp # cimport gives us access to NumPy's C API

cdef extern from "swim.h":
    ctypedef enum SWIM_RegionType:
        NotApplicable,
        DuanAndZhang2006,
        ZhangAndDuan2005,
        HoltenEtAl2014,
        WagnerEtAl2002
    double SWIM_getGibbsFreeEnergy(double t, double p, SWIM_RegionType region);
    double SWIM_getEnthalpy(double t, double p, SWIM_RegionType region);
    double SWIM_getEntropy(double t, double p, SWIM_RegionType region);
    double SWIM_getHeatCapacity(double t, double p, SWIM_RegionType region);
    double SWIM_getDcpDt(double t, double p, SWIM_RegionType region);
    double SWIM_getVolume(double t, double p, SWIM_RegionType region);
    double SWIM_getDvDt(double t, double p, SWIM_RegionType region);
    double SWIM_getDvDp(double t, double p, SWIM_RegionType region);
    double SWIM_getD2vDt2(double t, double p, SWIM_RegionType region);
    double SWIM_getD2vDtDp(double t, double p, SWIM_RegionType region);
    double SWIM_getD2vDp2(double t, double p, SWIM_RegionType region);

# Python instance of C enum variable
eos_type = NotApplicable

# here is the "wrapper" signature
def cy_SWIM_aqueous_identifier():
    result = 'Std-H2O-Int-Model'
    return result
def cy_SWIM_aqueous_calib_identifier():
    return cy_SWIM_aqueous_identifier()

def cy_SWIM_aqueous_name():
    result = 'SWIM'
    return result
def cy_SWIM_aqueous_calib_name():
    return cy_SWIM_aqueous_name()

def cy_SWIM_aqueous_formula():
    result = 'H2O'
    return result
def cy_SWIM_aqueous_calib_formula():
    return cy_SWIM_aqueous_formula()

def cy_SWIM_aqueous_mw():
    result = 18.01528
    return result
def cy_SWIM_aqueous_calib_mw():
    return cy_SWIM_aqueous_mw()

def cy_SWIM_aqueous_elements():
    np_array = np.zeros(106)
    np_array[1] = 2.0
    np_array[8] = 1.0
    return np_array
def cy_SWIM_aqueous_calib_elements():
    return cy_SWIM_aqueous_elements()

def cy_SWIM_aqueous_g(double t, double p):
    result = SWIM_getGibbsFreeEnergy(<double> t, <double> p, eos_type)
    return result
def cy_SWIM_aqueous_calib_g(double t, double p):
    return cy_SWIM_aqueous_g(t, p)

def cy_SWIM_aqueous_dgdt(double t, double p):
    result = -SWIM_getEntropy(<double> t, <double> p, eos_type)
    return result
def cy_SWIM_aqueous_calib_dgdt(double t, double p):
    return cy_SWIM_aqueous_dgdt(t, p)

def cy_SWIM_aqueous_dgdp(double t, double p):
    result = SWIM_getVolume(<double> t, <double> p, eos_type)
    return result
def cy_SWIM_aqueous_calib_dgdp(double t, double p):
    return cy_SWIM_aqueous_dgdp(t, p)

def cy_SWIM_aqueous_d2gdt2(double t, double p):
    result = -SWIM_getHeatCapacity(<double> t, <double> p, eos_type)/t
    return result
def cy_SWIM_aqueous_calib_d2gdt2(double t, double p):
    return cy_SWIM_aqueous_d2gdt2(t, p)

def cy_SWIM_aqueous_d2gdtdp(double t, double p):
    result = SWIM_getDvDt(<double> t, <double> p, eos_type)
    return result
def cy_SWIM_aqueous_calib_d2gdtdp(double t, double p):
    return cy_SWIM_aqueous_d2gdtdp(t, p)

def cy_SWIM_aqueous_d2gdp2(double t, double p):
    result = SWIM_getDvDp(<double> t, <double> p, eos_type)
    return result
def cy_SWIM_aqueous_calib_d2gdp2(double t, double p):
    return cy_SWIM_aqueous_d2gdp2(t, p)

def cy_SWIM_aqueous_d3gdt3(double t, double p):
    result  = -SWIM_getDcpDt(<double> t, <double> p, eos_type)/t
    result += SWIM_getHeatCapacity(<double> t, <double> p, eos_type)/t/t
    return result
def cy_SWIM_aqueous_calib_d3gdt3(double t, double p):
    return cy_SWIM_aqueous_d3gdt3(t, p)

def cy_SWIM_aqueous_d3gdt2dp(double t, double p):
    result = SWIM_getD2vDt2(<double> t, <double> p, eos_type)
    return result
def cy_SWIM_aqueous_calib_d3gdt2dp(double t, double p):
    return cy_SWIM_aqueous_d3gdt2dp(t, p)

def cy_SWIM_aqueous_d3gdtdp2(double t, double p):
    result = SWIM_getD2vDtDp(<double> t, <double> p, eos_type)
    return result
def cy_SWIM_aqueous_calib_d3gdtdp2(double t, double p):
    return cy_SWIM_aqueous_d3gdtdp2(t, p)

def cy_SWIM_aqueous_d3gdp3(double t, double p):
    result = SWIM_getD2vDp2(<double> t, <double> p, eos_type)
    return result
def cy_SWIM_aqueous_calib_d3gdp3(double t, double p):
    return cy_SWIM_aqueous_d3gdp3(t, p)

def cy_SWIM_aqueous_s(double t, double p):
    result = SWIM_getEntropy(<double> t, <double> p, eos_type)
    return result
def cy_SWIM_aqueous_calib_s(double t, double p):
    return cy_SWIM_aqueous_s(t, p)

def cy_SWIM_aqueous_v(double t, double p):
    result = SWIM_getVolume(<double> t, <double> p, eos_type)
    return result
def cy_SWIM_aqueous_calib_v(double t, double p):
    return cy_SWIM_aqueous_v(t, p)

def cy_SWIM_aqueous_cv(double t, double p):
    result = t*cy_SWIM_aqueous_v(t, p)*(cy_SWIM_aqueous_alpha(t, p))**2/cy_SWIM_aqueous_beta(t, p)
    return cy_SWIM_aqueous_cp(t, p) - result
def cy_SWIM_aqueous_calib_cv(double t, double p):
    return cy_SWIM_aqueous_cv(t, p)

def cy_SWIM_aqueous_cp(double t, double p):
    result = SWIM_getHeatCapacity(<double> t, <double> p, eos_type)
    return result
def cy_SWIM_aqueous_calib_cp(double t, double p):
    return cy_SWIM_aqueous_cp(t, p)

def cy_SWIM_aqueous_dcpdt(double t, double p):
    result = SWIM_getDcpDt(<double> t, <double> p, eos_type)
    return result
def cy_SWIM_aqueous_calib_dcpdt(double t, double p):
    return cy_SWIM_aqueous_dcpdt(t, p)

def cy_SWIM_aqueous_alpha(double t, double p):
    result = cy_SWIM_aqueous_d2gdtdp(t, p)/cy_SWIM_aqueous_v(t, p)
    return result
def cy_SWIM_aqueous_calib_alpha(double t, double p):
    return cy_SWIM_aqueous_alpha(t, p)

def cy_SWIM_aqueous_beta(double t, double p):
    result = -cy_SWIM_aqueous_d2gdp2(t, p)/cy_SWIM_aqueous_v(t, p)
    return result
def cy_SWIM_aqueous_calib_beta(double t, double p):
    return cy_SWIM_aqueous_beta(t, p)

def cy_SWIM_aqueous_K(double t, double p):
    result = 1.0/cy_SWIM_aqueous_beta(t, p)
    return result
def cy_SWIM_aqueous_calib_K(double t, double p):
    return cy_SWIM_aqueous_K(t, p)

def cy_SWIM_aqueous_Kp(double t, double p):
    result  = (cy_SWIM_aqueous_d2gdp2(t, p)/cy_SWIM_aqueous_v(t, p))**2
    result -= cy_SWIM_aqueous_d3gdp3(t, p)/cy_SWIM_aqueous_v(t, p)
    return result
def cy_SWIM_aqueous_calib_Kp(double t, double p):
    return cy_SWIM_aqueous_Kp(t, p)

# Methods below only apply to calibration models

def cy_SWIM_aqueous_get_param_number():
    return 1
def cy_SWIM_aqueous_get_param_names():
    result = []
    result.append('EOS/0-auto/1-DZ2006/2-ZD2005/3-Holten/4-Wagner')
    return result
def cy_SWIM_aqueous_get_param_units():
    result = []
    result.append('None')
    return result
def cy_SWIM_aqueous_get_param_values():
    np_array = np.array([eos_type])
    return np_array
def cy_SWIM_aqueous_set_param_values(np_array):
    global eos_type
    n = len(np_array)
    assert n == 1, 'Specify one parameter.'
    assert np_array[0] >= 0 and np_array[0] <5, '0 <= value <= 4'
    eos_type = np_array[0]
    return True
def cy_SWIM_aqueous_get_param_value(int index):
    assert index == 0, 'Only permitted index value is zero.'
    result = eos_type
    return result
def cy_SWIM_aqueous_set_param_value(int index, int value):
    global eos_type
    assert index == 0, 'Only permitted index value is zero.'
    assert isinstance(value, int), 'Value must be an integer'
    assert value >= 0 and value < 5, '0 <= value <= 4'
    eos_type = value
    return True
def cy_SWIM_aqueous_dparam_g(double t, double p, int index):
    return 0
def cy_SWIM_aqueous_dparam_dgdt(double t, double p, int index):
    return 0
def cy_SWIM_aqueous_dparam_dgdp(double t, double p, int index):
    return 0
def cy_SWIM_aqueous_dparam_d2gdt2(double t, double p, int index):
    return 0
def cy_SWIM_aqueous_dparam_d2gdtdp(double t, double p, int index):
    return 0
def cy_SWIM_aqueous_dparam_d2gdp2(double t, double p, int index):
    return 0
def cy_SWIM_aqueous_dparam_d3gdt3(double t, double p, int index):
    return 0
def cy_SWIM_aqueous_dparam_d3gdt2dp(double t, double p, int index):
    return 0
def cy_SWIM_aqueous_dparam_d3gdtdp2(double t, double p, int index):
    return 0
def cy_SWIM_aqueous_dparam_d3gdp3(double t, double p, int index):
    return 0
