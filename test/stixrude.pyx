import numpy as np
cimport numpy as cnp # cimport gives us access to NumPy's C API

# here we just replicate the function signature from the header
cdef extern from "Forsterite_stixrude_calc.h":
    const char *Forsterite_stixrude_identifier();
    const char *Forsterite_stixrude_name();
    const char *Forsterite_stixrude_formula();
    const double Forsterite_stixrude_mw();
    const double *Forsterite_stixrude_elements();
    double Forsterite_stixrude_g(double t, double p)
    double Forsterite_stixrude_dgdt(double t, double p)
    double Forsterite_stixrude_dgdp(double t, double p)
    double Forsterite_stixrude_d2gdt2(double t, double p)
    double Forsterite_stixrude_d2gdtdp(double t, double p)
    double Forsterite_stixrude_d2gdp2(double t, double p)
    double Forsterite_stixrude_d3gdt3(double t, double p)
    double Forsterite_stixrude_d3gdt2dp(double t, double p)
    double Forsterite_stixrude_d3gdtdp2(double t, double p)
    double Forsterite_stixrude_d3gdp3(double t, double p)
    double Forsterite_stixrude_s(double t, double p)
    double Forsterite_stixrude_v(double t, double p)
    double Forsterite_stixrude_cv(double t, double p)
    double Forsterite_stixrude_cp(double t, double p)
    double Forsterite_stixrude_dcpdt(double t, double p)
    double Forsterite_stixrude_alpha(double t, double p)
    double Forsterite_stixrude_beta(double t, double p)
    double Forsterite_stixrude_K(double t, double p)
    double Forsterite_stixrude_Kp(double t, double p)

# here is the "wrapper" signature
def cy_Forsterite_stixrude_identifier():
    result = <bytes> Forsterite_stixrude_identifier()
    return result.decode('UTF-8')
def cy_Forsterite_stixrude_name():
    result = <bytes> Forsterite_stixrude_name()
    return result.decode('UTF-8')
def cy_Forsterite_stixrude_formula():
    result = <bytes> Forsterite_stixrude_formula()
    return result.decode('UTF-8')
def cy_Forsterite_stixrude_mw():
    result = Forsterite_stixrude_mw()
    return result
def cy_Forsterite_stixrude_elements():
    cdef const double *e = Forsterite_stixrude_elements()
    np_array = np.zeros(106)
    for i in range(0,106):
        np_array[i] = e[i]
    return np_array
def cy_Forsterite_stixrude_g(double t, double p):
    result = Forsterite_stixrude_g(<double> t, <double> p)
    return result
def cy_Forsterite_stixrude_dgdt(double t, double p):
    result = Forsterite_stixrude_dgdt(<double> t, <double> p)
    return result
def cy_Forsterite_stixrude_dgdp(double t, double p):
    result = Forsterite_stixrude_dgdp(<double> t, <double> p)
    return result
def cy_Forsterite_stixrude_d2gdt2(double t, double p):
    result = Forsterite_stixrude_d2gdt2(<double> t, <double> p)
    return result
def cy_Forsterite_stixrude_d2gdtdp(double t, double p):
    result = Forsterite_stixrude_d2gdtdp(<double> t, <double> p)
    return result
def cy_Forsterite_stixrude_d2gdp2(double t, double p):
    result = Forsterite_stixrude_d2gdp2(<double> t, <double> p)
    return result
def cy_Forsterite_stixrude_d3gdt3(double t, double p):
    result = Forsterite_stixrude_d3gdt3(<double> t, <double> p)
    return result
def cy_Forsterite_stixrude_d3gdt2dp(double t, double p):
    result = Forsterite_stixrude_d3gdt2dp(<double> t, <double> p)
    return result
def cy_Forsterite_stixrude_d3gdtdp2(double t, double p):
    result = Forsterite_stixrude_d3gdtdp2(<double> t, <double> p)
    return result
def cy_Forsterite_stixrude_d3gdp3(double t, double p):
    result = Forsterite_stixrude_d3gdp3(<double> t, <double> p)
    return result
def cy_Forsterite_stixrude_s(double t, double p):
    result = Forsterite_stixrude_s(<double> t, <double> p)
    return result
def cy_Forsterite_stixrude_v(double t, double p):
    result = Forsterite_stixrude_v(<double> t, <double> p)
    return result
def cy_Forsterite_stixrude_cv(double t, double p):
    result = Forsterite_stixrude_cv(<double> t, <double> p)
    return result
def cy_Forsterite_stixrude_cp(double t, double p):
    result = Forsterite_stixrude_cp(<double> t, <double> p)
    return result
def cy_Forsterite_stixrude_dcpdt(double t, double p):
    result = Forsterite_stixrude_dcpdt(<double> t, <double> p)
    return result
def cy_Forsterite_stixrude_alpha(double t, double p):
    result = Forsterite_stixrude_alpha(<double> t, <double> p)
    return result
def cy_Forsterite_stixrude_beta(double t, double p):
    result = Forsterite_stixrude_beta(<double> t, <double> p)
    return result
def cy_Forsterite_stixrude_K(double t, double p):
    result = Forsterite_stixrude_K(<double> t, <double> p)
    return result
def cy_Forsterite_stixrude_Kp(double t, double p):
    result = Forsterite_stixrude_Kp(<double> t, <double> p)
    return result
