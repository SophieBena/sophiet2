'''
The Coder Templates module provides templates that suport code generation in coder.py classes.

'''
__all__ = ['create_calib_c_template',
           'create_calib_h_template',
           'create_calib_pyx_template',
           'create_code_for_born_functions',
           'create_code_for_debye_function',
           'create_code_for_dh_functions',
           'create_fast_c_template',
           'create_fast_h_template',
           'create_fast_pyx_template',
           'create_pyxbld_template',
           'create_redundant_calib_TV_template',
           'create_redundant_function_template',
           'create_soln_calc_template',
           'create_soln_calib_code_template',
           'create_soln_calib_extra_template',
           'create_soln_calib_include_template',
           'create_soln_calib_pyx_template',
           'create_soln_calib_template',
           'create_soln_deriv_template',
           'create_soln_fast_code_template',
           'create_soln_fast_include_template',
           'create_soln_fast_pyx_template',
           'create_soln_pyxbld_template',
           'create_soln_redun_template',
           'create_soln_std_state_include_template',
           'create_ordering_gaussj_template',
           'create_ordering_code_template',
           'create_complx_soln_calc_template',
           'create_complx_soln_calib_template',
           'create_speciation_code_template',
           'create_speciation_ordering_code_template'
           ]

##############################
# Generic External Functions #
##############################

def create_calib_c_template(language='C'):
    """
    Retrieves template for a C function file that implements calibration
    mode calculation of model functions for a specific phase instance.

    The calibration functions expose the model parameters using getters and
    setters, which otherwise are fixed as constants in the *fast* routines,
    allowing the compiler to optimize the resulting code.

    The user does not normally call this function directly.

    Parameters
    ----------
    language : string
        Language syntax for generated code, ("C" is the C99 programming
        language)

    Returns
    -------
    string :
        The template string.
    """

    if language == 'C':
        return _create_calib_c_template_c()
    elif language == 'C++':
        return _create_calib_c_template_cpp()
    else:
        raise NotImplementedError('Language not implemented.')
        return ""

def create_calib_h_template(language='C'):
    """
    Retrieves template for an include file that implements calibration mode
    calculation of model functions for a specific phase instance.

    The calibration functions expose the model parameters using getters and
    setters, which otherwise are fixed as constants in the *fast* routines,
    allowing the compiler to optimize the resulting code.

    The user does not normally call this function directly.

    Parameters
    ----------
    language : string
        Language syntax for generated code, ("C" is the C99 programming
        language)

    Returns
    -------
    string :
        The template string.
    """

    if language == 'C':
        return _create_calib_h_template_c()
    elif language == 'C++':
        return _create_calib_h_template_cpp()
    else:
        raise NotImplementedError('Language not implemented.')
        return ""

def create_calib_pyx_template(language='C'):
    """
    Retrieves calib code template for cython pyx file.

    Parameters
    ----------
    language : string
        Language syntax for generated code, ("C" is the C99 programming
        language)

    Returns
    -------
    string :
        The template string.
    """

    if language == 'C':
        return _create_calib_pyx_template_c()
    elif language == 'C++':
        return _create_calib_pyx_template_cpp()
    else:
        raise NotImplementedError('Language not implemented.')
        return ""

def create_code_for_born_functions(language='C'):
    """
    Retrieves code that provides a reference to the Born functions.

    Parameters
    ----------
    language : string
        Language syntax for generated code, ("C" is the C99 programming
        language)

    Returns
    -------
    string :
        The template string.
    """

    if language == 'C':
        return _create_code_for_born_functions_c()
    elif language == 'C++':
        return _create_code_for_born_functions_cpp()
    else:
        raise NotImplementedError('Language not implemented.')
        return ""

def create_code_for_debye_function(language='C'):
    """
    Retrieves a block of code that provides an implementation of the Debye
    function.

    Parameters
    ----------
    language : string
        Language syntax for generated code, ("C" is the C99 programming
        language)

    Returns
    -------
    string :
        The template string.
    """

    if language == 'C':
        return _create_code_for_debye_function_c()
    elif language == 'C++':
        return _create_code_for_debye_function_cpp()
    else:
        raise NotImplementedError('Language not implemented.')
        return ""

def create_code_for_dh_functions(language='C'):
    """
    Retrieves code that provides a reference to the Debye-Hückel
    solvent functions.

    Parameters
    ----------
    language : string
        Language syntax for generated code, ("C" is the C99 programming
        language)

    Returns
    -------
    string :
        The template string.
    """

    if language == 'C':
        return _create_code_for_dh_functions_c()
    elif language == 'C++':
        return _create_code_for_dh_functions_cpp()
    else:
        raise NotImplementedError('Language not implemented.')
        return ""

def create_fast_c_template(language='C'):
    """
    Retrieves template for a C function file that implements fast
    calculation of model functions for a specific phase instance.

    The user does not normally call this function directly.

    Parameters
    ----------
    language : string
        Language syntax for generated code, ("C" is the C99 programming
        language)

    Returns
    -------
    string :
        The template string.
    """

    if language == 'C':
        return _create_fast_c_template_c()
    elif language == 'C++':
        return _create_fast_c_template_cpp()
    else:
        raise NotImplementedError('Language not implemented.')
        return ""

def create_fast_h_template(language='C'):
    """
    Retrieves template for an include file that implements fast calculation
    of model functions for a specific phase instance.

    The user does not normally call this function directly.

    Parameters
    ----------
    language : string
        Language syntax for generated code, ("C" is the C99 programming
        language)

    Returns
    -------
    string :
        The template string.
    """

    if language == 'C':
        return _create_fast_h_template_c()
    elif language == 'C++':
        return _create_fast_h_template_cpp()
    else:
        raise NotImplementedError('Language not implemented.')
        return ""

def create_fast_pyx_template(language='C'):
    """
    Retrieves fast code template for cython pyx file.

    Parameters
    ----------
    language : string
        Language syntax for generated code, ("C" is the C99 programming language)

    Returns
    -------
    string :
        The template string.
    """

    if language == 'C':
        return _create_fast_pyx_template_c()
    elif language == 'C++':
        return _create_fast_pyx_template_cpp()
    else:
        raise NotImplementedError('Language not implemented.')
        return ""

def create_pyxbld_template(language='C'):
    """
    Retrieves template for cython pyxbld file.

    Parameters
    ----------
    language : string
        Language syntax for generated code, ("C" is the C99 programming
        language)

    Returns
    -------
    string :
        The template string.
    """

    if language == 'C':
        return _create_pyxbld_template_c()
    elif language == 'C++':
        return _create_pyxbld_template_cpp()
    else:
        raise NotImplementedError('Language not implemented.')
        return ""

def create_redundant_calib_TV_template(language='C'):
    """
    Retrieves template for redundant thermodynamic functions associated with
    Helmholtz models.

    Parameters
    ----------
    language : string
        Language syntax for generated code, ("C" is the C99 programming
        language)

    Returns
    -------
    string :
        The template string.
    """

    if language == 'C':
        return _create_redundant_calib_TV_template_c()
    elif language == 'C++':
        return _create_redundant_calib_TV_template_cpp()
    else:
        raise NotImplementedError('Language not implemented.')
        return ""

def create_redundant_function_template(language='C', model_type='TP'):
    """
    Retrieves template for redundant thermodynamic functions.

    Parameters
    ----------
    language : string
        Language syntax for generated code, ("C" is the C99 programming
        language)
    model_type: string
        Potential type, either Gibbs free energy ('TP') or Helmholtz free
        energy ('TV')

    Returns
    -------
    string :
        The template string.
    """

    if language == 'C':
        return _create_redundant_function_template_c(model_type=model_type)
    elif language == 'C++':
        return _create_redundant_function_template_cpp(model_type=model_type)
    else:
        raise NotImplementedError('Language not implemented.')
        return ""

def create_soln_calc_template(language='C'):
    """
    Retrieves template for solution functions.

    Parameters
    ----------
    language : string
        Language syntax for generated code, ("C" is the C99 programming
        language)

    Returns
    -------
    string :
        The template string.
    """

    if language == 'C':
        return _create_soln_calc_template_c()
    elif language == 'C++':
        return _create_soln_calc_template_cpp()
    else:
        raise NotImplementedError('Language not implemented.')
        return ""

def create_soln_calib_code_template(language='C'):
    """
    Retrieves template for calibration code solution template.

    Parameters
    ----------
    language : string
        Language syntax for generated code, ("C" is the C99 programming
        language)

    Returns
    -------
    string :
        The template string.
    """

    if language == 'C':
        return _create_soln_calib_code_template_c()
    elif language == 'C++':
        return _create_soln_calib_code_template_cpp()
    else:
        raise NotImplementedError('Language not implemented.')
        return ""

def create_soln_calib_extra_template(language='C'):
    """
    Retrieves template for solution calibration functions that retrieve and set
    values of model parameters.

    Parameters
    ----------
    language : string
        Language syntax for generated code, ("C" is the C99 programming
        language)

    Returns
    -------
    string :
        The template string.
    """

    if language == 'C':
        return _create_soln_calib_extra_template_c()
    elif language == 'C++':
        return _create_soln_calib_extra_template_cpp()
    else:
        raise NotImplementedError('Language not implemented.')
        return ""

def create_soln_calib_include_template(language='C'):
    """
    Retrieves template for calibration include file solution template.

    Parameters
    ----------
    language : string
        Language syntax for generated code, ("C" is the C99 programming
        language)

    Returns
    -------
    string :
        The template string.
    """

    if language == 'C':
        return _create_soln_calib_include_template_c()
    elif language == 'C++':
        return _create_soln_calib_include_template_cpp()
    else:
        raise NotImplementedError('Language not implemented.')
        return ""

def create_soln_calib_pyx_template(language='C'):
    """
    Retrieves calib solution code template for cython pyx file.

    Parameters
    ----------
    language : string
        Language syntax for generated code, ("C" is the C99 programming language)

    Returns
    -------
    string :
        The template string.
    """

    if language == 'C':
        return _create_soln_calib_pyx_template_c()
    elif language == 'C++':
        return _create_soln_calib_pyx_template_cpp()
    else:
        raise NotImplementedError('Language not implemented.')
        return ""

def create_soln_calib_template(language='C'):
    """
    Retrieves template for solution calibration functions.

    Parameters
    ----------
    language : string
        Language syntax for generated code, ("C" is the C99 programming
        language)

    Returns
    -------
    string :
        The template string.
    """

    if language == 'C':
        return _create_soln_calib_template_c()
    elif language == 'C++':
        return _create_soln_calib_template_cpp()
    else:
        raise NotImplementedError('Language not implemented.')
        return ""

def create_soln_deriv_template(language='C'):
    """
    Retrieves template for solution derivative functions.

    Parameters
    ----------
    language : string
        Language syntax for generated code, ("C" is the C99 programming
        language)

    Returns
    -------
    string :
        The template string.
    """

    if language == 'C':
        return _create_soln_deriv_template_c()
    elif language == 'C++':
        return _create_soln_deriv_template_cpp()
    else:
        raise NotImplementedError('Language not implemented.')
        return ""

def create_soln_fast_code_template(language='C'):
    """
    Retrieves template for fast code solution template.

    Parameters
    ----------
    language : string
        Language syntax for generated code, ("C" is the C99 programming
        language)

    Returns
    -------
    string :
        The template string.
    """

    if language == 'C':
        return _create_soln_fast_code_template_c()
    elif language == 'C++':
        return _create_soln_fast_code_template_cpp()
    else:
        raise NotImplementedError('Language not implemented.')
        return ""

def create_soln_fast_include_template(language='C'):
    """
    Retrieves template for fast include file solution template.

    Parameters
    ----------
    language : string
        Language syntax for generated code, ("C" is the C99 programming
        language)

    Returns
    -------
    string :
        The template string.
    """

    if language == 'C':
        return _create_soln_fast_include_template_c()
    elif language == 'C++':
        return _create_soln_fast_include_template_cpp()
    else:
        raise NotImplementedError('Language not implemented.')
        return ""

def create_soln_fast_pyx_template(language='C'):
    """
    Retrieves fast solution code template for cython pyx file.

    Parameters
    ----------
    language : string
        Language syntax for generated code, ("C" is the C99 programming language)

    Returns
    -------
    string :
        The template string.
    """

    if language == 'C':
        return _create_soln_fast_pyx_template_c()
    elif language == 'C++':
        return _create_soln_fast_pyx_template_cpp()
    else:
        raise NotImplementedError('Language not implemented.')
        return ""

def create_soln_pyxbld_template(language='C'):
    """
    Retrieves template for cython solution pyxbld file.

    Parameters
    ----------
    language : string
        Language syntax for generated code, ("C" is the C99 programming
        language)

    Returns
    -------
    string :
        The template string.
    """

    if language == 'C':
        return _create_soln_pyxbld_template_c()
    elif language == 'C++':
        return _create_soln_pyxbld_template_cpp()
    else:
        raise NotImplementedError('Language not implemented.')
        return ""

def create_soln_redun_template(language='C'):
    """
    Retrieves template for solution redundant functions.

    Parameters
    ----------
    language : string
        Language syntax for generated code, ("C" is the C99 programming
        language)

    Returns
    -------
    string :
        The template string.
    """

    if language == 'C':
        return _create_soln_redun_template_c()
    elif language == 'C++':
        return _create_soln_redun_template_cpp()
    else:
        raise NotImplementedError('Language not implemented.')
        return ""

def create_soln_std_state_include_template(language='C'):
    """
    Retrieves template for standard state properties include template.

    Parameters
    ----------
    language : string
        Language syntax for generated code, ("C" is the C99 programming
        language)

    Returns
    -------
    string :
        The template string.
    """

    if language == 'C':
        return _create_soln_std_state_include_template_c()
    elif language == 'C++':
        return _create_soln_std_state_include_template_cpp()
    else:
        raise NotImplementedError('Language not implemented.')
        return ""

def create_ordering_gaussj_template(language='C'):
    """
    Retrieves template for generation of guassj code.

    Parameters
    ----------
    language : string
        Language syntax for generated code, ("C" is the C99 programming
        language)

    Returns
    -------
    string :
        The template string.
    """

    if language == 'C':
        return _create_ordering_gaussj_template_c()
    elif language == 'C++':
        return _create_ordering_gaussj_template_cpp()
    else:
        raise NotImplementedError('Language not implemented.')
        return ""

def create_ordering_code_template(language='C'):
    """
    Retrieves template for generation of ordering code.

    Parameters
    ----------
    language : string
        Language syntax for generated code, ("C" is the C99 programming
        language)

    Returns
    -------
    string :
        The template string.
    """

    if language == 'C':
        return _create_ordering_code_template_c()
    elif language == 'C++':
        return _create_ordering_code_template_cpp()
    else:
        raise NotImplementedError('Language not implemented.')
        return ""

def create_complx_soln_calc_template(language='C'):
    """
    Retrieves template for generation of solution property code.

    Parameters
    ----------
    language : string
        Language syntax for generated code, ("C" is the C99 programming
        language)

    Returns
    -------
    string :
        The template string.
    """

    if language == 'C':
        return _create_complx_soln_calc_template_c()
    elif language == 'C++':
        return _create_complx_soln_calc_template_cpp()
    else:
        raise NotImplementedError('Language not implemented.')
        return ""

def create_complx_soln_calib_template(language='C'):
    """
    Retrieves template for generation of solution property code.

    Parameters
    ----------
    language : string
        Language syntax for generated code, ("C" is the C99 programming
        language)

    Returns
    -------
    string :
        The template string.
    """

    if language == 'C':
        return _create_complx_soln_calib_template_c()
    elif language == 'C++':
        return _create_complx_soln_calib_template_cpp()
    else:
        raise NotImplementedError('Language not implemented.')
        return ""

def create_speciation_code_template(language='C'):
    """
    Retrieves template for generation of speciation code.

    Parameters
    ----------
    language : string
        Language syntax for generated code, ("C" is the C99 programming
        language)

    Returns
    -------
    string :
        The template string.
    """

    if language == 'C':
        return _create_speciation_code_template_c()
    elif language == 'C++':
        return _create_speciation_code_template_cpp()
    else:
        raise NotImplementedError('Language not implemented.')
        return ""

def create_speciation_ordering_code_template(language='C'):
    """
    Retrieves template for generation of speciation ordering code.

    Parameters
    ----------
    language : string
        Language syntax for generated code, ("C" is the C99 programming
        language)

    Returns
    -------
    string :
        The template string.
    """

    if language == 'C':
        return _create_speciation_ordering_code_template_c()
    elif language == 'C++':
        return _create_speciation_ordering_code_template_cpp()
    else:
        raise NotImplementedError('Language not implemented.')
        return ""

##############################
# Standard State C templates #
##############################

def _create_calib_c_template_c():
    """
    C language implementation of create_calib_c_template()

    Returns
    -------
    string :
        The template string.
    """

    return """\

static char *identifier = "{git_identifier}";
{parameter_init_block}

{include_calc_h}
{include_calib_h}

const char *{phase}_{module}_calib_identifier(void) {{
    return identifier;
}}

const char *{phase}_{module}_calib_name(void) {{
    return "{phase}";
}}

const char *{phase}_{module}_calib_formula(void) {{
    return "{formula}";
}}

const double {phase}_{module}_calib_mw(void) {{
    return {mw};
}}

static const double elmformula[106] = {{
        {elmvector[0]},{elmvector[1]},{elmvector[2]},{elmvector[3]},{elmvector[4]},{elmvector[5]},
        {elmvector[6]},{elmvector[7]},{elmvector[8]},{elmvector[9]},{elmvector[10]},{elmvector[11]},
        {elmvector[12]},{elmvector[13]},{elmvector[14]},{elmvector[15]},{elmvector[16]},{elmvector[17]},
        {elmvector[18]},{elmvector[19]},{elmvector[20]},{elmvector[21]},{elmvector[22]},{elmvector[23]},
        {elmvector[24]},{elmvector[25]},{elmvector[26]},{elmvector[27]},{elmvector[28]},{elmvector[29]},
        {elmvector[30]},{elmvector[31]},{elmvector[32]},{elmvector[33]},{elmvector[34]},{elmvector[35]},
        {elmvector[36]},{elmvector[37]},{elmvector[38]},{elmvector[39]},{elmvector[40]},{elmvector[41]},
        {elmvector[42]},{elmvector[43]},{elmvector[44]},{elmvector[45]},{elmvector[46]},{elmvector[47]},
        {elmvector[48]},{elmvector[49]},{elmvector[50]},{elmvector[51]},{elmvector[52]},{elmvector[53]},
        {elmvector[54]},{elmvector[55]},{elmvector[56]},{elmvector[57]},{elmvector[58]},{elmvector[59]},
        {elmvector[60]},{elmvector[61]},{elmvector[62]},{elmvector[63]},{elmvector[64]},{elmvector[65]},
        {elmvector[66]},{elmvector[67]},{elmvector[68]},{elmvector[69]},{elmvector[70]},{elmvector[71]},
        {elmvector[72]},{elmvector[73]},{elmvector[74]},{elmvector[75]},{elmvector[76]},{elmvector[77]},
        {elmvector[78]},{elmvector[79]},{elmvector[80]},{elmvector[81]},{elmvector[82]},{elmvector[83]},
        {elmvector[84]},{elmvector[85]},{elmvector[86]},{elmvector[87]},{elmvector[88]},{elmvector[89]},
        {elmvector[90]},{elmvector[91]},{elmvector[92]},{elmvector[93]},{elmvector[94]},{elmvector[95]},
        {elmvector[96]},{elmvector[97]},{elmvector[98]},{elmvector[99]},{elmvector[100]},{elmvector[101]},
        {elmvector[102]},{elmvector[103]},{elmvector[104]},{elmvector[105]}
    }};

const double *{phase}_{module}_calib_elements(void) {{
    return elmformula;
}}

double {phase}_{module}_calib_g(double T, double P) {{
    return {module}_g(T, P);
}}

double {phase}_{module}_calib_dgdt(double T, double P) {{
    return {module}_dgdt(T, P);
}}

double {phase}_{module}_calib_dgdp(double T, double P) {{
    return {module}_dgdp(T, P);
}}

double {phase}_{module}_calib_d2gdt2(double T, double P) {{
    return {module}_d2gdt2(T, P);
}}

double {phase}_{module}_calib_d2gdtdp(double T, double P) {{
    return {module}_d2gdtdp(T, P);
}}

double {phase}_{module}_calib_d2gdp2(double T, double P) {{
    return {module}_d2gdp2(T, P);
}}

double {phase}_{module}_calib_d3gdt3(double T, double P) {{
    return {module}_d3gdt3(T, P);
}}

double {phase}_{module}_calib_d3gdt2dp(double T, double P) {{
    return {module}_d3gdt2dp(T, P);
}}

double {phase}_{module}_calib_d3gdtdp2(double T, double P) {{
    return {module}_d3gdtdp2(T, P);
}}

double {phase}_{module}_calib_d3gdp3(double T, double P) {{
    return {module}_d3gdp3(T, P);
}}

double {phase}_{module}_calib_s(double T, double P) {{
    return {module}_s(T, P);
}}

double {phase}_{module}_calib_v(double T, double P) {{
    return {module}_v(T, P);
}}

double {phase}_{module}_calib_cv(double T, double P) {{
    return {module}_cv(T, P);
}}

double {phase}_{module}_calib_cp(double T, double P) {{
    return {module}_cp(T, P);
}}

double {phase}_{module}_calib_dcpdt(double T, double P) {{
    return {module}_dcpdt(T, P);
}}

double {phase}_{module}_calib_alpha(double T, double P) {{
    return {module}_alpha(T, P);
}}

double {phase}_{module}_calib_beta(double T, double P) {{
    return {module}_beta(T, P);
}}

double {phase}_{module}_calib_K(double T, double P) {{
    return {module}_K(T, P);
}}

double {phase}_{module}_calib_Kp(double T, double P) {{
    return {module}_Kp(T, P);
}}

int {phase}_{module}_get_param_number(void) {{
    return {module}_get_param_number();
}}

const char **{phase}_{module}_get_param_names(void) {{
    return {module}_get_param_names();
}}

const char **{phase}_{module}_get_param_units(void) {{
    return {module}_get_param_units();
}}

void {phase}_{module}_get_param_values(double **values) {{
    {module}_get_param_values(values);
}}

int {phase}_{module}_set_param_values(double *values) {{
    return {module}_set_param_values(values);
}}

double {phase}_{module}_get_param_value(int index) {{
    return {module}_get_param_value(index);
}}

int {phase}_{module}_set_param_value(int index, double value) {{
    return {module}_set_param_value(index, value);
}}

double {phase}_{module}_dparam_g(double T, double P, int index) {{
    return {module}_dparam_g(T, P, index);
}}

double {phase}_{module}_dparam_dgdt(double T, double P, int index) {{
    return {module}_dparam_dgdt(T, P, index);
}}

double {phase}_{module}_dparam_dgdp(double T, double P, int index) {{
    return {module}_dparam_dgdp(T, P, index);
}}

double {phase}_{module}_dparam_d2gdt2(double T, double P, int index) {{
    return {module}_dparam_d2gdt2(T, P, index);
}}

double {phase}_{module}_dparam_d2gdtdp(double T, double P, int index) {{
    return {module}_dparam_d2gdtdp(T, P, index);
}}

double {phase}_{module}_dparam_d2gdp2(double T, double P, int index) {{
    return {module}_dparam_d2gdp2(T, P, index);
}}

double {phase}_{module}_dparam_d3gdt3(double T, double P, int index) {{
    return {module}_dparam_d3gdt3(T, P, index);
}}

double {phase}_{module}_dparam_d3gdt2dp(double T, double P, int index) {{
    return {module}_dparam_d3gdt2dp(T, P, index);
}}

double {phase}_{module}_dparam_d3gdtdp2(double T, double P, int index) {{
    return {module}_dparam_d3gdtdp2(T, P, index);
}}

double {phase}_{module}_dparam_d3gdp3(double T, double P, int index) {{
    return {module}_dparam_d3gdp3(T, P, index);
}}

\
"""

def _create_calib_h_template_c():
    """
    C language implementation of create_calib_h_template()

    Returns
    -------
    string :
        The template string.
    """

    return """\

const char *{phase}_{module}_calib_identifier(void);
const char *{phase}_{module}_calib_name(void);
const char *{phase}_{module}_calib_formula(void);
const double {phase}_{module}_calib_mw(void);
const double *{phase}_{module}_calib_elements(void);

double {phase}_{module}_calib_g(double T, double P);
double {phase}_{module}_calib_dgdt(double T, double P);
double {phase}_{module}_calib_dgdp(double T, double P);
double {phase}_{module}_calib_d2gdt2(double T, double P);
double {phase}_{module}_calib_d2gdtdp(double T, double P);
double {phase}_{module}_calib_d2gdp2(double T, double P);
double {phase}_{module}_calib_d3gdt3(double T, double P);
double {phase}_{module}_calib_d3gdt2dp(double T, double P);
double {phase}_{module}_calib_d3gdtdp2(double T, double P);
double {phase}_{module}_calib_d3gdp3(double T, double P);

double {phase}_{module}_calib_s(double T, double P);
double {phase}_{module}_calib_v(double T, double P);
double {phase}_{module}_calib_cv(double T, double P);
double {phase}_{module}_calib_cp(double T, double P);
double {phase}_{module}_calib_dcpdt(double T, double P);
double {phase}_{module}_calib_alpha(double T, double P);
double {phase}_{module}_calib_beta(double T, double P);
double {phase}_{module}_calib_K(double T, double P);
double {phase}_{module}_calib_Kp(double T, double P);

int {phase}_{module}_get_param_number(void);
const char **{phase}_{module}_get_param_names(void);
const char **{phase}_{module}_get_param_units(void);
void {phase}_{module}_get_param_values(double **values);
int {phase}_{module}_set_param_values(double *values);
double {phase}_{module}_get_param_value(int index);
int {phase}_{module}_set_param_value(int index, double value);

double {phase}_{module}_dparam_g(double T, double P, int index);
double {phase}_{module}_dparam_dgdt(double T, double P, int index);
double {phase}_{module}_dparam_dgdp(double T, double P, int index);
double {phase}_{module}_dparam_d2gdt2(double T, double P, int index);
double {phase}_{module}_dparam_d2gdtdp(double T, double P, int index);
double {phase}_{module}_dparam_d2gdp2(double T, double P, int index);
double {phase}_{module}_dparam_d3gdt3(double T, double P, int index);
double {phase}_{module}_dparam_d3gdt2dp(double T, double P, int index);
double {phase}_{module}_dparam_d3gdtdp2(double T, double P, int index);
double {phase}_{module}_dparam_d3gdp3(double T, double P, int index);

\
"""

def _create_calib_pyx_template_c():
    """
    C language implementation of create_calib_pyx_template()

    Returns
    -------
    string :
        The template string.
    """

    return """\
# Cython numpy wrapper code for arrays is taken from:
# http://gael-varoquaux.info/programming/cython-example-of-exposing-c-computed-arrays-in-python-without-data-copies.html
# Author: Gael Varoquaux, BSD license

# Declare the prototype of the C functions
cdef extern from "{phase}_{module}_calib.h":
    const char *{phase}_{module}_calib_identifier();
    const char *{phase}_{module}_calib_name();
    const char *{phase}_{module}_calib_formula();
    const double {phase}_{module}_calib_mw();
    const double *{phase}_{module}_calib_elements();
    double {phase}_{module}_calib_g(double t, double p);
    double {phase}_{module}_calib_dgdt(double t, double p);
    double {phase}_{module}_calib_dgdp(double t, double p);
    double {phase}_{module}_calib_d2gdt2(double t, double p);
    double {phase}_{module}_calib_d2gdtdp(double t, double p);
    double {phase}_{module}_calib_d2gdp2(double t, double p);
    double {phase}_{module}_calib_d3gdt3(double t, double p);
    double {phase}_{module}_calib_d3gdt2dp(double t, double p);
    double {phase}_{module}_calib_d3gdtdp2(double t, double p);
    double {phase}_{module}_calib_d3gdp3(double t, double p);
    double {phase}_{module}_calib_s(double t, double p);
    double {phase}_{module}_calib_v(double t, double p);
    double {phase}_{module}_calib_cv(double t, double p);
    double {phase}_{module}_calib_cp(double t, double p);
    double {phase}_{module}_calib_dcpdt(double t, double p);
    double {phase}_{module}_calib_alpha(double t, double p);
    double {phase}_{module}_calib_beta(double t, double p);
    double {phase}_{module}_calib_K(double t, double p);
    double {phase}_{module}_calib_Kp(double t, double p);
    int {phase}_{module}_get_param_number();
    const char **{phase}_{module}_get_param_names();
    const char **{phase}_{module}_get_param_units();
    void {phase}_{module}_get_param_values(double **values);
    int {phase}_{module}_set_param_values(double *values);
    double {phase}_{module}_get_param_value(int index);
    int {phase}_{module}_set_param_value(int index, double value);
    double {phase}_{module}_dparam_g(double t, double p, int index);
    double {phase}_{module}_dparam_dgdt(double t, double p, int index);
    double {phase}_{module}_dparam_dgdp(double t, double p, int index);
    double {phase}_{module}_dparam_d2gdt2(double t, double p, int index);
    double {phase}_{module}_dparam_d2gdtdp(double t, double p, int index);
    double {phase}_{module}_dparam_d2gdp2(double t, double p, int index);
    double {phase}_{module}_dparam_d3gdt3(double t, double p, int index);
    double {phase}_{module}_dparam_d3gdt2dp(double t, double p, int index);
    double {phase}_{module}_dparam_d3gdtdp2(double t, double p, int index);
    double {phase}_{module}_dparam_d3gdp3(double t, double p, int index);

from libc.stdlib cimport malloc, free
from cpython cimport PyObject, Py_INCREF
import ctypes

# Import the Python-level symbols of numpy
import numpy as np

# Import the C-level symbols of numpy
cimport numpy as np

# Numpy must be initialized. When using numpy from C or Cython you must
# _always_ do that, or you will have segfaults
np.import_array()

# here is the "wrapper" signature
def {prefix}_{phase}_{module}_calib_identifier():
    result = <bytes> {phase}_{module}_calib_identifier()
    return result.decode('UTF-8')
def {prefix}_{phase}_{module}_calib_name():
    result = <bytes> {phase}_{module}_calib_name()
    return result.decode('UTF-8')
def {prefix}_{phase}_{module}_calib_formula():
    result = <bytes> {phase}_{module}_calib_formula()
    return result.decode('UTF-8')
def {prefix}_{phase}_{module}_calib_mw():
    result = {phase}_{module}_calib_mw()
    return result
def {prefix}_{phase}_{module}_calib_elements():
    cdef const double *e = {phase}_{module}_calib_elements()
    np_array = np.zeros(106)
    for i in range(0,106):
        np_array[i] = e[i]
    return np_array
def {prefix}_{phase}_{module}_calib_g(double t, double p):
    result = {phase}_{module}_calib_g(<double> t, <double> p)
    return result
def {prefix}_{phase}_{module}_calib_dgdt(double t, double p):
    result = {phase}_{module}_calib_dgdt(<double> t, <double> p)
    return result
def {prefix}_{phase}_{module}_calib_dgdp(double t, double p):
    result = {phase}_{module}_calib_dgdp(<double> t, <double> p)
    return result
def {prefix}_{phase}_{module}_calib_d2gdt2(double t, double p):
    result = {phase}_{module}_calib_d2gdt2(<double> t, <double> p)
    return result
def {prefix}_{phase}_{module}_calib_d2gdtdp(double t, double p):
    result = {phase}_{module}_calib_d2gdtdp(<double> t, <double> p)
    return result
def {prefix}_{phase}_{module}_calib_d2gdp2(double t, double p):
    result = {phase}_{module}_calib_d2gdp2(<double> t, <double> p)
    return result
def {prefix}_{phase}_{module}_calib_d3gdt3(double t, double p):
    result = {phase}_{module}_calib_d3gdt3(<double> t, <double> p)
    return result
def {prefix}_{phase}_{module}_calib_d3gdt2dp(double t, double p):
    result = {phase}_{module}_calib_d3gdt2dp(<double> t, <double> p)
    return result
def {prefix}_{phase}_{module}_calib_d3gdtdp2(double t, double p):
    result = {phase}_{module}_calib_d3gdtdp2(<double> t, <double> p)
    return result
def {prefix}_{phase}_{module}_calib_d3gdp3(double t, double p):
    result = {phase}_{module}_calib_d3gdp3(<double> t, <double> p)
    return result
def {prefix}_{phase}_{module}_calib_s(double t, double p):
    result = {phase}_{module}_calib_s(<double> t, <double> p)
    return result
def {prefix}_{phase}_{module}_calib_v(double t, double p):
    result = {phase}_{module}_calib_v(<double> t, <double> p)
    return result
def {prefix}_{phase}_{module}_calib_cv(double t, double p):
    result = {phase}_{module}_calib_cv(<double> t, <double> p)
    return result
def {prefix}_{phase}_{module}_calib_cp(double t, double p):
    result = {phase}_{module}_calib_cp(<double> t, <double> p)
    return result
def {prefix}_{phase}_{module}_calib_dcpdt(double t, double p):
    result = {phase}_{module}_calib_dcpdt(<double> t, <double> p)
    return result
def {prefix}_{phase}_{module}_calib_alpha(double t, double p):
    result = {phase}_{module}_calib_alpha(<double> t, <double> p)
    return result
def {prefix}_{phase}_{module}_calib_beta(double t, double p):
    result = {phase}_{module}_calib_beta(<double> t, <double> p)
    return result
def {prefix}_{phase}_{module}_calib_K(double t, double p):
    result = {phase}_{module}_calib_K(<double> t, <double> p)
    return result
def {prefix}_{phase}_{module}_calib_Kp(double t, double p):
    result = {phase}_{module}_calib_Kp(<double> t, <double> p)
    return result
def {prefix}_{phase}_{module}_get_param_number():
    result = {phase}_{module}_get_param_number()
    return result
def {prefix}_{phase}_{module}_get_param_names():
    cdef const char **names = {phase}_{module}_get_param_names()
    n = {phase}_{module}_get_param_number()
    result = []
    for i in range(0,n):
        entry = <bytes> names[i]
        result.append(entry.decode('UTF-8'))
    return result
def {prefix}_{phase}_{module}_get_param_units():
    cdef const char **units = {phase}_{module}_get_param_units()
    n = {phase}_{module}_get_param_number()
    result = []
    for i in range(0,n):
        entry = <bytes> units[i]
        result.append(entry.decode('UTF-8'))
    return result
def {prefix}_{phase}_{module}_get_param_values():
    n = {phase}_{module}_get_param_number()
    cdef double *m = <double *>malloc(n*sizeof(double))
    {phase}_{module}_get_param_values(&m)
    np_array = np.zeros(n)
    for i in range(n):
        np_array[i] = m[i]
    free(m)
    return np_array
def {prefix}_{phase}_{module}_set_param_values(np_array):
    n = len(np_array)
    cdef double *m = <double *>malloc(n*sizeof(double))
    for i in range(np_array.size):
        m[i] = np_array[i]
    result = {phase}_{module}_set_param_values(m);
    free(m)
    return result
def {prefix}_{phase}_{module}_get_param_value(int index):
    result = {phase}_{module}_get_param_value(<int> index)
    return result
def {prefix}_{phase}_{module}_set_param_value(int index, double value):
    result = {phase}_{module}_set_param_value(<int> index, <double> value)
    return result
def {prefix}_{phase}_{module}_dparam_g(double t, double p, int index):
    result = {phase}_{module}_dparam_g(<double> t, <double> p, <int> index)
    return result
def {prefix}_{phase}_{module}_dparam_dgdt(double t, double p, int index):
    result = {phase}_{module}_dparam_dgdt(<double> t, <double> p, <int> index)
    return result
def {prefix}_{phase}_{module}_dparam_dgdp(double t, double p, int index):
    result = {phase}_{module}_dparam_dgdp(<double> t, <double> p, <int> index)
    return result
def {prefix}_{phase}_{module}_dparam_d2gdt2(double t, double p, int index):
    result = {phase}_{module}_dparam_d2gdt2(<double> t, <double> p, <int> index)
    return result
def {prefix}_{phase}_{module}_dparam_d2gdtdp(double t, double p, int index):
    result = {phase}_{module}_dparam_d2gdtdp(<double> t, <double> p, <int> index)
    return result
def {prefix}_{phase}_{module}_dparam_d2gdp2(double t, double p, int index):
    result = {phase}_{module}_dparam_d2gdp2(<double> t, <double> p, <int> index)
    return result
def {prefix}_{phase}_{module}_dparam_d3gdt3(double t, double p, int index):
    result = {phase}_{module}_dparam_d3gdt3(<double> t, <double> p, <int> index)
    return result
def {prefix}_{phase}_{module}_dparam_d3gdt2dp(double t, double p, int index):
    result = {phase}_{module}_dparam_d3gdt2dp(<double> t, <double> p, <int> index)
    return result
def {prefix}_{phase}_{module}_dparam_d3gdtdp2(double t, double p, int index):
    result = {phase}_{module}_dparam_d3gdtdp2(<double> t, <double> p, <int> index)
    return result
def {prefix}_{phase}_{module}_dparam_d3gdp3(double t, double p, int index):
    result = {phase}_{module}_dparam_d3gdp3(<double> t, <double> p, <int> index)
    return result
\
"""

def _create_code_for_born_functions_c():
    """
    C language implementation of create_code_for_born_functions()

    Returns
    -------
    string :
        The template string.
    """
    return """\
double born_B(double t, double p);
double born_Q(double t, double p);
double born_N(double t, double p);
double born_U(double t, double p);
double born_Y(double t, double p);
double born_X(double t, double p);
double born_dUdT(double t, double p);
double born_dUdP(double t, double p);
double born_dNdT(double t, double p);
double born_dNdP(double t, double p);
double born_dXdT(double t, double p);
double gSolvent(double t, double p);
double DgSolventDt(double t, double p);
double DgSolventDp(double t, double p);
double D2gSolventDt2(double t, double p);
double D2gSolventDtDp(double t, double p);
double D2gSolventDp2(double t, double p);
double D3gSolventDt3(double t, double p);
double D3gSolventDt2Dp(double t, double p);
double D3gSolventDtDp2(double t, double p);
double D3gSolventDp3(double t, double p);
double D4gSolventDt4(double t, double p);
\
"""

def _create_code_for_debye_function_c():
    """
    C language implementation of create_code_for_debye_function()

    Returns
    -------
    string :
        The template string.
    """
    return """\
#include <float.h>
#include <assert.h>

static double chebvalat(double x) {
    double c[17] = {
        2.707737068327440945 / 2.0, 0.340068135211091751, -0.12945150184440869e-01, 0.7963755380173816e-03,
        -0.546360009590824e-04, 0.39243019598805e-05, -0.2894032823539e-06, 0.217317613962e-07, -0.16542099950e-08,
        0.1272796189e-09, -0.987963460e-11, 0.7725074e-12, -0.607797e-13, 0.48076e-14, -0.3820e-15, 0.305e-16, -0.24e-17
    };
    double x2 = 2 * x;
    double c0 = c[17-2];
    double c1 = c[17-1];
    for (int i=3; i<18; i++) {
        double tmp = c0;
        c0 = c[17-i] - c1;
        c1 = tmp + c1 * x2;
    }
    return c0 + c1 * x;
}

static double Debye(double x) {
    //
    // returns D_3(x) = 3/x^3\int_0^t t^3/(e^t - 1) dt
    //
    
    double val_infinity = 19.4818182068004875;
    double sqrt_eps = sqrt(DBL_EPSILON);
    double log_eps = log(DBL_EPSILON);
    double xcut = -log_eps;

    //Check for negative x (was returning zero)
    assert(x >= 0.);

    if (x < (2.0*sqrt(2.0)*sqrt_eps)) return 1.0 - 3.0*x/8.0 + x*x/20.0;
    else if (x <= 4.0) {
        double t = x*x/8.0 - 1.0;
        double c = chebvalat(t);
        return c - 0.375*x;
    } else if (x < -(log(2.0)+log_eps)) {
        int nexp = (int)(floor(xcut / x));
        double ex = exp(-x);
        double xk = nexp * x;
        double rk = nexp;
        double sum = 0.0;
        for (int i=nexp; i>0; i--) {
            double xk_inv = 1.0/xk;
            sum *= ex;
            sum += (((6.0*xk_inv + 6.0)*xk_inv + 3.0)*xk_inv + 1.0)/rk;
            rk -= 1.0;
            xk -= x;
        }
        return val_infinity / (x * x * x) - 3.0 * sum * ex;
    } else if (x < xcut) {
        double x3 = x*x*x;
        double sum = 6.0 + 6.0*x + 3.0*x*x + x3;
        return (val_infinity - 3.0*sum*exp(-x))/x3;
    } else return ((val_infinity/x)/x)/x;
}

\
"""

def _create_code_for_dh_functions_c():
    """
    C language implementation of create_code_for_dh_functions()

    Returns
    -------
    string :
        The template string.
    """
    return """\
double Agamma(double t, double p);
double dAgammaDt(double t, double p);
double dAgammaDp(double t, double p);
double d2AgammaDt2(double t, double p);
double d2AgammaDtDp(double t, double p);
double d2AgammaDp2(double t, double p);
double d3AgammaDt3(double t, double p);
double d3AgammaDt2Dp(double t, double p);
double d3AgammaDtDp2(double t, double p);
double d3AgammaDp3(double t, double p);

double Bgamma(double t, double p);
double dBgammaDt(double t, double p);
double dBgammaDp(double t, double p);
double d2BgammaDt2(double t, double p);
double d2BgammaDtDp(double t, double p);
double d2BgammaDp2(double t, double p);
double d3BgammaDt3(double t, double p);
double d3BgammaDt2Dp(double t, double p);
double d3BgammaDtDp2(double t, double p);
double d3BgammaDp3(double t, double p);

double AsubG(double t, double p);
double AsubH(double t, double p);
double AsubJ(double t, double p);
double AsubV(double t, double p);
double AsubKappa(double t, double p);
double AsubEx(double t, double p);
double BsubG(double t, double p);
double BsubH(double t, double p);
double BsubJ(double t, double p);
double BsubV(double t, double p);
double BsubKappa(double t, double p);
double BsubEx(double t, double p);
\
"""

def _create_fast_c_template_c():
    """
    C language implementation of create_fast_c_template()

    Returns
    -------
    string :
        The template string.
    """

    return """\

static const char *identifier = "{git_identifier}";
{parameter_init_block}
{include_calc_h}
const char *{phase}_{module}_identifier(void) {{
    return identifier;
}}

const char *{phase}_{module}_name(void) {{
    return "{phase}";
}}

const char *{phase}_{module}_formula(void) {{
    return "{formula}";
}}

const double {phase}_{module}_mw(void) {{
    return {mw};
}}

static const double elmformula[106] = {{
        {elmvector[0]},{elmvector[1]},{elmvector[2]},{elmvector[3]},{elmvector[4]},{elmvector[5]},
        {elmvector[6]},{elmvector[7]},{elmvector[8]},{elmvector[9]},{elmvector[10]},{elmvector[11]},
        {elmvector[12]},{elmvector[13]},{elmvector[14]},{elmvector[15]},{elmvector[16]},{elmvector[17]},
        {elmvector[18]},{elmvector[19]},{elmvector[20]},{elmvector[21]},{elmvector[22]},{elmvector[23]},
        {elmvector[24]},{elmvector[25]},{elmvector[26]},{elmvector[27]},{elmvector[28]},{elmvector[29]},
        {elmvector[30]},{elmvector[31]},{elmvector[32]},{elmvector[33]},{elmvector[34]},{elmvector[35]},
        {elmvector[36]},{elmvector[37]},{elmvector[38]},{elmvector[39]},{elmvector[40]},{elmvector[41]},
        {elmvector[42]},{elmvector[43]},{elmvector[44]},{elmvector[45]},{elmvector[46]},{elmvector[47]},
        {elmvector[48]},{elmvector[49]},{elmvector[50]},{elmvector[51]},{elmvector[52]},{elmvector[53]},
        {elmvector[54]},{elmvector[55]},{elmvector[56]},{elmvector[57]},{elmvector[58]},{elmvector[59]},
        {elmvector[60]},{elmvector[61]},{elmvector[62]},{elmvector[63]},{elmvector[64]},{elmvector[65]},
        {elmvector[66]},{elmvector[67]},{elmvector[68]},{elmvector[69]},{elmvector[70]},{elmvector[71]},
        {elmvector[72]},{elmvector[73]},{elmvector[74]},{elmvector[75]},{elmvector[76]},{elmvector[77]},
        {elmvector[78]},{elmvector[79]},{elmvector[80]},{elmvector[81]},{elmvector[82]},{elmvector[83]},
        {elmvector[84]},{elmvector[85]},{elmvector[86]},{elmvector[87]},{elmvector[88]},{elmvector[89]},
        {elmvector[90]},{elmvector[91]},{elmvector[92]},{elmvector[93]},{elmvector[94]},{elmvector[95]},
        {elmvector[96]},{elmvector[97]},{elmvector[98]},{elmvector[99]},{elmvector[100]},{elmvector[101]},
        {elmvector[102]},{elmvector[103]},{elmvector[104]},{elmvector[105]}
    }};

const double *{phase}_{module}_elements(void) {{
    return elmformula;
}}

double {phase}_{module}_g(double T, double P) {{
    return {module}_g(T, P);
}}

double {phase}_{module}_dgdt(double T, double P) {{
    return {module}_dgdt(T, P);
}}

double {phase}_{module}_dgdp(double T, double P) {{
    return {module}_dgdp(T, P);
}}

double {phase}_{module}_d2gdt2(double T, double P) {{
    return {module}_d2gdt2(T, P);
}}

double {phase}_{module}_d2gdtdp(double T, double P) {{
    return {module}_d2gdtdp(T, P);
}}

double {phase}_{module}_d2gdp2(double T, double P) {{
    return {module}_d2gdp2(T, P);
}}

double {phase}_{module}_d3gdt3(double T, double P) {{
    return {module}_d3gdt3(T, P);
}}

double {phase}_{module}_d3gdt2dp(double T, double P) {{
    return {module}_d3gdt2dp(T, P);
}}

double {phase}_{module}_d3gdtdp2(double T, double P) {{
    return {module}_d3gdtdp2(T, P);
}}

double {phase}_{module}_d3gdp3(double T, double P) {{
    return {module}_d3gdp3(T, P);
}}

double {phase}_{module}_s(double T, double P) {{
    return {module}_s(T, P);
}}

double {phase}_{module}_v(double T, double P) {{
    return {module}_v(T, P);
}}

double {phase}_{module}_cv(double T, double P) {{
    return {module}_cv(T, P);
}}

double {phase}_{module}_cp(double T, double P) {{
    return {module}_cp(T, P);
}}

double {phase}_{module}_dcpdt(double T, double P) {{
    return {module}_dcpdt(T, P);
}}

double {phase}_{module}_alpha(double T, double P) {{
    return {module}_alpha(T, P);
}}

double {phase}_{module}_beta(double T, double P) {{
    return {module}_beta(T, P);
}}

double {phase}_{module}_K(double T, double P) {{
    return {module}_K(T, P);
}}

double {phase}_{module}_Kp(double T, double P) {{
    return {module}_Kp(T, P);
}}

\
"""

def _create_fast_h_template_c():
    """
    C language implementation of create_fast_h_template()

    Returns
    -------
    string :
        The template string.
    """

    return """\

const char *{phase}_{module}_identifier(void);
const char *{phase}_{module}_name(void);
const char *{phase}_{module}_formula(void);
const double {phase}_{module}_mw(void);
const double *{phase}_{module}_elements(void);

double {phase}_{module}_g(double T, double P);
double {phase}_{module}_dgdt(double T, double P);
double {phase}_{module}_dgdp(double T, double P);
double {phase}_{module}_d2gdt2(double T, double P);
double {phase}_{module}_d2gdtdp(double T, double P);
double {phase}_{module}_d2gdp2(double T, double P);
double {phase}_{module}_d3gdt3(double T, double P);
double {phase}_{module}_d3gdt2dp(double T, double P);
double {phase}_{module}_d3gdtdp2(double T, double P);
double {phase}_{module}_d3gdp3(double T, double P);

double {phase}_{module}_s(double T, double P);
double {phase}_{module}_v(double T, double P);
double {phase}_{module}_cv(double T, double P);
double {phase}_{module}_cp(double T, double P);
double {phase}_{module}_dcpdt(double T, double P);
double {phase}_{module}_alpha(double T, double P);
double {phase}_{module}_beta(double T, double P);
double {phase}_{module}_K(double T, double P);
double {phase}_{module}_Kp(double T, double P);

\
"""

def _create_fast_pyx_template_c():
    """
    C language implementation of create_fast_pyx_template()

    Returns
    -------
    string :
        The template string.
    """

    return """\
import numpy as np
cimport numpy as cnp # cimport gives us access to NumPy's C API

# here we just replicate the function signature from the header
cdef extern from "{phase}_{module}_calc.h":
    const char *{phase}_{module}_identifier();
    const char *{phase}_{module}_name();
    const char *{phase}_{module}_formula();
    const double {phase}_{module}_mw();
    const double *{phase}_{module}_elements();
    double {phase}_{module}_g(double t, double p)
    double {phase}_{module}_dgdt(double t, double p)
    double {phase}_{module}_dgdp(double t, double p)
    double {phase}_{module}_d2gdt2(double t, double p)
    double {phase}_{module}_d2gdtdp(double t, double p)
    double {phase}_{module}_d2gdp2(double t, double p)
    double {phase}_{module}_d3gdt3(double t, double p)
    double {phase}_{module}_d3gdt2dp(double t, double p)
    double {phase}_{module}_d3gdtdp2(double t, double p)
    double {phase}_{module}_d3gdp3(double t, double p)
    double {phase}_{module}_s(double t, double p)
    double {phase}_{module}_v(double t, double p)
    double {phase}_{module}_cv(double t, double p)
    double {phase}_{module}_cp(double t, double p)
    double {phase}_{module}_dcpdt(double t, double p)
    double {phase}_{module}_alpha(double t, double p)
    double {phase}_{module}_beta(double t, double p)
    double {phase}_{module}_K(double t, double p)
    double {phase}_{module}_Kp(double t, double p)

# here is the "wrapper" signature
def {prefix}_{phase}_{module}_identifier():
    result = <bytes> {phase}_{module}_identifier()
    return result.decode('UTF-8')
def {prefix}_{phase}_{module}_name():
    result = <bytes> {phase}_{module}_name()
    return result.decode('UTF-8')
def {prefix}_{phase}_{module}_formula():
    result = <bytes> {phase}_{module}_formula()
    return result.decode('UTF-8')
def {prefix}_{phase}_{module}_mw():
    result = {phase}_{module}_mw()
    return result
def {prefix}_{phase}_{module}_elements():
    cdef const double *e = {phase}_{module}_elements()
    np_array = np.zeros(106)
    for i in range(0,106):
        np_array[i] = e[i]
    return np_array
def {prefix}_{phase}_{module}_g(double t, double p):
    result = {phase}_{module}_g(<double> t, <double> p)
    return result
def {prefix}_{phase}_{module}_dgdt(double t, double p):
    result = {phase}_{module}_dgdt(<double> t, <double> p)
    return result
def {prefix}_{phase}_{module}_dgdp(double t, double p):
    result = {phase}_{module}_dgdp(<double> t, <double> p)
    return result
def {prefix}_{phase}_{module}_d2gdt2(double t, double p):
    result = {phase}_{module}_d2gdt2(<double> t, <double> p)
    return result
def {prefix}_{phase}_{module}_d2gdtdp(double t, double p):
    result = {phase}_{module}_d2gdtdp(<double> t, <double> p)
    return result
def {prefix}_{phase}_{module}_d2gdp2(double t, double p):
    result = {phase}_{module}_d2gdp2(<double> t, <double> p)
    return result
def {prefix}_{phase}_{module}_d3gdt3(double t, double p):
    result = {phase}_{module}_d3gdt3(<double> t, <double> p)
    return result
def {prefix}_{phase}_{module}_d3gdt2dp(double t, double p):
    result = {phase}_{module}_d3gdt2dp(<double> t, <double> p)
    return result
def {prefix}_{phase}_{module}_d3gdtdp2(double t, double p):
    result = {phase}_{module}_d3gdtdp2(<double> t, <double> p)
    return result
def {prefix}_{phase}_{module}_d3gdp3(double t, double p):
    result = {phase}_{module}_d3gdp3(<double> t, <double> p)
    return result
def {prefix}_{phase}_{module}_s(double t, double p):
    result = {phase}_{module}_s(<double> t, <double> p)
    return result
def {prefix}_{phase}_{module}_v(double t, double p):
    result = {phase}_{module}_v(<double> t, <double> p)
    return result
def {prefix}_{phase}_{module}_cv(double t, double p):
    result = {phase}_{module}_cv(<double> t, <double> p)
    return result
def {prefix}_{phase}_{module}_cp(double t, double p):
    result = {phase}_{module}_cp(<double> t, <double> p)
    return result
def {prefix}_{phase}_{module}_dcpdt(double t, double p):
    result = {phase}_{module}_dcpdt(<double> t, <double> p)
    return result
def {prefix}_{phase}_{module}_alpha(double t, double p):
    result = {phase}_{module}_alpha(<double> t, <double> p)
    return result
def {prefix}_{phase}_{module}_beta(double t, double p):
    result = {phase}_{module}_beta(<double> t, <double> p)
    return result
def {prefix}_{phase}_{module}_K(double t, double p):
    result = {phase}_{module}_K(<double> t, <double> p)
    return result
def {prefix}_{phase}_{module}_Kp(double t, double p):
    result = {phase}_{module}_Kp(<double> t, <double> p)
    return result
\
"""

def _create_pyxbld_template_c():
    """
    C language implementation of create_pyxbld_template()

    Returns
    -------
    string :
        The template string.

    Notes
    -----
    setuptools.extension.Extension
    self,
    name,
    sources,
    include_dirs=None,
    define_macros=None,
    undef_macros=None,
    library_dirs=None,
    libraries=None,
    runtime_library_dirs=None,
    extra_objects=None,
    extra_compile_args=None,
    extra_link_args=None,
    export_symbols=None,
    swig_opts=None,
    depends=None,
    language=None,
    optional=None,
    **kw

    -O0, -O1, -O2, -O3, -Ofast, -Os, -Oz, -Og, -O, -O4
    Specify which optimization level to use:

        -O0 Means “no optimization”: this level compiles the fastest and
        generates the most debuggable code.
        -O1 Somewhere between -O0 and -O2.
        -O2 Moderate level of optimization which enables most optimizations.
        -O3 Like -O2, except that it enables optimizations that take longer
        to perform or that may generate larger code (in an attempt to make
        the program run faster).
        -Ofast Enables all the optimizations from -O3 along with other
        agressive optimizations that may violate strict compliance with
        language standards.
        -Os Like -O2 with extra optimizations to reduce code size.
        -Oz Like -Os (and thus -O2), but reduces code size further.
        -Og Like -O1. In future versions, this option might disable
        different optimizations in order to improve debuggability.
        -O Equivalent to -O2.
        -O4 and higher
    """

    return """\
import numpy

#            module name specified by `%%cython_pyximport` magic
#            |        just `modname + ".pyx"`
#            |        |
def make_ext(modname, pyxfilename):
    from setuptools.extension import Extension
    return Extension(modname,
                     sources=[pyxfilename, '{file_to_compile}'],
                     include_dirs=['.', numpy.get_include()],
                     extra_compile_args=['-O3'],
                     libraries=['gsl', 'swimdew'],
                     library_dirs=['/usr/local/lib'],
                     runtime_library_dirs=['/usr/local/lib'])
"""

def _create_redundant_calib_TV_template_c():
    """
    C language implementation of create_redundant_calib_TV_template()

    Returns
    -------
    string :
        The template string.
    """

    return """\

static double {module}_dparam_g(double T, double P, int index) {{
    {module}_solve_V(T, P);
    double dAdz = {module}_dparam_a(T, V, index);
    return dAdz + P*V;
}}

static double {module}_dparam_dgdt(double T, double P, int index) {{
    {module}_solve_V(T, P);
    double dAdTdz = {module}_dparam_dadt(T, V, index);
    return dAdTdz;
}}

static double {module}_dparam_dgdp(double T, double P, int index) {{
    {module}_solve_V(T, P);
    return 0.0; /* V; */
}}

static double {module}_dparam_d2gdt2(double T, double P, int index) {{
    {module}_solve_V(T, P);
    double d2AdTdV = {module}_d2adtdv(T, V);
    double d2AdV2 = {module}_d2adv2(T, V);
    double d2AdT2dz = {module}_dparam_d2adt2(T, V, index);
    double d2AdTdVdz = {module}_dparam_d2adtdv(T, V, index);
    double d2AdV2dz = {module}_dparam_d2adv2(T, V, index);
    /* return d2AdT2 - d2AdTdV*d2AdTdV/d2AdV2; */
    return d2AdT2dz - 2.0*d2AdTdV*d2AdTdVdz/d2AdV2 + d2AdTdV*d2AdTdV*d2AdV2dz/d2AdV2/d2AdV2;
}}

static double {module}_dparam_d2gdtdp(double T, double P, int index) {{
    {module}_solve_V(T, P);
    double d2AdTdV = {module}_d2adtdv(T, V);
    double d2AdV2 = {module}_d2adv2(T, V);
    double d2AdTdVdz = {module}_dparam_d2adtdv(T, V, index);
    double d2AdV2dz = {module}_dparam_d2adv2(T, V, index);
    /* return - d2AdTdV/d2AdV2; */
    return - d2AdTdVdz/d2AdV2 + d2AdTdV*d2AdV2dz/d2AdV2/d2AdV2;
}}

static double {module}_dparam_d2gdp2(double T, double P, int index) {{
    {module}_solve_V(T, P);
    double d2AdV2 = {module}_d2adv2(T, V);
    double d2AdV2dz = {module}_dparam_d2adv2(T, V, index);
    /* return - 1.0/d2AdV2; */
    return d2AdV2dz/d2AdV2/d2AdV2;
}}

static double {module}_dparam_d3gdt3(double T, double P, int index) {{
    {module}_solve_V(T, P);
    double d3AdT2dV = {module}_d3adt2dv(T, V);
    double d3AdTdV2 = {module}_d3adtdv2(T, V);
    double d3AdV3 = {module}_d3adv3(T,V);
    double d2AdTdV = {module}_d2adtdv(T, V);
    double d2AdV2 = {module}_d2adv2(T, V);
    double d3AdT3dz = {module}_dparam_d3adt3(T, V, index);
    double d3AdT2dVdz = {module}_dparam_d3adt2dv(T, V, index);
    double d3AdTdV2dz = {module}_dparam_d3adtdv2(T, V, index);
    double d3AdV3dz = {module}_dparam_d3adv3(T, V, index);
    double d2AdTdVdz = {module}_dparam_d2adtdv(T, V, index);
    double d2AdV2dz = {module}_dparam_d2adv2(T, V, index);
    double dVdT = - d2AdTdV/d2AdV2;
    double dVdTdz = - d2AdTdVdz/d2AdV2 + d2AdTdV*d2AdV2dz/d2AdV2/d2AdV2;
    double d2VdT2 = (-d3AdT2dV - 2.0*d3AdTdV2*dVdT - d3AdV3*dVdT*dVdT)/d2AdV2;
    double d2VdT2dz = (-d3AdT2dVdz - 2.0*d3AdTdV2dz*dVdT - 2.0*d3AdTdV2*dVdTdz
                        - d3AdV3dz*dVdT*dVdT - 2.0*d3AdV3*dVdT*dVdTdz)/d2AdV2
                    - (-d3AdT2dV - 2.0*d3AdTdV2*dVdT - d3AdV3*dVdT*dVdT)*d2AdV2dz/d2AdV2/d2AdV2;
    /* return d3AdT3 + 2.0*d3AdT2dV*dVdT + d3AdTdV2*dVdT*dVdT + d2AdTdV*d2VdT2; */
    return d3AdT3dz + 2.0*d3AdT2dVdz*dVdT + 2.0*d3AdT2dV*dVdTdz + d3AdTdV2dz*dVdT*dVdT
            + 2.0*d3AdTdV2*dVdT*dVdTdz + d2AdTdVdz*d2VdT2 + d2AdTdV*d2VdT2dz;
}}

static double {module}_dparam_d3gdt2dp(double T, double P, int index) {{
    {module}_solve_V(T, P);
    double d3AdT2dV = {module}_d3adt2dv(T, V);
    double d3AdTdV2 = {module}_d3adtdv2(T, V);
    double d3AdV3 = {module}_d3adv3(T,V);
    double d2AdTdV = {module}_d2adtdv(T, V);
    double d2AdV2 = {module}_d2adv2(T, V);
    double d3AdT2dVdz = {module}_dparam_d3adt2dv(T, V, index);
    double d3AdTdV2dz = {module}_dparam_d3adtdv2(T, V, index);
    double d3AdV3dz = {module}_dparam_d3adv3(T,V, index);
    double d2AdTdVdz = {module}_dparam_d2adtdv(T, V, index);
    double d2AdV2dz = {module}_dparam_d2adv2(T, V, index);
    double dVdT = -d2AdTdV/d2AdV2;
    double dVdTdz = - d2AdTdVdz/d2AdV2 + d2AdTdV*d2AdV2dz/d2AdV2/d2AdV2;
    double dVdP = -1.0/d2AdV2;
    double dVdPdz = d2AdV2dz/d2AdV2/d2AdV2;
    double d2VdTdP = (-d3AdTdV2*dVdP - d3AdV3*dVdT*dVdP)/d2AdV2;
    double d2VdTdPdz = (-d3AdTdV2dz*dVdP -d3AdTdV2*dVdPdz
            - d3AdV3dz*dVdT*dVdP - d3AdV3*dVdTdz*dVdP - d3AdV3*dVdT*dVdPdz)/d2AdV2
            - (-d3AdTdV2*dVdP - d3AdV3*dVdT*dVdP)*d2AdV2dz/d2AdV2/d2AdV2;
    /* return d3AdT2dV*dVdP + d3AdTdV2*dVdT*dVdP + d2AdTdV*d2VdTdP; */
    return d3AdT2dVdz*dVdP + d3AdT2dV*dVdPdz + d3AdTdV2dz*dVdT*dVdP + d3AdTdV2*dVdTdz*dVdP
        + d3AdTdV2*dVdT*dVdPdz + d2AdTdVdz*d2VdTdP + d2AdTdV*d2VdTdPdz;
}}

static double {module}_dparam_d3gdtdp2(double T, double P, int index) {{
    {module}_solve_V(T, P);
    double d3AdTdV2 = {module}_d3adtdv2(T, V);
    double d3AdV3 = {module}_d3adv3(T,V);
    double d2AdTdV = {module}_d2adtdv(T, V);
    double d2AdV2 = {module}_d2adv2(T, V);
    double d3AdTdV2dz = {module}_dparam_d3adtdv2(T, V, index);
    double d3AdV3dz = {module}_dparam_d3adv3(T, V, index);
    double d2AdTdVdz = {module}_dparam_d2adtdv(T, V, index);
    double d2AdV2dz = {module}_dparam_d2adv2(T, V, index);
    double dVdT = -d2AdTdV/d2AdV2;
    double dVdTdz = - d2AdTdVdz/d2AdV2 + d2AdTdV*d2AdV2dz/d2AdV2/d2AdV2;
    /* return (d3AdTdV2 + d3AdV3*dVdT)/d2AdV2/d2AdV2; */
    return (d3AdTdV2dz + d3AdV3dz*dVdT + d3AdV3*dVdTdz)/d2AdV2/d2AdV2
        - 2.0*(d3AdTdV2 + d3AdV3*dVdT)*d2AdV2dz/d2AdV2/d2AdV2/d2AdV2;
}}

static double {module}_dparam_d3gdp3(double T, double P, int index) {{
    {module}_solve_V(T, P);
    double d3AdV3 = {module}_d3adv3(T,V);
    double d2AdV2 = {module}_d2adv2(T, V);
    double d3AdV3dz = {module}_dparam_d3adv3(T, V, index);
    double d2AdV2dz = {module}_dparam_d2adv2(T, V, index);
    double dVdP = -1.0/d2AdV2;
    double dVdPdz = d2AdV2dz/d2AdV2/d2AdV2;
    /* return d3AdV3*dVdP/d2AdV2/d2AdV2; */
    return d3AdV3dz*dVdP/d2AdV2/d2AdV2 + d3AdV3*dVdPdz/d2AdV2/d2AdV2
        - 2.0*d3AdV3*dVdP*d2AdV2dz/d2AdV2/d2AdV2/d2AdV2;
}}

\
"""

def _create_redundant_function_template_c(model_type='TP'):
    """
    C language implementation of create_redundant_function_template()

    Parameters
    ----------
    model_type: string
        Potential type, either Gibbs free energy ('TP') or Helmholtz free
        energy ('TV')

    Returns
    -------
    string :
        The template string.
    """

    if model_type == 'TP':
        return """\

static double {module}_s(double T, double P) {{
    double result = -{module}_dgdt(T, P);
    return result;
}}

static double {module}_v(double T, double P) {{
    double result = {module}_dgdp(T, P);
    return result;
}}

static double {module}_cv(double T, double P) {{
    double result = -T*{module}_d2gdt2(T, P);
    double dvdt = {module}_d2gdtdp(T, P);
    double dvdp = {module}_d2gdp2(T, P);
    result += T*dvdt*dvdt/dvdp;
    return result;
}}

static double {module}_cp(double T, double P) {{
    double result = -T*{module}_d2gdt2(T, P);
    return result;
}}

static double {module}_dcpdt(double T, double P) {{
    double result = -T*{module}_d3gdt3(T, P) - {module}_d2gdt2(T, P);
    return result;
}}

static double {module}_alpha(double T, double P) {{
    double result = {module}_d2gdtdp(T, P)/{module}_dgdp(T, P);
    return result;
}}

static double {module}_beta(double T, double P) {{
    double result = -{module}_d2gdp2(T, P)/{module}_dgdp(T, P);
    return result;
}}

static double {module}_K(double T, double P) {{
    double result = -{module}_dgdp(T, P)/{module}_d2gdp2(T, P);
    return result;
}}

static double {module}_Kp(double T, double P) {{
    double result = {module}_dgdp(T, P);
    result *= {module}_d3gdp3(T, P);
    result /= pow({module}_d2gdp2(T, P), 2.0);
    return result - 1.0;
}}

\
"""
    elif model_type == 'TV':
        return """\

static double V = {v_initial_guess};
static void {module}_solve_V(double T, double P) {{
    static double Told = 0.0;
    static double Pold = 0.0;
    // Update if *either* T or P changes (was only if both changed)
    // if ((T != Told) && (P != Pold)) {{
    if ((T != Told) || (P != Pold)) {{
        Told = T;
        Pold = P;
        double f = 0.0;
        int iter = 0;
        do {{
            f = -{module}_dadv(T, V) - P;
            double df = -{module}_d2adv2(T, V);
            if (df == 0.0) break;
            V -= f/df;
            if (V <= 0.0) V = 0.001;
            else if (V >= 4.0*{v_initial_guess}) V = 4.0*{v_initial_guess};
            iter++;
        }} while ((fabs(f) > 0.001) && (iter < 200));
    }}
}}

static double {module}_g(double T, double P) {{
    {module}_solve_V(T, P);
    double A = {module}_a(T, V);
    return A + P*V;
}}

static double {module}_dgdt(double T, double P) {{
    {module}_solve_V(T, P);
    double dAdT = {module}_dadt(T, V);
    return dAdT;
}}

static double {module}_dgdp(double T, double P) {{
    {module}_solve_V(T, P);
    return V;
}}

static double {module}_d2gdt2(double T, double P) {{
    {module}_solve_V(T, P);
    double d2AdT2 = {module}_d2adt2(T, V);
    double d2AdTdV = {module}_d2adtdv(T, V);
    double d2AdV2 = {module}_d2adv2(T, V);
    return d2AdT2 - d2AdTdV*d2AdTdV/d2AdV2;
}}

static double {module}_d2gdtdp(double T, double P) {{
    {module}_solve_V(T, P);
    double d2AdTdV = {module}_d2adtdv(T, V);
    double d2AdV2 = {module}_d2adv2(T, V);
    return - d2AdTdV/d2AdV2;
}}

static double {module}_d2gdp2(double T, double P) {{
    {module}_solve_V(T, P);
    double d2AdV2 = {module}_d2adv2(T, V);
    return - 1.0/d2AdV2;
}}

static double {module}_d3gdt3(double T, double P) {{
    {module}_solve_V(T, P);
    double d3AdT3 = {module}_d3adt3(T, V);
    double d3AdT2dV = {module}_d3adt2dv(T, V);
    double d3AdTdV2 = {module}_d3adtdv2(T, V);
    double d3AdV3 = {module}_d3adv3(T,V);
    double d2AdTdV = {module}_d2adtdv(T, V);
    double d2AdV2 = {module}_d2adv2(T, V);
    double dVdT = -d2AdTdV/d2AdV2;
    double d2VdT2 = (-d3AdT2dV - 2.0*d3AdTdV2*dVdT - d3AdV3*dVdT*dVdT)/d2AdV2;
    return d3AdT3 + 2.0*d3AdT2dV*dVdT + d3AdTdV2*dVdT*dVdT + d2AdTdV*d2VdT2;
}}

static double {module}_d3gdt2dp(double T, double P) {{
    {module}_solve_V(T, P);
    double d3AdT2dV = {module}_d3adt2dv(T, V);
    double d3AdTdV2 = {module}_d3adtdv2(T, V);
    double d3AdV3 = {module}_d3adv3(T,V);
    double d2AdTdV = {module}_d2adtdv(T, V);
    double d2AdV2 = {module}_d2adv2(T, V);
    double dVdT = -d2AdTdV/d2AdV2;
    double dVdP = -1.0/d2AdV2;
    double d2VdTdP = (-d3AdTdV2*dVdP - d3AdV3*dVdT*dVdP)/d2AdV2;
    return d3AdT2dV*dVdP + d3AdTdV2*dVdT*dVdP + d2AdTdV*d2VdTdP;
}}

static double {module}_d3gdtdp2(double T, double P) {{
    {module}_solve_V(T, P);
    double d3AdTdV2 = {module}_d3adtdv2(T, V);
    double d3AdV3 = {module}_d3adv3(T,V);
    double d2AdTdV = {module}_d2adtdv(T, V);
    double d2AdV2 = {module}_d2adv2(T, V);
    double dVdT = -d2AdTdV/d2AdV2;
    return (d3AdTdV2 + d3AdV3*dVdT)/d2AdV2/d2AdV2;
}}

static double {module}_d3gdp3(double T, double P) {{
    {module}_solve_V(T, P);
    double d3AdV3 = {module}_d3adv3(T,V);
    double d2AdV2 = {module}_d2adv2(T, V);
    double dVdP = -1.0/d2AdV2;
    return d3AdV3*dVdP/d2AdV2/d2AdV2;
}}

static double {module}_s(double T, double P) {{
    double result = -{module}_dgdt(T, P);
    return result;
}}

static double {module}_v(double T, double P) {{
    double result = {module}_dgdp(T, P);
    return result;
}}

static double {module}_cv(double T, double P) {{
    double result = -T*{module}_d2gdt2(T, P);
    double dvdt = {module}_d2gdtdp(T, P);
    double dvdp = {module}_d2gdp2(T, P);
    result += T*dvdt*dvdt/dvdp;
    return result;
}}

static double {module}_cp(double T, double P) {{
    double result = -T*{module}_d2gdt2(T, P);
    return result;
}}

static double {module}_dcpdt(double T, double P) {{
    double result = -T*{module}_d3gdt3(T, P) - {module}_d2gdt2(T, P);
    return result;
}}

static double {module}_alpha(double T, double P) {{
    double result = {module}_d2gdtdp(T, P)/{module}_dgdp(T, P);
    return result;
}}

static double {module}_beta(double T, double P) {{
    double result = -{module}_d2gdp2(T, P)/{module}_dgdp(T, P);
    return result;
}}

static double {module}_K(double T, double P) {{
    double result = -{module}_dgdp(T, P)/{module}_d2gdp2(T, P);
    return result;
}}

static double {module}_Kp(double T, double P) {{
    double result = {module}_dgdp(T, P);
    result *= {module}_d3gdp3(T, P);
    result /= pow({module}_d2gdp2(T, P), 2.0);
    return result - 1.0;
}}

\
"""
    else:
        print ("Unsupported model_type: ", model_type)
        return ""

###############################
# Simple Solution C templates #
###############################

def _create_soln_calc_template_c():
    """
    C language implementation of create_soln_calc_template()

    Returns
    -------
    string :
        The template string.
    """

    return """\

static double {module}_{func}(double T, double P, double n[{number_components}]) {{
{moles_assign}
    double result;
    {g_code}
    return result;
}}
    \
    """

def _create_soln_calib_code_template_c():
    """
    C language implementation of create_soln_calib_code_template()

    Returns
    -------
    string :
        The template string.
    """

    return """\

static char *identifier = "{git_identifier}";
{code_block_two}

#include "{module}_calc.h"
#include "{module}_calib.h"

const char *{phase}_{module}_calib_identifier(void) {{
    return identifier;
}}

const char *{phase}_{module}_calib_name(void) {{
    return "{phase}";
}}

char *{phase}_{module}_calib_formula(double T, double P, double n[{number_components}]) {{
{code_block_five}
}}

double *{phase}_{module}_calib_conv_elm_to_moles(double *e) {{
{code_block_six}
}}

int {phase}_{module}_calib_test_moles(double *n) {{
{code_block_seven}
}}

const char *{phase}_{module}_calib_endmember_name(int index) {{
    return (*endmember[index].name)();
}}

const char *{phase}_{module}_calib_endmember_formula(int index) {{
    return (*endmember[index].formula)();
}}

const double {phase}_{module}_calib_endmember_mw(int index) {{
    return (*endmember[index].mw)();
}}

const double *{phase}_{module}_calib_endmember_elements(int index) {{
    return (*endmember[index].elements)();
}}

double {phase}_{module}_calib_endmember_mu0(int index, double t, double p) {{
    return (*endmember[index].mu0)(t, p);
}}

double {phase}_{module}_calib_endmember_dmu0dT(int index, double t, double p) {{
    return (*endmember[index].dmu0dT)(t, p);
}}

double {phase}_{module}_calib_endmember_dmu0dP(int index, double t, double p) {{
    return (*endmember[index].dmu0dP)(t, p);
}}

double {phase}_{module}_calib_endmember_d2mu0dT2(int index, double t, double p) {{
    return (*endmember[index].d2mu0dT2)(t, p);
}}

double {phase}_{module}_calib_endmember_d2mu0dTdP(int index, double t, double p) {{
    return (*endmember[index].d2mu0dTdP)(t, p);
}}

double {phase}_{module}_calib_endmember_d2mu0dP2(int index, double t, double p) {{
    return (*endmember[index].d2mu0dP2)(t, p);
}}

double {phase}_{module}_calib_endmember_d3mu0dT3(int index, double t, double p) {{
    return (*endmember[index].d3mu0dT3)(t, p);
}}

double {phase}_{module}_calib_endmember_d3mu0dT2dP(int index, double t, double p) {{
    return (*endmember[index].d3mu0dT2dP)(t, p);
}}

double {phase}_{module}_calib_endmember_d3mu0dTdP2(int index, double t, double p) {{
    return (*endmember[index].d3mu0dTdP2)(t, p);
}}

double {phase}_{module}_calib_endmember_d3mu0dP3(int index, double t, double p) {{
    return (*endmember[index].d3mu0dP3)(t, p);
}}

double {phase}_{module}_calib_g(double T, double P, double n[{number_components}]) {{
    return {module}_g(T, P, n);
}}

void {phase}_{module}_calib_dgdn(double T, double P, double n[{number_components}], double result[{number_components}]) {{
    {module}_dgdn(T, P, n, result);
}}

void {phase}_{module}_calib_d2gdn2(double T, double P, double n[{number_components}], double result[{number_components}]) {{
    {module}_d2gdn2(T, P, n, result);
}}

void {phase}_{module}_calib_d3gdn3(double T, double P, double n[{number_components}], double result[{number_components}]) {{
    {module}_d3gdn3(T, P, n, result);
}}

double {phase}_{module}_calib_dgdt(double T, double P, double n[{number_components}]) {{
    return {module}_dgdt(T, P, n);
}}

void {phase}_{module}_calib_d2gdndt(double T, double P, double n[{number_components}], double result[{number_components}]) {{
    {module}_d2gdndt(T, P, n, result);
}}

void {phase}_{module}_calib_d3gdn2dt(double T, double P, double n[{number_components}], double result[{number_components}]) {{
    {module}_d3gdn2dt(T, P, n, result);
}}

void {phase}_{module}_calib_d4gdn3dt(double T, double P, double n[{number_components}], double result[{number_components}]) {{
    {module}_d4gdn3dt(T, P, n, result);
}}

double {phase}_{module}_calib_dgdp(double T, double P, double n[{number_components}]) {{
    return {module}_dgdp(T, P, n);
}}

void {phase}_{module}_calib_d2gdndp(double T, double P, double n[{number_components}], double result[{number_components}]) {{
    {module}_d2gdndp(T, P, n, result);
}}

void {phase}_{module}_calib_d3gdn2dp(double T, double P, double n[{number_components}], double result[{number_components}]) {{
    {module}_d3gdn2dp(T, P, n, result);
}}

void {phase}_{module}_calib_d4gdn3dp(double T, double P, double n[{number_components}], double result[{number_components}]) {{
    {module}_d4gdn3dp(T, P, n, result);
}}

double {phase}_{module}_calib_d2gdt2(double T, double P, double n[{number_components}]) {{
    return {module}_d2gdt2(T, P, n);
}}

void {phase}_{module}_calib_d3gdndt2(double T, double P, double n[{number_components}], double result[{number_components}]) {{
    {module}_d3gdndt2(T, P, n, result);
}}

void {phase}_{module}_calib_d4gdn2dt2(double T, double P, double n[{number_components}], double result[{number_components}]) {{
    {module}_d4gdn2dt2(T, P, n, result);
}}

void {phase}_{module}_calib_d5gdn3dt2(double T, double P, double n[{number_components}], double result[{number_components}]) {{
    {module}_d5gdn3dt2(T, P, n, result);
}}

double {phase}_{module}_calib_d2gdtdp(double T, double P, double n[{number_components}]) {{
    return {module}_d2gdtdp(T, P, n);
}}

void {phase}_{module}_calib_d3gdndtdp(double T, double P, double n[{number_components}], double result[{number_components}]) {{
    {module}_d3gdndtdp(T, P, n, result);
}}

void {phase}_{module}_calib_d4gdn2dtdp(double T, double P, double n[{number_components}], double result[{number_components}]) {{
    {module}_d4gdn2dtdp(T, P, n, result);
}}

void {phase}_{module}_calib_d5gdn3dtdp(double T, double P, double n[{number_components}], double result[{number_components}]) {{
    {module}_d5gdn3dtdp(T, P, n, result);
}}

double {phase}_{module}_calib_d2gdp2(double T, double P, double n[{number_components}]) {{
    return {module}_d2gdp2(T, P, n);
}}

void {phase}_{module}_calib_d3gdndp2(double T, double P, double n[{number_components}], double result[{number_components}]) {{
    {module}_d3gdndp2(T, P, n, result);
}}

void {phase}_{module}_calib_d4gdn2dp2(double T, double P, double n[{number_components}], double result[{number_components}]) {{
    {module}_d4gdn2dp2(T, P, n, result);
}}

void {phase}_{module}_calib_d5gdn3dp2(double T, double P, double n[{number_components}], double result[{number_components}]) {{
    {module}_d5gdn3dp2(T, P, n, result);
}}

double {phase}_{module}_calib_d3gdt3(double T, double P, double n[{number_components}]) {{
    return {module}_d3gdt3(T, P, n);
}}

void {phase}_{module}_calib_d4gdndt3(double T, double P, double n[{number_components}], double result[{number_components}]) {{
    {module}_d4gdndt3(T, P, n, result);
}}

void {phase}_{module}_calib_d5gdn2dt3(double T, double P, double n[{number_components}], double result[{number_components}]) {{
    {module}_d5gdn2dt3(T, P, n, result);
}}

void {phase}_{module}_calib_d6gdn3dt3(double T, double P, double n[{number_components}], double result[{number_components}]) {{
    {module}_d6gdn3dt3(T, P, n, result);
}}

double {phase}_{module}_calib_d3gdt2dp(double T, double P, double n[{number_components}]) {{
    return {module}_d3gdt2dp(T, P, n);
}}

void {phase}_{module}_calib_d4gdndt2dp(double T, double P, double n[{number_components}], double result[{number_components}]) {{
    {module}_d4gdndt2dp(T, P, n, result);
}}

void {phase}_{module}_calib_d5gdn2dt2dp(double T, double P, double n[{number_components}], double result[{number_components}]) {{
    {module}_d5gdn2dt2dp(T, P, n, result);
}}

void {phase}_{module}_calib_d6gdn3dt2dp(double T, double P, double n[{number_components}], double result[{number_components}]) {{
    {module}_d6gdn3dt2dp(T, P, n, result);
}}

double {phase}_{module}_calib_d3gdtdp2(double T, double P, double n[{number_components}]) {{
    return {module}_d3gdtdp2(T, P, n);
}}

void {phase}_{module}_calib_d4gdndtdp2(double T, double P, double n[{number_components}], double result[{number_components}]) {{
    {module}_d4gdndtdp2(T, P, n, result);
}}

void {phase}_{module}_calib_d5gdn2dtdp2(double T, double P, double n[{number_components}], double result[{number_components}]) {{
    {module}_d5gdn2dtdp2(T, P, n, result);
}}

void {phase}_{module}_calib_d6gdn3dtdp2(double T, double P, double n[{number_components}], double result[{number_components}]) {{
    {module}_d6gdn3dtdp2(T, P, n, result);
}}

double {phase}_{module}_calib_d3gdp3(double T, double P, double n[{number_components}]) {{
    return {module}_d3gdp3(T, P, n);
}}

void {phase}_{module}_calib_d4gdndp3(double T, double P, double n[{number_components}], double result[{number_components}]) {{
    {module}_d4gdndp3(T, P, n, result);
}}

void {phase}_{module}_calib_d5gdn2dp3(double T, double P, double n[{number_components}], double result[{number_components}]) {{
    {module}_d5gdn2dp3(T, P, n, result);
}}

void {phase}_{module}_calib_d6gdn3dp3(double T, double P, double n[{number_components}], double result[{number_components}]) {{
    {module}_d6gdn3dp3(T, P, n, result);
}}

double {phase}_{module}_calib_s(double T, double P, double n[{number_components}]) {{
    return {module}_s(T, P, n);
}}

double {phase}_{module}_calib_v(double T, double P, double n[{number_components}]) {{
    return {module}_v(T, P, n);
}}

double {phase}_{module}_calib_cv(double T, double P, double n[{number_components}]) {{
    return {module}_cv(T, P, n);
}}

double {phase}_{module}_calib_cp(double T, double P, double n[{number_components}]) {{
    return {module}_cp(T, P, n);
}}

double {phase}_{module}_calib_dcpdt(double T, double P, double n[{number_components}]) {{
    return {module}_dcpdt(T, P, n);
}}

double {phase}_{module}_calib_alpha(double T, double P, double n[{number_components}]) {{
    return {module}_alpha(T, P, n);
}}

double {phase}_{module}_calib_beta(double T, double P, double n[{number_components}]) {{
    return {module}_beta(T, P, n);
}}

double {phase}_{module}_calib_K(double T, double P, double n[{number_components}]) {{
    return {module}_K(T, P, n);
}}

double {phase}_{module}_calib_Kp(double T, double P, double n[{number_components}]) {{
    return {module}_Kp(T, P, n);
}}

int {phase}_{module}_get_param_number(void) {{
    return {module}_get_param_number();
}}

const char **{phase}_{module}_get_param_names(void) {{
    return {module}_get_param_names();
}}

const char **{phase}_{module}_get_param_units(void) {{
    return {module}_get_param_units();
}}

void {phase}_{module}_get_param_values(double **values) {{
    {module}_get_param_values(values);
}}

int {phase}_{module}_set_param_values(double *values) {{
    return {module}_set_param_values(values);
}}

double {phase}_{module}_get_param_value(int index) {{
    return {module}_get_param_value(index);
}}

int {phase}_{module}_set_param_value(int index, double value) {{
    return {module}_set_param_value(index, value);
}}

double {phase}_{module}_dparam_g(double T, double P, double n[{number_components}], int index) {{
    return {module}_dparam_g(T, P, n, index);
}}

double {phase}_{module}_dparam_dgdt(double T, double P, double n[{number_components}], int index) {{
    return {module}_dparam_dgdt(T, P, n, index);
}}

double {phase}_{module}_dparam_dgdp(double T, double P, double n[{number_components}], int index) {{
    return {module}_dparam_dgdp(T, P, n, index);
}}

double {phase}_{module}_dparam_d2gdt2(double T, double P, double n[{number_components}], int index) {{
    return {module}_dparam_d2gdt2(T, P, n, index);
}}

double {phase}_{module}_dparam_d2gdtdp(double T, double P, double n[{number_components}], int index) {{
    return {module}_dparam_d2gdtdp(T, P, n, index);
}}

double {phase}_{module}_dparam_d2gdp2(double T, double P, double n[{number_components}], int index) {{
    return {module}_dparam_d2gdp2(T, P, n, index);
}}

double {phase}_{module}_dparam_d3gdt3(double T, double P, double n[{number_components}], int index) {{
    return {module}_dparam_d3gdt3(T, P, n, index);
}}

double {phase}_{module}_dparam_d3gdt2dp(double T, double P, double n[{number_components}], int index) {{
    return {module}_dparam_d3gdt2dp(T, P, n, index);
}}

double {phase}_{module}_dparam_d3gdtdp2(double T, double P, double n[{number_components}], int index) {{
    return {module}_dparam_d3gdtdp2(T, P, n, index);
}}

double {phase}_{module}_dparam_d3gdp3(double T, double P, double n[{number_components}], int index) {{
    return {module}_dparam_d3gdp3(T, P, n, index);
}}

void {phase}_{module}_dparam_dgdn(double T, double P, double n[{number_components}], int index, double result[{number_components}]) {{
    {module}_dparam_dgdn(T, P, n, index, result);
}}

\
"""

def _create_soln_calib_extra_template_c():
    """
    C language implementation of create_soln_calib_extra_template()

    Returns
    -------
    string :
        The template string.
    """

    return """\

static int {module}_get_param_number(void) {{
    return {number_params};
}}

static const char *paramNames[{number_params}] = {names_params};
static const char *paramUnits[{number_params}] = {units_params};

static const char **{module}_get_param_names(void) {{
    return paramNames;
}}

static const char **{module}_get_param_units(void) {{
    return paramUnits;
}}

static void {module}_get_param_values(double **values) {{
{code_block_one}
}}

static int {module}_set_param_values(double *values) {{
{code_block_two}
    return 1;
}}

static double {module}_get_param_value(int index) {{
    double result = 0.0;
    switch (index) {{
{code_block_three}
    default:
        break;
    }}
    return result;
}}

static int {module}_set_param_value(int index, double value) {{
    int result = 1;
    switch (index) {{
{code_block_four}
    default:
        result = 0;
        break;
    }}
    return result;
}}

\
"""

def _create_soln_calib_include_template_c():
    """
    C language implementation of create_soln_fcalib_include_template()

    Returns
    -------
    string :
        The template string.
    """

    return """\

const char *{phase}_{module}_calib_identifier(void);
const char *{phase}_{module}_calib_name(void);
char *{phase}_{module}_calib_formula(double T, double P, double n[{number_components}]);
double *{phase}_{module}_calib_conv_elm_to_moles(double *e);
int {phase}_{module}_calib_test_moles(double *n);

const char *{phase}_{module}_calib_endmember_name(int index);
const char *{phase}_{module}_calib_endmember_formula(int index);
const double {phase}_{module}_calib_endmember_mw(int index);
const double *{phase}_{module}_calib_endmember_elements(int index);
double {phase}_{module}_calib_endmember_mu0(int index, double t, double p);
double {phase}_{module}_calib_endmember_dmu0dT(int index, double t, double p);
double {phase}_{module}_calib_endmember_dmu0dP(int index, double t, double p);
double {phase}_{module}_calib_endmember_d2mu0dT2(int index, double t, double p);
double {phase}_{module}_calib_endmember_d2mu0dTdP(int index, double t, double p);
double {phase}_{module}_calib_endmember_d2mu0dP2(int index, double t, double p);
double {phase}_{module}_calib_endmember_d3mu0dT3(int index, double t, double p);
double {phase}_{module}_calib_endmember_d3mu0dT2dP(int index, double t, double p);
double {phase}_{module}_calib_endmember_d3mu0dTdP2(int index, double t, double p);
double {phase}_{module}_calib_endmember_d3mu0dP3(int index, double t, double p);

double {phase}_{module}_calib_g(double T, double P, double n[{number_components}]);
double {phase}_{module}_calib_dgdt(double T, double P, double n[{number_components}]);
double {phase}_{module}_calib_dgdp(double T, double P, double n[{number_components}]);
double {phase}_{module}_calib_d2gdt2(double T, double P, double n[{number_components}]);
double {phase}_{module}_calib_d2gdtdp(double T, double P, double n[{number_components}]);
double {phase}_{module}_calib_d2gdp2(double T, double P, double n[{number_components}]);
double {phase}_{module}_calib_d3gdt3(double T, double P, double n[{number_components}]);
double {phase}_{module}_calib_d3gdt2dp(double T, double P, double n[{number_components}]);
double {phase}_{module}_calib_d3gdtdp2(double T, double P, double n[{number_components}]);
double {phase}_{module}_calib_d3gdp3(double T, double P, double n[{number_components}]);

void {phase}_{module}_calib_dgdn(double T, double P, double n[{number_components}], double result[{number_components}]);
void {phase}_{module}_calib_d2gdndt(double T, double P, double n[{number_components}], double result[{number_components}]);
void {phase}_{module}_calib_d2gdndp(double T, double P, double n[{number_components}], double result[{number_components}]);
void {phase}_{module}_calib_d3gdndt2(double T, double P, double n[{number_components}], double result[{number_components}]);
void {phase}_{module}_calib_d3gdndtdp(double T, double P, double n[{number_components}], double result[{number_components}]);
void {phase}_{module}_calib_d3gdndp2(double T, double P, double n[{number_components}], double result[{number_components}]);
void {phase}_{module}_calib_d4gdndt3(double T, double P, double n[{number_components}], double result[{number_components}]);
void {phase}_{module}_calib_d4gdndt2dp(double T, double P, double n[{number_components}], double result[{number_components}]);
void {phase}_{module}_calib_d4gdndtdp2(double T, double P, double n[{number_components}], double result[{number_components}]);
void {phase}_{module}_calib_d4gdndp3(double T, double P, double n[{number_components}], double result[{number_components}]);

void {phase}_{module}_calib_d2gdn2(double T, double P, double n[{number_components}], double result[{number_components}]);
void {phase}_{module}_calib_d3gdn2dt(double T, double P, double n[{number_components}], double result[{number_components}]);
void {phase}_{module}_calib_d3gdn2dp(double T, double P, double n[{number_components}], double result[{number_components}]);
void {phase}_{module}_calib_d4gdn2dt2(double T, double P, double n[{number_components}], double result[{number_components}]);
void {phase}_{module}_calib_d4gdn2dtdp(double T, double P, double n[{number_components}], double result[{number_components}]);
void {phase}_{module}_calib_d4gdn2dp2(double T, double P, double n[{number_components}], double result[{number_components}]);
void {phase}_{module}_calib_d5gdn2dt3(double T, double P, double n[{number_components}], double result[{number_components}]);
void {phase}_{module}_calib_d5gdn2dt2dp(double T, double P, double n[{number_components}], double result[{number_components}]);
void {phase}_{module}_calib_d5gdn2dtdp2(double T, double P, double n[{number_components}], double result[{number_components}]);
void {phase}_{module}_calib_d5gdn2dp3(double T, double P, double n[{number_components}], double result[{number_components}]);

void {phase}_{module}_calib_d3gdn3(double T, double P, double n[{number_components}], double result[{number_components}]);
void {phase}_{module}_calib_d4gdn3dt(double T, double P, double n[{number_components}], double result[{number_components}]);
void {phase}_{module}_calib_d4gdn3dp(double T, double P, double n[{number_components}], double result[{number_components}]);
void {phase}_{module}_calib_d5gdn3dt2(double T, double P, double n[{number_components}], double result[{number_components}]);
void {phase}_{module}_calib_d5gdn3dtdp(double T, double P, double n[{number_components}], double result[{number_components}]);
void {phase}_{module}_calib_d5gdn3dp2(double T, double P, double n[{number_components}], double result[{number_components}]);
void {phase}_{module}_calib_d6gdn3dt3(double T, double P, double n[{number_components}], double result[{number_components}]);
void {phase}_{module}_calib_d6gdn3dt2dp(double T, double P, double n[{number_components}], double result[{number_components}]);
void {phase}_{module}_calib_d6gdn3dtdp2(double T, double P, double n[{number_components}], double result[{number_components}]);
void {phase}_{module}_calib_d6gdn3dp3(double T, double P, double n[{number_components}], double result[{number_components}]);

double {phase}_{module}_calib_s(double T, double P, double n[{number_components}]);
double {phase}_{module}_calib_v(double T, double P, double n[{number_components}]);
double {phase}_{module}_calib_cv(double T, double P, double n[{number_components}]);
double {phase}_{module}_calib_cp(double T, double P, double n[{number_components}]);
double {phase}_{module}_calib_dcpdt(double T, double P, double n[{number_components}]);
double {phase}_{module}_calib_alpha(double T, double P, double n[{number_components}]);
double {phase}_{module}_calib_beta(double T, double P, double n[{number_components}]);
double {phase}_{module}_calib_K(double T, double P, double n[{number_components}]);
double {phase}_{module}_calib_Kp(double T, double P, double n[{number_components}]);


int {phase}_{module}_get_param_number(void);
const char **{phase}_{module}_get_param_names(void);
const char **{phase}_{module}_get_param_units(void);
void {phase}_{module}_get_param_values(double **values);
int {phase}_{module}_set_param_values(double *values);
double {phase}_{module}_get_param_value(int index);
int {phase}_{module}_set_param_value(int index, double value);

double {phase}_{module}_dparam_g(double T, double P, double n[{number_components}], int index);
double {phase}_{module}_dparam_dgdt(double T, double P, double n[{number_components}], int index);
double {phase}_{module}_dparam_dgdp(double T, double P, double n[{number_components}], int index);
double {phase}_{module}_dparam_d2gdt2(double T, double P, double n[{number_components}], int index);
double {phase}_{module}_dparam_d2gdtdp(double T, double P, double n[{number_components}], int index);
double {phase}_{module}_dparam_d2gdp2(double T, double P, double n[{number_components}], int index);
double {phase}_{module}_dparam_d3gdt3(double T, double P, double n[{number_components}], int index);
double {phase}_{module}_dparam_d3gdt2dp(double T, double P, double n[{number_components}], int index);
double {phase}_{module}_dparam_d3gdtdp2(double T, double P, double n[{number_components}], int index);
double {phase}_{module}_dparam_d3gdp3(double T, double P, double n[{number_components}], int index);

void {phase}_{module}_dparam_dgdn(double T, double P, double n[{number_components}], int index, double result[{number_components}]);

\
"""

def _create_soln_calib_pyx_template_c():
    """
    C language implementation of create_soln_calib_pyx_template()

    Returns
    -------
    string :
        The template string.
    """

    return """\
# Cython numpy wrapper code for arrays is taken from:
# http://gael-varoquaux.info/programming/cython-example-of-exposing-c-computed-arrays-in-python-without-data-copies.html
# Author: Gael Varoquaux, BSD license

# Declare the prototype of the C functions
cdef extern from "{phase}_{module}_calib.h":
    const char *{phase}_{module}_calib_identifier();
    const char *{phase}_{module}_calib_name();
    char *{phase}_{module}_calib_formula(double T, double P, double n[{number_components}]);
    double *{phase}_{module}_calib_conv_elm_to_moles(double *e);
    int {phase}_{module}_calib_test_moles(double *n);
    const char *{phase}_{module}_calib_endmember_name(int index);
    const char *{phase}_{module}_calib_endmember_formula(int index);
    const double {phase}_{module}_calib_endmember_mw(int index);
    const double *{phase}_{module}_calib_endmember_elements(int index);
    double {phase}_{module}_calib_endmember_mu0(int index, double t, double p);
    double {phase}_{module}_calib_endmember_dmu0dT(int index, double t, double p);
    double {phase}_{module}_calib_endmember_dmu0dP(int index, double t, double p);
    double {phase}_{module}_calib_endmember_d2mu0dT2(int index, double t, double p);
    double {phase}_{module}_calib_endmember_d2mu0dTdP(int index, double t, double p);
    double {phase}_{module}_calib_endmember_d2mu0dP2(int index, double t, double p);
    double {phase}_{module}_calib_endmember_d3mu0dT3(int index, double t, double p);
    double {phase}_{module}_calib_endmember_d3mu0dT2dP(int index, double t, double p);
    double {phase}_{module}_calib_endmember_d3mu0dTdP2(int index, double t, double p);
    double {phase}_{module}_calib_endmember_d3mu0dP3(int index, double t, double p);
    double {phase}_{module}_calib_g(double t, double p, double n[{number_components}])
    double {phase}_{module}_calib_dgdt(double t, double p, double n[{number_components}])
    double {phase}_{module}_calib_dgdp(double t, double p, double n[{number_components}])
    double {phase}_{module}_calib_d2gdt2(double t, double p, double n[{number_components}])
    double {phase}_{module}_calib_d2gdtdp(double t, double p, double n[{number_components}])
    double {phase}_{module}_calib_d2gdp2(double t, double p, double n[{number_components}])
    double {phase}_{module}_calib_d3gdt3(double t, double p, double n[{number_components}])
    double {phase}_{module}_calib_d3gdt2dp(double t, double p, double n[{number_components}])
    double {phase}_{module}_calib_d3gdtdp2(double t, double p, double n[{number_components}])
    double {phase}_{module}_calib_d3gdp3(double t, double p, double n[{number_components}])
    double {phase}_{module}_calib_s(double t, double p, double n[{number_components}])
    double {phase}_{module}_calib_v(double t, double p, double n[{number_components}])
    double {phase}_{module}_calib_cv(double t, double p, double n[{number_components}])
    double {phase}_{module}_calib_cp(double t, double p, double n[{number_components}])
    double {phase}_{module}_calib_dcpdt(double t, double p, double n[{number_components}])
    double {phase}_{module}_calib_alpha(double t, double p, double n[{number_components}])
    double {phase}_{module}_calib_beta(double t, double p, double n[{number_components}])
    double {phase}_{module}_calib_K(double t, double p, double n[{number_components}])
    double {phase}_{module}_calib_Kp(double t, double p, double n[{number_components}])

    void {phase}_{module}_calib_dgdn(double T, double P, double n[{number_components}], double result[{number_components}])
    void {phase}_{module}_calib_d2gdndt(double T, double P, double n[{number_components}], double result[{number_components}])
    void {phase}_{module}_calib_d2gdndp(double T, double P, double n[{number_components}], double result[{number_components}])
    void {phase}_{module}_calib_d3gdndt2(double T, double P, double n[{number_components}], double result[{number_components}])
    void {phase}_{module}_calib_d3gdndtdp(double T, double P, double n[{number_components}], double result[{number_components}])
    void {phase}_{module}_calib_d3gdndp2(double T, double P, double n[{number_components}], double result[{number_components}])
    void {phase}_{module}_calib_d4gdndt3(double T, double P, double n[{number_components}], double result[{number_components}])
    void {phase}_{module}_calib_d4gdndt2dp(double T, double P, double n[{number_components}], double result[{number_components}])
    void {phase}_{module}_calib_d4gdndtdp2(double T, double P, double n[{number_components}], double result[{number_components}])
    void {phase}_{module}_calib_d4gdndp3(double T, double P, double n[{number_components}], double result[{number_components}])

    void {phase}_{module}_calib_d2gdn2(double T, double P, double n[{number_components}], double result[{number_components}])
    void {phase}_{module}_calib_d3gdn2dt(double T, double P, double n[{number_components}], double result[{number_components}])
    void {phase}_{module}_calib_d3gdn2dp(double T, double P, double n[{number_components}], double result[{number_components}])
    void {phase}_{module}_calib_d4gdn2dt2(double T, double P, double n[{number_components}], double result[{number_components}])
    void {phase}_{module}_calib_d4gdn2dtdp(double T, double P, double n[{number_components}], double result[{number_components}])
    void {phase}_{module}_calib_d4gdn2dp2(double T, double P, double n[{number_components}], double result[{number_components}])
    void {phase}_{module}_calib_d5gdn2dt3(double T, double P, double n[{number_components}], double result[{number_components}])
    void {phase}_{module}_calib_d5gdn2dt2dp(double T, double P, double n[{number_components}], double result[{number_components}])
    void {phase}_{module}_calib_d5gdn2dtdp2(double T, double P, double n[{number_components}], double result[{number_components}])
    void {phase}_{module}_calib_d5gdn2dp3(double T, double P, double n[{number_components}], double result[{number_components}])

    void {phase}_{module}_calib_d3gdn3(double T, double P, double n[{number_components}], double result[{number_components}])
    void {phase}_{module}_calib_d4gdn3dt(double T, double P, double n[{number_components}], double result[{number_components}])
    void {phase}_{module}_calib_d4gdn3dp(double T, double P, double n[{number_components}], double result[{number_components}])
    void {phase}_{module}_calib_d5gdn3dt2(double T, double P, double n[{number_components}], double result[{number_components}])
    void {phase}_{module}_calib_d5gdn3dtdp(double T, double P, double n[{number_components}], double result[{number_components}])
    void {phase}_{module}_calib_d5gdn3dp2(double T, double P, double n[{number_components}], double result[{number_components}])
    void {phase}_{module}_calib_d6gdn3dt3(double T, double P, double n[{number_components}], double result[{number_components}])
    void {phase}_{module}_calib_d6gdn3dt2dp(double T, double P, double n[{number_components}], double result[{number_components}])
    void {phase}_{module}_calib_d6gdn3dtdp2(double T, double P, double n[{number_components}], double result[{number_components}])
    void {phase}_{module}_calib_d6gdn3dp3(double T, double P, double n[{number_components}], double result[{number_components}])

    int {phase}_{module}_get_param_number()
    const char **{phase}_{module}_get_param_names()
    const char **{phase}_{module}_get_param_units()
    void {phase}_{module}_get_param_values(double **values)
    int {phase}_{module}_set_param_values(double *values)
    double {phase}_{module}_get_param_value(int index)
    int {phase}_{module}_set_param_value(int index, double value)

    double {phase}_{module}_dparam_g(double t, double p, double n[{number_components}], int index);
    double {phase}_{module}_dparam_dgdt(double t, double p, double n[{number_components}], int index);
    double {phase}_{module}_dparam_dgdp(double t, double p, double n[{number_components}], int index);
    double {phase}_{module}_dparam_d2gdt2(double t, double p, double n[{number_components}], int index);
    double {phase}_{module}_dparam_d2gdtdp(double t, double p, double n[{number_components}], int index);
    double {phase}_{module}_dparam_d2gdp2(double t, double p, double n[{number_components}], int index);
    double {phase}_{module}_dparam_d3gdt3(double t, double p, double n[{number_components}], int index);
    double {phase}_{module}_dparam_d3gdt2dp(double t, double p, double n[{number_components}], int index);
    double {phase}_{module}_dparam_d3gdtdp2(double t, double p, double n[{number_components}], int index);
    double {phase}_{module}_dparam_d3gdp3(double t, double p, double n[{number_components}], int index);

    void {phase}_{module}_dparam_dgdn(double t, double p, double n[{number_components}], int index, double result[{number_components}])

from libc.stdlib cimport malloc, free
from cpython cimport PyObject, Py_INCREF
import ctypes

# Import the Python-level symbols of numpy
import numpy as np

# Import the C-level symbols of numpy
cimport numpy as np

# Numpy must be initialized. When using numpy from C or Cython you must
# _always_ do that, or you will have segfaults
np.import_array()

# here is the "wrapper" signature
def {prefix}_{phase}_{module}_calib_identifier():
    result = <bytes> {phase}_{module}_calib_identifier()
    return result.decode('UTF-8')
def {prefix}_{phase}_{module}_calib_name():
    result = <bytes> {phase}_{module}_calib_name()
    return result.decode('UTF-8')
def {prefix}_{phase}_{module}_calib_formula(double t, double p, np_array):
    cdef double *m = <double *>malloc(len(np_array)*sizeof(double))
    for i in range(np_array.size):
        m[i] = np_array[i]
    r = {phase}_{module}_calib_formula(<double> t, <double> p, <double *> m)
    result = <bytes> r
    free (m)
    free (r)
    return result.decode('UTF-8')
def {prefix}_{phase}_{module}_calib_conv_elm_to_moles(np_array):
    cdef double *e = <double *>malloc(len(np_array)*sizeof(double))
    for i in range(np_array.size):
        e[i] = np_array[i]
    r = {phase}_{module}_calib_conv_elm_to_moles(<double *> e)
    result = []
    for i in range({number_components}):
        result.append(r[i])
    free (e)
    free (r)
    return np.array(result)
def {prefix}_{phase}_{module}_calib_conv_elm_to_tot_moles(np_array):
    cdef double *e = <double *>malloc(len(np_array)*sizeof(double))
    for i in range(np_array.size):
        e[i] = np_array[i]
    r = {phase}_{module}_calib_conv_elm_to_moles(<double *> e)
    result = 0.0
    for i in range({number_components}):
        result += r[i]
    free (e)
    free (r)
    return result
def {prefix}_{phase}_{module}_calib_conv_elm_to_tot_grams(np_array):
    cdef double *e = <double *>malloc(len(np_array)*sizeof(double))
    for i in range(np_array.size):
        e[i] = np_array[i]
    r = {phase}_{module}_calib_conv_elm_to_moles(<double *> e)
    result = 0.0
    for i in range({number_components}):
        mw = {prefix}_{phase}_{module}_calib_endmember_mw(i)
        result += r[i]*mw
    free (e)
    free (r)
    return result
def {prefix}_{phase}_{module}_calib_conv_moles_to_tot_moles(np_array):
    result = 0.0
    for i in range({number_components}):
        result += np_array[i]
    return result
def {prefix}_{phase}_{module}_calib_conv_moles_to_mole_frac(np_array):
    result = np.zeros({number_components})
    sum = np.sum(np_array)
    for i in range({number_components}):
        result[i] += np_array[i]/sum
    return result
def {prefix}_{phase}_{module}_calib_conv_moles_to_elm(np_array):
    result = np.zeros(106)
    for i in range({number_components}):
        end = {prefix}_{phase}_{module}_calib_endmember_elements(i)
        for j in range(0,106):
            result[j] += np_array[i]*end[j]
    return result
def {prefix}_{phase}_{module}_calib_test_moles(np_array):
    cdef double *n = <double *>malloc(len(np_array)*sizeof(double))
    for i in range(np_array.size):
        n[i] = np_array[i]
    result = {phase}_{module}_calib_test_moles(<double *> n)
    free (n)
    return False if result == 0 else True

def {prefix}_{phase}_{module}_calib_endmember_number():
    return {number_components}
def {prefix}_{phase}_{module}_calib_endmember_name(int index):
    assert index in range(0,{number_components}), "index out of range"
    result = {phase}_{module}_calib_endmember_name(index);
    return result.decode('UTF-8')
def {prefix}_{phase}_{module}_calib_endmember_formula(int index):
    assert index in range(0,{number_components}), "index out of range"
    result = {phase}_{module}_calib_endmember_formula(index);
    return result.decode('UTF-8')
def {prefix}_{phase}_{module}_calib_endmember_mw(int index):
    assert index in range(0,{number_components}), "index out of range"
    return {phase}_{module}_calib_endmember_mw(index);
def {prefix}_{phase}_{module}_calib_endmember_elements(int index):
    assert index in range(0,{number_components}), "index out of range"
    r = {phase}_{module}_calib_endmember_elements(index);
    result = []
    for i in range(0,106):
        result.append(r[i])
    return np.array(result)

def {prefix}_{phase}_{module}_calib_species_number():
    return {number_species}
def {prefix}_{phase}_{module}_calib_species_name(int index):
    assert index in range(0,{number_species}), "index out of range"
    result = {phase}_{module}_calib_endmember_name(index);
    return result.decode('UTF-8')
def {prefix}_{phase}_{module}_calib_species_formula(int index):
    assert index in range(0,{number_species}), "index out of range"
    result = {phase}_{module}_calib_endmember_formula(index);
    return result.decode('UTF-8')
def {prefix}_{phase}_{module}_calib_species_mw(int index):
    assert index in range(0,{number_species}), "index out of range"
    return {phase}_{module}_calib_endmember_mw(index);
def {prefix}_{phase}_{module}_calib_species_elements(int index):
    assert index in range(0,{number_species}), "index out of range"
    r = {phase}_{module}_calib_endmember_elements(index);
    result = []
    for i in range(0,106):
        result.append(r[i])
    return np.array(result)

def {prefix}_{phase}_{module}_calib_endmember_mu0(int index, double t, double p):
    assert index in range(0,{number_components}), "index out of range"
    return {phase}_{module}_calib_endmember_mu0(index, <double> t, <double> p);
def {prefix}_{phase}_{module}_calib_endmember_dmu0dT(int index, double t, double p):
    assert index in range(0,{number_components}), "index out of range"
    return {phase}_{module}_calib_endmember_dmu0dT(index, <double> t, <double> p);
def {prefix}_{phase}_{module}_calib_endmember_dmu0dP(int index, double t, double p):
    assert index in range(0,{number_components}), "index out of range"
    return {phase}_{module}_calib_endmember_dmu0dP(index, <double> t, <double> p);
def {prefix}_{phase}_{module}_calib_endmember_d2mu0dT2(int index, double t, double p):
    assert index in range(0,{number_components}), "index out of range"
    return {phase}_{module}_calib_endmember_d2mu0dT2(index, <double> t, <double> p);
def {prefix}_{phase}_{module}_calib_endmember_d2mu0dTdP(int index, double t, double p):
    assert index in range(0,{number_components}), "index out of range"
    return {phase}_{module}_calib_endmember_d2mu0dTdP(index, <double> t, <double> p);
def {prefix}_{phase}_{module}_calib_endmember_d2mu0dP2(int index, double t, double p):
    assert index in range(0,{number_components}), "index out of range"
    return {phase}_{module}_calib_endmember_d2mu0dP2(index, <double> t, <double> p);
def {prefix}_{phase}_{module}_calib_endmember_d3mu0dT3(int index, double t, double p):
    assert index in range(0,{number_components}), "index out of range"
    return {phase}_{module}_calib_endmember_d3mu0dT3(index, <double> t, <double> p);
def {prefix}_{phase}_{module}_calib_endmember_d3mu0dT2dP(int index, double t, double p):
    assert index in range(0,{number_components}), "index out of range"
    return {phase}_{module}_calib_endmember_d3mu0dT2dP(index, <double> t, <double> p);
def {prefix}_{phase}_{module}_calib_endmember_d3mu0dTdP2(int index, double t, double p):
    assert index in range(0,{number_components}), "index out of range"
    return {phase}_{module}_calib_endmember_d3mu0dTdP2(index, <double> t, <double> p);
def {prefix}_{phase}_{module}_calib_endmember_d3mu0dP3(int index, double t, double p):
    assert index in range(0,{number_components}), "index out of range"
    return {phase}_{module}_calib_endmember_d3mu0dP3(index, <double> t, <double> p);

def {prefix}_{phase}_{module}_calib_g(double t, double p, np_array):
    cdef double *m = <double *>malloc(len(np_array)*sizeof(double))
    for i in range(np_array.size):
        m[i] = np_array[i]
    result = {phase}_{module}_calib_g(<double> t, <double> p, <double *> m)
    free (m)
    return result
def {prefix}_{phase}_{module}_calib_dgdt(double t, double p, np_array):
    cdef double *m = <double *>malloc(len(np_array)*sizeof(double))
    for i in range(np_array.size):
        m[i] = np_array[i]
    result = {phase}_{module}_calib_dgdt(<double> t, <double> p, <double *> m)
    free (m)
    return result
def {prefix}_{phase}_{module}_calib_dgdp(double t, double p, np_array):
    cdef double *m = <double *>malloc(len(np_array)*sizeof(double))
    for i in range(np_array.size):
        m[i] = np_array[i]
    result = {phase}_{module}_calib_dgdp(<double> t, <double> p, <double *> m)
    free (m)
    return result
def {prefix}_{phase}_{module}_calib_d2gdt2(double t, double p, np_array):
    cdef double *m = <double *>malloc(len(np_array)*sizeof(double))
    for i in range(np_array.size):
        m[i] = np_array[i]
    result = {phase}_{module}_calib_d2gdt2(<double> t, <double> p, <double *> m)
    free (m)
    return result
def {prefix}_{phase}_{module}_calib_d2gdtdp(double t, double p, np_array):
    cdef double *m = <double *>malloc(len(np_array)*sizeof(double))
    for i in range(np_array.size):
        m[i] = np_array[i]
    result = {phase}_{module}_calib_d2gdtdp(<double> t, <double> p, <double *> m)
    free (m)
    return result
def {prefix}_{phase}_{module}_calib_d2gdp2(double t, double p, np_array):
    cdef double *m = <double *>malloc(len(np_array)*sizeof(double))
    for i in range(np_array.size):
        m[i] = np_array[i]
    result = {phase}_{module}_calib_d2gdp2(<double> t, <double> p, <double *> m)
    free (m)
    return result
def {prefix}_{phase}_{module}_calib_d3gdt3(double t, double p, np_array):
    cdef double *m = <double *>malloc(len(np_array)*sizeof(double))
    for i in range(np_array.size):
        m[i] = np_array[i]
    result = {phase}_{module}_calib_d3gdt3(<double> t, <double> p, <double *> m)
    free (m)
    return result
def {prefix}_{phase}_{module}_calib_d3gdt2dp(double t, double p, np_array):
    cdef double *m = <double *>malloc(len(np_array)*sizeof(double))
    for i in range(np_array.size):
        m[i] = np_array[i]
    result = {phase}_{module}_calib_d3gdt2dp(<double> t, <double> p, <double *> m)
    free (m)
    return result
def {prefix}_{phase}_{module}_calib_d3gdtdp2(double t, double p, np_array):
    cdef double *m = <double *>malloc(len(np_array)*sizeof(double))
    for i in range(np_array.size):
        m[i] = np_array[i]
    result = {phase}_{module}_calib_d3gdtdp2(<double> t, <double> p, <double *> m)
    free (m)
    return result
def {prefix}_{phase}_{module}_calib_d3gdp3(double t, double p, np_array):
    cdef double *m = <double *>malloc(len(np_array)*sizeof(double))
    for i in range(np_array.size):
        m[i] = np_array[i]
    result = {phase}_{module}_calib_d3gdp3(<double> t, <double> p, <double *> m)
    free (m)
    return result
def {prefix}_{phase}_{module}_calib_s(double t, double p, np_array):
    cdef double *m = <double *>malloc(len(np_array)*sizeof(double))
    for i in range(np_array.size):
        m[i] = np_array[i]
    result = {phase}_{module}_calib_s(<double> t, <double> p, <double *> m)
    free (m)
    return result
def {prefix}_{phase}_{module}_calib_v(double t, double p, np_array):
    cdef double *m = <double *>malloc(len(np_array)*sizeof(double))
    for i in range(np_array.size):
        m[i] = np_array[i]
    result = {phase}_{module}_calib_v(<double> t, <double> p, <double *> m)
    free (m)
    return result
def {prefix}_{phase}_{module}_calib_cv(double t, double p, np_array):
    cdef double *m = <double *>malloc(len(np_array)*sizeof(double))
    for i in range(np_array.size):
        m[i] = np_array[i]
    result = {phase}_{module}_calib_cv(<double> t, <double> p, <double *> m)
    free (m)
    return result
def {prefix}_{phase}_{module}_calib_cp(double t, double p, np_array):
    cdef double *m = <double *>malloc(len(np_array)*sizeof(double))
    for i in range(np_array.size):
        m[i] = np_array[i]
    result = {phase}_{module}_calib_cp(<double> t, <double> p, <double *> m)
    free (m)
    return result
def {prefix}_{phase}_{module}_calib_dcpdt(double t, double p, np_array):
    cdef double *m = <double *>malloc(len(np_array)*sizeof(double))
    for i in range(np_array.size):
        m[i] = np_array[i]
    result = {phase}_{module}_calib_dcpdt(<double> t, <double> p, <double *> m)
    free (m)
    return result
def {prefix}_{phase}_{module}_calib_alpha(double t, double p, np_array):
    cdef double *m = <double *>malloc(len(np_array)*sizeof(double))
    for i in range(np_array.size):
        m[i] = np_array[i]
    result = {phase}_{module}_calib_alpha(<double> t, <double> p, <double *> m)
    free (m)
    return result
def {prefix}_{phase}_{module}_calib_beta(double t, double p, np_array):
    cdef double *m = <double *>malloc(len(np_array)*sizeof(double))
    for i in range(np_array.size):
        m[i] = np_array[i]
    result = {phase}_{module}_calib_beta(<double> t, <double> p, <double *> m)
    free (m)
    return result
def {prefix}_{phase}_{module}_calib_K(double t, double p, np_array):
    cdef double *m = <double *>malloc(len(np_array)*sizeof(double))
    for i in range(np_array.size):
        m[i] = np_array[i]
    result = {phase}_{module}_calib_K(<double> t, <double> p, <double *> m)
    free (m)
    return result
def {prefix}_{phase}_{module}_calib_Kp(double t, double p, np_array):
    cdef double *m = <double *>malloc(len(np_array)*sizeof(double))
    for i in range(np_array.size):
        m[i] = np_array[i]
    result = {phase}_{module}_calib_Kp(<double> t, <double> p, <double *> m)
    free (m)
    return result

def {prefix}_{phase}_{module}_calib_dgdn(double t, double p, np_array):
    cdef double *m = <double *>malloc(len(np_array)*sizeof(double))
    for i in range(np_array.size):
        m[i] = np_array[i]
    cdef double *r = <double *>malloc({number_components}*sizeof(double))
    {phase}_{module}_calib_dgdn(<double> t, <double> p, <double *> m, <double *> r)
    result = []
    for i in range({number_components}):
        result.append(r[i])
    free (m)
    free (r)
    return np.array(result)
def {prefix}_{phase}_{module}_calib_d2gdndt(double t, double p, np_array):
    cdef double *m = <double *>malloc(len(np_array)*sizeof(double))
    for i in range(np_array.size):
        m[i] = np_array[i]
    cdef double *r = <double *>malloc({number_components}*sizeof(double))
    {phase}_{module}_calib_d2gdndt(<double> t, <double> p, <double *> m, <double *> r)
    result = []
    for i in range({number_components}):
        result.append(r[i])
    free (m)
    free (r)
    return np.array(result)
def {prefix}_{phase}_{module}_calib_d2gdndp(double t, double p, np_array):
    cdef double *m = <double *>malloc(len(np_array)*sizeof(double))
    for i in range(np_array.size):
        m[i] = np_array[i]
    cdef double *r = <double *>malloc({number_components}*sizeof(double))
    {phase}_{module}_calib_d2gdndp(<double> t, <double> p, <double *> m, <double *> r)
    result = []
    for i in range({number_components}):
        result.append(r[i])
    free (m)
    free (r)
    return np.array(result)
def {prefix}_{phase}_{module}_calib_d3gdndt2(double t, double p, np_array):
    cdef double *m = <double *>malloc(len(np_array)*sizeof(double))
    for i in range(np_array.size):
        m[i] = np_array[i]
    cdef double *r = <double *>malloc({number_components}*sizeof(double))
    {phase}_{module}_calib_d3gdndt2(<double> t, <double> p, <double *> m, <double *> r)
    result = []
    for i in range({number_components}):
        result.append(r[i])
    free (m)
    free (r)
    return np.array(result)
def {prefix}_{phase}_{module}_calib_d3gdndtdp(double t, double p, np_array):
    cdef double *m = <double *>malloc(len(np_array)*sizeof(double))
    for i in range(np_array.size):
        m[i] = np_array[i]
    cdef double *r = <double *>malloc({number_components}*sizeof(double))
    {phase}_{module}_calib_d3gdndtdp(<double> t, <double> p, <double *> m, <double *> r)
    result = []
    for i in range({number_components}):
        result.append(r[i])
    free (m)
    free (r)
    return np.array(result)
def {prefix}_{phase}_{module}_calib_d3gdndp2(double t, double p, np_array):
    cdef double *m = <double *>malloc(len(np_array)*sizeof(double))
    for i in range(np_array.size):
        m[i] = np_array[i]
    cdef double *r = <double *>malloc({number_components}*sizeof(double))
    {phase}_{module}_calib_d3gdndp2(<double> t, <double> p, <double *> m, <double *> r)
    result = []
    for i in range({number_components}):
        result.append(r[i])
    free (m)
    free (r)
    return np.array(result)
def {prefix}_{phase}_{module}_calib_d4gdndt3(double t, double p, np_array):
    cdef double *m = <double *>malloc(len(np_array)*sizeof(double))
    for i in range(np_array.size):
        m[i] = np_array[i]
    cdef double *r = <double *>malloc({number_components}*sizeof(double))
    {phase}_{module}_calib_d4gdndt3(<double> t, <double> p, <double *> m, <double *> r)
    result = []
    for i in range({number_components}):
        result.append(r[i])
    free (m)
    free (r)
    return np.array(result)
def {prefix}_{phase}_{module}_calib_d4gdndt2dp(double t, double p, np_array):
    cdef double *m = <double *>malloc(len(np_array)*sizeof(double))
    for i in range(np_array.size):
        m[i] = np_array[i]
    cdef double *r = <double *>malloc({number_components}*sizeof(double))
    {phase}_{module}_calib_d4gdndt2dp(<double> t, <double> p, <double *> m, <double *> r)
    result = []
    for i in range({number_components}):
        result.append(r[i])
    free (m)
    free (r)
    return np.array(result)
def {prefix}_{phase}_{module}_calib_d4gdndtdp2(double t, double p, np_array):
    cdef double *m = <double *>malloc(len(np_array)*sizeof(double))
    for i in range(np_array.size):
        m[i] = np_array[i]
    cdef double *r = <double *>malloc({number_components}*sizeof(double))
    {phase}_{module}_calib_d4gdndtdp2(<double> t, <double> p, <double *> m, <double *> r)
    result = []
    for i in range({number_components}):
        result.append(r[i])
    free (m)
    free (r)
    return np.array(result)
def {prefix}_{phase}_{module}_calib_d4gdndp3(double t, double p, np_array):
    cdef double *m = <double *>malloc(len(np_array)*sizeof(double))
    for i in range(np_array.size):
        m[i] = np_array[i]
    cdef double *r = <double *>malloc({number_components}*sizeof(double))
    {phase}_{module}_calib_d4gdndp3(<double> t, <double> p, <double *> m, <double *> r)
    result = []
    for i in range({number_components}):
        result.append(r[i])
    free (m)
    free (r)
    return np.array(result)

def {prefix}_{phase}_{module}_calib_d2gdn2(double t, double p, np_array):
    cdef double *m = <double *>malloc(len(np_array)*sizeof(double))
    for i in range(np_array.size):
        m[i] = np_array[i]
    ndim = int({number_components}*({number_components}-1)/2 + {number_components})
    cdef double *r = <double *>malloc(ndim*sizeof(double))
    {phase}_{module}_calib_d2gdn2(<double> t, <double> p, <double *> m, <double *> r)
    result = []
    for i in range(ndim):
        result.append(r[i])
    free (m)
    free (r)
    return np.array(result)
def {prefix}_{phase}_{module}_calib_d3gdn2dt(double t, double p, np_array):
    cdef double *m = <double *>malloc(len(np_array)*sizeof(double))
    for i in range(np_array.size):
        m[i] = np_array[i]
    ndim = int({number_components}*({number_components}-1)/2 + {number_components})
    cdef double *r = <double *>malloc(ndim*sizeof(double))
    {phase}_{module}_calib_d3gdn2dt(<double> t, <double> p, <double *> m, <double *> r)
    result = []
    for i in range(ndim):
        result.append(r[i])
    free (m)
    free (r)
    return np.array(result)
def {prefix}_{phase}_{module}_calib_d3gdn2dp(double t, double p, np_array):
    cdef double *m = <double *>malloc(len(np_array)*sizeof(double))
    for i in range(np_array.size):
        m[i] = np_array[i]
    ndim = int({number_components}*({number_components}-1)/2 + {number_components})
    cdef double *r = <double *>malloc(ndim*sizeof(double))
    {phase}_{module}_calib_d3gdn2dp(<double> t, <double> p, <double *> m, <double *> r)
    result = []
    for i in range(ndim):
        result.append(r[i])
    free (m)
    free (r)
    return np.array(result)
def {prefix}_{phase}_{module}_calib_d4gdn2dt2(double t, double p, np_array):
    cdef double *m = <double *>malloc(len(np_array)*sizeof(double))
    for i in range(np_array.size):
        m[i] = np_array[i]
    ndim = int({number_components}*({number_components}-1)/2 + {number_components})
    cdef double *r = <double *>malloc(ndim*sizeof(double))
    {phase}_{module}_calib_d4gdn2dt2(<double> t, <double> p, <double *> m, <double *> r)
    result = []
    for i in range(ndim):
        result.append(r[i])
    free (m)
    free (r)
    return np.array(result)
def {prefix}_{phase}_{module}_calib_d4gdn2dtdp(double t, double p, np_array):
    cdef double *m = <double *>malloc(len(np_array)*sizeof(double))
    for i in range(np_array.size):
        m[i] = np_array[i]
    ndim = int({number_components}*({number_components}-1)/2 + {number_components})
    cdef double *r = <double *>malloc(ndim*sizeof(double))
    {phase}_{module}_calib_d4gdn2dtdp(<double> t, <double> p, <double *> m, <double *> r)
    result = []
    for i in range(ndim):
        result.append(r[i])
    free (m)
    free (r)
    return np.array(result)
def {prefix}_{phase}_{module}_calib_d4gdn2dp2(double t, double p, np_array):
    cdef double *m = <double *>malloc(len(np_array)*sizeof(double))
    for i in range(np_array.size):
        m[i] = np_array[i]
    ndim = int({number_components}*({number_components}-1)/2 + {number_components})
    cdef double *r = <double *>malloc(ndim*sizeof(double))
    {phase}_{module}_calib_d4gdn2dp2(<double> t, <double> p, <double *> m, <double *> r)
    result = []
    for i in range(ndim):
        result.append(r[i])
    free (m)
    free (r)
    return np.array(result)
def {prefix}_{phase}_{module}_calib_d5gdn2dt3(double t, double p, np_array):
    cdef double *m = <double *>malloc(len(np_array)*sizeof(double))
    for i in range(np_array.size):
        m[i] = np_array[i]
    ndim = int({number_components}*({number_components}-1)/2 + {number_components})
    cdef double *r = <double *>malloc(ndim*sizeof(double))
    {phase}_{module}_calib_d5gdn2dt3(<double> t, <double> p, <double *> m, <double *> r)
    result = []
    for i in range(ndim):
        result.append(r[i])
    free (m)
    free (r)
    return np.array(result)
def {prefix}_{phase}_{module}_calib_d5gdn2dt2dp(double t, double p, np_array):
    cdef double *m = <double *>malloc(len(np_array)*sizeof(double))
    for i in range(np_array.size):
        m[i] = np_array[i]
    ndim = int({number_components}*({number_components}-1)/2 + {number_components})
    cdef double *r = <double *>malloc(ndim*sizeof(double))
    {phase}_{module}_calib_d5gdn2dt2dp(<double> t, <double> p, <double *> m, <double *> r)
    result = []
    for i in range(ndim):
        result.append(r[i])
    free (m)
    free (r)
    return np.array(result)
def {prefix}_{phase}_{module}_calib_d5gdn2dtdp2(double t, double p, np_array):
    cdef double *m = <double *>malloc(len(np_array)*sizeof(double))
    for i in range(np_array.size):
        m[i] = np_array[i]
    ndim = int({number_components}*({number_components}-1)/2 + {number_components})
    cdef double *r = <double *>malloc(ndim*sizeof(double))
    {phase}_{module}_calib_d5gdn2dtdp2(<double> t, <double> p, <double *> m, <double *> r)
    result = []
    for i in range(ndim):
        result.append(r[i])
    free (m)
    free (r)
    return np.array(result)
def {prefix}_{phase}_{module}_calib_d5gdn2dp3(double t, double p, np_array):
    cdef double *m = <double *>malloc(len(np_array)*sizeof(double))
    for i in range(np_array.size):
        m[i] = np_array[i]
    ndim = int({number_components}*({number_components}-1)/2 + {number_components})
    cdef double *r = <double *>malloc(ndim*sizeof(double))
    {phase}_{module}_calib_d5gdn2dp3(<double> t, <double> p, <double *> m, <double *> r)
    result = []
    for i in range(ndim):
        result.append(r[i])
    free (m)
    free (r)
    return np.array(result)

def {prefix}_{phase}_{module}_calib_d3gdn3(double t, double p, np_array):
    cdef double *m = <double *>malloc(len(np_array)*sizeof(double))
    for i in range(np_array.size):
        m[i] = np_array[i]
    nc = {number_components}
    ndim = int(nc*(nc+1)*(nc+2)/6)
    cdef double *r = <double *>malloc(ndim*sizeof(double))
    {phase}_{module}_calib_d3gdn3(<double> t, <double> p, <double *> m, <double *> r)
    result = []
    for i in range(ndim):
        result.append(r[i])
    free (m)
    free (r)
    return np.array(result)
def {prefix}_{phase}_{module}_calib_d4gdn3dt(double t, double p, np_array):
    cdef double *m = <double *>malloc(len(np_array)*sizeof(double))
    for i in range(np_array.size):
        m[i] = np_array[i]
    nc = {number_components}
    ndim = int(nc*(nc+1)*(nc+2)/6)
    cdef double *r = <double *>malloc(ndim*sizeof(double))
    {phase}_{module}_calib_d4gdn3dt(<double> t, <double> p, <double *> m, <double *> r)
    result = []
    for i in range(ndim):
        result.append(r[i])
    free (m)
    free (r)
    return np.array(result)
def {prefix}_{phase}_{module}_calib_d4gdn3dp(double t, double p, np_array):
    cdef double *m = <double *>malloc(len(np_array)*sizeof(double))
    for i in range(np_array.size):
        m[i] = np_array[i]
    nc = {number_components}
    ndim = int(nc*(nc+1)*(nc+2)/6)
    cdef double *r = <double *>malloc(ndim*sizeof(double))
    {phase}_{module}_calib_d4gdn3dp(<double> t, <double> p, <double *> m, <double *> r)
    result = []
    for i in range(ndim):
        result.append(r[i])
    free (m)
    free (r)
    return np.array(result)
def {prefix}_{phase}_{module}_calib_d5gdn3dt2(double t, double p, np_array):
    cdef double *m = <double *>malloc(len(np_array)*sizeof(double))
    for i in range(np_array.size):
        m[i] = np_array[i]
    nc = {number_components}
    ndim = int(nc*(nc+1)*(nc+2)/6)
    cdef double *r = <double *>malloc(ndim*sizeof(double))
    {phase}_{module}_calib_d5gdn3dt2(<double> t, <double> p, <double *> m, <double *> r)
    result = []
    for i in range(ndim):
        result.append(r[i])
    free (m)
    free (r)
    return np.array(result)
def {prefix}_{phase}_{module}_calib_d5gdn3dtdp(double t, double p, np_array):
    cdef double *m = <double *>malloc(len(np_array)*sizeof(double))
    for i in range(np_array.size):
        m[i] = np_array[i]
    nc = {number_components}
    ndim = int(nc*(nc+1)*(nc+2)/6)
    cdef double *r = <double *>malloc(ndim*sizeof(double))
    {phase}_{module}_calib_d5gdn3dtdp(<double> t, <double> p, <double *> m, <double *> r)
    result = []
    for i in range(ndim):
        result.append(r[i])
    free (m)
    free (r)
    return np.array(result)
def {prefix}_{phase}_{module}_calib_d5gdn3dp2(double t, double p, np_array):
    cdef double *m = <double *>malloc(len(np_array)*sizeof(double))
    for i in range(np_array.size):
        m[i] = np_array[i]
    nc = {number_components}
    ndim = int(nc*(nc+1)*(nc+2)/6)
    cdef double *r = <double *>malloc(ndim*sizeof(double))
    {phase}_{module}_calib_d5gdn3dp2(<double> t, <double> p, <double *> m, <double *> r)
    result = []
    for i in range(ndim):
        result.append(r[i])
    free (m)
    free (r)
    return np.array(result)
def {prefix}_{phase}_{module}_calib_d6gdn3dt3(double t, double p, np_array):
    cdef double *m = <double *>malloc(len(np_array)*sizeof(double))
    for i in range(np_array.size):
        m[i] = np_array[i]
    nc = {number_components}
    ndim = int(nc*(nc+1)*(nc+2)/6)
    cdef double *r = <double *>malloc(ndim*sizeof(double))
    {phase}_{module}_calib_d6gdn3dt3(<double> t, <double> p, <double *> m, <double *> r)
    result = []
    for i in range(ndim):
        result.append(r[i])
    free (m)
    free (r)
    return np.array(result)
def {prefix}_{phase}_{module}_calib_d6gdn3dt2dp(double t, double p, np_array):
    cdef double *m = <double *>malloc(len(np_array)*sizeof(double))
    for i in range(np_array.size):
        m[i] = np_array[i]
    nc = {number_components}
    ndim = int(nc*(nc+1)*(nc+2)/6)
    cdef double *r = <double *>malloc(ndim*sizeof(double))
    {phase}_{module}_calib_d6gdn3dt2dp(<double> t, <double> p, <double *> m, <double *> r)
    result = []
    for i in range(ndim):
        result.append(r[i])
    free (m)
    free (r)
    return np.array(result)
def {prefix}_{phase}_{module}_calib_d6gdn3dtdp2(double t, double p, np_array):
    cdef double *m = <double *>malloc(len(np_array)*sizeof(double))
    for i in range(np_array.size):
        m[i] = np_array[i]
    nc = {number_components}
    ndim = int(nc*(nc+1)*(nc+2)/6)
    cdef double *r = <double *>malloc(ndim*sizeof(double))
    {phase}_{module}_calib_d6gdn3dtdp2(<double> t, <double> p, <double *> m, <double *> r)
    result = []
    for i in range(ndim):
        result.append(r[i])
    free (m)
    free (r)
    return np.array(result)
def {prefix}_{phase}_{module}_calib_d6gdn3dp3(double t, double p, np_array):
    cdef double *m = <double *>malloc(len(np_array)*sizeof(double))
    for i in range(np_array.size):
        m[i] = np_array[i]
    nc = {number_components}
    ndim = int(nc*(nc+1)*(nc+2)/6)
    cdef double *r = <double *>malloc(ndim*sizeof(double))
    {phase}_{module}_calib_d6gdn3dp3(<double> t, <double> p, <double *> m, <double *> r)
    result = []
    for i in range(ndim):
        result.append(r[i])
    free (m)
    free (r)
    return np.array(result)


def {prefix}_{phase}_{module}_get_param_number():
    result = {phase}_{module}_get_param_number()
    return result
def {prefix}_{phase}_{module}_get_param_names():
    cdef const char **names = {phase}_{module}_get_param_names()
    n = {phase}_{module}_get_param_number()
    result = []
    for i in range(0,n):
        entry = <bytes> names[i]
        result.append(entry.decode('UTF-8'))
    return result
def {prefix}_{phase}_{module}_get_param_units():
    cdef const char **units = {phase}_{module}_get_param_units()
    n = {phase}_{module}_get_param_number()
    result = []
    for i in range(0,n):
        entry = <bytes> units[i]
        result.append(entry.decode('UTF-8'))
    return result
def {prefix}_{phase}_{module}_get_param_values():
    n = {phase}_{module}_get_param_number()
    cdef double *m = <double *>malloc(n*sizeof(double))
    {phase}_{module}_get_param_values(&m)
    np_array = np.zeros(n)
    for i in range(n):
        np_array[i] = m[i]
    free(m)
    return np_array
def {prefix}_{phase}_{module}_set_param_values(np_array):
    n = len(np_array)
    cdef double *m = <double *>malloc(n*sizeof(double))
    for i in range(np_array.size):
        m[i] = np_array[i]
    result = {phase}_{module}_set_param_values(m);
    free(m)
    return result
def {prefix}_{phase}_{module}_get_param_value(int index):
    result = {phase}_{module}_get_param_value(<int> index)
    return result
def {prefix}_{phase}_{module}_set_param_value(int index, double value):
    result = {phase}_{module}_set_param_value(<int> index, <double> value)
    return result

def {prefix}_{phase}_{module}_dparam_g(double t, double p, np_array, int index):
    cdef double *m = <double *>malloc(len(np_array)*sizeof(double))
    for i in range(np_array.size):
        m[i] = np_array[i]
    result = {phase}_{module}_dparam_g(<double> t, <double> p, <double *> m, <int> index)
    free(m)
    return result
def {prefix}_{phase}_{module}_dparam_dgdt(double t, double p, np_array, int index):
    cdef double *m = <double *>malloc(len(np_array)*sizeof(double))
    for i in range(np_array.size):
        m[i] = np_array[i]
    result = {phase}_{module}_dparam_dgdt(<double> t, <double> p, <double *> m, <int> index)
    free(m)
    return result
def {prefix}_{phase}_{module}_dparam_dgdp(double t, double p, np_array, int index):
    cdef double *m = <double *>malloc(len(np_array)*sizeof(double))
    for i in range(np_array.size):
        m[i] = np_array[i]
    result = {phase}_{module}_dparam_dgdp(<double> t, <double> p, <double *> m, <int> index)
    free(m)
    return result
def {prefix}_{phase}_{module}_dparam_d2gdt2(double t, double p, np_array, int index):
    cdef double *m = <double *>malloc(len(np_array)*sizeof(double))
    for i in range(np_array.size):
        m[i] = np_array[i]
    result = {phase}_{module}_dparam_d2gdt2(<double> t, <double> p, <double *> m, <int> index)
    free(m)
    return result
def {prefix}_{phase}_{module}_dparam_d2gdtdp(double t, double p, np_array, int index):
    cdef double *m = <double *>malloc(len(np_array)*sizeof(double))
    for i in range(np_array.size):
        m[i] = np_array[i]
    result = {phase}_{module}_dparam_d2gdtdp(<double> t, <double> p, <double *> m, <int> index)
    free(m)
    return result
def {prefix}_{phase}_{module}_dparam_d2gdp2(double t, double p, np_array, int index):
    cdef double *m = <double *>malloc(len(np_array)*sizeof(double))
    for i in range(np_array.size):
        m[i] = np_array[i]
    result = {phase}_{module}_dparam_d2gdp2(<double> t, <double> p, <double *> m, <int> index)
    free(m)
    return result
def {prefix}_{phase}_{module}_dparam_d3gdt3(double t, double p, np_array, int index):
    cdef double *m = <double *>malloc(len(np_array)*sizeof(double))
    for i in range(np_array.size):
        m[i] = np_array[i]
    result = {phase}_{module}_dparam_d3gdt3(<double> t, <double> p, <double *> m, <int> index)
    free(m)
    return result
def {prefix}_{phase}_{module}_dparam_d3gdt2dp(double t, double p, np_array, int index):
    cdef double *m = <double *>malloc(len(np_array)*sizeof(double))
    for i in range(np_array.size):
        m[i] = np_array[i]
    result = {phase}_{module}_dparam_d3gdt2dp(<double> t, <double> p, <double *> m, <int> index)
    free(m)
    return result
def {prefix}_{phase}_{module}_dparam_d3gdtdp2(double t, double p, np_array, int index):
    cdef double *m = <double *>malloc(len(np_array)*sizeof(double))
    for i in range(np_array.size):
        m[i] = np_array[i]
    result = {phase}_{module}_dparam_d3gdtdp2(<double> t, <double> p, <double *> m, <int> index)
    free(m)
    return result
def {prefix}_{phase}_{module}_dparam_d3gdp3(double t, double p, np_array, int index):
    cdef double *m = <double *>malloc(len(np_array)*sizeof(double))
    for i in range(np_array.size):
        m[i] = np_array[i]
    result = {phase}_{module}_dparam_d3gdp3(<double> t, <double> p, <double *> m, <int> index)
    free(m)
    return result

def {prefix}_{phase}_{module}_dparam_dgdn(double t, double p, np_array, int index):
    cdef double *m = <double *>malloc(len(np_array)*sizeof(double))
    cdef double *r = <double *>malloc(len(np_array)*sizeof(double))
    for i in range(np_array.size):
        m[i] = np_array[i]
    {phase}_{module}_dparam_dgdn(<double> t, <double> p, <double *> m, <int> index, <double *> r)
    r_np_array = np.zeros(len(np_array))
    for i in range(r_np_array.size):
        r_np_array[i] = r[i]
    free(m)
    free(r)
    return r_np_array

\
"""

def _create_soln_calib_template_c():
    """
    C language implementation of create_soln_calib_template()

    Returns
    -------
    string :
        The template string.
    """

    return """\

static double {module}_dparam_{func}(double T, double P, double n[{number_components}], int index) {{
{moles_assign}
    double result = 0.0;
    switch (index) {{
{switch_code}
    default:
        break;
    }}
        return result;
}}
\
"""

def _create_soln_deriv_template_c():
    """
    C language implementation of create_soln_deriv_template()

    Returns
    -------
    string :
        The template string.
    """

    return """\

static void {module}_{func}(double T, double P, double n[{number_components}], double result[{number_components}]) {{
{moles_assign}
{derivative_code}
}}
    \
    """

def _create_soln_fast_code_template_c():
    """
    C language implementation of create_soln_fast_code_template()

    Returns
    -------
    string :
        The template string.
    """

    return """\

static const char *identifier = "{git_identifier}";
{code_block_one}

#include "{module}_calc.h"

const char *{phase}_{module}_identifier(void) {{
    return identifier;
}}

const char *{phase}_{module}_name(void) {{
    return "{phase}";
}}

char *{phase}_{module}_formula(double T, double P, double n[{number_components}]) {{
{code_block_five}
}}

double *{phase}_{module}_conv_elm_to_moles(double *e) {{
{code_block_six}
}}

int {phase}_{module}_test_moles(double *n) {{
{code_block_seven}
}}

const char *{phase}_{module}_endmember_name(int index) {{
    return (*endmember[index].name)();
}}

const char *{phase}_{module}_endmember_formula(int index) {{
    return (*endmember[index].formula)();
}}

const double {phase}_{module}_endmember_mw(int index) {{
    return (*endmember[index].mw)();
}}

const double *{phase}_{module}_endmember_elements(int index) {{
    return (*endmember[index].elements)();
}}

double {phase}_{module}_endmember_mu0(int index, double t, double p) {{
    return (*endmember[index].mu0)(t, p);
}}

double {phase}_{module}_endmember_dmu0dT(int index, double t, double p) {{
    return (*endmember[index].dmu0dT)(t, p);
}}

double {phase}_{module}_endmember_dmu0dP(int index, double t, double p) {{
    return (*endmember[index].dmu0dP)(t, p);
}}

double {phase}_{module}_endmember_d2mu0dT2(int index, double t, double p) {{
    return (*endmember[index].d2mu0dT2)(t, p);
}}

double {phase}_{module}_endmember_d2mu0dTdP(int index, double t, double p) {{
    return (*endmember[index].d2mu0dTdP)(t, p);
}}

double {phase}_{module}_endmember_d2mu0dP2(int index, double t, double p) {{
    return (*endmember[index].d2mu0dP2)(t, p);
}}

double {phase}_{module}_endmember_d3mu0dT3(int index, double t, double p) {{
    return (*endmember[index].d3mu0dT3)(t, p);
}}

double {phase}_{module}_endmember_d3mu0dT2dP(int index, double t, double p) {{
    return (*endmember[index].d3mu0dT2dP)(t, p);
}}

double {phase}_{module}_endmember_d3mu0dTdP2(int index, double t, double p) {{
    return (*endmember[index].d3mu0dTdP2)(t, p);
}}

double {phase}_{module}_endmember_d3mu0dP3(int index, double t, double p) {{
    return (*endmember[index].d3mu0dP3)(t, p);
}}


double {phase}_{module}_g(double T, double P, double n[{number_components}]) {{
    return {module}_g(T, P, n);
}}

void {phase}_{module}_dgdn(double T, double P, double n[{number_components}], double result[{number_components}]) {{
    {module}_dgdn(T, P, n, result);
}}

void {phase}_{module}_d2gdn2(double T, double P, double n[{number_components}], double result[{number_components}]) {{
    {module}_d2gdn2(T, P, n, result);
}}

void {phase}_{module}_d3gdn3(double T, double P, double n[{number_components}], double result[{number_components}]) {{
    {module}_d3gdn3(T, P, n, result);
}}

double {phase}_{module}_dgdt(double T, double P, double n[{number_components}]) {{
    return {module}_dgdt(T, P, n);
}}

void {phase}_{module}_d2gdndt(double T, double P, double n[{number_components}], double result[{number_components}]) {{
    {module}_d2gdndt(T, P, n, result);
}}

void {phase}_{module}_d3gdn2dt(double T, double P, double n[{number_components}], double result[{number_components}]) {{
    {module}_d3gdn2dt(T, P, n, result);
}}

void {phase}_{module}_d4gdn3dt(double T, double P, double n[{number_components}], double result[{number_components}]) {{
    {module}_d4gdn3dt(T, P, n, result);
}}

double {phase}_{module}_dgdp(double T, double P, double n[{number_components}]) {{
    return {module}_dgdp(T, P, n);
}}

void {phase}_{module}_d2gdndp(double T, double P, double n[{number_components}], double result[{number_components}]) {{
    {module}_d2gdndp(T, P, n, result);
}}

void {phase}_{module}_d3gdn2dp(double T, double P, double n[{number_components}], double result[{number_components}]) {{
    {module}_d3gdn2dp(T, P, n, result);
}}

void {phase}_{module}_d4gdn3dp(double T, double P, double n[{number_components}], double result[{number_components}]) {{
    {module}_d4gdn3dp(T, P, n, result);
}}

double {phase}_{module}_d2gdt2(double T, double P, double n[{number_components}]) {{
    return {module}_d2gdt2(T, P, n);
}}

void {phase}_{module}_d3gdndt2(double T, double P, double n[{number_components}], double result[{number_components}]) {{
    {module}_d3gdndt2(T, P, n, result);
}}

void {phase}_{module}_d4gdn2dt2(double T, double P, double n[{number_components}], double result[{number_components}]) {{
    {module}_d4gdn2dt2(T, P, n, result);
}}

void {phase}_{module}_d5gdn3dt2(double T, double P, double n[{number_components}], double result[{number_components}]) {{
    {module}_d5gdn3dt2(T, P, n, result);
}}

double {phase}_{module}_d2gdtdp(double T, double P, double n[{number_components}]) {{
    return {module}_d2gdtdp(T, P, n);
}}

void {phase}_{module}_d3gdndtdp(double T, double P, double n[{number_components}], double result[{number_components}]) {{
    {module}_d3gdndtdp(T, P, n, result);
}}

void {phase}_{module}_d4gdn2dtdp(double T, double P, double n[{number_components}], double result[{number_components}]) {{
    {module}_d4gdn2dtdp(T, P, n, result);
}}

void {phase}_{module}_d5gdn3dtdp(double T, double P, double n[{number_components}], double result[{number_components}]) {{
    {module}_d5gdn3dtdp(T, P, n, result);
}}

double {phase}_{module}_d2gdp2(double T, double P, double n[{number_components}]) {{
    return {module}_d2gdp2(T, P, n);
}}

void {phase}_{module}_d3gdndp2(double T, double P, double n[{number_components}], double result[{number_components}]) {{
    {module}_d3gdndp2(T, P, n, result);
}}

void {phase}_{module}_d4gdn2dp2(double T, double P, double n[{number_components}], double result[{number_components}]) {{
    {module}_d4gdn2dp2(T, P, n, result);
}}

void {phase}_{module}_d5gdn3dp2(double T, double P, double n[{number_components}], double result[{number_components}]) {{
    {module}_d5gdn3dp2(T, P, n, result);
}}

double {phase}_{module}_d3gdt3(double T, double P, double n[{number_components}]) {{
    return {module}_d3gdt3(T, P, n);
}}

void {phase}_{module}_d4gdndt3(double T, double P, double n[{number_components}], double result[{number_components}]) {{
    {module}_d4gdndt3(T, P, n, result);
}}

void {phase}_{module}_d5gdn2dt3(double T, double P, double n[{number_components}], double result[{number_components}]) {{
    {module}_d5gdn2dt3(T, P, n, result);
}}

void {phase}_{module}_d6gdn3dt3(double T, double P, double n[{number_components}], double result[{number_components}]) {{
    {module}_d6gdn3dt3(T, P, n, result);
}}

double {phase}_{module}_d3gdt2dp(double T, double P, double n[{number_components}]) {{
    return {module}_d3gdt2dp(T, P, n);
}}

void {phase}_{module}_d4gdndt2dp(double T, double P, double n[{number_components}], double result[{number_components}]) {{
    {module}_d4gdndt2dp(T, P, n, result);
}}

void {phase}_{module}_d5gdn2dt2dp(double T, double P, double n[{number_components}], double result[{number_components}]) {{
    {module}_d5gdn2dt2dp(T, P, n, result);
}}

void {phase}_{module}_d6gdn3dt2dp(double T, double P, double n[{number_components}], double result[{number_components}]) {{
    {module}_d6gdn3dt2dp(T, P, n, result);
}}

double {phase}_{module}_d3gdtdp2(double T, double P, double n[{number_components}]) {{
    return {module}_d3gdtdp2(T, P, n);
}}

void {phase}_{module}_d4gdndtdp2(double T, double P, double n[{number_components}], double result[{number_components}]) {{
    {module}_d4gdndtdp2(T, P, n, result);
}}

void {phase}_{module}_d5gdn2dtdp2(double T, double P, double n[{number_components}], double result[{number_components}]) {{
    {module}_d5gdn2dtdp2(T, P, n, result);
}}

void {phase}_{module}_d6gdn3dtdp2(double T, double P, double n[{number_components}], double result[{number_components}]) {{
    {module}_d6gdn3dtdp2(T, P, n, result);
}}

double {phase}_{module}_d3gdp3(double T, double P, double n[{number_components}]) {{
    return {module}_d3gdp3(T, P, n);
}}

void {phase}_{module}_d4gdndp3(double T, double P, double n[{number_components}], double result[{number_components}]) {{
    {module}_d4gdndp3(T, P, n, result);
}}

void {phase}_{module}_d5gdn2dp3(double T, double P, double n[{number_components}], double result[{number_components}]) {{
    {module}_d5gdn2dp3(T, P, n, result);
}}

void {phase}_{module}_d6gdn3dp3(double T, double P, double n[{number_components}], double result[{number_components}]) {{
    {module}_d6gdn3dp3(T, P, n, result);
}}

double {phase}_{module}_s(double T, double P, double n[{number_components}]) {{
    return {module}_s(T, P, n);
}}

double {phase}_{module}_v(double T, double P, double n[{number_components}]) {{
    return {module}_v(T, P, n);
}}

double {phase}_{module}_cv(double T, double P, double n[{number_components}]) {{
    return {module}_cv(T, P, n);
}}

double {phase}_{module}_cp(double T, double P, double n[{number_components}]) {{
    return {module}_cp(T, P, n);
}}

double {phase}_{module}_dcpdt(double T, double P, double n[{number_components}]) {{
    return {module}_dcpdt(T, P, n);
}}

double {phase}_{module}_alpha(double T, double P, double n[{number_components}]) {{
    return {module}_alpha(T, P, n);
}}

double {phase}_{module}_beta(double T, double P, double n[{number_components}]) {{
    return {module}_beta(T, P, n);
}}

double {phase}_{module}_K(double T, double P, double n[{number_components}]) {{
    return {module}_K(T, P, n);
}}

double {phase}_{module}_Kp(double T, double P, double n[{number_components}]) {{
    return {module}_Kp(T, P, n);
}}

\
"""

def _create_soln_fast_include_template_c():
    """
    C language implementation of create_soln_fast_include_template()

    Returns
    -------
    string :
        The template string.
    """

    return """\

const char *{phase}_{module}_identifier(void);
const char *{phase}_{module}_name(void);
char *{phase}_{module}_formula(double T, double P, double n[{number_components}]);
double *{phase}_{module}_conv_elm_to_moles(double *e);
int {phase}_{module}_test_moles(double *n);

const char *{phase}_{module}_endmember_name(int index);
const char *{phase}_{module}_endmember_formula(int index);
const double {phase}_{module}_endmember_mw(int index);
const double *{phase}_{module}_endmember_elements(int index);
double {phase}_{module}_endmember_mu0(int index, double t, double p);
double {phase}_{module}_endmember_dmu0dT(int index, double t, double p);
double {phase}_{module}_endmember_dmu0dP(int index, double t, double p);
double {phase}_{module}_endmember_d2mu0dT2(int index, double t, double p);
double {phase}_{module}_endmember_d2mu0dTdP(int index, double t, double p);
double {phase}_{module}_endmember_d2mu0dP2(int index, double t, double p);
double {phase}_{module}_endmember_d3mu0dT3(int index, double t, double p);
double {phase}_{module}_endmember_d3mu0dT2dP(int index, double t, double p);
double {phase}_{module}_endmember_d3mu0dTdP2(int index, double t, double p);
double {phase}_{module}_endmember_d3mu0dP3(int index, double t, double p);

double {phase}_{module}_g(double T, double P, double n[{number_components}]);
double {phase}_{module}_dgdt(double T, double P, double n[{number_components}]);
double {phase}_{module}_dgdp(double T, double P, double n[{number_components}]);
double {phase}_{module}_d2gdt2(double T, double P, double n[{number_components}]);
double {phase}_{module}_d2gdtdp(double T, double P, double n[{number_components}]);
double {phase}_{module}_d2gdp2(double T, double P, double n[{number_components}]);
double {phase}_{module}_d3gdt3(double T, double P, double n[{number_components}]);
double {phase}_{module}_d3gdt2dp(double T, double P, double n[{number_components}]);
double {phase}_{module}_d3gdtdp2(double T, double P, double n[{number_components}]);
double {phase}_{module}_d3gdp3(double T, double P, double n[{number_components}]);

void {phase}_{module}_dgdn(double T, double P, double n[{number_components}], double result[{number_components}]);
void {phase}_{module}_d2gdndt(double T, double P, double n[{number_components}], double result[{number_components}]);
void {phase}_{module}_d2gdndp(double T, double P, double n[{number_components}], double result[{number_components}]);
void {phase}_{module}_d3gdndt2(double T, double P, double n[{number_components}], double result[{number_components}]);
void {phase}_{module}_d3gdndtdp(double T, double P, double n[{number_components}], double result[{number_components}]);
void {phase}_{module}_d3gdndp2(double T, double P, double n[{number_components}], double result[{number_components}]);
void {phase}_{module}_d4gdndt3(double T, double P, double n[{number_components}], double result[{number_components}]);
void {phase}_{module}_d4gdndt2dp(double T, double P, double n[{number_components}], double result[{number_components}]);
void {phase}_{module}_d4gdndtdp2(double T, double P, double n[{number_components}], double result[{number_components}]);
void {phase}_{module}_d4gdndp3(double T, double P, double n[{number_components}], double result[{number_components}]);

void {phase}_{module}_d2gdn2(double T, double P, double n[{number_components}], double result[{number_components}]);
void {phase}_{module}_d3gdn2dt(double T, double P, double n[{number_components}], double result[{number_components}]);
void {phase}_{module}_d3gdn2dp(double T, double P, double n[{number_components}], double result[{number_components}]);
void {phase}_{module}_d4gdn2dt2(double T, double P, double n[{number_components}], double result[{number_components}]);
void {phase}_{module}_d4gdn2dtdp(double T, double P, double n[{number_components}], double result[{number_components}]);
void {phase}_{module}_d4gdn2dp2(double T, double P, double n[{number_components}], double result[{number_components}]);
void {phase}_{module}_d5gdn2dt3(double T, double P, double n[{number_components}], double result[{number_components}]);
void {phase}_{module}_d5gdn2dt2dp(double T, double P, double n[{number_components}], double result[{number_components}]);
void {phase}_{module}_d5gdn2dtdp2(double T, double P, double n[{number_components}], double result[{number_components}]);
void {phase}_{module}_d5gdn2dp3(double T, double P, double n[{number_components}], double result[{number_components}]);

void {phase}_{module}_d3gdn3(double T, double P, double n[{number_components}], double result[{number_components}]);
void {phase}_{module}_d4gdn3dt(double T, double P, double n[{number_components}], double result[{number_components}]);
void {phase}_{module}_d4gdn3dp(double T, double P, double n[{number_components}], double result[{number_components}]);
void {phase}_{module}_d5gdn3dt2(double T, double P, double n[{number_components}], double result[{number_components}]);
void {phase}_{module}_d5gdn3dtdp(double T, double P, double n[{number_components}], double result[{number_components}]);
void {phase}_{module}_d5gdn3dp2(double T, double P, double n[{number_components}], double result[{number_components}]);
void {phase}_{module}_d6gdn3dt3(double T, double P, double n[{number_components}], double result[{number_components}]);
void {phase}_{module}_d6gdn3dt2dp(double T, double P, double n[{number_components}], double result[{number_components}]);
void {phase}_{module}_d6gdn3dtdp2(double T, double P, double n[{number_components}], double result[{number_components}]);
void {phase}_{module}_d6gdn3dp3(double T, double P, double n[{number_components}], double result[{number_components}]);

double {phase}_{module}_s(double T, double P, double n[{number_components}]);
double {phase}_{module}_v(double T, double P, double n[{number_components}]);
double {phase}_{module}_cv(double T, double P, double n[{number_components}]);
double {phase}_{module}_cp(double T, double P, double n[{number_components}]);
double {phase}_{module}_dcpdt(double T, double P, double n[{number_components}]);
double {phase}_{module}_alpha(double T, double P, double n[{number_components}]);
double {phase}_{module}_beta(double T, double P, double n[{number_components}]);
double {phase}_{module}_K(double T, double P, double n[{number_components}]);
double {phase}_{module}_Kp(double T, double P, double n[{number_components}]);

\
"""

def _create_soln_fast_pyx_template_c():
    """
    C language implementation of create_soln_fast_pyx_template()

    Returns
    -------
    string :
        The template string.
    """

    return """\
# Cython numpy wrapper code for arrays is taken from:
# http://gael-varoquaux.info/programming/cython-example-of-exposing-c-computed-arrays-in-python-without-data-copies.html
# Author: Gael Varoquaux, BSD license

# Declare the prototype of the C functions
cdef extern from "{phase}_{module}_calc.h":
    const char *{phase}_{module}_identifier();
    const char *{phase}_{module}_name();
    char *{phase}_{module}_formula(double T, double P, double n[{number_components}]);
    double *{phase}_{module}_conv_elm_to_moles(double *e);
    int {phase}_{module}_test_moles(double *n);
    const char *{phase}_{module}_endmember_name(int index);
    const char *{phase}_{module}_endmember_formula(int index);
    const double {phase}_{module}_endmember_mw(int index);
    const double *{phase}_{module}_endmember_elements(int index);
    double {phase}_{module}_endmember_mu0(int index, double t, double p);
    double {phase}_{module}_endmember_dmu0dT(int index, double t, double p);
    double {phase}_{module}_endmember_dmu0dP(int index, double t, double p);
    double {phase}_{module}_endmember_d2mu0dT2(int index, double t, double p);
    double {phase}_{module}_endmember_d2mu0dTdP(int index, double t, double p);
    double {phase}_{module}_endmember_d2mu0dP2(int index, double t, double p);
    double {phase}_{module}_endmember_d3mu0dT3(int index, double t, double p);
    double {phase}_{module}_endmember_d3mu0dT2dP(int index, double t, double p);
    double {phase}_{module}_endmember_d3mu0dTdP2(int index, double t, double p);
    double {phase}_{module}_endmember_d3mu0dP3(int index, double t, double p);
    double {phase}_{module}_g(double t, double p, double n[{number_components}])
    double {phase}_{module}_dgdt(double t, double p, double n[{number_components}])
    double {phase}_{module}_dgdp(double t, double p, double n[{number_components}])
    double {phase}_{module}_d2gdt2(double t, double p, double n[{number_components}])
    double {phase}_{module}_d2gdtdp(double t, double p, double n[{number_components}])
    double {phase}_{module}_d2gdp2(double t, double p, double n[{number_components}])
    double {phase}_{module}_d3gdt3(double t, double p, double n[{number_components}])
    double {phase}_{module}_d3gdt2dp(double t, double p, double n[{number_components}])
    double {phase}_{module}_d3gdtdp2(double t, double p, double n[{number_components}])
    double {phase}_{module}_d3gdp3(double t, double p, double n[{number_components}])
    double {phase}_{module}_s(double t, double p, double n[{number_components}])
    double {phase}_{module}_v(double t, double p, double n[{number_components}])
    double {phase}_{module}_cv(double t, double p, double n[{number_components}])
    double {phase}_{module}_cp(double t, double p, double n[{number_components}])
    double {phase}_{module}_dcpdt(double t, double p, double n[{number_components}])
    double {phase}_{module}_alpha(double t, double p, double n[{number_components}])
    double {phase}_{module}_beta(double t, double p, double n[{number_components}])
    double {phase}_{module}_K(double t, double p, double n[{number_components}])
    double {phase}_{module}_Kp(double t, double p, double n[{number_components}])

    void {phase}_{module}_dgdn(double T, double P, double n[{number_components}], double result[{number_components}]);
    void {phase}_{module}_d2gdndt(double T, double P, double n[{number_components}], double result[{number_components}]);
    void {phase}_{module}_d2gdndp(double T, double P, double n[{number_components}], double result[{number_components}]);
    void {phase}_{module}_d3gdndt2(double T, double P, double n[{number_components}], double result[{number_components}]);
    void {phase}_{module}_d3gdndtdp(double T, double P, double n[{number_components}], double result[{number_components}]);
    void {phase}_{module}_d3gdndp2(double T, double P, double n[{number_components}], double result[{number_components}]);
    void {phase}_{module}_d4gdndt3(double T, double P, double n[{number_components}], double result[{number_components}]);
    void {phase}_{module}_d4gdndt2dp(double T, double P, double n[{number_components}], double result[{number_components}]);
    void {phase}_{module}_d4gdndtdp2(double T, double P, double n[{number_components}], double result[{number_components}]);
    void {phase}_{module}_d4gdndp3(double T, double P, double n[{number_components}], double result[{number_components}]);

    void {phase}_{module}_d2gdn2(double T, double P, double n[{number_components}], double result[{number_components}]);
    void {phase}_{module}_d3gdn2dt(double T, double P, double n[{number_components}], double result[{number_components}]);
    void {phase}_{module}_d3gdn2dp(double T, double P, double n[{number_components}], double result[{number_components}]);
    void {phase}_{module}_d4gdn2dt2(double T, double P, double n[{number_components}], double result[{number_components}]);
    void {phase}_{module}_d4gdn2dtdp(double T, double P, double n[{number_components}], double result[{number_components}]);
    void {phase}_{module}_d4gdn2dp2(double T, double P, double n[{number_components}], double result[{number_components}]);
    void {phase}_{module}_d5gdn2dt3(double T, double P, double n[{number_components}], double result[{number_components}]);
    void {phase}_{module}_d5gdn2dt2dp(double T, double P, double n[{number_components}], double result[{number_components}]);
    void {phase}_{module}_d5gdn2dtdp2(double T, double P, double n[{number_components}], double result[{number_components}]);
    void {phase}_{module}_d5gdn2dp3(double T, double P, double n[{number_components}], double result[{number_components}]);

    void {phase}_{module}_d3gdn3(double T, double P, double n[{number_components}], double result[{number_components}]);
    void {phase}_{module}_d4gdn3dt(double T, double P, double n[{number_components}], double result[{number_components}]);
    void {phase}_{module}_d4gdn3dp(double T, double P, double n[{number_components}], double result[{number_components}]);
    void {phase}_{module}_d5gdn3dt2(double T, double P, double n[{number_components}], double result[{number_components}]);
    void {phase}_{module}_d5gdn3dtdp(double T, double P, double n[{number_components}], double result[{number_components}]);
    void {phase}_{module}_d5gdn3dp2(double T, double P, double n[{number_components}], double result[{number_components}]);
    void {phase}_{module}_d6gdn3dt3(double T, double P, double n[{number_components}], double result[{number_components}]);
    void {phase}_{module}_d6gdn3dt2dp(double T, double P, double n[{number_components}], double result[{number_components}]);
    void {phase}_{module}_d6gdn3dtdp2(double T, double P, double n[{number_components}], double result[{number_components}]);
    void {phase}_{module}_d6gdn3dp3(double T, double P, double n[{number_components}], double result[{number_components}]);

from libc.stdlib cimport malloc, free
from cpython cimport PyObject, Py_INCREF
import ctypes

# Import the Python-level symbols of numpy
import numpy as np

# Import the C-level symbols of numpy
cimport numpy as np

# Numpy must be initialized. When using numpy from C or Cython you must
# _always_ do that, or you will have segfaults
np.import_array()

# here is the "wrapper" signature
def {prefix}_{phase}_{module}_identifier():
    result = <bytes> {phase}_{module}_identifier()
    return result.decode('UTF-8')
def {prefix}_{phase}_{module}_name():
    result = <bytes> {phase}_{module}_name()
    return result.decode('UTF-8')
def {prefix}_{phase}_{module}_formula(double t, double p, np_array):
    cdef double *m = <double *>malloc(len(np_array)*sizeof(double))
    for i in range(np_array.size):
        m[i] = np_array[i]
    r = {phase}_{module}_formula(<double> t, <double> p, <double *> m)
    result = <bytes> r
    free (m)
    free (r)
    return result.decode('UTF-8')
def {prefix}_{phase}_{module}_conv_elm_to_moles(np_array):
    cdef double *e = <double *>malloc(len(np_array)*sizeof(double))
    for i in range(np_array.size):
        e[i] = np_array[i]
    r = {phase}_{module}_conv_elm_to_moles(<double *> e)
    result = []
    for i in range({number_components}):
        result.append(r[i])
    free(r)
    return np.array(result)
def {prefix}_{phase}_{module}_test_moles(np_array):
    cdef double *n = <double *>malloc(len(np_array)*sizeof(double))
    for i in range(np_array.size):
        n[i] = np_array[i]
    result = {phase}_{module}_test_moles(<double *> n)
    return False if result == 0 else True

def {prefix}_{phase}_{module}_endmember_number():
    return {number_components}
def {prefix}_{phase}_{module}_endmember_name(int index):
    assert index in range(0,{number_components}), "index out of range"
    result = {phase}_{module}_endmember_name(index);
    return result.decode('UTF-8')
def {prefix}_{phase}_{module}_endmember_formula(int index):
    assert index in range(0,{number_components}), "index out of range"
    result = {phase}_{module}_endmember_formula(index);
    return result.decode('UTF-8')
def {prefix}_{phase}_{module}_endmember_mw(int index):
    assert index in range(0,{number_components}), "index out of range"
    return {phase}_{module}_endmember_mw(index);
def {prefix}_{phase}_{module}_endmember_elements(int index):
    assert index in range(0,{number_components}), "index out of range"
    r = {phase}_{module}_endmember_elements(index);
    result = []
    for i in range(0,106):
        result.append(r[i])
    return np.array(result)

def {prefix}_{phase}_{module}_species_number():
    return {number_species}
def {prefix}_{phase}_{module}_species_name(int index):
    assert index in range(0,{number_species}), "index out of range"
    result = {phase}_{module}_endmember_name(index);
    return result.decode('UTF-8')
def {prefix}_{phase}_{module}_species_formula(int index):
    assert index in range(0,{number_species}), "index out of range"
    result = {phase}_{module}_endmember_formula(index);
    return result.decode('UTF-8')
def {prefix}_{phase}_{module}_species_mw(int index):
    assert index in range(0,{number_species}), "index out of range"
    return {phase}_{module}_endmember_mw(index);
def {prefix}_{phase}_{module}_species_elements(int index):
    assert index in range(0,{number_species}), "index out of range"
    r = {phase}_{module}_endmember_elements(index);
    result = []
    for i in range(0,106):
        result.append(r[i])
    return np.array(result)

def {prefix}_{phase}_{module}_endmember_mu0(int index, double t, double p):
    assert index in range(0,{number_components}), "index out of range"
    return {phase}_{module}_endmember_mu0(index, <double> t, <double> p);
def {prefix}_{phase}_{module}_endmember_dmu0dT(int index, double t, double p):
    assert index in range(0,{number_components}), "index out of range"
    return {phase}_{module}_endmember_dmu0dT(index, <double> t, <double> p);
def {prefix}_{phase}_{module}_endmember_dmu0dP(int index, double t, double p):
    assert index in range(0,{number_components}), "index out of range"
    return {phase}_{module}_endmember_dmu0dP(index, <double> t, <double> p);
def {prefix}_{phase}_{module}_endmember_d2mu0dT2(int index, double t, double p):
    assert index in range(0,{number_components}), "index out of range"
    return {phase}_{module}_endmember_d2mu0dT2(index, <double> t, <double> p);
def {prefix}_{phase}_{module}_endmember_d2mu0dTdP(int index, double t, double p):
    assert index in range(0,{number_components}), "index out of range"
    return {phase}_{module}_endmember_d2mu0dTdP(index, <double> t, <double> p);
def {prefix}_{phase}_{module}_endmember_d2mu0dP2(int index, double t, double p):
    assert index in range(0,{number_components}), "index out of range"
    return {phase}_{module}_endmember_d2mu0dP2(index, <double> t, <double> p);
def {prefix}_{phase}_{module}_endmember_d3mu0dT3(int index, double t, double p):
    assert index in range(0,{number_components}), "index out of range"
    return {phase}_{module}_endmember_d3mu0dT3(index, <double> t, <double> p);
def {prefix}_{phase}_{module}_endmember_d3mu0dT2dP(int index, double t, double p):
    assert index in range(0,{number_components}), "index out of range"
    return {phase}_{module}_endmember_d3mu0dT2dP(index, <double> t, <double> p);
def {prefix}_{phase}_{module}_endmember_d3mu0dTdP2(int index, double t, double p):
    assert index in range(0,{number_components}), "index out of range"
    return {phase}_{module}_endmember_d3mu0dTdP2(index, <double> t, <double> p);
def {prefix}_{phase}_{module}_endmember_d3mu0dP3(int index, double t, double p):
    assert index in range(0,{number_components}), "index out of range"
    return {phase}_{module}_endmember_d3mu0dP3(index, <double> t, <double> p);

def {prefix}_{phase}_{module}_g(double t, double p, np_array):
    cdef double *m = <double *>malloc(len(np_array)*sizeof(double))
    for i in range(np_array.size):
        m[i] = np_array[i]
    result = {phase}_{module}_g(<double> t, <double> p, <double *> m)
    free (m)
    return result
def {prefix}_{phase}_{module}_dgdt(double t, double p, np_array):
    cdef double *m = <double *>malloc(len(np_array)*sizeof(double))
    for i in range(np_array.size):
        m[i] = np_array[i]
    result = {phase}_{module}_dgdt(<double> t, <double> p, <double *> m)
    free (m)
    return result
def {prefix}_{phase}_{module}_dgdp(double t, double p, np_array):
    cdef double *m = <double *>malloc(len(np_array)*sizeof(double))
    for i in range(np_array.size):
        m[i] = np_array[i]
    result = {phase}_{module}_dgdp(<double> t, <double> p, <double *> m)
    free (m)
    return result
def {prefix}_{phase}_{module}_d2gdt2(double t, double p, np_array):
    cdef double *m = <double *>malloc(len(np_array)*sizeof(double))
    for i in range(np_array.size):
        m[i] = np_array[i]
    result = {phase}_{module}_d2gdt2(<double> t, <double> p, <double *> m)
    free (m)
    return result
def {prefix}_{phase}_{module}_d2gdtdp(double t, double p, np_array):
    cdef double *m = <double *>malloc(len(np_array)*sizeof(double))
    for i in range(np_array.size):
        m[i] = np_array[i]
    result = {phase}_{module}_d2gdtdp(<double> t, <double> p, <double *> m)
    free (m)
    return result
def {prefix}_{phase}_{module}_d2gdp2(double t, double p, np_array):
    cdef double *m = <double *>malloc(len(np_array)*sizeof(double))
    for i in range(np_array.size):
        m[i] = np_array[i]
    result = {phase}_{module}_d2gdp2(<double> t, <double> p, <double *> m)
    free (m)
    return result
def {prefix}_{phase}_{module}_d3gdt3(double t, double p, np_array):
    cdef double *m = <double *>malloc(len(np_array)*sizeof(double))
    for i in range(np_array.size):
        m[i] = np_array[i]
    result = {phase}_{module}_d3gdt3(<double> t, <double> p, <double *> m)
    free (m)
    return result
def {prefix}_{phase}_{module}_d3gdt2dp(double t, double p, np_array):
    cdef double *m = <double *>malloc(len(np_array)*sizeof(double))
    for i in range(np_array.size):
        m[i] = np_array[i]
    result = {phase}_{module}_d3gdt2dp(<double> t, <double> p, <double *> m)
    free (m)
    return result
def {prefix}_{phase}_{module}_d3gdtdp2(double t, double p, np_array):
    cdef double *m = <double *>malloc(len(np_array)*sizeof(double))
    for i in range(np_array.size):
        m[i] = np_array[i]
    result = {phase}_{module}_d3gdtdp2(<double> t, <double> p, <double *> m)
    free (m)
    return result
def {prefix}_{phase}_{module}_d3gdp3(double t, double p, np_array):
    cdef double *m = <double *>malloc(len(np_array)*sizeof(double))
    for i in range(np_array.size):
        m[i] = np_array[i]
    result = {phase}_{module}_d3gdp3(<double> t, <double> p, <double *> m)
    free (m)
    return result
def {prefix}_{phase}_{module}_s(double t, double p, np_array):
    cdef double *m = <double *>malloc(len(np_array)*sizeof(double))
    for i in range(np_array.size):
        m[i] = np_array[i]
    result = {phase}_{module}_s(<double> t, <double> p, <double *> m)
    free (m)
    return result
def {prefix}_{phase}_{module}_v(double t, double p, np_array):
    cdef double *m = <double *>malloc(len(np_array)*sizeof(double))
    for i in range(np_array.size):
        m[i] = np_array[i]
    result = {phase}_{module}_v(<double> t, <double> p, <double *> m)
    free (m)
    return result
def {prefix}_{phase}_{module}_cv(double t, double p, np_array):
    cdef double *m = <double *>malloc(len(np_array)*sizeof(double))
    for i in range(np_array.size):
        m[i] = np_array[i]
    result = {phase}_{module}_cv(<double> t, <double> p, <double *> m)
    free (m)
    return result
def {prefix}_{phase}_{module}_cp(double t, double p, np_array):
    cdef double *m = <double *>malloc(len(np_array)*sizeof(double))
    for i in range(np_array.size):
        m[i] = np_array[i]
    result = {phase}_{module}_cp(<double> t, <double> p, <double *> m)
    free (m)
    return result
def {prefix}_{phase}_{module}_dcpdt(double t, double p, np_array):
    cdef double *m = <double *>malloc(len(np_array)*sizeof(double))
    for i in range(np_array.size):
        m[i] = np_array[i]
    result = {phase}_{module}_dcpdt(<double> t, <double> p, <double *> m)
    free (m)
    return result
def {prefix}_{phase}_{module}_alpha(double t, double p, np_array):
    cdef double *m = <double *>malloc(len(np_array)*sizeof(double))
    for i in range(np_array.size):
        m[i] = np_array[i]
    result = {phase}_{module}_alpha(<double> t, <double> p, <double *> m)
    free (m)
    return result
def {prefix}_{phase}_{module}_beta(double t, double p, np_array):
    cdef double *m = <double *>malloc(len(np_array)*sizeof(double))
    for i in range(np_array.size):
        m[i] = np_array[i]
    result = {phase}_{module}_beta(<double> t, <double> p, <double *> m)
    free (m)
    return result
def {prefix}_{phase}_{module}_K(double t, double p, np_array):
    cdef double *m = <double *>malloc(len(np_array)*sizeof(double))
    for i in range(np_array.size):
        m[i] = np_array[i]
    result = {phase}_{module}_K(<double> t, <double> p, <double *> m)
    free (m)
    return result
def {prefix}_{phase}_{module}_Kp(double t, double p, np_array):
    cdef double *m = <double *>malloc(len(np_array)*sizeof(double))
    for i in range(np_array.size):
        m[i] = np_array[i]
    result = {phase}_{module}_Kp(<double> t, <double> p, <double *> m)
    free (m)
    return result

def {prefix}_{phase}_{module}_dgdn(double t, double p, np_array):
    cdef double *m = <double *>malloc(len(np_array)*sizeof(double))
    for i in range(np_array.size):
        m[i] = np_array[i]
    cdef double *r = <double *>malloc({number_components}*sizeof(double))
    {phase}_{module}_dgdn(<double> t, <double> p, <double *> m, <double *> r)
    result = []
    for i in range({number_components}):
        result.append(r[i])
    free (m)
    free (r)
    return np.array(result)
def {prefix}_{phase}_{module}_d2gdndt(double t, double p, np_array):
    cdef double *m = <double *>malloc(len(np_array)*sizeof(double))
    for i in range(np_array.size):
        m[i] = np_array[i]
    cdef double *r = <double *>malloc({number_components}*sizeof(double))
    {phase}_{module}_d2gdndt(<double> t, <double> p, <double *> m, <double *> r)
    result = []
    for i in range({number_components}):
        result.append(r[i])
    free (m)
    free (r)
    return np.array(result)
def {prefix}_{phase}_{module}_d2gdndp(double t, double p, np_array):
    cdef double *m = <double *>malloc(len(np_array)*sizeof(double))
    for i in range(np_array.size):
        m[i] = np_array[i]
    cdef double *r = <double *>malloc({number_components}*sizeof(double))
    {phase}_{module}_d2gdndp(<double> t, <double> p, <double *> m, <double *> r)
    result = []
    for i in range({number_components}):
        result.append(r[i])
    free (m)
    free (r)
    return np.array(result)
def {prefix}_{phase}_{module}_d3gdndt2(double t, double p, np_array):
    cdef double *m = <double *>malloc(len(np_array)*sizeof(double))
    for i in range(np_array.size):
        m[i] = np_array[i]
    cdef double *r = <double *>malloc({number_components}*sizeof(double))
    {phase}_{module}_d3gdndt2(<double> t, <double> p, <double *> m, <double *> r)
    result = []
    for i in range({number_components}):
        result.append(r[i])
    free (m)
    free (r)
    return np.array(result)
def {prefix}_{phase}_{module}_d3gdndtdp(double t, double p, np_array):
    cdef double *m = <double *>malloc(len(np_array)*sizeof(double))
    for i in range(np_array.size):
        m[i] = np_array[i]
    cdef double *r = <double *>malloc({number_components}*sizeof(double))
    {phase}_{module}_d3gdndtdp(<double> t, <double> p, <double *> m, <double *> r)
    result = []
    for i in range({number_components}):
        result.append(r[i])
    free (m)
    free (r)
    return np.array(result)
def {prefix}_{phase}_{module}_d3gdndp2(double t, double p, np_array):
    cdef double *m = <double *>malloc(len(np_array)*sizeof(double))
    for i in range(np_array.size):
        m[i] = np_array[i]
    cdef double *r = <double *>malloc({number_components}*sizeof(double))
    {phase}_{module}_d3gdndp2(<double> t, <double> p, <double *> m, <double *> r)
    result = []
    for i in range({number_components}):
        result.append(r[i])
    free (m)
    free (r)
    return np.array(result)
def {prefix}_{phase}_{module}_d4gdndt3(double t, double p, np_array):
    cdef double *m = <double *>malloc(len(np_array)*sizeof(double))
    for i in range(np_array.size):
        m[i] = np_array[i]
    cdef double *r = <double *>malloc({number_components}*sizeof(double))
    {phase}_{module}_d4gdndt3(<double> t, <double> p, <double *> m, <double *> r)
    result = []
    for i in range({number_components}):
        result.append(r[i])
    free (m)
    free (r)
    return np.array(result)
def {prefix}_{phase}_{module}_d4gdndt2dp(double t, double p, np_array):
    cdef double *m = <double *>malloc(len(np_array)*sizeof(double))
    for i in range(np_array.size):
        m[i] = np_array[i]
    cdef double *r = <double *>malloc({number_components}*sizeof(double))
    {phase}_{module}_d4gdndt2dp(<double> t, <double> p, <double *> m, <double *> r)
    result = []
    for i in range({number_components}):
        result.append(r[i])
    free (m)
    free (r)
    return np.array(result)
def {prefix}_{phase}_{module}_d4gdndtdp2(double t, double p, np_array):
    cdef double *m = <double *>malloc(len(np_array)*sizeof(double))
    for i in range(np_array.size):
        m[i] = np_array[i]
    cdef double *r = <double *>malloc({number_components}*sizeof(double))
    {phase}_{module}_d4gdndtdp2(<double> t, <double> p, <double *> m, <double *> r)
    result = []
    for i in range({number_components}):
        result.append(r[i])
    free (m)
    free (r)
    return np.array(result)
def {prefix}_{phase}_{module}_d4gdndp3(double t, double p, np_array):
    cdef double *m = <double *>malloc(len(np_array)*sizeof(double))
    for i in range(np_array.size):
        m[i] = np_array[i]
    cdef double *r = <double *>malloc({number_components}*sizeof(double))
    {phase}_{module}_d4gdndp3(<double> t, <double> p, <double *> m, <double *> r)
    result = []
    for i in range({number_components}):
        result.append(r[i])
    free (m)
    free (r)
    return np.array(result)

def {prefix}_{phase}_{module}_d2gdn2(double t, double p, np_array):
    cdef double *m = <double *>malloc(len(np_array)*sizeof(double))
    for i in range(np_array.size):
        m[i] = np_array[i]
    ndim = int({number_components}*({number_components}-1)/2 + {number_components})
    cdef double *r = <double *>malloc(ndim*sizeof(double))
    {phase}_{module}_d2gdn2(<double> t, <double> p, <double *> m, <double *> r)
    result = []
    for i in range(ndim):
        result.append(r[i])
    free (m)
    free (r)
    return np.array(result)
def {prefix}_{phase}_{module}_d3gdn2dt(double t, double p, np_array):
    cdef double *m = <double *>malloc(len(np_array)*sizeof(double))
    for i in range(np_array.size):
        m[i] = np_array[i]
    ndim = int({number_components}*({number_components}-1)/2 + {number_components})
    cdef double *r = <double *>malloc(ndim*sizeof(double))
    {phase}_{module}_d3gdn2dt(<double> t, <double> p, <double *> m, <double *> r)
    result = []
    for i in range(ndim):
        result.append(r[i])
    free (m)
    free (r)
    return np.array(result)
def {prefix}_{phase}_{module}_d3gdn2dp(double t, double p, np_array):
    cdef double *m = <double *>malloc(len(np_array)*sizeof(double))
    for i in range(np_array.size):
        m[i] = np_array[i]
    ndim = int({number_components}*({number_components}-1)/2 + {number_components})
    cdef double *r = <double *>malloc(ndim*sizeof(double))
    {phase}_{module}_d3gdn2dp(<double> t, <double> p, <double *> m, <double *> r)
    result = []
    for i in range(ndim):
        result.append(r[i])
    free (m)
    free (r)
    return np.array(result)
def {prefix}_{phase}_{module}_d4gdn2dt2(double t, double p, np_array):
    cdef double *m = <double *>malloc(len(np_array)*sizeof(double))
    for i in range(np_array.size):
        m[i] = np_array[i]
    ndim = int({number_components}*({number_components}-1)/2 + {number_components})
    cdef double *r = <double *>malloc(ndim*sizeof(double))
    {phase}_{module}_d4gdn2dt2(<double> t, <double> p, <double *> m, <double *> r)
    result = []
    for i in range(ndim):
        result.append(r[i])
    free (m)
    free (r)
    return np.array(result)
def {prefix}_{phase}_{module}_d4gdn2dtdp(double t, double p, np_array):
    cdef double *m = <double *>malloc(len(np_array)*sizeof(double))
    for i in range(np_array.size):
        m[i] = np_array[i]
    ndim = int({number_components}*({number_components}-1)/2 + {number_components})
    cdef double *r = <double *>malloc(ndim*sizeof(double))
    {phase}_{module}_d4gdn2dtdp(<double> t, <double> p, <double *> m, <double *> r)
    result = []
    for i in range(ndim):
        result.append(r[i])
    free (m)
    free (r)
    return np.array(result)
def {prefix}_{phase}_{module}_d4gdn2dp2(double t, double p, np_array):
    cdef double *m = <double *>malloc(len(np_array)*sizeof(double))
    for i in range(np_array.size):
        m[i] = np_array[i]
    ndim = int({number_components}*({number_components}-1)/2 + {number_components})
    cdef double *r = <double *>malloc(ndim*sizeof(double))
    {phase}_{module}_d4gdn2dp2(<double> t, <double> p, <double *> m, <double *> r)
    result = []
    for i in range(ndim):
        result.append(r[i])
    free (m)
    free (r)
    return np.array(result)
def {prefix}_{phase}_{module}_d5gdn2dt3(double t, double p, np_array):
    cdef double *m = <double *>malloc(len(np_array)*sizeof(double))
    for i in range(np_array.size):
        m[i] = np_array[i]
    ndim = int({number_components}*({number_components}-1)/2 + {number_components})
    cdef double *r = <double *>malloc(ndim*sizeof(double))
    {phase}_{module}_d5gdn2dt3(<double> t, <double> p, <double *> m, <double *> r)
    result = []
    for i in range(ndim):
        result.append(r[i])
    free (m)
    free (r)
    return np.array(result)
def {prefix}_{phase}_{module}_d5gdn2dt2dp(double t, double p, np_array):
    cdef double *m = <double *>malloc(len(np_array)*sizeof(double))
    for i in range(np_array.size):
        m[i] = np_array[i]
    ndim = int({number_components}*({number_components}-1)/2 + {number_components})
    cdef double *r = <double *>malloc(ndim*sizeof(double))
    {phase}_{module}_d5gdn2dt2dp(<double> t, <double> p, <double *> m, <double *> r)
    result = []
    for i in range(ndim):
        result.append(r[i])
    free (m)
    free (r)
    return np.array(result)
def {prefix}_{phase}_{module}_d5gdn2dtdp2(double t, double p, np_array):
    cdef double *m = <double *>malloc(len(np_array)*sizeof(double))
    for i in range(np_array.size):
        m[i] = np_array[i]
    ndim = int({number_components}*({number_components}-1)/2 + {number_components})
    cdef double *r = <double *>malloc(ndim*sizeof(double))
    {phase}_{module}_d5gdn2dtdp2(<double> t, <double> p, <double *> m, <double *> r)
    result = []
    for i in range(ndim):
        result.append(r[i])
    free (m)
    free (r)
    return np.array(result)
def {prefix}_{phase}_{module}_d5gdn2dp3(double t, double p, np_array):
    cdef double *m = <double *>malloc(len(np_array)*sizeof(double))
    for i in range(np_array.size):
        m[i] = np_array[i]
    ndim = int({number_components}*({number_components}-1)/2 + {number_components})
    cdef double *r = <double *>malloc(ndim*sizeof(double))
    {phase}_{module}_d5gdn2dp3(<double> t, <double> p, <double *> m, <double *> r)
    result = []
    for i in range(ndim):
        result.append(r[i])
    free (m)
    free (r)
    return np.array(result)

def {prefix}_{phase}_{module}_d3gdn3(double t, double p, np_array):
    cdef double *m = <double *>malloc(len(np_array)*sizeof(double))
    for i in range(np_array.size):
        m[i] = np_array[i]
    nc = {number_components}
    ndim = int(nc*(nc+1)*(nc+2)/6)
    cdef double *r = <double *>malloc(ndim*sizeof(double))
    {phase}_{module}_d3gdn3(<double> t, <double> p, <double *> m, <double *> r)
    result = []
    for i in range(ndim):
        result.append(r[i])
    free (m)
    free (r)
    return np.array(result)
def {prefix}_{phase}_{module}_d4gdn3dt(double t, double p, np_array):
    cdef double *m = <double *>malloc(len(np_array)*sizeof(double))
    for i in range(np_array.size):
        m[i] = np_array[i]
    nc = {number_components}
    ndim = int(nc*(nc+1)*(nc+2)/6)
    cdef double *r = <double *>malloc(ndim*sizeof(double))
    {phase}_{module}_d4gdn3dt(<double> t, <double> p, <double *> m, <double *> r)
    result = []
    for i in range(ndim):
        result.append(r[i])
    free (m)
    free (r)
    return np.array(result)
def {prefix}_{phase}_{module}_d4gdn3dp(double t, double p, np_array):
    cdef double *m = <double *>malloc(len(np_array)*sizeof(double))
    for i in range(np_array.size):
        m[i] = np_array[i]
    nc = {number_components}
    ndim = int(nc*(nc+1)*(nc+2)/6)
    cdef double *r = <double *>malloc(ndim*sizeof(double))
    {phase}_{module}_d4gdn3dp(<double> t, <double> p, <double *> m, <double *> r)
    result = []
    for i in range(ndim):
        result.append(r[i])
    free (m)
    free (r)
    return np.array(result)
def {prefix}_{phase}_{module}_d5gdn3dt2(double t, double p, np_array):
    cdef double *m = <double *>malloc(len(np_array)*sizeof(double))
    for i in range(np_array.size):
        m[i] = np_array[i]
    nc = {number_components}
    ndim = int(nc*(nc+1)*(nc+2)/6)
    cdef double *r = <double *>malloc(ndim*sizeof(double))
    {phase}_{module}_d5gdn3dt2(<double> t, <double> p, <double *> m, <double *> r)
    result = []
    for i in range(ndim):
        result.append(r[i])
    free (m)
    free (r)
    return np.array(result)
def {prefix}_{phase}_{module}_d5gdn3dtdp(double t, double p, np_array):
    cdef double *m = <double *>malloc(len(np_array)*sizeof(double))
    for i in range(np_array.size):
        m[i] = np_array[i]
    nc = {number_components}
    ndim = int(nc*(nc+1)*(nc+2)/6)
    cdef double *r = <double *>malloc(ndim*sizeof(double))
    {phase}_{module}_d5gdn3dtdp(<double> t, <double> p, <double *> m, <double *> r)
    result = []
    for i in range(ndim):
        result.append(r[i])
    free (m)
    free (r)
    return np.array(result)
def {prefix}_{phase}_{module}_d5gdn3dp2(double t, double p, np_array):
    cdef double *m = <double *>malloc(len(np_array)*sizeof(double))
    for i in range(np_array.size):
        m[i] = np_array[i]
    nc = {number_components}
    ndim = int(nc*(nc+1)*(nc+2)/6)
    cdef double *r = <double *>malloc(ndim*sizeof(double))
    {phase}_{module}_d5gdn3dp2(<double> t, <double> p, <double *> m, <double *> r)
    result = []
    for i in range(ndim):
        result.append(r[i])
    free (m)
    free (r)
    return np.array(result)
def {prefix}_{phase}_{module}_d6gdn3dt3(double t, double p, np_array):
    cdef double *m = <double *>malloc(len(np_array)*sizeof(double))
    for i in range(np_array.size):
        m[i] = np_array[i]
    nc = {number_components}
    ndim = int(nc*(nc+1)*(nc+2)/6)
    cdef double *r = <double *>malloc(ndim*sizeof(double))
    {phase}_{module}_d6gdn3dt3(<double> t, <double> p, <double *> m, <double *> r)
    result = []
    for i in range(ndim):
        result.append(r[i])
    free (m)
    free (r)
    return np.array(result)
def {prefix}_{phase}_{module}_d6gdn3dt2dp(double t, double p, np_array):
    cdef double *m = <double *>malloc(len(np_array)*sizeof(double))
    for i in range(np_array.size):
        m[i] = np_array[i]
    nc = {number_components}
    ndim = int(nc*(nc+1)*(nc+2)/6)
    cdef double *r = <double *>malloc(ndim*sizeof(double))
    {phase}_{module}_d6gdn3dt2dp(<double> t, <double> p, <double *> m, <double *> r)
    result = []
    for i in range(ndim):
        result.append(r[i])
    free (m)
    free (r)
    return np.array(result)
def {prefix}_{phase}_{module}_d6gdn3dtdp2(double t, double p, np_array):
    cdef double *m = <double *>malloc(len(np_array)*sizeof(double))
    for i in range(np_array.size):
        m[i] = np_array[i]
    nc = {number_components}
    ndim = int(nc*(nc+1)*(nc+2)/6)
    cdef double *r = <double *>malloc(ndim*sizeof(double))
    {phase}_{module}_d6gdn3dtdp2(<double> t, <double> p, <double *> m, <double *> r)
    result = []
    for i in range(ndim):
        result.append(r[i])
    free (m)
    free (r)
    return np.array(result)
def {prefix}_{phase}_{module}_d6gdn3dp3(double t, double p, np_array):
    cdef double *m = <double *>malloc(len(np_array)*sizeof(double))
    for i in range(np_array.size):
        m[i] = np_array[i]
    nc = {number_components}
    ndim = int(nc*(nc+1)*(nc+2)/6)
    cdef double *r = <double *>malloc(ndim*sizeof(double))
    {phase}_{module}_d6gdn3dp3(<double> t, <double> p, <double *> m, <double *> r)
    result = []
    for i in range(ndim):
        result.append(r[i])
    free (m)
    free (r)
    return np.array(result)
\
"""

def _create_soln_pyxbld_template_c():
    """
    C language implementation of create_pyxbld_template()

    Returns
    -------
    string :
        The template string.

    Notes
    -----
    setuptools.extension.Extension
    self,
    name,
    sources,
    include_dirs=None,
    define_macros=None,
    undef_macros=None,
    library_dirs=None,
    libraries=None,
    runtime_library_dirs=None,
    extra_objects=None,
    extra_compile_args=None,
    extra_link_args=None,
    export_symbols=None,
    swig_opts=None,
    depends=None,
    language=None,
    optional=None,
    **kw

    -O0, -O1, -O2, -O3, -Ofast, -Os, -Oz, -Og, -O, -O4
    Specify which optimization level to use:

        -O0 Means “no optimization”: this level compiles the fastest and
        generates the most debuggable code.
        -O1 Somewhere between -O0 and -O2.
        -O2 Moderate level of optimization which enables most optimizations.
        -O3 Like -O2, except that it enables optimizations that take longer
        to perform or that may generate larger code (in an attempt to make
        the program run faster).
        -Ofast Enables all the optimizations from -O3 along with other
        agressive optimizations that may violate strict compliance with
        language standards.
        -Os Like -O2 with extra optimizations to reduce code size.
        -Oz Like -Os (and thus -O2), but reduces code size further.
        -Og Like -O1. In future versions, this option might disable
        different optimizations in order to improve debuggability.
        -O Equivalent to -O2.
        -O4 and higher
    """

    return """\
import numpy

#            module name specified by `%%cython_pyximport` magic
#            |        just `modname + ".pyx"`
#            |        |
def make_ext(modname, pyxfilename):
    from setuptools.extension import Extension
    return Extension(modname,
                     sources=[pyxfilename, {files_to_compile}],
                     include_dirs=['.', numpy.get_include()], extra_compile_args=['-O3'])
\
"""

def _create_soln_redun_template_c():
    """
    C language implementation of create_soln_deriv_template()

    Returns
    -------
    string :
        The template string.
    """

    return """\

static double {module}_s(double T, double P, double n[{number_components}]) {{
    double result = -{module}_dgdt(T, P, n);
    return result;
}}

static double {module}_v(double T, double P, double n[{number_components}]) {{
    double result = {module}_dgdp(T, P, n);
    return result;
}}

static double {module}_cv(double T, double P, double n[{number_components}]) {{
    double result = -T*{module}_d2gdt2(T, P, n);
    double dvdt = {module}_d2gdtdp(T, P, n);
    double dvdp = {module}_d2gdp2(T, P, n);
    result += T*dvdt*dvdt/dvdp;
    return result;
}}

static double {module}_cp(double T, double P, double n[{number_components}]) {{
    double result = -T*{module}_d2gdt2(T, P, n);
    return result;
}}

static double {module}_dcpdt(double T, double P, double n[{number_components}]) {{
    double result = -T*{module}_d3gdt3(T, P, n) - {module}_d2gdt2(T, P, n);
    return result;
}}

static double {module}_alpha(double T, double P, double n[{number_components}]) {{
    double result = {module}_d2gdtdp(T, P, n)/{module}_dgdp(T, P, n);
    return result;
}}

static double {module}_beta(double T, double P, double n[{number_components}]) {{
    double result = -{module}_d2gdp2(T, P, n)/{module}_dgdp(T, P, n);
    return result;
}}

static double {module}_K(double T, double P, double n[{number_components}]) {{
    double result = -{module}_dgdp(T, P, n)/{module}_d2gdp2(T, P, n);
    return result;
}}

static double {module}_Kp(double T, double P, double n[{number_components}]) {{
    double result = {module}_dgdp(T, P, n);
    result *= {module}_d3gdp3(T, P, n);
    result /= pow({module}_d2gdp2(T, P, n), 2.0);
    return result - 1.0;
}}

\
"""

def _create_soln_std_state_include_template_c():
    """
    C language implementation of create_soln_std_state_include_template()

    Returns
    -------
    string :
        The template string.

    Notes
    -----
    Structure looks like:
    typedef struct _endmembers {
        const char *(*name) (void);
        const char *(*formula) (void);
        const double (*mw) (void);
        double (*mu0) (double t, double p);
        double (*dmu0dT) (double t, double p);
        double (*dmu0dP) (double t, double p);
        double (*d2mu0dT2) (double t, double p);
        double (*d2mu0dTdP) (double t, double p);
        double (*d2mu0dP2) (double t, double p);
        double (*d3mu0dT3) (double t, double p);
        double (*d3mu0dT2dP) (double t, double p);
        double (*d3mu0dTdP2) (double t, double p);
        double (*d3mu0dP3) (double t, double p);
    } Endmembers;
    static Endmembers endmember[] = {
        {
            Albite_berman_name,
            Albite_berman_formula,
            Albite_berman_mw,
            Albite_berman_g,
            Albite_berman_dgdt,
            Albite_berman_dgdp,
            Albite_berman_d2gdt2,
            Albite_berman_d2gdtdp,
            Albite_berman_d2gdp2,
            Albite_berman_d3gdt3,
            Albite_berman_d3gdt2dp,
            Albite_berman_d3gdtdp2,
            Albite_berman_d3gdp3
        },
        {
            Anorthite_berman_name,
            Anorthite_berman_formula,
            Anorthite_berman_mw,
            Anorthite_berman_g,
            Anorthite_berman_dgdt,
            Anorthite_berman_dgdp,
            Anorthite_berman_d2gdt2,
            Anorthite_berman_d2gdtdp,
            Anorthite_berman_d2gdp2,
            Anorthite_berman_d3gdt3,
            Anorthite_berman_d3gdt2dp,
            Anorthite_berman_d3gdtdp2,
            Anorthite_berman_d3gdp3
        },
        {
            Sanidine_berman_name,
            Sanidine_berman_formula,
            Sanidine_berman_mw,
            Sanidine_berman_g,
            Sanidine_berman_dgdt,
            Sanidine_berman_dgdp,
            Sanidine_berman_d2gdt2,
            Sanidine_berman_d2gdtdp,
            Sanidine_berman_d2gdp2,
            Sanidine_berman_d3gdt3,
            Sanidine_berman_d3gdt2dp,
            Sanidine_berman_d3gdtdp2,
            Sanidine_berman_d3gdp3
        }
    };
    static int nc = (sizeof endmember / sizeof(_endmembers));
    """

    return """\

#include <stdlib.h>
#include <stdio.h>

{code_block_four}

typedef struct _endmembers {{
  const char *(*name) (void);
  const char *(*formula) (void);
  const double (*mw) (void);
  const double *(*elements) (void);
  double (*mu0) (double t, double p);
  double (*dmu0dT) (double t, double p);
  double (*dmu0dP) (double t, double p);
  double (*d2mu0dT2) (double t, double p);
  double (*d2mu0dTdP) (double t, double p);
  double (*d2mu0dP2) (double t, double p);
  double (*d3mu0dT3) (double t, double p);
  double (*d3mu0dT2dP) (double t, double p);
  double (*d3mu0dTdP2) (double t, double p);
  double (*d3mu0dP3) (double t, double p);
}} Endmembers;

static Endmembers endmember[] = {{
{code_block_three}
}};
static int nc = (sizeof endmember / sizeof(struct _endmembers));

static const double R=8.3143;

\
"""

def _create_ordering_gaussj_template_c():
    """
    C language implementation of create_ordering_gaussj_template()

    Returns
    -------
    string :
        The template string.
    """

    return """\

#define SWAP(a,b) {{temp=(a);(a)=(b);(b)=temp;}}

static void gaussj(double a[{NS}][{NS}]) {{
    int indxc[{NS}], indxr[{NS}], ipiv[{NS}];
    int i, icol = -1, irow = -1, j, k, l,ll;
    double big, dum, pivinv, temp;

    for (j=0; j<{NS}; j++) ipiv[j]=0;
    for (i=0; i<{NS}; i++) {{
        big=0.0;
        for (j=0; j<{NS}; j++)
            if (ipiv[j] != 1)
                for (k=0; k<{NS}; k++) {{
                    if (ipiv[k] == 0) {{
                        if (fabs(a[j][k]) >= big) {{
                            big = fabs(a[j][k]);
                            irow = j;
                            icol = k;
                        }}
                    }} else if (ipiv[k] > 1) return;
                }}
        ++(ipiv[icol]);
        if (irow != icol) {{
            for (l=0; l<{NS}; l++) SWAP(a[irow][l],a[icol][l])
                }}
        indxr[i] = irow;
        indxc[i] = icol;
        if (a[icol][icol] == 0.0) return;
        pivinv = 1.0/a[icol][icol];
        a[icol][icol] = 1.0;
        for (l=0; l<{NS}; l++) a[icol][l] *= pivinv;
        for (ll=0; ll<{NS}; ll++)
            if (ll != icol) {{
                dum = a[ll][icol];
                a[ll][icol] = 0.0;
                for (l=0; l<{NS}; l++) a[ll][l] -= a[icol][l]*dum;
            }}
    }}
    for (l=({NS}-1); l>=0; l--) {{
        if (indxr[l] != indxc[l])
            for (k=0; k<{NS}; k++)
                SWAP(a[k][indxr[l]],a[k][indxc[l]]);
    }}
}}
"""

def _create_ordering_code_template_c():
    """
    C language implementation of create_ordering_code_template()

    Returns
    -------
    string :
        The template string.
    """

    return """\

static void order_s(double T, double P, double n[{NC}], double s[{NS}],
    double invd2gds2[{NS}][{NS}]) {{
    double dgds[{NS}], sOld[{NS}];
    int i, j, iter = 0;
{ORDER_CODE_BLOCK_ZERO}
    do {{
        double deltaS[{NS}], steplength;
{ORDER_CODE_BLOCK_TWO}
        for (i=0; i<{NS}; i++) sOld[i] = s[i];

{ORDER_CODE_BLOCK_THREE}
        for (i=0; i<{NS}; i++) {{
            for(j=0; j<{NS}; j++) s[i] += - invd2gds2[i][j]*dgds[j];
            deltaS[i] = s[i] - sOld[i];
        }}

        steplength = 2.0;

{ORDER_CODE_BLOCK_FOUR}
        if (steplength < 1.0) for (i=0; i<{NS}; i++)
            s[i] = sOld[i] + steplength*deltaS[i];
        printf("\\n");
        printf("iter: %3d, step: %13.6e \\n", iter, steplength);
        for (i=0; i<{NS}; i++) {{
            printf(".s%d. %13.6e %13.6e %13.6e %13.6e \\n", i+1, sOld[i], s[i], deltaS[i], dgds[i]);
        }}
        iter++;
    }} while (({ORDER_CODE_BLOCK_ONE}) && (iter < {MAX_ITER}));
}}

static void order_dsdn(double T, double P, double n[{NC}], double s[{NS}],
    double invd2gds2[{NS}][{NS}], double dsdn[{NS}][{NC}]) {{
{ORDER_CODE_BLOCK_FIVE}
    int i,j,k;
    for (i=0; i<{NS}; i++) {{
        for (j=0; j<{NC}; j++) {{
            dsdn[i][j] = 0.0;
            for (k=0; k<{NS}; k++) dsdn[i][j] += - invd2gds2[i][k]*d2gdnds[j][k];
        }}
    }}
}}

static void order_dsdt(double T, double P, double n[{NC}], double s[{NS}],
    double invd2gds2[{NS}][{NS}], double dsdt[{NS}]) {{
    int i,j;
{ORDER_CODE_BLOCK_SIX}
    for (i=0; i<{NS}; i++) {{
        dsdt[i] = 0.0;
        for (j=0; j<{NS}; j++) dsdt[i] += - invd2gds2[i][j]*d2gdsdt[j];
    }}
}}

static void order_dsdp(double T, double P, double n[{NC}], double s[{NS}],
    double invd2gds2[{NS}][{NS}], double dsdp[{NS}]) {{
    int i,j;
{ORDER_CODE_BLOCK_SEVEN}
    for (i=0; i<{NS}; i++) {{
        dsdp[i] = 0.0;
        for (j=0; j<{NS}; j++) dsdp[i] += - invd2gds2[i][j]*d2gdsdp[j];
    }}
}}

static void order_d2sdn2(double T, double P, double n[{NC}], double s[{NS}],
    double invd2gds2[{NS}][{NS}], double d2sdn2[{NS}][{NC}][{NC}]) {{
    double dsdn[{NS}][{NC}], temp[{NS}];
    int i, j, k, l, m, o;
{ORDER_CODE_BLOCK_EIGHT}
    /* compute dsdn matrix */
    for (i=0; i<{NS}; i++) {{
        for (j=0; j<{NC}; j++) {{
        dsdn[i][j] = 0.0;
            for (k=0; k<{NS}; k++) dsdn[i][j] += - invd2gds2[i][k]*d2gdnds[j][k];
        }}
    }}

    /* compute dsdn2 cube */
    for (i=0; i<{NS}; i++) {{
        for (j=0; j<{NC}; j++) {{
            for (k=0; k<{NC}; k++) {{
                for (l=0; l<{NS}; l++) {{
                    temp[l] = d3gdn2ds[j][k][l];
                    for (m=0; m<{NS}; m++) {{
                        temp[l] += d3gdnds2[j][l][m]*dsdn[m][k]
                        + d3gdnds2[k][l][m]*dsdn[m][j];
                        for (o=0; o<{NS}; o++)
                            temp[l] += d3gds3[l][m][o]*dsdn[m][j]*dsdn[o][k];
                    }}
                }}
                d2sdn2[i][j][k] = 0.0;
                for (l=0; l<{NS}; l++) d2sdn2[i][j][k] += - invd2gds2[i][l]*temp[l];
            }}
        }}
    }}
}}

static void order_d2sdndt(double T, double P, double n[{NC}], double s[{NS}],
    double invd2gds2[{NS}][{NS}], double d2sdndt[{NS}][{NC}]) {{
    double dsdn[{NS}][{NC}], dsdt[{NS}], temp[{NS}];
    int i, j, k, l, m;
{ORDER_CODE_BLOCK_NINE}

    /* compute dsdn matrix */
    for (i=0; i<{NS}; i++) {{
        for (j=0; j<{NC}; j++) {{
            dsdn[i][j] = 0.0;
            for (k=0; k<{NS}; k++) dsdn[i][j] += - invd2gds2[i][k]*d2gdnds[j][k];
        }}
    }}

    /* compute dsdt vector */
    for (i=0; i<{NS}; i++) {{
        dsdt[i] = 0.0;
        for (j=0; j<{NS}; j++) dsdt[i] += - invd2gds2[i][j]*d2gdsdt[j];
    }}

    /* compute d2sdndt matrix */
    for (i=0; i<{NS}; i++) {{
        for (j=0; j<{NC}; j++) {{
            for (k=0; k<{NS}; k++) {{
                temp[k] = d3gdndsdt[j][k];
                for (l=0; l<{NS}; l++) {{
                    temp[k] += d3gdnds2[j][k][l]*dsdt[l] + d3gds2dt[k][l]*dsdn[l][j];
                    for (m=0; m<{NS}; m++) temp[k] += d3gds3[k][l][m]*dsdn[l][j]*dsdt[m];
                }}
            }}
            d2sdndt[i][j] = 0.0;
            for (k=0; k<{NS}; k++) d2sdndt[i][j] += - invd2gds2[i][k]*temp[k];
        }}
    }}
}}

static void order_d2sdndp(double T, double P, double n[{NC}], double s[{NS}],
    double invd2gds2[{NS}][{NS}], double d2sdndp[{NS}][{NC}]) {{
    double dsdn[{NS}][{NC}], dsdp[{NS}], temp[{NS}];
    int i, j, k, l, m;
{ORDER_CODE_BLOCK_TEN}

    /* compute dsdn matrix */
    for (i=0; i<{NS}; i++) {{
        for (j=0; j<{NC}; j++) {{
            dsdn[i][j] = 0.0;
            for (k=0; k<{NS}; k++) dsdn[i][j] += - invd2gds2[i][k]*d2gdnds[j][k];
        }}
    }}

    /* compute dsdp vector */
    for (i=0; i<{NS}; i++) {{
        dsdp[i] = 0.0;
        for (j=0; j<{NS}; j++) dsdp[i] += - invd2gds2[i][j]*d2gdsdp[j];
    }}

    /* compute d2sdndp matrix */
    for (i=0; i<{NS}; i++) {{
        for (j=0; j<{NC}; j++) {{
            for (k=0; k<{NS}; k++) {{
                temp[k] = d3gdndsdp[j][k];
                for (l=0; l<{NS}; l++) {{
                    temp[k] += d3gdnds2[j][k][l]*dsdp[l] + d3gds2dp[k][l]*dsdn[l][j];
                    for (m=0; m<{NS}; m++) temp[k] += d3gds3[k][l][m]*dsdn[l][j]*dsdp[m];
                }}
            }}
            d2sdndp[i][j] = 0.0;
            for (k=0; k<{NS}; k++) d2sdndp[i][j] += - invd2gds2[i][k]*temp[k];
        }}
    }}
}}

static void order_d2sdt2(double T, double P, double n[{NC}], double s[{NS}],
    double invd2gds2[{NS}][{NS}], double d2sdt2[{NS}]) {{
    double dsdt[{NS}], temp[{NS}];
    int i, j, k, l;
{ORDER_CODE_BLOCK_ELEVEN}

    /* compute dsdt vector */
    for (i=0; i<{NS}; i++) {{
        dsdt[i] = 0.0;
        for (j=0; j<{NS}; j++) dsdt[i] += - invd2gds2[i][j]*d2gdsdt[j];
    }}

    /* compute d2sdt2 vector */
    for (i=0; i<{NS}; i++) {{
        for (j=0; j<{NS}; j++) {{
            temp[j] = d3gdsdt2[j];
            for (k=0; k<{NS}; k++) {{
                temp[j] +=  2.0*d3gds2dt[j][k]*dsdt[k];
                for (l=0; l<{NS}; l++) temp[j] += d3gds3[j][k][l]*dsdt[k]*dsdt[l];
            }}
        }}
        d2sdt2[i] = 0.0;
        for (j=0; j<{NS}; j++) d2sdt2[i] += - invd2gds2[i][j]*temp[j];
    }}
}}

static void order_d2sdtdp(double T, double P, double n[{NC}], double s[{NS}],
    double invd2gds2[{NS}][{NS}], double d2sdtdp[{NS}]) {{
    double dsdt[{NS}], dsdp[{NS}], temp[{NS}];
    int i, j, k, l;
{ORDER_CODE_BLOCK_TWELVE}

    /* compute dsdt vector */
    for (i=0; i<{NS}; i++) {{
        dsdt[i] = 0.0;
        for (j=0; j<{NS}; j++) dsdt[i] += - invd2gds2[i][j]*d2gdsdt[j];
    }}

    /* compute dsdp vector */
    for (i=0; i<{NS}; i++) {{
        dsdp[i] = 0.0;
        for (j=0; j<{NS}; j++) dsdp[i] += - invd2gds2[i][j]*d2gdsdp[j];
    }}

    /* compute d2sdtdp vector */
    for (i=0; i<{NS}; i++) {{
        for (j=0; j<{NS}; j++) {{
            temp[j] = d3gdsdtdp[j];
            for (k=0; k<{NS}; k++) {{
                temp[j] += d3gds2dt[j][k]*dsdp[k] + d3gds2dp[j][k]*dsdt[k];
                for (l=0; l<{NS}; l++) temp[j] += d3gds3[j][k][l]*dsdt[k]*dsdp[l];
            }}
        }}
        d2sdtdp[i] = 0.0;
        for (j=0; j<{NS}; j++) d2sdtdp[i] += - invd2gds2[i][j]*temp[j];
    }}
}}

static void order_d2sdp2(double T, double P, double n[{NC}], double s[{NS}],
    double invd2gds2[{NS}][{NS}], double d2sdp2[{NS}]) {{
    double dsdp[{NS}], temp[{NS}];
    int i, j, k, l;
{ORDER_CODE_BLOCK_THIRTEEN}

    /* compute dsdp vector */
    for (i=0; i<{NS}; i++) {{
        dsdp[i] = 0.0;
        for (j=0; j<{NS}; j++) dsdp[i] += - invd2gds2[i][j]*d2gdsdp[j];
    }}

    /* compute d2sdp2 vector */
    for (i=0; i<{NS}; i++) {{
        for (j=0; j<{NS}; j++) {{
            temp[j] = d3gdsdp2[j];
            for (k=0; k<{NS}; k++) {{
                temp[j] +=  2.0*d3gds2dp[j][k]*dsdp[k];
                for (l=0; l<{NS}; l++) temp[j] += d3gds3[j][k][l]*dsdp[k]*dsdp[l];
            }}
        }}
        d2sdp2[i] = 0.0;
        for (j=0; j<{NS}; j++) d2sdp2[i] += - invd2gds2[i][j]*temp[j];
    }}
}}

    \
    """

def _create_complx_soln_calc_template_c():
    """
    C language implementation of create_soln_calc_template()

    Returns
    -------
    string :
        The template string.
    """

    return """\

static double *retrieveGuess(double T, double P, double n[{number_components}]) {{
{moles_assign}
    static double s[{number_ordering}];
    static double nOld[{number_components}];
    static int reset = 1;
    for (int i=0; i<{number_ordering}; i++) if (isnan(s[i])) reset = 1;
    for (int i=0; i<{number_components}; i++) if (nOld[i] != n[i]) {{
        reset = 1;
        break;
    }}
    if (reset) {{
{order_initial_guess}
        for (int i=0; i<{number_components}; i++) nOld[i] = n[i];
        reset = 0;
    }}
    return s;
}}

static double {module}_g(double T, double P, double n[{number_components}]) {{
{moles_assign}
    double invd2gds2[{number_ordering}][{number_ordering}];
    double *s = retrieveGuess(T, P, n);
    order_s(T, P, n, s, invd2gds2);
{order_assign}
    double result = {g_code};
    return result;
}}

static double {module}_dgdt(double T, double P, double n[{number_components}]) {{
{moles_assign}
    double invd2gds2[{number_ordering}][{number_ordering}];
    double *s = retrieveGuess(T, P, n);
    order_s(T, P, n, s, invd2gds2);
{order_assign}
    double result = {dgdt_code};
    return result;
}}

static double {module}_dgdp(double T, double P, double n[{number_components}]) {{
{moles_assign}
    double invd2gds2[{number_ordering}][{number_ordering}];
    double *s = retrieveGuess(T, P, n);
    order_s(T, P, n, s, invd2gds2);
{order_assign}
    double result = {dgdp_code};
    return result;
}}

static void {module}_dgdn(double T, double P, double n[{number_components}],
    double result[{number_components}]) {{
    {moles_assign}
    double invd2gds2[{number_ordering}][{number_ordering}];
    double dgdn[{number_components}];
    int i;
    double *s = retrieveGuess(T, P, n);
    order_s(T, P, n, s, invd2gds2);
{order_assign}
{dgdn_code}
    for (i=0; i<{number_components}; i++) result[i] = dgdn[i];
}}

static double {module}_d2gdt2(double T, double P, double n[{number_components}]) {{
{moles_assign}
    double invd2gds2[{number_ordering}][{number_ordering}];
    double dsdt[{number_ordering}];
    double d2gds2[{number_ordering}][{number_ordering}], d2gdsdt[{number_ordering}];
    int i,j;
    double *s = retrieveGuess(T, P, n);
    order_s(T, P, n, s, invd2gds2);
    order_dsdt(T, P, n, s, invd2gds2, dsdt);
{order_assign}
{fillD2GDS2}
{fillD2GDSDT}
    double result = {d2gdt2_code};
    for (i=0; i<{number_ordering}; i++) {{
        result += 2.0*d2gdsdt[i]*dsdt[i];
        for (j=0; j<{number_ordering}; j++) result += d2gds2[i][j]*dsdt[i]*dsdt[j];
    }}
    return result;
}}

static double {module}_d2gdtdp(double T, double P, double n[{number_components}]) {{
{moles_assign}
    double invd2gds2[{number_ordering}][{number_ordering}];
    double dsdt[{number_ordering}], dsdp[{number_ordering}];
    double d2gds2[{number_ordering}][{number_ordering}], d2gdsdt[{number_ordering}],
           d2gdsdp[{number_ordering}];
    int i,j;
    double *s = retrieveGuess(T, P, n);
    order_s(T, P, n, s, invd2gds2);
    order_dsdt(T, P, n, s, invd2gds2, dsdt);
    order_dsdp(T, P, n, s, invd2gds2, dsdp);
{order_assign}
{fillD2GDS2}
{fillD2GDSDT}
{fillD2GDSDP}
    double result = {d2gdtdp_code};
    for (i=0; i<{number_ordering}; i++) {{
        result += d2gdsdt[i]*dsdp[i] + d2gdsdp[i]*dsdt[i];
        for (j=0; j<{number_ordering}; j++) result += d2gds2[i][j]*dsdt[i]*dsdp[j];
    }}
    return result;
}}

static void {module}_d2gdndt(double T, double P, double n[{number_components}],
    double result[{number_components}]) {{
    {moles_assign}
    double invd2gds2[{number_ordering}][{number_ordering}];
    double dsdt[{number_ordering}], dsdn[{number_ordering}][{number_components}];
    double d2gdnds[{number_components}][{number_ordering}],
           d2gdndt[{number_components}], d2gdsdt[{number_ordering}],
           d2gds2[{number_ordering}][{number_ordering}];
    int i,k,l;
    double *s = retrieveGuess(T, P, n);
    order_s(T, P, n, s, invd2gds2);
    order_dsdt(T, P, n, s, invd2gds2, dsdt);
    order_dsdn(T, P, n, s, invd2gds2, dsdn);
{order_assign}
{fillD2GDNDS}
{fillD2GDNDT}
{fillD2GDS2}
{fillD2GDSDT}
    for (i=0; i<{number_components}; i++) {{
        result[i] = d2gdndt[i];
        for (k=0; k<{number_ordering}; k++) {{
            result[i] += d2gdnds[i][k]*dsdt[k] + d2gdsdt[k]*dsdn[k][i];
            for (l=0; l<{number_ordering}; l++) result[i] += d2gds2[k][l]*dsdt[k]*dsdn[l][i];
        }}
    }}
}}

static double {module}_d2gdp2(double T, double P, double n[{number_components}]) {{
{moles_assign}
    double invd2gds2[{number_ordering}][{number_ordering}];
    double dsdp[{number_ordering}];
    double d2gds2[{number_ordering}][{number_ordering}],d2gdsdp[{number_ordering}];
    int i,j;
    double *s = retrieveGuess(T, P, n);
    order_s(T, P, n, s, invd2gds2);
    order_dsdp(T, P, n, s, invd2gds2, dsdp);
{order_assign}
{fillD2GDS2}
{fillD2GDSDP}
    double result = {d2gdp2_code};
    for (i=0; i<{number_ordering}; i++) {{
        result += 2.0*d2gdsdp[i]*dsdp[i];
        for (j=0; j<{number_ordering}; j++) result += d2gds2[i][j]*dsdp[i]*dsdp[j];
    }}
    return result;
}}

static void {module}_d2gdndp(double T, double P, double n[{number_components}],
    double result[{number_components}]) {{
    {moles_assign}
    double invd2gds2[{number_ordering}][{number_ordering}];
    double dsdp[{number_ordering}], dsdn[{number_ordering}][{number_components}];
    double d2gdnds[{number_components}][{number_ordering}],
           d2gdndp[{number_components}], d2gdsdp[{number_ordering}],
           d2gds2[{number_ordering}][{number_ordering}];
    int i,k,l;
    double *s = retrieveGuess(T, P, n);
    order_s(T, P, n, s, invd2gds2);
    order_dsdp(T, P, n, s, invd2gds2, dsdp);
    order_dsdn(T, P, n, s, invd2gds2, dsdn);
{order_assign}
{fillD2GDNDS}
{fillD2GDNDP}
{fillD2GDS2}
{fillD2GDSDP}
    for (i=0; i<{number_components}; i++) {{
        result[i] = d2gdndp[i];
        for (k=0; k<{number_ordering}; k++) {{
            result[i] += d2gdnds[i][k]*dsdp[k] + d2gdsdp[k]*dsdn[k][i];
            for (l=0; l<{number_ordering}; l++) result[i] += d2gds2[k][l]*dsdp[k]*dsdn[l][i];
        }}
    }}
}}

static void {module}_d2gdn2(double T, double P, double n[{number_components}],
    double result[{number_symmetric_hessian_terms}]) {{
    {moles_assign}
    double invd2gds2[{number_ordering}][{number_ordering}];
    double dsdn[{number_ordering}][{number_components}];
    double d2gdn2[{number_components}][{number_components}],
           d2gdnds[{number_components}][{number_ordering}],
           d2gds2[{number_ordering}][{number_ordering}];
    int i,j,k,l,m;
    double *s = retrieveGuess(T, P, n);
    order_s(T, P, n, s, invd2gds2);
    order_dsdn(T, P, n, s, invd2gds2, dsdn);
{order_assign}
{fillD2GDN2}
{fillD2GDNDS}
{fillD2GDS2}
    m = 0;
    for (i=0; i<{number_components}; i++) {{
        for (j=i; j<{number_components}; j++) {{
            result[m] = d2gdn2[i][j];
            for (k=0; k<{number_ordering}; k++) {{
                result[m] += d2gdnds[i][k]*dsdn[k][j] + d2gdnds[j][k]*dsdn[k][i];
                for (l=0; l<{number_ordering}; l++) result[m] += d2gds2[k][l]*dsdn[k][i]*dsdn[l][j];
            }}
            m += 1;
        }}
    }}
}}

static double {module}_d3gdt3(double T, double P, double n[{number_components}]) {{
{moles_assign}
    double invd2gds2[{number_ordering}][{number_ordering}];
    double dsdt[{number_ordering}], d2sdt2[{number_ordering}];
    double d2gds2[{number_ordering}][{number_ordering}], d2gdsdt[{number_ordering}],
           d3gds3[{number_ordering}][{number_ordering}][{number_ordering}],
           d3gds2dt[{number_ordering}][{number_ordering}],
           d3gdsdt2[{number_ordering}];
    int i,j,k;
    double *s = retrieveGuess(T, P, n);
    order_s(T, P, n, s, invd2gds2);
    order_dsdt(T, P, n, s, invd2gds2, dsdt);
    order_d2sdt2(T, P, n, s, invd2gds2, d2sdt2);
{order_assign}
{fillD2GDS2}
{fillD2GDSDT}
{fillD3GDS3}
{fillD3GDS2DT}
{fillD3GDSDT2}
    double result = {d3gdt3_code};
    for (i=0; i<{number_ordering}; i++) {{
        result += 3.0*d3gdsdt2[i]*dsdt[i] + 3.0*d2gdsdt[i]*d2sdt2[i];
        for (j=0; j<{number_ordering}; j++) {{
            result += 3.0*d2gds2[i][j]*dsdt[i]*d2sdt2[j]
                    + 3.0*d3gds2dt[i][j]*dsdt[i]*dsdt[j];
            for (k=0; k<{number_ordering}; k++) result += d3gds3[i][j][k]*dsdt[i]*dsdt[j]*dsdt[k];
        }}
    }}
    return result;
}}

static double {module}_d3gdt2dp(double T, double P, double n[{number_components}]) {{
{moles_assign}
    double invd2gds2[{number_ordering}][{number_ordering}];
    double dsdt[{number_ordering}], dsdp[{number_ordering}];
    double d2sdt2[{number_ordering}], d2sdtdp[{number_ordering}];
    double d2gds2[{number_ordering}][{number_ordering}],
           d2gdsdt[{number_ordering}], d2gdsdp[{number_ordering}],
           d3gds3[{number_ordering}][{number_ordering}][{number_ordering}],
           d3gds2dt[{number_ordering}][{number_ordering}],
           d3gds2dp[{number_ordering}][{number_ordering}],
           d3gdsdt2[{number_ordering}], d3gdsdtdp[{number_ordering}];
    int i,j,k;
    double *s = retrieveGuess(T, P, n);
    order_s(T, P, n, s, invd2gds2);
    order_dsdt(T, P, n, s, invd2gds2, dsdt);
    order_dsdp(T, P, n, s, invd2gds2, dsdp);
    order_d2sdt2(T, P, n, s, invd2gds2, d2sdt2);
    order_d2sdtdp(T, P, n, s, invd2gds2, d2sdtdp);
{order_assign}
{fillD2GDS2}
{fillD2GDSDT}
{fillD2GDSDP}
{fillD3GDS3}
{fillD3GDS2DT}
{fillD3GDS2DP}
{fillD3GDSDT2}
{fillD3GDSDTDP}
    double result = {d3gdt2dp_code};
    for (i=0; i<{number_ordering}; i++) {{
        result += d3gdsdt2[i]*dsdp[i] + 2.0*d2gdsdt[i]*d2sdtdp[i]
                + d2gdsdp[i]*d2sdt2[i] + 2.0*d3gdsdtdp[i]*dsdt[i];
        for (j=0; j<{number_ordering}; j++) {{
            result += 2.0*d3gds2dt[i][j]*dsdt[i]*dsdp[j]
                    + d2gds2[i][j]*d2sdt2[i]*dsdp[j]
                    + 2.0*d2gds2[i][j]*dsdt[i]*d2sdtdp[j]
                    + d3gds2dp[i][j]*dsdt[i]*dsdt[j];
            for (k=0; k<{number_ordering}; k++) result += d3gds3[i][j][k]*dsdt[i]*dsdt[j]*dsdp[k];
        }}
    }}
    return result;
}}

static void {module}_d3gdndt2(double T, double P, double n[{number_components}],
    double result[{number_components}]) {{
    {moles_assign}
    double invd2gds2[{number_ordering}][{number_ordering}];
    double dsdt[{number_ordering}], dsdn[{number_ordering}][{number_components}];
    double d2sdt2[{number_ordering}], d2sdndt[{number_ordering}][{number_components}];
    double d2gdnds[{number_components}][{number_ordering}],
           d2gdndt[{number_components}], d2gdsdt[{number_ordering}],
           d2gds2[{number_ordering}][{number_ordering}],
           d3gdnds2[{number_components}][{number_ordering}][{number_ordering}],
           d3gdndsdt[{number_components}][{number_ordering}],
           d3gdndt2[{number_components}], d3gdsdt2[{number_ordering}],
           d3gds3[{number_ordering}][{number_ordering}][{number_ordering}],
           d3gds2dt[{number_ordering}][{number_ordering}];
    int i,j,k,l;
    double *s = retrieveGuess(T, P, n);
    order_s(T, P, n, s, invd2gds2);
    order_dsdt(T, P, n, s, invd2gds2, dsdt);
    order_dsdn(T, P, n, s, invd2gds2, dsdn);
    order_d2sdt2(T, P, n, s, invd2gds2, d2sdt2);
    order_d2sdndt(T, P, n, s, invd2gds2, d2sdndt);
{order_assign}
{fillD2GDNDS}
{fillD2GDNDT}
{fillD2GDS2}
{fillD2GDSDT}
{fillD3GDNDS2}
{fillD3GDNDSDT}
{fillD3GDNDT2}
{fillD3GDS3}
{fillD3GDS2DT}
{fillD3GDSDT2}
    for (i=0; i<{number_components}; i++) {{
        for (j=0,result[i]=d3gdndt2[i]; j<{number_ordering}; j++) {{
            result[i] += d3gdsdt2[j]*dsdn[j][i] + 2.0*d2gdsdt[j]*d2sdndt[j][i]
                       + 2.0*d3gdndsdt[i][j]*dsdt[j] + d2gdnds[i][j]*d2sdt2[j];
            for (k=0; k<{number_ordering}; k++) {{
                result[i] += d3gdnds2[i][j][k]*dsdt[j]*dsdt[k]
                           + 2.0*d2gds2[j][k]*dsdt[j]*d2sdndt[k][i]
                           + 2.0*d3gds2dt[j][k]*dsdn[j][i]*dsdt[k]
                           + d2gds2[j][k]*dsdn[j][i]*d2sdt2[k];
                for (l=0; l<{number_ordering}; l++)
                    result[i] += d3gds3[j][k][l]*dsdn[j][i]*dsdt[k]*dsdt[l];
            }}
        }}
    }}
}}

static double {module}_d3gdtdp2(double T, double P, double n[{number_components}]) {{
{moles_assign}
    double invd2gds2[{number_ordering}][{number_ordering}];
    double dsdt[{number_ordering}], dsdp[{number_ordering}];
    double d2sdp2[{number_ordering}], d2sdtdp[{number_ordering}];
    double d2gds2[{number_ordering}][{number_ordering}], d2gdsdt[{number_ordering}],
           d2gdsdp[{number_ordering}],
           d3gds3[{number_ordering}][{number_ordering}][{number_ordering}],
           d3gds2dt[{number_ordering}][{number_ordering}],
           d3gds2dp[{number_ordering}][{number_ordering}],
           d3gdsdt2[{number_ordering}], d3gdsdtdp[{number_ordering}],
           d3gdsdp2[{number_ordering}];
    int i,j,k;
    double *s = retrieveGuess(T, P, n);
    order_s(T, P, n, s, invd2gds2);
    order_dsdt(T, P, n, s, invd2gds2, dsdt);
    order_dsdp(T, P, n, s, invd2gds2, dsdp);
    order_d2sdp2(T, P, n, s, invd2gds2, d2sdp2);
    order_d2sdtdp(T, P, n, s, invd2gds2, d2sdtdp);
{order_assign}
{fillD2GDS2}
{fillD2GDSDT}
{fillD2GDSDP}
{fillD3GDS3}
{fillD3GDS2DT}
{fillD3GDS2DP}
{fillD3GDSDT2}
{fillD3GDSDTDP}
{fillD3GDSDP2}
    double result = {d3gdtdp2_code};
    for (i=0; i<{number_ordering}; i++) {{
        result += 2.0*d3gdsdtdp[i]*dsdp[i] + d2gdsdt[i]*d2sdp2[i]
                + 2.0*d2gdsdp[i]*d2sdtdp[i] + d3gdsdp2[i]*dsdt[i];
        for (j=0; j<{number_ordering}; j++) {{
            result += 2.0*d3gds2dp[i][j]*dsdt[i]*dsdp[j]
                    + d2gds2[i][j]*dsdt[i]*d2sdp2[j]
                    + 2.0*d2gds2[i][j]*d2sdtdp[i]*dsdp[j]
                    + d3gds2dt[i][j]*dsdp[i]*dsdp[j];
            for (k=0; k<{number_ordering}; k++) result += d3gds3[i][j][k]*dsdt[i]*dsdp[j]*dsdp[k];
        }}
    }}
    return result;
}}

static void {module}_d3gdndtdp(double T, double P, double n[{number_components}],
    double result[{number_components}]) {{
    {moles_assign}
    double invd2gds2[{number_ordering}][{number_ordering}];
    double dsdt[{number_ordering}], dsdp[{number_ordering}], dsdn[{number_ordering}][{number_components}];
    double d2sdtdp[{number_ordering}], d2sdndt[{number_ordering}][{number_components}];
    double d2sdndp[{number_ordering}][{number_components}];
    double d2gdnds[{number_components}][{number_ordering}],
           d2gds2[{number_ordering}][{number_ordering}],
           d2gdsdt[{number_ordering}], d2gdsdp[{number_ordering}],
           d3gdnds2[{number_components}][{number_ordering}][{number_ordering}],
           d3gdndsdt[{number_components}][{number_ordering}],
           d3gdndsdp[{number_components}][{number_ordering}],
           d3gdndtdp[{number_components}],
           d3gds3[{number_ordering}][{number_ordering}][{number_ordering}],
           d3gds2dt[{number_ordering}][{number_ordering}],
           d3gdsdtdp[{number_ordering}],
           d3gds2dp[{number_ordering}][{number_ordering}];
    int i,j,k,l;
    double *s = retrieveGuess(T, P, n);
    order_s(T, P, n, s, invd2gds2);
    order_dsdt(T, P, n, s, invd2gds2, dsdt);
    order_dsdn(T, P, n, s, invd2gds2, dsdn);
    order_d2sdtdp(T, P, n, s, invd2gds2, d2sdtdp);
    order_d2sdndt(T, P, n, s, invd2gds2, d2sdndt);
{order_assign}
{fillD2GDNDS}
{fillD2GDS2}
{fillD2GDSDT}
{fillD2GDSDP}
{fillD3GDNDS2}
{fillD3GDNDSDT}
{fillD3GDNDSDP}
{fillD3GDNDTDP}
{fillD3GDS3}
{fillD3GDS2DT}
{fillD3GDSDTDP}
{fillD3GDS2DP}
    for (i=0; i<{number_components}; i++) {{
        for (j=0,result[i]=d3gdndtdp[i]; j<{number_ordering}; j++) {{
            result[i] += d3gdsdtdp[j]*dsdn[j][i] + d2gdsdt[j]*d2sdndp[j][i]
                       + d3gdndsdt[i][j]*dsdp[j] + d2gdnds[i][j]*d2sdtdp[j]
                       + d3gdndsdp[i][j]*dsdt[j] + d2gdsdp[j]*d2sdndt[j][i];
            for (k=0; k<{number_ordering}; k++) {{
                result[i] += d3gdnds2[i][j][k]*dsdt[j]*dsdp[k]
                           + d2gds2[j][k]*dsdt[j]*d2sdndp[k][i]
                           + d2gds2[j][k]*dsdp[j]*d2sdndt[k][i]
                           + d3gds2dt[j][k]*dsdn[j][i]*dsdp[k]
                           + d3gds2dp[j][k]*dsdn[j][i]*dsdt[k]
                           + d2gds2[j][k]*dsdn[j][i]*d2sdtdp[k];
                for (l=0; l<{number_ordering}; l++)
                    result[i] += d3gds3[j][k][l]*dsdn[j][i]*dsdt[k]*dsdp[l];
            }}
        }}
    }}
}}

static void {module}_d3gdn2dt(double T, double P, double n[{number_components}],
    double result[{number_symmetric_hessian_terms}]) {{
    {moles_assign}
    double invd2gds2[{number_ordering}][{number_ordering}];
    double dsdn[{number_ordering}][{number_components}], dsdt[{number_ordering}];
    double d2sdndt[{number_ordering}][{number_components}];
    double d2sdn2[{number_ordering}][{number_components}][{number_components}];
    double d2gdnds[{number_components}][{number_ordering}],
           d2gds2[{number_ordering}][{number_ordering}], d2gdsdt[{number_ordering}],
           d3gdn2ds[{number_components}][{number_components}][{number_ordering}],
           d3gdn2dt[{number_components}][{number_components}],
           d3gdnds2[{number_components}][{number_ordering}][{number_ordering}],
           d3gdndsdt[{number_components}][{number_ordering}],
           d3gds3[{number_ordering}][{number_ordering}][{number_ordering}],
           d3gds2dt[{number_ordering}][{number_ordering}];
    int i,j,k,l,m,o;
    double *s = retrieveGuess(T, P, n);
    order_s(T, P, n, s, invd2gds2);
    order_dsdn(T, P, n, s, invd2gds2, dsdn);
    order_dsdt(T, P, n, s, invd2gds2, dsdt);
    order_d2sdn2(T, P, n, s, invd2gds2, d2sdn2);
    order_d2sdndt(T, P, n, s, invd2gds2, d2sdndt);
{order_assign}
{fillD2GDNDS}
{fillD2GDS2}
{fillD2GDSDT}
{fillD3GDN2DS}
{fillD3GDN2DT}
{fillD3GDNDS2}
{fillD3GDNDSDT}
{fillD3GDS3}
{fillD3GDS2DT}
    o = 0;
    for (i=0; i<{number_components}; i++) {{
        for (j=i; j<{number_components}; j++) {{
            result[o] = d3gdn2dt[i][j];
            for (k=0; k<{number_ordering}; k++) {{
                result[o] += d3gdn2ds[i][j][k]*dsdt[k]
                              + d3gdndsdt[i][k]*dsdn[k][j]
                              + d3gdndsdt[j][k]*dsdn[k][i]
                              + d2gdsdt[k]*d2sdn2[k][i][j]
                              + d2gdnds[i][k]*d2sdndt[k][j]
                              + d2gdnds[j][k]*d2sdndt[k][i];
                for (l=0; l<{number_ordering}; l++) {{
                    result[o] += d3gdnds2[i][k][l]*dsdn[k][j]*dsdt[l]
                                  + d3gdnds2[j][k][l]*dsdn[k][i]*dsdt[l]
                                  + d2gds2[k][l]*d2sdn2[k][i][j]*dsdt[l]
                                  + d3gds2dt[k][l]*dsdn[k][i]*dsdn[l][j]
                                  + d2gds2[k][l]*dsdn[k][i]*d2sdndt[l][j]
                                  + d2gds2[k][l]*dsdn[k][j]*d2sdndt[l][i];
                    for (m=0; m<{number_ordering}; m++)
                        result[o] += d3gds3[k][l][m]*dsdn[k][i]*dsdn[l][j]*dsdt[m];
                }}
            }}
            o += 1;
        }}
    }}
}}

static double {module}_d3gdp3(double T, double P, double n[{number_components}]) {{
{moles_assign}
    double invd2gds2[{number_ordering}][{number_ordering}];
    double dsdp[{number_ordering}], d2sdp2[{number_ordering}];
    double d2gds2[{number_ordering}][{number_ordering}], d2gdsdp[{number_ordering}],
           d3gds3[{number_ordering}][{number_ordering}][{number_ordering}],
           d3gds2dp[{number_ordering}][{number_ordering}],
           d3gdsdp2[{number_ordering}];
    int i,j,k;
    double *s = retrieveGuess(T, P, n);
    order_s(T, P, n, s, invd2gds2);
    order_dsdp(T, P, n, s, invd2gds2, dsdp);
    order_d2sdp2(T, P, n, s, invd2gds2, d2sdp2);
{order_assign}
{fillD2GDS2}
{fillD2GDSDP}
{fillD3GDS3}
{fillD3GDS2DP}
{fillD3GDSDP2}
    double result = {d3gdp3_code};
    for (i=0; i<{number_ordering}; i++) {{
        result += 3.0*d3gdsdp2[i]*dsdp[i] + 3.0*d2gdsdp[i]*d2sdp2[i];
        for (j=0; j<{number_ordering}; j++) {{
            result += 3.0*d2gds2[i][j]*dsdp[i]*d2sdp2[j]
                    + 3.0*d3gds2dp[i][j]*dsdp[i]*dsdp[j];
            for (k=0; k<{number_ordering}; k++) result += d3gds3[i][j][k]*dsdp[i]*dsdp[j]*dsdp[k];
        }}
    }}
    return result;
}}

static void {module}_d3gdndp2(double T, double P, double n[{number_components}],
    double result[{number_components}]) {{
    {moles_assign}
    double invd2gds2[{number_ordering}][{number_ordering}];
    double dsdp[{number_ordering}], dsdn[{number_ordering}][{number_components}];
    double d2sdp2[{number_ordering}], d2sdndp[{number_ordering}][{number_components}];
    double d2gdnds[{number_components}][{number_ordering}],
           d2gdndp[{number_components}], d3gdndp2[{number_components}],
           d2gds2[{number_ordering}][{number_ordering}], d2gdsdp[{number_ordering}],
           d3gdnds2[{number_components}][{number_ordering}][{number_ordering}],
           d3gdndsdp[{number_components}][{number_ordering}],
           d3gds3[{number_ordering}][{number_ordering}][{number_ordering}],
           d3gds2dp[{number_ordering}][{number_ordering}],
           d3gdsdp2[{number_ordering}];
    int i,j,k,l;
    double *s = retrieveGuess(T, P, n);
    order_s(T, P, n, s, invd2gds2);
    order_dsdp(T, P, n, s, invd2gds2, dsdp);
    order_dsdn(T, P, n, s, invd2gds2, dsdn);
    order_d2sdp2(T, P, n, s, invd2gds2, d2sdp2);
    order_d2sdndp(T, P, n, s, invd2gds2, d2sdndp);
{order_assign}
{fillD2GDNDS}
{fillD2GDNDP}
{fillD2GDS2}
{fillD2GDSDP}
{fillD3GDNDS2}
{fillD3GDNDSDP}
{fillD3GDNDP2}
{fillD3GDS3}
{fillD3GDS2DP}
{fillD3GDSDP2}
    for (i=0; i<{number_components}; i++) {{
        for (j=0,result[i]=d3gdndp2[i]; j<{number_ordering}; j++) {{
            result[i] += d3gdsdp2[j]*dsdn[j][i] + 2.0*d2gdsdp[j]*d2sdndp[j][i]
                       + 2.0*d3gdndsdp[i][j]*dsdp[j] + d2gdnds[i][j]*d2sdp2[j];
            for (k=0; k<{number_ordering}; k++) {{
                result[i] += d3gdnds2[i][j][k]*dsdp[j]*dsdp[k]
                           + 2.0*d2gds2[j][k]*dsdp[j]*d2sdndp[k][i]
                           + 2.0*d3gds2dp[j][k]*dsdn[j][i]*dsdp[k]
                           + d2gds2[j][k]*dsdn[j][i]*d2sdp2[k];
                for (l=0; l<{number_ordering}; l++)
                    result[i] += d3gds3[j][k][l]*dsdn[j][i]*dsdp[k]*dsdp[l];
            }}
        }}
    }}
}}

static void {module}_d3gdn2dp(double T, double P, double n[{number_components}],
    double result[{number_symmetric_hessian_terms}]) {{
    {moles_assign}
    double invd2gds2[{number_ordering}][{number_ordering}];
    double dsdn[{number_ordering}][{number_components}], dsdp[{number_ordering}];
    double d2sdndp[{number_ordering}][{number_components}];
    double d2sdn2[{number_ordering}][{number_components}][{number_components}];
    double d2gdnds[{number_components}][{number_ordering}],
           d2gds2[{number_ordering}][{number_ordering}], d2gdsdp[{number_ordering}],
           d3gdn2ds[{number_components}][{number_components}][{number_ordering}],
           d3gdn2dp[{number_components}][{number_components}],
           d3gdnds2[{number_components}][{number_ordering}][{number_ordering}],
           d3gdndsdp[{number_components}][{number_ordering}],
           d3gds3[{number_ordering}][{number_ordering}][{number_ordering}],
           d3gds2dp[{number_ordering}][{number_ordering}];
    int i,j,k,l,m,o;
    double *s = retrieveGuess(T, P, n);
    order_s(T, P, n, s, invd2gds2);
    order_dsdn(T, P, n, s, invd2gds2, dsdn);
    order_dsdp(T, P, n, s, invd2gds2, dsdp);
    order_d2sdn2(T, P, n, s, invd2gds2, d2sdn2);
    order_d2sdndp(T, P, n, s, invd2gds2, d2sdndp);
{order_assign}
{fillD2GDNDS}
{fillD2GDS2}
{fillD2GDSDP}
{fillD3GDN2DS}
{fillD3GDN2DP}
{fillD3GDNDS2}
{fillD3GDNDSDP}
{fillD3GDS3}
{fillD3GDS2DP}
    o = 0;
    for (i=0; i<{number_components}; i++) {{
        for (j=i; j<{number_components}; j++) {{
            result[o] = d3gdn2dp[i][j];
            for (k=0; k<{number_ordering}; k++) {{
                result[o] += d3gdn2ds[i][j][k]*dsdp[k]
                              + d3gdndsdp[i][k]*dsdn[k][j]
                              + d3gdndsdp[j][k]*dsdn[k][i]
                              + d2gdsdp[k]*d2sdn2[k][i][j]
                              + d2gdnds[i][k]*d2sdndp[k][j]
                              + d2gdnds[j][k]*d2sdndp[k][i];
                for (l=0; l<{number_ordering}; l++) {{
                    result[o] += d3gdnds2[i][k][l]*dsdn[k][j]*dsdp[l]
                                  + d3gdnds2[j][k][l]*dsdn[k][i]*dsdp[l]
                                  + d2gds2[k][l]*d2sdn2[k][i][j]*dsdp[l]
                                  + d3gds2dp[k][l]*dsdn[k][i]*dsdn[l][j]
                                  + d2gds2[k][l]*dsdn[k][i]*d2sdndp[l][j]
                                  + d2gds2[k][l]*dsdn[k][j]*d2sdndp[l][i];
                    for (m=0; m<{number_ordering}; m++)
                        result[o] += d3gds3[k][l][m]*dsdn[k][i]*dsdn[l][j]*dsdp[m];
                }}
            }}
            o += 1;
        }}
    }}
}}

static void {module}_d3gdn3(double T, double P, double n[{number_components}],
    double result[{number_symmetric_tensor_terms}]) {{
    {moles_assign}
    double invd2gds2[{number_ordering}][{number_ordering}];
    double dsdn[{number_ordering}][{number_components}];
    double d2sdn2[{number_ordering}][{number_components}][{number_components}];
    double d2gdnds[{number_components}][{number_ordering}],
           d2gds2[{number_ordering}][{number_ordering}],
           d3gdn3[{number_components}][{number_components}][{number_components}],
           d3gdn2ds[{number_components}][{number_components}][{number_ordering}],
           d3gdnds2[{number_components}][{number_ordering}][{number_ordering}],
           d3gds3[{number_ordering}][{number_ordering}][{number_ordering}];
    int i,j,k,l,m,q,o;
    double *s = retrieveGuess(T, P, n);
    order_s(T, P, n, s, invd2gds2);
    order_dsdn(T, P, n, s, invd2gds2, dsdn);
    order_d2sdn2(T, P, n, s, invd2gds2, d2sdn2);
{order_assign}
{fillD2GDNDS}
{fillD2GDS2}
{fillD3GDN3}
{fillD3GDN2DS}
{fillD3GDNDS2}
{fillD3GDS3}
    o = 0;
    for (i=0; i<{number_components}; i++) {{
        for (j=i; j<{number_components}; j++) {{
            for (k=j; k<{number_components}; k++) {{
                result[o] = d3gdn3[i][j][k];
                for (l=0; l<{number_ordering}; l++) {{
                    result[o] += d3gdn2ds[i][j][l]*dsdn[l][k]
                               + d3gdn2ds[j][k][l]*dsdn[l][i]
                               + d3gdn2ds[k][i][l]*dsdn[l][j];
                    result[o] += d2gdnds[i][l]*d2sdn2[l][j][k]
                               + d2gdnds[j][l]*d2sdn2[l][i][k]
                               + d2gdnds[k][l]*d2sdn2[l][j][i];
                    for (m=0; m<{number_ordering}; m++) {{
                        result[o] += d3gdnds2[i][l][m]*dsdn[l][j]*dsdn[m][k]
                                   + d3gdnds2[j][l][m]*dsdn[l][k]*dsdn[m][i]
                                   + d3gdnds2[k][l][m]*dsdn[l][i]*dsdn[m][j];
                        result[o] += d2gds2[l][m]*d2sdn2[l][j][k]*dsdn[m][i]
                                   + d2gds2[l][m]*d2sdn2[l][i][k]*dsdn[m][j]
                                   + d2gds2[l][m]*d2sdn2[l][i][j]*dsdn[m][k];
                        for (q=0; q<{number_ordering}; q++)
                            result[o] += d3gds3[l][m][q]*dsdn[l][i]*dsdn[m][j]*dsdn[q][k];
                    }}
                }}
                o += 1;
            }}
        }}
    }}
}}

static void {module}_d4gdndt3(double T, double P, double n[{number_components}],
    double result[{number_components}]) {{
    int i;
    for (i=0; i<{number_components}; i++) result[i] = 0.0;
}}

static void {module}_d4gdndt2dp(double T, double P, double n[{number_components}],
    double result[{number_components}]) {{
    int i;
    for (i=0; i<{number_components}; i++) result[i] = 0.0;
}}

static void {module}_d4gdndtdp2(double T, double P, double n[{number_components}],
    double result[{number_components}]) {{
    int i;
    for (i=0; i<{number_components}; i++) result[i] = 0.0;
}}

static void {module}_d4gdndp3(double T, double P, double n[{number_components}],
    double result[{number_components}]) {{
    int i;
    for (i=0; i<{number_components}; i++) result[i] = 0.0;
}}

static void {module}_d4gdn2dt2(double T, double P, double n[{number_components}],
    double result[{number_symmetric_hessian_terms}]) {{
    int i;
    for (i=0; i<{number_symmetric_hessian_terms}; i++) result[i] = 0.0;
}}

static void {module}_d4gdn2dtdp(double T, double P, double n[{number_components}],
    double result[{number_symmetric_hessian_terms}]) {{
    int i;
    for (i=0; i<{number_symmetric_hessian_terms}; i++) result[i] = 0.0;
}}

static void {module}_d4gdn2dp2(double T, double P, double n[{number_components}],
    double result[{number_symmetric_hessian_terms}]) {{
    int i;
    for (i=0; i<{number_symmetric_hessian_terms}; i++) result[i] = 0.0;
}}

static void {module}_d5gdn2dt3(double T, double P, double n[{number_components}],
    double result[{number_symmetric_hessian_terms}]) {{
    int i;
    for (i=0; i<{number_symmetric_hessian_terms}; i++) result[i] = 0.0;
}}

static void {module}_d5gdn2dt2dp(double T, double P, double n[{number_components}],
    double result[{number_symmetric_hessian_terms}]) {{
    int i;
    for (i=0; i<{number_symmetric_hessian_terms}; i++) result[i] = 0.0;
}}

static void {module}_d5gdn2dtdp2(double T, double P, double n[{number_components}],
    double result[{number_symmetric_hessian_terms}]) {{
    int i;
    for (i=0; i<{number_symmetric_hessian_terms}; i++) result[i] = 0.0;
}}

static void {module}_d5gdn2dp3(double T, double P, double n[{number_components}],
    double result[{number_symmetric_hessian_terms}]) {{
    int i;
    for (i=0; i<{number_symmetric_hessian_terms}; i++) result[i] = 0.0;
}}

static void {module}_d4gdn3dt(double T, double P, double n[{number_components}],
    double result[{number_symmetric_tensor_terms}]) {{
    int i;
    for (i=0; i<{number_symmetric_tensor_terms}; i++) result[i] = 0.0;
}}

static void {module}_d4gdn3dp(double T, double P, double n[{number_components}],
    double result[{number_symmetric_tensor_terms}]) {{
    int i;
    for (i=0; i<{number_symmetric_tensor_terms}; i++) result[i] = 0.0;
}}

static void {module}_d5gdn3dt2(double T, double P, double n[{number_components}],
    double result[{number_symmetric_tensor_terms}]) {{
    int i;
    for (i=0; i<{number_symmetric_tensor_terms}; i++) result[i] = 0.0;
}}

static void {module}_d5gdn3dtdp(double T, double P, double n[{number_components}],
    double result[{number_symmetric_tensor_terms}]) {{
    int i;
    for (i=0; i<{number_symmetric_tensor_terms}; i++) result[i] = 0.0;
}}

static void {module}_d5gdn3dp2(double T, double P, double n[{number_components}],
    double result[{number_symmetric_tensor_terms}]) {{
    int i;
    for (i=0; i<{number_symmetric_tensor_terms}; i++) result[i] = 0.0;
}}

static void {module}_d6gdn3dt3(double T, double P, double n[{number_components}],
    double result[{number_symmetric_tensor_terms}]) {{
    int i;
    for (i=0; i<{number_symmetric_tensor_terms}; i++) result[i] = 0.0;
}}

static void {module}_d6gdn3dt2dp(double T, double P, double n[{number_components}],
    double result[{number_symmetric_tensor_terms}]) {{
    int i;
    for (i=0; i<{number_symmetric_tensor_terms}; i++) result[i] = 0.0;
}}

static void {module}_d6gdn3dtdp2(double T, double P, double n[{number_components}],
    double result[{number_symmetric_tensor_terms}]) {{
    int i;
    for (i=0; i<{number_symmetric_tensor_terms}; i++) result[i] = 0.0;
}}

static void {module}_d6gdn3dp3(double T, double P, double n[{number_components}],
    double result[{number_symmetric_tensor_terms}]) {{
    int i;
    for (i=0; i<{number_symmetric_tensor_terms}; i++) result[i] = 0.0;
}}

    \
    """

def _create_complx_soln_calib_template_c():
    """
    C language implementation of create_complx_soln_calib_template()

    Returns
    -------
    string :
        The template string.
    """

    return """\

{extra_ordering_code}

static double {module}_dparam_{func}(double T, double P, double n[{number_components}], int index) {{
{moles_assign}
    double invd2gds2[{number_ordering}][{number_ordering}];
    double *s = retrieveGuess(T, P, n);
    order_s(T, P, n, s, invd2gds2);
{order_assign}
    double result = 0.0;
    switch (index) {{
{switch_code}
    default:
        break;
    }}
        return result;
}}
\
"""

def _create_speciation_code_template_c():
    """
    C language implementation of create_speciation_code_template()

    Returns
    -------
    string :
        The template string.
    """

    return """\
#include <stdlib.h>
#include <stdio.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_multiroots.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_cblas.h>

extern int nnlsWithConstraintMatrix(double **a, int m, int n, double *b, double *x, double *rnorm,
    double *w, double *zz, int *index, int debug);

#undef MAX
#undef MIN

#define MAX(a,b)  ((a) > (b) ? (a) : (b))
#define MIN(a,b)  ((a) < (b) ? (a) : (b))

#define SQUARE(x) ((x)*(x))

#ifndef TRUE
#define TRUE 1
#endif

#ifndef FALSE
#define FALSE 0
#endif

struct system_params {{
    gsl_vector *e;
    gsl_vector *lnQ;
    gsl_matrix *R;
    gsl_matrix *Cb;
    gsl_matrix *Cs;
}};

static int system_f(const gsl_vector *y, void *params, gsl_vector *f) {{
  gsl_vector *e   = ((struct system_params *) params)->e;
  gsl_vector *lnQ = ((struct system_params *) params)->lnQ;
  gsl_matrix *R   = ((struct system_params *) params)->R;
  gsl_matrix *Cb  = ((struct system_params *) params)->Cb;
  gsl_matrix *Cs  = ((struct system_params *) params)->Cs;
  int ne = e->size;
  int ns = lnQ->size;
  int nv = y->size-1;
  double nT = exp(y->data[y->size-1]);

  gsl_vector *prod = gsl_vector_alloc(lnQ->size);

  for (int i=0; i<ns; i++) {{
    double sum = gsl_vector_get(lnQ, i);
    for (int j=0; j<nv; j++) sum += gsl_matrix_get(R, i, j)*y->data[j];
    gsl_vector_set(prod, i, exp(sum));
  }}

  for (int i=0; i<ne; i++) {{
    double sum = 0.0;
    for (int j=0; j<ne; j++) sum += gsl_matrix_get(Cb, j, i)*exp(y->data[j]);
    for (int j=0; j<ns; j++) sum += gsl_matrix_get(Cs, j, i)*prod->data[j];
    gsl_vector_set(f, i, nT*sum - e->data[i]);
  }}

  double Xsum = 0.0;
  for (int i=0; i<ne; i++) Xsum += exp(y->data[i]);
  for (int i=0; i<ns; i++) Xsum += prod->data[i];
  gsl_vector_set(f, ne, Xsum - 1.0);

  gsl_vector_free(prod);
  return GSL_SUCCESS;
}}

static int system_df(const gsl_vector *y, void *params, gsl_matrix *J) {{
  gsl_vector *e   = ((struct system_params *) params)->e;
  gsl_vector *lnQ = ((struct system_params *) params)->lnQ;
  gsl_matrix *R   = ((struct system_params *) params)->R;
  gsl_matrix *Cb  = ((struct system_params *) params)->Cb;
  gsl_matrix *Cs  = ((struct system_params *) params)->Cs;
  int ne = e->size;
  int ns = lnQ->size;
  int nv = y->size-1;
  double nT = exp(y->data[y->size-1]);

  gsl_vector *prod = gsl_vector_alloc(lnQ->size);

  for (int i=0; i<ns; i++) {{
    double sum = gsl_vector_get(lnQ, i);
    for (int j=0; j<nv; j++) sum += gsl_matrix_get(R, i, j)*y->data[j];
    gsl_vector_set(prod, i, exp(sum));
  }}

  for (int i=0; i<ne; i++) {{
    for (int j=0; j<ne; j++) {{
      double sum = gsl_matrix_get(Cb, i, j)*exp(y->data[j]);
      for (int k=0; k<ns; k++) {{
        sum +=  gsl_matrix_get(R, k, i)*prod->data[k]*gsl_matrix_get(Cs, k, j);
      }}
      gsl_matrix_set (J, i, j, nT*sum);
    }}
    double sum = 0.0;
    for (int j=0; j<ne; j++) sum += gsl_matrix_get(Cb, j, i)*exp(y->data[j]);
    for (int j=0; j<ns; j++) sum += gsl_matrix_get(Cs, j, i)*prod->data[j];
    gsl_matrix_set(J, i, ne, nT*sum);
  }}

  for (int i=0; i<ne; i++) {{
    double Xsum = exp(y->data[i]);
    for (int j=0; j<ns; j++) Xsum += gsl_matrix_get(R, j, i)*prod->data[j];
    gsl_matrix_set(J, ne, i, Xsum);
  }}
  gsl_matrix_set(J, ne, ne, 0.0);

  gsl_vector_free(prod);
  return GSL_SUCCESS;
}}

static int system_fdf(const gsl_vector *y, void *params, gsl_vector *f, gsl_matrix *J) {{
  gsl_vector *e   = ((struct system_params *) params)->e;
  gsl_vector *lnQ = ((struct system_params *) params)->lnQ;
  gsl_matrix *R   = ((struct system_params *) params)->R;
  gsl_matrix *Cb  = ((struct system_params *) params)->Cb;
  gsl_matrix *Cs  = ((struct system_params *) params)->Cs;
  int ne = e->size;
  int ns = lnQ->size;
  int nv = y->size-1;
  double nT = exp(y->data[y->size-1]);

  gsl_vector *prod = gsl_vector_alloc(lnQ->size);

  for (int i=0; i<ns; i++) {{
    double sum = gsl_vector_get(lnQ, i);
    for (int j=0; j<nv; j++) sum += gsl_matrix_get(R, i, j)*y->data[j];
    gsl_vector_set(prod, i, exp(sum));
  }}

  for (int i=0; i<ne; i++) {{
    double sum = 0.0;
    for (int j=0; j<ne; j++) sum += gsl_matrix_get(Cb, j, i)*exp(y->data[j]);
    for (int j=0; j<ns; j++) sum += gsl_matrix_get(Cs, j, i)*prod->data[j];
    gsl_vector_set(f, i, nT*sum - e->data[i]);
    gsl_matrix_set(J, i, ne, nT*sum);

    for (int j=0; j<ne; j++) {{
      double sum = gsl_matrix_get(Cb, i, j)*exp(y->data[j]);
      for (int k=0; k<ns; k++) {{
        sum +=  gsl_matrix_get(R, k, i)*prod->data[k]*gsl_matrix_get(Cs, k, j);
      }}
      gsl_matrix_set (J, i, j, nT*sum);
    }}
  }}

  double Xsum = 0.0;
  for (int i=0; i<ne; i++) Xsum += exp(y->data[i]);
  for (int i=0; i<ns; i++) Xsum += prod->data[i];
  gsl_vector_set(f, ne, Xsum - 1.0);

  for (int i=0; i<ne; i++) {{
    double Xsum = exp(y->data[i]);
    for (int j=0; j<ns; j++) Xsum += gsl_matrix_get(R, j, i)*prod->data[j];
    gsl_matrix_set(J, ne, i, Xsum);
  }}
  gsl_matrix_set(J, ne, ne, 0.0);

  gsl_vector_free(prod);
  return GSL_SUCCESS;
}}

static void print_state (size_t iter, gsl_multiroot_fdfsolver *s) {{
  printf ("iter = %3lu\\n", iter);
  for (size_t i=0; i<s->x->size; i+=8) {{
    printf ("Solution:\\n");
    for (size_t j=i*8; j<MIN(i*8+8, s->x->size); j++) printf(" %10.3e", gsl_vector_get (s->x, j));
    printf ("\\nFunction:\\n");
    for (size_t j=i*8; j<MIN(i*8+8, s->x->size); j++) printf(" %10.3e", gsl_vector_get (s->f, j));
    printf ("\\n");
  }}
}}

#define ne {number_components}
#define ns {number_non_basis}

static void compute_lnQ(double t, double p, struct system_params *params) {{
    // compute mu's and get lnQ from R: lnQ = -(mu0_s - np.matmul(R,mu0_b))/(8.3143*t)
    double mu0_b[ne];
    for (int i=0; i<ne; i++) mu0_b[i] = endmember[i].mu0(t, p);
    for (int i=0; i<ns; i++) {{
        double sum = endmember[i+ne].mu0(t, p);
        for (int j=0; j<ne; j++) {{
            sum -= gsl_matrix_get(params->R, i, j)*mu0_b[j];
        }}
        gsl_vector_set(params->lnQ, i, -sum/(8.3143*t));
    }}
}}

static void gaussj(double a[ns][ns]);

static int speciate(double t, double p, double*e, double *b, double *s,
    double invd2gds2[ns][ns], int debug) {{
    static struct system_params params;
    static int first = TRUE;
    static gsl_vector *lnX;
    static double tOld, pOld, eOld[ne], bOld[ne], sOld[ns], invd2gds2Old[ns][ns];
    if (!first) {{
        int same  = (fabs(t-tOld) < 1.0e-2);
            same &= (fabs(p-pOld) < 1.0e-2);
        for (int i=0; i<ne; i++) same &= (fabs(e[i]-eOld[i]) < 1.0e-7);
        if (same) {{
            for (int i=0; i<ne; i++) {{
                b[i] = bOld[i];
            }}
            for (int i=0; i<ns; i++) {{
                s[i] = sOld[i];
                for (int j=0; j<ns; j++) invd2gds2[i][j] = invd2gds2Old[i][j];
            }}
            return TRUE;
        }}
    }} else {{
        first = FALSE;
        params.e   = gsl_vector_alloc(ne);
        params.lnQ = gsl_vector_alloc(ns);
        params.R   = gsl_matrix_alloc(ns, ne);
        params.Cb  = gsl_matrix_alloc(ne, ne);
        params.Cs  = gsl_matrix_alloc(ns, ne);
        lnX = gsl_vector_alloc(ne+1);
        //
        // Fill Cb, Cs matrices for this case
{fill_Cb}
{fill_Cs}
        gsl_matrix *CbInv = gsl_matrix_alloc(ne, ne);
        gsl_permutation *perm  = gsl_permutation_alloc(ne);
        gsl_permutation_init(perm);
        gsl_linalg_LU_invert(params.Cb, perm, CbInv);
        gsl_permutation_free(perm);
        gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, params.Cs, CbInv, 0.0, params.R);
        gsl_matrix_free(CbInv);
    }}
    // Fill lnQ for this case using R
    // compute mu's and get lnQ from R: lnQ = -(mu0_s - np.matmul(R,mu0_b))/(8.3143*t)
    compute_lnQ(t, p, &params);
    if (debug) {{
        for (int i=0; i<ns; i++) printf("lnQ[%d] = %10.3e\\n", i,
            gsl_vector_get(params.lnQ, i));
    }}
    // Now compute an initial guess
    double **aNNLS = (double **) malloc((size_t) ns*sizeof(double *));
    double  *bNNLS = (double *)  malloc((size_t) ns*sizeof(double));
    double  *xNNLS = (double *)  malloc((size_t) (ne+ns)*sizeof(double));
    double  *wNNLS = (double *)  malloc((size_t) (ne+ns)*sizeof(double));
    double *zzNNLS = (double *)  malloc((size_t) ns*sizeof(double));
    int *indexNNLS = (int *)     malloc((size_t) (ne+ns)*sizeof(int));
    double rNorm=0.0;
    for (int i=0; i<ns; i++) {{
        bNNLS[i] = -gsl_vector_get(params.lnQ, i);
        aNNLS[i] = (double *) malloc((size_t) (ne+ns)*sizeof(double));
        for (int j=0; j<ne; j++) {{
            aNNLS[i][j]    = -gsl_matrix_get(params.R, i, j);
        }}
        for (int j=0; j<ns; j++) {{
            aNNLS[i][ne+j] = (i == j) ? 1.0 : 0.0;
        }}
    }}
    int success = nnlsWithConstraintMatrix(aNNLS, ns, ne+ns, bNNLS, xNNLS, &rNorm,
        wNNLS, zzNNLS, indexNNLS, debug);
    if (!success) printf("Speciate error. Failure in initial guess routine (NNLS).\\n");
    if (debug) {{
        printf ("NNLS solution: (rNorm = %10.3e)\\n", rNorm);
        for (int i=0; i<ne+ns; i++) printf (" %10.2e", xNNLS[i]);
        printf ("\\n");
    }}
    for (int i=0; i<ne; i++) gsl_vector_set(lnX, i, -xNNLS[i]);
    for (int i=0; i<ns; i++) free(aNNLS[i]);
    free (aNNLS); free(bNNLS); free(xNNLS); free(wNNLS); free(zzNNLS); free(indexNNLS);

    // compute mu's and get lnQ from R: lnQ = -(mu0_s - np.matmul(R,mu0_b))/(8.3143*t)
    compute_lnQ(t, p, &params);
    if (debug) {{
        for (int i=0; i<ns; i++) printf("lnQ[%d] = %10.3e\\n", i,
            gsl_vector_get(params.lnQ, i));
    }}

    // set e
    double nT = 0.0;
    for (int i=0; i<ne; i++) {{
        gsl_vector_set(params.e, i, e[i]);
        nT += e[i];
    }}
    gsl_vector_set(lnX, ne, log(nT));

    // options:
    // gsl_multiroot_fdfsolver_hybridsj
    // gsl_multiroot_fdfsolver_hybridj
    // gsl_multiroot_fdfsolver_newton
    // gsl_multiroot_fdfsolver_gnewton
    const gsl_multiroot_fdfsolver_type *Tgsl = {gsl_multiroot_method};
    gsl_multiroot_fdfsolver *soln = gsl_multiroot_fdfsolver_alloc (Tgsl, ne+1);
    int status;
    size_t iter = 0;

    gsl_multiroot_function_fdf f = {{&system_f, &system_df, &system_fdf, ne+1, &params}};
    gsl_multiroot_fdfsolver_set (soln, &f, lnX);

    if (debug) print_state (iter, soln);

    do {{
        iter++;
        status = gsl_multiroot_fdfsolver_iterate (soln);
        if (debug) print_state (iter, soln);
        if (status) break;
        gsl_vector *x  = gsl_multiroot_fdfsolver_root(soln);
        gsl_vector *f  = gsl_multiroot_fdfsolver_f(soln);
        gsl_vector *dx = gsl_multiroot_fdfsolver_dx(soln);
        if (debug) {{
            for (int i=0; i<(ne+1); i++) {{
                printf ("%3d f = %13.6g x = %13.6g dx = %13.6g\\n",
                   i, gsl_vector_get(f, i), gsl_vector_get(x, i), gsl_vector_get(dx, i));
            }}
        }}
        // status = gsl_multiroot_test_residual (f, 1.e-7);
        status = gsl_multiroot_test_delta (dx, x, 1.e-7, 1.e-10);
    }} while (status == GSL_CONTINUE && iter < 50);

    // Outputs:
    // gsl_vector * gsl_multiroot_fdfsolver_root(const gsl_multiroot_fdfsolver * soln)
    // gsl_vector * gsl_multiroot_fdfsolver_f(const gsl_multiroot_fdfsolver * soln)
    // gsl_vector * gsl_multiroot_fdfsolver_dx(const gsl_multiroot_fdfsolver * soln)
    // gsl_multiroot_test_delta(const gsl_vector * dx, const gsl_vector * x, double epsabs, double epsrel)

    if (debug) {{
        printf ("status = %s\\n", gsl_strerror (status));
        print_state(iter, soln);
    }}

    gsl_vector *lnQ = params.lnQ;
    gsl_matrix *R   = params.R;
    for (int i=0; i<ne; i++) b[i] = exp(gsl_vector_get (soln->x, i));
    for (int i=0; i<ns; i++) {{
        double sum = gsl_vector_get(lnQ, i);
        for (int j=0; j<ne; j++) sum += gsl_matrix_get(R, i, j)*gsl_vector_get (soln->x, j);
        s[i] = exp(sum);
    }}

{fill_invd2gds2}

    if (debug) {{
        double eT = 0.0;
        for (int i=0; i<ne; i++) {{
            printf ("e[%d]  = %g\\n", i, e[i]);
            eT += e[i];
        }}
        printf ("eT    = %g\\n", eT);
        double Xsum = 0.0;
        for (int i=0; i<ne; i++) {{
            printf ("Xb[%d] = %g\\n", i, exp(gsl_vector_get (soln->x, i)));
            Xsum += exp(gsl_vector_get (soln->x, i));
        }}
        gsl_vector *lnQ = params.lnQ;
        gsl_matrix *R   = params.R;
        for (int i=0; i<ns; i++) {{
          double sum = gsl_vector_get(lnQ, i);
          for (int j=0; j<ne; j++) sum += gsl_matrix_get(R, i, j)*gsl_vector_get (soln->x, j);
          printf ("Xs[%d] = %g\\n", i, exp(sum));
          Xsum += exp(sum);
        }}
        printf ("Xsum  = %g\\n", Xsum);
    }}

    tOld = t;
    pOld = p;
    for (int i=0; i<ne; i++) {{
        bOld[i] = b[i];
        eOld[i] = e[i];
    }}
    for (int i=0; i<ns; i++) {{
        sOld[i] = s[i];
        for (int j=0; j<ns; j++) invd2gds2Old[i][j] = invd2gds2[i][j];
    }}
    gsl_multiroot_fdfsolver_free (soln);
    return TRUE;
}}
\
"""

def _create_speciation_ordering_code_template_c():
    """
    C language implementation of create_speciation_ordering_code_template()

    Returns
    -------
    string :
        The template string.
    """

    return """\

static void order_dsdn(double T, double P, double n[{NC}], double b[{NC}],
    double s[{NS}], double invd2gds2[{NS}][{NS}], double dsdn[{NS}][{NC}]) {{
{ORDER_CODE_BLOCK_FIVE}
    int i,j,k;
    for (i=0; i<{NS}; i++) {{
        for (j=0; j<{NC}; j++) {{
            dsdn[i][j] = 0.0;
            for (k=0; k<{NS}; k++) dsdn[i][j] += - invd2gds2[i][k]*d2gdnds[j][k];
        }}
    }}
}}

static void order_dsdt(double T, double P, double n[{NC}], double b[{NC}],
    double s[{NS}], double invd2gds2[{NS}][{NS}], double dsdt[{NS}]) {{
    int i,j;
{ORDER_CODE_BLOCK_SIX}
    for (i=0; i<{NS}; i++) {{
        dsdt[i] = 0.0;
        for (j=0; j<{NS}; j++) dsdt[i] += - invd2gds2[i][j]*d2gdsdt[j];
    }}
}}

static void order_dsdp(double T, double P, double n[{NC}], double b[{NC}],
    double s[{NS}], double invd2gds2[{NS}][{NS}], double dsdp[{NS}]) {{
    int i,j;
{ORDER_CODE_BLOCK_SEVEN}
    for (i=0; i<{NS}; i++) {{
        dsdp[i] = 0.0;
        for (j=0; j<{NS}; j++) dsdp[i] += - invd2gds2[i][j]*d2gdsdp[j];
    }}
}}

static void order_d2sdn2(double T, double P, double n[{NC}], double b[{NC}],
    double s[{NS}], double invd2gds2[{NS}][{NS}], double d2sdn2[{NS}][{NC}][{NC}]) {{
    double dsdn[{NS}][{NC}], temp[{NS}];
    int i, j, k, l, m, o;
{ORDER_CODE_BLOCK_EIGHT}
    /* compute dsdn matrix */
    for (i=0; i<{NS}; i++) {{
        for (j=0; j<{NC}; j++) {{
        dsdn[i][j] = 0.0;
            for (k=0; k<{NS}; k++) dsdn[i][j] += - invd2gds2[i][k]*d2gdnds[j][k];
        }}
    }}

    /* compute dsdn2 cube */
    for (i=0; i<{NS}; i++) {{
        for (j=0; j<{NC}; j++) {{
            for (k=0; k<{NC}; k++) {{
                for (l=0; l<{NS}; l++) {{
                    temp[l] = d3gdn2ds[j][k][l];
                    for (m=0; m<{NS}; m++) {{
                        temp[l] += d3gdnds2[j][l][m]*dsdn[m][k]
                        + d3gdnds2[k][l][m]*dsdn[m][j];
                        for (o=0; o<{NS}; o++)
                            temp[l] += d3gds3[l][m][o]*dsdn[m][j]*dsdn[o][k];
                    }}
                }}
                d2sdn2[i][j][k] = 0.0;
                for (l=0; l<{NS}; l++) d2sdn2[i][j][k] += - invd2gds2[i][l]*temp[l];
            }}
        }}
    }}
}}

static void order_d2sdndt(double T, double P, double n[{NC}], double b[{NC}],
    double s[{NS}], double invd2gds2[{NS}][{NS}], double d2sdndt[{NS}][{NC}]) {{
    double dsdn[{NS}][{NC}], dsdt[{NS}], temp[{NS}];
    int i, j, k, l, m;
{ORDER_CODE_BLOCK_NINE}

    /* compute dsdn matrix */
    for (i=0; i<{NS}; i++) {{
        for (j=0; j<{NC}; j++) {{
            dsdn[i][j] = 0.0;
            for (k=0; k<{NS}; k++) dsdn[i][j] += - invd2gds2[i][k]*d2gdnds[j][k];
        }}
    }}

    /* compute dsdt vector */
    for (i=0; i<{NS}; i++) {{
        dsdt[i] = 0.0;
        for (j=0; j<{NS}; j++) dsdt[i] += - invd2gds2[i][j]*d2gdsdt[j];
    }}

    /* compute d2sdndt matrix */
    for (i=0; i<{NS}; i++) {{
        for (j=0; j<{NC}; j++) {{
            for (k=0; k<{NS}; k++) {{
                temp[k] = d3gdndsdt[j][k];
                for (l=0; l<{NS}; l++) {{
                    temp[k] += d3gdnds2[j][k][l]*dsdt[l] + d3gds2dt[k][l]*dsdn[l][j];
                    for (m=0; m<{NS}; m++) temp[k] += d3gds3[k][l][m]*dsdn[l][j]*dsdt[m];
                }}
            }}
            d2sdndt[i][j] = 0.0;
            for (k=0; k<{NS}; k++) d2sdndt[i][j] += - invd2gds2[i][k]*temp[k];
        }}
    }}
}}

static void order_d2sdndp(double T, double P, double n[{NC}], double b[{NC}],
    double s[{NS}], double invd2gds2[{NS}][{NS}], double d2sdndp[{NS}][{NC}]) {{
    double dsdn[{NS}][{NC}], dsdp[{NS}], temp[{NS}];
    int i, j, k, l, m;
{ORDER_CODE_BLOCK_TEN}

    /* compute dsdn matrix */
    for (i=0; i<{NS}; i++) {{
        for (j=0; j<{NC}; j++) {{
            dsdn[i][j] = 0.0;
            for (k=0; k<{NS}; k++) dsdn[i][j] += - invd2gds2[i][k]*d2gdnds[j][k];
        }}
    }}

    /* compute dsdp vector */
    for (i=0; i<{NS}; i++) {{
        dsdp[i] = 0.0;
        for (j=0; j<{NS}; j++) dsdp[i] += - invd2gds2[i][j]*d2gdsdp[j];
    }}

    /* compute d2sdndp matrix */
    for (i=0; i<{NS}; i++) {{
        for (j=0; j<{NC}; j++) {{
            for (k=0; k<{NS}; k++) {{
                temp[k] = d3gdndsdp[j][k];
                for (l=0; l<{NS}; l++) {{
                    temp[k] += d3gdnds2[j][k][l]*dsdp[l] + d3gds2dp[k][l]*dsdn[l][j];
                    for (m=0; m<{NS}; m++) temp[k] += d3gds3[k][l][m]*dsdn[l][j]*dsdp[m];
                }}
            }}
            d2sdndp[i][j] = 0.0;
            for (k=0; k<{NS}; k++) d2sdndp[i][j] += - invd2gds2[i][k]*temp[k];
        }}
    }}
}}

static void order_d2sdt2(double T, double P, double n[{NC}], double b[{NC}],
    double s[{NS}], double invd2gds2[{NS}][{NS}], double d2sdt2[{NS}]) {{
    double dsdt[{NS}], temp[{NS}];
    int i, j, k, l;
{ORDER_CODE_BLOCK_ELEVEN}

    /* compute dsdt vector */
    for (i=0; i<{NS}; i++) {{
        dsdt[i] = 0.0;
        for (j=0; j<{NS}; j++) dsdt[i] += - invd2gds2[i][j]*d2gdsdt[j];
    }}

    /* compute d2sdt2 vector */
    for (i=0; i<{NS}; i++) {{
        for (j=0; j<{NS}; j++) {{
            temp[j] = d3gdsdt2[j];
            for (k=0; k<{NS}; k++) {{
                temp[j] +=  2.0*d3gds2dt[j][k]*dsdt[k];
                for (l=0; l<{NS}; l++) temp[j] += d3gds3[j][k][l]*dsdt[k]*dsdt[l];
            }}
        }}
        d2sdt2[i] = 0.0;
        for (j=0; j<{NS}; j++) d2sdt2[i] += - invd2gds2[i][j]*temp[j];
    }}
}}

static void order_d2sdtdp(double T, double P, double n[{NC}], double b[{NC}],
    double s[{NS}], double invd2gds2[{NS}][{NS}], double d2sdtdp[{NS}]) {{
    double dsdt[{NS}], dsdp[{NS}], temp[{NS}];
    int i, j, k, l;
{ORDER_CODE_BLOCK_TWELVE}

    /* compute dsdt vector */
    for (i=0; i<{NS}; i++) {{
        dsdt[i] = 0.0;
        for (j=0; j<{NS}; j++) dsdt[i] += - invd2gds2[i][j]*d2gdsdt[j];
    }}

    /* compute dsdp vector */
    for (i=0; i<{NS}; i++) {{
        dsdp[i] = 0.0;
        for (j=0; j<{NS}; j++) dsdp[i] += - invd2gds2[i][j]*d2gdsdp[j];
    }}

    /* compute d2sdtdp vector */
    for (i=0; i<{NS}; i++) {{
        for (j=0; j<{NS}; j++) {{
            temp[j] = d3gdsdtdp[j];
            for (k=0; k<{NS}; k++) {{
                temp[j] += d3gds2dt[j][k]*dsdp[k] + d3gds2dp[j][k]*dsdt[k];
                for (l=0; l<{NS}; l++) temp[j] += d3gds3[j][k][l]*dsdt[k]*dsdp[l];
            }}
        }}
        d2sdtdp[i] = 0.0;
        for (j=0; j<{NS}; j++) d2sdtdp[i] += - invd2gds2[i][j]*temp[j];
    }}
}}

static void order_d2sdp2(double T, double P, double n[{NC}], double b[{NC}],
    double s[{NS}], double invd2gds2[{NS}][{NS}], double d2sdp2[{NS}]) {{
    double dsdp[{NS}], temp[{NS}];
    int i, j, k, l;
{ORDER_CODE_BLOCK_THIRTEEN}

    /* compute dsdp vector */
    for (i=0; i<{NS}; i++) {{
        dsdp[i] = 0.0;
        for (j=0; j<{NS}; j++) dsdp[i] += - invd2gds2[i][j]*d2gdsdp[j];
    }}

    /* compute d2sdp2 vector */
    for (i=0; i<{NS}; i++) {{
        for (j=0; j<{NS}; j++) {{
            temp[j] = d3gdsdp2[j];
            for (k=0; k<{NS}; k++) {{
                temp[j] +=  2.0*d3gds2dp[j][k]*dsdp[k];
                for (l=0; l<{NS}; l++) temp[j] += d3gds3[j][k][l]*dsdp[k]*dsdp[l];
            }}
        }}
        d2sdp2[i] = 0.0;
        for (j=0; j<{NS}; j++) d2sdp2[i] += - invd2gds2[i][j]*temp[j];
    }}
}}

\
"""

#######################
# C++ implementations #
#######################

def _create_calib_c_template_cpp():
    """
    C++ language implementation of create_calib_c_template()

    Returns
    -------
    string :
        The template string.
    """

    return ""

def _create_calib_h_template_cpp():
    """
    C++ language implementation of create_calib_h_template()

    Returns
    -------
    string :
        The template string.
    """

    return ""

def _create_calib_pyx_template_cpp():
    """
    C++ language implementation of create_calib_pyx_template()

    Returns
    -------
    string :
        The template string.
    """

    return ""

def _create_code_for_born_functions_cpp():
    """
    C++ language implementation of create_code_for_born_functions()

    Returns
    -------
    string :
        The template string.
    """

    return ""

def _create_code_for_debye_function_cpp():
    """
    C++ language implementation of create_code_for_debye_function()

    Returns
    -------
    string :
        The template string.
    """

    return ""

def _create_code_for_dh_functions_cpp():
    """
    C++ language implementation of create_code_for_dh_functions()

    Returns
    -------
    string :
        The template string.
    """

    return ""

def _create_fast_c_template_cpp():
    """
    C++ language implementation of create_fast_c_template()

    Returns
    -------
    string :
        The template string.
    """

    return ""

def _create_fast_h_template_cpp():
    """
    C++ language implementation of create_fast_h_template()

    Returns
    -------
    string :
        The template string.
    """

    return ""

def _create_fast_pyx_template_cpp():
    """
    C++ language implementation of create_fast_pyx_template()

    Returns
    -------
    string :
        The template string.
    """

    return ""

def _create_pyxbld_template_cpp():
    """
    C++ language implementation of create_pyxbld_template()

    Returns
    -------
    string :
        The template string.
    """

    return ""

def _create_redundant_calib_TV_template_cpp():
    """
    C++ language implementation of create_redundant_calib_TV_template()

    Returns
    -------
    string :
        The template string.
    """

    return ""

def _create_redundant_function_template_cpp(model_type='TP'):
    """
    C++ language implementation of create_redundant_function_template()

    Parameters
    ----------
    model_type: string
        Potential type, either Gibbs free energy ('TP') or Helmholtz free
        energy ('TV')

    Returns
    -------
    string :
        The template string.
    """

    return ""

def _create_soln_calc_template_cpp():
    """
    C++ language implementation of create_soln_calc_template()

    Returns
    -------
    string :
        The template string.
    """

    return ""

def _create_soln_calib_code_template_cpp():
    """
    C++ language implementation of create_soln_calib_code_template()

    Returns
    -------
    string :
        The template string.
    """

    return ""

def _create_soln_calib_extra_template_cpp():
    """
    C++ language implementation of create_soln_calib_extra_template()

    Returns
    -------
    string :
        The template string.
    """

    return ""

def _create_soln_calib_include_template_cpp():
    """
    C++ language implementation of create_soln_calib_include_template()

    Returns
    -------
    string :
        The template string.
    """

    return ""

def _create_soln_calib_pyx_template_cpp():
    """
    C++ language implementation of create_soln_calib_pyx_template()

    Returns
    -------
    string :
        The template string.
    """

    return ""

def _create_soln_calib_template_cpp():
    """
    C++ language implementation of create_soln_calib_template()

    Returns
    -------
    string :
        The template string.
    """

    return ""

def _create_soln_deriv_template_cpp():
    """
    C++ language implementation of create_soln_deriv_template()

    Returns
    -------
    string :
        The template string.
    """

    return ""

def _create_soln_fast_code_template_cpp():
    """
    C++ language implementation of create_soln_fast_code_template()

    Returns
    -------
    string :
        The template string.
    """

    return ""

def _create_soln_fast_include_template_cpp():
    """
    C++ language implementation of create_soln_fast_include_template()

    Returns
    -------
    string :
        The template string.
    """

    return ""

def _create_soln_fast_pyx_template_cpp():
    """
    C++ language implementation of create_soln_fast_pyx_template()

    Returns
    -------
    string :
        The template string.
    """

    return ""

def _create_soln_pyxbld_template_cpp():
    """
    C++ language implementation of create_pyxbld_template()

    Returns
    -------
    string :
        The template string.
    """

    return ""

def _create_soln_redun_template_cpp():
    """
    C++ language implementation of create_soln_deriv_template()

    Returns
    -------
    string :
        The template string.
    """

    return ""

def _create_soln_std_state_include_template_cpp():
    """
    C++ language implementation of create_soln_std_state_include_template()

    Returns
    -------
    string :
        The template string.
    """

    return ""

def _create_ordering_gaussj_template_cpp():
    """
    C++ language implementation of create_ordering_gaussj_template()

    Returns
    -------
    string :
        The template string.
    """

    return ""

def _create_ordering_code_template_cpp():
    """
    C++ language implementation of create_ordering_code_template()

    Returns
    -------
    string :
        The template string.
    """

    return ""

def _create_complx_soln_calc_template_cpp():
    """
    C++ language implementation of create_complx_soln_calc_template()

    Returns
    -------
    string :
        The template string.
    """

    return ""

def _create_complx_soln_calib_template_cpp():
    """
    C++ language implementation of create_complx_soln_calib_template()

    Returns
    -------
    string :
        The template string.
    """

    return ""

def _create_speciation_code_template_cpp():
    """
    C++ language implementation of create_speciation_code_template()

    Returns
    -------
    string :
        The template string.
    """

    return ""

def _create_speciation_ordering_code_template_cpp():
    """
    C++ language implementation of create_speciation_ordering_code_template()

    Returns
    -------
    string :
        The template string.
    """

    return ""
