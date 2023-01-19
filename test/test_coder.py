import collections
from typing import NamedTuple

import sympy as sym

from thermoengine import model, coder
from thermoengine.coder import Debye as db
from thermoengine.test_utils import dict_values_are_close


class TestCoderPhaseCreation:
    FORSTERITE_PARAMS_STIXRUDE = {
        'a0': -2055.0 * 1000.0, 'n': 7.0, 'v0': 43.6 / 10.0, 'k00': 128.0 * 10000.0,
        'k0p': 4.2, 'theta0': 809.0, 'gamma0': 0.99, 'q': 2.1, 'refS': 0.0,
        'R': 8.314472, 'T_r': 300.00, 'V_r': 43.6 / 10.0}

    def test_should_create_Forsterite_endmem_from_Stixrude(self, tmp_path):
        phase_model = self._generate_Stixrude_pure_phase_model(model_name='stixrude_tmp')
        self._setup_working_dir(tmp_path)
        result = self._create_phase_module(
            phase_model, phase_name='Forsterite', formula='Mg(2)Si(1)O(4)',
            params=self.FORSTERITE_PARAMS_STIXRUDE)

        gibbs_surf = GibbsSurface('stixrude_tmp', 'Forsterite')
        surface_props = gibbs_surf.eval_surface_and_derivs(T=1000.0, P=1e4)

        EXPECTED_SURFACE_PROPS = {
            'Cp': 174.82540409193032,
            'Cv': 169.15509024550337,
            'G': -2147814.9365130505,
            'K': 1180799.9403863563,
            'Kp': 4.483865674080948,
            'S': 274.7395550094689,
            'V': 4.414370884328149,
            'alpha': 3.298230307060583e-05,
            'beta': 8.468835115902878e-07,
            'd2GdP2': -3.7384579159817474e-06,
            'd2GdTdP': 0.0001455961183729693,
            'd3GdP3': 1.7362129128106605e-11,
            'd3GdT2dP': 4.705474115146101e-08,
            'd3GdT3': 0.00015573012556530512,
            'd3GdTdP2': -8.244761179623636e-10,
            'dCpdT': 0.01909527852662521,
            'dGdP': 4.414370884328149,
            'dGdT': -274.7395550094689}

        assert dict_values_are_close(surface_props, EXPECTED_SURFACE_PROPS)

    def _setup_working_dir(self, tmp_path):
        # Avoid makeing a subdirectory
        # Make a working sub-directory and move down into the directory.  This is done so that generated files will not clash between alternate phase_model configurations.
        # model_working_dir = 'working'
        # d = tmp_path / model_working_dir
        # d.mkdir()
        # os.chdir(d)
        # print('current dir = ', os.curdir)
        # subprocess.run(['mkdir -p', model_working_dir], shell=True)
        # subprocess.run(['cd', model_working_dir], shell=True)
        pass

    def _create_phase_module(self, phase_model, phase_name=None,
                             formula=None, params=None):
        phase_model.set_reference_origin(Vr=params['V_r'])
        result = phase_model.create_code_module(phase=phase_name, formula=formula, params=params)
        return result


    def _generate_Stixrude_pure_phase_model(self, model_name='stixrude'):
        phase_model = coder.StdStateModel(model_type='TV')
        T = phase_model.get_symbol_for_t()
        V = phase_model.get_symbol_for_v()
        Tr = phase_model.get_symbol_for_tr()
        Vr = phase_model.get_symbol_for_vr()

        # ### Define model expressions applicable over all of T,P space
        # An expression for the Gibbs free energy, $G(T,P)$ or the Helmholtz energy $A(T,V)$ is constructed.  The expression may have multiple parts.  Often the heat capacity function is postulated, then integrated to yield expressions for the entahlpy, entropy, and in combination the energy potential. Then, an equation of state (EOS) is adopted and that term is integrated in pressure or volume and added to the heat capacity integrals. This proceedure is follwed here.
        # #### (1) Helmholtz free energy
        # Declare parameters of the Stixrude standard state model:
        a0, n, v0, k00, k0p, theta0, gamma0, q, refS, R = sym.symbols('a0 n v0 k00 k0p theta0 gamma0 q refS R')
        params = [('a0', 'J/m', a0), ('n', '', n), ('v0', 'J/bar-m', v0), ('k00', 'bar', k00), ('k0p', '', k0p),
                  ('theta0', 'K', theta0), ('gamma0', '', gamma0), ('q', '', q), ('refS', 'J/K-m', refS),
                  ('R', 'J/K-m', R)]

        # Define the Debye temperature:
        c1 = sym.S(9) * k00 * v0
        c2 = k0p / sym.S(2) - sym.S(2)
        c5 = sym.S(3) * gamma0
        c7 = c5 * (-sym.S(2) + sym.S(6) * gamma0 - sym.S(3) * q)
        f = (v0 / V) ** (sym.S(2) / sym.S(3)) / sym.S(2) - sym.S(1) / sym.S(2)
        d0 = theta0 * (sym.S(1) + c7 * f * f + sym.S(2) * c5 * f) ** (sym.S(1) / sym.S(2))
        d0

        # Define the Debye Helmholtz free energy ...
        # db(x) returns the Debye *integral* with upper limit *x*
        x = d0 / T
        A_db = n * R * T * (sym.S(3) * sym.ln(sym.S(1) - sym.exp(-x)) - db(x))

        # ... and from that the quasiharmonic approximation to the Helmholtz energy ...
        A_quasi = A_db - A_db.subs(T, Tr)
        # ... and finally the Stixrude phase_model expression for the Helmholtz free energy:
        A = a0 + c1 * f * f * (sym.S(1) / sym.S(2) + c2 * f) + A_quasi
        # ... and add this expression to the phase_model
        phase_model.add_expression_to_model(A, params)

        # Check that db(x) is actually returning D_3(x)
        #
        # In general
        #
        # $$
        #     \frac{d}{dx} D_n(x) = \frac{n}{e^x - 1} - n\frac{D_n(x)}{x}
        # $$

        db(T).diff(T)

        # ## Code Print the phase_model, compile the code and link a Python module
        # Name the phase_model class

        phase_model.set_module_name(model_name)
        phase_model.set_include_debye_code(include=True)
        return phase_model


class TestCoderPhaseUsage:

    def test_should_predict_Forsterite_endmem_properties_from_Stixrude(self, tmp_path):
        gibbs_surf = GibbsSurface('stixrude', 'Forsterite')
        surface_props = gibbs_surf.eval_surface_and_derivs(T=1000.0, P=1e4)

        Gval = gibbs_surf.G(T=1000.0, P=1e4)
        assert Gval == surface_props['G']

        EXPECTED_SURFACE_PROPS = {
            'Cp': 174.82540409193032,
            'Cv': 169.15509024550337,
            'G': -2147814.9365130505,
            'K': 1180799.9403863563,
            'Kp': 4.483865674080948,
            'S': 274.7395550094689,
            'V': 4.414370884328149,
            'alpha': 3.298230307060583e-05,
            'beta': 8.468835115902878e-07,
            'd2GdP2': -3.7384579159817474e-06,
            'd2GdTdP': 0.0001455961183729693,
            'd3GdP3': 1.7362129128106605e-11,
            'd3GdT2dP': 4.705474115146101e-08,
            'd3GdT3': 0.00015573012556530512,
            'd3GdTdP2': -8.244761179623636e-10,
            'dCpdT': 0.01909527852662521,
            'dGdP': 4.414370884328149,
            'dGdT': -274.7395550094689}

        assert dict_values_are_close(surface_props, EXPECTED_SURFACE_PROPS)

    def _verify_phase_attributes(self, phase=None):
        print(phase.cy_Forsterite_stixrude_identifier())
        print(phase.cy_Forsterite_stixrude_name())
        print(phase.cy_Forsterite_stixrude_formula())
        print(phase.cy_Forsterite_stixrude_mw())
        print(phase.cy_Forsterite_stixrude_elements())

    def _eval_phase_properties(self, phase=None, t=None, p=None):
        fmt = "{0:<10.10s} {1:13.6e} {2:<10.10s}"
        print(fmt.format('G', phase.cy_Forsterite_stixrude_g(t, p), 'J/m'))
        print(fmt.format('dGdT', phase.cy_Forsterite_stixrude_dgdt(t, p), 'J/K-m'))
        print(fmt.format('dGdP', phase.cy_Forsterite_stixrude_dgdp(t, p), 'J/bar-m'))
        print(fmt.format('d2GdP2', phase.cy_Forsterite_stixrude_d2gdt2(t, p), 'J/K^2-m'))
        print(fmt.format('d2GdTdP', phase.cy_Forsterite_stixrude_d2gdtdp(t, p), 'J/K-bar-m'))
        print(fmt.format('d2GdP2', phase.cy_Forsterite_stixrude_d2gdp2(t, p), 'J/bar^2-m'))
        print(fmt.format('d3GdT3', phase.cy_Forsterite_stixrude_d3gdt3(t, p), 'J/K^3-m'))
        print(fmt.format('d3GdT2dP', phase.cy_Forsterite_stixrude_d3gdt2dp(t, p), 'J/K^2-bar-m'))
        print(fmt.format('d3GdTdP2', phase.cy_Forsterite_stixrude_d3gdtdp2(t, p), 'J/K-bar^2-m'))
        print(fmt.format('d3GdP3', phase.cy_Forsterite_stixrude_d3gdp3(t, p), 'J/bar^3-m'))
        print(fmt.format('S', phase.cy_Forsterite_stixrude_s(t, p), 'J/K-m'))
        print(fmt.format('V', phase.cy_Forsterite_stixrude_v(t, p), 'J/bar-m'))
        print(fmt.format('Cv', phase.cy_Forsterite_stixrude_cv(t, p), 'J/K-m'))
        print(fmt.format('Cp', phase.cy_Forsterite_stixrude_cp(t, p), 'J/K-m'))
        print(fmt.format('dCpdT', phase.cy_Forsterite_stixrude_dcpdt(t, p), 'J/K^2-m'))
        print(fmt.format('alpha', phase.cy_Forsterite_stixrude_alpha(t, p), '1/K'))
        print(fmt.format('beta', phase.cy_Forsterite_stixrude_beta(t, p), '1/bar'))
        print(fmt.format('K', phase.cy_Forsterite_stixrude_K(t, p), 'bar'))
        print(fmt.format('Kp', phase.cy_Forsterite_stixrude_Kp(t, p), ''))

    def _compare_phase_props_with_standard(self, phase=None, t=None, p=None):
        stixrudeDB = model.Database(database="Stixrude")

        abbrv = ""
        for full_name, abbrv in zip(stixrudeDB.phase_info.phase_name, stixrudeDB.phase_info.abbrev):
            if full_name == 'Forsterite':
                break
        refModel = stixrudeDB.get_phase(abbrv)

        import math
        fmt = "{0:<10.10s} {1:13.6e} {2:13.6e} {3:13.6e} {4:6.2f}%"
        fmts = "{0:<10.10s} {1:13.6e}"
        x = phase.cy_Forsterite_stixrude_g(t, p)
        y = refModel.gibbs_energy(t, p)
        print(fmt.format('G', x, y, x - y, 100.0 * math.fabs((x - y) / y)))
        x = phase.cy_Forsterite_stixrude_dgdt(t, p)
        y = -refModel.entropy(t, p)
        print(fmt.format('dGdT', x, y, x - y, 100.0 * math.fabs((x - y) / y)))
        x = phase.cy_Forsterite_stixrude_dgdp(t, p)
        y = refModel.volume(t, p)
        print(fmt.format('dGdP', x, y, x - y, 100.0 * math.fabs((x - y) / y)))
        x = phase.cy_Forsterite_stixrude_d2gdt2(t, p)
        print(fmts.format('d2GdT2', x))
        x = phase.cy_Forsterite_stixrude_d2gdtdp(t, p)
        print(fmts.format('d2GdTdP', x))
        x = phase.cy_Forsterite_stixrude_d2gdp2(t, p)
        print(fmts.format('d2GdP2', x))
        x = phase.cy_Forsterite_stixrude_d3gdt3(t, p)
        print(fmts.format('d3GdT3', x))
        x = phase.cy_Forsterite_stixrude_d3gdt2dp(t, p)
        print(fmts.format('d3GdT2dP', x))
        x = phase.cy_Forsterite_stixrude_d3gdtdp2(t, p)
        print(fmts.format('d3GdTdP2', x))
        x = phase.cy_Forsterite_stixrude_d3gdp3(t, p)
        print(fmts.format('d3GdP3', x))
        x = phase.cy_Forsterite_stixrude_s(t, p)
        y = refModel.entropy(t, p)
        print(fmt.format('S', x, y, x - y, 100.0 * math.fabs((x - y) / y)))
        x = phase.cy_Forsterite_stixrude_v(t, p)
        y = refModel.volume(t, p)
        print(fmt.format('V', x, y, x - y, 100.0 * math.fabs((x - y) / y)))
        x = phase.cy_Forsterite_stixrude_cv(t, p)
        print(fmts.format('Cv', x))
        x = phase.cy_Forsterite_stixrude_cp(t, p)
        y = refModel.heat_capacity(t, p)
        print(fmt.format('Cp', x, y, x - y, 100.0 * math.fabs((x - y) / y)))
        x = phase.cy_Forsterite_stixrude_dcpdt(t, p)
        print(fmts.format('dCpdT', x))
        x = phase.cy_Forsterite_stixrude_alpha(t, p)
        print(fmts.format('alpha', x))
        x = phase.cy_Forsterite_stixrude_beta(t, p)
        print(fmts.format('beta', x))
        x = phase.cy_Forsterite_stixrude_K(t, p)
        print(fmts.format('K', x))
        x = phase.cy_Forsterite_stixrude_Kp(t, p)
        print(fmts.format('Kp', x))





class GibbsSurface:
    _pkg_name: str
    _phase_name: str
    _methods: NamedTuple

    def __init__(self, pkg_name: str, phase_name: str):
        self._init_coder_module(pkg_name, phase_name)
        self._init_gibbs_methods()

    def _init_coder_module(self, pkg_name, phase_name):
        self._pkg_name = pkg_name
        self._phase_name = phase_name

        self._phase_model = coder.import_coder_phase(pkg_name)
        self._method_prefix = self._get_coder_method_prefix(pkg_name, phase_name)

    def _get_coder_method_prefix(self, pkg_name, phase_name):
        return 'cy_' + phase_name + '_' + pkg_name + '_'

    def _init_gibbs_methods(self):
        def get_method(name):
            return getattr(self._phase_model, self._method_prefix + name)

        Methods = collections.namedtuple(
            'methods',
            ['G', 'dGdT', 'dGdP', 'd2GdT2', 'd2GdTdP', 'd2GdP2',
             'd3GdT3', 'd3GdT2dP', 'd3GdTdP2', 'd3GdP3',
             'S', 'V', 'Cv', 'Cp', 'dCpdT',
             'alpha', 'beta', 'K', 'Kp'])

        self._methods = Methods(G=get_method('g'),
                                dGdT=get_method('dgdt'),
                                dGdP=get_method('dgdp'),
                                d2GdT2=get_method('d2gdt2'),
                                d2GdTdP=get_method('d2gdtdp'),
                                d2GdP2=get_method('d2gdp2'),
                                d3GdT3=get_method('d3gdt3'),
                                d3GdT2dP=get_method('d3gdt2dp'),
                                d3GdTdP2=get_method('d3gdtdp2'),
                                d3GdP3=get_method('d3gdp3'),
                                S=get_method('s'),
                                V=get_method('v'),
                                Cv=get_method('cv'),
                                Cp=get_method('cp'),
                                dCpdT=get_method('dcpdt'),
                                alpha=get_method('alpha'),
                                beta=get_method('beta'),
                                K=get_method('K'),
                                Kp=get_method('Kp'))

    def eval_surface_and_derivs(self, T, P):
        surface_props = {}
        surface_props['G'] = self.G(T, P)
        surface_props['dGdT'] = self.dGdT(T, P)
        surface_props['dGdP'] = self.dGdP(T, P)
        surface_props['d2GdP2'] = self.d2GdT2(T, P)
        surface_props['d2GdTdP'] = self.d2GdTdP(T, P)
        surface_props['d2GdP2'] = self.d2GdP2(T, P)
        surface_props['d3GdT3'] = self.d3GdT3(T, P)
        surface_props['d3GdT2dP'] = self.d3GdT2dP(T, P)
        surface_props['d3GdTdP2'] = self.d3GdTdP2(T, P)
        surface_props['d3GdP3'] = self.d3GdP3(T, P)
        surface_props['S'] = self.S(T, P)
        surface_props['V'] = self.V(T, P)
        surface_props['Cv'] = self.Cv(T, P)
        surface_props['Cp'] = self.Cp(T, P)
        surface_props['dCpdT'] = self.dCpdT(T, P)
        surface_props['alpha'] = self.alpha(T, P)
        surface_props['beta'] = self.beta(T, P)
        surface_props['K'] = self.K(T, P)
        surface_props['Kp'] = self.Kp(T, P)

        return surface_props

    def G(self, T, P):
        return self._methods.G(T, P)

    def dGdT(self, T, P):
        return self._methods.dGdT(T, P)

    def dGdP(self, T, P):
        return self._methods.dGdP(T, P)

    def d2GdT2(self, T, P):
        return self._methods.d2GdT2(T, P)

    def d2GdTdP(self, T, P):
        return self._methods.d2GdTdP(T, P)

    def d2GdP2(self, T, P):
        return self._methods.d2GdP2(T, P)

    def d3GdT3(self, T, P):
        return self._methods.d3GdT3(T, P)

    def d3GdT2dP(self, T, P):
        return self._methods.d3GdT2dP(T, P)

    def d3GdTdP2(self, T, P):
        return self._methods.d3GdTdP2(T, P)

    def d3GdP3(self, T, P):
        return self._methods.d3GdP3(T, P)

    def S(self, T, P):
        return self._methods.S(T, P)

    def V(self, T, P):
        return self._methods.V(T, P)

    def Cv(self, T, P):
        return self._methods.Cv(T, P)

    def Cp(self, T, P):
        return self._methods.Cp(T, P)

    def dCpdT(self, T, P):
        return self._methods.dCpdT(T, P)

    def alpha(self, T, P):
        return self._methods.alpha(T, P)

    def beta(self, T, P):
        return self._methods.beta(T, P)

    def K(self, T, P):
        return self._methods.K(T, P)

    def Kp(self, T, P):
        return self._methods.Kp(T, P)