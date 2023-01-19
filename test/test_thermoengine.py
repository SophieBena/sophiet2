from __future__ import annotations  # Enable Python 4 type hints in Python 3

from pytest import raises

import numpy as np
import pandas as pd
from thermoengine import model
import thermoengine as thermo
from thermoengine.test_utils import are_close, are_roughly_close

DEFAULT_RESULT = True


class TestDatabase:
    def test_should_get_phase_model_from_database(self):
        modelDB = thermo.model.Database(database='Stixrude')
        assert modelDB.get_phase('Ol').abbrev == 'Ol'
        assert not modelDB.get_phase('Ol').abbrev =='Cpx'

    def test_should_raise_invalid_phase_symbol_exception(self):
        modelDB = thermo.model.Database(database='Stixrude')
        InvalidPhaseSymbol = thermo.model.Database.InvalidPhaseSymbol
        with raises(InvalidPhaseSymbol):
            modelDB.get_phase('NotASymbol')

    def test_should_recognize_special_phase_symbols_as_valid(self):
        modelDB = thermo.model.Database()
        assert 'Liq' in modelDB.phases
        assert 'H2O' in modelDB.phases
        assert 'MtlL' in modelDB.phases
        assert 'MtlS' in modelDB.phases

    def test_should_delay_init_of_phase_until_needed(self):
        modelDB = thermo.model.Database(database='Berman')
        assert 'Ol' not in modelDB.phases

        modelDB.get_phase('Ol')
        assert 'Ol' in modelDB.phases


class QuartzCpMeas:
    QZ_CP_MEAS = pd.DataFrame(
        [['low', 401.828, 53.754],
         ['pre-lambda', 700.569, 69.605],
         ['lambda', 843.839, 90.730],
         ['post-lambda', 899.940, 68.266]],
        columns=['meas', 'T', 'Cp']).set_index('meas')

    def get_normal_measurements(self):
        exp_data = self.QZ_CP_MEAS.drop(index='lambda')
        return {'T': exp_data['T'],
                'Cp': exp_data['Cp'],
                'regime': exp_data.index}

    def get_lambda_transition_measurement(self):
        lambda_trans_data = self.QZ_CP_MEAS
        return {'T': lambda_trans_data.loc['lambda', 'T'],
                'Cp': lambda_trans_data.loc['lambda', 'Cp'],
                'regime': 'lambda'}


class TestQuartzCpMeas(QuartzCpMeas):

    def test_should_retrieve_normal_measurements(self):
        assert 'lambda' not in self.get_normal_measurements()['regime']

    def test_should_retrieve_lambda_transition_measurements(self):
        assert self.get_lambda_transition_measurement()['regime'] is 'lambda'


class TestPurePhases:
    def test_should_predict_quartz_Cp(self):
        assert are_close(*self._compare_quartz_Cp_model_and_data(database='Berman'), rel_tol=0.01)
        assert are_close(*self._compare_quartz_Cp_model_and_data(database='HollandAndPowell'), rel_tol=0.01)
        assert are_roughly_close(*self._compare_quartz_Cp_model_and_data(database='Stixrude'),
                                 rel_tol_bounds=[0.01, 0.2])

    def test_should_predict_quartz_lambda_transition_temperature(self):
        T_measured = QuartzCpMeas().get_lambda_transition_measurement()['T']
        lambda_trans = self.LambdaTransition(T_measured, 'Qz')

        assert are_close(lambda_trans.find_transition(database='HollandAndPowell'), T_measured, abs_tol=10)
        assert are_close(lambda_trans.find_transition(database='Stixrude'), T_measured, abs_tol=10)
        assert are_close(lambda_trans.find_transition(database='Berman'), T_measured, abs_tol=10)

    class LambdaTransition():
        def __init__(self, T_lambda_measured, phase_symbol, P=1.0):
            self.T_lambda_measured = T_lambda_measured
            self.phase_symbol = phase_symbol
            self.P = P

        def _get_phase(self, database, phasename):
            modelDB = model.Database(database=database)
            phs = modelDB.get_phase(phasename)
            return phs

        def find_transition(self, database='Berman', dT_bracket=100):
            phase = self._get_phase(database, self.phase_symbol)
            T_lambda_modeled = self._find_max_Cp_Temperature(phase, self.T_lambda_measured, dT_bracket)
            return T_lambda_modeled

        def _find_max_Cp_Temperature(self, Qz, T_lambda_meas, dT_bracket):
            T = T_lambda_meas + np.linspace(-dT_bracket, +dT_bracket, 101)
            Cp_mod = Qz.heat_capacity(T, self.P)
            ind_max = np.argmax(Cp_mod)
            T_lambda_mod = T[ind_max]
            return T_lambda_mod

    def _compare_quartz_Cp_model_and_data(self, database='Berman'):
        Qz = self._get_phase(database, 'Qz')
        Cp_meas = QuartzCpMeas().get_normal_measurements()
        Cp_modeled = Qz.heat_capacity(Cp_meas['T'], 1.0)
        return Cp_modeled, Cp_meas['Cp']

    def _get_phase(self, database, phasename):
        modelDB = model.Database(database=database)
        phs = modelDB.get_phase(phasename)
        return phs


class TestPhaseTransitions:
    def xtest_should_predict_mantle_phase_transitions_MgSiO3_Stixrude(self):
        assert DEFAULT_RESULT

    def xtest_should_predict_mantle_phase_transitions_SiO2_Stixrude(self):
        assert DEFAULT_RESULT

    def xtest_should_predict_phase_diagram_Al2SiO5_Berman(self):
        assert DEFAULT_RESULT

    def xtest_should_predict_clapeyron_slope_for_Brg_Ppv_transition(self):
        assert DEFAULT_RESULT


class TestMelting:
    def xtest_should_produce_basalt_by_partially_melting_peridotite(self):
        assert DEFAULT_RESULT

    def xtest_should_recover_binary_phase_loop_for_olivine(self):
        assert DEFAULT_RESULT


class TestUtils:
    def test_roughly_close_bounds(self):
        assert are_roughly_close(10, 13, rel_tol_bounds=[0.01, 0.5])