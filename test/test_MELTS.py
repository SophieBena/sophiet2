from __future__ import annotations # Enable Python 4 type hints in Python 3

from typing import List

import numpy as np
import pandas as pd
from scipy import optimize as opt

import thermoengine as thermo
import thermoengine.model
import thermoengine.phases
from thermoengine.test_utils import are_close


class TestLiquidCrystallization:

    # TODO: This test became flaky suddenly and requires addressing
    #  (see task list for fix ideas)
    modelDB = thermo.model.Database(database='Berman')

    def test_should_first_precipitate_spinel_from_MORB_liquid(self):
        elem_comp = self._get_anhydrous_Cr_rich_MORB_comp()
        phs_sys = self.modelDB.get_phases(['Fsp', 'Ol', 'Cpx', 'SplS', 'Liq'])
        T_precip, state = self._find_liquidus_by_cooling(
            elem_comp, phs_sys, P=1e3, T0=1900, dT=10)
        other_crystal_names = self._get_crystal_phase_names(state, exclude='Spinel')

        assert all([state.tot_grams_phase(nm)==0 for nm in other_crystal_names])
        assert state.tot_grams_phase('Spinel') > 0

    def _get_crystal_phase_names(self, state, exclude=None):
        phase_names = list(state.phase_d.keys())
        other_crystal_names = phase_names.copy()
        other_crystal_names.remove('Liquid')
        if exclude is not None:
            other_crystal_names.remove(exclude)

        return other_crystal_names

    def test_should_produce_valid_MORB_elemental_composition(self):
        compDB = thermoengine.model.GeoCompDB()
        MORB_comp = compDB.get_syscomp('MORB', H2O=0)
        MORB_elem_comp = MORB_comp.mol_comp(components='elems').iloc[0]
        MORB_elem_comp /= MORB_elem_comp.sum()

        assert are_close(MORB_elem_comp, self._get_anhydrous_Cr_rich_MORB_comp(), abs_tol=1e-6)
        assert are_close(MORB_elem_comp, self._get_anhydrous_MORB_comp_legacy(), abs_tol=1e-2)
        assert not are_close(MORB_elem_comp, self._get_anhydrous_MORB_comp_legacy(), abs_tol=1e-3)

    def test_should_crystallize_Cr_spinel_at_MORB_liquidus_temperature(self):
        # This is for anhydrous MORB
        T_MORB_LIQUIDUS = 1500.0 #TODO replace with actually valid number
        CR_PARTITIONING_MIN = 100
        CR_SPINEL_CATION_MIN = 0.1

        elem_comp = self._get_anhydrous_Cr_rich_MORB_comp()
        phs_sys = self.modelDB.get_phases(['SplS', 'Liq'])

        T_precip, state = self._find_liquidus_by_cooling(elem_comp, phs_sys, P=1e3, T0=1900)
        Cr_frac_spl = self._get_spinel_Cr_abundance(state)

        CATION_FRAC = 3/7

        assert Cr_frac_spl / CATION_FRAC > CR_SPINEL_CATION_MIN
        assert Cr_frac_spl / elem_comp.Cr > CR_PARTITIONING_MIN
        # assert are_close(T_precip, T_MORB_LIQUIDUS, abs_tol=100)

    def _get_spinel_Cr_abundance(self, state):
        spinel = state.phase_d['Spinel']
        mol = spinel['moles']
        mol /= mol.sum()
        CHROMITE_IND = 0
        FORMULA_ATOM_NUM = 7
        Cr_frac_spl = 2 * mol[CHROMITE_IND] / FORMULA_ATOM_NUM
        return Cr_frac_spl

    def _find_liquidus_by_cooling(self, elem_comp, phs_sys: List[thermoengine.phases],
                                  P=1, T0=1600.0, dT=1, TOL=1e-3):
        equil = thermo.equilibrate.Equilibrate(elem_comp.index, phs_sys)
        state = equil.execute(T0, P, bulk_comp=elem_comp.values,
                              stats=True)
        T = T0
        while (state.affinities('Spinel') > TOL):
            T -= dT
            state = equil.execute(T, P, bulk_comp=elem_comp.values, stats=True)
        T_precip = T
        return T_precip, state

    def test_should_find_spinel_as_only_stable_crystal_phase_at_MORB_liquidus(self):
        pass

    def xtest_should_first_precipitate_Cr_spinel_from_MORB_liquid_legacy(self):
        T0 = 1600
        dT = 50
        P = 1000.0

        elem_comp = self._get_anhydrous_Cr_rich_MORB_comp()

        phasenames = ['Fsp', 'Ol', 'Cpx', 'SplS', 'Liq']
        phs_sys = self._get_phases('Berman', phasenames)

        equil = equilibrate.Equilibrate(elem_comp.index, phs_sys)

        def spinel_precip_fun(T):
            return self._calc_spinel_precip_progress(T, P, elem_comp, equil)

        root_results = opt.root_scalar(spinel_precip_fun, x0=T0, x1=T0-dT)
        T_precip = root_results.root

        state = equil.execute(T_precip, P, bulk_comp=elem_comp.values,
                              stats=True)
        phase_keys = state.phase_d.keys()

        # other_phases = [iphs_key if iphs_key in ['Liquid', 'Spinel']
        #                 else None for iphs_key in phase_keys]
        other_phases = [iphs_key for iphs_key in phase_keys
                        if iphs_key not in ['Liquid', 'Spinel']]

        other_masses = np.array([state.tot_grams_phase(iphs)
                                 for iphs in other_phases])
        other_affinities = np.array([state.affinities(iphs)
                                     for iphs in other_phases])

        assert are_close(state.affinities('Spinel'), 0)
        assert are_close(other_masses, 0)
        assert np.all(other_affinities > 0)

    def _calc_spinel_precip_progress(self, T, P, elem_comp, equil,
                                     mass_wt=3e4, stats=False):
        state = equil.execute(T, P, bulk_comp=elem_comp.values, stats=stats)
        A = state.affinities('Spinel')
        m = state.tot_grams_phase('Spinel')
        return mass_wt*m - A

    def _get_anhydrous_MORB_comp_legacy(self):
        elem_comp = pd.Series(
            data=[0.022203, 2.792676, 0.085513, 0.225782,
                  0.346014, 0.810195, 0.001127, 0.000637,
                  0.222007, 0.012641, 0.000559, 0.0,
                  0.116788, 0.0, 0.0],
            index=['H', 'O', 'Na', 'Mg', 'Al', 'Si', 'P',
                   'K', 'Ca', 'Ti', 'Cr', 'Mn', 'Fe',
                   'Co', 'Ni'])
        return elem_comp/elem_comp.sum()

    def _get_anhydrous_Cr_rich_MORB_comp(self):
        # NOTE: Cr is 10x more abundant than in legacy composition
        elem_comp = pd.Series(
            data=[0.00000000e+00, 6.05391526e-01, 1.91441309e-02, 4.23814692e-02,
                  7.74635199e-02, 1.76612430e-01, 2.52351146e-04, 1.42589847e-04,
                  4.96616279e-02, 2.82999095e-03, 1.25200848e-03, 0.00000000e+00,
                  2.48683560e-02, 0.00000000e+00, 0.00000000e+00],
            index=['H', 'O', 'Na', 'Mg', 'Al', 'Si', 'P', 'K', 'Ca', 'Ti',
                   'Cr', 'Mn', 'Fe', 'Co', 'Ni'])
        return elem_comp/elem_comp.sum()