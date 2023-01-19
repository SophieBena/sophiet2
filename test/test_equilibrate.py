from __future__ import annotations # Enable Python 4 type hints in Python 3

from pytest import raises
import time

import numpy as np


# API maintained for equilibrate package
from thermoengine.equilibrate import PhaseLibrary, GibbsMinimizer, BaseGibbsMinimizer, System

import thermoengine as thermo
from thermoengine import samples as smp
from thermoengine.samples import SampleMesh, Assemblage



from thermoengine.const import units
from thermoengine.core import UnorderedList

from thermoengine.test_utils import are_close, are_roughly_close, dict_values_are_close

CR_SPINEL_CATION_MIN = 0.1

DEFAULT_RESULT = True

bermanDB = thermo.model.Database(database='Berman')
stixrudeDB = thermo.model.Database(database='Stixrude')
MgO_SiO2_phase_library = PhaseLibrary(stixrudeDB.get_phases(
    ['Per',  # MgO
     'Fo', 'MgWds', 'MgRwd',  # Mg2SiO4
     'En', 'cEn', 'hpcEn', 'Maj', 'MgAki', 'MgPrv', 'MgPpv',  # MgSiO3
     'Qz', 'Coe', 'Sti', 'Seif']  # SiO2
))
all_stixrude_phases = stixrudeDB.get_all_phases()
stixrude_phs_library = PhaseLibrary(all_stixrude_phases)

BSE_comp = thermo.OxideMolComp(SiO2=38.80, Al2O3=1.79, FeO=5.75, MgO=50.83,
                                   CaO=2.58, Na2O=0.25)

class TestEquilUtils:
    def test_should_define_unordered_list(self):
        assert UnorderedList(['Ol', 'MgWds', 'MgRwd']) == ['Ol', 'MgWds', 'MgRwd']
        assert UnorderedList(['Ol', 'MgWds', 'MgRwd']) == ['Ol', 'MgRwd', 'MgWds']
        assert not UnorderedList(['Ol', 'Ol', 'MgRwd']) == ['Ol', 'MgRwd']
        assert UnorderedList(['Ol', 'Ol', 'MgRwd']) == ['Ol', 'MgRwd', 'Ol']


class TestSystemEquilibrate:

    def test_should_store_and_compare_environmental_conditions(self):
        system = thermo.equilibrate.System(T=1000.0, P=1)
        assert system == thermo.equilibrate.System(T=1000, P=1)
        assert not system == thermo.equilibrate.System(T=0, P=0)

    def test_should_store_and_compare_composition(self):
        system = thermo.equilibrate.System(
            comp=thermo.OxideMolComp(MgO=0.5, SiO2=0.5))

        assert system == thermo.equilibrate.System(
            comp=thermo.OxideMolComp(MgO=0.5, SiO2=0.5))
        assert not system == thermo.equilibrate.System(
            comp=thermo.OxideMolComp(MgO=0.1, SiO2=0.9))
        assert not system == thermo.equilibrate.System(
            comp=thermo.OxideMolComp(Al2O3=1))

    def test_should_define_phase_library(self):
        phs_lib = PhaseLibrary(
            stixrudeDB.get_phases(['Fo', 'MgRwd', 'MgWds']))

        assert phs_lib.available_phase_abbrevs == ['Fo', 'MgRwd', 'MgWds']
        assert phs_lib.available_phase_abbrevs == ['MgRwd', 'Fo', 'MgWds']
        assert not phs_lib.available_phase_abbrevs == ['MgRwd', 'Fo']

    def test_should_find_stable_mantle_olivine_polymorphs(self):
        phs_lib = PhaseLibrary(
            stixrudeDB.get_phases(['Fo', 'MgRwd', 'MgWds']))
        system = thermo.equilibrate.System(
            T=1673, comp=thermo.OxideMolComp(MgO=2, SiO2=1), phase_library=phs_lib,
            equil_method=GibbsMinimizer.GRIDDED)

        assert system.update(P=10*units.GPA).stable_phase_names == ['Fo']
        assert system.update(P=15*units.GPA).stable_phase_names == ['MgWds']
        assert system.update(P=20*units.GPA).stable_phase_names == ['MgRwd']
        assert not system.stable_phase_names == ['Fo']

    def test_should_reveal_miscibility_gap_for_olivine_polymorph(self):
        Wds_only = PhaseLibrary([stixrudeDB.get_phase('Wds')])
        MgFeSiO4 = thermo.OxideMolComp(MgO=1, FeO=1, SiO2=1)
        system = System(comp=MgFeSiO4, phase_library=Wds_only)
        T_crit = 496.1
        dT = 5
        assert not system.update(T=1500).spans_miscibility_gap
        assert not system.update(T=T_crit+dT).spans_miscibility_gap
        assert system.update(T=T_crit-dT).spans_miscibility_gap

    def test_should_locate_critical_pt_for_olivine_polymorph_at_desired_resolution(self):
        phs_lib = PhaseLibrary([stixrudeDB.get_phase('Wds')])
        MgFeSiO4 = thermo.OxideMolComp(MgO=1, FeO=1, SiO2=1)
        T_crit = 496.1

        system = System(T=T_crit-1, comp=MgFeSiO4, phase_library=phs_lib,
                        equil_method=GibbsMinimizer.FIXEDGRID)

        system_hires = System(T=T_crit - 1, comp=MgFeSiO4, phase_library=phs_lib,
                              equil_method=GibbsMinimizer.FIXEDGRID,
                              options={'grid_spacing':0.01})

        assert not system.spans_miscibility_gap
        assert system_hires.spans_miscibility_gap

    def test_should_reveal_miscibility_gap_for_ternary_feldspar(self):
        Fsp_only = PhaseLibrary([bermanDB.get_phase('Fsp')])
        Ab = thermo.ElemMolComp(K=0, Na=1, Ca=0, Al=1, Si=3, O=8)
        An = thermo.ElemMolComp(K=0, Na=0, Ca=1, Al=2, Si=2, O=8)
        Or = thermo.ElemMolComp(K=1, Na=0, Ca=0, Al=1, Si=3, O=8)
        system = System(T=600, phase_library=Fsp_only,
                        equil_method=GibbsMinimizer.FIXEDGRID)

        endmem = Ab
        stable_binary = Ab+An
        unstable_binary = An+Or
        unstable_ternary = Ab+An+Or

        # TODO(lone_samp_stable): simplify equilibrate if only a single sample
        #  is available to match sys comp. Need to just return that samp
        #  as stable
        system.update(comp=endmem)
        assert not system.update(comp=endmem).spans_miscibility_gap
        assert not system.update(comp=stable_binary).spans_miscibility_gap
        assert system.update(comp=unstable_binary).spans_miscibility_gap
        assert system.update(comp=unstable_ternary).spans_miscibility_gap

    def test_should_resolve_exsolved_samples_without_grid_aliasing_for_ternary_feldspar(self):
        Fsp_only = PhaseLibrary([bermanDB.get_phase('Fsp')])
        Ab = thermo.ElemMolComp(K=0, Na=1, Ca=0, Al=1, Si=3, O=8)
        An = thermo.ElemMolComp(K=0, Na=0, Ca=1, Al=2, Si=2, O=8)
        Or = thermo.ElemMolComp(K=1, Na=0, Ca=0, Al=1, Si=3, O=8)
        system = System(T=600, phase_library=Fsp_only, options={'grid_spacing': 1/10})

        stable_binary = Ab + An
        unstable_binary = An+Or
        unstable_ternary = Ab+An+Or

        assert system.update(comp=stable_binary).resolves_any_exsolved_samples()
        assert system.update(comp=unstable_binary).resolves_any_exsolved_samples()
        assert system.update(comp=unstable_ternary).resolves_any_exsolved_samples()

    def test_should_resolve_single_phase_without_grid_aliasing_for_ternary_garnet(self):
        Grt_only = PhaseLibrary([bermanDB.get_phase('Grt')])
        Alm = thermo.ElemMolComp(Fe=3, Mg=0, Ca=0, Al=2, Si=3, O=12)
        Grs = thermo.ElemMolComp(Fe=0, Mg=0, Ca=3, Al=2, Si=3, O=12)
        Prp = thermo.ElemMolComp(Fe=0, Mg=3, Ca=0, Al=2, Si=3, O=12)

        system = System(T=600, phase_library=Grt_only, options={'grid_spacing': 1/4})

        stable_ternary = Alm+Grs+Prp
        system.update(comp=stable_ternary)
        system.update(comp=stable_ternary)
        assert system.update(comp=stable_ternary).resolves_any_exsolved_samples()

    def test_should_distinguish_miscibility_gap_from_multiple_pure_phases(self):
        Mg2SiO4, SiO2 = thermo.OxideMolComp(MgO=2, SiO2=1), thermo.OxideMolComp(SiO2=1)
        system = System(T=500, phase_library=MgO_SiO2_phase_library)
        pure_phases = Mg2SiO4 + SiO2
        assert not system.update(comp=pure_phases).spans_miscibility_gap

    def test_should_find_olivine_stable_for_all_Fe_contents_at_1_bar(self):
        """
        data source: Stixrude & Lithgow-Bertelloni (2011), Fig. 10
        """
        phs_lib = PhaseLibrary(stixrudeDB.get_phases(
            ['Ol', 'Rwd', 'Wds']))
        Mg2SiO4 = thermo.OxideMolComp(MgO=2, SiO2=1)
        Fe2SiO4 = thermo.OxideMolComp(FeO=2, SiO2=1)
        endmems = np.array([Mg2SiO4,Fe2SiO4])

        system = thermo.equilibrate.System(T=1673, P=1, phase_library=phs_lib)

        assert system.update(comp=endmems.dot([0.99, 0.01])).stable_phase_names == ['Ol']
        assert system.update(comp=endmems.dot([0.5, 0.5])).stable_phase_names == ['Ol']
        assert system.update(comp=endmems.dot([0.01, 0.99])).stable_phase_names == ['Ol']

    def test_should_select_equilibration_method(self):
        system = thermo.equilibrate.System(equil_method=GibbsMinimizer.GRIDDED)
        assert system.equil_method == GibbsMinimizer.GRIDDED
        assert not system.equil_method == GibbsMinimizer.FIXEDGRID

    def test_should_find_PX_dependent_stable_solid_solution_phase(self):
        """
        data source: Stixrude & Lithgow-Bertelloni (2011), Fig. 10
        """
        phs_lib = PhaseLibrary(stixrudeDB.get_phases(
            ['Ol', 'Rwd', 'Wds']))
        Mg2SiO4 = thermo.OxideMolComp(MgO=2, SiO2=1)
        Fe2SiO4 = thermo.OxideMolComp(FeO=2, SiO2=1)
        endmems = np.array([Mg2SiO4,Fe2SiO4])


        system = thermo.equilibrate.System(T=1673, P=1, phase_library=phs_lib)
        # Mg-rich olivine polymorphs
        assert system.update(P=10*units.GPA, comp=endmems.dot([0.99, 0.01])
                             ).stable_phase_names == ['Ol']
        assert system.update(P=15*units.GPA, comp=endmems.dot([0.99, 0.01])
                             ).stable_phase_names == ['Wds']
        assert system.update(P=20*units.GPA, comp=endmems.dot([0.99, 0.01])
                             ).stable_phase_names == ['Rwd']

        # Mg/Fe mixed olivine polymorphs
        assert system.update(P=8*units.GPA, comp=endmems.dot([0.5, 0.5])
                             ).stable_phase_names == ['Ol']
        assert system.update(P=12*units.GPA, comp=endmems.dot([0.5, 0.5])
                             ).stable_phase_names == ['Rwd']

        # Fe-rich olivine polymorphs
        assert system.update(P=5*units.GPA, comp=endmems.dot([0.01, 0.99])
                             ).stable_phase_names == ['Ol']
        assert system.update(P=8*units.GPA, comp=endmems.dot([0.01, 0.99])
                             ).stable_phase_names == ['Rwd']

    def test_should_find_mantle_assemblages(self):
        system = System(T=1600, P=1, equil_method=GibbsMinimizer.FIXEDGRID,
                        options={'grid_spacing':1/10,
                                 'use_sample_pool':True},
                        phase_library=stixrude_phs_library)

        assert system.update(P=5 * units.GPA, comp=BSE_comp).stable_phase_names == UnorderedList(
            ['Ol', 'Opx', 'Cpx', 'Grt'])
        # assert system.update(P=10 * units.GPA).stable_phase_names == UnorderedList(
        #     ['Ol', 'hpCpx', 'Cpx', 'Grt'])
        assert system.update(P=17 * units.GPA).stable_phase_names == UnorderedList(
            ['Wds', 'Grt'])
        # assert system.update(P=24 * units.GPA).stable_phase_names == UnorderedList(
        #     ['CaPrv', 'PrvS', 'Fp', 'Grt'])
        assert system.update(P=100 * units.GPA).stable_phase_names == UnorderedList(
            ['PrvS', 'CaPrv', 'Fp', 'CfS'])

        assert system.update(P=136 * units.GPA).stable_phase_names == UnorderedList(
            ['PpvS', 'CaPrv', 'Fp', 'CfS'])

    def _test_should_rapidly_calculate_mantle_assem(self):
        # TODO(speedup): optimize new metastable calculation to allow rapid computation
        system = System(T=1600, P=1,
                        options={'grid_spacing': 1/8,
                                 'use_sample_pool': True},
                        phase_library=stixrude_phs_library,
                        equil_method=GibbsMinimizer.GRIDDED)

        t_start = time.time()
        system.update(P=5 * units.GPA, comp=BSE_comp)
        t_end = time.time()
        runtime_secs = t_end - t_start

        assert runtime_secs < 4

    def test_should_maintain_system_composition_during_fixed_grid_equilibration(self):
        system = System(T=1600, P=1, equil_method=GibbsMinimizer.FIXEDGRID,
                        options={'grid_spacing': 1 / 10, 'use_sample_pool':True},
                        phase_library=stixrude_phs_library)

        system.update(P=5 * units.GPA, comp=BSE_comp)


        assert system.stable_assemblage.total_comp == system.comp

    def test_should_maintain_system_composition_during_gridded_equilibration(self):
        system = System(T=1600, P=1, equil_method=GibbsMinimizer.GRIDDED,
                        options={'grid_spacing': 1 / 10, 'use_sample_pool':True},
                        phase_library=stixrude_phs_library)

        system.update(P=5 * units.GPA, comp=BSE_comp)

        # Ol,Grt,Opx,Cpx
        # 0.61173,0.12664,0.05706,0.20457

        assert system.stable_assemblage.total_comp == system.comp

class TestGibbsMinimizer:
    def _test_should_set_sample_library_to_near_endmember_comps_if_nonlinear_method(self):
        """
        Deactivated test for deleted NonlinearMinimizer,
        resurrect test & code from repo when ready to resume development
        """
        #TODO: TEST DEACTIVATED until development on the nonlinear
        #      minimizer resumes. Nonlinear Minimizer has been deleted to remove cruft
        #      resurrect this test if/when relevant.

        PrvS, Ol, Qz = stixrudeDB.get_phases(['PrvS', 'Ol', 'Qz'])

        phs_lib = PhaseLibrary([PrvS])
        minimizer = GibbsMinimizer(method=GibbsMinimizer.NONLINEAR,
                                   phase_library=phs_lib)
        self._assert_tests_each_endmember(minimizer.initial_sample_pool,
                                          endmem_num=3)  #

        phs_lib = PhaseLibrary([PrvS, Ol, Qz])
        minimizer = GibbsMinimizer(method=GibbsMinimizer.NONLINEAR,
                                   phase_library=phs_lib)
        endmem_total = 3 + 2 + 1  # (MgSi, FeSi, AlAl)O3 + (Mg, Fe)2SiO4 + SiO2
        assert len(minimizer.initial_sample_pool) == endmem_total

    def _assert_tests_each_endmember(self, sample_library, endmem_num=1):
        for samp, iX_endmem in zip(sample_library, np.eye(endmem_num)):
            assert are_close(samp.X, iX_endmem, abs_tol=1e-5)

    def test_should_update_sample_PT_conditions(self):
        comp = thermo.OxideMolComp(SiO2=1, FeO=0.2, MgO=1.8)
        minimizer = GibbsMinimizer(phase_library=stixrude_phs_library)
        minimizer.equilibrate(T=1600, P=5*units.GPA, comp=comp)

        assert minimizer.stable_assemblage[0].T == 1600
        assert minimizer.current_sample_pool[0].T == 1600

    def test_should_equilibrate_multiple_times_using_current_sample_pool(self):
        phs_lib = PhaseLibrary(stixrudeDB.get_phases(
            ['Per', 'Wus',  # MgO
             'Ol','Rwd','Wds',  # MgSiO3
             'Qz', 'Coe', 'Sti', 'Seif']  # SiO2
        ))
        comp = thermo.OxideMolComp(SiO2=1, FeO=0.2, MgO=1.8)
        minimizer = GibbsMinimizer(phase_library=phs_lib,
                                   method=GibbsMinimizer.GRIDDED,
                                   use_sample_pool=True)
        minimizer.equilibrate(T=1600, P=5*units.GPA, comp=comp)
        minimizer.equilibrate(T=1600, P=5 * units.GPA, comp=comp)

    def test_should_maintain_bulk_composition_during_equilibration(self):
        comp = thermo.OxideMolComp(SiO2=1, FeO=0.2, MgO=1.8)
        minimizer = GibbsMinimizer(phase_library=stixrude_phs_library)
        minimizer.equilibrate(T=1600, P=5 * units.GPA, comp=comp)

        assert minimizer.stable_assemblage.total_comp == comp

    def test___should_retrieve_stable_assemblage_after_equilibration(self):
        """Test to be deleted as it tests very little at this point"""
        minimizer = GibbsMinimizer(
            phase_library=MgO_SiO2_phase_library)

        assert (minimizer.stable_assemblage is None) or \
               (len(minimizer.stable_assemblage)==0)

        assert minimizer.equilibrate(
            T=1000,P=1,comp=thermo.OxideMolComp(SiO2=1, MgO=2)).stable_assemblage is not None
        # assert minimizer.stable_assemblage == Assemblage(minimizer.current_sample_pool)

    def test_should_equilibrate_by_chemical_exchange_on_stable_binary_solution_sample(self):
        phs = stixrudeDB.get_phase('Ol')
        sample = smp.SampleMaker.get_sample(phs, X=[0.5, 0.5], T=1673, P=1)

        samp_equil = sample.equilibrate_by_chemical_exchange(
            sample.chem_potential(X=[0.99, 0.01]))
        assert are_close(samp_equil.X, [0.99, 0.01], abs_tol=1e-4)

        samp_equil = sample.equilibrate_by_chemical_exchange(
            sample.chem_potential(X=[0.3, 0.7]))
        assert are_close(samp_equil.X, [0.3, 0.7], abs_tol=1e-4)

    def test_should_raise_missing_phase_library(self):
        with raises(GibbsMinimizer.MissingPhaseLibraryError):
            GibbsMinimizer()

    def test_should_raise_composition_not_viable_if_comp_outside_samp_library(self):
        minimizer = GibbsMinimizer(phase_library=PhaseLibrary(
            stixrudeDB.get_phases(['Per', 'Fo', 'En', 'Qz'])))

        with raises(minimizer.CompositionNotViable):
            assert minimizer.find_stable_phase(
                T=500, P=1, comp=thermo.OxideMolComp(Al2O3=1))

    def test_should_filter_phases_by_composition(self):
        phs_lib = PhaseLibrary(stixrudeDB.get_phases(
            ['Per',  # MgO
             'Fo', 'MgWds', 'MgRwd',  # Mg2SiO4
             'En', 'cEn', 'hpcEn', 'Maj', 'MgAki', 'MgPrv', 'MgPpv',  # MgSiO3
             'Qz', 'Coe', 'Sti', 'Seif']  # SiO2
        ))
        minimizer = GibbsMinimizer(method=GibbsMinimizer.FIXEDGRID,
                                   phase_library=phs_lib)
        SiO2 = thermo.OxideMolComp(SiO2=1)

        minimizer.init_sample_pool(T=300, P=1, comp=SiO2)
        assert len(minimizer.initial_sample_pool)==len(phs_lib.available_phases)
        assert len(minimizer.filter_phases_by_comp(SiO2).current_sample_pool) == 4

    def test_should_filter_allowable_samples_by_system_comp(self):
        minimizer = GibbsMinimizer(method=GibbsMinimizer.FIXEDGRID,
                                   phase_library=stixrude_phs_library)


        MgSiO3 = thermo.OxideMolComp(MgO=1, SiO2=1)
        minimizer.init_sample_pool(T=1000, P=1, comp=MgSiO3)

        assert UnorderedList(minimizer._sample_library.elems) == [
            'Mg', 'Si', 'O']

    def test_should_get_valid_minimizer_methods(self):
        methods = UnorderedList(GibbsMinimizer.get_valid_methods())
        assert methods == ['GibbsMinimizer.GRIDDED', 'GibbsMinimizer.FIXEDGRID',
                           'GibbsMinimizer.NONLINEAR']

    def test_should_raise_invalid_gibbs_minimizer_type(self):
        with raises(GibbsMinimizer.InvalidMinimizerMethod):
            GibbsMinimizer(method='gridded', phase_library=MgO_SiO2_phase_library)

    def test_should_find_single_stable_pure_phase_with_desired_composition(self):
        phs_lib = PhaseLibrary(stixrudeDB.get_phases(
            ['Per',  # MgO
             'Fo', 'MgWds', 'MgRwd',  # Mg2SiO4
             'En', 'cEn', 'hpcEn', 'Maj', 'MgAki', 'MgPrv', 'MgPpv',  # MgSiO3
             'Qz', 'Coe', 'Sti', 'Seif']  # SiO2
        ))
        minimizer = GibbsMinimizer(method=GibbsMinimizer.GRIDDED,
                                   phase_library=phs_lib)

        assert minimizer.find_stable_phase(
            T=500, P=1, comp=thermo.OxideMolComp(MgO=2, SiO2=1)).name == 'Fo'
        assert minimizer.find_stable_phase(
            T=500, P=1, comp=thermo.OxideMolComp(SiO2=1)).name == 'Qz'
        assert minimizer.find_stable_phase(
            T=500, P=1, comp=thermo.OxideMolComp(MgO=1, SiO2=1)).name == 'En'
        assert minimizer.find_stable_phase(
            T=500, P=75*units.GPA,
            comp=thermo.OxideMolComp(MgO=1, SiO2=1)).name == 'MgPrv'

    def test_should_find_pure_phase_assemblage_with_desired_composition(self):
        minimizer = GibbsMinimizer(method=GibbsMinimizer.GRIDDED,
                                   phase_library=MgO_SiO2_phase_library,
                                   )
        comp_Mg2SiO4, comp_SiO2 = thermo.OxideMolComp(MgO=2, SiO2=1), thermo.OxideMolComp(SiO2=1)

        minimizer.equilibrate(T=500, P=1, comp=comp_Mg2SiO4 + comp_SiO2)
        assert UnorderedList([samp.name for samp in minimizer.stable_assemblage]) == ['Fo', 'Qz']

        minimizer.equilibrate(T=500, P=75*units.GPA, comp=comp_Mg2SiO4)
        assert UnorderedList([samp.name for samp in minimizer.stable_assemblage]) == ['MgPrv', 'Per']

    def test_should_find_binary_miscibility_gap_using_grid_of_fixed_compositions(self):
        phs = stixrudeDB.get_phase('Wds')
        minimizer = GibbsMinimizer(method=GibbsMinimizer.FIXEDGRID,
                                   phase_library=PhaseLibrary([phs]),
                                   grid_spacing=0.02, use_sample_pool=True)
        MgFeSiO4 = thermo.OxideMolComp(MgO=1, FeO=1, SiO2=1)

        minimizer.equilibrate(T=490, P=1, comp=MgFeSiO4)
        assert len(minimizer.equilibrate(T=1500, P=1, comp=MgFeSiO4).stable_assemblage) == 1
        assert len(minimizer.equilibrate(T=500, P=1, comp=MgFeSiO4).stable_assemblage) == 1

        assert len(minimizer.equilibrate(T=490, P=1, comp=MgFeSiO4).stable_assemblage) == 2

    def test_should_find_ternary_miscibility_gap_using_grid_of_fixed_comps(self):

        Fsp = bermanDB.get_phase('Fsp')
        Ab = thermo.ElemMolComp(K=0, Na=1, Ca=0, Al=1, Si=3, O=8)
        An = thermo.ElemMolComp(K=0, Na=0, Ca=1, Al=2, Si=2, O=8)
        Or = thermo.ElemMolComp(K=1, Na=0, Ca=0, Al=1, Si=3, O=8)
        comp = (Ab + An + Or)*(1/3)
        minimizer = GibbsMinimizer(
            method=GibbsMinimizer.FIXEDGRID,
            phase_library=PhaseLibrary(phase_models=[Fsp]),
            grid_spacing=0.1,
            use_sample_pool=True
        )


        assert len(minimizer.equilibrate(T=300, P=1, comp=comp).stable_assemblage)==3
        [samp.X for samp in minimizer.stable_assemblage]




class TestNearlyStablePhases:
    def test_should_resolve_soln_phases_to_few_mol_percent(self):
        phases = stixrudeDB.get_phases(['PrvS','Fp','Qz'])
        spacing0 = 0.1
        minimizer = GibbsMinimizer(method=GibbsMinimizer.GRIDDED,
                                   phase_library=PhaseLibrary(phases),
                                   grid_spacing=spacing0)
        # init_sample_library = minimizer.sample_library.samples.copy()
        mixedMgFeSiO4 = thermo.OxideMolComp(MgO=.8*2, FeO=.2*2, SiO2=1)
        minimizer.equilibrate(T=1673, P=25*units.GPA, comp=mixedMgFeSiO4)

        stable_samp = minimizer.stable_assemblage.get_subset_for_phase('PrvS')
        samp_grid = minimizer.sample_library.get_subset_for_phase('PrvS')

        local_grid_spacing = SampleMesh.get_local_grid_spacing(
            X0=stable_samp.sample_endmem_comps,
            Xsamples=samp_grid.sample_endmem_comps)

        assert are_close(local_grid_spacing, spacing0/2)

    def test_should_return_affinities_that_agree_with_stable_assemblage_for_solutions(self):
        phs_lib = PhaseLibrary(stixrudeDB.get_phases(['Ol']))
        Mg2SiO4 = thermo.OxideMolComp(MgO=2, SiO2=1)
        Fe2SiO4 = thermo.OxideMolComp(FeO=2, SiO2=1)

        system = thermo.equilibrate.System(
            T=1673, P=1, phase_library=phs_lib,
            options={'grid_spacing': 1/10, 'use_sample_pool':True},
            equil_method=GibbsMinimizer.FIXEDGRID)

        comp = 0.51*Mg2SiO4 + 0.49*Fe2SiO4
        system.update(comp=comp)

        zero_affs = system.nearly_stable_samples.affinities < 1e-3
        num_zero_aff_phases = np.sum(zero_affs)
        assert num_zero_aff_phases == len(system.stable_assemblage)

    def test_should_return_accurate_affinities_for_polymorphs(self):
        pure_phases = PhaseLibrary(stixrudeDB.get_phases(
            ['Qz', 'Coe', 'Sti', 'Seif']  # SiO2
        ))

        minimizer = GibbsMinimizer(method=GibbsMinimizer.FIXEDGRID,
                                   phase_library=pure_phases)

        comp_SiO2 = thermo.OxideMolComp(SiO2=1)
        T = 1000

        minimizer.equilibrate(T=T, P=1, comp=comp_SiO2)
        zero_affs = minimizer.sample_library.affinities < 1e-3
        num_zero_aff_phases = np.sum(zero_affs)
        assert num_zero_aff_phases == 1

    def test_should_find_consistent_equil_for_range_of_resolutions(self):
        # TODO NEEDS SPEEDUP
        system = System(T=1600, P=1, options={'grid_spacing': 1/13},
                        phase_library=stixrude_phs_library,
                        equil_method=GibbsMinimizer.FIXEDGRID)
                        # options{'aff_thresh':1e1})

        # assert system.update(P=5 * units.GPA, comp=BSE_comp).stable_phase_names == UnorderedList(
        #     ['Ol', 'Opx', 'Cpx', 'Grt'])

        system.update(P=5 * units.GPA, comp=BSE_comp)

        assert system.nearly_stable_phase_names == UnorderedList(
            ['Ol', 'Opx', 'Cpx', 'Grt'])

    def _test_should_find_accurate_refined_equil_for_range_of_resolutions(self):
        # TODO fix metastability for refined Gridded minimization method
        system = System(T=1600, P=1, options={'grid_spacing': 1/13},
                        phase_library=stixrude_phs_library,
                        equil_method=GibbsMinimizer.GRIDDED)
                        # options{'aff_thresh':1e1})

        system.update(P=5 * units.GPA, comp=BSE_comp)
        assert system.nearly_stable_phase_names == UnorderedList(
            ['Ol', 'Opx', 'Cpx', 'Grt'])

    def test_should_output_equilibrated_system_summary(self):
        # TODO (sorted system summary for consistency): alphabetize by phase,
        #  then sort by comp, then amount...
        system = System(T=1600, P=100*units.GPA,
                        options={'grid_spacing': 1 / 10, 'use_sample_pool':True},
                        equil_method=GibbsMinimizer.FIXEDGRID,
                        phase_library=stixrude_phs_library)
        system.update(comp=BSE_comp)

        output = '\n'.join((
            'T = 1600.00 K, P = 100.0 GPa',
            'CaPerovskite     amt:   0.052778 mol',
            'Perovskite       amt:   0.730709 mol',
            ' | ----------- endmember mol frac ------------ |',
            ' |   MgPerovskite [MgSiO3]          0.9000 mol |',
            ' |   FePerovskite [FeSiO3]          0.1000 mol |',
            ' |   AlPerovskite [Al2O3]           0.0000 mol |',
            'CaFerritePhase   amt:   0.058424 mol',
            ' | ----------- endmember mol frac ------------ |',
            ' |           MgCF [MgAl2O4]         0.6549 mol |',
            ' |           FeCF [FeAl2O4]         0.1000 mol |',
            ' |           NaCF [NaAlSiO4]        0.2451 mol |',
            'Ferropericlase   amt:   0.158089 mol',
            ' | ----------- endmember mol frac ------------ |',
            ' |      Periclase [MgO]             0.8978 mol |',
            ' |       Wuestite [FeO]             0.1022 mol |',
        ))


        print()
        print(system.summary)

        assert str(system.summary) == output

    def test_should_find_metastability_and_zero_affinities_for_stable_phases(self):
        system = System(T=1600, P=100 * units.GPA,
                        options={'grid_spacing': 1 / 10, 'use_sample_pool': True},
                        equil_method=GibbsMinimizer.FIXEDGRID,
                        phase_library=stixrude_phs_library)
        system.update(comp=BSE_comp)

        minimizer = system.gibbs_minimizer

        # TODO: ensure that combined phases retain zero aff and proper
        #  metastability.
        assert np.all(minimizer.stable_assemblage.affinities==0)
        assert np.all(minimizer.stable_assemblage.metastability==True)
        # samp_lib = minimizer.sample_library
        # assert samp_

    def test_should_get_nearly_stable_phases_for_gridded_soln_phases(self):
        soln_phases = PhaseLibrary(stixrudeDB.get_phases(
            ['Ol', 'Wds', 'Rwd']  # (Mg,Fe)2SiO4
        ))

        minimizer = GibbsMinimizer(method=GibbsMinimizer.FIXEDGRID,
                                   grid_spacing=0.1,
                                   phase_library=soln_phases)

        comp_MgFeSiO4 = thermo.OxideMolComp(MgO=0.6*2, FeO=0.4*2, SiO2=1)
        T = 1673

        # TODO: create nearly stable phase names function that ignores repeats
        # Below bottom edge of binary loop near
        minimizer.equilibrate(T=T, P=10.5 * units.GPA, comp=comp_MgFeSiO4)
        assert UnorderedList(minimizer.stable_assemblage.names) == ['Ol']
        assert minimizer.nearly_stable_phase_names == ['Ol', 'Rwd', 'Wds']

        # Below eutectic inside binary loop
        minimizer.equilibrate(T=T, P=11.8 * units.GPA, comp=comp_MgFeSiO4)
        assert UnorderedList(minimizer.stable_assemblage.names) == ['Ol', 'Rwd']
        assert minimizer.nearly_stable_phase_names == ['Ol', 'Wds', 'Rwd']

        # Above eutectic
        minimizer.equilibrate(T=T, P=12.0 * units.GPA, comp=comp_MgFeSiO4)
        assert UnorderedList(minimizer.stable_assemblage.names) == ['Wds']
        assert minimizer.nearly_stable_phase_names == ['Ol', 'Wds', 'Rwd']

        # inside binary loop
        minimizer.equilibrate(T=T, P=13.0 * units.GPA, comp=comp_MgFeSiO4)
        assert UnorderedList(minimizer.stable_assemblage.names) == ['Wds', 'Rwd']
        assert minimizer.nearly_stable_phase_names == ['Ol','Wds', 'Rwd']
        #
        # # above top edge of binary loop
        # minimizer.equilibrate(T=T, P=14.0 * units.GPA, comp=comp_MgFeSiO4)
        # assert UnorderedList(minimizer.stable_assemblage.names) == ['Rwd']
        # assert UnorderedList(minimizer.nearly_stable_phases.names) == ['Wds', 'Rwd']

        # in stable phase region
        minimizer.equilibrate(T=T, P=18.0 * units.GPA, comp=comp_MgFeSiO4)
        assert UnorderedList(minimizer.stable_assemblage.names) == ['Rwd']
        assert minimizer.nearly_stable_phase_names == ['Rwd', 'Wds']

    def test_should_track_stable_and_nearly_stable_phase_changes_for_polymorphs(self):
        pure_phases = PhaseLibrary(stixrudeDB.get_phases(
            ['Fo', 'MgWds', 'MgRwd']  # Mg2SiO4
        ))
        minimizer = GibbsMinimizer(method=GibbsMinimizer.FIXEDGRID,
                                   grid_spacing=0.2,
                                   phase_library=pure_phases)

        comp_Mg2SiO4 = thermo.OxideMolComp(MgO=2, SiO2=1)
        T = 1673


        minimizer.equilibrate(T=T, P=12*units.GPA, comp=comp_Mg2SiO4)
        assert minimizer.stable_assemblage.names == ['Fo']
        assert UnorderedList(minimizer.nearly_stable_phases.names) == ['Fo', 'MgWds']

        minimizer.equilibrate(T=T, P=16*units.GPA, comp=comp_Mg2SiO4)
        assert minimizer.stable_assemblage.names == ['MgWds']
        assert UnorderedList(minimizer.nearly_stable_phases.names) == ['Fo', 'MgWds', 'MgRwd']

        minimizer.equilibrate(T=T, P=22 * units.GPA, comp=comp_Mg2SiO4)
        assert minimizer.stable_assemblage.names == ['MgRwd']
        assert UnorderedList(minimizer.nearly_stable_phases.names) == ['MgWds', 'MgRwd']

    def _test_should_return_stable_sample_even_if_no_metastable_samples_are_possible(self):
        Fsp_only = PhaseLibrary([bermanDB.get_phase('Fsp')])
        Ab = thermo.ElemMolComp(K=0, Na=1, Ca=0, Al=1, Si=3, O=8)
        system = System(T=600, phase_library=Fsp_only)

        endmem = Ab
        system.update(comp=endmem)

    def test_should_raise_composition_not_viable(self):
        pure_phases = PhaseLibrary(stixrudeDB.get_phases(
            ['Fo', 'Per']  # Mg2SiO4, MgO
        ))
        minimizer = GibbsMinimizer(method=GibbsMinimizer.FIXEDGRID,
                                   grid_spacing=0.2,
                                   phase_library=pure_phases)

        comp_SiO2 = thermo.OxideMolComp(SiO2=1)
        comp_MgSiO3 = thermo.OxideMolComp(MgO=1, SiO2=1)

        with raises(BaseGibbsMinimizer.CompositionNotViable):
            minimizer.equilibrate(T=1000, P=1, comp=comp_SiO2)

        with raises(BaseGibbsMinimizer.CompositionNotViable):
            minimizer.equilibrate(T=1000, P=1, comp=comp_MgSiO3)

    def test_should_only_get_nearly_stable_phases_with_viable_comp(self):
        pure_phases = PhaseLibrary(stixrudeDB.get_phases(
            ['MgRwd', 'MgPrv', 'Per', 'Wus']  # Mg2SiO4 = MgSiO3 + MgO
        ))
        minimizer = GibbsMinimizer(method=GibbsMinimizer.FIXEDGRID,
                                   grid_spacing=0.2,
                                   phase_library=pure_phases)

        comp_Mg2SiO4 = thermo.OxideMolComp(MgO=2, SiO2=1)
        T = 1000

        P_trans = 23.75

        minimizer.equilibrate(T=T, P=(P_trans - 1) * units.GPA, comp=comp_Mg2SiO4)
        assert 'Wus' not in minimizer.nearly_stable_phases.names

    def test_should_distinguish_polymorphs_within_nearly_stable_phases(self):
        pure_phases = PhaseLibrary(stixrudeDB.get_phases(
            ['MgRwd', 'MgPrv', 'Per', 'En']  # Mg2SiO4 = MgSiO3 + MgO
        ))
        minimizer = GibbsMinimizer(method=GibbsMinimizer.FIXEDGRID,
                                   grid_spacing=0.2,
                                   phase_library=pure_phases)

        comp_Mg2SiO4 = thermo.OxideMolComp(MgO=2, SiO2=1)
        T = 1000

        P_trans = 23.75

        minimizer.equilibrate(T=T, P=(P_trans - 1) * units.GPA, comp=comp_Mg2SiO4)
        assert 'En' not in minimizer.nearly_stable_phases.names

    def test_should_get_nearly_stable_phases_for_pure_phase_decomposition_rxn(self):
        pure_phases = PhaseLibrary(stixrudeDB.get_phases(
            ['MgRwd', 'MgPrv', 'Per']  # Mg2SiO4 = MgSiO3 + MgO
        ))
        minimizer = GibbsMinimizer(method=GibbsMinimizer.FIXEDGRID,
                                   grid_spacing=0.2,
                                   phase_library=pure_phases)

        comp_Mg2SiO4 = thermo.OxideMolComp(MgO=2, SiO2=1)
        T = 1000

        P_trans = 23.75

        minimizer.equilibrate(T=T, P=(P_trans-0.5) * units.GPA, comp=comp_Mg2SiO4)
        assert UnorderedList(minimizer.stable_assemblage.names) == ['MgRwd']
        assert UnorderedList(minimizer.nearly_stable_phases.names) == ['MgRwd', 'MgPrv', 'Per']

    def test_should_get_nearly_stable_phases_for_pure_phase_combination_rxn(self):
        pure_phases = PhaseLibrary(stixrudeDB.get_phases(
            ['MgRwd', 'MgPrv', 'Per']  # Mg2SiO4 = MgSiO3 + MgO
        ))
        minimizer = GibbsMinimizer(method=GibbsMinimizer.GRIDDED,
                                   grid_spacing=0.2,
                                   phase_library=pure_phases)

        comp_Mg2SiO4 = thermo.OxideMolComp(MgO=2, SiO2=1)
        T = 1000

        P_trans = 23.75

        minimizer.equilibrate(T=T, P=(P_trans+1)*units.GPA, comp=comp_Mg2SiO4)
        assert UnorderedList(minimizer.stable_assemblage.names) == ['MgPrv', 'Per']
        assert UnorderedList(minimizer.nearly_stable_phases.names) == ['MgRwd', 'MgPrv', 'Per']

    def test_should_find_multiple_nearly_stable_phases_for_olivine_polymorphs(self):
        phs_lib = PhaseLibrary(
            stixrudeDB.get_phases(['Fo', 'MgWds', 'MgRwd']))
        system = thermo.equilibrate.System(
            T=1673, phase_library=phs_lib,
            equil_method=GibbsMinimizer.FIXEDGRID, affinity_thresh=1e3)

        system.update(P=16 * units.GPA, comp=thermo.OxideMolComp(MgO=2, SiO2=1), )
        assert system.stable_phase_names == ['MgWds']
        assert UnorderedList(system.nearly_stable_phases.names) == ['MgWds', 'MgRwd', 'Fo']

    def test_should_retain_coexisting_samps_for_immiscible_soln_phases(self):
        soln_phases = PhaseLibrary(stixrudeDB.get_phases(
            ['PrvS']
        ))

        minimizer = GibbsMinimizer(method=GibbsMinimizer.FIXEDGRID,
                                   grid_spacing=0.1,
                                   phase_library=soln_phases)

        MgSiO3 = thermo.OxideMolComp(MgO=1, SiO2=1)
        Al2O3 = thermo.OxideMolComp(Al2O3=1)
        comp = 0.5*MgSiO3 + 0.5*Al2O3

        minimizer.equilibrate(T=2273, P=100 * units.GPA, comp=comp)

        assert UnorderedList(minimizer.stable_assemblage.names) == ['PrvS', 'PrvS']
        assert minimizer.nearly_stable_phases == minimizer.stable_assemblage

    def test_should_get_metastable_samps_for_immiscible_soln_phase(self):
        soln_phases = PhaseLibrary(stixrudeDB.get_phases(
            ['PrvS', 'Sti', 'MgCf']
        ))

        minimizer = GibbsMinimizer(method=GibbsMinimizer.FIXEDGRID,
                                   grid_spacing=0.1,
                                   phase_library=soln_phases, aff_thresh=2e3)

        MgSiO3 = thermo.OxideMolComp(MgO=1, SiO2=1)
        Al2O3 = thermo.OxideMolComp(Al2O3=1)
        comp = 0.5*MgSiO3 + 0.5*Al2O3

        minimizer.equilibrate(T=2273, P=100 * units.GPA, comp=comp)

        assert UnorderedList(minimizer.stable_assemblage.names) == ['MgCf', 'Sti']
        assert UnorderedList(minimizer.nearly_stable_phases.names) == [
            'MgCf', 'Sti'] + ['PrvS', 'PrvS']

    def test_should_get_all_metastable_samps_with_small_aff_for_immiscible_soln_phase(self):
        soln_phases = PhaseLibrary(stixrudeDB.get_phases(
            ['PrvS', 'Sti', 'MgCf']
        ))

        minimizer = GibbsMinimizer(method=GibbsMinimizer.FIXEDGRID,
                                   grid_spacing=0.1,
                                   phase_library=soln_phases, aff_thresh=2e3)

        MgSiO3 = thermo.OxideMolComp(MgO=1, SiO2=1)
        Al2O3 = thermo.OxideMolComp(Al2O3=1)
        comp = 0.5*MgSiO3 + 0.5*Al2O3

        minimizer.equilibrate(T=2273, P=100 * units.GPA, comp=comp)

        affs = minimizer.sample_library.affinities
        min_nonzero_aff = np.min(affs[affs > 0])
        best_metastable_samp_ind = np.where(np.abs(affs - min_nonzero_aff) < 1e-3)[0]
        best_metastable_samples = minimizer.sample_library.get_subset(best_metastable_samp_ind)
        assert np.all(best_metastable_samples.metastability)


    def test_should_get_asymmetric_metastable_samps_for_immiscible_soln_phase(self):
        soln_phases = PhaseLibrary(stixrudeDB.get_phases(
            ['PrvS', 'Sti', 'MgCf']
        ))

        minimizer = GibbsMinimizer(method=GibbsMinimizer.FIXEDGRID,
                                   grid_spacing=0.1,
                                   phase_library=soln_phases, aff_thresh=2e3)

        MgSiO3 = thermo.OxideMolComp(MgO=1, SiO2=1)
        Al2O3 = thermo.OxideMolComp(Al2O3=1)
        comp = 0.6*MgSiO3 + 0.4*Al2O3

        minimizer.equilibrate(T=2273, P=100 * units.GPA, comp=comp)

        assert UnorderedList(minimizer.stable_assemblage.names) == ['MgCf', 'Sti', 'PrvS']

        prv_samps = minimizer._sample_library.get_subset_for_phase('PrvS')
        assert np.sum(prv_samps.metastability) <= 3

        assert UnorderedList(minimizer.nearly_stable_phases.names) == [
            'PrvS', 'MgCf', 'Sti'] + ['PrvS', 'PrvS']

    def test_should_remove_redundant_pure_endmembers_when_soln_phase_present(self):
        soln_phases = PhaseLibrary(stixrudeDB.get_phases(
            ['MgCf', 'CfS']
        ))

        minimizer = GibbsMinimizer(method=GibbsMinimizer.FIXEDGRID,
                                   grid_spacing=0.1,
                                   phase_library=soln_phases, aff_thresh=2e3)

        MgAl2O4 = thermo.OxideMolComp(MgO=1, Al2O3=1)
        minimizer.equilibrate(T=2273, P=100 * units.GPA, comp=MgAl2O4)

        assert UnorderedList(minimizer.stable_assemblage.names) == ['CfS']

    def test_should_calc_equil_when_syscomp_only_allows_pure_endmember_for_soln_phase(self):
        soln_phases = PhaseLibrary(stixrudeDB.get_phases(
            ['PrvS', 'Sti', 'CfS']
        ))

        minimizer = GibbsMinimizer(method=GibbsMinimizer.FIXEDGRID,
                                   grid_spacing=0.1,
                                   phase_library=soln_phases, aff_thresh=2e3)

        MgSiO3 = thermo.OxideMolComp(MgO=1, SiO2=1)
        Al2O3 = thermo.OxideMolComp(Al2O3=1)
        comp = 0.6 * MgSiO3 + 0.4 * Al2O3

        minimizer.equilibrate(T=2273, P=100 * units.GPA, comp=comp)

        assert UnorderedList(minimizer.stable_assemblage.names) == ['PrvS', 'Sti', 'CfS']

    def test_should_calc_metastable_assem_when_syscomp_only_allows_pure_endmember_for_soln_phase(self):
        soln_phases = PhaseLibrary(stixrudeDB.get_phases(
            ['PpvS', 'Sti', 'CfS']
        ))

        minimizer = GibbsMinimizer(method=GibbsMinimizer.FIXEDGRID,
                                   grid_spacing=0.1,
                                   phase_library=soln_phases, aff_thresh=2e3)

        MgAl2SiO6 = thermo.OxideMolComp(MgO=1, Al2O3=1, SiO2=1)
        minimizer.equilibrate(T=2273, P=130 * units.GPA, comp=MgAl2SiO6)

        assert UnorderedList(minimizer.stable_assemblage.names) == ['PpvS']
        assert UnorderedList(minimizer.nearly_stable_phases.names) == ['PpvS', 'Sti', 'CfS']

    def test_should_calc_nearly_stable_phases_for_ternary_soln(self):
        soln_phases = PhaseLibrary(stixrudeDB.get_phases(
            ['PrvS', 'Sti', 'CfS']
        ))

        minimizer = GibbsMinimizer(method=GibbsMinimizer.FIXEDGRID,
                                   grid_spacing=0.1,
                                   phase_library=soln_phases, aff_thresh=2e3)

        MgSiO3 = thermo.OxideMolComp(MgO=1, SiO2=1)
        FeSiO3 = thermo.OxideMolComp(FeO=1, SiO2=1)
        Al2O3 = thermo.OxideMolComp(Al2O3=1)
        comp = 0.4*MgSiO3 + 0.1*FeSiO3 + 0.5*Al2O3

        minimizer.equilibrate(T=2273, P=100 * units.GPA, comp=comp)

        assert UnorderedList(minimizer.stable_assemblage.names) == ['Sti', 'CfS']
        assert UnorderedList(minimizer.nearly_stable_phase_names) == [
            'Sti', 'CfS', 'PrvS']



class PhaseAccuracy:


    def xtest_should_show_Wds_critical_pt_at_1000K(self):
        # store info on macroscopic vs microscopic W exchange energy
        pass


