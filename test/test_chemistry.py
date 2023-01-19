from __future__ import annotations # Enable Python 4 type hints in Python 3

import numpy as np
import pandas as pd
from pytest import raises

import thermoengine as thermo
from thermoengine import chemistry, chem_library
from thermoengine.core import UnorderedList
from thermoengine.chemistry import OxideMolComp, OxideWtComp, ElemMolComp, Oxides, OxideWt

from thermoengine.test_utils import are_close, are_roughly_close, dict_values_are_close

from typing import List, Set, Dict, Tuple, Optional

DEFAULT_RESULT = True

class TestCrystalSites:
    # Crystal Sites
    def test_should_store_potential_occupants_of_crystal_site(self):
        site = chemistry._CrystalSite(['Mg', 'Fe2+', 'Al'])
        assert site.potential_occupants == ['Mg', 'Fe2+', 'Al']
        assert site.potential_occupants == ['Fe2+', 'Mg', 'Al']
        assert not site.potential_occupants == ['Fe2+', 'Al']

    def test_should_store_current_occupant_of_crystal_site(self):
        site_occupants = ['Mg', 'Fe2+', 'Al']
        default_first_element = site_occupants[0]

        site = chemistry._CrystalSite(site_occupants)
        assert site.occupancy == default_first_element

        site.occupancy = 'Fe2+'
        assert site.occupancy == 'Fe2+'

    def test_should_store_multiple_occupants_on_crystal_site(self):
        site_occupants = ['Ti4+', 'Fe2+']
        site = chemistry._CrystalSite(site_occupants)

        assert site.potential_occupants == ['Ti4+', 'Fe2+']

    def test_should_get_elemental_composition_from_site(self):
        site = chemistry._CrystalSite(['Mg'])
        assert site.composition == {'Mg':1}

        site = chemistry._CrystalSite(['Fe2+'])
        assert site.composition == {'Fe2+': 1}

    def test_should_get_element_composition_from_site_with_multiple_elements(self):
        site = chemistry._CrystalSite(['Fe2+', 'Mg'])
        site.occupancy = 'Fe2+'
        assert site.composition == {'Fe2+': 1, 'Mg': 0}

        site.occupancy = 'Mg'
        assert site.composition == {'Fe2+': 0, 'Mg': 1}

    def test_should_get_elemental_compositions_from_site_with_partial_occupancy(self):
        site = chemistry._CrystalSite(['Al', 'Mg'])
        site.occupancy = {'Mg': 0.5, 'Al':0.5}

        assert site.composition == {'Al': 0.5, 'Mg': 0.5}

    def test_should_verify_occupancy_functionality_when_given_string_and_dict(self):
        site = chemistry._CrystalSite(['Fe2+', 'Mg'])

        site.occupancy = 'Fe2+'
        composition_using_str = site.composition

        site.occupancy = {'Fe2+': 1, 'Mg': 0}
        composition_using_dict = site.composition

        assert composition_using_dict == composition_using_str

    def test_should_notify_if_attempting_to_set_unknown_site_occupant(self):
        site = chemistry._CrystalSite(['Mg', 'Fe2+', 'Al'])
        invalid_site = 'purple cow'

        UnknownSiteOccupantError = chemistry._Occupancy.UnknownSiteOccupantError
        with raises(UnknownSiteOccupantError):
            site.occupancy = invalid_site

    # Site Element Groups
    def test_should_create_site_element_group(self):
        elem_group = chemistry._SiteElementGroup(['Mg','Fe2+','Al'])
        assert elem_group == ['Mg', 'Fe2+', 'Al']
        assert elem_group == ['Fe2+', 'Al', 'Mg']
        assert not elem_group == ['Fe2+']

    def test_should_get_site_element_properties(self):
        elem_group = chemistry._SiteElementGroup(['Mg', 'Fe2+', 'Fe3+', 'Al'])

        assert elem_group.element_types == ['Mg', 'Fe', 'Fe', 'Al']
        assert not elem_group.charges == [0, 0, 0, 0]
        assert elem_group.charges == [2, 2, 3, 3]

    # Site Elements
    def test_should_get_site_element_given_element_symbol(self):
        InvalidElementSymbol = chemistry._SiteElement.InvalidElementSymbol
        assert chemistry._SiteElement.get_element('Mg').name == 'Mg'
        assert chemistry._SiteElement.get_element('Fe2+').name == 'Fe2+'

        with raises(InvalidElementSymbol):
            chemistry._SiteElement.get_element('purple cow')

    def test_should_store_element_type_for_site_elements(self):
        Fe3 = chemistry._SiteElement.get_element('Fe3+')
        Fe2 = chemistry._SiteElement.get_element('Fe2+')
        Mg = chemistry._SiteElement.get_element('Mg')
        assert Fe3.element_type == Fe2.element_type
        assert not Fe3.element_type == Mg.element_type

    def test_should_get_realistic_elements_given_symbol(self):
        Mg = chemistry._SiteElement.get_element('Mg')
        K = chemistry._SiteElement.get_element('K')
        assert not Mg.charge == K.charge

        Fe2 = chemistry._SiteElement.get_element('Fe2+')
        Fe3 = chemistry._SiteElement.get_element('Fe3+')
        assert Fe3.charge == 3
        assert Fe2.charge == 2

    def test_should_store_charge_of_site_elements(self):
        site_occupant = chemistry._SiteElement('Mg', +2, 'Mg')
        assert site_occupant.name == 'Mg'
        assert site_occupant.element_type == 'Mg'
        assert site_occupant.charge == +2

    # Site Element Library
    def test_should_create_site_element_library(self):
        ELEMENT_LIBRARY = chem_library._SiteElementLibrary.get_library()
        assert ELEMENT_LIBRARY['O']['name'] == 'O'
        assert ELEMENT_LIBRARY['O']['charge'] == -2
        assert ELEMENT_LIBRARY['Fe2+']['name'] == 'Fe2+'
        assert ELEMENT_LIBRARY['Fe3+']['element_type'] == 'Fe'
        assert ELEMENT_LIBRARY['Fe2+']['element_type'] == 'Fe'

class TestOrderedCrystals:
    def test_should_store_group_of_crystal_sites_for_ordered_crystal(self):
        crystal = chemistry._OrderedCrystal()
        crystal.add_site('SiteA', ['Mg', 'Fe2+', 'Al'])
        crystal.add_site('SiteB', ['Ca', 'Fe2+'])
        assert crystal.get_potential_site_occupants('SiteA') == ['Mg', 'Fe2+', 'Al']
        assert crystal.get_potential_site_occupants('SiteB') == ['Ca', 'Fe2+']

    def test_should_notify_if_getting_empty_crystal_elemental_abundances(self):
        empty_crystal = chemistry._OrderedCrystal()
        with raises(chemistry._OrderedCrystal.EmptyCrystalError):
            empty_crystal.get_elemental_abundances()


    def test_should_calculate_elemental_abundances_for_ordered_crystal(self):
        crystal = chemistry._OrderedCrystal()
        crystal.add_site('SiteA', ['Mg', 'Fe2+', 'Al'])
        crystal.add_site('SiteB', ['Ca', 'Fe2+'])
        composition = {'Mg': 1, 'Al': 0, 'Ca': 1, 'Fe': 0}
        assert crystal.get_elemental_abundances() == composition

    def test_should_represent_ordered_crystal(self):
        # spinel_inv = chemistry._OrderedCrystal()
        spinel_norm = chemistry._OrderedCrystal()

        spinel_norm.add_site('Tet', ['Mg', 'Al'])
        spinel_norm.add_site('Oct', ['Mg', 'Al'], multiplicity=2, occupancy='Al')
        spinel_norm.add_site('O', ['O'], multiplicity=4)

        elemental_abundances = spinel_norm.get_elemental_abundances()
        assert elemental_abundances == {'Al': 2, 'Mg': 1, 'O': 4}

    def test_should_represent_ordered_crystal_with_multivalent_atoms(self):
        magnetite = chemistry._OrderedCrystal()
        # NOTE: magnetite is inverse spinel not normal

        magnetite.add_site('Tet', ['Mg', 'Al', 'Fe2+', 'Fe3+'])
        magnetite.add_site('Oct', ['Mg', 'Al', 'Fe2+', 'Fe3+'], multiplicity=2)
        magnetite.add_site('O', ['O'], multiplicity=4)

        magnetite.set_site_occupants({'Oct': 'Fe2+'})
        magnetite.set_site_occupants({'Tet': 'Fe3+'})
        elemental_abundances = magnetite.get_elemental_abundances()
        assert elemental_abundances == {'Al': 0, 'Fe':3, 'Mg': 0, 'O': 4}

        import copy
        magnetite_inverse = copy.deepcopy(magnetite)

        # assert magnetite_inverse==magnetite
        # TODO: create simple method for copying an ordered crystal


    def test_should_notify_if_missing_crystal_site_is_requested(self):
        crystal = chemistry._OrderedCrystal()

        MissingSiteError = chemistry._OrderedCrystal.MissingSiteError
        missing_site = 'SiteC'
        with raises(MissingSiteError):
            crystal._get_site(missing_site)

    def _test_should_recognize_crystal_sites_in_site_formula(self):
        import thermo.chemistry as chemistry
        site_formula = 'Tet_1[Mg,Fe]Oct_2[Al,Cr,Ti,Fe]O_4[O]'
        # sites = chemistry._parse_site_formula(site_formula)
        # assert sites.number == 3

class TestComp:
    def test_should_define_simple_elemental_compositions_by_chemical_formula(self):
        ox = thermo.ElemMolComp.get_by_formula('MgO')
        assert ox.comp == {'Mg':1,'O':1}

        assert thermo.ElemMolComp.get_by_formula('K2O') == thermo.ElemMolComp(K=2, O=1)
        assert thermo.ElemMolComp.get_by_formula('Al2O3') == {'Al': 2, 'O':3}
        assert thermo.ElemMolComp.get_by_formula('Mg2SiO4') == {'Mg': 2, 'Si':1, 'O': 4}

    def test_should_compare_approximately_equal_compositions(self):
        assert thermo.ElemMolComp(Si=1, O=2) == thermo.ElemMolComp(Si=1.0000001, O=2)
        assert not thermo.ElemMolComp(Si=1, O=2) == thermo.ElemMolComp(Si=1.01, O=2)

    def test__should_compare_comp_objects_enabling_ordering(self):
        H = thermo.ElemMolComp(H=1)
        Fe = thermo.ElemMolComp(Fe=1)

        assert H.sort_index < Fe.sort_index
        assert H < Fe
        assert Oxides.SiO2 < Oxides.Fe2O3


    def test_should_retrieve_oxides(self):
        assert Oxides.SiO2 == {'Si':1, 'O':2}
        assert Oxides.Cr2O3 == {'Cr':2, 'O':3}

    def test_should_add_elemental_compositions(self):
        assert (
                thermo.ElemMolComp.get_by_formula('MgO') +
                thermo.ElemMolComp.get_by_formula('SiO2') ==
                thermo.ElemMolComp.get_by_formula('MgSiO3')
        )

    def test_should_add_oxides(self):
        assert (Oxides.MgO + Oxides.SiO2 ==
                thermo.ElemMolComp.get_by_formula('MgSiO3'))
        assert (2 * Oxides.MgO + Oxides.SiO2 ==
                thermo.ElemMolComp.get_by_formula('Mg2SiO4'))

    def test_should_add_oxide_comps(self):
        MgO = thermo.OxideMolComp(MgO=1)
        SiO2 = thermo.OxideMolComp(SiO2=1)
        MgSiO3 = thermo.OxideMolComp(MgO=1, SiO2=1)

        assert MgO + SiO2 == MgSiO3

    def test_should_compare_oxide_comps_with_no_shared_components(self):
        MgO, SiO2 = thermo.OxideMolComp(MgO=1), thermo.OxideMolComp(SiO2=1)
        assert MgO == thermo.OxideMolComp(MgO=1)
        assert not MgO == SiO2

    def test_should_calc_linear_combination_of_oxide_comps(self):
        MgO, SiO2, Mg2SiO4 = thermo.OxideMolComp(MgO=1), thermo.OxideMolComp(SiO2=1), thermo.OxideMolComp(MgO=2, SiO2=1)
        Mg2SiO4_elems = thermo.ElemMolComp.get_by_formula('Mg2SiO4')
        components = np.array([MgO, SiO2])

        assert 2*MgO + SiO2 == Mg2SiO4
        assert MgO*2 + SiO2 == Mg2SiO4
        assert components.dot([2, 1]) == Mg2SiO4
        assert components.dot([2, 1]) == Mg2SiO4_elems

    def test_should_represent_empty_oxide_comp_as_default(self):
        empty_comp = thermo.OxideMolComp()

        assert empty_comp.data_is_empty
        assert empty_comp == {}

    def test__should_get_non_default_data_from_oxide_comp(self):
        comp = thermo.OxideMolComp(MgO=1, SiO2=1)
        assert {'MgO':1, 'SiO2':1} == comp.data

    def test__should_get_non_default_data_from_elem_comp(self):
        comp = thermo.ElemMolComp(**{'Mg': 1, 'Si': 1, 'O': 3})
        assert {'Mg':1, 'Si': 1, 'O': 3} == comp.data

    def test_should_get_nonzero_values_from_oxide_comp(self):
        comp = thermo.OxideMolComp(MgO=1, SiO2=1)
        assert not np.any(comp.values==0)

    def test_should_get_nonzero_values_from_elem_comp(self):
        comp = thermo.ElemMolComp(**{'Mg': 1, 'Si': 1, 'O': 3})
        assert not np.any(comp.values==0)

    def test_should_get_all_values_from_oxide_comp(self):
        comp = thermo.OxideMolComp(MgO=1, SiO2=1)
        assert np.any(comp.all_values==0)

    def test_should_get_all_values_from_elem_comp(self):
        comp = thermo.ElemMolComp(**{'Mg': 1, 'Si': 1, 'O': 3})
        assert np.any(comp.all_values == 0)

    def test_should_get_all_components_from_oxide_comp(self):
        comp = thermo.OxideMolComp(MgO=1, SiO2=1)
        assert len(comp.all_components) > len(comp.components)

    def test_should_get_all_components_from_elem_comp(self):
        comp = thermo.ElemMolComp(**{'Mg': 1, 'Si': 1, 'O': 3})
        assert len(comp.all_components) > len(comp.components)

    def test_should_get_nonzero_components_from_oxide_comp(self):
        comp = thermo.OxideMolComp(MgO=1, SiO2=1)
        all_data = pd.Series(comp.all_data)

        assert np.all(all_data[comp.components]>0)
        assert np.all(all_data[comp.zero_components]==0)

        assert UnorderedList(comp.all_components) == (
                list(comp.components) + list(comp.zero_components) )

    def test_should_get_nonzero_components_from_elem_comp(self):
        comp = thermo.ElemMolComp(**{'Mg': 1, 'Si': 1, 'O': 3})
        all_data = pd.Series(comp.all_data)

        assert np.all(all_data[comp.components]>0)
        assert np.all(all_data[comp.zero_components]==0)

        assert UnorderedList(comp.all_components) == (
                list(comp.components) + list(comp.zero_components) )

    def test_should_define_composition_by_oxide_molar_abundance(self):
        comp = thermo.OxideMolComp(MgO=1, SiO2=1)
        assert comp == {'MgO':1, 'SiO2':1}

    def test_should_compare_normalized_oxide_compositions(self):
        comp = thermo.OxideMolComp(MgO=1, SiO2=1)
        assert comp == {'MgO':0.5, 'SiO2':0.5}
        assert comp == {'MgO':3, 'SiO2':3}
        assert not comp == {'MgO':2, 'SiO2':0.5}

    def test_should_compare_normalized_elemental_compositions(self):
        comp = thermo.ElemMolComp(**{'Mg': 2, 'Si': 1, 'O': 4})
        assert comp == {'Mg': 2, 'Si': 1, 'O': 4}
        assert comp == {'Mg': 2/7, 'Si': 1/7, 'O': 4/7}

    def test_should_calculate_elemental_comp_from_oxides(self):
        assert thermo.OxideMolComp(Al2O3=1).elem_comp == {'Al': 2, 'O': 3}
        assert thermo.OxideMolComp(MgO=1, SiO2=1).elem_comp == {'Mg': 1, 'Si': 1, 'O': 3}

    def test_should_set_simple_oxide_comp_by_wt_or_mols(self):
        assert OxideWtComp(SiO2=100) == OxideMolComp(SiO2=1)
        assert OxideMolComp(MgO=1, SiO2=1) == OxideWtComp(
            MgO=40.3044, SiO2=60.0848)

    def test_should_compare_natural_comps_by_wt_or_mols(self):
        BSE = OxideWtComp(SiO2=44.95, TiO2=0.158, Al2O3=3.52,
                          Cr2O3=0.385, MnO=0.131, FeO=7.97, NiO=0.252,
                          MgO=39.50, CaO=2.79, Na2O=0.298, K2O=0.023,
                          P2O5=0.015)

        assert BSE == OxideMolComp(
            SiO2=38.8, TiO2=0.10256129117905899, Al2O3=1.7904987989927128,
            Cr2O3=0.13137471668468576, FeO=5.753338927767198, MnO=0.09577731951497395,
            MgO=50.82896713297678, NiO=0.17494113572627076, CaO=2.5802839102636868,
            Na2O=0.2493668785846322, K2O=0.012663821802104725, P2O5=0.0054807409820552275)

        assert BSE == OxideMolComp(
            SiO2=38.8, TiO2=0.1026, Al2O3=1.7905,
            Cr2O3=0.1314, FeO=5.7533, MnO=0.09578,
            MgO=50.829, NiO=0.1749, CaO=2.5803,
            Na2O=0.2494, K2O=0.01266, P2O5=0.0055)

        assert OxideMolComp(
            SiO2=38.8, TiO2=0.1026, Al2O3=1.7905,
            Cr2O3=0.1314, FeO=5.7533, MnO=0.09578,
            MgO=50.829, NiO=0.1749, CaO=2.5803,
            Na2O=0.2494, K2O=0.01266, P2O5=0.0055) == BSE
        # assert BSE == OxideMolComp(SiO2=38.80, Al2O3=1.79, FeO=5.75, MgO=50.83, CaO=2.58, Na2O=0.25)






