import numpy as np
import pandas as pd

import thermoengine as thermo
from thermoengine import samples as smp
from thermoengine.samples import SampleMaker, SampleLibrary, Assemblage, MonophaseAssemblage, SampleMesh

from thermoengine import UnorderedList
from thermoengine.test_utils import are_close, all_arrays_are_close


modelDB = thermo.model.Database()
stixrudeDB = thermo.model.Database(database='Stixrude')
XTOL = 1e-5


class TestPhaseSample:
    def test___should_get_sample_with_sample_maker(self):
        Qz, Cpx, Fsp = modelDB.get_phases(['Qz', 'Cpx', 'Fsp'])

        smp._PurePhaseSample(Qz)
        assert (smp.SampleMaker.get_sample(Qz) ==
                smp._PurePhaseSample(Qz))
        assert (smp.SampleMaker.get_sample(Cpx) ==
                smp._SolutionPhaseSample(Cpx))
        assert (smp.SampleMaker.get_fixed_comp_sample(Fsp) ==
                smp._FixedCompPhaseSample(Fsp))

    def test_should_return_pure_phase_sample_properties_using_atomic_basis(self):
        Qz = modelDB.get_phase('Qz')
        atom_num = 3
        qz_samp = smp.SampleMaker.get_sample(Qz)
        assert qz_samp.G == Qz.gibbs_energy(qz_samp.T, qz_samp.P) / atom_num

    def test_should_return_solution_phase_sample_properties_using_atomic_basis(self):
        Fsp = modelDB.get_phase('Fsp')
        atom_num = 13
        fsp_samp = smp.SampleMaker.get_sample(Fsp, X=[0.5, 0.5, 0])
        assert fsp_samp.G == Fsp.gibbs_energy(fsp_samp.T, fsp_samp.P, mol=fsp_samp.X) / atom_num
        assert are_close(
            fsp_samp.chem_potential(),
            Fsp.chem_potential(fsp_samp.T, fsp_samp.P, mol=fsp_samp.X) / atom_num)

    def test_should_store_fixed_comp(self):
        qz_samp = smp.SampleMaker.get_sample(modelDB.get_phase('Qz'))
        fsp_samp = smp.SampleMaker.get_sample(modelDB.get_phase('Fsp'))
        fixed_fsp_samp = smp.SampleMaker.get_fixed_comp_sample(modelDB.get_phase('Fsp'))

        assert qz_samp.fixed_comp
        assert not fsp_samp.fixed_comp
        assert fixed_fsp_samp.fixed_comp

    def test_should_differentiate_solution_samples_and_fixed_comp_samples(self):
        Fsp = modelDB.get_phase('Fsp')
        fsp_samp = smp.SampleMaker.get_sample(Fsp)
        fixed_fsp_samp = smp.SampleMaker.get_fixed_comp_sample(Fsp)

        assert not fixed_fsp_samp == fsp_samp

    def test_should_compare_approx_sample_comps(self):
        Fsp = modelDB.get_phase('Fsp')
        X0 = np.array([0.8, 0.1, 0.1])
        dX = 1e-7*np.array([+1, -1, 0])

        samp = SampleMaker.get_fixed_comp_sample(Fsp, X=X0)
        samp_almost_equal = SampleMaker.get_fixed_comp_sample(Fsp, X=X0+dX)

        assert samp == samp_almost_equal

    def test_should_output_sample_summary(self):
        Fsp = modelDB.get_phase('Fsp')
        samp = smp.SampleMaker.get_fixed_comp_sample(Fsp, X=[0.7,0.3,0])

        value = '\n'.join((
            'Feldspar',
            '| ----------- endmember mol frac ------------ |',
            '|         albite [NaAlSi3O8]       0.7000 mol |',
            '|      anorthite [CaAl2Si2O8]      0.3000 mol |',
            '|       sanidine [KAlSi3O8]        0.0000 mol |',
        ))

        print()
        print(samp.summary)

        assert str(samp.summary) == value

class TestSampleMesh:
    def test_should_find_correct_number_of_neighboring_mesh_points(self):
        neighbors_2dim = SampleMesh.calc_neighboring_mesh_pts(
            X0=[1 / 2, 1 / 2], spacing=1 / 4)
        neighbors_3dim = SampleMesh.calc_neighboring_mesh_pts(
            X0=[1 / 3, 1 / 3, 1 / 3], spacing=1 / 6)
        neighbors_4dim = SampleMesh.calc_neighboring_mesh_pts(
            X0=[1 / 4, 1 / 4, 1 / 4, 1 / 4], spacing=1 / 8)

        assert len(neighbors_2dim) == 2
        assert len(neighbors_3dim) == 6
        assert len(neighbors_4dim) == 12

    def test_should_find_neighboring_mesh_points(self):
        neighbors_2dim = SampleMesh.calc_neighboring_mesh_pts(
            X0=[1 / 2, 1 / 2], spacing=0.1)
        assert are_close(neighbors_2dim,
                         [[0.6, 0.4],
                          [0.4, 0.6]])

        neighbors_3dim = SampleMesh.calc_neighboring_mesh_pts(
            X0=[1 / 3, 1 / 3, 1 / 3], spacing=1 / 6)
        assert are_close(neighbors_3dim,
                         [[3/6, 1/6, 1/3],
                          [1/6, 3/6, 1/3],
                          [3/6, 1/3, 1/6],
                          [1/6, 1/3, 3/6],
                          [1/3, 3/6, 1/6],
                          [1/3, 1/6, 3/6]])

    def test_should_only_return_viable_neighboring_mesh_points(self):
        endmem_neighbors_3dim = SampleMesh.calc_neighboring_mesh_pts(
            X0=[1, 0, 0], spacing=1 / 2)
        assert len(endmem_neighbors_3dim) == 2

        endmem_neighbors_4dim = SampleMesh.calc_neighboring_mesh_pts(
            X0=[0, 1, 0, 0], spacing=1 / 2)
        assert len(endmem_neighbors_4dim) == 3

        bounding_binary_neighbors_3dim = SampleMesh.calc_neighboring_mesh_pts(
            X0=[1 / 2, 1 / 2, 0], spacing=1 / 4)
        assert len(bounding_binary_neighbors_3dim) == 4

        bounding_ternary_neighbors_4dim = SampleMesh.calc_neighboring_mesh_pts(
            X0=[1 / 3, 1 / 3, 1 / 3, 0], spacing=1 / 6)
        assert len(bounding_ternary_neighbors_4dim) == 9

    def test_should_refine_local_mesh_for_interior_pts(self):
        X0 = [1/2, 1/4, 1/4]
        Xneighbors0 = np.array(
            [[3/4, 0/4, 1/4],
             [1/4, 2/4, 1/4],
             [3/4, 1/4, 0.],
             [1/4, 1/4, 2/4],
             [2/4, 2/4, 0/4],
             [2/4, 0/4, 2/4]])
        Xlocal0 = np.vstack((X0, Xneighbors0))
        midpt_count = 6 + 6

        local_mesh_3dim = SampleMesh.refine_local_mesh(
            Xlocal0=Xlocal0, spacing0=1 / 4)

        assert len(local_mesh_3dim) == midpt_count

    def test_should_refine_local_mesh_for_corner_pts(self):
        X0 = [1, 0, 0]
        Xneighbors0 = np.array(
            [[1/2, 1/2, 0],
             [1/2, 0, 1/2]])
        Xlocal0 = np.vstack((X0, Xneighbors0))
        midpt_count = 2 + 1

        local_mesh_3dim = SampleMesh.refine_local_mesh(
            Xlocal0=Xlocal0, spacing0=1 / 2)

        assert len(local_mesh_3dim) == midpt_count

    def test_should_refine_local_mesh_for_edge_pts(self):
        X0 = [1/2, 1/2, 0]
        Xneighbors0 = np.array(
            [[1, 0, 0],
             [0, 1, 0],
             [0, 1/2, 1/2],
             [1/2, 0, 1/2]])
        Xlocal0 = np.vstack((X0, Xneighbors0))
        midpt_count = 4 + 3

        local_mesh_3dim = SampleMesh.refine_local_mesh(
            Xlocal0=Xlocal0, spacing0=1 / 2)

        assert len(local_mesh_3dim) == midpt_count

    def test_should_prevent_duplicated_samples_in_assemblage(self):
        X0_samples = [
            [0.2, 0.1, 0.6, 0., 0.1],
            [0.3, 0.1, 0.5, 0., 0.1]]

        Xlocal0, Xgrid_refined = SampleMesh.refine_mesh_for_multiple_samples(
            X0_samples, spacing0=0.1)

        assert len(Xlocal0) < 32
        assert len(Xgrid_refined) < 116

        assert len(Xlocal0) == 27
        assert len(Xgrid_refined) == 101

        Xlocal0_uniq = np.unique(np.round(Xlocal0, decimals=3), axis=0)
        Xgrid_refined_uniq = np.unique(np.round(Xgrid_refined, decimals=3), axis=0)

        assert are_close(np.sort(Xlocal0), np.sort(Xlocal0_uniq))
        assert are_close(np.sort(Xgrid_refined), np.sort(Xgrid_refined_uniq))

    def test_should_build_binary_mesh_with_desired_spacing(self):
        binary = SampleMesh.build_sample_mesh(ndim=2, spacing=0.25)
        binary_expected = np.array([[0,1],[.25,.75],[.5,.5],[.75,.25],[1,0]])

        assert are_close(np.sort(binary), np.sort(binary_expected))

    def test_should_build_ternary_mesh_with_desired_spacing(self):
        ternary = SampleMesh.build_sample_mesh(ndim=3, spacing=0.25)
        ternary_expected = np.array(
            [[0.  , 0.  , 1.  ],
             [0.25, 0.  , 0.75],
             [0.5 , 0.  , 0.5 ],
             [0.75, 0.  , 0.25],
             [1.  , 0.  , 0.  ],
             [0.  , 0.25, 0.75],
             [0.25, 0.25, 0.5 ],
             [0.5 , 0.25, 0.25],
             [0.75, 0.25, 0.  ],
             [0.  , 0.5 , 0.5 ],
             [0.25, 0.5 , 0.25],
             [0.5 , 0.5 , 0.  ],
             [0.  , 0.75, 0.25],
             [0.25, 0.75, 0.  ],
             [0.  , 1.  , 0.  ]])

        assert are_close(ternary, ternary_expected)

    def test_should_build_high_dim_mesh_with_desired_spacing(self):
        mesh = SampleMesh.build_sample_mesh(ndim=5, spacing=0.1)

        uniq_mesh = np.unique(mesh, axis=0)
        X_missing = [0.2, 0.1, 0.6, 0.1, 0.]

        assert mesh.shape == uniq_mesh.shape
        assert X_missing in mesh.tolist()
        assert len(mesh) == 1001

    def test_should_determine_similar_mol_comps_equiv_within_tol(self):
        X1 = np.array([0.07692308, 0.07692308, 0.61538462, 0.07692308, 0.15384615])
        X2 = np.array([0.077, 0.077, 0.615, 0.077, 0.154])
        assert SampleMesh.mol_comps_equiv(X1, X2, TOL=1e-3)

    def test_should_determine_mol_comps_with_diff_dimensions_not_equiv(self):
        assert not SampleMesh.mol_comps_equiv([0], [0,0])

class TestSampleComps:
    def test___should_compare_samples_for_equivalence(self):
        fsp = modelDB.get_phase('Fsp')
        kspar_samp = smp.SampleMaker.get_sample(fsp, X=[0.5, 0, 0.5])

        assert kspar_samp == smp.SampleMaker.get_sample(fsp, X=[0.5, 0, 0.5])
        assert not kspar_samp == smp.SampleMaker.get_sample(fsp, X=[0.5, 0, 0.5], P = 1e3)
        assert not kspar_samp == smp.SampleMaker.get_sample(fsp, X=[0.75, 0, 0.25])

    def test_should_store_individual_elemental_composition_of_solution_phases_sample(self):
        fsp = modelDB.get_phase('Fsp')
        kspar_samp = smp.SampleMaker.get_sample(fsp, X=[0.5, 0, 0.5])

        assert kspar_samp.comp == thermo.ElemMolComp(Na=0.5, Ca=0, K=0.5, Al=1, Si=3, O=8)

    def test___should_store_compositions_for_set_of_solution_phase_samples(self):
        fsp = modelDB.get_phase('Fsp')
        fsp_samp = [smp.SampleMaker.get_sample(fsp, X=[1, 0, 0]),
                    smp.SampleMaker.get_sample(fsp, X=[0, 1, 0]),
                    smp.SampleMaker.get_sample(fsp, X=[0, 0, 1])]
        fsp_samp_comp = UnorderedList([samp.comp for samp in fsp_samp])
        assert fsp_samp_comp == [
            thermo.ElemMolComp(K=0, Na=1, Ca=0, Al=1, Si=3, O=8),
            thermo.ElemMolComp(K=0, Na=0, Ca=1, Al=2, Si=2, O=8),
            thermo.ElemMolComp(K=1, Na=0, Ca=0, Al=1, Si=3, O=8)]

    def test_should_get_samples_with_endmember_comp_for_pure_phase(self):
        Qz_samps = smp.SampleMaker.get_sample_endmembers(modelDB.get_phase('Qz'))
        assert len(Qz_samps) == 1
        assert are_close([samp.X for samp in Qz_samps], [1], abs_tol=XTOL)

    def test_should_get_samples_with_endmember_comps_for_solution_phase(self):
        Fsp_samps = smp.SampleMaker.get_sample_endmembers(modelDB.get_phase('Fsp'))
        X_samps = np.vstack([samp.X for samp in Fsp_samps])
        assert len(Fsp_samps) == 3
        assert are_close(X_samps, [[1, 0, 0], [0, 1, 0], [0, 0, 1]],
                         abs_tol=XTOL)

    def test_should_get_sample_grid_over_compositional_binary(self):
        assert len(smp.SampleMaker.get_sample_grid(
            modelDB.get_phase('Bt'), grid_spacing=1/5)) == 6

    def test_should_get_sample_grid_over_compositional_ternary(self):
        ternary_samps = smp.SampleMaker.get_sample_grid(
            modelDB.get_phase('Fsp'), grid_spacing=1/4)
        assert len(ternary_samps) == 15
        assert np.unique([samp.X for samp in ternary_samps], axis=0).shape[0] == 15


    def test___should_build_single_sample_comp_grid_for_pure_phase(self):
        pure_comp_grid = smp.SampleMaker.build_comp_grid(1, 1)
        assert len(pure_comp_grid) == 1
        assert pure_comp_grid == 1


    def test___should_build_binary_composition_grid_with_desired_spacing(self):
        binary_grid = smp.SampleMaker.build_comp_grid(2, 0.1)
        X0 = np.sort(binary_grid[:, 1])
        assert are_close(np.diff(X0), 0.1)

    def test___should_build_ternary_composition_grid_with_desired_spacing(self):
        ternary_grid = smp.SampleMaker.build_comp_grid(3, 0.2)
        endmember_grids = [np.unique(Xi) for Xi in ternary_grid.T]
        assert are_close(np.diff(endmember_grids), 0.2)

    def test___should_build_Ndimensional_composition_grid_with_desired_spacing(self):
        ND_grid = smp.SampleMaker.build_comp_grid(5, 1/3)
        endmember_grids = [np.unique(Xi) for Xi in ND_grid.T]
        assert are_close(np.diff(endmember_grids), 1/3)

    def test___should_build_composition_grid_with_normalized_comps(self):
        ternary_grid = smp.SampleMaker.build_comp_grid(3, 0.2)
        assert np.all(ternary_grid.sum(axis=1) == 1)

        ND_grid = smp.SampleMaker.build_comp_grid(5, 1 / 4)
        assert np.all(ND_grid.sum(axis=1) == 1)

    def test___should_build_ternary_grid_with_desired_resolution(self):
        ternary_grid_1 = smp.SampleMaker.build_comp_grid(3, 1)
        assert len(ternary_grid_1) == 3

        ternary_grid_2 = smp.SampleMaker.build_comp_grid(3, 1/2)
        assert len(ternary_grid_2)== 6

        ternary_grid_3 = smp.SampleMaker.build_comp_grid(3, 1/3)
        assert len(ternary_grid_3) == 10

    def test___should_build_quaternary_grid_with_desired_resolution(self):
        ternary_grid_1 = smp.SampleMaker.build_comp_grid(4, 1)
        assert len(ternary_grid_1) == 4

        ternary_grid_2 = smp.SampleMaker.build_comp_grid(4, 1/2)
        assert len(ternary_grid_2)== 10

class TestSampleLibrary:
    def test_should_store_copy_of_initialized_samples(self):
        samp_lib_empty = SampleLibrary([])
        assert len(samp_lib_empty) == 0

        samp_list = SampleMaker.get_sample_grid(modelDB.get_phase('Grt'))
        samp_lib = SampleLibrary(samp_list)

        assert len(samp_lib) > 0
        assert samp_lib[0] == samp_list[0]
        assert samp_lib[0] is not samp_list[0]

    def test_should_store_original_samples_if_not_safe_copy(self):
        original_samp_list = SampleMaker.get_sample_grid(modelDB.get_phase('Grt'))
        init_size = len(original_samp_list)
        samp_lib = SampleLibrary(original_samp_list, safe_copy=False)

        assert samp_lib[0] is original_samp_list[0]

        samp_lib.extend([SampleMaker.get_sample(modelDB.get_phase('Qz'))])
        assert len(samp_lib) == init_size + 1
        assert len(original_samp_list) == init_size

    def test_should_provide_slices_that_are_still_sample_libraries(self):
        samp_list = SampleMaker.get_sample_grid(modelDB.get_phase('Grt'))
        samp_lib = SampleLibrary(samp_list)

        assert type(samp_lib) == SampleLibrary
        assert type(samp_lib[:5]) == SampleLibrary

    def test_should_return_elem_comps_for_every_sample(self):
        phases = modelDB.get_phases(['En','Qz','Fa','Grt'])
        samples = [SampleMaker.get_sample(phs) for phs in phases]
        samp_lib = SampleLibrary(samples)

        elem_comps = pd.DataFrame([
            [0.600000, 0.2, 0.0, 0.200000, 0.000000],
            [0.666667, 0.0, 0.0, 0.333333, 0.000000],
            [0.571429, 0.0, 0.0, 0.142857, 0.285714],
            [0.600000, 0.0, 0.1, 0.150000, 0.150000]
        ], index=['En', 'Qz', 'Fa', 'Grt'],
            columns=['O', 'Mg', 'Al', 'Si', 'Fe'])

        assert are_close(samp_lib.elem_comps, elem_comps)

    def test_should_get_subset_retaining_sample_affinities(self):
        samp_list = [SampleMaker.get_sample(phs)
                     for phs in stixrudeDB.get_phases(['Fo','MgWds','MgRwd'])]
        samp_lib = SampleLibrary(samp_list)

        affinities = np.array([0, 1e4, 2e4])
        samp_lib.update_affinities(affinities)
        ind_subset = [0,1]

        assert are_close(samp_lib.affinities, affinities)
        assert are_close(samp_lib.get_subset(ind_subset).affinities,
                         affinities[ind_subset])

    def test_should_get_subset_retaining_sample_amounts(self):
        samp_list = [SampleMaker.get_sample(phs)
                     for phs in stixrudeDB.get_phases(['Fo','MgWds','MgRwd'])]
        samp_lib = SampleLibrary(samp_list)

        amounts = np.array([0.0, 1.1, 2.2])
        samp_lib.update_amounts(amounts)
        ind_subset = [0,1]

        assert are_close(samp_lib.amounts, amounts)
        assert are_close(samp_lib.get_subset(ind_subset).amounts,
                         amounts[ind_subset])

    # TODO: nearly stable samples
    # TODO: nearly stable phases
class TestAssemblage:
    def test_should_create_assemblage_from_list_of_samples(self):
        samples = [smp.SampleMaker.get_sample(phs)
                   for phs in modelDB.get_phases(['Qz','Ol','Cpx','Fsp'])]
        assem = smp.Assemblage(samples)

        assert assem == smp.Assemblage(samples)
        assert not assem == smp.Assemblage()

    def test_should_retrieve_sample_properties_from_assemblage(self):

        samples = [SampleMaker.get_sample(phs)
                   for phs in modelDB.get_phases(['Qz', 'Ol', 'Cpx', 'Fsp'])]
        assem = smp.Assemblage(samples)


        assert are_close(assem.energies,
                         [-327577.07980718615, -280701.47311117913,
                          -344672.3421547144, -328799.95127329574])

        assert are_close(assem.sample_amounts, [1, 1, 1, 1])

        assert all_arrays_are_close(assem.sample_endmem_comps,
                                    [np.array([1.]),
                                     np.array([1., 0., 0., 0., 0., 0.]),
                                     np.array([1., 0., 0., 0., 0., 0., 0.]),
                                     np.array([1., 0., 0.])])

    def test_should_calculate_total_comp_of_samples(self):
        samples = [SampleMaker.get_sample(phs)
                   for phs in modelDB.get_phases(['Qz', 'Fo'])]

        assem = smp.Assemblage(samples)
        assem.update_amounts([2, 1])

        Qz_comp = thermo.ElemMolComp(Si=1, O=2).normalize()
        Fo_comp = thermo.ElemMolComp(Mg=2, Si=1, O=4).normalize()
        assert assem.total_comp == 2*Qz_comp + Fo_comp

    def test_should_get_samples_by_phase(self):
        Fsp = modelDB.get_phase('Fsp')
        Qz = modelDB.get_phase('Qz')
        Qz_sample = SampleMaker.get_sample(Qz)
        Fsp_samples = [SampleMaker.get_sample(Fsp, X=[0,1,0]),
                       SampleMaker.get_sample(Fsp, X=[0,0,1])]
        assem = Assemblage(Fsp_samples+[Qz_sample])

        assert assem.get_subset_for_phase('Fsp') == Assemblage(Fsp_samples)
        assert assem.get_subset_for_phase('Qz') == Assemblage([Qz_sample])

    def test_should_filter_duplicate_endmember_phases(self):
        Fsp = modelDB.get_phase('Fsp')
        Fsp_samples = SampleMaker.get_sample_endmembers(Fsp)

        endmembers = modelDB.get_phases(['Ab','An','Sa'])
        endmem_samples = [SampleMaker.get_sample(phs) for phs in endmembers]

        assem = Assemblage(Fsp_samples+endmem_samples)
        assem_filtered = assem.remove_redundant_endmembers()

        for endmem in endmem_samples:
            assert endmem not in assem_filtered

        assert assem_filtered == Assemblage(Fsp_samples)

    def test_should_segregate_resolved_samples_into_separate_assemblages(self):
        res = 0.1
        Fsp = modelDB.get_phase('Fsp')
        resolved_samples = [
            SampleMaker.get_sample(Fsp, X=[1, 0, 0]),
            SampleMaker.get_sample(Fsp, X=[0, 1, 0]),
            SampleMaker.get_sample(Fsp, X=[0, 0, 1])]
        resolved_assem_groups = (
            MonophaseAssemblage(resolved_samples).segregate_resolved_samples(res))
        assert len(resolved_assem_groups) == 3

    def test_should_group_unresolved_samples_into_single_assemblage(self):
        res = 0.1
        Fsp = modelDB.get_phase('Fsp')
        unresolved_samples = [
            SampleMaker.get_sample(Fsp, X=[1, 0, 0]),
            SampleMaker.get_sample(Fsp, X=[.95, .05, 0]),
            SampleMaker.get_sample(Fsp, X=[.95, 0, .05])]
        unresolved_assem_groups = (
            MonophaseAssemblage(unresolved_samples).segregate_resolved_samples(res))
        assert len(unresolved_assem_groups) == 1

    def test_should_separate_partly_resolved_samples_into_diff_assemblages(self):
        res = 0.1
        Fsp = modelDB.get_phase('Fsp')
        partly_resolved_samples = [
            SampleMaker.get_sample(Fsp, X=[.5, .5, 0]),
            SampleMaker.get_sample(Fsp, X=[.55, .45, 0]),
            SampleMaker.get_sample(Fsp, X=[0, 0, 1])]
        partly_resolved_assem_groups = (
            MonophaseAssemblage(partly_resolved_samples).segregate_resolved_samples(res))
        assert len(partly_resolved_assem_groups) == 2

    def test_should_segregate_unresolved_sample_clusters_into_diff_assemblages(self):
        res = 0.1
        Fsp = modelDB.get_phase('Fsp')
        unresolved_clustered_samples = [
            SampleMaker.get_sample(Fsp, X=[.5, .5, 0]),
            SampleMaker.get_sample(Fsp, X=[.55, .45, 0]),
            SampleMaker.get_sample(Fsp, X=[0, .3, .7]),
            SampleMaker.get_sample(Fsp, X=[0, .35, .65])
        ]
        resolved_assem_groups = (
            MonophaseAssemblage(unresolved_clustered_samples).segregate_resolved_samples(res))
        assert len(resolved_assem_groups) == 2
