from __future__ import absolute_import, print_function, division
from builtins import object
import numpy as np
from scipy.misc import derivative
import pytest

from thermoengine import chem
from thermoengine import phases
from thermoengine import model

class TestRxn():

    def test_pure_phase_rxn_affinity(self, modelDB):
        # Test Ferrosilite breakdown reaction
        # should occur between 1 and 1.5 GPa
        phase_symbols = ['Fa','Qz','Fs']
        endmember_ids = [0, 0, 0]
        rxn_coefs = [+1, +1, -2]

        rxn = modelDB.get_rxn(phase_symbols, endmember_ids, rxn_coefs)

        T = 1000.0 - 273.15

        P0 = 1
        affinity = rxn.affinity(T, P0)
        affinity_1GPa = rxn.affinity(T, 1e4)
        affinity_1_5GPa = rxn.affinity(T, 1.5e4)

        assert affinity_1GPa > 0, (
            'Ferrosilite breakdown is favored at low pressure, so affinity should be positive 1.0 GPa.'
        )
        assert affinity_1_5GPa < 0, (
            'Ferrosilite breakdown is inhibited at high pressure, so affinity should be negative at 1.5 GPa.'
        )

        phase_objs = rxn.phases
        phase_objs[0].chem_potential(T, P0)

        chem_potentials = [phase_obj.chem_potential(T, P0, endmember=0) for phase_obj in phase_objs]

        chem_potential_sum = np.sum(np.hstack(chem_potentials)*rxn_coefs)
        affinity_sum = - chem_potential_sum

        TOL = 1e-10
        assert np.abs(affinity_sum-affinity) < TOL, (
            'affinity calculation does not agree with manual chem potential sum.'
        )

    def test_pure_solution_rxn_affinity(self, modelDB):
        # Test Albite breakdown reaction
        #    Albite   =  Jadeite  + Quartz
        #  NaAlSi3O8  = NaAlSi2O6 +  SiO2
        # should occur around ~2 GPa at 1000 degrees C
        phase_Fsp = modelDB.get_phase('Fsp')
        phase_Cpx = modelDB.get_phase('Cpx')

        # NOTE: This is fairly awkward but correct
        albite_mol_oxide_comp = chem.format_mol_oxide_comp(
            {'Na2O':0.5,'Al2O3':0.5,'SiO2':3})
        jadeite_mol_oxide_comp = chem.format_mol_oxide_comp(
            {'Na2O':0.5,'Al2O3':0.5,'SiO2':2})

        Cpx_endmember_ind, Cpx_endmember_comp = (
            phase_Cpx.get_endmember_ind(jadeite_mol_oxide_comp,
                                        get_endmember_comp=True))

        Fsp_endmember_ind, Fsp_endmember_comp = (
            phase_Fsp.get_endmember_ind(albite_mol_oxide_comp,
                                        get_endmember_comp=True))


        phase_symbols = ['Fsp','Cpx','Qz']
        rxn_coefs = [-1, +1, +1]
        endmember_ids = [Fsp_endmember_ind,
                         Cpx_endmember_ind, 0]

        rxn = modelDB.get_rxn(phase_symbols, endmember_ids,
                              rxn_coefs)


        mols = {}
        mols['Fsp'] = Fsp_endmember_comp
        mols['Cpx'] = Cpx_endmember_comp

        T = 1000

        affinity_1GPa = rxn.affinity(T, 1e4, mols=mols)
        affinity_3GPa = rxn.affinity(T, 3e4, mols=mols)

        assert affinity_1GPa < 0, (
            'Albite breakdown is inhibited at low pressure, so affinity should be negative at 1.0 GPa.'
        )
        assert affinity_3GPa > 0, (
            'Albite breakdown is favored at high pressure, so affinity should be positive at 3 GPa.'
        )

    def test_get_balanced_rxns(self, modelDB):
        relevant_phases = ['Cpx', 'Grt', 'Ol', 'Opx', 'Bt',
                           'Crn', 'Rt', 'Fsp', 'Cam',
                           'Zo', 'Crd', 'Ms', 'Qz', 'Ky', 'Sil',
                           'SplS', 'Mll', 'Nph', 'Ilm']

        phase_assemblage = modelDB.get_assemblage(relevant_phases)
        # assert False
        # endmember_comp_matrix = phase_assemblage.get_endmember_comp_matrix()




        # assert False
        # iphs
class TestAssemblage():
    def test_gibbs_energy_pure(self, modelDB):
        phase_symbols = ['Qz', 'En', 'Fo']
        assemblage = modelDB.get_assemblage(phase_symbols)
        T = 2100.0 - 273.15
        P = 15000.0

        G_all = assemblage.gibbs_energy_all(T, P)

    def test_chem_potential_solution(self, modelDB):
        phase_symbols = ['Grt', 'Cpx']
        assemblage = modelDB.get_assemblage(phase_symbols)
        T = 2100.0 - 273.15
        P = 15000.0

        mols = {}
        iphs=assemblage.phases[0]

        # assert False
        # iphs

class TestSolutionPhase():
    def test_all_endmember_chem_potentials(self, modelDB):

        T = [2100.0 - 273.15, 2500-273.15, 3000-273.15]
        P = 15000.0

        solution_phases = modelDB.phase_obj['solution']


        mu_endmembers = []
        for iphase_key in solution_phases:
            iphase = solution_phases[iphase_key]
            endmember_num = iphase.endmember_num
            # imu_endmembers = np.zeros(endmember_num)
            imu_endmembers = []

            for ijendmem in range(endmember_num):
                ijmol = np.zeros(endmember_num)
                ijmol[ijendmem] = 1.0
                ijmu = iphase.chem_potential(T, P, mol=ijmol)
                if len(ijmu) != len(ijmol):
                    assert False
                else:
                    print(len(ijmu))

                # imu_endmembers[ijendmem] = ijmu
                imu_endmembers.append(ijmu)

            # imu_endmembers[-1][:] = np.nan
            mu_endmembers.append(np.array(imu_endmembers))

        # assert False

        assert np.all(~np.hstack([
            np.isnan(imu).ravel() for imu in mu_endmembers])), (
                'Endmember evaluations of the chemical '
                'potential should not produce NaN values.'
            )

    def test_specific_endmember_chem_potentials(self, modelDB):

        T0 = 2100.0 - 273.15
        T = T0
        P = 15000.0
        self._specific_endmember_chem_potentials(modelDB, T, P)

        T = [T0, T0+100, T0+200, T0+300, T0+400, T0+500]
        P = 15000.0
        self._specific_endmember_chem_potentials(modelDB, T, P)

    def _specific_endmember_chem_potentials(self, modelDB, T, P):

        NPT = np.max((np.array(T).size,np.array(P).size))

        solution_phases = modelDB.phase_obj['solution']


        # Loop over solution phases and get chem potential of first endmember for every phase
        mu_values = []
        for iphase_key in solution_phases:
            iphase = solution_phases[iphase_key]
            endmember_num = iphase.endmember_num
            endmember=0
            # imu_endmembers = np.zeros(endmember_num)


            imol = np.ones(endmember_num)
            imu = iphase.chem_potential(T, P, mol=imol, endmember=endmember)

            assert np.array(imu).size==NPT,(
                'The output of chemical potential for a defined endmember '
                'must match number of PT points'
            )

            mu_values.append(imu)

        mu_values = np.array(mu_values)


        # assert mu_values.ndim==1, (
        #     'Chemical potential with endmember flag should return '
        #     'only a single float.'
        # )
        assert np.all(~np.hstack([
            np.isnan(imu).ravel() for imu in mu_values])), (
                'Endmember evaluations of the chemical '
                'potential should not produce NaN values.'
            )

class TestPurePhase():
    def test_enthalpy(self, modelDB, TOL=1e-3):
        obj = modelDB.get_phase('Fo')
        T = 2100.0 - 273.15
        P = 15000.0


        H_est = obj.gibbs_energy(T, P) + T*obj.entropy(T, P)
        H_act = obj.enthalpy(T, P)
        H_err = (H_est-H_act)/H_act
        err = H_err

        assert np.abs(err) < TOL, 'Enthalpy is not consistent with gibbs energy and entropy'

    def test_entropy(self, modelDB, TOL=1e-3):
        obj = modelDB.get_phase('Fo')
        T = 2100.0 - 273.15
        P = 15000.0

        fun = lambda T: obj.gibbs_energy(T, P)

        S_est = -derivative(fun, T)
        S_act = obj.entropy(T, P)
        S_err = (S_est-S_act)/S_act
        err = S_err

        assert np.abs(err) < TOL, 'Entropy is not consistent with dG/dT'

    def test_volume(self, modelDB, TOL=1e-3):
        obj = modelDB.get_phase('Fo')
        T = 2100.0 - 273.15
        P = 15000.0

        fun = lambda P: obj.gibbs_energy(T, P)

        V_est = derivative(fun, P)
        V_act = obj.volume(T, P)
        V_err = (V_est-V_act)/V_act
        err = V_err

        assert np.abs(err) < TOL, 'Volume is not consistent with dG/dP'

    def test_heat_capacity(self, modelDB, TOL=1e-3):
        obj = modelDB.get_phase('Fo')
        T = 2100.0 - 273.15
        P = 15000.0

        fun = lambda T: obj.entropy(T, P)

        Cp_est = T*derivative(fun, T)
        Cp_act = obj.heat_capacity(T, P)
        Cp_err = (Cp_est-Cp_act)/Cp_act
        err = Cp_err

        assert np.abs(err) < TOL, 'Heat capacity is not consistent with dS/dT'

    def test_heat_capacity_Tderiv(self, modelDB, TOL=1e-3):
        obj = modelDB.get_phase('Fo')
        T = 2100.0 - 273.15
        P = 15000.0

        fun = lambda T: obj.heat_capacity(T, P)

        dCpdT_est = derivative(fun, T)
        dCpdT_act = obj.heat_capacity(T, P, deriv={'dT':1})
        dCpdT_err = (dCpdT_est-dCpdT_act)/dCpdT_act
        err = dCpdT_err

        assert np.abs(err) < TOL, 'Heat capacity temp. deriv is not consistent with dCp/dT'

    def test_volume_Tderiv(self, modelDB, TOL=1e-3):
        obj = modelDB.get_phase('Fo')
        T = 2100.0 - 273.15
        P = 15000.0

        fun = lambda T: obj.volume(T, P)

        dVdT_est = derivative(fun, T)
        dVdT_act = obj.volume(T, P, deriv={'dT':1})
        dVdT_err = (dVdT_est-dVdT_act)/dVdT_act

        err = dVdT_err

        assert np.abs(err) < TOL, 'Volume temp. deriv is not consistent with dV/dT'

    def test_volume_Pderiv(self, modelDB, TOL=1e-3):
        obj = modelDB.get_phase('Fo')
        T = 2100.0 - 273.15
        P = 15000.0

        fun = lambda P: obj.volume(T, P)

        dVdP_est = derivative(fun, P)
        dVdP_act = obj.volume(T, P, deriv={'dP':1})
        dVdP_err = (dVdP_est-dVdP_act)/dVdP_act
        err = dVdP_err

        assert np.abs(err) < TOL, 'Volume press. deriv is not consistent with dV/dP'

    @pytest.mark.xfail(reason='Bug in Stixrude PhaseObjC')
    def test_PhaseOBJC_invalid_Stixrude_calc(self, modelDB_Stix, TOL=1e-3):
        """
        NOTE that this only fails for Stixrude database.
        """
        obj = modelDB_Stix.get_phase('Fo')
        T = np.float64(2100.0 - 273.15)
        P = np.float64(15000.0)

        T0 = T.copy()
        P0 = P.copy()

        eval_start = obj.volume(T, P, deriv={'dP':1, 'dT':1})
        dVdT_fun = lambda P: obj.volume(T, P, deriv={'dT':1})

        # d2VdTdP_act

        d2VdTdP_est = derivative(dVdT_fun, P)
        eval_end = obj.volume(T, P, deriv={'dP':1, 'dT':1})

        assert np.abs(eval_start - eval_end) < 1e-12, 'Repeating the same calculation (d2V/dTdP of Fo) yields different answers.'

    # def test_volume_mixed_deriv(self, modelDB_Stix, TOL=1e-3):
    #     obj = modelDB_Stix.get_phase('Fo')
    #     T = np.float64(2100.0 - 273.15)
    #     P = np.float64(15000.0)

#
    #     T0 = T.copy()
    #     P0 = P.copy()
#
    #     eval_start = obj.volume(T, P, deriv={'dP':1, 'dT':1})
    #     dVdT_fun = lambda P: obj.volume(T, P, deriv={'dT':1})
#
    #     # d2VdTdP_act
#
    #     d2VdTdP_est = derivative(dVdT_fun, P)
    #     eval_end = obj.volume(T, P, deriv={'dP':1, 'dT':1})
#
    #     assert False
#
#
#
    #     dVdP_fun = lambda T: obj.volume(T, P, deriv={'dP':1})
#
    #     d2VdTdP_err = d2VdTdP_est/chicken - 1
    #     err = d2VdTdP_err
#
    #     assert False
    #     assert np.abs(err) < TOL, 'Mixed partial volume deriv (T-first) is not consistent with d2V/dTdP'
#
    #     d2VdPdT_est = derivative(dVdP_fun, T)
    #     d2VdPdT_err = d2VdPdT_est/d2VdTdP_act - 1
    #     err = d2VdPdT_err
    #     assert np.abs(err) < TOL, 'Mixed partial volume deriv (P-first) is not consistent with d2V/dTdP'


# def test_dvdp_dt(t,p):
#     d2vdtdp_est = derivative(dvdp, t, args=(True,))
#     d2vdtdp_act = obj.volume(t,p, deriv={'dT':1, 'dP':1})
#     d2vdtdp_err = (d2vdtdp_est-d2vdtdp_act)*100.0/d2vdtdp_act
#     print ("d2VdTDp {0:10.6f} % error, est: {1:15.6e} act: {2:15.6e}".format(d2vdtdp_err, d2vdtdp_est, d2vdtdp_act))

# def test_dvdp_dp(t,p):
#     d2vdp2_est = derivative(dvdp, p, args=(False,))
#     d2vdp2_act = obj.volume(t,p, deriv={'dP':2})
#     d2vdp2_err = (d2vdp2_est-d2vdp2_act)*100.0/d2vdp2_act
#     print ("d2VdP2  {0:10.6f} % error, est: {1:15.6e} act: {2:15.6e}".format(d2vdp2_err, d2vdp2_est, d2vdp2_act))
#     modelDBStix = model.Database(database='Stixrude')
# obj = modelDBStix.get_phase('Fo')
