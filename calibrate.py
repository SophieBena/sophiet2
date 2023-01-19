"""The phases subpackage implements a Python interface with the Phase
Objective-C classes as well as the infrastructure for pure phase thermodynamic calibration.

"""

from thermoengine import core
from thermoengine import model
from thermoengine import chem

import numpy as np
import scipy as sp
import pandas as pd
from os import path
import re
import copy
from pandas import DataFrame

from collections import OrderedDict

# specialized numerical imports
from scipy import stats, random, interpolate as interp, optimize as optim
import numdifftools as nd

RGAS = 8.3144598
# __all__ = ['Database','Data','ParamModel','RxnModel']
__all__ = ['Database']

#===================================================
class Database:
    """
    Calibration database model object.


    Parameters
    ----------
    rxn_data: pandas df
        experimental data input
    modelDB: str
        choose thermodynamic model database (e.g., 'Berman' is a valid input)
    reaction_library:
    TOLsvd: int
    Ndraw: int
    ortho_scale: int
        ratio that dictates level of rxn complexity (orthogonality/simplicity)
    TOL: int
        complexity of rxns
    ignore_kinetics: boolean, default False
    contaminant_phases: str list
    rxn_trans_typ: ['logisitic']
    Tscl_energy: 1.0

    Returns
    -------
    rxn_svd: array of ints
        matrix of valid, lineraly independent reactions

    TODO
    ----
    - get trusted data?
    - priors?

    """

    RXN_FACTORS = ['flux', 'seed', 'contaminant']
    RXN_SCL_VALUES = [1, 1, 1]
    RXN_DEFAULT_VALUES = [1, 1, 0]
    PHASE_PARAM_LOOKUP = {'V0': 'V', 'S0': 'S', 'dH0': 'delta H'}

    T0 = 300.0
    P0 = 1.0

    def __init__(self, rxn_data, modelDB=None,reaction_library=None,
                 ignore_kinetics=False, contaminant_phases=None,
                 rxn_trans_typ='logistic', TOLsvd=1e-4, Ndraw=10,
                 ortho_scale=15, TOL=1e-10, phase_priors=None, rxn_priors=None,
                 ref_energy_phases=None):

        if modelDB is None:
            modelDB = model.Database()

        self._init_rxns(rxn_data, modelDB, reaction_library, Ndraw, TOLsvd, ortho_scale, TOL)

        self._init_params(phase_priors, rxn_priors, ref_energy_phases)

        self.modelDB = modelDB

        self.ignore_kinetics = ignore_kinetics
        self.contaminant_phases = contaminant_phases
        self.rxn_trans_typ = rxn_trans_typ

    def _init_rxns(self, rxn_data, modelDB, reaction_library, Ndraw, TOLsvd, ortho_scale, TOL):
        rxn_coefs = None
        rxn_eqns = None #Need to revisit this
        phases = None
        if reaction_library is None:
            #endmember_ids = [0,0]

            #rxns = [modelDB.get_rxn(irxn_prop['phases'], endmember_ids,
                                    #irxn_prop['coefs'])
                    #for irxn_prop in rxn_data.rxn_props]
            #rxn_eqns = [irxn_prop['eqn'] for irxn_prop in rxn_data.rxn_props]

            #phases = []
            #for irxn in rxns:
                #phases.extend(irxn.phases)

            #phases = np.unique(phases)

            phase_symbols=chem.get_phase_symbols(rxn_data)
            rxn_svd_props = chem.calc_reaction_svd(phase_symbols, TOLsvd=TOLsvd)
            rxn_svd = rxn_svd_props['rxn_svd']
            Nbasis=len(rxn_svd)
            wtcoefs, costs, rxn_coefs_raw, wtcoefs_ortho = chem.get_rxns(rxn_svd, Ndraw=Ndraw, ortho_scale=ortho_scale, Nbasis=Nbasis, TOL=TOL)
            rxn_coefs = rxn_coefs_raw.copy()
            (np.place(rxn_coefs, abs(rxn_coefs)< 1e-2, 0))
            #rxns =

            #endmember_ids = np.arange(0,len(rxn_coefs[0]))
            phases = phase_symbols

        else:
            assert False, 'reaction_library is not None, need to implement user defined set of reactions'

        self.reaction_library = reaction_library
        self.rxn_data = rxn_data
        self.rxn_coefs = rxn_coefs
        self.rxn_eqns = rxn_eqns
        self.rxns = rxns
        self.phases = phases

        #self.endmember_ids = endmember_ids
        #TK
    def _init_params(self, phase_priors, rxn_priors,
                     ref_energy_phases):

        self._param_values = []
        self._param_names = []
        self._param_scales = []

        self._set_ref_energy_phases(ref_energy_phases)

        self._init_rxn_params()
        self._init_phase_params()

        N = len(self._param_names)
        self._param_values = np.array(self._param_values)
        self._param_names = np.array(self._param_names)
        self._param_scales = np.array(self._param_scales)
        self._free_params = np.tile(False, (N,))

        self._set_phase_priors(phase_priors, ref_energy_phases)
        self._set_rxn_priors(rxn_priors)

    def _set_ref_energy_phases(self, ref_energy_phases):
        phases = self.phases
        #TK
        if phases is None:
            print('phases is None')
        else:
            phase_symbols = np.array(
                [iphs.abbrev for iphs in phases])
            print('phases is good to go')

        #phase_symbols = np.array(
            #[iphs.abbrev for iphs in phases])

        if ref_energy_phases is None:
            ref_energy_phases = [phase_symbols[0]]

        else:
            assert np.all([ref_phs in phase_symbols
                           for ref_phs in ref_energy_phases]),(
                               'The ref_energy_phases provided '
                               'are not valid.'
                               )

        H0_ref_phases = []
        for phase_name in ref_energy_phases:
            iref_phase, = phases[phase_symbols==phase_name]
            H0_ref_phases.append(iref_phase)

        self._ref_energy_phases = H0_ref_phases

    def _get_ref_energy(self, phase):
        ref_energy_phases = self._ref_energy_phases

        assert len(ref_energy_phases)==1, (
            'Currently, only a single ref_energy_phase is implimented. Must work to implement composition-dependent enthalpy ref.'
        )

        ref_phase = ref_energy_phases[0]

        H0, = ref_phase.get_param_values(param_names=['delta H'])
        return H0

    def _set_phase_priors(self, phase_priors, ref_energy_phases):
        phases = self.phases
        phase_symbols = np.array(
            [iphs.abbrev for iphs in phases])

        # from IPython import embed; embed(); import ipdb; ipdb.set_trace()

        prior_name = []
        prior_avg = []
        prior_error = []

        for ind, iphs in enumerate(phases):
            isym = iphs.abbrev
            imask = phase_priors['phase']==isym
            if np.any(imask):
                ipriors = phase_priors[imask]
                for iname, iavg, ierror in zip(
                    ipriors['param_name'],
                    ipriors['param_value'],
                    ipriors['param_error']):

                    if iname=='H0':
                        iparam_name = 'dH0_P'+str(int(ind))
                        # Adjust dH0 value here

                        iavg = iavg
                    else:
                        iparam_name = iname+'_P'+str(int(ind))

                    prior_name.append(iparam_name)
                    prior_avg.append(iavg)
                    prior_error.append(ierror)

        self._prior_name = np.array(prior_name)
        self._prior_avg = np.array(prior_avg)
        self._prior_error = np.array(prior_error)

    def _set_rxn_priors(self, rxn_priors):
        pass

    def _append_param(self, param_name, param_value, scale=1):
        self._param_names.append(param_name)
        self._param_values.append(param_value/scale)
        self._param_scales.append(scale)

    def _init_rxn_params(self):
        rxns = self.rxns

        for factor, scl_val in zip(self.RXN_FACTORS,
                                   self.RXN_SCL_VALUES):
            param_basename = 'alpha_'+factor+'_R'
            self._append_param(param_basename+'*', 0.0)

            for ind, irxn in enumerate(rxns):
                self._append_param(param_basename+str(ind), 0.0)

    def _init_phase_params(self):
        T0 = self.T0
        P0 = self.P0

        def get_set_phase_param(param_name):
            if param_name.startswith('dH0_P'):
                H0, iphs = self._get_phase_param(param_name, return_phase=True)
                H0_ref = self._get_ref_energy(iphs)
                # val0 = (H0 - H0_ref)/1e3
                val0 = (H0 - H0_ref)
                scale = np.abs(val0)
            else:
                val0 = self._get_phase_param(param_name)
                scale = np.abs(val0)

            self._append_param(param_name, val0, scale=scale)

        for ind, iphs in enumerate(self.phases):
            get_set_phase_param('V0_P'+str(ind))
            get_set_phase_param('S0_P'+str(ind))
            if iphs not in self._ref_energy_phases:
                get_set_phase_param('dH0_P'+str(ind))

    def _split_param_name(self, param_name):
        ind_sep = param_name.rfind('_')
        param_typ, phs_key = param_name[0:ind_sep], param_name[ind_sep+1:]
        return param_typ, phs_key

    def _get_phase_param(self, param_name, return_phase=False):
        param_typ, phs_key = self._split_param_name(param_name)

        assert phs_key[0]=='P', (
            'This parameter must be a phase parameter.'
        )

        phs_ind = int(phs_key[1:])
        iphs = self.phases[phs_ind]
        value = iphs.get_param_values(
            param_names=self.PHASE_PARAM_LOOKUP[param_typ])[0]

        if return_phase:
            return value, iphs
        else:
            return value

    def _update_phase_params(self):
        phase_param_names = ['delta H', 'V', 'S']
        calib_param_basenames = ['dH0_P', 'V0_P', 'S0_P']

        def _add_set_params(basename, phase, calib_names, calib_values,
                           phase_param_names, phase_param_values):
            phase_param_lookup = {'V0':'V', 'S0':'S', 'dH0':'delta H'}

            mask = [iname.startswith(basename) for iname in calib_names]
            if np.any(mask):
                val = calib_values[mask]
            else:
                val = None

            if val is not None and basename=='dH0':
                H0_ref = self._get_ref_energy(phase)
                # H0 = 1e3*(val + H0_ref)
                H0 = val + H0_ref
                val = H0

            if val is not None:
                phase_param_values.append(val)
                phase_param_names.append(phase_param_lookup[basename])

            pass

        for ind, iphase in enumerate(self.phases):
            icalib_param_names = self.get_param_group(kind='phase', id=ind)

            # icalib_param_values = self.param_values(
            #     param_group=icalib_param_names, scale_params=False)
            icalib_param_values = self.param_values(
                param_group=icalib_param_names, scale_params=True)

            iphase_param_names = []
            iphase_param_values = []

            _add_set_params('V0', iphase,
                            icalib_param_names, icalib_param_values,
                            iphase_param_names, iphase_param_values)
            _add_set_params('S0', iphase,
                            icalib_param_names, icalib_param_values,
                            iphase_param_names, iphase_param_values)
            _add_set_params('dH0', iphase,
                            icalib_param_names, icalib_param_values,
                            iphase_param_names, iphase_param_values)

            iphase.set_param_values(param_names=iphase_param_names,
                                    param_values=iphase_param_values)

        pass
    #===========================
    def _get_param_group_index(self, param_group):
        if param_group is None:
            param_group = self.get_param_group(free=True)
        else:
            param_group = np.array(param_group)

        loc = np.zeros(len(param_group), dtype=int)
        for ind, name in enumerate(param_group):
            iloc, = np.where(self._param_names==name)
            loc[ind] = iloc

        return loc

    def get_param_group(self, kind='all', id=None,
                        base=None, free=None):
        """
        kind: ['all', 'phase', 'rxn']
        id: [None, '*', int]
        base: [None, str]
        free: [None, True, False]

        """

        param_names = self._param_names

        def _get_param_mask(symbol, id):
            base = '.*_'+symbol
            if id is None:
                mask = np.array(
                    [re.match(base+'[\*0-9]*', iname) is not None
                     for iname in param_names])
            elif id == '*':
                mask = np.array(
                    [re.match(base+'\*', iname) is not None
                     for iname in param_names])
            else:
                mask = np.array(
                    [re.match(base+str(int(id)), iname) is not None
                     for iname in param_names])

            return mask

        if kind == 'all':
            mask = np.tile(True, param_names.size)
        elif kind == 'phase':
            mask = _get_param_mask('P', id)
        elif kind == 'rxn':
            mask = _get_param_mask('R', id)
        else:
            assert False, (
                'That is not a valid param_group kind.'
            )

        if free is None:
            pass
        elif free:
            mask = mask & self._free_params
        else:
            mask = mask & ~self._free_params

        if base is not None:
            mask = mask & np.array([name.startswith(base)
                                    for name in param_names])

        param_group = param_names[mask]
        return param_group

    def scale_params(self, param_values, param_group=None):
        param_scales = self.param_scales(param_group=param_group)
        scaled_param_values = param_values*param_scales
        return scaled_param_values

    def unscale_params(self, scaled_param_values, param_group=None):
        param_scales = self.param_scales(param_group=param_group)
        param_values = param_values/param_scales
        return param_values

    def param_names(self, param_group=None):
        ind = self._get_param_group_index(param_group)
        return self._param_names[ind]

    def param_scales(self, param_group=None):
        ind = self._get_param_group_index(param_group)
        return self._param_scales[ind]

    def param_values(self, param_group=None, scale_params=False):
        ind = self._get_param_group_index(param_group)
        param_values = self._param_values[ind]

        if scale_params:
            return self.scale_params(param_values, param_group=param_group)
        else:
            return param_values

    def param_errors(self, param_group=None):
        ind = self._get_param_group_index(param_group)
        return self._param_errors[ind]

    def set_param_values(self, param_values, param_group=None):
        ind = self._get_param_group_index(param_group)

        self._param_values[ind] = param_values
        self._update_phase_params()
        pass

    def add_free_params(self, param_group):
        ind = self._get_param_group_index(param_group)
        self._free_params[ind] = True

    def del_free_params(self, param_group):
        ind = self._get_param_group_index(param_group)
        self._free_params[ind] = False
        pass
    #===========================
    def rxn_affinity(self):
        rxns = self.rxns
        rxn_data = self.rxn_data

        P = rxn_data.conditions['P']
        T = rxn_data.conditions['T']
        rxn_id = rxn_data.rxn['rxn_id']

        rxn_affinity = np.zeros(len(P))
        for ind, (irxn_id, iT, iP) in enumerate(zip(rxn_id, T, P)):
            irxn = rxns[irxn_id]
            rxn_affinity[ind] = irxn.affinity(iT, iP)

        return rxn_affinity

    def rxn_affinity_error(self):
        rxns = self.rxns
        rxn_data = self.rxn_data

        P = rxn_data.conditions['P']
        T = rxn_data.conditions['T']
        P_err = rxn_data.conditions['P_err']
        T_err = rxn_data.conditions['T_err']
        rxn_id = rxn_data.rxn['rxn_id']

        affinity_err = np.zeros(len(P))

        # from IPython import embed; embed(); import ipdb; ipdb.set_trace()
        for ind, (irxn_id, iT, iP, iTerr, iPerr) in enumerate(zip(rxn_id, T, P, T_err, P_err)):
            irxn = rxns[irxn_id]
            irxn_vol = irxn.volume(iT, iP)
            irxn_entropy = irxn.entropy(iT, iP)

            affinity_err[ind] = np.sqrt(
                (irxn_vol*iPerr)**2 + (irxn_entropy*iTerr)**2
            )


        return affinity_err

    def rxn_affinity_thresh(self):
        rxns = self.rxns
        rxn_data = self.rxn_data

        P = rxn_data.conditions['P']
        T = rxn_data.conditions['T']
        P_err = rxn_data.conditions['P_err']
        T_err = rxn_data.conditions['T_err']
        rxn_id = rxn_data.rxn['rxn_id']

        affinity_thresh = np.zeros(len(P))

        if self.ignore_kinetics:
            return affinity_thresh


        # # get universal reaction parameters
        # alpha_t_all   = param_d['alpha_t_rxn_all']
        # alpha_T_all   = param_d['alpha_T_rxn_all']
        # alpha_H2O_all = param_d['alpha_H2O_rxn_all']

        # for ind,(rxn_eqn, rxn_obj) in enumerate(zip(rxn_eqn_l,rxn_l)):
        #     msk_rxn = dat_d['rxn']==rxn_eqn
        #     idGrxn_a = rxn_obj.get_rxn_gibbs_energy(dat_d['T'][msk_rxn],
        #                                             dat_d['P'][msk_rxn],
        #                                             peratom=True )
        #     # idVrxn_a =  rxn_obj.get_rxn_volume(dat_d['T'][msk_rxn],
        #     #                                    dat_d['P'][msk_rxn],
        #     #                                    peratom=True )
        #     # idSrxn_a =  rxn_obj.get_rxn_entropy(dat_d['T'][msk_rxn],
        #     #                                     dat_d['P'][msk_rxn],
        #     #                                     peratom=True )


        #     # get reaction-specific parameters
        #     alpha_0_rxn = param_d['alpha_0_rxn'+str(ind)]
        #     dalpha_t_rxn = param_d['dalpha_t_rxn'+str(ind)]
        #     dalpha_T_rxn = param_d['dalpha_T_rxn'+str(ind)]
        #     dalpha_H2O_rxn = param_d['dalpha_H2O_rxn'+str(ind)]

        #     logGth_a[msk_rxn] = alpha_0_rxn \
        #         + (alpha_t_all+dalpha_t_rxn)*dat_d['time'][msk_rxn] \
        #         + (alpha_T_all+dalpha_T_rxn)*dat_d['T'][msk_rxn] \
        #         + (alpha_H2O_all+dalpha_H2O_rxn)*dat_d['water'][msk_rxn]

        #     dGrxn_a[msk_rxn] = idGrxn_a
        #     # dVrxn_a[msk_rxn] = idVrxn_a
        #     # dSrxn_a[msk_rxn] = idSrxn_a

        # Gth_a = Gthresh_scl*np.exp(logGth_a)
        # # sigG_a = np.sqrt((dat_d['Perr']*dVrxn_a)**2+(dat_d['Terr']*dSrxn_a)**2)

        # loglk_a = np.zeros(Gth_a.size)


        #     loglk_a[msk_rxndir] = iloglk_a


        # # from IPython import embed; embed(); import ipdb; ipdb.set_trace()
        # for ind, (irxn_id, iT, iP, iTerr, iPerr) in enumerate(zip(rxn_id, T, P, T_err, P_err)):
        #     irxn = rxns[irxn_id]
        #     irxn_vol = irxn.volume(iT, iP)
        #     irxn_entropy = irxn.entropy(iT, iP)

        #     affinity_err[ind] = np.sqrt(
        #         (irxn_vol*iPerr)**2 + (irxn_entropy*iTerr)**2
        #     )

        return affinity_thresh

    def eval_model_costfun(self, param_values, param_group=None,
                           full_output=False, kind='logistic'):

        N = len(param_values)

        # log_prior = -0.5*((param_values-np.array([29,-125, -3.5]))/10)**2
        # log_prior = -0.5*((param_values-np.array([5.147, 4.412, 4.983]))/0.1)**2

        # log_prior = -0.5*((param_values-np.array([1,1,1]))/0.02)**2
        prior = np.sign(param_values)*np.ones(param_values.shape)
        log_prior = -0.5*((param_values-prior)/0.02)**2

        # param_values = 5*np.exp(param_values*1e-3)
        self.set_param_values(param_values, param_group=param_group)

        # What about trust values?

        affinity = self.rxn_affinity()
        affinity_err = self.rxn_affinity_error()
        affinity_thresh = self.rxn_affinity_thresh()

        rxn_dir = self.rxn_data.rxn['rxn_dir']
        log_like = Stats.logprob_rxn_dir(rxn_dir, affinity, affinity_err,
                                         affinity_thresh, kind=kind)
        # log_prior = -0.5*((param_values-np.array([29,-125, -3.5]))/1)**2

        log_like_tot = np.sum(log_like)
        log_prior_tot = np.sum(log_prior)

        # log_prior = self.eval_log_prior()

        cost_val = np.hstack((-log_like, -log_prior))
        cost_tot = - log_like_tot - log_prior_tot

        if full_output:
            output = OrderedDict()
            output['cost_tot'] = cost_tot
            output['cost_val'] = cost_val
            output['log_like'] = log_like
            output['log_prior'] = log_prior

            output['log_prior_tot'] = log_prior_tot
            output['log_like_tot'] = log_like_tot

            output['affinity'] = affinity
            output['affinity_err'] = affinity_err
            output['affinity_thresh'] = affinity_thresh
            return output
        else:
            return cost_tot

    def fit_model(self, param_group=None, full_output=False,
                  kind='logistic', method='Nelder-Mead'):

        # Extract only trustworthy data
        # self._dat_trust_d = self.extract_trust_data()

        params0 = self.param_values(param_group=param_group)
        # params0 = np.log(params0/5)/1e-3

        model_cost0 = self.eval_model_costfun(
            params0, param_group=param_group,
            kind=kind, full_output=True)
        # print(params0)
        # print(model_cost0)

        # param0_unscale_a = self.get_param_values(free_params)
        # param0_a = self.unscale_params(param0_unscale_a, free_params)
        # param0_tbl = self.get_param_table(param_nm_a=free_params)

        # Precalculate approx Gibbs energy uncertainties
        # sigG_trust_a = self.propagate_data_errors(param0_unscale_a,
        #                                           free_params=free_params)
        # self._sigG_trust_a = sigG_trust_a

        # costfun = lambda params: self.eval_model_costfun_scl(
        #     params0, free_params=free_params)

        # lnprob_f = lambda param_a: -costfun(param_a)

        costfun = lambda params: self.eval_model_costfun(
            params, param_group=param_group, kind=kind, full_output=False)

        method = 'Nelder-Mead'
        # method = 'BFGS'
        result = optim.minimize(costfun, params0, method=method,
                                options={'disp':True, 'maxiter':1e4})

        # set best-fit value
        self.set_param_values(result.x, param_group=param_group)

        return result

        # def shift_one_param(shift,ind,mu_a=result.x,costfun=costfun):
        #     param_a = np.copy(mu_a)
        #     param_a[ind] += shift
        #     return costfun(param_a)

        # # Create re-scaled-shifted function for hessian
        # mu_a = result.x
        # cost0 = costfun(mu_a)
        # delx_param_scl = np.zeros(mu_a.shape)
        # dcost_target=1
        # for ind,param in enumerate(mu_a):
        #     del0 = 1e-2
        #     idelcostfun = lambda dx, ind=ind,target=dcost_target: \
        #         shift_one_param(dx,ind)-cost0-dcost_target
        #     delx = optim.fsolve(idelcostfun,del0)
        #     delx_param_scl[ind] = np.abs(delx)


        # norm_costfun = lambda dx_a, shift_scl=delx_param_scl,\
        #     mu_a=mu_a,costfun=costfun: costfun(dx_a*shift_scl+mu_a)


        # curv_scl_a = delx_param_scl*self.get_param_scl_values(free_params)
        # scl_mat_a = core.make_scale_matrix(curv_scl_a)

        # Hnorm_fun = nd.Hessian(norm_costfun,step = 1e-2)
        # Hnorm_a = Hnorm_fun(np.zeros(mu_a.shape))

        # covnorm_a = np.linalg.pinv(Hnorm_a)

        # cov_a = covnorm_a*scl_mat_a

        # try:
        #     err_a = np.sqrt(np.diag(cov_a))
        #     # print(cov_a)
        #     err_scl_a = core.make_scale_matrix(err_a)
        #     corr_a = cov_a/err_scl_a
        # except:
        #     err_a = None
        #     corr_a = None

        # # MCMC
        # # ndim = len(free_params)
        # # nwalkers = 10*ndim
        # # walker_pos0_a = [result['x'] + 1e-4*np.random.randn(ndim)
        # #                  for i in range(nwalkers)]
        # # sampler = emcee.EnsembleSampler(nwalkers, ndim,lnprob_f)
        # # sampler.run_mcmc(walker_pos0_a,10)


        # # from IPython import embed; embed(); import ipdb; ipdb.set_trace()

        # paramf_unscale_a = self.scale_params(result.x, free_params)
        # self.set_param_values(paramf_unscale_a,free_params)
        # self._fit_result_d = result

        # if full_output:

        #     output_d = {}
        #     output_d['err_a'] = err_a
        #     output_d['corr_a'] = corr_a

        #     self._mu_a = paramf_unscale_a
        #     self._err_a = err_a
        #     self._corr_a = corr_a

        #     for key in self._param_error_d:
        #         self._param_error_d[key] = np.nan

        #     for key,err_val in zip(free_params,err_a):
        #         self._param_error_d[key] = err_val

        #     model_cost_d = self.eval_model_costfun(paramf_unscale_a,
        #                                            free_params=free_params,
        #                                            full_output=True)
        #     param_tbl = self.get_param_table(param_nm_a=free_params)
        #     param_all_tbl = self.get_param_table(typ='all')

        #     output_d['free_params'] = free_params

        #     output_d['costval0'] = model_cost0_d['cost_val']
        #     output_d['costdata0_df'] = model_cost0_d['cost_data_df']
        #     output_d['param0_tbl'] = param0_tbl

        #     output_d['costval'] = model_cost_d['cost_val']
        #     output_d['costdata_df'] = model_cost_d['cost_data_df']
        #     output_d['prior_df'] = model_cost_d['prior_df']
        #     output_d['param_tbl'] = param_tbl
        #     output_d['param_all_tbl'] = param_all_tbl

        #     output_d['result'] = result

        #     # output_d['param_d'] = copy.copy(self._param_d)
        #     # output_d['param0_a'] = param0_a
        #     # output_d['paramf_a'] = result.x
        #     # output_d['param0_unscl_a'] = param0_unscale_a
        #     # output_d['paramf_unscl_a'] = paramf_unscale_a

        #     return output_d

        pass
#===================================================
class Stats:
    @classmethod
    def logprior_fun(cls, x, kind='studentt', dof=5):
        if kind == 'studentt':
            # Variance of student's t distribution is slightly larger than a
            # normal (depending on dof). Thus the relative residual x must be
            # scaled down to match the desired standard deviation.
            const = np.sqrt(1.0*dof/(dof-2))
            log_prob = stats.t.logpdf(x/const,dof)
        elif (kind == 'normal') | (kind == 'erf'):
            log_prob = stats.norm.log_pdf(x)

        return log_prob

    @classmethod
    def rxn_trans_fun(self, x, kind='logistic'):
        if kind == 'logistic':
            const = 1.8138 # pi/sqrt(3)
            # F_a = 1.0/(1+np.exp(-const*x))
            # Special optimized version of logistic function
            prob = sp.special.expit(const*x)
        elif (kind == 'normal') | (kind == 'erf'):
            const = 0.70711 # 1/sqrt(2)
            prob = 0.5*(1+sp.special.erf(const*x))

        return prob

    @classmethod
    def rxn_logtrans_fun(cls, x, kind='logistic'):
        if kind == 'logistic':
            const = 1.8138 # pi/sqrt(3)
            # Special optimized version of logistic function
            log_prob = -np.logaddexp(0,-const*x)
        elif (kind == 'normal') | (kind == 'erf'):
            const = 0.70711 # 1/sqrt(2)
            # Special optimized version of log-cdf for normal distribution
            log_prob = sp.special.log_ndtr(x)

        return log_prob

    @classmethod
    def logprob_rxn_dir(cls, rxn_dir, affinity, affinity_err, affinity_thresh,
                        kind='logistic'):
        """
        rxndir = ['FWD', 'REV', 'NC', 'FWD?', 'REV?', 'NC?']
        """

        shp = affinity.shape

        x_fwd = (affinity-affinity_thresh)/affinity_err
        x_rev = -(affinity+affinity_thresh)/affinity_err

        ones = np.ones(shp)
        zeros = np.zeros(shp)
        log_prob = np.zeros(shp)

        log_prob[rxn_dir=='FWD'] = cls.rxn_logtrans_fun(
            x_fwd[rxn_dir=='FWD'], kind=kind)

        log_prob[rxn_dir=='REV'] = cls.rxn_logtrans_fun(
            x_rev[rxn_dir=='REV'], kind=kind)

        log_prob[rxn_dir=='BIASED'] = 0.0

        log_prob[rxn_dir=='FWD?'] = sp.special.logsumexp(
            np.vstack((
                cls.rxn_logtrans_fun(x_rev[rxn_dir=='FWD?'], kind=kind),
                zeros[rxn_dir=='FWD?'])), axis=0,
            b=np.vstack((-ones[rxn_dir=='FWD?'], +ones[rxn_dir=='FWD?'])))

        log_prob[rxn_dir=='REV?'] = sp.special.logsumexp(
            np.vstack((
                cls.rxn_logtrans_fun(x_fwd[rxn_dir=='REV?'], kind=kind),
                zeros[rxn_dir=='REV?'])), axis=0,
            b=np.vstack((-ones[rxn_dir=='REV?'], +ones[rxn_dir=='REV?'])))

        log_prob[rxn_dir=='NC'] = sp.special.logsumexp(
            np.vstack((
                cls.rxn_logtrans_fun(x_fwd[rxn_dir=='REV?'], kind=kind),
                cls.rxn_logtrans_fun(x_rev[rxn_dir=='REV?'], kind=kind),
                zeros[rxn_dir=='REV?'])), axis=0,
            b=np.vstack((-ones[rxn_dir=='REV?'], -ones[rxn_dir=='REV?'], +ones[rxn_dir=='REV?'])))

        log_prob[np.isnan(log_prob)] = -np.inf
        # from IPython import embed; embed(); import ipdb; ipdb.set_trace()

        return log_prob
#===================================================
class Database_OLD:
    def __init__(self, rxn_data, thermoDB=None, ignore_kinetics=False,
                 contaminant_phases=None, rxn_trans_typ='logistic',
                 Tscl_energy=1.0):

        if thermoDB is None:
            thermoDB = model.Database()

        self.rxn_data = rxn_data
        self.thermoDB = thermoDB
        self.ignore_kinetics = ignore_kinetics
        self.contaminant_phases = contaminant_phases
        self.rxn_trans_typ = rxn_trans_typ


        # self.datadir = DATADIR
        self.Escl = 3.0/2*RGAS*Tscl_energy


        # self.Gthresh_scl = 3.0/2*Rgas*Tthresh_scl
        # self.Gthresh_scl = self.Escl

        # dat_df, rxn_d_l, phasesym_l = self.filter_phase_rev_data(
        #     raw_df, mask_phs_l=mask_phs_l)
        # self._dat_trust_d = self.extract_trust_data()

        # self.load_exp_prior_data()


        self._param_d = {}
        self.init_model_phases(phasesym_l)
        self.init_model_rxns(rxn_d_l)
        self.init_exp_priors()

        self.init_param_scl()

        param_error_d = self._param_d.copy()
        for key in param_error_d:
            param_error_d[key] = np.nan

        self._param_error_d = param_error_d

        param_name_a = np.array(list(self._param_d.keys()))
        self._free_param_a = param_name_a

        self.load_exp_prior_data()

        # sigG_trust_a = self.propagate_data_errors(param0_unscale_a)
        # self._sigG_trust_a = sigG_trust_a
        self._sigG_trust_a = None

    def load_exp_prior_data( self ):
        datadir = self.datadir
        filenm = 'ExpPriorData.xlsx'


        parentpath = path.dirname(__file__)
        pathnm = path.join(parentpath,datadir,filenm)
        exp_prior_data_df = pd.read_excel(pathnm,sheetname=None)
        # Cycle through sheets, each representing a different parameter
        all_prior_df = pd.DataFrame()
        for paramnm in exp_prior_data_df:
            if (paramnm[0]=='<') & (paramnm[-1]=='>'):
                # This sheet provides metadata (such as references) rather than
                # actual priors
                continue

            data_df = exp_prior_data_df[paramnm]
            param_s = pd.Series(np.tile(paramnm,data_df.shape[0]))
            data_df['Param'] = param_s
            prior_data_df = data_df[['Phase','Abbrev','Param','Trust',
                                     'Data','Error','Ref']]

            all_prior_df = all_prior_df.append(prior_data_df)

            # phssym_a = data_df['Abbrev']
            # val_a = data_df['Data']
            # err_a = data_df['Error']
            # trust_a = data_df['Trust']
            # refID_a = data_df['RefID']

            # for sym in phssym_a:
            # prior_d = {'refID':refID_a,'phase_sym':phssym_a,'value':val_a,
            #            'error':err_a,'trust':trust_a}


        self.all_prior_df = all_prior_df

        pass

    def init_exp_priors(self):
        all_prior_df = self.all_prior_df
        exp_prior_df = pd.DataFrame()
        for ind,phs in enumerate(self.phs_key):
            prior_dat_phs_df = all_prior_df[all_prior_df['Abbrev'] == phs].copy()
            param_name_s = pd.Series(prior_dat_phs_df['Param']+str(ind))
            prior_dat_phs_df['Param'] = param_name_s
            exp_prior_df = exp_prior_df.append(prior_dat_phs_df,ignore_index=True)
            # print(prior_dat_phs_df)

        self.exp_prior_df = exp_prior_df

        pass

    def init_param_scl(self):
        param_name_a = np.array(list(self._param_d.keys()))
        # Define the scale for each type of parameter
        S0_param_keys_a = np.sort(param_name_a[np.char.startswith(param_name_a,'S0_phs')])
        V0_param_keys_a = np.sort(param_name_a[np.char.startswith(param_name_a,'V0_phs')])
        dH0_param_keys_a = np.sort(param_name_a[np.char.startswith(param_name_a,'dH0_phs')])

        S0_scl_atom = np.mean(self.get_param_values(S0_param_keys_a)/self.Natom_a)
        V0_scl_atom = np.mean(self.get_param_values(V0_param_keys_a)/self.Natom_a)
        # dH0_a = self.get_param_values(dH0_param_keys_a)/self.Natom_a
        # dH0_scl_atom = 3./2*Rgas*1e1
        dH0_scl_atom = self.Escl

        # alpha_T_scl = 1.0/1000 # 1/K
        alpha_T_scl = 1.0/1.0 # 1/K
        alpha_t_scl = 1.0
        alpha_H2O_scl = 1.0
        alpha_0_scl = 1.0

        param_scl_d = {}
        for ind,phs in enumerate(self.phase_l):
            Natom=phs.props_d['Natom']
            param_scl_d['S0_phs'+str(ind)] = S0_scl_atom*Natom
            param_scl_d['V0_phs'+str(ind)] = V0_scl_atom*Natom
            param_scl_d['dH0_phs'+str(ind)] = dH0_scl_atom*Natom



        param_scl_d['alpha_t_rxn_all']   =alpha_t_scl
        param_scl_d['alpha_T_rxn_all']   =alpha_T_scl
        param_scl_d['alpha_H2O_rxn_all'] =alpha_H2O_scl

        # logGth_rxn = log(3/2*Rgas*1) + alpha_i*X_i
        rxn_l = self.rxn_l
        for ind,rxn in enumerate(rxn_l):
            # self._param_d['logGth0_rxn'+str(ind)] = -
            param_scl_d['alpha_0_rxn'+str(ind)]   = alpha_0_scl
            param_scl_d['dalpha_t_rxn'+str(ind)]   = alpha_t_scl
            param_scl_d['dalpha_T_rxn'+str(ind)]   = alpha_T_scl
            param_scl_d['dalpha_H2O_rxn'+str(ind)] = alpha_H2O_scl

        self.param_scl_d = param_scl_d

        pass

    def init_model_phases( self, phasesym_l ):
        phase_l = [ self.thermoDB.new_phase(phasesym) for phasesym in phasesym_l]
        Natom_a = [phs.props_d['Natom'] for phs in phase_l]

        self._phasesym_l = phasesym_l
        self.phase_l = phase_l
        self.Natom_a = Natom_a
        self.init_std_state_params()

    def init_std_state_params(self):
        phase_l = self.phase_l
        phasesym_l = []
        for ind,phs in enumerate(phase_l):
            iparam_a = phs.get_param_values( param_names=['S','V','delta H'] )
            phasesym_l.append(phs.props_d['abbrev'])
            self._param_d['S0_phs'+str(ind)] = iparam_a[0]
            self._param_d['V0_phs'+str(ind)] = iparam_a[1]
            self._param_d['dH0_phs'+str(ind)] = iparam_a[2]

        pass

    def init_model_rxns(self, rxn_d_l):
        rxn_obj_l = []
        rxn_eqn_l = []
        for rxn_d in rxn_d_l:
            rxn_obj = self.thermoDB.new_rxn( rxn_d['reac_l'], rxn_d['prod_l'] )
            rxn_obj_l.append(rxn_obj)
            rxn_eqn_l.append(rxn_d['rxn_eqn'])

        self.rxn_l = rxn_obj_l
        self._rxn_eqn_l = rxn_eqn_l
        self.init_rxn_params()

    def init_rxn_params(self,dTthresh0=1.0):
        self._param_d['alpha_t_rxn_all']   = 0.0
        self._param_d['alpha_T_rxn_all']   = 0.0
        self._param_d['alpha_H2O_rxn_all'] = 0.0

        # logGth_rxn = log(3/2*Rgas*1) + alpha_i*X_i
        rxn_l = self.rxn_l
        for ind,rxn in enumerate(rxn_l):
            # self._param_d['logGth0_rxn'+str(ind)] = -
            self._param_d['alpha_0_rxn'+str(ind)]   = 0.0
            self._param_d['dalpha_t_rxn'+str(ind)]   = 0.0
            self._param_d['dalpha_T_rxn'+str(ind)]   = 0.0
            self._param_d['dalpha_H2O_rxn'+str(ind)] = 0.0

        pass

    @property
    def param_d(self):
        return copy.copy(self._param_d)

    @property
    def rxn_key(self):
        return pd.Series(self._rxn_eqn_l)

    @property
    def phs_key(self):
        return pd.Series(self._phasesym_l)

    @property
    def param_d(self):
        return copy.copy(self._param_d)
    #########

    def get_param_names(self, typ='all'):
        # typ = ['all','free','rxn','phs','rxnadj','rxnall']
        param_names_a = np.array(list(self._param_d.keys()))

        if typ is not None:
            if typ == 'all':
                param_names_typ_a = param_names_a
            elif typ == 'free':
                param_names_typ_a = self._free_param_a
            elif typ == 'rxn':
                param_names_typ_a = \
                    param_names_a[np.char.rfind(param_names_a,'_rxn')>=0]
            elif typ == 'rxnadj':
                param_names_typ_a = \
                    param_names_a[np.char.startswith(param_names_a,'dalpha')]
            elif typ == 'rxnall':
                param_names_typ_a = \
                    param_names_a[np.char.rfind(param_names_a,'_rxn_all')>=0]
            elif typ == 'phs':
                param_names_typ_a = \
                    param_names_a[np.char.rfind(param_names_a,'_phs')>=0]
            else:
                assert False, typ+' is not a valid param typ for get_param_names.'

        # Finally, sort params into sensible order
        # msk_phs_a = np.char.find(param_names_typ_a,'_phs')>=0
        # msk_rxn_all_a = np.char.find(param_names_typ_a,'_rxn_all')>=0
        # msk_rxn_a = np.char.find(param_names_typ_a,'_rxn')>=0
        msk_phs_a = np.array([istr.find('_phs')>=0 for istr in param_names_typ_a])
        msk_rxn_all_a = np.array([istr.find('_rxn_all')>=0 for istr in param_names_typ_a])
        msk_rxn_a = np.array([istr.find('_rxn')>=0 for istr in param_names_typ_a])

        msk_rxn_a = msk_rxn_a*~msk_rxn_all_a
        msk_other_a = ~np.any((msk_phs_a,msk_rxn_all_a,msk_rxn_a),axis=0)

        param_names_sort_a = []
        if(np.any(msk_phs_a)):
            param_names_sort_a.extend(np.sort(param_names_typ_a[msk_phs_a]))
        if(np.any(msk_rxn_all_a)):
            param_names_sort_a.extend(np.sort(param_names_typ_a[msk_rxn_all_a]))
        if(np.any(msk_rxn_a)):
            param_names_sort_a.extend(np.sort(param_names_typ_a[msk_rxn_a]))
        if(np.any(msk_other_a)):
            param_names_sort_a.extend(np.sort(param_names_typ_a[msk_other_a]))

        param_names_sort_a = np.array(param_names_sort_a)
        return param_names_sort_a

    def get_param_values(self, param_key_a, typ=None):
        if typ is not None:
            param_key_a = self.get_param_names(typ=typ)

        param_a = []
        for key in param_key_a:
            param_a.append(self._param_d[key])

        return np.array(param_a)

    def get_param_errors(self, param_key_a, typ=None):
        if typ is not None:
            param_key_a = self.get_param_names(typ=typ)

        error_a = []
        for key in param_key_a:
            error_a.append(self._param_error_d[key])

        return np.array(error_a)

    def get_param_scl_values(self, param_key_a, typ=None):
        if typ is not None:
            param_key_a = self.get_param_names(typ=typ)

        param_scl_a = []
        for key in param_key_a:
            param_scl_a.append(self.param_scl_d[key])

        return np.array(param_scl_a)

    def get_param_table(self, param_nm_a=None, typ='all'):
        if param_nm_a is None:
            param_nm_a = self.get_param_names(typ=typ)

        param_val_a = self.get_param_values(param_nm_a)
        param_scl_a = self.get_param_scl_values(param_nm_a)
        scaled_param_a = self.unscale_params(param_val_a, param_nm_a)
        err_a = self.get_param_errors(param_nm_a)

        param_tbl_d = {'name':param_nm_a,
                       'value':param_val_a,
                       'error':err_a,
                       'scale':param_scl_a,
                       'scaled value':scaled_param_a}
        param_tbl_df = pd.DataFrame(param_tbl_d,columns=['name','value','error',
                                                         'scale','scaled value'])
        return param_tbl_df

    def set_param_values(self, param_val_a, param_key_a ):
        # print(param_val_a)
        for key, val in zip(param_key_a, param_val_a):
            self._param_d[key] = val
            if key.rfind('phs') >=0:
                self.set_phaseobj_param(val,key)

        pass

    def set_phaseobj_param(self, param_val, param_key):
        phsid = 'phs'
        loc = param_key.rfind(phsid)
        phs_ind = int(param_key[loc+len(phsid):])
        iphs = self.phase_l[phs_ind]
        # print(param_key)
        # print(param_val)

        if param_key.startswith('S0'):
            iphs.set_param_values(param_names=['S'],param_values=[param_val])
            # print(iphs.get_param_values(param_names=['S']))

        elif param_key.startswith('V0'):
            iphs.set_param_values(param_names=['V'],param_values=[param_val])
            # print(iphs.get_param_values(param_names=['V']))

        elif param_key.startswith('dH0'):
            initval= iphs.get_param_values(param_names=['delta H'])
            initgibbs = iphs.get_gibbs_energy(300,1)

            iphs.set_param_values(param_names=['delta H'],param_values=[param_val])
            setval = iphs.get_param_values(param_names=['delta H'])

            setgibbs = iphs.get_gibbs_energy(300,1)
            # print(initval,setval,setval-initval)
            # print('%%%%%')
            # print(initgibbs,setgibbs,setgibbs-initgibbs)
            # print('============')

        else:
            assert False, param_key+' is not a valid phase parameter name.'


        pass

    def add_free_params(self, new_free_param_a, typ=None):
        if typ is not None:
            new_free_param_a = self.get_param_names(typ=typ)

        free_param_a = np.hstack((self._free_param_a,new_free_param_a))
        self._free_param_a = np.unique(free_param_a)
        pass

    def del_free_params(self, fix_param_a, typ=None ):
        if typ is not None:
            fix_param_a = self.get_param_names(typ=typ)

        free_param_a = self._free_param_a
        self._free_param_a = np.setdiff1d(free_param_a, fix_param_a)
        pass

    #########

    def extract_trust_data(self):
        """
        Convert units
        """
        data_df = self.data_df

        trust_msk = data_df['Trust']=='Yes'

        # print(data_df[trust_msk])
        # test_df = data_df[trust_msk].set_index('Rxn')
        # print(test_df)

        # NEED TO CAST as FLOAT!!
        num_dat_df = data_df[['P','P_err','T','T_err',
                              'equil_time']][trust_msk].astype(np.float)

        P_a = 1e3*num_dat_df['P'].values
        Perr_a = 1e3*num_dat_df['P_err'].values
        T_a = 273.15+num_dat_df['T'].values
        Terr_a = num_dat_df['T_err'].values
        # time_a = np.log10(num_dat_df['equil_time'].values) # NOTE: use log time
        time_a = num_dat_df['equil_time'].values # NOTE: use log time

        pubid_a = data_df['PubID'][trust_msk].values
        run_num_a = data_df['Run Num.'][trust_msk].values

        rxn_dir_a = data_df['rxn_dir'][trust_msk].values
        rxn_a = data_df['rxn_studied'][trust_msk].values

        # NOTE: need to fix Water flag to flux_amt
        water_a = data_df['Water'][trust_msk].values
        water_flag_a = 1.0*np.ones(water_a.shape)
        water_flag_a[water_a=='dry'] = 0.0
        water_flag_a[water_a=='trace'] = 0.0

        data_trust_d = {}
        data_trust_d['pubid'] = pubid_a
        data_trust_d['run_num'] = run_num_a
        data_trust_d['P'] = P_a
        data_trust_d['Perr'] = Perr_a
        data_trust_d['T'] = T_a
        data_trust_d['Terr'] = Terr_a
        data_trust_d['time'] = time_a
        data_trust_d['rxn_dir'] = rxn_dir_a
        data_trust_d['rxn'] = rxn_a
        data_trust_d['water'] = water_flag_a

        return data_trust_d

    def fit_model(self, free_params=None, full_output=False,
                  method='Nelder-Mead'):
        if free_params is None:
            free_params = self.get_param_names(typ='free')

        # if not np.any([startswith(fix_param,'dH0') for fix_param in fix_param_a]):
        #     assert False, 'User MUST fix some of the Std. State enthalpy params, dH0_X, relative to the elements.'

        # Extract only trustworthy data
        self._dat_trust_d = self.extract_trust_data()

        param0_unscale_a = self.get_param_values(free_params)
        param0_a = self.unscale_params(param0_unscale_a, free_params)
        param0_tbl = self.get_param_table(param_nm_a=free_params)


        # Precalculate approx Gibbs energy uncertainties
        sigG_trust_a = self.propagate_data_errors(param0_unscale_a,
                                                  free_params=free_params)
        self._sigG_trust_a = sigG_trust_a

        model_cost0_d = self.eval_model_costfun(param0_unscale_a,
                                                free_params=free_params,
                                                full_output=True)

        costfun = lambda param_a: self.eval_model_costfun_scl(param_a,
                                                              free_params=free_params)

        lnprob_f = lambda param_a: -costfun(param_a)

        # costfun = self.eval_model_costfun_scl( param_free_vals_scl_a,
        #                                       free_params=free_params )

        # method = 'Nelder-Mead'
        # method = 'BFGS'
        result = optim.minimize(costfun,param0_a,method=method,
                                options={'disp':True,'maxiter':1e4})

        def shift_one_param(shift,ind,mu_a=result.x,costfun=costfun):
            param_a = np.copy(mu_a)
            param_a[ind] += shift
            return costfun(param_a)

        # Create re-scaled-shifted function for hessian
        mu_a = result.x
        cost0 = costfun(mu_a)
        delx_param_scl = np.zeros(mu_a.shape)
        dcost_target=1
        for ind,param in enumerate(mu_a):
            del0 = 1e-2
            idelcostfun = lambda dx, ind=ind,target=dcost_target: \
                shift_one_param(dx,ind)-cost0-dcost_target
            delx = optim.fsolve(idelcostfun,del0)
            delx_param_scl[ind] = np.abs(delx)


        norm_costfun = lambda dx_a, shift_scl=delx_param_scl,\
            mu_a=mu_a,costfun=costfun: costfun(dx_a*shift_scl+mu_a)


        curv_scl_a = delx_param_scl*self.get_param_scl_values(free_params)
        scl_mat_a = core.make_scale_matrix(curv_scl_a)

        Hnorm_fun = nd.Hessian(norm_costfun,step = 1e-2)
        Hnorm_a = Hnorm_fun(np.zeros(mu_a.shape))

        covnorm_a = np.linalg.pinv(Hnorm_a)

        cov_a = covnorm_a*scl_mat_a

        try:
            err_a = np.sqrt(np.diag(cov_a))
            # print(cov_a)
            err_scl_a = core.make_scale_matrix(err_a)
            corr_a = cov_a/err_scl_a
        except:
            err_a = None
            corr_a = None

        # MCMC
        # ndim = len(free_params)
        # nwalkers = 10*ndim
        # walker_pos0_a = [result['x'] + 1e-4*np.random.randn(ndim)
        #                  for i in range(nwalkers)]
        # sampler = emcee.EnsembleSampler(nwalkers, ndim,lnprob_f)
        # sampler.run_mcmc(walker_pos0_a,10)


        # from IPython import embed; embed(); import ipdb; ipdb.set_trace()

        paramf_unscale_a = self.scale_params(result.x, free_params)
        self.set_param_values(paramf_unscale_a,free_params)
        self._fit_result_d = result

        if full_output:

            output_d = {}
            output_d['err_a'] = err_a
            output_d['corr_a'] = corr_a

            self._mu_a = paramf_unscale_a
            self._err_a = err_a
            self._corr_a = corr_a

            for key in self._param_error_d:
                self._param_error_d[key] = np.nan

            for key,err_val in zip(free_params,err_a):
                self._param_error_d[key] = err_val

            model_cost_d = self.eval_model_costfun(paramf_unscale_a,
                                                   free_params=free_params,
                                                   full_output=True)
            param_tbl = self.get_param_table(param_nm_a=free_params)
            param_all_tbl = self.get_param_table(typ='all')

            output_d['free_params'] = free_params

            output_d['costval0'] = model_cost0_d['cost_val']
            output_d['costdata0_df'] = model_cost0_d['cost_data_df']
            output_d['param0_tbl'] = param0_tbl

            output_d['costval'] = model_cost_d['cost_val']
            output_d['costdata_df'] = model_cost_d['cost_data_df']
            output_d['prior_df'] = model_cost_d['prior_df']
            output_d['param_tbl'] = param_tbl
            output_d['param_all_tbl'] = param_all_tbl

            output_d['result'] = result

            # output_d['param_d'] = copy.copy(self._param_d)
            # output_d['param0_a'] = param0_a
            # output_d['paramf_a'] = result.x
            # output_d['param0_unscl_a'] = param0_unscale_a
            # output_d['paramf_unscl_a'] = paramf_unscale_a

            return output_d

        pass

    def posterior_draw(self, Ndraw=1):
        param_free = self.get_param_names(typ='free')
        mu_a = self._mu_a
        err_a = self._err_a
        corr_a = self._corr_a
        err_scl_mat_a = core.make_scale_matrix(err_a)

        cov_a = err_scl_mat_a*corr_a

        pdraw_a = np.squeeze(random.multivariate_normal(mu_a, cov_a, Ndraw))
        return pdraw_a

    def posterior_rxn_bound(self, Tlims, conf_lvl=0.68, Ndraw=100,
                            Nsamp=30, sampfac=3, convert_units=True):
        Nrxn = len(self.rxn_key)
        Nsamptot = np.round(Nsamp*sampfac)

        # T_bound_draw_a = np.zeros((Nrxn,Ndraw,Nsamptot))
        P_bound_draw_a = np.zeros((Nrxn,Ndraw,Nsamptot))
        T_a = np.linspace(Tlims[0],Tlims[1],Nsamptot)
        free_param_nm_a= self.get_param_names(typ='free')
        PT_triple_draw_a = np.zeros((2,Ndraw))

        for ind in range(Ndraw):
            pdraw_a = self.posterior_draw()
            self.set_param_values(pdraw_a,free_param_nm_a)

            for irxn,(rxn_eqn,rxn_obj) in enumerate(zip(self.rxn_key,self.rxn_l)):
                iTP_bound_a = rxn_obj.trace_rxn_bound(Tlims=Tlims,Nsamp=Nsamp)
                fun = interp.interp1d(iTP_bound_a[0],iTP_bound_a[1])
                Pbnd_a = fun(T_a)

                # T_bound_draw_a[irxn,ind,:] = T_a
                P_bound_draw_a[irxn,ind,:] = Pbnd_a

            T_tp,P_tp = self.rxn_l[0].get_simultaneous_rxn_cond(self.rxn_l[1])
            PT_triple_draw_a[0,ind] = T_tp
            PT_triple_draw_a[1,ind] = P_tp


        if convert_units:
            T_a -= 273.15
            P_bound_draw_a/=1e3
            PT_triple_draw_a[0] -= 273.15
            PT_triple_draw_a[1] /= 1e3

        posterior_rxn_d = {}
        posterior_rxn_d['T_a'] = T_a
        posterior_rxn_d['P_bound_draw_a'] = P_bound_draw_a
        posterior_rxn_d['PT_triple_draw_a'] = PT_triple_draw_a

        self.calc_rxn_bound_conf_lvl(posterior_rxn_d,conf_lvl=conf_lvl)

        return posterior_rxn_d

    def calc_rxn_bound_conf_lvl(self,posterior_rxn_d, conf_lvl=0.68):
        T_a = posterior_rxn_d['T_a']
        P_bound_draw_a = posterior_rxn_d['P_bound_draw_a']
        PT_triple_draw_a = posterior_rxn_d['PT_triple_draw_a']

        rxn_conf_bnd_a=np.percentile(P_bound_draw_a,
                                     [50-0.5*100*conf_lvl,50+0.5*100*conf_lvl],
                                     axis=1)
        PT_triple_mean_a = np.mean(PT_triple_draw_a,axis=1)
        PT_triple_cov_a = np.cov(PT_triple_draw_a)

        posterior_rxn_d['rxn_conf_bnd_a'] = rxn_conf_bnd_a
        posterior_rxn_d['PT_triple_mean_a'] = PT_triple_mean_a
        posterior_rxn_d['PT_triple_cov_a'] = PT_triple_cov_a
        pass

    def scale_params(self, param_vals_scl_a, param_names_a):
        param_scl_a = self.get_param_scl_values(param_names_a)
        param_vals_a = param_vals_scl_a*param_scl_a
        return param_vals_a

    def unscale_params(self, param_vals_a, param_names_a):
        param_scl_a = self.get_param_scl_values(param_names_a)
        param_vals_scl_a = param_vals_a/param_scl_a
        return param_vals_scl_a

    def eval_model_costfun_scl(self, param_free_vals_scl_a, free_params=None):
        if free_params is None:
            free_params = self.get_param_names(typ='free')

        # param_free_scl_a = thermoDB_mod.get_param_scl_values(free_params)
        # param_free_vals_a = param_free_vals_scl_a*param_free_scl_a

        param_free_vals_a = self.scale_params(param_free_vals_scl_a,free_params)

        cost_fun = self.eval_model_costfun(param_free_vals_a,
                                           free_params=free_params)

        return cost_fun

    def propagate_data_errors(self, param_free_vals_a, free_params=None):

        if free_params is None:
            free_params = self.get_param_names(typ='free')

        if param_free_vals_a is None:
            param_free_vals_a = self.get_param_values(free_params)



        # Extract only trustworthy data
        # if self._dat_trust_d is None:
        #     self._dat_trust_d = self.extract_trust_data()

        self.set_param_values(param_free_vals_a,free_params)
        # self.set_free_param_values(param_free_vals_a)

        data_df = self.data_df
        rxn_eqn_l = self._rxn_eqn_l
        rxn_l = self.rxn_l
        param_d = self._param_d

        dat_d = self._dat_trust_d
        # dat_d = self.extract_trust_data()

        # print(P_a)
        # print(P_a.reset_index())

        Ndat = dat_d['P'].size
        dVrxn_a = np.zeros(Ndat)
        dSrxn_a = np.zeros(Ndat)

        for ind,(rxn_eqn, rxn_obj) in enumerate(zip(rxn_eqn_l,rxn_l)):
            msk_rxn = dat_d['rxn']==rxn_eqn
            idVrxn_a =  rxn_obj.get_rxn_volume(dat_d['T'][msk_rxn],
                                               dat_d['P'][msk_rxn],
                                               peratom=True )
            idSrxn_a =  rxn_obj.get_rxn_entropy(dat_d['T'][msk_rxn],
                                                dat_d['P'][msk_rxn],
                                                peratom=True )

            dVrxn_a[msk_rxn] = idVrxn_a
            dSrxn_a[msk_rxn] = idSrxn_a

        sigG_trust_a = np.sqrt((dat_d['Perr']*dVrxn_a)**2+(dat_d['Terr']*dSrxn_a)**2)
        return sigG_trust_a

    def eval_model_costfun(self, param_free_vals_a, free_params=None,
                           full_output=False):

        if free_params is None:
            free_params = self.get_param_names(typ='free')


        # Extract only trustworthy data
        # if self._dat_trust_d is None:
        #     self._dat_trust_d = self.extract_trust_data()

        self.set_param_values(param_free_vals_a,free_params)
        # self.set_free_param_values(param_free_vals_a)

        # Try to use precalculated approx Gibbs energy uncertainties
        sigG_a = self._sigG_trust_a
        if sigG_a is None:
            sigG_a = self.propagate_data_errors(param_free_vals_a,
                                                free_params=free_params)
            self._sigG_trust_a = sigG_a

        Gthresh_scl = self.Gthresh_scl
        data_df = self.data_df
        rxn_eqn_l = self._rxn_eqn_l
        rxn_l = self.rxn_l
        param_d = self._param_d

        dat_d = self._dat_trust_d
        # dat_d = self.extract_trust_data()

        # print(P_a)
        # print(P_a.reset_index())

        Ndat = dat_d['P'].size
        dGrxn_a = np.zeros(Ndat)
        logGth_a = np.zeros(Ndat)

        # dVrxn_a = np.zeros(Ndat)
        # dSrxn_a = np.zeros(Ndat)

        # get universal reaction parameters
        alpha_t_all   = param_d['alpha_t_rxn_all']
        alpha_T_all   = param_d['alpha_T_rxn_all']
        alpha_H2O_all = param_d['alpha_H2O_rxn_all']

        for ind,(rxn_eqn, rxn_obj) in enumerate(zip(rxn_eqn_l,rxn_l)):
            msk_rxn = dat_d['rxn']==rxn_eqn
            idGrxn_a = rxn_obj.get_rxn_gibbs_energy(dat_d['T'][msk_rxn],
                                                    dat_d['P'][msk_rxn],
                                                    peratom=True )
            # idVrxn_a =  rxn_obj.get_rxn_volume(dat_d['T'][msk_rxn],
            #                                    dat_d['P'][msk_rxn],
            #                                    peratom=True )
            # idSrxn_a =  rxn_obj.get_rxn_entropy(dat_d['T'][msk_rxn],
            #                                     dat_d['P'][msk_rxn],
            #                                     peratom=True )


            # get reaction-specific parameters
            alpha_0_rxn = param_d['alpha_0_rxn'+str(ind)]
            dalpha_t_rxn = param_d['dalpha_t_rxn'+str(ind)]
            dalpha_T_rxn = param_d['dalpha_T_rxn'+str(ind)]
            dalpha_H2O_rxn = param_d['dalpha_H2O_rxn'+str(ind)]

            logGth_a[msk_rxn] = alpha_0_rxn \
                + (alpha_t_all+dalpha_t_rxn)*dat_d['time'][msk_rxn] \
                + (alpha_T_all+dalpha_T_rxn)*dat_d['T'][msk_rxn] \
                + (alpha_H2O_all+dalpha_H2O_rxn)*dat_d['water'][msk_rxn]

            dGrxn_a[msk_rxn] = idGrxn_a
            # dVrxn_a[msk_rxn] = idVrxn_a
            # dSrxn_a[msk_rxn] = idSrxn_a

        Gth_a = Gthresh_scl*np.exp(logGth_a)
        # sigG_a = np.sqrt((dat_d['Perr']*dVrxn_a)**2+(dat_d['Terr']*dSrxn_a)**2)

        loglk_a = np.zeros(Gth_a.size)

        for irxndir in self.rxn_dir_opt:
            msk_rxndir = dat_d['rxn_dir']==irxndir
            iNobs = np.sum(msk_rxndir==True)
            if iNobs == 0:
                continue

            if irxndir=='FWD':
                iloglk_a = self.logprob_fwd(dGrxn_a[msk_rxndir],
                                              Gth_a[msk_rxndir],
                                              sigG_a[msk_rxndir],
                                              rxn_trans_typ=self.rxn_trans_typ)
            elif irxndir=='REV':
                iloglk_a = self.logprob_rev(dGrxn_a[msk_rxndir],
                                              Gth_a[msk_rxndir],
                                              sigG_a[msk_rxndir],
                                              rxn_trans_typ=self.rxn_trans_typ)
            elif irxndir=='FWD?':
                iloglk_a = self.logprob_fwdq(dGrxn_a[msk_rxndir],
                                               Gth_a[msk_rxndir],
                                               sigG_a[msk_rxndir],
                                               rxn_trans_typ=self.rxn_trans_typ)
            if irxndir=='REV?':
                iloglk_a = self.logprob_revq(dGrxn_a[msk_rxndir],
                                               Gth_a[msk_rxndir],
                                               sigG_a[msk_rxndir],
                                               rxn_trans_typ=self.rxn_trans_typ)
            elif irxndir=='NC':
                iloglk_a = self.logprob_nc(dGrxn_a[msk_rxndir],
                                             Gth_a[msk_rxndir],
                                             sigG_a[msk_rxndir],
                                             rxn_trans_typ=self.rxn_trans_typ)

            loglk_a[msk_rxndir] = iloglk_a


        msk_zeroprob_a = np.isinf(loglk_a) & (loglk_a<0)

        logprior_a = self.eval_log_prior()
        # logprob_a[msk_zeroprob_a] = -160
        cost_val_a = np.hstack((-loglk_a,-logprior_a))
        cost_val = np.sum(cost_val_a)

        if full_output:
            model_cost_d = {}
            model_cost_d['cost_val'] = cost_val

            logprior_a,prior_df = self.eval_log_prior(full_output=True)

            cost_data_df = pd.DataFrame(dat_d)
            cost_data_df['zeroprob'] = pd.Series(msk_zeroprob_a)
            cost_data_df['log_lk'] = pd.Series(loglk_a)
            # cost_data_df['PubID'] = pubid_a
            # cost_data_df['run_num'] = run_num_a
            # cost_data_df['cost_fun'] = cost_fun_a
            # cost_data_df['rxn_dir'] = rxn_dir_a
            cost_data_df['Gth'] = pd.Series(Gth_a)
            cost_data_df['sigG'] = pd.Series(sigG_a)
            cost_data_df['dGrxn'] = pd.Series(dGrxn_a)
            cost_data_df['relG'] = pd.Series(dGrxn_a/sigG_a)

            model_cost_d['cost_data_df'] = cost_data_df
            model_cost_d['log_prior'] = pd.Series(logprior_a)
            model_cost_d['prior_df'] = prior_df
            return model_cost_d
        else:
            return cost_val

    def eval_log_prior(self,typ='studentt',dof=5,full_output=False):
        paramnm_s = self.exp_prior_df['Param']
        trust_s = self.exp_prior_df['Trust']
        val_data_s = self.exp_prior_df['Data']
        err_s = self.exp_prior_df['Error']

        log_prior_a = np.zeros(val_data_s.shape[0])
        val_model_a = np.zeros(val_data_s.shape[0])
        resid_a = np.zeros(val_data_s.shape[0])

        for ind,(paramnm,trust,val_data,err) in \
                enumerate(zip(paramnm_s,trust_s,val_data_s,err_s)):

            val_mod = self.param_d[paramnm]
            val_model_a[ind] = val_mod
            x = (val_mod-val_data)/err
            resid_a[ind] = x
            if trust=='Yes':
                log_prior_a[ind] = self.logprior_fun(x)
            else:
                log_prior_a[ind] = 0.0


        if full_output:
            # Get new dataframe by copything columns
            prior_df = self.exp_prior_df[['Param','Abbrev']].copy()
            prior_df['Data'] = self.exp_prior_df['Data']
            prior_df['Error'] = self.exp_prior_df['Error']
            prior_df['Model'] =  pd.Series(val_model_a)
            prior_df['Resid'] =  pd.Series(resid_a)
            prior_df['Trust'] = self.exp_prior_df['Trust']
            prior_df['log_prior'] =  pd.Series(log_prior_a)
            return log_prior_a, prior_df
        else:
            return log_prior_a
#===================================================
class PhaseStabilityData:
    def __init__(self, exp_data, phase_wt_comp, phase_symbol_key, modelDB,
                    T_units='K', P_units='bars'):
        # self._prune_missing_phases(exp_dataphase_wt_comp, modelDB)
        self._init_exps(exp_data, phase_wt_comp, phase_symbol_key, modelDB,
                            T_units, P_units)

    def _init_exps(self, exp_data, phs_wt_comp, phase_symbol_key, modelDB,
                    T_units, P_units):
        phase_symbol_key = chem.LEPR_phase_symbols #Jenna added this; to get
        # this function to work properly, it seems like we should take output
        # the phase_symbol_key as an input, and create a comprehensive lists
        # of these lepr phase symbols?
        phase_stability_exps=[]
        exp_index_invalid = []

        for exp_idx in exp_data.index:
            iP = exp_data.loc[exp_idx, 'P']
            iT = exp_data.loc[exp_idx, 'T']

            phase_symbols = []
            phase_names = []
            phase_wt_oxide_comp = []

            for key in phs_wt_comp:
                iphs_wt_comp = phs_wt_comp[key]
                if exp_idx in iphs_wt_comp.index:
                    phase_names.append(key)
                    phase_symbols.append(phase_symbol_key[key])
                    phase_wt_oxide_comp.append(iphs_wt_comp.loc[exp_idx].values)

            phase_wt_oxide_comp = np.array(phase_wt_oxide_comp)

            try:
                phase_stability_exp = PhaseStabilityExp(
                    iP, iT, 0, phase_symbols, phase_wt_oxide_comp, modelDB,
                    T_units, P_units)
                phase_stability_exps.append(phase_stability_exp)

            except:
                exp_index_invalid.append(exp_idx)

        self._exp_index_invalid = exp_index_invalid
        self._phase_stability_exps = phase_stability_exps

        #phase_stability_info = OrderedDict()
        #phase_stability_info['phase_symbols'] = self._phase_symbols = phase_symbols
        #phase_stability_info['phases'] = phases
        #phase_stability_info['phase_num'] = phase_num
        #phase_stability_info['phase_wt_oxide_comp'] = phase_wt_oxide_comp
        #phase_stability_info['phase_mol_oxide_comp'] = phase_mol_oxide_comp
        #phase_stability_info['phase_mol_endmem_comp'] = phase_mol_endmem_comp

    def calc_equil_rxn_affinities(self):
        phase_stability_exps = self._phase_stability_exps

        affinities = []

        for iexp in phase_stability_exps:
            iaffinity = iexp.calc_equil_rxn_affinities()
            affinities.append(iaffinity)

        return affinities

    # def _prune_missing_phases(exp_data, phase_wt_comp, modelDB):
    #     self._modelDB = modelDB
    #     prune_phases = []
    #
    #     for phase_sym in phase_wt_comp:
    #         try:
    #             phase = modelDB.get_phase(phase_sym)
    #
    #         except:
    #             prune_phases.append(phase_sym)

    #     for idx in exp_data.index:




#===================================================
class PhaseStabilityExp:
    """
    """

    def __init__(self, T, P, wt_oxide_comp,
                 phase_symbols, phase_wt_oxide_comp, modelDB,
                 T_err=None, P_err=None, wt_oxide_comp_err=None,
                 phase_wt_oxide_comp_err=None,
                 fO2_buffer=None, fO2_offset=0, fO2_err=None,
                 T_units='K', P_units='bars'):

        self._init_exp_cond(P, T, P_err, T_err,
                            fO2_buffer, fO2_offset, fO2_err, T_units, P_units)
        self._init_phase_info(phase_symbols, phase_wt_oxide_comp, modelDB)
        self._init_equil_rxns()
            # calculate rxns with svd
        #self._init_absent_rxns()

    def _init_exp_cond(self, T, P, T_err, P_err,
                       fO2_buffer, fO2_offset, fO2_err, T_units, P_units):
        if T_units == 'C':
            T = T + 273
        if P_units == 'GPa':
            P = P * 10000
            P_err = P_err * 10000

        self._P = P
        self._T = T

        self._P_err = P_err
        self._T_err = T_err

        self._fO2_buffer = fO2_buffer
        self._fO2_offset = fO2_offset
        self._fO2_err = fO2_err

    def _init_phase_info(self, phase_symbols, phase_wt_oxide_comp, modelDB):
        phase_wt_oxide_comp = np.array(phase_wt_oxide_comp)

        phase_num = len(phase_symbols)
        oxide_num = len(chem.OXIDE_ORDER)
        assert phase_wt_oxide_comp.shape[0]==phase_num, (
            'phase_wt_oxide_comp must define compositions '
            'for every phase in phase_symbols.'
        )
        assert phase_wt_oxide_comp.shape[1]==oxide_num, (
            'phase_wt_oxide_comp must define composition for every oxide '
            'in standard order.'
        )

        phases = [modelDB.get_phase(phs_sym) for phs_sym in phase_symbols]
        phase_mol_oxide_comp = chem.wt_to_mol_oxide(phase_wt_oxide_comp)

        phase_mol_endmem_comp = {}

        for phs_sym, phase, mol_oxide_comp in zip(
            phase_symbols, phases, phase_mol_oxide_comp):

            if phs_sym not in modelDB.phase_obj['pure']:
                endmem_comp = phase.calc_endmember_comp(
                    mol_oxide_comp, method='intrinsic', output_residual=False,
                    normalize=True)
            else:
                continue

            phase_mol_endmem_comp[phs_sym] = np.array(endmem_comp)


        self._modelDB = modelDB
        self._phase_symbols = phase_symbols
        self._phases = phases
        self._phase_num = phase_num
        self._phase_wt_oxide_comp = phase_wt_oxide_comp
        self._phase_mol_oxide_comp = phase_mol_oxide_comp
        self._phase_mol_endmem_comp = phase_mol_endmem_comp

    def _init_equil_rxns(self):
        phases = self._phases
        phase_symbols = self._phase_symbols
        modelDB = self._modelDB

        rxn_svd_props = chem.calc_reaction_svd(phase_symbols, TOLsvd=1e-4, modelDB=modelDB)
        equil_rxn_coefs = rxn_svd_props['rxn_svd']
        Nbasis=len(equil_rxn_coefs)
        rxn_endmember_name = rxn_svd_props['all_endmember_name']
        rxn_phase_symbols = rxn_svd_props['all_phase_symbol']
        endmember_ids = rxn_svd_props['all_endmember_id']

        equil_rxns = []
        for irxn_coefs in equil_rxn_coefs:
            irxn = modelDB.get_rxn(rxn_phase_symbols, endmember_ids, irxn_coefs, coefs_per_atom=True)

            equil_rxns.append(irxn)

        equil_rxn_num = len(equil_rxns)

        all_endmem_names = rxn_endmember_name

        self._all_endmem_names = all_endmem_names
        self._equil_rxn_num = equil_rxn_num
        self._equil_rxns = equil_rxns
        self._equil_rxn_coefs = equil_rxn_coefs

    def calc_equil_rxn_affinities(self):
        phases = self._phases
        equil_rxns = self._equil_rxns
        equil_rxn_num = self._equil_rxn_num
        T = self._T
        P = self._P
        phase_mol_endmem_comp = self._phase_mol_endmem_comp
        # print(T)
        # print(P)
        # print(phase_mol_endmem_comp)

        rxn_affinities = np.zeros(equil_rxn_num)
        for ind, rxn in enumerate(equil_rxns):
            rxn_affinities[ind] = rxn.affinity(T, P, mols=phase_mol_endmem_comp)

        return rxn_affinities

    def calc_absent_rxn_affinities(self):
        raise NotImplimented()

#===================================================
class RxnData:
    """
    TODO:
        - filter trusted data
        - get_data() - return only trusted info?
    """
    RESULT_OPTS = ['+','-','=','+?','-?']

    ERROR_DEFAULT = {'P':5, 'T':10, 'mol':1}
    ERROR_TYPE = {'P':'%', 'T':'absolute', 'mol':'%'}

    def __init__(self, input_data_sheets, error_default=None):
        self._init_data_tables()
        self.load_data(input_data_sheets)

    def _init_data_tables(self):
        self.reference = self._init_reference_data()
        self.setup = self._init_setup_data()
        self.conditions = self._init_conditions_data()
        self.rxn = self._init_rxn_data()
        self.comp = self._init_comp_data()
        pass

    def _init_reference_data(self):
        reference = pd.DataFrame(columns=[
            'pub_id', 'authors', 'date', 'title'])
        return reference

    def _init_setup_data(self):
        setup = pd.DataFrame(columns=[
            'pub_id','run_id',
            'total_mass', 'contaminant_phases', 'contact_geom',
            'container_aspect_ratio', 'grain_size', 'single_xtal_size',
            'other_phase_frac', 'flux_type', 'flux_amt',
            'init_reac_present', 'init_prod_present'])
        return setup

    def _init_conditions_data(self):
        conditions = pd.DataFrame(columns=[
            'P', 'P_err', 'T', 'T_err', 'equil_time', 'trust_conditions'])
        return conditions

    def _init_rxn_data(self):
        rxn = pd.DataFrame(columns=[
            'rxn_id', 'rxn_studied', 'init_rxn_progress', 'results', 'rxn_dir'])
        return rxn

    def _init_comp_data(self):
        comp = OrderedDict()
        return comp

    def load_data(self, input_data):
        self.input_data = input_data
        self._load_reference_data(input_data)
        self._load_setup_data(input_data)
        self._load_conditions_data(input_data)
        self._load_rxn_data(input_data)
        pass

    def _load_reference_data(self, input_data):
        reference = self.reference
        reference_cols = self.reference.columns
        ref_data, ref_units = self._read_data_sheet(
            input_data['reference'][reference_cols])
        self.reference = reference.append(ref_data)[reference_cols]
        pass

    def _load_setup_data(self, input_data):
        setup = self.setup
        setup_cols = setup.columns
        exp_conditions_cols = ['pub_id','run_id',
            'container_aspect_ratio']

        rxn_cols = ['contact_geom', 'total_mass',
                    'other_phase_frac', 'contaminant_phases',
                    'grain_size', 'single_xtal_size', 'flux_type', 'flux_amt',
                    'init_reac_present', 'init_prod_present']

        isetup = pd.concat([
            input_data['exp_conditions'][exp_conditions_cols],
            input_data['rxn'][rxn_cols]], axis=1)

        isetup_data, isetup_units = self._read_data_sheet(isetup)
        self.setup = setup.append(isetup_data)[setup_cols]
        pass

    def _load_conditions_data(self, input_data):
        conditions = self.conditions
        conditions_cols = conditions.columns
        iconditions_data, iconditions_units = self._read_data_sheet(
            input_data['exp_conditions'][conditions_cols])

        # Convert C to K
        if iconditions_units['T'] == 'C':
            iconditions_data['T'] += 273.15

        if iconditions_units['P'] == 'kbar':
            iconditions_data['P'] *= 1e3
            iconditions_data['P_err'] *= 1e3

        iconditions_data['trust_conditions'].fillna(value='Yes', inplace=True)

        iconditions_data = self._set_default_error(iconditions_data)
        self.conditions = conditions.append(iconditions_data)
        pass

    def _set_default_error(self, conditions_data):
        ERROR_DEFAULT = self.ERROR_DEFAULT
        ERROR_TYPE = self.ERROR_TYPE


        for key in ['T', 'P']:
            default_val = ERROR_DEFAULT[key]
            default_typ = ERROR_TYPE[key]

            # from IPython import embed; embed(); import ipdb; ipdb.set_trace()
            mask = conditions_data[key+'_err'].isnull()
            # print(mask)

            if default_typ=='absolute':
                ierr = default_val
            elif default_typ=='%':
                ival = conditions_data.loc[mask, key]
                ival = np.maximum(ival, 1.0)
                ierr = np.abs(ival*default_val/100)

            conditions_data.loc[mask, key+'_err'] = ierr

        return conditions_data

    def _load_rxn_data(self, input_data):
        rxn = self.rxn
        rxn_cols = rxn.columns
        rxn_input_cols = ['rxn_studied', 'init_rxn_progress', 'results']

        rxn_input_data, rxn_input_units = self._read_data_sheet(
            input_data['rxn'][rxn_input_cols])

        # Redetermine rxn directions by combining ALL rxns reported
        all_rxn_input_data = rxn[rxn_input_cols].append(rxn_input_data)

        rxn_props, phases, rxn_ids = self._init_rxns(all_rxn_input_data)
        rxn_dir = self._infer_rxn_dir(all_rxn_input_data, rxn_props, rxn_ids)

        self.rxn = pd.concat((pd.Series(rxn_ids, name='rxn_id'),
                              all_rxn_input_data,
                              pd.Series(rxn_dir, name='rxn_dir')), axis=1)

        # self.rxn = rxn.append(irxn)[rxn.columns]
        self.rxn_num = len(rxn_props)
        self.rxn_props = rxn_props
        self.phases = phases
        pass

    def _load_comp_data(self, comp, input_data):
        pass

    def _init_phase_comp(self, component):
        phase_comp = pd.DataFrame(columns=[
            component, component+'_err', 'init_'+component,
            'init_'+component+'_err'])
        self.comp[component] = phase_comp
        pass

    def _read_data_sheet(self, sheet):
        units = sheet.loc[0]
        data = sheet.loc[1:].reset_index(drop=True)
        return data, units

    def _init_rxns(self, rxn):
        rxn_eqns = rxn['rxn_studied'].unique()

        rxn_props = []
        all_phases = []
        rxn_ids = np.zeros(rxn.shape[0], dtype=int)
        for ind, irxn_eqn in enumerate(rxn_eqns):
            imask = rxn['rxn_studied']==irxn_eqn
            rxn_ids[imask] = ind

            iphs, icoef = self.extract_rxn_coefs(irxn_eqn)
            all_phases.append(iphs)
            irxn_prop = OrderedDict({'phases': iphs, 'coefs': icoef, 'eqn':irxn_eqn})
            rxn_props.append(irxn_prop)

        phases = np.unique(all_phases)

        return rxn_props, phases, rxn_ids

    def _infer_rxn_dir(self, rxn, rxn_props, rxn_ids):
        results = rxn['results']
        rxn_dir = []
        for iresult, irxn_id in zip(results, rxn_ids):
            iphases, ichanges = self.extract_rxn_results(iresult)
            ichange_dir = np.nan*np.ones(len(ichanges))

            irxn_prop = rxn_props[irxn_id]
            irxn_phases = irxn_prop['phases']
            irxn_coefs = irxn_prop['coefs']

            irxn_dir = np.zeros(len(ichanges))
            irxn_certain = np.tile(False, len(ichanges))
            for ind, (ijphs, ijchange) in enumerate(zip(iphases, ichanges)):
                ijdir = np.nan
                ijcertain = False
                if ijchange == '+':
                    ijdir = +1
                    ijcertain = True
                elif ijchange == '-':
                    ijdir = -1
                    ijcertain = True
                elif ijchange == '+?':
                    ijdir = +1
                    ijcertain = False
                elif ijchange == '-?':
                    ijdir = -1
                    ijcertain = False
                elif (ijchange=='='):
                    ijdir = 0
                    ijcertain = True
                elif (ijchange=='=?'):
                    ijdir = 0
                    ijcertain = False
                else:
                    ijdir = np.nan
                    ijcertain = False

                ijdir *= irxn_coefs[np.where(irxn_phases==ijphs)]
                irxn_dir[ind] = ijdir
                irxn_certain[ind] = ijcertain

            if np.all(irxn_dir==-1):
                if np.any(irxn_certain):
                    rxn_dir.append('REV')
                else:
                    rxn_dir.append('REV?')

            elif np.all(irxn_dir==+1):
                if np.any(irxn_certain):
                    rxn_dir.append('FWD')
                else:
                    rxn_dir.append('FWD?')

            elif np.all(irxn_dir==0):
                rxn_dir.append('NC')

            else:
                rxn_dir.append('BIASED')

        rxn_dir = np.array(rxn_dir)

        return rxn_dir






        # for irxn_props in rxn_props:
        #     irxn_
        #     imask = df['rxn_studied']==irxn_eqn
        #     iresults = pd.DataFrame([pd.Series(df[imask]['results'].str.startswith(iphs),name=iphs)
        #                          for iphs in irxn.phase_symbols]).T
        #
        #     irxn_dir = np.dot(irxn.rxn_coefs,
        #                       np.array([df[imask]['results'].str.startswith(iphs)
        #                                 for iphs in irxn.phase_symbols]))
        #     rxn_dir[imask] = irxn_dir
        #
        # df['rxn_dir'] = rxn_dir
        return None

    @classmethod
    def read_phase_rev_data(cls, filenm, sheetname=None):
        # Read file and concatenate all sheets
        data_d = pd.read_excel(filenm,sheetname=sheetname)
        # try to concatenate multiple sheets (if present)
        try:
            raw_df = pd.concat(data_d,ignore_index=True)

        except:
            raw_df = data_d

        return raw_df

    @classmethod
    # Determine which values are actually bounds
    def detect_bound(self, colnm, df):
        msk_lo = df[colnm].astype(np.object).str.startswith('>').fillna(value=False).astype(bool)
        msk_hi = df[colnm].astype(np.object).str.startswith('<').fillna(value=False).astype(bool)

        # NOTE: It is crucial that bound is a series (not a numpy array), otherwise msk indexing will fail
        bound_ser = pd.Series(np.tile('',msk_lo.size))
        bound_ser[msk_lo] = 'lower'
        bound_ser[msk_hi] = 'upper'


        bound_df = pd.DataFrame()
        bound_df[colnm+'_Bound'] =bound_ser

        bound_df[colnm] = df[colnm].copy()
        bound_df.loc[msk_lo,colnm] = df.loc[msk_lo,colnm].astype(np.object).str.slice(start=1).astype(np.float)
        bound_df.loc[msk_hi,colnm] = df.loc[msk_hi,colnm].astype(np.object).str.slice(start=1).astype(np.float)

        return bound_df

    def filter_phase_rev_data(self, raw_df, mask_phs_l=None):
        P_df = self.detect_bound('P',raw_df)
        T_df = self.detect_bound('T',raw_df)
        time_df = self.detect_bound('equil_time',raw_df)

        rxn_df = pd.DataFrame()

        # rxn_df['RxnPhases'] = raw_df[['rxn_studied']].applymap(get_reaction_str)
        rxn_df['RxnPhases'] = raw_df[['rxn_studied']].applymap(model.Database._get_reaction_phase_str)


        # Set Rxn equation as the most common string representaton in the raw database
        RxnPhases_uniq = rxn_df['RxnPhases'].unique()
        rxn_df['rxn_studied'] = raw_df['rxn_studied']
        for rxn_phs_str in RxnPhases_uniq:
            this_reaction = raw_df['rxn_studied'][rxn_df['RxnPhases']==rxn_phs_str]
            this_reaction = this_reaction.str.strip()
            # Store eqn only as most common variant
            #NOTE iloc crucial to obtain just value (not series object)
            rxn_df.loc[rxn_df['RxnPhases']==rxn_phs_str,'rxn_studied'] = this_reaction.mode().iloc[0]


        rxn_d_l = []
        phs_l = []
        rxn_eqn_uniq = rxn_df['rxn_studied'].unique()
        for rxn_eqn_str in rxn_eqn_uniq:
            rxn_d = model.Database.parse_rxn( rxn_eqn_str )
            #Rewrite Rxn equation using adopted rxn direction
            rxn_df.loc[rxn_df['rxn_studied']==rxn_eqn_str,'rxn_studied'] = rxn_d['rxn_eqn']

            # Remove masked phases
            if mask_phs_l is not None:
                curr_rxn_phs_l = []
                curr_rxn_phs_l.extend(rxn_d['reac_l'])
                curr_rxn_phs_l.extend(rxn_d['prod_l'])
                disallowed_phs_l = np.intersect1d( curr_rxn_phs_l, mask_phs_l )
                if len(disallowed_phs_l)>0:
                    continue

            phs_l.extend( rxn_d['reac_l'] )
            phs_l.extend( rxn_d['prod_l'] )
            rxn_d_l.append( rxn_d )


        # Remove masked phases
        trust_ser = raw_df['trust'].fillna(value='Yes')
        if mask_phs_l is not None:
            # Remove from phase list
            phs_l = np.setdiff1d( phs_l, mask_phs_l )

            # set Trust variable to No for rxns involving masked phase
            for mask_phs in mask_phs_l:
                trust_ser[rxn_df['RxnPhases'].str.contains(mask_phs)] = 'No'



        phs_uniq_l = np.unique( phs_l )

        # Determine Reaction Direction: ['FWD','REV','NC','FWD?','REV?','INV']
        rxn_dir_l = []
        for result, rxn_phs_str in zip(raw_df['results'],rxn_df['RxnPhases']):
            rxn_dir_l.append(model.Database._get_rxn_dir(rxn_phs_str, result))



        # rxn_uniq = rxn_df['Rxn'].unique()
        # print(rxn_uniq)
        rxn_df['rxn_dir'] = pd.Series(rxn_dir_l)
        # rxn_df['results'] = raw_df['results']

        # result = rxn_df[['Rxn']].applymap(get_reaction_phase_str)
        # print(get_reaction_phases(result.loc[0,'Rxn']))

        # print(result)
        #  print(result['Rxn'].unique())
        #  rxn_df = pd.concat((rxn_df,pd.DataFrame(result['Rxn'].tolist(),columns=['Reac_l','Prod_l'])),axis=1)

        # print(result.tolist())
        # print(pd.DataFrame(result,columns=['Reac','Prod']))

        dat_df = pd.concat((raw_df['pub_id'],raw_df['device'],
                            raw_df['run_id'],
                            time_df,raw_df['flux_amt'],
                            P_df,raw_df['P_err'],T_df,raw_df['T_err'],
                            rxn_df,trust_ser),axis=1)

        return dat_df, rxn_d_l, phs_uniq_l

    @classmethod
    def extract_rxn_coefs(cls, rxn_eqn_str):
        eqn_split = re.split('=', rxn_eqn_str)
        assert len(eqn_split) - 1 == 1, (
            'rxn_eqn_str must have exactly one = sign'
            )
        reac_str, prod_str = eqn_split
        reac_str = str.strip(reac_str)
        prod_str = str.strip(prod_str)

        def extract_coefs(eqn_side_str):
            eqn_terms_str = list(filter(None, re.split('[ +=]', eqn_side_str)))

            coefs = []
            phases_str = []

            for iterm_str in eqn_terms_str:
                match = re.match('^[0-9]', iterm_str)
                if not match:
                    icoef = 1.0
                else:
                    icoef = float(match.group())

                coefs.append(icoef)
                phases_str.append(re.split('^[0-9]', iterm_str)[-1])

            coefs = np.array(coefs)

            return phases_str, coefs

        phases_str_reac, coefs_reac = extract_coefs(reac_str)
        phases_str_prod, coefs_prod = extract_coefs(prod_str)

        phases_str = np.hstack((phases_str_reac, phases_str_prod))
        coefs = np.hstack((-coefs_reac, +coefs_prod))

        return phases_str, coefs

    @classmethod
    def extract_rxn_results(cls, result):
        result_terms = [str.strip(iresult) for iresult in str.split(result,';')]

        phases = []
        changes = []

        for iterm in result_terms:
            term_split = str.split(iterm, ' ')
            assert len(term_split)==2, (
                'phase and result symbols must be separated by a space.'
            )
            iphs, ichange = term_split

            assert ichange in cls.RESULT_OPTS, (
                'result is invalid. Every result must be selected from ' + str(cls.RESULT_OPTS)
            )

            phases.append(iphs)
            changes.append(ichange)

        return phases, changes
#===================================================
class ParamModel:
    pass
#===================================================
