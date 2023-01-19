import numpy as np
import re
import scipy as sp
from scipy import stats
from scipy import interpolate as interp
from scipy import random
# import emcee
import numdifftools as nd
import os
import pandas as pd
import ctypes
import copy
from scipy import optimize as optim
from ctypes import cdll
from ctypes import util
from rubicon.objc import ObjCClass, NSObject, objc_method
cdll.LoadLibrary(util.find_library('phaseobjc'))

def make_scale_matrix(array):
    scl_mat_a = np.dot(np.expand_dims(array,-1),
                       np.expand_dims(array,0))
    return scl_mat_a

def fill_array( var1, var2 ):
    """
    fix fill_array such that it returns two numpy arrays of equal size

    use numpy.full_like

    """
    var1_a = np.asarray( var1 )
    var2_a = np.asarray( var2 )

    if var1_a.shape==():
        var1_a = np.asarray( [var1] )
    if var2_a.shape==():
        var2_a = np.asarray( [var2] )

    # Begin try/except block to handle all cases for filling an array
    while True:
        try:
            assert var1_a.shape == var2_a.shape
            break
        except: pass
        try:
            var1_a = np.full_like( var2_a, var1_a )
            break
        except: pass
        try:
            var2_a = np.full_like( var1_a, var2_a )
            break
        except: pass

        # If none of the cases properly handle it, throw error
        assert False, 'var1 and var2 must both be equal shape or size=1'

    return var1_a, var2_a

def double_vector_to_array(vec):
    size = vec.size
    array = np.empty(size)
    for i in range(size):
        val = vec.valueAtIndex_(i)
        array[i] = val

    return array

def array_to_double_vector(array):
    doublevec_cls = ObjCClass('DoubleVector')
    # vec = (ctypes.c_double*array.size)()
    # ctypes.cast(vec, ctypes.POINTER(ctypes.c_double))
    vec = doublevec_cls.alloc().initWithSize_( array.size )
    vec_pointer = vec.pointerToDouble()
    for ind, val in enumerate(array):
        vec_pointer[ind] = val

    return vec

def double_matrix_to_array(mat):
    Nrow, Ncol = mat.rowSize, mat.colSize
    array = np.empty((Nrow,Ncol))
    for i in range(Nrow):
        for j in range(Ncol):
            val = mat.valueAtRowIndex_andColIndex(i,j)
            array[i,j] = val

    return array

def double_vector_to_array(vec):
    size = vec.size
    array = np.empty(size)
    for i in range(size):
        val = vec.valueAtIndex_(i)
        array[i] = val

    return array
#===================================================
class ThermoDBModel:
    def __init__(self, raw_df, thermoDB=None, mask_phs_l=['fSil'],
                 rxn_trans_typ='logistic', Tscl_energy=1.0):

        self.datadir = 'data'
        Rgas = 8.3144598
        self.Escl = 3.0/2*Rgas*Tscl_energy
        # self.Gthresh_scl = 3.0/2*Rgas*Tthresh_scl
        self.Gthresh_scl = self.Escl

        if thermoDB is None:
            thermoDB = ThermoDB()

        self.thermoDB = thermoDB
        self.rxn_trans_typ = rxn_trans_typ

        dat_df, rxn_d_l, phasesym_l = self.filter_phase_rev_data( raw_df, mask_phs_l=mask_phs_l )
        self.data_df = dat_df
        self._dat_trust_d = self.extract_trust_data()

        self.load_exp_prior_data()

        self.rxn_dir_opt = ['FWD','REV','NC','FWD?','REV?']


        self._param_d = {}
        self.init_model_phases( phasesym_l )
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


        parentpath = os.path.dirname(__file__)
        pathnm = os.path.join(parentpath,datadir,filenm)
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
        Rgas = 8.3144598
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
        num_dat_df = data_df[['P_kbar','Perr_kbar','T_C','Terr_C',
                              'time_hr']][trust_msk].astype(np.float)

        P_a = 1e3*num_dat_df['P_kbar'].values
        Perr_a = 1e3*num_dat_df['Perr_kbar'].values
        T_a = 273.15+num_dat_df['T_C'].values
        Terr_a = num_dat_df['Terr_C'].values
        # time_a = np.log10(num_dat_df['time_hr'].values) # NOTE: use log time
        time_a = num_dat_df['time_hr'].values # NOTE: use log time

        pubid_a = data_df['PubID'][trust_msk].values
        run_num_a = data_df['Run Num.'][trust_msk].values

        rxn_dir_a = data_df['RxnDir'][trust_msk].values
        rxn_a = data_df['Rxn'][trust_msk].values
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
        scl_mat_a = make_scale_matrix(curv_scl_a)

        Hnorm_fun = nd.Hessian(norm_costfun,step = 1e-2)
        Hnorm_a = Hnorm_fun(np.zeros(mu_a.shape))

        covnorm_a = np.linalg.pinv(Hnorm_a)

        cov_a = covnorm_a*scl_mat_a

        try:
            err_a = np.sqrt(np.diag(cov_a))
            # print(cov_a)
            err_scl_a = make_scale_matrix(err_a)
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
        err_scl_mat_a = make_scale_matrix(err_a)

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

    @classmethod
    def read_phase_rev_data(cls, filenm,sheetname=None):
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
    def detect_bound(cls, colnm, df):
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

    @classmethod
    def filter_phase_rev_data( cls, raw_df, mask_phs_l=None):
        P_df = cls.detect_bound('P_kbar',raw_df)
        T_df = cls.detect_bound('T_C',raw_df)
        time_df = cls.detect_bound('time_hr',raw_df)

        rxn_df = pd.DataFrame()

        # rxn_df['RxnPhases'] = raw_df[['Reaction_Studied']].applymap(get_reaction_str)
        rxn_df['RxnPhases'] = raw_df[['Reaction_Studied']].applymap(ThermoDB._get_reaction_phase_str)


        # Set Rxn equation as the most common string representaton in the raw database
        RxnPhases_uniq = rxn_df['RxnPhases'].unique()
        rxn_df['Rxn'] = raw_df['Reaction_Studied']
        for rxn_phs_str in RxnPhases_uniq:
            this_reaction = raw_df['Reaction_Studied'][rxn_df['RxnPhases']==rxn_phs_str]
            this_reaction = this_reaction.str.strip()
            # Store eqn only as most common variant
            #NOTE iloc crucial to obtain just value (not series object)
            rxn_df.loc[rxn_df['RxnPhases']==rxn_phs_str,'Rxn'] = this_reaction.mode().iloc[0]


        rxn_d_l = []
        phs_l = []
        rxn_eqn_uniq = rxn_df['Rxn'].unique()
        for rxn_eqn_str in rxn_eqn_uniq:
            rxn_d = ThermoDB.parse_rxn( rxn_eqn_str )
            #Rewrite Rxn equation using adopted rxn direction
            rxn_df.loc[rxn_df['Rxn']==rxn_eqn_str,'Rxn'] = rxn_d['rxn_eqn']

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
        trust_ser = raw_df['Trust'].fillna(value='Yes')
        if mask_phs_l is not None:
            # Remove from phase list
            phs_l = np.setdiff1d( phs_l, mask_phs_l )

            # set Trust variable to No for rxns involving masked phase
            for mask_phs in mask_phs_l:
                trust_ser[rxn_df['RxnPhases'].str.contains(mask_phs)] = 'No'



        phs_uniq_l = np.unique( phs_l )

        # Determine Reaction Direction: ['FWD','REV','NC','FWD?','REV?','INV']
        rxn_dir_l = []
        for result, rxn_phs_str in zip(raw_df['Results'],rxn_df['RxnPhases']):
            rxn_dir_l.append(ThermoDB._get_rxn_dir(rxn_phs_str, result))



        # rxn_uniq = rxn_df['Rxn'].unique()
        # print(rxn_uniq)
        rxn_df['RxnDir'] = pd.Series(rxn_dir_l)
        # rxn_df['Results'] = raw_df['Results']

        # result = rxn_df[['Rxn']].applymap(get_reaction_phase_str)
        # print(get_reaction_phases(result.loc[0,'Rxn']))

        # print(result)
        #  print(result['Rxn'].unique())
        #  rxn_df = pd.concat((rxn_df,pd.DataFrame(result['Rxn'].tolist(),columns=['Reac_l','Prod_l'])),axis=1)

        # print(result.tolist())
        # print(pd.DataFrame(result,columns=['Reac','Prod']))

        dat_df = pd.concat((raw_df['PubID'],raw_df['Aparatus'],raw_df['Run Num.'],
                            time_df,raw_df['Water'],
                            P_df,raw_df['Perr_kbar'],T_df,raw_df['Terr_C'],
                            rxn_df,trust_ser),axis=1)

        return dat_df, rxn_d_l, phs_uniq_l

    @classmethod
    def rxn_trans_fun(self, x_a,  rxn_trans_typ='logistic'):
        if rxn_trans_typ == 'logistic':
            const = 1.8138 # pi/sqrt(3)
            # F_a = 1.0/(1+np.exp(-const*x_a))
            # Special optimized version of logistic function
            F_a = sp.special.expit(const*x_a)
        elif (rxn_trans_typ == 'normal') | (rxn_trans_typ == 'erf'):
            const = 0.70711 # 1/sqrt(2)
            F_a = 0.5*(1+sp.special.erf(const*x_a))

        return F_a

    @classmethod
    def rxn_logtrans_fun(cls, x_a, rxn_trans_typ='logistic'):
        if rxn_trans_typ == 'logistic':
            const = 1.8138 # pi/sqrt(3)
            # Special optimized version of logistic function
            logF_a = -np.logaddexp(0,-const*x_a)
        elif (rxn_trans_typ == 'normal') | (rxn_trans_typ == 'erf'):
            const = 0.70711 # 1/sqrt(2)
            # Special optimized version of log-cdf for normal distribution
            logF_a = sp.special.log_ndtr(x_a)

        return logF_a

    @classmethod
    def logprior_fun(cls, x_a, typ='studentt', dof=5):
        if typ == 'studentt':
            # Variance of student's t distribution is slightly larger than a
            # normal (depending on dof). Thus the relative residual x must be
            # scaled down to match the desired standard deviation.
            const = np.sqrt(1.0*dof/(dof-2))
            logprob_a = stats.t.logpdf(x_a/const,dof)
        elif (typ == 'normal') | (rxn_trans_typ == 'erf'):
            logprob_a = stats.norm.log_pdf(x_a)

        return logprob_a

    @classmethod
    def logprob_fwd(cls, dGrxn_a,Gth_a,sigG_a, rxn_trans_typ='logistic'):
        x_a = -(dGrxn_a+Gth_a)/sigG_a
        logprob_a = cls.rxn_logtrans_fun(x_a,rxn_trans_typ=rxn_trans_typ)
        return logprob_a

    @classmethod
    def logprob_rev(cls, dGrxn_a,Gth_a,sigG_a, rxn_trans_typ='logistic'):
        x_a = +(dGrxn_a-Gth_a)/sigG_a
        logprob_a = cls.rxn_logtrans_fun(x_a,rxn_trans_typ=rxn_trans_typ)
        return logprob_a

    @classmethod
    def logprob_fwdq(cls, dGrxn_a,Gth_a,sigG_a, rxn_trans_typ='logistic'):
        N=dGrxn_a.size
        logPrev_a = cls.logprob_rev(dGrxn_a,Gth_a,sigG_a,
                                    rxn_trans_typ=rxn_trans_typ)
        logterms_a = np.vstack((np.zeros(N),logPrev_a))
        scl_a      = np.vstack((np.ones(N),-np.ones(N)))
        logprob_a = sp.misc.logsumexp(logterms_a,axis=0,b=scl_a)
        return logprob_a

    @classmethod
    def logprob_revq(cls, dGrxn_a,Gth_a,sigG_a, rxn_trans_typ='logistic'):
        N=dGrxn_a.size
        logPfwd_a = cls.logprob_fwd(dGrxn_a,Gth_a,sigG_a,
                                    rxn_trans_typ=rxn_trans_typ)
        logterms_a = np.vstack((np.zeros(N),logPfwd_a))
        scl_a      = np.vstack((np.ones(N),-np.ones(N)))
        logprob_a = sp.misc.logsumexp(logterms_a,axis=0,b=scl_a)
        return logprob_a

    @classmethod
    def logprob_nc(cls, dGrxn_a,Gth_a,sigG_a, rxn_trans_typ='logistic'):
        N=dGrxn_a.size
        logPfwd_a = cls.logprob_fwd(dGrxn_a,Gth_a,sigG_a,
                                    rxn_trans_typ=rxn_trans_typ)
        logPrev_a = cls.logprob_rev(dGrxn_a,Gth_a,sigG_a,
                                    rxn_trans_typ=rxn_trans_typ)
        logterms_a = np.vstack((np.zeros(N),logPfwd_a,logPrev_a))
        scl_a      = np.vstack((np.ones(N),-np.ones(N),-np.ones(N)))
        logprob_a = sp.misc.logsumexp(logterms_a,axis=0,b=scl_a)
        logprob_a[np.isnan(logprob_a)] = -np.inf
        return logprob_a

class ThermoDB:
    def __init__( self, database='Berman', exp_prior=''):

        self.database = database
        self.datadir = 'data'
        self.filebasenm_pure = 'PurePhases.csv'

        self.load_full_pure_phase_list()
        self.load_pure_phase_data()

        self.init_pure_database()

    def load_full_pure_phase_list( self, filenm='PurePhaseList.csv' ):
        parentpath = os.path.dirname(__file__)
        pathname = os.path.join(parentpath,self.datadir,filenm)
        try:
            all_purephases_df = pd.read_csv(pathname)
        except:
            assert False,'The '+filenm+' file cannot be found. '\
                'It is needed to define the standard phase abbreviations.'

        self.all_purephases_df = all_purephases_df
        self.all_purephases_pathnm = pathname
        pass

    def load_pure_phase_data( self ):
        filebasenm = self.filebasenm_pure
        datadir = self.datadir
        database = self.database

        parentpath = os.path.dirname(__file__)
        datafile = database+filebasenm
        pathnm = os.path.join(parentpath,datadir,datafile)
        purephaselib_df = pd.read_csv(pathnm)

        # Verify that all phases in library are part of full pure phase list
        # [phasesym for phasesym in purephaselib_df['Abbrev']]

        abbrev_valid = purephaselib_df['Abbrev'].isin(self.all_purephases_df['Abbrev'])
        assert abbrev_valid.all(), 'The pure phase library defined in "'+pathnm\
            +'" contains some invalid phase abbreviations, shown below: \n\n'\
            +str(purephaselib_df[~abbrev_valid])\
            +'\n\nCheck that the abbreviations conform to the list given in: "'\
            +self.all_purephases_pathnm+'"'


        self.purephaselib = pathnm
        self.purephaselib_df = purephaselib_df

        pass


    def init_pure_database( self, fixH2O=True ):
        purephaselib_df = self.purephaselib_df

        # Create and store class for each phase
        purephase_props_d = {}
        purephase_cls_d = {}
        purephase_obj_d = {}
        for indx, purephase in purephaselib_df.iterrows():
            abbrev = purephase['Abbrev']
            classnm = purephase['ClassName']
            try:

                phase_cls = ObjCClass(classnm+self.database)
                # phs = PurePhase(classnm, database=self.database)
            except:
                assert False, classnm+' is not a valid ClassName for the '\
                    +self.database+' database.'

            phase_obj = PurePhase(phase_cls, abbrev)

            # Create temporary phase object to determine fixed phase properties
            # props_d = phase_obj.props_d
            # props_d = self.get_phase_props(phase_cls)
            purephase_props_d[abbrev] = phase_obj.props_d
            purephase_cls_d[abbrev] = phase_cls
            purephase_obj_d[abbrev] = phase_obj


        if fixH2O:
            H2O_phase_cls = ObjCClass('GenericH2O')
            H2O_phase_obj = PurePhase(H2O_phase_cls, 'H2O', calib=False)
            purephase_props_d['H2O'] = H2O_phase_obj.props_d
            purephase_cls_d['H2O'] = H2O_phase_cls
            purephase_obj_d['H2O'] = H2O_phase_obj

        self.purephase_props_d = purephase_props_d
        self.purephase_cls_d = purephase_cls_d
        self.purephase_obj_d = purephase_obj_d
        pass

    def enable_gibbs_energy_reference_state(self):
        # call method on any phase class (automatically applied to all)
        key0 = next(iter(self.purephase_cls_d))
        phs0_cls = self.purephase_cls_d[key0]
        phs0_cls.enableGibbsFreeEnergyReferenceStateUsed()
        pass

    def disable_gibbs_energy_reference_state(self):
        # call method on any phase class (automatically applied to all)
        key0 = next(iter(self.purephase_cls_d))
        phs0_cls = self.purephase_cls_d[key0]
        phs0_cls.disableGibbsFreeEnergyReferenceStateUsed()
        pass

    def get_phase_cls(self, phasesym_l):
        phase_cls_l = []
        for phasesym in phasesym_l:
            try:
                phase_cls_l.append(self.purephase_cls_d[phasesym])
            except:
                assert False, '"'+phasesym+'" is not a valid phase abbreviation. Try again.'

        return phase_cls_l

    def get_phase_obj(self, phasesym_l):
        phase_obj_l = []
        for phasesym in phasesym_l:
            try:
                phase_obj_l.append(self.purephase_obj_d[phasesym])
            except:
                assert False, '"'+phasesym+'" is not a valid phase abbreviation. Try again.'

        return phase_obj_l

    def new_phase(self, phasesym ):
        # phase_cls = self.get_phase_cls([phasesym])[0]
        # phase = PurePhase(phase_cls, phasesym)
        phase = self.purephase_obj_d[phasesym]
        return phase

    def new_assemblage(self, phasesym_l ):
        phase_obj_l = []
        for phs_sym in self.purephase_obj_d:
            if phs_sym in phasesym_l:
                phase_obj_l.append(self.purephase_obj_d[phs_sym])

        # phase_cls_l = self.get_phase_cls(phasesym_l)
        # assemblage = PhaseAssemblage( phase_cls_l, phasesym_l )
        assemblage = PhaseAssemblage( phase_obj_l, phasesym_l, obj_is_class=False )
        return assemblage

    def new_rxn(self, reac_phasesym_l, prod_phasesym_l ):

        # reac_phase_cls_l = self.get_phase_cls(reac_phasesym_l)
        # prod_phase_cls_l = self.get_phase_cls(prod_phasesym_l)
        # rxn = PhaseRxn( reac_phase_cls_l, reac_phasesym_l,
        #                prod_phase_cls_l, prod_phasesym_l  )
        reac_phase_obj_l = self.get_phase_obj(reac_phasesym_l)
        prod_phase_obj_l = self.get_phase_obj(prod_phasesym_l)
        rxn = PhaseRxn( reac_phase_obj_l, reac_phasesym_l,
                       prod_phase_obj_l, prod_phasesym_l, obj_is_class=False )
        return rxn

    @classmethod
    def parse_rxn( cls, rxn_eqn_str, rxn_result_str=None, sort=True ):
        rxn_phs_str, rxn_eqn_str = cls._get_reaction_phase_str( rxn_eqn_str, sort=sort, full_output=True)
        reac_l, prod_l = cls._get_reaction_phases( rxn_phs_str )

        if rxn_result_str is not None:
            rxn_dir = cls._get_rxn_dir(rxn_phs_str, rxn_result_str)
        else:
            rxn_dir = None

        rxn_d = {}
        rxn_d['rxn_eqn'] = rxn_eqn_str
        rxn_d['rxn_phs_str'] = rxn_phs_str
        rxn_d['reac_l'] = reac_l
        rxn_d['prod_l'] = prod_l
        rxn_d['rxn_dir'] = rxn_dir

        return rxn_d

    @classmethod
    def _parse_rxn_result(cls, result_str):
        # Remove surrounding whitespace
        result_str = result_str.strip()
        if result_str == 'NC':
            phs_l = None
            obs_l = result_str
            return phs_l, obs_l

        parse_result = re.compile(r'\w+\s*[+-?]+\s*')
        phs_result_l = parse_result.findall(result_str)
        phs_a = []
        obs_a = []
        for iphs_result in phs_result_l:
            ires = iphs_result.strip().split()
            phs_a.append(ires[0])
            obs_a.append(ires[1])

        phs_a = np.array(phs_a)
        obs_a = np.array(obs_a)

        return phs_a, obs_a

    @classmethod
    def _get_rxn_dir(cls, rxn_phs_str, result):
        # Determine Reaction Direction: ['FWD','REV','NC','FWD?','REV?','INV']
        # print('rxn_phs_str = ',rxn_phs_str)
        # print('results = ',result)
        # print('===========')

        reac_l, prod_l = cls._get_reaction_phases(rxn_phs_str)
        phs_a, obs_a = cls._parse_rxn_result(result)

        if phs_a is None:
            rxn_dir = 'NC'
            return rxn_dir

        fwd_cnt = 0
        fwd_maybe_cnt = 0
        rev_cnt = 0
        rev_maybe_cnt = 0

        for reac in reac_l:
            if reac not in phs_a:
                continue
            obs = obs_a[phs_a==reac]
            if obs == '-':
                fwd_cnt += 1
            elif obs == '-?':
                fwd_maybe_cnt +=1
            elif obs == '+':
                rev_cnt += 1
            elif obs == '+?':
                rev_maybe_cnt += 1
            elif obs == '?':
                # Do nothing
                pass
            else:
                # print('*****************')
                # print('reac = ',reac)
                # print('prod = ',prod)
                # print('obs = ', obs)
                # print('prod_l = ',prod_l)
                # print('reac_l = ',reac_l)
                # print('phs_a = ',phs_a)
                # print('obs_a = ',obs_a)
                # print('*****************')
                assert False, obs + ' is not a supported Results code.'

        for prod in prod_l:
            if prod not in phs_a:
                continue
            obs = obs_a[phs_a==prod]
            if obs == '+':
                fwd_cnt += 1
            elif obs == '+?':
                fwd_maybe_cnt +=1
            elif obs == '-':
                rev_cnt += 1
            elif obs == '-?':
                rev_maybe_cnt += 1
            elif obs == '?':
                # Do nothing
                pass
            else:
                # print('*****************')
                # print('reac = ',reac)
                # print('prod = ',prod)
                # print('obs = ', obs)
                # print('prod_l = ',prod_l)
                # print('reac_l = ',reac_l)
                # print('phs_a = ',phs_a)
                # print('obs_a = ',obs_a)
                # print('*****************')
                assert False, obs + ' is not a supported Results code.'

        # Reaction Direction Options: ['FWD','REV','NC','FWD?','REV?','INV']
        if (fwd_cnt==0)&(fwd_maybe_cnt==0)&(rev_cnt==0)&(rev_maybe_cnt==0):
            rxn_dir = 'INV'
        if fwd_maybe_cnt > 0:
            rxn_dir = 'FWD?'
        if rev_maybe_cnt > 0:
            rxn_dir = 'REV?'
        if (fwd_maybe_cnt > 0) & (rev_maybe_cnt > 0):
            rxn_dir = 'NC'
        if fwd_cnt > 0:
            rxn_dir = 'FWD'
        if rev_cnt > 0:
            rxn_dir = 'REV'
        if (fwd_cnt > 0) & (rev_cnt > 0):
            rxn_dir = 'INV'

        return rxn_dir

    @classmethod
    def _get_reaction_phases(cls, rxn_phs_str):
        reac_str, prod_str = str.split(rxn_phs_str,':')
        reac_l = reac_str.split(',')
        prod_l = prod_str.split(',')
        return reac_l, prod_l

    @classmethod
    def _split_phases(cls, phs_combo_str):
        return [re.sub('^[0-9]','',phs.strip()) \
                for phs in phs_combo_str.split('+')]

    @classmethod
    def _get_reaction_phase_str(cls, rxn_eqn_str, sort=True, full_output=False ):
        reac_str, prod_str = str.split(rxn_eqn_str,'=')
        reac_str = str.strip(reac_str)
        prod_str = str.strip(prod_str)

        reac_l = cls._split_phases(reac_str)
        prod_l = cls._split_phases(prod_str)

        if sort:
            reac_l.sort()
            prod_l.sort()

            first_phs = [reac_l[0],prod_l[0]]
            first_phs.sort()
            if first_phs[0] == prod_l[0]:
                reac_l, prod_l = prod_l, reac_l
                rxn_eqn_str = prod_str + ' = ' + reac_str

        rxn_phs_str = ''
        for phs in reac_l[:-1]:
            rxn_phs_str += phs+','

        rxn_phs_str += reac_l[-1]
        rxn_phs_str += ':'
        for phs in prod_l[:-1]:
            rxn_phs_str += phs+','

        rxn_phs_str += prod_l[-1]

        if full_output:
            return rxn_phs_str, rxn_eqn_str
        else:
            return rxn_phs_str

class PhaseRxn:
    def __init__( self, reac_obj_l, reac_sym_l,
                 prod_obj_l, prod_sym_l, obj_is_class=False, TOL=1e-8 ):
        # has a Phase assemblage internally
        self._reac_assemblage = PhaseAssemblage( reac_obj_l, reac_sym_l,
                                                obj_is_class=obj_is_class )
        self._prod_assemblage = PhaseAssemblage( prod_obj_l, prod_sym_l,
                                                obj_is_class=obj_is_class )
        self.init_rxn(TOL)
        pass

    def init_rxn( self, TOL, TOL_COEF=1e-2 ):
        props_reac_d = self.reac_assemblage.props_d
        props_prod_d = self.prod_assemblage.props_d
        reac_abbrev_l = props_reac_d['abbrev_all_l']
        prod_abbrev_l = props_prod_d['abbrev_all_l']


        elem_formula_all_reac_a = props_reac_d['elem_formula_all_a']
        elem_formula_all_prod_a = props_prod_d['elem_formula_all_a']


        elem_formula_all_a = np.vstack((elem_formula_all_reac_a,
                                        elem_formula_all_prod_a))

        atomnum = 0.5*np.sum(elem_formula_all_a)

        # ind_nonzero_a = ~np.all(elem_formula_all_a==0,axis=1)
        ind_nonzero_a = ~np.all(elem_formula_all_a==0,axis=0)

        stoich_a = elem_formula_all_a[:,ind_nonzero_a]
        Nphase = stoich_a.shape[0]

        #######################
        #  Explore multiple possible reactions
        #######################

        #    U,s,Vh = np.linalg.svd(stoich_a)
        #    print('U')
        #    print(U)
        #    print('s')
        #    print(s)
        #    print('Vh')
        #    print(Vh)
        #    print('stoich')
        #    print(stoich_a)
        #    print('--------------')


        #    coefs_a = np.zeros((Nphase,Nphase))
        #    resid_a = np.zeros(Nphase)

        #    for i in range(Nphase):
        #        coefs_a = np.roll(coefs_a,1,axis=1)
        #        stoich_a = np.roll(stoich_a,1,axis=0)
        #        # print(stoich_a)
        #        xbar = np.dot(np.linalg.pinv(stoich_a[0:-1].T),stoich_a[-1])
        #        icoefs_a = np.hstack((xbar,-1))/np.min(np.abs(xbar))
        #        iresid = np.sqrt(np.sum((np.dot(stoich_a[0:-1].T,xbar)-stoich_a[-1])**2))
        #        resid_a[i] = iresid
        #        # print(icoefs_a)
        #        coefs_a[i] = icoefs_a
        #        # print('====')

        #    print('resid_a')
        #    print(resid_a)
        #    print('residmax')
        #    print(np.max(resid_a))

        #    decimals=8

        #    for ind in range(Nphase):
        #        if coefs_a[ind,0] <0:
        #            coefs_a[ind] *= -1

        #    print('max scale')
        #    coef_scl_a = np.round(np.vstack([coef/np.max(np.abs(coef)) for coef in coefs_a]),decimals=decimals+1)
        #    print(coef_scl_a)
        #    coef_scl2_a = np.round(np.vstack([coef/np.min(np.abs(coef[coef!=0])) for coef in coef_scl_a]),decimals=decimals)
        #    print(coef_scl2_a)

        #    ceofs_a = coef_scl2_a

        #    print(self.assemblage.phasesym_l)
        #    # print([key for key in props_d])
        #    # print(props_d['name_all_l'])
        #    coefs_a = np.round(coefs_a,decimals=6)
        #    print(coefs_a)
        #    # print(np.unique([tuple(row) for row in coefs_a]))
        #    coefs_uniq_a = np.vstack({tuple(row) for row in coefs_a})
        #    print('coefs uniq')
        #    print(coefs_uniq_a)
        # print(coefs_uniq_a[0]-coefs_uniq_a[1])
        # print(coefs_uniq_a[0]+coefs_uniq_a[1])
        # print(1./21*(coefs_uniq_a[0]-coefs_uniq_a[1])+coefs_uniq_a[2])

        # print(np.dot(stoich_a[0:-1].T,xbar))
        # print(stoich_a[-1])
        # resid=np.dot(stoich_a[0:-1].T,xbar)-stoich_a[-1]
        # print(resid)

        # U,s,Vh = np.linalg.svd(stoich_a)
        # print('U')
        # print(U)
        # print('s')
        # print(s)
        # print('Vh')
        # print(Vh)
        # print('stoich')
        # print(stoich_a)

        # print('Vh*stoich_a')
        # # print(np.dot(Vh.T,stoich_a))
        # print(np.dot(stoich_a,Vh.T))

        # print('=============')

        # Least Squares solution to a*x = b
        a = stoich_a[0:-1].T
        b = stoich_a[-1]

        results = np.linalg.lstsq(a,b)

        coef_a = results[0]
        resid = results[1]
        # print(results)

        assert len(resid)==1, 'This set of phases does not represent a valid reaction.'
        assert resid<TOL, 'This set of phases does not represent a complete reaction (large residuals).'

        coef_a = np.hstack((coef_a,-1))

        mincoef = np.min(np.abs(coef_a))
        maxcoef = np.max(np.abs(coef_a))
        assert mincoef/maxcoef > TOL, 'This set of phases includes additional phases that do not participate in the reaction.'

        coef_a /= mincoef

        # self._rxn_coef_a = coef_a
        reac_coef_a = +coef_a[coef_a>0]
        prod_coef_a = -coef_a[coef_a<0]

        rxn_eqn_str = ''
        for coef,reac in zip(reac_coef_a,reac_abbrev_l):
            coef_rnd = np.round(coef,0)
            if np.abs(coef_rnd-coef)<TOL_COEF:
                coef = int(coef_rnd)
            rxn_eqn_str += str(coef)+' '+reac+' '

        rxn_eqn_str += '='
        for coef,prod in zip(prod_coef_a,prod_abbrev_l):
            coef_rnd = np.round(coef,0)
            if np.abs(coef_rnd-coef)<TOL_COEF:
                coef = int(coef_rnd)
            rxn_eqn_str += ' '+str(coef)+' '+prod

        self.rxn_eqn_str = rxn_eqn_str
        self._reac_rxn_coef_a = reac_coef_a
        self._prod_rxn_coef_a = prod_coef_a

        self.reac_abbrev_l = reac_abbrev_l
        self.prod_abbrev_l = prod_abbrev_l
        self._rxn_atomnum = atomnum

        pass

    def __eq__(self, other):
        return (self.reac_assemblage == other.reac_assemblage) & \
            (self.prod_assemblage == other.prod_assemblage)

    @property
    def reac_assemblage(self):
        return self._reac_assemblage

    @property
    def prod_assemblage(self):
        return self._prod_assemblage

    @property
    def reac_rxn_coef_a(self):
        return self._reac_rxn_coef_a

    @property
    def prod_rxn_coef_a(self):
        return self._prod_rxn_coef_a

    @property
    def rxn_atomnum(self):
        return self._rxn_atomnum

    def get_rxn_gibbs_energy( self, T_a, P_a, peratom=False ):
        dG_rxn_a = self._calc_rxn_change( 'get_gibbs_energy_all', T_a, P_a,
                                         peratom=peratom )
        return dG_rxn_a

    def get_rxn_entropy( self, T_a, P_a, peratom=False ):
        dS_rxn_a = self._calc_rxn_change('get_entropy_all', T_a, P_a,
                                        peratom=peratom )
        return dS_rxn_a

    def get_rxn_enthalpy( self, T_a, P_a, peratom=False ):
        dH_rxn_a = self._calc_rxn_change('get_enthalpy_all', T_a, P_a,
                                         peratom=peratom )
        return dH_rxn_a

    def get_rxn_volume( self, T_a, P_a, peratom=False ):
        dV_rxn_a = self._calc_rxn_change('get_volume_all', T_a, P_a,
                                         peratom=peratom )
        return dV_rxn_a

    def get_rxn_bound( self, T=None, P=None, init_guess=None ):
        Tstd = 300.0
        Pstd = 1.0

        assert (T is not None) or (P is not None), \
            'Both T and P are None; Must define either temperature or pressure.'

        if T is not None:
            if init_guess is None:
                init_guess = Pstd

            Gfun = lambda P, T=T: self.get_rxn_gibbs_energy(T,P)
            dGfun = lambda P, T=T: self.get_rxn_volume(T,P)

        else:
            if init_guess is None:
                init_guess = Tstd

            Gfun = lambda T, P=P: self.get_rxn_gibbs_energy(T,P)
            dGfun = lambda T, P=P: -self.get_rxn_entropy(T,P)

        value_bound = optim.newton(Gfun,init_guess,fprime=dGfun)

        return value_bound

    def get_clapeyron_slope( self, T, P ):
        dTdP = self.get_rxn_volume(T,P)/self.get_rxn_entropy(T,P)
        return dTdP

    def trace_rxn_bound( self, Tlims=None, Plims=None, init_guess=None,
                        Nsamp=30 ):

        assert (Tlims is not None) or (Plims is not None), \
            'Both T and P are None; Must define either temperature or pressure.'

        if Tlims is not None:
            xval_a = np.linspace( Tlims[0], Tlims[1], Nsamp )

            bound_fun = lambda T, init_guess: self.get_rxn_bound\
                (T=T, init_guess=init_guess)

            dGdx_fun = lambda x, y: -self.get_rxn_entropy(x,y)
            dGdy_fun = lambda x, y: +self.get_rxn_volume (x,y)

        else:
            xval_a = np.linspace( Plims[0], Plims[1], Nsamp )

            bound_fun = lambda P, init_guess: self.get_rxn_bound\
                (P=P, init_guess=init_guess)

            dGdx_fun = lambda x, y: +self.get_rxn_volume (y,x)
            dGdy_fun = lambda x, y: -self.get_rxn_entropy(y,x)

        dx = xval_a[1]-xval_a[0]
        yval_a = np.zeros( xval_a.shape )

        for ind,x in enumerate(xval_a):
            y = bound_fun(x,init_guess)
            dy_guess = dx*dGdx_fun(x=x,y=y)/dGdy_fun(x=x,y=y)
            yval_a[ind] = y
            init_guess = y + dy_guess


        if Tlims is not None:
            TP_a = np.vstack((xval_a,yval_a))
        else:
            TP_a = np.vstack((yval_a,xval_a))

        return TP_a

    def get_rxn_stability( self, TP_bound_a, other_rxn_l, TOL=1e-6 ):
        # NOTE: other_rxn_l must be a list of reactions of the same composition
        peratom=True

        self._reac_assemblage
        G_rxn_a = self.get_rxn_gibbs_energy( TP_bound_a[0], TP_bound_a[1],
                                            peratom=peratom )

        assert np.all(np.abs(G_rxn_a)<TOL), \
            'Energy diff. must be less than TOL.'

        N_TP = G_rxn_a.size
        N_rxn = len(other_rxn_l)

        G_rxn_other_a = np.zeros((N_rxn,N_TP))

        # Calculate gibbs energy of rxn for all possible reactions

        for ind,rxn in enumerate(other_rxn_l):
            iG_rxn_other_a = rxn.get_rxn_gibbs_energy( TP_bound_a[0],
                                                      TP_bound_a[1],
                                                      peratom=peratom )
            # if other rxn phase assemblage is a subset of the current
            # assemblage, then it is a valid competing reaction

            if (rxn._reac_assemblage.issubset(self._reac_assemblage) | \
                rxn._reac_assemblage.issubset(self._prod_assemblage) ):
                G_rxn_other_a[ind] = iG_rxn_other_a

            elif (rxn._prod_assemblage.issubset(self._reac_assemblage) | \
                  rxn._prod_assemblage.issubset(self._prod_assemblage) ):
                # Store negative energy as reverse reaction energy change
                G_rxn_other_a[ind] = -iG_rxn_other_a

            else:
                G_rxn_other_a[ind] = 0.0

        # print(G_rxn_other_a)


        # If any rxn is energetically favored, then current reaction is not
        # stable
        stable_a = ~np.any(G_rxn_other_a < -TOL, axis=0)
        return stable_a

    def get_simultaneous_rxn_cond( self, other_rxn, Pinit=1.0, TOL=1e-5 ):
        P = Pinit

        while True:
            Tbnd1 = self.get_rxn_bound(P=P)
            Tbnd2 = other_rxn.get_rxn_bound(P=P)

            dTdP1 = self.get_clapeyron_slope(Tbnd1,P)
            dTdP2 = other_rxn.get_clapeyron_slope(Tbnd2,Pinit)

            if np.abs(np.log(Tbnd1/Tbnd2)) < TOL:
                T =  0.5*(Tbnd1+Tbnd2)
                break

            dP = - (Tbnd1 - Tbnd2)/(dTdP1-dTdP2)
            P += dP

        T = float(T)
        P = float(P)

        return T, P

    def _calc_reac_value( self, method_name, T_a, P_a, peratom=False ):
        reac_method = getattr(self.reac_assemblage, method_name)
        val_phs_a = reac_method( T_a, P_a )
        val_a = np.dot(self.reac_rxn_coef_a, val_phs_a)

        if peratom:
            val_a /= self.rxn_atomnum

        return val_a

    def _calc_prod_value( self, method_name, T_a, P_a, peratom=False ):
        prod_method = getattr(self.prod_assemblage, method_name)
        val_phs_a = prod_method( T_a, P_a )
        val_a = np.dot(self.prod_rxn_coef_a, val_phs_a)

        if peratom:
            val_a /= self.rxn_atomnum

        return val_a

    def _calc_rxn_change( self, method_name, T_a, P_a, peratom=False ):
        val_prod_a = self._calc_prod_value( method_name, T_a, P_a, peratom=peratom )
        val_reac_a = self._calc_reac_value( method_name, T_a, P_a, peratom=peratom )
        val_rxn_a = val_prod_a-val_reac_a
        return val_rxn_a

class PhaseAssemblage:
    def __init__( self, phase_obj_l, phasesym_l, obj_is_class=False ):
        indsort = np.argsort(phasesym_l)
        phase_obj_l = np.array(phase_obj_l)
        phasesym_l = np.array(phasesym_l)
        # N = len(indsort)
        # for ind in
        # print(indsort)
        # print(phase_cls_l)
        # print(phasesym_l[indsort])
        self.set_phase_assemblage(phase_obj_l[indsort],phasesym_l[indsort],
                                  obj_is_class=obj_is_class)
        pass

    def __eq__(self, other):
        return self._phase_l == other._phase_l

    def __lt__(self, other):
        if len(self._phase_l) < len(other._phase_l):
            return True

        return self._phase_l[0] < other._phase_l[0]

    def __gt__(self, other):
        if len(self._phase_l) > len(other._phase_l):
            return True

        return self._phase_l[0] > other._phase_l[0]

    def issubset(self, other):
        is_member_l = [phase in other._phase_l for phase in self._phase_l]
        return np.all(is_member_l)

    def set_phase_assemblage( self, phase_obj_l, phasesym_l,
                             obj_is_class=False ):
        if obj_is_class:
            phase_cls_l = phase_obj_l
            phase_obj_l = [PurePhase(phase_cls,phasesym) for \
                           (phase_cls,phasesym) in zip(phase_cls_l,phasesym_l)]

        self._phasesym_l = phasesym_l
        self._phase_l = phase_obj_l

        props_d = {}
        props_d['formula_all_l'] = [ x.props_d['formula'] for x in phase_obj_l]
        props_d['name_all_l'] = [ x.props_d['name'] for x in phase_obj_l ]
        props_d['abbrev_all_l'] = [ x.props_d['abbrev'] for x in phase_obj_l ]
        props_d['mw_all_a'] = np.array([ x.props_d['mw'] for x in phase_obj_l])
        props_d['elem_formula_all_a'] = \
            np.vstack([ x.props_d['elem_formula_a'] for x in phase_obj_l ] )
        props_d['elem_sym_a'] = phase_obj_l[0].props_d['elem_sym_a']

        self._props_d = props_d

        pass

    @property
    def phasesym_l(self):
        return self._phasesym_l

    @property
    def phase_l(self):
        return self._phase_l

    @property
    def props_d(self):
        return self._props_d

    def get_gibbs_energy_all( self, T_a, P_a ):
        return np.vstack([ x.get_gibbs_energy(T_a,P_a) for x in self.phase_l ])

    def get_enthalpy_all( self, T_a, P_a ):
        return np.vstack([ x.get_enthalpy(T_a,P_a) for x in self.phase_l ])

    def get_entropy_all( self, T_a, P_a ):
        return np.vstack([ x.get_entropy(T_a,P_a) for x in self.phase_l ])

    def get_heat_capacity_all( self, T_a, P_a ):
        return np.vstack([ x.get_heat_capacity(T_a,P_a) for x in self.phase_l ])

    def get_dCp_dT_all( self, T_a, P_a ):
        return np.vstack([ x.get_dCp_dT(T_a,P_a) for x in self.phase_l ])

    def get_volume_all( self, T_a, P_a ):
        return np.vstack([ x.get_volume(T_a,P_a) for x in self.phase_l ])

    def get_dV_dT_all( self, T_a, P_a ):
        return np.vstack([ x.get_dV_dT(T_a,P_a) for x in self.phase_l ])

    def get_dV_dP_all( self, T_a, P_a ):
        return np.vstack([ x.get_dV_dP(T_a,P_a) for x in self.phase_l ])

    def get_d2V_dT2_all( self, T_a, P_a ):
        return np.vstack([ x.get_d2V_dT2(T_a,P_a) for x in self.phase_l ])

    def get_d2V_dTdP_all( self, T_a, P_a ):
        return np.vstack([ x.get_d2V_dTdP(T_a,P_a) for x in self.phase_l ])

    def get_d2V_dP2_all( self, T_a, P_a ):
        return np.vstack([ x.get_d2V_dP2(T_a,P_a) for x in self.phase_l ])

class PurePhase:

    def __init__(self, phase_cls, abbrev, calib=True):
        if abbrev=='H2O':
            calib=False

        self._init_phase_props(phase_cls,abbrev, calib=calib)
        pass

    def _init_phase_props(self, phase_cls, abbrev, calib=True):
        phase_obj = phase_cls.alloc().init()

        self._phase_obj = phase_obj

        props_d = {}
        props_d['class_name'] = phase_cls.name
        props_d['abbrev'] = abbrev
        props_d['formula'] = phase_obj.phaseFormula
        props_d['name'] = phase_obj.phaseName
        props_d['mw'] = phase_obj.mw

        elem_formula_a = double_vector_to_array(phase_obj.formulaAsElementArray)
        props_d['elem_formula_a'] = elem_formula_a
        props_d['entropy_Robie1979_a'] = phase_obj.entropyFromRobieEtAl1979

        Natom = np.sum(elem_formula_a)
        Nelem = elem_formula_a.size
        props_d['Natom'] = Natom
        # NOTE: Why is Hg missing (represented as None)?
        elem_sym_a = np.array([phase_cls.elementNameFromAtomicNumber_(i) \
                               for i in range(Nelem)])

        props_d['elem_sym_a'] = elem_sym_a
        # props_d['supports_param_calib'] =  phase_obj.supportsParameterCalibration()

        if calib:
            supports_calib = phase_obj.supportsParameterCalibration()
            props_d['supports_calib'] = supports_calib

            Nparam = phase_obj.getNumberOfFreeParameters()
            props_d['param_num'] = Nparam

            param_names_NSArray = phase_obj.getArrayOfNamesOfFreeParameters()
            param_names_a = [param_names_NSArray.objectAtIndex_(i) for i in range(Nparam)]
            props_d['param_names_a'] =  param_names_a

            param_units_a = np.array([phase_obj.getUnitsForParameterName_(key) \
                                      for key in param_names_a])
            props_d['param_units_a'] =  param_units_a


            # Store initial param values (read from ObjC class)
            props_d['param0_a'] = np.array([phase_obj.getValueForParameterName_(key) \
                                            for key in param_names_a])


        self._props_d = props_d
        pass

    def __eq__(self, other):
        return self._props_d['abbrev'] == other._props_d['abbrev']

    def __lt__(self, other):
        return self._props_d['abbrev'] <= other._props_d['abbrev']

    def __gt__(self, other):
        return self._props_d['abbrev'] >= other._props_d['abbrev']

    @property
    def props_d(self):
        return self._props_d

    @property
    def phase_obj(self):
        return self._phase_obj

    def convert_elements_to_moles_of_phase(self, elements):
        elements_vec = array_to_double_vector(elements)
        return self._phase_obj.convertElementsToMolesOfPhase_(elements_vec)

    def convert_elements_to_mass_of_phase(self, elements):
        elements_vec = array_to_double_vector(elements)
        return self._phase_obj.convertElementsToMassOfPhase_(elements_vec)

    def set_ref_state(self, Tr=298.15, Pr=1.0, Trl=298.15):
        self._phase_obj.setTr_(Tr)
        self._phase_obj.setPr_(Pr)
        self._phase_obj.setTrl_(Trl)
        pass

    def get_gibbs_energy( self, T_a, P_a ):
        T_a, P_a = fill_array(T_a, P_a)

        G_a = np.zeros(T_a.shape)
        for ind,(T,P) in enumerate(zip(T_a,P_a)):
            G_a[ind] = self._phase_obj.getGibbsFreeEnergyFromT_andP_( T, P )

        return G_a

    def get_enthalpy( self, T_a, P_a ):
        T_a, P_a = fill_array(T_a, P_a)

        H_a = np.zeros(T_a.shape)
        for ind,(T,P) in enumerate(zip(T_a,P_a)):
            H_a[ind] = self._phase_obj.getEnthalpyFromT_andP_( T, P )

        return H_a

    def get_entropy( self, T_a, P_a ):
        T_a, P_a = fill_array(T_a, P_a)

        S_a = np.zeros(T_a.shape)
        for ind,(T,P) in enumerate(zip(T_a,P_a)):
            S_a[ind] = self._phase_obj.getEntropyFromT_andP_( T, P )

        return S_a

    def get_heat_capacity( self, T_a, P_a ):
        T_a, P_a = fill_array(T_a, P_a)

        Cp_a = np.zeros(T_a.shape)
        for ind,(T,P) in enumerate(zip(T_a,P_a)):
            Cp_a[ind] = self._phase_obj.getHeatCapacityFromT_andP_( T, P )

        return Cp_a

    def get_dCp_dT( self, T_a, P_a ):
        T_a, P_a = fill_array(T_a, P_a)

        dCp_dT_a = np.zeros(T_a.shape)
        for ind,(T,P) in enumerate(zip(T_a,P_a)):
            dCp_dT_a[ind] = self._phase_obj.getDcpDtFromT_andP_( T, P )

        return dCp_dT_a

    def get_volume( self, T_a, P_a ):
        T_a, P_a = fill_array(T_a, P_a)

        V_a = np.zeros(T_a.shape)
        for ind,(T,P) in enumerate(zip(T_a,P_a)):
            V_a[ind] = self._phase_obj.getVolumeFromT_andP_( T, P )

        return V_a

    def get_dV_dT( self, T_a, P_a ):
        T_a, P_a = fill_array(T_a, P_a)

        dV_dT_a = np.zeros(T_a.shape)
        for ind,(T,P) in enumerate(zip(T_a,P_a)):
            dV_dT_a[ind] = self._phase_obj.getDvDtFromT_andP_( T, P )

        return dV_dT_a

    def get_dV_dP( self, T_a, P_a ):
        T_a, P_a = fill_array(T_a, P_a)

        dV_dP_a = np.zeros(T_a.shape)
        for ind,(T,P) in enumerate(zip(T_a,P_a)):
            dV_dP_a[ind] = self._phase_obj.getDvDpFromT_andP_( T, P )

        return dV_dP_a

    def get_d2V_dT2( self, T_a, P_a ):
        T_a, P_a = fill_array(T_a, P_a)

        d2V_dT2_a = np.zeros(T_a.shape)
        for ind,(T,P) in enumerate(zip(T_a,P_a)):
            d2V_dT2_a[ind] = self._phase_obj.getD2vDt2FromT_andP_( T, P )

        return d2V_dT2_a

    def get_d2V_dTdP( self, T_a, P_a ):
        T_a, P_a = fill_array(T_a, P_a)

        d2V_dTdP_a = np.zeros(T_a.shape)
        for ind,(T,P) in enumerate(zip(T_a,P_a)):
            d2V_dTdP_a[ind] = self._phase_obj.getD2vDtDpFromT_andP_( T, P )

        return d2V_dTdP_a

    def get_d2V_dP2( self, T_a, P_a ):
        T_a, P_a = fill_array(T_a, P_a)

        d2V_dP2_a = np.zeros(T_a.shape)
        for ind,(T,P) in enumerate(zip(T_a,P_a)):
            d2V_dP2_a[ind] = self._phase_obj.getD2vDp2FromT_andP_( T, P )

        return d2V_dP2_a

    def get_chem_pot( self, T_a, P_a ):
        T_a, P_a = fill_array(T_a, P_a)

        chem_pot_a = np.zeros(T_a.shape)
        for ind,(T,P) in enumerate(zip(T_a,P_a)):
            chem_pot_a[ind] = self._phase_obj.getChemicalPotentialFromT_andP_( T, P )

        assert False, 'Is chem_pot an Objc vector? '+\
            'If so, the documentation needs to be updated'
        return chem_pot_a

    def get_param_units( self, param_names=[], all_params =False ):
        if all_params:
            param_names = self.get_param_names()
        elif type(param_names) == str:
            param_names = [param_names]

        units_l = []

        for ind,key in enumerate(param_names):
            unit = self._phase_obj.getUnitsForParameterName_(key)
            units_l.append(unit)

        return units_l

    def get_param_values( self, param_names=[], all_params =False ):
        if all_params:
            param_names = self.get_parameter_names()
        elif type(param_names) == str:
            param_names = [param_names]

        values_a = np.zeros(len(param_names))

        for ind,key in enumerate(param_names):
            value = self._phase_obj.getValueForParameterName_(key)
            values_a[ind] = value

        return values_a

    def set_param_values( self, param_names=[], param_values=[]  ):
        assert len(param_names)==len(param_values), \
            'param_names and param_values must have the same length'

        for name, value in zip(param_names, param_values):
            self._phase_obj.setParameterName_tovalue_(name,value)

        pass


    # def get_dG_dw( self, mol_a, T, P ):
    #     assert False, '
    #     mol_vec = array_to_double_vector(mol_a)
    #     print(mol_vec)
    #     dG_dw_vec = self._phase_obj.getDgDwFromMolesOfComponents_andT_andP_( mol_vec, T, P )
    #     print(dG_dw_vec)
    #     dG_dw_a = double_vector_to_array( dG_dw_vec )

    #     return dG_dw_a

    # def get_dchempot_dw
    # getChemicalPotentialDerivativesForParameterArray:(NSArray *)array usingMolesOfComponents:(double *)m andT:(double)t andP:

    # def get_dchempot_dw( self, T_a, P_a ):
    #     self._phase_obj.getChemicalPotentialDerivativesForParameter_usingMolesOfComponents_andT_andP_('S',m,T,P)
