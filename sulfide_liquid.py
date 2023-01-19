"""The sulfide_liquid package implements a coder-like solution phase
Python interface to the SulfLiq thermochemical model of Kress. See
https://gitlab.com/ENKI-portal/sulfliq for source code to the SulfLiq
package.

"""
import numpy as np
from thermoengine import chem
import SulfLiq

sl = SulfLiq.pySulfLiq()


def cy_SulfLiq_sulfide_liquid_calib_identifier():
    return "Version_1_0_0"
def cy_SulfLiq_sulfide_liquid_calib_name():
    return "Sulfide Liquid"

def cy_SulfLiq_sulfide_liquid_calib_formula(t, p, mol):
    result = ""
    for i in range(0, mol.shape[0]):
        result += sl.getCompFormula(i)
        result += str(round(mol[i], 3))
    return result
def cy_SulfLiq_sulfide_liquid_calib_conv_elm_to_moles(np_array):
    result = []
    iO  = np.where(chem.PERIODIC_ORDER == 'O')[0][0]
    iS  = np.where(chem.PERIODIC_ORDER == 'S')[0][0]
    iFe = np.where(chem.PERIODIC_ORDER == 'Fe')[0][0]
    iNi = np.where(chem.PERIODIC_ORDER == 'Ni')[0][0]
    iCu = np.where(chem.PERIODIC_ORDER == 'Cu')[0][0]
    result.append(np_array[iO])
    result.append(np_array[iS])
    result.append(np_array[iFe])
    result.append(np_array[iNi])
    result.append(np_array[iCu])
    return np.array(result)
def cy_SulfLiq_sulfide_liquid_calib_conv_elm_to_tot_moles(np_array):
    result = 0.0
    iO  = np.where(chem.PERIODIC_ORDER == 'O')[0][0]
    iS  = np.where(chem.PERIODIC_ORDER == 'S')[0][0]
    iFe = np.where(chem.PERIODIC_ORDER == 'Fe')[0][0]
    iNi = np.where(chem.PERIODIC_ORDER == 'Ni')[0][0]
    iCu = np.where(chem.PERIODIC_ORDER == 'Cu')[0][0]
    result += np_array[iO]
    result += np_array[iS]
    result += np_array[iFe]
    result += np_array[iNi]
    result += np_array[iCu]
    return result
def cy_SulfLiq_sulfide_liquid_calib_conv_elm_to_tot_grams(np_array):
    result = 0.0
    ind = np.where(chem.PERIODIC_ORDER == 'O')[0][0]
    result += np_array[ind]*chem.PERIODIC_WEIGHTS[ind]
    ind = np.where(chem.PERIODIC_ORDER == 'S')[0][0]
    result += np_array[ind]*chem.PERIODIC_WEIGHTS[ind]
    ind = np.where(chem.PERIODIC_ORDER == 'Fe')[0][0]
    result += np_array[ind]*chem.PERIODIC_WEIGHTS[ind]
    ind = np.where(chem.PERIODIC_ORDER == 'Ni')[0][0]
    result += np_array[ind]*chem.PERIODIC_WEIGHTS[ind]
    ind = np.where(chem.PERIODIC_ORDER == 'Cu')[0][0]
    result += np_array[ind]*chem.PERIODIC_WEIGHTS[ind]
    return result
def cy_SulfLiq_sulfide_liquid_calib_conv_moles_to_tot_moles(np_array):
    return np.sum(np_array)
def cy_SulfLiq_sulfide_liquid_calib_conv_moles_to_mole_frac(np_array):
    total = np.sum(np_array)
    return np_array/total
def cy_SulfLiq_sulfide_liquid_calib_conv_moles_to_elm(np_array):
    result = np.zeros(106)
    iO  = np.where(chem.PERIODIC_ORDER == 'O')[0][0]
    iS  = np.where(chem.PERIODIC_ORDER == 'S')[0][0]
    iFe = np.where(chem.PERIODIC_ORDER == 'Fe')[0][0]
    iNi = np.where(chem.PERIODIC_ORDER == 'Ni')[0][0]
    iCu = np.where(chem.PERIODIC_ORDER == 'Cu')[0][0]
    result[iO] = np_array[0]
    result[iS] = np_array[1]
    result[iFe] = np_array[2]
    result[iNi] = np_array[3]
    result[iCu] = np_array[4]
    return result
def cy_SulfLiq_sulfide_liquid_calib_test_moles(np_array):
    result = True
    for x in np_array:
        if x < 0.0:
            result = False
    return result

def cy_SulfLiq_sulfide_liquid_calib_endmember_number():
    return 5
def cy_SulfLiq_sulfide_liquid_calib_endmember_name(index):
    assert index in range(0,5), "index out of range"
    if index == 0:
        return "O"
    elif index == 1:
        return "S"
    elif index == 2:
        return "Fe"
    elif index == 3:
        return "Ni"
    else:
        return "Cu"
def cy_SulfLiq_sulfide_liquid_calib_endmember_formula(index):
    assert index in range(0,5), "index out of range"
    if index == 0:
        return "O"
    elif index == 1:
        return "S"
    elif index == 2:
        return "Fe"
    elif index == 3:
        return "Ni"
    else:
        return "Cu"
def cy_SulfLiq_sulfide_liquid_calib_endmember_mw(index):
    assert index in range(0,5), "index out of range"
    iO  = np.where(chem.PERIODIC_ORDER == 'O')[0][0]
    iS  = np.where(chem.PERIODIC_ORDER == 'S')[0][0]
    iFe = np.where(chem.PERIODIC_ORDER == 'Fe')[0][0]
    iNi = np.where(chem.PERIODIC_ORDER == 'Ni')[0][0]
    iCu = np.where(chem.PERIODIC_ORDER == 'Cu')[0][0]
    if index == 0:
        return chem.PERIODIC_WEIGHTS[iO]
    elif index == 1:
        return chem.PERIODIC_WEIGHTS[iS]
    elif index == 2:
        return chem.PERIODIC_WEIGHTS[iFe]
    elif index == 3:
        return chem.PERIODIC_WEIGHTS[iNi]
    else:
        return chem.PERIODIC_WEIGHTS[iCu]
def cy_SulfLiq_sulfide_liquid_calib_endmember_elements(index):
    assert index in range(0,5), "index out of range"
    result = np.zeros(106)
    iO  = np.where(chem.PERIODIC_ORDER == 'O')[0][0]
    iS  = np.where(chem.PERIODIC_ORDER == 'S')[0][0]
    iFe = np.where(chem.PERIODIC_ORDER == 'Fe')[0][0]
    iNi = np.where(chem.PERIODIC_ORDER == 'Ni')[0][0]
    iCu = np.where(chem.PERIODIC_ORDER == 'Cu')[0][0]
    if index == 0:
        result[iO] = 1.0
    elif index == 1:
        result[iS] = 1.0
    elif index == 2:
        result[iFe] = 1.0
    elif index == 3:
        result[iNi] = 1.0
    else:
        result[iCu] = 1.0
    return result
def cy_SulfLiq_sulfide_liquid_calib_species_number():
    return sl.getNspec()
def cy_SulfLiq_sulfide_liquid_calib_species_name(index):
    assert index in range(0,sl.getNspec()), "index out of range"
    result = sl.getSpecFormula(index)
    return result
def cy_SulfLiq_sulfide_liquid_calib_species_formula(index):
    assert index in range(0,sl.getNspec()), "index out of range"
    result = sl.getSpecFormula(index)
    return result
def cy_SulfLiq_sulfide_liquid_calib_species_mw(index):
    assert index in range(0,15), "index out of range"
    iO  = np.where(chem.PERIODIC_ORDER == 'O')[0][0]
    iS  = np.where(chem.PERIODIC_ORDER == 'S')[0][0]
    iFe = np.where(chem.PERIODIC_ORDER == 'Fe')[0][0]
    iNi = np.where(chem.PERIODIC_ORDER == 'Ni')[0][0]
    iCu = np.where(chem.PERIODIC_ORDER == 'Cu')[0][0]
    if index == 0:
        return chem.PERIODIC_WEIGHTS[iO]
    elif index == 1:
        return chem.PERIODIC_WEIGHTS[iS]
    elif index == 2:
        return chem.PERIODIC_WEIGHTS[iFe]
    elif index == 3:
        return chem.PERIODIC_WEIGHTS[iNi]
    elif index == 4:
        return chem.PERIODIC_WEIGHTS[iCu]
    elif index == 5:
        return chem.PERIODIC_WEIGHTS[iFe] + chem.PERIODIC_WEIGHTS[iO]
    elif index == 6:
        return chem.PERIODIC_WEIGHTS[iFe] + 1.5*chem.PERIODIC_WEIGHTS[iO]
    elif index == 7:
        return chem.PERIODIC_WEIGHTS[iNi] + chem.PERIODIC_WEIGHTS[iO]
    elif index == 8:
        return chem.PERIODIC_WEIGHTS[iFe] + chem.PERIODIC_WEIGHTS[iS]
    elif index == 9:
        return chem.PERIODIC_WEIGHTS[iNi] + chem.PERIODIC_WEIGHTS[iS]
    elif index == 10:
        return chem.PERIODIC_WEIGHTS[iCu] + chem.PERIODIC_WEIGHTS[iS]/2.0
    elif index == 11:
        return chem.PERIODIC_WEIGHTS[iFe] + chem.PERIODIC_WEIGHTS[iO] + chem.PERIODIC_WEIGHTS[iS]
    elif index == 12:
        return chem.PERIODIC_WEIGHTS[iNi]/4.0 + chem.PERIODIC_WEIGHTS[iO] + chem.PERIODIC_WEIGHTS[iS]/4.0
    elif index == 13:
        return chem.PERIODIC_WEIGHTS[iCu] + chem.PERIODIC_WEIGHTS[iS]
    elif index == 14:
        return chem.PERIODIC_WEIGHTS[iCu] + chem.PERIODIC_WEIGHTS[iO]
def cy_SulfLiq_sulfide_liquid_calib_species_elements(index):
    assert index in range(0,15), "index out of range"
    iO  = np.where(chem.PERIODIC_ORDER == 'O')[0][0]
    iS  = np.where(chem.PERIODIC_ORDER == 'S')[0][0]
    iFe = np.where(chem.PERIODIC_ORDER == 'Fe')[0][0]
    iNi = np.where(chem.PERIODIC_ORDER == 'Ni')[0][0]
    iCu = np.where(chem.PERIODIC_ORDER == 'Cu')[0][0]
    result = np.zeros(106)
    if index == 0:
        result[iO] = 1.0
    elif index == 1:
        result[iS] = 1.0
    elif index == 2:
        result[iFe] = 1.0
    elif index == 3:
        result[iNi] = 1.0
    elif index == 4:
        result[iCu] = 1.0
    elif index == 5:
        result[iFe] = 1.0
        result[iO] = 1.0
    elif index == 6:
        result[iFe] = 1.0
        result[iO] = 1.5
    elif index == 7:
        result[iNi] = 1.0
        result[iO] = 1.0
    elif index == 8:
        result[iFe] = 1.0
        result[iS] = 1.0
    elif index == 9:
        result[iNi] = 1.0
        result[iS] = 1.0
    elif index == 10:
        result[iCu] = 1.0
        result[iS] = 0.5
    elif index == 11:
        result[iFe] = 1.0
        result[iO] = 1.0
        result[iS] = 1.0
    elif index == 12:
        result[iNi] = 0.25
        result[iO] = 1.0
        result[iS] = 0.25
    elif index == 13:
        result[iCu] = 1.0
        result[iS] = 1.0
    elif index == 14:
        result[iCu] = 1.0
        result[iO] = 1.0
    return result

def cy_SulfLiq_sulfide_liquid_calib_endmember_mu0(index, t, p):
    assert index in range(0,5), "index out of range"
    sl.setTK(t)
    sl.setPa(p*1.0e5)
    return sl.getMu0(index)
def cy_SulfLiq_sulfide_liquid_calib_endmember_dmu0dT(index, t, p):
    assert index in range(0,5), "index out of range"
    eps = t*np.sqrt(np.finfo(float).eps)
    sl.setTK(t)
    sl.setPa(p*1.0e5)
    f = sl.getMu0(index)
    sl.setTK(t+eps)
    ff = sl.getMu0(index)
    return (ff-f)/eps
def cy_SulfLiq_sulfide_liquid_calib_endmember_dmu0dP(index, t, p):
    assert index in range(0,5), "index out of range"
    eps = p*1.0e5*np.sqrt(np.finfo(float).eps)
    sl.setTK(t)
    sl.setPa(p*1.0e5)
    f = sl.getMu0(index)
    sl.setPa(p*1.0e5+eps)
    ff = sl.getMu0(index)
    return (ff-f)/eps
def cy_SulfLiq_sulfide_liquid_calib_endmember_d2mu0dT2(index, t, p):
    assert index in range(0,5), "index out of range"
    eps = t*np.sqrt(np.finfo(float).eps)
    sl.setTK(t)
    sl.setPa(p*1.0e5)
    f = sl.getMu0(index)
    sl.setTK(t+eps)
    ff = sl.getMu0(index)
    sl.setTK(t-eps)
    fb = sl.getMu0(index)
    return (ff-2.0*f+fb)/eps/eps
def cy_SulfLiq_sulfide_liquid_calib_endmember_d2mu0dTdP(index, t, p):
    assert index in range(0,5), "index out of range"
    epsT = t*np.sqrt(np.finfo(float).eps)
    epsP = p*1.0e5*np.sqrt(np.finfo(float).eps)
    sl.setTK(t)
    sl.setPa(p*1.0e5)
    f = sl.getMu0(index)
    sl.setTK(t+epsT)
    sl.setPa(p*1.0e5+epsP)
    ff = sl.getMu0(index)
    sl.setTK(t-epsT)
    sl.setPa(p*1.0e5-epsP)
    fb = sl.getMu0(index)
    return (ff-2.0*f+fb)/epsT/epsP
def cy_SulfLiq_sulfide_liquid_calib_endmember_d2mu0dP2(index, t, p):
    assert index in range(0,5), "index out of range"
    eps = p*1.0e5*np.sqrt(np.finfo(float).eps)
    sl.setTK(t)
    sl.setPa(p*1.0e5)
    f = sl.getMu0(index)
    sl.setPa(p*1.0e5+eps)
    ff = sl.getMu0(index)
    sl.setPa(p*1.0e5-eps)
    fb = sl.getMu0(index)
    return (ff-2.0*f+fb)/eps/eps

def cy_SulfLiq_sulfide_liquid_calib_endmember_d3mu0dT3(index, t, p):
    assert index in range(0,5), "index out of range"
    return 0.0
def cy_SulfLiq_sulfide_liquid_calib_endmember_d3mu0dT2dP(index, t, p):
    assert index in range(0,5), "index out of range"
    return 0.0
def cy_SulfLiq_sulfide_liquid_calib_endmember_d3mu0dTdP2(index, t, p):
    assert index in range(0,5), "index out of range"
    return 0.0
def cy_SulfLiq_sulfide_liquid_calib_endmember_d3mu0dP3(index, t, p):
    assert index in range(0,5), "index out of range"
    return 0.0

def cy_SulfLiq_sulfide_liquid_calib_g(t, p, mol):
    sl.setComps(mol.tolist())
    sl.setTK(t)
    sl.setPa(p*1.0e5)
    result = sl.getGibbs()
    return result
def cy_SulfLiq_sulfide_liquid_calib_dgdt(t, p, mol):
    sl.setComps(mol.tolist())
    sl.setTK(t)
    sl.setPa(p*1.0e5)
    result = -sl.getEntropy()
    return result
def cy_SulfLiq_sulfide_liquid_calib_dgdp(t, p, mol):
    sl.setComps(mol.tolist())
    sl.setTK(t)
    sl.setPa(p*1.0e5)
    result = sl.getVolume()*1.0e5 # m^3 -> J/bar
    return result
def cy_SulfLiq_sulfide_liquid_calib_d2gdt2(t, p, mol):
    sl.setComps(mol.tolist())
    sl.setTK(t)
    sl.setPa(p*1.0e5)
    result = -sl.getCp()/t
    return result
def cy_SulfLiq_sulfide_liquid_calib_d2gdtdp(t, p, mol):
    sl.setComps(mol.tolist())
    sl.setTK(t)
    sl.setPa(p*1.0e5)
    result = sl.getdVdT()*1.0e5 # m^3/K -> J/bar/K
    return result
def cy_SulfLiq_sulfide_liquid_calib_d2gdp2(t, p, mol):
    sl.setComps(mol.tolist())
    sl.setTK(t)
    sl.setPa(p*1.0e5)
    result = sl.getdVdP()*1.0e10 # m^3/Pa -> J/bar^2
    return result
def cy_SulfLiq_sulfide_liquid_calib_d3gdt3(t, p, mol):
    sl.setComps(mol.tolist())
    sl.setTK(t)
    sl.setPa(p*1.0e5)
    result = sl.getCp()/t/t - sl.getdCpdT()/t
    return result
def cy_SulfLiq_sulfide_liquid_calib_d3gdt2dp(t, p, mol):
    sl.setComps(mol.tolist())
    sl.setTK(t)
    sl.setPa(p*1.0e5)
    result = sl.getd2VdT2() *1.0e5 # m^3/K^2 -> J/bar/K^2
    return result
def cy_SulfLiq_sulfide_liquid_calib_d3gdtdp2(t, p, mol):
    sl.setComps(mol.tolist())
    sl.setTK(t)
    sl.setPa(p*1.0e5)
    result = sl.getd2VdTdP() *1.0e10 # m^3/K/Pa -> J/K/bar^2
    return result
def cy_SulfLiq_sulfide_liquid_calib_d3gdp3(t, p, mol):
    sl.setComps(mol.tolist())
    sl.setTK(t)
    sl.setPa(p*1.0e5)
    result = sl.getd2VdP2() *1.0e15 # m^3/Pa^2 -> J/bar^3
    return result
def cy_SulfLiq_sulfide_liquid_calib_s(t, p, mol):
    sl.setComps(mol.tolist())
    sl.setTK(t)
    sl.setPa(p*1.0e5)
    result = sl.getEntropy()
    return result
def cy_SulfLiq_sulfide_liquid_calib_v(t, p, mol):
    sl.setComps(mol.tolist())
    sl.setTK(t)
    sl.setPa(p*1.0e5)
    result = sl.getVolume() *1.0e5 # m^3 -> J/bar
    return result
def cy_SulfLiq_sulfide_liquid_calib_cv(t, p, mol):
    sl.setComps(mol.tolist())
    sl.setTK(t)
    sl.setPa(p*1.0e5)
    dvdt = sl.getdVdT()
    result = sl.getCp() + t*dvdt*dvdt/sl.getdVdP()
    return result
def cy_SulfLiq_sulfide_liquid_calib_cp(t, p, mol):
    sl.setComps(mol.tolist())
    sl.setTK(t)
    sl.setPa(p*1.0e5)
    result = sl.getCp()
    return result
def cy_SulfLiq_sulfide_liquid_calib_dcpdt(t, p, mol):
    sl.setComps(mol.tolist())
    sl.setTK(t)
    sl.setPa(p*1.0e5)
    result = sl.getdCpdT()
    return result
def cy_SulfLiq_sulfide_liquid_calib_alpha(t, p, mol):
    sl.setComps(mol.tolist())
    sl.setTK(t)
    sl.setPa(p*1.0e5)
    result = sl.getdVdT()/sl.getVolume()
    return result
def cy_SulfLiq_sulfide_liquid_calib_beta(t, p, mol):
    sl.setComps(mol.tolist())
    sl.setTK(t)
    sl.setPa(p*1.0e5)
    result = -sl.getdVdP()*1.0e5/sl.getVolume() # 1/Pa -> 1/bar
    return result
def cy_SulfLiq_sulfide_liquid_calib_K(t, p, mol):
    sl.setComps(mol.tolist())
    sl.setTK(t)
    sl.setPa(p*1.0e5)
    result = -sl.getVolume()/sl.getdVdP()/1.0e5 # Pa -> bar
    return result
def cy_SulfLiq_sulfide_liquid_calib_Kp(t, p, mol):
    sl.setComps(mol.tolist())
    sl.setTK(t)
    sl.setPa(p*1.0e5)
    dvdp = sl.getdVdP()
    result = sl.getVolume()*sl.getd2VdP2()/dvdp/dvdp - 1.0
    return result

def cy_SulfLiq_sulfide_liquid_calib_dgdn(t, p, mol):
    sl.setComps(mol.tolist())
    sl.setTK(t)
    sl.setPa(p*1.0e5)
    result = sl.getdGdm()
    return np.array(result)
def cy_SulfLiq_sulfide_liquid_calib_d2gdndt(t, p, mol):
    sl.setComps(mol.tolist())
    sl.setTK(t)
    sl.setPa(p*1.0e5)
    result = -np.array(sl.getdSdm())
    return result
def cy_SulfLiq_sulfide_liquid_calib_d2gdndp(t, p, mol):
    sl.setComps(mol.tolist())
    sl.setTK(t)
    sl.setPa(p*1.0e5)
    result = np.array(sl.getdVdm())*1.0e5 # m^3 -> J/bar
    return result
def cy_SulfLiq_sulfide_liquid_calib_d3gdndt2(t, p, mol):
    sl.setComps(mol.tolist())
    sl.setTK(t)
    sl.setPa(p*1.0e5)
    result = -np.array(sl.getdCpdm())/t
    return result
def cy_SulfLiq_sulfide_liquid_calib_d3gdndtdp(t, p, mol):
    sl.setComps(mol.tolist())
    sl.setTK(t)
    sl.setPa(p*1.0e5)
    result = np.array(sl.getd2VdmdT())*1.0e5 # m^3/K -> J/bar/K
    return result
def cy_SulfLiq_sulfide_liquid_calib_d3gdndp2(t, p, mol):
    sl.setComps(mol.tolist())
    sl.setTK(t)
    sl.setPa(p*1.0e5)
    result = np.array(sl.getd2VdmdP())*1.0e10 # m^3/Pa -> J/bar^2
    return result
def cy_SulfLiq_sulfide_liquid_calib_d4gdndt3(t, p, mol):
    return np.zeros(mol.shape)
def cy_SulfLiq_sulfide_liquid_calib_d4gdndt2dp(t, p, mol):
    return np.zeros(mol.shape)
def cy_SulfLiq_sulfide_liquid_calib_d4gdndtdp2(t, p, mol):
    return np.zeros(mol.shape)
def cy_SulfLiq_sulfide_liquid_calib_d4gdndp3(t, p, mol):
    return np.zeros(mol.shape)

def cy_SulfLiq_sulfide_liquid_calib_d2gdn2(t, p, mol):
    sl.setComps(mol.tolist())
    sl.setTK(t)
    sl.setPa(p*1.0e5)
    result = np.array(sl.getd2Gdm2())
    return result[np.triu_indices(mol.shape[0])]
def cy_SulfLiq_sulfide_liquid_calib_d3gdn2dt(t, p, mol):
    sl.setComps(mol.tolist())
    sl.setTK(t)
    sl.setPa(p*1.0e5)
    result = -np.array(sl.getd2Sdm2())
    return result[np.triu_indices(mol.shape[0])]
def cy_SulfLiq_sulfide_liquid_calib_d3gdn2dp(t, p, mol):
    sl.setComps(mol.tolist())
    sl.setTK(t)
    sl.setPa(p*1.0e5)
    result = np.array(sl.getd2Vdm2())*1.0e5 # m^3-> J/bar
    return result[np.triu_indices(mol.shape[0])]

def cy_SulfLiq_sulfide_liquid_calib_d4gdn2dt2(t, p, mol):
    n = mol.shape[0]
    return np.zeros(int(n*(n-1)/2+n))
def cy_SulfLiq_sulfide_liquid_calib_d4gdn2dtdp(t, p, mol):
    n = mol.shape[0]
    return np.zeros(int(n*(n-1)/2+n))
def cy_SulfLiq_sulfide_liquid_calib_d4gdn2dp2(t, p, mol):
    n = mol.shape[0]
    return np.zeros(int(n*(n-1)/2+n))
def cy_SulfLiq_sulfide_liquid_calib_d5gdn2dt3(t, p, mol):
    n = mol.shape[0]
    return np.zeros(int(n*(n-1)/2+n))
def cy_SulfLiq_sulfide_liquid_calib_d5gdn2dt2dp(t, p, mol):
    n = mol.shape[0]
    return np.zeros(int(n*(n-1)/2+n))
def cy_SulfLiq_sulfide_liquid_calib_d5gdn2dtdp2(t, p, mol):
    n = mol.shape[0]
    return np.zeros(int(n*(n-1)/2+n))
def cy_SulfLiq_sulfide_liquid_calib_d5gdn2dp3(t, p, mol):
    n = mol.shape[0]
    return np.zeros(int(n*(n-1)/2+n))

def cy_SulfLiq_sulfide_liquid_calib_d3gdn3(t, p, mol):
    n = mol.shape[0]
    return np.zeros(int(n*(n+1)*(n+2)/6))
def cy_SulfLiq_sulfide_liquid_calib_d4gdn3dt(t, p, mol):
    n = mol.shape[0]
    return np.zeros(int(n*(n+1)*(n+2)/6))
def cy_SulfLiq_sulfide_liquid_calib_d4gdn3dp(t, p, mol):
    n = mol.shape[0]
    return np.zeros(int(n*(n+1)*(n+2)/6))
def cy_SulfLiq_sulfide_liquid_calib_d5gdn3dt2(t, p, mol):
    n = mol.shape[0]
    return np.zeros(int(n*(n+1)*(n+2)/6))
def cy_SulfLiq_sulfide_liquid_calib_d5gdn3dtdp(t, p, mol):
    n = mol.shape[0]
    return np.zeros(int(n*(n+1)*(n+2)/6))
def cy_SulfLiq_sulfide_liquid_calib_d5gdn3dp2(t, p, mol):
    n = mol.shape[0]
    return np.zeros(int(n*(n+1)*(n+2)/6))
def cy_SulfLiq_sulfide_liquid_calib_d6gdn3dt3(t, p, mol):
    n = mol.shape[0]
    return np.zeros(int(n*(n+1)*(n+2)/6))
def cy_SulfLiq_sulfide_liquid_calib_d6gdn3dt2dp(t, p, mol):
    n = mol.shape[0]
    return np.zeros(int(n*(n+1)*(n+2)/6))
def cy_SulfLiq_sulfide_liquid_calib_d6gdn3dtdp2(t, p, mol):
    n = mol.shape[0]
    return np.zeros(int(n*(n+1)*(n+2)/6))
def cy_SulfLiq_sulfide_liquid_calib_d6gdn3dp3(t, p, mol):
    n = mol.shape[0]
    return np.zeros(int(n*(n+1)*(n+2)/6))

def cy_SulfLiq_sulfide_liquid_get_param_number():
    return 0
def cy_SulfLiq_sulfide_liquid_get_param_names():
    return [""]
def cy_SulfLiq_sulfide_liquid_get_param_units():
    return [""]
def cy_SulfLiq_sulfide_liquid_get_param_values():
    return np.zeros(1)
def cy_SulfLiq_sulfide_liquid_set_param_values(np_array):
    return True
def cy_SulfLiq_sulfide_liquid_get_param_value(index):
    return True
def cy_SulfLiq_sulfide_liquid_set_param_value(index, value):
    return True

def cy_SulfLiq_sulfide_liquid_dparam_g(t, p, mol, index):
    return np.zeros(1)
def cy_SulfLiq_sulfide_liquid_dparam_dgdt(t, p, mol, index):
    return np.zeros(1)
def cy_SulfLiq_sulfide_liquid_dparam_dgdp(t, p, mol, index):
    return np.zeros(1)
def cy_SulfLiq_sulfide_liquid_dparam_d2gdt2(t, p, mol, index):
    return np.zeros(1)
def cy_SulfLiq_sulfide_liquid_dparam_d2gdtdp(t, p, mol, index):
    return np.zeros(1)
def cy_SulfLiq_sulfide_liquid_dparam_d2gdp2(t, p, mol, index):
    return np.zeros(1)
def cy_SulfLiq_sulfide_liquid_dparam_d3gdt3(t, p, mol, index):
    return np.zeros(1)
def cy_SulfLiq_sulfide_liquid_dparam_d3gdt2dp(t, p, mol, index):
    return np.zeros(1)
def cy_SulfLiq_sulfide_liquid_dparam_d3gdtdp2(t, p, mol, index):
    return np.zeros(1)
def cy_SulfLiq_sulfide_liquid_dparam_d3gdp3(t, p, mol, index):
    return np.zeros(1)

def cy_SulfLiq_sulfide_liquid_dparam_dgdn(t, p, mol, index):
    return np.zeros(len(mol))
