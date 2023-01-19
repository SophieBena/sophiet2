from __future__ import annotations  # Enable Python 4 type hints in Python 3
import abc
import collections
import copy
from functools import wraps
from typing import Type, List, Optional, Tuple
from typing import NewType, TypeVar, NamedTuple
from collections import namedtuple

import numpy as np
import pandas as pd
from numpy.typing import NDArray
import thermoengine.phases
from scipy import optimize as opt

from thermoengine import phases as phs
from thermoengine import chem, chemistry
from thermoengine import Phase, ElemMolComp

from thermoengine.core import UnorderedList

from dataclasses import dataclass

__all__ = ['SampleMaker']


# Phase = NewType('Phase', thermoengine.phases.Phase)
class SampleSummary:
    """Summary report for composition of phase sample."""

    def __init__(self, phase_sample: _PhaseSample):
        self.init_title(phase_sample)
        self.init_mol_comp(phase_sample)

        self.overview = '\n'.join([self.title] + self.mol_comp)

    def init_title(self, samp: _PhaseSample):
        self.title = samp.phase_name

    def init_mol_comp(self, samp: _PhaseSample):
        mol_comp = []
        mol_comp.append('| ----------- endmember mol frac ------------ |')
        for name, formula, X in zip(
                samp.endmember_names, samp.endmember_formulas, samp.X):
            formula = '[' + formula + ']'
            mol_comp.append(f'|{name:>15} {formula:17} {X:.4f} mol |')

        self.mol_comp = mol_comp

    def __str__(self):
        return self.overview


class _PhaseSample(abc.ABC):
    phase_model: phs.Phase
    X: np.ndarray
    T: float
    P: float
    fixed_comp: bool

    comp: chemistry.ElemMolComp
    G: float
    aff: float
    amount: float
    metastable: bool
    atom_num: int
    name: str
    endmember_num: int

    def __init__(self, phase_model: phs.Phase,
                 X: List[float] = None, T: float = 1000,
                 P: float = 1, amount: float = 1, aff: float = 0,
                 metastable: bool = False, XTOL: float = 1e-6):
        self._init_phase_model(phase_model)
        self._init_conditions(T, P, X)
        self.G = self.gibbs_energy()
        self.aff = aff
        self.amount = amount
        self.metastable = metastable
        self._XTOL = XTOL

    def _init_conditions(self, T, P, X):
        self.T = T
        self.P = P
        if X is None:
            self._X = np.eye(self.phase_model.endmember_num)[0]
        else:
            X = np.array(X)
            assert X.size == self.endmember_num, \
                'dimensions of X do not match number of endmembers.'
            self._X = X
        self.comp = self.get_molar_elem_comp()

    def update_conditions(self, T, P):
        self.T = T
        self.P = P
        self.G = self.gibbs_energy()
        self.aff = 0
        self.amount = 1

    @property
    def phase_model(self) -> phs.Phase:
        return self._phase_model

    def _init_phase_model(self, model):
        self._phase_model = model
        self.is_solution = model.phase_type == 'solution'
        self.fixed_comp = not self.is_solution
        self.name = model.abbrev
        self.endmember_num = model.endmember_num
        atom_num = model.props['atom_num']
        if not np.isscalar(atom_num):
            atom_num = atom_num[0]
        self.atom_num = atom_num

    @property
    def endmember_names(self) -> List[str]:
        return self.phase_model.endmember_names

    @property
    def endmember_formulas(self) -> List[str]:
        return self.phase_model.props['formula']

    @property
    def endmember_elem_comp(self) -> pd.DataFrame:
        return self.phase_model.endmember_elem_comp

    @property
    def phase_name(self) -> str:
        return self.phase_model.phase_name

    @property
    def X(self):
        return self._X

    @property
    def summary(self) -> SampleSummary:
        return SampleSummary(self)

    def get_molar_elem_comp(self):
        elem_comp = self.phase_model.props['element_comp']
        elem_mask = np.any(elem_comp > 0, axis=0)
        mol_elems = self.X.dot(elem_comp[:, elem_mask])
        elems = chem.PERIODIC_ORDER[elem_mask]
        return chemistry.ElemMolComp(**dict(zip(elems, mol_elems)))

    def modify_sample_conditions(self, T=None, P=None, X=None, amount=None,
                                 metastable=None):
        if T is None:
            T = self.T
        if P is None:
            P = self.P
        if X is None:
            X = self.X
        if amount is None:
            amount = self.amount
        if metastable is None:
            metastable = self.metastable

        return SampleMaker.get_sample(self.phase_model, X=X, T=T, P=P,
                                      amount=amount, metastable=metastable)

    def optionalX(func):
        @wraps(func)
        def check_input(*args, **kwargs):
            self = args[0]
            X = kwargs['X'] if ('X' in kwargs) else None
            if X is None:
                X = self.X

            result = func(self, X=X)
            return result

        return check_input

    @optionalX
    def gibbs_energy(self, X: Optional[np.ndarray] = None):
        return self._calc_gibbs_energy(X) / self.atom_num

    @optionalX
    def chem_potential(self, X=None):
        return self._calc_chem_potential(X) / self.atom_num

    @abc.abstractmethod
    def _calc_gibbs_energy(self, X):
        pass

    @abc.abstractmethod
    def _calc_chem_potential(self, X):
        pass

    ## variable_comp
    def equilibrate_by_chemical_exchange(self, chem_potential):
        self.chem_potential_ref = chem_potential
        self.Xref = self.X
        result = opt.minimize(self.affinity_perturbation,
                              x0=np.zeros(self.X.size))

        return SampleMaker.get_sample(
            self.phase_model, X=self._calc_perturbed_comp(result.x),
            T=self.T, P=self.P)

    def _calc_perturbed_comp(self, dlogX):
        X = np.exp(dlogX) * self.Xref
        return X / X.sum()

    def affinity_perturbation(self, dlogX):
        X = self._calc_perturbed_comp(dlogX)
        return self.affinity(X)

    def affinity(self, X):
        dmu = self.chem_potential(X=X) - self.chem_potential_ref
        return np.dot(dmu, X)

    def __eq__(self, other):
        same_model = self.phase_model == other.phase_model
        same_fixed_comp = self.fixed_comp == other.fixed_comp
        same_model_config = same_model and same_fixed_comp

        same_T = np.allclose(self.T, other.T)
        same_P = np.allclose(self.P, other.P)
        same_X = SampleMesh.mol_comps_equiv(self.X, other.X,
                                            TOL=self._XTOL)

        same_conditions = same_T and same_P and same_X

        return same_model_config and same_conditions

    def __lt__(self, other):
        return self.name < other.name


class SampleMaker:
    XTOL = 1e-8

    @classmethod
    def get_sample(cls, phase_model: Phase,
                   X: List[float] = None, T: float = 1000,
                   P: float = 1, amount: float = 1,
                   aff: float = 0, metastable: bool = False) -> _PhaseSample:
        """
        Create PhaseSample with given composition and environmental conditions.

        Parameters
        ----------
        phase_model
            Phase model providing gibbs surface.
        X
            Molar composition of phase in terms of endmember components.
        T
            Temperature in Kelvin.
        P
            Pressure in bars.

        Returns
        -------
        PhaseSample
            Phase Sample with desired composition and environmental variables.

        """

        if phase_model.phase_type == 'solution':
            return _SolutionPhaseSample(phase_model=phase_model, X=X, T=T, P=P,
                                        amount=amount, aff=aff, metastable=metastable)
        else:
            return _PurePhaseSample(phase_model=phase_model, X=X, T=T, P=P,
                                    amount=amount, aff=aff, metastable=metastable)

    @classmethod
    def copy_sample(cls, sample: _PhaseSample) -> _PhaseSample:
        if type(sample) is _FixedCompPhaseSample:
            return cls.get_fixed_comp_sample(sample.phase_model, sample.X,
                                             sample.T, sample.P, sample.amount,
                                             sample.aff, sample.metastable)
        else:
            return cls.get_sample(sample.phase_model, sample.X,
                                  sample.T, sample.P, sample.amount,
                                  sample.aff, sample.metastable)

    @classmethod
    def get_fixed_comp_sample(cls, phase_model: phs.Phase,
                              X: List[float] = None, T: float = 1000,
                              P: float = 1, amount: float = 1,
                              aff: float = 0, metastable: bool = False) -> _PhaseSample:

        if phase_model.phase_type == 'solution':
            return _FixedCompPhaseSample(phase_model=phase_model, X=X, T=T, P=P,
                                         amount=amount, aff=aff, metastable=metastable)
        else:
            return _PurePhaseSample(phase_model=phase_model, X=X, T=T, P=P,
                                    amount=amount, aff=aff, metastable=metastable)

    @classmethod
    def get_sample_endmembers(cls, phase_model: phs.Phase, T: float = 1000,
                              P: float = 1, amount: float = 1) -> List[_PhaseSample]:

        X_endmems = cls._get_phase_endmember_comps(phase_model)
        return cls._get_sample_set_from_X_grid(phase_model, X_endmems, T, P, amount)

    @classmethod
    def get_sample_grid(cls, phase_model: phs.Phase, grid_spacing: float = 0.2,
                        T: float = 1000, P: float = 1, amount: float = 1) -> List[_PhaseSample]:

        Nendmem = phase_model.endmember_num
        Xgrid = cls.build_comp_grid(Nendmem, grid_spacing)
        return cls._get_fixed_sample_set_from_X_grid(phase_model, Xgrid, T, P, amount)

    @classmethod
    def _get_fixed_sample_set_from_X_grid(cls, phase_model, Xgrid, T, P, amount) -> List[_PhaseSample]:
        sample_set = []
        [sample_set.append(SampleMaker.get_fixed_comp_sample(
            phase_model, X=X, T=T, P=P, amount=amount))
            for X in Xgrid]
        return sample_set

    @classmethod
    def _get_sample_set_from_X_grid(cls, phase_model, Xgrid, T, P, amount):
        sample_set = []
        [sample_set.append(SampleMaker.get_sample(phase_model, X=X, T=T, P=P, amount=amount))
         for X in Xgrid]
        return sample_set

    @classmethod
    def build_comp_grid(cls, Nendmem, grid_spacing, rel_tol=1e-2):
        if Nendmem == 1:
            return np.array([1])

        return SampleMesh.build_sample_mesh(ndim=Nendmem, spacing=grid_spacing)

    @classmethod
    def _append_final_normalized_endmem(cls, Xcoor_independent_norm):
        Xcoor_sum = Xcoor_independent_norm.sum(axis=1)
        Xcoor = np.hstack((Xcoor_independent_norm, 1 - Xcoor_sum[:, np.newaxis]))
        return Xcoor

    @classmethod
    def _filter_results_by_normalization_constraint(cls, Xcoor_independent):
        Xcoor_sum = Xcoor_independent.sum(axis=1)
        Xcoor_independent_norm = Xcoor_independent[Xcoor_sum <= 1]
        return Xcoor_independent_norm

    @classmethod
    def _build_mesh_of_independent_endmems(cls, Nendmem, iXgrid):
        Xgrid_set = [iXgrid for n in range(Nendmem)]
        Xi_values = np.meshgrid(*Xgrid_set[:-1])
        Xcoor_independent = np.vstack([Xi.ravel() for Xi in Xi_values])
        return Xcoor_independent.T

    @classmethod
    def _build_endmem_grid(cls, grid_spacing):
        inv_spacing = 1 / grid_spacing
        division_num = int(np.round(inv_spacing))
        iXgrid = np.linspace(0, 1, division_num + 1)
        return iXgrid

    @classmethod
    def _get_phase_endmember_comps(cls, phs):
        endmem_num = phs.endmember_num
        X_endmems = np.eye(endmem_num) + cls.XTOL
        return X_endmems


class _PurePhaseSample(_PhaseSample):
    def _calc_gibbs_energy(self, X):
        return self.phase_model.gibbs_energy(self.T, self.P)

    def _calc_chem_potential(self, X):
        return self.phase_model.chem_potential(self.T, self.P)


class _SolutionPhaseSample(_PhaseSample):
    def _calc_gibbs_energy(self, X):
        return self.phase_model.gibbs_energy(self.T, self.P, mol=X)

    def _calc_chem_potential(self, X):
        return self.phase_model.chem_potential(self.T, self.P, mol=X)


class SampleMesh:
    """
    Enables mesh refinement for gridding of solution phase
    """
    RTOL = 1e-2

    def __init__(self):
        pass

    @classmethod
    def calc_neighboring_mesh_pts(cls, X0: NDArray,
                                  spacing: float) -> List[NDArray]:
        X_neighbors = cls._calc_neighboring_mesh_coor(X0, spacing)
        return cls._filter_viable_mesh_points(X_neighbors, spacing)

    @classmethod
    def refine_local_mesh(cls, Xlocal0: NDArray,
                          spacing0: float) -> List[NDArray]:
        """Double local mesh resolution by returning mesh midpoints"""
        Xlocal0 = np.array(Xlocal0)
        ind_pairs = cls._find_neighbor_pairs(Xlocal0, spacing0)
        Xnew_mesh_pts = [0.5 * (Xlocal0[i] + Xlocal0[j])
                         for i, j in ind_pairs]
        return Xnew_mesh_pts

    @classmethod
    def refine_mesh_for_multiple_samples(cls, X0_samples: List[NDArray], spacing0: float) \
            -> Tuple[List[NDArray], List[NDArray]]:
        spacing = spacing0 / 2
        decimals_spacing = int(np.ceil(-np.log10(spacing)))

        Xlocal0 = []
        Xgrid_refined = []
        for iX0 in X0_samples:
            iXlocal0 = cls.calc_neighboring_mesh_pts(
                X0=iX0, spacing=spacing0)
            iXgrid_refined = cls.refine_local_mesh(
                iXlocal0 + [iX0], spacing0)

            Xlocal0.extend(iXlocal0)
            Xgrid_refined.extend(iXgrid_refined)

        Xlocal0 = cls._get_roughly_unique_rows(Xlocal0, decimals_spacing + 1)
        Xgrid_refined = cls._get_roughly_unique_rows(Xgrid_refined, decimals_spacing + 2)

        return Xlocal0, Xgrid_refined

    @classmethod
    def _get_roughly_unique_rows(cls, a: NDArray, decimals: int):
        a_round = np.round(a, decimals=decimals)
        a_uniq = np.unique(a_round, axis=0)

        uniq_ind = []
        for a_row in a_uniq:
            row_matches = np.all(a_round == a_row, axis=1)
            ind_match = np.argmax(row_matches)
            uniq_ind.append(ind_match)

        return [a[ind] for ind in uniq_ind]

    @classmethod
    def _find_neighbor_pairs(cls, Xlocal0: List[NDArray],
                             spacing0: float) -> List[NDArray]:
        Npts = len(Xlocal0)
        dX_pairs = (Xlocal0[np.newaxis, :, :] - Xlocal0[:, np.newaxis, :])
        rel_dist_sqr = np.sum((dX_pairs / spacing0) ** 2, axis=2)
        neighbor_pairs = rel_dist_sqr <= 2 * (1 + cls.RTOL)
        ind_pairs = []
        for i in range(Npts):
            for j in range(i + 1, Npts):
                if neighbor_pairs[i, j]:
                    ind_pairs.append([i, j])
        return ind_pairs

    @classmethod
    def _calc_neighboring_mesh_coor(cls, X0: NDArray,
                                    spacing: float) -> List[NDArray]:
        ndim = len(X0)
        iendmems, jendmems = np.triu_indices(ndim, +1)
        X_neighbors = []
        for idim, jdim in zip(iendmems, jendmems):
            dX = np.zeros(ndim)
            dX[idim] = +spacing
            dX[jdim] = -spacing
            X_neighbors.append(X0 + dX)
            X_neighbors.append(X0 - dX)
        return X_neighbors

    @classmethod
    def _filter_viable_mesh_points(cls, X_neighbors: List[NDArray],
                                   spacing: float) -> List[NDArray]:
        XTOL = cls.RTOL * spacing
        X_within_bounds = []
        for Xi in X_neighbors:
            if np.all(Xi > -XTOL) and np.all(Xi < 1 + XTOL):
                Xi[Xi < 0] = 0
                Xi[Xi > 1] = 1
                X_within_bounds.append(Xi)

        return X_within_bounds

    @classmethod
    def get_local_grid_spacing(cls, X0: NDArray, Xsamples: List[NDArray]):
        dX = Xsamples - X0
        dist_sqr = np.sum(dX ** 2, axis=1)
        min_dist_sqr = dist_sqr[dist_sqr > 0].min()
        local_grid_spacing = np.sqrt(min_dist_sqr / 2)
        return local_grid_spacing

    @classmethod
    def build_sample_mesh(cls, ndim: int, spacing: float):
        """
        Create sample mesh with desired grid spacing

        Parameters
        ----------
        ndim: Number of dimensions (components) in the phase
        spacing: sample resolution for mesh

        Returns
        -------
        Xmesh: array of sample compositions (rows) that fully grid the space

        """
        division_num = int(np.round(1 / spacing))

        ind_all_combos = cls._get_all_index_combos(division_num, ndim)
        valid_combos = ind_all_combos.sum(axis=1) == division_num
        X_mesh = ind_all_combos[valid_combos] / division_num
        return X_mesh

    @classmethod
    def _get_all_index_combos(cls, division_num, ndim):
        ind_grid = np.arange(division_num + 1)
        ind_grid_set = [ind_grid for n in range(ndim)]
        ind_combos = np.meshgrid(*ind_grid_set)
        ind_all_combos = np.vstack([inds.ravel() for inds in ind_combos]).T
        return ind_all_combos

    @classmethod
    def mol_comps_equiv(cls, X1, X2, TOL=1e-5):
        if len(X1) != len(X2):
            return False

        mol_tot = np.sum(X2)
        return np.allclose(X1, X2, atol=TOL * mol_tot)


class _FixedCompPhaseSample(_SolutionPhaseSample):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.fixed_comp = True


class SampleLibrary(collections.MutableSequence):
    """
    Stores collection of potential samples to determine relative stability


    Attributes
    ----------
    energies : NDArray[float]
        array of energies for each sample in assemblage

    affinties : NDArray[float]
        array of affinities for each phase

    amounts :  NDArray[float]
        array of molar amounts of every phase
    """

    def __init__(self, samples: Optional[List[_PhaseSample]] = None,
                 safe_copy=True):
        """
        Create new assemblage with desired phase samples

        Parameters
        ----------
        samples :
            list of phase samples for inclusion in assemblage
        """
        if samples is None:
            self._samples = []
        else:
            if safe_copy:
                self._samples = [SampleMaker.copy_sample(samp) for samp in samples]
            else:
                self._samples = [samp for samp in samples]

    def update_conditions(self, T: float = None, P: float = None):
        for samp in self.samples:
            samp.update_conditions(T=T, P=P)

    @property
    def energies(self) -> NDArray[float]:
        return np.array([samp.G for samp in self.samples])

    @property
    def affinities(self) -> NDArray[float]:
        return np.array([samp.aff for samp in self.samples])

    @property
    def amounts(self) -> NDArray[float]:
        return np.array([samp.amount for samp in self.samples])

    @property
    def metastability(self) -> NDArray[bool]:
        return np.array([samp.metastable for samp in self.samples])

    @property
    def names(self) -> NDArray[str]:
        return np.array([samp.name for samp in self.samples])

    @property
    def elem_comps(self) -> pd.DataFrame:
        elem_comp_data = [samp.comp.normalize().all_data for samp in self.samples]
        elem_comps = pd.DataFrame(elem_comp_data, index=self.names)
        drop_mask = (elem_comps == 0).all(axis=0)
        elem_comps = elem_comps.loc[:, ~drop_mask]
        return elem_comps

    @property
    def elems(self) -> List[str]:
        return list(self.elem_comps.columns)

    @property
    def sample_endmem_comps(self) -> NDArray[NDArray[float]]:
        """Return array of molar endmember compositions for every sample"""
        return np.array([samp.X for samp in self.samples])

    @property
    def sample_names(self) -> NDArray[str]:
        return np.array([samp.name for samp in self.samples])

    @property
    def unique_phase_names(self) -> NDArray[str]:
        return np.unique(self.sample_names)

    @property
    def unique_solution_phase_names(self):
        """non-repeating list of solution phase names from assemblage"""
        return [samp.name for samp in self.samples if samp.endmember_num > 1]

    @property
    def multi_sample_phase_names(self) -> NDArray[str]:
        """
        sample phase names for phases that appear multiple times in assemblage
        """
        phase_sample_count = np.array([np.sum(self.sample_names == name)
                                       for name in self.unique_phase_names])
        return self.unique_phase_names[phase_sample_count > 1]

    @property
    def samples(self) -> List[_PhaseSample]:
        return self._samples

    @property
    def pure_samples(self):
        """list of only pure samples from assemblage"""
        return [samp for samp in self.samples if samp.endmember_num == 1]

    @property
    def solution_samples(self):
        """list of only solution samples from assemblage"""
        return [samp for samp in self.samples if samp.endmember_num > 1]

    def update_affinities(self, affinities):
        for samp, aff in zip(self.samples, affinities):
            samp.aff = aff

    def update_amounts(self, amounts):
        for samp, amt in zip(self.samples, amounts):
            samp.amount = amt

    def update_metastability(self, is_metastable: List[bool]):
        for samp, metastable in zip(self.samples, is_metastable):
            samp.metastable = metastable

    def get_subset(self, ind, safe_copy=False):
        return type(self)([self.samples[i] for i in ind], safe_copy=safe_copy)

    def get_nearly_stable_samples(self, aff_thresh: float = 1e2):
        """all nearly stable samples including neighbors in composition space"""
        return type(self)([samp for samp in self.samples
                           if samp.aff < aff_thresh], safe_copy=False)

    def get_nearly_stable_phases(self, aff_thresh, ATOL=1e-3):
        """nearly stable samples with single metastable sample for each phase"""
        metastable_phases = []
        for name in self.unique_phase_names:
            phase_subset = self.get_subset_for_phase(name)
            affinities = phase_subset.affinities
            min_aff = np.min(affinities)
            ind_metastable = np.where(affinities - min_aff <= ATOL)[0]
            metastable_samples = [
                samp for samp in phase_subset.get_subset(ind_metastable)
                if samp.aff <= aff_thresh]

            metastable_phases.extend(metastable_samples)


            # phase_subset.aff
            # icomp = phase_subset.sample_endmem_comps
            # iaffinities = phase_subset.affinities
            #
            # Npts = len(icomp)
            # dX_pairs = (icomp[np.newaxis, :, :] - icomp[:, np.newaxis, :])
            # rel_dist_sqr = np.sum(dX_pairs** 2, axis=2)
            # dX_sqr_neigh = np.sort(rel_dist_sqr, axis=1)[:, 1]
            #
            # imetastable_samp = []
            # for ind in range(Npts):
            #     jdX = dX_sqr_neigh[ind]
            #     jrel_dist_sqr = rel_dist_sqr[ind]
            #     jneigh_mask = jrel_dist_sqr<= jdX*(1.01)
            #
            #     del_aff = iaffinities[jneigh_mask] - iaffinities[ind]
            #
            #
            #     if np.all(del_aff > -1e-3):
            #         print('hello')
            #         # add to metastable point list
            #
            #
            # for isamp in phase_subset:
            #     isamp

            # nearly_stable_subset = phase_subset[phase_subset.affinities < aff_thresh]

            # ind_metastable = np.argmin(phase_subset.affinities)
            # metastable_samp = phase_subset.samples[ind_metastable]
            # if metastable_samp.aff < aff_thresh:
            #     metastable_phases.append(metastable_samp)

        return type(self)(metastable_phases, safe_copy=False)

    def remove_redundant_endmembers(self) -> SampleLibrary:
        """filter out redundant pure samples represented as solution endmembers"""
        endmember_names = self._get_solution_endmember_names()
        filtered_samples = self._get_only_non_endmember_pure_samples(endmember_names)
        filtered_samples += self.solution_samples
        return SampleLibrary(filtered_samples)

    def get_subset_for_phase(self, phase_name: str) -> MonophaseAssemblage:
        """
        Get assemblage of all samples for desired phase

        Parameters
        ----------
        phase_name :
            Desired phase name

        Returns
        -------
        Assemblage of desired phase samples
        """
        mask = [name == phase_name for name in self.sample_names]
        ind = np.where(mask)[0]
        return MonophaseAssemblage([self.samples[i] for i in ind])

    def _get_only_non_endmember_pure_samples(self, endmember_names):
        filtered_samples = []
        for pure_samp in self.pure_samples:
            pure_name = pure_samp.phase_model.phase_name.lower()
            if pure_name not in endmember_names:
                filtered_samples.append(pure_samp)
        return filtered_samples

    def _get_solution_endmember_names(self) -> List[_PhaseSample]:
        endmember_names = []
        for samp in self.solution_samples:
            endmember_names.extend([
                nm.lower() for nm in samp.phase_model.endmember_names])
        return endmember_names

    def __eq__(self, other):
        if other is None:
            other = SampleLibrary()

        return sorted(self.samples) == sorted(other.samples)

    def __getitem__(self, index) -> SampleLibrary:
        ''' retrieves an item by its index'''
        if isinstance(index, slice):
            return type(self)(self.samples[index])

        return self.samples[index]

    def __setitem__(self, key, value):
        ''' set the item at index, key, to value '''
        self._samples[key] = value

    def __delitem__(self, key):
        ''' removes the item at index, key '''
        del self._samples[key]

    def __len__(self):
        return len(self.samples)

    def insert(self, key, value):
        ''' add an item, value, at index, key. '''
        self._samples.insert(key, value)




class Assemblage(SampleLibrary):
    """
    Stores state of coexisting phase samples

    Provides enhanced list for PhaseSample objects that form a
    coexisting phase assemblage. Provides attributes that track
    the state of each phase sample.

    NOTE: each sample is defined on a normalized atomic basis,
      representing 1 mol of atoms
    * The amounts of each phase thus represent the number of total moles
      of atoms present in that phase

    Attributes
    ----------
    samples : List[PhaseSample]
        list of phase samples present in assemblage

    total_energy : float
        total energy of phase assemblage given amounts/comp of each phase
    total_comp : chemistry.Comp
        total composition of assemblage


    sample_amounts : NDArray[float]
        list of molar amounts of each sample on a molar atomic basis
    sample_endmem_comps : NDArray[float]
        list of endmember compositions (X-values) for each sample
    sample_names : NDArray[str]
        list of sample names

    """

    @property
    def total_energy(self) -> float:
        return np.dot(self.sample_amounts, self.energies)

    @property
    def total_comp(self) -> chemistry.Comp:
        sample_elem_comps = self.elem_comps.T
        total_elem_comp = sample_elem_comps.dot(self.sample_amounts)
        return ElemMolComp(**total_elem_comp)

    @property
    def sample_amounts(self) -> NDArray[float]:
        return np.array([samp.amount for samp in self.samples])

    def remove_redundant_endmembers(self):
        return Assemblage(super().remove_redundant_endmembers())

    def get_subset_for_phase(self, phase_name: str) -> MonophaseAssemblage:
        """
        Get assemblage of all samples for desired phase

        Parameters
        ----------
        phase_name :
            Desired phase name

        Returns
        -------
        Assemblage of desired phase samples
        """
        mask = [name == phase_name for name in self.sample_names]
        ind = np.where(mask)[0]
        return MonophaseAssemblage([self.samples[i] for i in ind], safe_copy=False)

    def __eq__(self, other):
        if other is None:
            other = Assemblage()

        return sorted(self.samples) == sorted(other.samples)


class MonophaseAssemblage(Assemblage):
    """Special case of an assemblage with all coexisting samples from same phase"""

    def all_samples_resolved(self, resolution: float, rel_tol: float = 1e-2):
        """returns whether all samples are resolved relative to given resolution for a single phase assemblage"""
        sample_pairs_resolved = self._sample_pairs_resolved(resolution, rel_tol)
        uniq_sample_pairs_resolved = self._get_unique_pairs(sample_pairs_resolved)
        return np.all(uniq_sample_pairs_resolved)

    def segregate_resolved_samples(self, resolution: float,
                                   rel_tol: float = 1e-2) -> List[MonophaseAssemblage]:
        """
        Split monophase assemblage into a list of resolved sample clusters

        Parameters
        ----------
        resolution
            desired molar resolution defining indistinguishable sample clusters
        rel_tol
            relative tolerance on sample resolution

        Returns
        -------
        assem_groups

        """
        sample_pairs_resolved = self._sample_pairs_resolved(resolution, rel_tol)

        if self._all_samples_resolved(sample_pairs_resolved):
            return [MonophaseAssemblage([samp], safe_copy=False) for samp in self.samples]

        assem_pool = MonophaseAssemblage(self.samples, safe_copy=False)

        assem_groups = []
        while (len(assem_pool) > 0):
            sample_cluster = self._extract_sample_cluster_from_pool(assem_pool, resolution, rel_tol)
            assem_groups.append(sample_cluster)

        return assem_groups

    def get_mixed_subset_for_phase(self) -> MonophaseAssemblage:
        """
        Mix  all samples for desired phase and get as single phase assemblage

        Returns
        -------
        single phase Assemblage obtained by mixing selected phase samples
        """
        total_amount = self.sample_amounts.sum()

        sample_frac = self.sample_amounts / total_amount
        X_avg = np.dot(sample_frac, self.sample_endmem_comps)

        samp = self.samples[0]
        mixed_assem = MonophaseAssemblage([samp.modify_sample_conditions(
            X=X_avg, amount=total_amount)])
        return mixed_assem

    def _all_samples_resolved(self, sample_pairs_resolved) -> bool:
        return np.all(self._get_unique_pairs(sample_pairs_resolved))

    def _sample_pairs_resolved(self, resolution: float, rel_tol: float):
        Xsamps = self.sample_endmem_comps
        dX_pairs = np.swapaxes(Xsamps[:, :, np.newaxis] - Xsamps.T, 1, 2)
        axis_resolved = np.abs(dX_pairs) > resolution * (1 + rel_tol)
        sample_pairs_resolved = np.any(axis_resolved, axis=2)
        return sample_pairs_resolved

    def _get_unique_pairs(self, sample_pairs_resolved):
        ind_pairs = np.triu_indices(n=len(self.samples), k=1)
        uniq_sample_pairs_resolved = sample_pairs_resolved[ind_pairs]
        return uniq_sample_pairs_resolved

    def _extract_sample_cluster_from_pool(self, assem_pool: MonophaseAssemblage,
                                          resolution: float, rel_tol: float):
        resolved_pairs_in_pool = assem_pool._sample_pairs_resolved(
            resolution, rel_tol)
        in_cluster = ~resolved_pairs_in_pool[0]
        ind_cluster = np.where(in_cluster)[0]
        cluster_samples = [assem_pool[ind] for ind in ind_cluster]
        [assem_pool.remove(isamp) for isamp in cluster_samples]
        return MonophaseAssemblage(cluster_samples)
