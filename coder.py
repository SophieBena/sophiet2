'''
The Coder module provides code generation for symbolically defined thermodynamics.

This module implements a Python interface that enables the symbolic 
specification of the thermodynamic properties of a phase with code generation 
of the resulting model. Code may be produced for either fast computation of 
results or model parameter calibration. 

Code generation options:
    * Fast computation optimization: fixed parameters
    * Model calibration optimization: variable parameters

Classes are defined to implement standard state properties of pure components, 
thermodynamic properties of simple solutions, and thermodynamic properties of
complex solid solutions that exhibit both convergent and non-convergent cation 
order-disorder.

Thermodynamic properties are implemented for:
    * Pure phases at standard state
    * Simple solutions (asymmetric regular)
    * Complex solutions (convergent & non-convergent cation order-disorder)

By default, models are formulated in terms of the Gibbs free energy, with 
independent variables temperature, pressure, and (for solutions) mole numbers of 
components. Alternatively, models may be implemented in terms of the Helmholtz 
free energy, with independent variables temperature, volume, and mole numbers of 
components.

Model formulation options:
    * Gibbs Free Energy G(T, P, mols) [default]
    * Helmholtz Free Energy A(T, V, mols)


'''
from importlib import import_module

import numpy as np
import sympy as sym
import pandas as pd
import re
import time
import pyximport
import json
from os import path
from itertools import product

from pyximport import install
from sympy.printing.c import C99CodePrinter
from sympy.core.function import Function, UndefinedFunction
from numbers import Number

import thermoengine.chemistry
from thermoengine import core

import thermoengine.coder_templates as tpl

__all__ = ['StdStateModel', 'SimpleSolnModel', 'Debye', 'B', 'Y', 'Q',
'X', 'U', 'N', 'dXdT', 'dUdT', 'dUdP', 'dNdT', 'dNdP', 'Agamma', 'Bgamma',
'AsubG', 'AsubH', 'AsubJ', 'AsubV', 'AsubKappa', 'AsubEx',
'BsubG', 'BsubH', 'BsubJ', 'BsubV', 'BsubKappa', 'BsubEx', 'Expression', 'f_mu_e', 'gSolvent',
'Implicit_Function', 'ComplexSolnModel', 'Ordering_Functions', 'Parameter', 
'SpeciationSolnModel', 'SubCodePrinter']

DATADIR = 'data/coder'
AQUEOUSDIR = 'aqueous'


#####################################################
# Class to extend SymPy code printer                #
# Deals with powers and user function substitutions #
#####################################################

class SubCodePrinter(C99CodePrinter):
    """
    Subclass of the C99 code printer class in SymPy

    Extension expands :math:`x^2`, :math:`x^3`, and :math:`x^4` as multiplication rather than 
    pow(x, (double)n).
    
    The second implements a version of user function derivative printing 
    specified to endmember chemical potential derivatives called externally 
    from a C structure. (See implementation below.)
    
    Attributes
    ----------
    forIndex
    nBasis
    sum_result
    
    """
    
    def __init__(self, settings=None, nBasis=0, forIndex='i', 
        sum_result='result'):
        super().__init__(settings=settings)
        self._nBasis = nBasis
        self._forIndex = forIndex
        self._sum_result = sum_result
        self._print_subs = None
    
    @property
    def nBasis(self):
        """
        Number of basis species used in printing derivatives of standard state 
        chemical potentials
        
        See _print_derivative method for the function mu_s to illustrate 
        functionality. nBasis is the number of basis species in an instance
        of the SpeciationSolnModel class.
        
        Returns:
            str
        """
        return self._nBasis

    @nBasis.setter
    def nBasis(self, nBasis):
        self._nBasis = nBasis
    
    @property
    def forIndex(self):
        """
        Index variable used in for loop generation for printing derivatives of
        standard state chemical potentials
        
        See _print_derivative method for the function mu_s to illustrate 
        functionality. forIndex is the string equivalent of the index variable 
        (SymPy class Idx) used in summation expressions that are fed to the
        SpeciationSolnModel class.  
        
        Returns:
            [type] -- [description]
        """
        return self._forIndex

    @forIndex.setter
    def forIndex(self, forIndex):
        self._forIndex = forIndex

    @property
    def sum_result(self):
        """
        Name of assignment variable used when printing summations 
        
        See _print_Sum method to illustrate functionality. When an instance of
        the SymPy class Sum is printed, the result is assigned to a variable of
        this name.
        
        Returns:
            str
        """
        return self._sum_result

    @sum_result.setter
    def sum_result(self, sum_result):
        self._sum_result = sum_result

    @property
    def print_subs(self):
        """
        A list of tuples of expression substitutions used in doprint()
        
        Each tuple is of the form (sympy expression to be replaced, sympy 
        expression to be substituted). Replacement is done in expressions 
        passed to the doprint() method.
        
        Returns:
            list of tuples
        """
        return self._print_subs

    @print_subs.setter
    def print_subs(self, print_subs):
        self._print_subs = print_subs

    def doprint(self, expr, assign_to=None):
        if self.print_subs is not None:
            try:
                expr = expr.subs(self.print_subs)
            except AttributeError:
                pass
        return super(SubCodePrinter,self).doprint(expr, assign_to=assign_to)

    def _print_Pow(self, expr):
        if expr.exp.is_integer and expr.exp > 0 and expr.exp <= 4:
            result = ')*('.join([self._print(expr.base) for i in range(expr.exp)])
            return '((' + result + '))'
        else:
            return super()._print_Pow(expr)

    def _print_Sum(self, expr):
        ind = str(expr.limits[0][0])
        low = str(expr.limits[0][1])
        high = str(expr.limits[0][2])
        func = expr.function
        result =  '{\n'
        result += 'double sum = 0.0;\n'
        result += 'for (int '+ind+'='+low+'; '+ind+'<='+high+'; '+ind+'++) {\n'
        result += 'sum += ' + self.doprint(func) + ';\n'
        result += '}\n'
        result += self._sum_result + ' += sum;\n'
        result += '}\n'
        return result

    def _print_Derivative(self, expr):
        function, *vars = expr.args
        number_of_derivatives = len(expr.args) - 1
        
        if function.func.__name__[0:4] == 'mu_s':
            if number_of_derivatives == 1:
                derivative_string = repr(vars[0][0])
                derivative_order  = '' if vars[0][1] == 1 else str(vars[0][1])
                result = ('(*endmember['+str(self.nBasis)+'+'+self.forIndex+'-1].d' 
                          + derivative_order + 'mu0d' 
                          + derivative_string + derivative_order + ')(T, P)')
            elif number_of_derivatives == 2:
                derivative_string_2 = repr(vars[0][0])
                derivative_order_2  = '' if vars[0][1] == 1 else str(vars[0][1])
                derivative_string_1 = repr(vars[1][0])
                derivative_order_1  = '' if vars[1][1] == 1 else str(vars[1][1])
                derivative_total    = str(vars[0][1]+vars[1][1])
                result = ('(*endmember['+str(self.nBasis)+'+'+self.forIndex+'-1].d' 
                          + derivative_total + 'mu0d' 
                          + derivative_string_1 + derivative_order_1 +'d' 
                          + derivative_string_2 + derivative_order_2 + ')(T, P)')
            else:
                result = ''

        elif function.func.__name__[0:2] == 'mu':
            function_string_index = (
                int(sym.srepr(function).split("'")[1][2:]) - 1)
            if number_of_derivatives == 1:
                derivative_string = repr(vars[0][0])
                derivative_order  = '' if vars[0][1] == 1 else str(vars[0][1])
                result = ('(*endmember[' + str(function_string_index) + '].d' 
                          + derivative_order + 'mu0d' 
                          + derivative_string + derivative_order + ')(T, P)')
            elif number_of_derivatives == 2:
                derivative_string_2 = repr(vars[0][0])
                derivative_order_2  = '' if vars[0][1] == 1 else str(vars[0][1])
                derivative_string_1 = repr(vars[1][0])
                derivative_order_1  = '' if vars[1][1] == 1 else str(vars[1][1])
                derivative_total    = str(vars[0][1]+vars[1][1])
                result = ('(*endmember[' + str(function_string_index) + '].d' 
                          + derivative_total + 'mu0d' 
                          + derivative_string_1 + derivative_order_1 +'d' 
                          + derivative_string_2 + derivative_order_2 + ')(T, P)')
            else:
                result = ''
        
        elif (len(function.func.__name__) >= 6 and 
            function.func.__name__[1:6] == 'gamma'):
            if number_of_derivatives == 1:
                derivative_string = repr(vars[0][0]).lower()
                derivative_order  = '' if vars[0][1] == 1 else str(vars[0][1])
                result = ('d' + derivative_order + function.func.__name__ + 'd' 
                    + derivative_string + derivative_order + '(T, P)')
            elif number_of_derivatives == 2:
                derivative_string_2 = repr(vars[0][0]).lower()
                derivative_order_2  = '' if vars[0][1] == 1 else str(vars[0][1])
                derivative_string_1 = repr(vars[1][0]).lower()
                derivative_order_1  = '' if vars[1][1] == 1 else str(vars[1][1])
                derivative_total    = str(vars[0][1]+vars[1][1])
                result = ('d' + derivative_total + function.func.__name__ + 'D' 
                    + derivative_string_1 + derivative_order_1 +'D' 
                    + derivative_string_2 + derivative_order_2 + '(T, P)')
            else:
                result = ''

        elif (len(function.func.__name__) >= 8 and 
            function.func.__name__[0:9] == 'gSolvent'):
            if number_of_derivatives == 1:
                derivative_string = repr(vars[0][0]).lower()
                derivative_order  = '' if vars[0][1] == 1 else str(vars[0][1])
                result = ('D' + derivative_order + function.func.__name__ + 'D' 
                    + derivative_string + derivative_order + '(T, P)')
            elif number_of_derivatives == 2:
                derivative_string_2 = repr(vars[0][0]).lower()
                derivative_order_2  = '' if vars[0][1] == 1 else str(vars[0][1])
                derivative_string_1 = repr(vars[1][0]).lower()
                derivative_order_1  = '' if vars[1][1] == 1 else str(vars[1][1])
                derivative_total    = str(vars[0][1]+vars[1][1])
                result = ('D' + derivative_total + function.func.__name__ + 'D' 
                    + derivative_string_1 + derivative_order_1 +'D' 
                    + derivative_string_2 + derivative_order_2 + '(T, P)')
            else:
                result = ''

        else:
            if (not isinstance(type(function), UndefinedFunction) or 
                not all(isinstance(i, Symbol) for i in vars)):
                return super()._print_Derivative(expr)
        return result

#####################################################
# SymPy sub-classes implementing the Born Functions #
#####################################################

class dUdT(sym.Function):
    """
    SymPy Function package extension: Third partial derivative of the Born 
    function with respect to pressure and twice with respect to temperature

    Note that :math:`U = \\frac{{{\\partial ^2}B}}{{\\partial T\\partial P}}`.
    """
    nargs = (1,2)

class dUdP(sym.Function):
    """
    SymPy Function package extension: Third partial derivative of the Born 
    function with respect to temperature and twice with respect to pressure

    Note that :math:`U = \\frac{{{\\partial ^2}B}}{{\\partial T\\partial P}}`.  
    Identical to *dNdT*.
    """
    nargs = (1,2)

class dNdT(sym.Function):
    """
    SymPy Function package extension: Third partial derivative of the Born 
    function with respect to temperature and twice with respect to pressure

    Note that :math:`N = \\frac{{{\\partial ^2}B}}{{\\partial {P^2}}}`.
    """
    nargs = (1,2)

class dNdP(sym.Function):
    """
    SymPy Function package extension: Third partial derivative of the Born 
    function with respect to pressure

    Note that :math:`N = \\frac{{{\\partial ^2}B}}{{\\partial {P^2}}}`.
    """
    nargs = (1,2)

class dXdT(sym.Function):
    """
    SymPy Function package extension: Third partial derivative of the Born 
    function with respect to temperature

    Note that :math:`X = \\frac{{{\\partial ^2}B}}{{\\partial {T^2}}}`.   
    """
    nargs = (1,2)

    def fdiff(self, argindex=1):
        """
        Provides a definition for the derivative of the function.
        """
        return sym.S.Zero

class X(sym.Function):
    """
    SymPy Function package extension: Second partial derivative of the Born 
    function with respect to temperature

    Definition: :math:`X = \\frac{1}{\\varepsilon }\\left[ {{{\\left( {\\frac{{{\\partial ^2}\\ln \\varepsilon }}{{\\partial {T^2}}}} \\right)}_P} - \\left( {\\frac{{\\partial \\ln \\varepsilon }}{{\\partial T}}} \\right)_P^2} \\right]` 
    
    Note that :math:`X = \\frac{{{\\partial ^2}B}}{{\\partial {T^2}}}`.
    """
    nargs = (1, 2)

    def fdiff(self, argindex=1):
        """
        Provides a definition for the derivative of the function.
        """
        T,P = self.args
        if argindex == 1:
            return dXdT(T,P)
        elif argindex == 2:
            return dUdT(T,P)
        raise ArgumentIndexError(self, argindex)

class N(sym.Function):
    """
    SymPy Function package extension: Second partial derivative of the Born 
    function with respect to pressure

    Definition: :math:`N = \\frac{1}{\\varepsilon }\\left[ {{{\\left( {\\frac{{{\\partial ^2}\\ln \\varepsilon }}{{\\partial {P^2}}}} \\right)}_T} - \\left( {\\frac{{\\partial \\ln \\varepsilon }}{{\\partial P}}} \\right)_T^2} \\right]` 
    
    Note that :math:`N = \\frac{{{\\partial ^2}B}}{{\\partial {P^2}}}`.
    """
    nargs = (1, 2)

    def fdiff(self, argindex=1):
        """
        Provides a definition for the derivative of the function.
        """
        T,P = self.args
        if argindex == 1:
            return dNdT(T,P)
        elif argindex == 2:
            return dNdP(T,P)
        raise ArgumentIndexError(self, argindex)

class U(sym.Function):
    """
    SymPy Function package extension: Second partial derivative of the Born 
    function with respect to temperature and pressure

    Definition: :math:`U = \\frac{1}{\\varepsilon }\\left[ {\\left( {\\frac{{{\\partial ^2}\\ln \\varepsilon }}{{\\partial T\\partial P}}} \\right) - {{\\left( {\\frac{{\\partial \\ln \\varepsilon }}{{\\partial T}}} \\right)}_P}{{\\left( {\\frac{{\\partial \\ln \\varepsilon }}{{\\partial P}}} \\right)}_T}} \\right]` 
    
    Note that :math:`U = \\frac{{{\\partial ^2}B}}{{\\partial T\\partial P}}`.
    """
    nargs = (1, 2)

    def fdiff(self, argindex=1):
        """
        Provides a definition for the derivative of the function.
        """    
        T,P = self.args
        if argindex == 1:
            return dUdT(T,P)
        elif argindex == 2:
            return dUdP(T,P)
        raise ArgumentIndexError(self, argindex)

class Y(sym.Function):
    """
    SymPy Function package extension: Temperature derivative of the Born function

    Definition: :math:`Y = \\frac{1}{\\varepsilon }{\\left( {\\frac{{\\partial \\ln \\varepsilon }}{{\\partial T}}} \\right)_P}` 
    
    Note that :math:`Y = \\frac{{\\partial B}}{{\\partial T}}`.
    """
    nargs = (1, 2)
    
    def fdiff(self, argindex=1):
        """
        Provides a definition for the derivative of the function.
        """
        T,P = self.args
        if argindex == 1:
            return X(T,P)
        elif argindex == 2:
            return U(T,P)
        raise ArgumentIndexError(self, argindex)

class Q(sym.Function):
    """
    SymPy Function package extension: Pressure derivative of the Born function

    Definition: :math:`Q = \\frac{1}{\\varepsilon }{\\left( {\\frac{{\\partial \\ln \\varepsilon }}{{\\partial P}}} \\right)_T}` 
    
    Note that :math:`Q = \\frac{{\\partial B}}{{\\partial P}}`.
    """
    nargs = (1, 2)
    
    def fdiff(self, argindex=1):
        """
        Provides a definition for the derivative of the function.
        """
        T,P = self.args
        if argindex == 1:
            return U(T,P)
        elif argindex == 2:
            return N(T,P)
        raise ArgumentIndexError(self, argindex)
    
class B(sym.Function):
    """
    SymPy Function package extension: Born function

    Definition: :math:`B = - \\frac{1}{\\varepsilon }`, where 
    :math:`\\varepsilon` is the dielectric constant of water
    """
    nargs = (1, 2)
    
    def fdiff(self, argindex=1):
        """
        Provides a definition for the derivative of the function.
        """
        T,P = self.args
        if argindex == 1:
            return Y(T,P)
        elif argindex == 2:
            return Q(T,P)
        raise ArgumentIndexError(self, argindex)

#############################################################
# SymPy sub-classes implementing the Debye-Hückel Functions #
#############################################################

class AsubEx(sym.Function):
    """
    SymPy Function package extension: Temperature derivative of AsubV 

    Definition: :math:`{A_{Ex} } = {\\left( {\\frac{{\\partial {A_V}}}{{\\partial T}}} \\right)_P}`
    """
    nargs = (1, 2)

class AsubKappa(sym.Function):
    """
    SymPy Function package extension: Pressure derivative of AsubV 

    Definition: :math:`{A_\\kappa } = {\\left( {\\frac{{\\partial {A_V}}}{{\\partial P}}} \\right)_T}`
    """
    nargs = (1, 2)

class AsubJ(sym.Function):
    """
    SymPy Function package extension: Temperature derivative of AsubH 

    Definition: :math:`{A_J} = {\\left( {\\frac{{\\partial {A_H}}}{{\\partial T}}} \\right)_P}`
    """
    nargs = (1, 2)

class AsubV(sym.Function):
    """
    SymPy Function package extension: Pressure derivative of AsubG 

    Definition: :math:`{A_V} = {\\left( {\\frac{{\\partial {A_G}}}{{\\partial P}}} \\right)_T}`
    """
    nargs = (1, 2)

    def fdiff(self, argindex=1):
        """
        Provides a definition for the derivative of the function.
        """
        T,P = self.args
        if argindex == 1:
            return AsubEx(T,P)
        elif argindex == 2:
            return AsubKappa(T,P)
        raise ArgumentIndexError(self, argindex)

class AsubH(sym.Function):
    """
    SymPy Function package extension: Enthalpy equivalent of the Debye-Hückel 
    function Agamma 

    Definition: :math:`{A_H} = {A_G} - T{\\left( {\\frac{{\\partial {A_G}}}{{\\partial T}}} \\right)_P}`
    
    Pressure derivative: :math:`\\frac{{\\partial {A_H}}}{{\\partial P}} = \\frac{{\\partial {A_G}}}{{\\partial P}} - T\\frac{{{\\partial ^2}{A_G}}}{{\\partial T\\partial P}}` 
    
    or, :math:`\\frac{{\\partial {A_H}}}{{\\partial P}} = {A_V} - T{A_{Ex}}`
    """
    nargs = (1, 2)

    def fdiff(self, argindex=1):
        """
        Provides a definition for the derivative of the function.
        """
        T,P = self.args
        if argindex == 1:
            return AsubJ(T,P)
        elif argindex == 2:
            return AsubV(T,P) - T*AsubEx(T,P)
        raise ArgumentIndexError(self, argindex)

class AsubG(sym.Function):
    """
    SymPy Function package extension: Gibbs free energy equivalent of the
    Debye-Hückel function Agamma 

    Definition: :math:`{A_G} =  - 2\\left( {\\ln 2} \\right)RT{A_\\gamma }`

    """
    nargs = (1, 2)

    def fdiff(self, argindex=1):
        """
        Provides a definition for the derivative of the function.
        """
        T,P = self.args
        if argindex == 1:
            return -(AsubG(T,P) - AsubH(T,P))/T
        elif argindex == 2:
            return AsubV(T,P)
        raise ArgumentIndexError(self, argindex)

class Agamma(sym.Function):
    """
    SymPy Function package extension: Debye-Hückel function Agamma 
    
    Temperature and pressure derivatives to third order are coded for this
    function.
    """
    nargs = (1, 2)

class BsubEx(sym.Function):
    """
    SymPy Function package extension: Temperature derivative of BsubV 

    Definition: :math:`{B_{Ex} } = {\\left( {\\frac{{\\partial {B_V}}}{{\\partial T}}} \\right)_P}`
    """
    nargs = (1, 2)

class BsubKappa(sym.Function):
    """
    SymPy Function package extension: Pressure derivative of BsubV 

    Definition: :math:`{B_\\kappa } = {\\left( {\\frac{{\\partial {B_V}}}{{\\partial P}}} \\right)_T}`
    """
    nargs = (1, 2)

class BsubJ(sym.Function):
    """
    SymPy Function package extension: Temperature derivative of BsubH 

    Definition: :math:`{B_J} = {\\left( {\\frac{{\\partial {B_H}}}{{\\partial T}}} \\right)_P}`
    """
    nargs = (1, 2)

class BsubV(sym.Function):
    """
    SymPy Function package extension: Pressure derivative of BsubG 

    Definition: :math:`{B_V} = {\\left( {\\frac{{\\partial {B_G}}}{{\\partial P}}} \\right)_T}`
    """
    nargs = (1, 2)

    def fdiff(self, argindex=1):
        """
        Provides a definition for the derivative of the function.
        """
        T,P = self.args
        if argindex == 1:
            return BsubEx(T,P)
        elif argindex == 2:
            return BsubKappa(T,P)
        raise ArgumentIndexError(self, argindex)

class BsubH(sym.Function):
    """
    SymPy Function package extension: Enthalpy equivalent of the Debye-Hückel 
    function Bgamma 

    Definition: :math:`{B_H} = {B_G} - T{\\left( {\\frac{{\\partial {B_G}}}{{\\partial T}}} \\right)_P}`
    
    Pressure derivative: :math:`\\frac{{\\partial {B_H}}}{{\\partial P}} = \\frac{{\\partial {B_G}}}{{\\partial P}} - T\\frac{{{\\partial ^2}{B_G}}}{{\\partial T\\partial P}}` 
    
    or, :math:`\\frac{{\\partial {B_H}}}{{\\partial P}} = {B_V} - T{B_{Ex}}`
    """
    nargs = (1, 2)

    def fdiff(self, argindex=1):
        """
        Provides a definition for the derivative of the function.
        """
        T,P = self.args
        if argindex == 1:
            return BsubJ(T,P)
        elif argindex == 2:
            return BsubV(T,P) - T*BsubEx(T,P)
        raise ArgumentIndexError(self, argindex)

class BsubG(sym.Function):
    """
    SymPy Function package extension: Gibbs free energy equivalent of the
    Debye-Hückel function Bgamma 

    Definition: :math:`{B_G} =  - 2\\left( {\\ln 2} \\right)RT{B_\\gamma }`
    """
    nargs = (1, 2)

    def fdiff(self, argindex=1):
        """
        Provides a definition for the derivative of the function.
        """
        T,P = self.args
        if argindex == 1:
            return -(BsubG(T,P) - BsubH(T,P))/T
        elif argindex == 2:
            return BsubV(T,P)
        raise ArgumentIndexError(self, argindex)

class Bgamma(sym.Function):
    """
    SymPy Function package extension: Debye-Hückel function Bgamma 
    
    Temperature and pressure derivatives to third order are coded for this
    function.
    """
    nargs = (1, 2)

#########################################################
# SymPy sub-classes implementing the g solvent function #
#########################################################

class gSolvent(sym.Function):
    """
    SymPy Function package extension: g solvent function 
    
    Temperature and pressure derivatives to second order are coded for 
    this function.
    """
    nargs = (1, 2)

#####################################################
# SymPy sub-classes implementing the Debye integral #
#####################################################

class Debye(Function):
    """
    SymPy Function package extension: Debye function integral

    Defines symbolic derivatives using a recurrence relation.
    """

    @classmethod
    def eval(cls, x, n=None):
        if x is sym.S.Zero:
            return sym.S.One

    def fdiff(self, argindex=1):
        """
        Provides a definition for the derivative of the function.
        """
        x = self.args[0]

        if len(self.args) == 1:
            # if called without a value of n,  assume n=3
            n = 3
            if argindex == 1:
                return n/(sym.exp(x)-1) - n*Debye(x)/x
        else:
            n = self.args[1]
            if argindex == 1:
                return n/(sym.exp(x)-1) - n*Debye(x,n)/x
       
        raise ArgumentIndexError(self, argindex)

####################################################################
# Classes to store and organize information for the public classes #
####################################################################

class Parameter:
    """
    Class holds properties of a model parameter.

    Parameters
    ----------
    name : str
        Name of the parameter
    units : str
        Units of the parameter
    expression : sympy.core.symbol.Symbol
        Instance of a SymPy symbol class equivalent to name

    Attributes
    ----------
    expression 
    name 
    units

    Notes
    -----
    Class is utilized principally to store model parameters and model variables
    used by the Expression class to construct expressions for the Gibbs free
    energy or Helmholtz energy.

    """
    def __init__(self, name, units, expression):
        assert isinstance(name, str)
        assert isinstance(units, str)
        assert isinstance(expression, sym.Symbol)
        self._name = name
        self._units = units
        self._expression = expression 

    @property
    def name(self):
        """
        Name of the parameter
        
        Returns
        -------
        String      
        """
        return self._name

    @property
    def units(self):
        """
        Units of the parameter
        
        Returns
        -------
        String      
        """
        return self._units
    
    @property
    def expression(self):
        """
        Instance of a SymPy symbol class equivalent to property name
        
        Returns
        -------
        sympy.core.symbol.Symbol        
        """
        return self._expression

class Expression:
    """
    Class holds properties of a model expression.

    Parameters
    ----------
    expression : sympy.core.symbol.Symbol
        Instance of a SymPy symbol class that contains the model expression
    params : array of instances of the Parameter class
        Parameter class instances provided are relevant to the expression
    exp_type : str
        'unrestricted' if the expression applies to the whole of (T,P) or (T,V)
        space
        
        'restricted' if the expression applies to a portion of (T,P) or (T,V)
        space, specified by lower_limits and upper_limits
    lower_limit_T : Number
        Lower limit of applicability of expression in T units
    lower_limit_PorV : Number
        Lower limit of applicability of expression in P or V units
    upper_limit_T : Number
        Upper limit of applicability of expression in T units
    upper_limit_PorV : Number
        Upper limit of applicability of expression in P or V units

    Attributes
    ----------
    exp_type  
    expression  
    lower_limit_PorV  
    lower_limit_T  
    params  
    upper_limit_PorV  
    upper_limit_T  

    Notes
    -----
    A generic expression class for defining a thermodynamic model. Multiple
    expressions may be used to construct the model.
    """
    def __init__(self, expression, params, exp_type='unrestricted', 
        lower_limit_T=None, lower_limit_PorV=None, 
        upper_limit_T=None, upper_limit_PorV=None):
        assert exp_type in ['unrestricted', 'restricted']
        assert (isinstance(x, Parameter) for x in params)

        self._exp_type = exp_type
        self._expression = expression
        self._params = []
        for param in params:
            self._params.append(param)
        self._lower_limit_T = lower_limit_T
        self._lower_limit_PorV = lower_limit_PorV
        self._upper_limit_T = upper_limit_T
        self._upper_limit_PorV = upper_limit_PorV

    @property
    def exp_type(self):
        """
        Type of expression
        
        - 'unrestricted' if the expression applies to the whole of (T,P) or (T,V)
          space         
        - 'restricted' if the expression applies to a portion of (T,P) or (T,V)
          space, specified by lower_limits and upper_limits
        
        Returns
        -------
        String      
        """
        return self._exp_type
    
    @property
    def expression(self):
        """
        Instance of a SymPy symbol class that contains the model expression
        
        Returns
        -------
        sympy.core.symbol.Symbol
        """
        return self._expression
    
    @property
    def params(self):
        """
        Array of instances of the Parameter class
        
        Returns
        -------
        Parameter class instances relevant to the expression
        """
        return self._params
    
    @property
    def lower_limit_T(self):
        """
        Lower limit of applicability of expression in temperature (K)
        
        Returns
        -------
        Float
        """
        return self._lower_limit_T
    
    @property
    def lower_limit_PorV(self):
        """
        Lower limit of applicability of expression in pressure (bars) or volume (J/bar)
        
        Returns
        -------
        Float
        """
        return self._lower_limit_PorV
    
    @property
    def upper_limit_T(self):
        """
        Upper limit of applicability of expression in temperature (K)
        
        Returns
        -------
        Float   
        """
        return self._upper_limit_T

    @property
    def upper_limit_PorV(self):
        """
        Upper limit of applicability of expression in pressure (bars) or volume (J/bar)
            
        Returns
        -------
        Float   
        """
        return self._upper_limit_PorV

class Implicit_Function:
    """
    Class holds properties of an Implicit Function

    Parameters
    ----------
    implicit_function : sympy.core.symbol.Symbol
        A SymPy expression for an implicit function in two independent (T,P or 
        T,V) variables and one dependent variable (*f*). The function is equal 
        to  zero; e.g., :math:`log(f) - f = 0`. This expression must be defined in 
        terms of known parameters and Tr, Pr, T, P.
    dep_variable : sympy.core.symbol.Symbol
        A SymPy symbol for the dependent variable *f*
    initial_guess : sympy.core.symbol.Symbol
        A SymPy expression that initializes *f* in the iterative routine. This 
        expression must be defined in terms of known parameters and Tr, Pr, 
        T, P.

    Attributes
    ----------
    function
    guess   
    variable 

    Notes
    -----
    Implicit functions are solved prior to evaluating model expressions for the
    Gibbs free energy or Helmholtz energy.  They must be differentiable with 
    respect to model parameters T and P, or T and V.  Typically, implicit 
    functions are used to support expressions for equations of state, e.g. 
    the Birch Murnaghan expression, which is implicit in V and T, when coupled 
    with a reference pressure expression for the Gibbs free energy, which is 
    explicit in T and P, requires an implicit function that solves for V given
    an input value of P.  
    """
    def __init__(self, implicit_function, dep_variable, initial_guess):
        assert isinstance(dep_variable, sym.Function)
        self._function = implicit_function
        self._variable = dep_variable
        self._guess = initial_guess

    @property
    def function(self):
        """
        A SymPy expression for an implicit function in two independent (T,P or 
        T,V) variables and one dependent variable (*f*). The function is equal 
        to zero; e.g., :math:`log(f) - f = 0`. This expression must be defined in 
        terms of known parameters and Tr, Pr, T, P.
        
        Returns
        -------
        sympy.core.symbol.Symbol
        """
        return self._function
    
    @property
    def variable(self):
        """
        A SymPy symbol for the dependent variable *f*.
    
        Returns
        -------
        sympy.core.symbol.Symbol
        """
        return self._variable
    
    @property
    def guess(self):
        """
        A SymPy expression that initializes *f* in the iterative routine. This 
        expression must be defined in terms of known parameters and Tr, Pr, 
        T, P.
        
        Returns
        -------
        sympy.core.symbol.Symbol
        """
        return self._guess

##################
# Public classes #
##################

class StdStateModel:
    """
    Class creates representation of standard state properties of pure phase.

    Parameters
    ----------
    model_type : str
        Model type of 'TP' implies that expressions are of the form of the Gibbs 
        free energy.
        
        Model type of 'TV' implies that expressions are of the form of the 
        Helmholtz energy. This option is infrequently used.

    Attributes
    ----------
    a_list  
    expression_parts  
    g_list  
    implicit_functions  
    module  
    params  
    printer  
    variables
    
    Methods
    -------
    add_expression_to_model
    create_calc_h_file 
    create_calib_h_file
    create_code_module 
    get_berman_std_state_database
    get_include_born_code
    get_include_debye_code
    get_model_param_names
    get_model_param_symbols
    get_model_param_units
    get_module_name
    get_symbol_for_p
    get_symbol_for_pr
    get_symbol_for_t
    get_symbol_for_tr
    get_symbol_for_v
    get_symbol_for_vr
    parse_formula
    set_include_born_code
    set_include_debye_code
    set_module_name
    set_reference_origin

    Notes
    -----
    Class for creating a representation of the standard state properties of a 
    stoichiometric phase.  The principal use of the class is to construct a 
    model symbolically and to generate code that implements the model for model 
    parameter calibration and thermodynamic property calculation.

    The class supports standard state models that involve integrals of the 
    Debye function, implicit functions for equations of state, and models that
    require Born functions for estimation of the dielectric properties of 
    water.

    """
    def __init__(self, model_type='TP'):
        function_d = {"Debye":"Debye", "B":"born_B", "Y":"born_Y", "Q":"born_Q",
        "X":"born_X", "U":"born_U", "N":"born_N", "dXdT":"born_dXdT",
        "dUdT":"born_dUdT", "dNdT":"born_dNdT", "dUdP":"born_dUdP",
        "dNdP":"born_dNdP", "gSolvent":"gSolvent"}
        self._printer = SubCodePrinter(settings={"user_functions":function_d})
        self._g_list = [('g',0,0), ('dgdt',1,0), ('dgdp',0,1), ('d2gdt2',2,0), 
                        ('d2gdtdp',1,1), ('d2gdp2',0,2), ('d3gdt3',3,0), 
                        ('d3gdt2dp',2,1), ('d3gdtdp2',1,2), ('d3gdp3',0,3)]
        self._a_list = [('a',0,0), ('dadt',1,0), ('dadv',0,1), ('d2adt2',2,0), 
                        ('d2adtdv',1,1), ('d2adv2',0,2), ('d3adt3',3,0), 
                        ('d3adt2dv',2,1), ('d3adtdv2',1,2), ('d3adv3',0,3)]
        self._module = 'untitled'
        self._params = []
        self._variables = []
        self._expression_parts = []
        self._implicit_functions = []
        self._model_type = model_type
        self._T_r = 298.15 # Kelvins
        self._P_r = 1.0    # bars
        self._V_r = 2.5    # J/bar
        assert (model_type in ['TP', 'TV'])
        if model_type == 'TP':
            self._params.append(Parameter('T_r', 'K', sym.symbols('T_r')))
            self._params.append(Parameter('P_r', 'bar', sym.symbols('P_r')))
            self._variables.append(Parameter('T', 'K', sym.symbols('T')))
            self._variables.append(Parameter('P', 'bar', sym.symbols('P')))
        elif model_type == 'TV':
            self._params.append(Parameter('T_r', 'K', sym.symbols('T_r')))
            self._params.append(Parameter('V_r', 'J/bar', sym.symbols('V_r')))
            self._variables.append(Parameter('T', 'K', sym.symbols('T')))
            self._variables.append(Parameter('V', 'J/bar', sym.symbols('V')))
        else:
            print ('ERROR: Unsupported model_type. Must be either "TP" or "TV"')
        self._elements = pd.read_csv(path.join(path.dirname(__file__), DATADIR, 
            'elements.csv'))
        self._entropies = pd.read_csv(path.join(path.dirname(__file__), DATADIR, 
            'entropies.csv'))
        self._berman_db = None
        self._include_debye_code = False
        self._include_born_code = False
        self.calc_h = None
        self.calib_h = None
        self.fast_h_template = None
        self.calib_h_template = None
        self.fast_c_template = None 
        self.calib_c_template = None
        return

    @property
    def a_list(self):
        """
        Array of tuples used internally by the class
        
        Each tuple contains the following entries:   
    
        - [0] str - Function name for derivatives of the Helmholtz energy         
        - [1] int - Order of temperature derivative      
        - [2] int - Order of volume derivative
          
        Returns
        -------
        Array of tuples
        """
        return self._a_list
    
    @a_list.setter
    def a_list(self, a_list):
        self._a_list = a_list

    @property
    def g_list(self):
        """
        Array of tuples used internally by the class
        
        Each tuple contains the following entries:

        - [0] str - Function name for derivatives of the Gibbs free energy      
        - [1] int - Order of temperature derivative     
        - [2] int - Order of pressure derivative
        
        Generally for internal use by the class

        Returns
        -------
        Array of tuples
        """     
        return self._g_list
    
    @g_list.setter
    def g_list(self, g_list):
        self._g_list = g_list

    @property
    def expression_parts(self):
        """     
        Array of Expression class instances  

        A model is built by layering expression instances. Each instance may
        be applicable over all or a portion of the domain space of variables.       

        Returns
        -------
        Array of Expression class instances, [str,...]
        """ 
        return self._expression_parts
    
    @expression_parts.setter
    def expression_parts(self, expression_parts):
        self._expression_parts = expression_parts

    @property
    def implicit_functions(self):
        """     
        Array of Implicit_Function class instances

        A model may require that one or more implicit function expressions be 
        satisfied prior to evaluation of the model expressions. 
   
        Returns
        -------
        Array of Implicit_Function class instances, [str,...]
        """ 
        return self._implicit_functions
    
    @implicit_functions.setter
    def implicit_functions(self, implicit_functions):
        self._implicit_functions = implicit_functions

    @property
    def module(self):
        """     
        Name of module
                
        Returns
        -------
        Name of module (str)
        """
        return self._module

    @module.setter
    def module(self, module):
        assert isinstance(module, str)
        self._module = module

    @property
    def params(self):
        """     
        Array of Parameter class instances  

        Parameters are calibrated quantities that define the model, like T_r 
        and P_r.
        
        Returns
        -------
        Array of Parameter class instances, [Parameter,...]
        """
        return self._params
    
    @params.setter
    def params(self, params):
        self._params = params

    @property
    def printer(self):
        """
        Instance of a SymPy code printer class used in code generation
        
        Returns
        -------
        C99CodePrinter
        """
        return self._printer
    
    @printer.setter
    def printer(self, printer):
        self._printer = printer

    @property
    def variables(self):
        """
        Array of Parameter class instances  
        
        This property returns or sets an array of variables, either T, P, T_r and P_r 
        (if model_type is ’TP’) or T, V, T_r and V_r (if model type is equal to ’TV’). 
        Array elements are encapsulated as instances of the Parameter class object. 
        This is simply done as a convenience wrapper. Note that SymPy symbols that are 
        suitable for developing model expressions using these variables are probably 
        more easily accessed using the methods:
        
        - get_symbol_for_t()    
        - get_symbol_for_p()
        
        and so on.  
                
        Returns
        -------
        Array of Parameter class instances, [Parameter,...]
        """
        return self._variables
    
    @variables.setter
    def variables(self, variables):
        self._variables = variables

    def get_symbol_for_t(self):
        """
        Retrieves SymPy symbol for temperature.

        Returns
        -------
        T : sympy.core.symbol.Symbol
            SymPy symbol for temperature
        """
        for x in self.variables:
            if x.name == 'T':
                return x.expression
        return None

    def get_symbol_for_p(self):
        """
        Retrieves SymPy symbol for pressure.

        Returns
        -------
        P : sympy symbol
            SymPy symbol for pressure
        """
        for x in self.variables:
            if x.name == 'P':
                return x.expression
        return None

    def get_symbol_for_v(self):
        """
        Retrieves SymPy symbol for volume.

        Returns
        -------
        V : sympy symbol
            SymPy symbol for volume
        """
        for x in self.variables:
            if x.name == 'V':
                return x.expression
        return None

    def get_symbol_for_tr(self):
        """
        Retrieves SymPy symbol for reference temperature.

        Returns
        -------
        Tr : sympy symbol
             SymPy symbol for reference temperature
        """
        for x in self.params:
            if x.name == 'T_r':
                return x.expression
        return None

    def get_symbol_for_pr(self):
        """
        Retrieves SymPy symbol for reference pressure.

        Returns
        -------
        Pr : sympy symbol
             SymPy symbol for reference pressure
        """
        for x in self.params:
            if x.name == 'P_r':
                return x.expression
        return None

    def get_symbol_for_vr(self):
        """
        Retrieves SymPy symbol for reference volume.

        Returns
        -------
        Vr : sympy symbol
             SymPy symbol for reference volume
        """
        for x in self.params:
            if x.name == 'V_r':
                return x.expression
        return None

    def get_module_name(self):
        """
        Retrieves module name.

        The module name is used in creating file names for coded functions.
        Maintained for backward compatibility.

        Returns
        -------
        string : str
            The name of the module.
        """
        return self._module

    def set_module_name(self, module_s):
        """
        Retrieves module name.

        The module name is used in creating file names for coded functions.
        Maintained for backward compatibility.

        Parameters
        ----------
        module_s : str
            The name of the module
        """
        self._module = module_s

    def set_reference_origin(self, Tr=298.15, Pr=1.0, Vr=2.5):
        """
        Sets the reference conditions for the model.

        Tr must be specified along with one of Pr or Vr.

        Parameters
        ----------
        Tr : float
             Reference temperature for the model (in Kelvins)
        Pr : float
             Reference pressure for the model (in bars)
        Vr : float
             Reference volume of the model (in J/bar)
        """
        self._T_r = Tr
        self._P_r = Pr
        self._V_r = Vr

    def get_include_debye_code(self):
        """
        Retrieves a boolean flag specifying whether a block of code 
        implementing the Debye integral function will be generated.

        Returns
        -------
        boolean : bool
                  True or False (default)
        """
        return self._include_debye_code

    def set_include_debye_code(self,include=False):
        """
        Sets a boolean flag controlling the inclusion of a block of code 
        implementing the Debye integral function.

        Parameters
        ----------
        include : bool
                  True or False
        """
        if include:
            self._include_debye_code = True
        else:
            self._include_debye_code = False

    def get_include_born_code(self):
        """
        Retrieves a boolean flag specifying whether code implementing the Born 
        functions will be linked to the generated module.

        Returns
        -------
        boolean : bool
                  True or False (default)
        """
        return self._include_born_code

    def set_include_born_code(self,include=False):
        """
        Sets a boolean flag controlling the inclusion of code implementing the 
        Born functions.

        Parameters
        ----------
        include : bool
                  True or False
        """
        if include:
            self._include_born_code = True
        else:
            self._include_born_code = False

    def get_model_param_names(self):
        """
        Retrieves a list of strings of model parameter names.

        Returns
        -------
        result : list of str
            Names of model parameters
        """
        result = []
        for x in self.params:
            result.append(x.name)
        return result

    def get_model_param_units(self):
        """
        Retrieves a list of strings of model parameter units.

        Returns
        -------
        result : list of str
            Units of model parameters
        """
        result = []
        for x in self.params:
            result.append(x.units)
        return result

    def get_model_param_symbols(self):
        """
        Retrieves a list of SymPy symbols for model parameters.

        Returns
        -------
        result : list of sym.Symbol
            SymPy symbols for model parameters
        """
        result = []
        for x in self.params:
            result.append(x.expression)
        return result

    def add_expression_to_model(self, expression, params, 
        exp_type='unrestricted',
        lower_limits=(None, None), upper_limits=(None, None), 
        implicit_functions=None, extend_restricted_functions=True):
        """
        Adds an expression and associated parameters to the model.

        Adds an expression for the Gibbs or Helmholtz energy to the 
        standard state model along with a description of expression parameters.

        Parameters
        ----------
        expression : sympy.core.symbol.Symbol
            A SymPy expression for the Gibbs free energy (if model_type is 'TP')
            or the Helmholtz energy (if model_type is 'TV').
            
            The expression may contain an implicit variable, *f*, whose value
            is a function of T,P or T,V.  The value of *f* is determined 
            numerically by solving the implicit_function expression (defined 
            below) at runtime using Newton's method.

        params : An array of tuples
            Structure (string, string, SymPy expression):
            
            - [0] str - Name of the parameter
            
            - [1] str - Units of parameter
            
            - [2] sympy.core.symbol.Symbol - SymPy symbol for the parameter

        exp_type : str, default='unrestricted'
            - 'unrestricted' - Expression is applicable over the whole of T,P or T,V space.
            
            - 'restricted' - Expression applies only between the specified 
              lower_limits and upper_limits.

        lower_limits : tuple
            A tuple of SymPy expressions defining the inclusive lower (T,P) or 
            (T,V) limit of applicability of the expression. Used only if 
            exp_type is set to 'restricted'.

        upper_limits : tuple
            A tuple of SymPy expressions defining the inclusive upper (T,P) or 
            (T,V) limit of applicability of the expression. Used only if 
            exp_type is set to 'restricted'.

        implicit_functions : array of tuples
            A tuple element contains three parts:
             
            - [0] sympy.core.symbol.Symbol - SymPy expression for an implicit 
              function in two independent (T,P or T,V) variables and one 
              dependent variable (*f*). The function is equal to zero.     
            
            - [1] sympy.core.symbol.Symbol - SymPy symbol for the dependent 
              variable *f*. 
            
            - [2] sympy.core.symbol.Symbol - SymPy expression that initializes *f* 
              in the iterative routine. This expression must be defined in 
              terms of known parameters and Tr, Pr, T, P.

        extend_restricted_functions : bool, default=True
            A boolean that controls whether a "restricted" *exp_type* is 
            extended beyond its *upper_limits*. By default, temperature and 
            pressure/volume derivatives of *expression* are evaluated at the
            *upper_limits*, and these constants are added as linear functions of
            *T* and *P/V* applicable and restricted to conditions *above* the 
            *upper_limits*. You can disable the default behavior by setting 
            *extend_restricted_functions* to *False*.  Note that, in general, the
            default behavior is desired as it allows entropic and volumetric 
            contributions developed over the restricted domain of *expression* 
            to contribute to the energy potential beyond the domain. Such 
            behavior is consistent with SymPy Piecewise functions. Special 
            circumstances, however, may require the default behavior to be
            disabled. 
        """
        l_params = []
        for (x,y,z) in params:
            l_params.append(Parameter(x,y,z))
        for x in l_params:
            self._params.append(x)

        l_exp = Expression(expression, l_params, 
            exp_type=exp_type, 
            lower_limit_T=lower_limits[0], lower_limit_PorV=lower_limits[1],
            upper_limit_T=upper_limits[0], upper_limit_PorV=upper_limits[1])
        self._expression_parts.append(l_exp)

        if implicit_functions:
            for (x,y,z) in implicit_functions:
                self._implicit_functions.append(Implicit_Function(x,y,z))

        if not extend_restricted_functions:
            return 

        if upper_limits[0] == None and upper_limits[1] == None:
            return

        elif upper_limits[0] != None and upper_limits[1] != None:
            if self._model_type == 'TP':
                T = self.get_symbol_for_t()
                P = self.get_symbol_for_p()
                s = -expression.diff(T).subs(T,upper_limits[0]).subs(
                    P,upper_limits[1])
                v =  expression.diff(P).subs(T,upper_limits[0]).subs(
                    P,upper_limits[1])
                g =  -(T-upper_limits[0])*s + (P-upper_limits[1])*v
                self.expression_parts.append(
                    Expression(g, l_params, exp_type='restricted',
                        lower_limit_T=upper_limits[0], 
                        lower_limit_PorV=upper_limits[1],
                        upper_limit_T=None,
                        upper_limit_PorV=None)
                    )
            elif self._model_type == 'TV':
                T = self.get_symbol_for_t()
                V = self.get_symbol_for_v()
                s = -expression.diff(T).subs(T,upper_limits[0]).subs(
                    V,upper_limits[1])
                p = -expression.diff(V).subs(T,upper_limits[0]).subs(
                    V,upper_limits[1])
                a =  -(T-upper_limits[0])*s - (V-upper_limits[1])*p
                self.expression_parts.append(
                    Expression(a, l_params, exp_type='restricted', 
                        lower_limit_T=upper_limits[0], 
                        lower_limit_PorV=upper_limits[1], 
                        upper_limit_T=None, 
                        upper_limit_PorV=None)
                    )
        elif upper_limits[0] != None:
            if self._model_type == 'TP':
                T = self.get_symbol_for_t()
                s = -expression.diff(T).subs(T,upper_limits[0])
                g =  -(T-upper_limits[0])*s
                self.expression_parts.append(
                    Expression(g, l_params, exp_type='restricted', 
                        lower_limit_T=upper_limits[0], 
                        lower_limit_PorV=upper_limits[1], 
                        upper_limit_T=None, 
                        upper_limit_PorV=None)
                    )
            elif self._model_type == 'TV':
                T = self.get_symbol_for_t()
                s = -expression.diff(T).subs(T,upper_limits[0])
                a =  -(T-upper_limits[0])*s
                self.expression_parts.append(
                    Expression(a, l_params, exp_type='restricted', 
                        lower_limit_T=upper_limits[0], 
                        lower_limit_PorV=upper_limits[1], 
                        upper_limit_T=None, 
                        upper_limit_PorV=None)
                    )
        elif upper_limits[1] != None:
            if self._model_type == 'TP':
                P = self.get_symbol_for_p()
                v =  expression.diff(P).subs(P,upper_limits[1])
                g =  (P-upper_limits[1])*v
                self.expression_parts.append(
                    Expression(g, l_params, exp_type='restricted', 
                        lower_limit_T=upper_limits[0], 
                        lower_limit_PorV=upper_limits[1], 
                        upper_limit_T=None, 
                        upper_limit_PorV=None)
                    )
            elif self._model_type == 'TV':
                V = self.get_symbol_for_v()
                p = -expression.diff(V).subs(V,upper_limits[1])
                a = -(V-upper_limits[1])*p
                self.expression_parts.append(
                    Expression(a, l_params, exp_type='restricted', 
                        lower_limit_T=upper_limits[0], 
                        lower_limit_PorV=upper_limits[1], 
                        upper_limit_T=None, 
                        upper_limit_PorV=None)
                    )

    def create_calc_h_file(self, language='C', module_type='fast'):
        """
        Creates an include file implementing model calculations.

        Note that this include file contains code for functions that implement 
        the model for the generic case. It is meant to be included into a file 
        that implements a specific parameterized instance of the model. See 
        create_fast_code_module().

        The user does not normally call this function directly.

        Parameters
        ----------
        language : str
            Language syntax for generated code. ("C" is the C99 programming 
            language.)
        module_type : str
            Generate code that executes "fast" but does not expose hooks for 
            model parameter calibration. Alternately, generate code suitable for 
            "calib"ration of parameters in the model, which executes more slowly 
            and exposes additional functions that allow setting of parameters 
            and generation of derivatives of thermodynamic functions with 
            respect to model parameters. 
        """
        assert language == "C"
        s = "#include <math.h>\n\n"
        if self._include_debye_code:
            s += tpl.create_code_for_debye_function(language)
        if self._include_born_code:
            s += tpl.create_code_for_born_functions(language)
            s += "\n"
        if self._model_type == 'TP':
            T = self.get_symbol_for_t()
            P = self.get_symbol_for_p()
            implicit_calls = ""
            # list of lists [][], rows denote temperature derivatives, columns, 
            # pressure derivatives
            # each element contains a tuple: (function, symbol for function, 
            # sympy expression for function,
            # sympy expression for fully substituted function)
            self.global_subs = []
            if len(self._implicit_functions) > 0:
                for index,imp in enumerate(self._implicit_functions):
                    fn       = imp.function
                    var      = imp.variable
                    var_init = imp.guess
                    var_str  = "var" + str(index)
                    var_sym  = sym.symbols(var_str)
                    fn_subs  = fn.subs(var,var_sym)
                    if module_type == 'fast':
                        s += ("static double " + var_str + " = " 
                            + self.printer.doprint(var_init) + ";\n")
                    elif module_type == 'calib':
                        s += "static double " + var_str + " = -1.0;\n"
                    implicit_calls += ("    " + self.module + "_" + var_str 
                        + "(T, P);\n")
                    for oT in range(0,4):
                        self.global_subs.append([])
                        for oP in range(0,4):
                            if oT+oP == 0:
                                # definition of the primary implicit variable
                                self.global_subs[0].append(
                                    (var, var_str, var_str, var))
                            elif oT+oP <= 3:
                                var_diff = sym.diff(var,T,oT,P,oP)
                                fn_diff = sym.solve(sym.diff(fn,T,oT,P,oP),
                                    var_diff)[0]
                                var_sym_diff = (sym.symbols("d"+str(oT+oP)+"var"
                                    +str(index)+"d"+str(oT)+"Td"+str(oP)+"P"))
                                fn_sym_diff = fn_diff
                                ss_list = []
                                for oo in self.global_subs:
                                    for (aa,bb,cc,dd) in oo:
                                        ss_list.insert(0,(aa,cc))
                                fn_sym_diff = fn_sym_diff.subs(ss_list)
                                self.global_subs[oT].append(
                                    (var_diff, fn_diff, var_sym_diff, 
                                        fn_sym_diff))
                    s += ("static void " + self.module + "_" + var_str 
                        + "(double T, double P) {\n")
                    s += "    static double Told = 0.0;\n"
                    s += "    static double Pold = 0.0;\n"
                    if module_type == 'calib':
                        s += "    if (" + var_str + " == -1.0) {\n"
                        s += ("        " + var_str + " = " 
                            + self.printer.doprint(var_init) + ";\n")
                        s += "    }\n"
                    s += "    if ((T != Told) && (P != Pold)) {\n"
                    s += "        Told = T;\n"
                    s += "        Pold = P;\n"
                    s += "        double f = 0.0;\n"
                    s += "        int iter = 0;\n"
                    s += "        do {\n"
                    s += "            f = " + self.printer.doprint(fn_subs) + ";\n"
                    s += "            double df = " 
                    s += self.printer.doprint(sym.diff(fn_subs,var_sym)) + ";\n"
                    s += "            if (df == 0.0) break;\n"
                    s += "            " + var_str + " -= f/df;\n"
                    s += "            if (" + var_str + " <= 0.0) " 
                    s += var_str + " = 0.001;\n"
                    s += "            else if (" + var_str + " >= 2.0*V_TrPr) " 
                    s += var_str + " = 2.0*V_TrPr;\n"
                    s += "            iter++;\n"
                    s += "        } while ((fabs(f) > 0.001) && (iter < 200));\n"
                    s += "    }\n"
                    s += "}\n\n"
            for (f, oT, oP) in self.g_list:
                s += ("static double " + self.module + "_" + f 
                    + "(double T, double P) {\n")
                s += implicit_calls
                s += "    double result = 0.0;\n"
                #
                ss_list = []
                if len(self.global_subs) > 0:
                    ss_index = list(product([i for i in range(0,oT+1)],
                        [i for i in range(0,oP+1)]))
                    for (ooT,ooP) in ss_index:
                        if ooT+ooP > 0:
                            s += "    double " 
                            s += self.printer.doprint(self.global_subs[ooT][ooP][2]) 
                            s += " = "
                            s += self.printer.doprint(self.global_subs[ooT][ooP][3]) 
                            s += ";\n"
                        ss_list.insert(0,(self.global_subs[ooT][ooP][0], 
                            self.global_subs[ooT][ooP][2]))
                #
                inited = False
                for exp in self.expression_parts:
                    a = sym.diff(exp.expression,T,oT,P,oP).subs(ss_list)
                    if a == sym.S.Zero:
                        pass
                    elif exp.exp_type == 'unrestricted':
                        s += "    result += " + self.printer.doprint(a) + ";\n"
                        inited = True
                    elif exp.exp_type == 'restricted':
                        test_term = "    if("
                        separator = ""
                        for param in exp.params:
                            test_term += (separator + "((" 
                                       + self.printer.doprint(param.expression) 
                                       + ") != 0.0)")
                            separator = " || "
                        test_term += ") {\n"
                        s += test_term
                        test_term = "        if("
                        separator = ""
                        if exp.lower_limit_T != None:
                            test_term += (separator + "(T > (" 
                                + self.printer.doprint(exp.lower_limit_T) 
                                + "))")
                            separator = " && "
                        if exp.lower_limit_PorV != None:
                            test_term += (separator + "(P > (" 
                                + self.printer.doprint(exp.lower_limit_PorV) 
                                + "))")
                            separator = " && "
                        if exp.upper_limit_T != None:
                            test_term += (separator + "(T <= (" 
                                + self.printer.doprint(exp.upper_limit_T) 
                                + "))")
                            separator = " && "
                        if exp.upper_limit_PorV != None:
                            test_term += (separator + "(P <= (" 
                                + self.printer.doprint(exp.upper_limit_PorV) 
                                + "))")
                            separator = " && "
                        test_term += ") {\n"
                        s += test_term
                        s += "             result += " 
                        s += self.printer.doprint(a) + ";\n"
                        s += "        }\n"
                        s += "    }\n"
                        inited = True
                if not inited:
                    s += "    result += 0.0;\n"
                s += "    return result;\n}\n\n"
        elif self._model_type == 'TV':
            T = self.get_symbol_for_t()
            V = self.get_symbol_for_v()
            implicit_calls = ""
            # list of lists [][], rows denote temperature derivatives, columns, 
            # volume derivatives
            # each element contains a tuple: (function, symbol for function, 
            # sympy expression for function,
            # sympy expression for fully substituted function)
            self.global_subs = []
            if len(self._implicit_functions) > 0:
                for index,imp in enumerate(self._implicit_functions):
                    fn       = imp.function
                    var      = imp.variable 
                    var_init = imp.guess
                    var_str = "var" + str(index)
                    var_sym = sym.symbols(var_str)
                    fn_subs = fn.subs(var,var_sym)
                    if module_type == 'fast':
                        s += ("static double " + var_str + " = " 
                            + self.printer.doprint(var_init) + ";\n")
                    elif module_type == 'calib':
                        s += "static double " + var_str + " = -1.0;\n"
                    implicit_calls += ("    " + self.module + "_" + var_str 
                        + "(T, V);\n")
                    for oT in range(0,4):
                        self.global_subs.append([])
                        for oV in range(0,4):
                            if oT+oV == 0:
                                # definition of the primary implicit variable
                                self.global_subs[0].append((var, var_str, 
                                    var_str, var))
                            elif oT+oV <= 3:
                                var_diff = sym.diff(var,T,oT,V,oV)
                                fn_diff = sym.solve(sym.diff(fn,T,oT,V,oV),
                                    var_diff)[0]
                                var_sym_diff = sym.symbols("d"+str(oT+oV)+"var"
                                    +str(index)+"d"+str(oT)+"Td"+str(oV)+"V")
                                fn_sym_diff = fn_diff
                                ss_list = []
                                for oo in self.global_subs:
                                    for (aa,bb,cc,dd) in oo:
                                        ss_list.insert(0,(aa,cc))
                                fn_sym_diff = fn_sym_diff.subs(ss_list)
                                self.global_subs[oT].append((var_diff, fn_diff, 
                                    var_sym_diff, fn_sym_diff))
                    s += "static void " + self.module + "_" 
                    s += var_str + "(double T, double V) {\n"
                    s += "    static double Told = 0.0;\n"
                    s += "    static double Vold = 0.0;\n"
                    if module_type == 'calib':
                        s += "    if (" + var_str + " == -1.0) {\n"
                        s += "        " + var_str + " = " 
                        s += self.printer.doprint(var_init) + ";\n"
                        s += "    }\n"
                    s += "    if ((T != Told) && (V != Vold)) {\n"
                    s += "        Told = T;\n"
                    s += "        Vold = V;\n"
                    s += "        double f = 0.0;\n"
                    s += "        int iter = 0;\n"
                    s += "        do {\n"
                    s += "            f = " + self.printer.doprint(fn_subs) + ";\n"
                    s += "            double df = " 
                    s += self.printer.doprint(sym.diff(fn_subs,var_sym)) + ";\n"
                    s += "            if (df == 0.0) break;\n"
                    s += "            " + var_str + " -= f/df;\n"
                    s += "            if (" + var_str + " <= 0.0) " 
                    s += var_str + " = 0.001;\n"
                    s += "            else if (" + var_str + " >= 4.0*" 
                    s += self.printer.doprint(var_init) 
                    s +=                  +") " + var_str + " = 4.0*" 
                    s += self.printer.doprint(var_init) + ";\n"
                    s += "            iter++;\n"
                    s += "        } while ((fabs(f) > 0.001) && (iter < 200));\n"
                    s += "    }\n"
                    s += "}\n\n"
            for (f, oT, oV) in self.a_list:
                s += ("static double " + self.module + "_" + f 
                    + "(double T, double V) {\n")
                s += implicit_calls
                s += "    double result = 0.0;\n"
                #
                ss_list = []
                if len(self.global_subs) > 0:
                    ss_index = list(product([i for i in range(0,oT+1)],
                        [i for i in range(0,oV+1)]))
                    for (ooT,ooV) in ss_index:
                        if ooT+ooV > 0:
                            s += "    double " 
                            s += self.printer.doprint(self.global_subs[ooT][ooV][2]) 
                            s += " = "
                            s += self.printer.doprint(self.global_subs[ooT][ooV][3]) 
                            s += ";\n"
                        ss_list.insert(0,(self.global_subs[ooT][ooV][0], 
                            self.global_subs[ooT][ooV][2]))
                #
                inited = False
                for exp in self.expression_parts:
                    a = sym.diff(exp.expression,T,oT,V,oV).subs(ss_list)
                    if a == sym.S.Zero:
                        pass
                    elif exp.exp_type == 'unrestricted':
                        s += "    result += " + self.printer.doprint(a) + ";\n"
                        inited = True
                    elif exp.exp_type == 'restricted':
                        test_term = "    if("
                        separator = ""
                        for param in params:
                            test_term += (separator + "((" 
                                + self.printer.doprint(param.expression) 
                                + ") != 0.0)")
                            separator = " || "
                        test_term += ") {\n"
                        s += test_term
                        test_term = "        if("
                        separator = ""
                        if exp.lower_limit_T != None:
                            test_term += (separator + "(T > (" 
                                + self.printer.doprint(exp.lower_limit_T) 
                                + "))")
                            separator = " && "
                        if exp.lower_limit_PorV != None:
                            test_term += (separator + "(V > (" 
                                + self.printer.doprint(exp.lower_limit_PorV) 
                                + "))")
                            separator = " && "
                        if exp.upper_limit_T != None:
                            test_term += (separator + "(T <= (" 
                                + self.printer.doprint(exp.upper_limit_T) 
                                + "))")
                            separator = " && "
                        if exp.upper_limit_PorV != None:
                            test_term += (separator + "(V <= (" 
                                + self.printer.doprint(exp.upper_limit_PorV) 
                                + "))")
                            separator = " && "
                        test_term += ") {\n"
                        s += test_term
                        s += "             result += " + self.printer.doprint(a) 
                        s += ";\n"
                        s += "        }\n"
                        s += "    }\n"
                        inited = True
                if not inited:
                    s += "    result += 0.0;\n"
                s += "    return result;\n}\n\n"
        else:
            print ("Unsupported model_type: ", self._model_type)
            s = ""
        self.calc_h = (s + tpl.create_redundant_function_template(
            model_type=self._model_type).format(module=self.module,
            v_initial_guess="2.0"))

    def create_calib_h_file(self, language='C'):
        """
        Creates a C-code include file implementing model calibration functions.

        Parameters
        ----------
        language : string
            Language syntax for generated code. ("C" is the C99 programming 
            language.)

        Notes
        -----
        The calib_h include file contains code for functions that implement 
        the model for the generic case. It is meant to be included into a file 
        that implements a specific parameterized instance of the model. 
        See create_calib_code_module().

        The user does not normally call this function directly.
        """
        assert language == "C"
        s = "#include <math.h>\n\n"
        if self._model_type == 'TP':
            T = self.get_symbol_for_t()
            P = self.get_symbol_for_p()
            for (f, oT, oP) in self.g_list:
                s += ("static double " + self.module + "_dparam_" + f 
                    + "(double T, double P, int index) {\n")
                s += "    double result = 0.0;\n"
                #
                ss_list = []
                if len(self.global_subs) > 0:
                    ss_index = list(product([i for i in range(0,oT+1)],
                        [i for i in range(0,oP+1)]))
                    for (ooT,ooP) in ss_index:
                        if ooT+ooP > 0:
                            s += "    double " 
                            s += self.printer.doprint(self.global_subs[ooT][ooP][2]) 
                            s += " = "
                            s += self.printer.doprint(self.global_subs[ooT][ooP][3]) 
                            s += ";\n"
                        ss_list.insert(0,(self.global_subs[ooT][ooP][0], 
                            self.global_subs[ooT][ooP][2]))
                #
                s += "    switch (index) {\n"
                for i in range (0,len(self.params)):
                    s += ("    case " + str(i) + ": /* " 
                        + self.printer.doprint(self.params[i].expression) 
                        + " */ \n")
                    #
                    inited = False
                    for exp in self.expression_parts:
                        a = sym.diff(exp.expression,T,oT,P,oP,
                            self.params[i].expression,1).subs(ss_list)
                        if a == sym.S.Zero:
                            pass
                        elif exp.exp_type == 'unrestricted':
                            s += "        result += " 
                            s += self.printer.doprint(a) + ";\n"
                            inited = True
                        elif exp.exp_type == 'restricted':
                            test_term = "        if("
                            separator = ""
                            for param in exp.params:
                                test_term += (separator + "((" 
                                    + self.printer.doprint(param.expression) 
                                    + ") != 0.0)")
                                separator = " || "
                            test_term += ") {\n"
                            s += test_term
                            test_term = "            if("
                            separator = ""
                            if exp.lower_limit_T != None:
                                test_term += (separator + "(T > (" 
                                    + self.printer.doprint(exp.lower_limit_T) 
                                    + "))")
                                separator = " && "
                            if exp.lower_limit_PorV != None:
                                test_term += (separator + "(P > (" 
                                    + self.printer.doprint(exp.lower_limit_PorV) 
                                    + "))")
                                separator = " && "
                            if exp.upper_limit_T != None:
                                test_term += (separator + "(T <= (" 
                                    + self.printer.doprint(exp.upper_limit_T) 
                                    + "))")
                                separator = " && "
                            if exp.upper_limit_PorV != None:
                                test_term += (separator + "(P <= (" 
                                    + self.printer.doprint(exp.upper_limit_PorV) 
                                    + "))")
                                separator = " && "
                            test_term += ") {\n"
                            s += test_term
                            s += "                 result += " 
                            s += self.printer.doprint(a) + ";\n"
                            s += "            }\n"
                            s += "        }\n"
                            inited = True
                    #
                    if not inited:
                        s += "        result += 0.0;\n"
                    s += "        break;\n"
                s += "    }\n"
                s += "    return result;\n}\n\n"
        elif self._model_type == 'TV':
            T = self.get_symbol_for_t()
            V = self.get_symbol_for_v()
            for (f, oT, oV) in self.a_list:
                s += ("static double " + self.module + "_dparam_" + f 
                    + "(double T, double V, int index) {\n")
                s += "    double result = 0.0;\n"
                #
                ss_list = []
                if len(self.global_subs) > 0:
                    ss_index = list(product([i for i in range(0,oT+1)],
                        [i for i in range(0,oV+1)]))
                    for (ooT,ooV) in ss_index:
                        if ooT+ooV > 0:
                            s += "    double " 
                            s += self.printer.doprint(self.global_subs[ooT][ooV][2])
                            s += " = "
                            s += self.printer.doprint(self.global_subs[ooT][ooV][3]) 
                            s += ";\n"
                        ss_list.insert(0,(self.global_subs[ooT][ooV][0], 
                            self.global_subs[ooT][ooV][2]))
                #
                s += "    switch (index) {\n"
                for i in range (0,len(self.params)):
                    s += "    case " + str(i) + ": /* " 
                    s += self.printer.doprint(self.params[i].expression) 
                    s += " */ \n"
                    #
                    inited = False
                    for exp in self.expression_parts:
                        a = sym.diff(exp.expression,T,oT,V,oV,
                            self.params[i].expression,1).subs(ss_list)
                        if a == sym.S.Zero:
                            pass
                        elif exp.exp_type == 'unrestricted':
                            s += "        result += " + self.printer.doprint(a) 
                            s += ";\n"
                            inited = True
                        elif exp.exp_type == 'restricted':
                            test_term = "        if("
                            separator = ""
                            for param in params:
                                test_term += (separator + "((" 
                                    + self.printer.doprint(param.expression) 
                                    + ") != 0.0)")
                                separator = " || "
                            test_term += ") {\n"
                            s += test_term
                            test_term = "            if("
                            separator = ""
                            if exp.lower_limit_T != None:
                                test_term += (separator + "(T > (" 
                                    + self.printer.doprint(exp.lower_limit_T) 
                                    + "))")
                                separator = " && "
                            if exp.lower_limit_PorV != None:
                                test_term += (separator + "(V > (" 
                                    + self.printer.doprint(exp.lower_limit_PorV) 
                                    + "))")
                                separator = " && "
                            if exp.upper_limit_T != None:
                                test_term += (separator + "(T <= (" 
                                    + self.printer.doprint(exp.upper_limit_T) 
                                    + "))")
                                separator = " && "
                            if exp.upper_limit_PorV != None:
                                test_term += (separator + "(V <= (" 
                                    + self.printer.doprint(exp.upper_limit_PorV) 
                                    + "))")
                                separator = " && "
                            test_term += ") {\n"
                            s += test_term
                            s += "                 result += " 
                            s += self.printer.doprint(a) + ";\n"
                            s += "            }\n"
                            s += "        }\n"
                            inited = True
                    #
                    if not inited:
                        s += "        result += 0.0;\n"
                    s += "        break;\n"
                s += "    }\n"
                s += "    return result;\n}\n\n"
            s += tpl.create_redundant_calib_TV_template().format(
                module=self.module)
        else:
            print ("Unsupported model_type: ", self._model_type)
            s = ""
        
        # Additional functions
        s += "static int " + self.module + "_get_param_number(void) {\n"
        s += "    return " + str(len(self.params)) + ";\n"
        s += "}\n\n"
        #
        s += "static const char *paramNames["  + str(len(self.params)) + "] = {"
        separator = " "
        for param in self.params:
            s += separator + '"' + self.printer.doprint(param.expression) + '"'
            separator = ", "
        s += "  };\n\n"
        #
        s += "static const char *paramUnits["  + str(len(self.params)) + "] = {"
        separator = " "
        for param in self.params:
            s += separator + '"' + param.units + '"'
            separator = ", "
        s += "  };\n\n"
        #
        s += "static const char **" + self.module + "_get_param_names(void) {\n"
        s += "    return paramNames;\n"
        s += "}\n\n"
        #
        s += "static const char **" + self.module + "_get_param_units(void) {\n"
        s += "    return paramUnits;\n"
        s += "}\n\n"
        #
        s += "static void " + self.module + "_get_param_values(double **values) {\n"
        for i in range(0, len(self.params)):
            s += "    (*values)[" + str(i) + "] = " 
            s += self.printer.doprint(self.params[i].expression) + ";\n"
        s += "}\n\n"
        #
        s += "static int " + self.module + "_set_param_values(double *values) {\n"
        for i in range(0, len(self.params)):
            s += "    " + self.printer.doprint(self.params[i].expression) 
            s += "= values[" + str(i) + "];\n"
        s += "    return 1;\n"
        s += "}\n\n"
        #
        s += "static double " + self.module + "_get_param_value(int index) {\n"
        s += "    double result = 0.0;\n"
        s += "    switch (index) {\n"
        for i in range(0, len(self.params)):
            s += "    case " + str(i) + ":\n"
            s += "        result = " 
            s += self.printer.doprint(self.params[i].expression) + ";\n"
            s += "        break;\n"
        s += "     default:\n"
        s += "         break;\n"
        s += "    }\n"
        s += "    return result;\n"
        s += "}\n\n"
        #
        s += "static int " + self.module 
        s += "_set_param_value(int index, double value) {\n"
        s += "    int result = 1;\n"
        s += "    switch (index) {\n"
        for i in range(0, len(self.params)):
            s += "    case " + str(i) + ":\n"
            s += "        " + self.printer.doprint(self.params[i].expression) 
            s += " = value;\n"
            s += "        break;\n"
        s += "     default:\n"
        s += "         break;\n"
        s += "    }\n"
        s += "    return result;\n"
        s += "}\n\n"
        #
        self.calib_h = s

    def parse_formula(self, formula_string="Si(1)O(2)"):
        """
        Parses a chemical formula, and returns a molecular weight and an array of
        elemental concentrations.

        Parameters
        ----------
        formula_string : str 
            Formula of compound specified in the standard form
            SiO2 -> Si(1)O(2)

        Returns
        -------
        mw, elmvector : tuple
            - [0] mw is a float containing the molecular weight in grams.
                    
            - [1] elmvector is a numpy array of length 120 containing the mole numbers 
              of each element in the compound.
        """
        formula = formula_string.title().replace('(',',').replace(')',',').split(',')
        mw = 0.0
        elmvector = np.zeros(120)
        for i in range(0,len(formula)-1,2):
            element = self._elements.loc[self._elements['Abbrv'] == formula[i]]
            ind = element.index.values[0]
            mw += element.MW.values[0]*float(formula[i+1])
            elmvector[ind] = float(formula[i+1])
        return mw, elmvector
    
    def get_berman_std_state_database(self, identifier=None, extend_defs=False):
        """
        Retrieves priors from the thermodynamic database of Berman (1988).

        Parameters
        ----------
        identifier : int
            Value may be None or an integer in the range (0, length of Berman 
            database).
        extend_defs : bool
            False: The standard set of Berman parameters are retrieved.  
            
            True: Additional parameters are computed, including 'K' and 
            'K_P', the bulk modulus and its pressure derivative.

        Returns
        -------
        result : variable
            If identifier is None, then the function returns an array of tuples, 
            one for each phase in the database:
            
            - [0] database index  
            
            - [1] name of the phase  
            
            - [2] formula of the phase, in standard notation  
            
            If identifier is an integer, then the function returns a dictionary 
            with keys, corresponding parameter names, and values corresponding 
            to parameter values.
        
        """
        if self._berman_db is None:
            params = ['H_TrPr', 'S_TrPr', 'k0', 'k1', 'k2', 'k3', 'V_TrPr', 
            'v1', 'v2', 'v3', 'v4', 'l1', 'l2', 'k_lambda', 'T_lambda_Pr', 
            'T_lambda_ref', 'H_t', 'd0', 'd1', 'd2', 'd3', 'd4', 'd5', 'T_D', 
            'T_D_ref', 'T_r', 'P_r']
            self._berman_db = pd.read_json(path.join(path.dirname(__file__), 
                DATADIR, 'berman_1988.json'))
            self._berman_db['T_r'] = 298.15 
            self._berman_db['P_r'] = 1.0
            if extend_defs:
                params.append('K')
                self._berman_db['K'] = self._berman_db.apply(
                    lambda row: -1.0/row['v1'] if row['v1'] != 0 else 1000000.0, 
                    axis=1)
                params.append('K_P')
                self._berman_db['K_P'] = self._berman_db.apply(
                    lambda row: (2.0*row['v2']/row['v1']/row['v1']-1.0) \
                    if row['v1'] != 0 and row['v2'] != 0 else 4.0, axis=1)
            self._berman_db = self._berman_db[['Phase', 'Formula'] + params]
            self._berman_db.fillna(0, inplace=True)
        if identifier is None:
            result = []
            for i in range (0, len(self._berman_db)):
                row = self._berman_db.iloc[i,:]
                result.append((i, row.Phase.title(), 
                    row.Formula.title().replace("(1)","").replace("(", "").replace(")","")))
        elif identifier >= 0 and identifier < len(self._berman_db):
            row = self._berman_db.to_dict('records')[identifier]
            result = row
        else:
            result = None
        return result

    def create_code_module(self, phase="Quartz", formula="Si(1)O(2)", params={}, 
        identifier=None, prefix="cy", module_type="fast", silent=False, 
        language='C'):
        """
        Creates and writes an include and code file for a model instance.

        Parameters
        ----------
        phase : str
            Model instance title (e.g., phase name).  Used to name the generated 
            function. Cannot contain blank spaces or special characters; underscore 
            ("_") is permitted. Convention capitalizes the first letter and the 
            letter following an underscore ("_") character.
        formula : str
            Chemical formula of the model instance, in standard notation.
            Standard notation is of the form: Element Symbol followed by 
            parentheses enclosing a number, which may be decimalized.  E.g., 
            SiO2 is written as Si(1)O(2); CaMg1/2Ti1/2AlO6 is written as 
            Ca(1)Mg(0.5)Ti(0.5)Al(1)Si(1)O(6); CaMg(CO3)2 is written as 
            Ca(1)Mg(1)C(2)O(6).
        params : dict
            Parameter values for the model instance.
            The keys of this dictionary are validated against parameter symbols 
            stored for the model.
        identifier : str
            A unique identifier for the model instance.
            Defaults to local date and time when module is created (rounded to 
            the second).
        prefix : str
            Prefix to function names for Python bindings, e.g., 
            {prefix}_{phase}_{module}_g(T,P).
        module_type : str
            Generate code that executes "fast" but does not expose hooks for 
            model parameter calibration. Alternately, generate code suitable for 
            "calib"ration of parameters in the model, which executes more slowly 
            and exposes additional functions that allow setting of parameters 
            and generation of derivatives of thermodynamic functions with 
            respect to model parameters. 
        silent : bool
            Print (True) or do not print (False) status messages.
        language : string
            Language syntax for generated code. ("C" is the C99 programming 
            language; "C++" is the C++ programming language.)

        Returns
        -------
        result : Boolean
                 True if module is succesfully generated, False if some error occurred 
        """
        assert language == "C"
        if not re.match("^[a-zA-Z0-9_]*$", phase):
            print ("Error: ", phase, 
                " is only allowed to have characters a-z, A-Z, 0-9, and _")
            return False
        okay = True
        for param in self.params:
            if not param.name in params:
                print ("Error: paramter key ", param.name, 
                    " is missing from passed model parameter dictionary")
                okay = False
        if not okay:
            print ("Error in params dictionary (see above)")
            return False
        (mw, elmvector) = self.parse_formula(formula)
        if mw == 0.0:
            print ("Error: parsing of the formula ", formula, 
                " retunred a molecular weight of zero")
            return False
        module = self.get_module_name()
        if module == 'untitled':
            print ("Error: Please set module name before called this function")
            return False

        if self.calc_h is None:
            if not silent:
                print ("Creating (once only) generic fast model code file string")
            self.create_calc_h_file(language, module_type)

        if module_type == 'fast':
            if self.fast_h_template is None:
                if not silent:
                    print ("Creating (once only) generic model fast code template include file string")
                self.fast_h_template = tpl.create_fast_h_template()
            if self.fast_c_template is None:
                if not silent:
                    print ("Creating (once only) generic model fast code template code file string")
                self.fast_c_template = tpl.create_fast_c_template()
        elif module_type == 'calib':
            if self.calib_h_template is None:
                if not silent:
                    print ("Creating (once only) generic model calib code template include file string")
                self.calib_h_template = tpl.create_calib_h_template()
            if self.calib_c_template is None:
                if not silent:
                    print ("Creating (once only) generic model calib code template code file string")
                self.calib_c_template = tpl.create_calib_c_template()
            if self.calib_h is None:
                if not silent:
                    print ("Creating (once only) generic calib model code file string")
                self.create_calib_h_file()
        else:
            print ("Error: module_type must be set to either 'fast' or 'calib'")
            return False


        if not silent:
            print ("Creating include file ...")
        h_file = ""
        h_file_name = ""
        if module_type == 'fast':
            h_file = self.fast_h_template.format(module=module, phase=phase)
            h_file_name = phase + "_" + module + "_calc.h"
        elif module_type == 'calib':
            h_file = self.calib_h_template.format(module=module, phase=phase)
            h_file_name = phase + "_" + module + "_calib.h"
        if not silent:
            print ("... done!")

        if not silent:
            print ("Creating code file ...")
        if identifier is None:
            identifier = time.asctime(time.localtime(time.time()))
        formula_condensed = formula.title().replace("(1)","").replace("(", "").replace(")","") 
        param_block = "\n"
        for param in self.params:
            if module_type == 'fast':
                param_block += "static const double " + param.name 
                param_block += " = " + str(params[param.name]) + ";\n"
            elif module_type == 'calib':
                param_block += "static double " + param.name 
                param_block += " = " + str(params[param.name]) + ";\n"
        (param_names, param_values) = zip(*params.items())
        c_file = ""
        c_file_name = ""
        if module_type == 'fast':
            c_file = self.fast_c_template.format( module=module, phase=phase, 
                formula=formula_condensed, mw=mw, elmvector=elmvector, 
                parameter_init_block=param_block, git_identifier=identifier,
                include_calc_h=self.calc_h)
            c_file_name = phase + "_" + module + "_calc.c"
        elif module_type == 'calib':
            c_file = self.calib_c_template.format( module=module, phase=phase, 
                formula=formula_condensed, mw=mw, elmvector=elmvector, 
                parameter_init_block=param_block, git_identifier=identifier,
                include_calc_h=self.calc_h, include_calib_h=self.calib_h)
            c_file_name = phase + "_" + module + "_calib.c"
        if not silent:
            print ("... done")

        if not silent:
            print ("Writing include file to working directory ...")
        with open(h_file_name, 'w') as f:
            f.write(h_file)

        if not silent:
            print ("Writing code file to working directory ...")
        with open(c_file_name, 'w') as f:
            f.write(c_file)

        if not silent:
            print ("Writing pyxbld file to working directory ...")
        pyxbld_template = tpl.create_pyxbld_template()
        pyxbld_file = pyxbld_template.format(file_to_compile=c_file_name)
        pyxbld_file_name = module + ".pyxbld"
        with open(pyxbld_file_name, 'w') as f:
            f.write(pyxbld_file)

        if not silent:
            print ("writing pyx file to working directory ...")
        pyx_template = ""
        if module_type == "fast":
            pyx_template = tpl.create_fast_pyx_template()
        elif module_type == "calib":
            pyx_template = tpl.create_calib_pyx_template()
        pyx_file = pyx_template.format(prefix=prefix, phase=phase, module=module)
        pyx_file_name = module + ".pyx"
        with open(pyx_file_name, 'w') as f:
            f.write(pyx_file)

        # pyximport.install(pyximport=True, pyimport=False, build_dir=None, 
        # build_in_temp=True, setup_args=None, reload_support=False, 
        # load_py_module_on_import_failure=False, inplace=False, 
        # language_level=None)

        if not silent:
            print ("Compiling code and Python bindings ...")
        pyximport.install(language_level=3)

        if not silent:
            print ("Success! Import the module named ", module)
        return True

class SimpleSolnModel:
    """
    Class creates a model of the thermodynamic properties of a simple solution.

    Parameters
    ----------
    nc : int
        Number of thermodynamic components in the solution
    model_type : str
        Model type of 'TP' implies that expressions are of the form of the Gibbs 
        free energy.
        
        Model type of 'TV' implies that expressions are of the form of the 
        Helmholtz energy. This option is infrequently used.

    Attributes
    ----------
    a_list  
    conversion_string
    d2n_g_list  
    d3n_g_list  
    dn_g_list  
    expression_parts
    formula_string  
    g_list  
    implicit_functions  
    module  
    mu  
    n  
    nc  
    nT  
    params  
    printer  
    test_string  
    variables 
    
    Methods
    -------
    add_expression_to_model
    create_code_module
    create_conversion_code_block
    create_formula_code_block
    create_soln_calc_h_file
    create_soln_calib_h_file
    create_test_code_block
    get_include_dh_code
    get_model_param_names
    get_model_param_symbols
    get_model_param_units
    get_symbol_for_p
    get_symbol_for_pr
    get_symbol_for_t
    get_symbol_for_tr
    get_symbol_for_v
    get_symbol_for_vr
    set_include_dh_code
    set_reference_origin
    symmetric_index_from_2d_array
    symmetric_index_from_3d_array
    
    Notes
    -----
    This class is for creating a representation of the thermodynamic properties 
    of a simple solution phase.  The principal use of the class is to construct 
    a model symbolically and to generate code that implements the model for 
    model parameter calibration and thermodynamic property calculation.

    A simple solution is one that does not contain implicit variables that need
    be determined via solution of conditions of homogeneous equilibrium. 
    Examples of simple solutions include regular solutions, asymmetric binary 
    and ternary Margules expansions, and high order Taylor series expansions of 
    the Gibbs free energy of mixing.  Examples of solutions not compatible with
    this class are aqueous or liquid solutions involving complex formation 
    (speciation) or solid solutions that include cation-ordering or composition
    dependent symmetry-breaking phase transitions.

    """
    def __init__(self, nc=None, model_type='TP'):
        assert nc
        self._nc = nc
        self._g_list = [('g',0,0), ('dgdt',1,0), ('dgdp',0,1), ('d2gdt2',2,0), 
                        ('d2gdtdp',1,1), ('d2gdp2',0,2), ('d3gdt3',3,0), 
                        ('d3gdt2dp',2,1), ('d3gdtdp2',1,2), ('d3gdp3',0,3)]
        self._a_list = [('a',0,0), ('dadt',1,0), ('dadv',0,1), ('d2adt2',2,0), 
                        ('d2adtdv',1,1), ('d2adv2',0,2), ('d3adt3',3,0), 
                        ('d3adt2dv',2,1), ('d3adtdv2',1,2), ('d3adv3',0,3)]
        self._dn_g_list = ['dgdn', 'd2gdndt', 'd2gdndp', 'd3gdndt2', 
        'd3gdndtdp', 'd3gdndp2', 'd4gdndt3', 'd4gdndt2dp', 'd4gdndtdp2', 
        'd4gdndp3']
        self._d2n_g_list = ['d2gdn2', 'd3gdn2dt', 'd3gdn2dp', 'd4gdn2dt2', 
        'd4gdn2dtdp', 'd4gdn2dp2', 'd5gdn2dt3', 'd5gdn2dt2dp', 'd5gdn2dtdp2', 
        'd5gdn2dp3']
        self._d3n_g_list = ['d3gdn3', 'd4gdn3dt', 'd4gdn3dp', 'd5gdn3dt2', 
        'd5gdn3dtdp', 'd5gdn3dp2', 'd6gdn3dt3', 'd6gdn3dt2dp', 'd6gdn3dtdp2', 
        'd6gdn3dp3']
        self._module = 'untitled'
        self._params = []
        self._variables = []
        self._expression_parts = []
        self._implicit_functions = []
        self._model_type = model_type
        self._T_r = 298.15 # Kelvins
        self._P_r = 1.0    # bars
        self._V_r = 2.5    # J/bar
        self._include_dh_code = False
        assert (model_type in ['TP', 'TV'])
        if model_type == 'TP':
            self._params.append(Parameter('T_r', 'K', sym.symbols('T_r')))
            self._params.append(Parameter('P_r', 'bar', sym.symbols('P_r')))
            self._variables.append(Parameter('T', 'K', sym.symbols('T')))
            self._variables.append(Parameter('P', 'bar', sym.symbols('P')))
        elif model_type == 'TV':
            self._params.append(Parameter('T_r', 'K', sym.symbols('T_r')))
            self._params.append(Parameter('V_r', 'J/bar', sym.symbols('V_r')))
            self._variables.append(Parameter('T', 'K', sym.symbols('T')))
            self._variables.append(Parameter('V', 'J/bar', sym.symbols('V')))
        self.calc_h = None
        self.calib_h = None
        self.fast_h_template = None
        self.calib_h_template = None
        self.fast_c_template = None 
        self.calib_c_template = None
        self._formula_string = None
        self._conversion_string = None
        self._test_string = None

        component_string = ''
        endmember_dict = { "Agamma":"Agamma", "Bgamma":"Bgamma",
        "AsubG":"AsubG", "AsubH":"AsubH", "AsubV":"AsubV", "AsubJ":"AsubJ",
        "AsubKappa":"AsubKappa", "AsubEx":"AsubEx",
        "BsubG":"BsubG", "BsubH":"BsubH", "BsubV":"BsubV", "BsubJ":"BsubJ",
        "BsubKappa":"BsubKappa", "BsubEx":"BsubEx" }
        ss_list = []
        T = self.get_symbol_for_t()
        P = self.get_symbol_for_p()
        for i in range(1,nc+1):
            component_string += 'n' + str(i) + ' '
            ss_string = 'mu' + str(i)
            ss_list.append(sym.Function(ss_string)(T,P))
            endmember_dict[ss_string] = '(*endmember[' + str(i-1) + '].mu0)'
        self._n = sym.Matrix(list(sym.symbols(component_string)))
        self._nT = (sym.ones(1,nc) * self._n)[0]
        self._mu = sym.Matrix(ss_list)
        
        self._printer = SubCodePrinter(
            settings={'user_functions':endmember_dict})

        self._2d_symm_index = []
        for i in range(1,nc+1):
            for j in range (i,nc+1):
                self._2d_symm_index.append((i,j))

        self._3d_symm_index = []
        for i in range (1,nc+1):
            for j in range (i,nc+1):
                for k in range (j,nc+1):
                    self._3d_symm_index.append((i,j,k))
        
        return

    @property
    def nc(self):
        """
        Number of thermodynamic components in the solution

        Returns
        -------
        Number of components (int)  
        """ 
        return self._nc

    @property
    def n(self):
        """
        1-d matrix of SymPy symbols for the mole numbers of thermodynamic 
        components in the model
                
        Returns
        -------
        SymPy Matrix object (sympy.Matrix)
        """
        return self._n
    
    @property
    def nT(self):
        """
        SymPy expression for the total number of moles of all thermodynamic
        components in the solution (expressed as elements of n)
                
        Returns
        -------
        SymPy symbol object (sympy.Symbol)
        """
        return self._nT
    
    @property
    def mu(self):
        """
        1-d matrix of SymPy symbols for the chemical potentials of thermodynamic 
        components in the model

        Returns
        -------
        SymPy Matrix object (sympy.Matrix)
        """
        return self._mu 

    @property
    def a_list(self):
        """
        Array of tuples used internally by the class
    
        Each tuple contains the following entries:  
    
        - [0] str - Function name for derivatives of the Helmholtz energy
        - [1] int - Order of temperature derivative
        - [2] int - Order of volume derivative      
    
        Generally for internal use by the class

        Returns
        -------
        Array of tuples
        """
        return self._a_list
    
    @a_list.setter
    def a_list(self, a_list):
        self._a_list = a_list

    @property
    def dn_g_list(self):
        """
        List of strings that identify first order compositional derivatives produced 
        by the class

        The create_code_module(...) method generates temperature, pressure and
        compositional derivatives of the Gibbs free energy of solution. This property 
        returns a list of all the first order compositional derivatives callable from 
        the generated module. Some of these derivatives will evaluate to zero if the 
        minimal_deriv_set method parameter is initialized to True.

        Returns
        -------
        List of strings, [str,...]                     
        """
        return self._dn_g_list
    
    @property
    def d2n_g_list(self):
        """
        List of strings that identify second order compositional derivatives produced 
        by the class

        The create_code_module(...) method generates temperature, pressure, and
        compositional derivatives of the Gibbs free energy of solution. This property 
        returns a list of all the second order compositional derivatives callable from 
        the generated module. The majority of these derivatives will evaluate to zero 
        if the minimal_deriv_set method parameter is initialized to True.
        
        Returns
        -------
        List of strings, [str,...]
        """
        return self._d2n_g_list
    
    @property
    def d3n_g_list(self):
        """
        List of strings that identify third order compositional derivatives produced 
        by the class

        The create_code_module(...) method generates temperature, pressure, and
        compositional derivatives of the Gibbs free energy of solution. This property 
        returns a list of all the third order compositional derivatives callable from 
        the generated module. Most of these derivatives will evaluate to zero if the 
        minimal_deriv_set method parameter is initialized to True.
        
        Returns
        -------
        List of strings, [str,...]
        """ 
        return self._d3n_g_list

    @property
    def g_list(self):
        """
        Array of tuples used internally by the class
        
        Each tuple contains the following entries:  
        
        - [0] str - Function name for derivatives of the Gibbs free energy
        - [1] int - Order of temperature derivative     
        - [2] int - Order of pressure derivative    
            
        Generally for internal use by the class
                
        Returns
        -------
        Array of tuples
        """
        return self._g_list
    
    @g_list.setter
    def g_list(self, g_list):
        self._g_list = g_list

    @property
    def expression_parts(self):
        """
        Array of Expression class instances

        Builds a model by layering expression instances. Each instance may be 
        applicable over all or a portion of the domain space of variables.
        
        Returns
        -------
        Array of Expression class instances, [Expression,...]
        """
        return self._expression_parts
    
    @expression_parts.setter
    def expression_parts(self, expression_parts):
        self._expression_parts = expression_parts

    @property
    def implicit_functions(self):
        """
        Array of ImplicitFunction class instances
    
        A model may require that one or more implicit function expressions be 
        satisfied prior to evaluation of the model expressions.
        
        Returns
        -------
        Array of ImplicitFunction class instances, [ImplicitFunction,...]
        """
        return self._implicit_functions
    
    @implicit_functions.setter
    def implicit_functions(self, implicit_functions):
        self._implicit_functions = implicit_functions

    @property
    def module(self):
        """
        Name of module

        Returns
        -------
        Name of module (str)
        """
        return self._module

    @module.setter
    def module(self, module):
        assert isinstance(module, str)
        self._module = module

    @property
    def params(self):
        """
        Array of Parameter class instances

        Parameters are calibrated quantities that define the model, like T_r and P_r.

        Returns
        -------
        Array of Parameter class instances, [Parameter,...] 
        """
        return self._params
    
    @params.setter
    def params(self, params):
        self._params = params

    @property
    def printer(self):
        """
        Instance of a SymPy code printer class used in code generation

        Returns
        -------
        C99CodePrinter
        """
        return self._printer
    
    @printer.setter
    def printer(self, printer):
        self._printer = printer

    @property
    def variables(self):
        """
        Array of Parameter class instances
        
        This property returns or sets an array of variables, either T, P, T_r and P_r 
        (if model_type is ’TP’) or T, V, T_r and V_r  (if model type is ’TV’).  Array 
        elements are encapsulated as instances of the Parameter class object. This is 
        done simply as a convenience wrapper. Note that SymPy symbols that are suitable 
        for developing model expressions using these variables are probably more easily 
        accessed using the methods:
        
        - get_symbol_for_t()
        - get_symbol_for_p()
            
        and so on 
            
        Returns
        -------
        Array of Parameter class instances, [Parameter,...]
        """
        return self._variables
    
    @variables.setter
    def variables(self, variables):
        self._variables = variables

    @property
    def formula_string(self):
        """
        String that describes how the chemical formula of a phase will be displayed

        The formula string consists of a concatenation of descriptors that allow the 
        formula of the phase to be generated from its bulk composition. An example is:
        'Ca[Ca]Na[Na]K[K]Al[Al]Si[Si]O8', suitable for describing the compositions of 
        feldspars in the system NaAlSi3O8 - CaAl2Si2O8 - KAlSi3O8. The terms in square
        brackets within the string must contain standard element symbols. These bracketed
        quantities are replaced by the actual number of moles of the indicated element
        present in one mole of the solution when the formula of the phase is requested 
        in the generated module. Other than the bracketed quantities, the content of the 
        string is arbitrary.
        
        Returns
        -------
        str     
        """
        return self._formula_string

    @formula_string.setter
    def formula_string(self, formula_string):
        self._formula_string = formula_string

    @property
    def conversion_string(self):
        """
        List of strings that describe how to assign mole numbers to phase components

        Each string in the list equates the bracketed component number to a mathematical 
        expression involving mole numbers of elements in the bulk composition of the 
        solution. For example, ['[0]=[Na]', '[1]=[Ca]', '[2]=[K]'], which says that all 
        the sodium in the solution defines the moles of the 1st 
        comment (components are indexed starting with 0), all the calcium the 2nd 
        component, and all the potassium the 3rd. 
        
        The righthand side of each expression in the list may be more complex, 
        for example [Ca]-[Na]/2+[K], if appropriate. The conversion strings are used to 
        convert bulk composition of the solution when expressed as elemental 
        concentrations to concentrations of the endmember components.

        Returns
        -------
        List of strings, [str,...]
        """
        return self._conversion_string
    
    @conversion_string.setter
    def conversion_string(self, conversion_string):
        self._conversion_string = conversion_string

    @property
    def test_string(self):
        """
        List of strings that describe feasible limits on the mole numbers of phase 
        components

        Each string in the list constrains the bracketed component number using a 
        mathematical relation, e.g. ['[0] > 0.0', '[1] > 0.0', '[2] > 0.0'], 
        which says that the first endmember component (Components are indexed starting 
        with 0) must always have a value greater than 0, and similarly for the 2nd and 
        3rd components. The mathematical expressions may be complex, as '[0]+[1] > 0.0', 
        or '[0] - [1] < 1/2', as required. The test strings are used to constrain 
        feasible values of component mole numbers when the bulk composition of the 
        solution is set.

        Returns
        -------
        List of strings, [str,...]
        """
        return self._test_string
    
    @test_string.setter
    def test_string(self, test_string):
        self._test_string = test_string

    def get_symbol_for_t(self):
        """
        Retrieves SymPy symbol for temperature.

        Returns
        -------
        T : sympy.core.symbol.Symbol
            SymPy symbol for temperature
        """
        for x in self.variables:
            if x.name == 'T':
                return x.expression
        return None

    def get_symbol_for_p(self):
        """
        Retrieves SymPy symbol for pressure.

        Returns
        -------
        P : sympy symbol
            SymPy symbol for pressure
        """
        for x in self.variables:
            if x.name == 'P':
                return x.expression
        return None

    def get_symbol_for_v(self):
        """
        Retrieves SymPy symbol for volume.

        Returns
        -------
        V : sympy symbol
            SymPy symbol for volume
        """
        for x in self.variables:
            if x.name == 'V':
                return x.expression
        return None

    def get_symbol_for_tr(self):
        """
        Retrieves SymPy symbol for reference temperature.

        Returns
        -------
        Tr : sympy symbol
             SymPy symbol for reference temperature
        """
        for x in self.params:
            if x.name == 'T_r':
                return x.expression
        return None

    def get_symbol_for_pr(self):
        """
        Retrieves SymPy symbol for reference pressure.

        Returns
        -------
        Pr : sympy symbol
             SymPy symbol for reference pressure
        """
        for x in self.params:
            if x.name == 'P_r':
                return x.expression
        return None

    def get_symbol_for_vr(self):
        """
        Retrieves SymPy symbol for reference volume.

        Returns
        -------
        Vr : sympy symbol
             SymPy symbol for reference volume
        """
        for x in self.params:
            if x.name == 'V_r':
                return x.expression
        return None

    def set_reference_origin(self, Tr=298.15, Pr=1.0, Vr=2.5):
        """
        Sets the reference conditions for the model.

        Tr must be specified along with one of Pr or Vr.

        Parameters
        ----------
        Tr : float
             Reference temperature for the model (in Kelvins)
        Pr : float
             Reference pressure for the model (in bars)
        Vr : float
             Reference volume of the model (in J/bar)
        """
        self._T_r = Tr
        self._P_r = Pr
        self._V_r = Vr

    def symmetric_index_from_2d_array(self, elm=(0,0)):
        """
        Given a 2D array index, retrieves a 1D symmetric storage index.

        Parameters
        ----------
        elm : tuple of int
            Two element tuple containing the array index
            
            Values must be in the range 1 to nc.

        Returns
        -------
        Index of 1D compact array corresponding to the specified element of the
        symmetric array
        """
        assert len(elm) == 2
        assert isinstance(elm[0], int)
        assert isinstance(elm[1], int)
        assert elm[0] in range(1,self.nc+1)
        assert elm[1] in range(1,self.nc+1)
        return self._2d_symm_index.index(tuple(sorted(elm)))

    def symmetric_index_from_3d_array(self, elm=(0,0,0)):
        """
        Given a 3D array index, retrieves a 1D symmetric storage index.

        Parameters
        ----------
        elm : tuple of int
            Three element tuple containing the array index
            
            Values must be in the range 1 to nc.

        Returns
        -------
        Index of 1D compact array corresponding to the specified element of the
        symmetric array
        """
        assert len(elm) == 3
        assert isinstance(elm[0], int)
        assert isinstance(elm[1], int)
        assert isinstance(elm[2], int)
        assert elm[0] in range(1,self.nc+1)
        assert elm[1] in range(1,self.nc+1)
        assert elm[2] in range(1,self.nc+1)
        return self._3d_symm_index.index(tuple(sorted(elm)))

    def get_model_param_names(self):
        """
        Retrieves a list of strings of model parameter names.

        Returns
        -------
        result : list of str
            Names of model parameters
        """
        result = []
        for x in self.params:
            result.append(x.name)
        return result

    def get_model_param_units(self):
        """
        Retrieves a list of strings of model parameter units.

        Returns
        -------
        result : list of str
            Units of model parameters
        """
        result = []
        for x in self.params:
            result.append(x.units)
        return result

    def get_model_param_symbols(self):
        """
        Retrieves a list of SymPy symbols for model parameters.

        Returns
        -------
        result : list of sym.Symbol
            SymPy symbols for model parameters
        """
        result = []
        for x in self.params:
            result.append(x.expression)
        return result

    def get_include_dh_code(self):
        """
        Retrieves boolean flag specifying whether code implementing the Debye-Hückel 
        functions will be included in generated code.

        Returns
        -------
        boolean :
                  True or False (default)
        """
        return self._include_dh_code

    def set_include_dh_code(self,include=False):
        """
        Sets a boolean flag controlling the inclusion code implementing the 
        Debye-Hückel functions.

        Parameters
        ----------
        include : bool
                  True or False
        """
        if include:
            self._include_dh_code = True
        else:
            self._include_dh_code = False

    def add_expression_to_model(self, expression, params, 
        exp_type='unrestricted',
        lower_limits=(None, None), upper_limits=(None, None), 
        implicit_functions=None):
        """
        Adds an expression and associated parameters to the model.

        Adds an expression for the Gibbs or Helmholtz energy to the 
        solution model along with a description of expression parameters.

        Parameters
        ----------
        expression : sympy.core.symbol.Symbol
            A SymPy expression for the Gibbs free energy (if model_type is 'TP')
            or the Helmholtz energy (if model_type is 'TV').
            
            The expression may contain an implicit variable, *f*, whose value
            is a function of T,P or T,V.  The value of *f* is determined 
            numerically by solving the implicit_function expression (defined 
            below) at runtime using Newton's method.

        params : An array of tuples
            Structure (string, string, SymPy expression):
            
            - [0] str - Name of the parameter
            
            - [1] str - Units of parameter
            
            - [2] sympy.core.symbol.Symbol - SymPy symbol for the parameter

        exp_type : str, default='unrestricted'
            - 'unrestricted' - Expression is applicable over the whole of T,P or T,V space.
            
            - 'restricted' - Expression applies only between the specified 
              lower_limits and upper_limits.

        lower_limits : tuple, default (None, None)
            A tuple of SymPy expressions defining the inclusive lower (T,P) or 
            (T,V) limit of applicability of the expression. Used only if 
            exp_type is set to 'restricted'.

        upper_limits : tuple, default (None, None)
            A tuple of SymPy expressions defining the inclusive upper (T,P) or 
            (T,V) limit of applicability of the expression. Used only if 
            exp_type is set to 'restricted'.
            
        implicit_functions : Array of tuples
            A tuple element contains three parts:
              
            - [0] sympy.core.symbol.Symbol - SymPy expression for an implicit 
              function in two independent (T,P or T,V) variables and one 
              dependent variable (*f*). The function is equal to zero.
                 
            - [1] sympy.core.symbol.Symbol - SymPy symbol for the dependent 
              variable *f*.
                 
            - [2] sympy.core.symbol.Symbol - SymPy expression that initializes *f* 
              in the iterative routine. This expression must be defined in 
              terms of known parameters and Tr, Pr, T, P.
        """
        l_params = []
        for (x,y,z) in params:
            l_params.append(Parameter(x,y,z))
        l_exp = Expression(expression, l_params, 
            exp_type=exp_type, 
            lower_limit_T=lower_limits[0], lower_limit_PorV=lower_limits[1],
            upper_limit_T=upper_limits[0], upper_limit_PorV=upper_limits[1])

        for x in l_params:
            self._params.append(x)
        self._expression_parts.append(l_exp)

        if implicit_functions:
            for (x,y,z) in implicit_functions:
                self._implicit_functions.append(Implicit_Function(x,y,z))

        if upper_limits[0] == None and upper_limits[1] == None:
            return
        elif upper_limits[0] != None and upper_limits[1] != None:
            if self._model_type == 'TP':
                T = self.get_symbol_for_t()
                P = self.get_symbol_for_p()
                s = -expression.diff(T).subs(T,upper_limits[0]).subs(
                    P,upper_limits[1])
                v =  expression.diff(P).subs(T,upper_limits[0]).subs(
                    P,upper_limits[1])
                g =  -(T-upper_limits[0])*s + (P-upper_limits[1])*v
                self.expression_parts.append(
                    Expression(g, l_params, exp_type='restricted',
                        lower_limit_T=upper_limits[0], 
                        lower_limit_PorV=upper_limits[1],
                        upper_limit_T=None,
                        upper_limit_PorV=None)
                    )
            elif model_type == 'TV':
                T = self.get_symbol_for_t()
                V = self.get_symbol_for_v()
                s = -expression.diff(T).subs(T,upper_limits[0]).subs(
                    V,upper_limits[1])
                p = -expression.diff(V).subs(T,upper_limits[0]).subs(
                    V,upper_limits[1])
                a =  -(T-upper_limits[0])*s - (V-upper_limits[1])*p
                self.expression_parts.append(
                    Expression(a, l_params, exp_type='restricted', 
                        lower_limit_T=upper_limits[0], 
                        lower_limit_PorV=upper_limits[1], 
                        upper_limit_T=None, 
                        upper_limit_PorV=None)
                    )
        elif upper_limits[0] != None:
            if self._model_type == 'TP':
                T = self.get_symbol_for_t()
                s = -expression.diff(T).subs(T,upper_limits[0])
                g =  -(T-upper_limits[0])*s
                self.expression_parts.append(
                    Expression(g, l_params, exp_type='restricted', 
                        lower_limit_T=upper_limits[0], 
                        lower_limit_PorV=upper_limits[1], 
                        upper_limit_T=None, 
                        upper_limit_PorV=None)
                    )
            elif model_type == 'TV':
                T = self.get_symbol_for_t()
                s = -expression.diff(T).subs(T,upper_limits[0])
                a =  -(T-upper_limits[0])*s
                self.expression_parts.append(
                    Expression(a, l_params, exp_type='restricted', 
                        lower_limit_T=upper_limits[0], 
                        lower_limit_PorV=upper_limits[1], 
                        upper_limit_T=None, 
                        upper_limit_PorV=None)
                    )
        elif upper_limits[1] != None:
            if self._model_type == 'TP':
                P = self.get_symbol_for_p()
                v =  expression.diff(P).subs(P,upper_limits[1])
                g =  (P-upper_limits[1])*v
                self.expression_parts.append(
                    Expression(g, l_params, exp_type='restricted', 
                        lower_limit_T=upper_limits[0], 
                        lower_limit_PorV=upper_limits[1], 
                        upper_limit_T=None, 
                        upper_limit_PorV=None)
                    )
            elif model_type == 'TV':
                V = self.get_symbol_for_v()
                p = -expression.diff(V).subs(V,upper_limits[1])
                a = -(V-upper_limits[1])*p
                self.expression_parts.append(
                    Expression(a, l_params, exp_type='restricted', 
                        lower_limit_T=upper_limits[0], 
                        lower_limit_PorV=upper_limits[1], 
                        upper_limit_T=None, 
                        upper_limit_PorV=None)
                    )

    def create_soln_calc_h_file(self, language='C', module_type='fast',
        minimal_deriv_set=False):
        """
        Creates an include file implementing model calculations.

        Note that this include file contains code for functions that implement 
        the model for the generic case. It is meant to be included into a file 
        that implements a specific parameterized instance of the model. See 
        create_code_module().

        The user does not normally call this function directly.

        Parameters
        ----------
        language : string
            Language syntax for generated code. ("C" is the C99 programming 
            language.)
        module_type : string
            Generate code that executes "fast" but does not expose hooks for 
            model parameter calibration. Alternately, generate code suitable for 
            "calib"ration of parameters in the model, which executes more slowly 
            and exposes additional functions that allow setting of parameters 
            and generation of derivatives of thermodynamic functions with 
            respect to model parameters.
        minimal_deriv_set : bool
            Generate a minimal set of compositional derivatives: dgdn, d2gdndt, 
            d2gdndp, d3gdndt2, d3gdndtdp, d3gdndp2, d4gdndt3, d4gdndt2dp, 
            d4gdndtdp2, d4gdndp3, d2gdn2, d3gdn2dt, d3gdn2dp, d3gdn3.  This is 
            the subset of derivatives currently required for solution phases 
            that are imported into the phases module. Remaining derivatives are
            returned with values of zero. 
        """
        assert language == "C"
        # Code for assigning mole numbers to components
        c = self.nc
        moles_assign_text = ''
        for i in range(0,c):
            moles_assign_text += ('    double n' + str(i+1) + ' = n[' 
                + str(i) + '];')
            if i < c-1:
                moles_assign_text += '\n'

        gen_expr = []
        s = '#include <math.h>\n\n'
        if self._include_dh_code:
            s += tpl.create_code_for_dh_functions(language)
            s += "\n"
        solution_calc_template = tpl.create_soln_calc_template()
        solution_derivative_template = tpl.create_soln_deriv_template()
        n = self.n
        T = self.get_symbol_for_t()
        P = self.get_symbol_for_p()
        for idx, (f, oT, oP) in enumerate(self.g_list):
            #print (idx, f, oT, oP)
            # Potential and T, P derivatives
            a = sym.S.Zero
            for exp in self.expression_parts:
                assert exp.exp_type == 'unrestricted'
                a += exp.expression
            a = sym.diff(a, T, oT, P, oP)
            gen_expr.append(a)
            s += solution_calc_template.format(module=self.module, func=f, 
                g_code=self.printer.doprint(a, assign_to='result'), 
                number_components=c, moles_assign=moles_assign_text)
            # First compositional derivatives
            derivative_code_text = ''
            for j in range(0,c):
                if minimal_deriv_set and idx > 5:
                    derivative_code_text += ('    result[' + str(j) + '] = ' 
                        + '0.0' + ';\n')
                else:
                    #print (idx, f, oT, oP, j)
                    derivative_code_text += self.printer.doprint(a.diff(n[j]), 
                        assign_to='result[' + str(j) + ']') + '\n'
            s += solution_derivative_template.format(module=self.module, 
                func=self.dn_g_list[idx], number_components=c, 
                derivative_code=derivative_code_text, 
                moles_assign=moles_assign_text)
            # Second compositional derivatives
            derivative_code_text = ''
            l = 0
            for j in range(0,c):
                for k in range(j,c):
                    if minimal_deriv_set and idx > 2:
                        derivative_code_text += ('    result[' + str(l) + '] = ' 
                            + '0.0' + ';\n')
                    else:
                        #print (idx, f, oT, oP, j, k)
                        derivative_code_text += self.printer.doprint(
                            a.diff(n[j]).diff(n[k]), 
                            assign_to='result[' + str(l) + ']') + '\n'
                    l += 1
            s += solution_derivative_template.format(module=self.module, 
                func=self.d2n_g_list[idx], number_components=c, 
                derivative_code=derivative_code_text, 
                moles_assign=moles_assign_text)
            # Third compositional derivatives
            derivative_code_text = ''
            m = 0
            for j in range(0,c):
                for k in range(j,c):
                    for l in range(k,c):
                        if minimal_deriv_set and idx > 0:
                            derivative_code_text += ('    result[' + str(m) 
                                + '] = ' + '0.0' + ';\n')
                        else:
                            #print (idx, f, oT, oP, j, k, l)
                            derivative_code_text += self.printer.doprint(
                                    a.diff(n[j]).diff(n[k]).diff(n[l]),
                                    assign_to='result[' + str(m) + ']') + '\n'
                        m += 1
            s += solution_derivative_template.format(module=self.module, 
                func=self.d3n_g_list[idx], number_components=c, 
                derivative_code=derivative_code_text, 
                moles_assign=moles_assign_text)

        self._g_matrix = sym.Matrix(gen_expr)
        convenience_template = tpl.create_soln_redun_template()
        s += convenience_template.format(module=self.module, 
            number_components=c)
        return s

    def create_soln_calib_h_file(self, language='C'):
        """
        Creates an include file implementing model parameter derivatives.

        Note that this include file contains code for functions that implement 
        the model for the generic case. It is meant to be included into a file 
        that implements a specific parameterized instance of the model. See 
        create_code_module().

        The user does not normally call this function directly.

        Parameters
        ----------
        language : string
            Language syntax for generated code. ("C" is the C99 programming 
            language.)
        """
        assert language == "C"
        # Code for assigning mole numbers to components
        c = self.nc
        moles_assign_text = ''
        for i in range(0,c):
            moles_assign_text += ('    double n' + str(i+1) + ' = n[' 
                + str(i) + '];')
            if i < c-1:
                moles_assign_text += '\n'

        solution_calib_template = tpl.create_soln_calib_template()
        symparam = self.get_model_param_symbols()
        G_param_jac = sym.Matrix(self._g_matrix).jacobian(symparam)
        s = '#include <math.h>\n\n'

        for j in range(0,len(self.g_list)):
            G_jac_list = [self.printer.doprint(
                G_param_jac[j,i], assign_to='result') for i in range(
                0, len(self.params)) ]
    
            switch_code_text = ''
            for i in range(0, len(self.params)):
                switch_code_text += '    case ' + str(i) + ':\n'
                switch_code_text += '        ' + G_jac_list[i] + '\n'
                switch_code_text += '        break;\n'
            
            s += solution_calib_template.format(module=self.module, 
                number_components=c, func=self.g_list[j][0], 
                switch_code=switch_code_text, moles_assign=moles_assign_text)

        s += ('static void ' + self.module + '_dparam_dgdn(double T, double P, '
            +'double n[' + str(c) + '], int index, double result[' + str(c) 
            + ']) {\n')
        s += moles_assign_text
        s += '\n'
        switch_code_text = '    switch (index) {\n'
        n = self.n
        for i in range(0, len(self.params)):
            a = sym.S.Zero
            for exp in self.expression_parts:
                a += exp.expression
            a = a.diff(symparam[i])
            switch_code_text += '    case ' + str(i) + ':\n'
            for j in range(0,c):
                switch_code_text += ('        ' + 
                    self.printer.doprint(a.diff(n[j]), 
                        assign_to='result[' + str(j) + ']') + '\n')
            switch_code_text += '        break;\n'
        s += switch_code_text
        s += '    default:\n'
        s += '        break;\n'
        s += '    }\n'
        s += '}\n'

        value_params=[self.printer.doprint(symparam[i]) for i in range(0, 
            len(self.params))]
        code_block_one_text = ''
        code_block_two_text = ''
        code_block_three_text = ''
        code_block_four_text = ''
        for i in range(0,len(value_params)):
            code_block_one_text   += ('    (*values)[' + str(i) + '] = ' 
                + value_params[i] + ';\n')
            code_block_two_text   += ('    ' + value_params[i] + ' = values[' 
                + str(i) + '];\n')
            code_block_three_text += ('    case ' + str(i) + ':\n' 
                + '        result = ' + value_params[i] + ';\n' 
                + '        break;\n')
            code_block_four_text  += ('    case ' + str(i) + ':\n' 
                + '        ' + value_params[i] + ' = value;\n' 
                + '        break;\n')

        name_params = self.get_model_param_names()
        unit_params = self.get_model_param_units()
        extra_template = tpl.create_soln_calib_extra_template()
        s += extra_template.format(module=self.module, 
            number_params=len(self.params), 
            names_params=json.dumps(name_params).replace(
                '[', '{').replace(']', '}'), 
            units_params=json.dumps(unit_params).replace(
                '[', '{').replace(']', '}'), 
            code_block_one=code_block_one_text, 
            code_block_two=code_block_two_text, 
            code_block_three=code_block_three_text, 
            code_block_four=code_block_four_text)

        return s

    def create_formula_code_block(self):
        """
        Creates a block of code that computes a formula for this phase.

        The formula is formatted as specified in the property formula_string.

        The code block has access to the variables: T, P, and n, where
        n is a numpy array containing mole numbers of endmember components.

        Returns
        -------
        text : str
               A block of C code that renders the formula

        Notes
        -----
        The formula string looks like: Ca[Ca]Na[Na]K[K]Al[Al]Si[Si]O8

        """
        if self.formula_string is not None:
            elm = list(core.chem.PERIODIC_ORDER)
            format_str = ''
            a = self._formula_string.split(']')
            ne = 0
            for index, b in enumerate(a):
                a[index] = b.split('[')
                format_str += a[index][0]
                if len(a[index]) > 1:
                    a[index][1] = elm.index(a[index][1])
                    format_str += '%5.3f'
                    ne += 1
            text  = '    double sum, elm[' + str(ne) + '];\n'
            for i in range(0,self._nc):
                text += '    const double *end' + str(i) + ' = '
                text += '(*endmember[' + str(i) + '].elements)();\n'
            text += '    int i;\n'
            text += '    const char *fmt = "' + format_str + '";\n'
            text += '    char *result = (char *) malloc(' 
            text += str(len(format_str)+1) + '*sizeof(char));\n'
            text += '    for (i=0, sum=0.0; i<nc; i++) sum += n[i];\n'
            text += '    if (sum == 0.0) return result;\n'
            ne = 0
            format_out = ''
            format_dlm = ''
            for b in a:
                if len(b) > 1:
                    format_out += format_dlm + 'elm[' + str(ne) + ']'
                    format_dlm = ', '
                    j = b[1]
                    text += '    elm[' + str(ne) + '] = '
                    delim = ''
                    for i in range(0,self._nc):
                        text += delim + 'end' + str(i) + '[' + str(j) + ']'
                        text += '*n[' + str(i) + ']/sum'
                        delim = ' + '
                    text += ';\n'
                    ne += 1
            text += '    sprintf(result, fmt, ' + format_out + ');\n'
            text += '    return result;\n'
        else:
            text  = '    char *result = (char *)malloc(sizeof(char));\n'
            text += '    return result;\n'
        return text

    def create_conversion_code_block(self):
        """
        Creates a block of code that computes moles of endmember components
        from moles of elements according to a specified recipe.

        The recipe is read from the property conversion_string.

        The code block has access to the variables: T, P, and e, where
        e is a numpy array containing mole numbers of elements.

        Returns
        -------
        text : str
               A block of C code that renders the conversion

        Notes
        -----
        The recipe looks like: ['[0]=[Na]', '[1]=[Ca]', '[2]=[K]'],
        a list whose elements describe how to assign moles to each endmember.
        The entries in the list may contain algebraic combinations of elements 
        [Na] and moles of endmember components [0]. 

        """
        if self.conversion_string is not None:
            elm = list(core.chem.PERIODIC_ORDER)
            text  = '    double *n = (double *) malloc(' + str(self._nc)
            text += '*sizeof(double));\n'
            for entry in self._conversion_string:
                text += '    '
                b = entry.split(']')
                for c in b:
                    d = c.split('[')
                    if len(d) > 1:
                        text += d[0]
                        if d[1].isnumeric():
                            text += 'n[' + d[1] + ']'
                        else:
                            text += 'e[' + str(elm.index(d[1])) + ']'
                    else:
                        text += d[0]
                text += ';\n'
        else:
            text  = '    double *n = (double *) calloc(' + str(self._nc)
            text += '*sizeof(double));\n'
        text += '    return n;\n'
        return text

    def create_test_code_block(self):
        """
        Creates a block of code that tests moles of endmember components
        according to a specified recipe to ensure that values are viable.

        The recipe is read from the property test_string.

        The code block has access to the variables: T, P, and n, where
        n is a numpy array containing moles of endmember components.

        Returns
        -------
        text : str
               A block of C code that renders the testing

        Notes
        -----
        The recipe looks like: ['[0] > 0.0', '[1] > 0.0', '[2] > 0.0'],
        a list whose elements describe how to assign moles to each endmember.
        The entries in the list may contain algebraic combinations of moles of 
        endmember components [0], and the length of the list (i.e., the 
        number of tests) is unbounded. 

        """
        text  = '    int result = 1;\n'
        if self.test_string is not None:
            for entry in self._test_string:
                text += '    result &= ('
                b = entry.split(']')
                for c in b:
                    d = c.split('[')
                    if len(d) > 1:
                        text += d[0]
                        if d[1].isnumeric():
                            text += 'n[' + d[1] + ']'
                    else:
                        text += d[0]
                text += ');\n'
        text += '    return result;\n'
        return text

    def create_code_module(self, phase="Feldspar", params={}, endmembers=[],
        identifier=None, prefix="cy", module_type="fast", silent=False, 
        language='C', minimal_deriv_set=False):
        """
        Creates include and code file for a model instance.

        Parameters
        ----------
        phase : str
            Model instance title (e.g., phase name).  Used to name the generated 
            function. Cannot contain blank spaces or special characters; underscore 
            ("_") is permitted. Convention capitalizes the first letter and letter 
            following an underscore ("_") character.
        params : dict
            Parameter values for the model instance.
            The keys of this dictionary are validated against parameter symbols 
            stored for the model.
        endmembers : list of str
            A list of prefixes for standard state property functions for the 
            endmember components of this solution. E.g., "Albite_berman" will
            be used to call functions with names like Albite_berman_g(...)
            If the standard state property routines are coded by the 
            StdStateModel Class in this module, all required functions will be
            generated and they will automatically be compliant with this naming 
            convention.
        identifier : str
            A unique identifier for the model instance.
            Defaults to local when module is created (rounded to the second).
        prefix : str
            Prefix to function names for Python bindings, e.g., 
            {prefix}_{phase}_{module}_g(T,P)
        module_type : str
            Generate code that executes "fast", but does not expose hooks for 
            model parameter calibration. Alternately, generate code suitable for 
            "calib"ration of parameters in the model, which executes more slowly 
            and exposes additional functions that allow setting of parameters 
            and generation of derivatives of thermodynamic functions with 
            respect to model parameters. 
        silent : bool
            Do not print status messages.
        language : str
            Language syntax for generated code. ("C" is the C99 programming 
            language.)
        minimal_deriv_set : bool
            Generate a minimal set of compositional derivatives: dgdn, d2gdndt, 
            d2gdndp, d3gdndt2, d3gdndtdp, d3gdndp2, d4gdndt3, d4gdndt2dp, 
            d4gdndtdp2, d4gdndp3, d2gdn2, d3gdn2dt, d3gdn2dp, d3gdn3.  This is 
            the subset of derivatives currently required for solution phases 
            that are imported into the phases module. Remaining derivatives are
            returned with values of zero.

        Returns
        -------
        result : Boolean
                 True if module is succesfully generated, False if some error occurred 
        """
        assert language == "C"
        assert isinstance(phase, str)
        if not re.match("^[a-zA-Z0-9_]*$", phase):
            print ("Error: ", phase, 
                " is only allowed to have characters a-z, A-Z, 0-9, and _")
            return False
        assert len(params) > 0
        okay = True
        for x in self.params:
            if not x.name in params.keys():
                print ("Error: paramter key ", x.name, 
                    " is missing from passed model parameter dictionary")
                okay = False
        if not okay:
            print ("Error in params dictionary (see above)")
            return False
        module = self.module
        if module == 'untitled':
            print ("Error: Please set module name before called this function")
            return False
        assert len(endmembers) == self.nc
        for x in endmembers:
            assert isinstance(x, str)
        if identifier is None:
            identifier = time.asctime(time.localtime(time.time()))
        
        if not silent:
            print ("Creating generic fast model code file string")
        calc_h_file = self.create_soln_calc_h_file(language=language, 
            module_type=module_type, minimal_deriv_set=minimal_deriv_set)
        if not silent:
            print ("Writing include file to working directory ...")
        with open(module+'_calc.h', 'w') as f:
            f.write(calc_h_file)
        
        if module_type == 'fast':
            if not silent:
                print ("Creating (once only) generic model fast code" 
                    + "template include file string")
            fast_h_template = tpl.create_soln_fast_include_template(
                language=language)
            if not silent:
                print ("Creating (once only) generic model fast code" 
                    + "template code file string")
            fast_c_template = tpl.create_soln_fast_code_template(
                language=language)
            ss_xx = ''
        elif module_type == 'calib':
            if not silent:
                print ("Creating (once only) generic model calib code" 
                    + "template include file string")
            calib_h_template = tpl.create_soln_calib_include_template(
                language=language)
            if not silent:
                print ("Creating (once only) generic model calib code" 
                    + "template code file string")
            calib_c_template = tpl.create_soln_calib_code_template(
                language=language)
            if not silent:
                print ("Creating generic calib model code file string")
            calib_h_file = self.create_soln_calib_h_file(language=language)
            if not silent:
                print ("Writing include file to working directory ...")
            with open(module+'_calib.h', 'w') as f:
                f.write(calib_h_file)
            ss_xx = '_calib'
        else:
            print ("Error: module_type must be set to either 'fast' or 'calib'")
            return False
        
        if not silent:
            print ("Creating code blocks for standard state properties.")
        code_block_three_text = ''
        code_block_four_text = ''
        for i in range(0,len(endmembers)):
            code_block_three_text += '  {\n'
            code_block_three_text += '    ' + endmembers[i] + ss_xx + '_name,\n'
            code_block_three_text += '    ' + endmembers[i] + ss_xx + '_formula,\n'
            code_block_three_text += '    ' + endmembers[i] + ss_xx + '_mw,\n'
            code_block_three_text += '    ' + endmembers[i] + ss_xx + '_elements,\n'
            code_block_three_text += '    ' + endmembers[i] + ss_xx + '_g,\n'
            code_block_three_text += '    ' + endmembers[i] + ss_xx + '_dgdt,\n'
            code_block_three_text += '    ' + endmembers[i] + ss_xx + '_dgdp,\n'
            code_block_three_text += '    ' + endmembers[i] + ss_xx + '_d2gdt2,\n'
            code_block_three_text += '    ' + endmembers[i] + ss_xx + '_d2gdtdp,\n'
            code_block_three_text += '    ' + endmembers[i] + ss_xx + '_d2gdp2,\n'
            code_block_three_text += '    ' + endmembers[i] + ss_xx + '_d3gdt3,\n'
            code_block_three_text += '    ' + endmembers[i] + ss_xx + '_d3gdt2dp,\n'
            code_block_three_text += '    ' + endmembers[i] + ss_xx + '_d3gdtdp2,\n'
            code_block_three_text += '    ' + endmembers[i] + ss_xx + '_d3gdp3\n'
            code_block_three_text += '  },\n'
            if module_type == 'fast':
                code_block_four_text += ('#include "' + endmembers[i] 
                    + '_calc.h"\n')
            elif module_type == 'calib':
                code_block_four_text += ('#include "' + endmembers[i] 
                    + '_calib.h"\n')
        code_block_five_text  = self.create_formula_code_block()
        code_block_six_text = self.create_conversion_code_block()
        code_block_seven_text = self.create_test_code_block()

        code_block_one_text = ''
        code_block_two_text = ''
        sym_params = self.get_model_param_symbols()
        name_params = self.get_model_param_names()
        for i in range(0, len(self.params)):
            symbol = self.printer.doprint(sym_params[i])
            value = str(params[name_params[i]])
            code_block_one_text += ('static const double ' + symbol + ' = ' 
                + value + ';\n')
            code_block_two_text += ('static double ' + symbol + ' = ' 
                + value + ';\n')

        c_file = ''
        h_file = ''
        symparam = self.get_model_param_symbols()
        if module_type == 'fast':
            if not silent:
                print ("Creating fast code and include files ...")
            std_state_h_template = tpl.create_soln_std_state_include_template(
                language=language)
            c_file += std_state_h_template.format(
                code_block_three=code_block_three_text, 
                code_block_four=code_block_four_text)
            c_file += fast_c_template.format(module=module, 
                phase=phase, 
                param_names=[ self.printer.doprint(symparam[i]) for i in range(
                    0, len(params)) ], param_values=params, 
                git_identifier=identifier, 
                number_components=self.nc, 
                code_block_one=code_block_one_text,
                code_block_five=code_block_five_text,
                code_block_six=code_block_six_text,
                code_block_seven=code_block_seven_text)
            h_file += fast_h_template.format(module=module, 
                phase=phase, number_components=self.nc)
            h_file_name = phase + '_' + module + '_calc.h'
            if not silent:
                print ("Writing include file to working directory ...")
            with open(h_file_name, 'w') as f:
                f.write(h_file)
            c_file_name = phase + '_' + module + '_calc.c'
            if not silent:
                print ("Writing code file to working directory ...")
            with open(c_file_name, 'w') as f:
                f.write(c_file)
        elif module_type == 'calib':
            if not silent:
                print ("Creating calib code and include files ...")
            std_state_h_template = tpl.create_soln_std_state_include_template(
                language=language)
            c_file += std_state_h_template.format(
                code_block_three=code_block_three_text, 
                code_block_four=code_block_four_text)
            c_file += calib_c_template.format(
                module=module, 
                phase=phase, 
                param_names=[ self.printer.doprint(symparam[i]) for i in range(
                    0, len(params)) ], param_values=params, 
                git_identifier=identifier, 
                number_components=self.nc, 
                code_block_two=code_block_two_text,
                code_block_five=code_block_five_text,
                code_block_six=code_block_six_text,
                code_block_seven=code_block_seven_text)
            h_file += calib_h_template.format(module=module, 
                phase=phase, number_components=self.nc)
            h_file_name = phase + '_' + module + '_calib.h'
            if not silent:
                print ("Writing include file to working directory ...")
            with open(h_file_name, 'w') as f:
                f.write(h_file)
            c_file_name = phase + '_' + module + '_calib.c'
            if not silent:
                print ("Writing code file to working directory ...")
            with open(c_file_name, 'w') as f:
                f.write(c_file)
        if not silent:
            print ("... done")

        if not silent:
            print ("Writing pyxbld file to working directory ...")
        pyxbld_template = tpl.create_soln_pyxbld_template()
        pyx_file_list = "'" + c_file_name + "'"
        for x in endmembers:
            pyx_file_list += ','
            if module_type == 'fast':
                pyx_file_list += "'" + x + "_calc.c'"
            elif module_type == 'calib':
                pyx_file_list += "'" + x + "_calib.c'"
        pyxbld_file = pyxbld_template.format(files_to_compile=pyx_file_list)
        pyxbld_file_name = module + ".pyxbld"
        with open(pyxbld_file_name, 'w') as f:
            f.write(pyxbld_file)

        if not silent:
            print ("writing pyx file to working directory ...")
        pyx_template = ""
        if module_type == "fast":
            pyx_template = tpl.create_soln_fast_pyx_template()
        elif module_type == "calib":
            pyx_template = tpl.create_soln_calib_pyx_template()
        pyx_file = pyx_template.format(prefix=prefix, phase=phase, 
            module=module, number_components=self.nc, number_species=self.nc)
        pyx_file_name = module + ".pyx"
        with open(pyx_file_name, 'w') as f:
            f.write(pyx_file)

        if not silent:
            print ("Compiling code and Python bindings ...")
        pyximport.install(language_level=3)

        if not silent:
            print ("Success! Import the module named ", module)
        return True

class Ordering_Functions:
    """
    Class holds properties of implicit functions that describe homogeneous 
    equilibria via cation ordering or speciation

    Parameters
    ----------
    ordering_functions : list[sympy.core.symbol.Symbol]
        A list of SymPy expressions defining a system of implicit multi-variable 
        functions having 2+nc independent variables (T,P, and moles of endmember 
        components) and one or more dependent variables (ordering parameters). 
        The number of implicit equations in the system must equal the number of 
        ordering parameters. The implicit function expressions must all be equal 
        to zero; e.g., :math:`log(f) - f = 0`. The list of expressions must be 
        defined in terms of known parameters and Tr, Pr, T, P, n.
    dep_variables : list[sympy.core.symbol.Symbol]
        A list of SymPy symbols for the dependent variables, i.e., the ordering
        parameters or species mole fractions
    initial_guesses : list[sympy.core.symbol.Symbol]
        SymPy expressions that initialize the ordering parameters in the 
        iterative routines that solve the system. These expression must be 
        defined in terms of known solution parameters and Tr, Pr, T, P, n.
    bounds : sympy.logic.boolalg.And
        An instance of the And class, as output from reduce_inequalities(), 
        which holds a logical expression that embodies bound constraints on 
        the ordering parameters
    nullValues : list[sympy.core.symbol.Symbol]
        SymPy expressions that declare null values of the ordering parameters.
        These values wil not be altered during the homogeneous equilibrium 
        calculation.

    Attributes
    ----------
    bounds
    function
    guess
    nullValues
    variable

    Notes
    -----
    This system of equations is solved for the ordering parameters prior to 
    evaluating model expressions for the Gibbs free energy and its derivatives.  
    The equations must be differentiable with respect to model parameters T and 
    P.  
    """
    def __init__(self, ordering_functions, dep_variables, initial_guesses,
        dep_bounds, nullValues):
        self._function = ordering_functions
        self._variable = dep_variables
        self._guess = initial_guesses
        self._bounds = dep_bounds
        self._nullValues = nullValues

    @property
    def function(self):
        """
        Storage for the ordering_functions parameter

        Returns
        -------
        sympy.core.symbol.Symbol
        """
        return self._function
    
    @property
    def variable(self):
        """
        Storage for the dep_variables parameter    

        Returns
        -------
        sympy.core.symbol.Symbol
        """
        return self._variable
    
    @property
    def guess(self):
        """
        Storage for the initial_guesses parameter

        Returns
        -------
        sympy.core.symbol.Symbol
        """
        return self._guess

    @property
    def bounds(self):
        """
        Storage for the bounds parameter

        Returns
        -------
        sympy.logic.boolalg.And
        """
        return self._bounds

    @property
    def nullValues(self):
        """
        Null values of the order parameters
        
        Returns:
        sympy.core.symbol.Symbol
        """
        return self._nullValues
    

class ComplexSolnModel(SimpleSolnModel):
    """
    Class creates a model of the thermodynamic properties of a complex solution.

    Inherits all methods and functionality of the Simple Solution Model class.

    Parameters
    ----------
    nc : int
        Number of thermodynamic components in the solution
    ns : int    
        Number of ordering parameters in the solution
    nw : int    
        Number of species in the solution
    model_type : str
        Model type of 'TP' implies that expressions are of the form of the Gibbs 
        free energy.
        
        Model type of 'TV' implies that expressions are of the form of the 
        Helmholtz energy. This option is infrequently used.

    Attributes
    ----------
    ns
    nw
    order_sub
    s
    verbose
    de_minimis
        
    Methods
    -------
    add_expression_to_model
    create_code_module
    create_ordering_code
    create_soln_calc_h_file
    create_soln_calib_h_file
    eval_asymmetric_regular_param
    eval_endmember
    eval_regular_param
    eval_ternary_param
    exchange_indices
    taylor_expansion
    
    Notes
    -----
    This class is for creating a representation of the thermodynamic properties 
    of a complex solution phase.  The principal use of the class is to construct 
    a model symbolically and to generate code that implements the model for 
    model parameter calibration and thermodynamic property calculation.

    A complex solution is one that does not contain implicit variables that need
    be determined via solution of conditions of homogeneous equilibrium. 
    Examples of complex solutions include aqueous or liquid solutions involving 
    complex formation (speciation) or solid solutions that include 
    cation-ordering or composition dependent symmetry-breaking phase transitions.

    """
    def __init__(self, nc=None, ns=None, nw=None, model_type='TP'):
        super().__init__(nc=nc, model_type=model_type)
        nw = nc if nw is None else nw
        assert nw >= nc, 'Parameter nw must be equal to or exceed nc.'
        self._ns = ns
        self._nw = nw
        self._ordering_functions = None
        self._verbose = False
        self._de_minimis = False

        if ns is not None:
            component_string = ''
            if ns > 1:
                for i in range(1,ns+1):
                    component_string += 's' + str(i) + ' '
                self._s = sym.Matrix(list(sym.symbols(component_string)))
            else:
                self._s = sym.Matrix([sym.symbols('s1')])

            self._order_sub = []
            T = self.get_symbol_for_t()
            P = self.get_symbol_for_p()
            for i in range(1,ns+1):
                v = self._s[i-1]
                f = sym.Function('S'+str(i))(T, P, self.n)
                self._order_sub.append([v,f])

        else:
            self._s = None
            self._order_sub = None

    @property
    def ns(self):
        """
        Number of ordering parameters in the solution 

        Returns
        -------
        Number of ordering parameters in the solution (int)
        """
        return self._ns

    @property
    def s(self):
        """
        1-d matrix of SymPy symbols for the ordering parameters in the model
        
        Returns
        -------
        SymPy Matrix object (sympy.Matrix)
        """
        return self._s

    @property
    def order_sub(self):
        """
        Deprecated
        """
        return self._order_sub

    @property
    def nw(self):
        """
        Number of species in the solution

        Returns
        -------
        Number of species in the solution (int)
        """
        return self._nw

    @property
    def verbose(self):
        """
        Prints verbose progress during execution of *create_soln_calc_h_file*
        and *create_soln_calib_h_file* methods
        
        Initialized to False.

        Returns
        -------
        Verbosity (boolean)
        """
        return self._verbose

    @verbose.setter
    def verbose(self, verbose):
        self._verbose = verbose

    @property
    def de_minimis(self):
        """
        Construct the minimal number of thermodynamic quanitities 
        
        Generaly set for fastest compilation and for calibration.  Computes 
        only g, dgdt, dgdp, dgdn, d2gdndt, d2gdndp, d2gdn2, dgdw, and d2gdndw.
        Default is False.
        
        Returns
        -------
        Boolean
        """
        return self._de_minimis

    @de_minimis.setter
    def de_minimis(self, de_minimis):
        self._de_minimis = de_minimis

    def taylor_expansion(self, order=2, mod_type='non-convergent'):
        """
        Contructs a Taylor expansion of the non-configurational Gibbs free 
        energy of specified order using the SymPy symbols in the list var.

        Parameters
        ----------
        order : int
            Order of the Taylor expansion. Second order is the default.
            The maximum order is 4.
        mod_type : str
            1) 'convergent': order parameter expansion only has symmetric terms
            2) 'non-convergent': order parameter expansion has all terms
        
        Returns
        -------
        output : tuple
            1) Number of terms in the expansion
            2) Taylor expansion as a SymPy expression
            3) List of SymPy symbols for energetic terms in the expansion
            4) List of SymPy expressions multiplying energetic terms in the
               expansion

        """
        nc = self.nc
        ns = self.ns
        assert order >= 1, 'Minimum order of Taylor expansion is 1.' 
        assert order <= 4, 'Maximum order of Taylor expansion is 4.'
        assert nc >= 2,'There must be at least two components in the solution.'
        assert ns >= 1,'There must be at least one ordering parameter '+ \
            'in this solution.'
        assert mod_type == 'convergent' or mod_type == 'non-convergent', \
            'mod_type must be one of "convergent" or "non-convergent".'
        asymmetric = True if mod_type == 'non-convergent' else False

        n = self.n
        s = self.s
        nTot = self.nT
        var = []
        for i in range(1,nc):
            var.append(n[i]/nTot) #form r variables; nc-1 of these
        for i in range(0,ns):
            var.append(s[i])

        taylor_coeff = [sym.symbols('G0')]
        taylor_terms = [sym.S.One]
        taylor = taylor_coeff[sym.S.Zero]
        count = 1
        for i in range(1,nc):
            taylor_coeff.append(sym.symbols('Gr'+str(i)))
            taylor += taylor_coeff[count]*var[i-1]
            taylor_terms.append(var[i-1])
            count += 1
        if asymmetric:
            for i in range(1,ns+1):
                taylor_coeff.append(sym.symbols('Gs'+str(i)))
                taylor += taylor_coeff[count]*var[nc-2+i]
                taylor_terms.append(var[nc-2+i])
                count += 1
        if order > 1:
            for i in range(1,nc):
                for j in range(i,nc):
                    taylor_coeff.append(sym.symbols('Grr'+str(i)+str(j)))
                    taylor += taylor_coeff[count]*var[i-1]*var[j-1]
                    taylor_terms.append(var[i-1]*var[j-1])
                    count += 1
                if asymmetric:
                    for j in range(1,ns+1):
                        taylor_coeff.append(sym.symbols('Grs'+str(i)+str(j)))
                        taylor += taylor_coeff[count]*var[i-1]*var[nc-2+j]
                        taylor_terms.append(var[i-1]*var[nc-2+j])
                        count += 1
            for i in range(1,ns+1):
                for j in range(i,ns+1):
                    taylor_coeff.append(sym.symbols('Gss'+str(i)+str(j)))
                    taylor += taylor_coeff[count]*var[nc-2+i]*var[nc-2+j]
                    taylor_terms.append(var[nc-2+i]*var[nc-2+j])
                    count += 1
            if order > 2:
                for i in range(1,nc):
                    for j in range(i,nc):
                        for k in range(j,nc):
                            taylor_coeff.append(sym.symbols(
                                'Grrr'+str(i)+str(j)+str(k)))
                            taylor += taylor_coeff[count]* \
                                      var[i-1]*var[j-1]*var[k-1]
                            taylor_terms.append(var[i-1]*var[j-1]*var[k-1])
                            count += 1
                        if asymmetric:
                            for k in range(1,ns+1):
                                taylor_coeff.append(sym.symbols(
                                    'Grrs'+str(i)+str(j)+str(k)))
                                taylor += taylor_coeff[count]* \
                                          var[i-1]*var[j-1]*var[nc-2+k]
                                taylor_terms.append(var[i-1]*var[j-1]*var[nc-2+k])
                                count += 1
                    for j in range(1,ns+1):
                        for k in range(j,ns+1):
                            taylor_coeff.append(sym.symbols(
                                'Grss'+str(i)+str(j)+str(k)))
                            taylor += taylor_coeff[count]* \
                                      var[i-1]*var[nc-2+j]*var[nc-2+k]
                            taylor_terms.append(
                                var[i-1]*var[nc-2+j]*var[nc-2+k])
                            count += 1
                if asymmetric:
                    for i in range(1,ns+1):
                        for j in range(i,ns+1):
                            for k in range(j,ns+1):
                                taylor_coeff.append(sym.symbols(
                                    'Gsss'+str(i)+str(j)+str(k)))
                                taylor += taylor_coeff[count]* \
                                          var[nc-2+i]*var[nc-2+j]*var[nc-2+k]
                                taylor_terms.append(
                                    var[nc-2+i]*var[nc-2+j]*var[nc-2+k])
                                count += 1
                if order > 3:
                    for i in range(1,nc):
                        for j in range(i,nc):
                            for k in range(j,nc):
                                for l in range(k,nc):
                                    taylor_coeff.append(sym.symbols(
                                        'Grrrr'+str(i)+str(j)+str(k)+str(l)))
                                    taylor += taylor_coeff[count]* \
                                      var[i-1]*var[j-1]*var[k-1]*var[l-1]
                                    taylor_terms.append(
                                        var[i-1]*var[j-1]*var[k-1]*var[l-1])
                                    count += 1
                                if asymmetric:
                                    for l in range(1,ns+1):
                                        taylor_coeff.append(sym.symbols(
                                            'Grrrs'+str(i)+str(j)+str(k)+str(l)))
                                        taylor += taylor_coeff[count]* \
                                           var[i-1]*var[j-1]*var[k-1]*var[nc-2+l]
                                        taylor_terms.append(
                                            var[i-1]*var[j-1]*var[k-1]*var[nc-2+l])
                                        count += 1
                            for k in range(1,ns+1):
                                for l in range(k,ns+1):
                                    taylor_coeff.append(sym.symbols(
                                        'Grrss'+str(i)+str(j)+str(k)+str(l)))
                                    taylor += taylor_coeff[count]* \
                                        var[i-1]*var[j-1]*var[nc-2+k]*var[nc-2+l]
                                    taylor_terms.append(
                                        var[i-1]*var[j-1]*var[nc-2+k]*var[nc-2+l])
                                    count += 1
                        if asymmetric:
                            for j in range(1,ns+1):
                                for k in range(j,ns+1):
                                    for l in range(k,ns+1):
                                        taylor_coeff.append(sym.symbols(
                                            'Grsss'+str(i)+str(j)+str(k)+str(l)))
                                        taylor += taylor_coeff[count]* \
                                        var[i-1]*var[nc-2+j]*var[nc-2+k]*var[nc-2+l]
                                        taylor_terms.append(
                                            var[i-1]*var[nc-2+j]*var[nc-2+k]*var[nc-2+l])
                                        count += 1
                    for i in range(1,ns+1):
                        for j in range(i,ns+1):
                            for k in range(j,ns+1):
                                for l in range(k,ns+1):
                                    taylor_coeff.append(sym.symbols(
                                        'Gssss'+str(i)+str(j)+str(k)+str(l)))
                                    taylor += taylor_coeff[count]* \
                                        var[nc-2+i]*var[nc-2+j]*var[nc-2+k]*var[nc-2+l]
                                    taylor_terms.append(
                                        var[nc-2+i]*var[nc-2+j]*var[nc-2+k]*var[nc-2+l])
                                    count += 1
        return (count,taylor,taylor_coeff,taylor_terms)

    def eval_endmember(self, n_val, s_val, taylor):
        """
        Evaluates a Taylor expansion of the molar non-configurational Gibbs free
        energy at the specified composition.

        Parameters
        ----------
        n_val : []
            List of SymPy constants providing composition in terms of solution 
            components
            
            Must have length nc
        s_val : []
            List of SymPy constants providing values of ordering parameters for 
            the specified composition
            
            Must have length ns
        taylor : str
            A SymPy expression for Taylor expansion of the non-configurational
            Gibbs free energy interms of the elements of n and s

        Returns
        -------
        output : str
            A SymPy expression for the evaluated Taylor series
        """
        assert len(n_val) == self.nc, 'n_val must have length nc.'
        assert len(s_val) == self.ns, 's_val must have length ns.'

        nc = self.nc
        ns = self.ns
        n = self.n
        s = self.s

        sub_list = []
        for i in range(0,nc):
            sub_list.append((n[i],n_val[i]))
        for i in range(0,ns):
            sub_list.append((s[i],s_val[i]))
        return taylor.subs(sub_list)

    def eval_regular_param(self, nA_val, sA_val, nB_val, sB_val, 
        taylor):
        """
        Evaluates a Taylor expansion of the molar non-configurational Gibbs free
        energy equivalent to a regular solution parameter for the A-B join.

        Parameters
        ----------
        nA_val : []
            List of SymPy constants providing composition in terms of solution 
            components of endmember A
            
            Must have length nc
        sA_val : []
            List of SymPy constants providing values of ordering parameters for 
            the specified composition of endmember A
            
            Must have length ns
        nB_val : []
            List of SymPy constants providing composition in terms of solution 
            components of endmember B
            
            Must have length nc
        sB_val : []
            List of SymPy constants providing values of ordering parameters for 
            the specified composition of endmember B
            
            Must have length ns
        taylor : str
            A SymPy expression for Taylor expansion of the non-configurational
            Gibbs free energy interms of the elements of n and s

        Returns
        -------
        output : str
            A SymPy expression for the evaluated Taylor series

        Notes
        -----
        join A-B ==>  4 ( G(A/2+B/2) - G(A)/2 - G(B)/2 )

        """
        assert len(nA_val) == self.nc, 'nA_val must have length nc.'
        assert len(nB_val) == self.nc, 'nB_val must have length nc.'
        assert len(sA_val) == self.ns, 'sA_val must have length ns.'
        assert len(sB_val) == self.ns, 'sB_val must have length ns.'

        nc = self.nc
        ns = self.ns
        n = self.n
        s = self.s

        sub_A_list = []
        sub_B_list = []
        sub_AB_list = []
        for i in range(0,nc):
            sub_A_list.append((n[i],nA_val[i]))
            sub_B_list.append((n[i],nB_val[i]))
            sub_AB_list.append((n[i],(nA_val[i]+nB_val[i])/sym.S(2)))
        for i in range(0,ns):
            sub_A_list.append((s[i],sA_val[i]))
            sub_B_list.append((s[i],sB_val[i]))
            sub_AB_list.append((s[i],(sA_val[i]+sB_val[i])/sym.S(2)))
        gA = taylor.subs(sub_A_list)
        gB = taylor.subs(sub_B_list)
        gAB = taylor.subs(sub_AB_list)

        return 4*gAB - 2*gA - 2*gB

    def eval_asymmetric_regular_param(self, nA_val, sA_val, nB_val, sB_val, 
        taylor):
        """
        Evaluates a Taylor expansion of the molar non-configurational Gibbs free
        energy equivalent to an asymmetric regular solution parameter.

        This method yields a quantity DW, related to the asymmetric regular 
        solution parameters, W112 and W122, as W + DW and W - DW:

        dW = (27/2)G(2A/3+B/3) - 12G(A/2+B/2) - 3G(A) + 3G(B)/2

        See method eval_regular_param to obtain W.

        Parameters
        ----------
        nA_val : []
            List of SymPy constants providing composition in terms of solution 
            components of endmember A
            
            Must have length nc
        sA_val : []
            List of SymPy constants providing values of ordering parameters for 
            the specified composition of endmember A
            
            Must have length ns
        nB_val : []
            List of SymPy constants providing composition in terms of solution 
            components of endmember B
            
            Must have length nc
        sB_val : []
            List of SymPy constants providing values of ordering parameters for 
            the specified composition of endmember B
            
            Must have length ns
        taylor : str
            A SymPy expression for Taylor expansion of the non-configurational
            Gibbs free energy interms of the elements of n and s

        Returns
        -------
        output : str
            A SymPy expression for the evaluated Taylor series

        Notes
        -----
        join A-B ==>  (1/3)(2/3)(W + dW(2/3-1/3)) =  G(2A/3+B/3)-2G(A)/3-G(B)/3

        """
        assert len(nA_val) == self.nc, 'nA_val must have length nc.'
        assert len(nB_val) == self.nc, 'nB_val must have length nc.'
        assert len(sA_val) == self.ns, 'sA_val must have length ns.'
        assert len(sB_val) == self.ns, 'sB_val must have length ns.'

        nc = self.nc
        ns = self.ns
        n = self.n
        s = self.s

        sub_A_list = []
        sub_B_list = []
        sub_AB_list = []
        sub_2AB_list = []
        for i in range(0,nc):
            sub_A_list.append((n[i],nA_val[i]))
            sub_B_list.append((n[i],nB_val[i]))
            sub_AB_list.append((n[i],(nA_val[i]+nB_val[i])/sym.S(2)))
            sub_2AB_list.append(
                (n[i],(sym.S(2)*nA_val[i]+nB_val[i])/sym.S(3)))
        for i in range(0,ns):
            sub_A_list.append((s[i],sA_val[i]))
            sub_B_list.append((s[i],sB_val[i]))
            sub_AB_list.append((s[i],(sA_val[i]+sB_val[i])/sym.S(2)))
            sub_2AB_list.append(
                (s[i],(sym.S(2)*sA_val[i]+sB_val[i])/sym.S(3)))
        gA = taylor.subs(sub_A_list)
        gB = taylor.subs(sub_B_list)
        gAB = taylor.subs(sub_AB_list)
        g2AB = taylor.subs(sub_2AB_list)

        return sym.Rational(27,2)*g2AB - 12*gAB - 3*gA + sym.Rational(3,2)*gB

    def eval_ternary_param(self, nA_val, sA_val, nB_val, sB_val, nC_val, sC_val, 
        taylor):
        """
        Evaluates a Taylor expansion of the molar non-configurational Gibbs free energy equivalent to a strict ternary solution parameter.  

        WT = 27G(A/3+B/3+C/3) - 12G(A/2+B/2) - 12G(A/2,C/2) - 12G(B/2,C/2) + 3G(A) + 
        3G(B) + 3G(C)

        Parameters
        ----------
        nA_val : []
            List of SymPy constants providing composition in terms of solution 
            components of endmember A
            
            Must have length nc
        sA_val : []
            List of SymPy constants providing values of ordering parameters for 
            the specified composition of endmember A
            
            Must have length ns
        nB_val : []
            List of SymPy constants providing composition in terms of solution 
            components of endmember B
            
            Must have length nc
        sB_val : []
            List of SymPy constants providing values of ordering parameters for 
            the specified composition of endmember B
            
            Must have length ns
        nC_val : []
            List of SymPy constants providing composition in terms of solution 
            components of endmember C
            
            Must have length nc
        sC_val : []
            List of SymPy constants providing values of ordering parameters for 
            the specified composition of endmember C
            
            Must have length ns
        taylor : str
            A SymPy expression for Taylor expansion of the non-configurational
            Gibbs free energy interms of the elements of n and s

        Returns
        -------
        output : str
            A SymPy expression for the evaluated Taylor series

        Notes
        -----
        A-B-C ==> (1/27) W(A,B,C) + (1/9) W(A,B) + (1/9) W(A,C) + (1/9) W(B,C)
               = G(A/3+B/3+C/3) - G(A)/3 - G(B)/3 - G(C)/3

        """
        assert len(nA_val) == self.nc, 'nA_val must have length nc.'
        assert len(nB_val) == self.nc, 'nB_val must have length nc.'
        assert len(nC_val) == self.nc, 'nC_val must have length nc.'
        assert len(sA_val) == self.ns, 'sA_val must have length ns.'
        assert len(sB_val) == self.ns, 'sB_val must have length ns.'
        assert len(sC_val) == self.ns, 'sC_val must have length ns.'

        nc = self.nc
        ns = self.ns
        n = self.n
        s = self.s

        sub_A_list = []
        sub_B_list = []
        sub_C_list = []
        sub_AB_list = []
        sub_AC_list = []
        sub_BC_list = []
        sub_ABC_list = []
        for i in range(0,nc):
            sub_A_list.append((n[i],nA_val[i]))
            sub_B_list.append((n[i],nB_val[i]))
            sub_C_list.append((n[i],nC_val[i]))
            sub_AB_list.append((n[i],(nA_val[i]+nB_val[i])/sym.S(2)))
            sub_AC_list.append((n[i],(nA_val[i]+nC_val[i])/sym.S(2)))
            sub_BC_list.append((n[i],(nB_val[i]+nC_val[i])/sym.S(2)))
            sub_ABC_list.append(
                (n[i],(nA_val[i]+nB_val[i]+nC_val[i])/sym.S(3)))
        for i in range(0,ns):
            sub_A_list.append((s[i],sA_val[i]))
            sub_B_list.append((s[i],sB_val[i]))
            sub_C_list.append((s[i],sC_val[i]))
            sub_AB_list.append((s[i],(sA_val[i]+sB_val[i])/sym.S(2)))
            sub_AC_list.append((s[i],(sA_val[i]+sC_val[i])/sym.S(2)))
            sub_BC_list.append((s[i],(sB_val[i]+sC_val[i])/sym.S(2)))
            sub_ABC_list.append(
                (s[i],(sA_val[i]+sB_val[i]+sC_val[i])/sym.S(3)))
        gA = taylor.subs(sub_A_list)
        gB = taylor.subs(sub_B_list)
        gC = taylor.subs(sub_C_list)
        gAB = taylor.subs(sub_AB_list)
        gAC = taylor.subs(sub_AC_list)
        gBC = taylor.subs(sub_BC_list)
        gABC = taylor.subs(sub_ABC_list)

        return 27*gABC - 12*gAB - 12*gAC - 12*gBC + 3*gA + 3*gB + 3*gC
    
    def add_expression_to_model(self, expression, params, 
        ordering_functions=None):
        """
        Adds an expression and associated parameters to the model.

        Parameters
        ----------
        expression : sympy.core.symbol.Symbol
            A SymPy expression for the Gibbs free energy (if model_type is 'TP')
            or the Helmholtz energy (if model_type is 'TV').
            
            The expression may contain an implicit variable, *f*, whose value
            is a function of T,P or T,V.  The value of *f* is determined 
            numerically by solving the implicit_function expression (defined 
            below) at runtime using Newton's method.

        params : An array of tuples
            Structure (string, string, SymPy expression):
            
            - [0] str - Name of the parameter
            
            - [1] str - Units of parameter
            
            - [2] sympy.core.symbol.Symbol - SymPy symbol for the parameter

        ordering_functions : tuple
            A tuple element contains four or five parts:
              
            - [0] [sympy.core.symbol.Symbol] - List of SymPy expressions for a
              system of implicit functions in independent variables T,P,n, and one 
              or more dependent variables (ordering parameters). Each of the functions
              in the system must evaluate to zero.
                 
            - [1] [sympy.core.symbol.Symbol] - List of SymPy symbols for the 
              dependent variables in the specified system.
                 
            - [2] [sympy.core.symbol.Symbol] - List of SymPy expressions that 
              initialize the ordering parameters in the iterative routine that
              solve the system of implicit functions. 

            - [3] [sympy.logic.boolalg.And] - An instance of the And class, as 
              output from reduce_inequalities(), that holds a logical expression
              that embodies bound constraints on the ordering parameters

            - [4] [sympy.core.symbol.Symbol] - List of Sympy symbols for null values
              of the ordering parameters 
        """
        
        super(ComplexSolnModel, self).add_expression_to_model(expression, 
            params)

        if ordering_functions:
            if len(ordering_functions) == 4:
                (x,y,z,v) = ordering_functions
                self._ordering_functions = Ordering_Functions(x,y,z,v,None)
            elif len(ordering_functions) == 5:
                (x,y,z,v,w) = ordering_functions
                self._ordering_functions = Ordering_Functions(x,y,z,v,w)
            else:
                assert True, 'Ordering argument must be a tuple of length 4 of 5.'

    def create_ordering_code(self, moles_assign_text, order_assign_text):
        """
        Generates ordering function code for a complex solution model.

        Parameters
        ----------
        moles_assign_text : str
            Text string containing C code for assigning moles of each 
            component (n1, n2, ..., nc) from a vector of moles, n.

        order_assign_text : str
            Text string containing C code for assigning values of each
            ordering variable (s1, s2, ..., ss) from a vector, s

        Returns
        -------
        output : str
            String of code that implements ordering functions and derivatives.

        Notes
        -----
        The user does not normally call this function directly.

        """
        nc = self.nc
        n  = self.n
        ns = self.ns
        s  = self.s
        sFunc = self._ordering_functions
        T = self.get_symbol_for_t()
        P = self.get_symbol_for_p()

        iter_max = 200

        if self.verbose:
            print ("... ... ordering code")
        cb_0 = moles_assign_text + '    double '
        separator = ''
        for i in range(0,ns):
            cb_0 += separator + 's' + str(i+1)
            separator = ', '
        cb_0 += ';\n'
        cb_1 = ''
        cb_2 = ''
        for i in range(0,ns):
            cb_2 += '        s' + str(i+1) + ' = s[' + str(i) + '];\n'
        cb_3 = ''
        cb_4 = ''
        if sFunc.nullValues is not None:
            for i in range(0,ns):
                cb_4 += '        if (sOld[' + str(i) + '] == '
                cb_4 += self.printer.doprint(sFunc.nullValues[i])
                cb_4 += ') deltaS[' + str(i) + '] = 0.0;\n'
        cb_4 += '        do { \n'
        cb_4 += '            steplength /= 2.0;\n'
        
        moles_assign_text += order_assign_text
        separator = ''
        for i in range(0,ns):
            # (fabs(sNew[0]-sOld[0]) > 10.0*DBL_EPSILON) ||
            cb_1 += separator + '(fabs(s[' + str(i) + ']-sOld[' + str(i) + \
                    ']) > 10.0*DBL_EPSILON)'
            separator = ' || '
            # dgds[0] = DGDS0;
            cb_2 += '        dgds[' + str(i) + '] = '
            cb_2 += self.printer.doprint(sFunc.function[i]) + ';\n'
            # invd2gds2[0][0] = D2GDS0S0;
            for j in range(i,ns):
                cb_3 += '        invd2gds2[' + str(i) + '][' + str(j) + '] = '
                a = sym.diff(sFunc.function[i], s[j], 1)
                cb_3 += self.printer.doprint(a) + ';\n'
                if i < j:
                    cb_3 += '        invd2gds2[' + str(j) + '][' + str(i) + \
                            '] = invd2gds2[' + str(i) + '][' + str(j) + '];\n'
            ###########################################
            # bound constraints on ordering parameters#
            ###########################################
            cb_4 += '            s' + str(i+1) + ' = sOld[' + str(i) \
                   + ']+ steplength*deltaS[' + str(i) + '];\n'

        cb_3 += '\n'
        if ns == 1:
            cb_3 += '        invd2gds2[0][0] = 1.0/invd2gds2[0][0];\n'
        else:
            cb_3 += '        gaussj(invd2gds2);\n'
        cb_4 += '        } while(!(' + self.printer.doprint(sFunc.bounds) \
              + ') && (steplength > DBL_EPSILON));\n'
        cb_4 += '        if (steplength <= DBL_EPSILON) iter = 200;\n'

        cb_5  = '    double d2gdnds['+str(nc)+']['+str(ns)+'];\n'
        cb_5 += moles_assign_text
        cb_6  = '    double d2gdsdt['+str(ns)+'];\n'
        cb_6 += moles_assign_text
        cb_7  = '    double d2gdsdp['+str(ns)+'];\n'
        cb_7 += moles_assign_text
        cb_8  = '    double d2gdnds['+str(nc)+']['+str(ns)+'],d3gdn2ds['+ \
                str(nc)+']['+str(nc)+']['+str(ns)+'],d3gdnds2['+str(nc)+']['+ \
                str(ns)+']['+str(ns)+'],d3gds3['+str(ns)+']['+str(ns)+']['+ \
                str(ns)+'];\n'
        cb_8 += moles_assign_text
        cb_9  = '    double d2gdnds['+str(nc)+']['+str(ns)+'],d3gdndsdt['+ \
                str(nc)+']['+str(ns)+'], d3gdnds2['+str(nc)+']['+str(ns)+']['+ \
                str(ns)+'],d2gdsdt['+str(ns)+'],d3gds2dt['+str(ns)+']['+ \
                str(ns)+'],d3gds3['+str(ns)+']['+str(ns)+']['+str(ns)+'];\n'
        cb_9 += moles_assign_text

        cb_10 = '    double d2gdnds['+str(nc)+']['+str(ns)+'],d3gdndsdp['+ \
                str(nc)+']['+str(ns)+'],d3gdnds2['+str(nc)+']['+str(ns)+']['+ \
                str(ns)+'],d2gdsdp['+str(ns)+'],d3gds3['+str(ns)+']['+str(ns)+ \
                ']['+str(ns)+'],d3gds2dp['+str(ns)+']['+str(ns)+'];\n'
        cb_10 += moles_assign_text
        cb_11  = '    double d2gdsdt['+str(ns)+'],d3gdsdt2['+str(ns)+ \
                 '],d3gds2dt['+ str(ns)+']['+str(ns)+'],d3gds3['+str(ns)+ \
                 ']['+str(ns)+']['+ str(ns)+'];\n'
        cb_11 += moles_assign_text
        cb_12  = '    double d2gdsdt['+str(ns)+'],d2gdsdp['+str(ns)+ \
                 '],d3gdsdtdp['+str(ns)+'],d3gds2dt['+str(ns)+']['+str(ns)+ \
                 '],d3gds2dp['+str(ns)+']['+str(ns)+'],d3gds3['+str(ns)+']['+ \
                 str(ns)+']['+str(ns)+'];\n'
        cb_12 += moles_assign_text
        cb_13  = '    double d2gdsdp['+str(ns)+'],d3gdsdp2['+str(ns)+ \
                 '],d3gds2dp['+str(ns)+']['+str(ns)+'],d3gds3['+str(ns)+']['+ \
                 str(ns)+ ']['+str(ns)+'];\n'
        cb_13 += moles_assign_text

        dmin = self.de_minimis
        if self.verbose:
            print ("... ... c["+str(nc)+"]: ", end="")
        for i in range(0,nc):
            if self.verbose:
                print (".", end="")
            for j in range(0,ns):
                # d2gdnds[{NC}][{NS}]
                a = '    d2gdnds[' + str(i) + '][' + str(j) + '] = '
                a += self.printer.doprint(
                    sym.diff(sFunc.function[j], n[i])) + ';\n'
                cb_5  += a
                cb_8  += a
                cb_9  += a
                cb_10 += a
                # d3gdndsdt[{NC}][{NS}]
                a = '    d3gdndsdt[' + str(i) + '][' + str(j) + '] = '
                a += '0' if dmin else self.printer.doprint(
                    sym.diff(sFunc.function[j], n[i], T))
                a += ';\n'
                cb_9 += a
                # d3gdndsdp[{NC}][{NS}]
                a = '    d3gdndsdp[' + str(i) + '][' + str(j) + '] = '
                a += '0' if dmin else self.printer.doprint(
                    sym.diff(sFunc.function[j], n[i], P))
                a += ';\n'
                cb_10 += a
                for k in range(i,nc):
                    # d3gdn2ds[{NC}][{NC}][{NS}]
                    a = '    d3gdn2ds[' + str(i) + '][' + str(k) + '][' \
                      + str(j) + '] = '
                    a += '0' if dmin else self.printer.doprint(
                        sym.diff(sFunc.function[j], n[i], n[k]))
                    a += ';\n'
                    cb_8 += a
                    if k > i:
                        cb_8 += '    d3gdn2ds[' + str(k) + '][' + str(i) + \
                                '][' + str(j) + '] = d3gdn2ds[' + str(i) + \
                                '][' + str(k) + '][' + str(j) + '];\n'
                for k in range(j,ns):
                    # d3gdnds2[{NC}][{NS}][{NS}]
                    a = '    d3gdnds2[' + str(i) + '][' + str(j) + '][' \
                      + str(k) + '] = '
                    a += '0' if dmin else self.printer.doprint(
                        sym.diff(sFunc.function[j], n[i], s[k]))
                    a += ';\n'
                    cb_8  += a
                    cb_9  += a
                    cb_10 += a
                    if k > j:
                        cb_8 += '    d3gdnds2[' + str(i) + '][' + str(k) + \
                                '][' + str(j) + '] = d3gdnds2[' + str(i) + \
                                '][' + str(j) + '][' + str(k) + '];\n'
                        cb_9 += '    d3gdnds2[' + str(i) + '][' + str(k) + \
                                '][' + str(j) + '] = d3gdnds2[' + str(i) + \
                                '][' + str(j) + '][' + str(k) + '];\n'
                        cb_10 += '    d3gdnds2[' + str(i) + '][' + str(k) + \
                                 '][' + str(j) + '] = d3gdnds2[' + str(i) + \
                                 '][' + str(j) + '][' + str(k) + '];\n'

        if self.verbose:
            print ("")
            print ("... ... s["+str(ns)+"]: ", end="")
        for i in range(0,ns):
            print (".", end="")
            # d2gdsdt[{NS}]
            a = '    d2gdsdt[' + str(i) + '] = '
            a += self.printer.doprint(sym.diff(sFunc.function[i], T)) + ';\n'
            cb_6  += a
            cb_9  += a
            cb_11 += a
            cb_12 += a
            # d2gdsdp[{NS}]
            a = '    d2gdsdp[' + str(i) + '] = '
            a += self.printer.doprint(sym.diff(sFunc.function[i], P)) + ';\n'
            cb_7  += a
            cb_10 += a
            cb_12 += a
            cb_13 += a
            # d3gdsdt2[{NS}]
            a = '    d3gdsdt2[' + str(i) + '] = '
            a += '0' if dmin else self.printer.doprint(
                sym.diff(sFunc.function[i], T, 2))
            a += ';\n'
            cb_11 += a
            # d3gdsdtdp[{NS}]
            a = '    d3gdsdtdp[' + str(i) + '] = '
            a += '0' if dmin else self.printer.doprint(
                sym.diff(sFunc.function[i], T, P))
            a += ';\n'
            cb_12 += a
            # d3gdsdp2[{NS}]
            a = '    d3gdsdp2[' + str(i) + '] = '
            a += '0' if dmin else self.printer.doprint(
                sym.diff(sFunc.function[i], P, 2))
            a += ';\n'
            cb_13 += a
            for j in range(i,ns):
                # d3gds2dt[{NS}][{NS}]
                a  = '    d3gds2dt[' + str(i) + '][' + str(j) + '] = '
                a += '0' if dmin else self.printer.doprint(
                    sym.diff(sFunc.function[i], s[j], T))
                a += ';\n'
                cb_9  += a
                cb_11 += a
                cb_12 += a
                if j > i:
                    cb_9 += '    d3gds2dt[' + str(j) + '][' + str(i) + '] = ' \
                          + 'd3gds2dt[' + str(i) + '][' + str(j) + '];\n'
                    cb_11 += '    d3gds2dt[' + str(j) + '][' + str(i) + '] = ' \
                           + 'd3gds2dt[' + str(i) + '][' + str(j) + '];\n'
                    cb_12 += '    d3gds2dt[' + str(j) + '][' + str(i) + '] = ' \
                           + 'd3gds2dt[' + str(i) + '][' + str(j) + '];\n'
                # d3gds2dp[{NS}][{NS}]
                a = '    d3gds2dp[' + str(i) + '][' + str(j) + '] = '
                a += '0' if dmin else self.printer.doprint(
                    sym.diff(sFunc.function[i], s[j], P))
                a += ';\n'
                cb_10 += a
                cb_12 += a
                cb_13 += a
                if j > i:
                    cb_10 += '    d3gds2dp[' + str(j) + '][' + str(i) + '] = ' \
                           + 'd3gds2dp[' + str(i) + '][' + str(j) + '];\n'
                    cb_12 += '    d3gds2dp[' + str(j) + '][' + str(i) + '] = ' \
                           + 'd3gds2dp[' + str(i) + '][' + str(j) + '];\n'
                    cb_13 += '    d3gds2dp[' + str(j) + '][' + str(i) + '] = ' \
                           + 'd3gds2dp[' + str(i) + '][' + str(j) + '];\n'
                for k in range(j,ns):
                    # d3gds3[{NS}][{NS}][{NS}]
                    a  = '    d3gds3[' + str(i) + '][' + str(j) + '][' + \
                         str(k) + '] = '
                    a += '0' if dmin else self.printer.doprint(
                        sym.diff(sFunc.function[i], s[j], s[k]))
                    a += ';\n'
                    a += self.exchange_indices(i,j,k,'d3gds3')
                    cb_8  += a
                    cb_9  += a
                    cb_10 += a
                    cb_11 += a
                    cb_12 += a
                    cb_13 += a
        if self.verbose:
            print ("")
            print ("... ... fill template and return")
        return tpl.create_ordering_code_template().format(
            NS=ns, NC=nc, MAX_ITER=iter_max,
            ORDER_CODE_BLOCK_ZERO=cb_0,
            ORDER_CODE_BLOCK_ONE=cb_1,
            ORDER_CODE_BLOCK_TWO=cb_2,
            ORDER_CODE_BLOCK_THREE=cb_3,
            ORDER_CODE_BLOCK_FOUR=cb_4,
            ORDER_CODE_BLOCK_FIVE=cb_5,
            ORDER_CODE_BLOCK_SIX=cb_6,
            ORDER_CODE_BLOCK_SEVEN=cb_7,
            ORDER_CODE_BLOCK_EIGHT=cb_8,
            ORDER_CODE_BLOCK_NINE=cb_9,
            ORDER_CODE_BLOCK_TEN=cb_10,
            ORDER_CODE_BLOCK_ELEVEN=cb_11,
            ORDER_CODE_BLOCK_TWELVE=cb_12,
            ORDER_CODE_BLOCK_THIRTEEN=cb_13)

    def exchange_indices(self, i, j, k, var):
        """
        Reassigns indices in a derivative tensor to generate dependent entries.

        Parameters
        ----------
        i,j,k : int
            Indices, i <= j <= k, of a tensor derivative

        var : str
            Name of tensor

        Returns
        -------
        output : str
            Code fragment (in C) that assigns name[i][j][k] to redundant 
            elements of the tensor, e.g. name[k][j][i]

        Notes
        -----
        The user does not normally call this function directly.

        """
        s = ''
        if i == j == k:
            return s
        elif i == j:
            # i,k,j and k,i,j must equal i,j,k
            s += '    '+var+'['+str(i)+']['+str(k)+']['+str(j)+'] = ' \
                       +var+'['+str(i)+']['+str(j)+']['+str(k)+'];\n'
            s += '    '+var+'['+str(k)+']['+str(i)+']['+str(j)+'] = ' \
                       +var+'['+str(i)+']['+str(j)+']['+str(k)+'];\n'
            return s
        elif j == k:
            # j,i,k and j,k,i must equal i,j,k
            s += '    '+var+'['+str(j)+']['+str(i)+']['+str(k)+'] = ' \
                       +var+'['+str(i)+']['+str(j)+']['+str(k)+'];\n'
            s += '    '+var+'['+str(j)+']['+str(k)+']['+str(i)+'] = ' \
                       +var+'['+str(i)+']['+str(j)+']['+str(k)+'];\n'
            return s
        else:
            # j,i,k and k,i,j and j,k,i and k,j,i must equal i,j,k
            s += '    '+var+'['+str(j)+']['+str(i)+']['+str(k)+'] = ' \
                       +var+'['+str(i)+']['+str(j)+']['+str(k)+'];\n'
            s += '    '+var+'['+str(k)+']['+str(i)+']['+str(j)+'] = ' \
                       +var+'['+str(i)+']['+str(j)+']['+str(k)+'];\n'
            s += '    '+var+'['+str(j)+']['+str(k)+']['+str(i)+'] = ' \
                       +var+'['+str(i)+']['+str(j)+']['+str(k)+'];\n'
            s += '    '+var+'['+str(k)+']['+str(j)+']['+str(i)+'] = ' \
                       +var+'['+str(i)+']['+str(j)+']['+str(k)+'];\n'
            return s

    def create_soln_calc_h_file(self, language='C', module_type='fast', 
        minimal_deriv_set=True):
        """
        Creates an include file implementing model calculations.

        Note that this include file contains code for functions that implement 
        the model for the generic case. It is meant to be included into a file 
        that implements a specific parameterized instance of the model. See 
        create_code_module().

        The user does not normally call this function directly.

        Parameters
        ----------
        language : string
            Language syntax for generated code. ("C" is the C99 programming 
            language.)
        module_type : string
            Generate code that executes "fast" but does not expose hooks for 
            model parameter calibration. Alternately, generate code suitable for 
            "calib"ration of parameters in the model, which executes more slowly 
            and exposes additional functions that allow setting of parameters 
            and generation of derivatives of thermodynamic functions with 
            respect to model parameters. 
        minimal_deriv_set : bool
            Generate a minimal set of compositional derivatives: dgdn, d2gdndt, 
            d2gdndp, d3gdndt2, d3gdndtdp, d3gdndp2, d4gdndt3, d4gdndt2dp, 
            d4gdndtdp2, d4gdndp3, d2gdn2, d3gdn2dt, d3gdn2dp, d3gdn3.  This is 
            the subset of derivatives currently required for solution phases 
            that are imported into the phases module. Remaining derivatives are
            returned with values of zero. 
        """
        assert language == "C"
        # Code for assigning mole numbers to components
        if self.verbose:
            print ("... entering create_soln_calc_h_file")
        nc = self.nc
        moles_assign_text = ''
        for i in range(0,nc):
            moles_assign_text += ('    double n' + str(i+1) + ' = n[' 
                + str(i) + '];')
            if i < nc-1:
                moles_assign_text += '\n'

        gen_expr = []
        w = '#include <math.h>\n\n'
        if self._include_dh_code:
            w += tpl.create_code_for_dh_functions(language)
            w += "\n"

        if self.verbose:
            print ("... generating ordering method")
        ns = self.ns
        order_assign_text = ''
        for i in range(0,ns):
            order_assign_text += ('    double s' + str(i+1) + ' = s[' 
                + str(i) + '];')
            if i < ns-1:
                order_assign_text += '\n'
        if ns > 1:
            w += tpl.create_ordering_gaussj_template(language).format(
                NS=ns)
        if ns > 0:
            w += '#include <float.h>\n\n'
            w += self.create_ordering_code(moles_assign_text+'\n',
                                           order_assign_text+'\n')

        n = self.n
        s = self.s 
        T = self.get_symbol_for_t()
        P = self.get_symbol_for_p()
        if self.verbose:
            print ("... generating derivatives of G")
        for idx, (f, oT, oP) in enumerate(self.g_list):
            # Potential and T, P derivatives
            a = sym.S.Zero
            for exp in self.expression_parts:
                assert exp.exp_type == 'unrestricted'
                a += exp.expression
            a = sym.diff(a, T, oT, P, oP)
            gen_expr.append(a)
        self._g_matrix = sym.Matrix(gen_expr)

        sFunc = self._ordering_functions
        order_initial_guess_text = ''
        separator = '        '
        for i in range(0,ns):
            order_initial_guess_text += separator
            order_initial_guess_text += 's[' + str(i) + '] = '
            order_initial_guess_text += self.printer.doprint(sFunc.guess[i])
            order_initial_guess_text += ';\n'

        complx_soln_calc_template = tpl.create_complx_soln_calc_template()
        if self.verbose:
            print ("... generating derivative strings for template")
        G = self.expression_parts[0].expression
        dgdn_code_text = ''
        d2gdndt_fill_text = ''
        d2gdndp_fill_text = ''
        d2gdn2_fill_text = ''
        d2gdnds_fill_text = ''
        d3gdn3_fill_text = ''
        d3gdn2dt_fill_text = ''
        d3gdn2dp_fill_text = ''
        d3gdn2ds_fill_text = ''
        d3gdndt2_fill_text = ''
        d3gdndtdp_fill_text = ''
        d3gdndp2_fill_text = ''
        d3gdndsdt_fill_text = ''
        d3gdndsdp_fill_text = ''
        d3gdnds2_fill_text = ''
        if self.verbose:
            print ("... ... c["+str(nc)+"]: ", end="")
        for i in range(0,nc):
            if self.verbose:
                print (".", end="")
            dgdn_code_text += '    dgdn['+str(i)+'] = '
            dgdn_code_text += self.printer.doprint(G.diff(n[i]))
            dgdn_code_text += ';\n'
            d2gdndt_fill_text += '    d2gdndt['+str(i)+'] = '
            d2gdndt_fill_text += self.printer.doprint(G.diff(T).diff(n[i]))
            d2gdndt_fill_text += ';\n'
            d2gdndp_fill_text += '    d2gdndp['+str(i)+'] = '
            d2gdndp_fill_text += self.printer.doprint(G.diff(P).diff(n[i]))
            d2gdndp_fill_text += ';\n'
            d3gdndt2_fill_text += '    d3gdndt2['+str(i)+'] = '
            d3gdndt2_fill_text += "0" if self.de_minimis else self.printer.doprint(
                G.diff(T,2).diff(n[i]))
            d3gdndt2_fill_text += ';\n'
            d3gdndtdp_fill_text += '    d3gdndtdp['+str(i)+'] = '
            d3gdndtdp_fill_text += "0" if self.de_minimis else self.printer.doprint(
                G.diff(T).diff(P).diff(n[i]))
            d3gdndtdp_fill_text += ';\n'
            d3gdndp2_fill_text += '    d3gdndp2['+str(i)+'] = '
            d3gdndp2_fill_text += "0" if self.de_minimis else self.printer.doprint(
                G.diff(P,2).diff(n[i]))
            d3gdndp2_fill_text += ';\n'
            for j in range(i,nc):
                d2gdn2_fill_text += '    d2gdn2['+str(i)+']['+str(j)+'] = '
                d2gdn2_fill_text += self.printer.doprint(G.diff(n[i]).diff(n[j]))
                d2gdn2_fill_text += ';\n'
                d3gdn2dt_fill_text += '    d3gdn2dt['+str(i)+']['+str(j)+'] = '
                d3gdn2dt_fill_text += "0" if self.de_minimis else self.printer.doprint(
                    G.diff(T).diff(n[i]).diff(n[j]))
                d3gdn2dt_fill_text += ';\n'
                d3gdn2dp_fill_text += '    d3gdn2dp['+str(i)+']['+str(j)+'] = '
                d3gdn2dp_fill_text += "0" if self.de_minimis else self.printer.doprint(
                    G.diff(P).diff(n[i]).diff(n[j]))
                d3gdn2dp_fill_text += ';\n'
                if j > i:
                    d2gdn2_fill_text += '    d2gdn2['+str(j)+']['+str(i)+'] = '
                    d2gdn2_fill_text +=     'd2gdn2['+str(i)+']['+str(j)+'];\n'
                    d3gdn2dt_fill_text += '    d3gdn2dt['+str(j)+']['+str(i)+'] = '
                    d3gdn2dt_fill_text +=     'd3gdn2dt['+str(i)+']['+str(j)+'];\n'
                    d3gdn2dp_fill_text += '    d3gdn2dp['+str(j)+']['+str(i)+'] = '
                    d3gdn2dp_fill_text +=     'd3gdn2dp['+str(i)+']['+str(j)+'];\n'
                for k in range(j,nc):
                    d3gdn3_fill_text += '    d3gdn3['+str(i)+']['+str(j)+'][' \
                                      +str(k)+'] = '
                    d3gdn3_fill_text += "0" if self.de_minimis else self.printer.doprint(
                        G.diff(n[i]).diff(n[j]).diff(n[k]))
                    d3gdn3_fill_text += ';\n'
                    d3gdn3_fill_text += self.exchange_indices(i,j,k,'d3gdn3')
                for k in range(0,ns):
                    d3gdn2ds_fill_text += '    d3gdn2ds['+str(i)+']['+str(j)+'][' \
                                        +str(k)+'] = '
                    d3gdn2ds_fill_text += "0" if self.de_minimis else self.printer.doprint(
                        G.diff(n[i]).diff(n[j]).diff(s[k]))
                    d3gdn2ds_fill_text += ';\n'
                    if j > i:
                        d3gdn2ds_fill_text += '    d3gdn2ds['+str(j)+'][' \
                                            +str(i)+']['+str(k)+'] = '
                        d3gdn2ds_fill_text +=     'd3gdn2ds['+str(i)+'][' \
                                            +str(j)+']['+str(k)+'];\n'

            for j in range(0,ns):
                d2gdnds_fill_text += '    d2gdnds['+str(i)+']['+str(j)+'] = '
                d2gdnds_fill_text += self.printer.doprint(G.diff(n[i]).diff(s[j]))
                d2gdnds_fill_text += ';\n'
                d3gdndsdt_fill_text += '    d3gdndsdt['+str(i)+']['+str(j)+'] = '
                d3gdndsdt_fill_text += "0" if self.de_minimis else self.printer.doprint(
                    G.diff(T).diff(n[i]).diff(s[j]))
                d3gdndsdt_fill_text += ';\n'
                d3gdndsdp_fill_text += '    d3gdndsdp['+str(i)+']['+str(j)+'] = '
                d3gdndsdp_fill_text += "0" if self.de_minimis else self.printer.doprint(
                    G.diff(P).diff(n[i]).diff(s[j]))
                d3gdndsdp_fill_text += ';\n'
                for k in range(0,ns):
                    d3gdnds2_fill_text += '    d3gdnds2['+str(i)+']['+str(j) \
                                        +']['+str(k)+'] = '
                    d3gdnds2_fill_text += "0" if self.de_minimis else self.printer.doprint(
                        G.diff(n[i]).diff(s[j]).diff(s[k]))
                    d3gdnds2_fill_text += ';\n'
                    if k > j:
                        d3gdnds2_fill_text += '    d3gdnds2['+str(i)+'][' \
                                            +str(k)+']['+str(j)+'] = '
                        d3gdnds2_fill_text +=     'd3gdnds2['+str(i)+'][' \
                                            +str(j)+']['+str(k)+'];\n'
        if self.verbose:
            print ("")
            print ("... ... s["+str(ns)+"]: ", end="")
        d2gdsdt_fill_text = ''
        d2gdsdp_fill_text = ''
        d2gds2_fill_text = ''
        d3gdsdt2_fill_text = ''
        d3gdsdtdp_fill_text = ''
        d3gdsdp2_fill_text = ''
        d3gds2dt_fill_text = ''
        d3gds2dp_fill_text = ''
        d3gds3_fill_text = ''
        for i in range(0,ns):
            if self.verbose:
                print (".", end="")
            d2gdsdt_fill_text += '    d2gdsdt['+str(i)+'] = '
            d2gdsdt_fill_text += self.printer.doprint(G.diff(T).diff(s[i]))
            d2gdsdt_fill_text += ';\n'
            d2gdsdp_fill_text += '    d2gdsdp['+str(i)+'] = '
            d2gdsdp_fill_text += self.printer.doprint(G.diff(P).diff(s[i]))
            d2gdsdp_fill_text += ';\n'
            d3gdsdt2_fill_text += '    d3gdsdt2['+str(i)+'] = '
            d3gdsdt2_fill_text += "0" if self.de_minimis else self.printer.doprint(
                G.diff(T,2).diff(s[i]))
            d3gdsdt2_fill_text += ';\n'
            d3gdsdtdp_fill_text += '    d3gdsdtdp['+str(i)+'] = '
            d3gdsdtdp_fill_text += "0" if self.de_minimis else self.printer.doprint(
                G.diff(T).diff(P).diff(s[i]))
            d3gdsdtdp_fill_text += ';\n'
            d3gdsdp2_fill_text += '    d3gdsdp2['+str(i)+'] = '
            d3gdsdp2_fill_text += "0" if self.de_minimis else self.printer.doprint(
                G.diff(P,2).diff(s[i]))
            d3gdsdp2_fill_text += ';\n'
            for j in range(i,ns):
                d2gds2_fill_text += '    d2gds2['+str(i)+']['+str(j)+'] = '
                d2gds2_fill_text += self.printer.doprint(G.diff(s[i]).diff(s[j]))
                d2gds2_fill_text += ';\n'
                d3gds2dt_fill_text += '    d3gds2dt['+str(i)+']['+str(j)+'] = '
                d3gds2dt_fill_text += "0" if self.de_minimis else self.printer.doprint(
                    G.diff(T).diff(s[i]).diff(s[j]))
                d3gds2dt_fill_text += ';\n'
                d3gds2dp_fill_text += '    d3gds2dp['+str(i)+']['+str(j)+'] = '
                d3gds2dp_fill_text += "0" if self.de_minimis else self.printer.doprint(
                    G.diff(P).diff(s[i]).diff(s[j]))
                d3gds2dp_fill_text += ';\n'
                if j > i:
                    d2gds2_fill_text += '    d2gds2['+str(j)+']['+str(i)+'] = '
                    d2gds2_fill_text +=     'd2gds2['+str(i)+']['+str(j)+'];\n'
                    d3gds2dt_fill_text += '    d3gds2dt['+str(j)+']['+str(i)+'] = '
                    d3gds2dt_fill_text +=     'd3gds2dt['+str(i)+']['+str(j)+'];\n'
                    d3gds2dp_fill_text += '    d3gds2dp['+str(j)+']['+str(i)+'] = '
                    d3gds2dp_fill_text +=     'd3gds2dp['+str(i)+']['+str(j)+'];\n'
                for k in range(j,ns):
                    d3gds3_fill_text += '    d3gds3['+str(j)+']['+str(i)+ \
                                        ']['+str(k)+'] = '
                    d3gds3_fill_text += "0" if self.de_minimis else self.printer.doprint(
                        G.diff(s[i]).diff(s[j]).diff(s[k]))
                    d3gds3_fill_text += ';\n'
                    d3gds3_fill_text += self.exchange_indices(i,j,k,'d3gds3')

        if self.verbose:
            print ("")
            print ("... generating string using template format")
        w += complx_soln_calc_template.format(
            module=self.module,
            number_components=nc,
            number_ordering=ns,
            number_symmetric_hessian_terms=int(nc+nc*(nc-1)/2),
            number_symmetric_tensor_terms=int(nc*(nc+1)*(nc+2)/6), 
            moles_assign=moles_assign_text,
            order_assign=order_assign_text,
            order_initial_guess=order_initial_guess_text,
            g_code=self.printer.doprint(G),
            dgdt_code=self.printer.doprint(G.diff(T)),
            dgdp_code=self.printer.doprint(G.diff(P)),
            dgdn_code=dgdn_code_text,
            d2gdt2_code="0" if self.de_minimis else self.printer.doprint(
                G.diff(T,2)),
            d2gdtdp_code="0" if self.de_minimis else self.printer.doprint(
                G.diff(T).diff(P)),
            d2gdp2_code="0" if self.de_minimis else self.printer.doprint(
                G.diff(P,2)),
            d3gdt3_code="0" if self.de_minimis else self.printer.doprint(
                G.diff(T,3)),
            d3gdt2dp_code="0" if self.de_minimis else self.printer.doprint(
                G.diff(T,2).diff(P)),
            d3gdtdp2_code="0" if self.de_minimis else self.printer.doprint(
                G.diff(T).diff(P,2)),
            d3gdp3_code="0" if self.de_minimis else self.printer.doprint(
                G.diff(P,3)),
            fillD2GDNDT=d2gdndt_fill_text,
            fillD2GDNDP=d2gdndp_fill_text,
            fillD2GDN2=d2gdn2_fill_text,
            fillD2GDNDS=d2gdnds_fill_text,
            fillD2GDSDT=d2gdsdt_fill_text,
            fillD2GDSDP=d2gdsdp_fill_text,
            fillD2GDS2=d2gds2_fill_text,
            fillD3GDN3=d3gdn3_fill_text,
            fillD3GDN2DT=d3gdn2dt_fill_text,
            fillD3GDN2DP=d3gdn2dp_fill_text,
            fillD3GDN2DS=d3gdn2ds_fill_text,
            fillD3GDNDT2=d3gdndt2_fill_text,
            fillD3GDNDTDP=d3gdndtdp_fill_text,
            fillD3GDNDP2=d3gdndp2_fill_text,
            fillD3GDNDSDT=d3gdndsdt_fill_text,
            fillD3GDNDSDP=d3gdndsdp_fill_text,
            fillD3GDNDS2=d3gdnds2_fill_text,
            fillD3GDSDT2=d3gdsdt2_fill_text,
            fillD3GDSDTDP=d3gdsdtdp_fill_text,
            fillD3GDSDP2=d3gdsdp2_fill_text,
            fillD3GDS2DT=d3gds2dt_fill_text,
            fillD3GDS2DP=d3gds2dp_fill_text,
            fillD3GDS3=d3gds3_fill_text
            )

        if self.verbose:
            print ("... generating convenience code")
        convenience_template = tpl.create_soln_redun_template()
        w += convenience_template.format(module=self.module, 
            number_components=nc)
        if self.verbose:
            print ("... exiting create_soln_calc_h_file")
        return w

    def create_soln_calib_h_file(self, language='C'):
        """
        Creates an include file implementing model parameter derivatives.

        Note that this include file contains code for functions that implement 
        the model for the generic case. It is meant to be included into a file 
        that implements a specific parameterized instance of the model. See 
        create_code_module().

        The user does not normally call this function directly.

        Parameters
        ----------
        language : string
            Language syntax for generated code. ("C" is the C99 programming 
            language.)
        """
        assert language == "C"
        SpMdClass = True if self.__class__.__name__ == 'SpeciationSolnModel' \
            else False
        if self.verbose:
            print ("... entering create_soln_calib_h_file")

        # Code for assigning mole numbers to components
        c = self.nc
        moles_assign_text = ''
        for i in range(0,c):
            moles_assign_text += ('    double n' + str(i+1) + ' = n[' 
                + str(i) + '];')
            if i < c-1:
                moles_assign_text += '\n'

        ns = self.ns

        sFunc = self._ordering_functions
        order_initial_guess_text = ''
        separator = '        '
        for i in range(0,ns):
            order_initial_guess_text += separator
            order_initial_guess_text += 's[' + str(i) + '] = '
            order_initial_guess_text += self.printer.doprint(sFunc.guess[i])
            order_initial_guess_text += ';\n'

        order_assign_text = ''
        for i in range(0,ns):
            order_assign_text += ('    double s' + str(i+1) + ' = s[' 
                + str(i) + '];')
            if i < ns-1:
                order_assign_text += '\n'

        solution_calib_template = tpl.create_complx_soln_calib_template()
        symparam = self.get_model_param_symbols()
        G_param_jac = sym.Matrix(self._g_matrix).jacobian(symparam)

        w = '#include <math.h>\n'

        eCB = ''
        np = len(self.params)
        G = self.expression_parts[0].expression
        n = self.n
        s = self.s 
        T = self.get_symbol_for_t()
        P = self.get_symbol_for_p()

        dmin = self.de_minimis
        
        ##########################
        # Block of code for dgdw #
        ##########################
        if self.verbose:
            print ("... ( 1/18) computing and writing code for dgdw")
        j = 0
        G_jac_list = [self.printer.doprint(G_param_jac[j,i]) for i in range(
                0, np)]

        switch_code_text = ''
        for i in range(0, len(self.params)):
            switch_code_text += '    case ' + str(i) + ':\n'
            switch_code_text += '        result += ' + G_jac_list[i] + ';\n'
            switch_code_text += '        break;\n'
            
        w += solution_calib_template.format(module=self.module, 
            number_components=c, func=self.g_list[j][0], 
            switch_code=switch_code_text, moles_assign=moles_assign_text,
            number_ordering=ns, 
            order_initial_guess=order_initial_guess_text,
            order_assign=order_assign_text,
            extra_ordering_code=eCB)

        #######################
        # Extra ordering code #
        # dsdw                #
        #######################
        if self.verbose:
            print ("... ( 2/18) computing and writing code for dsdw")
        eCB = 'static void order_dsdw(double T, double P,' \
            + ' double n['+str(c)+'],' \
            + (' double b['+str(c)+'],' if SpMdClass else '') \
            + ' double s['+str(ns)+'],' \
            + ' double invd2gds2['+str(ns)+']['+str(ns)+'],' \
            + ' double dsdw['+str(ns)+']['+str(np)+']) {\n'
        eCB += '    int i,j,k;\n'
        eCB += '    double d2gdsdw['+str(ns)+']['+str(np)+'];\n'
        eCB += moles_assign_text + '\n'
        if SpMdClass:
            for i in range(0,c):
                eCB += '    double b' + str(i+1) + ' = b[' + str(i) + '];\n'
        eCB += order_assign_text + '\n'
        for k in range(0,ns):
            f = sFunc.function[k]
            for l in range (0,np):
                eCB += '    d2gdsdw['+str(k)+']['+str(l) \
                     + '] = ' + self.printer.doprint(f.diff(symparam[l])) \
                     + ';\n'
        eCB += '    for (i=0; i<'+str(np)+'; i++) {\n'
        eCB += '        for (j=0; j<'+str(ns)+'; j++) {\n'
        eCB += '            dsdw[j][i] = 0.0;\n'
        eCB += '            for (k=0; k<'+str(ns)+'; k++) dsdw[j][i] += ' \
             + '- invd2gds2[j][k]*d2gdsdw[k][i];\n'
        eCB += '        }\n'
        eCB += '    }\n'
        eCB += '}'

        #############################
        # Block of code for d2gdtdw #
        #############################
        if self.verbose:
            print ("... ( 3/18) computing and writing code for d2gdtdw")
        extraCb  = '    double d2gdsdt['+str(ns)+'];\n'
        extraCb += '    double d2gdsdw['+str(ns)+']['+str(np)+'];\n'
        extraCb += '    double d2gds2['+str(ns)+']['+str(ns)+'];\n'
        extraCb += '    double dsdt['+str(ns)+'];\n'
        extraCb += '    double dsdw['+str(ns)+']['+str(np)+'];\n'
        extraCb += '    int i,j;\n'
        extraCb += '    order_dsdt(T, P, n, s, invd2gds2, dsdt);\n'
        extraCb += '    order_dsdw(T, P, n, s, invd2gds2, dsdw);\n'
        extraCa  = ''
        for k in range(0,ns):
            extraCa += '    d2gdsdt['+str(k)+'] = '
            extraCa += '0' if dmin else self.printer.doprint(
                G.diff(T).diff(s[k]))
            extraCa += ';\n'
            for l in range(0,ns):
                extraCa += '    d2gds2['+str(k)+']['+str(l)+'] = '
                extraCa += '0' if dmin else self.printer.doprint(
                    G.diff(s[k]).diff(s[l]))
                extraCa += ';\n'
            for l in range(0,np):
                extraCa += '    d2gdsdw['+str(k)+']['+str(l)+'] = '
                extraCa += '0' if dmin else self.printer.doprint(
                    G.diff(s[k]).diff(symparam[l]))
                extraCa += ';\n'

        j += 1
        if dmin:
            G_jac_list = ['0' for i in range(0, np)]
        else:
            G_jac_list = [self.printer.doprint(
                G_param_jac[j,i]) for i in range(0, np)]

        switch_code_text = ''
        for i in range(0, len(self.params)):
            switch_code_text += '    case ' + str(i) + ':\n'
            switch_code_text += '        result += ' + G_jac_list[i] + ';\n'
            switch_code_text += '        for (i=0; i<'+str(ns)+'; i++) {\n'
            switch_code_text += '            result += d2gdsdt[i]*dsdw[i][' \
                +str(i)+'] + d2gdsdw[i]['+str(i)+']*dsdt[i];\n'
            switch_code_text += '            for (j=0; j<'+str(ns)+'; j++) ' \
                + 'result += d2gds2[i][j]*dsdt[i]*dsdw[j]['+str(i)+'];\n'
            switch_code_text += '        }\n'
            switch_code_text += '        break;\n'
            
        w += solution_calib_template.format(module=self.module, 
            number_components=c, func=self.g_list[j][0], 
            switch_code=switch_code_text, moles_assign=moles_assign_text,
            number_ordering=ns, 
            order_initial_guess=order_initial_guess_text,
            order_assign=extraCb+order_assign_text+extraCa,
            extra_ordering_code=eCB)

        #############################
        # Block of code for d2gdpdw #
        #############################
        if self.verbose:
            print ("... ( 4/18) computing and writing code for d2gdpdw")
        extraCb  = '    double d2gdsdp['+str(ns)+'];\n'
        extraCb += '    double d2gdsdw['+str(ns)+']['+str(np)+'];\n'
        extraCb += '    double d2gds2['+str(ns)+']['+str(ns)+'];\n'
        extraCb += '    double dsdp['+str(ns)+'];\n'
        extraCb += '    double dsdw['+str(ns)+']['+str(np)+'];\n'
        extraCb += '    int i,j;\n'
        extraCb += '    order_dsdp(T, P, n, s, invd2gds2, dsdp);\n'
        extraCb += '    order_dsdw(T, P, n, s, invd2gds2, dsdw);\n'
        extraCa  = ''
        for k in range(0,ns):
            extraCa += '    d2gdsdp['+str(k)+'] = '
            extraCa += '0' if dmin else self.printer.doprint(
                G.diff(P).diff(s[k]))
            extraCa += ';\n'
            for l in range(0,ns):
                extraCa += '    d2gds2['+str(k)+']['+str(l)+'] = '
                extraCa += '0' if dmin else self.printer.doprint(
                    G.diff(s[k]).diff(s[l]))
                extraCa += ';\n'
            for l in range(0,np):
                extraCa += '    d2gdsdw['+str(k)+']['+str(l)+'] = '
                extraCa += '0' if dmin else self.printer.doprint(
                    G.diff(s[k]).diff(symparam[l]))
                extraCa += ';\n'

        j += 1
        if dmin:
            G_jac_list = ['0' for i in range(0, np)]
        else:
            G_jac_list = [self.printer.doprint(
                G_param_jac[j,i]) for i in range(0, np)]

        switch_code_text = ''
        for i in range(0, len(self.params)):
            switch_code_text += '    case ' + str(i) + ':\n'
            switch_code_text += '        result += ' + G_jac_list[i] + ';\n'
            switch_code_text += '        for (i=0; i<'+str(ns)+'; i++) {\n'
            switch_code_text += '            result += d2gdsdp[i]*dsdw[i][' \
                +str(i)+'] + d2gdsdw[i]['+str(i)+']*dsdp[i];\n'
            switch_code_text += '            for (j=0; j<'+str(ns)+'; j++) ' \
                + 'result += d2gds2[i][j]*dsdp[i]*dsdw[j]['+str(i)+'];\n'
            switch_code_text += '        }\n'
            switch_code_text += '        break;\n'
            
        w += solution_calib_template.format(module=self.module, 
            number_components=c, func=self.g_list[j][0], 
            switch_code=switch_code_text, moles_assign=moles_assign_text,
            number_ordering=ns, 
            order_initial_guess=order_initial_guess_text,
            order_assign=extraCb+order_assign_text+extraCa,
            extra_ordering_code='')

        #############################
        # Block of code for d2gdndw #
        #############################
        if self.verbose:
            print ("... ( 5/18) computing and writing code for d2gdndw")
        extraCb  = '    double d2gdsdn['+str(ns)+']['+str(c)+'];\n'
        extraCb += '    double d2gdsdw['+str(ns)+']['+str(np)+'];\n'
        extraCb += '    double d2gds2['+str(ns)+']['+str(ns)+'];\n'
        extraCb += '    double dsdn['+str(ns)+']['+str(c)+'];\n'
        extraCb += '    double dsdw['+str(ns)+']['+str(np)+'];\n'
        extraCb += '    int i,j;\n'
        extraCb += '    order_dsdn(T, P, n, s, invd2gds2, dsdn);\n'
        extraCb += '    order_dsdw(T, P, n, s, invd2gds2, dsdw);\n'
        extraCa  = ''
        for k in range(0,ns):
            for l in range(0,c):
                extraCa += '    d2gdsdn['+str(k)+']['+str(l)+'] = ' \
                    + self.printer.doprint(G.diff(s[k]).diff(n[l])) + ';\n'
            for l in range(0,ns):
                extraCa += '    d2gds2['+str(k)+']['+str(l)+'] = ' \
                         + self.printer.doprint(G.diff(s[k]).diff(s[l])) + ';\n'
            for l in range(0,np):
                extraCa += '    d2gdsdw['+str(k)+']['+str(l)+'] = ' \
                  + self.printer.doprint(G.diff(s[k]).diff(symparam[l])) + ';\n'

        w += ('static void ' + self.module + '_dparam_dgdn(double T, double P, '
            +'double n[' + str(c) + '], int index, double result[' + str(c) 
            + ']) {\n')
        w += moles_assign_text + '\n'
        w += '    double invd2gds2['+str(ns)+']['+str(ns)+'];\n'
        w += '    double *s = retrieveGuess(T, P, n);\n'
        w += '    order_s(T, P, n, s, invd2gds2);\n'
        w += extraCb
        w += order_assign_text + '\n'
        w += extraCa
        switch_code_text = '    switch (index) {\n'
        for i in range(0, len(self.params)):
            a = self.expression_parts[0].expression.diff(symparam[i])
            switch_code_text += '    case ' + str(i) + ':\n'
            for jjj in range(0,c):
                switch_code_text += ('        result['+str(jjj)+'] = ' 
                    + self.printer.doprint(a.diff(n[jjj])) + ';\n')

                switch_code_text += '        for (i=0; i<'+str(ns)+'; i++) {\n'
                switch_code_text += '            result['+str(jjj)+']' \
                    + '+= d2gdsdn[i]['+str(jjj)+']*dsdw[i]['+str(i)+']' \
                    + ' + d2gdsdw[i]['+str(i)+']*dsdn[i]['+str(jjj)+'];\n'
                switch_code_text += '            for (j=0; j<'+str(ns)+'; j++) ' \
                    + 'result['+str(jjj)+'] += d2gds2[i][j]*dsdn[i]['+str(jjj)+']' \
                    + '*dsdw[j]['+str(i)+'];\n'
                switch_code_text += '        }\n'


            switch_code_text += '        break;\n'
        w += switch_code_text
        w += '    default:\n'
        w += '        break;\n'
        w += '    }\n'
        w += '}\n'

        #######################
        # Extra ordering code #
        # d2sdtdw             #
        #######################
        if self.verbose:
            print ("... ( 6/18) computing and writing code for d2sdtdw")
        eCB = 'static void order_d2sdtdw(double T, double P,' \
            + ' double n['+str(c)+'],' \
            + (' double b['+str(c)+'],' if SpMdClass else '') \
            + ' double s['+str(ns)+'],' \
            + ' double invd2gds2['+str(ns)+']['+str(ns)+'],' \
            + ' double d2sdtdw['+str(ns)+']['+str(np)+']) {\n'
        eCB += '    double dsdt['+str(ns)+'], dsdw['+str(ns)+']['+str(np)+'];\n'
        eCB += '    double temp['+str(ns)+'];\n'
        eCB += '    int i,j,k,l,ll;\n'
        eCB += '    double d2gdsdt['+str(ns)+'];\n'
        eCB += '    double d2gdsdw['+str(ns)+']['+str(np)+'];\n'
        eCB += '    double d3gdsdtdw['+str(ns)+']['+str(np)+'];\n'
        eCB += '    double d3gds2dt['+str(ns)+']['+str(ns)+'];\n'
        eCB += '    double d3gds2dw['+str(ns)+']['+str(ns)+']['+str(np)+'];\n'
        eCB += '    double d3gds3['+str(ns)+']['+str(ns)+']['+str(ns)+'];\n'
        eCB += moles_assign_text + '\n'
        if SpMdClass:
            for i in range(0,c):
                eCB += '    double b' + str(i+1) + ' = b[' + str(i) + '];\n'
        eCB += order_assign_text + '\n'
        for k in range(0,ns):
            f = sFunc.function[k]
            eCB += '    d2gdsdt['+str(k)+'] = '
            eCB += '0' if dmin else self.printer.doprint(f.diff(T))
            eCB += ';\n'
            for l in range(0,np):
                eCB += '    d2gdsdw['+str(k)+']['+str(l) + '] = '
                eCB += '0' if dmin else self.printer.doprint(
                    f.diff(symparam[l])) 
                eCB += ';\n'
                eCB += '    d3gdsdtdw['+str(k)+']['+str(l) + '] = '
                eCB += '0' if dmin else self.printer.doprint(
                    f.diff(symparam[l]).diff(T)) 
                eCB += ';\n'
            for l in range(k,ns):
                eCB += '    d3gds2dt['+str(k)+']['+str(l)+'] = '
                eCB += '0' if dmin else self.printer.doprint(
                    f.diff(s[l]).diff(T)) 
                eCB += ';\n'
                if l > k:
                    eCB += '    d3gds2dt['+str(l)+']['+str(k)+'] = '
                    eCB += 'd3gds2dt['+str(k)+']['+str(l)+'];\n'
                for ll in range(0,np):
                    eCB += '    d3gds2dw['+str(k)+']['+str(l)+']['+str(ll)+'] = '
                    eCB += '0' if dmin else self.printer.doprint(
                        f.diff(s[l]).diff(symparam[ll]))
                    eCB += ';\n'
                    if l > k:
                        eCB += '    d3gds2dw['+str(l)+']['+str(k)+']['+str(ll) \
                              +'] = d3gds2dw['+str(k)+']['+str(l)+']['+str(ll) \
                             + '];\n'
                for ll in range(l,ns):
                    eCB += '    d3gds3['+str(k)+']['+str(l)+']['+str(ll)+'] = '
                    eCB += '0' if dmin else self.printer.doprint(
                        f.diff(s[l]).diff(s[ll])) 
                    eCB += ';\n'
                    eCB += self.exchange_indices(k,l,ll,'d3gds3')
        eCB += '    for (i=0; i<'+str(ns)+'; i++) {\n'
        eCB += '        dsdt[i] = 0.0;\n'
        eCB += '        for (j=0; j<'+str(ns)+'; j++) dsdt[i] += '
        eCB +=            '- invd2gds2[i][j]*d2gdsdt[j];\n'
        eCB += '    }\n'
        eCB += '    for (ll=0; ll<'+str(np)+'; ll++) {\n'
        eCB += '        for (i=0; i<'+str(ns)+'; i++) {\n'
        eCB += '            dsdw[i][ll] = 0.0;\n'
        eCB += '            for (j=0; j<'+str(ns)+'; j++) '
        eCB +=              'dsdw[i][ll] += - invd2gds2[i][j]*d2gdsdw[j][ll];\n'
        eCB += '        }\n'
        eCB += '    }\n'
        eCB += '    for (ll=0; ll<'+str(np)+'; ll++) {\n'
        eCB += '        for (i=0; i<'+str(ns)+'; i++) {\n'
        eCB += '            for (j=0; j<'+str(ns)+'; j++) {\n'
        eCB += '                temp[j] = d3gdsdtdw[j][ll];\n'
        eCB += '                for (k=0; k<'+str(ns)+'; k++) {\n'
        eCB += '                    temp[j] += d3gds2dt[j][k]*dsdw[k][ll] '
        eCB +=                              '+ d3gds2dw[j][k][ll]*dsdt[k];\n'
        eCB += '                    for (l=0; l<'+str(ns)+'; l++) temp[j] += '
        eCB +=                         'd3gds3[j][k][l]*dsdt[k]*dsdw[l][ll];\n'
        eCB += '                }\n'
        eCB += '            }\n'
        eCB += '            d2sdtdw[i][ll] = 0.0;\n'
        eCB += '            for (j=0; j<'+str(ns)+'; j++) '
        eCB +=              'd2sdtdw[i][ll] += - invd2gds2[i][j]*temp[j];\n'
        eCB += '        }\n'
        eCB += '    }\n'
        eCB += '}\n'

        ##############################
        # Block of code for d3gdt2dw #
        ##############################
        if self.verbose:
            print ("... ( 7/18) computing and writing code for d3gdt2dw")
        extraCb  = '    double d2gdsdt['+str(ns)+'];\n'
        extraCb += '    double d2gdsdw['+str(ns)+']['+str(np)+'];\n'
        extraCb += '    double d2gds2['+str(ns)+']['+str(ns)+'];\n'
        extraCb += '    double dsdt['+str(ns)+'];\n'
        extraCb += '    double dsdw['+str(ns)+']['+str(np)+'];\n'
        extraCb += '    double d2sdt2['+str(ns)+'];\n'
        extraCb += '    double d2sdtdw['+str(ns)+']['+str(np)+'];\n'
        extraCb += '    double d3gds3['+str(ns)+']['+str(ns)+']['+str(ns)+'];\n'
        extraCb += '    double d3gds2dt['+str(ns)+']['+str(ns)+'];\n'
        extraCb += '    double d3gds2dw['+str(ns)+']['+str(ns)+']['+str(np)+'];\n'
        extraCb += '    double d3gdsdt2['+str(ns)+'];\n'
        extraCb += '    double d3gdsdtdw['+str(ns)+']['+str(np)+'];\n'
        extraCb += '    int i,j,k;\n'
        extraCb += '    order_dsdt(T, P, n, s, invd2gds2, dsdt);\n'
        extraCb += '    order_dsdw(T, P, n, s, invd2gds2, dsdw);\n'
        extraCb += '    order_d2sdt2(T, P, n, s, invd2gds2, d2sdt2);\n'
        extraCb += '    order_d2sdtdw(T, P, n, s, invd2gds2, d2sdtdw);\n'
        extraCa  = ''
        for k in range(0,ns):
            extraCa += '    d2gdsdt['+str(k)+'] = '
            extraCa += '0' if dmin else self.printer.doprint(
                G.diff(T).diff(s[k])) 
            extraCa += ';\n'
            extraCa += '    d3gdsdt2['+str(k)+'] = '
            extraCa += '0' if dmin else self.printer.doprint(
                G.diff(T,2).diff(s[k])) 
            extraCa += ';\n'
            for l in range(k,ns):
                extraCa += '    d2gds2['+str(k)+']['+str(l)+'] = '
                extraCa += '0' if dmin else self.printer.doprint(
                    G.diff(s[k]).diff(s[l])) 
                extraCa += ';\n'
                extraCa += '    d3gds2dt['+str(k)+']['+str(l)+'] = '
                extraCa += '0' if dmin else self.printer.doprint(
                    G.diff(s[k]).diff(s[l]).diff(T))
                extraCa += ';\n'
                if l > k:
                    extraCa += '    d2gds2['+str(k)+']['+str(l)+'] = '
                    extraCa +=    ' d2gds2['+str(l)+']['+str(k)+'];\n'
                    extraCa += '    d3gds2dt['+str(k)+']['+str(l)+'] = '
                    extraCa +=    ' d3gds2dt['+str(l)+']['+str(k)+'];\n'
                for ll in range(0,np):
                    extraCa += '    d3gds2dw['+str(k)+']['+str(l)+']['+str(ll)+'] = '
                    extraCa += '0' if dmin else self.printer.doprint(
                        G.diff(s[k]).diff(s[l]).diff(symparam[ll]))
                    extraCa += ';\n'
                    if l > k:
                        extraCa += '    d3gds2dw['+str(k)+']['+str(l)+']['+str(ll)+'] = '
                        extraCa +=    ' d3gds2dw['+str(l)+']['+str(k)+']['+str(ll)+'];\n'
                for ll in range(l,ns):
                    extraCa += '    d3gds3['+str(k)+']['+str(l)+']['+str(ll)+'] = '
                    extraCa += '0' if dmin else  self.printer.doprint(
                        G.diff(s[k]).diff(s[l]).diff(s[ll]))
                    extraCa += ';\n'
                    extraCa += self.exchange_indices(k,l,ll,'d3gds3')
            for l in range(0,np):
                extraCa += '    d2gdsdw['+str(k)+']['+str(l)+'] = '
                extraCa += '0' if dmin else self.printer.doprint(
                    G.diff(s[k]).diff(symparam[l])) 
                extraCa += ';\n'
                extraCa += '    d3gdsdtdw['+str(k)+']['+str(l)+'] = '
                extraCa += '0' if dmin else self.printer.doprint(
                    G.diff(s[k]).diff(symparam[l]).diff(T))
                extraCa += ';\n'

        j += 1
        if dmin:
            G_jac_list = ['0' for i in range(0, np)]
        else:
            G_jac_list = [self.printer.doprint(
                G_param_jac[j,i]) for i in range(0, np)]

        switch_code_text = ''
        for i in range(0, len(self.params)):
            switch_code_text += '    case ' + str(i) + ':\n'
            switch_code_text += '        result += ' + G_jac_list[i] + ';\n'

            switch_code_text += '        for (i=0; i<'+str(ns)+'; i++) {\n'
            switch_code_text += '            result += ' \
            'd3gdsdt2[i]*dsdw[i]['+str(i)+'] ' \
            + '+ 2.0*d2gdsdt[i]*d2sdtdw[i]['+str(i)+'] ' \
            + '+ d2gdsdw[i]['+str(i)+']*d2sdt2[i] ' \
            + '+ 2.0*d3gdsdtdw[i]['+str(i)+']*dsdt[i];\n'
            switch_code_text += '            for (j=0; j<'+str(ns)+'; j++) {\n'
            switch_code_text += '                result += ' \
            + '2.0*d3gds2dt[i][j]*dsdt[i]*dsdw[j]['+str(i)+']' \
            + '+ d2gds2[i][j]*d2sdt2[i]*dsdw[j]['+str(i)+']' \
            + '+ 2.0*d2gds2[i][j]*dsdt[i]*d2sdtdw[j]['+str(i)+']' \
            + '+ d3gds2dw[i][j]['+str(i)+']*dsdt[i]*dsdt[j];\n'
            switch_code_text += '                for (k=0; k<1; k++) ' \
            'result += d3gds3[i][j][k]*dsdt[i]*dsdt[j]*dsdw[k]['+str(i)+'];\n' 
            switch_code_text += '            }\n'
            switch_code_text += '        }\n'

            switch_code_text += '        break;\n'

        w += solution_calib_template.format(module=self.module, 
            number_components=c, func=self.g_list[j][0], 
            switch_code=switch_code_text, moles_assign=moles_assign_text,
            number_ordering=ns, 
            order_initial_guess=order_initial_guess_text,
            order_assign=extraCb+order_assign_text+extraCa,
            extra_ordering_code=eCB)

        #######################
        # Extra ordering code #
        # d2sdpdw             #
        #######################
        if self.verbose:
            print ("... ( 8/18) computing and writing code for d2sdpdw")
        eCB = 'static void order_d2sdpdw(double T, double P,' \
            + ' double n['+str(c)+'],' \
            + (' double b['+str(c)+'],' if SpMdClass else '') \
            + ' double s['+str(ns)+'],' \
            + ' double invd2gds2['+str(ns)+']['+str(ns)+'],' \
            + ' double d2sdpdw['+str(ns)+']['+str(np)+']) {\n'
        eCB += '    double dsdp['+str(ns)+'], dsdw['+str(ns)+']['+str(np)+'];\n'
        eCB += '    double temp['+str(ns)+'];\n'
        eCB += '    int i,j,k,l,ll;\n'
        eCB += '    double d2gdsdp['+str(ns)+'];\n'
        eCB += '    double d2gdsdw['+str(ns)+']['+str(np)+'];\n'
        eCB += '    double d3gdsdpdw['+str(ns)+']['+str(np)+'];\n'
        eCB += '    double d3gds2dp['+str(ns)+']['+str(ns)+'];\n'
        eCB += '    double d3gds2dw['+str(ns)+']['+str(ns)+']['+str(np)+'];\n'
        eCB += '    double d3gds3['+str(ns)+']['+str(ns)+']['+str(ns)+'];\n'
        eCB += moles_assign_text + '\n'
        if SpMdClass:
            for i in range(0,c):
                eCB += '    double b' + str(i+1) + ' = b[' + str(i) + '];\n'
        eCB += order_assign_text + '\n'
        for k in range(0,ns):
            f = sFunc.function[k]
            eCB += '    d2gdsdp['+str(k)+'] = '
            eCB += '0' if dmin else  self.printer.doprint(f.diff(P)) 
            eCB += ';\n'
            for l in range(0,np):
                eCB += '    d2gdsdw['+str(k)+']['+str(l) + '] = '
                eCB += '0' if dmin else self.printer.doprint(
                    f.diff(symparam[l])) 
                eCB += ';\n'
                eCB += '    d3gdsdpdw['+str(k)+']['+str(l) + '] = '
                eCB += '0' if dmin else self.printer.doprint(
                    f.diff(symparam[l]).diff(P))
                eCB += ';\n'
            for l in range(k,ns):
                eCB += '    d3gds2dp['+str(k)+']['+str(l)+'] = '
                eCB += '0' if dmin else self.printer.doprint(
                    f.diff(s[l]).diff(P))
                eCB += ';\n'
                if l > k:
                    eCB += '    d3gds2dp['+str(l)+']['+str(k)+'] = '
                    eCB += 'd3gds2dp['+str(k)+']['+str(l)+'];\n'
                for ll in range(0,np):
                    eCB += '    d3gds2dw['+str(k)+']['+str(l)+']['+str(ll)+'] = '
                    eCB += '0' if dmin else self.printer.doprint(
                        f.diff(s[l]).diff(symparam[ll]))
                    eCB += ';\n'
                    if l > k:
                        eCB += '    d3gds2dw['+str(l)+']['+str(k)+']['+str(ll) \
                              +'] = d3gds2dw['+str(k)+']['+str(l)+']['+str(ll) \
                             + '];\n'
                for ll in range(l,ns):
                    eCB += '    d3gds3['+str(k)+']['+str(l)+']['+str(ll)+'] = '
                    eCB += '0' if dmin else self.printer.doprint(
                        f.diff(s[l]).diff(s[ll])) 
                    eCB += ';\n'
                    eCB += self.exchange_indices(k,l,ll,'d3gds3')
        eCB += '    for (i=0; i<'+str(ns)+'; i++) {\n'
        eCB += '        dsdp[i] = 0.0;\n'
        eCB += '        for (j=0; j<'+str(ns)+'; j++) dsdp[i] += '
        eCB +=            '- invd2gds2[i][j]*d2gdsdp[j];\n'
        eCB += '    }\n'
        eCB += '    for (ll=0; ll<'+str(np)+'; ll++) {\n'
        eCB += '        for (i=0; i<'+str(ns)+'; i++) {\n'
        eCB += '            dsdw[i][ll] = 0.0;\n'
        eCB += '            for (j=0; j<'+str(ns)+'; j++) '
        eCB +=              'dsdw[i][ll] += - invd2gds2[i][j]*d2gdsdw[j][ll];\n'
        eCB += '        }\n'
        eCB += '    }\n'
        eCB += '    for (ll=0; ll<'+str(np)+'; ll++) {\n'
        eCB += '        for (i=0; i<'+str(ns)+'; i++) {\n'
        eCB += '            for (j=0; j<'+str(ns)+'; j++) {\n'
        eCB += '                temp[j] = d3gdsdpdw[j][ll];\n'
        eCB += '                for (k=0; k<'+str(ns)+'; k++) {\n'
        eCB += '                    temp[j] += d3gds2dp[j][k]*dsdw[k][ll] '
        eCB +=                              '+ d3gds2dw[j][k][ll]*dsdp[k];\n'
        eCB += '                    for (l=0; l<'+str(ns)+'; l++) temp[j] += '
        eCB +=                         'd3gds3[j][k][l]*dsdp[k]*dsdw[l][ll];\n'
        eCB += '                }\n'
        eCB += '            }\n'
        eCB += '            d2sdpdw[i][ll] = 0.0;\n'
        eCB += '            for (j=0; j<'+str(ns)+'; j++) '
        eCB +=              'd2sdpdw[i][ll] += - invd2gds2[i][j]*temp[j];\n'
        eCB += '        }\n'
        eCB += '    }\n'
        eCB += '}\n'

        ###############################
        # Block of code for d3gdtdpdw #
        ###############################
        if self.verbose:
            print ("... ( 9/18) computing and writing code for d3gdtdpdw")
        extraCb  = '    double d2gdsdt['+str(ns)+'];\n'
        extraCb += '    double d2gdsdp['+str(ns)+'];\n'
        extraCb += '    double d2gdsdw['+str(ns)+']['+str(np)+'];\n'
        extraCb += '    double d2gds2['+str(ns)+']['+str(ns)+'];\n'
        extraCb += '    double dsdt['+str(ns)+'];\n'
        extraCb += '    double dsdp['+str(ns)+'];\n'
        extraCb += '    double dsdw['+str(ns)+']['+str(np)+'];\n'
        extraCb += '    double d2sdtdp['+str(ns)+'];\n'
        extraCb += '    double d2sdtdw['+str(ns)+']['+str(np)+'];\n'
        extraCb += '    double d2sdpdw['+str(ns)+']['+str(np)+'];\n'
        extraCb += '    double d3gds3['+str(ns)+']['+str(ns)+']['+str(ns)+'];\n'
        extraCb += '    double d3gds2dt['+str(ns)+']['+str(ns)+'];\n'
        extraCb += '    double d3gds2dp['+str(ns)+']['+str(ns)+'];\n'
        extraCb += '    double d3gds2dw['+str(ns)+']['+str(ns)+']['+str(np)+'];\n'
        extraCb += '    double d3gdsdtdp['+str(ns)+'];\n'
        extraCb += '    double d3gdsdtdw['+str(ns)+']['+str(np)+'];\n'
        extraCb += '    double d3gdsdpdw['+str(ns)+']['+str(np)+'];\n'
        extraCb += '    int i,j,k;\n'
        extraCb += '    order_dsdt(T, P, n, s, invd2gds2, dsdt);\n'
        extraCb += '    order_dsdp(T, P, n, s, invd2gds2, dsdp);\n'
        extraCb += '    order_dsdw(T, P, n, s, invd2gds2, dsdw);\n'
        extraCb += '    order_d2sdtdp(T, P, n, s, invd2gds2, d2sdtdp);\n'
        extraCb += '    order_d2sdtdw(T, P, n, s, invd2gds2, d2sdtdw);\n'
        extraCb += '    order_d2sdpdw(T, P, n, s, invd2gds2, d2sdpdw);\n'
        extraCa  = ''
        for k in range(0,ns):
            extraCa += '    d2gdsdt['+str(k)+'] = '
            extraCa += '0' if dmin else self.printer.doprint(
                G.diff(T).diff(s[k]))
            extraCa += ';\n'
            extraCa += '    d2gdsdp['+str(k)+'] = '
            extraCa += '0' if dmin else self.printer.doprint(
                G.diff(P).diff(s[k]))
            extraCa += ';\n'
            extraCa += '    d3gdsdtdp['+str(k)+'] = '
            extraCa += '0' if dmin else self.printer.doprint(
                G.diff(T).diff(P).diff(s[k]))
            extraCa += ';\n'
            for l in range(k,ns):
                extraCa += '    d2gds2['+str(k)+']['+str(l)+'] = '
                extraCa += '0' if dmin else self.printer.doprint(
                    G.diff(s[k]).diff(s[l]))
                extraCa += ';\n'
                extraCa += '    d3gds2dt['+str(k)+']['+str(l)+'] = '
                extraCa += '0' if dmin else self.printer.doprint(
                    G.diff(s[k]).diff(s[l]).diff(T))
                extraCa += ';\n'
                extraCa += '    d3gds2dp['+str(k)+']['+str(l)+'] = '
                extraCa += '0' if dmin else self.printer.doprint(
                    G.diff(s[k]).diff(s[l]).diff(P))
                extraCa += ';\n'
                if l > k:
                    extraCa += '    d2gds2['+str(k)+']['+str(l)+'] = '
                    extraCa +=    ' d2gds2['+str(l)+']['+str(k)+'];\n'
                    extraCa += '    d3gds2dt['+str(k)+']['+str(l)+'] = '
                    extraCa +=    ' d3gds2dt['+str(l)+']['+str(k)+'];\n'
                    extraCa += '    d3gds2dp['+str(k)+']['+str(l)+'] = '
                    extraCa +=    ' d3gds2dp['+str(l)+']['+str(k)+'];\n'
                for ll in range(0,np):
                    extraCa += '    d3gds2dw['+str(k)+']['+str(l)+']['+str(ll)+'] = '
                    extraCa += '0' if dmin else self.printer.doprint(
                        G.diff(s[k]).diff(s[l]).diff(symparam[ll]))
                    extraCa += ';\n'
                    if l > k:
                        extraCa += '    d3gds2dw['+str(k)+']['+str(l)+']['+str(ll)+'] = '
                        extraCa +=    ' d3gds2dw['+str(l)+']['+str(k)+']['+str(ll)+'];\n'
                for ll in range(l,ns):
                    extraCa += '    d3gds3['+str(k)+']['+str(l)+']['+str(ll)+'] = '
                    extraCa += '0' if dmin else self.printer.doprint(
                        G.diff(s[k]).diff(s[l]).diff(s[ll]))
                    extraCa += ';\n'
                    extraCa += self.exchange_indices(k,l,ll,'d3gds3')
            for l in range(0,np):
                extraCa += '    d2gdsdw['+str(k)+']['+str(l)+'] = '
                extraCa += '0' if dmin else self.printer.doprint(
                    G.diff(s[k]).diff(symparam[l]))
                extraCa += ';\n'
                extraCa += '    d3gdsdtdw['+str(k)+']['+str(l)+'] = '
                extraCa += '0' if dmin else self.printer.doprint(
                    G.diff(s[k]).diff(symparam[l]).diff(T))
                extraCa += ';\n'
                extraCa += '    d3gdsdpdw['+str(k)+']['+str(l)+'] = '
                extraCa += '0' if dmin else self.printer.doprint(
                    G.diff(s[k]).diff(symparam[l]).diff(P))
                extraCa += ';\n'

        j += 1
        if dmin:
            G_jac_list = ['0' for i in range(0, np)]
        else:
            G_jac_list = [self.printer.doprint(
                G_param_jac[j,i]) for i in range(0, np)]

        switch_code_text = ''
        for i in range(0, len(self.params)):
            switch_code_text += '    case ' + str(i) + ':\n'
            switch_code_text += '        result += ' + G_jac_list[i] + ';\n'

            switch_code_text += '        for (i=0; i<'+str(ns)+'; i++) {\n'
            switch_code_text += '            result += ' \
            'd3gdsdtdp[i]*dsdw[i]['+str(i)+'] ' \
            + '+ d2gdsdt[i]*d2sdpdw[i]['+str(i)+'] ' \
            + '+ d2gdsdp[i]*d2sdtdw[i]['+str(i)+'] ' \
            + '+ d2gdsdw[i]['+str(i)+']*d2sdtdp[i] ' \
            + '+ d3gdsdtdw[i]['+str(i)+']*dsdp[i] ' \
            + '+ d3gdsdpdw[i]['+str(i)+']*dsdt[i];\n'
            switch_code_text += '            for (j=0; j<'+str(ns)+'; j++) {\n'
            switch_code_text += '                result += ' \
            + 'd3gds2dt[i][j]*dsdp[i]*dsdw[j]['+str(i)+']' \
            + '+ d3gds2dp[i][j]*dsdt[i]*dsdw[j]['+str(i)+']' \
            + '+ d2gds2[i][j]*d2sdtdp[i]*dsdw[j]['+str(i)+']' \
            + '+ d2gds2[i][j]*dsdp[i]*d2sdtdw[j]['+str(i)+']' \
            + '+ d2gds2[i][j]*dsdt[i]*d2sdpdw[j]['+str(i)+']' \
            + '+ d3gds2dw[i][j]['+str(i)+']*dsdt[i]*dsdp[j];\n'
            switch_code_text += '                for (k=0; k<1; k++) ' \
            'result += d3gds3[i][j][k]*dsdt[i]*dsdp[j]*dsdw[k]['+str(i)+'];\n' 
            switch_code_text += '            }\n'
            switch_code_text += '        }\n'

            switch_code_text += '        break;\n'

        w += solution_calib_template.format(module=self.module, 
            number_components=c, func=self.g_list[j][0], 
            switch_code=switch_code_text, moles_assign=moles_assign_text,
            number_ordering=ns, 
            order_initial_guess=order_initial_guess_text,
            order_assign=extraCb+order_assign_text+extraCa,
            extra_ordering_code=eCB)

        ##############################
        # Block of code for d3gdp2dw #
        ##############################
        if self.verbose:
            print ("... (10/18) computing and writing code for d3gdp2dw")
        extraCb  = '    double d2gdsdp['+str(ns)+'];\n'
        extraCb += '    double d2gdsdw['+str(ns)+']['+str(np)+'];\n'
        extraCb += '    double d2gds2['+str(ns)+']['+str(ns)+'];\n'
        extraCb += '    double dsdp['+str(ns)+'];\n'
        extraCb += '    double dsdw['+str(ns)+']['+str(np)+'];\n'
        extraCb += '    double d2sdp2['+str(ns)+'];\n'
        extraCb += '    double d2sdpdw['+str(ns)+']['+str(np)+'];\n'
        extraCb += '    double d3gds3['+str(ns)+']['+str(ns)+']['+str(ns)+'];\n'
        extraCb += '    double d3gds2dp['+str(ns)+']['+str(ns)+'];\n'
        extraCb += '    double d3gds2dw['+str(ns)+']['+str(ns)+']['+str(np)+'];\n'
        extraCb += '    double d3gdsdp2['+str(ns)+'];\n'
        extraCb += '    double d3gdsdpdw['+str(ns)+']['+str(np)+'];\n'
        extraCb += '    int i,j,k;\n'
        extraCb += '    order_dsdp(T, P, n, s, invd2gds2, dsdp);\n'
        extraCb += '    order_dsdw(T, P, n, s, invd2gds2, dsdw);\n'
        extraCb += '    order_d2sdp2(T, P, n, s, invd2gds2, d2sdp2);\n'
        extraCb += '    order_d2sdpdw(T, P, n, s, invd2gds2, d2sdpdw);\n'
        extraCa  = ''
        for k in range(0,ns):
            extraCa += '    d2gdsdp['+str(k)+'] = '
            extraCa += '0' if dmin else self.printer.doprint(
                G.diff(P).diff(s[k]))
            extraCa += ';\n'
            extraCa += '    d3gdsdp2['+str(k)+'] = '
            extraCa += '0' if dmin else self.printer.doprint(
                G.diff(P,2).diff(s[k]))
            extraCa += ';\n'
            for l in range(k,ns):
                extraCa += '    d2gds2['+str(k)+']['+str(l)+'] = '
                extraCa += '0' if dmin else self.printer.doprint(
                    G.diff(s[k]).diff(s[l]))
                extraCa += ';\n'
                extraCa += '    d3gds2dp['+str(k)+']['+str(l)+'] = '
                extraCa += '0' if dmin else self.printer.doprint(
                    G.diff(s[k]).diff(s[l]).diff(P))
                extraCa += ';\n'
                if l > k:
                    extraCa += '    d2gds2['+str(k)+']['+str(l)+'] = '
                    extraCa +=    ' d2gds2['+str(l)+']['+str(k)+'];\n'
                    extraCa += '    d3gds2dp['+str(k)+']['+str(l)+'] = '
                    extraCa +=    ' d3gds2dp['+str(l)+']['+str(k)+'];\n'
                for ll in range(0,np):
                    extraCa += '    d3gds2dw['+str(k)+']['+str(l)+']['+str(ll)+'] = '
                    extraCa += '0' if dmin else self.printer.doprint(
                        G.diff(s[k]).diff(s[l]).diff(symparam[ll]))
                    extraCa += ';\n'
                    if l > k:
                        extraCa += '    d3gds2dw['+str(k)+']['+str(l)+']['+str(ll)+'] = '
                        extraCa +=    ' d3gds2dw['+str(l)+']['+str(k)+']['+str(ll)+'];\n'
                for ll in range(l,ns):
                    extraCa += '    d3gds3['+str(k)+']['+str(l)+']['+str(ll)+'] = '
                    extraCa += '0' if dmin else self.printer.doprint(
                        G.diff(s[k]).diff(s[l]).diff(s[ll]))
                    extraCa += ';\n'
                    extraCa += self.exchange_indices(k,l,ll,'d3gds3')
            for l in range(0,np):
                extraCa += '    d2gdsdw['+str(k)+']['+str(l)+'] = '
                extraCa += '0' if dmin else self.printer.doprint(
                    G.diff(s[k]).diff(symparam[l]))
                extraCa += ';\n'
                extraCa += '    d3gdsdpdw['+str(k)+']['+str(l)+'] = '
                extraCa += '0' if dmin else self.printer.doprint(
                    G.diff(s[k]).diff(symparam[l]).diff(P))
                extraCa += ';\n'

        j += 1
        if dmin:
            G_jac_list = ['0' for i in range(0, np)]
        else:
            G_jac_list = [self.printer.doprint(
                G_param_jac[j,i]) for i in range(0, np)]

        switch_code_text = ''
        for i in range(0, len(self.params)):
            switch_code_text += '    case ' + str(i) + ':\n'
            switch_code_text += '        result += ' + G_jac_list[i] + ';\n'

            switch_code_text += '        for (i=0; i<'+str(ns)+'; i++) {\n'
            switch_code_text += '            result += ' \
            'd3gdsdp2[i]*dsdw[i]['+str(i)+'] ' \
            + '+ 2.0*d2gdsdp[i]*d2sdpdw[i]['+str(i)+'] ' \
            + '+ d2gdsdw[i]['+str(i)+']*d2sdp2[i] ' \
            + '+ 2.0*d3gdsdpdw[i]['+str(i)+']*dsdp[i];\n'
            switch_code_text += '            for (j=0; j<'+str(ns)+'; j++) {\n'
            switch_code_text += '                result += ' \
            + '2.0*d3gds2dp[i][j]*dsdp[i]*dsdw[j]['+str(i)+']' \
            + '+ d2gds2[i][j]*d2sdp2[i]*dsdw[j]['+str(i)+']' \
            + '+ 2.0*d2gds2[i][j]*dsdp[i]*d2sdpdw[j]['+str(i)+']' \
            + '+ d3gds2dw[i][j]['+str(i)+']*dsdp[i]*dsdp[j];\n'
            switch_code_text += '                for (k=0; k<1; k++) ' \
            'result += d3gds3[i][j][k]*dsdp[i]*dsdp[j]*dsdw[k]['+str(i)+'];\n' 
            switch_code_text += '            }\n'
            switch_code_text += '        }\n'

            switch_code_text += '        break;\n'

        w += solution_calib_template.format(module=self.module, 
            number_components=c, func=self.g_list[j][0], 
            switch_code=switch_code_text, moles_assign=moles_assign_text,
            number_ordering=ns, 
            order_initial_guess=order_initial_guess_text,
            order_assign=extraCb+order_assign_text+extraCa,
            extra_ordering_code='')

        #######################
        # Extra ordering code #
        # d3sdt2dw, d3sdt3    #
        #######################
        if self.verbose:
            print ("... (11/18) computing and writing code for d3sdt2dw, d3sdt3")
        eCB = 'static void order_d3sdt2dw(double T, double P,' \
            + ' double n['+str(c)+'],' \
            + (' double b['+str(c)+'],' if SpMdClass else '') \
            + ' double s['+str(ns)+'],' \
            + ' double invd2gds2['+str(ns)+']['+str(ns)+'],' \
            + ' double d3sdt2dw['+str(ns)+']['+str(np)+']) {\n'
        eCB += '    double dsdt['+str(ns)+'];\n'
        eCB += '    double dsdw['+str(ns)+']['+str(np)+'];\n'
        eCB += '    double d2sdt2['+str(ns)+'];\n'
        eCB += '    double d2sdtdw['+str(ns)+']['+str(np)+'];\n'
        eCB += '    order_dsdt(T, P, n, s, invd2gds2, dsdt);\n'
        eCB += '    order_dsdw(T, P, n, s, invd2gds2, dsdw);\n'
        eCB += '    order_d2sdt2(T, P, n, s, invd2gds2, d2sdt2);\n'
        eCB += '    order_d2sdtdw(T, P, n, s, invd2gds2, d2sdtdw);\n'
        eCB += '    int i,j,k,l,m,ll;\n'
        eCB += '    double d4gdsdt2dw['+str(ns)+']['+str(np)+'];\n'
        eCB += '    double d4gds2dtdw['+str(ns)+']['+str(ns)+']['+str(np)+'];\n'
        eCB += '    double d4gds2dt2['+str(ns)+']['+str(ns)+'];\n'
        eCB += '    double d4gds3dt['+str(ns)+']['+str(ns)+']['+str(ns)+'];\n'
        eCB += '    double d4gds3dw['+str(ns)+']['+str(ns)+']['+str(ns)+']['+str(np)+'];\n'
        eCB += '    double d4gds4['+str(ns)+']['+str(ns)+']['+str(ns)+']['+str(ns)+'];\n'
        eCB += '    double d3gds2dt['+str(ns)+']['+str(ns)+'];\n'
        eCB += '    double d3gds2dw['+str(ns)+']['+str(ns)+']['+str(np)+'];\n'
        eCB += '    double d3gds3['+str(ns)+']['+str(ns)+']['+str(ns)+'];\n'
        eCB += '    double temp['+str(ns)+'];\n'
        eCB += moles_assign_text + '\n'
        if SpMdClass:
            for i in range(0,c):
                eCB += '    double b' + str(i+1) + ' = b[' + str(i) + '];\n'
        eCB += order_assign_text + '\n'
        for k in range(0,ns):
            f = sFunc.function[k]
            for l in range(0,np):
                eCB += '    d4gdsdt2dw['+str(k)+']['+str(l)+'] = '
                eCB += '0' if dmin else self.printer.doprint(
                        f.diff(symparam[l]).diff(T,2))
                eCB += ';\n'
            for l in range(0,ns):
                eCB += '    d3gds2dt['+str(k)+']['+str(l)+'] = '
                eCB += '0' if dmin else self.printer.doprint(
                    f.diff(s[l]).diff(T))
                eCB += ';\n'
                eCB += '    d4gds2dt2['+str(k)+']['+str(l)+'] = '
                eCB += '0' if dmin else self.printer.doprint(
                    f.diff(s[l]).diff(T,2))
                eCB += ';\n'
                for kk in range(0,np):
                    eCB += '    d3gds2dw['+str(k)+']['+str(l)+']['+str(kk)+'] = '
                    eCB += '0' if dmin else self.printer.doprint(
                            f.diff(s[l]).diff(symparam[kk]))
                    eCB += ';\n'
                    eCB += '    d4gds2dtdw['+str(k)+']['+str(l)+']['+str(kk)+'] = '
                    eCB += '0' if dmin else self.printer.doprint(
                            f.diff(s[l]).diff(symparam[kk]).diff(T))
                    eCB += ';\n'
                for kk in range(0,ns):
                    eCB += '    d3gds3['+str(k)+']['+str(l)+']['+str(kk)+'] = '
                    eCB += '0' if dmin else self.printer.doprint(
                        f.diff(s[l]).diff(s[kk]))
                    eCB += ';\n'
                    eCB += '    d4gds3dt['+str(k)+']['+str(l)+']['+str(kk)+'] = '
                    eCB += '0' if dmin else self.printer.doprint(
                            f.diff(s[l]).diff(s[kk]).diff(T))
                    eCB += ';\n'
                    for ll in range(0,np):
                        eCB += '    d4gds3dw['+str(k)+']['+str(l)+'][' \
                        +str(kk)+']['+str(ll)+'] = '
                        eCB += '0' if dmin else self.printer.doprint(
                            f.diff(s[l]).diff(s[kk]).diff(symparam[ll]))
                        eCB += ';\n'
                    for ll in range(0,ns):
                        eCB += '    d4gds4['+str(k)+']['+str(l)+'][' \
                        +str(kk)+']['+str(ll)+'] = '
                        eCB += '0' if dmin else self.printer.doprint(
                            f.diff(s[l]).diff(s[kk]).diff(s[ll]))
                        eCB += ';\n'
        eCB += '    for (ll=0; ll<'+str(np)+'; ll++) {\n'
        eCB += '        for (i=0; i<'+str(ns)+'; i++) {\n'
        eCB += '            for (j=0; j<'+str(ns)+'; j++) {\n'
        eCB += '                temp[j] = d4gdsdt2dw[j][ll];\n'
        eCB += '                for (k=0; k<'+str(ns)+'; k++) {\n'
        eCB += '                    temp[j] += d4gds2dt2[j][k]*dsdw[k][ll];\n'
        eCB += '                    temp[j] += 2.0*d4gds2dtdw[j][k][ll]*dsdt[k];\n'
        eCB += '                    temp[j] += 2.0*d3gds2dt[j][k]*d2sdtdw[k][ll];\n'
        eCB += '                    temp[j] += d3gds2dw[j][k][ll]*d2sdt2[k];\n'
        eCB += '                    for (l=0; l<'+str(ns)+'; l++) {\n' 
        eCB += '                        temp[j] += 2.0*d4gds3dt[j][k][l]*dsdt[k]*dsdw[l][ll];\n'
        eCB += '                        temp[j] += d4gds3dw[j][k][l][ll]*dsdt[k]*dsdt[l];\n'
        eCB += '                        temp[j] += 2.0*d3gds3[j][k][l]*d2sdtdw[k][ll]'
        eCB +=                          '*dsdt[l];\n'
        eCB += '                        temp[j] += d3gds3[j][k][l]*d2sdt2[k]'
        eCB +=                          '*dsdw[l][ll];\n'
        eCB += '                        for (m=0; m<'+str(ns)+'; m++) {\n'
        eCB += '                            temp[j] += d4gds4[j][k][l][m]*'
        eCB +=                              'dsdt[k]*dsdt[l]*dsdw[m][ll];\n'
        eCB += '                        }\n'
        eCB += '                    }\n'
        eCB += '                }\n'
        eCB += '            }\n'
        eCB += '            d3sdt2dw[i][ll] = 0.0;\n'
        eCB += '            for (j=0; j<'+str(ns)+'; j++) '
        eCB +=                'd3sdt2dw[i][ll] += - invd2gds2[i][j]*temp[j];\n'
        eCB += '        }\n'
        eCB += '    }\n'
        eCB += '}\n'

        eCB += 'static void order_d3sdt3(double T, double P,' \
            +  ' double n['+str(c)+'],' \
            +  (' double b['+str(c)+'],' if SpMdClass else '') \
            +  ' double s['+str(ns)+'],' \
            +  ' double invd2gds2['+str(ns)+']['+str(ns)+'],' \
            +  ' double d3sdt3['+str(ns)+']) {\n'
        eCB += '    double dsdt['+str(ns)+'];\n'
        eCB += '    double d2sdt2['+str(ns)+'];\n'
        eCB += '    order_dsdt(T, P, n, s, invd2gds2, dsdt);\n'
        eCB += '    order_d2sdt2(T, P, n, s, invd2gds2, d2sdt2);\n'
        eCB += '    int i,j,k,l,m;\n'
        eCB += '    double d4gdsdt3['+str(ns)+'];\n'
        eCB += '    double d4gds2dt2['+str(ns)+']['+str(ns)+'];\n'
        eCB += '    double d4gds3dt['+str(ns)+']['+str(ns)+']['+str(ns)+'];\n'
        eCB += '    double d4gds4['+str(ns)+']['+str(ns)+']['+str(ns)+']['+str(ns)+'];\n'
        eCB += '    double d3gds2dt['+str(ns)+']['+str(ns)+'];\n'
        eCB += '    double d3gds3['+str(ns)+']['+str(ns)+']['+str(ns)+'];\n'
        eCB += '    double temp['+str(ns)+'];\n'
        eCB += moles_assign_text + '\n'
        if SpMdClass:
            for i in range(0,c):
                eCB += '    double b' + str(i+1) + ' = b[' + str(i) + '];\n'
        eCB += order_assign_text + '\n'
        for k in range(0,ns):
            f = sFunc.function[k]
            for l in range(0,ns):
                eCB += '    d3gds2dt['+str(k)+']['+str(l)+'] = '
                eCB += '0' if dmin else self.printer.doprint(
                    f.diff(s[l]).diff(T))
                eCB += ';\n'
                eCB += '    d4gds2dt2['+str(k)+']['+str(l)+'] = '
                eCB += '0' if dmin else self.printer.doprint(
                    f.diff(s[l]).diff(T,2))
                eCB += ';\n'
                for kk in range(0,ns):
                    eCB += '    d3gds3['+str(k)+']['+str(l)+']['+str(kk)+'] = '
                    eCB += '0' if dmin else self.printer.doprint(
                        f.diff(s[l]).diff(s[kk]))
                    eCB += ';\n'
                    eCB += '    d4gds3dt['+str(k)+']['+str(l)+']['+str(kk)+'] = '
                    eCB += '0' if dmin else self.printer.doprint(
                            f.diff(s[l]).diff(s[kk]).diff(T))
                    eCB += ';\n'
                    for ll in range(0,ns):
                        eCB += '    d4gds4['+str(k)+']['+str(l)+'][' \
                        +str(kk)+']['+str(ll)+'] = '
                        eCB += '0' if dmin else self.printer.doprint(
                            f.diff(s[l]).diff(s[kk]).diff(s[ll]))
                        eCB += ';\n'
        eCB += '    for (i=0; i<'+str(ns)+'; i++) {\n'
        eCB += '        for (j=0; j<'+str(ns)+'; j++) {\n'
        eCB += '            temp[j] = d4gdsdt3[j];\n'
        eCB += '            for (k=0; k<'+str(ns)+'; k++) {\n'
        eCB += '                temp[j] += 3.0*d4gds2dt2[j][k]*dsdt[k];\n'
        eCB += '                temp[j] += 3.0*d3gds2dt[j][k]*d2sdt2[k];\n'
        eCB += '                for (l=0; l<'+str(ns)+'; l++) {\n' 
        eCB += '                    temp[j] += 3.0*d4gds3dt[j][k][l]*dsdt[k]*dsdt[l];\n'
        eCB += '                    temp[j] += 3.0*d3gds3[j][k][l]*d2sdt2[k]*dsdt[l];\n'
        eCB += '                    for (m=0; m<'+str(ns)+'; m++) {\n'
        eCB += '                        temp[j] += d4gds4[j][k][l][m]*'
        eCB +=                          'dsdt[k]*dsdt[l]*dsdt[m];\n'
        eCB += '                    }\n'
        eCB += '                }\n'
        eCB += '            }\n'
        eCB += '        }\n'
        eCB += '        d3sdt3[i] = 0.0;\n'
        eCB += '        for (j=0; j<'+str(ns)+'; j++) '
        eCB +=            'd3sdt3[i] += - invd2gds2[i][j]*temp[j];\n'
        eCB += '    }\n'
        eCB += '}\n'

        ##############################
        # Block of code for d4gdt3dw #
        ##############################
        if self.verbose:
            print ("... (12/18) computing and writing code for d4gdt3dw")
        extraCb  = '    double d4gdsdt3['+str(ns)+'];\n'
        extraCb += '    double d4gds2dt2['+str(ns)+']['+str(ns)+'];\n'
        extraCb += '    double d4gdsdt2dw['+str(ns)+']['+str(np)+'];\n'
        extraCb += '    double d4gds3dt['+str(ns)+']['+str(ns)+']['+str(ns)+'];\n'
        extraCb += '    double d4gds2dtdw['+str(ns)+']['+str(ns)+']['+str(np)+'];\n'
        extraCb += '    double d4gds3dw['+str(ns)+']['+str(ns)+']['+str(ns)+']['+str(np)+'];\n'
        extraCb += '    double d4gds4['+str(ns)+']['+str(ns)+']['+str(ns)+']['+str(ns)+'];\n'
        extraCb += '    double d3gdsdt2['+str(ns)+'];\n'
        extraCb += '    double d3gds2dt['+str(ns)+']['+str(ns)+'];\n'
        extraCb += '    double d3gdsdtdw['+str(ns)+']['+str(np)+'];\n'
        extraCb += '    double d3gds2dw['+str(ns)+']['+str(ns)+']['+str(np)+'];\n'
        extraCb += '    double d3gds3['+str(ns)+']['+str(ns)+']['+str(ns)+'];\n'
        extraCb += '    double d2gdsdt['+str(ns)+'];\n'
        extraCb += '    double d2gdsdw['+str(ns)+']['+str(np)+'];\n'
        extraCb += '    double d2gds2['+str(ns)+']['+str(ns)+'];\n'
        extraCb += '    double dsdt['+str(ns)+'];\n'
        extraCb += '    double dsdw['+str(ns)+']['+str(np)+'];\n'
        extraCb += '    double d2sdt2['+str(ns)+'];\n'
        extraCb += '    double d2sdtdw['+str(ns)+']['+str(np)+'];\n'
        extraCb += '    double d3sdt3['+str(ns)+'];\n'
        extraCb += '    double d3sdt2dw['+str(ns)+']['+str(np)+'];\n'
        extraCb += '    int i,j,k,l;\n'
        extraCb += '    order_dsdt(T, P, n, s, invd2gds2, dsdt);\n'
        extraCb += '    order_dsdw(T, P, n, s, invd2gds2, dsdw);\n'
        extraCb += '    order_d2sdt2(T, P, n, s, invd2gds2, d2sdt2);\n'
        extraCb += '    order_d2sdtdw(T, P, n, s, invd2gds2, d2sdtdw);\n'
        extraCb += '    order_d3sdt3(T, P, n, s, invd2gds2, d3sdt3);\n'
        extraCb += '    order_d3sdt2dw(T, P, n, s, invd2gds2, d3sdt2dw);\n'

        for k in range(0,ns):
            extraCa  = '    d4gdsdt3['+str(k)+'] = '
            extraCa += '0' if dmin else self.printer.doprint(
                G.diff(T,3).diff(s[k]))
            extraCa += ';\n'
            extraCa += '    d3gdsdt2['+str(k)+'] = '
            extraCa += '0' if dmin else self.printer.doprint(
                G.diff(T,2).diff(s[k]))
            extraCa += ';\n'
            extraCa += '    d2gdsdt['+str(k)+'] = '
            extraCa += '0' if dmin else self.printer.doprint(
                G.diff(T).diff(s[k]))
            extraCa += ';\n'
            for l in range(0,np):
                extraCa += '    d4gdsdt2dw['+str(k)+']['+str(l)+'] = '
                extraCa += '0' if dmin else self.printer.doprint(
                    G.diff(T,2).diff(s[k]).diff(symparam[l]))
                extraCa += ';\n'
                extraCa += '    d3gdsdtdw['+str(k)+']['+str(l)+'] = '
                extraCa += '0' if dmin else self.printer.doprint(
                    G.diff(T).diff(s[k]).diff(symparam[l]))
                extraCa += ';\n'
                extraCa += '    d2gdsdw['+str(k)+']['+str(l)+'] = '
                extraCa += '0' if dmin else self.printer.doprint(
                    G.diff(s[k]).diff(symparam[l]))
                extraCa += ';\n'
            for l in range(0,ns):
                extraCa += '    d4gds2dt2['+str(k)+']['+str(l)+'] = '
                extraCa += '0' if dmin else self.printer.doprint(
                    G.diff(T,2).diff(s[k]).diff(s[l]))
                extraCa += ';\n'
                extraCa += '    d3gds2dt['+str(k)+']['+str(l)+'] = '
                extraCa += '0' if dmin else self.printer.doprint(
                    G.diff(T).diff(s[k]).diff(s[l]))
                extraCa += ';\n'
                extraCa += '    d2gds2['+str(k)+']['+str(l)+'] = '
                extraCa += '0' if dmin else self.printer.doprint(
                    G.diff(s[k]).diff(s[l]))
                extraCa += ';\n'
                for kk in range(0,ns):
                    extraCa += '    d4gds3dt['+str(k)+']['+str(l)+']['+str(kk)+'] = '
                    extraCa += '0' if dmin else self.printer.doprint(
                        G.diff(T).diff(s[k]).diff(s[l]).diff(s[kk]))
                    extraCa += ';\n'
                    extraCa += '    d3gds3['+str(k)+']['+str(l)+']['+str(kk)+'] = '
                    extraCa += '0' if dmin else self.printer.doprint(
                        G.diff(s[k]).diff(s[l]).diff(s[kk]))
                    extraCa += ';\n'
                    for ll in range(0,ns):
                        extraCa += '    d4gds4['+str(k)+']['+str(l)+'][' \
                                 +str(kk)+']['+str(ll)+'] = '
                        extraCa += '0' if dmin else self.printer.doprint(
                            G.diff(s[k]).diff(s[l]).diff(s[kk]).diff(s[ll]))
                        extraCa += ';\n'
                    for ll in range(0,np):
                        extraCa += '    d4gds3dw['+str(k)+']['+str(l)+'][' \
                                 +str(kk)+']['+str(ll)+'] = '
                        extraCa += '0' if dmin else self.printer.doprint(
                            G.diff(s[k]).diff(s[l]).diff(s[kk]).diff(symparam[ll]))
                        extraCa += ';\n'
                for kk in range(0,np):
                    extraCa += '    d4gds2dtdw['+str(k)+']['+str(l)+']['+str(kk)+'] = '
                    extraCa += '0' if dmin else self.printer.doprint(
                        G.diff(T).diff(s[k]).diff(s[l]).diff(symparam[kk]))
                    extraCa += ';\n'
                    extraCa += '    d3gds2dw['+str(k)+']['+str(l)+']['+str(kk)+'] = '
                    extraCa += '0' if dmin else self.printer.doprint(
                        G.diff(s[k]).diff(s[l]).diff(symparam[kk]))
                    extraCa += ';\n'

        j += 1
        if dmin:
            G_jac_list = ['0' for i in range(0, np)]
        else:
            G_jac_list = [self.printer.doprint(
                G_param_jac[j,i]) for i in range(0, np)]

        switch_code_text = ''
        for i in range(0, len(self.params)):
            switch_code_text += '    case ' + str(i) + ':\n'
            switch_code_text += '        result += ' + G_jac_list[i] + ';\n'

            switch_code_text += '        for (i=0; i<'+str(ns)+'; i++) {\n'
            switch_code_text += '            result += ' \
            + '3.0*d4gdsdt2dw[i]['+str(i)+']*dsdt[i] + ' \
            + '    d4gdsdt3[i]*dsdw[i]['+str(i)+'] + ' \
            + '3.0*d3gdsdt2[i]*d2sdtdw[i]['+str(i)+'] + ' \
            + '3.0*d3gdsdtdw[i]['+str(i)+']*d2sdt2[i] + ' \
            + '3.0*d2gdsdt[i]*d3sdt2dw[i]['+str(i)+'] + ' \
            + '    d2gdsdw[i]['+str(i)+']*d3sdt3[i];\n'
            switch_code_text += '            for (j=0; j<'+str(ns)+'; j++) {\n'
            switch_code_text += '                result += ' \
            + '3.0*d4gds2dt2[i][j]*dsdt[i]*dsdw[j]['+str(i)+'] + ' \
            + '3.0*d4gds2dtdw[i][j]['+str(i)+']*dsdt[i]*dsdt[j] + ' \
            + '6.0*d3gds2dt[i][j]*d2sdtdw[i]['+str(i)+']*dsdt[j] + ' \
            + '3.0*d3gds2dt[i][j]*d2sdt2[i]*dsdw[j]['+str(i)+'] + ' \
            + '3.0*d3gds2dw[i][j]['+str(i)+']*d2sdt2[i]*dsdt[j] + ' \
            + '3.0*d2gds2[i][j]*d2sdt2[i]*d2sdtdw[j]['+str(i)+'];\n'
            switch_code_text += '                for (k=0; k<'+str(ns)+'; k++) {\n'
            switch_code_text += '                    result += ' \
            + '3.0*d4gds3dt[i][j][k]*dsdt[i]*dsdt[j]*dsdw[k]['+str(i)+'] + ' \
            + '    d4gds3dw[i][j][k]['+str(i)+']*dsdt[i]*dsdt[j]*dsdt[k] + ' \
            + '3.0*d3gds3[i][j][k]*d2sdt2[i]*dsdt[j]*dsdw[k]['+str(i)+'] + ' \
            + '3.0*d3gds3[i][j][k]*d2sdtdw[i]['+str(i)+']*dsdt[j]*dsdt[k];\n'
            switch_code_text += '                    for (l=0; l<'+str(ns)+'; l++) ' \
            + 'result += d4gds4[i][j][k][l]*dsdt[i]*dsdt[j]*dsdt[k]*dsdw[l][' \
            + str(i)+'];\n'
            switch_code_text += '                }\n'
            switch_code_text += '            }\n'
            switch_code_text += '        }\n'

            switch_code_text += '        break;\n'

        w += solution_calib_template.format(module=self.module, 
            number_components=c, func=self.g_list[j][0], 
            switch_code=switch_code_text, moles_assign=moles_assign_text,
            number_ordering=ns, 
            order_initial_guess=order_initial_guess_text,
            order_assign=extraCb+order_assign_text+extraCa,
            extra_ordering_code=eCB)

        #######################
        # Extra ordering code #
        # d3sdtdpdw, d3sdt2dp #
        #######################
        if self.verbose:
            print ("... (13/18) computing and writing code for d3sdtdpdw, d3sdt2dp")
        eCB = 'static void order_d3sdtdpdw(double T, double P,' \
            + ' double n['+str(c)+'],' \
            + (' double b['+str(c)+'],' if SpMdClass else '') \
            + ' double s['+str(ns)+'],' \
            + ' double invd2gds2['+str(ns)+']['+str(ns)+'],' \
            + ' double d3sdtdpdw['+str(ns)+']['+str(np)+']) {\n'
        eCB += '    double dsdt['+str(ns)+'];\n'
        eCB += '    double dsdp['+str(ns)+'];\n'
        eCB += '    double dsdw['+str(ns)+']['+str(np)+'];\n'
        eCB += '    double d2sdtdp['+str(ns)+'];\n'
        eCB += '    double d2sdtdw['+str(ns)+']['+str(np)+'];\n'
        eCB += '    double d2sdpdw['+str(ns)+']['+str(np)+'];\n'
        eCB += '    order_dsdt(T, P, n, s, invd2gds2, dsdt);\n'
        eCB += '    order_dsdp(T, P, n, s, invd2gds2, dsdp);\n'
        eCB += '    order_dsdw(T, P, n, s, invd2gds2, dsdw);\n'
        eCB += '    order_d2sdtdp(T, P, n, s, invd2gds2, d2sdtdp);\n'
        eCB += '    order_d2sdtdw(T, P, n, s, invd2gds2, d2sdtdw);\n'
        eCB += '    order_d2sdpdw(T, P, n, s, invd2gds2, d2sdpdw);\n'
        eCB += '    int i,j,k,l,m,ll;\n'
        eCB += '    double d4gdsdtdpdw['+str(ns)+']['+str(np)+'];\n'
        eCB += '    double d4gds2dtdw['+str(ns)+']['+str(ns)+']['+str(np)+'];\n'
        eCB += '    double d4gds2dpdw['+str(ns)+']['+str(ns)+']['+str(np)+'];\n'
        eCB += '    double d4gds2dtdp['+str(ns)+']['+str(ns)+'];\n'
        eCB += '    double d4gds3dt['+str(ns)+']['+str(ns)+']['+str(ns)+'];\n'
        eCB += '    double d4gds3dp['+str(ns)+']['+str(ns)+']['+str(ns)+'];\n'
        eCB += '    double d4gds3dw['+str(ns)+']['+str(ns)+']['+str(ns)+']['+str(np)+'];\n'
        eCB += '    double d4gds4['+str(ns)+']['+str(ns)+']['+str(ns)+']['+str(ns)+'];\n'
        eCB += '    double d3gds2dt['+str(ns)+']['+str(ns)+'];\n'
        eCB += '    double d3gds2dp['+str(ns)+']['+str(ns)+'];\n'
        eCB += '    double d3gds2dw['+str(ns)+']['+str(ns)+']['+str(np)+'];\n'
        eCB += '    double d3gds3['+str(ns)+']['+str(ns)+']['+str(ns)+'];\n'
        eCB += '    double temp['+str(ns)+'];\n'
        eCB += moles_assign_text + '\n'
        if SpMdClass:
            for i in range(0,c):
                eCB += '    double b' + str(i+1) + ' = b[' + str(i) + '];\n'
        eCB += order_assign_text + '\n'
        for k in range(0,ns):
            f = sFunc.function[k]
            for l in range(0,np):
                eCB += '    d4gdsdtdpdw['+str(k)+']['+str(l)+'] = '
                eCB += '0' if dmin else self.printer.doprint(
                        f.diff(symparam[l]).diff(T).diff(P))
                eCB += ';\n'
            for l in range(0,ns):
                eCB += '    d3gds2dt['+str(k)+']['+str(l)+'] = '
                eCB += '0' if dmin else self.printer.doprint(
                    f.diff(s[l]).diff(T))
                eCB += ';\n'
                eCB += '    d3gds2dp['+str(k)+']['+str(l)+'] = '
                eCB += '0' if dmin else self.printer.doprint(
                    f.diff(s[l]).diff(P))
                eCB += ';\n'
                eCB += '    d4gds2dtdp['+str(k)+']['+str(l)+'] = '
                eCB += '0' if dmin else self.printer.doprint(
                    f.diff(s[l]).diff(T).diff(P))
                eCB += ';\n'
                for kk in range(0,np):
                    eCB += '    d3gds2dw['+str(k)+']['+str(l)+']['+str(kk)+'] = '
                    eCB += '0' if dmin else self.printer.doprint(
                            f.diff(s[l]).diff(symparam[kk]))
                    eCB += ';\n'
                    eCB += '    d4gds2dtdw['+str(k)+']['+str(l)+']['+str(kk)+'] = '
                    eCB += '0' if dmin else self.printer.doprint(
                            f.diff(s[l]).diff(symparam[kk]).diff(T))
                    eCB += ';\n'
                    eCB += '    d4gds2dpdw['+str(k)+']['+str(l)+']['+str(kk)+'] = '
                    eCB += '0' if dmin else self.printer.doprint(
                            f.diff(s[l]).diff(symparam[kk]).diff(P))
                    eCB += ';\n'
                for kk in range(0,ns):
                    eCB += '    d3gds3['+str(k)+']['+str(l)+']['+str(kk)+'] = '
                    eCB += '0' if dmin else self.printer.doprint(
                        f.diff(s[l]).diff(s[kk]))
                    eCB += ';\n'
                    eCB += '    d4gds3dt['+str(k)+']['+str(l)+']['+str(kk)+'] = '
                    eCB += '0' if dmin else self.printer.doprint(
                            f.diff(s[l]).diff(s[kk]).diff(T))
                    eCB += ';\n'
                    eCB += '    d4gds3dp['+str(k)+']['+str(l)+']['+str(kk)+'] = '
                    eCB += '0' if dmin else self.printer.doprint(
                            f.diff(s[l]).diff(s[kk]).diff(P))
                    eCB += ';\n'
                    for ll in range(0,np):
                        eCB += '    d4gds3dw['+str(k)+']['+str(l)+'][' \
                             +str(kk)+']['+str(ll)+'] = ' 
                        eCB += '0' if dmin else self.printer.doprint(
                            f.diff(s[l]).diff(s[kk]).diff(symparam[ll]))
                        eCB += ';\n'
                    for ll in range(0,ns):
                        eCB += '    d4gds4['+str(k)+']['+str(l)+'][' \
                        +str(kk)+']['+str(ll)+'] = ' 
                        eCB += '0' if dmin else self.printer.doprint(
                            f.diff(s[l]).diff(s[kk]).diff(s[ll]))
                        eCB += ';\n'
        eCB += '    for (ll=0; ll<'+str(np)+'; ll++) {\n'
        eCB += '        for (i=0; i<'+str(ns)+'; i++) {\n'
        eCB += '            for (j=0; j<'+str(ns)+'; j++) {\n'
        eCB += '                temp[j] = d4gdsdtdpdw[j][ll];\n'
        eCB += '                for (k=0; k<'+str(ns)+'; k++) {\n'
        eCB += '                    temp[j] += d4gds2dtdp[j][k]*dsdw[k][ll];\n'
        eCB += '                    temp[j] += d4gds2dtdw[j][k][ll]*dsdp[k];\n'
        eCB += '                    temp[j] += d4gds2dpdw[j][k][ll]*dsdt[k];\n'
        eCB += '                    temp[j] += d3gds2dt[j][k]*d2sdpdw[k][ll];\n'
        eCB += '                    temp[j] += d3gds2dp[j][k]*d2sdtdw[k][ll];\n'
        eCB += '                    temp[j] += d3gds2dw[j][k][ll]*d2sdtdp[k];\n'
        eCB += '                    for (l=0; l<'+str(ns)+'; l++) {\n' 
        eCB += '                        temp[j] += d4gds3dt[j][k][l]*dsdp[k]*dsdw[l][ll];\n'
        eCB += '                        temp[j] += d4gds3dp[j][k][l]*dsdt[k]*dsdw[l][ll];\n'
        eCB += '                        temp[j] += d4gds3dw[j][k][l][ll]*dsdt[k]*dsdp[l];\n'
        eCB += '                        temp[j] += d3gds3[j][k][l]*d2sdpdw[k][ll]'
        eCB +=                          '*dsdt[l];\n'
        eCB += '                        temp[j] += d3gds3[j][k][l]*d2sdtdw[k][ll]'
        eCB +=                          '*dsdp[l];\n'
        eCB += '                        temp[j] += d3gds3[j][k][l]*d2sdtdp[k]'
        eCB +=                          '*dsdw[l][ll];\n'
        eCB += '                        for (m=0; m<'+str(ns)+'; m++) {\n'
        eCB += '                            temp[j] += d4gds4[j][k][l][m]*'
        eCB +=                              'dsdt[k]*dsdp[l]*dsdw[m][ll];\n'
        eCB += '                        }\n'
        eCB += '                    }\n'
        eCB += '                }\n'
        eCB += '            }\n'
        eCB += '            d3sdtdpdw[i][ll] = 0.0;\n'
        eCB += '            for (j=0; j<'+str(ns)+'; j++) '
        eCB +=                'd3sdtdpdw[i][ll] += - invd2gds2[i][j]*temp[j];\n'
        eCB += '        }\n'
        eCB += '    }\n'
        eCB += '}\n'

        eCB += 'static void order_d3sdt2dp(double T, double P,' \
             + ' double n['+str(c)+'],' \
             + (' double b['+str(c)+'],' if SpMdClass else '') \
             + ' double s['+str(ns)+'],' \
             + ' double invd2gds2['+str(ns)+']['+str(ns)+'],' \
             + ' double d3sdt2dp['+str(ns)+']) {\n'
        eCB += '    double dsdt['+str(ns)+'];\n'
        eCB += '    double dsdp['+str(ns)+'];\n'
        eCB += '    double d2sdt2['+str(ns)+'];\n'
        eCB += '    double d2sdtdp['+str(ns)+'];\n'
        eCB += '    order_dsdt(T, P, n, s, invd2gds2, dsdt);\n'
        eCB += '    order_dsdp(T, P, n, s, invd2gds2, dsdp);\n'
        eCB += '    order_d2sdt2(T, P, n, s, invd2gds2, d2sdt2);\n'
        eCB += '    order_d2sdtdp(T, P, n, s, invd2gds2, d2sdtdp);\n'
        eCB += '    int i,j,k,l,m;\n'
        eCB += '    double d4gdsdt2dp['+str(ns)+'];\n'
        eCB += '    double d4gds2dtdp['+str(ns)+']['+str(ns)+'];\n'
        eCB += '    double d4gds2dt2['+str(ns)+']['+str(ns)+'];\n'
        eCB += '    double d4gds3dt['+str(ns)+']['+str(ns)+']['+str(ns)+'];\n'
        eCB += '    double d4gds3dp['+str(ns)+']['+str(ns)+']['+str(ns)+'];\n'
        eCB += '    double d4gds4['+str(ns)+']['+str(ns)+']['+str(ns)+']['+str(ns)+'];\n'
        eCB += '    double d3gds2dt['+str(ns)+']['+str(ns)+'];\n'
        eCB += '    double d3gds2dp['+str(ns)+']['+str(ns)+'];\n'
        eCB += '    double d3gds3['+str(ns)+']['+str(ns)+']['+str(ns)+'];\n'
        eCB += '    double temp['+str(ns)+'];\n'
        eCB += moles_assign_text + '\n'
        if SpMdClass:
            for i in range(0,c):
                eCB += '    double b' + str(i+1) + ' = b[' + str(i) + '];\n'
        eCB += order_assign_text + '\n'
        for k in range(0,ns):
            f = sFunc.function[k]
            eCB += '    d4gdsdt2dp['+str(k)+'] = '
            eCB += '0' if dmin else self.printer.doprint(
                f.diff(P).diff(T,2))
            eCB += ';\n'
            for l in range(0,ns):
                eCB += '    d3gds2dt['+str(k)+']['+str(l)+'] = '
                eCB += '0' if dmin else self.printer.doprint(
                    f.diff(s[l]).diff(T))
                eCB += ';\n'
                eCB += '    d4gds2dt2['+str(k)+']['+str(l)+'] = '
                eCB += '0' if dmin else self.printer.doprint(
                    f.diff(s[l]).diff(T,2))
                eCB += ';\n'
                eCB += '    d3gds2dp['+str(k)+']['+str(l)+'] = '
                eCB += '0' if dmin else self.printer.doprint(
                    f.diff(s[l]).diff(P))
                eCB += ';\n'
                eCB += '    d4gds2dtdp['+str(k)+']['+str(l)+'] = '
                eCB += '0' if dmin else self.printer.doprint(
                    f.diff(s[l]).diff(P).diff(T))
                eCB += ';\n'
                for kk in range(0,ns):
                    eCB += '    d3gds3['+str(k)+']['+str(l)+']['+str(kk)+'] = '
                    eCB += '0' if dmin else self.printer.doprint(
                        f.diff(s[l]).diff(s[kk]))
                    eCB += ';\n'
                    eCB += '    d4gds3dt['+str(k)+']['+str(l)+']['+str(kk)+'] = '
                    eCB += '0' if dmin else self.printer.doprint(
                            f.diff(s[l]).diff(s[kk]).diff(T))
                    eCB += ';\n'
                    eCB += '    d4gds3dp['+str(k)+']['+str(l)+']['+str(kk)+'] = '
                    eCB += '0' if dmin else self.printer.doprint(
                            f.diff(s[l]).diff(s[kk]).diff(P))
                    eCB += ';\n'
                    for ll in range(0,ns):
                        eCB += '    d4gds4['+str(k)+']['+str(l)+'][' \
                             +str(kk)+']['+str(ll)+'] = ' 
                        eCB += '0' if dmin else self.printer.doprint(
                            f.diff(s[l]).diff(s[kk]).diff(s[ll]))
                        eCB += ';\n'
        eCB += '    for (i=0; i<'+str(ns)+'; i++) {\n'
        eCB += '        for (j=0; j<'+str(ns)+'; j++) {\n'
        eCB += '            temp[j] = d4gdsdt2dp[j];\n'
        eCB += '            for (k=0; k<'+str(ns)+'; k++) {\n'
        eCB += '                temp[j] += d4gds2dt2[j][k]*dsdp[k];\n'
        eCB += '                temp[j] += 2.0*d4gds2dtdp[j][k]*dsdt[k];\n'
        eCB += '                temp[j] += 2.0*d3gds2dt[j][k]*d2sdtdp[k];\n'
        eCB += '                temp[j] += d3gds2dp[j][k]*d2sdt2[k];\n'
        eCB += '                for (l=0; l<'+str(ns)+'; l++) {\n' 
        eCB += '                    temp[j] += 2.0*d4gds3dt[j][k][l]*dsdt[k]*dsdp[l];\n'
        eCB += '                    temp[j] += d4gds3dp[j][k][l]*dsdt[k]*dsdt[l];\n'
        eCB += '                    temp[j] += 2.0*d3gds3[j][k][l]*d2sdtdp[k]*dsdt[l];\n'
        eCB += '                    temp[j] += d3gds3[j][k][l]*d2sdt2[k]*dsdp[l];\n'
        eCB += '                    for (m=0; m<'+str(ns)+'; m++) {\n'
        eCB += '                        temp[j] += d4gds4[j][k][l][m]*'
        eCB +=                          'dsdt[k]*dsdt[l]*dsdp[m];\n'
        eCB += '                    }\n'
        eCB += '                }\n'
        eCB += '            }\n'
        eCB += '        }\n'
        eCB += '        d3sdt2dp[i] = 0.0;\n'
        eCB += '        for (j=0; j<'+str(ns)+'; j++) '
        eCB +=            'd3sdt2dp[i] += - invd2gds2[i][j]*temp[j];\n'
        eCB += '    }\n'
        eCB += '}\n'

        ################################
        # Block of code for d4gdt2dpdw #
        ################################
        if self.verbose:
            print ("... (14/18) computing and writing code for d4gdt2dpdw")
        extraCb  = '    double d4gdsdt2dp['+str(ns)+'];\n'
        extraCb += '    double d4gds2dt2['+str(ns)+']['+str(ns)+'];\n'
        extraCb += '    double d4gds2dtdp['+str(ns)+']['+str(ns)+'];\n'
        extraCb += '    double d4gdsdt2dw['+str(ns)+']['+str(np)+'];\n'
        extraCb += '    double d4gdsdtdpdw['+str(ns)+']['+str(np)+'];\n'
        extraCb += '    double d4gds3dt['+str(ns)+']['+str(ns)+']['+str(ns)+'];\n'
        extraCb += '    double d4gds3dp['+str(ns)+']['+str(ns)+']['+str(ns)+'];\n'
        extraCb += '    double d4gds2dtdw['+str(ns)+']['+str(ns)+']['+str(np)+'];\n'
        extraCb += '    double d4gds2dpdw['+str(ns)+']['+str(ns)+']['+str(np)+'];\n'
        extraCb += '    double d4gds3dw['+str(ns)+']['+str(ns)+']['+str(ns)+']['+str(np)+'];\n'
        extraCb += '    double d4gds4['+str(ns)+']['+str(ns)+']['+str(ns)+']['+str(ns)+'];\n'
        extraCb += '    double d3gdsdt2['+str(ns)+'];\n'
        extraCb += '    double d3gdsdtdp['+str(ns)+'];\n'
        extraCb += '    double d3gds2dt['+str(ns)+']['+str(ns)+'];\n'
        extraCb += '    double d3gds2dp['+str(ns)+']['+str(ns)+'];\n'
        extraCb += '    double d3gdsdtdw['+str(ns)+']['+str(np)+'];\n'
        extraCb += '    double d3gdsdpdw['+str(ns)+']['+str(np)+'];\n'
        extraCb += '    double d3gds2dw['+str(ns)+']['+str(ns)+']['+str(np)+'];\n'
        extraCb += '    double d3gds3['+str(ns)+']['+str(ns)+']['+str(ns)+'];\n'
        extraCb += '    double d2gdsdt['+str(ns)+'];\n'
        extraCb += '    double d2gdsdp['+str(ns)+'];\n'
        extraCb += '    double d2gdsdw['+str(ns)+']['+str(np)+'];\n'
        extraCb += '    double d2gds2['+str(ns)+']['+str(ns)+'];\n'
        extraCb += '    double dsdt['+str(ns)+'];\n'
        extraCb += '    double dsdp['+str(ns)+'];\n'
        extraCb += '    double dsdw['+str(ns)+']['+str(np)+'];\n'
        extraCb += '    double d2sdt2['+str(ns)+'];\n'
        extraCb += '    double d2sdtdp['+str(ns)+'];\n'
        extraCb += '    double d2sdtdw['+str(ns)+']['+str(np)+'];\n'
        extraCb += '    double d2sdpdw['+str(ns)+']['+str(np)+'];\n'
        extraCb += '    double d3sdt2dp['+str(ns)+'];\n'
        extraCb += '    double d3sdt2dw['+str(ns)+']['+str(np)+'];\n'
        extraCb += '    double d3sdtdpdw['+str(ns)+']['+str(np)+'];\n'
        extraCb += '    int i,j,k,l;\n'
        extraCb += '    order_dsdt(T, P, n, s, invd2gds2, dsdt);\n'
        extraCb += '    order_dsdp(T, P, n, s, invd2gds2, dsdp);\n'
        extraCb += '    order_dsdw(T, P, n, s, invd2gds2, dsdw);\n'
        extraCb += '    order_d2sdt2(T, P, n, s, invd2gds2, d2sdt2);\n'
        extraCb += '    order_d2sdtdp(T, P, n, s, invd2gds2, d2sdtdp);\n'
        extraCb += '    order_d2sdtdw(T, P, n, s, invd2gds2, d2sdtdw);\n'
        extraCb += '    order_d2sdpdw(T, P, n, s, invd2gds2, d2sdpdw);\n'
        extraCb += '    order_d3sdt2dp(T, P, n, s, invd2gds2, d3sdt2dp);\n'
        extraCb += '    order_d3sdt2dw(T, P, n, s, invd2gds2, d3sdt2dw);\n'
        extraCb += '    order_d3sdtdpdw(T, P, n, s, invd2gds2, d3sdtdpdw);\n'

        for k in range(0,ns):
            extraCa  = '    d4gdsdt2dp['+str(k)+'] = '
            extraCa += '0' if dmin else self.printer.doprint(
                G.diff(T,2).diff(P).diff(s[k]))
            extraCa += ';\n'
            extraCa += '    d3gdsdt2['+str(k)+'] = '
            extraCa += '0' if dmin else self.printer.doprint(
                G.diff(T,2).diff(s[k]))
            extraCa += ';\n'
            extraCa += '    d3gdsdtdp['+str(k)+'] = '
            extraCa += '0' if dmin else self.printer.doprint(
                G.diff(T).diff(P).diff(s[k]))
            extraCa += ';\n'
            extraCa += '    d2gdsdt['+str(k)+'] = '
            extraCa += '0' if dmin else self.printer.doprint(
                G.diff(T).diff(s[k]))
            extraCa += ';\n'
            extraCa += '    d2gdsdp['+str(k)+'] = '
            extraCa += '0' if dmin else self.printer.doprint(
                G.diff(P).diff(s[k]))
            extraCa += ';\n'
            for l in range(0,np):
                extraCa += '    d4gdsdt2dw['+str(k)+']['+str(l)+'] = '
                extraCa += '0' if dmin else self.printer.doprint(
                    G.diff(T,2).diff(s[k]).diff(symparam[l]))
                extraCa += ';\n'
                extraCa += '    d4gdsdtdpdw['+str(k)+']['+str(l)+'] = '
                extraCa += '0' if dmin else self.printer.doprint(
                    G.diff(T).diff(P).diff(s[k]).diff(symparam[l]))
                extraCa += ';\n'
                extraCa += '    d3gdsdtdw['+str(k)+']['+str(l)+'] = '
                extraCa += '0' if dmin else self.printer.doprint(
                    G.diff(T).diff(s[k]).diff(symparam[l]))
                extraCa += ';\n'
                extraCa += '    d3gdsdpdw['+str(k)+']['+str(l)+'] = '
                extraCa += '0' if dmin else self.printer.doprint(
                    G.diff(P).diff(s[k]).diff(symparam[l]))
                extraCa += ';\n'
                extraCa += '    d2gdsdw['+str(k)+']['+str(l)+'] = '
                extraCa += '0' if dmin else self.printer.doprint(
                    G.diff(s[k]).diff(symparam[l]))
                extraCa += ';\n'
            for l in range(0,ns):
                extraCa += '    d4gds2dt2['+str(k)+']['+str(l)+'] = '
                extraCa += '0' if dmin else self.printer.doprint(
                    G.diff(T,2).diff(s[k]).diff(s[l]))
                extraCa += ';\n'
                extraCa += '    d4gds2dtdp['+str(k)+']['+str(l)+'] = '
                extraCa += '0' if dmin else self.printer.doprint(
                    G.diff(T).diff(P).diff(s[k]).diff(s[l]))
                extraCa += ';\n'
                extraCa += '    d3gds2dt['+str(k)+']['+str(l)+'] = '
                extraCa += '0' if dmin else self.printer.doprint(
                    G.diff(T).diff(s[k]).diff(s[l]))
                extraCa += ';\n'
                extraCa += '    d3gds2dp['+str(k)+']['+str(l)+'] = '
                extraCa += '0' if dmin else self.printer.doprint(
                    G.diff(P).diff(s[k]).diff(s[l]))
                extraCa += ';\n'
                extraCa += '    d2gds2['+str(k)+']['+str(l)+'] = '
                extraCa += '0' if dmin else self.printer.doprint(
                    G.diff(s[k]).diff(s[l]))
                extraCa += ';\n'
                for kk in range(0,ns):
                    extraCa += '    d4gds3dt['+str(k)+']['+str(l)+']['+str(kk)+'] = '
                    extraCa += '0' if dmin else self.printer.doprint(
                        G.diff(T).diff(s[k]).diff(s[l]).diff(s[kk]))
                    extraCa += ';\n'
                    extraCa += '    d4gds3dp['+str(k)+']['+str(l)+']['+str(kk)+'] = '
                    extraCa += '0' if dmin else self.printer.doprint(
                        G.diff(P).diff(s[k]).diff(s[l]).diff(s[kk]))
                    extraCa += ';\n'
                    extraCa += '    d3gds3['+str(k)+']['+str(l)+']['+str(kk)+'] = '
                    extraCa += '0' if dmin else self.printer.doprint(
                        G.diff(s[k]).diff(s[l]).diff(s[kk]))
                    extraCa += ';\n'
                    for ll in range(0,ns):
                        extraCa += '    d4gds4['+str(k)+']['+str(l)+'][' \
                                 +str(kk)+']['+str(ll)+'] = '
                        extraCa += '0' if dmin else self.printer.doprint(
                            G.diff(s[k]).diff(s[l]).diff(s[kk]).diff(s[ll]))
                        extraCa += ';\n'
                    for ll in range(0,np):
                        extraCa += '    d4gds3dw['+str(k)+']['+str(l)+'][' \
                                 +str(kk)+']['+str(ll)+'] = '
                        extraCa += '0' if dmin else self.printer.doprint(
                            G.diff(s[k]).diff(s[l]).diff(s[kk]).diff(symparam[ll]))
                        extraCa += ';\n'
                for kk in range(0,np):
                    extraCa += '    d4gds2dtdw['+str(k)+']['+str(l)+']['+str(kk)+'] = '
                    extraCa += '0' if dmin else self.printer.doprint(
                        G.diff(T).diff(s[k]).diff(s[l]).diff(symparam[kk]))
                    extraCa += ';\n'
                    extraCa += '    d4gds2dpdw['+str(k)+']['+str(l)+']['+str(kk)+'] = '
                    extraCa += '0' if dmin else self.printer.doprint(
                        G.diff(P).diff(s[k]).diff(s[l]).diff(symparam[kk]))
                    extraCa += ';\n'
                    extraCa += '    d3gds2dw['+str(k)+']['+str(l)+']['+str(kk)+'] = '
                    extraCa += '0' if dmin else self.printer.doprint(
                        G.diff(s[k]).diff(s[l]).diff(symparam[kk]))
                    extraCa += ';\n'

        j += 1
        if dmin:
            G_jac_list = ['0' for i in range(0, np)]
        else:
            G_jac_list = [self.printer.doprint(
                G_param_jac[j,i]) for i in range(0, np)]

        switch_code_text = ''
        for i in range(0, len(self.params)):
            switch_code_text += '    case ' + str(i) + ':\n'
            switch_code_text += '        result += ' + G_jac_list[i] + ';\n'

            switch_code_text += '        for (i=0; i<'+str(ns)+'; i++) {\n'
            switch_code_text += '            result += ' \
            + '2.0*d4gdsdtdpdw[i]['+str(i)+']*dsdt[i] + ' \
            + '    d4gdsdt2dw[i]['+str(i)+']*dsdp[i] + ' \
            + '    d4gdsdt2dp[i]*dsdw[i]['+str(i)+'] + ' \
            + '    d3gdsdt2[i]*d2sdpdw[i]['+str(i)+'] + ' \
            + '2.0*d3gdsdtdp[i]*d2sdtdw[i]['+str(i)+'] + ' \
            + '2.0*d3gdsdtdw[i]['+str(i)+']*d2sdtdp[i] + ' \
            + '    d3gdsdpdw[i]['+str(i)+']*d2sdt2[i] + ' \
            + '2.0*d2gdsdt[i]*d3sdtdpdw[i]['+str(i)+'] + ' \
            + '    d2gdsdp[i]*d3sdt2dw[i]['+str(i)+'] + ' \
            + '    d2gdsdw[i]['+str(i)+']*d3sdt2dp[i];\n'
            switch_code_text += '            for (j=0; j<'+str(ns)+'; j++) {\n'
            switch_code_text += '                result += ' \
            + '    d4gds2dt2[i][j]*dsdp[i]*dsdw[j]['+str(i)+'] + ' \
            + '2.0*d4gds2dtdp[i][j]*dsdt[i]*dsdw[j]['+str(i)+'] + ' \
            + '2.0*d4gds2dtdw[i][j]['+str(i)+']*dsdt[i]*dsdp[j] + ' \
            + '    d4gds2dpdw[i][j]['+str(i)+']*dsdt[i]*dsdt[j] + ' \
            + '2.0*d3gds2dt[i][j]*d2sdtdw[i]['+str(i)+']*dsdp[j] + ' \
            + '2.0*d3gds2dt[i][j]*d2sdpdw[i]['+str(i)+']*dsdt[j] + ' \
            + '2.0*d3gds2dt[i][j]*d2sdtdp[i]*dsdw[j]['+str(i)+'] + ' \
            + '2.0*d3gds2dp[i][j]*d2sdtdw[i]['+str(i)+']*dsdt[j] + ' \
            + '    d3gds2dp[i][j]*d2sdt2[i]*dsdw[j]['+str(i)+'] + ' \
            + '    d3gds2dw[i][j]['+str(i)+']*d2sdt2[i]*dsdp[j] + ' \
            + '2.0*d3gds2dw[i][j]['+str(i)+']*d2sdtdp[i]*dsdt[j] + ' \
            + '2.0*d2gds2[i][j]*d2sdtdp[i]*d2sdtdw[j]['+str(i)+'] + ' \
            + '    d2gds2[i][j]*d2sdt2[i]*d2sdpdw[j]['+str(i)+'];\n'
            switch_code_text += '                for (k=0; k<'+str(ns)+'; k++) {\n'
            switch_code_text += '                    result += ' \
            + '2.0*d4gds3dt[i][j][k]*dsdp[i]*dsdt[j]*dsdw[k]['+str(i)+'] + ' \
            + '    d4gds3dp[i][j][k]*dsdt[i]*dsdt[j]*dsdw[k]['+str(i)+'] + ' \
            + '    d4gds3dw[i][j][k]['+str(i)+']*dsdt[i]*dsdt[j]*dsdp[k] + ' \
            + '    d3gds3[i][j][k]*d2sdt2[i]*dsdp[j]*dsdw[k]['+str(i)+'] + ' \
            + '2.0*d3gds3[i][j][k]*d2sdtdp[i]*dsdt[j]*dsdw[k]['+str(i)+'] + ' \
            + '2.0*d3gds3[i][j][k]*d2sdtdw[i]['+str(i)+']*dsdt[j]*dsdp[k] + ' \
            + '    d3gds3[i][j][k]*d2sdpdw[i]['+str(i)+']*dsdt[j]*dsdt[k];\n'
            switch_code_text += '                    for (l=0; l<'+str(ns)+'; l++) ' \
            + 'result += d4gds4[i][j][k][l]*dsdt[i]*dsdt[j]*dsdp[k]*dsdw[l][' \
            + str(i)+'];\n'
            switch_code_text += '                }\n'
            switch_code_text += '            }\n'
            switch_code_text += '        }\n'

            switch_code_text += '        break;\n'

        w += solution_calib_template.format(module=self.module, 
            number_components=c, func=self.g_list[j][0], 
            switch_code=switch_code_text, moles_assign=moles_assign_text,
            number_ordering=ns, 
            order_initial_guess=order_initial_guess_text,
            order_assign=extraCb+order_assign_text+extraCa,
            extra_ordering_code=eCB)

        #######################
        # Extra ordering code #
        # d3sdp2dw, d3sdtdp2  #
        #######################
        if self.verbose:
            print ("... (15/18) computing and writing code for d3sdp2dw, d3sdtdp2")
        eCB = 'static void order_d3sdp2dw(double T, double P,' \
            + ' double n['+str(c)+'],' \
            + (' double b['+str(c)+'],' if SpMdClass else '') \
            + ' double s['+str(ns)+'],' \
            + ' double invd2gds2['+str(ns)+']['+str(ns)+'],' \
            + ' double d3sdp2dw['+str(ns)+']['+str(np)+']) {\n'
        eCB += '    double dsdp['+str(ns)+'];\n'
        eCB += '    double dsdw['+str(ns)+']['+str(np)+'];\n'
        eCB += '    double d2sdp2['+str(ns)+'];\n'
        eCB += '    double d2sdpdw['+str(ns)+']['+str(np)+'];\n'
        eCB += '    order_dsdp(T, P, n, s, invd2gds2, dsdp);\n'
        eCB += '    order_dsdw(T, P, n, s, invd2gds2, dsdw);\n'
        eCB += '    order_d2sdp2(T, P, n, s, invd2gds2, d2sdp2);\n'
        eCB += '    order_d2sdpdw(T, P, n, s, invd2gds2, d2sdpdw);\n'
        eCB += '    int i,j,k,l,m,ll;\n'
        eCB += '    double d4gdsdp2dw['+str(ns)+']['+str(np)+'];\n'
        eCB += '    double d4gds2dpdw['+str(ns)+']['+str(ns)+']['+str(np)+'];\n'
        eCB += '    double d4gds2dp2['+str(ns)+']['+str(ns)+'];\n'
        eCB += '    double d4gds3dp['+str(ns)+']['+str(ns)+']['+str(ns)+'];\n'
        eCB += '    double d4gds3dw['+str(ns)+']['+str(ns)+']['+str(ns)+']['+str(np)+'];\n'
        eCB += '    double d4gds4['+str(ns)+']['+str(ns)+']['+str(ns)+']['+str(ns)+'];\n'
        eCB += '    double d3gds2dp['+str(ns)+']['+str(ns)+'];\n'
        eCB += '    double d3gds2dw['+str(ns)+']['+str(ns)+']['+str(np)+'];\n'
        eCB += '    double d3gds3['+str(ns)+']['+str(ns)+']['+str(ns)+'];\n'
        eCB += '    double temp['+str(ns)+'];\n'
        eCB += moles_assign_text + '\n'
        if SpMdClass:
            for i in range(0,c):
                eCB += '    double b' + str(i+1) + ' = b[' + str(i) + '];\n'
        eCB += order_assign_text + '\n'
        for k in range(0,ns):
            f = sFunc.function[k]
            for l in range(0,np):
                eCB += '    d4gdsdp2dw['+str(k)+']['+str(l)+'] = '
                eCB += '0' if dmin else self.printer.doprint(
                        f.diff(symparam[l]).diff(P,2))
                eCB += ';\n'
            for l in range(0,ns):
                eCB += '    d3gds2dp['+str(k)+']['+str(l)+'] = '
                eCB += '0' if dmin else self.printer.doprint(
                    f.diff(s[l]).diff(P))
                eCB += ';\n'
                eCB += '    d4gds2dp2['+str(k)+']['+str(l)+'] = '
                eCB += '0' if dmin else self.printer.doprint(
                    f.diff(s[l]).diff(P,2))
                eCB += ';\n'
                for kk in range(0,np):
                    eCB += '    d3gds2dw['+str(k)+']['+str(l)+']['+str(kk)+'] = '
                    eCB += '0' if dmin else self.printer.doprint(
                            f.diff(s[l]).diff(symparam[kk]))
                    eCB += ';\n'
                    eCB += '    d4gds2dpdw['+str(k)+']['+str(l)+']['+str(kk)+'] = '
                    eCB += '0' if dmin else self.printer.doprint(
                            f.diff(s[l]).diff(symparam[kk]).diff(P))
                    eCB += ';\n'
                for kk in range(0,ns):
                    eCB += '    d3gds3['+str(k)+']['+str(l)+']['+str(kk)+'] = '
                    eCB += '0' if dmin else self.printer.doprint(
                        f.diff(s[l]).diff(s[kk]))
                    eCB += ';\n'
                    eCB += '    d4gds3dp['+str(k)+']['+str(l)+']['+str(kk)+'] = '
                    eCB += '0' if dmin else self.printer.doprint(
                            f.diff(s[l]).diff(s[kk]).diff(P))
                    eCB += ';\n'
                    for ll in range(0,np):
                        eCB += '    d4gds3dw['+str(k)+']['+str(l)+'][' \
                             +str(kk)+']['+str(ll)+'] = '
                        eCB += '0' if dmin else self.printer.doprint(
                            f.diff(s[l]).diff(s[kk]).diff(symparam[ll]))
                        eCB += ';\n'
                    for ll in range(0,ns):
                        eCB += '    d4gds4['+str(k)+']['+str(l)+'][' \
                             +str(kk)+']['+str(ll)+'] = '
                        eCB += '0' if dmin else self.printer.doprint(
                            f.diff(s[l]).diff(s[kk]).diff(s[ll]))
                        eCB += ';\n'
        eCB += '    for (ll=0; ll<'+str(np)+'; ll++) {\n'
        eCB += '        for (i=0; i<'+str(ns)+'; i++) {\n'
        eCB += '            for (j=0; j<'+str(ns)+'; j++) {\n'
        eCB += '                temp[j] = d4gdsdp2dw[j][ll];\n'
        eCB += '                for (k=0; k<'+str(ns)+'; k++) {\n'
        eCB += '                    temp[j] += d4gds2dp2[j][k]*dsdw[k][ll];\n'
        eCB += '                    temp[j] += 2.0*d4gds2dpdw[j][k][ll]*dsdp[k];\n'
        eCB += '                    temp[j] += 2.0*d3gds2dp[j][k]*d2sdpdw[k][ll];\n'
        eCB += '                    temp[j] += d3gds2dw[j][k][ll]*d2sdp2[k];\n'
        eCB += '                    for (l=0; l<'+str(ns)+'; l++) {\n' 
        eCB += '                        temp[j] += 2.0*d4gds3dp[j][k][l]*dsdp[k]*dsdw[l][ll];\n'
        eCB += '                        temp[j] += d4gds3dw[j][k][l][ll]*dsdp[k]*dsdp[l];\n'
        eCB += '                        temp[j] += 2.0*d3gds3[j][k][l]*d2sdpdw[k][ll]'
        eCB +=                          '*dsdp[l];\n'
        eCB += '                        temp[j] += d3gds3[j][k][l]*d2sdp2[k]'
        eCB +=                          '*dsdw[l][ll];\n'
        eCB += '                        for (m=0; m<'+str(ns)+'; m++) {\n'
        eCB += '                            temp[j] += d4gds4[j][k][l][m]*'
        eCB +=                              'dsdp[k]*dsdp[l]*dsdw[m][ll];\n'
        eCB += '                        }\n'
        eCB += '                    }\n'
        eCB += '                }\n'
        eCB += '            }\n'
        eCB += '            d3sdp2dw[i][ll] = 0.0;\n'
        eCB += '            for (j=0; j<'+str(ns)+'; j++) '
        eCB +=                'd3sdp2dw[i][ll] += - invd2gds2[i][j]*temp[j];\n'
        eCB += '        }\n'
        eCB += '    }\n'
        eCB += '}\n'

        eCB += 'static void order_d3sdtdp2(double T, double P,' \
             + ' double n['+str(c)+'],' \
             + (' double b['+str(c)+'],' if SpMdClass else '') \
             + ' double s['+str(ns)+'],' \
             + ' double invd2gds2['+str(ns)+']['+str(ns)+'],' \
             + ' double d3sdtdp2['+str(ns)+']) {\n'
        eCB += '    double dsdt['+str(ns)+'];\n'
        eCB += '    double dsdp['+str(ns)+'];\n'
        eCB += '    double d2sdp2['+str(ns)+'];\n'
        eCB += '    double d2sdtdp['+str(ns)+'];\n'
        eCB += '    order_dsdt(T, P, n, s, invd2gds2, dsdt);\n'
        eCB += '    order_dsdp(T, P, n, s, invd2gds2, dsdp);\n'
        eCB += '    order_d2sdp2(T, P, n, s, invd2gds2, d2sdp2);\n'
        eCB += '    order_d2sdtdp(T, P, n, s, invd2gds2, d2sdtdp);\n'
        eCB += '    int i,j,k,l,m;\n'
        eCB += '    double d4gdsdtdp2['+str(ns)+'];\n'
        eCB += '    double d4gds2dtdp['+str(ns)+']['+str(ns)+'];\n'
        eCB += '    double d4gds2dp2['+str(ns)+']['+str(ns)+'];\n'
        eCB += '    double d4gds3dt['+str(ns)+']['+str(ns)+']['+str(ns)+'];\n'
        eCB += '    double d4gds3dp['+str(ns)+']['+str(ns)+']['+str(ns)+'];\n'
        eCB += '    double d4gds4['+str(ns)+']['+str(ns)+']['+str(ns)+']['+str(ns)+'];\n'
        eCB += '    double d3gds2dt['+str(ns)+']['+str(ns)+'];\n'
        eCB += '    double d3gds2dp['+str(ns)+']['+str(ns)+'];\n'
        eCB += '    double d3gds3['+str(ns)+']['+str(ns)+']['+str(ns)+'];\n'
        eCB += '    double temp['+str(ns)+'];\n'
        eCB += moles_assign_text + '\n'
        if SpMdClass:
            for i in range(0,c):
                eCB += '    double b' + str(i+1) + ' = b[' + str(i) + '];\n'
        eCB += order_assign_text + '\n'
        for k in range(0,ns):
            f = sFunc.function[k]
            eCB += '    d4gdsdtdp2['+str(k)+'] = '
            eCB += '0' if dmin else self.printer.doprint(
                f.diff(P,2).diff(T))
            eCB += ';\n'
            for l in range(0,ns):
                eCB += '    d3gds2dt['+str(k)+']['+str(l)+'] = '
                eCB += '0' if dmin else self.printer.doprint(
                    f.diff(s[l]).diff(T))
                eCB += ';\n'
                eCB += '    d4gds2dp2['+str(k)+']['+str(l)+'] = '
                eCB += '0' if dmin else self.printer.doprint(
                    f.diff(s[l]).diff(P,2))
                eCB += ';\n'
                eCB += '    d3gds2dp['+str(k)+']['+str(l)+'] = '
                eCB += '0' if dmin else self.printer.doprint(
                    f.diff(s[l]).diff(P))
                eCB += ';\n'
                eCB += '    d4gds2dtdp['+str(k)+']['+str(l)+'] = '
                eCB += '0' if dmin else self.printer.doprint(
                    f.diff(s[l]).diff(P).diff(T))
                eCB += ';\n'
                for kk in range(0,ns):
                    eCB += '    d3gds3['+str(k)+']['+str(l)+']['+str(kk)+'] = '
                    eCB += '0' if dmin else self.printer.doprint(
                        f.diff(s[l]).diff(s[kk]))
                    eCB += ';\n'
                    eCB += '    d4gds3dt['+str(k)+']['+str(l)+']['+str(kk)+'] = '
                    eCB += '0' if dmin else self.printer.doprint(
                            f.diff(s[l]).diff(s[kk]).diff(T))
                    eCB += ';\n'
                    eCB += '    d4gds3dp['+str(k)+']['+str(l)+']['+str(kk)+'] = '
                    eCB += '0' if dmin else self.printer.doprint(
                            f.diff(s[l]).diff(s[kk]).diff(P))
                    eCB += ';\n'
                    for ll in range(0,ns):
                        eCB += '    d4gds4['+str(k)+']['+str(l)+'][' \
                             +str(kk)+']['+str(ll)+'] = '
                        eCB += '0' if dmin else self.printer.doprint(
                            f.diff(s[l]).diff(s[kk]).diff(s[ll]))
                        eCB += ';\n'
        eCB += '    for (i=0; i<'+str(ns)+'; i++) {\n'
        eCB += '        for (j=0; j<'+str(ns)+'; j++) {\n'
        eCB += '            temp[j] = d4gdsdtdp2[j];\n'
        eCB += '            for (k=0; k<'+str(ns)+'; k++) {\n'
        eCB += '                temp[j] += d4gds2dp2[j][k]*dsdt[k];\n'
        eCB += '                temp[j] += 2.0*d4gds2dtdp[j][k]*dsdp[k];\n'
        eCB += '                temp[j] += 2.0*d3gds2dp[j][k]*d2sdtdp[k];\n'
        eCB += '                temp[j] += d3gds2dt[j][k]*d2sdp2[k];\n'
        eCB += '                for (l=0; l<'+str(ns)+'; l++) {\n' 
        eCB += '                    temp[j] += 2.0*d4gds3dp[j][k][l]*dsdt[k]*dsdp[l];\n'
        eCB += '                    temp[j] += d4gds3dt[j][k][l]*dsdp[k]*dsdp[l];\n'
        eCB += '                    temp[j] += 2.0*d3gds3[j][k][l]*d2sdtdp[k]*dsdp[l];\n'
        eCB += '                    temp[j] += d3gds3[j][k][l]*d2sdp2[k]*dsdt[l];\n'
        eCB += '                    for (m=0; m<'+str(ns)+'; m++) {\n'
        eCB += '                        temp[j] += d4gds4[j][k][l][m]*'
        eCB +=                          'dsdp[k]*dsdp[l]*dsdt[m];\n'
        eCB += '                    }\n'
        eCB += '                }\n'
        eCB += '            }\n'
        eCB += '        }\n'
        eCB += '        d3sdtdp2[i] = 0.0;\n'
        eCB += '        for (j=0; j<'+str(ns)+'; j++) '
        eCB +=            'd3sdtdp2[i] += - invd2gds2[i][j]*temp[j];\n'
        eCB += '    }\n'
        eCB += '}\n'

        ################################
        # Block of code for d4gdtdp2dw #
        ################################
        if self.verbose:
            print ("... (16/18) computing and writing code for d4gdtdp2dw")
        extraCb  = '    double d4gdsdtdp2['+str(ns)+'];\n'
        extraCb += '    double d4gds2dp2['+str(ns)+']['+str(ns)+'];\n'
        extraCb += '    double d4gds2dtdp['+str(ns)+']['+str(ns)+'];\n'
        extraCb += '    double d4gdsdp2dw['+str(ns)+']['+str(np)+'];\n'
        extraCb += '    double d4gdsdtdpdw['+str(ns)+']['+str(np)+'];\n'
        extraCb += '    double d4gds3dp['+str(ns)+']['+str(ns)+']['+str(ns)+'];\n'
        extraCb += '    double d4gds3dt['+str(ns)+']['+str(ns)+']['+str(ns)+'];\n'
        extraCb += '    double d4gds2dpdw['+str(ns)+']['+str(ns)+']['+str(np)+'];\n'
        extraCb += '    double d4gds2dtdw['+str(ns)+']['+str(ns)+']['+str(np)+'];\n'
        extraCb += '    double d4gds3dw['+str(ns)+']['+str(ns)+']['+str(ns)+']['+str(np)+'];\n'
        extraCb += '    double d4gds4['+str(ns)+']['+str(ns)+']['+str(ns)+']['+str(ns)+'];\n'
        extraCb += '    double d3gdsdp2['+str(ns)+'];\n'
        extraCb += '    double d3gdsdtdp['+str(ns)+'];\n'
        extraCb += '    double d3gds2dp['+str(ns)+']['+str(ns)+'];\n'
        extraCb += '    double d3gds2dt['+str(ns)+']['+str(ns)+'];\n'
        extraCb += '    double d3gdsdpdw['+str(ns)+']['+str(np)+'];\n'
        extraCb += '    double d3gdsdtdw['+str(ns)+']['+str(np)+'];\n'
        extraCb += '    double d3gds2dw['+str(ns)+']['+str(ns)+']['+str(np)+'];\n'
        extraCb += '    double d3gds3['+str(ns)+']['+str(ns)+']['+str(ns)+'];\n'
        extraCb += '    double d2gdsdp['+str(ns)+'];\n'
        extraCb += '    double d2gdsdt['+str(ns)+'];\n'
        extraCb += '    double d2gdsdw['+str(ns)+']['+str(np)+'];\n'
        extraCb += '    double d2gds2['+str(ns)+']['+str(ns)+'];\n'
        extraCb += '    double dsdp['+str(ns)+'];\n'
        extraCb += '    double dsdt['+str(ns)+'];\n'
        extraCb += '    double dsdw['+str(ns)+']['+str(np)+'];\n'
        extraCb += '    double d2sdp2['+str(ns)+'];\n'
        extraCb += '    double d2sdtdp['+str(ns)+'];\n'
        extraCb += '    double d2sdpdw['+str(ns)+']['+str(np)+'];\n'
        extraCb += '    double d2sdtdw['+str(ns)+']['+str(np)+'];\n'
        extraCb += '    double d3sdtdp2['+str(ns)+'];\n'
        extraCb += '    double d3sdp2dw['+str(ns)+']['+str(np)+'];\n'
        extraCb += '    double d3sdtdpdw['+str(ns)+']['+str(np)+'];\n'
        extraCb += '    int i,j,k,l;\n'
        extraCb += '    order_dsdp(T, P, n, s, invd2gds2, dsdp);\n'
        extraCb += '    order_dsdt(T, P, n, s, invd2gds2, dsdt);\n'
        extraCb += '    order_dsdw(T, P, n, s, invd2gds2, dsdw);\n'
        extraCb += '    order_d2sdp2(T, P, n, s, invd2gds2, d2sdp2);\n'
        extraCb += '    order_d2sdtdp(T, P, n, s, invd2gds2, d2sdtdp);\n'
        extraCb += '    order_d2sdpdw(T, P, n, s, invd2gds2, d2sdpdw);\n'
        extraCb += '    order_d2sdtdw(T, P, n, s, invd2gds2, d2sdtdw);\n'
        extraCb += '    order_d3sdtdp2(T, P, n, s, invd2gds2, d3sdtdp2);\n'
        extraCb += '    order_d3sdp2dw(T, P, n, s, invd2gds2, d3sdp2dw);\n'
        extraCb += '    order_d3sdtdpdw(T, P, n, s, invd2gds2, d3sdtdpdw);\n'

        for k in range(0,ns):
            extraCa  = '    d4gdsdtdp2['+str(k)+'] = '
            extraCa += '0' if dmin else self.printer.doprint(
                G.diff(P,2).diff(T).diff(s[k]))
            extraCa += ';\n'
            extraCa += '    d3gdsdp2['+str(k)+'] = '
            extraCa += '0' if dmin else self.printer.doprint(
                G.diff(P,2).diff(s[k]))
            extraCa += ';\n'
            extraCa += '    d3gdsdtdp['+str(k)+'] = '
            extraCa += '0' if dmin else self.printer.doprint(
                G.diff(T).diff(P).diff(s[k]))
            extraCa += ';\n'
            extraCa += '    d2gdsdp['+str(k)+'] = '
            extraCa += '0' if dmin else self.printer.doprint(
                G.diff(P).diff(s[k]))
            extraCa += ';\n'
            extraCa += '    d2gdsdt['+str(k)+'] = '
            extraCa += '0' if dmin else self.printer.doprint(
                G.diff(T).diff(s[k]))
            extraCa += ';\n'
            for l in range(0,np):
                extraCa += '    d4gdsdp2dw['+str(k)+']['+str(l)+'] = '
                extraCa += '0' if dmin else self.printer.doprint(
                    G.diff(P,2).diff(s[k]).diff(symparam[l]))
                extraCa += ';\n'
                extraCa += '    d4gdsdtdpdw['+str(k)+']['+str(l)+'] = '
                extraCa += '0' if dmin else self.printer.doprint(
                    G.diff(T).diff(P).diff(s[k]).diff(symparam[l]))
                extraCa += ';\n'
                extraCa += '    d3gdsdpdw['+str(k)+']['+str(l)+'] = '
                extraCa += '0' if dmin else self.printer.doprint(
                    G.diff(P).diff(s[k]).diff(symparam[l]))
                extraCa += ';\n'
                extraCa += '    d3gdsdtdw['+str(k)+']['+str(l)+'] = '
                extraCa += '0' if dmin else self.printer.doprint(
                    G.diff(T).diff(s[k]).diff(symparam[l]))
                extraCa += ';\n'
                extraCa += '    d2gdsdw['+str(k)+']['+str(l)+'] = '
                extraCa += '0' if dmin else self.printer.doprint(
                    G.diff(s[k]).diff(symparam[l]))
                extraCa += ';\n'
            for l in range(0,ns):
                extraCa += '    d4gds2dp2['+str(k)+']['+str(l)+'] = '
                extraCa += '0' if dmin else self.printer.doprint(
                    G.diff(P,2).diff(s[k]).diff(s[l]))
                extraCa += ';\n'
                extraCa += '    d4gds2dtdp['+str(k)+']['+str(l)+'] = '
                extraCa += '0' if dmin else self.printer.doprint(
                    G.diff(T).diff(P).diff(s[k]).diff(s[l]))
                extraCa += ';\n'
                extraCa += '    d3gds2dp['+str(k)+']['+str(l)+'] = '
                extraCa += '0' if dmin else self.printer.doprint(
                    G.diff(P).diff(s[k]).diff(s[l]))
                extraCa += ';\n'
                extraCa += '    d3gds2dt['+str(k)+']['+str(l)+'] = '
                extraCa += '0' if dmin else self.printer.doprint(
                    G.diff(T).diff(s[k]).diff(s[l]))
                extraCa += ';\n'
                extraCa += '    d2gds2['+str(k)+']['+str(l)+'] = '
                extraCa += '0' if dmin else self.printer.doprint(
                    G.diff(s[k]).diff(s[l]))
                extraCa += ';\n'
                for kk in range(0,ns):
                    extraCa += '    d4gds3dp['+str(k)+']['+str(l)+']['+str(kk)+'] = '
                    extraCa += '0' if dmin else self.printer.doprint(
                        G.diff(P).diff(s[k]).diff(s[l]).diff(s[kk]))
                    extraCa += ';\n'
                    extraCa += '    d4gds3dt['+str(k)+']['+str(l)+']['+str(kk)+'] = '
                    extraCa += '0' if dmin else self.printer.doprint(
                        G.diff(T).diff(s[k]).diff(s[l]).diff(s[kk]))
                    extraCa += ';\n'
                    extraCa += '    d3gds3['+str(k)+']['+str(l)+']['+str(kk)+'] = '
                    extraCa += '0' if dmin else self.printer.doprint(
                        G.diff(s[k]).diff(s[l]).diff(s[kk]))
                    extraCa += ';\n'
                    for ll in range(0,ns):
                        extraCa += '    d4gds4['+str(k)+']['+str(l)+'][' \
                                 +str(kk)+']['+str(ll)+'] = '
                        extraCa += '0' if dmin else self.printer.doprint(
                            G.diff(s[k]).diff(s[l]).diff(s[kk]).diff(s[ll]))
                        extraCa += ';\n'
                    for ll in range(0,np):
                        extraCa += '    d4gds3dw['+str(k)+']['+str(l)+'][' \
                                 +str(kk)+']['+str(ll)+'] = '
                        extraCa += '0' if dmin else self.printer.doprint(
                            G.diff(s[k]).diff(s[l]).diff(s[kk]).diff(symparam[ll]))
                        extraCa += ';\n'
                for kk in range(0,np):
                    extraCa += '    d4gds2dpdw['+str(k)+']['+str(l)+']['+str(kk)+'] = '
                    extraCa += '0' if dmin else self.printer.doprint(
                        G.diff(P).diff(s[k]).diff(s[l]).diff(symparam[kk]))
                    extraCa += ';\n'
                    extraCa += '    d4gds2dtdw['+str(k)+']['+str(l)+']['+str(kk)+'] = '
                    extraCa += '0' if dmin else self.printer.doprint(
                        G.diff(T).diff(s[k]).diff(s[l]).diff(symparam[kk]))
                    extraCa += ';\n'
                    extraCa += '    d3gds2dw['+str(k)+']['+str(l)+']['+str(kk)+'] = '
                    extraCa += '0' if dmin else self.printer.doprint(
                        G.diff(s[k]).diff(s[l]).diff(symparam[kk]))
                    extraCa += ';\n'

        j += 1
        if dmin:
            G_jac_list = ['0' for i in range(0, np)]
        else:
            G_jac_list = [self.printer.doprint(
                G_param_jac[j,i]) for i in range(0, np)]

        switch_code_text = ''
        for i in range(0, len(self.params)):
            switch_code_text += '    case ' + str(i) + ':\n'
            switch_code_text += '        result += ' + G_jac_list[i] + ';\n'

            switch_code_text += '        for (i=0; i<'+str(ns)+'; i++) {\n'
            switch_code_text += '            result += ' \
            + '2.0*d4gdsdtdpdw[i]['+str(i)+']*dsdp[i] + ' \
            + '    d4gdsdp2dw[i]['+str(i)+']*dsdt[i] + ' \
            + '    d4gdsdtdp2[i]*dsdw[i]['+str(i)+'] + ' \
            + '    d3gdsdp2[i]*d2sdtdw[i]['+str(i)+'] + ' \
            + '2.0*d3gdsdtdp[i]*d2sdpdw[i]['+str(i)+'] + ' \
            + '2.0*d3gdsdpdw[i]['+str(i)+']*d2sdtdp[i] + ' \
            + '    d3gdsdtdw[i]['+str(i)+']*d2sdp2[i] + ' \
            + '2.0*d2gdsdp[i]*d3sdtdpdw[i]['+str(i)+'] + ' \
            + '    d2gdsdt[i]*d3sdp2dw[i]['+str(i)+'] + ' \
            + '    d2gdsdw[i]['+str(i)+']*d3sdtdp2[i];\n'
            switch_code_text += '            for (j=0; j<'+str(ns)+'; j++) {\n'
            switch_code_text += '                result += ' \
            + '    d4gds2dp2[i][j]*dsdt[i]*dsdw[j]['+str(i)+'] + ' \
            + '2.0*d4gds2dtdp[i][j]*dsdp[i]*dsdw[j]['+str(i)+'] + ' \
            + '2.0*d4gds2dpdw[i][j]['+str(i)+']*dsdp[i]*dsdt[j] + ' \
            + '    d4gds2dtdw[i][j]['+str(i)+']*dsdp[i]*dsdp[j] + ' \
            + '2.0*d3gds2dp[i][j]*d2sdpdw[i]['+str(i)+']*dsdt[j] + ' \
            + '2.0*d3gds2dp[i][j]*d2sdtdw[i]['+str(i)+']*dsdp[j] + ' \
            + '2.0*d3gds2dp[i][j]*d2sdtdp[i]*dsdw[j]['+str(i)+'] + ' \
            + '2.0*d3gds2dt[i][j]*d2sdpdw[i]['+str(i)+']*dsdp[j] + ' \
            + '    d3gds2dt[i][j]*d2sdp2[i]*dsdw[j]['+str(i)+'] + ' \
            + '    d3gds2dw[i][j]['+str(i)+']*d2sdp2[i]*dsdt[j] + ' \
            + '2.0*d3gds2dw[i][j]['+str(i)+']*d2sdtdp[i]*dsdp[j] + ' \
            + '2.0*d2gds2[i][j]*d2sdtdp[i]*d2sdpdw[j]['+str(i)+'] + ' \
            + '    d2gds2[i][j]*d2sdp2[i]*d2sdtdw[j]['+str(i)+'];\n'
            switch_code_text += '                for (k=0; k<'+str(ns)+'; k++) {\n'
            switch_code_text += '                    result += ' \
            + '2.0*d4gds3dp[i][j][k]*dsdt[i]*dsdp[j]*dsdw[k]['+str(i)+'] + ' \
            + '    d4gds3dt[i][j][k]*dsdp[i]*dsdp[j]*dsdw[k]['+str(i)+'] + ' \
            + '    d4gds3dw[i][j][k]['+str(i)+']*dsdt[i]*dsdp[j]*dsdp[k] + ' \
            + '    d3gds3[i][j][k]*d2sdp2[i]*dsdt[j]*dsdw[k]['+str(i)+'] + ' \
            + '2.0*d3gds3[i][j][k]*d2sdtdp[i]*dsdp[j]*dsdw[k]['+str(i)+'] + ' \
            + '2.0*d3gds3[i][j][k]*d2sdpdw[i]['+str(i)+']*dsdp[j]*dsdt[k] + ' \
            + '    d3gds3[i][j][k]*d2sdtdw[i]['+str(i)+']*dsdp[j]*dsdp[k];\n'
            switch_code_text += '                    for (l=0; l<'+str(ns)+'; l++) ' \
            + 'result += d4gds4[i][j][k][l]*dsdp[i]*dsdp[j]*dsdt[k]*dsdw[l][' \
            + str(i)+'];\n'
            switch_code_text += '                }\n'
            switch_code_text += '            }\n'
            switch_code_text += '        }\n'

            switch_code_text += '        break;\n'

        w += solution_calib_template.format(module=self.module, 
            number_components=c, func=self.g_list[j][0], 
            switch_code=switch_code_text, moles_assign=moles_assign_text,
            number_ordering=ns, 
            order_initial_guess=order_initial_guess_text,
            order_assign=extraCb+order_assign_text+extraCa,
            extra_ordering_code=eCB)

        #######################
        # Extra ordering code #
        # d3sdp3              #
        #######################
        if self.verbose:
            print ("... (17/18) computing and writing code for d3sdp3")
        eCB  = 'static void order_d3sdp3(double T, double P,' \
            +  ' double n['+str(c)+'],' \
            +  (' double b['+str(c)+'],' if SpMdClass else '') \
            +  ' double s['+str(ns)+'],' \
            +  ' double invd2gds2['+str(ns)+']['+str(ns)+'],' \
            +  ' double d3sdp3['+str(ns)+']) {\n'
        eCB += '    double dsdp['+str(ns)+'];\n'
        eCB += '    double d2sdp2['+str(ns)+'];\n'
        eCB += '    order_dsdp(T, P, n, s, invd2gds2, dsdp);\n'
        eCB += '    order_d2sdp2(T, P, n, s, invd2gds2, d2sdp2);\n'
        eCB += '    int i,j,k,l,m;\n'
        eCB += '    double d4gdsdp3['+str(ns)+'];\n'
        eCB += '    double d4gds2dp2['+str(ns)+']['+str(ns)+'];\n'
        eCB += '    double d4gds3dp['+str(ns)+']['+str(ns)+']['+str(ns)+'];\n'
        eCB += '    double d4gds4['+str(ns)+']['+str(ns)+']['+str(ns)+']['+str(ns)+'];\n'
        eCB += '    double d3gds2dp['+str(ns)+']['+str(ns)+'];\n'
        eCB += '    double d3gds3['+str(ns)+']['+str(ns)+']['+str(ns)+'];\n'
        eCB += '    double temp['+str(ns)+'];\n'
        eCB += moles_assign_text + '\n'
        if SpMdClass:
            for i in range(0,c):
                eCB += '    double b' + str(i+1) + ' = b[' + str(i) + '];\n'
        eCB += order_assign_text + '\n'
        for k in range(0,ns):
            f = sFunc.function[k]
            for l in range(0,ns):
                eCB += '    d3gds2dp['+str(k)+']['+str(l)+'] = '
                eCB += '0' if dmin else self.printer.doprint(
                    f.diff(s[l]).diff(P))
                eCB += ';\n'
                eCB += '    d4gds2dp2['+str(k)+']['+str(l)+'] = '
                eCB += '0' if dmin else self.printer.doprint(
                    f.diff(s[l]).diff(P,2))
                eCB += ';\n'
                for kk in range(0,ns):
                    eCB += '    d3gds3['+str(k)+']['+str(l)+']['+str(kk)+'] = '
                    eCB += '0' if dmin else self.printer.doprint(
                        f.diff(s[l]).diff(s[kk]))
                    eCB += ';\n'
                    eCB += '    d4gds3dp['+str(k)+']['+str(l)+']['+str(kk)+'] = '
                    eCB += '0' if dmin else self.printer.doprint(
                            f.diff(s[l]).diff(s[kk]).diff(P))
                    eCB += ';\n'
                    for ll in range(0,ns):
                        eCB += '    d4gds4['+str(k)+']['+str(l)+'][' \
                             +str(kk)+']['+str(ll)+'] = '
                        eCB += '0' if dmin else self.printer.doprint(
                            f.diff(s[l]).diff(s[kk]).diff(s[ll]))
                        eCB += ';\n'
        eCB += '    for (i=0; i<'+str(ns)+'; i++) {\n'
        eCB += '        for (j=0; j<'+str(ns)+'; j++) {\n'
        eCB += '            temp[j] = d4gdsdp3[j];\n'
        eCB += '            for (k=0; k<'+str(ns)+'; k++) {\n'
        eCB += '                temp[j] += 3.0*d4gds2dp2[j][k]*dsdp[k];\n'
        eCB += '                temp[j] += 3.0*d3gds2dp[j][k]*d2sdp2[k];\n'
        eCB += '                for (l=0; l<'+str(ns)+'; l++) {\n' 
        eCB += '                    temp[j] += 3.0*d4gds3dp[j][k][l]*dsdp[k]*dsdp[l];\n'
        eCB += '                    temp[j] += 3.0*d3gds3[j][k][l]*d2sdp2[k]*dsdp[l];\n'
        eCB += '                    for (m=0; m<'+str(ns)+'; m++) {\n'
        eCB += '                        temp[j] += d4gds4[j][k][l][m]*'
        eCB +=                          'dsdp[k]*dsdp[l]*dsdp[m];\n'
        eCB += '                    }\n'
        eCB += '                }\n'
        eCB += '            }\n'
        eCB += '        }\n'
        eCB += '        d3sdp3[i] = 0.0;\n'
        eCB += '        for (j=0; j<'+str(ns)+'; j++) '
        eCB +=            'd3sdp3[i] += - invd2gds2[i][j]*temp[j];\n'
        eCB += '    }\n'
        eCB += '}\n'

        ##############################
        # Block of code for d4gdp3dw #
        ##############################
        if self.verbose:
            print ("... (18/18) computing and writing code for d4gdp3dw")
        extraCb  = '    double d4gdsdp3['+str(ns)+'];\n'
        extraCb += '    double d4gds2dp2['+str(ns)+']['+str(ns)+'];\n'
        extraCb += '    double d4gdsdp2dw['+str(ns)+']['+str(np)+'];\n'
        extraCb += '    double d4gds3dp['+str(ns)+']['+str(ns)+']['+str(ns)+'];\n'
        extraCb += '    double d4gds2dpdw['+str(ns)+']['+str(ns)+']['+str(np)+'];\n'
        extraCb += '    double d4gds3dw['+str(ns)+']['+str(ns)+']['+str(ns)+']['+str(np)+'];\n'
        extraCb += '    double d4gds4['+str(ns)+']['+str(ns)+']['+str(ns)+']['+str(ns)+'];\n'
        extraCb += '    double d3gdsdp2['+str(ns)+'];\n'
        extraCb += '    double d3gds2dp['+str(ns)+']['+str(ns)+'];\n'
        extraCb += '    double d3gdsdpdw['+str(ns)+']['+str(np)+'];\n'
        extraCb += '    double d3gds2dw['+str(ns)+']['+str(ns)+']['+str(np)+'];\n'
        extraCb += '    double d3gds3['+str(ns)+']['+str(ns)+']['+str(ns)+'];\n'
        extraCb += '    double d2gdsdp['+str(ns)+'];\n'
        extraCb += '    double d2gdsdw['+str(ns)+']['+str(np)+'];\n'
        extraCb += '    double d2gds2['+str(ns)+']['+str(ns)+'];\n'
        extraCb += '    double dsdp['+str(ns)+'];\n'
        extraCb += '    double dsdw['+str(ns)+']['+str(np)+'];\n'
        extraCb += '    double d2sdp2['+str(ns)+'];\n'
        extraCb += '    double d2sdpdw['+str(ns)+']['+str(np)+'];\n'
        extraCb += '    double d3sdp3['+str(ns)+'];\n'
        extraCb += '    double d3sdp2dw['+str(ns)+']['+str(np)+'];\n'
        extraCb += '    int i,j,k,l;\n'
        extraCb += '    order_dsdp(T, P, n, s, invd2gds2, dsdp);\n'
        extraCb += '    order_dsdw(T, P, n, s, invd2gds2, dsdw);\n'
        extraCb += '    order_d2sdp2(T, P, n, s, invd2gds2, d2sdp2);\n'
        extraCb += '    order_d2sdpdw(T, P, n, s, invd2gds2, d2sdpdw);\n'
        extraCb += '    order_d3sdp3(T, P, n, s, invd2gds2, d3sdp3);\n'
        extraCb += '    order_d3sdp2dw(T, P, n, s, invd2gds2, d3sdp2dw);\n'

        for k in range(0,ns):
            extraCa  = '    d4gdsdp3['+str(k)+'] = '
            extraCa += '0' if dmin else self.printer.doprint(
                G.diff(P,3).diff(s[k]))
            extraCa += ';\n'
            extraCa += '    d3gdsdp2['+str(k)+'] = '
            extraCa += '0' if dmin else self.printer.doprint(
                G.diff(P,2).diff(s[k]))
            extraCa += ';\n'
            extraCa += '    d2gdsdp['+str(k)+'] = '
            extraCa += '0' if dmin else self.printer.doprint(
                G.diff(P).diff(s[k]))
            extraCa += ';\n'
            for l in range(0,np):
                extraCa += '    d4gdsdp2dw['+str(k)+']['+str(l)+'] = '
                extraCa += '0' if dmin else self.printer.doprint(
                    G.diff(P,2).diff(s[k]).diff(symparam[l]))
                extraCa += ';\n'
                extraCa += '    d3gdsdpdw['+str(k)+']['+str(l)+'] = '
                extraCa += '0' if dmin else self.printer.doprint(
                    G.diff(P).diff(s[k]).diff(symparam[l]))
                extraCa += ';\n'
                extraCa += '    d2gdsdw['+str(k)+']['+str(l)+'] = '
                extraCa += '0' if dmin else self.printer.doprint(
                    G.diff(s[k]).diff(symparam[l]))
                extraCa += ';\n'
            for l in range(0,ns):
                extraCa += '    d4gds2dp2['+str(k)+']['+str(l)+'] = '
                extraCa += '0' if dmin else self.printer.doprint(
                    G.diff(P,2).diff(s[k]).diff(s[l]))
                extraCa += ';\n'
                extraCa += '    d3gds2dp['+str(k)+']['+str(l)+'] = '
                extraCa += '0' if dmin else self.printer.doprint(
                    G.diff(P).diff(s[k]).diff(s[l]))
                extraCa += ';\n'
                extraCa += '    d2gds2['+str(k)+']['+str(l)+'] = '
                extraCa += '0' if dmin else self.printer.doprint(
                    G.diff(s[k]).diff(s[l]))
                extraCa += ';\n'
                for kk in range(0,ns):
                    extraCa += '    d4gds3dp['+str(k)+']['+str(l)+']['+str(kk)+'] = '
                    extraCa += '0' if dmin else self.printer.doprint(
                        G.diff(P).diff(s[k]).diff(s[l]).diff(s[kk]))
                    extraCa += ';\n'
                    extraCa += '    d3gds3['+str(k)+']['+str(l)+']['+str(kk)+'] = '
                    extraCa += '0' if dmin else self.printer.doprint(
                        G.diff(s[k]).diff(s[l]).diff(s[kk]))
                    extraCa += ';\n'
                    for ll in range(0,ns):
                        extraCa += '    d4gds4['+str(k)+']['+str(l)+'][' \
                                 +str(kk)+']['+str(ll)+'] = '
                        extraCa += '0' if dmin else self.printer.doprint(
                            G.diff(s[k]).diff(s[l]).diff(s[kk]).diff(s[ll]))
                        extraCa += ';\n'
                    for ll in range(0,np):
                        extraCa += '    d4gds3dw['+str(k)+']['+str(l)+'][' \
                                 +str(kk)+']['+str(ll)+'] = '
                        extraCa += '0' if dmin else self.printer.doprint(
                            G.diff(s[k]).diff(s[l]).diff(s[kk]).diff(symparam[ll]))
                        extraCa += ';\n'
                for kk in range(0,np):
                    extraCa += '    d4gds2dpdw['+str(k)+']['+str(l)+']['+str(kk)+'] = '
                    extraCa += '0' if dmin else self.printer.doprint(
                        G.diff(P).diff(s[k]).diff(s[l]).diff(symparam[kk]))
                    extraCa += ';\n'
                    extraCa += '    d3gds2dw['+str(k)+']['+str(l)+']['+str(kk)+'] = '
                    extraCa += '0' if dmin else self.printer.doprint(
                        G.diff(s[k]).diff(s[l]).diff(symparam[kk]))
                    extraCa += ';\n'

        j += 1
        if dmin:
            G_jac_list = ['0' for i in range(0, np)]
        else:
            G_jac_list = [self.printer.doprint(
                G_param_jac[j,i]) for i in range(0, np)]

        switch_code_text = ''
        for i in range(0, len(self.params)):
            switch_code_text += '    case ' + str(i) + ':\n'
            switch_code_text += '        result += ' + G_jac_list[i] + ';\n'

            switch_code_text += '        for (i=0; i<'+str(ns)+'; i++) {\n'
            switch_code_text += '            result += ' \
            + '3.0*d4gdsdp2dw[i]['+str(i)+']*dsdp[i] + ' \
            + '    d4gdsdp3[i]*dsdw[i]['+str(i)+'] + ' \
            + '3.0*d3gdsdp2[i]*d2sdpdw[i]['+str(i)+'] + ' \
            + '3.0*d3gdsdpdw[i]['+str(i)+']*d2sdp2[i] + ' \
            + '3.0*d2gdsdp[i]*d3sdp2dw[i]['+str(i)+'] + ' \
            + '    d2gdsdw[i]['+str(i)+']*d3sdp3[i];\n'
            switch_code_text += '            for (j=0; j<'+str(ns)+'; j++) {\n'
            switch_code_text += '                result += ' \
            + '3.0*d4gds2dp2[i][j]*dsdp[i]*dsdw[j]['+str(i)+'] + ' \
            + '3.0*d4gds2dpdw[i][j]['+str(i)+']*dsdp[i]*dsdp[j] + ' \
            + '6.0*d3gds2dp[i][j]*d2sdpdw[i]['+str(i)+']*dsdp[j] + ' \
            + '3.0*d3gds2dp[i][j]*d2sdp2[i]*dsdw[j]['+str(i)+'] + ' \
            + '3.0*d3gds2dw[i][j]['+str(i)+']*d2sdp2[i]*dsdp[j] + ' \
            + '3.0*d2gds2[i][j]*d2sdp2[i]*d2sdpdw[j]['+str(i)+'];\n'
            switch_code_text += '                for (k=0; k<'+str(ns)+'; k++) {\n'
            switch_code_text += '                    result += ' \
            + '3.0*d4gds3dp[i][j][k]*dsdp[i]*dsdp[j]*dsdw[k]['+str(i)+'] + ' \
            + '    d4gds3dw[i][j][k]['+str(i)+']*dsdp[i]*dsdp[j]*dsdp[k] + ' \
            + '3.0*d3gds3[i][j][k]*d2sdp2[i]*dsdp[j]*dsdw[k]['+str(i)+'] + ' \
            + '3.0*d3gds3[i][j][k]*d2sdpdw[i]['+str(i)+']*dsdp[j]*dsdp[k];\n'
            switch_code_text += '                    for (l=0; l<'+str(ns)+'; l++) ' \
            + 'result += d4gds4[i][j][k][l]*dsdp[i]*dsdp[j]*dsdp[k]*dsdw[l][' \
            + str(i)+'];\n'
            switch_code_text += '                }\n'
            switch_code_text += '            }\n'
            switch_code_text += '        }\n'

            switch_code_text += '        break;\n'

        w += solution_calib_template.format(module=self.module, 
            number_components=c, func=self.g_list[j][0], 
            switch_code=switch_code_text, moles_assign=moles_assign_text,
            number_ordering=ns, 
            order_initial_guess=order_initial_guess_text,
            order_assign=extraCb+order_assign_text+extraCa,
            extra_ordering_code=eCB)

        ########
        # Done #
        ########

        if self.verbose:
            print ("... completing extra code block")
        value_params=[self.printer.doprint(symparam[i]) for i in range(0, 
            len(self.params))]
        code_block_one_text = ''
        code_block_two_text = ''
        code_block_three_text = ''
        code_block_four_text = ''
        for i in range(0,len(value_params)):
            code_block_one_text   += ('    (*values)[' + str(i) + '] = ' 
                + value_params[i] + ';\n')
            code_block_two_text   += ('    ' + value_params[i] + ' = values[' 
                + str(i) + '];\n')
            code_block_three_text += ('    case ' + str(i) + ':\n' 
                + '        result = ' + value_params[i] + ';\n' 
                + '        break;\n')
            code_block_four_text  += ('    case ' + str(i) + ':\n' 
                + '        ' + value_params[i] + ' = value;\n' 
                + '        break;\n')

        name_params = self.get_model_param_names()
        unit_params = self.get_model_param_units()
        extra_template = tpl.create_soln_calib_extra_template()
        w += extra_template.format(module=self.module, 
            number_params=len(self.params), 
            names_params=json.dumps(name_params).replace(
                '[', '{').replace(']', '}'), 
            units_params=json.dumps(unit_params).replace(
                '[', '{').replace(']', '}'), 
            code_block_one=code_block_one_text, 
            code_block_two=code_block_two_text, 
            code_block_three=code_block_three_text, 
            code_block_four=code_block_four_text)

        if self.verbose:
            print ("... exiting create_soln_calib_h_file")
        return w

    def create_code_module(self, phase="Feldspar", params={}, endmembers=[],
        identifier=None, prefix="cy", module_type="fast", silent=False, 
        language='C', add_code_to_access_order_paramater=False,
        minimal_deriv_set=False):
        """
        Creates include and code file for a model instance.

        See documentation for superclass.

        Parameters
        ----------
        add_code_to_access_order_paramater : bool
            Generates additional functions to output values of order parameters
            as a function of composition, temperature and pressure
        """
        super(ComplexSolnModel,self).create_code_module(phase, params,
            endmembers, identifier, prefix, module_type, silent, language,
            minimal_deriv_set)
        if add_code_to_access_order_paramater:
            if not silent:
                print ('Generating additional methods for order parameter ' +
                    'access.')
            h_file_name  = phase + '_' + self.module
            h_file_name += '_calib.h' if module_type == 'calib' else '.h'
            c_file_name  = phase + '_' + self.module 
            c_file_name += '_calib.c' if module_type == 'calib' else '.c'
            pyx_file_name = self.module + '.pyx'
            nc = self.nc
            ns = self.ns

            func_signature  = 'void ' + phase + '_' + self.module + '_'
            func_signature += 'order_params(double T, double P, '
            func_signature += 'double n[' + str(nc) + '], '
            func_signature += 'double result[' + str(ns) + '])'

            with open(h_file_name, 'a') as f:
                f.write(func_signature+';\n')
            with open(c_file_name, 'a') as f:
                f.write(func_signature+' {\n')
                fc  = '    double invd2gds2['+str(ns)+']['+str(ns)+'];\n'
                fc += '    int i;\n'
                fc += '    double *s = retrieveGuess(T, P, n);\n'
                fc += '    order_s(T, P, n, s, invd2gds2);\n'
                fc += '    for (i=0; i<'+str(ns)+'; i++) result[i] = s[i];\n'
                fc += '}\n'
                f.write(fc)
            with open(pyx_file_name, 'a') as f:
                fc  = 'cdef extern from "' + phase + '_'
                fc += self.module + '_'
                fc += 'calib.h' if module_type == 'calib' else '.h'
                fc += '":\n'
                fc += '    ' + func_signature + ';\n'
                fc += '\n'
                fc += 'def ' + prefix + '_' + phase + '_' + self.module 
                fc += '_order_params(double t, double p, np_array):\n'
                fc += '    cdef double *m = <double *>malloc(' + str(nc) \
                    + '*sizeof(double))\n'
                fc += '    cdef double *r = <double *>malloc(' + str(ns) \
                    + '*sizeof(double))\n'
                fc += '    for i in range('+str(nc)+'):\n'
                fc += '        m[i] = np_array[i]\n'
                fc += '    ' + phase + '_' + self.module + '_'
                fc += 'order_params(<double> t, <double> p, <double *> m, '
                fc += '<double *> r)\n'
                #fc += '    print(r[0])\n'
                fc += '    r_np_array = np.zeros('+str(ns)+')\n'
                fc += '    for i in range(0,'+str(ns)+'):\n'
                fc += '        r_np_array[i] = r[i]\n'
                fc += '    free(m)\n'
                fc += '    free(r)\n'
                fc += '    return r_np_array\n'
                f.write(fc)
        if not silent:
            print ('Done!')

class f_mu_e(sym.Function):
    """
    A placeholder class that creates a SymPy function that accepts three 
    arguments and returns a value of zero/
    
    This class is used by the SpeciationSolnModel class to initial functions 
    that return excess chemical potentials for an ideal solution.
    
    Extends:
        sym.Function
    
    Variables:
        nargs {number} -- Number of arguments to the function
    """
    nargs = 3
    @classmethod
    def eval(cls, x, y, z):
        return sym.S.Zero

class SpeciationSolnModel(ComplexSolnModel):
    """
    Class creates a model of the thermodynamic properties of a speciated solution.

    Inherits all methods and functionality of the Complex Solution Model class.

    Parameters
    ----------
    nc : int
        Number of thermodynamic components in the solution
    nb : int    
        Number of basis species in the solution (usually, nb == nc)
    ns : int    
        Number of non-basis species in the solution
    Cb : numpy 2-d array
        An array of nb rows and ne columns that maps the composition of basis 
        species to moles of endmember components. Often, this matrix is the 
        identity matrix, but that condition is not required. Cb must be an 
        invertible matrix.
    Cs : numpy 2-d array
        An array of ns rows and ne columns that maps the composition of 
        non-basis species to moles of endmember components.
    R : SymPy Matrix, (ns,nb)
        A SymPy Matrix with ns rows and nb columns containing reaction 
        coefficients that map basis to non-basis species. This matrix need be
        provided only if an expression added to the model utilizes a SymPy
        indexed symbol (e.g., the *r* property of this class).  Otherwise, 
        the quanity is unreferenced. 
    model_type : (str, opt ...)
        Type of speciation model. Default is ('ideal') implying ideal mixing 
        between basis and non-basis species. 
        Acceptable alternatives are:
        
        - ('debye-huckel-limit', []) - Non-ideal interactions are described by 
          the Debye-Hückel limiting law and the mole fraction -> molality 
          conversion term. The second element of the tuple is a list of charges
          associated with the basis plus non-basis species; anions have negative
          values, cations have positive values, neutral species have zero values.  
        
        - ('debye-huckel', [], []) - Non-ideal interactions are described 
          by the full version of Debye-Hückel theory (including the denominator 
          term) and the mole fraction -> molality conversion term. The second 
          element of the tuple is a list of charges associated with the basis 
          plus non-basis species. The third term in the tuple is a list of 
          constants associated with the azero contribution to the Debye-Hückel
          term for each basis and non-basis species.
        
        - ('debye-huckel-ext', [], [], []) - Non-ideal interactions are 
          described by the full version of Debye-Hückel theory (including the 
          denominator term) plus the mole fraction -> molality conversion term 
          and an extended term. The second element of the tuple is a list of 
          charges associated with the basis plus non-basis species. The third 
          term in the tuple is a list of constants associated with the azero 
          contribution to the Debye-Hückel term for each basis and non-basis 
          species. The fourth term is a list of constant coefficients that pre-
          multiply the species molality and define an extended term contribution
          to the excess Gibbs energy. 
    debug_print : bool
        Print debug messages from root solvers in generated C code
    multiroot_method : str
        Method to use in GNU Scientific Library (GSL) solver for homogeneous
        equilibrium speciation. Options are:
        
        gsl_multiroot_fdfsolver_hybridsj
        
        gsl_multiroot_fdfsolver_hybridj
        
        gsl_multiroot_fdfsolver_newton
        
        gsl_multiroot_fdfsolver_gnewton (default)

    Attributes
    ----------
    a_list (inherited from super class)
    b_list
    Cb
    Cs
    debug_print
    expression_parts (inherited from super class)
    g_list (inherited from super class)
    i
    implicit_functions (inherited from super class)
    j
    k
    module (inherited from super class)
    mu (overrides super class)
    multiroot_method
    mu_e
    mu_ex
    mu_s
    n (inherited from super class)
    nb
    nc (inherited from super class)
    ns
    nT (inherited from super class)
    params (inherited from super class)
    printer (inherited from super class)
    r
    s (inherited from super class)
    R
    Rg
    variables (inherited from super class)
 
    Notes
    -----
    The property mu declared by the SimpleSolnModel subclass is modified by 
    this class to return a 1-d SymPy Matrix of symbols for endmember chemical 
    potentials of all solution species, nb+ns, numbered as 1 ... nb+ns

    This class is for creating a representation of the thermodynamic properties 
    of a speciated solution phase.  The principal use of the class is to 
    construct a model symbolically, and to generate code that implements the 
    model for model parameter calibration and thermodynamic property calculation.

    A speciated solution is one that contains implicit concentration variables 
    that need to be determined via solution of conditions of homogeneous 
    equilibrium.  Examples of speciated solutions include aqueous solutions with
    complexes or gas solutions with multiple species.

    Model specification:
    
    - A model may be contructed in a manner similar to the SimpleSolnModel 
      class by adding an expression for the total Gibbs free energy of solution 
      using the add_expression_to_model() method.  This methoid of construction is
      suitable if the specification is compact and does not involve iteration over 
      a common expression. Alternatively,
    
    - A model may be formed by adding a basis species and non-basis species
      contribution, where the non-basis species contribution is specified as an
      implicit sum over a parameterized exprression of fixed mathematical form. 
      A model of this kind is constructed by adding multiple terms using the 
      add_expression_to_model() method. Typically, the terms are:
        - A SymPy expression for the ideal Gibbs free energy of solution 
          contribution of the basis species. This expression should include 
          standard state, and configurational contributions
        
        - A SymPy expression for the ideal Gibbs free energy of solution 
          contribution of the non-basis species. This sympy expression 
          represents a term that involves an implicit summation over non-basis 
          species. This expression may be parameterized using SymPy 
          IndexedBase-class symbols and a sympy function describing the standard 
          state chemical potential of the non-basis species.
        
        - A SymPy expression for the non-ideal Gibbs free energy of solution. 
          The expression may be parameterized using an arbitrary number of SymPy 
          IndexedBase-class symbols.
    """
    def __init__(self, nc=None, nb=None, ns=None, Cb=None, Cs=None, R=None,
        model_type=('ideal'), debug_print=False, 
        multiroot_method="gsl_multiroot_fdfsolver_gnewton"):
        assert Cb is not None, 'Cb must be set.'
        assert Cs is not None, 'Cs must be set.'
        super().__init__(nc=nc, ns=ns, nw=nc+ns)
        nb = nc if nb is None else nb
        ns = 0 if ns is None else ns
        assert nb >= nc, 'Parameter nb must be equal to or exceed nc.'
        assert ns > 0, 'Parameter ns must be greater than zero.'
        self._nb = nb
        self._ns = ns
        self._Rgas = sym.symbols('Rgas')

        # Update standard state chemical potentials to include species
        ss_list = []
        ns_list = []
        T = self.get_symbol_for_t()
        P = self.get_symbol_for_p()

        n = self.n
        for i in range(1, ns+1):
            ss_string = 'mu' + str(i+nc)
            ss_func = sym.Function(ss_string)(T, P)
            ss_list.append(ss_func)
            self._printer.known_functions[ss_string] = \
            '(*endmember[' + str(i+nc-1) + '].mu0)'
        self._mu = self._mu.row_insert(nc, sym.Matrix(ss_list))

        ss_list = []
        for i in range(1, nb+ns+1):
            ss_string = 'mu_e' + str(i)
            ss_list.append(f_mu_e(T, P, n))
        self._mu_ex = sym.Matrix(ss_list)
        self._mu_e = f_mu_e(T, P, n)

        self._Cb = np.copy(Cb)
        self._Cs = np.copy(Cs)
        self._R = R

        # Create basis species mole fractions
        print_subs = []
        self._b_list = []
        nT = self.nT
        s = self.s
        for i in range(0,nc):
            entry = n[i]/nT
            for j in range(0,ns):
                entry -= R[j,i]*s[j]
            self._b_list.append(entry)
            print_subs.append((entry, sym.symbols('b'+str(i+1))))
        self.printer.print_subs = print_subs

        # Create index variables and functions for loop constructions
        self._mu_s = sym.Function('mu_s')(T, P)
        self._i = sym.Idx('i', (1, nb))
        self._j = sym.Idx('j', (1, ns))
        self._k = sym.Idx('k', (1, nb+ns))
        ss_list = []
        for i in range(1, nc+1):
            ss_list.append(sym.IndexedBase('r'+str(i)))
        self._r = sym.Matrix(ss_list)

        self._printer.known_functions['mu_s'] = \
            '(*endmember[' + str(nb) + '+j-1].mu0)'
        self._printer.nBasis = nb
        self._printer.forIndex = 'j'

        self._debug_print = debug_print
        self._multiroot_method = multiroot_method

    @property
    def b_list(self):
        """
        Symbolic definitions of basis species

        A list of SymPy expressions that define the mole fractions of basis
        species (the elements of b_list) in terms of moles of endmember 
        components and mole fractions of non-basis species (elements of s)
        
        Returns:
            list of SymPy expressions
        """
        return self._b_list

    @property
    def Cb(self):
        """
        Stoichiometric array mapping basis species to elements
        
        An array of nb rows and ne columns that maps the composition of basis 
        species to moles of endmember components. Often, this matrix is the 
        identity matrix, but that condition is not required. Cb must be an 
        invertible matrix.
        
        Returns:
            Numpy 2d array
        """
        return self._Cb

    @property
    def Cs(self):
        """
        Stoichiometric array mapping non-basis species to elements
        
        [description]
        
        Returns:
            [type] -- [description]
        """
        return self._Cs

    @property
    def debug_print(self):
        """
        Print debug messages from root solvers in generated C code
    
        Passed to code generation routines to flag generation of output related
        to the functioning of root finding and related numerical methods for 
        speciation.

        Returns:
            bool
        """
        return self._debug_print

    @property
    def i(self):
        """
        Index variable for looping over basis species
        
        A SymPy instance of the class Idx that is used to loop over terms in
        the Gibbs free energy that describe properties of basis species. The 
        range of the index is 1 to nb+1. Use this symbol as an index for 
        summation or product terms in the potential.
        
        Returns:
            SymPy.Idx
        """
        return self._j

    @property
    def j(self):
        """
        Index variable for looping over non-basis species
        
        A SymPy instance of the class Idx that is used to loop over terms in
        the Gibbs free energy that describe properties of non-basis species. 
        The range of the index is 1 to ns+1. Use this symbol as an index for 
        summation or product terms in the potential.
        
        Returns:
            SymPy.Idx
        """
        return self._j

    @property
    def k(self):
        """
        Index variable for looping over all species
        
        A SymPy instance of the class Idx that is used to loop over terms in
        the Gibbs free energy that describe properties of species. The range
        of the index is 1 to nb+ns+1. Use this symbol as an index for summation
        or product terms in the potential.
        
        Returns:
            SymPy.Idx
        """
        return self._j

    @property
    def mu(self):
        """
        1-d matrix of SymPy symbols for the chemical potentials of basis and 
        non-basis species in the model

        Notes
        -----
        The property definition is changed from that of the underlying Simple 
        Solution superclass
                
        Returns
        -------
        SymPy Matrix object (sympy.Matrix)
        """
        return self._mu 

    @property
    def multiroot_method(self):
        """
        Method to use in GNU Scientific Library (GSL) multiroot solver 
        
        The gsl_multiroot_fdfsolver_type. Options are:  
        
        "gsl_multiroot_fdfsolver_hybridsj"
        
        "gsl_multiroot_fdfsolver_hybridj"
        
        "gsl_multiroot_fdfsolver_newton"
        
        "gsl_multiroot_fdfsolver_gnewton" (default)
        
        Returns:
            str
        """
        return self._multiroot_method

    @property
    def mu_e(self):
        """
        Generic excess chemical potential of a non-basis species
        
        A SymPy instance of the class Function that is used to represent the
        value of the excess chemical potential of a non-basis species in an 
        expression that loops species properties.  In such a loop, values
        of mu_e are replaced with specific entries of mu_ex[...] when the code 
        is printed. The instance is a function of T, P, and n.

        By default, for an ideal solution this function returns a value of zero.
        
        Returns:
            SymPy function
        """
        return self._mu_e

    @property
    def mu_ex(self):
        """
        A SymPy matrix of SymPy functions that calculate the excess chemical
        potential of each species in solution

        The functions must have arguments T, P and n, where T is temperature, 
        P is pressure and n is a vector of mole numbers of basis species in 
        solution.  The matrix has shape (nb+ns, 1)

        By default, for an ideal solution each function returns a value of zero.
        
        Returns:
            SymPy Matrix
        """
        return self._mu_ex

    @property
    def mu_s(self):
        """
        Generic standard state chemical potenial of a non-basis species
        
        A SymPy instance of the class Function that is used to represent the
        value of the standard state Gibbs free energy of a non-basis species
        in an expression that loops species properties.  In such a loop, values
        of mu_s are replaced with specific entries of mu[...] when the code is
        printed. The instance is a function of T and P.
        
        Returns:
            SymPy.Function
        """
        return self._mu_s

    @property
    def nb(self):
        """
        Number of basis species
        
        The number of basis spcies is generally equal to the number of 
        thermodynamic components (nc) required to describe the solution phase.
        
        Returns:
            int
        """
        return self._nb

    @property
    def ns(self):
        """
        Number of non-basis species
        
        The number of non-basis species is strictly positive but unbounded.
        Concentrations of non-basis species are determined by solution of mass
        action expressions corresponding to conditions of homogeneous 
        equilibrium.
        
        Returns:
            int
        """
        return self._ns

    @property
    def r(self):
        """
        Vector of SymPy symbols for indexed stoichiometric reaction coefficients 

        A SymPy Matrix of length *nc* whose elements are instances of the 
        SymPy class IndexedBase. These symbols are used to parameterize model
        Gibbs free energy expressions that require coefficients that map the
        stoichiometry of basis species to that of non-basis species; i.e.,
        (moles of non-basis species) = r[0]*(moles of 1st basis species) + ...
        r[nc-1]*(moles of nc-1 basis species).  Numerical values for these 
        symbols will be assigned by letting r[...] denote rows of the numpy 
        matrix R, where R is given by Cs (Cb)^-1. Cb is an array of nb rows and 
        ne columns (ne is the number of elements; generally ne = nb) that maps 
        the composition of basis species to moles of endmember components. Cs is
        an array of ns rows and ne columns that maps the composition of 
        non-basis species to moles of endmember components. Numpy matrices Cb 
        and Cs are provided to the class at the code generation step when 
        calling the create_code_module() method.

        Returns:
            SymPy.Matrix -- length nb
        """
        return self._r

    @property
    def R(self):
        """
        SymPy Matrix of stoichiometric reaction coefficients
        
        A SymPy Matrix of shape (ns, nb) whose elements are constants that map
        the stoichiometry of basis species to that of non-basis species. The 
        constants are SymPy compatible rational numbers. The columns of R are
        assigned to the elements of r at code generation.  See the description 
        of the r property for further explanation as to how R is defined.

        R need be specified only if one or more of the expressions added to the
        model contain a SymPy IndexedBase variable.
        
        Returns:
            SymPy Matrix (ns,nb)
        """
        return self._R

    @property
    def Rgas(self):
        """
        Sympy symbol for the Universal gas constant
        
        The universal gas constant used as a parameter in constructing expresions
        for the Gibbs free energy. Must be given units and a parameter value.
        
        Returns:
            SymPy Symbol
        """
        return self._Rgas

    def create_ordering_code(self, moles_assign_text, order_assign_text):
        """
        Generates ordering function code for a speciation solution model.

        Parameters
        ----------
        moles_assign_text : str
            Text string containing C code for assigning moles of each 
            component (n1, n2, ..., nc) from a vector of moles, n.

        order_assign_text : str
            Text string containing C code for assigning values of each
            ordering variable (s1, s2, ..., ss) from a vector, s

        Returns
        -------
        output : str
            String of code that implements ordering functions and derivatives.

        Notes
        -----
        The user does not normally call this function directly.

        """
        nc = self.nc
        n  = self.n
        ns = self.ns
        s  = self.s
        sFunc = self._ordering_functions
        T = self.get_symbol_for_t()
        P = self.get_symbol_for_p()

        for i in range(0,ns):
            moles_assign_text += '    double b'+str(i+1)+' = b['+str(i)+'];\n'
        moles_assign_text += order_assign_text
        
        cb_5  = '    double d2gdnds['+str(nc)+']['+str(ns)+'];\n'
        cb_5 += moles_assign_text
        cb_6  = '    double d2gdsdt['+str(ns)+'];\n'
        cb_6 += moles_assign_text
        cb_7  = '    double d2gdsdp['+str(ns)+'];\n'
        cb_7 += moles_assign_text
        cb_8  = '    double d2gdnds['+str(nc)+']['+str(ns)+'],d3gdn2ds['+ \
                str(nc)+']['+str(nc)+']['+str(ns)+'],d3gdnds2['+str(nc)+']['+ \
                str(ns)+']['+str(ns)+'],d3gds3['+str(ns)+']['+str(ns)+']['+ \
                str(ns)+'];\n'
        cb_8 += moles_assign_text
        cb_9  = '    double d2gdnds['+str(nc)+']['+str(ns)+'],d3gdndsdt['+ \
                str(nc)+']['+str(ns)+'], d3gdnds2['+str(nc)+']['+str(ns)+']['+ \
                str(ns)+'],d2gdsdt['+str(ns)+'],d3gds2dt['+str(ns)+']['+ \
                str(ns)+'],d3gds3['+str(ns)+']['+str(ns)+']['+str(ns)+'];\n'
        cb_9 += moles_assign_text

        cb_10 = '    double d2gdnds['+str(nc)+']['+str(ns)+'],d3gdndsdp['+ \
                str(nc)+']['+str(ns)+'],d3gdnds2['+str(nc)+']['+str(ns)+']['+ \
                str(ns)+'],d2gdsdp['+str(ns)+'],d3gds3['+str(ns)+']['+str(ns)+ \
                ']['+str(ns)+'],d3gds2dp['+str(ns)+']['+str(ns)+'];\n'
        cb_10 += moles_assign_text
        cb_11  = '    double d2gdsdt['+str(ns)+'],d3gdsdt2['+str(ns)+ \
                 '],d3gds2dt['+ str(ns)+']['+str(ns)+'],d3gds3['+str(ns)+ \
                 ']['+str(ns)+']['+ str(ns)+'];\n'
        cb_11 += moles_assign_text
        cb_12  = '    double d2gdsdt['+str(ns)+'],d2gdsdp['+str(ns)+ \
                 '],d3gdsdtdp['+str(ns)+'],d3gds2dt['+str(ns)+']['+str(ns)+ \
                 '],d3gds2dp['+str(ns)+']['+str(ns)+'],d3gds3['+str(ns)+']['+ \
                 str(ns)+']['+str(ns)+'];\n'
        cb_12 += moles_assign_text
        cb_13  = '    double d2gdsdp['+str(ns)+'],d3gdsdp2['+str(ns)+ \
                 '],d3gds2dp['+str(ns)+']['+str(ns)+'],d3gds3['+str(ns)+']['+ \
                 str(ns)+ ']['+str(ns)+'];\n'
        cb_13 += moles_assign_text

        dmin = self.de_minimis
        for i in range(0,nc):
            for j in range(0,ns):
                # d2gdnds[{NC}][{NS}]
                a = '    d2gdnds[' + str(i) + '][' + str(j) + '] = '
                a += self.printer.doprint(
                    sym.diff(sFunc.function[j], n[i])) + ';\n'
                cb_5  += a
                cb_8  += a
                cb_9  += a
                cb_10 += a
                # d3gdndsdt[{NC}][{NS}]
                a = '    d3gdndsdt[' + str(i) + '][' + str(j) + '] = '
                a += '0' if dmin else self.printer.doprint(
                    sym.diff(sFunc.function[j], n[i], T))
                a += ';\n'
                cb_9 += a
                # d3gdndsdp[{NC}][{NS}]
                a = '    d3gdndsdp[' + str(i) + '][' + str(j) + '] = '
                a += '0' if dmin else self.printer.doprint(
                    sym.diff(sFunc.function[j], n[i], P))
                a += ';\n'
                cb_10 += a
                for k in range(i,nc):
                    # d3gdn2ds[{NC}][{NC}][{NS}]
                    a = '    d3gdn2ds[' + str(i) + '][' + str(k) + '][' \
                      + str(j) + '] = '
                    a += '0' if dmin else self.printer.doprint(
                        sym.diff(sFunc.function[j], n[i], n[k]))
                    a += ';\n'
                    cb_8 += a
                    if k > i:
                        cb_8 += '    d3gdn2ds[' + str(k) + '][' + str(i) + \
                                '][' + str(j) + '] = d3gdn2ds[' + str(i) + \
                                '][' + str(k) + '][' + str(j) + '];\n'
                for k in range(j,ns):
                    # d3gdnds2[{NC}][{NS}][{NS}]
                    a = '    d3gdnds2[' + str(i) + '][' + str(j) + '][' \
                      + str(k) + '] = '
                    a += '0' if dmin else self.printer.doprint(
                        sym.diff(sFunc.function[j], n[i], s[k]))
                    a += ';\n'
                    cb_8  += a
                    cb_9  += a
                    cb_10 += a
                    if k > j:
                        cb_8 += '    d3gdnds2[' + str(i) + '][' + str(k) + \
                                '][' + str(j) + '] = d3gdnds2[' + str(i) + \
                                '][' + str(j) + '][' + str(k) + '];\n'
                        cb_9 += '    d3gdnds2[' + str(i) + '][' + str(k) + \
                                '][' + str(j) + '] = d3gdnds2[' + str(i) + \
                                '][' + str(j) + '][' + str(k) + '];\n'
                        cb_10 += '    d3gdnds2[' + str(i) + '][' + str(k) + \
                                 '][' + str(j) + '] = d3gdnds2[' + str(i) + \
                                 '][' + str(j) + '][' + str(k) + '];\n'

        for i in range(0,ns):
            # d2gdsdt[{NS}]
            a = '    d2gdsdt[' + str(i) + '] = '
            a += self.printer.doprint(sym.diff(sFunc.function[i], T)) + ';\n'
            cb_6  += a
            cb_9  += a
            cb_11 += a
            cb_12 += a
            # d2gdsdp[{NS}]
            a = '    d2gdsdp[' + str(i) + '] = '
            a += self.printer.doprint(sym.diff(sFunc.function[i], P)) + ';\n'
            cb_7  += a
            cb_10 += a
            cb_12 += a
            cb_13 += a
            # d3gdsdt2[{NS}]
            a = '    d3gdsdt2[' + str(i) + '] = '
            a += '0' if dmin else self.printer.doprint(
                sym.diff(sFunc.function[i], T, 2))
            a += ';\n'
            cb_11 += a
            # d3gdsdtdp[{NS}]
            a = '    d3gdsdtdp[' + str(i) + '] = '
            a += '0' if dmin else self.printer.doprint(
                sym.diff(sFunc.function[i], T, P))
            a += ';\n'
            cb_12 += a
            # d3gdsdp2[{NS}]
            a = '    d3gdsdp2[' + str(i) + '] = '
            a += '0' if dmin else self.printer.doprint(
                sym.diff(sFunc.function[i], P, 2))
            a += ';\n'
            cb_13 += a
            for j in range(i,ns):
                # d3gds2dt[{NS}][{NS}]
                a  = '    d3gds2dt[' + str(i) + '][' + str(j) + '] = '
                a += '0' if dmin else self.printer.doprint(
                    sym.diff(sFunc.function[i], s[j], T))
                a += ';\n'
                cb_9  += a
                cb_11 += a
                cb_12 += a
                if j > i:
                    cb_9 += '    d3gds2dt[' + str(j) + '][' + str(i) + '] = ' \
                          + 'd3gds2dt[' + str(i) + '][' + str(j) + '];\n'
                    cb_11 += '    d3gds2dt[' + str(j) + '][' + str(i) + '] = ' \
                           + 'd3gds2dt[' + str(i) + '][' + str(j) + '];\n'
                    cb_12 += '    d3gds2dt[' + str(j) + '][' + str(i) + '] = ' \
                           + 'd3gds2dt[' + str(i) + '][' + str(j) + '];\n'
                # d3gds2dp[{NS}][{NS}]
                a = '    d3gds2dp[' + str(i) + '][' + str(j) + '] = '
                a += '0' if dmin else self.printer.doprint(
                    sym.diff(sFunc.function[i], s[j], P))
                a += ';\n'
                cb_10 += a
                cb_12 += a
                cb_13 += a
                if j > i:
                    cb_10 += '    d3gds2dp[' + str(j) + '][' + str(i) + '] = ' \
                           + 'd3gds2dp[' + str(i) + '][' + str(j) + '];\n'
                    cb_12 += '    d3gds2dp[' + str(j) + '][' + str(i) + '] = ' \
                           + 'd3gds2dp[' + str(i) + '][' + str(j) + '];\n'
                    cb_13 += '    d3gds2dp[' + str(j) + '][' + str(i) + '] = ' \
                           + 'd3gds2dp[' + str(i) + '][' + str(j) + '];\n'
                for k in range(j,ns):
                    # d3gds3[{NS}][{NS}][{NS}]
                    a  = '    d3gds3[' + str(i) + '][' + str(j) + '][' + \
                         str(k) + '] = '
                    a += '0' if dmin else self.printer.doprint(
                        sym.diff(sFunc.function[i], s[j], s[k]))
                    a += ';\n'
                    a += self.exchange_indices(i,j,k,'d3gds3')
                    cb_8  += a
                    cb_9  += a
                    cb_10 += a
                    cb_11 += a
                    cb_12 += a
                    cb_13 += a
        return tpl.create_speciation_ordering_code_template().format(
            NS=ns, NC=nc,
            ORDER_CODE_BLOCK_FIVE=cb_5,
            ORDER_CODE_BLOCK_SIX=cb_6,
            ORDER_CODE_BLOCK_SEVEN=cb_7,
            ORDER_CODE_BLOCK_EIGHT=cb_8,
            ORDER_CODE_BLOCK_NINE=cb_9,
            ORDER_CODE_BLOCK_TEN=cb_10,
            ORDER_CODE_BLOCK_ELEVEN=cb_11,
            ORDER_CODE_BLOCK_TWELVE=cb_12,
            ORDER_CODE_BLOCK_THIRTEEN=cb_13)

    def create_code_module(self, phase="IdealGas", params={}, endmembers=[],
        identifier=None, prefix="cy", module_type="fast", 
        silent=False, language='C', minimal_deriv_set=False, debug=0):
        """
        Creates include and code file for a model instance.

        Parameters
        ----------
        phase : str
            Model instance title (e.g., phase name).  Used to name the generated 
            function. Cannot contain blanks or special characters; underscore 
            ("_") is permitted. Convention capitalizes the first letter and letter 
            following a "_" character.
        params : dict
            Parameter values for the model instance.
            The keys of this dictionary are validated against parameter symbols 
            stored for the model.
        endmembers : list of str
            A list of prefixes for standard state property functions for the 
            endmember components of this solution. e.g., "Albite_berman" will
            be used to call functions with names like Albite_berman_g(...)
            If the standard state property routines are coded by the 
            StdStateModel Class in this module, all required functions will be
            generated and they will automatically be compliant with this naming 
            convention. There must be nc+ns entries in this list, i.e.,
            list elements must include basis + non-basis species in an order
            that is internally consistent with specification of the G function.
        identifier : str
            A unique identifier for the model instance.
            Defaults to date/time when module is created (rounded to the second).
        prefix : str
            Prefix to function names for python bindings, e.g., 
            {prefix}_{phase}_{module}_g(T,P)
        module_type : str
            Generate code that executes "fast", but does not expose hooks for 
            model parameter calibration. Alternately, generate code suitable for 
            "calib"ration of parameters in the model, which executes more slowly 
            and exposes additional functions that allow setting of parameters 
            and generation of derivatives of thermodynamic functions with 
            respect to model parameters. 
        silent : bool
            Do not print status messages.
        language : str
            Language syntax for generated code, ("C" is the C99 programming 
            language.)
        minimal_deriv_set : bool
            Generate a minimal set of compositional derivatives: dgdn, d2gdndt, 
            d2gdndp, d3gdndt2, d3gdndtdp, d3gdndp2, d4gdndt3, d4gdndt2dp, 
            d4gdndtdp2, d4gdndp3, d2gdn2, d3gdn2dt, d3gdn2dp, d3gdn3.  This is 
            the subset of derivatives currently required for solution phases 
            that are imported into the phases module. Remaining derivatives are
            returned with values of zero.
        debug : int
            Level of debugging output generated by module:
            0 : None
            1 : Informational messages
            2 : Informational plus normal debugging messages
            3 : Verbose informational and debugging messages

        Returns
        -------
        result : Boolean
                 True if module is succesfully generated, False if some error 
                 occurred 

        Notes
        -----
        This method overrides that of the super class.
        """
        nc = self.nc
        nb = self.nb
        ns = self.ns

        # Generate Simple Solution code files
        success = super(SpeciationSolnModel, self).create_code_module(
            phase=phase, 
            params=params, 
            endmembers=endmembers[:nc],
            identifier=identifier, 
            prefix=prefix, 
            module_type=module_type, 
            silent=silent, 
            language=language,
            minimal_deriv_set=minimal_deriv_set)

        #
        # Updates to fast calc_h file
        if not silent:
            print ("Updating code blocks for speciation code.")

        # Use this text to search for code injection points
        x = '#include <math.h>\n\n'
        # Generate code to inject after these lines
        Cb_text = ''
        Cs_text = ''
        for j in range(0,nc):
            for i in range(0,nb):
                Cb_text += '        gsl_matrix_set(params.Cb, ' + str(i) + \
                    ',' + str(j) + ',' + str(self.Cb[i,j]) + ');\n'
            for i in range(0,ns):
                Cs_text += '        gsl_matrix_set(params.Cs, ' + str(i) + \
                    ',' + str(j) + ',' + str(self.Cs[i,j]) + ');\n'

        speciation_template = tpl.create_speciation_code_template(language)
        # generate fill_invd2gds2 for speciation insert
        moles_assign_text  = '    double T = t;\n    double P = p;\n'
        for i in range(0,nc):
            moles_assign_text += ('    double n' + str(i+1) + ' = e[' 
                + str(i) + '];')
            moles_assign_text += '\n'
            moles_assign_text += ('    double b' + str(i+1) + ' = b[' 
                + str(i) + '];')
            moles_assign_text += '\n'
        invd2gds2_text = moles_assign_text + '    double '
        separator = ''
        for i in range(0,ns):
            invd2gds2_text += separator + 's' + str(i+1)
            separator = ', '
        invd2gds2_text += ';\n'
        for i in range(0,ns):
            invd2gds2_text += '    s' + str(i+1) + ' = s[' + str(i) + '];\n'
        s = self.s
        sFunc = self._ordering_functions
        for i in range(0,ns):
            # invd2gds2[0][0] = D2GDS0S0;
            for j in range(i,ns):
                invd2gds2_text += '    invd2gds2[' 
                invd2gds2_text += str(i) + '][' + str(j) + '] = '
                a = sym.diff(sFunc.function[i], s[j], 1)
                invd2gds2_text += self.printer.doprint(a) + ';\n'
                if i < j:
                    invd2gds2_text += '    invd2gds2[' + str(j) + '][' 
                    invd2gds2_text += str(i) + '] = invd2gds2[' + str(i) 
                    invd2gds2_text += '][' + str(j) + '];\n'
        invd2gds2_text += '\n'
        if ns == 1:
            invd2gds2_text += '    invd2gds2[0][0] = 1.0/invd2gds2[0][0];\n'
        else:
            invd2gds2_text += '    gaussj(invd2gds2);'
        code_to_inject = speciation_template.format(number_components=nc, 
            number_non_basis=ns, 
            gsl_multiroot_method=self.multiroot_method,
            fill_Cb=Cb_text, 
            fill_Cs=Cs_text,
            fill_invd2gds2=invd2gds2_text)

        # Retreive the file to modify
        module = self.module
        with open(module+'_calc.h', 'r') as f:
            fold = f.read()
            f.close()
        fnew = fold.replace(x, x + code_to_inject)
        fold = fnew

        # Use this text to search for code injection points
        x = '    order_s(T, P, n, s, invd2gds2);\n'
        # Add code to call speciation routine
        code_to_add = """    double b[{nc}];
    speciate(T, P, n, b, s, invd2gds2, {debug});
""".format(debug=(1 if debug > 0 else 0),nc=nc)
        for i in range(0,nc):
            code_to_add += ('    double b' + str(i+1) + ' = b[' 
                + str(i) + '];')
            code_to_add += '\n'
        fnew = fold.replace(x, code_to_add)
        fold = fnew

        #order parameter derivatives
        fnew = fold.replace('order_dsdt(T, P, n, s, invd2gds2, dsdt)',
                            'order_dsdt(T, P, n, b, s, invd2gds2, dsdt)')
        fold = fnew
        fnew = fold.replace('order_dsdp(T, P, n, s, invd2gds2, dsdp)',
                            'order_dsdp(T, P, n, b, s, invd2gds2, dsdp)')
        fold = fnew
        fnew = fold.replace('order_dsdn(T, P, n, s, invd2gds2, dsdn)',
                            'order_dsdn(T, P, n, b, s, invd2gds2, dsdn)')
        fold = fnew
        fnew = fold.replace('order_d2sdt2(T, P, n, s, invd2gds2, d2sdt2)',
                            'order_d2sdt2(T, P, n, b, s, invd2gds2, d2sdt2)')
        fold = fnew
        fnew = fold.replace('order_d2sdtdp(T, P, n, s, invd2gds2, d2sdtdp)',
                            'order_d2sdtdp(T, P, n, b, s, invd2gds2, d2sdtdp)')
        fold = fnew
        fnew = fold.replace('order_d2sdp2(T, P, n, s, invd2gds2, d2sdp2)',
                            'order_d2sdp2(T, P, n, b, s, invd2gds2, d2sdp2)')
        fold = fnew
        fnew = fold.replace('order_d2sdndt(T, P, n, s, invd2gds2, d2sdndt)',
                            'order_d2sdndt(T, P, n, b, s, invd2gds2, d2sdndt)')
        fold = fnew
        fnew = fold.replace('order_d2sdndp(T, P, n, s, invd2gds2, d2sdndp)',
                            'order_d2sdndp(T, P, n, b, s, invd2gds2, d2sdndp)')
        fold = fnew
        fnew = fold.replace('order_d2sdn2(T, P, n, s, invd2gds2, d2sdn2)',
                            'order_d2sdn2(T, P, n, b, s, invd2gds2, d2sdn2)')
        fold = fnew

        # Write modified code 
        with open(module+'_calc.h', 'w') as f:
            f.write(fnew)
            f.close()

        #
        # Update the _calib.h files if it was generated
        if module_type == 'calib':
            with open(module+'_calib.h', 'r') as f:
                fold = f.read()
                f.close()
            fnew = fold.replace(x, code_to_add)
            fold = fnew
            #order parameter derivatives
            fnew = fold.replace('order_dsdt(T, P, n, s, invd2gds2, dsdt)',
                                'order_dsdt(T, P, n, b, s, invd2gds2, dsdt)')
            fold = fnew
            fnew = fold.replace('order_dsdp(T, P, n, s, invd2gds2, dsdp)',
                                'order_dsdp(T, P, n, b, s, invd2gds2, dsdp)')
            fold = fnew
            fnew = fold.replace('order_dsdn(T, P, n, s, invd2gds2, dsdn)',
                                'order_dsdn(T, P, n, b, s, invd2gds2, dsdn)')
            fold = fnew
            fnew = fold.replace('order_dsdw(T, P, n, s, invd2gds2, dsdw)',
                                'order_dsdw(T, P, n, b, s, invd2gds2, dsdw)')
            fold = fnew
            fnew = fold.replace('order_d2sdt2(T, P, n, s, invd2gds2, d2sdt2)',
                                'order_d2sdt2(T, P, n, b, s, invd2gds2, d2sdt2)')
            fold = fnew
            fnew = fold.replace('order_d2sdtdp(T, P, n, s, invd2gds2, d2sdtdp)',
                                'order_d2sdtdp(T, P, n, b, s, invd2gds2, d2sdtdp)')
            fold = fnew
            fnew = fold.replace('order_d2sdtdw(T, P, n, s, invd2gds2, d2sdtdw)',
                                'order_d2sdtdw(T, P, n, b, s, invd2gds2, d2sdtdw)')
            fold = fnew
            fnew = fold.replace('order_d2sdp2(T, P, n, s, invd2gds2, d2sdp2)',
                                'order_d2sdp2(T, P, n, b, s, invd2gds2, d2sdp2)')
            fold = fnew
            fnew = fold.replace('order_d2sdpdw(T, P, n, s, invd2gds2, d2sdpdw)',
                                'order_d2sdpdw(T, P, n, b, s, invd2gds2, d2sdpdw)')
            fold = fnew
            fnew = fold.replace('order_d2sdndt(T, P, n, s, invd2gds2, d2sdndt)',
                                'order_d2sdndt(T, P, n, b, s, invd2gds2, d2sdndt)')
            fold = fnew
            fnew = fold.replace('order_d2sdndp(T, P, n, s, invd2gds2, d2sdndp)',
                                'order_d2sdndp(T, P, n, b, s, invd2gds2, d2sdndp)')
            fold = fnew
            fnew = fold.replace('order_d2sdn2(T, P, n, s, invd2gds2, d2sdn2)',
                                'order_d2sdn2(T, P, n, b, s, invd2gds2, d2sdn2)')
            fold = fnew
            fnew = fold.replace('order_d3sdt3(T, P, n, s, invd2gds2, d3sdt3)',
                                'order_d3sdt3(T, P, n, b, s, invd2gds2, d3sdt3)')
            fold = fnew
            fnew = fold.replace('order_d3sdt2dp(T, P, n, s, invd2gds2, d3sdt2dp)',
                                'order_d3sdt2dp(T, P, n, b, s, invd2gds2, d3sdt2dp)')
            fold = fnew
            fnew = fold.replace('order_d3sdt2dw(T, P, n, s, invd2gds2, d3sdt2dw)',
                                'order_d3sdt2dw(T, P, n, b, s, invd2gds2, d3sdt2dw)')
            fold = fnew
            fnew = fold.replace('order_d3sdtdp2(T, P, n, s, invd2gds2, d3sdtdp2)',
                                'order_d3sdtdp2(T, P, n, b, s, invd2gds2, d3sdtdp2)')
            fold = fnew
            fnew = fold.replace('order_d3sdtdpdw(T, P, n, s, invd2gds2, d3sdtdpdw)',
                                'order_d3sdtdpdw(T, P, n, b, s, invd2gds2, d3sdtdpdw)')
            fold = fnew
            fnew = fold.replace('order_d3sdp3(T, P, n, s, invd2gds2, d3sdp3)',
                                'order_d3sdp3(T, P, n, b, s, invd2gds2, d3sdp3)')
            fold = fnew
            fnew = fold.replace('order_d3sdp2dw(T, P, n, s, invd2gds2, d3sdp2dw)',
                                'order_d3sdp2dw(T, P, n, b, s, invd2gds2, d3sdp2dw)')
            fold = fnew
            with open(module+'_calib.h', 'w') as f:
                f.write(fnew)
                f.close()

        #
        # Updates to calc_c file
        if not silent:
            print ("Updating code blocks for standard state properties.")

        code_block_find_text = ''
        code_block_repl_text = ''
        incl_block_find_text = ''
        incl_block_repl_text = ''
        ss_xx = '' if module_type == 'fast' else '_calib'
        for i,x in enumerate(endmembers):
            code_block_text  = '  {\n'
            code_block_text += '    ' + x + ss_xx + '_name,\n'
            code_block_text += '    ' + x + ss_xx + '_formula,\n'
            code_block_text += '    ' + x + ss_xx + '_mw,\n'
            code_block_text += '    ' + x + ss_xx + '_elements,\n'
            code_block_text += '    ' + x + ss_xx + '_g,\n'
            code_block_text += '    ' + x + ss_xx + '_dgdt,\n'
            code_block_text += '    ' + x + ss_xx + '_dgdp,\n'
            code_block_text += '    ' + x + ss_xx + '_d2gdt2,\n'
            code_block_text += '    ' + x + ss_xx + '_d2gdtdp,\n'
            code_block_text += '    ' + x + ss_xx + '_d2gdp2,\n'
            code_block_text += '    ' + x + ss_xx + '_d3gdt3,\n'
            code_block_text += '    ' + x + ss_xx + '_d3gdt2dp,\n'
            code_block_text += '    ' + x + ss_xx + '_d3gdtdp2,\n'
            code_block_text += '    ' + x + ss_xx + '_d3gdp3\n'
            code_block_text += '  },\n'
            if module_type == 'fast':
                incl_block_text = ('#include "' + x + '_calc.h"\n')
            elif module_type == 'calib':
                incl_block_text = ('#include "' + x + '_calib.h"\n')
            if i < nc:
                code_block_find_text += code_block_text
                incl_block_find_text += incl_block_text
            code_block_repl_text += code_block_text
            incl_block_repl_text += incl_block_text

        # Open file and replace code
        with open(phase + '_' + module + ss_xx + '.c', 'r') as f:
            fold = f.read()
            f.close()
        fnew = fold.replace(code_block_find_text, code_block_repl_text)
        fold = fnew
        fnew = fold.replace(incl_block_find_text, incl_block_repl_text)
        fold = fnew
        x   = 'static int nc = (sizeof endmember / sizeof(struct _endmembers));'
        rx  = 'static int nc = ' + str(nc) + ';\n' 
        rx += 'static int ns = ' + str(ns) + ';\n'
        fnew = fold.replace(x, rx)
        fold=fnew

        # Add R matrix constants for loop code
        add_R_code = False
        for exp in self.expression_parts:
            add_R_code |= exp.expression.has(sym.IndexedBase)
        if add_R_code:
            assert self.R is not None, \
                'R must be specified because an expression contains an ' + \
                'IndexedBase SymPy symbol.'
            assert self.R.shape == (ns,nb), \
                'R has shape ' + str(self.R.shape) + \
                'when it should have shape ' + str((ns,nb)) + '.'
            x = 'static const double R=8.3143;'
            rx = ''
            for i in range(1,nb+1):
                rx += 'static const double r'
                rx += str(i) + '[' + str(ns+1) + '] = {'
                delim = ' 0, '
                for j in range(0,ns):
                    rx += delim + self.printer.doprint(self.R[j,i-1])
                    delim = ', '
                rx += ' };\n'
            fnew = fold.replace(x, rx)

        # Write modified code 
        with open(phase + '_' + module + ss_xx + '.c', 'w') as f:
            f.write(fnew)
            f.close()

        #
        # Updates to pyxbld file
        if not silent:
            print ("Updating _pyxbld file.")
        with open(module + '.pyxbld', 'r') as f:
            fold = f.read()
            f.close()
        x   = "=['-O3']"
        rx  = "=['-O3'],\n"
        rx += "                     libraries=['gsl', 'speciation'],\n"
        rx += "                     library_dirs=['/usr/local/lib'],\n"
        rx += "                     runtime_library_dirs=['/usr/local/lib']"
        fnew = fold.replace(x, rx)
        fold = fnew

        c_file_name = phase + '_' + module 
        if module_type == 'calib':
            c_file_name += '_calib.c'
        else:
            c_file_name += '_calc.c'
        new_pyx_file_list = "'" + c_file_name + "'"
        old_pyx_file_list = "'" + c_file_name + "'"
        for i,x in enumerate(endmembers):
            new_pyx_file_list += ','
            old_pyx_file_list += ',' if i < nc else ''
            if module_type == 'fast':
                new_pyx_file_list += "'" + x + "_calc.c'"
                old_pyx_file_list += "'" + x + "_calc.c'" if i < nc else ''
            elif module_type == 'calib':
                new_pyx_file_list += "'" + x + "_calib.c'"
                old_pyx_file_list += "'" + x + "_calib.c'" if i < nc else ''
        fnew = fold.replace(old_pyx_file_list, new_pyx_file_list)

        with open(module + '.pyxbld', 'w') as f:
            f.write(fnew)
            f.close()

        if not silent:
            print ("Updates completed.")

        return success

def import_coder_phase(module_name):
    from pyximport import install
    from importlib import import_module
    install()
    coder_phase = import_module(module_name)
    return coder_phase
