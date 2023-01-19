# Notes
* These are notes in attemts to create an auto-installing conda package for thermoengine
* All info for environment is stored in environment.yml file in top level directory

"""
conda env update -f environment.yml
"""

"""
conda env remove -n thermoengine
conda env create -n thermoengine -f environment.yml
conda activate thermoengine
"""


* tests all pass except 1 in the Equilibrate/7a-K-Quartz-Equilibrate.ipynb notebook


"""
======================================== test session starts =========================================
platform darwin -- Python 3.8.5, pytest-6.2.2, py-1.10.0, pluggy-0.13.1
rootdir: /Users/aswolf/Documents/projects/ENKI/ThermoEngineASW/Notebooks
plugins: cov-2.11.1, nbval-0.9.6
collected 18 items

Equilibrate/7a-K-Quartz-Equilibrate.ipynb ................F.                                   [100%]Coverage.py warning: Module thermoengine was never imported. (module-not-imported)


============================================== FAILURES ==============================================

______________________ Equilibrate/7a-K-Quartz-Equilibrate.ipynb::Cell 16 ______________________
Notebook cell execution failed
Cell 16: Cell execution caused an exception

Input:
state = equil.execute(t-20.0, p, state=state, debug=0, stats=True)
state.print_state()

Traceback:

---------------------------------------------------------------------------
ValueError                                Traceback (most recent call last)
<ipython-input-18-5341e00c0562> in <module>
----> 1 state = equil.execute(t-20.0, p, state=state, debug=0, stats=True)
      2 state.print_state()

~/opt/anaconda3/envs/thermoengine/lib/python3.8/site-packages/thermoengine/equilibrate.py in execute(self, t, p, bulk_comp, state, con_deltaNNO, debug, stats)
   3429                         else:
   3430                             self._print_matrix_repr(P_nz)
-> 3431                     A, df, Q1, Q2, R11 = self._compute_a_and_qr(t, p, state,
   3432                         P_nz, debug)
   3433                     if debug > 0:

~/opt/anaconda3/envs/thermoengine/lib/python3.8/site-packages/thermoengine/equilibrate.py in _compute_a_and_qr(self, t, p, state, P_nz, debug)
   2602                 print ('A matrix has reduced row space by:', row-row_rank)
   2603             df = col - row_rank
-> 2604         R, Q = sp.linalg.rq(A, mode='full')
   2605         if debug > 1:
   2606             print ('RQ = A, R matrix', R.shape)

~/opt/anaconda3/envs/thermoengine/lib/python3.8/site-packages/scipy/linalg/decomp_qr.py in rq(a, overwrite_a, lwork, mode, check_finite)
    389
    390     if check_finite:
--> 391         a1 = numpy.asarray_chkfinite(a)
    392     else:
    393         a1 = numpy.asarray(a)

~/opt/anaconda3/envs/thermoengine/lib/python3.8/site-packages/numpy/lib/function_base.py in asarray_chkfinite(a, dtype, order)
    486     a = asarray(a, dtype=dtype, order=order)
    487     if a.dtype.char in typecodes['AllFloat'] and not np.isfinite(a).all():
--> 488         raise ValueError(
    489             "array must not contain infs or NaNs")
    490     return a

ValueError: array must not contain infs or NaNs






====================================== short test summary info =======================================
FAILED Equilibrate/7a-K-Quartz-Equilibrate.ipynb::Cell 16
=================================== 1 failed, 17 passed in 18.41s ====================================
"""
