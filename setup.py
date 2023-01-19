#!/usr/bin/env python
""" file: setup.py
    modified: Mark S. Ghiorso, OFM Research
    date: June 12, 2017, rev June 27, 2017, rev cython Dec 19, 2019

    description: Distutils installer script for thermoengine.
"""
from setuptools import setup, find_packages
from setuptools.extension import Extension
from Cython.Build import cythonize
import numpy

from sys import platform
if platform == "linux" or platform == "linux2":
    from distutils import sysconfig
elif platform == "darwin":
    pass
elif platform == "win32":
    pass

extensions = [
    Extension(
        "thermoengine.aqueous",
        sources=["thermoengine/aqueous/aqueous.pyx",
        "thermoengine/aqueous/swim.c",
        "thermoengine/aqueous/born.c",
        "thermoengine/aqueous/duanzhang.c",
        "thermoengine/aqueous/holten.c",
        "thermoengine/aqueous/wagner.c",
        "thermoengine/aqueous/zhangduan.c",
        "thermoengine/aqueous/FreeSteam2.1/b23.c",
        "thermoengine/aqueous/FreeSteam2.1/backwards.c",
        "thermoengine/aqueous/FreeSteam2.1/bounds.c",
        "thermoengine/aqueous/FreeSteam2.1/common.c",
        "thermoengine/aqueous/FreeSteam2.1/derivs.c",
        "thermoengine/aqueous/FreeSteam2.1/region1.c",
        "thermoengine/aqueous/FreeSteam2.1/region2.c",
        "thermoengine/aqueous/FreeSteam2.1/region3.c",
        "thermoengine/aqueous/FreeSteam2.1/region4.c",
        "thermoengine/aqueous/FreeSteam2.1/solver2.c",
        "thermoengine/aqueous/FreeSteam2.1/steam.c",
        "thermoengine/aqueous/FreeSteam2.1/steam_ph.c",
        "thermoengine/aqueous/FreeSteam2.1/steam_ps.c",
        "thermoengine/aqueous/FreeSteam2.1/steam_pT.c",
        "thermoengine/aqueous/FreeSteam2.1/steam_pu.c",
        "thermoengine/aqueous/FreeSteam2.1/steam_pv.c",
        "thermoengine/aqueous/FreeSteam2.1/steam_Ts.c",
        "thermoengine/aqueous/FreeSteam2.1/steam_Tx.c",
        "thermoengine/aqueous/FreeSteam2.1/surftens.c",
        "thermoengine/aqueous/FreeSteam2.1/thcond.c",
        "thermoengine/aqueous/FreeSteam2.1/viscosity.c",
        "thermoengine/aqueous/FreeSteam2.1/zeroin.c"],
        include_dirs=['./thermoengine/aqueous', './thermoengine/aqueous/FreeSteam2.1', numpy.get_include()],
        extra_compile_args=['-O3', '-DNPY_NO_DEPRECATED_API=NPY_1_7_API_VERSION'],
        libraries=['gsl'],
        library_dirs=['/usr/local/lib'],
        runtime_library_dirs=['/usr/local/lib']
    ),
]

if platform == "linux" or platform == "linux2":
    sysconfig.get_config_vars()['CC'] = 'clang'
    sysconfig.get_config_vars()['CXX'] = 'clang++'
    sysconfig.get_config_vars()['CCSHARED'] = '-fPIC'
    sysconfig.get_config_vars()['LDSHARED'] = 'clang -shared'

def readme():
    with open('README.rst') as f:
        return f.read()

setup(
      name='thermoengine',
      version='1.0',
      description='Principal Python package for ENKI thermodynamics modules',
      long_description=readme(),
      url='http://gitlab.com/enki-portal/ThermoEngine',
      author='Aaron S. Wolf; Mark S. Ghiorso',
      author_email='aswolf@umich.edu, ghiorso@ofm-research.org',
      license='GNU AFFERO GENERAL PUBLIC LICENSE Version 3',
      packages=[
          'thermoengine',
      ],
      ext_modules = cythonize(extensions),
      include_package_data=True,
      install_requires=[
         'coverage',
         'cython',
         'deprecation',
         'ipykernel',
         'jupyter',
         'jupyterlab',
         'matplotlib',
         'nbval',
         'numdifftools',
         'numpy>=1.21.2',
         'openpyxl',
         'pandas',
         'pytest',
         'pytest-cov',
         'scipy',
         'seaborn',
         'statsmodels',
         'sympy',
         'cmake>=3.14.0',
         'SulfLiq',
         'elasticsearch',
         'elasticsearch-dsl',
      ] + ['rubicon-objc>=0.4.2' if platform == 'darwin' else 'rubicon-objc==0.2.10'],
      classifiers=[
        'Development Status :: 3 - Beta',
        'Environment :: Plugins',
        'Intended Audience :: Developers',
        'Intended Audience :: Education',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: BSD License',
        'Natural Language :: English',
        'Operating System :: OS Independent',
        'Programming Language :: Python',
        'Topic :: Scientific/Engineering',
        'Topic :: Scientific/Engineering :: Chemistry',
        'Topic :: Scientific/Engineering :: Mathematics',
        'Topic :: Scientific/Engineering :: Physics',
        'Topic :: Software Development :: Libraries :: Python Modules'
      ],
      zip_safe=False,
      test_suite='nose.collector',
      tests_require=['nose']
)
