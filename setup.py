'''
Description: setup.py for allskyScene package
Author: Hejun Xie
Date: 2022-05-18 17:08:52
LastEditors: Hejun Xie
LastEditTime: 2022-05-20 01:03:31
'''


from distutils.core import setup
from Cython.Build import cythonize
import numpy

import Cython.Compiler.Options
Cython.Compiler.Options.annotate = True


DISTNAME = "ZJU_AERO"
AUTHOR = 'Hejun Xie - Zhejiang University'
AUTHOR_EMAIL = 'hejun.xie@zju.edu.cn'
URL = "https://github.com/Usami-Renko/allskyScene"
LICENCE = 'MIT'
PYTHON_REQUIRES = ">=3.7"
INSTALL_REQUIRES = ['numpy', 'scipy', 'xarray', 'h5py', 'cartopy', 'proplot', 'matplotlib', 'satpy']
DESCRIPTION = "Allsky quick overview"
LONG_DESCRIPTION = """For allsky quick overview.
"""
PLATFORMS = ["Linux"]
CLASSIFIERS = [
    'Development Status :: 1 - Planning',
    'Intended Audience :: Science/Research',
    'Intended Audience :: Developers',
    'License :: OSI Approved :: MIT License',
    'Programming Language :: Python',
    'Programming Language :: Python :: 3',
    'Programming Language :: Python :: 3.7',
    'Topic :: Scientific/Engineering',
    'Topic :: Scientific/Engineering :: Atmospheric Science',
    'Operating System :: POSIX :: Linux']
PACKAGES = ['allskyScene']

setup(
  name = 'allskyScene',
  version = '0.1.0',
  author       = AUTHOR,
  license      = LICENCE,
  author_email = AUTHOR_EMAIL,
  description  = DESCRIPTION,
  long_description  = LONG_DESCRIPTION,
  python_requires   = PYTHON_REQUIRES,
  install_requires  = INSTALL_REQUIRES,
  url               = URL,
  platforms         = PLATFORMS,
  classifiers       = CLASSIFIERS,
  packages          = PACKAGES,
  ext_modules = cythonize('./allskyScene/monotonic.pyx', compiler_directives={'language_level': "3"}),
  include_dirs=[numpy.get_include()]
)
