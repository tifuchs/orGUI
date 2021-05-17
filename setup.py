# -*- coding: utf-8 -*-
###############################################################################
# Copyright (c) 2020 Timo Fuchs, Olaf Magnussen all rights reserved
#
# This software was developed during the PhD work of Timo Fuchs,
# within the group of Olaf Magnussen. Usage within the group is hereby granted.
###############################################################################
"""Module descripiton

"""
__author__ = "Timo Fuchs"
__copyright__ = "Copyright 2020, Timo Fuchs, Olaf Magnussen all rights reserved"
__credits__ = []
__license__ = "all rights reserved"
__version__ = "1.0.0"
__maintainer__ = "Timo Fuchs"
__email__ = "fuchs@physik.uni-kiel.de"

import os
import logging

logger = logging.getLogger('orgui')

#from distutils.core import setup

packages = ['orgui',"orgui.app"]
#package_data = {'datautils.xrayutils' : ['element_colors.txt','water_scattering.dat','water_density_molecularformfactor.dat']}

dry_run = False

if dry_run:
	# DRY_RUN implies actions which do not require dependencies, like NumPy
	try:
		from setuptools import setup
		logger.info("Use setuptools.setup")
	except ImportError:
		from distutils.core import setup
		logger.info("Use distutils.core.setup")
else:
	try:
		from setuptools import setup
	except ImportError:
		from numpy.distutils.core import setup
		logger.info("Use numpy.distutils.setup")


setup(name='orgui', version='1.0.0',
      description='User interface for angle calculations for HESXRD',
      long_description='FIXME',
      packages=packages,
      entry_points = {
        'console_scripts': ['orGUI=orgui.main:main']
      },
      author="Timo Fuchs",
      author_email="fuchs@physik.uni-kiel.de",
      url='FIXME',
      license='All rights reserved',
      install_requires=[
      'numpy >= 1.12',
      'scipy >= 1.0',
      'matplotlib >= 1.3',
      'fabio >= 0.7',
      'silx >= 0.9',
      'pyFAI >= 0.16'],
      classifiers=[
          'Topic :: Scientific/Engineering',
          'Development Status :: 3 - Alpha',
          'Operating System :: OS Independent',
          'Programming Language :: Python ']
      )
