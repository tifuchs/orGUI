# -*- coding: utf-8 -*-
# /*##########################################################################
#
# Copyright (c) 2020-2024 Timo Fuchs
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.
#
# ###########################################################################*/
"""Module descripiton

"""
__author__ = "Timo Fuchs"
__copyright__ = "Copyright 2020-2024, Timo Fuchs"
__credits__ = []
__license__ = "MIT License"
__version__ = "1.0.0"
__maintainer__ = "Timo Fuchs"
__email__ = "fuchs@physik.uni-kiel.de"

import os
import logging

logger = logging.getLogger('orgui')

#from distutils.core import setup

packages = ['orgui',"orgui.app", "orgui.backend", "orgui.backend.beamline", 
            "orgui.resources", "orgui.datautils", "orgui.datautils.xrayutils",
            "orgui.datautils.xrayutils.unitcells"]
package_data = {'orgui.resources' : ['icons/*.svg','icons/*.png'],
                "orgui.datautils.xrayutils.unitcells" : ['*.bul', '*.sur', '*.xtal', '*.xpr']}
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
      package_data=package_data,
      author="Timo Fuchs",
      author_email="fuchs@physik.uni-kiel.de",
      url='FIXME',
      license='All rights reserved',
      install_requires=[
      'numpy >= 1.12',
      'scipy >= 1.0',
      'matplotlib >= 1.3',
      'fabio >= 0.7',
      'silx >= 1.1',
      'pyFAI >= 0.16',
      'pytz >= 2022'],
      classifiers=[
          'Topic :: Scientific/Engineering',
          'Development Status :: 3 - Alpha',
          'Operating System :: OS Independent',
          'Programming Language :: Python ']
      )
