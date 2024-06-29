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
__author__ = "Timo Fuchs"
__copyright__ = "Copyright 2020-2024 Timo Fuchs"
__license__ = "MIT License"
__version__ = "1.0.0"
__maintainer__ = "Timo Fuchs"
__email__ = "fuchs@physik.uni-kiel.de"

import unittest

from .. import CTRcalc, CTRplotutil
from ... import util

import numpy as np
import os


class TestReadSXRDCrystal(unittest.TestCase):

    xpr_file = """E = 68.00000 keV
# UnitCell relaxations
0000 occupancy = 1.00000 +- nan
return
Coherent 1.00000   1.00000 2.00000 3.00000 4.00000 5.00000 6.00000 7.00000 8.00000 9.00000 10.00000 11.00000 0.00000
2.7748 2.7748 3.9242 90.0000 90.0000 90.0000
Name   x/frac     y/frac     z/frac     iDW     oDW      occup
00  Pt  (0.00000 +- nan)  (0.00000 +- nan)  (0.00249 +- 0.00007)  (0.4871 +- nan)  (0.5663 +- nan)  (1.0000 +- nan)
01  Pt  (0.50000 +- nan)  (0.50000 +- nan)  (0.51591 +- 0.00010)  (0.7070 +- nan)  (0.8379 +- nan)  (1.0000 +- nan)

# UnitCell adsorbates
0001 occupancy = 1.00000 +- nan
return
Coherent 0.50000   1.00000 0.00000 0.00000 0.00000 1.00000 0.00000 0.00000 0.00000 1.00000 0.00000 0.00000 0.00000
Coherent 0.50000   0.00000 -1.00000 0.00000 1.00000 0.00000 0.00000 0.00000 0.00000 1.00000 0.00000 0.00000 0.00000
2.7748 2.7748 3.9242 90.0000 90.0000 90.0000
Name   x/frac     y/frac     z/frac     iDW     oDW      occup
00  O  (0.00000 +- nan)  (0.50000 +- nan)  (1.00100 +- nan)  (0.4350 +- nan)  (0.4350 +- nan)  (0.0000 +- nan)

# UnitCell bulk
return
Coherent 1.00000   1.00000 0.00000 0.00000 0.00000 1.00000 0.00000 0.00000 0.00000 1.00000 0.00000 0.00000 0.00000
2.7748 2.7748 3.9242 90.0000 90.0000 90.0000
Name   x/frac     y/frac     z/frac     iDW     oDW      occup
00  Pt  (0.00000 +- nan)  (0.00000 +- nan)  (-1.00000 +- nan)  (0.4350 +- nan)  (0.4350 +- nan)  (1.0000 +- nan)
01  Pt  (0.50000 +- nan)  (0.50000 +- nan)  (-0.50000 +- nan)  (0.4350 +- nan)  (0.4350 +- nan)  (1.0000 +- nan)
"""

    xpr_file_orig = """E = 68.00000 keV
# UnitCell relaxations
0000 occupancy = 1.00000 +- nan
return
Coherent 1.00000   1.00000 0.00000 0.00000 0.00000 1.00000 0.00000 0.00000 0.00000 1.00000 0.00000 0.00000 0.00000
2.7748 2.7748 3.9242 90.0000 90.0000 90.0000
Name   x/frac     y/frac     z/frac     iDW     oDW      occup
00  Pt  (0.00000 +- nan)  (0.00000 +- nan)  (0.00249 +- 0.00007)  (0.4871 +- nan)  (0.5663 +- nan)  (1.0000 +- nan)
01  Pt  (0.50000 +- nan)  (0.50000 +- nan)  (0.51591 +- 0.00010)  (0.7070 +- nan)  (0.8379 +- nan)  (1.0000 +- nan)

# UnitCell adsorbates
0001 occupancy = 1.00000 +- nan
return
Coherent 0.50000   1.00000 0.00000 0.00000 0.00000 1.00000 0.00000 0.00000 0.00000 1.00000 0.00000 0.00000 0.00000
Coherent 0.50000   0.00000 -1.00000 0.00000 1.00000 0.00000 0.00000 0.00000 0.00000 1.00000 0.00000 0.00000 0.00000
2.7748 2.7748 3.9242 90.0000 90.0000 90.0000
Name   x/frac     y/frac     z/frac     iDW     oDW      occup
00  O  (0.00000 +- nan)  (0.50000 +- nan)  (1.00100 +- nan)  (0.4350 +- nan)  (0.4350 +- nan)  (0.0000 +- nan)

# UnitCell bulk
return
Coherent 1.00000   1.00000 0.00000 0.00000 0.00000 1.00000 0.00000 0.00000 0.00000 1.00000 0.00000 0.00000 0.00000
2.7748 2.7748 3.9242 90.0000 90.0000 90.0000
Name   x/frac     y/frac     z/frac     iDW     oDW      occup
00  Pt  (0.00000 +- nan)  (0.00000 +- nan)  (-1.00000 +- nan)  (0.4350 +- nan)  (0.4350 +- nan)  (1.0000 +- nan)
01  Pt  (0.50000 +- nan)  (0.50000 +- nan)  (-0.50000 +- nan)  (0.4350 +- nan)  (0.4350 +- nan)  (1.0000 +- nan)
"""

    def testFromStr(self):
        xtal = CTRcalc.SXRDCrystal.fromStr(TestReadSXRDCrystal.xpr_file)
        self.assertIsInstance(xtal['relaxations'], CTRcalc.UnitCell)
        with self.assertRaises(ValueError):
            uc = xtal['notexisting']
        self.assertEqual(str(xtal), TestReadSXRDCrystal.xpr_file)
        
        self.assertTrue(np.array_equal(xtal['relaxations'].coherentDomainMatrix[0], 
            np.array([[1.00000, 2.00000, 3.00000, 10.00000], 
                      [4.00000, 5.00000, 6.00000, 11.00000],
                      [7.00000, 8.00000, 9.00000, 0.00000]])))
    
    def testFromFile(self):
        fp = os.path.split(__file__)[0]
        xtal = CTRcalc.SXRDCrystal.fromFile(os.path.join(fp,"testdata/0V12_calculated.xpr"))
        self.assertEqual(str(xtal), TestReadSXRDCrystal.xpr_file_orig)
        
    
    
class TestCTRcalculationNumPy(unittest.TestCase):
    def setUp(self):
        fp = os.path.split(__file__)[0]
        self.xtal_unitcells = CTRcalc.SXRDCrystal.fromFile(os.path.join(fp,"testdata/0V12_calculated.xpr"))
        self.CTRs = CTRplotutil.CTRCollection.fromANAROD(os.path.join(fp,"testdata/0V12_calculated.dat"), RODexport=True)
        pt100 = CTRcalc.UnitCell([3.9242, 3.9242, 3.9242] ,[90.0000 ,90.0000, 90.0000])
        self.xtal_unitcells.setGlobalReferenceUnitCell(pt100,util.z_rotation(np.deg2rad(45.)))
        CTRcalc.HAS_NUMBA_ACCEL = False
        
    def testStructureFactorEqual(self):
        calc_CTRs = self.CTRs.generateCollectionFromXtal(self.xtal_unitcells)

        for calc, reference in zip(calc_CTRs, self.CTRs):
            self.assertTrue(np.allclose(calc.sfI, reference.sfI, rtol=1e-02))
        
        
class TestCTRcalculationNumba(unittest.TestCase):
    def setUp(self):
        fp = os.path.split(__file__)[0]
        self.xtal_unitcells = CTRcalc.SXRDCrystal.fromFile(os.path.join(fp,"testdata/0V12_calculated.xpr"))
        self.CTRs = CTRplotutil.CTRCollection.fromANAROD(os.path.join(fp,"testdata/0V12_calculated.dat"), RODexport=True)
        pt100 = CTRcalc.UnitCell([3.9242, 3.9242, 3.9242] ,[90.0000 ,90.0000, 90.0000])
        self.xtal_unitcells.setGlobalReferenceUnitCell(pt100,util.z_rotation(np.deg2rad(45.)))
        if hasattr(CTRcalc, "_CTRcalc_accel"):
            CTRcalc.HAS_NUMBA_ACCEL = True
        else:
            raise Exception("Can not perform Numba tests: _CTRcalc_accel library not imported. Is Numba installed?")
            
    def testStructureFactorEqual(self):
        calc_CTRs = self.CTRs.generateCollectionFromXtal(self.xtal_unitcells)

        for calc, reference in zip(calc_CTRs, self.CTRs):
            self.assertTrue(np.allclose(calc.sfI, reference.sfI, rtol=1e-02))
        

