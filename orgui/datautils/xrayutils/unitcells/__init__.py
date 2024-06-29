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

import os
import glob

_uc_path = os.path.dirname(__file__)

_files_bulk = glob.glob(os.path.join(_uc_path,"*.bul"))
_files_surf = glob.glob(os.path.join(_uc_path,"*.sur"))
_files_xtal = glob.glob(os.path.join(_uc_path,"*.xtal"))

availablebulk = [os.path.splitext(os.path.basename(bf))[0] for bf in _files_bulk]
availablesur = [os.path.splitext(os.path.basename(bf))[0] for bf in _files_surf]
availablextal = [os.path.splitext(os.path.basename(bf))[0] for bf in _files_xtal]

BULFILES = dict( (os.path.splitext(os.path.basename(bf))[0].lower(), bf) for bf in _files_bulk)
SURFILES = dict( (os.path.splitext(os.path.basename(bf))[0].lower(), bf) for bf in _files_surf)
XTALFILES = dict( (os.path.splitext(os.path.basename(bf))[0].lower(), bf) for bf in _files_xtal)


def unitcell(name: str):
    from .. import CTRcalc
    name = name.lower()
    if name in BULFILES:
        return CTRcalc.UnitCell.fromBULfile(BULFILES[name])
    elif name in SURFILES:
        return CTRcalc.UnitCell.fromSURfile(SURFILES[name])
    else:
        raise ValueError("Unit cell %s is not available." % name)


def crystal(name: str):
    from .. import CTRcalc
    name = name.lower()
    if name in XTALFILES:
        return CTRcalc.SXRDCrystal.fromFile(XTALFILES[name])
    else:
        raise ValueError("Crystal %s is not available." % name)