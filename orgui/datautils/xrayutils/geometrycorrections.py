# /*##########################################################################
#
# Copyright (c) 2026 Timo Fuchs
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
"""Backwards-compatible alias of :mod:`orgui.datautils.xrayutils.corrections.geometry`.

The z-axis geometrical corrections moved into the
:mod:`orgui.datautils.xrayutils.corrections` package, which collects every
factor between detector counts and a structure factor. This module was
released under its own name, so it stays importable and re-exports the same
objects. New code should import from
:mod:`orgui.datautils.xrayutils.corrections.geometry`.
"""

from .corrections.geometry import (  # noqa: F401
    AREA,
    REFLECTIVITY_ROCKING,
    ROCKING,
    STATIONARY,
    area_correction,
    lorentz_factor,
    lorentz_reflectivity_rocking_scan,
    lorentz_rocking_scan,
    lorentz_stationary,
    rod_interception,
)

__all__ = [
    "AREA",
    "REFLECTIVITY_ROCKING",
    "ROCKING",
    "STATIONARY",
    "area_correction",
    "lorentz_factor",
    "lorentz_reflectivity_rocking_scan",
    "lorentz_rocking_scan",
    "lorentz_stationary",
    "rod_interception",
]
