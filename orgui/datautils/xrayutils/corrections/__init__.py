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
r"""Everything that turns detector counts into a structure factor.

One home for the correction factors of

.. math::

    I = \Phi_0 \frac{r_e^2 A \lambda^2}{A_u^2}\,
        |F_{hkl}|^2 \, P \, \eta \, C_\mathrm{det}

so that a rocking scan, a stationary area-detector measurement, a
reflectivity curve and a reciprocal-space reconstruction of the same sample
are corrected by the same code and land on the same scale. The measurement
equation is E. Vlieg, *J. Appl. Cryst.* **30** (1997) 532, in the
two-dimensional-detector form of J. Drnec *et al.*,
*J. Appl. Cryst.* **47** (2014) 365; :mod:`~.measurement` documents it in
full.

The split between the modules is by *what a factor depends on*, which is also
what makes each of them testable on its own:

:mod:`~.geometry`
    Diffractometer angles only. The z-axis table of the ANA/ROD manual:
    Lorentz factors, rod interception, the geometric area factor.
:mod:`~.beamprofile`
    The vertical profile of the incident beam and its integrals over a
    finite sample.
:mod:`~.activearea`
    The illuminated active surface area :math:`A`, in square meter, in both
    the slit-limited and the beam-limited case.
:mod:`~.detector`
    Per-pixel factors of a detector image: solid angle and polarization.
:mod:`~.normalization`
    Counting time and monitor.
:mod:`~.roi`
    Reducing per-pixel factors onto a summed region of interest.
:mod:`~.measurement`
    Which of the above apply to which kind of scan, and the reduction to
    :math:`|F_{hkl}|^2` and to absolute reflectivity.

Everything here is physics: arrays and scalars in, arrays out. Nothing in
this package reads a scan object, a configuration file or a GUI widget.
Resolving *which* corrections an experiment asked for, and pulling counters
out of a beamline scan, belongs to :mod:`orgui.backend.scans` and to the
application layer.

For backwards compatibility ``orgui.datautils.xrayutils.geometrycorrections``
and ``orgui.datautils.xrayutils.beamprofile`` remain importable and re-export
:mod:`~.geometry` and :mod:`~.beamprofile` unchanged.
"""

from . import (  # noqa: F401
    activearea,
    beamprofile,
    detector,
    geometry,
    measurement,
    normalization,
    roi,
)

__all__ = [
    "activearea",
    "beamprofile",
    "detector",
    "geometry",
    "measurement",
    "normalization",
    "roi",
]
