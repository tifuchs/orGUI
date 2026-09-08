# /*##########################################################################
#
# Copyright (c) 2020-2025 Timo Fuchs
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
__copyright__ = "Copyright 2020-2025 Timo Fuchs"
__license__ = "MIT License"
__version__ = "1.3.0"
__maintainer__ = "Timo Fuchs"
__email__ = "tfuchs@cornell.edu"


import os
import matplotlib.pyplot as plt
import matplotlib.figure as mplfig
import matplotlib.transforms as mtransforms
import scipy.interpolate as interp
import numpy as np
import math
import json
import warnings
from collections import OrderedDict
from dataclasses import dataclass
from typing import Literal
from .CTRcalc import SXRDCrystal
from .. import util


@dataclass(frozen=True, eq=False)
class PolarizationReduction:
    """Describe polarization already represented by a CTR measurement.

    :param float s_fraction:
        Incoherent incident s-polarization fraction in the local Renaud basis.
        One is pure s and zero is pure p.
    :param str outgoing:
        Analysed outgoing channel, ``"s"`` or ``"p"``, or ``"unanalysed"``.
    :param polarization_factor:
        Optional conventional pointwise intensity factor P. The array is
        copied, converted to float, and stored read-only.
    :type polarization_factor: array-like or None
    :raises ValueError:
        If the mixture, outgoing channel, or supplied factor is invalid.

    Instances compare by value but are intentionally unhashable because the
    pointwise factor is an array.
    """

    s_fraction: float
    outgoing: Literal["s", "p", "unanalysed"]
    polarization_factor: np.ndarray | None = None
    __hash__ = None

    def __post_init__(self):
        if isinstance(self.s_fraction, bool | np.bool_):
            raise ValueError("s_fraction must be a finite scalar in [0, 1].")
        fraction = np.asarray(self.s_fraction)
        if fraction.ndim != 0:
            raise ValueError("s_fraction must be a finite scalar in [0, 1].")
        fraction = float(fraction)
        if not np.isfinite(fraction) or not 0.0 <= fraction <= 1.0:
            raise ValueError("s_fraction must be a finite scalar in [0, 1].")
        if self.outgoing not in {"s", "p", "unanalysed"}:
            raise ValueError('outgoing must be "s", "p", or "unanalysed".')
        object.__setattr__(self, "s_fraction", fraction)

        if self.polarization_factor is not None:
            factor = np.array(
                self.polarization_factor, dtype=np.float64, copy=True
            )
            if not np.all(np.isfinite(factor)) or np.any(factor <= 0.0):
                raise ValueError(
                    "polarization_factor must contain finite, strictly positive P."
                )
            factor.setflags(write=False)
            object.__setattr__(self, "polarization_factor", factor)

    def __eq__(self, other):
        if not isinstance(other, PolarizationReduction):
            return NotImplemented
        if (
            self.s_fraction != other.s_fraction
            or self.outgoing != other.outgoing
        ):
            return False
        if self.polarization_factor is None or other.polarization_factor is None:
            return (
                self.polarization_factor is None
                and other.polarization_factor is None
            )
        return np.array_equal(
            self.polarization_factor, other.polarization_factor
        )


@dataclass(frozen=True, eq=False)
class MeasurementReduction:
    """Describe the stored quantity and its polarization reduction.

    :param str quantity:
        ``"structure_factor"`` for stored F or ``"reflectivity"`` for
        corrected field-intensity ratio R.
    :param PolarizationReduction polarization:
        Polarization provenance, or ``None`` when it is unknown.

    Instances compare by value but are intentionally unhashable because their
    polarization metadata may own an array.
    """

    quantity: Literal["structure_factor", "reflectivity"] = "structure_factor"
    polarization: PolarizationReduction | None = None
    __hash__ = None

    def __post_init__(self):
        if self.quantity not in {"structure_factor", "reflectivity"}:
            raise ValueError(
                'quantity must be "structure_factor" or "reflectivity".'
            )
        if self.polarization is not None and not isinstance(
            self.polarization, PolarizationReduction
        ):
            raise TypeError(
                "polarization must be a PolarizationReduction or None."
            )

    def __eq__(self, other):
        if not isinstance(other, MeasurementReduction):
            return NotImplemented
        return (
            self.quantity == other.quantity
            and self.polarization == other.polarization
        )


@dataclass(frozen=True)
class CTRScanGeometry:
    """Record the Vlieg z-mode scan rule for one measured CTR.

    :param str fixed:
        ``"in"`` for fixed incidence, ``"out"`` for fixed exit, or ``"eq"``
        for equal incident and exit angles.
    :param float angle:
        Fixed angle in rad for ``"in"`` and ``"out"``. Equal-angle scans
        require ``None``.
    :param bool mirrorx:
        Select the negative-delta scattering branch.
    """

    fixed: Literal["in", "out", "eq"]
    angle: float | None = None
    mirrorx: bool = False

    def __post_init__(self):
        if self.fixed not in {"in", "out", "eq"}:
            raise ValueError('fixed must be "in", "out", or "eq".')
        if not isinstance(self.mirrorx, bool | np.bool_):
            raise TypeError("mirrorx must be boolean.")
        object.__setattr__(self, "mirrorx", bool(self.mirrorx))
        if self.fixed == "eq":
            if self.angle is not None:
                raise ValueError('angle must be None when fixed="eq".')
            return
        if self.angle is None:
            raise ValueError(
                'angle must be a finite scalar in (0, pi/2] rad for "in"/"out".'
            )
        if isinstance(self.angle, bool | np.bool_):
            raise ValueError(
                'angle must be a finite scalar in (0, pi/2] rad for "in"/"out".'
            )
        angle = np.asarray(self.angle)
        if angle.ndim != 0:
            raise ValueError(
                'angle must be a finite scalar in (0, pi/2] rad for "in"/"out".'
            )
        angle = float(angle)
        if not np.isfinite(angle) or not 0.0 < angle <= np.pi / 2.0:
            raise ValueError(
                'angle must be a finite scalar in (0, pi/2] rad for "in"/"out".'
            )
        object.__setattr__(self, "angle", angle)


_DEFAULT_REDUCTION = object()
_ANGLE_DTYPE = np.dtype(
    [
        ("alpha", "f8"),
        ("delta", "f8"),
        ("gamma", "f8"),
        ("omega", "f8"),
        ("chi", "f8"),
        ("phi", "f8"),
    ]
)
_QUANTITY_YLABELS = {
    "structure_factor": "Structure factor / arb. units",
    "reflectivity": "Reflectivity / dimensionless",
}


def _signed_sqrt(values):
    """Map signed intensity values reversibly to signed amplitudes."""
    values = np.asarray(values, dtype=np.float64)
    return np.sign(values) * np.sqrt(np.abs(values))


def _signed_sqrt_uncertainty(intensity, uncertainty):
    """Return the stable transformed interval half-width."""
    intensity = np.asarray(intensity, dtype=np.float64)
    uncertainty = np.asarray(uncertainty, dtype=np.float64)
    if uncertainty.shape != intensity.shape:
        raise ValueError("Intensity uncertainties must match the CTR data shape.")

    finite_intensity = np.isfinite(intensity)
    invalid_uncertainty = finite_intensity & (
        ~np.isfinite(uncertainty) | (uncertainty <= 0.0)
    )
    if np.any(invalid_uncertainty):
        raise ValueError(
            "Finite intensities require finite, strictly positive uncertainties."
        )

    transformed = np.full_like(intensity, np.nan, dtype=np.float64)
    magnitude = np.abs(intensity[finite_intensity])
    sigma = uncertainty[finite_intensity]
    strong = magnitude > sigma
    result = np.empty_like(magnitude)
    strong_ratio = sigma[strong] / magnitude[strong]
    result[strong] = (sigma[strong] / np.sqrt(magnitude[strong])) / (
        np.sqrt(1.0 + strong_ratio) + np.sqrt(1.0 - strong_ratio)
    )
    weak_ratio = magnitude[~strong] / sigma[~strong]
    result[~strong] = 0.5 * np.sqrt(sigma[~strong]) * (
        np.sqrt(1.0 + weak_ratio) + np.sqrt(1.0 - weak_ratio)
    )
    transformed[finite_intensity] = result
    return transformed


def _nx_text(value):
    """Return one NeXus scalar attribute as text."""
    value = np.asarray(value)
    if value.size != 1:
        raise ValueError("Expected one scalar NeXus metadata value.")
    value = value.reshape(-1)[0]
    if isinstance(value, bytes):
        return value.decode("utf-8")
    return str(value)


def _nx_bool(value):
    """Return one NeXus scalar attribute as bool."""
    value = _nx_text(value).lower()
    if value in {"true", "1"}:
        return True
    if value in {"false", "0"}:
        return False
    raise ValueError("NeXus mirrorx metadata must be boolean.")


def _select_reduction(reduction, selection):
    """Return reduction metadata after applying a point selection."""
    polarization = reduction.polarization
    if polarization is None or polarization.polarization_factor is None:
        return reduction
    selected = PolarizationReduction(
        polarization.s_fraction,
        polarization.outgoing,
        polarization.polarization_factor[selection],
    )
    return MeasurementReduction(reduction.quantity, selected)


def _reduction_for_import(reduction, selection):
    """Select a full-table polarization factor for one imported rod."""
    if reduction is _DEFAULT_REDUCTION:
        return reduction
    polarization = reduction.polarization
    if polarization is None or polarization.polarization_factor is None:
        return reduction
    factor = polarization.polarization_factor
    if factor.shape == np.asarray(selection).shape:
        return _select_reduction(reduction, selection)
    return reduction


def _load_ctr_file(filename, reduction):
    """Load an orGUI CTR table and its metadata header."""
    metadata = {}
    with open(filename, encoding="utf-8") as stream:
        for raw_line in stream:
            line = raw_line.strip()
            if not line:
                continue
            if not line.startswith("#"):
                break
            comment = line[1:].strip()
            if ":" in comment:
                key, value = comment.split(":", 1)
                metadata[key.strip().lower()] = value.strip()

    if metadata.get("orgui_ctr_schema") != "2":
        raise ValueError(
            "orGUI CTR files require '# orgui_ctr_schema: 2'; use "
            "CTRCollection.fromANAROD for legacy files."
        )
    if "columns" not in metadata:
        raise ValueError("orGUI CTR files require a columns header.")

    table = np.atleast_2d(np.loadtxt(filename))
    column_names = metadata["columns"].lower().split()
    if len(column_names) != table.shape[1]:
        raise ValueError(
            "orGUI CTR columns metadata must match the numerical table."
        )
    if (
        len(column_names) not in {5, 6}
        or column_names[:3] != ["h", "k", "l"]
        or (
            len(column_names) == 6
            and column_names[5] != "polarization_factor"
        )
    ):
        raise ValueError(
            "orGUI CTR files require H, K, L, value, and uncertainty columns "
            "with an optional trailing polarization_factor."
        )
    if reduction is not _DEFAULT_REDUCTION:
        return table, reduction

    quantity = metadata.get("quantity")
    s_fraction = metadata.get("s_fraction")
    outgoing = metadata.get("outgoing")
    if quantity is None:
        raise ValueError("orGUI CTR files require quantity metadata.")
    if (s_fraction is None) != (outgoing is None):
        raise ValueError(
            "orGUI CTR metadata must provide both s_fraction and outgoing."
        )

    polarization = None
    if s_fraction is not None:
        factor = table[:, 5] if len(column_names) == 6 else None
        polarization = PolarizationReduction(
            float(s_fraction), outgoing, factor
        )
    return table, MeasurementReduction(
        quantity, polarization
    )


def _calculate_angles_zmode(
    h,
    k,
    l_values,
    vliegangles,
    fixedangle,
    fixed="in",
    chi=0.0,
    phi=0.0,
    *,
    hkl_transform=None,
    **keyargs,
):
    """Return round-trip-validated Vlieg z-mode angle records in rad."""
    h, k, l_values = np.broadcast_arrays(
        np.asarray(h, dtype=np.float64),
        np.asarray(k, dtype=np.float64),
        np.asarray(l_values, dtype=np.float64),
    )
    reference_hkl = np.vstack((h.ravel(), k.ravel(), l_values.ravel()))
    calculator_hkl = reference_hkl
    if hkl_transform is not None:
        transform = np.asarray(hkl_transform, dtype=np.float64)
        if transform.shape != (3, 3) or not np.all(np.isfinite(transform)):
            raise ValueError("hkl_transform must be a finite 3-by-3 matrix.")
        if np.linalg.matrix_rank(transform) != 3:
            raise ValueError("hkl_transform must be nonsingular.")
        calculator_hkl = transform @ reference_hkl

    positions = np.asarray(
        vliegangles.anglesZmode(
            calculator_hkl,
            fixedangle,
            fixed=fixed,
            chi=chi,
            phi=phi,
            **keyargs,
        ),
        dtype=np.float64,
    )
    if positions.ndim == 1:
        positions = positions.reshape(1, -1)
    if positions.shape != (reference_hkl.shape[1], 6):
        raise ValueError("Vlieg angle solver returned an unexpected array shape.")
    records = np.core.records.fromarrays(positions.T, dtype=_ANGLE_DTYPE)
    if hasattr(vliegangles, "anglesToHkl"):
        reconstructed = np.vstack(
            vliegangles.anglesToHkl(
                *(records[name] for name in records.dtype.names)
            )
        )
        if not np.allclose(
            reconstructed, calculator_hkl, rtol=1e-7, atol=1e-8
        ):
            raise ValueError(
                "Calculated Vlieg angles do not reproduce every requested H, K, "
                "and L."
            )
    return records


# don't init CTRFigure directly, use ctrfigure instead
class CTRFigure(mplfig.Figure):
    def __init__(self, **figargs):
        super().__init__(**figargs)
        self.data = OrderedDict()
        self.xlabels = "L / r.l.u."
        self.ylabels = "Structure factor / arb. units"
        self._ylabels_explicit = False
        self.xlim = None  # [0,5]
        self.ylim = None  # [3e1,2e3]
        self.wspace = 0
        self.hspace = 0.05

    def settings(self, **settings):
        if "ylabels" in settings:
            self._ylabels_explicit = True
        self.__dict__.update(settings)

    def addCTR(self, ctr, *plotargs, **keyargs):
        if ctr.ctr_id in self.data:
            existing_quantity = self.data[ctr.ctr_id][0][0].reduction.quantity
            if existing_quantity != ctr.reduction.quantity:
                raise ValueError(
                    "Datasets stored as different quantities cannot share one "
                    "CTR plot axis without explicit conversion."
                )
            self.data[ctr.ctr_id].append([ctr, plotargs, keyargs])
        else:
            self.data[ctr.ctr_id] = [[ctr, plotargs, keyargs]]

    def addCollection(self, collection):
        for rod in collection:
            self.addCTR(rod, *collection.plotsett, **collection.plotkeyargs)

    def generateCTRplot(self, cols=2, maxrelerr=None, **keyargs):
        rows = math.ceil(float(len(self.data)) / float(cols))
        self.add_subplot(rows, cols, 1)
        [self.add_subplot(rows, cols, i) for i in range(2, rows * cols + 1)]
        self.axes_ctr_id = []
        self.axes_ctr_hk = []
        share = keyargs.get("share_axes", True)
        self.to_sharex = np.array([share for i in range(len(self.axes))])
        self.to_sharey = np.array([share for i in range(len(self.axes))])
        if keyargs.get("sort_hk", False):
            data = sorted(self.data, key=lambda x: np.sum(np.abs(x[0])))
        else:
            data = self.data
        offset = 0
        panel_quantities = {}
        for i, ctrkey in enumerate(data):
            i += offset
            while i in keyargs.get("skip_panel", []):
                i += 1
                offset += 1
                self.axes_ctr_id.append(i)
                self.axes_ctr_hk.append(i)
            ax = self.axes[i]
            ctrname, ctrid = ctrkey
            panel_quantities[i] = self.data[ctrkey][0][0].reduction.quantity
            self.axes_ctr_id.append(ctrkey)
            self.axes_ctr_hk.append(ctrname)

            def formatMillerName(h):
                h = round(h, 3)
                if not h % 1:
                    h = int(h)
                if h < 0:
                    return rf"\overline{{{abs(h)}}}"
                else:
                    return rf"{abs(h)}"

            if "rodLabelL" in keyargs:
                ctrstr = rf"$( \, {formatMillerName(ctrname[0])} \, {formatMillerName(ctrname[1])} \, \ell \, )$"  # noqa: E501
            else:
                ctrstr = rf"$( \,{formatMillerName(ctrname[0])} \, {formatMillerName(ctrname[1])} \,)$"  # noqa: E501

            if hasattr(self, "rodlabelsize"):
                size = self.rodlabelsize
            else:
                size = None

            if hasattr(self, "rodlabelweight"):
                weight = self.rodlabelweight
            else:
                weight = "bold"

            if "rodlabel" in keyargs:
                if keyargs["rodlabel"] == "bottom":
                    ax.text(
                        keyargs.get("labelxpos", 0.02),
                        keyargs.get("labelypos", 0.02),
                        ctrstr,
                        horizontalalignment="left",
                        verticalalignment="bottom",
                        size=size,
                        transform=ax.transAxes,
                        weight=weight,
                    )
                if keyargs["rodlabel"] == "topright":
                    ax.text(
                        keyargs.get("labelxpos", 0.95),
                        keyargs.get("labelypos", 0.95),
                        ctrstr,
                        horizontalalignment="right",
                        verticalalignment="top",
                        size=size,
                        transform=ax.transAxes,
                        weight=weight,
                    )
            else:
                ax.text(
                    keyargs.get("labelxpos", 0.05),
                    keyargs.get("labelypos", 0.95),
                    ctrstr,
                    horizontalalignment="left",
                    verticalalignment="top",
                    size=size,
                    transform=ax.transAxes,
                    weight=weight,
                )
            # txt.set_bbox(dict(facecolor='white', alpha=1., edgecolor=))
            xlim = None
            ylim = None
            for ctr, plotargs, keyargs_c in self.data[ctrkey]:
                xlim = keyargs_c.pop("xlim", xlim)
                ylim = keyargs_c.pop("ylim", ylim)
                if ctr.isWithError:
                    if maxrelerr is not None:
                        relerr = ctr.err / ctr.sfI
                        mask = relerr > maxrelerr
                        ax.errorbar(
                            ctr.l[~mask],
                            ctr.sfI[~mask],
                            yerr=ctr.err[~mask],
                            **keyargs_c,
                        )
                        if "plotcaps" in keyargs:
                            errbar = ax.errorbar(
                                ctr.l[mask],
                                ctr.sfI[mask],
                                yerr=ctr.err[mask],
                                **keyargs_c,
                            )
                            _, _, barlinecols = errbar.lines
                            barlinecols[0].set_visible(False)
                    else:
                        ax.errorbar(
                            ctr.l, ctr.sfI, yerr=ctr.err, **keyargs_c
                        )  # ,fmt='.',elinewidth=0.5,capsize=1. ,errorevery=1,color='k',zorder=1  # noqa: E501
                else:
                    ax.plot(ctr.l, ctr.sfI, *plotargs, **keyargs_c)
            if self.data[ctrkey][0][0].difference:
                ax.set_yscale("linear")
            else:
                ax.set_yscale("log")
            if xlim is not None:
                self.to_sharex[i] = False
                ax.set_xlim(xlim)
            elif self.xlim is not None:
                ax.set_xlim(self.xlim)

            if ylim is not None:
                self.to_sharey[i] = False
                ax.set_ylim(ylim)
            elif self.xlim is not None:
                ax.set_ylim(self.ylim)

        mixed_quantities = len(set(panel_quantities.values())) > 1
        if mixed_quantities:
            self.to_sharey[:] = False

        if np.any(self.to_sharey):
            shareyaxes = np.array(self.axes)[self.to_sharey]
            for axy in shareyaxes:
                if axy is shareyaxes[0]:
                    continue
                axy.sharey(shareyaxes[0])

        if np.any(self.to_sharex):
            sharexaxes = np.array(self.axes)[self.to_sharex]
            for axx in sharexaxes:
                if axx is sharexaxes[0]:
                    continue
                axx.sharex(sharexaxes[0])

        axes = np.reshape(self.axes, (rows, cols))
        [a.tick_params(axis="x", labelbottom=False) for a in axes[:-1, :].flat]
        [a.tick_params(axis="y", labelleft=False) for a in axes[:, 1:].flat]
        if not self._ylabels_explicit and mixed_quantities:
            for index, quantity in panel_quantities.items():
                self.axes[index].set_ylabel(_QUANTITY_YLABELS[quantity])
                self.axes[index].tick_params(axis="y", labelleft=True)
            if isinstance(self.xlabels, list):
                [a.set_xlabel(xlbl) for xlbl, a in zip(self.xlabels, axes[-1, :])]
            else:
                [a.set_xlabel(self.xlabels) for a in axes[-1, :]]
        elif isinstance(self.xlabels, list):
            [a.set_ylabel(ylbl) for ylbl, a in zip(self.ylabels, axes[:, 0])]
            [a.set_xlabel(xlbl) for xlbl, a in zip(self.xlabels, axes[-1, :])]
        elif keyargs.get("sharexyLabels", False):
            if keyargs.get("sharexyLabels") == "x":
                axes[-1, 0].set_xlabel(self.xlabels)
                [a.set_ylabel(self.ylabels) for a in axes[:, 0]]
                xy_lablel = "x"
            elif keyargs.get("sharexyLabels") == "y":
                axes[0, 0].set_ylabel(self.ylabels)
                [a.set_xlabel(self.xlabels) for a in axes[-1, :]]
                xy_lablel = "y"
            else:
                axes[0, 0].set_ylabel(
                    self.ylabels
                )  # create labels for the tight_layout() call below
                axes[-1, 0].set_xlabel(self.xlabels)
                xy_lablel = "xy"
        else:
            if not self._ylabels_explicit and panel_quantities:
                quantity = next(iter(panel_quantities.values()))
                self.ylabels = _QUANTITY_YLABELS[quantity]
            [a.set_ylabel(self.ylabels) for a in axes[:, 0]]
            [a.set_xlabel(self.xlabels) for a in axes[-1, :]]

        self.tight_layout()
        self.subplots_adjust(wspace=self.wspace, hspace=self.hspace)

        if keyargs.get("sharexyLabels", False):
            if "x" in xy_lablel:
                avepos_width = 0.5 * (self.subplotpars.left + self.subplotpars.right)
                transform_xlabel = mtransforms.blended_transform_factory(
                    self.transFigure, mtransforms.IdentityTransform()
                )
                axes[-1, 0].xaxis.label.set_transform(transform_xlabel)
                axes[-1, 0].xaxis.label.set_x(avepos_width)
            if "y" in xy_lablel:
                avepos_height = 0.5 * (self.subplotpars.bottom + self.subplotpars.top)
                transform_ylabel = mtransforms.blended_transform_factory(
                    mtransforms.IdentityTransform(), self.transFigure
                )
                axes[0, 0].yaxis.label.set_transform(transform_ylabel)
                axes[0, 0].yaxis.label.set_y(avepos_height)

    def get_ctr_ax(self, hk, ident=None):
        hk = (float(hk[0]), float(hk[1]))
        if ident is None:
            idx = self.axes_ctr_hk.index(hk)
        else:
            idx = self.axes_ctr_id.index((hk, ident))
        return self.axes[idx]

    def set_ctr_xlim(self, hk, xlim, ident=None):
        hk = (float(hk[0]), float(hk[1]))
        if self.axes:
            ax = self.get_ctr_ax(hk, ident)
            if len(ax.get_shared_x_axes().get_siblings(ax)) > 1:
                raise RuntimeError("CTR axes limits must be set before plot creation.")
            else:
                ax.set_xlim(xlim)
        else:
            if ident is None:
                for ids in self.data:
                    if ids[0] == hk:
                        self.data[ids][0][2]["xlim"] = xlim
                        break
                else:
                    raise KeyError(f"CTR {str(hk)} not found in figure data")
            else:
                self.data[(hk, ident)][0][2]["xlim"] = xlim

    def set_ctr_ylim(self, hk, ylim, ident=None):
        hk = (float(hk[0]), float(hk[1]))
        if self.axes:
            ax = self.get_ctr_ax(hk, ident)
            if len(ax.get_shared_y_axes().get_siblings(ax)) > 1:
                raise RuntimeError("CTR axes limits must be set before plot creation.")
            else:
                ax.set_ylim(ylim)
        else:
            if ident is None:
                for ids in self.data:
                    if ids[0] == hk:
                        self.data[ids][0][2]["ylim"] = ylim
                        break
                else:
                    raise KeyError(f"CTR {str(hk)} not found in figure data")
            else:
                self.data[(hk, ident)][0][2]["ylim"] = ylim


class CTR:
    ctrtype = 100

    optional_counters = ["bgI", "ctrI", "croi_pix", "bgroi_pix", "weight"]

    def __init__(  # noqa: E741
        self,
        hk,
        l=None,  # noqa: E741
        sfI=None,
        err=None,
        phi=None,
        *,
        reduction=_DEFAULT_REDUCTION,
        scan_geometry=None,
        **keyargs,
    ):
        """Create one measured CTR dataset.

        :param sequence hk: In-plane reference indices in r.l.u.
        :param array-like l: Pointwise reference L in r.l.u.
        :param array-like sfI: Stored F or R values, selected by ``reduction``.
        :param array-like err: Optional uncertainties in the stored quantity.
        :param array-like phi: Optional structure-factor phases in deg.
        :param MeasurementReduction reduction:
            Measurement reduction metadata. Omission supplies the legacy
            structure-factor default with unknown polarization provenance.
            Explicit ``None`` is invalid.
        :param CTRScanGeometry scan_geometry:
            Optional z-mode scan rule; fixed angles are in rad.
        """
        self.hk = tuple(hk)
        h, k = hk
        self.l = np.ascontiguousarray(l)
        self.sfI = np.ascontiguousarray(sfI)
        if err is not None:
            self.err = np.ascontiguousarray(err)
        else:
            self.err = None
        self.ctrtype = CTR.ctrtype
        CTR.ctrtype += 1
        self.withErr = True
        self.phi = phi
        self.difference = False
        self.weight = 1
        if "name" in keyargs:
            self.name = keyargs["name"]
        else:
            self.name = "default"
        self.harr = np.full_like(self.l, h)
        self.karr = np.full_like(self.l, k)
        if reduction is _DEFAULT_REDUCTION:
            reduction = MeasurementReduction()
        self.reduction = reduction
        self.scan_geometry = scan_geometry

    @property
    def reduction(self):
        """Measurement reduction describing the values stored in ``sfI``."""
        return self._reduction

    @reduction.setter
    def reduction(self, reduction):
        if reduction is None or not isinstance(reduction, MeasurementReduction):
            raise TypeError(
                f"{self!r}: reduction must be a MeasurementReduction; "
                "use MeasurementReduction() to reset it."
            )
        polarization = reduction.polarization
        if reduction.quantity == "reflectivity" and self.phi is not None:
            raise ValueError(
                f"{self!r}: reflectivity data cannot carry structure-factor "
                "phase information."
            )
        if polarization is not None:
            factor = polarization.polarization_factor
            if factor is not None and factor.shape != self.l.shape:
                raise ValueError(
                    f"{self!r}: polarization_factor must have the same shape "
                    "as the CTR data."
                )
            if reduction.quantity == "reflectivity" and factor is not None:
                raise ValueError(
                    f"{self!r}: reflectivity data cannot carry a "
                    "polarization_factor; stored R must already be corrected."
                )
        self._reduction = reduction

    @property
    def scan_geometry(self):
        """Optional z-mode scan rule for this measured CTR."""
        return self._scan_geometry

    @scan_geometry.setter
    def scan_geometry(self, geometry):
        if geometry is not None and not isinstance(geometry, CTRScanGeometry):
            raise TypeError(
                f"{self!r}: scan_geometry must be a CTRScanGeometry or None."
            )
        self._scan_geometry = geometry

    def toNXdict(self):
        """Export CTR data with reference HKL and six-circle angles in rad."""
        reduction_group = {
            "@NX_class": "NXcollection",
            "@quantity": self.reduction.quantity,
        }
        if self.reduction.polarization is not None:
            polarization = self.reduction.polarization
            polarization_group = {
                "@NX_class": "NXcollection",
                "@s_fraction": polarization.s_fraction,
                "@outgoing": polarization.outgoing,
            }
            if polarization.polarization_factor is not None:
                polarization_group["polarization_factor"] = (
                    polarization.polarization_factor
                )
                polarization_group["@polarization_factor_unit"] = (
                    "dimensionless"
                )
            reduction_group["polarization"] = polarization_group

        nxdict = {
            "@NX_class": "NXdata",
            "@orgui_ctr_schema": 2,
            "sixc_angles": {"@NX_class": "NXpositioner", "@unit": "rad"},
            "hkl": {
                "@NX_class": "NXcollection",
                "h": self.harr,
                "k": self.karr,
                "l": self.l,
                "@unit": "r.l.u.",
            },
            "counters": {"@NX_class": "NXdetector", "structurefactor": self.sfI},
            "@signal": "counters/structurefactor",
            "@axes": "hkl/l",
            "@title": repr(self),
            "@name": self.name,
            "@difference": self.difference,
            "measurement_reduction": reduction_group,
        }
        if self.scan_geometry is not None:
            geometry_group = {
                "@NX_class": "NXcollection",
                "@fixed": self.scan_geometry.fixed,
                "@mirrorx": self.scan_geometry.mirrorx,
            }
            if self.scan_geometry.angle is not None:
                geometry_group["@angle"] = self.scan_geometry.angle
                geometry_group["@angle_unit"] = "rad"
            nxdict["scan_geometry"] = geometry_group
        for cnter in CTR.optional_counters:
            if hasattr(self, cnter):
                nxdict["counters"][cnter] = getattr(self, cnter)

        if self.err is not None:
            nxdict["counters"]["structurefactor_errors"] = self.err

        if self.phi is not None:
            nxdict["counters"]["phase"] = self.phi

        if hasattr(self, "angles"):
            for ang in self.angles.dtype.fields:
                nxdict["sixc_angles"][ang] = self.angles[ang]

        return nxdict

    @classmethod
    def fromNXdict(cls, nxdict, *, angle_units=None):
        """Load CTR data, storing six-circle angles internally in radians.

        :param dict nxdict: Payload from :meth:`toNXdict` or an external source.
        :param str angle_units:
            Explicit input units, ``"rad"`` or ``"deg"``. By default, schema
            2 and later use the declared ``sixc_angles/@unit``. Unversioned
            payloads retain their numeric angles as radians: old orGUI writers
            incorrectly labeled radians as degrees. Use ``"deg"`` explicitly
            for genuine external degree-valued legacy payloads.
        :raises ValueError: If angle units or angle-array shapes are invalid.

        H and K are preserved as stored. Historical files with K erroneously
        copied from H cannot be repaired without independent index information.
        """
        angle_group = nxdict.get("sixc_angles", {})
        if angle_units is None:
            angle_units = (
                angle_group.get("@unit", "rad")
                if int(nxdict.get("@orgui_ctr_schema", 0)) >= 2
                else "rad"
            )
        if isinstance(angle_units, bytes):
            angle_units = angle_units.decode("ascii")
        if angle_units not in ("rad", "deg"):
            raise ValueError('angle_units must be "rad" or "deg".')
        h = nxdict["hkl"]["h"]
        k = nxdict["hkl"]["k"]
        l = nxdict["hkl"]["l"]  # noqa: E741
        sfI = nxdict["counters"]["structurefactor"]
        err = nxdict["counters"].get("structurefactor_errors", None)
        phi = nxdict["counters"].get("phase", None)
        name = nxdict.get("@name", "default")

        reduction_group = nxdict.get("measurement_reduction")
        if reduction_group is None:
            reduction = MeasurementReduction()
        else:
            quantity = _nx_text(
                reduction_group.get("@quantity", "structure_factor")
            )
            polarization_group = reduction_group.get("polarization")
            if polarization_group is None:
                polarization = None
            else:
                polarization = PolarizationReduction(
                    float(
                        np.asarray(
                            polarization_group["@s_fraction"]
                        ).reshape(-1)[0]
                    ),
                    _nx_text(polarization_group["@outgoing"]),
                    polarization_group.get("polarization_factor"),
                )
            reduction = MeasurementReduction(quantity, polarization)

        geometry_group = nxdict.get("scan_geometry")
        if geometry_group is None:
            scan_geometry = None
        else:
            fixed = _nx_text(geometry_group["@fixed"])
            angle = geometry_group.get("@angle")
            if angle is not None:
                angle = float(np.asarray(angle).reshape(-1)[0])
                units = _nx_text(geometry_group.get("@angle_unit", "rad"))
                if units == "deg":
                    angle = float(np.deg2rad(angle))
                elif units != "rad":
                    raise ValueError(
                        'scan_geometry angle units must be "rad" or "deg".'
                    )
            scan_geometry = CTRScanGeometry(
                fixed=fixed,
                angle=angle,
                mirrorx=_nx_bool(geometry_group.get("@mirrorx", False)),
            )

        ctr = cls(
            (h[0], k[0]),
            l,
            sfI,
            err,
            phi,
            name=name,
            reduction=reduction,
            scan_geometry=scan_geometry,
        )
        ctr.harr = h
        ctr.karr = k

        if "@difference" in nxdict:
            ctr.difference = nxdict["@difference"]

        for cnter in CTR.optional_counters:
            if cnter in nxdict["counters"]:
                setattr(ctr, cnter, nxdict["counters"][cnter])

        angles = []
        angles_names = []
        for ang in angle_group:
            if not ang.startswith("@"):
                values = np.asarray(angle_group[ang], dtype=np.float64)
                if values.shape != ctr.l.shape:
                    raise ValueError(
                        f"Vlieg angle {ang!r} must have the same shape as L."
                    )
                angles.append(np.deg2rad(values) if angle_units == "deg" else values)
                angles_names.append(ang)
        if angles:
            dt = np.dtype([(ang, "f8") for ang in angles_names])
            angles = np.core.records.fromarrays(angles, dtype=dt)
            ctr.angles = angles
        return ctr

    def getPlotLabel(self):
        if self.plotlabel in self:
            return self.plotlabel
        else:
            return self.name

    def _select_points(self, selection):
        """Apply one point selection to every aligned CTR field."""
        old_size = self.l.size
        selected_reduction = _select_reduction(self.reduction, selection)
        for attribute in (
            "l",
            "harr",
            "karr",
            "sfI",
            "err",
            "phi",
            "angles",
            "bgI",
            "ctrI",
            "croi_pix",
            "bgroi_pix",
        ):
            value = getattr(self, attribute, None)
            if isinstance(value, np.ndarray) and value.ndim > 0:
                if value.shape[0] == old_size:
                    setattr(self, attribute, value[selection])
        self.reduction = selected_reduction

    def convertToF(self, excludeInvalid=True):
        """Convert corrected signed intensity to signed amplitude in place.

        Input values are intensities in F-squared arbitrary units. The central
        value becomes ``sign(I) * sqrt(abs(I))`` in structure-factor arbitrary
        units. If uncertainties are present, each becomes half the interval
        between the signed-square-root transforms of ``I - sigma_I`` and
        ``I + sigma_I``. This effective symmetric uncertainty is finite at
        zero and is not an exact Normal uncertainty near zero.

        :param bool excludeInvalid:
            Remove points whose central intensity is not finite. Finite
            negative and zero intensities are always retained.
        :raises ValueError:
            If called on reflectivity, or if a retained intensity has a
            nonfinite, nonpositive, or misaligned uncertainty.
        """
        self._require_structure_factor("intensity-to-amplitude conversion")
        intensity = np.asarray(self.sfI, dtype=np.float64)
        mask = (
            np.isfinite(intensity)
            if excludeInvalid
            else np.ones_like(intensity, dtype=np.bool_)
        )
        amplitude = _signed_sqrt(intensity)
        uncertainty = None
        if self.isWithError:
            uncertainty = _signed_sqrt_uncertainty(intensity, self.err)

        auxiliary = {}
        for name in ("bgI", "ctrI"):
            if hasattr(self, name):
                auxiliary[name] = _signed_sqrt(getattr(self, name))

        self.sfI = amplitude
        if uncertainty is not None:
            self.err = uncertainty
        for name, values in auxiliary.items():
            setattr(self, name, values)
        self._select_points(mask)

    # in degrees
    def setPhase(self, phi):
        self._require_structure_factor("structure-factor phase assignment")
        self.phi = phi

    def setWithError(self, err):
        self.withErr = err

    @property
    def isWithPhase(self):
        if not isinstance(self.phi, np.ndarray):
            return False
        else:
            return True

    @property
    def isWithError(self):
        if not isinstance(self.err, np.ndarray):
            return False
        return self.withErr

    def getComplexSF(self):
        self._require_structure_factor("complex structure factors")
        if not self.isWithPhase:
            raise Exception(f"{repr(self)}:\nNo phase informaion available.")
        return self.sfI * np.exp(1j * np.deg2rad(self.phi))

    @property
    def ctr_id(self):
        return tuple(np.around(self.hk, 2)), self.ctrtype

    def _require_structure_factor(self, operation):
        """Reject this CTR when an operation assumes structure factors."""
        if self.reduction.quantity != "structure_factor":
            raise ValueError(
                f"{self!r}: {operation} supports structure-factor data only; "
                "reflectivity requires an explicit supported workflow."
            )

    def setToDefaultID(self):
        self.ctrtype = 0

    def generateDifference(self, other):
        self._require_structure_factor("CTR differences")
        other._require_structure_factor("CTR differences")
        otherinter = interp.interp1d(other.l, other.sfI)
        self.sfI -= otherinter(self.l)
        self.difference = True
        return self

    def meanSF(self, lowerL, upperL):
        upper = np.nanargmin(np.abs(self.l - upperL))
        lower = np.nanargmin(np.abs(self.l - lowerL))
        return np.nanmean(self.sfI[lower:upper])

    def __imul__(self, valOrArray):
        self.sfI *= valOrArray
        if isinstance(self.err, np.ndarray):
            self.err *= valOrArray
        return self

    def __iadd__(self, valOrArray):
        self.sfI += valOrArray
        # raise NotImplementedError()
        # if isinstance(self.err,np.ndarray):
        #    self.err += valOrArray
        return self

    def cut(self, lower, upper, invert=False):
        """Restricts the CTR to the selected lower and upper index

        If invert is True, will remove the data within the selected range.

        If invert is set to 'insertNAN', the sfI of the lower index will be set to nan.
        This is useful for CTR plotting to interrupt the lines at this point.
        But should never be used for any CTR that is supposed to be used for
        further computations!

        """
        if invert:
            mask = np.ones_like(self.l, dtype=np.bool_)
            if invert == "insertNAN":
                mask[lower + 1 : upper] = False
                self.sfI[lower] = np.nan
            else:
                mask[lower:upper] = False

        else:
            mask = slice(lower, upper)
        self._select_points(mask)

    def cutToL(self, lowerL, upperL, invert=False):
        """Restricts the CTR to the selected lowerL and upperL.

        If invert is True, will remove the data within the selected range.

        If invert is set to 'insertNAN', the sfI of the lower index will be set to nan.
        This is useful for CTR plotting to interrupt the lines at this point.
        But should never be used for any CTR that is supposed to be used for
        further computations!
        """
        upper = np.nanargmin(np.abs(self.l - upperL))
        lower = np.nanargmin(np.abs(self.l - lowerL))
        self.cut(lower, upper, invert)

    def cutToROIfile(self, jsonfile):
        with open(jsonfile) as f:
            udict = json.load(f)

        rois = udict["ROI"]["roidict"]
        del rois["ICR"]
        mask = np.zeros_like(self.l, dtype=np.bool_)

        for roikey in rois:
            fr = rois[roikey]["from"]
            to = rois[roikey]["to"]
            lower = np.nanargmin(np.abs(self.l - fr))
            upper = np.nanargmin(np.abs(self.l - to))
            mask[lower:upper] = 1.0
        self._select_points(mask)

    def get_scale(self, xtal, omitErrors=False, lognorm=False):
        self._require_structure_factor("kinematical crystal scaling")
        if not hasattr(self, "err") or omitErrors:
            err = None
        else:
            err = self.err
        F_cryst = np.abs(xtal.F(self.harr, self.karr, self.l))
        if lognorm:
            return util.get_scale_logchi2(F_cryst, self.sfI)
        else:
            return util.get_scale_chi2(F_cryst, self.sfI, err)

    def scaleToXtal(self, xtal, omitErrors=False, lognorm=False):
        self._require_structure_factor("kinematical crystal scaling")
        self.__imul__(self.get_scale(xtal, omitErrors, lognorm))

    def calcAnglesZmode(
        self,
        vliegangles,
        fixedangle=np.deg2rad(0.1),
        fixed="in",
        chi=0.0,
        phi=0.0,
        *,
        hkl_transform=None,
        **keyargs,
    ):
        """Calculate and store Vlieg z-mode angles in radians.

        :param HKLVlieg.VliegAngles vliegangles:
            Angle calculator configured for the target lattice.
        :param float fixedangle:
            Fixed incidence or exit angle in rad.
        :param str fixed:
            Angle constraint: ``"in"``, ``"out"``, or ``"eq"``.
        :param float chi:
            Fixed chi angle in rad.
        :param float phi:
            Fixed phi angle in rad.
        :param hkl_transform:
            Optional finite, nonsingular 3-by-3 matrix mapping CTR reference
            HKL to the calculator's HKL (both in r.l.u.). The default leaves
            HKL unchanged. For a bulk-lattice calculator, pass the crystal's
            ``uc_bulk.refHKLTransform``.
        :returns:
            Structured records containing alpha, delta, gamma, omega, chi,
            and phi in rad.
        :rtype: numpy.recarray
        """
        try:
            self.angles = _calculate_angles_zmode(
                self.harr,
                self.karr,
                self.l,
                vliegangles,
                fixedangle,
                fixed=fixed,
                chi=chi,
                phi=phi,
                hkl_transform=hkl_transform,
                **keyargs,
            )
        except (TypeError, ValueError) as error:
            raise type(error)(f"{self!r}: {error}") from error
        return self.angles

    def toArray(self, mode=None):
        if mode is not None:
            data = np.empty((6, self.l.size))
            data[5] = mode
        else:
            data = np.empty((5, self.l.size))
        data[0] = self.harr
        data[1] = self.karr
        data[2] = self.l
        data[3] = self.sfI
        if self.err is not None:
            data[4] = self.err
        elif self.phi is not None:
            data[4] = self.phi
        else:
            data[4] = np.nan

        return data.T

    # returns a list of CTRs!!!
    @staticmethod
    def fromANAROD(
        filenameOrArray,
        RODexport=False,
        *,
        reduction=_DEFAULT_REDUCTION,
        scan_geometry=None,
    ):
        warnings.warn(
            "CTR.fromANAROD is deprecated, use CTRCollection.fromANAROD instead!",
            DeprecationWarning,
        )
        if not isinstance(filenameOrArray, np.ndarray):
            filenameOrArray = np.atleast_2d(
                np.loadtxt(filenameOrArray, skiprows=1)
            )
        else:
            filenameOrArray = np.atleast_2d(filenameOrArray)
        rods = np.unique(filenameOrArray[:, :2], axis=0)
        CTRs = []
        for hk in rods:
            rodmask = np.logical_and(
                filenameOrArray[:, 0] == hk[0], filenameOrArray[:, 1] == hk[1]
            )
            rod = filenameOrArray[rodmask]
            l = rod[:, 2]  # noqa: E741
            sfI = rod[:, 3]
            rod_reduction = _reduction_for_import(reduction, rodmask)
            if RODexport:
                ctr = CTR(
                    tuple(hk),
                    l,
                    sfI,
                    reduction=rod_reduction,
                    scan_geometry=scan_geometry,
                )
                ctr.setPhase(rod[:, 4])
                CTRs.append(ctr)
            else:
                err = rod[:, 4] if rod.shape[1] > 4 else None
                CTRs.append(
                    CTR(
                        tuple(hk),
                        l,
                        sfI,
                        err,
                        reduction=rod_reduction,
                        scan_geometry=scan_geometry,
                    )
                )

        return CTRCollection(CTRs)

    @classmethod
    def fromArray(
        cls,
        array,
        RODexport=False,
        *,
        reduction=_DEFAULT_REDUCTION,
        scan_geometry=None,
    ):
        h = array[:, 0][0]
        k = array[:, 1][0]
        l = array[:, 2]  # noqa: E741
        sfI = array[:, 3]
        if RODexport:
            err = None
            phase = array[:, 4]
        else:
            err = array[:, 4] if array.shape[1] > 4 else None
            phase = None
        return cls(
            [h, k],
            l,
            sfI,
            err,
            phase,
            reduction=reduction,
            scan_geometry=scan_geometry,
        )

    def millerIdentifier(self):
        h, k = self.hk
        idstr = f"{round(h):d}_{round((h % 1) * 100):02d}_{round(k):d}_{round((k % 1) * 100):02d}"  # noqa: E501
        return idstr

    def rollbackHK(self):
        idstrsp = self.name.split("_")
        h = float(idstrsp[0] + "." + idstrsp[1])
        k = float(idstrsp[2] + "." + idstrsp[3])
        self.hk = h, k

    def generateAverage(self, step_size=None, **kwargs):
        """Creates a new averaged CTR.

        args: step_size:  step size along l
              nbins : number of bins along l
              overlap (float from 0 to 1): define how much the first and last bins exceed the data range.
                                           This can solve issues with edge data points
        provide either step_size or nbins

        """  # noqa: E501
        if (
            self.reduction != MeasurementReduction()
            or self.scan_geometry is not None
            or hasattr(self, "angles")
        ):
            raise NotImplementedError(
                f"{self!r}: averaging CTR geometry or measurement-reduction "
                "metadata is not supported."
            )
        overlap = kwargs.get("overlap", 0.25)

        lmax = np.amax(self.l)
        lmin = np.amin(self.l)
        size_max = self.l.size

        l_range = lmax - lmin

        if "nbins" in kwargs:
            nbins = int(kwargs.get("nbins"))
            #
            l_full_range = l_range / (1.0 - (2 * overlap) / nbins)
            step = l_full_range / nbins
            l_first_bin = lmin - step * overlap
            bin_edges = l_first_bin + step * np.arange(nbins + 1)

        elif step_size is not None:
            nbins = int(np.floor(l_range / abs(step_size))) + 1
            l_full_range = l_range / (1.0 - (2 * overlap) / nbins)

            nbins = int(np.floor(l_full_range / abs(step_size))) + 1

            l_first_bin = lmin - step_size * overlap
            bin_edges = l_first_bin + step_size * np.arange(nbins + 1)
        else:
            nbins = size_max + 1

            l_full_range = l_range / (1.0 - (2 * overlap) / nbins)
            step = l_full_range / nbins

            l_first_bin = lmin - step * overlap
            bin_edges = l_first_bin + step * np.arange(nbins + 1)

        l_cntr = np.zeros(nbins)

        indexes = np.digitize(self.l, bin_edges)

        if np.any(indexes == 0) or np.any(indexes == nbins + 1):
            raise Exception(
                "bin edges were chosen incorrectly. This is probably a bug."
            )

        indexes -= 1

        weights = np.zeros_like(l_cntr)
        np.add.at(weights, indexes, 1.0)

        np.add.at(l_cntr, indexes, self.l)
        l_cntr /= weights

        I = np.zeros_like(l_cntr)  # noqa: E741
        np.add.at(I, indexes, self.sfI)
        I /= weights  # noqa: E741

        if hasattr(self, "err") and self.err is not None:
            Ierr = np.zeros_like(l_cntr)
            np.add.at(Ierr, indexes, self.err**2)
            Ierr = np.sqrt(Ierr)
            Ierr /= weights
        else:
            Ierr = None

        mask = np.logical_and(weights != 0, np.isfinite(I))

        l_masked = l_cntr[mask]
        I_masked = I[mask]
        if Ierr is not None:
            Ierr_masked = Ierr[mask]
        else:
            Ierr_masked = None
        weights_masked = weights[mask]

        newctr = CTR(self.hk, l_masked, I_masked, Ierr_masked)
        newctr.contributions = weights_masked
        return newctr

    def __repr__(self):
        # return "<CTR %s ctrtype %s at %016X>" % (tuple(np.around(self.hk,2)), self.ctrtype , id(self))  # noqa: E501
        return (
            f"<CTR<{self.name}> {tuple(np.around(self.hk, 2))} ctrtype {self.ctrtype}>"  # noqa: E501
        )


class CTRCollection(list):
    def __init__(self, iterableCTRs=[], **kwargs):
        super().__init__(iterableCTRs)
        self.plotsett = ()
        self.plotkeyargs = {"linestyle": "-", "marker": ".", "color": "k", "zorder": 1}
        # self._updaterodtypes()
        self.name = kwargs.get("name", "CTRs")

    def setPlotSettings(self, *settings, **keyargs):
        self.plotsett = settings
        self.plotkeyargs = keyargs

    def getReprList(self):
        return [repr(rod) for rod in self]

    def getHKList(self):
        return [rod.hk for rod in self]

    def setAllToDefaultID(self):
        for rod in self:
            rod.setToDefaultID()

    def _require_structure_factors(self, operation):
        """Reject collections containing reflectivity for legacy operations."""
        for rod in self:
            rod._require_structure_factor(operation)

    def deleteRod(self, key):
        if isinstance(key, tuple):
            for rod in self:
                if rod.hk == key:
                    self.remove(rod)

    def __and__(self, other):
        coll = CTRCollection()
        for rod in self:
            if repr(rod) in other.getReprList():
                coll.append(rod)
        # coll._updaterodtypes()
        return coll

    def convertToF(self, excludeInvalid=True):
        """Convert every member's signed intensity to amplitude in place."""
        for ctr in self:
            ctr.convertToF(excludeInvalid=excludeInvalid)

    def calcAnglesZmode(
        self,
        vliegangles,
        fixedangle=np.deg2rad(0.1),
        fixed="in",
        chi=0.0,
        phi=0.0,
        *,
        hkl_transform=None,
        **keyargs,
    ):
        """Calculate and store Vlieg z-mode angles for every CTR.

        All angle arguments and returned angle records are in rad.

        :param HKLVlieg.VliegAngles vliegangles:
            Angle calculator configured for the target lattice.
        :param float fixedangle:
            Fixed incidence or exit angle in rad.
        :param str fixed:
            Angle constraint: ``"in"``, ``"out"``, or ``"eq"``.
        :param float chi:
            Fixed chi angle in rad.
        :param float phi:
            Fixed phi angle in rad.
        :param hkl_transform:
            Optional reference-to-calculator HKL matrix, forwarded to each
            :meth:`CTR.calcAnglesZmode`. The default leaves HKL unchanged.
        :returns:
            One structured angle-record array per CTR.
        :rtype: list[numpy.recarray]
        """
        return [
            ctr.calcAnglesZmode(
                vliegangles,
                fixedangle=fixedangle,
                fixed=fixed,
                chi=chi,
                phi=phi,
                hkl_transform=hkl_transform,
                **keyargs,
            )
            for ctr in self
        ]

    def generateAverage(self, step_size=None, **kwargs):
        """Creates a new CTRCollection with CTRs, which were individually
        averaged.

        args: step_size:  step size along l
              nbins : number of bins along l
              overlap (float from 0 to 1): define how much the first and last bins exceed the data range.
                                           This can solve issues with edge data points
        provide either step_size or nbins

        """  # noqa: E501
        coll = CTRCollection(name="AVE: " + self.name)
        for ctr in self:
            coll.append(ctr.generateAverage(step_size, **kwargs))
        return coll

    def generateDifferenceCollection(self, other, sortby="repr"):
        self._require_structure_factors("CTR differences")
        coll = CTRCollection()
        if isinstance(other, CTRCollection):
            other._require_structure_factors("CTR differences")
            for rod in self:
                if sortby == "repr":
                    try:
                        idx = other.getReprList().index(repr(rod))
                    except ValueError:
                        continue
                else:
                    try:
                        idx = other.getHKList().index(rod.hk)
                    except ValueError:
                        continue
                coll.append(rod.generateDifference(other[idx]))
            return coll
        elif isinstance(other, SXRDCrystal):
            for ctr in self:
                l = np.copy(ctr.l)  # noqa: E741
                h, k = ctr.hk
                harr = np.full_like(l, h)
                karr = np.full_like(l, k)
                F_theo = np.abs(other.F(harr, karr, l))
                scale = np.prod(ctr.sfI / F_theo) ** (1 / ctr.l.size)
                F_theo *= scale
                diffCTR = CTR((h, k), l, ctr.sfI - F_theo, name=ctr.name)
                diffCTR.difference = True
                coll.append(diffCTR)
            return coll
        else:
            raise NotImplementedError(
                f"Can not generate difference collection for type {type(other)}"
            )

    def addScalar(self, scalar):
        for rod in self:
            rod += scalar
        return self

    def generateCollectionFromXtal(self, xtal, samples=None, lrange=None):
        CTRs = CTRCollection()
        if samples is not None or lrange is not None:
            if lrange is not None:
                for ctr in self:
                    h, k = ctr.hk
                    l = np.linspace(lrange[0], lrange[1], samples)  # noqa: E741
                    h = np.full_like(l, h)
                    k = np.full_like(l, k)
                    F = xtal.F(h, k, l)
                    CTRs.append(CTR(ctr.hk, l, np.abs(F), phi=np.angle(F)))
                return CTRs
            else:
                for ctr in self:
                    h, k = ctr.hk
                    l = np.linspace(np.amin(ctr.l), np.amax(ctr.l), samples)  # noqa: E741
                    h = np.full_like(l, h)
                    k = np.full_like(l, k)
                    F = xtal.F(h, k, l)
                    CTRs.append(CTR(ctr.hk, l, np.abs(F), phi=np.angle(F)))
                return CTRs
        else:
            for ctr in self:
                F = xtal.F(ctr.harr, ctr.karr, ctr.l)
                CTRs.append(CTR(ctr.hk, ctr.l, np.abs(F), phi=np.angle(F)))
            return CTRs

    def toANAROD(self, filename, mode=-3):
        self._require_structure_factors("ANAROD structure-factor export")
        if self.__getitem__(0).isWithError:
            header = "H  K  L  F_HKL  errorF  mode"
        else:
            header = "H  K  L  F_HKL  phi"

        data_combined = np.vstack([ctr.toArray(mode) for ctr in self])

        np.savetxt(filename, data_combined, header=header, fmt="%.5f")

    @staticmethod
    def fromCTRFile(filename, **kwargs):
        """Load a metadata-aware orGUI CTR file.

        The file uses whitespace-separated H, K, L, value, and uncertainty
        columns followed by an optional ``polarization_factor``. It deliberately
        omits ANAROD's unused mode column and adds a schema-2 comment header.
        Header metadata is applied to every rod, while the pointwise factor is
        split with the corresponding rod values.

        :param path-like filename: File to load.
        :param str name: Optional collection name.
        :param MeasurementReduction reduction:
            Optional explicit override for the stored reduction metadata.
        :param CTRScanGeometry scan_geometry:
            Optional scan geometry applied to every imported rod.
        :returns: Loaded CTR collection.
        :rtype: CTRCollection
        """
        reduction = kwargs.pop("reduction", _DEFAULT_REDUCTION)
        scan_geometry = kwargs.pop("scan_geometry", None)
        name = kwargs.pop("name", os.path.basename(filename))
        if kwargs:
            unexpected = ", ".join(sorted(kwargs))
            raise TypeError(f"Unexpected CTR file option(s): {unexpected}")
        table, reduction = _load_ctr_file(filename, reduction)
        return CTRCollection.fromANAROD(
            table,
            reduction=reduction,
            scan_geometry=scan_geometry,
            name=name,
        )

    @staticmethod
    def fromANAROD(filenameOrArray, RODexport=False, **kwargs):
        reduction = kwargs.pop("reduction", _DEFAULT_REDUCTION)
        scan_geometry = kwargs.pop("scan_geometry", None)
        if not isinstance(filenameOrArray, np.ndarray):
            name = kwargs.get("name", os.path.basename(filenameOrArray))
            filenameOrArray = np.atleast_2d(
                np.loadtxt(filenameOrArray, skiprows=1)
            )
        else:
            name = kwargs.get("name", "Array_CTR_import")
            filenameOrArray = np.atleast_2d(filenameOrArray)
        rods = np.unique(filenameOrArray[:, :2], axis=0)
        CTRs = []
        for hk in rods:
            rodmask = np.logical_and(
                filenameOrArray[:, 0] == hk[0], filenameOrArray[:, 1] == hk[1]
            )
            rod = filenameOrArray[rodmask]
            l = rod[:, 2]  # noqa: E741
            sfI = rod[:, 3]
            rod_reduction = _reduction_for_import(reduction, rodmask)
            if RODexport:
                ctr = CTR(
                    tuple(hk),
                    l,
                    sfI,
                    reduction=rod_reduction,
                    scan_geometry=scan_geometry,
                )
                ctr.setPhase(rod[:, 4])
                CTRs.append(ctr)
            else:
                err = rod[:, 4] if rod.shape[1] > 4 else None
                CTRs.append(
                    CTR(
                        tuple(hk),
                        l,
                        sfI,
                        err,
                        reduction=rod_reduction,
                        scan_geometry=scan_geometry,
                    )
                )

        return CTRCollection(CTRs, name=name)

    def toNXdict(self):
        nxdict = {
            "@title": f"{self.name}",
            "@NX_class": "NXentry",
            "@name": f"{self.name}",
        }
        for ctr in self:
            d = ctr.toNXdict()
            nxdict[d["@title"]] = d
            nxdict["@default"] = d["@title"]

        return nxdict

    @classmethod
    def fromNXdict(cls, nxdict, *, angle_units=None):
        """Load a collection, forwarding input angle units to every CTR.

        :param str angle_units:
            ``"rad"``, ``"deg"``, or ``None``; see :meth:`CTR.fromNXdict`
            for schema-aware defaults and legacy radian compatibility.
        """
        ctrs = []
        for dt in nxdict:
            if not dt.startswith("@"):
                ctr = CTR.fromNXdict(nxdict[dt], angle_units=angle_units)
                ctrs.append(ctr)
        name = nxdict["@name"]

        CTRs = cls(ctrs, name=name)
        return CTRs

    def get_flat(self):
        h = []
        k = []
        l = []  # noqa: E741
        F = []
        for ctr in self:
            h.append(ctr.harr)
            k.append(ctr.karr)
            l.append(ctr.l)
            F.append(ctr.sfI)
        return (
            np.concatenate(h),
            np.concatenate(k),
            np.concatenate(l),
        ), np.concatenate(F)

    def get_err_flat(self):
        err = []
        for ctr in self:
            err.append(ctr.err)
        return np.concatenate(err)

    def get_scale(self, xtal, omitErrors=False, lognorm=False):
        self._require_structure_factors("kinematical crystal scaling")
        hkl, F = self.get_flat()
        F_cryst = np.abs(xtal.F(*hkl))
        if omitErrors:
            err = None
        else:
            err = self.get_err_flat()
        if lognorm:
            return util.get_scale_logchi2(F_cryst, F)
        else:
            return util.get_scale_chi2(F_cryst, F, err)

    def scaleToXtal(self, xtal, individual=True, omitErrors=False, lognorm=False):
        self._require_structure_factors("kinematical crystal scaling")
        if individual:
            for ctr in self:
                ctr.scaleToXtal(xtal, lognorm, omitErrors)
        else:
            self.__imul__(self.get_scale(xtal, omitErrors, lognorm))

    def __imul__(self, valOrArray):
        self._require_structure_factors("collection scaling")
        for rod in self:
            rod *= valOrArray
        return self

    def __getitem__(self, key):
        if isinstance(key, tuple):
            for rod in self:
                if rod.hk == key:
                    return rod
            else:
                raise KeyError(
                    f"<{type(self).__name__}: {self.name}> CTR indices {str(key)} not found."  # noqa: E501
                )
        return super().__getitem__(key)

    def __repr__(self):
        s = f"<{type(self).__name__}: {self.name}\n"
        s += super().__repr__()
        s += ">"
        return s


def ctrfigure(**figargs):
    figargs["FigureClass"] = CTRFigure
    fig = plt.figure(**figargs)
    return fig


if __name__ == "__main__":
    fig = ctrfigure(figsize=(12, 8))
    ctrlist = CTR.fromANAROD("CTRfit/data_in/0V17/CH4977_0V17_0001.dat")
    [ctr.setWithError(False) for ctr in ctrlist]
    [ctr.setToDefaultID() for ctr in ctrlist]
    [fig.addCTR(r, linestyle="", marker=".", color="k", zorder=1) for r in ctrlist]
    ctrlist = CTR.fromANAROD("CTRfit/data_in/0V47/CH4977_0V47_0001.dat")
    [ctr.setWithError(False) for ctr in ctrlist]
    [ctr.setToDefaultID() for ctr in ctrlist]
    [fig.addCTR(r, linestyle="", marker=".", color="g", zorder=1) for r in ctrlist]

    ctrlist = CTR.fromANAROD("CTRfit/data_in/0V72/CH4977_0V72_0001.dat")
    [ctr.setWithError(False) for ctr in ctrlist]
    [ctr.setToDefaultID() for ctr in ctrlist]
    [fig.addCTR(r, linestyle="", marker=".", color="b", zorder=1) for r in ctrlist]

    ctrlist = CTR.fromANAROD("CTRfit/data_in/0V02/CH4977_0V02_0001.dat")
    [ctr.setWithError(False) for ctr in ctrlist]
    [ctr.setToDefaultID() for ctr in ctrlist]
    [fig.addCTR(r, linestyle="", marker=".", color="y", zorder=1) for r in ctrlist]

    fig.settings(wspace=0.05, hspace=0, ylabels="|F| / a.u.", ylim=[1e1, 1e3])
    fig.generateCTRplot()
    fig.show()
