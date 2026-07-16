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
__version__ = "1.2.0"
__maintainer__ = "Timo Fuchs"
__email__ = "tfuchs@cornell.edu"

import numpy as np
from .. import util
import re
# random.seed(45)


from .CTRutil import _ensure_contiguous, next_skip_comment, LinearFitFunctions
from .CTRoptics import combine_profiles

from .CTRuc import UnitCell, ctr_accel_enabled
from .CTRdistributions import (
    DEFAULT_TAIL_PROBABILITY,
    PoissonProfile,
    SkellamProfile,
)
from .CTRstacking import (
    LayerCycle,
    LayerState,
    LayerTransition,
    resolve_upper_start,
)


def _cyclic_layer_order(layer_ids, below_layer):
    layer_ids = np.asarray(layer_ids)
    if below_layer == -1:
        start = 0
    else:
        matches = np.flatnonzero(layer_ids == below_layer)
        if matches.size == 0:
            raise ValueError(
                f"Layer {below_layer} is not present in cyclic layer sequence {layer_ids.tolist()}."  # noqa: E501
            )
        start = (matches[0] + 1) % len(layer_ids)
    indices = np.roll(np.arange(len(layer_ids)), -start)
    return layer_ids[indices], indices


def _unwrapped_layer_positions(layer_positions, order):
    positions = np.asarray(layer_positions)[order].astype(np.float64)
    if positions.size == 0:
        return positions
    positions = np.copy(positions)
    for i in range(1, len(positions)):
        while positions[i] <= positions[i - 1]:
            positions[i] += 1.0
    return positions


def _translate_domains(layers, height):
    for layer in layers:
        offset = height / layer.a[2]
        for matrix in layer.coherentDomainMatrix:
            matrix[2, 3] += offset


def _move_domains(layers, delta_height):
    for layer in layers:
        offset = delta_height / layer.a[2]
        for matrix in layer.coherentDomainMatrix:
            matrix[2, 3] += offset


class _LayerStackingMixin:
    def _initialize_layer_stacking(self, kwargs):
        self.layer_behavior = kwargs.get(
            "layer_behavior",
            kwargs.get("layer_behaviour", "select"),
        )
        transition = kwargs.get("layer_transition")
        if transition is not None and not isinstance(transition, LayerTransition):
            transition = LayerTransition(transition)
        self.layer_transition = transition
        self._start_layer_number = -1.0
        self.start_layer_number = -1.0

    @property
    def layer_behavior(self):
        """Return whether cyclic layer selection is applied."""
        return self._layer_behavior

    @layer_behavior.setter
    def layer_behavior(self, behavior):
        if behavior not in ("ignore", "select"):
            raise ValueError(
                f"layer_behavior must be 'ignore' or 'select', not {behavior!r}"
            )
        self._layer_behavior = behavior

    @property
    def layer_behaviour(self):
        """Return the legacy spelling of :attr:`layer_behavior`."""
        return self.layer_behavior

    @layer_behaviour.setter
    def layer_behaviour(self, behavior):
        self.layer_behavior = behavior

    def _wyckoff_target_unitcells(self, kwargs):
        if hasattr(self, "uc_top") and hasattr(self, "uc_bottom"):
            if "unitcell" not in kwargs:
                raise ValueError(
                    "Missing unit cell name. Provide unit cell name as kwarg "
                    "'unitcell'"
                )
            unitcell = kwargs["unitcell"]
            if isinstance(unitcell, list | tuple):
                return [self[ucn] for ucn in unitcell]
            return [self[unitcell]]
        return [self.unitcell]

    def _forward_wyckoff_call(self, method_name, *args, **kwargs):
        parameters = [
            getattr(unitcell, method_name)(*args, **kwargs)
            for unitcell in self._wyckoff_target_unitcells(kwargs)
        ]
        return parameters if len(parameters) > 1 else parameters[0]

    def addWyckoffParameter(
        self,
        site_id,
        variable,
        limits=(-np.inf, np.inf),
        absolute_limits=None,
        **kwargs,
    ):
        """Add an absolute Wyckoff variable parameter on contained unit cells.

        :param str site_id:
            Wyckoff site identifier.
        :param str variable:
            Wyckoff variable name, for example ``"u"``.
        :param tuple limits:
            Absolute variable limits in parent fractional units.
        :param tuple absolute_limits:
            Deprecated alias for ``limits`` retained for compatibility.
        :param kwargs:
            For ``EpitaxyInterface``, provide ``unitcell="top"``,
            ``unitcell="bottom"``, or a list of these names.
        :returns:
            Created parameter or list of parameters.
        """
        return self._forward_wyckoff_call(
            "addWyckoffParameter",
            site_id,
            variable,
            limits=limits,
            absolute_limits=absolute_limits,
            **kwargs,
        )

    def addWyckoffParameters(
        self,
        site_id,
        variables=None,
        limits=(-np.inf, np.inf),
        absolute_limits=None,
        **kwargs,
    ):
        """Add absolute Wyckoff variable parameters on contained unit cells.

        :param str site_id:
            Wyckoff site identifier.
        :param variables:
            Iterable of Wyckoff variable names. If ``None``, fit all variables.
        :param tuple limits:
            Absolute variable limits in parent fractional units, or a mapping
            from variable names to absolute limits.
        :param absolute_limits:
            Deprecated alias for ``limits`` retained for compatibility.
        :param kwargs:
            For ``EpitaxyInterface``, provide ``unitcell="top"``,
            ``unitcell="bottom"``, or a list of these names.
        :returns:
            Created parameters. A list of lists is returned when multiple
            contained unit cells are selected.
        """
        return self._forward_wyckoff_call(
            "addWyckoffParameters",
            site_id,
            variables=variables,
            limits=limits,
            absolute_limits=absolute_limits,
            **kwargs,
        )

    def addWyckoffShift(
        self,
        site_id,
        axis,
        limits=(-np.inf, np.inf),
        absolute_limits=None,
        **kwargs,
    ):
        """Add a relative symmetry-site shift along a parent direction.

        The parameter is a displacement from the representative symmetry-site
        coordinate, not an absolute parent coordinate.

        :param str site_id:
            Wyckoff site identifier.
        :param str axis:
            Parent conventional representative coordinate, one of ``"x"``,
            ``"y"``, or ``"z"``.
        :param tuple limits:
            Delta fit limits in parent fractional units.
        :param tuple absolute_limits:
            Optional absolute parent-coordinate bounds converted to relative
            shift bounds. The fitted value remains a relative displacement.
        :param kwargs:
            For ``EpitaxyInterface``, provide ``unitcell="top"``,
            ``unitcell="bottom"``, or a list of these names.
        :returns:
            Created parameter or list of parameters.
        """
        return self._forward_wyckoff_call(
            "addWyckoffShift",
            site_id,
            axis,
            limits=limits,
            absolute_limits=absolute_limits,
            **kwargs,
        )

    def addWyckoffShifts(
        self,
        site_id,
        axes=("x", "y", "z"),
        limits=(-np.inf, np.inf),
        absolute_limits=None,
        **kwargs,
    ):
        """Add relative symmetry-site shifts along parent directions.

        Every parameter is a displacement from the representative
        symmetry-site coordinate, not an absolute parent coordinate.

        :param str site_id:
            Wyckoff site identifier.
        :param axes:
            Parent conventional representative coordinate axes.
        :param tuple limits:
            Delta fit limits in parent fractional units.
        :param absolute_limits:
            Optional absolute parent-coordinate bounds converted to relative
            shift bounds. The fitted values remain relative displacements.
        :param kwargs:
            For ``EpitaxyInterface``, provide ``unitcell="top"``,
            ``unitcell="bottom"``, or a list of these names.
        :returns:
            Created parameters. A list of lists is returned when multiple
            contained unit cells are selected.
        """
        return self._forward_wyckoff_call(
            "addWyckoffShifts",
            site_id,
            axes=axes,
            limits=limits,
            absolute_limits=absolute_limits,
            **kwargs,
        )

    @property
    def start_layer_number(self):
        """Return the first cyclic layer created by this object."""
        return self._start_layer_number

    @start_layer_number.setter
    def start_layer_number(self, below_layer):
        if self.layer_behavior == "ignore":
            below_layer = -1
        layer_order, order = _cyclic_layer_order(self._layer_ids, below_layer)
        self._start_layer_number = layer_order[0]
        self._set_layer_order(layer_order, order)
        self._basis_created = np.full_like(self.basis, np.nan)

    @property
    def layer_cycle(self):
        """Return the ordered local structural-layer cycle."""
        return LayerCycle(self._layer_ids)

    @property
    def layer_state(self):
        """Return the top structural-layer state of this object."""
        return LayerState(self.layer_cycle, self.end_layer_number)

    @property
    def stacking_height_absolute(self):
        """Return the nominal height passed to the object above."""
        return self.height_absolute

    @property
    def stacking_loc_absolute(self):
        """Return the nominal reference location passed upward."""
        return self.loc_absolute

    def _set_start_layer(self, start_layer):
        matches = np.flatnonzero(self._layer_ids == start_layer)
        if matches.size == 0:
            raise ValueError(
                f"Layer {start_layer} is not present in cyclic layer sequence {self._layer_ids.tolist()}."  # noqa: E501
            )
        start = matches[0]
        indices = np.roll(np.arange(len(self._layer_ids)), -start)
        layer_order = self._layer_ids[indices]
        self._start_layer_number = layer_order[0]
        self._set_layer_order(layer_order, indices)
        self._basis_created = np.full_like(self.basis, np.nan)

    def stack_on(self, below_loc, below_height, below_layer=-1, below_state=None):
        """Generate this object's layer order and height from the object below.

        :param float below_loc:
            Absolute reference location of the object below in Angstrom.
        :param float below_height:
            Absolute top height of the object below in Angstrom.
        :param float below_layer:
            Top cyclic layer identifier of the object below.
        """
        if below_state is None:
            below_state = LayerState(self.layer_cycle, below_layer)
        if self.layer_behavior == "ignore":
            self.start_layer_number = -1
        else:
            start_layer = resolve_upper_start(
                below_state, self.layer_cycle, self.layer_transition
            )
            self._set_start_layer(start_layer)
        self.set_below(below_loc, below_height)

    def _stacking_metadata_to_str(self):
        if self.layer_transition is None:
            return ""
        pairs = ", ".join(
            "{} = {}".format(*pair) for pair in self.layer_transition.mapping.items()
        )
        return f"layer_transition: {pairs}\n"


def _parse_layer_transition(string):
    for line in string.splitlines():
        if line.strip().lower().startswith("layer_transition:"):
            values = line.split(":", 1)[1]
            mapping = {}
            for pair in values.split(","):
                lower, upper = pair.split("=")
                mapping[float(lower)] = float(upper)
            return mapping
    return None


def _parse_float_metadata(string, name, default):
    prefix = name.lower() + ":"
    equals_prefix = name.lower() + "="
    for line in string.splitlines():
        normalized = line.strip().lower().replace(" ", "")
        if normalized.startswith(prefix):
            return float(normalized.split(":", 1)[1])
        if normalized.startswith(equals_prefix):
            return float(normalized.split("=", 1)[1])
    return default


class EpitaxyInterface(_LayerStackingMixin, LinearFitFunctions):
    parameterOrder = "Width/cells Skew/cells"

    parameterLookup = {"W": 0, "S": 1}

    avail_types = ["skellam"]

    parameterLookup_inv = dict(map(reversed, parameterLookup.items()))

    def __init__(self, uc_top, uc_bottom, type="skellam", **kwargs):
        """
        sigma_calc : automatic selection of number of unitcells for calculation of interface
                     this will shift the location of the interface.
        fixed_ucs : create interface with fixed number of unitcells. Location of interface
                    will only change if skew != 0
                    total number of uc will be fixed_ucs*2 +1



        """  # noqa: E501
        super().__init__()
        profile = kwargs.pop("profile", None)
        if isinstance(type, SkellamProfile):
            profile = type
            type = "skellam"
        if profile is not None and not isinstance(profile, SkellamProfile):
            raise TypeError("EpitaxyInterface profile must be SkellamProfile")
        self.type = type
        self.profile = profile
        self.sigma_calc = kwargs.get("sigma_calc", 3)
        self.fixed_ucs = kwargs.get("fixed_ucs", False)
        self.set_ucs(uc_top, uc_bottom, **kwargs)
        if profile is None:
            self.basis = np.array([0.0, 0.0])
        else:
            self.basis = np.array([profile.width, profile.asymmetry])
        self._basis_created = np.array([np.nan, np.nan])
        self.basis_0 = np.array([0.0, 0.0])
        self.errors = None
        if "name" in kwargs:
            self.name = kwargs["name"]
        else:
            self.name = "unnamed"

        self.below_loc = 0.0
        self.below_H = 0.0
        self.below_layer = -1.0
        self._initialize_layer_stacking(kwargs)

    def set_ucs(self, uc_top, uc_bottom, **kwargs):
        if not np.all(uc_top.layers == uc_bottom.layers):
            raise ValueError(
                "Top and bottom layer numbers must be equal. "
                f"Is Top: {uc_top.layers}, Bottom: {uc_bottom.layers}"
            )
        self.uc_top = uc_top
        self.uc_bottom = uc_bottom
        self.reference_uc = self.uc_bottom

        self.uc_top.setReferenceUnitCell(
            self.uc_bottom, kwargs.get("rot", np.identity(3))
        )

        self.uc_layers_top = self.uc_top.split_in_layers()
        self.uc_layers_bottom = self.uc_bottom.split_in_layers()

        self._layer_ids = np.array(list(self.uc_layers_top))
        self._top_layers_base = [self.uc_layers_top[uc] for uc in self._layer_ids]
        self._bottom_layers_base = [self.uc_layers_bottom[uc] for uc in self._layer_ids]
        self._layerpos_base = np.array(
            [self.uc_top.layerpos[i] for i in self._layer_ids]
        )
        self.top_layers = [self.uc_layers_top[uc] for uc in self.uc_layers_top]
        self.bottom_layers = [self.uc_layers_bottom[uc] for uc in self.uc_layers_bottom]
        self.layerpos = np.copy(self._layerpos_base)

    @property
    def layers(self):
        """Return the cyclic layer identifiers."""
        return np.copy(self._layer_ids)

    @property
    def uc_area(self):
        """Return the lower interface unit-cell area in Angstrom squared.

        The lower unit cell defines the canonical lateral cell of the
        interface structure factor.
        """
        return self.uc_bottom.uc_area

    def _set_layer_order(self, layer_order, order):
        self.layer_order = np.asarray(layer_order)
        self._layer_order_indices = np.asarray(order)
        self.top_layers = [self._top_layers_base[i] for i in order]
        self.bottom_layers = [self._bottom_layers_base[i] for i in order]
        self.layerpos = _unwrapped_layer_positions(self._layerpos_base, order)

    def setReferenceUnitCell(self, uc, rotMatrix=np.identity(3)):
        """Set the reference frame on the interface and generated layers.

        :param UnitCell uc:
            Unit cell defining input reciprocal lattice units.
        :param numpy.ndarray rotMatrix:
            Optional 3-by-3 rotation from the reference frame into the
            interface crystal frame.
        """
        self.uc_top.setReferenceUnitCell(uc, rotMatrix)
        self.uc_bottom.setReferenceUnitCell(uc, rotMatrix)
        for layer in self._top_layers_base:
            layer.setReferenceUnitCell(uc, rotMatrix)
        for layer in self._bottom_layers_base:
            layer.setReferenceUnitCell(uc, rotMatrix)
        self.reference_uc = uc

    def setEnergy(self, E):
        """Set X-ray energy for the source and generated unit cells.

        :param float E:
            X-ray energy in eV.
        """
        self.E = E
        self.uc_top.setEnergy(E)
        self.uc_bottom.setEnergy(E)
        for layer in self._top_layers_base:
            layer.setEnergy(E)
        for layer in self._bottom_layers_base:
            layer.setEnergy(E)

    def set_below(self, loc, height):
        """Place the interface at the generated height of the object below."""
        self.below_loc = loc
        self.below_H = height
        self.createInterfaceCells()

    @property
    def loc_absolute(self):
        if np.any(self._basis_created != self.basis):
            self.createInterfaceCells()
        return self._loc_absolute

    @property
    def height_absolute(self):
        if np.any(self._basis_created != self.basis):
            self.createInterfaceCells()
        upper_layer = self.top_layers[-1]
        pos = upper_layer.coherentDomainMatrix[-1][2, 3]
        strain = upper_layer.coherentDomainMatrix[-1][2, 2]
        layer_id = self.layer_order[-1]
        layer_space = np.diff(self.layerpos, append=self.layerpos[0] + 1)
        H = (
            pos + strain * (self.uc_top.layerpos[layer_id] + layer_space[-1])
        ) * upper_layer.a[2]
        return H

    @property
    def pos_absolute(self):
        if np.any(self._basis_created != self.basis):
            self.createInterfaceCells()
        lower_layer = self.top_layers[0]
        matrix = lower_layer.coherentDomainMatrix[0]
        layer_id = self.layer_order[0]
        H = (
            matrix[2, 3] + matrix[2, 2] * self.uc_top.layerpos[layer_id]
        ) * lower_layer.a[2]
        return H

    @pos_absolute.setter
    def pos_absolute(self, pos):
        if np.any(self._basis_created != self.basis):
            self.createInterfaceCells()
        delta_height = pos - self.pos_absolute
        _move_domains(self.top_layers, delta_height)
        _move_domains(self.bottom_layers, delta_height)
        self.below_H = pos
        self._loc_absolute = self._loc_absolute_ref + self.below_H

    @property
    def end_layer_number(self):
        if np.any(self._basis_created != self.basis):
            self.createInterfaceCells()
        return self._end_layer_number

    @property
    def stacking_height_absolute(self):
        """Return the physical upper support height in Angstrom.

        Components stacked above this interface start at this height, so the
        interface alone owns every layer in its generated support.
        """
        return self.height_absolute

    @property
    def stacking_loc_absolute(self):
        """Return the nominal interface boundary location in Angstrom.

        A Film uses this location as the origin of its requested total width,
        while :attr:`stacking_height_absolute` supplies the physical support
        it must not duplicate.
        """
        return self.loc_absolute

    @property
    def layer_state(self):
        """Return the terminal layer state of the physical upper support."""
        return super().layer_state

    def createInterfaceCells(self):
        n_layers = len(self.uc_top.layers)
        if self.type == "skellam":
            tail_probability = (
                self.profile.tail_probability
                if self.profile is not None
                else DEFAULT_TAIL_PROBABILITY
            )
            profile = SkellamProfile(self.basis[0], self.basis[1], tail_probability)
            mu1, mu2 = profile.parameters(n_layers)

            loc = mu1 - mu2
            loc_int = int(round(loc / n_layers, 0)) * n_layers

            if self.fixed_ucs:
                uc_number = self.fixed_ucs
                unitcells = np.arange(-uc_number, uc_number) + loc_int
            elif self.profile is not None:
                support_low, support_high = profile.support(n_layers)
                first = (support_low // n_layers) * n_layers
                stop = (support_high // n_layers + 1) * n_layers
                unitcells = np.arange(first, stop)
                uc_number = -first
            else:
                sigma = self.basis[0] * n_layers
                uc_number = (
                    int(np.ceil((self.sigma_calc * sigma) / n_layers)) + 1
                ) * n_layers
                uc_number = int(uc_number)
                unitcells = np.arange(-uc_number, uc_number) + loc_int
            assert unitcells.size % len(self.top_layers) == 0

            probability_top = profile.occupancy(unitcells, n_layers).reshape(
                (-1, n_layers)
            )
            probability_bottom = 1.0 - probability_top
            sharp_top = (
                (unitcells >= loc).astype(np.float64).reshape((-1, n_layers))
            )
            # The upper material owns the complete generated support because
            # Film starts only above ``stacking_height_absolute``.  The lower
            # material remains a signed correction to the semi-infinite bulk
            # that already occupies the lower side of the nominal boundary.
            occupancy_top = probability_top
            occupancy_bottom = probability_bottom - (1.0 - sharp_top)

            a3_top = self.top_layers[0].a[2]
            a3_bottom = self.bottom_layers[0].a[2]

            for i, (uc_t, uc_b) in enumerate(zip(self.top_layers, self.bottom_layers)):
                uc_t.coherentDomainMatrix = []
                uc_t.coherentDomainOccupancy = np.ascontiguousarray(occupancy_top.T[i])
                uc_b.coherentDomainMatrix = []
                uc_b.coherentDomainOccupancy = np.ascontiguousarray(
                    occupancy_bottom.T[i]
                )

            mat_0 = np.vstack((np.identity(3).T, np.array([0, 0, 0]))).T

            ratio_top = a3_bottom / a3_top
            ratio_bottom = 1 / ratio_top
            h = 0.0

            for p_t, p_b in zip(probability_top, probability_bottom):
                top_strains = p_t + ratio_top * p_b
                bottom_strains = ratio_bottom * p_t + p_b
                physical_top_index = np.flatnonzero(
                    self._layer_order_indices == n_layers - 1
                )[0]
                cell_height = a3_top * top_strains[physical_top_index]

                for i, (uc_t, uc_b) in enumerate(
                    zip(self.top_layers, self.bottom_layers)
                ):
                    mat_top_i = np.copy(mat_0)
                    top_strain_and_h = top_strains[i]
                    mat_top_i[2, 2] = top_strain_and_h
                    relative_layer_position = self.layerpos[i] - self.layerpos[0]
                    top_layer_offset = (
                        relative_layer_position
                        - (self.uc_top.layerpos[self.layer_order[i]])
                    )
                    mat_top_i[2, 3] = h / a3_top + top_strain_and_h * top_layer_offset

                    uc_t.coherentDomainMatrix.append(mat_top_i)

                    mat_bottom_i = np.copy(mat_0)
                    bottom_strain_and_h = bottom_strains[i]
                    mat_bottom_i[2, 2] = bottom_strain_and_h
                    bottom_layer_offset = (
                        relative_layer_position
                        - (self.uc_bottom.layerpos[self.layer_order[i]])
                    )
                    mat_bottom_i[2, 3] = (
                        h / a3_bottom + bottom_strain_and_h * bottom_layer_offset
                    )
                    # h_bottom += bottom_strain_and_h
                    uc_b.coherentDomainMatrix.append(mat_bottom_i)
                h += cell_height
                # h_bottom += bottom_strain_and_h
                # h_top += top_strain_and_h

            loc_rescaled = loc - unitcells[0]
            uc_no_loc = int(np.floor(loc_rescaled)) // n_layers
            layer_no_loc = int(np.floor(loc_rescaled)) % n_layers
            loc_remainder = (loc_rescaled % n_layers) % 1
            loc_mat = self.top_layers[layer_no_loc].coherentDomainMatrix[uc_no_loc]
            self._loc_absolute_ref = (
                loc_mat[2, 3] * a3_top + loc_remainder * loc_mat[2, 2] * a3_top
            )
            translation = self.below_H - self._loc_absolute_ref
            self._loc_absolute = self.below_H
            _translate_domains(self.top_layers, translation)
            _translate_domains(self.bottom_layers, translation)
            self._end_layer_number = self.layer_order[-1]
            # self._loc_absolute = self._loc_absolute_ref + self.below_H

            self._basis_created = np.copy(self.basis)
        else:
            raise NotImplementedError(f"{self.type} is not a valid interface model")

    def F_uc(self, h, k, l):  # noqa: E741
        """Return the interface structure factor in electrons.

        Upper- and lower-material amplitudes are first converted to area
        densities using their own :attr:`UnitCell.uc_area`, added, and then
        multiplied by the lower unit-cell area. The returned amplitude
        therefore corresponds to one lower lateral unit cell.

        :param numpy.ndarray h:
            Reference-frame reciprocal coordinate in r.l.u.
        :param numpy.ndarray k:
            Reference-frame reciprocal coordinate in r.l.u.
        :param numpy.ndarray l:
            Reference-frame reciprocal coordinate in r.l.u.
        :returns:
            Complex interface amplitude in electrons.
        :rtype: numpy.ndarray
        """
        if np.any(self._basis_created != self.basis):
            self.createInterfaceCells()
        if ctr_accel_enabled():
            h, k, l = _ensure_contiguous(h, k, l, testOnly=False, astype=np.float64)  # noqa: E741
        F = np.zeros_like(l, dtype=np.complex128)
        for uc_t, uc_b in zip(self.top_layers, self.bottom_layers):
            F += self.uc_area * uc_t.F_uc(h, k, l) / uc_t.uc_area
            F += self.uc_area * uc_b.F_uc(h, k, l) / uc_b.uc_area
        return F

    def zDensity_G(self, z, h, k):
        if np.any(self._basis_created != self.basis):
            self.createInterfaceCells()
        rho = np.zeros_like(z, dtype=np.complex128)
        for uc_t, uc_b in zip(self.top_layers, self.bottom_layers):
            rho += uc_t.zDensity_G(z, h, k)
            rho += uc_b.zDensity_G(z, h, k)
        return rho

    def optical_profile(self):
        """Return the combined homogeneous optical profile of the interface.

        :returns:
            C-contiguous ``(N, 3)`` array with columns ``z``, ``delta``, and
            ``beta``.
        :rtype: numpy.ndarray
        """
        if np.any(self._basis_created != self.basis):
            self.createInterfaceCells()
        profiles = [uc.optical_profile() for uc in self.top_layers]
        profiles.extend(uc.optical_profile() for uc in self.bottom_layers)
        return combine_profiles(*profiles)

    def addFitParameter(self, indexarray, limits=(-np.inf, np.inf), **kwarg):
        """to assign multiple unitcells with the same fitparameter, provide list of
        unitcell names as kwarg `unitcell`
        """
        if len(np.array(indexarray).shape) < 2:
            return super().addFitParameter(indexarray, limits, **kwarg)
        if "unitcell" not in kwarg:
            raise ValueError(
                "Missing unit cell name. Provide unit cell name as kwarg 'unitcell'"
            )
        if isinstance(kwarg["unitcell"], list):
            fp = []
            for ucn in kwarg["unitcell"]:
                fp.append(self[ucn].addFitParameter(indexarray, limits, **kwarg))
            return fp
        else:
            return self[kwarg["unitcell"]].addFitParameter(indexarray, limits, **kwarg)

    def addRelParameter(self, indexarray, factors, limits=(-np.inf, np.inf), **kwarg):
        """to assign multiple unitcells with the same fitparameter, provide list of
        unitcell names as kwarg `unitcell`
        """
        if len(np.array(indexarray).shape) < 2:
            return super().addRelParameter(indexarray, factors, limits, **kwarg)
        if "unitcell" not in kwarg:
            raise ValueError(
                "Missing unit cell name. Provide unit cell name as kwarg 'unitcell'"
            )
        if isinstance(kwarg["unitcell"], list):
            fp = []
            for ucn in kwarg["unitcell"]:
                fp.append(
                    self[ucn].addRelParameter(indexarray, factors, limits, **kwarg)
                )
            return fp
        else:
            return self[kwarg["unitcell"]].addRelParameter(
                indexarray, factors, limits, **kwarg
            )

    def getStartParamAndLimits(self, force_recalculate=False):
        # if self.basis_0 is None:
        #    self.basis_0 = np.copy(self.basis)
        x0, lower, upper = super().getStartParamAndLimits(
            force_recalculate
        )  # absolute and relative
        top_x0, top_lower, top_upper = self.uc_top.getStartParamAndLimits(
            force_recalculate
        )
        bottom_x0, bottom_lower, bottom_upper = self.uc_bottom.getStartParamAndLimits(
            force_recalculate
        )
        return (
            np.concatenate([x0, top_x0, bottom_x0]),
            np.concatenate([lower, top_lower, bottom_lower]),
            np.concatenate([upper, top_upper, bottom_upper]),
        )

    def setFitParameters(self, x):
        abs_rel_no = len(self.parameters["absolute"]) + len(self.parameters["relative"])
        fp_top_no = len(self.uc_top.fitparnames)
        fp_bottom_no = len(self.uc_bottom.fitparnames)
        super().setFitParameters(x[:abs_rel_no])
        if self.profile is not None:
            self.profile = SkellamProfile(
                self.basis[0],
                self.basis[1],
                self.profile.tail_probability,
            )
        self.uc_top.setFitParameters(x[abs_rel_no : abs_rel_no + fp_top_no])
        self.uc_bottom.setFitParameters(
            x[abs_rel_no + fp_top_no : abs_rel_no + fp_top_no + fp_bottom_no]
        )

    def setLimits(self, lim):
        abs_rel_no = len(self.parameters["absolute"]) + len(self.parameters["relative"])
        fp_top_no = len(self.uc_top.fitparnames)
        fp_bottom_no = len(self.uc_bottom.fitparnames)
        super().setLimits(lim[:abs_rel_no])
        self.uc_top.setLimits(lim[abs_rel_no : abs_rel_no + fp_top_no])
        self.uc_bottom.setLimits(
            lim[abs_rel_no + fp_top_no : abs_rel_no + fp_top_no + fp_bottom_no]
        )

    def setFitErrors(self, errors):
        abs_rel_no = len(self.parameters["absolute"]) + len(self.parameters["relative"])
        fp_top_no = len(self.uc_top.fitparnames)
        fp_bottom_no = len(self.uc_bottom.fitparnames)
        super().setFitErrors(errors[:abs_rel_no])
        self.uc_top.setFitErrors(errors[abs_rel_no : abs_rel_no + fp_top_no])
        self.uc_bottom.setFitErrors(
            errors[abs_rel_no + fp_top_no : abs_rel_no + fp_top_no + fp_bottom_no]
        )

    def getFitErrors(self):
        err = super().getFitErrors()
        err_t = self.uc_top.getFitErrors()
        err_b = self.uc_bottom.getFitErrors()
        return np.concatenate([err, err_t, err_b])

    @property
    def fitparnames(self):
        return (
            super().fitparnames + self.uc_top.fitparnames + self.uc_bottom.fitparnames
        )

    @property
    def priors(self):
        return super().priors + self.uc_top.priors + self.uc_bottom.priors

    def parametersToDict(self):
        d = super().parametersToDict()
        d["unitcells"] = {}
        d["unitcells"]["top"] = self.uc_top.parametersToDict()
        d["unitcells"]["bottom"] = self.uc_bottom.parametersToDict()
        return d

    def clearParameters(self):
        super().clearParameters()
        self.uc_top.clearParameters()
        self.uc_bottom.clearParameters()

    def parametersFromDict(self, d, override_values=True):
        self.uc_top.parametersFromDict(d["unitcells"]["top"], override_values)
        self.uc_bottom.parametersFromDict(d["unitcells"]["bottom"], override_values)
        super().parametersFromDict(d, override_values)

    def updateFromParameters(self):
        """Update basis from the values stored in the Parameters"""
        self.uc_top.updateFromParameters()
        self.uc_bottom.updateFromParameters()
        super().updateFromParameters()
        if self.profile is not None:
            self.profile = SkellamProfile(
                self.basis[0],
                self.basis[1],
                self.profile.tail_probability,
            )

    def __getitem__(self, uc_name_or_index):
        if isinstance(uc_name_or_index, str):
            if uc_name_or_index.lower() in ["top", "t", "upper", self.uc_top.name]:
                return self.uc_top
            elif uc_name_or_index.lower() in [
                "bottom",
                "b",
                "lower",
                self.uc_bottom.name,
            ]:
                return self.uc_bottom
            else:
                raise KeyError(
                    f"No unit cell {uc_name_or_index} in EpitaxyInterface {self.name}"
                )
        else:
            raise ValueError(f"must be str, not {type(uc_name_or_index)}")

    def parameter_list(self):
        return (
            super().parameter_list()
            + self.uc_top.parameter_list()
            + self.uc_bottom.parameter_list()
        )

    @classmethod
    def fromStr(cls, string):
        with util.StringIO(string) as f:
            # parse header
            line = next_skip_comment(f).split()
            if line[0].lower() != "type":
                raise ValueError(
                    "You must specify a epitaxy type in line 1."
                    f" Available are {EpitaxyInterface.avail_types}"
                )
            if line[1].lower() not in EpitaxyInterface.avail_types:
                raise ValueError(
                    f"Expitaxy type {line[1]} is not valid."
                    f" Must be one of {EpitaxyInterface.avail_types}"
                )
            ep_type = line[1].lower()

            statistics = dict()
            line = next_skip_comment(f)
            while "Width" in line or "=" in line:  # parameter header or statistics line
                if "=" in line:
                    try:
                        splitted = [n.split(",") for n in line.split("=")]
                        splitted = [item for sublist in splitted for item in sublist]
                        for i in range(0, len(splitted), 2):
                            statistics[splitted[i].strip()] = float(splitted[i + 1])
                    except Exception:
                        print(f"Cannot read statistics string: {line}")
                line = next_skip_comment(f)
            # epitaxy parameters
            sline = line.split()
            if "+-" in sline:
                params = re.findall(r"\(([^)]+)", line)
                params_array = np.array(
                    [np.array(p.split("+-"), dtype=np.float64) for p in params]
                ).T
                basis = params_array[0]
                errors = params_array[1]
            else:
                basis = np.array(sline, dtype=np.float64)
                errors = None

        # very explicit searching for the lines containing TopUnitCell and BottomUnitCell:  # noqa: E501
        sp_str = string.splitlines()
        top_pos = -1
        bottom_pos = -1
        for i, l in enumerate(sp_str):  # noqa: E741
            if top_pos == -1:
                if "TopUnitCell" in l:
                    top_pos = i  # found it, and save line number
            if bottom_pos == -1:
                if "BottomUnitCell" in l:
                    bottom_pos = i  # found it, and save line number
            if top_pos != -1 and bottom_pos != -1:
                break
        else:
            msg = "Cannot create EpitaxyInterface. "
            if top_pos < 0:
                msg += "No TopUnitCell provided. "
            if bottom_pos < 0:
                msg += "No BottomUnitCell provided."
            raise ValueError(msg)

        classname_top, top_name = sp_str[top_pos].split(maxsplit=1)
        classname_bottom, bottom_name = sp_str[bottom_pos].split(maxsplit=1)

        assert classname_top == "TopUnitCell"
        assert classname_bottom == "BottomUnitCell"

        if top_pos < bottom_pos:
            top_uc_str = "\n".join(sp_str[top_pos + 1 : bottom_pos])
            bottom_uc_str = "\n".join(sp_str[bottom_pos + 1 :])
        else:
            top_uc_str = "\n".join(sp_str[top_pos + 1 :])
            bottom_uc_str = "\n".join(sp_str[bottom_pos + 1 : top_pos])

        uc_top = UnitCell.fromStr(top_uc_str)
        uc_bottom = UnitCell.fromStr(bottom_uc_str)

        uc_top.name = top_name
        uc_bottom.name = bottom_name

        tail_probability = _parse_float_metadata(
            string, "tail_probability", DEFAULT_TAIL_PROBABILITY
        )
        profile = SkellamProfile(basis[0], basis[1], tail_probability)
        epit = cls(
            uc_top,
            uc_bottom,
            ep_type,
            profile=profile,
            layer_transition=_parse_layer_transition(string),
        )
        epit.statistics = statistics
        epit.basis = basis
        epit.basis_0 = np.copy(basis)
        epit.errors = errors
        return epit

    def toStr(self, showErrors=True):
        """Serialize the interface as plain text.

        :param bool showErrors:
            Include propagated interface and nested unit-cell errors.
        :returns:
            Plain-text interface representation.
        :rtype: str
        """
        s = f"type {self.type}"
        s += (
            "\n"
            + EpitaxyInterface.parameterOrder
            + "\n"
            + self.epitToStr(showErrors=showErrors)
        )
        if (
            self.profile is not None
            and self.profile.tail_probability != DEFAULT_TAIL_PROBABILITY
        ):
            s += f"\ntail_probability: {self.profile.tail_probability:.12g}"
        metadata = self._stacking_metadata_to_str()
        if metadata:
            s += "\n" + metadata.rstrip()
        s += "\n\n"
        s += f"TopUnitCell {self.uc_top.name}\n"
        s += self.uc_top.toStr(showErrors=showErrors) + "\n\n"
        s += f"BottomUnitCell {self.uc_bottom.name}\n"
        s += self.uc_bottom.toStr(showErrors=showErrors) + "\n"
        return s

    def __repr__(self):
        return self.toStr()

    def epitToStr(self, showErrors=True):
        param = self.basis
        if (self.errors is not None) and showErrors:
            errors = self.errors
            l = []  # noqa: E741
            for p, err in zip(param, errors):
                l.append(f"({p:.5f} +- {err:.5f})")
            return "   ".join(l)
        else:
            l = []  # noqa: E741
            for p in param:
                l.append(f"{p:.5f} ")
            return "   ".join(l)


class Film(_LayerStackingMixin, LinearFitFunctions):
    parameterOrder = "Width/layers"

    parameterLookup = {"W": 0}

    parameterLookup_inv = dict(map(reversed, parameterLookup.items()))

    def __init__(self, unitcell, **kwargs):
        super().__init__()
        self.type = type
        self.set_ucs(unitcell, **kwargs)
        self.basis = np.array([0.0])
        self._basis_created = np.array([np.nan])
        self.basis_0 = np.array([0.0])
        self.errors = None
        if "name" in kwargs:
            self.name = kwargs["name"]
        else:
            self.name = "unnamed"

        self.below_loc = 0.0
        self.below_H = 0.0
        self.below_layer = -1.0
        self._initialize_layer_stacking(kwargs)

    def set_ucs(self, unitcell, **kwargs):
        self.unitcell = unitcell
        self.uc_layers = self.unitcell.split_in_layers()
        self._layer_ids = np.array(list(self.uc_layers))
        self._layer_ucs_base = [self.uc_layers[uc] for uc in self._layer_ids]
        self._layerpos_base = np.array(
            [self.unitcell.layerpos[i] for i in self._layer_ids]
        )
        self.layer_ucs = list(self._layer_ucs_base)
        self.layerpos = np.copy(self._layerpos_base)

    @property
    def layers(self):
        """Return the cyclic layer identifiers."""
        return np.copy(self._layer_ids)

    @property
    def uc_area(self):
        """Return the Film lateral unit-cell area in Angstrom squared."""
        return self.unitcell.uc_area

    def _set_layer_order(self, layer_order, order):
        self.layer_order = np.asarray(layer_order)
        self._layer_order_indices = np.asarray(order)
        self.layer_ucs = [self._layer_ucs_base[i] for i in order]
        self.layerpos = _unwrapped_layer_positions(self._layerpos_base, order)

    def setReferenceUnitCell(self, uc, rotMatrix=np.identity(3)):
        """Set the reference frame on the Film and every generated layer.

        :param UnitCell uc:
            Unit cell defining input reciprocal lattice units.
        :param numpy.ndarray rotMatrix:
            Optional 3-by-3 rotation from the reference frame into the Film
            crystal frame.
        """
        self.unitcell.setReferenceUnitCell(uc, rotMatrix)
        for layer in self._layer_ucs_base:
            layer.setReferenceUnitCell(uc, rotMatrix)

    def setEnergy(self, E):
        """Set X-ray energy for the Film and generated layers.

        :param float E:
            X-ray energy in eV.
        """
        self.E = E
        self.unitcell.setEnergy(E)
        for layer in self._layer_ucs_base:
            layer.setEnergy(E)

    @property
    def loc_absolute(self):
        return self.pos_absolute

    def set_below(self, loc, height):
        """Place Film above support while retaining its nominal width origin.

        :param float loc:
            Nominal lower boundary in Angstrom, used as the origin of
            :attr:`basis` width.
        :param float height:
            Physical upper support height in Angstrom. Generated unstrained
            Film layers begin here.
        """
        self.below_loc = loc
        self.below_H = height
        self.createLayers()

    @property
    def height_absolute(self):
        if np.any(self._basis_created != self.basis):
            self.createLayers()
        if self._layers_to_create == 0:
            return self.below_H
        upper_layer_id = self.end_layer_number
        upper_layer = self.uc_layers[upper_layer_id]
        idx = self.layer_ucs.index(upper_layer)
        pos = upper_layer.coherentDomainMatrix[-1][2, 3] * upper_layer.a[2]
        strain = upper_layer.coherentDomainMatrix[-1][2, 2]
        layerpos = self.unitcell.layerpos[upper_layer_id]
        layer_space = np.diff(self.layerpos, append=self.layerpos[0] + 1)
        H = pos + strain * (layerpos + layer_space[idx]) * upper_layer.a[2]
        return H

    @property
    def pos_absolute(self):
        if np.any(self._basis_created != self.basis):
            self.createLayers()
        if self._layers_to_create == 0:
            return self.below_H
        lower_layer = self.layer_ucs[0]
        matrix = lower_layer.coherentDomainMatrix[0]
        layer_id = self.layer_order[0]
        H = (
            matrix[2, 3] + matrix[2, 2] * self.unitcell.layerpos[layer_id]
        ) * lower_layer.a[2]
        return H

    @pos_absolute.setter
    def pos_absolute(self, pos):
        if np.any(self._basis_created != self.basis):
            self.createLayers()
        delta_height = pos - self.pos_absolute
        _move_domains(self.layer_ucs, delta_height)
        self.below_H = pos

    @property
    def end_layer_number(self):
        if np.any(self._basis_created != self.basis):
            self.createLayers()
        return self._end_layer_number

    def createLayers(self):
        """Create unstrained Film layers above lower-component support.

        The Film width in :attr:`basis` is measured from ``below_loc``. When
        the component below has already generated support up to ``below_H``,
        only the remaining width is represented by Film layers.
        """
        n_layers_in_uc = len(self.unitcell.layers)
        scaled_width = self.basis[0] - n_layers_in_uc * (
            (self.below_H - self.below_loc) / self.unitcell.a[2]
        )
        if scaled_width < 0 and not np.isclose(scaled_width, 0.0):
            raise ValueError(
                "Effective film width is shorter than lower-component support"
            )
        layers_to_create = int(round(scaled_width, 0))
        if layers_to_create < 0:
            raise ValueError(
                "Effective film width is shorter than lower-component support"
            )
        for i, uc in enumerate(self.layer_ucs):
            uc.coherentDomainMatrix = []
            uc.coherentDomainOccupancy = []

        self._layers_to_create = layers_to_create
        if layers_to_create == 0:
            self._end_layer_number = self.layer_order[-1]
            self._basis_created = np.copy(self.basis)
            return

        mat_0 = np.vstack((np.identity(3).T, np.array([0, 0, 0]))).T
        strain = self.unitcell.coherentDomainMatrix[0][2, 2]
        occup = self.unitcell.coherentDomainOccupancy[0]

        for layer_index in range(layers_to_create):
            order_index = layer_index % n_layers_in_uc
            cycle_index = layer_index // n_layers_in_uc
            uc = self.layer_ucs[order_index]
            mat_i = np.copy(mat_0)
            layer_id = self.layer_order[order_index]
            relative_layer_position = (
                cycle_index + self.layerpos[order_index] - self.layerpos[0]
            )
            layer_offset = relative_layer_position - self.unitcell.layerpos[layer_id]

            mat_i[2, 2] = strain
            mat_i[2, 3] = layer_offset * strain

            uc.coherentDomainMatrix.append(mat_i)
            uc.coherentDomainOccupancy.append(occup)

        upper_layer = uc.basis[0, 7]
        self._end_layer_number = upper_layer
        _translate_domains(self.layer_ucs, self.below_H)

        self._basis_created = np.copy(self.basis)

    def F_uc(self, h, k, l):  # noqa: E741
        """Return the Film structure factor in electrons.

        The result is the unnormalized sum over all generated Film layers and
        corresponds to one lateral Film unit cell.

        :param numpy.ndarray h:
            Reference-frame reciprocal coordinate in r.l.u.
        :param numpy.ndarray k:
            Reference-frame reciprocal coordinate in r.l.u.
        :param numpy.ndarray l:
            Reference-frame reciprocal coordinate in r.l.u.
        :returns:
            Complex Film amplitude in electrons.
        :rtype: numpy.ndarray
        """
        if np.any(self._basis_created != self.basis):
            self.createLayers()
        if ctr_accel_enabled():
            h, k, l = _ensure_contiguous(h, k, l, testOnly=False, astype=np.float64)  # noqa: E741
        F = np.zeros_like(l, dtype=np.complex128)
        for uc in self.layer_ucs:
            F += uc.F_uc(h, k, l)
        return F

    def zDensity_G(self, z, h, k):
        if np.any(self._basis_created != self.basis):
            self.createLayers()
        rho = np.zeros_like(z, dtype=np.complex128)
        for uc in self.layer_ucs:
            rho += uc.zDensity_G(z, h, k)
        return rho

    def optical_profile(self):
        """Return the combined homogeneous optical profile of the Film.

        :returns:
            C-contiguous ``(N, 3)`` array with columns ``z``, ``delta``, and
            ``beta``.
        :rtype: numpy.ndarray
        """
        if np.any(self._basis_created != self.basis):
            self.createLayers()
        return combine_profiles(*(uc.optical_profile() for uc in self.layer_ucs))

    def addFitParameter(self, indexarray, limits=(-np.inf, np.inf), **kwarg):
        """to assign multiple unitcells with the same fitparameter, provide list of
        unitcell names as kwarg `unitcell`
        """
        if len(np.array(indexarray).shape) < 2:
            return super().addFitParameter(indexarray, limits, **kwarg)

        return self.unitcell.addFitParameter(indexarray, limits, **kwarg)

    def addRelParameter(self, indexarray, factors, limits=(-np.inf, np.inf), **kwarg):
        """to assign multiple unitcells with the same fitparameter, provide list of
        unitcell names as kwarg `unitcell`
        """
        if len(np.array(indexarray).shape) < 2:
            return super().addRelParameter(indexarray, factors, limits, **kwarg)
        return self.unitcell.addRelParameter(indexarray, factors, limits, **kwarg)

    def getStartParamAndLimits(self, force_recalculate=False):
        # if self.basis_0 is None:
        #    self.basis_0 = np.copy(self.basis)
        x0, lower, upper = super().getStartParamAndLimits(
            force_recalculate
        )  # absolute and relative
        uc_x0, uc_lower, uc_upper = self.unitcell.getStartParamAndLimits(
            force_recalculate
        )
        return (
            np.concatenate([x0, uc_x0]),
            np.concatenate([lower, uc_lower]),
            np.concatenate([upper, uc_upper]),
        )

    def setFitParameters(self, x):
        abs_rel_no = len(self.parameters["absolute"]) + len(self.parameters["relative"])
        fp_no = len(self.unitcell.fitparnames)
        super().setFitParameters(x[:abs_rel_no])
        self.unitcell.setFitParameters(x[abs_rel_no : abs_rel_no + fp_no])

    def setLimits(self, lim):
        abs_rel_no = len(self.parameters["absolute"]) + len(self.parameters["relative"])
        fp_no = len(self.unitcell.fitparnames)
        super().setLimits(lim[:abs_rel_no])
        self.unitcell.setLimits(lim[abs_rel_no : abs_rel_no + fp_no])

    def setFitErrors(self, errors):
        abs_rel_no = len(self.parameters["absolute"]) + len(self.parameters["relative"])
        fp_no = len(self.unitcell.fitparnames)
        super().setFitErrors(errors[:abs_rel_no])
        self.unitcell.setFitErrors(errors[abs_rel_no : abs_rel_no + fp_no])

    def getFitErrors(self):
        err = super().getFitErrors()
        err_uc = self.unitcell.getFitErrors()
        return np.concatenate([err, err_uc])

    @property
    def fitparnames(self):
        return super().fitparnames + self.unitcell.fitparnames

    @property
    def priors(self):
        return super().priors + self.unitcell.priors

    def parametersToDict(self):
        d = super().parametersToDict()
        d["unitcells"] = {}
        d["unitcells"]["unitcell"] = self.unitcell.parametersToDict()
        return d

    def clearParameters(self):
        super().clearParameters()
        self.unitcell.clearParameters()

    def parametersFromDict(self, d, override_values=True):
        self.unitcell.parametersFromDict(d["unitcells"]["unitcell"], override_values)
        super().parametersFromDict(d, override_values)

    def updateFromParameters(self):
        """Update basis from the values stored in the Parameters"""
        self.unitcell.updateFromParameters()
        super().updateFromParameters()

    def __getitem__(self, uc_name_or_index):
        if isinstance(uc_name_or_index, str):
            if uc_name_or_index.lower() in ["uc", "unitcell", self.uc_top.name]:
                return self.unitcell
            else:
                raise KeyError(
                    f"No unit cell {uc_name_or_index} in EpitaxyInterface {self.name}"
                )
        else:
            raise ValueError(f"must be str, not {type(uc_name_or_index)}")

    def parameter_list(self):
        return super().parameter_list() + self.unitcell.parameter_list()

    @classmethod
    def fromStr(cls, string):
        with util.StringIO(string) as f:
            # parse header
            line = next_skip_comment(f).split()
            # if line[0].lower() != 'type':
            #    raise ValueError("You must specify a epitaxy type in line 1."
            #    " Available are %s" % EpitaxyInterface.avail_types)
            # if line[1].lower() not in EpitaxyInterface.avail_types:
            #    raise ValueError("Expitaxy type %s is not valid."
            #    " Must be one of %s" % (line[1], EpitaxyInterface.avail_types))
            # ep_type = line[1].lower()

            statistics = dict()
            line = next_skip_comment(f)
            while "Width" in line or "=" in line:  # parameter header or statistics line
                if "=" in line:
                    try:
                        splitted = [n.split(",") for n in line.split("=")]
                        splitted = [item for sublist in splitted for item in sublist]
                        for i in range(0, len(splitted), 2):
                            statistics[splitted[i].strip()] = float(splitted[i + 1])
                    except Exception:
                        print(f"Cannot read statistics string: {line}")
                line = next_skip_comment(f)
            # epitaxy parameters
            sline = line.split()
            if "+-" in sline:
                params = re.findall(r"\(([^)]+)", line)
                params_array = np.array(
                    [np.array(p.split("+-"), dtype=np.float64) for p in params]
                ).T
                basis = params_array[0]
                errors = params_array[1]
            else:
                basis = np.array(sline, dtype=np.float64)
                errors = None

        # very explicit searching for the lines containing TopUnitCell and BottomUnitCell:  # noqa: E501
        sp_str = string.splitlines()
        uc_pos = -1
        for i, l in enumerate(sp_str):  # noqa: E741
            if uc_pos == -1:
                if "UnitCell" in l:
                    uc_pos = i  # found it, and save line number
            if uc_pos != -1:
                break
        else:
            msg = "Cannot create Film. "
            if uc_pos < 0:
                msg += "No UnitCell provided. "
            raise ValueError(msg)

        uc_classname, uc_name = sp_str[uc_pos].split(maxsplit=1)

        assert uc_classname == "UnitCell"
        uc_str = "\n".join(sp_str[uc_pos + 1 :])

        uc = UnitCell.fromStr(uc_str)

        uc.name = uc_name

        film = cls(uc, layer_transition=_parse_layer_transition(string))
        film.statistics = statistics
        film.basis = basis
        film.basis_0 = np.copy(basis)
        film.errors = errors
        return film

    def toStr(self, showErrors=True):
        """Serialize the Film as plain text.

        :param bool showErrors:
            Include propagated Film and nested unit-cell errors.
        :returns:
            Plain-text Film representation.
        :rtype: str
        """
        # s = "type %s" % self.type
        s = "\n" + Film.parameterOrder + "\n" + self.filmToStr(showErrors=showErrors)
        metadata = self._stacking_metadata_to_str()
        if metadata:
            s += "\n" + metadata.rstrip()
        s += "\n\n"
        s += f"UnitCell {self.unitcell.name}\n"
        s += self.unitcell.toStr(showErrors=showErrors) + "\n"
        return s

    def __repr__(self):
        return self.toStr()

    def filmToStr(self, showErrors=True):
        param = self.basis
        if (self.errors is not None) and showErrors:
            errors = self.errors
            l = []  # noqa: E741
            for p, err in zip(param, errors):
                l.append(f"({p:.5f} +- {err:.5f})")
            return "   ".join(l)
        else:
            l = []  # noqa: E741
            for p in param:
                l.append(f"{p:.5f} ")
            return "   ".join(l)


class PoissonSurface(_LayerStackingMixin, LinearFitFunctions):
    """Signed Poisson growth or etching at a nominal Film boundary.

    For the current API, ``W`` is the signed Poisson mean in structural layers
    and ``offset`` is a deterministic height offset. Positive ``W`` models
    growth, negative ``W`` models dissolution or etching, and
    ``offset = -W`` preserves the original mean Film height.

    Legacy files using ``Width/layers deltaW/layers`` remain readable with
    their historical absolute-width semantics.
    """

    parameterOrder = "W/layers offset/layers"

    parameterLookup = {"W": 0, "offset": 1, "deltaW": 1}

    parameterLookup_inv = {0: "W", 1: "offset"}

    def __init__(self, unitcell, **kwargs):
        super().__init__()
        profile = kwargs.pop("profile", None)
        if profile is not None and not isinstance(profile, PoissonProfile):
            raise TypeError("PoissonSurface profile must be PoissonProfile")
        self.profile = profile
        self._legacy_absolute_width = kwargs.pop(
            "legacy_absolute_width", profile is None
        )
        self.type = type
        self.set_ucs(unitcell, **kwargs)
        if profile is None:
            self.basis = np.array([0.0, 0.0])
        else:
            self.basis = np.array([profile.mean_change, profile.offset])
        self._basis_created = np.array([np.nan, np.nan])
        self.basis_0 = np.array([0.0, 0.0])
        self.errors = None
        if "name" in kwargs:
            self.name = kwargs["name"]
        else:
            self.name = "unnamed"

        self.below_loc = 0.0
        self.below_H = 0.0
        self.below_layer = -1.0
        self._initialize_layer_stacking(kwargs)

    def set_ucs(self, unitcell, **kwargs):
        self.unitcell = unitcell
        self.uc_layers = self.unitcell.split_in_layers()
        self._layer_ids = np.array(list(self.uc_layers))
        self._layer_ucs_base = [self.uc_layers[uc] for uc in self._layer_ids]
        self._layerpos_base = np.array(
            [self.unitcell.layerpos[i] for i in self._layer_ids]
        )
        self.layer_ucs = list(self._layer_ucs_base)
        self.layerpos = np.copy(self._layerpos_base)

    @property
    def layers(self):
        """Return the cyclic layer identifiers."""
        return np.copy(self._layer_ids)

    @property
    def uc_area(self):
        """Return the surface lateral unit-cell area in Angstrom squared."""
        return self.unitcell.uc_area

    def _set_layer_order(self, layer_order, order):
        self.layer_order = np.asarray(layer_order)
        self._layer_order_indices = np.asarray(order)
        self.layer_ucs = [self._layer_ucs_base[i] for i in order]
        self.layerpos = _unwrapped_layer_positions(self._layerpos_base, order)

    def setReferenceUnitCell(self, uc, rotMatrix=np.identity(3)):
        """Set the reference frame on the surface and generated layers.

        :param UnitCell uc:
            Unit cell defining input reciprocal lattice units.
        :param numpy.ndarray rotMatrix:
            Optional 3-by-3 rotation from the reference frame into the surface
            crystal frame.
        """
        self.unitcell.setReferenceUnitCell(uc, rotMatrix)
        for layer in self._layer_ucs_base:
            layer.setReferenceUnitCell(uc, rotMatrix)

    def setEnergy(self, E):
        """Set X-ray energy for the surface and generated layers.

        :param float E:
            X-ray energy in eV.
        """
        self.E = E
        self.unitcell.setEnergy(E)
        for layer in self._layer_ucs_base:
            layer.setEnergy(E)

    @property
    def loc_absolute(self):
        return self.pos_absolute

    def set_below(self, loc, height):
        self.below_loc = loc
        self.below_H = height
        self.createLayers()

    @property
    def height_absolute(self):
        if np.any(self._basis_created != self.basis):
            self.createLayers()
        upper_layer_id = self.end_layer_number
        upper_layer = self.uc_layers[upper_layer_id]
        idx = self.layer_ucs.index(upper_layer)
        pos = upper_layer.coherentDomainMatrix[-1][2, 3] * upper_layer.a[2]
        strain = upper_layer.coherentDomainMatrix[-1][2, 2]
        layerpos = self.unitcell.layerpos[upper_layer_id]
        layer_space = np.diff(self.layerpos, append=self.layerpos[0] + 1)
        H = pos + strain * (layerpos + layer_space[idx]) * upper_layer.a[2]
        return H

    @property
    def pos_absolute(self):
        if np.any(self._basis_created != self.basis):
            self.createLayers()
        lower_layer = self.layer_ucs[0]
        matrix = lower_layer.coherentDomainMatrix[0]
        layer_id = self.layer_order[0]
        H = (
            matrix[2, 3] + matrix[2, 2] * self.unitcell.layerpos[layer_id]
        ) * lower_layer.a[2]
        return H

    @pos_absolute.setter
    def pos_absolute(self, pos):
        if np.any(self._basis_created != self.basis):
            self.createLayers()
        delta_height = pos - self.pos_absolute
        _move_domains(self.layer_ucs, delta_height)
        self.below_H = pos

    @property
    def end_layer_number(self):
        if np.any(self._basis_created != self.basis):
            self.createLayers()
        return self._end_layer_number

    @property
    def stacking_height_absolute(self):
        """Return the expected surface height in Angstrom."""
        if self._legacy_absolute_width:
            return self.height_absolute
        return self.mean_height_absolute

    @property
    def stacking_loc_absolute(self):
        """Return the nominal boundary location in Angstrom."""
        return self.stacking_height_absolute

    @property
    def mean_height_absolute(self):
        """Return the expected surface height in Angstrom."""
        layer_height = self.unitcell.a[2] / len(self.unitcell.layers)
        return self.below_H + layer_height * (self.basis[0] + self.basis[1])

    def createLayers(self):
        """Create layer domains and their Poisson surface occupancies."""
        n_layers_in_uc = len(self.unitcell.layers)
        tail_probability = (
            self.profile.tail_probability
            if self.profile is not None
            else DEFAULT_TAIL_PROBABILITY
        )
        if self._legacy_absolute_width:
            delta_width = self.basis[1]
            if delta_width < 0:
                raise ValueError("legacy deltaW must be greater than or equal to zero")
            profile = PoissonProfile(
                mean_change=delta_width,
                tail_probability=tail_probability,
            )
            mean_width = self.basis[0] - n_layers_in_uc * (
                (self.below_H - self.below_loc) / self.unitcell.a[2]
            )
            scaled_width = mean_width + delta_width
            layers_to_create = int(round(scaled_width, 0))
            if layers_to_create <= 0:
                raise ValueError("Effective film width <= 0. Cannot create layers")
            layer_numbers = np.arange(layers_to_create)
            if delta_width == 0:
                layer_occupancy = np.ones(layers_to_create)
            else:
                minimum_width = mean_width - delta_width
                layer_occupancy = profile.occupancy(layer_numbers - minimum_width)
        else:
            profile = PoissonProfile(
                mean_change=self.basis[0],
                offset=self.basis[1],
                tail_probability=tail_probability,
            )
            support_low, support_high = profile.support()
            layer_numbers = np.arange(support_low, support_high + 1)
            layer_occupancy = profile.correction(layer_numbers)
            represented = np.abs(layer_occupancy) > tail_probability
            layer_numbers = layer_numbers[represented]
            layer_occupancy = layer_occupancy[represented]
            layers_to_create = len(layer_numbers)
        for i, uc in enumerate(self.layer_ucs):
            uc.coherentDomainMatrix = []
            uc.coherentDomainOccupancy = []

        mat_0 = np.vstack((np.identity(3).T, np.array([0, 0, 0]))).T
        strain = self.unitcell.coherentDomainMatrix[0][2, 2]
        occup = self.unitcell.coherentDomainOccupancy[0]

        layer_occupancy *= occup

        for layer_index, layer_number in enumerate(layer_numbers):
            order_index = layer_number % n_layers_in_uc
            cycle_index = layer_number // n_layers_in_uc
            uc = self.layer_ucs[order_index]
            mat_i = np.copy(mat_0)
            layer_id = self.layer_order[order_index]
            relative_layer_position = (
                cycle_index + self.layerpos[order_index] - self.layerpos[0]
            )
            layer_offset = relative_layer_position - self.unitcell.layerpos[layer_id]

            mat_i[2, 2] = strain
            mat_i[2, 3] = layer_offset * strain

            uc.coherentDomainMatrix.append(mat_i)
            uc.coherentDomainOccupancy.append(layer_occupancy[layer_index])

        if layers_to_create:
            upper_layer = uc.basis[0, 7]
            self._end_layer_number = upper_layer
        else:
            self._end_layer_number = self.start_layer_number
        _translate_domains(self.layer_ucs, self.below_H)

        self._basis_created = np.copy(self.basis)

    def F_uc(self, h, k, l):  # noqa: E741
        """Return the Poisson surface correction in electrons.

        Positive occupancies add grown material and negative occupancies remove
        etched material relative to the sharp Film boundary. The amplitude
        corresponds to one lateral surface unit cell.

        :param numpy.ndarray h:
            Reference-frame reciprocal coordinate in r.l.u.
        :param numpy.ndarray k:
            Reference-frame reciprocal coordinate in r.l.u.
        :param numpy.ndarray l:
            Reference-frame reciprocal coordinate in r.l.u.
        :returns:
            Complex surface-correction amplitude in electrons.
        :rtype: numpy.ndarray
        """
        if np.any(self._basis_created != self.basis):
            self.createLayers()
        if ctr_accel_enabled():
            h, k, l = _ensure_contiguous(h, k, l, testOnly=False, astype=np.float64)  # noqa: E741
        F = np.zeros_like(l, dtype=np.complex128)
        for uc in self.layer_ucs:
            if uc.coherentDomainMatrix:
                F += uc.F_uc(h, k, l)
        return F

    def zDensity_G(self, z, h, k):
        if np.any(self._basis_created != self.basis):
            self.createLayers()
        rho = np.zeros_like(z, dtype=np.complex128)
        for uc in self.layer_ucs:
            if uc.coherentDomainMatrix:
                rho += uc.zDensity_G(z, h, k)
        return rho

    def optical_profile(self):
        """Return the homogeneous optical profile of the Poisson surface.

        :returns:
            C-contiguous ``(N, 3)`` array with columns ``z``, ``delta``, and
            ``beta``.
        :rtype: numpy.ndarray
        """
        if np.any(self._basis_created != self.basis):
            self.createLayers()
        return combine_profiles(
            *(uc.optical_profile() for uc in self.layer_ucs if uc.coherentDomainMatrix)
        )

    def addFitParameter(self, indexarray, limits=(-np.inf, np.inf), **kwarg):
        """to assign multiple unitcells with the same fitparameter, provide list of
        unitcell names as kwarg `unitcell`
        """
        if len(np.array(indexarray).shape) < 2:
            return super().addFitParameter(indexarray, limits, **kwarg)

        return self.unitcell.addFitParameter(indexarray, limits, **kwarg)

    def addRelParameter(self, indexarray, factors, limits=(-np.inf, np.inf), **kwarg):
        """to assign multiple unitcells with the same fitparameter, provide list of
        unitcell names as kwarg `unitcell`
        """
        if len(np.array(indexarray).shape) < 2:
            return super().addRelParameter(indexarray, factors, limits, **kwarg)
        return self.unitcell.addRelParameter(indexarray, factors, limits, **kwarg)

    def getStartParamAndLimits(self, force_recalculate=False):
        # if self.basis_0 is None:
        #    self.basis_0 = np.copy(self.basis)
        x0, lower, upper = super().getStartParamAndLimits(
            force_recalculate
        )  # absolute and relative
        uc_x0, uc_lower, uc_upper = self.unitcell.getStartParamAndLimits(
            force_recalculate
        )
        return (
            np.concatenate([x0, uc_x0]),
            np.concatenate([lower, uc_lower]),
            np.concatenate([upper, uc_upper]),
        )

    def setFitParameters(self, x):
        abs_rel_no = len(self.parameters["absolute"]) + len(self.parameters["relative"])
        fp_no = len(self.unitcell.fitparnames)
        super().setFitParameters(x[:abs_rel_no])
        if self.profile is not None:
            self.profile = PoissonProfile(
                mean_change=self.basis[0],
                offset=self.basis[1],
                tail_probability=self.profile.tail_probability,
            )
        self.unitcell.setFitParameters(x[abs_rel_no : abs_rel_no + fp_no])

    def setLimits(self, lim):
        abs_rel_no = len(self.parameters["absolute"]) + len(self.parameters["relative"])
        fp_no = len(self.unitcell.fitparnames)
        super().setLimits(lim[:abs_rel_no])
        self.unitcell.setLimits(lim[abs_rel_no : abs_rel_no + fp_no])

    def setFitErrors(self, errors):
        abs_rel_no = len(self.parameters["absolute"]) + len(self.parameters["relative"])
        fp_no = len(self.unitcell.fitparnames)
        super().setFitErrors(errors[:abs_rel_no])
        self.unitcell.setFitErrors(errors[abs_rel_no : abs_rel_no + fp_no])

    def getFitErrors(self):
        err = super().getFitErrors()
        err_uc = self.unitcell.getFitErrors()
        return np.concatenate([err, err_uc])

    @property
    def fitparnames(self):
        return super().fitparnames + self.unitcell.fitparnames

    @property
    def priors(self):
        return super().priors + self.unitcell.priors

    def parametersToDict(self):
        d = super().parametersToDict()
        d["unitcells"] = {}
        d["unitcells"]["unitcell"] = self.unitcell.parametersToDict()
        return d

    def clearParameters(self):
        super().clearParameters()
        self.unitcell.clearParameters()

    def parametersFromDict(self, d, override_values=True):
        self.unitcell.parametersFromDict(d["unitcells"]["unitcell"], override_values)
        super().parametersFromDict(d, override_values)

    def updateFromParameters(self):
        """Update basis from the values stored in the Parameters"""
        self.unitcell.updateFromParameters()
        super().updateFromParameters()
        if self.profile is not None:
            self.profile = PoissonProfile(
                mean_change=self.basis[0],
                offset=self.basis[1],
                tail_probability=self.profile.tail_probability,
            )

    def __getitem__(self, uc_name_or_index):
        if isinstance(uc_name_or_index, str):
            if uc_name_or_index.lower() in [
                "uc",
                "unitcell",
                self.unitcell.name.lower(),
            ]:
                return self.unitcell
            else:
                raise KeyError(
                    f"No unit cell {uc_name_or_index} in PoissonSurface {self.name}"
                )
        else:
            raise ValueError(f"must be str, not {type(uc_name_or_index)}")

    def parameter_list(self):
        return super().parameter_list() + self.unitcell.parameter_list()

    @classmethod
    def fromStr(cls, string):
        with util.StringIO(string) as f:
            # parse header
            line = next_skip_comment(f).split()
            # if line[0].lower() != 'type':
            #    raise ValueError("You must specify a epitaxy type in line 1."
            #    " Available are %s" % EpitaxyInterface.avail_types)
            # if line[1].lower() not in EpitaxyInterface.avail_types:
            #    raise ValueError("Expitaxy type %s is not valid."
            #    " Must be one of %s" % (line[1], EpitaxyInterface.avail_types))
            # ep_type = line[1].lower()

            statistics = dict()
            line = next_skip_comment(f)
            while "Width" in line or "=" in line:  # parameter header or statistics line
                if "=" in line:
                    try:
                        splitted = [n.split(",") for n in line.split("=")]
                        splitted = [item for sublist in splitted for item in sublist]
                        for i in range(0, len(splitted), 2):
                            statistics[splitted[i].strip()] = float(splitted[i + 1])
                    except Exception:
                        print(f"Cannot read statistics string: {line}")
                line = next_skip_comment(f)
            # epitaxy parameters
            sline = line.split()
            if "+-" in sline:
                params = re.findall(r"\(([^)]+)", line)
                params_array = np.array(
                    [np.array(p.split("+-"), dtype=np.float64) for p in params]
                ).T
                basis = params_array[0]
                errors = params_array[1]
            else:
                basis = np.array(sline, dtype=np.float64)
                errors = None

        # very explicit searching for the lines containing TopUnitCell and BottomUnitCell:  # noqa: E501
        sp_str = string.splitlines()
        uc_pos = -1
        for i, l in enumerate(sp_str):  # noqa: E741
            if uc_pos == -1:
                if "UnitCell" in l:
                    uc_pos = i  # found it, and save line number
            if uc_pos != -1:
                break
        else:
            msg = "Cannot create Film. "
            if uc_pos < 0:
                msg += "No UnitCell provided. "
            raise ValueError(msg)

        uc_classname, uc_name = sp_str[uc_pos].split(maxsplit=1)

        assert uc_classname == "UnitCell"
        uc_str = "\n".join(sp_str[uc_pos + 1 :])

        uc = UnitCell.fromStr(uc_str)

        uc.name = uc_name

        tail_probability = _parse_float_metadata(
            string, "tail_probability", DEFAULT_TAIL_PROBABILITY
        )
        legacy_parameter_names = any("deltaW" in line for line in string.splitlines())
        legacy_absolute_width = legacy_parameter_names
        if legacy_absolute_width:
            profile = PoissonProfile(
                mean_change=basis[1],
                tail_probability=tail_probability,
            )
        elif legacy_parameter_names:
            # Migrate the short-lived nominal-cursor format where W was zero
            # and deltaW was a positive, mean-preserving roughness parameter.
            profile = PoissonProfile(
                mean_change=basis[1],
                offset=-basis[1],
                tail_probability=tail_probability,
            )
            basis = np.array([profile.mean_change, profile.offset])
        else:
            profile = PoissonProfile(
                mean_change=basis[0],
                offset=basis[1],
                tail_probability=tail_probability,
            )
        film = cls(
            uc,
            profile=profile,
            layer_transition=_parse_layer_transition(string),
            legacy_absolute_width=legacy_absolute_width,
        )
        film.statistics = statistics
        film.basis = basis
        film.basis_0 = np.copy(basis)
        film.errors = errors
        return film

    def toStr(self, showErrors=True):
        """Serialize the Poisson surface as plain text.

        :param bool showErrors:
            Include propagated surface and nested unit-cell errors.
        :returns:
            Plain-text Poisson-surface representation.
        :rtype: str
        """
        # s = "type %s" % self.type
        if self._legacy_absolute_width:
            parameter_order = "Width/layers deltaW/layers"
        else:
            parameter_order = PoissonSurface.parameterOrder
        s = "\n" + parameter_order + "\n" + self.filmToStr(showErrors=showErrors)
        if (
            self.profile is not None
            and self.profile.tail_probability != DEFAULT_TAIL_PROBABILITY
        ):
            s += f"\ntail_probability = {self.profile.tail_probability:.12g}"
        metadata = self._stacking_metadata_to_str()
        if metadata:
            s += "\n" + metadata.rstrip()
        s += "\n\n"
        s += f"UnitCell {self.unitcell.name}\n"
        s += self.unitcell.toStr(showErrors=showErrors) + "\n"
        return s

    def __repr__(self):
        return self.toStr()

    def filmToStr(self, showErrors=True):
        param = self.basis
        if (self.errors is not None) and showErrors:
            errors = self.errors
            l = []  # noqa: E741
            for p, err in zip(param, errors):
                l.append(f"({p:.5f} +- {err:.5f})")
            return "   ".join(l)
        else:
            l = []  # noqa: E741
            for p in param:
                l.append(f"{p:.5f} ")
            return "   ".join(l)
