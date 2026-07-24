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

import numpy as np
from .. import util
import warnings
import os

# random.seed(45)
import errno

from .CTRutil import ParameterType, Parameter, next_skip_comment

from .CTRuc import WaterModel, UnitCell
from .CTRfilm import EpitaxyInterface, Film, PoissonSurface
from .CTRoptics import (
    add_structural_to_sampled_profile,
    combine_profiles,
    simplify_profile,
    solve_wavefield,
    stratify_profile,
    top_layer_spacing,
)
from .CTRdistributions import PoissonProfile, SkellamProfile  # noqa: F401
from .CTRstacking import (  # noqa: F401
    LayerCycle,
    LayerState,
    LayerTransition,
)


class SXRDCrystal:
    """Compose bulk and surface amplitudes on one reference lateral cell.

    Every component ``F_uc`` method returns an unnormalized structure factor
    in electrons for that component's own lateral unit cell. This class
    converts each amplitude to the selected reference area using
    ``reference_uc.uc_area / component.uc_area`` before adding it.

    The bulk unit cell is the default reciprocal-coordinate and area reference.
    No illuminated sample area, detector response, or experimental scale
    factor is included.
    """

    def __init__(self, uc_bulk, *uc_surface, **keyargs):
        """Create a surface X-ray diffraction crystal model.

        :param UnitCell uc_bulk:
            Semi-infinite bulk unit cell. It is also the default reciprocal
            coordinate and lateral-area reference.
        :param uc_surface:
            Ordered surface, interface, Film, or water components.
        :param UnitCell reference_uc:
            Optional reference unit cell. Defaults to ``uc_bulk``.
        :param numpy.ndarray reference_rotation:
            Optional 3-by-3 rotation from the reference crystal frame into each
            component's frame. Defaults to the identity matrix.
        :param numpy.ndarray stacking:
            Integer stacking levels for the surface components.
        :param float atten:
            Dimensionless bulk attenuation exponent per reference-cell
            out-of-plane repeat.
        """
        self.uc_bulk = uc_bulk
        self.uc_surface_list = list(uc_surface)
        self.enable_uc_stacking = keyargs.get("enable_stacking", False)
        if not self.enable_uc_stacking:
            for uc in self.uc_surface_list:
                stackable = Film | PoissonSurface | EpitaxyInterface | WaterModel
                if isinstance(uc, stackable):
                    self.enable_uc_stacking = True
                    break

        self.uc_stacking = keyargs.get(
            "stacking", np.arange(len(self.uc_surface_list))[::-1]
        )
        order = np.argsort(self.uc_stacking)
        self.uc_surface_list_ordered = np.array(self.uc_surface_list)[order]
        self.uc_stacking_ordered = self.uc_stacking[order]
        self.domains = [
            [(np.identity(3, dtype=np.float64), 1)] for i in self.uc_surface_list
        ]
        self.atten = keyargs.get("atten", 0.01)
        self.specularRes = keyargs.get("spec_res", 0.0)
        self.name = keyargs.get("name", "xtal")

        if not self.uc_surface_list:
            self.weights = np.array([], dtype=np.float64)
        elif self.enable_uc_stacking:
            # Stacked objects are additive material slices, not alternative
            # surface domains. Their complete occupancies must therefore be
            # included unless an explicit crystal weight is later assigned.
            self.weights = np.ones(len(self.uc_surface_list), dtype=np.float64)
        else:
            self.weights = np.ones(len(self.uc_surface_list), dtype=np.float64) / len(
                self.uc_surface_list
            )
        self.weights_0 = np.copy(self.weights)
        self._weights_parvalues = None
        self.werrors = None
        self._werrors_parvalues = None
        self.parameters = {"coupled": [], "weight": [], "domain": []}
        self._parIdNo = 0
        self.fit_metadata_cache = None
        reference_uc = keyargs.get("reference_uc", self.uc_bulk)
        reference_rotation = keyargs.get("reference_rotation", np.identity(3))
        self.setGlobalReferenceUnitCell(reference_uc, reference_rotation)

    def parametersToDict(self):
        d = dict()
        d["weights_0"] = self.weights_0
        d["_parIdNo"] = self._parIdNo

        for p in self.parameters:
            d[p] = dict()
            for i, param in enumerate(self.parameters[p]):
                d[p][f"{i}_{param.name}"] = param.asdict()

        d["unitcells"] = dict()
        for uc in self.uc_surface_list:
            d["unitcells"][uc.name] = uc.parametersToDict()
        d["unitcells"]["bulk"] = self.uc_bulk.parametersToDict()
        return d

    def clearParameters(self):
        for p in self.parameters:
            self.parameters[p] = []
        self._parIdNo = 0
        self.fit_metadata_cache = None
        self._weights_parvalues = None
        self._werrors_parvalues = None
        for uc in self.uc_surface_list:
            uc.clearParameters()
        self.uc_bulk.clearParameters()

    def parametersFromDict(self, d, override_values=True):
        for p in self.parameters:
            self.parameters[p] = []

        self.weights_0 = d["weights_0"].astype(np.float64)
        self._parIdNo = int(d["_parIdNo"])

        for p in self.parameters:
            if p in d:
                for dkey in sorted(d[p].keys()):
                    self.parameters[p].append(Parameter(**d[p][dkey]))

        for uc in d["unitcells"]:
            self[uc].parametersFromDict(d["unitcells"][uc], override_values)

        if override_values:
            try:
                self.updateFromParameters()
            except Exception as e:
                warnings.warn(
                    f"Can not apply weight fit parameter values to this SXRDCrystal: {e}"  # noqa: E501
                )
        else:
            warnings.warn(
                "Weight values not updated from Parameter values. This can cause a mismatch between weight and fitparameter values!"  # noqa: E501
            )

    def updateFromParameters(self):
        self.weights = np.copy(self.weights_0)
        self.werrors = np.full_like(self.weights, np.nan)
        no_errors = False
        for par in self.parameters["weight"]:
            if par.value is not None:
                self.weights[par.indices] += par.factors * par.value
                self.werrors[par.indices] = np.nan_to_num(
                    self.werrors[par.indices], nan=0.0
                )
            else:
                raise ValueError(
                    f"Can not set weight values from parameters. Value of Parameter {par.name} is None."  # noqa: E501
                )
            if par.error is not None:
                self.werrors[par.indices] += par.factors * par.error
            else:
                no_errors = True
        self._weights_parvalues = np.copy(self.weights)
        if not no_errors:
            self._werrors_parvalues = np.copy(self.werrors)

    def setGlobalReferenceUnitCell(self, uc_reference, rotMatrix=np.identity(3)):
        """Set the crystal-wide coordinate and lateral-area reference.

        Reciprocal coordinates supplied to structure-factor methods are
        interpreted in ``uc_reference`` reciprocal lattice units. The same
        unit cell supplies the reference area used by :meth:`F`.

        This method propagates the coordinate transform to bulk, source unit
        cells, and all generated Film, interface, and surface layer cells.

        :param UnitCell uc_reference:
            Unit cell defining reciprocal lattice units and ``uc_area``.
        :param numpy.ndarray rotMatrix:
            Optional 3-by-3 rotation from the reference crystal frame into the
            component crystal frames.
        """
        self.reference_uc = uc_reference
        self.reference_area = uc_reference.uc_area
        self.reference_rotation = np.asarray(rotMatrix, dtype=np.float64)
        self.uc_bulk.setReferenceUnitCell(uc_reference, rotMatrix)
        for uc in self.uc_surface_list:
            uc.setReferenceUnitCell(uc_reference, rotMatrix)

    def F_surf(self, harray, karray, Larray):
        """Return the combined surface amplitude in electrons.

        Each component amplitude is converted from its own lateral unit cell to
        the configured reference lateral cell. The bulk contribution is not
        included.

        :param numpy.ndarray harray:
            Reference-frame reciprocal coordinate in r.l.u.
        :param numpy.ndarray karray:
            Reference-frame reciprocal coordinate in r.l.u.
        :param numpy.ndarray Larray:
            Reference-frame reciprocal coordinate in r.l.u.
        :returns:
            Complex surface amplitude in electrons per reference lateral cell.
        :rtype: numpy.ndarray
        """
        F = np.zeros_like(Larray, dtype=np.complex128)
        for uc, weight in zip(self.uc_surface_list, self.weights):
            area_scale = self.reference_area / uc.uc_area
            F += area_scale * weight * uc.F_uc(harray, karray, Larray)
        return F

    def apply_stacking(self):
        """Generate object positions and cyclic layers from bottom to top."""
        self.enable_uc_stacking = True
        height_new = 0.0
        i = 0
        layer_state_new = getattr(self.uc_bulk, "layer_state", None)
        layer_number_new = layer_state_new.layer if layer_state_new is not None else -1
        loc_new = 0.0
        component_new = self.uc_bulk
        for stacking_level in np.unique(self.uc_stacking_ordered):
            below_height = height_new
            below_layer = layer_number_new
            below_loc = loc_new
            below_component = component_new
            while self.uc_stacking_ordered[i] == stacking_level:
                uc = self.uc_surface_list_ordered[i]
                if hasattr(uc, "stack_on"):
                    stack_kwargs = {"below_state": layer_state_new}
                    if isinstance(uc, PoissonSurface):
                        stack_kwargs["below_component"] = below_component
                    uc.stack_on(
                        below_loc,
                        below_height,
                        below_layer,
                        **stack_kwargs,
                    )
                else:
                    uc.start_layer_number = below_layer
                    uc.pos_absolute = below_height
                if hasattr(uc, "stacking_height_absolute"):
                    height_new = uc.stacking_height_absolute
                else:
                    height_new = uc.height_absolute
                if hasattr(uc, "stacking_loc_absolute"):
                    loc_new = uc.stacking_loc_absolute
                else:
                    loc_new = uc.loc_absolute
                layer_number_new = uc.end_layer_number
                layer_state_new = getattr(uc, "layer_state", None)
                component_new = uc
                i += 1
                if i == len(self.uc_stacking_ordered):
                    return

    def F(self, harray, karray, Larray):
        """Return the complete crystal structure factor in electrons.

        The returned amplitude corresponds to one lateral unit cell of
        :attr:`reference_uc`. For every component ``j``, the raw amplitude in
        electrons is scaled by

        ``reference_uc.uc_area / component.uc_area``.

        Component weights, coherent-domain occupancies, and attenuation are
        dimensionless. No illuminated footprint or experimental intensity
        scale is included; calculated intensity is proportional to
        ``abs(F)**2``.

        :param numpy.ndarray harray:
            Reference-frame reciprocal coordinate in r.l.u.
        :param numpy.ndarray karray:
            Reference-frame reciprocal coordinate in r.l.u.
        :param numpy.ndarray Larray:
            Reference-frame reciprocal coordinate in r.l.u.
        :returns:
            Complex crystal amplitude in electrons per reference lateral cell.
        :rtype: numpy.ndarray
        """
        bulk_scale = self.reference_area / self.uc_bulk.uc_area
        F = bulk_scale * self.uc_bulk.F_bulk(harray, karray, Larray, self.atten)
        hkl = np.vstack((harray, karray, Larray))
        if self.enable_uc_stacking:
            self.apply_stacking()

        for uc, weight, domains in zip(
            self.uc_surface_list, self.weights, self.domains
        ):
            area_scale = self.reference_area / uc.uc_area
            for matrix, occup in domains:
                hkl_n = np.dot(matrix, hkl)
                F += area_scale * occup * weight * uc.F_uc(hkl_n[0], hkl_n[1], hkl_n[2])
        return F

    def setDomain(self, uc_no, domains):
        """
        should be a list of tuples with:
        (rotmatrix, occupancy)
        """
        self.domains[uc_no] = domains

    def addRelParameter(
        self,
        index0_or_name,
        indexneg_or_name=None,
        initial=None,
        limits=(-np.inf, np.inf),
        prior=None,
        **keyargs,
    ):
        """
        sets a new fit parameter for the weight w of a unit cell (Legacy API!, consider using addWeightParameter)

        F_xtal = w * F_uc1

        if indexneg_or_name is not None the structure factor will be calculated
        to

        F_xtal = w * F_uc1 +  (1 - w) * F_uc2

        Parameters
        ----------
        index0_or_name : str or int
            index or name of the unit cell 1 (with F_uc1)
        indexneg_or_name : default: None, str or int
            index or name of the unit cell 2 (with F_uc2)
        initial: float
            initial value of w
        limits: list
            fit limits of w
        """  # noqa: E501

        index0 = self.getUcIndex(index0_or_name)

        if initial is None:
            initial = self.weights[index0]

        self.weights_0[index0] = 0.0
        self.weights[index0] = initial
        if indexneg_or_name is not None:
            indexneg = self.getUcIndex(indexneg_or_name)
            self.weights_0[indexneg] = 1.0
            self.weights[indexneg] = 1.0 - initial
            names = np.array([index0, indexneg], dtype=np.intp)
            factors = np.array([1.0, -1.0])
        else:
            names = np.array([index0])
            factors = np.array([1.0])

        keyargs["limits"] = limits
        keyargs["prior"] = prior

        return self.addWeightParameter(names, factors, **keyargs)

    def addWeightParameter(self, indices_or_names, factors, **keyargs):
        """
        sets a new fit parameter for the weight w of a unit cell

        Parameters
        ----------
        index0_or_name : str or int
            index or name of the unit cell 1 (with F_uc1)
        indexneg_or_name : default: None, str or int
            index or name of the unit cell 2 (with F_uc2)
        initial: float
            initial value of w
        limits: list
            fit limits of w
        """

        prior = keyargs.get("prior", None)
        name = keyargs.get("name", f"{self.name}_weight_{self._parIdNo}")
        limits = keyargs.get("limits", (-np.inf, np.inf))
        initial = keyargs.get("initial", None)
        basis_init = keyargs.get("basis_init", None)

        indexarray = np.array(
            [self.getUcIndex(ucn) for ucn in indices_or_names], dtype=np.intp
        )

        factors = np.array(factors, dtype=np.float64)

        if basis_init is not None:
            self.weights_0[indexarray] = basis_init

        if initial is not None:
            self.weights[indexarray] = self.weights_0[indexarray] + factors * initial
        par = Parameter(
            name, indexarray, ParameterType.RELATIVE, limits, prior, keyargs, factors
        )
        self.parameters["weight"].append(par)

        self._parIdNo += 1

        self.fit_metadata_cache = None

        return par

    def addFitParameter(self, parameter, limits=(-np.inf, np.inf), **keyargs):
        """Add a coupled fit parameter to crystal components.

        ``parameter`` maps component names to dictionaries containing
        ``atoms`` and ``par`` index selections. An optional ``factors`` entry
        creates a relative parameter, and ``settings`` is forwarded to the
        component fit-parameter constructor.

        Example::

            {
                "film": {
                    "atoms": (1, 2, 3),
                    "par": ("z", "oDW", "iDW"),
                    "factors": (1, -1, 1),
                    "settings": {"name": "film_relaxation"},
                }
            }

        :param dict parameter:
            Component-specific fit-parameter definitions.
        :param tuple limits:
            Lower and upper fit limits.
        :param keyargs:
            Additional coupled-parameter settings.
        :returns:
            The created coupled parameter.
        :rtype: Parameter
        """
        prior = keyargs.get("prior", None)
        name = keyargs.get("name", f"{self.name}_par_{self._parIdNo}")
        keyargs["name"] = name
        keyargs["prior"] = prior

        ucnames = np.array(list(parameter.keys()))
        ucindices = np.array([self.getUcIndex(ucn) for ucn in ucnames])
        sidx = np.argsort(ucindices)

        ucindices = ucindices[sidx]
        ucnames = ucnames[sidx]

        if not isinstance(parameter[ucnames[0]], dict):  #  legacy API fix!
            for uc_identifier in ucnames:
                atoms, par = parameter[uc_identifier]
                parameter[uc_identifier] = {"atoms": atoms, "par": par}

        # children = []

        partype = ParameterType(0)
        for uc_identifier in ucnames:
            par = parameter[uc_identifier]
            settings = par.get("settings", dict())
            settings.update(keyargs)
            if "factors" in par:  # relative parameters
                partype |= ParameterType.RELATIVE
                self[uc_identifier].addRelParameter(
                    (par["atoms"], par["par"]), par["factors"], limits, **settings
                )
                # children.append(fp)
            else:
                partype |= ParameterType.ABSOLUTE
                self[uc_identifier].addFitParameter(
                    (par["atoms"], par["par"]), limits, **settings
                )
                # children.append(fp)

        par = Parameter(name, ucindices, partype, limits, prior, keyargs)
        # par.children = children

        self.parameters["coupled"].append(par)
        self._parIdNo += 1

        self.fit_metadata_cache = None

        return par

    def addWyckoffParameter(self, parameter, limits=(-np.inf, np.inf), **keyargs):
        """Add a coupled symmetry-preserving Wyckoff parameter.

        ``parameter`` maps crystal component names to Wyckoff selections. A
        selection can be ``(site_id, variable)`` or a dictionary with
        ``site_id`` and ``variable`` entries. For ``EpitaxyInterface``
        components, either use a key ``(component, unitcell)`` or provide
        ``unitcell="top"``, ``unitcell="bottom"``, or a list in the selection
        dictionary.

        The coupled value and limits are absolute. This method distributes the
        request to :meth:`UnitCell.addWyckoffParameter` on each selected unit
        cell and links the resulting absolute parameters by name.

        :param dict parameter:
            Component-specific Wyckoff variable selections.
        :param tuple limits:
            Shared absolute limits in fractional units.
        :param keyargs:
            Additional coupled-parameter settings such as ``name`` or
            ``prior``.
        :returns:
            Created coupled parameter.
        :rtype:
            Parameter
        """
        return self._add_wyckoff_coupled_parameter(
            parameter,
            kind="coordinate",
            value_key="variable",
            limits=limits,
            **keyargs,
        )

    def addWyckoffShift(self, parameter, limits=(-np.inf, np.inf), **keyargs):
        """Add coupled relative shifts along parent-cell directions.

        ``parameter`` maps crystal component names to Wyckoff shift
        selections. A selection can be ``(site_id, axis)`` or a dictionary with
        ``site_id`` and ``axis`` entries. For ``EpitaxyInterface`` components,
        either use a key ``(component, unitcell)`` or provide
        ``unitcell="top"``, ``unitcell="bottom"``, or a list in the selection
        dictionary.

        :param dict parameter:
            Component-specific Wyckoff shift selections.
        :param tuple limits:
            Shared relative-shift limits in parent fractional units.
        :param keyargs:
            Additional coupled-parameter settings such as ``name`` or
            ``prior``.
        :returns:
            Created coupled parameter.
        :rtype:
            Parameter
        """
        return self._add_wyckoff_coupled_parameter(
            parameter,
            kind="site_displacement",
            value_key="axis",
            limits=limits,
            **keyargs,
        )

    def _add_wyckoff_coupled_parameter(
        self,
        parameter,
        kind,
        value_key,
        limits=(-np.inf, np.inf),
        **keyargs,
    ):
        if not parameter:
            raise ValueError("At least one Wyckoff target must be provided.")

        prior = keyargs.get("prior", None)
        name = keyargs.get("name", f"{self.name}_par_{self._parIdNo}")
        keyargs["name"] = name
        keyargs["prior"] = prior

        targets = []
        for component_key, selection in parameter.items():
            targets.extend(
                (component_key, *target)
                for target in self._wyckoff_targets(
                    component_key,
                    selection,
                    value_key,
                )
            )

        component_indices = []
        for component_key, component_index, unitcell, site_id, value in targets:
            settings = dict(keyargs)
            if kind == "coordinate":
                unitcell.addWyckoffParameter(
                    site_id,
                    value,
                    limits=limits,
                    **settings,
                )
            else:
                unitcell.addWyckoffShift(
                    site_id,
                    value,
                    limits=limits,
                    **settings,
                )
            component_indices.append(component_index)

        coupled_settings = dict(keyargs)
        coupled_settings["wyckoff"] = {
            "kind": kind,
            "value_kind": "absolute" if kind == "coordinate" else "delta",
        }

        par = Parameter(
            name,
            np.asarray(component_indices, dtype=np.intp),
            ParameterType.RELATIVE,
            limits,
            prior,
            coupled_settings,
        )
        self.parameters["coupled"].append(par)
        self._parIdNo += 1
        self.fit_metadata_cache = None
        return par

    def _wyckoff_targets(self, component_key, selection, value_key):
        internal = None
        if isinstance(component_key, tuple):
            if len(component_key) != 2:
                raise ValueError(
                    "Tuple Wyckoff component keys must be "
                    "(component, unitcell)."
                )
            component_key, internal = component_key

        if isinstance(selection, dict):
            site_id = selection["site_id"]
            value = selection[value_key]
            internal = selection.get("unitcell", internal)
        else:
            site_id, value = selection

        component = self[component_key]
        component_index = self.getUcIndex(component_key)
        if isinstance(internal, list | tuple):
            for internal_name in internal:
                yield (
                    component_index,
                    self._wyckoff_target_unitcell(component, internal_name),
                    site_id,
                    value,
                )
        else:
            yield (
                component_index,
                self._wyckoff_target_unitcell(component, internal),
                site_id,
                value,
            )

    @staticmethod
    def _wyckoff_target_unitcell(component, internal):
        if isinstance(component, UnitCell):
            if internal is not None:
                raise ValueError("UnitCell components do not use 'unitcell'.")
            return component
        if isinstance(component, Film | PoissonSurface):
            if internal is not None:
                return component[internal]
            return component.unitcell
        if isinstance(component, EpitaxyInterface):
            if internal is None:
                raise ValueError(
                    "EpitaxyInterface Wyckoff parameters require "
                    "unitcell='top' or unitcell='bottom'."
                )
            return component[internal]
        if internal is not None and hasattr(component, "__getitem__"):
            return component[internal]
        if hasattr(component, "unitcell"):
            return component.unitcell
        raise TypeError(
            f"Component {component!r} does not expose Wyckoff-bearing unit cells."
        )

    def getSurfaceBasis(self):
        return np.concatenate([uc.basis for uc in self if isinstance(uc, UnitCell)])

    def getSurfaceDWConstraintEnable(self):
        return np.concatenate(
            [uc.dw_increase_constraint for uc in self if isinstance(uc, UnitCell)]
        )

    def fitparameterList(self, verbosity=2):
        s = "## Parameters\n"
        for ucn in self.getUcNames():
            uc = self[ucn]
            s += f"# {uc.__class__.__name__} {uc.name}\n"
            s += uc.fitparameterList(verbosity) + "\n"

        s += "## UnitCell Weights\n"
        s += "{:10} {:10} {:10} {:10} {:10}\n".format(
            "uc1", "occup", "uc2", "occup", "Limits"
        )
        for [index0, indexneg, value], lim in zip(
            self.weightparameters, self.weightparameterLimits
        ):
            uc1name = self.uc_surface_list[index0].name
            uc1value = value
            if indexneg is not None:
                uc2name = self.uc_surface_list[indexneg].name
                uc2value = 1.0 - value
            else:
                uc2name = ""
                uc2value = ""
            s += f"{str(uc1name):10} {str(uc1value):10} {str(uc2name):10} {str(uc2value):10} {str(lim):10}\n"  # noqa: E501

        return s

    def __uc_fitparnames(self):
        names = []
        names += self.uc_bulk.fitparnames
        for uc in self.uc_surface_list:
            names += uc.fitparnames
        return names

    @property
    def fitparnames(self):
        if self.fit_metadata_cache is None:
            self.getStartParamAndLimits()  # will generate metadata, if not yet done so.
        names = [
            par.name
            for par in self.parameters["coupled"]
            + self.parameters["weight"]
            + self.parameters["domain"]
        ]
        names += self.__uc_fitparnames()
        names = list(np.array(names)[self.fit_metadata_cache["par_mask"]])
        return names

    @property
    def priors(self):
        if self.fit_metadata_cache is None:
            self.getStartParamAndLimits()  # will generate metadata, if not yet done so.
        priorlist = []
        for par in (
            self.parameters["coupled"]
            + self.parameters["weight"]
            + self.parameters["domain"]
        ):
            if par.prior is None:
                if tuple(par.limits) == (-np.inf, np.inf):
                    raise Exception(
                        f"Infinite range given as prior for parameter {par.name}. You must provide an explicit prior."  # noqa: E501
                    )
                else:
                    priorlist.append(
                        tuple(par.limits)
                    )  # has to be converted to correct distribution in CTROptimizer
            else:
                priorlist.append(par.prior)  # real prior distribution

        priorlist += self.uc_bulk.priors
        for uc in self.uc_surface_list:
            priorlist += uc.priors
        priorlist = list(np.array(priorlist)[self.fit_metadata_cache["par_mask"]])
        return priorlist

    def parameter_list(self):
        if self.fit_metadata_cache is None:
            self.getStartParamAndLimits()  # will generate metadata, if not yet done so.
        l = (  # noqa: E741
            self.parameters["coupled"]
            + self.parameters["weight"]
            + self.parameters["domain"]
        )
        l += self.uc_bulk.parameter_list()  # noqa: E741
        for uc in self.uc_surface_list:
            l += uc.parameter_list()  # noqa: E741
        l = list(np.array(l)[self.fit_metadata_cache["par_mask"]])  # noqa: E741
        return l

    def getInitialParameters(self, force_recalculate=False):
        p, _, _ = self.getStartParamAndLimits(force_recalculate)
        return p

    def getStartParamAndLimits(self, force_recalculate=False):
        x0 = []
        lower = []
        higher = []
        try:
            mismatch = not np.allclose(
                self._weights_parvalues, self.weights, equal_nan=True
            )
        except Exception:
            mismatch = True
        recalculate = force_recalculate or mismatch

        number_xtalpar = sum(len(p) for p in self.parameters.values())

        par_names = []
        for par in self.parameters["weight"]:
            if recalculate or par.value is None:
                x0.append(
                    np.mean(
                        (self.weights[par.indices] - self.weights_0[par.indices])
                        / par.factors
                    )
                )
            else:
                x0.append(par.value)
            lower.append(par.limits[0])
            higher.append(par.limits[1])
            par_names.append(par.name)

        for par in self.parameters["domain"]:
            raise NotImplementedError("Domain fitting not yet implemented")
            lower.append(par.limits[0])
            higher.append(par.limits[1])
            par_names.append(par.name)

        uc_numberpars = []
        p, low, high = self.uc_bulk.getStartParamAndLimits()
        x0.extend(p)
        lower.extend(low)
        higher.extend(high)

        uc_numberpars.append(len(p))

        for uc in self.uc_surface_list:
            p, low, high = uc.getStartParamAndLimits()
            uc_numberpars.append(len(p))
            x0.extend(p)
            lower.extend(low)
            higher.extend(high)

        par_names += self.__uc_fitparnames()

        number_coupled = len(self.parameters["coupled"])

        par_names = np.array(par_names)
        uc_par_idx = np.arange(par_names.size) + number_coupled
        uc_par_mask = np.ones(par_names.size, dtype=np.bool_)

        x0_coupled = []
        lower_coupled = []
        higher_coupled = []
        idx = 0
        for par in self.parameters["coupled"]:
            paridx = par_names == par.name
            x0_coupled.append(np.mean(np.array(x0)[paridx]))
            uc_par_idx[paridx] = idx
            uc_par_mask[paridx] = False

            lower_coupled.append(par.limits[0])
            higher_coupled.append(par.limits[1])
            idx += 1

        uc_par_idx_trunc = uc_par_idx[uc_par_mask]
        uc_par_idx_trunc_idx = np.arange(uc_par_idx_trunc.size) + number_coupled
        uc_par_idx[uc_par_mask] = uc_par_idx_trunc_idx

        par_idx_sortarray = np.concatenate(
            (np.arange(number_coupled), uc_par_idx)
        ).astype(np.intp)

        self.fit_metadata_cache = {
            "par_idx_sortarray": par_idx_sortarray,
            "par_mask": np.concatenate(
                (np.ones(number_coupled, dtype=np.bool_), uc_par_mask)
            ),
            "number_xtalpar": number_xtalpar,
            "number_coupled": number_coupled,
            "uc_numberpars": uc_numberpars,
        }

        x0 = np.array(x0)[uc_par_mask]
        lower = np.array(lower)[uc_par_mask]
        higher = np.array(higher)[uc_par_mask]

        return (
            np.concatenate((x0_coupled, x0)),
            np.concatenate((lower_coupled, lower)),
            np.concatenate((higher_coupled, higher)),
        )

    """
    def _readFitparValues(self):

        updates fitp_values from fitparameters

        reads parameter values in self.fitparameters and writes
        the mean value to the corresponding entry in fitp_values

        Returns
        -------
        None.


        values = np.zeros(len(self.fitp_values))
        count = np.zeros(len(self.fitp_values))
        for uc_pars in self.fitparameters:
            for p in uc_pars:
                values[int(p[0])] += p[2]
                count[int(p[0])] += 1.

        self.fitp_values = values/count
    """

    def setParameters(self, x):
        if self.fit_metadata_cache is None:
            self.getStartParamAndLimits()  # will generate metadata, if not yet done so.
        x = np.asarray(x)
        x_r = x[
            self.fit_metadata_cache["par_idx_sortarray"]
        ]  # reorder fitpars, including coupled absolute Wyckoff parameters
        x_ucs = x_r[self.fit_metadata_cache["number_xtalpar"] :]
        uc_numberpars = self.fit_metadata_cache["uc_numberpars"]

        lower = 0
        for uc, noparam in zip([self.uc_bulk] + self.uc_surface_list, uc_numberpars):
            uc.setFitParameters(x_ucs[lower : noparam + lower])
            lower = noparam + lower

        for i, par in enumerate(self.parameters["coupled"]):
            par.value = x[i]

        self.weights = np.copy(self.weights_0)
        idx = self.fit_metadata_cache["number_coupled"]
        for par in self.parameters["weight"]:
            self.weights[par.indices] += par.factors * x_r[idx]
            par.value = x_r[idx]
            idx += 1

        self._weights_parvalues = np.copy(self.weights)

    def setLimits(self, lim):
        """fit bounds.

        lim shape: (n, 2)

        """
        if self.fit_metadata_cache is None:
            self.getStartParamAndLimits()  # will generate metadata, if not yet done so.
        lim = np.asarray(lim)
        x_r = lim[
            self.fit_metadata_cache["par_idx_sortarray"]
        ]  # reorder limits, including coupled absolute Wyckoff parameters
        x_ucs = x_r[self.fit_metadata_cache["number_xtalpar"] :]
        uc_numberpars = self.fit_metadata_cache["uc_numberpars"]

        lower = 0
        for uc, noparam in zip([self.uc_bulk] + self.uc_surface_list, uc_numberpars):
            uc.setLimits(x_ucs[lower : noparam + lower])
            lower = noparam + lower

        for i, par in enumerate(self.parameters["coupled"]):
            if lim[i][0] > lim[i][1]:
                raise ValueError(
                    "Upper fit bound of coupled par must be larger than the lower bound."  # noqa: E501
                )
            par.limits = (lim[i][0], lim[i][1])

        idx = self.fit_metadata_cache["number_coupled"]
        for par in self.parameters["weight"]:
            if x_r[idx][0] > x_r[idx][1]:
                raise ValueError(
                    "Upper fit bound of weight pars must be larger than the lower bound."  # noqa: E501
                )
            par.limits = (x_r[idx][0], x_r[idx][1])
            idx += 1

    def setFitErrors(self, errors):
        if self.fit_metadata_cache is None:
            self.getStartParamAndLimits()  # will generate metadata, if not yet done so.
        errors = np.asarray(errors)
        errors_r = errors[
            self.fit_metadata_cache["par_idx_sortarray"]
        ]  # reorder fitpars
        errors_ucs = errors_r[self.fit_metadata_cache["number_xtalpar"] :]
        uc_numberpars = self.fit_metadata_cache["uc_numberpars"]

        lower = 0
        for uc, noparam in zip([self.uc_bulk] + self.uc_surface_list, uc_numberpars):
            uc.setFitErrors(errors_ucs[lower : noparam + lower])
            lower = noparam + lower

        for i, par in enumerate(self.parameters["coupled"]):
            par.error = errors[i]

        self.werrors = np.full_like(self.weights, np.nan)
        idx = 0
        for par in self.parameters["weight"]:
            self.werrors[par.indices] = np.nan_to_num(
                self.werrors[par.indices], nan=0.0
            )
            self.werrors[par.indices] += par.factors * errors_r[idx]
            par.error = errors_r[idx]
            idx += 1
        self._werrors_parvalues = np.copy(self.werrors)

    def getFitErrors(self):
        if self.fit_metadata_cache is None:
            self.getStartParamAndLimits()  # will generate metadata, if not yet done so.
        if self.werrors is None:
            raise ValueError(
                "No errors for SXRDCrystal UnitCell weights have been set."
            )
        try:
            mismatch = not np.allclose(
                self._werrors_parvalues, self.werrors, equal_nan=True
            )
        except Exception:
            mismatch = True
        err0 = []

        sum(len(p) for p in self.parameters.values())

        par_names = []
        for par in self.parameters["weight"]:
            if mismatch:
                err0.append(np.mean(self.werrors[par.indices] / np.abs(par.factors)))
            else:
                err0.append(par.error)
            par_names.append(par.name)

        for par in self.parameters["domain"]:
            raise NotImplementedError("Domain fitting not yet implemented")
            # lower.append(par.limits[0]); higher.append(par.limits[1])
            # par_names.append(par.name)

        # uc_numberpars = []
        p = self.uc_bulk.getFitErrors()
        err0.extend(p)

        # uc_numberpars.append(len(p))

        for uc in self.uc_surface_list:
            p = uc.getFitErrors()
            # uc_numberpars.append(len(p))
            err0.extend(p)

        par_names += self.__uc_fitparnames()

        par_names = np.array(par_names)

        err0_coupled = []
        for par in self.parameters["coupled"]:
            paridx = par_names == par.name
            err0_coupled.append(np.mean(np.array(err0)[paridx]))

        err = np.concatenate((err0_coupled, err0))[self.fit_metadata_cache["par_mask"]]

        return err

    def setEnergy(self, E):
        for uc in self.uc_surface_list:
            uc.setEnergy(E)
        self.uc_bulk.setEnergy(E)

    def zDensity_G(self, z, h, k):
        if self.enable_uc_stacking:
            self.apply_stacking()
        rho = self.uc_bulk.zDensity_G_asbulk(z, h, k)
        for uc, w in zip(self.uc_surface_list, self.weights):
            rho += uc.zDensity_G(z, h, k) * w
        return rho

    def optical_profile(self):
        """Return the combined homogeneous optical profile of this crystal.

        Crystal surface weights are applied directly to ``delta`` and ``beta``.
        Water is sampled on the top atomistic layer grid. Without water, one
        vacuum row is appended one top-layer spacing above the highest atomic
        layer.

        :returns:
            C-contiguous ``(N, 3)`` array with columns ``z`` in Angstrom,
            ``delta``, and ``beta``.
        :rtype: numpy.ndarray
        :raises NotImplementedError:
            If a component does not expose an atomistic optical profile.
        """
        if self.enable_uc_stacking:
            self.apply_stacking()
        profiles = [self.uc_bulk.optical_profile_asbulk()]
        water_components = []
        for component, weight in zip(self.uc_surface_list, self.weights):
            if not hasattr(component, "optical_profile"):
                raise NotImplementedError(
                    f"{type(component).__name__} has no atomistic optical profile."
                )
            if isinstance(component, WaterModel):
                water_components.append((component, weight))
            else:
                profile = component.optical_profile().copy()
                profile[:, 1:] *= weight
                profiles.append(profile)
        structural_profile = combine_profiles(*profiles)
        dz = top_layer_spacing(structural_profile)
        if not water_components:
            vacuum = np.array(
                [[structural_profile[-1, 0] + dz, 0.0, 0.0]], dtype=np.float64
            )
            return np.ascontiguousarray(
                np.concatenate((structural_profile, vacuum), axis=0)
            )
        water_profiles = []
        z_origin = structural_profile[-1, 0]
        for component, weight in water_components:
            profile = component.optical_profile(
                z_step=dz, z_origin=z_origin
            ).copy()
            profile[:, 1:] *= weight
            water_profiles.append(profile)
        return add_structural_to_sampled_profile(
            structural_profile, *water_profiles
        )

    def simplified_optical_profile(
        self, delta_tolerance=1e-9, beta_tolerance=None
    ):
        """Return a thickness-conserving simplified optical profile.

        :param float delta_tolerance:
            Maximum delta range within a merged finite layer.
        :param float beta_tolerance:
            Maximum beta range within a merged finite layer. Defaults to the
            delta tolerance.
        :returns:
            Simplified ``(N, 3)`` optical profile.
        :rtype: numpy.ndarray
        """
        return simplify_profile(
            self.optical_profile(), delta_tolerance, beta_tolerance
        )

    def stratified_optical_profile(
        self, delta_tolerance=1e-9, beta_tolerance=None
    ):
        """Return optical media and their physical interface positions.

        Profile z coordinates are centers of homogeneous media. Interfaces
        are centered between adjacent samples before optional simplification;
        merged media retain those original outer boundaries.

        :param float delta_tolerance:
            Maximum delta range within a merged finite layer. Pass ``None``
            to retain every sampled medium.
        :param float beta_tolerance:
            Maximum beta range within a merged finite layer. Defaults to the
            delta tolerance when simplification is requested.
        :returns:
            Layer-center optical constants and substrate-to-ambient edges.
        :rtype: CTRoptics.StratifiedProfile
        """
        return stratify_profile(
            self.optical_profile(), delta_tolerance, beta_tolerance
        )

    def wavefield(
        self,
        alpha,
        polarization="s",
        delta_tolerance=1e-9,
        beta_tolerance=None,
    ):
        """Return the unperturbed layered electric field.

        The field follows the Renaud convention with downward and upward
        amplitudes ``A_plus`` and ``A_minus``. Slabs are ordered from the
        incident medium towards the substrate.

        :param float alpha:
            Glancing incidence angle in degrees inside the incident medium.
        :param str polarization:
            ``"s"`` or ``"p"``.
        :param float delta_tolerance:
            Maximum delta range used to simplify adjacent layers. Pass
            ``None`` to disable simplification.
        :param float beta_tolerance:
            Optional beta simplification tolerance. Defaults to the delta
            tolerance when simplification is requested.
        :returns:
            Layer amplitudes, normal wavevectors, and sampled electric field.
        :rtype: CTRoptics.Wavefield
        """
        profile = self.stratified_optical_profile(
            delta_tolerance, beta_tolerance
        )
        return solve_wavefield(
            profile.values,
            self.uc_bulk._E,
            alpha,
            polarization,
            boundaries=profile.boundaries,
        )

    def specular_reflectivity(
        self,
        alpha,
        polarization="s",
        delta_tolerance=1e-9,
        beta_tolerance=None,
    ):
        """Return scalar optical specular reflectivity.

        :param float or numpy.ndarray alpha:
            Glancing incidence angle or angles in degrees inside the incident
            medium.
        :param str polarization:
            ``"s"``, ``"p"``, or ``"unpolarized"``. Unpolarized intensity
            is the equal incoherent average of s and p reflectivity.
        :param float delta_tolerance:
            Maximum delta range used to simplify adjacent layers. Pass
            ``None`` to disable simplification.
        :param float beta_tolerance:
            Optional beta simplification tolerance. Defaults to the delta
            tolerance when simplification is requested.
        :returns:
            Reflectivity with the same scalar/array shape as ``alpha``.
        :rtype: float or numpy.ndarray
        """
        if polarization not in {"s", "p", "unpolarized"}:
            raise ValueError(
                "polarization must be 's', 'p', or 'unpolarized'."
            )
        alpha_array = np.asarray(alpha, dtype=np.float64)
        scalar = alpha_array.ndim == 0
        angles = np.atleast_1d(alpha_array)
        profile = self.stratified_optical_profile(
            delta_tolerance, beta_tolerance
        )

        def reflectivity(pol):
            return np.asarray(
                [
                    abs(
                        solve_wavefield(
                            profile.values,
                            self.uc_bulk._E,
                            angle,
                            pol,
                            boundaries=profile.boundaries,
                        ).r_S
                    )
                    ** 2
                    for angle in angles
                ],
                dtype=np.float64,
            )

        if polarization == "unpolarized":
            result = 0.5 * (reflectivity("s") + reflectivity("p"))
        else:
            result = reflectivity(polarization)
        if scalar:
            return float(result[0])
        return result.reshape(alpha_array.shape)

    def toRODStr(self):
        s = f"E = {self.uc_bulk._E * 1e-3:.5f} keV\n"
        for i, w, uc in zip(self.uc_stacking, self.weights, self.uc_surface_list):
            s += f"# {uc.__class__.__name__} {uc.name}\n"
            s += f"{i:04d} {w:.5f}\n"
            s += uc.toRODStr() + "\n"

        s += "# UnitCell bulk\n" + self.uc_bulk.toRODStr()
        return s

    def toStr(self, showErrors=True):
        """Serialize the complete crystal model as plain text.

        :param bool showErrors:
            Include propagated component, atom, and weight errors.
        :returns:
            Plain-text crystal representation.
        :rtype: str
        """
        s = f"E = {self.uc_bulk._E * 1e-3:.5f} keV\n"
        for i, (no, w, uc) in enumerate(
            zip(self.uc_stacking, self.weights, self.uc_surface_list)
        ):
            s += f"# {uc.__class__.__name__} {uc.name}\n"
            if showErrors and self.werrors is not None:
                s += f"{no:04d} occupancy = {w:.5f} +- {self.werrors[i]:.5f}\n"
            else:
                s += f"{no:04d} occupancy = {w:.5f}\n"
            s += uc.toStr(showErrors=showErrors) + "\n"

        s += "# UnitCell bulk\n" + self.uc_bulk.toStr(showErrors=showErrors)
        return s

    @staticmethod
    def fromStr(string, atten=0.01):
        weights = []
        werrors = []
        uc_suface = []
        uc_stacking = []

        strio = util.StringIO(string)
        Estr = next_skip_comment(strio).split("=")[1].split("keV")[0]
        E = float(Estr) * 1e3
        string = strio.read()
        errors = False
        splstring = string.split("#")
        unitCells = splstring[1:]
        for i in range(
            len(splstring) - 1, 0, -1
        ):  # ugly fix for commenting of line with #
            lineprefix = splstring[i - 1][
                splstring[i - 1].rfind("\n") :
            ]  # get contents of line before '#'
            if "//" in lineprefix:
                del unitCells[
                    i - 1
                ]  # index of unitCells is already shifted by 1 in "unitCells = splstring[1:]"  # noqa: E501

        for ucstr in unitCells[:-1]:
            strio = util.StringIO(ucstr)
            line1 = next_skip_comment(strio)
            classname, name = line1.split(maxsplit=1)
            name = name.strip()
            line2 = next_skip_comment(strio).split()
            uc_stacking.append(int(line2[0]))
            if "+-" in line2:
                errors = True
                weights.append(float(line2[line2.index("+-") - 1]))
                werrors.append(float(line2[line2.index("+-") + 1]))
            elif "=" in line2:
                weights.append(float(line2[line2.index("=") + 1]))
                werrors.append(float("nan"))
            else:
                weights.append(float(line2[1]))
                werrors.append(float("nan"))
            if classname == "UnitCell":
                uc = UnitCell.fromStr(strio.read())
            elif classname == "WaterModel":
                uc = WaterModel.fromStr(strio.read())
            elif classname == "EpitaxyInterface":
                uc = EpitaxyInterface.fromStr(strio.read())
            elif classname == "Film":
                uc = Film.fromStr(strio.read())
            elif classname == "PoissonSurface":
                uc = PoissonSurface.fromStr(strio.read())
            else:
                raise NotImplementedError(f"class name not understood: {classname}")
            uc.name = name
            uc.setEnergy(E)
            uc_suface.append(uc)
        bulkstrio = util.StringIO(unitCells[-1])
        next_skip_comment(bulkstrio)
        uc_bulk = UnitCell.fromStr(bulkstrio.read())
        uc_bulk.setEnergy(E)
        uc_stacking = np.array(uc_stacking)
        xtal = SXRDCrystal(uc_bulk, *uc_suface, atten=atten, stacking=uc_stacking)
        xtal.weights = np.array(weights)
        xtal.weights_0 = np.copy(weights)
        xtal.werrors = np.array(werrors) if errors else None
        return xtal

    def __getitem__(self, uc_name_or_index):
        if isinstance(uc_name_or_index, int):
            if uc_name_or_index == -1:
                return self.uc_bulk
            return self.uc_surface_list[uc_name_or_index]
        elif isinstance(uc_name_or_index, str):
            if uc_name_or_index == "bulk":
                return self.uc_bulk
            namelist = [uc.name for uc in self.uc_surface_list]
            idx = namelist.index(uc_name_or_index)
            return self.uc_surface_list[idx]
        else:
            raise ValueError(f"must be str or int, not {type(uc_name_or_index)}")

    def getUcIndex(self, uc_name_or_index):
        if isinstance(uc_name_or_index, int | np.integer):
            return int(uc_name_or_index)
        elif isinstance(uc_name_or_index, str):
            if uc_name_or_index == "bulk":
                return -1
            namelist = [uc.name for uc in self.uc_surface_list]
            idx = namelist.index(uc_name_or_index)
            return idx
        else:
            raise ValueError(f"must be str or int, not {type(uc_name_or_index)}")

    def plot3d(
        self,
        ucx=2,
        ucy=2,
        ucz=-3,
        dwon=False,
        occuon=False,
        figure=None,
        translate=np.array([0.0, 0.0, 0.0]),
        backend="auto",
        radius_scale=1.0,
        **keyargs,
    ):
        """Plot the bulk and surface unit cells in one 3D viewer.

        Coordinates and radii passed to either backend are in Angstrom.
        Assign the result or terminate an incremental call with a semicolon to
        prevent Jupyter from also rendering the returned view as a new output.

        :param int ucx: Number of cells along the first lattice direction.
        :param int ucy: Number of cells along the second lattice direction.
        :param int ucz: Number of bulk cells along the third lattice direction.
        :param bool dwon:
            Sample atom positions using the stored Debye-Waller disorder.
        :param bool occuon:
            Randomly omit atoms according to their occupancies.
        :param figure:
            Existing Mayavi figure or py3Dmol view to extend.
        :param numpy.ndarray translate:
            Fractional translation retained for API compatibility.
        :param str backend:
            ``"auto"``, ``"mayavi"``, or ``"py3dmol"``.
        :param float radius_scale:
            Positive dimensionless multiplier applied to covalent radii.
        :returns: The selected backend's shared figure or view.
        """
        keyargs["_defer_update"] = True

        for mat in self.uc_bulk.coherentDomainMatrix:
            figure = self.uc_bulk.plot3d(
                ucx,
                ucy,
                ucz,
                dwon,
                occuon,
                figure,
                mat,
                backend=backend,
                radius_scale=radius_scale,
                **keyargs,
            )

        for uc in self.uc_surface_list:
            for mat in uc.coherentDomainMatrix:
                figure = uc.plot3d(
                    ucx,
                    ucy,
                    1,
                    dwon,
                    occuon,
                    figure,
                    mat,
                    backend=backend,
                    radius_scale=radius_scale,
                    **keyargs,
                )

        if getattr(figure, "_orgui_plot3d_backend", None) == "py3dmol":
            figure.zoomTo()
            if getattr(figure, "uniqueid", None) is not None:
                figure.update()
        return figure

    def getUcNames(self):
        return [uc.name for uc in self.uc_surface_list] + ["bulk"]

    def _ipython_key_completions_(self):
        return self.getUcNames()

    def toFile(self, filename):
        """Write a plain-text ``.xtal`` or ``.xpr`` model.

        ``.xtal`` files store fitted values without propagated errors.
        ``.xpr`` files additionally store propagated errors for weights,
        CTRfilm parameters, and atom parameters.

        :param str filename:
            Output path ending in ``.xtal`` or ``.xpr``.
        :raises ValueError:
            If the filename has another extension.
        """
        extension = os.path.splitext(filename)[1].lower()
        if extension not in (".xtal", ".xpr"):
            raise ValueError("CTR model filename must end in .xtal or .xpr")
        with open(filename, "w") as f:
            f.write(self.toStr(showErrors=extension == ".xpr"))

    def toRODFile(self, filename):
        with open(filename, "w") as f:
            f.write(self.toRODStr())

    @staticmethod
    def fromFile(filename):
        """Read a plain-text ``.xtal`` or ``.xpr`` crystal model.

        Values and propagated errors are reconstructed from the text. Fit
        parameter definitions remain in the optional companion ``.h5`` file
        used by the existing fitting workflow.

        :param str filename:
            Model path, or basename for automatic ``.xpr``/``.xtal`` lookup.
        :returns:
            Reconstructed crystal model.
        :rtype: SXRDCrystal
        """
        f, ext = os.path.splitext(filename)
        if ext == "":
            for fext in [".xpr", ".xtal"]:
                if os.path.isfile(filename + fext):
                    with open(filename + fext) as f:
                        fstr = f.read()
                    xtal = SXRDCrystal.fromStr(fstr)
                    xtal.name = os.path.basename(filename + fext)
                    break
            else:
                raise FileNotFoundError(
                    errno.ENOENT, os.strerror(errno.ENOENT), filename
                )
            if os.path.isfile(filename + ".h5"):
                try:
                    from silx.io import dictdump

                    fitdict = dictdump.h5todict(filename + ".h5")
                    xtal.parametersFromDict(fitdict)
                except Exception as e:
                    warnings.warn(
                        f"Fitsettings for SXRDCrystal {os.path.basename(filename)} seem to exist, but cannot be read: {e}"  # noqa: E501
                    )
        else:
            with open(filename) as f:
                fstr = f.read()
            xtal = SXRDCrystal.fromStr(fstr)
            xtal.name = os.path.basename(filename)
        return xtal

    def __repr__(self):
        return self.toStr()
