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
import xraydb
import warnings
import json

# random.seed(45)
from collections import OrderedDict
import dataclasses
from dataclasses import dataclass, field
import typing
import enum
from abc import ABC


special_elementcolors = {"Pt2+": (1.0, 0.75, 0.0), "Pt4+": (0.0, 0.0, 1.0)}

special_formfactors = OrderedDict()
# special_eDensity = OrderedDict([('H2O' , eDens_water_fun )])
special_eDensity = OrderedDict()
__db = xraydb.get_xraydb()
if hasattr(__db, "atomic_symbols"):  # XrayDB <= 4.5
    special_onset = len(__db.atomic_symbols)
else:
    __elems = __db.get_cache("elements")
    special_onset = len([e.element for e in __elems])
__db.close()

special_numbers = OrderedDict(
    zip(special_formfactors.keys(), np.arange(len(special_formfactors)) + special_onset)
)


class ParameterType(enum.IntFlag):
    UNDEFINED = enum.auto()
    ABSOLUTE = enum.auto()
    RELATIVE = enum.auto()


@dataclass
class Parameter:
    name: str
    indices: np.ndarray
    kind: ParameterType
    limits: tuple = (-np.inf, np.inf)
    prior: typing.Any = None
    settings: dict = field(default_factory=dict)
    factors: np.ndarray = None
    # children: list = field(default_factory=list)
    value: float = None
    error: float = None

    def __post_init__(self):  # ensures correct data types
        for data_field in dataclasses.fields(self):
            value = getattr(self, data_field.name)
            if data_field.name == "indices":
                idx = np.array(value, dtype=np.intp)
                if len(idx.shape) == 2:
                    idx = tuple(idx)
                setattr(self, data_field.name, idx)
            elif data_field.name == "factors":
                if value is not None:
                    factors = np.array(value, dtype=np.float64)
                    setattr(self, data_field.name, factors)
            elif data_field.name in ["prior"]:
                pass  # ToDo convert to distribution here!
            else:
                if value is not None:
                    setattr(self, data_field.name, data_field.type(value))

    def asdict(self):
        d = dataclasses.asdict(self)
        prior_a = np.array(d["prior"])
        if not np.issubdtype(prior_a.dtype, np.number):
            d["prior"] = None
        # d['children'] = [c.asdict() for c in d['children']]
        for p in ["indices_t", "factors_t"]:
            d.pop(p, None)
        return d


class LinearFitFunctions(ABC):
    parameterOrder = ""
    parameterLookup = {}
    parameterLookup_inv = dict(map(reversed, parameterLookup.items()))

    def __init__(self):
        self.basis = np.array([])
        self.parameters = {"absolute": [], "relative": []}
        self.basis_0 = np.array([])
        self._basis_parvalues = None
        self.errors = None
        self._errors_parvalues = None

    def parametersToDict(self):
        d = dict()
        d["basis_0"] = self.basis_0
        for p in self.parameters:
            d[p] = dict()
            for i, param in enumerate(self.parameters[p]):
                d[p][f"{i}_{param.name}"] = param.asdict()
        return d

    def clearParameters(self):
        for p in self.parameters:
            self.parameters[p] = []
        self._basis_parvalues = None
        self._errors_parvalues = None

    def parametersFromDict(self, d, override_values=True):
        for p in self.parameters:
            self.parameters[p] = []

        self.basis_0 = d["basis_0"].astype(np.float64)
        for p in self.parameters:
            for dkey in sorted(d[p].keys()):
                self.parameters[p].append(Parameter(**d[p][dkey]))

        if override_values:
            try:
                self.updateFromParameters()
            except Exception as e:
                warnings.warn(
                    f"Can not apply fit parameter values to this {self.__class__.__name__}: {e}"  # noqa: E501
                )
        else:
            warnings.warn(
                "Basis values not updated from Parameter values. This can cause a mismatch between basis and fitparameter values!"  # noqa: E501
            )

    def updateFromParameters(self):
        """Update basis and propagated errors from stored fit parameters.

        Absolute parameters replace their selected basis values. Relative
        parameters add ``factors * value`` to ``basis_0``. Stored parameter
        errors are propagated with the absolute value of their factors.
        """
        no_errors = False
        self.basis = np.copy(self.basis_0)
        self.errors = np.full_like(self.basis, np.nan)
        for par in self.parameters["absolute"]:
            if par.value is not None:
                self.basis[par.indices] = par.value
            else:
                raise ValueError(
                    f"Can not set basis values from parameters. Value of Parameter {par.name} is None."  # noqa: E501
                )
            if par.error is not None:
                self.errors[par.indices] = par.error
            else:
                no_errors = True
        for par in self.parameters["relative"]:
            if par.value is not None:
                self.basis[par.indices] += par.factors * par.value
            else:
                raise ValueError(
                    f"Can not set basis values from parameters. Value of Parameter {par.name} is None."  # noqa: E501
                )
            if par.error is not None:
                self.errors[par.indices] = np.nan_to_num(
                    self.errors[par.indices], nan=0.0
                )
                self.errors[par.indices] += np.abs(par.factors) * par.error
            else:
                no_errors = True

        self._basis_parvalues = np.copy(self.basis)
        if not no_errors:
            self._errors_parvalues = np.copy(self.errors)

    def addFitParameter(self, indexarray, limits=(-np.inf, np.inf), **keyargs):
        par = indexarray
        try:
            parindexes = np.array(
                [
                    p if isinstance(p, np.integer | int) else self.parameterLookup[p]
                    for p in np.atleast_1d(par)
                ],
                dtype=np.intp,
            )
        except Exception as e:
            raise ValueError(
                f"Invalid parameter name {par}., Parameters for class {self.__class__.__name__} are {list(self.parameterLookup.keys())}"  # noqa: E501
            ) from e

        if "name" in keyargs:
            name = keyargs["name"]
        else:
            parameternames = tuple([self.parameterLookup_inv[n] for n in parindexes])
            name = self.name + " " + " ".join([p for p in parameternames])

        if name in self.fitparnames:
            raise ValueError(f"Absolute fit parameter {name} already exists.")

        indexarray = parindexes
        try:
            curr_val = self.basis[indexarray]
        except IndexError as e:
            raise ValueError("Invalid parameter indices.") from e

        if not limits[0] <= np.mean(curr_val) <= limits[1]:
            raise ValueError(
                f"start parameter is not within limits! Is not: {limits[0]} <= {np.mean(self.basis[indexarray])} <= {limits[1]}"  # noqa: E501
            )
        prior = keyargs.get("prior", None)
        keyargs["name"] = name
        keyargs["prior"] = prior

        par = Parameter(
            name, indexarray, ParameterType.ABSOLUTE, limits, prior, keyargs
        )
        self.parameters["absolute"].append(par)
        return par

    def addRelParameter(self, indexarray, factors, limits=(-np.inf, np.inf), **keyargs):
        par = indexarray
        try:
            parindexes = np.array(
                [
                    p if isinstance(p, np.integer | int) else self.parameterLookup[p]
                    for p in np.atleast_1d(par)
                ],
                dtype=np.intp,
            )
        except Exception as e:
            raise ValueError(
                f"Invalid parameter name {par}., Parameters for class {self.__class__.__name__} are {list(self.parameterLookup.keys())}"  # noqa: E501
            ) from e

        if "name" in keyargs:
            name = keyargs["name"]
        else:
            parameternames = tuple(
                [self.parameterLookup_inv[n] + "_r" for n in parindexes]
            )
            name = self.name + " " + " ".join([p for p in parameternames])

        if name in self.fitparnames:
            raise ValueError(f"Relative fit parameter {name} already exists.")

        prior = keyargs.get("prior", None)
        keyargs["name"] = name
        keyargs["prior"] = prior

        factors = np.atleast_1d(factors)
        try:
            curr_val = self.basis[parindexes]
        except IndexError as e:
            raise ValueError("Invalid parameter indices.") from e

        if curr_val.shape != factors.shape:
            raise ValueError(
                "Number of basis parameters does not match number of factors."
            )

        par = Parameter(
            name, parindexes, ParameterType.RELATIVE, limits, prior, keyargs, factors
        )
        self.parameters["relative"].append(par)
        return par

    def getInitialParameters(self, force_recalculate=False):
        x0, _, _ = self.getStartParamAndLimits(force_recalculate)
        return x0

    def getStartParamAndLimits(self, force_recalculate=False):
        # if self.basis_0 is None:
        #    self.basis_0 = np.copy(self.basis)
        try:
            mismatch = not np.allclose(
                self._basis_parvalues, self.basis, equal_nan=True
            )
        except Exception:
            mismatch = True
        recalculate = force_recalculate or mismatch
        abspar = self.parameters["absolute"]
        relpar = self.parameters["relative"]

        x0 = []
        lower = []
        upper = []
        for par in abspar:
            if recalculate or par.value is None:
                x0.append(np.mean(self.basis[par.indices]))
            else:
                x0.append(par.value)
            lower.append(par.limits[0])
            upper.append(par.limits[1])

        for par in relpar:
            if recalculate or par.value is None:
                x0.append(
                    np.mean(
                        (self.basis[par.indices] - self.basis_0[par.indices])
                        / par.factors
                    )
                )
            else:
                x0.append(par.value)
            lower.append(par.limits[0])
            upper.append(par.limits[1])

        return np.array(x0), np.array(lower), np.array(upper)

    def fitparameterList(self, verbosity=3):
        if verbosity > 1:
            outstr = "{}".format("# Direct fit parameters:\n")
        else:
            outstr = ""

        if len(self.parameters["absolute"]) > 0:
            outstr += "{:10} {:20} {:10} {:10}".format(
                "Id", "parameters", "Value", "Limits"
            )

            for par in self.parameters["absolute"]:
                outstr += "\n" + self.fitparToStr(par)
        elif verbosity > 2:
            outstr += "# No direct fit parameters"

        outstr += "\n{}".format("# Relative fit parameters:\n")

        if len(self.parameters["relative"]) > 0:
            outstr += "{:10} {:20} {:10} {:10} {:10}".format(
                "Id", "parameters", "Value", "Factors", "Limits"
            )

            for par in self.parameters["relative"]:
                outstr += "\n" + self.fitparToStr(par)
        elif verbosity > 2:
            outstr += "# No rel. parameters"

        return outstr

    def fitparToStr(self, par):
        if self.basis_0 is None:
            basis_0 = np.copy(self.basis)
        else:
            basis_0 = self.basis_0

        if par.kind & ParameterType.ABSOLUTE:
            val = np.mean(self.basis[par.indices])
            parameternames = tuple([self.parameterLookup_inv[n] for n in par.indices])
            return (
                f"{par.name:10} {str(parameternames):15} {val:20} {str(par.limits):10}"
            )
        elif par.kind & ParameterType.RELATIVE:
            val = np.mean(
                (self.basis[par.indices] - basis_0[par.indices]) / par.factors
            )
            parameternames = tuple(
                [self.parameterLookup_inv[n] + "_r" for n in par.indices]
            )

            return f"{par.name:10} {str(parameternames):20} {val:10} {str(par.factors):10} {str(par.limits):10}"  # noqa: E501
        else:
            raise ValueError(
                f"Unvalid parameter type {par.kind} for class {self.__class__.__name__}"
            )

    def setFitParameters(self, x):
        x_0 = x[: len(self.parameters["absolute"])]
        x_r = x[len(self.parameters["absolute"]) :]
        self.basis[:] = self.basis_0
        for val, par in zip(x_0, self.parameters["absolute"]):
            self.basis[par.indices] = val
            par.value = val
        for val, par in zip(x_r, self.parameters["relative"]):
            self.basis[par.indices] += par.factors * val
            par.value = val
        self._basis_parvalues = np.copy(self.basis)

    def setLimits(self, lim):
        x_0 = lim[: len(self.parameters["absolute"])]
        x_r = lim[len(self.parameters["absolute"]) :]
        for val, par in zip(x_0, self.parameters["absolute"]):
            if val[0] > val[1]:
                raise ValueError("Upper fit bound must be larger than the lower bound.")
            par.limits = (val[0], val[1])
        for val, par in zip(x_r, self.parameters["relative"]):
            if val[0] > val[1]:
                raise ValueError("Upper fit bound must be larger than the lower bound.")
            par.limits = (val[0], val[1])

    def setFitErrors(self, errors):
        self.errors = np.full_like(self.basis, np.nan)
        err_0 = errors[: len(self.parameters["absolute"])]
        err_r = errors[len(self.parameters["absolute"]) :]
        for val, par in zip(err_0, self.parameters["absolute"]):
            self.errors[par.indices] = val
            par.error = val
        for val, par in zip(err_r, self.parameters["relative"]):
            self.errors[par.indices] = np.nan_to_num(self.errors[par.indices], nan=0.0)
            self.errors[par.indices] += np.abs(par.factors) * val
            par.error = val
        self._errors_parvalues = np.copy(self.errors)

    def getFitErrors(self):
        if self.errors is None:
            raise ValueError("No errors have been set.")
        try:
            mismatch = not np.allclose(
                self._errors_parvalues, self.errors, equal_nan=True
            )
        except Exception:
            mismatch = True
        abspar = self.parameters["absolute"]
        relpar = self.parameters["relative"]
        err0 = []
        for par in abspar:
            if mismatch:
                err0.append(np.mean(self.errors[par.indices]))
            else:
                err0.append(par.error)

        for par in relpar:
            if mismatch:
                err0.append(np.mean(self.errors[par.indices] / np.abs(par.factors)))
            else:
                err0.append(par.error)

        err0 = np.array(err0)
        if np.any(~np.isfinite(err0)):
            raise ValueError("Some errors are non-finite or not set.")
        return err0

    def parameter_list(self):
        return self.parameters["absolute"] + self.parameters["relative"]

    @property
    def fitparnames(self):
        # pars = [self.fitparToStr(i,True) for i in range(len(self.fitparameters))] + [self.fitparToStr(i + len(self.fitparameters),True) for i in range(len(self.relfitparam))]  # noqa: E501
        pars = []
        for p_kind in self.parameters:
            pars += [p.name for p in self.parameters[p_kind]]
        return pars

    @property
    def priors(self):
        priorlist = []
        pars = []
        for p_kind in self.parameters:
            pars += [p for p in self.parameters[p_kind]]
        for par in pars:
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
        return priorlist


def _ensure_contiguous(*arrays, testOnly=False, astype=None):
    if testOnly:
        for a in arrays:
            if not a.flags["C_CONTIGUOUS"]:
                raise ValueError(
                    "ndArray is not contiguous. Use np.ascontiguousarray() to convert the array"  # noqa: E501
                )
            if astype is not None and not a.dtype == astype:
                raise ValueError(
                    f"Wrong data type of ndArray. Should be {astype}, but is {a.dtype}"
                )
        return arrays
    else:
        a_c = []
        for a in arrays:
            if not a.flags["C_CONTIGUOUS"]:
                warnings.warn(
                    "h,k,l is not contiguous. This causes calculation speed loss. Use np.ascontiguousarray() to convert the arrays"  # noqa: E501
                )
                a_c.append(np.ascontiguousarray(a, dtype=astype))
            elif astype is not None and not a.dtype == astype:
                warnings.warn(
                    f"Not matching data type of ndArray. Should be {astype}, but is {a.dtype}. This causes calculation speed loss."  # noqa: E501
                )
                a_c.append(np.ascontiguousarray(a, dtype=astype))
            else:
                a_c.append(a)
        return a_c


def next_skip_comment(it, comment=("//", "return")):
    while True:
        line = next(it)
        if not line.startswith(comment):
            line = line.rsplit(comment[0])[0].strip()
            if line:
                break
    return line


def generate_surface_termination_cells(
    surface_supercell,
    film_cycle,
    name_template=None,
):
    """Generate one complete surface-slab cell per Film termination.

    The returned cells are suitable for any roughness model that selects a
    surface structure from the primitive Film stacking cycle.  Each variant is
    first shifted with ``UnitCell.affine_layer_transform`` so the requested
    primitive-cycle member is the uppermost internal layer.  All atoms are
    then assigned that one termination identifier with
    ``UnitCell.as_surface_termination``.  Internal slab planes remain encoded
    by their fractional z coordinates and may be fitted independently.

    :param UnitCell surface_supercell:
        Surface slab produced by ``UnitCell.supercell``.  Its number of
        structural layers must be an integer multiple of the Film-cycle
        length.
    :param film_cycle:
        A Film, UnitCell, LayerCycle, or iterable of primitive Film layer
        identifiers, ordered from bottom to top.
    :param str name_template:
        Optional ``str.format`` template for generated names.  The fields
        ``name`` (the source-cell name) and ``layer`` are available.  The
        default is ``"{name}_termination_<layer>"``.
    :returns:
        Mapping from primitive Film layer identifier to an independent,
        complete termination unit cell.
    :rtype: dict
    :raises ValueError:
        If either cycle is empty, contains duplicate identifiers, or the slab
        layer count is not an integer multiple of the Film cycle.
    :raises TypeError:
        If the supplied surface cell does not provide the required layered
        transformation methods.
    """
    required_methods = ("affine_layer_transform", "as_surface_termination")
    if any(not hasattr(surface_supercell, method) for method in required_methods):
        raise TypeError(
            "surface_supercell must be a layered UnitCell with affine and "
            "surface-termination transforms"
        )

    cycle = getattr(film_cycle, "layer_cycle", film_cycle)
    cycle = getattr(cycle, "layers", cycle)
    film_layers = tuple(cycle)
    if not film_layers:
        raise ValueError("Film termination cycle cannot be empty")
    if len(set(film_layers)) != len(film_layers):
        raise ValueError("Film termination cycle identifiers must be unique")

    internal_layers = tuple(surface_supercell.layer_cycle.layers)
    if not internal_layers:
        raise ValueError("Surface supercell layer cycle cannot be empty")
    if len(internal_layers) % len(film_layers):
        raise ValueError(
            "Surface slab layer count must be an integer multiple of the "
            "Film layer-cycle length"
        )

    terminations = {}
    for termination_index, termination in enumerate(film_layers):
        source_index = list(
            range(termination_index, len(internal_layers), len(film_layers))
        )[-1]
        shift = len(internal_layers) - 1 - source_index
        if name_template is None:
            layer_label = f"{termination:g}" if isinstance(
                termination, int | float | np.integer | np.floating
            ) else str(termination)
            name = f"{surface_supercell.name}_termination_{layer_label}"
        else:
            name = name_template.format(
                name=surface_supercell.name,
                layer=termination,
            )
        transformed = surface_supercell.affine_layer_transform(
            [0, 0, shift],
            name=name,
        )
        origin = max(transformed.layerpos.values())
        terminations[termination] = transformed.as_surface_termination(
            termination,
            name=name,
            origin=origin,
        )
    return terminations


SQRT2pi = np.sqrt(2 * np.pi)


def _eDensity_old(z, element, deltaZ, Qpara, deltaPara, E):
    atom_parameters = xraydb.read_Waasmaier(element)
    atom_parameters["exponents"] /= (4 * np.pi) ** 2
    dispersion = xraydb.f1_chantler(element, E) + 1j * xraydb.f2_chantler(element, E)

    rho_c = ((atom_parameters["c"] + dispersion) / (SQRT2pi * deltaZ)) * np.exp(
        -0.5 * ((deltaPara**2) * (Qpara**2) + (z**2) / (deltaZ**2))
    )

    rho_0 = np.zeros_like(z)
    for a, b in zip(atom_parameters["scale"], atom_parameters["exponents"]):
        exp_dpara = b + 0.5 * (deltaPara**2)
        exp_dz = b + 0.5 * (deltaZ**2)
        rho_0 += (a / (np.sqrt(4.0 * np.pi * exp_dz))) * np.exp(
            -exp_dpara * (Qpara**2) - ((z**2) / (4 * exp_dz))
        )
        # rho_0 += (1./(np.sqrt(4.*np.pi*exp_dz))) * np.exp(- ((z**2)/(4*exp_dz))  )

    return rho_c + rho_0


def DWtoDisorder(dw):
    return np.sqrt(dw / (8 * np.pi**2))


# returns a,b,c for Q = 4pi/lambda * sin th (instead of s)
def readWaasmaier(element):
    xraydb_t = xraydb.get_xraydb()
    wtab = xraydb_t.tables["Waasmaier"]

    row = xraydb_t.query(wtab)
    if isinstance(element, int):
        row = row.filter(wtab.c.atomic_number == element).one()
    else:
        row = row.filter(wtab.c.ion == element.title()).one()
    # if len(row) > 0:
    #    row = row[0]

    c = row.offset
    a = json.loads(row.scale)
    b = np.array(json.loads(row.exponents)) / (4 * np.pi) ** 2
    xraydb_t.close()
    return np.concatenate((a, b, [c]))


# incorrect for ions!!
def readDispersion(element, E):
    return xraydb.f1_chantler(atomic_number(element), E), xraydb.f2_chantler(
        atomic_number(element), E
    )


def atomic_number(elementname):
    if elementname in special_numbers:
        return special_numbers[elementname]
    elif elementname[-1] == "+" or elementname[-1] == "-":
        return xraydb.atomic_number(elementname[:-2])
    else:
        return xraydb.atomic_number(elementname)


def estimateDispersionCompound(name, E):
    f1 = 0.0
    f2 = 0.0
    for atom, frac in xraydb.chemparse(name).items():
        f1 += xraydb.f1_chantler(atom, E) * frac
        f2 += xraydb.f2_chantler(atom, E) * frac
    return f1, f2
