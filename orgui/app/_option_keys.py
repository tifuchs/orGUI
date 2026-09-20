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
"""Option-dictionary key names, and the legacy spellings they replaced.

The integration options and the advanced region-of-interest options are
passed around as plain dictionaries, including by user batch scripts. Their
keys were originally spelled in mixed case (``solidAngle``,
``DetectorInclination``, ``sizeX``) while everything those options are stored
next to uses ``snake_case``. The names are now ``snake_case`` throughout, and
the old spellings keep working for one deprecation cycle.

Both directions are covered, and they need different mechanisms:

* a script that *passes* an old key is handled by :func:`canonical_options`,
  which rewrites the mapping on the way in. Rewriting on read is not enough,
  because ``for key in ddict`` and ``**ddict`` never reach ``__getitem__``.
* code that *reads* an old key off a returned dictionary is handled by
  :class:`LegacyKeyDict`, which resolves the alias on lookup.

Each legacy key warns twice on first use: a :exc:`DeprecationWarning` for
tooling, and one log record, because :exc:`DeprecationWarning` is invisible by
default and these dictionaries are mostly written in batch scripts whose only
output is a log.
"""

import logging
import warnings

logger = logging.getLogger(__name__)

__all__ = [
    "LEGACY_KEYS",
    "LegacyKeyDict",
    "canonical_key",
    "canonical_options",
]

#: Legacy option-dictionary key -> the name that replaced it.
LEGACY_KEYS = {
    # Integration options
    "solidAngle": "solid_angle",
    # Advanced region-of-interest options
    "DetectorInclination": "detector_inclination",
    "ProjectSampleSize": "project_sample_size",
    "sizeX": "sample_size_x",
    "sizeY": "sample_size_y",
    "sizeZ": "sample_size_z",
    "xoffset": "offset_x",
    "yoffset": "offset_y",
    "FittedBackground": "fitted_background",
    "FittedBackgroundOrder": "fitted_background_order",
}

#: Legacy keys already reported, so a loop over a scan warns once, not once
#: per frame.
_reported = set()


def _report(legacy):
    """Warn once per legacy key, to both the warnings system and the log."""
    message = (
        f"The integration option key {legacy!r} is deprecated; use "
        f"{LEGACY_KEYS[legacy]!r} instead. The old spelling will be removed "
        f"in a future release."
    )
    warnings.warn(message, DeprecationWarning, stacklevel=3)
    if legacy not in _reported:
        _reported.add(legacy)
        logger.warning(message)


def canonical_key(key):
    """Current spelling of one option key.

    :param str key: A current or legacy key.
    :returns: The current spelling; unknown keys are returned unchanged, so
        that a caller passing an option this version does not know about gets
        the same "ignored" behaviour as before rather than an error.
    :rtype: str
    """
    if key in LEGACY_KEYS:
        _report(key)
        return LEGACY_KEYS[key]
    return key


def canonical_options(mapping):
    """Rewrite an option mapping onto the current key spellings.

    Applied at the entry of every setter, because a mapping is consumed by
    iteration and unpacking as often as by lookup.

    :param mapping: Option dictionary, possibly using legacy keys. ``None``
        is treated as empty.
    :returns: A new plain :class:`dict` with current keys. A legacy key and
        its replacement in the same mapping is an error rather than a silent
        precedence rule.
    :rtype: dict
    :raises ValueError: If a key is given under both spellings.
    """
    result = {}
    for key, value in dict(mapping or {}).items():
        name = canonical_key(key)
        if name in result:
            raise ValueError(
                f"option {name!r} was given twice, once as {key!r}; pass only "
                f"one spelling"
            )
        result[name] = value
    return result


class LegacyKeyDict(dict):
    """Option dictionary that still answers to the legacy key spellings.

    Returned by the option getters so that existing code and scripts reading
    ``options["solidAngle"]`` keep working, with a deprecation warning. Only
    lookup is aliased: iteration, :meth:`keys` and unpacking expose the
    current names, which is what makes the old spellings disappear from
    anything that round-trips a whole dictionary.
    """

    def __getitem__(self, key):
        return super().__getitem__(canonical_key(key))

    def __contains__(self, key):
        return super().__contains__(canonical_key(key))

    def get(self, key, default=None):
        return super().get(canonical_key(key), default)

    def pop(self, key, *default):
        return super().pop(canonical_key(key), *default)
