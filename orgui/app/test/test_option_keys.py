"""Regression tests for the option-key deprecation shim.

The integration and region-of-interest option dictionaries are part of the
scripting surface -- ``startup_setup.py``-style batch scripts write them by
hand -- so renaming their keys has to keep the old spellings working for a
deprecation cycle. See :mod:`orgui.app._option_keys`.
"""

import logging
import warnings

import pytest

from orgui.app import _option_keys
from orgui.app._option_keys import (
    LEGACY_KEYS,
    LegacyKeyDict,
    canonical_key,
    canonical_options,
)


@pytest.fixture(autouse=True)
def _forget_reported():
    """Each test sees a fresh once-per-key log state."""
    _option_keys._reported.clear()
    yield
    _option_keys._reported.clear()


def test_a_current_key_passes_through_silently():
    """The common path must not warn, or every run would be noisy."""
    with warnings.catch_warnings():
        warnings.simplefilter("error", DeprecationWarning)
        assert canonical_key("solid_angle") == "solid_angle"
        assert canonical_options({"mask": True}) == {"mask": True}


def test_a_legacy_key_is_translated_and_warns():
    """Old scripts keep working, loudly enough to be noticed."""
    with pytest.deprecated_call():
        assert canonical_key("solidAngle") == "solid_angle"


def test_every_legacy_key_maps_to_a_distinct_current_name():
    """A mapping collision would silently merge two options."""
    assert len(set(LEGACY_KEYS.values())) == len(LEGACY_KEYS)
    # No legacy key may also be a current name, or translation would loop.
    assert not (set(LEGACY_KEYS) & set(LEGACY_KEYS.values()))


def test_an_unknown_key_is_left_alone():
    """Options this version does not know stay ignored, not fatal.

    ``set_integration_options`` has always skipped keys it has no branch for,
    which is what lets an older orGUI read a newer configuration.
    """
    assert canonical_key("something_new") == "something_new"
    assert canonical_options({"something_new": 1}) == {"something_new": 1}


def test_the_whole_mapping_is_rewritten_on_the_way_in():
    """Setters consume mappings by iteration, which no alias can intercept."""
    with pytest.deprecated_call():
        got = canonical_options(
            {"solidAngle": True, "sizeX": 5e-4, "mask": False}
        )
    assert got == {"solid_angle": True, "sample_size_x": 5e-4, "mask": False}
    assert list(got) == ["solid_angle", "sample_size_x", "mask"]


def test_both_spellings_at_once_is_an_error():
    """Silently preferring one would hide a real inconsistency."""
    with pytest.raises(ValueError, match="given twice"):
        with pytest.deprecated_call():
            canonical_options({"solidAngle": True, "solid_angle": False})


def test_reading_a_legacy_key_off_a_returned_dict_still_works():
    """Existing consumers index the getter's result directly."""
    options = LegacyKeyDict({"solid_angle": True, "sample_size_x": 5e-4})

    with pytest.deprecated_call():
        assert options["solidAngle"] is True
    with pytest.deprecated_call():
        assert options.get("sizeX") == 5e-4
    with pytest.deprecated_call():
        assert "solidAngle" in options
    assert options["solid_angle"] is True


def test_iterating_a_returned_dict_exposes_only_current_names():
    """A whole-dictionary round trip must not carry the old names along."""
    options = LegacyKeyDict({"solid_angle": True})

    assert list(options) == ["solid_angle"]
    assert list(options.keys()) == ["solid_angle"]
    assert dict(**options) == {"solid_angle": True}


def test_the_log_warns_once_per_key(caplog):
    """DeprecationWarning is invisible by default; the log is not.

    A scan loop would otherwise emit one record per frame.
    """
    with caplog.at_level(logging.WARNING, logger="orgui.app._option_keys"):
        with pytest.deprecated_call():
            canonical_key("solidAngle")
            canonical_key("solidAngle")
            canonical_key("sizeX")

    records = [r for r in caplog.records if "deprecated" in r.message]
    assert len(records) == 2, [r.message for r in records]
