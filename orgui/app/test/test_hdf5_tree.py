"""Regressions for browsing malformed SPEC measurement groups."""

import sys
from unittest.mock import Mock

import h5py
import numpy as np
import pytest
from silx.gui import qt
from silx.gui.hdf5 import Hdf5TreeView, NexusSortFilterProxyModel

from orgui.app.qutils import RobustHdf5TreeModel


@pytest.fixture(scope="session")
def qapp():
    """Keep a Qt application alive for the tree views."""
    return qt.QApplication.instance() or qt.QApplication([])


def _child(model, parent, name):
    for row in range(model.rowCount(parent)):
        index = model.index(row, 0, parent)
        if model.nodeFromIndex(index).basename == name:
            return index
    raise AssertionError(f"Missing tree child: {name!r}")


@pytest.mark.parametrize("extra_column", [False, True])
def test_expand_spec_measurement_with_empty_label(
    qapp, tmp_path, monkeypatch, extra_column
):
    """A trailing blank label must not crash painting or sorting the tree."""
    path = tmp_path / "scan.spec"
    content = (
        "#F scan.spec\n#E 1455180875\n#D Thu Feb 11 09:54:35 2016\n"
        "#C regression fixture\n#O0 phi\n\n"
        "#S 1 ascan phi 0 1 1 1\n#D Thu Feb 11 09:55:20 2016\n"
        f"#P0 0\n#N {3 if extra_column else 2}\n#L phi  ic1  \n"
        + ("0 10 100\n1 20 200\n\n" if extra_column else "0 10\n1 20\n\n")
    )
    path.write_text(content, encoding="utf-8", newline="\n")
    original = path.read_bytes()
    view = Hdf5TreeView()
    model = RobustHdf5TreeModel(view, ownFiles=True)
    proxy = NexusSortFilterProxyModel(view)
    proxy.setSourceModel(model)
    view.setModel(proxy)
    view.setSortingEnabled(True)
    errors = Mock()
    monkeypatch.setattr(sys, "excepthook", errors)
    model.insertFile(str(path))
    root = model.index(0, 0)
    try:
        scan = _child(model, root, "1.1")
        measurement = _child(model, scan, "measurement")
        unnamed = _child(model, measurement, "")
        assert model.rowCount(measurement) == 3
        assert model.data(unnamed) == "Unnamed SPEC column"
        assert model.flags(unnamed) == qt.Qt.NoItemFlags
        assert model.data(unnamed, model.H5PY_OBJECT_ROLE) is None
        assert "empty label" in model.data(unnamed, qt.Qt.ToolTipRole)

        for column in range(model.columnCount()):
            index = model.index(unnamed.row(), column, measurement)
            proxy_index = proxy.mapFromSource(index)
            for role in (
                qt.Qt.DisplayRole, qt.Qt.DecorationRole, qt.Qt.ToolTipRole,
                model.H5PY_ITEM_ROLE, model.H5PY_OBJECT_ROLE,
            ):
                proxy.data(proxy_index, role)

        for order in (qt.Qt.AscendingOrder, qt.Qt.DescendingOrder):
            proxy.sort(0, order)
            for index in (root, scan, measurement):
                view.expand(proxy.mapFromSource(index))
            view.show()
            qapp.processEvents()
            view.collapse(proxy.mapFromSource(measurement))
            qapp.processEvents()
        errors.assert_not_called()

        counter = _child(model, measurement, "ic1")
        assert model.flags(counter) & qt.Qt.ItemIsEnabled
        np.testing.assert_array_equal(
            model.data(counter, model.H5PY_OBJECT_ROLE)[()], [10, 20]
        )
        # Retain even the inaccessible key in the original adapter; this is
        # only a GUI workaround, not a rewrite of measurements or labels.
        group = model.data(measurement, model.H5PY_OBJECT_ROLE)
        assert list(group.keys()) == ["phi", "ic1", ""]
        assert path.read_bytes() == original
    finally:
        view.close()
        model.removeIndex(root)
        view.deleteLater()


def test_named_hdf5_dataset_keeps_normal_tree_behavior(qapp, tmp_path):
    """Ordinary HDF5 counters retain their names, values and selectable flags."""
    path = tmp_path / "scan.h5"
    with h5py.File(path, "w") as file:
        file["measurement/ic1"] = [10, 20]
    model = RobustHdf5TreeModel(ownFiles=True)
    model.insertFile(str(path))
    root = model.index(0, 0)
    try:
        measurement = _child(model, root, "measurement")
        counter = _child(model, measurement, "ic1")
        assert model.data(counter) == "ic1"
        assert model.flags(counter) & qt.Qt.ItemIsSelectable
        assert model.data(counter, qt.Qt.DecorationRole) is not None
        np.testing.assert_array_equal(
            model.data(counter, model.H5PY_OBJECT_ROLE)[()], [10, 20]
        )
    finally:
        model.removeIndex(root)
