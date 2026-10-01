"""Bounded HDF5 views and assembly of reduced rocking results."""

import numpy as np


class _RowDataset:
    def __init__(self, dataset, rows):
        self.dataset = dataset
        self.rows = rows
        self.shape = (rows.stop - rows.start,) + dataset.shape[1:]
        self.ndim = dataset.ndim
        self.attrs = dataset.attrs

    def __getitem__(self, key):
        # Reading a tile first also supports the old [:, 0][()] idiom.
        # This allocation is bounded by the chosen curve count, never R*F.
        return self.dataset[self.rows][key]


class _RowGroup:
    def __init__(self, group, rows):
        self.group = group
        self.rows = rows
        self.attrs = group.attrs

    def __contains__(self, key):
        return key in self.group

    def __getitem__(self, key):
        dataset = self.group[key]
        return _RowDataset(dataset, self.rows) if dataset.ndim else dataset


def _read_roi_tile(group, rows):
    result = {}
    for name, item in group.items():
        if hasattr(item, "items"):
            result[name] = _read_roi_tile(item, rows)
        else:
            result[name] = item[rows] if item.ndim else item[()]
    return result


def _join_reduced_tiles(tiles):
    """Join nested scalar-per-curve outputs, retaining NeXus attributes."""
    first = tiles[0]
    if isinstance(first, dict):
        return {
            key: _join_reduced_tiles([tile[key] for tile in tiles]) for key in first
        }
    if isinstance(first, np.ndarray):
        return np.concatenate(tiles, axis=0)
    return first
