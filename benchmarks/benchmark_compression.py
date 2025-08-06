# Copyright (c) 2023 European Synchrotron Radiation Facility
# Adapted for orGUI by Timo Fuchs, Copyright (c) 2025 Timo Fuchs
# License: MIT https://opensource.org/license/mit/
from __future__ import annotations

import json
import os
import sys
import tempfile
import time
from typing import Optional, NamedTuple

import h5py
import numpy


# Set affinity and env. var. before importing hdf5plugin

if len(sys.argv) >= 2:
    NCPU = int(sys.argv[1])
    print(f"Set affinity: Using {NCPU} CPU(s)")
    os.sched_setaffinity(0, list(range(0, NCPU, 4)))
    os.environ["OPENMP_NUM_THREADS"] = str(NCPU)
    os.environ["BLOSC_NTHREADS"] = str(NCPU)
else:
    NCPU = os.cpu_count()
    print(f"Default affinity: Using {NCPU} CPU(s)")

# Directory where to run read/write benchmarks
DIRECTORY = "."

DATABASE_PATH = './database_example.h5'

from silx.io.dictdump import dicttonx ,nxtodict


import hdf5plugin


class Result(NamedTuple):
    """Store benchmark result"""

    raw_nbytes: int
    compressed_nbytes: int
    write_duration: float
    read_duration: float

    compression_rate = property(
        lambda self: self.raw_nbytes / self.compressed_nbytes)
    write_speed = property(
        lambda self: (self.raw_nbytes / 1024**2) / self.write_duration,
        doc="Unit: MB/sec")
    read_speed = property(
        lambda self: (self.raw_nbytes / 1024**2) / self.read_duration,
        doc="Unit: MB/sec")

    
def benchmark(
    data: dict,
    data_nbytes : int,
    directory: Optional[str] = None,
    **kwargs
) -> Result:
    """Run benchmark for given conditions

    :param data: Dataset to use
    :param directory: Directory where to write HDF5 file for benchmark
    :param kwargs: Arguments passed to h5py.Group.create_dataset
    """
    if directory is None:
        tmpdir = tempfile.TemporaryDirectory()
        dirname = tmpdir.name
    else:
        dirname = directory

    filename = os.path.join(dirname, 'hdf5plugin_benchmark.h5')

    if os.path.exists(filename):
        os.remove(filename)

    # Compression
    start_write_time = time.perf_counter()
    dicttonx(data, filename, create_dataset_args=kwargs)
        
        
    #with h5py.File(filename, "w") as h5file:
    #    dataset = h5file.create_dataset(
    #        "data", shape=data.shape, dtype=data.dtype, **kwargs)
    #    
    #    dataset[:] = data
    #    dataset.flush()
    write_duration = time.perf_counter() - start_write_time

    # Decompression
    start_time = time.perf_counter()
    alldata = nxtodict(filename)
    #with h5py.File(filename, "r") as h5file:
    #    dataset = h5file["data"]
    #    start_time = time.perf_counter()
    #    read_data = dataset[()]
    #    read_duration = time.perf_counter() - start_time
    #    storage_size = dataset.id.get_storage_size()
    #    chunks = dataset.chunks
    #
    read_duration = time.perf_counter() - start_time
    storage_size = os.path.getsize(filename)
    #storage_size = dataset.id.get_storage_size()
    os.remove(filename)

    return Result(data_nbytes, storage_size, write_duration, read_duration)

    
DEFAULT_FILTERS = {  # Filters available with h5py/libhdf5
    "Raw": None,
    "GZip": "gzip",
    "LZF": "lzf",
}

LOSSLESS_FILTERS = {
    "BZip2": hdf5plugin.BZip2(),
    "LZ4": hdf5plugin.LZ4(),
    "ZStd": hdf5plugin.Zstd(),
}

BITSHUFFLE_FILTERS = {
    "Bitshuffle-lz4": hdf5plugin.Bitshuffle(cname='lz4'),
    "Bitshuffle-zstd": hdf5plugin.Bitshuffle(cname='zstd'),
}
    
BLOSC_FILTERS = {}
for cname in ('lz4', 'blosclz', 'lz4', 'lz4hc', 'snappy', 'zlib', 'zstd'):
    for shuffle_name, shuffle in [('NoShuffle', hdf5plugin.Blosc.NOSHUFFLE),
                                  ('Shuffle', hdf5plugin.Blosc.SHUFFLE),
                                  ('BitShuffle', hdf5plugin.Blosc.BITSHUFFLE)]:
        for clevel in [5]: #(1, 3, 5, 9):
            BLOSC_FILTERS[f"Blosc-{cname}-{shuffle_name}-{clevel}"] = hdf5plugin.Blosc(
                cname=cname, clevel=clevel, shuffle=shuffle)


BLOSC2_FILTERS = {}
for cname in ('lz4', 'blosclz', 'lz4', 'lz4hc', 'zlib', 'zstd'):
    for filters_name, filters in [('NoFilter', hdf5plugin.Blosc2.NOFILTER),
                                ('Shuffle', hdf5plugin.Blosc2.SHUFFLE),
                                ('BitShuffle', hdf5plugin.Blosc2.BITSHUFFLE)]:
        for clevel in [5]: # (1, 3, 5, 9):
            BLOSC2_FILTERS[f"Blosc2-{cname}-{filters_name}-{clevel}"] = hdf5plugin.Blosc2(
                cname=cname, clevel=clevel, filters=filters)


FILTERS = {
    **DEFAULT_FILTERS,
    **LOSSLESS_FILTERS,
    **BITSHUFFLE_FILTERS,
    **BLOSC_FILTERS,
    **BLOSC2_FILTERS,
}


print("Read benchmark data...")
data_size = os.path.getsize('database_example.h5')
data = nxtodict('database_example.h5')
print("Run benchmarks:")

results = {}
for name, compression in FILTERS.items():
    print(f"- {name}")
    results[name] = benchmark(
        data,
        data_size,
        directory=DIRECTORY,
        compression=compression,
    )._asdict()


print("Write results")
with open("benchmark.json", "w") as f:
    json.dump(
        dict(
            config=dict(
                ncpu=NCPU,
                affinity=os.cpu_count(),
                openmp_num_threads=os.environ.get("OPENMP_NUM_THREADS", "unset"),
                blosc_nthreads=os.environ.get("BLOSC_NTHREADS", "unset"),
                hdf5plugin=hdf5plugin.get_config().build_config._asdict()
            ),
            results=results,
        ),
        f,
        indent=2,
        default=lambda obj: type(obj).__name__,
    )
