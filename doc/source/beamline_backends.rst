Beamline Backends
=================

Beamline backends tell orGUI how to open a scan, read its motor/counter
metadata, and find the detector images. They are needed because beamlines often
store the same physical information under different names or in different file
layouts.

Most users only need to select an existing backend. New backend files are useful
when data from another beamline or a changed beamline file format should be
loaded into orGUI.

Using A Backend
---------------

There are two normal ways to select a backend in GUI mode.

Config File
~~~~~~~~~~~

Add a ``[backend]`` section to the orGUI config file. To load a backend Python
file:

.. code-block:: ini

   [backend]
   file = ./backend/P212_backend.py

The path is resolved relative to the config file. When the config is loaded,
orGUI loads the Python file, finds its ``Scan`` subclass, adds it to the
backend selector, selects it, and disables auto-detection.

Alternatively, select a backend that is already registered in orGUI:

.. code-block:: ini

   [backend]
   beamtime = id31_default

GUI Selection
~~~~~~~~~~~~~

The same selection can be made in the ``Scan data`` dock described in
:doc:`gui_overview`.

1. Open the ``Scan data`` window.
2. Disable ``auto detect`` in the lower backend controls.
3. Choose an existing backend from the ``Backend`` drop-down, or click ``...``
   to load a backend Python file.
4. Select the scan file and scan number.
5. Click the open-scan button or double click a compatible scan in the file
   tree.

Auto-detection is convenient when the file metadata can identify a known
beamtime. Disable it when using a custom backend file or when the automatic
choice is wrong for the current data.

What A Backend Does
-------------------

A backend normalizes beamline-specific data into the conventions expected by
orGUI. The examples in ``examples/backend`` show the usual structure:

* read scan metadata from a Nexus/HDF5/SPEC-like file
* map beamline motor names to orGUI motor attributes
* identify the scan axis
* locate the detector image files
* return one image at a time as a numpy array
* expose counters that should be copied into integrated output files

The base class is ``orgui.backend.scans.Scan``. A backend file loaded through
the GUI must define exactly one subclass of ``Scan``.

Required Scan Attributes
------------------------

The backend constructor should populate these attributes:

``axisname``
   Name of the scanned axis. The most fully supported scan axes are ``"th"``
   and ``"mu"``.

``axis``
   Numpy array with the scan-axis values. Values are user-facing motor values
   in degrees.

``th``
   Theta motor values in degrees.

``omega``
   Internal omega-convention values in degrees. In the current orGUI
   diffractometer convention, this is usually ``omega = -th``.

``mu``
   Incidence-angle motor values in degrees.

``title``
   Human-readable scan title.

``name``
   Optional unique scan identifier used as a key in output files.

Required Scan Methods
---------------------

``__init__(hdffilepath_orNode=None, scanno=None)``
   Read the selected file or already-open HDF5 node and populate scan metadata.
   The examples handle both a file path and a node selected in the GUI file
   browser.

``parse_h5_node(cls, node)``
   Return a dictionary with at least ``scanno`` and ``name`` for a scan node
   double clicked in the Nexus tree.

``__len__()``
   Return the number of images or points in the scan.

``get_raw_img(i)``
   Return an image object for image number ``i``. The returned object must have
   an ``img`` attribute containing the image as a numpy array.

``__getitem__(i)``
   Optional but convenient. The examples implement this by returning
   ``get_raw_img(i)``.

Image Objects
-------------

The minimal image contract is small:

.. code-block:: python

   class MyImage:
       def __init__(self, filename):
           with fabio.open(filename) as img:
               self.img = img.data.astype(np.float64)
               self.header = img.header
           self.motors = {}
           self.counters = {}

Only ``img`` is required for image display and integration. ``motors`` and
``counters`` are optional dictionaries for additional metadata.

Auxiliary Counters
------------------

Backends can expose counters that should be copied into the integration output
database by overriding the ``auxillary_counters`` property. The spelling is
intentional and must match the API name:

.. code-block:: python

   @property
   def auxillary_counters(self):
       return ["exposure_time", "time", "diode"]

After integration, orGUI looks for attributes with these names on the ``Scan``
object and copies matching values into the database.

Minimal Backend Skeleton
------------------------

This skeleton shows the shape of a custom backend. Real backends usually need
beamline-specific path handling and motor-name conversion, as shown in
``examples/backend/P212_backend.py`` and ``examples/backend/CHESS_QM2.py``.

.. code-block:: python

   import os

   import fabio
   import numpy as np
   import silx.io

   from orgui.backend.scans import Scan


   class ExampleImage:
       def __init__(self, filename):
           with fabio.open(filename) as img:
               self.img = img.data.astype(np.float64)
               self.header = img.header
           self.motors = {}
           self.counters = {}


   class ExampleBackend(Scan):
       def __init__(self, hdffilepath_orNode=None, scanno=None):
           if hdffilepath_orNode is None:
               return

           if isinstance(hdffilepath_orNode, str):
               filepath = os.path.abspath(hdffilepath_orNode)
               h5file = silx.io.open(filepath)
               close_file = True
           else:
               filepath = hdffilepath_orNode.local_filename
               h5file = hdffilepath_orNode.file
               close_file = False

           try:
               scan_name = str(scanno)
               scan = h5file[scan_name]

               self.name = scan_name
               self.title = scan["title"][()]

               # Replace these paths and motor names with beamline-specific
               # metadata parsing.
               measurement = scan["measurement"]
               self.th = np.asarray(measurement["th"], dtype=np.float64)
               self.omega = -self.th
               self.mu = np.asarray(measurement["mu"], dtype=np.float64)

               self.axisname = "th"
               self.axis = self.th
               self.nopoints = len(self.axis)

               image_folder = os.path.join(os.path.dirname(filepath), "images")
               self.image_paths = [
                   os.path.join(image_folder, f"image_{i:05d}.cbf")
                   for i in range(self.nopoints)
               ]
           finally:
               if close_file:
                   h5file.close()

       @classmethod
       def parse_h5_node(cls, node):
           return {
               "scanno": int(node.basename),
               "name": node.local_name,
           }

       def __len__(self):
           return self.nopoints

       def get_raw_img(self, i):
           return ExampleImage(self.image_paths[i])

       def __getitem__(self, i):
           return self.get_raw_img(i)

Example Backend Patterns
------------------------

``examples/backend/P212_backend.py``
   Reads a P21.2-style HDF5 scan, finds available detector cameras, maps
   beamline motors to six-circle equivalents, derives ``th`` and ``mu``, and
   constructs detector image paths from detector metadata.

``examples/backend/CHESS_QM2.py``
   Reads CHESS QM2-style HDF5 metadata, maps scan motors to orGUI motor
   attributes, searches beamline image folders for matching Pilatus CBF images,
   and parses useful values from the Pilatus header.

Both examples illustrate the same principle: keep beamline-specific naming and
file-layout logic inside the backend, and expose a normalized ``Scan`` object to
the rest of orGUI.

Implementation API
------------------

For the formal base classes and method signatures, see
:ref:`backend-api`.
