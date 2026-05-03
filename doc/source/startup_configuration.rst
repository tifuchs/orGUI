Startup And Configuration
=========================

This page describes the normal graphical user interface (GUI) startup path and
the configuration file used to initialize orGUI. It is intended for users who
already know the beamline, detector calibration, sample lattice, and scan file
format they want to use.

Start orGUI with a configuration file from a terminal:

.. code-block:: bash

   orGUI examples/config_minimal

The same configuration can also be loaded from the GUI with
``Config -> Load config``.

Configuration File Overview
---------------------------

The config file is an INI-style text file. Section names are written in square
brackets, and values are written as ``key = value`` pairs. Relative file paths,
such as ``./poni_files/P3_100_calib.poni``, are resolved relative to the
directory containing the config file.

The main sections are:

``[Machine]``
   X-ray source, detector geometry, detector size, pixel size, and optionally a
   pyFAI ``.poni`` calibration file.

``[Lattice]``
   Crystal lattice and, optionally, a crystal/unit-cell file containing the
   atomic basis used for structure-factor based Bragg reflection and CTR
   calculations.

``[Diffractometer]``
   Sample and instrument angles used to connect detector pixels to reciprocal
   space.

``[backend]``
   Optional scan-loading backend selection. This can point to a backend Python
   file or select a registered beamtime/backend id.

``[Settings]``
   GUI startup behavior and output database compression.

Full Configuration Example
--------------------------

This example combines the commonly used settings from the bundled example
files. It uses a pyFAI PONI file for detector geometry, a CIF file for the
crystal unit cell, a backend file for scan loading, and explicit GUI settings.

.. code-block:: ini

   [Machine]
   # Detector calibration from pyFAI.
   # When this is set, the PONI file supplies the detector geometry and
   # wavelength used by orGUI.
   poni = ./poni_files/P3_100_calib.poni

   # Fallback detector/source values. These are useful when no PONI file is
   # given. If a PONI file is given, the loaded PONI geometry takes precedence.
   # Energy is in keV.
   E = 77.0
   # Sample-detector distance is in m.
   SDD = 0.78
   # Detector pixel size is in m.
   pixelsize = 172e-6
   # Detector size is in pixels.
   sizex = 1475
   sizey = 1679
   # Direct-beam position is in detector pixels.
   cpx = 731.0
   cpy = 1587.856

   [Lattice]
   # Crystal/unit-cell file used for structure-factor based Bragg reflection
   # and CTR calculations. CIF files require ASE support.
   crystal = ./crystal_files/Pt_mp-126_symmetrized.cif

   # Refraction correction parameter. This is usually delta in n = 1 - delta.
   # Set to 0 if no refraction correction should be applied.
   refractionindex = 1.1415e-06

   # Optional lattice override. If these values are present, they replace the
   # lattice constants read from the crystal file while preserving relative
   # atomic coordinates.
   # Lengths are in Angstrom; angles are in degrees.
   # a1 = 3.9242
   # a2 = 3.9242
   # a3 = 3.9242
   # alpha1 = 90.0
   # alpha2 = 90.0
   # alpha3 = 90.0

   [Diffractometer]
   # Rotates the diffractometer around the beam direction, in degrees.
   azimuthal_reference = 90

   # Polarization correction settings.
   polarization_axis = 0
   polarization_factor = 1

   # User-facing diffractometer angles are in degrees.
   mu = 0.1
   chi = 0
   phi = 0

   [backend]
   # Either load a backend file...
   file = ./backend/P212_backend.py

   # ...or select a registered backend/beamtime id instead.
   # beamtime = id31_default

   [Settings]
   # Automatically load all images to create max/sum images when opening a scan.
   autoload = True

   # Output database compression. Raw is the most portable option.
   # Other filters require the corresponding HDF5/plugin support.
   compression = Raw

Machine Section
---------------

The simplest detector setup gives the detector distance, pixel size, detector
shape, and direct-beam pixel position directly in the config file. These values
are enough for basic detector geometry when no PONI file is used.

If ``poni`` is set, orGUI loads detector geometry from that pyFAI file. In that
case, the PONI file provides the detector geometry and wavelength used by the
loaded detector model. The explicit ``E``, ``SDD``, ``pixelsize``, ``sizex``,
``sizey``, ``cpx``, and ``cpy`` entries remain useful documentation and
fallback values, but the PONI calibration is the active geometry.

Common units in this section:

* ``E``: keV
* ``SDD``: m
* ``pixelsize``: m
* ``sizex`` / ``sizey``: pixels
* ``cpx`` / ``cpy``: pixels

Lattice Section
---------------

The ``[Lattice]`` section defines the crystal lattice and, optionally, the
atomic basis. There are two common levels of detail:

* Lattice constants only: orGUI can show all integer Bragg positions for that
  lattice, but it does not know which reflections are forbidden by the basis.
* Crystal/unit-cell file: orGUI reads atomic coordinates and can calculate
  structure-factor based Bragg reflections and available CTRs.

Lattice constants use Angstrom and degrees:

.. code-block:: ini

   [Lattice]
   a1 = 2.7740
   a2 = 2.7740
   a3 = 6.7960
   alpha1 = 90.0
   alpha2 = 90.0
   alpha3 = 120.0
   refractionindex = 1.1415e-06

For structure-factor based calculations, add ``crystal = ...``:

.. code-block:: ini

   [Lattice]
   crystal = ./crystal_files/Pt_mp-126_symmetrized.cif
   refractionindex = 1.1415e-06

If both ``crystal`` and explicit lattice constants are present, orGUI keeps the
relative atomic coordinates from the crystal file and overrides the lattice
constants with the values in the config.

Diffractometer Section
----------------------

The ``[Diffractometer]`` section sets the instrument orientation used by the
detector-to-reciprocal-space transform. These values are user-facing GUI/config
values and are written in degrees.

.. code-block:: ini

   [Diffractometer]
   azimuthal_reference = 90
   polarization_axis = 0
   polarization_factor = 1
   mu = 0.1
   chi = 0
   phi = 0

``azimuthal_reference`` rotates the diffractometer around the beam direction.
``mu`` is the angle of incidence. ``chi`` and ``phi`` set the sample
orientation circles. The polarization settings are used when polarization
correction is enabled during integration.

Backend Section
---------------

The backend tells orGUI how to read scan metadata and image locations for a
beamline or file format. A config can either load a backend Python file:

.. code-block:: ini

   [backend]
   file = ./backend/P212_backend.py

or select a registered backend id:

.. code-block:: ini

   [backend]
   beamtime = id31_default

When ``beamtime`` is used, orGUI selects that backend in the scan-data panel and
disables backend auto-detection. If neither option is set, the GUI can still use
the backend selector and ``auto detect`` controls in the ``Scan data`` dock.

Settings Section
----------------

``autoload`` controls whether orGUI automatically loads all images to create a
maximum and summed image when a scan is opened. This is convenient for
interactive browsing, but can be slow for large scans.

``compression`` selects the compression filter for new output database files.
``Raw`` is always the safest portable choice. Other choices, such as
``Blosc-lz4-Shuffle-5``, require matching HDF5/plugin support in the runtime
environment.

Crystal File Options
--------------------

The ``crystal`` key in ``[Lattice]`` can refer to different ways of describing
the same basic information: a lattice and, when available, the atomic basis of a
crystal unit cell.

CIF Files
~~~~~~~~~

CIF support uses the ASE package. A simple reference looks like:

.. code-block:: ini

   [Lattice]
   crystal = ./crystal_files/Pt_mp-126_symmetrized.cif
   refractionindex = 1.1415e-06

CIF is a good choice when the structure is available from crystallographic
databases or other structure-building tools.

orGUI ``.bul`` Files
~~~~~~~~~~~~~~~~~~~~

``.bul`` files are orGUI-specific unit-cell files. They are also compatible
with bulk files from the surface diffraction modeling software
`ANA-ROD <https://www.esrf.fr/computing/scientific/joint_projects/ANA-ROD/index.html>`_.
The ESRF ANA-ROD page describes ANA, AVE, and ROD as surface-diffraction
programs originally written by Elias Vlieg and later maintained at ESRF.

Example config reference:

.. code-block:: ini

   [Lattice]
   crystal = ./crystal_files/Pt100.bul
   refractionindex = 1.1415e-06

Minimal ``.bul`` content:

.. code-block:: text

   // Pt(100) larger unit cell for compatibility
   return
   3.9242 3.9242 3.9242 90.0000 90.0000 90.0000
   Pt     0.00000     0.00000     -1.00000  0.4353  0.4353  1.0000
   Pt     0.50000     0.50000     -1.00000  0.4353  0.4353  1.0000

Bundled ``.bul`` files can also be referenced by name, for example:

.. code-block:: ini

   [Lattice]
   crystal = Pt100

orGUI ``.xtal`` Files
~~~~~~~~~~~~~~~~~~~~~

``.xtal`` files are orGUI-specific crystal model files. When used for GUI
orientation and reflection calculations, orGUI reads the bulk unit cell from
the file.

Example config reference:

.. code-block:: ini

   [Lattice]
   crystal = ./crystal_files/Pt100_with_surface_model.xtal
   refractionindex = 1.1415e-06

Excerpt of a bundled ``.xtal`` example:

.. code-block:: text

   E = 68.00000 keV
   # UnitCell bulk
   return
   Coherent 1.00000   1.00000 0.00000 0.00000 0.00000 1.00000 0.00000 0.00000 0.00000 1.00000
   2.7748 2.7748 3.9242 90.0000 90.0000 90.0000
   Pt     0.00000     0.00000     -1.00000  0.4350  0.4350  1.0000
   Pt     0.50000     0.50000     -0.50000  0.4350  0.4350  1.0000
