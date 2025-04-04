*********
Changelog
*********

1.2.0 (2024-12-03)
##################

Release that adds 3 major improvements: 

#. Rocking scan integration is now properly supported with fast calculation speed and visulalization of the ROIs
#. Backends can now be loaded as standalone python file. No reinstalling of orGUI is required to add new backends.
#. hdf5 file handling has been improved to enable reloading of changed hdf5 files 

In addition, there are many more smaller improvements and changes. See the detailed changelog. 


Added
=====

* New backend style. Backends are now a standalone python file, which include a class that inherits scans.Scan. Only changing this single file is sufficient to add new backends (080642ca0d81e362dbff0fd1ff029df872cac6df, 9f1bfc83c14ffdbff9f51b5cd529059fddbfbb7e)
* add DESY P21.2 backend (080642ca0d81e362dbff0fd1ff029df872cac6df)
* add CHESS QM2 backend (1a21fbb19128a1924d8337b3d2a20ea262a56b29)  
* enable loading backends from config file (5c9f9c80ba0da4204c099aad8373d6450c8f95ef)
* U matrix editing: add manual Euler rotation adjustment (570a83ad34215d9a4e265388251a1810622890fd)
* Add batch processing example (001313244543b8988b9180813089dacf67ee506b)
* anglesToHKL is now fully numpy parallel (09776674ea491e49db76c5ce5c825528d3907971)
* NewROI class RectangleBgROI, that handles center and background ROIs (5a5efeb3d7d53ba7849fbe6d945e8efc33fb7b2c)
* Visulalization of rocking scan ROIs (829e328edc195ec9d93fc9e10ab21c58a25ea4c5)
* add PyOpenGL to optional dependencies (d3eeb1885f052ff94bcffe1915d4a8659ff99b31)

Changed
=======
* add new rocking scan integration method, that avoids duplicated file reads and leads to massive speedup. (77600563a309618fd57985c34765958e4d66b6ef)
* rocking scan integration now accessible in ROI integration panel (a629c91463b70ae0be7714039330105fd30a3d76)
* Reworked refreshing of hdf5 files (0451171a9dd2cdb9713e8a46fe4bc42494f320c5)
* HDF5_USE_FILE_LOCKING can now explicitly be set by user (a5f836050dbb31c1e1cb023a5ddc1530d5dc8d2f) 


1.1.2 (2024-11-03)
##################

Bugfix version that fixes a bug that prevents installing orGUI on certain systems. Full install target advertised on pypi is now available. More minor bug fixes, see Fixed list.

Added
=====

* pyproject.toml - add full target (37df17dc815700f0c4893d7df5f5921f5c9bb90f)

Fixed
=====

* pyproject.toml - change directory name of version path to fix a bug blocking installation on case sensitive UNIX systems (f06fe9ce36d25f51c1b837ec5db3ff17614df45e)
* Qt - replace deprecated .desktop() API for screen size determination (8ae7c187de44b0f3704a06a57968bb1213776348)
* Changes to ensure compatibility with the new Kiel diffraction setup (d532920b94aad0b00706fe20a319079a5381d51e)
* Improved error handling of 'Generate scan from images' function (13a7615ecb9f80421a7ef34d934a9e437e01f05e, c8e530c5eca27d6a193ee472c45d67a2068ee828)
* legacy id31 backend: Add 'potv' counter for fast potential measurement to backend (ead0c8061f029aec9d8565bdb3099dd59ccfc49e)
* plot3d, pos_cart_all: correct lookup of atom radius, bugfix in translate argument (88ddb0e40449ebe6f1127091eb419bbac1dfd73f, 888cdd92e51ff9c69c995dc10e72f1131cfe527b)

1.1.1 (2024-08-06)
##################

Version that fixes a critical bug in the config file reading when a poni file is given.

Fixed
=====

* Energy set from poni file if provided in config file (434fae447351f626c280ab3f5c4b94dd4124b813)

1.1.0 (2024-08-05)
##################

This is a release, that reworks the configuration widgets. 
It is aimed to correctly show all configuration options in the GUI, which were previously only available through the config file.

orGUI now also has a proper logo.

Added
=====

* Config files in example folder to illustrate the multitude of config options
* logo and application icon
* This CHANGELOG
* Enable optional loading of atomic coordinates and lattice from `ASE <https://wiki.fysik.dtu.dk/ase/>`_ supported files. These include:
  
  * Cif files
  * VASP / CONTCAR files
  * many more ... see `files supported by ASE <https://wiki.fysik.dtu.dk/ase/ase/io/io.html>`_

Changed
=======

* Rework of machine parameter widget: (4bc83799e082c7ef8178e017b24c373dc0bec97c)

  * Uses pyFAI DetectorSelectorDialog for detector selection 
  * Uses pyFAI geometry widget to display detector geometry as in pyFAI or Fit2D.
  * Add button to directly load poni file.
  
* Crystal parameters widget:

  * Add button to link lattice parameters upon editing. (51f21be749ab7d16d7ae2ae883b5edd40489ecda)
  * Add button to show the coordinates of atoms in the unit cell
  
* Detector size can now only be changed explicitly in machine parameter widget. (24bdffa41e2c61f3f57dc61cb2ccd11298daa4f7) 
  
* Splashscreen: add version number (c36c05ea421239da28bf49d4adff7f0654f40f90)

Fixed
=====

* Removed unnecessary config warning upon startup, if no config file was provided with the start. (c36c05ea421239da28bf49d4adff7f0654f40f90)
* Fix QLayout double assignment warning. (8950c9edd8841891f0f9ba15999345865d3419f3)
* Explicitly close file handles in universal scan loader after data was read (0a9f5d8e8828bdc0523a2cd9673d04b45c634bd7) 
* loading of xtal and bulk files from relative file path
  

1.0.1 (2024-07-29)
##################

First release version on PyPi, project availabe under `orGUI <https://pypi.org/project/orGUI/>`_.
This version aimes to fix some minor issues from the first release and provides a minimal documentation in form of a README and a diffractometer image showing the geometry.

Added
=====

* Add README
* Add diffractometer image in help menu (bb4fd6724f7ec040cea1a764afa79cce5afb32cc)

Changed
=======

* Enable start of *orGUI* without config file. 
* Replace old `setup.py` build system by `pyproject.toml`.

Fixed
=====

* Fix license statement in about dialog to show the `MIT License` instead of `All rights reserved`. (bb4fd6724f7ec040cea1a764afa79cce5afb32cc)
* Crash if hkl reflection is searched for static roi, but no scan is selected (6f6546ddeb9dd2a28ab11963f9d6d7bd694a78c8) 
* Enable loading of poni file with relative file path in config file. (1cdae8219ab3124d28231bb2282e2475729ffec8)
* Rename deprecated ``np.product`` -> ``np.prod`` (38348f11dc0787413b5b1513b7bc9bcf93fda22e)


1.0.0 (2024-07-02)
##################

* First public release on Zenodo with the DOI `https://doi.org/10.5281/zenodo.12592486 <https://doi.org/10.5281/zenodo.12592486>`_.
* Added MIT License
* Merged with `datautils`

  * Backends

    * datautils.xrayutils.id31_tools -> orgui.backend.beamline.id31_tools
    * datautils.xrayutils.P212_tools -> orgui.backend.beamline.P212_tools
    * datautils.xrayutils.ID31DiffractLinTilt -> orgui.backend.beamline.ID31DiffractLinTilt
    * datautils.xrayutils.fio_reader -> orgui.backend.beamline.fio_reader

  * Physics / Diffraction caculations

    * datautils.util -> orgui.datautils.util
    * datautils.xrayutils.unitcells -> orgui.datautils.xrayutils.unitcells (and .bul files therein)
    * datautils.xrayutils.test -> orgui.datautils.xrayutils.test (and datautils test code therein)
    * datautils.xrayutils.CTRcalc -> orgui.datautils.xrayutils.CTRcalc
    * datautils.xrayutils._CTRcalc_accel -> orgui.datautils.xrayutils._CTRcalc_accel
    * datautils.xrayutils.CTRopt -> orgui.datautils.xrayutils.CTRopt
    * datautils.xrayutils.CTRplotutil -> orgui.datautils.xrayutils.CTRplotutil
    * datautils.xrayutils.DetectorCalibration -> orgui.datautils.xrayutils.DetectorCalibration
    * datautils.xrayutils.element_data -> orgui.datautils.xrayutils.element_data
    * datautils.xrayutils.HKLVlieg -> orgui.datautils.xrayutils.HKLVlieg
    * datautils.xrayutils.ReciprocalNavigation -> orgui.datautils.xrayutils.ReciprocalNavigation

* *orGUI* is now a standalone package, that has only publicly available dependencies.

0.9-alpha (2024-06-29)
######################

* Last inoffical version before merge with the internal dependency `datautils`.
* Last version under `All rights reserved`.