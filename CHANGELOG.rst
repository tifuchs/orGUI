*********
Changelog
*********

1.1.0a (2024-08-03)
###################

This is a release, that reworks the configuration widgets. 
It is aimed to correctly show all configuration options in the GUI, which were previously only available through the config file.

Added
=====

* This CHANGELOG

Changed
=======

* Rework of machine parameter widget: (4bc83799e082c7ef8178e017b24c373dc0bec97c)
  * Uses pyFAI DetectorSelectorDialog for detector selection 
  * Uses pyFAI geometry widget to display detector geometry as in pyFAI or Fit2D.
  * Add button to directly load poni file.
  
* Crystal parameters widget:
  * Add button to link lattice parameters upon editing. (51f21be749ab7d16d7ae2ae883b5edd40489ecda)
  
* Detector size can now only be changed explicitly in machine parameter widget. (24bdffa41e2c61f3f57dc61cb2ccd11298daa4f7) 
  
* Splashscreen: add version number (c36c05ea421239da28bf49d4adff7f0654f40f90)

Fixed
=====

* Removed unnecessary config warning upon startup, if no config file was provided with the start. (c36c05ea421239da28bf49d4adff7f0654f40f90)
* Fix QLayout double assignment warning. (8950c9edd8841891f0f9ba15999345865d3419f3)
* Explicitly close file handles in universal scan loader after data was read (0a9f5d8e8828bdc0523a2cd9673d04b45c634bd7) 
  

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
* Merged with `datautils`:
  * Backends:
    * datautils.xrayutils.id31_tools :raw-html:`&rarr;` orgui.backend.beamline.id31_tools
    * datautils.xrayutils.P212_tools :raw-html:`&rarr;` orgui.backend.beamline.P212_tools
    * datautils.xrayutils.ID31DiffractLinTilt :raw-html:`&rarr;` orgui.backend.beamline.ID31DiffractLinTilt
    * datautils.xrayutils.fio_reader :raw-html:`&rarr;` orgui.backend.beamline.fio_reader
  * Physics / Diffraction caculations:
    * datautils.util :raw-html:`&rarr;` orgui.datautils.util
    * datautils.xrayutils.unitcells :raw-html:`&rarr;` orgui.datautils.xrayutils.unitcells (and .bul files therein)
    * datautils.xrayutils.test :raw-html:`&rarr;` orgui.datautils.xrayutils.test (and datautils test code therein)
	* datautils.xrayutils.CTRcalc :raw-html:`&rarr;` orgui.datautils.xrayutils.CTRcalc
	* datautils.xrayutils._CTRcalc_accel :raw-html:`&rarr;` orgui.datautils.xrayutils._CTRcalc_accel
	* datautils.xrayutils.CTRopt :raw-html:`&rarr;` orgui.datautils.xrayutils.CTRopt
	* datautils.xrayutils.CTRplotutil :raw-html:`&rarr;` orgui.datautils.xrayutils.CTRplotutil
	* datautils.xrayutils.DetectorCalibration :raw-html:`&rarr;` orgui.datautils.xrayutils.DetectorCalibration
	* datautils.xrayutils.element_data :raw-html:`&rarr;` orgui.datautils.xrayutils.element_data
	* datautils.xrayutils.HKLVlieg :raw-html:`&rarr;` orgui.datautils.xrayutils.HKLVlieg
	* datautils.xrayutils.ReciprocalNavigation :raw-html:`&rarr;` orgui.datautils.xrayutils.ReciprocalNavigation
* *orGUI* is now a standalone package, that has only publicly available dependencies.

0.9-alpha (2024-06-29)
######################

* Last inoffical version before merge with the internal dependency `datautils`.
* Last version under `All rights reserved`.