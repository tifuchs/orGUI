# Changelog

This is the changelog for the software orGUI, written by Timo Fuchs


## [1.4.0] (2026-01-20)

[ef10ee8](https://github.com/tifuchs/orGUI/commit/ef10ee8e805c13b60bb51361c38868da928d7116)...[eb2d82f](https://github.com/tifuchs/orGUI/commit/eb2d82f796a5f34798d7f78f7794ef73b9231332)

This is a feature release with lots of additional bug fixes.

- Added a new dialog which automatically calculates the position of the
  Bragg peak, which helps to improve the U matrix the most and adds it
  to the list of reference reflections.
- Added a Bragg peak search, which searches for the center of mass in a
  small section of multiple images. This can be used to quickly find the
  image and position of a Bragg peak. To be used when an estimate of the
  Bragg position is already given (such as with the new Bragg peak
  dialog)
- Auxilliary counters (potential, temperature, \...) are now also
  integrated in all integration functions.
- A static background image can now be subtracted from the detector
  images. This can be used to subtract static sample environment
  scattering such as water or air (beamline background) scattering.
- Multiple scans can be concatenated and treated as a single scan. This
  is done to avoid high-intensity sample positions which would lead to
  damage to the detector (Thanks to Finn!).
- Positions of Bragg reflections can now be integrated as rocking scans.
  This enables fast conventional single crystal X-ray difffraction
  analysis of Bragg reflection structure factor. Strain with respect to
  the reference lattice can be added for epitaxially grown thin films.
- 2D ROI size can now be dynamically adjusted. Currently, the ROI size
  can be corected for
  - Peak size increase due to scattered beam - detector inclination
  - Peak size increase due to the illuminated sample size. Currently,
    the sample is approximated as a cuboid (size to be set in Advanced
    ROI options dialog) and its outlines are projected onto the
    detector.
- The orGUI hdf5 data base can now be compressed using
  [hdf5plugin](https://github.com/silx-kit/hdf5plugin). This compression
  is not visible to the user, but reduces the file size by a factor of
  2x - 3x. Calculation speed loss is possible and depends on the OS,
  available libraries , etc. See [Issue
  #16](https://github.com/tifuchs/orGUI/issues/16) for an overview.
  Opening compressed data bases with external software requires the user
  to install and import
  [hdf5plugin](https://github.com/silx-kit/hdf5plugin).

### Added

- [ReflectionSelector] add automatic peak search (2d image peak and rocking peak) and automatic Bragg reflection add function ([985f9ca](https://github.com/tifuchs/orGUI/commit/985f9ca803eead6487d68a3341c6652ff39b72cc))
- [ReflectionSelector] add automatic Bragg reflection search dialog ([4fe31db](https://github.com/tifuchs/orGUI/commit/4fe31dbc1a181b9d37561dd16bb488c1b8c6253e))
- [Database] add optional compression, add compression benchmarks ([e468efd](https://github.com/tifuchs/orGUI/commit/e468efdbfcba6bd13a2791674388f84274999a39))
- [RoIntegr] add aux counter integration ([2aa1cc9](https://github.com/tifuchs/orGUI/commit/2aa1cc9ee9e9a4b580562c0df6dd83ece1a4c40a))
- Add static background image subtraction ([b9777f9](https://github.com/tifuchs/orGUI/commit/b9777f90aed18ecad12c6a86c6b5fc6e1a109b08))
- Add social media logo ([e900408](https://github.com/tifuchs/orGUI/commit/e900408ba5c756744ef3a56ac2f5fceaa3c4479c))
- OrGUI.py - add 'select all' functionality when removing curves in the integrated data plot ([cf26b34](https://github.com/tifuchs/orGUI/commit/cf26b348ed726aaa0d51e9960f02f05a520ca39b))
- [CTRplot] add CTR average functions ([915c2d5](https://github.com/tifuchs/orGUI/commit/915c2d53c8e825314d8fbaae856f320c3d10bd10))
- Added functionality to concatenate scans ([2aa59b1](https://github.com/tifuchs/orGUI/commit/2aa59b1233d4bbd8bb4aaaa21622d59ad3da9368))
- Add Bragg reflection rocking integration, add variable ROI sizes ([798269f](https://github.com/tifuchs/orGUI/commit/798269fa23edfd7008f08b17d3d8cccb8f20147d))
- Add missing ROIutils.py ([871f486](https://github.com/tifuchs/orGUI/commit/871f486ea5544ae6d7b487938bfadd38b156de7a))
- Add title showing hkl and s to rocking peak integr graph ([bf42ad6](https://github.com/tifuchs/orGUI/commit/bf42ad67ebe17e0a532074577d80bc518883adbc))
- [config] add compression and autoload setting to config file ([598b558](https://github.com/tifuchs/orGUI/commit/598b55853a135332634edf6b98464209758bf19c))
- Add setuptools_scm for automatic version numbering, update .gitignore ([0c33c71](https://github.com/tifuchs/orGUI/commit/0c33c7194e5a748f8ab8b86d72d4668022329b6e))
- Add cliff.toml for automatic CHANGELOG creation, replace CHANGELOG ([a155cdb](https://github.com/tifuchs/orGUI/commit/a155cdb23fbef6245fba8cf4fb42016f53a805b2))

### Changed

- [orGUI.py] - integrated curves can now be hidden/removed ([9c0a93e](https://github.com/tifuchs/orGUI/commit/9c0a93e39e95a773ef2a5b26a10326548d057e93))
- [ROIRoIntegr] drop requirement of rocking_integrate for s, H_0 and H_1 ([d2400f7](https://github.com/tifuchs/orGUI/commit/d2400f72cb19fd14c8a68acdfd0d6744e62d95b3))
- [Recip Space Nav dock] automatically expand reflection table to its ma size ([02e842a](https://github.com/tifuchs/orGUI/commit/02e842a967d7ed56d53348753af8a42599cfb341))
- [integration] Normalize croibg by the number of available pixels ([c126d5c](https://github.com/tifuchs/orGUI/commit/c126d5c2b2f5ae56189e61c784cfee9585365ab5))
- [RoIntegr] force integrated intensities positive, bugfix in error propagation ([cfbf202](https://github.com/tifuchs/orGUI/commit/cfbf2029d714131c538c4efe3874f3ed5b79eb74))
- Rework integration: Corrections are applied after background subtraction, add scaling of bg image as method of background subtraction, croi, bgroi are now raw counts, add Corrections and bg image as counter in database ([9067dbb](https://github.com/tifuchs/orGUI/commit/9067dbbce4722fe470ea884a8c57e02ba4a3f40f))
- Integrate functions return a status message to handle errors in batch processing ([152a713](https://github.com/tifuchs/orGUI/commit/152a7137f36d46672a616c6a62b85542f85f1d19))
- [Database] Hdf5Model set ownfile=False to avoid multiple open file handles, improve handling of harddrive faults and general I/O errors ([e682e56](https://github.com/tifuchs/orGUI/commit/e682e562fa660f4de4ed05c70ef5b6394e8c90ca))
- Backends.py / id31_tools.py - adjustments to correctly read and process beamtime ch7856 datasets ([d90b3ff](https://github.com/tifuchs/orGUI/commit/d90b3ffeb3984ea8d8dddf11930134d22e8e35f3))
- Merge branch 'master' of https://github.com/tifuchs/orGUI ([fddbe5e](https://github.com/tifuchs/orGUI/commit/fddbe5eeabf4b4df239146bdd17eac3bbd53cc92))
- [CTR] raise more descripitve error when CTR is not in CTRcollection ([a8586ed](https://github.com/tifuchs/orGUI/commit/a8586ed5ab237527e8d5fe458e2ad0c3cd4eee7a))
- Backends.py - further improved ch7856 scan opening backend function ([dfcd034](https://github.com/tifuchs/orGUI/commit/dfcd03439a9cb8a8cec94c71fcdbbb645688bfc4))
- OrGUI.py - optimize plotting of rocking scan curves: scans over panel edges (zero intensity in every point) are not plotted anymore and hiding + removing of curves is fixed ([a3a35fb](https://github.com/tifuchs/orGUI/commit/a3a35fbc3d7838f7029ee4ae20b0aa21ee3e4bac))
- ReciprocalNavigation.py - hotfix of allowed Bragg peak searching function ([9f3ce20](https://github.com/tifuchs/orGUI/commit/9f3ce205eb5b73fd625d84f1cd11da3649212de5))
- QReflectionSelector.py - typo ([5e22e16](https://github.com/tifuchs/orGUI/commit/5e22e166ab053c52fa8b99458cd1249378217790))
- Revert "ReciprocalNavigation.py - hotfix of allowed Bragg peak searching function" ([96aeb74](https://github.com/tifuchs/orGUI/commit/96aeb7432969ce5804008be8b7db849ca109a24f))
- ReciprocalNavigation.py - temporary fix of incorrect Qmax ([5195d34](https://github.com/tifuchs/orGUI/commit/5195d3441173449fc81692eff2d5134ed0fd068a))
- QReflectionSelector.py - catch errors if Bragg lists are empty ([77d24da](https://github.com/tifuchs/orGUI/commit/77d24dad49803b15dfb12aeee880b3d3866f2867))
- OrGUI.py - catch exception in rocking_extraction() ([109c18e](https://github.com/tifuchs/orGUI/commit/109c18e199d1456cc7a50d9df0ef97e4c423410d))
- Orgui.py - small improvements to 'load segemented scan' functionality ([37fc6b7](https://github.com/tifuchs/orGUI/commit/37fc6b7b4b46df099cf76c9086ad954c81bd730c))
- Improvements to scan concatenation function ([64a9351](https://github.com/tifuchs/orGUI/commit/64a93514214847e818c85b741e34bef9d47e10ca))
- Merge pull request #32 from tifuchs/ro-Bragg ([bc77481](https://github.com/tifuchs/orGUI/commit/bc774812da79075f2d5ae9905826f30a75638736))
- [UI] make background image action checkable - allows to disable the image subtraction ([3785b65](https://github.com/tifuchs/orGUI/commit/3785b650d85dbea0f75876d74d28437b0204ca4f))
- Merge branch 'master' into concatenate-partial-scans ([e97da08](https://github.com/tifuchs/orGUI/commit/e97da08bf7e4922b8444bcee1ceebb5c965a7854))
- Merge pull request #35 from tifuchs/concatenate-partial-scans ([c3b6a73](https://github.com/tifuchs/orGUI/commit/c3b6a73cc6ca27dd7f4cc7427a15b524bab72965))
- Update deprecated pyFAI dependencies and functions ([24a4600](https://github.com/tifuchs/orGUI/commit/24a4600332cf4c61754317d23105278b48246dec))

### Fixed

- [RoIntegr] Bugfix: apply mask to s and hkl to handle missing data ([a9829fb](https://github.com/tifuchs/orGUI/commit/a9829fb957d4da457e1bc2dd91f599202a92cb3c))
- [Database] Bugfix: fix gui lock when database has pending operations and main thread is performing database operations by allowing Qt to process events in wait loop ([1822c35](https://github.com/tifuchs/orGUI/commit/1822c353c918fc822a0039fe26f0abcb39a7ada0))
- [RoPeakIntegr] Bugfix: normalize F2 to integration range ([cda3540](https://github.com/tifuchs/orGUI/commit/cda3540077d21eba800551f8efb16a42e93bc4cc))
- [RoIntegr] Fix crash when scan is active which is not in database. ([30b6dd0](https://github.com/tifuchs/orGUI/commit/30b6dd0b92663bde6acd2ba5ffa2ab04f0beefc6))
- [UBcalc] Fix loading of refraction index from config file ([285cbc1](https://github.com/tifuchs/orGUI/commit/285cbc1410f74bbff199567a3f7d45e609c78f64))
- [RoPeakIntegr] Fix peak position interpolation ([cbc907e](https://github.com/tifuchs/orGUI/commit/cbc907e12d46bc5cd629f66aee237bf77eab0473))
- [RoPeakIntegr] Fix fit along rod by using intersect point s as abscissa instead of axis ([eb04e01](https://github.com/tifuchs/orGUI/commit/eb04e0104f1a24933c1397d814bf6976bc1ff1cd))
- [CTRuc] Fix loading of xpr files, fix layered_O water model ([e7d1d00](https://github.com/tifuchs/orGUI/commit/e7d1d005f10dd46d01eab0f0f23de3b8889c68e8))
- [CTRplot] fix typo ([17543ab](https://github.com/tifuchs/orGUI/commit/17543ab8b3b6ddb6ace46b394e9b61fe752d9d32))
- OrGUI.py - bugfix in sum image creation ([6b82ccb](https://github.com/tifuchs/orGUI/commit/6b82ccbd02de704c002982c2e42ef56db39e8e1c))
- [interlacedScan] fix crash if data is missing or corrupted ([ab6abcf](https://github.com/tifuchs/orGUI/commit/ab6abcfd18e825647db742ff452763315f4d5abe))
- OrGUI.py - fix uncaught exception in integrdata Plot function ([e0ca4fe](https://github.com/tifuchs/orGUI/commit/e0ca4fe320cc833032520ea1cb79498fb49f0eaf))
- [CTRuc] fix crash in calculation of cartesian position of atom in UnitCell ([774dfe3](https://github.com/tifuchs/orGUI/commit/774dfe34d7672e93830ebd810624bc876181f046))
- Fix loading of Bragg reflection rocking scans in roIntegr window ([0469f61](https://github.com/tifuchs/orGUI/commit/0469f61eb26c5503491791771e52a69282307a43))
- [CTRuc] fix loading of unit cell from file with single atom ([baf1d12](https://github.com/tifuchs/orGUI/commit/baf1d121d1c997ab9a9951f6c44acb36c6474aa8))
- [RoIntegr] fix saving of x y coordinates ([26221e1](https://github.com/tifuchs/orGUI/commit/26221e19313a727c7e1db39fb1c59b27d2f91660))
- [Ro1DIntegr] fix / remove integration range normalization for F2 calculation ([88403f3](https://github.com/tifuchs/orGUI/commit/88403f33641be594f46c10789705105901da2faf))
- Fix commit links in clffi.toml template ([eb2d82f](https://github.com/tifuchs/orGUI/commit/eb2d82f796a5f34798d7f78f7794ef73b9231332))

### Removed

- InterlacedScanLoader.py - delete function that was defined twice ([2b77efc](https://github.com/tifuchs/orGUI/commit/2b77efc13d33a788ae2846de29b459f10a8a0e63))
- [ROIoptions] remove size adjustment of static ROI, remove unused scale factor in ROIoptions dialog. ([4ad37ba](https://github.com/tifuchs/orGUI/commit/4ad37bae01d9d69e5b64025f5851fac5b43e4392))

## [1.3.0] (2025-07-22)

[dd31b00](https://github.com/tifuchs/orGUI/commit/dd31b000c631d1df3752d39ad2f71414322a49c0)...[ef10ee8](https://github.com/tifuchs/orGUI/commit/ef10ee8e805c13b60bb51361c38868da928d7116)

This is a feature release with major improvements and bugfixes (\> 40
commits), while still preserving mostly backward compatibility. One
potential incompatibility is the file format of saved rocking scans,
which had to be changed for more efficient saving of the data (see the
Changed section).

- Simplisitic rocking scan integration feature is now available.
  Conventional rocking scans (\$phi\$- or \$omega\$-scans) and
  reflectivity rocking scans are supported.

  - New tab in main GUI called \`Rocking scan integrate\`:
    - Select rocking scan in database, which was previously ROI
      integrated, right click and select show as rocking
      scan to display the data in the new tab.
    - 1D ROIs can be added in the Modify ROIs area. These
      are set relative to the calculated position of the reflection.
      They are shown in the graph as vertical lines.
    - ROI positons can be modified by drag+drop. ROIs which are selected
      as anchor are fixed when using the Fit between
      anchors feature.
    - Fit between anchors linearly interpolates between
      ROI anchor positions
    - Lorentz and footprint corections can be applied.
    - Saving/Loading of ROI locations

  This process is very manual at the moment, but will be improved in the
  future.

- Image integration is now thread-parallel and jit compiled while
  releasing the gil in numba nopython nogil mode, leading to a speedup
  of up to 10x (depending on available thread count).

- Refactored CTR calculation code, which is now split internally into
  the following modules:

  - CTRcalc: Main entry point for the CTR calculation, the end user
    should only use this top level API. API should be completely
    backward compatible.
    - CTRutil: reading of constants (form factors etc.), Fit
      functionality such as FitMixins, Parameter class, \...
    - CTRuc: fundamental structure factor calculation of simple unit
      cells and water models
    - CTRfilm: Code to construct coherent thin films on a substrate.
      Includes interfacial models such as Skellam interfacial roughness

Thanks to Finn for the rocking extraction script, changing the default
scan name, fixing the reflections view bug and the maintenance of the
ID31 backend.

### Added

- Add layer idx to UnitCell, split_to_layers function, EpitaxyInterface now uses splitted UnitCells, add EpitaxyInteface.toStr and fromStr ([dff4a41](https://github.com/tifuchs/orGUI/commit/dff4a411ae81b8033adc02214a671fade2552978))
- Add fit code to EpitaxyInterface using common LinearFitFunctions mixin class ([08a3637](https://github.com/tifuchs/orGUI/commit/08a3637f19020f9c0ceeff22d25f8ed571731ac6))
- Add Film class, stacking of layers in SXRDcrystal, F and eDensity of test examples correct ([50dc46c](https://github.com/tifuchs/orGUI/commit/50dc46c6a736502051a1c230d6956ab1b9814446))
- Add set_arci_result to numpy CTRopt, fit bugfixes ([905ab66](https://github.com/tifuchs/orGUI/commit/905ab66bcabea4f82200a0500ac649403360e2a3))
- [Qutils] add DataRangeSlider ([87f887a](https://github.com/tifuchs/orGUI/commit/87f887af725c47d90e7d297e715506313c8a4bc8))
- Add rocking scan tab, scolling through ro scans possible, CHANGE RO SCAN DATABASE STRUCTURE ([ca08c5e](https://github.com/tifuchs/orGUI/commit/ca08c5e21be4947f2eb6a1746244d4e29a5ab7b6))
- Add ROI modification layout, adding, deleting, modify, area correnctions dialog ([32da858](https://github.com/tifuchs/orGUI/commit/32da8586fa731f482704d319a649e7832e374908))
- Add Ro scan integration, fit anchor along rod ([a472e85](https://github.com/tifuchs/orGUI/commit/a472e858fecbbd4e00b6fde0ce94ef10bf24934e))
- Add button to create new empty database ([13000fd](https://github.com/tifuchs/orGUI/commit/13000fd5c06fc5ada9cb8157ca18c0ae44ed2eff))
- Add saving/loading of roi locations, many bugfixes ([475506d](https://github.com/tifuchs/orGUI/commit/475506d2bf9ed86ac1a9615993ab570dadbcbb43))
- Backends.py - added recent ESRF beamtime CH7149 ([599fdbe](https://github.com/tifuchs/orGUI/commit/599fdbe311de02cf8fb765ac377e1f2207889f07))
- Id31_tools.py - added fast potential counter 'scaled_potv2f' as auxilliary counter ([341fed2](https://github.com/tifuchs/orGUI/commit/341fed2cf6b799f86df03be5647bc72a7aea6f16))
- Add script to perform rocking scan extractions ([cfc307d](https://github.com/tifuchs/orGUI/commit/cfc307d719de4eabced3d3a0f54003b1aa3b5505))
- Add QM2 backend for loading files in tiff folder, change phi name ([47a862f](https://github.com/tifuchs/orGUI/commit/47a862fefa2948381680b5b50444e76d7d0b9ae4))
- Add rocking scan CTR integration ([3033893](https://github.com/tifuchs/orGUI/commit/3033893a61aa0244742a6a20cd7143208445a0f7))

### Changed

- Merge branch 'master' of https://cau-git.rz.uni-kiel.de/ieap/solid/ec/members/t_fuchs/orgui into thinfilm ([2d9ec33](https://github.com/tifuchs/orGUI/commit/2d9ec33bcd25ae830f39e22e82a77024d7b1f859))
- Id31_tools.py - include new fast potential data in aux counters ([e7203b8](https://github.com/tifuchs/orGUI/commit/e7203b87ff877031b852370eeb9e9749832cf792))
- Merge branch 'master' of https://cau-git.rz.uni-kiel.de/ieap/solid/ec/members/t_fuchs/orgui ([e41f843](https://github.com/tifuchs/orGUI/commit/e41f843d54824bc35fea4b0e99747f2f98d4659b))
- Merge branch 'master' of https://cau-git.rz.uni-kiel.de/ieap/solid/ec/members/t_fuchs/orgui into thinfilm ([0969ee9](https://github.com/tifuchs/orGUI/commit/0969ee9bec419ba74d41fd99a82a250d6836e418))
- [CTRuc] pos_cart, plot3d, pos_cart_all now are calculated takng into account the domain transformations ([7758dc8](https://github.com/tifuchs/orGUI/commit/7758dc8757bc2cda5040bec03befa8319c8924f2))
- Merge branch 'master' of https://cau-git.rz.uni-kiel.de/ieap/solid/ec/members/t_fuchs/orgui into thinfilm ([9a0c604](https://github.com/tifuchs/orGUI/commit/9a0c6048f2eecdd8a04505dd8c92d29bab3f7ed2))
- [rocking_int] save each rocking scan in its own hdf5 path in database ([224398a](https://github.com/tifuchs/orGUI/commit/224398a5b94451b58d5b883368433fc6de4407ea))
- Save rocking curves and integration ranges in 2d arrays for speedup, add auto zoom ([5403779](https://github.com/tifuchs/orGUI/commit/540377951d4980c0fb18b776039cb9246c72702c))
- Run roi rocking integration thread parallel using nogil and nopython mode of numba jit: speedup of factor 10 with many rois ([dbdd016](https://github.com/tifuchs/orGUI/commit/dbdd016606f167f15d8f2907dd8e2ddbe56350af))
- Ro integr switching anchors update cureve slider value ([8fba6f1](https://github.com/tifuchs/orGUI/commit/8fba6f11559842576e4b33074efb712d94c85bb2))
- Merge branch 'ro-peak-integr' of https://github.com/tifuchs/orGUI into ro-peak-integr ([b92c801](https://github.com/tifuchs/orGUI/commit/b92c801ce5690676976712875ca61411478874c9))
- Sort images according to suffix number ([4a2785a](https://github.com/tifuchs/orGUI/commit/4a2785a4fe9135ec808bd7b2f3a742557df7f21b))
- QScanSelector.py - more information in ddict in _onLoadScan ([309f567](https://github.com/tifuchs/orGUI/commit/309f5672b553f8c4088b24f41afde66c3f4a8faf))
- OrGUI.py - further reduce number of plots created by repeated rocking integration ([872cbb0](https://github.com/tifuchs/orGUI/commit/872cbb0f7d45bbcc204d84315e6f1174951fcfba))
- Merge branch 'master' of https://github.com/tifuchs/orGUI into ro-peak-integr ([88dfdcf](https://github.com/tifuchs/orGUI/commit/88dfdcfaaaf4b64dc41ca7dd637412cc5bf816e6))
- [QUBcalculator] calcRefection adjusts angle range to data angle range in fscan, addesses #6 ([6771ba4](https://github.com/tifuchs/orGUI/commit/6771ba4802495e344010d725f1f25c1916048b1f))
- Update CHANGELOG.rst for release 1.3.0 ([ecdac8a](https://github.com/tifuchs/orGUI/commit/ecdac8a4f0188f26fb33e09b16279d88d550ee23))
- Bump version to 1.3.0, update file headers, implement license formatting according to PEP 639 ([ef10ee8](https://github.com/tifuchs/orGUI/commit/ef10ee8e805c13b60bb51361c38868da928d7116))

### Fixed

- [DetCal] Bugfix of detector binning if number of pixels is divisible ([c2737f8](https://github.com/tifuchs/orGUI/commit/c2737f827d79e463cb94a851582c5c6e1199c9ba))
- [CTRuc] fix loading of xpr file with errors, save layerpos and layer_behaviour in xtal, xpr file ([cc49330](https://github.com/tifuchs/orGUI/commit/cc49330c151bb10dd2b0d0c0fb3ec45cab10352c))
- [CTRuc] Bugfix: layer_behaviour can now be read with '_', add better error message, fix // comment in layer_behaviour line ([109eae5](https://github.com/tifuchs/orGUI/commit/109eae513783a14981eae188e720b8ebd2eb7a44))
- QReflectionSelector.py - fix bug that deletes list of set reflections when toggling 'view - reference reflections' option ([9260ed5](https://github.com/tifuchs/orGUI/commit/9260ed549c1d45dc65b25c61ed8a1031285017e7))
- [DetCal] Bugfix: calculation of detector Q range fixed ([702bb31](https://github.com/tifuchs/orGUI/commit/702bb3191ba0c1047f26461be5b8344c3c3dc93f))
- Rocking integration: populate C_arr ([0e0959c](https://github.com/tifuchs/orGUI/commit/0e0959c4d47c21949efb64c84d63ec49dad7f1f5))
- Bugfix database.py - disable renaming of 'rois' group ([bd677c9](https://github.com/tifuchs/orGUI/commit/bd677c99b73583e3836b8efa6df5d721260a2ac0))

## [1.2.0] (2024-12-03)

[524579a](https://github.com/tifuchs/orGUI/commit/524579a1718dd07f6129b7ea1d301731d09fb768)...[dd31b00](https://github.com/tifuchs/orGUI/commit/dd31b000c631d1df3752d39ad2f71414322a49c0)


Release that adds 3 major improvements:

1.  Rocking scan integration is now properly supported with fast
    calculation speed and visulalization of the ROIs
2.  Backends can now be loaded as standalone python file. No
    reinstalling of orGUI is required to add new backends.
3.  hdf5 file handling has been improved to enable reloading of changed
    hdf5 files

In addition, there are many more smaller improvements and changes. See
the detailed changelog.

### Added

- Add QM2 / CHESS standalone backend ([1a21fbb](https://github.com/tifuchs/orGUI/commit/1a21fbb19128a1924d8337b3d2a20ea262a56b29))
- Add manual Euler rotation adjustment ([570a83a](https://github.com/tifuchs/orGUI/commit/570a83ad34215d9a4e265388251a1810622890fd))
- Add user provided fscan.name used as database key. Add batch processing example ([0013132](https://github.com/tifuchs/orGUI/commit/001313244543b8988b9180813089dacf67ee506b))
- Added rocking scan integration tab into ROI integration area ([a0452ce](https://github.com/tifuchs/orGUI/commit/a0452ceabf162ddb4776281a2aa2b1614759061d))
- Add ROI class RectangleBgROI, that handles center and background, enable drag+drop editing of static ROI ([5a5efeb](https://github.com/tifuchs/orGUI/commit/5a5efeb3d7d53ba7849fbe6d945e8efc33fb7b2c))
- Add angles, desired hkl and more descriptive name to saved rocking scan data ([5879ac7](https://github.com/tifuchs/orGUI/commit/5879ac72ce47f28e83f29a62475325e953140fd5))
- Add CTRfilm.InterfaceEpitaxy with Skellam distr, Refactor CTRcalc ([e2be0c4](https://github.com/tifuchs/orGUI/commit/e2be0c43cc7ee613a516cf09e7fedd50ce6853f4))
- Added notebook to test hdf5 file locking behaviour of orGUI ([4b63d68](https://github.com/tifuchs/orGUI/commit/4b63d6882d3f381dd02791de44571ff01d28ab5d))
- Add HDF5_USE_FILE_LOCKING explicitly in argparser. Display setting of hdf5 locking upon startup. ([a5f8360](https://github.com/tifuchs/orGUI/commit/a5f836050dbb31c1e1cb023a5ddc1530d5dc8d2f))
- Add PyOpenGL to optional dependencies ([d3eeb18](https://github.com/tifuchs/orGUI/commit/d3eeb1885f052ff94bcffe1915d4a8659ff99b31))

### Changed

- Merge branch 'master' into 'backend-rework' ([41cb06c](https://github.com/tifuchs/orGUI/commit/41cb06cd55fd77e38036cccd23de1417091fcaee))
- Merge branch 'backend-rework' into 'master' ([9f1bfc8](https://github.com/tifuchs/orGUI/commit/9f1bfc83c14ffdbff9f51b5cd529059fddbfbb7e))
- UniversalScanLoader.py - adapt class to reworked 'Scan' backend class ([3c36f67](https://github.com/tifuchs/orGUI/commit/3c36f677ee08eb490112301f52480d715b9a7495))
- Merge branch 'master' into optimize_rocking_scan_integration ([f5137fd](https://github.com/tifuchs/orGUI/commit/f5137fd45dca75e17ac63577f39ed677899ed59d))
- Enable loading of backend from config file, add example config file, addresses #1 ([5c9f9c8](https://github.com/tifuchs/orGUI/commit/5c9f9c80ba0da4204c099aad8373d6450c8f95ef))
- Merge branch 'optimize_rocking_scan_integration' into 'master' ([7760056](https://github.com/tifuchs/orGUI/commit/77600563a309618fd57985c34765958e4d66b6ef))
- OrGUI.py - reduce number of plots created by rocking scan integration to max. 30 ([e920d19](https://github.com/tifuchs/orGUI/commit/e920d1923db1078136938ffc2d94a40533ccb401))
- QScanSelector - change text in ROI selection tab ([a629c91](https://github.com/tifuchs/orGUI/commit/a629c91463b70ae0be7714039330105fd30a3d76))
- QUbcalculator.calcReflection now accepts hkl ndarray ([3cc8053](https://github.com/tifuchs/orGUI/commit/3cc8053f84a33153c7da2b5cc2833aaf9c5c0342))
- Make anglesToHKL fully numpy parallel ([0977667](https://github.com/tifuchs/orGUI/commit/09776674ea491e49db76c5ce5c825528d3907971))
- Enable showing of ro scan ROIs, update ro scan options to general H_1*s + H_0 scans, add fast parallel rocking roi, angles and hkl calculation ([829e328](https://github.com/tifuchs/orGUI/commit/829e328edc195ec9d93fc9e10ab21c58a25ea4c5))
- Cosmetic changes to rocking scan integration: ([a72a312](https://github.com/tifuchs/orGUI/commit/a72a312bdafa1bef02ca9b3cb987781b0e1d3ea0))
- Rocking scan: Set auto h and v roi size to median of differences ([704fc27](https://github.com/tifuchs/orGUI/commit/704fc2741143bd459332d10fa30f01d4ac3b2adb))
- QScanSelector - change hdf5 tree refresh (adapted from silx view) ([61e2209](https://github.com/tifuchs/orGUI/commit/61e2209e6be3d683e75edb0a92428f0aa467f42b))
- Main - disable hdf5 file locking when starting orGUI ([c64489c](https://github.com/tifuchs/orGUI/commit/c64489c4eb29319dcaff98aa05b2eb2df0650375))
- Pyproject.toml - change directory name of version path to fix a bug blocking installation on case sensitive UNIX systems ([9172eed](https://github.com/tifuchs/orGUI/commit/9172eedaf7b4698026375d50a6b37ab8d2915b38))
- Updated tests for file locking/synchronization of hdf5 ([bb359b4](https://github.com/tifuchs/orGUI/commit/bb359b499befe129c1c60423cbed078ab11ca5d0))
- QScanSelector.py - rework of hdf5 refresh by creating a new HDF5TreeView object ([0451171](https://github.com/tifuchs/orGUI/commit/0451171a9dd2cdb9713e8a46fe4bc42494f320c5))
- Merge branch '2-make-rocking-scan-integration-more-user-friendly' into 'master' ([bbd640f](https://github.com/tifuchs/orGUI/commit/bbd640f321f871cfc05dec64b5054493369948e5))
- Merge branch 'master' of https://cau-git.rz.uni-kiel.de/ieap/solid/ec/members/t_fuchs/orgui into hdf5_handling/locking ([85eae51](https://github.com/tifuchs/orGUI/commit/85eae515e9470418c691c3813844d023183ae2a0))
- Merge branch 'hdf5_handling/locking' into 'master' ([b147a88](https://github.com/tifuchs/orGUI/commit/b147a883225488a9a3bdc7ee1cf48aa25df13720))
- Update CHANGELOG.rst for version 1.2.0 ([c5647e3](https://github.com/tifuchs/orGUI/commit/c5647e350721160854a8a64b24a9bd7925c02a20))
- Bump version to 1.2.0 ([dd31b00](https://github.com/tifuchs/orGUI/commit/dd31b000c631d1df3752d39ad2f71414322a49c0))

### Fixed

- Fix typo ([be2bbb7](https://github.com/tifuchs/orGUI/commit/be2bbb76d3be6c131593a30868b83ed80090589e))
- Fix legacy auxillary_counter in ro scan integration ([c415ddc](https://github.com/tifuchs/orGUI/commit/c415ddc31655d2d327b194775980d4e29bdcac81))
- QScanSelector - fixed a bug that would cause a crash when searching for the Ewald intersect position in rocking scan integration ([e4a0299](https://github.com/tifuchs/orGUI/commit/e4a02999b8beee9f8a1022e7853e6b58bde3dbec))
- Fix crash when closing spec or fio files ([cf60616](https://github.com/tifuchs/orGUI/commit/cf60616b5b6fd749c42c341bae4ce15d642fa7e9))

### Removed

- OrGUI.py - delete old rocking extraction functionality from 'data processing' menu ([f99ed7c](https://github.com/tifuchs/orGUI/commit/f99ed7c42b57086c6afa471318feef209bee80d5))
- Remove confusing fallback lattice vectors message ([94228b4](https://github.com/tifuchs/orGUI/commit/94228b4ee20965c649c349205feeee0053d01eb5))

## [1.1.2] (2024-11-04)

[4d2621c](https://github.com/tifuchs/orGUI/commit/4d2621c78a9cce3f98a9c42ae1d61b4de71ea155)...[524579a](https://github.com/tifuchs/orGUI/commit/524579a1718dd07f6129b7ea1d301731d09fb768)

Bugfix version that fixes a bug that prevents installing orGUI on
certain systems. Full install target advertised on pypi is now
available. More minor bug fixes, see Fixed list.

### Added

- Legacy id31 backend: Add 'potv' counter for fast potential measurement to backend ([ead0c80](https://github.com/tifuchs/orGUI/commit/ead0c8061f029aec9d8565bdb3099dd59ccfc49e))
- Add experimental rocking extaction function - avoids opening images multiple times but currently still worse runtime ([dc1ba12](https://github.com/tifuchs/orGUI/commit/dc1ba120cce1a9c1363eb32488fe13d880efca7a))
- [HKL] add universal coordinate transform function ([c7ba52c](https://github.com/tifuchs/orGUI/commit/c7ba52c50e0e34443c264f6fa0aa194f0191a367))
- Add full install target ([37df17d](https://github.com/tifuchs/orGUI/commit/37df17dc815700f0c4893d7df5f5921f5c9bb90f))

### Changed

- Changes to ensure compatibility with the new Kiel diffraction setup: ([d532920](https://github.com/tifuchs/orGUI/commit/d532920b94aad0b00706fe20a319079a5381d51e))
- Improved error handling of 'Generate scan from images' function ([13a7615](https://github.com/tifuchs/orGUI/commit/13a7615ecb9f80421a7ef34d934a9e437e01f05e))
- Improved error handling of 'Generate scan from images' function (continuation) ([c8e530c](https://github.com/tifuchs/orGUI/commit/c8e530c5eca27d6a193ee472c45d67a2068ee828))
- Moved all ROI operations from integrate_beta to rocking_extract_beta -> massive runtime improvement ([b54424b](https://github.com/tifuchs/orGUI/commit/b54424bcc08fc5d512989f87dfa3ec1ee48365f5))
- Apply new backend API to old backends ([5f13a77](https://github.com/tifuchs/orGUI/commit/5f13a77aeb090486b2c7208f1678f0b7f1e403e2))
- Merge branch 'master' of cau-git.rz.uni-kiel.de:ieap/solid/ec/members/t_fuchs/orgui into backend-rework ([2502039](https://github.com/tifuchs/orGUI/commit/25020397d207352f8e47f9320258eb8464cf4271))
- Bump version to 1.2.0a ([4c31a41](https://github.com/tifuchs/orGUI/commit/4c31a4182d2f412798c5b26efb9eaf98adc28ce5))
- Enable loading of backend file from GUI, add standalone P21.2 backend ([080642c](https://github.com/tifuchs/orGUI/commit/080642ca0d81e362dbff0fd1ff029df872cac6df))
- Qt - replace deprecated .desktop() API for screen size determination ([8ae7c18](https://github.com/tifuchs/orGUI/commit/8ae7c187de44b0f3704a06a57968bb1213776348))
- Pyproject.toml - change directory name of version path to fix a bug blocking installation on case sensitive UNIX systems ([f06fe9c](https://github.com/tifuchs/orGUI/commit/f06fe9ce36d25f51c1b837ec5db3ff17614df45e))
- Update CHANGELOG.rst for version 1.1.2 ([3fa618c](https://github.com/tifuchs/orGUI/commit/3fa618c2d8d0d178a629aa99f78d97e3287ccbb3))
- Bump version to 1.1.2 ([524579a](https://github.com/tifuchs/orGUI/commit/524579a1718dd07f6129b7ea1d301731d09fb768))

### Fixed

- Bugfix in database creation and removal of dataset2 ([b8f52ca](https://github.com/tifuchs/orGUI/commit/b8f52ca608fbdd5557eea2c8b51e952820276884))
- OrGUI.py - bugfix in rocking inegtation dataset creation and code cleanup ([22df1eb](https://github.com/tifuchs/orGUI/commit/22df1eb3ec7a86a3c282793ffdbad74da6951f67))
- [CTRcalc] Bugfix: 3dplot correct lookup of atom radius ([88ddb0e](https://github.com/tifuchs/orGUI/commit/88ddb0e40449ebe6f1127091eb419bbac1dfd73f))
- [CTRcalc] Bugfix: plot3d, pos_cart_all: translate arg replace np.matrix.A1 by np.ndarray.flatten since translate is ndarray ([888cdd9](https://github.com/tifuchs/orGUI/commit/888cdd92e51ff9c69c995dc10e72f1131cfe527b))

### Removed

- OrGUI.py - remove old rocking scan function and tidy up code ([d4e2015](https://github.com/tifuchs/orGUI/commit/d4e2015474fdd64a37272098cd6646b86f2c8923))

## [1.1.1] (2024-08-07)

[ef0ac24](https://github.com/tifuchs/orGUI/commit/ef0ac24af01de0a14ac6ee79d74c85a2a906afc5)...[4d2621c](https://github.com/tifuchs/orGUI/commit/4d2621c78a9cce3f98a9c42ae1d61b4de71ea155)

Version that fixes a critical bug in the config file reading when a poni
file is given.

### Changed

- Bump version to 1.1.1 ([d6c84a3](https://github.com/tifuchs/orGUI/commit/d6c84a3c3feb7fc4c6c09e47aa543477c436df42))
- Update CHANGELOG.rst for 1.1.1 ([4d2621c](https://github.com/tifuchs/orGUI/commit/4d2621c78a9cce3f98a9c42ae1d61b4de71ea155))

### Fixed

- Fix setting of energy when poni file was given in config file ([434fae4](https://github.com/tifuchs/orGUI/commit/434fae447351f626c280ab3f5c4b94dd4124b813))

## [1.1.0] (2024-08-05)

[3e19f24](https://github.com/tifuchs/orGUI/commit/3e19f2416c4e5f501db57685028d45016cccc8ad)...[ef0ac24](https://github.com/tifuchs/orGUI/commit/ef0ac24af01de0a14ac6ee79d74c85a2a906afc5)

This is a release, that reworks the configuration widgets. It is aimed
to correctly show all configuration options in the GUI, which were
previously only available through the config file.

orGUI now also has a proper logo.

- Enable optional loading of atomic coordinates and lattice from
  [ASE](https://wiki.fysik.dtu.dk/ase/) supported files. These include:
  - Cif files
  - VASP / CONTCAR files
  - many more \... see [files supported by
    ASE](https://wiki.fysik.dtu.dk/ase/ase/io/io.html)

### Added

- Crystal config widget: add lattice parameter scaling options ([51f21be](https://github.com/tifuchs/orGUI/commit/51f21be749ab7d16d7ae2ae883b5edd40489ecda))
- Add CHANGELOG ([461c111](https://github.com/tifuchs/orGUI/commit/461c111de426e38c183b83029772dc13eb7cc257))
- Add new logo and application icon, bump version for release 1.1.0 ([5b6a6d5](https://github.com/tifuchs/orGUI/commit/5b6a6d5c79fc42a4dff2aa35f7e4354b95ea1b96))

### Changed

- Bump version to 1.1.0a ([401368c](https://github.com/tifuchs/orGUI/commit/401368c496da5b0da87b1019ca293c139310208c))
- Merge branch 'master' into origin/config-rework ([2489f51](https://github.com/tifuchs/orGUI/commit/2489f5107ddbb75809eae6ab1d436be99b079e99))
- Reorder machine parameter widget, add load poni file button, rework machine parameter config API ([4bc8379](https://github.com/tifuchs/orGUI/commit/4bc83799e082c7ef8178e017b24c373dc0bec97c))
- Universal scan loader: close file handles after data was read ([0a9f5d8](https://github.com/tifuchs/orGUI/commit/0a9f5d8e8828bdc0523a2cd9673d04b45c634bd7))
- Suppress runtime warning in Snell law when angle is smaller than critical angle ([346547f](https://github.com/tifuchs/orGUI/commit/346547f848acafb1f3d8161019a3a5c9f2cc4a1b))
- Merge branch 'config-rework' into 'master' ([ce7ac6f](https://github.com/tifuchs/orGUI/commit/ce7ac6f29cc0d2d33266c9684e1de1dc635f194f))
- Machine parameter widget: make detector select button smaller, add border and title to detector geometry ([b0f7b19](https://github.com/tifuchs/orGUI/commit/b0f7b19f02b7f329b33d59c9562dabf639f2dffd))
- Enable loading of ASEsupported file types, Fix loading minimal config file, add button to show atoms in unit cell ([51d857c](https://github.com/tifuchs/orGUI/commit/51d857c5ef12c3381d136513b8c0b719781681d3))
- Rework example config files, now showing multiple options of config files ([210604f](https://github.com/tifuchs/orGUI/commit/210604f63b32fe2cdc367dbd96fd24df45a05a7e))
- Delete deprecated config in root directory, now config only in examples folder ([3ff9b9d](https://github.com/tifuchs/orGUI/commit/3ff9b9d487faa867cbeb9b5c5bcd2b7553a5aaa5))
- Update CHANGELOG.rst for 1.1.0 ([ef0ac24](https://github.com/tifuchs/orGUI/commit/ef0ac24af01de0a14ac6ee79d74c85a2a906afc5))

### Fixed

- Fix detector loading order, now starts gui ([956d652](https://github.com/tifuchs/orGUI/commit/956d6527bb5942f386671a1a206330f50b1bf677))
- Fix QLayout double assignment warning ([8950c9e](https://github.com/tifuchs/orGUI/commit/8950c9edd8841891f0f9ba15999345865d3419f3))
- Fix CHANGELOG syntax errors ([873501d](https://github.com/tifuchs/orGUI/commit/873501d561f77c68404b52dd08cc43f44d20192a))

### Removed

- Remove detector size from scan creator widgets, should now only be edited in machine config ([24bdffa](https://github.com/tifuchs/orGUI/commit/24bdffa41e2c61f3f57dc61cb2ccd11298daa4f7))
- Remove config warning on startup, add version number to splashscreen ([c36c05e](https://github.com/tifuchs/orGUI/commit/c36c05ea421239da28bf49d4adff7f0654f40f90))

## [1.0.1] (2024-07-31)

[4f09f6a](https://github.com/tifuchs/orGUI/commit/4f09f6aea2abc092e5b370ad444e0d7e5441f64a)...[3e19f24](https://github.com/tifuchs/orGUI/commit/3e19f2416c4e5f501db57685028d45016cccc8ad)

First release version on PyPi, project availabe under
[orGUI](https://pypi.org/project/orGUI/). This version aimes to fix some
minor issues from the first release and provides a minimal documentation
in form of a README and a diffractometer image showing the geometry.

### Added

- Add README and pyproject.toml ([3289878](https://github.com/tifuchs/orGUI/commit/3289878b71f0ea9e79fe2e4dfa9b154afb29f7dc))
- Add diffractometer image to README ([42f0081](https://github.com/tifuchs/orGUI/commit/42f00813635876263b966ae0ea83fec17f3ffe34))
- Add link to Busing Levy paper ([1dd4057](https://github.com/tifuchs/orGUI/commit/1dd4057b553e5172d822a64d9749fadcc34ad8be))
- Add DetectorDialog and pyFAI config widget ([08e20d8](https://github.com/tifuchs/orGUI/commit/08e20d85860d46aa09e7e8ab68bd9d2ebd4f7970))
- Add getting started section to README ([bb6fb2c](https://github.com/tifuchs/orGUI/commit/bb6fb2c861549f44f86e9185deaca78ad4c4ac35))

### Changed

- Update About statement, dependency fixes and include diffractometer image ([bb4fd67](https://github.com/tifuchs/orGUI/commit/bb4fd6724f7ec040cea1a764afa79cce5afb32cc))
- Adjust diffractometer image size ([d231004](https://github.com/tifuchs/orGUI/commit/d231004aa524e9d00f7cd4295213c79eef2f33c8))
- Correct phi-circle rotation axis ([5ab7611](https://github.com/tifuchs/orGUI/commit/5ab76111bce35f8cc6d94c2788ca0b26cf77df92))
- Correct Latex code ([6d31682](https://github.com/tifuchs/orGUI/commit/6d31682f9293b41b600a3574ba577ce42e6a71f5))
- HKLVlieg.py - np.product -> np.prod ([38348f1](https://github.com/tifuchs/orGUI/commit/38348f11dc0787413b5b1513b7bc9bcf93fda22e))

### Fixed

- Fix crash if hkl reflection is seached for static roi and no scan is selected ([6f6546d](https://github.com/tifuchs/orGUI/commit/6f6546ddeb9dd2a28ab11963f9d6d7bd694a78c8))
- Fix README typo and add clarification about poni files ([3e19f24](https://github.com/tifuchs/orGUI/commit/3e19f2416c4e5f501db57685028d45016cccc8ad))

### Removed

- Remove deprecated setup.py, Bugfix: enable loading of poni file relative to config file location, enable start of orGUI without config ([1cdae82](https://github.com/tifuchs/orGUI/commit/1cdae8219ab3124d28231bb2282e2475729ffec8))

## [1.0.0] (2024-06-29)

[fd93969](https://github.com/tifuchs/orGUI/commit/fd939694709cba6371affe51e047dc50f61ab748)...[4f09f6a](https://github.com/tifuchs/orGUI/commit/4f09f6aea2abc092e5b370ad444e0d7e5441f64a)

- First public release on Zenodo with the DOI
  <https://doi.org/10.5281/zenodo.12592486>.
- Added MIT License
- Merged with *datautils*
  - Backends
    - datautils.xrayutils.id31_tools -\>
      orgui.backend.beamline.id31_tools
    - datautils.xrayutils.P212_tools -\>
      orgui.backend.beamline.P212_tools
    - datautils.xrayutils.ID31DiffractLinTilt -\>
      orgui.backend.beamline.ID31DiffractLinTilt
    - datautils.xrayutils.fio_reader -\>
      orgui.backend.beamline.fio_reader
  - Physics / Diffraction caculations
    - datautils.util -\> orgui.datautils.util
    - datautils.xrayutils.unitcells -\>
      orgui.datautils.xrayutils.unitcells (and .bul files therein)
    - datautils.xrayutils.test -\> orgui.datautils.xrayutils.test (and
      datautils test code therein)
    - datautils.xrayutils.CTRcalc -\> orgui.datautils.xrayutils.CTRcalc
    - datautils.xrayutils.\_CTRcalc_accel -\>
      orgui.datautils.xrayutils.\_CTRcalc_accel
    - datautils.xrayutils.CTRopt -\> orgui.datautils.xrayutils.CTRopt
    - datautils.xrayutils.CTRplotutil -\>
      orgui.datautils.xrayutils.CTRplotutil
    - datautils.xrayutils.DetectorCalibration -\>
      orgui.datautils.xrayutils.DetectorCalibration
    - datautils.xrayutils.element_data -\>
      orgui.datautils.xrayutils.element_data
    - datautils.xrayutils.HKLVlieg -\>
      orgui.datautils.xrayutils.HKLVlieg
    - datautils.xrayutils.ReciprocalNavigation -\>
      orgui.datautils.xrayutils.ReciprocalNavigation
- *orGUI* is now a standalone package, that has only publicly available
  dependencies.

### Added

- Added a universal backend to load raw images as a scan into orGUI ([01e49bd](https://github.com/tifuchs/orGUI/commit/01e49bdb4734d9bb021154397c0896fd05a9d1a3))
- [scan loader backend] Add license ([4f09f6a](https://github.com/tifuchs/orGUI/commit/4f09f6aea2abc092e5b370ad444e0d7e5441f64a))

### Changed

- MERGE WITH DATAUTILS ([23253ce](https://github.com/tifuchs/orGUI/commit/23253ce989de52d245515795747e49c2988dffe0))
- Deleted old code - universalScanLoader.py ([6d97730](https://github.com/tifuchs/orGUI/commit/6d977301561ad09a351d5b5b7036360d76055d51))
- Include universal scan creation functionality into orGUI.py ([382d918](https://github.com/tifuchs/orGUI/commit/382d918f352d926c6bdb1a72bb8d186ddea3f1e6))
- Merge branch 'master' of cau-git.rz.uni-kiel.de:ieap/solid/ec/members/t_fuchs/orgui ([59b93df](https://github.com/tifuchs/orGUI/commit/59b93df5bf4d62cf26fce766c0a9ce651c36e141))

## [0.9-alpha] (2024-06-29)

- Last inoffical version before merge with the internal dependency
  *datautils*.
- Last version under *All rights reserved*.

### Added

- Add gitignore ([5a0365e](https://github.com/tifuchs/orGUI/commit/5a0365e5fa93de9a5d9a6dcc779778aa1e83366a))
- Add gitignore ([aeb956c](https://github.com/tifuchs/orGUI/commit/aeb956c9ede7444c95084f3010611c333bb4eb79))
- [packaging] added setup.py, orGUI can now be started as shell script ([7838377](https://github.com/tifuchs/orGUI/commit/7838377fc5280ae52d8306e367a4194d94d70c73))
- Add simulation mode: dummy scan, add licence statements ([bfec4d9](https://github.com/tifuchs/orGUI/commit/bfec4d9064085678b3197eb8e1263cda62983823))
- Add ROI integration along trajectory in reciprocal space ([86c3971](https://github.com/tifuchs/orGUI/commit/86c3971304f0db4684df0227265e335ca3e26065))
- Add crystal name in config file ([6f41c09](https://github.com/tifuchs/orGUI/commit/6f41c0960092727f8ba14eee0f94cbe0a558b6b0))
- Add basic angle-of-incidence scan capability: mu-scan ([b53e539](https://github.com/tifuchs/orGUI/commit/b53e539ce2623f840fbf36c0f00bc2ff2cc907a5))
- Add integration database ([4d9b5fd](https://github.com/tifuchs/orGUI/commit/4d9b5fd9b31867184ce35929d4a96aa1d46e5d7c))
- Add static ROI timescan integration ([843a0c6](https://github.com/tifuchs/orGUI/commit/843a0c6f1b16d0c3d05021e05ae50841448f3be9))
- Add beamtime id select box for P21 backend ([6a87369](https://github.com/tifuchs/orGUI/commit/6a8736914dc03d71bd1e7c3cc980b72c207ab32d))
- Add ArrayTableDialog for editing of arrays; Add button and viewer to exclude images from max/sum image ([d973d23](https://github.com/tifuchs/orGUI/commit/d973d23e2976c0e1bfc08dc6d88420628905d245))
- Add splash widget during orGUI loading ([2ff08f4](https://github.com/tifuchs/orGUI/commit/2ff08f4c6965b9608e9f02233c157c1bfdc37501))
- Add qutils.py for error messages in QMessageBox, where the python error traceback is in the detail text ([6ea625c](https://github.com/tifuchs/orGUI/commit/6ea625ca6801292026f192847b6f06327649b260))
- Add ipython console widget for debugging and orGUI scripting ([787a56c](https://github.com/tifuchs/orGUI/commit/787a56c8e50f1b7c4d6476e4f73e5d0e7d3ee011))
- [examples] add example poni file, used by config_id31_EBS ([ad171e9](https://github.com/tifuchs/orGUI/commit/ad171e94ee37de45f535ef29beae2bd1509c4353))
- Add U Matrix edit dialog, currently without functionality, some bugfixes ([d430397](https://github.com/tifuchs/orGUI/commit/d430397b0cfa9e79f0176371836a7f5d6a4275fc))
- Add enable OpenGL rendering action to plot toolbars ([7386161](https://github.com/tifuchs/orGUI/commit/738616112bdd4cdff2f4ad765f1a884af2501569))
- Add machine params view action ([9ed9b4a](https://github.com/tifuchs/orGUI/commit/9ed9b4ab144cd20d6e64093ec7ff3edee1500b7f))
- Add reflection select dialog, add mu scan as simulation scan ([a8cbff7](https://github.com/tifuchs/orGUI/commit/a8cbff71d2ff791fd4ba9d20ee963b0fcee6634e))
- Add default backends and toggle button for selecting the backend/beamtime id ([ba158cc](https://github.com/tifuchs/orGUI/commit/ba158cc7565c04d1fe11eb7139c0bb212fd7ac3a))
- Add pytz as dependency ([337bbc7](https://github.com/tifuchs/orGUI/commit/337bbc73968ddc64aabf17676a46acc2bce0c2fe))
- [UB] add calculate miscut function ([fd93969](https://github.com/tifuchs/orGUI/commit/fd939694709cba6371affe51e047dc50f61ab748))

### Changed

- Initial commit ([6157e95](https://github.com/tifuchs/orGUI/commit/6157e9581bf2d9b5dea457007a5a41d9257fcec7))
- Initial bliss scan / nexus file support ([b9dfe4a](https://github.com/tifuchs/orGUI/commit/b9dfe4a756b13a8bdb25a7ac17f6162bca533812))
- Initial configfile support, reads from ./config ([abd14b0](https://github.com/tifuchs/orGUI/commit/abd14b08659ba5536b502a4587f957596518aa0d))
- Put in option for poni file caibration, pyFAI not yet in use!, example config ([e119397](https://github.com/tifuchs/orGUI/commit/e119397a0c43c2e075713db90d96143fa7412b1b))
- Make readConfig for pyFAI ([62c2b31](https://github.com/tifuchs/orGUI/commit/62c2b313169a31548617c474c4488b7697fbf015))
- Config looks good, also HKL calculation seem to work with pyFAi ([cc7bc8a](https://github.com/tifuchs/orGUI/commit/cc7bc8a30c15278f1c77023a9fd9c6d52a8c5c56))
- First version pyFAI, has to be tested ([4aca32b](https://github.com/tifuchs/orGUI/commit/4aca32b5caad6f8c3055abab6aa43ab7df442731))
- Load config, but some bug: shows save dialog instead of load ([d26a605](https://github.com/tifuchs/orGUI/commit/d26a6052a273dc96bc668c25bd11d6c63b7b246a))
- Fit of lattice parameters ([3d7a348](https://github.com/tifuchs/orGUI/commit/3d7a34896e61119749db60b71176ddc3b13a1955))
- [11006959] snapshot after beamtime ([183a75c](https://github.com/tifuchs/orGUI/commit/183a75cf278bdf4e873cb7f9265b396ed80ee3d4))
- Merge branch 'CrudeP212' into 'master' ([483e19a](https://github.com/tifuchs/orGUI/commit/483e19a0a3573ec4c8926c30b09eb1c77fc2147d))
- Major rework: drop spec support, now only working with silx and bliss ([53f1fdb](https://github.com/tifuchs/orGUI/commit/53f1fdbb079799a03f0006f95c0400bcc2bd816b))
- Now also acccepts only one reflection for ub calculation in z-mode ([ef5550e](https://github.com/tifuchs/orGUI/commit/ef5550e96d48d6718c6bbba38abda02761971835))
- [Navigation] calcuation of available CTRs ([7122e77](https://github.com/tifuchs/orGUI/commit/7122e77680907d1fd49f37f1c3d9ede0622460a8))
- [navigation] saving of available CTRs; some bugfixes and error reporting ([7f97a83](https://github.com/tifuchs/orGUI/commit/7f97a83fdcfa1a0d58dd9520cee167e0fc0b6d9c))
- Improved the avaliable CTR save option and added a save-dialogue ([b619ada](https://github.com/tifuchs/orGUI/commit/b619ada54e0e5dec5c534c34762493d3d7c0e38a))
- Reenable legacy ch5523 bliss support ([6d2e011](https://github.com/tifuchs/orGUI/commit/6d2e01170911785a1a0d3e246f036086da713287))
- Toggle view of reference reflections ([edc622c](https://github.com/tifuchs/orGUI/commit/edc622c56cebe363b4eaccd76e0cab9ed4f14c43))
- Moved crystal parameters to a dialog ([14f8165](https://github.com/tifuchs/orGUI/commit/14f816541d60f348fb89069f6e60fa490c15de7d))
- Crystal bulk unit cell import, new crystal dialog ([2d37ed4](https://github.com/tifuchs/orGUI/commit/2d37ed402d8c029c1cae6ffd6c6aa3204c35b4ac))
- Show Bragg reflections ([c461215](https://github.com/tifuchs/orGUI/commit/c46121562c79b757fc0d5d8ffe245ad7c06b4bc7))
- ROI integration working. Reorder GUI ([a198274](https://github.com/tifuchs/orGUI/commit/a198274d371d68e982a7b6d032b26b5b28ba8960))
- Dialog for setting cpu count, prev and next icons in image slider ([4e8d09f](https://github.com/tifuchs/orGUI/commit/4e8d09fe65e1fb91f558eaf52c73bf52d6ef08cb))
- Integrated curves plot, add integration correction factors ([1419573](https://github.com/tifuchs/orGUI/commit/1419573a6fa1e9cc6bc45c8d2bc85d8fe2af51d0))
- Enable export of angles of allowed Bragg reflections ([880c427](https://github.com/tifuchs/orGUI/commit/880c427c1ca72747f70ee9887708220209ecce46))
- Update popup messages for Bragg reflection calculation. Fix crash when selecting showROI without scan ([702d9e4](https://github.com/tifuchs/orGUI/commit/702d9e421fe9cf2514af67dfa0532411ec316433))
- Scan loading now handled in dedicated backend file ([03643d1](https://github.com/tifuchs/orGUI/commit/03643d1724c21c9fffba9fa1d784bc703eb6c0d6))
- Mask persists when changing images ([495c19b](https://github.com/tifuchs/orGUI/commit/495c19b601141f78314d7544ffbfa46e026ace28))
- Reorder Max/sum/alpha buttons into toolbar ([4601033](https://github.com/tifuchs/orGUI/commit/46010332da6a84817ee3fc2b4a87a4fc5641c31a))
- Static scan support, i.e. without diffractometer motor movement ([8f30d7e](https://github.com/tifuchs/orGUI/commit/8f30d7e440c52caaf8c052b4e2fbec195677878a))
- Make splash screen prettier by setting transformation to smooth ([03ca642](https://github.com/tifuchs/orGUI/commit/03ca6422c0c42f3ea669a63ac0e16bff0eb72e9f))
- Complete functionality of edit U dialog. Add user defined defaults of U in backend ([e0a6203](https://github.com/tifuchs/orGUI/commit/e0a620361621a5ea7b69da851b4957bb9c22adc3))
- Update for datautils.HKLVlieg API change ([13a52f0](https://github.com/tifuchs/orGUI/commit/13a52f0ed17604245cdbeba8dc8ff2c6c4aa0543))
- Rework reflection selector view ([d7af74d](https://github.com/tifuchs/orGUI/commit/d7af74d37e2009a5411cb00f4462ed64ba7e5219))
- Replace hkl reflection search by toolbar, fix UncaughtHook, center window on screen ([5305bae](https://github.com/tifuchs/orGUI/commit/5305bae5aad264326182c250f51ba150593aed1c))
- Enable static hkl reflection search ([ebcbac6](https://github.com/tifuchs/orGUI/commit/ebcbac6cf2e7f6ffd2a86649738d3c064cfc8b0f))
- Correct a bug in refraction correction before calling intersectEwald function ([b3dc3c6](https://github.com/tifuchs/orGUI/commit/b3dc3c6924b5cb4540931ef0b44c5965cfdf2558))
- Saving and loading of UB matrix as ascii, manual U edit updates plot ([25f0761](https://github.com/tifuchs/orGUI/commit/25f07611452dd76647f0e0047038e68e059e8b34))
- Separate search for hkl reflection function ([e1b5755](https://github.com/tifuchs/orGUI/commit/e1b5755fce35a6d9c3d41bf31e30c027204bab50))
- Enable U calculation in time scans, use correct omega, if provided ([5b75220](https://github.com/tifuchs/orGUI/commit/5b75220ffc281503c389d87c7957037e0d795cc2))
- Update dependencies ([d56a9ab](https://github.com/tifuchs/orGUI/commit/d56a9abb17e5d94509302b66cda2abfb591fb5c3))
- Replace orGUI.py ([0cf55c9](https://github.com/tifuchs/orGUI/commit/0cf55c9f2058d96e207a392e560f38d93f94347a))
- Merge branch 'rocking_scan_integration' into 'master' ([0fc2566](https://github.com/tifuchs/orGUI/commit/0fc2566dbe504b2456de1e34160d226d6bbf8be6))
- Rocking scan integration: Use GUI settings for mask, solidangle etc., add selection dialog for reflection, add more error handling, Bugfix: race condition with async database writes ([7076b8f](https://github.com/tifuchs/orGUI/commit/7076b8f49e58c52c3c0ee64a1968023f83278b4b))

### Fixed

- Fix hkl api change, now also save xmirror ([30f76d1](https://github.com/tifuchs/orGUI/commit/30f76d13029828a0a0cd3afea58447eb984cd03c))
- Remove duplicate notifications in machine and crystal params widget during startup ([19991a1](https://github.com/tifuchs/orGUI/commit/19991a1d9e657232e3dcb8e0f16aea9e29db05a9))
- Calc U with > 2 reflections working again. Multithreaded reading of images ([97bf694](https://github.com/tifuchs/orGUI/commit/97bf6949e3e5cd70f26849d48f2e9614c081b08a))
- Fix program crash if mask is not selected ([d54f403](https://github.com/tifuchs/orGUI/commit/d54f4034754830f30944cde6cc3837bf43751302))
- Reflection select on dclick, Bugfix: ROI masking, Add title and aux counters to database, Add ROI offset ([fd150f1](https://github.com/tifuchs/orGUI/commit/fd150f11860a0063873d777ef08b70d4694563d2))
- Fix missing square in number of pix ratio in uncertainty progression ([bf90b2b](https://github.com/tifuchs/orGUI/commit/bf90b2bc027646df83cd8b7f952ddcf5509e61ed))
- [scan loading] Bugfix: enable image loading from menu if auto load is deselected ([4349d01](https://github.com/tifuchs/orGUI/commit/4349d01ab6a42e1930d47e89f89e0eda5fb45a64))
- Display image enable in taskbar button, fix load array if only one element in ArrayEditWidget ([c0506fb](https://github.com/tifuchs/orGUI/commit/c0506fb734226d32d2b58040d722334e8c5b707c))
- Completely clear Bragg reflection cache when setting new reflections ([c77773b](https://github.com/tifuchs/orGUI/commit/c77773ba79c9bea3086aa78fa4efbef55ca88866))
- Enable chi,phi in calculation of possible Bragg reflections / CTRs ([ab18d04](https://github.com/tifuchs/orGUI/commit/ab18d04bd658d0cae9e675b6523cb3fea0590886))
- [UBCalc] fix API change in datautils.DetectorCalibration.rangedelgam ([7a5c344](https://github.com/tifuchs/orGUI/commit/7a5c344e20dd45397bee8d2c75f962eb357f39d0))
- Exclude CTR reflections from the mirrored detector plane from calculations by checking the maximum and minimum available Q of the detector ([3bb97f2](https://github.com/tifuchs/orGUI/commit/3bb97f205528046cfc6d6512cc5267a817027ce6))
- Fix critical bug in SIXC angle view if reflections are deleted ([2333a6f](https://github.com/tifuchs/orGUI/commit/2333a6fbf4c0054ff2e5c633c07fba7401e8b1d3))
- Convert float to int as splashscreen scaling argument ([bdae0a3](https://github.com/tifuchs/orGUI/commit/bdae0a3e52ed31858c6c4cebb101d45efc6aa790))
- Bugfix numpy>1.20: rename np.float -> np.float64 ([6b016d6](https://github.com/tifuchs/orGUI/commit/6b016d6469e192171b0d4153a1c479426792146e))

### Removed

- Remove old reflection spinbox edit, add toolbar buttons instead ([7854679](https://github.com/tifuchs/orGUI/commit/7854679ff96344e40fd9d88334672309f62be219))
- Remove array index widget from ArrayTableWidget, as only 1D or 2D arrays are used ([25796e2](https://github.com/tifuchs/orGUI/commit/25796e21f103a2a6f9e4dd8544e941bd7886605d))

[1.4.0]: https://github.com/tifuchs/orGUI/compare/1.3.0..v1.4.0
[1.3.0]: https://github.com/tifuchs/orGUI/compare/1.2.0..1.3.0
[1.2.0]: https://github.com/tifuchs/orGUI/compare/1.1.2..1.2.0
[1.1.2]: https://github.com/tifuchs/orGUI/compare/1.1.1..1.1.2
[1.1.1]: https://github.com/tifuchs/orGUI/compare/1.1.0..1.1.1
[1.1.0]: https://github.com/tifuchs/orGUI/compare/1.0.1..1.1.0
[1.0.1]: https://github.com/tifuchs/orGUI/compare/1.0.0..1.0.1
[1.0.0]: https://github.com/tifuchs/orGUI/compare/0.9-alpha..1.0.0
[0.9-alpha]: https://github.com/tifuchs/orGUI/tree/0.9-alpha

<!-- generated by git-cliff -->
