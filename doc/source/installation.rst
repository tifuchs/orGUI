Installation
============

orGUI can be installed from PyPI with its full optional dependency set:

.. code-block:: bash

   pip install orGUI[full]

For a smaller installation with only the core dependencies, use:

.. code-block:: bash

   pip install orGUI

The application is started from the command line with:

.. code-block:: bash

   orGUI

To start orGUI with a configuration file, pass the file path as the first
argument:

.. code-block:: bash

   orGUI <configfile>

Machine and crystal parameters are available from the GUI through
``Config -> Machine parameters`` and ``Config -> Crystal parameters``.
Example configuration files are distributed with the source code release.

Data Loading
------------

orGUI is designed to load scans stored as Nexus files or as SPEC files with
separately stored detector images. The file browser on the left side of the
main window displays the scans found in the selected file, and a scan is opened
by double-clicking it.

Beamlines often use different counter names, metadata conventions, and image
locations. orGUI therefore uses backend modules for beamline-specific loading.
Current backends support ID31 at ESRF and P21.2 at DESY. Support for another
beamline usually requires adding a backend in ``orgui/backend`` that normalizes
the scan metadata into orGUI conventions. See :doc:`beamline_backends` for
backend usage from the GUI or config file and for the backend implementation
structure.

For simple image stacks without scan metadata, use
``File -> Generate scan from images``.
