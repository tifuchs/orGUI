Acceleration Backends
=====================

orGUI has two independent acceleration selectors: one for ROI image
integration and one for CTR structure-factor calculations. Both selectors use
the same process-global syntax. Set ``ORGUI_ACCEL_BACKEND`` before importing
the affected module, before starting ``orGUI``, or before launching a batch
script.

ROI Integration
---------------

ROI integration uses the C++ backend by default. If the C++ extension is not
available and no backend is requested, ROI integration falls back to NumPy.
This is the normal path for the GUI and for batch integrations that call
``orgui.app._roi_sum_accel``.

To select the ROI backend globally for the current process, set
``ORGUI_ACCEL_BACKEND``:

.. code-block:: bash

   # Default: C++ ROI summing and C++ fitted-background support.
   export ORGUI_ACCEL_BACKEND=cpp

   # Optional: Numba ROI summing for strip-background integration.
   export ORGUI_ACCEL_BACKEND=numba

   # Disable acceleration and use the explicit NumPy fallback paths.
   export ORGUI_ACCEL_BACKEND=numpy

Valid values are ``cpp``, ``numba``, and ``numpy``. If the variable is unset,
``cpp`` is used when available. ``numpy`` means no ROI acceleration backend is
selected; the application uses its explicit NumPy fallback integration code.
The Numba backend is imported only when ``ORGUI_ACCEL_BACKEND=numba`` is set
or when it is selected explicitly in Python:

.. code-block:: python

   import orgui.app._roi_sum_accel as roi_sum

   roi_sum.set_accel_backend("cpp")
   roi_sum.set_accel_backend("numba")
   roi_sum.set_accel_backend("numpy")  # disable acceleration

The Numba ROI backend mirrors the C++ strip-background and image accumulation
functions:

* ``processImage_Carr``
* ``processImage_bg_Carr``
* ``calcMaxSum``
* ``calcBgSub``
* ``calcMaxSum_bg``

Polynomial fitted-background integration is C++ only. When the C++ extension
is available, ``processImage_polybg_Carr`` and ``interpolate_polybg_croi``
remain available from ``orgui.app._roi_sum_accel`` even if
``ORGUI_ACCEL_BACKEND=numba`` selects Numba for the strip-background counters.

CTR Structure Factors
---------------------

CTR calculations default to the C++ backend when it is available. If the C++
extension is unavailable and no backend is requested, they run the explicit
NumPy structure-factor code path. Numba is never selected or imported
automatically.

To select the CTR backend globally for the current process, set
``ORGUI_ACCEL_BACKEND``:

.. code-block:: bash

   export ORGUI_ACCEL_BACKEND=cpp
   export ORGUI_ACCEL_BACKEND=numba
   export ORGUI_ACCEL_BACKEND=numpy

Alternatively, select the backend in Python:

.. code-block:: python

   from orgui.datautils.xrayutils import CTRuc

   CTRuc.set_accel_backend("cpp")    # default when available
   CTRuc.set_accel_backend("numba")  # optional Numba backend
   CTRuc.set_accel_backend("numpy")  # disable acceleration

The active backend is stored in ``CTRuc.CTR_ACCEL_BACKEND``. To check whether
the optional Numba backend can be loaded, use:

.. code-block:: python

   CTRuc.ctr_numba_accel_available()

Calling ``CTRuc.ctr_numba_accel_available()`` imports the optional Numba CTR
module when it has not already been probed. ``CTRuc.set_accel_backend("numba")``
raises ``ValueError`` if the optional Numba backend cannot be imported.
