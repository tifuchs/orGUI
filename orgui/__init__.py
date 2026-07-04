# /*##########################################################################
#
# Copyright (c) 2020-2026 Timo Fuchs
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.
#
# ###########################################################################*/
"""Module descripiton"""

__author__ = "Timo Fuchs"
__copyright__ = "Copyright 2020-2026 Timo Fuchs"
__credits__ = []
__license__ = "MIT License"
__maintainer__ = "Timo Fuchs"
__email__ = "tfuchs@cornell.edu"

__all__ = ["main", "logger_settings", "BUILD_CONFIG", "get_build_config"]

import logging

logger = logging.getLogger(__name__)


def _get_version():
    try:
        from setuptools_scm import get_version

        return get_version(root="..", relative_to=__file__)
    except Exception:
        pass

    try:
        from importlib.metadata import version

        return version("orGUI")
    except Exception:
        pass

    try:
        from ._version import version

        return version
    except Exception:
        logger.info("orGUI is not installed. Version number will be incorrect!")
        return "1.4.1-unknown"


try:
    __version__ = _get_version()
except Exception:
    logger.exception("Cannot determine orGUI version")
    __version__ = "1.4.1-unknown"


def get_build_config():
    """Return Meson build settings recorded in the installed package.

    :returns:
        A dictionary with compiler, host, and native optimization settings.
        Source-tree imports before a Meson build return ``{"available": False}``.
    :rtype: dict
    """
    try:
        from ._build_config import BUILD_CONFIG

        return dict(BUILD_CONFIG, available=True)
    except Exception:
        return {"available": False}


BUILD_CONFIG = get_build_config()
