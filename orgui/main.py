# -*- coding: utf-8 -*-
# /*##########################################################################
#
# Copyright (c) 2020-2024 Timo Fuchs
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
"""Module descripiton

"""
__author__ = "Timo Fuchs"
__copyright__ = "Copyright 2020-2024 Timo Fuchs"
__credits__ = []
__license__ = "MIT License"
from . import __version__
__maintainer__ = "Timo Fuchs"
__email__ = "fuchs@physik.uni-kiel.de"

import os
import sys
import datetime
from argparse import ArgumentParser

description = """Load th scans of HESXRD reciprocal space mappings from 
single crystal surfaces. Allows maximum image calculation and setting of
an orientation matrix for that scan. The angles for a reciprocal space
coordinate can be subsequently calculated.  
"""

epilog = """The reflection list allows calculation of binned reciprocal
space using the HESXRD backends of binoculars."""

usage = "orGUI [options] configfile"

defaultconfigfile = os.path.expanduser("~/orgui")

def main():
    parser = ArgumentParser(usage=usage, description=description, epilog=epilog)
    parser.add_argument("configfile", metavar="FILE", 
                        help="configuration file, will use ~/orgui otherwise", 
                        nargs='?',default=defaultconfigfile)
    parser.add_argument("--opengl", "--gl", dest="opengl",
                        action="store_true",
                        default=False,
                        help="Enable OpenGL rendering (else matplotlib is used)")

    options = parser.parse_args()
    
    import silx
    from silx.gui import qt
    
    if options.opengl:
        silx.config.DEFAULT_PLOT_BACKEND = "opengl"
    
    if os.path.isfile(options.configfile) or options.configfile == defaultconfigfile:
        app = qt.QApplication(sys.argv)
        app.setApplicationName("orGUI")
        from .resources import getQicon
        app.setWindowIcon(getQicon("orguiicon"))
        from .resources import getSplashScreen
        pixmap = getSplashScreen(__version__)
        desktopWidget = app.desktop()
        screenGeometry = desktopWidget.screenGeometry()
        splashpm = pixmap.scaledToHeight(int(screenGeometry.height()/3), qt.Qt.SmoothTransformation)
        splash = qt.QSplashScreen(splashpm)
        splash.show()
        configfile = options.configfile
        if options.configfile == defaultconfigfile and not os.path.isfile(options.configfile):
            configfile = None
        
        splash.showMessage("jit compile libraries, this may take a while on first start...", qt.Qt.AlignLeft | qt.Qt.AlignBottom)
        from .datautils.xrayutils import CTRcalc, CTRplotutil
        splash.showMessage("load orGUI libraries", qt.Qt.AlignLeft | qt.Qt.AlignBottom)
        splash.raise_()
        from orgui.app.orGUI import orGUI, UncaughtHook
        splash.showMessage("starting orGUI", qt.Qt.AlignLeft | qt.Qt.AlignBottom)
        splash.raise_()
        qt_exception_hook = UncaughtHook()
        mainWindow = orGUI(configfile)
        qt_exception_hook.set_orgui(mainWindow)
        splash.finish(mainWindow)

        mainWindow.show()
        current_screen = mainWindow.window().windowHandle().screen()
        qr = mainWindow.frameGeometry()
        qr.moveCenter(current_screen.geometry().center())
        mainWindow.move(qr.topLeft())
        
        return app.exec()
    else:
        raise Exception("%s is no file" % options.configfile)

