# -*- coding: utf-8 -*-
###############################################################################
# Copyright (c) 2020 Timo Fuchs, Olaf Magnussen all rights reserved
#
# This software was developed during the PhD work of Timo Fuchs,
# within the group of Olaf Magnussen. Usage within the group is hereby granted.
###############################################################################
"""Module descripiton

"""
__author__ = "Timo Fuchs"
__copyright__ = "Copyright 2020, Timo Fuchs, Olaf Magnussen all rights reserved"
__credits__ = []
__license__ = "all rights reserved"
__version__ = "1.0.0"
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

    if os.path.isfile(options.configfile):
        app = qt.QApplication(sys.argv)
        app.setApplicationName("orGUI")
        from .resources import getSplashScreen
        pixmap = getSplashScreen()
        
        desktopWidget = app.desktop()
        screenGeometry = desktopWidget.screenGeometry()
        splashpm = pixmap.scaledToHeight(screenGeometry.height()/3, qt.Qt.SmoothTransformation)
        splash = qt.QSplashScreen(splashpm)
        splash.show()
        splash.showMessage("jit compile libraries", qt.Qt.AlignVCenter | qt.Qt.AlignRight)
        from datautils.xrayutils import CTRcalc, CTRplotutil
        splash.showMessage("load orGUI", qt.Qt.AlignVCenter | qt.Qt.AlignRight)
        from orgui.app.orGUI import orGUI, UncaughtHook
        qt_exception_hook = UncaughtHook()
        mainWindow = orGUI(options.configfile)
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

