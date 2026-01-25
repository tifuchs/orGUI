# -*- coding: utf-8 -*-
# /*##########################################################################
#
# Copyright (c) 2020-2025 Timo Fuchs
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
__copyright__ = "Copyright 2020-2025 Timo Fuchs"
__credits__ = []
__license__ = "MIT License"
from . import __version__
__maintainer__ = "Timo Fuchs"
__email__ = "tfuchs@cornell.edu"

import os
import sys
import datetime
from argparse import ArgumentParser
import logging 
from . import logger_settings

logger = logging.getLogger(__name__)

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
    parser.add_argument("--nogui", "--cli", "--headless", dest="cli",
                        action="store_true",
                        default=False,
                        help="Start an orGUI IPython shell (no graphical user interface)")
    parser.add_argument("--opengl", "--gl", dest="opengl",
                        action="store_true",
                        default=False,
                        help="Enable OpenGL rendering (else matplotlib is used)")
    group = parser.add_argument_group('HDF5 file locking', "Sets HDF5_USE_FILE_LOCKING variable. "
                        "Setting to False can resolve some file read issues "
                        "but can also cause file corruption! "
                        "This will not override a manually set environment variable. ")
    locking_parser = group.add_mutually_exclusive_group(required=False)
    locking_parser.add_argument('--hdflocking', '-l', dest='locking', action='store_true', help="HDF5_USE_FILE_LOCKING=True (default)")
    locking_parser.add_argument('--no-hdflocking', '-nl', dest='locking', action='store_false', help="HDF5_USE_FILE_LOCKING=False")
    
    parser.set_defaults(locking=True)

    options = parser.parse_args()
    
    if "HDF5_USE_FILE_LOCKING" not in os.environ:
        if options.locking:
            os.environ["HDF5_USE_FILE_LOCKING"] = "TRUE"
        else:
            os.environ["HDF5_USE_FILE_LOCKING"] = "FALSE"
    if str(os.environ["HDF5_USE_FILE_LOCKING"]).lower() in ['y', 'true', '1', 'yes']:
        logger.info("HDF5_USE_FILE_LOCKING=%s" % os.environ["HDF5_USE_FILE_LOCKING"])
        logger.info("Save hdf5 mode, if you encounter hdf5 file read issues, try starting orGUI with the option -nl")
    else:
        logger.info("HDF5_USE_FILE_LOCKING=%s" % os.environ["HDF5_USE_FILE_LOCKING"])
        logger.info("No hdf5 locking. This is potentially dangerous and can cause file corruption. Especially if orGUI crashes.")
    
    if options.cli:
        logger_settings.set_logging_context('cli')
        _start_CLI(options)
    else:
        logger_settings.set_logging_context('gui')
        _start_GUI(options)

        
def _start_CLI(options):
    
    os.environ["QT_QPA_PLATFORM"] = "offscreen" # maybe use minimal instead
    os.environ["QT_LOGGING_RULES"] = "*.warning=false"
    
    if os.path.isfile(options.configfile) or options.configfile == defaultconfigfile:
        from IPython.terminal.embed import InteractiveShellEmbed
        from IPython.terminal.prompts import Prompts, Token
        
        class OrGUIPrompts(Prompts):
            def in_prompt_tokens(self):
                return [
                    (Token.Prompt.Mode, self.vi_mode()),
                    (
                        Token.Prompt.LineNumber,
                        self.shell.prompt_line_number_format.format(
                            line=1, rel_line=-self.current_line()
                        ),
                    ),
                    (Token.Prompt, "orCLI In["),
                    (Token.PromptNum, str(self.shell.execution_count)),
                    (Token.Prompt, ']: '),
                ]

            def out_prompt_tokens(self):
                return [
                    (Token.OutPrompt, '     Out['),
                    (Token.OutPromptNum, str(self.shell.execution_count)),
                    (Token.OutPrompt, ']: '),
                ]

        custom_banner = f"""orCLI {__version__} console  - the command line interface to orGUI
Available variables:
orgui : top level gui
ub : gui for UB matrix and angle calculations 
"""
        ipshell = InteractiveShellEmbed(banner2=custom_banner )
        ipshell.prompts = OrGUIPrompts(ipshell)

        ipshell.enable_gui('qt')
        from silx.gui import qt
        app = qt.QApplication(sys.argv)
        app.setApplicationName("orGUI")
        configfile = options.configfile
        if options.configfile == defaultconfigfile and not os.path.isfile(options.configfile):
            configfile = None
        logger.info("jit compile libraries, this may take a while on first start...")
        from .datautils.xrayutils import CTRcalc, CTRplotutil
        logger.info("load orGUI libraries")
        from orgui.app.orGUI import orGUI, UncaughtHook
        logger.info("starting orGUI")
        mainWindow = orGUI(configfile)
        app.aboutToQuit.connect(mainWindow.database.close)
        ipshell(local_ns={'app': app, 'orgui': mainWindow, 'ub': mainWindow.ubcalc})
        app.quit()


def _start_GUI(options):
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
        screenGeometry = app.primaryScreen().availableGeometry()
        splashpm = pixmap.scaledToHeight(int(screenGeometry.height()/3), qt.Qt.SmoothTransformation)
        splash = qt.QSplashScreen(splashpm)
        splash.show()
        configfile = options.configfile
        if options.configfile == defaultconfigfile and not os.path.isfile(options.configfile):
            configfile = None
        
        splash.showMessage("jit compile libraries, this may take a while on first start...", qt.Qt.AlignLeft | qt.Qt.AlignBottom)
        logger.info("jit compile libraries, this may take a while on first start...")
        from .datautils.xrayutils import CTRcalc, CTRplotutil
        splash.showMessage("load orGUI libraries", qt.Qt.AlignLeft | qt.Qt.AlignBottom)
        logger.info("load orGUI libraries")
        splash.raise_()
        from orgui.app.orGUI import orGUI, UncaughtHook
        splash.showMessage("starting orGUI", qt.Qt.AlignLeft | qt.Qt.AlignBottom)
        logger.info("starting orGUI")
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
        app.aboutToQuit.connect(mainWindow.database.close)
        return app.exec()
    else:
        raise Exception("%s is no file" % options.configfile)