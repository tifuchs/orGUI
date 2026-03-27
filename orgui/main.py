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
import runpy
import argparse
from argparse import ArgumentParser
import logging 
from . import logger_utils

def existing_file(path):
    if not os.path.isfile(path):
        raise argparse.ArgumentTypeError(f"File '{path}' does not exist")
    return path

def writable_file(path: str):
    """Argparse type: ensure the file can be created or written."""
    try:
        # Try opening the file in append mode; creates it if needed
        with open(path, "a"):
            pass
    except Exception as e:
        raise argparse.ArgumentTypeError(
            f"Cannot write to '{path}': {e.strerror}"
        ) from e
    return path



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
    parser.add_argument("-i", "--input", type=existing_file, 
        help="Path to a batch script file, will be executed with access to the orgui API, orGUI will quit after execution")
    parser.add_argument("--keep-running", dest='keeprunning', action='store_true', default=False,
        help="set to keep orGUI running after completion of the batch script")
    parser.add_argument("--opengl", "--gl", dest="opengl",
                        action="store_true",
                        default=False,
                        help="Enable OpenGL rendering (else matplotlib is used)")
    parser.add_argument("--logfile", type=writable_file, 
                        help="Path for a log file (default level: INFO)")
    parser.add_argument("--errorlog", type=writable_file, 
                        help="Path for a separate error log file")
    
    parser.add_argument("--cpus", type=int, 
                        help="Number of threads used for computation and file reads (will use SLURM_CPUS_ON_NODE or NSLOTS if available, or try to detect otherwise)")
    
    group = parser.add_argument_group('HDF5 file locking', "Sets HDF5_USE_FILE_LOCKING variable. "
                        "Setting to False can resolve some file read issues "
                        "but can also cause file corruption! "
                        "This will not override a manually set environment variable. ")
    locking_parser = group.add_mutually_exclusive_group(required=False)
    locking_parser.add_argument('--hdflocking', '-l', dest='locking', action='store_true', help="HDF5_USE_FILE_LOCKING=True (default)")
    locking_parser.add_argument('--no-hdflocking', '-nl', dest='locking', action='store_false', help="HDF5_USE_FILE_LOCKING=False")
    
    parser.set_defaults(locking=True)

    options = parser.parse_args()
    
    # set up log files first
    
    root_logger = logging.getLogger() # use root logger, i.e. log everything
    
    formatter = logging.Formatter(
        "%(levelname)s:%(asctime)s:%(name)s:%(message)s",
        datefmt="%Y-%m-%dT%H:%M:%S"
    )

    console_stdout = logging.StreamHandler(sys.stdout)
    console_stdout.setLevel(logging.INFO)
    console_stdout.addFilter(lambda r: r.levelno < logging.WARNING)
    console_stdout.setFormatter(formatter)
    root_logger.addHandler(console_stdout)

    console_stderr = logging.StreamHandler(sys.stderr)
    console_stderr.setLevel(logging.WARNING)
    console_stderr.setFormatter(formatter)
    root_logger.addHandler(console_stderr)

    root_logger.setLevel(logging.WARNING)
    
    orgui_logger = logging.getLogger() # use root logger for now, i.e. log everything
    
    orgui_logger.setLevel(logging.INFO)
    
    
    if options.logfile:
        # use root logger, i.e. log everything
        formatter = logging.Formatter(
            "%(levelname)s:%(asctime)s:%(name)s:%(message)s",
            datefmt="%Y-%m-%dT%H:%M:%S"
        )
        logfile_handler = logging.FileHandler(options.logfile)
        logfile_handler.setLevel(logging.INFO)
        logfile_handler.setFormatter(formatter)
        orgui_logger.addHandler(logfile_handler)
    
    if options.errorlog:
        
        formatter = logging.Formatter(
            "%(levelname)s:%(asctime)s:%(name)s:Line:%(lineno)d:%(funcName)s():%(message)s",
            datefmt="%Y-%m-%dT%H:%M:%S"
        )
        errorlog_handler = logging.FileHandler(options.errorlog)
        errorlog_handler.setLevel(logging.ERROR)
        errorlog_handler.setFormatter(formatter)
        orgui_logger.addHandler(errorlog_handler)
    
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
    
    if options.cpus:
        ncpus = options.cpus
    else:
        if 'SLURM_CPUS_ON_NODE' in os.environ:
            ncpus = int(os.environ['SLURM_CPUS_ON_NODE'])
            logger.info('Detected SLURM_CPUS_ON_NODE=%s (on SLURM).' % ncpus)
        elif 'NSLOTS' in os.environ:
            ncpus = int(os.environ['NSLOTS'])
            logger.info('Detected NSLOTS=%s. (on SGE).' % ncpus)
        else:
            ncpus = int(min(os.cpu_count(), 16)) if os.cpu_count() is not None else 1
            logger.info('Detected ncpus=%s. (capped to 16)' % ncpus)
    logger.info('Using ncpus=%s. (capped to 16)' % ncpus)

    
    if options.cli:
        logger_utils.set_logging_context('cli')
        _start_CLI(options, ncpus)
        
    else:
        logger_utils.set_logging_context('gui')
        _start_GUI(options, ncpus)

        
def _start_CLI(options, ncpu):
    
    os.environ["QT_QPA_PLATFORM"] = "minimal" # "offscreen" # maybe use minimal instead
    # os.environ["QT_LOGGING_RULES"] = "*.warning=false"
    
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
        logger.info("loading orGUI")
        mainWindow = orGUI(configfile)
        mainWindow.numberthreads = ncpu
        app.aboutToQuit.connect(mainWindow.database.close)
        
        namespace = {'app': app, 'orgui': mainWindow, 'ub': mainWindow.ubcalc}

        if options.input:
            logger.info("Run batch script %s" % options.input)
            try:
                runpy.run_path(options.input, init_globals=namespace)
            except:
                mainWindow.database.close()
                app.quit()
                raise
                
            logger.info("Completed batch script %s" % options.input)
            if not options.keeprunning:
                logger.info("All tasks completed. Application will now exit.")
                mainWindow.database.close()
                app.quit()
                return
        logger.info("starting orGUI")
        ipshell(local_ns=namespace)
        app.quit()


def _start_GUI(options, ncpu):
    from silx.gui import qt
    import silx
    
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
        splash.showMessage("loading orGUI", qt.Qt.AlignLeft | qt.Qt.AlignBottom)
        logger.info("loading orGUI")
        splash.raise_()
        qt_exception_hook = UncaughtHook()
        mainWindow = orGUI(configfile)
        mainWindow.numberthreads = ncpu
        qt_exception_hook.set_orgui(mainWindow)
        splash.finish(mainWindow)

        mainWindow.show()
        current_screen = mainWindow.window().windowHandle().screen()
        qr = mainWindow.frameGeometry()
        qr.moveCenter(current_screen.geometry().center())
        mainWindow.move(qr.topLeft())
        app.aboutToQuit.connect(mainWindow.database.close)
        namespace = {'app': app, 'orgui': mainWindow, 'ub': mainWindow.ubcalc}
        if options.input:
            logger.info("Run batch script %s" % options.input)
            try:
                runpy.run_path(options.input, init_globals=namespace)
            except:
                raise
            logger.info("Completed batch script %s" % options.input)
            if not options.keeprunning:
                logger.info("All tasks completed. Application will now exit.")
                mainWindow.database.close()
                app.quit()
                return
        logger.info("starting orGUI")
        return app.exec()
    else:
        raise Exception("%s is no file" % options.configfile)
