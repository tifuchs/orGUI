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
        from orgui.app.orGUI import main
        return main(options.configfile)
    else:
        raise Exception("%s is no file" % options.configfile)

