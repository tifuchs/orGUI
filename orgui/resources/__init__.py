import os
import glob

from silx.gui import qt

_iconpath = os.path.join(os.path.dirname(__file__), "icons")

def getQicon(name: str):

    return qt.QIcon(os.path.join(_iconpath, name))