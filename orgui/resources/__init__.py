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
import os
import glob

from silx.gui import qt

_iconpath = os.path.join(os.path.dirname(__file__), "icons")

def getQicon(name: str):

    return qt.QIcon(os.path.join(_iconpath, name))
    
def getSplashScreen(version_number=None):
    pixmap = qt.QPixmap(os.path.join(_iconpath, "logo"))
    if version_number is None:
        return pixmap
    else:
        painter = qt.QPainter(pixmap)
        font = qt.QFont()
        font.setPixelSize(30)
        painter.setFont(font)
        painter.drawText(720, 403, "version %s" % version_number)
        return pixmap 
    
def getDiffractometerPath():
    return os.path.join(_iconpath, "diffractometer_v3.png")
    
def getPath(name):
    return os.path.join(_iconpath, name)