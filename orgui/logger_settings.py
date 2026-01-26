# -*- coding: utf-8 -*-
# /*##########################################################################
#
# Copyright (c) 2026 Timo Fuchs
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
__copyright__ = "Copyright 2026 Timo Fuchs"
__credits__ = []
__license__ = "MIT License"
__maintainer__ = "Timo Fuchs"
__email__ = "tfuchs@cornell.edu"

import os
import sys
import logging

from . import __version__

from silx.gui import qt
from silx.gui import icons
import numpy as np
import traceback


# QMetaObject.invokeMethod(
#             QApplication.instance(),
#             lambda: show_messagebox(level, msg),
#             Qt.QueuedConnection
#         )


def messagebox_detailed_message(parent, title, text, detailed_text, icon, buttons=qt.QMessageBox.Ok):
    diag = qt.QMessageBox(icon, title, text, buttons, parent)
    if detailed_text != "":
        diag.setDetailedText(detailed_text)
    return diag.exec()

def critical_detailed_message(parent, title, text, detailed_text, buttons=qt.QMessageBox.Ok):
    return messagebox_detailed_message(parent, title, text, detailed_text, qt.QMessageBox.Critical, buttons=buttons)
    
def warning_detailed_message(parent, title, text, detailed_text, buttons=qt.QMessageBox.Ok):
    return messagebox_detailed_message(parent, title, text, detailed_text, qt.QMessageBox.Warning, buttons=buttons)

def information_detailed_message(parent, title, text, detailed_text, buttons=qt.QMessageBox.Ok):
    return messagebox_detailed_message(parent, title, text, detailed_text, qt.QMessageBox.Information, buttons=buttons)


logger = logging.getLogger("orgui")
logger.setLevel(logging.INFO)

_LOGGING_CONTEXT = 'None'

def get_logging_context():
    global _LOGGING_CONTEXT
    return _LOGGING_CONTEXT

def set_logging_context(context):
    global _LOGGING_CONTEXT, logger
    if context == _LOGGING_CONTEXT:
        return
    for h in logger.handlers[:]:
        logger.removeHandler(h)
        h.close()
    
    if context.lower() == 'gui':
        handler = MessageBoxHandler()
        formatter = logging.Formatter("%(asctime)s: %(levelname)s: %(message)s")
        handler.setFormatter(formatter)
        logger.addHandler(handler)
        _LOGGING_CONTEXT = 'gui'
    elif context.lower() == 'cli':
        handler = CLIExceptionHandler()
        formatter = logging.Formatter("%(asctime)s: %(levelname)s: %(message)s")
        handler.setFormatter(formatter)
        logger.addHandler(handler)
        _LOGGING_CONTEXT = 'cli'
    else:
        raise ValueError('Logging context %s is unknown.' % context)
        

class CLIExceptionHandler(logging.Handler):
    def emit(self, record):
        if record.levelno < logging.ERROR:
            return # super().emit(record)
        msg = self.format(record)

        if record.levelno >= logging.ERROR:
            if record.exc_info:
                exc_type, exc_value, exc_tb = record.exc_info
                raise exc_value.with_traceback(exc_tb)
            else:
                raise RuntimeError(msg)

    
class MessageBoxHandler(logging.Handler):
    def emit(self, record):
        show_dialog = getattr(record, "show_dialog", False) # only show dialog when explicitly requested
        if not show_dialog:
            return
        msg = self.format(record)
        defaulttitle =  f"{record.levelname} [{record.module}:{record.funcName}:{record.lineno}]"
        title = getattr(record, "title", defaulttitle)
        description = getattr(record, "description", "")
        detailed_text = getattr(record, "detailed_text", "")
        dialog_level = getattr(record, "dialog_level", record.levelno)
        show_dialog = getattr(record, "show_dialog", False)
        parent = getattr(record, "parent", None)
        if description != "":
            msg += '\n' + description
            
        if dialog_level >= logging.ERROR:
            if record.exc_info:
                exc_text = "".join(traceback.format_exception(*record.exc_info))
                detailed_text += exc_text
            critical_detailed_message(parent, title, msg, detailed_text)
        elif dialog_level >= logging.WARNING:
            if record.exc_info:
                exc_text = "".join(traceback.format_exception(*record.exc_info))
                detailed_text += exc_text
            warning_detailed_message(parent, title, msg, detailed_text)
        elif dialog_level >= logging.INFO:
            if record.exc_info:
                exc_text = "".join(traceback.format_exception(*record.exc_info))
                detailed_text += exc_text
            information_detailed_message(parent, title, msg, detailed_text)





