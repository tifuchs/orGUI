# -*- coding: utf-8 -*-
###############################################################################
# Copyright (c) 2020 Timo Fuchs all rights reserved
#
# This software was developed during the PhD work of Timo Fuchs,
# within the group of Olaf Magnussen. Usage within the group is hereby granted.
###############################################################################
"""Module descripiton

"""
__author__ = "Timo Fuchs"
__copyright__ = "Copyright 2022, Timo Fuchs all rights reserved"
__credits__ = []
__license__ = "all rights reserved"
__version__ = "1.0.0"
__maintainer__ = "Timo Fuchs"
__email__ = "fuchs@physik.uni-kiel.de"

from silx.gui import qt


def messagebox_detailed_message(parent, title, text, detailed_text, icon, buttons=qt.QMessageBox.Ok):
    diag = qt.QMessageBox(icon, title, text, buttons, parent)
    diag.setDetailedText(detailed_text)
    return diag.exec()

def critical_detailed_message(parent, title, text, detailed_text, buttons=qt.QMessageBox.Ok):
    return messagebox_detailed_message(parent, title, text, detailed_text, qt.QMessageBox.Critical, buttons=buttons)
    
def warning_detailed_message(parent, title, text, detailed_text, buttons=qt.QMessageBox.Ok):
    return messagebox_detailed_message(parent, title, text, detailed_text, qt.QMessageBox.Warning, buttons=buttons)

def information_detailed_message(parent, title, text, detailed_text, buttons=qt.QMessageBox.Ok):
    return messagebox_detailed_message(parent, title, text, detailed_text, qt.QMessageBox.Information, buttons=buttons)




