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


import numpy as np


class h5_Image:
    
    def __init__(self, data):
        
        self.img = data 
        self.motors = dict()
        self.counters = dict()


class SimulationThScan():

    def __init__(self, detshape, ommin, ommax, points):
        
        self.shape = detshape
        
        self.omega = np.linspace(ommin,ommax,points)
        self.th = -1*self.omega
        
        self.nopoints = points
        
        self.images = np.zeros((points,*detshape))
        
    def __len__(self):
        return self.nopoints
        
    def get_raw_img(self, i):
        return h5_Image(self.images[i])
        
    def set_raw_img(self, i, data):
        self.images[i] = data
        
        