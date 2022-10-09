# -*- coding: utf-8 -*-
###############################################################################
# Copyright (c) 2021 Timo Fuchs, Olaf Magnussen all rights reserved
#
# This software was developed during the PhD work of Timo Fuchs,
# within the group of Olaf Magnussen. Usage within the group is hereby granted.
###############################################################################
"""Module descripiton

"""
__author__ = "Timo Fuchs"
__copyright__ = "Copyright 2021, Timo Fuchs, Olaf Magnussen all rights reserved"
__credits__ = []
__license__ = "all rights reserved"
__version__ = "1.0.0"
__maintainer__ = "Timo Fuchs"
__email__ = "fuchs@physik.uni-kiel.de"

import numpy as np


class Scan():

    def __init__(self,hdffilepath_orNode=None, scanno=None):
        self.axisname = None # either "th" or "mu"
        self.axis = None # value of either "th" or "mu"
        self.th = 0. # or self.mu, depending on scanaxis
        self.omega = 0. # = -1*th
        self.title = "generic_scan"
        # for mu-scan you must provide a value for omega/theta 

    def __len__(self):
        raise NotImplementedError()

    def get_raw_img(self, i):
        raise NotImplementedError()
        


class h5_Image:
    
    def __init__(self, data):
        
        self.img = data 
        self.motors = dict()
        self.counters = dict()


class SimulationScan(Scan):

    def __init__(self, detshape, axismin, axismax, points, axis='th', fixed=0.):
        
        self.shape = detshape
        self.axisname = axis
        self.axis = np.linspace(axismin,axismax,points)
        if axis == 'th':
            self.th = self.axis
            self.omega = -1*self.th
            self.mu = fixed
        elif axis == 'mu':
            self.mu = self.axis
            self.th = fixed
            self.omega = -1*self.th
        else:
            raise ValueError("%s is not an implemented scan axis." % axis)
        self.nopoints = points
        
        self.images = np.zeros((points,*detshape))
        self.title = "sim ascan %s %s %s %s" % (self.axisname,axismin,axismax,points)
        
    def __len__(self):
        return self.nopoints
        
    def get_raw_img(self, i):
        return h5_Image(self.images[i])
        
    def set_raw_img(self, i, data): #for intensity simulation in the future.
        self.images[i] = data
        
        

