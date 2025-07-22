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
__author__ = "Timo Fuchs"
__copyright__ = "Copyright 2020-2025 Timo Fuchs"
__license__ = "MIT License"
__version__ = "1.3.0"
__maintainer__ = "Timo Fuchs"
__email__ = "tfuchs@cornell.edu"


import numpy as np
from abc import ABC, abstractmethod

class Scan(ABC):
    """Base class for every backend.

    All methods in this class must be populated for a working backend. 
    """
    
    @abstractmethod
    def __init__(self,hdffilepath_orNode=None, scanno=None):
        """Constructor of the class. 
        
        hdffilepath_orNode is an either an object selected from the file browser,
        where data scan be accessed, or a file path. Example on how to open 
        
        .. code-block:: python
            if isinstance(hdffilepath_orNode, str): # check if this is a file path or already an opened h5node
                with silx.io.open(hdffilepath_orNode) as f: # file path, so we have to open the file manually
                    self.th = f[scanno ,'th'] # read data from the file, see below for the required entries
                    ...
            else:       # file already open, so just read the data from the file
                hdffilepath = hdffilepath_orNode.local_filename # to get the filename
                f = hdffilepath_orNode.file
                self.th = f[scanno ,'th']
                ...
        
        A backend must populate the following entries
        *  axisname : 'th' or 'mu', defines the scan axis
        *  axis : data value of the axis (should be a numpy array)
        *  title : a meaningful title of the scan (e.g. "ascan th 0 90 90 1")
        *  th : motor value (as array)
        *  om : = -th (just populate the value)
        *  mu : motor value
        with actual data 
        """
        self.axisname = None # either "th" or "mu", this defines the scan axis
        self.axis = None # value of either "th" or "mu"
        
        self.th = 0. # or self.mu, depending on scanaxis
        self.omega = -self.th # = -1*th
        self.title = "generic_scan"
        # for mu-scan you must provide a value for omega/theta
        
        # optional: provide a unique identifier, which is used as key in h5 database to identify 
        # the dataset
        # self.name = "identifier"
    
    @property
    def auxillary_counters(self):
        """Optional: provide a list of counters or motor names, that should 
        be copied into the orGUI data base for further processing.
        
        e.g. return ['exposure_time', 'elapsed_time','time', 'srcur', 'mondio', 'epoch']
        
        after each integration, orGUI will search for these counter names in the Scan
        object and copy the entries into the database.
        """
        return [] 

    @classmethod    
    @abstractmethod
    def parse_h5_node(cls, node):
        pass

    @abstractmethod
    def __len__(self):
        """returns the number of entries in the scan.
        """
        raise NotImplementedError()

    @abstractmethod
    def get_raw_img(self, i):
        """This should return a populated image class such as h5_Image (see below)
        
        Only the image data is required as numpy array, accessible as h5_image.img
        """
        raise NotImplementedError()
        


class h5_Image:
    
    def __init__(self, data):
        """Only the image data is required as numpy array.
        motors and counters do not need to be populated.
        """
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
    
    def parse_h5_node(cls, node): # unused
        pass
        

