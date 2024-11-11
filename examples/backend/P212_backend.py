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
__author__ = "Timo Fuchs"
__copyright__ = "Copyright 2020-2024 Timo Fuchs"
__license__ = "MIT License"
__version__ = "1.2.0"
__maintainer__ = "Timo Fuchs"
__email__ = "fuchs@physik.uni-kiel.de"


import numpy as np
import datetime
import re
import warnings
import fabio
import h5py
import os
import re
from datetime import datetime
import json

import silx.io
from silx.io import dictdump

from orgui.backend.scans import Scan


            


class Varex_Image:
    def __init__(self, filename, dark=None):
        
        with fabio.open(filename) as img:
            self.img = img.data.astype(np.float64)
            self.header = img.header
        self.motors = dict() # not jet available
        self.counters = dict() # not jet available
        if dark is not None:
            self.img -= dark
            # crude repair of images
            #if np.mean(self.img) < 0:
            #    self.img += 72.2
            
class P3_Image:
    def __init__(self, filename):
        with fabio.open(filename) as img:
            self.img = img.data.astype(np.float64)
            self.header = img.header
        self.motors = dict() # not jet available
        self.counters = dict() # not jet available
        
        
# This is the backend, must implement `Scan`!
class P212_backend_2024(Scan):
    def __init__(self, hdffilepath_orNode=None, scanno=None):
        self.dark = None

        #self.scanname = scanname
        self.cameras = {}
        data = None
        
        if hdffilepath_orNode is None:
            return
        self.hdffilepath_orNode=hdffilepath_orNode
        
        ## --- dump all data for the specified scan into `data` 

        if isinstance(hdffilepath_orNode, str):
            self.hdffilepath_orNode = os.path.abspath(hdffilepath_orNode)
            filepath,filename = os.path.split(self.hdffilepath_orNode)
            filename_noext,_ = os.path.splitext(filepath)
            #filename_noext = filename.split('.')[0]
            self.filename_base = filename_noext #filename_noext[:filename_noext.rfind('_')]
            with silx.io.open(self.hdffilepath_orNode) as f:
                #print([d for d in f])
                for d in f:
                    scanno_s, subscanno = d.split('.')
                    if int(scanno_s) == scanno:
                        break
                else:
                    raise IOError("Scan number %s not found in file" % scanno)
                self.scanname = '.'.join((scanno_s,subscanno))

                data = dictdump.h5todict(f,self.scanname) # here is the actual loading 
        else:
            hdffilepath = hdffilepath_orNode.local_filename
            filepath,filename = os.path.split(hdffilepath)
            filename_noext,_ = os.path.splitext(filepath)
            self.filename_base = filename_noext
            
            f = hdffilepath_orNode.file
            for d in f:
                scanno_s, subscanno = d.split('.')
                if int(scanno_s) == scanno:
                    break
            else:
                raise IOError("Scan number %s not found in file" % scanno)
            self.scanname = '.'.join((scanno_s,subscanno))

            data = dictdump.h5todict(f,self.scanname) # here is the actual loading 
        self.name = self.scanname
        self.title = data['title'][()]
        # find available cameras
        
        possibleCameras = ['Pilatus', 'Varex_1', 'Varex_2', 'Varex_3', 'Varex_4']
        self.defaultCamera = 'Pilatus'
        
        parameter_section = data['instrument']['parameter']
        measurement_section = data['measurement']
        
        for cam in possibleCameras:
            if cam in parameter_section:
                self.cameras[cam] = {}
        
        # parse image file location
        for cam in self.cameras:
            
            self.cameras[cam] = json.loads(parameter_section[cam][()]) # detector configuration
            loc = self.cameras[cam]['Filedir']
            loc = loc.replace('/gpfs/current/', '')
            self.cameras[cam]['Filedir_current'] = os.path.join(self.filename_base , '../..', loc) # assumes is in fio directory
            self.cameras[cam]['imagenumbers'] = np.atleast_1d(measurement_section[cam])
        
        # now all data for the scan in the hdf5 file is in the variable data 
        
        
        # motor conversion table:
        motor_names =     ['idrz1', 'idrx2', 'idry1', 'idry2', 'idtx2', 'idty2', 'idtz2']
        sixc_equivalent = ['omega',   'chi', 'alpha',   'phi',   'say',   'sax',   'saz']
        sixc_sign =       [     -1,       1,      -1,       -1,      1,      -1,       1]
        fastsweep_names = ['(start)', '(end)']
        
        scanned_motors_fast = {}
        for mot in motor_names:
            ext = fastsweep_names[0]
            if mot + ext in measurement_section:
                mot_start = measurement_section[mot + ext]
                ext = fastsweep_names[1]
                mot_end = measurement_section[mot + ext]
                scanned_motors_fast[mot] = (mot_start + mot_end) / 2 # disregard range for now
        
        scanned_motors = {}
        for mot in motor_names:
            if mot in measurement_section:
                scanned_motors[mot] = measurement_section[mot]
        
        for mot, sixc_mot, sign in zip(motor_names, sixc_equivalent, sixc_sign):
            try:
                val = sign * np.asarray(parameter_section[mot][()], dtype=np.float64) # values at start of scan 
                setattr(self, sixc_mot, val)
            except:
                pass # can fail
                
            if mot in scanned_motors:
                val = sign * scanned_motors[mot] # values during scan
                setattr(self, sixc_mot, val)
                
            if mot in scanned_motors_fast:
                val = sign * scanned_motors_fast[mot] # values during continuous scan
                setattr(self, sixc_mot, val)
                
                
        self.time = measurement_section['timestamp']
        if 'dt' in measurement_section:
            self.exposure_time = measurement_section['dt']
        else:
            self.exposure_time = np.gradient(self.time)
            
         
        self.th = -1*self.omega
        self.mu = self.alpha
                    
        if scanned_motors_fast:
            mot = list(scanned_motors_fast.keys())[0] # use only first, this might be changed in the future
            idx = motor_names.index(mot)
            sixc_angle = sixc_equivalent[idx]
            if sixc_angle == 'omega':
                self.axisname = 'th'
            elif sixc_angle == 'alpha':
                self.axisname = 'mu'
            else:
                self.axisname = sixc_angle
        elif scanned_motors:
            mot = list(scanned_motors.keys())[0] # use only first, this might be changed in the future
            idx = motor_names.index(mot)
            sixc_angle = sixc_equivalent[idx]
            if sixc_angle == 'omega':
                self.axisname = 'th'
            elif sixc_angle == 'alpha':
                self.axisname = 'mu'
            else:
                self.axisname = sixc_angle
        else:
            self.axisname = 'time'
            
        self.axis = np.atleast_1d(getattr(self,self.axisname))
        
        self.nopoints = self.axis.size
        
        
    
    @property
    def auxillary_counters(self):
        """Optional: provide a list of counters or motor names, that should 
        be copied into the orGUI data base for further processing.
        
        e.g. return ['exposure_time', 'elapsed_time','time', 'srcur', 'mondio', 'epoch']
        
        after each integration, orGUI will search for these counter names in the Scan
        object and copy the entries into the database.
        """
        return ['beckvolt_1', 'beckvolt_2', 'eh3_entrance', 'oh2_diode1', 'oh2_diode2', 'petracurrent', 'timestamp'] 
    
    @classmethod
    def parse_h5_node(cls, obj):
        """parse scan number from the selected h5node.
        
        obj.basename is the name of the object double cleicked 
        in the NEXUS panel. Usually this includes the scan number.
        
        Here: the nodes are named as in SPEC files, i.e. 1.1, 1.2, 1.3, ... 
        """
        ddict = dict() 
        scanname = obj.basename
        scanno, subscanno = scanname.split('.')
        ddict['scanno'] = int(scanno)
        ddict['name'] = obj.local_name
        return ddict
        
    # returns single image, attempts to use the default camera!
    def get_raw_img(self,img, **kwargs):
        if 'camera' in kwargs:
            cam = kwargs['camera']
            if cam not in self.cameras:
                raise ValueError("No images were recorded with camera %s" % cam)
        elif self.defaultCamera in self.cameras:
            cam = self.defaultCamera
        elif self.cameras:
            cam = list(self.cameras.keys())[0] # use first image in dict (dict is ordered as of python 3.7)
        else:
            raise ValueError("No images were recorded.")
            
        folder = self.cameras[cam]['Filedir_current']
        filenumber = self.cameras[cam]['imagenumbers'][img]
        imagename = self.cameras[cam]['Filepattern'] % int(filenumber)
        imagepath = os.path.join(folder, imagename)
        
        if cam == 'Pilatus':
            return P3_Image(imagepath)
        else:
            return Varex_Image(imagepath)

    
    def __getitem__(self,key):
        return self.get_raw_img(key)
    
    def __len__(self):
        return self.nopoints


