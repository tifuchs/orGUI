# -*- coding: utf-8 -*-
# /*##########################################################################
#
# Copyright (c) 2024 Timo Fuchs
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
__copyright__ = "Copyright 2024 Timo Fuchs"
__license__ = "MIT License"
__version__ = "1.2.0"
__maintainer__ = "Timo Fuchs"
__email__ = "tfuchs@cornell.edu"


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


class P3_Image:
    def __init__(self, filename):
        with fabio.open(filename) as img:
            self.img = img.data.astype(np.float64)
            self.header = img.header
        header = PilatusHeader(self.header['_array_data.header_contents'])
        self.motors = dict() # not jet available
        self.counters = dict() # not jet available
        self.counters['Exposure_time'] = header.header_dict['Exposure_time']
        self.counters['Exposure_period'] = header.header_dict['Exposure_period']
        self.counters['datetime'] = header.get_date_time()
        self.header_dict = header.header_dict
        
# This is the backend, must implement `Scan`!
class QM2_backend_2024(Scan):
    def __init__(self, hdffilepath_orNode=None, scanno=None):
        data = None
        if hdffilepath_orNode is None:
            return
        self.hdffilepath_orNode=hdffilepath_orNode
        
        ## --- dump all data for the specified scan into `data` 

        if isinstance(hdffilepath_orNode, str):
            self.hdffilepath_orNode = os.path.abspath(hdffilepath_orNode)
            self.filepath,self.filename = os.path.split(self.hdffilepath_orNode)
            filename_noext,_ = os.path.splitext(self.filepath)
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

                data = f[self.scanname]
        else:
            hdffilepath = hdffilepath_orNode.local_filename
            self.filepath,self.filename = os.path.split(hdffilepath)
            filename_noext,_ = os.path.splitext(self.filepath)
            self.filename_base = filename_noext
            
            
            f = hdffilepath_orNode.file
            for d in f:
                scanno_s, subscanno = d.split('.')
                if int(scanno_s) == scanno:
                    break
            else:
                raise IOError("Scan number %s not found in file" % scanno)
            self.scanname = '.'.join((scanno_s,subscanno))

            data = f[self.scanname]
            
        self.name = "%s.1" % scanno
        self.title = data['title'][()]
        
        positioners_section = data['instrument']['positioners']
        measurement_section = data['measurement']
        

        # motor conversion table:
        motor_names =     ['phi2', 'chi', 'th', 'phi']
        sixc_equivalent = ['omega',   'chi', 'alpha', 'omega']
        sixc_sign =       [     -1,       1,      1,  1.]
        sxic_offset =     [      0,    -270,      0,   0]
        
        
        scanned_motors = {}
        for mot in motor_names:
            if mot in measurement_section:
                scanned_motors[mot] = measurement_section[mot]
        
        for mot, sixc_mot, sign, offset in zip(motor_names, sixc_equivalent, sixc_sign, sxic_offset):
            try:
                val = sign * np.asarray(positioners_section[mot][()] + offset, dtype=np.float64) # values at start of scan 
                setattr(self, sixc_mot, val)
            except:
                pass # can fail
                
            if mot in scanned_motors:
                val = sign * scanned_motors[mot] + offset # values during scan
                setattr(self, sixc_mot, val)
                
        
        if 'Time' in measurement_section:
            self.time = measurement_section['Time'][()]
        else:
            self.time = np.cumsum(measurement_section['Seconds'][()]) # this is just an approximation

        self.th = -1*self.omega
        self.mu = self.alpha
                    
        if scanned_motors:
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

        
        for c in self.auxillary_counters:
            try:
                val = np.asarray(measurement_section[c][()], dtype=np.float64) # values at start of scan 
                setattr(self, c, val)
            except:
                try:
                    val = np.asarray(positioners_section[c][()], dtype=np.float64) # values at start of scan 
                    setattr(self, c, val)
                except:
                    pass
        
        # search for scan folder in raw6M
        
        foldername_scan = self.filename + ("_%03i" % scanno)
        root = os.path.join(self.filepath, 'raw6M', self.filename)

        if os.path.isdir(root):
            #from IPython import embed; embed()
            subfolders = [f.path for f in os.scandir(root) if f.is_dir()]
            sssfolders = []
            for sf in subfolders:
                for ssf in os.scandir(sf):
                    if ssf.is_dir():
                        sssfolders.extend([f.path for f in os.scandir(ssf) if f.is_dir()])
            
            for sssf in sssfolders:
                directory, fname = os.path.split(sssf)
                if fname == foldername_scan:
                    break
            else: # search in tiff folder instead
                root_tiff = os.path.join(self.filepath, 'tiffs')
                subfolders = [f.path for f in os.scandir(root_tiff) if f.is_dir()]
                for sssf in subfolders:
                    directory, fname = os.path.split(sssf)
                    if fname == foldername_scan:
                        break
                else:
                    raise ValueError("Cannot find image folder '%s' for scan number %s" % (foldername_scan, scanno))
        else:
            root_tiff = os.path.join(self.filepath, 'tiffs')
            subfolders = [f.path for f in os.scandir(root_tiff) if f.is_dir()]
            for sssf in subfolders:
                directory, fname = os.path.split(sssf)
                if fname == foldername_scan:
                    break
            else:
                raise ValueError("Cannot find image folder '%s' for scan number %s" % (foldername_scan, scanno))
        #from IPython import embed; embed()
        self.imagefolder = os.path.abspath(sssf)
        self.image_prefix = self.filename + ("_PIL10_%03i_" % scanno)
        self.image_suffix = ".cbf"
        
        
    
    @property
    def auxillary_counters(self):
        """Optional: provide a list of counters or motor names, that should 
        be copied into the orGUI data base for further processing.
        
        e.g. return ['exposure_time', 'elapsed_time','time', 'srcur', 'mondio', 'epoch']
        
        after each integration, orGUI will search for these counter names in the Scan
        object and copy the entries into the database.
        """
        return ['diode', 'ic1', 'ic2', 'emon', 'nemon', 'pemon'] 
    
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
        imagepath = os.path.join(self.imagefolder, self.image_prefix + ("%05i" % img) + self.image_suffix)
        try:
            return P3_Image(imagepath)
        except IOError:
            imagepath = os.path.join(self.imagefolder, self.image_prefix + ("%03i" % img) + self.image_suffix)
            return P3_Image(imagepath)
        


    
    def __getitem__(self,key):
        return self.get_raw_img(key)
    
    def __len__(self):
        return self.nopoints


#
# ==========================================================================
# == pilatus_header.py
# == Parsing detector and experiment parameters from Pilatus file headers
# == 19.09.2011
# == Marcus Mueller (marcus.mueller@dectris.com)
# == Dectris Ltd.
#
# ==========================================================================

SPACE_EQUIVALENT_CHARACTERS = '#:=,()'
NON_OPTIONAL_KEYWORDS = {
#key: (pattern, [value_indeces], type)
    'Detector_identifier': ('Detector ', [slice(1, None)], str),
    'Pixel_size': ('Pixel_size', [1, 4], float),
    'Silicon': ('Silicon', [3], float),
    'Exposure_time': ('Exposure_time', [1], float),
    'Exposure_period': ('Exposure_period', [1], float),
    'Tau': ('Tau', [1], float),
    'Count_cutoff': ('Count_cutoff', [1], int),
    'Threshold_setting': ('Threshold_setting', [1], float),
    'Gain_setting': ('Gain_setting', [1, 2], str),
    'N_excluded_pixels': ('N_excluded_pixels', [1], int),
    'Excluded_pixels': ('Excluded_pixels', [1], str),
    'Flat_field': ('Flat_field', [1], str),
    'Trim_file': ('Trim_file', [1], str),
    'Image_path': ('Image_path', [1], str),
}
OPTIONAL_KEYWORDS = {
    'Wavelength': ('Wavelength', [1], float),
    'Energy_range': ('Energy_range', [1, 2], float),
    'Detector_distance': ('Detector_distance', [1], float),
    'Detector_Voffset': ('Detector_Voffset', [1], float),
    'Beam_xy': ('Beam_xy', [1, 2], float),
    'Beam_x': ('Beam_xy', [1], float),
    'Beam_y': ('Beam_xy', [2], float),
    'Flux': ('Flux', [1], str),
    'Filter_transmission': ('Filter_transmission', [1], float),
    'Start_angle': ('Start_angle', [1], float),
    'Angle_increment': ('Angle_increment', [1], float),
    'Detector_2theta': ('Detector_2theta', [1], float),
    'Polarization': ('Polarization', [1], float),
    'Alpha': ('Alpha', [1], float),
    'Kappa': ('Kappa', [1], float),
    'Phi': ('Phi', [1], float),
    'Phi_increment': ('Phi_increment', [1], float),
    'Chi': ('Chi', [1], float),
    'Chi_increment': ('Chi_increment', [1], float),
    'Oscillation_axis': ('Oscillation_axis', [slice(1, None)], str),
    'N_oscillations': ('N_oscillations', [1], int),
    'Start_position': ('Start_position', [1], float),
    'Position_increment': ('Position_increment', [1], float),
    'Shutter_time': ('Shutter_time', [1], float),
}
ALL_KEYWORDS = {}
ALL_KEYWORDS.update(NON_OPTIONAL_KEYWORDS)
ALL_KEYWORDS.update(OPTIONAL_KEYWORDS)

class PilatusHeader(object):
    """
    Class for parsing contents of a Pilatus cbf header from a minimal
    or full cbf file.
    Parsing the Pilatus cbf header populates the header_dict dictionary,
    using Pilatus cbf header keywords as keys.
    """
    def __init__(self, rawheader):
        self.rawheader = rawheader
        self.header_dict = {}
        self.header_lines = rawheader.splitlines()
        #self.read_header_lines()
        self.parse_header()
        
    def _rawheader(self):
        return self.rawheader
    
    def has_pilatus_cbf_convention(self):
        # Check for the _array_data.header_convention data item
        pilatus_header_pattern = re.compile(
            r'''_array_data.header_convention +["']?(SLS|PILATUS)'''
            r'''_\d+(\.?\d*)*["']?''')
        return bool(pilatus_header_pattern.search(self._rawheader()))
    
    def read_header_lines(self):
        """
        Populate the self.header_lines list.
        """
        contents_pattern = re.compile(
            r'''_array_data.header_contents\s+;.*?;''',
            re.DOTALL)
        contents_match = contents_pattern.search(self._rawheader())
        assert contents_match is not None
        self.header_lines = contents_match.group().splitlines()
        
    def _spaced_header_lines(self):
        """
        Return header_lines with all space equivalent charecters converted
        to space.
        """
        spaced_header_lines = []
        for line in self.header_lines:
            for space_equivalent in SPACE_EQUIVALENT_CHARACTERS:
                line = line.replace(space_equivalent, ' ')
            spaced_header_lines.append(line)
        return spaced_header_lines
    
    def parse_header(self):
        """
        Populate self.header_dict with contents of Pilatus cbf header
        """
        #assert self.has_pilatus_cbf_convention()
        if len(self.header_lines) == 0:
            self.read_header_lines()
        # parse the header lines
        for key, (pattern, valueindices, datatype) in ALL_KEYWORDS.items():
            for line in self._spaced_header_lines():
                if pattern in line:
                    values = []
                    for idx in valueindices:
                        try:
                        # handle multiple or single values
                            if isinstance(idx, slice):
                                values += line.split()[idx]
                            else:
                                values.append(line.split()[idx])
                        except IndexError:
                            print ('No value at index %d on header line:\n'
                                '%s' % (idx, line))
                    value = self._datatype_handling(values, key, datatype)
                    if value is not None:
                        self.header_dict[key] = value
    
    def _datatype_handling(self, values, key, datatype):
        # handle rare cases of value "not set"
        if datatype is float and values[0] == 'not':
            # NON_OPTIONAL_KEYWORDS should always have value, at least NaN
            if key in NON_OPTIONAL_KEYWORDS:
                return float('NaN')
            else:
                return None
        # do the conversion for standard cases
        if len(values) == 1:
            values = datatype(values[0])
        else:
            if datatype is str:
                values = ' '.join(values)
            else:
                values = tuple([datatype(v) for v in values])
        return values
    
    def get_beam_xy_mm(self, factor=1000):
        return tuple([n * size * factor for n, size in zip(
            self.header_dict['Beam_xy'], self.header_dict['Pixel_size'])])
    
    def get_date_time(self):
        """
        Return date and time of image acquistion.
        Works for format of current camserver versions
        2011-06-04T04:57:02.976
        or format of old camserver versions
        2011/Sep/12 09:21:27.252
        """
        date_time_pattern = re.compile(
            r'(\d{4}-\d{2}-\d{2}T|\d{4}/\D+/\d{2} )\d{2}:\d{2}:\d{2}.\d+')
        return date_time_pattern.search(self._rawheader()).group()
    
    def get_time(self):
        time_pattern = re.compile(r'\d{2}:\d{2}:\d{2}.\d+')
        return time_pattern.search(self.get_date_time()).group()
    
    date_time = property(get_date_time)
    time = property(get_time)
        
