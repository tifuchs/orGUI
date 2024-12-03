# -*- coding: utf-8 -*-
# /*##########################################################################
#
# Copyright (c) 2024 Finn Schroeter
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
__author__ = "Finn Schroeter"
__credits__ = ['Finn Schroeter']
__copyright__ = "Copyright 2024 Finn Schroeter"
__license__ = "MIT License"
__version__ = "1.0.0"

import re
import os
import fabio
import numpy as np

from .scans import h5_Image

class ImportImagesScan():

    def __init__(self, imgpath):
        self.filename = imgpath

        self.inpath = self.find_files()

        if self.inpath == None:
            return

        with fabio.open(self.inpath[0] + self.inpath[1][0]) as fabf:
            img_data = fabf.data
            self.FramesPerFile = fabf.nframes
        self.shape = (img_data.shape[0], img_data.shape[1])

        if self.FramesPerFile > 1:
            with fabio.open(self.inpath[0] + self.inpath[1][-1]) as last_file:
                self.FramesLastFile = last_file.nframes
        else:
            self.FramesLastFile = self.FramesPerFile

        self.images = np.zeros((1,*self.shape))

        self.axisname = ''
        self.axis = [0]
        self.th = 0
        self.omega = 0
        self.mu = 0

        self.title = "manually loaded scan" 

        self.nopoints = (len(self.inpath[1])-1)*self.FramesPerFile + self.FramesLastFile


    def find_files(self):
        re_str = re.compile(r'_\d+' + os.path.splitext(self.filename)[1]) #define search string. It matches a file source with syntax name_0000i.extension -> may need to be adapted 
        selected_directory = os.path.dirname(os.path.abspath(self.filename))
        filenames = ''.join(os.listdir(selected_directory))

        found_scanfiles = re_str.findall(filenames) #list of found filenames (suffix only)
        #found_scannrs = [e[1:-4] for e in found_scanfiles] #only the scan numbers, eg. 00015
        if re_str.findall(self.filename) == []:
            print('Could not load images')
            return
        else:
            suffix = re_str.findall(self.filename)[0]
            imagePrefix = self.filename.removesuffix(suffix)
            return [imagePrefix,found_scanfiles]

    def set_axis(self,axismin,axismax,axis,fixedAxisValue):
        self.axis = np.linspace(axismin,axismax,self.nopoints)
        self.axisname = axis
        if axis == 'th':
            self.th = self.axis
            self.omega = -1*self.th
            self.mu = fixedAxisValue

        
    def __len__(self):
        return self.nopoints
        
    def get_raw_img(self, i):
        if self.FramesPerFile > 1:
            index = i // self.FramesPerFile
            frame = i % self.FramesPerFile
            with fabio.open(self.inpath[0] + self.inpath[1][index]) as fabf:
                img_data = fabf.get_frame(frame).data
        else:
            with fabio.open(self.inpath[0] + self.inpath[1][i]) as fabf:
                img_data = fabf.data
        return h5_Image(img_data)
        
    @classmethod
    def parse_h5_node(cls, obj):
        ddict = dict()
        return ddict

    @property
    def auxillary_counters(self):
        return [] 