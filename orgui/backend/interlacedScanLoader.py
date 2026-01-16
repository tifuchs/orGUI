# -*- coding: utf-8 -*-
# /*##########################################################################
#
# Copyright (c) 2024-2025 Finn Schroeter
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
__copyright__ = "Copyright 2024-2025 Finn Schroeter"
__license__ = "MIT License"
__version__ = "1.3.0"

import re
import os
import fabio
import numpy as np

from .scans import h5_Image

class InterlacedScan():

    def __init__(self, scansegments):
        self.subscans = scansegments

        scan0 = self.subscans[0]
        lensum = 0
        th_area = []

        for i in scansegments:
            lensum += i.nopoints
            th_area.append(i.th)

        self.axisname = scan0.axisname # todo: check if all scans share the same axis
        
        self.th = np.concatenate(th_area) # todo: add option to sort th entries
        self.axis = np.concatenate(th_area)
        self.omega = -1*self.th
        self.title = "interlaced scan"

        self.mu = scan0.mu # todo: check if all scans share the same mu

        self.nopoints = lensum

    def get_raw_img(self, i):
        len_previous = 0
        for k in self.subscans:
            if i < k.nopoints+len_previous:
                return k.get_raw_img(i-len_previous)
            else:
                len_previous += k.nopoints

    @property
    def auxillary_counters(self):
        #todo: combine aux counters of all scan segments 
        #return ['current', 'potential', 'exposure_time', 'elapsed_time','time', 'srcur', 'mondio', 'epoch','scaled_potv2f'] 
        return [] 

    def __len__(self):
        return self.nopoints
        
    @classmethod
    def parse_h5_node(cls, obj):
        ddict = dict()
