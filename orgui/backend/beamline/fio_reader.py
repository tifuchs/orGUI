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
__version__ = "1.0.0"
__maintainer__ = "Timo Fuchs"
__email__ = "fuchs@physik.uni-kiel.de"


import numpy as np
import os
import warnings
import traceback

dtypeConverter = {'STRING' : 'U32',
                  'DOUBLE' : 'f8',
                  'FLOAT' : 'f4',
                  'INTEGER' : 'i8',
                  'BOOLEAN' : '?'}



class FioFile(object):
    
    def __init__(self,filepath):
        _,filename = os.path.split(filepath)
        fnowithsuffix = filename.split('_')[-1]
        self.scanno = int(fnowithsuffix.split('.')[0])
        
        with open(filepath,'r') as fiof:
            
            prev = 0
            while(True):
                line = fiof.readline()
                if line.startswith('!'):
                    prev = fiof.tell()
                    continue
                if line.startswith('%c'):
                    self.comment = ''
                    line = fiof.readline()
                    while( not line.startswith('%') and not line.startswith('!')):
                        self.comment += line
                        prev = fiof.tell()
                        line = fiof.readline()
                if line.startswith('%p'):
                    self.parameters_str = ''
                    line = fiof.readline()
                    while( not line.startswith('%') and not line.startswith('!')):
                        self.parameters_str += line
                        prev = fiof.tell()
                        line = fiof.readline()
                if line.startswith('%d'):
                    self.datacols = []
                    self.names = []
                    self.dtypes = []
                    line = fiof.readline()
                    while(line.startswith(' Col') ):
                        splitline = line.split()
                        name = splitline[-2]
                        self.names.append(name) 
                        dtype = dtypeConverter[splitline[-1]]
                        self.dtypes.append(dtype)
                        self.datacols.append((name,dtype))
                        prev = fiof.tell()
                        line = fiof.readline()
                    fiof.seek(prev)
                    break
                
            # special cases:
            if "idrz1(encoder)" in self.names:
                idx = self.names.index("idrz1(encoder)")
                self.names[idx] = "idrz1"
            
            self.data = np.loadtxt(fiof,dtype={'names' : tuple(self.names), 'formats' : tuple(self.dtypes)},comments="!")
            #print(self.data)
        self.parameter = {}
        
        # make parameter section nicer, a bit dangerous like this:
        try:
            for line in self.parameters_str.splitlines():
                param , value = line.split(' = ')
                self.parameter[param] = value
        except Exception:
            warnings.warn("Cannot parse parameter section : %s" % traceback.format_exc())
            

        
