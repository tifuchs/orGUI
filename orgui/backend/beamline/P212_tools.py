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
from PyMca5.PyMcaIO import EdfFile
from PyMca5.PyMcaIO import PilatusCBF
from PyMca5.PyMcaIO import specfile
import scipy.ndimage
import scipy.interpolate
import scipy.optimize
from matplotlib.colors import LinearSegmentedColormap
import datetime
import re
sqrt3 = np.sqrt(3)
import warnings
import fabio
import h5py
import os
import re
from datetime import datetime
from .fio_reader import FioFile
#from fio_reader import FioFile
import silx.io
from silx.io import dictdump

# not in use!!!
class P3_Image:
    def __init__(self, filename,edf=False):
        if edf:
            image = fabio.open(filename)
            #edf = EdfFile.EdfFile(filename, 'r')
            #self.header = edf.GetHeader(0)
            self.header = image.header
            #self.header = edf.GetHeader(0)
            self.motors = dict()
            motor_mne = self.header['motor_mne'].split()
            motor_pos = self.header['motor_pos'].split()
            for i in range(len(motor_mne)):
                self.motors[motor_mne[i]] = float(motor_pos[i])
            self.counters = dict()
            counter_mne = self.header['counter_mne'].split()
            counter_pos = self.header['counter_pos'].split()
            for i in range(len(counter_mne)):
                self.counters[counter_mne[i]] = float(counter_pos[i])
            self.img = image.data.astype(np.float64)
            #self.img = np.array(edf.GetData(0), dtype=float)
        else:
            try:
                image = fabio.open(filename)
                self.header = image.header
                self.img = image.data.astype(np.float64)
            except Exception as e:
                warnings.warn("Cannot open image, ok only during test run!!! %s" % e)
                self.header = ""
                #self.img = np.zeros((2880, 2880))


            

            self.motors = dict() # not available in cbf
            self.counters = dict() # not available in cbf

            


class PE_Image:
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
        self.img = np.fliplr(self.img.T)

class H5Fastsweep(object):
    def __init__(self, hdffilepath_orNode=None, scanno=None):
    

        #self.scanname = scanname
        self.cameras = []
        data = None
        
        if hdffilepath_orNode is None:
            return
        self.hdffilepath_orNode=hdffilepath_orNode    

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

                data = dictdump.h5todict(f,self.scanname)
                
                    
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

            data = dictdump.h5todict(f,self.scanname)
    
        self.th = data['measurement']['idrz1(encoder)']
        dth = np.mean(np.diff(self.th))
        self.th += dth
        self.th = self.th[1:-1] #last image is missing, first image is dark file
        
        self.th -= 180.
        
        #meandt = np.mean(np.diff(datafile.data['time']))
        self.exposure_time = np.full_like(self.th,1.) #np.diff(datafile.data['time'],append=(datafile.data['time'][-1] + meandt))
        self.time = data['measurement']['time'][1:-1]
        #self.th = -1*self.th  # orientation of th is inversed
        self.omega = -1*self.th
        
        
        
        #self.exposure_time = np.full_like(self.th,self.exposure_time)
        self.directFolder = True
        self.imagenames = data['measurement']['filename']
        #from IPython import embed; embed()
        parsedscanname = re.compile('=\"([^\"]*)\"').findall(str(data['title']))
        
        directory = parsedscanname[1]
        _,folder = os.path.split(directory)
        imagefolder = os.path.join(os.path.split(filepath)[0],folder)
        
        self.filenames = [os.path.join(imagefolder,imgname) for imgname in self.imagenames ]
        with fabio.open(self.filenames[0]) as fiof:
            self.dark = fiof.data
        
        self.filenames = self.filenames[1:-1]
        #self.set_image_folder(self.path+ "/" + detector )
        self.nopoints = self.th.size
        self.axisname = 'th'
        self.axis = getattr(self,self.axisname)
   

        
    def set_th_offset(self,offset):
        self.th += offset
        self.omega = -1*self.th

    """
    def set_image_folder(self,path_to_folder):
        #self.filenames = [None]*len(self.th)
        self.filenames = []
        basefolder = path_to_folder.split('/')[-1]
        self.imageprefix = basefolder.split('_')[0]
        if self.directFolder:
            for i in range(self.th.size):
                self.filenames.append("%s/%s_%05i.%s" % (path_to_folder,\
                    self.imageprefix,i+self.first_imageno,self.imagesuffix))
        else:
            for i in range(self.th.size):
                self.filenames.append("%s/%s/%s_%05.cbf" % (path_to_folder,\
                    self.filename_base,self.filename_base,self.first_imageno+i))
    """        
    def get_raw_P3_img(self,img):
        return PE_Image(self.filenames[img],self.dark)
        
    # returns single default image!
    def get_raw_img(self,img):
        return PE_Image(self.filenames[img],self.dark)
    
    def get_P3_img(self,img):
        imgdata = PE_Image(self.filenames[img],self.dark)
        imgdata.counters['Time'] = self.time[img]
        imgdata.counters['TrigTime'] = self.exposure_time[img]
        imgdata.counters['th'] = self.th[img]
        imgdata.counters['om'] = self.omega[img]
        if self.current is not None:
            imgdata.counters['srcur'] = self.current[img]
        return imgdata
    
    def __getitem__(self,key):
        return self.get_P3_img(key)
    
    def __len__(self):
        return self.nopoints


# only log file existing
class XRRScan(object):
    def __init__(self,logfilepath,img_prefix="xx",scanno=None):
        self.path,fname = os.path.split(logfilepath)
        #fnamewotsuff, suffix = fname.split('.')
        #prefix, nostr = fnamewotsuff.split('_')
        #if scanno is None:
        #    scanno = int(nostr)
        #fiofilecorrectedscanno = "%s/%s_%05i.%s" % (self.path,prefix,scanno,suffix) 
        #self.dtime = []
        #self.current = [] # sr current or some other monitor is missing!
        with open(logfilepath,'r') as logfile:
            line = logfile.readline()
            self.header = []
            while line.startswith('#'): # in header
                self.header.append(line.strip('#').strip())
                lineno = logfile.tell()
                line = logfile.readline()
            logfile.seek(lineno)
            data = np.genfromtxt(logfile)
        self.mu = -1*data[:,1]
        self.th = np.full_like( self.mu,float( self.header[4].split()[-1] ))
            
        
        #datafile = FioFile(fiofilecorrectedscanno)
        #self.th = datafile.data['idrz1']
        # ********** temporary fix *************************************
        self.th -= 80.0
        # **************************************************************
        #imgfilepath, img_foldername = os.path.split(img_prefix)
        #imagefoldernumber = int(img_foldername.split('_')[0])
        #voltagelog = np.genfromtxt("%s/%i_voltage.log" % (imgfilepath,imagefoldernumber))
        # file not saved !!! shit!
        
        #meandt = np.mean(np.diff(datafile.data['time']))
        
        self.exposure_time = np.full_like( self.mu,0.5) # for all XRRs in Dec 2019 beamtime
        
        self.time = data[:,0]
        #self.th = -1*self.th  # orientation of th is inversed
        self.omega = -1*self.th
        self.current = data[:,2]
        
        
        #self.exposure_time = np.full_like(self.th,self.exposure_time)
        self.directFolder = True
        #self.imagenames = datafile.data['filename']

        self.imageprefix,fmiddle = img_prefix, "00001.tif"
        self.imagesuffix = fmiddle.split('.')[-1]
        self.first_imageno = 1
        
        #self.set_image_folder(self.path+ "/" + detector )
        self.nopoints = self.th.size
        
        self.set_image_folder(img_prefix)

    def set_th_offset(self,offset):
        self.th += offset
        self.omega = -1*self.th
        
    def set_image_folder(self,path_to_folder):
        #self.filenames = [None]*len(self.th)
        self.filenames = []
        basefolder = path_to_folder.split('/')[-1]
        self.imageprefix = basefolder.split('_')[0]
        if self.directFolder:
            for i in range(self.th.size):
                self.filenames.append("%s/%s_%05i.%s" % (path_to_folder,\
                    self.imageprefix,i+self.first_imageno,self.imagesuffix))
        else:
            for i in range(self.th.size):
                self.filenames.append("%s/%s/%s_%05.cbf" % (path_to_folder,\
                    self.filename_base,self.filename_base,self.first_imageno+i))
            
    def get_raw_P3_img(self,img):
        return P3_Image(self.filenames[img])
        
    # returns single default image!
    def get_raw_img(self,img):
        return P3_Image(self.filenames[img])
    
    def get_P3_img(self,img):
        imgdata = P3_Image(self.filenames[img])
        imgdata.counters['Time'] = self.time[img]
        imgdata.counters['TrigTime'] = self.exposure_time[img]
        imgdata.counters['th'] = self.th[img]
        imgdata.counters['om'] = self.omega[img]
        if self.current is not None:
            imgdata.counters['srcur'] = self.current[img]
        if self.mu is not None:
            imgdata.counters['mu'] = self.mu[img]
        return imgdata
    
    def __getitem__(self,key):
        return self.get_P3_img(key)
    
    def __len__(self):
        return self.nopoints    

      
class FioFastsweep(object):
    def __init__(self,fiofilepath,img_prefix="xx",scanno=None):
        self.path,fname = os.path.split(fiofilepath)
        fnamewotsuff, suffix = fname.split('.')
        prefix, nostr = fnamewotsuff.split('_')
        if scanno is None:
            scanno = int(nostr)
        fiofilecorrectedscanno = "%s/%s_%05i.%s" % (self.path,prefix,scanno,suffix) 
        #self.dtime = []
        #self.current = [] # sr current or some other monitor is missing!
        datafile = FioFile(fiofilecorrectedscanno)
        self.th = datafile.data['idrz1']
        # ********** temporary fix *************************************
        self.th -= 80.0
        # **************************************************************
        #imgfilepath, img_foldername = os.path.split(img_prefix)
        #imagefoldernumber = int(img_foldername.split('_')[0])
        #voltagelog = np.genfromtxt("%s/%i_voltage.log" % (imgfilepath,imagefoldernumber))
        # file not saved !!! shit!
        
        meandt = np.mean(np.diff(datafile.data['time']))
        self.exposure_time = np.diff(datafile.data['time'],append=(datafile.data['time'][-1] + meandt))
        self.time = datafile.data['time']
        #self.th = -1*self.th  # orientation of th is inversed
        self.omega = -1*self.th
        
        
        
        #self.exposure_time = np.full_like(self.th,self.exposure_time)
        self.directFolder = True
        self.imagenames = datafile.data['filename']
        try:
            self.imageprefix,fmiddle = self.imagenames[0].split('_')
            self.imagesuffix = fmiddle.split('.')[-1]
            print(self.imageprefix)
        except:
            self.imageprefix,fmiddle = img_prefix, "00001.tif"
            self.imagesuffix = fmiddle.split('.')[-1]
        self.first_imageno = 1
        self.current = None
        #self.set_image_folder(self.path+ "/" + detector )
        self.nopoints = self.th.size
        
        
   

        
    def set_th_offset(self,offset):
        self.th += offset
        self.omega = -1*self.th

        
    def set_image_folder(self,path_to_folder):
        #self.filenames = [None]*len(self.th)
        self.filenames = []
        basefolder = path_to_folder.split('/')[-1]
        self.imageprefix = basefolder.split('_')[0]
        if self.directFolder:
            for i in range(self.th.size):
                self.filenames.append("%s/%s_%05i.%s" % (path_to_folder,\
                    self.imageprefix,i+self.first_imageno,self.imagesuffix))
        else:
            for i in range(self.th.size):
                self.filenames.append("%s/%s/%s_%05.cbf" % (path_to_folder,\
                    self.filename_base,self.filename_base,self.first_imageno+i))
            
    def get_raw_P3_img(self,img):
        return P3_Image(self.filenames[img])
        
    # returns single default image!
    def get_raw_img(self,img):
        return P3_Image(self.filenames[img])
    
    def get_P3_img(self,img):
        imgdata = P3_Image(self.filenames[img])
        imgdata.counters['Time'] = self.time[img]
        imgdata.counters['TrigTime'] = self.exposure_time[img]
        imgdata.counters['th'] = self.th[img]
        imgdata.counters['om'] = self.omega[img]
        if self.current is not None:
            imgdata.counters['srcur'] = self.current[img]
        return imgdata
    
    def __getitem__(self,key):
        return self.get_P3_img(key)
    
    def __len__(self):
        return self.nopoints
        
class FioFastsweep_original(object):
    def __init__(self,fiofilepath,scanno=None):
        self.path,fname = os.path.split(fiofilepath)
        fnamewotsuff, suffix = fname.split('.')
        prefix, nostr = fnamewotsuff.split('_')
        if scanno is None:
            scanno = int(nostr)
        fiofilecorrectedscanno = "%s/%s_%05i.%s" % (self.path,prefix,scanno,suffix) 
        #self.dtime = []
        #self.current = [] # sr current or some other monitor is missing!
        datafile = FioFile(fiofilecorrectedscanno)
        self.th = datafile.data['idrz1']
        meandt = np.mean(np.diff(datafile.data['time']))
        self.exposure_time = np.diff(datafile.data['time'],append=(datafile.data['time'][-1] + meandt))
        self.time = datafile.data['time']
        #self.th = -1*self.th  # orientation of th is inversed
        self.omega = -1*self.th
        #self.exposure_time = np.full_like(self.th,self.exposure_time)
        self.directFolder = True
        self.imagenames = datafile.data['filename']
        print(self.imagenames)
        self.imageprefix,fmiddle = self.imagenames[0].split('_')
        self.imagesuffix = fmiddle.split('.')[-1]
        self.first_imageno = 1
        self.current = None
        #self.set_image_folder(self.path+ "/" + detector )
        self.nopoints = self.th.size
            

        
    def set_th_offset(self,offset):
        self.th += offset
        self.omega = -1*self.th

        
    def set_image_folder(self,path_to_folder):
        #self.filenames = [None]*len(self.th)
        self.filenames = []
        if self.directFolder:
            for i in range(self.th.size):
                self.filenames.append("%s/%s_%05i.%s" % (path_to_folder,\
                    self.imageprefix,i+self.first_imageno,self.imagesuffix))
        else:
            for i in range(self.th.size):
                self.filenames.append("%s/%s/%s_%05.cbf" % (path_to_folder,\
                    self.filename_base,self.filename_base,self.first_imageno+i))
            
    def get_raw_P3_img(self,img):
        return P3_Image(self.filenames[img])
        
    # returns single default image!
    def get_raw_img(self,img):
        return P3_Image(self.filenames[img])
    
    def get_P3_img(self,img):
        imgdata = P3_Image(self.filenames[img])
        imgdata.counters['Time'] = self.time[img]
        imgdata.counters['TrigTime'] = self.exposure_time[img]
        imgdata.counters['th'] = self.th[img]
        imgdata.counters['om'] = self.omega[img]
        if self.current is not None:
            imgdata.counters['srcur'] = self.current[img]
        return imgdata
    
    def __getitem__(self,key):
        return self.get_P3_img(key)
    
    def __len__(self):
        return self.nopoints
            
    
        
        
class CrudeThScan(object):
    def __init__(self,logfilepath,detector='PE1',darkfile=None):
        self.path,_ = os.path.split(logfilepath)
        self.th = []
        self.dtime = []
        self.current = []
        self.dark = None
        with open(logfilepath,'r') as f:
            for l in f:
                if l.startswith('#'):
                    if 'Exposure time' in l:
                        lsp = l.split()
                        self.exposure_time = float(lsp[2][:-1])
                    continue
                lspl = l.split()
                self.th.append(float(lspl[0]))
                self.current.append(float(lspl[2]))
                self.dtime.append(datetime.strptime(lspl[3] + ' ' + lspl[4],"%Y-%m-%d %H:%M:%S.%f"))
            self.th = np.array(self.th)
            self.current = np.array(self.current)
            self.time = []
            for t in self.dtime:
                self.time.append(t.day*24*3600 + t.hour*3600 + t.minute*60 + t.second + t.microsecond*1e-6)
            self.time = np.array(self.time)
        #self.th = -1*self.th  # orientation of th is inversed
        self.omega = -1*self.th
        self.exposure_time = np.full_like(self.th,self.exposure_time)
        self.directFolder = True
        self.filename_base = 'pe'
        self.first_imageno = 1
        self.set_image_folder(self.path+ "/" + detector )
        self.nopoints = self.th.size
        
        if darkfile is not None:
            with fabio.open(darkfile) as img:
                self.dark = img.data.astype(np.float64)
            

        
    def set_th_offset(self,offset):
        self.th += offset
        self.omega = -1*self.th

        
    def set_image_folder(self,path_to_folder):
        #self.filenames = [None]*len(self.th)
        self.filenames = []
        if self.directFolder:
            for i in range(self.th.size):
                self.filenames.append("%s/%s%05i.tif.gz" % (path_to_folder,\
                    self.filename_base,self.first_imageno+i))
        else:
            for i in range(self.th.size):
                self.filenames.append("%s/%s/%s_%05.cbf" % (path_to_folder,\
                    self.filename_base,self.filename_base,self.first_imageno+i))
            
    def get_raw_PE_img(self,img):
        return PE_Image(self.filenames[img],self.dark)
        
    # returns single default image!
    def get_raw_img(self,img):
        return PE_Image(self.filenames[img],self.dark)
    
    def get_PE_img(self,img):
        imgdata = PE_Image(self.filenames[img],self.dark)
        imgdata.counters['Time'] = self.time[img]
        imgdata.counters['TrigTime'] = self.exposure_time[img]
        imgdata.counters['th'] = self.th[img]
        imgdata.counters['om'] = self.omega[img]
        if self.current is not None:
            imgdata.counters['srcur'] = self.current[img]
        return imgdata
    
    def __getitem__(self,key):
        return self.get_PE_img(key)
    
    def __len__(self):
        return self.nopoints
        
    
