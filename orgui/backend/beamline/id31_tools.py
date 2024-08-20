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
import scipy.interpolate
import datetime
import re
import warnings
import traceback
import fabio
import h5py
import os
import copy
from silx.io import dictdump
import silx.io

from .ID31DiffractLinTilt import ID31DiffractLinTilt


class P3_Image:
    def __init__(self, filename,edf=True):
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
            image = fabio.open(filename)
            #cbf = PilatusCBF.PilatusCBF(filename)
            self.header = image.header
            #self.header = cbf.getInfo('') # only detector related information
            self.motors = dict() # not available in cbf
            self.counters = dict() # not available in cbf
            self.img = image.data.astype(np.float64)
            #del cbf
class h5_Image:
    
    def __init__(self, data):
        
        self.img = data 
        self.motors = dict()
        self.counters = dict()

class MPX_Image:
    def __init__(self, filename):
        edf = EdfFile.EdfFile(filename, 'r')
        self.header = edf.GetHeader(0)
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
        self.img = np.array(edf.GetData(0))


class PE_Image:
    def __init__(self, filename):
        edf = EdfFile.EdfFile(filename, 'r')
        self.header = edf.GetHeader(0)
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

        self.img = np.array(edf.GetData(0))
        # detector was mounted rotated by 90deg to the right
        self.img = np.rot90(self.img, 3)
        
    
    

# currently only set up for th scans in TOMO session, cbf fileformat
class Fastscan(object):
    def __init__(self,fastscan_specfile, scanno):
        id31_fastscan_spec = specfile.Specfile(fastscan_specfile)
        scan = id31_fastscan_spec[scanno-1]
        filename_line = scan.header('C next_image_file')[0]
        self.filename_base, next_image_name = filename_line.split('/')[-2:]
        
        imagnostr = next_image_name.split('.')[0].split('_')[-1]
        imagnostr_clean = ''.join(list(filter( str.isdigit ,imagnostr)))
        self.directFolder = False
        if imagnostr_clean != imagnostr:
            warnings.warn(("Imageno formatting doesn't follow the fastscan numbering scheme.\n"
                           "Fall back to direct folder selection. (please select folder manually!!)")
                          ,UserWarning)
            self.directFolder = True
            cutposition = next_image_name.find(imagnostr_clean)
            self.filename_base = next_image_name[:cutposition]
        self.first_imageno = int(imagnostr_clean)
        
        #self.first_imageno = int( next_image_name.split('.')[0].split('_')[-1] )
        try:
            self.starttime = datetime.datetime.strptime(scan.header('D')[0][3:]\
                                                    , "%a %b %d %H:%M:%S %Y") 
        except Exception as e:
            warnings.warn("can not read starttime, don't know why yet!. Error: %s" % e)
        self.sr_status = dict(zip(scan.header('UMI0')[0].split()[1:],\
                        re.split(r'\s{2,}', scan.header('UMI1')[0])[1:]))
        scan = id31_fastscan_spec[scanno]
        self.th = np.mean([scan.datacol('th_UpPos'),\
              scan.datacol('th_DownPos')],axis=0)
        self.omega = -1*self.th
        self.exposure_time = scan.datacol('TrigTime')
        commands = scan.command().split()
        self.nopoints = int(commands[-4])
        totaltime = (self.nopoints * float(commands[2]))/float(commands[-3])
        self.time = np.linspace(0,totaltime,self.nopoints)
        self.current = None
        
    def set_th_offset(self,offset):
        self.th += offset
        self.omega = -1*self.th

        
    def set_image_folder(self,path_to_folder):
        #self.filenames = [None]*len(self.th)
        self.filenames = []
        if self.directFolder:
            for i in range(self.th.size):
                self.filenames.append("%s/%s%04i.cbf" % (path_to_folder,\
                    self.filename_base,self.first_imageno+i))
        else:
            for i in range(self.th.size):
                self.filenames.append("%s/%s/%s_%04i.cbf" % (path_to_folder,\
                    self.filename_base,self.filename_base,self.first_imageno+i))
            
    def get_raw_p3_img(self,img):
        return P3_Image(self.filenames[img],False)
    
    # returns single default image!
    def get_raw_img(self,img):
        return P3_Image(self.filenames[img],False)
        
    def get_p3_img(self,img):
        imgdata = P3_Image(self.filenames[img],False)
        imgdata.counters['Time'] = self.time[img]
        imgdata.counters['TrigTime'] = self.exposure_time[img]
        imgdata.counters['th'] = self.th[img]
        imgdata.counters['om'] = self.omega[img]
        if self.current is not None:
            imgdata.counters['srcur'] = self.current[img]
        return imgdata

    # for normalization
    # example: 0.004213 mA/s was ok for 16 bunch, top up mode (Feb 2018) 
    def estimateRingCurrent(self,lossPerSecond):
        self.current = float(self.sr_status['Current']) -\
        lossPerSecond*self.time
        return self.current

    def get_p3_img_recal(self,img):
        imgdata = self.get_p3_img(img)
        imgdata.img = np.multiply(imgdata.img, id31_cal_matrix)
        return imgdata
    
    def __getitem__(self,key):
        return self.get_p3_img(key)
    
    def __len__(self):
        return self.nopoints

class BlissScan_EBS(Fastscan):
    def __init__(self, hdffilepath_orNode=None, scanno=None, loadimg=True, saveh5data=False, **keyargs):
        #self.scanname = scanname
        self.cameras = []
        data_1 = None
        data_2 = None
        self.loadimg=loadimg
        excludenames = None if self.loadimg else ['p3']
        
        if hdffilepath_orNode is None:
            return
        self.hdffilepath_orNode=hdffilepath_orNode    

        if isinstance(hdffilepath_orNode, str):
            filepath,filename = os.path.split(hdffilepath_orNode)
            _,filename_noext = os.path.split(filepath)
            #filename_noext = filename.split('.')[0]
            self.filename_base = filename_noext #filename_noext[:filename_noext.rfind('_')]
            with silx.io.open(hdffilepath_orNode) as f:
                #print([d for d in f])
                for d in f:
                    scansuffix = d.split('_')[-1]
                    scanname_nosuffix = '_'.join(d.split('_')[:-1])
                    scanno_s, subscanno = scansuffix.split('.')
                    if int(scanno_s) == scanno:
                        break
                else:
                    raise IOError("Scan number %s not found in file" % scanno)
                self.scanno1 = str(scanno) + '.' + '1'
                self.scanno2 = str(scanno) + '.' + '2'
                self.scanname_1 = scanname_nosuffix + '_' + self.scanno1
                self.scanname_2 = scanname_nosuffix + '_' + self.scanno2
                self.scanname = self.scanname_1
                if self.scanname_1 in f:
                    data_1 = dictdump.h5todict(f,self.scanname_1, exclude_names=excludenames)
                    if 'p3' in f[self.scanname_1]['measurement']:
                        self.cameras.append('p3')
                        self.nopoints = f[self.scanname_1]['measurement']['p3'].shape[0]
                    if self.scanname_2 in f:
                        data_2 = dictdump.h5todict(f,self.scanname_2, exclude_names=excludenames)
                else:
                    data_1 = dictdump.h5todict(f,self.scanno1, exclude_names=excludenames)                
                    if 'p3' in f[self.scanno1]['measurement']:
                        self.cameras.append('p3')
                        self.nopoints = f[self.scanno1]['measurement']['p3'].shape[0]
                    if self.scanno2 in f:
                        data_2 = dictdump.h5todict(f,self.scanno2, exclude_names=excludenames)
                    
        else:
            hdffilepath = hdffilepath_orNode.local_filename
            filepath,filename = os.path.split(hdffilepath)
            _,filename_noext = os.path.split(filepath)
            self.filename_base = filename_noext
            
            f = hdffilepath_orNode.file
            for d in f:
                scansuffix = d.split('_')[-1]
                scanname_nosuffix = '_'.join(d.split('_')[:-1])
                scanno_s, subscanno = scansuffix.split('.')
                
                if int(scanno_s) == scanno:
                    break
            else:
                raise IOError("Scan number %s not found in file" % scanno)
            self.scanno1 = str(scanno) + '.' + '1'
            self.scanno2 = str(scanno) + '.' + '2'
            self.scanname_1 = scanname_nosuffix + '_' + self.scanno1
            self.scanname_2 = scanname_nosuffix + '_' + self.scanno2
            self.scanname = self.scanname_1
            
            if self.scanname_1 in f:
                data_1 = dictdump.h5todict(f,self.scanname_1, exclude_names=excludenames)
                if 'p3' in f[self.scanname_1]['measurement']:
                    self.cameras.append('p3')
                    self.nopoints = f[self.scanname_1]['measurement']['p3'].shape[0]
                if 'mpx' in f[self.scanname_1]['measurement']:
                    self.cameras.append('mpx')
                    self.nopoints = f[self.scanname_1]['measurement']['mpx'].shape[0]
                if self.scanname_2 in f:
                    data_2 = dictdump.h5todict(f,self.scanname_2, exclude_names=excludenames)
            else:
                data_1 = dictdump.h5todict(f,self.scanno1, exclude_names=excludenames)
                if 'p3' in f[self.scanno1]['measurement']:
                    self.cameras.append('p3')
                    self.nopoints = f[self.scanno1]['measurement']['p3'].shape[0]
                if 'mpx' in f[self.scanno1]['measurement']:
                    self.cameras.append('mpx')
                    self.nopoints = f[self.scanno1]['measurement']['mpx'].shape[0]
                if self.scanno2 in f:
                    data_2 = dictdump.h5todict(f,self.scanno2, exclude_names=excludenames)
        if saveh5data:
            self.data_1 = data_1
            self.data_2 = data_2
        
        if 'title' in data_1:
            self.title = data_1['title']
        
        self.positioners = data_1['instrument']['positioners']
        if 'p3' in data_1['measurement']:
            self.cameras.append('p3')
            self.nopoints = data_1['measurement']['p3'].shape[0]
            if self.loadimg:
                print('load p3 images')
                self.p3 = data_1['measurement']['p3'][()]
            
        if 'mpx' in data_1['measurement']:
            self.mpx = data_1['measurement']['mpx'][()]
            self.nopoints = data_1['measurement']['mpx'].shape[0]
            self.cameras.append('mpx')        
        
        if 'th' in data_1['measurement']:
            self.th = data_1['measurement']['th'][:self.nopoints]
            if 'th_trig' in data_1['measurement'] and 'th_delta' in data_1['measurement']:
                self.th = data_1['measurement']['th_trig'][:self.nopoints] + data_1['measurement']['th_delta'][:self.nopoints] / 2
            self.axisname = 'th'
            self.mu = self.positioners['mu']

        elif 'uth' in data_1['measurement']:
            self.th = data_1['measurement']['uth'][:self.nopoints]
            self.axisname = 'th'
            self.mu = self.positioners['mu'] *-1
            
        elif 'mu' in data_1['measurement']:
            self.mu = data_1['measurement']['mu'][:self.nopoints]
            self.axisname = 'mu'
            self.th = self.positioners['th']
        
        elif 'linai' in data_1['measurement']:
            lintomu = ID31DiffractLinTilt()
            if 'muoffset' in keyargs:
                lintomu.config['muoffset'] = keyargs['muoffset']
            self.linai = data_1['measurement']['linai'][:self.nopoints]
            if 'linai_trig' in data_1['measurement'] and 'linai_delta' in data_1['measurement']:
                self.linai = data_1['measurement']['linai_trig'][:self.nopoints] + data_1['measurement']['linai_delta'][:self.nopoints] / 2
            self.mu = lintomu.linai_to_mu(self.linai)
            self.axisname = 'mu'
            self.th = self.positioners['th']
        else:
            self.axisname = 'time'
            self.th = self.positioners['th']
            self.mu = self.positioners['mu']
            
        if 'srcur' in data_1['measurement']:
            self.srcur = data_1['measurement']['srcur'][:self.nopoints]
        else:
            self.srcur = np.ones(self.nopoints)
            
        if 'timer_trig' in data_1['measurement']:
            self.time = data_1['measurement']['timer_trig'][:self.nopoints]
        elif 'elapsed_time' in data_1['measurement']:
            self.time = data_1['measurement']['elapsed_time'][:self.nopoints]
        else:
            raise IOError("Cannot find time counter in scan")
        
        if 'epoch_trig' in data_1['measurement']:
            self.epoch = data_1['measurement']['epoch_trig'][:self.nopoints]
        elif 'epoch' in data_1['measurement']:
            self.epoch = data_1['measurement']['epoch'][:self.nopoints]
        else:
            raise IOError("Cannot find epoch counter in scan")
        
        if data_2 is not None:
            if 'mondio' in data_2['measurement']:
                try:
                    mondiointer = scipy.interpolate.interp1d(data_2['measurement']['epoch'],data_2['measurement']['mondio'])
                    self.mondio = mondiointer(self.epoch)
                except ValueError as e0:
                    warnings.warn('Cannot interpolate mondio in scan %s out of dataset %s \n%s' % (self.scanno2, filename, traceback.format_exc()))
                
            if 'potential' in data_2['measurement']:
                try:
                    potentialinter = scipy.interpolate.interp1d(data_2['measurement']['epoch'],data_2['measurement']['potential'])
                    self.potential = potentialinter(self.epoch)
                except ValueError as e1:
                    warnings.warn('Cannot interpolate potential in scan %s out of dataset %s \n%s' % (self.scanno2, filename, traceback.format_exc()))
                
            if 'current' in data_2['measurement']:
                try:
                    currentinter = scipy.interpolate.interp1d(data_2['measurement']['epoch'],data_2['measurement']['current'])
                    self.current = currentinter(self.epoch)
                except ValueError as e2:
                    warnings.warn('Cannot interpolate current in scan %s out of dataset %s \n%s' % (self.scanno2, filename, traceback.format_exc()))
        else: 
            if 'potential' in data_1['measurement']:
                self.potential = data_1['measurement']['potential'][:self.nopoints]
            #else:
            #    raise IOError("Didn't find potential counter.")
            if 'current' in data_1['measurement']:
                self.current = data_1['measurement']['current'][:self.nopoints]
            #else:
            #    raise IOError("Didn't find current counter.")
        if 'mondio' in data_1['measurement']:
            self.mondio = data_1['measurement']['mondio'][:self.nopoints]
        
        if not hasattr(self,'mondio'):
            self.mondio = np.ones(self.nopoints)
        self.mondio = self.mondio/np.mean(self.mondio)
        self.srcur = self.srcur/np.mean(self.srcur)
        if 'timer_delta' in data_1['measurement']:
            self.exposure_time = data_1['measurement']['timer_delta'][:self.nopoints]
        elif 'sec' in data_1['measurement']:
            self.exposure_time = data_1['measurement']['sec'][:self.nopoints]
        else:
            raise IOError("Cannot find exposure time in scan")

        self.relative_exposure = self.exposure_time / np.mean(self.exposure_time)
        
        self.offsetindex = 0
        self.scandatapoints = self.nopoints
        self.imageno = np.arange(self.nopoints)
        
        if hasattr(self,"th"):
            self.omega = -1*self.th
            
        self.axis = getattr(self,self.axisname)
            
    # returns single default image!
    def get_raw_img(self,img):
        if not hasattr(self, 'p3'):
            if isinstance(self.hdffilepath_orNode, str):
                with silx.io.h5py_utils.File(self.hdffilepath_orNode, 'r') as f:
                    if self.scanname_1 in f:
                        data_1 = f[self.scanname_1]
                    else:
                        data_1 = f[self.scanno1]
                    image = data_1['measurement']['p3'][img][()]
            else:
                f = self.hdffilepath_orNode.file
                if self.scanname_1 in f:
                    data_1 = f[self.scanname_1]
                else:
                    data_1 = f[self.scanno1]
                image = data_1['measurement']['p3'][img][()]
        else:
            image = self.p3[img]
        return h5_Image(image)
            
    def slice(self,startno,endno):
        if startno > endno:
            startno, endno = endno, startno
 
        fscan = BlissScan_EBS()
        
        for dat in self.__dict__:
            if isinstance(self.__dict__[dat], np.ndarray):
                if self.__dict__[dat].ndim > 0:
                    fscan.__dict__[dat] = copy.deepcopy(self.__dict__[dat][startno:endno])
                else:
                    fscan.__dict__[dat] = copy.deepcopy(self.__dict__[dat])
            else:
                fscan.__dict__[dat] = copy.deepcopy(self.__dict__[dat])
                
        if self.loadimg==False:    
            with silx.io.h5py_utils.File(self.hdffilepath_orNode, 'r') as f:
                data_1 = f[self.scanname_1]
                fscan.p3 = data_1['measurement']['p3'][startno:endno][()]
        
        fscan.offsetindex = copy.deepcopy(startno)
        fscan.nopoints = copy.deepcopy(endno - startno)
        
        return fscan
        

    def get_p3_img(self,img):
        #img = img - self.offsetindex
        imgdata = self.get_raw_img(img)
        imgdata.counters['Time'] = self.time[img]
        imgdata.counters['exposure'] = self.exposure_time[img]
        imgdata.counters['TrigTime'] = self.relative_exposure[img]
        imgdata.counters['imageno'] = self.imageno[img]
        if hasattr(self,"th"):
            if self.axisname == 'th':
                imgdata.counters['th'] = self.th[img]
                imgdata.counters['om'] = self.omega[img]
            else:
                imgdata.counters['th'] = self.th
                imgdata.counters['om'] = self.omega
        #if hasattr(self,"fixed_th"):        
        #    imgdata.counters['th'] = self.fixed_th
        #    imgdata.counters['om'] = self.fixed_th*-1
        if hasattr(self,"mu"):
            if self.axisname == 'mu':
                imgdata.counters['mu'] = self.mu[img]
            else:
                imgdata.counters['mu'] = self.mu
        if self.srcur is not None:
            imgdata.counters['srcur'] = self.srcur[img]
        
        if self.mondio is not None:
            imgdata.counters['mondio'] = self.mondio[img]
            
        return imgdata                   
    
# TODO: change subclassing! this is terrible, but works for the moment
class BlissScan(Fastscan):
    def __init__(self,hdffilepath, scanname):
        hdffilepath = os.path.abspath(hdffilepath)
        filepath,filename = os.path.split(hdffilepath)
        self.filepath = filepath
        #print(filepath,filename)
        self.scanname = scanname
        _,filename_noext = os.path.split(filepath)
        #filename_noext = filename.split('.')[0]
        self.filename_base = filename_noext #filename_noext[:filename_noext.rfind('_')]
        #print(self.filename_base)
        self.cameras = []
        #self.fixed_th = None
        #self.th = None
        #self.mu = None
        with silx.io.h5py_utils.File(hdffilepath,'r') as f:
            if scanname not in f:
                raise Exception("Scan name %s does not exist in file" % scanname)
            scangroup = f[scanname]
            command = scangroup['title'][()].decode('utf-8')
            self.title = command
            measgroup = scangroup['measurement']
            if command.startswith("cscan"):
                for key in measgroup['th']:
                    try:
                        data = measgroup['th'][key][:]
                        dataname = key.split(':')[1]
                    except ValueError:
                        continue
                    setattr(self,dataname,data)
                    
                self.th = self.th_mean
                self.current = None
                self.nopoints = self.th.size
                self.time = self.time_mean
                self.exposure_time = self.time_delta
                
                if "p3" in measgroup:
                    self.cameras.append('p3')
                    
                if 'mpx' in measgroup:
                    self.cameras.append('mpx')
                self.valid = np.zeros(self.th.size,dtype=bool)
                
                self.axisname = 'th'
                self.axis = self.th
            else:
                scanmotorgroup = None
                for key in measgroup:
                    if key.startswith('group'):
                        scanmotorgroup = measgroup[key]
                        break
                if scanmotorgroup is not None:
                    for key in scanmotorgroup:
                        try:
                            data = scanmotorgroup[key][:]
                            dataname = key.split(':')[1]
                            self.axisname = dataname
                            self.axis = data
                        except ValueError:
                            continue
                    
                    setattr(self,dataname,data)
                    
                for key in measgroup['timer']:
                    try:
                        data = measgroup['timer'][key][:]
                        dataname = key.split(':')[1]
                    except ValueError:
                        continue
                    setattr(self,dataname,data)
                    
                if not hasattr(self, "axisname"):
                    self.axisname = 'time'
                    self.axis = self.elapsed_time
                    
                positioner_group = measgroup['instrument/positioners']
                if 'th' in positioner_group:
                    self.th = np.full_like(self.srcur,positioner_group['th'][()])
                
                #self.current = self.srcur
                self.nopoints = self.srcur.size
                self.time = self.elapsed_time
                self.exposure_time = np.full_like(self.srcur.size, float(command.split(' ')[-1]) )
        
                if "p3" in measgroup:
                    self.cameras.append('p3')
                    
                if 'mpx' in measgroup:
                    self.cameras.append('mpx')
                self.valid = np.zeros(self.srcur.size,dtype=bool)
                
                #self.cameras = measgroup['musst'].keys()
            if hasattr(self,"th"):
                self.omega = -1*self.th
            
            
                #  To implement:
                #self.starttime = datetime.datetime.strptime(scan.header('D')[0][3:]\
                #                                        , "%a %b %d %H:%M:%S %Y")
            for item in self.__dict__:
                if isinstance(self.__dict__[item],np.ndarray):
                    self.__dict__[item] = np.ma.array(self.__dict__[item])
            
            for cam in self.cameras:
                imgfolder = os.path.join(self.filepath,cam)
                if os.path.isdir(imgfolder):
                    self.set_image_folder(imgfolder)
                    

            

    def set_image_folder(self,path_to_folder):
        #self.filenames = [None]*len(self.th)
        self.valid = np.zeros(len(self),dtype=bool)
        self.filenames = []
        camera = self.cameras[0] # p3
        for i in range(len(self)):
            name = "%s/%s_%s_%s_%04i.edf.gz" % (path_to_folder,\
                self.filename_base,camera,self.scanname,i)
            if os.path.isfile(name):
                self.filenames.append(name)
                self.valid[i] = True
            else:
                self.valid[i] = False
        if (~self.valid).any():
            for item in self.__dict__:
                
                if isinstance(self.__dict__[item],np.ndarray) and self.__dict__[item].shape == self.valid.shape:
                    self.__dict__[item] = self.__dict__[item][self.valid]
        self.nopoints = int(self.valid.sum())

    def get_p3_img(self,img):
        imgdata = P3_Image(self.filenames[img],False)
        imgdata.counters['Time'] = self.time[img]
        imgdata.counters['TrigTime'] = self.exposure_time[img]
        if hasattr(self,"th"): 
            imgdata.counters['th'] = self.th[img]
            imgdata.counters['om'] = self.omega[img]
        #if hasattr(self,"fixed_th"):        
        #    imgdata.counters['th'] = self.fixed_th
        #    imgdata.counters['om'] = self.fixed_th*-1
        if hasattr(self,"mu"):       
            imgdata.counters['mu'] = self.mu[img]
        if self.current is not None:
            imgdata.counters['srcur'] = self.srcur[img]
        return imgdata                    

