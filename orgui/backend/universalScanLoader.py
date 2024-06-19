import re
import os
import fabio
import numpy as np

from .scans import h5_Image

class ImportImagesScan():

    def __init__(self, imgpath):
        self.filename = imgpath

        self.inpath = self.find_files()

        img_data = fabio.open(self.inpath[0] + self.inpath[1][0])
        self.shape = (img_data.data.shape[0], img_data.data.shape[1])

        self.FramesPerFile = img_data.nframes
        if self.FramesPerFile > 1:
            last_file = fabio.open(self.inpath[0] + self.inpath[1][-1])
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
        suffix = re_str.findall(self.filename)[0]
        imagePrefix = self.filename.removesuffix(suffix)

        #found_scannrs = [e[1:-4] for e in found_scanfiles] #only the scan numbers, eg. 00015

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
            img_data = fabio.open(self.inpath[0] + self.inpath[1][index]).get_frame(frame)
        else:
            img_data = fabio.open(self.inpath[0] + self.inpath[1][i])

        return h5_Image(img_data.data)
        
    def set_raw_img(self, i, data): #for intensity simulation in the future.
        self.images[i] = data
        

class ImportImagesScanOld():

    def __init__(self, detshape, axismin, axismax, points, axis='th', fixed=0.,imgpath=[0,0],framesperfile=1):
        
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
        self.title = "manually loaded scan %s %s %s %s" % (self.axisname,axismin,axismax,points)
        self.impath = imgpath

        self.framesperfile = framesperfile
        self.images = np.zeros((1,*detshape))
        
    def __len__(self):
        return self.nopoints
        
    def get_raw_img(self, i):
        if self.framesperfile > 1:
            index = i // self.framesperfile
            frame = i % self.framesperfile
            img_data = fabio.open(self.impath[0] + self.impath[1][index]).get_frame(frame)
        else:
            img_data = fabio.open(self.impath[0] + self.impath[1][i])

        return h5_Image(img_data.data)
        
    def set_raw_img(self, i, data): #for intensity simulation in the future.
        self.images[i] = data

    def find_files(self):
        re_str = re.compile(r'_\d+' + os.path.splitext(self.filename)[1]) #define search string. It matches a file source with syntax name_0000i.extension -> may need to be adapted 
        selected_directory = os.path.dirname(os.path.abspath(self.filename))
        filenames = ''.join(os.listdir(selected_directory))

        found_scanfiles = re_str.findall(filenames) #list of found filenames (suffix only)
        suffix = re_str.findall(self.filename)[0]
        imagePrefix = self.filename.removesuffix(suffix)

        #found_scannrs = [e[1:-4] for e in found_scanfiles] #only the scan numbers, eg. 00015
        img_data = fabio.open(imagePrefix + found_scanfiles[0]) # load first found scan image to get nr of pixels

        return imagePrefix, found_scanfiles
    
    def set_metadata(self, detshape, frames):
        self.detshape = detshape
        self.frames = frames
