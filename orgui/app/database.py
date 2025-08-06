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

import sys
import numpy as np
import silx.gui.hdf5
from silx.gui import icons
from silx.gui.data import DataViewerFrame
from silx.gui.hdf5.Hdf5TreeModel import Hdf5TreeModel
from silx.gui import qt
import tempfile
import silx.io.h5py_utils
from silx.io import nxdata
import silx.io.utils
import h5py

import datetime
import os
import traceback
import time

from ..datautils.xrayutils import unitcells
from ..datautils.xrayutils import HKLVlieg, CTRcalc
from ..datautils.xrayutils import DetectorCalibration
from .. import resources


from silx.io.dictdump import dicttonx ,nxtodict

DEFAULT_FILTERS = {  # Filters available with h5py/libhdf5
    "Raw": None,
    "GZip": "gzip",
    "LZF": "lzf",
}

FILTERS = {
    **DEFAULT_FILTERS
}

try:
    import hdf5plugin

    LOSSLESS_FILTERS = {
        "BZip2": hdf5plugin.BZip2(),
        "LZ4": hdf5plugin.LZ4(),
        "ZStd": hdf5plugin.Zstd(),
    }
    FILTERS.update(**LOSSLESS_FILTERS)

    BITSHUFFLE_FILTERS = {
        "Bitshuffle-lz4": hdf5plugin.Bitshuffle(cname='lz4'),
        "Bitshuffle-zstd": hdf5plugin.Bitshuffle(cname='zstd'),
    }
    FILTERS.update(**BITSHUFFLE_FILTERS)
        
    BLOSC_FILTERS = {}
    for cname in ('lz4', 'blosclz', 'lz4', 'lz4hc', 'snappy', 'zlib', 'zstd'):
        for shuffle_name, shuffle in [('NoShuffle', hdf5plugin.Blosc.NOSHUFFLE),
                                      ('Shuffle', hdf5plugin.Blosc.SHUFFLE),
                                      ('BitShuffle', hdf5plugin.Blosc.BITSHUFFLE)]:
            for clevel in [5]: #(1, 3, 5, 9):
                BLOSC_FILTERS[f"Blosc-{cname}-{shuffle_name}-{clevel}"] = hdf5plugin.Blosc(
                    cname=cname, clevel=clevel, shuffle=shuffle)
    FILTERS.update(**BLOSC_FILTERS)

    BLOSC2_FILTERS = {}
    for cname in ('lz4', 'blosclz', 'lz4', 'lz4hc', 'zlib', 'zstd'):
        for filters_name, filters in [('NoFilter', hdf5plugin.Blosc2.NOFILTER),
                                    ('Shuffle', hdf5plugin.Blosc2.SHUFFLE),
                                    ('BitShuffle', hdf5plugin.Blosc2.BITSHUFFLE)]:
            for clevel in [5]: # (1, 3, 5, 9):
                BLOSC2_FILTERS[f"Blosc2-{cname}-{filters_name}-{clevel}"] = hdf5plugin.Blosc2(
                    cname=cname, clevel=clevel, filters=filters)
    FILTERS.update(**BLOSC2_FILTERS)

except Exception as e:
    print("Unable to load hdf5plugin compression filters:")
    traceback.print_exc()

class ConfigData(qt.QObject):
    def __init__(self, config=None):
        self.detector = DetectorCalibration.Detector2D_SXRD()
        
        pass
    
    #@classmethod
    def readConfig(self, filename):
        config = configparser.ConfigParser()
        config.read(configfile)
        
        machine = config['Machine']
        lattice = config['Lattice']
        diffrac = config['Diffractometer']
            
        self.azimuth = np.deg2rad(diffrac.getfloat('azimuthal_reference',0))
        self.polaxis = np.deg2rad(diffrac.getfloat('polarization_axis',90))
        self.polfactor = diffrac.getfloat('polarization_factor',0)
            
        sdd = machine.getfloat('SDD',0.729) #m
        E =  machine.getfloat('E',78.0) #keV
        pixelsize = machine.getfloat('pixelsize',172e-6) #m
        cpx = machine.getfloat('cpx',731)
        cpy = machine.getfloat('cpy',1587)
        cp = [cpx,cpy]

        self.mu = np.deg2rad(diffrac.getfloat('mu',0.05))
        self.chi = np.deg2rad(diffrac.getfloat('chi',0.0))
        self.phi = np.deg2rad(diffrac.getfloat('phi',0.0))
            

        a1 = lattice.getfloat('a1',-1)
        a2 = lattice.getfloat('a2',-1)
        a3 = lattice.getfloat('a3',-1)
        alpha1 = lattice.getfloat('alpha1',-1)
        alpha2 = lattice.getfloat('alpha2',-1)
        alpha3 = lattice.getfloat('alpha3',-1)
        self.n = 1 - lattice.getfloat('refractionindex',0.0)
            
        lat = np.array([a1,a2,a3])
            
        latticeoverride = True
        latangle = np.array([alpha1,alpha2,alpha3])
        if np.any(lat < 0.) or np.any(latangle < 0):
            latticeoverride = False
            a1 = a2 = a3 = 1.
            alpha1 = alpha2 = alpha3 = 90.
            
        
        self.crystal = CTRcalc.UnitCell([a1,a2,a3],[alpha1,alpha2,alpha3])
        self.crystal.addAtom('Pt',[0.,0.,0.],0.1,0.1,1.)
        self.crystal.setEnergy(E*1e3)
        
        self.ubCal = HKLVlieg.UBCalculator(self.crystal,E)
        self.ubCal.defaultU()
        self.angles = HKLVlieg.VliegAngles(self.ubCal)
        
        if 'crystal' in lattice:
            idx = self.crystalparams.crystalComboBox.findText(lattice['crystal'],qt.Qt.MatchFixedString)
            if idx == -1:
                qt.QMessageBox.warning(self,"Did not find crystal","Can not find crystal <%s> \nException occured during read of configfile %s,\nException:\n%s" % (unitcellsconfigfile,traceback.format_exc()))
            else:
                self.crystalparams.crystalComboBox.setCurrentIndex(idx)
                self.crystalparams.onSwitchCrystal(idx)
                #self.crystal = self.crystalparams.getCrystal()

        if latticeoverride:
            print('foo')
            self.crystal.setLattice([a1,a2,a3],[alpha1,alpha2,alpha3])
        

        if 'poni' in machine:
            if machine['poni']:
                self.detectorCal.load(machine['poni'])
                self.ubCal.setLambda(self.detectorCal.get_wavelength()*1e10)
                
            else:
                self.detectorCal.setFit2D(sdd*1e3,cpx,cpy,pixelX=pixelsize*1e6, pixelY=pixelsize*1e6)
                self.detectorCal.set_wavelength(self.ubCal.getLambda()*1e-10)
                self.detectorCal.detector.shape = (2880,2880) # Perkin 
                
        else:
            self.detectorCal.setFit2D(sdd*1e3,cpx,cpy,pixelX=pixelsize*1e6, pixelY=pixelsize*1e6)
            self.detectorCal.set_wavelength(self.ubCal.getLambda()*1e-10)
            self.detectorCal.detector.shape = (2880,2880)
            
        self.detectorCal.setAzimuthalReference(self.azimuth)
        self.detectorCal.setPolarization(self.polaxis,self.polfactor)
        
        fit2dCal = self.detectorCal.getFit2D()

        paramlist = [self.ubCal.getEnergy(),self.mu,fit2dCal['directDist']/1e3,
                     fit2dCal['pixelX']*1e-6,[fit2dCal['centerX'],fit2dCal['centerY']],self.polaxis
                     ,self.polfactor,self.azimuth,self.chi,self.phi]
        self.crystalparams.setValues(self.crystal,self.n)
        self.machineParams.setValues(paramlist)
        
        


class DataBase(qt.QMainWindow):
    
    compression = FILTERS['Raw']
    
    sigChangeRockingScan = qt.Signal(object)
    
    def __init__(self, plot, parent=None):
        qt.QMainWindow.__init__(self, parent)
        self.nxfile = None
        self.filedialogdir = os.getcwd()
        self.plot = plot
        
        self.view = silx.gui.hdf5.Hdf5TreeView()
        self.view.setSortingEnabled(True)
        
        #self.view.addContextMenuCallback(self.nexus_treeview_callback)
        #self.hdf5model = self.view.findHdf5TreeModel()
        
        customModel = Hdf5TreeModel
        customModel.NAME_COLUMN = 0
        customModel.DESCRIPTION_COLUMN = 1
        customModel.VALUE_COLUMN = 2
        customModel.SHAPE_COLUMN = 3
        customModel.TYPE_COLUMN = 4
        customModel.NODE_COLUMN = 5
        customModel.LINK_COLUMN = 6
        
        self.hdf5model = customModel(self.view)
        self.view.setModel(self.hdf5model)

        self.dataviewer = DataViewerFrame.DataViewerFrame()
        self.dataviewerDialog = qt.QDialog(self)
        dvlayout = qt.QVBoxLayout()
        dvlayout.addWidget(self.dataviewer)
        self.dataviewerDialog.setLayout(dvlayout)
        self.dataviewerDialog.setModal(False)
        
        self.view.setExpandsOnDoubleClick(False)
        self.hdf5model.setFileMoveEnabled(True)
        #self.__treeModelSorted = silx.gui.hdf5.NexusSortFilterProxyModel(self.view)
        #self.__treeModelSorted.setSourceModel(self.hdf5model)
        #self.__treeModelSorted.sort(0, qt.Qt.AscendingOrder)
        #self.__treeModelSorted.setSortCaseSensitivity(qt.Qt.CaseInsensitive)

        #self.hdfTreeView.setModel(self.__treeModelSorted)
        
        self.view.doubleClicked.connect(self._onNEXUSDoubleClicked)
        self.view.addContextMenuCallback(self.nexus_treeview_callback)
        
        self.temp_directory = tempfile.TemporaryDirectory(dir=os.path.join(os.getcwd()))
        tempfilepath = os.path.join(self.temp_directory.name,"orgui_database.h5")
        self.createNewDBFile(tempfilepath)
        #self.add_nxdict(gauss)
        #self.add_nxdict(gauss)
        self.setCentralWidget(self.view)
        
        toolbar = qt.QToolBar("Database toolbar",self)
        
        
        
        newDatabaseAct = toolbar.addAction(resources.getQicon("document-nx-new"),"Create new orgui database")
        newDatabaseAct.triggered.connect(self.onNewDatabase)
        
        loadDatabaseAct = toolbar.addAction(icons.getQIcon("document-open"),"Open orgui database")
        loadDatabaseAct.triggered.connect(self.onOpenDatabase)
        
        savenewact = toolbar.addAction(icons.getQIcon("layer-nx"),"Select orgui database location")
        savenewact.triggered.connect(self.onSaveNewDBFile)
                
        saveact = toolbar.addAction(icons.getQIcon("document-save"),"Save orgui database")
        saveact.triggered.connect(self.onSaveDBFile)
        
        self.addToolBar(toolbar)
        
    def _onNEXUSDoubleClicked(self,index): # ToDo add try except with popup message!
        nodes = list(self.view.selectedH5Nodes())
        if len(nodes) > 0:
            obj = nodes[0]
            if obj.ntype is h5py.Dataset:
                roi_node = self.get_roinode(obj.h5py_object)
                if roi_node is not None:
                    self.plot_signal_callback(roi_node, obj)
                    return
            self.plot_default(obj.h5py_object)
                
    def plot_default(self, h5py_object):
        try:
            nxdat = nxdata.get_default(h5py_object)
            #print(nxdat.axes, nxdat.signal, nxdat.title)
            #print(nxdat.signal_name, nxdat.errors, nxdat.axes_names)
            if nxdat is not None and len(nxdat.axes) == 1:
                self.plot.addCurve(nxdat.axes[0], nxdat.signal, legend=nxdat.title, xlabel=nxdat.axes_names[0], ylabel=nxdat.signal_name, yerror=nxdat.errors)
        except Exception as e:
            traceback.print_exc()
            print("Cannot plot data: %s" % e)

    def get_roinode(self, obj):
        if silx.io.utils.get_h5_class(obj) is None:
            return None
            
        while(obj.name != '/'):
            meta = obj.attrs.get('orgui_meta', False)
            if meta and meta == 'roi':
                return obj
            obj = obj.parent
            
    def get_scannode(self, obj):
        if silx.io.utils.get_h5_class(obj) is None:
            return None
            
        while(obj.name != '/'):
            meta = obj.attrs.get('orgui_meta', False)
            if meta and meta == 'scan':
                return obj
            obj = obj.parent
            
    def nexus_treeview_callback(self,event):
        objects = list(event.source().selectedH5Nodes())
        obj = objects[0]  # for single selection
        menu = event.menu()
        action = qt.QAction("Refresh", menu)
        action.triggered.connect(lambda:  self.hdf5model.synchronizeH5pyObject(obj))
        menu.addAction(action)
        action = qt.QAction("details", menu)
        action.triggered.connect(lambda:  self.view_data_callback(obj))
        menu.addAction(action)
        
        if obj.ntype is h5py.Dataset:
            roi_node = self.get_roinode(obj.h5py_object)
            if roi_node is not None:
                action = qt.QAction("plot %s" % obj.name, menu)
                action.triggered.connect(lambda:  self.plot_signal_callback(roi_node, obj))
                menu.addAction(action)
        
        if obj.ntype is h5py.Group:
            meta = obj.h5py_object.attrs.get('orgui_meta', False)
            if meta and 'roi' in meta:
                if 'rocking' in meta:
                    action = qt.QAction("Show in rocking integration", menu)
                    action.triggered.connect(lambda:  self.onShowRoIntegrate(obj))
                    menu.addAction(action)
                menu.addSeparator()
                action = qt.QAction("delete", menu)
                action.triggered.connect(lambda:  self.delete_node(obj.h5py_object))
                menu.addAction(action)
                
            elif meta and 'rocking' in meta:
                action = qt.QAction("Show in rocking integration", menu)
                action.triggered.connect(lambda: self.onShowRoIntegrate(obj))
                menu.addAction(action)
                action = qt.QAction("rename", menu)
                action.triggered.connect(lambda:  self.onRenameNode(obj.h5py_object))
                menu.addAction(action)
                menu.addSeparator()
                action = qt.QAction("delete", menu)
                action.triggered.connect(lambda:  self.delete_node(obj.h5py_object))
                menu.addAction(action)
                
            if meta and 'scan' in meta:
                menu.addSeparator()
                action = qt.QAction("delete", menu)
                action.triggered.connect(lambda:  self.onDeleteScan(obj.h5py_object))
                menu.addAction(action)
                

            
        """    
        if obj.ntype is h5py.File:
            action = qt.QAction("remove", menu)
            action.triggered.connect(lambda:  self._onCloseFile())
            menu.addAction(action)
        """
    def plot_signal_callback(self, roi_node, dataset):
        try:
            nxdat = nxdata.get_default(roi_node)
            data = dataset.h5py_object[()]
            #print(nxdat.axes, nxdat.signal, nxdat.title)
            #print(nxdat.signal_name, nxdat.errors, nxdat.axes_names)
            if nxdat is not None and len(nxdat.axes) == 1:
                self.plot.addCurve(nxdat.axes[0], data, legend=nxdat.title + "_" + dataset.name, xlabel=nxdat.axes_names[0], ylabel=dataset.name)
        except Exception as e:
            traceback.print_exc()
            print("Cannot plot data: %s" % e)
        
    def view_data_callback(self,obj):
        self.dataviewer.setData(obj)
        self.dataviewerDialog.open()
        
    def onSaveNewDBFile(self):
        fileTypeDict = {'NEXUS Files (*.h5)': '.h5', 'All files (*)': '' }
        fileTypeFilter = ""
        for f in fileTypeDict:
            fileTypeFilter += f + ";;"
            
        filename, filetype = qt.QFileDialog.getSaveFileName(self,"Select database location",
                                                  self.filedialogdir,
                                                  fileTypeFilter[:-2])
        if filename == '':
            return
        
        self.filedialogdir = os.path.splitext(filename)[0]
        #filename += fileTypeDict[filetype]
        self.saveNewDBFile(filename)
        
    def onNewDatabase(self):
        fileTypeDict = {'NEXUS Files (*.h5)': '.h5', 'All files (*)': '' }
        fileTypeFilter = ""
        for f in fileTypeDict:
            fileTypeFilter += f + ";;"
            
        filename, filetype = qt.QFileDialog.getSaveFileName(self,"Select database location",
                                                  self.filedialogdir,
                                                  fileTypeFilter[:-2])
        if filename == '':
            return
        self.filedialogdir = os.path.splitext(filename)[0]
        self.createNewDBFile(filename)

    def onSaveDBFile(self):
        fileTypeDict = {'NEXUS Files (*.h5)': '.h5', 'All files (*)': '' }
        fileTypeFilter = ""
        for f in fileTypeDict:
            fileTypeFilter += f + ";;"
            
        filename, filetype = qt.QFileDialog.getSaveFileName(self,"Save database",
                                                  self.filedialogdir,
                                                  fileTypeFilter[:-2])
        if filename == '':
            return
        
        self.filedialogdir = os.path.splitext(filename)[0]
        #filename += fileTypeDict[filetype]
        self.saveDBFile(filename)
        
    def onOpenDatabase(self):
        fileTypeDict = {'NEXUS Files (*.h5)': '.h5', 'All files (*)': '' }
        fileTypeFilter = ""
        for f in fileTypeDict:
            fileTypeFilter += f + ";;"
            
        filename, filetype = qt.QFileDialog.getOpenFileName(self,"Save database",
                                                  self.filedialogdir,
                                                  fileTypeFilter[:-2])
        if filename == '':
            return
        
        self.filedialogdir = os.path.splitext(filename)[0]
        #filename += fileTypeDict[filetype]
        try:
            self.openDBFile(filename)
        except:
            msgbox = qt.QMessageBox(qt.QMessageBox.Critical,'Cannot open db file', 
            'Cannot open db file %s.' % filename,
            qt.QMessageBox.Ok, self)
            msgbox.setDetailedText(traceback.format_exc())
            clickedbutton = msgbox.exec()
        
    def saveNewDBFile(self, filename):
        alldata = nxtodict(self.nxfile)
        self.createNewDBFile(filename, alldata)
        
    def saveDBFile(self, filename):
        alldata = nxtodict(self.nxfile)
        dicttonx(alldata, filename, create_dataset_args={'compression' : self.compression})
        
    def createNewDBFile(self, filename, datadict=None):
        if self.nxfile is not None:
            self.close()

        fileattrs = {"@NX_class": u"NXroot",
                     "@creator": u"orGUI version %s" % __version__,
                     "@file_name": str(os.path.basename(filename)),
                     "@file_time": datetime.datetime.utcnow().isoformat()}
        if datadict is None:
            datadict = fileattrs
        else:
            datadict.update(fileattrs)
        try:
            dicttonx(datadict, filename, create_dataset_args={'compression' : self.compression})
            self.openDBFile(filename)
        except:
            msgbox = qt.QMessageBox(qt.QMessageBox.Critical,'Cannot create db file', 
            'Cannot create file %s.' % filename,
            qt.QMessageBox.Ok, self)
            msgbox.setDetailedText(traceback.format_exc())
            clickedbutton = msgbox.exec()
        
    def openDBFile(self, filename):
        if self.nxfile is not None:
            self.close()
        self.nxfile = silx.io.h5py_utils.File(filename,'a')
        self._filepath = filename
        while(self.hdf5model.hasPendingOperations()):
            qt.QApplication.processEvents()
            time.sleep(0.01)
        self.hdf5model.insertH5pyObject(self.nxfile)
        self.view.expandToDepth(0)

    def add_nxdict(self, nxentry, update_mode='add', h5path='/'):
        if self.nxfile is None:
            raise Exception("No database file open.")
        dicttonx(nxentry, self.nxfile, h5path=h5path, update_mode='add', create_dataset_args={'compression' : self.compression})
        while(self.hdf5model.hasPendingOperations()):
            qt.QApplication.processEvents()
            time.sleep(0.01)
        self.hdf5model.synchronizeH5pyObject(self.nxfile)
        self.view.expandToDepth(0)
        
    def onDeleteScan(self, obj):
        btn = qt.QMessageBox.question(self,"Delete scan?","Are you sure that you want to delete %s from the orgui database?" % obj.name.split("/")[-1])
        if btn == qt.QMessageBox.Yes:
            self.delete_node(obj)
            
    def onShowRoIntegrate(self, obj):
        meta = obj.h5py_object.attrs.get('orgui_meta', False)
        h5_obj = obj.h5py_object
        while(h5_obj.name != '/'): # search for rocking scan group
            meta = h5_obj.attrs.get('orgui_meta', False)
            if meta and meta == 'rocking':
                break
            h5_obj = h5_obj.parent
        else:
            msgbox = qt.QMessageBox(qt.QMessageBox.Critical, 'Invalid rocking scan', 
            'Not a rocking scan: %s.' % (dir_name),
            qt.QMessageBox.Ok, self)
            clickedbutton = msgbox.exec()
            return # invalid dataset
        
        dir_name = h5_obj.name
        self.sigChangeRockingScan.emit(dir_name)
        
    def delete_node(self, obj):
        basename = obj.name.split("/")[-1]
        objpar = obj.parent
        del objpar[basename]
        while(self.hdf5model.hasPendingOperations()):
            qt.QApplication.processEvents()
            time.sleep(0.01)
        self.hdf5model.synchronizeH5pyObject(self.nxfile)
        self.view.expandToDepth(0)
        
    def onRenameNode(self, obj):
        basename = obj.name.split("/")[-1]
        newname, success = qt.QInputDialog.getText(self,"Rename NEXUS node",
                               "New name:", qt.QLineEdit.EchoMode.Normal,
                               basename)
        if success and newname != '':
            self.rename_node(obj, newname)
                               
    def rename_node(self, obj, newname):
        basename = obj.name.split("/")[-1]
        objpar = obj.parent
        objpar.move(basename, newname)
        while(self.hdf5model.hasPendingOperations()):
            qt.QApplication.processEvents()
            time.sleep(0.01)
        self.hdf5model.synchronizeH5pyObject(self.nxfile)
        self.view.expandToDepth(0)
        
        
    def close(self):
        if self.nxfile is not None:
            timer = 1
            if self.hdf5model.hasPendingOperations():
                qt.QApplication.processEvents()
            while(self.hdf5model.hasPendingOperations() and timer < 6000):
                time.sleep(0.01)
                if not( timer % 1000):
                    qt.QApplication.processEvents()
                timer += 1
            if timer == 6000:
                raise Exception('Timeout on hdf5 model operation, This is probably a bug, or a very long writing operation occurs, please report if this is a long writing opertion')
            self.hdf5model.removeH5pyObject(self.nxfile)
            timer = 1
            while(self.hdf5model.hasPendingOperations() and timer < 6000):
                time.sleep(0.01)
                if not( timer % 1000):
                    qt.QApplication.processEvents()
                timer += 1
            if timer == 6000:
                raise Exception('Timeout on hdf5 model operation, This is probably a bug, or a very long writing operation occurs, please report if this is a long writing opertion')
            try:
                self.nxfile.close()
            except RuntimeError:
                traceback.print_exc()
                print('Closing of database file failed. The database file might be corrupted!')
            self.nxfile = None
            if hasattr(self, "temp_directory"):
                del self.temp_directory
        
    
    
