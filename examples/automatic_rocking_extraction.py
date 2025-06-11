import pyFAI
import os

def loadScan(scanselector,scanno):
    # load scan
    ddict = dict()
    ddict['event'] = "loadScan"
    ddict['file'] = scanselector.pathedit.text()
    ddict['scanno'] = scanno
    ddict['name'] = os.path.splitext(os.path.basename(scanselector.pathedit.text()))[0] + '.' + str(scanno)
    ddict['beamtime'] = 'ch7149'
    scanselector.sigScanChanged.emit(ddict)

def loadPoni(orgui2,ponipath):
    # load poni
    az = pyFAI.load(ponipath)
    #self._detectorDialog.selectDetector(az.detector)
    orgui2.ubcalc.machineParams.set_Xray_source({'wavelength' : az.get_wavelength()*1e10})
    model = orgui2.ubcalc.machineParams.geometryTabs.geometryModel()
    model.lockSignals()
    model.distance().setValue(az.get_dist())
    model.poni1().setValue(az.get_poni1())
    model.poni2().setValue(az.get_poni2())
    model.rotation1().setValue(az.get_rot1())
    model.rotation2().setValue(az.get_rot2())
    model.rotation3().setValue(az.get_rot3())
    #model.wavelength().setValue(az.get_wavelength())
    model.unlockSignals()
    orgui2.ubcalc.machineParams.set_detector(az.detector)
    orgui2.ubcalc.machineParams._onAnyValueChanged()

def editROI(scanselector,smax=3.,h0=[0.,0.,0.],hsize=5.,delS=0.05):
    # edit roi
    scanselector.roscanMaxS.setValue(smax)
    scanselector.roscanDeltaS.setValue(delS)

    scanselector.ro_H_0_dialog.set_hkl(h0)
    scanselector.ro_H_1_dialog.set_hkl([0.,0.,1.])

    scanselector.hsize.setValue(hsize)
    scanselector.vsize.setValue(6.)
    #scanselector.autoROIVsize = qt.QCheckBox("auto v")

def integrate(orgui2):
    #integrate
    refldict = orgui2.get_rocking_coordinates()
    if orgui2.scanSelector.intersS1Act.isChecked():
        intersect = 1
    elif orgui2.scanSelector.intersS2Act.isChecked():
        intersect = 2
    else:
        intersect = 1 # default
    mask = refldict['mask_%s' % intersect]
    xy = refldict['xy_%s' % intersect][mask]
    refldict['angles'] = refldict['angles_%s' % intersect][mask]
    roi_keys = orgui2.intbkgkeys_rocking(refldict)
    hkl_del_gam = orgui2.getStaticROIparams(xy)
    orgui2.rocking_integrate(xy, roi_keys, hkl_del_gam, refldict)


def full_rocking(orgui2,scan1,scan2):
    #orgui2.ubcalc.machineParams.azimbox.setValue(270.2)
    orgui2.scanSelector.pathedit.setText(r"E:\ch7149\RAW_DATA\waterCTR_100mM\ch7149_waterCTR_100mM.h5")
    orgui2.scanSelector.scannoBox.setValue(1)
    orgui2.scanSelector.btid.setCurrentIndex(orgui2.scanSelector.btid.findText('ch7149'))
    orgui2.scanSelector._onLoadScan()

    ponipanel = 'path1'
    poninew = 'path2'
    poniedit = 'path3'

    destination = 'path4' + str(scan1) + str(scan2) + '.h5'

    loadScan(orgui2.scanSelector,scan1)
    loadPoni(orgui2,poninew)
    editROI(orgui2.scanSelector,smax=3.,h0=[0.,0.,0.],hsize=3.,delS=0.025)
    integrate(orgui2)

    loadPoni(orgui2,poniedit)
    editROI(orgui2.scanSelector,smax=3.,h0=[0.,0.,3.],hsize=2.,delS=0.025)
    integrate(orgui2)

    editROI(orgui2.scanSelector,smax=2.5,h0=[0.,0.,6.],hsize=1.,delS=0.025)
    integrate(orgui2)

    loadScan(orgui2.scanSelector,scan2)
    loadPoni(orgui2,ponipanel)
    editROI(orgui2.scanSelector,smax=0.9,h0=[0.,0.,4.295],hsize=2.,delS=0.025)
    integrate(orgui2)

    orgui2.database.saveDBFile(destination)
