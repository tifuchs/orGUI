# -*- coding: utf-8 -*-
###############################################################################
# Copyright (c) 2021 Timo Fuchs, Olaf Magnussen all rights reserved
#
# This software was developed during the PhD work of Timo Fuchs,
# within the group of Olaf Magnussen. Usage within the group is hereby granted.
###############################################################################
"""Module descripiton

"""
__author__ = "Timo Fuchs"
__copyright__ = "Copyright 2021, Timo Fuchs, Olaf Magnussen all rights reserved"
__credits__ = []
__license__ = "all rights reserved"
__version__ = "1.0.0"
__maintainer__ = "Timo Fuchs"
__email__ = "fuchs@physik.uni-kiel.de"

import numpy as np

from .scans import Scan
from datetime import datetime
import pytz

# --- parse h5node name and return the scan number and name ---

def parseCH5523(obj):
    ddict = dict() 
    scanname = obj.local_name
    scanno = int(scanname.split('_')[-1])
    ddict['scanno'] = scanno
    ddict['name'] = obj.local_name.strip('/')
    return ddict

def parseP212H5(obj):
    ddict = dict() 
    scanname = obj.basename
    scanno, subscanno = scanname.split('.')
    ddict['scanno'] = int(scanno)
    ddict['name'] = obj.local_name
    return ddict
    
def parseID31Bliss(obj):
    ddict = dict() 
    scanname = obj.local_name
    scansuffix = scanname.split('_')[-1]
    scanname_nosuffix = '_'.join(scanname.split('_')[:-1])
    scanno, subscanno = scansuffix.split('.')
    ddict['scanno'] = int(scanno)
    ddict['name'] = obj.local_name
    return ddict
    
# orgui will search for these counters in the Scan object and copy them into the database, if available
auxillary_counters = ['current', 'potential', 'exposure_time', 'elapsed_time','time', 'srcur', 'mondio', 'epoch']
    
# assign the name parser to the beamtime identifiers:
             
scannoConverter = {'ch5523': parseCH5523,
             '20190017': parseP212H5,
             'ch5700': parseID31Bliss,
             '20200028': parseP212H5,
             'ch5918' : parseID31Bliss
             }
    
    
# --- dates of the beamtimes. Used to identify the beamtimes from the files ---

beamtimes = {'ch5523': (datetime(2018, 9, 22), datetime(2018, 10, 5)),
             '20190017': (datetime(2019, 12, 8), datetime(2019, 12, 24)),
             'ch5700': (datetime(2020, 11, 10), datetime(2020, 11, 27)),
             '20200028': (datetime(2021, 4, 27), datetime(2021, 5, 10)),
             'ch5918' : (datetime(2021, 7, 18), datetime(2021, 8, 1))
             }

def localize(dt):
    if dt.tzinfo is None or dt.tzinfo.utcoffset(dt) is None:
        grenobletime = pytz.timezone('Europe/Paris')
        return grenobletime.localize(dt)
    else:
        return dt
             
def getBeamtimeId(dt):
    for bt in beamtimes:
        start, end = beamtimes[bt]
        if localize(start) <= localize(dt) <= localize(end):
            return bt
    else:
        raise Exception("Didn't find matching beamtime for date %s" % dt.ctime())
 
            
# add actual backends here, which perform the file reads:
# They must implement scans.Scan

from datautils.xrayutils.id31_tools import BlissScan_EBS, Fastscan, BlissScan
from datautils.xrayutils.P212_tools import H5Fastsweep

fscans = {'ch5523': BlissScan,
             '20190017': H5Fastsweep, #probably doesn't work since image names are not saved
             'ch5700': BlissScan_EBS,
             '20200028': H5Fastsweep,
             'ch5918' : BlissScan_EBS
             }


def openScan(btid, ddict):
    fscancls = fscans[btid]
    
    # ideally, now the scan should only be opened with:
    # fscan = fscancls(ddict['file'], ddict['scanno'])
    # but this doesn't always work. So handle special cases here
    
    if btid == 'ch5523':
        fscan = fscancls(ddict['file'],ddict['name'])
    elif btid == '20190017' or btid == '20200028':
        fscan = fscancls(ddict['file'],ddict['scanno'])
    elif btid == 'ch5700' or btid == 'ch5918':
        if 'node' in ddict:
            fscan = fscancls(ddict['node'],ddict['scanno'], loadimg=False)
        else:
            fscan = fscancls(ddict['file'],ddict['scanno'], loadimg=False)
    else:
        try:
            fscan = fscancls(ddict['file'], ddict['scanno'])
        except Exception:
            raise ValueError("Did not find matching scan in backends for beamtime id %s" % btid)
    return fscan 


