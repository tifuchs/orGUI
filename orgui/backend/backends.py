# -*- coding: utf-8 -*-
# /*##########################################################################
#
# Copyright (c) 2020-2026 Timo Fuchs
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
__copyright__ = "Copyright 2020-2026 Timo Fuchs"
__license__ = "MIT License"
__version__ = "1.3.0"
__maintainer__ = "Timo Fuchs"
__email__ = "tfuchs@cornell.edu"

import numpy as np

from .scans import Scan
from datetime import datetime
import pytz
import logging
logger = logging.getLogger(__name__)

# --- parse h5node name and return the scan number and name ---

    
# --- dates of the beamtimes. Used to identify the beamtimes from the files ---

# This enables the beamtime auto detect functionality.
# However, you can always override this setting by unchecking 
# the box in the lower left corner in the gui and manually selecting
# a backend

default_beamtime = 'id31_default' 

beamtimes = {'ch5523': (datetime(2018, 9, 22), datetime(2018, 10, 5)),
             '20190017': (datetime(2019, 12, 8), datetime(2019, 12, 24)),
             'ch5700': (datetime(2020, 11, 10), datetime(2020, 11, 27)),
             '20200028': (datetime(2021, 4, 27), datetime(2021, 5, 10)),
             'ch5918' : (datetime(2021, 7, 18), datetime(2021, 8, 1)),
             'P212_default' : (datetime(1902, 7, 18), datetime(1903, 8, 1)),
             '20230930' : (datetime(2024, 7, 16), datetime(2024, 7, 25)),
             'ch7149' : (datetime(2025, 2, 17), datetime(2025, 2, 25)),
             'ch7856' : (datetime(2025, 9, 20), datetime(2025, 9, 30))#,
             #'id31_default' : (datetime(2024, 8, 2), datetime(2500, 1, 1)) # all data not found in this list automatically detected as ID31 data.
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
        return default_beamtime
        #raise Exception("Didn't find matching beamtime for date %s" % dt.ctime())
 
            
# add actual backends here, which perform the file reads:
# They must implement scans.Scan

from .beamline.id31_tools import BlissScan_EBS, BlissScan
from .beamline.P212_tools import H5Fastsweep

fscans = {'ch5523': BlissScan,
             '20190017': H5Fastsweep, #probably doesn't work since image names are not saved
             'ch5700': BlissScan_EBS,
             '20200028': H5Fastsweep,
             'ch5918' : BlissScan_EBS,
             'P212_default' : H5Fastsweep,
             'ch7149' : BlissScan_EBS,
             'ch7856' : BlissScan_EBS,
             'id31_default' : BlissScan_EBS
             }


def openScan(btid, ddict):
    fscancls = fscans[btid]
    
    # ideally, now the scan should only be opened with:
    # fscan = fscancls(ddict['file'], ddict['scanno'])
    # but this doesn't always work. So handle special cases here
    
    if btid == 'ch5523':
        fscan = fscancls(ddict['file'],ddict['name'])
        
        if ddict['name'].startswith('ascan') and 'Pt111_3' in ddict['file']:
            if fscan.axisname == 'mu':
                 mu = fscan.mu - 0.055 # misalignment! 
                 fscan.axis = mu
                 fscan.mu = mu
                 logger.warning("Correct mu misalignment 0.055 deg,  Pt111_3")

    elif btid == 'ch7149':
        if 'node' in ddict:
            fscan = fscancls(ddict['node'],ddict['scanno'], loadimg=False, muoffset=-0.544424)
        else:
            fscan = fscancls(ddict['file'],ddict['scanno'], loadimg=False, muoffset=-0.544424)
        
    elif btid == '20190017' or btid == '20200028' or btid == 'P212_default':
        fscan = fscancls(ddict['file'],ddict['scanno'])

    elif btid == 'ch5700' or btid == 'ch5918' or btid == 'id31_default':
        if 'node' in ddict:
            fscan = fscancls(ddict['node'],ddict['scanno'], loadimg=False)
        else:
            fscan = fscancls(ddict['file'],ddict['scanno'], loadimg=False)

    elif btid == 'ch7856':
        if ('pt110_prep2_naoh' in ddict['file']) and 'node' in ddict:
            fscan = fscancls(ddict['node'],ddict['scanno'], loadimg=False, muoffset=-7.381881105095823)
        elif ('pt110_prep2_naoh' in ddict['file']) and 'node' not in ddict:
            fscan = fscancls(ddict['file'],ddict['scanno'], loadimg=False, muoffset=-7.381881105095823)
        elif ('pt110_prep3_CsOH' in ddict['file'] or 'pt110_prep1_pt110_prep1_CsOH_CsClO4' in ddict['file']) and 'node' in ddict:
            fscan = fscancls(ddict['node'],ddict['scanno'], loadimg=False, muoffset=-7.254956854)
        elif ('pt110_prep3_CsOH' in ddict['file'] or 'pt110_prep1_pt110_prep1_CsOH_CsClO4' in ddict['file']) and 'node' not in ddict:
            fscan = fscancls(ddict['file'],ddict['scanno'], loadimg=False, muoffset=-7.254956854)
        elif ('Pt111' in ddict['file'] or 'pt111_prep1_CsOH' in ddict['file'] or 'pt111_prepSHIT_NaOH' in ddict['file'] or 'pt111_bs' in ddict['file']) and 'node' in ddict:
            fscan = fscancls(ddict['node'],ddict['scanno'], loadimg=False, muoffset=-5.944448934945881)
        elif ('Pt111' in ddict['file'] or 'pt111_prep1_CsOH' in ddict['file'] or 'pt111_prepSHIT_NaOH' in ddict['file'] or 'pt111_bs' in ddict['file']) and 'node' not in ddict:
            fscan = fscancls(ddict['file'],ddict['scanno'], loadimg=False, muoffset=-5.944448934945881)
        elif 'node' in ddict:
            fscan = fscancls(ddict['node'],ddict['scanno'], loadimg=False, muoffset=0)
        else:
            fscan = fscancls(ddict['file'],ddict['scanno'], loadimg=False, muoffset=0)

    else:
        try:
            fscan = fscancls(ddict['file'], ddict['scanno'])
        except Exception:
            raise ValueError("Did not find matching scan in backends for beamtime id %s" % btid)
    return fscan 


