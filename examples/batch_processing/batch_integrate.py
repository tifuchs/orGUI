# 
# This is currently a workaround and requries a bit of
# tinkering with the code of orGUI
# 
# Generally the idea is to start this script from the 
# orGUI console using the command
# 
# %run -i batch_integrate.py
# 
# this exposes the whole source code of the orGUI instance
# in form of the variables 'orgui' and 'ub' to the script.
# Both can be used to do essentially everything the gui can do.
# However, the source code was of course never intented to 
# be directly accessed by the user. So unintended results may 
# occur. It is recommended to first test the code in the console
# and subsequently transfer it to a script.

# This is an example how to select a new scannumber, and perform
# a hklscan integration (which is probably the most common usecase)

# All settings (like the backend, ub, etc.) must be set in 
# the gui and the main file opened before starting the script

orgui.autoLoadAct.setChecked(False) # disable automatic calculation of max/sum image

scans = [63, 64, 65, 66, 67, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81]
H_0 = [4, 0, 0] # you can also change these for each scannumber, or extract multiple hkl scans
H_1 = [0, 0, 1]
# xy_static = [100, 100] # pixel cordinates for static roi 

for i, sc in enumerate(scans):
    print("Integrate %s / %s: scan number %s" % (i+1, len(scans),sc)) # give status update
    orgui.scanSelector.scannoBox.setValue(sc) # select new scan number
    orgui.scanSelector._onLoadScan() # actually load scan
    
    # for hkl scan
    orgui.scanSelector.scanstab.setCurrentIndex(0) # select hkl integration (index 0)
    for j, h0 in enumerate(H_0):
        orgui.scanSelector.H_0[j].setValue(h0) # set H_0 value
    for j, h1 in enumerate(H_1):
        orgui.scanSelector.H_1[j].setValue(h1) # set H_1 value
        
    # or access static roi via
    #orgui.scanSelector.scanstab.setCurrentIndex(0) # select static integration (index 1)
    #for j, coord in enumerate(xy_static):
    #    orgui.scanSelector.xy_static[j].setValue(coord) # set xy value
    
    orgui.integrateROI() # perform integration
    

 