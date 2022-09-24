import numpy as np
from datautils.xrayutils import HKLVlieg
from datautils import util

def defaultU_GID(ubCalculator):
    """Sets a default U matrix for Grazing incidence surface diffraction geometry.
    
    Geometry:
    L along z axis in alpha frame
    K along x axis in alpha frame
    for chi = 0°, phi = 0°, omega = 0° (omega = -theta)
    
    Aligns crystal L axis along omega rotation axis. 
    Alpha tilts this rotation axis. Typically alpha is equal to the angle of incidence.

    """
    # Compute the two reflections' reciprical lattice vectors in the
    # cartesian crystal frame (hc = B @ hkl)
    h1c = ubCalculator.lattice.reciprocalVectorCart([0.,0.,1.]).flatten() # for hkl = (0, 0, 1)
    h2c = ubCalculator.lattice.reciprocalVectorCart([0.,1.,0.]).flatten() # for hkl = (0, 1, 0)
    
    # Calculate the sample rotation matrices
    
    chi1 = np.deg2rad(0.)
    phi1 = np.deg2rad(0.)
    omega1 = np.deg2rad(0.)
    
    chi2 = np.deg2rad(0.)
    phi2 = np.deg2rad(0.)
    omega2 = np.deg2rad(0.)
    
    # [ALPHA, DELTA, GAMMA, OMEGA, CHI, PHI]
    pos1 = np.array([None, None, None, omega1, chi1, phi1])
    pos2 = np.array([None, None, None, omega2, chi2, phi2])
    
    _,_,_, OM1, CHI1, PHI1 = HKLVlieg.createVliegMatrices(pos1) # rotation matrices
    _,_,_, OM2, CHI2, PHI2 = HKLVlieg.createVliegMatrices(pos2) # rotation matrices
    
    # define reference directions in alpha frame
    
    Qalp1 = np.array([0.,0.,1.]) * np.linalg.norm(h1c) # reference 1:  z direction
    Qalp2 = np.array([1.,0.,0.]) * np.linalg.norm(h2c) # reference 2:  x direction
    
    # Transform Qalp in Qphi
    # hint: matrix inverse of rotation matrices is the transpose (T)
    Qphi1 = PHI1.T @ CHI1.T @ OM1.T @ Qalp1
    Qphi2 = PHI2.T @ CHI2.T @ OM2.T @ Qalp2
    
    # 1: primary vector, 2: secondary vector
    
    U, stats = HKLVlieg.UBCalculator.calc_U_from_vectors(Qphi1, Qphi2, h1c, h2c)
    
    ubCalculator.setU(U)
    
    return ubCalculator

def defaultU_TSD(ubCalculator):
    """Sets a default U matrix for Transmission Surface Diffraction geometry.
    
    Geometry:
    L along y axis in laboratory frame (along beam direction) at omega = 0°
    K along z axis in laboratory frame
    for chi = 90°
    phi = 0°
    alpha = 0°
    
    In TSD mode, the crystal L axis points towards the beam direction at one 
    specific omega value as described here:
    
    Transmission Surface Diffraction for Operando Studies of Heterogeneous Interfaces
    Finn Reikowski, Tim Wiegmann, Jochim Stettner, Jakub Drnec, Veijo Honkimäki, Fouad Maroun, Philippe Allongue, and Olaf M. Magnussen
    The Journal of Physical Chemistry Letters 2017 8 (5), 1067-1071
    DOI: 10.1021/acs.jpclett.7b00332

    """
    # Compute the two reflections' reciprical lattice vectors in the
    # cartesian crystal frame (hc = B @ hkl)
    h1c = ubCalculator.lattice.reciprocalVectorCart([0.,0.,1.]).flatten() # for hkl = (0, 0, 1)
    h2c = ubCalculator.lattice.reciprocalVectorCart([0.,1.,0.]).flatten() # for hkl = (0, 1, 0)
    
    # Calculate the sample rotation matrices
    
    chi1 = np.deg2rad(90.)
    phi1 = np.deg2rad(0.)
    omega1 = np.deg2rad(0.)
    alpha1 = np.deg2rad(0.)
    
    chi2 = np.deg2rad(90.)
    phi2 = np.deg2rad(0.)
    omega2 = np.deg2rad(0.)
    alpha2 = np.deg2rad(0.)
    
    # [ALPHA, DELTA, GAMMA, OMEGA, CHI, PHI]
    pos1 = np.array([alpha1, None, None, omega1, chi1, phi1])
    pos2 = np.array([alpha2, None, None, omega2, chi2, phi2])
    
    A1,_,_, OM1, CHI1, PHI1 = HKLVlieg.createVliegMatrices(pos1) # rotation matrices
    A2,_,_, OM2, CHI2, PHI2 = HKLVlieg.createVliegMatrices(pos2) # rotation matrices
    
    # define reference directions in laboratory frame
    
    Qlab1 = np.array([0.,1.,0.]) * np.linalg.norm(h1c) # reference 1:  y direction
    Qlab2 = np.array([0.,0.,1.]) * np.linalg.norm(h2c) # reference 2:  z direction
    
    # Transform Qlab in Qphi
    # hint: matrix inverse of rotation matrices is the transpose (T)
    Qphi1 = PHI1.T @ CHI1.T @ OM1.T @ A1.T @ Qlab1
    Qphi2 = PHI2.T @ CHI2.T @ OM2.T @ A2.T @ Qlab2
    
    # 1: primary vector, 2: secondary vector
    
    U, stats = HKLVlieg.UBCalculator.calc_U_from_vectors(Qphi1, Qphi2, h1c, h2c)
    
    ubCalculator.setU(U)
    
    return ubCalculator
    
def manualU(ubCalculator):
    U = np.array([[1., 0., 0.],
                  [0., 1., 0.],
                  [0., 0., 1.]])
    # could also do something like 
    # U = util.z_rotation(np.deg2rad(45.))
    ubCalculator.setU(U)
    return ubCalculator


u_defaults = {
    'Grazing incidence' : defaultU_GID,
    'Transmission' : defaultU_TSD,
    'manual U' : manualU
}
