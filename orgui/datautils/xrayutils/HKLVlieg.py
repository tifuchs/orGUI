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
import numpy.linalg as LA
from .. import util
import scipy.optimize as opt
import warnings
from scipy.spatial.transform import Rotation

# needs to be multiplied with K = 2*pi/lambda to get Qphi
def calculate_q_phi(pos,K=1.):
    [ALPHA, DELTA, GAMMA, OMEGA, CHI, PHI] = createVliegMatrices(pos)

    u1a = (DELTA @ GAMMA - ALPHA.T) @ np.array([0.,K,0.])
    u1p = PHI.T @ CHI.T @ OMEGA.T @ u1a
    return u1p


"""
pos = [ALPHA, DELTA, GAMMA, OMEGA, CHI, PHI] (angles)
accepts also 1d-arrays of the angles in pos
""" 
def createVliegMatrices(pos):

    ALPHA = None if pos[0] is None else calcALPHA(pos[0])
    DELTA = None if pos[1] is None else calcDELTA(pos[1])
    GAMMA = None if pos[2] is None else calcGAMMA(pos[2])
    OMEGA = None if pos[3] is None else calcOMEGA(pos[3])
    CHI = None if pos[4] is None else calcCHI(pos[4])
    PHI = None if pos[5] is None else calcPHI(pos[5])
    return ALPHA, DELTA, GAMMA, OMEGA, CHI, PHI

def calcALPHA(alpha): # alpha = mu if GID geometry and there is no miscut 
    if isinstance(alpha,np.ndarray):
        return util.x_rotationArray(alpha)
    return util.x_rotation(alpha)


def calcDELTA(delta):
    if isinstance(delta,np.ndarray):
        return util.z_rotationArray(-delta)
    return util.z_rotation(-delta)


def calcGAMMA(gamma):
    if isinstance(gamma,np.ndarray):
        return util.x_rotationArray(gamma)
    return util.x_rotation(gamma)


def calcOMEGA(omega): 
    if isinstance(omega,np.ndarray):
        return util.z_rotationArray(-omega)
    return util.z_rotation(-omega)


def calcCHI(chi):
    if isinstance(chi,np.ndarray):
        return util.y_rotationArray(chi)
    return util.y_rotation(chi)


def calcPHI(phi):
    if isinstance(phi,np.ndarray):
        return util.x_rotationArray(phi)
    return util.x_rotation(phi)

def calcSIGMA(sigma):
    if isinstance(sigma,np.ndarray):
        return util.y_rotationArray(sigma)
    return util.y_rotation(sigma)

def calcTAU(tau):
    if isinstance(tau,np.ndarray):
        return util.y_rotationArray(-tau)
    return util.y_rotation(-tau)


def primBeamAngles(pos):
    [alpha,delta,gamma,omega,chi,phi] = pos
    gamma_p = np.arcsin( np.cos(alpha)*np.sin(gamma) + np.sin(alpha)*np.cos(delta)*np.cos(gamma) )
    delta_p = np.arcsin( (np.sin(delta)*np.cos(gamma))/np.cos(gamma_p) )
    return [alpha,delta_p,gamma_p,omega,chi,phi]

def vliegDiffracAngles(pos_p):
    [alpha,delta_p,gamma_p,omega,chi,phi] = pos_p
    gamma = np.arcsin( np.cos(alpha)*np.sin(gamma_p) - np.sin(alpha)*np.cos(delta_p)*np.cos(gamma_p) )
    delta = np.arcsin( (np.sin(delta_p)*np.cos(gamma_p))/np.cos(gamma) )
    return [alpha,delta,gamma,omega,chi,phi]
    

"""
transforms alpha and gamma into the crystal frame using 
snellius refraction law

angles smaller than 0 are treated as below sample horizon, i.e. no refraction 

"""
def crystalAngles(pos,refraction_index):
    pos = np.array(pos)
    if len(pos.shape) == 1:
        pos[0] = crystalAngles_singleArray(pos[0], refraction_index) #np.arccos(np.cos(pos[0]) / refraction_index)
        pos[2] = crystalAngles_singleArray(pos[2], refraction_index) #np.arccos(np.cos(pos[2]) / refraction_index)
        #pos[np.isnan(pos)] = 0. 
    else:
        pos[:,0] = crystalAngles_singleArray(pos[:,0], refraction_index) #np.arccos(np.cos(pos[:,0]) / refraction_index)
        pos[:,2] = crystalAngles_singleArray(pos[:,2], refraction_index) #np.arccos(np.cos(pos[:,2]) / refraction_index)
        #pos[np.isnan(pos)] = 0. 
    return pos

def crystalAngles_singleArray(angle,refraction_index):
    if isinstance(angle,np.ndarray):
        sign = np.sign(angle)
        mask = sign > 0.
        angle[mask] = np.arccos(np.cos(angle[mask]) / refraction_index)
        angle[np.isnan(angle)] = 0.
        #angle *= sign
    else:
        sign = np.sign(angle)
        if sign < 0.:
            return angle
        angle = np.arccos(np.cos(angle) / refraction_index)
        if np.isnan(angle):
            angle = 0.
    return angle
	
def vacAngles(pos,refraction_index):
    pos = np.array(pos)
    if len(pos.shape) == 1:
        pos[0] = vacAngles_singleArray(pos[0], refraction_index) #np.arccos(np.cos(pos[0]) * refraction_index)
        pos[2] = vacAngles_singleArray(pos[2], refraction_index) #np.arccos(np.cos(pos[2]) * refraction_index)
        #pos[np.isnan(pos)] = 0. 
    else:
        pos[:,0] = vacAngles_singleArray(pos[:,0], refraction_index) #np.arccos(np.cos(pos[:,0]) * refraction_index)
        pos[:,2] = vacAngles_singleArray(pos[:,2], refraction_index) #np.arccos(np.cos(pos[:,2]) * refraction_index)
        #pos[np.isnan(pos)] = 0. 
    return pos
    
def vacAngles_singleArray(angle,refraction_index):
    if isinstance(angle,np.ndarray):
        sign = np.sign(angle)
        mask = sign > 0.
        angle[mask] = np.arccos(np.cos(angle[mask]) * refraction_index)
        angle = np.nan_to_num(angle,nan=0.0)
        #angle *= sign
    else:
        sign = np.sign(angle)
        if sign < 0.:
            return angle
        angle = np.arccos(np.cos(angle) * refraction_index)
        angle = np.nan_to_num(angle,nan=0.0)
        #angle *= sign
    return angle

def printPos(pos,phichi=False):
    pos = np.rad2deg(pos)
    if phichi:
        print("alp=%.2f, del=%.2f, gam=%.2f, om=%.2f, phi=%.2f, chi=%.2f" % tuple(pos))
    else:
        print("alp=%.2f, del=%.2f, gam=%.2f, om=%.2f" % tuple(np.array((pos))[:-2]) )
        
def strPos(pos,phichi=False):
    pos = np.rad2deg(pos)
    if phichi:
        spos = "alp=%.2f, del=%.2f, gam=%.2f, om=%.2f, phi=%.2f, chi=%.2f" % tuple(pos)
    else:
        spos = "alp=%.2f, del=%.2f, gam=%.2f, om=%.2f" % tuple(np.array((pos))[:-2])
    return spos
    
def strPos_prim(pos,phichi=False):
    spos = "Vlieg: "
    spos += strPos(pos,phichi)
    spos += "  pb:"
    spos += strPos(primBeamAngles(pos),phichi)
    return spos
    
def printPos_prim(pos,phichi=False):
    print("Vlieg angles:")
    printPos(pos,phichi)
    print("angles ref prim beam")
    printPos(primBeamAngles(pos),phichi)

def spec_pa(ub):
    print(str(ub.lattice))    
    
 
"""

"""

class Lattice(object):
    """The coordiante system of the lattice is defined by first calculating the\
    reciprocal lattice's vector lengths and angles. After that, the\
    definition of the reciprocal lattice coordiante system as in\
    Busing, W. R. & Levy, H. A. Angle calculations for 3- and 4-circle\
    X-ray and neutron diffractometers. Acta Crystallogr. 22, 457–464 (1967).\
    is used:
    
    The x axis parallel to b1,
    the y axis in the plane of bl and b2, 
    and the z axis perpendicular to that plane
    
    The factor 2pi is included in the reciprocal lattice vectors.
    """
    def __init__(self,a,alpha):
        self.setLattice(a,alpha)
    
    def setLattice(self,a,alpha):
        """Sets real space lattice from lengths and angles of lattice vectors.
        
        :param np.ndarray a: lattice vector lengths in Angstrom
        :param np.ndarray alpha: angles between lattice vectors in degrees
        """
        self._alpha = np.deg2rad(alpha)
        self._a = np.asarray(a)
        self._b, self._beta = Lattice._calcReciprocalLattice(self._a, self._alpha)
        self.volume = (np.prod(self._a) *
               np.sqrt(1 + 2 * np.cos(self._alpha[0]) * np.cos(self._alpha[1]) * np.cos(self._alpha[2]) -
               np.cos(self._alpha[0]) ** 2 - np.cos(self._alpha[1]) ** 2 - np.cos(self._alpha[2]) ** 2))
        self.uc_area = self._a[0]*self._a[1]*np.sin(self._alpha[2]) #unit cell area
        
        self.B_mat = Lattice.reciprocalMatrix(self._a,self._alpha,self._b,self._beta)
        self.B_mat_inv = LA.inv(self.B_mat)
        self.R_mat = 2*np.pi * np.ascontiguousarray(self.B_mat_inv.T) 
        self.RealspaceMatrix = self.R_mat # old: Lattice.realspaceMatrix(self._a,self._alpha,self._b,self._beta)
        self.R_mat_inv = LA.inv(self.R_mat)
        
    #def setReciprocalLattice(self, b, beta):
        """Sets reciprocal space lattice.
        
        :param np.ndarray b: reciprocal lattice vector lengths in Angstrom-1. Must include the factor 2pi.
        :param np.ndarray beta: angles between lattice vectors in degrees
        """
    
    @property
    def a(self):
        return self._a
    
    @a.setter
    def a(self, a):
        self._a = np.asarray(a)
        self._b, self._beta = Lattice._calcReciprocalLattice(self._a, self._alpha)
        self.volume = (np.prod(self._a) *
               np.sqrt(1 + 2 * np.cos(self._alpha[0]) * np.cos(self._alpha[1]) * np.cos(self._alpha[2]) -
               np.cos(self._alpha[0]) ** 2 - np.cos(self._alpha[1]) ** 2 - np.cos(self._alpha[2]) ** 2))
               
        self.uc_area = self._a[0]*self._a[1]*np.sin(self._alpha[2]) #unit cell area
        
        self.B_mat = Lattice.reciprocalMatrix(self._a,self._alpha,self._b,self._beta)
        self.B_mat_inv = LA.inv(self.B_mat)
        self.R_mat = 2*np.pi * np.ascontiguousarray(self.B_mat_inv.T) 
        self.RealspaceMatrix = self.R_mat # old: Lattice.realspaceMatrix(self._a,self._alpha,self._b,self._beta)
        self.R_mat_inv = LA.inv(self.R_mat)

    @property
    def alpha(self):
        """Angle between real space lattice vectors.
        
        in radians.
        """
        return self._alpha
    
    @alpha.setter
    def alpha(self, alpha):
        self._alpha = np.asarray(alpha)
        self._b, self._beta = Lattice._calcReciprocalLattice(self._a, self._alpha)
        self.volume = (np.product(self._a) *
               np.sqrt(1 + 2 * np.cos(self._alpha[0]) * np.cos(self._alpha[1]) * np.cos(self._alpha[2]) -
               np.cos(self._alpha[0]) ** 2 - np.cos(self._alpha[1]) ** 2 - np.cos(self._alpha[2]) ** 2))
        self.uc_area = self._a[0]*self._a[1]*np.sin(self._alpha[2]) #unit cell area
        
        self.B_mat = Lattice.reciprocalMatrix(self._a,self._alpha,self._b,self._beta)
        self.B_mat_inv = LA.inv(self.B_mat)
        self.R_mat = 2*np.pi * np.ascontiguousarray(self.B_mat_inv.T) 
        self.RealspaceMatrix = self.R_mat # old: Lattice.realspaceMatrix(self._a,self._alpha,self._b,self._beta)
        self.R_mat_inv = LA.inv(self.R_mat)
    
    @property
    def beta(self):
        """Angle between reciprocal lattice vectors.
        
        in radians.
        """
        return self._beta
        
    @property
    def b(self):
        return self._b

    def getLatticeParameters(self):
        """Legacy. Use Lattice.a and Lattice.alpha instead
        """
        return self._a, self._alpha, self._b, self._beta
    
    def getB(self):
        """Legacy. Use Lattice.B_mat instead
        """
        return self.B_mat
    
    #in rad, shape of hkl must be either (3,) or (3,n)
    def get2ThetaFromHKL(self,hkl,energy):
        """Calculates scattering angle 2theta for given hkl.
        
        From |B_mat @ H| = |Q| = (4pi / lambda) sin(theta) 
        """
        wavelength = 12.39842 / energy
        hkl = np.array(hkl)
        if len(hkl.shape) == 1:
            G = LA.norm(self.reciprocalVectorCart(hkl).T)
        else:
            G = LA.norm(self.reciprocalVectorCart(hkl).T,axis=1)
        return 2*np.arcsin((G*wavelength)/(4*np.pi))

    
    # calculates atomic positions from fractional coordinates xyz_frac 
    # to cartesian coordinates in Angstroms  
    def directVectorCart(self,xyz_frac):
        """Calculates real space vector in cartesian coordinates in crystal frame\
        from fractional coordinates xyz_frac.
        
        See definition of the reciprocalMatrix B_mat. The 

        in Angstrom
        """
        return self.R_mat @ np.asarray(xyz_frac).T
    
    def reciprocalVectorCart(self,hkl):
        """Calculates reciprocal vector in cartesian coordinates in crystal frame\
        from reciprocal lattice units hkl.
        
        in Angstrom-1
        """
        return self.B_mat @ np.asarray(hkl).T
        
    # calculates reciprocal vector in cartesian coordinates 
    # from lattice units hkl
    def getReciprocalVectorCart(self,hkl):
        """deprecated. use reciprocalVectorCart instead.
        """
        warnings.warn("getReciprocalVectorCart is deprecated. Use reciprocalVectorCart instead.", FutureWarning)
        return self.B_mat @ np.asarray(hkl).T
    
    @staticmethod
    def _calcReciprocalLattice(a, alpha):
        beta = np.empty(3)
        
        beta[0] = np.arccos((np.cos(alpha[1]) * np.cos(alpha[2]) - np.cos(alpha[0])) /
                           (np.sin(alpha[1]) * np.sin(alpha[2])))

        beta[1] = np.arccos((np.cos(alpha[0]) * np.cos(alpha[2]) - np.cos(alpha[1])) /
                                   (np.sin(alpha[0]) * np.sin(alpha[2])))

        beta[2] = np.arccos((np.cos(alpha[0]) * np.cos(alpha[1]) - np.cos(alpha[2])) /
                                   (np.sin(alpha[0]) * np.sin(alpha[1])))
        
        b = np.empty(3)
        
        volume = (np.product(a) *
          np.sqrt(1 + 2 * np.cos(alpha[0]) * np.cos(alpha[1]) * np.cos(alpha[2]) -
               np.cos(alpha[0]) ** 2 - np.cos(alpha[1]) ** 2 - np.cos(alpha[2]) ** 2))
        
        b[0] = 2 * np.pi * a[1] * a[2] * np.sin(alpha[0]) / volume
        b[1] = 2 * np.pi * a[0] * a[2] * np.sin(alpha[1]) / volume
        b[2] = 2 * np.pi * a[0] * a[1] * np.sin(alpha[2]) / volume
        
        return b, beta
    
    @staticmethod
    def reciprocalMatrix(a,alpha,b,beta):
        """Creates a matrix that transforms reciprocal lattice coordinates into cartesian coordinates.
        h along x, k 
        
        """
        return np.array([(b[0],     b[1]*np.cos(beta[2]),   b[2]*np.cos(beta[1])),
                         (0.,        b[1]*np.sin(beta[2]),   -b[2]*np.sin(beta[1])*np.cos(alpha[0])),
                         (0.,        0.,                     2.*np.pi/a[2])])
    @staticmethod 
    def realspaceMatrix(a,alpha,b,beta):
        """Creates a matrix that transforms lattice coordinates into cartesian coordinates.
        
        legacy version! This matrix is not the reciprocal to reciprocalMatrix! 
        
        """
        return np.array([(a[0],     a[1]*np.cos(alpha[2]),   a[2]*np.cos(alpha[1])),
                         (0.,        a[1]*np.sin(alpha[2]),   -a[2]*np.sin(alpha[1])*np.cos(beta[0])),
                         (0.,        0.,                     2.*np.pi/b[2])])

    
    
    def __str__(self):
        name = "Lattice:\nReal space: \t" + str(self._a) + " / " + str(np.rad2deg(self._alpha)) + "\n"
        name += "Reciprocal space: \t" + str(self._b) + " / " + str(np.rad2deg(self._beta))
        return name


    def __repr__(self):
        name = "Lattice:\nReal space: \t" + str(self._a) + " / " + str(np.rad2deg(self._alpha)) + "\n"
        name += "Reciprocal space: \t" + str(self._b) + " / " + str(np.rad2deg(self._beta))
        return name

class Crystal(Lattice):
    def __init__(self, *args, **kwargs):
        warnings.warn("Crystal is deprecated. Use Lattice instead.", FutureWarning)
        super().__init__(*args, **kwargs)

class UBCalculator():
    
    def __init__(self,lattice, energy):
        self._U = None
        self._UB = None
        self.setLattice(lattice)
        self.setEnergy(energy)
    
    def setLattice(self, lattice):
        self.lattice = lattice
        if self._U is not None:
            self._UB = self._U @ self.lattice.B_mat
        
    def setEnergy(self,energy):
        """in keV
        
        """
        self._energy = energy
        self._lambda = 12.39842 / energy
        self._K = (2*np.pi)/self._lambda
        
    def setLambda(self,lmbda):
        """in Angstrom
        
        """
        self._energy = 12.39842 / lmbda
        self._lambda = lmbda
        self._K = (2*np.pi)/self._lambda
    
    def getLattice(self):
        """depricated. use UBCalculator.lattice instead.
        """
        return self.lattice
    
    def getEnergy(self):
        return self._energy
    
    def getLambda(self):
        return self._lambda
    
    def getK(self):
        return self._K
    
    def setPrimaryReflection(self,pos,hkl):
        self._primary = (pos,hkl)
        
    def setSecondayReflection(self,pos,hkl):
        self._secondary = (pos,hkl)
        
    # l in z-direction
    def defaultU(self):
        """For compatibility of legacy code. Use one of the explicit defaults for future developments instead. 
        """
        self.defaultU_GID()
        """
        TwoTheta1 = self.lattice.get2ThetaFromHKL([0,0,1],self._energy)
        TwoTheta2 = self.lattice.get2ThetaFromHKL([0,1,0],self._energy)
        self.setPrimaryReflection([TwoTheta1/2.,0.,TwoTheta1/2.,0,0.,0.],[0,0,1])
        self.setSecondayReflection([0.,TwoTheta2,0.,0.,0.,0.],[0,1,0])
        self.calculateU()
        """
        
    def alignU_lab(self, hkl, pos, orientation):
        """Rotate current U matrix in the shortest way,
        so that the reflection ´hkl´ with the angles ´pos´
        points along ´orientation´ in lab coordinates.
        
        :param hkl np.ndarray: Miller indices of reflection 
        :param pos np.ndarray: sixc diffractometer angles [ALPHA, DELTA, GAMMA, OMEGA, CHI, PHI] in rad
        :param orientation np.ndarray: direction to rotate in laboratory coordinates
        """
        ALPHA, _, _, OMEGA, CHI, PHI = createVliegMatrices(pos)
        
        Vlab = ALPHA @ OMEGA @ CHI @ PHI
        Qphi = np.asarray(self.getUB()) @ hkl
        Qlab = Vlab @ Qphi
        
        Qlab_norm = Qlab / np.linalg.norm(Qlab)
        orientation_norm = orientation / np.linalg.norm(orientation)
        
        rotation = util.rotation_matrix_from_vectors(Qlab_norm, orientation_norm)
        
        Unew = Vlab.T @ rotation @ Vlab @ np.asarray(self.getU())
        
        self.setU(Unew)
        
    def alignU_alpha(self, hkl, pos, orientation):
        """Rotate current U matrix in the shortest way,
        so that the reflection ´hkl´ with the angles ´pos´
        points along ´orientation´ in lab coordinates.
        
        :param hkl np.ndarray: Miller indices of reflection 
        :param pos np.ndarray: sixc diffractometer angles [ALPHA, DELTA, GAMMA, OMEGA, CHI, PHI] in rad
        :param orientation np.ndarray: direction to rotate in laboratory coordinates
        """
        _, _, _, OMEGA, CHI, PHI = createVliegMatrices(pos)
        
        Valp = OMEGA @ CHI @ PHI
        Qphi = np.asarray(self.getUB()) @ hkl
        Qalp = Valp @ Qphi
        
        Qalp_norm = Qalp / np.linalg.norm(Qalp)
        orientation_norm = orientation / np.linalg.norm(orientation)
        
        rotation = util.rotation_matrix_from_vectors(Qalp_norm, orientation_norm)
        
        Unew = Valp.T @ rotation @ Valp @ np.asarray(self.getU())
        
        self.setU(Unew)
    
    def defaultU_GID(self):
        """Sets a default U matrix for Grazing incidence surface diffraction geometry.
        
        Geometry:
        L along z axis in alpha frame
        K along x axis in alpha frame
        for chi = 0°, phi = 0°
        
        Aligns crystal L axis along omega rotation axis. 
        Alpha tilts this rotation axis. Typically alpha is equal to the angle of incidence.

        """
        # Compute the two reflections' reciprical lattice vectors in the
        # cartesian crystal frame (hc = B @ hkl)
        h1c = self.lattice.reciprocalVectorCart([0.,0.,1.]).flatten() # for hkl = (0, 0, 1)
        h2c = self.lattice.reciprocalVectorCart([0.,1.,0.]).flatten() # for hkl = (0, 1, 0)
        
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
        
        _,_,_, OM1, CHI1, PHI1 = createVliegMatrices(pos1) # rotation matrices
        _,_,_, OM2, CHI2, PHI2 = createVliegMatrices(pos2) # rotation matrices
        
        # define reference directions in alpha frame
        
        Qalp1 = np.array([0.,0.,1.]) * np.linalg.norm(h1c) # reference 1:  z direction
        Qalp2 = np.array([1.,0.,0.]) * np.linalg.norm(h2c) # reference 2:  x direction
        
        # Transform Qalp in Qphi
        # hint: matrix inverse of rotation matrices is the transpose (T)
        Qphi1 = PHI1.T @ CHI1.T @ OM1.T @ Qalp1
        Qphi2 = PHI2.T @ CHI2.T @ OM2.T @ Qalp2
        
        # 1: primary vector, 2: secondary vector
        
        U, stats = UBCalculator.calc_U_from_vectors(Qphi1, Qphi2, h1c, h2c)
        
        self.setU(U)
        
    def defaultU_TSD(self):
        """Sets a default U matrix for Transmission Surface Diffraction geometry.
        
        Geometry:
        L along y axis in laboratory frame (along beam direction) at omega = 0°
        K along z axis in laboratory frame
        for phi = 90°
        chi = 0°
        
        In TSD mode, the crystal L axis points towards the beam direction at one 
        specific omega value as described here:
        
        Transmission Surface Diffraction for Operando Studies of Heterogeneous Interfaces
        Finn Reikowski, Tim Wiegmann, Jochim Stettner, Jakub Drnec, Veijo Honkimäki, Fouad Maroun, Philippe Allongue, and Olaf M. Magnussen
        The Journal of Physical Chemistry Letters 2017 8 (5), 1067-1071
        DOI: 10.1021/acs.jpclett.7b00332

        """
        # Compute the two reflections' reciprical lattice vectors in the
        # cartesian crystal frame (hc = B @ hkl)
        h1c = self.lattice.reciprocalVectorCart([0.,0.,1.]).flatten() # for hkl = (0, 0, 1)
        h2c = self.lattice.reciprocalVectorCart([0.,1.,0.]).flatten() # for hkl = (0, 1, 0)
        
        # Calculate the sample rotation matrices
        
        chi1 = np.deg2rad(00.)
        phi1 = np.deg2rad(90.)
        omega1 = np.deg2rad(0.)
        alpha1 = np.deg2rad(0.)
        
        chi2 = np.deg2rad(0.)
        phi2 = np.deg2rad(90.)
        omega2 = np.deg2rad(0.)
        alpha2 = np.deg2rad(0.)
        
        # [ALPHA, DELTA, GAMMA, OMEGA, CHI, PHI]
        pos1 = np.array([alpha1, None, None, omega1, chi1, phi1])
        pos2 = np.array([alpha2, None, None, omega2, chi2, phi2])
        
        A1,_,_, OM1, CHI1, PHI1 = createVliegMatrices(pos1) # rotation matrices
        A2,_,_, OM2, CHI2, PHI2 = createVliegMatrices(pos2) # rotation matrices
        
        # define reference directions in laboratory frame
        
        Qlab1 = np.array([0.,1.,0.]) * np.linalg.norm(h1c) # reference 1:  y direction
        Qlab2 = np.array([0.,0.,1.]) * np.linalg.norm(h2c) # reference 2:  z direction
        
        # Transform Qlab in Qphi
        # hint: matrix inverse of rotation matrices is the transpose (T)
        Qphi1 = PHI1.T @ CHI1.T @ OM1.T @ A1.T @ Qlab1
        Qphi2 = PHI2.T @ CHI2.T @ OM2.T @ A2.T @ Qlab2
        
        # 1: primary vector, 2: secondary vector
        
        U, stats = UBCalculator.calc_U_from_vectors(Qphi1, Qphi2, h1c, h2c)
        
        self.setU(U)
        
        
    def zmodeUSingleRefl(self,pos,hkl):
        """
        assumes L pointing in z-direction

        Parameters
        ----------
        pos : TYPE
            pos = [ALPHA, DELTA, GAMMA, OMEGA, CHI, PHI] (angles)
        hkl : TYPE
            DESCRIPTION.

        Returns
        -------
        None.

        """
        TwoTheta1 = self.lattice.get2ThetaFromHKL([0,0,1],self._energy)
        self.setPrimaryReflection([TwoTheta1/2.,0.,TwoTheta1/2.,0,0.,0.],[0,0,1])
        self.setSecondayReflection(pos,hkl)
        self.calculateU()
        """
        U = np.asarray(self.getU())
        U[0,2] = 0
        U[1,2] = 0
        U[2,0] = 0
        U[2,1] = 0
        eigenvalue = np.linalg.eigvals(U) # make matrix unitary
        evreal = eigenvalue[~np.iscomplex(eigenvalue)]
        if len(evreal.shape) > 0:
            evreal = evreal[0]
        U /= np.real(evreal)  
        self.setU(U)
        """
    
    @staticmethod
    def calc_U_from_vectors(Q_phi_1, Q_phi_2, Q_c_1, Q_c_2):
        """Calculate the orientation matrix ´U´ from two reference reflections.
        
        Calculates a orthogonal matrix U that closely satisfies:

        Q_phi_1 = U @ Q_c_1
        Q_phi_2 = U @ Q_c_2
        
        The procedure is described in
        Busing, W. R. & Levy, H. A. Angle calculations for 3- and 4-circle X-ray and neutron diffractometers. Acta Crystallogr. 22, 457–464 (1967).
        
        Q_phi_1, Q_phi_2 are the momentum transfer vector in the phi frame 
            (attached to the diffractometer at all angles = 0)
        
        Q_c_1, Q_c_2 are the momentum transfer vectors in the crystal system
        
        reflection vectors Q_phi_1, Q_c_1 are the primary vectors, i.e. they define a fixed reference,
        and the secondary vectors Q_phi_2, Q_c_2 only lead to a rotation around the primary vector.
        """
        
        # Normalize vectors as we are only interested in orientations
        
        u_phi_1 = Q_phi_1 / np.linalg.norm(Q_phi_1) # could also normalize after doing cross products...
        u_phi_2 = Q_phi_2 / np.linalg.norm(Q_phi_2)
        
        h_c_1 = Q_c_1 / np.linalg.norm(Q_c_1)
        h_c_2 = Q_c_2 / np.linalg.norm(Q_c_2)
        
        vector_norm_mismatch_1 = np.linalg.norm(Q_phi_1) / np.linalg.norm(Q_c_1)
        vector_norm_mismatch_2 = np.linalg.norm(Q_phi_2) / np.linalg.norm(Q_c_2)
        
        # Create modified, orthogonal unit vectors t_c and t_phi in crystal and phi frame, respectively
        
        # t_c_1 parallel to h_c_1, t_c_2 lies in the plane of h_c_1 and h_c_2, and t_c_3 is perpendicular to this plane.
        t_c_1 = h_c_1
        t_c_3 = np.cross(h_c_1, h_c_2)
        t_c_2 = np.cross(t_c_3, t_c_1)
        
        vector_angle_c = np.arcsin(np.linalg.norm(t_c_3))
        
        t_p_1 = u_phi_1
        t_p_3 = np.cross(u_phi_1, u_phi_2)
        t_p_2 = np.cross(t_p_3, t_p_1)
        
        vector_angle_p = np.arcsin(np.linalg.norm(t_p_3))
        
        vector_angle_mismatch = np.abs(vector_angle_c - vector_angle_p)
        vector_angle = np.mean([vector_angle_c, vector_angle_p])
        
        if np.rad2deg(vector_angle) < 5.:
            raise Exception("Angle between reference vectors < 5°, this would result in a bad orientation matrix")
        
        t_c_3 /= np.linalg.norm(t_c_3)
        t_c_2 /= np.linalg.norm(t_c_2)
        
        t_p_3 /= np.linalg.norm(t_p_3)
        t_p_2 /= np.linalg.norm(t_p_2)
        
        Tc = np.column_stack([t_c_1, t_c_2, t_c_3])
        Tp = np.column_stack([t_p_1, t_p_2, t_p_3])
        
        U = Tp @ Tc.T
        
        stats = {'angle' : vector_angle,
                 'angle_mismatch' : vector_angle_mismatch,
                 'norm_mismatch_1' : vector_norm_mismatch_1,
                 'norm_mismatch_2' : vector_norm_mismatch_2}
        
        return U, stats
    
    def calculateU(self):
        """Calculate the orientation matrix ´U´ from the two reference reflections.
        
        The procedure is described in
        Busing, W. R. & Levy, H. A. Angle calculations for 3- and 4-circle X-ray and neutron diffractometers. Acta Crystallogr. 22, 457–464 (1967).

        """
        ppos, phkl =  self._primary
        spos, shkl =  self._secondary
        
        # Compute the two reflections' reciprical lattice vectors in the
        # cartesian crystal frame (hc = B * hkl)
        h1c = self.lattice.reciprocalVectorCart(phkl).flatten()
        h2c = self.lattice.reciprocalVectorCart(shkl).flatten()
        
        # Calculate observed vectors in the phi frame. 
        u1p = calculate_q_phi(ppos, self._K).flatten()
        u2p = calculate_q_phi(spos, self._K).flatten()
        
        U, stats = UBCalculator.calc_U_from_vectors(u1p, u2p, h1c, h2c)
        
        self.setU(U)
        return self._U

    def getU(self):
        return self._U 
        
    def getUB(self):
        return self._UB

    def getUmB(self):
        B = self.lattice.B_mat
        return  self._U @ B
    
    def setU(self,U):
        U = np.asarray(U).reshape((3,3))
        self._U = np.ascontiguousarray(U)
        self._UB = np.ascontiguousarray(self._U @ self.lattice.getB())
        
    # this fits the U matrix and also fits the lattice constants 
    # you can either scale the lattice constants equally or fit them 
    # individually
    # mode: either 'scale' or 'indivdual'
    def refineULattice(self,hkl,angles,mode='scale',rod=None,factor=100.):
        qphi = []
        for pos in angles:
            qphi.append(calculate_q_phi(pos,self._K).T)
        qphi = np.array(qphi).T
        # p[0], p[1], p[2]: rotation angles
        
        weights = np.ones(angles.shape[0])
        
        if rod is not None:
            weights[np.all(hkl[:,:2] == rod,axis=1)] = factor
        
        #a,alpha,_,_ = self.lattice.getLatticeParameters()
        a = self.lattice.a
        alpha = np.rad2deg(self.lattice.alpha)
        a = np.array(a)
        if mode=='scale':
            def Chi2(p):
                self.lattice.setLattice(a*p[3],alpha)
                UBnew = util.x_rotation(p[0]) @ util.y_rotation(p[1]) @ util.z_rotation(p[2]) @ self._U @ self.lattice.B_mat                
                hklnew = (np.linalg.inv(UBnew) @ qphi).T
                return np.sum(LA.norm(hkl - hklnew,axis=1) * weights)
            res = opt.minimize(Chi2,[0,0,0,1.])
        else:
            
            def Chi2(p):
                self.lattice.setLattice(a*np.array(p[3:]),alpha)
                UBnew = util.x_rotation(p[0]) @ util.y_rotation(p[1]) @ util.z_rotation(p[2]) @ self._U @ self.lattice.B_mat
                hklnew = (np.linalg.inv(UBnew) @ qphi).T
                return np.sum(LA.norm(hkl - hklnew,axis=1) * weights)
            res = opt.minimize(Chi2,[0,0,0,1.,1.,1.])
        #print(res)
        U = util.x_rotation(res.x[0]) @ util.y_rotation(res.x[1]) @ util.z_rotation(res.x[2]) @ self._U
        self.setU(U)
    
    
    def refineU(self,hkl,angles,allowPhiChi_opt=False,rod=None,factor=100.):
        qphi = []
        for pos in angles:
            qphi.append(calculate_q_phi(pos,self.getK()).T)
        qphi = np.array(qphi).T
        # p[0], p[1], p[2]: rotation angles
        
        weights = np.ones(angles.shape[0])
        
        if rod is not None:
            weights[np.all(hkl[:,:2] == rod,axis=1)] = factor
        
        if allowPhiChi_opt:
            def Chi2(p):
                qphi = []
                UBnew = util.x_rotation(p[0]) @ util.y_rotation(p[1]) @ util.z_rotation(p[2]) @ self._U @ self.lattice.B_mat
                #if input("Type exit") == "exit": raise Exception("bla")
                for pos in angles:
                    pos[4] = p[3]
                    pos[5] = p[4]
                    qphi.append(calculate_q_phi(pos,self.getK()).T)
                qphi = np.array(qphi).T
                
                hklnew = (np.linalg.inv(UBnew) @ qphi).T
                #print(np.sum(LA.norm(hkl - hklnew,axis=1)))
                return np.sum(LA.norm(hkl - hklnew,axis=1) * weights )
            res = opt.minimize(Chi2,[0,0,0,0,0])
        else:
            
            def Chi2(p):
                UBnew = util.x_rotation(p[0]) @ util.y_rotation(p[1]) @ util.z_rotation(p[2]) @ self._U @ self.lattice.B_mat
                hklnew = (np.linalg.inv(UBnew) @ qphi).T
                return np.sum(LA.norm(hkl - hklnew,axis=1) * weights)
                
            res = opt.minimize(Chi2,[0,0,0])
        U = util.x_rotation(res.x[0]) @ util.y_rotation(res.x[1]) @ util.z_rotation(res.x[2]) @ self._U
        self.setU(U)

    def bruteForceU(self,hkl,angles):
        qphi = []
        for pos in angles:
            qphi.append(calculate_q_phi(pos,self.getK()).T)
        qphi = np.array(qphi).T
        # p[0], p[1], p[2]: rotation angles
        #I = np.matrix([[1.,0.,0.],[0.,1.,0.],[0.,0.,1.]])
        def Chi2(p):
            UBnew = util.x_rotation(p[0]) @ util.y_rotation(p[1]) @ util.z_rotation(p[2]) @ self.lattice.B_mat
            hklnew = (np.linalg.inv(UBnew) * qphi).T
            #print(hklnew)
            #print(hkl)
            #print(np.sum(LA.norm(hkl - hklnew,axis=1)))
            return np.sum(LA.norm(hkl - hklnew,axis=1))
        res = opt.minimize(Chi2,[0,0,0])
        U = util.x_rotation(res.x[0]) @ util.y_rotation(res.x[1]) @ util.z_rotation(res.x[2])
        self.setU(U)
    
    def __str__(self):
        pstr = 'E = ' + str(self.getEnergy()) + ' keV, lambda = ' +str(self.getLambda()) + "\n"
        pstr += str(self.lattice)
        if hasattr(self, "_primary"):
            ppos, phkl =  self._primary
            pstr += '\nprimary reflection (or0): %s\n%s\n' % (phkl , strPos_prim(ppos) )
        if hasattr(self, "_secondary"):
            spos, shkl =  self._secondary 
            pstr += 'secondary reflection (or1): %s\n%s\n' % (shkl , strPos_prim(spos) )
        return pstr
    
    def __repr__(self):
        pstr = 'E = ' + str(self.getEnergy()) + ' keV, lambda = ' +str(self.getLambda()) + "\n"
        pstr += str(self.lattice)
        if hasattr(self, "_primary"):
            ppos, phkl =  self._primary
            pstr += '\nprimary reflection (or0): %s\n%s\n' % (phkl , strPos_prim(ppos) )
        if hasattr(self, "_secondary"):
            spos, shkl =  self._secondary 
            pstr += 'secondary reflection (or1): %s\n%s\n' % (shkl , strPos_prim(spos) )
        return pstr
        
    
    


class VliegAngles():
    def __init__(self,ubCalculator):
        self._ubCalculator = ubCalculator
        
    """
    Returns the hkl values of a aingle detector frame. 
    phi,chi,alpha,omega are fixed, gamma and delta are 1-d arrays
    The calculation is optimized with numpy
    Only kinematical calculation!
    """
    def anglesToHklDetector(self,alpha,delta,gamma,omega,chi,phi):
        [ALPHA, DELTA, GAMMA, OMEGA, CHI, PHI] = createVliegMatrices([alpha,delta,gamma,omega,chi,phi])
        hkl = np.empty((gamma.size,delta.size,3))
        K = self._ubCalculator.getK()
        UBi = np.linalg.inv(self._ubCalculator.getUB())
        ALPHAi = ALPHA.T
        OMEGAi = OMEGA.T
        CHIi = CHI.T
        PHIi = PHI.T
        for i in range(gamma.size):
            #calculate ( DELTA * GAMMA - ALPHA**-1 ) * K_lab = Q_alpha
            DEL_GAM_minALP = np.matmul((np.matmul(DELTA,GAMMA[i]) - ALPHAi), np.array([0.,K,0.]) ).T
            #calculate UBi * PHIi * CHIi * OMEGAi * Q_alpha 
            hkl[i] = np.matmul(UBi,np.matmul(PHIi,np.matmul(CHIi,np.matmul(OMEGAi,DEL_GAM_minALP)))).T
        return hkl[:,:,0], hkl[:,:,1], hkl[:,:,2] # h k l
        
    def anglesToHklDetector_mesh(self,alpha,delta,gamma,omega,chi,phi):
        warnings.warn("anglesToHklDetector_mesh is deprecated and will be removed in the future. Use anglesToHkl instead.", FutureWarning)
        return self.anglesToHkl(alpha,delta,gamma,omega,chi,phi)

    
    def anglesToHkl(self,alpha,delta,gamma,omega,chi,phi):
        delta = np.asarray(delta)
        gamma = np.asarray(gamma)
        assert delta.shape == gamma.shape
        
        [_, _, _, OMEGA, CHI, PHI] = createVliegMatrices([None,None,None,omega,chi,phi])
        #hkl = np.empty((*shape,3))
        K = self._ubCalculator.getK()
        UBi = np.linalg.inv(self._ubCalculator.getUB())

        OMEGAi = OMEGA.T
        CHIi = CHI.T
        PHIi = PHI.T

        Qalp = self.QAlpha(alpha,delta,gamma)
        shape = Qalp.shape
        Qalp = np.ascontiguousarray(Qalp.reshape((-1,3)).T)
        #calculate UBi * PHIi * CHIi * OMEGAi * Q_alpha 

        hkl = (UBi @ PHIi @ CHIi @ OMEGAi) @ Qalp
        #hkl = np.matmul(UBi,np.matmul(PHIi,np.matmul(CHIi,np.matmul(OMEGAi,DEL_GAM_minALP.T)))).T.reshape((*shape,3))
        hkl = hkl.T.reshape(shape)
        return hkl[...,0], hkl[...,1], hkl[...,2] # h k l

    # only for single points
    def __anglesToHkl(self,pos):
        """
        Returns hkl from pos object in radians.
        Only kinematical calculation!
        """
        pos = np.array(pos)
        if len(pos.shape) == 1:
            return np.linalg.inv(self._ubCalculator.getUB()) @ calculate_q_phi(pos,self._ubCalculator.getK())
        else:
            qphi = []
            for p in pos:
                qphi.append(calculate_q_phi(p,self._ubCalculator.getK()).T)
            qphi = np.array(qphi).T
            return (np.linalg.inv(self._ubCalculator.getUB()) @ qphi).T

            
    def QAlpha(self,alpha,delta,gamma):
        delta = np.asarray(delta)
        gamma = np.asarray(gamma)
        shape = delta.shape
        assert delta.shape == gamma.shape
        Qxyz = np.empty((*shape,3),dtype=np.float64)
        
        cosgam = np.cos(gamma) 
        Qxyz[...,0] = - np.sin(-delta)*cosgam
        Qxyz[...,1] = np.cos(-delta)*cosgam - np.cos(alpha)
        Qxyz[...,2] = np.sin(gamma) + np.sin(alpha)
        
        Qxyz *= self._ubCalculator.getK()
        return Qxyz
        
        
    def QAlphaDetector(self,alpha,delta,gamma):
        [ALPHA, DELTA, GAMMA, OMEGA, CHI, PHI] = createVliegMatrices([alpha,delta,gamma,None,None,None])
        K = self._ubCalculator.getK()
        Qxyz = np.empty((gamma.size,delta.size,3))
        ALPHAi = ALPHA.T
        for i in range(gamma.size):
            #calculate ( DELTA * GAMMA - ALPHA**-1 ) * K_lab = Q_alpha
            Qxyz[i] = np.matmul((np.matmul(DELTA,GAMMA[i]) - ALPHAi), np.array([0.,K,0.]) )
        return Qxyz[:,:,0], Qxyz[:,:,1], Qxyz[:,:,2] # Qx Qy Qz
    
    def anglesZmode(self,hkl,fixedangle,fixed='in',chi=0,phi=0,**keyargs):
        """
        Calculates the diffractometer angles for the z-mode of the
        the 6-circle diffractometer as described in Lohmeier and Vlieg 1993.
        
        Currently only fully functioning in forward scattering.
        
        Parameters
        ----------
        hkl: array of shape (3,n) or array_like of shape (3,)
            reciprocal space coordinates in reciprocal lattice units 
            of the lattice specified in the given UBCalculator
            
        fixedangle: float
            value of the fixed angle in rad, 
            has no effect if equal angle constraint is given.
            
        fixed: str (default: 'in')
            specifies one of three possible angle constraints:
            'in': angle of incidence (alpha) fixed
            'out': exit angle fixed (gamma)
            'eq': equal incident and exit angles. 
                This typically is used to measure the 
                specular reflectivity or the specular CTR
                The value of fixedangle has no effect, pass any number!
        
        chi, phi: float (default: 0)
            values of the fixed inner sample circles chi and phi in rad
            
        **keyargs:
            'mirrorx' : bool (default : False)
            indicates that the reflections with negative delta should be
            calculated. 
            
        Returns
        -------
        out: diffractometer angles
        alpha, delta, gamma, omega, chi, phi in rad
        if hkl was of shape (3,):
            6-tuple of the angles
        if hkl was of shape (3,n):
            ndarray of shape (n,6) with the angles
        """
        [ALPHA, DELTA, GAMMA, OMEGA, CHI, PHI] = createVliegMatrices([None,None,None,None,chi,phi])
        hkl = np.array(hkl)
        K = self._ubCalculator.getK()
        Hphi = np.matmul(self._ubCalculator.getUB(),hkl)
        
        Homega = CHI @ PHI @ Hphi
        if fixed == 'in':
            alpha = fixedangle
            gamma = np.arcsin(Homega[2]/K - np.sin(alpha) )
        elif fixed == 'out':
            gamma = fixedangle
            alpha = np.arcsin(Homega[2]/K - np.sin(gamma) )            
        elif fixed == 'eq':
            gamma = alpha = np.arcsin(Homega[2]/(2*K))
        else:
            raise Exception("No valid angle constraint given. Should be one of 'in', 'out' or 'eq'")
        if len(hkl.shape) > 1:
            accos_arg = (1. - np.sum(Homega**2,axis=0) / (2*K**2) + np.sin(gamma)*np.sin(alpha)) * \
                              (1 / (np.cos(gamma)*np.cos(alpha)))
            mask = np.logical_or(accos_arg > 0.99, accos_arg < -0.99) # handle numerical precision issue close to arccos(1). 
            accos_arg[mask] = np.round(accos_arg[mask], np.finfo(np.float64).precision - 2)
            delta = np.arccos(accos_arg)

            gamma = np.full_like(delta,gamma) if np.array(gamma).size < 2 else gamma
            alpha = np.full_like(delta,alpha) if np.array(alpha).size < 2 else alpha
            chi = np.full_like(delta,chi)
            phi = np.full_like(delta,phi)
        else:
            accos_arg = (1. - np.sum(Homega**2) / (2*K**2) + np.sin(gamma)*np.sin(alpha)) * \
                              (1 / (np.cos(gamma)*np.cos(alpha)))
            if accos_arg > 0.99 or accos_arg < -0.99:
                accos_arg = np.round(accos_arg, np.finfo(np.float64).precision - 2)
            delta = np.arccos(accos_arg)
        if keyargs.get('mirrorx',False):
            delta = - delta # solution since delta is calculated using arccos above
        omega = np.arctan2((Homega[1]*np.sin(delta)*np.cos(gamma) - Homega[0]*(np.cos(delta)*np.cos(gamma) - np.cos(alpha))),
                          (Homega[0]*np.sin(delta)*np.cos(gamma) + Homega[1]*(np.cos(delta)*np.cos(gamma) - np.cos(alpha))))
        
        # omega now in range +-pi
        # convert to 0 - 2*pi, since most diffractometers work from 0 to 360 degrees
        #omega = (omega + np.pi) % np.pi
        if len(hkl.shape) > 1:
            return np.vstack((alpha, delta, gamma, omega, chi, phi)).T
        else:
            return alpha, delta, gamma, omega, chi, phi
        
    def intersectLineEwald(self, H_0, H_1, alpha,omega,phi=0.,chi=0.,**keyargs):
        shape = H_0.shape
        assert len(H_0.shape) == 2 or len(H_0.shape) == 1
        #assert H_0.shape == H_1.shape
        
        [ALPHA, DELTA, GAMMA, OMEGA, CHI, PHI] = createVliegMatrices([alpha, None, None, omega, chi, phi])
        #rotmatrices = [np.asarray(mat) for mat in rotmatrices]
        #[ALPHA, DELTA, GAMMA, OMEGA, CHI, PHI] = rotmatrices
        K = self._ubCalculator.getK()
        ub = np.asarray(self._ubCalculator.getUB())

        Vmat = OMEGA @ (CHI @ (PHI @ ub)) # this can be super expensive!
        Vmat_H_0 = np.squeeze(Vmat @ H_0[...,np.newaxis])
        Vmat_H_1 = np.squeeze(Vmat @ H_1[...,np.newaxis])
        
        # np.swapaxes(ALPHA,1,2) transposes matrices inside the stack, works also with a single matrix
        C_vec = (np.swapaxes(ALPHA,-1,-2)[...,1] * K) + Vmat_H_0
        
        C_vec_sqr = np.sum(C_vec**2, axis=-1)
        Vmat_H_1_sqr = np.sum(Vmat_H_1**2, axis=-1)
        Vmat_H_1_C = np.sum(Vmat_H_1*C_vec, axis=-1)
        
        term1 = Vmat_H_1_C / Vmat_H_1_sqr

        with warnings.catch_warnings(): # will be nan if there is no solution!
            warnings.simplefilter("ignore")
            sqrtTerm = np.sqrt(term1**2 + (K**2 - C_vec_sqr) / Vmat_H_1_sqr)
            
        
        s1 = -(term1 + sqrtTerm)
        s2 = -(term1 - sqrtTerm)
        
        #from IPython import embed; embed()
        
        HKL_1 = H_1*s1[...,np.newaxis] + H_0
        HKL_2 = H_1*s2[...,np.newaxis] + H_0
        
        if keyargs.get('Qalpha',False):
            Qalpha1 = np.squeeze(Vmat @ HKL_1[...,np.newaxis])
            Qalpha2 = np.squeeze(Vmat @ HKL_2[...,np.newaxis])
            return np.concatenate((HKL_1,Qalpha1,s1[...,np.newaxis]),axis=-1), np.concatenate((HKL_2,Qalpha2,s2[...,np.newaxis]),axis=-1)
        else:
            return HKL_1, HKL_2
        
    def anglesIntersectLineEwald(self, H_0, H_1, alpha,omega,phi=0.,chi=0.,**keyargs):
        HKL_Q1, HKL_Q2 = self.intersectLineEwald(H_0, H_1, alpha,omega,phi,chi, Qalpha=True)
        K = self._ubCalculator.getK()
        HKL_1 = HKL_Q1[...,:3]
        HKL_2 = HKL_Q2[...,:3]
        
        Qa_1 = HKL_Q1[...,3:-1] / K
        Qa_2 = HKL_Q2[...,3:-1] / K
        
        sinalpha = np.sin(alpha)
        gamma_1 = np.arcsin(Qa_1[...,2] - sinalpha)
        gamma_2 = np.arcsin(Qa_2[...,2] - sinalpha)
        
        delta_1 = np.arctan2(Qa_1[...,0], Qa_1[...,1] + np.cos(alpha))
        delta_2 = np.arctan2(Qa_2[...,0], Qa_2[...,1] + np.cos(alpha))
        
        #gamma_1 = np.arctan2((Qa_1[...,2] - sinalpha) * np.sin(delta_1), Qa_1[...,0])
        #gamma_2 = np.arctan2((Qa_2[...,2] - sinalpha) * np.sin(delta_2), Qa_2[...,0])
        if keyargs.get('Qalpha',False):
            return np.concatenate((HKL_1,delta_1[...,np.newaxis], gamma_1[...,np.newaxis],HKL_Q1[...,-1][...,np.newaxis]),axis=-1),\
                   np.concatenate((HKL_2,delta_2[...,np.newaxis], gamma_2[...,np.newaxis],HKL_Q2[...,-1][...,np.newaxis]),axis=-1),\
                   Qa_1*K, Qa_2*K
        else:
            return np.concatenate((HKL_1,delta_1[...,np.newaxis], gamma_1[...,np.newaxis],HKL_Q1[...,-1][...,np.newaxis]),axis=-1),\
                   np.concatenate((HKL_2,delta_2[...,np.newaxis], gamma_2[...,np.newaxis],HKL_Q2[...,-1][...,np.newaxis]),axis=-1)
        
    def hkIntersect(self,rod,alpha,omega,phi=0.,chi=0.):
        hk = np.asarray(rod)
        if len(hk.shape) > 1:
            H_0 = np.concatenate((hk,np.zeros((hk.shape[0],1))),axis=-1)
        else:
            H_0 = np.array([*hk,0.])
        H_1 = np.array([0.,0.,1.])
        
        hkl_del_gam_1, hkl_del_gam_2 = self.anglesIntersectLineEwald(H_0, H_1, alpha,omega,phi,chi)
        
        L1 = hkl_del_gam_1[...,2]
        L2 = hkl_del_gam_2[...,2]
        delta1 = hkl_del_gam_1[...,3]
        delta2 = hkl_del_gam_2[...,3]
        gam1 = hkl_del_gam_1[...,4]
        gam2 = hkl_del_gam_2[...,4]
        
        return (L1,gam1,delta1), (L2,gam2,delta2)

    def coordinatesAlpha(self, xyz_rel, omega,phi=0.,chi=0.):
        """Calculates the coordinates xyz_rel in the alpha frame (tilted sample).
        
        :param xyz_rel: real space coordinates in lattice units
        :type xyz_rel: ndarray shape (n,3)
        
        """
        xyz_rel = np.atleast_2d(np.asarray(xyz_rel))
        shape = xyz_rel.shape
    
        OMEGA = calcOMEGA(omega)
        PHI = calcPHI(phi)
        CHI = calcCHI(chi)
        U = self._ubCalculator.getU()
        R = self._ubCalculator.lattice.RealspaceMatrix
        xyz_alpha = (OMEGA @ CHI @ PHI @ U @ R) @ xyz_rel.T
        return (xyz_alpha.T).reshape(shape)
        
    def coordinatesLab(self, xyz_rel,alpha, omega,phi=0.,chi=0.):
        """Calculates the coordinates xyz_rel in the laboratory frame (y along beam).
        
        :param xyz_rel: real space coordinates in lattice units
        :type xyz_rel: ndarray shape (n,3)
        
        """
        xyz_rel = np.atleast_2d(np.asarray(xyz_rel))
        shape = xyz_rel.shape
        
        ALPHA = calcALPHA(alpha)
        xyz_alpha = self.coordinatesAlpha(xyz_rel, omega, phi, chi)
        return ((ALPHA @ xyz_alpha.T).T).reshape(shape)
        
    def anglesOrientationAlpha(self, xyz_rel, xyz_direction):
        xyz_rel = np.atleast_2d(np.asarray(xyz_rel))
        xyz_direction = np.atleast_2d(np.asarray(xyz_direction))
        shape = xyz_rel.shape
        assert xyz_rel.shape == xyz_direction.shape
        
        xyz_rel = xyz_rel.reshape((-1,3))
        xyz_direction = xyz_direction.reshape((-1,3))
        
        om = np.empty(xyz_rel.shape[0])
        chi = np.empty(xyz_rel.shape[0])
        phi = np.empty(xyz_rel.shape[0])
        
        U = np.array(self._ubCalculator.getU())
        R = np.array(self._ubCalculator.lattice.RealspaceMatrix)
        UR = U @ R
        
        for i, (xyz_r, xyz_d) in enumerate(zip(xyz_rel, xyz_direction)):
            urx = UR @ xyz_r
            rotmat = util.rotation_matrix_from_vectors(urx, xyz_d)
            rot = Rotation.from_matrix(rotmat)
            rotangles = rot.as_euler('xyz')
            phi[i] = rotangles[0]
            chi[i] = rotangles[1] 
            om[i] = -rotangles[2] 
        
        return om.reshape(shape[:-1]), chi.reshape(shape[:-1]), phi.reshape(shape[:-1])
        
        
        
        
    def getGeometryCorrection(self):
        return GeometryCorrection(self)

# only for phi/omega scans, partially zmode
# deprecated! can be removed in release!
class GeometryCorrection():
    def __init__(self,vliegangles):
        self._angles = vliegangles
    
    def lorentzFactor(self,delta,beta_in,gamma):
        return 1./(np.sin(delta)*np.cos(beta_in)*np.cos(gamma))
    
    def polarization(self,delta,gamma,alpha,fraction_horiz=1.):
        P_hor = 1. - (np.sin(alpha)*np.cos(delta)*np.cos(gamma) + np.cos(alpha)*np.sin(gamma))**2
        P_vert = 1. - (np.sin(delta)**2)*(np.cos(gamma)**2) 
        return fraction_horiz*P_hor + (1.-fraction_horiz)*P_vert
    
    # without footprint correction
    def activeSurfaceArea(self,delta,alpha,beta_in):
        return 1./(np.sin(delta)*np.cos(alpha-beta_in))
    
    def correctionZmode(self,hkl,fixedangle,fixed='in',polarization_horiz=1.):
        alpha, delta, gamma, omega, chi, phi = self._angles.anglesZmode(hkl,fixedangle,fixed)
        P = self.polarization(delta,gamma,alpha,polarization_horiz)
        #L_phi = self.lorentzFactor(delta,alpha,gamma)
        Carea = self.activeSurfaceArea(delta,alpha,alpha)
        return P*Carea
    
    def correctDatasetZmode(self,hkl,I,fixedangle,fixed='in',polarization_horiz=1.):
        corr = np.empty_like(I)
        for i in range(I.size):
            corr[i] = self.correctionZmode(hkl[i],fixedangle,fixed,polarization_horiz)
        corr /= np.mean(corr)
        return I/corr
    
    def correctionFactorZmode(self,alpha,delta,gamma,polarization_horiz=1.):
        delta = np.abs(delta)
        corr = np.empty_like(delta)
        for i in range(delta.shape[0]):
            P = self.polarization(delta[i],gamma[i],alpha,polarization_horiz)
            Carea = self.activeSurfaceArea(delta[i],alpha,alpha)
            corr[i] = (P*Carea)
        return corr
    
    def correctImageZmode(self,intensity,alpha,delta,gamma,polarization_horiz=1.):
        I = np.copy(intensity)
        delta = np.abs(delta)
        
        for i in range(delta.shape[0]):
            #print(i)
            P = self.polarization(delta[i],gamma[i],alpha,polarization_horiz)
            Carea = self.activeSurfaceArea(delta[i],alpha,alpha)
            #print(I[i].shape)
            #print((P*Carea).shape)
            I[i] = I[i]/(P*Carea)
        return I
    
    def applyImageZmode(self,intensity,alpha,delta,gamma,polarization_horiz=1.):
        I = np.copy(intensity)
        delta = np.abs(delta)
        
        for i in range(delta.shape[0]):
            #print(i)
            P = self.polarization(delta[i],gamma[i],alpha,polarization_horiz)
            Carea = self.activeSurfaceArea(delta[i],alpha,alpha)
            #print(I[i].shape)
            #print((P*Carea).shape)
            I[i] = I[i]*(P*Carea)
        return I
        
        

#pos = [ALPHA, DELTA, GAMMA, OMEGA, CHI, PHI] (angles)

if __name__ == "__main__":

    pt111 = Lattice([ 2.7748 , 2.7748 , 6.7969],[  90. ,  90. , 120.])
    """
    basis = [[1.,0.,0.,0.],[1.,2./3.,1./3.,1./3.],[1.,1./3.,2./3.,2./3.]]
    h = 1.
    k = -1
    for l in range(10):
        print("[%s , %s, %s]: F = %s" % (h,k,l,np.absolute(pt111.F_hkl([h,k,l],basis))))
    """
    ub = UBCalculator(pt111,69.971)
    ub.defaultU()
    
    #ub.setPrimaryReflection(np.deg2rad([0.15,4.2168,1.354,0,0.,0.]),[1.,0.,1.])
    #ub.setSecondayReflection(np.deg2rad([0.15,4.211,2.839,+32.43 + 28.12,0.,0.]),[0.,1.,2.])
    #ub.calculateU()
    
    
    angles = VliegAngles(ub)
    
    
    
    #delta = np.linspace(-np.pi/8,np.pi/8,1100)
    #gamma = np.linspace(-np.pi/8,np.pi/8,1600)
    h , k , l = angles.anglesToHklDetector(0.15,delta,gamma,0.1,0,0)
    


