# -*- coding: utf-8 -*-

"""ID31 beamline code.

Used for this software with permission from Jakub Drnec

"""

import numpy as np

class ID31DiffractLinTilt(object):
    """
    ID31DiffractLinTilt: transforms linear motors into tilt motor (angular)

    Some angular motors are mechanically moved through linear translations.
    
    muoffset has to be set for each experiment since it is dependent on
    the diffractometer alignment!
    
    Adapted version from the ID31 beamline code.
    """
    def __init__(self, *args, **kwargs):
        self.config = {}
        self.config['a'] = 938
        self.config['b'] = 400
        self.config['muoffset'] = -0.06442416994811659

    @property
    def a(self):
        return self.config.get("a", float)

    @property
    def b(self):
        return self.config.get("b", float)
        
    @property
    def muoffset(self):
        return self.config.get("muoffset", float)

    def c(self, a=None, b=None):
        a = self.a if a is None else a
        b = self.b if b is None else b
        return np.sqrt(np.square(a) + np.square(b))

    def d(self, a=None, b=None):
        a = self.a if a is None else a
        b = self.b if b is None else b
        return np.arctan(b / a)

    def calc_from_real(self, positions_dict):
        a, b = self.a, self.b
        c, d = self.c(a, b), self.d(a, b)
        a2, c2 = np.square(a), np.square(c)
        linear = positions_dict["linear"]
        bc2 = np.square(b + linear)
        tilt = np.arccos((a2 + c2 - bc2) / (2 * a * c)) - d
        return dict(tilt=np.rad2deg(tilt))

    def calc_to_real(self, positions_dict):
        a, b = self.a, self.b
        c, d = self.c(a, b), self.d(a, b)
        a2, c2 = np.square(a), np.square(c)
        tilt = np.deg2rad(positions_dict["tilt"])
        bc = np.sqrt(a2 + c2 - 2 * a * c * np.cos(tilt + d))
        linear = bc - b
        return dict(linear=linear)
        
    def linai_to_mu(self,linai):
        positions = {'linear' : linai}
        tilt = self.calc_from_real(positions)
        return tilt["tilt"] - self.muoffset
        
    def mu_to_linai(self,mu):
        ai = mu + self.muoffset #np.linspace(-1.0,1,11) + muoffset
        tilts = {'tilt' : ai}
        pos = self.calc_to_real(tilts)
        return pos['linear']
