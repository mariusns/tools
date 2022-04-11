'''

	@author: Marius N. Sletten
	masonlib._material
'''
''' Import a minimum since it's a library, also make sure of same functionality atleast on  python ~2.7.10 and ~3.4 '''

from __future__ import print_function, division
from numpy import pi, sin, tan,sqrt, asarray, zeros, float64, complex128

''' Basic data structure '''
class struct (object):
    
    name = ''
    def __init__(self, name=''):
        self.name = name
    def fields(self):
        return self.__dict__.keys()
    def __str__(self):
        s = ''
        for f in self.fields():        
            s += f + ': ' + str(getattr(self, f)) + '\n'
        return s
    
    def view(self):
        print(self)

''' Mechanical material object, struct super class '''
class mechanic(struct):
    
    r = 0
    l = 0
    c = 0
    rho = 0
    c33D = 0    
    Qc33 = 0
    
    def __init__(self,name=''):
        self.name = name

    
    def calc(self):
        ''' Calculate internal relations '''
        
        if self.c33D !=0:
            # Calc soundspeed using c33d 
            try:
                if self.Qc33 != 0 and self.c33D.imag == 0:
                    self.c33D = self.c33D*(1 + 1j/self.Qc33)
                self.c = sqrt(self.c33D/self.rho)
            except ZeroDivisionError:
                pass        
        else:  
            try:
                self.c = self.c*(1 + 1j/(2*self.Qc33))
            except ZeroDivisionError:
                pass


        self.A = pi*self.r**2;
        self.Z0s = self.rho*self.c;
        self.Z0 = self.Z0s*self.A;
        
        try:
            self.Al = self.A/self.l
        except ZeroDivisionError:
            self.Al = 0
        
    def kxl(self,w):
        ''' The complex wavenumber length for the element '''
        return w/self.c*self.l

    def impedance(self,w, el=False):
        """Charateristic layer impedance, mechanical"""

        kxl = self.kxl(w)
        Za = 1j*self.Z0x*tan(kxl/2) 
        Zb = self.Z0x/(1j*sin(kxl))
        return Za, Zb


''' Piezoelectric material, mechanic super class '''
class piezoelement(mechanic):
    ''' Use Q_epsS or tand for dielectric loss '''
    Q_epsS = 0
    tandeps = 0
    eps_33S_rel = 0
    ''' Use kt if e33 is 0 '''
    e_33 = 0
    kt = 0
    eps_0 = 8.8541878176e-12
    C_0 = 0
    c = 0
    def __init__(self, name=''):
        self.name = name
    
    def calc(self):
        ''' The first part is calculated by it's  super ??'''
        
        if self.c33D !=0:
            # Calc soundspeed using c33d 
            try:
                if self.Qc33 != 0 and self.c33D.imag == 0:
                    self.c33D = self.c33D*(1 + 1j/self.Qc33)
                    
                self.c = sqrt(self.c33D/self.rho)
                
            except ZeroDivisionError:
                pass
        else:  
            try:
                self.c = self.c*(1 + 1j/(2*self.Qc33))
            except ZeroDivisionError:
                pass
               
        self.A = pi*self.r**2;
        self.Z0s = self.rho*self.c;
        self.Z0 = self.Z0s*self.A;
        try:
            self.Al = self.A/self.l
        except ZeroDivisionError:
            self.Al = 0
        
        self.eps_33S = self.eps_33S_rel*self.eps_0*(1 - 1j*self.tandeps)
        ''' Set e33 or kt '''
        if self.kt > 0:
            self.e_33 = self.kt * sqrt(self.eps_33S*self.c33D)
        else:
            self.kt = self.e_33 / sqrt(self.eps_33S*self.c33D)
            
        self.phi = self.e_33*self.Al
        self.phi2 = self.phi**2
        self.C_0 = self.eps_33S*self.Al

''' Material constant converters '''


''' Helper functions for conversion using Young modulus and Poisson's'''

def stiffness(Y, q, Q=(0,0,0)):
    ''' 
        Youngs / Poisson's to c_D for non-piezo material 
    
    '''
    cD = zeros(3, dtype=complex128)
    
    cD[0] = Y/(1-q)*(1 + q/(1-2*q))    
    cD[1] = q*Y/((1+q)*(1-2*q))
    cD[2] = Y/(2*(1 + q))
    
    for i in range(3):
        try:
            cD[i] = cD[i]*(1 +1j/Q[i])
        except ZeroDivisionError:
            pass

    return cD

def speedofsound(Y,q,rho, Q=(0,0,0)):
    ''' 
        This gives the compressional and shear speed 
    '''
    
    cD = stiffness(Y, q, Q = Q)
    c = zeros(2, dtype=complex128)
    c[0] = sqrt(cD[0]/rho)
    c[1] = sqrt(cD[2]/rho)
        
    return c
