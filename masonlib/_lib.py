'''
	@author: Marius N. Sletten
	
	masonlib._lib
'''

''' Import a minimum since it's a library, also make sure of same functionality atleast on  python ~2.7.10 and ~3.4 '''
from __future__ import print_function, division
from numpy import tan, sin, cos, pi, arange, sqrt, isscalar,\
                    zeros, float64, complex128, ndarray, asarray
from scipy import fftpack
from ._material import struct

def layer(w,m, Lr, name=''):
    '''
        Calculate properties for mechanical layer  defined by the material m
        
        returns
        
        Z = mechanical impedance
        H23vv = mech to mech transferer function
    '''
    
    out = struct(name)
    out.Za, out.Zb = m.impedance(w)
    
    ''' Impedance - Is Lr and previous layer or impedance'''
    if not isscalar(Lr):
        Zr = Lr.Z
        H23vv = Lr.H23vv
    else: 
        Zr = Lr
        H23vv = 1.0
        
    Z0 = m.Z0
    Z2 = Zr/Z0*tan(m.kxl(w))
    out.Z = (1 + 1j * Z0)/(1 + 1j * Z2) * Zr
    out.H23vv = out.Zb/(out.Za + out.Zb + Zr)*H23vv
    
    return out

def transmitter(w, m,  L1, L2, name='transmitter'):
    ''' 
        m is a piezo material.
        Put everything you want on the out-struct 
        
        w = wave number array
        m = piezo material
        L1 = backing layers
        L2 = front layers
        
        returns:
        
        Zt = Electrical impedance
        H23Vv = Volt to particle velocity
        
    '''
    
    
    
    Za, Zb = m.impedance(w)
    
    Za1 = Za + L1.Z
    Za2 = Za + L2.Z    
    
    ZC0 = 1/(1j*w*m.C_0)
    
    Zar = 1/(1/Za1 + 1/Za2);
    Zbar = (Zb + Zar)/m.phi2
    
    ''' Transferer function on the front '''
    
    Zcb = Zb + 1j*m.phi2/(w*m.C_0)
    Zm2 = (1 + Zcb/Zar)*Za2
    
    ''' Output '''
    
    out = struct(name)
    out.Zt = (ZC0*(1 - ZC0/Zbar))
    out.H23Vv = m.phi/Zm2*L2.H23vv
    
    return out

def receiver(w, m, Zr1, ZL, name='receiver'):
    '''
        w = wavenumber
        m = piezo material
        Zr1 = front layers impedance
        ZL = electrical output load
        
        returns 
        Zml = mechanical input impedance
        H45FV = mech to voltage transferer function
    '''
    
    out = struct(name)
    
    Za, Zb = m.impedance(w)
    Za1 = Za + Zr1
    ZC0 = 1/(1j*w*m.C_0)
    
    Zcb = Zb - m.phi2 * ZC0
    Zcl = m.phi2 * (1/ZC0 + 1/ZL)**-1
    Zcbl = Zcb + Zcl
    
    out.Zml = Za + Za1*Zcbl/(Za1 +Zcbl)
    out.H45FV = Za1*Zcl/(m.phi*(Za+Zcbl)*Za1 + Za*Zcbl)
    
    return out

def sourcesensitivity(w, source, c, z = 1.0, name='Sv'):
        
    
    out = struct(name)
    k = w/c
    ''' Frequency domain diffration correction '''
    
    out.Hdiff = source.H23Vv*diffcorr(k, source.r, z, N=256)
    return out

def diffcorr(k,a,z,N=256):
    '''
    k, wavenumber 
    a, aperture radi
    z, distance
    
    k or z is a vector, not both
    '''

    ''' Vectors and constants '''
    theta = pi*arange(0,N)/(2*N);
    kz = k*z;
    
    ''' Simplifications and lambda functions '''
    sint2 = sin(theta)**2;
    
    if isscalar(z):
        Lt = sqrt(z**2 + (2*a*cos(theta))**2);
    else:
        Lt = zeros((z.shape[0],N))
        for i in range(z.shape[0]):
            Lt[i] = sqrt(z[i]**2 + (2*a*cos(theta))**2);
   
    trap = lambda X: 2/N*(sum(X[0:]) + 0.5*(X[0] + X[-1]));
    Xc = lambda X: cos(X)*sint2;
    Xd = lambda X: sin(X)*sint2;
    
    Hd = zeros(kz.shape[0], dtype = complex128);
    if not isscalar(k):
    
        for ii in range(k.shape[0]):
            ''' Integration variable '''
            rho = k[ii]*Lt;
            C = trap(Xc(rho));
            D = trap(Xd(rho));
    
            ''' Real and complex parts '''
            A = C * cos(kz[ii]) + D * sin(kz[ii]);
            B = C * sin(kz[ii]) - D * cos(kz[ii]);
    
            Hd[ii] = 1 - A - 1j*B;
    else:
        ''' Integration variable '''
        
        for ii in range(kz.shape[0]):
            rho = k*Lt[ii]
            C = trap(Xc(rho))
            D = trap(Xd(rho))
            ''' Real and complex parts '''
            A = C * cos(kz[ii]) + D * sin(kz[ii]);
            B = C * sin(kz[ii]) - D * cos(kz[ii]);
    
            Hd[ii] = 1 - A - 1j*B;    
    
    return Hd

def sos_water(T, Pg=1, components=False):
    
    '''
        An approximation of sound speed in fresh water. 
        0.05% accuracy for 0<T<100 [C] and 1<Pg<2 [bar] 
    
    '''
    
    t = T/100.
    
    a0 = 1402.7
    a1 = 488.
    a2 = -482.
    a3 = 153.

    b0 = 15.9
    b1 = 2.8
    b2 = 2.4
    
    if isinstance(T, ndarray):
        p = zeros(T.shape[0],dtype=float64)
        c = zeros(T.shape[0], dtype=float64)
        
        p[:] = (b0 + b1*t + b2*t**2)*Pg/100
        c[:] = a0 + a1*t +a2*t**2 +a3*t**3
    else:
        p = (b0 + b1*t + b2*t**2)*Pg/100
        c = a0 + a1*t +a2*t**2 +a3*t**3
        
    if components:
        return  asarray((c, p))
    else:
        return c + p
    
    
