'''

    @author: Marius N. Sletten
    
    masonlib.stack
'''

''' Import a minimum since it's a library, also make sure of same functionality atleast on  python ~2.7.10 and ~3.4 '''
from __future__ import print_function, division
from numpy import *
import . as ml
from matplotlib import pyplot as pp

def stack_materials(a):
    
    ''' 
        First  define the mechanical building blocks including:
        piezoelectric material
        front and backing materials
        radiation medias
    '''
    
    piezo = ml.piezoelement('Unknown_material')
    piezo.r = a
    piezo.l = 0
    piezo.rho = 0
    piezo.c33D = 0
    piezo.Qc33 = 0
    piezo.kt = 0
    piezo.eps_33S
    piezo.tandeps = 0
    piezo.calc()
    
    ''' Front and backing layers '''
    
    mat1 = ml.mechanic('mat1')
    mat1.r = a
    mat1.l = 0
    mat1.rho = 0
    mat1.c = 0
    mat1.Qc33 = 0
    mat1.calc()
    
    mat2 = ml.mechanic('mat2')
    mat2.r = a
    mat2.l = 0
    mat2.rho = 0
    mat2.c = 0
    mat2.Qc33 = 0
    mat2.calc()
    
    m = ml.struct('materials')
    m.piezo = piezo
    m.mat1 = mat1
    m.mat2 = mat2
    
    
    return m
    
def stack_layout(w, T = 0):
    ''' 
        Together with the stack layout  a frequency interval 
        has to be defined  in kHz unit
        
        also radiation medium has to be defined
        
    '''
    a = 0
    A = pi*a**2
    
    ''' Radiation into water at T degree C '''
    rho_m = 0
    c_m = ml.sos_water(T)
    Z_m = c_m*rho_m*A
    
    
    mat = stack_materials(a=a)
    
    ''' Front ''' 

    Lf = ml.layer(w, mat.mat1, Z_m)
    
    ''' Back '''
    
    Lb = ml.layer(w, mat.mat2, Z_m)
   
    ''' Source transducer '''
    source = ml.transmsitter(w, mat.piezo, Lb, Lf)

    return source

    
def main():
    
    ''' Run parameters '''
    scale = 0 
    start = 0
    stop = 0
    inc = 0
    T = 0
    c_m = ml.sos_water(T)
    
    ''' Frequency vector '''
    f = arange(start, stop, inc, dtype=float64)
    w = 2*pi*f*scale
    
    
    ''' Calculate source transducer resonse '''
    source = stack_layout(w, T)
    
    ''' Have a look at source sensitivity '''
    Sv = ml.sourcesensitivity(w, source, c_m)
    
    fig, ax = pp.subplots(nrows=2, sharex=True)
    fig.suptitle('Source sensitivity modul and phase')
    ax[0].plot(f, abs(Sv)*scale)
    ax[1].plot(f, arctan2(Sv.imag, Sv.real))
    ax[0].grid(True)
    ax[1].grid(True)
    ax[1].set_xlabel('[kHz]')
    pp.show()
    
    
    
    
if __name__ =='__main__':
    main()
