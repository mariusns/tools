'''
	@author: Marius N. Sletten
	masonlib.__init__
	
'''
''' Import a minimum since this a library, also make sure of same functionality atleast on  python ~2.7.10 and ~3.4 '''

from ._material import struct, mechanic, piezoelement, stiffness, speedofsound 
from ._lib import layer, transmitter, receiver, sourcesensitivity, sos_water
