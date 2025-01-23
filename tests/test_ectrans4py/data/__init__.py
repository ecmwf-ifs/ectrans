import numpy
import os
_here = os.path.abspath(os.path.dirname(__file__))

lon_number_by_lat = numpy.load(os.path.join(_here, 'lon_number_by_lat.npy'))
zonal_wavenumbers = numpy.load(os.path.join(_here, 'zonal_wavenumbers.npy'))

antwrp1300 = {
    'sp'       : numpy.load(os.path.join(_here, 'antwrp1300-s1t@sp.npy')),
    'sp2gp'    : numpy.load(os.path.join(_here, 'antwrp1300-s1t@sp2gp.npy')),
    }
tl149_c24 = {
    'sp'       : numpy.load(os.path.join(_here, 'tl149-c24-s1t@sp.npy')),
    'sp2gp'    : numpy.load(os.path.join(_here, 'tl149-c24-s1t@sp2gp.npy')),
    }
