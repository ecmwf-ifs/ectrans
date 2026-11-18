# (C) Copyright 2025- Meteo-France.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.

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
