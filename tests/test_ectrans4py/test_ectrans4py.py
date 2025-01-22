from unittest import main, TestCase
import numpy
from . import data
import ectrans4py

ectrans4py.init_env()

KNUMMAXRESOL = 10


class TestLAM(TestCase):

    gpdims = {'X':54,
              'Y':48,
              'X_CIzone':43,
              'Y_CIzone':37,
              'X_resolution':1300.0,
              'Y_resolution':1300.0}
    truncation = {'in_X':26,
                  'in_Y':23}
    spectra_data_sizes = (2592, 1968)

    def test_etrans_inq(self):
        spectra_data_sizes = ectrans4py.etrans_inq4py(
            self.gpdims['X'],
            self.gpdims['Y'],
            self.gpdims['X_CIzone'],
            self.gpdims['Y_CIzone'],
            self.truncation['in_X'],
            self.truncation['in_Y'],
            KNUMMAXRESOL,
            self.gpdims['X_resolution'],
            self.gpdims['Y_resolution'])
        self.assertEqual(spectra_data_sizes, self.spectra_data_sizes)


class TestGlobal(TestCase):

    gpdims = {'lat_number':150,
              'lon_number_by_lat':data.lon_number_by_lat}
    truncation = {'max':148}
    spectral_data_sizes = (
            33052,
            11175,
            data.zonal_wavenumbers)

    def test_trans_inq4py(self):
        spectral_data_sizes = ectrans4py.trans_inq4py(
            self.gpdims['lat_number'],
            self.truncation['max'],
            len(self.gpdims['lon_number_by_lat']),
            numpy.array(self.gpdims['lon_number_by_lat']),
            KNUMMAXRESOL)
        self.assertEqual(spectral_data_sizes[0:2], self.spectral_data_sizes[0:2])  # dimensions
        numpy.testing.assert_array_equal(spectral_data_sizes[2], self.spectral_data_sizes[2])  # zonal_wavenumbers
