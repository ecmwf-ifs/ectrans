from unittest import main, TestCase
import numpy
from . import data
import ectrans4py
import platform

system = platform.system()
if system == "Linux":
    ectrans4py.init_env(unlimited_stack=True)
elif system == "Darwin":
    ectrans4py.init_env(unlimited_stack=False)
else:
    raise NotImplementedError("ectrans4py does not support Windows")

KNUMMAXRESOL = 10
EPSILON = 1e-10


class ArraysAlmostEqual(object):

    def assert_arrays_diff_under_epsilon(self, x, y):
        diff = x - y
        diffmax = abs(diff.max())
        diffmin = abs(diff.min())
        self.assertTrue(diffmax < EPSILON, "diffmax is {}".format(diffmax))
        self.assertTrue(diffmin < EPSILON, "diffmin is {}".format(diffmin))


class TestLAM(TestCase, ArraysAlmostEqual):

    gpdims = {'X':54,
              'Y':48,
              'X_CIzone':43,
              'Y_CIzone':37,
              'X_resolution':1300.0,
              'Y_resolution':1300.0}
    truncation = {'in_X':26,
                  'in_Y':23}
    spectral_data_sizes = (2592, 1968)
    spdata = data.antwrp1300['sp']
    gpdata = data.antwrp1300['sp2gp']

    def test_etrans_inq(self):
        spectral_data_sizes = ectrans4py.etrans_inq4py(
            self.gpdims['X'],
            self.gpdims['Y'],
            self.gpdims['X_CIzone'],
            self.gpdims['Y_CIzone'],
            self.truncation['in_X'],
            self.truncation['in_Y'],
            KNUMMAXRESOL,
            self.gpdims['X_resolution'],
            self.gpdims['Y_resolution'])
        self.assertEqual(spectral_data_sizes, self.spectral_data_sizes)

    def test_sp2gp(self):
        gpdata = ectrans4py.sp2gp_lam4py(
            self.gpdims['X'],
            self.gpdims['Y'],
            self.gpdims['X_CIzone'],
            self.gpdims['Y_CIzone'],
            self.truncation['in_X'],
            self.truncation['in_Y'],
            KNUMMAXRESOL,
            len(self.spdata.flatten()),
            False,  # no derivatives
            False, # spectral_coeff_order != 'model',
            self.gpdims['X_resolution'],
            self.gpdims['Y_resolution'],
            self.spdata.flatten())[0]
        self.assert_arrays_diff_under_epsilon(gpdata, gpdata.flatten())
    
    def test_gp2sp(self):
        spdata = ectrans4py.gp2sp_lam4py(
            self.spectral_data_sizes[1],
            self.gpdims['X'],
            self.gpdims['Y'],
            self.gpdims['X_CIzone'],
            self.gpdims['Y_CIzone'],
            self.truncation['in_X'],
            self.truncation['in_Y'],
            KNUMMAXRESOL,
            self.gpdims['X_resolution'],
            self.gpdims['Y_resolution'],
            False,  # spectral_coeff_order != 'model',
            self.gpdata.flatten())
        self.assert_arrays_diff_under_epsilon(spdata, spdata.flatten())

class TestGlobal(TestCase, ArraysAlmostEqual):

    gpdims = {'lat_number':150,
              'lon_number_by_lat':data.lon_number_by_lat}
    truncation = {'max':148}
    spectral_data_sizes = (
            33052,
            11175,
            data.zonal_wavenumbers)
    spdata = data.tl149_c24['sp']
    gpdata_latlon = data.tl149_c24['sp2gp']

    # Pack latlon gridded data to reduced grid
    gpdata = numpy.zeros((sum(gpdims['lon_number_by_lat'])))
    offset = 0
    for i in range(gpdims['lat_number']):
        nlon = gpdims['lon_number_by_lat'][i]
        gpdata[offset:offset+nlon] = gpdata_latlon[i,:nlon]
        offset += nlon

    def test_get_legendre_assets(self):
        nspec = sum([self.truncation['max'] + 2 - im for im in range(self.truncation['max']+1)])
        knmeng, weights, polys = ectrans4py.get_legendre_assets(
            self.gpdims['lat_number'],
            self.truncation['max'],
            len(self.gpdims['lon_number_by_lat']),
            nspec,
            self.gpdims['lon_number_by_lat'],
            KNUMMAXRESOL
        )
        weights_sum = sum(weights)
        # The sum of the Gaussian weights should be equal to 1.0
        self.assertTrue(abs(weights_sum - 1.0) < EPSILON, f"sum of weights is {weights_sum}")

    def test_trans_inq4py(self):
        spectral_data_sizes = ectrans4py.trans_inq4py(
            self.gpdims['lat_number'],
            self.truncation['max'],
            len(self.gpdims['lon_number_by_lat']),
            self.gpdims['lon_number_by_lat'],
            KNUMMAXRESOL)
        self.assertEqual(spectral_data_sizes[0:2], self.spectral_data_sizes[0:2])  # dimensions
        numpy.testing.assert_array_equal(spectral_data_sizes[2], self.spectral_data_sizes[2])  # zonal_wavenumbers
    
    def test_sp2gp(self):
        gpdata = ectrans4py.sp2gp_gauss4py(
            self.gpdims['lat_number'],
            self.truncation['max'],
            KNUMMAXRESOL,
            sum(self.gpdims['lon_number_by_lat']),
            len(self.gpdims['lon_number_by_lat']),
            self.gpdims['lon_number_by_lat'],
            len(self.spdata),
            False,  # no derivatives
            False, # spectral_coeff_order != 'model',
            self.spdata)[0]
        self.assert_arrays_diff_under_epsilon(self.gpdata, gpdata)
    
    def test_gp2sp(self):
        spdata = ectrans4py.gp2sp_gauss4py(
            self.spectral_data_sizes[1] * 2,  # *2 for complex coefficients
            self.gpdims['lat_number'],
            self.truncation['max'],
            KNUMMAXRESOL,
            len(self.gpdims['lon_number_by_lat']),
            self.gpdims['lon_number_by_lat'],
            len(self.gpdata),
            False,  # spectral_coeff_order != 'model',
            self.gpdata)
        self.assert_arrays_diff_under_epsilon(self.spdata, spdata)

