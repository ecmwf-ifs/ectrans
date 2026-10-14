# (C) Copyright 2026- ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.

import platform
from types import SimpleNamespace
from unittest import TestCase

import numpy as np
import pytest

import ectrans4py

mpi4py = pytest.importorskip("mpi4py")
from mpi4py import MPI

TRUNCATION = 79
NUM_LATS = 2 * (TRUNCATION + 1)
LONS_PER_LAT = [20 + 4 * i for i in range(NUM_LATS // 2)]
LONS_PER_LAT = np.array(LONS_PER_LAT + LONS_PER_LAT[::-1])
NFLD = 20

# A second resolution, set up next to the first one: another truncation, on a regular Gaussian
# grid. The maximum number of resolutions leaves room for a third one, never set up.
TRUNCATION_2 = 47
NUM_LATS_2 = 72
LONS_PER_LAT_2 = np.full(NUM_LATS_2, 2 * NUM_LATS_2, dtype=np.int64)
KMAX_RESOL = 3

# Earth radius used by ecTrans (the default of SETUP_TRANS0's PRAD)
EARTH_RADIUS = 6371229.0

# Tolerances. Data movement and minima/maxima are exact and compared bit for bit; the rest is
# compared relative to the size of the expected values. EPS: machine epsilon of the build.
EPS = np.finfo(ectrans4py._REAL).eps
# One spherical harmonic, transformed once: rounding only, except that winds -> vorticity
# leaves spurious high-degree coefficients of up to 2e-13 in double precision
ANALYTIC_TOL = max(1.0e-12, 1000 * EPS)
# Spectral -> grid -> spectral is inexact on the reduced grid, which drops the zonal wavenumbers
# above each latitude's cut-off near the poles (up to 7e-10 at T79); rounding dominates in sp
ROUND_TRIP_TOL = max(1.0e-8, 1000 * EPS)
# Adjoint identities: the threshold of ecTrans's own adjoint tests (tests/trans/adjoint)
ADJOINT_TOL = 20000 * EPS

# Shared transform geometry, populated once by setUpModule
context = SimpleNamespace()
# The same for the second resolution
context_2 = SimpleNamespace()


def _approx_eq(a, b, tol=1.0e-4):
    return np.all(np.abs(a - b) <= tol)


def _from_which_rank(rank, nfld):
    return np.ones(nfld, dtype=np.int64)


def _initialise_global_spectral_field(rank, nfld, num_spec2_glob):
    # Only rank 0 will have non-zero, meaningful data in the global spectral field
    spectral_field_glob = np.zeros((nfld, num_spec2_glob), dtype=ectrans4py._REAL)
    if rank == 0:
        spectral_field_glob[:, 0] = 1.0 # (m = 0, n = 0) real part is 1.0
    return spectral_field_glob



def _initialise_local_spectral_field(rank, nfld, num_spec2_glob, num_spec2_loc):
    spectral_field_glob = _initialise_global_spectral_field(rank, nfld, num_spec2_glob)
    return ectrans4py.dist_spec4py(num_spec2_glob, num_spec2_loc, nfld,
                                   _from_which_rank(rank, nfld), spectral_field_glob)



def _inquire_geometry(geom, kresol, truncation, lons_per_lat):
    """Store in geom the sizes of resolution kresol, set up with the given truncation and grid
    (number of points of each latitude row, north to south), and the coordinates of the points
    and coefficients of its global fields."""
    num_lats = len(lons_per_lat)
    (num_grid_points_loc, _, num_spec2_loc, num_grid_points_glob, num_spec2_glob,
     inquired_truncation, _, sin_lats, gauss_weights) = \
        ectrans4py.trans_inq4py(kresol, num_lats, 0, 0, np.array([], dtype=np.int64), 0)
    geom.kresol = kresol
    geom.truncation = truncation
    geom.lons_per_lat = lons_per_lat
    geom.num_grid_points_loc = num_grid_points_loc
    geom.num_spec2_loc = num_spec2_loc
    geom.num_grid_points_glob = num_grid_points_glob
    geom.num_spec2_glob = num_spec2_glob
    geom.inquired_truncation = inquired_truncation
    geom.sin_lats = sin_lats
    geom.gauss_weights = gauss_weights

    # Grid points: latitude rows from north to south, each starting at longitude 0
    row = np.repeat(np.arange(num_lats), lons_per_lat)
    geom.sin_lat = geom.sin_lats[row]
    # (1 - sin)(1 + sin) avoids the cancellation of 1 - sin^2 near the poles
    geom.cos_lat = np.sqrt((1.0 - geom.sin_lat) * (1.0 + geom.sin_lat))
    geom.lon = np.concatenate([2.0 * np.pi * np.arange(nlon) / nlon for nlon in lons_per_lat])

    # Spectral coefficients: m = 0..T, then n = m..T, each a (real, imaginary) pair. Weights for
    # inner products and norms: 2 for m > 0 (a coefficient also stands for its conjugate), 0 for
    # the imaginary part of m = 0 (a real field has none), 1 for the real part of m = 0
    zonal, total = np.array([(m, n) for m in range(truncation + 1)
                             for n in range(m, truncation + 1)]).T
    imaginary = np.tile([False, True], zonal.size)
    geom.zonal_wavenumber = np.repeat(zonal, 2)
    geom.total_wavenumber = np.repeat(total, 2)
    geom.spectral_weights = np.where(geom.zonal_wavenumber > 0, 2.0,
                                     np.where(imaginary, 0.0, 1.0))


def _spectral_index(geom, m, n):
    """Index of the real part of coefficient (m, n) in a global spectral field."""
    return np.flatnonzero((geom.zonal_wavenumber == m) & (geom.total_wavenumber == n))[0]


def _random_field(shape, seed):
    """Uniform random values in [-1, 1) in working precision, the same on every rank."""
    return np.random.default_rng(seed).uniform(-1.0, 1.0, shape).astype(ectrans4py._REAL)


def _random_grid_field(geom, nfld, seed):
    return _random_field((nfld, geom.num_grid_points_glob), seed)


def _random_spectral_field(geom, nfld, seed, zero_mean=False):
    """Random global spectral field of real fields (no imaginary part for m = 0). With
    zero_mean, the (0, 0) coefficient is zero too, as for vorticity and divergence."""
    field = _random_field((nfld, geom.num_spec2_glob), seed)
    field[:, geom.spectral_weights == 0.0] = 0.0
    if zero_mean:
        field[:, geom.total_wavenumber == 0] = 0.0
    return field


def _distribute_spectral(geom, field_glob):
    nfld = field_glob.shape[0]
    return ectrans4py.dist_spec4py(geom.num_spec2_glob, geom.num_spec2_loc, nfld,
                                   _from_which_rank(context.my_rank, nfld), field_glob,
                                   geom.kresol)


def _gather_spectral(geom, field_loc):
    nfld = field_loc.shape[0]
    return ectrans4py.gath_spec4py(geom.num_spec2_glob, geom.num_spec2_loc, nfld,
                                   _from_which_rank(context.my_rank, nfld), field_loc,
                                   geom.kresol)


def _distribute_grid(geom, field_glob):
    nfld = field_glob.shape[0]
    return ectrans4py.dist_grid4py(geom.num_grid_points_glob, geom.num_grid_points_loc, nfld,
                                   _from_which_rank(context.my_rank, nfld), field_glob,
                                   geom.kresol)


def _gather_grid(geom, field_loc):
    nfld = field_loc.shape[0]
    return ectrans4py.gath_grid4py(geom.num_grid_points_glob, geom.num_grid_points_loc, nfld,
                                   _from_which_rank(context.my_rank, nfld), field_loc,
                                   geom.kresol)


def _on_rank_0(compute):
    """compute() on rank 0, which holds the gathered fields, broadcast so that all ranks agree."""
    return MPI.COMM_WORLD.bcast(compute() if context.my_rank == 0 else None, root=0)


def _relative_error(actual, expected):
    """Largest difference relative to the largest expected magnitude."""
    expected = np.asarray(expected, dtype=np.float64)
    return float(np.max(np.abs(np.asarray(actual, dtype=np.float64) - expected))
                 / np.max(np.abs(expected)))


def _grid_dot(*pairs):
    """Grid-point inner product summed over pairs of local fields: unweighted sum over all grid
    points, as SCALPRODGP of tests/trans/adjoint."""
    local = sum(np.sum(a.astype(np.float64) * b.astype(np.float64)) for a, b in pairs)
    return MPI.COMM_WORLD.allreduce(float(local), op=MPI.SUM)


def _spectral_dot(geom, *pairs):
    """Spectral inner product summed over pairs of local fields, weighted as SCALPRODSP of
    tests/trans/adjoint."""
    gathered = [(_gather_spectral(geom, a), _gather_spectral(geom, b)) for a, b in pairs]
    return _on_rank_0(lambda: float(sum(np.sum(geom.spectral_weights * a.astype(np.float64)
                                               * b.astype(np.float64)) for a, b in gathered)))


def _adjoint_error(lhs, rhs):
    """Relative difference between <F x, y> and <x, F* y>, as in tests/trans/adjoint."""
    return abs(lhs - rhs) / abs(lhs)


def _first_task_only(outputs):
    """outputs on the first task, which receives the norms, and nothing on the others (where
    the average of gpnorm_trans4py is undefined)."""
    return outputs if context.my_rank == 0 else ()


def _binding_calls(geom):
    """Every distributed binding that takes KRESOL, as (name, call) pairs: call(*kresol) runs
    the binding on fixed inputs of resolution geom, with *kresol appended to its arguments
    (none: the default resolution), and returns its outputs."""
    nfld = 2
    root = np.ones(nfld, dtype=np.int64)
    spectral_glob = _random_spectral_field(geom, nfld, seed=25, zero_mean=True)
    grid_glob = _random_grid_field(geom, nfld, seed=26)
    spectral_loc = _distribute_spectral(geom, spectral_glob)
    grid_loc = _distribute_grid(geom, grid_glob)
    spec2_glob, spec2_loc = geom.num_spec2_glob, geom.num_spec2_loc
    gp_glob, gp_loc = geom.num_grid_points_glob, geom.num_grid_points_loc
    sizes = (spec2_loc, gp_loc, nfld)
    return [
        ("dist_spec4py", lambda *k: ectrans4py.dist_spec4py(
            spec2_glob, spec2_loc, nfld, root, spectral_glob, *k)),
        ("gath_spec4py", lambda *k: ectrans4py.gath_spec4py(
            spec2_glob, spec2_loc, nfld, root, spectral_loc, *k)),
        ("dist_grid4py", lambda *k: ectrans4py.dist_grid4py(
            gp_glob, gp_loc, nfld, root, grid_glob, *k)),
        ("gath_grid4py", lambda *k: ectrans4py.gath_grid4py(
            gp_glob, gp_loc, nfld, root, grid_loc, *k)),
        ("inv_trans_scalar_dist4py", lambda *k: ectrans4py.inv_trans_scalar_dist4py(
            *sizes, spectral_loc, *k)),
        ("inv_trans_scalar_ders_dist4py", lambda *k: ectrans4py.inv_trans_scalar_ders_dist4py(
            *sizes, spectral_loc, *k)),
        ("dir_trans_scalar_dist4py", lambda *k: ectrans4py.dir_trans_scalar_dist4py(
            *sizes, grid_loc, *k)),
        ("inv_trans_uv_dist4py", lambda *k: ectrans4py.inv_trans_uv_dist4py(
            *sizes, spectral_loc, spectral_loc, *k)),
        ("dir_trans_uv_dist4py", lambda *k: ectrans4py.dir_trans_uv_dist4py(
            *sizes, grid_loc, grid_loc, *k)),
        ("inv_trans_scalar_dist4py_ad", lambda *k: ectrans4py.inv_trans_scalar_dist4py_ad(
            *sizes, grid_loc, *k)),
        ("inv_trans_scalar_ders_dist4py_ad",
         lambda *k: ectrans4py.inv_trans_scalar_ders_dist4py_ad(
             *sizes, grid_loc, grid_loc, grid_loc, *k)),
        ("dir_trans_scalar_dist4py_ad", lambda *k: ectrans4py.dir_trans_scalar_dist4py_ad(
            *sizes, spectral_loc, *k)),
        ("inv_trans_uv_dist4py_ad", lambda *k: ectrans4py.inv_trans_uv_dist4py_ad(
            *sizes, grid_loc, grid_loc, *k)),
        ("dir_trans_uv_dist4py_ad", lambda *k: ectrans4py.dir_trans_uv_dist4py_ad(
            *sizes, spectral_loc, spectral_loc, *k)),
        ("specnorm4py", lambda *k: _first_task_only(ectrans4py.specnorm4py(
            spec2_loc, nfld, spectral_loc, *k))),
        ("gpnorm_trans4py", lambda *k: _first_task_only(ectrans4py.gpnorm_trans4py(
            gp_loc, nfld, grid_loc, *k))),
    ]


def _differing_outputs(names, outputs, other_outputs):
    """Names of the bindings whose outputs (an array or a tuple of arrays) differ in any bit on
    any task."""
    def as_tuple(output):
        return output if isinstance(output, tuple) else (output,)

    differing = [name for name, output, other in zip(names, outputs, other_outputs)
                 if not all(np.array_equal(a, b) for a, b in zip(as_tuple(output),
                                                                as_tuple(other)))]
    return sorted(set(sum(MPI.COMM_WORLD.allgather(differing), [])))


def setUpModule():
    system = platform.system()
    if system == "Linux":
        ectrans4py.init_env(unlimited_stack=True)
    elif system == "Darwin":
        ectrans4py.init_env(unlimited_stack=False)
    else:
        raise NotImplementedError("ectrans4py does not support Windows")

    comm = MPI.COMM_WORLD
    n_rank = comm.Get_size()
    my_rank = comm.Get_rank()
    if comm.Get_rank() == 0:
        print(f"Testing with {comm.Get_size()} MPI processes")
    ectrans4py.mpl_init4py()

    # KPRGPNS=n_rank, KPRGPEW=1, KPRTRW=n_rank, LDEQ_REGIONS=False, KMAX_RESOL, LDMPOFF=False
    ectrans4py.setup_trans0_4py(n_rank, 1, n_rank, False, KMAX_RESOL, False)

    # KSMAX=TRUNCATION, KDGL=num_lats, KSLOEN=num_lats, KLOEN=lons_per_lat, LDSPLIT=True
    # LDUSEFLT=False
    kresol = ectrans4py.setup_trans_4py(TRUNCATION, NUM_LATS, NUM_LATS, LONS_PER_LAT, True, False)

    context.my_rank = my_rank
    context.n_rank = n_rank
    _inquire_geometry(context, kresol, TRUNCATION, LONS_PER_LAT)

    # Outputs of the first resolution before the second one is set up
    context.reference_outputs = [call() for _, call in _binding_calls(context)]

    kresol_2 = ectrans4py.setup_trans_4py(TRUNCATION_2, NUM_LATS_2, NUM_LATS_2, LONS_PER_LAT_2,
                                          True, False)
    _inquire_geometry(context_2, kresol_2, TRUNCATION_2, LONS_PER_LAT_2)


def tearDownModule():
    ectrans4py.mpl_end4py()


class DistributedTestCase(TestCase):
    # Transform geometry the tests run on: the first resolution set up by setUpModule (the
    # classes at the end of the module run the same tests on the second one)
    geom = context


class TestInquiry(DistributedTestCase):
    def test_sizes(self):
        geom = self.geom
        # Every grid point and every spectral coefficient is held by exactly one rank
        num_grid_points = MPI.COMM_WORLD.allreduce(geom.num_grid_points_loc, op=MPI.SUM)
        num_spec2 = MPI.COMM_WORLD.allreduce(geom.num_spec2_loc, op=MPI.SUM)

        assert geom.inquired_truncation == geom.truncation
        assert geom.num_grid_points_glob == np.sum(geom.lons_per_lat)
        assert geom.num_spec2_glob == (geom.truncation + 1) * (geom.truncation + 2)
        assert num_grid_points == geom.num_grid_points_glob
        assert num_spec2 == geom.num_spec2_glob

    def test_latitudes(self):
        geom = self.geom
        # Sines of the Gaussian latitudes, north to south, and Gaussian weights summing to one
        # (each weight rounded to working precision)
        assert np.all(np.abs(geom.sin_lats) < 1.0)
        assert np.all(np.diff(geom.sin_lats) < 0.0)
        weights_sum = np.sum(geom.gauss_weights, dtype=np.float64)
        assert abs(weights_sum - 1.0) <= len(geom.lons_per_lat) * EPS

    def test_sin_lats_double_precision(self):
        geom = self.geom
        # The sines of the Gaussian latitudes come in double precision whatever the precision of
        # the build: the Gauss-Legendre nodes, north to south. numpy's nodes agree with them to
        # within an ulp or two; rounded to single precision, they would be off by up to 3e-8
        nodes, _ = np.polynomial.legendre.leggauss(len(geom.lons_per_lat))
        assert geom.sin_lats.dtype == np.float64
        assert np.max(np.abs(geom.sin_lats - nodes[::-1])) <= 100 * np.finfo(np.float64).eps


class TestInverse(DistributedTestCase):
    def test_scalar(self):
        spectral_field_loc = _initialise_local_spectral_field(context.my_rank, NFLD,
                                                              context.num_spec2_glob,
                                                              context.num_spec2_loc)
        gridded_field_loc = ectrans4py.inv_trans_scalar_dist4py(context.num_spec2_loc,
                                                                context.num_grid_points_loc, NFLD,
                                                                spectral_field_loc)
        gridded_field_glob = ectrans4py.gath_grid4py(context.num_grid_points_glob,
                                                     context.num_grid_points_loc, NFLD,
                                                     _from_which_rank(context.my_rank, NFLD),
                                                     gridded_field_loc)

        if context.my_rank == 0:
            range_val = float(np.max(gridded_field_glob) - np.min(gridded_field_glob))
            assert range_val <= 1.0e-4

    def test_scalar_ders(self):
        geom = self.geom
        # Y(n=1, m=0) = sqrt(3) sin(lat) and the real part of Y(n=1, m=1), sqrt(6) cos(lat)
        # cos(lon), with their N-S derivatives (1/a) d/dlat and E-W derivatives
        # 1/(a cos(lat)) d/dlon, a being the Earth's radius
        spectral_field_glob = np.zeros((2, geom.num_spec2_glob), dtype=ectrans4py._REAL)
        spectral_field_glob[0, _spectral_index(geom, 0, 1)] = 1.0
        spectral_field_glob[1, _spectral_index(geom, 1, 1)] = 1.0
        fields_loc = ectrans4py.inv_trans_scalar_ders_dist4py(
            geom.num_spec2_loc, geom.num_grid_points_loc, 2,
            _distribute_spectral(geom, spectral_field_glob), geom.kresol)
        value, ns_der, ew_der = [_gather_grid(geom, field) for field in fields_loc]

        def error():
            sin_lat, cos_lat, lon = geom.sin_lat, geom.cos_lat, geom.lon
            a = EARTH_RADIUS
            return max(
                _relative_error(value, [np.sqrt(3) * sin_lat,
                                        np.sqrt(6) * cos_lat * np.cos(lon)]),
                _relative_error(ns_der, [np.sqrt(3) * cos_lat / a,
                                         -np.sqrt(6) * sin_lat * np.cos(lon) / a]),
                _relative_error(ew_der, [np.zeros_like(lon), -np.sqrt(6) * np.sin(lon) / a]))

        assert _on_rank_0(error) <= ANALYTIC_TOL

    def test_uv(self):
        geom = self.geom
        # Solid-body rotation: vorticity 2 sin(lat) / a, i.e. 2 / (sqrt(3) a) times
        # Y(n=1, m=0), and no divergence give the winds u = cos(lat), v = 0
        vorticity_glob = np.zeros((1, geom.num_spec2_glob), dtype=ectrans4py._REAL)
        vorticity_glob[0, _spectral_index(geom, 0, 1)] = 2.0 / (np.sqrt(3.0) * EARTH_RADIUS)
        divergence_glob = np.zeros((1, geom.num_spec2_glob), dtype=ectrans4py._REAL)
        u_loc, v_loc = ectrans4py.inv_trans_uv_dist4py(
            geom.num_spec2_loc, geom.num_grid_points_loc, 1,
            _distribute_spectral(geom, vorticity_glob),
            _distribute_spectral(geom, divergence_glob), geom.kresol)
        u_glob, v_glob = _gather_grid(geom, u_loc), _gather_grid(geom, v_loc)

        def error():
            return max(_relative_error(u_glob[0], geom.cos_lat), float(np.max(np.abs(v_glob))))

        assert _on_rank_0(error) <= ANALYTIC_TOL


class TestDirect(DistributedTestCase):
    def test_scalar(self):
        gridded_field_glob = np.ones((NFLD, context.num_grid_points_glob), dtype=ectrans4py._REAL)

        gridded_field_loc = ectrans4py.dist_grid4py(context.num_grid_points_glob,
                                                    context.num_grid_points_loc, NFLD,
                                                    _from_which_rank(context.my_rank, NFLD),
                                                    gridded_field_glob)

        spectral_field_loc = ectrans4py.dir_trans_scalar_dist4py(context.num_spec2_loc,
                                                                 context.num_grid_points_loc, NFLD,
                                                                 gridded_field_loc)

        spectral_field_glob = ectrans4py.gath_spec4py(context.num_spec2_glob,
                                                      context.num_spec2_loc, NFLD,
                                                      _from_which_rank(context.my_rank, NFLD),
                                                      spectral_field_loc)

        if context.my_rank == 0:
            assert _approx_eq(spectral_field_glob[:,0], 1.0)
            assert _approx_eq(spectral_field_glob[:,1:], 0.0)

    def test_uv(self):
        geom = self.geom
        # The solid-body rotation of TestInverse.test_uv: u = cos(lat) and v = 0 give the
        # vorticity 2 / (sqrt(3) a) times Y(n=1, m=0) and no divergence
        u_glob = geom.cos_lat[np.newaxis, :].astype(ectrans4py._REAL)
        v_glob = np.zeros_like(u_glob)
        vorticity_loc, divergence_loc = ectrans4py.dir_trans_uv_dist4py(
            geom.num_spec2_loc, geom.num_grid_points_loc, 1,
            _distribute_grid(geom, u_glob), _distribute_grid(geom, v_glob), geom.kresol)
        vorticity_glob = _gather_spectral(geom, vorticity_loc)
        divergence_glob = _gather_spectral(geom, divergence_loc)

        def error():
            expected = np.zeros((1, geom.num_spec2_glob))
            expected[0, _spectral_index(geom, 0, 1)] = 2.0 / (np.sqrt(3.0) * EARTH_RADIUS)
            return max(_relative_error(vorticity_glob, expected),
                       float(np.max(np.abs(divergence_glob)) / np.max(expected)))

        assert _on_rank_0(error) <= ANALYTIC_TOL


class TestRoundTrip(DistributedTestCase):
    def test_scalar(self):
        geom = self.geom
        # spectral -> grid -> spectral
        spectral_field_glob = _random_spectral_field(geom, NFLD, seed=1)
        gridded_field_loc = ectrans4py.inv_trans_scalar_dist4py(
            geom.num_spec2_loc, geom.num_grid_points_loc, NFLD,
            _distribute_spectral(geom, spectral_field_glob), geom.kresol)
        spectral_field_loc = ectrans4py.dir_trans_scalar_dist4py(
            geom.num_spec2_loc, geom.num_grid_points_loc, NFLD, gridded_field_loc, geom.kresol)
        round_trip_glob = _gather_spectral(geom, spectral_field_loc)

        def error():
            return _relative_error(round_trip_glob, spectral_field_glob)

        assert _on_rank_0(error) <= ROUND_TRIP_TOL

    def test_uv(self):
        geom = self.geom
        # vorticity, divergence -> u, v -> vorticity, divergence
        vorticity_glob = _random_spectral_field(geom, NFLD, seed=2, zero_mean=True)
        divergence_glob = _random_spectral_field(geom, NFLD, seed=3, zero_mean=True)
        u_loc, v_loc = ectrans4py.inv_trans_uv_dist4py(
            geom.num_spec2_loc, geom.num_grid_points_loc, NFLD,
            _distribute_spectral(geom, vorticity_glob),
            _distribute_spectral(geom, divergence_glob), geom.kresol)
        vorticity_loc, divergence_loc = ectrans4py.dir_trans_uv_dist4py(
            geom.num_spec2_loc, geom.num_grid_points_loc, NFLD, u_loc, v_loc, geom.kresol)
        round_trip_glob = [_gather_spectral(geom, vorticity_loc),
                           _gather_spectral(geom, divergence_loc)]

        def error():
            return _relative_error(round_trip_glob, [vorticity_glob, divergence_glob])

        assert _on_rank_0(error) <= ROUND_TRIP_TOL

    def test_scalar_ders(self):
        geom = self.geom
        # The E-W and N-S derivatives of a velocity potential chi are the winds u and v of a
        # purely divergent flow: transformed back, they have no vorticity and the divergence
        # -n (n + 1) chi / a^2
        chi_glob = _random_spectral_field(geom, NFLD, seed=4)
        _, ns_der_loc, ew_der_loc = ectrans4py.inv_trans_scalar_ders_dist4py(
            geom.num_spec2_loc, geom.num_grid_points_loc, NFLD,
            _distribute_spectral(geom, chi_glob), geom.kresol)
        vorticity_loc, divergence_loc = ectrans4py.dir_trans_uv_dist4py(
            geom.num_spec2_loc, geom.num_grid_points_loc, NFLD, ew_der_loc, ns_der_loc,
            geom.kresol)
        vorticity_glob = _gather_spectral(geom, vorticity_loc)
        divergence_glob = _gather_spectral(geom, divergence_loc)

        def error():
            n = geom.total_wavenumber
            laplacian = -n * (n + 1) / EARTH_RADIUS ** 2 * chi_glob.astype(np.float64)
            return _relative_error([vorticity_glob, divergence_glob],
                                   [np.zeros_like(laplacian), laplacian])

        assert _on_rank_0(error) <= ROUND_TRIP_TOL


class TestDistributeGather(DistributedTestCase):
    # Field j is distributed from, and gathered back to, MPL rank j mod (number of ranks) + 1.
    # The global array of a rank holds only its own fields, in its first rows (as the ecTrans
    # DIST_* and GATH_* routines expect). Distributing and gathering only move data, so the
    # round trip reproduces the global arrays bit for bit (GATH_SPEC returns no imaginary part
    # for m = 0, so the spectral fields have none).

    def _own_fields(self, field_glob):
        """Rank owning each field, and field_glob with the rows of this rank's fields first and
        the other rows zero."""
        owners = np.arange(NFLD, dtype=np.int64) % context.n_rank + 1
        num_own = np.count_nonzero(owners == context.my_rank + 1)
        own_field_glob = np.zeros_like(field_glob)
        own_field_glob[:num_own] = field_glob[:num_own]
        return owners, own_field_glob

    def test_dist_spec(self):
        geom = self.geom
        owners, spectral_field_glob = self._own_fields(
            _random_spectral_field(geom, NFLD, seed=5 + context.my_rank))
        spectral_field_loc = ectrans4py.dist_spec4py(geom.num_spec2_glob, geom.num_spec2_loc,
                                                     NFLD, owners, spectral_field_glob,
                                                     geom.kresol)
        gathered_glob = ectrans4py.gath_spec4py(geom.num_spec2_glob, geom.num_spec2_loc, NFLD,
                                                owners, spectral_field_loc, geom.kresol)

        identical = bool(np.array_equal(gathered_glob, spectral_field_glob))
        assert MPI.COMM_WORLD.allreduce(identical, op=MPI.LAND)

    def test_dist_grid(self):
        geom = self.geom
        owners, gridded_field_glob = self._own_fields(
            _random_grid_field(geom, NFLD, seed=5 + context.my_rank))
        gridded_field_loc = ectrans4py.dist_grid4py(geom.num_grid_points_glob,
                                                    geom.num_grid_points_loc, NFLD, owners,
                                                    gridded_field_glob, geom.kresol)
        gathered_glob = ectrans4py.gath_grid4py(geom.num_grid_points_glob,
                                                geom.num_grid_points_loc, NFLD, owners,
                                                gridded_field_loc, geom.kresol)

        identical = bool(np.array_equal(gathered_glob, gridded_field_glob))
        assert MPI.COMM_WORLD.allreduce(identical, op=MPI.LAND)


class TestNorm(DistributedTestCase):
    def test_specnorm(self):
        geom = self.geom
        # Square root of the weighted sum of squares of the coefficients (weights of
        # _spectral_dot), returned on rank 0
        spectral_field_glob = _random_spectral_field(geom, NFLD, seed=7)
        norms = ectrans4py.specnorm4py(geom.num_spec2_loc, NFLD,
                                       _distribute_spectral(geom, spectral_field_glob),
                                       geom.kresol)

        def error():
            squares = geom.spectral_weights * spectral_field_glob.astype(np.float64) ** 2
            expected = np.sqrt(np.sum(squares, axis=1))
            return float(np.max(np.abs(norms - expected) / expected))

        # Sums of squares in working precision over T + 1 wavenumbers n, then T + 1 m
        assert _on_rank_0(error) <= 2 * (geom.truncation + 1) * EPS

    def test_gpnorm_trans(self):
        geom = self.geom
        # Average (Gaussian weight times mean of each latitude row, summed over the rows),
        # minimum and maximum, returned on rank 0
        gridded_field_glob = _random_grid_field(geom, NFLD, seed=8)
        average, minimum, maximum = ectrans4py.gpnorm_trans4py(
            geom.num_grid_points_loc, NFLD, _distribute_grid(geom, gridded_field_glob),
            geom.kresol)

        def errors():
            rows = np.split(gridded_field_glob.astype(np.float64),
                            np.cumsum(geom.lons_per_lat)[:-1], axis=1)
            expected = sum(weight * row.mean(axis=1)
                           for weight, row in zip(geom.gauss_weights.astype(np.float64), rows))
            exact = (np.array_equal(minimum, np.min(gridded_field_glob, axis=1))
                     and np.array_equal(maximum, np.max(gridded_field_glob, axis=1)))
            average_error = np.max(np.abs(average - expected)) / np.max(np.abs(gridded_field_glob))
            return bool(exact), float(average_error)

        min_max_exact, average_error = _on_rank_0(errors)
        assert min_max_exact
        # The latitude contributions are summed in working precision, one rounding per row
        assert average_error <= len(geom.lons_per_lat) * EPS


class TestAdjoint(DistributedTestCase):
    # <F x, y> = <x, F* y> for each transform F and its adjoint F*, with the inner products of
    # ecTrans's adjoint tests (_grid_dot, _spectral_dot)

    def test_inv_trans_scalar(self):
        geom = self.geom
        x_loc = _distribute_spectral(geom, _random_spectral_field(geom, NFLD, seed=9))
        y_loc = _distribute_grid(geom, _random_grid_field(geom, NFLD, seed=10))
        inv_x_loc = ectrans4py.inv_trans_scalar_dist4py(
            geom.num_spec2_loc, geom.num_grid_points_loc, NFLD, x_loc, geom.kresol)
        adj_y_loc = ectrans4py.inv_trans_scalar_dist4py_ad(
            geom.num_spec2_loc, geom.num_grid_points_loc, NFLD, y_loc, geom.kresol)

        lhs = _grid_dot((inv_x_loc, y_loc))
        rhs = _spectral_dot(geom, (x_loc, adj_y_loc))
        assert _adjoint_error(lhs, rhs) <= ADJOINT_TOL

    def test_inv_trans_scalar_ders(self):
        geom = self.geom
        # The derivatives are smaller than the values by a factor of about the Earth's radius,
        # so their seeds are scaled up by it to weigh as much in the inner product
        x_loc = _distribute_spectral(geom, _random_spectral_field(geom, NFLD, seed=11))
        y_loc, y_ns_loc, y_ew_loc = [
            _distribute_grid(geom, scale * _random_grid_field(geom, NFLD, seed))
            for scale, seed in ((1.0, 12), (EARTH_RADIUS, 13), (EARTH_RADIUS, 14))]
        value_loc, ns_der_loc, ew_der_loc = ectrans4py.inv_trans_scalar_ders_dist4py(
            geom.num_spec2_loc, geom.num_grid_points_loc, NFLD, x_loc, geom.kresol)
        adj_y_loc = ectrans4py.inv_trans_scalar_ders_dist4py_ad(
            geom.num_spec2_loc, geom.num_grid_points_loc, NFLD, y_loc, y_ns_loc, y_ew_loc,
            geom.kresol)

        lhs = _grid_dot((value_loc, y_loc), (ns_der_loc, y_ns_loc), (ew_der_loc, y_ew_loc))
        rhs = _spectral_dot(geom, (x_loc, adj_y_loc))
        assert _adjoint_error(lhs, rhs) <= ADJOINT_TOL

    def test_dir_trans_scalar(self):
        geom = self.geom
        u_loc = _distribute_grid(geom, _random_grid_field(geom, NFLD, seed=15))
        v_loc = _distribute_spectral(geom, _random_spectral_field(geom, NFLD, seed=16))
        dir_u_loc = ectrans4py.dir_trans_scalar_dist4py(
            geom.num_spec2_loc, geom.num_grid_points_loc, NFLD, u_loc, geom.kresol)
        adj_v_loc = ectrans4py.dir_trans_scalar_dist4py_ad(
            geom.num_spec2_loc, geom.num_grid_points_loc, NFLD, v_loc, geom.kresol)

        lhs = _spectral_dot(geom, (dir_u_loc, v_loc))
        rhs = _grid_dot((u_loc, adj_v_loc))
        assert _adjoint_error(lhs, rhs) <= ADJOINT_TOL

    def test_inv_trans_uv(self):
        geom = self.geom
        vorticity_loc, divergence_loc = [
            _distribute_spectral(geom, _random_spectral_field(geom, NFLD, seed, zero_mean=True))
            for seed in (17, 18)]
        y_u_loc, y_v_loc = [_distribute_grid(geom, _random_grid_field(geom, NFLD, seed))
                            for seed in (19, 20)]
        u_loc, v_loc = ectrans4py.inv_trans_uv_dist4py(
            geom.num_spec2_loc, geom.num_grid_points_loc, NFLD, vorticity_loc, divergence_loc,
            geom.kresol)
        adj_vorticity_loc, adj_divergence_loc = ectrans4py.inv_trans_uv_dist4py_ad(
            geom.num_spec2_loc, geom.num_grid_points_loc, NFLD, y_u_loc, y_v_loc, geom.kresol)

        lhs = _grid_dot((u_loc, y_u_loc), (v_loc, y_v_loc))
        rhs = _spectral_dot(geom, (vorticity_loc, adj_vorticity_loc),
                            (divergence_loc, adj_divergence_loc))
        assert _adjoint_error(lhs, rhs) <= ADJOINT_TOL

    def test_dir_trans_uv(self):
        geom = self.geom
        u_loc, v_loc = [_distribute_grid(geom, _random_grid_field(geom, NFLD, seed))
                        for seed in (21, 22)]
        y_vorticity_loc, y_divergence_loc = [
            _distribute_spectral(geom, _random_spectral_field(geom, NFLD, seed, zero_mean=True))
            for seed in (23, 24)]
        vorticity_loc, divergence_loc = ectrans4py.dir_trans_uv_dist4py(
            geom.num_spec2_loc, geom.num_grid_points_loc, NFLD, u_loc, v_loc, geom.kresol)
        adj_u_loc, adj_v_loc = ectrans4py.dir_trans_uv_dist4py_ad(
            geom.num_spec2_loc, geom.num_grid_points_loc, NFLD, y_vorticity_loc,
            y_divergence_loc, geom.kresol)

        lhs = _spectral_dot(geom, (vorticity_loc, y_vorticity_loc),
                            (divergence_loc, y_divergence_loc))
        rhs = _grid_dot((u_loc, adj_u_loc), (v_loc, adj_v_loc))
        assert _adjoint_error(lhs, rhs) <= ADJOINT_TOL


class TestInputArrays(DistributedTestCase):
    def test_grid_point_inputs_unchanged(self):
        # The bindings pass grid-point inputs to ecTrans in place, without a copy, and must leave
        # them bit for bit as they were
        geom = self.geom
        sizes = (geom.num_spec2_loc, geom.num_grid_points_loc, NFLD)
        gridded_field_loc = _distribute_grid(geom, _random_grid_field(geom, NFLD, seed=27))
        original = gridded_field_loc.copy()
        ectrans4py.dir_trans_scalar_dist4py(*sizes, gridded_field_loc, geom.kresol)
        ectrans4py.inv_trans_scalar_dist4py_ad(*sizes, gridded_field_loc, geom.kresol)
        ectrans4py.gath_grid4py(geom.num_grid_points_glob, geom.num_grid_points_loc, NFLD,
                                _from_which_rank(context.my_rank, NFLD), gridded_field_loc,
                                geom.kresol)
        ectrans4py.gpnorm_trans4py(geom.num_grid_points_loc, NFLD, gridded_field_loc, geom.kresol)

        unchanged = gridded_field_loc.tobytes() == original.tobytes()
        assert MPI.COMM_WORLD.allreduce(unchanged, op=MPI.LAND)


# The tests above on the second resolution: every binding is passed its KRESOL, and runs at
# that resolution's truncation and on its grid. (The test_scalar of TestInverse and TestDirect
# use the first resolution and the default KRESOL, and are not repeated.)


class TestInquiryResolution2(TestInquiry):
    geom = context_2


class TestInverseResolution2(TestInverse):
    geom = context_2
    test_scalar = None


class TestDirectResolution2(TestDirect):
    geom = context_2
    test_scalar = None


class TestRoundTripResolution2(TestRoundTrip):
    geom = context_2


class TestDistributeGatherResolution2(TestDistributeGather):
    geom = context_2


class TestNormResolution2(TestNorm):
    geom = context_2


class TestAdjointResolution2(TestAdjoint):
    geom = context_2


class TestInputArraysResolution2(TestInputArrays):
    geom = context_2


class TestResolutions(TestCase):
    def test_default_resolution(self):
        # The first resolution set up has the tag 1, and without KRESOL every binding runs on
        # it, bit for bit as with KRESOL=1
        calls = _binding_calls(context)
        names = [name for name, _ in calls]
        assert context.kresol == 1 and context_2.kresol == 2
        assert _differing_outputs(names, [call() for _, call in calls],
                                  [call(1) for _, call in calls]) == []

    def test_first_resolution_unchanged(self):
        # After the second resolution has been set up and used, the outputs of the first one
        # are bit-identical to those before the second one was set up
        for _, call in _binding_calls(context_2):
            call(context_2.kresol)
        calls = _binding_calls(context)
        names = [name for name, _ in calls]
        assert _differing_outputs(names, context.reference_outputs,
                                  [call() for _, call in calls]) == []

    def test_unknown_resolution(self):
        # Every binding raises error code 101 for a KRESOL that is not a resolution set up: 0
        # and KMAX_RESOL + 1, out of range, and KMAX_RESOL, never set up. (trans_inq4py takes
        # KRESOL <= 0 to set up a serial resolution from its other arguments.)
        for kresol in (0, KMAX_RESOL, KMAX_RESOL + 1):
            for _, call in _binding_calls(context):
                with pytest.raises(RuntimeError, match="Error code 101"):
                    call(kresol)
            if kresol > 0:
                with pytest.raises(RuntimeError, match="Error code 101"):
                    ectrans4py.trans_inq4py(kresol, NUM_LATS, 0, 0, np.array([], dtype=np.int64),
                                            0)
