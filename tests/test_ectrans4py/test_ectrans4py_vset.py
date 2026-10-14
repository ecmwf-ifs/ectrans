# (C) Copyright 2026- ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.

# ectrans4py with the spectral fields distributed over V-sets as well as over wave sets (the 2-D
# spectral decomposition): the *_vset4py bindings. ecTrans takes its processor decomposition when
# it is set up, once per process, so this module runs in an MPI job of its own:
#     mpirun -np 4 python -m pytest -k test_ectrans4py_vset
# With an even number of tasks NPROC, the spectral coefficients are distributed over NPROC / 2
# wave sets and 2 V-sets, with an odd number over NPROC wave sets and a single V-set. The
# environment variable ECTRANS4PY_TEST_NPRTRW sets another number of wave sets NPRTRW, a divisor
# of NPROC (the number of V-sets is NPROC / NPRTRW).

import os
import platform
from types import SimpleNamespace
from unittest import TestCase

import numpy as np
import pytest

import ectrans4py

mpi4py = pytest.importorskip("mpi4py")
from mpi4py import MPI

TRUNCATION = 21
NUM_LATS = 32
LONS_PER_LAT = np.full(NUM_LATS, 2 * NUM_LATS, dtype=np.int64)
# Number of global fields: not a multiple of the number of V-sets, which then hold different
# numbers of fields
NFLD = 7

# Earth radius used by ecTrans (the default of SETUP_TRANS0's PRAD)
EARTH_RADIUS = 6371229.0

# Tolerances, relative to the size of the expected values. EPS: machine epsilon of the build.
EPS = np.finfo(ectrans4py._REAL).eps
# On a regular quadratic Gaussian grid, transforms are exact up to rounding: of a single
# spherical harmonic, there and back, and with another map from the fields to the V-sets, which
# changes the order of some of ecTrans's sums (by up to about 10 EPS, when measured)
TOL = 1000 * EPS
# Adjoint identities: the threshold of ecTrans's own adjoint tests (tests/trans/adjoint)
ADJOINT_TOL = 20000 * EPS

# Maps from the global fields to the V-sets (see _vset_map)
VSET_MAPS = ("round-robin", "blocked", "last")

# Transform geometry and decomposition, set up once by setUpModule
context = SimpleNamespace()


def _vset_map(kind):
    """V-set (1-based) of each global field: round-robin (field j in V-set j mod NPRTRV + 1),
    blocked (consecutive fields in the same V-set, as the IFS distributes its levels) or last
    (all fields in the last V-set, the other V-sets holding none)."""
    nprtrv = context.nprtrv
    if kind == "round-robin":
        kvset = np.arange(NFLD) % nprtrv + 1
    elif kind == "blocked":
        kvset = np.arange(NFLD) * nprtrv // NFLD + 1
    else:
        kvset = np.full(NFLD, nprtrv)
    return kvset.astype(np.int64)


def _num_local_fields(kvset):
    """Number of fields of this task's V-set, the KFLDL of the *_vset4py bindings"""
    return int(np.count_nonzero(kvset == context.mysetv))


def _sizes(kvset):
    """The leading arguments of the *_vset4py transforms, for the map kvset: KSPEC2, KGPTOT,
    KFLDL, KFLDG and KVSETSC or KVSETUV"""
    return (context.num_spec2_loc, context.num_grid_points_loc, _num_local_fields(kvset), NFLD,
            kvset)


def _root(nfld):
    """MPL rank (1-based) of the task holding the global fields: the first"""
    return np.ones(nfld, dtype=np.int64)


def _distribute_spectral(field_glob, kvset):
    return ectrans4py.dist_spec_vset4py(context.num_spec2_glob, context.num_spec2_loc, NFLD,
                                        _num_local_fields(kvset), _root(NFLD), kvset, field_glob)


def _gather_spectral(field_loc, kvset):
    return ectrans4py.gath_spec_vset4py(context.num_spec2_glob, context.num_spec2_loc, NFLD,
                                        _num_local_fields(kvset), _root(NFLD), kvset, field_loc)


def _distribute_grid(field_glob):
    return ectrans4py.dist_grid4py(context.num_grid_points_glob, context.num_grid_points_loc,
                                   NFLD, _root(NFLD), field_glob)


def _random_field(shape, seed):
    """Uniform random values in [-1, 1) in working precision, the same on every rank."""
    return np.random.default_rng(seed).uniform(-1.0, 1.0, shape).astype(ectrans4py._REAL)


def _random_spectral_field(seed, zero_mean=False):
    """Random global spectral field of real fields (no imaginary part for m = 0). With
    zero_mean, the (0, 0) coefficient is zero too, as for vorticity and divergence."""
    field = _random_field((NFLD, context.num_spec2_glob), seed)
    field[:, context.spectral_weights == 0.0] = 0.0
    if zero_mean:
        field[:, context.total_wavenumber == 0] = 0.0
    return field


def _random_grid_field(seed):
    return _random_field((NFLD, context.num_grid_points_glob), seed)


def _on_rank_0(compute):
    """compute() on rank 0, which holds the gathered fields, broadcast so that all ranks agree."""
    return MPI.COMM_WORLD.bcast(compute() if MPI.COMM_WORLD.Get_rank() == 0 else None, root=0)


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


def _spectral_dot(kvset, *pairs):
    """Spectral inner product summed over pairs of local fields of the map kvset, weighted as
    SCALPRODSP of tests/trans/adjoint."""
    gathered = [(_gather_spectral(a, kvset), _gather_spectral(b, kvset)) for a, b in pairs]
    return _on_rank_0(lambda: float(sum(np.sum(context.spectral_weights * a.astype(np.float64)
                                               * b.astype(np.float64)) for a, b in gathered)))


def _adjoint_error(lhs, rhs):
    """Relative difference between <F x, y> and <x, F* y>, as in tests/trans/adjoint."""
    return abs(lhs - rhs) / abs(lhs)


def _binding_calls():
    """Every V-set binding, as (name, call) pairs: call(*kresol) runs the binding on fixed
    inputs, with *kresol appended to its arguments (none: the default resolution), and returns
    its outputs."""
    kvset = _vset_map("round-robin")
    sizes = _sizes(kvset)
    spectral_glob = _random_spectral_field(seed=30, zero_mean=True)
    spectral_loc = _distribute_spectral(spectral_glob, kvset)
    grid_loc = _distribute_grid(_random_grid_field(seed=31))
    spec2_glob, spec2_loc, num_loc = (context.num_spec2_glob, context.num_spec2_loc,
                                      _num_local_fields(kvset))
    return [
        ("trans_inq_vset4py", lambda *k: ectrans4py.trans_inq_vset4py(*k)),
        ("dist_spec_vset4py", lambda *k: ectrans4py.dist_spec_vset4py(
            spec2_glob, spec2_loc, NFLD, num_loc, _root(NFLD), kvset, spectral_glob, *k)),
        ("gath_spec_vset4py", lambda *k: ectrans4py.gath_spec_vset4py(
            spec2_glob, spec2_loc, NFLD, num_loc, _root(NFLD), kvset, spectral_loc, *k)),
        ("inv_trans_scalar_vset4py", lambda *k: ectrans4py.inv_trans_scalar_vset4py(
            *sizes, spectral_loc, *k)),
        ("inv_trans_scalar_ders_vset4py", lambda *k: ectrans4py.inv_trans_scalar_ders_vset4py(
            *sizes, spectral_loc, *k)),
        ("dir_trans_scalar_vset4py", lambda *k: ectrans4py.dir_trans_scalar_vset4py(
            *sizes, grid_loc, *k)),
        ("inv_trans_uv_vset4py", lambda *k: ectrans4py.inv_trans_uv_vset4py(
            *sizes, spectral_loc, spectral_loc, *k)),
        ("dir_trans_uv_vset4py", lambda *k: ectrans4py.dir_trans_uv_vset4py(
            *sizes, grid_loc, grid_loc, *k)),
        ("inv_trans_scalar_vset4py_ad", lambda *k: ectrans4py.inv_trans_scalar_vset4py_ad(
            *sizes, grid_loc, *k)),
        ("inv_trans_scalar_ders_vset4py_ad",
         lambda *k: ectrans4py.inv_trans_scalar_ders_vset4py_ad(
             *sizes, grid_loc, grid_loc, grid_loc, *k)),
        ("dir_trans_scalar_vset4py_ad", lambda *k: ectrans4py.dir_trans_scalar_vset4py_ad(
            *sizes, spectral_loc, *k)),
        ("inv_trans_uv_vset4py_ad", lambda *k: ectrans4py.inv_trans_uv_vset4py_ad(
            *sizes, grid_loc, grid_loc, *k)),
        ("dir_trans_uv_vset4py_ad", lambda *k: ectrans4py.dir_trans_uv_vset4py_ad(
            *sizes, spectral_loc, spectral_loc, *k)),
        ("specnorm_vset4py", lambda *k: ectrans4py.specnorm_vset4py(
            spec2_loc, NFLD, num_loc, kvset, spectral_loc, *k)),
    ]


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
    nprtrw = int(os.environ.get("ECTRANS4PY_TEST_NPRTRW",
                                n_rank // 2 if n_rank % 2 == 0 else n_rank))
    # (ecTrans would abort the MPI job)
    if nprtrw < 1 or n_rank % nprtrw != 0 or n_rank // nprtrw > NFLD:
        raise ValueError(f"{nprtrw} wave sets on {n_rank} tasks: the number of wave sets must "
                         f"divide the number of tasks, into at most {NFLD} V-sets")
    ectrans4py.mpl_init4py()

    # KPRGPNS=n_rank, KPRGPEW=1, KPRTRW=nprtrw, LDEQ_REGIONS=False, KMAX_RESOL=1, LDMPOFF=False
    ectrans4py.setup_trans0_4py(n_rank, 1, nprtrw, False, 1, False)

    # KSMAX=TRUNCATION, KDGL=NUM_LATS, KSLOEN=NUM_LATS, KLOEN=LONS_PER_LAT, LDSPLIT=True,
    # LDUSEFLT=False
    context.kresol = ectrans4py.setup_trans_4py(TRUNCATION, NUM_LATS, NUM_LATS, LONS_PER_LAT,
                                                True, False)
    (context.num_grid_points_loc, _, context.num_spec2_loc, context.num_grid_points_glob,
     context.num_spec2_glob, *_) = ectrans4py.trans_inq4py(
        context.kresol, NUM_LATS, 0, 0, np.array([], dtype=np.int64), 0)
    context.nprtrw = nprtrw
    context.nprtrv = n_rank // nprtrw
    context.kprtrw, context.mysetw, context.mysetv = ectrans4py.trans_inq_vset4py(context.kresol)

    # Spectral coefficients: m = 0..T, then n = m..T, each a (real, imaginary) pair. Weights for
    # inner products and norms: 2 for m > 0 (a coefficient also stands for its conjugate), 0 for
    # the imaginary part of m = 0 (a real field has none), 1 for the real part of m = 0
    zonal, total = np.array([(m, n) for m in range(TRUNCATION + 1)
                             for n in range(m, TRUNCATION + 1)]).T
    imaginary = np.tile([False, True], zonal.size)
    context.total_wavenumber = np.repeat(total, 2)
    context.spectral_weights = np.where(np.repeat(zonal, 2) > 0, 2.0,
                                        np.where(imaginary, 0.0, 1.0))


def tearDownModule():
    ectrans4py.mpl_end4py()


class TestDecomposition(TestCase):
    def test_vset_position(self):
        # The tasks' wave sets and V-sets: each pair once, on the NPRTRW x NPRTRV grid
        positions = MPI.COMM_WORLD.allgather((context.mysetw, context.mysetv))
        layout = f"{context.nprtrw} wave sets and {context.nprtrv} V-sets"

        assert context.kprtrw == context.nprtrw, layout
        assert sorted(positions) == [(w, v) for w in range(1, context.nprtrw + 1)
                                     for v in range(1, context.nprtrv + 1)], layout

    def test_local_fields(self):
        # A task holds the fields of its V-set, in the order of the global fields, and the tasks
        # of a V-set hold every spectral coefficient once between them. The values of the
        # distributed fields give their global field and coefficient (exactly, in any precision).
        kvset = _vset_map("round-robin")
        coefficients = np.arange(context.num_spec2_glob)
        index_glob = 1000 * np.arange(NFLD)[:, np.newaxis] + coefficients
        index_loc = _distribute_spectral(index_glob.astype(ectrans4py._REAL), kvset)
        index_loc = index_loc.astype(np.int64)
        own_fields = np.flatnonzero(kvset == context.mysetv)
        local_coefficients = index_loc % 1000
        layout_ok = bool(np.all(index_loc // 1000 == own_fields[:, np.newaxis])
                         and np.all(local_coefficients == local_coefficients[0]))
        held = MPI.COMM_WORLD.allgather((context.mysetv, local_coefficients[0]))

        assert MPI.COMM_WORLD.allreduce(layout_ok, op=MPI.LAND)
        for vset in range(1, context.nprtrv + 1):
            held_by_vset = np.concatenate([c for v, c in held if v == vset])
            assert np.array_equal(np.sort(held_by_vset), coefficients)


class TestDistributeGather(TestCase):
    def test_round_trip(self):
        # Distributing and gathering only move data: bit for bit the global fields, whatever the
        # map from fields to V-sets (GATH_SPEC returns no imaginary part for m = 0, so the fields
        # have none)
        spectral_field_glob = _random_spectral_field(seed=1)
        gathered = []
        for kind in VSET_MAPS:
            kvset = _vset_map(kind)
            gathered.append(_gather_spectral(_distribute_spectral(spectral_field_glob, kvset),
                                             kvset))

        assert _on_rank_0(lambda: all(np.array_equal(field, spectral_field_glob)
                                      for field in gathered))


class TestInverse(TestCase):
    def test_constant(self):
        # The coefficient (n = 0, m = 0) alone, 1, gives the constant field 1, whatever the V-set
        # of the field
        spectral_field_glob = np.zeros((NFLD, context.num_spec2_glob), dtype=ectrans4py._REAL)
        spectral_field_glob[:, 0] = 1.0
        kvset = _vset_map("round-robin")
        gridded_field_loc = ectrans4py.inv_trans_scalar_vset4py(
            *_sizes(kvset), _distribute_spectral(spectral_field_glob, kvset), context.kresol)
        error = float(np.max(np.abs(gridded_field_loc.astype(np.float64) - 1.0)))

        assert MPI.COMM_WORLD.allreduce(error, op=MPI.MAX) <= TOL


class TestRoundTrip(TestCase):
    def test_scalar(self):
        # spectral -> grid -> spectral
        kvset = _vset_map("round-robin")
        spectral_field_glob = _random_spectral_field(seed=2)
        gridded_field_loc = ectrans4py.inv_trans_scalar_vset4py(
            *_sizes(kvset), _distribute_spectral(spectral_field_glob, kvset), context.kresol)
        spectral_field_loc = ectrans4py.dir_trans_scalar_vset4py(
            *_sizes(kvset), gridded_field_loc, context.kresol)
        round_trip_glob = _gather_spectral(spectral_field_loc, kvset)

        assert _on_rank_0(lambda: _relative_error(round_trip_glob, spectral_field_glob)) <= TOL

    def test_uv(self):
        # vorticity, divergence -> u, v -> vorticity, divergence
        kvset = _vset_map("round-robin")
        vorticity_glob = _random_spectral_field(seed=3, zero_mean=True)
        divergence_glob = _random_spectral_field(seed=4, zero_mean=True)
        u_loc, v_loc = ectrans4py.inv_trans_uv_vset4py(
            *_sizes(kvset), _distribute_spectral(vorticity_glob, kvset),
            _distribute_spectral(divergence_glob, kvset), context.kresol)
        vorticity_loc, divergence_loc = ectrans4py.dir_trans_uv_vset4py(
            *_sizes(kvset), u_loc, v_loc, context.kresol)
        round_trip_glob = [_gather_spectral(vorticity_loc, kvset),
                           _gather_spectral(divergence_loc, kvset)]

        assert _on_rank_0(lambda: _relative_error(round_trip_glob,
                                                  [vorticity_glob, divergence_glob])) <= TOL

    def test_scalar_ders(self):
        # The E-W and N-S derivatives of a velocity potential chi are the winds u and v of a
        # purely divergent flow: transformed back, they have no vorticity and the divergence
        # -n (n + 1) chi / a^2
        kvset = _vset_map("round-robin")
        chi_glob = _random_spectral_field(seed=5)
        _, ns_der_loc, ew_der_loc = ectrans4py.inv_trans_scalar_ders_vset4py(
            *_sizes(kvset), _distribute_spectral(chi_glob, kvset), context.kresol)
        vorticity_loc, divergence_loc = ectrans4py.dir_trans_uv_vset4py(
            *_sizes(kvset), ew_der_loc, ns_der_loc, context.kresol)
        vorticity_glob = _gather_spectral(vorticity_loc, kvset)
        divergence_glob = _gather_spectral(divergence_loc, kvset)

        def error():
            n = context.total_wavenumber
            laplacian = -n * (n + 1) / EARTH_RADIUS ** 2 * chi_glob.astype(np.float64)
            return _relative_error([vorticity_glob, divergence_glob],
                                   [np.zeros_like(laplacian), laplacian])

        assert _on_rank_0(error) <= TOL


class TestNorm(TestCase):
    def test_specnorm_vset(self):
        # The norms of all global fields, from the fields of each V-set, on the first task:
        # square roots of the weighted sums of squares of the coefficients (weights of
        # _spectral_dot), whatever the map from fields to V-sets
        spectral_field_glob = _random_spectral_field(seed=6)
        norms = []
        for kind in VSET_MAPS:
            kvset = _vset_map(kind)
            norms.append(ectrans4py.specnorm_vset4py(
                context.num_spec2_loc, NFLD, _num_local_fields(kvset), kvset,
                _distribute_spectral(spectral_field_glob, kvset), context.kresol))

        def error():
            squares = context.spectral_weights * spectral_field_glob.astype(np.float64) ** 2
            expected = np.sqrt(np.sum(squares, axis=1))
            return max(float(np.max(np.abs(norm - expected) / expected)) for norm in norms)

        # Sums of squares in working precision over T + 1 wavenumbers n, then T + 1 m
        assert _on_rank_0(error) <= 2 * (TRUNCATION + 1) * EPS


class TestAdjoint(TestCase):
    # <F x, y> = <x, F* y> for each transform F and its adjoint F*, with the inner products of
    # ecTrans's adjoint tests (_grid_dot, _spectral_dot)

    def test_inv_trans_scalar(self):
        kvset = _vset_map("round-robin")
        x_loc = _distribute_spectral(_random_spectral_field(seed=7), kvset)
        y_loc = _distribute_grid(_random_grid_field(seed=8))
        inv_x_loc = ectrans4py.inv_trans_scalar_vset4py(*_sizes(kvset), x_loc, context.kresol)
        adj_y_loc = ectrans4py.inv_trans_scalar_vset4py_ad(*_sizes(kvset), y_loc, context.kresol)

        lhs = _grid_dot((inv_x_loc, y_loc))
        rhs = _spectral_dot(kvset, (x_loc, adj_y_loc))
        assert _adjoint_error(lhs, rhs) <= ADJOINT_TOL

    def test_inv_trans_scalar_ders(self):
        # The derivatives are smaller than the values by a factor of about the Earth's radius,
        # so their seeds are scaled up by it to weigh as much in the inner product
        kvset = _vset_map("round-robin")
        x_loc = _distribute_spectral(_random_spectral_field(seed=9), kvset)
        y_loc, y_ns_loc, y_ew_loc = [
            _distribute_grid(scale * _random_grid_field(seed))
            for scale, seed in ((1.0, 10), (EARTH_RADIUS, 11), (EARTH_RADIUS, 12))]
        value_loc, ns_der_loc, ew_der_loc = ectrans4py.inv_trans_scalar_ders_vset4py(
            *_sizes(kvset), x_loc, context.kresol)
        adj_y_loc = ectrans4py.inv_trans_scalar_ders_vset4py_ad(
            *_sizes(kvset), y_loc, y_ns_loc, y_ew_loc, context.kresol)

        lhs = _grid_dot((value_loc, y_loc), (ns_der_loc, y_ns_loc), (ew_der_loc, y_ew_loc))
        rhs = _spectral_dot(kvset, (x_loc, adj_y_loc))
        assert _adjoint_error(lhs, rhs) <= ADJOINT_TOL

    def test_dir_trans_scalar(self):
        kvset = _vset_map("round-robin")
        u_loc = _distribute_grid(_random_grid_field(seed=13))
        v_loc = _distribute_spectral(_random_spectral_field(seed=14), kvset)
        dir_u_loc = ectrans4py.dir_trans_scalar_vset4py(*_sizes(kvset), u_loc, context.kresol)
        adj_v_loc = ectrans4py.dir_trans_scalar_vset4py_ad(*_sizes(kvset), v_loc, context.kresol)

        lhs = _spectral_dot(kvset, (dir_u_loc, v_loc))
        rhs = _grid_dot((u_loc, adj_v_loc))
        assert _adjoint_error(lhs, rhs) <= ADJOINT_TOL

    def test_inv_trans_uv(self):
        kvset = _vset_map("round-robin")
        vorticity_loc, divergence_loc = [
            _distribute_spectral(_random_spectral_field(seed, zero_mean=True), kvset)
            for seed in (15, 16)]
        y_u_loc, y_v_loc = [_distribute_grid(_random_grid_field(seed)) for seed in (17, 18)]
        u_loc, v_loc = ectrans4py.inv_trans_uv_vset4py(
            *_sizes(kvset), vorticity_loc, divergence_loc, context.kresol)
        adj_vorticity_loc, adj_divergence_loc = ectrans4py.inv_trans_uv_vset4py_ad(
            *_sizes(kvset), y_u_loc, y_v_loc, context.kresol)

        lhs = _grid_dot((u_loc, y_u_loc), (v_loc, y_v_loc))
        rhs = _spectral_dot(kvset, (vorticity_loc, adj_vorticity_loc),
                            (divergence_loc, adj_divergence_loc))
        assert _adjoint_error(lhs, rhs) <= ADJOINT_TOL

    def test_dir_trans_uv(self):
        kvset = _vset_map("round-robin")
        u_loc, v_loc = [_distribute_grid(_random_grid_field(seed)) for seed in (19, 20)]
        y_vorticity_loc, y_divergence_loc = [
            _distribute_spectral(_random_spectral_field(seed, zero_mean=True), kvset)
            for seed in (21, 22)]
        vorticity_loc, divergence_loc = ectrans4py.dir_trans_uv_vset4py(
            *_sizes(kvset), u_loc, v_loc, context.kresol)
        adj_u_loc, adj_v_loc = ectrans4py.dir_trans_uv_vset4py_ad(
            *_sizes(kvset), y_vorticity_loc, y_divergence_loc, context.kresol)

        lhs = _spectral_dot(kvset, (vorticity_loc, y_vorticity_loc),
                            (divergence_loc, y_divergence_loc))
        rhs = _grid_dot((u_loc, adj_u_loc), (v_loc, adj_v_loc))
        assert _adjoint_error(lhs, rhs) <= ADJOINT_TOL


class TestFieldMaps(TestCase):
    # The same global fields, distributed over the V-sets by different maps, give the same
    # results up to rounding. With the last map, the tasks of the V-sets without fields still
    # transform the grid-point values of all fields (and their derivatives, for the adjoint of
    # INV_TRANS with LDSCDERS).

    def test_spectral_to_grid(self):
        # The inverse transforms and the adjoints of the direct transforms, on every task
        spectral_field_glob = _random_spectral_field(seed=23, zero_mean=True)
        outputs = []
        for kind in VSET_MAPS:
            kvset = _vset_map(kind)
            sizes = _sizes(kvset)
            x_loc = _distribute_spectral(spectral_field_glob, kvset)
            outputs.append(
                [ectrans4py.inv_trans_scalar_vset4py(*sizes, x_loc, context.kresol)]
                + list(ectrans4py.inv_trans_scalar_ders_vset4py(*sizes, x_loc, context.kresol))
                + list(ectrans4py.inv_trans_uv_vset4py(*sizes, x_loc, x_loc, context.kresol))
                + [ectrans4py.dir_trans_scalar_vset4py_ad(*sizes, x_loc, context.kresol)]
                + list(ectrans4py.dir_trans_uv_vset4py_ad(*sizes, x_loc, x_loc,
                                                          context.kresol)))
        error = max(_relative_error(field, reference) for other in outputs[1:]
                    for field, reference in zip(other, outputs[0]))

        assert MPI.COMM_WORLD.allreduce(error, op=MPI.MAX) <= TOL

    def test_grid_to_spectral(self):
        # The direct transforms and the adjoints of the inverse transforms, gathered on the
        # first task
        grid_loc = _distribute_grid(_random_grid_field(seed=24))
        outputs = []
        for kind in VSET_MAPS:
            kvset = _vset_map(kind)
            sizes = _sizes(kvset)
            fields_loc = (
                [ectrans4py.dir_trans_scalar_vset4py(*sizes, grid_loc, context.kresol)]
                + list(ectrans4py.dir_trans_uv_vset4py(*sizes, grid_loc, grid_loc,
                                                       context.kresol))
                + [ectrans4py.inv_trans_scalar_vset4py_ad(*sizes, grid_loc, context.kresol),
                   ectrans4py.inv_trans_scalar_ders_vset4py_ad(*sizes, grid_loc, grid_loc,
                                                               grid_loc, context.kresol)]
                + list(ectrans4py.inv_trans_uv_vset4py_ad(*sizes, grid_loc, grid_loc,
                                                          context.kresol)))
            outputs.append([_gather_spectral(field, kvset) for field in fields_loc])

        assert _on_rank_0(lambda: max(_relative_error(field, reference) for other in outputs[1:]
                                      for field, reference in zip(other, outputs[0]))) <= TOL


class TestResolution(TestCase):
    def test_unknown_resolution(self):
        # Every V-set binding raises error code 101 for a KRESOL that is not a resolution set up
        for kresol in (0, context.kresol + 1):
            for _, call in _binding_calls():
                with pytest.raises(RuntimeError, match="Error code 101"):
                    call(kresol)

    def test_default_resolution(self):
        # Without KRESOL, every V-set binding runs on the first resolution, bit for bit as with
        # KRESOL=1
        def as_tuple(output):
            return output if isinstance(output, tuple) else (output,)

        differing = [name for name, call in _binding_calls()
                     if not all(np.array_equal(a, b)
                                for a, b in zip(as_tuple(call()), as_tuple(call(1))))]

        assert context.kresol == 1
        assert sum(MPI.COMM_WORLD.allgather(differing), []) == []
