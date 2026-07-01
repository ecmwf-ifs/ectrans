"""Distributed-memory ectrans4py smoke test (self-contained, no external data).

Run:  mpirun -n 2 python tests/test_ectrans4py/mpi_distributed.py   (also -n 1)

Exercises the MPI interface end to end on a small analytic Gaussian grid:
  mpl_init4py -> setup_trans0_4py (processor grid) -> setup_trans_4py (LDSPLIT)
  -> trans_inq4py(KRESOL>0)  [the unified inquiry, keyed by the resolution handle]
  -> dist_spec4py / inv/dir_trans_scalar_dist4py / gath_grid4py round trip.

Checks: truncation & global sizes; local grid sizes sum to the global size across
tasks; the (0,0) spectral coefficient maps to a constant grid field; and the
grid -> dir -> inv -> grid round trip is the identity.
"""
import sys

import numpy as np
from mpi4py import MPI

import ectrans4py

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

ectrans4py.init_env(unlimited_stack=False)
ectrans4py.mpl_init4py()

# field-data dtype matches the loaded build precision (float32 for sp, float64 for dp)
REAL = getattr(ectrans4py, "_REAL", np.float64)

# --- small analytic full Gaussian grid: NDGL=32 latitudes, 64 longitudes, T21 ---
NDGL = 32
KSMAX = 21
KLOEN = np.full(NDGL, 2 * NDGL, dtype=np.int64)
NGPTOTG = int(KLOEN.sum())

# processor grid (N-S split), distributed
ectrans4py.setup_trans0_4py(size, 1, size, False, 10)
kresol = ectrans4py.setup_trans_4py(KSMAX, NDGL, NDGL, KLOEN, True, False)

(kgptot, kspec, kspec2, kgptotg, kspec2g, ksmax, knmeng, pmu, pgw) = \
    ectrans4py.trans_inq4py(kresol, NDGL, KSMAX, NDGL, KLOEN, 10)

ok = True


def check(cond, msg):
    global ok
    if not cond:
        ok = False
        if rank == 0:
            print(f"[dist] FAIL: {msg}")


check(ksmax == KSMAX, f"ksmax {ksmax} != {KSMAX}")
check(kgptotg == NGPTOTG, f"kgptotg {kgptotg} != {NGPTOTG}")
check(np.all(np.abs(pmu) <= 1.0 + 1e-12), "pmu out of [-1,1]")
check(np.all(np.diff(pmu) > 0) or np.all(np.diff(pmu) < 0), "pmu not monotone")
# local grid-point counts partition the global grid
check(comm.allreduce(kgptot, op=MPI.SUM) == kgptotg, "local kgptot do not sum to global")

# --- (0,0) coefficient -> constant grid field --------------------------------
nfld = 1
kfrom = np.ones(nfld, dtype=np.int64)
specg = np.zeros((nfld, kspec2g), dtype=REAL)
if rank == 0:
    specg[0, 0] = 1.0                       # (m=0, n=0) real part -> constant field
sloc = ectrans4py.dist_spec4py(kspec2g, kspec2, nfld, kfrom, specg)
g = ectrans4py.inv_trans_scalar_dist4py(kspec2, kgptot, nfld, sloc)
gg = ectrans4py.gath_grid4py(kgptotg, kgptot, nfld, kfrom, g)
if rank == 0:
    rng = float(np.max(gg) - np.min(gg))
    check(rng < 1e-4, f"(0,0) coeff did not give a constant field (range {rng:.2e})")

# --- grid -> dir -> inv -> grid identity -------------------------------------
s2 = ectrans4py.dir_trans_scalar_dist4py(kspec2, kgptot, nfld, g)
g2 = ectrans4py.inv_trans_scalar_dist4py(kspec2, kgptot, nfld, s2)
err = comm.allreduce(float(np.max(np.abs(g2 - g))), op=MPI.MAX)
check(err < 1e-4, f"grid->dir->inv round trip error {err:.2e}")

if rank == 0:
    print(f"[dist] size={size} kgptot={kgptot} kspec2={kspec2} kgptotg={kgptotg} "
          f"ksmax={ksmax} pgw_sum={float(pgw.sum()):.6f} rt_err={err:.2e}")
    print("[dist] PASS" if ok else "[dist] FAIL")

ectrans4py.mpl_end4py()
if rank == 0 and not ok:
    sys.exit(1)
