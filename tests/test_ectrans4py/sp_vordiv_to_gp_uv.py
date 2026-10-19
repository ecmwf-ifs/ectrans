import numpy as np
import ectrans4py

# field-data dtype matches the loaded build precision (float32 for sp, float64 for dp)
REAL = getattr(ectrans4py, "_REAL", np.float64)

# --- T21 octahedral reduced Gaussian grid ---
truncation = 639
nlat = (truncation + 1) * 2
lons_per_lat = [20 + 4 * i for i in range(nlat // 2)]
lons_per_lat = np.array(lons_per_lat + lons_per_lat[::-1], dtype=np.int64)
n_grid_points = sum(lons_per_lat)

ectrans4py.setup_trans0_4py(1, 1, 1, False, 1, True)
kresol = ectrans4py.setup_trans_4py(truncation, nlat, nlat, lons_per_lat, True, False)

kgptot, kspec, kspec2, kgptotg, kspec2g, trunc_ectrans, knmeng, pmu, pgw = ectrans4py.trans_inq4py(
    kresol, nlat, truncation, nlat, lons_per_lat, 1
)

assert truncation == trunc_ectrans, f"truncation {truncation} != {trunc_ectrans}"
assert np.all(np.abs(pmu) <= 1.0 + 1e-12), "pmu out of [-1,1]"
assert np.all(np.diff(pmu) > 0) or np.all(np.diff(pmu) < 0), "pmu not monotone"
assert kgptot == n_grid_points, f"kgptot {kgptot} != n_grid_points {n_grid_points}"
assert kgptotg == n_grid_points, f"kgptotg {kgptotg} != n_grid_points {n_grid_points}"

# --- (0,0) coefficient -> constant grid field --------------------------------
nfld = 100
vor_spec = np.zeros((nfld, kspec2g), dtype=REAL)
div_spec = np.zeros((nfld, kspec2g), dtype=REAL)
# (m=0, n=0) real part -> constant field
vor_spec[0, 0] = 1.0
div_spec[0, 0] = 1.0

print("Computing inverse transform to grid...")

u_grid, v_grid = ectrans4py.inv_trans_uv_dist4py(kspec2, kgptot, nfld, vor_spec, div_spec)

print(kgptotg)
print(u_grid.shape, v_grid.shape)
# rng = float(np.max(u_grid) - np.min(u_grid))
# check(rng < 1e-4, f"(0,0) coeff did not give a constant field (range {rng:.2e})")

# err = float(np.max(np.abs(g2 - g)))
# check(err < 1e-4, f"grid->dir->inv round trip error {err:.2e}")

# print(f"[dist] kgptot={kgptot} kspec2={kspec2} kgptotg={kgptotg} "
#         f"truncation={truncation} pgw_sum={float(pgw.sum()):.6f} rt_err={err:.2e}")
# print("[dist] PASS" if ok else "[dist] FAIL")

