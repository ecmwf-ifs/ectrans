#include <unordered_map>
#include <vector>
#include <algorithm>
#include <memory>
#include <fftw3.h>

typedef enum fftType_t {
  FFT_R2C = 0x2a,  // Real to complex (interleaved)
  FFT_C2R = 0x2c,  // Complex (interleaved) to real
  FFT_D2Z = 0x6a,  // Double to double-complex (interleaved)
  FFT_Z2D = 0x6c,  // Double-complex (interleaved) to double
} fftType;

namespace {
struct Double {
  using real = double;
  using cmplx = fftw_complex;
  using plan_type = fftw_plan;
};
struct Float {
  using real = float;
  using cmplx = fftwf_complex;
  using plan_type = fftwf_plan;
};

template <class Type, fftType Direction> class fft_plan {
  using real = typename Type::real;
  using cmplx = typename Type::cmplx;
  using plan_type = typename Type::plan_type;

public:
  void exec(real *data_real, cmplx *data_complex) const {
    real *data_real_l = &data_real[offset];
    cmplx *data_complex_l = &data_complex[offset / 2];
    if constexpr (Direction == FFT_R2C)
      fftwf_execute_dft_r2c(*plan_ptr, data_real_l, data_complex_l);
    else if constexpr (Direction == FFT_C2R)
      fftwf_execute_dft_c2r(*plan_ptr, data_complex_l, data_real_l);
    else if constexpr (Direction == FFT_D2Z)
      fftw_execute_dft_r2c(*plan_ptr, data_real_l, data_complex_l);
    else if constexpr (Direction == FFT_Z2D)
      fftw_execute_dft_c2r(*plan_ptr, data_complex_l, data_real_l);
  }
  fft_plan(plan_type plan_, int64_t offset_)
      : plan_ptr(new plan_type{plan_},
                   [](auto ptr) {
                     if constexpr (Direction == FFT_R2C || Direction == FFT_C2R)
                       fftwf_destroy_plan(*ptr);
                     else if constexpr (Direction == FFT_D2Z || Direction == FFT_Z2D)
                       fftw_destroy_plan(*ptr);
                     delete ptr;
                   }),
        offset(offset_) {}

private:
  std::shared_ptr<plan_type> plan_ptr;
  int64_t offset;
};

struct cache_key {
  int resol_id;
  int kfield;
  bool operator==(const cache_key &other) const {
    return resol_id == other.resol_id && kfield == other.kfield;
  }
  cache_key(int resol_id_, int kfield_)
      : resol_id(resol_id_), kfield(kfield_) {}
};
} // namespace

template <> struct std::hash<cache_key> {
  std::size_t operator()(const cache_key &k) const {
    return k.resol_id * 10000 + k.kfield;
  }
};

namespace {

// kfield -> handles
template <class Type, fftType Direction> auto &get_fft_plan_cache() {
  static std::unordered_map<cache_key,
                            std::vector<fft_plan<Type, Direction>>>
      fftPlansCache;
  return fftPlansCache;
}

// // kfield -> ptrs
// template <class Type, hipfftType Direction> auto &get_ptr_cache() {
//   using real = typename Type::real;
//   using cmplx = typename Type::cmplx;
//   static std::unordered_map<cache_key, std::pair<real *, cmplx *>> ptrCache;
//   return ptrCache;
// }

template <typename Cache>
void erase_resol_from_cache(Cache &cache, int resol_id) {
  // Note that in C++20 this could also be std::erase_if
  int erased = 0;
  for (auto it = cache.begin(); it != cache.end();) {
    if (it->first.resol_id == resol_id) {
      it = cache.erase(it);
      ++erased;
    } else
      ++it;
  }
}
template <class Type, fftType Direction>
void erase_from_caches(int resol_id) {
  erase_resol_from_cache(get_fft_plan_cache<Type, Direction>(), resol_id);
  // erase_resol_from_cache(get_ptr_cache<Type, Direction>(), resol_id);
}

template <class Type, fftType Direction>
std::vector<fft_plan<Type, Direction>> plan_all(int resol_id, int kfield, int *loens,
                                                   int nfft, int64_t *offsets) {
  using cmplx = typename Type::cmplx;
  using plan_type = typename Type::plan_type;
  int flags = FFTW_ESTIMATE | FFTW_NO_SIMD;

  auto key = cache_key{resol_id, kfield};
  auto &fftPlansCache = get_fft_plan_cache<Type, Direction>();
  auto fftPlans = fftPlansCache.find(key);
  if (fftPlans == fftPlansCache.end()) {
    // the fft plans do not exist yet
    cmplx *dummy;
    if constexpr (Direction == FFT_R2C || Direction == FFT_C2R)
      dummy = fftwf_alloc_complex((size_t)1);
    else
      dummy = fftw_alloc_complex((size_t)1);

    std::vector<fft_plan<Type, Direction>> newPlans;
    newPlans.reserve(nfft);
    for (int i = 0; i < nfft; ++i) {
      int nloen = loens[i];

      plan_type plan;
      int dist = offsets[i + 1] - offsets[i];
      int embed[] = {1};

      if constexpr (Direction == FFT_R2C)
        plan = fftwf_plan_many_dft_r2c(
          1, &nloen, kfield, (float*)dummy, embed, 1, dist, dummy, embed, 1, dist / 2, flags
        );
      else if constexpr (Direction == FFT_C2R)
        plan = fftwf_plan_many_dft_c2r(
          1, &nloen, kfield, dummy, embed, 1, dist / 2, (float*)dummy, embed, 1, dist, flags
        );
      else if constexpr (Direction == FFT_D2Z)
        plan = fftw_plan_many_dft_r2c(
          1, &nloen, kfield, (double*)dummy, embed, 1, dist, dummy, embed, 1, dist / 2, flags
        );
      else if constexpr (Direction == FFT_Z2D)
        plan = fftw_plan_many_dft_c2r(
          1, &nloen, kfield, dummy, embed, 1, dist / 2, (double*)dummy, embed, 1, dist, flags
        );
      newPlans.emplace_back(plan, kfield * offsets[i]);
    }
    fftPlansCache.insert({key, newPlans});
  }
  return fftPlansCache.find(key)->second;
}
} // namespace

extern "C" {
void execute_dir_fft_float(float *data_real, fftwf_complex *data_complex,
                           int resol_id, int kfield, int *loens, int64_t *offsets, int nfft,
                           void *growing_allocator) {
  auto plans = plan_all<Float, FFT_R2C>(resol_id, kfield, loens, nfft, offsets);
  for (auto &plan : plans)
    plan.exec(data_real, data_complex);
}
void execute_inv_fft_float(fftwf_complex *data_complex, float *data_real,
                           int resol_id, int kfield, int *loens, int64_t *offsets, int nfft,
                           void *growing_allocator) {
  auto plans = plan_all<Float, FFT_C2R>(resol_id, kfield, loens, nfft, offsets);
  for (auto &plan : plans)
    plan.exec(data_real, data_complex);
}
void execute_dir_fft_double(double *data_real,
                            fftw_complex *data_complex, int resol_id, int kfield,
                            int *loens, int64_t *offsets, int nfft,
                            void *growing_allocator) {
  auto plans = plan_all<Double, FFT_D2Z>(resol_id, kfield, loens, nfft, offsets);
  for (auto &plan : plans)
    plan.exec(data_real, data_complex);
}
void execute_inv_fft_double(fftw_complex *data_complex,
                            double *data_real, int resol_id, int kfield, int *loens,
                            int64_t *offsets, int nfft, void *growing_allocator) {
  auto plans = plan_all<Double, FFT_Z2D>(resol_id, kfield, loens, nfft, offsets);
  for (auto &plan : plans)
    plan.exec(data_real, data_complex);
}

void clean_fft(int resol_id) {
  erase_from_caches<Float, FFT_R2C>(resol_id);
  erase_from_caches<Float, FFT_C2R>(resol_id);
  erase_from_caches<Double, FFT_D2Z>(resol_id);
  erase_from_caches<Double, FFT_Z2D>(resol_id);
}
}
