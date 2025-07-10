#include "hicfft.h"

#include <algorithm>
#include <memory>

#include "growing_allocator.h"

#define fftSafeCall(err) __fftSafeCall(err, __FILE__, __LINE__)

// __global__ void debug(int varId, int N, HIP_DATA_TYPE_COMPLEX *x) {
//     for (int i = 0; i < N; i++)
//     {
//         HIP_DATA_TYPE_COMPLEX a = x[i];
//         double b = (double)a.x;
//         double c = (double)a.y;
//         if (varId == 0) printf("GPU: input[%d]=(%2.4f,%2.4f)\n",i+1,b,c);
//         if (varId == 1) printf("GPU: output[%d]=(%2.4f,%2.4f)\n",i+1,b,c);
//     }
// }

// __global__ void debugFloat(int varId, int N, HIP_DATA_TYPE_REAL *x) {
//     for (int i = 0; i < N; i++)
//     {
//         double a = (double)x[i];
//         if (varId == 0) printf("GPU: input[%d]=%2.4f\n",i+1,a);
//         if (varId == 1) printf("GPU: output[%d]=%2.4f\n",i+1,a);
//     }
// }

namespace {
struct Double {
  using real = double;
  using cmplx = hipfftDoubleComplex;
};
struct Float {
  using real = float;
  using cmplx = hipfftComplex;
};

template <class Type, hipfftType Direction> class hicfft_plan {
  using real = typename Type::real;
  using cmplx = typename Type::cmplx;

public:
  void exec(real *data_real, cmplx *data_complex) const {
    real *data_real_l = &data_real[offset];
    cmplx *data_complex_l = &data_complex[offset / 2];
    if constexpr (Direction == HIPFFT_R2C)
      fftSafeCall(hipfftExecR2C(*handle_ptr, data_real_l, data_complex_l));
    else if constexpr (Direction == HIPFFT_C2R)
      fftSafeCall(hipfftExecC2R(*handle_ptr, data_complex_l, data_real_l));
    else if constexpr (Direction == HIPFFT_D2Z)
      fftSafeCall(hipfftExecD2Z(*handle_ptr, data_real_l, data_complex_l));
    else if constexpr (Direction == HIPFFT_Z2D)
      fftSafeCall(hipfftExecZ2D(*handle_ptr, data_complex_l, data_real_l));
  }
  void set_stream(hipStream_t stream) {
    fftSafeCall(hipfftSetStream(*handle_ptr, stream));
  }
  hicfft_plan(hipfftHandle handle_, int64_t offset_)
      : handle_ptr(new hipfftHandle{handle_},
                   [](auto ptr) {
                     fftSafeCall(hipfftDestroy(*ptr));
                     delete ptr;
                   }),
        offset(offset_) {}

private:
  std::shared_ptr<hipfftHandle> handle_ptr;
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
template <class Type, hipfftType Direction> auto &get_fft_plan_cache() {
  static std::unordered_map<cache_key,
                            std::vector<hicfft_plan<Type, Direction>>>
      fftPlansCache;
  return fftPlansCache;
}
// kfield -> graphs
template <class Type, hipfftType Direction> auto &get_graph_cache() {
  static std::unordered_map<cache_key, std::shared_ptr<hipGraphExec_t>>
      graphCache;
  return graphCache;
}
// kfield -> ptrs
template <class Type, hipfftType Direction> auto &get_ptr_cache() {
  using real = typename Type::real;
  using cmplx = typename Type::cmplx;
  static std::unordered_map<cache_key, std::pair<real *, cmplx *>> ptrCache;
  return ptrCache;
}

template <class Type, hipfftType Direction>
void free_fft_graph_cache(float *, size_t) {
  get_graph_cache<Type, Direction>().clear();
  get_ptr_cache<Type, Direction>().clear();
}
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
template <class Type, hipfftType Direction>
void erase_from_caches(int resol_id) {
  erase_resol_from_cache(get_fft_plan_cache<Type, Direction>(), resol_id);
  erase_resol_from_cache(get_graph_cache<Type, Direction>(), resol_id);
  erase_resol_from_cache(get_ptr_cache<Type, Direction>(), resol_id);
}

template <class Type, hipfftType Direction>
std::vector<hicfft_plan<Type, Direction>> plan_all(int resol_id, int kfield, int *loens,
                                                   int nfft, int64_t *offsets) {
  static constexpr bool is_forward =
      Direction == HIPFFT_R2C || Direction == HIPFFT_D2Z;

  auto key = cache_key{resol_id, kfield};
  auto &fftPlansCache = get_fft_plan_cache<Type, Direction>();
  auto fftPlans = fftPlansCache.find(key);
  if (fftPlans == fftPlansCache.end()) {
    // the fft plans do not exist yet
    std::vector<hicfft_plan<Type, Direction>> newPlans;
    newPlans.reserve(nfft);
    for (int i = 0; i < nfft; ++i) {
      int nloen = loens[i];

      hipfftHandle plan;
      int dist = offsets[i + 1] - offsets[i];
      int embed[] = {1};
      fftSafeCall(hipfftPlanMany(
          &plan, 1, &nloen, embed, 1, is_forward ? dist : dist / 2, embed, 1,
          is_forward ? dist / 2 : dist, Direction, abs(kfield)));
      newPlans.emplace_back(plan, abs(kfield) * offsets[i]);
    }
    fftPlansCache.insert({key, newPlans});
  }
  return fftPlansCache.find(key)->second;
}

template <class Type, hipfftType Direction>
void run_group_graph(typename Type::real *data_real,
                     typename Type::cmplx *data_complex, int resol_id, int kfield, int *loens,
                     int64_t *offsets, int nfft, void *growing_allocator) {

  growing_allocator_register_free_c(growing_allocator,
                                    free_fft_graph_cache<Type, Direction>);

  // if the pointers are changed, we need to update the graph
  auto &ptrCache = get_ptr_cache<Type, Direction>();     // kfield -> ptrs
  auto &graphCache = get_graph_cache<Type, Direction>(); // kfield -> graphs

  auto key = cache_key{resol_id, kfield};
  auto ptrs = ptrCache.find(key);
  if (ptrs != ptrCache.end() && (ptrs->second.first != data_real ||
                                 ptrs->second.second != data_complex)) {
    // the plan is cached, but the pointers are not correct. we remove and
    // delete the graph, but we keep the FFT plans, if this happens more often,
    // we should cache this...
    std::cout << "WARNING FFT: POINTER CHANGE --> THIS MIGHT BE SLOW"
              << std::endl;
    graphCache.erase(key);
    ptrCache.erase(key);
  }

  auto graph = graphCache.find(key);
  if (graph == graphCache.end()) {
    // this graph does not exist yet
    auto plans =
        plan_all<Type, Direction>(resol_id, kfield, loens, nfft, offsets);

    // create a temporary stream
    hipStream_t stream;
    HIC_CHECK(hipStreamCreate(&stream));

    for (auto &plan : plans) // set the streams
      plan.set_stream(stream);

    // now create the graph
    hipGraph_t new_graph;
    hipGraphCreate(&new_graph, 0);
    for (auto &plan : plans) {
      HIC_CHECK(hipStreamBeginCapture(stream, hipStreamCaptureModeGlobal));
      plan.exec(data_real, data_complex);
      hipGraph_t my_graph;
      HIC_CHECK(hipStreamEndCapture(stream, &my_graph));
      hipGraphNode_t my_node;
      HIC_CHECK(
          hipGraphAddChildGraphNode(&my_node, new_graph, nullptr, 0, my_graph));
    }
    hipGraphExec_t instance;
    HIC_CHECK(hipGraphInstantiate(&instance, new_graph, NULL, NULL, 0));
    HIC_CHECK(hipStreamDestroy(stream));
    HIC_CHECK(hipGraphDestroy(new_graph));

    graphCache.insert({key, std::shared_ptr<hipGraphExec_t>(
                                new hipGraphExec_t{instance}, [](auto ptr) {
                                  HIC_CHECK(hipGraphExecDestroy(*ptr));
                                  delete ptr;
                                })});
    ptrCache.insert({key, std::make_pair(data_real, data_complex)});
  }

  HIC_CHECK(hipGraphLaunch(*graphCache.at(key), 0));
  HIC_CHECK(hipDeviceSynchronize());
}

template <class Type, hipfftType Direction>
void run_group(typename Type::real *data_real,
               typename Type::cmplx *data_complex, int resol_id, int kfield, int *loens,
               int64_t *offsets, int nfft, void *growing_allocator) {
  auto plans = plan_all<Type, Direction>(resol_id, kfield, loens, nfft, offsets);

  for (auto &plan : plans)
    plan.exec(data_real, data_complex);
  HIC_CHECK(hipDeviceSynchronize());
}
} // namespace

extern "C" {
#ifdef USE_GRAPHS_FFT
#define RUN run_group_graph
#else
#define RUN run_group
#endif
void execute_dir_fft_float(float *data_real, hipfftComplex *data_complex,
                           int resol_id, int kfield, int *loens, int64_t *offsets, int nfft,
                           void *growing_allocator) {
  RUN<Float, HIPFFT_R2C>(data_real, data_complex, resol_id, kfield, loens, offsets, nfft,
                         growing_allocator);
}
void execute_inv_fft_float(hipfftComplex *data_complex, float *data_real,
                           int resol_id, int kfield, int *loens, int64_t *offsets, int nfft,
                           void *growing_allocator) {
  RUN<Float, HIPFFT_C2R>(data_real, data_complex, resol_id, kfield, loens, offsets, nfft,
                         growing_allocator);
}
void execute_dir_fft_double(double *data_real,
                            hipfftDoubleComplex *data_complex, int resol_id, int kfield,
                            int *loens, int64_t *offsets, int nfft,
                            void *growing_allocator) {
  RUN<Double, HIPFFT_D2Z>(data_real, data_complex, resol_id, kfield, loens,
                          offsets, nfft, growing_allocator);
}
void execute_inv_fft_double(hipfftDoubleComplex *data_complex,
                            double *data_real, int resol_id, int kfield, int *loens,
                            int64_t *offsets, int nfft, void *growing_allocator) {
  RUN<Double, HIPFFT_Z2D>(data_real, data_complex, resol_id, kfield, loens, offsets, nfft,
                          growing_allocator);
}
#undef RUN

void clean_fft(int resol_id) {
  erase_from_caches<Float, HIPFFT_R2C>(resol_id);
  erase_from_caches<Float, HIPFFT_C2R>(resol_id);
  erase_from_caches<Double, HIPFFT_D2Z>(resol_id);
  erase_from_caches<Double, HIPFFT_Z2D>(resol_id);
}
}
