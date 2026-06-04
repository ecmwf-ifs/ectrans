#include "hicfft_rocfft.h"

#include <algorithm>
#include <memory>
#include <unordered_map>
#include <vector>

#include "growing_allocator.h"

#define fftSafeCall(err) __fftSafeCall(err, __FILE__, __LINE__)

// Lazy initialization of rocfft library
static void ensure_rocfft_initialized() {
  static bool initialized = false;
  if (!initialized) {
    fftSafeCall(rocfft_setup());
    initialized = true;
  }
}

namespace {
struct Double {
  using real = double;
  using cmplx = double2;
};
struct Float {
  using real = float;
  using cmplx = float2;
};

template <class Type, rocfft_transform_type Direction> class hicfft_plan {
  using real = typename Type::real;
  using cmplx = typename Type::cmplx;

public:
  void exec(real *data_real, cmplx *data_complex) const {
    real *data_real_l = &data_real[offset];
    cmplx *data_complex_l = &data_complex[offset / 2];

    void* in_buffer[1];
    void* out_buffer[1];

    // Out-of-place transforms: separate input and output buffers
    if constexpr (Direction == rocfft_transform_type_real_forward) {
      // R2C: input real, output complex
      in_buffer[0] = (void*)data_real_l;
      out_buffer[0] = (void*)data_complex_l;
    } else if constexpr (Direction == rocfft_transform_type_real_inverse) {
      // C2R: input complex, output real
      in_buffer[0] = (void*)data_complex_l;
      out_buffer[0] = (void*)data_real_l;
    }

    fftSafeCall(rocfft_execute(*plan_ptr, in_buffer, out_buffer, *exec_info_ptr));
  }

  void set_stream(hipStream_t stream) {
    fftSafeCall(rocfft_execution_info_set_stream(*exec_info_ptr, stream));
  }

  hicfft_plan(rocfft_plan plan_, rocfft_execution_info exec_info_, int64_t offset_)
      : plan_ptr(new rocfft_plan{plan_},
                   [](auto ptr) {
                     fftSafeCall(rocfft_plan_destroy(*ptr));
                     delete ptr;
                   }),
        exec_info_ptr(new rocfft_execution_info{exec_info_},
                      [](auto ptr) {
                        fftSafeCall(rocfft_execution_info_destroy(*ptr));
                        delete ptr;
                      }),
        offset(offset_) {}

private:
  std::shared_ptr<rocfft_plan> plan_ptr;
  std::shared_ptr<rocfft_execution_info> exec_info_ptr;
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
template <class Type, rocfft_transform_type Direction> auto &get_fft_plan_cache() {
  static std::unordered_map<cache_key,
                            std::vector<hicfft_plan<Type, Direction>>>
      fftPlansCache;
  return fftPlansCache;
}
// kfield -> graphs
template <class Type, rocfft_transform_type Direction> auto &get_graph_cache() {
  static std::unordered_map<cache_key, std::shared_ptr<hipGraphExec_t>>
      graphCache;
  return graphCache;
}
// kfield -> ptrs
template <class Type, rocfft_transform_type Direction> auto &get_ptr_cache() {
  using real = typename Type::real;
  using cmplx = typename Type::cmplx;
  static std::unordered_map<cache_key, std::pair<real *, cmplx *>> ptrCache;
  return ptrCache;
}
template <class Type, rocfft_transform_type Direction>
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

template <class Type, rocfft_transform_type Direction>
void erase_from_caches(int resol_id) {
  erase_resol_from_cache(get_fft_plan_cache<Type, Direction>(), resol_id);
  erase_resol_from_cache(get_graph_cache<Type, Direction>(), resol_id);
  erase_resol_from_cache(get_ptr_cache<Type, Direction>(), resol_id);
}

template <class Type, rocfft_transform_type Direction>
std::vector<hicfft_plan<Type, Direction>> plan_all(int resol_id, int kfield, int *loens,
                                                   int nfft, int64_t *offsets) {
  ensure_rocfft_initialized();

  static constexpr bool is_forward =
      Direction == rocfft_transform_type_real_forward;

  auto key = cache_key{resol_id, kfield};
  auto &fftPlansCache = get_fft_plan_cache<Type, Direction>();
  auto fftPlans = fftPlansCache.find(key);

  if (fftPlans == fftPlansCache.end()) {
    // the fft plans do not exist yet
    std::vector<hicfft_plan<Type, Direction>> newPlans;
    newPlans.reserve(nfft);

    for (int i = 0; i < nfft; ++i) {
      int nloen = loens[i];

      rocfft_plan plan;
      rocfft_plan_description desc = nullptr;

      // Setup plan description for batched FFT
      fftSafeCall(rocfft_plan_description_create(&desc));

      int dist = offsets[i + 1] - offsets[i];
      size_t istride = 1;
      size_t ostride = 1;
      size_t idist = is_forward ? dist : dist / 2;
      size_t odist = is_forward ? dist / 2 : dist;

      fftSafeCall(rocfft_plan_description_set_data_layout(
          desc,
          is_forward ? rocfft_array_type_real : rocfft_array_type_hermitian_interleaved,
          is_forward ? rocfft_array_type_hermitian_interleaved : rocfft_array_type_real,
          nullptr, nullptr,  // offsets
          1, &istride, idist,
          1, &ostride, odist));

      size_t lengths[1] = {(size_t)nloen};

      fftSafeCall(rocfft_plan_create(
          &plan,
          rocfft_placement_notinplace,
          Direction,
          sizeof(typename Type::real) == 8 ? rocfft_precision_double : rocfft_precision_single,
          1,  // dimensions
          lengths,
          kfield,  // number of transforms (batch)
          desc));

      fftSafeCall(rocfft_plan_description_destroy(desc));

      // Create execution info for this plan
      rocfft_execution_info exec_info;
      fftSafeCall(rocfft_execution_info_create(&exec_info));

      // Query work buffer size and allocate if needed
      size_t work_buf_size = 0;
      fftSafeCall(rocfft_plan_get_work_buffer_size(plan, &work_buf_size));
      if (work_buf_size > 0) {
        void* work_buf = nullptr;
        HIC_CHECK(hipMalloc(&work_buf, work_buf_size));
        fftSafeCall(rocfft_execution_info_set_work_buffer(exec_info, work_buf, work_buf_size));
      }

      newPlans.emplace_back(plan, exec_info, kfield * offsets[i]);
    }
    fftPlansCache.insert({key, newPlans});
  }
  return fftPlansCache.find(key)->second;
}

template <class Type, rocfft_transform_type Direction>
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
    HIC_CHECK(hipStreamBeginCapture(stream, hipStreamCaptureModeGlobal));
    for (auto &plan : plans) {
      plan.exec(data_real, data_complex);
    }
    hipGraph_t my_graph;
    HIC_CHECK(hipStreamEndCapture(stream, &my_graph));
    hipGraphExec_t instance;
    HIC_CHECK(hipGraphInstantiate(&instance, my_graph, NULL, NULL, 0));
    HIC_CHECK(hipStreamDestroy(stream));

    graphCache.insert({key, std::shared_ptr<hipGraphExec_t>(
                                new hipGraphExec_t{instance}, [](auto ptr) {
                                  HIC_CHECK(hipGraphExecDestroy(*ptr));
                                  delete ptr;
                                })});
    ptrCache.insert({key, std::make_pair(data_real, data_complex)});
  }

  /* running in stream 0 */
  HIC_CHECK(hipGraphLaunch(*graphCache.at(key), 0));
  HIC_CHECK(hipStreamSynchronize(0));
}

template <class Type, rocfft_transform_type Direction>
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

void execute_dir_rocfft_float(float *data_real, float2 *data_complex,
                               int resol_id, int kfield, int *loens, int64_t *offsets, int nfft,
                               void *growing_allocator) {
  RUN<Float, rocfft_transform_type_real_forward>(data_real, data_complex, resol_id, kfield, loens, offsets, nfft,
                         growing_allocator);
}

void execute_inv_rocfft_float(float2 *data_complex, float *data_real,
                               int resol_id, int kfield, int *loens, int64_t *offsets, int nfft,
                               void *growing_allocator) {
  RUN<Float, rocfft_transform_type_real_inverse>(data_real, data_complex, resol_id, kfield, loens, offsets, nfft,
                         growing_allocator);
}

void execute_dir_rocfft_double(double *data_real,
                                double2 *data_complex, int resol_id, int kfield,
                                int *loens, int64_t *offsets, int nfft,
                                void *growing_allocator) {
  RUN<Double, rocfft_transform_type_real_forward>(data_real, data_complex, resol_id, kfield, loens,
                          offsets, nfft, growing_allocator);
}

void execute_inv_rocfft_double(double2 *data_complex,
                                double *data_real, int resol_id, int kfield, int *loens,
                                int64_t *offsets, int nfft, void *growing_allocator) {
  RUN<Double, rocfft_transform_type_real_inverse>(data_real, data_complex, resol_id, kfield, loens, offsets, nfft,
                          growing_allocator);
}
#undef RUN

void clean_rocfft(int resol_id) {
  erase_from_caches<Float, rocfft_transform_type_real_forward>(resol_id);
  erase_from_caches<Float, rocfft_transform_type_real_inverse>(resol_id);
  erase_from_caches<Double, rocfft_transform_type_real_forward>(resol_id);
  erase_from_caches<Double, rocfft_transform_type_real_inverse>(resol_id);
}
}

