#include "hicfft.h"

#define fftSafeCall(err) __fftSafeCall(err, __FILE__, __LINE__)

#ifdef TRANS_SINGLE
typedef float DATA_TYPE;
typedef hipfftComplex HIP_DATA_TYPE_COMPLEX;
typedef hipfftReal HIP_DATA_TYPE_REAL;
#define fftExecDir hipfftExecR2C
#define fftExecInv hipfftExecC2R
#else
typedef double DATA_TYPE;
typedef hipfftDoubleComplex HIP_DATA_TYPE_COMPLEX;
typedef hipfftDoubleReal HIP_DATA_TYPE_REAL;
#define fftExecDir hipfftExecD2Z
#define fftExecInv hipfftExecZ2D
#endif

__global__ void debug(int varId, int N, HIP_DATA_TYPE_COMPLEX *x) {
    for (int i = 0; i < N; i++)
    {
        HIP_DATA_TYPE_COMPLEX a = x[i];
        double b = (double)a.x;
        double c = (double)a.y;
        if (varId == 0) printf("GPU: input[%d]=(%2.4f,%2.4f)\n",i+1,b,c);
        if (varId == 1) printf("GPU: output[%d]=(%2.4f,%2.4f)\n",i+1,b,c);
    }
}

__global__ void debugFloat(int varId, int N, HIP_DATA_TYPE_REAL *x) {
    for (int i = 0; i < N; i++)
    {
        double a = (double)x[i];
        if (varId == 0) printf("GPU: input[%d]=%2.4f\n",i+1,a);
        if (varId == 1) printf("GPU: output[%d]=%2.4f\n",i+1,a);
    }
}

extern "C"
void
hicfft_execute_plan_(int ISIGNp, int N, DATA_TYPE *data_in_host, DATA_TYPE *data_out_host, long *iplan)
{
    HIP_DATA_TYPE_COMPLEX *data_in = reinterpret_cast<HIP_DATA_TYPE_COMPLEX*>(data_in_host);
    HIP_DATA_TYPE_COMPLEX *data_out = reinterpret_cast<HIP_DATA_TYPE_COMPLEX*>(data_out_host);
    hipfftHandle* PLANp = reinterpret_cast<hipfftHandle*>(iplan);
    hipfftHandle plan = *PLANp;
    int ISIGN = ISIGNp;

    /*if (hipDeviceSynchronize() != hipSuccess){
        fprintf(stderr, "GPU runtime error: Failed to synchronize\n");
        return;
    }*/

    if( ISIGN== -1 ){
        fftSafeCall(fftExecDir(plan, (HIP_DATA_TYPE_REAL*)data_in, data_out));
    }
    else if( ISIGN== 1){
        fftSafeCall(fftExecInv(plan, data_in, (HIP_DATA_TYPE_REAL*)data_out));
    }
    else {
        abort();
    }

    if (hipDeviceSynchronize() != hipSuccess){
        fprintf(stderr, "GPU runtime error: Failed to synchronize\n");
        return;
    }

}

namespace {
struct Double {
  using real = double;
  using cmplx = hipfftDoubleComplex;
};
struct Float {
  using real = float;
  using cmplx = hipfftComplex;
};
}
template <class Type, hipfftType Direction>
void execute_fft(typename Type::real *data_real, typename Type::cmplx *data_complex,
        int kfield, int *loens, int *offsets, int nfft) {
  using real = typename Type::real;
  using cmplx = typename Type::cmplx;

  /* static std::unordered_map<int, void *> allocationCache; // nloens -> ptr */
  static std::unordered_map<int, std::vector<hipfftHandle>> fftPlansCache; // kfield -> handles
  static std::unordered_map<int, hipGraphExec_t> graphCache; // kfield -> graphs

  // if the pointers are changed, we need to update the graph
  static std::unordered_map<int, std::pair<real *, cmplx *>> ptrCache; // kfield -> ptrs

  auto ptrs = ptrCache.find(kfield);
  if (ptrs != ptrCache.end() && (
              ptrs->second.first != data_real || ptrs->second.second != data_complex)) {
      // the plan is cached, but the pointers are not correct. we remove and delete the graph,
      // but we keep the FFT plans, if this happens more often, we should cache this...
      std::cout << "WARNING: POINTER CHANGE --> THIS MIGHT BE SLOW" << std::endl;
      HIC_CHECK(hipGraphExecDestroy(graphCache[kfield]));
      graphCache.erase(kfield);
      ptrCache.erase(kfield);
  }

  auto graph = graphCache.find(kfield);
  if (graph == graphCache.end()) {
      // this graph does not exist yet

      auto fftPlans = fftPlansCache.find(kfield);
      if (fftPlans == fftPlansCache.end()) {
          // the fft plans do not exist yet
          std::vector<hipfftHandle> newPlans;
          newPlans.resize(nfft);
          for (int i = 0; i < nfft; ++i) {
            int nloen = loens[i];

            hipfftHandle plan;
            fftSafeCall(hipfftCreate(&plan));
            int dist = 1;
            int embed[] = {1};
            fftSafeCall(hipfftPlanMany(&plan, 1, &nloen, embed, kfield, dist, embed,
                                      kfield, dist, Direction, kfield));
            newPlans[i] = plan;
          }
          fftPlansCache.insert({kfield, newPlans});
      }
      fftPlans = fftPlansCache.find(kfield);

      // create a temporary stream
      hipStream_t stream;
      HIC_CHECK(hipStreamCreate(&stream));

      for (auto &plan : fftPlans->second) // set the streams
        fftSafeCall(hipfftSetStream(plan, stream));

      // now create the graph
      hipGraph_t new_graph;
      hipGraphCreate(&new_graph, 0);
      for (int i = 0; i < nfft; ++i) {
        int offset = offsets[i];
        real *data_real_l = &data_real[kfield * offset];
        cmplx *data_complex_l = &data_complex[kfield * offset / 2];
        HIC_CHECK(hipStreamBeginCapture(stream, hipStreamCaptureModeGlobal));
        if constexpr(Direction == HIPFFT_R2C)
          fftSafeCall(hipfftExecR2C(fftPlans->second[i], data_real_l, data_complex_l));
        else if constexpr(Direction == HIPFFT_C2R)
          fftSafeCall(hipfftExecC2R(fftPlans->second[i], data_complex_l, data_real_l));
        else if constexpr(Direction == HIPFFT_D2Z)
          fftSafeCall(hipfftExecD2Z(fftPlans->second[i], data_real_l, data_complex_l));
        else if constexpr(Direction == HIPFFT_Z2D)
          fftSafeCall(hipfftExecZ2D(fftPlans->second[i], data_complex_l, data_real_l));
        hipGraph_t my_graph;
        HIC_CHECK(hipStreamEndCapture(stream, &my_graph));
        hipGraphNode_t my_node;
        HIC_CHECK(hipGraphAddChildGraphNode(&my_node, new_graph, nullptr, 0, my_graph));
      }
      hipGraphExec_t instance;
      HIC_CHECK(hipGraphInstantiate(&instance, new_graph, NULL, NULL, 0));
      HIC_CHECK(hipStreamDestroy(stream));
      HIC_CHECK(hipGraphDestroy(new_graph));

      graphCache.insert({kfield, instance});
      ptrCache.insert({kfield, std::make_pair(data_real, data_complex)});
  }

  HIC_CHECK(hipGraphLaunch(graphCache.at(kfield), 0));
      /* for (int i = 0; i < nfft; ++i) { */
        /* int nloen = loens[i]; */

        /* hipfftHandle plan; */
        /* fftSafeCall(hipfftCreate(&plan)); */
        /* int dist = 1; */
        /* int embed[] = {1}; */
        /* fftSafeCall(hipfftPlanMany(&plan, 1, &nloen, embed, kfield, dist, embed, */
                                  /* kfield, dist, Direction, kfield)); */
        /* int offset = offsets[i]; */
        /* real *data_real_l = &data_real[kfield * offset]; */
        /* cmplx *data_complex_l = &data_complex[kfield * offset / 2]; */
        /* if (Direction == HIPFFT_R2C) */
          /* fftSafeCall(hipfftExecR2C(plan, data_real_l, data_complex_l)) */
        /* else */
          /* fftSafeCall(hipfftExecC2R(plan, data_complex_l, data_real_l)); */
        /* fftSafeCall(hipfftDestroy(plan)); */
      /* } */
  HIC_CHECK(hipDeviceSynchronize());
}

template <class Type, hipfftType Direction>
void execute_fft_new(typename Type::real *data_real, typename Type::cmplx *data_complex,
        int kfield, int *loens, int *offsets, int nfft) {
  using real = typename Type::real;
  using cmplx = typename Type::cmplx;

  /* static std::unordered_map<int, void *> allocationCache; // nloens -> ptr */
  static std::unordered_map<int, std::vector<hipfftHandle>> fftPlansCache; // kfield -> handles
  static std::unordered_map<int, hipGraphExec_t> graphCache; // kfield -> graphs

  // if the pointers are changed, we need to update the graph
  static std::unordered_map<int, std::pair<real *, cmplx *>> ptrCache; // kfield -> ptrs

  auto ptrs = ptrCache.find(kfield);
  if (ptrs != ptrCache.end() && (
              ptrs->second.first != data_real || ptrs->second.second != data_complex)) {
      // the plan is cached, but the pointers are not correct. we remove and delete the graph,
      // but we keep the FFT plans, if this happens more often, we should cache this...
      std::cout << "WARNING: POINTER CHANGE --> THIS MIGHT BE SLOW" << std::endl;
      HIC_CHECK(hipGraphExecDestroy(graphCache[kfield]));
      graphCache.erase(kfield);
      ptrCache.erase(kfield);
  }

  auto graph = graphCache.find(kfield);
  if (graph == graphCache.end()) {
      // this graph does not exist yet

      auto fftPlans = fftPlansCache.find(kfield);
      if (fftPlans == fftPlansCache.end()) {
          // the fft plans do not exist yet
          std::vector<hipfftHandle> newPlans;
          newPlans.resize(nfft);
          for (int i = 0; i < nfft; ++i) {
            int nloen = loens[i];

            hipfftHandle plan;
            fftSafeCall(hipfftCreate(&plan));
            int dist = offsets[i+1] - offsets[i];
            int embed[] = {1};
            fftSafeCall(hipfftPlanMany(&plan, 1, &nloen, embed, 1, dist, embed,
                                      1, dist / 2, Direction, kfield));
            newPlans[i] = plan;
          }
          fftPlansCache.insert({kfield, newPlans});
      }
      fftPlans = fftPlansCache.find(kfield);

      // create a temporary stream
      hipStream_t stream;
      HIC_CHECK(hipStreamCreate(&stream));

      for (auto &plan : fftPlans->second) // set the streams
        fftSafeCall(hipfftSetStream(plan, stream));

      // now create the graph
      hipGraph_t new_graph;
      hipGraphCreate(&new_graph, 0);
      for (int i = 0; i < nfft; ++i) {
        int offset = offsets[i];
        real *data_real_l = &data_real[kfield * offset];
        cmplx *data_complex_l = &data_complex[kfield * offset / 2];
        HIC_CHECK(hipStreamBeginCapture(stream, hipStreamCaptureModeGlobal));
        if constexpr(Direction == HIPFFT_R2C)
          fftSafeCall(hipfftExecR2C(fftPlans->second[i], data_real_l, data_complex_l));
        else if constexpr(Direction == HIPFFT_C2R)
          fftSafeCall(hipfftExecC2R(fftPlans->second[i], data_complex_l, data_real_l));
        else if constexpr(Direction == HIPFFT_D2Z)
          fftSafeCall(hipfftExecD2Z(fftPlans->second[i], data_real_l, data_complex_l));
        else if constexpr(Direction == HIPFFT_Z2D)
          fftSafeCall(hipfftExecZ2D(fftPlans->second[i], data_complex_l, data_real_l));
        hipGraph_t my_graph;
        HIC_CHECK(hipStreamEndCapture(stream, &my_graph));
        hipGraphNode_t my_node;
        HIC_CHECK(hipGraphAddChildGraphNode(&my_node, new_graph, nullptr, 0, my_graph));
      }
      hipGraphExec_t instance;
      HIC_CHECK(hipGraphInstantiate(&instance, new_graph, NULL, NULL, 0));
      HIC_CHECK(hipStreamDestroy(stream));
      HIC_CHECK(hipGraphDestroy(new_graph));

      graphCache.insert({kfield, instance});
      ptrCache.insert({kfield, std::make_pair(data_real, data_complex)});
  }

  HIC_CHECK(hipGraphLaunch(graphCache.at(kfield), 0));
      /* for (int i = 0; i < nfft; ++i) { */
        /* int nloen = loens[i]; */

        /* hipfftHandle plan; */
        /* HIPFFT_CHECK(hipfftCreate(&plan)); */
        /* int dist = 1; */
        /* int embed[] = {1}; */
        /* HIPFFT_CHECK(hipfftPlanMany(&plan, 1, &nloen, embed, kfield, dist, embed, */
                                  /* kfield, dist, Direction, kfield)); */
        /* int offset = offsets[i]; */
        /* real *data_real_l = &data_real[kfield * offset]; */
        /* cmplx *data_complex_l = &data_complex[kfield * offset / 2]; */
        /* if (Direction == HIPFFT_R2C) */
          /* HIPFFT_CHECK(hipfftExecR2C(plan, data_real_l, data_complex_l)) */
        /* else */
          /* HIPFFT_CHECK(hipfftExecC2R(plan, data_complex_l, data_real_l)); */
        /* HIPFFT_CHECK(hipfftDestroy(plan)); */
      /* } */
  HIC_CHECK(hipDeviceSynchronize());
}


extern "C" {
void execute_dir_fft_float(float *data_real, hipfftComplex *data_complex,
        int kfield, int *loens, int *offsets, int nfft) {
    //execute_fft<Float, HIPFFT_R2C>(data_real, data_complex, kfield, loens, offsets, nfft);
    execute_fft_new<Float, HIPFFT_R2C>(data_real, data_complex, kfield, loens, offsets, nfft);
}
void execute_inv_fft_float(hipfftComplex *data_complex, float *data_real,
        int kfield, int *loens, int *offsets, int nfft) {
    execute_fft<Float, HIPFFT_C2R>(data_real, data_complex, kfield, loens, offsets, nfft);
}
void execute_dir_fft_double(double *data_real, hipfftDoubleComplex *data_complex,
        int kfield, int *loens, int *offsets, int nfft) {
    //execute_fft<Double, HIPFFT_D2Z>(data_real, data_complex, kfield, loens, offsets, nfft);
    execute_fft_new<Double, HIPFFT_D2Z>(data_real, data_complex, kfield, loens, offsets, nfft);
}
void execute_inv_fft_double(hipfftDoubleComplex *data_complex, double *data_real,
        int kfield, int *loens, int *offsets, int nfft) {
    execute_fft<Double, HIPFFT_Z2D>(data_real, data_complex, kfield, loens, offsets, nfft);
}
}

