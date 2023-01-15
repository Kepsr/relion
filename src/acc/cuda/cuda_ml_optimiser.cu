#undef ALTCPU
#include <sys/time.h>
#include <stdio.h>
#include <time.h>
#include <math.h>
#include <ctime>
#include <vector>
#include <iostream>
#include "src/ml_optimiser.h"
#include <cuda_runtime.h>
#include <curand.h>
#include <curand_kernel.h>

#include "src/acc/acc_ptr.h"
#include "src/acc/acc_projector.h"
#include "src/acc/acc_backprojector.h"
#include "src/acc/acc_backprojector_impl.h"
#include "src/acc/acc_projector_plan.h"
#include "src/acc/cuda/cuda_benchmark_utils.h"
#include "src/acc/cuda/cuda_kernels/helper.cuh"
#include "src/acc/cuda/cuda_kernels/diff2.cuh"
#include "src/acc/cuda/cuda_kernels/wavg.cuh"
#include "src/acc/cuda/cuda_mem_utils.h"
#include "src/acc/cuda/cuda_fft.h"
#include "src/acc/data_types.h"
#include "src/complex.h"
#include "src/helix.h"
#include "src/error.h"
#include <fstream>
#include "src/parallel.h"
#include <signal.h>
#include <map>

#ifdef CUDA_FORCESTL
#include "src/acc/cuda/cuda_utils_stl.cuh"
#else
#include "src/acc/cuda/cuda_utils_cub.cuh"
#endif

#include "src/acc/utilities.h"
#include "src/acc/utilities_impl.h"

#include "src/acc/acc_ml_optimiser.h"
#include "src/acc/cuda/cuda_ml_optimiser.h"
#include "src/acc/acc_helper_functions.h"
#include "src/acc/acc_ml_optimiser_impl.h"

// -------------------------------  Some explicit template instantiations
template __global__ void CudaKernels::cuda_kernel_translate2D<XFLOAT>(
    XFLOAT const*, XFLOAT*, int, int, int, int, int
);

template __global__ void CudaKernels::cuda_kernel_translate3D<XFLOAT>(
    XFLOAT const*, XFLOAT*, int, int, int, int, int, int, int
);

template __global__ void cuda_kernel_multi<XFLOAT>(
    XFLOAT const*, XFLOAT*, XFLOAT, int
);

template __global__ void CudaKernels::cuda_kernel_multi<XFLOAT>(
    XFLOAT*, XFLOAT, int
);

template __global__ void cuda_kernel_multi<XFLOAT>(
    XFLOAT const*, XFLOAT const*, XFLOAT*, XFLOAT, int
);

// ----------------------------------------------------------------------

// High-level CUDA objects

size_t MlDeviceBundle::checkFixedSizedObjects(int shares) {
    int devCount;
    HANDLE_ERROR(cudaGetDeviceCount(&devCount));
    if (device_id >= devCount) {
        CRITICAL(ERR_GPUID);
    }

    HANDLE_ERROR(cudaSetDevice(device_id));

    size_t free(0), total(0);
    DEBUG_HANDLE_ERROR(cudaMemGetInfo(&free, &total));
    float margin(1.05);
    const size_t BoxLimit = pow(free / (margin * 2.5 * sizeof(XFLOAT) * ((float) shares)),(1 / 3.0)) / ((float) baseMLO->mymodel.padding_factor);
    // size_t BytesNeeded = ((float) shares) * margin * 2.5 * sizeof(XFLOAT) * pow((baseMLO->mymodel.ori_size * baseMLO->mymodel.padding_factor),3);
    return BoxLimit;
}

void MlDeviceBundle::setupFixedSizedObjects() {
    int devCount;
    HANDLE_ERROR(cudaGetDeviceCount(&devCount));
    if (device_id >= devCount) {
        //std::cerr << " using device_id=" << device_id << " (device no. " << device_id+1 << ") which is higher than the available number of devices=" << devCount << std::endl;
        CRITICAL(ERR_GPUID);
    } else {
        HANDLE_ERROR(cudaSetDevice(device_id));
    }

    // Can we pre-generate projector plan and corresponding euler matrices for all particles
    generateProjectionPlanOnTheFly = baseMLO->do_skip_align || baseMLO->do_skip_rotate || baseMLO->do_auto_refine || baseMLO->mymodel.orientational_prior_mode != NOPRIOR;

    const unsigned nr_proj  = baseMLO->mymodel.PPref.size();
    const unsigned nr_bproj = baseMLO->wsum_model.BPref.size();

    projectors.resize(nr_proj);
    backprojectors.resize(nr_bproj);

    /*======================================================
                  PROJECTOR AND BACKPROJECTOR
    ======================================================*/

    for (int i = 0; i < nr_proj; i++) {
        const auto& pp = baseMLO->mymodel.PPref[i];
        projectors[i].setMdlDim(
            pp.data.xdim, pp.data.ydim, pp.data.zdim,
            pp.data.yinit, pp.data.zinit,
            pp.r_max, pp.padding_factor
        );
        projectors[i].initMdl(pp.data.data);
    }

    for (int i = 0; i < nr_bproj; i++) {
        const auto& bp = baseMLO->wsum_model.BPref[i];
        backprojectors[i].setMdlDim(
            bp.data.xdim, bp.data.ydim, bp.data.zdim,
            bp.data.yinit, bp.data.zinit,
            bp.r_max, bp.padding_factor
        );
        backprojectors[i].initMdl();
    }

    /*======================================================
                        CUSTOM ALLOCATOR
    ======================================================*/

    int memAlignmentSize;
    cudaDeviceGetAttribute(&memAlignmentSize, cudaDevAttrTextureAlignment, device_id);
    allocator = new CudaCustomAllocator(0, memAlignmentSize);
}

void MlDeviceBundle::setupTunableSizedObjects(size_t allocationSize) {
    int devCount;
    HANDLE_ERROR(cudaGetDeviceCount(&devCount));
    if (device_id >= devCount) {
        //std::cerr << " using device_id=" << device_id << " (device no. " << device_id+1 << ") which is higher than the available number of devices=" << devCount << std::endl;
        CRITICAL(ERR_GPUID);
    } else {
        HANDLE_ERROR(cudaSetDevice(device_id));
    }

    /*======================================================
                        CUSTOM ALLOCATOR
    ======================================================*/
    #ifdef DEBUG_CUDA
    printf("DEBUG: Total GPU allocation size set to %zu MB on device id %d.\n", allocationSize / (1000*1000), device_id);
    #endif
    #ifndef CUDA_NO_CUSTOM_ALLOCATION
    allocator->resize(allocationSize);
    #endif

    /*======================================================
                        PROJECTION PLAN
    ======================================================*/

    const unsigned nr_models = baseMLO->mymodel.nr_classes;
    coarseProjectionPlans.resize(nr_models, allocator);

    if (generateProjectionPlanOnTheFly) return;
    // Otherwise, if doing predefined projector plan:
    for (int iclass = 0; iclass < nr_models; iclass++) {
        // Skip insignificant classes
        if (baseMLO->mymodel.pdf_class[iclass] <= 0.0) continue;
        std::vector<int> exp_pointer_dir_nonzeroprior;
        std::vector<int> exp_pointer_psi_nonzeroprior;
        std::vector<RFLOAT> exp_directions_prior;
        std::vector<RFLOAT> exp_psi_prior;

        const long unsigned itrans_max = baseMLO->sampling.NrTranslationalSamplings() - 1;
        const long unsigned nr_idir    = baseMLO->sampling.NrDirections();
        const long unsigned nr_ipsi    = baseMLO->sampling.NrPsiSamplings();

        coarseProjectionPlans[iclass].setup(
            baseMLO->sampling,
            exp_directions_prior,
            exp_psi_prior,
            exp_pointer_dir_nonzeroprior,
            exp_pointer_psi_nonzeroprior,
            nullptr,  // Mcoarse_significant
            baseMLO->mymodel.pdf_class,
            baseMLO->mymodel.pdf_direction,
            nr_idir, nr_ipsi,
            0 /*idir_min*/, nr_idir - 1 /*idir_max*/,
            0 /*ipsi_min*/, nr_ipsi - 1 /*ipsi_max*/,
            0 /*itrans_min*/, itrans_max,
            0 /*current_oversampling*/, 1 /*nr_oversampled_rot*/,
            iclass,
            true,  // coarse
            !IS_NOT_INV,
            baseMLO->do_skip_align,
            baseMLO->do_skip_rotate,
            baseMLO->mymodel.orientational_prior_mode
        );
    }
};

void MlOptimiserCuda::resetData() {
    int devCount;
    HANDLE_ERROR(cudaGetDeviceCount(&devCount));
    if (device_id >= devCount) {
        //std::cerr << " using device_id=" << device_id << " (device no. " << device_id+1 << ") which is higher than the available number of devices=" << devCount << std::endl;
        CRITICAL(ERR_GPUID);
    } else {
        HANDLE_ERROR(cudaSetDevice(device_id));
    }

    classStreams.resize(baseMLO->mymodel.nr_classes, nullptr);
    for (cudaStream_t& stream: classStreams) {
        HANDLE_ERROR(cudaStreamCreate(&stream));
        // HANDLE_ERROR(cudaStreamCreateWithFlags(&stream, cudaStreamNonBlocking));
    }

    transformer1.clear();
    transformer2.clear();
};

void MlOptimiserCuda::doThreadExpectationSomeParticles(int thread_id) {
    #ifdef TIMING
    // Only time one thread
    if (thread_id == 0)
    baseMLO->timer.tic(baseMLO->TIMING_ESP_THR);
    #endif
	// CTOC(cudaMLO->timer, "interParticle");

    int devCount;
    HANDLE_ERROR(cudaGetDeviceCount(&devCount));
    if (device_id >= devCount) {
        // std::cerr << " using device_id=" << device_id << " (device no. " << device_id+1 << ") which is higher than the available number of devices=" << devCount << std::endl;
        CRITICAL(ERR_GPUID);
    } else {
        DEBUG_HANDLE_ERROR(cudaSetDevice(device_id));
    }
    // std::cerr << " calling on device " << device_id << std::endl;
    // put mweight allocation here
    size_t first_ipart = 0, last_ipart = 0;

    while (baseMLO->exp_ipart_ThreadTaskDistributor->getTasks(first_ipart, last_ipart)) {
        CTICTOC(timer, "oneTask", ({
        for (long unsigned ipart = first_ipart; ipart <= last_ipart; ipart++) {
            #ifdef TIMING
            // Only time one thread
            if (thread_id == 0)
            baseMLO->timer.tic(baseMLO->TIMING_ESP_DIFF2_A);
            #endif

            AccPtrFactory<> ptrFactory (allocator, cudaStreamPerThread);
            accDoExpectationOneParticle<MlOptimiserCuda>(this, baseMLO->exp_my_first_part_id + ipart, thread_id, ptrFactory);

        }
        }))
    }

    // CTIC(cudaMLO->timer, "interParticle");
    // exit(0);

    #ifdef TIMING
    // Only time one thread
    if (thread_id == 0)
    baseMLO->timer.toc(baseMLO->TIMING_ESP_THR);
    #endif
}
