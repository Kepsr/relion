#undef ALTCPU
#include <cuda_runtime.h>
#include "src/ml_optimiser.h"
#include "src/acc/acc_ptr.h"
#include "src/acc/acc_projector.h"
#include "src/acc/acc_projector_plan.h"
#include "src/acc/acc_backprojector.h"
#include "src/acc/cuda/cuda_settings.h"
#include "src/acc/cuda/cuda_fft.h"
#include "src/acc/cuda/cuda_kernels/cuda_device_utils.cuh"

#ifdef CUDA_FORCESTL
#include "src/acc/cuda/cuda_utils_stl.cuh"
#else
#include "src/acc/cuda/cuda_utils_cub.cuh"
#endif

#include "src/acc/utilities.h"
#include "src/acc/acc_helper_functions.h"
#include "src/acc/cuda/cuda_kernels/BP.cuh"
#include "src/macros.h"
#include "src/error.h"

#include "src/acc/acc_ml_optimiser.h"
#include "src/acc/cuda/cuda_ml_optimiser.h"
#include "src/acc/acc_helper_functions.h"


#include "src/acc/acc_helper_functions_impl.h"

#include <cuda_runtime.h>
#include "src/acc/cuda/cuda_settings.h"
#include "src/acc/cuda/cuda_helper_functions.cuh"
#include "src/acc/cuda/cuda_kernels/BP.cuh"
#include "src/macros.h"
#include "src/error.h"

/*
 * Maps weights to a decoupled indexing of translations and orientations
 */
void mapWeights(
        unsigned long orientation_start,
        XFLOAT *mapped_weights,
        unsigned orientation_num,
        unsigned long idxArr_start,
        unsigned long idxArr_end,
        unsigned translation_num,
        XFLOAT *weights,
        long unsigned *rot_idx,
        long unsigned *trans_idx,
        unsigned long current_oversampling)
{

    for (long unsigned i = 0; i < orientation_num*translation_num; i++)
        mapped_weights[i] = -999.;

    for (long unsigned i = idxArr_start; i < idxArr_end; i++)
        mapped_weights[ (rot_idx[i]-orientation_start) * translation_num + trans_idx[i] ]= weights[i];
}

void generateEulerMatrices(
    XFLOAT padding_factor,
    ProjectionParams &ProjectionData,
    XFLOAT *eulers,
    bool inverse
) {

    for (long int i = 0; i < ProjectionData.rots.size(); i++) {
        // TODO In a sense we're doing degrees just to do radians here.
        // The only place the degree value is actually used is in the metadata assignment.

        RFLOAT alpha = radians(ProjectionData.rots[i]);
        RFLOAT beta  = radians(ProjectionData.tilts[i]);
        RFLOAT gamma = radians(ProjectionData.psis[i]);

        RFLOAT ca, sa, cb, sb, cg, sg;
        sincos(alpha, &sa, &ca);
        sincos(beta,  &sb, &cb);
        sincos(gamma, &sg, &cg);

        RFLOAT cc = cb * ca;
        RFLOAT cs = cb * sa;
        RFLOAT sc = sb * ca;
        RFLOAT ss = sb * sa;

        if (inverse) {
            eulers[9 * i + 0] = ( cg * cc - sg * sa) ;// * padding_factor; //00
            eulers[9 * i + 1] = (-sg * cc - cg * sa) ;// * padding_factor; //10
            eulers[9 * i + 2] = ( sc )               ;// * padding_factor; //20
            eulers[9 * i + 3] = ( cg * cs + sg * ca) ;// * padding_factor; //01
            eulers[9 * i + 4] = (-sg * cs + cg * ca) ;// * padding_factor; //11
            eulers[9 * i + 5] = ( ss )               ;// * padding_factor; //21
            eulers[9 * i + 6] = (-cg * sb )          ;// * padding_factor; //02
            eulers[9 * i + 7] = ( sg * sb )          ;// * padding_factor; //12
            eulers[9 * i + 8] = ( cb )               ;// * padding_factor; //22
        } else {
            eulers[9 * i + 0] = ( cg * cc - sg * sa) ;// * padding_factor; //00
            eulers[9 * i + 1] = ( cg * cs + sg * ca) ;// * padding_factor; //01
            eulers[9 * i + 2] = (-cg * sb )          ;// * padding_factor; //02
            eulers[9 * i + 3] = (-sg * cc - cg * sa) ;// * padding_factor; //10
            eulers[9 * i + 4] = (-sg * cs + cg * ca) ;// * padding_factor; //11
            eulers[9 * i + 5] = ( sg * sb )          ;// * padding_factor; //12
            eulers[9 * i + 6] = ( sc )               ;// * padding_factor; //20
            eulers[9 * i + 7] = ( ss )               ;// * padding_factor; //21
            eulers[9 * i + 8] = ( cb )               ;// * padding_factor; //22
        }
    }
}

__global__ void cuda_kernel_allweights_to_mweights(
    unsigned long *d_iorient,
    XFLOAT *d_allweights, XFLOAT *d_mweights,
    unsigned long orientation_num, unsigned long translation_num
) {
    size_t idx = blockIdx.x * WEIGHT_MAP_BLOCK_SIZE + threadIdx.x;
    if (idx < orientation_num * translation_num)
        d_mweights[d_iorient[idx / translation_num] * translation_num + idx % translation_num] =
                d_allweights[idx / translation_num  * translation_num + idx % translation_num];
}

size_t findThresholdIdxInCumulativeSum(CudaGlobalPtr<XFLOAT> &data, XFLOAT threshold) {
    int grid_size = ceil((float) (data.getSize() - 1) / (float) FIND_IN_CUMULATIVE_BLOCK_SIZE);
    if (grid_size == 0) return 0;
    CudaGlobalPtr<size_t> idx(1, data.getStream(), data.getAllocator());
    idx[0] = 0;
    idx.put_on_device();

    cuda_kernel_find_threshold_idx_in_cumulative
    <<<grid_size, FIND_IN_CUMULATIVE_BLOCK_SIZE, 0, data.getStream()>>>
    (~data, threshold, data.getSize() - 1, ~idx);

    idx.cp_to_host();
    DEBUG_HANDLE_ERROR(cudaStreamSynchronize(data.getStream()));

    return idx[0];
}

void runCollect2jobs(
    dim3 grid_dim,
    XFLOAT *oo_otrans_x,          // otrans-size -> make const
    XFLOAT *oo_otrans_y,          // otrans-size -> make const
    XFLOAT *oo_otrans_z,          // otrans-size -> make const
    XFLOAT *myp_oo_otrans_x2y2z2, // otrans-size -> make const
    XFLOAT *weights,
    XFLOAT significant_weight,
    XFLOAT sum_weight,
    unsigned long nr_trans,
    unsigned long nr_oversampled_trans,
    unsigned long nr_oversampled_rot,
    int oversamples, bool skip_rots,
    XFLOAT *p_weights,
    XFLOAT *p_thr_wsum_prior_offsetx_class,
    XFLOAT *p_thr_wsum_prior_offsety_class,
    XFLOAT *p_thr_wsum_prior_offsetz_class,
    XFLOAT *p_thr_wsum_sigma2_offset,
    size_t *rot_idx, size_t *trans_idx, size_t *jobOrigin, size_t *jobExtent,
    bool data_is_3D
) {
    if (data_is_3D) {
        size_t shared_buffer = sizeof(XFLOAT) * SUMW_BLOCK_SIZE * 5; // x+y+z+myp+weights
        cuda_kernel_collect2jobs<true><<<grid_dim, SUMW_BLOCK_SIZE, shared_buffer>>> (
            oo_otrans_x,          // otrans-size -> make const
            oo_otrans_y,          // otrans-size -> make const
            oo_otrans_z,          // otrans-size -> make const
            myp_oo_otrans_x2y2z2, // otrans-size -> make const
            weights, significant_weight, sum_weight,
            nr_trans, nr_oversampled_trans, nr_oversampled_rot,
            oversamples, skip_rots, p_weights,
            p_thr_wsum_prior_offsetx_class,
            p_thr_wsum_prior_offsety_class,
            p_thr_wsum_prior_offsetz_class,
            p_thr_wsum_sigma2_offset,
            rot_idx, trans_idx, jobOrigin, jobExtent
        );
    } else {
        size_t shared_buffer = sizeof(XFLOAT) * SUMW_BLOCK_SIZE * 4; // x+y+myp+weights
        cuda_kernel_collect2jobs<false><<<grid_dim, SUMW_BLOCK_SIZE, shared_buffer>>>
        (
            oo_otrans_x,          // otrans-size -> make const
            oo_otrans_y,          // otrans-size -> make const
            oo_otrans_z,          // otrans-size -> make const
            myp_oo_otrans_x2y2z2, // otrans-size -> make const
            weights, significant_weight, sum_weight,
            nr_trans, nr_oversampled_trans, nr_oversampled_rot,
            oversamples, skip_rots, p_weights,
            p_thr_wsum_prior_offsetx_class,
            p_thr_wsum_prior_offsety_class,
            p_thr_wsum_prior_offsetz_class,
            p_thr_wsum_sigma2_offset,
            rot_idx, trans_idx, jobOrigin, jobExtent
        );
    }
}

//  void windowFourierTransform2(
//      XFLOAT *d_in_real,  XFLOAT *d_in_imag,
//      XFLOAT *d_out_real, XFLOAT *d_out_imag,
//      unsigned iX, unsigned iY, unsigned iZ,  // Input dimensions
//      unsigned oX, unsigned oY, unsigned oZ,  // Output dimensions
//      cudaStream_t stream
//  ) {
//      if (iX > 1 && iY / 2 + 1 != iX)
//          REPORT_ERROR("windowFourierTransform ERROR: the Fourier transform should be of an image with equal sizes in all dimensions!");
//
//      if (oY == iX)
//          REPORT_ERROR("windowFourierTransform ERROR: there is a one-to-one map between input and output!");
//
//      cudaMemInit<XFLOAT>(d_out_real, 0, (size_t) oX * oY * oZ, stream);
//      cudaMemInit<XFLOAT>(d_out_imag, 0, (size_t) oX * oY * oZ, stream);
//
//      if (oY > iX) {
//          long int max_r2 = (iX - 1) * (iX - 1);
//
//          unsigned grid_dim = ceil((float)(iX*iY*iZ) / (float) WINDOW_FT_BLOCK_SIZE);
//          cuda_kernel_window_fourier_transform<true>
//          <<<grid_dim, WINDOW_FT_BLOCK_SIZE, 0, stream>>>
//          (
//              d_in_real,  d_in_imag,
//              d_out_real, d_out_imag,
//              iX, iY, iZ, iX * iY,  // Input dimensions
//              oX, oY, oZ, oX * oY,  // Output dimensions
//              iX * iY * iZ,
//              max_r2
//          );
//      } else {
//          unsigned grid_dim = ceil((float)(oX*oY*oZ) / (float) WINDOW_FT_BLOCK_SIZE);
//          cuda_kernel_window_fourier_transform<false>
//          <<<grid_dim, WINDOW_FT_BLOCK_SIZE, 0, stream>>>
//          (
//              d_in_real,  d_in_imag,
//              d_out_real, d_out_imag,
//              iX, iY, iZ, iX * iY,  // Input dimensions
//              oX, oY, oZ, oX * oY,  // Output dimensions
//              oX*oY*oZ
//          );
//      }
//  }

