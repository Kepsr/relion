#ifndef CUDA_HELPER_KERNELS_CUH_
#define CUDA_HELPER_KERNELS_CUH_

#include <cuda_runtime.h>
#include <curand.h>
#include <curand_kernel.h>
#include <vector>
#include <iostream>
#include <fstream>
#include "src/acc/cuda/cuda_settings.h"
#include "src/acc/cuda/cuda_kernels/cuda_device_utils.cuh"
#include "src/acc/acc_projector.h"
#include "src/acc/acc_projectorkernel_impl.h"

template<typename T>
__global__ void cuda_kernel_weights_exponent_coarse(
    T *g_pdf_orientation, bool *g_pdf_orientation_zeros,
    T *g_pdf_offset, bool *g_pdf_offset_zeros,
    T *g_weights, T g_min_diff2,
    int nr_coarse_orient, int nr_coarse_trans, int max_idx
) {
    const int i = blockIdx.x * SUMW_BLOCK_SIZE + threadIdx.x;
    if (i >= max_idx) return;
    const int itrans = i % nr_coarse_trans;
    const int iorient = (i - itrans) / nr_coarse_trans;
    const T diff2 = g_weights[i];
    if (diff2 < g_min_diff2 || g_pdf_orientation_zeros[iorient] || g_pdf_offset_zeros[itrans])
        g_weights[i] = -99e99;  // large negative number
    else
        g_weights[i] = g_pdf_orientation[iorient] + g_pdf_offset[itrans] + g_min_diff2 - diff2;
}

template<typename T>
__global__ void cuda_kernel_exponentiate(
    T *g_array, T add, size_t size
) {
    const int i = threadIdx.x + blockIdx.x * blockDim.x;
    if (i >= size) return;
    T a = g_array[i] + add;
    #ifdef ACC_DOUBLE_PRECISION
    g_array[i] = a < -700.0 ? 0.0 : exp(a);
    #else
    g_array[i] = a < -88.f ? 0.f : expf(a);
    #endif
}

template<bool DATA3D>
__global__ void cuda_kernel_collect2jobs(
    XFLOAT *g_oo_otrans_x,          // otrans-size -> make const
    XFLOAT *g_oo_otrans_y,          // otrans-size -> make const
    XFLOAT *g_oo_otrans_z,          // otrans-size -> make const
    XFLOAT *g_myp_oo_otrans_x2y2z2, // otrans-size -> make const
    XFLOAT *g_i_weights,
    XFLOAT op_significant_weight,    // TODO Put in const
    XFLOAT op_sum_weight,            // TODO Put in const
    int   coarse_trans,
    int   oversamples_trans,
    int   oversamples_orient,
    int   oversamples,
    bool  do_ignore_pdf_direction,
    XFLOAT *g_o_weights,
    XFLOAT *g_thr_wsum_prior_offsetx_class,
    XFLOAT *g_thr_wsum_prior_offsety_class,
    XFLOAT *g_thr_wsum_prior_offsetz_class,
    XFLOAT *g_thr_wsum_sigma2_offset,
    unsigned long * d_rot_idx,
    unsigned long * d_trans_idx,
    unsigned long * d_job_idx,
    unsigned long * d_job_num
) {
    extern __shared__ XFLOAT buffer[];
    XFLOAT* s_o_weights 				   = buffer;
    XFLOAT* s_thr_wsum_sigma2_offset       = buffer + SUMW_BLOCK_SIZE;
    XFLOAT* s_thr_wsum_prior_offsetx_class = buffer + SUMW_BLOCK_SIZE * 2;
    XFLOAT* s_thr_wsum_prior_offsety_class = buffer + SUMW_BLOCK_SIZE * 3;
    XFLOAT* s_thr_wsum_prior_offsetz_class = 0;
    if (DATA3D)
    s_thr_wsum_prior_offsetz_class = buffer + SUMW_BLOCK_SIZE * 4;

    s_o_weights             [threadIdx.x] = 0.0;
    s_thr_wsum_sigma2_offset[threadIdx.x]  = 0.0;

    s_thr_wsum_prior_offsetx_class[threadIdx.x] = 0.0;
    s_thr_wsum_prior_offsety_class[threadIdx.x] = 0.0;
    if (DATA3D)
    s_thr_wsum_prior_offsety_class[threadIdx.x] = 0.0;

    long int pos = d_job_idx[blockIdx.x] + threadIdx.x;  // pos is thread-resolved
    int job_size = d_job_num[blockIdx.x];

    const int pass_num = ceilfracf(job_size, SUMW_BLOCK_SIZE);
    __syncthreads();
    for (int pass = 0; pass < pass_num && pass * SUMW_BLOCK_SIZE + threadIdx.x < job_size; pass++, pos += SUMW_BLOCK_SIZE) {
        // loop the available warps enough to complete all translations for this orientation
        // if there is a translation that needs to be done still for this thread
            // index of comparison
            long int iy = d_trans_idx[pos];              // ...and its own trans...

            XFLOAT weight = g_i_weights[pos];
            if( weight >= op_significant_weight ) //TODO Might be slow (divergent threads)
                weight /= op_sum_weight;
            else
                weight = (XFLOAT)0.0;

            s_o_weights[threadIdx.x] 					+= weight;
            s_thr_wsum_sigma2_offset[threadIdx.x]       += weight * g_myp_oo_otrans_x2y2z2[iy];
            s_thr_wsum_prior_offsetx_class[threadIdx.x] += weight *          g_oo_otrans_x[iy];
            s_thr_wsum_prior_offsety_class[threadIdx.x] += weight *          g_oo_otrans_y[iy];
            if (DATA3D)
            s_thr_wsum_prior_offsetz_class[threadIdx.x] += weight *          g_oo_otrans_z[iy];
    }

    // Reduction of all translations this orientation
    __syncthreads();
    for (int j = SUMW_BLOCK_SIZE / 2; j > 0; j /= 2) {
        if (threadIdx.x < j) {
            s_thr_wsum_sigma2_offset      [threadIdx.x] += s_thr_wsum_sigma2_offset      [threadIdx.x + j];
            s_o_weights                   [threadIdx.x] += s_o_weights                   [threadIdx.x + j];
            s_thr_wsum_prior_offsetx_class[threadIdx.x] += s_thr_wsum_prior_offsetx_class[threadIdx.x + j];
            s_thr_wsum_prior_offsety_class[threadIdx.x] += s_thr_wsum_prior_offsety_class[threadIdx.x + j];
            if (DATA3D)
            s_thr_wsum_prior_offsetz_class[threadIdx.x] += s_thr_wsum_prior_offsetz_class[threadIdx.x + j];
        }
        __syncthreads();
    }

    g_o_weights                   [blockIdx.x] = s_o_weights[0];
    g_thr_wsum_sigma2_offset      [blockIdx.x] = s_thr_wsum_sigma2_offset[0];
    g_thr_wsum_prior_offsetx_class[blockIdx.x] = s_thr_wsum_prior_offsetx_class[0];
    g_thr_wsum_prior_offsety_class[blockIdx.x] = s_thr_wsum_prior_offsety_class[0];
    if (DATA3D)
    g_thr_wsum_prior_offsetz_class[blockIdx.x] = s_thr_wsum_prior_offsetz_class[0];
}

__global__ void cuda_kernel_exponentiate_weights_fine(
    XFLOAT *g_pdf_orientation, bool *g_pdf_orientation_zeros,
    XFLOAT *g_pdf_offset, bool *g_pdf_offset_zeros,
    XFLOAT *g_weights,
    XFLOAT min_diff2,
    int oversamples_orient, int oversamples_trans,
    unsigned long * d_rot_id, unsigned long * d_trans_idx,
    unsigned long * d_job_idx, unsigned long * d_job_num,
    long int job_num
);

__global__ void cuda_kernel_initRND(unsigned long seed, curandState *States);

__global__ void cuda_kernel_RNDnormalDitributionComplexWithPowerModulation2D(
    acc::Complex *Image, curandState *States, long int xdim, XFLOAT* spectra
);

__global__ void cuda_kernel_RNDnormalDitributionComplexWithPowerModulation3D(
    acc::Complex *Image, curandState *States, long int xdim, long int ydim, XFLOAT* spectra
);

// Basically just AccDataTypes::Image<XFLOAT>
template <typename T, int D>
struct image_view {

    T* data;
    long int size;
    long int dim[D];
    long int init[D];

};

__global__ void cuda_kernel_softMaskOutsideMap(
    image_view<XFLOAT, 3> vol, bool do_Mnoise,
    XFLOAT radius, XFLOAT radius_p, XFLOAT cosine_width
);

__global__ void cuda_kernel_softMaskBackgroundValue(
    image_view<XFLOAT, 3> vol,
    XFLOAT radius, XFLOAT radius_p, XFLOAT cosine_width,
    XFLOAT *g_sum, XFLOAT *g_sum_bg
);

__global__ void cuda_kernel_cosineFilter(
    image_view<XFLOAT, 3> vol,
    bool do_snoise, XFLOAT *noise,
    XFLOAT radius, XFLOAT radius_p, XFLOAT cosine_width,
    XFLOAT sum_bg_total
);

__global__ void cuda_kernel_initOrientations(RFLOAT *pdfs, XFLOAT *pdf_orientation, bool *pdf_orientation_zeros, size_t sz);

//----------------------------------------------------------------------------
namespace CudaKernels {

template <typename T>
__global__ void cuda_kernel_translate2D(
    T const * g_image_in, T * g_image_out,
    int image_size, int xdim, int ydim, int dx, int dy
) {
    const int i = threadIdx.x + blockIdx.x * blockDim.x;
    if (i >= image_size) return;
    const int x = i % xdim;
    const int xp = x + dx;
    const int yp = (i - x) / xdim + dy;
    if (yp < 0 || xp < 0 || yp >= ydim || xp >= xdim) return;
    const int j = yp * xdim + xp;
    // If displacement is negative, j could be negative
    if (j >= 0 && j < image_size)
        g_image_out[j] = g_image_in[i];
}

template <typename T>
__global__ void cuda_kernel_translate3D(
    T const * g_image_in, T * g_image_out,
    int image_size, int xdim, int ydim, int zdim,
    int dx, int dy, int dz
) {
    const int i = threadIdx.x + blockIdx.x * blockDim.x;
    if (i >= image_size) return;
    const int xydim = xdim * ydim;
    const int xy = i % xydim;
    const int zp = i / xydim + dz;
    const int yp = xy / xdim + dy;
    const int xp = xy % xdim + dx;
    if (zp < 0 || yp < 0 || xp < 0 || zp >= zdim || yp >= ydim || xp >= xdim) return;
    const int j = zp * xydim +  yp * xdim + xp;
    if (j >= 0 && j < image_size) // if displacement is negative, new_pixel could be less than 0
        g_image_out[j] = g_image_in[i];
}

}

//----------------------------------------------------------------------------
//__global__ void cuda_kernel_selfTranslate2D(	XFLOAT * g_image_in,
//												XFLOAT * g_image_out,
//												int image_size,
//												int xdim,
//												int ydim, //not used
//												int dx,
//												int dy);
//
//__global__ void cuda_kernel_selfTranslate3D(	XFLOAT * g_image_in,
//												XFLOAT * g_image_out,
//												int image_size,
//												int xdim,
//												int ydim,
//												int zdim, //not used
//												int dx,
//												int dy,
//												int dz);
//----------------------------------------------------------------------------
//__global__ void cuda_kernel_powerClass2D(	acc::Complex * g_image,
//											XFLOAT * g_spectrum,
//											int image_size,
//											int spectrum_size,
//											int xdim,
//											int ydim,
//											int res_limit,
//											XFLOAT * g_highres_Xi2);
//
//__global__ void cuda_kernel_powerClass3D(	acc::Complex * g_image,
//											XFLOAT * g_spectrum,
//											int image_size,
//											int spectrum_size,
//											int xdim,
//											int ydim,
//											int zdim,
//											int res_limit,
//											XFLOAT * g_highres_Xi2);

//----------------------------------------------------------------------------
__global__ void cuda_kernel_probRatio(  XFLOAT * d_Mccf,
                                        XFLOAT * d_Mpsi,
                                        XFLOAT * d_Maux,
                                        XFLOAT * d_Mmean,
                                        XFLOAT * d_Mstddev,
                                        int image_size,
                                        XFLOAT normfft,
                                        XFLOAT sum_ref_under_circ_mask,
                                        XFLOAT sum_ref2_under_circ_mask,
                                        XFLOAT expected_Pratio,
                                        int NpsiThisBatch,
                                        int startPsi,
                                        int totalPsis);

__global__ void cuda_kernel_rotateOnly(
    acc::Complex * d_Faux,
    XFLOAT psi, AccProjectorKernel projector, int startPsi
);

__global__ void cuda_kernel_rotateAndCtf(
    acc::Complex * d_Faux, XFLOAT const * d_ctf,
    XFLOAT psi, AccProjectorKernel projector, int startPsi = 0
);

/// In-place pointwise multiply complex array A by B, after conjugating A
__global__ void cuda_kernel_convol_A(acc::Complex * d_A, acc::Complex const * d_B, int image_size);

/// In-place pointwise multiply complex array A by B, after conjugating A. Write to C.
__global__ void cuda_kernel_convol_A(acc::Complex * d_A, acc::Complex const * d_B, acc::Complex * d_C, int image_size);

/// In-place pointwise multiply many complex arrays A by a single B after conjugating A
__global__ void cuda_kernel_batch_convol_A(acc::Complex * d_A, acc::Complex const * d_B, int image_size);

/// Pointwise multiply many complex arrays A by a single B, after conjugating A
__global__ void cuda_kernel_batch_convol_A(acc::Complex * d_A, acc::Complex const * d_B, acc::Complex * d_C, int image_size);

/// In-place pointwise multiply complex array A by B, after conjugating B
__global__ void cuda_kernel_convol_B(acc::Complex * d_A, acc::Complex const * d_B, int image_size);

/// In-place pointwise multiply complex array A by B, after conjugating B. Write to C.
__global__ void cuda_kernel_convol_B(acc::Complex * d_A, acc::Complex const * d_B, acc::Complex * d_C, int image_size);

/// In-place pointwise multiply many complex arrays A by a single one B, after conjugating B
__global__ void cuda_kernel_batch_convol_B(acc::Complex * d_A, acc::Complex const * d_B, int image_size);

/*
 * Multiply scalar array A by a scalar S
 *
 *  OUT[i] = A[i]*S
 */
template <typename T>
__global__ void cuda_kernel_multi(T const* A, T *OUT, T S, int image_size) {
    const int i = threadIdx.x + blockIdx.x * blockDim.x;
    if (i < image_size) OUT[i] = A[i] * S;
}

namespace CudaKernels {
/*
 * In place multiplies scalar array A by a scalar S
 *
 *  A[i] = A[i]*S
 */
template <typename T>
__global__ void cuda_kernel_multi(T *A, T S, int image_size) {
    const int i = threadIdx.x + blockIdx.x * blockDim.x;
    if (i < image_size) A[i] *= S;
}

}

/*
 * Multiply scalar array A by scalar array B and a scalar S, pixel-by-pixel
 *
 *  OUT[i] = A[i]*B[i]*S
 */
template <typename T>
__global__ void cuda_kernel_multi(T const* A, T const* B, T* OUT, T S, int image_size) {
    const int i = threadIdx.x + blockIdx.x * blockDim.x;
    if (i < image_size) OUT[i] = A[i] * B[i] * S;
}

__global__ void cuda_kernel_finalizeMstddev(XFLOAT *Mstddev, XFLOAT *aux, XFLOAT S, int image_size);

/// Square values in-place
__global__ void cuda_kernel_square(XFLOAT *A, int image_size);

/*
 * Casts on device so we can copy_to_host directly into a multidimarray.
 */
template <typename T1, typename T2 >
__global__ void cuda_kernel_cast(T1 const* IN, T2 *OUT, int size) {
    const int i = threadIdx.x + blockIdx.x * blockDim.x;
    if (i < size) OUT[i] = IN[i];
}

template<bool do_highpass>
__global__ void cuda_kernel_frequencyPass(
    acc::Complex *A, long int ori_size,
    size_t Xdim, size_t Ydim, size_t Zdim,
    XFLOAT edge_low, XFLOAT edge_width, XFLOAT edge_high,
    XFLOAT angpix, int image_size
) {
    const int i = threadIdx.x + blockIdx.x * blockDim.x;
    if (i >= image_size) return;

    const int z = i / (Xdim * Ydim);
    const int xy = (i - z * Xdim * Ydim);
    const int y = xy / Xdim;

    const int xp = xy - y * Xdim;
    const int yp = y < Xdim ? y : y - Ydim;
    const int zp = z < Xdim ? z : z - Zdim;

    const RFLOAT r2 = xp * xp + yp * yp + zp * zp;
    const RFLOAT res = sqrt(r2) / ori_size;
    if (do_highpass) {
        if (res < edge_low) {
            // highpass => lows are dead
            A[i].x = 0.0;
            A[i].y = 0.0;
        } else if (res < edge_high) {
            // highpass => medium lows are almost dead
            const XFLOAT factor = 0.5 * (1.0 - cos((res - edge_low) * PI / edge_width));
            A[i].x *= factor;
            A[i].y *= factor;
        }
    } else {  // lowpass
        if (res > edge_high) {
            // lowpass => highs are dead
            A[i].x = 0.0;
            A[i].y = 0.0;
        } else if (res > edge_low) {
            //lowpass => medium highs are almost dead
            const XFLOAT factor = 0.5 + (1.0 + cos((res - edge_low) * PI / edge_width));
            A[i].x *= factor;
            A[i].y *= factor;
        }
    }
}

template<bool DATA3D>
__global__ void cuda_kernel_powerClass(
    acc::Complex const* g_image, XFLOAT* g_spectrum,
    int image_size, int spectrum_size,
    int xdim, int ydim, int zdim,
    int res_limit,
    XFLOAT * g_highres_Xi2
) {
    __shared__ XFLOAT s_highres_Xi2[POWERCLASS_BLOCK_SIZE];

    s_highres_Xi2[threadIdx.x] = 0.0;

    const int xydim = xdim * ydim;
    const int i = threadIdx.x + blockIdx.x * POWERCLASS_BLOCK_SIZE;
    if (i < image_size) {
        XFLOAT d;
        bool coords_in_range;
        if (DATA3D) {
                  int z =  i / xydim;
            const int xy = i % xydim;
                  int y =  xy / xdim;
            const int x =  xy % xdim;
            y = y < xdim ? y : y - ydim;
            z = z < xdim ? z : z - zdim;
            d = x * x + y * y + z * z;
            coords_in_range = x != 0 || y >= 0.0 || z >= 0.f;
        } else {
            const int x = i % xdim;
                  int y = (i - x) / xdim;
            y = y < xdim ? y : y - ydim;
            d = x * x + y * y;
            coords_in_range = x != 0 || y >= 0.0;
        }

        #if defined(ACC_DOUBLE_PRECISION)
        const int ires = __double2int_rn(sqrt(d));
        #else
        const int ires = __float2int_rn(sqrtf(d));
        #endif
        if (ires > 0.f && ires < spectrum_size && coords_in_range) {
            const XFLOAT normFaux = g_image[i].x * g_image[i].x + g_image[i].y * g_image[i].y;
            cuda_atomic_add(g_spectrum + ires, normFaux);
            if (ires >= res_limit)
                s_highres_Xi2[threadIdx.x] = normFaux;
        }
    }

    // Reduce the higres_Xi2-values for all threads.
    // (I tried a straight atomic-write: for 128 threads it was ~3x slower)
    __syncthreads();
    for (int j = POWERCLASS_BLOCK_SIZE / 2; j > 0; j /= 2) {
        if (threadIdx.x < j)
            s_highres_Xi2[threadIdx.x] += s_highres_Xi2[threadIdx.x + j];
        __syncthreads();
    }
    if (threadIdx.x == 0)
        cuda_atomic_add(g_highres_Xi2, s_highres_Xi2[0]);
}

template<bool invert>
__global__ void cuda_kernel_make_eulers_2D(
    XFLOAT *alphas, XFLOAT *eulers, unsigned orientation_num
);

template<bool invert,bool doL, bool doR>
__global__ void cuda_kernel_make_eulers_3D(
    XFLOAT *alphas, XFLOAT *betas, XFLOAT *gammas, XFLOAT *eulers,
    unsigned orientation_num,
    XFLOAT *L, XFLOAT *R
);

#define INIT_VALUE_BLOCK_SIZE 512
template <typename T>
__global__ void cuda_kernel_init_complex_value(
    T *data, XFLOAT value, size_t size
) {
    const size_t i = blockIdx.x * INIT_VALUE_BLOCK_SIZE + threadIdx.x;
    if (i >= size) return;
    data[i].x = value;
    data[i].y = value;
}

template <typename T>
__global__ void cuda_kernel_init_value(
    T *data, T value, size_t size
) {
    const size_t i = blockIdx.x * INIT_VALUE_BLOCK_SIZE + threadIdx.x;
    if (i < size)
        data[i] = value;
}

#define WEIGHT_MAP_BLOCK_SIZE 512
__global__ void cuda_kernel_allweights_to_mweights(
    unsigned long* d_iorient,
    XFLOAT* d_allweights, XFLOAT* d_mweights,
    unsigned long orientation_num, unsigned long translation_num,
    int block_size
);

#define OVER_THRESHOLD_BLOCK_SIZE 512
template <typename T>
__global__ void cuda_kernel_array_over_threshold(
    T const* data, bool *passed, T threshold, size_t size
) {
    const size_t i = blockIdx.x * OVER_THRESHOLD_BLOCK_SIZE + threadIdx.x;
    if (i >= size) return;
    passed[i] = data[i] >= threshold;
}

#define FIND_IN_CUMULATIVE_BLOCK_SIZE 512
template <typename T>
__global__ void cuda_kernel_find_threshold_idx_in_cumulative(
    T const* data, T threshold, size_t size_m1 /*data size minus 1*/, size_t *idx
) {
    const size_t i = blockIdx.x * FIND_IN_CUMULATIVE_BLOCK_SIZE + threadIdx.x;
    if (i < size_m1 && data[i] <= threshold && threshold < data[i + 1])
        idx[0] = i + 1;
}

#define WINDOW_FT_BLOCK_SIZE 128
template<bool check_max_r2>
__global__ void cuda_kernel_window_fourier_transform(
    XFLOAT const* g_in_real,  XFLOAT const* g_in_imag,
    XFLOAT *g_out_real, XFLOAT *g_out_imag,
    unsigned iX, unsigned iY, unsigned iZ, unsigned iYX,  // Input dimensions
    unsigned oX, unsigned oY, unsigned oZ, unsigned oYX,  // Output dimensions
    unsigned max_idx,
    unsigned max_r2 = 0
) {
    const unsigned n = threadIdx.x + WINDOW_FT_BLOCK_SIZE * blockIdx.x;
    if (n >= max_idx) return;
    const long int image_offset = oX * oY * oZ * blockIdx.y;
    const int iXiY = iX * iY;
    const int oXoY = oX * oY;
    int i, k, ip, jp, kp;
    if (check_max_r2) {
        k = n / iXiY;
        i = n % iXiY / iX;

        kp = k < iX ? k : k - iZ;
        ip = i < iX ? i : i - iY;
        jp = n % iX;

        if (kp * kp + ip * ip + jp * jp > max_r2) return;
    } else {
        k = n / oXoY;
        i = n % oXoY / oX;

        kp = k < oX ? k : k - oZ;
        ip = i < oX ? i : i - oY;
        jp = n % oX;
    }
    const int ipi = ip < 0 ? ip + iY : ip;
    const int kpi = kp < 0 ? kp + iZ : kp;
    const int ipo = ip < 0 ? ip + oY : ip;
    const int kpo = kp < 0 ? kp + oZ : kp;
    const int ikpi = kpi * iYX + ipi * iX + jp + image_offset;
    const int ikpo = kpo * oYX + ipo * oX + jp + image_offset;
    g_out_real[ikpo] = g_in_real[ikpi];
    g_out_imag[ikpo] = g_in_imag[ikpi];
}

#define WINDOW_FT_BLOCK_SIZE 128
template<bool check_max_r2>
__global__ void cuda_kernel_window_fourier_transform(
    acc::Complex const* g_in, acc::Complex *g_out,
    size_t iX, size_t iY, size_t iZ, size_t iYX, //Input dimensions
    size_t oX, size_t oY, size_t oZ, size_t oYX, //Output dimensions
    size_t max_idx,
    size_t max_r2 = 0
) {
    const size_t n = threadIdx.x + WINDOW_FT_BLOCK_SIZE * blockIdx.x;
    if (n >= max_idx) return;
    const size_t oOFF = oX * oY * oZ * blockIdx.y;
    const size_t iOFF = iX * iY * iZ * blockIdx.y;
    long int i, k, ip, jp, kp;
    if (check_max_r2) {
        k = n / (iX * iY);
        i = n % (iX * iY) / iX;
        kp = k < iX ? k : k - iZ;
        ip = i < iX ? i : i - iY;
        jp = n % iX;
        if (kp * kp + ip * ip + jp * jp > max_r2) return;
    } else {
        k = n / (oX * oY);
        i = n % (oX * oY) / oX;
        kp = k < oX ? k : k - oZ;
        ip = i < oX ? i : i - oY;
        jp = n % oX;
    }
    const long int idxi = (kp < 0 ? kp + iZ : kp) * iYX + (ip < 0 ? ip + iY : ip) * iX + jp;
    const long int idxo = (kp < 0 ? kp + oZ : kp) * oYX + (ip < 0 ? ip + oY : ip) * oX + jp;
    g_out[idxo + oOFF] = g_in[idxi + iOFF];
}

#define NEAREST_NEIGHBOUR 0
#define TRILINEAR 1
__global__ void cuda_kernel_griddingCorrect(RFLOAT *vol, int interpolator, RFLOAT rrval, RFLOAT r_min_nn,
                                            size_t iX, size_t iY, size_t iZ);

template <typename T>
__global__ void cuda_kernel_window_transform(
    T * d_in, T * d_out,
    int iszX, int iszY, int iszZ,  // Input dimensions
    int oftX, int oftY, int oftZ, int oszX, int oszY, int oszZ  // Output dimensions
) {
    const int idx = blockIdx.x * blockDim.x + threadIdx.x;
    const int idy = blockIdx.y * blockDim.y + threadIdx.y;
    const int idz = blockIdx.z * blockDim.z + threadIdx.z;
    if (idx >= oszX || idy >= oszY || idz >= oszZ) return;
    d_out[idz * oszX * oszY + idy * oszX + idx] = 
        idx >= oftX && idx < oftX + iszX &&
        idy >= oftY && idy < oftY + iszY &&
        idz >= oftZ && idz < oftZ + iszZ ?
        d_in[(idz - oftZ) * iszX * iszY + (idy - oftY) * iszX + (idx - oftX)] : 0.0;
}

template <typename T>
__device__ inline void swap(T& a, T& b) {
    const T c = a;
    a = b;
    b = c;
}

__device__ inline float floordivf(float a, float b) {
    return floorf(a / b);
}

template <typename T>
__global__ void cuda_kernel_centerFFT_2D(
    T *img_in, int image_size, int xdim, int ydim, int xshift, int yshift
) {
    const int i = threadIdx.x + blockIdx.x * blockDim.x;
    if (i >= image_size / 2) return;
    const long int image_offset = image_size * blockIdx.y;
    const int y = floordivf(i, xdim);
    const int x = i % xdim;	// also = i - y*xdim, but this depends on y having been calculated, i.e. serial evaluation
    const int xp = (x + xshift + xdim) % xdim;
    const int yp = (y + yshift + ydim) % ydim;
    const int j = yp * xdim + xp;
    swap(img_in[image_offset + i], img_in[image_offset + j]);
}

template __global__ void cuda_kernel_centerFFT_2D<double>(double*, int, int, int, int, int);
template __global__ void cuda_kernel_centerFFT_2D<float>(float*, int, int, int, int, int);

template <typename T>
__global__ void cuda_kernel_centerFFT_3D(
    T *img_in,
    int image_size,
    int xdim, int ydim, int zdim,
    int xshift, int yshift, int zshift
) {
    const int i = threadIdx.x + blockIdx.x * blockDim.x;
    if (i >= image_size / 2) return;
    const long int image_offset = image_size * blockIdx.y;
    const int xydim = xdim * ydim;
    const int z = floordivf(i, xydim);
    const int xy = i % xydim;
    const int y = floordivf(xy, xdim);
    const int x = xy % xdim;
    const int xp = (x + xshift + xdim) % xdim;
    const int yp = (y + yshift + ydim) % ydim;
    const int zp = (z + zshift + zdim) % zdim;
    const int j = zp * xydim + yp * xdim + xp;
    swap(img_in[image_offset + i], img_in[image_offset + j]);
}

template __global__ void cuda_kernel_centerFFT_3D<double>(double*, int, int, int, int, int, int, int);
template __global__ void cuda_kernel_centerFFT_3D<float>(float*, int, int, int, int, int, int, int);

template <typename T>
__global__ void cuda_kernel_centerFFTbySign(
    T *img_in, int xdim, int ydim, int zdim
) {
    const int x = threadIdx.x + blockIdx.x * blockDim.x;
    const int y = threadIdx.y + blockIdx.y * blockDim.y;
    const int z = threadIdx.z + blockIdx.z * blockDim.z;
    const int i = z * xdim * ydim + y * xdim + x;
    if (x >= xdim || y >= ydim || z >= zdim) return;
    if ((x ^ y ^ z) & 1 == 0) return;
    img_in[i].x = -img_in[i].x;
    img_in[i].y = -img_in[i].y;
}

template __global__ void cuda_kernel_centerFFTbySign<double2>(double2*, int, int, int);
template __global__ void cuda_kernel_centerFFTbySign<float2>(float2*, int, int, int);

#if !defined(__CUDA_ARCH__) || __CUDA_ARCH__ >= 600
#else
__device__ double atomicAdd(double* address, double val) {
    using ull = unsigned long long int;
    ull* address_as_ull = (ull*) address;
    ull old = *address_as_ull, assumed;
    do {
        assumed = old;
        old = atomicCAS(address_as_ull, assumed,
                        __double_as_longlong(val + __longlong_as_double(assumed)));
    } while (assumed != old);
    return __longlong_as_double(old);
}
#endif

template <typename T>
__global__ void cuda_kernel_calcPowerSpectrum(
    T const* dFaux, int padoridim, T *ddata, int data_sz, RFLOAT *dpower_spectrum, RFLOAT *dcounter,
    int max_r2, int min_r2, RFLOAT normfft, RFLOAT padding_factor, RFLOAT weight,
    RFLOAT const* dfourier_mask, int fx, int fy, int fz, bool do_fourier_mask
) {
    const int idx = blockIdx.x * blockDim.x + threadIdx.x;
    const int idy = blockIdx.y * blockDim.y + threadIdx.y;
    const int idz = blockIdx.z * blockDim.z + threadIdx.z;
    const int XSIZE = padoridim / 2 + 1;
    if (idx >= XSIZE || idy >= padoridim || idz >= padoridim) return;
    const int dx = data_sz / 2 + 1;
    const int dxy = blockDim.z == 1 ? 0 : data_sz * dx;
    const int jp = idx;
    const int ip = idy < XSIZE ? idy : idy - padoridim;
    const int kp = idz < XSIZE ? idz : idz - padoridim;
    const int r2 = kp * kp + ip * ip + jp * jp;
    // The Fourier Transforms are all "normalised" for 2D transforms of size = ori_size x ori_size
    if (r2 > max_r2) return;
    // if (do_fourier_mask) weight = FFTW::elem(*fourier_mask, round(ip/padding_factor), round(jp/padding_factor), round(kp/padding_factor));
    if (do_fourier_mask) {
        int lkp = round(kp / padding_factor);
        int lip = round(ip / padding_factor);
        int ljp = round(jp / padding_factor);
        lkp = lkp < 0 ? lkp + fz : lkp;
        lip = lip < 0 ? lip + fy : lip;
        weight = dfourier_mask[lkp * fy * fx + lip * fx + ljp];
    }
    // Set data array
    T val = dFaux[idz * XSIZE * padoridim + idy * XSIZE + idx];
    val.x *= normfft * weight;
    val.y *= normfft * weight;
    // data.elem(ip, jp, kp) = weight * direct::elem(Faux, i, j, k) * normfft;
    const long int n = (kp + data_sz / 2) * dxy + (ip + data_sz / 2) * dx + jp;
    ddata[n] = val;

    // Calculate power spectrum
    const int ires = round(sqrt((RFLOAT) r2) / padding_factor);
    // Factor two because of two-dimensionality of the complex plane
    const RFLOAT norm = (val.x * val.x + val.y * val.y);
    atomicAdd(dpower_spectrum + ires, norm / 2.0);    // direct::elem(power_spectrum, ires) += norm(data.elem(ip, jp, kp)) / 2.0;
    atomicAdd(dcounter        + ires, weight);  // direct::elem(counter, ires) += weight;
    // Apply high pass filter of the reference only after calculating the power spectrum
    val.x = val.y = 0.0;
    if (r2 <= min_r2)
        ddata[n] = val;  // data.elem(ip, jp, kp) = 0;
}

template __global__ void cuda_kernel_calcPowerSpectrum(
    double2 const*, int, double2*, int, RFLOAT*, RFLOAT*, int, int, RFLOAT, RFLOAT, RFLOAT,
    RFLOAT const*, int, int, int, bool
);

template __global__ void cuda_kernel_calcPowerSpectrum(
    float2 const*, int, float2*, int, RFLOAT*, RFLOAT*, int, int, RFLOAT, RFLOAT, RFLOAT,
    RFLOAT const*, int, int, int, bool
);

__global__ void cuda_kernel_updatePowerSpectrum(RFLOAT *dcounter, RFLOAT *dpower_spectrum, int sz);
#endif /* CUDA_HELPER_KERNELS_CUH_ */
