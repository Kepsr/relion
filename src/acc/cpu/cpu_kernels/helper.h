#ifndef HELPER_KERNELS_H_
#define HELPER_KERNELS_H_

#include <vector>
#include <iostream>
#include <fstream>
#include "src/macros.h"
#include "src/acc/cpu/cpu_settings.h"
#include "src/acc/cpu/cpu_kernels/cpu_utils.h"
#include "src/acc/acc_projector.h"
#include "src/acc/acc_projectorkernel_impl.h"

#ifdef __INTEL_COMPILER
#define ALWAYS_INLINE_GCC inline
#else
#define ALWAYS_INLINE_GCC __attribute__((always_inline)) inline
#endif

#ifdef ACC_DOUBLE_PRECISION
#define ACC_SINCOS( theta, sin_ptr, cos_ptr ) sincos(theta, sin_ptr, cos_ptr)
inline double raised_cos(double theta) { return 0.5 * (1.0 + cos(theta)); }
#else
#define ACC_SINCOS( theta, sin_ptr, cos_ptr ) sincosf(theta, sin_ptr, cos_ptr)
inline float raised_cos(float theta) { return 0.5 * (1.0 + cos(theta)); }
#endif

namespace CpuKernels {

template<typename T>
ALWAYS_INLINE_GCC void weights_exponent_coarse(
    T *g_pdf_orientation, bool *g_pdf_orientation_zeros,
    T *g_pdf_offset, bool *g_pdf_offset_zeros,
    T *g_weights, T g_min_diff2,
    unsigned long nr_coarse_orient,
    unsigned long nr_coarse_trans,
    size_t max_idx
) {
    for (size_t idx = 0; idx < max_idx; idx++) {
        const unsigned long  itrans = idx % nr_coarse_trans;
        const unsigned long  iorient = (idx - itrans) / nr_coarse_trans;
        T& w = g_weights[idx];
        if (w < g_min_diff2 || g_pdf_orientation_zeros[iorient] || g_pdf_offset_zeros[itrans])
            w = std::numeric_limits<T>::lowest();
        else
            w = g_pdf_orientation[iorient] + g_pdf_offset[itrans] + g_min_diff2 - w;
    }
}


template<typename T>
ALWAYS_INLINE_GCC void exponentiate(T *g_array, T add, size_t size) {
    for (size_t idx = 0; idx < size; idx++) {
        const T a = g_array[idx] + add;
        #ifdef ACC_DOUBLE_PRECISION
        g_array[idx] = a < -700.0 ? 0.0 : exp(a);
        #else
        g_array[idx] = a < -88.f ? 0.f : expf(a);
        #endif
    }
}

template<bool DATA3D>
ALWAYS_INLINE_GCC void collect2jobs(
    int grid_size, int block_size,
    XFLOAT *g_oo_otrans_x,          // otrans-size -> make const
    XFLOAT *g_oo_otrans_y,          // otrans-size -> make const
    XFLOAT *g_oo_otrans_z,          // otrans-size -> make const
    XFLOAT *g_myp_oo_otrans_x2y2z2, // otrans-size -> make const
    XFLOAT *g_i_weights,
    XFLOAT op_significant_weight,    // TODO Put in const
    XFLOAT op_sum_weight,            // TODO Put in const
    unsigned long coarse_trans,
    unsigned long oversamples_trans, unsigned long oversamples_orient,
    unsigned long oversamples,
    bool do_ignore_pdf_direction,
    XFLOAT *g_o_weights,
    XFLOAT *g_thr_wsum_prior_offsetx_class,
    XFLOAT *g_thr_wsum_prior_offsety_class,
    XFLOAT *g_thr_wsum_prior_offsetz_class,
    XFLOAT *g_thr_wsum_sigma2_offset,
    unsigned long *d_rot_idx, unsigned long *d_trans_idx,
    unsigned long *d_job_idx, unsigned long *d_job_num
) {
    // block id
    for (int bid = 0; bid < grid_size; bid++) {

        XFLOAT s_o_weights[block_size];
        XFLOAT s_thr_wsum_sigma2_offset[block_size];;
        XFLOAT s_thr_wsum_prior_offsetx_class[block_size];
        XFLOAT s_thr_wsum_prior_offsety_class[block_size];
        XFLOAT s_thr_wsum_prior_offsetz_class[block_size];

        unsigned long pos = d_job_idx[bid];
        unsigned long job_size = d_job_num[bid];

        int pass_num = ceilfracf(job_size, block_size);

        for(int tid = 0; tid < block_size; tid++) {
            s_o_weights[tid]                    = 0.0;
            s_thr_wsum_sigma2_offset[tid]       = 0.0;
            s_thr_wsum_prior_offsetx_class[tid] = 0.0;
            s_thr_wsum_prior_offsety_class[tid] = 0.0;
            if (DATA3D)
            s_thr_wsum_prior_offsety_class[tid] = 0.0;
        }

        for (int pass = 0; pass < pass_num; pass++, pos+=block_size) {
            // loop the available warps enough to complete all translations for this orientation
            for (int tid = 0; tid < block_size; tid++) {
                if (pass * block_size + tid >= job_size) continue;
                // if there is a translation that needs to be done still for this thread
                // index of comparison
                long int iy = d_trans_idx[pos + tid];              // ...and its own trans...

                XFLOAT weight = g_i_weights[pos + tid];
                if (weight >= op_significant_weight)
                    // TODO Might be slow (divergent threads)
                    weight /= op_sum_weight;
                else
                    weight = 0.0;

                s_o_weights[tid]                    += weight;
                s_thr_wsum_prior_offsetx_class[tid] += weight *          g_oo_otrans_x[iy];
                s_thr_wsum_prior_offsety_class[tid] += weight *          g_oo_otrans_y[iy];
                s_thr_wsum_sigma2_offset[tid]       += weight * g_myp_oo_otrans_x2y2z2[iy];
            }
        }

        for (int tid = 1; tid < block_size; tid++) {
            s_o_weights[0]                    += s_o_weights[tid];
            s_thr_wsum_sigma2_offset[0]       += s_thr_wsum_sigma2_offset[tid];
            s_thr_wsum_prior_offsetx_class[0] += s_thr_wsum_prior_offsetx_class[tid];
            s_thr_wsum_prior_offsety_class[0] += s_thr_wsum_prior_offsety_class[tid];
            if (DATA3D)
            s_thr_wsum_prior_offsetz_class[0] += s_thr_wsum_prior_offsetz_class[tid];
        }
        g_o_weights[bid]			        = s_o_weights[0];
        g_thr_wsum_sigma2_offset[bid]       = s_thr_wsum_sigma2_offset[0];
        g_thr_wsum_prior_offsetx_class[bid] = s_thr_wsum_prior_offsetx_class[0];
        g_thr_wsum_prior_offsety_class[bid] = s_thr_wsum_prior_offsety_class[0];
        if (DATA3D)
        g_thr_wsum_prior_offsetz_class[bid] = s_thr_wsum_prior_offsetz_class[0];
    }
}

void exponentiate_weights_fine(
    XFLOAT *g_pdf_orientation, bool *g_pdf_orientation_zeros,
    XFLOAT *g_pdf_offset, bool *g_pdf_offset_zeros,
    XFLOAT *g_weights, XFLOAT min_diff2,
    unsigned long oversamples_orient, unsigned long oversamples_trans,
    unsigned long *d_rot_id, unsigned long *d_trans_idx,
    unsigned long *d_job_idx, unsigned long *d_job_num,
    long int job_num
);

void RNDnormalDitributionComplexWithPowerModulation2D(acc::Complex* Image, size_t xdim, XFLOAT *spectra);
void RNDnormalDitributionComplexWithPowerModulation3D(acc::Complex* Image, size_t xdim, size_t ydim, XFLOAT *spectra);

void softMaskBackgroundValue(
    size_t block_dim, size_t block_size,
    XFLOAT *vol, size_t vol_size,
    long int xdim, long int ydim, long int zdim,
    long int xinit, long int yinit, long int zinit,
    XFLOAT radius, XFLOAT radius_p, XFLOAT cosine_width,
    XFLOAT *g_sum, XFLOAT *g_sum_bg
);

void cosineFilter(
    size_t block_dim, size_t block_size,
    XFLOAT *vol, size_t vol_size,
    long int xdim, long int ydim, long int zdim,
    long int xinit, long int yinit, long int zinit,
    bool do_Mnoise, XFLOAT *noise,
    XFLOAT radius, XFLOAT radius_p, XFLOAT cosine_width, XFLOAT sum_bg_total
);

//----------------------------------------------------------------------------

template <typename T>
void cpu_translate2D(
    T * g_image_in, T * g_image_out, size_t image_size,
    int xdim, int ydim /*not used*/, int dx, int dy
);

template <typename T>
void cpu_translate3D(
    T * g_image_in, T * g_image_out, size_t image_size,
    int xdim, int ydim, int zdim /*not used*/,
    int dx, int dy, int dz
);

//----------------------------------------------------------------------------
template <typename T>
void centerFFT_2D(
    size_t batch_size, size_t pixel_start, size_t pixel_end,
    T *img_in, size_t image_size, int xdim, int ydim, int xshift, int yshift
);

template <typename T>
void centerFFT_3D(
    size_t batch_size, size_t pixel_start, size_t pixel_end,
    T *img_in, size_t image_size,
    int xdim, int ydim, int zdim,
    int xshift, int yshift, int zshift
);
//----------------------------------------------------------------------------
/*void probRatio( int       blockIdx_x,
                int       threadIdx_x,
                XFLOAT *d_Mccf,
                XFLOAT *d_Mpsi,
                XFLOAT *d_Maux,
                XFLOAT *d_Mmean,
                XFLOAT *d_Mstddev,
                size_t image_size,
                XFLOAT  normfft,
                XFLOAT  sum_ref_under_circ_mask,
                XFLOAT sum_ref2_under_circ_mask,
                XFLOAT expected_Pratio,
                int NpsiThisBatch,
                int startPsi,
                int totalPsis);

void rotateOnly(int              blockIdx_x,
                int              blockIdx_y,
                int              threadIdx_x,
                acc::Complex     *d_Faux,
                XFLOAT           psi,
                AccProjectorKernel &projector,

void rotateAndCtf(  int              blockIdx_x,
                    int              blockIdx_y,
                    int              threadIdx_x,
                    acc::Complex     *d_Faux,
                    XFLOAT          *d_ctf,
                    XFLOAT           psi,
                    AccProjectorKernel &projector,
                    int              startPsi = 0);

|*
 * Multiplies complex array A (in-place) by B, pixel-by-pixel, after conjugating A
 *|
void convol_A(  int          blockIdx_x,
                int          threadIdx_x,
                acc::Complex *d_A,
                acc::Complex *d_B,
                size_t       image_size);

|*
 * Multiplies complex array A (in-place) by B, pixel-by-pixel, after conjugating A, writes to C
 *|
void convol_A(  int          blockIdx_x,
                int          threadIdx_x,
                acc::Complex *d_A,
                acc::Complex *d_B,
                acc::Complex *d_C,
                size_t       image_size);

|*
 * Multiplies many complex arrays A (in-place) by a single B, pixel-by-pixel, after conjugating A
 *|
void batch_convol_A(int           blockIdx_x,
                    int           threadIdx_x,
                    acc::Complex  *d_A,
                    acc::Complex  *d_B,
                    size_t        image_size);

|*
* Multiplies many complex arrays A (not in-place) by a single B, pixel-by-pixel, after conjugating A
*|
void batch_convol_A(int          blockIdx_x,
                    int          threadIdx_x,
                    acc::Complex *d_A,
                    acc::Complex *d_B,
                    acc::Complex *d_C,
                    size_t       image_size);

|*
 * Multiplies complex array A (in-place) by B, pixel-by-pixel, after conjugating B
 *|
void convol_B(  int          blockIdx_x,
                int          threadIdx_x,
                acc::Complex *d_A,
                acc::Complex *d_B,
                size_t       image_size);

|*
 * Multiplies complex array A (in-place) by B, pixel-by-pixel, after conjugating B, writes to C
 *|
void convol_B(  int       blockIdx_x,
                int       threadIdx_x,
                acc::Complex *d_A,
                acc::Complex *d_B,
                acc::Complex *d_C,
                size_t    image_size);
|*
 * Multiplies many complex arrays A (in-place) by a single one B, pixel-by-pixel, after conjugating B
 *|
void batch_convol_B(int           blockIdx_x,
                    int           threadIdx_x,
                    acc::Complex  *d_A,
                    acc::Complex  *d_B,
                    size_t        image_size);
|*
 * Multiplies scalar array A by a scalar S
 *
 *  OUT[i] = A[i]*S
 *|
template <typename T>
void cpu_kernel_multi( T   *A,
            T   *OUT,
            T    S,
            size_t   image_size);
*/
/*
 * In place multiplies scalar array A by a scalar S
 *
 *  A[i] = A[i]*S
 */
template <typename T>
void cpu_kernel_multi(T *A, T S, size_t image_size);
/*
 * Multiplies scalar array A by scalar array B and a scalar S, pixel-by-pixel
 *
 *  OUT[i] = A[i]*B[i]*S
 */
template <typename T>
void cpu_kernel_multi(T *A, T *B, T *OUT, T S, size_t image_size);
/*
void finalizeMstddev(   int       blockIdx_x,
                        int       threadIdx_x,
                        XFLOAT   *Mstddev,
                        XFLOAT   *aux,
                        XFLOAT    S,
                        size_t       image_size);

|*
 * In place squares array in place
 *
 *  A[i] = A[i]*A[i]
 *|
void square(int       blockIdx_x,
            int       threadIdx_x,
            XFLOAT   *A,
            size_t       image_size);
*/
/*
 * Casts on device so we can copy_to_host directly into a multidimarray.
 *
template <typename T1, typename T2 >
void cast(  int blockIdx_x,
            int threadIdx_x,
            T1 *IN,
            T2 *OUT,
            size_t size)
{
    size_t pixel = (size_t)threadIdx_x + (size_t)blockIdx_x*(size_t)BLOCK_SIZE;
    if(pixel<size)
        OUT[pixel] = IN[pixel];
}
*/

template<bool do_highpass>
ALWAYS_INLINE_GCC void kernel_frequencyPass(
    size_t grid_size, size_t block_size,
    acc::Complex *A, long int ori_size,
    size_t Xdim, size_t Ydim, size_t Zdim,
    XFLOAT edge_low, XFLOAT edge_width, XFLOAT edge_high,
    XFLOAT angpix, size_t image_size
) {
    #ifdef DEBUG_CUDA
    if (grid_size * block_size > std::numeric_limits<int>::max())
        CHECK_INDEX_DEBUG_FATAL("kernel_frequencyPass:  grid_size*(size_t)block_size > (size_t)std::numeric_limits<int>::max()");
    #endif
    // TODO - why not a single loop over image_size pixels?
    for (size_t blk = 0; blk < grid_size; blk++)
    for (size_t tid = 0; tid < block_size; tid++) {
        size_t texel = tid + blk * block_size;

        int z = texel / (Xdim * Ydim);
        int xy = texel - z * Xdim * Ydim;
        int y = xy / Xdim;

        int xp = xy - y * Xdim;

        int zp = z < Xdim ? z : z - Zdim;
        int yp = y < Xdim ? y : y - Ydim;

        RFLOAT r2 = xp * xp + yp * yp + zp * zp;

        RFLOAT res;
        if (texel < image_size) {
            res = sqrt(r2) / ori_size;

            if (do_highpass) {
                if (res < edge_low) {
                    // highpass => lows are dead
                    A[texel].x = 0.0;
                    A[texel].y = 0.0;
                } else if (res < edge_high) {
                    // highpass => medium lows are almost dead
                    const XFLOAT factor = 0.5 - 0.5 * cos(PI * (res - edge_low) / edge_width);
                    A[texel].x *= factor;
                    A[texel].y *= factor;
                }
            } else {
                if (res > edge_high) {
                    // lowpass => highs are dead
                    A[texel].x = 0.0;
                    A[texel].y = 0.0;
                } else if (res > edge_low) {
                    // lowpass => medium highs are almost dead
                    const XFLOAT factor = 0.5 + 0.5 * cos(PI * (res - edge_low) / edge_width);
                    A[texel].x *= factor;
                    A[texel].y *= factor;
                }
            }
        }
    }
}

template<bool DATA3D>
ALWAYS_INLINE_GCC void powerClass(
    size_t gridSize,
    acc::Complex *g_image, XFLOAT *g_spectrum,
    size_t image_size, size_t spectrum_size,
    int xdim, int ydim, int zdim,
    int res_limit, XFLOAT *g_highres_Xi2
) {
    #ifdef DEBUG_CUDA
    if (gridSize * POWERCLASS_BLOCK_SIZE > std::numeric_limits<int>::max())
        CHECK_INDEX_DEBUG_FATAL("kernel_frequencyPass:  gridSize*(size_t)POWERCLASS_BLOCK_SIZE > (size_t)std::numeric_limits<int>::max()");
    #endif
    for (size_t bid = 0; bid < gridSize; bid++) {

        XFLOAT normFaux;
        int x, y, xy, d;
        int xydim = xdim * ydim;
        bool coords_in_range = true;

        XFLOAT s_highres_Xi2[POWERCLASS_BLOCK_SIZE];
        for (int tid = 0; tid < POWERCLASS_BLOCK_SIZE; tid++)
            s_highres_Xi2[tid] = (XFLOAT)0.;

        for (size_t tid = 0; tid < POWERCLASS_BLOCK_SIZE; tid++) {
            size_t voxel = tid + bid * POWERCLASS_BLOCK_SIZE;
            if (voxel >= image_size) continue;
            if (DATA3D) {
                int z =  voxel / xydim;
                xy = voxel % xydim;
                y =  xy / xdim;
                x =  xy % xdim;
                y = y < xdim ? y : y - ydim;
                z = z < xdim ? z : z - zdim;
                d  = x * x + y * y + z * z;
                coords_in_range = x != 0 || y >= 0.f || z >= 0.f;
            } else {
                x = voxel % xdim;
                y = (voxel - x) / xdim;
                y = y < xdim ? y : y - ydim;
                d  = x * x + y * y;
                coords_in_range = x != 0 || y >= 0.f;
            }

            #if defined(ACC_DOUBLE_PRECISION)
            size_t ires = sqrt(d) + 0.5;
            #else
            size_t ires = sqrtf(d) + 0.5f;
            #endif
            if (ires < spectrum_size && coords_in_range) {
                normFaux = g_image[voxel].x * g_image[voxel].x + g_image[voxel].y * g_image[voxel].y;
                g_spectrum[ires] += normFaux;
                if (ires>=res_limit)
                    s_highres_Xi2[tid] += normFaux;
            }
        }

        for (int tid = 1; tid < POWERCLASS_BLOCK_SIZE; tid++)
            s_highres_Xi2[0] += s_highres_Xi2[tid];

        g_highres_Xi2[0] += s_highres_Xi2[0];
    }
}

ALWAYS_INLINE_GCC void translatePixel(
    int x, int y, XFLOAT tx, XFLOAT ty,
    XFLOAT &real, XFLOAT &imag, XFLOAT &tReal, XFLOAT &tImag
) {
    XFLOAT s, c;
    ACC_SINCOS(x * tx + y * ty, &s, &c);
    tReal = c * real - s * imag;
    tImag = c * imag + s * real;
}

ALWAYS_INLINE_GCC void translatePixel(
    int x, int y, int z, XFLOAT tx, XFLOAT ty, XFLOAT tz,
    XFLOAT &real, XFLOAT &imag, XFLOAT &tReal, XFLOAT &tImag
) {
    XFLOAT s, c;
    ACC_SINCOS(x * tx + y * ty + z * tz, &s, &c);
    tReal = c * real - s * imag;
    tImag = c * imag + s * real;
}

// sincos lookup table optimization. Function translatePixel calls
// sincos(x*tx + y*ty). We precompute 2D lookup tables for x and y directions.
// The first dimension is x or y pixel index, and the second dimension is x or y
// translation index. Since sin(a+B) = sin(A) * cos(B) + cos(A) * sin(B), and
// cos(A+B) = cos(A) * cos(B) - sin(A) * sin(B), we can use lookup table to
// compute sin(x*tx + y*ty) and cos(x*tx + y*ty).
ALWAYS_INLINE_GCC void computeSincosLookupTable2D(
    unsigned long  trans_num,
    XFLOAT  *trans_x, XFLOAT  *trans_y,
    int      xSize, int      ySize,
    XFLOAT  *sin_x, XFLOAT  *cos_x,
    XFLOAT  *sin_y, XFLOAT  *cos_y
) {
    for (unsigned long i = 0; i < trans_num; i++) {
        const XFLOAT tx = trans_x[i];
        for (int x = 0; x < xSize; x++) {
           unsigned long index = i * xSize + x;
            ACC_SINCOS(x * tx, sin_x + index, cos_x + index);
        }
        const XFLOAT ty = trans_y[i];
        for (int y = 0; y < ySize; y++) {
            unsigned long index = i * ySize + y;
            ACC_SINCOS(y * ty, sin_y + index, cos_y + index);
        }
    }
}

ALWAYS_INLINE_GCC void computeSincosLookupTable3D(
    unsigned long  trans_num,
    XFLOAT  *trans_x, XFLOAT  *trans_y, XFLOAT  *trans_z,
    int      xSize, int      ySize, int      zSize,
    XFLOAT  *sin_x, XFLOAT  *cos_x,
    XFLOAT  *sin_y, XFLOAT  *cos_y,
    XFLOAT  *sin_z, XFLOAT  *cos_z
) {
    for (unsigned long i = 0; i < trans_num; i++) {
        const XFLOAT tx = trans_x[i];
        for (int x = 0; x < xSize; x++) {
           unsigned long index = i * xSize + x;
            ACC_SINCOS(x * tx, sin_x + index, cos_x + index);
        }
        const XFLOAT ty = trans_y[i];
        for (int y = 0; y < ySize; y++) {
            unsigned long index = i * ySize + y;
            ACC_SINCOS(y * ty, sin_y + index, cos_y + index);
        }
        const XFLOAT tz = trans_z[i];
        for (int z = 0; z < zSize; z++) {
            unsigned long index = i * zSize + z;
            ACC_SINCOS(z * tz, sin_z + index, cos_z + index);
        }
    }
}

template<bool invert>
void cpu_kernel_make_eulers_2D(
    size_t grid_size, size_t block_size,
    XFLOAT *alphas, XFLOAT *eulers,
    unsigned long orientation_num
);

template<bool invert,bool doL, bool doR>
void cpu_kernel_make_eulers_3D(
    size_t grid_size, size_t block_size,
    XFLOAT *alphas, XFLOAT *betas, XFLOAT *gammas, XFLOAT *eulers,
    unsigned long orientation_num,
    XFLOAT *L, XFLOAT *R
);

}

#endif /* HELPER_KERNELS_H_ */
