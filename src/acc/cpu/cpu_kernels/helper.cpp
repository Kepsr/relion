#include "src/acc/cpu/cuda_stubs.h"

#include "src/acc/acc_ptr.h"
#include "src/acc/acc_projector.h"
#include "src/acc/acc_backprojector.h"
#include "src/acc/acc_projector_plan.h"
#include "src/acc/cpu/cpu_benchmark_utils.h"
#include "src/acc/cpu/cpu_helper_functions.h"
#include "src/acc/cpu/cpu_kernels/helper.h"
#include "src/acc/cpu/cpu_kernels/diff2.h"
#include "src/acc/cpu/cpu_kernels/wavg.h"
#include "src/acc/cpu/cpu_kernels/BP.h"
#include "src/acc/utilities.h"
#include "src/acc/data_types.h"

#include "src/acc/acc_helper_functions.h"

#include "src/acc/cpu/cpu_kernels/cpu_utils.h"

namespace CpuKernels {

/*
 * This draft of a kernel assumes input that has jobs which have a single orientation and sequential translations within each job.
 *
 */
void exponentiate_weights_fine(
    XFLOAT *g_pdf_orientation, bool *g_pdf_orientation_zeros,
    XFLOAT *g_pdf_offset, bool *g_pdf_offset_zeros,
    XFLOAT *g_weights, XFLOAT min_diff2,
    unsigned long oversamples_orient, unsigned long oversamples_trans,
    unsigned long *d_rot_id, unsigned long *d_trans_idx,
    unsigned long *d_job_idx, unsigned long *d_job_num,
    long int job_num
) {
    for (long int jobid = 0; jobid < job_num; jobid++) {
        const long int pos = d_job_idx[jobid];
        // index of comparison
        const long int ix = d_rot_id   [pos];   // each thread gets its own orient...
              long int iy = d_trans_idx[pos];   // ...and its starting trans...
        const long int in = d_job_num  [jobid]; // ...AND the number of translations to go through
        for (long int itrans = 0; itrans < in; itrans++, iy++) {
            const int c_itrans = (iy - iy % oversamples_trans) / oversamples_trans;
            auto& w = g_weights[pos + itrans];
            w = w < min_diff2 || g_pdf_orientation_zeros[ix] || g_pdf_offset_zeros[c_itrans] ?
                std::numeric_limits<XFLOAT>::lowest() :  g_pdf_orientation[ix] + g_pdf_offset[c_itrans] + min_diff2 - w;
        }
    }
}

void RNDnormalDitributionComplexWithPowerModulation2D(
    acc::Complex* Image, size_t xdim, XFLOAT *spectrum
) {
    const size_t dim2 = (xdim - 1) * 2;
    const size_t size = xdim * dim2;
    for (size_t i = 0; i < size; i++) {
        size_t y = i / xdim;  // fftshift in one of two dims;
        size_t x = i % xdim;
        if (y >= xdim) y -= dim2;
        const int ires = sqrtf(x * x + y * y);
        if (ires < xdim) {
            Image[i].x = rnd_gaus(0.0, spectrum[ires]);
            Image[i].y = rnd_gaus(0.0, spectrum[ires]);
        } else {
            Image[i].x = 0.0;
            Image[i].y = 0.0;
        }
    }
}

void RNDnormalDitributionComplexWithPowerModulation3D(
    acc::Complex* Image, size_t xdim, size_t ydim, XFLOAT *spectrum
) {
    const int xydim = xdim * ydim;
    const int dim2 = (xdim - 1) * 2;  // assuming square input images (particles)
    const int size = xdim * dim2;
    for (int i = 0; i < size; i++) {
        int z = i / xydim;
        int y = i - z / ydim;
        int x = i % xdim;
        // fftshift in two of three dims;
        if (z >= xdim) z -= dim2;					
        if (y >= xdim) y -= dim2;
        const int ires = sqrtf(x * x + y * y + z * z);
        if (ires < xdim) {
            Image[i].x = rnd_gaus(0.0, spectrum[ires]);
            Image[i].y = rnd_gaus(0.0, spectrum[ires]);
        } else {
            Image[i].x = 0;
            Image[i].y = 0;
        }
    }
}

void softMaskBackgroundValue(
    size_t block_dim, size_t block_size,
    XFLOAT* vol, size_t vol_size,
    long int xdim, long int ydim, long int zdim,
    long int xinit, long int yinit, long int zinit,
    XFLOAT radius, XFLOAT radius_p, XFLOAT cosine_width,
    XFLOAT* g_sum, XFLOAT* g_sum_bg
) {
    const int xydim = xdim * ydim;
    const size_t n_passes = ceilfracf(vol_size, block_size * block_dim);
    for (size_t bid = 0; bid < block_dim;  bid++)
    for (size_t tid = 0; tid < block_size; tid++) {
        size_t i = bid * block_size * n_passes + tid;
        for (size_t pass = 0; pass < n_passes; pass++, i += block_size) {
            // loop through all translations for this orientation
            if (i >= vol_size) break;
            const XFLOAT img_pixels = vol[i];
            const int z =  i / xydim         - zinit;
            const int y = (i % xydim) / xdim - yinit;
            const int x = (i % xydim) % xdim - xinit;
            const XFLOAT r = hypot((double) x, y, z);
            if (r < radius) {
                continue;
            } else if (r > radius_p) {
                g_sum[tid]    += 1.0;
                g_sum_bg[tid] += img_pixels;
            } else {
                const XFLOAT factor = raised_cos((radius_p - r) / cosine_width * M_PI);
                g_sum[tid]    += factor;
                g_sum_bg[tid] += factor * img_pixels;
            }
        }
    }
}

void cosineFilter(
    size_t block_dim, size_t block_size,
    XFLOAT *vol, size_t vol_size,
    long int xdim, long int ydim, long int zdim,
    long int xinit, long int yinit, long int zinit,
    bool do_noise, XFLOAT *noise,
    XFLOAT radius, XFLOAT radius_p,
    XFLOAT cosine_width, XFLOAT bg_value
) {
    const size_t n_passes = ceilfracf(vol_size, block_size * block_dim);
    for (size_t bid = 0; bid < block_dim;  bid++)
    for (size_t tid = 0; tid < block_size; tid++) {
        size_t i = bid * block_size * n_passes + tid;
        XFLOAT defVal = bg_value;
        for (size_t pass = 0; pass < n_passes; pass++, i += block_size) {
            // loop the available warps enough to complete all translations for this orientation
            if (i >= vol_size) break;
            const int xydim = xdim * ydim;
            const int z =  i / xydim          - zinit;
            const int y = (i % xydim) / xdim  - yinit;
            const int x = (i % xydim) % xdim  - xinit;
            const XFLOAT r = hypot((double) x, y, z);
            if (do_noise)
                defVal = noise[i];
            if (r < radius) {
                continue;
            } else if (r > radius_p) {
                vol[i] = defVal;
            } else {
                const XFLOAT factor = raised_cos((radius_p - r) / cosine_width * M_PI);
                vol[i] = vol[i] * (1 - factor) + defVal * factor;
            }
        }
    }
}

template <typename T>
ALWAYS_INLINE_GCC void cpu_translate2D(
    T *g_image_in, T *g_image_out,
    size_t image_size, int xdim, int ydim, int dx, int dy
) {
    #ifdef DEBUG_CUDA
    if (image_size > std::numeric_limits<int>::max())
        AccPtr<T>::HandleDebugInformational("cpu_translate2D: image_size > std::numeric_limits<int>::max()", __FILE__, __LINE__);
    #endif
    for (size_t i = 0; i < image_size; i++) {
        const int x = i % xdim;
        const int y = (i - x) / xdim;
        const int xp = x + dx;
        const int yp = y + dy;
        if (yp < 0 || xp < 0 || yp >= ydim || xp >= xdim) continue;
        const size_t j = yp * xdim + xp;
        if (j >= 0 && j < image_size)  // if displacement is negative, j could be less than 0
            g_image_out[j] = g_image_in[i];
    }
}

template <typename T>
ALWAYS_INLINE_GCC void cpu_translate3D(
    T * g_image_in, T*  g_image_out, size_t image_size,
    int  xdim, int  ydim, int  zdim,
    int  dx, int  dy, int  dz
) {
    #ifdef DEBUG_CUDA
    if (image_size > std::numeric_limits<int>::max())
        AccPtr<T>::HandleDebugInformational("cpu_translate3D: image_size > std::numeric_limits<int>::max()", __FILE__, __LINE__);
    #endif
    const int xydim = xdim * ydim;
    for (size_t i = 0; i < image_size; i++) {
        const int z =  i / xydim;
        const int zp = z + dz;
        const int xy = i % xydim;
        const int y =  xy / xdim;
        const int yp = y + dy;
        const int x =  xy % xdim;
        const int xp = x + dx;
        if (zp < 0 || yp < 0 || xp < 0 || zp >= zdim || yp >= ydim || xp >= xdim) continue;
        const size_t j = zp * xydim + yp * xdim + xp;
        if (j >= 0 && j < image_size) // if displacement is negative, j could be less than 0
            g_image_out[j] = g_image_in[i];
    }
}

template <typename T>
ALWAYS_INLINE_GCC void centerFFT_2D(
    size_t batch_size, size_t pixel_start, size_t pixel_end,
    T *img_in, size_t image_size,
    int xdim, int ydim, int xshift, int yshift
) {
    #ifdef DEBUG_CUDA
    if (image_size > std::numeric_limits<int>::max())
        AccPtr<T>::HandleDebugInformational("centerFFT_2D: image_size > std::numeric_limits<int>::max()", __FILE__, __LINE__);
    if (image_size * batch_size > std::numeric_limits<int>::max())
        AccPtr<T>::HandleDebugInformational("centerFFT_2D: image_size*batch_size > std::numeric_limits<int>::max()", __FILE__, __LINE__);
    if (pixel_end > image_size)
        AccPtr<T>::HandleDebugInformational("centerFFT_2D: pixel_end > image_size", __FILE__, __LINE__);
    #endif
    for (int batch = 0; batch < batch_size; batch++) {
        const size_t image_offset = image_size * batch;
    for (size_t i = pixel_start; i < pixel_end; i++) {
        const int y = floorf((XFLOAT) i / (XFLOAT) xdim);
        const int x = i % xdim;    // also = i - y*xdim, but this depends on y having been calculated, i.e. serial evaluation
        const int yp = (y + yshift + ydim) % ydim;
        const int xp = (x + xshift + xdim) % xdim;
        const size_t j = yp * xdim + xp;
        std::swap(img_in[image_offset + i], img_in[image_offset + j]);
    }}
}

template void centerFFT_2D<float>
(size_t batch_size, size_t pixel_start, size_t pixel_end,
float *img_in, size_t image_size,
int xdim, int ydim, int xshift, int yshift);

template void centerFFT_2D<double>
(size_t batch_size, size_t pixel_start, size_t pixel_end,
double *img_in, size_t image_size,
int xdim, int ydim, int xshift, int yshift);

template <typename T>
ALWAYS_INLINE_GCC void centerFFT_3D(
    size_t batch_size, size_t pixel_start, size_t pixel_end,
    T *img_in, size_t image_size,
    int xdim, int ydim, int zdim,
    int xshift, int yshift, int zshift
) {
    #ifdef DEBUG_CUDA
    if (image_size > std::numeric_limits<int>::max())
        AccPtr<T>::HandleDebugInformational("centerFFT_3D: image_size > std::numeric_limits<int>::max()", __FILE__, __LINE__);
    if (image_size * batch_size > std::numeric_limits<int>::max())
        AccPtr<T>::HandleDebugInformational("centerFFT_3D: image_size*batch_size > std::numeric_limits<int>::max()", __FILE__, __LINE__);
    if (pixel_end > image_size)
        AccPtr<T>::HandleDebugInformational("centerFFT_3D: pixel_end > image_size", __FILE__, __LINE__);
    #endif
    const int xydim = xdim * ydim;
    for (int batch = 0; batch < batch_size; batch++) {
        const size_t image_offset = image_size * batch;
    for (size_t i = pixel_start; i < pixel_end; i++) {
        const int z = floorf((XFLOAT) i / (XFLOAT) xydim);
        const int xy = i % xydim;
        const int y = floorf((XFLOAT) xy / (XFLOAT) xdim);
        const int x = xy % xdim;
        const int xp = (x + xshift + xdim) % xdim;
        const int yp = (y + yshift + ydim) % ydim;
        const int zp = (z + zshift + zdim) % zdim;
        const size_t j = zp * xydim + yp * xdim + xp;
        std::swap(img_in[image_offset + i], img_in[image_offset + j]);
    }}
}

template void centerFFT_3D<float>
(size_t batch_size, size_t pixel_start, size_t pixel_end,
float *img_in, size_t image_size,
int xdim, int ydim, int zdim,
int xshift, int yshift, int zshift);

template void centerFFT_3D<double>
(size_t batch_size, size_t pixel_start, size_t pixel_end,
double *img_in, size_t image_size,
int xdim, int ydim, int zdim,
int xshift, int yshift, int zshift);

/* TODO - if create optimized CPU version of autopicker
 * All these functions need to be converted to use internal loops rather than
 * block and thread indices to operate like other active functions seen in this file
void probRatio( int       blockIdx_x,
                int       threadIdx_x,
                XFLOAT   *d_Mccf,
                XFLOAT   *d_Mpsi,
                XFLOAT   *d_Maux,
                XFLOAT   *d_Mmean,
                XFLOAT   *d_Mstddev,
                size_t       image_size,
                XFLOAT    normfft,
                XFLOAT    sum_ref_under_circ_mask,
                XFLOAT    sum_ref2_under_circ_mask,
                XFLOAT    expected_Pratio,
                int       NpsiThisBatch,
                int       startPsi,
                int       totalPsis)
{
    |* PLAN TO:
     *
     * 1) Pre-filter
     * 		d_Mstddev[i] = 1 / (2*d_Mstddev[i])   ( if d_Mstddev[pixel] > 1E-10 )
     * 		d_Mstddev[i] = 1    				  ( else )
     *
     * 2) Set
     * 		sum_ref2_under_circ_mask /= 2.
     *
     * 3) Total expression becomes
     * 		diff2 = ( exp(k) - 1.f ) / (expected_Pratio - 1.f)
     * 	  where
     * 	  	k = (normfft * d_Maux[pixel] + d_Mmean[pixel] * sum_ref_under_circ_mask)*d_Mstddev[i] + sum_ref2_under_circ_mask
     *
     *|

    size_t pixel = (size_t)threadIdx_x + (size_t)blockIdx_x*(size_t)PROBRATIO_BLOCK_SIZE;
    if(pixel<image_size)
    {
        XFLOAT Kccf = d_Mccf[pixel];
        XFLOAT Kpsi =(XFLOAT)-1.0;
        for(int psi = 0; psi < NpsiThisBatch; psi++ )
        {
            XFLOAT diff2 = normfft * d_Maux[pixel + image_size*psi];
            diff2 += d_Mmean[pixel] * sum_ref_under_circ_mask;

    //		if (d_Mstddev[pixel] > (XFLOAT)1E-10)
            diff2 *= d_Mstddev[pixel];
            diff2 += sum_ref2_under_circ_mask;

            #if defined(ACC_DOUBLE_PRECISION)
            diff2 = exp(-diff2 / 2.); // exponentiate to reflect the Gaussian error model. sigma=1 after normalization, 0.4=1/sqrt(2pi)
            #else
            diff2 = expf(-diff2 / 2.f);
            #endif

            // Store fraction of (1 - probability-ratio) wrt  (1 - expected Pratio)
            diff2 = (diff2 - (XFLOAT)1.0) / (expected_Pratio - (XFLOAT)1.0);
            if (diff2 > Kccf)
            {
                Kccf = diff2;
                Kpsi = (startPsi + psi)*(360/totalPsis);
            }
        }
        d_Mccf[pixel] = Kccf;
        if (Kpsi >= 0.)
            d_Mpsi[pixel] = Kpsi;
    }
}

void rotateOnly(int              blockIdx_x,
                int              blockIdx_y,
                int              threadIdx_x,
                acc::Complex     *d_Faux,
                XFLOAT           psi,
                AccProjectorKernel &projector,
                int              startPsi
               )
{
    int proj = blockIdx_y;
    size_t image_size=(size_t)projector.imgX*(size_t)projector.imgY;
    size_t pixel = (size_t)threadIdx_x + (size_t)blockIdx_x*(size_t)BLOCK_SIZE;
    if(pixel<image_size)
    {
        int y = floorfracf(pixel,projector.imgX);
        int x = pixel % projector.imgX;

        if (y > projector.maxR)
        {
            if (y >= projector.imgY - projector.maxR)
                y = y - projector.imgY;
            else
                x = projector.maxR;
        }

        XFLOAT sa, ca;
        ACC_SINCOS((proj+startPsi)*psi, &sa, &ca);

        acc::Complex val;

        projector.project2Dmodel(	 x,y,
                                     ca,
                                    -sa,
                                     sa,
                                     ca,
                                     val.x,val.y);

        long int out_pixel = proj*image_size + pixel;

        d_Faux[out_pixel].x =val.x;
        d_Faux[out_pixel].y =val.y;
    }
}

void rotateAndCtf(  int              blockIdx_x,
                    int              blockIdx_y,
                    int              threadIdx_x,
                    acc::Complex     *d_Faux,
                    XFLOAT          *d_ctf,
                    XFLOAT           psi,
                    AccProjectorKernel &projector,
                    int       startPsi
                )
{
    int proj = blockIdx_y;
    size_t image_size=(size_t)projector.imgX*(size_t)projector.imgY;
    size_t pixel = (size_t)threadIdx_x + (size_t)blockIdx_x*(size_t)BLOCK_SIZE;
    if(pixel<image_size)
    {
        int y = floorfracf(pixel,projector.imgX);
        int x = pixel % projector.imgX;

        if (y > projector.maxR)
        {
            if (y >= projector.imgY - projector.maxR)
                y = y - projector.imgY;
            else
                x = projector.maxR;
        }

        XFLOAT sa, ca;
        ACC_SINCOS((proj+startPsi)*psi, &sa, &ca);

        acc::Complex val;

        projector.project2Dmodel(	 x,y,
                                     ca,
                                    -sa,
                                     sa,
                                     ca,
                                     val.x,val.y);

        long int out_pixel = proj*image_size + pixel;

        d_Faux[out_pixel].x =val.x*d_ctf[pixel];
        d_Faux[out_pixel].y =val.y*d_ctf[pixel];

    }
}


void convol_A(  int           blockIdx_x,
                int           threadIdx_x,
                acc::Complex  *d_A,
                acc::Complex  *d_B,
                size_t           image_size)
{
    size_t pixel = (size_t)threadIdx_x + (size_t)blockIdx_x*(size_t)BLOCK_SIZE;
    if(pixel<image_size)
    {
        XFLOAT tr =   d_A[pixel].x;
        XFLOAT ti = - d_A[pixel].y;
        d_A[pixel].x =   tr*d_B[pixel].x - ti*d_B[pixel].y;
        d_A[pixel].y =   ti*d_B[pixel].x + tr*d_B[pixel].y;
    }
}

void convol_A(  int          blockIdx_x,
                int          threadIdx_x,
                acc::Complex *d_A,
                acc::Complex *d_B,
                acc::Complex *d_C,
                size_t          image_size)
{
    size_t pixel = (size_t)threadIdx_x + (size_t)blockIdx_x*(size_t)BLOCK_SIZE;
    if(pixel<image_size)
    {
        XFLOAT tr =   d_A[pixel].x;
        XFLOAT ti = - d_A[pixel].y;
        d_C[pixel].x =   tr*d_B[pixel].x - ti*d_B[pixel].y;
        d_C[pixel].y =   ti*d_B[pixel].x + tr*d_B[pixel].y;
    }
}

void batch_convol_A(int           blockIdx_x,
                    int           blockIdx_y,
                    int           threadIdx_x,
                    acc::Complex  *d_A,
                    acc::Complex  *d_B,
                    size_t           image_size)
{
    size_t pixel = (size_t)threadIdx_x + (size_t)blockIdx_x*(size_t)BLOCK_SIZE;
    int A_off = blockIdx_y * image_size;
    if(pixel<image_size)
    {
        XFLOAT tr =   d_A[pixel + A_off].x;
        XFLOAT ti = - d_A[pixel + A_off].y;
        d_A[pixel + A_off].x =   tr*d_B[pixel].x - ti*d_B[pixel].y;
        d_A[pixel + A_off].y =   ti*d_B[pixel].x + tr*d_B[pixel].y;
    }
}

void batch_convol_A(int           blockIdx_x,
                    int           blockIdx_y,
                    int           threadIdx_x,
                    acc::Complex  *d_A,
                    acc::Complex  *d_B,
                    acc::Complex  *d_C,
                    size_t           image_size)
{
    size_t pixel = (size_t)threadIdx_x + (size_t)blockIdx_x*(size_t)BLOCK_SIZE;
    int A_off = blockIdx_y*image_size;
    if(pixel<image_size)
    {
        XFLOAT tr =   d_A[pixel + A_off].x;
        XFLOAT ti = - d_A[pixel + A_off].y;
        d_C[pixel + A_off].x =   tr*d_B[pixel].x - ti*d_B[pixel].y;
        d_C[pixel + A_off].y =   ti*d_B[pixel].x + tr*d_B[pixel].y;
    }
}

void convol_B(  int          blockIdx_x,
                int          threadIdx_x,
                acc::Complex *d_A,
                acc::Complex *d_B,
                size_t          image_size)
{
    size_t pixel =  (size_t)threadIdx_x + (size_t)blockIdx_x*(size_t)BLOCK_SIZE;
    if(pixel<image_size)
    {
        XFLOAT tr = d_A[pixel].x;
        XFLOAT ti = d_A[pixel].y;
        d_A[pixel].x =   tr*d_B[pixel].x + ti*d_B[pixel].y;
        d_A[pixel].y =   ti*d_B[pixel].x - tr*d_B[pixel].y;
    }
}

void convol_B(  int           blockIdx_x,
                int           threadIdx_x,
                acc::Complex  *d_A,
                acc::Complex  *d_B,
                acc::Complex  *d_C,
                size_t           image_size)
{
    size_t pixel = (size_t)threadIdx_x + (size_t)blockIdx_x*(size_t)BLOCK_SIZE;
    if(pixel<image_size)
    {
        XFLOAT tr = d_A[pixel].x;
        XFLOAT ti = d_A[pixel].y;
        d_C[pixel].x =   tr*d_B[pixel].x + ti*d_B[pixel].y;
        d_C[pixel].y =   ti*d_B[pixel].x - tr*d_B[pixel].y;
    }
}

void batch_convol_B(int          blockIdx_x,
                    int          blockIdx_y,
                    int          threadIdx_x,
                    acc::Complex *d_A,
                    acc::Complex *d_B,
                    size_t          image_size)
{
    long int pixel = (size_t)threadIdx_x + (size_t)blockIdx_x*(size_t)BLOCK_SIZE;
    int A_off = blockIdx_y*image_size;
    if(pixel<image_size)
    {
        XFLOAT tr = d_A[pixel + A_off].x;
        XFLOAT ti = d_A[pixel + A_off].y;
        d_A[pixel + A_off].x =   tr*d_B[pixel].x + ti*d_B[pixel].y;
        d_A[pixel + A_off].y =   ti*d_B[pixel].x - tr*d_B[pixel].y;
    }
}
*/
template <typename T>
void cpu_kernel_multi(T *A, T *OUT, T  S, size_t     image_size) {
    #ifdef DEBUG_CUDA
    if (image_size < 0)
        AccPtr<T>::HandleDebugInformational("cpu_kernel_multi:  image_size < 0", __FILE__, __LINE__);
    #endif
    for (size_t i = 0; i < image_size; i ++)
        OUT[i] = A[i] * S;
}

template <typename T>
void cpu_kernel_multi( T *A, T  S, size_t     image_size) {
    #ifdef DEBUG_CUDA
    if (image_size < 0)
        AccPtr<T>::HandleDebugInformational("cpu_kernel_multi2:  image_size < 0", __FILE__, __LINE__);
    #endif
    for (size_t i = 0; i < image_size; i ++)
        A[i] *= S;
}

template <typename T>
void cpu_kernel_multi(
    T *A, T *B, T *OUT, T  S, size_t     image_size
) {
    #ifdef DEBUG_CUDA
    if (image_size < 0)
        AccPtr<T>::HandleDebugInformational("cpu_kernel_multi3:  image_size < 0", __FILE__, __LINE__);
    #endif
    for (size_t i = 0; i < image_size; i ++)
        OUT[i] = A[i] * B[i] * S;
}
/*
void batch_multi(   int     blockIdx_x,
                    int     blockIdx_y,
                    int     threadIdx_x,
                    XFLOAT *A,
                    XFLOAT *B,
                    XFLOAT *OUT,
                    XFLOAT  S,
                    size_t     image_size)
{
    sise_t pixel = (size_t)threadIdx_x + (size_t)blockIdx_x*(size_t)BLOCK_SIZE;
    if(pixel<image_size)
        OUT[pixel + blockIdx_y*image_size] = A[pixel + blockIdx_y*image_size]*B[pixel + blockIdx_y*image_size]*S;
}
 */
/* TODO - CPU-optimized autopicker
void finalizeMstddev(   int     blockIdx_x,
                        int     threadIdx_x,
                        XFLOAT *Mstddev,
                        XFLOAT *aux,
                        XFLOAT  S,
                        size_t     image_size)
{
    int size_t = (size_t)threadIdx_x + (size_t)blockIdx_x*(size_t)BLOCK_SIZE;
    if(pixel<image_size)
    {
        XFLOAT temp = Mstddev[pixel] + S * aux[pixel];
        if(temp > 0)
            Mstddev[pixel] = sqrt(temp);
        else
            Mstddev[pixel] = 0;
    }
}

void square(int     blockIdx_x,
            int     threadIdx_x,
            XFLOAT *A,
            size_t     image_size)
{
    size_t pixel = (size_t)threadIdx_x + (size_t)blockIdx_x*(size_t)BLOCK_SIZE;
    if(pixel<image_size)
        A[pixel] = A[pixel]*A[pixel];
}
*/

template<bool invert>
ALWAYS_INLINE_GCC void cpu_kernel_make_eulers_2D(
    size_t grid_size, size_t block_size, XFLOAT *alphas, XFLOAT *eulers, unsigned long orientation_num
) {
    #ifdef DEBUG_CUDA
    if (grid_size * block_size > std::numeric_limits<int>::max())
        AccPtr<XFLOAT>::HandleDebugInformational("cpu_kernel_make_eulers_2D: grid_size*block_size > std::numeric_limits<int>::max()", __FILE__, __LINE__);
    #endif
    for (size_t blockIdx_x = 0;  blockIdx_x  < grid_size;  blockIdx_x++)
    for (size_t threadIdx_x = 0; threadIdx_x < block_size; threadIdx_x++) {
        const unsigned long i = blockIdx_x * block_size + threadIdx_x;  // Orientation id
        if (i >= orientation_num) break;
        const XFLOAT a = radians(invert ? -alphas[i] : alphas[i]);
        XFLOAT ca, sa;
        ACC_SINCOS(a, &sa, &ca);

        eulers[9 * i + 0] =  ca;
        eulers[9 * i + 1] = -sa;
        eulers[9 * i + 2] =   0;
        eulers[9 * i + 3] =  sa;
        eulers[9 * i + 4] =  ca;
        eulers[9 * i + 5] =   0;
        eulers[9 * i + 6] =   0;
        eulers[9 * i + 7] =   0;
        eulers[9 * i + 8] =   1;
    }
}

template<bool invert, bool doL, bool doR>
ALWAYS_INLINE_GCC void cpu_kernel_make_eulers_3D(
    size_t grid_size, size_t block_size,
    XFLOAT *alphas, XFLOAT *betas, XFLOAT *gammas, XFLOAT *eulers,
    unsigned long orientation_num,
    XFLOAT *L, XFLOAT *R
) {
    #ifdef DEBUG_CUDA
    if (grid_size * block_size > std::numeric_limits<int>::max())
        AccPtr<XFLOAT>::HandleDebugInformational("cpu_kernel_make_eulers_3D: grid_size*block_size > std::numeric_limits<int>::max()", __FILE__, __LINE__);
    #endif
    for (size_t blockIdx_x  = 0; blockIdx_x  < grid_size;  blockIdx_x++)
    for (size_t threadIdx_x = 0; threadIdx_x < block_size; threadIdx_x++) {

        const unsigned long i = blockIdx_x * block_size + threadIdx_x;  // Orientation id

        if (i >= orientation_num) return;

        const XFLOAT a = radians(alphas[i]);
        const XFLOAT b = radians(betas [i]);
        const XFLOAT g = radians(gammas[i]);

        XFLOAT ca, sa, cb, sb, cg, sg;
        ACC_SINCOS(a, &sa, &ca);
        ACC_SINCOS(b, &sb, &cb);
        ACC_SINCOS(g, &sg, &cg);

        const XFLOAT cc = cb * ca;
        const XFLOAT cs = cb * sa;
        const XFLOAT sc = sb * ca;
        const XFLOAT ss = sb * sa;

        XFLOAT A[9], B[9];

        A[0] = ( cg * cc - sg * sa);//00
        A[1] = ( cg * cs + sg * ca);//01
        A[2] = (-cg * sb )         ;//02
        A[3] = (-sg * cc - cg * sa);//10
        A[4] = (-sg * cs + cg * ca);//11
        A[5] = ( sg * sb )         ;//12
        A[6] = ( sc )              ;//20
        A[7] = ( ss )              ;//21
        A[8] = ( cb )              ;//22

        if (doR) {
            std::fill_n(B, 9, 0.f);
            for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++)
            for (int k = 0; k < 3; k++)
                B[i * 3 + j] += A[i * 3 + k] * R[k * 3 + j];
        } else {
            std::copy_n(A, 9, B);
        }

        if (doL) {
            if (doR)
                std::copy_n(B, 9, A);

            std::fill_n(B, 9, 0.f);
            for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++)
            for (int k = 0; k < 3; k++)
                B[i * 3 + j] += L[i * 3 + k] * A[k * 3 + j];
        }

        if (invert) {
            if (doL) {
                // this could have anisotropy, so inverse != transpose!!!
                const XFLOAT det
                        = B[0] * (B[4] * B[8] - B[7] * B[5])
                        - B[1] * (B[3] * B[8] - B[6] * B[5])
                        + B[2] * (B[3] * B[7] - B[6] * B[4]);

                eulers[9 * i + 0] = (B[4] * B[8] - B[7] * B[5]) / det;
                eulers[9 * i + 1] = (B[7] * B[2] - B[1] * B[8]) / det;
                eulers[9 * i + 2] = (B[1] * B[5] - B[4] * B[2]) / det;
                eulers[9 * i + 3] = (B[5] * B[6] - B[8] * B[3]) / det;
                eulers[9 * i + 4] = (B[8] * B[0] - B[2] * B[6]) / det;
                eulers[9 * i + 5] = (B[2] * B[3] - B[5] * B[0]) / det;
                eulers[9 * i + 6] = (B[3] * B[7] - B[6] * B[4]) / det;
                eulers[9 * i + 7] = (B[6] * B[1] - B[0] * B[7]) / det;
                eulers[9 * i + 8] = (B[0] * B[4] - B[3] * B[1]) / det;
            } else {
                eulers[9 * i + 0] = B[0];//00
                eulers[9 * i + 1] = B[3];//01
                eulers[9 * i + 2] = B[6];//02
                eulers[9 * i + 3] = B[1];//10
                eulers[9 * i + 4] = B[4];//11
                eulers[9 * i + 5] = B[7];//12
                eulers[9 * i + 6] = B[2];//20
                eulers[9 * i + 7] = B[5];//21
                eulers[9 * i + 8] = B[8];//22
            }
        } else {
            std::copy_n(B, 9, eulers + 9 * i);
        }
    }
}

}

/// Some explicit template instantiations
template void CpuKernels::cpu_translate2D<XFLOAT>
(XFLOAT*, XFLOAT*, size_t, int, int, int, int);

template void CpuKernels::cpu_translate3D<XFLOAT>
(XFLOAT*, XFLOAT*, size_t, int, int, int, int, int, int);

template void CpuKernels::cpu_kernel_multi<XFLOAT>
(XFLOAT*, XFLOAT, size_t);

#define INSTANTIATE( A, B, C ) \
    template void CpuKernels::cpu_kernel_make_eulers_3D<A, B, C> \
    (size_t, size_t, XFLOAT*, XFLOAT*, XFLOAT*, XFLOAT*, unsigned long, XFLOAT*, XFLOAT*);

INSTANTIATE(0, 0, 0)
INSTANTIATE(0, 0, 1)
INSTANTIATE(0, 1, 0)
INSTANTIATE(0, 1, 1)
INSTANTIATE(1, 0, 0)
INSTANTIATE(1, 0, 1)
INSTANTIATE(1, 1, 0)
INSTANTIATE(1, 1, 1)

#undef INSTANTIATE

template void CpuKernels::cpu_kernel_make_eulers_2D<true>
(size_t, size_t, XFLOAT*, XFLOAT*, unsigned long);

template void CpuKernels::cpu_kernel_make_eulers_2D<false>
(size_t, size_t, XFLOAT*, XFLOAT*, unsigned long);
