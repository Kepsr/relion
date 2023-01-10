#include "src/acc/settings.h"
#include "src/acc/cuda/cuda_kernels/cuda_device_utils.cuh"
#include "src/acc/cuda/cuda_kernels/helper.cuh"
#include "src/acc/cuda/cuda_settings.h"

#include <curand.h>
#include <curand_kernel.h>
#include "src/acc/cuda/cub/cub.cuh"

// Macros, because we can't call __host__ functions from __global__ functions.
#if defined(ACC_DOUBLE_PRECISION)
#define RAISED_COSPI(x) (0.5  + 0.5  * cospi(x))
#else
#define RAISED_COSPI(x) (0.5f + 0.5f * cospif(x))
#endif

// It may surprise you,
// but floating-point division is faster than integer division.
__device__ double floordiv(double x, double y) {
    return floor(x / y);
}

namespace device {

    __device__ inline double hypot(int x, int y, int z) {
        return sqrt(double(x * x + y * y + z * z));
    }

    __device__ inline double hypot(int x, int y) {
        return sqrt(double(x * x + y * y));
    }

    __device__ constexpr inline float radians(float theta) {
        return theta * (float) PI / (float) 180.0;
    }

    __device__ constexpr inline double radians(double theta) {
        return theta * (double) PI / (double) 180.0;
    }

}

/// Needed explicit template instantiations
template __global__ void cuda_kernel_make_eulers_2D<true> (XFLOAT*, XFLOAT*, unsigned);
template __global__ void cuda_kernel_make_eulers_2D<false>(XFLOAT*, XFLOAT*, unsigned);

template __global__ void cuda_kernel_make_eulers_3D<true,  true,  true> (XFLOAT*, XFLOAT*, XFLOAT *, XFLOAT*, unsigned, XFLOAT*, XFLOAT*);
template __global__ void cuda_kernel_make_eulers_3D<true,  true,  false>(XFLOAT*, XFLOAT*, XFLOAT *, XFLOAT*, unsigned, XFLOAT*, XFLOAT*);
template __global__ void cuda_kernel_make_eulers_3D<true,  false, true> (XFLOAT*, XFLOAT*, XFLOAT *, XFLOAT*, unsigned, XFLOAT*, XFLOAT*);
template __global__ void cuda_kernel_make_eulers_3D<true,  false, false>(XFLOAT*, XFLOAT*, XFLOAT *, XFLOAT*, unsigned, XFLOAT*, XFLOAT*);
template __global__ void cuda_kernel_make_eulers_3D<false, true,  true> (XFLOAT*, XFLOAT*, XFLOAT *, XFLOAT*, unsigned, XFLOAT*, XFLOAT*);
template __global__ void cuda_kernel_make_eulers_3D<false, true,  false>(XFLOAT*, XFLOAT*, XFLOAT *, XFLOAT*, unsigned, XFLOAT*, XFLOAT*);
template __global__ void cuda_kernel_make_eulers_3D<false, false, true> (XFLOAT*, XFLOAT*, XFLOAT *, XFLOAT*, unsigned, XFLOAT*, XFLOAT*);
template __global__ void cuda_kernel_make_eulers_3D<false, false, false>(XFLOAT*, XFLOAT*, XFLOAT *, XFLOAT*, unsigned, XFLOAT*, XFLOAT*);

/*
 * This draft of a kernel assumes input that has jobs which have a single orientation and sequential translations within each job.
 *
 */
__global__ void cuda_kernel_exponentiate_weights_fine(
    XFLOAT *g_pdf_orientation, bool *g_pdf_orientation_zeros,
    XFLOAT *g_pdf_offset, bool *g_pdf_offset_zeros,
    XFLOAT *g_weights,
    XFLOAT min_diff2,
    int oversamples_orient, int oversamples_trans,
    unsigned long * d_rot_id, unsigned long * d_trans_idx,
    unsigned long * d_job_idx, unsigned long * d_job_num,
    long int job_num
) {
    const long int i = blockIdx.x * SUMW_BLOCK_SIZE + threadIdx.x;
    if (i >= job_num) return;
    const long int j = d_job_idx[i];  // index of comparison
    const long int in = d_job_num[i]; // Number of translations to go through
    const long int ix = d_rot_id[j];   // Orient
          long int iy = d_trans_idx[j];  // Starting trans
    for (int itrans = 0; itrans < in; itrans++, iy++) {
        const int c_itrans = (iy - iy % oversamples_trans) / oversamples_trans;

        auto& w = g_weights[j + itrans];
        if (w < min_diff2 || g_pdf_orientation_zeros[ix] || g_pdf_offset_zeros[c_itrans])
            w = -99e99;  // large negative number
        else
            w = g_pdf_orientation[ix] + g_pdf_offset[c_itrans] + min_diff2 - w;
    }
}

__global__ void cuda_kernel_initRND(unsigned long seed, curandState *States) {
    const int i = blockIdx.x * RND_BLOCK_SIZE + threadIdx.x;
    curand_init(seed, i, 0, States + i);
}

#if defined(ACC_DOUBLE_PRECISION)
#define CURAND_NORMAL2( x ) curand_normal2_double(x)
#else
#define CURAND_NORMAL2( x ) curand_normal2(x)
#endif

__global__ void cuda_kernel_RNDnormalDitributionComplexWithPowerModulation2D(
    acc::Complex* Image, curandState *States, long int xdim, XFLOAT* spectra
) {
    const int i = blockIdx.x * RND_BLOCK_SIZE + threadIdx.x;
    // curand_init(1234, pixel, 0, States + i);
    const int dim2 = (xdim - 1) * 2;  // Assuming square images (particles)
    const int size = xdim * dim2;
    const int passes = size / (RND_BLOCK_NUM * RND_BLOCK_SIZE) + 1;
    for (int pass = 0, j = i; pass != passes && j != size; pass++, j += RND_BLOCK_NUM * RND_BLOCK_SIZE) {
        const int x = j % xdim;
              int y = j / xdim;
        // fftshift
        if (y >= xdim) y -= dim2;
        const int ires = rintf(hypotf(x, y));
        const XFLOAT scale = ires < xdim ? spectra[ires] : 0.0;
        Image[j] = CURAND_NORMAL2(States + i) * scale;
    }
}

__global__ void cuda_kernel_RNDnormalDitributionComplexWithPowerModulation3D(
    acc::Complex* Image, curandState *States, long int xdim, long int ydim, XFLOAT* spectrum
) {
    const int i = blockIdx.x * RND_BLOCK_SIZE + threadIdx.x;
    // curand_init(1234, j, 0, States + i);
    const int dim2 = (xdim - 1) * 2;  // Assuming square images (particles)
    const int size = xdim * dim2 * dim2;
    const int xydim = xdim * ydim;
    const int passes = size / (RND_BLOCK_NUM * RND_BLOCK_SIZE) + 1;
    for (int pass = 0, j = i; pass != passes && j != size; pass++, j += RND_BLOCK_NUM * RND_BLOCK_SIZE) {
        const int x = j % xdim;
              int y = j - (j / xdim);
              int z = j / xydim;
        // fftshift
        if (z >= xdim) z -= dim2;
        if (y >= xdim) y -= dim2;
        const int ires = rintf(device::hypot(x, y, z));
        const XFLOAT scale = ires < xdim ? spectrum[ires] : 0.0;
        Image[j] = CURAND_NORMAL2(States + i) * scale;
    }
}

#undef CURAND_NORMAL2

// __global__ void cuda_kernel_exponentiate_weights_fine2(
//     XFLOAT *g_pdf_orientation, XFLOAT *g_pdf_offset, XFLOAT *g_weights,
//     XFLOAT avg_diff2,
//     int oversamples_orient, int oversamples_trans,
//     unsigned long * d_rot_id, unsigned long * d_trans_idx,
//     unsigned long * d_job_idx, unsigned long * d_job_num,
//     long int job_num
// ) {
//     const long int i = threadIdx.x + blockIdx.x * SUMW_BLOCK_SIZE;
//     if (i >= job_num) return;
//     const long int j = d_job_idx[i];  // Index of comparison
//     const long int in = d_job_num[i];
//     for (int itrans = 0; itrans < in; itrans++) {
//         const XFLOAT a = g_weights[j + itrans] + avg_diff2;
//         #if defined(ACC_DOUBLE_PRECISION)
//         g_weights[j + itrans] = a < -700.0 ? 0.0 : exp(a);
//         #else
//         if (a < -88.)
//         g_weights[j + itrans] = a < -88.0 ? 0.0 : expf(a);
//         #endif
//     }
// }

__global__ void cuda_kernel_softMaskOutsideMap(
    image_view<XFLOAT, 3> vol, bool do_Mnoise,
    XFLOAT radius, XFLOAT radius_p, XFLOAT cosine_width
) {
    __shared__ XFLOAT     img_pixels[SOFTMASK_BLOCK_SIZE];
    __shared__ XFLOAT    partial_sum[SOFTMASK_BLOCK_SIZE];
    __shared__ XFLOAT partial_sum_bg[SOFTMASK_BLOCK_SIZE];

    long int texel_pass_num = ceilfracf(vol.size, SOFTMASK_BLOCK_SIZE);

    partial_sum   [threadIdx.x] = 0.0;
    partial_sum_bg[threadIdx.x] = 0.0;

    const long int xydim = vol.dim[0] * vol.dim[1];
    if (do_Mnoise) {
        for (int pass = 0, texel = threadIdx.x; pass < texel_pass_num; pass++, texel += SOFTMASK_BLOCK_SIZE) {
            // loop the available warps enough to complete all translations for this orientation
            if (texel < vol.size) {
                img_pixels[threadIdx.x] = __ldg(vol.data + texel);
                XFLOAT x, y, z;
                z = floordiv(texel, xydim);
                y = floordiv(texel - xydim * z, vol.dim[0]);
                x = texel - xydim * z - vol.dim[0] * y;
                z -= vol.init[2];
                y -= vol.init[1];
                x -= vol.init[0];
                const XFLOAT r = device::hypot(x, y, z);
                if (r < radius) {
                    continue;
                } else if (r > radius_p) {
                    partial_sum   [threadIdx.x] += 1.0;
                    partial_sum_bg[threadIdx.x] += img_pixels[threadIdx.x];
                } else {
                    const XFLOAT factor = RAISED_COSPI((radius_p - r) / cosine_width);
                    partial_sum   [threadIdx.x] += factor;
                    partial_sum_bg[threadIdx.x] += factor * img_pixels[threadIdx.x];
                }
            }
        }
    }

    // Parallel sum reduce
    __syncthreads();
    for (int j = SOFTMASK_BLOCK_SIZE / 2; j > 0; j /= 2) {
        if (threadIdx.x < j) {
            partial_sum   [threadIdx.x] += partial_sum   [threadIdx.x + j];
            partial_sum_bg[threadIdx.x] += partial_sum_bg[threadIdx.x + j];
        }
        __syncthreads();
    }

    const XFLOAT sum_bg_total = partial_sum_bg[0] / partial_sum[0];

    int texel = threadIdx.x;
    for (int pass = 0; pass < texel_pass_num; pass++, texel += SOFTMASK_BLOCK_SIZE) {
        // loop the available warps enough to complete all translations for this orientation
        if (texel < vol.size) {
            img_pixels[threadIdx.x] = __ldg(vol.data + texel);
            XFLOAT x, y, z;
            z = floordiv(texel, xydim);
            y = floordiv(texel - xydim * z, vol.dim[0]);
            x = texel - xydim * z - vol.dim[0] * y;
            z -= vol.init[2];
            y -= vol.init[1];
            x -= vol.init[0];
            const XFLOAT r = device::hypot(x, y, z);
            if (r < radius) {
                continue;
            } else if (r > radius_p) {
                img_pixels[threadIdx.x] = sum_bg_total;
            } else {
                const XFLOAT factor = RAISED_COSPI((radius_p - r) / cosine_width);
                img_pixels[threadIdx.x] = img_pixels[threadIdx.x] * (1 - factor) + sum_bg_total * factor;
            }
            vol.data[texel] = img_pixels[threadIdx.x];
        }
    }
}

__global__ void cuda_kernel_softMaskBackgroundValue(
    image_view<XFLOAT, 3> vol,
    XFLOAT radius, XFLOAT radius_p, XFLOAT cosine_width,
    XFLOAT *g_sum, XFLOAT *g_sum_bg
) {
    __shared__ XFLOAT     img_pixels[SOFTMASK_BLOCK_SIZE];
    __shared__ XFLOAT    partial_sum[SOFTMASK_BLOCK_SIZE];
    __shared__ XFLOAT partial_sum_bg[SOFTMASK_BLOCK_SIZE];

    const long int texel_pass_num = ceilfracf(vol.size, SOFTMASK_BLOCK_SIZE * gridDim.x);
    int texel = blockIdx.x * SOFTMASK_BLOCK_SIZE * texel_pass_num + threadIdx.x;

    partial_sum   [threadIdx.x] = 0.0;
    partial_sum_bg[threadIdx.x] = 0.0;

    const long int xydim = vol.dim[0] * vol.dim[1];
    // loop the available warps enough to complete all translations for this orientation
    for (int pass = 0; pass < texel_pass_num; pass++, texel += SOFTMASK_BLOCK_SIZE) {
        if (texel < vol.size) {
            img_pixels[threadIdx.x] = __ldg(vol.data + texel);
            const int z = texel / xydim              - vol.init[2];
            const int y = texel % xydim / vol.dim[0] - vol.init[1];
            const int x = texel % xydim % vol.dim[0] - vol.init[0];
            const XFLOAT r = device::hypot(x, y, z);
            if (r < radius) {
                continue;
            } else if (r > radius_p) {
                partial_sum   [threadIdx.x] += 1.0;
                partial_sum_bg[threadIdx.x] += img_pixels[threadIdx.x];
            } else {
                const XFLOAT factor = RAISED_COSPI((radius_p - r) / cosine_width);
                partial_sum   [threadIdx.x] += factor;
                partial_sum_bg[threadIdx.x] += factor * img_pixels[threadIdx.x];
            }
        }
    }

    cuda_atomic_add(g_sum    + threadIdx.x, partial_sum   [threadIdx.x]);
    cuda_atomic_add(g_sum_bg + threadIdx.x, partial_sum_bg[threadIdx.x]);
}


__global__ void cuda_kernel_cosineFilter(
    image_view<XFLOAT, 3> vol,
    bool do_noise, XFLOAT *noise,
    XFLOAT radius, XFLOAT radius_p, XFLOAT cosine_width, XFLOAT bg_value
) {
    __shared__ XFLOAT img_pixels[SOFTMASK_BLOCK_SIZE];

    const long int xydim = vol.dim[0] * vol.dim[1];
    // loop the available warps enough to complete all translations for this orientation
    const long int texel_pass_num = ceilfracf(vol.size, SOFTMASK_BLOCK_SIZE * gridDim.x);
    int texel = blockIdx.x * SOFTMASK_BLOCK_SIZE * texel_pass_num + threadIdx.x;
    for (int pass = 0; pass < texel_pass_num && texel < vol.size; pass++, texel += SOFTMASK_BLOCK_SIZE) {
        img_pixels[threadIdx.x] = __ldg(vol.data + texel);
        const int z = texel / xydim              - vol.init[2];
        const int y = texel % xydim / vol.dim[0] - vol.init[1];
        const int x = texel % xydim % vol.dim[0] - vol.init[0];
        const XFLOAT r = device::hypot(x, y, z);
        const XFLOAT def = do_noise ? noise[texel] : bg_value;
        if (r < radius) {
            continue;
        } else if (r > radius_p) {
            img_pixels[threadIdx.x] = def;
        } else {
            const XFLOAT factor = RAISED_COSPI((radius_p - r) / cosine_width);
            img_pixels[threadIdx.x] = img_pixels[threadIdx.x] * (1 - factor) + def * factor;

        }
        vol.data[texel] = img_pixels[threadIdx.x];
    }
}

__global__ void cuda_kernel_probRatio(
    XFLOAT * d_Mccf, XFLOAT * d_Mpsi, XFLOAT * d_Maux, XFLOAT * d_Mmean, XFLOAT * d_Mstddev,
    int image_size,
    XFLOAT normfft, XFLOAT sum_ref_under_circ_mask, XFLOAT sum_ref2_under_circ_mask, XFLOAT expected_Pratio,
    int NpsiThisBatch, int startPsi, int totalPsis
) {
    /* PLAN TO:
     *
     * 1) Pre-filter
     * 		d_Mstddev[i] = 1 / (2* d_Mstddev[i])   ( if d_Mstddev[pixel] > 1E-10 )
     * 		d_Mstddev[i] = 1    				  ( else )
     *
     * 2) Set
     * 		sum_ref2_under_circ_mask /= 2.
     *
     * 3) Total expression becomes
     * 		diff2 = ( exp(k) - 1.f ) / (expected_Pratio - 1.f)
     * 	  where
     * 	  	k = (normfft * d_Maux[pixel] + d_Mmean[pixel] * sum_ref_under_circ_mask)* d_Mstddev[i] + sum_ref2_under_circ_mask
     *
     */
    const int pixel = threadIdx.x + blockIdx.x * (int) PROBRATIO_BLOCK_SIZE;
    if (pixel < image_size) {
        XFLOAT Kccf = d_Mccf[pixel];
        XFLOAT Kpsi = -1.0;
        for (int psi = 0; psi < NpsiThisBatch; psi++) {
            XFLOAT diff2
                = normfft * d_Maux[pixel + image_size * psi]
                + d_Mmean[pixel] * sum_ref_under_circ_mask;
            // if (d_Mstddev[pixel] > 1E-10)
            diff2 = diff2 * d_Mstddev[pixel] + sum_ref2_under_circ_mask;
            // exponentiate to reflect the Gaussian error model. sigma=1 after normalization, 0.4=1/sqrt(2pi)
            #if defined(ACC_DOUBLE_PRECISION)
            diff2 = exp(-diff2 / 2.0);
            #else
            diff2 = expf(-diff2 / 2.f);
            #endif
            // Store fraction of (1 - probability-ratio) wrt (1 - expected Pratio)
            diff2 = (diff2 - 1.f) / (expected_Pratio - 1.f);
            if (diff2 > Kccf) {
                Kccf = diff2;
                Kpsi = (startPsi + psi) * (360 / totalPsis);
            }
        }
        d_Mccf[pixel] = Kccf;
        if (Kpsi >= 0.0)
            d_Mpsi[pixel] = Kpsi;
    }
}

__global__ void cuda_kernel_rotateOnly(
    acc::Complex * d_Faux,
    XFLOAT psi,
    AccProjectorKernel projector,
    int startPsi
) {
    const int proj = blockIdx.y;
    const int image_size = projector.imgX * projector.imgY;
    const int i = threadIdx.x + blockIdx.x * blockDim.x;
    if (i >= image_size) return;
    int y = floorfracf(i, projector.imgX);
    int x = i % projector.imgX;
    if (y > projector.maxR) {
        if (y >= projector.imgY - projector.maxR) {
            y -= projector.imgY;
        } else {
            x = projector.maxR;
        }
    }

    XFLOAT sa, ca;
    sincos((proj + startPsi) * psi, &sa, &ca);
    acc::Complex z;
    projector.project2Dmodel(x, y, ca, -sa, sa, ca, z.x, z.y);
    const long int j = proj * image_size + i;
    d_Faux[j] = z;
}

__global__ void cuda_kernel_rotateAndCtf(
    acc::Complex * d_Faux,
    XFLOAT const * d_ctf,
    XFLOAT psi,
    AccProjectorKernel projector,
    int startPsi
){
    const int proj = blockIdx.y;
    const int image_size = projector.imgX * projector.imgY;
    const int i = threadIdx.x + blockIdx.x * blockDim.x;
    if (i >= image_size) return;
    int y = floorfracf(i, projector.imgX);
    int x = i % projector.imgX;
    if (y > projector.maxR) {
        if (y >= projector.imgY - projector.maxR) {
            y -= projector.imgY;
        } else {
            x = projector.maxR;
        }
    }

    XFLOAT sa, ca;
    sincos((proj + startPsi) * psi, &sa, &ca);
    acc::Complex z;
    projector.project2Dmodel(x, y, ca, -sa, sa, ca, z.x, z.y);
    const long int j = proj * image_size + i;
    d_Faux[j].x = z.x * d_ctf[i];
    d_Faux[j].y = z.y * d_ctf[i];
}

__global__ void cuda_kernel_convol_A(
    acc::Complex * d_A, acc::Complex const * d_B, int image_size
) {
    const int i = threadIdx.x + blockIdx.x * blockDim.x;
    if (i >= image_size) return;
    // conj(d_A) * d_B
    const XFLOAT real = d_A[i].x * d_B[i].x + d_A[i].y * d_B[i].y;
    const XFLOAT imag = d_A[i].x * d_B[i].y - d_A[i].y * d_B[i].x;
    d_A[i].x = real;
    d_A[i].y = imag;
}

__global__ void cuda_kernel_convol_A(
    acc::Complex * d_A, acc::Complex const * d_B, acc::Complex * d_C, int image_size
) {
    const int i = threadIdx.x + blockIdx.x * blockDim.x;
    if (i >= image_size) return;
    // conj(d_A) * d_B
    d_C[i].x = d_A[i].x * d_B[i].x + d_A[i].y * d_B[i].y;
    d_C[i].y = d_A[i].x * d_B[i].y - d_A[i].y * d_B[i].x;
}

__global__ void cuda_kernel_batch_convol_A(
    acc::Complex * d_A, acc::Complex const * d_B, int image_size
) {
    const int i = threadIdx.x + blockIdx.x * blockDim.x;
    if (i >= image_size) return;
    const int j = blockIdx.y * image_size;
    const XFLOAT imag = d_A[i + j].x * d_B[i].x + d_A[i + j].y * d_B[i].y;
    const XFLOAT real = d_A[i + j].x * d_B[i].y - d_A[i + j].y * d_B[i].x;
    d_A[i + j].x = real;
    d_A[i + j].y = imag;
}

__global__ void cuda_kernel_batch_convol_A(
    acc::Complex * d_A, acc::Complex const * d_B, acc::Complex * d_C, int image_size
) {
    const int i = threadIdx.x + blockIdx.x * blockDim.x;
    if (i >= image_size) return;
    const int j = blockIdx.y * image_size;
    d_C[i + j].x = d_A[i + j].x * d_B[i].x + d_A[i + j].y * d_B[i].y;
    d_C[i + j].y = d_A[i + j].x * d_B[i].y - d_A[i + j].y * d_B[i].x;
}

__global__ void cuda_kernel_convol_B(
    acc::Complex * d_A, acc::Complex const * d_B, int image_size
) {
    const int i = threadIdx.x + blockIdx.x * blockDim.x;
    if (i >= image_size) return;
    const XFLOAT real = d_A[i].x * d_B[i].x + d_A[i].y * d_B[i].y;
    const XFLOAT imag = d_A[i].y * d_B[i].x - d_A[i].x * d_B[i].y;
    d_A[i].x = real;
    d_A[i].y = imag;
}

__global__ void cuda_kernel_convol_B(
    acc::Complex * d_A, acc::Complex const * d_B, acc::Complex * d_C, int image_size
) {
    const int i = threadIdx.x + blockIdx.x * blockDim.x;
    if (i >= image_size) return;
    d_C[i].x = d_A[i].x * d_B[i].x + d_A[i].y * d_B[i].y;
    d_C[i].y = d_A[i].y * d_B[i].x - d_A[i].x * d_B[i].y;
}

__global__ void cuda_kernel_batch_convol_B(
    acc::Complex * d_A, acc::Complex const * d_B, int image_size
) {
    const long int i = threadIdx.x + blockIdx.x * blockDim.x;
    const int j = blockIdx.y * image_size;
    if (i >= image_size) return;
    const XFLOAT real = d_A[i + j].x * d_B[i].x + d_A[i + j].y * d_B[i].y;
    const XFLOAT imag = d_A[i + j].y * d_B[i].x - d_A[i + j].x * d_B[i].y;
    d_A[i + j].x = real;
    d_A[i + j].y = imag;
}

__global__ void cuda_kernel_batch_multi(
    XFLOAT *A, XFLOAT *B, XFLOAT *OUT, XFLOAT S, int image_size
) {
    const int i = threadIdx.x + blockIdx.x * blockDim.x;
    if (i >= image_size) return;
    const int j = blockIdx.y * image_size;
    OUT[i + j] = A[i + j] * S * B[i + j];
}

__global__ void cuda_kernel_finalizeMstddev(
    XFLOAT *Mstddev, XFLOAT *aux, XFLOAT S, int image_size
) {
    const int i = threadIdx.x + blockIdx.x * blockDim.x;
    if (i >= image_size) return;
    const XFLOAT x = Mstddev[i] + S * aux[i];
    Mstddev[i] = x > 0 ? sqrt(x) : 0;
}

__global__ void cuda_kernel_square(XFLOAT *A, int image_size) {
    const int i = threadIdx.x + blockIdx.x * blockDim.x;
    if (i >= image_size) return;
    A[i] = A[i] * A[i];
}

template<bool invert>
__global__ void cuda_kernel_make_eulers_2D(
    XFLOAT *alphas, XFLOAT *eulers, unsigned orientation_num
) {
    // Orientation id
    const unsigned i = threadIdx.x + blockIdx.x * blockDim.x;
    if (i >= orientation_num) return;
    XFLOAT ca, sa;
    const XFLOAT a = device::radians(invert ? -alphas[i] : alphas[i]);
    #ifdef ACC_DOUBLE_PRECISION
    sincos(a, &sa, &ca);
    #else
    sincosf(a, &sa, &ca);
    #endif
    eulers[9 * i + 0] = +ca;
    eulers[9 * i + 1] = +sa;
    eulers[9 * i + 2] = 0;
    eulers[9 * i + 3] = -sa;
    eulers[9 * i + 4] = +ca;
    eulers[9 * i + 5] = 0;
    eulers[9 * i + 6] = 0;
    eulers[9 * i + 7] = 0;
    eulers[9 * i + 8] = 1;
}

template<bool invert, bool doL, bool doR>
__global__ void cuda_kernel_make_eulers_3D(
        XFLOAT *alphas, XFLOAT *betas, XFLOAT *gammas, XFLOAT *eulers,
        unsigned orientation_num,
        XFLOAT *L, XFLOAT *R
) {
    const unsigned i = blockIdx.x * blockDim.x + threadIdx.x;  // Orientation id
    if (i >= orientation_num) return;

    const XFLOAT a = device::radians(alphas[i]);
    const XFLOAT b = device::radians(betas [i]);
    const XFLOAT g = device::radians(gammas[i]);

    XFLOAT ca, sa, cb, sb, cg, sg;
    #ifdef ACC_DOUBLE_PRECISION
    sincos(a, &sa, &ca);
    sincos(b, &sb, &cb);
    sincos(g, &sg, &cg);
    #else
    sincosf(a, &sa, &ca);
    sincosf(b, &sb, &cb);
    sincosf(g, &sg, &cg);
    #endif

    const XFLOAT cc = cb * ca;
    const XFLOAT cs = cb * sa;
    const XFLOAT sc = sb * ca;
    const XFLOAT ss = sb * sa;

    XFLOAT A[9], B[9];
    A[0] = cg * cc - sg * sa;   // 00
    A[1] = cg * cs + sg * ca;   // 01
    A[2] = -cg * sb;            // 02
    A[3] = -sg * cc - cg * sa;  // 10
    A[4] = -sg * cs + cg * ca;  // 11
    A[5] = sg * sb;             // 12
    A[6] = sc;                  // 20
    A[7] = ss;                  // 21
    A[8] = cb;                  // 22    

    if (doR) {
        for (int i = 0; i < 9; i++)
            B[i] = 0.0;
        for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
        for (int k = 0; k < 3; k++)
            B[i * 3 + j] += A[i * 3 + k] * R[k * 3 + j];
    } else {
        for (int i = 0; i < 9; i++)
            B[i] = A[i];
    }

    if (doL) {
        if (doR) {
            for (int i = 0; i < 9; i++)
                A[i] = B[i];
        }

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
            eulers[9 * i + 0] = B[0];  // 00
            eulers[9 * i + 1] = B[3];  // 01
            eulers[9 * i + 2] = B[6];  // 02
            eulers[9 * i + 3] = B[1];  // 10
            eulers[9 * i + 4] = B[4];  // 11
            eulers[9 * i + 5] = B[7];  // 12
            eulers[9 * i + 6] = B[2];  // 20
            eulers[9 * i + 7] = B[5];  // 21
            eulers[9 * i + 8] = B[8];  // 22
        }
    } else {
        for (int j = 0; j != 9; j++)
            eulers[9 * i + j] = B[j];
    }
}

__global__ void cuda_kernel_allweights_to_mweights(
    unsigned long * d_iorient, XFLOAT * d_allweights, XFLOAT * d_mweights,
    unsigned long orientation_num, unsigned long translation_num, int block_size
) {
    const size_t i = blockIdx.x * block_size + threadIdx.x;
    if (i >= orientation_num * translation_num) return;
    const long quot = i / translation_num;
    const long rem  = i % translation_num;
    d_mweights[d_iorient[quot] * translation_num + rem] = d_allweights[i];
}

__global__ void cuda_kernel_initOrientations(
    RFLOAT *pdfs, XFLOAT *pdf_orientation, bool *pdf_orientation_zeros, size_t sz
) {
    const int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= sz) return;
    pdf_orientation_zeros[i] = pdfs[i] == 0.0;
    pdf_orientation[i]       = pdfs[i] == 0.0 ? 0.0 : log(pdfs[i]);
}

__global__ void cuda_kernel_griddingCorrect(
    RFLOAT *vol, int interpolator, RFLOAT rrval, RFLOAT r_min_nn,
    size_t iX, size_t iY, size_t iZ
) {
    const int idx = blockIdx.x * blockDim.x + threadIdx.x;
    const int idy = blockIdx.y * blockDim.y + threadIdx.y;
    const int idz = blockIdx.z * blockDim.z + threadIdx.z;
    if (idx >= iX || idy >= iY || idz >= iZ) return;
    const RFLOAT r = device::hypot(idx - iX / 2, idy - iY / 2, idz - iZ / 2);
    if (r <= 0.0) return;
    const RFLOAT theta = r / rrval;
    const RFLOAT sinc = sin(PI * theta) / (PI * theta);
    if (interpolator == NEAREST_NEIGHBOUR && r_min_nn == 0.0)
        vol[idz * iX * iY + idy * iX + idx] /= sinc;
    else if (interpolator == TRILINEAR || interpolator == NEAREST_NEIGHBOUR && r_min_nn > 0)
        vol[idz * iX * iY + idy * iX + idx] /= sinc * sinc;
}

__global__ void cuda_kernel_updatePowerSpectrum(
    RFLOAT *dcounter, RFLOAT *dpower_spectrum, int sz
) {
    const int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= sz) return;
    if (dcounter[i] < 1.0)
        dpower_spectrum[i] = 0.0;
    else
        dpower_spectrum[i] /= dcounter[i];
}
