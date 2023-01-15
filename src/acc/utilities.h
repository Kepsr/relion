#ifndef ACC_UTILITIES_H_
#define ACC_UTILITIES_H_

#include "src/acc/cpu/cpu_kernels/wavg.h"
#include "src/acc/cpu/cpu_kernels/helper.h"
#include "src/acc/cpu/cpu_helper_functions.h"

#include "src/acc/acc_ptr.h"
#include "src/acc/data_types.h"
#include "src/error.h"
#ifdef CUDA
#include "src/acc/cuda/cuda_kernels/helper.cuh"
#include "src/acc/cuda/cuda_kernels/wavg.cuh"
#include "src/acc/cuda/cuda_kernels/diff2.cuh"
#include "src/acc/cuda/cuda_fft.h"
#else
#include "src/acc/cpu/cpu_kernels/helper.h"
#include "src/acc/cpu/cpu_kernels/wavg.h"
#include "src/acc/cpu/cpu_kernels/diff2.h"
#endif

void dump_array(char *name, bool *ptr, size_t size);
void dump_array(char *name, int *ptr, size_t size);
void dump_array(char *name, size_t *ptr, size_t size);
void dump_array(char *name, float *ptr, size_t size);
void dump_complex_array(char *name, acc::Complex *ptr, size_t size);
void dump_complex_array(char *name, Complex *ptr, size_t size);
void dump_double_array(char *name, float *ptr, float *ptr2, size_t size);
void dump_triple_array(char *name, float *ptr, float *ptr2, float *ptr3, size_t size);
void dump_array(char *name, double *ptr, size_t size);
void dump_double_array(char *name, double *ptr, double *ptr2, size_t size);
void dump_triple_array(char *name, double *ptr, double *ptr2, double *ptr3, size_t size);

namespace AccUtilities {

template <typename T>
static T getSumOnDevice(AccPtr<T, acc::cpu> &ptr) {
    #ifdef DEBUG_CUDA
    if (ptr.getSize() == 0)
        printf("DEBUG_ERROR: getSumOnDevice called with pointer of zero size.\n");
    if (!ptr.getHostPtr())
        printf("DEBUG_ERROR: getSumOnDevice called with null device pointer.\n");
    #endif
    return std::accumulate(ptr.getHostPtr(), ptr.getHostPtr() + ptr.getSize(), T(0));
}

template <typename T>
static T getMinOnDevice(AccPtr<T, acc::cpu> &ptr) {
    #ifdef DEBUG_CUDA
    if (ptr.getSize() == 0)
        printf("DEBUG_ERROR: getMinOnDevice called with pointer of zero size.\n");
    if (!ptr.getHostPtr())
        printf("DEBUG_ERROR: getMinOnDevice called with null device pointer.\n");
    #endif
    return CpuKernels::getMin<T>(ptr.getAccPtr(), ptr.getSize());
}

template <typename T>
static T getMaxOnDevice(AccPtr<T, acc::cpu> &ptr) {
    #ifdef DEBUG_CUDA
    if (ptr.getSize() == 0)
        printf("DEBUG_ERROR: getMaxOnDevice called with pointer of zero size.\n");
    if (!ptr.getHostPtr())
        printf("DEBUG_ERROR: getMaxOnDevice called with null device pointer.\n");
    #endif
    return CpuKernels::getMax<T>(ptr.getAccPtr(), ptr.getSize());
}

template <typename T>
static std::pair<size_t, T> getArgMinOnDevice(AccPtr<T, acc::cpu> &ptr) {
    #ifdef DEBUG_CUDA
    if (ptr.getSize() == 0)
        printf("DEBUG_ERROR: getArgMinOnDevice called with pointer of zero size.\n");
    if (!ptr.getHostPtr())
        printf("DEBUG_ERROR: getArgMinOnDevice called with null device pointer.\n");
    #endif
    return CpuKernels::getArgMin<T>(ptr.getAccPtr(), ptr.getSize());
}

template <typename T>
static std::pair<size_t, T> getArgMaxOnDevice(AccPtr<T, acc::cpu> &ptr) {
    #ifdef DEBUG_CUDA
    if (ptr.getSize() == 0)
        printf("DEBUG_ERROR: getArgMaxOnDevice called with pointer of zero size.\n");
    if (!ptr.getHostPtr())
        printf("DEBUG_ERROR: getArgMaxOnDevice called with null device pointer.\n");
    #endif
    return CpuKernels::getArgMax<T>(ptr.getAccPtr(), ptr.getSize());
}

template <typename T>
static void filterGreaterZeroOnDevice(AccPtr<T, acc::cpu> &in, AccPtr<T, acc::cpu> &out) {
    // Determine size of output array
    const size_t filt_size = std::count_if(
        in.getHostPtr(), in.getHostPtr() + in.getSize(),
        [] (const T& x) -> bool { return x > 0.0; });
    #ifdef DEBUG_CUDA
    if (filt_size == 0)
        in.HandleDebugFatal("filterGreaterZeroOnDevice - No filtered values greater than 0.\n", __FILE__, __LINE__);
    #endif
    out.resizeHost(filt_size);
    // Populate output array
    for (size_t i = 0, j = 0; i < in.getSize(); i++)
        if (in.getHostPtr()[i] > (T) 0.0)
            out.getHostPtr()[j++] = in.getHostPtr()[i];
}

template <typename T>
static void sortOnDevice(AccPtr<T, acc::cpu> &in, AccPtr<T, acc::cpu> &out) {
    T const* const src = in.getAccPtr();
    T* const dest = out.getHostPtr();
    const size_t n = in.getSize();
    std::copy_n(src, n, dest);
    std::sort(dest, dest + n);
}

template <typename T>
static void scanOnDevice(AccPtr<T, acc::cpu> &in, AccPtr<T, acc::cpu> &out) {
    T sum = 0.0;
    for (size_t i = 0; i < in.getSize(); i++) {
        sum += in.getHostPtr()[i];
        out.getHostPtr()[i] = sum;
    }
}

static void TranslateAndNormCorrect(
    MultidimArray<RFLOAT> &img_in, AccPtr<XFLOAT> &img_out,
    XFLOAT normcorr, RFLOAT xOff, RFLOAT yOff, RFLOAT zOff,
    bool DATA3D
);

static void softMaskBackgroundValue(
    int inblock_dim, int inblock_size,
    XFLOAT *vol, Image<RFLOAT> &img,
    XFLOAT radius, XFLOAT radius_p, XFLOAT cosine_width,
    XFLOAT *g_sum, XFLOAT *g_sum_bg
);

static void cosineFilter(
    int inblock_dim, int inblock_size,
    XFLOAT *vol, long int vol_size,
    long int xdim, long int ydim, long int zdim,
    long int xinit, long int yinit, long int zinit,
    bool do_Mnoise,
    XFLOAT radius, XFLOAT radius_p, XFLOAT cosine_width, XFLOAT sum_bg_total
);

template <acc::Type>
void powerClass(
    int grid_dimensions, int block_dimensions,
    acc::Complex *g_image, XFLOAT *g_spectrum,
    size_t image_size, size_t spectrum_size,
    int xdim, int ydim, int zdim,
    int res_limit, XFLOAT *g_highres_Xi2, bool DATA3D
);

template <>
inline void powerClass<acc::cpu>(
    int grid_dimensions, int block_dimensions,
    acc::Complex *g_image, XFLOAT *g_spectrum,
    size_t image_size, size_t spectrum_size,
    int xdim, int ydim, int zdim,
    int res_limit, XFLOAT *g_highres_Xi2, bool DATA3D
) {
    if (DATA3D)
        CpuKernels::powerClass<true>(
            grid_dimensions, g_image, g_spectrum,
            image_size, spectrum_size,
            xdim, ydim, zdim,
            res_limit, g_highres_Xi2
        );
    else
        CpuKernels::powerClass<false>(
            grid_dimensions, g_image, g_spectrum,
            image_size, spectrum_size,
            xdim, ydim, zdim,
            res_limit, g_highres_Xi2
        );
}

template <acc::Type accType>
void acc_make_eulers_2D(
    int grid_size, int block_size, cudaStream_t stream,
    XFLOAT *alphas, XFLOAT *eulers,
    unsigned long orientation_num, bool invert
);

template <>
inline void acc_make_eulers_2D<acc::cpu>(
    int grid_size, int block_size, cudaStream_t stream,
    XFLOAT *alphas, XFLOAT *eulers,
    unsigned long orientation_num, bool invert
) {
    if (invert)
        CpuKernels::cpu_kernel_make_eulers_2D<true>
            (grid_size, block_size, alphas, eulers, orientation_num);
    else
        CpuKernels::cpu_kernel_make_eulers_2D<false>
            (grid_size, block_size, alphas, eulers, orientation_num);
}

template <acc::Type>
void acc_make_eulers_3D(
    int grid_size, int block_size, cudaStream_t stream,
    XFLOAT *alphas, XFLOAT *betas, XFLOAT *gammas, XFLOAT *eulers,
    unsigned long orientation_num,
    XFLOAT *L, XFLOAT *R, bool doL, bool doR, bool invert
);

template <>
inline void acc_make_eulers_3D<acc::cpu>(
    int grid_size, int block_size, cudaStream_t stream,
    XFLOAT *alphas, XFLOAT *betas, XFLOAT *gammas, XFLOAT *eulers,
    unsigned long orientation_num,
    XFLOAT *L, XFLOAT *R, bool doL, bool doR, bool invert
) {
    #define MAKE_EULERS( A, B, C ) \
        CpuKernels::cpu_kernel_make_eulers_3D<A, B, C> \
        (grid_size, block_size, alphas, betas, gammas, eulers, orientation_num, L, R);
    const int bits = (int) doL << 2 | (int) doR << 1 | (int) invert;
    switch (bits) {
        case 0: MAKE_EULERS( 0, 0, 0 ); break;
        case 1: MAKE_EULERS( 0, 0, 1 ); break;
        case 2: MAKE_EULERS( 0, 1, 0 ); break;
        case 3: MAKE_EULERS( 0, 1, 1 ); break;
        case 4: MAKE_EULERS( 1, 0, 0 ); break;
        case 5: MAKE_EULERS( 1, 0, 1 ); break;
        case 6: MAKE_EULERS( 1, 1, 0 ); break;
        case 7: MAKE_EULERS( 1, 1, 1 ); break;
    }
    #undef MAKE_EULERS
}

template <acc::Type>
void frequencyPass(
    int grid_size, int block_size, cudaStream_t stream,
    acc::Complex *A, long int ori_size,
    size_t Xdim, size_t Ydim, size_t Zdim,
    XFLOAT edge_low, XFLOAT edge_width, XFLOAT edge_high,
    XFLOAT angpix, size_t image_size, bool do_highpass
);

template <acc::Type>
void kernel_wavg(
    int block_sz, XFLOAT *g_eulers,
    AccProjectorKernel &projector,
    unsigned long image_size, unsigned long orientation_num,
    XFLOAT *g_img_real, XFLOAT *g_img_imag,
    XFLOAT *g_trans_x, XFLOAT *g_trans_y, XFLOAT *g_trans_z,
    XFLOAT* g_weights, XFLOAT* g_ctfs,
    XFLOAT *g_wdiff2s_parts, XFLOAT *g_wdiff2s_AA, XFLOAT *g_wdiff2s_XA,
    unsigned long translation_num,
    XFLOAT weight_norm, XFLOAT significant_weight, XFLOAT part_scale,
    cudaStream_t stream,
    bool CTFPREMULTIPLIED, bool REFCTF, bool REF3D, bool DATA3D
);

#ifdef CUDA
#define INIT_VALUE_BLOCK_SIZE 512

template <typename T>
static T getSumOnDevice(AccPtr<T, acc::cuda> &ptr) {
    return CudaKernels::getSumOnDevice<T>(ptr);
}

template <typename T>
static T getMinOnDevice(AccPtr<T, acc::cuda> &ptr) {
    return CudaKernels::getMinOnDevice<T>(ptr);
}

template <typename T>
static T getMaxOnDevice(AccPtr<T, acc::cuda> &ptr) {
    return CudaKernels::getMaxOnDevice<T>(ptr);
}

template <typename T>
static std::pair<size_t, T> getArgMinOnDevice(AccPtr<T, acc::cuda> &ptr) {
    return CudaKernels::getArgMinOnDevice<T>(ptr);
}

template <typename T>
static std::pair<size_t, T> getArgMaxOnDevice(AccPtr<T, acc::cuda> &ptr) {
    return CudaKernels::getArgMaxOnDevice<T>(ptr);
}

template <typename T>
static void filterGreaterZeroOnDevice(AccPtr<T, acc::cuda> &in, AccPtr<T, acc::cuda> &out) {
    CudaKernels::device_greater_than<T> closure (0);
    CudaKernels::filterOnDevice(in, out, closure);
}

template <typename T>
static void sortOnDevice(AccPtr<T, acc::cuda> &in, AccPtr<T, acc::cuda> &out) {
    CudaKernels::sortOnDevice(in, out);
}

template <typename T>
static void scanOnDevice(AccPtr<T, acc::cuda> &in, AccPtr<T, acc::cuda> &out) {
    CudaKernels::scanOnDevice(in, out);
}

template <>
inline void powerClass<acc::cuda>(
    int grid_dimensions, int block_dimensions,
    acc::Complex *g_image, XFLOAT *g_spectrum,
    size_t image_size, size_t spectrum_size,
    int xdim, int ydim, int zdim,
    int res_limit, XFLOAT *g_highres_Xi2, bool DATA3D
) {
    if (DATA3D)
        cuda_kernel_powerClass<1>
        <<<dim3(grid_dimensions), block_dimensions, 0, 0>>>(
            g_image, g_spectrum, image_size, spectrum_size,
            xdim, ydim, zdim,
            res_limit, g_highres_Xi2
        );
    else
        cuda_kernel_powerClass<0>
        <<<dim3(grid_dimensions), block_dimensions, 0, 0>>>(
            g_image, g_spectrum, image_size, spectrum_size,
            xdim, ydim, zdim,
            res_limit, g_highres_Xi2
        );
}

template <>
inline void acc_make_eulers_2D<acc::cuda>(
    int grid_size, int block_size, cudaStream_t stream,
    XFLOAT *alphas, XFLOAT *eulers,
    unsigned long orientation_num, bool invert
) {
    if (invert)
        cuda_kernel_make_eulers_2D<1><<<grid_size, block_size, 0, stream>>>
            (alphas, eulers, orientation_num);
    else 
        cuda_kernel_make_eulers_2D<0><<<grid_size, block_size, 0, stream>>>
            (alphas, eulers, orientation_num);
}

template <>
inline void acc_make_eulers_3D<acc::cuda>(
    int grid_size, int block_size, cudaStream_t stream,
    XFLOAT *alphas, XFLOAT *betas, XFLOAT *gammas, XFLOAT *eulers,
    unsigned long orientation_num,
    XFLOAT *L, XFLOAT *R, bool doL, bool doR, bool invert
) {
    #define MAKE_EULERS( A, B, C ) \
        cuda_kernel_make_eulers_3D<A, B, C> \
        <<<grid_size, block_size, 0, stream>>> \
        (alphas, betas, gammas, eulers, orientation_num, L, R)
    const int bits = (int) doL << 2 | (int) doR << 1 | (int) invert;
    switch (bits) {
        case 0: MAKE_EULERS( 0, 0, 0 ); break;
        case 1: MAKE_EULERS( 0, 0, 1 ); break;
        case 2: MAKE_EULERS( 0, 1, 0 ); break;
        case 3: MAKE_EULERS( 0, 1, 1 ); break;
        case 4: MAKE_EULERS( 1, 0, 0 ); break;
        case 5: MAKE_EULERS( 1, 0, 1 ); break;
        case 6: MAKE_EULERS( 1, 1, 0 ); break;
        case 7: MAKE_EULERS( 1, 1, 1 ); break;
    }
    #undef MAKE_EULERS
}

template <typename T>
void InitComplexValue(AccPtr<T, acc::cuda> &data, XFLOAT value) {
    const int grid_size = ceil((float) data.getSize() / (float) INIT_VALUE_BLOCK_SIZE);
    cuda_kernel_init_complex_value<T><<<grid_size, INIT_VALUE_BLOCK_SIZE, 0, data.getStream()>>>(
        data.getAccPtr(), value, data.getSize(), INIT_VALUE_BLOCK_SIZE
    );
}

template <typename T>
void InitValue(AccPtr<T, acc::cuda> &data, T value, size_t size) {
    const int grid_size = ceil((float) size / (float) INIT_VALUE_BLOCK_SIZE);
    cuda_kernel_init_value<T><<<grid_size, INIT_VALUE_BLOCK_SIZE, 0, data.getStream()>>>(
        data.getAccPtr(), value, size, INIT_VALUE_BLOCK_SIZE);
    LAUNCH_HANDLE_ERROR(cudaGetLastError());
}

template <>
inline void frequencyPass<acc::cuda>(
    int grid_size, int block_size,
    cudaStream_t stream,
    acc::Complex *A, long int ori_size,
    size_t Xdim, size_t Ydim, size_t Zdim,
    XFLOAT edge_low, XFLOAT edge_width, XFLOAT edge_high,
    XFLOAT angpix, size_t image_size, bool do_highpass
) {
    if (do_highpass) 
        cuda_kernel_frequencyPass
        <<<dim3(grid_size), block_size, 0, stream>>>(
            A, ori_size, Xdim, Ydim, Zdim,
            angpix, image_size, device_high_pass(edge_low, edge_high, edge_width)
        );
    else
        cuda_kernel_frequencyPass
        <<<dim3(grid_size), block_size, 0, stream>>>(
            A, ori_size, Xdim, Ydim, Zdim,
            angpix, image_size, device_low_pass(edge_low, edge_high, edge_width)
        );
}

template <>
inline void kernel_wavg<acc::cuda>(
    int block_sz, XFLOAT *g_eulers,
    AccProjectorKernel &projector,
    unsigned long image_size, unsigned long orientation_num,
    XFLOAT *g_img_real, XFLOAT *g_img_imag,
    XFLOAT *g_trans_x, XFLOAT *g_trans_y, XFLOAT *g_trans_z,
    XFLOAT* g_weights, XFLOAT* g_ctfs,
    XFLOAT *g_wdiff2s_parts, XFLOAT *g_wdiff2s_AA, XFLOAT *g_wdiff2s_XA,
    unsigned long translation_num,
    XFLOAT weight_norm, XFLOAT significant_weight, XFLOAT part_scale,
    cudaStream_t stream,
    bool CTFPREMULTIPLIED, bool REFCTF, bool REF3D, bool DATA3D
) {
    // We only want as many blocks as there are chunks of orientations to be treated
    // within the same block (this is done to reduce memory loads in the kernel).
    // ceil((float) orientation_num / (float) REF_GROUP_SIZE);
    #define KERNEL( A, B, C, D ) cuda_kernel_wavg<A, B, C, D> \
    <<<dim3(orientation_num), block_sz, (3 * block_sz + 9) * sizeof(XFLOAT), stream>>>( \
        block_sz, g_eulers, projector, image_size, orientation_num, \
        g_img_real, g_img_imag, \
        g_trans_x, g_trans_y, g_trans_z, \
        g_weights, g_ctfs, \
        g_wdiff2s_parts, g_wdiff2s_AA, g_wdiff2s_XA, \
        translation_num, weight_norm, significant_weight, part_scale \
    )
    const int bits = (int) CTFPREMULTIPLIED << 3 | (int) REFCTF << 2 | (int) REF3D << 1 | (int) DATA3D << 0;
    switch (bits) {
        case  0: KERNEL( 0, 0, 0, 0 ); break;
        case  1: KERNEL( 0, 0, 0, 1 ); break;
        case  2: KERNEL( 0, 0, 1, 0 ); break;
        case  3: KERNEL( 0, 0, 1, 1 ); break;
        case  4: KERNEL( 0, 1, 0, 0 ); break;
        case  5: KERNEL( 0, 1, 0, 1 ); break;
        case  6: KERNEL( 0, 1, 1, 0 ); break;
        case  7: KERNEL( 0, 1, 1, 1 ); break;
        case  8: KERNEL( 1, 0, 0, 0 ); break;
        case  9: KERNEL( 1, 0, 0, 1 ); break;
        case 10: KERNEL( 1, 0, 1, 0 ); break;
        case 11: KERNEL( 1, 0, 1, 1 ); break;
        case 12: KERNEL( 1, 1, 0, 0 ); break;
        case 13: KERNEL( 1, 1, 0, 1 ); break;
        case 14: KERNEL( 1, 1, 1, 0 ); break;
        case 15: KERNEL( 1, 1, 1, 1 ); break;
    }
    #undef KERNEL
}

template <typename T>
void kernel_weights_exponent_coarse(
    size_t num_classes,
    AccPtr<T, acc::cuda> &g_pdf_orientation, AccPtr<bool, acc::cuda> &g_pdf_orientation_zeros,
    AccPtr<T, acc::cuda> &g_pdf_offset, AccPtr<bool, acc::cuda> &g_pdf_offset_zeros,
    AccPtr<T, acc::cuda> &g_Mweight,
    T g_min_diff2,
    size_t nr_coarse_orient, size_t  nr_coarse_trans
) {
    const size_t block_num = ceilf(((double) nr_coarse_orient * nr_coarse_trans * num_classes) / (double) SUMW_BLOCK_SIZE);
    cuda_kernel_weights_exponent_coarse<T>
    <<<block_num, SUMW_BLOCK_SIZE, 0, g_Mweight.getStream()>>>(
        g_pdf_orientation.getAccPtr(), g_pdf_orientation_zeros.getAccPtr(),
        g_pdf_offset.getAccPtr(), g_pdf_offset_zeros.getAccPtr(),
        g_Mweight.getAccPtr(),
        g_min_diff2,
        nr_coarse_orient, nr_coarse_trans,
        nr_coarse_orient * nr_coarse_trans * num_classes
    );
}

template <typename T>
void kernel_exponentiate(AccPtr<T, acc::cuda> &array, T add) {
    const int Dg = ceilf((float) array.getSize() / BLOCK_SIZE);
    cuda_kernel_exponentiate<T>
    <<<Dg, BLOCK_SIZE, 0, array.getStream()>>>
        (array.getAccPtr(), add, array.getSize());
}
#endif

template <typename T>
void InitComplexValue(AccPtr<T, acc::cpu> &data, XFLOAT value) {
    for (size_t i = 0; i < data.getSize(); i++) {
        data.getHostPtr()[i].x = value;
        data.getHostPtr()[i].y = value;
    }
}

template <typename T>
void InitValue(AccPtr<T, acc::cpu> &data, T value, size_t size) {
    std::fill_n(data.getHostPtr(), size, value);
}

template <typename T, acc::Type accType>
void InitValue(AccPtr<T, accType> &data, T value) {
    InitValue(data, value, data.getSize());
}

void initOrientations(AccPtr<RFLOAT> &pdfs, AccPtr<XFLOAT> &pdf_orientation, AccPtr<XFLOAT> &pdf_orientation_zeros);

namespace CpuKernels {

void centerFFT_2D(
    int grid_size, int batch_size, int block_size, cudaStream_t stream,
    XFLOAT *img_in, size_t image_size,
    int xdim, int ydim, int xshift, int yshift
);

void centerFFT_2D(
    int grid_size, int batch_size, int block_size,
    XFLOAT *img_in, size_t image_size,
    int xdim, int ydim, int xshift, int yshift
);

void centerFFT_3D(
    int grid_size, int batch_size, int block_size, cudaStream_t stream,
    XFLOAT *img_in, size_t image_size,
    int xdim, int ydim, int zdim, int xshift, int yshift, int zshift
);

#ifdef ALTCPU
template <typename T, acc::Type accType>
void translate(
    int block_size, AccDataTypes::Image<T, accType> &in, AccDataTypes::Image<T, accType> &out,
    int dx, int dy, int dz = 0
) {
    if (in.getAccPtr() == out.getAccPtr())
        CRITICAL(ERRUNSAFEOBJECTREUSE);
    if (in.is3D()) {
        ::CpuKernels::cpu_translate3D<T>(
            in(), out(), in.getxyz(), in.getx(), in.gety(), in.getz(), dx, dy, dz
        );
    } else {
        ::CpuKernels::cpu_translate2D<T>(
            in(), out(), in.getxyz(), in.getx(), in.gety(), dx, dy
        );
    }
}

template <typename T, acc::Type accType>
void multiply(
    int block_dimensions, AccDataTypes::Image<T, accType> &ptr, T value
) {
    ::CpuKernels::cpu_kernel_multi<T>(ptr.getAccPtr(), value, ptr.getSize());
}

template <typename T>
void multiply(
    int grid_dimensions, int block_dimensions, cudaStream_t stream, T *array, T value, size_t size
) {
    ::CpuKernels::cpu_kernel_multi<T>(array, value, size);
}

template <bool REF3D, bool DATA3D, int block_sz, int eulers_per_block, int prefetch_fraction>
void diff2_coarse(
    unsigned long grid_size, int block_size, XFLOAT *g_eulers,
    XFLOAT *trans_x, XFLOAT *trans_y, XFLOAT *trans_z,
    XFLOAT *g_real, XFLOAT *g_imag,
    AccProjectorKernel projector,
    XFLOAT *g_corr, XFLOAT *g_diff2s,
    unsigned long translation_num, unsigned long image_size,
    cudaStream_t stream
) {
    #if 1
    ::CpuKernels::diff2_coarse<REF3D, DATA3D, block_sz, eulers_per_block, prefetch_fraction>(
        grid_size, g_eulers, trans_x, trans_y, trans_z,
        g_real, g_imag,
        projector, g_corr, g_diff2s,
        translation_num, image_size
    );
    #else
    if (DATA3D) {
        ::CpuKernels::diff2_coarse_3D<eulers_per_block>(
            grid_size, g_eulers, trans_x, trans_y, trans_z,
            g_real, g_imag,
            projector, g_corr, g_diff2s,
            translation_num, image_size
        );
    } else {
        ::CpuKernels::diff2_coarse_2D<REF3D, eulers_per_block>(
            grid_size, g_eulers,
            trans_x, trans_y, trans_z,
            g_real, g_imag,
            projector, g_corr, g_diff2s,
            translation_num, image_size
        );
    }
    #endif
}

template <bool REF3D, bool DATA3D, int block_sz>
void diff2_CC_coarse(
    unsigned long grid_size, int block_size,
    XFLOAT *g_eulers,
    XFLOAT *g_imgs_real, XFLOAT *g_imgs_imag,
    XFLOAT *g_trans_x, XFLOAT *g_trans_y, XFLOAT *g_trans_z,
    AccProjectorKernel projector,
    XFLOAT *g_corr_img, XFLOAT *g_diff2s,
    unsigned long translation_num, unsigned long image_size,
    XFLOAT exp_local_sqrtXi2,
    cudaStream_t stream
) {
    if (DATA3D) {
        ::CpuKernels::diff2_CC_coarse_3D(
            grid_size,
            g_eulers, g_imgs_real, g_imgs_imag,
            g_trans_x, g_trans_y, g_trans_z,
            projector,
            g_corr_img, g_diff2s,
            translation_num, image_size, exp_local_sqrtXi2
        );
    } else {
        ::CpuKernels::diff2_CC_coarse_2D<REF3D>(
            grid_size,
            g_eulers, g_imgs_real, g_imgs_imag,
            g_trans_x, g_trans_y,
            projector,
            g_corr_img, g_diff2s,
            translation_num, image_size, exp_local_sqrtXi2
        );
    }
}

template <bool REF3D, bool DATA3D, int block_sz, int chunk_sz>
void diff2_fine(
    unsigned long grid_size, int block_size,
    XFLOAT *g_eulers, XFLOAT *g_imgs_real, XFLOAT *g_imgs_imag,
    XFLOAT *trans_x, XFLOAT *trans_y, XFLOAT *trans_z,
    AccProjectorKernel projector,
    XFLOAT *g_corr_img, XFLOAT *g_diff2s,
    unsigned long image_size, XFLOAT sum_init,
    unsigned long orientation_num, unsigned long translation_num,
    unsigned long todo_blocks,
    unsigned long *d_rot_idx, unsigned long *d_trans_idx,
    unsigned long *d_job_idx, unsigned long *d_job_num,
    cudaStream_t stream
) {
    // In these non-CC kernels, g_corr_img is effectively an adjusted MinvSigma2
    if (DATA3D)
        ::CpuKernels::diff2_fine_3D(
            grid_size,
            g_eulers, g_imgs_real, g_imgs_imag,
            trans_x, trans_y, trans_z,
            projector,
            g_corr_img, g_diff2s,
            image_size, sum_init,
            orientation_num, translation_num,
            todo_blocks,  // significant_num,
            d_rot_idx, d_trans_idx, d_job_idx, d_job_num
        );
    else
        ::CpuKernels::diff2_fine_2D<REF3D>(
            grid_size,
            g_eulers, g_imgs_real, g_imgs_imag,
            trans_x, trans_y, trans_z,
            projector,
            g_corr_img, g_diff2s,
            image_size, sum_init,
            orientation_num, translation_num,
            todo_blocks,  // significant_num,
            d_rot_idx, d_trans_idx, d_job_idx, d_job_num
        );
}

template <bool REF3D, bool DATA3D, int block_sz,int chunk_sz>
void diff2_CC_fine(
    unsigned long grid_size, int block_size,
    XFLOAT *g_eulers, XFLOAT *g_imgs_real, XFLOAT *g_imgs_imag,
    XFLOAT *g_trans_x, XFLOAT *g_trans_y, XFLOAT *g_trans_z,
    AccProjectorKernel &projector,
    XFLOAT *g_corr_img,
    XFLOAT *g_diff2s,
    unsigned long image_size,
    XFLOAT sum_init, XFLOAT exp_local_sqrtXi2,
    unsigned long orientation_num, unsigned long translation_num,
    unsigned long todo_blocks,
    unsigned long *d_rot_idx, unsigned long *d_trans_idx,
    unsigned long *d_job_idx, unsigned long *d_job_num,
    cudaStream_t stream
) {
    if (DATA3D)
        ::CpuKernels::diff2_CC_fine_3D(
            grid_size,
            g_eulers, g_imgs_real, g_imgs_imag,
            g_trans_x, g_trans_y, g_trans_z,
            projector,
            g_corr_img, g_diff2s,
            image_size,
            sum_init,
            exp_local_sqrtXi2,
            orientation_num, translation_num,
            todo_blocks,
            d_rot_idx, d_trans_idx,
            d_job_idx, d_job_num
        );
    else
        ::CpuKernels::diff2_CC_fine_2D<REF3D>(
            grid_size,
            g_eulers, g_imgs_real, g_imgs_imag,
            g_trans_x, g_trans_y,
            projector,
            g_corr_img, g_diff2s,
            image_size, sum_init,
            exp_local_sqrtXi2,
            orientation_num, translation_num,
            todo_blocks,
            d_rot_idx, d_trans_idx,
            d_job_idx, d_job_num
        );
}
#endif

void kernel_exponentiate_weights_fine(
    int grid_size, int block_size,
    XFLOAT *g_pdf_orientation, XFLOAT *g_pdf_offset, XFLOAT *g_weights,
    unsigned long  oversamples_orient, unsigned long  oversamples_trans,
    unsigned long *d_rot_id, unsigned long *d_trans_idx,
    unsigned long *d_job_idx, unsigned long *d_job_num,
    long int job_num,
    cudaStream_t stream
);

};

namespace GpuKernels {

void centerFFT_2D(
    int grid_size, int batch_size, int block_size, cudaStream_t stream,
    XFLOAT *img_in, size_t image_size,
    int xdim, int ydim, int xshift, int yshift
);

void centerFFT_2D(
    int grid_size, int batch_size, int block_size,
    XFLOAT *img_in, size_t image_size,
    int xdim, int ydim, int xshift, int yshift
);

void centerFFT_3D(
    int grid_size, int batch_size, int block_size, cudaStream_t stream,
    XFLOAT *img_in, size_t image_size,
    int xdim, int ydim, int zdim, int xshift, int yshift, int zshift
);

#ifdef CUDA
template <typename T, acc::Type accType>
void multiply(int block_dimensions, AccDataTypes::Image<T, accType> &ptr, T value) {
    const int grid_dimensions = ceilf((float) ptr.getSize() / (float) block_dimensions);
    CudaKernels::cuda_kernel_multi<T>
    <<<grid_dimensions, block_dimensions, 0, ptr.getStream()>>>
        (ptr.getAccPtr(), value, ptr.getSize());
}

template <typename T>
void multiply(
    int grid_dimensions, int block_dimensions, cudaStream_t stream, T *array, T value, size_t size
) {
    CudaKernels::cuda_kernel_multi<T><<<grid_dimensions, block_dimensions, 0, stream>>>
        (array, value, size);
}

template <typename T, acc::Type accType>
void translate(
    int block_size, AccDataTypes::Image<T, accType> &in, AccDataTypes::Image<T, accType> &out,
    int dx, int dy, int dz = 0
) {
    if (in.getAccPtr() == out.getAccPtr())
        CRITICAL(ERRUNSAFEOBJECTREUSE);
    const int Dg = ceilf((float) in.getxyz() / (float) block_size);
    if (in.is3D()) {
        CudaKernels::cuda_kernel_translate3D<T><<<Dg, block_size, 0, in.getStream()>>>(
            in(), out(), in.getxyz(), in.getx(), in.gety(), in.getz(), dx, dy, dz
        );
    } else {
        CudaKernels::cuda_kernel_translate2D<T><<<Dg, block_size, 0, in.getStream()>>>(
            in(), out(), in.getxyz(), in.getx(), in.gety(), dx, dy
        );
    }
}

template <bool REF3D, bool DATA3D, int block_sz, int eulers_per_block, int prefetch_fraction>
void diff2_coarse(
    unsigned long grid_size, int block_size, XFLOAT *g_eulers,
    XFLOAT *trans_x, XFLOAT *trans_y, XFLOAT *trans_z,
    XFLOAT *g_real, XFLOAT *g_imag,
    AccProjectorKernel projector,
    XFLOAT *g_corr, XFLOAT *g_diff2s,
    unsigned long translation_num, unsigned long image_size,
    cudaStream_t stream
) {
    cuda_kernel_diff2_coarse<REF3D, DATA3D, block_sz, eulers_per_block, prefetch_fraction>
    <<<grid_size, block_size, 0, stream>>>(
        g_eulers, trans_x, trans_y, trans_z,
        g_real, g_imag,
        projector, g_corr, g_diff2s,
        translation_num, image_size
    );
}

template <bool REF3D, bool DATA3D, int block_sz>
void diff2_CC_coarse(
    unsigned long grid_size, int block_size,
    XFLOAT *g_eulers,
    XFLOAT *g_imgs_real, XFLOAT *g_imgs_imag,
    XFLOAT *g_trans_x, XFLOAT *g_trans_y, XFLOAT *g_trans_z,
    AccProjectorKernel projector,
    XFLOAT *g_corr_img, XFLOAT *g_diff2s,
    unsigned long translation_num, unsigned long image_size,
    XFLOAT exp_local_sqrtXi2,
    cudaStream_t stream
) {
    cuda_kernel_diff2_CC_coarse<REF3D, DATA3D, block_sz>
    <<<dim3(grid_size, translation_num), block_size, 0, stream>>>(
        g_eulers, g_imgs_real, g_imgs_imag,
        g_trans_x, g_trans_y, g_trans_z,
        projector,
        g_corr_img, g_diff2s,
        translation_num, image_size, exp_local_sqrtXi2
    );
}

template <bool REF3D, bool DATA3D, int block_sz, int chunk_sz>
void diff2_fine(
    unsigned long grid_size, int block_size,
    XFLOAT *g_eulers, XFLOAT *g_imgs_real, XFLOAT *g_imgs_imag,
    XFLOAT *trans_x, XFLOAT *trans_y, XFLOAT *trans_z,
    AccProjectorKernel projector,
    XFLOAT *g_corr_img, XFLOAT *g_diff2s,
    unsigned long image_size, XFLOAT sum_init,
    unsigned long orientation_num, unsigned long translation_num,
    unsigned long todo_blocks,
    unsigned long *d_rot_idx, unsigned long *d_trans_idx,
    unsigned long *d_job_idx, unsigned long *d_job_num,
    cudaStream_t stream
) {
    cuda_kernel_diff2_fine<REF3D, DATA3D, block_sz, chunk_sz>
    <<<grid_size, block_size, 0, stream>>>(
        g_eulers, g_imgs_real, g_imgs_imag,
        trans_x, trans_y, trans_z,
        projector,
        g_corr_img,    // in these non-CC kernels this is effectively an adjusted MinvSigma2
        g_diff2s,
        image_size, sum_init,
        orientation_num, translation_num,
        todo_blocks, //significant_num,
        d_rot_idx, d_trans_idx, d_job_idx, d_job_num
    );
}

template <bool REF3D, bool DATA3D, int block_sz,int chunk_sz>
void diff2_CC_fine(
    unsigned long grid_size, int block_size,
    XFLOAT *g_eulers, XFLOAT *g_imgs_real, XFLOAT *g_imgs_imag,
    XFLOAT *g_trans_x, XFLOAT *g_trans_y, XFLOAT *g_trans_z,
    AccProjectorKernel &projector,
    XFLOAT *g_corr_img,
    XFLOAT *g_diff2s,
    unsigned long image_size,
    XFLOAT sum_init, XFLOAT exp_local_sqrtXi2,
    unsigned long orientation_num, unsigned long translation_num,
    unsigned long todo_blocks,
    unsigned long *d_rot_idx, unsigned long *d_trans_idx,
    unsigned long *d_job_idx, unsigned long *d_job_num,
    cudaStream_t stream
) {
    dim3 block_dim = grid_size;
    cuda_kernel_diff2_CC_fine<REF3D, DATA3D, block_sz, chunk_sz>
    <<<block_dim,block_size,0,stream>>>(
        g_eulers, g_imgs_real, g_imgs_imag,
        g_trans_x, g_trans_y, g_trans_z,
        projector,
        g_corr_img, g_diff2s,
        image_size, sum_init,
        exp_local_sqrtXi2,
        orientation_num, translation_num,
        todo_blocks,
        d_rot_idx, d_trans_idx,
        d_job_idx, d_job_num
    );
}
#endif

void kernel_exponentiate_weights_fine(
    int grid_size, int block_size,
    XFLOAT *g_pdf_orientation, XFLOAT *g_pdf_offset, XFLOAT *g_weights,
    unsigned long  oversamples_orient, unsigned long  oversamples_trans,
    unsigned long *d_rot_id, unsigned long *d_trans_idx,
    unsigned long *d_job_idx, unsigned long *d_job_num,
    long int job_num,
    cudaStream_t stream
);

};

#ifdef ALTCPU
template <>
inline void frequencyPass<acc::cpu>(
    int grid_size, int block_size,
    cudaStream_t stream,
    acc::Complex *A, long int ori_size,
    size_t Xdim, size_t Ydim, size_t Zdim,
    XFLOAT edge_low, XFLOAT edge_width, XFLOAT edge_high,
    XFLOAT angpix, size_t image_size, bool do_highpass
) {
    if (do_highpass)
        ::CpuKernels::kernel_frequencyPass<1>(
            grid_size, block_size,
            A, ori_size, Xdim, Ydim, Zdim,
            edge_low, edge_width, edge_high,
            angpix, image_size
        );
    else
        ::CpuKernels::kernel_frequencyPass<0>(
            grid_size, block_size,
            A, ori_size, Xdim, Ydim, Zdim,
            edge_low, edge_width, edge_high,
            angpix, image_size
        );
}

template <>
inline void kernel_wavg<acc::cpu>(
    int block_sz, XFLOAT *g_eulers,
    AccProjectorKernel &projector,
    unsigned long image_size, unsigned long orientation_num,
    XFLOAT *g_img_real, XFLOAT *g_img_imag,
    XFLOAT *g_trans_x, XFLOAT *g_trans_y, XFLOAT *g_trans_z,
    XFLOAT* g_weights, XFLOAT* g_ctfs,
    XFLOAT *g_wdiff2s_parts, XFLOAT *g_wdiff2s_AA, XFLOAT *g_wdiff2s_XA,
    unsigned long translation_num,
    XFLOAT weight_norm, XFLOAT significant_weight, XFLOAT part_scale,
    cudaStream_t stream,
    bool CTFPREMULTIPLIED, bool REFCTF, bool REF3D, bool DATA3D
) {
    if (DATA3D) {
        #define KERNEL( A, B ) ::CpuKernels::wavg_3D<A, B>( \
            g_eulers, projector, image_size, orientation_num, \
            g_img_real, g_img_imag, \
            g_trans_x, g_trans_y, g_trans_z, \
            g_weights, g_ctfs, \
            g_wdiff2s_parts, g_wdiff2s_AA, g_wdiff2s_XA, \
            translation_num, weight_norm, significant_weight, part_scale)
        if (CTFPREMULTIPLIED) {
            if (REFCTF) KERNEL( 1, 1 ); else KERNEL( 1, 0 );
        } else {
            if (REFCTF) KERNEL( 0, 1 ); else KERNEL( 0, 0 );
        }
        #undef KERNEL
    } else {
        #define KERNEL( A, B, C ) ::CpuKernels::wavg_ref3D<A, B, C>( \
            g_eulers, projector, image_size, orientation_num, \
            g_img_real, g_img_imag, \
            g_trans_x, g_trans_y, g_trans_z, \
            g_weights, g_ctfs, \
            g_wdiff2s_parts, g_wdiff2s_AA, g_wdiff2s_XA, \
            translation_num, weight_norm, significant_weight, part_scale)
        if (CTFPREMULTIPLIED) {
            if (REFCTF) {
                if (REF3D) KERNEL( 1, 1, 1 ); else KERNEL( 1, 1, 0 );
            } else {
                if (REF3D) KERNEL( 1, 0, 1 ); else KERNEL( 1, 0, 0 );
            }
        } else {
            if (REFCTF) {
                if (REF3D) KERNEL( 0, 1, 1 ); else KERNEL( 0, 1, 0 );
            } else {
                if (REF3D) KERNEL( 0, 0, 1 ); else KERNEL( 0, 0, 0 );
            }
        }
        #undef KERNEL
    }
}
#endif

template <typename T>
void kernel_weights_exponent_coarse(
    size_t num_classes,
    AccPtr<T, acc::cpu> &g_pdf_orientation, AccPtr<bool, acc::cpu> &g_pdf_orientation_zeros,
    AccPtr<T, acc::cpu> &g_pdf_offset, AccPtr<bool, acc::cpu> &g_pdf_offset_zeros,
    AccPtr<T, acc::cpu> &g_Mweight,
    T g_min_diff2,
    size_t nr_coarse_orient, size_t  nr_coarse_trans
) {
    const size_t block_num = ceilf(((double) nr_coarse_orient * nr_coarse_trans * num_classes) / (double) SUMW_BLOCK_SIZE);
    ::CpuKernels::weights_exponent_coarse(
        g_pdf_orientation.getAccPtr(), g_pdf_orientation_zeros.getAccPtr(),
        g_pdf_offset.getAccPtr(), g_pdf_offset_zeros.getAccPtr(),
        g_Mweight.getAccPtr(),
        g_min_diff2,
        nr_coarse_orient, nr_coarse_trans,
        nr_coarse_orient * nr_coarse_trans * num_classes
    );
}

template <typename T>
void kernel_exponentiate(AccPtr<T, acc::cpu> &array, T add) {
    const int Dg = ceilf((float) array.getSize() / BLOCK_SIZE);
    ::CpuKernels::exponentiate<T>
        (array.getAccPtr(), add, array.getSize());
}

};

#endif //ACC_UTILITIES_H_
