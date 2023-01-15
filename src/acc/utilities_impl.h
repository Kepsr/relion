#ifndef ACC_UTILITIES_IMPL_H_
#define ACC_UTILITIES_IMPL_H_

#include "src/acc/acc_ptr.h"
#include "src/acc/data_types.h"
#include "src/acc/acc_helper_functions.h"
#include "src/acc/utilities.h"
#ifdef CUDA
#include "src/acc/cuda/cuda_kernels/helper.cuh"
#include "src/acc/cuda/cuda_kernels/wavg.cuh"
#include "src/acc/cuda/cuda_kernels/diff2.cuh"
#include "src/acc/cuda/cuda_fft.h"
#include "src/acc/cuda/cuda_ml_optimiser.h"
#else
#include "src/acc/cpu/cpu_kernels/helper.h"
#include "src/acc/cpu/cpu_kernels/wavg.h"
#include "src/acc/cpu/cpu_kernels/diff2.h"
#endif

void dump_array(char *name, bool *ptr, size_t size) {
    int count = 0;
    FILE *fp = fopen(name, "w");
    fprintf(fp, "Array size:  %ld\n", size);
    for (size_t i = 0; i < size; i++) {
        fprintf(fp, "%d, ", ptr[i]);
        count++;
        if (count > 10) {
            fprintf(fp, "\n");
            count = 0;
        }
    }
    fprintf(fp, "\n");
    fflush(fp);
    fclose(fp);
}

void dump_array(char *name, int *ptr, size_t size) {
    int count = 0;
    FILE *fp = fopen(name, "w");
    fprintf(fp, "Array size:  %ld\n", size);
    for (size_t i = 0; i < size; i++) {
        fprintf(fp, "%d, ", ptr[i]);
        count++;
        if (count > 10) {
            fprintf(fp, "\n");
            count = 0;
        }
    }
    fprintf(fp, "\n");
    fflush(fp);
    fclose(fp);
}

void dump_array(char *name, size_t *ptr, size_t size) {
    int count = 0;
    FILE *fp = fopen(name, "w");
    fprintf(fp, "Array size:  %ld\n", size);
    for (size_t i = 0; i < size; i++) {
        fprintf(fp, "%zu, ", ptr[i]);
        count++;
        if (count > 10) {
            fprintf(fp, "\n");
            count = 0;
        }
    }
    fprintf(fp, "\n");
    fflush(fp);
    fclose(fp);
}

void dump_array(char *name, float *ptr, size_t size) {
    int count = 0;
    FILE *fp = fopen(name, "w");
    fprintf(fp, "Array size:  %ld\n", size);
    for (size_t i = 0; i < size; i++) {
        fprintf(fp, "%f, ", ptr[i]);
        count++;
        if (count > 10) {
            fprintf(fp, "\n");
            count = 0;
        }
    }
    fprintf(fp, "\n");
    fflush(fp);
    fclose(fp);
}

void dump_complex_array(char *name, acc::Complex *ptr, size_t size) {
    int count = 0;
    FILE *fp = fopen(name, "w");
    fprintf(fp, "Array size:  %ld\n", size);
    for (size_t i = 0; i < size; i++) {
        fprintf(fp, "%f,%f, ", ptr[i].x, ptr[i].y);
        count++;
        if (count > 10) {
            fprintf(fp, "\n");
            count = 0;
        }
    }
    fprintf(fp, "\n");
    fflush(fp);
    fclose(fp);
}

void dump_complex_array(char *name, Complex *ptr, size_t size) {
    int count = 0;
    FILE *fp = fopen(name, "w");
    fprintf(fp, "Array size:  %ld\n", size);
    for (size_t i = 0; i < size; i++) {
        fprintf(fp, "%f,%f, ", ptr[i].real, ptr[i].imag);
        count++;
        if (count > 10) {
            fprintf(fp, "\n");
            count = 0;
        }
    }
    fprintf(fp, "\n");
    fflush(fp);
    fclose(fp);
}

void dump_double_array(char *name, float *ptr, float *ptr2, size_t size) {
    int count = 0;
    FILE *fp = fopen(name, "w");
    fprintf(fp, "Array size:  %ld\n", size);
    for (size_t i = 0; i < size; i++) {
        fprintf(fp, "%f,%f, ", ptr[i], ptr2[i]);
        count++;
        if (count > 10) {
            fprintf(fp, "\n");
            count = 0;
        }
    }
    fprintf(fp, "\n");
    fflush(fp);
    fclose(fp);
}

void dump_triple_array(char *name, float *ptr, float *ptr2, float *ptr3, size_t size) {
    int count = 0;
    FILE *fp = fopen(name, "w");
    fprintf(fp, "Array size:  %ld\n", size);
    for (size_t i = 0; i < size; i++) {
        fprintf(fp, "%f,%f,%f, ", ptr[i], ptr2[i], ptr3[i]);
        count++;
        if (count > 10) {
            fprintf(fp, "\n");
            count = 0;
        }
    }
    fprintf(fp, "\n");
    fflush(fp);
    fclose(fp);
}

void dump_array(char *name, double *ptr, size_t size) {
    int count = 0;
    FILE *fp = fopen(name, "w");
    fprintf(fp, "Array size:  %ld\n", size);
    for (size_t i = 0; i < size; i++) {
        fprintf(fp, "%f, ", ptr[i]);
        count++;
        if (count > 10) {
            fprintf(fp, "\n");
            count = 0;
        }
    }
    fprintf(fp, "\n");
    fflush(fp);
    fclose(fp);
}

void dump_double_array(char *name, double *ptr, double *ptr2, size_t size) {
    int count = 0;
    FILE *fp = fopen(name, "w");
    fprintf(fp, "Array size:  %ld\n", size);
    for (size_t i = 0; i < size; i++) {
        fprintf(fp, "%f,%f, ", ptr[i], ptr2[i]);
            count++;
        if (count > 10) {
            fprintf(fp, "\n");
            count = 0;
        }
    }
    fprintf(fp, "\n");
    fflush(fp);
    fclose(fp);
}

void dump_triple_array(char *name, double *ptr, double *ptr2, double *ptr3, size_t size) {
    int count = 0;
    FILE *fp = fopen(name, "w");
    fprintf(fp, "Array size:  %ld\n", size);
    for (size_t i = 0; i < size; i++) {
        fprintf(fp, "%f,%f,%f, ", ptr[i], ptr2[i], ptr3[i]);
        count++;
        if (count > 10) {
            fprintf(fp, "\n");
            count = 0;
        }
    }
    fprintf(fp, "\n");
    fflush(fp);
    fclose(fp);
}

namespace AccUtilities {

#ifdef ALTCPU
void makeNoiseImage(
    XFLOAT sigmaFudgeFactor, MultidimArray<RFLOAT> &sigmaNoiseSpectra,
    long int seed, MlOptimiserCpu *accMLO, AccPtr<XFLOAT> &RandomImage, bool is3D
) {
    // Different MPI-distributed subsets may otherwise have different instances of the random noise below,
    // because work is on an on-demand basis and therefore variable with the timing of distinct nodes...
    // Have the seed based on the part_id, so that each particle has a different instant of the noise
    init_random_generator(seed);

    // Make a holder for the spectral profile and copy to the GPU
    // AccDataTypes::Image<XFLOAT> NoiseSpectra(sigmaNoiseSpectra, ptrFactory);
    AccPtr<XFLOAT> NoiseSpectra = RandomImage.make<XFLOAT>(sigmaNoiseSpectra.size());
    NoiseSpectra.allAlloc();

    for (int n = 0; n < sigmaNoiseSpectra.size(); n++)
        NoiseSpectra.getHostPtr()[n] = sqrt(sigmaFudgeFactor * sigmaNoiseSpectra.data[n]);

    // Create noise image with the correct spectral profile
    if (is3D) {
        ::CpuKernels::RNDnormalDitributionComplexWithPowerModulation3D(
            accMLO->transformer1.fouriers.getAccPtr(),
            accMLO->transformer1.sizef[0], accMLO->transformer1.sizef[1],
            NoiseSpectra.getAccPtr());
    } else {
        ::CpuKernels::RNDnormalDitributionComplexWithPowerModulation2D(
            accMLO->transformer1.fouriers.getAccPtr(),
            accMLO->transformer1.sizef[0],
            NoiseSpectra.getAccPtr());
    }

    // Transform to real-space, to get something which look like
    // the particle image without actual signal (a particle)
    accMLO->transformer1.backward();

    // Copy the randomized image to A separate device-array, so that the
    // transformer can be used to set up the actual particle image
    std::copy_n(accMLO->transformer1.reals.getHostPtr(),
                RandomImage.getSize(),
                RandomImage.getHostPtr());
}

void TranslateAndNormCorrect(
    MultidimArray<RFLOAT> &img_in, AccPtr<XFLOAT, acc::cpu> &img_out,
    XFLOAT normcorr, RFLOAT xOff, RFLOAT yOff, RFLOAT zOff, bool DATA3D
) {
    const size_t N = img_in.size();
    // Temporary array because translate is out-of-place
    AccPtr<XFLOAT> temp = img_out.make_like(img_in.data, N);
    // Apply the norm_correction term
    if (normcorr != 1) {
        ::CpuKernels::cpu_kernel_multi<XFLOAT>
        (temp.getAccPtr(), normcorr, N);
    }
    // LAUNCH_PRIVATE_ERROR(cudaGetLastError(),accMLO->errorStatus);
    if (temp.getAccPtr() == img_out.getAccPtr())
        CRITICAL(ERRUNSAFEOBJECTREUSE);
    if (DATA3D)
        ::CpuKernels::cpu_translate3D<XFLOAT>
        (temp.getAccPtr(), img_out.getAccPtr(),
        img_in.xdim * img_in.ydim * img_in.zdim,
        img_in.xdim, img_in.ydim, img_in.zdim, xOff, yOff, zOff);
    else
        ::CpuKernels::cpu_translate2D<XFLOAT>
        (temp.getAccPtr(), img_out.getAccPtr(),
        img_in.xdim * img_in.ydim * img_in.zdim,
        img_in.xdim, img_in.ydim, xOff, yOff);
}

void normalizeAndTransformImage(
    AccPtr<XFLOAT, acc::cpu> &img_in, MultidimArray<Complex> &img_out, MlOptimiserCpu *accMLO,
    size_t xSize, size_t ySize, size_t zSize
) {
    img_in.cpOnAcc(accMLO->transformer1.reals.getDevicePtr());
    runCenterFFT(
        accMLO->transformer1.reals,
        (int) accMLO->transformer1.sizer[0],
        (int) accMLO->transformer1.sizer[1],
        (int) accMLO->transformer1.sizer[2],
        false
    );
    accMLO->transformer1.reals.streamSync();
    accMLO->transformer1.forward();
    accMLO->transformer1.fouriers.streamSync();

    size_t FMultiBsize = ((int) ceilf((float) accMLO->transformer1.fouriers.getSize() * 2 / (float) BLOCK_SIZE));
    CpuKernels::multiply<XFLOAT>(
        FMultiBsize, BLOCK_SIZE, accMLO->transformer1.fouriers.getStream(),
        (XFLOAT*) accMLO->transformer1.fouriers.getAccPtr(),
        (XFLOAT) 1 / (XFLOAT) accMLO->transformer1.reals.getSize(),
        accMLO->transformer1.fouriers.getSize() * 2
    );
    // LAUNCH_PRIVATE_ERROR(cudaGetLastError(),accMLO->errorStatus);

    AccPtr<acc::Complex> d_Fimg = img_in.make<acc::Complex>(xSize * ySize * zSize);
    d_Fimg.allAlloc();
    accMLO->transformer1.fouriers.streamSync();
    windowFourierTransform2(
        accMLO->transformer1.fouriers,
        d_Fimg,
        accMLO->transformer1.sizef[0], accMLO->transformer1.sizef[1], accMLO->transformer1.sizef[2], //Input dimensions
        xSize, ySize, zSize  //Output dimensions
    );
    accMLO->transformer1.fouriers.streamSync();

    d_Fimg.cpToHost();
    d_Fimg.streamSync();
    img_out.resize(xSize, ySize, zSize);
    for (unsigned long i = 0; i < img_out.size(); i++) {
        img_out.data[i].real = d_Fimg.getHostPtr()[i].x;
        img_out.data[i].imag = d_Fimg.getHostPtr()[i].y;
    }
}

void softMaskBackgroundValue(
    AccDataTypes::Image<XFLOAT, acc::cpu> &vol,
    XFLOAT radius, XFLOAT radius_p, XFLOAT cosine_width,
    AccPtr<XFLOAT, acc::cpu> &g_sum, AccPtr<XFLOAT, acc::cpu> &g_sum_bg
) {
    const int grid_dimensions = 128;  // TODO: set balanced (hardware-dependent?)
    ::CpuKernels::softMaskBackgroundValue(
        grid_dimensions, SOFTMASK_BLOCK_SIZE,
        vol.ptr.getAccPtr(), vol.getxyz(),
        vol.getx(), vol.gety(), vol.getz(),
        vol.getx() / 2, vol.gety() / 2, vol.getz() / 2,
        radius, radius_p,
        cosine_width, g_sum.getAccPtr(), g_sum_bg.getAccPtr()
    );
}

void cosineFilter(
    AccDataTypes::Image<XFLOAT, acc::cpu> &vol,
    bool do_Mnoise, AccDataTypes::Image<XFLOAT, acc::cpu> Noise,
    XFLOAT radius, XFLOAT radius_p, XFLOAT cosine_width,
    XFLOAT sum_bg_total
) {
    const int grid_dimensions = 128;  // TODO: set balanced (hardware-dependent?)
    ::CpuKernels::cosineFilter(
        grid_dimensions, SOFTMASK_BLOCK_SIZE,
        vol.ptr.getAccPtr(), vol.getxyz(),
        vol.getx(), vol.gety(), vol.getz(),
        vol.getx() / 2, vol.gety() / 2, vol.getz() / 2,
        !do_Mnoise, Noise.ptr.getAccPtr(),
        radius, radius_p,
        cosine_width, sum_bg_total
    );
}
#endif

#ifdef CUDA
void makeNoiseImage(
    XFLOAT sigmaFudgeFactor, MultidimArray<RFLOAT> &sigmaNoiseSpectra,
    long int seed, MlOptimiserCuda *accMLO, AccPtr<XFLOAT> &RandomImage, bool is3D
) {
    // Different MPI-distributed subsets may otherwise have different instances of the random noise below,
    // because work is on an on-demand basis and therefore variable with the timing of distinct nodes...
    // Have the seed based on the part_id, so that each particle has a different instant of the noise
    init_random_generator(seed);

    // Make a holder for the spectral profile and copy to the GPU
    // AccDataTypes::Image<XFLOAT> NoiseSpectra(sigmaNoiseSpectra, ptrFactory);
    AccPtr<XFLOAT> NoiseSpectra = RandomImage.make<XFLOAT>(sigmaNoiseSpectra.size());
    NoiseSpectra.allAlloc();
    for (int n = 0; n < sigmaNoiseSpectra.size(); n++)
        NoiseSpectra.getHostPtr()[n] = sqrt(sigmaFudgeFactor * sigmaNoiseSpectra.data[n]);

    // Set up states to seeda and run randomization on the GPU
    // AccDataTypes::Image<curandState > RandomStates(RND_BLOCK_NUM*RND_BLOCK_SIZE,ptrFactory);
    AccPtr<curandState> RandomStates = RandomImage.make<curandState>(RND_BLOCK_NUM * RND_BLOCK_SIZE);
    RandomStates.deviceAlloc();

    NoiseSpectra.cpToDevice();
    NoiseSpectra.streamSync();
    LAUNCH_PRIVATE_ERROR(cudaGetLastError(), accMLO->errorStatus);

    // Initialize randomization by particle ID, like on the CPU-side
    cuda_kernel_initRND<<<RND_BLOCK_NUM, RND_BLOCK_SIZE>>>
        (seed, RandomStates.getAccPtr());
    LAUNCH_PRIVATE_ERROR(cudaGetLastError(), accMLO->errorStatus);

    // Create noise image with the correct spectral profile
    if (is3D) {
        cuda_kernel_RNDnormalDitributionComplexWithPowerModulation3D
        <<<RND_BLOCK_NUM, RND_BLOCK_SIZE>>>(
            accMLO->transformer1.fouriers.getAccPtr(),
            RandomStates.getAccPtr(),
            accMLO->transformer1.sizef[0],
            accMLO->transformer1.sizef[1],
            NoiseSpectra.getAccPtr());
    } else {
        cuda_kernel_RNDnormalDitributionComplexWithPowerModulation2D
        <<<RND_BLOCK_NUM, RND_BLOCK_SIZE>>>(
            accMLO->transformer1.fouriers.getAccPtr(),
            RandomStates.getAccPtr(),
            accMLO->transformer1.sizef[0],
            NoiseSpectra.getAccPtr());
    }
    LAUNCH_PRIVATE_ERROR(cudaGetLastError(), accMLO->errorStatus);

    // Transform to real-space, to get something which look like
    // the particle image without actual signal (a particle)
    accMLO->transformer1.backward();

    // Copy the randomized image to A separate device-array, so that the
    // transformer can be used to set up the actual particle image
    accMLO->transformer1.reals.cpOnDevice(RandomImage.getAccPtr());
    // cudaMLO->transformer1.reals.streamSync();
}

void TranslateAndNormCorrect(
    MultidimArray<RFLOAT> &img_in, AccPtr<XFLOAT, acc::cuda> &img_out,
    XFLOAT normcorr, RFLOAT xOff, RFLOAT yOff, RFLOAT zOff, bool DATA3D
) {
    const size_t N = img_in.size();
    // Temporary array because translate is out-of-place
    AccPtr<XFLOAT> temp = img_out.make_like(img_in.data, N);
    // Apply the norm_correction term
    const int grid_dimensions = ceilf((float) N / (float) BLOCK_SIZE);
    if (normcorr != 1) {
        ::CudaKernels::cuda_kernel_multi<XFLOAT>
        <<<grid_dimensions, BLOCK_SIZE, 0, temp.getStream()>>>
        (temp.getAccPtr(), normcorr, N);
    }
    // LAUNCH_PRIVATE_ERROR(cudaGetLastError(),accMLO->errorStatus);

    if (temp.getAccPtr() == img_out.getAccPtr())
        CRITICAL(ERRUNSAFEOBJECTREUSE);
    if (DATA3D)
        ::CudaKernels::cuda_kernel_translate3D<XFLOAT>
        <<<grid_dimensions, BLOCK_SIZE, 0, temp.getStream()>>>
        (temp.getAccPtr(), img_out.getAccPtr(),
        img_in.xdim * img_in.ydim * img_in.zdim,
        img_in.xdim, img_in.ydim, img_in.zdim, xOff, yOff, zOff);
    else
        ::CudaKernels::cuda_kernel_translate2D<XFLOAT>
        <<<grid_dimensions, BLOCK_SIZE, 0, temp.getStream()>>>
        (temp.getAccPtr(), img_out.getAccPtr(),
        img_in.xdim * img_in.ydim * img_in.zdim,  // BUG? should be no zdim?
        img_in.xdim, img_in.ydim, xOff, yOff);
    // LAUNCH_PRIVATE_ERROR(cudaGetLastError(), accMLO->errorStatus);
}

void normalizeAndTransformImage(
    AccPtr<XFLOAT, acc::cuda> &img_in, MultidimArray<Complex> &img_out, MlOptimiserCuda *accMLO,
    size_t xSize, size_t ySize, size_t zSize
) {
    img_in.cpOnAcc(accMLO->transformer1.reals.getDevicePtr());
    runCenterFFT(
        accMLO->transformer1.reals,
        (int) accMLO->transformer1.sizer[0],
        (int) accMLO->transformer1.sizer[1],
        (int) accMLO->transformer1.sizer[2],
        false
    );
    accMLO->transformer1.reals.streamSync();
    accMLO->transformer1.forward();
    accMLO->transformer1.fouriers.streamSync();

    size_t FMultiBsize = ((int) ceilf((float) accMLO->transformer1.fouriers.getSize() * 2 / (float) BLOCK_SIZE));
    GpuKernels::multiply(
        FMultiBsize, BLOCK_SIZE, accMLO->transformer1.fouriers.getStream(),
        (XFLOAT*) accMLO->transformer1.fouriers.getAccPtr(),
        (XFLOAT) 1 / (XFLOAT) accMLO->transformer1.reals.getSize(),
        accMLO->transformer1.fouriers.getSize() * 2
    );
    // LAUNCH_PRIVATE_ERROR(cudaGetLastError(),accMLO->errorStatus);

    AccPtr<acc::Complex> d_Fimg = img_in.make<acc::Complex>(xSize * ySize * zSize);
    d_Fimg.allAlloc();
    accMLO->transformer1.fouriers.streamSync();
    windowFourierTransform2(
        accMLO->transformer1.fouriers,
        d_Fimg,
        accMLO->transformer1.sizef[0], accMLO->transformer1.sizef[1], accMLO->transformer1.sizef[2], //Input dimensions
        xSize, ySize, zSize  //Output dimensions
    );
    accMLO->transformer1.fouriers.streamSync();

    d_Fimg.cpToHost();
    d_Fimg.streamSync();
    img_out.resize(xSize, ySize, zSize);
    for (unsigned long i = 0; i < img_out.size(); i++) {
        img_out.data[i].real = d_Fimg.getHostPtr()[i].x;
        img_out.data[i].imag = d_Fimg.getHostPtr()[i].y;
    }
}

void softMaskBackgroundValue(
    AccDataTypes::Image<XFLOAT, acc::cuda> &vol,
    XFLOAT radius, XFLOAT radius_p, XFLOAT cosine_width,
    AccPtr<XFLOAT, acc::cuda> &g_sum, AccPtr<XFLOAT, acc::cuda> &g_sum_bg
) {
    const int grid_dimensions = 128;  // TODO: set balanced (hardware-dependent?)
    cuda_kernel_softMaskBackgroundValue
    <<<grid_dimensions, SOFTMASK_BLOCK_SIZE, 0, vol.ptr.getStream()>>>(
        {vol.ptr.getAccPtr(), vol.getxyz(),
         vol.getx(),     vol.gety(),     vol.getz(),
         vol.getx() / 2, vol.gety() / 2, vol.getz() / 2},
        radius, radius_p, cosine_width, g_sum.getAccPtr(), g_sum_bg.getAccPtr()
    );
}

void cosineFilter(
    AccDataTypes::Image<XFLOAT, acc::cuda> &vol,
    bool do_Mnoise, AccDataTypes::Image<XFLOAT, acc::cuda> Noise,
    XFLOAT radius, XFLOAT radius_p, XFLOAT cosine_width,
    XFLOAT sum_bg_total
) {
    const int grid_dimensions = 128;  // TODO: set balanced (hardware-dependent?)
    cuda_kernel_cosineFilter
    <<<grid_dimensions, SOFTMASK_BLOCK_SIZE, 0, vol.ptr.getStream()>>>(
        {vol.ptr.getAccPtr(), vol.getxyz(),
         vol.getx(),     vol.gety(),     vol.getz(),
         vol.getx() / 2, vol.gety() / 2, vol.getz() / 2},
        !do_Mnoise, Noise.ptr.getAccPtr(),
        radius, radius_p,
        cosine_width, sum_bg_total
    );
}
#endif

void initOrientations(
    AccPtr<RFLOAT, acc::cpu> &pdfs, AccPtr<XFLOAT, acc::cpu> &pdf_orientation,
    AccPtr<bool, acc::cpu> &pdf_orientation_zeros
) {
    for (int i = 0; i < pdfs.getSize(); i++) {
        const RFLOAT x = pdfs.getHostPtr()[i];
        pdf_orientation      .getHostPtr()[i] = x == 0.0 ? 0.0 : log(x);
        pdf_orientation_zeros.getHostPtr()[i] = x == 0.0;
    }
}

#ifdef CUDA
void initOrientations(
    AccPtr<RFLOAT, acc::cuda> &pdfs, AccPtr<XFLOAT, acc::cuda> &pdf_orientation,
    AccPtr<bool, acc::cuda> &pdf_orientation_zeros
) {
    const int Db = 512, Dg = ceil(pdfs.getSize() / (float) Db);
    cuda_kernel_initOrientations<<<Dg, Db, 0, pdfs.getStream()>>>
    (pdfs.getAccPtr(), pdf_orientation.getAccPtr(),
     pdf_orientation_zeros.getAccPtr(), pdfs.getSize());
    LAUNCH_HANDLE_ERROR(cudaGetLastError());
}
#endif

#ifdef ALTCPU
namespace CpuKernels {

void centerFFT_2D(
    int grid_size, int batch_size, int block_size, cudaStream_t stream,
    XFLOAT *img_in, size_t image_size, int xdim, int ydim, int xshift, int yshift
) {
    ::CpuKernels::centerFFT_2D<XFLOAT>(
        batch_size, 0, image_size / 2,
        img_in, image_size, xdim, ydim, xshift, yshift);
}

void centerFFT_2D(
    int grid_size, int batch_size, int block_size,
    XFLOAT *img_in, size_t image_size, int xdim, int ydim, int xshift, int yshift
) {
    ::CpuKernels::centerFFT_2D<XFLOAT>(
        batch_size, 0, image_size / 2,
        img_in, image_size, xdim, ydim, xshift, yshift);
}

void centerFFT_3D(
    int grid_size, int batch_size, int block_size, cudaStream_t stream,
    XFLOAT *img_in, size_t image_size, int xdim, int ydim, int zdim, int xshift, int yshift, int zshift
) {
    ::CpuKernels::centerFFT_3D<XFLOAT>(
        batch_size, 0, image_size / 2,
        img_in, image_size, xdim, ydim, zdim, xshift, yshift, zshift);
}

void kernel_exponentiate_weights_fine(
    XFLOAT *g_pdf_orientation, bool *g_pdf_orientation_zeros,
    XFLOAT *g_pdf_offset, bool *g_pdf_offset_zeros,
    XFLOAT *g_weights, XFLOAT min_diff2,
    unsigned long oversamples_orient, unsigned long oversamples_trans,
    unsigned long *d_rot_id, unsigned long *d_trans_idx,
    unsigned long *d_job_idx, unsigned long *d_job_num,
    long int job_num, cudaStream_t stream
) {
    ::CpuKernels::exponentiate_weights_fine(
        g_pdf_orientation, g_pdf_orientation_zeros,
        g_pdf_offset, g_pdf_offset_zeros,
        g_weights, min_diff2,
        oversamples_orient, oversamples_trans,
        d_rot_id, d_trans_idx,
        d_job_idx, d_job_num, job_num
    );
}

};
#endif

#ifdef CUDA
namespace GpuKernels {

void centerFFT_2D(
    int grid_size, int batch_size, int block_size, cudaStream_t stream,
    XFLOAT *img_in, size_t image_size, int xdim, int ydim, int xshift, int yshift
) {
    const dim3 blocks (grid_size, batch_size);
    cuda_kernel_centerFFT_2D<<<blocks, block_size, 0, stream>>>(
        img_in, image_size, xdim, ydim, xshift, yshift);
}

void centerFFT_2D(
    int grid_size, int batch_size, int block_size,
    XFLOAT *img_in, size_t image_size, int xdim, int ydim, int xshift, int yshift
) {
    const dim3 blocks (grid_size, batch_size);
    cuda_kernel_centerFFT_2D<<<blocks, block_size>>>(
        img_in, image_size, xdim, ydim, xshift, yshift);
}

void centerFFT_3D(
    int grid_size, int batch_size, int block_size, cudaStream_t stream,
    XFLOAT *img_in, size_t image_size, int xdim, int ydim, int zdim, int xshift, int yshift, int zshift
) {
    const dim3 blocks (grid_size, batch_size);
    cuda_kernel_centerFFT_3D<<<blocks, block_size, 0, stream>>>(
        img_in, image_size, xdim, ydim, zdim, xshift, yshift, zshift);
}

void kernel_exponentiate_weights_fine(
    XFLOAT *g_pdf_orientation, bool *g_pdf_orientation_zeros,
    XFLOAT *g_pdf_offset, bool *g_pdf_offset_zeros,
    XFLOAT *g_weights, XFLOAT min_diff2,
    unsigned long oversamples_orient, unsigned long oversamples_trans,
    unsigned long *d_rot_id, unsigned long *d_trans_idx,
    unsigned long *d_job_idx, unsigned long *d_job_num,
    long int job_num, cudaStream_t stream
) {
    const size_t block_num = ceil((double) job_num / (double) SUMW_BLOCK_SIZE);
    cuda_kernel_exponentiate_weights_fine<<<block_num, SUMW_BLOCK_SIZE, 0, stream>>>(
        g_pdf_orientation, g_pdf_orientation_zeros,
        g_pdf_offset, g_pdf_offset_zeros,
        g_weights, min_diff2,
        oversamples_orient, oversamples_trans,
        d_rot_id, d_trans_idx,
        d_job_idx, d_job_num, job_num
    );
}

};
#endif

};

#ifdef CUDA
void run_padTranslatedMap(
    RFLOAT *d_in, RFLOAT *d_out,
    size_t isX, size_t ieX, size_t isY, size_t ieY, size_t isZ, size_t ieZ,  // Input dimensions
    size_t osX, size_t oeX, size_t osY, size_t oeY, size_t osZ, size_t oeZ,  // Output dimensions
    cudaStream_t stream
) {
    const size_t iszX = ieX - isX + 1;
    const size_t iszY = ieY - isY + 1;
    const size_t iszZ = ieZ - isZ + 1;
    const size_t oszX = oeX - osX + 1;
    const size_t oszY = oeY - osY + 1;
    const size_t oszZ = oeZ - osZ + 1;
    if (iszX == oszX && iszY == oszY && iszZ == oszZ) {
        cudaCpyDeviceToDevice(d_in, d_out, iszX*iszY*iszZ, stream);
    } else {
        const dim3 block_dim (16, 4, 2);
        const dim3 grid_dim (ceil(oszX / (float) block_dim.x),
                             ceil(oszY / (float) block_dim.y),
                             ceil(oszZ / (float) block_dim.z));
        cuda_kernel_window_transform<RFLOAT><<<grid_dim, block_dim, 0, stream>>>(
            d_in, d_out,
            iszX, iszY, iszZ, // Input dimensions
            isX - osX, isY - osY, isZ - osZ, oszX, oszY, oszZ  // Output dimensions
        );
        LAUNCH_HANDLE_ERROR(cudaGetLastError());
    }
}

void run_griddingCorrect(
    RFLOAT *vol, int interpolator, 
    RFLOAT rrval, RFLOAT r_min_nn,
    size_t iX, size_t iY, size_t iZ
) {
    const dim3 Db (32, 4, 2);
    const dim3 Dg (ceil(iX / (float) Db.x),
                   ceil(iY / (float) Db.y),
                   ceil(iZ / (float) Db.z));
    cuda_kernel_griddingCorrect<<<Dg, Db>>>(vol, interpolator, rrval, r_min_nn, iX, iY, iZ);
    LAUNCH_HANDLE_ERROR(cudaGetLastError());
}

void run_CenterFFTbySign(
    Complex *img_in, int xSize, int ySize, int zSize, cudaStream_t stream
) {
    const dim3 Db (32, 4, 2);
    const dim3 Dg (ceil(xSize / (float) Db.x),
                   ceil(ySize / (float) Db.y),
                   ceil(zSize / (float) Db.z));
    using rfloat2 = std::conditional<std::is_same<RFLOAT, double>::value, double2, float2>::type;
    cuda_kernel_centerFFTbySign<<<Dg, Db, 0, stream>>>
    ((rfloat2*) img_in, xSize, ySize, zSize);
    LAUNCH_HANDLE_ERROR(cudaGetLastError());
}
#endif

#endif //ACC_UTILITIES_H_
