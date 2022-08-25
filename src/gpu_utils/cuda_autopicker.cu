/**
 * WARNING:
 * This file (src/gpu_utils/cuda_autopicker.cu) is not referenced by any other,
 * and is a doublet of src/acc/cuda/cuda_autopicker.cu
 *
 * ('Doublet' in the sense that it implements functions
 * of the same name as those implemented in src/acc/cuda/cuda_autopicker.cu,
 * but those functions may not behave exactly the same.)
*/
#include <sys/time.h>
#include <stdio.h>
#include <time.h>
#include <math.h>
#include <ctime>
#include <iostream>
#include <fstream>
#include <cuda_runtime.h>
#include <signal.h>
#include "src/gpu_utils/cuda_autopicker.h"

#include "src/gpu_utils/cuda_mem_utils.h"
#include "src/gpu_utils/cuda_projector.h"
#include "src/gpu_utils/cuda_settings.h"
#include "src/gpu_utils/cuda_benchmark_utils.h"
#include "src/gpu_utils/cuda_helper_functions.cuh"
#include "src/gpu_utils/cuda_fft.h"

#include "src/macros.h"
#include "src/error.h"

#ifdef CUDA_FORCESTL
#include "src/gpu_utils/cuda_utils_stl.cuh"
#else
#include "src/gpu_utils/cuda_utils_cub.cuh"
#endif

// Z-score of x given mean mu and standard deviation sigma
inline RFLOAT Z(RFLOAT x, RFLOAT mu, RFLOAT sigma) {
    return (x - mu) / sigma;
}

AutoPickerCuda::AutoPickerCuda(AutoPicker *basePicker, int dev_id, const char * timing_fnm):
    node(NULL), basePckr(basePicker),
    allocator(new CudaCustomAllocator(0, 1)),
    micTransformer(0, allocator),
    cudaTransformer1(0, allocator),
    #ifdef TIMING_FILES
    timer(timing_fnm),
    #endif
    cudaTransformer2(0, allocator) {

    cudaProjectors.resize(basePckr->Mrefs.size());
    have_warned_batching = false;
    /*======================================================
                        DEVICE SETTINGS
    ======================================================*/
    device_id = dev_id;
    int devCount;
    HANDLE_ERROR(cudaGetDeviceCount(&devCount));

    if (dev_id >= devCount) {
        // std::cerr << " using device_id=" << dev_id << " (device no. " << dev_id + 1 << ") which is higher than the available number of devices=" << devCount << std::endl;
        CRITICAL(ERR_GPUID);
    } else {
        HANDLE_ERROR(cudaSetDevice(dev_id));
    }
};

AutoPickerCuda::AutoPickerCuda(AutoPickerMpi *basePicker, int dev_id, const char *timing_fnm):
    basePckr(basePicker),
    allocator(new CudaCustomAllocator(0, 1)),
    micTransformer(0, allocator),
    cudaTransformer1(0, allocator),
    #ifdef TIMING_FILES
    timer(timing_fnm),
    #endif
    cudaTransformer2(0, allocator) {

    node = basePicker->getNode();
    /// BUG: class "MpiNode" has no member "isMaster"
    basePicker->verb = node->isMaster() ? 1 : 0;

    cudaProjectors.resize(basePckr->Mrefs.size());
    have_warned_batching = false;
    /*======================================================
                        DEVICE SETTINGS
    ======================================================*/
    device_id = dev_id;
    int devCount;
    HANDLE_ERROR(cudaGetDeviceCount(&devCount));

    if (dev_id >= devCount) {
        // std::cerr << " using device_id=" << dev_id << " (device no. " << dev_id + 1 << ") which is higher than the available number of devices=" << devCount << std::endl;
        CRITICAL(ERR_GPUID);
    } else {
        HANDLE_ERROR(cudaSetDevice(dev_id));
    }
};

void AutoPickerCuda::run() {
    long int my_first_micrograph, my_last_micrograph, my_nr_micrographs;
    if (node != NULL) {
        // Each node does part of the work
        divide_equally(basePckr->fn_micrographs.size(), node->size, node->rank, my_first_micrograph, my_last_micrograph);
    } else {
        my_first_micrograph = 0;
        my_last_micrograph = basePckr->fn_micrographs.size() - 1;
    }
    my_nr_micrographs = my_last_micrograph - my_first_micrograph + 1;

    int barstep;
    if (basePckr->verb > 0) {
        std::cout << " Autopicking ..." << std::endl;
        init_progress_bar(my_nr_micrographs);
        barstep = std::max(1, (int) my_nr_micrographs / 60);
    }

    if (!basePckr->do_read_fom_maps) {
        CTICTOC(timer, "setupProjectors", ({
        for (int iref = 0; iref < basePckr->Mrefs.size(); iref++) {
            cudaProjectors[iref].setMdlDim(
                basePckr->PPref[iref].data.xdim,
                basePckr->PPref[iref].data.ydim,
                basePckr->PPref[iref].data.zdim,
                basePckr->PPref[iref].data.yinit,
                basePckr->PPref[iref].data.zinit,
                basePckr->PPref[iref].r_max,
                basePckr->PPref[iref].padding_factor
            );
            cudaProjectors[iref].initMdl(&basePckr->PPref[iref].data.data[0]);
        }
        }))
    }

    FileName fn_olddir = "";

    for (long int imic = my_first_micrograph; imic <= my_last_micrograph; imic++) {
        if (basePckr->verb > 0 && imic % barstep == 0)
            progress_bar(imic);


        // Check new-style outputdirectory exists and make it if not!
        FileName fn_dir = basePckr->getOutputRootName(basePckr->fn_micrographs[imic]);
        fn_dir = fn_dir.beforeLastOf("/");
        if (fn_dir != fn_olddir) {
            // Make a Particles directory
            system(("mkdir -p " + fn_dir).c_str());
            fn_olddir = fn_dir;
        }
        #ifdef TIMING
        basePckr->timer.tic(basePckr->TIMING_A5);
        #endif
        autoPickOneMicrograph(basePckr->fn_micrographs[imic], imic);
    }
    /// BUG: Should this toc be inside the above loop?
        #ifdef TIMING
        basePckr->timer.toc(basePckr->TIMING_A5);
        #endif
    if (basePckr->verb > 0)
        progress_bar(my_nr_micrographs);

    cudaDeviceReset();

}

void AutoPickerCuda::calculateStddevAndMeanUnderMask(
    CudaGlobalPtr<CUDACOMPLEX> &d_Fmic, CudaGlobalPtr<CUDACOMPLEX> &d_Fmic2, CudaGlobalPtr<CUDACOMPLEX> &d_Fmsk,
    int nr_nonzero_pixels_mask,
    CudaGlobalPtr<XFLOAT> &d_Mstddev, CudaGlobalPtr<XFLOAT> &d_Mmean,
    size_t x, size_t y, size_t mic_size, size_t workSize
) {
    cudaTransformer2.setSize(workSize, workSize, 1);

    deviceInitValue(d_Mstddev, (XFLOAT) 0.0);

    RFLOAT normfft = (RFLOAT) (mic_size * mic_size) / (RFLOAT) nr_nonzero_pixels_mask;

    CudaGlobalPtr<CUDACOMPLEX> d_Fcov(d_Fmic.getAllocator());
    d_Fcov.device_alloc(d_Fmic.getSize());

    CTICTOC(timer, "PRE-multi_0", ({
    int Bsize = ceilf((float) d_Fmic.size / (float) BLOCK_SIZE);
    cuda_kernel_convol_B<<<Bsize, BLOCK_SIZE>>>(
        ~d_Fmic, ~d_Fmsk, ~d_Fcov, d_Fmic.getSize()
    );
    LAUNCH_HANDLE_ERROR(cudaGetLastError());
    }))

    CTICTOC(timer, "PRE-window_0", ({
    windowFourierTransform2(
        d_Fcov,
        cudaTransformer2.fouriers,
        x, y, 1,
        workSize / 2 + 1, workSize, 1
    );
    }))

    CTICTOC(timer, "PRE-Transform_0", ({ cudaTransformer2.backward(); }))

    Bsize = ceilf((float) cudaTransformer2.reals.size / (float) BLOCK_SIZE);
    cuda_kernel_multi<<<Bsize, BLOCK_SIZE>>>(
        cudaTransformer2.reals.d_ptr,
        cudaTransformer2.reals.d_ptr,
        (XFLOAT) normfft,
        cudaTransformer2.reals.size
    );
    LAUNCH_HANDLE_ERROR(cudaGetLastError());

    CTICTOC(timer, "PRE-multi_1", ({
    cuda_kernel_multi<<<Bsize, BLOCK_SIZE>>>(
        cudaTransformer2.reals.d_ptr,
        cudaTransformer2.reals.d_ptr,
        d_Mstddev.d_ptr,
        (XFLOAT) -1,
        cudaTransformer2.reals.size
    );
    LAUNCH_HANDLE_ERROR(cudaGetLastError());
    }))

    CTICTOC(timer, "PRE-CenterFFT_0", ({
    runCenterFFT(
        cudaTransformer2.reals,
        (int) cudaTransformer2.xSize,
        (int) cudaTransformer2.ySize,
        false,
        1
    );
    }))

    cudaTransformer2.reals.cp_on_device(d_Mmean); //TODO remove the need for this

    CTICTOC(timer, "PRE-multi_2", ({
    Bsize = ((int) ceilf((float) d_Fmsk.size / (float) BLOCK_SIZE));
    cuda_kernel_convol_A<<<Bsize, BLOCK_SIZE>>>(
        ~d_Fmsk, ~d_Fmic2, ~d_Fcov, d_Fmsk.size
    );
    LAUNCH_HANDLE_ERROR(cudaGetLastError());
    }))

    CTICTOC(timer, "PRE-window_1", ({
    windowFourierTransform2(
        d_Fcov,
        cudaTransformer2.fouriers,
        x, y, 1,
        workSize / 2 + 1, workSize, 1
    );
    }))

    CTICTOC(timer, "PRE-Transform_1", ({ cudaTransformer2.backward(); }))

    CTICTOC(timer, "PRE-multi_3", ({
    Bsize = ceilf((float) d_Mstddev.size / (float) BLOCK_SIZE);
    cuda_kernel_finalizeMstddev<<<Bsize, BLOCK_SIZE>>>(
        d_Mstddev.d_ptr,
        cudaTransformer2.reals.d_ptr,
        normfft,
        d_Mstddev.size
    );
    LAUNCH_HANDLE_ERROR(cudaGetLastError());
    }))

    CTICTOC(timer, "PRE-CenterFFT_1", ({
    runCenterFFT(d_Mstddev, (int) workSize, (int) workSize, false, 1);
    }))

}

void AutoPickerCuda::autoPickOneMicrograph(FileName &fn_mic, long int imic) {
    Image<RFLOAT> Imic;
    MultidimArray<Complex> Faux, Faux2, Fmic;
    MultidimArray<RFLOAT> Maux, Mstddev, Mccf_best, Mpsi_best, Fctf, Mccf_best_combined;
    MultidimArray<int> Mclass_best_combined;

    CudaGlobalPtr<XFLOAT> d_Mccf_best(basePckr->workSize * basePckr->workSize, allocator);
    CudaGlobalPtr<XFLOAT> d_Mpsi_best(basePckr->workSize * basePckr->workSize, allocator);
    d_Mccf_best.device_alloc();
    d_Mpsi_best.device_alloc();

    // Always use the same random seed
    init_random_generator(basePckr->random_seed + imic);

    RFLOAT sum_ref_under_circ_mask, sum_ref2_under_circ_mask;
    int my_skip_side = basePckr->autopick_skip_side + basePckr->particle_size / 2;

    int Npsi = 360 / basePckr->psi_sampling;

    int min_distance_pix = round(basePckr->min_particle_distance / basePckr->angpix);
    XFLOAT scale = (XFLOAT) basePckr->workSize / (XFLOAT) basePckr->micrograph_size;

    // Read in the micrograph
    {
    ifdefTIMING(TicToc tt (basePckr->timer.tic, basePckr->TIMING_A6);)
    CTICTOC(timer, "readMicrograph", ({ Imic.read(fn_mic); }))
    CTICTOC(timer, "setXmippOrigin_0", ({ Imic().setXmippOrigin(); }))
    }

    // Let's just check the square size again....
    RFLOAT my_xsize = XSIZE(Imic());
    RFLOAT my_ysize = YSIZE(Imic());
    RFLOAT my_size = std::max(my_xsize, my_ysize);

    if (
        my_xsize != basePckr->micrograph_xsize ||
        my_ysize != basePckr->micrograph_ysize ||
        my_size  != basePckr->micrograph_size
    ) {
        Imic().printShape();
        std::cerr << " micrograph_size= " << basePckr->micrograph_size << " micrograph_xsize= " << basePckr->micrograph_xsize << " micrograph_ysize= " << basePckr->micrograph_ysize << std::endl;
        REPORT_ERROR("AutoPicker::autoPickOneMicrograph ERROR: No differently sized micrographs are allowed in one run, sorry you will have to run separately for each size...");
    }

    if (!basePckr->do_read_fom_maps) {
        CTICTOC(timer, "setSize_micTr", ({
        micTransformer.setSize(basePckr->micrograph_size, basePckr->micrograph_size, 1, 1);
        }))

        CTICTOC(timer, "setSize_cudaTr", ({
        cudaTransformer1.setSize(basePckr->workSize,basePckr->workSize, 1, Npsi, FFTW_BACKWARD);
        }))
    }
    HANDLE_ERROR(cudaDeviceSynchronize());

    if (cudaTransformer1.batchSize.size() > 1 && !have_warned_batching) {
        have_warned_batching = true;
        std::cerr << std::endl << "*-----------------------------WARNING------------------------------------------------*" << std::endl;
        std::cerr              << "With the current settings the GPU memory is imposing a soft limit on your performace," << std::endl;
        std::cerr              << "since one or more micrographs has to use (at least " << cudaTransformer1.batchSize.size() << ") batches of orientations to " << std::endl;
        std::cerr              << "achieve the total requested " << Npsi << " orientations. Consider using" << std::endl;
        std::cerr              << "\t higher --ang" << std::endl;
        std::cerr              << "\t harder --shrink" << std::endl;
        std::cerr              << "\t higher --lowpass with --shrink 0" << std::endl;
        std::cerr              << "*------------------------------------------------------------------------------------*" << std::endl;
    }

    {
    ifdefTIMING(TicToc tt (basePckr->timer, basePckr->TIMING_A7);)
    const auto stats = [&] () -> Stats<RFLOAT> {
        CTICTOC(timer, "computeStats", ({
        // Set mean to zero and stddev to 1 to prevent numerical problems with one-sweep stddev calculations....
        return computeStats(Imic());
        }))}();
        avg0 = stats.avg;
        stddev0 = stats.stddev;
    }

    CTICTOC(timer, "middlePassFilter", ({
    for (long int n = 0; n < Imic().size(); n++) {
        // Remove pixel values that are too far away from the mean
        if (abs(Z(Imic()[n], avg0, stddev0)) > basePckr->outlier_removal_zscore)
            Imic()[n] = avg0;

        Imic()[n] = Z(Imic()[n], avg0, stddev0);
    }
    }))

    if (basePckr->micrograph_xsize != basePckr->micrograph_ysize) {
        // Window non-square micrographs to be a square with the largest side
        CTICTOC(timer, "rewindow", ({ rewindow(Imic, basePckr->micrograph_size); }))

        // Fill region outside the original window with white Gaussian noise to prevent all-zeros in Mstddev
        CTICTOC(timer, "gaussNoiseOutside", ({
        FOR_ALL_ELEMENTS_IN_ARRAY2D(Imic(), i, j) {
            if (
                j < Xmipp::init(basePckr->micrograph_ysize) ||
                j > Xmipp::last(basePckr->micrograph_ysize) ||
                i < Xmipp::init(basePckr->micrograph_xsize) ||
                i > Xmipp::last(basePckr->micrograph_xsize)
            ) {
                Imic().elem(i, j) = rnd_gaus(0.0, 1.0);
            }
        }
        }))
    }

    #ifdef TIMING
    basePckr->timer.tic(basePckr->TIMING_A8);
    #endif
    CTICTOC(timer, "CTFread", ({
    // Read in the CTF information if needed
    if (basePckr->do_ctf) {
        // Search for this micrograph in the metadata table
        for (long int index : basePckr->MDmic) {
            FileName fn_tmp = basePckr->MDmic.getValue(EMDL::MICROGRAPH_NAME);
            if (fn_tmp == fn_mic) {
                CTF ctf = CTF(basePckr->MDmic);
                Fctf = CtfHelper::getFftwImage(
                    ctf,
                    basePckr->workSize / 2 + 1, basePckr->workSize, 
                    basePckr->micrograph_size, basePckr->micrograph_size,
                    basePckr->angpix,
                    NULL,  // No ObservationModel
                    false, false, basePckr->intact_ctf_first_peak, true
                );
                break;
            }
        }
    }
    }))
    #ifdef TIMING
    basePckr->timer.toc(basePckr->TIMING_A8);
    #endif

    #ifdef TIMING
    basePckr->timer.tic(basePckr->TIMING_A9);
    #endif
    CTICTOC(timer, "mccfResize", ({
    Mccf_best.resize(basePckr->workSize,basePckr->workSize);
    }))

    CTICTOC(timer, "mpsiResize", ({
    Mpsi_best.resize(basePckr->workSize,basePckr->workSize);
    }))
    #ifdef TIMING
    basePckr->timer.toc(basePckr->TIMING_A9);
    #endif

    CudaGlobalPtr<CUDACOMPLEX> d_Fmic(allocator);
    CudaGlobalPtr<XFLOAT> d_Mmean(allocator);
    CudaGlobalPtr<XFLOAT> d_Mstddev(allocator);

    #ifdef TIMING
    basePckr->timer.tic(basePckr->TIMING_B1);
    #endif
    RFLOAT normfft = (RFLOAT) (basePckr->micrograph_size * basePckr->micrograph_size) / (RFLOAT) basePckr->nr_pixels_circular_mask;
    if (basePckr->do_read_fom_maps) {
        CTICTOC(timer, "readFromFomMaps_0", ({
        FileName fn_tmp = basePckr->getOutputRootName(fn_mic) + "_" + basePckr->fn_out + "_stddevNoise.spi";
        Image<RFLOAT> It;
        It.read(fn_tmp);
        Mstddev = It();
        }))
    } else {
        /*
         * Squared difference FOM:
         * Sum ( (X-mu)/sig  - A )^2 =
         *  = Sum((X-mu)/sig)^2 - 2 Sum (A*(X-mu)/sig) + Sum(A)^2
         *  = (1/sig^2)*Sum(X^2) - (2*mu/sig^2)*Sum(X) + (mu^2/sig^2)*Sum(1) - (2/sig)*Sum(AX) + (2*mu/sig)*Sum(A) + Sum(A^2)
         *
         * However, the squared difference with an "empty" ie all-zero reference is:
         * Sum ( (X-mu)/sig)^2
         *
         * The ratio of the probabilities thereby becomes:
         * P(ref) = 1/sqrt(2pi) * exp (( (X-mu)/sig  - A )^2 / -2 )   // assuming sigma = 1!
         * P(zero) = 1/sqrt(2pi) * exp (( (X-mu)/sig )^2 / -2 )
         *
         * P(ref)/P(zero) = exp(( (X-mu)/sig  - A )^2 / -2) / exp ( ( (X-mu)/sig )^2 / -2)
         *                = exp( (- (2/sig)*Sum(AX) + (2*mu/sig)*Sum(A) + Sum(A^2)) / - 2 )
         *
         *                Therefore, I do not need to calculate (X-mu)/sig beforehand!!!
         *
         */

        CTICTOC(timer, "Imic_insert", ({
        for (int i = 0; i < Imic().size(); i++)
            micTransformer.reals[i] = (XFLOAT) Imic().data[i];
        micTransformer.reals.cp_to_device();
        }))

        CTICTOC(timer, "runCenterFFT_0", ({
        runCenterFFT(micTransformer.reals, micTransformer.xSize, micTransformer.ySize, true, 1);
        }))

        CTICTOC(timer, "FourierTransform_0", ({
        micTransformer.forward();
        int FMultiBsize = ceilf((float) micTransformer.fouriers.getSize() * 2 / (float) BLOCK_SIZE);
        cuda_kernel_multi<<<FMultiBsize, BLOCK_SIZE>>>(
            (XFLOAT*) ~micTransformer.fouriers,
            (XFLOAT) 1 / (XFLOAT) micTransformer.reals.getSize(),
            micTransformer.fouriers.getSize() * 2
        );
        LAUNCH_HANDLE_ERROR(cudaGetLastError());
        }))

        if (basePckr->highpass > 0.0) {
            CTICTOC(timer, "highpass", ({
            micTransformer.fouriers.streamSync();
            lowPassFilterMapGPU(
                micTransformer.fouriers,
                (size_t) 1,
                micTransformer.yFSize, micTransformer.xFSize,
                XSIZE(Imic()),
                basePckr->lowpass, basePckr->highpass,
                basePckr->angpix,
                2,
                true // false = lowpass, true = highpass
            );
            micTransformer.fouriers.streamSync();
            micTransformer.backward();
            micTransformer.reals.streamSync();
            }))
        }

        CTICTOC(timer, "F_cp", ({
        CudaGlobalPtr<CUDACOMPLEX> Ftmp(allocator);
        Ftmp.setSize(micTransformer.fouriers.getSize());
        Ftmp.device_alloc();
        micTransformer.fouriers.cp_on_device(Ftmp);
        }))

        // Also calculate the FFT of the squared micrograph
        CTICTOC(timer, "SquareImic", ({

        cuda_kernel_square<<<FMultiBsize, BLOCK_SIZE>>>(
            ~micTransformer.reals, micTransformer.reals.getSize()
        );
        LAUNCH_HANDLE_ERROR(cudaGetLastError());
        }))

        CTICTOC(timer, "FourierTransform_1", ({

        micTransformer.forward();
        cuda_kernel_multi<<<FMultiBsize, BLOCK_SIZE>>>(
            (XFLOAT*) ~micTransformer.fouriers,
            (XFLOAT) 1 / (XFLOAT) micTransformer.reals.getSize(),
            micTransformer.fouriers.getSize() * 2
        );
        LAUNCH_HANDLE_ERROR(cudaGetLastError());
        }))

        // The following calculate mu and sig under the solvent area at every position in the micrograph
        CTICTOC(timer, "calculateStddevAndMeanUnderMask", ({

        d_Mstddev.device_alloc(basePckr->workSize * basePckr->workSize);
        d_Mmean.device_alloc(basePckr->workSize * basePckr->workSize);


        /// TODO: Do this only once further up in scope
        CudaGlobalPtr<CUDACOMPLEX> d_Fmsk(basePckr->Finvmsk.size(), allocator);
        for (int i = 0; i< d_Fmsk.size ; i++) {
            d_Fmsk[i].x = basePckr->Finvmsk.data[i].real;
            d_Fmsk[i].y = basePckr->Finvmsk.data[i].imag;
        }
        d_Fmsk.put_on_device();
        d_Fmsk.streamSync();

        calculateStddevAndMeanUnderMask(
            Ftmp, micTransformer.fouriers,
            d_Fmsk, basePckr->nr_pixels_circular_invmask,
            d_Mstddev, d_Mmean,
            micTransformer.xFSize, micTransformer.yFSize,
            basePckr->micrograph_size, basePckr->workSize
        );


        /// TODO: remove this
        d_Mstddev.host_alloc();
        d_Mstddev.cp_to_host();
        d_Mstddev.streamSync();

        Mstddev.resizeNoCp(basePckr->workSize, basePckr->workSize);

        /// TODO: put this in a kernel
        for (int i = 0; i < d_Mstddev.size ; i ++) {
            Mstddev.data[i] = d_Mstddev[i];
            if (d_Mstddev[i] > (XFLOAT) 1E-10) {
                d_Mstddev[i] = 1 / d_Mstddev[i];
            } else {
                d_Mstddev[i] = 1;
            }
        }

        d_Mstddev.cp_to_device();
        d_Mstddev.streamSync();

        }))

        // From now on use downsized Fmic, as the cross-correlation with the references can be done at lower resolution
        CTICTOC(timer, "windowFourierTransform_0", ({

        d_Fmic.setSize((basePckr->workSize / 2 + 1) * basePckr->workSize);
        d_Fmic.device_alloc();
        windowFourierTransform2(
            Ftmp, d_Fmic,
            basePckr->micrograph_size / 2 + 1, basePckr->micrograph_size, 1, // Input dimensions
            basePckr->workSize        / 2 + 1, basePckr->workSize,        1  // Output dimensions
        );
        }))

        if (basePckr->do_write_fom_maps) {
            CTICTOC(timer, "writeToFomMaps", ({
            // TMP output
            FileName fn_tmp = basePckr->getOutputRootName(fn_mic) + "_" + basePckr->fn_out + "_stddevNoise.spi";
            Image<RFLOAT> It;
            It() = Mstddev;
            It.write(fn_tmp);
            }))
        }
    }

    // Now start looking for the peaks of all references
    // Clear the output vector with all peaks
    std::vector<Peak> peaks;
    CTICTOC(timer, "initPeaks", ({ peaks.clear(); }))
    #ifdef TIMING
    basePckr->timer.toc(basePckr->TIMING_B1);
    #endif

    if (basePckr->autopick_helical_segments) {
        if (basePckr->do_read_fom_maps) {
            FileName fn_ccf = basePckr->getOutputRootName(fn_mic) + "_" + basePckr->fn_out + "_combinedCCF.spi";
            Image<RFLOAT> It_float;
            It_float.read(fn_ccf);
            Mccf_best_combined = It_float();

            FileName fn_class = basePckr->getOutputRootName(fn_mic) + "_" + basePckr->fn_out + "_combinedCLASS.spi";
            It_int.read(fn_class);
            Image<int> It_int;
            Mclass_best_combined = It_int();
        } else {
            Mccf_best_combined.clear();
            Mccf_best_combined.resize(basePckr->workSize, basePckr->workSize);
            Mccf_best_combined.initConstant(-99.0e99);

            Mclass_best_combined.clear();
            Mclass_best_combined.resize(basePckr->workSize, basePckr->workSize);
            Mclass_best_combined.initConstant(-1);
        }
    }

    CudaGlobalPtr<XFLOAT> d_ctf(Fctf.size(), allocator);
    if (basePckr->do_ctf) {
        for (int i = 0; i < d_ctf.size; i++)
            d_ctf[i] = Fctf.data[i];
        d_ctf.put_on_device();
    }

    for (int iref = 0; iref < basePckr->Mrefs.size(); iref++) {

        CTICTOC(timer, "OneReference", ({
        RFLOAT expected_Pratio; // the expectedFOM for this (ctf-corrected) reference
        if (basePckr->do_read_fom_maps) {
            #ifdef TIMING
            basePckr->timer.tic(basePckr->TIMING_B2);
            #endif
            if (!basePckr->autopick_helical_segments) {
                CTICTOC(timer, "readFromFomMaps", ({
                FileName fn_tmp;
                Image<RFLOAT> It;

                fn_tmp.compose(basePckr->getOutputRootName(fn_mic) + "_" + basePckr->fn_out + "_ref", iref, "_bestCCF.spi");
                It.read(fn_tmp);
                Mccf_best = It();
                expected_Pratio = It.MDMainHeader.getValue(EMDL::IMAGE_STATS_MAX);  // Retrieve expected_Pratio from the header of the image

                fn_tmp.compose(basePckr->getOutputRootName(fn_mic) + "_" + basePckr->fn_out + "_ref", iref, "_bestPSI.spi");
                It.read(fn_tmp);
                Mpsi_best = It();
                }))
            }
            #ifdef TIMING
            basePckr->timer.toc(basePckr->TIMING_B2);
            #endif
        } else {
            #ifdef TIMING
            basePckr->timer.tic(basePckr->TIMING_B3);
            #endif
            CTICTOC(timer, "mccfInit", ({
            deviceInitValue(d_Mccf_best, (XFLOAT) -LARGE_NUMBER);
            }))
            CudaProjectorKernel projKernel = CudaProjectorKernel::makeKernel(
                cudaProjectors[iref],
                (int) basePckr->workSize / 2 + 1,
                (int) basePckr->workSize,
                1, // Zdim, always 1 in autopicker.
                (int) basePckr->workSize / 2 + 1 -1 // ?!?
            );

            int FauxStride = (basePckr->workSize / 2 + 1) * basePckr->workSize;

            #ifdef TIMING
            basePckr->timer.tic(basePckr->TIMING_B4);
            #endif
            CTICTOC(timer, "SingleProjection", ({
            dim3 blocks((int) ceilf((float) FauxStride / (float) BLOCK_SIZE), 1);
            if (basePckr->do_ctf) {
                cuda_kernel_rotateAndCtf<<<blocks, BLOCK_SIZE>>>(
                    ~cudaTransformer1.fouriers, ~d_ctf, 0, projKernel, 0
                );
            } else {
                cuda_kernel_rotateOnly<<<blocks, BLOCK_SIZE>>>(
                    ~cudaTransformer1.fouriers, 0, projKernel, 0
                );
            }
            LAUNCH_HANDLE_ERROR(cudaGetLastError());
            }))
            #ifdef TIMING
            basePckr->timer.toc(basePckr->TIMING_B4);
            #endif
            /*
            *    FIRST PSI WAS USED FOR PREP CALCS - THIS IS NOW A DEDICATED SECTION
            *    -------------------------------------------------------------------
            */

            CTICTOC(timer, "PREP_CALCS", ({

            #ifdef TIMING
            basePckr->timer.tic(basePckr->TIMING_B5);
            #endif
            // Sjors 20 April 2016: The calculation for sum_ref_under_circ_mask, etc below needs to be done on original micrograph_size!
            CTICTOC(timer, "windowFourierTransform_FP", ({
            windowFourierTransform2(
                cudaTransformer1.fouriers,
                micTransformer.fouriers,
                basePckr->workSize / 2 + 1,        basePckr->workSize,        1, // Input dimensions
                basePckr->micrograph_size / 2 + 1, basePckr->micrograph_size, 1  // Output dimensions
            );
            }))

            CTICTOC(timer, "inverseFourierTransform_FP", ({
            micTransformer.backward();
            }))

            CTICTOC(timer, "runCenterFFT_FP", ({
            runCenterFFT(
                micTransformer.reals,
                (int) micTransformer.xSize, (int) micTransformer.ySize,
                false, 1
            );
            }))

            micTransformer.reals.cp_to_host();

            Maux.resizeNoCp(basePckr->micrograph_size, basePckr->micrograph_size);

            micTransformer.reals.streamSync();
            for (int i = 0; i < micTransformer.reals.size ; i ++)
                Maux.data[i] = micTransformer.reals[i];

            CTICTOC(timer, "setXmippOrigin_FP_0", ({ Maux.setXmippOrigin(); })
            /// TODO: check whether I need CenterFFT(Maux, false)
            // Sjors 20 Apr 2016: checked, somehow not needed.

            sum_ref_under_circ_mask = 0.0;
            sum_ref2_under_circ_mask = 0.0;
            RFLOAT suma2 = 0.0;
            RFLOAT sumn = 1.0;

            MultidimArray<RFLOAT> Mctfref(basePckr->particle_size, basePckr->particle_size);
            CTICTOC(timer, "setXmippOrigin_FP_1", ({
            Mctfref.setXmippOrigin();
            }))

            CTICTOC(timer, "suma_FP", ({
            FOR_ALL_ELEMENTS_IN_ARRAY2D(Mctfref, i, j) {
                // only loop over smaller Mctfref, but take values from large Maux!
                if (i * i + j * j < basePckr->particle_radius2) {
                    const auto &x = Maux.elem(i, j);
                    suma2 += x * x;
                    suma2 += 2.0 * x * rnd_gaus(0.0, 1.0);
                    sum_ref_under_circ_mask += x;
                    sum_ref2_under_circ_mask += x * x;
                    sumn += 1.0;
                }
            }
            sum_ref_under_circ_mask /= sumn;
            sum_ref2_under_circ_mask /= sumn;
            expected_Pratio = exp(suma2 / (2.0 * sumn));
            }))

            }))

            CTICTOC(timer, "AllPsi", ({
            int startPsi = 0;
            // for all batches
            for (int psiIter = 0; psiIter < cudaTransformer1.batchIters; psiIter++) {
                // psi-batches for possible memory-limits

                CTICTOC(timer, "Projection", ({
                dim3 blocks((int) ceilf((float) FauxStride / (float) BLOCK_SIZE), cudaTransformer1.batchSize[psiIter]);
                if (basePckr->do_ctf) {
                    cuda_kernel_rotateAndCtf<<<blocks, BLOCK_SIZE>>>(
                        ~cudaTransformer1.fouriers, ~d_ctf,
                        radians(basePckr->psi_sampling),
                        projKernel, startPsi
                    );
                } else {
                    cuda_kernel_rotateOnly<<<blocks, BLOCK_SIZE>>>(
                        ~cudaTransformer1.fouriers,
                        radians(basePckr->psi_sampling),
                        projKernel, startPsi
                    );
                }
                LAUNCH_HANDLE_ERROR(cudaGetLastError());
                }))

                // Now multiply template and micrograph to calculate the cross-correlation
                CTICTOC(timer, "convol", ({
                dim3 blocks2(
                    (int) ceilf((float) FauxStride / (float) BLOCK_SIZE),
                    cudaTransformer1.batchSize[psiIter]
                );
                cuda_kernel_batch_convol_A<<<blocks2, BLOCK_SIZE>>>(
                    cudaTransformer1.fouriers.d_ptr, d_Fmic.d_ptr, FauxStride
                );
                LAUNCH_HANDLE_ERROR(cudaGetLastError());
                }))

                CTICTOC(timer, "CudaInverseFourierTransform_1", ({
                cudaTransformer1.backward();
                HANDLE_ERROR(cudaDeviceSynchronize());
                }))


                CTICTOC(timer, "runCenterFFT_1", ({
                runCenterFFT(
                    cudaTransformer1.reals,
                    (int) cudaTransformer1.xSize,
                    (int) cudaTransformer1.ySize,
                    false,
                    cudaTransformer1.batchSize[psiIter]
                );
                }))
                // Calculate ratio of prabilities P(ref)/P(zero)
                // Keep track of the best values and their corresponding iref and psi
                // ------------------------------------------------------------------
                // So now we already had precalculated: Mdiff2 = 1/sig*Sum(X^2) - 2/sig*Sum(X) + mu^2/sig*Sum(1)
                // Still to do (per reference): - 2/sig*Sum(AX) + 2*mu/sig*Sum(A) + Sum(A^2)
                CTICTOC(timer, "probRatio", ({
                HANDLE_ERROR(cudaDeviceSynchronize());
                dim3 PR_blocks(ceilf((float)(cudaTransformer1.reals.size/cudaTransformer1.batchSize[psiIter])/(float)PROBRATIO_BLOCK_SIZE));
                cuda_kernel_probRatio<<<PR_blocks, PROBRATIO_BLOCK_SIZE>>>(
                    d_Mccf_best.d_ptr, d_Mpsi_best.d_ptr,
                    cudaTransformer1.reals.d_ptr,
                    d_Mmean.d_ptr, d_Mstddev.d_ptr,
                    cudaTransformer1.reals.size / cudaTransformer1.batchSize[0],
                    (XFLOAT) -2 * normfft,
                    (XFLOAT) 2 * sum_ref_under_circ_mask,
                    (XFLOAT) sum_ref2_under_circ_mask,
                    (XFLOAT) expected_Pratio,
                    cudaTransformer1.batchSize[psiIter],
                    startPsi, Npsi
                );
                LAUNCH_HANDLE_ERROR(cudaGetLastError());
                startPsi += cudaTransformer1.batchSize[psiIter];
                }))

            } // end for psi-batches
            }))
            #ifdef TIMING
            basePckr->timer.toc(basePckr->TIMING_B6);
            #endif

            #ifdef TIMING
            basePckr->timer.tic(basePckr->TIMING_B7);
            #endif
            CTICTOC(timer, "output", ({
            d_Mccf_best.cp_to_host();
            d_Mpsi_best.cp_to_host();
            d_Mccf_best.streamSync();
            for (int i = 0; i < Mccf_best.size(); i++) {
                Mccf_best.data[i] = d_Mccf_best[i];
                Mpsi_best.data[i] = d_Mpsi_best[i];
            }
            }))

            if (basePckr->do_write_fom_maps && !basePckr->autopick_helical_segments) {
                CTICTOC(timer, "writeFomMaps", ({
                // TMP output
                FileName fn_tmp;
                Image<RFLOAT> It;
                It() = Mccf_best;
                // Store expected_Pratio in the header of the image..
                It.MDMainHeader.setValue(EMDL::IMAGE_STATS_MAX, expected_Pratio);  // Store expected_Pratio in the header of the image
                fn_tmp.compose(basePckr->getOutputRootName(fn_mic) + "_" + basePckr->fn_out + "_ref", iref, "_bestCCF.spi");
                It.write(fn_tmp);

                It() = Mpsi_best;
                fn_tmp.compose(basePckr->getOutputRootName(fn_mic) + "_" + basePckr->fn_out + "_ref", iref, "_bestPSI.spi");
                It.write(fn_tmp);
                }))

            }
            #ifdef TIMING
            basePckr->timer.toc(basePckr->TIMING_B7);
            #endif
            #ifdef TIMING
            basePckr->timer.toc(basePckr->TIMING_B3);
            #endif
        }


        /// TODO: FIX HELICAL SEGMENTS SUPPORT
        if (basePckr->autopick_helical_segments) {
            if (!basePckr->do_read_fom_maps) {
                // Combine Mccf_best and Mpsi_best from all refs
                for (long int n = 0; n < Mccf_best.size(); n++) {
                    RFLOAT new_ccf = Mccf_best[n];
                    RFLOAT old_ccf = Mccf_best_combined[n];
                    if (new_ccf > old_ccf) {
                        Mccf_best_combined[n] = new_ccf;
                        Mclass_best_combined[n] = iref;
                    }
                }
            }
        } else {
            #ifdef TIMING
            basePckr->timer.tic(basePckr->TIMING_B8);
            #endif
            // Now that we have Mccf_best and Mpsi_best, get the peaks
            std::vector<Peak> my_ref_peaks;
            CTICTOC(timer, "setXmippOriginX3", ({
            Mstddev.setXmippOrigin();
            Mccf_best.setXmippOrigin();
            Mpsi_best.setXmippOrigin();
            }))

            CTICTOC(timer, "peakSearch", ({
            basePckr->peakSearch(Mccf_best, Mpsi_best, Mstddev, iref, my_skip_side, my_ref_peaks, scale);
            }))

            CTICTOC(timer, "peakPrune", ({
            basePckr->prunePeakClusters(my_ref_peaks, min_distance_pix, scale);
            }))

            CTICTOC(timer, "peakInsert", ({
            // append the peaks of this reference to all the other peaks
            peaks.insert(peaks.end(), my_ref_peaks.begin(), my_ref_peaks.end());
            }))

            #ifdef TIMING
            basePckr->timer.toc(basePckr->TIMING_B8);
            #endif
        }
        }))
    }

    if (basePckr->autopick_helical_segments) {
        RFLOAT thres = basePckr->min_fraction_expected_Pratio;
        int peak_r_min = 1;
        std::vector<ccfPeak> ccf_peak_list;
        std::vector<std::vector<ccfPeak>> tube_coord_list, tube_track_list;
        std::vector<RFLOAT> tube_len_list;
        MultidimArray<RFLOAT> Mccfplot;

        Mccf_best_combined.setXmippOrigin();
        Mclass_best_combined.setXmippOrigin();
        basePckr->pickCCFPeaks(
            Mccf_best_combined, Mclass_best_combined, thres, peak_r_min,
            (basePckr->particle_diameter / basePckr->angpix),
            ccf_peak_list, Mccfplot, my_skip_side, scale
        );
        basePckr->extractHelicalTubes(
            ccf_peak_list, tube_coord_list, tube_len_list, tube_track_list,
            basePckr->particle_diameter / basePckr->angpix, basePckr->helical_tube_curvature_factor_max,
            basePckr->min_particle_distance / basePckr->angpix,
            basePckr->helical_tube_diameter / basePckr->angpix, scale
        );
        basePckr->exportHelicalTubes(
            Mccf_best_combined, Mccfplot, Mclass_best_combined,
            tube_coord_list, tube_track_list, tube_len_list,
            fn_mic, basePckr->fn_out,
            basePckr->particle_diameter / basePckr->angpix,
            basePckr->helical_tube_length_min / basePckr->angpix,
            my_skip_side, scale
        );

        if (basePckr->do_write_fom_maps) {

            Image<RFLOAT> It_float;
            It_float() = Mccf_best_combined;
            It_float.write(basePckr->getOutputRootName(fn_mic) + "_" + basePckr->fn_out + "_combinedCCF.spi");

            Image<int> It_int;
            It_int() = Mclass_best_combined;
            It_int.write(basePckr->getOutputRootName(fn_mic) + "_" + basePckr->fn_out + "_combinedCLASS.spi");

        }

        if (basePckr->do_write_fom_maps || basePckr->do_read_fom_maps) {
            Image<RFLOAT> It;
            It() = Mccfplot;
            It.write(basePckr->getOutputRootName(fn_mic) + "_" + basePckr->fn_out + "_combinedPLOT.spi");
        }
    } else {
        #ifdef TIMING
        basePckr->timer.tic(basePckr->TIMING_B9);
        #endif
        // Now that we have done all references, prune the list again...
        CTICTOC(timer, "finalPeakPrune", ({
        basePckr->prunePeakClusters(peaks, min_distance_pix, scale);
        }))

        // And remove all too close neighbours
        basePckr->removeTooCloselyNeighbouringPeaks(peaks, min_distance_pix, scale);

        // Write out a STAR file with the coordinates
        MetaDataTable MDout;
        for (int ipeak = 0; ipeak < peaks.size(); ipeak++) {
            MDout.addObject();
            MDout.setValue(EMDL::IMAGE_COORD_X, (RFLOAT) peaks[ipeak].x / scale);
            MDout.setValue(EMDL::IMAGE_COORD_Y, (RFLOAT) peaks[ipeak].y / scale);
            MDout.setValue(EMDL::PARTICLE_CLASS,         peaks[ipeak].ref + 1); // start counting at 1
            MDout.setValue(EMDL::PARTICLE_AUTOPICK_FOM,  peaks[ipeak].fom);
            MDout.setValue(EMDL::ORIENT_PSI,             peaks[ipeak].psi);
        }
        FileName fn_tmp = basePckr->getOutputRootName(fn_mic) + "_" + basePckr->fn_out + ".star";
        MDout.write(fn_tmp);
        #ifdef TIMING
        basePckr->timer.toc(basePckr->TIMING_B9);
        #endif
    }
}
