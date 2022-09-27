/***************************************************************************
 *
 * Author: "Sjors H.W. Scheres"
 * MRC Laboratory of Molecular Biology
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * This complete copyright notice must be included in any revised version of the
 * source code. Additional authorship citations may be added, but existing
 * author citations must be preserved.
 ***************************************************************************/

// #define DEBUG_CHECKSIZES
// #define DEBUG_HELICAL_ORIENTATIONAL_SEARCH
// #define PRINT_GPU_MEM_INFO
// #define DEBUG_BODIES

#include <sys/time.h>
#include <stdio.h>
#include <time.h>
#include <math.h>
#include <ctime>
#include <iostream>
#include <string>
#include <fstream>
#include "src/macros.h"
#include "src/error.h"
#include "src/ml_optimiser.h"
#include "src/jaz/ctf_helper.h"
#ifdef CUDA
    #include "src/acc/cuda/cuda_ml_optimiser.h"
    #include <nvToolsExt.h>
    #include <cuda_profiler_api.h>
#endif
#ifdef ALTCPU
    #include <atomic>
    #include <tbb/tbb.h>
    #include <tbb/parallel_for.h>
    #define TBB_PREVIEW_GLOBAL_CONTROL 1
    #include <tbb/global_control.h>
    #include "src/acc/cpu/cpu_ml_optimiser.h"
#endif

#define NR_CLASS_MUTEXES 5

// Some global thread management variables
static pthread_mutex_t global_mutex2[NR_CLASS_MUTEXES] = { PTHREAD_MUTEX_INITIALIZER };
static pthread_mutex_t global_mutex = PTHREAD_MUTEX_INITIALIZER;
Barrier * global_barrier;
ThreadManager * global_ThreadManager;

template <typename T>
MultidimArray<T>& copyXYZ(MultidimArray<T> &dest, const MultidimArray<T> &src, long int n, long int offset = 0) {
    dest.resize(n, n, n);
    for (long int k = 0; k < n; k++)
    for (long int j = 0; j < n; j++)
    for (long int i = 0; i < n; i++)
    direct::elem(dest, i, j, k) = direct::elem(src, i, j, k + offset);
    return dest;
}

template <typename T>
MultidimArray<T>& copyXY(MultidimArray<T> &dest, const MultidimArray<T> &src, long int n, long int offset = 0) {
    dest.resize(n, n);
    for (long int j = 0; j < n; j++)
    for (long int i = 0; i < n; i++)
    direct::elem(dest, i, j) = direct::elem(src, i, j, offset);
    return dest;
}


template <typename T>
MultidimArray<T>& copyXYZ_fourier(MultidimArray<T> &FT, long int n, long int rank) {
    FT.resize(n / 2 + 1, n, rank == 3 ? n : 1);
    FT = -1;
    FOR_ALL_ELEMENTS_IN_FFTW_TRANSFORM(FT) {
        const int ires = round(euclid(ip, jp, kp));
        // Exclude points beyond ires, and exclude and half (y<0) of the x=0 column that is stored twice in FFTW
        // exclude lowest-resolution points
        // TODO: better check for volume_refine, but the same still seems to hold...
        // Half of the yz plane (either ip<0 or kp<0 is redundant at jp==0)
        if (ires < n / 2 + 1 && (jp != 0 || ip >= 0)) {
            direct::elem(FT, i, j, k) = ires;
        }
    }
    return FT;
}

Complex *get_shifted_image(
    int img_id, int exp_nr_oversampled_trans,
    int exp_nr_images, int itrans, int exp_itrans_min,  int iover_trans,
    std::vector<std::vector<MultidimArray<Complex>>> &exp_local_Fimgs_shifted,
    const std::vector<MultidimArray<Complex>> &global_fftshifts_ab_coarse,
    const std::vector<MultidimArray<Complex>> &global_fftshifts_ab_current,
    const std::vector<MultidimArray<Complex>> &global_fftshifts_ab2_coarse,
    const std::vector<MultidimArray<Complex>> &global_fftshifts_ab2_current,
    bool do_shifts_onthefly, bool do_skip_align,
    bool do_helical_refine, bool ignore_helical_symmetry,
    const std::vector<RFLOAT> &oversampled_translations_x,
    const std::vector<RFLOAT> &oversampled_translations_y,
    const std::vector<RFLOAT> &oversampled_translations_z,
    int data_dim, int ori_size,
    int exp_current_oversampling,
    const MultidimArray<RFLOAT> &exp_metadata,
    long int offset,
    float strict_highres_exp,
    int current_size, int coarse_size, long int ydim,
    MultidimArray<Complex> &Fimg_otfshift,
    TabSine &tab_sin, TabCosine &tab_cos
) {
    if (do_shifts_onthefly) {
        // Calculate shifted image on-the-fly to save replicating memory in multi-threaded jobs.
        // Feb01,2017 - Shaoda, on-the-fly shifts in helical reconstuctions (2D and 3D)
        auto &first_Fimg_shifted = exp_local_Fimgs_shifted[img_id][0];
        if (do_helical_refine && !ignore_helical_symmetry) {
            RFLOAT xshift = 0.0, yshift = 0.0, zshift = 0.0;
            xshift = oversampled_translations_x[exp_current_oversampling == 0 ? 0 : iover_trans];
            yshift = oversampled_translations_y[exp_current_oversampling == 0 ? 0 : iover_trans];
            if (data_dim == 3)
            zshift = oversampled_translations_z[exp_current_oversampling == 0 ? 0 : iover_trans];

            const RFLOAT rot_deg  = direct::elem(exp_metadata, offset, METADATA_ROT);
            const RFLOAT tilt_deg = direct::elem(exp_metadata, offset, METADATA_TILT);
            const RFLOAT psi_deg  = direct::elem(exp_metadata, offset, METADATA_PSI);
            transformCartesianAndHelicalCoords(
                xshift, yshift, zshift,
                rot_deg, tilt_deg, psi_deg,
                data_dim, HELICAL_TO_CART_COORDS
            );

            bool use_coarse_size =
                exp_current_oversampling == 0 && ydim == coarse_size ||
                exp_current_oversampling >  0 && strict_highres_exp > 0.0;
            shiftImageInFourierTransformWithTabSincos(
                first_Fimg_shifted,
                Fimg_otfshift,
                (RFLOAT) ori_size,
                use_coarse_size ? coarse_size : current_size,
                tab_sin, tab_cos,
                xshift, yshift, zshift
            );
        } else {
            Complex *AB;
            if (exp_current_oversampling == 0) {
                AB = (ydim == coarse_size ?
                    global_fftshifts_ab_coarse : global_fftshifts_ab_current
                )[itrans].data;
            } else {
                int iitrans = itrans * exp_nr_oversampled_trans + iover_trans;
                AB = (strict_highres_exp > 0.0 ?
                    global_fftshifts_ab2_coarse : global_fftshifts_ab2_current
                )[iitrans].data;
            }
            for (long int n = 0; n < first_Fimg_shifted.size(); n++) {
                Complex A = AB[n];
                Complex X = first_Fimg_shifted[n];
                Fimg_otfshift[n] = Complex(
                    A.real * X.real - A.imag * X.imag,  // A dot conj X
                    A.real * X.imag + A.imag * X.real   // A dot (i conj X)
                );
            }
        }
        return Fimg_otfshift.data;
    } else {
        long int ishift = img_id * exp_nr_oversampled_trans * exp_nr_images +
            (itrans - exp_itrans_min) * exp_nr_oversampled_trans + iover_trans;
        if (do_skip_align)
            ishift = img_id;
        #ifdef DEBUG_CHECKSIZES
        if (ishift >= exp_local_Fimgs_shifted.size()) {
            std::cerr<< "ishift= "<<ishift<<" exp_local_Fimgs_shifted.size()= "<< exp_local_Fimgs_shifted.size() <<std::endl;
            std::cerr << " itrans= " << itrans << std::endl;
            std::cerr << " img_id= " << img_id << std::endl;
            std::cerr << " exp_nr_oversampled_trans= " << exp_nr_oversampled_trans << " exp_nr_trans= " << exp_nr_trans << " iover_trans= " << iover_trans << std::endl;
            REPORT_ERROR("ishift >= exp_local_Fimgs_shifted.size()");
        }
        #endif
        return exp_local_Fimgs_shifted[img_id][ishift].data;
    }
}

template <typename T>
T tryget(const MetaDataTable &mdt, EMDL::EMDLabel emdl_index, T defaultval, long object_id = -1) {
    try {
        return mdt.getValue<T>(emdl_index, object_id);
    } catch (const char* errmsg) {
        return defaultval;
    }
}

/** ========================== Threaded parallelization of expectation ===== */

void globalThreadExpectationSomeParticles(ThreadArgument &thArg) {
    MlOptimiser *MLO = (MlOptimiser*) thArg.workClass;

    try {
        #ifdef CUDA
        if (MLO->do_gpu)
            ((MlOptimiserCuda*) MLO->cudaOptimisers[thArg.thread_id])->doThreadExpectationSomeParticles(thArg.thread_id);
        else
        #endif
        MLO->doThreadExpectationSomeParticles(thArg.thread_id);
    } catch (RelionError XE) {
        RelionError *gE = new RelionError(XE.msg, XE.file, XE.line);
        gE->msg = XE.msg;
        MLO->threadException = gE;
    }
}


/** ========================== I/O operations  =========================== */

void MlOptimiser::usage() {
    parser.writeUsage(std::cout);
}

void MlOptimiser::read(int argc, char **argv, int rank) {
    // #define DEBUG_READ

    parser.setCommandLine(argc, argv);

    if (checkParameter(argc, argv, "--continue")) {
        // Do this before reading in the data.star file below!
        do_preread_images   = checkParameter(argc, argv, "--preread_images");
        do_parallel_disc_io = !checkParameter(argc, argv, "--no_parallel_disc_io");

        parser.addSection("Continue options");
        FileName fn_in = parser.getOption("--continue", "_optimiser.star file of the iteration after which to continue");
        // Read in previously calculated parameters
        if (!fn_in.empty())
            read(fn_in, rank);

        // And look for additional command-line options...
        parseContinue(argc, argv);
    } else {
        // Start a new run from scratch
        parseInitial(argc, argv);
    }
}

void MlOptimiser::parseContinue(int argc, char **argv) {
    #ifdef DEBUG
    std::cerr << "Entering parseContinue" << std::endl;
    #endif

    int general_section = parser.addSection("General options");
    // Not all parameters are accessible here...
    FileName fn_out_new = parser.getOption("--o", "Output rootname", "OLD_ctX");
    if (fn_out_new == "OLD_ctX" || fn_out_new == fn_out ) {
        fn_out += "_ct" + integerToString(iter);
    } else {
        fn_out = fn_out_new;
    }

    do_force_converge =  parser.checkOption("--force_converge", "Force an auto-refinement run to converge immediately upon continuation.");

    // For multi-body refinement
    bool fn_body_masks_was_empty = fn_body_masks == "None";
    std::string fnt;
    fnt = parser.getOption("--multibody_masks", "STAR file with masks and metadata for multi-body refinement", "OLD");
    if (fnt != "OLD")
        fn_body_masks = fnt;
    // Don't use _ctXX at start of a multibody refinement
    if (fn_body_masks_was_empty && !fn_body_masks.empty())
        fn_out = parser.getOption("--o", "Output rootname", "run");

    // Also allow change of padding...
    fnt = parser.getOption("--pad", "Oversampling factor for the Fourier transforms of the references", "OLD");
    if (fnt != "OLD") {
        if (textToInteger(fnt) != mymodel.padding_factor) {
            if (mymodel.nr_bodies > 1)
                REPORT_ERROR("ERROR: cannot change padding factor in a continuation of a multi-body refinement...");
            mymodel.padding_factor = textToInteger(fnt);
            // Re-initialise the model to get the right padding factors in the PPref vectors
            mymodel.initialise();
        }
    }

    // Is this a new multi-body refinement?
    do_initialise_bodies = fn_body_masks_was_empty && fn_body_masks != "None";

    if (do_initialise_bodies) {
        ini_high = textToFloat(parser.getOption("--ini_high", "Resolution (in Angstroms) to which to limit refinement in the first iteration ", "-1"));
        mymodel.norm_body_mask_overlap = parser.checkOption("--multibody_norm_overlap", "Overlapping regions between bodies are normalized. This reduces memory requirements.");
    }
    do_reconstruct_subtracted_bodies = parser.checkOption("--reconstruct_subtracted_bodies", "Use this flag to perform reconstructions with the subtracted images in multi-body refinement");

    fnt = parser.getOption("--iter", "Maximum number of iterations to perform", "OLD");
    if (fnt != "OLD")
        nr_iter = textToInteger(fnt);

    fnt = parser.getOption("--tau2_fudge", "Regularisation parameter (values higher than 1 give more weight to the data)", "OLD");
    if (fnt != "OLD")
        mymodel.tau2_fudge_factor = textToFloat(fnt);

    auto_ignore_angle_changes = parser.checkOption("--auto_ignore_angles", "In auto-refinement, update angular sampling regardless of changes in orientations for convergence. This makes convergence faster.");
    auto_resolution_based_angles= parser.checkOption("--auto_resol_angles", "In auto-refinement, update angular sampling based on resolution-based required sampling. This makes convergence faster.");
    allow_coarser_samplings = parser.checkOption("--allow_coarser_sampling", "In 2D/3D classification, allow coarser angular and translational samplings if accuracies are bad (typically in earlier iterations.");

    // Solvent flattening
    do_solvent |= parser.checkOption("--flatten_solvent", "Switch on masking on the references?", "OLD");

    // Check whether the mask has changed
    fnt = parser.getOption("--solvent_mask", "User-provided mask for the references", "OLD");
    if (fnt != "OLD")
        fn_mask = fnt;

    // Check whether the secondary mask has changed
    fnt = parser.getOption("--solvent_mask2", "User-provided secondary mask", "OLD");
    if (fnt != "OLD")
        fn_mask2 = fnt;

    // These are still experimental; so not in the optimiser.star yet.
    fn_lowpass_mask = parser.getOption("--lowpass_mask", "User-provided mask for low-pass filtering", "None");
    lowpass = textToFloat(parser.getOption("--lowpass", "User-provided cutoff for region specified above", "0"));

    // Check whether tau2-spectrum has changed
    fnt = parser.getOption("--tau", "STAR file with input tau2-spectrum (to be kept constant)", "OLD");
    if (fnt != "OLD")
        fn_tau = fnt;

    // Check whether particle diameter has changed
    fnt = parser.getOption("--particle_diameter", "Diameter of the circular mask that will be applied to the experimental images (in Angstroms)", "OLD");
    if (fnt != "OLD")
        particle_diameter = textToFloat(fnt);

    // SGD stuff
    fnt = parser.getOption("--sgd_ini_iter", "Number of initial SGD iterations", "OLD");
    if (fnt != "OLD")
        sgd_ini_iter = textToInteger(fnt);

    fnt = parser.getOption("--sgd_fin_iter", "Number of final SGD iterations", "OLD");
    if (fnt != "OLD")
        sgd_fin_iter = textToInteger(fnt);

    fnt = parser.getOption("--sgd_inbetween_iter", "Number of SGD iterations between the initial and final ones", "OLD");
    if (fnt != "OLD")
        sgd_inbetween_iter = textToInteger(fnt);

    fnt = parser.getOption("--sgd_ini_resol", "Resolution cutoff during the initial SGD iterations (A)", "OLD");
    if (fnt != "OLD")
        sgd_ini_resol = textToFloat(fnt);

    fnt = parser.getOption("--sgd_fin_resol", "Resolution cutoff during the final SGD iterations (A)", "OLD");
    if (fnt != "OLD")
        sgd_fin_resol = textToFloat(fnt);

    fnt = parser.getOption("--sgd_ini_subset", "Mini-batch size during the initial SGD iterations", "OLD");
    if (fnt != "OLD")
        sgd_ini_subset_size = textToInteger(fnt);

    fnt = parser.getOption("--sgd_fin_subset", "Mini-batch size during the final SGD iterations", "OLD");
    if (fnt != "OLD")
        sgd_fin_subset_size = textToInteger(fnt);

    fnt = parser.getOption("--sgd_stepsize", "Step size parameter for SGD updates", "OLD");
    if (fnt != "OLD")
        sgd_stepsize = textToInteger(fnt);

    fnt = parser.getOption("--mu", "Momentum parameter for SGD updates", "OLD");
    if (fnt != "OLD")
        mu = textToFloat(fnt);

    fnt = parser.getOption("--sgd_write_iter", "Write out model every so many iterations in SGD", "OLD");
    if (fnt != "OLD")
        write_every_sgd_iter = textToInteger(fnt);

    fnt = parser.getOption("--relax_sym", "The symmetry to be relaxed", "OLD");
    if (fnt != "OLD") {
        sampling.fn_sym_relax = fnt;
    }

    do_join_random_halves = parser.checkOption("--join_random_halves", "Join previously split random halves again (typically to perform a final reconstruction).");

    // ORIENTATIONS
    int orientations_section = parser.addSection("Orientations");

    fnt = parser.getOption("--oversampling", "Adaptive oversampling order to speed-up calculations (0=no oversampling, 1=2x, 2=4x, etc)", "OLD");
    if (fnt != "OLD")
        adaptive_oversampling = textToInteger(fnt);

    // Check whether angular sampling has changed
    // Do not do this for auto_refine, but make sure to do this when initialising multi-body refinement!
    if (!do_auto_refine || do_initialise_bodies) {
        directions_have_changed = false;
        fnt = parser.getOption("--healpix_order", "Healpix order for the angular sampling rate on the sphere (before oversampling): hp2=15deg, hp3=7.5deg, etc", "OLD");
        if (fnt != "OLD") {
            int _order = textToInteger(fnt);
            if (_order != sampling.healpix_order) {
                directions_have_changed = true;
                sampling.healpix_order = _order;
            }
        }

        fnt = parser.getOption("--psi_step", "Angular sampling (before oversampling) for the in-plane angle (default=10deg for 2D, hp sampling for 3D)", "OLD");
        if (fnt != "OLD")
            sampling.psi_step = textToFloat(fnt);

        fnt = parser.getOption("--offset_range", "Search range for origin offsets (in pixels)", "OLD");
        if (fnt != "OLD") {
            sampling.offset_range = textToFloat(fnt);
            sampling.offset_range *= mymodel.pixel_size; // sampling.offset_range is in Angstroms, but command line in pixels!
        }

        fnt = parser.getOption("--offset_step", "Sampling rate for origin offsets (in pixels)", "OLD");
        if (fnt != "OLD") {
            sampling.offset_step = textToFloat(fnt);
            sampling.offset_step *= mymodel.pixel_size; // sampling.offset_step is in Angstroms, but command line in pixels!
        }
    }

    fnt = parser.getOption("--auto_local_healpix_order", "Minimum healpix order (before oversampling) from which auto-refine procedure will use local searches", "OLD");
    if (fnt != "OLD")
        autosampling_hporder_local_searches = textToInteger(fnt);

    // Check whether the prior mode changes
    RFLOAT _sigma_rot, _sigma_tilt, _sigma_psi, _sigma_off;
    int _mode;
    fnt = parser.getOption("--sigma_ang", "Stddev on all three Euler angles for local angular searches (of +/- 3 stddev)", "OLD");
    if (fnt != "OLD") {
        mymodel.orientational_prior_mode = PRIOR_ROTTILT_PSI;
        mymodel.sigma2_rot = mymodel.sigma2_tilt = mymodel.sigma2_psi = textToFloat(fnt) * textToFloat(fnt);
    }
    fnt = parser.getOption("--sigma_rot", "Stddev on the first Euler angle for local angular searches (of +/- 3 stddev)", "OLD");
    if (fnt != "OLD") {
        mymodel.orientational_prior_mode = PRIOR_ROTTILT_PSI;
        mymodel.sigma2_rot = textToFloat(fnt) * textToFloat(fnt);
    }
    fnt = parser.getOption("--sigma_tilt", "Stddev on the first Euler angle for local angular searches (of +/- 3 stddev)", "OLD");
    if (fnt != "OLD") {
        mymodel.orientational_prior_mode = PRIOR_ROTTILT_PSI;
        mymodel.sigma2_tilt = textToFloat(fnt) * textToFloat(fnt);
    }
    fnt = parser.getOption("--sigma_psi", "Stddev on the in-plane angle for local angular searches (of +/- 3 stddev)", "OLD");
    if (fnt != "OLD") {
        mymodel.orientational_prior_mode = PRIOR_ROTTILT_PSI;
        mymodel.sigma2_psi = textToFloat(fnt) * textToFloat(fnt);
    }
    fnt = parser.getOption("--sigma_off", "Stddev. on the translations", "OLD");
    if (fnt != "OLD") {
        mymodel.sigma2_offset = textToFloat(fnt) * textToFloat(fnt);
    }
    fnt = parser.getOption("--helical_inner_diameter", "Inner diameter of helical tubes in Angstroms (for masks of helical references and particles)", "OLD");
    if (fnt != "OLD") {
        helical_tube_inner_diameter = textToFloat(fnt);
    }
    fnt = parser.getOption("--helical_outer_diameter", "Outer diameter of helical tubes in Angstroms (for masks of helical references and particles)", "OLD");
    if (fnt != "OLD") {
        helical_tube_outer_diameter = textToFloat(fnt);
    }
    fnt = parser.getOption("--perturb", "Perturbation factor for the angular sampling (0=no perturb; 0.5=perturb)", "OLD");
    if (fnt != "OLD") {
        sampling.perturbation_factor = textToFloat(fnt);
    }

    // do_skip_align should normally be false
    do_skip_align = parser.checkOption("--skip_align", "Skip orientational assignment (only classify)?");

    // do_skip_rotate should normally be false...
    do_skip_rotate = parser.checkOption("--skip_rotate", "Skip rotational assignment (only translate and classify)?");

    do_bimodal_psi = parser.checkOption("--bimodal_psi", "Do bimodal searches of psi angle?"); // 7 Oct 2015 - Shaoda, bimodal psi

    do_skip_maximization = parser.checkOption("--skip_maximize", "Skip maximization step (only write out data.star file)?");

    int corrections_section = parser.addSection("Corrections");

    do_ctf_padding = parser.checkOption("--pad_ctf", "Perform CTF padding to treat CTF aliaising better?");
    if (do_ctf_padding)
        REPORT_ERROR("--pad_ctf currently disabled.");

    // Can also switch the following option OFF
    do_scale_correction |=  parser.checkOption("--scale",    "Switch on intensity-scale corrections on image groups",  "OLD");
    do_scale_correction &= !parser.checkOption("--no_scale", "Switch off intensity-scale corrections on image groups", "OLD");

    // Can also switch the following option OFF
    do_norm_correction |=  parser.checkOption("--norm",    "Switch on normalisation-error correction",  "OLD");
    do_norm_correction &= !parser.checkOption("--no_norm", "Switch off normalisation-error correction", "OLD");

    int computation_section = parser.addSection("Computation");

    x_pool = textToInteger(parser.getOption("--pool", "Number of images to pool for each thread task", "1"));
    nr_threads = textToInteger(parser.getOption("--j", "Number of threads to run in parallel (only useful on multi-core machines)", "1"));
    do_parallel_disc_io = !parser.checkOption("--no_parallel_disc_io", "Do NOT let parallel (MPI) processes access the disc simultaneously (use this option with NFS)");
    combine_weights_thru_disc = !parser.checkOption("--dont_combine_weights_via_disc", "Send the large arrays of summed weights through the MPI network, instead of writing large files to disc");
    do_shifts_onthefly = parser.checkOption("--onthefly_shifts", "Calculate shifted images on-the-fly, do not store precalculated ones in memory");
    do_preread_images  = parser.checkOption("--preread_images", "Use this to let the leader process read all particles into memory. Be careful you have enough RAM for large data sets!");
    fn_scratch = parser.getOption("--scratch_dir", "If provided, particle stacks will be copied to this local scratch disk prior to refinement.", "");
    keep_free_scratch_Gb = textToFloat(parser.getOption("--keep_free_scratch", "Space available for copying particle stacks (in Gb)", "10"));
    do_reuse_scratch = parser.checkOption("--reuse_scratch", "Re-use data on scratchdir, instead of wiping it and re-copying all data. This works only when ALL particles have already been cached.");
    keep_scratch = parser.checkOption("--keep_scratch", "Don't remove scratch after convergence. Following jobs that use EXACTLY the same particles should use --reuse_scratch.");

    #ifdef ALTCPU
    do_cpu = parser.checkOption("--cpu", "Use intel vectorisation implementation for CPU");
    #else
    do_cpu = false;
    #endif

    failsafe_threshold = textToInteger(parser.getOption("--failsafe_threshold", "Maximum number of particles permitted to be drop, due to zero sum of weights, before exiting with an error (GPU only).", "40"));

    do_gpu = parser.checkOption("--gpu", "Use available gpu resources for some calculations");
    gpu_ids = parser.getOption("--gpu", "Device ids for each MPI-thread","default");
    #ifndef CUDA
    if (do_gpu) {
        std::cerr << "+ WARNING : Relion was compiled without CUDA of at least version 7.0 - you do NOT have support for GPUs" << std::endl;
        do_gpu = false;
    }
    #endif
    double temp_reqSize = textToDouble(parser.getOption("--free_gpu_memory", "GPU device memory (in Mb) to leave free after allocation.", "0"));
    if (!do_zero_mask)
        temp_reqSize += 100;
    temp_reqSize *= 1000 * 1000;
    if (temp_reqSize < 0) {
        REPORT_ERROR("Invalid free_gpu_memory value.");
    } else {
        requested_free_gpu_memory = temp_reqSize;
    }

    // only allow switching ON solvent_fsc, not off
    do_phase_random_fsc |= parser.checkOption("--solvent_correct_fsc", "Correct FSC curve for the effects of the solvent mask?");
    verb = textToInteger(parser.getOption("--verb", "Verbosity (1=normal, 0=silent)", "1"));

    int expert_section = parser.addSection("Expert options");

    fnt = parser.getOption("--strict_highres_exp", "Resolution limit (in Angstrom) to restrict probability calculations in the expectation step", "OLD");
    if (fnt != "OLD")
        strict_highres_exp = textToFloat(fnt);

    do_trust_ref_size = parser.checkOption("--trust_ref_size", "Trust the pixel and box size of the input reference; by default the program will die if these are different from the first optics group of the data");

    // Debugging/analysis/hidden stuff
    do_map = !checkParameter(argc, argv, "--no_map");
    minres_map = textToInteger(getParameter(argc, argv, "--minres_map", "5"));
    gridding_nr_iter = textToInteger(getParameter(argc, argv, "--gridding_iter", "10"));
    debug1 = textToFloat(getParameter(argc, argv, "--debug1", "0."));
    debug2 = textToFloat(getParameter(argc, argv, "--debug2", "0."));
    debug3 = textToFloat(getParameter(argc, argv, "--debug3", "0."));
    do_bfactor = checkParameter(argc, argv, "--bfactor");
    // Read in initial sigmaNoise spectrum
    fn_sigma = getParameter(argc, argv, "--sigma","");
    sigma2_fudge = textToFloat(getParameter(argc, argv, "--sigma2_fudge", "1."));
    do_acc_currentsize_despite_highres_exp = checkParameter(argc, argv, "--accuracy_current_size");
    do_sequential_halves_recons  = checkParameter(argc, argv, "--sequential_halves_recons");
    do_always_join_random_halves = checkParameter(argc, argv, "--always_join_random_halves");
    do_use_all_data = checkParameter(argc, argv, "--use_all_data");
    do_always_cc  = checkParameter(argc, argv, "--always_cc");
    do_only_sample_tilt  = checkParameter(argc, argv, "--only_sample_tilt");
    minimum_angular_sampling = textToFloat(getParameter(argc, argv, "--minimum_angular_sampling", "0"));
    maximum_angular_sampling = textToFloat(getParameter(argc, argv, "--maximum_angular_sampling", "0"));
    asymmetric_padding = parser.checkOption("--asymmetric_padding", "", "false", true);
    maximum_significants = textToInteger(parser.getOption("--maxsig", "Maximum number of poses & translations to consider", "-1"));
    skip_gridding = parser.checkOption("--skip_gridding", "Skip gridding in the M step");
    nr_iter_max = textToInteger(parser.getOption("--auto_iter_max", "In auto-refinement, stop at this iteration.", "999"));
    debug_split_random_half = textToInteger(getParameter(argc, argv, "--debug_split_random_half", "0"));

    do_print_metadata_labels = false;
    do_print_symmetry_ops = false;
    #ifdef DEBUG
    std::cerr << "Leaving parseContinue" << std::endl;
    #endif

}

void MlOptimiser::parseInitial(int argc, char **argv) {
    #ifdef DEBUG_READ
    std::cerr<<"MlOptimiser::parseInitial Entering "<<std::endl;
    #endif

    // Read/initialise mymodel and sampling from a STAR file
    FileName fn_model = getParameter(argc, argv, "--model", "None");
    if (fn_model != "None") {
        mymodel.read(fn_model);
    }
    // Read in the sampling information from a _sampling.star file
    FileName fn_sampling = getParameter(argc, argv, "--sampling", "None");
    if (fn_sampling != "None") {
        sampling.read(fn_sampling);
    }

    // General optimiser I/O stuff
    int general_section = parser.addSection("General options");
    fn_data = parser.getOption("--i", "Input images (in a star-file)", "");
    fn_out = parser.getOption("--o", "Output rootname", "");
    nr_iter = textToInteger(parser.getOption("--iter", "Maximum number of iterations to perform", "50"));
    mymodel.tau2_fudge_factor = textToFloat(parser.getOption("--tau2_fudge", "Regularisation parameter (values higher than 1 give more weight to the data)", "1"));
    mymodel.nr_classes = textToInteger(parser.getOption("--K", "Number of references to be refined", "1"));
    particle_diameter = textToFloat(parser.getOption("--particle_diameter", "Diameter of the circular mask that will be applied to the experimental images (in Angstroms)", "-1"));
    do_zero_mask = parser.checkOption("--zero_mask","Mask surrounding background in particles to zero (by default the solvent area is filled with random noise)");
    do_solvent = parser.checkOption("--flatten_solvent", "Perform masking on the references as well?");
    fn_mask = parser.getOption("--solvent_mask", "User-provided mask for the references (default is to use spherical mask with particle_diameter)", "None");
    fn_mask2 = parser.getOption("--solvent_mask2", "User-provided secondary mask (with its own average density)", "None");
    fn_lowpass_mask = parser.getOption("--lowpass_mask", "User-provided mask for low-pass filtering", "None");
    lowpass = textToFloat(parser.getOption("--lowpass", "User-provided cutoff for region specified above", "0"));
    fn_tau = parser.getOption("--tau", "STAR file with input tau2-spectrum (to be kept constant)", "None");
    fn_local_symmetry = parser.getOption("--local_symmetry", "Local symmetry description file containing list of masks and their operators", "None");
    do_split_random_halves = parser.checkOption("--split_random_halves", "Refine two random halves of the data completely separately");
    low_resol_join_halves = textToFloat(parser.getOption("--low_resol_join_halves", "Resolution (in Angstrom) up to which the two random half-reconstructions will not be independent to prevent diverging orientations","-1"));

    // Initialisation
    int init_section = parser.addSection("Initialisation");
    fn_ref = parser.getOption("--ref", "Image, stack or star-file with the reference(s). (Compulsory for 3D refinement!)", "None");
    is_3d_model = parser.checkOption("--denovo_3dref", "Make an initial 3D model from randomly oriented 2D particles");
    mymodel.sigma2_offset = textToFloat(parser.getOption("--offset", "Initial estimated stddev for the origin offsets (in Angstroms)", "10"));
    mymodel.sigma2_offset *= mymodel.sigma2_offset;

    // Perform cross-product comparison at first iteration
    do_firstiter_cc = parser.checkOption("--firstiter_cc", "Perform CC-calculation in the first iteration (use this if references are not on the absolute intensity scale)");
    ini_high = textToFloat(parser.getOption("--ini_high", "Resolution (in Angstroms) to which to limit refinement in the first iteration ", "-1"));

    // Set the orientations
    int orientations_section = parser.addSection("Orientations");
    // Move these to sampling
    adaptive_oversampling = textToInteger(parser.getOption("--oversampling", "Adaptive oversampling order to speed-up calculations (0=no oversampling, 1=2x, 2=4x, etc)", "1"));
    sampling.healpix_order = textToInteger(parser.getOption("--healpix_order", "Healpix order for the angular sampling (before oversampling) on the (3D) sphere: hp2=15deg, hp3=7.5deg, etc", "2"));
    sampling.psi_step = textToFloat(parser.getOption("--psi_step", "Sampling rate (before oversampling) for the in-plane angle (default=10deg for 2D, hp sampling for 3D)", "-1"));
    sampling.limit_tilt = textToFloat(parser.getOption("--limit_tilt", "Limited tilt angle: positive for keeping side views, negative for keeping top views", "-91"));

    std::string sym_ = parser.getOption("--sym", "Symmetry group", "c1");

    // Check if a comma-separated list was provided
    if (sym_.find(",") != std::string::npos) {
        std::stringstream ss(sym_);
        std::string item;
        while (std::getline(ss, item, ','))
            fn_multi_sym.push_back(item);
    } else {
        sampling.fn_sym = sym_;
    }

    // Check for relax_symmetry option
    std::string sym_relax_ = parser.getOption("--relax_sym", "Symmetry to be relaxed", "");
    sampling.fn_sym_relax = sym_relax_;

    sampling.offset_range = textToFloat(parser.getOption("--offset_range", "Search range for origin offsets (in pixels)", "6"));
    sampling.offset_step = textToFloat(parser.getOption("--offset_step", "Sampling rate (before oversampling) for origin offsets (in pixels)", "2"));
    // 19 Jun 2015 - Shaoda, Helical refinement
    sampling.helical_offset_step = textToFloat(parser.getOption("--helical_offset_step", "Sampling rate (before oversampling) for offsets along helical axis (in Angstroms)", "-1"));
    sampling.perturbation_factor = textToFloat(parser.getOption("--perturb", "Perturbation factor for the angular sampling (0=no perturb; 0.5=perturb)", "0.5"));
    do_auto_refine = parser.checkOption("--auto_refine", "Perform 3D auto-refine procedure?");
    autosampling_hporder_local_searches = textToInteger(parser.getOption("--auto_local_healpix_order", "Minimum healpix order (before oversampling) from which autosampling procedure will use local searches", "4"));
    parser.setSection(orientations_section);
    RFLOAT _sigma_ang = textToFloat(parser.getOption("--sigma_ang", "Stddev on all three Euler angles for local angular searches (of +/- 3 stddev)", "-1"));
    RFLOAT _sigma_rot = textToFloat(parser.getOption("--sigma_rot", "Stddev on the first Euler angle for local angular searches (of +/- 3 stddev)", "-1"));
    RFLOAT _sigma_tilt = textToFloat(parser.getOption("--sigma_tilt", "Stddev on the second Euler angle for local angular searches (of +/- 3 stddev)", "-1"));
    RFLOAT _sigma_psi = textToFloat(parser.getOption("--sigma_psi", "Stddev on the in-plane angle for local angular searches (of +/- 3 stddev)", "-1"));

    if (_sigma_ang > 0.0) {
        mymodel.orientational_prior_mode = PRIOR_ROTTILT_PSI;
        // the sigma-values for the orientational prior are in model (and not in sampling) because one might like to estimate them
        // from the data by calculating weighted sums of all angular differences: therefore it needs to be in wsum_model and thus in mymodel.
        mymodel.sigma2_rot = mymodel.sigma2_tilt = mymodel.sigma2_psi = _sigma_ang * _sigma_ang;
    } else if (_sigma_rot > 0.0 || _sigma_tilt > 0.0 || _sigma_psi > 0.0) {
        mymodel.orientational_prior_mode = PRIOR_ROTTILT_PSI;
        mymodel.sigma2_rot  = _sigma_rot  > 0.0 ? _sigma_rot  * _sigma_rot  : 0.0;
        mymodel.sigma2_tilt = _sigma_tilt > 0.0 ? _sigma_tilt * _sigma_tilt : 0.0;
        mymodel.sigma2_psi  = _sigma_psi  > 0.0 ? _sigma_psi  * _sigma_psi  : 0.0;
    } else {
        // default
        // Very small to force the algorithm to take the current orientation
        if (!sym_relax_.empty()) {
            mymodel.orientational_prior_mode = PRIOR_ROTTILT_PSI;
            _sigma_ang = 0.0033;
            mymodel.sigma2_rot = mymodel.sigma2_tilt = mymodel.sigma2_psi = _sigma_ang * _sigma_ang;
        } else {
            mymodel.orientational_prior_mode = NOPRIOR;
            mymodel.sigma2_rot = mymodel.sigma2_tilt = mymodel.sigma2_psi = 0.0;
        }
    }
    do_skip_align = parser.checkOption("--skip_align", "Skip orientational assignment (only classify)?");
    do_skip_rotate = parser.checkOption("--skip_rotate", "Skip rotational assignment (only translate and classify)?");
    do_bimodal_psi = parser.checkOption("--bimodal_psi", "Do bimodal searches of psi angle?"); // Oct07,2015 - Shaoda, bimodal psi
    do_skip_maximization = false;

    // Helical reconstruction
    int helical_section = parser.addSection("Helical reconstruction (in development...0)");
    do_helical_refine = parser.checkOption("--helix", "Perform 3D classification or refinement for helices?");
    ignore_helical_symmetry = parser.checkOption("--ignore_helical_symmetry", "Ignore helical symmetry?");
    mymodel.helical_nr_asu = textToInteger(parser.getOption("--helical_nr_asu", "Number of new helical asymmetric units (asu) per box (1 means no helical symmetry is present)", "1"));
    helical_twist_initial = textToFloat(parser.getOption("--helical_twist_initial", "Helical twist (in degrees, positive values for right-handedness)", "0."));
    mymodel.helical_twist_min = textToFloat(parser.getOption("--helical_twist_min", "Minimum helical twist (in degrees, positive values for right-handedness)", "0."));
    mymodel.helical_twist_max = textToFloat(parser.getOption("--helical_twist_max", "Maximum helical twist (in degrees, positive values for right-handedness)", "0."));
    mymodel.helical_twist_inistep = textToFloat(parser.getOption("--helical_twist_inistep", "Initial step of helical twist search (in degrees)", "0."));
    helical_rise_initial = textToFloat(parser.getOption("--helical_rise_initial", "Helical rise (in Angstroms)", "0."));
    mymodel.helical_rise_min = textToFloat(parser.getOption("--helical_rise_min", "Minimum helical rise (in Angstroms)", "0."));
    mymodel.helical_rise_max = textToFloat(parser.getOption("--helical_rise_max", "Maximum helical rise (in Angstroms)", "0."));
    mymodel.helical_rise_inistep = textToFloat(parser.getOption("--helical_rise_inistep", "Initial step of helical rise search (in Angstroms)", "0."));
    helical_nstart = textToInteger(parser.getOption("--helical_nstart", "N-number for the N-start helix (only useful for rotational priors)", "1"));
    helical_z_percentage = textToFloat(parser.getOption("--helical_z_percentage", "This box length along the center of Z axis contains good information of the helix. Important in imposing and refining symmetry", "0.3"));
    helical_tube_inner_diameter = textToFloat(parser.getOption("--helical_inner_diameter", "Inner diameter of helical tubes in Angstroms (for masks of helical references and particles)", "-1."));
    helical_tube_outer_diameter = textToFloat(parser.getOption("--helical_outer_diameter", "Outer diameter of helical tubes in Angstroms (for masks of helical references and particles)", "-1."));
    do_helical_symmetry_local_refinement = parser.checkOption("--helical_symmetry_search", "Perform local refinement of helical symmetry?");
    helical_sigma_distance = textToFloat(parser.getOption("--helical_sigma_distance", "Sigma of distance along the helical tracks", "-1."));
    helical_keep_tilt_prior_fixed = parser.checkOption("--helical_keep_tilt_prior_fixed", "Keep helical tilt priors fixed (at 90 degrees) in global angular searches?");
    if (ignore_helical_symmetry) {
        mymodel.helical_nr_asu = 1; // IMPORTANT !
        do_helical_symmetry_local_refinement = false;
        helical_twist_initial = mymodel.helical_twist_min = mymodel.helical_twist_max = mymodel.helical_twist_inistep = 0.0;
        helical_rise_initial = mymodel.helical_rise_min = mymodel.helical_rise_max = mymodel.helical_rise_inistep = 0.0;
        helical_z_percentage = 0.0;
    }
    mymodel.initialiseHelicalParametersLists(helical_twist_initial, helical_rise_initial);
    mymodel.is_helix = do_helical_refine;
    if (mymodel.helical_rise_min > mymodel.helical_rise_max)
        std::swap(mymodel.helical_rise_min, mymodel.helical_rise_max);
    if (mymodel.helical_twist_min > mymodel.helical_twist_max)
        std::swap(mymodel.helical_twist_min, mymodel.helical_twist_max);
    helical_fourier_mask_resols = parser.getOption("--helical_exclude_resols", "Resolutions (in A) along helical axis to exclude from refinement (comma-separated pairs, e.g. 50-5)", "");
    fn_fourier_mask = parser.getOption("--fourier_mask", "Originally-sized, FFTW-centred image with Fourier mask for Projector", "None");

    // CTF, norm, scale, bfactor correction etc.
    int corrections_section = parser.addSection("Corrections");
    do_ctf_correction = parser.checkOption("--ctf", "Perform CTF correction?");
    do_ctf_padding = parser.checkOption("--pad_ctf", "Perform CTF padding to treat CTF aliaising better?");
    if (do_ctf_padding)
        REPORT_ERROR("--pad_ctf currently disabled.");
    intact_ctf_first_peak = parser.checkOption("--ctf_intact_first_peak", "Ignore CTFs until their first peak?");
    refs_are_ctf_corrected = parser.checkOption("--ctf_corrected_ref", "Have the input references been CTF-amplitude corrected?");
    ctf_phase_flipped = parser.checkOption("--ctf_phase_flipped", "Have the data been CTF phase-flipped?");
    only_flip_phases = parser.checkOption("--only_flip_phases", "Only perform CTF phase-flipping? (default is full amplitude-correction)");
    do_norm_correction = parser.checkOption("--norm", "Perform normalisation-error correction?");
    do_scale_correction = parser.checkOption("--scale", "Perform intensity-scale corrections on image groups?");
    // Allow switching off norm and scale (which is on by default in the GUI)
    if (parser.checkOption("--no_norm", "Switch off normalisation-error correction?"))
        do_norm_correction = false;
    if (parser.checkOption("--no_scale", "Switch off intensity-scale corrections on image groups?"))
        do_scale_correction = false;

    // SGD stuff
    int sgd_section = parser.addSection("Stochastic Gradient Descent");
    do_sgd = parser.checkOption("--sgd", "Perform stochastic gradient descent instead of default expectation-maximization");
    do_avoid_sgd = parser.checkOption("--stochastic_em", "Perform stochastic EM instead of SGD to avoid patent problems for initial model generation by commercial users");
    // Stochastic EM is implemented as a variant of SGD, though it is really a different algorithm!
    if (do_avoid_sgd) do_sgd = true;
    sgd_ini_iter = textToInteger(parser.getOption("--sgd_ini_iter", "Number of initial SGD iterations", "50"));
    sgd_fin_iter = textToInteger(parser.getOption("--sgd_fin_iter", "Number of final SGD iterations", "50"));
    sgd_inbetween_iter = textToInteger(parser.getOption("--sgd_inbetween_iter", "Number of SGD iterations between the initial and final ones", "200"));
    sgd_ini_resol = textToInteger(parser.getOption("--sgd_ini_resol", "Resolution cutoff during the initial SGD iterations (A)", "35"));
    sgd_fin_resol = textToInteger(parser.getOption("--sgd_fin_resol", "Resolution cutoff during the final SGD iterations (A)", "15"));
    sgd_ini_subset_size = textToInteger(parser.getOption("--sgd_ini_subset", "Mini-batch size during the initial SGD iterations", "100"));
    sgd_fin_subset_size = textToInteger(parser.getOption("--sgd_fin_subset", "Mini-batch size during the final SGD iterations", "500"));
    mu = textToFloat(parser.getOption("--mu", "Momentum parameter for SGD updates", "0.9"));
    sgd_stepsize = textToFloat(parser.getOption("--sgd_stepsize", "Step size parameter for SGD updates", "0.5"));
    sgd_sigma2fudge_ini = textToFloat(parser.getOption("--sgd_sigma2fudge_initial", "Initial factor by which the noise variance will be multiplied for SGD (not used if halftime is negative)", "8"));
    sgd_sigma2fudge_halflife = textToInteger(parser.getOption("--sgd_sigma2fudge_halflife", "Initialise SGD with 8x higher noise-variance, and reduce with this half-life in # of particles (default is keep normal variance)", "-1"));
    do_sgd_skip_anneal = parser.checkOption("--sgd_skip_anneal", "By default, multiple references are annealed during the in_between iterations. Use this option to switch annealing off");
    write_every_sgd_iter = textToInteger(parser.getOption("--sgd_write_iter", "Write out model every so many iterations in SGD (default is writing out all iters)", "1"));

    // Computation stuff
    // The number of threads is always read from the command line
    int computation_section = parser.addSection("Computation");
    x_pool = textToInteger(parser.getOption("--pool", "Number of images to pool for each thread task", "1"));
    nr_threads = textToInteger(parser.getOption("--j", "Number of threads to run in parallel (only useful on multi-core machines)", "1"));
    combine_weights_thru_disc = !parser.checkOption("--dont_combine_weights_via_disc", "Send the large arrays of summed weights through the MPI network, instead of writing large files to disc");
    do_shifts_onthefly = parser.checkOption("--onthefly_shifts", "Calculate shifted images on-the-fly, do not store precalculated ones in memory");
    do_parallel_disc_io = !parser.checkOption("--no_parallel_disc_io", "Do NOT let parallel (MPI) processes access the disc simultaneously (use this option with NFS)");
    do_preread_images  = parser.checkOption("--preread_images", "Use this to let the leader process read all particles into memory. Be careful you have enough RAM for large data sets!");
    fn_scratch = parser.getOption("--scratch_dir", "If provided, particle stacks will be copied to this local scratch disk prior to refinement.", "");
    keep_free_scratch_Gb = textToFloat(parser.getOption("--keep_free_scratch", "Space available for copying particle stacks (in Gb)", "10"));
    do_reuse_scratch = parser.checkOption("--reuse_scratch", "Re-use data on scratchdir, instead of wiping it and re-copying all data.");
    keep_scratch = parser.checkOption("--keep_scratch", "Don't remove scratch after convergence. Following jobs that use EXACTLY the same particles should use --reuse_scratch.");
    do_fast_subsets = parser.checkOption("--fast_subsets", "Use faster optimisation by using subsets of the data in the first 15 iterations");
    #ifdef ALTCPU
    do_cpu = parser.checkOption("--cpu", "Use intel vectorisation implementation for CPU");
    #else
    do_cpu = false;
    #endif

    do_gpu = parser.checkOption("--gpu", "Use available gpu resources for some calculations");
    gpu_ids = parser.getOption("--gpu", "Device ids for each MPI-thread","default");
    #ifndef CUDA
    if (do_gpu) {
        std::cerr << "+ WARNING : Relion was compiled without CUDA of at least version 7.0 - you do NOT have support for GPUs" << std::endl;
        do_gpu = false;
    }
    #endif
    double temp_reqSize = textToDouble(parser.getOption("--free_gpu_memory", "GPU device memory (in Mb) to leave free after allocation.", "0"));
    if (!do_zero_mask) { temp_reqSize += 100; }
    temp_reqSize *= 1000 * 1000;
    if (temp_reqSize < 0) {
        REPORT_ERROR("Invalid free_gpu_memory value.");
    } else {
        requested_free_gpu_memory =  temp_reqSize;
    }

    // Expert options
    int expert_section = parser.addSection("Expert options");
    mymodel.padding_factor = textToFloat(parser.getOption("--pad", "Oversampling factor for the Fourier transforms of the references", "2"));
    ref_angpix = textToFloat(parser.getOption("--ref_angpix", "Pixel size (in A) for the input reference (default is to read from header)", "-1."));
    mymodel.interpolator = parser.checkOption("--NN", "Perform nearest-neighbour instead of linear Fourier-space interpolation?") ? NEAREST_NEIGHBOUR : TRILINEAR;
    mymodel.r_min_nn = textToInteger(parser.getOption("--r_min_nn", "Minimum number of Fourier shells to perform linear Fourier-space interpolation", "10"));
    verb = textToInteger(parser.getOption("--verb", "Verbosity (1=normal, 0=silent)", "1"));
    random_seed = textToInteger(parser.getOption("--random_seed", "Number for the random seed generator", "-1"));
    max_coarse_size = textToInteger(parser.getOption("--coarse_size", "Maximum image size for the first pass of the adaptive sampling approach", "-1"));
    adaptive_fraction = textToFloat(parser.getOption("--adaptive_fraction", "Fraction of the weights to be considered in the first pass of adaptive oversampling ", "0.999"));
    width_mask_edge = textToInteger(parser.getOption("--maskedge", "Width of the soft edge of the spherical mask (in pixels)", "5"));
    fix_sigma_noise = parser.checkOption("--fix_sigma_noise", "Fix the experimental noise spectra?");
    fix_sigma_offset = parser.checkOption("--fix_sigma_offset", "Fix the stddev in the origin offsets?");
    incr_size = textToInteger(parser.getOption("--incr_size", "Number of Fourier shells beyond the current resolution to be included in refinement", "10"));
    do_print_metadata_labels = parser.checkOption("--print_metadata_labels", "Print a table with definitions of all metadata labels, and exit");
    do_print_symmetry_ops = parser.checkOption("--print_symmetry_ops", "Print all symmetry transformation matrices, and exit");
    strict_highres_exp = textToFloat(parser.getOption("--strict_highres_exp", "High resolution limit (in Angstrom) to restrict probability calculations in the expectation step", "-1"));
    strict_lowres_exp = textToFloat(parser.getOption("--strict_lowres_exp", "Low resolution limit (in Angstrom) to restrict probability calculations in the expectation step", "-1"));
    dont_raise_norm_error = parser.checkOption("--dont_check_norm", "Skip the check whether the images are normalised correctly");
    do_always_cc  = parser.checkOption("--always_cc", "Perform CC-calculation in all iterations (useful for faster denovo model generation?)");
    do_phase_random_fsc = parser.checkOption("--solvent_correct_fsc", "Correct FSC curve for the effects of the solvent mask?");
    do_skip_maximization = parser.checkOption("--skip_maximize", "Skip maximization step (only write out data.star file)?");
    failsafe_threshold = textToInteger(parser.getOption("--failsafe_threshold", "Maximum number of particles permitted to be handled by fail-safe mode, due to zero sum of weights, before exiting with an error (GPU only).", "40"));
    do_external_reconstruct = parser.checkOption("--external_reconstruct", "Perform the reconstruction step outside relion_refine, e.g. for learned priors?)");
    nr_iter_max = textToInteger(parser.getOption("--auto_iter_max", "In auto-refinement, stop at this iteration.", "999"));
    auto_ignore_angle_changes = parser.checkOption("--auto_ignore_angles", "In auto-refinement, update angular sampling regardless of changes in orientations for convergence. This makes convergence faster.");
    auto_resolution_based_angles= parser.checkOption("--auto_resol_angles", "In auto-refinement, update angular sampling based on resolution-based required sampling. This makes convergence faster.");
    allow_coarser_samplings = parser.checkOption("--allow_coarser_sampling", "In 2D/3D classification, allow coarser angular and translational samplings if accuracies are bad (typically in earlier iterations.");
    do_trust_ref_size = parser.checkOption("--trust_ref_size", "Trust the pixel and box size of the input reference; by default the program will die if these are different from the first optics group of the data");
    ///////////////// Special stuff for first iteration (only accessible via CL, not through readSTAR ////////////////////

    // When reading from the CL: always start at iteration 1 and subset 1
    iter = 0;
    // When starting from CL: always calculate initial sigma_noise
    do_calculate_initial_sigma_noise = true;
    // Start average norm correction at 1!
    mymodel.avg_norm_correction = 1.0;
    // Always initialise the PDF of the directions
    directions_have_changed = true;

    // Only reconstruct and join random halves are only available when continuing an old run
    do_join_random_halves = false;

    // For auto-sampling and convergence check
    nr_iter_wo_resol_gain = 0;
    nr_iter_wo_large_hidden_variable_changes = 0;
    current_changes_optimal_classes = 9999999;
    current_changes_optimal_offsets = 999.0;
    current_changes_optimal_orientations = 999.0;
    smallest_changes_optimal_classes = 9999999;
    smallest_changes_optimal_offsets = 999.0;
    smallest_changes_optimal_orientations = 999.0;
    acc_rot = acc_trans = 999.0;

    best_resol_thus_far = 1./999.0;
    has_converged = false;
    has_high_fsc_at_limit = false;
    has_large_incr_size_iter_ago = 0;
    do_initialise_bodies = false;

    // By default, start with nr_bodies to 1
    mymodel.nr_bodies = 1;
    fn_body_masks = "None";

    // Debugging/analysis/hidden stuff
    do_map = !checkParameter(argc, argv, "--no_map");
    minres_map = textToInteger(getParameter(argc, argv, "--minres_map", "5"));
    do_bfactor = checkParameter(argc, argv, "--bfactor");
    gridding_nr_iter = textToInteger(getParameter(argc, argv, "--gridding_iter", "10"));
    debug1 = textToFloat(getParameter(argc, argv, "--debug1", "0"));
    debug2 = textToFloat(getParameter(argc, argv, "--debug2", "0"));
    debug3 = textToFloat(getParameter(argc, argv, "--debug3", "0"));
    // Read in initial sigmaNoise spectrum
    fn_sigma = getParameter(argc, argv, "--sigma", "");
    do_calculate_initial_sigma_noise = fn_sigma.empty();
    sigma2_fudge = textToFloat(getParameter(argc, argv, "--sigma2_fudge", "1"));
    do_acc_currentsize_despite_highres_exp = checkParameter(argc, argv, "--accuracy_current_size");
    do_sequential_halves_recons  = checkParameter(argc, argv, "--sequential_halves_recons");
    do_always_join_random_halves = checkParameter(argc, argv, "--always_join_random_halves");
    do_use_all_data = checkParameter(argc, argv, "--use_all_data");
    do_only_sample_tilt  = checkParameter(argc, argv, "--only_sample_tilt");
    minimum_angular_sampling = textToFloat(getParameter(argc, argv, "--minimum_angular_sampling", "0"));
    maximum_angular_sampling = textToFloat(getParameter(argc, argv, "--maximum_angular_sampling", "0"));
    asymmetric_padding = parser.checkOption("--asymmetric_padding", "", "false", true);
    maximum_significants = textToInteger(parser.getOption("--maxsig", "Maximum number of poses & translations to consider", "-1"));
    skip_gridding = parser.checkOption("--skip_gridding", "Skip gridding in the M step");
    debug_split_random_half = textToInteger(getParameter(argc, argv, "--debug_split_random_half", "0"));

    #ifdef DEBUG_READ
    std::cerr<<"MlOptimiser::parseInitial Done"<<std::endl;
    #endif
}

void MlOptimiser::read(FileName fn_in, int rank, bool do_prevent_preread) {
    // #define DEBUG_READ
    #ifdef DEBUG_READ
    std::cerr<<"MlOptimiser::readStar entering ..."<<std::endl;
    #endif

    if (rank == 0)
        std::cout << " Reading in optimiser.star ..." << std::endl;

    // Open input file
    std::ifstream in(fn_in.data(), std::ios_base::in);
    if (in.fail())
        REPORT_ERROR( (std::string) "MlOptimiser::readStar: File " + fn_in + " cannot be read." );

    MetaDataTable MD;

    // Read general stuff
    FileName fn_model, fn_model2, fn_sampling;
    MD.readStar(in, "optimiser_general");
    in.close();

    const long int i = MD.size() - 1;
    try {
        fn_out = MD.getValue<std::string>(EMDL::OPTIMISER_OUTPUT_ROOTNAME, i);
        fn_model = MD.getValue<std::string>(EMDL::OPTIMISER_MODEL_STARFILE, i);
        fn_data = MD.getValue<std::string>(EMDL::OPTIMISER_DATA_STARFILE, i);
        fn_sampling = MD.getValue<std::string>(EMDL::OPTIMISER_SAMPLING_STARFILE, i);
        iter = MD.getValue<int>(EMDL::OPTIMISER_ITERATION_NO, i);
        nr_iter = MD.getValue<int>(EMDL::OPTIMISER_NR_ITERATIONS, i);
        do_split_random_halves = MD.getValue<bool>(EMDL::OPTIMISER_DO_SPLIT_RANDOM_HALVES, i);
        low_resol_join_halves = MD.getValue<RFLOAT>(EMDL::OPTIMISER_LOWRES_JOIN_RANDOM_HALVES, i);
        adaptive_oversampling = MD.getValue<int>(EMDL::OPTIMISER_ADAPTIVE_OVERSAMPLING, i);
        adaptive_fraction = MD.getValue<RFLOAT>(EMDL::OPTIMISER_ADAPTIVE_FRACTION, i);
        random_seed = MD.getValue<int>(EMDL::OPTIMISER_RANDOM_SEED, i);
        particle_diameter = MD.getValue<RFLOAT>(EMDL::OPTIMISER_PARTICLE_DIAMETER, i);
        width_mask_edge = MD.getValue<int>(EMDL::OPTIMISER_WIDTH_MASK_EDGE, i);
        do_zero_mask = MD.getValue<bool>(EMDL::OPTIMISER_DO_ZERO_MASK, i);
        do_solvent = MD.getValue<bool>(EMDL::OPTIMISER_DO_SOLVENT_FLATTEN, i);
        fn_mask = MD.getValue<std::string>(EMDL::OPTIMISER_SOLVENT_MASK_NAME, i);
        fn_mask2 = MD.getValue<std::string>(EMDL::OPTIMISER_SOLVENT_MASK2_NAME, i);
        fn_tau = MD.getValue<std::string>(EMDL::OPTIMISER_TAU_SPECTRUM_NAME, i);
        max_coarse_size = MD.getValue<int>(EMDL::OPTIMISER_MAX_COARSE_SIZE, i);
        strict_highres_exp = MD.getValue<RFLOAT>(EMDL::OPTIMISER_HIGHRES_LIMIT_EXP, i);
        incr_size = MD.getValue<int>(EMDL::OPTIMISER_INCR_SIZE, i);
        do_map = MD.getValue<bool>(EMDL::OPTIMISER_DO_MAP, i);
        do_auto_refine = MD.getValue<bool>(EMDL::OPTIMISER_DO_AUTO_REFINE, i);
        autosampling_hporder_local_searches = MD.getValue<int>(EMDL::OPTIMISER_AUTO_LOCAL_HP_ORDER, i);
        nr_iter_wo_resol_gain = MD.getValue<int>(EMDL::OPTIMISER_NR_ITER_WO_RESOL_GAIN, i);
        best_resol_thus_far = MD.getValue<RFLOAT>(EMDL::OPTIMISER_BEST_RESOL_THUS_FAR, i);
        nr_iter_wo_large_hidden_variable_changes = MD.getValue<int>(EMDL::OPTIMISER_NR_ITER_WO_HIDDEN_VAR_CHANGES, i);
        do_skip_align = MD.getValue<bool>(EMDL::OPTIMISER_DO_SKIP_ALIGN, i);
        // do_skip_rotate = MD.getValue<bool>(EMDL::OPTIMISER_DO_SKIP_ROTATE, i);
        acc_rot = MD.getValue<RFLOAT>(EMDL::OPTIMISER_ACCURACY_ROT, i);
        current_changes_optimal_orientations = MD.getValue<RFLOAT>(EMDL::OPTIMISER_CHANGES_OPTIMAL_ORIENTS, i);
        current_changes_optimal_offsets = MD.getValue<RFLOAT>(EMDL::OPTIMISER_CHANGES_OPTIMAL_OFFSETS, i);
        current_changes_optimal_classes = MD.getValue<RFLOAT>(EMDL::OPTIMISER_CHANGES_OPTIMAL_CLASSES, i);
        smallest_changes_optimal_orientations = MD.getValue<RFLOAT>(EMDL::OPTIMISER_SMALLEST_CHANGES_OPT_ORIENTS, i);
        smallest_changes_optimal_offsets = MD.getValue<RFLOAT>(EMDL::OPTIMISER_SMALLEST_CHANGES_OPT_OFFSETS, i);
        smallest_changes_optimal_classes = MD.getValue<RFLOAT>(EMDL::OPTIMISER_SMALLEST_CHANGES_OPT_CLASSES, i);
        has_converged = MD.getValue<bool>(EMDL::OPTIMISER_HAS_CONVERGED, i);
        has_high_fsc_at_limit = MD.getValue<bool>(EMDL::OPTIMISER_HAS_HIGH_FSC_AT_LIMIT, i);
        has_large_incr_size_iter_ago = MD.getValue<bool>(EMDL::OPTIMISER_HAS_LARGE_INCR_SIZE_ITER_AGO, i);
        do_norm_correction = MD.getValue<bool>(EMDL::OPTIMISER_DO_CORRECT_NORM, i);
        do_scale_correction = MD.getValue<bool>(EMDL::OPTIMISER_DO_CORRECT_SCALE, i);
        do_ctf_correction = MD.getValue<bool>(EMDL::OPTIMISER_DO_CORRECT_CTF, i);
        intact_ctf_first_peak = MD.getValue<bool>(EMDL::OPTIMISER_IGNORE_CTF_UNTIL_FIRST_PEAK, i);
        ctf_phase_flipped = MD.getValue<bool>(EMDL::OPTIMISER_DATA_ARE_CTF_PHASE_FLIPPED, i);
        only_flip_phases = MD.getValue<bool>(EMDL::OPTIMISER_DO_ONLY_FLIP_CTF_PHASES, i);
        refs_are_ctf_corrected = MD.getValue<bool>(EMDL::OPTIMISER_REFS_ARE_CTF_CORRECTED, i);
        fix_sigma_noise = MD.getValue<bool>(EMDL::OPTIMISER_FIX_SIGMA_NOISE, i);
        fix_sigma_offset = MD.getValue<bool>(EMDL::OPTIMISER_FIX_SIGMA_OFFSET, i);
        nr_pool = MD.getValue<int>(EMDL::OPTIMISER_MAX_NR_POOL, i);
    } catch (const char* errmsg) {
        REPORT_ERROR("MlOptimiser::readStar: incorrect optimiser_general table");
    }

    // Backwards compatibility with RELION-1.4
    fn_local_symmetry                    = tryget<std::string>(MD, EMDL::OPTIMISER_LOCAL_SYMMETRY_FILENAME,          "None", i);
    do_helical_refine                    = tryget<bool>       (MD, EMDL::OPTIMISER_DO_HELICAL_REFINE,                 false, i);
    ignore_helical_symmetry              = tryget<bool>       (MD, EMDL::OPTIMISER_IGNORE_HELICAL_SYMMETRY,           false, i);
    helical_twist_initial                = tryget<RFLOAT>     (MD, EMDL::OPTIMISER_HELICAL_TWIST_INITIAL,              0.0, i);
    helical_rise_initial                 = tryget<RFLOAT>     (MD, EMDL::OPTIMISER_HELICAL_RISE_INITIAL,               0.0, i);
    helical_z_percentage                 = tryget<RFLOAT>     (MD, EMDL::OPTIMISER_HELICAL_Z_PERCENTAGE,               0.3, i);
    helical_nstart                       = tryget<int>        (MD, EMDL::OPTIMISER_HELICAL_NSTART,                     1, i);
    helical_tube_inner_diameter          = tryget<RFLOAT>     (MD, EMDL::OPTIMISER_HELICAL_TUBE_INNER_DIAMETER,       -1.0, i);
    helical_tube_outer_diameter          = tryget<RFLOAT>     (MD, EMDL::OPTIMISER_HELICAL_TUBE_OUTER_DIAMETER,       -1.0, i);
    do_helical_symmetry_local_refinement = tryget<bool>       (MD, EMDL::OPTIMISER_HELICAL_SYMMETRY_LOCAL_REFINEMENT, false, i);
    helical_sigma_distance               = tryget<RFLOAT>     (MD, EMDL::OPTIMISER_HELICAL_SIGMA_DISTANCE,            -1.0, i);
    helical_keep_tilt_prior_fixed        = tryget<bool>       (MD, EMDL::OPTIMISER_HELICAL_KEEP_TILT_PRIOR_FIXED,     false, i);

    // New SGD (13 Feb 2018)
    do_sgd                   = tryget<bool>       (MD, EMDL::OPTIMISER_DO_SGD,                    false, i);
    do_avoid_sgd             = tryget<bool>       (MD, EMDL::OPTIMISER_DO_STOCHASTIC_EM,          false, i);
    sgd_ini_iter             = tryget<int>        (MD, EMDL::OPTIMISER_SGD_INI_ITER,              50, i);
    sgd_fin_iter             = tryget<int>        (MD, EMDL::OPTIMISER_SGD_FIN_ITER,              50, i);
    sgd_inbetween_iter       = tryget<int>        (MD, EMDL::OPTIMISER_SGD_INBETWEEN_ITER,        200, i);
    sgd_ini_resol            = tryget<RFLOAT>     (MD, EMDL::OPTIMISER_SGD_INI_RESOL,             35.0, i);
    sgd_fin_resol            = tryget<RFLOAT>     (MD, EMDL::OPTIMISER_SGD_FIN_RESOL,             15.0, i);
    sgd_ini_subset_size      = tryget<int>        (MD, EMDL::OPTIMISER_SGD_INI_SUBSET_SIZE,       100, i);
    sgd_fin_subset_size      = tryget<int>        (MD, EMDL::OPTIMISER_SGD_FIN_SUBSET_SIZE,       500, i);
    mu                       = tryget<RFLOAT>     (MD, EMDL::OPTIMISER_SGD_MU,                    0.9, i);
    sgd_sigma2fudge_ini      = tryget<RFLOAT>     (MD, EMDL::OPTIMISER_SGD_SIGMA2FUDGE_INI,       8.0, i);
    sgd_sigma2fudge_halflife = tryget<long>       (MD, EMDL::OPTIMISER_SGD_SIGMA2FUDGE_HALFLIFE, -1, i);
    do_sgd_skip_anneal       = tryget<bool>       (MD, EMDL::OPTIMISER_SGD_SKIP_ANNNEAL,         false, i);
    subset_size              = tryget<long>       (MD, EMDL::OPTIMISER_SGD_SUBSET_SIZE,          -1, i);
    sgd_stepsize             = tryget<RFLOAT>     (MD, EMDL::OPTIMISER_SGD_STEPSIZE,              0.5, i);
    write_every_sgd_iter     = tryget<int>        (MD, EMDL::OPTIMISER_SGD_WRITE_EVERY_SUBSET,    1, i);
    fn_body_masks            = tryget<std::string>(MD, EMDL::BODY_STAR_FILE,                     "None", i);
    do_phase_random_fsc      = tryget<bool>       (MD, EMDL::OPTIMISER_DO_SOLVENT_FSC,           false, i);
    do_fast_subsets          = tryget<bool>       (MD, EMDL::OPTIMISER_FAST_SUBSETS,             false, i);
    do_external_reconstruct  = tryget<bool>       (MD, EMDL::OPTIMISER_DO_EXTERNAL_RECONSTRUCT,  false, i);

    // Backwards compatibility with RELION 3.0
    try {
        acc_trans = MD.getValue<RFLOAT>(EMDL::OPTIMISER_ACCURACY_TRANS_ANGSTROM, i);
    } catch (const char* errmsg) { try {
        acc_trans = MD.getValue<RFLOAT>(EMDL::OPTIMISER_ACCURACY_TRANS, i);
    } catch (const char* errmsg) {
        REPORT_ERROR("MlOptimiser::readStar::ERROR no accuracy of translations defined!");
    } }

    try {
        fn_fourier_mask = MD.getValue<std::string>(EMDL::OPTIMISER_FOURIER_MASK, i);
    } catch (const char* errmsg) {
        fn_fourier_mask = "None";
    }

    if (do_split_random_halves) {
        try {
            fn_model2 = MD.getValue<std::string>(EMDL::OPTIMISER_MODEL_STARFILE2, i);
        } catch (const char* errmsg) {
            REPORT_ERROR("MlOptimiser::readStar: splitting data into two random halves, but rlnModelStarFile2 not found in optimiser_general table");
        }
        if (fn_model2.empty())
            REPORT_ERROR("MlOptimiser::readStar: splitting data into two random halves, but rlnModelStarFile2 is empty. Probably you specified an optimiser STAR file generated with --force_converge. You cannot perform continuation or subtraction from this file. Please use one from the previous iteration.");
        try {
            strict_lowres_exp = MD.getValue<RFLOAT>(EMDL::OPTIMISER_LOWRES_LIMIT_EXP, i);
        } catch (const char* errmsg) {
            strict_lowres_exp = -1.0;
        }
    }

    // Initialise some stuff for first-iteration only (not relevant here...0)
    do_calculate_initial_sigma_noise = false;
    do_average_unaligned = false;
    do_generate_seeds = false;
    do_firstiter_cc = false;
    ini_high = 0;

    // Initialise some of the other, hidden or debugging stuff
    minres_map = 5;
    do_bfactor = false;
    gridding_nr_iter = 10;
    debug1 = debug2 = debug3 = 0.0;

    // Then read in sampling, mydata and mymodel stuff
    // If do_preread_images: when not do_parallel_disc_io: only the leader reads all images into RAM; otherwise: everyone reads in images into RAM
    #ifdef DEBUG_READ
    std::cerr<<"MlOptimiser::readStar before data."<<std::endl;
    #endif
    bool do_preread = do_preread_images && !do_prevent_preread && (do_parallel_disc_io || rank == 0);
    bool is_helical_segment = do_helical_refine || mymodel.ref_dim == 2 && helical_tube_outer_diameter > 0.0;
    mydata.read(fn_data, false, false, do_preread, is_helical_segment);

    #ifdef DEBUG_READ
    std::cerr << "MlOptimiser::readStar before model." << std::endl;
    #endif
    if (do_split_random_halves) {
        if (debug_split_random_half == 1) {
            mymodel.read(fn_model);
        } else if (debug_split_random_half == 2) {
            mymodel.read(fn_model2);
        } else if (rank % 2 == 1) {
            mymodel.read(fn_model);
        } else {
            mymodel.read(fn_model2);
        }
    } else {
        mymodel.read(fn_model);
    }
    // Set up the bodies in the model, if this is a continuation of a multibody refinement (otherwise this is done in initialiseGeneral)
    if (fn_body_masks != "None") {
        mymodel.initialiseBodies(fn_body_masks, fn_out, false, rank); // also_initialise_rest = false

        if (mymodel.nr_bodies != mydata.nr_bodies)
            REPORT_ERROR("ERROR: Unequal number of bodies in model.star and data.star files!");
    }

    #ifdef DEBUG_READ
    std::cerr<<"MlOptimiser::readStar before sampling."<<std::endl;
    #endif
    sampling.read(fn_sampling);

    #ifdef DEBUG_READ
    std::cerr<<"MlOptimiser::readStar done."<<std::endl;
    #endif
}


void MlOptimiser::write(bool do_write_sampling, bool do_write_data, bool do_write_optimiser, bool do_write_model, int random_subset) {
    if (subset_size > 0 && (iter % write_every_sgd_iter) != 0 && iter != nr_iter)
        return;

    FileName fn_tmp, fn_model, fn_model2, fn_data, fn_sampling;
    FileName fn_root = iter > -1 ? FileName::compose(fn_out + "_it", iter, "", 3) : fn_out;
    // fn_root2 is used to write out the model and optimiser, and adds a subset number in SGD
    FileName fn_root2 = fn_root;
    bool do_write_bild = !do_skip_align && !do_skip_rotate && !do_sgd;

    // First write "main" STAR file with all information from this run
    // Do this for random_subset==0 and random_subset==1
    if (do_write_optimiser && random_subset < 2) {
        fn_tmp = fn_root2 + "_optimiser.star";
        std::ofstream fh((fn_tmp).c_str(), std::ios::out);
        if (!fh) REPORT_ERROR((std::string) "MlOptimiser::write: Cannot write file: " + fn_tmp);

        // Write the command line as a comment in the header
        fh << "# RELION optimiser; version " << g_RELION_VERSION <<std::endl;
        fh << "# ";
        parser.writeCommandLine(fh);

        if (do_split_random_halves && !do_join_random_halves) {
            fn_model  = fn_root2 + "_half1_model.star";
            fn_model2 = fn_root2 + "_half2_model.star";
        } else {
            fn_model = fn_root2 + "_model.star";
        }
        fn_data     = fn_root + "_data.star";
        fn_sampling = fn_root + "_sampling.star";

        MetaDataTable MD;
        MD.isList = true;
        MD.name = "optimiser_general";
        MD.addObject();
        const long int i = MD.size() - 1;
        MD.setValue(EMDL::OPTIMISER_OUTPUT_ROOTNAME, fn_out, i);
        if (do_split_random_halves) {
            MD.setValue(EMDL::OPTIMISER_MODEL_STARFILE, fn_model, i);
            MD.setValue(EMDL::OPTIMISER_MODEL_STARFILE2, fn_model2, i);
        } else {
            MD.setValue(EMDL::OPTIMISER_MODEL_STARFILE, fn_model, i);
        }
        MD.setValue(EMDL::OPTIMISER_DATA_STARFILE, fn_data, i);
        MD.setValue(EMDL::OPTIMISER_SAMPLING_STARFILE, fn_sampling, i);
        MD.setValue(EMDL::OPTIMISER_ITERATION_NO, iter, i);
        MD.setValue(EMDL::OPTIMISER_NR_ITERATIONS, nr_iter, i);
        MD.setValue(EMDL::OPTIMISER_DO_SPLIT_RANDOM_HALVES, do_split_random_halves, i);
        MD.setValue(EMDL::OPTIMISER_LOWRES_JOIN_RANDOM_HALVES, low_resol_join_halves, i);
        MD.setValue(EMDL::OPTIMISER_ADAPTIVE_OVERSAMPLING, adaptive_oversampling, i);
        MD.setValue(EMDL::OPTIMISER_ADAPTIVE_FRACTION, adaptive_fraction, i);
        MD.setValue(EMDL::OPTIMISER_RANDOM_SEED, random_seed, i);
        MD.setValue(EMDL::OPTIMISER_PARTICLE_DIAMETER, particle_diameter, i);
        MD.setValue(EMDL::OPTIMISER_WIDTH_MASK_EDGE, width_mask_edge, i);
        MD.setValue(EMDL::OPTIMISER_DO_ZERO_MASK, do_zero_mask, i);
        MD.setValue(EMDL::OPTIMISER_DO_SOLVENT_FLATTEN, do_solvent, i);
        MD.setValue(EMDL::OPTIMISER_DO_SOLVENT_FSC, do_phase_random_fsc, i);
        MD.setValue(EMDL::OPTIMISER_SOLVENT_MASK_NAME, fn_mask, i);
        MD.setValue(EMDL::OPTIMISER_SOLVENT_MASK2_NAME, fn_mask2, i);
        MD.setValue(EMDL::BODY_STAR_FILE, fn_body_masks, i);
        MD.setValue(EMDL::OPTIMISER_TAU_SPECTRUM_NAME, fn_tau, i);
        MD.setValue(EMDL::OPTIMISER_MAX_COARSE_SIZE, max_coarse_size, i);
        MD.setValue(EMDL::OPTIMISER_HIGHRES_LIMIT_EXP, strict_highres_exp, i);
        MD.setValue(EMDL::OPTIMISER_LOWRES_LIMIT_EXP, strict_lowres_exp, i);
        MD.setValue(EMDL::OPTIMISER_INCR_SIZE, incr_size, i);
        MD.setValue(EMDL::OPTIMISER_DO_MAP, do_map, i);
        MD.setValue(EMDL::OPTIMISER_FAST_SUBSETS, do_fast_subsets, i);
        MD.setValue(EMDL::OPTIMISER_DO_EXTERNAL_RECONSTRUCT, do_external_reconstruct, i);
        MD.setValue(EMDL::OPTIMISER_DO_SGD, do_sgd, i);
        MD.setValue(EMDL::OPTIMISER_DO_STOCHASTIC_EM, do_avoid_sgd, i);
        MD.setValue(EMDL::OPTIMISER_SGD_INI_ITER, sgd_ini_iter, i);
        MD.setValue(EMDL::OPTIMISER_SGD_FIN_ITER, sgd_fin_iter, i);
        MD.setValue(EMDL::OPTIMISER_SGD_INBETWEEN_ITER, sgd_inbetween_iter, i);
        MD.setValue(EMDL::OPTIMISER_SGD_INI_RESOL, sgd_ini_resol, i);
        MD.setValue(EMDL::OPTIMISER_SGD_FIN_RESOL, sgd_fin_resol, i);
        MD.setValue(EMDL::OPTIMISER_SGD_INI_SUBSET_SIZE, sgd_ini_subset_size, i);
        MD.setValue(EMDL::OPTIMISER_SGD_FIN_SUBSET_SIZE, sgd_fin_subset_size, i);
        MD.setValue(EMDL::OPTIMISER_SGD_MU, mu, i);
        MD.setValue(EMDL::OPTIMISER_SGD_SIGMA2FUDGE_INI, sgd_sigma2fudge_ini, i);
        MD.setValue(EMDL::OPTIMISER_SGD_SIGMA2FUDGE_HALFLIFE, sgd_sigma2fudge_halflife, i);
        MD.setValue(EMDL::OPTIMISER_SGD_SKIP_ANNNEAL, do_sgd_skip_anneal, i);
        MD.setValue(EMDL::OPTIMISER_SGD_SUBSET_SIZE, subset_size, i);
        MD.setValue(EMDL::OPTIMISER_SGD_WRITE_EVERY_SUBSET, write_every_sgd_iter, i);
        MD.setValue(EMDL::OPTIMISER_SGD_STEPSIZE, sgd_stepsize, i);
        MD.setValue(EMDL::OPTIMISER_DO_AUTO_REFINE, do_auto_refine, i);
        MD.setValue(EMDL::OPTIMISER_AUTO_LOCAL_HP_ORDER, autosampling_hporder_local_searches, i);
        MD.setValue(EMDL::OPTIMISER_NR_ITER_WO_RESOL_GAIN, nr_iter_wo_resol_gain, i);
        MD.setValue(EMDL::OPTIMISER_BEST_RESOL_THUS_FAR,best_resol_thus_far, i);
        MD.setValue(EMDL::OPTIMISER_NR_ITER_WO_HIDDEN_VAR_CHANGES, nr_iter_wo_large_hidden_variable_changes, i);
        MD.setValue(EMDL::OPTIMISER_DO_SKIP_ALIGN, do_skip_align, i);
        MD.setValue(EMDL::OPTIMISER_DO_SKIP_ROTATE, do_skip_rotate, i);
        MD.setValue(EMDL::OPTIMISER_ACCURACY_ROT, acc_rot, i);
        MD.setValue(EMDL::OPTIMISER_ACCURACY_TRANS_ANGSTROM, acc_trans, i);
        MD.setValue(EMDL::OPTIMISER_CHANGES_OPTIMAL_ORIENTS, current_changes_optimal_orientations, i);
        MD.setValue(EMDL::OPTIMISER_CHANGES_OPTIMAL_OFFSETS, current_changes_optimal_offsets, i);
        MD.setValue(EMDL::OPTIMISER_CHANGES_OPTIMAL_CLASSES, current_changes_optimal_classes, i);
        MD.setValue(EMDL::OPTIMISER_SMALLEST_CHANGES_OPT_ORIENTS, smallest_changes_optimal_orientations, i);
        MD.setValue(EMDL::OPTIMISER_SMALLEST_CHANGES_OPT_OFFSETS, smallest_changes_optimal_offsets, i);
        MD.setValue(EMDL::OPTIMISER_SMALLEST_CHANGES_OPT_CLASSES, smallest_changes_optimal_classes, i);
        MD.setValue(EMDL::OPTIMISER_LOCAL_SYMMETRY_FILENAME, fn_local_symmetry, i);
        MD.setValue(EMDL::OPTIMISER_DO_HELICAL_REFINE, do_helical_refine, i);
        MD.setValue(EMDL::OPTIMISER_IGNORE_HELICAL_SYMMETRY, ignore_helical_symmetry, i);
        MD.setValue(EMDL::OPTIMISER_FOURIER_MASK, fn_fourier_mask, i);
        MD.setValue(EMDL::OPTIMISER_HELICAL_TWIST_INITIAL, helical_twist_initial, i);
        MD.setValue(EMDL::OPTIMISER_HELICAL_RISE_INITIAL, helical_rise_initial, i);
        MD.setValue(EMDL::OPTIMISER_HELICAL_Z_PERCENTAGE, helical_z_percentage, i);
        MD.setValue(EMDL::OPTIMISER_HELICAL_NSTART, helical_nstart, i);
        MD.setValue(EMDL::OPTIMISER_HELICAL_TUBE_INNER_DIAMETER, helical_tube_inner_diameter, i);
        MD.setValue(EMDL::OPTIMISER_HELICAL_TUBE_OUTER_DIAMETER, helical_tube_outer_diameter, i);
        MD.setValue(EMDL::OPTIMISER_HELICAL_SYMMETRY_LOCAL_REFINEMENT, do_helical_symmetry_local_refinement, i);
        MD.setValue(EMDL::OPTIMISER_HELICAL_SIGMA_DISTANCE, helical_sigma_distance, i);
        MD.setValue(EMDL::OPTIMISER_HELICAL_KEEP_TILT_PRIOR_FIXED, helical_keep_tilt_prior_fixed, i);
        MD.setValue(EMDL::OPTIMISER_HAS_CONVERGED, has_converged, i);
        MD.setValue(EMDL::OPTIMISER_HAS_HIGH_FSC_AT_LIMIT, has_high_fsc_at_limit, i);
        MD.setValue(EMDL::OPTIMISER_HAS_LARGE_INCR_SIZE_ITER_AGO, has_large_incr_size_iter_ago, i);
        MD.setValue(EMDL::OPTIMISER_DO_CORRECT_NORM, do_norm_correction, i);
        MD.setValue(EMDL::OPTIMISER_DO_CORRECT_SCALE, do_scale_correction, i);
        MD.setValue(EMDL::OPTIMISER_DO_CORRECT_CTF, do_ctf_correction, i);
        MD.setValue(EMDL::OPTIMISER_IGNORE_CTF_UNTIL_FIRST_PEAK, intact_ctf_first_peak, i);
        MD.setValue(EMDL::OPTIMISER_DATA_ARE_CTF_PHASE_FLIPPED, ctf_phase_flipped, i);
        MD.setValue(EMDL::OPTIMISER_DO_ONLY_FLIP_CTF_PHASES, only_flip_phases, i);
        MD.setValue(EMDL::OPTIMISER_REFS_ARE_CTF_CORRECTED, refs_are_ctf_corrected, i);
        MD.setValue(EMDL::OPTIMISER_FIX_SIGMA_NOISE, fix_sigma_noise, i);
        MD.setValue(EMDL::OPTIMISER_FIX_SIGMA_OFFSET, fix_sigma_offset, i);
        MD.setValue(EMDL::OPTIMISER_MAX_NR_POOL, nr_pool, i);

        MD.write(fh);
    }

    // Then write the mymodel to file
    if (do_write_model)
    mymodel.write(
        do_split_random_halves && !do_join_random_halves ? fn_root2 + "_half" + integerToString(random_subset) :
                                                           fn_root2,
        sampling, do_write_bild
    );

    // And write the mydata to file
    if (do_write_data)
        mydata.write(fn_root);

    // And write the sampling object
    if (do_write_sampling)
        sampling.write(fn_root);
}

/** ========================== Initialisation  =========================== */

void MlOptimiser::initialise() {
    #ifdef DEBUG
    std::cerr << "MlOptimiser::initialise Entering" << std::endl;
    #endif

    if (do_gpu) {
        #ifdef CUDA
        int devCount;
        HANDLE_ERROR(cudaGetDeviceCount(&devCount));

        cudaDeviceProp deviceProp;
        int compatibleDevices(0);
        // Send device count seen by this follower
        HANDLE_ERROR(cudaGetDeviceCount(&devCount));
        for (int i = 0; i < devCount; i++) {
            HANDLE_ERROR(cudaGetDeviceProperties(&deviceProp, i));
            if (
                deviceProp.major > CUDA_CC_MAJOR ||
                deviceProp.major == CUDA_CC_MAJOR && deviceProp.minor >= CUDA_CC_MINOR
            ) {
                compatibleDevices += 1;
            } else {
                // std::cout << "Found a " << deviceProp.name << " GPU with compute-capability " << deviceProp.major << "." << deviceProp.minor << std::endl;
            }
        }
        if (compatibleDevices == 0) {
            REPORT_ERROR("You have no GPUs compatible with RELION (CUDA-capable and compute-capability >= 3.5");
        } else if (compatibleDevices != devCount) {
            std::cerr << "WARNING : at least one of your GPUs is not compatible with RELION (CUDA-capable and compute-capability >= 3.5)" << std::endl;
        }

        std::vector <std::vector<std::string>> allThreadIDs;
        untangleDeviceIDs(gpu_ids, allThreadIDs);

        // Sequential initialisation of GPUs on all ranks
        bool fullAutomaticMapping = true;
        bool semiAutomaticMapping = true;
        if (allThreadIDs[0].empty() || !std::isdigit(*gpu_ids.begin())) {
            std::cout << "gpu-ids not specified, threads will automatically be mapped to devices (incrementally)." << std::endl;
        } else {
            fullAutomaticMapping = false;
            if (allThreadIDs[0].size() != nr_threads) {
                std::cout << " Will distribute threads over devices ";
                for (int j = 0; j < allThreadIDs[0].size(); j++)
                    std::cout << " "  << allThreadIDs[0][j];
                std::cout  << std::endl;
            } else {
                semiAutomaticMapping = false;
                std::cout << " Using explicit indexing to assign devices ";
                for (int j = 0; j < allThreadIDs[0].size(); j++)
                    std::cout << " "  << allThreadIDs[0][j];
                std::cout << std::endl;
            }
        }

        for (int i = 0; i < nr_threads; i++) {
            int dev_id = semiAutomaticMapping ? (
                fullAutomaticMapping ?
                // Sjors: hack to make use of several cards; will only work if all MPI followers are on the same node!
                // Bjorn: Better hack
                devCount * i / nr_threads :
                textToInteger(allThreadIDs[0][i % allThreadIDs[0].size()].c_str())
                // textToInteger(allThreadIDs[0][allThreadIDs[0].size() * i / nr_threads].c_str());
            ) : textToInteger(allThreadIDs[0][i].c_str());  // not semiAutomatic => explicit


            std::cout << " Thread " << i << " mapped to device " << dev_id << std::endl;

            //Only make a new bundle of not existing on device
            int bundleId(-1);

            for (int j = 0; j < cudaDevices.size(); j++)
                if (cudaDevices[j] == dev_id)
                    bundleId = j;

            if (bundleId == -1) {
                bundleId = cudaDevices.size();
                cudaDevices.push_back(dev_id);
            }

            cudaOptimiserDeviceMap.push_back(bundleId);
        }
        mymodel.do_gpu = do_gpu;
        #else
        REPORT_ERROR("GPU usage requested, but RELION was compiled without CUDA support");
        #endif
    }

    initialiseGeneral();

    initialiseWorkLoad();

    #ifdef MKLFFT
    // Enable multi-threaded FFTW
    int success = fftw_init_threads();
    if (0 == success)
        REPORT_ERROR("Multithreaded FFTW failed to initialize");

    // And allow plans before expectation to run using allowed
    // number of threads
    fftw_plan_with_nthreads(nr_threads);
    #endif

    if (!fn_sigma.empty()) {
        // Read in sigma_noise spetrum from file DEVELOPMENTAL!!! FOR DEBUGGING ONLY....
        MetaDataTable MDsigma;
        MDsigma.read(fn_sigma);
        int idx;
        for (long int i : MDsigma) {
            idx = MDsigma.getValue<int>(EMDL::SPECTRAL_IDX, i);
            RFLOAT val = MDsigma.getValue<RFLOAT>(EMDL::MLMODEL_SIGMA2_NOISE, i);
            if (idx < Xsize(mymodel.sigma2_noise[0]))
                mymodel.sigma2_noise[0].elem(idx) = val;
        }
        if (idx < Xsize(mymodel.sigma2_noise[0]) - 1) {
            if (verb > 0)
                std::cout
                << " WARNING: provided sigma2_noise-spectrum has fewer entries ("
                << idx + 1 << ") than needed (" << Xsize(mymodel.sigma2_noise[0])
                << "). Set rest to zero..." << std::endl;
        }

        mydata.getNumberOfImagesPerGroup(mymodel.nr_particles_per_group);
        for (int igroup = 0; igroup < mymodel.nr_groups; igroup++) {
            // Use the same spectrum for all classes
            mymodel.sigma2_noise[igroup] =  mymodel.sigma2_noise[0];
            // We set wsum_model.sumw_group as in calculateSumOfPowerSpectraAndAverageImage
            wsum_model.sumw_group[igroup] = mymodel.nr_particles_per_group[igroup];
        }
    } else if (do_calculate_initial_sigma_noise || do_average_unaligned) {
        MultidimArray<RFLOAT> Mavg;

        // Calculate initial sigma noise model from power_class spectra of the individual images
        calculateSumOfPowerSpectraAndAverageImage(Mavg);

        // Abort through the pipeline_control system
        if (pipeline_control_check_abort_job())
            exit(RELION_EXIT_ABORTED);

        // Set sigma2_noise and Iref from averaged poser spectra and Mavg
        setSigmaNoiseEstimatesAndSetAverageImage(Mavg);
    }

    // First low-pass filter the initial references
    if (iter == 0)
        initialLowPassFilterReferences();

    // Initialise the data_versus_prior ratio to get the initial current_size right
    if (iter == 0 && !do_initialise_bodies)
        mymodel.initialiseDataVersusPrior(fix_tau); // fix_tau was set in initialiseGeneral

    // Check minimum group size of 10 particles
    if (verb > 0) {
        bool do_warn = false;
        for (int igroup = 0; igroup < mymodel.nr_groups; igroup++) {
            if (mymodel.nr_particles_per_group[igroup] < 10) {
                std:: cout << "WARNING: There are only " << mymodel.nr_particles_per_group[igroup] << " particles in group " << igroup + 1 << std::endl;
                do_warn = true;
            }
        }
        if (do_warn) {
            std:: cout << "WARNING: You may want to consider joining some micrographs into larger groups to obtain more robust noise estimates. " << std::endl;
            std:: cout << "         You can do so by using the same rlnMicrographName label for particles from multiple different micrographs in the input STAR file. " << std::endl;
        }
    }

    // Write out initial mymodel
    write(DONT_WRITE_SAMPLING, DO_WRITE_DATA, DO_WRITE_OPTIMISER, DO_WRITE_MODEL, 0);


    #ifdef DEBUG
    std::cerr<<"MlOptimiser::initialise Done"<<std::endl;
    #endif
}

void MlOptimiser::checkMask(FileName &_fn_mask, int solvent_nr, int rank) {
    int ref_box_size = Xsize(mymodel.Iref[0]);

    Image<RFLOAT> Isolvent;
    Isolvent.read(_fn_mask);
    Isolvent().setXmippOrigin();
    RFLOAT mask_pixel_size = Isolvent.header.getValue<RFLOAT>(EMDL::IMAGE_SAMPLINGRATE_X, Isolvent.header.size() - 1);

    bool need_new_mask = false;
    if (fabs(mask_pixel_size-mymodel.pixel_size) > 0.001) {
        need_new_mask = true;

        if (verb > 0) {
            std::cerr << " + WARNING: solvent mask pixel size: " << mask_pixel_size
            << " is not the same as the reference pixel size: " << mymodel.pixel_size << std::endl;
            std::cerr << " + WARNING: re-scaling the mask... " << std::endl;
        }

        if (rank == 0) {
            // Only the leader writes out the new mask
            int rescale_size = round(Xsize(Isolvent()) * mask_pixel_size / mymodel.pixel_size);
            rescale_size += rescale_size % 2;  // Ensure divisibility by 2 (+= rather than -=)
            resizeMap(Isolvent(), rescale_size);
            Isolvent.setSamplingRateInHeader(mymodel.pixel_size);
        }
    }

    if (Xsize(Isolvent()) != ref_box_size) {
        need_new_mask = true;

        if (verb > 0) {
            std::cerr << " + WARNING: solvent mask box size: " << Xsize(Isolvent())
            << " is not the same as the reference box size: " << ref_box_size << std::endl;
            std::cerr << " + WARNING: re-windowing the mask... " << std::endl;
        }

        if (rank == 0) {
            // Only the leader writes out the new mask
            Isolvent() = Isolvent().setXmippOrigin().windowed(
                Xmipp::init(ref_box_size), Xmipp::last(ref_box_size),
                Xmipp::init(ref_box_size), Xmipp::last(ref_box_size),
                Xmipp::init(ref_box_size), Xmipp::last(ref_box_size)
            ).setXmippOrigin();
        }
    }

    const auto range = minmax(Isolvent());
    if (range.first < 0.0 || range.second > 1.0) {
        need_new_mask = true;

        if (verb > 0) {
            std::cerr << " + WARNING: solvent mask minimum: " << range.first
            << " or maximum: " << range.second << " are outside the [0,1] range." << std::endl;
            std::cerr << " + WARNING: thresholding the mask value to [0,1] range ... " << std::endl;
        }

        for (auto &x : Isolvent())
            if (x < 0.0) { x = 0.0; } else if (x > 1.0) { x = 1.0; }
    }

    if (need_new_mask) {
        // everyone should know about the new mask
        _fn_mask = fn_out + "_solvent" + (solvent_nr == 1 ? "_mask" : integerToString(solvent_nr)) + ".mrc";
        if (rank == 0) Isolvent.write(_fn_mask);
    }

    return;
}

void MlOptimiser::initialiseGeneral(int rank) {

    #ifdef DEBUG
    std::cerr << "Entering initialiseGeneral" << std::endl;
    #endif

    #ifdef TIMING
    // DIFFF = timer.setNew("difff");
    TIMING_EXP    = timer.setNew("expectation");
    TIMING_EXP_1  = timer.setNew("expectation_1");
    TIMING_EXP_1a = timer.setNew("expectation_1a");
    TIMING_EXP_2  = timer.setNew("expectation_2");
    TIMING_EXP_3  = timer.setNew("expectation_3");
    TIMING_EXP_4  = timer.setNew("expectation_4");
    TIMING_EXP_4a = timer.setNew("expectation_4a");
    TIMING_EXP_4b = timer.setNew("expectation_4b");
    TIMING_EXP_4c = timer.setNew("expectation_4c");
    TIMING_EXP_4d = timer.setNew("expectation_4d");
    TIMING_EXP_5  = timer.setNew("expectation_5");
    TIMING_EXP_6  = timer.setNew("expectation_6");
    TIMING_EXP_7  = timer.setNew("expectation_7");
    TIMING_EXP_8  = timer.setNew("expectation_8");
    TIMING_EXP_9  = timer.setNew("expectation_9");

    TIMING_EXP_METADATA =  timer.setNew(" - EXP: metadata shuffling");
    TIMING_EXP_CHANGES =   timer.setNew(" - EXP: monitor changes hidden variables");
    TIMING_MAX =           timer.setNew("maximization");
    TIMING_SOLVFLAT =      timer.setNew("flatten solvent");
    TIMING_UPDATERES =     timer.setNew("update resolution");
    TIMING_RECONS =        timer.setNew("reconstruction");
    TIMING_ESP =           timer.setNew("expectationSomeParticles");
    TIMING_ESP_THR =       timer.setNew("doThreadExpectationSomeParticles");
    TIMING_ESP_ONEPART =   timer.setNew("expectationOneParticle (thr0)");
    TIMING_ESP_ONEPARTN =  timer.setNew("expectationOneParticle (thrN)");
    TIMING_ESP_INI=        timer.setNew(" - EOP: initialise memory");
    TIMING_ESP_FT =        timer.setNew(" - EOP: getFourierTransforms");
    TIMING_ESP_PREC1 =     timer.setNew(" - EOP: precalcShifts1");
    TIMING_ESP_PREC2 =     timer.setNew(" - EOP: precalcShifts2");
    TIMING_ESP_DIFF1 =     timer.setNew(" - EOP: getAllSquaredDifferences1");
    TIMING_ESP_DIFF2 =     timer.setNew(" - EOP: getAllSquaredDifferences2");
    TIMING_ESP_DIFF2_A =   timer.setNew(" - EOP: getD2: A");
    TIMING_ESP_DIFF2_B =   timer.setNew(" - EOP: getD2: B");
    TIMING_ESP_DIFF2_C =   timer.setNew(" - EOP: getD2: C");
    TIMING_ESP_DIFF2_D =   timer.setNew(" - EOP: getD2: D");
    TIMING_ESP_DIFF2_E =   timer.setNew(" - EOP: getD2: E");
    TIMING_DIFF_PROJ =     timer.setNew(" -  - EOPdiff2: project");
    TIMING_DIFF_SHIFT =    timer.setNew(" -  - EOPdiff2: shift");
    TIMING_DIFF2_GETSHIFT =timer.setNew(" -  - EOPdiff2: get shifted img");
    TIMING_DIFF_DIFF2 =    timer.setNew(" -  - EOPdiff2: diff2");
    TIMING_ESP_WEIGHT1 =   timer.setNew(" - EOP: convertDiff2ToWeights1");
    TIMING_ESP_WEIGHT2 =   timer.setNew(" - EOP: convertDiff2ToWeights2");
    TIMING_WEIGHT_EXP =    timer.setNew(" -  - EOPweight: exp");
    TIMING_WEIGHT_SORT =   timer.setNew(" -  - EOPweight: sort");
    TIMING_ESP_WSUM =      timer.setNew(" - EOP: storeWeightedSums");
    TIMING_ESP_PRECW =     timer.setNew(" -  - EOPwsum: precalcShiftsW");
    TIMING_WSUM_PROJ =     timer.setNew(" -  - EOPwsum: project");
    TIMING_WSUM_GETSHIFT = timer.setNew(" -  - EOPwsum: get shifted img");
    TIMING_WSUM_DIFF2 =    timer.setNew(" -  - EOPwsum: diff2");
    TIMING_WSUM_LOCALSUMS =timer.setNew(" -  - EOPwsum: localsums");
    TIMING_WSUM_SUMSHIFT = timer.setNew(" -  - EOPwsum: shiftimg");
    TIMING_WSUM_BACKPROJ = timer.setNew(" -  - EOPwsum: backproject");

    TIMING_ITER_HELICALREFINE = timer.setNew("iterate:  helicalRefinement");
    TIMING_ITER_WRITE         = timer.setNew("iterate:  writeOutput");
    TIMING_ITER_LOCALSYM      = timer.setNew("iterate:  ApplyLocalSymmetry");

    TIMING_EXTRA1= timer.setNew(" -extra1");
    TIMING_EXTRA2= timer.setNew(" -extra2");
    TIMING_EXTRA3= timer.setNew(" -extra3");

    RCT_1 = timer.setNew(" RcT1_BPrefRecon ");
    RCT_2 = timer.setNew(" RcT2_BroadCast ");
    RCT_3 = timer.setNew(" RcT3_maxiOthers ");
    RCT_4 = timer.setNew(" RcT4_monitorHidden ");
    RCT_5 = timer.setNew(" RcT5_sums ");
    RCT_6 = timer.setNew(" RcT6_updatePdf ");
    RCT_7 = timer.setNew(" RcT7_updateNoise ");
    RCT_8 = timer.setNew(" RcT8_initials ");
    #endif

    // Check for errors in the command-line option
    if (parser.checkForErrors(verb))
        REPORT_ERROR("Errors encountered on the command line (see above), exiting...");

        #ifdef RELION_SINGLE_PRECISION
        if (verb > 0)
            std::cout << " Running CPU instructions in single precision. Runs might not be exactly reproducible." << std::endl;
        #else
        if (verb > 0)
            std::cout << " Running CPU instructions in double precision. " << std::endl;
        #endif

    // print symmetry operators or metadata labels before doing anything else...
    if (do_print_symmetry_ops) {
        if (verb > 0) {
            SymList SL;
            SL.writeDefinition(std::cout, sampling.symmetryGroup());
        }
        exit(0);
    }

    if (do_print_metadata_labels) {
        if (verb > 0)
            EMDL::printDefinitions(std::cout);
        exit(0);
    }

    if (fn_data.empty() || fn_out.empty()) {
        REPORT_ERROR("ERROR: provide both --i and --o arguments");
    }

    // For safeguarding the gold-standard separation
    my_halfset = (debug_split_random_half > 0) ? debug_split_random_half : -1;

    // Check if output directory exists
    FileName fn_dir = fn_out.beforeLastOf("/");
    if (!exists(fn_dir))
        REPORT_ERROR("ERROR: output directory does not exist!");

    // Just die if trying to use accelerators and skipping alignments
    if (do_skip_align && (do_gpu || do_cpu))
        REPORT_ERROR("ERROR: you cannot use accelerators when skipping alignments");

    if (do_always_cc)
        do_calculate_initial_sigma_noise = false;


    if (do_shifts_onthefly && (do_gpu || do_cpu)) {
        std::cerr << "WARNING: --onthefly_shifts cannot be combined with --cpu or --gpu, setting do_shifts_onthefly to false" << std::endl;
        do_shifts_onthefly = false;
    }

    // If we are not continuing an old run, now read in the data and the reference images
    if (iter == 0) {
        // Read in the experimental image metadata
        // If do_preread_images: only the leader reads all images into RAM
        bool do_preread = do_preread_images ? do_parallel_disc_io || rank == 0 : false;
        bool is_helical_segment = do_helical_refine || mymodel.ref_dim == 2 && helical_tube_outer_diameter > 0.0;
        int myverb = rank == 0;
        mydata.read(fn_data, true, false, do_preread, is_helical_segment, myverb); // true means ignore original particle name

        // Read in the reference(s) and initialise mymodel
        int refdim = fn_ref == "denovo" ? 3 : 2;
        mymodel.initialiseFromImages(fn_ref, is_3d_model, mydata,
                do_average_unaligned, do_generate_seeds, refs_are_ctf_corrected, ref_angpix, do_sgd, do_trust_ref_size, (rank==0));

    }

    if (mymodel.nr_classes > 1 && do_split_random_halves)
        REPORT_ERROR("ERROR: One cannot use --split_random_halves with more than 1 reference... You could first classify, and then refine each class separately using --random_halves.");

    if (do_join_random_halves && !do_split_random_halves)
        REPORT_ERROR("ERROR: cannot join random halves because they were not split in the previous run");

    // Check all images have the same image_size, otherwise disable non-parallel disc I/O
    if (!mydata.obsModel.allBoxSizesIdentical() && !do_parallel_disc_io)
        REPORT_ERROR("ERROR: non-parallel disc I/O is not implemented when multiple different box sizes are present in the data set. Sorry....");

    // Local symmetry operators
    fn_local_symmetry_masks.clear();
    fn_local_symmetry_operators.clear();
    if (fn_local_symmetry != "None")
        readRelionFormatMasksAndOperators(fn_local_symmetry, fn_local_symmetry_masks, fn_local_symmetry_operators, mymodel.pixel_size, false);

    // For multi-body refinement: Read in the masks in the input STAR file, add a soft-edge to them, and write them to disc with the standard name
    if (do_initialise_bodies) {

        if (verb > 0)
            std::cout << " + Initialising multi-body refinement ..." << std::endl;

        if (mymodel.nr_classes > 1)
            REPORT_ERROR("ERROR: One cannot use multiple classes with multi-body refinement!");
        if (do_sgd)
            REPORT_ERROR("ERROR: One cannot use SGD with multi-body refinement!");
        if (do_helical_refine)
            REPORT_ERROR("ERROR: One cannot use helical symmetry with multi-body refinement!");
        if (!do_split_random_halves)
            REPORT_ERROR("ERROR: One has to use split random halves with multi-body refinement!");

        // This reads the masks, calculates com_bodies and orient_bodies
        mymodel.initialiseBodies(fn_body_masks, fn_out, true, rank);
        mymodel.writeBildFileBodies(fn_out + "_bodies.bild");

        // For multi-body refinement: expand the MetaDataTables with orientations for all bodies
        mydata.initialiseBodies(mymodel.nr_bodies);

        nr_iter_wo_resol_gain = 0;
        nr_iter_wo_large_hidden_variable_changes = 0;
        current_changes_optimal_classes = 9999999;
        current_changes_optimal_offsets = 999.0;
        current_changes_optimal_orientations = 999.0;
        smallest_changes_optimal_classes = 9999999;
        smallest_changes_optimal_offsets = 999.0;
        smallest_changes_optimal_orientations = 999.0;

        if (autosampling_hporder_local_searches > sampling.healpix_order) {
            mymodel.sigma2_rot = mymodel.sigma2_tilt = mymodel.sigma2_psi = 0.0;
        }

        // We're not using pdf_directions anyway (only in wsum for bild file), so reset the sizes of the vectors just in case
        directions_have_changed = true;

        // TMP?
        //do_norm_correction = false;

        // don't perturb angles anymore
        //sampling.perturbation_factor = 0.0;
        // Only do a single pass through the local-search orientations
        //adaptive_oversampling = 0;

        // Start at iteration 1 again
        iter = 0;

    } else if (fn_body_masks.empty()) {
        mymodel.nr_bodies = 1;
    }

    if (mymodel.nr_bodies > 1) {
        // This creates a rotation matrix for (rot,tilt,psi) = (0,90,0)
        // It will be used to make all Abody orientation matrices relative to (0,90,0) instead of the more logical (0,0,0)
        // This is useful, as psi-priors are ill-defined around tilt=0, as rot becomes the same as -psi!!
        A_rot90 = rotation3DMatrix(-90.0, 'Y', false);
        A_rot90T = A_rot90.transpose();
    }

    if (fn_fourier_mask != "None") {
        // Used also for continuations...
        Image<RFLOAT> Itmp;
        Itmp.read(fn_fourier_mask);
        helical_fourier_mask = Itmp();
    }

    // 09 Jun 2015 - Shaoda, Helical refinement
    if (do_helical_refine) {
        if (fn_fourier_mask == "None" && !helical_fourier_mask_resols.empty()) {
            std::vector<std::string> resols = split(helical_fourier_mask_resols, ",");
            if (resols.size() % 2 == 1) REPORT_ERROR("Provide an even number of start-end resolutions for --fourier_exclude_resols");
            std::vector<RFLOAT> resols_end, resols_start;
            for (int nshell = 0; nshell < resols.size() / 2; nshell++) {
                resols_start.push_back(textToFloat(resols[2 * nshell    ]));
                resols_end  .push_back(textToFloat(resols[2 * nshell + 1]));
            }


            MultidimArray<RFLOAT> tmpmsk = generateBinaryHelicalFourierMask(
                Xsize(mymodel.Iref[0]), Ysize(mymodel.Iref[0]), Zsize(mymodel.Iref[0]) / 2 + 1,
                resols_start, resols_end, mymodel.pixel_size
            );
            tmpmsk.setXmippOrigin();
            // make a 2-pixel soft edge of fourier mask
            autoMask(tmpmsk, helical_fourier_mask, 0.5, 0.0, 2.0, false, nr_threads);

            // Save the mask to the output directory
            fn_fourier_mask = fn_out + "fourier_mask.mrc";
            if (verb > 0) {
                Image<RFLOAT> Itmp;
                Itmp() = helical_fourier_mask;
                // use this also for continuation
                Itmp.write(fn_fourier_mask);
            }

        }


        if (mymodel.nr_bodies != 1)
            REPORT_ERROR("ERROR: cannot do multi-body refinement for helices!");

        if (do_shifts_onthefly && !ignore_helical_symmetry && verb > 0) {
            std::cerr << " WARNING: On-the-fly shifts slow down helical reconstructions with CPUs considerably. "
                    << "Enable this option only if limited RAM causes trouble (e.g. too large segment boxes used or in 3D sub-tomogram averaging). "
                    << std::endl;
        }

        // Set particle diameter to 90% the box size if user does not give this parameter
        if (particle_diameter < 0.0) {
            if (verb > 0)
                std::cout << " Automatically set particle diameter to 90% the box size for 3D helical reconstruction." << std::endl;

            particle_diameter = (RFLOAT) mymodel.ori_size;
            if ((RFLOAT) mymodel.ori_size * 0.05 < width_mask_edge)
                particle_diameter -= 2.0 * width_mask_edge;
            particle_diameter *= 0.90 * mymodel.pixel_size;
        }

        if (ignore_helical_symmetry) { mymodel.helical_nr_asu = 1; }

        if (!do_helical_symmetry_local_refinement || ignore_helical_symmetry) {
            mymodel.helical_twist_min = mymodel.helical_twist_max = helical_twist_initial;
            mymodel.helical_rise_min = mymodel.helical_rise_max = helical_rise_initial;
        }

        if (mymodel.ref_dim == 3) {
            checkParametersFor3DHelicalReconstruction(
                ignore_helical_symmetry,
                do_helical_symmetry_local_refinement,
                mymodel.helical_nr_asu,
                helical_rise_initial,
                mymodel.helical_rise_min, mymodel.helical_rise_max,
                helical_twist_initial,
                mymodel.helical_twist_min, mymodel.helical_twist_max,
                mymodel.ori_size,
                mymodel.pixel_size,
                helical_z_percentage,
                particle_diameter,
                helical_tube_inner_diameter, helical_tube_outer_diameter
            );
        }
    } else {
        if (do_helical_symmetry_local_refinement)
            REPORT_ERROR("ERROR: cannot do local refinement of helical parameters for non-helical segments!");
    }
    if (mymodel.ref_dim == 2 && helical_tube_outer_diameter > 0.0 && particle_diameter < helical_tube_outer_diameter)
        REPORT_ERROR("ERROR: 2D classification: Helical tube diameter should be smaller than particle diameter!");

    if (do_always_join_random_halves)
        std::cout << " Joining half-reconstructions at each iteration: this is a developmental option to test sub-optimal FSC usage only! " << std::endl;

    // If fn_tau is provided, read in the tau spectrum
    fix_tau = false;
    if (fn_tau != "None") {
        fix_tau = true;
        mymodel.readTauSpectrum(fn_tau, verb);
    }

    if (do_auto_refine) {
        nr_iter = nr_iter_max;
        has_fine_enough_angular_sampling = false;
        has_converged = false;

        if (mymodel.tau2_fudge_factor > 1. && verb > 0) {
            std::cerr << " WARNING: Using tau2_fudge of " <<mymodel.tau2_fudge_factor << ", which will lead to inflated resolution estimates during refinement." << std::endl;
            std::cerr << " WARNING: This option will most likely lead to overfitting and therefore difficult to interpret maps: proceed with caution and know the risks!" << std::endl;
            std::cerr << " WARNING: You have to run postprocessing afterwards to get a reliable, gold-standard resolution estimate! " << std::endl;
        }

        if (iter == 0 && sampling.healpix_order >= autosampling_hporder_local_searches) {
            mymodel.orientational_prior_mode = PRIOR_ROTTILT_PSI;
            sampling.is_3D = (mymodel.ref_dim == 3);
            RFLOAT rottilt_step = sampling.getAngularSampling(adaptive_oversampling);
            mymodel.sigma2_rot = mymodel.sigma2_tilt = mymodel.sigma2_psi = 2. * 2. * rottilt_step * rottilt_step;
        }

        // Check whether we had converged already
        // Jobs often fail in the last iteration, these lines below allow restarting of only the last iteration
        if (do_auto_refine && iter > 0) {
            if (do_force_converge) {
                has_converged = true;
                do_join_random_halves = true;
                // In the last iteration, include all data until Nyquist
                do_use_all_data = true;
            }
        }

    }

    // Initialise the sampling object (sets prior mode and fills translations and rotations inside sampling object)
    // May06,2015 - Shaoda & Sjors, initialise for helical translations
    bool do_local_searches_helical = do_auto_refine && do_helical_refine && sampling.healpix_order >= autosampling_hporder_local_searches;

    if (iter == 0) {
        // Sjors 26Jul2019: the sampling.offset_step and range are in Angstroms, but inputs from command line are in pixels!
        // For continuation runs, this is done in parseContinue, but for new refinements this still needs to be done.
        // Use the pixel size from the model for this!
        sampling.offset_range *= mymodel.pixel_size;
        sampling.offset_step *= mymodel.pixel_size;
    }
    sampling.initialise(
        mymodel.ref_dim, mymodel.data_dim == 3, do_gpu, verb > 0,
        do_local_searches_helical, do_helical_refine && !ignore_helical_symmetry,
        helical_rise_initial, helical_twist_initial
    );

    // Now that sampling is initialised, also modify sigma2_rot for the helical refinement
    if (do_local_searches_helical && !ignore_helical_symmetry && iter == 0) {
        // Aug20,2015 - Shaoda, Helical refinement
        RFLOAT rottilt_step = sampling.getAngularSampling(adaptive_oversampling);
        mymodel.sigma2_rot = getHelicalSigma2Rot(helical_rise_initial, helical_twist_initial, sampling.helical_offset_step, rottilt_step, mymodel.sigma2_rot);
    }

    if (particle_diameter < 0.0)
        particle_diameter = (mymodel.ori_size - width_mask_edge) * mymodel.pixel_size;

    // For do_average_unaligned, always use initial low_pass filter
    if (do_average_unaligned && ini_high < 0.0) {
        // By default, use 0.07 dig.freq. low-pass filter
        // See S.H.W. Scheres (2010) Meth Enzym.
        ini_high = 1.0 / mymodel.getResolution(round(0.07 * mymodel.ori_size));
    }

    // For skipped alignments
    // Also do not perturb this orientation, nor do oversampling or priors
    // Switch off on-the-fly shifts, as that will not work when skipping alignments! (it isn't necessary anyway in that case)
    if (do_skip_align || do_skip_rotate) {
        mymodel.orientational_prior_mode = NOPRIOR;
        adaptive_oversampling = 0;
        sampling.perturbation_factor = 0.0;
        sampling.random_perturbation = 0.0;
        sampling.addOneOrientation(0.,0.,0., true);
        directions_have_changed = true;
    } else if (do_only_sample_tilt) {
        std::cout << " Only sampling tilt, keep rot and psi fixed." << std::endl;

        mymodel.orientational_prior_mode = NOPRIOR;
        adaptive_oversampling = 0;
        sampling.perturbation_factor = 0.0;
        sampling.random_perturbation = 0.0;

        sampling.directions_ipix.clear();
        sampling.rot_angles.clear();
        sampling.tilt_angles.clear();
        RFLOAT rot = 0.0, psi = 0.0, tilt;
        int ipix = 0;
        for (tilt = -180.0; tilt < 180.0; tilt+= sampling.getAngularSampling(), ipix++) {
            sampling.rot_angles.push_back(rot);
            sampling.tilt_angles.push_back(tilt);
            sampling.directions_ipix.push_back(ipix);
        }
        sampling.psi_angles.clear();
        sampling.psi_angles.push_back(psi);
        directions_have_changed = true;
    }
    if (do_skip_align) {
        RFLOAT dummy = 0.0;
        sampling.addOneTranslation(dummy, dummy, dummy, true);
        do_shifts_onthefly = false; // on-the-fly shifts are incompatible with do_skip_align!
    }
    if (do_bimodal_psi && mymodel.sigma2_psi > 0.0 && verb > 0)
        std::cout << " Using bimodal priors on the psi angle..." << std::endl;

    // Resize the pdf_direction arrays to the correct size and fill with an even distribution
    if (directions_have_changed)
        mymodel.initialisePdfDirection(sampling.NrDirections());

    // Initialise the wsum_model according to the mymodel
    wsum_model.initialise(mymodel, sampling.symmetryGroup(), asymmetric_padding, skip_gridding);

    // Initialise sums of hidden variable changes
    // In later iterations, this will be done in updateOverallChangesInHiddenVariables
    sum_changes_optimal_orientations = 0.0;
    sum_changes_optimal_offsets = 0.0;
    sum_changes_optimal_classes = 0.0;
    sum_changes_count = 0.0;

    if (mymodel.data_dim == 3) {

        // TODO: later do norm correction?!
        // Don't do norm correction for volume averaging at this stage....
        do_norm_correction = false;

        if ((!do_helical_refine || ignore_helical_symmetry) && !do_cpu && !do_gpu) // For 3D helical sub-tomogram averaging, either is OK, so let the user decide
            do_shifts_onthefly = true; // save RAM for volume data (storing all shifted versions would take a lot!)

        if (do_skip_align)
            do_shifts_onthefly = false; // on-the-fly shifts are incompatible with do_skip_align!
        // getMetaAndImageData is not made for passing multiple volumes!
        do_parallel_disc_io = true;
    }

    // Tabulated sine and cosine values (for 2D helical segments / 3D helical sub-tomogram averaging with on-the-fly shifts)
    if (do_shifts_onthefly && do_helical_refine && !ignore_helical_symmetry) {
        tab_sin.initialise(100000);
        tab_cos.initialise(100000);
    }

    // Skip scale correction if there are no groups
    if (mymodel.nr_groups == 1) { do_scale_correction = false; }

    // Check for rlnReconstructImageName in the data.star file. If it is present, set do_use_reconstruct_images to true
    do_use_reconstruct_images = mydata.MDimg.containsLabel(EMDL::IMAGE_RECONSTRUCT_NAME);
    if (do_use_reconstruct_images && verb > 0)
        std::cout <<" Using rlnReconstructImageName from the input data.star file!" << std::endl;

    // For new thread-parallelization: each thread does 1 particle, so nr_pool=nr_threads
    nr_pool = x_pool * nr_threads;

    if (do_sgd) {
        #ifndef ALLOW_CTF_IN_SGD
        do_ctf_correction = false;
        if (verb > 0) {
            std::cerr << " + Skipping CTF-modulation in SGD, as mentioned in Claim 1 of patent US10,282,513B2." << std::endl;
            std::cerr << " + Note that the output map will not be CTF-corrrected, and this should be specified for subsequent refinement." << std::endl;
        }
        #endif

        sgd_inires_pix = mymodel.getPixelFromResolution(1.0 / sgd_ini_resol);
        sgd_finres_pix = mymodel.getPixelFromResolution(1.0 / sgd_fin_resol);
        // for continuation jobs (iter>0): could do some more iterations as specified by nr_iter
        nr_iter = sgd_ini_iter + sgd_fin_iter + sgd_inbetween_iter;
    } else {
        subset_size = -1;
        mu = 0.0;
    }

    if (do_fast_subsets) {
        if (nr_iter < 20)
            REPORT_ERROR("ERROR: when using --fast_subsets you have to perform at least 20 iterations!");
        if (do_auto_refine)
            REPORT_ERROR("ERROR: you cannot use --fast_subsets together with --auto_refine");
    }

    // Check mask angpix, boxsize and [0,1] compliance right away.
    if (fn_mask != "None") checkMask(fn_mask, 1, rank);
    if (fn_mask2 != "None") checkMask(fn_mask2, 2, rank);

    // Write out unmasked 2D class averages
    do_write_unmasked_refs = (mymodel.ref_dim == 2);

    #ifdef DEBUG
    std::cerr << "Leaving initialiseGeneral" << std::endl;
    #endif

}

void MlOptimiser::initialiseWorkLoad() {

    // Note, this function is overloaded in ml_optimiser_mpi (where random_seed is only set by the leader and then send to all followers!)

    // Randomise the order of the particles
    if (random_seed == -1) { random_seed = time(NULL); }
    // Also randomize random-number-generator for perturbations on the angles
    init_random_generator(random_seed);

    if (do_split_random_halves && debug_split_random_half > 0) {

        // Split the data into two random halves
        mydata.divideParticlesInRandomHalves(random_seed, do_helical_refine);
        // rank=0 will work on subset 2, because rank%2==0
        my_halfset = debug_split_random_half;
    }

    divide_equally(mydata.numberOfParticles(), 1, 0, my_first_particle_id, my_last_particle_id);

    // Now copy particle stacks to scratch if needed
    if (!fn_scratch.empty() && !do_preread_images) {
        mydata.setScratchDirectory(fn_scratch, do_reuse_scratch, 1);

        if (!do_reuse_scratch) {
            mydata.prepareScratchDirectory(fn_scratch);
            bool also_do_ctfimage = (mymodel.data_dim == 3 && do_ctf_correction);
            mydata.copyParticlesToScratch(1, true, also_do_ctfimage, keep_free_scratch_Gb);
        }
    }
}

void MlOptimiser::calculateSumOfPowerSpectraAndAverageImage(MultidimArray<RFLOAT> &Mavg, bool myverb) {

    #ifdef DEBUG_INI
    std::cerr<<"MlOptimiser::calculateSumOfPowerSpectraAndAverageImage Entering"<<std::endl;
    #endif

    int barstep, my_nr_particles = my_last_particle_id - my_first_particle_id + 1;

    // Initialise Mavg
    if (mydata.is_3D) {
        Mavg.initZeros(mymodel.ori_size, mymodel.ori_size, mymodel.ori_size);
    } else {
        Mavg.initZeros(mymodel.ori_size, mymodel.ori_size);
    }
    Mavg.setXmippOrigin();

    if (my_nr_particles < 1) {
        // Master doesn't do anything here...
        // But still set Mavg the right size for AllReduce later on
        return;
    }

    if (myverb > 0) {
        std::cout << " Estimating initial noise spectra " << std::endl;
        init_progress_bar(my_nr_particles);
        barstep = std::max(1, my_nr_particles / 60);
    }

    // Only open stacks once and then read multiple images
    fImageHandler hFile;

    // Note the loop over the particles (part_id) is MPI-parallelized
    int nr_particles_done = 0;
    FileName fn_img;
    // For spectrum calculation: recycle the transformer (so do not call getSpectrum all the time)
    FourierTransformer transformer;
    MetaDataTable MDimg;

    // Start reconstructions at ini_high or 0.07 digital frequencies....
    wsum_model.current_size = ini_high <= 0.0 ?
        1.0 / mymodel.getResolution(round(0.07 * mymodel.ori_size)) :
        mymodel.getPixelFromResolution(1.0 / ini_high);
    wsum_model.initZeros();

    for (long int part_id_sorted = my_first_particle_id; part_id_sorted <= my_last_particle_id; part_id_sorted++, nr_particles_done++) {
        const long int part_id = mydata.sorted_idx[part_id_sorted];
        for (int img_id = 0; img_id < mydata.numberOfImagesInParticle(part_id); img_id++) {
            long int group_id = mydata.getGroupId(part_id, img_id);
            int optics_group = mydata.getOpticsGroup(part_id, img_id);
            RFLOAT my_pixel_size = mydata.getOpticsPixelSize(optics_group);
            int my_image_size = mydata.getOpticsImageSize(optics_group);

            // Extract the relevant MetaDataTable row from MDimg
            MDimg = mydata.getMetaDataImage(part_id, img_id);

            // Read image from disc
            Image<RFLOAT> img;
            if (do_preread_images && do_parallel_disc_io) {
                img() = mydata.particles[part_id].images[img_id].img;
            } else {
                try {
                    fn_img = mydata.getImageNameOnScratch(part_id, img_id);
                    if (!do_parallel_disc_io) {
                        // When not doing parallel disk IO,
                        // only those MPI processes running on the same node as the leader have scratch.
                        FileName fn_stack;
                        long int dump;
                        fn_img.decompose(dump, fn_stack);
                        if (!exists(fn_stack))
                            fn_img = MDimg.getValue<std::string>(EMDL::IMAGE_NAME, MDimg.size() - 1);
                    }
                } catch (const char *errmsg) {
                    fn_img = MDimg.getValue<std::string>(EMDL::IMAGE_NAME, MDimg.size() - 1);
                }

                FileName fn_stack;
                long int dump;
                fn_img.decompose(dump, fn_stack);
                FileName fn_open_stack = "";
                if (fn_stack != fn_open_stack) {
                    hFile.openFile(fn_stack, WRITE_READONLY);
                    fn_open_stack = fn_stack;
                }
                img.readFromOpenFile(fn_img, hFile, -1, false);
                img().setXmippOrigin();
            }

            // 24 May 2015 - Shaoda & Sjors, Helical refinement
            RFLOAT psi_deg = 0.0, tilt_deg = 0.0;
            bool is_helical_segment = do_helical_refine || mymodel.ref_dim == 2 && helical_tube_outer_diameter > 0.0;
            const long int i = MDimg.size() - 1;
            if (is_helical_segment) {
                try {
                    psi_deg = MDimg.getValue<RFLOAT>(EMDL::ORIENT_PSI_PRIOR, i);
                } catch (const char* errmsg) { try {
                    psi_deg = MDimg.getValue<RFLOAT>(EMDL::ORIENT_PSI, i);
                } catch (const char* errmsg) {
                    REPORT_ERROR("ml_optimiser.cpp::calculateSumOfPowerSpectraAndAverageImage: Psi priors of helical segments are missing!");
                } }

                try {
                    tilt_deg = MDimg.getValue<RFLOAT>(EMDL::ORIENT_TILT_PRIOR, i);
                } catch (const char* errmsg) { try {
                    tilt_deg = MDimg.getValue<RFLOAT>(EMDL::ORIENT_TILT, i);
                } catch (const char* errmsg) {
                    REPORT_ERROR("ml_optimiser.cpp::calculateSumOfPowerSpectraAndAverageImage: Tilt priors of helical segments are missing!");
                } }
            }

            // Check that the average in the noise area is approximately zero and the stddev is one
            if (!dont_raise_norm_error && verb > 0) {
                // NEW METHOD
                RFLOAT sphere_radius_pix, cyl_radius_pix;
                cyl_radius_pix = helical_tube_outer_diameter / (2.0 * my_pixel_size);
                sphere_radius_pix = particle_diameter / (2.0 * my_pixel_size);
                const auto bg_stats = calculateBackgroundAvgStddev(img, round(sphere_radius_pix), is_helical_segment, cyl_radius_pix, tilt_deg, psi_deg);

                // Average should be close to 0, i.e. max +/-50% of stddev...
                // Stddev should be close to 1, i.e. larger than 0.5 and smaller than 2)
                if (abs(bg_stats.avg / bg_stats.stddev) > 0.5 || bg_stats.stddev < 0.5 || bg_stats.stddev > 2.0) {
                    std::cerr << " fn_img= " << fn_img << " bg_avg= " << bg_stats.avg << " bg_stddev= " << bg_stats.stddev << std::flush;
                    if (is_helical_segment) {
                        std::cerr << " tube_bg_radius= " << cyl_radius_pix << " psi_deg= " << psi_deg << " tilt_deg= " << tilt_deg << " (this is a particle from a helix)" << std::flush;
                    } else {
                        std::cerr << " bg_radius= " << sphere_radius_pix << std::flush;
                    }
                    std::cerr << std::endl;
                    std::cerr << "WARNING: It appears that these images have not been normalised to an average background value of 0 and a stddev value of 1. \n \
                            Note that the average and stddev values for the background are calculated: \n \
                            (1) for single particles: outside a circle with the particle diameter \n \
                            (2) for helical segments: outside a cylinder (tube) with the helical tube diameter \n \
                            You can use the relion_preprocess program to normalise your images \n \
                            If you are sure you have normalised the images correctly (also see the RELION Wiki), you can switch off this warning message using the --dont_check_norm command line option" <<std::endl;
                    dont_raise_norm_error = true;
                }
            }


            // Apply a similar softMask as below (assume zero translations)
            if (do_zero_mask) {
                // 24 May 2015 - Shaoda & Sjors, Helical refinement
                if (is_helical_segment) {
                    softMaskOutsideMapForHelix(
                        img(), psi_deg, tilt_deg,
                        particle_diameter / (2.0 * my_pixel_size),
                        helical_tube_outer_diameter / (2.0 * my_pixel_size),
                        width_mask_edge
                    );
                } else {
                    softMaskOutsideMap(img(), particle_diameter / (2.0 * my_pixel_size), width_mask_edge);
                }
            }

            // Keep track of the average image (only to correct power spectra, no longer for initial references!)

            // Rescale img() onto Mavg, as optics_groups may have different box sizes and pixel sizes...
            // a) rescale to same pixel size
            if (fabs(my_pixel_size - mymodel.pixel_size) > 0.0001) {
                int rescalesize = round(Xsize(img()) * my_pixel_size / mymodel.pixel_size);
                // Enforce divisibility by 2
                rescalesize += rescalesize % 2;
                resizeMap(img(), rescalesize);
            }
            // b) window to same box size
            img().setXmippOrigin();
            if (fabs(Xsize(img()) - mymodel.ori_size) > 0) {
                switch (mymodel.data_dim) {
                    case 2:
                    img() = img().windowed(
                        Xmipp::init(mymodel.ori_size), Xmipp::last(mymodel.ori_size),
                        Xmipp::init(mymodel.ori_size), Xmipp::last(mymodel.ori_size)
                    );
                    break;
                    case 3:
                    img() = img().windowed(
                        Xmipp::init(mymodel.ori_size), Xmipp::last(mymodel.ori_size),
                        Xmipp::init(mymodel.ori_size), Xmipp::last(mymodel.ori_size),
                        Xmipp::init(mymodel.ori_size), Xmipp::last(mymodel.ori_size)
                    );
                    break;
                }
            }
            Mavg += img();

            // Calculate the power spectrum of this particle
            int spectral_size = mymodel.ori_size / 2 + 1;
            MultidimArray<RFLOAT> ind_spectrum = MultidimArray<RFLOAT>::zeros(spectral_size);
            MultidimArray<RFLOAT> count        = MultidimArray<RFLOAT>::zeros(spectral_size);
            // recycle the same transformer for all images
            MultidimArray<Complex> &Faux = transformer.FourierTransform(img());

            FOR_ALL_ELEMENTS_IN_FFTW_TRANSFORM(Faux) {
                long int idx = round(euclid(ip, jp, kp));
                if (idx < spectral_size) {
                    ind_spectrum.elem(idx) += norm(direct::elem(Faux, i, j, k));
                    count.elem(idx) += 1.0;
                }
            }
            ind_spectrum /= count;

            // Resize the power_class spectrum to the correct size and keep sum
            wsum_model.sigma2_noise[group_id] += ind_spectrum;
            wsum_model.sumw_group[group_id] += 1.0;

            // When doing SGD, only take the first sgd_ini_subset_size*mymodel.nr_classes images to calculate the initial reconstruction
            if (fn_ref == "None" && !(do_sgd && part_id < sgd_ini_subset_size * mymodel.nr_classes)) {

                MultidimArray<RFLOAT> Fctf, Fweight;

                // Make sure MPI and sequential behave exactly the same
                init_random_generator(random_seed + part_id);
                // Randomize the initial orientations for initial reference generation at this step....
                // TODO: this is not an even angular distribution....
                RFLOAT rot  = mymodel.ref_dim == 2 ? 0.0 : rnd_unif() * 360.0;
                RFLOAT tilt = mymodel.ref_dim == 2 ? 0.0 : rnd_unif() * 180.0;
                RFLOAT psi  = rnd_unif() * 360.0;
                int iclass  = rnd_unif() * mymodel.nr_classes;
                if (iclass == mymodel.nr_classes)
                    iclass = mymodel.nr_classes - 1;
                if (iclass >= mymodel.nr_classes) {
                    // Should not happen but without this some people get errors in Set2DFourierTransform
                    // TODO: investigate
                    std::cerr << "WARNING: numerical issue in initial class assignment. Your result is NOT compromised but please report this to our issue tracker.\n";
                    std::cerr << "         iclass = " << iclass << " nr_classes = " << mymodel.nr_classes << " sizeof(RFLOAT) = " << sizeof(RFLOAT) << std::endl;
                    iclass = mymodel.nr_classes - 1;
                }
                Matrix2D<RFLOAT> A = Euler::angles2matrix(rot, tilt, psi);

                // At this point anisotropic magnification shouldn't matter
                // Also: dont apply scaleDifference, as img() was rescaled to mymodel.ori_size and mymodel.pixel_size
                // if (mydata.obsModel.hasMagMatrices) { A *= mydata.obsModel.anisoMag(optics_group); }
                //A *= mydata.obsModel.scaleDifference(optics_group, mymodel.ori_size, mymodel.pixel_size);
                // Construct initial references from random subsets
                MultidimArray<Complex> Fimg = windowFourierTransform(Faux, wsum_model.current_size);
                Fctf.resize(CenterFFTbySign(Fimg));
                Fctf = 1.0;

                // Apply CTF if necessary (skip this for subtomograms!)
                if (do_ctf_correction && mymodel.data_dim != 3) {
                    CTF ctf = CtfHelper::makeCTF(MDimg, &mydata.obsModel, 0); // This MDimg only contains one particle!
                    Fctf = CtfHelper::getFftwImage(
                        ctf, Xsize(Fctf), Ysize(Fctf),
                        mymodel.ori_size, mymodel.ori_size, mymodel.pixel_size,
                        &mydata.obsModel,
                        ctf_phase_flipped, only_flip_phases, intact_ctf_first_peak, true, do_ctf_padding
                    );

                    for (long int n = 0; n < Fimg.size(); n++) {
                        Fimg[n] *= Fctf[n];
                        Fctf[n] *= Fctf[n];
                    }
                }
                wsum_model.BPref[iclass].set2DFourierTransform(Fimg, A, &Fctf);
            }
        }

        if (myverb > 0 && nr_particles_done % barstep == 0) {
            progress_bar(nr_particles_done);
            // Abort through the pipeline_control system
            if (pipeline_control_check_abort_job())
                break;
        }
    }


    // Clean up the fftw object completely
    // This is something that needs to be done manually, as among multiple threads only one of them may actually do this
    transformer.cleanup();

    if (myverb > 0)
        progress_bar(my_nr_particles);

    #ifdef DEBUG_INI
    std::cerr << "MlOptimiser::calculateSumOfPowerSpectraAndAverageImage Leaving" << std::endl;
    #endif

}

void MlOptimiser::setSigmaNoiseEstimatesAndSetAverageImage(MultidimArray<RFLOAT> &Mavg) {

    #ifdef DEBUG_INI
    std::cerr << "MlOptimiser::setSigmaNoiseEstimatesAndSetAverageImage Entering" << std::endl;
    #endif

    // First calculate average image
    RFLOAT total_sum = 0.0;
    for (int igroup = 0; igroup < mymodel.nr_groups; igroup++) {
        mymodel.nr_particles_per_group[igroup] = round(wsum_model.sumw_group[igroup]);
        total_sum += wsum_model.sumw_group[igroup];
    }
    Mavg /= total_sum;

    if (fn_ref == "None") {
        for (int iclass = 0; iclass < mymodel.nr_classes * mymodel.nr_bodies; iclass++) {

            MultidimArray<RFLOAT> tau2;
            mymodel.Iref[iclass] = wsum_model.BPref[iclass].reconstruct(gridding_nr_iter, false, tau2);
            // 2D projection data were CTF-corrected, subtomograms were not
            refs_are_ctf_corrected = mymodel.data_dim != 3;
        }
    }

    // Calculate sigma2_noise estimates as average of power class spectra, and subtract power spectrum of the average image from that
    if (do_calculate_initial_sigma_noise) {
        // Calculate power spectrum of the average image
        auto spect = getSpectrum(Mavg, POWER_SPECTRUM) / 2.0;  // because of 2-dimensionality of the complex plane
        spect.resize(mymodel.sigma2_noise[0]);

        // Set noise spectra, once for each group
        for (int igroup = 0; igroup < wsum_model.nr_groups; igroup++) {
            // Factor 2 because of 2-dimensionality of the complex plane
            if (wsum_model.sumw_group[igroup] > 0.0) {
                //std::cerr << " igroup= " << igroup << " wsum_model.sigma2_noise[igroup].sum()= " << wsum_model.sigma2_noise[igroup].sum() << " wsum_model.sumw_group[igroup]= " << wsum_model.sumw_group[igroup] << std::endl;
                mymodel.sigma2_noise[igroup] = wsum_model.sigma2_noise[igroup] / (2.0 * wsum_model.sumw_group[igroup]);

                // Now subtract power spectrum of the average image from the average power spectrum of the individual images
                mymodel.sigma2_noise[igroup] -= spect;

                // Convert any negative sigma2_noise values to nearby positive value
                for (long int n = 0; n < (spect).size(); n++) {
                    RFLOAT curr = mymodel.sigma2_noise[igroup][n];
                    RFLOAT prev = mymodel.sigma2_noise[igroup][n - 1];
                    if (curr < 0.0) {
                        if (n - 1 >= 0 && prev > 0.0) {
                            // Use the previous value (if positive)
                            curr = prev;
                        } else {
                            // Try to find the next positive value
                            for (int nn = n + 1; true; nn++) {
                                if (nn > Xsize(mymodel.sigma2_noise[igroup])) {
                                    std::cerr << " igroup= " << igroup << " n= " << n << " mymodel.sigma2_noise[igroup]= " << mymodel.sigma2_noise[igroup] << std::endl;
                                    REPORT_ERROR("BUG! cannot find positive values in sigma2_noise spectrum");
                                }
                                RFLOAT next = mymodel.sigma2_noise[igroup][nn];
                                if (next > 0.0) {
                                    curr = next;
                                    break;
                                }
                            }
                        }
                    }
                }
            } else {
                // no particles in this group...
                mymodel.sigma2_noise[igroup].initZeros();
            }
        }
    }
    #ifdef DEBUG_INI
        std::cerr<<"MlOptimiser::setSigmaNoiseEstimatesAndSetAverageImage Leaving"<<std::endl;
    #endif
}

void MlOptimiser::initialLowPassFilterReferences() {
    if (ini_high <= 0.0) return;
    // Make a soft (raised cosine) filter in Fourier space to prevent artefacts in real-space
    // The raised cosine goes through 0.5 at the filter frequency and has a width of width_mask_edge fourier pixels
    RFLOAT radius = mymodel.ori_size * mymodel.pixel_size / ini_high;
    radius -= WIDTH_FMASK_EDGE / 2.0;
    RFLOAT radius_p = radius + WIDTH_FMASK_EDGE;
    FourierTransformer transformer;
    for (int iclass = 0; iclass < mymodel.nr_classes; iclass++) {
        auto &reference = mymodel.Iref[iclass];
        MultidimArray<Complex> Faux = transformer.FourierTransform(reference);
        FOR_ALL_ELEMENTS_IN_FFTW_TRANSFORM(Faux) {
            RFLOAT r = euclid(ip, jp, kp);
            if (r < radius) {
                continue;
            } else if (r > radius_p) {
                direct::elem(Faux, i, j, k) = 0.0;
            } else {
                direct::elem(Faux, i, j, k) *= 0.5 - 0.5 * cos(PI * (radius_p - r) / WIDTH_FMASK_EDGE);
            }
        }
        reference = transformer.inverseFourierTransform(Faux);
    }
}

/** ========================== EM-Iteration  ================================= */

void MlOptimiser::iterateSetup() {

    // Make a barrier where all working threads wait
    global_barrier = new Barrier(nr_threads - 1);

    // Create threads to start working
    global_ThreadManager = new ThreadManager(nr_threads, this);

    // Set up the thread task distributors for the particles and the orientations (will be resized later on)
    exp_ipart_ThreadTaskDistributor = new ThreadTaskDistributor(nr_threads, 1);

}
void MlOptimiser::iterateWrapUp() {

    // delete barrier, threads and task distributors
    delete global_barrier;
    delete global_ThreadManager;
    delete exp_ipart_ThreadTaskDistributor;

    // Delete volatile space on scratch
    if (!keep_scratch)
        mydata.deleteDataOnScratch();

    #ifdef MKLFFT
    fftw_cleanup_threads();
    #endif
}

void MlOptimiser::iterate() {

    if (do_split_random_halves && debug_split_random_half == 0)
        REPORT_ERROR("ERROR: Cannot split data into random halves without using MPI! For debugging ONLY, use --debug_split_random_half 1 (or 2)");

    // launch threads etc
    iterateSetup();

    // Update the current resolution and image sizes, and precalculate resolution pointers
    // The rest of the time this will be done after maximization and before writing output files,
    // so that current resolution is in the output files of the current iteration
    updateCurrentResolution();

    /*
    // If we're doing a restart from subsets, then do not increment the iteration number in the restart!
    if (subset > 0) {

        iter--;
        std::cerr << " iter= " << iter << std::endl;
    }
    */

    bool has_already_reached_convergence = false;
    for (iter = iter + 1; iter <= nr_iter; iter++) {
        #ifdef DEBUG_HELICAL_ORIENTATIONAL_SEARCH
        // 18 May 2015 - Shaoda & Sjors, Helical refinement (orientational searches)
        std::cerr << "\n\n"
        << "MlOptimiser::iterate()\n"
        << "iter = " << iter << ", do_helical_refine == "
        << (do_helical_refine ? "true" : "false")
        << std::endl;
        #endif

        {
        ifdefTIMING(TicToc tt (timer, TIMING_EXP);)

        // Update subset_size
        updateSubsetSize();

        // Randomly take different subset of the particles each time we do a new "iteration" in SGD
        if (random_seed > 0) {
            mydata.randomiseParticlesOrder(random_seed+iter, do_split_random_halves,  subset_size < mydata.numberOfParticles() );
        } else if (verb > 0) {
            std::cerr << " WARNING: skipping randomisation of particle order because random_seed equals zero..." << std::endl;
        }

        if (do_auto_refine) {
            // Check whether we have converged by now
            // If we have, set do_join_random_halves and do_use_all_data for the next iteration
            checkConvergence();
        }

        expectation();


        // Sjors & Shaoda Apr 2015
        // This function does enforceHermitianSymmetry, applyHelicalSymmetry and applyPointGroupSymmetry sequentially.
        // First it enforces Hermitian symmetry to the back-projected Fourier 3D matrix.
        // Then helical symmetry is applied in Fourier space. It does rise and twist for all asymmetrical units in Fourier space.
        // Finally it applies point group symmetry (such as Cn, ...0).
        // DEBUG
        if (verb > 0) {
            if (do_helical_refine && !ignore_helical_symmetry) {
                if (mymodel.helical_nr_asu > 1)
                    std::cout << " Applying helical symmetry from the last iteration for all asymmetrical units in Fourier space..." << std::endl;
                if (iter > 1 && do_helical_symmetry_local_refinement) {
                    std::cout << " Refining helical symmetry in real space..." << std::endl;
                    std::cout << " Applying refined helical symmetry in real space..." << std::endl;
                } else {
                    std::cout << " Applying helical symmetry from the last iteration in real space..." << std::endl;
                }
            }
        }
        symmetriseReconstructions();

        }

        {
        ifdefTIMING(TicToc tt (timer, TIMING_MAX);)

        if (do_skip_maximization) {
            // Only write data.star file and break from the iteration loop
            write(DONT_WRITE_SAMPLING, DO_WRITE_DATA, DONT_WRITE_OPTIMISER, DONT_WRITE_MODEL, 0);
            break;
        }

        maximization();
        }

        {
        ifdefTIMING(TicToc tt (timer, TIMING_ITER_LOCALSYM);)
        // Apply local symmetry according to a list of masks and their operators
        applyLocalSymmetryForEachRef();
        }

        {
        ifdefTIMING(TicToc tt (timer, TIMING_ITER_HELICALREFINE);)
        // Shaoda 26 Jul 2015
        // Helical symmetry refinement and imposition of real space helical symmetry.
        if (do_helical_refine && mymodel.ref_dim == 3) {
            if (!ignore_helical_symmetry) {
                makeGoodHelixForEachRef();
            }

            if (!do_skip_align && !do_skip_rotate) {
                updatePriorsForHelicalReconstruction(
                    mydata.MDimg,
                    helical_sigma_distance * (RFLOAT) mymodel.ori_size,
                    mymodel.helical_rise,
                    mymodel.helical_twist,
                    helical_nstart,
                    mymodel.data_dim == 3,
                    do_auto_refine,
                    mymodel.sigma2_rot,
                    mymodel.sigma2_tilt,
                    mymodel.sigma2_psi,
                    mymodel.sigma2_offset,
                    helical_keep_tilt_prior_fixed,
                    verb
                );
            }
        }

        // Directly use fn_out, without "_it" specifier, so unmasked refs will be overwritten at every iteration
        if (do_write_unmasked_refs)
            mymodel.write(fn_out+"_unmasked", sampling, false, true);

        }

        {
        ifdefTIMING(TicToc tt (timer, TIMING_SOLVFLAT);)
        // Apply masks to the reference images
        // At the last iteration, do not mask the map for validation purposes
        if (do_solvent && !has_converged)
            solventFlatten();
        }

        {
        ifdefTIMING(TicToc tt (timer, TIMING_UPDATERES);)
        // Re-calculate the current resolution, do this before writing to get the correct values in the output files
        updateCurrentResolution();
        }

        {
        ifdefTIMING(TicToc tt (timer, TIMING_ITER_WRITE);)
        // Write output files
        write(DO_WRITE_SAMPLING, DO_WRITE_DATA, DO_WRITE_OPTIMISER, DO_WRITE_MODEL, 0);
        }

        ifdefTIMING(if (verb > 0) timer.printTimes(false);)

    }

    // Delete threads etc
    iterateWrapUp();

}

void MlOptimiser::expectation() {

    // #define DEBUG_EXP
    #ifdef DEBUG_EXP
    std::cerr << "Entering expectation" << std::endl;
    #endif

    #ifdef MKLFFT
    // Allow parallel FFTW execution
    fftw_plan_with_nthreads(nr_threads);
    #endif

    // Initialise some stuff
    // A. Update current size (may have been changed to ori_size in autoAdjustAngularSampling) and resolution pointers
    updateImageSizeAndResolutionPointers();

    // B. Initialise Fouriertransform, set weights in wsum_model to zero, initialise AB-matrices for FFT-phase shifts, etc
    expectationSetup();

    #ifdef DEBUG_EXP
    std::cerr << "Expectation: done setup" << std::endl;
    #endif

    // C. Calculate expected minimum angular errors (only for 3D refinements)
    // And possibly update orientational sampling automatically
    if (!do_cc() && !do_skip_align && !do_only_sample_tilt && !do_sgd) {
        // Set the exp_metadata (but not the exp_imagedata which is not needed for calculateExpectedAngularErrors)
        int n_trials_acc = std::min(
            mymodel.ref_dim == 3 && mymodel.data_dim != 3 ? 100 : 10,
            (int) mydata.numberOfParticles()
        );
        getMetaAndImageDataSubset     (0, n_trials_acc - 1, false);
        calculateExpectedAngularErrors(0, n_trials_acc - 1);
    }

    // D. Update the angular sampling (all nodes except leader)
    if (
        iter > 1 && (do_auto_refine || do_sgd) ||
        mymodel.nr_classes > 1 && allow_coarser_samplings
    ) {
        updateAngularSampling();
    }

    // E. Check whether everything fits into memory
    expectationSetupCheckMemory(verb);

    // F. Precalculate AB-matrices for on-the-fly shifts
    // Use tabulated sine and cosine values instead for 2D helical segments / 3D helical sub-tomogram averaging with on-the-fly shifts
    if (do_shifts_onthefly && (!do_helical_refine || ignore_helical_symmetry) && (!do_sgd || iter <= 1))
        precalculateABMatrices();


    #ifdef DEBUG_EXP
    std::cerr << "Expectation: done setupCheckMemory" << std::endl;
    #endif

    #ifdef CUDA
    /************************************************************************/
    // GPU memory setup

    if (do_gpu) {
        for (int i = 0; i < cudaDevices.size(); i++) {
            MlDeviceBundle *b = new MlDeviceBundle(this);
            b->setDevice(cudaDevices[i]);
            b->setupFixedSizedObjects();
            accDataBundles.push_back((void*) b);
        }

        std::vector<unsigned> threadcountOnDevice(accDataBundles.size(),0);

        for (int i = 0; i < cudaOptimiserDeviceMap.size(); i++) {
            std::stringstream didSs;
            didSs << "RRt" << i;
            MlOptimiserCuda *b = new MlOptimiserCuda(this, (MlDeviceBundle*) accDataBundles[cudaOptimiserDeviceMap[i]], didSs.str().c_str());
            b->resetData();
            cudaOptimisers.push_back((void*) b);
            threadcountOnDevice[cudaOptimiserDeviceMap[i]] ++;
        }

        int devCount;
        HANDLE_ERROR(cudaGetDeviceCount(&devCount));
        HANDLE_ERROR(cudaDeviceSynchronize());
        for (int i = 0; i < accDataBundles.size(); i++) {
            if (((MlDeviceBundle*) accDataBundles[i])->device_id >= devCount || ((MlDeviceBundle*) accDataBundles[i])->device_id < 0 ) {
                // std::cerr << " using device_id=" << ((MlDeviceBundle*) accDataBundles[i])->device_id << " (device no. " << ((MlDeviceBundle*) accDataBundles[i])->device_id + 1 << ") which is not within the available device range" << devCount << std::endl;
                CRITICAL(ERR_GPUID);
            } else {
                HANDLE_ERROR(cudaSetDevice(((MlDeviceBundle*)accDataBundles[i])->device_id));
            }

            size_t free, total, allocationSize;
            HANDLE_ERROR(cudaMemGetInfo( &free, &total ));

            size_t required_free = requested_free_gpu_memory + GPU_THREAD_MEMORY_OVERHEAD_MB * 1000 * 1000 * threadcountOnDevice[i];

            if (free < required_free) {
                printf("WARNING: Ignoring required free GPU memory amount of %zu MB, due to space insufficiency.\n", required_free / 1000000);
                allocationSize = (double) free * 0.7;
            } else {
                allocationSize = free - required_free;
            }

            if (allocationSize < 200000000)
                printf("WARNING: The available space on the GPU after initialization (%zu MB) might be insufficient for the expectation step.\n", allocationSize/1000000);

            #ifdef PRINT_GPU_MEM_INFO
            printf("INFO: Projector model size %dx%dx%d\n", (int) mymodel.PPref[0].data.xdim, (int) mymodel.PPref[0].data.ydim, (int) mymodel.PPref[0].data.zdim);
            printf("INFO: Free memory for Custom Allocator of device bundle %d is %d MB\n", i, (int) ( ((float)allocationSize)/1000000.0 ) );
            #endif

            ((MlDeviceBundle*) accDataBundles[i])->setupTunableSizedObjects(allocationSize);
        }
    }
    #endif
    #ifdef ALTCPU
    if (do_cpu) {
        unsigned nr_classes = mymodel.PPref.size();
        // Allocate Array of complex arrays for this class
        if (posix_memalign((void **) &mdlClassComplex, MEM_ALIGN, nr_classes * sizeof(std::complex<XFLOAT> *)))
            CRITICAL(RAMERR);

        // Set up XFLOAT complex array shared by all threads for each class
        for (int iclass = 0; iclass < nr_classes; iclass++) {
            int mdlX = mymodel.PPref[iclass].data.xdim;
            int mdlY = mymodel.PPref[iclass].data.ydim;
            int mdlZ = mymodel.PPref[iclass].data.zdim;
            size_t mdlXYZ = mdlZ == 0 ? (size_t) mdlX * (size_t) mdlY : (size_t) mdlX * (size_t) mdlY * (size_t) mdlZ;

            try {
                mdlClassComplex[iclass] = new std::complex<XFLOAT>[mdlXYZ];
            } catch (std::bad_alloc &ba) {
                CRITICAL(RAMERR);
            }

            std::complex<XFLOAT> *pData = mdlClassComplex[iclass];

            // Copy results into complex number array
            for (size_t i = 0; i < mdlXYZ; i++) {
                std::complex<XFLOAT> arrayval(
                    (XFLOAT) mymodel.PPref[iclass].data.data[i].real,
                    (XFLOAT) mymodel.PPref[iclass].data.data[i].imag
                );
                pData[i] = arrayval;
            }
        }

        MlDataBundle *b = new MlDataBundle();
        b->setup(this);
        accDataBundles.push_back((void*) b);
    }
    #endif // ALTCPU
    /************************************************************************/

    #ifdef MKLFFT
    // Single-threaded FFTW execution for code inside parallel processing loop
    fftw_plan_with_nthreads(1);
    #endif

    // Now perform real expectation over all particles
    // Use local parameters here, as also done in the same overloaded function in MlOptimiserMpi

    long int my_nr_particles = subset_size > 0 ? subset_size : mydata.numberOfParticles();
    int barstep = std::max(1, (int) (my_nr_particles / 60));
    long int prev_barstep = 0;
    long int my_first_part_id = 0.0;
    long int my_last_part_id = my_nr_particles - 1;
    long int nr_particles_done = 0;
    if (verb > 0) {
        if (do_sgd) {
            std::cout << (do_avoid_sgd ? " Stochastic Expectation Maximisation iteration " : " Stochastic Gradient Descent iteration ") << iter << " of " << nr_iter;
        } else {
            std::cout << " Expectation iteration " << iter;
            if (!do_auto_refine)
                std::cout << " of " << nr_iter;
            if (my_nr_particles < mydata.numberOfParticles())
                std::cout << " (with " << my_nr_particles << " particles)";
        }
        std::cout << std::endl;
        init_progress_bar(my_nr_particles);
    }

    while (nr_particles_done < my_nr_particles) {

        long int my_pool_first_part_id = my_first_part_id + nr_particles_done;
        long int my_pool_last_part_id = std::min(my_last_part_id, my_pool_first_part_id + nr_pool - 1);

        {
        ifdefTIMING(TicToc tt (timer, TIMING_EXP_METADATA);)
        // Get the metadata for these particles
        getMetaAndImageDataSubset(my_pool_first_part_id, my_pool_last_part_id, !do_parallel_disc_io);
        }

        // Abort through the pipeline_control system
        if (pipeline_control_check_abort_job())
            exit(RELION_EXIT_ABORTED);

        // perform the actual expectation step on several particles
        expectationSomeParticles(my_pool_first_part_id, my_pool_last_part_id);

        {
        ifdefTIMING(TicToc tt (timer, TIMING_EXP_CHANGES);)
        // Also monitor the changes in the optimal orientations and classes
        monitorHiddenVariableChanges(my_pool_first_part_id, my_pool_last_part_id);
        }

        {
        ifdefTIMING(TicToc tt (timer, TIMING_EXP_METADATA);)
        // Set the metadata for these particles
        setMetaDataSubset(my_pool_first_part_id, my_pool_last_part_id);
        }

        nr_particles_done += my_pool_last_part_id - my_pool_first_part_id + 1;
        if (verb > 0 && nr_particles_done - prev_barstep > barstep) {
            prev_barstep = nr_particles_done;
            progress_bar(nr_particles_done);
        }

    }

    if (verb > 0)
        progress_bar(my_nr_particles);

    #ifdef CUDA
    if (do_gpu) {
        for (int i = 0; i < accDataBundles.size(); i++) {
            MlDeviceBundle *b = (MlDeviceBundle*) accDataBundles[i];
            b->syncAllBackprojects();

            for (int j = 0; j < b->backprojectors.size(); j++) {
                unsigned long s = wsum_model.BPref[j].data.size();
                XFLOAT *reals = new XFLOAT[s];
                XFLOAT *imags = new XFLOAT[s];
                XFLOAT *weights = new XFLOAT[s];

                b->backprojectors[j].getMdlData(reals, imags, weights);

                for (unsigned long n = 0; n < s; n++) {
                    wsum_model.BPref[j].data.data[n].real += (RFLOAT) reals[n];
                    wsum_model.BPref[j].data.data[n].imag += (RFLOAT) imags[n];
                    wsum_model.BPref[j].weight.data[n] += (RFLOAT) weights[n];
                }

                delete [] reals;
                delete [] imags;
                delete [] weights;

                b->projectors[j].clear();
                b->backprojectors[j].clear();
            }

            for (int j = 0; j < b->coarseProjectionPlans.size(); j++)
                b->coarseProjectionPlans[j].clear();
        }

        for (int i = 0; i < cudaOptimisers.size(); i++)
            delete (MlOptimiserCuda*) cudaOptimisers[i];

        cudaOptimisers.clear();

        for (int i = 0; i < accDataBundles.size(); i++) {

            ((MlDeviceBundle*) accDataBundles[i])->allocator->syncReadyEvents();
            ((MlDeviceBundle*) accDataBundles[i])->allocator->freeReadyAllocs();

            #ifdef DEBUG_CUDA
            if (((MlDeviceBundle*) accDataBundles[i])->allocator->getNumberOfAllocs() != 0) {
                printf("DEBUG_ERROR: Non-zero allocation count encountered in custom allocator between iterations.\n");
                ((MlDeviceBundle*) accDataBundles[i])->allocator->printState();
                fflush(stdout);
                CRITICAL(ERR_CANZ);
            }

        #endif
        }

        for (int i = 0; i < accDataBundles.size(); i++)
            delete (MlDeviceBundle*) accDataBundles[i];

        accDataBundles.clear();
    }
    #endif
    #ifdef ALTCPU
    if (do_cpu) {
        MlDataBundle* b = (MlDataBundle*) accDataBundles[0];

        for (int j = 0; j < b->backprojectors.size(); j++) {
            unsigned long s = wsum_model.BPref[j].data.size();
            XFLOAT *reals = NULL;
            XFLOAT *imags = NULL;
            XFLOAT *weights = NULL;

            b->backprojectors[j].getMdlDataPtrs(reals, imags, weights);

            for (unsigned long n = 0; n < s; n++) {
                wsum_model.BPref[j].data.data[n].real += (RFLOAT) reals[n];
                wsum_model.BPref[j].data.data[n].imag += (RFLOAT) imags[n];
                wsum_model.BPref[j].weight.data[n] += (RFLOAT) weights[n];
            }

            b->projectors[j].clear();
            b->backprojectors[j].clear();
        }

        for (int j = 0; j < b->coarseProjectionPlans.size(); j++)
            b->coarseProjectionPlans[j].clear();

        delete b;
        accDataBundles.clear();

        // Now clean up
        unsigned nr_classes = mymodel.nr_classes;
        for (int iclass = 0; iclass < nr_classes; iclass++) {
            delete [] mdlClassComplex[iclass];
        }
        free(mdlClassComplex);

        tbbCpuOptimiser.clear();
    }
    #endif
    #ifdef  MKLFFT
    // Allow parallel FFTW execution to continue now that we are outside the parallel
    // portion of expectation
    fftw_plan_with_nthreads(nr_threads);
    #endif

    // Clean up some memory
    for (int iclass = 0; iclass < mymodel.nr_classes; iclass++)
        mymodel.PPref[iclass].data.clear();

    #ifdef DEBUG_EXP
    std::cerr << "Expectation: done " << std::endl;
    #endif

}


void MlOptimiser::expectationSetup() {
    #ifdef DEBUG
    std::cerr << "Entering expectationSetup" << std::endl;
    #endif

    // Re-initialise the random seed, because with a noisy_mask, inside the previous iteration different timings of different MPI nodes may have given rise to different number of calls to ran1
    // Use the iteration number so that each iteration has a different random seed
    init_random_generator(random_seed + iter);

    // Reset the random perturbation for this sampling
    sampling.resetRandomlyPerturbedSampling();

    // Initialise Projectors and fill vector with power_spectra for all classes
    MultidimArray<RFLOAT> *my_fourier_mask = Xsize(helical_fourier_mask) > 0 ? &helical_fourier_mask : NULL;
    mymodel.setFourierTransformMaps(!fix_tau, nr_threads, strict_lowres_exp, my_fourier_mask);

    // Initialise all weighted sums to zero
    wsum_model.initZeros();

    // If we're doing SGD with gradual decrease of sigma2_fudge: calculate current fudge-factor here
    if (do_sgd && sgd_sigma2fudge_halflife > 0) {
        RFLOAT NN = (RFLOAT)(iter * subset_size);
        RFLOAT f = NN / (NN + sgd_sigma2fudge_halflife);
        // new sigma2_fudge = f * 1.0  +  (1 - f) * sgd_ini_sigma2fudge
        sigma2_fudge = f + (1.0 - f) * sgd_sigma2fudge_ini;
    }
}

void MlOptimiser::expectationSetupCheckMemory(int myverb) {

    std::vector<int> pointer_dir_nonzeroprior, pointer_psi_nonzeroprior;
    std::vector<RFLOAT> directions_prior, psi_prior;
    if (mymodel.orientational_prior_mode != NOPRIOR) {
        // First select one random direction and psi-angle for selectOrientationsWithNonZeroPriorProbability
        // This is to get an idea how many non-zero probabilities there will be
        RFLOAT ran_rot, ran_tilt, ran_psi;
        if (mymodel.nr_bodies > 1) {
            ran_rot = ran_psi = 0.0;
            ran_tilt = 90.0;
        } else {
            int randir = (int) (rnd_unif() * sampling.NrDirections());
            int ranpsi = (int) (rnd_unif() * sampling.NrPsiSamplings());
            sampling.getDirection(randir, ran_rot, ran_tilt);
            sampling.getPsiAngle(ranpsi, ran_psi);
        }
        // Calculate local searches for these angles
        // 4 Jun 2015 - Shaoda & Sjors, bimodal psi searches for helices
        if (do_helical_refine && mymodel.ref_dim == 3) {
            bool do_auto_refine_local_searches =
                do_auto_refine && sampling.healpix_order >= autosampling_hporder_local_searches;
            bool do_classification_local_searches =
                !do_auto_refine && mymodel.orientational_prior_mode == PRIOR_ROTTILT_PSI &&
                mymodel.sigma2_rot > 0.0 && mymodel.sigma2_tilt > 0.0 && mymodel.sigma2_psi > 0.0;
            bool do_local_angular_searches = do_auto_refine_local_searches || do_classification_local_searches;
            sampling.selectOrientationsWithNonZeroPriorProbabilityFor3DHelicalReconstruction(
                ran_rot, ran_tilt, ran_psi,
                sqrt(mymodel.sigma2_rot), sqrt(mymodel.sigma2_tilt), sqrt(mymodel.sigma2_psi),
                pointer_dir_nonzeroprior, directions_prior,
                pointer_psi_nonzeroprior, psi_prior,
                do_local_angular_searches
            );
        } else if (mymodel.nr_bodies > 1) {
            sampling.selectOrientationsWithNonZeroPriorProbability(
                ran_rot, ran_tilt, ran_psi,
                sqrt(mymodel.sigma2_rot), sqrt(mymodel.sigma2_tilt), sqrt(mymodel.sigma2_psi),
                pointer_dir_nonzeroprior, directions_prior,
                pointer_psi_nonzeroprior, psi_prior,
                false, 3.0,
                mymodel.sigma_tilt_bodies[0], mymodel.sigma_psi_bodies[0]
            );
        } else {
            sampling.selectOrientationsWithNonZeroPriorProbability(ran_rot, ran_tilt, ran_psi,
                                    sqrt(mymodel.sigma2_rot), sqrt(mymodel.sigma2_tilt), sqrt(mymodel.sigma2_psi),
                                    pointer_dir_nonzeroprior, directions_prior, pointer_psi_nonzeroprior, psi_prior,
                                    ((do_bimodal_psi) && (mymodel.sigma2_psi > 0.0)) );
        }
    }


    if (myverb > 0) {
        // Calculate number of sampled hidden variables:
        int nr_ang_steps = ceil(PI * particle_diameter * mymodel.current_resolution);
        RFLOAT myresol_angstep = 360.0 / nr_ang_steps;
        std::cout << " CurrentResolution= " << 1.0 / mymodel.current_resolution << " Angstroms, which requires orientationSampling of at least " << myresol_angstep
                   << " degrees for a particle of diameter " << particle_diameter << " Angstroms" << std::endl;
        for (int oversampling = 0; oversampling <= adaptive_oversampling; oversampling++) {
            std::cout << " Oversampling= " << oversampling << " NrHiddenVariableSamplingPoints= " << mymodel.nr_classes * sampling.NrSamplingPoints(oversampling, &pointer_dir_nonzeroprior, &pointer_psi_nonzeroprior) << std::endl;
            if (!sampling.fn_sym_relax.empty())
                std::cout<<"Relaxing symmetry to "<<sampling.fn_sym_relax<<std::endl;
            int nr_orient = (do_only_sample_tilt) ? sampling.NrDirections(oversampling, &pointer_dir_nonzeroprior) : sampling.NrDirections(oversampling, &pointer_dir_nonzeroprior) * sampling.NrPsiSamplings(oversampling, &pointer_psi_nonzeroprior);
            if (do_skip_rotate || do_skip_align)
                nr_orient = 1;
            std::cout << " OrientationalSampling= " << sampling.getAngularSampling(oversampling) << " NrOrientations= "<< nr_orient <<std::endl;
            if ( (do_helical_refine) && (!ignore_helical_symmetry) )
                std::cout << " TranslationalSamplingAlongHelicalAxis= " << sampling.getHelicalTranslationalSampling(oversampling) << std::flush;
            int nr_trans = (do_skip_align) ? 1 : sampling.NrTranslationalSamplings(oversampling);
            std::cout << " TranslationalSampling= " << sampling.getTranslationalSampling(oversampling)
                    << " NrTranslations= " << nr_trans << std::endl;
            std::cout << "=============================" << std::endl;
        }
    }


    if (myverb > 1) {
        // Check whether things will fit into memory
        // Each RFLOAT takes 8 bytes, and their are mymodel.nr_classes references, express in Gb
        RFLOAT Gb = sizeof(RFLOAT) / (1024.0 * 1024.0 * 1024.0);
        // A. Calculate approximate size of the reference maps
        // Forward projector has complex data, backprojector has complex data and real weight
        RFLOAT mem_references = Gb * mymodel.nr_classes * (2 * mymodel.PPref[0].data.size() + 3 * wsum_model.BPref[0].data.size());
        // B. Weight vectors
        RFLOAT mem_pool = Gb * mymodel.nr_classes * sampling.NrSamplingPoints(adaptive_oversampling,
                &pointer_dir_nonzeroprior, &pointer_psi_nonzeroprior);
        // C. The original image data
        int nr_pix = mymodel.data_dim == 2 ? mymodel.current_size * mymodel.current_size : mymodel.current_size * mymodel.current_size * mymodel.current_size;
        mem_pool += Gb * nr_pix;
        if (!do_shifts_onthefly) {
            // D. All precalculated shifted images as well (both masked and unmasked)
            mem_pool += Gb * nr_pix * 2 * sampling.NrTranslationalSamplings(adaptive_oversampling);
        }
        // Estimate the rest of the program at 0.1 Gb?
        RFLOAT mem_rest = 0.1; // This one does NOT scale with nr_pool
        // Use tabulated sine and cosine values instead for 2D helical segments / 3D helical sub-tomogram averaging with on-the-fly shifts
        if (do_shifts_onthefly && (!do_helical_refine || ignore_helical_symmetry)) {
            // E. Store all AB-matrices
            mem_rest += Gb * nr_pix * sampling.NrTranslationalSamplings(adaptive_oversampling);
        }

        RFLOAT total_mem_Gb_exp = mem_references + nr_pool * mem_pool + mem_rest;
        // Each reconstruction has to store 1 extra complex array (Fconv) and 4 extra RFLOAT arrays (Fweight, Fnewweight. vol_out and Mconv in convoluteBlobRealSpace),
        // in adddition to the RFLOAT weight-array and the complex data-array of the BPref
        // That makes a total of 2*2 + 5 = 9 * a RFLOAT array of size BPref
        RFLOAT total_mem_Gb_max = Gb * 9 * wsum_model.BPref[0].data.size();

        std::cout << " Estimated memory for expectation  step > " << total_mem_Gb_exp << " Gb."<<std::endl;
        std::cout << " Estimated memory for maximization step > " << total_mem_Gb_max << " Gb."<<std::endl;
    }

    #ifdef DEBUG
    std::cerr << "Leaving expectationSetup" << std::endl;
    #endif

}

void MlOptimiser::precalculateABMatrices() {

    global_fftshifts_ab_coarse.clear();
    global_fftshifts_ab_current.clear();
    global_fftshifts_ab2_coarse.clear();
    global_fftshifts_ab2_current.clear();
    for (int optics_group = 0; optics_group < mydata.numberOfOpticsGroups(); optics_group++) {

        std::vector<MultidimArray<Complex> > dummy;
        global_fftshifts_ab_coarse.push_back(dummy);
        global_fftshifts_ab_current.push_back(dummy);
        global_fftshifts_ab2_coarse.push_back(dummy);
        global_fftshifts_ab2_current.push_back(dummy);

        RFLOAT my_pixel_size = mydata.getOpticsPixelSize(optics_group);

        // Set the global AB-matrices for the FFT phase-shifted images
        MultidimArray<Complex> Fab_current;
        if (mymodel.data_dim == 3) {
            Fab_current.resize(
                image_current_size[optics_group],
                image_current_size[optics_group],
                image_current_size[optics_group] / 2 + 1
            );
        } else {
            Fab_current.resize(
                image_current_size[optics_group],
                image_current_size[optics_group] / 2 + 1
            );
        }
        long int exp_nr_trans = sampling.NrTranslationalSamplings();
        std::vector<RFLOAT> oversampled_translations_x, oversampled_translations_y, oversampled_translations_z;
        // Note that do_shifts_onthefly is incompatible with do_skip_align because of the loop below
        for (long int itrans = 0; itrans < exp_nr_trans; itrans++) {
            // First get the non-oversampled translations as defined by the sampling object
            // 1 Feb 2017 - Shaoda, obsolete, helical reconstuctions never call this function

            // TODO: see how this works with multiple pixel sizes......
            sampling.getTranslationsInPixel(
                itrans, 0, my_pixel_size,
                oversampled_translations_x, oversampled_translations_y, oversampled_translations_z,
                do_helical_refine && !ignore_helical_symmetry
            ); // need getTranslations to add random_perturbation

            // Precalculate AB-matrices
            RFLOAT tmp_zoff = mymodel.data_dim == 2 ? 0.0 : oversampled_translations_z[0];
            getAbMatricesForShiftImageInFourierTransform(
                Fab_current, Fab_current,
                (RFLOAT) image_full_size[optics_group],
                oversampled_translations_x[0], oversampled_translations_y[0],
                tmp_zoff
            );

            const MultidimArray<Complex> Fab_coarse = windowFourierTransform(Fab_current, image_coarse_size[optics_group]);
            global_fftshifts_ab_coarse[optics_group].push_back(Fab_coarse);
            if (adaptive_oversampling == 0) {
                global_fftshifts_ab_current[optics_group].push_back(Fab_current);
            } else {
                // Then also loop over all its oversampled relatives
                // Then loop over all its oversampled relatives
                // Feb01,2017 - Shaoda, obsolete, helical reconstuctions never call this function
                sampling.getTranslationsInPixel(
                    itrans, adaptive_oversampling, my_pixel_size,
                    oversampled_translations_x, oversampled_translations_y, oversampled_translations_z,
                    do_helical_refine && !ignore_helical_symmetry
                );
                for (long int iover_trans = 0; iover_trans < oversampled_translations_x.size(); iover_trans++) {
                    // Shift through phase-shifts in the Fourier transform
                    // Note that the shift search range is centered around (exp_old_xoff, exp_old_yoff)

                    RFLOAT tmp_zoff = mymodel.data_dim == 2 ? 0.0 : oversampled_translations_z[iover_trans];
                    getAbMatricesForShiftImageInFourierTransform(Fab_current, Fab_current, (RFLOAT)image_full_size[optics_group], oversampled_translations_x[iover_trans], oversampled_translations_y[iover_trans], tmp_zoff);

                    global_fftshifts_ab2_current[optics_group].push_back(Fab_current);
                    if (strict_highres_exp > 0.0) {
                        const MultidimArray<Complex> Fab_coarse = windowFourierTransform(Fab_current, image_coarse_size[optics_group]);
                        global_fftshifts_ab2_coarse[optics_group].push_back(Fab_coarse);
                    }
                }
            }
        }

        #ifdef DEBUG_AB
        std::cerr << " global_fftshifts_ab_coarse[optics_group].size()= " << global_fftshifts_ab_coarse[optics_group].size() << " global_fftshifts_ab_current[optics_group].size()= " << global_fftshifts_ab_current[optics_group].size() << std::endl;
        std::cerr << " global_fftshifts_ab2_coarse[optics_group].size()= " << global_fftshifts_ab2_coarse[optics_group].size() << " global_fftshifts_ab2_current[optics_group].size()= " << global_fftshifts_ab2_current[optics_group].size() << std::endl;
        #endif
    }

}

void MlOptimiser::expectationSomeParticles(
    long int my_first_part_id, long int my_last_part_id
) {

    {
    ifdefTIMING(TicToc tt (timer, TIMING_ESP);)

    // #define DEBUG_EXPSOME
    #ifdef DEBUG_EXPSOME
    std::cerr << "Entering expectationSomeParticles..." << std::endl;
    #endif

    // Use global variables for thread visibility (before there were local ones for similar call in MPI version!)
    exp_my_first_part_id = my_first_part_id;
    exp_my_last_part_id  = my_last_part_id;

    // Make sure random division is always the same with the same seed
    if (do_generate_seeds && (
        iter == 1 && !do_firstiter_cc ||
        iter == 2 &&  do_firstiter_cc
    )) {
        // calculate the random class for these SomeParticles
        exp_random_class_some_particles.clear();
        for (long int part_id_sorted = my_first_part_id; part_id_sorted <= my_last_part_id; part_id_sorted++) {
            init_random_generator(random_seed + part_id_sorted);
            exp_random_class_some_particles.push_back(rand() % mymodel.nr_classes);
        }
    }

    // Only open/close stacks once
    fImageHandler hFile;

    // Store total number of particle images in this bunch of SomeParticles, and set translations and orientations for skip_align/rotate
    long int my_metadata_offset = 0;
    exp_imgs.clear();
    int metadata_offset = 0;
    for (long int part_id_sorted = my_first_part_id; part_id_sorted <= my_last_part_id; part_id_sorted++) {

        const long int part_id = mydata.sorted_idx[part_id_sorted];

        // If skipping alignment or rotations, then store the old translation and orientation for each particle
        // If we do local angular searches, get the previously assigned angles to center the prior
        bool do_clear = part_id_sorted == my_first_part_id;
        if (do_skip_align || do_skip_rotate) {
            // Also set the rotations
            RFLOAT old_rot  = direct::elem(exp_metadata, metadata_offset, METADATA_ROT);
            RFLOAT old_tilt = direct::elem(exp_metadata, metadata_offset, METADATA_TILT);
            RFLOAT old_psi  = direct::elem(exp_metadata, metadata_offset, METADATA_PSI);
            sampling.addOneOrientation(old_rot, old_tilt, old_psi, do_clear);
        } else if (do_only_sample_tilt) {
            if (do_clear) // only clear psi_angles for the first particle, as one psi-angle is stored for each particle!
                sampling.psi_angles.clear();
            RFLOAT old_psi = direct::elem(exp_metadata, metadata_offset, METADATA_PSI);
            sampling.psi_angles.push_back(old_psi);
        }
        if (do_skip_align) {
            // Rounded translations will be applied to the image upon reading,
            // set the unique translation in the sampling object to the fractional difference
            RFLOAT my_old_offset_x, my_old_offset_y, my_old_offset_z;
            RFLOAT rounded_offset_x, rounded_offset_y, rounded_offset_z;
            RFLOAT rot_deg, tilt_deg, psi_deg;
            my_old_offset_x = direct::elem(exp_metadata, metadata_offset, METADATA_XOFF);
            my_old_offset_y = direct::elem(exp_metadata, metadata_offset, METADATA_YOFF);
            rounded_offset_x = my_old_offset_x - round(my_old_offset_x);
            rounded_offset_y = my_old_offset_y - round(my_old_offset_y);
            if (mymodel.data_dim == 3) {
                my_old_offset_z = direct::elem(exp_metadata, metadata_offset, METADATA_ZOFF);
                rounded_offset_z = my_old_offset_z - round(my_old_offset_z);
            }
            if (do_helical_refine) {
                rot_deg  = direct::elem(exp_metadata, metadata_offset, METADATA_ROT);
                tilt_deg = direct::elem(exp_metadata, metadata_offset, METADATA_TILT);
                psi_deg  = direct::elem(exp_metadata, metadata_offset, METADATA_PSI);
            }

            // TODO: this will not work if pixel size is different for different images in one particle....
            RFLOAT my_pixel_size = mydata.getImagePixelSize(part_id, 0);
            sampling.addOneTranslation(rounded_offset_x * my_pixel_size, rounded_offset_y * my_pixel_size, rounded_offset_z * my_pixel_size,
                    do_clear, (do_helical_refine) && (!ignore_helical_symmetry), rot_deg, tilt_deg, psi_deg); // clear for first particle
        }

        // Store total number of images in this bunch of SomeParticles
        metadata_offset += mydata.numberOfImagesInParticle(part_id);

        // Sjors 7 March 2016 to prevent too high disk access... Read in all pooled images simultaneously
        // Don't do this for sub-tomograms to save RAM!
        if (do_parallel_disc_io && !do_preread_images && mymodel.data_dim != 3) {
            // Read in the actual image from disc, only open/close common stacks once
            // Read in all images, only open/close common stacks once
            for (int img_id = 0; img_id < mydata.numberOfImagesInParticle(part_id); img_id++, my_metadata_offset++) {

                // Get the filename
                FileName fn_img;
                try {
                    fn_img = mydata.getImageNameOnScratch(part_id, img_id);
                } catch (const char* errmsg) {
                    std::istringstream split (exp_fn_img);
                    for (int i = 0; i <= my_metadata_offset; i++)
                        getline(split, fn_img);
                }

                // Only open again a new stackname
                FileName fn_stack;
                long int dump;
                fn_img.decompose(dump, fn_stack);
                FileName fn_open_stack = "";
                if (fn_stack != fn_open_stack) {
                    hFile.openFile(fn_stack, WRITE_READONLY);
                    fn_open_stack = fn_stack;
                }
                Image<RFLOAT> img;
                #ifdef DEBUG_BODIES
                std::cerr << " fn_img= " << fn_img << " part_id= " << part_id << std::endl;
                #endif
                img.readFromOpenFile(fn_img, hFile, -1, false);
                img().setXmippOrigin();
                exp_imgs.push_back(img());
            }
        }
    }


    #ifdef DEBUG_EXPSOME
    std::cerr << " exp_my_first_part_id= " << exp_my_first_part_id << " exp_my_last_part_id= " << exp_my_last_part_id << std::endl;
    #endif
    if (!do_cpu) {
        // GPU and traditional CPU case - use RELION's built-in task manager to
        // process multiple particles at once
        exp_ipart_ThreadTaskDistributor->resize(my_last_part_id - my_first_part_id + 1, 1);
        exp_ipart_ThreadTaskDistributor->reset();
        global_ThreadManager->run(globalThreadExpectationSomeParticles);
    }
    #ifdef ALTCPU
    else {
        // "New" CPU case - use TBB's tasking system to process multiple
        // particles in parallel.  Like the GPU implementation, the lower-
        // level parallelism is implemented by compiler vectorization
        // (roughly equivalent to GPU "threads").
        std::atomic<int> tCount(0);

        // Set the size of the TBB thread pool for these particles
        tbb::global_control gc(tbb::global_control::max_allowed_parallelism, nr_threads);
        // process all passed particles in parallel
        //for (unsigned long i = my_first_part_id; i <= my_last_part_id; i++) {
        tbb::parallel_for (my_first_part_id, my_last_part_id + 1, [&](long int i) {
            CpuOptimiserType::reference ref = tbbCpuOptimiser.local();
            MlOptimiserCpu *cpuOptimiser = (MlOptimiserCpu *) ref;
            if (!cpuOptimiser) {
                cpuOptimiser = new MlOptimiserCpu(this, (MlDataBundle*)accDataBundles[0], "cpu_optimiser");
                cpuOptimiser->resetData();
                ref = cpuOptimiser;
                cpuOptimiser->thread_id = tCount.fetch_add(1);
            }

            cpuOptimiser->expectationOneParticle(i, cpuOptimiser->thread_id);
        });
        //}
    }
    #endif

    if (threadException) throw *threadException;

    }

}


void MlOptimiser::doThreadExpectationSomeParticles(int thread_id) {

    #ifdef TIMING
    // Only time one thread
    if (thread_id == 0)
    timer.tic(TIMING_ESP_THR);
    #endif

    size_t first_ipart = 0, last_ipart = 0;
    while (exp_ipart_ThreadTaskDistributor->getTasks(first_ipart, last_ipart)) {
        // #define DEBUG_EXPSOMETHR
        #ifdef DEBUG_EXPSOMETHR
        pthread_mutex_lock(&global_mutex);
        std::cerr << " thread_id= " << thread_id << " first_ipart= " << first_ipart << " last_ipart= " << last_ipart << std::endl;
        std::cerr << " exp_my_first_part_id= " << exp_my_first_part_id << " exp_my_last_part_id= " << exp_my_last_part_id << std::endl;
        pthread_mutex_unlock(&global_mutex);
        #endif

        for (long int ipart = first_ipart; ipart <= last_ipart; ipart++) {

            #ifdef TIMING
            // Only time one thread
            if (thread_id == 0) {
                timer.tic(TIMING_ESP_ONEPART);
            } else if (thread_id == nr_threads -1) {
                timer.tic(TIMING_ESP_ONEPARTN);
            }
            #endif

            expectationOneParticle(exp_my_first_part_id + ipart, thread_id);

            #ifdef TIMING
            // Only time one thread
            if (thread_id == 0) {
                timer.toc(TIMING_ESP_ONEPART);
            } else if (thread_id == nr_threads -1) {
                timer.toc(TIMING_ESP_ONEPARTN);
            }
            #endif

        }
    }

    #ifdef TIMING
    // Only time one thread
    if (thread_id == 0)
    timer.toc(TIMING_ESP_THR);
    #endif

}


void MlOptimiser::expectationOneParticle(long int part_id_sorted, int thread_id) {
    #ifdef TIMING
    if (part_id_sorted == exp_my_first_part_id)
    timer.tic(TIMING_ESP_INI);
    #endif

    const long int part_id = mydata.sorted_idx[part_id_sorted];

    // In the first iteration, multiple seeds will be generated
    // A single random class is selected for each pool of images, and one does not marginalise over the orientations
    // The optimal orientation is based on signal-product (rather than the signal-intensity sensitive Gaussian)
    // If do_firstiter_cc, then first perform a single iteration with K=1 and cross-correlation criteria, afterwards

    // Decide which classes to integrate over (for random class assignment in 1st iteration)
    int exp_iclass_min = 0;
    int exp_iclass_max = mymodel.nr_classes - 1;
    // low-pass filter again and generate the seeds
    if (do_generate_seeds) {
        if (iter == 1 && do_firstiter_cc) {
            // In first (CC) iter, use a single reference (and CC)
            exp_iclass_min = exp_iclass_max = 0;
        } else if (
            iter == 1 ||
            iter == 2 &&  do_firstiter_cc
        ) {
            // In second CC iter, or first iter without CC: generate the seeds
            // Now select a single random class
            // exp_part_id_sorted is already in randomized order (controlled by -seed)
            // WARNING: USING SAME iclass_min AND iclass_max FOR SomeParticles!!
            long int idx = part_id_sorted - exp_my_first_part_id;
            if (idx >= exp_random_class_some_particles.size())
                REPORT_ERROR("BUG: expectationOneParticle idx>random_class_some_particles.size()");
            exp_iclass_min = exp_iclass_max = exp_random_class_some_particles[idx];
        }
    }


    // This debug is a good one to step through the separate steps of the expectation to see where trouble lies....
    // #define DEBUG_ESP_MEM
    #ifdef DEBUG_ESP_MEM

    std::cerr << "Entering MlOptimiser::expectationOneParticle" << std::endl;
    std::cerr << " part_id= " << part_id << std::endl;
    if (thread_id == 0) {
        char c;
        std::cerr << "Before getFourierTransformsAndCtfs, press any key to continue... " << std::endl;
        std::cin >> c;
    }
    global_barrier->wait();
    #endif


    // Loop over all bodies of the multi-body refinement
    // Basically, subsequently align and store weighted sums for each body
    for (int ibody = 0; ibody < mymodel.nr_bodies; ibody++) {

        // Skip this body if keep_fixed_bodies[ibody] or if it's angular accuracy is worse than 1.5x the sampling rate
        if (mymodel.nr_bodies > 1 && mymodel.keep_fixed_bodies[ibody] > 0)
            continue;

        // Here define all kind of local arrays that will be needed
        std::vector<MultidimArray<Complex > > exp_Fimg, exp_Fimg_nomask;
        std::vector<std::vector<MultidimArray<Complex > > > exp_local_Fimgs_shifted, exp_local_Fimgs_shifted_nomask;
        std::vector<MultidimArray<RFLOAT> > exp_Fctf, exp_local_Fctf, exp_local_Minvsigma2;
        std::vector<int> exp_pointer_dir_nonzeroprior, exp_pointer_psi_nonzeroprior;
        std::vector<RFLOAT> exp_directions_prior, exp_psi_prior, exp_local_sqrtXi2;
        int exp_current_image_size, exp_current_oversampling;
        std::vector<RFLOAT> exp_highres_Xi2_img, exp_min_diff2;
        MultidimArray<RFLOAT> exp_Mweight;
        MultidimArray<bool> exp_Mcoarse_significant;
        // And from storeWeightedSums
        std::vector<RFLOAT> exp_sum_weight, exp_significant_weight, exp_max_weight;
        std::vector<Matrix1D<RFLOAT> > exp_old_offset, exp_prior;
        std::vector<RFLOAT> exp_wsum_norm_correction;
        std::vector<MultidimArray<RFLOAT> > exp_power_imgs;

        int my_nr_images = mydata.numberOfImagesInParticle(part_id);
        // Global exp_metadata array has metadata of all ori_particles. Where does my_ori_particle start?
        int metadata_offset = 0;
        for (long int iori = exp_my_first_part_id; iori <= exp_my_last_part_id; iori++) {
            if (iori == part_id_sorted)
                break;
            metadata_offset += mydata.numberOfImagesInParticle(mydata.sorted_idx[iori]);
        }

        // Resize vectors for all particles
        exp_power_imgs.resize(my_nr_images);
        exp_highres_Xi2_img.resize(my_nr_images);
        exp_Fimg.resize(my_nr_images);
        exp_Fimg_nomask.resize(my_nr_images);
        exp_Fctf.resize(my_nr_images);
        exp_old_offset.resize(my_nr_images);
        exp_prior.resize(my_nr_images);

        // Then calculate all Fourier Transform of masked and unmasked image and the CTF
        #ifdef TIMING
        if (part_id_sorted == exp_my_first_part_id) {
            timer.toc(TIMING_ESP_INI);
            timer.tic(TIMING_ESP_FT);
        }
        #endif
        getFourierTransformsAndCtfs(part_id, ibody, metadata_offset, exp_Fimg, exp_Fimg_nomask, exp_Fctf,
                exp_old_offset, exp_prior, exp_power_imgs, exp_highres_Xi2_img,
                exp_pointer_dir_nonzeroprior, exp_pointer_psi_nonzeroprior, exp_directions_prior, exp_psi_prior);

        #ifdef TIMING
        if (part_id_sorted == exp_my_first_part_id)
            timer.toc(TIMING_ESP_FT);
        #endif

        #ifdef DEBUG_ESP_MEM
        if (thread_id == 0) {
            char c;
            std::cerr << " my_nr_images= " << my_nr_images << " metadata_offset= " << metadata_offset << std::endl;
            std::cerr << "After getFourierTransformsAndCtfs, press any key to continue... " << std::endl;
            std::cin >> c;
        }
        global_barrier->wait();
        #endif

        // To deal with skipped alignments/rotations
        int exp_itrans_min, exp_itrans_max, exp_idir_min, exp_idir_max, exp_ipsi_min, exp_ipsi_max;
        if (do_skip_align) {
            exp_itrans_min = exp_itrans_max = part_id_sorted - exp_my_first_part_id;
        } else {
            exp_itrans_min = 0;
            exp_itrans_max = sampling.NrTranslationalSamplings() - 1;
        }
        if (do_skip_align || do_skip_rotate) {
            exp_idir_min = exp_idir_max = exp_ipsi_min = exp_ipsi_max =
                    part_id_sorted - exp_my_first_part_id;
        } else if (do_only_sample_tilt) {
            exp_idir_min = 0;
            exp_idir_max = sampling.NrDirections(0, &exp_pointer_dir_nonzeroprior) - 1;
            exp_ipsi_min = exp_ipsi_max = part_id_sorted - exp_my_first_part_id;
        } else {
            exp_idir_min = exp_ipsi_min = 0;
            exp_idir_max = sampling.NrDirections  (0, &exp_pointer_dir_nonzeroprior) - 1;
            exp_ipsi_max = sampling.NrPsiSamplings(0, &exp_pointer_psi_nonzeroprior) - 1;
        }

        // Initialise significant weight to minus one, so that all coarse sampling points will be handled in the first pass
        exp_significant_weight.resize(my_nr_images, -1.0);

        // Only perform a second pass when using adaptive oversampling
        int nr_sampling_passes = 1 + (adaptive_oversampling > 0);

        // Pass twice through the sampling of the entire space of rot, tilt and psi
        // The first pass uses a coarser angular sampling and possibly smaller FFTs than the second pass.
        // Only those sampling points that contribute to the highest x% of the weights in the first pass are oversampled in the second pass
        // Only those sampling points will contribute to the weighted sums in the third loop below
        for (int exp_ipass = 0; exp_ipass < nr_sampling_passes; exp_ipass++) {

            // Use coarse sampling in the first pass, oversampled one the second pass
            exp_current_oversampling = exp_ipass == 0 ? 0 : adaptive_oversampling;

            #ifdef DEBUG_ESP_MEM
            if (thread_id == 0) {
                char c;
                std::cerr << " exp_current_image_size= " << exp_current_image_size << " exp_current_oversampling= " << exp_current_oversampling << " nr_sampling_passes= " << nr_sampling_passes << std::endl;
                std::cerr << "Before getAllSquaredDifferences, use top to see memory usage and then press any key to continue... " << std::endl;
                std::cin >> c;
            }
            global_barrier->wait();
            #endif

            // Calculate the squared difference terms inside the Gaussian kernel for all hidden variables
            getAllSquaredDifferences(
                part_id, ibody, exp_ipass, exp_current_oversampling,
                metadata_offset, exp_idir_min, exp_idir_max, exp_ipsi_min, exp_ipsi_max,
                exp_itrans_min, exp_itrans_max, exp_iclass_min, exp_iclass_max, exp_min_diff2, exp_highres_Xi2_img,
                exp_Fimg, exp_Fctf, exp_Mweight, exp_Mcoarse_significant,
                exp_pointer_dir_nonzeroprior, exp_pointer_psi_nonzeroprior, exp_directions_prior, exp_psi_prior,
                exp_local_Fimgs_shifted, exp_local_Minvsigma2, exp_local_Fctf, exp_local_sqrtXi2
            );

            #ifdef DEBUG_ESP_MEM
            if (thread_id == 0) {
                char c;
                std::cerr << "After getAllSquaredDifferences, use top to see memory usage and then press any key to continue... " << std::endl;
                std::cin >> c;
            }
            global_barrier->wait();
            #endif

            // Now convert the squared difference terms to weights,
            // also calculate exp_sum_weight, and in case of adaptive oversampling also exp_significant_weight
            convertAllSquaredDifferencesToWeights(
                part_id, ibody, exp_ipass, exp_current_oversampling, metadata_offset,
                exp_idir_min, exp_idir_max, exp_ipsi_min, exp_ipsi_max,
                exp_itrans_min, exp_itrans_max, exp_iclass_min, exp_iclass_max,
                exp_Mweight, exp_Mcoarse_significant, exp_significant_weight,
                exp_sum_weight, exp_old_offset, exp_prior, exp_min_diff2,
                exp_pointer_dir_nonzeroprior, exp_pointer_psi_nonzeroprior, exp_directions_prior, exp_psi_prior
            );

        #ifdef DEBUG_ESP_MEM
        if (thread_id == 0) {
            char c;
            std::cerr << "After convertAllSquaredDifferencesToWeights, press any key to continue... " << std::endl;
            std::cin >> c;
        }
        global_barrier->wait();
        #endif

        }

        #ifdef RELION_TESTING
        std::string mode = do_gpu ? "gpu" : "cpu";
        std::cerr << " " << std::endl;
        std::cerr << " finished running diffs in  " << mode << " mode." << std::endl;
        Image<RFLOAT> tt;
        tt().resize(exp_current_image_size, exp_current_image_size);
        MultidimArray<Complex> Fimg1;
        Fimg1 = exp_local_Fimgs_shifted[0][0];
        FourierTransformer transformer;
        tt() = transformer.inverseFourierTransform(Fimg1);
        CenterFFT(tt(), -1);
        std::string fnm = mode + std::string("_out_shifted_image.mrc");
        tt.write(fnm);
        tt().resize(Xsize(Mresol_coarse[optics_group]), Ysize(Mresol_coarse[optics_group]));
        for (long int n = 0; n < tt().size(); n++) {
            tt()[n] = (RFLOAT) Mresol_coarse[optics_group][n];
        }
        fnm = mode + std::string("_out_mresol_coarse.mrc");
        tt.write(fnm);
        tt().resize(Xsize(Mresol_fine), Ysize(Mresol_fine[optics_group]));
        for (long int n = 0; n < (tt()).size(); n++) {
            tt()[n] = (RFLOAT) Mresol_fine[optics_group][n];
        }
        fnm = mode + std::string("_out_mresol_fine.mrc");
        tt.write(fnm);
        tt().resize(Xsize(exp_local_Fctfs[0]), Ysize(exp_local_Fctfs[0]));
        for (long int n = 0; n < (exp_local_Fctfs[0]).size(); n++) {
            tt()[n] = (RFLOAT) exp_local_Fctfs[0][n];
        }
        fnm = mode + std::string("_out_ctf.mrc");
        tt.write(fnm);
        // exp_local_Fctfs[0]), tt());
        // CenterFFT(tt(), -1);
        // fnm = mode + std::string("_out_ctf.mrc");
        // tt.write(fnm);
        fnm = mode + std::string("_out_10k_weights.txt");
        char *text = &fnm[0];
        freopen(text, "w", stdout);
        for (int n = 0; n < 10000; n++) {
            printf("%4.8f \n", exp_Mweight[n]); // << std::endl;
        }
        fclose(stdout);
         // exit(0);
        #endif
        // For the reconstruction step use mymodel.current_size!
        // as of 3.1 no longer needed?? CHECK!! exp_current_image_size = mymodel.current_size;

        #ifdef DEBUG_ESP_MEM
        if (thread_id == 0) {
            char c;
            std::cerr << "Before storeWeightedSums, press any key to continue... " << std::endl;
            std::cin >> c;
        }
        global_barrier->wait();
        #endif

        storeWeightedSums(
            part_id, ibody, exp_current_oversampling, metadata_offset,
            exp_idir_min, exp_idir_max, exp_ipsi_min, exp_ipsi_max,
            exp_itrans_min, exp_itrans_max, exp_iclass_min, exp_iclass_max,
            exp_min_diff2, exp_highres_Xi2_img, exp_Fimg, exp_Fimg_nomask, exp_Fctf,
            exp_power_imgs, exp_old_offset, exp_prior, exp_Mweight, exp_Mcoarse_significant,
            exp_significant_weight, exp_sum_weight, exp_max_weight,
            exp_pointer_dir_nonzeroprior, exp_pointer_psi_nonzeroprior, exp_directions_prior, exp_psi_prior,
            exp_local_Fimgs_shifted, exp_local_Fimgs_shifted_nomask, exp_local_Minvsigma2, exp_local_Fctf, exp_local_sqrtXi2
        );

        #ifdef RELION_TESTING
        // std::string mode;
        mode = do_gpu ? "gpu" : "cpu";
        std::cerr << " "<< std::endl;
        std::cerr << " finished running diffs in  " << mode << " mode." << std::endl;
        fnm = mode + std::string("_out_10k_weights_afterstore.txt");
        text = &fnm[0];
        freopen(text,"w",stdout);
        // Write the first 10k diffs to be sure
        for (int n = 0; n < 10000; n++) {
            //std::cout << exp_Mweight[n] << std::endl;
            printf("%4.8f \n", exp_Mweight[n]);
        }
        // For tests we want to exit now
        // if (iter == 2)
        // 	exit(0);

        #endif

        #ifdef DEBUG_ESP_MEM
        if (thread_id==0) {
            char c;
            std::cerr << "After storeWeightedSums, press any key to continue... " << std::endl;
            std::cin >> c;
        }
        global_barrier->wait();
        #endif

    }

    #ifdef DEBUG_BODIES
    if (part_id_sorted == round(debug1))
        exit(1);
    #endif

    #ifdef DEBUG_EXPSINGLE
    std::cerr << "Leaving expectationOneParticle..." << std::endl;
    #endif

}

void MlOptimiser::symmetriseReconstructions() {
    for (int ibody = 0; ibody < mymodel.nr_bodies; ibody++) {
        if (mymodel.nr_bodies > 1 && mymodel.keep_fixed_bodies[ibody] > 0)
            continue;

        for (int iclass = 0; iclass < mymodel.nr_classes; iclass++) {
            // either ibody or iclass can be larger than 0, never 2 at the same time!
            int ith_recons = mymodel.nr_bodies > 1 ? ibody : iclass;

            if (mymodel.pdf_class[iclass] > 0.0) {
                // Immediately after expectation process. Do rise and twist for all asymmetrical units in Fourier space
                // Also convert helical rise to pixels for BPref object
                wsum_model.BPref[ith_recons].enforceHermitianSymmetry();

                // Then apply helical and point group symmetry (order irrelevant?)
                if (mymodel.nr_bodies == 1)
                    wsum_model.BPref[ith_recons].applyHelicalSymmetry(
                        mymodel.helical_nr_asu,
                        mymodel.helical_twist[ith_recons],
                        mymodel.helical_rise[ith_recons] / mymodel.pixel_size
                    );

                // Always false if (size == 0)
                if (fn_multi_sym.size() > ith_recons) {
                    //Modify symmetry settings
                    wsum_model.BPref[ith_recons].SL.read_sym_file(fn_multi_sym[ith_recons]);
                    std::cerr << " Applying point symmetry " << fn_multi_sym[ith_recons] << " to body/class " << ith_recons << std::endl;
                }

                wsum_model.BPref[ith_recons].applyPointGroupSymmetry();
            }
        }
    }
    return;
}

void MlOptimiser::applyLocalSymmetryForEachRef() {

    if (
        fn_local_symmetry_masks.size() < 1 ||
        fn_local_symmetry_operators.size() < 1
    ) return;

    if (verb > 0)
        std::cout << " Applying local symmetry in real space according to " << fn_local_symmetry_operators.size() << " operators..." << std::endl;

    for (int ibody = 0; ibody < mymodel.nr_bodies; ibody++) {
        if (mymodel.nr_bodies > 1 && mymodel.keep_fixed_bodies[ibody] > 0)
            continue;

        for (int iclass = 0; iclass < mymodel.nr_classes; iclass++) {
            // either ibody or iclass can be larger than 0, never 2 at the same time!
            int ith_recons = mymodel.nr_bodies > 1 ? ibody : iclass;
            applyLocalSymmetry(mymodel.Iref[ith_recons], fn_local_symmetry_masks, fn_local_symmetry_operators);
        }
    }
}

void MlOptimiser::makeGoodHelixForEachRef() {
    if (!do_helical_refine || ignore_helical_symmetry || mymodel.ref_dim == 2)
        return;

    for (int ibody = 0; ibody < mymodel.nr_bodies; ibody++) {
        if (mymodel.nr_bodies > 1 && mymodel.keep_fixed_bodies[ibody] > 0)
            continue;

        for (int iclass = 0; iclass < mymodel.nr_classes; iclass++) {
            // either ibody or iclass can be larger than 0, never 2 at the same time!
            int ith_recons = mymodel.nr_bodies > 1 ? ibody : iclass;
            if (mymodel.pdf_class[iclass] > 0.0) {
                if (iter > 1 && do_helical_symmetry_local_refinement) {
                    localSearchHelicalSymmetry(
                        mymodel.Iref[ith_recons],
                        mymodel.pixel_size,
                        particle_diameter / 2.0,
                        helical_tube_inner_diameter / 2.0, helical_tube_outer_diameter / 2.0,
                        helical_z_percentage,
                        mymodel.helical_rise_min, mymodel.helical_rise_max,
                        mymodel.helical_rise_inistep, mymodel.helical_rise[iclass],
                        mymodel.helical_twist_min, mymodel.helical_twist_max,
                        mymodel.helical_twist_inistep, mymodel.helical_twist[iclass]
                    );
                }
                imposeHelicalSymmetryInRealSpace(
                    mymodel.Iref[ith_recons],
                    mymodel.pixel_size,
                    particle_diameter / 2.0,
                    helical_tube_inner_diameter / 2.0,
                    helical_tube_outer_diameter / 2.0,
                    helical_z_percentage,
                    mymodel.helical_rise[iclass],
                    mymodel.helical_twist[iclass],
                    width_mask_edge
                );
            }
        }
    }

    if (verb > 0) {
        outputHelicalSymmetryStatus(
            iter,
            helical_rise_initial,
            mymodel.helical_rise_min, mymodel.helical_rise_max,
            helical_twist_initial,
            mymodel.helical_twist_min, mymodel.helical_twist_max,
            do_helical_symmetry_local_refinement,
            mymodel.helical_rise, mymodel.helical_twist,
            0.0, 0.0, 0.0, 0.0, false,
            std::cout
        );
    }
    return;
}

void MlOptimiser::maximization() {

    if (verb > 0) {
        std::cout << " Maximization ..." << std::endl;
        init_progress_bar(mymodel.nr_classes);
    }

    // First reconstruct the images for each class
    // multi-body refinement will never get here, as it is only 3D auto-refine and that requires MPI!
    for (int iclass = 0; iclass < mymodel.nr_classes * mymodel.nr_bodies; iclass++) {
        {
        ifdefTIMING(TicToc tt (timer, RCT_1);)
        if (mymodel.pdf_class[iclass] > 0.0 || mymodel.nr_bodies > 1) {

            if (wsum_model.BPref[iclass].weight.sum() > Xmipp::epsilon) {
                MultidimArray<RFLOAT> Iref_old;
                if (do_sgd) { Iref_old = mymodel.Iref[iclass]; }

                wsum_model.BPref[iclass].updateSSNRarrays(
                    mymodel.tau2_fudge_factor,
                    mymodel.tau2_class[iclass],
                    mymodel.sigma2_class[iclass],
                    mymodel.data_vs_prior_class[iclass],
                    mymodel.fourier_coverage_class[iclass],
                    mymodel.fsc_halves_class[0],
                    false,
                    false
                );

                if (do_external_reconstruct) {
                    FileName fn_ext_root = iter > -1 ?
                        FileName::compose(fn_out + "_it", iter, "", 3) :
                        fn_out;
                    fn_ext_root =
                        FileName::compose(fn_ext_root + "_class", iclass + 1, "", 3);
                    mymodel.Iref[iclass] = wsum_model.BPref[iclass].externalReconstruct(
                        fn_ext_root,
                        mymodel.fsc_halves_class[iclass],
                        mymodel.tau2_class[iclass],
                        mymodel.sigma2_class[iclass],
                        mymodel.data_vs_prior_class[iclass],
                        do_join_random_halves || do_always_join_random_halves,
                        mymodel.tau2_fudge_factor,
                        1  // verbose
                    );
                } else {
                    mymodel.Iref[iclass] = wsum_model.BPref[iclass].reconstruct(
                        gridding_nr_iter,
                        do_map,
                        mymodel.tau2_class[iclass],
                        mymodel.tau2_fudge_factor,
                        wsum_model.pdf_class[iclass],
                        minres_map,
                        iclass == 0
                    );
                }

                if (do_sgd) {
                    // Use stochastic expectation maximisation, instead of SGD.
                    if (do_avoid_sgd) {
                        if (iter < sgd_ini_iter) {
                            for (long int n = 0; n < (mymodel.Iref[iclass]).size(); n++) {
                                mymodel.Iref[iclass][n] = std::max(0.0, mymodel.Iref[iclass][n]);
                            }
                        }
                        mymodel.Iref[iclass] = mymodel.Iref[iclass] - Iref_old;
                    }

                    // Now update formula: dV_kl^(n) = (mu) * dV_kl^(n-1) + (1-mu)*step_size*G_kl^(n)
                    // where G_kl^(n) is now in mymodel.Iref[iclass]!!!
                    for (long int n = 0; n < mymodel.Igrad[iclass].size(); n++)
                        mymodel.Igrad[iclass][n] = mu * mymodel.Igrad[iclass][n] +
                        (1.0 - mu) * sgd_stepsize * mymodel.Iref[iclass][n];

                    // update formula: V_kl^(n+1) = V_kl^(n) + dV_kl^(n)
                    mymodel.Iref[iclass] = Iref_old + mymodel.Igrad[iclass];

                    // #define DEBUG_SGD
                    #ifdef DEBUG_SGD
                    FileName fn_tmp = "grad_class" + integerToString(iclass) + ".spi";
                    Image<RFLOAT> It;
                    It()=mymodel.Igrad[iclass];
                    It.write(fn_tmp);
                    fn_tmp = "ref_class" + integerToString(iclass) + ".spi";
                    It() = mymodel.Iref[iclass];
                    It.write(fn_tmp);
                    #endif
                    // Enforce positivity?
                    // Low-pass filter according to current resolution??
                    // Some sort of regularisation may be necessary....?
                }
            }
        } else {
            (do_sgd ?
                mymodel.Igrad : // When doing SGD, keep the previous reference and re-initialise the gradient to zero
                mymodel.Iref    // When not doing SGD, initialise the reference to zero.
            )[iclass].initZeros();
        }
        }
        if (verb > 0)
            progress_bar(iclass);
    }

    {
    ifdefTIMING(TicToc tt (timer, RCT_3);)
    // Then perform the update of all other model parameters
    maximizationOtherParameters();
    }

    {
    ifdefTIMING(TicToc tt (timer, RCT_4);)
    // Keep track of changes in hidden variables
    updateOverallChangesInHiddenVariables();
    }

    if (verb > 0) progress_bar(mymodel.nr_classes);

}

void MlOptimiser::maximizationOtherParameters() {
    // Note that reconstructions are done elsewhere!
    #ifdef DEBUG
    std::cerr << "Entering maximizationOtherParameters" << std::endl;
    #endif

    RFLOAT sum_weight;
    {
    ifdefTIMING(TicToc tt (timer, RCT_5);)
    // Calculate total sum of weights, and average CTF for each class (for SSNR estimation)
    sum_weight = 0.0;
    for (int iclass = 0; iclass < mymodel.nr_classes; iclass++)
        sum_weight += wsum_model.pdf_class[iclass];

    // For multi-body refinement: it is possible we haven't done any bodies anymore, so sum_weight is zero
    // in that case we need to leave all parameters as they were
    if (sum_weight < Xmipp::epsilon) return;

    // Annealing of multiple-references in SGD
    if (do_sgd && !do_sgd_skip_anneal && mymodel.nr_classes > 1 && iter < sgd_ini_iter + sgd_inbetween_iter) {
        MultidimArray<RFLOAT> Iavg = MultidimArray<RFLOAT>::zeros(mymodel.Iref[0]);
        for (int iclass = 0; iclass < mymodel.nr_classes; iclass++) {
            Iavg += mymodel.Iref[iclass];
        }
        Iavg /= (RFLOAT) mymodel.nr_classes;

        int diffiter = iter - sgd_ini_iter;
        RFLOAT frac = RFLOAT(iter - sgd_ini_iter) / RFLOAT(sgd_inbetween_iter);
        for (int iclass = 0; iclass < mymodel.nr_classes; iclass++) {
            for (long int n = 0; n < (Iavg).size(); n++) {
                mymodel.Iref[iclass][n] *= frac;
                mymodel.Iref[iclass][n] += (1.0 - frac) * Iavg[n];
            }
        }
    }

    // Update average norm_correction, don't update norm corrections anymore for multi-body refinements!
    if (do_norm_correction  && mymodel.nr_bodies == 1) {
        mymodel.avg_norm_correction *= mu;
        mymodel.avg_norm_correction += (1.0 - mu) * wsum_model.avg_norm_correction / sum_weight;
    }

    // Don't update scales in maxCC or in multi-body refinement
    if (do_scale_correction && (iter != 1 || !do_firstiter_cc) && !do_always_cc && mymodel.nr_bodies <= 1) {
        for (int igroup = 0; igroup < mymodel.nr_groups; igroup++) {
            mymodel.scale_correction[igroup] *= mu;
            mymodel.scale_correction[igroup] += wsum_model.wsum_reference_power[igroup] > 0.0 ?
                (1.0 - mu) * wsum_model.wsum_signal_product[igroup] / wsum_model.wsum_reference_power[igroup] :
                (1.0 - mu);
        }

        // TODO! Avoid extremities in scale estimates, because they lead to catastrophic events and instabilities in refinement
        // Let's exclude anything bigger than 5x the median or smaller than 1/5 of the median...
        // Use the median instead of the mean, because it is much more robust to outliers.
        std::vector<RFLOAT> sorted = mymodel.scale_correction;
        std::sort(sorted.begin(), sorted.end());
        RFLOAT median = sorted[mymodel.nr_groups / 2];

        RFLOAT avg_scale_correction = 0.0, nr_part = 0.0;
        for (int igroup = 0; igroup < mymodel.nr_groups; igroup++) {

            if (mymodel.scale_correction[igroup] > 5.0 * median) {
                mymodel.scale_correction[igroup] = 5.0 * median;
            } else if (mymodel.scale_correction[igroup] < median / 5.0) {
                mymodel.scale_correction[igroup] = median / 5.0;
            }

            avg_scale_correction += (RFLOAT) mymodel.nr_particles_per_group[igroup] * mymodel.scale_correction[igroup];
            nr_part += (RFLOAT) mymodel.nr_particles_per_group[igroup];
        }

        // Constrain average scale_correction to one.
        avg_scale_correction /= nr_part;
        for (int igroup = 0; igroup < mymodel.nr_groups; igroup++) {
            mymodel.scale_correction[igroup] /= avg_scale_correction;
            // #define DEBUG_UPDATE_SCALE
            #ifdef DEBUG_UPDATE_SCALE
            if (verb > 0) {
                std::cerr << "Group " << igroup + 1 << ": scale_correction= " << mymodel.scale_correction[igroup] << std::endl;
                for (int i = 0; i < Xsize(wsum_model.wsum_reference_power[igroup]); i++)
                    if (wsum_model.wsum_reference_power[igroup](i) > 0.0)
                        std::cerr << " i= " << i << " XA= " << wsum_model.wsum_signal_product[igroup](i)
                                            << " A2= " << wsum_model.wsum_reference_power[igroup](i)
                                            << " XA/A2= " << wsum_model.wsum_signal_product[igroup](i) / wsum_model.wsum_reference_power[igroup](i) << std::endl;

            }
            #endif
        }
    }
    }

    {
    ifdefTIMING(TicToc tt (timer, RCT_6);)
    // Update model.pdf_class vector (for each k)
    for (int iclass = 0; iclass < mymodel.nr_classes; iclass++) {

        // Update pdf_class (for SGD: update with taking mu into account! For non-SGD: mu equals zero)
        mymodel.pdf_class[iclass] *= mu;
        mymodel.pdf_class[iclass] += (1.0 - mu) * wsum_model.pdf_class[iclass] / sum_weight;

        // for 2D also update priors of translations for each class!
        if (mymodel.ref_dim == 2) {
            if (wsum_model.pdf_class[iclass] > 0.0) {
                mymodel.prior_offset_class[iclass] *= mu;
                mymodel.prior_offset_class[iclass] += (1.0 - mu) * wsum_model.prior_offset_class[iclass] / sum_weight;
            } else {
                mymodel.prior_offset_class[iclass].initZeros();
            }
        }
    }

    for (int iclass = 0; iclass < mymodel.nr_classes * mymodel.nr_bodies; iclass++) {
        // Use sampling.NrDirections() to include all directions (also those with zero prior probability for any given image)
        if (!do_skip_align && !do_skip_rotate && !do_sgd) {
            for (int idir = 0; idir < sampling.NrDirections(); idir++) {
                mymodel.pdf_direction[iclass].elem(idir) *= mu;
                mymodel.pdf_direction[iclass].elem(idir) += (1.0 - mu) * wsum_model.pdf_direction[iclass].elem(idir) / sum_weight;
            }
        }
    }

    // Update sigma2_offset
    // Factor 2 because of the 2-dimensionality of the xy-plane
    if (!fix_sigma_offset && mymodel.nr_bodies == 1) {
        mymodel.sigma2_offset *= mu;
        mymodel.sigma2_offset += (1.0 - mu) * wsum_model.sigma2_offset /
        (sum_weight * (float) ((mymodel.data_dim == 3 ? 3 : 2) - (do_helical_refine && !ignore_helical_symmetry)));
    }

    // TODO: update estimates for sigma2_rot, sigma2_tilt and sigma2_psi!
    }

    {
    ifdefTIMING(TicToc tt (timer, RCT_7);)
    // Also refrain from updating sigma_noise after the first iteration with first_iter_cc!
    if (!fix_sigma_noise && !(iter == 1 && do_firstiter_cc || do_always_cc)) {
        for (int igroup = 0; igroup < mymodel.nr_groups; igroup++) {
            RFLOAT tsum = wsum_model.sigma2_noise[igroup].sum();
            if (tsum != 0) {
                // Factor 2 because of the 2-dimensionality of the complex-plane
                for (long int n = 0; n < (mymodel.sigma2_noise[igroup]).size(); n++) {
                    mymodel.sigma2_noise[igroup][n] *= mu;
                    mymodel.sigma2_noise[igroup][n] +=
                            (1.0 - mu) * wsum_model.sigma2_noise[igroup][n] /
                                (2.0 * wsum_model.sumw_group[igroup] * Npix_per_shell[n]);
                    // Watch out for all-zero sigma2 in case of CTF-premultiplication!
                    if (mydata.hasCtfPremultiplied())
                        mymodel.sigma2_noise[igroup][n] = std::max(mymodel.sigma2_noise[igroup][n], 1e-15);

                    // With unequal box sizes and pixel sizes in optics groups, some pixels in the 1D-spectra may contain zeros:
                    // in that case, set sigma2_noise to the value in the previous pixel.
                    if (mymodel.sigma2_noise[igroup][n] < 1e-14 && n > 0) {
                        mymodel.sigma2_noise[igroup][n] = mymodel.sigma2_noise[igroup][n - 1];
                    }
                }
            }
        }
    }
    }

    {
    ifdefTIMING(TicToc tt (timer, RCT_8);)
    // After the first iteration the references are always CTF-corrected
    if (do_ctf_correction)
        refs_are_ctf_corrected = true;

    // Some statistics to output
    mymodel.LL = wsum_model.LL;
    if (iter == 1 && do_firstiter_cc || do_always_cc)
        mymodel.LL /= sum_weight; // this now stores the average ccf
    mymodel.ave_Pmax = wsum_model.ave_Pmax / sum_weight;

    // After the first, special iteration, apply low-pass filter of -ini_high again
    if (iter == 1 && do_firstiter_cc) {
        initialLowPassFilterReferences();
        if (ini_high > 0.0) {
            // Adjust the tau2_class and data_vs_prior_class, because they were calculated on the unfiltered maps
            // This is merely a matter of having correct output in the model.star file (these values are not used in the calculations)
            RFLOAT radius = mymodel.ori_size * mymodel.pixel_size / ini_high;
            radius -= WIDTH_FMASK_EDGE / 2.0;
            RFLOAT radius_p = radius + WIDTH_FMASK_EDGE;

            for (int iclass = 0; iclass < mymodel.nr_classes; iclass++) {
                for (int rr = 0; rr < Xsize(mymodel.tau2_class[iclass]); rr++) {
                    RFLOAT r = (RFLOAT) rr;
                    if (r < radius) {
                        continue;
                    } else if (r > radius_p) {
                        direct::elem(mymodel.tau2_class[iclass], rr) = 0.0;
                        direct::elem(mymodel.data_vs_prior_class[iclass], rr) = 0.0;
                    } else {
                        RFLOAT raisedcos = 0.5 - 0.5 * cos(PI * (radius_p - r) / WIDTH_FMASK_EDGE);
                        direct::elem(mymodel.tau2_class[iclass], rr) *= raisedcos * raisedcos;
                        direct::elem(mymodel.data_vs_prior_class[iclass], rr) *= raisedcos * raisedcos;
                    }
                }
            }
        }

        if (do_generate_seeds && mymodel.nr_classes > 1) {
            // In the first CC-iteration only a single reference was used
            // Now copy this one reference to all K references, for seed generation in the second iteration
            for (int iclass = 1; iclass < mymodel.nr_classes; iclass++) {
                mymodel.tau2_class         [iclass] = mymodel.tau2_class[0];
                mymodel.data_vs_prior_class[iclass] = mymodel.data_vs_prior_class[0];
                mymodel.pdf_class          [iclass] = mymodel.pdf_class[0] / mymodel.nr_classes;
                mymodel.pdf_direction      [iclass] = mymodel.pdf_direction[0];
                mymodel.Iref               [iclass] = mymodel.Iref[0];
            }
            mymodel.pdf_class[0] /= mymodel.nr_classes;
        }
    }
    }

    #ifdef DEBUG
    std::cerr << "Leaving maximizationOtherParameters" << std::endl;
    #endif
}


void MlOptimiser::solventFlatten() {
    #ifdef DEBUG
    std::cerr << "Entering MlOptimiser::solventFlatten" << std::endl;
    #endif

    // If we're doing multibody refinement: don't do solvent flattening anymore. This is already done per body
    if (mymodel.nr_bodies > 1)
        return;

    // If we're doing SGD: enforce non-negativity during the first sgd_ini_iter iterations
    if (do_sgd && iter < sgd_ini_iter) {
        for (int iclass = 0; iclass < mymodel.nr_classes; iclass++) {
            for (long int n = 0; n < (mymodel.Iref[iclass]).size(); n++) {
                mymodel.Iref[iclass][n] = std::max(0.0, mymodel.Iref[iclass][n]);
            }
        }
    }

    // First read solvent mask from disc, or pre-calculate it
    Image<RFLOAT> Isolvent, Isolvent2, Ilowpass;
    Isolvent().resize(mymodel.Iref[0]);
    Isolvent().setXmippOrigin();
    Isolvent().initZeros();
    if (fn_mask == "None") {
        // 9 Jun 2015 - Shaoda, Helical refinement
        // Solvent flatten for helices has already been done in 'makeHelicalReferenceInRealSpace()'
        if (do_helical_refine && mymodel.ref_dim == 3) {
            if (ignore_helical_symmetry) {
                createCylindricalReference(
                    Isolvent(),
                    Xsize(mymodel.Iref[0]),
                    helical_tube_inner_diameter / mymodel.pixel_size,
                    helical_tube_outer_diameter / mymodel.pixel_size,
                    width_mask_edge
                );
            } else {
                FOR_ALL_ELEMENTS_IN_ARRAY3D(Isolvent(), i, j, k) {
                    Isolvent().elem(i, j, k) = 1.0;
                }
            }
        } else {
            RFLOAT radius = particle_diameter / (2.0 * mymodel.pixel_size);
            RFLOAT radius_p = radius + width_mask_edge;
            FOR_ALL_ELEMENTS_IN_ARRAY3D(Isolvent(), i, j, k) {
                RFLOAT r = euclid(i, j, k);
                Isolvent().elem(i, j, k) =
                r < radius   ? 1.0 :
                r > radius_p ? 0.0 :
                0.5 * (1 - cos(PI * (radius_p - r) / width_mask_edge));
            }
        }
    } else {
        Isolvent.read(fn_mask);
        Isolvent().setXmippOrigin();

        const auto range = minmax(Isolvent());
        if (range.first < 0.0 || range.second > 1.0)
            REPORT_ERROR("MlOptimiser::solventFlatten: ERROR solvent mask should contain values between 0 and 1 only...");
    }

    // Also read a second solvent mask if necessary
    if (fn_mask2 != "None") {
        Isolvent2.read(fn_mask2);
        Isolvent2().setXmippOrigin();
        if (!Isolvent2().sameShape(Isolvent()))
            REPORT_ERROR("MlOptimiser::solventFlatten ERROR: second solvent mask is of incorrect size.");
    }

    // Also read a lowpass mask if necessary
    if (fn_lowpass_mask != "None") {
        Ilowpass.read(fn_lowpass_mask);
        Ilowpass().setXmippOrigin();
        if (!Ilowpass().sameShape(Isolvent()))
            REPORT_ERROR("MlOptimiser::solventFlatten ERROR: second solvent mask is of incorrect size.");
    }

    for (int iclass = 0; iclass < mymodel.nr_classes; iclass++) {
        MultidimArray<RFLOAT> tmp;
        if (fn_lowpass_mask != "None") {
            tmp = mymodel.Iref[iclass] * Ilowpass();
            lowPassFilterMap(tmp, lowpass, mymodel.pixel_size);
        }

        // Then apply the expanded solvent mask to the map
        mymodel.Iref[iclass] *= Isolvent();  // this is the tight mask

        if (fn_lowpass_mask != "None")
            mymodel.Iref[iclass] += tmp;

        // Apply a second solvent mask if necessary
        // This may for example be useful to set the interior of icosahedral viruses to a constant density value that is higher than the solvent
        // Invert the solvent mask, so that an input mask can be given where 1 is the masked area and 0 is protein....
        if (fn_mask2 != "None")
            softMaskOutsideMap(mymodel.Iref[iclass], Isolvent2(), true);
    }
    #ifdef DEBUG
    std::cerr << "Leaving MlOptimiser::solventFlatten" << std::endl;
    #endif

}

void MlOptimiser::updateCurrentResolution() {
    #ifdef DEBUG
    std::cerr << "Entering MlOptimiser::updateCurrentResolution" << std::endl;
    #endif

    if (do_sgd && !do_split_random_halves) {
        // Do initial iterations with completely identical K references, 100-particle batch size, enforce non-negativity and 35A resolution limit
        if (iter < sgd_ini_iter) {
            mymodel.current_resolution = 1.0 / sgd_ini_resol;
        } else if (iter < sgd_ini_iter + sgd_inbetween_iter) {
            int newpixres = sgd_inires_pix + round((RFLOAT(iter - sgd_ini_iter) / RFLOAT(sgd_inbetween_iter)) * (sgd_finres_pix - sgd_inires_pix));
            mymodel.current_resolution = mymodel.getResolution(newpixres);
        } else {
            mymodel.current_resolution = 1.0 / sgd_fin_resol;
        }
    } else {
        RFLOAT best_current_resolution = 0.0;
        int nr_iter_wo_resol_gain_sum_bodies = 0;
        for (int ibody = 0; ibody < mymodel.nr_bodies; ibody++) {
            int maxres = 0;
            if (do_map) {
                // Set current resolution
                if (ini_high > 0.0 && (iter == 0 || iter == 1 && do_firstiter_cc)) {
                    maxres = round(mymodel.ori_size * mymodel.pixel_size / ini_high);
                } else {
                    // Calculate at which resolution shell the data_vs_prior drops below 1
                    int ires;
                    for (int iclass = 0; iclass < mymodel.nr_classes; iclass++) {
                        int iclass_body = (mymodel.nr_bodies > 1) ? ibody: iclass;
                        for (ires = 1; ires < mymodel.ori_size/2; ires++) {
                            if (direct::elem(mymodel.data_vs_prior_class[iclass_body], ires) < 1.0)
                                break;
                        }
                        // Subtract one shell to be back on the safe side
                        ires--;

                        if (do_split_random_halves && do_auto_refine) {
                            // Let's also try and check from the high-res side. Sometimes phase-randomisation gives artefacts
                            int ires2;
                            for (ires2 = mymodel.ori_size/2-1; ires2 >= ires; ires2--) {
                                if (direct::elem(mymodel.data_vs_prior_class[iclass_body], ires2) > 1.0)
                                    break;
                            }
                            if (ires2 > ires + 3) {
                                if (verb > 0) {
                                    float higher = mymodel.getResolutionAngstrom(ires2);
                                    float lower  = mymodel.getResolutionAngstrom(ires);
                                    if (mymodel.nr_bodies > 1) {
                                        std::cerr << " WARNING: For the " << ibody+1 << "th body:" << std::endl;
                                        std::cerr << " WARNING: FSC dipped below 0.5 and rose again. Using higher resolution of "
                                                  << higher << " A, instead of " << lower << " A." << std::endl;
                                        std::cerr << "          This is not necessarily a bad thing. Often it is caused by too tight masks." << std::endl;
                                    }
                                }
                                ires = ires2;
                            }
                        }
                        if (ires > maxres) { maxres = ires; }
                    }

                    // Never allow smaller maxres than minres_map
                    maxres = std::max(maxres, minres_map);
                }
            } else {
                // If we are not doing MAP-estimation, set maxres to Nyquist
                maxres = mymodel.ori_size / 2;
            }
            RFLOAT newres = mymodel.getResolution(maxres);

            // best resolution over all bodies
            best_current_resolution = std::max(best_current_resolution, newres);

            // Check whether resolution improved, if not increase nr_iter_wo_resol_gain
            // if (newres <= best_resol_thus_far)
            if (newres <= mymodel.current_resolution + 0.0001) {
                // Add 0.0001 to avoid problems due to rounding error
                nr_iter_wo_resol_gain_sum_bodies++;
            } else {
                nr_iter_wo_resol_gain = 0;
            }

            // Store best resolution thus far (but no longer do anything with it anymore...0)
            if (newres > best_resol_thus_far)
                best_resol_thus_far = newres;

        }

        // Set the new resolution to be the highest resolution over all bodies
        mymodel.current_resolution = best_current_resolution;

        if (nr_iter_wo_resol_gain_sum_bodies == mymodel.nr_bodies)
            nr_iter_wo_resol_gain++;
    }

    #ifdef DEBUG
    std::cerr << "Leaving MlOptimiser::updateCurrentResolution" << std::endl;
    #endif

}

void MlOptimiser::updateImageSizeAndResolutionPointers() {

    // Increment the current_size
    // If we are far from convergence (in the initial stages of refinement) take steps of 25% the image size
    // Do this whenever the FSC at the current_size is larger than 0.2, but NOT when this is in combination with very low Pmax values,
    // in the latter case, over-marginalisation may lead to spuriously high FSCs (2 smoothed maps may look very similar at high-res: all zero!)
    //
    int maxres = mymodel.getPixelFromResolution(mymodel.current_resolution) + (
        mymodel.ave_Pmax > 0.1 && has_high_fsc_at_limit ?
            round(mymodel.ori_size * 0.125) :
            incr_size  // If we are near our resolution limit, use incr_size (by default 10 shells)
    );

    // Go back from resolution shells (i.e. radius) to image size, which are BTW always even...
    mymodel.current_size = maxres * 2;

    // Go all the way because resolution increase may be substantial
    if (do_use_all_data)
        mymodel.current_size = mymodel.ori_size;

    // current_size can never be larger than ori_size:
    mymodel.current_size = std::min(mymodel.current_size, mymodel.ori_size);

    // The current size is also used in wsum_model (in unpacking)
    wsum_model.current_size = mymodel.current_size;

    // Calculate number of pixels per resolution shell
    Npix_per_shell.initZeros(mymodel.ori_size / 2 + 1);
    MultidimArray<RFLOAT> aux (
        mymodel.ori_size / 2 + 1, mymodel.ori_size,
        mymodel.data_dim == 3 ? mymodel.ori_size : 1
    );

    FOR_ALL_ELEMENTS_IN_FFTW_TRANSFORM(aux) {
        int ires = round(euclid(ip, jp, kp));
        // TODO: better check for volume_refine, but the same still seems to hold... Half of the yz plane (either ip<0 or kp<0 is redundant at jp==0)
        // Exclude points beyond Xsize(Npix_per_shell), and exclude half of the x=0 column that is stored twice in FFTW
        if (ires < mymodel.ori_size / 2 + 1 && (jp != 0 || ip >= 0))
            Npix_per_shell.elem(ires) += 1;
    }

    // Also set sizes for the images in all optics groups
    const int nr_optics_groups = mydata.numberOfOpticsGroups();
    image_coarse_size .resize(nr_optics_groups);
    image_current_size.resize(nr_optics_groups);
    image_full_size   .resize(nr_optics_groups);
    Mresol_fine       .resize(nr_optics_groups);
    Mresol_coarse     .resize(nr_optics_groups);
    for (int optics_group = 0; optics_group < nr_optics_groups; optics_group++) {

        const RFLOAT my_pixel_size = mydata.getOpticsPixelSize(optics_group);
        const int    my_image_size = mydata.getOpticsImageSize(optics_group);
        const RFLOAT remap_sizes   = (my_pixel_size * my_image_size) / (mymodel.pixel_size * mymodel.ori_size);

        image_full_size[optics_group] = my_image_size;
        // Remap from model size to mysize, and keep even!
        image_current_size[optics_group] = 2 * ceil(0.5 * remap_sizes * mymodel.current_size);
        // Current size can never become bigger than original image size for this optics_group!
        image_current_size[optics_group] = std::min(my_image_size, image_current_size[optics_group]);

        int my_max_coarse_size = max_coarse_size > 0 ? remap_sizes * max_coarse_size : image_full_size[optics_group];

        // Update coarse_size
        if (strict_highres_exp > 0.0) {
            // Strictly limit the coarse size to the one corresponding to strict_highres_exp
            image_coarse_size[optics_group] = 2 * round((RFLOAT) (remap_sizes * mymodel.ori_size * mymodel.pixel_size / strict_highres_exp));
        } else if (adaptive_oversampling > 0.0) {
            // Dependency of coarse_size on the angular sampling used in the first pass
            const RFLOAT rotated_distance  = sampling.getAngularSampling() * PI * particle_diameter / 360.0;
            const RFLOAT keepsafe_factor   = mymodel.ref_dim == 3 ? 1.2 : 1.5;
            const RFLOAT coarse_resolution = rotated_distance / keepsafe_factor;
            // Note coarse_size should be even-valued!
            image_coarse_size[optics_group] = 2 * ceil(remap_sizes * mymodel.pixel_size * mymodel.ori_size / coarse_resolution);
            // Coarse size can never be larger than max_coarse_size
            image_coarse_size[optics_group] = std::min(my_max_coarse_size, image_coarse_size[optics_group]);
        } else {
            image_coarse_size[optics_group] = image_current_size[optics_group];
        }

        // Coarse_size can never become bigger than current_size
        image_coarse_size[optics_group] = std::min(image_current_size[optics_group], image_coarse_size[optics_group]);

        /// Also update the resolution pointers here
        copyXYZ_fourier(Mresol_fine[optics_group],   image_current_size[optics_group], mymodel.data_dim);
        copyXYZ_fourier(Mresol_coarse[optics_group], image_coarse_size[optics_group],  mymodel.data_dim);

        // #define DEBUG_MRESOL
        #ifdef DEBUG_MRESOL
        Image<RFLOAT> img;
        img().resize(Xsize(Mresol_fine[optics_group]), Ysize(Mresol_fine[optics_group]));
        for (long int n = 0; n < img().size(); n++) {
            img()[n] = (RFLOAT) Mresol_fine[optics_group][n];
        }
        img.write("Mresol_fine.mrc");
        img().resize(Xsize(Mresol_coarse[optics_group]), Ysize(Mresol_coarse[optics_group]));
        for (long int n = 0; n < img().size(); n++) {
            img()[n] = (RFLOAT) Mresol_coarse[optics_group][n];
        }
        img.write("Mresol_coarse.mrc");
        #endif

        #ifdef DEBUG
        std::cerr << " current_size= " << mymodel.current_size << " optics_group= " << optics_group << " image_current_size= " << image_current_size[optics_group] << " image_coarse_size= " << image_coarse_size[optics_group] << " current_resolution= " << mymodel.current_resolution << std::endl;
        #endif

    }
}


void MlOptimiser::getFourierTransformsAndCtfs(
    long int part_id, int ibody, int metadata_offset,
    std::vector<MultidimArray<Complex> > &exp_Fimg,
    std::vector<MultidimArray<Complex> > &exp_Fimg_nomask,
    std::vector<MultidimArray<RFLOAT> > &exp_Fctf,
    std::vector<Matrix1D<RFLOAT> > &exp_old_offset,
    std::vector<Matrix1D<RFLOAT> > &exp_prior,
    std::vector<MultidimArray<RFLOAT> > &exp_power_img,
    std::vector<RFLOAT> &exp_highres_Xi2_img,
    std::vector<int> &exp_pointer_dir_nonzeroprior,
    std::vector<int> &exp_pointer_psi_nonzeroprior,
    std::vector<RFLOAT> &exp_directions_prior,
    std::vector<RFLOAT> &exp_psi_prior
) {

    FourierTransformer transformer;
    for (int img_id = 0; img_id < mydata.numberOfImagesInParticle(part_id); img_id++) {
        Image<RFLOAT> img, rec_img;
        MultidimArray<RFLOAT> Fctf;
        Matrix2D<RFLOAT> Aori;
        Matrix1D<RFLOAT> my_projected_com(mymodel.data_dim), my_refined_ibody_offset(mymodel.data_dim);

        // To which group do I belong?
        int group_id = mydata.getGroupId(part_id, img_id);
        // What is my optics group?
        int optics_group = mydata.getOpticsGroup(part_id, img_id);
        RFLOAT my_pixel_size = mydata.getOpticsPixelSize(optics_group);
        int my_image_size = mydata.getOpticsImageSize(optics_group);

        // metadata offset for this image in the particle
        int my_metadata_offset = metadata_offset + img_id;

        // Get the norm_correction (for multi-body refinement: still use the one from the consensus refinement!)
        RFLOAT normcorr = direct::elem(exp_metadata, my_metadata_offset, METADATA_NORM);

        // Safeguard against gold-standard separation
        if (do_split_random_halves) {
            int halfset = direct::elem(exp_metadata, my_metadata_offset, METADATA_NR_SIGN);
            if (halfset != my_halfset) {
                std::cerr << "BUG!!! halfset= " << halfset << " my_halfset= " << my_halfset << " part_id= " << part_id << std::endl;
                REPORT_ERROR("BUG! Mixing gold-standard separation!!!!");
            }
        }

        // Get the old offsets and the priors on the offsets
        // Sjors 5mar18: it is very important that my_old_offset has baseMLO->mymodel.data_dim and not just (3), as transformCartesianAndHelicalCoords will give different results!!!
        Matrix1D<RFLOAT> my_old_offset(mymodel.data_dim), my_prior(mymodel.data_dim), my_old_offset_ori;

        my_old_offset[0] = direct::elem(exp_metadata, my_metadata_offset, METADATA_XOFF);
        my_prior[0]      = direct::elem(exp_metadata, my_metadata_offset, METADATA_XOFF_PRIOR);
        my_old_offset[1] = direct::elem(exp_metadata, my_metadata_offset, METADATA_YOFF);
        my_prior[1]      = direct::elem(exp_metadata, my_metadata_offset, METADATA_YOFF_PRIOR);
        if (mymodel.data_dim == 3) {
        my_old_offset[2] = direct::elem(exp_metadata, my_metadata_offset, METADATA_ZOFF);
        my_prior[2]      = direct::elem(exp_metadata, my_metadata_offset, METADATA_ZOFF_PRIOR);
        }
        if (mymodel.nr_bodies > 1) {

            // 17 May 2017: Shift image to the projected COM for this body!
            // Aori is the original transformation matrix of the consensus refinement
            Aori = Euler::angles2matrix(
                direct::elem(exp_metadata, my_metadata_offset, METADATA_ROT),
                direct::elem(exp_metadata, my_metadata_offset, METADATA_TILT),
                direct::elem(exp_metadata, my_metadata_offset, METADATA_PSI)
            );
            my_projected_com = Aori * mymodel.com_bodies[ibody];
            // This will have made my_projected_com of size 3 again! resize to mymodel.data_dim
            my_projected_com.resize(mymodel.data_dim);


            #ifdef DEBUG_BODIES
            if (part_id == round(debug1)) {
                std::cerr << "ibody: " << ibody + 1 << " projected COM: " << my_projected_com[0] << " , " << my_projected_com[1] << std::endl;
                std::cerr << "ibody: " << ibody + 1 << " consensus offset: " << my_old_offset[0] << " , " << my_old_offset[1] << std::endl;
            }
            #endif

            // Subtract the projected COM offset, to position this body in the center
            // Also keep the my_old_offset in my_old_offset_ori
            my_old_offset_ori = my_old_offset;
            my_old_offset -= my_projected_com;

            // Also get refined offset for this body
            const int N = METADATA_LINE_LENGTH_BEFORE_BODIES + ibody * METADATA_NR_BODY_PARAMS;
            for (int i = 0; i < mymodel.data_dim; ++i) {
                my_refined_ibody_offset[i] = direct::elem(exp_metadata, my_metadata_offset, N + 3 + i);
            }

            // For multi-body refinement: set the priors of the translations to zero (i.e. everything centred around consensus offset)
            my_prior.initZeros();

            #ifdef DEBUG_BODIES
            if (part_id == round(debug1)) {
                std::cerr << "ibody: " << ibody + 1 << " refined x,y= " << direct::elem(exp_metadata, my_metadata_offset, icol_xoff)
                        << "  , " << direct::elem(exp_metadata, my_metadata_offset, icol_yoff) << std::endl;
                std::cerr << "FINAL translation ibody: " << ibody + 1 << " : " << my_old_offset[0] << " , " << my_old_offset[1] << std::endl;
            }
            #endif

        }

        // Uninitialised priors were set to 999.
        for (int i = 0; i < mymodel.data_dim; i++) {
            if (my_prior[i] > 998.99 && my_prior[i] < 999.01) {
                my_prior[i] = 0.0;
            }
        }

        // Orientational priors
        if (mymodel.nr_bodies > 1) {

            // Centre local searches around the orientation from the previous iteration, this one goes with overall sigma2_ang
            // On top of that, apply prior on the deviation from (0,0,0) with mymodel.sigma_tilt_bodies[ibody] and mymodel.sigma_psi_bodies[ibody]
            const int N = METADATA_LINE_LENGTH_BEFORE_BODIES + ibody * METADATA_NR_BODY_PARAMS;
            const RFLOAT prior_rot  = direct::elem(exp_metadata, my_metadata_offset, N + 0);
            const RFLOAT prior_tilt = direct::elem(exp_metadata, my_metadata_offset, N + 1);
            const RFLOAT prior_psi  = direct::elem(exp_metadata, my_metadata_offset, N + 2);
            sampling.selectOrientationsWithNonZeroPriorProbability(
                prior_rot, prior_tilt, prior_psi,
                sqrt(mymodel.sigma2_rot), sqrt(mymodel.sigma2_tilt), sqrt(mymodel.sigma2_psi),
                exp_pointer_dir_nonzeroprior, exp_directions_prior,
                exp_pointer_psi_nonzeroprior, exp_psi_prior,
                false, 3.0,
                mymodel.sigma_tilt_bodies[ibody],
                mymodel.sigma_psi_bodies[ibody]
            );

        } else if (mymodel.orientational_prior_mode != NOPRIOR && !do_skip_align && !do_skip_rotate && !do_only_sample_tilt) {
            // First try if there are some fixed prior angles
            // For multi-body refinements, ignore the original priors and get the refined residual angles from the previous iteration
            RFLOAT prior_rot            = direct::elem(exp_metadata, my_metadata_offset, METADATA_ROT_PRIOR);
            RFLOAT prior_tilt           = direct::elem(exp_metadata, my_metadata_offset, METADATA_TILT_PRIOR);
            RFLOAT prior_psi            = direct::elem(exp_metadata, my_metadata_offset, METADATA_PSI_PRIOR);
            RFLOAT prior_psi_flip_ratio = direct::elem(exp_metadata, my_metadata_offset, METADATA_PSI_PRIOR_FLIP_RATIO);
            RFLOAT prior_rot_flip_ratio = direct::elem(exp_metadata, my_metadata_offset, METADATA_ROT_PRIOR_FLIP_RATIO);  // Kthurber

            bool do_auto_refine_local_searches = do_auto_refine && sampling.healpix_order >= autosampling_hporder_local_searches;
            bool do_classification_local_searches = !do_auto_refine && mymodel.orientational_prior_mode == PRIOR_ROTTILT_PSI &&
                mymodel.sigma2_rot > 0.0 && mymodel.sigma2_tilt > 0.0 && mymodel.sigma2_psi > 0.0;
            bool do_local_angular_searches = do_auto_refine_local_searches || do_classification_local_searches;

            // If there were no defined priors (i.e. their values were 999.0), then use the "normal" angles
            if (abs(prior_rot - 999.0) > 0.01)
                prior_rot = direct::elem(exp_metadata, my_metadata_offset, METADATA_ROT);
            if (abs(prior_tilt - 999.0) > 0.01)
                prior_tilt = direct::elem(exp_metadata, my_metadata_offset, METADATA_TILT);
            if (abs(prior_psi - 999.0) > 0.01)
                prior_psi = direct::elem(exp_metadata, my_metadata_offset, METADATA_PSI);
            if (abs(prior_psi_flip_ratio - 999.0) > 0.01)
                prior_psi_flip_ratio = 0.5;
            // Kthurber:
            if (abs(prior_rot_flip_ratio - 999.0) > 0.01)
                prior_rot_flip_ratio = 0.5;

            // Select only those orientations that have non-zero prior probability
            // Jun04,2015 - Shaoda & Sjors, bimodal psi searches for helices
            if (do_helical_refine && mymodel.ref_dim == 3) {
                sampling.selectOrientationsWithNonZeroPriorProbabilityFor3DHelicalReconstruction(
                    prior_rot, prior_tilt, prior_psi,
                    sqrt(mymodel.sigma2_rot), sqrt(mymodel.sigma2_tilt), sqrt(mymodel.sigma2_psi),
                    exp_pointer_dir_nonzeroprior, exp_directions_prior, exp_pointer_psi_nonzeroprior, exp_psi_prior,
                    do_local_angular_searches, prior_psi_flip_ratio, prior_rot_flip_ratio
                );
            } else {
                sampling.selectOrientationsWithNonZeroPriorProbability(
                    prior_rot, prior_tilt, prior_psi,
                    sqrt(mymodel.sigma2_rot), sqrt(mymodel.sigma2_tilt), sqrt(mymodel.sigma2_psi),
                    exp_pointer_dir_nonzeroprior, exp_directions_prior, exp_pointer_psi_nonzeroprior, exp_psi_prior,
                    do_bimodal_psi && mymodel.sigma2_psi > 0.0
                );
            }

            long int nr_orients = sampling.NrDirections(0, &exp_pointer_dir_nonzeroprior) * sampling.NrPsiSamplings(0, &exp_pointer_psi_nonzeroprior);
            if (nr_orients == 0) {
                std::cerr << " sampling.NrDirections()= " << sampling.NrDirections(0, &exp_pointer_dir_nonzeroprior)
                        << " sampling.NrPsiSamplings()= " << sampling.NrPsiSamplings(0, &exp_pointer_psi_nonzeroprior) << std::endl;
                REPORT_ERROR("Zero orientations fall within the local angular search. Increase the sigma-value(s) on the orientations!");
            }
        }

        // Get the image and recimg data
        if (do_parallel_disc_io) {
            // If all followers had preread images into RAM: get those now
            if (do_preread_images) {

                img() = mydata.particles[part_id].images[img_id].img;

            } else {

                // #define DEBUG_SIMULTANEOUS_READ
                #ifdef DEBUG_SIMULTANEOUS_READ
                // Read from disc
                FileName fn_img;
                std::istringstream split(exp_fn_img);
                for (int i = 0; i <= my_metadata_offset; i++)
                    getline(split, fn_img);

                img.read(fn_img);
                img().setXmippOrigin();

                // Check that this is the same as the image in exp_imgs vector
                Image<RFLOAT> diff;
                if (my_metadata_offset >= exp_imgs.size()) {
                    std::cerr << " my_metadata_offset= " << my_metadata_offset << " exp_imgs.size()= " << exp_imgs.size() << std::endl;
                    REPORT_ERROR("BUG: my_metadata_offset runs out of bounds!");
                }
                diff() = img() - exp_imgs[my_metadata_offset];
                if (diff().max() > 1e-6) {
                    std::cerr << "metadata_offset= " << metadata_offset << " fn_img=" << fn_img << " diff avg=" << diff().average() << std::endl;
                    diff.write("diff.spi");
                    diff() = exp_imgs[metadata_offset];
                    diff.write("preread.spi");
                    img.write("img.spi");
                    REPORT_ERROR("unequal pre-read images... BUG!");
                }
                #else
                if (mymodel.data_dim == 3) {

                    // Read sub-tomograms from disc in parallel (to save RAM in exp_imgs)
                    FileName fn_img;
                    try {
                        fn_img = mydata.getImageNameOnScratch(part_id, img_id);
                    } catch (const char *errmsg) {
                        std::istringstream split (exp_fn_img);
                        for (int i = 0; i <= my_metadata_offset; i++)
                            getline(split, fn_img);
                    }
                    img.read(fn_img);
                    img().setXmippOrigin();
                } else {
                    img() = exp_imgs[my_metadata_offset];
                }
                #endif
            } if (has_converged && do_use_reconstruct_images) {

                FileName fn_recimg;
                std::istringstream split2(exp_fn_recimg);
                // Get the right line in the exp_fn_img string
                for (int i = 0; i <= my_metadata_offset; i++)
                    getline(split2, fn_recimg);
                rec_img.read(fn_recimg);
                rec_img().setXmippOrigin();
            }
        } else {

            // Unpack the image from the imagedata
            const long int n = image_full_size[optics_group];
            if (mymodel.data_dim == 3) {
                // Only allow a single image per call of this function!!! nr_pool needs to be set to 1!!!!
                // This will save memory, as we'll need to store all translated images in memory....
                copyXYZ(img(), exp_imagedata, n).setXmippOrigin();

                if (has_converged && do_use_reconstruct_images) {
                    const int offset = do_ctf_correction ? 2 * n : n;
                    copyXYZ(rec_img(), exp_imagedata, n, offset).setXmippOrigin();
                }

            } else {

                copyXY(img(), exp_imagedata, n, my_metadata_offset).setXmippOrigin();

                if (has_converged && do_use_reconstruct_images) {
                    /// TODO: this will be WRONG for multi-image particles, but I guess that's not going to happen anyway...
                    const int my_nr_particles = exp_my_last_part_id + 1 - exp_my_first_part_id;
                    ////////////// TODO: think this through for no-threads here.....
                    copyXY(rec_img(), exp_imagedata, n, my_metadata_offset + my_nr_particles).setXmippOrigin();
                }

            }
        }

        // Apply the norm_correction term
        if (do_norm_correction) {
            // #define DEBUG_NORM
            #ifdef DEBUG_NORM
            if (normcorr < 0.001 || normcorr > 1000. || mymodel.avg_norm_correction < 0.001 || mymodel.avg_norm_correction > 1000.0) {
                std::cerr << " ** normcorr= " << normcorr << std::endl;
                std::cerr << " ** mymodel.avg_norm_correction= " << mymodel.avg_norm_correction << std::endl;
                std::cerr << " ** fn_img= " << fn_img << " part_id= " << part_id << " img_id= " << img_id << std::endl;
                int group_id = mydata.getGroupId(part_id);
                std::cerr << " ml_model.sigma2_noise[group_id]= " << mymodel.sigma2_noise[group_id] << " group_id= " << group_id <<std::endl;
                std::cerr << " img_id= " << img_id << std::endl;
                REPORT_ERROR("Very small or very big (avg) normcorr!");
            }
            #endif
            img() *= mymodel.avg_norm_correction / normcorr;
        }

        // Helical reconstruction: calculate old_offset in the system of coordinates of the helix, i.e. parallel & perpendicular, depending on psi-angle!
        // For helices do NOT apply old_offset along the direction of the helix!!
        Matrix1D<RFLOAT> my_old_offset_helix_coords;
        RFLOAT rot_deg  = direct::elem(exp_metadata, my_metadata_offset, METADATA_ROT);
        RFLOAT tilt_deg = direct::elem(exp_metadata, my_metadata_offset, METADATA_TILT);
        RFLOAT psi_deg  = direct::elem(exp_metadata, my_metadata_offset, METADATA_PSI);
        //std::cerr << " rot_deg= " << rot_deg << " tilt_deg= " << tilt_deg << " psi_deg= " << psi_deg << std::endl;
        if (do_helical_refine && !ignore_helical_symmetry) {
            // Calculate my_old_offset_helix_coords from my_old_offset and psi angle
            transformCartesianAndHelicalCoords(my_old_offset, my_old_offset_helix_coords, rot_deg, tilt_deg, psi_deg, CART_TO_HELICAL_COORDS);
            #ifdef DEBUG_HELICAL_ORIENTATIONAL_SEARCH
            // May 18, 2015 - Shaoda & Sjors - Helical refinement (orientational searches)
            std::cerr << "MlOptimiser::getFourierTransformsAndCtfs()" << std::endl;
            std::cerr << " Transform old Cartesian offsets to helical ones..." << std::endl;
            if (my_old_offset.size() == 2) {
                std::cerr << "  psi_deg = " << psi_deg << " degrees" << std::endl;
                std::cerr << "  old_offset(x, y) = (" << my_old_offset[0] << ", " << my_old_offset[1] << ")" << std::endl;
                std::cerr << "  old_offset_helix(r, p) = (" << my_old_offset_helix_coords[0] << ", " << my_old_offset_helix_coords[1] << ")" << std::endl;
            } else {
                std::cerr << "  psi_deg = " << psi_deg << " degrees, tilt_deg = " << tilt_deg << " degrees"<< std::endl;
                std::cerr << "  old_offset(x, y, z) = (" << my_old_offset[0] << ", " << my_old_offset[1] << ", " << my_old_offset[2] << ")" << std::endl;
                std::cerr << "  old_offset_helix(p1, p2, z) = (" << my_old_offset_helix_coords[0] << ", " << my_old_offset_helix_coords[1] << "," << my_old_offset_helix_coords[2] << ")" << std::endl;
            }
            #endif
            // We do NOT want to accumulate the offsets in the direction along the helix (which is X in the helical coordinate system!)
            // However, when doing helical local searches, we accumulate offsets
            // Do NOT accumulate offsets in 3D classification of helices
            if (!do_skip_align && !do_skip_rotate) {
                // TODO: check whether the following lines make sense
                bool do_auto_refine_local_searches = do_auto_refine && sampling.healpix_order >= autosampling_hporder_local_searches;
                bool do_classification_local_searches = !do_auto_refine && mymodel.orientational_prior_mode == PRIOR_ROTTILT_PSI &&
                    mymodel.sigma2_rot > 0.0 && mymodel.sigma2_tilt > 0.0 && mymodel.sigma2_psi > 0.0;
                bool do_local_angular_searches = (do_auto_refine_local_searches) || (do_classification_local_searches);
                if (!do_local_angular_searches) {
                           if (mymodel.data_dim == 2) {
                        my_old_offset_helix_coords[0] = 0.0;
                    } else if (mymodel.data_dim == 3) {
                        my_old_offset_helix_coords[2] = 0.0;
                    }
                }
            }
            #ifdef DEBUG_HELICAL_ORIENTATIONAL_SEARCH
            std::cerr << " Set r (translation along helical axis) to zero..." << std::endl;
            if (my_old_offset.size() == 2) {
                std::cerr << "  old_offset_helix(r, p) = (" << my_old_offset_helix_coords[0] << ", " << my_old_offset_helix_coords[1] << ")" << std::endl;
            } else {
                std::cerr << "  old_offset_helix(p1, p2, z) = (" << my_old_offset_helix_coords[0] << ", " << my_old_offset_helix_coords[1] << "," << my_old_offset_helix_coords[2] << ")" << std::endl;
            }
            #endif
            // Now re-calculate the my_old_offset in the real (or image) system of coordinate (rotate -psi angle)
            transformCartesianAndHelicalCoords(my_old_offset_helix_coords, my_old_offset, rot_deg, tilt_deg, psi_deg, HELICAL_TO_CART_COORDS);
            #ifdef DEBUG_HELICAL_ORIENTATIONAL_SEARCH
            std::cerr << " Transform helical offsets back to Cartesian ones..." << std::endl;
            if (my_old_offset.size() == 2)
                std::cerr << "  old_offset(x, y) = (" << my_old_offset[0] << ", " << my_old_offset[1] << ")" << std::endl;
            else
                std::cerr << "  old_offset(x, y, z) = (" << my_old_offset[0] << ", " << my_old_offset[1] << ", " << my_old_offset[2] << ")" << std::endl;
            #endif
        }

        for (auto &x : my_old_offset) { x = round(x); }
        img() = translate(img(), my_old_offset, DONT_WRAP);
        if (has_converged && do_use_reconstruct_images) {
            rec_img() = translate(rec_img(), my_old_offset, DONT_WRAP);
        }

        #ifdef DEBUG_HELICAL_ORIENTATIONAL_SEARCH
        if (do_helical_refine && !ignore_helical_symmetry) {
            std::cerr << " Apply (rounded) old offsets (r = 0, p) & (psi, tilt) for helices..." << std::endl;
            if (my_old_offset.size() == 2) {
                std::cerr << "  old_offset(x, y) = (" << my_old_offset[0] << ", " << my_old_offset[1] << ")" << std::endl;
            } else {
                std::cerr << "  old_offset(x, y, z) = (" << my_old_offset[0] << ", " << my_old_offset[1] << ", " << my_old_offset[2] << ")" << std::endl;
            }
            Image<RFLOAT> tt;
            tt = img;
            tt.write("selftranslated_helix.spi");
            tt.clear();
            std::cerr << " written selftranslated_helix.spi; press any key to continue..." << std::endl;
            std::string str;
            std::cin >> str;
        }
        #endif

        if (do_helical_refine && !ignore_helical_symmetry) {
            // Transform rounded Cartesian offsets to corresponding helical ones
            transformCartesianAndHelicalCoords(my_old_offset, my_old_offset_helix_coords, rot_deg, tilt_deg, psi_deg, CART_TO_HELICAL_COORDS);
            exp_old_offset[img_id] = my_old_offset_helix_coords;
        } else {
            // For multi-bodies: store only the old refined offset, not the constant consensus offset or the projected COM of this body
            exp_old_offset[img_id] = mymodel.nr_bodies > 1 ? my_refined_ibody_offset :
                                                             my_old_offset;  // Not doing helical refinement. Rounded Cartesian offsets are stored.
        }
        #ifdef DEBUG_HELICAL_ORIENTATIONAL_SEARCH
        if (do_helical_refine && !ignore_helical_symmetry) {
            if (exp_old_offset[img_id].size() == 2)
                std::cerr << "exp_old_offset = (" << exp_old_offset[img_id][0] << ", " << exp_old_offset[img_id][1] << ")" << std::endl;
            else
                std::cerr << "exp_old_offset = (" << exp_old_offset[img_id][0] << ", " << exp_old_offset[img_id][1] << ", " << exp_old_offset[img_id][2] << ")" << std::endl;
        }
        #endif
        // Also store priors on translations
        exp_prior[img_id] = my_prior;

        // #define DEBUG_SOFTMASK
        #ifdef DEBUG_SOFTMASK
        Image<RFLOAT>(img()).write("Fimg_unmasked.spi");
        std::cerr << "written Fimg_unmasked.spi; press any key to continue..." << std::endl;
        char c;
		// std::cin >> c;
        #endif

        // Always store FT of image without mask (to be used for the reconstruction)
        MultidimArray<RFLOAT> img_aux = has_converged && do_use_reconstruct_images ? rec_img() : img();
        MultidimArray<Complex> Faux = transformer.FourierTransform(img_aux);
        MultidimArray<Complex> Fimg = windowFourierTransform(Faux, image_current_size[optics_group]);
        CenterFFTbySign(Fimg);

        // Here apply the aberration corrections if necessary
        mydata.obsModel.demodulatePhase(optics_group, Fimg);
        mydata.obsModel.divideByMtf(optics_group, Fimg);
        exp_Fimg_nomask[img_id] = Fimg;

        MultidimArray<RFLOAT> Mnoise;
        bool is_helical_segment = do_helical_refine || mymodel.ref_dim == 2 && helical_tube_outer_diameter > 0.0;
        // For multibodies: have the mask radius equal to maximum radius within body mask plus the translational offset search range
        RFLOAT my_mask_radius = (mymodel.nr_bodies > 1 ?
            mymodel.max_radius_mask_bodies[ibody] + sampling.offset_range :
            particle_diameter / 2.0) / my_pixel_size;
        if (!do_zero_mask) {
            // Make a noisy background image with the same spectrum as the sigma2_noise

            // Different MPI-distributed subsets may otherwise have different instances of the random noise below,
            // because work is on an on-demand basis and therefore variable with the timing of distinct nodes...
            // Have the seed based on the part_id, so that each particle has a different instant of the noise
            init_random_generator(random_seed + part_id); // This only serves for exact reproducibility tests with 1.3-code...

            // Create noisy image for outside the mask
            Mnoise.resize(img());
            transformer.setReal(Mnoise);
            MultidimArray<Complex> &Fnoise = transformer.getFourier();

            // Remap mymodel.sigma2_noise[group_id] onto remapped_sigma2_noise for this images's size and angpix
            MultidimArray<RFLOAT> remapped_sigma2_noise = MultidimArray<RFLOAT>::zeros(Xsize(Mnoise) / 2 + 1);
            RFLOAT remap_image_sizes = (my_image_size * my_pixel_size) / (mymodel.ori_size * mymodel.pixel_size);
            for (long int i = 0; i < Xsize(mymodel.sigma2_noise[group_id]); i++) {
                int i_remap = round(remap_image_sizes * i);
                if (i_remap < Xsize(remapped_sigma2_noise))
                    direct::elem(remapped_sigma2_noise, i_remap) = direct::elem(mymodel.sigma2_noise[group_id], i);
            }

            // Fill Fnoise with random numbers, use power spectrum of the noise for its variance
            FOR_ALL_ELEMENTS_IN_FFTW_TRANSFORM(Fnoise) {
                int ires = round(sqrt((RFLOAT) (kp * kp + ip * ip + jp * jp)));
                if (ires >= 0 && ires < Xsize(remapped_sigma2_noise)) {
                    RFLOAT sigma = sqrt(sigma2_fudge * direct::elem(remapped_sigma2_noise, ires));
                    direct::elem(Fnoise, i, j, k).real = rnd_gaus(0., sigma);
                    direct::elem(Fnoise, i, j, k).imag = rnd_gaus(0., sigma);
                } else {
                    direct::elem(Fnoise, i, j, k) = 0.0;
                }
            }
            // Back to real space Mnoise
            transformer.inverseFourierTransform();
            Mnoise.setXmippOrigin();

            // 24 May 2014 - Shaoda & Sjors, Helical refinement
            if (is_helical_segment) {
                softMaskOutsideMapForHelix(
                    img(), psi_deg, tilt_deg, my_mask_radius,
                    helical_tube_outer_diameter / (2.0 * my_pixel_size), width_mask_edge, &Mnoise
                );
            } else {
                softMaskOutsideMap(img(), my_mask_radius, (RFLOAT) width_mask_edge, &Mnoise);
            }
        } else {
            // 24 May 2014 - Shaoda & Sjors, Helical refinement
            if (is_helical_segment) {
                softMaskOutsideMapForHelix(
                    img(), psi_deg, tilt_deg, my_mask_radius,
                    helical_tube_outer_diameter / (2.0 * my_pixel_size), width_mask_edge
                );
            } else {
                softMaskOutsideMap(img(), my_mask_radius, (RFLOAT)width_mask_edge);
            }
        }
        #ifdef DEBUG_SOFTMASK
        Image<RFLOAT>(img()).write("Fimg_masked.spi");
        std::cerr << "written Fimg_masked.spi; dying now..." << std::endl;
		// exit(0);
        #endif

        // Store the Fourier Transform of the image Fimg
        Faux = transformer.FourierTransform(img());

        // Store the power_class spectrum of the whole image (to fill sigma2_noise between current_size and ori_size
        if (image_current_size[optics_group] < image_full_size[optics_group]) {
            MultidimArray<RFLOAT> spectrum = MultidimArray<RFLOAT>::zeros(image_full_size[optics_group] / 2 + 1);
            RFLOAT highres_Xi2 = 0.0;
            FOR_ALL_ELEMENTS_IN_FFTW_TRANSFORM(Faux) {
                int ires = round(sqrt((RFLOAT) (kp * kp + ip * ip + jp * jp)));
                // Skip Hermitian pairs in the x==0 column

                if (ires > 0 && ires < image_full_size[optics_group] / 2 + 1 && jp != 0 || ip >= 0) {
                    RFLOAT normFaux = norm(direct::elem(Faux, i, j, k));
                    direct::elem(spectrum, ires) += normFaux;
                    // Store sumXi2 from current_size until ori_size
                    if (ires >= image_current_size[optics_group] / 2 + 1)
                        highres_Xi2 += normFaux;
                }
            }

            // Let's use .at() here instead of [] to check whether we go outside the vectors bounds
            exp_power_img[img_id] = spectrum;
            exp_highres_Xi2_img[img_id] = highres_Xi2;
        } else {
            exp_highres_Xi2_img[img_id] = 0.0;
        }

        // We never need any resolutions higher than current_size
        // So resize the Fourier transforms
        Fimg = windowFourierTransform(Faux, image_current_size[optics_group]);
        // Inside Projector and Backprojector the origin of the Fourier Transform is centered!
        CenterFFTbySign(Fimg);

        // Also perform aberration correction on the masked image (which will be used for alignment)
        mydata.obsModel.demodulatePhase(optics_group, Fimg);
        mydata.obsModel.divideByMtf(optics_group, Fimg);

        exp_Fimg[img_id] = Fimg;

        // Also store its CTF
        Fctf.resize(Fimg);

        // Now calculate the actual CTF
        if (do_ctf_correction) {
            if (mymodel.data_dim == 3) {
                Image<RFLOAT> Ictf;
                if (do_parallel_disc_io) {
                    // Read CTF-image from disc
                    FileName fn_ctf;
                    try {
                        fn_ctf = mydata.getImageNameOnScratch(part_id, img_id, true);
                    } catch (const char *errmsg) {
                        std::istringstream split (exp_fn_ctf);
                        // Get the right line in the exp_fn_img string
                        for (int i = 0; i <= my_metadata_offset; i++)
                            getline(split, fn_ctf);
                    }
                    Ictf.read(fn_ctf);
                } else {
                    // Unpack the CTF-image from the exp_imagedata array
                    Ictf().resize(image_full_size[optics_group], image_full_size[optics_group], image_full_size[optics_group]);
                    for (long int k = 0; k < Zsize(Ictf()); k++) \
                    for (long int j = 0; j < Ysize(Ictf()); j++) \
                    for (long int i = 0; i < Xsize(Ictf()); i++) {
                        direct::elem(Ictf(), i, j, k) = direct::elem(exp_imagedata, i, j, image_full_size[optics_group] + k);
                    }
                }

                // If there is a redundant half, get rid of it
                if (Xsize(Ictf()) == Ysize(Ictf())) {
                    // Set the CTF-image in Fctf
                    Ictf().setXmippOrigin();
                    FOR_ALL_ELEMENTS_IN_FFTW_TRANSFORM(Fctf) {
                        // Use negative kp,ip and jp indices, because the origin in the ctf_img lies half a pixel to the right of the actual center....
                        direct::elem(Fctf, i, j, k) = Ictf().elem(-ip, -jp, -kp);
                    }
                } else if (Xsize(Ictf()) == Ysize(Ictf()) / 2 + 1) {
                    // otherwise, just window the CTF to the current resolution
                    Fctf = windowFourierTransform(Ictf(), Ysize(Fctf));
                } else {
                    // if dimensions are neither cubical nor FFTW, stop
                    REPORT_ERROR("3D CTF volume must be either cubical or adhere to FFTW format!");
                }
            } else {

                // Get parameters that change per-particle from the exp_metadata
                CTF ctf = CtfHelper::makeCTF(
                    &mydata.obsModel, optics_group,
                    direct::elem(exp_metadata, my_metadata_offset, METADATA_CTF_DEFOCUS_U),
                    direct::elem(exp_metadata, my_metadata_offset, METADATA_CTF_DEFOCUS_V),
                    direct::elem(exp_metadata, my_metadata_offset, METADATA_CTF_DEFOCUS_ANGLE),
                    direct::elem(exp_metadata, my_metadata_offset, METADATA_CTF_BFACTOR),
                    direct::elem(exp_metadata, my_metadata_offset, METADATA_CTF_KFACTOR),
                    direct::elem(exp_metadata, my_metadata_offset, METADATA_CTF_PHASE_SHIFT)
                );

                Fctf = CtfHelper::getFftwImage(
                    ctf, Xsize(Fctf), Ysize(Fctf),
                    image_full_size[optics_group], image_full_size[optics_group], my_pixel_size,
                    &mydata.obsModel,
                    ctf_phase_flipped, only_flip_phases, intact_ctf_first_peak, true, do_ctf_padding
                );

            }

            // #define DEBUG_CTF_FFTW_IMAGE
            #ifdef DEBUG_CTF_FFTW_IMAGE
            Image<RFLOAT>(Fctf).write("relion_ctf.spi");
            std::cerr << "Written relion_ctf.spi, now exiting..." << std::endl;
            exit(1);
            #endif
            // #define DEBUG_GETCTF
            #ifdef DEBUG_GETCTF
            std::cerr << " intact_ctf_first_peak= " << intact_ctf_first_peak << std::endl;
            ctf.write(std::cerr);
            Image<RFLOAT> tmp;
            tmp() = Fctf;
            tmp.write("Fctf.spi");
            tmp() = ctf.getCenteredImage(
                mymodel.ori_size, mymodel.ori_size, mymodel.pixel_size,
                &mydata.obsModel, ctf_phase_flipped, only_flip_phases, intact_ctf_first_peak, true
            );
            tmp.write("Fctf_cen.spi");
            std::cerr << "Written Fctf.spi, Fctf_cen.spi. Press any key to continue..." << std::endl;
            char c;
            std::cin >> c;
            #endif
        } else {
            Fctf = 1.0;
        }

        // Store Fctf
        exp_Fctf[img_id] = Fctf;

        // If we're doing multibody refinement, now subtract projections of the other bodies from both the masked and the unmasked particle
        if (mymodel.nr_bodies > 1) {
            MultidimArray<Complex> Fsum_obody = MultidimArray<Complex>::zeros(Fimg);

            for (int obody = 0; obody < mymodel.nr_bodies; obody++) {
                if (obody != ibody) {
                    // Only subtract if other body is not this body....
                    // Get the right metadata
                    int ocol_rot  = 0 + METADATA_LINE_LENGTH_BEFORE_BODIES + obody * METADATA_NR_BODY_PARAMS;
                    int ocol_tilt = 1 + METADATA_LINE_LENGTH_BEFORE_BODIES + obody * METADATA_NR_BODY_PARAMS;
                    int ocol_psi  = 2 + METADATA_LINE_LENGTH_BEFORE_BODIES + obody * METADATA_NR_BODY_PARAMS;
                    int ocol_xoff = 3 + METADATA_LINE_LENGTH_BEFORE_BODIES + obody * METADATA_NR_BODY_PARAMS;
                    int ocol_yoff = 4 + METADATA_LINE_LENGTH_BEFORE_BODIES + obody * METADATA_NR_BODY_PARAMS;
                    int ocol_zoff = 5 + METADATA_LINE_LENGTH_BEFORE_BODIES + obody * METADATA_NR_BODY_PARAMS;
                    int ocol_norm = 6 + METADATA_LINE_LENGTH_BEFORE_BODIES + obody * METADATA_NR_BODY_PARAMS;

                    // Aresi is the residual orientation for this obody
                    Matrix2D<RFLOAT> Aresi = Euler::angles2matrix(
                        direct::elem(exp_metadata, my_metadata_offset, ocol_rot),
                        direct::elem(exp_metadata, my_metadata_offset, ocol_tilt),
                        direct::elem(exp_metadata, my_metadata_offset, ocol_psi)
                    );
                    // The real orientation to be applied is the obody transformation applied and the original one
                    Matrix2D<RFLOAT> Abody = Aori * (mymodel.orient_bodies[obody]).transpose() * A_rot90 * Aresi * mymodel.orient_bodies[obody];

                    // Apply anisotropic mag and scaling
                    if (mydata.obsModel.hasMagMatrices) { Abody *= mydata.obsModel.anisoMag(optics_group); }
                    Abody *= mydata.obsModel.scaleDifference(optics_group, mymodel.ori_size, mymodel.pixel_size);

                    // The following line gets the correct pointer to account for overlap in the bodies
                    int oobody = direct::elem(mymodel.pointer_body_overlap, ibody, obody);
                    // Get the FT of the projection in the right direction
                    auto FTo = mymodel.PPref[oobody].get2DFourierTransform(
                        Fimg.xdim, Fimg.ydim, Fimg.zdim, Abody);

                    #ifdef DEBUG_BODIES
                    if (part_id == round(debug1)) {
                        /*
                        for (int j = 0; j < Xsize(exp_metadata); j++)
                            std::cerr << " j= " << j << " direct::elem(exp_metadata, my_metadata_offset, j)= " << direct::elem(exp_metadata, my_metadata_offset, j) << std::endl;
                        Matrix2D<RFLOAT> B;
                        B = (mymodel.orient_bodies[obody]).transpose() * Aresi * mymodel.orient_bodies[obody];
                        std::cerr << " B= " << B << std::endl;
                        std::cerr << " Aresi= " << Aresi << std::endl;
                        std::cerr << " mymodel.orient_bodies[obody]= " << mymodel.orient_bodies[obody] << std::endl;
                        std::cerr << " Aori= " << Aori << std::endl;
                        std::cerr << " Abody= " << Abody << std::endl;
                        std::cerr << " obody= " << obody+1 << "ocol_rot= " << direct::elem(exp_metadata, my_metadata_offset, ocol_rot)
                                << " obody= " << obody+1 << "ocol_tilt= " << direct::elem(exp_metadata, my_metadata_offset, ocol_tilt)
                                << " obody= " << obody+1 << "ocol_psi= " << direct::elem(exp_metadata, my_metadata_offset, ocol_psi)
                                << " ocol_xoff= " << direct::elem(exp_metadata, my_metadata_offset, ocol_xoff)
                                << " ocol_yoff= " << direct::elem(exp_metadata, my_metadata_offset, ocol_yoff) << std::endl;
                        */
                        Faux = windowFourierTransform(FTo, mymodel.ori_size);
                        img() = transformer.inverseFourierTransform(Faux);
                        CenterFFT(img(), -1);
                        FileName fn_img = "unshifted.spi";
                        fn_img = fn_img.insertBeforeExtension("_ibody" + integerToString(ibody+1));
                        fn_img = fn_img.insertBeforeExtension("_obody" + integerToString(obody+1));
                        img.write(fn_img);
                        std::cerr << "written " << fn_img << std::endl;
                    }
                    #endif

                    // 17May2017: Body is centered at its own COM
                    // move it back to its place in the original particle image
                    Matrix1D<RFLOAT> other_projected_com(mymodel.data_dim);

                    // Projected COM for this body (using Aori, just like above for ibody and my_projected_com!!!)
                    other_projected_com = Aori * (mymodel.com_bodies[obody]);
                    // This will have made other_projected_com of size 3 again! resize to mymodel.data_dim
                    other_projected_com.resize(mymodel.data_dim);

                    // Do the exact same as was done for the ibody, but DONT round here, as later phaseShift applied to ibody below!!!
                    other_projected_com -= my_old_offset_ori;

                    #ifdef DEBUG_BODIES
                    if (part_id == round(debug1))
                        std::cerr << " obody: " << obody + 1 << " projected COM= " << other_projected_com.transpose() << std::endl;
                        std::cerr << " obody: " << obody + 1 << " refined (x,y)= " << direct::elem(exp_metadata, my_metadata_offset, ocol_xoff)
                            << "  , " << direct::elem(exp_metadata, my_metadata_offset, ocol_yoff) << std::endl;
                    #endif

                    // Subtract refined obody-displacement
                    other_projected_com[0] -= direct::elem(exp_metadata, my_metadata_offset, ocol_xoff);
                    other_projected_com[1] -= direct::elem(exp_metadata, my_metadata_offset, ocol_yoff);
                    if (mymodel.data_dim == 3)
                    other_projected_com[2] -= direct::elem(exp_metadata, my_metadata_offset, ocol_zoff);

                    // Add the my_old_offset=selfRound(my_old_offset_ori - my_projected_com) already applied to this image for ibody
                    other_projected_com += my_old_offset;

                    #ifdef DEBUG_BODIES
                    if (part_id == round(debug1)) {
                        std::cerr << " obody: " << obody + 1 << " APPLIED translation obody= " << other_projected_com.transpose() << std::endl;
                    }
                    #endif
                    shiftImageInFourierTransform(
                        FTo, (RFLOAT) mymodel.ori_size,
                        other_projected_com[0], other_projected_com[1], mymodel.data_dim == 3 ? other_projected_com[2] : 0
                    );

                    // Sum the Fourier transforms of all the obodies
                    Fsum_obody += FTo;

                }
            }

            // Now that we have all the summed projections of the obodies, apply CTF, masks etc
            // Apply the CTF to this reference projection
            if (do_ctf_correction) {

                if (mydata.obsModel.getCtfPremultiplied(optics_group)) {
                    Fsum_obody *= Fctf * Fctf;
                } else {
                    Fsum_obody *= Fctf;
                }

                // Also do phase modulation, for beam tilt correction and other asymmetric aberrations
                mydata.obsModel.demodulatePhase(optics_group, Fsum_obody, true); // true means do_modulate_instead
                mydata.obsModel.divideByMtf    (optics_group, Fsum_obody, true); // true means do_multiply_instead
            }

            // Subtract the other-body FT from the current image FT
            // First the unmasked one, which will be used for reconstruction
            // Only do this if the flag below is true. Otherwise, use the original particles for reconstruction
            if (do_reconstruct_subtracted_bodies) {
                exp_Fimg_nomask[img_id] -= Fsum_obody;
            }

            // For the masked one, have to mask outside the circular mask to prevent negative values outside the mask in the subtracted image!
            CenterFFTbySign(Fsum_obody);
            Faux = windowFourierTransform(Fsum_obody, image_full_size[optics_group]);
            img() = transformer.inverseFourierTransform(Faux);

            #ifdef DEBUG_BODIES
            if (part_id == round(debug1)) {
                fn_img = "shifted_beforemask.spi";
                fn_img = fn_img.insertBeforeExtension("_ibody" + integerToString(ibody+1));
                img.write(fn_img);
                std::cerr << "Written::: " << fn_img << std::endl;
            }
            #endif
            softMaskOutsideMap(img(), my_mask_radius, (RFLOAT) width_mask_edge);

            #ifdef DEBUG_BODIES
            if (part_id == round(debug1)) {
                fn_img = "shifted_aftermask.spi";
                fn_img = fn_img.insertBeforeExtension("_ibody" + integerToString(ibody+1));
                img.write(fn_img);
                std::cerr << "Written::: " << fn_img << std::endl;
            }
            #endif
            // And back to Fourier space now
            Faux = transformer.FourierTransform(img());
            Fsum_obody = windowFourierTransform(Faux, image_current_size[optics_group]);
            CenterFFTbySign(Fsum_obody);

            // Subtract the other-body FT from the masked exp_Fimgs
            exp_Fimg[img_id] -= Fsum_obody;

            // 23 Jul 17: NEW: as we haven't applied the (nonROUNDED!!)  my_refined_ibody_offset yet, do this now in the FourierTransform
            shiftImageInFourierTransform(
                exp_Fimg[img_id], (RFLOAT) image_full_size[optics_group],
                my_refined_ibody_offset[0], my_refined_ibody_offset[1], mymodel.data_dim == 3 ? my_refined_ibody_offset[2] : 0.0
            );
            shiftImageInFourierTransform(
                exp_Fimg_nomask[img_id], (RFLOAT)image_full_size[optics_group],
                my_refined_ibody_offset[0], my_refined_ibody_offset[1], mymodel.data_dim == 3 ? my_refined_ibody_offset[2] : 0.0
            );

            #ifdef DEBUG_BODIES
            if (part_id == round(debug1)) {
                Faux = windowFourierTransform(exp_Fimg, image_full_size[optics_group]);
                img() = transformer.inverseFourierTransform(Faux);
                CenterFFT(img(), -1);
                fn_img = "exp_Fimgs_subtracted.spi";
                fn_img = fn_img.insertBeforeExtension("_ibody" + integerToString(ibody+1));
                img.write(fn_img);
                std::cerr << "written " << fn_img << std::endl;
                Faux = windowFourierTransform(exp_Fimg_nomask[img_id], image_full_size[optics_group]);
                img() = transformer.inverseFourierTransform(Faux);
                CenterFFT(img(), -1);
                fn_img = "exp_Fimgs_nomask_subtracted.spi";
                fn_img = fn_img.insertBeforeExtension("_ibody" + integerToString(ibody+1));
                img.write(fn_img);
                std::cerr << "written " << fn_img << std::endl;
            }
            #endif
        } // end if mymodel.nr_bodies > 1

    } // end loop img_id


    transformer.clear();

#ifdef DEBUG
    std::cerr << " leaving getFourierTransformsAndCtfs..." << std::endl;
#endif

}

void MlOptimiser::precalculateShiftedImagesCtfsAndInvSigma2s(bool do_also_unmasked, bool is_for_store_wsums,
        long int part_id, int exp_current_oversampling, int metadata_offset,
        int exp_itrans_min, int exp_itrans_max,
        std::vector<MultidimArray<Complex>> &exp_Fimg,
        std::vector<MultidimArray<Complex>> &exp_Fimg_nomask,
        std::vector<MultidimArray<RFLOAT>> &exp_Fctf,
        std::vector<std::vector<MultidimArray<Complex>>> &exp_local_Fimgs_shifted,
        std::vector<std::vector<MultidimArray<Complex>>> &exp_local_Fimgs_shifted_nomask,
        std::vector<MultidimArray<RFLOAT>>&exp_local_Fctf,
        std::vector<RFLOAT> &exp_local_sqrtXi2,
        std::vector<MultidimArray<RFLOAT>>&exp_local_Minvsigma2)
{

    #ifdef TIMING
    if (part_id == mydata.sorted_idx[exp_my_first_part_id]) {
        timer.tic(
            do_also_unmasked ?              TIMING_ESP_PRECW :
            exp_current_oversampling == 0 ? TIMING_ESP_PREC1 :
                                            TIMING_ESP_PREC2
        );
    }
    #endif

    int exp_nr_images = mydata.numberOfImagesInParticle(part_id);
    int nr_shifts = do_shifts_onthefly || do_skip_align ? exp_nr_images : exp_nr_images * sampling.NrTranslationalSamplings(exp_current_oversampling);
    // Don't re-do if nothing has changed....

    // Use pre-sized vectors instead of push_backs!!
    // NO! Use .reserve()!!!  --JZ
    exp_local_Fimgs_shifted.resize(exp_nr_images);
    if (do_also_unmasked)
        exp_local_Fimgs_shifted_nomask.resize(exp_nr_images);
    exp_local_Minvsigma2.resize(exp_nr_images);
    exp_local_Fctf.resize(exp_nr_images);
    exp_local_sqrtXi2.resize(exp_nr_images);

    for (int img_id = 0, my_trans_image = 0; img_id < exp_nr_images; img_id++) {

        int group_id = mydata.getGroupId(part_id, img_id);
        int optics_group = mydata.getOpticsGroup(part_id, img_id);
        RFLOAT my_pixel_size = mydata.getOpticsPixelSize(optics_group);
        int my_image_size = mydata.getOpticsImageSize(optics_group);

        int my_metadata_offset = metadata_offset + img_id;

        int exp_current_image_size;
        if (is_for_store_wsums) {
            // Always use full size of images for weighted sums in reconstruction!
            exp_current_image_size = image_current_size[optics_group];
        } else if (strict_highres_exp > 0.0) {
            // Use smaller images in both passes and keep a maximum on coarse_size, just like in FREALIGN
            exp_current_image_size = image_coarse_size[optics_group];
        } else if (adaptive_oversampling > 0) {
            // Use smaller images in the first pass, larger ones in the second pass
            exp_current_image_size = (exp_current_oversampling == 0) ? image_coarse_size[optics_group] : image_current_size[optics_group];
        } else {
            exp_current_image_size = image_current_size[optics_group];
        }
        bool do_ctf_invsig = exp_local_Fctf.size() > 0 ? Ysize(exp_local_Fctf[0])  != exp_current_image_size : true; // size has changed
        bool do_masked_shifts = (do_ctf_invsig || nr_shifts != exp_local_Fimgs_shifted[img_id].size()); // size or nr_shifts has changed

        MultidimArray<Complex> Fimg;
        if (do_masked_shifts) {
            Fimg = windowFourierTransform(exp_Fimg[img_id], exp_current_image_size);
            exp_local_Fimgs_shifted[img_id].resize(nr_shifts);
        }
        MultidimArray<Complex> Fimg_nomask;
        if (do_also_unmasked) {
            Fimg_nomask = windowFourierTransform(exp_Fimg_nomask[img_id], exp_current_image_size);
            exp_local_Fimgs_shifted_nomask[img_id].resize(nr_shifts);
        }

        if (do_ctf_invsig) {
            // Also precalculate the sqrt of the sum of all Xi2
            // Could exp_current_image_size ever be different from mymodel.current_size?
            // Probably therefore do it here rather than in getFourierTransforms
            if (do_cc()) {
                // std::transform_reduce
                RFLOAT sumxi2 = 0.0;
                for (long int n = 0; n < Fimg.size(); n++) {
                    sumxi2 += norm(Fimg[n]);
                }
                // Normalised cross-correlation coefficient: divide by power of reference (power of image is a constant)
                exp_local_sqrtXi2[img_id] = sqrt(sumxi2);
            }

            // Also store downsized Fctfs
            // In the second pass of the adaptive approach this will have no effect,
            // since then exp_current_image_size will be the same as the size of exp_Fctfs
            exp_local_Fctf[img_id] = windowFourierTransform(exp_Fctf[img_id], exp_current_image_size);

            // Also prepare Minvsigma2
            if (mymodel.data_dim == 3) {
                exp_local_Minvsigma2[img_id].initZeros(Zsize(Fimg), Ysize(Fimg), Xsize(Fimg));
            } else {
                exp_local_Minvsigma2[img_id].initZeros(Ysize(Fimg), Xsize(Fimg));
            }

            // Map from model_size sigma2_noise array to my_image_size
            RFLOAT remap_image_sizes = (mymodel.ori_size * mymodel.pixel_size) / (my_image_size * my_pixel_size);
            int *myMresol = Ysize(Fimg) == image_coarse_size[optics_group] ? Mresol_coarse[optics_group].data : Mresol_fine[optics_group].data;
            // With group_id and relevant size of Fimg, calculate inverse of sigma^2 for relevant parts of Mresol
            for (long int n = 0; n < (exp_local_Minvsigma2[img_id]).size(); n++) {
                int ires = myMresol[n];
                int ires_remapped = round(remap_image_sizes * ires);
                // Exclude origin (ires==0) from the Probability-calculation
                // This way we are invariant to additive factors
                if (ires > 0 && ires_remapped < Xsize(mymodel.sigma2_noise[group_id]))
                    exp_local_Minvsigma2[img_id][n] = 1.0 / (sigma2_fudge * direct::elem(mymodel.sigma2_noise[group_id], ires_remapped));
            }
        }

        //Shifts are done on the fly on the gpu, if do_gpu || do_cpu, do_shifts_onthefly is always false!
        if (do_shifts_onthefly) {
            // Store a single, down-sized version of exp_Fimg[img_id] in exp_local_Fimgs_shifted[img_id]
            if (do_masked_shifts)
                exp_local_Fimgs_shifted[img_id][0] = Fimg;
            if (do_also_unmasked)
                exp_local_Fimgs_shifted_nomask[img_id][0] = Fimg_nomask;
            #ifdef DEBUG_HELICAL_ORIENTATIONAL_SEARCH
            std::cerr << " MlOptimiser::precalculateShiftedImagesCtfsAndInvSigma2s(): do_shifts_onthefly && !do_gpu" << std::endl;
            #endif
        } else if (!do_gpu && !do_cpu) {
            #ifdef DEBUG_HELICAL_ORIENTATIONAL_SEARCH
            Image<RFLOAT> img_save_ori    = Image<RFLOAT>::zeros(mymodel.data_dim == 3 ? mymodel.ori_size : 1, mymodel.ori_size, mymodel.ori_size);
            Image<RFLOAT> img_save_mask   = Image<RFLOAT>::zeros(mymodel.data_dim == 3 ? mymodel.ori_size : 1, mymodel.ori_size, mymodel.ori_size);
            Image<RFLOAT> img_save_nomask = Image<RFLOAT>::zeros(mymodel.data_dim == 3 ? mymodel.ori_size : 1, mymodel.ori_size, mymodel.ori_size);
            #endif
            // Store all translated variants of Fimg
            int my_trans_image = 0;
            for (long int itrans = exp_itrans_min; itrans <= exp_itrans_max; itrans++) {

                // First get the non-oversampled translations as defined by the sampling object
                std::vector<RFLOAT> oversampled_translations_x, oversampled_translations_y, oversampled_translations_z;
                // Jun01,2014 - Shaoda & Sjors, Helical refinement
                sampling.getTranslationsInPixel(
                    itrans, exp_current_oversampling, my_pixel_size,
                    oversampled_translations_x, oversampled_translations_y, oversampled_translations_z,
                    do_helical_refine && !ignore_helical_symmetry
                );
                #ifdef DEBUG_HELICAL_ORIENTATIONAL_SEARCH
                std::cerr << "MlOptimiser::precalculateShiftedImagesCtfsAndInvSigma2s(): Store all translated variants of Fimg" << std::endl;
                #endif
                // Then loop over all its oversampled relatives
                for (long int iover_trans = 0; iover_trans < oversampled_translations_x.size(); iover_trans++, my_trans_image++) {
                    // Helical reconstruction: rotate oversampled_translations_x[iover_trans] and oversampled_translations_y[iover_trans] according to rlnAnglePsi of this particle!
                    RFLOAT xshift = 0., yshift = 0., zshift = 0.0;

                    xshift = oversampled_translations_x[iover_trans];
                    yshift = oversampled_translations_y[iover_trans];
                    if (mymodel.data_dim == 3)
                    zshift = oversampled_translations_z[iover_trans];

                    if (do_helical_refine && !ignore_helical_symmetry) {
                        RFLOAT rot_deg  = direct::elem(exp_metadata, my_metadata_offset, METADATA_ROT);
                        RFLOAT tilt_deg = direct::elem(exp_metadata, my_metadata_offset, METADATA_TILT);
                        RFLOAT psi_deg  = direct::elem(exp_metadata, my_metadata_offset, METADATA_PSI);
                        #ifdef DEBUG_HELICAL_ORIENTATIONAL_SEARCH
                        std::cerr << "Helical xyz shift = (" << xshift << ", " << yshift << ", " << zshift << ")" << std::endl;
                        #endif
                        transformCartesianAndHelicalCoords(
                            xshift, yshift, zshift,
                            rot_deg, tilt_deg, psi_deg,
                            mymodel.data_dim, HELICAL_TO_CART_COORDS
                        );
                        #ifdef DEBUG_HELICAL_ORIENTATIONAL_SEARCH
                        std::cerr << "Cartesian xyz shift = (" << xshift << ", " << yshift << ", " << zshift << ")" << std::endl;
                        #endif
                    }

                    // Shift through phase-shifts in the Fourier transform
                    // Note that the shift search range is centered around (exp_old_xoff, exp_old_yoff)
                    if (do_masked_shifts) {
                        exp_local_Fimgs_shifted[img_id][my_trans_image].resize(Fimg);
                        shiftImageInFourierTransform(Fimg, exp_local_Fimgs_shifted[img_id][my_trans_image], (RFLOAT) mymodel.ori_size, xshift, yshift, zshift);
                        #ifdef DEBUG_HELICAL_ORIENTATIONAL_SEARCH
                        if (do_helical_refine && !ignore_helical_symmetry) {
                            // Shall we let 2D classification do this as well?
                            std::cerr << " Size of Fourier map (Z, Y, X) = "
                                    << Zsize(exp_local_Fimgs_shifted[img_id][my_trans_image]) << ", "
                                    << Ysize(exp_local_Fimgs_shifted[img_id][my_trans_image]) << ", "
                                    << Xsize(exp_local_Fimgs_shifted[img_id][my_trans_image]) << std::endl;
                            std::cerr << " mymodel.ori_size = " << mymodel.ori_size << std::endl;
                            Image<RFLOAT> tt;
                            FourierTransformer transformer;
                            tt().resize(
                                mymodel.data_dim == 3 ? mymodel.ori_size : 1, mymodel.ori_size,
                                mymodel.ori_size
                            );
                            MultidimArray<Complex> Faux = exp_local_Fimgs_shifted[img_id][my_trans_image];
                            MultidimArray<Complex> Fo = windowFourierTransform(Faux, mymodel.ori_size);
                            tt() = transformer.inverseFourierTransform(Fo);
                            CenterFFT(tt(), -1);
                            img_save_mask() += tt();
                            img_save_mask.write("translational_searches_mask_helix.spi");
                            std::cerr << " written translational_searches_mask_helix.spi; press any key to continue..." << std::endl;
                            std::string str;
                            std::cin >> str;
                        }
                        #endif
                    }
                    if (do_also_unmasked) {
                        exp_local_Fimgs_shifted_nomask[img_id][my_trans_image].resize(Fimg_nomask);
                        shiftImageInFourierTransform(Fimg_nomask, exp_local_Fimgs_shifted_nomask[img_id][my_trans_image], (RFLOAT)mymodel.ori_size, xshift, yshift, zshift);
                        #ifdef DEBUG_HELICAL_ORIENTATIONAL_SEARCH
                        if (do_helical_refine && !ignore_helical_symmetry) {
                            std::cerr << " Size of Fourier map (Z, Y, X) = "
                                    << Zsize(exp_local_Fimgs_shifted_nomask[img_id][my_trans_image]) << ", "
                                    << Ysize(exp_local_Fimgs_shifted_nomask[img_id][my_trans_image]) << ", "
                                    << Xsize(exp_local_Fimgs_shifted_nomask[img_id][my_trans_image]) << std::endl;
                            std::cerr << " mymodel.ori_size = " << mymodel.ori_size << std::endl;
                            Image<RFLOAT> tt;
                            FourierTransformer transformer;
                            tt().resize((mymodel.data_dim == 3) ? (mymodel.ori_size) : (1), mymodel.ori_size, mymodel.ori_size);
                            copMultidimArray<Complex> Faux = exp_local_Fimgs_shifted_nomask[img_id][my_trans_image];
                            copMultidimArray<Complex> Fo = windowFourierTransform(Faux, mymodel.ori_size);
                            tt() = transformer.inverseFourierTransform(Fo);
                            CenterFFT(tt(), -1);
                            img_save_nomask() += tt();
                            img_save_nomask.write("translational_searches_nomask_helix.spi");
                            std::cerr << " written translational_searches_nomask_helix.spi; press any key to continue..." << std::endl;
                            std::string str;
                            std::cin >> str;
                        }
                        #endif
                    }
                }
            }
        }
    }
    #ifdef TIMING
    if (part_id == mydata.sorted_idx[exp_my_first_part_id]) {
        timer.toc(
            do_also_unmasked ?              TIMING_ESP_PRECW :
            exp_current_oversampling == 0 ? TIMING_ESP_PREC1 :
                                            TIMING_ESP_PREC2
        );
    }
    #endif
}

bool MlOptimiser::isSignificantAnyImageAnyTranslation(
    long int iorient, int exp_itrans_min, int exp_itrans_max,
    MultidimArray<bool> &exp_Mcoarse_significant
) {

    long int exp_nr_trans = exp_itrans_max - exp_itrans_min + 1;
    for (long int ipart = 0; ipart < Ysize(exp_Mcoarse_significant); ipart++) {
        long int ihidden = iorient * exp_nr_trans;
        for (long int itrans = exp_itrans_min; itrans <= exp_itrans_max; itrans++, ihidden++) {
            #ifdef DEBUG_CHECKSIZES
            if (ihidden >= Xsize(exp_Mcoarse_significant)) {
                std::cerr << " ihidden= " << ihidden << " Xsize(exp_Mcoarse_significant)= " << Xsize(exp_Mcoarse_significant) << std::endl;
                std::cerr << " iorient= " << iorient << " itrans= " << itrans << " exp_nr_trans= " << exp_nr_trans << std::endl;
                REPORT_ERROR("ihidden > Xsize: ");
            }
            #endif
            if (direct::elem(exp_Mcoarse_significant, ipart, ihidden)) return true;
        }
    }
    return false;

}


void MlOptimiser::getAllSquaredDifferences(
    long int part_id, int ibody,
    int exp_ipass, int exp_current_oversampling, int metadata_offset,
    int exp_idir_min,   int exp_idir_max,   int exp_ipsi_min,   int exp_ipsi_max,
    int exp_itrans_min, int exp_itrans_max, int exp_iclass_min, int exp_iclass_max,
    std::vector<RFLOAT> &exp_min_diff2,
    std::vector<RFLOAT> &exp_highres_Xi2_img,
    std::vector<MultidimArray<Complex> > &exp_Fimg,
    std::vector<MultidimArray<RFLOAT> > &exp_Fctf,
    MultidimArray<RFLOAT> &exp_Mweight,
    MultidimArray<bool> &exp_Mcoarse_significant,
    std::vector<int> &exp_pointer_dir_nonzeroprior, std::vector<int> &exp_pointer_psi_nonzeroprior,
    std::vector<RFLOAT> &exp_directions_prior, std::vector<RFLOAT> &exp_psi_prior,
    std::vector<std::vector<MultidimArray<Complex> > > &exp_local_Fimgs_shifted,
    std::vector<MultidimArray<RFLOAT> > &exp_local_Minvsigma2,
    std::vector<MultidimArray<RFLOAT> > &exp_local_Fctf,
    std::vector<RFLOAT> &exp_local_sqrtXi2
) {

    #ifdef TIMING
    if (part_id == mydata.sorted_idx[exp_my_first_part_id])
        timer.tic(exp_ipass == 0 ? TIMING_ESP_DIFF1 : TIMING_ESP_DIFF2);
    #endif

    // #define DEBUG_GETALLDIFF2
    #ifdef DEBUG_GETALLDIFF2
    std::cerr << " ipass= " << exp_ipass << " exp_current_oversampling= " << exp_current_oversampling << std::endl;
    std::cerr << " sampling.NrPsiSamplings(exp_current_oversampling)= " << sampling.NrPsiSamplings(exp_current_oversampling) << std::endl;
    std::cerr << " sampling.NrTranslationalSamplings(exp_current_oversampling)= " << sampling.NrTranslationalSamplings(exp_current_oversampling) << std::endl;
    std::cerr << " sampling.NrSamplingPoints(exp_current_oversampling)= " << sampling.NrSamplingPoints(exp_current_oversampling) << std::endl;
    std::cerr << " sampling.oversamplingFactorOrientations(exp_current_oversampling)= "<<sampling.oversamplingFactorOrientations(exp_current_oversampling) << std::endl;
    std::cerr << " sampling.oversamplingFactorTranslations(exp_current_oversampling)= "<<sampling.oversamplingFactorTranslations(exp_current_oversampling) << std::endl;
    #endif

    // Initialise min_diff and exp_Mweight for this pass
    int exp_nr_images = mydata.numberOfImagesInParticle(part_id);
    long int exp_nr_dir = do_skip_align || do_skip_rotate ? 1 : sampling.NrDirections(0, &exp_pointer_dir_nonzeroprior);
    long int exp_nr_psi = do_skip_align || do_skip_rotate || do_only_sample_tilt ? 1 : sampling.NrPsiSamplings(0, &exp_pointer_psi_nonzeroprior);
    long int exp_nr_trans = do_skip_align ? 1 : sampling.NrTranslationalSamplings();
    long int exp_nr_oversampled_rot = sampling.oversamplingFactorOrientations(exp_current_oversampling);
    long int exp_nr_oversampled_trans = sampling.oversamplingFactorTranslations(exp_current_oversampling);

    exp_Mweight.resize(exp_nr_images, mymodel.nr_classes * exp_nr_dir * exp_nr_psi * exp_nr_trans * exp_nr_oversampled_rot * exp_nr_oversampled_trans);
    exp_Mweight = -999.0;
    if (exp_ipass == 0)
        exp_Mcoarse_significant.clear();

    exp_min_diff2.clear();
    exp_min_diff2.resize(exp_nr_images, LARGE_NUMBER);

    std::vector<MultidimArray<Complex> > dummy;
    std::vector<std::vector<MultidimArray<Complex> > > dummy2;

    precalculateShiftedImagesCtfsAndInvSigma2s(
        false, false, part_id, exp_current_oversampling, metadata_offset,
        exp_itrans_min, exp_itrans_max, exp_Fimg, dummy, exp_Fctf, exp_local_Fimgs_shifted, dummy2,
        exp_local_Fctf, exp_local_sqrtXi2, exp_local_Minvsigma2
    );

    // Loop only from exp_iclass_min to exp_iclass_max to deal with seed generation in first iteration
    for (int exp_iclass = exp_iclass_min; exp_iclass <= exp_iclass_max; exp_iclass++) {
        if (mymodel.pdf_class[exp_iclass] > 0.0) {
            // Local variables
            std::vector<RFLOAT> oversampled_rot, oversampled_tilt, oversampled_psi;
            std::vector<RFLOAT> oversampled_translations_x, oversampled_translations_y, oversampled_translations_z;
            RFLOAT *Minvsigma2;
            Matrix2D<RFLOAT> Abody, Aori;

            if (mymodel.nr_bodies > 1) {
                // ipart=0 because in multi-body refinement we do not do movie frames!
                RFLOAT rot_ori =  direct::elem(exp_metadata, metadata_offset, METADATA_ROT);
                RFLOAT tilt_ori = direct::elem(exp_metadata, metadata_offset, METADATA_TILT);
                RFLOAT psi_ori =  direct::elem(exp_metadata, metadata_offset, METADATA_PSI);
                Aori = Euler::angles2matrix(rot_ori, tilt_ori, psi_ori);
            }

            MultidimArray<Complex> Fref, Frefctf, Fimg_otfshift;
            Fref.resize(exp_local_Minvsigma2[0]);
            Frefctf.resize(exp_local_Minvsigma2[0]);
            if (do_shifts_onthefly)
                Fimg_otfshift.resize(Frefctf);

            for (long int idir = exp_idir_min, iorient = 0; idir <= exp_idir_max; idir++)
            for (long int ipsi = exp_ipsi_min; ipsi <= exp_ipsi_max; ipsi++, iorient++) {
                long int iorientclass = exp_iclass * exp_nr_dir * exp_nr_psi + iorient;

                // Get prior for this direction and skip calculation if prior==0
                RFLOAT pdf_orientation;
                if (do_skip_align || do_skip_rotate) {
                    #ifdef DEBUG_CHECKSIZES
                    if (exp_iclass >= mymodel.pdf_class.size()) {
                        std::cerr<< "exp_iclass= "<<exp_iclass<<" mymodel.pdf_class.size()= "<< mymodel.pdf_class.size() <<std::endl;
                        REPORT_ERROR("exp_iclass >= mymodel.pdf_class.size()");
                    }
                    #endif
                    pdf_orientation = mymodel.pdf_class[exp_iclass];
                } else if (mymodel.orientational_prior_mode == NOPRIOR) {
                    #ifdef DEBUG_CHECKSIZES
                    if (idir >= Xsize(mymodel.pdf_direction[exp_iclass])) {
                        std::cerr<< "idir= " << idir << " Xsize(mymodel.pdf_direction[exp_iclass])= " << Xsize(mymodel.pdf_direction[exp_iclass]) << std::endl;
                        REPORT_ERROR("idir >= mymodel.pdf_direction[exp_iclass].size()");
                    }
                    #endif
                    pdf_orientation = mymodel.pdf_direction[exp_iclass][idir];
                } else {
                    pdf_orientation = exp_directions_prior[idir] * exp_psi_prior[ipsi];
                }
                // In the first pass, always proceed
                // In the second pass, check whether one of the translations for this orientation had a significant weight in the first pass
                // if so, proceed with projecting the reference in that direction
                if ((exp_ipass == 0 || isSignificantAnyImageAnyTranslation(
                    iorientclass, exp_itrans_min, exp_itrans_max,
                    exp_Mcoarse_significant
                )) && pdf_orientation > 0.0) {
                    // Now get the oversampled (rot, tilt, psi) triplets
                    // This will be only the original (rot,tilt,psi) triplet in the first pass (exp_current_oversampling==0)
                    sampling.getOrientations(idir, ipsi, exp_current_oversampling, oversampled_rot, oversampled_tilt, oversampled_psi,
                            exp_pointer_dir_nonzeroprior, exp_directions_prior, exp_pointer_psi_nonzeroprior, exp_psi_prior);
                    // Loop over all oversampled orientations (only a single one in the first pass)
                    for (long int iover_rot = 0; iover_rot < exp_nr_oversampled_rot; iover_rot++) {
                        // loop over all images inside this particle
                        for (int img_id = 0; img_id < mydata.numberOfImagesInParticle(part_id); img_id++) {

                            int my_metadata_offset = metadata_offset + img_id;
                            RFLOAT my_pixel_size = mydata.getImagePixelSize(part_id, img_id);
                            int optics_group = mydata.getOpticsGroup(part_id, img_id);
                            bool ctf_premultiplied = mydata.obsModel.getCtfPremultiplied(optics_group);

                            // Get the Euler matrix
                            Matrix2D<RFLOAT> A = Euler::angles2matrix(
                                oversampled_rot[iover_rot],
                                oversampled_tilt[iover_rot],
                                oversampled_psi[iover_rot]
                            );

                            // Project the reference map (into Fref)
                            #ifdef TIMING
                            // Only time one thread, as I also only time one MPI process
                            if (part_id == mydata.sorted_idx[exp_my_first_part_id])
                                timer.tic(TIMING_DIFF_PROJ);
                            #endif

                            // For multi-body refinements, A are only 'residual' orientations, Abody is the complete Euler matrix
                            if (mymodel.nr_bodies > 1) {
                                Abody =  Aori * (mymodel.orient_bodies[ibody]).transpose() * A_rot90 * A * mymodel.orient_bodies[ibody];
                                if (mydata.obsModel.hasMagMatrices) { Abody *= mydata.obsModel.anisoMag(optics_group); }
                                Abody *= mydata.obsModel.scaleDifference(optics_group, mymodel.ori_size, mymodel.pixel_size);
                                
                                Fref = mymodel.PPref[ibody].get2DFourierTransform(
                                    exp_local_Minvsigma2[0].xdim, exp_local_Minvsigma2[0].ydim, exp_local_Minvsigma2[0].zdim, Abody);
                            } else {
                                if (mydata.obsModel.hasMagMatrices) { A *= mydata.obsModel.anisoMag(optics_group); }
                                A *= mydata.obsModel.scaleDifference(optics_group, mymodel.ori_size, mymodel.pixel_size);
                                Fref = mymodel.PPref[exp_iclass].get2DFourierTransform(
                                    exp_local_Minvsigma2[0].xdim, exp_local_Minvsigma2[0].ydim, exp_local_Minvsigma2[0].zdim, A);
                            }

                            #ifdef TIMING
                            // Only time one thread, as I also only time one MPI process
                            if (part_id == mydata.sorted_idx[exp_my_first_part_id])
                                timer.toc(TIMING_DIFF_PROJ);
                            #endif

                            Minvsigma2 = exp_local_Minvsigma2[img_id].data;

                            // Apply CTF to reference projection
                            if (do_ctf_correction && refs_are_ctf_corrected) {
                                if (ctf_premultiplied) {
                                    // TODO: ignore CTF until first peak of premultiplied CTF?
                                    for (long int n = 0; n < Fref.size(); n++) {
                                        Frefctf[n] = Fref[n] * exp_local_Fctf[img_id][n] * exp_local_Fctf[img_id][n];
                                    }
                                } else {
                                    for (long int n = 0; n < Fref.size(); n++) {
                                        Frefctf[n] = Fref[n] * exp_local_Fctf[img_id][n];
                                    }
                                }
                            } else {
                                Frefctf = Fref;
                            }

                            if (do_scale_correction) {
                                int group_id = mydata.getGroupId(part_id, img_id);
                                RFLOAT myscale = mymodel.scale_correction[group_id];
                                for (long int n = 0; n < (Frefctf).size(); n++) {
                                    Frefctf[n] *= myscale;
                                }
                            }

                            long int ihidden = iorientclass * exp_nr_trans;
                            for (long int itrans = exp_itrans_min; itrans <= exp_itrans_max; itrans++, ihidden++) {
                                #ifdef DEBUG_CHECKSIZES
                                if (exp_ipass > 0 && ihidden >= Xsize(exp_Mcoarse_significant)) {
                                    std::cerr << "ihidden= " << ihidden << " Xsize(exp_Mcoarse_significant)= " << Xsize(exp_Mcoarse_significant) << std::endl;
                                    REPORT_ERROR("ihidden >= Xsize(exp_Mcoarse_significant)");
                                }
                                #endif
                                // In the first pass, always proceed
                                // In the second pass, check whether this translations (&orientation) had a significant weight in the first pass
                                if (exp_ipass == 0 || direct::elem(exp_Mcoarse_significant, img_id, ihidden)) {
                                    // 1 Jun 2015 - Shaoda & Sjors, Helical refinement
                                    sampling.getTranslationsInPixel(
                                        itrans, exp_current_oversampling, my_pixel_size,
                                        oversampled_translations_x, oversampled_translations_y, oversampled_translations_z,
                                        do_helical_refine && !ignore_helical_symmetry
                                    );
                                    for (long int iover_trans = 0; iover_trans < exp_nr_oversampled_trans; iover_trans++) {
                                        #ifdef TIMING
                                        // Only time one thread, as I also only time one MPI process
                                        if (part_id == mydata.sorted_idx[exp_my_first_part_id])
                                            timer.tic(TIMING_DIFF2_GETSHIFT);
                                        #endif
                                        /// Now get the shifted image
                                        // Use a pointer to avoid copying the entire array again in this highly expensive loop
                                        Complex *Fimg_shift = get_shifted_image(
                                            img_id, exp_nr_oversampled_trans,
                                            exp_nr_images, itrans, exp_itrans_min, iover_trans,
                                            exp_local_Fimgs_shifted,
                                            global_fftshifts_ab_coarse[optics_group],
                                            global_fftshifts_ab_current[optics_group],
                                            global_fftshifts_ab2_coarse[optics_group],
                                            global_fftshifts_ab2_current[optics_group],
                                            do_shifts_onthefly, do_skip_align,
                                            do_helical_refine, ignore_helical_symmetry,
                                            oversampled_translations_x,
                                            oversampled_translations_y,
                                            oversampled_translations_z,
                                            mymodel.data_dim, mymodel.ori_size,
                                            exp_current_oversampling,
                                            exp_metadata,
                                            my_metadata_offset,
                                            strict_highres_exp,
                                            image_current_size[optics_group],
                                            image_coarse_size[optics_group],
                                            Ysize(Frefctf),
                                            Fimg_otfshift,
                                            tab_sin, tab_cos
                                        );
                                        #ifdef TIMING
                                        // Only time one thread, as I also only time one MPI process
                                        if (part_id == mydata.sorted_idx[exp_my_first_part_id])
                                            timer.toc(TIMING_DIFF2_GETSHIFT);
                                        #endif

                                        #ifdef TIMING
                                        // Only time one thread, as I also only time one MPI process
                                        if (part_id == mydata.sorted_idx[exp_my_first_part_id])
                                            timer.tic(TIMING_DIFF_DIFF2);
                                        #endif
                                        RFLOAT diff2;
                                        if (do_cc()) {
                                            // Do not calculate squared-differences, but signal product
                                            // Negative values because smaller is worse in this case
                                            diff2 = 0.0;
                                            RFLOAT suma2 = 0.0;
                                            for (long int n = 0; n < Frefctf.size(); n++) {
                                                Complex A = Frefctf[n];
                                                Complex X = Fimg_shift[n];
                                                diff2 -= A.real * X.real + A.imag * X.imag;  // A dot X
                                                suma2 += norm(A);
                                            }
                                            // Normalised cross-correlation coefficient: divide by power of reference (power of image is a constant)
                                            diff2 /= sqrt(suma2) * exp_local_sqrtXi2[img_id];
                                        } else {
                                            // Calculate the actual squared difference term of the Gaussian probability function
                                            // If current_size < mymodel.ori_size diff2 is initialised to the sum of
                                            // all |Xij|2 terms that lie between current_size and ori_size
                                            // Factor two because of factor 2 in division below, NOT because of 2-dimensionality of the complex plane!
                                            diff2 = exp_highres_Xi2_img[img_id] / 2.0;
                                            for (long int n = 0; n < (Frefctf).size(); n++) {
                                                Complex diff = Frefctf[n] - Fimg_shift[n];
                                                diff2 += (diff.real * diff.real + diff.imag * diff.imag) * 0.5 * Minvsigma2[n];
                                            }
                                        }
                                        #ifdef TIMING
                                        // Only time one thread, as I also only time one MPI process
                                        if (part_id == mydata.sorted_idx[exp_my_first_part_id])
                                            timer.toc(TIMING_DIFF_DIFF2);
                                        #endif

                                        // Store all diff2 in exp_Mweight
                                        long int ihidden_over = sampling.getPositionOversampledSamplingPoint(
                                            ihidden, exp_current_oversampling, iover_rot, iover_trans
                                        );
                                        //#define DEBUG_GETALLDIFF2
                                        #ifdef DEBUG_GETALLDIFF2
                                        pthread_mutex_lock(&global_mutex);
                                        if (itrans == exp_itrans_min && iover_trans == 0 && ipsi == exp_ipsi_min) {
                                        //if (ibody==1 && part_id == 0 && exp_ipass==0 && ihidden_over == 40217)
                                            //std::cerr << " iover_rot= "<<iover_rot << "exp_nr_oversampled_rot= " << exp_nr_oversampled_rot << " oversampled_rot[iover_rot]= " << oversampled_rot[iover_rot]
                                            //		  << " oversampled_tilt[iover_rot]= " << oversampled_tilt[iover_rot]
                                            //	      << " oversampled_psi[iover_rot]= " <<  oversampled_psi[iover_rot];
                                            angles_t angles = Euler::matrix2angles(A);
                                            std::cerr << " ihidden_over= " << ihidden_over << " diff2= " << diff2
                                                    << " rot= " << angles.rot
                                                    << " tilt= " << angles.tilt
                                                    << " psi= " << angles.psi
                                                    // non-oversampling correct only!!
                                                    << " x= " << oversampled_translations_x[0] << " y=" << oversampled_translations_y[0];
                                            // std::cerr << " A= " << A << std::endl;
                                            // angles_t angles = Euler::matrix2angles(Abody);
                                            // std::cerr << " Brot= " << angles.rot
                                            // 		<< " Btilt= " << angles.tilt
                                            // 		<< " Bpsi= " << angles.psi << std::endl;

                                            FourierTransformer transformer;
                                            MultidimArray<Complex> Fish;
                                            Fish.resize(exp_local_Minvsigma2[img_id]);
                                            for (long int n = 0; n < Fish.size(); n++) {
                                                Fish[n] = Fimg_shift[n];
                                            }
                                            Image<RFLOAT> tt;
                                            int exp_current_image_size;
                                            if (strict_highres_exp > 0.0)
                                                // Use smaller images in both passes and keep a maximum on coarse_size, just like in FREALIGN
                                                exp_current_image_size = image_coarse_size[optics_group];
                                            else if (adaptive_oversampling > 0)
                                                // Use smaller images in the first pass, larger ones in the second pass
                                                exp_current_image_size = (exp_current_oversampling == 0) ? image_coarse_size[optics_group] : image_current_size[optics_group];
                                            else
                                                exp_current_image_size = image_current_size[optics_group];
                                            if (mymodel.data_dim == 3)
                                                tt().resize(exp_current_image_size, exp_current_image_size, exp_current_image_size);
                                            else
                                                tt().resize(exp_current_image_size, exp_current_image_size);
                                            tt() = transformer.inverseFourierTransform(Fish);
                                            CenterFFT(tt(), -1);
                                            FileName fnt = "Fimg.spi";
                                            // fnt = FileName::compose("Fimg_shift1_i", ihidden_over, "spi");
                                            tt.write(fnt);

                                            tt() = transformer.inverseFourierTransform(Frefctf);
                                            CenterFFT(tt(), -1);
                                            fnt = "Fref.spi";
                                            // fnt = FileName::compose("Fref1_i", ihidden_over, "spi");
                                            tt.write(fnt);

                                            //for (int i = 0; i< mymodel.scale_correction.size(); i++)
                                            //	std::cerr << i << " scale="<<mymodel.scale_correction[i]<<std::endl;
                                            int group_id = mydata.getGroupId(part_id, img_id);
                                            RFLOAT myscale = mymodel.scale_correction[group_id];
                                            // std::cerr << " oversampled_rot[iover_rot]= " << oversampled_rot[iover_rot] << " oversampled_tilt[iover_rot]= " << oversampled_tilt[iover_rot] << " oversampled_psi[iover_rot]= " << oversampled_psi[iover_rot] << std::endl;
                                            // std::cerr << " group_id= " << group_id << " myscale= " << myscale <<std::endl;
                                            std::cerr << " itrans= "  << itrans  << " itrans * exp_nr_oversampled_trans +  iover_trans= " << itrans * exp_nr_oversampled_trans +  iover_trans << " ihidden= " << ihidden << std::endl;
                                            std::cerr << " part_id= " << part_id << " name= " << mydata.particles[part_id].name << std::endl;
                                            std::cerr << " img_id= "  << img_id  << " name= " << mydata.particles[part_id].images[img_id].name << std::endl;

                                            // std::cerr << " myrank= "<< myrank<<std::endl;
                                            // std::cerr << "Written Fimg_shift.spi and Fref.spi. Press any key to continue... part_id= " << part_id<< std::endl;
                                            char c;
                                            std::cin >> c;
                                            // exit(0);
                                        }
                                        pthread_mutex_unlock(&global_mutex);

                                        #endif
                                        //#define DEBUG_DIFF2_ISNAN
                                        #ifdef DEBUG_DIFF2_ISNAN
                                        if (std::isnan(diff2))
                                        {
                                            pthread_mutex_lock(&global_mutex);
                                            std::cerr <<" img_id= "<<img_id<<" name= "<< mydata.particles[part_id].images[img_id].name << std::endl;
                                            std::cerr << " exp_iclass= " << exp_iclass << std::endl;
                                            std::cerr << " diff2= " << diff2 << std::endl;
                                            std::cerr << " exp_highres_Xi2_img[img_id]= " << exp_highres_Xi2_img[img_id] << std::endl;
                                            std::cerr<< " exp_nr_oversampled_trans="<<exp_nr_oversampled_trans<<std::endl;
                                            std::cerr<< " exp_nr_oversampled_rot="<<exp_nr_oversampled_rot<<std::endl;
                                            std::cerr << " iover_rot= " << iover_rot << " iover_trans= " << iover_trans << " ihidden= " << ihidden << std::endl;
                                            std::cerr << " exp_current_oversampling= " << exp_current_oversampling << std::endl;
                                            std::cerr << " ihidden_over= " << ihidden_over << " Xsize(Mweight)= " << Xsize(exp_Mweight) << std::endl;
                                            std::cerr << " (mymodel.PPref[exp_iclass]).ori_size= " << (mymodel.PPref[exp_iclass]).ori_size << " (mymodel.PPref[exp_iclass]).r_max= " << (mymodel.PPref[exp_iclass]).r_max << std::endl;
                                            int group_id = mydata.getGroupId(part_id, img_id);
                                            std::cerr << " mymodel.scale_correction[group_id]= " << mymodel.scale_correction[group_id] << std::endl;
                                            if (std::isnan(mymodel.scale_correction[group_id]))
                                            {
                                                for (int i=0; i < mymodel.scale_correction.size(); i++)
                                                    std::cerr << " i= " << i << " mymodel.scale_correction[i]= " << mymodel.scale_correction[i] << std::endl;
                                            }
                                            std::cerr << " group_id= " << group_id << std::endl;
                                            Image<RFLOAT> It;
                                            std::cerr << "Frefctf shape= "; Frefctf.printShape(std::cerr);
                                            MultidimArray<Complex> Fish;
                                            Fish.resize(exp_local_Minvsigma2[img_id]);
                                            for (long int n = 0; n < Fish.size(); n++) {
                                                Fish[n] = Fimg_shift[n];
                                            }
                                            std::cerr << "Fimg_shift shape= "; (Fish).printShape(std::cerr);
                                            It()=exp_local_Fctf[img_id];
                                            It.write("exp_local_Fctf.spi");
                                            std::cerr << "written exp_local_Fctf.spi" << std::endl;
                                            FourierTransformer transformer;
                                            Image<RFLOAT> tt;
                                            int exp_current_image_size;
                                            if (strict_highres_exp > 0.0)
                                                // Use smaller images in both passes and keep a maximum on coarse_size, just like in FREALIGN
                                                exp_current_image_size = image_coarse_size[optics_group];
                                            else if (adaptive_oversampling > 0)
                                                // Use smaller images in the first pass, larger ones in the second pass
                                                exp_current_image_size = (exp_current_oversampling == 0) ? image_coarse_size[optics_group] : image_current_size[optics_group];
                                            else
                                                exp_current_image_size = image_current_size[optics_group];
                                            tt().resize(exp_current_image_size, exp_current_image_size);
                                            tt() = transformer.inverseFourierTransform(Fish);
                                            CenterFFT(tt(), -1);
                                            tt.write("Fimg_shift.spi");
                                            std::cerr << "written Fimg_shift.spi" << std::endl;
                                            FourierTransformer transformer2;
                                            tt().initZeros();
                                            tt() = transformer2.inverseFourierTransform(Frefctf);
                                            CenterFFT(tt(), -1);
                                            tt.write("Frefctf.spi");
                                            std::cerr << "written Frefctf.spi" << std::endl;
                                            FourierTransformer transformer3;
                                            tt().initZeros();
                                            tt() = transformer3.inverseFourierTransform(Fref);
                                            CenterFFT(tt(), -1);
                                            tt.write("Fref.spi");
                                            std::cerr << "written Fref.spi" << std::endl;
                                            std::cerr << " A= " << A << std::endl;
                                            std::cerr << "written Frefctf.spi" << std::endl;

                                            std::cerr << " exp_iclass= " << exp_iclass << std::endl;
                                            Fref = mymodel.PPref[exp_iclass].get2DFourierTransform(
                                                exp_local_Minvsigma2[img_id].xdim, exp_local_Minvsigma2[img_id].ydim, exp_local_Minvsigma2[img_id].zdim, A);
                                            tt() = transformer3.inverseFourierTransform(Fref);
                                            CenterFFT(tt(), -1);
                                            tt.write("Fref2.spi");
                                            std::cerr << "written Fref2.spi" << std::endl;
                                            Image<RFLOAT> Itt;
                                            Itt().resize(Xsize(mymodel.PPref[exp_iclass].data), Ysize(mymodel.PPref[exp_iclass].data), Zsize(mymodel.PPref[exp_iclass].data));
                                            for (long int n = 0; n < (Itt()).size(); n++) {
                                                Itt()[n] = abs(mymodel.PPref[exp_iclass].data[n]);
                                            }
                                            Itt.write("PPref_data.spi");
                                            REPORT_ERROR("diff2 is not a number");
                                            pthread_mutex_unlock(&global_mutex);
                                            exit(0);
                                        }
                                        #endif
                                        //#define DEBUG_VERBOSE
                                        #ifdef DEBUG_VERBOSE
                                        pthread_mutex_lock(&global_mutex);
                                        std::cout <<" name= "<< mydata.particles[part_id].images[img_id].name << " rot= " << oversampled_rot[iover_rot] << " tilt= "<< oversampled_tilt[iover_rot] << " psi= " << oversampled_psi[iover_rot] << std::endl;
                                        std::cout <<" name= "<< mydata.particles[part_id].images[img_id].name << " ihidden_over= " << ihidden_over << " diff2= " << diff2 << " exp_min_diff2= " << exp_min_diff2 << std::endl;
                                        pthread_mutex_unlock(&global_mutex);
                                        #endif
                                        #ifdef DEBUG_CHECKSIZES
                                        if (ihidden_over >= Xsize(exp_Mweight) )
                                        {
                                            std::cerr<< " exp_nr_oversampled_trans="<<exp_nr_oversampled_trans<<std::endl;
                                            std::cerr<< " exp_nr_oversampled_rot="<<exp_nr_oversampled_rot<<std::endl;
                                            std::cerr << " iover_rot= " << iover_rot << " iover_trans= " << iover_trans << " ihidden= " << ihidden << std::endl;
                                            std::cerr << " exp_current_oversampling= " << exp_current_oversampling << std::endl;
                                            std::cerr << " exp_itrans_min= " << exp_itrans_min <<" exp_nr_trans= " << exp_nr_trans << std::endl;
                                            std::cerr << " exp_itrans_max= " << exp_itrans_max << " iorientclass= " << iorientclass << " itrans= " << itrans << std::endl;
                                            std::cerr << " exp_nr_dir= " << exp_nr_dir << " exp_idir_min= " << exp_idir_min << " exp_idir_max= " << exp_idir_max << std::endl;
                                            std::cerr << " exp_nr_psi= " << exp_nr_psi << " exp_ipsi_min= " << exp_ipsi_min << " exp_ipsi_max= " << exp_ipsi_max << std::endl;
                                            std::cerr << " exp_iclass= " << exp_iclass << std::endl;
                                            std::cerr << " iorient= " << iorient << std::endl;
                                            std::cerr << " ihidden_over= " << ihidden_over << " Xsize(Mweight)= " << Xsize(exp_Mweight) << std::endl;
                                            REPORT_ERROR("ihidden_over >= Xsize(Mweight)");
                                        }
                                        #endif
                                        direct::elem(exp_Mweight, img_id, ihidden_over) = diff2;

                                        // Keep track of minimum of all diff2, only for the last image in this series
                                        if (diff2 < exp_min_diff2[img_id]) {
                                            exp_min_diff2[img_id] = diff2;
                                            /*
                                            if (part_id == 0)
                                            {
                                                std::cerr << " part_id= " << part_id << " ihidden_over= " << ihidden_over << " diff2= " << diff2
                                                << " x= " << oversampled_translations_x[iover_trans] << " y=" <<oversampled_translations_y[iover_trans]
                                                << " iover_trans= "<<iover_trans << "Xi2= " << exp_highres_Xi2_img[img_id] << " Minv_sigma2= " << exp_local_Minvsigma2[img_id][10]
                                                << " xsize= " << Xsize(Frefctf)
                                                << " Frefctf= " << Frefctf[10].real
                                                << " Fimgshift= " << (*(Fimg_shift + 10)).real
                                                << std::endl;
                                                }
                                                */
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    #ifdef TIMING
    if (part_id == mydata.sorted_idx[exp_my_first_part_id]) {
        timer.toc(exp_ipass == 0 ? TIMING_ESP_DIFF1 : TIMING_ESP_DIFF2);
    }
    #endif
}


void MlOptimiser::convertAllSquaredDifferencesToWeights(
    long int part_id, int ibody, int exp_ipass,
    int exp_current_oversampling, int metadata_offset,
    int exp_idir_min,   int exp_idir_max,   int exp_ipsi_min, int exp_ipsi_max,
    int exp_itrans_min, int exp_itrans_max, int exp_iclass_min, int exp_iclass_max,
    MultidimArray<RFLOAT> &exp_Mweight, MultidimArray<bool> &exp_Mcoarse_significant,
    std::vector<RFLOAT> &exp_significant_weight, std::vector<RFLOAT> &exp_sum_weight,
    std::vector<Matrix1D<RFLOAT> > &exp_old_offset, std::vector<Matrix1D<RFLOAT> > &exp_prior, std::vector<RFLOAT> &exp_min_diff2,
    std::vector<int> &exp_pointer_dir_nonzeroprior, std::vector<int> &exp_pointer_psi_nonzeroprior,
    std::vector<RFLOAT> &exp_directions_prior, std::vector<RFLOAT> &exp_psi_prior
) {

    #ifdef TIMING
    if (part_id == mydata.sorted_idx[exp_my_first_part_id]) {
        timer.tic(exp_ipass == 0 ? TIMING_ESP_WEIGHT1 : TIMING_ESP_WEIGHT2);
    }
    #endif

    // Convert the squared differences into weights
    // Note there is only one weight for each part_id, because a whole series of images is treated as one particle

    long int exp_nr_dir = (do_skip_align || do_skip_rotate) ? 1 : sampling.NrDirections(0, &exp_pointer_dir_nonzeroprior);
    long int exp_nr_psi = (do_skip_align || do_skip_rotate || do_only_sample_tilt) ? 1 : sampling.NrPsiSamplings(0, &exp_pointer_psi_nonzeroprior);
    long int exp_nr_trans = (do_skip_align) ? 1 : sampling.NrTranslationalSamplings();
    int exp_nr_images = mydata.numberOfImagesInParticle(part_id);
    long int exp_nr_oversampled_rot = sampling.oversamplingFactorOrientations(exp_current_oversampling);
    long int exp_nr_oversampled_trans = sampling.oversamplingFactorTranslations(exp_current_oversampling);

    // Initialise
    exp_sum_weight.clear();
    exp_sum_weight.resize(exp_nr_images, 0.0);

    RFLOAT my_sigma2_offset = mymodel.nr_bodies > 1 ?
        mymodel.sigma_offset_bodies[ibody] * mymodel.sigma_offset_bodies[ibody] :
        mymodel.sigma2_offset;

    // #define DEBUG_CONVERTDIFF2W
    #ifdef DEBUG_CONVERTDIFF2W
    RFLOAT max_weight = -1.0;
    RFLOAT opt_psi, opt_xoff, opt_yoff;
    int opt_iover_rot, opt_iover_trans, opt_ipsi, opt_itrans;
    long int opt_ihidden, opt_ihidden_over;
    #endif

    // loop over all images inside this particle
    for (int img_id = 0; img_id < exp_nr_images; img_id++) {
        int my_metadata_offset = metadata_offset + img_id;
        RFLOAT my_pixel_size = mydata.getImagePixelSize(part_id, img_id);

        RFLOAT exp_thisimage_sumweight = 0.0;
        RFLOAT old_offset_x, old_offset_y, old_offset_z;
        if (mymodel.nr_bodies > 1) {
            old_offset_x = old_offset_y = old_offset_z = 0.0;
        } else {
            old_offset_x = exp_old_offset[img_id][0];
            old_offset_y = exp_old_offset[img_id][1];
            if (mymodel.data_dim == 3)
            old_offset_z = exp_old_offset[img_id][2];
        }

        if (do_cc()) {
            // Binarize the squared differences array to skip marginalisation
            RFLOAT mymindiff2 = 99.e10;
            long int myminidx = -1;
            // Find the smallest element in this row of exp_Mweight
            for (long int i = 0; i < Xsize(exp_Mweight); i++) {

                RFLOAT cc = direct::elem(exp_Mweight, img_id, i);
                // ignore non-determined cc
                if (cc == -999.0)
                    continue;

                // just search for the maximum
                if (cc < mymindiff2) {
                    mymindiff2 = cc;
                    myminidx = i;
                }
            }
            // Set all except for the best hidden variable to zero and the smallest element to 1
            for (long int i = 0; i < Xsize(exp_Mweight); i++)
                direct::elem(exp_Mweight, i, img_id) = 0.0;

            direct::elem(exp_Mweight, myminidx, img_id) = 1.0;
            exp_thisimage_sumweight += 1.0;

        } else {
            // Extra normalization
            RFLOAT pdf_orientation_mean(0), pdf_offset_mean(0);
            unsigned long pdf_orientation_count(0), pdf_offset_count(0);
            for (int exp_iclass = exp_iclass_min; exp_iclass <= exp_iclass_max; exp_iclass++) {

                for (long int idir = exp_idir_min, iorient = 0; idir <= exp_idir_max; idir++)
                for (long int ipsi = exp_ipsi_min; ipsi <= exp_ipsi_max; ipsi++, iorient++) {
                    pdf_orientation_mean +=
                        do_skip_align || do_skip_rotate ? mymodel.pdf_class[exp_iclass] :
                        mymodel.orientational_prior_mode == NOPRIOR ? mymodel.pdf_direction[exp_iclass][idir] :
                        exp_directions_prior[idir] * exp_psi_prior[ipsi];
                    pdf_orientation_count++;
                }

                RFLOAT myprior_x, myprior_y, myprior_z;
                if (mymodel.nr_bodies > 1) {
                    myprior_x = myprior_y = myprior_z = 0.0;
                }
                else if (mymodel.ref_dim == 2 && !do_helical_refine) {
                    myprior_x = mymodel.prior_offset_class[exp_iclass][0];
                    myprior_y = mymodel.prior_offset_class[exp_iclass][1];
                } else {
                    myprior_x = exp_prior[img_id][0];
                    myprior_y = exp_prior[img_id][1];
                    if (mymodel.data_dim == 3)
                    myprior_z = exp_prior[img_id][2];
                }
                for (long int itrans = exp_itrans_min; itrans <= exp_itrans_max; itrans++) {
                    RFLOAT offset_x = old_offset_x + sampling.translations_x[itrans];
                    RFLOAT offset_y = old_offset_y + sampling.translations_y[itrans];
                    RFLOAT tdiff2 = 0.0;
                    if (!do_helical_refine || ignore_helical_symmetry || mymodel.data_dim == 3)
                        tdiff2 += (offset_x - myprior_x) * (offset_x - myprior_x);
                    tdiff2 += (offset_y - myprior_y) * (offset_y - myprior_y);
                    if (mymodel.data_dim == 3) {
                        RFLOAT offset_z = old_offset_z + sampling.translations_z[itrans];
                        if (!do_helical_refine || ignore_helical_symmetry)
                            tdiff2 += (offset_z - myprior_z) * (offset_z - myprior_z);
                    }
                    // As of version 3.1, sigma_offsets are in Angstroms!
                    tdiff2 *= my_pixel_size * my_pixel_size;

                    // P(offset|sigma2_offset)
                    // This is the probability of the offset, given the model offset and variance.
                    RFLOAT pdf_offset;
                    pdf_offset_mean += my_sigma2_offset < 0.0001 ? (tdiff2 > 0.0 ? 0.0 : 1.0) : exp(tdiff2 / (-2.0 * my_sigma2_offset)) / (2.0 * PI * my_sigma2_offset);
                    pdf_offset_count++;
                }
            }
            pdf_orientation_mean /= (RFLOAT) pdf_orientation_count;
            pdf_offset_mean      /= (RFLOAT) pdf_offset_count;
            // Loop from iclass_min to iclass_max to deal with seed generation in first iteration
            for (int exp_iclass = exp_iclass_min; exp_iclass <= exp_iclass_max; exp_iclass++) {

                // Make PdfOffset calculation much faster...
                RFLOAT myprior_x, myprior_y, myprior_z;
                if (mymodel.nr_bodies > 1) {
                    myprior_x = myprior_y = myprior_z = 0.0;
                } else if (mymodel.ref_dim == 2) {
                    myprior_x = mymodel.prior_offset_class[exp_iclass][0];
                    myprior_y = mymodel.prior_offset_class[exp_iclass][1];
                } else {
                    myprior_x = exp_prior[img_id][0];
                    myprior_y = exp_prior[img_id][1];
                    if (mymodel.data_dim == 3)
                    myprior_z = exp_prior[img_id][2];
                }
                for (long int idir = exp_idir_min, iorient = 0; idir <= exp_idir_max; idir++) {
                    for (long int ipsi = exp_ipsi_min; ipsi <= exp_ipsi_max; ipsi++, iorient++) {
                        long int iorientclass = exp_iclass * exp_nr_dir * exp_nr_psi + iorient;
                        RFLOAT pdf_orientation;

                        // Get prior for this direction
                        if (do_skip_align || do_skip_rotate) {
                            pdf_orientation = mymodel.pdf_class[exp_iclass];
                        } else if (mymodel.orientational_prior_mode == NOPRIOR) {
                            pdf_orientation = mymodel.pdf_direction[exp_iclass][idir];
                        } else {
                            // P(orientation) = P(idir|dir_prior) * P(ipsi|psi_prior)
                            // This is the probability of the orientation, given the gathered
                            // statistics of all assigned orientations of the dataset, since we
                            // are assigning a gaussian prior to all parameters.
                            pdf_orientation = exp_directions_prior[idir] * exp_psi_prior[ipsi];
                        }

                        if (pdf_orientation_mean != 0.0)
                            pdf_orientation /= pdf_orientation_mean;

                        // Loop over all translations
                        long int ihidden = iorientclass * exp_nr_trans;
                        for (long int itrans = exp_itrans_min; itrans <= exp_itrans_max; itrans++, ihidden++) {
                            // 18 May 2015 - Shaoda & Sjors - Helical refinement (translational searches)
                            // Calculate the vector length of myprior
                            RFLOAT mypriors_len2 = myprior_x * myprior_x + myprior_y * myprior_y;
                            if (mymodel.data_dim == 3)
                                mypriors_len2 += myprior_z * myprior_z;
                            // If it is doing helical refinement AND Cartesian vector myprior has a length > 0, transform the vector to its helical coordinates
                            if (do_helical_refine && !ignore_helical_symmetry && mypriors_len2 > 0.00001) {
                                RFLOAT rot_deg  = direct::elem(exp_metadata, my_metadata_offset, METADATA_ROT);
                                RFLOAT tilt_deg = direct::elem(exp_metadata, my_metadata_offset, METADATA_TILT);
                                RFLOAT psi_deg  = direct::elem(exp_metadata, my_metadata_offset, METADATA_PSI);
                                transformCartesianAndHelicalCoords(
                                    myprior_x, myprior_y, myprior_z,
                                    rot_deg, tilt_deg, psi_deg,
                                    mymodel.data_dim, CART_TO_HELICAL_COORDS
                                );
                            }
                            // (For helical refinement) Now offset, old_offset, sampling.translations and myprior are all in helical coordinates

                            // To speed things up, only calculate pdf_offset at the coarse sampling.
                            // That should not matter much, and that way one does not need to calculate all the OversampledTranslations
                            RFLOAT offset_x = old_offset_x + sampling.translations_x[itrans];
                            RFLOAT offset_y = old_offset_y + sampling.translations_y[itrans];
                            RFLOAT tdiff2 = 0.0;
                            if (!do_helical_refine || ignore_helical_symmetry || mymodel.data_dim == 3)
                                tdiff2 += (offset_x - myprior_x) * (offset_x - myprior_x);
                            tdiff2 += (offset_y - myprior_y) * (offset_y - myprior_y);
                            if (mymodel.data_dim == 3) {
                                RFLOAT offset_z = old_offset_z + sampling.translations_z[itrans];
                                if (!do_helical_refine || ignore_helical_symmetry)
                                    tdiff2 += (offset_z - myprior_z) * (offset_z - myprior_z);
                            }

                            // As of version 3.1, sigma_offsets are in Angstroms!
                            tdiff2 *= my_pixel_size * my_pixel_size;

                            // P(offset|sigma2_offset)
                            // This is the probability of the offset, given the model offset and variance.
                            RFLOAT pdf_offset = my_sigma2_offset < 0.0001 ? (tdiff2 > 0.0 ? 0.0 : 1.0) : exp(tdiff2 / (-2.0 * my_sigma2_offset)) / (2.0 * PI * my_sigma2_offset);

                            if (pdf_offset_mean > 0.0) { pdf_offset /= pdf_offset_mean; }

                            #ifdef TIMING
                            // Only time one thread, as I also only time one MPI process
                            if (part_id == mydata.sorted_idx[exp_my_first_part_id])
                                timer.tic(TIMING_WEIGHT_EXP);
                            #endif
                            // Now first loop over iover_rot, because that is the order in exp_Mweight as well
                            long int ihidden_over = ihidden * exp_nr_oversampled_rot * exp_nr_oversampled_trans;
                            for (long int iover_rot = 0; iover_rot < exp_nr_oversampled_rot; iover_rot++) {
                                // Then loop over iover_trans
                                for (long int iover_trans = 0; iover_trans < exp_nr_oversampled_trans; iover_trans++, ihidden_over++) {
                                    // Only exponentiate for determined values of exp_Mweight
                                    // (this is always true in the first pass, but not so in the second pass)
                                    // Only deal with this sampling point if its weight was significant
                                    if (direct::elem(exp_Mweight, img_id, ihidden_over) < 0.0) {
                                        direct::elem(exp_Mweight, img_id, ihidden_over) = 0.0;
                                    } else {
                                        // Set the weight base to the probability of the parameters given the prior
                                        RFLOAT weight = pdf_orientation * pdf_offset;
                                        RFLOAT diff2 = direct::elem(exp_Mweight, img_id, ihidden_over) - exp_min_diff2[img_id];
                                        // next line because of numerical precision of exp-function
                                        #ifdef RELION_SINGLE_PRECISION
                                        if (diff2 > 88.0)
                                            weight = 0.0;
                                        #else
                                        if (diff2 > 700.0)
                                            weight = 0.0;
                                        #endif
                                        // TODO: use tabulated exp function?
                                        else weight *= exp(-diff2);

                                        //std::cerr << "ihidden_over= "<<ihidden_over << " weight= " << weight << " diff2= " << diff2
                                        //		<< " pdf_orientation= " << pdf_orientation << " pdf_offset= " << pdf_offset<< std::endl;
                                        //#define DEBUG_PSIANGLE_PDISTRIBUTION
                                        #ifdef DEBUG_PSIANGLE_PDISTRIBUTION
                                        std::cout << ipsi * 360.0 / sampling.NrPsiSamplings() << " " << weight << std::endl;
                                        #endif
                                        // Store the weight
                                        direct::elem(exp_Mweight, img_id, ihidden_over) = weight;
                                        #ifdef DEBUG_CHECKSIZES
                                        if (std::isnan(weight)) {
                                            pthread_mutex_lock(&global_mutex);
                                            std::cerr<< "weight= "<<weight<<" is not a number! " <<std::endl;
                                            std::cerr << " exp_min_diff2= " << exp_min_diff2 << std::endl;
                                            std::cerr << " part_id= " << part_id << " img_id= "<< img_id << std::endl;
                                            std::cerr << " direct::elem(exp_Mweight, img_id, ihidden_over)= " << direct::elem(exp_Mweight, img_id, ihidden_over) << std::endl;
                                            REPORT_ERROR("weight is not a number");
                                            pthread_mutex_unlock(&global_mutex);
                                        }
                                        #endif
                                        // Keep track of sum and maximum of all weights for this particle
                                        // Later add all to exp_thisimage_sumweight, but inside this loop sum to local thisthread_sumweight first
                                        exp_thisimage_sumweight += weight;
                                    }
                                }
                            }
                            #ifdef TIMING
                            // Only time one thread, as I also only time one MPI process
                            if (part_id == mydata.sorted_idx[exp_my_first_part_id])
                                timer.toc(TIMING_WEIGHT_EXP);
                            #endif
                        }
                    }
                }
            }
        }

        //Store parameters for this image
        exp_sum_weight[img_id] = exp_thisimage_sumweight;

        // Check the sum of weights is not zero
        // On a Mac, the isnan function does not compile. Just uncomment the define statement, as this is merely a debugging statement
        // #define MAC_OSX
        #ifndef MAC_OSX
        if (exp_thisimage_sumweight == 0. || std::isnan(exp_thisimage_sumweight)) {
            std::cerr << " exp_thisimage_sumweight= " << exp_thisimage_sumweight << std::endl;
            Image<RFLOAT> It;
            It() = exp_Mweight;
            It.write("Mweight.spi");
            // It() = DEBUGGING_COPY_exp_Mweight;
            // It.write("Mweight_copy.spi");
            It().resize(exp_Mcoarse_significant);
            if (It().size() > 0) {
                for (long int n = 0; n < (It()).size(); n++) {
                    It()[n] = exp_Mcoarse_significant[n];
                }
                It.write("Mcoarse_significant.spi");
            }
            std::cerr << " part_id= " << part_id << std::endl;
            std::cerr << " img_id= " << img_id << std::endl;
            /*
            FourierTransformer transformer;
            MultidimArray<Complex> Faux = windowFourierTransform(exp_Fimg, mymodel.ori_size);
            It().resize(mymodel.ori_size, mymodel.ori_size);
            It() = transformer.inverseFourierTransform(Faux);
            CenterFFT(It(), -1);
            It.write("exp_Fimg.spi");
            std::cerr << "written exp_Fimgs.spi " << std::endl;
            */
            int group_id = mydata.getGroupId(part_id, img_id);
            std::cerr << " group_id= " << group_id << " mymodel.scale_correction[group_id]= " << mymodel.scale_correction[group_id] << std::endl;
            std::cerr << " exp_ipass= " << exp_ipass << std::endl;
            std::cerr << " sampling.NrDirections(0, true)= " << sampling.NrDirections()
                      << " sampling.NrDirections(0, false)= " << sampling.NrDirections(0, &exp_pointer_dir_nonzeroprior) << std::endl;
            std::cerr << " sampling.NrPsiSamplings(0, true)= " << sampling.NrPsiSamplings()
                      << " sampling.NrPsiSamplings(0, false)= " << sampling.NrPsiSamplings(0, &exp_pointer_psi_nonzeroprior) << std::endl;
            std::cerr << " mymodel.sigma2_noise[group_id]= " << mymodel.sigma2_noise[group_id] << std::endl;
            //std::cerr << " wsum_model.sigma2_noise[group_id]= " << wsum_model.sigma2_noise[group_id] << std::endl;
            if (do_norm_correction) {
                std::cerr << " mymodel.avg_norm_correction= " << mymodel.avg_norm_correction << std::endl;
                std::cerr << " wsum_model.avg_norm_correction= " << wsum_model.avg_norm_correction << std::endl;
            }

            std::cerr << "written out Mweight.spi" << std::endl;
            std::cerr << " exp_thisimage_sumweight= " << exp_thisimage_sumweight << std::endl;
            std::cerr << " exp_min_diff2[img_id]= " << exp_min_diff2[img_id] << std::endl;
            REPORT_ERROR("ERROR!!! zero sum of weights....");
        }
        #endif
    }

    // Initialise exp_Mcoarse_significant
    if (exp_ipass == 0)
        exp_Mcoarse_significant.resize(exp_nr_images, Xsize(exp_Mweight));

    // Now, for each image,  find the exp_significant_weight that encompasses adaptive_fraction of exp_sum_weight
    exp_significant_weight.clear();
    exp_significant_weight.resize(exp_nr_images, 0.0);

    for (int img_id = 0; img_id < exp_nr_images; img_id++) {

        int my_metadata_offset = metadata_offset + img_id;
        #ifdef TIMING
        if (part_id == mydata.sorted_idx[exp_my_first_part_id])
            timer.tic(TIMING_WEIGHT_SORT);
        #endif
        MultidimArray<RFLOAT> sorted_weight;
        // Get the relevant row for this particle
        exp_Mweight.getRow(img_id, sorted_weight);

        // Only select non-zero probabilities to speed up sorting
        long int np = 0;
        for (long int n = 0; n < sorted_weight.size(); n++) {
            if (sorted_weight[n] > 0.0) {
                sorted_weight[np++] = sorted_weight[n];
            }
        }
        sorted_weight.resize(np);

        // Sort from low to high values
        std::sort(sorted_weight.begin(), sorted_weight.end());

        #ifdef TIMING
        if (part_id == mydata.sorted_idx[exp_my_first_part_id])
            timer.toc(TIMING_WEIGHT_SORT);
        #endif
        RFLOAT frac_weight = 0.0;
        RFLOAT my_significant_weight;
        long int my_nr_significant_coarse_samples = 0;
        for (long int i = Xsize(sorted_weight) - 1; i >= 0; i--) {
            if (maximum_significants > 0 ) {
                if (my_nr_significant_coarse_samples < maximum_significants) {
                    if (exp_ipass == 0) { my_nr_significant_coarse_samples++; }
                    my_significant_weight = direct::elem(sorted_weight, i);
                }
            } else {
                if (exp_ipass == 0) { my_nr_significant_coarse_samples++; }
                my_significant_weight = direct::elem(sorted_weight, i);
            }
            frac_weight += direct::elem(sorted_weight, i);
            if (frac_weight > adaptive_fraction * exp_sum_weight[img_id])
                break;
        }

        #ifdef DEBUG_SORT
        // Check sorted array is really sorted
        RFLOAT prev = 0.0;
        for (long int n = 0; n < sorted_weight.size(); n++) {
            if (sorted_weight[n] < prev) {
                Image<RFLOAT>(sorted_weight * 10000).write("sorted_weight.spi");
                std::cerr << "written sorted_weight.spi" << std::endl;
                REPORT_ERROR("Error in sorting!");
            }
            prev = sorted_weight[n];
        }
        #endif

        if (exp_ipass == 0 && my_nr_significant_coarse_samples == 0) {
            std::cerr << " part_id= " << part_id << " img_id= " << img_id << " adaptive_fraction= " << adaptive_fraction << std::endl;
            std::cerr << " frac-weight= " << frac_weight << std::endl;
            std::cerr << " exp_sum_weight[img_id]= " << exp_sum_weight[img_id] << std::endl;
            std::cerr << " Xsize(exp_Mweight)= " << Xsize(exp_Mweight) << std::endl;
            Image<RFLOAT> It;
            It() = exp_Mweight;
            It() *= 10000;
            It.write("Mweight2.spi");
            std::cerr << "written Mweight2.spi" << std::endl;
            std::cerr << " np= " << np << std::endl;
            It() = sorted_weight;
            It() *= 10000;
            std::cerr << " Xsize(sorted_weight)= " << Xsize(sorted_weight) << std::endl;
            if (Xsize(sorted_weight) > 0) {
                It.write("sorted_weight.spi");
                std::cerr << "written sorted_weight.spi" << std::endl;
            }
            REPORT_ERROR("my_nr_significant_coarse_samples == 0");
        }

        if (exp_ipass == 0) {
            // Store nr_significant_coarse_samples for this particle
            // Don't do this for multibody, as it would be overwritten for each body,
            // and we also use METADATA_NR_SIGN in the new safeguard for the gold-standard separation
            if (mymodel.nr_bodies == 1)
                direct::elem(exp_metadata, my_metadata_offset, METADATA_NR_SIGN) = (RFLOAT) my_nr_significant_coarse_samples;

            // Keep track of which coarse samplings were significant were significant for this particle
            for (int ihidden = 0; ihidden < Xsize(exp_Mcoarse_significant); ihidden++) {
                direct::elem(exp_Mcoarse_significant, img_id, ihidden) = direct::elem(exp_Mweight, img_id, ihidden) >= my_significant_weight;
            }

        }
        exp_significant_weight[img_id] = my_significant_weight;
    }

    #ifdef TIMING
    if (part_id == mydata.sorted_idx[exp_my_first_part_id]) {
        timer.toc(exp_ipass == 0 ? TIMING_ESP_WEIGHT1 : TIMING_ESP_WEIGHT2);
    }
    #endif
}

template <typename T>
T update_and_remember(T &y, T z) {
    const T x = y;
    y = z;
    return x;
}

void MlOptimiser::storeWeightedSums(
    long int part_id, int ibody,
    int exp_current_oversampling, int metadata_offset,
    int exp_idir_min,   int exp_idir_max,   int exp_ipsi_min,   int exp_ipsi_max,
    int exp_itrans_min, int exp_itrans_max, int exp_iclass_min, int exp_iclass_max,
    std::vector<RFLOAT> &exp_min_diff2,
    std::vector<RFLOAT> &exp_highres_Xi2_img,
    std::vector<MultidimArray<Complex > > &exp_Fimg,
    std::vector<MultidimArray<Complex > > &exp_Fimg_nomask,
    std::vector<MultidimArray<RFLOAT> > &exp_Fctf,
    std::vector<MultidimArray<RFLOAT> > &exp_power_img,
    std::vector<Matrix1D<RFLOAT> > &exp_old_offset,
    std::vector<Matrix1D<RFLOAT> > &exp_prior,
    MultidimArray<RFLOAT> &exp_Mweight,
    MultidimArray<bool> &exp_Mcoarse_significant,
    std::vector<RFLOAT> &exp_significant_weight,
    std::vector<RFLOAT> &exp_sum_weight,
    std::vector<RFLOAT> &exp_max_weight,
    std::vector<int> &exp_pointer_dir_nonzeroprior, std::vector<int> &exp_pointer_psi_nonzeroprior,
    std::vector<RFLOAT> &exp_directions_prior, std::vector<RFLOAT> &exp_psi_prior,
    std::vector<std::vector<MultidimArray<Complex > > > &exp_local_Fimgs_shifted,
    std::vector<std::vector<MultidimArray<Complex > > > &exp_local_Fimgs_shifted_nomask,
    std::vector<MultidimArray<RFLOAT> > &exp_local_Minvsigma2,
    std::vector<MultidimArray<RFLOAT> > &exp_local_Fctf,
    std::vector<RFLOAT> &exp_local_sqrtXi2
) {
    #ifdef TIMING
    if (part_id == mydata.sorted_idx[exp_my_first_part_id])
        timer.tic(TIMING_ESP_WSUM);
    #endif

    int exp_nr_images = mydata.numberOfImagesInParticle(part_id);
    long int exp_nr_dir = (do_skip_align || do_skip_rotate) ? 1 : sampling.NrDirections(0, &exp_pointer_dir_nonzeroprior);
    long int exp_nr_psi = (do_skip_align || do_skip_rotate || do_only_sample_tilt) ? 1 : sampling.NrPsiSamplings(0, &exp_pointer_psi_nonzeroprior);
    long int exp_nr_trans = (do_skip_align) ? 1 : sampling.NrTranslationalSamplings();
    long int exp_nr_oversampled_rot = sampling.oversamplingFactorOrientations(exp_current_oversampling);
    long int exp_nr_oversampled_trans = sampling.oversamplingFactorTranslations(exp_current_oversampling);

    // Re-do below because now also want unmasked images AND if (stricht_highres_exp >0.0) then may need to resize
    precalculateShiftedImagesCtfsAndInvSigma2s(true, true, part_id, exp_current_oversampling, metadata_offset,
            exp_itrans_min, exp_itrans_max, exp_Fimg, exp_Fimg_nomask, exp_Fctf, exp_local_Fimgs_shifted, exp_local_Fimgs_shifted_nomask,
            exp_local_Fctf, exp_local_sqrtXi2, exp_local_Minvsigma2);

    // In doThreadPrecalculateShiftedImagesCtfsAndInvSigma2s() the origin of the exp_local_Minvsigma2s was omitted.
    // Set those back here
    for (int img_id = 0; img_id < exp_nr_images; img_id++) {
        int group_id = mydata.getGroupId(part_id, img_id);
        exp_local_Minvsigma2[img_id][0] = 1.0 / (sigma2_fudge * direct::elem(mymodel.sigma2_noise[group_id], 0));
    }

    // Initialise the maximum of all weights to a negative value
    exp_max_weight.clear();
    exp_max_weight.resize(exp_nr_images, -1.0);

    // For norm_correction and scale_correction of this particle
    std::vector<RFLOAT> exp_wsum_norm_correction;
    std::vector<RFLOAT> exp_wsum_scale_correction_XA, exp_wsum_scale_correction_AA;
    std::vector<RFLOAT> thr_wsum_signal_product_spectra, thr_wsum_reference_power_spectra;
    exp_wsum_norm_correction.resize(exp_nr_images, 0.0);

    // For scale_correction
    if (do_scale_correction) {
        exp_wsum_scale_correction_XA.resize(exp_nr_images);
        exp_wsum_scale_correction_AA.resize(exp_nr_images);
        thr_wsum_signal_product_spectra.resize(exp_nr_images);
        thr_wsum_reference_power_spectra.resize(exp_nr_images);
    }

    //Sigma2_noise estimation
    std::vector<MultidimArray<RFLOAT> > thr_wsum_sigma2_noise;
    // Wsum_sigma_noise2 is a 1D-spectrum for each img_id
    thr_wsum_sigma2_noise.resize(exp_nr_images);

    for (int img_id = 0; img_id < exp_nr_images; img_id++) {
        int optics_group = mydata.getOpticsGroup(part_id, img_id);
        thr_wsum_sigma2_noise[img_id].initZeros(image_full_size[optics_group]/2 + 1);
        if (do_scale_correction) {
            exp_wsum_scale_correction_XA    [img_id] = 0.0;
            exp_wsum_scale_correction_AA    [img_id] = 0.0;
            thr_wsum_signal_product_spectra [img_id] = 0.0;
            thr_wsum_reference_power_spectra[img_id] = 0.0;
        }
    }

    std::vector< RFLOAT> oversampled_rot, oversampled_tilt, oversampled_psi;
    std::vector<RFLOAT> oversampled_translations_x, oversampled_translations_y, oversampled_translations_z;
    Matrix2D<RFLOAT> A, Abody, Aori;
    MultidimArray<Complex > Fimg, Fref, Frefctf, Fimg_otfshift, Fimg_otfshift_nomask, Fimg_store_sgd;
    MultidimArray<RFLOAT> Minvsigma2, Mctf, Fweight;
    RFLOAT rot, tilt, psi;
    bool have_warned_small_scale = false;
    // Initialising... exp_Fimgs[0] has image_current_size[optics_group] (not coarse_size!)
    Fref.resize(exp_Fimg[0]);
    Frefctf.resize(exp_Fimg[0]);
    Fweight.resize(exp_Fimg[0]);
    Fimg.resize(exp_Fimg[0]);
    // Initialise Mctf to all-1 for if !do_ctf_corection
    Mctf.resize(exp_Fimg[0]);
    Mctf = 1.0;
    // Initialise Minvsigma2 to all-1 for if !do_map
    Minvsigma2.resize(exp_Fimg[0]);
    Minvsigma2 = 1.0;
    if (do_shifts_onthefly) {
        Fimg_otfshift.resize(Frefctf);
        Fimg_otfshift_nomask.resize(Frefctf);
    }
    if (do_sgd && ! do_avoid_sgd) {
        Fimg_store_sgd.resize(Frefctf);
    }


    if (mymodel.nr_bodies > 1) {
        RFLOAT rot_ori  = direct::elem(exp_metadata, metadata_offset, METADATA_ROT);
        RFLOAT tilt_ori = direct::elem(exp_metadata, metadata_offset, METADATA_TILT);
        RFLOAT psi_ori  = direct::elem(exp_metadata, metadata_offset, METADATA_PSI);
        Aori = Euler::angles2matrix(rot_ori, tilt_ori, psi_ori);
    }

    // Make local copies of weighted sums (except BPrefs, which are too big)
    // so that there are not too many mutex locks below
    std::vector<MultidimArray<RFLOAT> > thr_wsum_pdf_direction;
    std::vector<RFLOAT> thr_sumw_group, thr_wsum_pdf_class, thr_wsum_prior_offsetx_class, thr_wsum_prior_offsety_class;
    RFLOAT thr_wsum_sigma2_offset;
    // wsum_pdf_direction is a 1D-array (of length sampling.NrDirections()) for each class
    MultidimArray<RFLOAT> zeros = MultidimArray<RFLOAT>::zeros(sampling.NrDirections());
    thr_wsum_pdf_direction.resize(mymodel.nr_classes * mymodel.nr_bodies, zeros);
    // sumw_group is a RFLOAT for each group
    thr_sumw_group.resize(exp_nr_images, 0.0);
    // wsum_pdf_class is a RFLOAT for each class
    thr_wsum_pdf_class.resize(mymodel.nr_classes, 0.0);
    if (mymodel.ref_dim == 2) {
        thr_wsum_prior_offsetx_class.resize(mymodel.nr_classes, 0.0);
        thr_wsum_prior_offsety_class.resize(mymodel.nr_classes, 0.0);
    }
    // wsum_sigma2_offset is just a RFLOAT
    thr_wsum_sigma2_offset = 0.0;

    // Loop from iclass_min to iclass_max to deal with seed generation in first iteration
    for (int exp_iclass = exp_iclass_min;            exp_iclass <= exp_iclass_max; exp_iclass++)
    for (long int idir  = exp_idir_min, iorient = 0; idir <= exp_idir_max;         idir++)
    for (long int ipsi  = exp_ipsi_min;              ipsi <= exp_ipsi_max;         ipsi++, iorient++) {
        long int iorientclass = exp_iclass * exp_nr_dir * exp_nr_psi + iorient;

        // Only proceed if there was a significant coarsely sampled translation
        if (isSignificantAnyImageAnyTranslation(iorientclass, exp_itrans_min, exp_itrans_max, exp_Mcoarse_significant)) {

            // Now get the oversampled (rot, tilt, psi) triplets
            // This will be only the original (rot,tilt,psi) triplet if (adaptive_oversampling==0)
            sampling.getOrientations(
                idir, ipsi, adaptive_oversampling, oversampled_rot, oversampled_tilt, oversampled_psi,
                exp_pointer_dir_nonzeroprior, exp_directions_prior, exp_pointer_psi_nonzeroprior, exp_psi_prior
            );

            // The order of the looping here has changed for 3.1: different img_id have different optics_group and therefore different magnifications....
            for (int img_id = 0; img_id < exp_nr_images; img_id++) {
                int group_id = mydata.getGroupId(part_id, img_id);
                const int optics_group = mydata.getOpticsGroup(part_id, img_id);
                int my_metadata_offset = metadata_offset + img_id;
                RFLOAT my_pixel_size = mydata.getImagePixelSize(part_id, img_id);
                bool ctf_premultiplied = mydata.obsModel.getCtfPremultiplied(optics_group);

                // Loop over all oversampled orientations (only a single one in the first pass)
                for (long int iover_rot = 0; iover_rot < exp_nr_oversampled_rot; iover_rot++) {
                    rot  = oversampled_rot [iover_rot];
                    tilt = oversampled_tilt[iover_rot];
                    psi  = oversampled_psi [iover_rot];
                    // Get the Euler matrix
                    A = Euler::angles2matrix(rot, tilt, psi);


                    // For multi-body refinements, A are only 'residual' orientations, Abody is the complete Euler matrix
                    if (mymodel.nr_bodies > 1) {
                        Abody = Aori * (mymodel.orient_bodies[ibody]).transpose() * A_rot90 * A * mymodel.orient_bodies[ibody];
                        if (mydata.obsModel.hasMagMatrices) { Abody *= mydata.obsModel.anisoMag(optics_group); }
                        Abody *= mydata.obsModel.scaleDifference(optics_group, mymodel.ori_size, mymodel.pixel_size);
                    } else {
                        if (mydata.obsModel.hasMagMatrices) { A *= mydata.obsModel.anisoMag(optics_group); }
                        A *= mydata.obsModel.scaleDifference(optics_group, mymodel.ori_size, mymodel.pixel_size);
                    }

                    #ifdef TIMING
                    // Only time one thread, as I also only time one MPI process
                    if (part_id == mydata.sorted_idx[exp_my_first_part_id])
                        timer.tic(TIMING_WSUM_PROJ);
                    #endif
                    // Project the reference map (into Fref)
                    if (!do_skip_maximization) {
                        if (mymodel.nr_bodies > 1) {
                            Fref = mymodel.PPref[ibody].get2DFourierTransform(
                                exp_Fimg[0].xdim, exp_Fimg[0].ydim, exp_Fimg[0].zdim, Abody);
                        } else {
                            Fref = mymodel.PPref[exp_iclass].get2DFourierTransform(
                                exp_Fimg[0].xdim, exp_Fimg[0].ydim, exp_Fimg[0].zdim, A);
                        }
                    }
                    #ifdef TIMING
                    // Only time one thread, as I also only time one MPI process
                    if (part_id == mydata.sorted_idx[exp_my_first_part_id])
                        timer.toc(TIMING_WSUM_PROJ);
                    #endif
                    // Inside the loop over all translations and all part_id sum all shift Fimg's and their weights
                    // Then outside this loop do the actual backprojection
                    Fimg.initZeros();
                    Fweight.initZeros();

                    // This is an attempt to speed up illogically slow updates of wsum_sigma2_offset....
                    // It seems to make a big difference!
                    RFLOAT myprior_x, myprior_y, myprior_z, old_offset_z;
                    RFLOAT old_offset_x = exp_old_offset[img_id][0];
                    RFLOAT old_offset_y = exp_old_offset[img_id][1];
                    if (mymodel.ref_dim == 2 && mymodel.nr_bodies == 1) {
                        myprior_x = mymodel.prior_offset_class[exp_iclass][0];
                        myprior_y = mymodel.prior_offset_class[exp_iclass][1];
                    } else {
                        myprior_x = exp_prior[img_id][0];
                        myprior_y = exp_prior[img_id][1];
                        if (mymodel.data_dim == 3) {
                            myprior_z    = exp_prior[img_id][2];
                            old_offset_z = exp_old_offset[img_id][2];
                        }
                    }

                    if (!do_skip_maximization) {
                        if (do_map) { Minvsigma2 = exp_local_Minvsigma2[img_id]; }
                        // else Minvsigma2 was initialised to ones
                        // Apply CTF to reference projection
                        if (do_ctf_correction) {
                            Mctf = exp_local_Fctf[img_id];
                            if (refs_are_ctf_corrected) {
                                if (ctf_premultiplied) {
                                    for (long int n = 0; n < (Fref).size(); n++) {
                                        Frefctf[n] = Fref[n] * Mctf[n] * Mctf[n];
                                    }
                                } else {
                                    for (long int n = 0; n < (Fref).size(); n++) {
                                        Frefctf[n] = Fref[n] * Mctf[n];
                                    }
                                }
                            } else {
                                Frefctf = Fref;
                            }
                        } else {
                            // initialise because there are multiple particles and Mctf gets selfMultiplied for scale_correction
                            Mctf = 1.0;
                            Frefctf = Fref;
                        }
                        if (do_scale_correction) {
                            RFLOAT myscale = mymodel.scale_correction[group_id];
                            if (myscale > 10000.0) {
                                std::cerr << " rlnMicrographScaleCorrection= " << myscale << " group= " << group_id + 1 << std::endl;
                                REPORT_ERROR("ERROR: rlnMicrographScaleCorrection is very high. Did you normalize your data?");
                            }
                            else if (myscale < 0.001) {
                                if (!have_warned_small_scale) {
                                    std::cout << " WARNING: ignoring group " << group_id + 1
                                    << " with very small or negative scale (" << myscale << ");"
                                    << " Use larger groups for more stable scale estimates." << std::endl;
                                    have_warned_small_scale = true;
                                }
                                myscale = 0.001;
                            }
                            for (long int n = 0; n < (Frefctf).size(); n++) {
                                Frefctf[n] *= myscale;
                            }
                            // For CTF-terms in BP
                            Mctf *= myscale;
                        }
                    }

                    long int ihidden = iorientclass * exp_nr_trans;
                    for (long int itrans = exp_itrans_min, iitrans = 0; itrans <= exp_itrans_max; itrans++, ihidden++) {
                        // 1 Jun 2015 - Shaoda & Sjors, Helical refinement
                        sampling.getTranslationsInPixel(itrans, exp_current_oversampling, my_pixel_size, oversampled_translations_x, oversampled_translations_y, oversampled_translations_z,
                                (do_helical_refine) && (!ignore_helical_symmetry));
                        for (long int iover_trans = 0; iover_trans < exp_nr_oversampled_trans; iover_trans++, iitrans++) {
                            // Only deal with this sampling point if its weight was significant
                            long int ihidden_over = ihidden * exp_nr_oversampled_trans * exp_nr_oversampled_rot +
                                    iover_rot * exp_nr_oversampled_trans + iover_trans;
                            RFLOAT weight = direct::elem(exp_Mweight, img_id, ihidden_over);
                            // Only sum weights for non-zero weights
                            if (weight >= exp_significant_weight[img_id]) {
                                // Normalise the weight (do this after the comparison with exp_significant_weight!)
                                weight /= exp_sum_weight[img_id];

                                if (!do_skip_maximization) {

                                    #ifdef TIMING
                                    // Only time one thread, as I also only time one MPI process
                                    if (part_id == mydata.sorted_idx[exp_my_first_part_id])
                                        timer.tic(TIMING_WSUM_GETSHIFT);
                                    #endif

                                    /// Now get the shifted image
                                    // Use a pointer to avoid copying the entire array again in this highly expensive loop
                                    Complex *Fimg_shift, *Fimg_shift_nomask;
                                    if (!do_shifts_onthefly) {
                                        long int ishift = img_id * exp_nr_oversampled_trans * exp_nr_trans + iitrans;
                                        Fimg_shift        = exp_local_Fimgs_shifted       [img_id][ishift].data;
                                        Fimg_shift_nomask = exp_local_Fimgs_shifted_nomask[img_id][ishift].data;
                                    }
                                    else {
                                        // Feb01,2017 - Shaoda, on-the-fly shifts in helical reconstuctions (2D and 3D)
                                        if (do_helical_refine && !ignore_helical_symmetry) {
                                            RFLOAT xshift = 0.0, yshift = 0.0, zshift = 0.0;

                                            xshift = oversampled_translations_x[iover_trans];
                                            yshift = oversampled_translations_y[iover_trans];
                                            if (mymodel.data_dim == 3)
                                            zshift = oversampled_translations_z[iover_trans];

                                            RFLOAT rot_deg  = direct::elem(exp_metadata, my_metadata_offset, METADATA_ROT);
                                            RFLOAT tilt_deg = direct::elem(exp_metadata, my_metadata_offset, METADATA_TILT);
                                            RFLOAT psi_deg  = direct::elem(exp_metadata, my_metadata_offset, METADATA_PSI);
                                            transformCartesianAndHelicalCoords(
                                                xshift, yshift, zshift,
                                                rot_deg, tilt_deg, psi_deg,
                                                mymodel.data_dim,
                                                HELICAL_TO_CART_COORDS
                                            );
                                            // Fimg_shift
                                            shiftImageInFourierTransformWithTabSincos(
                                                exp_local_Fimgs_shifted[img_id][0],
                                                Fimg_otfshift,
                                                (RFLOAT) image_full_size[optics_group],
                                                image_current_size[optics_group],
                                                tab_sin, tab_cos,
                                                xshift, yshift, zshift
                                            );
                                            // Fimg_shift_nomask
                                            shiftImageInFourierTransformWithTabSincos(
                                                exp_local_Fimgs_shifted_nomask[img_id][0],
                                                Fimg_otfshift_nomask,
                                                (RFLOAT) image_full_size[optics_group],
                                                image_current_size[optics_group],
                                                tab_sin, tab_cos,
                                                xshift, yshift, zshift
                                            );
                                        } else {
                                            Complex *AB = (
                                                adaptive_oversampling == 0 ? global_fftshifts_ab_current :
                                                global_fftshifts_ab2_current
                                            )[optics_group][iitrans].data;
                                            for (long int n = 0; n < (exp_local_Fimgs_shifted[img_id][0]).size(); n++) {
                                                Complex A = AB[n];
                                                // Fimg_shift
                                                Complex X = exp_local_Fimgs_shifted[img_id][0][n];
                                                Fimg_otfshift[n] = Complex(
                                                    A.real * X.real - A.imag * X.imag, // A dot conj X
                                                    A.real * X.imag + A.imag * X.real  // A dot (i conj X)
                                                );
                                                // Fimg_shift_nomask
                                                Complex Y = exp_local_Fimgs_shifted_nomask[img_id][0][n];
                                                Fimg_otfshift_nomask[n] = Complex(
                                                    A.real * Y.real - A.imag * Y.imag,  // A dot conj Y
                                                    A.real * Y.imag + A.imag * Y.real   // A dot (i conj Y)
                                                );
                                            }
                                        }
                                        Fimg_shift = Fimg_otfshift.data;
                                        Fimg_shift_nomask = Fimg_otfshift_nomask.data;
                                    }
                                    #ifdef TIMING
                                    // Only time one thread, as I also only time one MPI process
                                    if (part_id == mydata.sorted_idx[exp_my_first_part_id]) {
                                        timer.toc(TIMING_WSUM_GETSHIFT);
                                        timer.tic(TIMING_WSUM_DIFF2);
                                    }
                                    #endif

                                    // Store weighted sum of squared differences for sigma2_noise estimation
                                    // Suggestion Robert Sinkovitz: merge difference and scale steps to make better use of cache
                                    for (long int n = 0; n < (Mresol_fine[optics_group]).size(); n++) {
                                        int ires = Mresol_fine[optics_group][n];
                                        if (ires > -1) {
                                            // Use FT of masked image for noise estimation!
                                            Complex A = Frefctf[n];
                                            Complex X = Fimg_shift[n];
                                            Complex diff = A - X;
                                            RFLOAT wdiff2 = weight * (diff.real * diff.real + diff.imag * diff.imag);
                                            // group-wise sigma2_noise
                                            thr_wsum_sigma2_noise[img_id][ires] += wdiff2;
                                            // For norm_correction
                                            exp_wsum_norm_correction[img_id] += wdiff2;
                                            if (do_scale_correction && direct::elem(mymodel.data_vs_prior_class[exp_iclass], ires) > 3.0) {
                                                RFLOAT AdotX = A.real * X.real + A.imag * X.imag;
                                                exp_wsum_scale_correction_XA[img_id] += weight * AdotX;
                                                RFLOAT AdotA = A.real * A.real + A.imag * A.imag;
                                                exp_wsum_scale_correction_AA[img_id] += weight * AdotA;
                                            }
                                        }
                                    }
                                    #ifdef TIMING
                                    // Only time one thread, as I also only time one MPI process
                                    if (part_id == mydata.sorted_idx[exp_my_first_part_id]) {
                                        timer.toc(TIMING_WSUM_DIFF2);
                                        timer.tic(TIMING_WSUM_LOCALSUMS);
                                    }
                                    #endif

                                    // Store sum of weights for this group
                                    thr_sumw_group[img_id] += weight;
                                    // Store weights for this class and orientation
                                    thr_wsum_pdf_class[exp_iclass] += weight;

                                    // The following goes MUCH faster than the original lines below....
                                    if (mymodel.ref_dim == 2) {
                                        thr_wsum_prior_offsetx_class[exp_iclass] += weight * my_pixel_size * (old_offset_x + oversampled_translations_x[iover_trans]);
                                        thr_wsum_prior_offsety_class[exp_iclass] += weight * my_pixel_size * (old_offset_y + oversampled_translations_y[iover_trans]);
                                    }
                                    // May18,2015 - Shaoda & Sjors, Helical refinement (translational searches)
                                    // Calculate the vector length of myprior
                                    RFLOAT mypriors_len2 = myprior_x * myprior_x + myprior_y * myprior_y;
                                    if (mymodel.data_dim == 3)
                                        mypriors_len2 += myprior_z * myprior_z;
                                    // If it is doing helical refinement AND Cartesian vector myprior has a length > 0, transform the vector to its helical coordinates
                                    if (do_helical_refine && !ignore_helical_symmetry && mypriors_len2 > 0.00001) {
                                        const RFLOAT rot_deg  = direct::elem(exp_metadata, my_metadata_offset, METADATA_ROT);
                                        const RFLOAT tilt_deg = direct::elem(exp_metadata, my_metadata_offset, METADATA_TILT);
                                        const RFLOAT psi_deg  = direct::elem(exp_metadata, my_metadata_offset, METADATA_PSI);
                                        transformCartesianAndHelicalCoords(
                                            myprior_x, myprior_y, myprior_z,
                                            rot_deg, tilt_deg, psi_deg,
                                            mymodel.data_dim, CART_TO_HELICAL_COORDS
                                        );
                                    }

                                    if (!do_helical_refine || ignore_helical_symmetry || mymodel.data_dim == 3) {
                                        RFLOAT diffx = myprior_x - old_offset_x - oversampled_translations_x[iover_trans];
                                        thr_wsum_sigma2_offset += weight * my_pixel_size * my_pixel_size * diffx * diffx;
                                    }
                                    RFLOAT diffy = myprior_y - old_offset_y - oversampled_translations_y[iover_trans];
                                    thr_wsum_sigma2_offset += weight * my_pixel_size * my_pixel_size * diffy * diffy;
                                    if (mymodel.data_dim == 3) {
                                        RFLOAT diffz = myprior_z - old_offset_z - oversampled_translations_z[iover_trans];
                                        if (!do_helical_refine || ignore_helical_symmetry)
                                            thr_wsum_sigma2_offset += weight * my_pixel_size * my_pixel_size * diffz * diffz;
                                    }

                                    // Store weight for this direction of this class
                                    if (do_skip_align || do_skip_rotate ) {
                                        // ignore pdf_direction
                                    } else if (mymodel.orientational_prior_mode == NOPRIOR) {
                                        thr_wsum_pdf_direction[exp_iclass][idir] += weight;
                                    } else {
                                        // In the case of orientational priors, get the original number of the direction back
                                        thr_wsum_pdf_direction[mymodel.nr_bodies > 1 ? ibody : exp_iclass][exp_pointer_dir_nonzeroprior[idir]] += weight;
                                    }

                                    #ifdef TIMING
                                    // Only time one thread, as I also only time one MPI process
                                    if (part_id == mydata.sorted_idx[exp_my_first_part_id]) {
                                        timer.toc(TIMING_WSUM_LOCALSUMS);
                                        timer.tic(TIMING_WSUM_SUMSHIFT);
                                    }
                                    #endif

                                    Complex *Fimg_store = [&] () {
                                        if (!do_sgd || do_avoid_sgd)
                                            return Fimg_shift_nomask;
                                        for (long int n = 0; n < Frefctf.size(); n++) {
                                            Fimg_store_sgd[n] = Fimg_shift_nomask[n] - Frefctf[n];
                                        }
                                        return Fimg_store_sgd.data;
                                    }();
                                    // #define DEBUG_BODIES2
                                    #ifdef DEBUG_BODIES2
                                    FourierTransformer transformer;
                                    MultidimArray<Complex> Ftt(Frefctf);
                                    for (long int n = 0; n < Ftt.size(); n++)
                                        Ftt[n] = Fimg_store[n];

                                    Image<RFLOAT> tt;
                                    tt().resize(exp_current_image_size, exp_current_image_size);
                                    tt() = transformer.inverseFourierTransform(Ftt);
                                    CenterFFT(tt(), -1);
                                    FileName fnt = "BPimg_body" + integerToString(ibody + 1, 1) + "_ihidden" + integerToString(ihidden_over) + ".spi";
                                    tt.write(fnt);
                                    Ftt = Frefctf;
                                    tt().resize(exp_current_image_size, exp_current_image_size);
                                    tt() = transformer.inverseFourierTransform(Ftt);
                                    CenterFFT(tt(), -1);
                                    fnt = "Fref_body" + integerToString(ibody + 1, 1) + "_ihidden" + integerToString(ihidden_over) + ".spi";
                                    tt.write(fnt);


                                    std::cerr << " rot= " << rot << " tilt= " << tilt << " psi= " << psi << std::endl;
                                    std::cerr << " itrans= " << itrans << " iover_trans= " << iover_trans << std::endl;
                                    std::cerr << " ihidden_over= " << ihidden_over << " weight= " << weight << std::endl;
                                    std::cerr << "written " << fnt <<std::endl;
                                    #endif

                                    // Store sum of weight*SSNR*Fimg in data and sum of weight*SSNR in weight
                                    // Use the FT of the unmasked image to back-project in order to prevent reconstruction artefacts! SS 25oct11
                                    if (ctf_premultiplied) {
                                        for (long int n = 0; n < Fimg.size(); n++) {
                                            RFLOAT ctf = Mctf[n];
                                            RFLOAT weightxinvsigma2 = weight * Minvsigma2[n];
                                            // now Fimg stores sum of all shifted w*Fimg
                                            Fimg[n] += Fimg_store[n] * weightxinvsigma2;
                                            // now Fweight stores sum of all w and multiply by CTF^2
                                            Fweight[n] += weightxinvsigma2 * ctf * ctf;
                                        }
                                    } else {
                                        for (long int n = 0; n < Fimg.size(); n++) {
                                            RFLOAT ctf = Mctf[n];
                                            RFLOAT weightxinvsigma2 = weight * ctf * Minvsigma2[n];
                                            // now Fimg stores sum of all shifted w*Fimg
                                            Fimg[n] += Fimg_store[n] * weightxinvsigma2;
                                            // now Fweight stores sum of all w
                                            // Note that CTF needs to be squared in Fweight, weightxinvsigma2 already contained one copy
                                            Fweight[n] += weightxinvsigma2 * ctf;

                                        }
                                    }

                                #ifdef TIMING
                                    // Only time one thread, as I also only time one MPI process
                                    if (part_id == mydata.sorted_idx[exp_my_first_part_id])
                                        timer.toc(TIMING_WSUM_SUMSHIFT);
                                #endif
                                }

                                // Keep track of max_weight and the corresponding optimal hidden variables
                                if (weight > exp_max_weight[img_id]) {
                                    // Store optimal image parameters
                                    exp_max_weight[img_id] = weight;

                                    // This is not necessary as rot, tilt and psi remain unchanged!
                                    // angles_t angles = Euler::matrix2angles(A);

                                    const int imx = mymodel.nr_bodies == 1 ? 0 : METADATA_LINE_LENGTH_BEFORE_BODIES + ibody * METADATA_NR_BODY_PARAMS;
                                    int icol_rot  = METADATA_ROT  + imx;
                                    int icol_tilt = METADATA_TILT + imx;
                                    int icol_psi  = METADATA_PSI  + imx;
                                    int icol_xoff = METADATA_XOFF + imx;
                                    int icol_yoff = METADATA_YOFF + imx;
                                    int icol_zoff = METADATA_ZOFF + imx;

                                    RFLOAT old_rot  = update_and_remember(direct::elem(exp_metadata, my_metadata_offset, icol_rot),  rot);
                                    RFLOAT old_tilt = update_and_remember(direct::elem(exp_metadata, my_metadata_offset, icol_tilt), tilt);
                                    RFLOAT old_psi  = update_and_remember(direct::elem(exp_metadata, my_metadata_offset, icol_psi),  psi);

                                    Matrix1D<RFLOAT> shifts(mymodel.data_dim);
                                    // include old_offsets for normal refinement (i.e. non multi-body)
                                    shifts[0] = exp_old_offset[img_id][0] + oversampled_translations_x[iover_trans];
                                    shifts[1] = exp_old_offset[img_id][1] + oversampled_translations_y[iover_trans];
                                    if (mymodel.data_dim == 3)
                                    shifts[2] = exp_old_offset[img_id][2] + oversampled_translations_z[iover_trans];

                                    #ifdef DEBUG_BODIES2
                                    std::cerr << ihidden_over << " weight= " << weight;
                                    std::cerr << " exp_old_offset= " << exp_old_offset[img_id].transpose() << std::endl;
                                    std::cerr << " SET: rot= " << rot << " tilt= " << tilt << " psi= " << psi;
                                    std::cerr << " xx-old= " << exp_old_offset[img_id][0];
                                    std::cerr << " yy-old= " << exp_old_offset[img_id][1];
                                    std::cerr << " add-xx= " << oversampled_translations_x[iover_trans];
                                    std::cerr << " add-yy= " << oversampled_translations_y[iover_trans];
                                    std::cerr << " xnew= " << shifts[0];
                                    std::cerr << " ynew= " << shifts[1];
                                    std::cerr << std::endl;
                                    #endif

                                    #ifdef DEBUG_HELICAL_ORIENTATIONAL_SEARCH
                                    std::cerr << "MlOptimiser::storeWeightedSums()" << std::endl;
                                    if (mymodel.data_dim == 2) {
                                        std::cerr << " exp_old_offset = (" << exp_old_offset[img_id][0] << ", " << exp_old_offset[img_id][1] << ")" << std::endl;
                                        std::cerr << " Oversampled trans = (" << oversampled_translations_x[iover_trans] << ", " << oversampled_translations_y[iover_trans] << ")" << std::endl;
                                        std::cerr << " shifts = (" << shifts[0] << ", " << shifts[1] << ")" << std::endl;
                                    } else {
                                        std::cerr << " exp_old_offset = (" << exp_old_offset[img_id][0] << ", " << exp_old_offset[img_id][1] << ", " << exp_old_offset[img_id][2] << ")" << std::endl;
                                        std::cerr << " Oversampled trans = (" << oversampled_translations_x[iover_trans] << ", " << oversampled_translations_y[iover_trans] << ", " << oversampled_translations_z[iover_trans] << ")" << std::endl;
                                        std::cerr << " shifts = (" << shifts[0] << ", " << shifts[1] << ", " << shifts[2] << ")" << std::endl;
                                    }
                                    #endif

                                    // Helical reconstruction: use oldpsi-angle to rotate back the exp_old_offset[0] + oversampled_translations_x[iover_trans] and
                                    if (do_helical_refine && !ignore_helical_symmetry) {
                                        // Bring xshift, yshift and zshift back to cartesian coords for outputting in the STAR file
                                        #ifdef DEBUG_HELICAL_ORIENTATIONAL_SEARCH
                                        std::cerr << "MlOptimiser::storeWeightedSums()" << std::endl;
                                        std::cerr << "Bring xy(z) shifts back to Cartesian coordinates for output in the STAR file" << std::endl;
                                        std::cerr << " itrans = " << itrans << ", iover_trans = " << iover_trans << std::endl;
                                        std::cerr << "  old_psi = " << old_psi << " degrees";
                                        if (shifts.size() == 2) {
                                            std::cerr << std::endl;
                                            std::cerr << "  Helical offsets (r, p) = (" << shifts[0] << ", " << shifts[1] << ")" << std::endl;
                                        } else {
                                            std::cerr << ", old_tilt = " << old_tilt << " degrees" << std::endl;
                                            std::cerr << "  Helical offsets (p1, p2, r) = (" << shifts[0] << ", " << shifts[1] << ", " << shifts[2] << ")" << std::endl;
                                        }
                                        #endif
                                        transformCartesianAndHelicalCoords(shifts, shifts, old_rot, old_tilt, old_psi, HELICAL_TO_CART_COORDS);
                                        #ifdef DEBUG_HELICAL_ORIENTATIONAL_SEARCH
                                        std::cerr << "  Cartesian offsets " << (
                                            shifts.size() == 2 ?
                                                "(x, y) = (" << shifts[0] << ", " << shifts[1] << ")" :
                                            "(x, y, z) = (" << shifts[0] << ", " << shifts[1] << ", " << shifts[2] << ")"
                                        ) << std::endl;
                                        #endif
                                    }

                                    direct::elem(exp_metadata, my_metadata_offset, icol_xoff) = shifts[0];
                                    direct::elem(exp_metadata, my_metadata_offset, icol_yoff) = shifts[1];
                                    if (mymodel.data_dim == 3)
                                    direct::elem(exp_metadata, my_metadata_offset, icol_zoff) = shifts[2];

                                    if (ibody == 0) {
                                        direct::elem(exp_metadata, my_metadata_offset, METADATA_CLASS) = (RFLOAT) exp_iclass + 1;
                                        direct::elem(exp_metadata, my_metadata_offset, METADATA_PMAX) = exp_max_weight[img_id];
                                    }
                                }
                            }
                        }
                    }
                    #ifdef RELION_TESTING
                    std::string fnm = std::string("cpu_out_exp_wsum_norm_correction.txt");
                    char *text = &fnm[0];
                    freopen(text, "w", stdout);
                    printf("%4.8f \n", exp_wsum_norm_correction);
                    fclose(stdout);
                    //----------
                    fnm = std::string("cpu_out_thr_wsum_sigma2_noise.txt");
                    text = &fnm[0];
                    freopen(text, "w", stdout);
                    for (long int n = 0; n < (Mresol_fine[optics_group]).size(); n++) {
                        printf("%4.8f \n", thr_wsum_sigma2_noise[0].data[n]);
                    }
                    fclose(stdout);
                    //----------
                    fnm = std::string("cpu_out_Fweights.txt");
                    text = &fnm[0];
                    freopen(text, "w", stdout);
                    for (int n = 0; n < 1000; n++) {
                        printf("%4.8f \n", *(Fweight.data + n * 60 + 50));
                    }
                    fclose(stdout);
                    #endif

                    if (!do_skip_maximization) {
                        #ifdef TIMING
                        // Only time one thread, as I also only time one MPI process
                        if (part_id == mydata.sorted_idx[exp_my_first_part_id])
                            timer.tic(TIMING_WSUM_BACKPROJ);
                        #endif
                        // Perform the actual back-projection.
                        // This is done with the sum of all (in-plane) shifted Fimg's
                        // Perform this inside a mutex
                        int my_mutex = exp_iclass % NR_CLASS_MUTEXES;
                        pthread_mutex_lock(&global_mutex2[my_mutex]);
                        if (mymodel.nr_bodies > 1) {
                            wsum_model.BPref[ibody].set2DFourierTransform(Fimg, Abody, &Fweight);
                        } else {
                            wsum_model.BPref[exp_iclass].set2DFourierTransform(Fimg, A, &Fweight);
                        }
                        pthread_mutex_unlock(&global_mutex2[my_mutex]);
                        #ifdef TIMING
                        // Only time one thread, as I also only time one MPI process
                        if (part_id == mydata.sorted_idx[exp_my_first_part_id])
                            timer.toc(TIMING_WSUM_BACKPROJ);
                        #endif
                    }
                }
            }
        }
    }

    // Extend norm_correction and sigma2_noise estimation to higher resolutions for all particles
    // Also calculate dLL for each particle and store in metadata
    RFLOAT thr_avg_norm_correction = 0.0;
    RFLOAT thr_sum_dLL = 0., thr_sum_Pmax = 0.0;

    // loop over all images inside this particle
    for (int img_id = 0; img_id < exp_nr_images; img_id++) {
        int group_id = mydata.getGroupId(part_id, img_id);
        int my_metadata_offset = metadata_offset + img_id;
        int optics_group = mydata.getOpticsGroup(part_id, img_id);
        RFLOAT my_pixel_size = mydata.getOpticsPixelSize(optics_group);
        int my_image_size = mydata.getOpticsImageSize(optics_group);

        // If the current images were smaller than the original size, fill the rest of wsum_model.sigma2_noise with the power_class spectrum of the images
        for (int ires = image_current_size[optics_group]/2 + 1; ires < image_full_size[optics_group] / 2 + 1; ires++) {
            direct::elem(thr_wsum_sigma2_noise[img_id], ires) += direct::elem(exp_power_img[img_id], ires);
            // Also extend the weighted sum of the norm_correction
            exp_wsum_norm_correction[img_id] += direct::elem(exp_power_img[img_id], ires);
        }

        // Store norm_correction
        // Multiply by old value because the old norm_correction term was already applied to the image
        // Don't do this for multi-body refinement, where one always uses the norm_correction from the consensus refinement
        if (do_norm_correction && mymodel.nr_bodies == 1) {
            RFLOAT old_norm_correction = direct::elem(exp_metadata, my_metadata_offset, METADATA_NORM);
            old_norm_correction /= mymodel.avg_norm_correction;
            // The factor two below is because exp_wsum_norm_correctiom is similar to sigma2_noise, which is the variance for the real/imag components
            // The variance of the total image (on which one normalizes) is twice this value!
            RFLOAT normcorr = old_norm_correction * sqrt(exp_wsum_norm_correction[img_id] * 2.0);
            thr_avg_norm_correction += normcorr;
            // Now set the new norm_correction in the relevant position of exp_metadata
            direct::elem(exp_metadata, my_metadata_offset, METADATA_NORM) = normcorr;

            // Print warning for strange norm-correction values
            if (!do_cc() && direct::elem(exp_metadata, my_metadata_offset, METADATA_NORM) > 10.0) {
                std::cout << " WARNING: norm_correction= "<< direct::elem(exp_metadata, my_metadata_offset, METADATA_NORM)
                        << " for particle " << part_id << " in group " << group_id + 1
                        << "; Are your groups large enough?  Or is the reference on the correct greyscale?" << std::endl;
            }

        }

        // Store weighted sums for scale_correction
        if (do_scale_correction) {
            // Divide XA by the old scale_correction and AA by the square of that, because was incorporated into Fctf
            exp_wsum_scale_correction_XA[img_id] /= mymodel.scale_correction[group_id];
            exp_wsum_scale_correction_AA [img_id]/= mymodel.scale_correction[group_id] * mymodel.scale_correction[group_id];

            thr_wsum_signal_product_spectra[img_id] += exp_wsum_scale_correction_XA[img_id];
            thr_wsum_reference_power_spectra[img_id] += exp_wsum_scale_correction_AA[img_id];
        }

        // Calculate DLL for each particle
        RFLOAT logsigma2 = 0.0;
        RFLOAT remap_image_sizes = (mymodel.ori_size * mymodel.pixel_size) / (my_image_size * my_pixel_size);
        for (long int n = 0; n < (Mresol_fine[optics_group]).size(); n++) {
            int ires = Mresol_fine[optics_group][n];
            int ires_remapped = round(remap_image_sizes * ires);
            // Note there is no sqrt in the normalisation term because of the 2-dimensionality of the complex-plane
            // Also exclude origin from logsigma2, as this will not be considered in the P-calculations
            if (ires > 0 && ires_remapped < Xsize(mymodel.sigma2_noise[group_id]))
                logsigma2 += log(2.0 * PI * direct::elem(mymodel.sigma2_noise[group_id], ires_remapped));
        }
        if (exp_sum_weight[img_id] == 0) {
            std::cerr << " part_id= " << part_id << std::endl;
            std::cerr << " img_id= " << img_id << std::endl;
            std::cerr << " exp_min_diff2[img_id]= " << exp_min_diff2[img_id]<< std::endl;
            std::cerr << " logsigma2= " << logsigma2 << std::endl;
            int group_id = mydata.getGroupId(part_id, img_id);
            std::cerr << " group_id= " << group_id << std::endl;
            std::cerr << " ml_model.scale_correction[group_id]= " << mymodel.scale_correction[group_id] << std::endl;
            std::cerr << " exp_significant_weight[img_id]= " << exp_significant_weight[img_id] << std::endl;
            std::cerr << " exp_max_weight[img_id]= " << exp_max_weight[img_id] << std::endl;
            std::cerr << " ml_model.sigma2_noise[group_id]= " << mymodel.sigma2_noise[group_id] << std::endl;
            REPORT_ERROR("ERROR: exp_sum_weight[img_id]==0");
        }
        RFLOAT dLL = do_cc() ?
            -exp_min_diff2[img_id] :
            log(exp_sum_weight[img_id]) - exp_min_diff2[img_id] - logsigma2;

        // Store dLL of each image in the output array, and keep track of total sum
        direct::elem(exp_metadata, my_metadata_offset, METADATA_DLL) = dLL;
        thr_sum_dLL += dLL;

        // Also store sum of Pmax
        thr_sum_Pmax += direct::elem(exp_metadata, my_metadata_offset, METADATA_PMAX);

    }

    // Now, inside a global_mutex, update the other weighted sums among all threads
    if (!do_skip_maximization) {
        pthread_mutex_lock(&global_mutex);
        for (int img_id = 0; img_id < exp_nr_images; img_id++) {
            long int igroup = mydata.getGroupId(part_id, img_id);
            int optics_group = mydata.getOpticsGroup(part_id, img_id);
            int my_image_size = mydata.getOpticsImageSize(optics_group);
            RFLOAT my_pixel_size = mydata.getOpticsPixelSize(optics_group);
            RFLOAT remap_image_sizes = (mymodel.ori_size * mymodel.pixel_size) / (my_image_size * my_pixel_size);
            for (long int i = 0; i < Xsize(thr_wsum_sigma2_noise[img_id]); i++) {
                int i_resam = round(i * remap_image_sizes);
                if (i_resam < Xsize(wsum_model.sigma2_noise[igroup])) {
                    direct::elem(wsum_model.sigma2_noise[igroup], i_resam) += direct::elem(thr_wsum_sigma2_noise[img_id], i);
                }
            }
            wsum_model.sumw_group[igroup] += thr_sumw_group[img_id];
            if (do_scale_correction) {
                wsum_model.wsum_signal_product [igroup] += thr_wsum_signal_product_spectra [img_id];
                wsum_model.wsum_reference_power[igroup] += thr_wsum_reference_power_spectra[img_id];
            }
        }
        for (int n = 0; n < mymodel.nr_classes; n++) {
            wsum_model.pdf_class[n] += thr_wsum_pdf_class[n];
            if (mymodel.ref_dim == 2) {
                wsum_model.prior_offset_class[n][0] += thr_wsum_prior_offsetx_class[n];
                wsum_model.prior_offset_class[n][1] += thr_wsum_prior_offsety_class[n];
            }
            #ifdef CHECKSIZES
            if (Xsize(wsum_model.pdf_direction[n]) != Xsize(thr_wsum_pdf_direction[n])) {
                std::cerr << " Xsize(wsum_model.pdf_direction[n])= " << Xsize(wsum_model.pdf_direction[n]) << " Xsize(thr_wsum_pdf_direction[n])= " << Xsize(thr_wsum_pdf_direction[n]) << std::endl;
                REPORT_ERROR("Xsize(wsum_model.pdf_direction[n]) != Xsize(thr_wsum_pdf_direction[n])");
            }
            #endif
        }
        for (int n = 0; n < mymodel.nr_classes * mymodel.nr_bodies; n++) {
            if (!do_skip_align && !do_skip_rotate)
                wsum_model.pdf_direction[n] += thr_wsum_pdf_direction[n];
        }
        wsum_model.sigma2_offset += thr_wsum_sigma2_offset;
        if (do_norm_correction && mymodel.nr_bodies == 1)
            wsum_model.avg_norm_correction += thr_avg_norm_correction;
        wsum_model.LL += thr_sum_dLL;
        wsum_model.ave_Pmax += thr_sum_Pmax;
        pthread_mutex_unlock(&global_mutex);
    }

    #ifdef TIMING
    if (part_id == mydata.sorted_idx[exp_my_first_part_id])
        timer.toc(TIMING_ESP_WSUM);
    #endif
}

/** Monitor the changes in the optimal translations, orientations and class assignments for some particles */
void MlOptimiser::monitorHiddenVariableChanges(long int my_first_part_id, long int my_last_part_id) {

    for (long int part_id_sorted = my_first_part_id, metadata_offset = 0; part_id_sorted <= my_last_part_id; part_id_sorted++) {

        const long int part_id = mydata.sorted_idx[part_id_sorted];
        for (int img_id = 0; img_id < mydata.numberOfImagesInParticle(part_id); img_id++, metadata_offset++) {

            const long int ori_img_id = mydata.particles[part_id].images[img_id].id;
            RFLOAT my_pixel_size = mydata.getImagePixelSize(part_id, img_id);

            for (int ibody = 0; ibody < mymodel.nr_bodies; ibody++) {

                if (mymodel.nr_bodies > 1 && mymodel.keep_fixed_bodies[ibody] > 0)
                    continue;

                RFLOAT old_rot, old_tilt, old_psi, old_xoff, old_yoff, old_zoff = 0.0;
                RFLOAT rot, tilt, psi, xoff, yoff, zoff = 0.0;
                int old_iclass, iclass;

                if (mymodel.nr_bodies > 1) {

                    // Old optimal parameters
                    old_rot  = mydata.MDbodies[ibody].getValue<RFLOAT>(EMDL::ORIENT_ROT,               ori_img_id);
                    old_tilt = mydata.MDbodies[ibody].getValue<RFLOAT>(EMDL::ORIENT_TILT,              ori_img_id);
                    old_psi  = mydata.MDbodies[ibody].getValue<RFLOAT>(EMDL::ORIENT_PSI,               ori_img_id);
                    old_xoff = mydata.MDbodies[ibody].getValue<RFLOAT>(EMDL::ORIENT_ORIGIN_X_ANGSTROM, ori_img_id);
                    old_yoff = mydata.MDbodies[ibody].getValue<RFLOAT>(EMDL::ORIENT_ORIGIN_Y_ANGSTROM, ori_img_id);
                    if (mymodel.data_dim == 3) {
                    old_zoff = mydata.MDbodies[ibody].getValue<RFLOAT>(EMDL::ORIENT_ORIGIN_Z_ANGSTROM, ori_img_id);
                    }
                    old_iclass = 0;

                    // New optimal parameters
                    const int imx = METADATA_LINE_LENGTH_BEFORE_BODIES + ibody * METADATA_NR_BODY_PARAMS;
                    rot  = direct::elem(exp_metadata, metadata_offset, imx + METADATA_ROT);
                    tilt = direct::elem(exp_metadata, metadata_offset, imx + METADATA_TILT);
                    psi  = direct::elem(exp_metadata, metadata_offset, imx + METADATA_PSI);
                    xoff = direct::elem(exp_metadata, metadata_offset, imx + METADATA_XOFF) * my_pixel_size;
                    yoff = direct::elem(exp_metadata, metadata_offset, imx + METADATA_YOFF) * my_pixel_size;
                    if (mymodel.data_dim == 3)
                    zoff = direct::elem(exp_metadata, metadata_offset, imx + METADATA_ZOFF) * my_pixel_size;
                    iclass = 0;

                } else {

                    // Old optimal parameters
                    old_rot  = mydata.MDimg.getValue<RFLOAT>(EMDL::ORIENT_ROT,               ori_img_id);
                    old_tilt = mydata.MDimg.getValue<RFLOAT>(EMDL::ORIENT_TILT,              ori_img_id);
                    old_psi  = mydata.MDimg.getValue<RFLOAT>(EMDL::ORIENT_PSI,               ori_img_id);
                    old_xoff = mydata.MDimg.getValue<RFLOAT>(EMDL::ORIENT_ORIGIN_X_ANGSTROM, ori_img_id);
                    old_yoff = mydata.MDimg.getValue<RFLOAT>(EMDL::ORIENT_ORIGIN_Y_ANGSTROM, ori_img_id);
                    if (mymodel.data_dim == 3)
                    old_zoff = mydata.MDimg.getValue<RFLOAT>(EMDL::ORIENT_ORIGIN_Z_ANGSTROM, ori_img_id);
                    old_iclass = mydata.MDimg.getValue<int>(EMDL::PARTICLE_CLASS, ori_img_id);

                    // New optimal parameters
                    rot  = direct::elem(exp_metadata, metadata_offset, METADATA_ROT);
                    tilt = direct::elem(exp_metadata, metadata_offset, METADATA_TILT);
                    psi  = direct::elem(exp_metadata, metadata_offset, METADATA_PSI);
                    xoff = direct::elem(exp_metadata, metadata_offset, METADATA_XOFF) * my_pixel_size;
                    yoff = direct::elem(exp_metadata, metadata_offset, METADATA_YOFF) * my_pixel_size;
                    if (mymodel.data_dim == 3)
                    zoff = direct::elem(exp_metadata, metadata_offset, METADATA_ZOFF) * my_pixel_size;
                    iclass = (int) direct::elem(exp_metadata, metadata_offset, METADATA_CLASS);

                }

                // Some orientational distance...
                sum_changes_optimal_orientations += sampling.calculateAngularDistance(rot, tilt, psi, old_rot, old_tilt, old_psi);
                sum_changes_optimal_offsets += (xoff - old_xoff) * (xoff - old_xoff) + (yoff - old_yoff) * (yoff - old_yoff) + (zoff - old_zoff) * (zoff - old_zoff);
                if (iclass != old_iclass)
                    sum_changes_optimal_classes += 1.0;
                sum_changes_count += 1.0;

            }

        }

    }

}

void MlOptimiser::updateOverallChangesInHiddenVariables() {

    // Calculate hidden variable changes
    if (sum_changes_count > 0.0) {
        current_changes_optimal_classes = sum_changes_optimal_classes / sum_changes_count;
        current_changes_optimal_orientations = sum_changes_optimal_orientations / sum_changes_count;
        current_changes_optimal_offsets = sqrt(sum_changes_optimal_offsets / (2. * sum_changes_count));
    } else {
        current_changes_optimal_classes = 0.0;
        current_changes_optimal_orientations = 0.0;
        current_changes_optimal_offsets = 0.0;
    }

    // Reset the sums
    sum_changes_optimal_classes = 0.0;
    sum_changes_optimal_orientations = 0.0;
    sum_changes_optimal_offsets = 0.0;
    sum_changes_count = 0.0;

    RFLOAT ratio_orient_changes = current_changes_optimal_orientations /  sampling.getAngularSampling(adaptive_oversampling);
    RFLOAT ratio_trans_changes = current_changes_optimal_offsets /  sampling.getTranslationalSampling(adaptive_oversampling);

    // Update nr_iter_wo_large_hidden_variable_changes if all three assignment types are within 3% of the smallest thus far
    // Or if changes in offsets or orientations are smaller than 40% of the current sampling
    if (1.03 * current_changes_optimal_classes >= smallest_changes_optimal_classes &&
        (ratio_trans_changes < 0.40 || 1.03 * current_changes_optimal_offsets >= smallest_changes_optimal_offsets) &&
        (ratio_orient_changes < 0.40 || 1.03 * current_changes_optimal_orientations >= smallest_changes_optimal_orientations) )
        nr_iter_wo_large_hidden_variable_changes++;
    else
        nr_iter_wo_large_hidden_variable_changes = 0;

    // Update smallest changes in hidden variables thus far
    if (current_changes_optimal_classes < smallest_changes_optimal_classes)
        smallest_changes_optimal_classes = round(current_changes_optimal_classes);
    if (current_changes_optimal_offsets < smallest_changes_optimal_offsets)
        smallest_changes_optimal_offsets = current_changes_optimal_offsets;
    if (current_changes_optimal_orientations < smallest_changes_optimal_orientations)
        smallest_changes_optimal_orientations = current_changes_optimal_orientations;

}


void MlOptimiser::calculateExpectedAngularErrors(long int my_first_part_id, long int my_last_part_id) {

    long int n_trials = 0;
    for (int part_id = my_first_part_id; part_id <= my_last_part_id; part_id++) {
        n_trials += mydata.numberOfImagesInParticle(part_id);

        if (mydata.numberOfImagesInParticle(part_id) > 1) {
            REPORT_ERROR("ERROR: calculateExpectedAngularErrors will not work for multiple images per particle from 3.1...");
        }
    }

    // Separate angular error estimate for each of the classes
    acc_rot = acc_trans = 999.0; // later minimum will be taken to find the best class...

    // P(X | X_1) / P(X | X_2) = exp ( |F_1 - F_2|^2 / (-2 sigma2) )
    // exp(-4.60517) = 0.01
    RFLOAT pvalue = 4.60517;
    //if (mymodel.data_dim == 3)
    //	pvalue *= 2.0;

    std::cout << " Estimating accuracies in the orientational assignment ... " << std::endl;
    int nr_particles = my_last_part_id - my_first_part_id + 1;
    init_progress_bar(nr_particles * mymodel.nr_classes * mymodel.nr_bodies);
    for (int iclass = 0; iclass < mymodel.nr_classes * mymodel.nr_bodies; iclass++) {

        // Don't do this for (almost) empty classes, but always for multi-body refinement
        if (mymodel.nr_bodies == 1 && mymodel.pdf_class[iclass] < 0.01) {
            mymodel.acc_rot  [iclass] = 999.0;
            mymodel.acc_trans[iclass] = 999.0;
            continue;
        }

        // Initialise the orientability arrays that will be written out in the model.star file
        // These are for the user's information only: nothing will be actually done with them
        #ifdef DEBUG_CHECKSIZES
        if (iclass >= (mymodel.orientability_contrib).size()) {
            std::cerr << "iclass= " << iclass << " (mymodel.orientability_contrib).size()= " << (mymodel.orientability_contrib).size() << std::endl;
            REPORT_ERROR("iclass >= (mymodel.orientability_contrib).size()");
        }
        #endif
        mymodel.orientability_contrib[iclass].initZeros(mymodel.ori_size / 2 + 1);

        RFLOAT acc_rot_class   = 0.0;
        RFLOAT acc_trans_class = 0.0;
        // Particles are already in random order, so just move from 0 to n_trials
        for (long int part_id_sorted = my_first_part_id, metadata_offset = 0; part_id_sorted <= my_last_part_id; part_id_sorted++) {

            const long int part_id = mydata.sorted_idx[part_id_sorted];
            for (int img_id = 0; img_id < mydata.numberOfImagesInParticle(part_id); img_id++, metadata_offset++) {

                int group_id = mydata.getGroupId(part_id, img_id);
                RFLOAT my_pixel_size = mydata.getImagePixelSize(part_id, img_id);
                const int optics_group = mydata.getOpticsGroup(part_id, img_id);
                bool ctf_premultiplied = mydata.obsModel.getCtfPremultiplied(optics_group);

                // Set current_image_size to the coarse_size to calculate expected angular errors
                int current_image_size = strict_highres_exp > 0.0 && !do_acc_currentsize_despite_highres_exp ?
                    image_coarse_size [optics_group] :  // Use smaller images in both passes and keep a maximum on coarse_size, just like in FREALIGN
                    image_current_size[optics_group];   // Use smaller images in the first pass, but larger ones in the second pass


                MultidimArray<RFLOAT> Fctf;
                // Get CTF for this particle
                if (do_ctf_correction) {
                    if (mymodel.data_dim == 3) {
                        // Read CTF-image from disc
                        FileName fn_ctf;
                        try {
                            fn_ctf = mydata.getImageNameOnScratch(part_id, img_id, true);
                        } catch (const char *errmsg) {
                            std::istringstream split(exp_fn_ctf);
                            // Get the right line in the exp_fn_img string
                            for (int i = 0; i <= metadata_offset; i++)
                                getline(split, fn_ctf);
                        }
                        auto Ictf = Image<RFLOAT>::from_filename(fn_ctf);
                        Fctf.resize(current_image_size / 2 + 1, current_image_size, current_image_size);

                        // If there is a redundant half, get rid of it
                        if (Xsize(Ictf()) == Ysize(Ictf())) {
                            // Set the CTF-image in Fctf
                            Ictf().setXmippOrigin();
                            FOR_ALL_ELEMENTS_IN_FFTW_TRANSFORM(Fctf) {
                                // Use negative kp, ip and jp indices, because the origin in the ctf_img lies half a pixel to the right of the actual center....
                                direct::elem(Fctf, i, j, k) = Ictf().elem(-ip, -jp, -kp);
                            }
                        } else if (Xsize(Ictf()) == Ysize(Ictf()) / 2 + 1) {
                            // otherwise, just window the CTF to the current resolution
                            Fctf = windowFourierTransform(Ictf(), Ysize(Fctf));
                        } else {
                            // if dimensions are neither cubical nor FFTW, stop
                            REPORT_ERROR("3D CTF volume must be either cubical or adhere to FFTW format!");
                        }
                    } else {

                        // Get parameters that change per-particle from the exp_metadata
                        CTF ctf = CtfHelper::makeCTF(
                            &mydata.obsModel, optics_group,
                            direct::elem(exp_metadata, metadata_offset, METADATA_CTF_DEFOCUS_U),
                            direct::elem(exp_metadata, metadata_offset, METADATA_CTF_DEFOCUS_V),
                            direct::elem(exp_metadata, metadata_offset, METADATA_CTF_DEFOCUS_ANGLE),
                            direct::elem(exp_metadata, metadata_offset, METADATA_CTF_BFACTOR),
                            direct::elem(exp_metadata, metadata_offset, METADATA_CTF_KFACTOR),
                            direct::elem(exp_metadata, metadata_offset, METADATA_CTF_PHASE_SHIFT)
                        );

                        Fctf = CtfHelper::getFftwImage(
                            ctf, current_image_size / 2 + 1, current_image_size,
                            image_full_size[optics_group], image_full_size[optics_group], my_pixel_size,
                            &mydata.obsModel,
                            ctf_phase_flipped, only_flip_phases, intact_ctf_first_peak, true, do_ctf_padding
                        );
                    }
                }

                // Search 2 times: ang and off
                // Don't estimate rotational accuracies if we're doing do_skip_rotate
                int imode_start = do_skip_rotate;
                for (int imode = imode_start; imode < 2; imode++) {
                    RFLOAT ang_error = 0.0;
                    RFLOAT sh_error = 0.0;
                    RFLOAT my_snr = 0.0;

                    // Search for ang_error and sh_error where there are at least 3-sigma differences!
                    // 13feb12: change for explicit probability at P=0.01
                    while (my_snr <= pvalue) {
                        // Graduallly increase the step size
                        const RFLOAT ang_step =
                            ang_error <  0.2 ? 0.05 :
                            ang_error <  1.0 ? 0.10 :
                            ang_error <  2.0 ? 0.20 :
                            ang_error <  5.0 ? 0.50 :
                            ang_error < 10.0 ? 1.00 :
                            ang_error < 20.0 ? 2.00 :
                                               5.00 ;
                        ang_error += ang_step;

                        const RFLOAT sh_step =
                            sh_error <  1.0 ? 0.1 :
                            sh_error <  2.0 ? 0.2 :
                            sh_error <  5.0 ? 0.5 :
                            sh_error < 10.0 ? 1.0 :
                                              2.0 ;
                        sh_error += sh_step;

                        // Prevent an endless while by putting boundaries on ang_error and sh_error
                        if (imode == 0 && ang_error > 30.0 || imode == 1 && sh_error > 10.0)
                            break;

                        // ori_img_id to keep exactly the same as in relion-3.0....
                        init_random_generator(random_seed + mydata.getOriginalImageId(part_id, img_id));

                        RFLOAT rot1  = direct::elem(exp_metadata, metadata_offset, METADATA_ROT);
                        RFLOAT tilt1 = direct::elem(exp_metadata, metadata_offset, METADATA_TILT);
                        RFLOAT psi1  = direct::elem(exp_metadata, metadata_offset, METADATA_PSI);
                        RFLOAT xoff1 = 0.0;
                        RFLOAT yoff1 = 0.0;
                        RFLOAT zoff1 = 0.0;

                        // Get the FT of the first image
                        Matrix2D<RFLOAT> A1 = Euler::angles2matrix(rot1, tilt1, psi1);
                        if (mydata.obsModel.hasMagMatrices) { A1 *= mydata.obsModel.anisoMag(optics_group); }
                        A1 *= mydata.obsModel.scaleDifference(optics_group, mymodel.ori_size, mymodel.pixel_size);
                        auto F1 = mymodel.PPref[iclass].get2DFourierTransform(
                            current_image_size / 2 + 1, current_image_size, mymodel.data_dim == 2 ? 1 : current_image_size, A1);

                        // Apply the angular or shift error
                        RFLOAT rot2  = rot1;
                        RFLOAT tilt2 = tilt1;
                        RFLOAT psi2  = psi1;
                        RFLOAT xshift = xoff1;
                        RFLOAT yshift = yoff1;
                        RFLOAT zshift = zoff1;

                        // Perturb psi or xoff , depending on the mode
                        if (imode == 0) {
                            if (mymodel.ref_dim == 3) {
                                // Randomly change rot, tilt or psi
                                RFLOAT ran = rnd_unif();
                                if (ran < 0.3333) {
                                    rot2 = rot1 + ang_error;
                                } else if (ran < 0.6667) {
                                    tilt2 = tilt1 + ang_error;
                                } else {
                                    psi2  = psi1 + ang_error;
                                }
                            } else {
                                psi2 = psi1 + ang_error;
                            }
                        } else {
                            // Randomly change xoff or yoff
                            RFLOAT ran = rnd_unif ();
                            if (mymodel.data_dim == 3) {
                                if (ran < 0.3333) {
                                    xshift = xoff1 + sh_error;
                                } else if (ran < 0.6667) {
                                    yshift = yoff1 + sh_error;
                                } else {
                                    zshift = zoff1 + sh_error;
                                }
                            } else {
                                if (ran < 0.5) {
                                    xshift = xoff1 + sh_error;
                                } else {
                                    yshift = yoff1 + sh_error;
                                }
                            }
                        }

                        // Get the FT of the second image
                        auto F2 = MultidimArray<Complex>::zeros(
                            current_image_size / 2 + 1,
                            current_image_size,
                            mymodel.data_dim == 3 ? current_image_size : 1
                        );

                        if (imode == 0) {
                            // Get new rotated version of reference
                            Matrix2D<RFLOAT> A2 = Euler::angles2matrix(rot2, tilt2, psi2);
                            if (mydata.obsModel.hasMagMatrices) { A2 *= mydata.obsModel.anisoMag(optics_group); }
                            A2 *= mydata.obsModel.scaleDifference(optics_group, mymodel.ori_size, mymodel.pixel_size);
                            F2 = mymodel.PPref[iclass].get2DFourierTransform(
                                current_image_size / 2 + 1, current_image_size,
                                mymodel.data_dim == 3 ? current_image_size : 1, A2);
                        } else {
                            // Get shifted version
                            shiftImageInFourierTransform(
                                F1, F2, (RFLOAT) image_full_size[optics_group],
                                -xshift, -yshift, -zshift
                            );
                        }

                        // Apply CTF to F1 and F2 if necessary
                        if (do_ctf_correction) {
                            #ifdef DEBUG_CHECKSIZES
                            if (!Fctf.sameShape(F1) || !Fctf.sameShape(F2)) {
                                std::cerr<<" Fctf: "; Fctf.printShape(std::cerr);
                                std::cerr<<" F1:   "; F1.printShape(std::cerr);
                                std::cerr<<" F2:   "; F2.printShape(std::cerr);
                                REPORT_ERROR("ERROR: Fctf has a different shape from F1 and F2");
                            }
                            #endif
                            for (long int n = 0; n < F1.size(); n++) {
                                F1[n] *= Fctf[n];
                                F2[n] *= Fctf[n];
                            }
                            if (ctf_premultiplied) {
                                for (long int n = 0; n < F1.size(); n++) {
                                    F1[n] *= Fctf[n];
                                    F2[n] *= Fctf[n];
                                }
                            }
                        }

                        RFLOAT remap_image_sizes =
                            (mymodel.ori_size * mymodel.pixel_size) /
                            (image_full_size[optics_group] * my_pixel_size);
                        const auto &myMresol = Ysize(F1) == image_coarse_size[optics_group] ?
                            Mresol_coarse[optics_group] : Mresol_fine[optics_group];
                        my_snr = 0.0;
                        for (long int n = 0; n < F1.size(); n++) {
                            int ires = myMresol[n];
                            int ires_remapped = round(remap_image_sizes * ires);
                            if (ires > 0 && ires_remapped < Xsize(mymodel.sigma2_noise[group_id])) {
                                my_snr += norm(F1[n] - F2[n]) / (2 * sigma2_fudge * mymodel.sigma2_noise[group_id].elem(ires_remapped));
                            }
                        }

                        // Only for the psi-angle and the translations, and only when my_prob < 0.01 calculate a histogram of the contributions at each resolution shell
                        if (my_snr > pvalue && imode == 0) {
                            for (long int n = 0; n < F1.size(); n++) {
                                int ires = myMresol[n];
                                int ires_remapped = round(remap_image_sizes * ires);
                                if (ires > 0 && ires_remapped < Xsize(mymodel.sigma2_noise[group_id]))
                                    mymodel.orientability_contrib[iclass].elem(ires_remapped) +=
                                        norm(F1[n] - F2[n]) /
                                        (2 * sigma2_fudge * mymodel.sigma2_noise[group_id].elem(ires_remapped));
                            }
                        }

                    }
                    if (imode == 0) {
                        acc_rot_class += ang_error;
                    } else if (imode == 1) {
                        acc_trans_class += my_pixel_size * sh_error; // now in Angstroms!
                    }
                }
            }

            progress_bar(n_trials*iclass + metadata_offset);
        }

        mymodel.acc_rot  [iclass] = acc_rot_class   / (RFLOAT) n_trials;
        mymodel.acc_trans[iclass] = acc_trans_class / (RFLOAT) n_trials;

        // Store normalised spectral contributions to orientability
        if (mymodel.orientability_contrib[iclass].sum() > 0.0)
            mymodel.orientability_contrib[iclass]   /= mymodel.orientability_contrib[iclass].sum();

        // Keep the orientational accuracy of the best class for the auto-sampling approach
        acc_rot   = std::min(mymodel.acc_rot  [iclass], acc_rot);
        acc_trans = std::min(mymodel.acc_trans[iclass], acc_trans);


        // Richard's formula with Greg's constant
        //RFLOAT b_orient = (acc_rot_class*acc_rot_class* particle_diameter*particle_diameter) / 3000.0;
        //std::cout << " + expected B-factor from the orientational errors = "
        //		<< b_orient<<std::endl;
        // B=8 PI^2 U^2
        //std::cout << " + expected B-factor from the translational errors = "
        //		<< 8 * PI * PI * mymodel.pixel_size * mymodel.pixel_size * acc_trans_class * acc_trans_class << std::endl;

    } // end loop iclass
    progress_bar(n_trials * mymodel.nr_classes * mymodel.nr_bodies);

    std::cout << " Auto-refine: Estimated accuracy angles= " << acc_rot << " degrees; offsets= " << acc_trans << " Angstroms" << std::endl;
    // Warn for inflated resolution estimates
    if (acc_rot > 10.0 && do_auto_refine) {
        std::cout << " Auto-refine: WARNING: The angular accuracy is worse than 10 degrees, so basically you cannot align your particles (yet)!" << std::endl;
        std::cout << " Auto-refine: WARNING: You probably need not worry if the accuracy improves during the next few iterations." << std::endl;
        std::cout << " Auto-refine: WARNING: However, if the problem persists it may lead to spurious FSC curves, so be wary of inflated resolution estimates..." << std::endl;
        std::cout << " Auto-refine: WARNING: Sometimes it is better to tune resolution yourself by adjusting T in a 3D-classification with a single class." << std::endl;
    }

}

void MlOptimiser::updateAngularSampling(bool myverb) {

    // For SGD: onyl update initial angular sampling after the sgd_ini_iter have passed
    if (do_sgd) {
        if (iter > sgd_ini_iter) {
            RFLOAT min_sampling = 360.0 / ceil(PI * particle_diameter * mymodel.current_resolution);
            RFLOAT old_rottilt_step = sampling.getAngularSampling(adaptive_oversampling);
            if (old_rottilt_step > min_sampling) {
                has_fine_enough_angular_sampling = false;

                int new_hp_order = sampling.healpix_order + 1;
                RFLOAT new_rottilt_step = 60.0 / exp2(new_hp_order + adaptive_oversampling);

                // Set the new sampling in the sampling-object
                sampling.setOrientations(new_hp_order, new_rottilt_step * exp2(adaptive_oversampling));

                // Resize the pdf_direction arrays to the correct size and fill with an even distribution
                mymodel.initialisePdfDirection(sampling.NrDirections());

                // Also reset the nr_directions in wsum_model
                wsum_model.nr_directions = mymodel.nr_directions;

                // Also resize and initialise wsum_model.pdf_direction for each class!
                for (int iclass = 0; iclass < mymodel.nr_classes * mymodel.nr_bodies; iclass++)
                    wsum_model.pdf_direction[iclass].initZeros(mymodel.nr_directions);
            }
        }
    } else if (mymodel.nr_classes > 1 && allow_coarser_samplings) {
        // For 2D/3D classification: use coarser angular and translational samplings when estimated accuracies are low...

        // A. Coarser rotational sampling
        // Stay a bit on the safe side: 80% of estimated accuracy
        if (mymodel.ref_dim == 3) {

            // If doing CC first iteration, there will not be a acc_rot yet: use minimum sampling based on resolution instead
            RFLOAT my_min_sampling = iter == 1 && do_firstiter_cc ?
                360.0 / ceil(PI * particle_diameter * mymodel.current_resolution) :
                acc_rot;

            // 3D classification
            int previous_healpix_order = sampling.healpix_order;
            // Always go down from original healpix order
            sampling.healpix_order = sampling.healpix_order_ori;
            bool is_decreased = false;
            while (sampling.getAngularSampling(adaptive_oversampling) < 0.8 * my_min_sampling) {
                sampling.healpix_order--;
                is_decreased = true;
            }
            // Now have healpix_order that gives coarser sampling than min_sampling, go one up to have finer sampling than min_sampling
            // Only do this is (is_decreased), as we don't want to go finer than the original one...
            if (is_decreased) { sampling.healpix_order++; }

            // Don't go beyond original sampling
            sampling.healpix_order = std::min(sampling.healpix_order, sampling.healpix_order_ori);
            if (sampling.healpix_order != previous_healpix_order) {
                sampling.setOrientations(sampling.healpix_order, sampling.getAngularSampling());

                // Resize the pdf_direction arrays to the correct size and fill with an even distribution
                mymodel.initialisePdfDirection(sampling.NrDirections());

                // Also reset the nr_directions in wsum_model
                wsum_model.nr_directions = mymodel.nr_directions;

                // Also resize and initialise wsum_model.pdf_direction for each class!
                for (int iclass = 0; iclass < mymodel.nr_classes * mymodel.nr_bodies; iclass++)
                    wsum_model.pdf_direction[iclass].initZeros(mymodel.nr_directions);

            }
        } else {
            // 2D classification
            RFLOAT new_psi_step = 0.8 * acc_rot * exp2(adaptive_oversampling);
            new_psi_step = std::max(new_psi_step, sampling.psi_step_ori);
            if (fabs(new_psi_step - sampling.psi_step) > 0.001 ) {
                sampling.setOrientations(-1, new_psi_step);
            }
        }

        // B. Coarser translational sampling
        // Stay a bit on the safe side: 80% of estimated accuracy
        RFLOAT new_offset_step = 0.8 * acc_trans * exp2(adaptive_oversampling);
        // Don't go coarser than the 95% of the offset_range (so at least 5 samplings are done)
        new_offset_step = std::min(new_offset_step, 0.95 * sampling.offset_range);
        // Don't go finer than the original offset_step!
        new_offset_step = std::max(new_offset_step, sampling.offset_step_ori);
        sampling.setTranslations(
            new_offset_step, sampling.offset_range, false,
            do_helical_refine && !ignore_helical_symmetry,
            sampling.helical_offset_step, helical_rise_initial, helical_twist_initial
        );

        // Print to screen
        if (myverb) {
            std::cout << " Coarser-sampling: Angular step= " << sampling.getAngularSampling(adaptive_oversampling) << " degrees." << std::endl;
            std::cout << " Coarser-sampling: Offset search range= " << sampling.offset_range << " Angstroms; offset step= " << sampling.getTranslationalSampling(adaptive_oversampling) << " Angstroms" << std::endl;
        }

    } else {

        if (!do_split_random_halves)
            REPORT_ERROR("MlOptimiser::updateAngularSampling: BUG! updating of angular sampling should only happen for gold-standard (auto-) refinements.");

        if (do_skip_rotate)
            REPORT_ERROR("ERROR: --skip_rotate can only be used in movie-frame refinement ...");

        // Only change the sampling if the resolution has not improved during the last 2 iterations
        // AND the hidden variables have not changed during the last 2 iterations
        RFLOAT old_rottilt_step = sampling.getAngularSampling(adaptive_oversampling);

        // If the angular accuracy and the necessary angular step for the current resolution is finer than the current angular step, make it finer.
        // But don't go to local search until it stabilises or look at change in angles?
        int nr_ang_steps = ceil(PI * particle_diameter * mymodel.current_resolution);
        RFLOAT myresol_angstep = 360.0 / nr_ang_steps;
        // But don't go down to local searches too early, i.e. at last exhaustive sampling first stabilise resolution
        bool do_proceed_resolution = auto_resolution_based_angles && myresol_angstep < old_rottilt_step && sampling.healpix_order + 1 != autosampling_hporder_local_searches ||
            nr_iter_wo_resol_gain >= MAX_NR_ITER_WO_RESOL_GAIN;

        const bool do_proceed_hidden_variables = auto_ignore_angle_changes || nr_iter_wo_large_hidden_variable_changes >= MAX_NR_ITER_WO_LARGE_HIDDEN_VARIABLE_CHANGES;

        // Only use a finer angular sampling if the angular accuracy is still above 75% of the estimated accuracy
        // If it is already below, nothing will change and eventually nr_iter_wo_resol_gain or nr_iter_wo_large_hidden_variable_changes will go above MAX_NR_ITER_WO_RESOL_GAIN
        if (do_proceed_resolution && do_proceed_hidden_variables ) {

            bool all_bodies_are_done = false;
            // For multi-body refinement: switch off those bodies that don't have high enough angular accuracy
            if (mymodel.nr_bodies > 1) {
                all_bodies_are_done = true;
                for (int ibody = 0; ibody < mymodel.nr_bodies; ibody++) {
                    // Stop multi-body refinements a bit earlier than normal ones: no 75%, but 90% of accuracy
                    // if has_converged: in the final iteration include all bodies again!
                    if (old_rottilt_step < 0.90 * mymodel.acc_rot[ibody] && !has_converged) {
                        if (myverb)
                            std::cout << " Body: " <<ibody << " with rotational accuracy of " << mymodel.acc_rot[ibody] << " will be kept fixed " << std::endl;
                        mymodel.keep_fixed_bodies[ibody] = 1;
                    } else {
                        all_bodies_are_done = false;
                    }
                }
            }

            if (
                all_bodies_are_done ||
                old_rottilt_step < maximum_angular_sampling && maximum_angular_sampling > 0.0 ||
                old_rottilt_step < 0.75 * acc_rot && (old_rottilt_step <= minimum_angular_sampling || minimum_angular_sampling <= 0.0)
                // Old rottilt step is already below 75% of estimated accuracy: have to stop refinement?
                // If a minimum_angular_sampling is given and we're not there yet, also just continue
            ) {
                // don't change angular sampling, as it is already fine enough
                has_fine_enough_angular_sampling = true;
            } else {
                has_fine_enough_angular_sampling = false;

                // A. Use translational sampling as suggested by acc_trans

                // Prevent very coarse translational samplings: max 1.5
                // Also stay a bit on the safe side with the translational sampling: 75% of estimated accuracy
                RFLOAT new_step = std::min(1.5, 0.75 * acc_trans) * exp2(adaptive_oversampling);

                // For subtomogram averaging: use at least half times previous step size
                if (mymodel.data_dim == 3) // TODO: check: this might just as well work for 2D data...
                    new_step = std::max(sampling.offset_step / 2.0, new_step);

                // Search ranges are five times the last observed changes in offsets
                // Only 3x for subtomogram averaging....
                RFLOAT new_range = mymodel.data_dim == 2 ? 5.0 * current_changes_optimal_offsets : 3 * current_changes_optimal_offsets;

                // New range can only become 30% bigger than the previous range (to prevent very slow iterations in the beginning)
                new_range = std::min(1.3 * sampling.offset_range, new_range);

                // Prevent too narrow searches: always at least 3x3 pixels in the coarse search
                new_range = std::max(new_range, 1.5 * new_step);

                // Also prevent too wide searches: that will lead to memory problems:
                // If steps size < 1/4th of search range, then decrease search range by 50%
                if (new_range > 4.0 * new_step) { new_range /= 2.0; }

                // If even that was not enough: use coarser step size and hope things will settle down later...
                if (new_range > 4.0 * new_step) { new_step = new_range / 4.0; }

                // 8 Jun 2015 Shaoda & Sjors, Helical refinement
                RFLOAT new_helical_offset_step = sampling.helical_offset_step;
                if (mymodel.ref_dim == 3) {
                    // 8 Jun 2015 Shaoda & Sjors, Helical refinement
                    // new_helical_offset_step /= 2.0;

                    // AFTER AUG17,2015
                    // ??? new_step ~= 1/4 * helical_offset_step ??? still divide helical_offset_step by 2 ! That is reasonable... Because it cannot happen (see above)
                    if (new_step < new_helical_offset_step) { new_helical_offset_step /= 2.0; }
                }

                // B. Use twice as fine angular sampling
                int new_hp_order;
                RFLOAT new_rottilt_step, new_psi_step;
                if (mymodel.ref_dim == 3) {
                    new_hp_order = sampling.healpix_order + 1;
                    new_rottilt_step = new_psi_step = 60.0 / exp2(new_hp_order + adaptive_oversampling);

                    // Set the new sampling in the sampling-object
                    sampling.setOrientations(new_hp_order, new_psi_step * exp2(adaptive_oversampling));

                    // Resize the pdf_direction arrays to the correct size and fill with an even distribution
                    mymodel.initialisePdfDirection(sampling.NrDirections());

                    // Also reset the nr_directions in wsum_model
                    wsum_model.nr_directions = mymodel.nr_directions;

                    // Also resize and initialise wsum_model.pdf_direction for each class!
                    for (int iclass = 0; iclass < mymodel.nr_classes * mymodel.nr_bodies; iclass++)
                        wsum_model.pdf_direction[iclass].initZeros(mymodel.nr_directions);

                } else if (mymodel.ref_dim == 2) {
                    sampling.psi_step /= 2.0;
                } else {
                    REPORT_ERROR("MlOptimiser::autoAdjustAngularSampling BUG: ref_dim should be two or three");
                }

                // 8 Jun 2015 Shaoda & Sjors, Helical refinement
                bool do_local_searches_helical = do_auto_refine && do_helical_refine &&
                    sampling.healpix_order >= autosampling_hporder_local_searches;

                // Don't go to coarse angular samplings. Then just keep doing as it was
                if (new_step > sampling.offset_step) {
                    new_step = sampling.offset_step;
                    new_range = sampling.offset_range;
                }

                sampling.setTranslations(new_step, new_range, do_local_searches_helical, (do_helical_refine) && (!ignore_helical_symmetry), new_helical_offset_step, helical_rise_initial, helical_twist_initial);

                // Reset iteration counters
                nr_iter_wo_resol_gain = 0;
                nr_iter_wo_large_hidden_variable_changes = 0;

                // Reset smallest changes hidden variables
                smallest_changes_optimal_classes = 9999999;
                smallest_changes_optimal_offsets = 999.0;
                smallest_changes_optimal_orientations = 999.0;

                // If the angular sampling is smaller than autosampling_hporder_local_searches, then use local searches of +/- 6 times the angular sampling
                if (mymodel.ref_dim == 3 && new_hp_order >= autosampling_hporder_local_searches) {
                    // Switch ON local angular searches
                    mymodel.orientational_prior_mode = PRIOR_ROTTILT_PSI;
                    mymodel.sigma2_rot = mymodel.sigma2_psi = 4.0 * new_rottilt_step * new_rottilt_step;
                    if (!do_helical_refine || !helical_keep_tilt_prior_fixed)
                        mymodel.sigma2_tilt = mymodel.sigma2_rot;

                    // 20 Aug 2015 - Shaoda, Helical refinement
                    if (do_helical_refine && !ignore_helical_symmetry)
                        mymodel.sigma2_rot = getHelicalSigma2Rot(
                            helical_rise_initial, helical_twist_initial,
                            sampling.helical_offset_step, new_rottilt_step, mymodel.sigma2_rot
                        );
                }
            }
        }

        // Print to screen
        if (myverb) {
            std::cout << " Auto-refine: Angular step= " << sampling.getAngularSampling(adaptive_oversampling) << " degrees; local searches= "
            << (mymodel.orientational_prior_mode == NOPRIOR ? "false" : "true") << std::endl;
            // 8 Jun 2015 Shaoda & Sjors, Helical refine
            if (do_helical_refine && !ignore_helical_symmetry) {
                std::cout << " Auto-refine: Helical refinement... Local translational searches along helical axis= "
                << (mymodel.ref_dim == 3 && do_auto_refine && sampling.healpix_order >= autosampling_hporder_local_searches ? "true" : "false")
                << std::endl;
            }
            std::cout << " Auto-refine: Offset search range= " << sampling.offset_range << " Angstroms; offset step= " << sampling.getTranslationalSampling(adaptive_oversampling) << " Angstroms";
            if (do_helical_refine && !ignore_helical_symmetry)
                std::cout << "; offset step along helical axis= " << sampling.getHelicalTranslationalSampling(adaptive_oversampling) << " pixels";
            std::cout << std::endl;
        }
    }
}

void MlOptimiser::updateSubsetSize(bool myverb) {
    // If we're doing cisTEM-like acceleration of refinement through subsets: set the subset size here
    long int old_subset_size = subset_size;
    if (do_fast_subsets) {
        long int min_parts_per_class = mymodel.ref_dim == 2 ? 100 : 1500;
        if (iter <= 5) {
            subset_size = min_parts_per_class * mymodel.nr_classes;
        } else if (iter <= 10) {
            subset_size = 3 * min_parts_per_class * mymodel.nr_classes;
        } else if (iter <= 15) {
            subset_size = std::max((double) (3 * min_parts_per_class * mymodel.nr_classes), 0.3 * mydata.numberOfParticles());
        } else {
            subset_size = -1;
        }
        if (subset_size > mydata.numberOfParticles())
            subset_size = -1;
    } else if (do_sgd) {
        // Do sgd_ini_iter iterations with completely identical K references, sigd_ini_subset_size, enforce non-negativity and sgd_ini_resol resolution limit
        if (iter < sgd_ini_iter) {
            subset_size = sgd_ini_subset_size;
        } else if (iter < sgd_ini_iter + sgd_inbetween_iter) {
            subset_size = sgd_ini_subset_size + round((RFLOAT(iter - sgd_ini_iter) / RFLOAT(sgd_inbetween_iter)) * sgd_fin_subset_size-sgd_ini_subset_size);
        } else {
            subset_size = sgd_fin_subset_size;
        }
        if (subset_size > mydata.numberOfParticles())
            subset_size = -1;
    }

    if (myverb && subset_size != old_subset_size)
        std::cout << " Setting subset size to " << subset_size << " particles" << std::endl;
}

void MlOptimiser::checkConvergence(bool myverb) {

    if (
        has_fine_enough_angular_sampling &&
        nr_iter_wo_resol_gain >= MAX_NR_ITER_WO_RESOL_GAIN &&
        (auto_ignore_angle_changes || nr_iter_wo_large_hidden_variable_changes >= MAX_NR_ITER_WO_LARGE_HIDDEN_VARIABLE_CHANGES)
    ) {
        has_converged = true;
        do_join_random_halves = true;
        // In the last iteration, include all data until Nyquist
        do_use_all_data = true;

        // For multibody refinement: reset all bodies to not-fixed (if they were originally) for the final iteration with all data to Nyquist
        if (mymodel.nr_bodies > 1) {
            for (int ibody = 0; ibody < mymodel.nr_bodies; ibody++) {
                mymodel.keep_fixed_bodies[ibody] =
                    mymodel.sigma_tilt_bodies  [ibody] < 0.001 &&
                    mymodel.sigma_psi_bodies   [ibody] < 0.001 &&
                    mymodel.sigma_offset_bodies[ibody] < 0.001;
            }
        }
    }

    if (myverb) {
        std::cout << " Auto-refine: Iteration= "<< iter<< std::endl;
        std::cout << " Auto-refine: Resolution= "<< 1.0 / mymodel.current_resolution<< " (no gain for " << nr_iter_wo_resol_gain << " iter) "<< std::endl;
        std::cout << " Auto-refine: Changes in angles= " << current_changes_optimal_orientations << " degrees; and in offsets= " << current_changes_optimal_offsets
        << " Angstroms (no gain for " << nr_iter_wo_large_hidden_variable_changes << " iter) "<< std::endl;

        if (has_converged) {
            std::cout << " Auto-refine: Refinement has converged, entering last iteration where two halves will be combined..." << std::endl;
            std::cout << " Auto-refine: The last iteration will use data to Nyquist frequency, which may take more CPU and RAM." << std::endl;
        }
    }

}

void MlOptimiser::setMetaDataSubset(long int first_part_id, long int last_part_id) {

    for (long int part_id_sorted = first_part_id, metadata_offset = 0; part_id_sorted <= last_part_id; part_id_sorted++) {

        const long int part_id = mydata.sorted_idx[part_id_sorted];
        for (int img_id = 0; img_id < mydata.numberOfImagesInParticle(part_id); img_id++, metadata_offset++) {

            const long int ori_img_id = mydata.particles[part_id].images[img_id].id;
            RFLOAT my_pixel_size = mydata.getImagePixelSize(part_id, img_id);

            mydata.MDimg.setValue(EMDL::ORIENT_PSI,                               direct::elem(exp_metadata, metadata_offset, METADATA_PSI),  ori_img_id);
            mydata.MDimg.setValue(EMDL::ORIENT_ORIGIN_X_ANGSTROM, my_pixel_size * direct::elem(exp_metadata, metadata_offset, METADATA_XOFF), ori_img_id);
            mydata.MDimg.setValue(EMDL::ORIENT_ORIGIN_Y_ANGSTROM, my_pixel_size * direct::elem(exp_metadata, metadata_offset, METADATA_YOFF), ori_img_id);
            if (mymodel.data_dim == 3)
            mydata.MDimg.setValue(EMDL::ORIENT_ORIGIN_Z_ANGSTROM, my_pixel_size * direct::elem(exp_metadata, metadata_offset, METADATA_ZOFF), ori_img_id);
            // SHWS: Upon request of Juha Huiskonen, 5 Apr 2016
            if (mymodel.ref_dim > 2) {
                mydata.MDimg.setValue(EMDL::ORIENT_ROT,  direct::elem(exp_metadata, metadata_offset, METADATA_ROT), ori_img_id);
                mydata.MDimg.setValue(EMDL::ORIENT_TILT, direct::elem(exp_metadata, metadata_offset, METADATA_TILT), ori_img_id);
            }
            mydata.MDimg.setValue(EMDL::PARTICLE_CLASS, (int) direct::elem(exp_metadata, metadata_offset, METADATA_CLASS), ori_img_id);
            mydata.MDimg.setValue(EMDL::PARTICLE_DLL,         direct::elem(exp_metadata, metadata_offset, METADATA_DLL),   ori_img_id);
            mydata.MDimg.setValue(EMDL::PARTICLE_PMAX,        direct::elem(exp_metadata, metadata_offset, METADATA_PMAX),  ori_img_id);
            mydata.MDimg.setValue(EMDL::PARTICLE_NR_SIGNIFICANT_SAMPLES, (int) direct::elem(exp_metadata, metadata_offset, METADATA_NR_SIGN), ori_img_id);
            mydata.MDimg.setValue(EMDL::IMAGE_NORM_CORRECTION, direct::elem(exp_metadata, metadata_offset, METADATA_NORM), ori_img_id);

            // For the moment, CTF, prior and transformation matrix info is NOT updated...
            RFLOAT prior_x = direct::elem(exp_metadata, metadata_offset, METADATA_XOFF_PRIOR);
            RFLOAT prior_y = direct::elem(exp_metadata, metadata_offset, METADATA_YOFF_PRIOR);
            if (prior_x < 999.0)
                mydata.MDimg.setValue(EMDL::ORIENT_ORIGIN_X_PRIOR_ANGSTROM, prior_x * my_pixel_size, ori_img_id);
            if (prior_y < 999.0)
                mydata.MDimg.setValue(EMDL::ORIENT_ORIGIN_Y_PRIOR_ANGSTROM, prior_y * my_pixel_size, ori_img_id);
            if (mymodel.data_dim == 3) {
                RFLOAT prior_z = direct::elem(exp_metadata, metadata_offset, METADATA_ZOFF_PRIOR);
                if (prior_z < 999.0)
                    mydata.MDimg.setValue(EMDL::ORIENT_ORIGIN_Z_PRIOR_ANGSTROM, prior_z * my_pixel_size, ori_img_id);
            }

            // For multi-body refinement
            if (mymodel.nr_bodies > 1) {
                for (int ibody = 0; ibody < mymodel.nr_bodies; ibody++) {
                    auto &mdt = mydata.MDbodies[ibody];
                    const int imx = METADATA_LINE_LENGTH_BEFORE_BODIES + ibody * METADATA_NR_BODY_PARAMS;
                    const RFLOAT rot  = direct::elem(exp_metadata, metadata_offset, imx + METADATA_ROT);
                    const RFLOAT tilt = direct::elem(exp_metadata, metadata_offset, imx + METADATA_TILT);
                    const RFLOAT psi  = direct::elem(exp_metadata, metadata_offset, imx + METADATA_PSI);
                    const RFLOAT xoff = direct::elem(exp_metadata, metadata_offset, imx + METADATA_XOFF);
                    const RFLOAT yoff = direct::elem(exp_metadata, metadata_offset, imx + METADATA_YOFF);
                    const RFLOAT zoff = direct::elem(exp_metadata, metadata_offset, imx + METADATA_ZOFF);
                    mdt.setValue(EMDL::ORIENT_ROT,               rot,                  ori_img_id);
                    mdt.setValue(EMDL::ORIENT_TILT,              tilt,                 ori_img_id);
                    mdt.setValue(EMDL::ORIENT_PSI,               psi,                  ori_img_id);
                    mdt.setValue(EMDL::ORIENT_ORIGIN_X_ANGSTROM, xoff * my_pixel_size, ori_img_id);
                    mdt.setValue(EMDL::ORIENT_ORIGIN_Y_ANGSTROM, yoff * my_pixel_size, ori_img_id);
                    if (mymodel.data_dim == 3)
                    mdt.setValue(EMDL::ORIENT_ORIGIN_Z_ANGSTROM, zoff * my_pixel_size, ori_img_id);
                }
            }

        }
    }
}

void MlOptimiser::getMetaAndImageDataSubset(long int first_part_id, long int last_part_id, bool do_also_imagedata) {

    // In case we're reading images here, only open stacks once and then read multiple images
    fImageHandler hFile;

    // Initialise filename strings if not reading imagedata here
    if (!do_also_imagedata) {
        exp_fn_img = "";
        exp_fn_ctf = "";
        exp_fn_recimg = "";
    }

    int nr_images = 0;
    for (long int part_id_sorted = first_part_id; part_id_sorted <= last_part_id; part_id_sorted++) {
        const long int part_id = mydata.sorted_idx[part_id_sorted];
        nr_images += mydata.numberOfImagesInParticle(part_id);
    }
    exp_metadata.initZeros(nr_images, METADATA_LINE_LENGTH_BEFORE_BODIES + mymodel.nr_bodies * METADATA_NR_BODY_PARAMS);

    // This assumes all images in first_part_id to last_part_id have the same image_size
    // If not, then do_also_imagedata will not work! Also warn during intialiseGeneral!
    int common_image_size = mydata.getOpticsImageSize(mydata.getOpticsGroup(first_part_id, 0));

    if (do_also_imagedata) {
        if (mymodel.data_dim == 3) {
            if (nr_images > 1)
                REPORT_ERROR("MlOptimiser::getMetaAndImageDataSubset ERROR: cannot get multiple images for 3D data!");
            exp_imagedata.resize(common_image_size, common_image_size, ((has_converged && do_use_reconstruct_images) + 1 + do_ctf_correction) * common_image_size);
        } else {
            exp_imagedata.resize(common_image_size, common_image_size, ((has_converged && do_use_reconstruct_images) + 1) * nr_images);
        }
    }

    for (long int part_id_sorted = first_part_id, metadata_offset = 0; part_id_sorted <= last_part_id; part_id_sorted++) {

        const long int part_id = mydata.sorted_idx[part_id_sorted];
        for (int img_id = 0; img_id < mydata.numberOfImagesInParticle(part_id); img_id++, metadata_offset++) {

            const long int ori_img_id = mydata.particles[part_id].images[img_id].id;
            RFLOAT my_pixel_size = mydata.getImagePixelSize(part_id, img_id);
            int my_image_size = mydata.getOpticsImageSize(mydata.getOpticsGroup(part_id, img_id));

            // Get the image names from the MDimg table
            FileName fn_img = [&] () {
                try {
                    return mydata.getImageNameOnScratch(part_id, img_id);
                } catch (const char *errmsg) {
                    return (FileName) mydata.MDimg.getValue<std::string>(EMDL::IMAGE_NAME, ori_img_id);
                }
            }();

            FileName fn_ctf = "";
            if (mymodel.data_dim == 3 && do_ctf_correction) {
                // Also read the CTF image from disc
                try {
                     fn_ctf = mydata.getImageNameOnScratch(part_id, img_id, true);
                } catch (const char *errmsg) { try {
                    fn_ctf = mydata.MDimg.getValue<std::string>(EMDL::CTF_IMAGE, ori_img_id);
                } catch (const char *errmsg) {
                    REPORT_ERROR("MlOptimiser::getMetaAndImageDataSubset ERROR: cannot find rlnCtfImage for 3D CTF correction!");
                } }
            }

            FileName fn_rec_img = "";
            if (has_converged && do_use_reconstruct_images) {
                fn_rec_img = mydata.MDimg.getValue<std::string>(EMDL::IMAGE_RECONSTRUCT_NAME, ori_img_id);
            }

            if (do_also_imagedata) {
                if (my_image_size != common_image_size)
                    REPORT_ERROR("ERROR: non-parallel disc I/O is not supported when images with different box sizes are present in the data set.");

                // First read the image from disc or get it from the preread images in the mydata structure
                Image<RFLOAT> img, rec_img;
                if (do_preread_images) {
                    img() = mydata.particles[part_id].images[img_id].img;
                } else {
                    // only open new stacks
                    FileName fn_stack, fn_open_stack = "";
                    long int dump;
                    fn_img.decompose(dump, fn_stack);
                    if (fn_stack != fn_open_stack) {
                        hFile.openFile(fn_stack, WRITE_READONLY);
                        fn_open_stack = fn_stack;
                    }
                    img.readFromOpenFile(fn_img, hFile, -1, false);
                    img().setXmippOrigin();
                }
                if (Xsize(img()) != Xsize(exp_imagedata) || Ysize(img()) != Ysize(exp_imagedata)) {
                    std::cerr << " fn_img= " << fn_img << " Xsize(img())= " << Xsize(img()) << " Ysize(img())= " << Ysize(img()) << std::endl;
                    std::cerr << " while Xsize(exp_imagedata)= " << Xsize(exp_imagedata) << " and Ysize(exp_imagedata)= " << Ysize(exp_imagedata) << std::endl;
                    REPORT_ERROR("MlOptimiser::getMetaAndImageDataSubset ERROR: incorrect image size");
                }
                if (has_converged && do_use_reconstruct_images) {
                    rec_img.read(fn_rec_img);
                    if (Xsize(rec_img()) != Xsize(exp_imagedata) || Ysize(rec_img()) != Ysize(exp_imagedata)) {
                        std::cerr << " fn_rec_img= " << fn_rec_img << " Xsize(rec_img())= " << Xsize(rec_img()) << " Ysize(rec_img())= " << Ysize(rec_img()) << std::endl;
                        REPORT_ERROR("MlOptimiser::getMetaAndImageDataSubset ERROR: incorrect reconstruct_image size");
                    }
                }
                if (mymodel.data_dim == 3) {

                    for (long int k = 0; k < Zsize(img()); k++)
                    for (long int j = 0; j < Ysize(img()); j++)
                    for (long int i = 0; i < Xsize(img()); i++) {
                        direct::elem(exp_imagedata, i, j, k) = direct::elem(img(), i, j, k);
                    }

                    if (do_ctf_correction) {
                        img.read(fn_ctf);
                        for (long int k = 0; k < Zsize(img()); k++)
                        for (long int j = 0; j < Ysize(img()); j++)
                        for (long int i = 0; i < Xsize(img()); i++) {
                            direct::elem(exp_imagedata, i, j, my_image_size + k) = direct::elem(img(), i, j, k);
                        }
                    }

                    if (has_converged && do_use_reconstruct_images) {
                        int offset = do_ctf_correction ? 2 * my_image_size : my_image_size;
                        for (long int k = 0; k < Zsize(img()); k++)
                        for (long int j = 0; j < Ysize(img()); j++)
                        for (long int i = 0; i < Xsize(img()); i++) {
                            direct::elem(exp_imagedata, i, j, offset + k) = direct::elem(rec_img(), i, j, k);
                        }
                    }

                } else {
                    for (long int j = 0; j < Ysize(img()); j++)
                    for (long int i = 0; i < Xsize(img()); i++) {
                        direct::elem(exp_imagedata, i, j, metadata_offset) = direct::elem(img(), i, j);
                    }

                    if (has_converged && do_use_reconstruct_images) {
                        for (long int j = 0; j < Ysize(rec_img()); j++)
                        for (long int i = 0; i < Xsize(rec_img()); i++) {
                            direct::elem(exp_imagedata, i, j, metadata_offset) = direct::elem(rec_img(), i, j);
                        }
                    }
                }
            } else {
                exp_fn_img += fn_img + "\n";
                if (!fn_ctf.empty())
                exp_fn_ctf += fn_ctf + "\n";
                if (!fn_rec_img.empty())
                exp_fn_recimg += fn_rec_img + "\n";
            }

            // Now get the metadata
            direct::elem(exp_metadata, metadata_offset, METADATA_ROT)  = mydata.MDimg.getValue<RFLOAT>(EMDL::ORIENT_ROT,  ori_img_id);
            direct::elem(exp_metadata, metadata_offset, METADATA_TILT) = mydata.MDimg.getValue<RFLOAT>(EMDL::ORIENT_TILT, ori_img_id);
            if (mymodel.data_dim == 3)
            direct::elem(exp_metadata, metadata_offset, METADATA_PSI)  = mydata.MDimg.getValue<RFLOAT>(EMDL::ORIENT_PSI,  ori_img_id);
            RFLOAT xoff_A = mydata.MDimg.getValue<RFLOAT>(EMDL::ORIENT_ORIGIN_X_ANGSTROM, ori_img_id);
            RFLOAT yoff_A = mydata.MDimg.getValue<RFLOAT>(EMDL::ORIENT_ORIGIN_Y_ANGSTROM, ori_img_id);
            RFLOAT zoff_A = mydata.MDimg.getValue<RFLOAT>(EMDL::ORIENT_ORIGIN_Z_ANGSTROM, ori_img_id);
            direct::elem(exp_metadata, metadata_offset, METADATA_XOFF) = xoff_A / my_pixel_size;
            direct::elem(exp_metadata, metadata_offset, METADATA_YOFF) = yoff_A / my_pixel_size;
            if (mymodel.data_dim == 3)
            direct::elem(exp_metadata, metadata_offset, METADATA_ZOFF) = zoff_A / my_pixel_size;

            direct::elem(exp_metadata, metadata_offset, METADATA_CLASS) = mydata.MDimg.getValue<RFLOAT>(EMDL::PARTICLE_CLASS, ori_img_id);
            direct::elem(exp_metadata, metadata_offset, METADATA_DLL)   = mydata.MDimg.getValue<RFLOAT>(EMDL::PARTICLE_DLL,   ori_img_id);
            direct::elem(exp_metadata, metadata_offset, METADATA_PMAX)  = mydata.MDimg.getValue<RFLOAT>(EMDL::PARTICLE_PMAX,  ori_img_id);

            // 5 July 2017: we do not need EMDL::PARTICLE_NR_SIGNIFICANT_SAMPLES for calculations. Send randomsubset instead!
            direct::elem(exp_metadata, metadata_offset, METADATA_NR_SIGN) = do_split_random_halves ?
                mydata.MDimg.getValue<RFLOAT>(EMDL::PARTICLE_RANDOM_SUBSET, ori_img_id) : mydata.MDimg.getValue<RFLOAT>(EMDL::PARTICLE_NR_SIGNIFICANT_SAMPLES, ori_img_id);
            try {
                direct::elem(exp_metadata, metadata_offset, METADATA_NORM) = mydata.MDimg.getValue<RFLOAT>(EMDL::IMAGE_NORM_CORRECTION, ori_img_id);
            } catch (const char* errmsg) {
                direct::elem(exp_metadata, metadata_offset, METADATA_NORM) = 1.0;
            }

            // If the priors are not set, set them to 999.0.
            const RFLOAT rot_prior  = tryget<RFLOAT>(mydata.MDimg, EMDL::ORIENT_ROT_PRIOR,  999.0, ori_img_id);
            const RFLOAT tilt_prior = tryget<RFLOAT>(mydata.MDimg, EMDL::ORIENT_TILT_PRIOR, 999.0, ori_img_id);
            const RFLOAT psi_prior  = tryget<RFLOAT>(mydata.MDimg, EMDL::ORIENT_PSI_PRIOR,  999.0, ori_img_id);
            direct::elem(exp_metadata, metadata_offset, METADATA_ROT_PRIOR)  = rot_prior;
            direct::elem(exp_metadata, metadata_offset, METADATA_TILT_PRIOR) = tilt_prior;
            if (mymodel.data_dim == 3)
            direct::elem(exp_metadata, metadata_offset, METADATA_PSI_PRIOR)  = psi_prior;

            // Not ideal: if tryget encounters an error, we will end up doing 999.0 * my_pixel_size / my_pixel size.
            xoff_A = tryget<RFLOAT>(mydata.MDimg, EMDL::ORIENT_ORIGIN_X_PRIOR_ANGSTROM, 999.0 * my_pixel_size, ori_img_id);
            yoff_A = tryget<RFLOAT>(mydata.MDimg, EMDL::ORIENT_ORIGIN_Y_PRIOR_ANGSTROM, 999.0 * my_pixel_size, ori_img_id);
            zoff_A = tryget<RFLOAT>(mydata.MDimg, EMDL::ORIENT_ORIGIN_Z_PRIOR_ANGSTROM, 999.0 * my_pixel_size, ori_img_id);
            direct::elem(exp_metadata, metadata_offset, METADATA_XOFF_PRIOR) = xoff_A / my_pixel_size;
            direct::elem(exp_metadata, metadata_offset, METADATA_YOFF_PRIOR) = yoff_A / my_pixel_size;
            if (mymodel.data_dim == 3)
            direct::elem(exp_metadata, metadata_offset, METADATA_ZOFF_PRIOR) = zoff_A / my_pixel_size;

            const RFLOAT psi_prior_flip_ratio = tryget<RFLOAT>(mydata.MDimg, EMDL::ORIENT_PSI_PRIOR_FLIP_RATIO, 999.0, ori_img_id);
            direct::elem(exp_metadata, metadata_offset, METADATA_PSI_PRIOR_FLIP_RATIO) = psi_prior_flip_ratio;

            // The following per-particle parameters are passed around through metadata
            // Note beamtilt is no longer part of this: it is now in the optics group
            if (do_ctf_correction) {
                long int mic_id = mydata.getMicrographId(part_id, img_id);  // Unused

                const RFLOAT DeltafU          = tryget<RFLOAT>(mydata.MDimg, EMDL::CTF_DEFOCUSU,      0.0,     ori_img_id);
                const RFLOAT DeltafV          = tryget<RFLOAT>(mydata.MDimg, EMDL::CTF_DEFOCUSV,      DeltafV, ori_img_id);
                const RFLOAT azimuthal_angle  = tryget<RFLOAT>(mydata.MDimg, EMDL::CTF_DEFOCUS_ANGLE, 0.0,     ori_img_id);
                const RFLOAT Bfac             = tryget<RFLOAT>(mydata.MDimg, EMDL::CTF_BFACTOR,       0.0,     ori_img_id);
                const RFLOAT kfac             = tryget<RFLOAT>(mydata.MDimg, EMDL::CTF_SCALEFACTOR,   1.0,     ori_img_id);
                const RFLOAT phase_shift      = tryget<RFLOAT>(mydata.MDimg, EMDL::CTF_PHASESHIFT,    0.0,     ori_img_id);

                direct::elem(exp_metadata, metadata_offset, METADATA_CTF_DEFOCUS_U)     = DeltafU;
                direct::elem(exp_metadata, metadata_offset, METADATA_CTF_DEFOCUS_V)     = DeltafV;
                direct::elem(exp_metadata, metadata_offset, METADATA_CTF_DEFOCUS_ANGLE) = azimuthal_angle;
                direct::elem(exp_metadata, metadata_offset, METADATA_CTF_BFACTOR)       = Bfac;
                direct::elem(exp_metadata, metadata_offset, METADATA_CTF_KFACTOR)       = kfac;
                direct::elem(exp_metadata, metadata_offset, METADATA_CTF_PHASE_SHIFT)   = phase_shift;

            }

            // For multi-body refinement
            if (mymodel.nr_bodies > 1) {
                for (int ibody = 0; ibody < mymodel.nr_bodies; ibody++) {
                    const RFLOAT rot  = tryget<RFLOAT>(mydata.MDbodies[ibody], EMDL::ORIENT_ROT,               0.0, ori_img_id);
                    const RFLOAT tilt = tryget<RFLOAT>(mydata.MDbodies[ibody], EMDL::ORIENT_TILT,              0.0, ori_img_id);
                    const RFLOAT psi  = tryget<RFLOAT>(mydata.MDbodies[ibody], EMDL::ORIENT_PSI,               0.0, ori_img_id);
                    const RFLOAT xoff = tryget<RFLOAT>(mydata.MDbodies[ibody], EMDL::ORIENT_ORIGIN_X_ANGSTROM, 0.0, ori_img_id);
                    const RFLOAT yoff = tryget<RFLOAT>(mydata.MDbodies[ibody], EMDL::ORIENT_ORIGIN_Y_ANGSTROM, 0.0, ori_img_id);
                    const RFLOAT zoff = tryget<RFLOAT>(mydata.MDbodies[ibody], EMDL::ORIENT_ORIGIN_Z_ANGSTROM, 0.0, ori_img_id);
                    const int imx = METADATA_LINE_LENGTH_BEFORE_BODIES + ibody * METADATA_NR_BODY_PARAMS;
                    direct::elem(exp_metadata, metadata_offset, imx + METADATA_ROT)  = rot;
                    direct::elem(exp_metadata, metadata_offset, imx + METADATA_TILT) = tilt;
                    direct::elem(exp_metadata, metadata_offset, imx + METADATA_PSI)  = psi;
                    direct::elem(exp_metadata, metadata_offset, imx + METADATA_XOFF) = xoff / my_pixel_size;
                    direct::elem(exp_metadata, metadata_offset, imx + METADATA_YOFF) = yoff / my_pixel_size;
                    direct::elem(exp_metadata, metadata_offset, imx + METADATA_ZOFF) = zoff / my_pixel_size;
                }
            }
        }
    }
}
