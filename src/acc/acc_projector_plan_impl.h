#include "src/acc/acc_projector_plan.h"
#include "src/acc/utilities.h"
#include "src/time.h"

// #define PP_TIMING
#ifdef PP_TIMING
    Timer timer;
    int TIMING_TOP        = timer.setNew("setup");
    int TIMING_SAMPLING   =	timer.setNew(" sampling");
    int TIMING_PRIOR      = timer.setNew("  prior");
    int TIMING_PROC_CALC  = timer.setNew("  procCalc");
    int TIMING_PROC       = timer.setNew("  proc");
    int TIMING_GEN        = timer.setNew("   genOri");
    int TIMING_PERTURB    = timer.setNew("   perturb");
    int TIMING_EULERS     = timer.setNew(" eulers");
    #define TIMING_TIC(id) timer.tic(id)
    #define TIMING_TOC(id) timer.toc(id)
#else
    #define TIMING_TIC(id)
    #define TIMING_TOC(id)
#endif

/// HACK: Imitate a context manager.
#define TICTOC(n, block) TIMING_TIC(n); block; TIMING_TOC(n);

void getOrientations(HealpixSampling &sampling, long int idir, long int ipsi, int oversampling_order,
    std::vector<RFLOAT > &my_rot, std::vector<RFLOAT > &my_tilt, std::vector<RFLOAT > &my_psi,
    std::vector<int> &pointer_dir_nonzeroprior, std::vector<RFLOAT> &directions_prior,
    std::vector<int> &pointer_psi_nonzeroprior, std::vector<RFLOAT> &psi_prior
) {
    my_rot.clear();
    my_tilt.clear();
    my_psi.clear();
    long int my_idir, my_ipsi;
    if (pointer_dir_nonzeroprior.size() > idir && pointer_psi_nonzeroprior.size() > ipsi) {
        // nonzeroprior vectors have been initialised, so use priors!
        my_idir = pointer_dir_nonzeroprior[idir];
        my_ipsi = pointer_psi_nonzeroprior[ipsi];
    } else {
        // no priors
        my_idir = idir;
        my_ipsi = ipsi;
    }

    if (oversampling_order == 0) {
        my_rot.push_back(sampling.rot_angles[my_idir]);
        my_tilt.push_back(sampling.tilt_angles[my_idir]);
        my_psi.push_back(sampling.psi_angles[my_ipsi]);
    } else if (!sampling.is_3D) {
        // for 2D sampling, only push back oversampled psi rotations
        sampling.pushbackOversampledPsiAngles(my_ipsi, oversampling_order, 0.0, 0.0, my_rot, my_tilt, my_psi);
    } else {
        // Set up oversampled grid for 3D sampling
        Healpix_Base HealPixOver(oversampling_order + sampling.healpix_order, NEST);
        int fact = HealPixOver.Nside() / sampling.healpix_base.Nside();
        // Get x, y and face for the original, coarse grid
        long int ipix = sampling.directions_ipix[my_idir];
        int x, y, face;
        sampling.healpix_base.nest2xyf(ipix, x, y, face);
        // Loop over the oversampled Healpix pixels on the fine grid
        for (int j = fact * y; j < fact * y + fact; ++j) {
            for (int i = fact * x; i < fact * x + fact; ++i) {
                long int overpix = HealPixOver.xyf2nest(i, j, face);
                // this one always has to be double (also for SINGLE_PRECISION CALCULATIONS) for call to external library
                double zz, phi;
                HealPixOver.pix2ang_z_phi(overpix, zz, phi);
                RFLOAT rot  = degrees(phi);
                RFLOAT tilt = degrees(acos(zz));

                // The geometrical considerations about the symmetry below require that rot = [-180,180] and tilt [0,180]
                sampling.checkDirection(rot, tilt);

                sampling.pushbackOversampledPsiAngles(my_ipsi, oversampling_order, rot, tilt, my_rot, my_tilt, my_psi);
            }
        }
    }
}

void AccProjectorPlan::setup(
    HealpixSampling &sampling,
    std::vector<RFLOAT> &directions_prior, std::vector<RFLOAT> &psi_prior,
    std::vector<int> &pointer_dir_nonzeroprior, std::vector<int> &pointer_psi_nonzeroprior,
    MultidimArray<bool> *Mcoarse_significant,
    std::vector<RFLOAT> &pdf_class,
    std::vector<MultidimArray<RFLOAT> > &pdf_direction,
    unsigned long nr_dir, unsigned long nr_psi,
    unsigned long idir_min, unsigned long idir_max,
    unsigned long ipsi_min, unsigned long ipsi_max,
    unsigned long itrans_min, unsigned long itrans_max,
    unsigned long current_oversampling, unsigned long nr_oversampled_rot,
    unsigned iclass,
    bool coarse, bool inverseMatrix, bool do_skip_align, bool do_skip_rotate,
    int orientational_prior_mode,
    Matrix<RFLOAT> &L_, Matrix<RFLOAT> &R_
) {
    TICTOC(TIMING_TOP, ({

    std::vector<RFLOAT> oversampled_rot, oversampled_tilt, oversampled_psi;

    AccPtr<XFLOAT> alphas =  eulers.make<XFLOAT>(nr_dir * nr_psi * nr_oversampled_rot * 9);
    AccPtr<XFLOAT> betas =   eulers.make<XFLOAT>(nr_dir * nr_psi * nr_oversampled_rot * 9);
    AccPtr<XFLOAT> gammas =  eulers.make<XFLOAT>(nr_dir * nr_psi * nr_oversampled_rot * 9);
    AccPtr<XFLOAT> perturb = eulers.make<XFLOAT>((size_t) 9);
    AccPtr<XFLOAT> adjustL = eulers.make<XFLOAT>((size_t) 9);
    AccPtr<XFLOAT> adjustR = eulers.make<XFLOAT>((size_t) 9);

    alphas.hostAlloc();
    betas .hostAlloc();
    gammas.hostAlloc();

    eulers.free();
    eulers.setSize(nr_dir * nr_psi * nr_oversampled_rot * 9);
    eulers.hostAlloc();

    iorientclasses.free();
    iorientclasses.setSize(nr_dir * nr_psi * nr_oversampled_rot);
    iorientclasses.hostAlloc();

    orientation_num = 0;

    auto L = Matrix<RFLOAT>::identity(3);
    auto R = Matrix<RFLOAT>::identity(3);

    bool doL = false, doR = false;
    RFLOAT myperturb = 0.0;

    if (L_.shape() == L.shape()) {
        doL = true;
        L = L.matmul(L_);
    }

    if (abs(sampling.random_perturbation) > 0.0) {
        myperturb = sampling.random_perturbation * sampling.getAngularSampling();
        if (sampling.is_3D) {
            R = Euler::angles2matrix(myperturb, myperturb, myperturb);
        }
        doR = true;
    }

    if (R_.shape() == R.shape()) {
        doR = true;
        R = R.matmul(R_);
    }

    TICTOC(TIMING_SAMPLING, ({

    for (long int idir = idir_min, iorient = 0; idir <= idir_max; idir++) {
        for (long int ipsi = ipsi_min, ipart = 0; ipsi <= ipsi_max; ipsi++, iorient++) {
            long int iorientclass = iclass * nr_dir * nr_psi + iorient;

            RFLOAT pdf_orientation;
            TICTOC(TIMING_PRIOR, ({
            // Get prior for this direction and skip calculation if prior == 0
            if (do_skip_align || do_skip_rotate) {
                pdf_orientation = pdf_class[iclass];
            } else if (orientational_prior_mode == NOPRIOR) {
                pdf_orientation = pdf_direction[iclass][idir];
            } else {
                pdf_orientation = directions_prior[idir] * psi_prior[ipsi];
            }
            }));

            // In the first pass, always proceed
            // In the second pass, check whether one of the translations for this orientation of any of the particles had a significant weight in the first pass
            // if so, proceed with projecting the reference in that direction

            bool do_proceed = false;

            TICTOC(TIMING_PROC_CALC, ({
            if (coarse && pdf_orientation > 0.0) {
                do_proceed = true;
            } else if (pdf_orientation > 0.0) {
                long int nr_trans = itrans_max - itrans_min + 1;
                for (long int ipart = 0; ipart < Ysize(*Mcoarse_significant); ipart++) {
                    long int ihidden = iorient * nr_trans;
                    for (long int itrans = itrans_min; itrans <= itrans_max; itrans++, ihidden++) {
                        if (direct::elem(*Mcoarse_significant, ipart, ihidden)) {
                            do_proceed = true;
                            break;
                        }
                    }
                }
            }
            }));

            TICTOC(TIMING_PROC, ({
            if (do_proceed) {
                // Now get the oversampled (rot, tilt, psi) triplets
                // This will be only the original (rot, tilt, psi) triplet in the first pass (sp.current_oversampling == 0)
                TICTOC(TIMING_GEN, ({
                getOrientations(
                    sampling, idir, ipsi, current_oversampling, oversampled_rot, oversampled_tilt, oversampled_psi,
                    pointer_dir_nonzeroprior, directions_prior, pointer_psi_nonzeroprior, psi_prior
                );
                }));

                // Loop over all oversampled orientations (only a single one in the first pass)
                for (long int iover_rot = 0; iover_rot < nr_oversampled_rot; iover_rot++, ipart++) {
                    if (sampling.is_3D) {
                        alphas.getHostPtr()[orientation_num] = oversampled_rot [iover_rot];
                        betas .getHostPtr()[orientation_num] = oversampled_tilt[iover_rot];
                        gammas.getHostPtr()[orientation_num] = oversampled_psi [iover_rot];
                    } else {
                        alphas.getHostPtr()[orientation_num] = oversampled_psi [iover_rot] + myperturb;
                    }

                    iorientclasses.getHostPtr()[orientation_num] = iorientclass;
                    orientation_num++;
                }
            }
            }));
        }
    }
    }));

    iorientclasses.resizeHost(orientation_num);
    iorientclasses.putOnDevice();

    eulers.resizeHost(orientation_num * 9);
    eulers.deviceAlloc();

    alphas.resizeHost(orientation_num);
    alphas.putOnDevice();

    if (sampling.is_3D) {
        betas.resizeHost(orientation_num);
        betas.putOnDevice();
        gammas.resizeHost(orientation_num);
        gammas.putOnDevice();
    }

    if (doL) {
        adjustL.hostAlloc();
        std::copy_n(L.data(), 9, adjustL.getHostPtr());
        adjustL.putOnDevice();
    }

    if (doR) {
        adjustR.hostAlloc();
        std::copy_n(R.data(), 9, adjustR.getHostPtr());
        adjustR.putOnDevice();
    }

    const int grid_size = ceil((float) orientation_num / (float) BLOCK_SIZE);

    if (sampling.is_3D) {
        AccUtilities::acc_make_eulers_3D<acc::type>(
            grid_size, BLOCK_SIZE, eulers.getStream(),
            alphas.getAccPtr(), betas.getAccPtr(), gammas.getAccPtr(), eulers.getAccPtr(),
            orientation_num,
            doL ? adjustL.getAccPtr() : nullptr,
            doR ? adjustR.getAccPtr() : nullptr,
            doL, doR, inverseMatrix
        );
    } else {
        AccUtilities::acc_make_eulers_2D<acc::type>(
            grid_size, BLOCK_SIZE, eulers.getStream(),
            alphas.getAccPtr(), eulers.getAccPtr(),
            orientation_num, inverseMatrix
        );
    }

    }));
}

void AccProjectorPlan::printTo(std::ostream &os) {
    // print
    os << "orientation_num = " << orientation_num << std::endl;
    os << "iorientclasses.getSize() = " << iorientclasses.getSize() << std::endl;
    os << std::endl << "iorientclasses\tiover_rots\teulers" << std::endl;

    for (int i = 0; i < iorientclasses.getSize(); i++) {
        os << iorientclasses.getHostPtr()[i] << "\t\t" << "\t";
        for (int j = 0; j < 9; j++)
            os << eulers.getHostPtr()[i * 9 + j] << "\t";
        os << std::endl;
    }
}

void AccProjectorPlan::clear() {
    orientation_num = 0;
    iorientclasses.free();
    iorientclasses.setSize(0);
    eulers.free();
    eulers.setSize(0);
    #ifdef PP_TIMING
    timer.printTimes(false);
    #endif
}
