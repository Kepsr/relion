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

#include <src/projector.h>
#include <src/backprojector.h>
#include <src/fftw.h>
#include <src/args.h>
#include <src/ctf.h>
#include <src/strings.h>
#include <src/funcs.h>
#include <src/memory.h>
#include <src/euler.h>
#include <src/time.h>
#include <src/metadata_table.h>
#include <src/ml_model.h>
#include <src/exp_model.h>
#include <src/healpix_sampling.h>
#include <src/jaz/obs_model.h>


class project_parameters {

    public:

    FileName fn_map, fn_ang, fn_out, fn_img, fn_model, fn_sym, fn_mask, fn_ang_simulate;
    RFLOAT rot, tilt, psi, xoff, yoff, zoff, angpix, maxres, stddev_white_noise, particle_diameter, ana_prob_range, ana_prob_step, sigma_offset;
    int padding_factor;
    int r_max, r_min_nn, interpolator, nr_uniform;
    bool do_only_one, do_ctf, do_ctf2, ctf_phase_flipped, do_ctf_intact_1st_peak, do_timing, do_add_noise, do_subtract_exp, do_ignore_particle_name, do_3d_rot;
    bool do_simulate;
    RFLOAT simulate_SNR;
    // I/O Parser
    IOParser parser;
    MlModel model;
    ObservationModel obsModel;

    void usage() { parser.writeUsage(std::cerr); }

    void read(int argc, char **argv) {
        parser.setCommandLine(argc, argv);

        int general_section = parser.addSection("Options");
        fn_map = parser.getOption("--i", "Input map to be projected");
        fn_out = parser.getOption("--o", "Rootname for output projections", "proj");
        do_ctf = parser.checkOption("--ctf", "Apply CTF to reference projections");
        ctf_phase_flipped = parser.checkOption("--ctf_phase_flip", "Flip phases of the CTF in the output projections");
        do_ctf_intact_1st_peak = parser.checkOption("--ctf_intact_first_peak", "Ignore CTFs until their first peak?");
        angpix = textToFloat(parser.getOption("--angpix", "Pixel size (in Angstroms)", "-1"));
        fn_mask = parser.getOption("--mask", "Mask that will be applied to the input map prior to making projections", "");
        fn_ang = parser.getOption("--ang", "Particle STAR file with orientations and CTF for multiple projections (if None, assume single projection)", "None");
        nr_uniform = textToInteger(parser.getOption("--nr_uniform", " OR get this many random samples from a uniform angular distribution", "-1"));
        sigma_offset = textToFloat(parser.getOption("--sigma_offset", "Apply Gaussian errors with this stddev to the XY-offsets", "0"));
        rot = textToFloat(parser.getOption("--rot", "First Euler angle (for a single projection)", "0"));
        tilt = textToFloat(parser.getOption("--tilt", "Second Euler angle (for a single projection)", "0"));
        psi = textToFloat(parser.getOption("--psi", "Third Euler angle (for a single projection)", "0"));
        xoff = textToFloat(parser.getOption("--xoff", "Origin X-offsets (in pixels) (for a single projection)", "0"));
        yoff = textToFloat(parser.getOption("--yoff", "Origin Y-offsets (in pixels) (for a single projection)", "0"));
        zoff = textToFloat(parser.getOption("--zoff", "Origin Z-offsets (in pixels) (for a single 3D rotation)", "0"));
        do_add_noise = parser.checkOption("--add_noise", "Add noise to the output projections (only with --ang)");
        stddev_white_noise = textToFloat(parser.getOption("--white_noise", "Standard deviation of added white Gaussian noise", "0"));
        fn_model = parser.getOption("--model_noise", "Model STAR file with power spectra for coloured Gaussian noise", "");
        do_subtract_exp = parser.checkOption("--subtract_exp", "Subtract projections from experimental images (in --ang)");
        do_ignore_particle_name = parser.checkOption("--ignore_particle_name", "Ignore the rlnParticleName column (in --ang)");
        do_only_one = (fn_ang == "None" && nr_uniform < 0);
        do_3d_rot = parser.checkOption("--3d_rot", "Perform 3D rotations instead of projection into 2D images");
        do_simulate = parser.checkOption("--simulate", "Simulate data with known ground-truth by subtracting signal and adding projection in random orientation.");
        simulate_SNR = textToFloat(parser.getOption("--adjust_simulation_SNR", "Relative SNR compared to input images for realistic simulation of data", "1."));
        fn_ang_simulate = parser.getOption("--ang_simulate", "STAR file with orientations for projections of realistic simulations (random from --ang STAR file by default)", "");

        maxres = textToFloat(parser.getOption("--maxres", "Maximum resolution (in Angstrom) to consider in Fourier space (default Nyquist)", "-1"));
        padding_factor = textToInteger(parser.getOption("--pad", "Padding factor", "2"));
        do_ctf2 = parser.checkOption("--ctf2", "Apply CTF*CTF to reference projections");
        interpolator = parser.checkOption("--NN", "Use nearest-neighbour instead of linear interpolation") ?
            NEAREST_NEIGHBOUR : TRILINEAR;

        // Hidden
        r_min_nn = textToInteger(getParameter(argc, argv, "--r_min_nn", "10"));

        if (do_simulate) {
            do_ctf = true;
        }

        // Check for errors in the command-line option
        if (parser.checkForErrors())
            REPORT_ERROR("Errors encountered on the command line (see above), exiting...");
    }

    void project() {
        MetaDataTable DFo, MDang, MDang_sim;
        Matrix2D<RFLOAT> A3D;
        FileName fn_expimg;

        MultidimArray<Complex > F3D, F2D, Fexpimg;
        MultidimArray<RFLOAT> Fctf, dummy;
        Image<RFLOAT> vol, img, expimg;
        FourierTransformer transformer, transformer_expimg;

        std::cout << " Reading map: " << fn_map << std::endl;
        vol.read(fn_map);
        std::cout << " Done reading map!" << std::endl;

        if (fn_mask != "") {
            Image<RFLOAT> msk;
            msk.read(fn_mask);
            if (!msk().sameShape(vol()))
                REPORT_ERROR("project ERROR: mask and map have different sizes!");
            FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(vol())
            DIRECT_MULTIDIM_ELEM(vol(), n) *= DIRECT_MULTIDIM_ELEM(msk(), n);
        }

        if (nr_uniform > 0) {
            std::cout << " Generating " << nr_uniform << " projections taken randomly from a uniform angular distribution ..." << std::endl;
            MDang.clear();
            randomize_random_generator();
            for (long int i = 0; i < nr_uniform; i++) {
                RFLOAT rot = rnd_unif() * 360.0;
                // tilt will be a random angle in the interval [0.0, 180.0]
                // sin(tilt) (which will thus be in the interval [0.0, 1.0]) 
                // must be greater than a random number from that same interval [0.0, 1.0]
                RFLOAT tilt; 
                do {
                    tilt = rnd_unif() * 180.0;
                } while (sin(radians(tilt)) <= rnd_unif());
                RFLOAT psi = rnd_unif() * 360.0;
                RFLOAT xoff = rnd_gaus(0.0, sigma_offset);
                RFLOAT yoff = rnd_gaus(0.0, sigma_offset);
                MDang.addObject();
                MDang.setValue(EMDL::ORIENT_ROT, rot);
                MDang.setValue(EMDL::ORIENT_TILT, tilt);
                MDang.setValue(EMDL::ORIENT_PSI, psi);
                MDang.setValue(EMDL::ORIENT_ORIGIN_X, xoff);
                MDang.setValue(EMDL::ORIENT_ORIGIN_Y, yoff);
                MDang.setValue(EMDL::IMAGE_OPTICS_GROUP, 1);
            }

            std::cout << " Setting default values for optics table, though CTFs are not used in the projections ... " << std::endl;
            MetaDataTable MDopt;
            MDopt.addObject();
            MDopt.setValue(EMDL::IMAGE_OPTICS_GROUP, 1);
            std::string name = "optics1";
            MDopt.setValue(EMDL::IMAGE_OPTICS_GROUP_NAME, name);
            MDopt.setValue(EMDL::CTF_VOLTAGE, 300.0);
            MDopt.setValue(EMDL::CTF_CS, 2.7);
            angpix = vol.MDMainHeader.getValue<RFLOAT>(EMDL::IMAGE_SAMPLINGRATE_X);
            MDopt.setValue(EMDL::IMAGE_PIXEL_SIZE, angpix);
            MDopt.setValue(EMDL::IMAGE_SIZE, Xsize(vol()));
            int mydim = do_3d_rot ? 3 : 2;
            MDopt.setValue(EMDL::IMAGE_DIMENSIONALITY, mydim);

            obsModel = ObservationModel(MDopt);
        } else if (!do_only_one) {
            std::cout << " Reading STAR file with all angles " << fn_ang << std::endl;
            ObservationModel::loadSafely(fn_ang, obsModel, MDang);
            std::cout << " Done reading STAR file!" << std::endl;

            if (do_simulate && fn_ang_simulate != "") {
                std::cout << " Reading STAR file with angles for simulated images " << fn_ang << std::endl;
                MDang_sim.read(fn_ang_simulate);
                std::cout << " Done reading STAR file with angles for simulated images!" << std::endl;
                if (MDang_sim.numberOfObjects() < MDang.numberOfObjects()) {
                    REPORT_ERROR("ERROR: STAR file with angles for simulated images has fewer entries than the input STAR file with all angles.");
                }
            }
        }

        if (angpix < 0.0) {
            if (!do_only_one) {
                // Get angpix from the first optics group in the obsModel
                angpix = obsModel.getPixelSize(0);
                std::cout << " + Using pixel size from the first optics group in the --ang STAR file: " << angpix << std::endl;
            } else {
                angpix = vol.samplingRateX();
                std::cerr << "WARNING: The pixel size (--angpix) was not specified." << std::endl;
                std::cerr << "         The value in the input image header (= " << angpix << ") is used instead." << std::endl;
            }
        }

        // Now that we have the size of the volume, check r_max
        r_max = maxres < 0.0 ? Xsize(vol()) : ceil(Xsize(vol()) * angpix / maxres);

        // Set right size of F2D and initialize to zero
        if (do_3d_rot) {
            img().resize(Zsize(vol()), Ysize(vol()), Xsize(vol()));
        } else {
            img().resize(Ysize(vol()), Xsize(vol()));
        }
        transformer.setReal(img());
        transformer.getFourierAlias(F2D);

        // Set up the projector
        int data_dim = (do_3d_rot) ? 3 : 2;
        Projector projector((int) Xsize(vol()), interpolator, padding_factor, r_min_nn, data_dim);
        projector.computeFourierTransformMap(vol(), dummy, 2* r_max);

        if (do_only_one) {
            Euler_rotation3DMatrix(rot, tilt, psi, A3D);
            F2D.initZeros();
            projector.get2DFourierTransform(F2D, A3D);
            if (abs(xoff) > 0.001 || abs(yoff) > 0.001 || do_3d_rot && abs(zoff) > 0.001) {
                Matrix1D<RFLOAT> shift(2);
                XX(shift) = -xoff;
                YY(shift) = -yoff;
                if (do_3d_rot) {
                    shift.resize(3);
                    ZZ(shift) = -zoff;
                    shiftImageInFourierTransform(F2D, F2D, Xsize(vol()), XX(shift), YY(shift), ZZ(shift));
                } else {
                    shiftImageInFourierTransform(F2D, F2D, Xsize(vol()), XX(shift), YY(shift));
                }
            }

            // 1 Feb 2017 - Shaoda, add white noise to 2D / 3D single images
            if (do_add_noise) {
                if (!(stddev_white_noise > 0.0) || fn_model != "")
                    REPORT_ERROR("ERROR: Only add --white_noise to a single image!");
                // fftw normalization and factor sqrt(2) for two-dimensionality of complex plane
                // TODO: sqrt(2) ??? Why ???
                stddev_white_noise /= data_dim == 3 ? Xsize(vol()) * Xsize(vol()) : Xsize(vol()) * sqrt(2);
                // Add white noise
                FOR_ALL_ELEMENTS_IN_FFTW_TRANSFORM(F2D) {
                    DIRECT_A3D_ELEM(F2D, k, i, j).real += rnd_gaus(0.0, stddev_white_noise);
                    DIRECT_A3D_ELEM(F2D, k, i, j).imag += rnd_gaus(0.0, stddev_white_noise);
                }
            }

            transformer.inverseFourierTransform(F2D, img());
            // Shift the image back to the center...
            CenterFFT(img(), false);
            img.setSamplingRateInHeader(angpix);
            img.write(fn_out);
            std::cout<<" Done writing "<<fn_out<<std::endl;
        } else {  // !do_only_one
            init_progress_bar(MDang.numberOfObjects());
            DFo.clear();
            rot = tilt = psi = xoff = yoff = zoff = 0.;

            // Can only add noise to multiple images
            // Feb 01,2017 - Shaoda, now we can add white noise to 2D / 3D single images
            if (do_add_noise) {
                if (fn_model != "") {
                    model.read(fn_model);
                } else if (stddev_white_noise > 0.) {
                    stddev_white_noise /= Xsize(vol()) * sqrt(2); // fftw normalization and factor sqrt(2) for two-dimensionality of complex plane
                } else {
                    REPORT_ERROR("ERROR: When adding noise provide either --model_noise or --white_noise");
                }
            }

            long int imgno = 0;
            long int max_imgno = MDang.numberOfObjects() - 1;
            FOR_ALL_OBJECTS_IN_METADATA_TABLE(MDang) {
                rot  = MDang.getValue<RFLOAT>(EMDL::ORIENT_ROT);
                tilt = MDang.getValue<RFLOAT>(EMDL::ORIENT_TILT);
                psi  = MDang.getValue<RFLOAT>(EMDL::ORIENT_PSI);
                xoff = MDang.getValue<RFLOAT>(EMDL::ORIENT_ORIGIN_X_ANGSTROM);
                yoff = MDang.getValue<RFLOAT>(EMDL::ORIENT_ORIGIN_Y_ANGSTROM);
                if (do_3d_rot)
                zoff = MDang.getValue<RFLOAT>(EMDL::ORIENT_ORIGIN_Z_ANGSTROM);

                xoff /= angpix;
                yoff /= angpix;
                zoff /= angpix;

                Euler_rotation3DMatrix(rot, tilt, psi, A3D);
                F2D.initZeros();
                projector.get2DFourierTransform(F2D, A3D);

                if (abs(xoff) > 0.001 || abs(yoff) > 0.001 || do_3d_rot && abs(zoff) > 0.001) {
                    Matrix1D<RFLOAT> shift(2);
                    XX(shift) = -xoff;
                    YY(shift) = -yoff;

                    if (do_3d_rot) {
                        shift.resize(3);
                        ZZ(shift) = -zoff;
                        shiftImageInFourierTransform(F2D, F2D, Xsize(vol()), XX(shift), YY(shift), ZZ(shift) );
                    } else {
                        shiftImageInFourierTransform(F2D, F2D, Xsize(vol()), XX(shift), YY(shift));
                    }
                }

                // Apply CTF if necessary
                if (do_ctf || do_ctf2) {
                    if (do_3d_rot) {
                        Image<RFLOAT> Ictf;
                        FileName fn_ctf = MDang.getValue<FileName>(EMDL::CTF_IMAGE);
                        Ictf.read(fn_ctf);

                        // Set the CTF-image in Fctf
                        Fctf.resize(F2D);

                        // If there is a redundant half, get rid of it
                        if (Xsize(Ictf()) == Ysize(Ictf())) {
                            Ictf().setXmippOrigin();
                            // Set the CTF-image in Fctf
                            FOR_ALL_ELEMENTS_IN_FFTW_TRANSFORM(Fctf) {
                                // Use negative kp,ip and jp indices, because the origin in the ctf_img lies half a pixel to the right of the actual center....
                                DIRECT_A3D_ELEM(Fctf, k, i, j) = A3D_ELEM(Ictf(), -kp, -ip, -jp);
                            }
                        } else if (Xsize(Ictf()) == Ysize(Ictf()) / 2 + 1) {
                            // otherwise, just window the CTF to the current resolution
                            windowFourierTransform(Ictf(), Fctf, Ysize(Fctf));
                        } else {
                            // if dimensions are neither cubical nor FFTW, stop
                            REPORT_ERROR("3D CTF volume must be either cubical or adhere to FFTW format!");
                        }
                    } else {
                        CTF ctf = CTF(MDang, &obsModel); // This MDimg only contains one particle!
                        Fctf.resize(F2D);
                        ctf.getFftwImage(Fctf, Xsize(vol()), Xsize(vol()), angpix, ctf_phase_flipped, false,  do_ctf_intact_1st_peak, true);
                    }

                    FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(F2D) {
                        DIRECT_MULTIDIM_ELEM(F2D, n) *= DIRECT_MULTIDIM_ELEM(Fctf, n);
                        if (do_ctf2)
                        DIRECT_MULTIDIM_ELEM(F2D, n) *= DIRECT_MULTIDIM_ELEM(Fctf, n);
                    }
                }

                // Apply Gaussian noise
                if (do_add_noise) {
                    if (fn_model != "") {
                        //// 23 May 2014: for preparation of 1.3 release: removed reading a exp_model, replaced by just reading MDang
                        // This does however mean that I no longer know mic_id of this image: replace by 0....
                        FileName fn_group;
                        if (MDang.containsLabel(EMDL::MLMODEL_GROUP_NAME)) {
                            fn_group = MDang.getValue<FileName>(EMDL::MLMODEL_GROUP_NAME);
                        } else {
                            if (MDang.containsLabel(EMDL::MICROGRAPH_NAME)) {
                                FileName fn_orig, fn_pre, fn_jobnr;
                                fn_orig = MDang.getValue<FileName>(EMDL::MICROGRAPH_NAME);
                                if (!decomposePipelineFileName(fn_orig, fn_pre, fn_jobnr, fn_group)) {
                                    fn_group = fn_orig; // Not a pipeline filename; use as is
                                }
                            } else {
                                REPORT_ERROR("ERROR: cannot find rlnGroupName or rlnMicrographName in the input --ang file...");
                            }
                        }
                        int my_mic_id = -1;
                        for (int mic_id = 0; mic_id < model.group_names.size(); mic_id++) {
                            if (fn_group == model.group_names[mic_id]) {
                                my_mic_id = mic_id;
                                break;
                            }
                        }
                        if (my_mic_id < 0)
                            REPORT_ERROR("ERROR: cannot find " + fn_group + " in the input model file...");

                        RFLOAT normcorr = 1.0;
                        if (MDang.containsLabel(EMDL::IMAGE_NORM_CORRECTION)) {
                            normcorr = MDang.getValue<RFLOAT>(EMDL::IMAGE_NORM_CORRECTION);
                        }

                        // Add coloured noise
                        FOR_ALL_ELEMENTS_IN_FFTW_TRANSFORM(F2D) {
                            int Nyquist = model.ori_size / 2;
                            int ires = std::min((int) round(sqrt((RFLOAT) (kp * kp + ip * ip + jp * jp))), Nyquist);   // at freqs higher than Nyquist: use last sigma2 value

                            RFLOAT sigma = sqrt(DIRECT_A1D_ELEM(model.sigma2_noise[my_mic_id], ires));
                            DIRECT_A3D_ELEM(F2D, k, i, j).real += rnd_gaus(0.0, sigma);
                            DIRECT_A3D_ELEM(F2D, k, i, j).imag += rnd_gaus(0.0, sigma);
                        }
                    } else {
                        // Add white noise
                        FOR_ALL_ELEMENTS_IN_FFTW_TRANSFORM(F2D) {
                            DIRECT_A3D_ELEM(F2D, k, i, j).real += rnd_gaus(0.0, stddev_white_noise);
                            DIRECT_A3D_ELEM(F2D, k, i, j).imag += rnd_gaus(0.0, stddev_white_noise);
                        }
                    }
                }

                img().initZeros();
                transformer.inverseFourierTransform(F2D, img());
                // Shift the image back to the center...
                CenterFFT(img(), false);

                // Subtract the projection from the corresponding experimental image
                if (do_subtract_exp || do_simulate) {
                    fn_expimg = MDang.getValue<FileName>(EMDL::IMAGE_NAME);
                    MDang.setValue(EMDL::IMAGE_ORI_NAME, fn_expimg); // Store fn_expimg in rlnOriginalParticleName
                    expimg.read(fn_expimg);
                    img() = expimg() - img();
                }

                // If we're simulating realistic images, then now add CTF-affected projection again
                if (do_simulate) {
                    // Take random orientation from the input STAR file is fn_ang_simulate is empty. Otherwise, use fn_ang_simulate
                    if (fn_ang_simulate == "") {
                        long int random_imgno = -1;
                        while (random_imgno < 0 || random_imgno > max_imgno) {
                            random_imgno = rnd_unif() * max_imgno;
                        }

                        rot  = MDang.getValue<RFLOAT>(EMDL::ORIENT_ROT,               random_imgno);
                        tilt = MDang.getValue<RFLOAT>(EMDL::ORIENT_TILT,              random_imgno);
                        psi  = MDang.getValue<RFLOAT>(EMDL::ORIENT_PSI,               random_imgno);
                        xoff = MDang.getValue<RFLOAT>(EMDL::ORIENT_ORIGIN_X_ANGSTROM, random_imgno);
                        yoff = MDang.getValue<RFLOAT>(EMDL::ORIENT_ORIGIN_Y_ANGSTROM, random_imgno);
                        if (do_3d_rot)
                        zoff = MDang.getValue<RFLOAT>(EMDL::ORIENT_ORIGIN_Z_ANGSTROM, random_imgno);

                        xoff /= angpix;
                        yoff /= angpix;
                        zoff /= angpix;
                    } else {
                        rot  = MDang_sim.getValue<RFLOAT>(EMDL::ORIENT_ROT,               imgno);
                        tilt = MDang_sim.getValue<RFLOAT>(EMDL::ORIENT_TILT,              imgno);
                        psi  = MDang_sim.getValue<RFLOAT>(EMDL::ORIENT_PSI,               imgno);
                        xoff = MDang_sim.getValue<RFLOAT>(EMDL::ORIENT_ORIGIN_X_ANGSTROM, imgno);
                        yoff = MDang_sim.getValue<RFLOAT>(EMDL::ORIENT_ORIGIN_Y_ANGSTROM, imgno);
                        if (do_3d_rot)
                        zoff = MDang_sim.getValue<RFLOAT>(EMDL::ORIENT_ORIGIN_Z_ANGSTROM, imgno);

                        xoff /= angpix;
                        yoff /= angpix;
                        zoff /= angpix;
                    }

                    Euler_rotation3DMatrix(rot, tilt, psi, A3D);
                    F2D.initZeros();
                    projector.get2DFourierTransform(F2D, A3D);

                    if (abs(xoff) > 0.001 || abs(yoff) > 0.001 || do_3d_rot && abs(zoff) > 0.001) {
                        Matrix1D<RFLOAT> shift(2);
                        XX(shift) = -xoff;
                        YY(shift) = -yoff;

                        if (do_3d_rot) {
                            shift.resize(3);
                            ZZ(shift) = -zoff;
                            shiftImageInFourierTransform(F2D, F2D, Xsize(vol()), XX(shift), YY(shift), ZZ(shift) );
                        } else {
                            shiftImageInFourierTransform(F2D, F2D, Xsize(vol()), XX(shift), YY(shift));
                        }
                    }

                    // Apply CTF
                    if (do_ctf || do_ctf2) {
                        if (do_3d_rot) {
                            Image<RFLOAT> Ictf;
                            FileName fn_ctf = MDang.getValue<FileName>(EMDL::CTF_IMAGE);
                            Ictf.read(fn_ctf);
                            Ictf().setXmippOrigin();

                            // If there is a redundant half, get rid of it
                            if (Xsize(Ictf()) == Ysize(Ictf())) {
                                Ictf().setXmippOrigin();
                                FOR_ALL_ELEMENTS_IN_FFTW_TRANSFORM(Fctf) {
                                    // Use negative kp,ip and jp indices, because the origin in the ctf_img lies half a pixel to the right of the actual center....
                                    DIRECT_A3D_ELEM(Fctf, k, i, j) = A3D_ELEM(Ictf(), -kp, -ip, -jp);
                                }
                            } else if (Xsize(Ictf()) == Ysize(Ictf()) / 2 + 1) {
                                // otherwise, just window the CTF to the current resolution
                                windowFourierTransform(Ictf(), Fctf, Ysize(Fctf));
                            } else {
                                // if dimensions are neither cubical nor FFTW, stop
                                REPORT_ERROR("3D CTF volume must be either cubical or adhere to FFTW format!");
                            }
                        } else {
                            CTF ctf = CTF(MDang, MDang, imgno);  // Repetition of MDang is redundant
                            Fctf.resize(F2D);
                            ctf.getFftwImage(Fctf, Xsize(vol()), Xsize(vol()), angpix, ctf_phase_flipped, false,  do_ctf_intact_1st_peak, true);
                        }

                        FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(F2D) {
                            DIRECT_MULTIDIM_ELEM(F2D, n) *= DIRECT_MULTIDIM_ELEM(Fctf, n);
                            if (do_ctf2)
                            DIRECT_MULTIDIM_ELEM(F2D, n) *= DIRECT_MULTIDIM_ELEM(Fctf, n);
                        }
                    }

                    expimg().initZeros();
                    transformer.inverseFourierTransform(F2D, expimg());
                    // Shift the image back to the center...
                    CenterFFT(expimg(), false);

                    // Modify the strength of the signal
                    if (fabs(simulate_SNR - 1.0) > 0.000001) {
                        expimg() *= simulate_SNR;
                    }

                    img() += expimg();

                }

                img.setSamplingRateInHeader(angpix);
                if (do_3d_rot) {
                    fn_img.compose(fn_out, imgno + 1,"mrc");
                    img.write(fn_img);
                } else {
                    // Write this particle to the stack on disc
                    // First particle: write stack in overwrite mode, from then on just append to it
                    fn_img.compose(imgno+1,fn_out+".mrcs");
                    if (imgno == 0) {
                        img.write(fn_img, -1, false, WRITE_OVERWRITE);
                    } else {
                        img.write(fn_img, -1, false, WRITE_APPEND);
                    }
                }

                // Set the image name to the output STAR file
                DFo.addObject();
                DFo.setObject(MDang.getObject());
                DFo.setValue(EMDL::IMAGE_NAME,fn_img);

                if (do_simulate) {
                    DFo.setValue(EMDL::ORIENT_ROT, rot);
                    DFo.setValue(EMDL::ORIENT_TILT, tilt);
                    DFo.setValue(EMDL::ORIENT_PSI, psi);
                    DFo.setValue(EMDL::ORIENT_ORIGIN_X_ANGSTROM, xoff * angpix);
                    DFo.setValue(EMDL::ORIENT_ORIGIN_Y_ANGSTROM, yoff * angpix);
                    if (do_3d_rot)
                        DFo.setValue(EMDL::ORIENT_ORIGIN_Z_ANGSTROM, zoff * angpix);
                }

                if (imgno % 60 == 0) progress_bar(imgno);
                imgno++;
            }
            progress_bar(MDang.numberOfObjects());

            // Write out STAR file with all information
            fn_img = fn_out + ".star";
            obsModel.save(DFo, fn_img);
            std::cout << " Done writing " << imgno << " images in " << fn_img << std::endl;

        }
    }
};


int main(int argc, char *argv[]) {
    time_config();
    project_parameters prm;

    try {
        prm.read(argc, argv);
        prm.project();
    } catch (RelionError XE) {
        //prm.usage();
        std::cerr << XE;
        return RELION_EXIT_FAILURE;
    }

    return RELION_EXIT_SUCCESS;
}
