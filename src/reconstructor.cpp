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
#include "src/reconstructor.h"
#include "src/jaz/ctf_helper.h"

void Reconstructor::read(int argc, char **argv) {
    parser.setCommandLine(argc, argv);

    int general_section = parser.addSection("General options");
    fn_sel = parser.getOption("--i", "Input STAR file with the projection images and their orientations", "");
    fn_out = parser.getOption("--o", "Name for output reconstruction","relion.mrc");
    fn_sym = parser.getOption("--sym", "Symmetry group", "c1");
    maxres = textToFloat(parser.getOption("--maxres", "Maximum resolution (in Angstrom) to consider in Fourier space (default Nyquist)", "-1"));
    padding_factor = textToFloat(parser.getOption("--pad", "Padding factor", "2"));
    image_path = parser.getOption("--img", "Optional: image path prefix", "");
    subset = textToInteger(parser.getOption("--subset", "Subset of images to consider (1: only reconstruct half1; 2: only half2; other: reconstruct all)", "-1"));
    chosen_class = textToInteger(parser.getOption("--class", "Consider only this class (-1: use all classes)", "-1"));
    angpix  = textToFloat(parser.getOption("--angpix", "Pixel size in the reconstruction (take from first optics group by default)", "-1"));

    int ctf_section = parser.addSection("CTF options");
    do_ctf = parser.checkOption("--ctf", "Apply CTF correction");
    intact_ctf_first_peak = parser.checkOption("--ctf_intact_first_peak", "Leave CTFs intact until first peak");
    ctf_phase_flipped = parser.checkOption("--ctf_phase_flipped", "Images have been phase flipped");
    only_flip_phases = parser.checkOption("--only_flip_phases", "Do not correct CTF-amplitudes, only flip phases");

    int ewald_section = parser.addSection("Ewald-sphere correction options");
    do_ewald = parser.checkOption("--ewald", "Correct for Ewald-sphere curvature (developmental)");
    mask_diameter  = textToFloat(parser.getOption("--mask_diameter", "Diameter (in A) of mask for Ewald-sphere curvature correction", "-1."));
    width_mask_edge = textToInteger(parser.getOption("--width_mask_edge", "Width (in pixels) of the soft edge on the mask", "3"));
    is_reverse = parser.checkOption("--reverse_curvature", "Try curvature the other way around");
    newbox = textToInteger(parser.getOption("--newbox", "Box size of reconstruction after Ewald sphere correction", "-1"));
    nr_sectors = textToInteger(parser.getOption("--sectors", "Number of sectors for Ewald sphere correction", "2"));
    skip_mask = parser.checkOption("--skip_mask", "Do not apply real space mask during Ewald sphere correction");
    skip_weighting = parser.checkOption("--skip_weighting", "Do not apply weighting during Ewald sphere correction");

    if (verb > 0 && do_ewald && mask_diameter < 0 && !(skip_mask && skip_weighting))
        REPORT_ERROR("To apply Ewald sphere correction (--ewald), you have to specify the mask diameter(--mask_diameter).");

    int helical_section = parser.addSection("Helical options");
    nr_helical_asu = textToInteger(parser.getOption("--nr_helical_asu", "Number of helical asymmetrical units", "1"));
    helical_rise = textToFloat(parser.getOption("--helical_rise", "Helical rise (in Angstroms)", "0."));
    helical_twist = textToFloat(parser.getOption("--helical_twist", "Helical twist (in degrees, + for right-handedness)", "0."));

    int expert_section = parser.addSection("Expert options");
    fn_sub = parser.getOption("--subtract","Subtract projections of this map from the images used for reconstruction", "");
    interpolator = parser.checkOption("--NN", "Use nearest-neighbour instead of linear interpolation before gridding correction") ?
        NEAREST_NEIGHBOUR : TRILINEAR;
    blob_radius = textToFloat(parser.getOption("--blob_r", "Radius of blob for gridding interpolation", "1.9"));
    blob_order = textToInteger(parser.getOption("--blob_m", "Order of blob for gridding interpolation", "0"));
    blob_alpha = textToFloat(parser.getOption("--blob_a", "Alpha-value of blob for gridding interpolation", "15"));
    iter = textToInteger(parser.getOption("--iter", "Number of gridding-correction iterations", "10"));
    ref_dim = textToInteger(parser.getOption("--refdim", "Dimension of the reconstruction (2D or 3D)", "3"));
    angular_error = textToFloat(parser.getOption("--angular_error", "Apply random deviations with this standard deviation (in degrees) to each of the 3 Euler angles", "0."));
    shift_error = textToFloat(parser.getOption("--shift_error", "Apply random deviations with this standard deviation (in Angstrom) to each of the 2 translations", "0."));
    do_fom_weighting = parser.checkOption("--fom_weighting", "Weight particles according to their figure-of-merit (_rlnParticleFigureOfMerit)");
    fn_fsc = parser.getOption("--fsc", "FSC-curve for regularized reconstruction", "");
    do_3d_rot = parser.checkOption("--3d_rot", "Perform 3D rotations instead of backprojections from 2D images");
    ctf_dim  = textToInteger(parser.getOption("--reconstruct_ctf", "Perform a 3D reconstruction from 2D CTF-images, with the given size in pixels", "-1"));
    do_reconstruct_ctf2 = parser.checkOption("--ctf2", "Reconstruct CTF^2 and then take the sqrt of that");
    skip_gridding = parser.checkOption("--skip_gridding", "Skip gridding part of the reconstruction");
    fn_debug = parser.getOption("--debug", "Rootname for debug reconstruction files", "");
    debug_ori_size =  textToInteger(parser.getOption("--debug_ori_size", "Rootname for debug reconstruction files", "1"));
    debug_size =  textToInteger(parser.getOption("--debug_size", "Rootname for debug reconstruction files", "1"));
    fn_noise = parser.getOption("--reconstruct_noise","Reconstruct noise using sigma2 values in this model STAR file", "");
    read_weights = parser.checkOption("--read_weights", "Developmental: read freq. weight files");
    do_debug = parser.checkOption("--write_debug_output", "Write out arrays with data and weight terms prior to reconstruct");
    do_external_reconstruct = parser.checkOption("--external_reconstruct", "Write out BP denominator and numerator for external_reconstruct program");
    verb = textToInteger(parser.getOption("--verb", "Verbosity", "1"));

    // Hidden
    r_min_nn = textToInteger(getParameter(argc, argv, "--r_min_nn", "10"));

    // Check for errors in the command-line option
    if (parser.checkForErrors())
        REPORT_ERROR("Errors encountered on the command line (see above), exiting...");
}

void Reconstructor::usage() {
    parser.writeUsage(std::cout);
}

void Reconstructor::initialise() {
    do_reconstruct_ctf = ctf_dim > 0;
    if (do_reconstruct_ctf) {
        do_ctf = false;
        padding_factor = 1.0;
    }

    do_ignore_optics = false;
    // Read MetaData file, which should have the image names and their angles!
    if (fn_debug.empty()) {
        ObservationModel::loadSafely(fn_sel, obsModel, DF, "particles", 0, false);
        if (obsModel.opticsMdt.empty()) {
            do_ignore_optics = true;
            DF.read(fn_sel);
        }
    }

    if (verb > 0 && (subset == 1 || subset == 2) && !DF.containsLabel(EMDL::PARTICLE_RANDOM_SUBSET)) {
        REPORT_ERROR("The rlnRandomSubset column is missing in the input STAR file.");
    }

    if (verb > 0 && (chosen_class >= 0) && !DF.containsLabel(EMDL::PARTICLE_CLASS)) {
        REPORT_ERROR("The rlnClassNumber column is missing in the input STAR file.");
    }

    randomize_random_generator();

    if (do_ewald) do_ctf = true;

    // Is this 2D or 3D data?
    data_dim = 2; // Initial default value

    if (!fn_noise.empty())
        model.read(fn_noise);

    // Get dimension of the images
    if (do_reconstruct_ctf) {
        output_boxsize = ctf_dim;
    } else {
        fn_img = DF.getValue<std::string>(EMDL::IMAGE_NAME, DF.size() - 1);

        if (!image_path.empty()) {
            fn_img = image_path + "/" + fn_img.substr(fn_img.find_last_of("/") + 1);
        }

        auto img0 = Image<RFLOAT>::from_filename(fn_img, false);
        output_boxsize = Xsize(img0());
        // When doing Ewald-curvature correction or when having optics groups: allow reconstructing smaller box than the input images (which should have large boxes!!)
        if ((do_ewald || !do_ignore_optics) && newbox > 0) {
            output_boxsize = newbox;
        }

        if (do_3d_rot) {
            data_dim = 3;
        } else {
            // If not specifically provided, we autodetect it
            if (do_ignore_optics) {
                data_dim = img0().getDim();
                std::cout << " + Taking data dimensions from the first image: " << data_dim << std::endl;
            } else {
                data_dim = obsModel.opticsMdt.getValue<int>(EMDL::IMAGE_DIMENSIONALITY, 0);
                std::cout << " + Taking data dimensions from the first optics group: " << data_dim << std::endl;
            }
        }
    }

    if (angpix < 0.0) {
        if (do_ignore_optics) {
            if (
                DF.containsLabel(EMDL::CTF_MAGNIFICATION) && DF.containsLabel(EMDL::CTF_DETECTOR_PIXEL_SIZE)
            ) {
                const long int i = DF.size() - 1;
                RFLOAT mag   = DF.getValue<RFLOAT>(EMDL::CTF_MAGNIFICATION, i);
                RFLOAT dstep = DF.getValue<RFLOAT>(EMDL::CTF_DETECTOR_PIXEL_SIZE, i);
                angpix = 10000.0 * dstep / mag;
                if (verb > 0)
                    std::cout << " + Using pixel size calculated from magnification and detector pixel size in the input STAR file: " << angpix << std::endl;
            } else {
                REPORT_ERROR("ERROR: cannot find pixel size in input STAR file, provide it using --angpix");
            }
        } else {
            angpix = obsModel.getPixelSize(0);
            std::cout << " + Taking angpix from the first optics group: " << angpix << std::endl;
        }
    }
    r_max = maxres < 0.0 ? -1 : ceil(output_boxsize * angpix / maxres);
}

void Reconstructor::run() {
    if (!fn_debug.empty()) {
        readDebugArrays();
    } else {
        initialise();
        backproject();
    }

    reconstruct();
}

void Reconstructor::readDebugArrays() {
    if (verb > 0)
        std::cout << " + Reading in the debug arrays ... " << std::endl;

    // We first read the image to set the data_dim automatically from backprojector data
    auto Ireal = Image<RFLOAT>::from_filename(fn_debug + "_data_real.mrc");
    data_dim = Ireal().getDim();

    backprojector = BackProjector(debug_ori_size, 3, fn_sym, interpolator, padding_factor, r_min_nn, blob_order, blob_radius, blob_alpha, data_dim, skip_gridding);

    backprojector.initialiseDataAndWeight(debug_size);
    if (verb > 0) {
        std::cout << " Size of data array: " ;
        backprojector.data.printShape();
        std::cout << " Size of weight array: " ;
        backprojector.weight.printShape();
    }

    Ireal().setXmippOrigin();
    Ireal().xinit = 0;

    if (verb > 0) {
        std::cout << " Size of reconstruction: " ;
        Ireal().printShape();
    }

    auto Iimag = Image<RFLOAT>::from_filename(fn_debug + "_data_imag.mrc");
    Iimag().setXmippOrigin();
    Iimag().xinit = 0;
    FOR_ALL_ELEMENTS_IN_ARRAY3D(backprojector.data, i, j, k) {
        backprojector.data.elem(i, j, k).real = Ireal().elem(i, j, k);
        backprojector.data.elem(i, j, k).imag = Iimag().elem(i, j, k);
    }

    auto It = Image<RFLOAT>::from_filename(fn_debug + "_weight.mrc");
    It().setXmippOrigin();
    It().xinit = 0;
    FOR_ALL_ELEMENTS_IN_ARRAY3D(It(), i, j, k) {
        backprojector.weight.elem(i, j, k) = It().elem(i, j, k);
    }
    output_boxsize = debug_ori_size;
}

void Reconstructor::backproject(int rank, int size) {
    if (!fn_sub.empty()) {
        projector = Projector(output_boxsize, interpolator, padding_factor, r_min_nn);
        auto sub = Image<RFLOAT>::from_filename(fn_sub);
        MultidimArray<RFLOAT> dummy;
        projector.computeFourierTransformMap(sub(), dummy, 2 * r_max);
    }

    backprojector = BackProjector(
        output_boxsize, ref_dim, fn_sym, interpolator,
        padding_factor, r_min_nn, blob_order,
        blob_radius, blob_alpha, data_dim, skip_gridding
    );
    backprojector.initZeros(2 * r_max);

    long int nr_parts = DF.size();
    long int barstep = std::max(1l, nr_parts / (size * 120));
    if (verb > 0) {
        std::cout << " + Back-projecting all images ..." << std::endl;
        time_config();
        init_progress_bar(nr_parts);
    }

    for (long int ipart = 0; ipart < nr_parts; ipart++) {
        if (ipart % size == rank)
            backprojectOneParticle(ipart);

        if (ipart % barstep == 0 && verb > 0)
            progress_bar(ipart);
    }

    if (verb > 0)
        progress_bar(nr_parts);
}

void Reconstructor::backprojectOneParticle(long int p) {
    RFLOAT rot, tilt, psi, r_ewald_sphere;
    FourierTransformer transformer;

    int randSubset = 0, classid = 0;
    randSubset = DF.getValue<int>(EMDL::PARTICLE_RANDOM_SUBSET, p);
    classid    = DF.getValue<int>(EMDL::PARTICLE_CLASS,         p);

    if (subset >= 1 && subset <= 2 && randSubset != subset)
        return;

    if (chosen_class >= 0 && chosen_class != classid)
        return;

    // Rotations
    if (ref_dim == 2) {
        rot = tilt = 0.0;
    } else {
        rot  = DF.getValue<RFLOAT>(EMDL::ORIENT_ROT,  p);
        tilt = DF.getValue<RFLOAT>(EMDL::ORIENT_TILT, p);
    }

    psi = DF.containsLabel(EMDL::ORIENT_PSI) ? DF.getValue<RFLOAT>(EMDL::ORIENT_PSI, p) : 0.0;

    if (angular_error > 0.0) {
        rot  += rnd_gaus(0.0, angular_error);
        tilt += rnd_gaus(0.0, angular_error);
        psi  += rnd_gaus(0.0, angular_error);
    }

    Matrix2D<RFLOAT> A3D = Euler::angles2matrix(rot, tilt, psi);

    // If we are considering Ewald sphere curvature, the mag. matrix
    // has to be provided to the backprojector explicitly
    // (to avoid creating an Ewald ellipsoid)
    int opticsGroup = -1;
    int myBoxSize = output_boxsize; // Without optics groups, the output box size is always the same as the one from the input images
    RFLOAT myPixelSize = angpix; // Without optics groups, the pixel size is always the same as the one from the input images
    bool ctf_premultiplied = false;
    if (!do_ignore_optics) {
        opticsGroup       = obsModel.getOpticsGroup(DF, p);
        myBoxSize         = obsModel.getBoxSize(opticsGroup);
        myPixelSize       = obsModel.getPixelSize(opticsGroup);
        ctf_premultiplied = obsModel.getCtfPremultiplied(opticsGroup);
        if (do_ewald && ctf_premultiplied)
            REPORT_ERROR("We cannot perform Ewald sphere correction on CTF premultiplied particles.");
        if (!do_ewald && obsModel.hasMagMatrices) {
            A3D *= obsModel.anisoMag(opticsGroup);
        }
        A3D *= obsModel.scaleDifference(opticsGroup, output_boxsize, angpix);
    }

    // Translations (either through phase-shifts or in real space
    Matrix1D<RFLOAT> trans = Matrix1D<RFLOAT>::zeros(data_dim);
    std::array<EMDL::EMDLabel, 3> origin_labels {
        EMDL::ORIENT_ORIGIN_X_ANGSTROM,
        EMDL::ORIENT_ORIGIN_Y_ANGSTROM,
        EMDL::ORIENT_ORIGIN_Z_ANGSTROM
    };
    for (int i = 0; i < data_dim; ++i) {
        trans[i] = DF.getValue<RFLOAT>(origin_labels[i], p);
    }

    if (shift_error > 0.0) {
        for (int i = 0; i < data_dim; ++i) { trans[i] += rnd_gaus(0.0, shift_error); }
    }

    // As of v3.1, shifts are in Angstroms in the STAR files, convert back to pixels here
    trans /= myPixelSize;

    RFLOAT fom;
    if (do_fom_weighting)
    fom = DF.getValue<RFLOAT>(EMDL::PARTICLE_FOM, p);

    // Use either translate OR shiftImageInFourierTransform!!
    // img() = translate(img(), trans, WRAP);

    MultidimArray<Complex> F2D;
    FileName fn_img;
    Image<RFLOAT> img;

    if (!do_reconstruct_ctf && fn_noise.empty()) {
        fn_img = DF.getValue<std::string>(EMDL::IMAGE_NAME, p);
        img.read(fn_img);
        img().setXmippOrigin();
        F2D = transformer.FourierTransform(img());
        CenterFFTbySign(F2D);

        if (abs(trans[0]) > 0.0 || abs(trans[1]) > 0.0 || abs(trans[2]) > 0.0) {
            // trans[2] is 0 in case data_dim=2
            shiftImageInFourierTransform(F2D, Xsize(img()), trans[0], trans[1], trans[2]);
        }
    } else {
        if (data_dim == 3) {
            F2D.resize(myBoxSize, myBoxSize, myBoxSize / 2 + 1);
        } else {
            F2D.resize(myBoxSize,            myBoxSize / 2 + 1);
        }
    }

    if (!fn_noise.empty()) {
        // TODO: Refactor code duplication from relion_project!
        FileName fn_group;
        const long int i = DF.size() - 1;
        if (DF.containsLabel(EMDL::MLMODEL_GROUP_NAME)) {
            fn_group = DF.getValue<std::string>(EMDL::MLMODEL_GROUP_NAME, i);
        } else if (DF.containsLabel(EMDL::MICROGRAPH_NAME)) {
            fn_group = DF.getValue<std::string>(EMDL::MICROGRAPH_NAME, i);
        } else {
            REPORT_ERROR("ERROR: cannot find rlnGroupName or rlnMicrographName in the input --i file...");
        }

        const auto it = std::find(model.group_names.begin(), model.group_names.end(), fn_group);
        if (it == model.group_names.end()) REPORT_ERROR("ERROR: cannot find " + fn_group + " in the input model file...");
        int imic = it - model.group_names.begin();

        RFLOAT normcorr = 1.0;
        if (DF.containsLabel(EMDL::IMAGE_NORM_CORRECTION))
        normcorr = DF.getValue<RFLOAT>(EMDL::IMAGE_NORM_CORRECTION, i);

        // Make coloured noise image
        FOR_ALL_ELEMENTS_IN_FFTW_TRANSFORM(F2D) {
            const int ires = std::min((int) round(euclid(ip, jp, kp)), myBoxSize / 2);
            // at freqs higher than Nyquist: use last sigma2 value
            const RFLOAT sigma = sqrt(direct::elem(model.sigma2_noise[imic], ires));
            direct::elem(F2D, i, j, k) += Complex(rnd_gaus(0.0, sigma), rnd_gaus(0.0, sigma));
        }
    }

    MultidimArray<RFLOAT> Fctf;
    Fctf.resize(F2D);
    Fctf = 1.0;
    MultidimArray<Complex> F2DP, F2DQ;

    // Apply CTF if necessary
    if (do_ctf || do_reconstruct_ctf) {

        // Also allow 3D CTF correction here
        if (data_dim == 3) {
            Image<RFLOAT> Ictf;
            FileName fn_ctf;
            try {
                fn_ctf = DF.getValue<std::string>(EMDL::CTF_IMAGE, p);
            } catch (const char *errmsg) {
                REPORT_ERROR("ERROR: cannot find rlnCtfImage for 3D CTF correction!");
            }
            Ictf.read(fn_ctf);

            // If there is a redundant half, get rid of it
            if (Xsize(Ictf()) == Ysize(Ictf())) {
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
            CTF ctf = do_ignore_optics ? CtfHelper::makeCTF(DF, p)
                                       : CtfHelper::makeCTF(DF, &obsModel, p);

            Fctf = CtfHelper::getFftwImage(
                ctf,
                Xsize(Fctf), Ysize(Fctf), myBoxSize, myBoxSize, myPixelSize,
                do_ignore_optics ? nullptr : &obsModel,
                ctf_phase_flipped, only_flip_phases,
                intact_ctf_first_peak, true
            );

            int opticsGroup;
            if (!do_ignore_optics) {
                opticsGroup = DF.getValue<int>(EMDL::IMAGE_OPTICS_GROUP, p) - 1;
                obsModel.demodulatePhase(opticsGroup, F2D, p);
                obsModel.divideByMtf    (opticsGroup, F2D, p);
            }

            // Ewald-sphere curvature correction
            if (do_ewald) {
                applyCTFPandCTFQ(
                    F2D,
                    ctf, do_ignore_optics ? nullptr : &obsModel, opticsGroup,
                    transformer, F2DP, F2DQ, skip_mask
                );

                if (!skip_weighting) {
                    // Also calculate W, store again in Fctf
                    CtfHelper::applyWeightEwaldSphereCurvature_noAniso(
                        ctf,
                        Fctf, myBoxSize, myBoxSize, myPixelSize,
                        do_ignore_optics ? nullptr : &obsModel,
                        opticsGroup,
                        mask_diameter
                    );
                }

                // Also calculate the radius of the Ewald sphere (in pixels)
                r_ewald_sphere = myBoxSize * myPixelSize / ctf.lambda;
            }
        }
    }

    // Subtract reference projection
    MultidimArray<Complex> Fsub;
    if (!fn_sub.empty()) {
        Fsub.resize(F2D);
        projector.get2DFourierTransform(Fsub, A3D);

        // Apply CTF if necessary
        if (do_ctf) { Fsub *= Fctf; }

        F2D -= Fsub;

        // Back-project difference image
        backprojector.set2DFourierTransform(F2D, A3D);
    } else {
        if (do_reconstruct_ctf) {
            for (long int n = 0; n < F2D.size(); n++) {
                F2D[n] = Fctf[n];
                if (do_reconstruct_ctf2) { F2D[n] *= Fctf[n]; }
                Fctf[n] = 1.0;
            }
        } else if (do_ewald) {
            Fctf *= Fctf;
        } else if (do_ctf) {
            // "Normal" reconstruction, multiply X by CTF, and W by CTF^2
            if (!ctf_premultiplied) { F2D *= Fctf; }
            Fctf *= Fctf;
        }

        // Do the following after squaring the CTFs!
        if (do_fom_weighting) {
            for (long int n = 0; n < F2D.size(); n++) {
                F2D[n]  *= fom;
                Fctf[n] *= fom;
            }
        }

        if (read_weights) {
            std::string fullName = DF.getValue<std::string>(EMDL::IMAGE_NAME, 0);
            std::string name = fullName.substr(fullName.find("@") + 1);

            if (!image_path.empty()) {
                name = image_path + "/" + name.substr(name.find_last_of("/") + 1);
            }

            std::string wghName = name;
            wghName = wghName.substr(0, wghName.find_last_of('.')) + "_weight.mrc";

            Image<RFLOAT> wgh;
            wgh.read(wghName);

            if (
                Fctf.ndim != wgh().ndim ||
                Fctf.zdim != wgh().zdim ||
                Fctf.ydim != wgh().ydim ||
                Fctf.xdim != wgh().xdim
            ) {
                REPORT_ERROR(wghName + " and " + name + " are of unequal size.\n");
            }

            for (long int n = 0; n < Fctf.ndim; n++)
            for (long int z = 0; z < Fctf.zdim; z++)
            for (long int y = 0; y < Fctf.ydim; y++)
            for (long int x = 0; x < Fctf.xdim; x++) {
                direct::elem(Fctf, x, y, z, n) *= direct::elem(wgh(), x, y, z, n);
            }
        }

        direct::elem(F2D, 0, 0) = 0.0;

        if (do_ewald) {
            Matrix2D<RFLOAT> magMat = !do_ignore_optics && obsModel.hasMagMatrices ?
                obsModel.getMagMatrix(opticsGroup) :
                Matrix2D<RFLOAT>::identity(2);

            backprojector.set2DFourierTransform(F2DP, A3D, &Fctf, r_ewald_sphere, +1.0, &magMat);
            backprojector.set2DFourierTransform(F2DQ, A3D, &Fctf, r_ewald_sphere, -1.0, &magMat);
        } else {
            backprojector.set2DFourierTransform(F2D, A3D, &Fctf);
        }
    }
}

void Reconstructor::reconstruct() {
    bool do_map = false;
    bool do_use_fsc = false;
    MultidimArray<RFLOAT> fsc;
    fsc.resize(output_boxsize / 2 + 1);

    if (!fn_fsc.empty()) {
        do_map = true;
        do_use_fsc = true;
        MetaDataTable MDfsc;
        MDfsc.read(fn_fsc);
        for (long int i : MDfsc) {
            int    idx = MDfsc.getValue<int>(EMDL::SPECTRAL_IDX, i);
            RFLOAT val = MDfsc.getValue<RFLOAT>(EMDL::MLMODEL_FSC_HALVES_REF, i);
            fsc(idx) = val;
        }
    }

    if (verb > 0)
        std::cout << " + Starting the reconstruction ..." << std::endl;

    backprojector.symmetrise(nr_helical_asu, helical_twist, helical_rise/angpix);

    Image<RFLOAT> vol;
    if (do_reconstruct_ctf) {

        vol().initZeros(ctf_dim, ctf_dim, ctf_dim);
        vol().setXmippOrigin();

        FOR_ALL_ELEMENTS_IN_ARRAY3D(vol(), i, j, k) {
            int jp = j;
            int ip = i;
            int kp = k;

            // for negative j's: use inverse
            if (j < 0) {
                jp = -j;
                ip = -i;
                kp = -k;
            }

            if (
                jp >= Xinit(backprojector.data) && jp <= Xlast(backprojector.data) &&
                ip >= Yinit(backprojector.data) && ip <= Ylast(backprojector.data) &&
                kp >= Zinit(backprojector.data) && kp <= Zlast(backprojector.data)
            ) {
                if (backprojector.weight.elem(ip, jp, kp) > 0.0) {
                    vol().elem(i, j, k) = backprojector.data.elem(ip, jp, kp) / backprojector.weight.elem(ip, jp, kp);
                    if (do_reconstruct_ctf2)
                        vol().elem(i, j, k) = sqrt(vol().elem(i, j, k));
                }
            }
        }
    } else {

        if (do_debug) {
            Image<RFLOAT> It;
            FileName fn_tmp = fn_out.withoutExtension();
            It().resize(backprojector.data);
            for (long int n = 0; n < It().size(); n++) {
                It()[n] = backprojector.data[n].real;
            }
            It.write(fn_tmp + "_data_real.mrc");
            for (long int n = 0; n < It().size(); n++) {
                It()[n] = backprojector.data[n].imag;
            }
            It.write(fn_tmp + "_data_imag.mrc");
            It() = backprojector.weight;
            It.write(fn_tmp + "_weight.mrc");
        }

        MultidimArray<RFLOAT> tau2, tmp;
        if (do_use_fsc) 
        backprojector.updateSSNRarrays(
            1.0, tau2, tmp, tmp, tmp, fsc, do_use_fsc, true
        );

        if (do_external_reconstruct) {
            FileName fn_root = fn_out.withoutExtension();
            vol() = backprojector.externalReconstruct(
                fn_root, tau2, tmp, tmp, tmp, false, 1.0, 1
            );
        } else {
            vol() = backprojector.reconstruct(iter, do_map, tau2);
        }
    }

    vol.setSamplingRateInHeader(angpix);
    vol.write(fn_out);
    if (verb > 0)
        std::cout << " + Done! Written output map in: " << fn_out << std::endl;
}

void Reconstructor::applyCTFPandCTFQ(
    MultidimArray<Complex> &Fin,
    CTF &ctf, ObservationModel *obsModel, int opticsGroup,
    FourierTransformer &transformer,
    MultidimArray<Complex> &outP, MultidimArray<Complex> &outQ, bool skip_mask
) {
    // FourierTransformer transformer;
    outP.resize(Fin);
    outQ.resize(Fin);
    float angle_step = 180.0 / nr_sectors;
    for (float angle = 0.0; angle < 180.0; angle += angle_step) {
        MultidimArray<Complex> CTFP(Fin), Fapp(Fin);
        MultidimArray<RFLOAT> Iapp(Ysize(Fin), Ysize(Fin));
        // Two passes: one for CTFP, one for CTFQ
        for (int ipass = 0; ipass < 2; ipass++) {
            bool is_my_positive = (ipass == 1) == is_reverse;

            // Get CTFP and multiply the Fapp with it
            CTFP = CtfHelper::getCTFPImage(
                ctf,
                Fin.xdim, Fin.ydim, Ysize(Fin), Ysize(Fin), angpix, 
                obsModel, opticsGroup,
                is_my_positive, angle
            );

            Fapp = Fin * CTFP; // element-wise complex multiplication!

            if (!skip_mask) {
                // inverse transform and mask out the particle....
                CenterFFTbySign(Fapp);
                Iapp = transformer.inverseFourierTransform(Fapp);

                softMaskOutsideMap(Iapp, round(mask_diameter / (angpix * 2.0)), (RFLOAT) width_mask_edge);

                // Re-box to a smaller size if necessary....
                if (newbox > 0 && newbox < Ysize(Fin)) {
                    Iapp.setXmippOrigin();
                    Iapp.window(
                        Xmipp::init(newbox), Xmipp::init(newbox),
                        Xmipp::last(newbox), Xmipp::last(newbox)
                    );

                }
                // Back into Fourier-space
                Fapp = transformer.FourierTransform(Iapp);  // std::move?
                CenterFFTbySign(Fapp);
            }

            // First time round: resize the output arrays
            if (ipass == 0 && fabs(angle) < Xmipp::epsilon) {
                outP.resize(Fapp);
                outQ.resize(Fapp);
            }

            // Now set back the right parts into outP (first pass) or outQ (second pass)
            float anglemin = angle + 90.0 - 0.5 * angle_step;
            float anglemax = angle + 90.0 + 0.5 * angle_step;

            // Angles larger than 180
            bool is_angle_reverse = anglemin >= 180.0;
            if (is_angle_reverse) {
                anglemin -= 180.0;
                anglemax -= 180.0;
            }

            bool PorQoutP = is_angle_reverse != (ipass == 0);
            auto &myCTFPorQ  = PorQoutP ? outP : outQ;
            auto &myCTFPorQb = PorQoutP ? outQ : outP;

            // Deal with sectors with the Y-axis in the middle of the sector...
            bool do_wrap_max = anglemin < 180.0 && anglemax > 180.0;
            if (do_wrap_max) { anglemax -= 180.0; }

            // Convert to radians
            anglemin = radians(anglemin);
            anglemax = radians(anglemax);
            FOR_ALL_ELEMENTS_IN_FFTW_TRANSFORM2D(CTFP) {
                RFLOAT theta = atan2(jp, ip);
                // Only take the relevant sector now...
                if (do_wrap_max) {
                    if (theta >= anglemin) {
                        direct::elem(myCTFPorQ, i, j) = direct::elem(Fapp, i, j);
                    } else if (theta < anglemax) {
                        direct::elem(myCTFPorQb, i, j) = direct::elem(Fapp, i, j);
                    }
                } else {
                    if (theta >= anglemin && theta < anglemax) {
                        direct::elem(myCTFPorQ, i, j) = direct::elem(Fapp, i, j);
                    }
                }
            }
        }
    }
}
