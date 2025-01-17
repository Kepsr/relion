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

#include "src/particle_subtractor.h"
#include "src/jaz/ctf_helper.h"


// Like modulo, but return the modulus instead of zero.
static inline int modulo_alt(int x, int y) {
    return (x - 1) % y + 1;
}

void ParticleSubtractor::read(int argc, char **argv) {
    parser.setCommandLine(argc, argv);
    int gen_section = parser.addSection("General options");
    fn_opt = parser.getOption("--i", "Name of optimiser.star file from refinement/classification to use for subtraction", "");
    fn_out = parser.getOption("--o", "Output directory name", "Subtract/");
    fn_msk = parser.getOption("--mask", "Name of the 3D mask with all density that should be kept, i.e. not subtracted", "");
    fn_sel = parser.getOption("--data", "Name of particle STAR file, in case not all particles from optimiser are to be used", "");
    ignore_class = parser.checkOption("--ignore_class", "Ignore the rlnClassNumber column in the particle STAR file.");
    fn_revert = parser.getOption("--revert", "Name of particle STAR file to revert. When this is provided, all other options are ignored.", "");
    do_ssnr = parser.checkOption("--ssnr", "Don't subtract, only calculate average spectral SNR in the images");

    int center_section = parser.addSection("Centering options");
    do_recenter_on_mask = parser.checkOption("--recenter_on_mask", "Use this flag to center the subtracted particles on projections of the centre-of-mass of the input mask");
    new_center.resize(3);
    XX(new_center) = textToInteger(parser.getOption("--center_x", "X-coordinate of 3D coordinate, which will be projected to center the subtracted particles.", "9999"));
    YY(new_center) = textToInteger(parser.getOption("--center_y", "Y-coordinate of 3D coordinate, which will be projected to center the subtracted particles.", "9999"));
    ZZ(new_center) = textToInteger(parser.getOption("--center_z", "Z-coordinate of 3D coordinate, which will be projected to center the subtracted particles.", "9999"));
    boxsize = textToInteger(parser.getOption("--new_box", "Output size of the subtracted particles", "-1"));

    verb = 1;
    // Check for errors in the command-line option
    if (parser.checkForErrors())
        REPORT_ERROR("Errors encountered on the command line (see above), exiting...");

    if (fn_opt.empty() == fn_revert.empty())
        REPORT_ERROR("Please specify only one of --i OR --revert");
}

void ParticleSubtractor::usage() {
    parser.writeUsage(std::cout);
}

void ParticleSubtractor::divideLabour(int _rank, int _size, long int &my_first, long int &my_last) {
    if (opt.do_split_random_halves) {
        if (size < 2) REPORT_ERROR("You need to run this program with at least 2 MPI processes for subtraction of 3D auto-refine jobs");
        opt.mydata.divideParticlesInRandomHalves(0, opt.do_helical_refine);

        int my_halfset = modulo_alt(_rank, 2);
        int mysize = my_halfset == 1 ? _size / 2 : _size / 2 + _size % 2;
        divide_equally(opt.mydata.numberOfParticles(my_halfset), mysize, _rank / 2, my_first, my_last);
        if (my_halfset == 2) {
            my_first += opt.mydata.numberOfParticles(1);
            my_last  += opt.mydata.numberOfParticles(1);
        }
    } else {
        divide_equally(opt.mydata.numberOfParticles(), _size, _rank, my_first, my_last);
    }
    //std::cerr << " rank= " << _rank << "size= " << _size << " my_first= " << my_first << " my_last= " << my_last << std::endl;
}

// Initialise some stuff after reading
void ParticleSubtractor::initialise(int _rank, int _size) {
    rank = _rank;
    size = _size;

    if (rank > 0) verb = 0;

    if (!fn_revert.empty())
        return;

    // Make directory for output particles
    if (fn_out[fn_out.length() - 1] != '/') fn_out += "/";
    if (verb > 0) {
        int res = system(("mkdir -p " + fn_out + "Particles").c_str());
    }

    opt.read(fn_opt, rank, true); // true means: prevent prereading all particle images
    nr_particles_in_optics_group.resize(opt.mydata.obsModel.opticsMdt.size(), 0);

    // Overwrite the particles STAR file with a smaller subset
    if (!fn_sel.empty()) {
        opt.mydata.clear();
        bool is_helical_segment = opt.do_helical_refine || opt.mymodel.ref_dim == 2 && opt.helical_tube_outer_diameter > 0.0;
        opt.mydata.read(fn_sel, false, false, false, is_helical_segment);
    }

    divideLabour(rank, size, my_first_part_id, my_last_part_id);

    Image<RFLOAT> Imask;
    if (!fn_msk.empty() && !do_ssnr) {
        if (verb > 0) std::cout << " + Reading in mask ... " << std::endl;
        // Mask stuff
        Imask.read(fn_msk);
        Imask().setXmippOrigin();

        const auto range = minmax(Imask());
        if (range.first < 0.0 || range.second > 1.0) {
            REPORT_ERROR("ERROR: the keep_inside mask has values outside the range [0,1]");
        }
    } else {
        Imask().reshape(opt.mymodel.Iref[0]);
        Imask().initZeros();
    }

    if (do_ssnr) {
        sum_S2.   initZeros(opt.mymodel.ori_size / 2 + 1);
        sum_N2.   initZeros(opt.mymodel.ori_size / 2 + 1);
        sum_count.initZeros(opt.mymodel.ori_size / 2 + 1);
    }

    if (do_recenter_on_mask) {
        Imask().centerOfMass(new_center);
        do_center = true;
    } else if (
        (int) XX(new_center) != 9999 && 
        (int) YY(new_center) != 9999 && 
        (int) ZZ(new_center) != 9999
    ) {
        do_center = true;
    } else {
        std::fill(new_center.begin(), new_center.end(), 0);
        do_center = false;
    }

    if (verb > 0 && (do_center || opt.fn_body_masks != "None")) {
        std::cout << " + The subtracted particles will be re-centred on projections of 3D-coordinate: ("
                  << XX(new_center) << " , " << YY(new_center) << " , " << ZZ(new_center) << ")" << std::endl;
    }

    if (opt.fn_body_masks != "None") {
        if (verb > 0) {
            std::cout << " + Initialising multi-body subtraction ..." << std::endl;
        }

        // This creates a rotation matrix for (rot,tilt,psi) = (0,90,0)
        // It will be used to make all Abody orientation matrices relative to (0,90,0) instead of the more logical (0,0,0)
        // This is useful, as psi-priors are ill-defined around tilt=0, as rot becomes the same as -psi!!
        A_rot90 = rotation3DMatrix(-90.0, 'Y', false);

        // Find out which body has the biggest overlap with the keepmask, use these orientations
        RFLOAT best_overlap = 0.0;
        subtract_body = -1;
        for (int ibody = 0; ibody < opt.mymodel.nr_bodies; ibody++) {
            auto &body = opt.mymodel.masks_bodies[ibody];
            if (!Imask().sameShape(body)) {
                Imask().printShape();
                body.printShape();
                REPORT_ERROR("ERROR: input mask is not of same shape as body masks.");
            }

            RFLOAT overlap = 0.0;
            for (long int n = 0; n < Imask().size(); n++)
                overlap += body[n] * Imask()[n];

            if (overlap > best_overlap) {
                best_overlap = overlap;
                subtract_body = ibody;
            }
        }

        if (subtract_body < 0) REPORT_ERROR("ERROR: input mask does not overlap with any of the bodies....");

        // Apply the inverse of the keepmask to all the mask_bodies
        for (int obody = 0; obody < opt.mymodel.nr_bodies; obody++) {
            int ii = direct::elem(opt.mymodel.pointer_body_overlap, subtract_body, obody);
            auto &body = opt.mymodel.masks_bodies[ii];
            for (long int n = 0; n < Imask().size(); n++) {
                body[n] *= 1.0 - Imask()[n];  // Expression templates
            }
        }
    } else {
        // For normal refinement/classification: just apply the inverse of keepmask to the references
        for (int iclass = 0; iclass < opt.mymodel.nr_classes; iclass++) {
            auto &class_ = opt.mymodel.Iref[iclass];
            if (!Imask().sameShape(class_)) {
                Imask().printShape();
                class_.printShape();
                REPORT_ERROR("ERROR: input mask is not of same shape as reference inside the optimiser.");
            }

            for (long int n = 0; n < Imask().size(); n++) {
                class_[n] *= 1.0 - Imask()[n];  // Expression templates
            }
        }
    }

    if (verb > 0) {
        std::cout << " + Calculating Fourier transforms of the maps ..." << std::endl;
    }

    // Now set up the Projectors inside the model
    opt.mymodel.setFourierTransformMaps(false); // false means ignore tau2_class

    // ensure even boxsize of subtracted images
    if (boxsize > 0) { boxsize -= boxsize % 2; }
}

void ParticleSubtractor::revert() {
    ObservationModel obsModel;
    MetaDataTable MD;

    ObservationModel::loadSafely(fn_revert, obsModel, MD);

    if (!MD.containsLabel(EMDL::IMAGE_ORI_NAME))
        REPORT_ERROR("The input STAR file does not contain the rlnImageOriginalName column.");

    if (!MD.containsLabel(EMDL::IMAGE_NAME))
        REPORT_ERROR("The input STAR file does not contain the rlnImageName column");

    // Swap image names
    for (long int i : MD) {
        const FileName f1 = MD.getValue<std::string>(EMDL::IMAGE_ORI_NAME, i);
        const FileName f2 = MD.getValue<std::string>(EMDL::IMAGE_NAME,     i);
        MD.setValue(EMDL::IMAGE_ORI_NAME, f2, i);
        MD.setValue(EMDL::IMAGE_NAME,     f1, i);
    }

    // Fix box size
    std::vector<bool> fixed_box_size(obsModel.numberOfOpticsGroups(), false);
    for (long int i : MD) {
        const int og = obsModel.getOpticsGroup(MD);
        if (fixed_box_size[og]) continue;

        FileName img_name = MD.getValue<std::string>(EMDL::IMAGE_NAME, i);
        FileName fn_img;
        long int dummy;
        img_name.decompose(dummy, fn_img);

        if (!exists(fn_img))
            REPORT_ERROR("Failed to read " + fn_img + " to determine the box size.");
        Image<RFLOAT> Ihead;
        Ihead.read(img_name, false, -1, nullptr, true);
        if (Xsize(Ihead()) != Ysize(Ihead()))
            REPORT_ERROR("Particle " + img_name + " is not square.");
        obsModel.setBoxSize(og, Xsize(Ihead()));
        obsModel.opticsMdt.setValue(EMDL::IMAGE_SIZE, Xsize(Ihead()), og);

        fixed_box_size[og] = true;
    }

    for (int i = 0; i < obsModel.numberOfOpticsGroups(); i++) {
        if (!fixed_box_size[i]) {
            std::cerr << "WARNING: could not determine the box size of optics group " << obsModel.getGroupName(i) << std::endl;
        } else {
            std::cout << "The box size of the optics group " << obsModel.getGroupName(i) << " was updated to " << obsModel.getBoxSize(i) << " px" << std::endl;
        }
    }

    obsModel.save(MD, fn_out + "original.star");
    std::cout << "Writen " << (fn_out + "original.star") << std::endl;
}

void ParticleSubtractor::run() {

    long int nr_parts = my_last_part_id - my_first_part_id + 1;
    long int barstep = std::max(1l, nr_parts / 120);
    if (verb > 0) {
        std::cout << (do_ssnr ?
            " + Calculating SNR for all particles ..." :
            " + Subtracting all particles ..."
        ) << std::endl;
        time_config();
        init_progress_bar(nr_parts);
    }

    MDimg_out.clear();
    //for (long int part_id_sorted = my_first_part_id, cc = 0; part_id_sorted <= my_last_part_id; part_id_sorted++, cc++)
    for (long int part_id_sorted = my_first_part_id, cc = 0; part_id_sorted <= my_last_part_id; part_id_sorted++, cc++) {

        //long int part_id = opt.mydata.sorted_idx[part_id_sorted];
        if (cc % barstep == 0) {
            if (pipeline_control_check_abort_job())
                exit(RELION_EXIT_ABORTED);
        }

        long int part_id = opt.mydata.sorted_idx[part_id_sorted];
        subtractOneParticle(part_id, 0, cc);

        if (cc % barstep == 0 && verb > 0) progress_bar(cc);
    }

    if (verb > 0) progress_bar(nr_parts);
}

void ParticleSubtractor::saveStarFile(int myrank) {

    if (do_ssnr) {

        // Only leader writes out the STAR file
        if (myrank == 0) {
            MetaDataTable MD;
            for (int ires = 0; ires < Xsize(sum_S2); ires++) {
                MD.addObject();
                MD.setValue(EMDL::SPECTRAL_IDX, ires, ires);
                MD.setValue(EMDL::RESOLUTION, opt.mymodel.getResolution(ires), ires);
                MD.setValue(EMDL::RESOLUTION_ANGSTROM, opt.mymodel.getResolutionAngstrom(ires), ires);
                MD.setValue(EMDL::MLMODEL_SSNR_REF, sum_S2.elem(ires) / sum_N2.elem(ires), ires);
                MD.setValue(EMDL::MLMODEL_TAU2_REF, sum_S2.elem(ires) / sum_count.elem(ires) , ires);
                MD.setValue(EMDL::MLMODEL_SIGMA2_NOISE, sum_N2.elem(ires) / sum_count.elem(ires) , ires);
            }
            std::cout << " Writing out STAR file with spectral SNR in: " << fn_out << "spectral_snr.star" << " ..." << std::endl;
            MD.write(fn_out+"spectral_snr.star");
        }

    } else {

        // Reset image size in optics table, if the images were rewindowed in a different box
        if (boxsize > 0) {
            for (long int i : opt.mydata.obsModel.opticsMdt) {
                opt.mydata.obsModel.opticsMdt.setValue(EMDL::IMAGE_SIZE, boxsize, i);
            }
        }

        // Remove origin prior columns if present, as we have re-centered.
        if (do_center || opt.fn_body_masks != "None") {
            if (MDimg_out.containsLabel(EMDL::ORIENT_ORIGIN_X_PRIOR_ANGSTROM)) {
                MDimg_out.deactivateLabel(EMDL::ORIENT_ORIGIN_X_PRIOR_ANGSTROM);
            }
            if (MDimg_out.containsLabel(EMDL::ORIENT_ORIGIN_Y_PRIOR_ANGSTROM)) {
                MDimg_out.deactivateLabel(EMDL::ORIENT_ORIGIN_Y_PRIOR_ANGSTROM);
            }
        }

        FileName fn_star;
        if (size == 0) {
            fn_star = fn_out + "particles_subtracted.star";
            MDimg_out.deactivateLabel(EMDL::IMAGE_ID);
        } else {
            fn_star = fn_out + "Particles/subtracted_rank" + integerToString(myrank) + "star";
        }
        opt.mydata.obsModel.save(MDimg_out, fn_star);

    }

    #ifdef DEBUG
    std::cout << "myrank = " << myrank << " size = " << size << " my_first = " << my_first << " my_last = " << my_last << " num_items = " << MD.size() << " writing to " << fn_star << std::endl;
    #endif
}

void ParticleSubtractor::combineStarFile(int myrank) {

    if (do_ssnr) return;

    if (myrank != 0) REPORT_ERROR("BUG: this function should only be called by leader!");

    MetaDataTable MD;
    for (int i = 1; i < size; i++) {
        FileName fn_star = fn_out + "Particles/subtracted_rank" + integerToString(i) + "star";
        MD.read(fn_star, "particles");
        MDimg_out.append(MD);
    }

    MDimg_out.sort(EMDL::IMAGE_ID);
    MDimg_out.deactivateLabel(EMDL::IMAGE_ID);
    opt.mydata.obsModel.save(MDimg_out, fn_out + "particles_subtracted.star");

    std::cout << " + Saved STAR file with " << MDimg_out.size()
              << " subtracted particles in " << fn_out << "particles_subtracted.star" << std::endl;

}

FileName ParticleSubtractor::getParticleName(
    long int imgno, int myrank, int optics_group
) {
    if (imgno_to_filename.find(imgno) != imgno_to_filename.end())
        return imgno_to_filename[imgno];

    if (optics_group == -1) {
        std::cerr << "rank = " << rank << " imgno = " << imgno << std::endl;
        REPORT_ERROR("Logic error: optics group must be specified to register a new entry");
    }

    nr_particles_in_optics_group[optics_group]++;

    // Now write out the image
    FileName fn_stack = fn_out + "Particles/subtracted";
    if (size > 1)
        fn_stack += "_rank" + integerToString(myrank + 1);

    const auto fn_img = opt.mymodel.data_dim == 3 ?
        FileName::compose(fn_stack, imgno + 1, "mrc") :
        FileName::compose(nr_particles_in_optics_group[optics_group], fn_stack + "_opticsgroup" + integerToString(optics_group + 1) + ".mrcs");

    imgno_to_filename[imgno] = fn_img;
    #ifdef DEBUG
    std::cout << "rank = " << rank << " imgno = " << imgno << " fn_img = " << fn_img << std::endl;
    #endif
    return fn_img;
}

void ParticleSubtractor::subtractOneParticle(
    long int part_id, long int imgno, long int counter
) {
    // Read the particle image
    long int ori_img_id = opt.mydata.particles[part_id].images[imgno].id;
    int optics_group = opt.mydata.getOpticsGroup(part_id, 0);
    auto img = Image<RFLOAT>::from_filename(opt.mydata.particles[part_id].images[0].name);
    img().setXmippOrigin();

    // Make sure gold-standard is adhered to!
    int my_subset = modulo_alt(rank, 2);
    if (opt.do_split_random_halves && my_subset != opt.mydata.getRandomSubset(part_id)) {
        std::cerr << " rank= " << rank << " part_id= " << part_id << " opt.mydata.getRandomSubset(part_id)= " << opt.mydata.getRandomSubset(part_id) << std::endl;
        REPORT_ERROR("BUG:: gold-standard separation of halves is broken!");
    }

    // Get the consensus class, orientational parameters and norm (if present)
    RFLOAT my_pixel_size = opt.mydata.getImagePixelSize(part_id, 0);
    RFLOAT remap_image_sizes = (opt.mymodel.ori_size * opt.mymodel.pixel_size) / (Xsize(img()) * my_pixel_size);
    Vector<RFLOAT> my_old_offset (3), my_residual_offset (3), centering_offset (3);
    Matrix<RFLOAT> Aori;
    RFLOAT xoff, yoff, zoff, mynorm, scale;
    int myclass = 0;
    if (!ignore_class && opt.mydata.MDimg.containsLabel(EMDL::PARTICLE_CLASS)) {
        myclass = opt.mydata.MDimg.getValue<int>(EMDL::PARTICLE_CLASS, ori_img_id);
        if (myclass > opt.mymodel.nr_classes) {
            std::cerr << "A particle belongs to class " << myclass << " while the number of classes in the optimiser.star is only " << opt.mymodel.nr_classes << "." << std::endl;
            REPORT_ERROR("Tried to subtract a non-existing class from a particle. If you have performed non-alignment Class3D after Refine3D and want to subtract a map from the Refine3D job, use the --ignore_class option.");
        }
        myclass--; // Count from zero instead of one
    }
    RFLOAT rot        = opt.mydata.MDimg.getValue<RFLOAT>(EMDL::ORIENT_ROT,               ori_img_id);
    RFLOAT tilt       = opt.mydata.MDimg.getValue<RFLOAT>(EMDL::ORIENT_TILT,              ori_img_id);
    RFLOAT psi        = opt.mydata.MDimg.getValue<RFLOAT>(EMDL::ORIENT_PSI,               ori_img_id);
    XX(my_old_offset) = opt.mydata.MDimg.getValue<RFLOAT>(EMDL::ORIENT_ORIGIN_X_ANGSTROM, ori_img_id);
    YY(my_old_offset) = opt.mydata.MDimg.getValue<RFLOAT>(EMDL::ORIENT_ORIGIN_Y_ANGSTROM, ori_img_id);
    if (opt.mymodel.data_dim == 3)
    ZZ(my_old_offset) = opt.mydata.MDimg.getValue<RFLOAT>(EMDL::ORIENT_ORIGIN_Z_ANGSTROM, ori_img_id);
    // As of v3.1, offsets are in Angstrom: convert back to pixels!
    my_old_offset /= my_pixel_size;

    // Apply the norm_correction term
    try {
        mynorm = opt.mydata.MDimg.getValue<RFLOAT>(EMDL::IMAGE_NORM_CORRECTION, ori_img_id);
    } catch (const char *errmsg) {
        mynorm = 1.0;
    }

    if (opt.do_norm_correction) 
        img() *= opt.mymodel.avg_norm_correction / mynorm;

    Vector<RFLOAT> my_projected_com (3), my_refined_ibody_offset (3);
    if (opt.fn_body_masks != "None") {
        // 17May2017: Shift image to the projected COM for this body!
        // Aori is the original transformation matrix of the consensus refinement
        Aori = Euler::angles2matrix(rot, tilt, psi);
        my_projected_com = matmul(Aori, opt.mymodel.com_bodies[subtract_body]);

        // Subtract the projected COM offset, to position this body in the center
        my_old_offset -= my_projected_com;
        my_residual_offset = my_old_offset;
        // Apply the old_offset (rounded to avoid interpolation errors)
        for (auto& x: my_old_offset) { x = round(x); }
        img() = translate(img(), my_old_offset, WRAP);
        // keep track of the differences between the rounded and the original offsets
        my_residual_offset -= my_old_offset;
    }

    // Now that the particle is centered (for multibody), get the FourierTransform of the particle
    FourierTransformer transformer;
    MultidimArray<Complex> Fimg = transformer.FourierTransform(img());
    CenterFFTbySign(Fimg);
    auto Fctf = MultidimArray<RFLOAT>::ones(Fimg.xdim, Fimg.ydim, Fimg.zdim, Fimg.ndim);

    if (opt.do_ctf_correction) {
        if (opt.mymodel.data_dim == 3) {
            FileName fn_ctf = opt.mydata.MDimg.getValue<std::string>(EMDL::CTF_IMAGE, ori_img_id);
            auto Ictf = Image<RFLOAT>::from_filename(fn_ctf);

            if (Xsize(Ictf()) == Ysize(Ictf())) {
                // If there is a redundant half, get rid of it
                Ictf().setXmippOrigin();
                FOR_ALL_ELEMENTS_IN_FFTW_TRANSFORM(Fctf) {
                    // Use negative indices, because the origin in the ctf_img lies half a pixel to the right of the actual center....
                    direct::elem(Fctf, i, j, k) = Ictf().elem(-ip, -jp, -kp);
                }
            } else if (Xsize(Ictf()) == Ysize(Ictf()) / 2 + 1) {
                // Otherwise, just window the CTF to the current resolution
                Fctf = windowFourierTransform(Ictf(), Ysize(Fctf));
            } else {
                // if dimensions are neither cubical nor FFTW, stop
                REPORT_ERROR("3D CTF volume must be either cubical or adhere to FFTW format!");
            }
        } else {
            CTF ctf = CtfHelper::makeCTF(opt.mydata.MDimg, &opt.mydata.obsModel, ori_img_id);
            Fctf = CtfHelper::getFftwImage(
                ctf,
                Xsize(Fctf), Ysize(Fctf), Xsize(img()), Ysize(img()), my_pixel_size,
                &opt.mydata.obsModel,
                opt.ctf_phase_flipped, false, opt.intact_ctf_first_peak, true
            );
        }
    }

    MultidimArray<Complex> Fsubtrahend = MultidimArray<Complex>::zeros(Fimg);

    if (opt.fn_body_masks != "None") {
        // For multi-body refinement
        Matrix<RFLOAT> Aresi_subtract;
        for (int obody = 0; obody < opt.mymodel.nr_bodies; obody++) {
            // Unlike getFourierTransformsAndCtfs, no check for ibody==obody: also subtract rest of subtract_body!

            Vector<RFLOAT> body_offset(3);
            RFLOAT body_rot  = opt.mydata.MDbodies[obody].getValue<RFLOAT>(EMDL::ORIENT_ROT,               ori_img_id);
            RFLOAT body_tilt = opt.mydata.MDbodies[obody].getValue<RFLOAT>(EMDL::ORIENT_TILT,              ori_img_id);
            RFLOAT body_psi  = opt.mydata.MDbodies[obody].getValue<RFLOAT>(EMDL::ORIENT_PSI,               ori_img_id);
            XX(body_offset)  = opt.mydata.MDbodies[obody].getValue<RFLOAT>(EMDL::ORIENT_ORIGIN_X_ANGSTROM, ori_img_id);
            YY(body_offset)  = opt.mydata.MDbodies[obody].getValue<RFLOAT>(EMDL::ORIENT_ORIGIN_Y_ANGSTROM, ori_img_id);
            if (opt.mymodel.data_dim == 3)
            ZZ(body_offset)  = opt.mydata.MDbodies[obody].getValue<RFLOAT>(EMDL::ORIENT_ORIGIN_Z_ANGSTROM, ori_img_id);

            // As of v3.1, offsets are in Angstrom: convert back to pixels!
            body_offset /= my_pixel_size;

            // Aresi is the residual orientation for this obody
            Matrix<RFLOAT> Aresi = Euler::angles2matrix(body_rot, body_tilt, body_psi);
            if (obody == subtract_body) { Aresi_subtract = Aresi; }
            // The real orientation to be applied is the obody transformation applied and the original one
            Matrix<RFLOAT> Abody = Aori
                .matmul(opt.mymodel.orient_bodies[obody].transpose())
                .matmul(A_rot90)
                .matmul(Aresi)
                .matmul(opt.mymodel.orient_bodies[obody]);
            // Apply anisotropic mag and scaling
            if (opt.mydata.obsModel.hasMagMatrices)
                Abody = Abody.matmul(opt.mydata.obsModel.anisoMag(optics_group));
            Abody *= opt.mydata.obsModel.scaleDifference(optics_group, opt.mymodel.ori_size, opt.mymodel.pixel_size);

            // The following line gets the correct pointer to account for overlap in the bodies
            int oobody = direct::elem(opt.mymodel.pointer_body_overlap, subtract_body, obody);
            // Get the FT of the projection in the right direction
            auto FTo = opt.mymodel.PPref[oobody].get2DFourierTransform(
                Fimg.xdim, Fimg.ydim, Fimg.zdim, Abody);

            // Body is centered at its own COM: move it back to its place in the original particle image

            // Projected COM for this body (using Aori, just like above for ibody and my_projected_com!!!)
            Vector<RFLOAT> other_projected_com = matmul(Aori, opt.mymodel.com_bodies[obody]);

            // Subtract refined obody-displacement
            other_projected_com -= body_offset;

            // Subtract the projected COM already applied to this image for ibody
            other_projected_com -= my_projected_com;

            shiftImageInFourierTransform(
                FTo, (RFLOAT) Xsize(img()),
                XX(other_projected_com), YY(other_projected_com), ZZ(other_projected_com)
            );

            // Sum the Fourier transforms of all the obodies
            Fsubtrahend += FTo;
        }

        // Set orientations back into the original RELION system of coordinates

        // Write out the rot, tilt, psi as the combination of Aori and Aresi!! So get rid of the rotations around the tilt=90 axes,
        Matrix<RFLOAT> Abody = Aori
            .matmul(opt.mymodel.orient_bodies[subtract_body].transpose())
            .matmul(A_rot90)
            .matmul(Aresi_subtract)
            .matmul(opt.mymodel.orient_bodies[subtract_body]);
        angles_t angles = Euler::matrix2angles(Abody);

        // Store the optimal orientations in the MDimg table
        opt.mydata.MDimg.setValue(EMDL::ORIENT_ROT,  angles.rot,  ori_img_id);
        opt.mydata.MDimg.setValue(EMDL::ORIENT_TILT, angles.tilt, ori_img_id);
        opt.mydata.MDimg.setValue(EMDL::ORIENT_PSI,  angles.psi,  ori_img_id);

        // Also get refined offset for this body
        XX(my_refined_ibody_offset) = opt.mydata.MDbodies[subtract_body].getValue<RFLOAT>(EMDL::ORIENT_ORIGIN_X_ANGSTROM, ori_img_id);
        YY(my_refined_ibody_offset) = opt.mydata.MDbodies[subtract_body].getValue<RFLOAT>(EMDL::ORIENT_ORIGIN_Y_ANGSTROM, ori_img_id);
        if (opt.mymodel.data_dim == 3)
        ZZ(my_refined_ibody_offset) = opt.mydata.MDbodies[subtract_body].getValue<RFLOAT>(EMDL::ORIENT_ORIGIN_Z_ANGSTROM, ori_img_id);
        // As of v3.1, offsets are in Angstrom: convert back to pixels!
        my_refined_ibody_offset /= my_pixel_size;

        // re-center to new_center
        my_residual_offset += my_refined_ibody_offset + matmul(Abody, opt.mymodel.com_bodies[subtract_body] - new_center * opt.mymodel.pixel_size / my_pixel_size);
    } else {
        // Normal 3D classification/refinement: get the projection in rot,tilt,psi for the corresponding class
        auto A3D_pure_rot = Euler::angles2matrix(rot, tilt, psi);

        // Apply anisotropic mag and scaling
        auto A3D = opt.mydata.obsModel.hasMagMatrices ?
            A3D_pure_rot :
            A3D_pure_rot.matmul(opt.mydata.obsModel.anisoMag(optics_group));
        A3D *= opt.mydata.obsModel.scaleDifference(optics_group, opt.mymodel.ori_size, opt.mymodel.pixel_size);
        Fsubtrahend = opt.mymodel.PPref[myclass].get2DFourierTransform(
            Fimg.xdim, Fimg.ydim, Fimg.zdim, A3D);

        // Shift in opposite direction as offsets in the STAR file
        shiftImageInFourierTransform(
            Fsubtrahend, (RFLOAT) Xsize(img()),
            -XX(my_old_offset), -YY(my_old_offset), -ZZ(my_old_offset)
        );

        if (do_center)
            // Re-center the output particle to a new centre...
            my_residual_offset = my_old_offset - matmul(A3D_pure_rot, new_center) * opt.mymodel.pixel_size / my_pixel_size;
    }

    // Apply the CTF to the subtrahend projection
    if (opt.do_ctf_correction) {
        if (opt.mydata.obsModel.getCtfPremultiplied(optics_group)) {
            for (long int n = 0; n < Fsubtrahend.size(); n++) {
                Fsubtrahend[n] *= Fctf[n] * Fctf[n];  // Expression templates
            }
        } else {
            Fsubtrahend *= Fctf;
        }

        // Also do phase modulation, for beam tilt correction and other asymmetric aberrations
        opt.mydata.obsModel.modulatePhase(optics_group, Fsubtrahend);
        opt.mydata.obsModel.multiplyByMtf(optics_group, Fsubtrahend);
    }

    if (opt.do_scale_correction) {
        int group_id = opt.mydata.getGroupId(part_id, 0);
        RFLOAT scale = opt.mymodel.scale_correction[group_id];
        Fsubtrahend *= scale;
    }

    // Do the actual subtraction
    Fimg -= Fsubtrahend;

    if (do_ssnr) {
        // Don't write out subtracted image,
        // only accumulate power of the signal (in Fsubtrahend) divided by the power of the noise (now in Fimg)
        FOR_ALL_ELEMENTS_IN_FFTW_TRANSFORM(Fimg) {
            long int idx = round(hypot((double) ip, jp, kp));
            int idx_remapped = round(remap_image_sizes * idx);
            if (idx_remapped < opt.mymodel.ori_size / 2 + 1) {
                RFLOAT S2 = norm(direct::elem(Fsubtrahend, i, j, k));
                RFLOAT N2 = norm(direct::elem(Fimg,        i, j, k));
                // division by two keeps the numbers similar to tau2 and sigma2_noise,
                // which are per real/imaginary component
                sum_S2.elem(idx_remapped) += S2 / 2.0;
                sum_N2.elem(idx_remapped) += N2 / 2.0;
                sum_count.elem(idx_remapped) += 1.0;
            }
        }
    } else {
        // And go finally back to real-space
        CenterFFTbySign(Fimg);
        img() = transformer.inverseFourierTransform(Fimg);

        if (do_center || opt.fn_body_masks != "None") {
            // Recenter the particles
            centering_offset = my_residual_offset;
            for (auto &x : centering_offset) { x = round(x); }
            my_residual_offset -= centering_offset;
            img() = translate(img(), centering_offset, WRAP);

            // Set the non-integer difference between the rounded centering offset and the actual offsets in the STAR file
            opt.mydata.MDimg.setValue(EMDL::ORIENT_ORIGIN_X_ANGSTROM, my_pixel_size * XX(my_residual_offset), ori_img_id);
            opt.mydata.MDimg.setValue(EMDL::ORIENT_ORIGIN_Y_ANGSTROM, my_pixel_size * YY(my_residual_offset), ori_img_id);
            if (opt.mymodel.data_dim == 3) {
                opt.mydata.MDimg.setValue(EMDL::ORIENT_ORIGIN_Z_ANGSTROM, my_pixel_size * ZZ(my_residual_offset), ori_img_id);
            }
        }

        // Rebox the image
        if (boxsize > 0) {
            int dimensionality = img().getDim();
            if (dimensionality == 2) {
                img() = img().windowed(
                    Xmipp::init(boxsize), Xmipp::last(boxsize),
                    Xmipp::init(boxsize), Xmipp::last(boxsize)
                );
            } else if (dimensionality == 3) {
                img() = img().windowed(
                    Xmipp::init(boxsize), Xmipp::last(boxsize),
                    Xmipp::init(boxsize), Xmipp::last(boxsize),
                    Xmipp::init(boxsize), Xmipp::last(boxsize)
                );
            }
        }

        // Now write out the image & set filenames in output metadatatable
        FileName fn_img = getParticleName(counter, rank, optics_group);
        opt.mydata.MDimg.setValue(EMDL::IMAGE_NAME, fn_img, ori_img_id);
        opt.mydata.MDimg.setValue(EMDL::IMAGE_ORI_NAME, opt.mydata.particles[part_id].images[0].name, ori_img_id);
        //Also set the original order in the input STAR file for later combination
        opt.mydata.MDimg.setValue(EMDL::IMAGE_ID, ori_img_id, ori_img_id);
        const long int i = MDimg_out.addObject();
        MDimg_out.setObject(opt.mydata.MDimg.getObject(ori_img_id), i);

        //printf("Writing: fn_orig = %s counter = %ld rank = %d optics_group = %d fn_img = %s SIZE = %d nr_particles_in_optics_group[optics_group] = %d\n", fn_orig.c_str(), counter, rank, optics_group+1, fn_img.c_str(), Xsize(img()), nr_particles_in_optics_group[optics_group]);
        img.setSamplingRateInHeader(my_pixel_size);
        if (opt.mymodel.data_dim == 3) {
            img.write(fn_img);
        } else {
            img.write(
                fn_img, -1, false, bool(nr_particles_in_optics_group[optics_group]) ?
                    WRITE_APPEND : WRITE_OVERWRITE
                    // If there are particles, enter append mode.
                    // Otherwise (if there are no particles), enter (over)write mode.
            );
        }
    }
}
