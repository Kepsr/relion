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
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * This complete copyright notice must be included in any revised version of the
 * source code. Additional authorship citations may be added, but existing
 * author citations must be preserved.
 ***************************************************************************/

#include "src/ml_model.h"

#ifdef MDL_TIMING
    Timer mdl_timer;
    int TIMING_MDL_1 = proj_timer.setNew("MDL_1");
#define TIMING_TOC(id) mdl_timer.toc(id)
#else
    #define TIMING_TIC(id)
    #define TIMING_TOC(id)
#endif


void MlModel::initialise(bool do_sgd) {

    // Auxiliary vector with relevant size in Fourier space
    const MultidimArray<RFLOAT> zeros = MultidimArray<RFLOAT>::zeros(ori_size / 2 + 1);

    // Now resize all relevant vectors
    Iref.resize(nr_classes * nr_bodies);
    masks_bodies           .resize(nr_bodies);
    com_bodies             .resize(nr_bodies);
    rotate_direction_bodies.resize(nr_bodies);
    orient_bodies          .resize(nr_bodies);
    sigma_tilt_bodies      .resize(nr_bodies, 0.0);
    sigma_psi_bodies       .resize(nr_bodies, 0.0);
    sigma_offset_bodies    .resize(nr_bodies, 0.0);
    keep_fixed_bodies      .resize(nr_bodies, 0);
    pointer_body_overlap   .resize(nr_bodies, nr_bodies);
    max_radius_mask_bodies .resize(nr_bodies, -1);
    pdf_class              .resize(nr_classes, 1.0 / (RFLOAT) nr_classes);
    pdf_direction          .resize(nr_classes * nr_bodies);
    group_names            .resize(nr_groups, "");
    sigma2_noise           .resize(nr_groups);
    nr_particles_per_group .resize(nr_groups);
    tau2_class             .resize(nr_classes * nr_bodies, zeros);
    fsc_halves_class       .resize(nr_classes * nr_bodies, zeros);
    sigma2_class           .resize(nr_classes * nr_bodies, zeros);
    data_vs_prior_class    .resize(nr_classes * nr_bodies, zeros);
    fourier_coverage_class .resize(nr_classes * nr_bodies, zeros);
    // TODO handle these two correctly.
    bfactor_correction.resize(nr_groups, 0.0);
    scale_correction  .resize(nr_groups, 1.0);

    acc_rot               .resize(nr_classes * nr_bodies, 0);
    acc_trans             .resize(nr_classes * nr_bodies, 0);
    estimated_resolution  .resize(nr_classes * nr_bodies, 0);
    total_fourier_coverage.resize(nr_classes * nr_bodies, 0);

    helical_twist.resize(nr_classes, 0);
    helical_rise .resize(nr_classes, 0);

    if (ref_dim == 2) {
        Matrix1D<RFLOAT> empty (2);
        prior_offset_class.resize(nr_classes * nr_bodies, empty);
    }
    // These arrays will be resized when they are filled
    orientability_contrib.resize(nr_classes * nr_bodies);

    Projector ref(ori_size, interpolator, padding_factor, r_min_nn, data_dim);
    PPref.clear();
    PPrefRank.clear();
    // Now fill the entire vector with instances of "ref"
    if (nr_classes != 1 && nr_bodies != 1)
        REPORT_ERROR("MlModel::initialise() - nr_bodies or nr_classes must be 1");
    PPref.resize(nr_classes * nr_bodies, ref);

    if (this->do_sgd = do_sgd)
        Igrad.resize(nr_classes);

}

// Reading from a file
void MlModel::read(FileName fn_in) {

    // Clear current model
    clear();

    // Open input file
    std::ifstream in(fn_in.data(), std::ios_base::in);
    if (in.fail())
        REPORT_ERROR((std::string) "MlModel::readStar: File " + fn_in + " cannot be read.");

    // Read general stuff
    MetaDataTable MDlog;
    MDlog.readStar(in, "model_general");
    const long int i = MDlog.size() - 1;
    try {
        ref_dim                  = MDlog.getValue<int>(EMDL::MLMODEL_DIMENSIONALITY, i);
        ori_size                 = MDlog.getValue<int>(EMDL::MLMODEL_ORIGINAL_SIZE, i);
        current_resolution       = MDlog.getValue<RFLOAT>(EMDL::MLMODEL_CURRENT_RESOLUTION, i);
        current_size             = MDlog.getValue<int>(EMDL::MLMODEL_CURRENT_SIZE, i);
        padding_factor           = MDlog.getValue<RFLOAT>(EMDL::MLMODEL_PADDING_FACTOR, i);
        interpolator             = MDlog.getValue<int>(EMDL::MLMODEL_INTERPOLATOR, i);
        r_min_nn                 = MDlog.getValue<int>(EMDL::MLMODEL_MINIMUM_RADIUS_NN_INTERPOLATION, i);
        pixel_size               = MDlog.getValue<RFLOAT>(EMDL::MLMODEL_PIXEL_SIZE, i);
        nr_classes               = MDlog.getValue<int>(EMDL::MLMODEL_NR_CLASSES, i);
        nr_groups                = MDlog.getValue<int>(EMDL::MLMODEL_NR_GROUPS, i);
        tau2_fudge_factor        = MDlog.getValue<RFLOAT>(EMDL::MLMODEL_TAU2_FUDGE_FACTOR, i);
        avg_norm_correction      = MDlog.getValue<RFLOAT>(EMDL::MLMODEL_NORM_CORRECTION_AVG, i);
        orientational_prior_mode = MDlog.getValue<int>(EMDL::MLMODEL_PRIOR_MODE, i);
        sigma2_rot               = MDlog.getValue<RFLOAT>(EMDL::MLMODEL_SIGMA_ROT, i);
        sigma2_tilt              = MDlog.getValue<RFLOAT>(EMDL::MLMODEL_SIGMA_TILT, i);
        sigma2_psi               = MDlog.getValue<RFLOAT>(EMDL::MLMODEL_SIGMA_PSI, i);
        LL                       = MDlog.getValue<RFLOAT>(EMDL::MLMODEL_LL, i);
        ave_Pmax                 = MDlog.getValue<RFLOAT>(EMDL::MLMODEL_AVE_PMAX, i);
    } catch (const char *errmsg) {
        REPORT_ERROR("MlModel::readStar: incorrect model_general table");
    }

    try {
        sigma2_offset = MDlog.getValue<RFLOAT>(EMDL::MLMODEL_SIGMA_OFFSET_ANGSTROM, i);
    } catch (const char *errmsg) { try {
        sigma2_offset = MDlog.getValue<RFLOAT>(EMDL::MLMODEL_SIGMA_OFFSET, i) * pixel_size;
    } catch (const char *errmsg) {
        REPORT_ERROR("MlModel::readStar: incorrect model_general table: cannot find sigma_offset");
    } }

    // Retain compability with model files written by Relion prior to 1.4

    #define TRY_MDLOG_GET(variable, T, label, defaultvalue) \
    try { \
        (variable) = MDlog.getValue<T>(label, i); \
    } catch (const char *errmsg) { \
        (variable) = (defaultvalue); \
    }

    TRY_MDLOG_GET(data_dim,  int,  EMDL::MLMODEL_DIMENSIONALITY_DATA, 2)
    TRY_MDLOG_GET(nr_bodies, int,  EMDL::MLMODEL_NR_BODIES,           1)
    TRY_MDLOG_GET(is_helix,  bool, EMDL::MLMODEL_IS_HELIX,            false)

    if (is_helix && nr_bodies != 1) {
        REPORT_ERROR("MlModel::readStar: incorrect nr_bodies for helix");
    }

    TRY_MDLOG_GET(helical_nr_asu, int, EMDL::MLMODEL_HELICAL_NR_ASU, 1)

    TRY_MDLOG_GET(helical_twist_min,     RFLOAT, EMDL::MLMODEL_HELICAL_TWIST_MIN,          0.0)
    TRY_MDLOG_GET(helical_twist_max,     RFLOAT, EMDL::MLMODEL_HELICAL_TWIST_MAX,          0.0)
    TRY_MDLOG_GET(helical_twist_inistep, RFLOAT, EMDL::MLMODEL_HELICAL_TWIST_INITIAL_STEP, 0.0)

    TRY_MDLOG_GET(helical_rise_min,      RFLOAT, EMDL::MLMODEL_HELICAL_RISE_MIN,           0.0)
    TRY_MDLOG_GET(helical_rise_max,      RFLOAT, EMDL::MLMODEL_HELICAL_RISE_MAX,           0.0)
    TRY_MDLOG_GET(helical_rise_inistep,  RFLOAT, EMDL::MLMODEL_HELICAL_RISE_INITIAL_STEP,  0.0)

    #undef TRY_MDLOG_GET

    // Treat classes or bodies (for multi-body refinement) in the same way...
    int nr_classes_bodies = nr_bodies > 1 ? nr_bodies : nr_classes;

    if (nr_classes > 1 && nr_bodies > 1)
        REPORT_ERROR("MlModel::readStar: nr_classes and nr_bodies cannot both be greater than 1.");

    // Take inverse again of current resolution:
    current_resolution = 1.0 / current_resolution;

    sigma2_offset *= sigma2_offset;
    sigma2_rot    *= sigma2_rot;
    sigma2_tilt   *= sigma2_tilt;
    sigma2_psi    *= sigma2_psi;

    // Resize vectors
    initialise();

    // Read classes
    MetaDataTable MDclass;
    MDclass.readStar(in, nr_bodies > 1 ? "model_bodies" : "model_classes");

    do_sgd = false;
    for (long int iclass : MDclass) {
        try {
            acc_trans[iclass] = MDclass.getValue<RFLOAT>(EMDL::MLMODEL_ACCURACY_TRANS_ANGSTROM, iclass);
        } catch (const char *errmsg) { try {
            acc_trans[iclass] = MDclass.getValue<RFLOAT>(EMDL::MLMODEL_ACCURACY_TRANS, iclass) * pixel_size;
        } catch (const char *errmsg) {
            REPORT_ERROR("MlModel::readStar: incorrect model_classes/bodies table: no acc_trans");
        } }

        FileName fn_ref;
        try {
            fn_ref          = MDclass.getValue<std::string>(EMDL::MLMODEL_REF_IMAGE, iclass);
            acc_rot[iclass] = MDclass.getValue<RFLOAT>(EMDL::MLMODEL_ACCURACY_ROT, iclass);
        } catch (const char *errmsg) {
            REPORT_ERROR("MlModel::readStar: incorrect model_classes/bodies table: no ref_image or acc_rot");
        }

        // Backwards compatibility
        if (MDlog.containsLabel(EMDL::MLMODEL_ESTIM_RESOL_REF)) {
            estimated_resolution[iclass] = MDlog.getValue<RFLOAT>(EMDL::MLMODEL_ESTIM_RESOL_REF, iclass);
        } else {
            estimated_resolution[iclass] = 0.0;
        }
        if (MDlog.containsLabel(EMDL::MLMODEL_FOURIER_COVERAGE_TOTAL_REF)) {
            total_fourier_coverage[iclass] = MDlog.getValue<RFLOAT>(EMDL::MLMODEL_FOURIER_COVERAGE_TOTAL_REF, iclass);
        } else {
            total_fourier_coverage[iclass] = 0.0;
        }

        if (ref_dim == 2) try {
            XX(prior_offset_class[iclass]) = MDclass.getValue<RFLOAT>(EMDL::MLMODEL_PRIOR_OFFX_CLASS, iclass);
            YY(prior_offset_class[iclass]) = MDclass.getValue<RFLOAT>(EMDL::MLMODEL_PRIOR_OFFY_CLASS, iclass);
        } catch (const char *errmsg) {
            REPORT_ERROR("MlModel::readStar: incorrect model_classes/bodies table: no offset priors for 2D classes");
        }
        if (iclass == 0 || nr_bodies == 1) try {
            // there is only one pdf_class for multibody, but multiple for classification!
            pdf_class[iclass] = MDclass.getValue<RFLOAT>(EMDL::MLMODEL_PDF_CLASS, iclass);
        } catch (const char *errmsg) {
            REPORT_ERROR("MlModel::readStar: incorrect model_classes table: no pdf_class");
        }
        if (is_helix) try {
            helical_rise [iclass] = MDclass.getValue<RFLOAT>(EMDL::MLMODEL_HELICAL_RISE, iclass);
            helical_twist[iclass] = MDclass.getValue<RFLOAT>(EMDL::MLMODEL_HELICAL_TWIST, iclass);
        } catch (const char *errmsg) {
            REPORT_ERROR("MlModel::readStar: incorrect helical parameters");
        }
        if (nr_bodies > 1) {
            keep_fixed_bodies[iclass] = MDclass.containsLabel(EMDL::BODY_KEEP_FIXED) ?
                MDclass.getValue<int>(EMDL::BODY_KEEP_FIXED, iclass) : 0;
        }

        // Read in actual reference image
        Iref[iclass] = Image<RFLOAT>::from_filename(fn_ref)().setXmippOrigin();

        // Check to see whether there is a SGD-gradient entry as well
        try {
            do_sgd = true;
            if (iclass == 0) Igrad.resize(nr_classes);
            Igrad[iclass] = Image<RFLOAT>::from_filename(
                MDclass.getValue<std::string>(EMDL::MLMODEL_SGD_GRADIENT_IMAGE, iclass)
            )();
        } catch (const char *errmsg) {}
    }

    // Read group stuff
    MetaDataTable MDgroup;
    MDgroup.readStar(in, "model_groups");
    // long int optics_group;
    try { for (long int i : MDgroup) {
            long int igroup = MDgroup.getValue<long int>(EMDL::MLMODEL_GROUP_NO, i);
            // Groups are indexed from 1.
            scale_correction      [igroup - 1] = MDgroup.getValue<RFLOAT>(EMDL::MLMODEL_GROUP_SCALE_CORRECTION, i);
            nr_particles_per_group[igroup - 1] = MDgroup.getValue<long>(EMDL::MLMODEL_GROUP_NR_PARTICLES, i);
            group_names           [igroup - 1] = MDgroup.getValue<std::string>(EMDL::MLMODEL_GROUP_NAME, i);
    } } catch (const char *errmsg) {
        REPORT_ERROR("MlModel::readStar: incorrect model_groups table");
    }

    // Read SSNR, noise reduction, tau2_class spectra for each class
    for (int iclass = 0; iclass < nr_classes_bodies; iclass++) {
        MetaDataTable MDsigma;
        MDsigma.readStar(in, (nr_bodies > 1 ? "model_body_" : "model_class_") + integerToString(iclass + 1));
        for (long int i : MDsigma) try {
            const int idx = MDsigma.getValue<int>(EMDL::SPECTRAL_IDX, i);
            data_vs_prior_class   [iclass](idx) = MDsigma.getValue<RFLOAT>(EMDL::MLMODEL_DATA_VS_PRIOR_REF,  i);
            tau2_class            [iclass](idx) = MDsigma.getValue<RFLOAT>(EMDL::MLMODEL_TAU2_REF,  i);
            fsc_halves_class      [iclass](idx) = MDsigma.getValue<RFLOAT>(EMDL::MLMODEL_FSC_HALVES_REF,  i);
            sigma2_class          [iclass](idx) = MDsigma.getValue<RFLOAT>(EMDL::MLMODEL_SIGMA2_REF,  i);
            // Backwards compatibility with STAR files without Fourier coverage
            try {
            fourier_coverage_class[iclass](idx) = MDsigma.getValue<RFLOAT>(EMDL::MLMODEL_FOURIER_COVERAGE_REF,  i);
            } catch (const char *errmsg) {
            fourier_coverage_class[iclass](idx) = 0.0;
            }
        } catch (const char *errmsg) {
            REPORT_ERROR("MlModel::readStar: incorrect table model_class/body_" + integerToString(iclass + 1));
        }
    }

    // Read sigma models for each group
    for (int igroup = 0; igroup < nr_groups; igroup++) {
        // Allow sigma2_noise with different sizes!
        sigma2_noise[igroup].resize(ori_size / 2 + 1);

        if (nr_particles_per_group[igroup] > 0) {
            MetaDataTable MDsigma;
            MDsigma.readStar(in, "model_group_" + integerToString(igroup + 1));
            for (long int i : MDsigma) try {
                const int idx             = MDsigma.getValue<int>(EMDL::SPECTRAL_IDX, i);
                sigma2_noise[igroup](idx) = MDsigma.getValue<RFLOAT>(EMDL::MLMODEL_SIGMA2_NOISE, i);
            } catch (const char *errmsg) {
                REPORT_ERROR("MlModel::readStar: incorrect table model_group_" + integerToString(igroup + 1));
            }
        } else {
            for (auto &x : sigma2_noise[igroup]) { x = 0.0; }
        }
    }

    // Read pdf_direction models for each class
    if (ref_dim == 3) {
        for (int iclass = 0; iclass < nr_classes_bodies; iclass++) {
            MetaDataTable MDclass;
            MDclass.readStar(in, (nr_bodies > 1 ? "model_pdf_orient_body_" : "model_pdf_orient_class_") + integerToString(iclass + 1));
            auto &directions = pdf_direction[iclass];
            directions.clear();
            std::vector<RFLOAT> vdirections;
            vdirections.reserve(MDclass.size());
            for (long int i : MDclass) try {
                vdirections.push_back(MDclass.getValue<RFLOAT>(EMDL::MLMODEL_PDF_ORIENT, i));
            } catch (const char *errmsg) {
                REPORT_ERROR("MlModel::readStar: incorrect table model_pdf_orient_class_" + integerToString(iclass + 1));
            }
            directions.resize(vdirections.size());
            for (long int i = 0; i < Xsize(directions); i++) {
                direct::elem(directions, i) = vdirections[i];
            }
            nr_directions = vdirections.size();
        }
    } else {
        // For 2D case, just fill pdf_direction with ones.
        for (int iclass = 0; iclass < nr_classes_bodies; iclass++) {
            auto &directions = pdf_direction[iclass];
            directions.clear();
            directions.resize(1);
            direct::elem(directions, 0) = 1.0;
        }
        nr_directions = 1;
    }
}

void MlModel::write(FileName fn_out, HealpixSampling &sampling, bool do_write_bild, bool only_write_images) {

    // Treat classes or bodies (for multi-body refinement) in the same way...
    int nr_classes_bodies = nr_bodies > 1 ? nr_bodies : nr_classes;
    // A. Write images
    if (ref_dim == 2) {
        Image<RFLOAT> img(Xsize(Iref[0]), Ysize(Iref[0]), 1, nr_classes_bodies);
        for (int iclass = 0; iclass < nr_classes_bodies; iclass++) {
            for (long int j = 0; j < Ysize(Iref[iclass]); j++)
            for (long int i = 0; i < Xsize(Iref[iclass]); i++) {
                direct::elem(img(), i, j, 0, iclass) = direct::elem(Iref[iclass], i, j);
            }
        }
        img.setSamplingRateInHeader(pixel_size);
        img.write(fn_out + "_" + (nr_bodies > 1 ? "bodies" : "classes") + ".mrcs");

        if (do_sgd) {
            for (int iclass = 0; iclass < nr_classes; iclass++) {
                for (long int j = 0; j < Ysize(Igrad[iclass]); j++)
                for (long int i = 0; i < Xsize(Igrad[iclass]); i++) {
                    direct::elem(img(), i, j, 0, iclass) = direct::elem(Igrad[iclass], i, j);
                }
            }
            img.write(fn_out + "_gradients.mrcs");
        }
    } else {
        // Set correct voxel size in the header
        for (int iclass = 0; iclass < nr_classes_bodies; iclass++) {
            Image<RFLOAT> img (Iref[iclass]);
            img.setSamplingRateInHeader(pixel_size);
            const auto fn = FileName::compose(fn_out + (nr_bodies > 1 ? "_body" : "_class"), iclass + 1, "mrc", 3);
            // apply the body mask for output to the user
            // No! That interferes with a clean continuation of multibody refinement, as ref will be masked 2x then!
            // img() *= masks_bodies[iclass];
            img.write(fn);
        }

        if (do_sgd) {
            for (int iclass = 0; iclass < nr_classes; iclass++) {
                const auto fn = FileName::compose(fn_out + "_grad", iclass + 1, "mrc", 3);
                Image<RFLOAT>(Igrad[iclass]).write(fn);
            }
        }

        if (do_write_bild) {
            // Also write out bild files with the orientational distribution of each class
            // Also write out angular distributions
            // Don't do this for bodies, only for classes!
            for (int iclass = 0; iclass < nr_classes_bodies; iclass++) {
                FileName fn_bild = FileName::compose(fn_out + (nr_bodies > 1 ? "_body" : "_class"), iclass + 1, "", 3) + "_angdist.bild";
                const RFLOAT offset = ori_size * pixel_size / 2.0;
                if (nr_bodies > 1) {
                    // 14 Jul 2017: rotations are all relative to (rot,tilt)=(0,90) to prevent problems with psi-prior around  tilt=0!
                    sampling.writeBildFileOrientationalDistribution(
                        pdf_direction[iclass], fn_bild, offset, offset,
                        &orient_bodies[iclass], &com_bodies[iclass]
                    );
                } else {
                    sampling.writeBildFileOrientationalDistribution(
                        pdf_direction[iclass], fn_bild, offset, offset
                    );
                }
            }
        }
    }

    if (only_write_images)
        return;

    // B. Write STAR file with metadata
    FileName fn_model = fn_out + "_model.star";
    std::ofstream fh (fn_model.c_str(), std::ios::out);
    if (!fh)
        REPORT_ERROR((std::string) "MlModel::write: Cannot write file: " + fn_model);

    // Write the output STAR file
    MetaDataTable MDlog;
    MDlog.isList = true;
    const long int i = MDlog.addObject();
    MDlog.name = "model_general";
    MDlog.setValue(EMDL::MLMODEL_DIMENSIONALITY, ref_dim, i);
    MDlog.setValue(EMDL::MLMODEL_DIMENSIONALITY_DATA, data_dim, i);
    MDlog.setValue(EMDL::MLMODEL_ORIGINAL_SIZE, ori_size, i);
    MDlog.setValue(EMDL::MLMODEL_CURRENT_RESOLUTION, 1.0 / current_resolution, i);
    MDlog.setValue(EMDL::MLMODEL_CURRENT_SIZE, current_size, i);
    MDlog.setValue(EMDL::MLMODEL_PADDING_FACTOR, padding_factor, i);
    MDlog.setValue(EMDL::MLMODEL_IS_HELIX, is_helix, i);
    if (is_helix) {
        MDlog.setValue(EMDL::MLMODEL_HELICAL_NR_ASU, helical_nr_asu, i);
        MDlog.setValue(EMDL::MLMODEL_HELICAL_TWIST_MIN, helical_twist_min, i);
        MDlog.setValue(EMDL::MLMODEL_HELICAL_TWIST_MAX, helical_twist_max, i);
        MDlog.setValue(EMDL::MLMODEL_HELICAL_TWIST_INITIAL_STEP, helical_twist_inistep, i);
        MDlog.setValue(EMDL::MLMODEL_HELICAL_RISE_MIN, helical_rise_min, i);
        MDlog.setValue(EMDL::MLMODEL_HELICAL_RISE_MAX, helical_rise_max, i);
        MDlog.setValue(EMDL::MLMODEL_HELICAL_RISE_INITIAL_STEP, helical_rise_inistep, i);
    }
    MDlog.setValue(EMDL::MLMODEL_INTERPOLATOR, interpolator, i);
    MDlog.setValue(EMDL::MLMODEL_MINIMUM_RADIUS_NN_INTERPOLATION, r_min_nn, i);
    MDlog.setValue(EMDL::MLMODEL_PIXEL_SIZE, pixel_size, i);
    MDlog.setValue(EMDL::MLMODEL_NR_CLASSES, nr_classes, i);
    MDlog.setValue(EMDL::MLMODEL_NR_BODIES, nr_bodies, i);
    MDlog.setValue(EMDL::MLMODEL_NR_GROUPS, nr_groups, i);
    MDlog.setValue(EMDL::MLMODEL_TAU2_FUDGE_FACTOR, tau2_fudge_factor, i);
    MDlog.setValue(EMDL::MLMODEL_NORM_CORRECTION_AVG, avg_norm_correction, i);
    MDlog.setValue(EMDL::MLMODEL_SIGMA_OFFSET_ANGSTROM, sqrt(sigma2_offset), i);
    MDlog.setValue(EMDL::MLMODEL_PRIOR_MODE, orientational_prior_mode, i);
    MDlog.setValue(EMDL::MLMODEL_SIGMA_ROT, sqrt(sigma2_rot), i);
    MDlog.setValue(EMDL::MLMODEL_SIGMA_TILT, sqrt(sigma2_tilt), i);
    MDlog.setValue(EMDL::MLMODEL_SIGMA_PSI, sqrt(sigma2_psi), i);
    MDlog.setValue(EMDL::MLMODEL_LL, LL, i);
    MDlog.setValue(EMDL::MLMODEL_AVE_PMAX, ave_Pmax, i);
    MDlog.write(fh);

    // Calculate resolutions and total Fourier coverages for each class
    calculateTotalFourierCoverage();

    // Write metadata and images for all classes
    MetaDataTable MDclass;
    const FileName fn_root = fn_out.beforeFirstOf("_it");
    MDclass.name = nr_bodies > 1 ? "model_bodies" : "model_classes";
    for (int iclass = 0; iclass < nr_classes_bodies; iclass++) {
        MDclass.addObject();
        Image<RFLOAT> Itmp;
        FileName fn_tmp;
        if (ref_dim == 2) {
            if (nr_bodies > 1) {
                fn_tmp = fn_out + "_bodies.mrcs";
                const auto fn_tmp2 = FileName::compose(fn_root + "_body", iclass + 1, "", 3) + "_mask.mrc";
                // class number from 1 to K!
            } else {
                fn_tmp = fn_out + "_classes.mrcs";
            }
            fn_tmp = FileName::compose(iclass + 1, fn_tmp);  // fn_tmp = integerToString(iclass) + "@" + fn_tmp;
        } else {
            if (nr_bodies > 1) {
                fn_tmp = FileName::compose(fn_out  + "_body", iclass + 1, "mrc", 3); // class number from 1 to K!
                const auto fn_tmp2 = FileName::compose(fn_root + "_body", iclass + 1, "", 3) + "_mask.mrc";
            } else {
                fn_tmp = FileName::compose(fn_out + "_class", iclass + 1, "mrc", 3); // class number from 1 to K!
            }
        }
        MDclass.setValue(EMDL::MLMODEL_REF_IMAGE, fn_tmp, iclass);

        if (do_sgd) {
            fn_tmp = ref_dim == 2 ?
                FileName::compose(iclass + 1, fn_out + "_gradients.mrcs") :
                FileName::compose(fn_out + "_grad", iclass + 1, "mrc", 3);
            MDclass.setValue(EMDL::MLMODEL_SGD_GRADIENT_IMAGE, fn_tmp, iclass);
        }

        // For multiple bodies: only star PDF_CLASS in the first one!
        int myclass = nr_bodies > 1 ? 0 : iclass; // for multi-body: just set iclass=0
        MDclass.setValue(EMDL::MLMODEL_PDF_CLASS, pdf_class[myclass], iclass);
        MDclass.setValue(EMDL::MLMODEL_ACCURACY_ROT, acc_rot[iclass], iclass);
        MDclass.setValue(EMDL::MLMODEL_ACCURACY_TRANS_ANGSTROM, acc_trans[iclass], iclass);
        MDclass.setValue(EMDL::MLMODEL_ESTIM_RESOL_REF, estimated_resolution[iclass], iclass);
        MDclass.setValue(EMDL::MLMODEL_FOURIER_COVERAGE_TOTAL_REF, total_fourier_coverage[iclass], iclass);
        if (nr_bodies > 1) {
            MDclass.setValue(EMDL::BODY_ROTATE_DIRECTION_X, XX(rotate_direction_bodies[iclass]), iclass);
            MDclass.setValue(EMDL::BODY_ROTATE_DIRECTION_Y, YY(rotate_direction_bodies[iclass]), iclass);
            MDclass.setValue(EMDL::BODY_ROTATE_DIRECTION_Z, ZZ(rotate_direction_bodies[iclass]), iclass);
            MDclass.setValue(EMDL::BODY_KEEP_FIXED, keep_fixed_bodies[iclass], iclass);
        }

        if (ref_dim == 2) {
            MDclass.setValue(EMDL::MLMODEL_PRIOR_OFFX_CLASS, XX(prior_offset_class[iclass]), iclass);
            MDclass.setValue(EMDL::MLMODEL_PRIOR_OFFY_CLASS, YY(prior_offset_class[iclass]), iclass);
        }

        if (is_helix) {
            MDclass.setValue(EMDL::MLMODEL_HELICAL_RISE, helical_rise[iclass], iclass);
            MDclass.setValue(EMDL::MLMODEL_HELICAL_TWIST, helical_twist[iclass], iclass);
        }
    }
    MDclass.write(fh);

    // Write radial_average of tau2_class and data_vs_prior_class for each reference
    for (int iclass = 0; iclass < nr_classes_bodies; iclass++) {
        MetaDataTable MDsigma;
        MDsigma.name = nr_bodies > 1 ?
            "model_body_"  + integerToString(iclass + 1) :
            "model_class_" + integerToString(iclass + 1) ;
        for (int i = 0; i < Xsize(tau2_class[iclass]); i++) {
            MDsigma.addObject();
            MDsigma.setValue(EMDL::SPECTRAL_IDX, i, i);
            MDsigma.setValue(EMDL::RESOLUTION, getResolution(i), i);
            MDsigma.setValue(EMDL::RESOLUTION_ANGSTROM, getResolutionAngstrom(i), i);
            MDsigma.setValue(EMDL::MLMODEL_DATA_VS_PRIOR_REF, data_vs_prior_class[iclass](i), i);
            MDsigma.setValue(EMDL::MLMODEL_FSC_HALVES_REF, fsc_halves_class[iclass](i), i);
            MDsigma.setValue(EMDL::MLMODEL_FOURIER_COVERAGE_REF, fourier_coverage_class[iclass](i), i);
            MDsigma.setValue(EMDL::MLMODEL_SIGMA2_REF, sigma2_class[iclass](i), i);
            MDsigma.setValue(EMDL::MLMODEL_TAU2_REF, tau2_class[iclass](i), i);
            // Only write orientabilities if they have been determined
            if (Xsize(orientability_contrib[iclass]) == Xsize(tau2_class[iclass]))
            MDsigma.setValue(EMDL::MLMODEL_ORIENTABILITY_CONTRIBUTION, orientability_contrib[iclass](i), i);
        }
        MDsigma.write(fh);
    }

    // Write scale-correction for all groups
    MetaDataTable MDgroup;
    MDgroup.name = "model_groups";
    for (long int igroup = 0; igroup < nr_groups; igroup++) {
        MDgroup.addObject();
        MDgroup.setValue(EMDL::MLMODEL_GROUP_NO, igroup + 1, igroup);
        MDgroup.setValue(EMDL::MLMODEL_GROUP_NAME, group_names[igroup], igroup);
        MDgroup.setValue(EMDL::MLMODEL_GROUP_NR_PARTICLES, nr_particles_per_group[igroup], igroup);
        MDgroup.setValue(EMDL::MLMODEL_GROUP_SCALE_CORRECTION, scale_correction[igroup], igroup);
    }
    MDgroup.write(fh);

    // Write sigma models for each group
    for (int igroup = 0; igroup < nr_groups; igroup++) {
        if (nr_particles_per_group[igroup] > 0) {
            MetaDataTable MDsigma;
            MDsigma.name = "model_group_" + integerToString(igroup + 1);
            for (int i = 0; i < Xsize(sigma2_noise[igroup]); i++) {
                MDsigma.addObject();
                // Some points in sigma2_noise arrays are never used...
                RFLOAT sigma = sigma2_noise[igroup](i);
                if (sigma > 0.0) {
                    MDsigma.setValue(EMDL::SPECTRAL_IDX, i, i);
                    MDsigma.setValue(EMDL::RESOLUTION, getResolution(i), i);
                    MDsigma.setValue(EMDL::MLMODEL_SIGMA2_NOISE, sigma, i);
                }
            }
            MDsigma.write(fh);
        }
    }

    // Write pdf_direction models for each class
    if (ref_dim == 3) {
        for (int iclass = 0; iclass < nr_classes_bodies; iclass++) {
            MetaDataTable MDclass;
            MDclass.name = nr_bodies > 1 ?
                "model_pdf_orient_body_"  + integerToString(iclass + 1) :
                "model_pdf_orient_class_" + integerToString(iclass + 1) ;
            for (RFLOAT x : pdf_direction[iclass]) {
                MDclass.setValue(EMDL::MLMODEL_PDF_ORIENT, x, MDclass.addObject());
            }
            MDclass.write(fh);
        }
    }
}

void  MlModel::readTauSpectrum(FileName fn_tau, int verb) {
    auto MDtau = MetaDataTable::from_filename(fn_tau);
    RFLOAT val;
    int idx;
    for (long int i : MDtau) {
        idx = MDtau.getValue<int>(EMDL::SPECTRAL_IDX, i);
        val = MDtau.getValue<RFLOAT>(EMDL::MLMODEL_TAU2_REF, i);
        if (idx < Xsize(tau2_class[0]))
            tau2_class[0](idx) = tau2_fudge_factor * val;
    }
    if (idx < Xsize(tau2_class[0]) - 1) {
        if (verb > 0) std::cerr << " Warning: provided tau2-spectrum has fewer entries ("<< idx + 1 << ") than needed (" << Xsize(tau2_class[0]) << "). Set rest to zero..." << std::endl;
    }
    // Use the same spectrum for all classes
    for (int iclass = 0; iclass < nr_classes; iclass++)
        tau2_class[iclass] = tau2_class[0];

}

// Reading images from disc
void MlModel::initialiseFromImages(
    FileName fn_ref, bool _is_3d_model, Experiment &_mydata,
    bool &do_average_unaligned, bool &do_generate_seeds, bool &refs_are_ctf_corrected,
    RFLOAT _ref_angpix, bool do_sgd, bool _do_trust_ref_size, bool verb
) {

    // Data dimensionality
    if (!_mydata.obsModel.opticsMdt.containsLabel(EMDL::IMAGE_DIMENSIONALITY)) {
        if (verb > 0) std::cerr << " WARNING: input particles STAR file does not have a column for image dimensionality, assuming 2D images ..." << std::endl;
        data_dim = 2;
    } else {
        data_dim = _mydata.obsModel.opticsMdt.getValue<int>(EMDL::IMAGE_DIMENSIONALITY, 0);
    }

    // Read references into memory
    Image<RFLOAT> img;
    if (fn_ref != "None") {
        // Read the references into memory
        do_average_unaligned = false;
        // If this is a STAR file, ignore nr_classes and read all references from this file
        if (fn_ref.isStarFile()) {
            auto MDref = MetaDataTable::from_filename(fn_ref, "model_classes");
            const long int i = MDref.size() - 1;
            try {
                FileName fn_tmp = MDref.getValue<std::string>(EMDL::MLMODEL_REF_IMAGE, i);
            } catch (const char *errmsg) {
                // if we did not find the meta-data label _rlnReferenceImage in a directed search, try more generally
                MDref.read(fn_ref);
            }
            try {
                FileName fn_tmp = MDref.getValue<std::string>(EMDL::MLMODEL_REF_IMAGE, i);
            } catch (const char* errmsg) {
                // if we still did not find the meta-data label _rlnReferenceImage, report an error
                REPORT_ERROR("When specifying a .star-file as --ref input, you need to have the _rlnReferenceImage field");
            }

            do_generate_seeds = false;
            // ignore nr_classes from the command line, use number of entries in STAR file
            nr_classes = 0;
            Iref.clear();
            Igrad.clear();
            for (long int i : MDref) {
                img.read(MDref.getValue<std::string>(EMDL::MLMODEL_REF_IMAGE, i));
                img().setXmippOrigin();
                if (_ref_angpix > 0.0) {
                    pixel_size = _ref_angpix;
                } else {
                    RFLOAT header_pixel_size = img.MDMainHeader.getValue<RFLOAT>(EMDL::IMAGE_SAMPLINGRATE_X, img.MDMainHeader.size() - 1);
                    if (nr_classes == 0) {
                        pixel_size = header_pixel_size;
                    } else {
                        if (fabs(header_pixel_size - pixel_size) > 0.001) {
                            REPORT_ERROR("MlModel::readImages ERROR: different models have different pixel sizes in their headers!");
                        }
                    }
                }

                ori_size = Xsize(img());
                ref_dim = img().getDim();
                Iref.push_back(img());
                if (do_sgd) {
                    img() *= 0.0;
                    Igrad.push_back(img());
                }
                nr_classes++;
            }
        } else {
            // For a single image, read this image as reference and set it in all nr_classes Irefs
            img.read(fn_ref);
            img().setXmippOrigin();
            if (_ref_angpix > 0.0) {
                pixel_size = _ref_angpix;
            } else {
                RFLOAT header_pixel_size = img.MDMainHeader.getValue<RFLOAT>(EMDL::IMAGE_SAMPLINGRATE_X, img.MDMainHeader.size() - 1);
                if (header_pixel_size <= 0) {
                    std::cerr << " header_pixel_size = " << header_pixel_size << std::endl;
                    REPORT_ERROR("MlModel::initialiseFromImages: Pixel size of reference image is not set!");
                }
                pixel_size = header_pixel_size;
            }
            ori_size = Xsize(img());
            ref_dim = img().getDim();
            if (ori_size != Xsize(img()) || ori_size != Ysize(img())) {
                std::cerr << " ori_size= " << ori_size << " Xsize(img())= " << Xsize(img()) << std::endl;
                REPORT_ERROR("MlOptimiser::read: size of reference image is not the same as the experimental images!");
            }
            Iref.clear();
            Igrad.clear();
            if (nr_bodies > 1) {
                for (int ibody = 0; ibody < nr_bodies; ibody++) {
                    Iref.push_back(img());
                    if (masks_bodies.size() <= ibody)
                        REPORT_ERROR("BUG: masks_bodies.size() < ibody. Did you initialise the body masks before reading the references?");
                }
            } else {
                for (int iclass = 0; iclass < nr_classes; iclass++) {
                    Iref.push_back(img());
                    if (do_sgd) {
                        img() *= 0.0;
                        Igrad.push_back(img());
                    }
                }
            }
            do_generate_seeds = nr_classes > 1;
        }
    }

    // Make sure that the model has the same box and pixel size as (the first optics group of) the data
    RFLOAT pixel_size_first_optics_group = _mydata.getOpticsPixelSize(0);
    int box_size_first_optics_group      = _mydata.getOpticsImageSize(0);

    if (fn_ref != "None") {

        if (
            fabs(pixel_size - pixel_size_first_optics_group) > 0.001 ||
            ori_size != box_size_first_optics_group
        ) {

            std::string mesg = "";
            if (fabs(pixel_size - pixel_size_first_optics_group) > 0.001) {
                mesg = " The reference pixel size is " + floatToString(pixel_size)
                     + " A/px, but the pixel size of the first optics group of the data is "
                     + floatToString(pixel_size_first_optics_group) + " A/px! \n";
            }
            if (ori_size != box_size_first_optics_group) {
                mesg += " The reference box size is " + integerToString(ori_size)
                     + " px, but the box size of the first optics group of the data is "
                     + integerToString(box_size_first_optics_group) + " px!\n";
            }

            if (!_do_trust_ref_size) {
                REPORT_ERROR("ERROR " + mesg + "\nIf you want to re-scale and/or re-box input particles into the pixel size and the box size of the reference, re-run the program with the --trust_ref_size option.");
            } else if (verb) {
                std::cerr << " WARNING " << mesg;
            }
        }

    } else {
        pixel_size = pixel_size_first_optics_group;
        ori_size = box_size_first_optics_group;

        // Calculate average of all unaligned images later on.
        do_average_unaligned = true;
        do_generate_seeds = false; // after SGD introduction, this is now done in the estimation of initial sigma2 step!
        refs_are_ctf_corrected = true;
        if (_is_3d_model || data_dim == 3) {
            ref_dim = 3;
            img().initZeros(ori_size, ori_size, ori_size);
        } else {
            ref_dim = 2;
            img().initZeros(ori_size, ori_size);
        }
        img().setXmippOrigin();
        Iref.clear();
        Igrad.clear();
        for (int iclass = 0; iclass < nr_classes; iclass++) {
            Iref.push_back(img());
            if (do_sgd)
            Igrad.push_back(img());
        }
    }


    // Set some group stuff
    nr_groups = _mydata.groups.size();
    sigma2_noise.resize(nr_groups, MultidimArray<RFLOAT>::zeros(ori_size / 2 + 1));

    initialise(do_sgd);

    // Now set the group names from the Experiment groups list
    for (int i = 0; i < nr_groups; i++)
        group_names[i] = _mydata.groups[i].name;

}

void MlModel::initialisePdfDirection(long long int newsize) {

    // If the pdf_direction were already filled (size!=0), and newsize=oldsize then leave them as they were
    // If they were still empty, or if the size changes, then initialise them with an even distribution
    for (int iclass = 0; iclass < nr_classes * nr_bodies; iclass++) {
        auto &direction = pdf_direction[iclass];
        long long int oldsize = direction.size();
        if (oldsize == 0 || oldsize != newsize) {
            direction.resize(newsize);
            direction = 1.0 / ((RFLOAT) nr_classes * newsize);
        }
    }
    nr_directions = newsize;

}

void MlModel::initialiseBodies(FileName fn_masks, FileName fn_root_out, bool also_initialise_rest, int rank) {
    auto MD = MetaDataTable::from_filename(fn_masks);
    if (!MD.containsLabel(EMDL::BODY_MASK_NAME))
        REPORT_ERROR("ERROR MlModel::initialiseBodyMasks: body-mask STAR file does not contain rlnBodyMaskName label.");

    const auto n = MD.size();
    masks_bodies           .resize(n);
    com_bodies             .resize(n);
    rotate_direction_bodies.resize(n);
    orient_bodies          .resize(n);
    sigma_tilt_bodies      .resize(n);
    sigma_psi_bodies       .resize(n);
    sigma_offset_bodies    .resize(n);
    keep_fixed_bodies      .resize(n);
    max_radius_mask_bodies .resize(n);
    std::vector<int> relatives_to;
    bool has_rotate_directions = false;
    nr_bodies = 0;
    for (long int i : MD) {
        FileName fn_mask = MD.getValue<std::string>(EMDL::BODY_MASK_NAME, i);
        auto Imask = Image<RFLOAT>::from_filename(fn_mask);
        const auto range = minmax(Imask());
        if (range.first < 0.0 || range.second > 1.0)
            REPORT_ERROR("ERROR: the mask " + fn_mask + " has values outside the range [0,1]");

        Imask().setXmippOrigin();
        masks_bodies[nr_bodies] = Imask();
        Imask.setSamplingRateInHeader(pixel_size);
        // For rotations, find center-of-mass (com)
        int mydim = Imask().getDim();
        Matrix1D<RFLOAT> com(mydim);
        Imask().centerOfMass(com);
        com_bodies[nr_bodies].resize(3);
        for (int i = 0; i < 3; i++) {
            com_bodies[nr_bodies][i] = i >= mydim ? 0.0 : round(com[i]);
            // Round to avoid interpolation artifacts in translate(Iref)
        }
        // Find maximum radius for mask around com
        int max_d2 = 0.0;
        FOR_ALL_ELEMENTS_IN_ARRAY3D(Imask(), i, j, k) {
            if (Imask().elem(i, j, k) > 0.05) {
                int d2 = (k - ZZ(com)) * (k - ZZ(com))
                       + (j - YY(com)) * (j - YY(com))
                       + (i - XX(com)) * (i - XX(com));
                if (d2 > max_d2) { max_d2 = d2; }
            }
        }
        max_radius_mask_bodies[nr_bodies] = ceil(pixel_size * sqrt((RFLOAT) max_d2));

        // Get which body to rotate relative to
        const int relative_to = MD.containsLabel(EMDL::BODY_ROTATE_RELATIVE_TO) ?
            MD.getValue<int>(EMDL::BODY_ROTATE_RELATIVE_TO, i) - 1 : -1;
            // numbering in STAR file starts with 1
        relatives_to.push_back(relative_to);

        if (
            MD.containsLabel(EMDL::BODY_ROTATE_DIRECTION_X) &&
            MD.containsLabel(EMDL::BODY_ROTATE_DIRECTION_Y) &&
            MD.containsLabel(EMDL::BODY_ROTATE_DIRECTION_Z)
        ) {
            has_rotate_directions = true;
            Matrix1D<RFLOAT> body_rotate_direction (3);
            XX(body_rotate_direction) = MD.getValue<RFLOAT>(EMDL::BODY_ROTATE_DIRECTION_X, i);
            YY(body_rotate_direction) = MD.getValue<RFLOAT>(EMDL::BODY_ROTATE_DIRECTION_Y, i);
            ZZ(body_rotate_direction) = MD.getValue<RFLOAT>(EMDL::BODY_ROTATE_DIRECTION_Z, i);
            rotate_direction_bodies.push_back(body_rotate_direction);
        }

        if (MD.containsLabel(EMDL::BODY_SIGMA_ANG)) {
            const RFLOAT sigma = MD.getValue<RFLOAT>(EMDL::BODY_SIGMA_ANG, i);
            sigma_tilt_bodies[nr_bodies] = sigma;
            sigma_psi_bodies [nr_bodies] = sigma;
        } else {
            if (!MD.containsLabel(EMDL::BODY_SIGMA_TILT) || !MD.containsLabel(EMDL::BODY_SIGMA_PSI))  // Logic error?
                REPORT_ERROR("ERROR: either provide rlnBodySigmaAngles OR provide rlnBodySigmaTilt and rlnBodySigmaPsi in the body STAR file.");
            sigma_tilt_bodies[nr_bodies] = MD.getValue<RFLOAT>(EMDL::BODY_SIGMA_TILT, i);
            sigma_psi_bodies [nr_bodies] = MD.getValue<RFLOAT>(EMDL::BODY_SIGMA_PSI,  i);
        }

        try {
            sigma_offset_bodies[nr_bodies] = MD.getValue<RFLOAT>(EMDL::BODY_SIGMA_OFFSET_ANGSTROM, i);
        } catch (const char *errmsg) { try {
            sigma_offset_bodies[nr_bodies] = MD.getValue<RFLOAT>(EMDL::BODY_SIGMA_OFFSET, i) * pixel_size;
        } catch (const char *errmsg) {
            REPORT_ERROR("ERROR: the body STAR file should contain a rlnBodySigmaOffsetAngst column for the prior on the offsets for each body");
        } }

        // Also write the mask with the standard name to disk
        fn_mask = FileName::compose(fn_root_out + "_body", nr_bodies + 1, "", 3) + "_mask.mrc";
        // body number from 1 to K!

        if (rank == 0)
            Imask.write(fn_mask);

        // update counter at the end!
        nr_bodies++;
    }

    // Now that we have the COMs, also get the orientation matrix and the direction of rotation for each body
    for (int ibody = 0; ibody < nr_bodies; ibody++) {
        if (relatives_to[ibody] >= 0) {
            // If another body was given in the input STAR file, rotate this body wrt the COM of the other body
            rotate_direction_bodies[ibody] = com_bodies[relatives_to[ibody]] - com_bodies[ibody];
        } else if (has_rotate_directions) {
            // If the rotation vector is specified directly, just use this one
        } else {
            // if no relative-bodies, nor explicit rotation directions are specified in the STAR file, then rotate relative to (0,0,0)
            rotate_direction_bodies[ibody] = -com_bodies[ibody];
        }

        rotate_direction_bodies[ibody].normalise();
        alignWithZ(-rotate_direction_bodies[ibody], orient_bodies[ibody], false);
    }

    if (also_initialise_rest) {
        if (Iref.size() != 1)
            REPORT_ERROR("BUG: at this point, there should only be a single reference!");

        for (int ibody = 1; ibody < nr_bodies; ibody++) {

            Iref.push_back(Iref[0]);
            tau2_class.push_back(tau2_class[0]);
            fsc_halves_class.push_back(fsc_halves_class[0]);
            sigma2_class.push_back(sigma2_class[0]);
            data_vs_prior_class.push_back(data_vs_prior_class[0]);
            fourier_coverage_class.push_back(fourier_coverage_class[0]);
            acc_rot.push_back(acc_rot[0]);
            acc_trans.push_back(acc_trans[0]);
            estimated_resolution.push_back(estimated_resolution[0]);
            total_fourier_coverage.push_back(total_fourier_coverage[0]);
            if (ref_dim == 2)
            prior_offset_class.push_back(prior_offset_class[0]);
            orientability_contrib.push_back(orientability_contrib[0]);
            PPref.push_back(PPref[0]);
            pdf_direction.push_back(pdf_direction[0]);

            // If all sigmas are zero, ignore this body in the refinement
            keep_fixed_bodies[ibody] =
                sigma_tilt_bodies  [ibody] < 0.001 &&
                sigma_psi_bodies   [ibody] < 0.001 &&
                sigma_offset_bodies[ibody] < 0.001;
        }

        // If provided a specific reference, re-set the corresponding Iref entry
        if (MD.containsLabel(EMDL::BODY_REFERENCE_NAME)) {
            for (long int i : MD) {
                const FileName fn_ref = MD.getValue<std::string>(EMDL::BODY_REFERENCE_NAME, i);
                if (fn_ref != "None")
                    Iref[i] = Image<RFLOAT>::from_filename(fn_ref)().setXmippOrigin();
            }
        }
    }

    // Find the overlap of the bodies, and extend the Iref, PPref and masks_bodies vectors
    pointer_body_overlap.resize(nr_bodies, nr_bodies);
    pointer_body_overlap_inv.resize(nr_bodies);

    // #define DEBUG_OVERLAP
    if (norm_body_mask_overlap) {
        const MultidimArray<RFLOAT> sum_mask = std::accumulate(
            masks_bodies.begin() + 1, masks_bodies.end(), masks_bodies.front());

        for (long int i = 0; i < Xsize(sum_mask); i++)
            if (direct::elem(sum_mask, i) > 1.0)
                for (int ibody = 0; ibody < nr_bodies; ibody++)
                    direct::elem(masks_bodies[ibody], i) /= direct::elem(sum_mask, i);

        for (int ibody = 0; ibody < nr_bodies; ibody++) {
            for (int obody = 0; obody < nr_bodies; obody++)
                direct::elem(pointer_body_overlap, ibody, obody) = obody;
            pointer_body_overlap_inv[ibody] = ibody;
        }

        #ifdef DEBUG_OVERLAP
        for (int ibody = 0; ibody < nr_bodies; ibody++) {
            const FileName fnt = "mask_ibody" + integerToString(ibody) + ".spi";
            Image<RFLOAT>(masks_bodies[ibody]).write(fnt);
            std::cerr << " PPref.size()= " << PPref.size() << std::endl;
        }
        #endif
    } else {
        for (int ibody = 0; ibody < nr_bodies; ibody++) {
            #ifdef DEBUG_OVERLAP
            const FileName fnt = "mask_ibody" + integerToString(ibody) + ".spi";
            Image<RFLOAT>(masks_bodies[ibody]).write(fnt);
            #endif
            for (int obody = 0; obody < nr_bodies; obody++) {
                if (ibody == obody) {
                    direct::elem(pointer_body_overlap, ibody, obody) = obody;
                    pointer_body_overlap_inv[obody] = obody;
                } else {
                    // Sum all the previously done obody masks to see whether there is also overlap with any of them
                    MultidimArray<RFLOAT> overlap_mask = masks_bodies[ibody];
                    for (int oldobody = 0; oldobody < obody; oldobody++) {
                        if (oldobody != ibody) {
                            const int ii = direct::elem(pointer_body_overlap, ibody, oldobody);
                            overlap_mask += masks_bodies[ii];
                        }
                    }
                    // Calculate the overlap between the sum of ibody and all the old obodies until now
                    overlap_mask *= masks_bodies[obody];  // element-wise multiplication
                    // If there is overlap, generate another PPref
                    if (overlap_mask.sum() > 0.0) {
                        // Calculate the mask that has the overlap subtracted from the obody mask
                        overlap_mask = masks_bodies[obody] - overlap_mask;
                        // set the right pointer in the 2D matrix
                        direct::elem(pointer_body_overlap, ibody, obody) = PPref.size();
                        //std::cerr << " ibody= " << ibody << " obody= " << obody << " overlap= " << overlap_mask.sum() << " icc= " << PPref.size() << std::endl;
                        // Extend the two vectors here!
                        PPref.push_back(PPref[obody]);
                        masks_bodies.push_back(overlap_mask);
                        // And keep track of which ibody this entry belonged to
                        pointer_body_overlap_inv.push_back(obody);

                        #ifdef DEBUG_OVERLAP
                        It() = overlap_mask;
                        fnt = "mask_ibody" + integerToString(ibody) + "_obody" + integerToString(obody) + "_overlap.spi";
                        It.write(fnt);
                        std::cerr << " PPref.size()= " << PPref.size() << std::endl;
                        #endif
                    } else {
                        // if there is no overlap: just point to the original obody
                        direct::elem(pointer_body_overlap, ibody, obody) = obody;
                    }
                }
            }
        }
    }
}


void MlModel::writeBildFileBodies(FileName fn_bild) {

    std::ofstream fh_bild (fn_bild.c_str(), std::ios::out);
    if (!fh_bild)
        REPORT_ERROR("HealpixSampling::writeBildFileOrientationalDistribution: cannot open " + fn_bild);

    const RFLOAT xcen = -Xinit(Iref[0]) * pixel_size;
    const RFLOAT ycen = -Yinit(Iref[0]) * pixel_size;
    const RFLOAT zcen = -Zinit(Iref[0]) * pixel_size;
    // Place a black sphere in the centre of the box
    fh_bild << ".color 0 0 0 " << std::endl;
    fh_bild << ".sphere " << xcen << " " << ycen << " " << zcen << " 3 "  << std::endl;
    for (int ibody = 0; ibody < nr_bodies; ibody++) {
        // Sample evenly colors from the rainbow
        RFLOAT r, g, b;
        HSL2RGB((RFLOAT) ibody / (RFLOAT) nr_bodies, 1.0, 0.5, r, g, b);
        fh_bild << ".color " << r << " " << g << " " << b << std::endl;

        const auto &com              = com_bodies[ibody];
        const auto &rotate_direction = rotate_direction_bodies[ibody];

        // Place a sphere at the centre-of-mass
        const RFLOAT x = XX(com) * pixel_size + pixel_size + xcen;
        const RFLOAT y = YY(com) * pixel_size + pixel_size + ycen;
        const RFLOAT z = ZZ(com) * pixel_size + pixel_size + zcen;
        // Add the center of the box to the coordinates
        fh_bild << ".sphere " << x << " " << y << " " << z << " 3 "  << std::endl;
        // Add a label
        fh_bild << ".cmov " << x + 5 << " " << y + 5 << " " << z + 5 << std::endl;
        fh_bild << "body " << ibody + 1 << std::endl;
        // Add an arrow for the direction of the rotation
        const RFLOAT length = 10.0;
        fh_bild << ".arrow " << x << " " << y << " " << z << " "
                << x + length * XX(rotate_direction) * pixel_size << " "
                << y + length * YY(rotate_direction) * pixel_size << " "
                << z + length * ZZ(rotate_direction) * pixel_size << " 1 " << std::endl;
    }
}

void MlModel::setFourierTransformMaps(
    bool update_tau2_spectra, int nr_threads, RFLOAT strict_lowres_exp,
    const MultidimArray<RFLOAT> *fourier_mask
) {

    int min_ires = -1;
    if (strict_lowres_exp > 0) {
        min_ires = round(pixel_size * ori_size / strict_lowres_exp);
        // std::cout << "MlModel::setFourierTransformMaps: strict_lowres_exp = " << strict_lowres_exp
        //           << " pixel_size = " << pixel_size << " ori_size = " << ori_size << " min_ires = " << min_ires << std::endl;;
    }

    // Note that PPref.size() can be bigger than nr_bodies in multi-body refinement, due to extra PPrefs needed for overlapping bodies
    // These only exist in PPref form, they are not needed for reconstructions, only for subtractions in getFourierTransformsAndCtfs
    bool do_heavy = true;
    for (int iclass = 0; iclass < PPref.size(); iclass++) {

        MultidimArray<RFLOAT> Irefp = [&] () {
            if (nr_bodies == 0) return Iref[iclass];

            // ibody deals with overlapping bodies here, as iclass can be larger than nr_bodies when bodies overlap,
            // but there are only nr_bodies Iref; ibody is the number of the original body (max nr_bodies)
            int ibody = pointer_body_overlap_inv[iclass];
            // Place each body with its center-of-mass in the center of the box
            return translate(Iref[ibody] * masks_bodies[iclass], -com_bodies[ibody], DONT_WRAP);
        }();

        if (PPrefRank.size() > 1)
            do_heavy = PPrefRank[iclass];

        MultidimArray<RFLOAT> dummy;
        PPref[iclass].computeFourierTransformMap(Irefp,
            update_tau2_spectra && iclass < nr_classes * nr_bodies ?
                tau2_class[iclass] : dummy,
            current_size, nr_threads, true, do_heavy, min_ires, fourier_mask, do_gpu
        );
    }
}

void MlModel::initialiseDataVersusPrior(bool fix_tau) {

    // Get total number of particles
    const RFLOAT nr_particles = std::accumulate(nr_particles_per_group.begin(), nr_particles_per_group.end(), 0.0);

    // Calculate average sigma2_noise over all image groups
    MultidimArray<RFLOAT> avg_sigma2_noise = MultidimArray<RFLOAT>::zeros(ori_size / 2 + 1);
    MultidimArray<RFLOAT> sum_parts        = MultidimArray<RFLOAT>::zeros(ori_size / 2 + 1);
    for (int igroup = 0; igroup < nr_particles_per_group.size(); igroup++) {
        avg_sigma2_noise += (RFLOAT) nr_particles_per_group[igroup] * sigma2_noise[igroup];
    }
    avg_sigma2_noise /= nr_particles;

    // Get the FT of all reference structures
    // The Fourier Transforms are all "normalised" for 2D transforms of size = ori_size x ori_size
    // And spectrum is squared, so ori_size*ori_size in the 3D case!
    RFLOAT normfft = ref_dim == 3 && data_dim == 2 ? (RFLOAT) (ori_size * ori_size) : 1.0;

    const int nr_classes_bodies = nr_classes * nr_bodies; // also set multiple bodies!
    for (int iclass = 0; iclass < nr_classes_bodies; iclass++) {
        // Initialise output arrays to correct size
        tau2_class[iclass].resize(ori_size / 2 + 1);

        // Get the power spectrum of the reference
        const auto spectrum = getSpectrum(Iref[iclass], POWER_SPECTRUM) * normfft / 2.0;
        // Factor two because of two-dimensionality of the complex plane
        // (just like sigma2_noise estimates, the power spectra should be divided by 2)

        // Update the tau2_class spectrum for this reference
        // This is only for writing out in the it000000_model.star file
        if (!fix_tau) {
            for (long int i = 0; i < Xsize(tau2_class[iclass]); i++) {
                direct::elem(tau2_class[iclass], i) = tau2_fudge_factor * direct::elem(spectrum, i);
            }
        }

        // Calculate data_vs_prior_class as spectral_nr_observations_per_class/sigma2_noise vs 1/tau2_class
        data_vs_prior_class[iclass].resize(ori_size / 2 + 1);
        if (nr_bodies > 1) {
            fsc_halves_class[iclass].initZeros(ori_size / 2 + 1);
        }

        for (long int i = 0; i < Xsize(tau2_class[iclass]); i++) {
            RFLOAT evidence = nr_particles * pdf_class[iclass] / direct::elem(avg_sigma2_noise, i);
            // empirical accounting for ratio of pixels in 3D shells compared to 2D shells
            if (ref_dim == 3 && i > 0)
                evidence /= 2.0 * (RFLOAT) i;
            RFLOAT prior = 1.0 /  direct::elem(tau2_class[iclass], i);
            RFLOAT myssnr = evidence / prior;
            direct::elem(data_vs_prior_class[iclass], i) = myssnr;
            // Also initialise FSC-halves here (...)
            // direct::elem(fsc_halves_class[iclass], i) = myssnr / (myssnr + 1);
        }
    }

}

void MlModel::initialiseHelicalParametersLists(RFLOAT _helical_twist, RFLOAT _helical_rise) {
    if (nr_classes < 1)
        REPORT_ERROR("MlModel::initialiseHelicalParametersLists  nr_classes is smaller than 1");
    helical_twist.resize(nr_classes);
    helical_rise.resize(nr_classes);
    for (int iclass = 0; iclass < nr_classes; iclass++) {
        helical_twist[iclass] = _helical_twist;
        helical_rise [iclass] = _helical_rise;
    }
}

void MlModel::calculateTotalFourierCoverage() {
    for (int iclass = 0; iclass < nr_classes * nr_bodies; iclass++) {
        int maxres = 0;
        for (int ires = 0; ires < Xsize(data_vs_prior_class[iclass]); ires++) {
            if (direct::elem(data_vs_prior_class[iclass], ires) < 1.0)
                break;
            maxres = ires;
        }
        int coverwindow = maxres * 2 - 1;

        estimated_resolution[iclass] = 1.0 / getResolution(maxres);
        total_fourier_coverage[iclass] = 0.0;
        RFLOAT count = 0;
        for (long int k = Xmipp::init(coverwindow); k <= Xmipp::last(coverwindow); k++)
        for (long int i = Xmipp::init(coverwindow); i <= Xmipp::last(coverwindow); i++)
        for (long int j = Xmipp::init(coverwindow); j <= Xmipp::last(coverwindow); j++) {
            int r = sqrt(RFLOAT(k * k + i * i + j * j));
            if (r <= maxres) {
                total_fourier_coverage[iclass] += direct::elem(fourier_coverage_class[iclass], r);
                count += 1.0;
            }
        }
        total_fourier_coverage[iclass] /= count;
    }

}


// MlWsumModel
void MlWsumModel::initialise(MlModel &_model, FileName fn_sym, bool asymmetric_padding, bool _skip_gridding) {
    pixel_size = _model.pixel_size;
    nr_classes = _model.nr_classes;
    nr_bodies = _model.nr_bodies;
    nr_groups = _model.nr_groups;
    nr_directions = _model.nr_directions;
    ref_dim = _model.ref_dim;
    data_dim = _model.data_dim;
    ori_size = _model.ori_size;
    pdf_class = _model.pdf_class;
    if (ref_dim == 2)
        prior_offset_class = _model.prior_offset_class;
    pdf_direction = _model.pdf_direction;
    sigma2_offset = _model.sigma2_offset;
    sigma2_noise = _model.sigma2_noise;
    sigma2_rot = _model.sigma2_rot;
    sigma2_tilt = _model.sigma2_tilt;
    sigma2_psi = _model.sigma2_psi;
    interpolator = _model.interpolator;
    r_min_nn = _model.r_min_nn;
    is_helix = _model.is_helix;
    helical_nr_asu = _model.helical_nr_asu;
    helical_twist_min = _model.helical_twist_min;
    helical_twist_max = _model.helical_twist_max;
    helical_twist_inistep = _model.helical_twist_inistep;
    helical_rise_min = _model.helical_rise_min;
    helical_rise_max = _model.helical_rise_max;
    helical_rise_inistep = _model.helical_rise_inistep;

    padding_factor = _model.padding_factor;
    if (asymmetric_padding)
        padding_factor ++;

    // Don't need forward projectors in MlWsumModel!
    PPref.clear();
    // Don't need scale_correction and bfactor_correction, keep wsum_signal_product and wsum_reference_power instead
    scale_correction.clear();
    bfactor_correction.clear();
    tau2_class.clear();
    data_vs_prior_class.clear();
    acc_rot.clear();
    acc_trans.clear();
    estimated_resolution.clear();
    total_fourier_coverage.clear();
    orientability_contrib.clear();

    helical_twist.resize(nr_classes);
    helical_rise.resize(nr_classes);
    for (int iclass = 0; iclass < nr_classes; iclass++) {
        helical_twist[iclass] = _model.helical_twist[iclass];
        helical_rise[iclass] = _model.helical_rise[iclass];
    }

    wsum_signal_product.resize(nr_groups);
    wsum_reference_power.resize(nr_groups);
    for (long int igroup = 0; igroup < nr_groups; igroup++) {
        wsum_signal_product[igroup] = 0.0;
        wsum_reference_power[igroup] = 0.0;
    }

    // Resize MlWsumModel-specific vectors
    BackProjector BP(
        ori_size, ref_dim, fn_sym, interpolator, padding_factor, r_min_nn,
        ML_BLOB_ORDER, ML_BLOB_RADIUS, ML_BLOB_ALPHA, data_dim, _skip_gridding
    );
    BPref.clear();
    BPref.resize(nr_classes * nr_bodies, BP); // also set multiple bodies
    sumw_group.resize(nr_groups);

}

void MlWsumModel::initZeros() {

    LL = 0.0;
    ave_Pmax = 0.0;
    sigma2_offset = 0.0;
    avg_norm_correction = 0.0;
    sigma2_rot = 0.0;
    sigma2_tilt = 0.0;
    sigma2_psi = 0.0;

    // Set all weighted sums to zero

    for (int iclass = 0; iclass < nr_classes * nr_bodies; iclass++) {
        BPref[iclass].initZeros(current_size);
        // Assume pdf_direction is already of the right size...
        pdf_direction[iclass].initZeros();
    }

    for (int iclass = 0; iclass < nr_classes; iclass++) {
        pdf_class[iclass] = 0.0;
        if (ref_dim == 2)
            prior_offset_class[iclass].initZeros();
    }

    // Initialise sigma2_noise spectra and sumw_group
    for (int igroup = 0; igroup < nr_groups; igroup++) {
        sumw_group[igroup] = 0.0;
        sigma2_noise[igroup].initZeros();
        wsum_signal_product[igroup] = 0.0;
        wsum_reference_power[igroup] = 0.0;
    }
}

// #define DEBUG_PACK
#ifdef DEBUG_PACK
#define MAX_PACK_SIZE 100000
#else
// Approximately 1024 * 1024 * 1024 / 8 / 2 ~ 0.5 Gb
#define MAX_PACK_SIZE 671010000
#endif

void MlWsumModel::pack(MultidimArray<RFLOAT> &packed) {
    unsigned long long packed_size = 0;
    int spectral_size = ori_size / 2 + 1;

    // for LL & avePmax & sigma2_offset & avg_norm_correction & sigma2_rot & sigma2_tilt & sigma2_psi
    packed_size += 7 ;
    // for all group-related stuff
    packed_size += nr_groups * spectral_size;
    // for sumw_group
    packed_size += 3 * nr_groups;
    // for all class-related stuff
    // data is complex: multiply by two!
    packed_size += nr_classes * nr_bodies * 2 * (unsigned long long) BPref[0].getSize();
    packed_size += nr_classes * nr_bodies *     (unsigned long long) BPref[0].getSize();
    packed_size += nr_classes * nr_bodies *     (unsigned long long) nr_directions;
    // for pdf_class
    packed_size += nr_classes;
    // for priors for each class
    if (ref_dim == 2)
        packed_size += nr_classes * 2;

    // Get memory for the packed array
    packed.clear();
    packed.resize(packed_size);

    // Start packing
    unsigned long long idx = 0;

    packed[idx++] = LL;
    packed[idx++] = ave_Pmax;
    packed[idx++] = sigma2_offset;
    packed[idx++] = avg_norm_correction;
    packed[idx++] = sigma2_rot;
    packed[idx++] = sigma2_tilt;
    packed[idx++] = sigma2_psi;

    for (int igroup = 0; igroup < nr_groups; igroup++) {
        for (auto &x : sigma2_noise[igroup]) { packed[idx++] = x; }
        sigma2_noise[igroup].clear();

        packed[idx++] = wsum_signal_product[igroup];
        packed[idx++] = wsum_reference_power[igroup];
        packed[idx++] = sumw_group[igroup];

    }
    for (int iclass = 0; iclass < nr_classes * nr_bodies; iclass++) {

        for (auto &x : BPref[iclass].data) {
            packed[idx++] = x.real;
            packed[idx++] = x.imag;
        }
        BPref[iclass].data.clear();

        for (auto &x : BPref[iclass].weight) {
            packed[idx++] = x;
        }
        BPref[iclass].weight.clear();
        for (auto &x : pdf_direction[iclass]) {
            packed[idx++] = x;
        }
    }
    for (int iclass = 0; iclass < nr_classes; iclass++) {
        pdf_direction[iclass].clear();

        packed[idx++] = pdf_class[iclass];

        if (ref_dim == 2) {
            packed[idx++] = XX(prior_offset_class[iclass]);
            packed[idx++] = YY(prior_offset_class[iclass]);
        }
    }
    #ifdef DEBUG_PACK
    std::cerr << " idx= " << idx << " packed_size= " << packed_size << std::endl;
    #endif

    // Just to check whether we went outside our memory...
    if (idx != packed_size) {
        std::cerr << "idx= " << idx << "packed_size= " << packed_size << std::endl;
        REPORT_ERROR("MlWsumModel::pack: idx != packed_size");
    }

}
void MlWsumModel::unpack(MultidimArray<RFLOAT> &packed) {

    unsigned long long idx = 0;
    int spectral_size = ori_size / 2 + 1;

    LL = packed[idx++];
    ave_Pmax = packed[idx++];
    sigma2_offset = packed[idx++];
    avg_norm_correction = packed[idx++];
    sigma2_rot = packed[idx++];
    sigma2_tilt = packed[idx++];
    sigma2_psi = packed[idx++];

    for (int igroup = 0; igroup < nr_groups; igroup++) {
        sigma2_noise[igroup].resize(spectral_size);
        for (auto &x : sigma2_noise[igroup]) {
            x = packed[idx++];
        }
        wsum_signal_product[igroup] = packed[idx++];
        wsum_reference_power[igroup] = packed[idx++];
        sumw_group[igroup] = packed[idx++];
    }

    for (int iclass = 0; iclass < nr_classes * nr_bodies; iclass++) {
        BPref[iclass].initialiseDataAndWeight(current_size);
        for (auto &x : BPref[iclass].data) {
            x.real = packed[idx++];
            x.imag = packed[idx++];
        }
        for (auto &x : BPref[iclass].weight) {
            x = packed[idx++];
        }
        pdf_direction[iclass].resize(nr_directions);
        for (auto &x : pdf_direction[iclass]) {
            x = packed[idx++];
        }
    }
    for (int iclass = 0; iclass < nr_classes; iclass++) {
        pdf_class[iclass] = packed[idx++];

        if (ref_dim == 2) {
            XX(prior_offset_class[iclass]) = packed[idx++];
            YY(prior_offset_class[iclass]) = packed[idx++];
        }
    }

    unsigned long long packed_size = packed.size();
    packed.clear();

    // Just to check whether we went outside our memory...
    if (idx != packed_size) {
        std::cerr << "idx= " << idx << " packed_size= " << packed_size << std::endl;
        REPORT_ERROR("MlWsumModel::unpack: idx != idx_stop-idx_start");
    }

}


void MlWsumModel::pack(MultidimArray<RFLOAT> &packed, int &piece, int &nr_pieces, bool do_clear) {


    // Determine size of the packed array
    unsigned long long nr_groups = sigma2_noise.size();
    unsigned long long nr_classes_bodies = BPref.size();
    unsigned long long nr_classes = pdf_class.size();
    unsigned long long spectral_size = ori_size / 2 + 1;
    unsigned long long packed_size = 0;
    unsigned long long idx_start, idx_stop;

    // for LL & avePmax & sigma2_offset & avg_norm_correction & sigma2_rot & sigma2_tilt & sigma2_psi
    packed_size += 7;
    // for group-related spectra
    packed_size += nr_groups * spectral_size; // sigma2_noise[spectral_size]
    // for sumw_group
    packed_size += 3 * nr_groups; // wsum_signal_product, wsum_reference_power, sumw_group
    // for all class-related stuff
    // data is complex: multiply by two!
    packed_size += nr_classes_bodies * 2 * (unsigned long long) BPref[0].getSize(); // BPref.data
    packed_size += nr_classes_bodies *     (unsigned long long) BPref[0].getSize(); // BPref.weight
    packed_size += nr_classes_bodies *     (unsigned long long) nr_directions; // pdf_directions
    // for pdf_class
    packed_size += nr_classes;
    // for priors for each class
    if (ref_dim == 2)
        packed_size += nr_classes * 2;

    if (piece < 0 && nr_pieces < 0) {
        // Special case: prevent making multiple pieces if input piece and nr_pieces are both negative
        idx_start = 0;
        idx_stop = packed_size;
    } else if (packed_size > MAX_PACK_SIZE) {
        idx_start = (unsigned long long) piece * MAX_PACK_SIZE;
        idx_stop = std::min(idx_start + MAX_PACK_SIZE, packed_size);
        nr_pieces = ceil((RFLOAT) packed_size / (RFLOAT) MAX_PACK_SIZE);
    } else {
        idx_start = 0;
        idx_stop = packed_size;
        nr_pieces = 1;
    }

    // increment piece so that pack will be called again
    piece++;
    // #define DEBUG_PACK
    #ifdef DEBUG_PACK
    std::cerr << " PACK: idx_start= " << idx_start << " idx_stop= " << idx_stop << " piece= " << piece << " nr_pieces= " << nr_pieces <<" packed_size= "<<packed_size<< std::endl;
    std::cerr << " nr_classes= " << nr_classes << " nr_groups= " << nr_groups << " packed_size= " << packed_size << std::endl;
    std::cerr << " sigma2_noise[0].size()= " << sigma2_noise[0].size() /*<< " wsum_signal_product_spectra[0].size()= " << wsum_signal_product[0].size() << " wsum_reference_power_spectra[0].size()= " << wsum_reference_power[0].size() */<< std::endl;
    std::cerr << " sigma2_noise.size()= " << sigma2_noise.size() << " wsum_signal_product_spectra.size()= " << wsum_signal_product.size() << " wsum_signal_product_spectra.size()= " << wsum_signal_product.size() << std::endl;
    std::cerr << " pdf_direction[0].size()= " << pdf_direction[0].size() << " pdf_direction.size()= " << pdf_direction.size()<<std::endl;
    #endif

    // Get memory for the packed array
    packed.clear();
    packed.resize(idx_stop - idx_start);

    unsigned long long idx = 0;
    unsigned long long ori_idx = 0;
    if (ori_idx >= idx_start && ori_idx < idx_stop) packed[idx++] = LL;
    ori_idx++;
    if (ori_idx >= idx_start && ori_idx < idx_stop) packed[idx++] = ave_Pmax;
    ori_idx++;
    if (ori_idx >= idx_start && ori_idx < idx_stop) packed[idx++] = sigma2_offset;
    ori_idx++;
    if (ori_idx >= idx_start && ori_idx < idx_stop) packed[idx++] = avg_norm_correction;
    ori_idx++;
    if (ori_idx >= idx_start && ori_idx < idx_stop) packed[idx++] = sigma2_rot;
    ori_idx++;
    if (ori_idx >= idx_start && ori_idx < idx_stop) packed[idx++] = sigma2_tilt;
    ori_idx++;
    if (ori_idx >= idx_start && ori_idx < idx_stop) packed[idx++] = sigma2_psi;
    ori_idx++;

    for (int igroup = 0; igroup < nr_groups; igroup++) {
        for (auto &x : sigma2_noise[igroup]) {
            if (ori_idx >= idx_start && ori_idx < idx_stop) { packed[idx++] = x; }
            ori_idx++;
        }
        if (idx == ori_idx && do_clear)
            sigma2_noise[igroup].clear();

        if (ori_idx >= idx_start && ori_idx < idx_stop) packed[idx++] = wsum_signal_product[igroup];
        ori_idx++;
        if (ori_idx >= idx_start && ori_idx < idx_stop) packed[idx++] = wsum_reference_power[igroup];
        ori_idx++;

        if (ori_idx >= idx_start && ori_idx < idx_stop) packed[idx++] = sumw_group[igroup];
        ori_idx++;

    }
    for (int iclass = 0; iclass < nr_classes_bodies; iclass++) {
        for (auto &x : BPref[iclass].data) {
            if (ori_idx >= idx_start && ori_idx < idx_stop) { packed[idx++] = x.real; }
            ori_idx++;
            if (ori_idx >= idx_start && ori_idx < idx_stop) { packed[idx++] = x.imag; }
            ori_idx++;
        }
        // Only clear after the whole array has been packed... i.e. not when we reached the pack_size halfway through
        if (idx == ori_idx && do_clear)
            BPref[iclass].data.clear();

        for (auto &x : BPref[iclass].weight) {
            if (ori_idx >= idx_start && ori_idx < idx_stop) { packed[idx++] = x; }
            ori_idx++;
        }
        if (idx == ori_idx && do_clear)
            BPref[iclass].weight.clear();

        for (auto &x : pdf_direction[iclass]) {
            if (ori_idx >= idx_start && ori_idx < idx_stop) { packed[idx++] = x; }
            ori_idx++;
        }
        if (idx == ori_idx && do_clear)
            pdf_direction[iclass].clear();
    }

    for (int iclass = 0; iclass < nr_classes; iclass++) {

        if (ori_idx >= idx_start && ori_idx < idx_stop)
            packed[idx++] = pdf_class[iclass];
        ori_idx++;

        if (ref_dim == 2) {
            if (ori_idx >= idx_start && ori_idx < idx_stop)
                packed[idx++] = XX(prior_offset_class[iclass]);
            ori_idx++;
            if (ori_idx >= idx_start && ori_idx < idx_stop)
                packed[idx++] = YY(prior_offset_class[iclass]);
            ori_idx++;
        }
    }
    #ifdef DEBUG_PACK
    std::cerr << " idx= " << idx << " packed_size= " << packed_size << std::endl;
    #endif

    // Just to check whether we went outside our memory...
    //std::cerr << " PACK piece= " << piece-1 << " nr_pieces= " << nr_pieces << " ori_idx= " << ori_idx<< " packed_size= " << packed_size << std::endl;
    //std::cerr << " PACK idx= " << idx << " idx_stop-idx_start= " << idx_stop-idx_start << " idx_start= " << idx_start << " idx_stop= " << idx_stop	<< std::endl;
    if (idx != idx_stop - idx_start) {
        std::cerr << "idx= " << idx << "ori_idx= " << ori_idx << " idx_start= " << idx_start << " idx_stop= " << idx_stop << " packed_size= " << packed_size << std::endl;
        REPORT_ERROR("MlWsumModel::pack: idx != idx_stop-idx_start");
    }

}

void MlWsumModel::unpack(MultidimArray<RFLOAT> &packed, int piece, bool do_clear) {

    int nr_groups = sigma2_noise.size();
    int nr_classes_bodies = BPref.size();
    int nr_classes = pdf_class.size();
    int spectral_size = ori_size / 2 + 1;
    unsigned long long idx_start;
    unsigned long long idx_stop;
    if (piece < 0) {
        // Special case: prevent making multiple pieces if input piece is negative
        idx_start = 0;
        idx_stop  = packed.size();
    } else {
        idx_start = (unsigned long long) piece * MAX_PACK_SIZE;
        idx_stop  = idx_start + (unsigned long long) packed.size();
    }
    unsigned long long ori_idx = 0;
    unsigned long long idx = 0;
    #ifdef DEBUG_PACK
    std::cerr << " UNPACK piece= " << piece << " idx_start= " << idx_start << " idx_stop= " << idx_stop << std::endl;
    #endif

    if (ori_idx >= idx_start && ori_idx < idx_stop) LL = packed[idx++];
    ori_idx++;
    if (ori_idx >= idx_start && ori_idx < idx_stop) ave_Pmax = packed[idx++];
    ori_idx++;
    if (ori_idx >= idx_start && ori_idx < idx_stop) sigma2_offset = packed[idx++];
    ori_idx++;
    if (ori_idx >= idx_start && ori_idx < idx_stop) avg_norm_correction = packed[idx++];
    ori_idx++;
    if (ori_idx >= idx_start && ori_idx < idx_stop) sigma2_rot = packed[idx++];
    ori_idx++;
    if (ori_idx >= idx_start && ori_idx < idx_stop) sigma2_tilt = packed[idx++];
    ori_idx++;
    if (ori_idx >= idx_start && ori_idx < idx_stop) sigma2_psi = packed[idx++];
    ori_idx++;

    for (int igroup = 0; igroup < nr_groups; igroup++) {
        if (idx == ori_idx)
            sigma2_noise[igroup].resize(spectral_size);
        for (auto &x : sigma2_noise[igroup]) {
            if (ori_idx >= idx_start && ori_idx < idx_stop) { x = packed[idx++]; }
            ori_idx++;
        }

        if (ori_idx >= idx_start && ori_idx < idx_stop)
            wsum_signal_product[igroup] = packed[idx++];
        ori_idx++;
        if (ori_idx >= idx_start && ori_idx < idx_stop)
            wsum_reference_power[igroup] = packed[idx++];
        ori_idx++;
        if (ori_idx >= idx_start && ori_idx < idx_stop)
            sumw_group[igroup] = packed[idx++];
        ori_idx++;
    }

    for (int iclass = 0; iclass < nr_classes_bodies; iclass++) {
        if (idx == ori_idx)
            BPref[iclass].initialiseDataAndWeight(current_size);
        for (auto &x : BPref[iclass].data) {
            if (ori_idx >= idx_start && ori_idx < idx_stop)
                x.real = packed[idx++];
            ori_idx++;

            if (ori_idx >= idx_start && ori_idx < idx_stop)
                x.imag = packed[idx++];
            ori_idx++;
            // x = Complex(re, im);
        }

        for (auto &x : BPref[iclass].weight) {
            if (ori_idx >= idx_start && ori_idx < idx_stop)
                x = packed[idx++];
            ori_idx++;
        }

        if (idx == ori_idx)
            pdf_direction[iclass].resize(nr_directions);
        for (auto &x : pdf_direction[iclass]) {
            if (ori_idx >= idx_start && ori_idx < idx_stop)
                x = packed[idx++];
            ori_idx++;
        }

    }

    for (int iclass = 0; iclass < nr_classes; iclass++) {
        if (ori_idx >= idx_start && ori_idx < idx_stop)
            pdf_class[iclass] = packed[idx++];
        ori_idx++;

        if (ref_dim == 2) {
            if (ori_idx >= idx_start && ori_idx < idx_stop)
                XX(prior_offset_class[iclass]) = packed[idx++];
            ori_idx++;
            if (ori_idx >= idx_start && ori_idx < idx_stop)
                YY(prior_offset_class[iclass]) = packed[idx++];
            ori_idx++;
        }
    }

    unsigned long long packed_size = packed.size();
    // Free memory
    if (do_clear)
        packed.clear();

    // Just to check whether we went outside our memory...
    //std::cerr << " UNPACK piece= " << piece << " idx= " << idx << " idx_stop-idx_start= " << idx_stop-idx_start << " idx_start= " << idx_start << " idx_stop= " << idx_stop	 << std::endl;
    if (idx != idx_stop - idx_start) {
        std::cerr << "idx= " << idx << "ori_idx= " << ori_idx << " idx_start= " << idx_start << " idx_stop= " << idx_stop << " packed_size= " << packed_size << std::endl;
        REPORT_ERROR("MlWsumModel::unpack: idx != idx_stop-idx_start");
    }
}
