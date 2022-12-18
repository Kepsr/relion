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

#include <src/args.h>
#include <src/ml_optimiser.h>
#include <src/jaz/obs_model.h>
#include "src/jaz/ctf_helper.h"
#include <stdlib.h>

class particle_reposition_parameters {

    public:

    FileName fn_in, fn_opt, fn_out, fn_dat, fn_odir;

    RFLOAT micrograph_background;
    int norm_radius;
    bool do_invert, do_ctf, do_subtract;
    ObservationModel obsModelMics;

    // I/O Parser
    IOParser parser;

    void usage() {
        parser.writeUsage(std::cerr);
    }

    void read(int argc, char **argv) {
        parser.setCommandLine(argc, argv);

        int general_section = parser.addSection("Options");

        fn_in  = parser.getOption("--i", "Input STAR file with rlnMicrographName's ");
        fn_out = parser.getOption("--o", "Output rootname, to be added to input micrograph names", "");
        fn_odir = parser.getOption("--odir", "Output directory (default is same as input micrographs directory", "");
        fn_opt = parser.getOption("--opt", "Optimiser STAR file with the 2D classes or 3D maps to be repositioned");
        fn_dat = parser.getOption("--data", "Data STAR file with selected particles (default is to use all particles)", "");
        micrograph_background = textToFloat(parser.getOption("--background", "The fraction of micrograph background noise in the output micrograph", "0.1"));
        do_invert= parser.checkOption("--invert", "Invert the contrast in the references?");

        do_ctf = parser.checkOption("--ctf", "Apply CTF for each particle to the references?");
        norm_radius = textToFloat(parser.getOption("--norm_radius", "Radius of the circle used for background normalisation (in pixels)", "-1"));
        do_subtract = parser.checkOption("--subtract", "Subtract repositioned micrographs from the input ones?");

        // Check for errors in the command-line option
        if (parser.checkForErrors())
            REPORT_ERROR("Errors encountered on the command line (see above), exiting...");
    }

    void run() {

        if (fn_out.empty() && fn_odir.empty())
            REPORT_ERROR("ERROR: You need to provide either --o or --odir");

        if (fn_odir.length() > 0 && fn_odir[fn_odir.length()-1] != '/') fn_odir += "/";

        int xdim, ydim, radius;
        MetaDataTable DFi, DFopt, MDmics_out;
        ObservationModel::loadSafely(fn_in, obsModelMics, DFi, "micrographs");

        MlOptimiser optimiser;
        optimiser.do_preread_images = false;

        optimiser.read(fn_opt);
        optimiser.mymodel.setFourierTransformMaps(false);

        // Use a user-provided subset of particles instead of all of them?
        if (!fn_dat.empty()) {
            std::cout <<" Reading data ..." << std::endl;
            MetaDataTable MDdata;
            MDdata.read(fn_dat);
            optimiser.mydata.MDimg = MDdata;
        }


        // Loop over all micrographs
        int barstep = std::max(1, (int) DFi.size() / 60);
        init_progress_bar(DFi.size());
        long int imgno = 0;
        FileName fn_prevdir = "";
        for (long int i : DFi) {

            FileName fn_mic = DFi.getValue<std::string>(EMDL::MICROGRAPH_NAME, i);
            FileName fn_mic_out = fn_out.empty() ? fn_mic : fn_mic.insertBeforeExtension("_" + fn_out);
            if (!fn_odir.empty()) {
                FileName fn_pre, fn_jobnr, fn_post;
                if (decomposePipelineFileName(fn_mic_out, fn_pre, fn_jobnr, fn_post)) {
                    fn_mic_out = fn_odir + fn_post;
                } else {
                    fn_mic_out = fn_odir + fn_mic_out;
                }
                FileName fn_onlydir = fn_mic_out.beforeLastOf("/");
                if (fn_onlydir != fn_prevdir) {
                    std::string command = " mkdir -p " + fn_onlydir;
                    int res = system(command.c_str());
                    fn_prevdir = fn_onlydir;
                }
            }

            FourierTransformer transformer;
            MetaDataTable MDcoord;
            // Read in the first micrograph
            auto Imic_in = Image<RFLOAT>::from_filename(fn_mic);
            Imic_in().setXmippOrigin();
            Image<RFLOAT> Imic_out (MultidimArray<RFLOAT>::zeros(Imic_in()));
            auto Imic_sum = MultidimArray<RFLOAT>::zeros(Imic_in());
            Imic_sum.setXmippOrigin();
            // Get mean and stddev of the input micrograph
            const auto micrograph_stats = computeStats(Imic_in());

            const int optics_group_mic = DFi.getValue<int>(EMDL::IMAGE_OPTICS_GROUP, i);
            RFLOAT mic_pixel_size = -1.0;
            for (int j = 0; j < obsModelMics.opticsMdt.size(); j++) {
                const int my_optics_group = obsModelMics.opticsMdt.getValue<int>(EMDL::IMAGE_OPTICS_GROUP, j);
                if (my_optics_group == optics_group_mic) {
                    mic_pixel_size = obsModelMics.opticsMdt.getValue<RFLOAT>(EMDL::MICROGRAPH_PIXEL_SIZE, j);
                    break;
                }
            }

            if (mic_pixel_size < 0.0)
                REPORT_ERROR("ERROR: could not find correct optics group in micrograph star file...");

            FileName fn_mic_pre, fn_mic_jobnr, fn_mic_post;
            decomposePipelineFileName(fn_mic, fn_mic_pre, fn_mic_jobnr, fn_mic_post);

            // Loop over all particles in the mydata.MDimg table
            bool found_one = false;
            for (long int part_id = 0; part_id < optimiser.mydata.numberOfParticles(); part_id++) {
                long int ori_img_id = optimiser.mydata.particles[part_id].images[0].id;
                int optics_group = optimiser.mydata.getOpticsGroup(part_id, 0);
                RFLOAT my_pixel_size = optimiser.mydata.getImagePixelSize(part_id, 0);
                int my_image_size = optimiser.mydata.getOpticsImageSize(optics_group);

                if (do_subtract && fabs(my_pixel_size - mic_pixel_size) > 1e-6)
                    REPORT_ERROR("ERROR: subtract code has only been validated with same pixel size for particles and micrographs... Sorry!");

                FileName fn_mic2 = optimiser.mydata.MDimg.getValue<std::string>(EMDL::MICROGRAPH_NAME, ori_img_id);
                FileName fn_mic2_pre, fn_mic2_jobnr, fn_mic2_post;
                decomposePipelineFileName(fn_mic2, fn_mic2_pre, fn_mic2_jobnr, fn_mic2_post);

                if (fn_mic2_post == fn_mic_post) {

                    found_one = true;

                    // Prepare transformer
                    MultidimArray<RFLOAT> Mref (my_image_size, my_image_size,
                        optimiser.mymodel.data_dim == 3 ? my_image_size : 1);

                    Matrix<RFLOAT> A;
                    Vector<RFLOAT> offsets (3);

                    const long int i = MDcoord.addObject();
                    MDcoord.setObject(optimiser.mydata.MDimg.getObject(ori_img_id), i);
                    MDcoord.setValue(EMDL::MICROGRAPH_NAME, fn_mic_out, i);

                    RFLOAT xcoord      = optimiser.mydata.MDimg.getValue<RFLOAT>(EMDL::IMAGE_COORD_X,            ori_img_id);
                    RFLOAT ycoord      = optimiser.mydata.MDimg.getValue<RFLOAT>(EMDL::IMAGE_COORD_Y,            ori_img_id);
                    RFLOAT zcoord = 0.0;
                    RFLOAT rot         = optimiser.mydata.MDimg.getValue<RFLOAT>(EMDL::ORIENT_ROT,               ori_img_id);
                    RFLOAT tilt        = optimiser.mydata.MDimg.getValue<RFLOAT>(EMDL::ORIENT_TILT,              ori_img_id);
                    RFLOAT psi         = optimiser.mydata.MDimg.getValue<RFLOAT>(EMDL::ORIENT_PSI,               ori_img_id);
                    XX(offsets) = optimiser.mydata.MDimg.getValue<RFLOAT>(EMDL::ORIENT_ORIGIN_X_ANGSTROM, ori_img_id);
                    YY(offsets) = optimiser.mydata.MDimg.getValue<RFLOAT>(EMDL::ORIENT_ORIGIN_Y_ANGSTROM, ori_img_id);
                    if (optimiser.mymodel.data_dim == 3) {
                    ZZ(offsets) = optimiser.mydata.MDimg.getValue<RFLOAT>(EMDL::ORIENT_ORIGIN_Z_ANGSTROM, ori_img_id);
                    zcoord = optimiser.mydata.MDimg.getValue<RFLOAT>(EMDL::IMAGE_COORD_Z,            ori_img_id);
                    } else {
                    ZZ(offsets) = zcoord = 0.0;
                    }

                    // Offsets in pixels
                    offsets /= my_pixel_size;

                    const int iclass = optimiser.mydata.MDimg.getValue<int>(EMDL::PARTICLE_CLASS, ori_img_id) - 1;

                    A = Euler::angles2matrix(rot, tilt, psi);
                    if (do_ctf) {
                        if (optimiser.mydata.obsModel.hasMagMatrices)
                            A = A.matmul(optimiser.mydata.obsModel.anisoMag(optics_group));
                        A *= optimiser.mydata.obsModel.scaleDifference(optics_group, optimiser.mymodel.ori_size, optimiser.mymodel.pixel_size);
                    }

                    // Get the 2D image (in its ori_size)
                    auto Fref = optimiser.mymodel.PPref[iclass].get2DFourierTransform(
                        my_image_size / 2 + 1, my_image_size,
                        optimiser.mymodel.data_dim == 3 ? my_image_size : 1, A);

                    shiftImageInFourierTransform(
                        Fref, my_image_size, -XX(offsets), -YY(offsets),
                        optimiser.mymodel.data_dim == 2 ? 0 : -ZZ(offsets)
                    );

                    if (do_ctf) {
                        MultidimArray<RFLOAT> Fctf;
                        Fctf.resize(Fref);

                        if (optimiser.mymodel.data_dim == 3) {

                            FileName fn_ctf = optimiser.mydata.MDimg.getValue<std::string>(EMDL::CTF_IMAGE, ori_img_id);
                            auto Ictf = Image<RFLOAT>::from_filename(fn_ctf);

                            // If there is a redundant half, get rid of it
                            if (Xsize(Ictf()) == Ysize(Ictf())) {
                                Ictf().setXmippOrigin();
                                FOR_ALL_ELEMENTS_IN_FFTW_TRANSFORM(Fctf) {
                                    // Use negative ip, jp, kp indices
                                    // because the origin in the ctf_img lies half a pixel to the right of the actual center.
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
                            const CTF ctf = CtfHelper::makeCTF(optimiser.mydata.MDimg, &optimiser.mydata.obsModel, ori_img_id);
                            Fctf = CtfHelper::getFftwImage(
                                ctf,
                                Xsize(Fctf), Ysize(Fctf), my_image_size, my_image_size, my_pixel_size,
                                &optimiser.mydata.obsModel,
                                optimiser.ctf_phase_flipped, false, optimiser.intact_ctf_first_peak, true
                            );
                        }

                        Fref *= optimiser.mydata.obsModel.getCtfPremultiplied(optics_group) ?
                            Fctf * Fctf :  // Expression templates
                            Fctf;

                        // Also do phase modulation, for beam tilt correction and other asymmetric aberrations
                        optimiser.mydata.obsModel.modulatePhase(optics_group, Fref);
                        optimiser.mydata.obsModel.multiplyByMtf(optics_group, Fref);

                    }

                    if (optimiser.do_scale_correction) {
                        const int group_id = optimiser.mydata.getGroupId(part_id, 0);
                        const RFLOAT scale = optimiser.mymodel.scale_correction[group_id];
                        Fref *= scale;
                    }

                    // Take inverse transform
                    Mref = transformer.inverseFourierTransform(Fref);
                    CenterFFT(Mref, -1);
                    Mref.setXmippOrigin();

                    int mic_image_size = ceil(my_image_size * my_pixel_size / mic_pixel_size);
                    MultidimArray<RFLOAT> Mpart_mic = Mref;
                    if (mic_image_size != my_image_size) {
                        resizeMap(Mpart_mic, mic_image_size);
                        Mpart_mic.setXmippOrigin();
                    }

                    // Image<RFLOAT>(Mpart_mic).write("It.spi");
                    // exit(1);

                    // To keep raw micrograph and reference projections on the same scale, need to re-obtain
                    // the multiplicative normalisation of the background area (outside circle) again

                    RFLOAT norm_factor = 1.0;
                    if (norm_radius > 0) {
                        Image<RFLOAT> Ipart;
                        Ipart().resize(Mpart_mic);
                        Ipart() = micrograph_stats.avg; // set areas outside the micrograph to average of micrograph (just like in preprocessing)
                        Imic_in().xinit = -round(xcoord);
                        Imic_in().yinit = -round(ycoord);
                        Imic_in().zinit = -round(zcoord);
                        FOR_ALL_ELEMENTS_IN_ARRAY3D(Mpart_mic, i, j, k) {
                            // check the particles do not go off the side
                            int ip = i - Xinit(Imic_in());
                            int jp = j - Yinit(Imic_in());
                            int kp = k - Zinit(Imic_in());
                            if (
                                ip >= 0 && ip < Xsize(Imic_in()) &&
                                jp >= 0 && jp < Ysize(Imic_in()) &&
                                kp >= 0 && kp < Zsize(Imic_in())
                            ) {
                                Ipart().elem(i, j, k) = Imic_in().elem(i, j, k);
                            }
                        }

                        RFLOAT psi_deg = 0.0, tilt_deg = 90.0;
                        if (optimiser.do_helical_refine) {
                            tilt_deg = optimiser.mydata.MDimg.getValue<RFLOAT>(EMDL::ORIENT_TILT_PRIOR, ori_img_id);
                            psi_deg  = optimiser.mydata.MDimg.getValue<RFLOAT>(EMDL::ORIENT_PSI_PRIOR,  ori_img_id);
                        }

                        norm_factor = calculateBackgroundAvgStddev(
                            Ipart, norm_radius, optimiser.do_helical_refine,
                            optimiser.helical_tube_outer_diameter / (2.0 * mic_pixel_size), tilt_deg, psi_deg
                        ).stddev;  // Could do with just calculating stddev.

                        // Apply the per-particle norm_correction term
                        if (optimiser.do_norm_correction) {
                            // TODO: check whether this is the right way around!!!
                            norm_factor *= optimiser.mydata.MDimg.getValue<RFLOAT>(EMDL::IMAGE_NORM_CORRECTION, ori_img_id)
                                         / optimiser.mymodel.avg_norm_correction;
                        }
                    }

                    // Reposition Mpart_mic back into the micrograph
                    Imic_out().xinit = -round(xcoord);
                    Imic_out().yinit = -round(ycoord);
                    Imic_out().zinit = -round(zcoord);
                    Imic_sum.xinit   = -round(xcoord);
                    Imic_sum.yinit   = -round(ycoord);
                    Imic_sum.zinit   = -round(zcoord);
                    radius = optimiser.particle_diameter / (2.0 * mic_pixel_size);
                    FOR_ALL_ELEMENTS_IN_ARRAY3D(Mpart_mic, i, j, k) {
                        long int idx = round(euclid(i, j, k));
                        if (idx < radius) {
                            // check the particles do not go off the side
                            const int ip = i - Xinit(Imic_sum);
                            const int jp = j - Yinit(Imic_sum);
                            const int kp = k - Zinit(Imic_sum);
                            if (
                                ip >= 0 && ip < Xsize(Imic_sum) &&
                                jp >= 0 && jp < Ysize(Imic_sum) &&
                                kp >= 0 && kp < Zsize(Imic_sum)
                            ) {
                                Imic_out().elem(i, j, k) += norm_factor * Mpart_mic.elem(i, j, k);
                                Imic_sum.elem(i, j, k) += 1.0;
                            }
                        }
                    }
                }
            }

            if (found_one) {
                for (long int n = 0; n < Imic_out().size(); n++) {
                    if (Imic_sum[n] > 0.0)
                        Imic_out()[n] /= Imic_sum[n];
                    if (do_invert)
                        Imic_out()[n] *= -1.0;
                    if (do_subtract) {
                        Imic_out()[n] = Imic_in()[n] - Imic_out()[n];
                    } else if (micrograph_background > 0.0) {
                        // normalize Imic_in on the fly
                        Imic_in()[n] -= micrograph_stats.avg;
                        Imic_in()[n] /= micrograph_stats.stddev;
                        // And add a precentage to Imic_out
                        Imic_out()[n] *= 1.0 - micrograph_background;
                        Imic_out()[n] += micrograph_background * Imic_in()[n];
                    }
                }

                // Write out the new micrograph
                Imic_out.write(fn_mic_out);

                const long int j = MDmics_out.addObject();
                MDmics_out.setObject(DFi.getObject(i), j);
                MDmics_out.setValue(EMDL::MICROGRAPH_NAME, fn_mic_out, j);

                // Also write out a STAR file with the particles used
                const FileName fn_coord_out = fn_mic_out.withoutExtension() + "_coord.star";
                MDcoord.write(fn_coord_out);
                MDcoord.clear();
            } else {
                const long int j = MDmics_out.addObject();
                MDmics_out.setObject(DFi.getObject(i), j);
            }

            if (imgno % barstep == 0) progress_bar(imgno);
            imgno++;

        } // end loop over input MetadataTable
        progress_bar(DFi.size());

        FileName fn_star_out = fn_odir + "micrographs_reposition.star";
        if (!fn_out.empty()) fn_star_out = fn_star_out.insertBeforeExtension("_" + fn_out);
        std::cout << "Writing out star file with the new micrographs: " << fn_star_out << std::endl;
        obsModelMics.save(MDmics_out, fn_star_out, "micrographs");
        std::cout << " Done!" << std::endl;
    }
};


int main(int argc, char *argv[]) {
    time_config();
    particle_reposition_parameters prm;

    try {
        prm.read(argc, argv);
        prm.run();
    } catch (RelionError XE) {
        // prm.usage();
        std::cerr << XE;
        return RELION_EXIT_FAILURE;
    }

    return RELION_EXIT_SUCCESS;
}
