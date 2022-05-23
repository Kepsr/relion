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

        if (fn_out == "" && fn_odir == "")
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
        if (fn_dat != "") {
            std::cout <<" Reading data ..." << std::endl;
            MetaDataTable MDdata;
            MDdata.read(fn_dat);
            optimiser.mydata.MDimg = MDdata;
        }


        // Loop over all micrographs
        int barstep = std::max(1, (int) DFi.numberOfObjects() / 60);
        init_progress_bar(DFi.numberOfObjects());
        long int imgno = 0;
        FileName fn_prevdir = "";
        FOR_ALL_OBJECTS_IN_METADATA_TABLE(DFi) {

            FileName fn_mic = DFi.getValue<FileName>(EMDL::MICROGRAPH_NAME);
            FileName fn_mic_out = fn_out == "" ? fn_mic : fn_mic.insertBeforeExtension("_" + fn_out);
            if (fn_odir != "") {
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
            Image<RFLOAT> Imic_in, Imic_out;
            // Read in the first micrograph
            Imic_in.read(fn_mic);
            Imic_in().setXmippOrigin();
            Imic_out().initZeros(Imic_in());
            MultidimArray<RFLOAT> Imic_sum = MultidimArray<RFLOAT>::zeros(Imic_in());
            Imic_sum.setXmippOrigin();
            // Get mean and stddev of the input micrograph
            // The min and max are not used
            Stats<RFLOAT> stats = Imic_in().computeStats();
            RFLOAT mean_mic = stats.avg;
            RFLOAT stddev_mic = stats.stddev;

            int optics_group_mic = DFi.getValue<int>(EMDL::IMAGE_OPTICS_GROUP);
            RFLOAT mic_pixel_size = -1.0;
            for (int i = 0; i < obsModelMics.opticsMdt.numberOfObjects(); i++) {
                int my_optics_group = obsModelMics.opticsMdt.getValue<int>(EMDL::IMAGE_OPTICS_GROUP);
                if (my_optics_group == optics_group_mic) {
                    mic_pixel_size = obsModelMics.opticsMdt.getValue<RFLOAT>(EMDL::MICROGRAPH_PIXEL_SIZE);
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

                FileName fn_mic2 = optimiser.mydata.MDimg.getValue<FileName>(EMDL::MICROGRAPH_NAME, ori_img_id);
                FileName fn_mic2_pre, fn_mic2_jobnr, fn_mic2_post;
                decomposePipelineFileName(fn_mic2, fn_mic2_pre, fn_mic2_jobnr, fn_mic2_post);

                if (fn_mic2_post == fn_mic_post) {

                    found_one = true;

                    // Prepare transformer
                    MultidimArray<Complex > Fref;
                    MultidimArray<RFLOAT> Mref;
                    if (optimiser.mymodel.data_dim == 3) {
                        Mref.resize(my_image_size, my_image_size, my_image_size);
                        Fref.resize(my_image_size, my_image_size, my_image_size / 2 + 1);
                    } else {
                        Mref.resize(my_image_size, my_image_size);
                        Fref.resize(my_image_size, my_image_size / 2 + 1);
                    }


                    RFLOAT rot, tilt, psi, xcoord = 0.0, ycoord = 0.0, zcoord = 0.0;
                    int iclass;
                    Matrix2D<RFLOAT> A;
                    Matrix1D<RFLOAT> offsets(3);


                    MDcoord.addObject();
                    MDcoord.setObject(optimiser.mydata.MDimg.getObject(ori_img_id));
                    MDcoord.setValue(EMDL::MICROGRAPH_NAME,fn_mic_out);

                    xcoord      = optimiser.mydata.MDimg.getValue<RFLOAT>(EMDL::IMAGE_COORD_X,            ori_img_id);
                    ycoord      = optimiser.mydata.MDimg.getValue<RFLOAT>(EMDL::IMAGE_COORD_Y,            ori_img_id);
                    rot         = optimiser.mydata.MDimg.getValue<RFLOAT>(EMDL::ORIENT_ROT,               ori_img_id);
                    tilt        = optimiser.mydata.MDimg.getValue<RFLOAT>(EMDL::ORIENT_TILT,              ori_img_id);
                    psi         = optimiser.mydata.MDimg.getValue<RFLOAT>(EMDL::ORIENT_PSI,               ori_img_id);
                    XX(offsets) = optimiser.mydata.MDimg.getValue<RFLOAT>(EMDL::ORIENT_ORIGIN_X_ANGSTROM, ori_img_id);
                    YY(offsets) = optimiser.mydata.MDimg.getValue<RFLOAT>(EMDL::ORIENT_ORIGIN_Y_ANGSTROM, ori_img_id);
                    if (optimiser.mymodel.data_dim == 3) {
                    ZZ(offsets) = optimiser.mydata.MDimg.getValue<RFLOAT>(EMDL::ORIENT_ORIGIN_Z_ANGSTROM, ori_img_id);
                    zcoord      = optimiser.mydata.MDimg.getValue<RFLOAT>(EMDL::IMAGE_COORD_Z,            ori_img_id);
                    } else {
                    ZZ(offsets) = zcoord = 0.0;
                    }

                    // Offsets in pixels
                    offsets /= my_pixel_size;

                    iclass = optimiser.mydata.MDimg.getValue<int>(EMDL::PARTICLE_CLASS, ori_img_id) - 1;

                    A = Euler::angles2matrix(rot, tilt, psi);
                    if (do_ctf) {
                        A = optimiser.mydata.obsModel.applyAnisoMag(A, optics_group);
                        A = optimiser.mydata.obsModel.applyScaleDifference(A, optics_group, optimiser.mymodel.ori_size, optimiser.mymodel.pixel_size);
                    }

                    // Get the 2D image (in its ori_size)
                    (optimiser.mymodel.PPref[iclass]).get2DFourierTransform(Fref, A);

                    if (optimiser.mymodel.data_dim == 2) {
                        shiftImageInFourierTransform(Fref, Fref, my_image_size, -XX(offsets), -YY(offsets));
                    } else {
                        shiftImageInFourierTransform(Fref, Fref, my_image_size, -XX(offsets), -YY(offsets), -ZZ(offsets));
                    }


                    if (do_ctf) {
                        MultidimArray<RFLOAT> Fctf;
                        Fctf.resize(Fref);

                        if (optimiser.mymodel.data_dim == 3) {
                            Image<RFLOAT> Ictf;

                            FileName fn_ctf = optimiser.mydata.MDimg.getValue<FileName>(EMDL::CTF_IMAGE, ori_img_id);
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
                                windowFourierTransform(Ictf(), Fctf, Ysize(Fctf));
                            } else {
                                // if dimensions are neither cubical nor FFTW, stop
                                REPORT_ERROR("3D CTF volume must be either cubical or adhere to FFTW format!");
                            }
                        } else {
                            CTF ctf = CTF(optimiser.mydata.MDimg, &optimiser.mydata.obsModel, ori_img_id);
                            Fctf = ctf.getFftwImage(
                                Xsize(Fctf), Ysize(Fctf), my_image_size, my_image_size, my_pixel_size,
                                optimiser.ctf_phase_flipped, false, optimiser.intact_ctf_first_peak, true
                            );
                        }

                        if (optimiser.mydata.obsModel.getCtfPremultiplied(optics_group)) {
                            Fref *= Fctf * Fctf;
                        } else {
                            Fref *= Fctf;
                        }

                        // Also do phase modulation, for beam tilt correction and other asymmetric aberrations
                        optimiser.mydata.obsModel.demodulatePhase(optics_group, Fref, true); // true means do_modulate_instead
                        optimiser.mydata.obsModel.divideByMtf(optics_group, Fref, true); // true means do_multiply_instead

                    }

                    if (optimiser.do_scale_correction) {
                        int group_id = optimiser.mydata.getGroupId(part_id, 0);
                        RFLOAT myscale = optimiser.mymodel.scale_correction[group_id];
                        for (auto &x : Fref) { x *= myscale; }
                    }

                    // Take inverse transform
                    transformer.inverseFourierTransform(Fref, Mref);
                    CenterFFT(Mref, false);
                    Mref.setXmippOrigin();

                    int mic_image_size = ceil(my_image_size * my_pixel_size / mic_pixel_size);
                    MultidimArray<RFLOAT> Mpart_mic = Mref;
                    if (mic_image_size != my_image_size) {
                        resizeMap(Mpart_mic, mic_image_size);
                        Mpart_mic.setXmippOrigin();
                    }

                    //Image<RFLOAT> It;
                    //It()=Mpart_mic;
                    //It.write("It.spi");
                    //exit(1);

                    // To keep raw micrograph and reference projections on the same scale, need to re-obtain
                    // the multiplicative normalisation of the background area (outside circle) again

                    RFLOAT norm_factor = 1.0;
                    if (norm_radius > 0) {
                        Image<RFLOAT> Ipart;
                        Ipart().resize(Mpart_mic);
                        Ipart().initConstant(mean_mic); // set areas outside the micrograph to average of micrograph (just like in preprocessing)
                        Imic_in().xinit = -round(xcoord);
                        Imic_in().yinit = -round(ycoord);
                        Imic_in().zinit = -round(zcoord);
                        FOR_ALL_ELEMENTS_IN_ARRAY3D(Mpart_mic) {
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
                    Imic_sum.xinit = -round(xcoord);
                    Imic_sum.yinit = -round(ycoord);
                    Imic_sum.zinit = -round(zcoord);
                    radius = optimiser.particle_diameter / (2.0 * mic_pixel_size);
                    FOR_ALL_ELEMENTS_IN_ARRAY3D(Mpart_mic) {
                        long int idx = round(sqrt(i * i + j * j + k * k));
                        if (idx < radius) {
                            // check the particles do not go off the side
                            int ip = i - Xinit(Imic_sum);
                            int jp = j - Yinit(Imic_sum);
                            int kp = k - Zinit(Imic_sum);
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
                        Imic_in()[n] -= mean_mic;
                        Imic_in()[n] /= stddev_mic;
                        // And add a precentage to Imic_out
                        Imic_out()[n] *= 1.0 - micrograph_background;
                        Imic_out()[n] += micrograph_background * Imic_in()[n];
                    }
                }

                // Write out the new micrograph
                Imic_out.write(fn_mic_out);

                MDmics_out.addObject();
                MDmics_out.setObject(DFi.getObject());
                MDmics_out.setValue(EMDL::MICROGRAPH_NAME, fn_mic_out);

                // Also write out a STAR file with the particles used
                FileName fn_coord_out = fn_mic_out.withoutExtension()+ "_coord.star";
                MDcoord.write(fn_coord_out);
                MDcoord.clear();
            } else {
                MDmics_out.addObject();
                MDmics_out.setObject(DFi.getObject());
            }

            if (imgno % barstep == 0) progress_bar(imgno);
            imgno++;

        } // end loop over input MetadataTable
        progress_bar(DFi.numberOfObjects());


        FileName fn_star_out = fn_odir + "micrographs_reposition.star";
        if (fn_out != "") fn_star_out = fn_star_out.insertBeforeExtension("_" + fn_out);
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
