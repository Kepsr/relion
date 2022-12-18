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

#include <src/image.h>
#include <src/funcs.h>
#include <src/ctf.h>
#include <src/args.h>
#include <src/error.h>
#include <src/euler.h>
#include <src/time.h>
#include <src/jaz/obs_model.h>

// TODO: set pixel sizes in the outputs

class stack_create_parameters {

    public:

    FileName fn_star, fn_root, fn_ext;
    MetaDataTable MD;
    // I/O Parser
    IOParser parser;
    bool do_spider, do_split_per_micrograph, do_apply_trans, do_apply_trans_only, do_ignore_optics, do_one_by_one;
    ObservationModel obsModel;

    void usage() { parser.writeUsage(std::cerr); }

    void read(int argc, char **argv) {
        parser.setCommandLine(argc, argv);

        int general_section = parser.addSection("General options");
        fn_star = parser.getOption("--i", "Input STAR file with the images (as rlnImageName) to be saved in a stack");
        fn_root = parser.getOption("--o", "Output rootname","output");
        do_spider  = parser.checkOption("--spider_format", "Write out in SPIDER stack format (by default MRC stack format)");
        do_split_per_micrograph = parser.checkOption("--split_per_micrograph", "Write out separate stacks for each micrograph (needs rlnMicrographName in STAR file)");
        do_apply_trans = parser.checkOption("--apply_transformation", "Apply the inplane-transformations (needs _rlnOriginX/Y and _rlnAnglePsi in STAR file) by real space interpolation");
        do_apply_trans_only = parser.checkOption("--apply_rounded_offsets_only", "Apply the rounded translations only (so-recentering without interpolation; needs _rlnOriginX/Y in STAR file)");
        do_ignore_optics = parser.checkOption("--ignore_optics", "Ignore optics groups. This allows you to read and write RELION 3.0 STAR files but does NOT allow you to convert 3.1 STAR files back to the 3.0 format.");
        do_one_by_one = parser.checkOption("--one_by_one", "Write particles one by one. This saves memory but can be slower.");

        if (do_apply_trans)
            std::cerr << "WARNING: --apply_transformation uses real space interpolation. It also invalidates CTF parameters (e.g. beam tilt & astigmatism). This can degrade the resolution. USE WITH CARE!!" << std::endl;

        fn_ext = (do_spider) ? ".spi" : ".mrcs";

        // Check for errors in the command-line option
        if (parser.checkForErrors())
                REPORT_ERROR("Errors encountered on the command line, exiting...");
    }

    void run() {
        if (do_ignore_optics && (do_apply_trans || do_apply_trans_only))
            REPORT_ERROR("ERROR: you cannot ignore optics and apply transformations");

        if (do_ignore_optics) {
            MD.read(fn_star);
        } else {
            ObservationModel::loadSafely(fn_star, obsModel, MD, "particles");
        }

        // Check for rlnImageName label
        if (!MD.containsLabel(EMDL::IMAGE_NAME))
            REPORT_ERROR("ERROR: Input STAR file does not contain the rlnImageName label. Aren't you reading RELION 3.1 STAR files with --ignore_optics?");

        if (do_split_per_micrograph && !MD.containsLabel(EMDL::MICROGRAPH_NAME))
            REPORT_ERROR("ERROR: Input STAR file does not contain the rlnMicrographName label");

        Image<RFLOAT> in;
        FileName fn_mic;
        std::vector<FileName> fn_mics;
        std::vector<int> mics_ndims;

        // First get number of images and their size
        bool is_first = true;
        int xdim, ydim, zdim;
        for (long int ndim : MD) {
            if (is_first) {
                const FileName fn_img = MD.getValue<std::string>(EMDL::IMAGE_NAME, ndim);
                in.read(fn_img);
                xdim = in().xdim;
                ydim = in().ydim;
                zdim = in().zdim;
                is_first = false;
            }

            if (do_split_per_micrograph) {
                fn_mic = MD.getValue<std::string>(EMDL::MICROGRAPH_NAME, ndim);
                bool have_found = false;
                for (int m = 0; m < fn_mics.size(); m++) {
                    if (fn_mic == fn_mics[m]) {
                        have_found = true;
                        mics_ndims[m]++;
                        break;
                    }
                }
                if (!have_found) {
                    fn_mics.push_back(fn_mic);
                    mics_ndims.push_back(1);
                }
            }
        }

        // If not splitting, just fill fn_mics and mics_ndim with one entry (to re-use loop below)
        if (!do_split_per_micrograph) {
            fn_mics.emplace_back();
            mics_ndims.push_back(MD.size());
        }

        // Loop over all micrographs
        int ndim = 0;
        for (int m = 0; m < fn_mics.size(); m++) {
            ndim = mics_ndims[m];
            fn_mic = fn_mics[m];

            Image<RFLOAT> out;

            if (!do_one_by_one) {
                // Resize the output image
                std::cout << "Resizing the output stack to "<< ndim << " images of size: " << xdim << "x" << ydim << "x" << zdim << std::endl;
                out().reshape(xdim, ydim, zdim, ndim);
                RFLOAT Gb = RFLOAT(out().size() * sizeof(RFLOAT)) / (1024.0 * 1024.0 * 1024.0);
                std::cout << "This will require " << Gb << "Gb of memory...." << std::endl;
                std::cout << "If this runs out of memory, please try --one_by_one." << std::endl;
            }

            const FileName fn_out = (do_split_per_micrograph ?
                // Remove any extensions from micrograph names....
                fn_root + "_" + fn_mic.withoutExtension() :
                fn_root) + fn_ext;

            // Make all necessary output directories
            if (fn_out.contains("/")) {
                const auto path = fn_out.beforeLastOf("/");
                std::string command = " mkdir -p " + path;
                system(command.c_str());
            }

            int n = 0;
            init_progress_bar(ndim);
            for (auto i : MD) {

                const FileName fn_mymic = do_split_per_micrograph ? MD.getValue<std::string>(EMDL::MICROGRAPH_NAME, i) : "";
                const int optics_group = MD.getValue<int>(EMDL::IMAGE_OPTICS_GROUP, i) - 1;
                const RFLOAT angpix = do_ignore_optics ? 1.0 : obsModel.getPixelSize(optics_group);

                if (fn_mymic == fn_mic) {

                    in.read(MD.getValue<std::string>(EMDL::IMAGE_NAME, i));

                    if (do_apply_trans || do_apply_trans_only) {
                        const RFLOAT ori_xoff = MD.getValue<RFLOAT>(EMDL::ORIENT_ORIGIN_X_ANGSTROM, i) / angpix;
                        const RFLOAT ori_yoff = MD.getValue<RFLOAT>(EMDL::ORIENT_ORIGIN_Y_ANGSTROM, i) / angpix;
                        const RFLOAT ori_psi  = MD.getValue<RFLOAT>(EMDL::ORIENT_PSI, i);

                        RFLOAT xoff, yoff, psi;
                        if (do_apply_trans_only) {
                            xoff = round(ori_xoff);
                            yoff = round(ori_yoff);
                            psi  = 0.0;
                        } else {
                            xoff = ori_xoff;
                            yoff = ori_yoff;
                            psi  = ori_psi;
                        }

                        // Apply the actual transformation
                        Matrix<RFLOAT> A = rotation2DMatrix(psi);
                        A.at(0, 2) = xoff * cos(radians(psi)) - yoff * sin(radians(psi));
                        A.at(1, 2) = yoff * cos(radians(psi)) + xoff * sin(radians(psi));
                        in() = applyGeometry(in(), A, IS_NOT_INV, DONT_WRAP);

                        MD.setValue(EMDL::ORIENT_ORIGIN_X_ANGSTROM, (ori_xoff - xoff) * angpix, i);
                        MD.setValue(EMDL::ORIENT_ORIGIN_Y_ANGSTROM, (ori_yoff - yoff) * angpix, i);
                        MD.setValue(EMDL::ORIENT_PSI, ori_psi - psi, i);
                    }
                    const auto fn_img = FileName::compose(n + 1, fn_out);
                    MD.setValue(EMDL::IMAGE_NAME, fn_img, i);

                    if (do_one_by_one) {
                        if (n == 0) {
                            in.write(fn_img, -1, false, WRITE_OVERWRITE);
                        } else {
                            in.write(fn_img, -1, true, WRITE_APPEND);
                        }
                    } else {
                        out().printShape();
                        in().printShape();
                        out().setImage(n, in());
                    }

                    if (++n % 100 == 0) progress_bar(n);
                }
            }
            progress_bar(ndim);

            if (!do_one_by_one)
                out.write(fn_out);
            std::cout << "Written out: " << fn_out << std::endl;
        }

        const auto fn_star = fn_root + ".star";
        if (do_ignore_optics) MD.write(fn_star);
        else obsModel.save(MD, fn_star, "particles");
        std::cout << "Written out: " << fn_star << std::endl;
        std::cout << "Done!" <<std::endl;
    }
};


int main(int argc, char *argv[]) {
    stack_create_parameters prm;

    try {
        prm.read(argc, argv);
        prm.run();
    } catch (RelionError XE) {
        std::cerr << XE;
        // prm.usage();
        return RELION_EXIT_FAILURE;
    }
    return RELION_EXIT_SUCCESS;
}
