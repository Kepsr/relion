#include "star_converter.h"

void StarConverter::convert_3p0_particlesTo_3p1(
    const MetaDataTable &in, MetaDataTable &outParticles, MetaDataTable &outOptics,
    std::string tablename, bool do_die_upon_error
) {
    int ver = in.getVersion();
    int curVer = MetaDataTable::getCurrentVersion();

    if (ver == curVer) {
        if (do_die_upon_error) {
            REPORT_ERROR_STR(
                "StarConverter::convert_3p0_particlesTo_3p1: Star file is already at version "
                    << curVer / 10000.0
            );
        } else {
            return;
        }
    } else if (ver > curVer) {
        if (do_die_upon_error) {
            REPORT_ERROR_STR(
                "StarConverter::convert_3p0_particlesTo_3p1: Star file is at version "
                << ver / 10000.0 << " - this is beyond the current version of Relion ("
                << curVer / 10000.0 << ")\n"
                << "You are either using an outdated copy of Relion, or the file is from the future.\n");
        } else {
            return;
        }
    }

    const int particleCount = in.numberOfObjects();

    std::vector<EMDL::EMDLabel> allOpticsLabels_double(0);

    allOpticsLabels_double.push_back(EMDL::CTF_Q0);
    allOpticsLabels_double.push_back(EMDL::IMAGE_BEAMTILT_X);
    allOpticsLabels_double.push_back(EMDL::IMAGE_BEAMTILT_Y);
    allOpticsLabels_double.push_back(EMDL::CTF_CS);
    allOpticsLabels_double.push_back(EMDL::CTF_VOLTAGE);

    allOpticsLabels_double.push_back(EMDL::CTF_DETECTOR_PIXEL_SIZE);
    allOpticsLabels_double.push_back(EMDL::CTF_MAGNIFICATION);

    std::vector<EMDL::EMDLabel> opticsLabels_double(0);

    for (int l = 0; l < allOpticsLabels_double.size(); l++) {
        if (in.labelExists(allOpticsLabels_double[l])) {
            opticsLabels_double.push_back(allOpticsLabels_double[l]);
        }
    }

    const int opticsLabelCount_double = opticsLabels_double.size();

    std::vector<std::vector<double>> groupValues_double(0);
    std::vector<int> opticsClasses(particleCount, -1);

    for (long int p = 0; p < particleCount; p++) {
        int foundGroup = -1;

        std::vector<double> curVals_double(opticsLabelCount_double);

        for (int l = 0; l < opticsLabelCount_double; l++) {
            curVals_double[l] = in.getValue<double>(opticsLabels_double[l], p);
        }

        for (int g = 0; g < groupValues_double.size(); g++) {
            bool groupGood = true;

            for (int l = 0; l < opticsLabelCount_double; l++) {
                if (curVals_double[l] != groupValues_double[g][l]) {
                    groupGood = false;
                    break;
                }
            }

            if (groupGood) {
                foundGroup = g;
                break;
            }
        }

        if (foundGroup >= 0) {
            opticsClasses[p] = foundGroup;
        } else {
            groupValues_double.push_back(curVals_double);
            opticsClasses[p] = groupValues_double.size() - 1;
        }
    }

    outParticles = in;

    for (EMDL::EMDLabel label : opticsLabels_double)
        outParticles.deactivateLabel(label);

    outParticles.addLabel(EMDL::IMAGE_OPTICS_GROUP);

    for (long int p = 0; p < particleCount; p++) {
        outParticles.setValue(EMDL::IMAGE_OPTICS_GROUP, opticsClasses[p] + 1, p);
    }

    // Determine the data type
    if (tablename.empty()) {
        tablename = in.containsLabel(EMDL::IMAGE_NAME) ? "particles" :
                    in.containsLabel(EMDL::MICROGRAPH_METADATA_NAME) ? "movies" :
                    "micrographs";
    }

    outParticles.setName(tablename);
    outParticles.setVersion(curVer);

    outOptics.setName("optics");
    outOptics.setVersion(curVer);
    outOptics.addLabel(EMDL::IMAGE_OPTICS_GROUP);
    outOptics.addLabel(EMDL::IMAGE_OPTICS_GROUP_NAME);

    for (EMDL::EMDLabel label : opticsLabels_double) {
        outOptics.addLabel(label);
    }

    for (int g = 0; g < groupValues_double.size(); g++) {
        outOptics.addObject();
        outOptics.setValue(EMDL::IMAGE_OPTICS_GROUP, g + 1, g);
        std::string mygroupname = "opticsGroup" + integerToString(g + 1);
        outOptics.setValue(EMDL::IMAGE_OPTICS_GROUP_NAME, mygroupname, g);

        for (int l = 0; l < opticsLabelCount_double; l++) {
            outOptics.setValue(opticsLabels_double[l], groupValues_double[g][l], g);
        }
    }

    // set IMAGE_PIXEL_SIZE/MICROGRAPH_PIXEL_SIZE instead of DETECTOR_PIXEL_SIZE and MAGNIFICATION
    // This does not do anything if DETECTOR_PIXEL_SIZE or MAGNIFICATION are not in the input STAR file
    unifyPixelSize(outOptics, tablename);

    if (tablename == "particles" || tablename.empty()) {
        // Make translations in Angstroms instead of in pixels
        translateOffsets(outParticles, outOptics);

        // Also read in one image for each optics group to set the image sizes in the outOptics table
        // Also set the image_size for each optics_group
        int nr_optics_groups_found = 0;
        int nr_optics_groups = groupValues_double.size();
        std::vector<bool> found_this_group;
        found_this_group.resize(nr_optics_groups, false);

        for (long int p = 0; p < particleCount; p++) {
            int g = opticsClasses[p];

            if (!found_this_group[g]) {

                FileName fn_img;
                try {
                    fn_img = outParticles.getValue<FileName>(EMDL::IMAGE_NAME, p);
                } catch (const char *errmsg) {
                    if (do_die_upon_error) {
                        REPORT_ERROR("BUG: cannot find name for particle...");
                    } else {
                        return;
                    }
                }

                try {
                    Image<double> img;
                    img.read(fn_img, false);  // false means read only header, skip real data
                    int image_size = img().xdim;

                    if (image_size % 2 != 0) {
                        REPORT_ERROR("ERROR: this program only works with even values for the image dimensions!");
                    }

                    if (image_size != img().ydim) {
                        REPORT_ERROR("ERROR: xsize != ysize: only squared images are allowed");
                    }

                    outOptics.setValue(EMDL::IMAGE_SIZE, image_size, g);
                    found_this_group[g] = true;
                    nr_optics_groups_found++;

                    if (img().zdim > 1) {
                        if (image_size != img().zdim) {
                            REPORT_ERROR("ERROR: xsize != zsize: only cube 3D images allowed");
                        }
                        outOptics.setValue(EMDL::IMAGE_DIMENSIONALITY, 3, g);
                    } else {
                        outOptics.setValue(EMDL::IMAGE_DIMENSIONALITY, 2, g);
                    }
                } catch (RelionError e) {
                    std::cerr << "Warning: " << fn_img << " not found.\n";
                    break;
                }
            }

            if (nr_optics_groups_found == nr_optics_groups) {
                break;
            }
        }

        if (nr_optics_groups_found != nr_optics_groups) {
            std::cerr << "Warning: Not all image files could be found.\n";
            std::cerr << "  Image sizes and dimensionalities will be missing from the star file.\n";
            std::cerr << "  Later steps (e.g. re-extraction, CtfRefine) can fail!!\n";
            std::cerr << "  Repeat this job after fixing the image paths.\n";

            //REPORT_ERROR("BUG: something went wrong with finding the optics groups...");
        }
    }

    return;
}

void StarConverter::unifyPixelSize(
    MetaDataTable& outOptics, const std::string &tablename
) {
    if (
        outOptics.containsLabel(EMDL::CTF_DETECTOR_PIXEL_SIZE) && 
        outOptics.containsLabel(EMDL::CTF_MAGNIFICATION)
    ) {
        for (int i = 0; i < outOptics.numberOfObjects(); i++) {

            double dstep = outOptics.getValue<double>(EMDL::CTF_DETECTOR_PIXEL_SIZE, i);
            double mag   = outOptics.getValue<double>(EMDL::CTF_MAGNIFICATION, i);

            double angpix = 10000 * dstep / mag;

            if (tablename == "particles" || tablename.empty()) {
                outOptics.setValue(EMDL::IMAGE_PIXEL_SIZE, angpix, i);
                // Do not set EMDL::MICROGRAPH_ORIGINAL_PIXEL_SIZE, because particles might have been down-sampled.
            } else if (tablename == "micrographs") {
                outOptics.setValue(EMDL::MICROGRAPH_PIXEL_SIZE,          angpix, i);
                outOptics.setValue(EMDL::MICROGRAPH_ORIGINAL_PIXEL_SIZE, angpix, i);
            }
        }

        outOptics.deactivateLabel(EMDL::CTF_DETECTOR_PIXEL_SIZE);
        outOptics.deactivateLabel(EMDL::CTF_MAGNIFICATION);
    }
}

void StarConverter::translateOffsets(
    MetaDataTable &outParticles, const MetaDataTable &optics
) {
    for (int i = 0; i < outParticles.numberOfObjects(); i++) {

        int og = outParticles.getValue<int>(EMDL::IMAGE_OPTICS_GROUP, i) - 1;
        double angpix = optics.getValue<double>(EMDL::IMAGE_PIXEL_SIZE, og);

        if (outParticles.containsLabel(EMDL::ORIENT_ORIGIN_X)) {
            double x = outParticles.getValue<double>(EMDL::ORIENT_ORIGIN_X, i);
            outParticles.setValue(EMDL::ORIENT_ORIGIN_X_ANGSTROM, x * angpix, i);
        }

        if (outParticles.containsLabel(EMDL::ORIENT_ORIGIN_Y)) {
            double y = outParticles.getValue<double>(EMDL::ORIENT_ORIGIN_Y, i);
            outParticles.setValue(EMDL::ORIENT_ORIGIN_Y_ANGSTROM, y * angpix, i);
        }

        if (outParticles.containsLabel(EMDL::ORIENT_ORIGIN_Z)) {
            double z = outParticles.getValue<double>(EMDL::ORIENT_ORIGIN_Z, i);
            outParticles.setValue(EMDL::ORIENT_ORIGIN_Z_ANGSTROM, z * angpix, i);
        }

        if (outParticles.containsLabel(EMDL::ORIENT_ORIGIN_X_PRIOR)) {
            double x = outParticles.getValue<double>(EMDL::ORIENT_ORIGIN_X_PRIOR, i);
            outParticles.setValue(EMDL::ORIENT_ORIGIN_X_PRIOR_ANGSTROM, x * angpix, i);
        }

        if (outParticles.containsLabel(EMDL::ORIENT_ORIGIN_Y_PRIOR)) {
            double y = outParticles.getValue<double> (EMDL::ORIENT_ORIGIN_Y_PRIOR, i);
            outParticles.setValue(EMDL::ORIENT_ORIGIN_Y_PRIOR_ANGSTROM, y * angpix, i);
        }

        if (outParticles.containsLabel(EMDL::ORIENT_ORIGIN_Z_PRIOR)) {
            double z = outParticles.getValue<double>(EMDL::ORIENT_ORIGIN_Z_PRIOR, i);
            outParticles.setValue(EMDL::ORIENT_ORIGIN_Z_PRIOR_ANGSTROM, z * angpix, i);
        }

        if (outParticles.containsLabel(EMDL::PARTICLE_HELICAL_TRACK_LENGTH)) {
            double d = outParticles.getValue<double>(EMDL::PARTICLE_HELICAL_TRACK_LENGTH, i);
            outParticles.setValue(EMDL::PARTICLE_HELICAL_TRACK_LENGTH_ANGSTROM, d * angpix, i);
        }
    }

    outParticles.deactivateLabel(EMDL::ORIENT_ORIGIN_X);
    outParticles.deactivateLabel(EMDL::ORIENT_ORIGIN_Y);
    outParticles.deactivateLabel(EMDL::ORIENT_ORIGIN_Z);
    outParticles.deactivateLabel(EMDL::ORIENT_ORIGIN_X_PRIOR);
    outParticles.deactivateLabel(EMDL::ORIENT_ORIGIN_Y_PRIOR);
    outParticles.deactivateLabel(EMDL::ORIENT_ORIGIN_Z_PRIOR);
    outParticles.deactivateLabel(EMDL::PARTICLE_HELICAL_TRACK_LENGTH);
}
