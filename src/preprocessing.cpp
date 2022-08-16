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
#include "src/preprocessing.h"
#include "src/jaz/ctf_helper.h"


// #define PREP_TIMING
#ifdef PREP_TIMING
Timer timer;
int TIMING_TOP              = timer.setNew("extractParticlesFromFieldOfView");
int TIMING_READ_COORD       = timer.setNew("readInCoordinateFile");
int TIMING_BIAS_CORRECT     = timer.setNew("biasCorrect");
int TIMING_EXTCT_FROM_FRAME = timer.setNew("extractParticlesFromOneFrame");
int TIMING_READ_IMG         = timer.setNew("-readImg");
int TIMING_WINDOW           = timer.setNew("-window");
int TIMING_BOUNDARY         = timer.setNew("-checkBoundary");
int TIMING_PRE_IMG_OPS      = timer.setNew("-performPerImageOperations");
int TIMING_NORMALIZE        = timer.setNew("--normalize");
int TIMING_INV_CONT         = timer.setNew("--invert_contrast");
int TIMING_COMP_STATS       = timer.setNew("--computeStats");
int TIMING_PER_IMG_OP_WRITE = timer.setNew("--write");
int TIMING_REST             = timer.setNew("-rest");
#define ifdefPREP_TIMING(statement) statement
#else
#define ifdefPREP_TIMING(statement)
#endif


void Preprocessing::read(int argc, char **argv, int rank) {
    parser.setCommandLine(argc, argv);
    int gen_section = parser.addSection("General options");
    fn_star_in = parser.getOption("--i", "The STAR file with all (selected) micrographs to extract particles from","");
    fn_coord_suffix = parser.getOption("--coord_suffix", "The suffix for the coordinate files, e.g. \"_picked.star\" or \".box\"","");
    fn_coord_dir = parser.getOption("--coord_dir", "The directory where the coordinate files are (default is same as micrographs)", "ASINPUT");
    fn_part_dir = parser.getOption("--part_dir", "Output directory for particle stacks", "Particles/");
    fn_part_star = parser.getOption("--part_star", "Output STAR file with all particles metadata", "");
    fn_data = parser.getOption("--reextract_data_star", "A _data.star file from a refinement to re-extract, e.g. with different binning or re-centered (instead of --coord_suffix)", "");
    keep_ctf_from_micrographs  = parser.checkOption("--keep_ctfs_micrographs", "By default, CTFs from fn_data will be kept. Use this flag to keep CTFs from input micrographs STAR file");
    do_reset_offsets = parser.checkOption("--reset_offsets", "reset the origin offsets from the input _data.star file to zero?");
    do_recenter = parser.checkOption("--recenter", "Re-center particle according to rlnOriginX/Y in --reextract_data_star STAR file");
    recenter_x = textToFloat(parser.getOption("--recenter_x", "X-coordinate (in pixel inside the reference) to recenter re-extracted data on", "0."));
    recenter_y = textToFloat(parser.getOption("--recenter_y", "Y-coordinate (in pixel inside the reference) to recenter re-extracted data on", "0."));
    recenter_z = textToFloat(parser.getOption("--recenter_z", "Z-coordinate (in pixel inside the reference) to recenter re-extracted data on", "0."));
    ref_angpix = textToFloat(parser.getOption("--ref_angpix", "Pixel size of the reference used for recentering. -1 uses the pixel size of particles.", "-1"));

    int extract_section = parser.addSection("Particle extraction");
    do_extract = parser.checkOption("--extract", "Extract all particles from the micrographs");
    extract_size = textToInteger(parser.getOption("--extract_size", "Size of the box to extract the particles in (in pixels)", "-1"));
    do_premultiply_ctf = parser.checkOption("--premultiply_ctf", "Premultiply the micrograph/frame with its CTF prior to particle extraction");
    premultiply_ctf_extract_size = textToInteger(parser.getOption("--premultiply_extract_size", "Size of the box to extract the particles in (in pixels) before CTF premultiplication", "-1"));
    if (premultiply_ctf_extract_size < 0)
        premultiply_ctf_extract_size = extract_size;
    do_ctf_intact_first_peak = parser.checkOption("--ctf_intact_first_peak", "When premultiplying with the CTF, leave frequencies intact until the first peak");
    do_phase_flip = parser.checkOption("--phase_flip", "Flip CTF-phases in the micrograph/frame prior to particle extraction");
    extract_bias_x  = textToInteger(parser.getOption("--extract_bias_x", "Bias in X-direction of picked particles (this value in pixels will be added to the coords)", "0"));
    extract_bias_y  = textToInteger(parser.getOption("--extract_bias_y", "Bias in Y-direction of picked particles (this value in pixels will be added to the coords)", "0"));
    only_extract_unfinished = parser.checkOption("--only_do_unfinished", "Extract only particles if the STAR file for that micrograph does not yet exist.");

    int perpart_section = parser.addSection("Particle operations");
    do_project_3d = parser.checkOption("--project3d", "Project sub-tomograms along Z to generate 2D particles");
    scale  = textToInteger(parser.getOption("--scale", "Re-scale the particles to this size (in pixels)", "-1"));
    window  = textToInteger(parser.getOption("--window", "Re-window the particles to this size (in pixels)", "-1"));
    do_normalise = parser.checkOption("--norm", "Normalise the background to average zero and stddev one");
    do_ramp = !parser.checkOption("--no_ramp", "Just subtract the background mean in the normalisation, instead of subtracting a fitted ramping background. ");
    bg_radius = textToInteger(parser.getOption("--bg_radius", "Radius of the circular mask that will be used to define the background area (in pixels)", "-1"));
    white_dust_stddev = textToFloat(parser.getOption("--white_dust", "Sigma-values above which white dust will be removed (negative value means no dust removal)","-1"));
    black_dust_stddev = textToFloat(parser.getOption("--black_dust", "Sigma-values above which black dust will be removed (negative value means no dust removal)","-1"));
    do_invert_contrast = parser.checkOption("--invert_contrast", "Invert the contrast in the input images");
    fn_operate_in = parser.getOption("--operate_on", "Use this option to operate on an input image stack ", "");
    fn_operate_out = parser.getOption("--operate_out", "Output name when operating on an input image stack", "preprocessed.mrcs");

    int helix_section = parser.addSection("Helix extraction");
    do_extract_helix = parser.checkOption("--helix", "Extract helical segments");
    helical_tube_outer_diameter = textToFloat(parser.getOption("--helical_outer_diameter", "Outer diameter of helical tubes in Angstroms (for masks of helical segments)", "-1."));
    do_extract_helical_tubes = parser.checkOption("--helical_tubes", "Extract helical segments from tube coordinates");
    helical_nr_asu = textToInteger(parser.getOption("--helical_nr_asu", "Number of helical asymmetrical units", "1"));
    helical_rise = textToFloat(parser.getOption("--helical_rise", "Helical rise (in Angstroms)", "0."));
    helical_bimodal_angular_priors = parser.checkOption("--helical_bimodal_angular_priors", "Add bimodal angular priors for helical segments");
    helical_cut_into_segments = parser.checkOption("--helical_cut_into_segments", "Cut helical tubes into segments");
    // Initialise verb for non-parallel execution
    verb = 1;
}

void Preprocessing::usage() {
    parser.writeUsage(std::cout);
}

void Preprocessing::initialise() {
    // Check for errors in the command-line option
    if (parser.checkForErrors())
        REPORT_ERROR("Errors encountered on the command line (see above), exiting...");

    if (!do_extract && fn_operate_in == "")
        REPORT_ERROR("Provide either --extract or --operate_on");

    // Make sure the output directory name ends with a '/'
    if (fn_part_dir[fn_part_dir.length() - 1] != '/')
        fn_part_dir += "/";

    // Set up which coordinate files to extract particles from (or to join STAR file for)
    if (do_extract) {
        if (verb > 0) {
            if (fn_star_in=="" || (fn_coord_suffix=="" && fn_data == ""))
                REPORT_ERROR("Preprocessing::initialise ERROR: please provide --i and either (--coord_suffix or --reextract_data_star) to extract particles");

            if (extract_size < 0)
                REPORT_ERROR("Preprocessing::initialise ERROR: please provide the size of the box to extract particle using --extract_size ");

            if (extract_size % 2 != 0)
                REPORT_ERROR("Preprocessing::initialise ERROR: only extracting to even-sized images is allowed in RELION...");
        }

        // Read in the micrographs STAR file
        ObservationModel::loadSafely(fn_star_in, obsModelMic, MDmics, "micrographs", verb);

        if (!MDmics.containsLabel(EMDL::MICROGRAPH_NAME))
            REPORT_ERROR("Preprocessing::initialise ERROR: Input micrograph STAR file has no rlnMicrographName column!");

        if (!obsModelMic.opticsMdt.containsLabel(EMDL::MICROGRAPH_PIXEL_SIZE))
            REPORT_ERROR("Preprocessing::initialise ERROR: Input micrograph STAR file has no rlnMicrographPixelSize column!");

        // Get the dimensionality of the micrographs

        // Read the header of the micrograph to see how many frames there are.
        Image<RFLOAT> Imic;

        MDmics.goToObject(0);
        FileName fn_mic = MDmics.getValue<std::string>(EMDL::MICROGRAPH_NAME);
        Imic.read(fn_mic, false, -1, false); // readData = false, select_image = -1, mapData= false, is_2D = true);
        Image<RFLOAT>::Dimensions dimensions = Imic.getDimensions();
             int xdim = dimensions.x;
             int ydim = dimensions.y;
             int zdim = dimensions.z;
        long int ndim = dimensions.n;
        dimensionality = 2 + (zdim > 1);

        mic_star_has_ctf = MDmics.containsLabel(EMDL::CTF_DEFOCUSU);

        if ((do_phase_flip || do_premultiply_ctf) && !mic_star_has_ctf)
            REPORT_ERROR("Preprocessing::initialise ERROR: No CTF information found in the input micrograph STAR-file");

        if (fn_data != "") {
            if (verb > 0) {
                std::cout << " + Re-extracting particles based on coordinates from input _data.star file " << std::endl;
                std::cout << " + " << fn_data << std::endl;
            }
            ObservationModel::loadSafely(fn_data, obsModelPart, MDimg, "particles", verb);
            data_star_has_ctf = MDimg.containsLabel(EMDL::CTF_DEFOCUSU);

            if (do_recenter && ref_angpix <= 0) {
                if (!obsModelPart.allPixelSizesIdentical())
                    REPORT_ERROR("The pixel sizes in the particle STAR file are not identical. Please specify the pixel size of the reference for re-centering as --ref_angpix.");
            }

            if (do_recenter && verb > 0) {
                std::cout << " + And re-centering particles based on refined coordinates in the _data.star file." << std::endl;
                if (fabs(recenter_x) > 0.0 || fabs(recenter_y) > 0.0 || fabs(recenter_z) > 0.0) {
                    if (ref_angpix > 0) {
                        std::cout << "   This uses " << ref_angpix << " A/px to convert the recentering coordinate from pixels to Angstrom." << std::endl;
                    } else {
                        std::cout << "   This assumes the particle pixel size is the same as the reference pixel size.\n   If this is not the case, please specify --ref_angpix." << std::endl;
                    }
                }
            }

        } else {
            data_star_has_ctf = false;
            // Make sure the coordinate file directory names end with a '/'
            if (fn_coord_dir != "ASINPUT" && fn_coord_dir[fn_coord_dir.length() - 1] != '/')
                fn_coord_dir += "/";

            // Loop over all micrographs in the input STAR file and warn of coordinate file or micrograph file do not exist
            if (verb > 0 && fn_data == "") {
                FileName fn_mic;
                for (long int i : MDmics) {
                    fn_mic = MDmics.getValue<std::string>(EMDL::MICROGRAPH_NAME);
                    FileName fn_coord = getCoordinateFileName(fn_mic);
                    if (!exists(fn_coord))
                        std::cerr << "Warning: coordinate file " << fn_coord << " does not exist..." << std::endl;
                    if (!exists(fn_mic))
                        std::cerr << "Warning: micrograph file " << fn_mic << " does not exist..." << std::endl;
                }
            }
        }

        if (do_phase_flip || do_premultiply_ctf) {
            if (!mic_star_has_ctf && !data_star_has_ctf) {
                REPORT_ERROR("Preprocessing:: ERROR: cannot phase flip or premultiply CTF without input CTF parameters in the micrograph or data STAR file");
            }
        }
    }

    if (do_extract || fn_operate_in != "") {
        // Check whether to do re-scaling
        do_rescale = scale > 0;
        if (do_rescale && scale % 2 != 0) {
            scale++;
            std::cerr << " Warning: only re-scaling to even-sized images is allowed in RELION, setting scale to: " << scale << std::endl;
        }

        // Check whether to do re-windowing
        do_rewindow = (window > 0);
        if (do_rewindow && window % 2 != 0) {
            window++;
            std::cerr << " Warning: only re-windowing to even-sized images is allowed in RELION, setting window to: " << window << std::endl;
        }

        // Check for bg_radius in case of normalisation
        if (do_normalise && bg_radius < 0)
            REPORT_ERROR("ERROR: please provide a radius for a circle that defines the background area when normalising...");

        // Extract helical segments
        if (do_extract_helix) {
            if (!do_extract || fn_operate_in != "")
                REPORT_ERROR("ERROR: cannot perform operations on helical segments other than extraction");
            if (do_normalise && helical_tube_outer_diameter < 0)
                REPORT_ERROR("ERROR: please provide a tube radius that defines the background area when normalising helical segments...");
        }
    }
}

void Preprocessing::run() {
    if (do_extract) {
        runExtractParticles();
    } else if (fn_operate_in != "") {
        runOperateOnInputFile();
    }

    if (verb > 0)
        std::cout << " Done preprocessing!" << std::endl;

    #ifdef PREP_TIMING
    timer.printTimes(false);
    #endif
}

void Preprocessing::joinAllStarFiles() {
    FileName fn_ostar;
    int og;
    std::cout << " Joining metadata of all particles from " << MDmics.numberOfObjects() << " micrographs in one STAR file..." << std::endl;

    long int imic = 0, ibatch = 0;
    MetaDataTable MDout, MDmicnames, MDbatch;
    for (long int current_object1 : MDmics) {
        // Micrograph filename
        FileName fn_mic = MDmics.getValue<std::string>(EMDL::MICROGRAPH_NAME);

        // Get the filename of the STAR file for just this micrograph
        FileName fn_star = getOutputFileNameRoot(fn_mic) + "_extract.star";

        if (fn_part_star != "") {
            if (exists(fn_star)) {
                MetaDataTable MDonestack;
                MDonestack.read(fn_star);

                if (MDout.numberOfObjects() > 0 && !MetaDataTable::compareLabels(MDout, MDonestack)) {
                    std::cout << "The STAR file " << fn_star << " contains a column not present in others. Missing values will be filled by default values (0 or empty string)" << std::endl;
                    MDout.addMissingLabels(MDonestack);
                    MDonestack.addMissingLabels(MDout);
                }
                MDout.append(MDonestack);
            }
        }

        imic++;
    }

    std::string optgroup_name;
    // Write out the joined star files
    if (fn_part_star != "") {
        // Get pixel size in original micrograph from obsModelMic, as this may no longer be present in obsModelPart
        std::map<std::string, RFLOAT> optics_group_mic_angpix;
        if (fn_data != "") {
            RFLOAT mic_angpix;
            for (long int i : obsModelMic.opticsMdt) {
                mic_angpix    = obsModelMic.opticsMdt.getValue<RFLOAT>     (EMDL::MICROGRAPH_PIXEL_SIZE);
                optgroup_name = obsModelMic.opticsMdt.getValue<std::string>(EMDL::IMAGE_OPTICS_GROUP_NAME);

                optics_group_mic_angpix.insert(std::make_pair(optgroup_name, mic_angpix));
            }
        }

        ObservationModel *myOutObsModel;
        myOutObsModel = fn_data == "" || keep_ctf_from_micrographs ? &obsModelMic : &obsModelPart;

        std::set<std::string> isOgPresent;
        for (long int i : MDout) {
            og = myOutObsModel->getOpticsGroup(MDout);
            std::string optgroup_name = myOutObsModel->opticsMdt.getValue<std::string>(EMDL::IMAGE_OPTICS_GROUP_NAME, og);
            isOgPresent.insert(optgroup_name);
        }

        RFLOAT my_angpix;
        // Set the (possibly rescale output_angpix and the output image size in the opticsMdt
        for (long int i : myOutObsModel->opticsMdt) {
            // Find the pixel size for the original micrograph
            if (fn_data.empty()) {
                my_angpix = obsModelMic.opticsMdt.getValue<RFLOAT>(EMDL::MICROGRAPH_PIXEL_SIZE);
            } else {
                std::string optgroup_name = myOutObsModel->opticsMdt.getValue<std::string>(EMDL::IMAGE_OPTICS_GROUP_NAME);
                if (optics_group_mic_angpix.count(optgroup_name) == 0) {
                    if (isOgPresent.count(optgroup_name) != 0) {
                        REPORT_ERROR("ERROR: optics group \"" + optgroup_name + "\" does not exist in micrograph STAR file...");
                    } else {
                        my_angpix = -1; // mark for deletion
                    }
                } else {
                    my_angpix = optics_group_mic_angpix[optgroup_name];
                }
            }

            if (do_rescale)
                my_angpix *= (RFLOAT)extract_size / (RFLOAT)scale;
            myOutObsModel->opticsMdt.setValue(EMDL::IMAGE_PIXEL_SIZE, my_angpix);

            if (do_rewindow) {
                myOutObsModel->opticsMdt.setValue(EMDL::IMAGE_SIZE, window);
            } else if (do_rescale) {
                myOutObsModel->opticsMdt.setValue(EMDL::IMAGE_SIZE, scale);
            } else {
                myOutObsModel->opticsMdt.setValue(EMDL::IMAGE_SIZE, extract_size);
            }

            myOutObsModel->opticsMdt.setValue(EMDL::IMAGE_DIMENSIONALITY, dimensionality);

            if (do_premultiply_ctf)
                myOutObsModel->opticsMdt.setValue(EMDL::OPTIMISER_DATA_ARE_CTF_PREMULTIPLIED, true);

            int igroup = myOutObsModel->opticsMdt.getValue<int>(EMDL::IMAGE_OPTICS_GROUP);
            if (my_angpix < 0) {
                std::cerr << "optics group \"" + optgroup_name + "\" will be removed because no extracted particle belong to it." << std::endl;
            } else {
                std::cout << " The pixel size of the extracted particles in optics group " << igroup << " is " << my_angpix << " Angstrom/pixel." << std::endl;
            }
        }

        if (fn_data.empty()) {
            myOutObsModel->opticsMdt.deactivateLabel(EMDL::MICROGRAPH_PIXEL_SIZE);
        }

        // Remove absent optics groups; After this, NOTHING should be done except for saving. obsModel's internal data structure is now corrupted!
        og = 0;
        while (og < myOutObsModel->opticsMdt.numberOfObjects()) {
            RFLOAT og_angpix = myOutObsModel->opticsMdt.getValue<RFLOAT>(EMDL::IMAGE_PIXEL_SIZE, og);
            if (og_angpix < 0) {
                myOutObsModel->opticsMdt.removeObject(og);
            } else {
                og++;
            }
        }

        ObservationModel::saveNew(MDout, myOutObsModel->opticsMdt, fn_part_star, "particles");
        std::cout << " Written out STAR file with " << MDout.numberOfObjects() << " particles in " << fn_part_star<< std::endl;
    }
}

void Preprocessing::runExtractParticles() {
    long int nr_mics = MDmics.numberOfObjects();

    int barstep;
    if (verb > 0) {
        std::cout << " Extracting particles from " << nr_mics << " micrographs ..." << std::endl;
        init_progress_bar(nr_mics);
        barstep = std::max(1, (int) nr_mics / 60);
    }
    MetaDataTable MDoutMics;  // during re-extraction we may not always use particles from all mics.
    FileName fn_mic, fn_olddir = "";
    long int imic = 0;
    bool micIsUsed;
    for (long int i : MDmics) {
        // Abort through the pipeline_control system
        if (pipeline_control_check_abort_job())
            exit(RELION_EXIT_ABORTED);

        fn_mic = MDmics.getValue<std::string>(EMDL::MICROGRAPH_NAME);
        int optics_group = obsModelMic.getOpticsGroup(MDmics);

        // Set the pixel size for this micrograph
        angpix = obsModelMic.getPixelSize(optics_group);
        // Also set the output_angpix (which could be rescaled)
        output_angpix = angpix;
        if (do_rescale) output_angpix *= (RFLOAT) extract_size / (RFLOAT) scale;

        // Check new-style outputdirectory exists and make it if not!
        FileName fn_dir = getOutputFileNameRoot(fn_mic);
        fn_dir = fn_dir.beforeLastOf("/");
        if (fn_dir != fn_olddir && !exists(fn_dir)) {
            // Make a Particles directory
            system(("mkdir -p " + fn_dir).c_str());
            fn_olddir = fn_dir;
        }

        if (verb > 0 && imic % barstep == 0) progress_bar(imic);

        {
        ifdefPREP_TIMING(TicToc tt (timer, TIMING_TOP);)
        micIsUsed = extractParticlesFromFieldOfView(fn_mic, imic);
        }

        if (micIsUsed) MDoutMics.addObject(MDmics.getObject(i));

        imic++;
    }

    MDmics = MDoutMics;
    if (verb > 0) progress_bar(fn_coords.size());

    // Now combine all metadata in a single STAR file
    joinAllStarFiles();
}

void Preprocessing::readCoordinates(FileName fn_coord, MetaDataTable &MD) {
    MD.clear();

    bool is_star = fn_coord.getExtension() == "star";
    bool is_box  = fn_coord.getExtension() == "box";

    if (is_star) {
        MD.read(fn_coord);
    } else {
        std::ifstream in (fn_coord.data(), std::ios_base::in);
        if (in.fail())
            REPORT_ERROR((std::string) "Preprocessing::readCoordinates ERROR: File " + fn_coord + " does not exists");

        // Start reading the ifstream at the top
        in.seekg(0);
        std::string line;
        for (int n = 0; getline(in, line, '\n'); n++) {
            std::vector<std::string> words;
            tokenize(line, words);

            if (is_box) {
                if (words.size() < 4)
                    REPORT_ERROR("Preprocessing::readCoordinates Unexpected number of words on data line of " + fn_coord);
                int num1 = textToInteger(words[0]);
                int num2 = textToInteger(words[1]);
                int num3 = textToInteger(words[2]);
                int num4 = textToInteger(words[3]);
                int xpos = num1 + num3 / 2;
                int ypos = num2 + num4 / 2;

                MD.addObject();
                MD.setValue(EMDL::IMAGE_COORD_X, (RFLOAT)xpos);
                MD.setValue(EMDL::IMAGE_COORD_Y, (RFLOAT)ypos);
            } else {
                // Try reading as plain ASCII....
                // The first line might be a special header (as in ximdisp or xmipp2 format)
                // But could also be a data line (as in plain text format)
                if (n == 0) {
                    int num1, num2, num3;
                    bool is_data = false;

                    // Ignore lines that do not have at least two integer numbers on it (at this point I do not know dimensionality yet....)
                    if (
                        words.size() > 1 &&
                        sscanf(words[0].c_str(), "%d", &num1) &&
                        sscanf(words[1].c_str(), "%d", &num2)
                    ) {
                        MD.addObject();
                        MD.setValue(EMDL::IMAGE_COORD_X, (RFLOAT) num1);
                        MD.setValue(EMDL::IMAGE_COORD_Y, (RFLOAT) num2);

                        // It could also be a X,Y,Z coordinate...
                        if (words.size() > 2 && sscanf(words[2].c_str(), "%d", &num3))
                            MD.setValue(EMDL::IMAGE_COORD_Z, (RFLOAT) num3);
                    }
                } else {
                    // All other lines contain x, y as first two entries, and possibly a Z-coordinate as well.
                    // Note the Z-coordinate will only be used for 3D micrographs (i.e. tomograms), so the third column for a 2D Ximdisp format will be ignored later on
                    if (words.size() < 2)
                        REPORT_ERROR("Preprocessing::readCoordinates Unexpected number of words on data line of " + fn_coord);

                    int num1 = textToInteger(words[0]);
                    int num2 = textToInteger(words[1]);
                    int num3;

                    MD.addObject();
                    MD.setValue(EMDL::IMAGE_COORD_X, (RFLOAT) num1);
                    MD.setValue(EMDL::IMAGE_COORD_Y, (RFLOAT) num2);
                    // It could also be a X,Y,Z coordinate...
                    if (words.size() > 2 && sscanf(words[2].c_str(), "%d", &num3))
                    MD.setValue(EMDL::IMAGE_COORD_Z, (RFLOAT)num3);
                }
            }
        }
    }
}

void Preprocessing::readHelicalCoordinates(FileName fn_mic, FileName fn_coord, MetaDataTable &MD) {
    MD.clear();

    // Get micrograph name and check it exists
    if (!exists(fn_mic)) {
        std::cerr << "WARNING: cannot find micrograph file " << fn_mic << std::endl;
        return;
    }

    // Read the header of the micrograph to see X and Y dimensions
    Image<RFLOAT> Imic;
    Imic.read(fn_mic, false, -1, false); // readData = false, select_image = -1, mapData= false, is_2D = true);

    Image<RFLOAT>::Dimensions dimensions = Imic.getDimensions();

    if (dimensions.n != 1) REPORT_ERROR("Preprocessing::readCoordinates ERROR: Extraction of helical segments - " + (std::string)(fn_mic) + " is a stack, not a 2D micrograph or 3D tomogram!");

    bool is_3D     = dimensions.z > 1;
    bool is_star   = fn_coord.getExtension() == "star";
    bool is_box    = fn_coord.getExtension() == "box";
    bool is_coords = fn_coord.getExtension() == "coords";
    if (!is_star && !is_box && !is_coords)
        REPORT_ERROR("Preprocessing::readCoordinates ERROR: Extraction of helical segments - Unknown file extension (RELION *.star, EMAN2 *.box and XIMDISP *.coords are supported).");
    if (is_3D && !is_star)
        REPORT_ERROR("Preprocessing::readCoordinates ERROR: Extraction of 3D helical subtomograms - Only STAR coordinate files are supported!");

    int total_segments, total_tubes;
    if (is_star) {
        // std::cerr << " DEBUG: Extracting helical segments / subtomograms from RELION STAR coordinate files..." << std::endl;
        if (do_extract_helical_tubes) {
            if (is_3D) REPORT_ERROR("Preprocessing::readCoordinates ERROR: Cannot extract 3D helical subtomograms from start-end coordinates!");
            convertHelicalTubeCoordsToMetaDataTable(fn_coord, MD, total_segments, total_tubes, helical_nr_asu, helical_rise, angpix, dimensions.x, dimensions.y, extract_size, helical_bimodal_angular_priors, helical_cut_into_segments);
        } else {
            convertHelicalSegmentCoordsToMetaDataTable(fn_coord, MD, total_segments, is_3D, dimensions.x, dimensions.y, dimensions.z, extract_size, helical_bimodal_angular_priors);
        }
    } else if (is_box) {
        if (do_extract_helical_tubes) {
            convertEmanHelicalTubeCoordsToMetaDataTable(fn_coord, MD, total_segments, total_tubes, helical_nr_asu, helical_rise, angpix, dimensions.x, dimensions.y, extract_size, helical_bimodal_angular_priors, helical_cut_into_segments);
        } else {
            convertEmanHelicalSegmentCoordsToMetaDataTable(fn_coord, MD, total_segments, total_tubes, angpix, dimensions.x, dimensions.y, extract_size, helical_bimodal_angular_priors);
        }
    } else if (is_coords) {
        if (do_extract_helical_tubes) {
            convertXimdispHelicalTubeCoordsToMetaDataTable(fn_coord, MD, total_segments, total_tubes, helical_nr_asu, helical_rise, angpix, dimensions.x, dimensions.y, extract_size, helical_bimodal_angular_priors, helical_cut_into_segments);
        } else {
            convertXimdispHelicalSegmentCoordsToMetaDataTable(fn_coord, MD, total_segments, total_tubes, angpix, dimensions.x, dimensions.y, extract_size, helical_bimodal_angular_priors);
        }
    } else {
        REPORT_ERROR("Preprocessing::readCoordinates ERROR: Extraction of helical segments - Unknown file extension (RELION *.star, EMAN2 *.box and XIMDISP *.coords are supported).");
    }
}

bool Preprocessing::extractParticlesFromFieldOfView(FileName fn_mic, long int imic) {
    // Name of the output particle stack

    FileName fn_pre, fn_jobnr, fn_post;
    decomposePipelineFileName(fn_mic, fn_pre, fn_jobnr, fn_post);
    FileName fn_output_img_root = fn_part_dir + fn_post.withoutExtension();
    FileName fn_oristack = fn_post.withoutExtension() + "_extract.mrcs";

    // Name of this micrographs STAR file
    FileName fn_star = fn_output_img_root + "_extract.star";

    if (exists(fn_star) && only_extract_unfinished) return true;

    MetaDataTable MDin, MDout;
    {
    ifdefPREP_TIMING(TicToc tt (timer, TIMING_READ_COORD);)
    // Read in the coordinates file
    if (fn_data != "") {
        // Search for this micrograph in the MDdata table
        MDin = getCoordinateMetaDataTable(fn_mic);
    } else {
        FileName fn_coord = getCoordinateFileName(fn_mic);
        if (!exists(fn_coord)) return false;
        if (do_extract_helix) {
            readHelicalCoordinates(fn_mic, fn_coord, MDin);
        } else {
            readCoordinates(fn_coord, MDin);
        }
    }
    }

    if (MDin.numberOfObjects() > 0) {
        {
        ifdefPREP_TIMING(TicToc tt (timer, TIMING_BIAS_CORRECT);)
        // Correct for bias in the picked coordinates
        if (abs(extract_bias_x) > 0 || abs(extract_bias_y) > 0) {
            for (long int i : MDin) {
                RFLOAT xcoor = MDin.getValue<RFLOAT>(EMDL::IMAGE_COORD_X)
                             + extract_bias_x;
                RFLOAT ycoor = MDin.getValue<RFLOAT>(EMDL::IMAGE_COORD_Y)
                             + extract_bias_y;
                MDin.setValue(EMDL::IMAGE_COORD_X, xcoor);
                MDin.setValue(EMDL::IMAGE_COORD_Y, ycoor);
            }
        }
        }

        // Warn for small groups
        int npos = MDin.numberOfObjects();
        if (npos < 10) {
            std::cerr << "Warning: There are only " << npos << " particles in micrograph " << fn_mic <<". Consider joining multiple micrographs into one group. "<< std::endl;
        }

        // Get micrograph name and check it exists
        if (!exists(fn_mic)) {
            std::cerr << "WARNING: Skipping " << fn_mic << ", which has " << npos << " particles, because cannot find the file..." << std::endl;
            return false;
        }

        if (dimensionality == 3) {
            do_ramp = false;
            if (do_phase_flip || do_premultiply_ctf) {
                REPORT_ERROR(std::string(__func__) + " ERROR: cannot do CTF premultiplication or phase flipping as dimensionality is not 2!");
            }
        }

        long int my_current_nr_images = 0;
        RFLOAT all_avg = 0;
        RFLOAT all_stddev = 0;
        RFLOAT all_minval =  LARGE_NUMBER;
        RFLOAT all_maxval = -LARGE_NUMBER;

        extractParticlesFromOneMicrograph(
            MDin, fn_mic, imic, fn_output_img_root, fn_oristack,
            my_current_nr_images, npos,
            all_avg, all_stddev, all_minval, all_maxval
        );

        MDout.append(MDin);
        // Keep track of total number of images extracted thus far
        my_current_nr_images += npos;

        /// BUG: No TIC!
        ifdefPREP_TIMING(timer.toc(TIMING_EXTCT_FROM_FRAME);)

        MDout.setName("images");
        MDout.write(fn_star);
        return true;
    } else {
        std::cerr << " WARNING: no particles on micrograph: " << fn_mic << std::endl;
        return false;
    }
}

// Actually extract particles. This can be from one micrograph
void Preprocessing::extractParticlesFromOneMicrograph(MetaDataTable &MD,
    FileName fn_mic, int imic,
    FileName fn_output_img_root, FileName fn_oristack, long int &my_current_nr_images, long int my_total_nr_images,
    RFLOAT &all_avg, RFLOAT &all_stddev, RFLOAT &all_minval, RFLOAT &all_maxval
) {

    bool MDin_has_optics_group = MD.containsLabel(EMDL::IMAGE_OPTICS_GROUP); // i.e. re-extracting
    bool MDin_has_beamtilt     = MD.containsLabel(EMDL::IMAGE_BEAMTILT_X) || MD.containsLabel(EMDL::IMAGE_BEAMTILT_Y);
    bool MDin_has_ctf          = MD.containsLabel(EMDL::CTF_DEFOCUSU);
    bool MDin_has_tiltgroup    = MD.containsLabel(EMDL::PARTICLE_BEAM_TILT_CLASS);
    int my_extract_size        = do_phase_flip || do_premultiply_ctf ? premultiply_ctf_extract_size : extract_size;

    Image<RFLOAT> Imic;
    RFLOAT my_angpix, mic_avg;
    {
    ifdefPREP_TIMING(TicToc tt (timer, TIMING_READ_IMG);)

    Imic.read(fn_mic);
    // Calculate average value in the micrograph, for filling empty region around large-box extraction for premultiplication with CTF
    mic_avg = average(Imic());

    }

    CTF ctf;
    ObservationModel *obsModel;
    int optics_group;
    if (mic_star_has_ctf || keep_ctf_from_micrographs) {
        ctf = CtfHelper::makeCTF(MDmics, &obsModelMic, imic);
        obsModel = &obsModelMic;
        optics_group = obsModelMic.getOpticsGroup(MDmics, imic);

        // Micrograph STAR file might not have a box size
        obsModelMic.setBoxSize(optics_group, my_extract_size);
        my_angpix = obsModelMic.opticsMdt.getValue<RFLOAT>(EMDL::MICROGRAPH_PIXEL_SIZE, optics_group);
    }

    // Now window all particles from the micrograph
    // Now do the actual phase flipping or CTF-multiplication
    MultidimArray<Complex> FT;
    FourierTransformer transformer;
    int ipos = 0;
    for (long int _ : MD) {
        RFLOAT dxpos = MD.getValue<RFLOAT>(EMDL::IMAGE_COORD_X);
        RFLOAT dypos = MD.getValue<RFLOAT>(EMDL::IMAGE_COORD_Y);
        long int xpos = (long int) dxpos;
        long int ypos = (long int) dypos;

        long int x0 = xpos + Xmipp::init(my_extract_size);
        long int xF = xpos + Xmipp::last(my_extract_size);
        long int y0 = ypos + Xmipp::init(my_extract_size);
        long int yF = ypos + Xmipp::last(my_extract_size);

        RFLOAT dzpos;
        long int zpos;
        long int z0, zF;
        if (dimensionality == 3) {
            dzpos = MD.getValue<RFLOAT>(EMDL::IMAGE_COORD_Z);
            zpos = (long int) dzpos;
            z0 = zpos + Xmipp::init(extract_size);
            zF = zpos + Xmipp::last(extract_size);
        }

        // Discard particles that are completely outside the micrograph and print a warning
        if (
                                     xF < 0 || x0 >= Xsize(Imic()) ||
                                     yF < 0 || y0 >= Ysize(Imic()) ||
            (dimensionality == 3 && (zF < 0 || z0 >= Zsize(Imic())))
        ) {
            std::cerr << " micrograph x,y,z,n-size= " << Xsize(Imic()) << " , " << Ysize(Imic()) << " , " << Zsize(Imic()) << " , " << Nsize(Imic()) << std::endl;
            std::cerr << " particle position= " << xpos << " , " << ypos;
            if (dimensionality == 3)
                std::cerr << " , " << zpos;
            std::cerr << std::endl;
            REPORT_ERROR("Preprocessing::extractParticlesFromOneFrame ERROR: particle" + integerToString(ipos + 1) + " lies completely outside micrograph " + fn_mic);
        }

        // Read per-particle CTF
        if (MDin_has_ctf && !keep_ctf_from_micrographs) {
            ctf = CtfHelper::makeCTF(MD, &obsModelPart, optics_group);
            obsModel = &obsModelPart;
            optics_group = obsModelPart.getOpticsGroup(MD);
            if (obsModelPart.getBoxSize(optics_group) != my_extract_size)
                obsModelPart.setBoxSize(optics_group, my_extract_size);
            my_angpix = obsModelPart.opticsMdt.getValue<RFLOAT>(EMDL::MICROGRAPH_PIXEL_SIZE, optics_group);
        }

        Image<RFLOAT> Ipart;
        {
        ifdefPREP_TIMING(TicToc tt (timer, TIMING_WINDOW);)
        // extract one particle in Ipart
        if (dimensionality == 3) {
            Imic().window(Ipart(), z0, y0, x0, zF, yF, xF);
        } else {
            Imic().window(Ipart(), y0, x0, yF, xF, mic_avg);
        }
        Ipart().setXmippOrigin();
        }

        // Premultiply the CTF of each particle, possibly in a bigger box (premultiply_ctf_extract_size)
        if (do_phase_flip || do_premultiply_ctf) {
            FT = transformer.FourierTransform(Ipart());  // std::move?

            // 190802 TAKANORI: The original code using getCTF was do_damping=false, but for consistency with Polishing, I changed it.
            // The boxsize in ObsModel has been updated above.
            // In contrast to Polish, we premultiply particle BEFORE down-sampling, so PixelSize in ObsModel is OK.
            // But we are doing this after extraction, so there is not much merit...
            MultidimArray<RFLOAT> Fctf = CtfHelper::getFftwImage(
                ctf,
                Xsize(FT), Ysize(FT),
                my_extract_size, my_extract_size, my_angpix,
                obsModel,
                false, do_phase_flip, do_ctf_intact_first_peak, true, false
                // do_abs, phase_flip, intact_first_peak, damping, padding
            );

            FT *= Fctf;

            Ipart() = transformer.inverseFourierTransform(FT);

            if (extract_size != premultiply_ctf_extract_size) {
                Ipart().window(
                    Xmipp::init(extract_size), Xmipp::init(extract_size),
                    Xmipp::last(extract_size), Xmipp::last(extract_size)
                );
            }
        }

        {
        ifdefPREP_TIMING(TicToc tt (timer, TIMING_BOUNDARY);)
        // Check boundaries: fill pixels outside the boundary with the nearest ones inside
        // This will create lines at the edges, rather than zeros
        Ipart().setXmippOrigin();

        // X-boundaries
        if (x0 < 0 || xF >= Xsize(Imic())) {
            FOR_ALL_ELEMENTS_IN_ARRAY3D(Ipart()) {
                if (i + xpos < 0) {
                    Ipart().elem(i, j, k) = Ipart().elem(-xpos, j, k);
                } else if (i + xpos >= Xsize(Imic())) {
                    Ipart().elem(i, j, k) = Ipart().elem(Xsize(Imic()) - xpos - 1, j, k);
                }
            }
        }

        // Y-boundaries
        if (y0 < 0 || yF >= Ysize(Imic())) {
            FOR_ALL_ELEMENTS_IN_ARRAY3D(Ipart()) {
                if (j + ypos < 0) {
                    Ipart().elem(i, j, k) = Ipart().elem(i, -ypos, k);
                } else if (j + ypos >= Ysize(Imic())) {
                    Ipart().elem(i, j, k) = Ipart().elem(i, Ysize(Imic()) - ypos - 1, k);
                }
            }
        }

        if (dimensionality == 3) {
            // Z-boundaries
            if (z0 < 0 || zF >= Zsize(Imic())) {
                FOR_ALL_ELEMENTS_IN_ARRAY3D(Ipart()) {
                    if (k + zpos < 0) {
                        Ipart().elem(i, j, k) = Ipart().elem(i, j, -zpos);
                    } else if (k + zpos >= Zsize(Imic())) {
                        Ipart().elem(i, j, k) = Ipart().elem(i, j, Zsize(Imic()) - zpos - 1);
                    }
                }
            }
        }

        // 2D projection of 3D sub-tomograms
        if (dimensionality == 3 && do_project_3d) {
            // Project the 3D sub-tomogram into a 2D particle again
            Image<RFLOAT> Iproj(Ysize(Ipart()), Xsize(Ipart()));
            Iproj().setXmippOrigin();
            for (long int k = 0; k < Zsize(Ipart()); k++)
            for (long int j = 0; j < Ysize(Ipart()); j++)
            for (long int i = 0; i < Xsize(Ipart()); i++) {
                direct::elem(Iproj(), i, j) += direct::elem(Ipart(), i, j, k);
            }
            Ipart = Iproj;
        }
        }

        // performPerImageOperations will also append the particle to the output stack in fn_stack
        // Jun24,2015 - Shaoda, extract helical segments
        RFLOAT tilt_deg, psi_deg;
        tilt_deg = psi_deg = 0.0;
        if (do_extract_helix) {
            // If priors do not exist, errors will occur in 'readHelicalCoordinates()'.
            tilt_deg = MD.getValue<RFLOAT>(EMDL::ORIENT_TILT_PRIOR);
            psi_deg  = MD.getValue<RFLOAT>(EMDL::ORIENT_PSI_PRIOR);
        }

        {
        ifdefPREP_TIMING(TicToc tt (timer, TIMING_PRE_IMG_OPS);)
        performPerImageOperations(
            Ipart, fn_output_img_root,
            my_current_nr_images + ipos, my_total_nr_images,
            tilt_deg, psi_deg,
            all_avg, all_stddev, all_minval, all_maxval
        );
        }

        {
        ifdefPREP_TIMING(TicToc tt (timer, TIMING_REST);)
        // Also store all the particles information in the STAR file
        FileName fn_img;
        if (Ipart().getDim() == 3) {
            fn_img.compose(fn_output_img_root, my_current_nr_images + ipos + 1, "mrc");
        } else {
            fn_img.compose(my_current_nr_images + ipos + 1, fn_output_img_root + ".mrcs"); // start image counting in stacks at 1!
        }
        MD.setValue(EMDL::IMAGE_NAME,      fn_img);
        MD.setValue(EMDL::MICROGRAPH_NAME, fn_mic);

        // Set the optics group for this particle to the one from the micrograph
        if (!MDin_has_optics_group) {
            int optics_group = MDmics.getValue<int>(EMDL::IMAGE_OPTICS_GROUP, imic);
            MD.setValue(EMDL::IMAGE_OPTICS_GROUP, optics_group);
        }

        // Also fill in the per-particle CTF parameters
        if (mic_star_has_ctf) {
            // Only set CTF parameters from the micrographs STAR file if the input STAR file did not contain it!
            if (!MDin_has_ctf || keep_ctf_from_micrographs) {
                RFLOAT maxres, fom;
                if (MDmics.containsLabel(EMDL::CTF_MAXRES)) {
                    maxres = MDmics.getValue<RFLOAT>(EMDL::CTF_MAXRES, imic);
                    MD.setValue(EMDL::CTF_MAXRES, maxres);
                }
                if (MDmics.containsLabel(EMDL::CTF_FOM)) {
                    fom = MDmics.getValue<RFLOAT>(EMDL::CTF_FOM, imic);
                    MD.setValue(EMDL::CTF_FOM, fom);
                }

                CtfHelper::write(ctf, MD);
            }

            // Only set beamtilt from the micrographs STAR file if the input STAR file did not contain it!
            if (!MDin_has_beamtilt || keep_ctf_from_micrographs) {
                if (MDmics.containsLabel(EMDL::IMAGE_BEAMTILT_X)) {
                    RFLOAT tilt_x = MDmics.getValue<RFLOAT>(EMDL::IMAGE_BEAMTILT_X, imic);
                    MD.setValue(EMDL::IMAGE_BEAMTILT_X, tilt_x);
                }
                if (MDmics.containsLabel(EMDL::IMAGE_BEAMTILT_Y)) {
                    RFLOAT tilt_y = MDmics.getValue<RFLOAT>(EMDL::IMAGE_BEAMTILT_Y, imic);
                    MD.setValue(EMDL::IMAGE_BEAMTILT_Y, tilt_y);
                }
            }

            // Copy rlnBeamTiltGroupName from the micrograph STAR file only when absent in the particle STAR file
            // for backwards compatibility with release 3.0
            if (!MDin_has_tiltgroup && MDmics.containsLabel(EMDL::PARTICLE_BEAM_TILT_CLASS)) {
                int tilt_class = MDmics.getValue<int>(EMDL::PARTICLE_BEAM_TILT_CLASS, imic);
                MD.setValue(EMDL::PARTICLE_BEAM_TILT_CLASS, tilt_class);
            }
        }

        }

        ipos++;
    }
}

void Preprocessing::runOperateOnInputFile() {
    Image<RFLOAT> Ipart, Iout;
    long int Nimg;

    if (fn_operate_in.isStarFile()) {
        REPORT_ERROR("ERROR: this functionality is no longer supported on STAR files. You can still operate on image stacks.");
    } else {
        // Read the header of the stack to see how many images there
        Iout.read(fn_operate_in, false);
        Nimg = Nsize(Iout());
    }

    RFLOAT all_avg = 0;
    RFLOAT all_stddev = 0;
    RFLOAT all_minval = LARGE_NUMBER;
    RFLOAT all_maxval = -LARGE_NUMBER;
    init_progress_bar(Nimg);
    int barstep = std::max(1, (int) Nimg / 120);
    for (long int i = 0; i < Nimg; i++) {
        FileName fn_tmp;

        // Read in individual miages from the stack
        Ipart.clear();
        Ipart.read(fn_operate_in, true, i);

        RFLOAT tilt_deg, psi_deg;
        tilt_deg = psi_deg = 0.0;
        performPerImageOperations(
            Ipart, fn_operate_out.withoutExtension(), i, Nimg,
            tilt_deg, psi_deg,
            all_avg, all_stddev, all_minval, all_maxval
        );

        // progress bar
        if (i % barstep == 0) progress_bar(i);

    }
    progress_bar(Nimg);

    std::cout << " Done writing to " << fn_operate_out << std::endl;
}

void Preprocessing::performPerImageOperations(
    Image<RFLOAT> &Ipart,
    FileName fn_output_img_root,
    long int image_nr,  long int nr_of_images,
    RFLOAT tilt_deg,    RFLOAT psi_deg,
    RFLOAT &all_avg,    RFLOAT &all_stddev,
    RFLOAT &all_minval, RFLOAT &all_maxval
) {

    Ipart().setXmippOrigin();

    if (do_rescale) rescale(Ipart, scale);

    if (do_rewindow) rewindow(Ipart, window);

    Ipart().setXmippOrigin();

    {
    ifdefPREP_TIMING(TicToc tt (timer, TIMING_NORMALIZE);)
    // 24 June 2015 - Shaoda, helical segments
    if (do_normalise) {
        RFLOAT bg_helical_radius = helical_tube_outer_diameter * 0.5 / angpix;
        if (do_rescale)
            bg_helical_radius *= scale / extract_size;
        normalise(
            Ipart, bg_radius, white_dust_stddev, black_dust_stddev, do_ramp,
            do_extract_helix, bg_helical_radius, tilt_deg, psi_deg
        );
    }
    }

    {
    ifdefPREP_TIMING(TicToc tt (timer, TIMING_INV_CONT);)
    if (do_invert_contrast) invert_contrast(Ipart);
    }

    // Calculate mean, stddev, min and max
    Stats<RFLOAT> stats = [&] () {
        ifdefPREP_TIMING(TicToc tt (timer, TIMING_COMP_STATS);)
        return computeStats(Ipart());
    }();

    if (Ipart().getDim() == 3) {
        Ipart.MDMainHeader.setValue(EMDL::IMAGE_STATS_MIN,    stats.min);
        Ipart.MDMainHeader.setValue(EMDL::IMAGE_STATS_MAX,    stats.max);
        Ipart.MDMainHeader.setValue(EMDL::IMAGE_STATS_AVG,    stats.avg);
        Ipart.MDMainHeader.setValue(EMDL::IMAGE_STATS_STDDEV, stats.stddev);
        Ipart.setSamplingRateInHeader(output_angpix);

        {
        ifdefPREP_TIMING(TicToc tt (timer, TIMING_PER_IMG_OP_WRITE);)
        // Write one mrc file for every subtomogram
        FileName fn_img;
        fn_img.compose(fn_output_img_root, image_nr + 1, "mrc");
        Ipart.write(fn_img);
        }
    } else {
        // Keep track of overall statistics
        all_minval = std::min(stats.min, all_minval);
        all_maxval = std::max(stats.max, all_maxval);
        all_avg	   += stats.avg;
        all_stddev += stats.stddev * stats.stddev;

        // Last particle: reset the min, max, avg and stddev values in the main header
        if (image_nr == nr_of_images - 1) {
            all_avg /= nr_of_images;
            all_stddev = sqrt(all_stddev / nr_of_images);
            Ipart.MDMainHeader.setValue(EMDL::IMAGE_STATS_MIN, all_minval);
            Ipart.MDMainHeader.setValue(EMDL::IMAGE_STATS_MAX, all_maxval);
            Ipart.MDMainHeader.setValue(EMDL::IMAGE_STATS_AVG, all_avg);
            Ipart.MDMainHeader.setValue(EMDL::IMAGE_STATS_STDDEV, all_stddev);
            Ipart.setSamplingRateInHeader(output_angpix);
        }

        {
        ifdefPREP_TIMING(TicToc tt (timer, TIMING_PER_IMG_OP_WRITE);)
        // Write this particle to the stack on disc
        // First particle: write stack in overwrite mode, from then on just append to it
        if (image_nr == 0) {
            Ipart.write(fn_output_img_root + ".mrcs", -1, nr_of_images > 1, WRITE_OVERWRITE);
        } else {
            Ipart.write(fn_output_img_root + ".mrcs", -1, false,            WRITE_APPEND);
        }
        }
    }
}

// Get the coordinate file from a given micrograph filename from MDdata
MetaDataTable Preprocessing::getCoordinateMetaDataTable(FileName fn_mic) {
    // Get the micropgraph name without the UNIQDATE string into fn_post, and only read that micrograph in the MDresult table
    FileName fn_pre, fn_jobnr, fn_post;
    decomposePipelineFileName(fn_mic, fn_pre, fn_jobnr, fn_post);

    MetaDataTable MDresult;

    // To upgrade to ObsModel, the following no longer works.
    // MDresult.read(fn_data, "particles", NULL, fn_post);
    for (long int i : MDimg) {
        FileName fn_mic_in_particle_star = MDimg.getValue<std::string>(EMDL::MICROGRAPH_NAME);
        if (fn_mic_in_particle_star.contains(fn_post))
            MDresult.addObject(MDimg.getObject(i));
    }

    RFLOAT mag2, dstep2, particle_angpix, rescale_fndata = 1.0;
    if (MDresult.numberOfObjects() > 0) {
        MDresult.goToObject(0);

        bool contains_xy = MDresult.containsLabel(EMDL::ORIENT_ORIGIN_X_ANGSTROM) && MDresult.containsLabel(EMDL::ORIENT_ORIGIN_Y_ANGSTROM);
        bool contains_z  = MDresult.containsLabel(EMDL::ORIENT_ORIGIN_Z_ANGSTROM);

        if (do_recenter && !contains_xy) {
            REPORT_ERROR("Preprocessing::initialise ERROR: input _data.star file does not contain rlnOriginX/YAngst for re-centering.");
        }

        RFLOAT xoff, yoff, zoff, rot, tilt, psi, xcoord, ycoord, zcoord, diffx, diffy, diffz;
        for (long int i : MDresult) {
            // Get the optics group for this particle
            /// FIXME: How can we guarantee optics group IDs are consistent among Part and Mic?
            int optics_group = obsModelMic.getOpticsGroup(MDresult);
            particle_angpix = obsModelPart.getPixelSize(optics_group);
            rescale_fndata = particle_angpix / angpix; // angpix is the pixel size of this micrograph

            // reset input offsets
            if (do_reset_offsets) {
                if (contains_xy) {
                    MDresult.setValue(EMDL::ORIENT_ORIGIN_X_ANGSTROM, 0.0);
                    MDresult.setValue(EMDL::ORIENT_ORIGIN_Y_ANGSTROM, 0.0);
                }
                if (contains_z) {
                    MDresult.setValue(EMDL::ORIENT_ORIGIN_Z_ANGSTROM, 0.0);
                }
            } else if (abs(rescale_fndata - 1.0) > 1e-6 || do_recenter) {
                // re-scale or re-center (irrelevant if do_reset_offsets)
                Matrix1D<RFLOAT> my_projected_center(3);
                my_projected_center.initZeros();

                xoff = MDresult.getValue<RFLOAT>(EMDL::ORIENT_ORIGIN_X_ANGSTROM);
                yoff = MDresult.getValue<RFLOAT>(EMDL::ORIENT_ORIGIN_Y_ANGSTROM);

                xoff /= particle_angpix; // Now in run_data's pixels
                yoff /= particle_angpix;

                if (do_recenter && (
                    fabs(recenter_x) > 0.0 ||
                    fabs(recenter_y) > 0.0 ||
                    fabs(recenter_z) > 0.0
                )) {
                    rot  = MDresult.getValue<RFLOAT>(EMDL::ORIENT_ROT);
                    tilt = MDresult.getValue<RFLOAT>(EMDL::ORIENT_TILT);
                    psi  = MDresult.getValue<RFLOAT>(EMDL::ORIENT_PSI);

                    // Project the center-coordinates
                    Matrix1D<RFLOAT> my_center(3);
                    Matrix2D<RFLOAT> A3D;
                    XX(my_center) = recenter_x; // in run_data's pixels
                    YY(my_center) = recenter_y;
                    ZZ(my_center) = recenter_z;
                    if (ref_angpix > 0)
                        my_center = my_center * (ref_angpix / particle_angpix);
                    A3D = Euler::angles2matrix(rot, tilt, psi);
                    my_projected_center = A3D * my_center;
                }

                xoff -= XX(my_projected_center);
                yoff -= YY(my_projected_center);
                xoff *= rescale_fndata; // now in micrograph's pixel
                yoff *= rescale_fndata;

                if (do_recenter) {
                    xcoord = MDresult.getValue<RFLOAT>(EMDL::IMAGE_COORD_X);
                    ycoord = MDresult.getValue<RFLOAT>(EMDL::IMAGE_COORD_Y);

                    diffx = xoff - round(xoff);
                    diffy = yoff - round(yoff);
                    xcoord -= round(xoff);
                    ycoord -= round(yoff);
                    xoff = diffx;
                    yoff = diffy;
                    MDresult.setValue(EMDL::IMAGE_COORD_X, xcoord);
                    MDresult.setValue(EMDL::IMAGE_COORD_Y, ycoord);
                }

                MDresult.setValue(EMDL::ORIENT_ORIGIN_X_ANGSTROM, angpix * xoff);
                MDresult.setValue(EMDL::ORIENT_ORIGIN_Y_ANGSTROM, angpix * yoff);
                if (contains_z) {
                    zoff = MDresult.getValue<RFLOAT>(EMDL::ORIENT_ORIGIN_Z_ANGSTROM);
                    zoff /= particle_angpix;
                    zoff -= ZZ(my_projected_center);
                    zoff *= rescale_fndata;

                    if (do_recenter) {
                        zcoord = MDresult.getValue<RFLOAT>(EMDL::IMAGE_COORD_Z);
                        diffz = zoff - round(zoff);
                        zcoord -= round(zoff);
                        zoff = diffz;
                        MDresult.setValue(EMDL::IMAGE_COORD_Z, zcoord);
                    }
                    MDresult.setValue(EMDL::ORIENT_ORIGIN_Z, angpix * zoff);
                }
            }
        }
    }
    return MDresult;
}

// Get the coordinate filename from the micrograph filename
FileName Preprocessing::getCoordinateFileName(FileName fn_mic) {
    FileName fn_pre, fn_jobnr, fn_post;
    decomposePipelineFileName(fn_mic, fn_pre, fn_jobnr, fn_post);
    FileName fn_coord = fn_coord_dir + fn_post.withoutExtension() + fn_coord_suffix;
    return fn_coord;
}

// Get the coordinate filename from the micrograph filename
FileName Preprocessing::getOutputFileNameRoot(FileName fn_mic) {
    FileName fn_pre, fn_jobnr, fn_post;
    decomposePipelineFileName(fn_mic, fn_pre, fn_jobnr, fn_post);
    FileName fn_part = fn_part_dir + fn_post.withoutExtension();
    return fn_part;
}
