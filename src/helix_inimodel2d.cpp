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


#include "helix_inimodel2d.h"
// #define DEBUG


void HelixAlignerModel::initialise(int nr_classes, int ydim, int xdim) {

    MultidimArray<RFLOAT> tmp = MultidimArray<RFLOAT>::zeros(ydim, xdim);
    tmp.setXmippOrigin();
    for (int iclass = 0; iclass < nr_classes; iclass++) {
        Aref.push_back(tmp);
        Asum.push_back(tmp);
        Asumw.push_back(tmp);
        pdf.push_back(0);
    }
    tmp.initZeros(ydim, ydim);
    for (int iclass = 0; iclass < nr_classes; iclass++) {
        Arec.push_back(tmp);
    }
}

void HelixAlignerModel::initZeroSums() {
    for (int iclass = 0; iclass < Asum.size(); iclass++) {
        Asum[iclass].initZeros();
        Asumw[iclass].initZeros();
        pdf[iclass] = 0.0;
    }
}

// Cleaning up
void HelixAlignerModel::clear() {
    Aref.clear();
    Arec.clear();
    Asum.clear();
    Asumw.clear();
    pdf.clear();
}

void HelixAligner::clear() {
    /// TODO: Clean up.
}

void HelixAligner::usage() {
    parser.writeUsage(std::cout);
}

void HelixAligner::parseInitial(int argc, char **argv) {

    parser.setCommandLine(argc, argv);

    // General optimiser I/O stuff
    int general_section = parser.addSection("General options");

    fn_out = parser.getOption("--o", "Output rootname","");
    fn_imgs = parser.getOption("--i", " STAR file with the input images and orientation parameters","");
    // deactivate fn_mics approach: never really worked...
    fn_mics = "";
    /*
    fn_mics = parser.getOption("--mic", "OR: STAR file with the input micrographs","");
    fn_coord_suffix = parser.getOption("--coord_suffix", "The suffix for the start-end coordinate files, e.g. \"_picked.star\" or \".box\"","");
    fn_coord_dir = parser.getOption("--coord_dir", "The directory where the coordinate files are (default is same as micrographs)", "ASINPUT");
    extract_width = textToInteger(parser.getOption("--extract_width", "Width (in pixels) of the images for the helices to be extracted ", "100"));
    */
    int param_section = parser.addSection("Parameters");
    crossover_distance = textToFloat(parser.getOption("--crossover_distance", "Distance in Angstroms between 2 cross-overs",""));
    nr_iter = textToInteger(parser.getOption("--iter", "Maximum number of iterations to perform", "10"));
    nr_classes = textToInteger(parser.getOption("--K", "Number of classes", "1"));
    angpix = textToFloat(parser.getOption("--angpix", "Pixel size in Angstroms (default take from STAR file)", "-1"));
    maxres = textToFloat(parser.getOption("--maxres", "Limit calculations to approximately this resolution in Angstroms", "-1"));
    max_shift_A = textToFloat(parser.getOption("--search_shift", "How many Angstroms to search translations perpendicular to helical axis?", "0"));
    max_rotate = textToFloat(parser.getOption("--search_angle", "How many degrees to search in-plane rotations?", "0"));
    step_rotate = textToFloat(parser.getOption("--step_angle", "The step size (in degrees) of the rotational searches", "1"));
    fn_inimodel = parser.getOption("--iniref", "An initial model to starting optimisation path", "");
    symmetry = textToInteger(parser.getOption("--sym", "Order of symmetry in the 2D xy-slice?", "1"));
    max_smear = textToInteger(parser.getOption("--smear", "Smear out each image along X to ensure continuity", "0"));
    random_seed = textToInteger(parser.getOption("--random_seed", "Random seed (default is with clock)", "-1"));
    search_size = textToInteger(parser.getOption("--search_size", "Search this many pixels up/down of the target downscaled size to fit best crossover distance", "5"));
    mask_diameter = textToFloat(parser.getOption("--mask_diameter", "The diameter (A) of a mask to be aplpied to the 2D reconstruction", "-1"));
    nr_threads = textToInteger(parser.getOption("--j", "Number of (openMP) threads", "1"));
    do_only_make_3d = parser.checkOption("--only_make_3d", "Take the iniref image, and create a 3D model from that without any alignment of the input images");

    verb = 1;

    if (parser.checkForErrors(verb))
        REPORT_ERROR("Errors encountered on the command line (see above), exiting...");
}


// Run multiple iterations


void HelixAligner::initialise() {

    // Randomise the order of the particles
    if (random_seed == -1) random_seed = time(NULL);
    // Also randomize random-number-generator for perturbations on the angles
    init_random_generator(random_seed);

    if (!fn_imgs.empty()) {
        // Get the image size
        MetaDataTable MD;
        MD.read(fn_imgs);
        FileName fn_img;
        Image<RFLOAT> img;
        const int i = MD.size() - 1;
        if (MD.containsLabel(EMDL::IMAGE_NAME)) {
            fn_img = MD.getValue<std::string>(EMDL::IMAGE_NAME, i);
        } else if (MD.containsLabel(EMDL::MLMODEL_REF_IMAGE)) {
            fn_img = MD.getValue<std::string>(EMDL::MLMODEL_REF_IMAGE, i);
        } else {
            REPORT_ERROR("ERROR: input STAR file does not contain rlnImageName or rlnReferenceImage!");
        }
        img.read(fn_img, false); // only read the header
        ori_size = img().xdim;
        if (img().xdim != img().ydim || img().zdim != 1)
            REPORT_ERROR("ERROR: only square 2D images are allowed.");

        // Get the pixel size
        if (MD.containsLabel(EMDL::CTF_MAGNIFICATION) && MD.containsLabel(EMDL::CTF_DETECTOR_PIXEL_SIZE)) {
            RFLOAT mag = MD.getValue<RFLOAT>(EMDL::CTF_MAGNIFICATION, i);
            RFLOAT dstep = MD.getValue<RFLOAT>(EMDL::CTF_DETECTOR_PIXEL_SIZE, i);
            RFLOAT my_angpix = 10000.0 * dstep / mag;
            std::cout << " Using pixel size from the input STAR file: " << my_angpix << std::endl;
            angpix = my_angpix;
        }

    } else if (!fn_mics.empty()) {
        // Read in the micrographs STAR file
        MDmics.read(fn_mics);

        // Get the pixel size
        if (MDmics.containsLabel(EMDL::CTF_MAGNIFICATION) && MDmics.containsLabel(EMDL::CTF_DETECTOR_PIXEL_SIZE)) {
            const long int i = MDmics.size() - 1;
            RFLOAT mag = MDmics.getValue<RFLOAT>(EMDL::CTF_MAGNIFICATION, i);
            RFLOAT dstep = MDmics.getValue<RFLOAT>(EMDL::CTF_DETECTOR_PIXEL_SIZE, i);
            RFLOAT my_angpix = 10000.0 * dstep / mag;
            std::cout << " Using pixel size from the input STAR file: " << my_angpix << std::endl;
            angpix = my_angpix;
        }

        // Make sure the coordinate file directory names end with a '/'
        if (fn_coord_dir != "ASINPUT" && fn_coord_dir[fn_coord_dir.length()-1] != '/')
            fn_coord_dir+="/";

        // Loop over all micrographs in the input STAR file and warn of coordinate file or micrograph file do not exist
        for (long int i : MDmics) {
            FileName fn_mic = MDmics.getValue<std::string>(EMDL::MICROGRAPH_NAME, i);
            FileName fn_pre, fn_jobnr, fn_post;
            decomposePipelineFileName(fn_mic, fn_pre, fn_jobnr, fn_post);
            FileName fn_coord = fn_coord_dir + fn_post.withoutExtension() + fn_coord_suffix;
            if (!exists(fn_coord))
                std::cerr << "Warning: coordinate file " << fn_coord << " does not exist..." << std::endl;
            if (!exists(fn_mic))
                std::cerr << "Warning: micrograph file " << fn_mic << " does not exist..." << std::endl;
        }

        ori_size = extract_width;
    } else if (do_only_make_3d && fn_inimodel != "") {
        Image<RFLOAT> img;
        img.read(fn_inimodel);
        img().setXmippOrigin();
        if (angpix < 0.0) {
            angpix = img.MDMainHeader.getValue<RFLOAT>(EMDL::IMAGE_SAMPLINGRATE_X, img.MDMainHeader.size() - 1);
            std::cout << " Using pixel size from the input file header: " << angpix << std::endl;
        }
        ori_size = img().xdim;

        // The 3D reconstruction
        float deg_per_pixel = 180.0 * angpix / crossover_distance;
        Image<RFLOAT> vol;
        vol().resize(ori_size, ori_size, ori_size);
        for (int k = 0; k < Zsize(vol()); k++) {
            float ang = deg_per_pixel * k;
            Matrix2D<RFLOAT> Arot = rotation2DMatrix(ang);

            MultidimArray<RFLOAT> Mrot = applyGeometry(img(), Arot, true, false);
            for (long int j = 0; j < Ysize(Mrot); j++)
            for (long int i = 0; i < Xsize(Mrot); i++) {
                direct::elem(vol(), i, j, k) = direct::elem(Mrot, i, j);
            }
        }
        vol.setSamplingRateInHeader(angpix);
        vol.write(fn_out + ".mrc");
        std::cout << " * Written " << fn_out << ".mrc" << std::endl;
        exit(RELION_EXIT_SUCCESS);
    } else {
        REPORT_ERROR("ERROR: provide --i, -mic, or --only_make_3d and --iniref");
    }

    if (angpix < 0.0) {
        REPORT_ERROR("ERROR: provide pixel size through --angpix or through the magnification and detectorpixel size in the input STAR file.");
    }

    if (maxres < 0.0 || maxres < 2.0 * angpix) {
        maxres = 2.0 * angpix;
        std::cout << " Setting maximum resolution to " << maxres << std::endl;
    }

    down_size = ori_size * angpix * 2.0 / maxres;

    // Make sure that the crossover distance is close to an integer times the (downsized) pixel size of the model!
    float best_fit = 1.;
    down_angpix = 0.;
    int best_size = 0;
    for (int delta_size = -search_size; delta_size <= search_size; delta_size += 2) {
        const int mysize = make_even(down_size + delta_size);
        if (mysize <= ori_size) {
            float myangpix = angpix * (float) ori_size / (float) mysize;
            float mydiv = 2.0 * crossover_distance / myangpix;
            // Also want even number of pixels in rectangle!
            float myfit = fmod(mydiv, 2);
            if (myfit > 1.0)
                myfit -= 2.0;
            myfit = fabs(myfit);
            if (myfit < best_fit) {
                best_fit = myfit;
                down_angpix = myangpix;
                best_size = mysize;
            }
            std::cout << " *   mydiv= " << mydiv << " myangpix= " << myangpix << " myfit= " << myfit << std::endl;
        }
    }
    std::cout <<     " *** best_angpix= " << down_angpix << " rectangles xsize= " << 2.0 * crossover_distance / down_angpix << std::endl;

    down_size = best_size;
    yrect  = make_even(round(ori_size * angpix / down_angpix));
    xrect  = round(2.0 * crossover_distance / down_angpix);
    model.initialise(nr_classes, yrect, xrect);
    max_shift = ceil(max_shift_A / down_angpix);
    mask_radius_pix = mask_diameter > 0 ? ceil(mask_diameter / (2.0 * down_angpix)) : yrect / 2 - 2;
    std::cout << " maxres= " << maxres << " angpix= " << angpix << " down_size= " << down_size << std::endl;
    std::cout << " xrect= " << xrect << " yrect= " << yrect << " down_angpix= " << down_angpix << std::endl;
    std::cout << " max_shift= " << max_shift << " mask_radius_pix= "<< mask_radius_pix<< std::endl;

    // Now read in all images
    if (fn_mics == "") readImages();
    else getHelicesFromMics();

    initialiseClasses();

}

// Read in all the images
void HelixAligner::readImages() {

    MD.read(fn_imgs);

    if (verb > 0) {
        std::cout << " Reading in all images ..." << std::endl;
        init_progress_bar(MD.size());
    }

    for (long int ipart : MD) {

        FileName fn_img;
        Image<RFLOAT> img;
        const long int i = MD.size() - 1;
        if (MD.containsLabel(EMDL::IMAGE_NAME)) {
            fn_img = MD.getValue<std::string>(EMDL::IMAGE_NAME, i);
        } else if (MD.containsLabel(EMDL::MLMODEL_REF_IMAGE)) {
            fn_img = MD.getValue<std::string>(EMDL::MLMODEL_REF_IMAGE, i);
        } else {
            REPORT_ERROR("ERROR: input STAR file does not contain rlnImageName or rlnReferenceImage!");
        }
        img.read(fn_img);
        img().setXmippOrigin();
        // Rethink this when expanding program to 3D!
        RFLOAT yoff = MD.containsLabel(EMDL::ORIENT_ORIGIN_Y_ANGSTROM) ? MD.getValue<RFLOAT>(EMDL::ORIENT_ORIGIN_Y_ANGSTROM, i) : 0.0;
        RFLOAT psi  = MD.containsLabel(EMDL::ORIENT_PSI)               ? MD.getValue<RFLOAT>(EMDL::ORIENT_PSI, i)               : 0.0;
        ori_psis.push_back(psi);
        ori_yoffs.push_back(yoff);
        // Apply the actual transformation
        Matrix2D<RFLOAT> A = rotation2DMatrix(psi);
        A.at(1, 2) = -yoff / angpix;
        img() = applyGeometry(img(), A, IS_INV, DONT_WRAP);

        Xrects.emplace_back();

        // Calculate all rotated versions
        if (ipart == 0) psis.clear();

        for (int iflip = 0; iflip < 2; iflip++) {
            for (RFLOAT ang = 0; ang <= max_rotate; ang += step_rotate) {
                RFLOAT myang = iflip == 1 ? ang + 180.0 : ang;
                Matrix2D<RFLOAT> Arot = rotation2DMatrix(myang);
                MultidimArray<RFLOAT> Irot = applyGeometry(img(), Arot, true, false);
                resizeMap(Irot, down_size);
                Irot.setXmippOrigin();
                Xrects[Xrects.size() - 1].push_back(Irot);

                if (ipart == 0) psis.push_back(myang);

                if (ang > 0.0) {
                    // Also rotate in the opposite direction
                    Irot = applyGeometry(img(), Arot, false, false);
                    resizeMap(Irot, down_size);
                    Irot.setXmippOrigin();
                    Xrects[Xrects.size()-1].push_back(Irot);

                    if (ipart == 0) psis.push_back(-myang);
                }
            }
        }

        if (verb > 0 && ipart % 50 == 0)
            progress_bar(ipart);

        // #define DEBUG_READIMAGES
        #ifdef DEBUG_READIMAGES
        const auto fnt = FileName::compose("helixnew", Xrects.size(),"spi",3);
        Image<RFLOAT>(Xrects[Xrects.size() - 1][3]).write(fnt);
        #endif
    }

    if (verb > 0)
        progress_bar(MD.size());

    #ifdef DEBUG
    std::cerr << "done readImages" << std::endl;
    #endif

}

void HelixAligner::getHelicesFromMics() {
    if (verb > 0) {
        std::cout << " Reading in all micrographs ..." << std::endl;
        init_progress_bar(MDmics.size());
    }

    // Loop over all micrographs in the input STAR file and warn of coordinate file or micrograph file do not exist
    // Surely imic should be part of the for loop?
    long int imic = 0;
    for (long int i : MDmics) {
        imic++;
        FileName fn_mic = MDmics.getValue<std::string>(EMDL::MICROGRAPH_NAME, i);
        FileName fn_pre, fn_jobnr, fn_post;
        decomposePipelineFileName(fn_mic, fn_pre, fn_jobnr, fn_post);
        FileName fn_coord = fn_coord_dir + fn_post.withoutExtension() + fn_coord_suffix;
        if (!exists(fn_mic) || !exists(fn_coord)) {
            if (!exists(fn_mic))
                std::cerr << "Warning: micrograph file " << fn_mic << " does not exist..." << std::endl;
            if (!exists(fn_coord))
                std::cerr << "Warning: coordinate file " << fn_coord << " does not exist..." << std::endl;
        } else {
            Image<RFLOAT> Imic;
            Imic.read(fn_mic);
            RFLOAT avg = average(Imic());

            // Read in the coordinate files
            MetaDataTable MDcoords;
            MDcoords.read(fn_coord);
            if (MDcoords.size() % 2 == 1) {
                std::cerr << " ERROR: odd number of entries in " << fn_coord << "! Skipping this micrograph... " << std::endl;
                continue;
            }

            // Get all start-end coordinate pairs
            std::vector<RFLOAT> x1_coord_list, y1_coord_list, x2_coord_list, y2_coord_list, pitch_list;
            RFLOAT xp, yp;
            for (long int i : MDcoords) {
                xp = MDcoords.getValue<RFLOAT>(EMDL::IMAGE_COORD_X, i);
                yp = MDcoords.getValue<RFLOAT>(EMDL::IMAGE_COORD_Y, i);
                if (i % 2 == 0) {
                    x1_coord_list.push_back(xp);
                    y1_coord_list.push_back(yp);
                } else {
                    x2_coord_list.push_back(xp);
                    y2_coord_list.push_back(yp);
                }
            }

            // Now extract the images: make all helices stand upright... Y becomes helical axis, X becomes helix width
            // For that we need to do interpolations...
            for (int ipair = 0; ipair < x1_coord_list.size(); ipair++) {
                std::vector<MultidimArray<RFLOAT>> dummy;
                Xrects.push_back(dummy);

                // Calculate all rotated versions
                int oldxsize, oldysize, oldsize;
                bool do_set_oldsize = true;
                // For all angles
                for (RFLOAT ang = 0.0; ang <= max_rotate; ang += step_rotate) {
                    RFLOAT x1,x2,y1,y2,xcen,ycen;
                    x1 = x1_coord_list[ipair];
                    x2 = x2_coord_list[ipair];
                    y1 = y1_coord_list[ipair];
                    y2 = y2_coord_list[ipair];
                    xcen = x1 + (x2 - x1) / 2;
                    ycen = y1 + (y2 - y1) / 2;

                    int xsize = floor(sqrt((x2 - x1) * (x2 - x1) + (y2 - y1) * (y2 - y1)));
                    RFLOAT phi = degrees(atan(RFLOAT(y2 - y1) / RFLOAT(x2 - x1)));
                    MultidimArray<RFLOAT> Ihelix;
                    Ihelix.resize(extract_width, xsize);
                    Ihelix.setXmippOrigin();

                    int nrots = ang > 0.0 ? 2 : 1;
                    // For for positive and negative rotations
                    for (int irot = 0; irot < nrots; irot++) {
                        Matrix2D<RFLOAT> Arot = rotation2DMatrix(irot == 0 ? phi + ang : phi - ang);
                        Arot(0, 2) = xcen;
                        Arot(1, 2) = ycen;

                        int m1, n1, m2, n2;
                        RFLOAT x, y, xp, yp;
                        RFLOAT wx, wy;

                        // Find center and limits of image
                        int cen_y  = (int) (Ysize(Ihelix) / 2);
                        int cen_x  = (int) (Xsize(Ihelix) / 2);
                        int cen_yp = (int) (Ysize(Imic()) / 2);
                        int cen_xp = (int) (Xsize(Imic()) / 2);
                        RFLOAT minxp  = 0;
                        RFLOAT minyp  = 0;
                        RFLOAT maxxp  = Xsize(Imic()) - 1;
                        RFLOAT maxyp  = Ysize(Imic()) - 1;
                        int Xdim   = Xsize(Imic());
                        int Ydim   = Ysize(Imic());

                        for (int i = 0; i < Ysize(Ihelix); i++) {
                            // Calculate position of the beginning of the row in the output image
                            x = -cen_x;
                            y = i - cen_y;

                            // Calculate this position in the input image according to the
                            // geometrical transformation
                            // they are related by
                            // coords_output(=x,y) = A * coords_input (=xp,yp)
                            xp = x * Arot(0, 0) + y * Arot(0, 1) + Arot(0, 2);
                            yp = x * Arot(1, 0) + y * Arot(1, 1) + Arot(1, 2);

                            for (int j = 0; j < Xsize(Ihelix); j++) {
                                bool interp;
                                RFLOAT tmp;

                                // If the point is outside the image, apply a periodic extension
                                // of the image, what exits by one side enters by the other
                                interp = true;
                                if (xp < minxp  || xp > maxxp) {
                                    interp = false;
                                }

                                if (yp < minyp  || yp > maxyp) {
                                    interp = false;
                                }

                                if (interp) {
                                    // Linear interpolation

                                    // Calculate the integer position in input image, be careful
                                    // that it is not the nearest but the one at the top left corner
                                    // of the interpolation square. Ie, (0.7,0.7) would give (0,0)
                                    // Calculate also weights for point m1+1,n1+1
                                    wx = xp;// + cen_xp;
                                    m1 = (int) wx;
                                    wx = wx - m1;
                                    m2 = m1 + 1;
                                    wy = yp;// + cen_yp;
                                    n1 = (int) wy;
                                    wy = wy - n1;
                                    n2 = n1 + 1;

                                    // Perform interpolation
                                    // if wx == 0 means that the rightest point is useless for this
                                    // interpolation, and even it might not be defined if m1=xdim-1
                                    // The same can be said for wy.
                                    tmp = (RFLOAT) ((1 - wy) * (1 - wx) * direct::elem(Imic(), m1, n1));

                                    if (m2 < Xdim)
                                        tmp += (RFLOAT) ((1 - wy) * wx * direct::elem(Imic(), m2, n1));

                                    if (n2 < Ydim) {
                                        tmp += (RFLOAT) (wy * (1 - wx) * direct::elem(Imic(), m1, n2));

                                        if (m2 < Xdim) {
                                            tmp += (RFLOAT) (wy * wx * direct::elem(Imic(), m2, n2));
                                        }
                                    }

                                    direct::elem(Ihelix, i, j) = tmp;

                                } else {
                                    direct::elem(Ihelix, i, j) = avg;
                                }

                                // Compute new point inside input image
                                xp += Arot(0, 0);
                                yp += Arot(1, 0);
                            }
                        }
                        // #define DEBUG_GETHELICESFROMMICS
                        #ifdef DEBUG_GETHELICESFROMMICS
                        const auto fntt = FileName::compose("helixnew1_beforedown", Xrects.size(), "spi", 3);
                        Image<RFLOAT>(Ihelix).write(fntt);
                        #endif

                        // Downscale if needed
                        MultidimArray<RFLOAT> Idown = Ihelix;
                        if (down_angpix > angpix) {
                            RFLOAT avg = average(Idown);

                            int oldxsize = Xsize(Idown);
                            int oldysize = Ysize(Idown);
                            int oldsize = oldxsize;

                            if (oldxsize != oldysize) {
                                oldsize = std::max(oldxsize, oldysize);
                                Idown.setXmippOrigin();
                                Idown.windowed(
                                    Xmipp::init(oldsize), Xmipp::init(oldsize),
                                    Xmipp::last(oldsize), Xmipp::last(oldsize),
                                    avg
                                );
                            }
                            int newsize = make_even(round(oldsize * angpix / down_angpix));
                            resizeMap(Idown, newsize);

                            if (oldxsize != oldysize) {
                                int newxsize = make_even(round(oldxsize * angpix / down_angpix));
                                int newysize = make_even(round(oldysize * angpix / down_angpix));
                                Idown.setXmippOrigin();
                                Idown.windowed(
                                    Xmipp::init(newysize), Xmipp::init(newxsize),
                                    Xmipp::last(newysize), Xmipp::last(newxsize)
                                );
                            }

                        }

                        // Ad-hoc image normalisation
                        const auto stats = computeStats(Idown);
                        Idown = (Idown - stats.avg) / -stats.stddev;  // Invert contrast

                        Xrects[Xrects.size() - 1].push_back(Idown);
                    }
                }

                if (verb > 0)
                    progress_bar(imic);

                // #define DEBUG_GETHELICESFROMMICS2
                #ifdef DEBUG_GETHELICESFROMMICS2
                const auto fnt = FileName::compose("helixnew1", Xrects.size(), "spi", 3);
                Image<RFLOAT>(Xrects[Xrects.size() - 1][1]).write(fnt);
                #endif
            }
        }
    }
    if (verb > 0)
        progress_bar(MDmics.size());
}

void HelixAligner::initialiseClasses() {
    #ifdef DEBUG
    std::cerr << "Entering initialiseClasses" << std::endl;
    #endif
    if (model.Aref.empty())
        REPORT_ERROR("BUG: non-initialised model!");

    if (verb > 0)
        std::cout << " Initialising reference(s) ..." << std::endl;

    if (fn_inimodel != "") {

        if (nr_classes > 1)
            REPORT_ERROR("ERROR: can only use initial reference for single-class!");

        Image<RFLOAT> Iref;
        Iref.read(fn_inimodel);
        resizeMap(Iref(), Ysize(model.Aref[0]));
        Iref().setXmippOrigin();
        std::cerr << " model.Arec.size()= " << model.Arec.size() << std::endl;
        model.Arec[0] = Iref();
        // Now project the reconstruction back out into the model.Aref[iclass]
        Projector PP(Ysize(model.Aref[0]), TRILINEAR, 2, 1, 1);
        // Set the FT of img inside the Projector
        MultidimArray<RFLOAT> dummy;
        PP.computeFourierTransformMap(Iref(), dummy, Ysize(model.Aref[0]), 1);

        // Calculate all projected lines
        for (int i = 0; i < Xsize(model.Aref[0]); i++) {
            FourierTransformer transformer;

            RFLOAT rot = (RFLOAT) i * 360.0 / (Xsize(model.Aref[0]));
            auto A2D = rotation2DMatrix(rot);
            auto myFline = PP.get2DFourierTransform(Ysize(model.Aref[0]) / 2 + 1, 1, 1, A2D);
            MultidimArray<RFLOAT> myline =
                transformer.inverseFourierTransform(myFline).setXmippOrigin();
            // Shift the image back to the center...
            CenterFFT(myline, false);
            for (int j = 0; j < Ysize(model.Aref[0]); j++)
                direct::elem(model.Aref[0], i, j) = direct::elem(myline, j);
        }

    #define DEBUGREC2D
    #ifdef DEBUGREC2D
        Image<RFLOAT> It;
        It() = model.Aref[0];
        It.write("after_reproject.spi");
    #endif

    } else {
        // Randomly position all particles along the X-direction
        model.initZeroSums();
        // Loop over all particles
        if (verb > 0)
            init_progress_bar(Xrects.size());

        for (int ipart = 0; ipart < Xrects.size(); ipart++) {
            // Set into a random class
            int myclass        = rnd_unif() * nr_classes;
            int random_xoffset = rnd_unif() * xrect;
            for (int i_smear = -max_smear; i_smear <= max_smear; i_smear++) {

                double smearw = max_smear == 0 ? 1 : gaussian1D((double) i_smear, (double) max_smear / 3);
                FOR_ALL_ELEMENTS_IN_ARRAY2D(Xrects[ipart][0], i, j) {
                    int ip = i + random_xoffset + i_smear;
                    while (ip < Xinit(model.Aref[myclass])) { ip += xrect; }
                    while (ip > Xlast(model.Aref[myclass])) { ip -= xrect; }

                    // this places the original image in the offset-translated center of the rectangle
                    model.Asum [myclass].elem(ip, j) += smearw * Xrects[ipart][0].elem(i, j);
                    model.Asumw[myclass].elem(ip, j) += smearw;

                    // This places the Y-flipped image at half a cross-over distance from the first one
                    int jp = -j;
                    if (jp >= Yinit(Xrects[ipart][0]) && jp <= Ylast(Xrects[ipart][0])) {
                        int iip = ip + xrect / 2;
                        while (iip < Xinit(model.Aref[myclass])) { iip += xrect; }
                        while (iip > Xlast(model.Aref[myclass])) { iip -= xrect; }
                        model.Asum [myclass].elem(iip, jp) += smearw * Xrects[ipart][0].elem(i, j);
                        model.Asumw[myclass].elem(iip, jp) += smearw;

                    }
                }
            }
            model.pdf[myclass] += 1.0;
            if (verb > 0)
                progress_bar(ipart);
        }

        if (verb > 0)
            progress_bar(Xrects.size());

        // After all images have been set, maximise the references in the model
        maximisation();
    }
    #ifdef DEBUG
    std::cerr << "Leaving initialiseClasses" << std::endl;
    #endif
}

void HelixAligner::expectationOneParticleNoFFT(long int ipart) {

    int twostarty = 2 * Xmipp::init(yrect);
    double maxccf = -100.;
    int best_class = -1;
    int best_k_rot = -1;
    int best_i_offset = -1;
    int best_j_offset = -1;
    for (int iclass = 0; iclass < nr_classes; iclass++) {

        for (int k_rot = 0; k_rot < Xrects[ipart].size(); k_rot++) {
            for (int j_offset = -max_shift; j_offset <= max_shift; j_offset++) {
            for (int i_offset = 0; i_offset < xrect; i_offset++) {
                double ccf_xa = 0;
                double ccf_x2 = 0;
                double ccf_a2 = 0;
                for (long int j = Yinit(Xrects[ipart][k_rot]); j <= Ylast(Xrects[ipart][k_rot]); j++) {

                    int jp = j + j_offset;
                    if (jp < -mask_radius_pix || jp > mask_radius_pix)
                        continue;

                    /*
                    while (jp < Yinit(model.Aref[iclass]))
                        jp += yrect;
                    while (jp > Ylast(model.Aref[iclass]))
                        jp -= yrect;
                    */

                    for (long int i = Xinit(Xrects[ipart][k_rot]); i <= Xlast(Xrects[ipart][k_rot]); i++) {
                        int ip = i + i_offset;
                        while (ip < Xinit(model.Aref[iclass])) { ip += xrect; }
                        while (ip > Xlast(model.Aref[iclass])) { ip -= xrect; }

                        // This places the Y-flipped image at half a cross-over distance from the first one
                        int jpp = -jp;
                        // Don't let the image run out of the height of the box
                        if (jpp >= Yinit(Xrects[ipart][k_rot]) && jpp <= Ylast(Xrects[ipart][k_rot])) {
                            int ipp = ip + xrect / 2;
                            while (ipp < Xinit(model.Aref[iclass])) { ipp += xrect; }
                            while (ipp > Xlast(model.Aref[iclass])) { ipp -= xrect; }

                            // this places the original image in the offset-translated center of the rectangle
                            ccf_xa += model.Aref[iclass].elem(ip, jp) * Xrects[ipart][k_rot].elem(i, j);
                            ccf_a2 += model.Aref[iclass].elem(ip, jp) * model.Aref[iclass].elem(ip, jp);

                            ccf_xa += model.Aref[iclass].elem(ipp, jpp) * Xrects[ipart][k_rot].elem(i, j);
                            ccf_a2 += model.Aref[iclass].elem(ipp, jpp) * model.Aref[iclass].elem(ipp, jpp);

                            ccf_x2 += 2.0 * Xrects[ipart][k_rot].elem(i, j) * Xrects[ipart][k_rot].elem(i, j);

                        }
                    }
                }

                double ccf = ccf_x2 > 0.0 && ccf_a2 > 0.0 ? ccf_xa / (sqrt(ccf_x2) * sqrt(ccf_a2)) : 0.0;

                // Find the best fit
                if (ccf > maxccf) {
                    maxccf = ccf;
                    best_class = iclass;
                    best_k_rot = k_rot;
                    best_i_offset = i_offset;
                    best_j_offset = j_offset;
                }
            }
            }
        }
    }


    if (maxccf < -1.0)
        REPORT_ERROR("BUG: not found maxccf!");

    // Now set the optimal Y-translations and rotations in the output STAR file
    RFLOAT psi = ori_psis[ipart] + psis[best_k_rot];
    RFLOAT yoff = ori_yoffs[ipart] + best_i_offset * down_angpix;

    MD.setValue(EMDL::ORIENT_ORIGIN_Y_ANGSTROM, yoff, ipart);
    MD.setValue(EMDL::ORIENT_PSI,               psi,  ipart);

    // Now add the image to that class reference
    // To ensure continuity in the reference: smear out every image along X
    #pragma omp critical
    {
        for (int i_smear = -max_smear; i_smear <= max_smear; i_smear++) {
            double smearw = max_smear < Xmipp::epsilon ? 1 : gaussian1D((double) i_smear, (double) max_smear / 3);
            FOR_ALL_ELEMENTS_IN_ARRAY2D(Xrects[ipart][best_k_rot], i, j) {

                int ip = i + best_i_offset + i_smear;
                while (ip < Xinit(model.Aref[best_class])) { ip += xrect; }
                while (ip > Xlast(model.Aref[best_class])) { ip -= xrect; }

                int jp = j + best_j_offset;
                while (jp < Yinit(model.Aref[best_class])) { jp += yrect; }
                while (jp > Ylast(model.Aref[best_class])) { jp -= yrect; }

                // this places the original image in the offset-translated center of the rectangle
                model.Asum [best_class].elem(ip, jp) += smearw * Xrects[ipart][best_k_rot].elem(i, j);
                model.Asumw[best_class].elem(ip, jp) += smearw;

                // This places the Y-flipped image at half a cross-over distance from the first one
                int jpp = -jp;
                if (jpp >= Yinit(Xrects[ipart][best_k_rot]) && jpp <= Ylast(Xrects[ipart][best_k_rot])) {
                    int ipp = ip + xrect / 2;
                    while (ipp > Xlast(model.Aref[best_class])) { ipp -= xrect; }
                    model.Asum [best_class].elem(ipp, jpp) += smearw * Xrects[ipart][best_k_rot].elem(i, j);
                    model.Asumw[best_class].elem(ipp, jpp) += smearw;
                }
            }
        }
        model.pdf[best_class] += 1.0;
    }
}

void HelixAligner::expectation() {

    // Initialise the wsum_model to zeros
    model.initZeroSums();

    if (verb > 0) { init_progress_bar(Xrects.size()); }

    #pragma omp parallel for num_threads(nr_threads)
    for (long int ipart = 0; ipart < Xrects.size(); ipart++) {

        expectationOneParticleNoFFT(ipart);

        if (ipart % nr_threads == 0)
            progress_bar(ipart);
    }

    progress_bar(Xrects.size());

}

void HelixAligner::maximisation() {
    #ifdef DEBUGREC2D
    Image<RFLOAT> It;
    It() = model.Asumw[0];
    It.write("Asumw.spi");
    It() = model.Asum[0];
    It.write("Asum.spi");
    #endif

    // Update the references
    double allsum = 0.;
    for (int iclass = 0; iclass < nr_classes; iclass++) {
        for (int j = 0; j < yrect; j++) {
        for (int i = 0; i < xrect; i++) {
            if (direct::elem(model.Asumw[iclass], i, j) > 0.0) {
                direct::elem(model.Aref[iclass], i, j) =  direct::elem(model.Asum[iclass], i, j) / direct::elem(model.Asumw[iclass], i, j);
            } else {
                direct::elem(model.Aref[iclass], i, j) = 0.0;
            }

            // Also store  sum of classes in Asum for writeOut
            direct::elem(model.Asum[iclass], i, j)  = direct::elem(model.Aref[iclass], i, j);
        }
        }
        allsum += model.pdf[iclass];

        reconstruct2D(iclass);
    }
    for (int iclass = 0; iclass < nr_classes; iclass++) {
        model.pdf[iclass] /= allsum;
    }

}

void HelixAligner::reconstruct2D(int iclass) {
    #ifdef DEBUG
    std::cerr << "Entering reconstruct2D" << std::endl;
    #endif
    #ifdef DEBUGREC2D
    Image<RFLOAT> It;
    It() = model.Aref[iclass];
    It.write("before_reproject.spi");
    #endif

    // Loop over the length of the helix to get the transforms of all 1D images
    std::vector<MultidimArray<Complex> > myFlines;
    for (int i = 0; i < Xsize(model.Aref[iclass]); i++) {
        MultidimArray<RFLOAT> myline(Ysize(model.Aref[iclass]));

        for (int j = 0; j < Ysize(model.Aref[iclass]); j++)
            direct::elem(myline, j) = direct::elem(model.Aref[iclass], i, j);

        CenterFFT(myline, true);
        FourierTransformer transformer;
        MultidimArray<Complex> &myFline = transformer.FourierTransform(myline);
        myFlines.push_back(myFline);
    }

    // Then reconstruct
    BackProjector BP(Ysize(model.Aref[iclass]), 2, "C1", TRILINEAR, 2, 1, 0, 1.9, 15, 1, false);
    BP.initialiseDataAndWeight(Ysize(model.Aref[iclass]));

    for (int j = 0; j < myFlines.size(); j++) {
        RFLOAT rot = (RFLOAT) j * 360.0 / Xsize(model.Aref[iclass]);
        Matrix2D<RFLOAT> A2D = rotation2DMatrix(rot);
        BP.set2DFourierTransform(myFlines[j], A2D);
    }
    MultidimArray<RFLOAT> tau2;
    model.Arec[iclass] = BP.reconstruct(10, false, tau2);

    if (symmetry > 1) {

        #ifdef DEBUGREC2D
        It() = model.Arec[iclass];
        resizeMap(It(), ori_size);
        It.write("rec_beforesym.spi");
        #endif

        MultidimArray<RFLOAT> Asum = model.Arec[iclass];
        for (int i = 1; i < symmetry; i++) {
            RFLOAT ang = i * 360.0 / (RFLOAT) symmetry;
            Matrix2D<RFLOAT> A2D = rotation2DMatrix(ang);
            Asum += applyGeometry(model.Arec[iclass], A2D, false, false);
        }
        model.Arec[iclass] = Asum / (RFLOAT) symmetry;
    }

    if (mask_diameter > 0.0) {
        RFLOAT pixel_radius = mask_diameter / (2.0 * down_angpix);
        softMaskOutsideMap(model.Arec[iclass], pixel_radius, 0.0);
    }

    #ifdef DEBUGREC2D
    It() = model.Arec[iclass];
    resizeMap(It(), ori_size);
    It.write("rec.spi");
    #endif

    // Now project the reconstruction back out into the model.Aref[iclass]
    Projector PP(Ysize(model.Aref[iclass]), TRILINEAR, 2, 1, 1);
    // Set the FT of img inside the Projector
    MultidimArray<RFLOAT> power_spectrum;
    PP.computeFourierTransformMap(model.Arec[iclass], power_spectrum, Ysize(model.Aref[iclass]), 1);

    // Calculate all projected lines
    for (int i = 0; i < myFlines.size(); i++) {

        const RFLOAT rot = (RFLOAT) i * 360.0 / Xsize(model.Aref[iclass]);
        const Matrix2D<RFLOAT> A2D = rotation2DMatrix(rot);
        myFlines[i] = PP.get2DFourierTransform(
            myFlines[i].xdim, myFlines[i].ydim, myFlines[i].zdim, A2D);
        FourierTransformer transformer;
        MultidimArray<RFLOAT> line = transformer.inverseFourierTransform(myFlines[i]);
        // Shift the image back to the center...
        CenterFFT(line, false);

        for (int j = 0; j < Ysize(model.Aref[iclass]); j++)
            direct::elem(model.Aref[iclass], i, j) = direct::elem(line, j);

    }
    #ifdef DEBUGREC2D
    It() = model.Aref[iclass];
    It.write("after_reproject.spi");
    #endif
    #ifdef DEBUG
    std::cerr << "Leaving reconstruct2D" << std::endl;
    #endif
}

static void write_subroutine(Image<RFLOAT> &img, const std::vector<MultidimArray<RFLOAT>> &v, const FileName &fn) {
    for (int iclass = 0; iclass < Nsize(img()); iclass++) {
        const auto &arr = v[iclass];
        for (long int j = 0; j < Ysize(arr); j++)
        for (long int i = 0; i < Xsize(arr); i++)
            direct::elem(img(), i, j, 0, iclass) = direct::elem(arr, i, j);
    }
    img.write(fn);
}

void HelixAligner::writeOut(int iter) {

    // std::cout << " **** Model for iteration " << iter << std::endl;

    #ifdef DEBUG
    for (int iclass = 0; iclass < nr_classes; iclass++) {
        FileName fn_class = fn_out + "_it" + integerToString(iter, 3) + "_class" + integerToString(iclass + 1, 3) + ".spi";
        Image<RFLOAT> Ic;
        Ic() = model.Aref[iclass];
        Ic.write(fn_class);
        std::cout << " * Written " << fn_class << std::endl;
        fn_class = fn_out + "_it" + integerToString(iter, 3) + "_class" + integerToString(iclass + 1, 3) + "_reconstructed.spi";
        Ic() = model.Arec[iclass];
        Ic.write(fn_class);
        std::cout << " * Written " << fn_class << std::endl;
    }
    #else
    FileName fn_iter = fn_out + "_it" + integerToString(iter, 3);
    MD.write(fn_iter + ".star");

    Image<RFLOAT> Aimg (xrect, yrect, 1, nr_classes);
    write_subroutine(Aimg, model.Aref, fn_iter + "_reprojections.mrcs");
    write_subroutine(Aimg, model.Asum, fn_iter + "_summed_classes.mrcs");

    Image<RFLOAT> Aimg2 (yrect, yrect, 1, nr_classes);
    write_subroutine(Aimg2, model.Arec, fn_iter + "_reconstructed.mrcs");
    #endif

    if (nr_classes > 1) {
        for (int iclass = 0; iclass < nr_classes; iclass++) {
            std:: cout << " * Fraction class " << iclass + 1 << " = " << model.pdf[iclass] << std::endl;
        }
    }

}

void HelixAligner::reconstruct3D() {

    for (int iclass = 0; iclass < nr_classes; iclass++) {
        FileName fn_class = fn_out + "_class" + integerToString(iclass + 1, 3) + "_projections.spi";
        Image<RFLOAT> Ic;
        Ic() = model.Aref[iclass];
        Ic.setSamplingRateInHeader(angpix);
        Ic.write(fn_class);
        std::cout << " * Written " << fn_class << std::endl;
        fn_class = fn_out + "_class" + integerToString(iclass + 1, 3) + "_rec2d.spi";
        Ic() = model.Arec[iclass];
        resizeMap(Ic(), ori_size);
        Ic.setSamplingRateInHeader(angpix);
        Ic.write(fn_class);
        MultidimArray<RFLOAT> Mori = Ic();
        std::cout << " * Written " << fn_class << std::endl;

        // The 3D reconstruction
        float deg_per_pixel = 180.0 * angpix / crossover_distance;
        Ic().resize(ori_size, ori_size, ori_size);
        for (int k = 0; k < Zsize(Ic()); k++) {
            float ang = deg_per_pixel * k;
            Matrix2D<RFLOAT> Arot = rotation2DMatrix(ang);

            MultidimArray<RFLOAT> Mrot = applyGeometry(Mori, Arot, true, false);
            for (long int j = 0; j < Ysize(Mrot); j++) \
            for (long int i = 0; i < Xsize(Mrot); i++) {
                direct::elem(Ic(), i, j, k) = direct::elem(Mrot, i, j);
            }
        }
        fn_class = fn_out + "_class" + integerToString(iclass + 1, 3) + "_rec3d.mrc";
        Ic.setSamplingRateInHeader(angpix);
        Ic.write(fn_class);
        std::cout << " * Written " << fn_class << std::endl;
    }

}

// Run multiple iterations
void HelixAligner::run() {

    // Write out the starting model as well
    writeOut(0);

    int decrease_smear = round((float) max_smear / (float) (nr_iter + 5));
    for (int iter = 1; iter <= nr_iter; iter++) {

        if (verb > 0) {
            std::cout << " Iteration " << iter <<" of " << nr_iter << std::endl;
            if (max_smear > 0)
                std::cout << "  = smearing references by " << max_smear << " downsampled pixels along helical axis " << std::endl;
        }

        expectation();

        maximisation();

        writeOut(iter);

        if (max_smear > 0)
            max_smear -= decrease_smear;

    }

    // Reconstruct the final solution in 3D
    reconstruct3D();
}
