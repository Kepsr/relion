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
#include <src/strings.h>
#include <src/funcs.h>
#include <src/memory.h>
#include <src/euler.h>
#include <src/image.h>
#include <src/projector.h>
#include <src/metadata_table.h>
#include <src/fftw.h>
#include <src/ctf.h>
#include <src/time.h>
#include <src/symmetries.h>



class tiltpair_plot_parameters {

    public:
    FileName fn_unt, fn_til, fn_eps, fn_sym;
    MetaDataTable MDu, MDt;
    RFLOAT exp_tilt, exp_beta, dist_from_alpha, dist_from_tilt, plot_max_tilt, plot_spot_radius;
    // I/O Parser
    IOParser parser;
    SymList SL;
    std::ofstream fh_eps;


    void usage() {
        parser.writeUsage(std::cerr);
    }

    void read(int argc, char **argv) {
        parser.setCommandLine(argc, argv);

        int general_section = parser.addSection("General options");
        fn_unt = parser.getOption("--u", "Input STAR file with untilted particles");
        fn_til = parser.getOption("--t", "Input STAR file with tilted particles");
        fn_eps = parser.getOption("--o", "Output EPS file ", "tiltpair.eps");
        fn_sym = parser.getOption("--sym", "Symmetry point group", "C1");
        exp_tilt = textToFloat(parser.getOption("--exp_tilt", "Choose symmetry operator that gives tilt angle closest to this value", "0."));
        exp_beta = textToFloat(parser.getOption("--exp_beta", "Choose symmetry operator that gives beta angle closest to this value", "0."));
        dist_from_alpha = textToFloat(parser.getOption("--dist_from_alpha", "Direction (alpha angle) of tilt axis from which to calculate distance", "0."));
        dist_from_tilt = textToFloat(parser.getOption("--dist_from_tilt", "Tilt angle from which to calculate distance", "0."));
        plot_max_tilt = textToFloat(parser.getOption("--max_tilt", "Maximum tilt angle to plot in the EPS file", "90."));
        plot_spot_radius = textToInteger(parser.getOption("--spot_radius", "Radius in pixels of the spots in the tiltpair plot", "3"));

        // Check for errors in the command-line option
        if (parser.checkForErrors())
            REPORT_ERROR("Errors encountered on the command line, exiting...");
    }

    void initialise() {

        // Get the MDs for both untilted and tilted particles
        MDu.read(fn_unt);
        MDt.read(fn_til);
        if (MDu.size() != MDt.size())
            REPORT_ERROR("Tiltpair plot ERROR: untilted and tilted STAR files have unequal number of entries.");

        // Get the symmetry point group
        int pgGroup, pgOrder;
        SL.isSymmetryGroup(fn_sym.removeDirectories(), pgGroup, pgOrder);
        SL.read_sym_file(fn_sym);

        // Make postscript header
        fh_eps.open(fn_eps.c_str(), std::ios::out);
        if (!fh_eps)
            REPORT_ERROR("Tiltpair plot ERROR: Cannot open " + fn_eps + " for output");

        fh_eps << "%%!PS-Adobe-2.0\n";
        fh_eps << "%% Creator: Tilt pair analysis \n";
        fh_eps << "%% Pages: 1\n";
        fh_eps << "0 setgray\n";
        fh_eps << "0.1 setlinewidth\n";
        // Draw circles on postscript: 250pixels=plot_max_tilt
        fh_eps << "300 400 83 0 360 arc closepath stroke\n";
        fh_eps << "300 400 167 0 360 arc closepath stroke\n";
        fh_eps << "300 400 250 0 360 arc closepath stroke\n";
        fh_eps << "300 150 newpath moveto 300 650 lineto stroke\n";
        fh_eps << "50 400 newpath moveto 550 400 lineto stroke\n";
    }

    void add_to_postscript(RFLOAT tilt_angle, RFLOAT alpha, RFLOAT beta) {

        RFLOAT rr = tilt_angle / plot_max_tilt * 250;
        // SINCOS?
        RFLOAT x = 300.0 + rr * cos(radians(alpha));
        RFLOAT y = 400.0 + rr * sin(radians(alpha));
        RFLOAT r, g, b;
        value_to_redblue_scale(abs(90.0 - beta), 0.0, 90.0, r, g, b);
        fh_eps << x << " " << y << " " << plot_spot_radius << " 0 360 arc closepath " << r << " " << g << " " << b << " setrgbcolor fill stroke\n";
    }

    void value_to_redblue_scale(RFLOAT val, RFLOAT minF, RFLOAT maxF, RFLOAT &r, RFLOAT &g, RFLOAT &b) {
        RFLOAT half = (maxF - minF) / 2.0;
        if (val < half) {
            r = val / half;
            b = 1.0;
        } else {
            b = (maxF - val) / half;
            r = 1.0;
        }
        g = 0.0;
    }

    RFLOAT check_symmetries(
        RFLOAT  rot1, RFLOAT  tilt1, RFLOAT  psi1,
        RFLOAT &rot2, RFLOAT &tilt2, RFLOAT &psi2
    ) {

        int imax = SL.SymsNo() + 1;
        Matrix2D<RFLOAT> L(4, 4), R(4, 4);  // A matrix from the list
        RFLOAT best_ang_dist = 3600;
        RFLOAT best_rot2, best_tilt2, best_psi2;

        for (int i = 0; i < imax; i++) {
            RFLOAT rot2p, tilt2p, psi2p;
            if (i == 0) {
                rot2p  = rot2;
                tilt2p = tilt2;
                psi2p  = psi2;
            } else {
                SL.get_matrices(i - 1, L, R);
                L.resize(3, 3); // Erase last row and column
                R.resize(3, 3); // as only the relative orientation is useful and not the translation
                angles_t angles = Euler::apply_transf(L, R, rot2, tilt2, psi2);
                rot2p = angles.rot; tilt2p = angles.tilt; psi2p = angles.psi;
            }

            RFLOAT ang_dist = check_tilt_pairs(rot1, tilt1, psi1, rot2p, tilt2p, psi2p);

            if (ang_dist < best_ang_dist) {
                best_ang_dist = ang_dist;
                best_rot2 = rot2p;
                best_tilt2 = tilt2p;
                best_psi2 = psi2p;
            }
        }

        rot2  = best_rot2;
        tilt2 = best_tilt2;
        psi2  = best_psi2;

        return best_ang_dist;
    }

    RFLOAT check_tilt_pairs(
        RFLOAT rot1, RFLOAT tilt1, RFLOAT psi1,
        RFLOAT &alpha, RFLOAT &tilt_angle, RFLOAT &beta
    ) {
        // Transformation matrices
        RFLOAT rot2 = alpha, tilt2 = tilt_angle, psi2 = beta;

        // Calculate the transformation from one setting to the second one.
        Matrix2D<RFLOAT> E1 = Euler::angles2matrix(psi1, tilt1, rot1);
        Matrix2D<RFLOAT> E2 = Euler::angles2matrix(psi2, tilt2, rot2).matmul(E1.inv());

        // Get the tilt angle (and its sine)
        RFLOAT ah = (E2(0, 0) + E2(1, 1) + E2(2, 2) - 1.0) / 2.0;
        if (Xmipp::gt(std::abs(ah), 1.0)) REPORT_ERROR("BUG: ah > 1");  // std::abs needed because abs comes from elsewhere
        RFLOAT tilt_angle_radians = acos(ah);
        tilt_angle = degrees(tilt_angle_radians);
        RFLOAT sine_tilt_angle = 2.0 * sin(tilt_angle_radians);

        Matrix1D<RFLOAT> axis;
        // Get the tilt axis direction in angles alpha and beta
        if (sine_tilt_angle > Xmipp::epsilon) {
            axis = Matrix1D<RFLOAT>({
                (E2(2, 1) - E2(1, 2)) / sine_tilt_angle,
                (E2(0, 2) - E2(2, 0)) / sine_tilt_angle,
                (E2(1, 0) - E2(0, 1)) / sine_tilt_angle});
        } else {
            axis = Matrix1D<RFLOAT>({0.0, 0.0, 1.0});
        }

        // Apply E1.inv() to the axis to get everyone in the same coordinate system again
        axis = matmul(E1.inv(), axis);

        // Convert to alpha and beta angle
        Euler::direction2angles(axis, alpha, beta);

        // Enforce positive beta: choose the other Euler angle combination to express the same direction
        if (beta < 0.0) {
            beta = -beta;
            alpha += 180.0;
        }

        // Let alpha go from 0 to 360 degrees
        alpha = wrap(alpha, 0.0, 360.0);

        // Return the value that needs to be optimized
        RFLOAT minimizer = 0.0;
        if (exp_beta < 999.0) { minimizer  = abs(beta       - exp_beta); }
        if (exp_tilt < 999.0) { minimizer += abs(tilt_angle - exp_tilt); }

        return minimizer;
    }

    void run() {
        for (long int iline : MDu) {

            // Read input data
            // RFLOAT best_tilt, best_alpha, best_beta;

            RFLOAT rot1  = MDu.getValue<RFLOAT>(EMDL::ORIENT_ROT,  iline);
            RFLOAT rot2  = MDt.getValue<RFLOAT>(EMDL::ORIENT_ROT,  iline);
            RFLOAT tilt1 = MDu.getValue<RFLOAT>(EMDL::ORIENT_TILT, iline);
            RFLOAT tilt2 = MDt.getValue<RFLOAT>(EMDL::ORIENT_TILT, iline);
            RFLOAT psi1  = MDu.getValue<RFLOAT>(EMDL::ORIENT_PSI,  iline);
            RFLOAT psi2  = MDt.getValue<RFLOAT>(EMDL::ORIENT_PSI,  iline);

            // Bring both angles to a normalized set
            rot1  = wrap(rot1,  -180.0, +180.0);
            tilt1 = wrap(tilt1, -180.0, +180.0);
            psi1  = wrap(psi1,  -180.0, +180.0);
            rot2  = wrap(rot2,  -180.0, +180.0);
            tilt2 = wrap(tilt2, -180.0, +180.0);
            psi2  = wrap(psi2,  -180.0, +180.0);

            // Apply rotations to find the minimum distance angles
            RFLOAT rot2p  = rot2;
            RFLOAT tilt2p = tilt2;
            RFLOAT psi2p  = psi2;
            RFLOAT distp  = check_symmetries(rot1, tilt1, psi1, rot2p, tilt2p, psi2p);

            // Calculate distance to user-defined point
            Matrix1D<RFLOAT> aux2 (4);
            // SINCOS?
            RFLOAT xp = dist_from_tilt * cos(radians(dist_from_alpha));
            RFLOAT yp = dist_from_tilt * sin(radians(dist_from_alpha));
            RFLOAT x = tilt2p * cos(radians(rot2p));
            RFLOAT y = tilt2p * sin(radians(rot2p));
            XX(aux2) = tilt2p;
            YY(aux2) = rot2p;
            ZZ(aux2) = psi2p;
            aux2[3] = euclid(xp - x, yp - y);
            add_to_postscript(tilt2p, rot2p, psi2p);
        }

        // Close the EPS file to write it to disk
        fh_eps << "showpage\n";
        fh_eps.close();
    }
};


int main(int argc, char *argv[]) {
    tiltpair_plot_parameters prm;

    try {
        prm.read(argc, argv);
        prm.initialise();
        prm.run();
    } catch (RelionError XE) {
        std::cerr << XE;
        // prm.usage();
        return RELION_EXIT_FAILURE;
    }
    return RELION_EXIT_SUCCESS;
}
