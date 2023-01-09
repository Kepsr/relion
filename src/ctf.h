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
/***************************************************************************
 *
 * Authors: Carlos Oscar S. Sorzano (coss@cnb.csic.es)
 *
 * Unidad de Bioinformatica of Centro Nacional de Biotecnologia , CSIC
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
 e You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
 * 02111-1307  USA
 *
 * All comments concerning this program package may be sent to the
 * e-mail address 'xmipp@cnb.csic.es'
 ***************************************************************************/

#ifndef _CTF_HH
#define _CTF_HH

#include <vector>
#include "src/complex.h"
#include "src/jaz/gravis/t2Vector.h"

class CTF {

    protected:

    // Different constants
    RFLOAT K1, K2, K3, K4, K5;

    // Astigmatism stored in symmetrical matrix form
    RFLOAT Axx, Axy, Ayy;

    // Azimuthal angle (radians)
    RFLOAT rad_azimuth;

    // defocus_average = (defocus_u + defocus_v)/2
    RFLOAT defocus_average;

    // defocus_deviation = (defocus_u - defocus_v)/2
    RFLOAT defocus_deviation;

    public:

    // Acceleration voltage (kilovolts)
    RFLOAT kV;

    // Defocus in U (in Angstroms).
    // Positive values are underfocused.
    RFLOAT DeltafU;

    // Defocus in V (in Angstroms).
    // Positive values are underfocused.
    RFLOAT DeltafV;

    // Azimuthal angle (between X and U) in degrees
    RFLOAT azimuthal_angle;

    // Electron wavelength (Angstroms)
    RFLOAT lambda;

    // Radius of the aperture (in microns)
    // RFLOAT aperture;

    // Spherical aberration (in millimeters).
    // Typical value 5.6
    RFLOAT Cs;

    // Chromatic aberration (in millimeters).
    // Typical value 2
    RFLOAT Ca;

    // Mean energy loss (in eV) due to interaction with sample.
    // Typical value 1
    RFLOAT espr;

    // Objective lens stability (deltaI/I) (ppm).
    // Typical value 1
    RFLOAT ispr;

    // Convergence cone semiangle (in mrad).
    // Typical value 0.5
    RFLOAT alpha;

    // Longitudinal mechanical displacement (Angstrom). Typical value 100
    RFLOAT DeltaF;

    // Transverse mechanical displacement (Angstrom). Typical value 3
    RFLOAT DeltaR;

    // Amplitude contrast. Typical values 0.07 for cryo, 0.2 for negative stain
    RFLOAT Q0;

    // B-factor fall-off
    RFLOAT Bfac;

    // Overall scale-factor of CTF
    RFLOAT scale;

    // Phase-shift from a phase-plate (in rad)
    RFLOAT phase_shift;

    /** Empty constructor. */
    CTF():
    kV(200), DeltafU(0), DeltafV(0), azimuthal_angle(0), phase_shift(0),
    Cs(0), Bfac(0), Q0(0), scale(1)
    {}

    CTF(
        RFLOAT defU, RFLOAT defV, RFLOAT defAng,
        RFLOAT voltage, RFLOAT Cs, RFLOAT Q0,
        RFLOAT Bfac = 0.0, RFLOAT scale = 1.0, RFLOAT phase_shift = 0.0
    ):  DeltafU(defU), DeltafV(defV), azimuthal_angle(defAng), 
        kV(voltage), Cs(Cs), Q0(Q0),
        scale(scale), phase_shift(phase_shift), Bfac(Bfac) {
        initialise();
    }

    void initialise();

    // Compute CTF at (U,V). Continuous frequencies
    inline RFLOAT operator () (
        RFLOAT X, RFLOAT Y,
        bool do_only_flip_phases = false,
        bool do_intact_until_first_peak = false, bool do_damping = true,
        double gammaOffset = 0.0, bool do_intact_after_first_peak = false
    ) const {

        const RFLOAT u2 = X * X + Y * Y;
        const RFLOAT u4 = u2 * u2;

        // if (u2>=ua2) return 0;
        // RFLOAT deltaf = getDeltaF(X, Y);
        // RFLOAT gamma = K1 * deltaf * u2 + K2 * u4 - K5 - K3 + gammaOffset;
        RFLOAT gamma = K1 * astigDefocus(X, Y) + K2 * u4 - K5 - K3 + gammaOffset;
        // Quadratic: xx + 2xy + yy

        RFLOAT CT = do_intact_until_first_peak && 2 * abs(gamma) < PI ||
                    do_intact_after_first_peak && 2 * abs(gamma) > PI ? 1 : -sin(gamma);

        if (do_damping)
            CT *= exp(K4 * u2);  // B-factor decay (K4 = -Bfac/4);

        if (do_only_flip_phases)
            CT = float(CT >= 0);

        CT *= scale;

        // SHWS 25-2-2019: testing a new idea to improve code stability
        // In order to prevent division by zero in GPU code, 
        // don't allow very small CTF values.
        if (fabs(CT) < 1e-8)
            CT = copysign(1e-8, CT);

        return CT;
    }

    RFLOAT getGamma(RFLOAT X, RFLOAT Y) const;

    // compute the local frequency of the ctf
    // (i.e. the radial slope of 'double gamma' in operator ())
    // -- deprecated, use getGammaGrad().length()
    RFLOAT getCtfFreq(RFLOAT X, RFLOAT Y);

    gravis::t2Vector<RFLOAT> getGammaGrad(RFLOAT X, RFLOAT Y) const;

    inline Complex getCTFP(RFLOAT X, RFLOAT Y, double gammaOffset = 0.0) const {
        const RFLOAT u2 = X * X + Y * Y;
        const RFLOAT u4 = u2 * u2;
        const RFLOAT gamma = K1 * astigDefocus(X, Y) + K2 * u4 - K5 - K3 + gammaOffset + PI / 2.0;
        return Complex::unit(gamma);
    }

    // Compute Deltaf at a given direction (no longer used by operator ())
    inline RFLOAT getDeltaF(RFLOAT X, RFLOAT Y) const {
        if (abs(X) < Xmipp::epsilon<RFLOAT>() && abs(Y) < Xmipp::epsilon<RFLOAT>())
            return 0;

        RFLOAT ellipsoid_ang = atan2(Y, X) - rad_azimuth;
        /**
         * For a derivation of this formula,
         * see Principles of Electron Optics p. 1380.
         * In particular, term defocus and twofold axial astigmatism
         * take into account that a1 and a2 are
         * the coefficient of the Zernike polynomials difference of defocus
         * at 0 and at 45 degrees.
         * In this case, a2 = 0.
         */
        return defocus_average + defocus_deviation * cos(2 * ellipsoid_ang);

    }

    std::vector<double> getK();

    double getAxx();

    double getAxy();

    double getAyy();

    inline RFLOAT getLambda() const { return lambda; }

    inline double astigDefocus(RFLOAT x, RFLOAT y) const {
        return Axx * x * x + 2 * Axy * x * y + Ayy * y * y;
    }

};

#endif
