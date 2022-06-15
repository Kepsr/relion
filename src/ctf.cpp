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
 * Authors:     Carlos Oscar S. Sorzano (coss@cnb.csic.es)
 *
 * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
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
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
 * 02111-1307  USA
 *
 *  All comments concerning this program package may be sent to the
 *  e-mail address 'xmipp@cnb.csic.es'
 ***************************************************************************/

#include "src/ctf.h"
#include "src/args.h"
#include "src/fftw.h"
#include <src/jaz/obs_model.h>
#include <src/jaz/gravis/t2Matrix.h>
#include "src/numerical_recipes.h" // For Pythag

using namespace gravis;

static void assert_square(int orixdim, int oriydim) {
    if (orixdim != oriydim) {
        REPORT_ERROR_STR(
            "CTF::getFftwImage: currently, symmetric aberrations are supported "
            << "only for square images.\n"
        );
    }
}

void CTF::setValues(
    RFLOAT _defU,    RFLOAT _defV,  RFLOAT _defAng,
    RFLOAT _voltage, RFLOAT _Cs,    RFLOAT _Q0,
    RFLOAT _Bfac,    RFLOAT _scale, RFLOAT _phase_shift
) {
    kV              = _voltage;
    DeltafU         = _defU;
    DeltafV         = _defV;
    azimuthal_angle = _defAng;
    Cs              = _Cs;
    Bfac            = _Bfac;
    scale           = _scale;
    Q0              = _Q0;
    phase_shift     = _phase_shift;

    initialise();
}

void CTF::initialise() {
    // Change units
    RFLOAT local_Cs = Cs * 1e7;
    RFLOAT local_kV = kV * 1e3;
    rad_azimuth = radians(azimuthal_angle);

    // Average focus and deviation
    defocus_average   = -(DeltafU + DeltafV) * 0.5;
    defocus_deviation = -(DeltafU - DeltafV) * 0.5;

    // lambda=h/sqrt(2*m*e*kV)
    //    h: Planck constant
    //    m: electron mass
    //    e: electron charge
    // lambda=0.387832/sqrt(kV * (1.+0.000978466*kV)); // Hewz: Angstroms
    // lambda=h/sqrt(2*m*e*kV)
    lambda = 12.2643247 / sqrt(local_kV * (1.0 + local_kV * 0.978466e-6));
    // See http://en.wikipedia.org/wiki/Electron_diffraction

    // Helpful constants
    // ICE: X(u)=-PI/2*deltaf(u)*lambda*u^2+PI/2*Cs*lambda^3*u^4
    //          = K1*deltaf(u)*u^2         +K2*u^4
    K1 = PI / 2 * 2 * lambda;
    K2 = PI / 2 * local_Cs * lambda * lambda * lambda;
    K3 = atan(Q0 / sqrt(1 - Q0 * Q0));
    K4 = -Bfac / 4.0;

    // Phase shift in radians
    K5 = radians(phase_shift);

    if (Q0 < 0.0 || Q0 > 1.0)
        REPORT_ERROR("CTF::initialise ERROR: AmplitudeContrast Q0 cannot be smaller than zero or larger than one!");

    if (abs(DeltafU) < 1e-6 && abs(DeltafV) < 1e-6 && abs(Q0) < 1e-6 && abs(Cs) < 1e-6)
        REPORT_ERROR("CTF::initialise: ERROR: CTF initialises to all-zero values. Was a correct STAR file provided?");

    // express astigmatism as a bilinear form:

    const double sin_az = sin(rad_azimuth);
    const double cos_az = cos(rad_azimuth);

    d2Matrix Q (cos_az, +sin_az, -sin_az, cos_az);
    d2Matrix Qt(cos_az, -sin_az, +sin_az, cos_az);
    d2Matrix D (-DeltafU, 0.0, 0.0, -DeltafV);

    d2Matrix A = Qt * D * Q;

    Axx = A(0, 0);
    Axy = A(0, 1);
    Ayy = A(1, 1);
}

RFLOAT CTF::getGamma(RFLOAT X, RFLOAT Y) const {

    RFLOAT u2 = X * X + Y * Y;
    RFLOAT u4 = u2 * u2;

    return K1 * (Axx * X * X + 2.0 * Axy * X * Y + Ayy * Y * Y) + K2 * u4 - K5 - K3;
}

RFLOAT CTF::getCtfFreq(RFLOAT X, RFLOAT Y) {
    RFLOAT u2 = X * X + Y * Y;
    RFLOAT u = sqrt(u2);

    RFLOAT deltaf = getDeltaF(X, Y);

    return 2.0 * K1 * deltaf * u + 4.0 * K2 * u * u * u;
}

t2Vector<RFLOAT> CTF::getGammaGrad(RFLOAT X, RFLOAT Y) const {

    RFLOAT u2 = X * X + Y * Y;
    // RFLOAT u4 = u2 * u2;

    // u4 = (X² + Y²)²
    // du4/dx = 2 (X² + Y²) 2 X = 4 (X³ + XY²) = 4 u2 X

    return t2Vector<RFLOAT>(
        2.0 * (K1 * Axx * X + K1 * Axy * Y + 2.0 * K2 * u2 * X),
        2.0 * (K1 * Ayy * Y + K1 * Axy * X + 2.0 * K2 * u2 * Y)
    );
}

// Generate a complete CTF Image ----------------------------------------------
MultidimArray<RFLOAT> CTF::getFftwImage(
    long int Xdim, long int Ydim, int orixdim, int oriydim,
    RFLOAT angpix, ObservationModel *obsModel, int opticsGroup,
    bool do_abs, bool do_only_flip_phases, bool do_intact_until_first_peak,
    bool do_damping, bool do_ctf_padding, bool do_intact_after_first_peak
) const {
    MultidimArray<RFLOAT> result(Xdim, Ydim);
    // Boxing the particle in a small box from the whole micrograph leads to loss of delocalised information (or aliaising in the CTF)
    // Here, calculate the CTF in a 2x larger box to support finer oscillations,
    // and then rescale the large CTF to simulate the effect of the windowing operation
    if (do_ctf_padding) {
        bool ctf_premultiplied = obsModel && obsModel->getCtfPremultiplied(opticsGroup);

        // two-fold padding, increased to 4-fold padding for pre-multiplied CTFs
        int orixdim_pad = 2 * orixdim;
        int oriydim_pad = 2 * oriydim;
        // TODO: Such a big box may not really be necessary...
        if (ctf_premultiplied) {
            orixdim_pad *= 2;
            oriydim_pad *= 2;
        }

        MultidimArray<RFLOAT> Fctf = getFftwImage(
            oriydim_pad / 2 + 1, oriydim_pad, orixdim_pad, oriydim_pad, angpix,
            obsModel,
            do_abs, do_only_flip_phases, do_intact_until_first_peak, do_damping, false, do_intact_after_first_peak
        );

        // From half to whole
        MultidimArray<RFLOAT> Mctf(oriydim_pad, orixdim_pad);
        Mctf.setXmippOrigin();
        for (int j = 0 ; j < Ysize(Fctf); j++) {
            // Don't take the middle row of the half-transform
            if (j != Ysize(Fctf) / 2) {
                int jp = j < Xsize(Fctf) ? j : j - Ysize(Fctf);
                // Don't take the last column from the half-transform
                for (int i = 0; i < Xsize(Fctf) - 1; i++) {
                    RFLOAT fctfij = direct::elem(Fctf, i, j);
                    if (ctf_premultiplied) {
                        Mctf.elem(+i, +jp) = fctfij * fctfij;
                        Mctf.elem(-i, -jp) = fctfij * fctfij;
                    } else {
                        Mctf.elem(+i, +jp) = fctfij;
                        Mctf.elem(-i, -jp) = fctfij;
                    }
                }
            }
        }

        resizeMap(Mctf, orixdim);

        Mctf.setXmippOrigin();
        // From whole to half
        for (int j = 0 ; j < Ysize(result); j++) {
            // Don't take the middle row of the half-transform
            if (j != Ysize(result) / 2) {
                int jp = j < Xsize(result) ? j : j - Ysize(result);
                // Don't take the last column from the half-transform
                for (int i = 0; i < Xsize(result) - 1; i++) {
                    // Make just one lookup on Mctf.data
                    RFLOAT mctfipj = Mctf.elem(i, jp);
                    if (ctf_premultiplied) {
                        // Constrain result[i, j] to the interval [0.0, 1.0]
                        direct::elem(result, i, j) =
                            mctfipj < 0.0 ? 0.0 :
                            mctfipj > 1.0 ? 1.0 :
                            sqrt(mctfipj);
                    } else {
                        direct::elem(result, i, j) = mctfipj;
                    }
                }
            }
        }
    } else {
        RFLOAT xs = (RFLOAT) orixdim * angpix;
        RFLOAT ys = (RFLOAT) oriydim * angpix;

        if (obsModel && obsModel->hasEvenZernike) {

            assert_square(orixdim, oriydim);

            if (obsModel->getBoxSize(opticsGroup) != orixdim) {
                REPORT_ERROR_STR(
                    "CTF::getFftwImage: requested output image size "
                    << orixdim
                    << " is not consistent with that in the optics group table "
                    << obsModel->getBoxSize(opticsGroup) << "\n"
                );
            }

            if (fabs(obsModel->getPixelSize(opticsGroup) - angpix) > 1e-4) {
                REPORT_ERROR_STR(
                    "CTF::getFftwImage: requested pixel size "
                    << angpix
                    << " is not consistent with that in the optics group table "
                    << obsModel->getPixelSize(opticsGroup) << "\n"
                );
            }

            const Image<RFLOAT>& gammaOffset = obsModel->getGammaOffset(opticsGroup, oriydim);

            for (int j = 0; j < result.ydim; j++)
            for (int i = 0; i < result.xdim; i++) {
                RFLOAT x = i / xs;
                RFLOAT y = (j <= result.ydim / 2 ? j : j - result.ydim) / ys;

                const int x0 = i;
                const int y0 = j <= result.ydim / 2 ? j : gammaOffset.data.ydim + j - result.ydim;

                if (obsModel) obsModel->magnify(x, y, obsModel->getMagMatrix(opticsGroup));
                RFLOAT t = getCTF(
                    x, y,
                    do_only_flip_phases, do_intact_until_first_peak,
                    do_damping, gammaOffset(y0, x0), do_intact_after_first_peak
                );
                direct::elem(result, i, j) = do_abs ? abs(t) : t;
            }
        } else {
            FOR_ALL_ELEMENTS_IN_FFTW_TRANSFORM2D(result) {
                RFLOAT x = (RFLOAT) ip / xs;
                RFLOAT y = (RFLOAT) jp / ys;

                if (obsModel) obsModel->magnify(x, y, obsModel->getMagMatrix(opticsGroup));
                RFLOAT t = getCTF(
                    x, y,
                    do_only_flip_phases, do_intact_until_first_peak,
                    do_damping, 0.0, do_intact_after_first_peak
                );
                direct::elem(result, i, j) = do_abs ? abs(t) : t;
            }
        }
    }
    return result;
}

// Generate a complete CTFP (complex) image (with sector along angle) ---------
MultidimArray<Complex> CTF::getCTFPImage(
    long int Xdim, long int Ydim, int orixdim, int oriydim, RFLOAT angpix,
    ObservationModel *obsModel, int opticsGroup,
    bool is_positive, float angle
) {
    if (angle < 0.0 || angle >= 360.0)
        REPORT_ERROR("CTF::getCTFPImage: angle should be in the interval [0,360)");

    MultidimArray<Complex> result(Xdim, Ydim);

    // Flip angles greater than 180 degrees
    if (angle >= 180.0) {
        angle -= 180.0;
        is_positive = !is_positive;
    }

    float anglerad = radians(angle);

    RFLOAT xs = (RFLOAT) orixdim * angpix;
    RFLOAT ys = (RFLOAT) oriydim * angpix;

    if (obsModel && obsModel->hasEvenZernike) {

        assert_square(orixdim, oriydim);

        const Image<RFLOAT> &gammaOffset = obsModel->getGammaOffset(opticsGroup, oriydim);

        if (gammaOffset.data.xdim < result.xdim || gammaOffset.data.ydim < result.ydim) {
            REPORT_ERROR_STR(
                "CTF::getFftwImage: size requested for output image "
                << "is greater than size of original image: "
                << result.xdim << "×" << result.ydim << " was requested, "
                << "but only " << gammaOffset.data.xdim << "×" << gammaOffset.data.ydim << " is available\n"
            );
        }

        // Why do we have i <= Ysize(result) / 2 here, but i < Xsize(result) below?
        for (int j = 0, jp = 0; j < Ysize(result); j++, jp = j <= Ysize(result) / 2 ? j : j - Ysize(result))
        for (int i = 0, ip = 0; i < Xsize(result); i++, ip = i) {
            RFLOAT x = (RFLOAT) ip / xs;
            RFLOAT y = (RFLOAT) jp / ys;
            RFLOAT myangle = x * x + y * y > 0 ? acos(y / Pythag(x, y)) : 0; // dot-product with Y-axis: (0, 1)
            const int x0 = i;
            const int y0 = j <= Ysize(result) / 2 ? jp : Ysize(gammaOffset.data) + jp;

            if (obsModel) obsModel->magnify(x, y, obsModel->getMagMatrix(opticsGroup));
            Complex ctfp = getCTFP(x, y, gammaOffset(y0, x0));
            if ((myangle >= anglerad) != is_positive) { ctfp = ctfp.conj(); }
            direct::elem(result, i, j) = ctfp;
        }
    } else {
        // FOR_ALL_ELEMENTS_IN_FFTW_TRANSFORM2D(result) {
        for (long int j = 0, jp = 0; j < Ysize(result); j++, jp = j < Xsize(result) ? j : j - Ysize(result)) \
        for (long int i = 0, ip = 0; i < Xsize(result); i++, ip = i) {
            // If i < Xsize(result), ip = i. Else, ip = i - Ysize(result).
            RFLOAT x = (RFLOAT) ip / xs;
            RFLOAT y = (RFLOAT) jp / ys;
            RFLOAT myangle = x * x + y * y > 0 ? acos(y / Pythag(x, y)) : 0; // dot-product with Y-axis: (0, 1)

            if (obsModel) obsModel->magnify(x, y, obsModel->getMagMatrix(opticsGroup));
            Complex ctfp = getCTFP(x, y);
            if ((myangle >= anglerad) != is_positive) { ctfp = ctfp.conj(); }
            direct::elem(result, i, j) = ctfp;
        }
    }
    // Special line along the vertical (Y-)axis, where FFTW stores both Friedel mates and Friedel symmetry needs to remain
    if (angle == 0.0) {
        for (int j = result.ydim / 2 + 1; j < result.ydim; j++) {
            direct::elem(result, 0, j) = conj(direct::elem(result, 0, result.ydim - j));
        }
    }
    return result;
}

MultidimArray<RFLOAT> CTF::getCenteredImage(
    long int Xdim, long int Ydim,
    RFLOAT Tm, ObservationModel *obsModel, int opticsGroup,
    bool do_abs, bool do_only_flip_phases, bool do_intact_until_first_peak,
    bool do_damping, bool do_intact_after_first_peak
) {
    MultidimArray<RFLOAT> result(Xdim, Ydim);
    result.setXmippOrigin();
    RFLOAT xs = (RFLOAT) Xsize(result) * Tm;
    RFLOAT ys = (RFLOAT) Ysize(result) * Tm;

    // Maybe use the absolute value
    auto f = do_abs ? [] (RFLOAT x) -> RFLOAT { return abs(x); }
                    : [] (RFLOAT x) -> RFLOAT { return x; };


    // Maybe apply magnification
    auto g = obsModel ? [] (ObservationModel *obsModel, int opticsGroup, RFLOAT &x, RFLOAT &y) -> void { obsModel->magnify(x, y, obsModel->getMagMatrix(opticsGroup)); }
                      : [] (ObservationModel *obsModel, int opticsGroup, RFLOAT &x, RFLOAT &y) -> void {};

    FOR_ALL_ELEMENTS_IN_ARRAY2D(result) {
        RFLOAT x = (RFLOAT) i / xs;
        RFLOAT y = (RFLOAT) j / ys;
        g(obsModel, opticsGroup, x, y);
        RFLOAT t = getCTF(
            x, y,
            do_only_flip_phases, do_intact_until_first_peak,
            do_damping, 0.0, do_intact_after_first_peak
        );
        result.elem(i, j) = f(t);
    }
    return result;
}

void CTF::get1DProfile(
    MultidimArray<RFLOAT> &result, RFLOAT angle, RFLOAT Tm,
    ObservationModel *obsModel, int opticsGroup,
    bool do_abs, bool do_only_flip_phases, bool do_intact_until_first_peak,
    bool do_damping, bool do_intact_after_first_peak
) {

    result.setXmippOrigin();
    RFLOAT xs = (RFLOAT) Xsize(result) * Tm; // assuming result is at the image size!

    auto f = do_abs ? [] (RFLOAT x) -> RFLOAT { return abs(x); }
                    : [] (RFLOAT x) -> RFLOAT { return x; };

    auto g = obsModel ? [] (ObservationModel *obsModel, int opticsGroup, RFLOAT x, RFLOAT y) -> void { obsModel->magnify(x, y, obsModel->getMagMatrix(opticsGroup)); }
                      : [] (ObservationModel *obsModel, int opticsGroup, RFLOAT x, RFLOAT y) -> void {};

    for (int i = Xinit(result); i <= Xlast(result); i++) {
        RFLOAT x = (RFLOAT) i * cos(radians(angle)) / xs;
        RFLOAT y = (RFLOAT) i * sin(radians(angle)) / xs;
        g(obsModel, opticsGroup, x, y);
        RFLOAT t = getCTF(
            x, y,
            do_only_flip_phases, do_intact_until_first_peak,
            do_damping, 0.0, do_intact_after_first_peak
        );
        result.elem(i) = f(t);
    }
}

void CTF::applyWeightEwaldSphereCurvature(
    MultidimArray<RFLOAT> &result, int orixdim, int oriydim,
    RFLOAT angpix, ObservationModel *obsModel, int opticsGroup, RFLOAT particle_diameter
) {
    RFLOAT xs = (RFLOAT) orixdim * angpix;
    RFLOAT ys = (RFLOAT) oriydim * angpix;

    Matrix2D<RFLOAT> M = obsModel && obsModel->hasMagMatrices ?
        obsModel->getMagMatrix(opticsGroup) :
        Matrix2D<RFLOAT>::identity(2);

    FOR_ALL_ELEMENTS_IN_FFTW_TRANSFORM2D(result) {
        RFLOAT xu = (RFLOAT) ip / xs;
        RFLOAT yu = (RFLOAT) jp / ys;

        RFLOAT x = M(0, 0) * xu + M(0, 1) * yu;
        RFLOAT y = M(1, 0) * xu + M(1, 1) * yu;

        const RFLOAT astigDefocus = Axx * x * x + 2.0 * Axy * x * y + Ayy * y * y;
        RFLOAT u2 = x * x + y * y;
        RFLOAT u4 = u2 * u2;
        RFLOAT gamma = K1 * astigDefocus + K2 * u4 - K5 - K3;

        RFLOAT deltaf = u2 > 0.0 ? std::abs(astigDefocus / u2) : 0.0;
        RFLOAT inv_d = sqrt(u2);
        RFLOAT aux = 2.0 * deltaf * lambda * inv_d / particle_diameter;
        RFLOAT A = aux > 1.0 ? 0.0 : (acos(aux) - aux * sqrt(1 - aux * aux)) * 2.0 / PI;
        // sin(acos(x)) is almost 50% slower than sqrt(1 - x * x)

        direct::elem(result, i, j) = 0.5 * (A * (2.0 * fabs(sin(gamma)) - 1.0) + 1.0);
        // Within RELION, sin(chi) is used rather than 2 * sin(chi).
        // Hence the 0.5 above to keep everything on the same scale.
    }
}

void CTF::applyWeightEwaldSphereCurvature_new(
    MultidimArray<RFLOAT>& result, int orixdim, int oriydim,
    RFLOAT angpix, ObservationModel *obsModel, int opticsGroup, RFLOAT particle_diameter
) {
    const int s = oriydim;
    const int half_s = s / 2 + 1;
    const double as = angpix * s;
    const double Dpx = particle_diameter / angpix;

    for (int yi = 0; yi <      s;  yi++)
    for (int xi = 0; xi < half_s; xi++) {
        double x = xi / as;
        double y = yi < half_s ? yi / as : (yi - s) / as;

        if (obsModel) obsModel->magnify(x, y, obsModel->getMagMatrix(opticsGroup));

        // shift of this frequency resulting from CTF:
        const t2Vector<RFLOAT> shift2D = RFLOAT(1.0 / (2 * angpix * PI)) * getGammaGrad(x, y);
        const double shift1D = 2.0 * shift2D.length();

        // Let two circles of equal size intersect at exactly two points.
        // Let alpha be the angle between
        // one intersection point,
        // the centre of either circle,
        // and the other intersection point.
        const double alpha = shift1D > Dpx ? 0.0 : 2.0 * acos(shift1D / Dpx);

        // Then the area of intersection between the two circles,
        // divided by the area of either circle will be:
        RFLOAT A = alpha == 0.0 ? 0.0 : (alpha - sin(alpha)) / PI;

        // abs. value of CTFR (no damping):
        const double ctf_val = getCTF(x, y, true, false, false, false, 0.0);

        direct::elem(result, xi, yi) = 0.5 * (A * (2.0 * ctf_val - 1.0) + 1.0);
        // Within RELION, sin(chi) is used rather than 2 * sin(chi).
        // Hence the 0.5 above to keep everything on the same scale.
    }
}

void CTF::applyWeightEwaldSphereCurvature_noAniso(
    MultidimArray <RFLOAT> &result, int orixdim, int oriydim,
    RFLOAT angpix, ObservationModel *obsModel, int opticsGroup, RFLOAT particle_diameter
) {
    RFLOAT xs = (RFLOAT) orixdim * angpix;
    RFLOAT ys = (RFLOAT) oriydim * angpix;
    FOR_ALL_ELEMENTS_IN_FFTW_TRANSFORM2D(result) {
        RFLOAT x = (RFLOAT) ip / xs;
        RFLOAT y = (RFLOAT) jp / ys;
        RFLOAT deltaf = fabs(getDeltaF(x, y));
        RFLOAT inv_d = Pythag(x, y);
        RFLOAT aux = 2.0 * deltaf * lambda * inv_d / particle_diameter;
        RFLOAT A = aux > 1.0 ? 0.0 : (acos(aux) - aux * sqrt(1 - aux * aux)) * 2.0 / PI;
        if (obsModel) obsModel->magnify(x, y, obsModel->getMagMatrix(opticsGroup));
        direct::elem(result, i, j) = 0.5 * (A * (2.0 * fabs(getCTF(x, y)) - 1.0) + 1.0);
        // Within RELION, sin(chi) is used rather than 2 * sin(chi).
        // Hence the 0.5 above to keep everything on the same scale.
    }
}

std::vector<double> CTF::getK() {
    // offset by one to maintain indices (K[1] = K1)
    return std::vector<double>{0, K1, K2, K3, K4, K5};
}

double CTF::getAxx() {
    return Axx;
}

double CTF::getAxy() {
    return Axy;
}

double CTF::getAyy() {
    return Ayy;
}
