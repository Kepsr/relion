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
#include "src/metadata_table.h"
#include <src/jaz/obs_model.h>
#include <src/jaz/gravis/t2Matrix.h>
#include "src/numerical_recipes.h" // For Pythag

using namespace gravis;

static void ensure_square(int orixdim, int oriydim) {
    if (orixdim != oriydim) {
        REPORT_ERROR_STR(
            "CTF::getFftwImage: currently, symmetric aberrations are supported "
            << "only for square images.\n"
        );
    }
}

/* Read -------------------------------------------------------------------- */
void CTF::readByGroup(
    const MetaDataTable &partMdt, ObservationModel* obs, long int particle
) {

    opticsGroup = 
        obs == 0 ? -1 :
        partMdt.getValue<int>(EMDL::IMAGE_OPTICS_GROUP, particle) - 1;

    readValue(EMDL::CTF_VOLTAGE,       kV,              200,     particle, opticsGroup, partMdt, obs);
    readValue(EMDL::CTF_DEFOCUSU,      DeltafU,         0,       particle, opticsGroup, partMdt, obs);
    readValue(EMDL::CTF_DEFOCUSV,      DeltafV,         DeltafU, particle, opticsGroup, partMdt, obs);
    readValue(EMDL::CTF_DEFOCUS_ANGLE, azimuthal_angle, 0,       particle, opticsGroup, partMdt, obs);
    readValue(EMDL::CTF_CS,            Cs,              0,       particle, opticsGroup, partMdt, obs);
    readValue(EMDL::CTF_BFACTOR,       Bfac,            0,       particle, opticsGroup, partMdt, obs);
    readValue(EMDL::CTF_SCALEFACTOR,   scale,           1,       particle, opticsGroup, partMdt, obs);
    readValue(EMDL::CTF_Q0,            Q0,              0,       particle, opticsGroup, partMdt, obs);
    readValue(EMDL::CTF_PHASESHIFT,    phase_shift,     0,       particle, opticsGroup, partMdt, obs);

    initialise();

    obsModel = obs;
}

/// TODO: Return dest
void CTF::readValue(
    EMDL::EMDLabel label, RFLOAT& dest, RFLOAT defaultVal,
    long int particle, int opticsGroup,
    const MetaDataTable& partMdt, const ObservationModel* obs
) {
    try {
        dest = partMdt.getValue<RFLOAT>(label, particle);
    } catch (const char *errmsg) { try {
        if (opticsGroup < 0) { throw "Negative optics group!"; }
        dest = obs->opticsMdt.getValue<RFLOAT>(label, opticsGroup);
    } catch (const char *errmsg) {
        dest = defaultVal;
    } }
}

void CTF::read(const MetaDataTable &MD1, const MetaDataTable &MD2, long int objectID) {

    try {
        kV = MD1.getValue<RFLOAT>(EMDL::CTF_VOLTAGE, objectID);
    } catch (const char *errmsg) { try {
        kV = MD2.getValue<RFLOAT>(EMDL::CTF_VOLTAGE, objectID);
    } catch (const char *errmsg) {
        kV = 200;
    } }

    try {
        DeltafU = MD1.getValue<RFLOAT>(EMDL::CTF_DEFOCUSU, objectID);
    } catch (const char *errmsg) { try {
        DeltafU = MD2.getValue<RFLOAT>(EMDL::CTF_DEFOCUSU, objectID);
    } catch (const char *errmsg) {
        DeltafU = 0;
    } }

    try {
        DeltafV = MD1.getValue<RFLOAT>(EMDL::CTF_DEFOCUSV, objectID);
    } catch (const char *errmsg) { try {
        DeltafV = MD2.getValue<RFLOAT>(EMDL::CTF_DEFOCUSV, objectID);
    } catch (const char *errmsg) {
        DeltafV = DeltafU;
    } }

    try {
        azimuthal_angle = MD1.getValue<RFLOAT>(EMDL::CTF_DEFOCUS_ANGLE, objectID);
    } catch (const char *errmsg) { try {
        azimuthal_angle = MD2.getValue<RFLOAT>(EMDL::CTF_DEFOCUS_ANGLE, objectID);
    } catch (const char *errmsg) {
        azimuthal_angle = 0;
    } }

    try {
        Cs = MD1.getValue<RFLOAT>(EMDL::CTF_CS, objectID);
    } catch (const char *errmsg) { try {
        Cs = MD2.getValue<RFLOAT>(EMDL::CTF_CS, objectID);
    } catch (const char *errmsg) {
        Cs = 0;
    } }

    try {
        Bfac = MD1.getValue<RFLOAT>(EMDL::CTF_BFACTOR, objectID);
    } catch (const char *errmsg) { try {
        Bfac = MD2.getValue<RFLOAT>(EMDL::CTF_BFACTOR, objectID);
    } catch (const char *errmsg) {
        Bfac = 0;
    } }

    try {
        scale = MD1.getValue<RFLOAT>(EMDL::CTF_SCALEFACTOR, objectID);
    } catch (const char *errmsg) { try {
        scale = MD2.getValue<RFLOAT>(EMDL::CTF_SCALEFACTOR, objectID);
    } catch (const char *errmsg) {
        scale = 1;
    } }

    try {
        Q0 = MD1.getValue<RFLOAT>(EMDL::CTF_Q0, objectID);
    } catch (const char *errmsg) { try {
        Q0 = MD2.getValue<RFLOAT>(EMDL::CTF_Q0, objectID);
    } catch (const char *errmsg) {
        Q0 = 0;
    } }

    try {
        phase_shift = MD1.getValue<RFLOAT>(EMDL::CTF_PHASESHIFT, objectID);
    } catch (const char *errmsg) { try {
        phase_shift = MD2.getValue<RFLOAT>(EMDL::CTF_PHASESHIFT, objectID);
    } catch (const char *errmsg) {
        phase_shift = 0;
    } }

    initialise();
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

void CTF::setValuesByGroup(
    ObservationModel *obs, int _opticsGroup,
    RFLOAT _defU, RFLOAT _defV, RFLOAT _defAng,
    RFLOAT _Bfac, RFLOAT _scale, RFLOAT _phase_shift
) {
    opticsGroup     = _opticsGroup;

    DeltafU         = _defU;
    DeltafV         = _defV;
    azimuthal_angle = _defAng;

    Bfac            = _Bfac;
    scale           = _scale;
    phase_shift     = _phase_shift;

    kV = obs->opticsMdt.getValue<RFLOAT>(EMDL::CTF_VOLTAGE, opticsGroup);
    Cs = obs->opticsMdt.getValue<RFLOAT>(EMDL::CTF_CS,      opticsGroup);
    Q0 = obs->opticsMdt.getValue<RFLOAT>(EMDL::CTF_Q0,      opticsGroup);

    initialise();

    obsModel = obs;
}

/* Read from 1 MetaDataTable ----------------------------------------------- */
void CTF::read(const MetaDataTable &MD) {
    MetaDataTable MDempty;
    MDempty.addObject(); // add one empty object
    read(MD, MDempty);
}

/** Write to an existing object in a MetaDataTable. */
void CTF::write(MetaDataTable &MD) {
    // From version-3.1 onwards: store kV, Cs, Q0 in optics table
    // MD.setValue(EMDL::CTF_VOLTAGE, kV);
    MD.setValue(EMDL::CTF_DEFOCUSU, DeltafU);
    MD.setValue(EMDL::CTF_DEFOCUSV, DeltafV);
    MD.setValue(EMDL::CTF_DEFOCUS_ANGLE, azimuthal_angle);
    // MD.setValue(EMDL::CTF_CS, Cs);
    MD.setValue(EMDL::CTF_BFACTOR, Bfac);
    MD.setValue(EMDL::CTF_SCALEFACTOR, scale);
    MD.setValue(EMDL::CTF_PHASESHIFT, phase_shift);
    // MD.setValue(EMDL::CTF_Q0, Q0);
}

/* Write ------------------------------------------------------------------- */
void CTF::write(std::ostream &out) {
    MetaDataTable MD;
    MD.addObject();
    write(MD);
    MD.write(out);
}

/* Initialise the CTF ------------------------------------------------------ */
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
    K4 = -Bfac / 4.;

    // Phase shift in radians
    K5 = radians(phase_shift);

    if (Q0 < 0.0 || Q0 > 1.0)
        REPORT_ERROR("CTF::initialise ERROR: AmplitudeContrast Q0 cannot be smaller than zero or larger than one!");

    if (abs(DeltafU) < 1e-6 && abs(DeltafV) < 1e-6 && abs(Q0) < 1e-6 && abs(Cs) < 1e-6)
        REPORT_ERROR("CTF::initialise: ERROR: CTF initialises to all-zero values. Was a correct STAR file provided?");

    // express astigmatism as a bilinear form:

    const double sin_az = sin(rad_azimuth);
    const double cos_az = cos(rad_azimuth);

    d2Matrix Q(cos_az, sin_az, -sin_az, cos_az);
    d2Matrix Qt(cos_az, -sin_az, sin_az, cos_az);
    d2Matrix D(-DeltafU, 0.0, 0.0, -DeltafV);

    d2Matrix A = Qt * D * Q;

    Axx = A(0, 0);
    Axy = A(0, 1);
    Ayy = A(1, 1);
}

RFLOAT CTF::getGamma(RFLOAT X, RFLOAT Y) const {
    if (obsModel != 0 && obsModel->hasMagMatrices) {
        const Matrix2D<RFLOAT>& M = obsModel->getMagMatrix(opticsGroup);
        RFLOAT XX = M(0,0) * X + M(0,1) * Y;
        RFLOAT YY = M(1,0) * X + M(1,1) * Y;

        X = XX;
        Y = YY;
    }

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
    if (obsModel != 0 && obsModel->hasMagMatrices) {
        const Matrix2D<RFLOAT>& M = obsModel->getMagMatrix(opticsGroup);
        RFLOAT XX = M(0, 0) * X + M(0, 1) * Y;
        RFLOAT YY = M(1, 0) * X + M(1, 1) * Y;

        X = XX;
        Y = YY;
    }

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
void CTF::getFftwImage(
    MultidimArray<RFLOAT> &result, int orixdim, int oriydim, RFLOAT angpix,
    bool do_abs, bool do_only_flip_phases, bool do_intact_until_first_peak,
    bool do_damping, bool do_ctf_padding, bool do_intact_after_first_peak
) const {
    // Boxing the particle in a small box from the whole micrograph leads to loss of delocalised information (or aliaising in the CTF)
    // Here, calculate the CTF in a 2x larger box to support finer oscillations,
    // and then rescale the large CTF to simulate the effect of the windowing operation
    if (do_ctf_padding) {
        bool ctf_premultiplied = false;
        if (obsModel != 0) {
            ctf_premultiplied = obsModel->getCtfPremultiplied(opticsGroup);
        }

        // two-fold padding, increased to 4-fold padding for pre-multiplied CTFs
        int orixdim_pad = 2 * orixdim;
        int oriydim_pad = 2 * oriydim;
        // TODO: Such a big box may not really be necessary...
        if (ctf_premultiplied) {
            orixdim_pad *= 2;
            oriydim_pad *= 2;
        }

        MultidimArray<RFLOAT> Fctf(oriydim_pad, orixdim_pad / 2 + 1);

        getFftwImage(
            Fctf, orixdim_pad, oriydim_pad, angpix, do_abs,
            do_only_flip_phases, do_intact_until_first_peak, do_damping, false, do_intact_after_first_peak
        );

        // From half to whole
        MultidimArray<RFLOAT> Mctf(oriydim_pad, orixdim_pad);
        Mctf.setXmippOrigin();
        for (int i = 0 ; i < YSIZE(Fctf); i++) {
            // Don't take the middle row of the half-transform
            if (i != YSIZE(Fctf) / 2) {
                int ip = i < XSIZE(Fctf) ? i : i - YSIZE(Fctf);
                // Don't take the last column from the half-transform
                for (int j = 0; j < XSIZE(Fctf) - 1; j++) {
                    // Make just one lookup on Fctf.data
                    RFLOAT fctfij = DIRECT_A2D_ELEM(Fctf, i, j);
                    if (ctf_premultiplied) {
                        A2D_ELEM(Mctf,  ip,  j) = fctfij * fctfij;
                        A2D_ELEM(Mctf, -ip, -j) = fctfij * fctfij;
                    } else {
                        A2D_ELEM(Mctf,  ip,  j) = fctfij;
                        A2D_ELEM(Mctf, -ip, -j) = fctfij;
                    }
                }
            }
        }

        resizeMap(Mctf, orixdim);

        Mctf.setXmippOrigin();
        // From whole to half
        for (int i = 0 ; i < YSIZE(result); i++) {
            // Don't take the middle row of the half-transform
            if (i != YSIZE(result) / 2) {
                int ip = (i < XSIZE(result)) ? i : i - YSIZE(result);
                // Don't take the last column from the half-transform
                for (int j = 0; j < XSIZE(result) - 1; j++) {
                    // Make just one lookup on Mctf.data
                    RFLOAT mctfipj = A2D_ELEM(Mctf, ip, j);
                    if (ctf_premultiplied) {
                        // Constrain result[i, j] to the interval [0.0, 1.0]
                        if (mctfipj < 0.0) {
                            DIRECT_A2D_ELEM(result, i, j) = 0.0;
                        } else if (mctfipj > 1.0) {
                            DIRECT_A2D_ELEM(result, i, j) = 1.0;
                        } else {
                            DIRECT_A2D_ELEM(result, i, j) = sqrt(mctfipj);
                        }
                    } else {
                        DIRECT_A2D_ELEM(result, i, j) = mctfipj;
                    }
                }
            }
        }
    } else {
        RFLOAT xs = (RFLOAT)orixdim * angpix;
        RFLOAT ys = (RFLOAT)oriydim * angpix;

        if (obsModel != 0 && obsModel->hasEvenZernike) {

            ensure_square(orixdim, oriydim);

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

            for (int i = 0; i < result.ydim; i++)
            for (int j = 0; j < result.xdim; j++) {
                RFLOAT x = j / xs;
                RFLOAT y = i <= result.ydim / 2 ? i / ys : (i - result.ydim) / ys;

                const int x0 = j;
                const int y0 = i <= result.ydim / 2 ? i : gammaOffset.data.ydim + i - result.ydim;

                DIRECT_A2D_ELEM(result, i, j) = getCTF(
                    x, y, do_abs, do_only_flip_phases, do_intact_until_first_peak,
                    do_damping, gammaOffset(y0, x0), do_intact_after_first_peak
                );
            }
        } else {
            FOR_ALL_ELEMENTS_IN_FFTW_TRANSFORM2D(result) {
                RFLOAT x = (RFLOAT)jp / xs;
                RFLOAT y = (RFLOAT)ip / ys;

                DIRECT_A2D_ELEM(result, i, j) = getCTF(
                    x, y, do_abs, do_only_flip_phases, do_intact_until_first_peak,
                    do_damping, 0.0, do_intact_after_first_peak
                );
            }
        }
    }
}

// Generate a complete CTFP (complex) image (with sector along angle) ---------
void CTF::getCTFPImage(
    MultidimArray<Complex> &result, int orixdim, int oriydim, RFLOAT angpix,
    bool is_positive, float angle
) {
    if (angle < 0.0 || angle >= 360.0)
        REPORT_ERROR("CTF::getCTFPImage: angle should be in the interval [0,360)");

    // Flip angles greater than 180 degrees
    if (angle >= 180.0) {
        angle -= 180.0;
        is_positive = !is_positive;
    }

    float anglerad = radians(angle);

    RFLOAT xs = (RFLOAT) orixdim * angpix;
    RFLOAT ys = (RFLOAT) oriydim * angpix;

    if (obsModel != 0 && obsModel->hasEvenZernike) {
        
        ensure_square(orixdim, oriydim);

        const Image<RFLOAT>& gammaOffset = obsModel->getGammaOffset(opticsGroup, oriydim);

        if (gammaOffset.data.xdim < result.xdim || gammaOffset.data.ydim < result.ydim) {
            REPORT_ERROR_STR(
                "CTF::getFftwImage: size requested for output image "
                << "is greater than size of original image: "
                << result.xdim << "×" << result.ydim << " was requested, "
                << "but only " << gammaOffset.data.xdim << "×" << gammaOffset.data.ydim << " is available\n"
            );
        }

        // Why do we have i <= YSIZE(result) / 2 here, but i < XSIZE(result) below?
        for (int i = 0, ip = 0; i < YSIZE(result); i++, ip = i <= YSIZE(result) / 2 ? i : i - YSIZE(result))
        for (int j = 0, jp = 0; j < XSIZE(result); j++, jp = j) {
            RFLOAT x = (RFLOAT)jp / xs;
            RFLOAT y = (RFLOAT)ip / ys;
            RFLOAT myangle = (x * x + y * y > 0) ? acos(y / Pythag(x, y)) : 0; // dot-product with Y-axis: (0,1)
            const int x0 = j;
            const int y0 = i <= YSIZE(result) / 2 ? ip : YSIZE(gammaOffset.data) + ip;

            DIRECT_A2D_ELEM(result, i, j) = getCTFP(
                x, y,
                (myangle >= anglerad ? is_positive : !is_positive),
                gammaOffset(y0, x0)
            );
        }
    } else {
        // FOR_ALL_ELEMENTS_IN_FFTW_TRANSFORM2D(result) {
        for (long int i = 0, ip = 0; i < YSIZE(result); i++, ip = i < XSIZE(result) ? i : i - YSIZE(result)) \
        for (long int j = 0, jp = 0; j < XSIZE(result); j++, jp = j) {
            // If i < XSIZE(result), ip = i. Else, ip = i - YSIZE(result).
            RFLOAT x = (RFLOAT)jp / xs;
            RFLOAT y = (RFLOAT)ip / ys;
            RFLOAT myangle = (x * x + y * y > 0) ? acos(y / Pythag(x, y)) : 0; // dot-product with Y-axis: (0,1)

            DIRECT_A2D_ELEM(result, i, j) = getCTFP(
                x, y,
                (myangle >= anglerad ? is_positive : !is_positive)
            );
        }
    }
    // Special line along the vertical (Y-)axis, where FFTW stores both Friedel mates and Friedel symmetry needs to remain
    if (angle == 0.0) {
        int dim = YSIZE(result);

        for (int i = dim / 2 + 1; i < dim; i++) {
            DIRECT_A2D_ELEM(result, i, 0) = conj(DIRECT_A2D_ELEM(result, dim - i, 0));
        }
    }
}

void CTF::getCenteredImage(
    MultidimArray<RFLOAT> &result, RFLOAT Tm,
    bool do_abs, bool do_only_flip_phases, bool do_intact_until_first_peak,
    bool do_damping, bool do_intact_after_first_peak
) {
    result.setXmippOrigin();
    RFLOAT xs = (RFLOAT)XSIZE(result) * Tm;
    RFLOAT ys = (RFLOAT)YSIZE(result) * Tm;

    FOR_ALL_ELEMENTS_IN_ARRAY2D(result) {
        RFLOAT x = (RFLOAT)j / xs;
        RFLOAT y = (RFLOAT)i / ys;
        A2D_ELEM(result, i, j) = getCTF(
            x, y, do_abs, do_only_flip_phases, do_intact_until_first_peak,
            do_damping, 0.0, do_intact_after_first_peak
        );
    }
}

void CTF::get1DProfile(
    MultidimArray<RFLOAT> &result, RFLOAT angle, RFLOAT Tm,
    bool do_abs, bool do_only_flip_phases, bool do_intact_until_first_peak,
    bool do_damping, bool do_intact_after_first_peak
) {

    result.setXmippOrigin();
    RFLOAT xs = (RFLOAT) XSIZE(result) * Tm; // assuming result is at the image size!

    FOR_ALL_ELEMENTS_IN_ARRAY1D(result) {
        RFLOAT x = (RFLOAT) i * cos(radians(angle)) / xs;
        RFLOAT y = (RFLOAT) i * sin(radians(angle)) / xs;
        A1D_ELEM(result, i) = getCTF(
            x, y, do_abs, do_only_flip_phases, do_intact_until_first_peak,
            do_damping, 0.0, do_intact_after_first_peak
        );
    }
}

void CTF::applyWeightEwaldSphereCurvature(
    MultidimArray<RFLOAT> &result, int orixdim, int oriydim,
    RFLOAT angpix, RFLOAT particle_diameter
) {
    RFLOAT xs = (RFLOAT) orixdim * angpix;
    RFLOAT ys = (RFLOAT) oriydim * angpix;

    Matrix2D<RFLOAT> M(2,2);

    if (obsModel != 0 && obsModel->hasMagMatrices) {
        M = obsModel->getMagMatrix(opticsGroup);
    } else {
        M.initIdentity();
    }

    FOR_ALL_ELEMENTS_IN_FFTW_TRANSFORM2D(result) {
        RFLOAT xu = (RFLOAT) jp / xs;
        RFLOAT yu = (RFLOAT) ip / ys;

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

        DIRECT_A2D_ELEM(result, i, j) = 0.5 * (A * (2.0 * fabs(sin(gamma)) - 1.0) + 1.0);
        // Within RELION, sin(chi) is used rather than 2 * sin(chi).
        // Hence the 0.5 above to keep everything on the same scale.
    }
}

void CTF::applyWeightEwaldSphereCurvature_new(
    MultidimArray<RFLOAT>& result, int orixdim, int oriydim,
    RFLOAT angpix, RFLOAT particle_diameter
) {
    const int s = oriydim;
    const int half_s = s / 2 + 1;
    const double as = angpix * s;
    const double Dpx = particle_diameter / angpix;

    for (int yi = 0; yi <      s;  yi++)
    for (int xi = 0; xi < half_s; xi++) {
        const double x = xi / as;
        const double y = yi < half_s ? yi / as : (yi - s) / as;

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

        DIRECT_A2D_ELEM(result, yi, xi) = 0.5 * (A * (2.0 * ctf_val - 1.0) + 1.0);
        // Within RELION, sin(chi) is used rather than 2 * sin(chi).
        // Hence the 0.5 above to keep everything on the same scale.
    }
}

void CTF::applyWeightEwaldSphereCurvature_noAniso(
    MultidimArray <RFLOAT> &result, int orixdim, int oriydim,
    RFLOAT angpix, RFLOAT particle_diameter
) {
    RFLOAT xs = (RFLOAT) orixdim * angpix;
    RFLOAT ys = (RFLOAT) oriydim * angpix;
    FOR_ALL_ELEMENTS_IN_FFTW_TRANSFORM2D(result) {
        RFLOAT x = (RFLOAT) jp / xs;
        RFLOAT y = (RFLOAT) ip / ys;
        RFLOAT deltaf = fabs(getDeltaF(x, y));
        RFLOAT inv_d = Pythag(x, y);
        RFLOAT aux = 2.0 * deltaf * lambda * inv_d / particle_diameter;
        RFLOAT A = aux > 1.0 ? 0.0 : (acos(aux) - aux * sqrt(1 - aux * aux)) * 2.0 / PI;
        DIRECT_A2D_ELEM(result, i, j) = 0.5 * (A * (2.0 * fabs(getCTF(x, y)) - 1.0) + 1.0);
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
