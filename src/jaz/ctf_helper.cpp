/***************************************************************************
 *
 * Author: "Jasenko Zivanov"
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

#include <src/jaz/ctf_helper.h>
#include <src/jaz/slice_helper.h>
#include <src/projector.h>
#include <src/jaz/img_proc/filter_helper.h>
#include <src/jaz/optimization/nelder_mead.h>
#include <src/jaz/gravis/t4Matrix.h>

using namespace gravis;

std::vector<CTF> CtfHelper::loadCtffind4(
    std::string path, int imageCount, 
    double voltage, double Cs, double Q0, double Bfac, double scale
) {
    /*
     example:
    # Output from CTFFind version 4.1.5, run on 2017-03-30 15:12:45
    # Input file: /beegfs/zivanov/tomograms/ts_05/frames/05_f32.mrc ; Number of micrographs: 1
    # Pixel size: 1.000 Angstroms ; acceleration voltage: 300.0 keV ; spherical aberration: 2.70 mm ; amplitude contrast: 0.07
    # Box size: 512 pixels ; min. res.: 30.0 Angstroms ; max. res.: 5.0 Angstroms ; min. def.: 5000.0 um; max. def. 50000.0 um
    # Columns: #1 - micrograph number; #2 - defocus 1 [Angstroms]; #3 - defocus 2; #4 - azimuth of astigmatism; #5 - additional phase shift [radians]; #6 - cross correlation; #7 - spacing (in Angstroms) up to which CTF rings were fit successfully
    1.000000 10295.926758 10012.275391 -38.856349 0.000000 0.030650 5.279412
    */
    
    std::vector<CTF> ctfs(imageCount);
    
    size_t ast = path.find_first_of('*');
    
    if (ast == std::string::npos) {
        std::ifstream file(path);
        int currImg = 0;
        
        char text[4096];

        while (file.getline(text, 4096)) {
            if (text[0] == '#') continue;
            
            std::stringstream line(text);

            ctfs[currImg] = setFromFile(line, voltage, Cs, Q0, Bfac, scale);
            currImg++;
            
            if (currImg >= imageCount) break;
        }
        
        if (currImg < imageCount) {
            REPORT_ERROR_STR("Insufficient number of CTFs found in " << path << ".\n"
                             << imageCount << " requested, " << currImg << " found.\n");
        }
    } else {
        std::string fnBase = path.substr(0, ast);
        std::string fnEnd = path.substr(ast + 1);
        
        for (int i = 0; i < imageCount; i++) {
            std::stringstream sts;
            sts << i;
            std::string fnm;
            sts >> fnm;
    
            std::string fn = fnBase + fnm + fnEnd;
            std::ifstream file(fn.c_str());
    
            if (!file.is_open()) {
                REPORT_ERROR("failed to open " + fn + '\n');
            }
            
            char text[4096];
    
            while (file.getline(text, 4096)) {
                if (text[0] == '#') continue;
    
                std::stringstream line(text);
    
                ctfs[i] = setFromFile(line, voltage, Cs, Q0, Bfac, scale);
            }
        }
    }
    
    return ctfs;
}

CTF CtfHelper::setFromFile(
    std::stringstream &line,
    double voltage, double Cs, double Q0, double Bfac, double scale
) {
    /*
    #1 - micrograph number;
    #2 - defocus 1 [Angstroms];
    #3 - defocus 2;
    #4 - azimuth of astigmatism;
    #5 - additional phase shift [radians];
    #6 - cross correlation;
    #7 - spacing (in Angstroms) up to which CTF rings were fit successfully
    */

    double imgNumber, defocus1, defocus2, azimuth, phaseShift, crossCorr, bestBefore;

    line >> imgNumber;
    line >> defocus1;
    line >> defocus2;
    line >> azimuth;
    line >> phaseShift;
    line >> crossCorr;
    line >> bestBefore;

    return CTF(defocus1, defocus2, azimuth, voltage, Cs, Q0, Bfac, scale, phaseShift);
}

void CtfHelper::setValuesByGroup(
    CTF &ctf, ObservationModel *obs, int opticsGroup,
    RFLOAT defU, RFLOAT defV, RFLOAT defAng,
    RFLOAT Bfac, RFLOAT scale, RFLOAT phase_shift
) {
    ctf.DeltafU         = defU;
    ctf.DeltafV         = defV;
    ctf.azimuthal_angle = defAng;

    ctf.Bfac            = Bfac;
    ctf.scale           = scale;
    ctf.phase_shift     = phase_shift;

    ctf.kV = obs->opticsMdt.getValue<RFLOAT>(EMDL::CTF_VOLTAGE, opticsGroup);
    ctf.Cs = obs->opticsMdt.getValue<RFLOAT>(EMDL::CTF_CS,      opticsGroup);
    ctf.Q0 = obs->opticsMdt.getValue<RFLOAT>(EMDL::CTF_Q0,      opticsGroup);

    ctf.initialise();
}

/* Read -------------------------------------------------------------------- */
void CtfHelper::readByGroup(
    CTF &ctf, const MetaDataTable &partMdt, ObservationModel *obsModel, long int particle
) {
    int opticsGroup = obsModel ? partMdt.getValue<int>(EMDL::IMAGE_OPTICS_GROUP, particle) - 1 : -1;

    ctf.kV              = readValue(EMDL::CTF_VOLTAGE,       200,         particle, opticsGroup, partMdt, obsModel);
    ctf.DeltafU         = readValue(EMDL::CTF_DEFOCUSU,      0,           particle, opticsGroup, partMdt, obsModel);
    ctf.DeltafV         = readValue(EMDL::CTF_DEFOCUSV,      ctf.DeltafU, particle, opticsGroup, partMdt, obsModel);
    ctf.azimuthal_angle = readValue(EMDL::CTF_DEFOCUS_ANGLE, 0,           particle, opticsGroup, partMdt, obsModel);
    ctf.Cs              = readValue(EMDL::CTF_CS,            0,           particle, opticsGroup, partMdt, obsModel);
    ctf.Bfac            = readValue(EMDL::CTF_BFACTOR,       0,           particle, opticsGroup, partMdt, obsModel);
    ctf.scale           = readValue(EMDL::CTF_SCALEFACTOR,   1,           particle, opticsGroup, partMdt, obsModel);
    ctf.Q0              = readValue(EMDL::CTF_Q0,            0,           particle, opticsGroup, partMdt, obsModel);
    ctf.phase_shift     = readValue(EMDL::CTF_PHASESHIFT,    0,           particle, opticsGroup, partMdt, obsModel);

    ctf.initialise();
}

CTF CtfHelper::makeCTF(const MetaDataTable &partMdt, ObservationModel *obs, long int particle) {
    CTF ctf;
    readByGroup(ctf, partMdt, obs, particle);
    return ctf;
}

CTF CtfHelper::makeCTF(const MetaDataTable &MD1, const MetaDataTable &MD2, long int objectID) {
    CTF ctf;
    read(ctf, MD1, MD2, objectID);
    return ctf;
}

CTF CtfHelper::makeCTF(
    ObservationModel *obs, int opticsGroup,
    RFLOAT defU, RFLOAT defV, RFLOAT defAng,
    RFLOAT Bfac, RFLOAT scale, RFLOAT phase_shift
) {
    CTF ctf;
    setValuesByGroup(ctf, obs, opticsGroup, defU, defV, defAng, Bfac, scale, phase_shift);
    return ctf;
}

// Read from a MetaDataTable
void CtfHelper::read(CTF &ctf, const MetaDataTable &MD) {
    MetaDataTable MDempty;
    MDempty.addObject(); // add one empty object
    read(ctf, MD, MDempty);
}

template <typename T>
T getMDT(EMDL::EMDLabel label, const MetaDataTable &mdt1, const MetaDataTable &mdt2, long int objectID, T defval) {
    try {
        return mdt1.getValue<T>(label, objectID);
    } catch (const char *errmsg) { try {
        return mdt2.getValue<T>(label, objectID);
    } catch (const char *errmsg) {
        return defval;
    } }
}

// Read parameters from MetaDataTables containing micrograph/particle information
void CtfHelper::read(CTF &ctf, const MetaDataTable &MD1, const MetaDataTable &MD2, long int objectID) {

    // Parameterse that MD1 does not contain, are tried to be read from MD2.
    ctf.kV              = getMDT<RFLOAT>(EMDL::CTF_VOLTAGE,       MD1, MD2, objectID, 200);
    ctf.DeltafU         = getMDT<RFLOAT>(EMDL::CTF_DEFOCUSU,      MD1, MD2, objectID, 0);
    ctf.DeltafV         = getMDT<RFLOAT>(EMDL::CTF_DEFOCUSV,      MD1, MD2, objectID, ctf.DeltafU);
    ctf.azimuthal_angle = getMDT<RFLOAT>(EMDL::CTF_DEFOCUS_ANGLE, MD1, MD2, objectID, 0);
    ctf.Cs              = getMDT<RFLOAT>(EMDL::CTF_CS,            MD1, MD2, objectID, 0);
    ctf.Bfac            = getMDT<RFLOAT>(EMDL::CTF_BFACTOR,       MD1, MD2, objectID, 0);
    ctf.scale           = getMDT<RFLOAT>(EMDL::CTF_SCALEFACTOR,   MD1, MD2, objectID, 1);
    ctf.Q0              = getMDT<RFLOAT>(EMDL::CTF_Q0,            MD1, MD2, objectID, 0);
    ctf.phase_shift     = getMDT<RFLOAT>(EMDL::CTF_PHASESHIFT,    MD1, MD2, objectID, 0);

    ctf.initialise();
}

RFLOAT CtfHelper::readValue(
    EMDL::EMDLabel label, RFLOAT defaultVal,
    long int particle, int opticsGroup,
    const MetaDataTable &partMdt, const ObservationModel* obs
) {
    try {
        return partMdt.getValue<RFLOAT>(label, particle);
    } catch (const char *errmsg) { try {
        if (opticsGroup < 0) { throw "Negative optics group!"; }
        if (!obs) { throw "No ObservationModel!"; }
        return obs->opticsMdt.getValue<RFLOAT>(label, opticsGroup);
    } catch (const char *errmsg) {
        return defaultVal;
    } }
}

/* Read -------------------------------------------------------------------- */

// Write to an existing object in a MetaDataTable
void CtfHelper::write(CTF &ctf, MetaDataTable &MD) {
    // For versions >= 3.1: store kV, Cs, Q0 in optics table
    // MD.setValue(EMDL::CTF_VOLTAGE, ctf.kV);
    MD.setValue(EMDL::CTF_DEFOCUSU, ctf.DeltafU);
    MD.setValue(EMDL::CTF_DEFOCUSV, ctf.DeltafV);
    MD.setValue(EMDL::CTF_DEFOCUS_ANGLE, ctf.azimuthal_angle);
    // MD.setValue(EMDL::CTF_CS, ctf.Cs);
    MD.setValue(EMDL::CTF_BFACTOR, ctf.Bfac);
    MD.setValue(EMDL::CTF_SCALEFACTOR, ctf.scale);
    MD.setValue(EMDL::CTF_PHASESHIFT, ctf.phase_shift);
    // MD.setValue(EMDL::CTF_Q0, ctf.Q0);
}

void CtfHelper::write(CTF &ctf, std::ostream &out) {
    MetaDataTable MD;
    MD.addObject();
    write(ctf, MD);
    MD.write(out);
}

static void assert_square(int orixdim, int oriydim) {
    if (orixdim != oriydim) {
        REPORT_ERROR_STR(
            "CTF::getFftwImage: currently, symmetric aberrations are supported "
            << "only for square images.\n"
        );
    }
}

// Generate a complete CTF Image ----------------------------------------------
MultidimArray<RFLOAT> CtfHelper::getFftwImage(
    const CTF &ctf,
    long int Xdim, long int Ydim, int orixdim, int oriydim,
    RFLOAT angpix, ObservationModel *obsModel, int opticsGroup,
    bool do_abs, bool do_only_flip_phases, bool do_intact_until_first_peak,
    bool do_damping, bool do_ctf_padding, bool do_intact_after_first_peak
) {
    MultidimArray<RFLOAT> result(Xdim, Ydim);
    // Boxing the particle in a small box from the whole micrograph leads to loss of delocalised information (or aliaising in the CTF)
    // Here, calculate the CTF in a 2x larger box to support finer oscillations,
    // and then rescale the large CTF to simulate the effect of the windowing operation
    if (do_ctf_padding) {
        bool ctf_premultiplied = obsModel && obsModel->getCtfPremultiplied(opticsGroup);

        // two-fold padding, increased to 4-fold padding for pre-multiplied CTFs
        int orixdim_pad = 2 * orixdim;
        int oriydim_pad = 2 * oriydim;
        // Such a big box might not be necessary
        if (ctf_premultiplied) {
            orixdim_pad *= 2;
            oriydim_pad *= 2;
        }

        MultidimArray<RFLOAT> Fctf = getFftwImage(
            ctf,
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
                RFLOAT t = ctf.getCTF(
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
                RFLOAT t = ctf.getCTF(
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

Complex _get_ctfp_(
    int ip, int jp, RFLOAT xs, RFLOAT ys,
    const CTF &ctf, ObservationModel *obsModel, int opticsGroup,
    float anglerad, bool is_positive, double gamma_offset
) {
    RFLOAT x = (RFLOAT) ip / xs;
    RFLOAT y = (RFLOAT) jp / ys;
    RFLOAT myangle = x * x + y * y > 0 ? acos(y / Pythag(x, y)) : 0; // dot-product with Y-axis: (0, 1)

    if (obsModel) obsModel->magnify(x, y, obsModel->getMagMatrix(opticsGroup));
    Complex ctfp = ctf.getCTFP(x, y, gamma_offset);
    return (myangle >= anglerad) == is_positive ? ctfp : ctfp.conj();
}

// Generate a complete CTFP (complex) image (with sector along angle) ---------
MultidimArray<Complex> CtfHelper::getCTFPImage(
    const CTF &ctf,
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
            const int x0 = i;
            const int y0 = j <= Ysize(result) / 2 ? jp : Ysize(gammaOffset.data) + jp;
            direct::elem(result, i, j) = _get_ctfp_(ip, jp, xs, ys, ctf, obsModel, opticsGroup, anglerad, is_positive, gammaOffset(y0, x0));
        }
    } else {
        // FOR_ALL_ELEMENTS_IN_FFTW_TRANSFORM2D(result) {
        for (long int j = 0, jp = 0; j < Ysize(result); j++, jp = j < Xsize(result) ? j : j - Ysize(result))
        for (long int i = 0, ip = 0; i < Xsize(result); i++, ip = i) {
            // If i < Xsize(result), ip = i. Else, ip = i - Ysize(result).
            direct::elem(result, i, j) = _get_ctfp_(ip, jp, xs, ys, ctf, obsModel, opticsGroup, anglerad, is_positive, 0.0);
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

MultidimArray<RFLOAT> CtfHelper::getCenteredImage(
    const CTF &ctf,
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
        RFLOAT t = ctf.getCTF(
            x, y,
            do_only_flip_phases, do_intact_until_first_peak,
            do_damping, 0.0, do_intact_after_first_peak
        );
        result.elem(i, j) = f(t);
    }
    return result;
}

MultidimArray<RFLOAT> CtfHelper::get1DProfile(
    const CTF &ctf, long int Xdim, long int Ydim, RFLOAT angle, RFLOAT Tm,
    ObservationModel *obsModel, int opticsGroup,
    bool do_abs, bool do_only_flip_phases, bool do_intact_until_first_peak,
    bool do_damping, bool do_intact_after_first_peak
) {
    MultidimArray<RFLOAT> result (Xdim, Ydim);
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
        RFLOAT t = ctf.getCTF(
            x, y,
            do_only_flip_phases, do_intact_until_first_peak,
            do_damping, 0.0, do_intact_after_first_peak
        );
        result.elem(i) = f(t);
    }
    return result;
}

void CtfHelper::applyWeightEwaldSphereCurvature(
    CTF &ctf,
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

        const RFLOAT astigDefocus = ctf.astigDefocus(x, y);
        RFLOAT u2 = x * x + y * y;
        RFLOAT u4 = u2 * u2;
        RFLOAT gamma = ctf.getGamma(x, y);  // XXX Computes u2, u4, astigDefocus() again internally

        RFLOAT deltaf = u2 > 0.0 ? std::abs(astigDefocus / u2) : 0.0;
        RFLOAT inv_d = sqrt(u2);
        RFLOAT aux = 2.0 * deltaf * ctf.getLambda() * inv_d / particle_diameter;
        RFLOAT A = aux > 1.0 ? 0.0 : (acos(aux) - aux * sqrt(1 - aux * aux)) * 2.0 / PI;

        direct::elem(result, i, j) = 0.5 * (A * (2.0 * fabs(sin(gamma)) - 1.0) + 1.0);
        // Within RELION, sin(chi) is used rather than 2 * sin(chi).
        // Hence the 0.5 above to keep everything on the same scale.
    }
}

void CtfHelper::applyWeightEwaldSphereCurvature_new(
    CTF &ctf,
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
        const t2Vector<RFLOAT> shift2D = RFLOAT(1.0 / (2 * angpix * PI)) * ctf.getGammaGrad(x, y);
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
        const double ctf_val = ctf.getCTF(x, y, true, false, false, false, 0.0);

        direct::elem(result, xi, yi) = 0.5 * (A * (2.0 * ctf_val - 1.0) + 1.0);
        // Within RELION, sin(chi) is used rather than 2 * sin(chi).
        // Hence the 0.5 above to keep everything on the same scale.
    }
}

void CtfHelper::applyWeightEwaldSphereCurvature_noAniso(
    CTF &ctf,
    MultidimArray <RFLOAT> &result, int orixdim, int oriydim,
    RFLOAT angpix, ObservationModel *obsModel, int opticsGroup, RFLOAT particle_diameter
) {
    RFLOAT xs = (RFLOAT) orixdim * angpix;
    RFLOAT ys = (RFLOAT) oriydim * angpix;
    FOR_ALL_ELEMENTS_IN_FFTW_TRANSFORM2D(result) {
        RFLOAT x = (RFLOAT) ip / xs;
        RFLOAT y = (RFLOAT) jp / ys;
        RFLOAT deltaf = fabs(ctf.getDeltaF(x, y));
        RFLOAT inv_d = Pythag(x, y);
        RFLOAT aux = 2.0 * deltaf * ctf.lambda * inv_d / particle_diameter;
        RFLOAT A = aux > 1.0 ? 0.0 : (acos(aux) - aux * sqrt(1 - aux * aux)) * 2.0 / PI;
        if (obsModel) obsModel->magnify(x, y, obsModel->getMagMatrix(opticsGroup));
        direct::elem(result, i, j) = 0.5 * (A * (2.0 * fabs(ctf.getCTF(x, y)) - 1.0) + 1.0);
        // Within RELION, sin(chi) is used rather than 2 * sin(chi).
        // Hence the 0.5 above to keep everything on the same scale.
    }
}
