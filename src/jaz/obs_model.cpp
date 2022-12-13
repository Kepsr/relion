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

#include "src/jaz/obs_model.h"
#include "src/jaz/stack_helper.h"
#include "src/jaz/img_proc/filter_helper.h"
#include "src/jaz/Fourier_helper.h"
#include "src/jaz/ctf/tilt_helper.h"
#include "src/jaz/math/Zernike.h"
#include "src/jaz/vtk_helper.h"
#include "src/jaz/io/star_converter.h"
#include "src/jaz/ctf_helper.h"

#include <src/backprojector.h>

#include <set>
#include <omp.h>

using namespace gravis;

void ObservationModel::loadSafely(
    std::string filename, ObservationModel& obsModel,
    MetaDataTable &particlesMdt, std::string tablename,
    int verb, bool do_die_upon_error
) {

    std::string mytablename;

    if (tablename == "discover") {
        if (particlesMdt.read(filename, "particles")) {
            mytablename = "particles";
        } else if (particlesMdt.read(filename, "micrographs")) {
            mytablename = "micrographs";
        } else if (particlesMdt.read(filename, "movies")) {
            mytablename = "movies";
        }
    } else {
        particlesMdt.read(filename, tablename);
        mytablename = tablename;
    }

    MetaDataTable opticsMdt;
    opticsMdt.read(filename, "optics");

    if (opticsMdt.empty()) {
        if (verb > 0) {
            std::cerr << "WARNING: " << filename << " seems to be from a previous version of Relion. Attempting conversion...\n";
            std::cerr << "         You should make sure metadata in the optics group table after conversion is correct.\n";
        }

        MetaDataTable oldMdt;
        oldMdt.read(filename);

        StarConverter::convert_3p0_particlesTo_3p1(oldMdt, particlesMdt, opticsMdt, mytablename, do_die_upon_error);
        if (!do_die_upon_error && opticsMdt.empty()) return;  // return an empty optics table if error was raised

        if (mytablename.empty() || mytablename == "discover") {
            particlesMdt.name =
                particlesMdt.containsLabel(EMDL::IMAGE_NAME)            ? "particles" :
                particlesMdt.containsLabel(EMDL::MICROGRAPH_MOVIE_NAME) ? "movies" :
                                                                          "micrographs";
        }
    }

    obsModel = ObservationModel(opticsMdt, do_die_upon_error);
    if (!do_die_upon_error && obsModel.opticsMdt.empty()) return;  // return an empty optics table if error was raised

    // make sure all optics groups are defined

    const std::vector<int> undefinedOptGroups = obsModel.findUndefinedOptGroups(particlesMdt);

    if (!undefinedOptGroups.empty()) {
        std::vector<std::string> v;
        v.reserve(undefinedOptGroups.size());
        for (int i : undefinedOptGroups)
            v.push_back(std::to_string(i));
        REPORT_ERROR("ERROR: The following optics groups were not defined in " + filename + ": " + join(v, ", "));
    }

    // make sure the optics groups appear in the right order (and rename them if necessary)
    if (!obsModel.opticsGroupsSorted()) {
        if (verb > 0) {
            std::cerr << "   - Warning: the optics groups in " << filename
                      << " are not in the right order - renaming them now" << std::endl;
        }

        obsModel.sortOpticsGroups(particlesMdt);
    }

    if (mytablename != "particles" && obsModel.opticsMdt.containsLabel(EMDL::IMAGE_PIXEL_SIZE)) {
        std::cerr << "WARNING: This is not a particle STAR file but contains rlnImagePixelSize column." << std::endl;
        if (!obsModel.opticsMdt.containsLabel(EMDL::MICROGRAPH_PIXEL_SIZE)) {
            std::cerr << "Pixel size in rlnImagePixelSize will be copied to rlnMicrographPixelSize column. Please make sure this is correct!" << std::endl;

            for (long int i : obsModel.opticsMdt) {
                const RFLOAT image_angpix = obsModel.opticsMdt.getValue<RFLOAT>(EMDL::IMAGE_PIXEL_SIZE, i);
                obsModel.opticsMdt.setValue(EMDL::MICROGRAPH_PIXEL_SIZE, image_angpix, i);
            }
        }
    }
}

void ObservationModel::saveNew(
    MetaDataTable &particlesMdt, MetaDataTable &opticsMdt,
    std::string filename, std::string tablename
) {
    std::string tmpfilename = filename + ".tmp";
    std::ofstream of(tmpfilename);

    opticsMdt.name = "optics";
    opticsMdt.write(of);

    particlesMdt.name = tablename;
    particlesMdt.write(of);

    std::rename(tmpfilename.c_str(), filename.c_str());
}

void ObservationModel::save(
    MetaDataTable &particlesMdt, std::string filename, std::string tablename
) {
    saveNew(particlesMdt, opticsMdt, filename, tablename);
}

ObservationModel::ObservationModel() {}

ObservationModel::ObservationModel(const MetaDataTable &_opticsMdt, bool do_die_upon_error):
opticsMdt(_opticsMdt),
angpix(_opticsMdt.size()),
lambda(_opticsMdt.size()),
Cs(_opticsMdt.size()),
boxSizes(_opticsMdt.size(), 0.0),
CtfPremultiplied(_opticsMdt.size(), false
) {

    if (
        !opticsMdt.containsLabel(EMDL::CTF_VOLTAGE) ||
        !opticsMdt.containsLabel(EMDL::CTF_CS) ||
        !opticsMdt.containsLabel(EMDL::IMAGE_PIXEL_SIZE) &&
        !opticsMdt.containsLabel(EMDL::MICROGRAPH_PIXEL_SIZE) &&
        !opticsMdt.containsLabel(EMDL::MICROGRAPH_ORIGINAL_PIXEL_SIZE)
    ) {
        if (do_die_upon_error) {
            REPORT_ERROR_STR(
                "ERROR: not all necessary variables defined in _optics.star file: "
                << "rlnPixelSize, rlnVoltage and rlnSphericalAberration. Make sure to convert older STAR files anew in version-3.1, "
                << "with relion_convert_star."
            );
        } else {
            opticsMdt.clear();
            return;
        }
    }

    // symmetrical high-order aberrations:
    hasEvenZernike = opticsMdt.containsLabel(EMDL::IMAGE_EVEN_ZERNIKE_COEFFS);
    evenZernikeCoeffs = std::vector<std::vector<double>>(opticsMdt.size(), std::vector<double>(0));
    gammaOffset = std::vector<std::map<int, Image<RFLOAT>>>(opticsMdt.size());

    // antisymmetrical high-order aberrations:
    hasOddZernike = opticsMdt.containsLabel(EMDL::IMAGE_ODD_ZERNIKE_COEFFS);
    oddZernikeCoeffs = std::vector<std::vector<double>>(opticsMdt.size(), std::vector<double>(0));
    phaseCorr = std::vector<std::map<int, Image<Complex>>>(opticsMdt.size());

    const bool hasTilt = opticsMdt.containsLabel(EMDL::IMAGE_BEAMTILT_X) ||
                         opticsMdt.containsLabel(EMDL::IMAGE_BEAMTILT_Y);

    // anisotropic magnification:
    hasMagMatrices = opticsMdt.containsLabel(EMDL::IMAGE_MAG_MATRIX_00) ||
                     opticsMdt.containsLabel(EMDL::IMAGE_MAG_MATRIX_01) ||
                     opticsMdt.containsLabel(EMDL::IMAGE_MAG_MATRIX_10) ||
                     opticsMdt.containsLabel(EMDL::IMAGE_MAG_MATRIX_11);

    magMatrices.resize(opticsMdt.size());

    hasBoxSizes = opticsMdt.containsLabel(EMDL::IMAGE_SIZE);

    if (opticsMdt.containsLabel(EMDL::IMAGE_OPTICS_GROUP_NAME)) {
        groupNames.resize(opticsMdt.size());
    }

    if (opticsMdt.containsLabel(EMDL::IMAGE_MTF_FILENAME)) {
        fnMtfs.resize(opticsMdt.size());
        mtfImage = std::vector<std::map<int, Image<RFLOAT>>>(opticsMdt.size());
    }
    if (opticsMdt.containsLabel(EMDL::MICROGRAPH_ORIGINAL_PIXEL_SIZE)) {
        originalAngpix.resize(opticsMdt.size());
    }

    for (int i = 0; i < opticsMdt.size(); i++) {
        angpix[i] =
            opticsMdt.containsLabel(EMDL::IMAGE_PIXEL_SIZE) ?
            opticsMdt.getValue<RFLOAT>(EMDL::IMAGE_PIXEL_SIZE, i) :
            opticsMdt.containsLabel(EMDL::MICROGRAPH_PIXEL_SIZE) ?
            opticsMdt.getValue<RFLOAT>(EMDL::MICROGRAPH_PIXEL_SIZE, i) :
            opticsMdt.getValue<RFLOAT>(EMDL::MICROGRAPH_ORIGINAL_PIXEL_SIZE, i);

        if (opticsMdt.containsLabel(EMDL::IMAGE_OPTICS_GROUP_NAME))
            groupNames[i] = opticsMdt.getValue<std::string>(EMDL::IMAGE_OPTICS_GROUP_NAME, i);
        if (opticsMdt.containsLabel(EMDL::IMAGE_MTF_FILENAME))
            fnMtfs[i] = opticsMdt.getValue<std::string>(EMDL::IMAGE_MTF_FILENAME, i);
        if (opticsMdt.containsLabel(EMDL::MICROGRAPH_ORIGINAL_PIXEL_SIZE))
            originalAngpix[i] = opticsMdt.getValue<RFLOAT>(EMDL::MICROGRAPH_ORIGINAL_PIXEL_SIZE, i);
        if (opticsMdt.containsLabel(EMDL::OPTIMISER_DATA_ARE_CTF_PREMULTIPLIED)) {
            CtfPremultiplied[i] = opticsMdt.getValue<bool>(EMDL::OPTIMISER_DATA_ARE_CTF_PREMULTIPLIED, i);
        }
        boxSizes[i] = opticsMdt.getValue<int>(EMDL::IMAGE_SIZE, i);

        const double kV = opticsMdt.getValue<double>(EMDL::CTF_VOLTAGE, i), V = kV * 1e3;
        lambda[i] = 12.2643247 / sqrt(V * (1.0 + V * 0.978466e-6));

        Cs[i] = opticsMdt.getValue<RFLOAT>(EMDL::CTF_CS, i);

        if (hasEvenZernike)
        evenZernikeCoeffs[i] = opticsMdt.getValue<std::vector<RFLOAT>>(EMDL::IMAGE_EVEN_ZERNIKE_COEFFS, i);

        if (hasOddZernike)
        oddZernikeCoeffs[i] = opticsMdt.getValue<std::vector<RFLOAT>>(EMDL::IMAGE_ODD_ZERNIKE_COEFFS, i);

        if (hasTilt) {
            const double tx = opticsMdt.getValue<double>(EMDL::IMAGE_BEAMTILT_X, i);
            const double ty = opticsMdt.getValue<double>(EMDL::IMAGE_BEAMTILT_Y, i);

            if (!hasOddZernike) {
                oddZernikeCoeffs[i] = std::vector<double>(6, 0.0);
            }

            TiltHelper::insertTilt(oddZernikeCoeffs[i], tx, ty, Cs[i], lambda[i]);
        }

        // Always keep a set of mag matrices
        // If none are defined, keep a set of identity matrices

        auto &magMatrix = magMatrices[i];
        magMatrix = Matrix2D<RFLOAT>::identity(2);

        // See if there is more than one MTF, for more rapid divideByMtf
        hasMultipleMtfs = std::adjacent_find(
            fnMtfs.begin(), fnMtfs.end(), std::not_equal_to<FileName>()
        ) != fnMtfs.end();

        if (hasMagMatrices) {
            magMatrix(0, 0) = opticsMdt.getValue<RFLOAT>(EMDL::IMAGE_MAG_MATRIX_00, i);
            magMatrix(0, 1) = opticsMdt.getValue<RFLOAT>(EMDL::IMAGE_MAG_MATRIX_01, i);
            magMatrix(1, 0) = opticsMdt.getValue<RFLOAT>(EMDL::IMAGE_MAG_MATRIX_10, i);
            magMatrix(1, 1) = opticsMdt.getValue<RFLOAT>(EMDL::IMAGE_MAG_MATRIX_11, i);
        }
    }

    hasOddZernike |= hasTilt;

}

MultidimArray<Complex> ObservationModel::predictObservation(
    const Projector &proj, const MetaDataTable &partMdt,
    long int particle, double angpix_ref,
    bool applyCtf, bool shiftPhases, bool applyShift, bool applyMtf, bool applyCtfPadding
) {

    const int s_ref = proj.ori_size;

    int opticsGroup = partMdt.getValue<int>(EMDL::IMAGE_OPTICS_GROUP, particle) - 1;

    if (!hasBoxSizes) {
        REPORT_ERROR_STR("ObservationModel::predictObservation: Unable to make a prediction without knowing the box size.\n");
    }

    const int s_out = boxSizes[opticsGroup];
    const int sh_out = s_out / 2 + 1;

    double xoff = partMdt.getValue<double>(EMDL::ORIENT_ORIGIN_X_ANGSTROM, particle) / angpix[opticsGroup];
    double yoff = partMdt.getValue<double>(EMDL::ORIENT_ORIGIN_Y_ANGSTROM, particle) / angpix[opticsGroup];

    double rot  = partMdt.getValue<double>(EMDL::ORIENT_ROT,  particle);
    double tilt = partMdt.getValue<double>(EMDL::ORIENT_TILT, particle);
    double psi  = partMdt.getValue<double>(EMDL::ORIENT_PSI,  particle);

    Matrix2D<RFLOAT> A3D = Euler::angles2matrix(rot, tilt, psi);
    if (hasMagMatrices) { A3D *= anisoMag(opticsGroup); }
    A3D *= scaleDifference(opticsGroup, s_ref, angpix_ref);

    auto pred = proj.get2DFourierTransform(sh_out, s_out, 1, A3D);

    if (applyShift) {
        shiftImageInFourierTransform(pred, s_out, s_out / 2 - xoff, s_out / 2 - yoff);
    }

    if (applyCtf) {
        const CTF ctf = CtfHelper::makeCTF(partMdt, this, particle);

        const auto ctfImg = CtfHelper::getFftwImage(
            ctf,
            sh_out, s_out, s_out, s_out, angpix[opticsGroup], this,
            false, false, false, true, applyCtfPadding
        );

        if (getCtfPremultiplied(opticsGroup))
            for (long int i = 0; i < pred.size(); i++)
                pred[i] *= ctfImg[i] * ctfImg[i];
        else
            for (long int i = 0; i < pred.size(); i++)
                pred[i] *= ctfImg[i];

    }

    if (
        shiftPhases &&
        opticsGroup < oddZernikeCoeffs.size() &&
        oddZernikeCoeffs[opticsGroup].size() > 0
    ) {
        pred *= getPhaseCorrection(opticsGroup, s_out).data;
    }

    if (applyMtf && opticsGroup < fnMtfs.size()) {
        pred *= getMtfImage(opticsGroup, s_out).data;
    }

    return pred;
}

Volume<t2Vector<Complex>> ObservationModel::predictComplexGradient(
    Projector &proj, const MetaDataTable &partMdt,
    long particle, double angpix_ref,
    bool applyCtf, bool shiftPhases, bool applyShift,
    bool applyMtf, bool applyCtfPadding
) {
    if (applyCtf || applyShift || applyCtfPadding) {
        REPORT_ERROR_STR(
            "ObservationModel::predictComplexGradient: "
            << "applyCtf and applyShift and applyCtfPadding are currently not supported\n"
        );
    }

    const int s_ref = proj.ori_size;

    int opticsGroup = partMdt.getValue<int>(EMDL::IMAGE_OPTICS_GROUP, particle) - 1;

    double xoff = partMdt.getValue<double>(EMDL::ORIENT_ORIGIN_X_ANGSTROM, particle) / angpix[opticsGroup];
    double yoff = partMdt.getValue<double>(EMDL::ORIENT_ORIGIN_Y_ANGSTROM, particle) / angpix[opticsGroup];

    double rot  = partMdt.getValue<double>(EMDL::ORIENT_ROT,  particle);
    double tilt = partMdt.getValue<double>(EMDL::ORIENT_TILT, particle);
    double psi  = partMdt.getValue<double>(EMDL::ORIENT_PSI,  particle);

    auto A3D = Euler::angles2matrix(rot, tilt, psi);
    if (hasMagMatrices) { A3D *= anisoMag(opticsGroup); }
    A3D *= scaleDifference(opticsGroup, s_ref, angpix_ref);

    const int s_out = boxSizes[opticsGroup];
    const int sh_out = s_out / 2 + 1;
    auto out = proj.projectGradient(sh_out, s_out, A3D);

    if (
        shiftPhases && opticsGroup < oddZernikeCoeffs.size() &&
        oddZernikeCoeffs[opticsGroup].size() > 0
    ) {
        const Image<Complex> &corr = getPhaseCorrection(opticsGroup, s_out);

        for (int y = 0; y < s_out;  y++)
        for (int x = 0; x < sh_out; x++) {
            out(x, y, 0).x *= corr(y, x);
            out(x, y, 0).y *= corr(y, x);
        }
    }

    if (applyMtf && fnMtfs.size() > opticsGroup) {
        const Image<RFLOAT> &mtf = getMtfImage(opticsGroup, s_out);

        for (int y = 0; y < s_out;  y++)
        for (int x = 0; x < sh_out; x++) {
            out(x, y, 0).x *= mtf(y, x);
            out(x, y, 0).y *= mtf(y, x);
        }
    }

    return out;
}

template <typename Operator>
void ObservationModel::operateByMtf(
    int opticsGroup, MultidimArray<Complex> &obsImage,
    bool do_correct_average_mtf, Operator assign  // e.g. Complex::operator *=
) {
    const int sh = obsImage.xdim, s = obsImage.ydim;
    // If there is only a single MTF and we are correcting for the average, then do nothing...
    if (do_correct_average_mtf && !hasMultipleMtfs) return;

    if (fnMtfs.size() > opticsGroup) {
        const auto &mtf    = getMtfImage(opticsGroup, s);
        const auto &avgmtf = getAverageMtfImage(s);

        if (do_correct_average_mtf) {
            for (int y = 0; y < s;  y++)
            for (int x = 0; x < sh; x++) {
                assign(obsImage.elem(y, x), mtf(y, x) / avgmtf(y, x));
            }
        } else {
            for (int y = 0; y < s;  y++)
            for (int x = 0; x < sh; x++) {
                assign(obsImage.elem(y, x), mtf(y, x));
            }
        }
    }
}

void ObservationModel::multiplyByMtf(
    int opticsGroup, MultidimArray<Complex> &obsImage, bool do_correct_average_mtf
) {
    operateByMtf(opticsGroup, obsImage, do_correct_average_mtf,
        [] (Complex &x, Complex y) { x *= y; });
}

void ObservationModel::divideByMtf(
    int opticsGroup, MultidimArray<Complex> &obsImage, bool do_correct_average_mtf
) {
    operateByMtf(opticsGroup, obsImage, do_correct_average_mtf,
        [] (Complex &x, Complex y) { x /= y; });
}

template <typename F>
void ObservationModel::operatePhase(
    int opticsGroup, MultidimArray<Complex> &obsImage, F f  // e.g. Complex::conj
) {
    if (
        opticsGroup >= oddZernikeCoeffs.size() ||
        oddZernikeCoeffs[opticsGroup].empty()
    ) return;

    const Image<Complex> &corr = getPhaseCorrection(opticsGroup, obsImage.ydim);

    // Precondition: obsImage and corr are 2D
    for (long int i = 0; i < obsImage.size(); i++) {
        obsImage[i] *= f(corr.data[i]);
    }
}

void ObservationModel::modulatePhase(
    int opticsGroup, MultidimArray<Complex> &obsImage
) {
    operatePhase(opticsGroup, obsImage, [] (Complex x) { return x; });
};

void ObservationModel::demodulatePhase(
    int opticsGroup, MultidimArray<Complex> &obsImage
) {
    operatePhase(opticsGroup, obsImage, [] (Complex x) { return x.conj(); });
}

double ObservationModel::angToPix(double a, int s, int opticsGroup) const {
    return s * angpix[opticsGroup] / a;
}

double ObservationModel::pixToAng(double p, int s, int opticsGroup) const {
    return s * angpix[opticsGroup] / p;
}

void ObservationModel::setPixelSize(int opticsGroup, RFLOAT newPixelSize) {
    if (opticsGroup < 0 || opticsGroup >= boxSizes.size()) {
        REPORT_ERROR("ObservationModel::setPixelSize: wrong opticsGroup");
    }

    angpix[opticsGroup] = newPixelSize;

    phaseCorr[opticsGroup].clear();
    gammaOffset[opticsGroup].clear();

    // mtfImage can be empty
    if (!mtfImage.empty())
        mtfImage[opticsGroup].clear();
}

double ObservationModel::getPixelSize(int opticsGroup) const {
    return angpix[opticsGroup];
}

std::vector<double> ObservationModel::getPixelSizes() const {
    return angpix;
}

double ObservationModel::getWavelength(int opticsGroup) const {
    return lambda[opticsGroup];
}

std::vector<double> ObservationModel::getWavelengths() const {
    return lambda;
}

double ObservationModel::getSphericalAberration(int opticsGroup) const {
    return Cs[opticsGroup];
}

std::vector<double> ObservationModel::getSphericalAberrations() const {
    return Cs;
}

void ObservationModel::setBoxSize(int opticsGroup, int newBoxSize) {
    if (opticsGroup < 0 || opticsGroup >= boxSizes.size()) {
        REPORT_ERROR("ObservationModel::setBoxSize: wrong opticsGroup");
    }

    boxSizes[opticsGroup] = newBoxSize;

    phaseCorr[opticsGroup].clear();
    gammaOffset[opticsGroup].clear();

    // mtfImage can be empty
    if (mtfImage.size() > 0)
        mtfImage[opticsGroup].clear();
}

int ObservationModel::getBoxSize(int opticsGroup) const {
    if (!hasBoxSizes) {
        REPORT_ERROR("ObservationModel::getBoxSize: box sizes not available. Make sure particle images are available before converting/importing STAR files from earlier versions of RELION.\n");
    }

    return boxSizes[opticsGroup];
}

void ObservationModel::getBoxSizes(std::vector<int>& sDest, std::vector<int>& shDest) const {
    if (!hasBoxSizes) {
        REPORT_ERROR("ObservationModel::getBoxSizes: box sizes not available. Make sure particle images are available before converting/importing STAR files from earlier versions of RELION.\n");
    }

     sDest.resize(boxSizes.size());
    shDest.resize(boxSizes.size());

    for (int i = 0; i < boxSizes.size(); i++) {
         sDest[i] = boxSizes[i];
        shDest[i] = boxSizes[i] / 2 + 1;
    }
}

Matrix2D<RFLOAT> ObservationModel::getMagMatrix(int opticsGroup) const {
    return magMatrices[opticsGroup];
}

void ObservationModel::setMagMatrix(int opticsGroup, const Matrix2D<RFLOAT> &M) {
    magMatrices[opticsGroup] = M;
}

std::vector<Matrix2D<RFLOAT>> ObservationModel::getMagMatrices() const {
    return magMatrices;
}

Matrix2D<RFLOAT> ObservationModel::anisoMag(int opticsGroup) const {
    const Matrix2D<RFLOAT> &magmatrix = magMatrices[opticsGroup];
    Matrix2D<RFLOAT> mag = Matrix2D<RFLOAT>::identity(3);
    mag(0, 0) = magmatrix(0, 0);
    mag(0, 1) = magmatrix(0, 1);
    mag(1, 0) = magmatrix(1, 0);
    mag(1, 1) = magmatrix(1, 1);
    return mag.inv();
}

int ObservationModel::getOpticsGroup(const MetaDataTable &particlesMdt, long int particle) const {
    try {
        return particlesMdt.getValue<int>(EMDL::IMAGE_OPTICS_GROUP, particle) - 1;
    } catch (const char *errmsg) {
        REPORT_ERROR("ObservationModel::getOpticsGroup: Failed to get optics group for particle #" + particle);
    }
}

bool ObservationModel::getCtfPremultiplied(int og) const {
    return og < CtfPremultiplied.size() && CtfPremultiplied[og];
}

void ObservationModel::setCtfPremultiplied(int og, bool val) {
    CtfPremultiplied[og] = val;
}

std::string ObservationModel::getGroupName(int og) {
    if (og < groupNames.size()) return groupNames[og];
    return std::to_string(og + 1);
}

bool ObservationModel::allPixelAndBoxSizesIdentical(const MetaDataTable &mdt) {
    int og0 = getOpticsGroup(mdt, 0);

    int boxSize0 = getBoxSize(og0);
    double angpix0 = getPixelSize(og0);

    const int pc = mdt.size();

    for (int p = 1; p < pc; p++) {
        int og = getOpticsGroup(mdt, p);

        if (og != og0) {
            int boxSize = getBoxSize(og);
            double angpix = getPixelSize(og);

            if (boxSize != boxSize0 || angpix != angpix0) {
                return false;
            }
        }
    }

    return true;
}

bool ObservationModel::containsGroup(const MetaDataTable &mdt, int group) {

    const int pc = mdt.size();
    for (int p = 0; p < pc; p++) {
        if (getOpticsGroup(mdt, p) == group) {
            return true;
        }
    }
    return false;
}

int ObservationModel::numberOfOpticsGroups() const {
    return opticsMdt.size();
}

bool ObservationModel::opticsGroupsSorted() const {
    for (int i = 0; i < opticsMdt.size(); i++) {
        if (opticsMdt.getValue<int>(EMDL::IMAGE_OPTICS_GROUP, i) != i + 1) {
            return false;
        }
    }
    return true;
}

std::vector<int> ObservationModel::findUndefinedOptGroups(const MetaDataTable &partMdt) const {

    std::set<int> definedGroups;
    for (int i = 0; i < opticsMdt.size(); i++) {
        definedGroups.insert(opticsMdt.getValue<int>(EMDL::IMAGE_OPTICS_GROUP, i));
    }

    std::vector<int> out;
    out.reserve(opticsMdt.size());

    for (long int i = 0; i < partMdt.size(); i++) {
        int og = partMdt.getValue<int>(EMDL::IMAGE_OPTICS_GROUP, i);
        if (definedGroups.find(og) == definedGroups.end()) {
            out.push_back(og);
        }
    }

    return out;
}

void ObservationModel::sortOpticsGroups(MetaDataTable &partMdt) {

    std::map<int, int> old2new;
    for (int i = 0; i < opticsMdt.size(); i++) {
        int og = opticsMdt.getValue<int>(EMDL::IMAGE_OPTICS_GROUP, i);
        old2new[og] = i + 1;
        opticsMdt.setValue(EMDL::IMAGE_OPTICS_GROUP, i + 1, i);
    }

    for (long int i = 0; i < partMdt.size(); i++) {
        int og = partMdt.getValue<int>(EMDL::IMAGE_OPTICS_GROUP, i);
        partMdt.setValue(EMDL::IMAGE_OPTICS_GROUP, old2new[og], i);
    }
}

std::vector<int> ObservationModel::getOptGroupsPresent(const MetaDataTable &partMdt) const {
    const int gc = opticsMdt.size();
    const long long int pc = partMdt.size();

    std::vector<bool> optGroupIsPresent(gc, false);
    for (long int p = 0; p < pc; p++) {
        int og = partMdt.getValue<int>(EMDL::IMAGE_OPTICS_GROUP, p);
        optGroupIsPresent[og - 1] = true;
    }

    std::vector<int> out;
    out.reserve(gc);
    for (int g = 0; g < gc; g++) {
        if (optGroupIsPresent[g]) {
            out.push_back(g);
        }
    }
    return out;
}

std::vector<std::pair<int, std::vector<int>>> ObservationModel::splitParticlesByOpticsGroup(const MetaDataTable &partMdt) const {

    std::vector<int> presentGroups = ObservationModel::getOptGroupsPresent(partMdt);
    const int pogc = presentGroups.size();
    const int ogc = opticsMdt.size();

    std::vector<int> groupToPresentGroup(ogc, -1);

    for (int pog = 0; pog < pogc; pog++) {
        const int og = presentGroups[pog];
        groupToPresentGroup[og] = pog;
    }

    std::vector<std::pair<int, std::vector<int>>> out(pogc);

    for (int pog = 0; pog < pogc; pog++) {
        out[pog] = std::make_pair(presentGroups[pog], std::vector<int>(0));
    }

    const long long int pc = partMdt.size();

    for (long int p = 0; p < pc; p++) {
        int og = partMdt.getValue<int>(EMDL::IMAGE_OPTICS_GROUP, p) - 1;
        int pog = groupToPresentGroup[og];
        out[pog].second.push_back(p);
    }

    return out;
}

const Image<RFLOAT>& ObservationModel::getMtfImage(int optGroup, int s) {
    #pragma omp critical(ObservationModel_getMtfImage)
    {
        if (mtfImage[optGroup].find(s) == mtfImage[optGroup].end()) {
            if (mtfImage[optGroup].size() > 100) {
                std::cerr << "Warning: " << (mtfImage[optGroup].size() + 1)
                          << " mtf images in cache for the same ObservationModel." << std::endl;
            }

            if (optGroup >= originalAngpix.size())
                REPORT_ERROR("For MTF correction, the rlnMicrographOriginalPixelSize column is necessary in the optics table.");

            MetaDataTable MDmtf;
            MDmtf.read(fnMtfs[optGroup]);
            MultidimArray<RFLOAT> mtf_resol, mtf_value;
            mtf_resol.resize(MDmtf.size());
            mtf_value.resize(mtf_resol);

            RFLOAT resol_inv_pixel = MDmtf.getValue<RFLOAT>(EMDL::RESOLUTION_INVPIXEL, 0);
            for (long int i : MDmtf) {
                direct::elem(mtf_resol, i) = resol_inv_pixel / originalAngpix[optGroup];  // resolution needs to be given in 1/Ang
                direct::elem(mtf_value, i) = MDmtf.getValue<RFLOAT>(EMDL::POSTPROCESS_MTF_VALUE, i);
                if (direct::elem(mtf_value, i) < 1e-10) {
                    std::cerr << " i= " << i <<  " mtf_value[i]= " << direct::elem(mtf_value, i) << std::endl;
                    REPORT_ERROR("ERROR: zero or negative values encountered in MTF curve: " + fnMtfs[optGroup]);
                }
            }

            // Calculate slope of resolution (in 1/A) per element in the MTF array, in order to interpolate below
            int i = MDmtf.size();
            RFLOAT res_per_elem = (direct::elem(mtf_resol, i - 1) - direct::elem(mtf_resol, 0)) / (RFLOAT) i;
            if (res_per_elem < 1e-10) REPORT_ERROR(" ERROR: the resolution in the MTF star file does not go up....");

            const int sh = s / 2 + 1;
            mtfImage[optGroup][s] = Image<RFLOAT>(s, sh);
            Image<RFLOAT> &img = mtfImage[optGroup][s];
            const double as = angpix[optGroup] * boxSizes[optGroup];

            for (int y = 0; y < s;  y++)
            for (int x = 0; x < sh; x++) {
                const double xx = x / as;  // logical X-coordinate in 1/A
                const double yy = y < sh ? y / as : (y - s) / as;  // logical Y-coordinate in 1/A

                RFLOAT res = sqrt(xx * xx + yy * yy);  // get resolution in 1/Ang
                int i_0 = floor(res / res_per_elem);
                RFLOAT mtf;
                if (i_0 <= 0) {
                    // Check array boundaries
                    mtf = direct::elem(mtf_value, 0);
                } else if (i_0 >= mtf_value.size() - 1) {
                    mtf = direct::elem(mtf_value,  mtf_value.size() - 1);
                } else {
                    // Linear interpolation
                    RFLOAT x_0 = direct::elem(mtf_resol, i_0);
                    RFLOAT y_0 = direct::elem(mtf_value, i_0);
                    RFLOAT x_1 = direct::elem(mtf_resol, i_0 + 1);
                    RFLOAT y_1 = direct::elem(mtf_value, i_0 + 1);
                    mtf = y_0 + (y_1 - y_0) * (res - x_0) / (x_1 - x_0);
                }
                img.data.elem(y, x) = mtf;
            }
        }
    }

    return mtfImage[optGroup][s];
}

const Image<RFLOAT>& ObservationModel::getAverageMtfImage(int s) {
    #pragma omp critical(ObservationModel_getAverageMtfImage)
    {
        if (avgMtfImage.find(s) == avgMtfImage.end()) {
            // get first mtfImage
            avgMtfImage[s] = getMtfImage(0, s);
            // Then add rest of optics groups
            for (int i = 1; i < mtfImage.size(); i++) {
                avgMtfImage[s].data += getMtfImage(i, s).data;
            }
            avgMtfImage[s].data /= (RFLOAT) mtfImage.size();
        }
    }

    return avgMtfImage[s];
}

const Image<Complex>& ObservationModel::getPhaseCorrection(int optGroup, int s) {
    #pragma omp critical(ObservationModel_getPhaseCorrection)
    {
        if (phaseCorr[optGroup].find(s) == phaseCorr[optGroup].end()) {
            if (phaseCorr[optGroup].size() > 100) {
                std::cerr << "Warning: " << (phaseCorr[optGroup].size() + 1)
                          << " phase shift images in cache for the same ObservationModel." << std::endl;
            }

            const int sh = s / 2 + 1;
            phaseCorr[optGroup][s] = Image<Complex>(s, sh);
            Image<Complex> &img = phaseCorr[optGroup][s];
            const double as = angpix[optGroup] * boxSizes[optGroup];
            const Matrix2D<RFLOAT> &M = magMatrices[optGroup];

            for (int y = 0; y < s;  y++)
            for (int x = 0; x < sh; x++) {

                double phase = 0.0;
                for (int i = 0; i < oddZernikeCoeffs[optGroup].size(); i++) {

                    const double xx0 = x / as;
                    const double yy0 = y < sh - 1 ? y / as : (y - s) / as;

                    const double xx = M(0, 0) * xx0 + M(0, 1) * yy0;
                    const double yy = M(1, 0) * xx0 + M(1, 1) * yy0;

                    Zernike::Z zmn = Zernike::Z::fromOddIndex(i);
                    phase += oddZernikeCoeffs[optGroup][i] * zmn.cart(xx, yy);

                }

                img.data.elem(y, x) = Complex::unit(phase);
            }
        }
    }

    return phaseCorr[optGroup][s];
}

const Image<RFLOAT>& ObservationModel::getGammaOffset(int optGroup, int s) {
    #pragma omp critical(ObservationModel_getGammaOffset)
    {
        if (gammaOffset[optGroup].find(s) == gammaOffset[optGroup].end()) {

            if (gammaOffset[optGroup].size() > 100) {
                std::cerr << "Warning: " << (gammaOffset[optGroup].size() + 1)
                          << " gamma offset images in cache for the same ObservationModel." << std::endl;
            }

            const int sh = s / 2 + 1;
            gammaOffset[optGroup][s] = Image<RFLOAT>(s, sh);
            Image<RFLOAT> &img = gammaOffset[optGroup][s];

            const double as = angpix[optGroup] * boxSizes[optGroup];
            const Matrix2D<RFLOAT> &M = magMatrices[optGroup];

            for (int y = 0; y < s;  y++)
            for (int x = 0; x < sh; x++) {

                double phase = 0.0;
                for (int i = 0; i < evenZernikeCoeffs[optGroup].size(); i++) {

                    const double xx0 = x / as;
                    const double yy0 = y < sh - 1 ? y / as : (y - s) / as;

                    const double xx = M(0, 0) * xx0 + M(0, 1) * yy0;
                    const double yy = M(1, 0) * xx0 + M(1, 1) * yy0;

                    Zernike::Z zmn = Zernike::Z::fromEvenIndex(i);
                    phase += evenZernikeCoeffs[optGroup][i] * zmn.cart(xx, yy);

                }

                img.data.elem(y, x) = phase;
            }
        }
    }

    return gammaOffset[optGroup][s];
}

double ObservationModel::scaleDifference(int opticsGroup, int s3D, double angpix3D) {
    return (boxSizes[opticsGroup] * angpix[opticsGroup]) / (s3D * angpix3D);
}

bool ObservationModel::containsAllColumnsNeededForPrediction(const MetaDataTable &partMdt) {
    return partMdt.containsLabel(EMDL::ORIENT_ORIGIN_X_ANGSTROM) &&
           partMdt.containsLabel(EMDL::ORIENT_ORIGIN_Y_ANGSTROM) &&
           partMdt.containsLabel(EMDL::ORIENT_ROT) &&
           partMdt.containsLabel(EMDL::ORIENT_TILT) &&
           partMdt.containsLabel(EMDL::ORIENT_PSI) &&
           partMdt.containsLabel(EMDL::PARTICLE_RANDOM_SUBSET);
}
