/***************************************************************************
 *
 * Author: "Sjors H.W. Scheres", "Takanori Nakane"
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

#include "src/micrograph_model.h"
#include "src/metadata_table.h"
#include "src/image.h"
#include "src/motioncorr_runner.h"
#include "src/renderEER.h"

// TODO: Think about first frame for local model

const RFLOAT Micrograph::NOT_OBSERVED = -9999;
const int ThirdOrderPolynomialModel::NUM_COEFFS_PER_DIM = 18;

int ThirdOrderPolynomialModel::getShiftAt(RFLOAT z, RFLOAT x, RFLOAT y, RFLOAT &shiftx, RFLOAT &shifty) const {
    const RFLOAT x2 = x * x, y2 = y * y, xy = x * y, z2 = z * z;
    const RFLOAT z3 = z * z2;

    shiftx = (coeffX[0]  * z + coeffX[1]  * z2 + coeffX[2]  * z3)
           + (coeffX[3]  * z + coeffX[4]  * z2 + coeffX[5]  * z3) * x
           + (coeffX[6]  * z + coeffX[7]  * z2 + coeffX[8]  * z3) * x2
           + (coeffX[9]  * z + coeffX[10] * z2 + coeffX[11] * z3) * y
           + (coeffX[12] * z + coeffX[13] * z2 + coeffX[14] * z3) * y2
           + (coeffX[15] * z + coeffX[16] * z2 + coeffX[17] * z3) * xy;

    shifty = (coeffY[0]  * z + coeffY[1]  * z2 + coeffY[2]  * z3)
           + (coeffY[3]  * z + coeffY[4]  * z2 + coeffY[5]  * z3) * x
           + (coeffY[6]  * z + coeffY[7]  * z2 + coeffY[8]  * z3) * x2
           + (coeffY[9]  * z + coeffY[10] * z2 + coeffY[11] * z3) * y
           + (coeffY[12] * z + coeffY[13] * z2 + coeffY[14] * z3) * y2
           + (coeffY[15] * z + coeffY[16] * z2 + coeffY[17] * z3) * xy;

    return 0;
}

MotionModel* ThirdOrderPolynomialModel::clone() const {
    return new ThirdOrderPolynomialModel(*this);
}

void ThirdOrderPolynomialModel::write(std::ostream &fh, std::string block_name) {
    MetaDataTable MD;
    MD.name = block_name;

    int coeff_idx = 0;

    // Write coeffX
    for (int i = 0; i < NUM_COEFFS_PER_DIM; i++, coeff_idx++) {
        MD.addObject();
        MD.setValue(EMDL::MICROGRAPH_MOTION_COEFFS_IDX, coeff_idx, coeff_idx);
        MD.setValue(EMDL::MICROGRAPH_MOTION_COEFF,      coeffX[i], coeff_idx);
    }

    // Write coeffY
    for (int i = 0; i < NUM_COEFFS_PER_DIM; i++, coeff_idx++) {
        MD.addObject();
        MD.setValue(EMDL::MICROGRAPH_MOTION_COEFFS_IDX, coeff_idx, coeff_idx);
        MD.setValue(EMDL::MICROGRAPH_MOTION_COEFF,      coeffY[i], coeff_idx);
    }

    MD.write(fh);
}

void ThirdOrderPolynomialModel::read(std::ifstream &fh, std::string block_name) {
    MetaDataTable MD;
    MD.readStar(fh, block_name);

    const int NUM_COEFFS = NUM_COEFFS_PER_DIM * 2;

    coeffX.resize(NUM_COEFFS_PER_DIM);
    std::fill(coeffX.begin(), coeffX.end(), 0);
    coeffY.resize(NUM_COEFFS_PER_DIM);
    std::fill(coeffY.begin(), coeffY.end(), 0);

    for (long int i : MD) {

        int idx;
        RFLOAT val;
        try {
            idx = MD.getValue<int>(EMDL::MICROGRAPH_MOTION_COEFFS_IDX, i);
            val = MD.getValue<RFLOAT>(EMDL::MICROGRAPH_MOTION_COEFF, i);
        } catch (const char *errmsg) {
            REPORT_ERROR("ThirdOrderPolynomialModel coefficients table: missing index or coefficients");
        }

        if (idx >= 0 && idx < NUM_COEFFS_PER_DIM) {
            coeffX[idx] = val;
        } else if (idx >= NUM_COEFFS_PER_DIM && idx < NUM_COEFFS) {
            coeffY[idx - NUM_COEFFS_PER_DIM] = val;
        } else {
            REPORT_ERROR("ThirdOrderPolynomialModel coefficients table: wrong index");
        }
    }

    if (MD.size() != NUM_COEFFS) {
        REPORT_ERROR("ThirdOrderPolynomialModel coefficients table: incomplete values");
    }
}

Micrograph::Micrograph(): ready(false), model(nullptr) { clearFields(); }

Micrograph::Micrograph(const Micrograph& m): ready(m.ready), model(m.model ? m.model->clone() : nullptr) {
    copyFieldsFrom(m);
}

Micrograph::Micrograph(FileName filename, FileName fnGain, RFLOAT binning, int eer_upsampling, int eer_grouping): ready(false), model(nullptr) {

    clearFields();

    if (filename.getExtension() == "star" && fnGain == "") {
        if (eer_upsampling > 0 || eer_grouping > 0)
            REPORT_ERROR("Micrograph::Micrograph: When reading STAR file, you shouldn't specify eer_upsampling and eer_grouping. They are read from the STAR file.");

        read(filename);
    } else {
        this->eer_upsampling = eer_upsampling;
        this->eer_grouping = eer_grouping;
        setMovie(filename, fnGain, binning);
    }

    ready = true;
}

Micrograph::~Micrograph() {
    delete model;
}

Micrograph& Micrograph::operator = (const Micrograph& m) {
    ready = m.ready;

    delete model;
    model = m.model ? m.model->clone() : nullptr;

    copyFieldsFrom(m);

    return *this;
}

void Micrograph::write(FileName filename) {
    checkReadyFlag("write");


    std::ofstream fh (filename.c_str());
    if (!fh) {
        REPORT_ERROR((std::string)"Micrograph::write: Cannot write file: " + filename);
    }

    MetaDataTable MD;
    MD.name = "general";
    MD.isList = true;
    const long int i = MD.addObject();
    MD.setValue(EMDL::IMAGE_SIZE_X, width, i);
    MD.setValue(EMDL::IMAGE_SIZE_Y, height, i);
    MD.setValue(EMDL::IMAGE_SIZE_Z, n_frames, i);
    MD.setValue(EMDL::MICROGRAPH_MOVIE_NAME, fnMovie, i);

    if (!fnGain.empty())
        MD.setValue(EMDL::MICROGRAPH_GAIN_NAME, fnGain, i);

    if (!fnDefect.empty())
        MD.setValue(EMDL::MICROGRAPH_DEFECT_FILE, fnDefect, i);

    MD.setValue(EMDL::MICROGRAPH_BINNING, binning, i);

    if (angpix != -1)
        MD.setValue(EMDL::MICROGRAPH_ORIGINAL_PIXEL_SIZE, angpix, i);

    if (dose_per_frame != -1)
        MD.setValue(EMDL::MICROGRAPH_DOSE_RATE, dose_per_frame, i);

    if (pre_exposure != -1)
        MD.setValue(EMDL::MICROGRAPH_PRE_EXPOSURE, pre_exposure, i);

    if (voltage != -1)
        MD.setValue(EMDL::CTF_VOLTAGE, voltage, i);

    MD.setValue(EMDL::MICROGRAPH_START_FRAME, first_frame, i); // 1-indexed

    if (EERRenderer::isEER(fnMovie)) {
        if (eer_upsampling > 0)
            MD.setValue(EMDL::MICROGRAPH_EER_UPSAMPLING, this->eer_upsampling, i);

        if (eer_grouping > 0)
            MD.setValue(EMDL::MICROGRAPH_EER_GROUPING, this->eer_grouping, i);
    }

    MD.setValue(EMDL::MICROGRAPH_MOTION_MODEL_VERSION, model ? model->getModelVersion() : (int) MOTION_MODEL_NULL, i);

    MD.write(fh);

    MD.clear();
    MD.name = "global_shift";
    for (int frame = 0; frame < n_frames; frame++) {
        MD.addObject();
        MD.setValue(EMDL::MICROGRAPH_FRAME_NUMBER, frame + 1, frame); // make 1-indexed
        MD.setValue(EMDL::MICROGRAPH_SHIFT_X, globalShiftX[frame], frame);
        MD.setValue(EMDL::MICROGRAPH_SHIFT_Y, globalShiftY[frame], frame);
    }
    MD.write(fh);

    if (model) {
        std::string block_name = "local_motion_model";
        model->write(fh, block_name);
    }

    MD.clear();
    MD.name = "hot_pixels";

    if (hotpixelX.size() != hotpixelY.size())
        REPORT_ERROR("Logic error: hotpixelX.size() != hotpixelY.size()");

    for (int i = 0, ilim = hotpixelX.size(); i < ilim; i++) {
        const long int j = MD.addObject();
        MD.setValue(EMDL::IMAGE_COORD_X, (RFLOAT) hotpixelX[i], j);
        MD.setValue(EMDL::IMAGE_COORD_Y, (RFLOAT) hotpixelY[i], j);
    }
    MD.write(fh);

    MD.clear();
    MD.name = "local_shift";

    int n_local_trajectory = localShiftX.size();
    if (
        n_local_trajectory != localShiftY.size() ||
        n_local_trajectory != patchX.size() ||
        n_local_trajectory != patchY.size() ||
        n_local_trajectory != patchZ.size()
    ) {
        REPORT_ERROR("Logic error: inconsistent local trajectory");
    }
    for (int i = 0; i < n_local_trajectory; i++) {
        const long int j = MD.addObject();
        MD.setValue(EMDL::MICROGRAPH_FRAME_NUMBER, (int) patchZ[i], j);
        MD.setValue(EMDL::IMAGE_COORD_X, patchX[i], j);
        MD.setValue(EMDL::IMAGE_COORD_Y, patchY[i], j);
        MD.setValue(EMDL::MICROGRAPH_SHIFT_X, localShiftX[i], j);
        MD.setValue(EMDL::MICROGRAPH_SHIFT_Y, localShiftY[i], j);
    }
    MD.write(fh);

}

FileName Micrograph::getGainFilename() const {
    return fnGain;
}

RFLOAT Micrograph::getBinningFactor() const {
    return binning;
}

FileName Micrograph::getMovieFilename() const {
    return fnMovie;
}

int Micrograph::getWidth() const {
    return width;
}

int Micrograph::getHeight() const {
    return height;
}

int Micrograph::getNframes() const {
    return n_frames;
}

int Micrograph::getEERUpsampling() const {
    return eer_upsampling;
}

int Micrograph::getEERGrouping() const {
    return eer_grouping;
}

void Micrograph::fillDefectAndHotpixels(MultidimArray<bool> &mask) const {
    checkReadyFlag("getShiftAt");

    mask.initZeros(height, width);

    bool fix_defect = (fnDefect != "");
    if (fnDefect.getExtension() == "txt" && MotioncorrRunner::detectSerialEMDefectText(fnDefect)) {
        std::cerr << "WARNING: The defect file specified in the micrograph metadata STAR file seems to be a SerialEM's defect file and not in the MotionCor2's format (x y w h). The defect file is ignored." << std::endl;
        fix_defect = false;
    }

    if (fix_defect) {
        MotioncorrRunner::fillDefectMask(mask, fnDefect);
    }

    if (hotpixelX.size() != hotpixelY.size()) {
        REPORT_ERROR("Logic error: hotpixelX.size() != hotpixelY.size()");
    }

    for (int i = 0, ilim = hotpixelX.size(); i < ilim; i++) {
        direct::elem(mask, hotpixelX[i], hotpixelY[i]) = true;
    }
}

int Micrograph::getShiftAt(RFLOAT frame, RFLOAT x, RFLOAT y, RFLOAT &shiftx, RFLOAT &shifty, bool use_local, bool normalise) const {
    checkReadyFlag("getShiftAt");

    if (normalise) {
        x = x / width - 0.5;
        y = y / height - 0.5;
    }

    if (globalShiftX[frame - 1] == NOT_OBSERVED || globalShiftX[frame - 1] == NOT_OBSERVED) {
        // Find the shift of the closest observed frame.
        // If the given 'frame' is unobserved due to initial frame truncation (--first_frame),
        // the output becomes zero. This is OK because the shift of the first observed frame
        // is zero by definition. So we don't have to search after the 'frame'.
        shiftx = shifty = 0;

        for (int i = frame - 1; i >= 0; i--) {
            if (globalShiftX[i] != NOT_OBSERVED && globalShiftY[i] != NOT_OBSERVED) {
                shiftx = globalShiftX[i];
                shifty = globalShiftY[i];
                break;
            }
        }
        return -1;
    }

    if (model && use_local) {
        // both frame and first_frame is 1 indexed
        model->getShiftAt(frame - first_frame, x, y, shiftx, shifty);
    } else {
        shiftx = 0;
        shifty = 0;
    }

    // frame is 1-indexed!
    shiftx += globalShiftX[frame - 1];
    shifty += globalShiftY[frame - 1];

    return 0;
}

void Micrograph::setGlobalShift(int frame, RFLOAT shiftx, RFLOAT shifty) {
    checkReadyFlag("setGlobalShift");

    if (frame <= 0 || frame > n_frames) {
        std::cout << "Frame: " << frame << " n_frames: " << n_frames << std::endl;
        REPORT_ERROR("Micrograph::setGlobalShift() frame out of range");
    }

    frame--; // frame is 1-indexed
    globalShiftX[frame] = shiftx;
    globalShiftY[frame] = shifty;
}

void Micrograph::read(FileName fn_in, bool read_hotpixels) {
    delete model;
    model = nullptr;

    // Clear current model
    clearFields();

    // Open input file
    std::ifstream in(fn_in.data(), std::ios_base::in);
    if (in.fail()) {
        REPORT_ERROR("MicrographModel::read: File " + fn_in + " cannot be read.");
    }

    MetaDataTable MDglobal, MDhot;

    // Read Image metadata
    MDglobal.readStar(in, "general");
    const long int i = MDglobal.size() - 1;

    try {
        width    = MDglobal.getValue<int>(EMDL::IMAGE_SIZE_X, i);
        height   = MDglobal.getValue<int>(EMDL::IMAGE_SIZE_Y, i);
        n_frames = MDglobal.getValue<int>(EMDL::IMAGE_SIZE_Z, i);
        fnMovie  = MDglobal.getValue<std::string>(EMDL::MICROGRAPH_MOVIE_NAME, i);
    } catch (const char *errmsg) {
        REPORT_ERROR("MicrographModel::read: insufficient general information in " + fn_in);
    }

    globalShiftX.resize(n_frames, NOT_OBSERVED);
    globalShiftY.resize(n_frames, NOT_OBSERVED);

    #define TRY_ASSIGN_FROM_MDT(assignmenttarget, table, T, label, defaultvalue) \
    try { (assignmenttarget) = (table).getValue<T>((label), i); } \
    catch (const char *errmsg) { (assignmenttarget) = (defaultvalue); }

    TRY_ASSIGN_FROM_MDT(fnGain,         MDglobal, std::string, EMDL::MICROGRAPH_GAIN_NAME,   "");
    TRY_ASSIGN_FROM_MDT(fnDefect,       MDglobal, std::string, EMDL::MICROGRAPH_DEFECT_FILE, "");
    TRY_ASSIGN_FROM_MDT(binning,        MDglobal, RFLOAT,      EMDL::MICROGRAPH_BINNING,     1.0);
    TRY_ASSIGN_FROM_MDT(angpix,         MDglobal, RFLOAT,      EMDL::MICROGRAPH_ORIGINAL_PIXEL_SIZE, -1);
    TRY_ASSIGN_FROM_MDT(pre_exposure,   MDglobal, RFLOAT,      EMDL::MICROGRAPH_PRE_EXPOSURE, -1);
    TRY_ASSIGN_FROM_MDT(dose_per_frame, MDglobal, RFLOAT,      EMDL::MICROGRAPH_DOSE_RATE, -1);
    TRY_ASSIGN_FROM_MDT(voltage,        MDglobal, RFLOAT,      EMDL::CTF_VOLTAGE, -1);
    TRY_ASSIGN_FROM_MDT(first_frame,    MDglobal, int,         EMDL::MICROGRAPH_START_FRAME, 1);  // 1-indexed

    if (EERRenderer::isEER(fnMovie)) {
    TRY_ASSIGN_FROM_MDT(eer_upsampling, MDglobal, int,      EMDL::MICROGRAPH_EER_UPSAMPLING, -1);
    TRY_ASSIGN_FROM_MDT(eer_grouping,   MDglobal, int,      EMDL::MICROGRAPH_EER_GROUPING,   -1);
    }

    #undef TRY_ASSIGN_FROM_MDT

    model = nullptr;
    try {
        const int model_version = MDglobal.getValue<int>(EMDL::MICROGRAPH_MOTION_MODEL_VERSION, i);
        if (model_version == MOTION_MODEL_THIRD_ORDER_POLYNOMIAL) {
            model = new ThirdOrderPolynomialModel();  // Where does this memory get deallocated?
        } else if (model_version == (int) MOTION_MODEL_NULL) {
            model = nullptr;
        } else {
            std::cerr << "Warning: Ignoring unknown motion model " << model_version << std::endl;
        }
    } catch (const char *errmsg) {
        std::cerr << "Warning: local motion model is absent in the micrograph star file." << std::endl;
    }

    if (model)
        model->read(in, "local_motion_model");

    // Read global shifts
    MDglobal.readStar(in, "global_shift");

    try {
        for (long int i : MDglobal) {
            const int frame = MDglobal.getValue<int>(EMDL::MICROGRAPH_FRAME_NUMBER, i);
            const RFLOAT shiftX = MDglobal.getValue<RFLOAT>(EMDL::MICROGRAPH_SHIFT_X, i);
            const RFLOAT shiftY = MDglobal.getValue<RFLOAT>(EMDL::MICROGRAPH_SHIFT_Y, i);
            // frame is 1-indexed!
            globalShiftX[frame - 1] = shiftX;
            globalShiftY[frame - 1] = shiftY;
        }
    } catch (const char *errmsg) {
        REPORT_ERROR("MicrographModel::read: incorrect global_shift table in " + fn_in);
    }

    if (read_hotpixels) {
        MDhot.readStar(in, "hot_pixels");
        try {
            for (long int i : MDhot) {
                const RFLOAT x = MDhot.getValue<RFLOAT>(EMDL::IMAGE_COORD_X, i);
                const RFLOAT y = MDhot.getValue<RFLOAT>(EMDL::IMAGE_COORD_Y, i);
                hotpixelX.push_back((int) x);
                hotpixelY.push_back((int) y);
            }
        } catch (const char *errmsg) {
            REPORT_ERROR("MicrographModel::read: incorrect hot_pixels table in " + fn_in);
        }
    }
}

void Micrograph::setMovie(FileName fnMovie, FileName fnGain, RFLOAT binning) {
    if (EERRenderer::isEER(fnMovie)) {
        EERRenderer renderer;
        renderer.read(fnMovie, eer_upsampling);
        width    = renderer.getWidth();
        height   = renderer.getHeight();
        n_frames = renderer.getNFrames() / eer_grouping;
    } else {
        Image<RFLOAT> Ihead;
        Ihead.read(fnMovie, false, -1, nullptr, true); // select_img -1, mmap false, is_2D true

        width    = Xsize(Ihead());
        height   = Ysize(Ihead());
        n_frames = Nsize(Ihead());
    }

    this->binning = binning;

    globalShiftX.resize(n_frames, NOT_OBSERVED);
    globalShiftY.resize(n_frames, NOT_OBSERVED);

    this->fnMovie = fnMovie;
    this->fnGain = fnGain;
}

void Micrograph::clearFields() {
    width = 0;
    height = 0;
    n_frames = 0;
    first_frame = 0;
    binning = 1;

    angpix = -1;
    voltage = -1;
    dose_per_frame = -1;
    pre_exposure = -1;

    eer_upsampling = -1;
    eer_grouping = -1;

    fnMovie = "";
    fnGain = "";
    fnDefect = "";

    hotpixelX.resize(0);
    hotpixelY.resize(0);

    globalShiftX.resize(0);
    globalShiftY.resize(0);

    localShiftX.resize(0);
    localShiftY.resize(0);
    localFitX.resize(0);
    localFitY.resize(0);
    patchX.resize(0);
    patchY.resize(0);
    patchZ.resize(0);
    patchW.resize(0);
    patchH.resize(0);
}

void Micrograph::copyFieldsFrom(const Micrograph& m) {
    width = m.width;
    height = m.height;
    n_frames = m.n_frames;
    first_frame = m.first_frame;
    binning = m.binning;

    angpix = m.angpix;
    voltage = m.voltage;
    dose_per_frame = m.dose_per_frame;
    pre_exposure = m.pre_exposure;

    eer_upsampling = m.eer_upsampling;
    eer_grouping = m.eer_grouping;

    fnMovie = m.fnMovie;
    fnGain = m.fnGain;
    fnDefect = m.fnDefect;

    hotpixelX = m.hotpixelX;
    hotpixelY = m.hotpixelY;

    globalShiftX = m.globalShiftX;
    globalShiftY = m.globalShiftY;

    localShiftX = m.localShiftX;
    localShiftY = m.localShiftY;
    localFitX = m.localFitX;
    localFitY = m.localFitY;
    patchX = m.patchX;
    patchY = m.patchY;
    patchZ = m.patchZ;
    patchW = m.patchW;
    patchH = m.patchH;
}

void Micrograph::checkReadyFlag(std::string origin) const {
    if (!ready) {
        REPORT_ERROR("Micrograph::"+origin+": instance not initialized.\n");
    }
}
