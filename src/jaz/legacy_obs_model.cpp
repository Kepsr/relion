#include "src/jaz/legacy_obs_model.h"
#include "src/jaz/stack_helper.h"
#include "src/jaz/img_proc/filter_helper.h"
#include "src/jaz/ctf_helper.h"
#include "src/jaz/Fourier_helper.h"

#include <src/backprojector.h>

LegacyObservationModel::LegacyObservationModel(): angpix(-1), anisoTilt(false) {}

LegacyObservationModel::LegacyObservationModel(double angpix, double Cs, double voltage): 
angpix(angpix), 
lambda(12.2643247 / sqrt(voltage * (1.0 + voltage * 0.978466e-6))),
Cs(Cs),
anisoTilt(false)
{}

void LegacyObservationModel::predictObservation(
    Projector& proj, const MetaDataTable &mdt, int particle,
    MultidimArray<Complex> &dest,
    bool applyCtf, bool applyTilt, bool applyShift
) const {

    const int s = proj.ori_size;
    const int sh = s / 2 + 1;

    double xoff = mdt.getValue<double>(EMDL::ORIENT_ORIGIN_X, particle);
    double yoff = mdt.getValue<double>(EMDL::ORIENT_ORIGIN_Y, particle);

    double rot  = mdt.getValue<double>(EMDL::ORIENT_ROT,  particle);
    double tilt = mdt.getValue<double>(EMDL::ORIENT_TILT, particle);
    double psi  = mdt.getValue<double>(EMDL::ORIENT_PSI,  particle);

    Matrix<RFLOAT> A3D = Euler::angles2matrix(rot, tilt, psi);

    dest = proj.get2DFourierTransform(sh, s, 1, A3D);

    if (applyShift) {
        shiftImageInFourierTransform(dest, s, s / 2 - xoff, s / 2 - yoff);
    }

    if (applyCtf) {
        CTF ctf = CtfHelper::makeCTF(mdt, particle);
        FilterHelper::modulate(dest, ctf, nullptr, -1, angpix);
    }

    if (applyTilt) {
        double tx = 0.0, ty = 0.0;
        tx = mdt.getValue<double>(EMDL::IMAGE_BEAMTILT_X, particle);
        ty = mdt.getValue<double>(EMDL::IMAGE_BEAMTILT_Y, particle);

        if (tx != 0.0 && ty != 0.0) {
            if (anisoTilt) {
                selfApplyBeamTilt(
                    dest, -tx, -ty, 
                    beamtilt_xx, beamtilt_xy, beamtilt_yy, 
                    lambda, Cs, angpix, s
                );
            } else {
                selfApplyBeamTilt(
                    dest, -tx, -ty, 
                    lambda, Cs, angpix, s
                );
            }
        }
    }
}


Image<Complex> LegacyObservationModel::predictObservation(
    Projector& proj, const MetaDataTable &mdt, int particle,
    bool applyCtf, bool applyTilt, bool applyShift
) const {
    MultidimArray<Complex> pred;
    predictObservation(proj, mdt, particle, pred, applyCtf, applyTilt, applyShift);
    return Image<Complex>(pred);
}

std::vector<Image<Complex>> LegacyObservationModel::predictObservations(
    Projector &proj, const MetaDataTable &mdt, int threads,
    bool applyCtf, bool applyTilt, bool applyShift
) const {
    const int pc = mdt.size();
    std::vector<Image<Complex>> out (pc);

    #pragma omp parallel for num_threads(threads)
    for (int p = 0; p < pc; p++) {
        out[p] = predictObservation(proj, mdt, p, applyCtf, applyTilt, applyShift);
    }

    return out;
}

void LegacyObservationModel::insertObservation(
    const Image<Complex> &img, BackProjector &bproj,
    const MetaDataTable &mdt, int particle,
    bool applyCtf, bool applyTilt, double shift_x, double shift_y
) {
    const int sh = img.data.xdim, s = img.data.ydim;

    const RFLOAT rot  = mdt.getValue<RFLOAT>(EMDL::ORIENT_ROT,  particle);
    const RFLOAT tilt = mdt.getValue<RFLOAT>(EMDL::ORIENT_TILT, particle);
    const RFLOAT psi  = mdt.getValue<RFLOAT>(EMDL::ORIENT_PSI,  particle);

    const Matrix<RFLOAT> A3D = Euler::angles2matrix(rot, tilt, psi);

    double tx = 0.0, ty = 0.0;
    tx = mdt.getValue<double>(EMDL::ORIENT_ORIGIN_X, particle) + shift_x;
    ty = mdt.getValue<double>(EMDL::ORIENT_ORIGIN_Y, particle) + shift_y;

    MultidimArray<Complex> F2D = img.data;

    shiftImageInFourierTransform(F2D, s, tx, ty);

    auto Fctf = MultidimArray<RFLOAT>::ones(F2D.xdim, F2D.ydim, F2D.zdim, F2D.ndim);

    if (applyCtf) {
        CTF ctf = CtfHelper::makeCTF(mdt, particle);
        Fctf = CtfHelper::getFftwImage(ctf, sh, s, s, s, angpix, nullptr, -1);

        for (long int n = 0; n < F2D.size(); n++) {
            F2D[n]  *= Fctf[n];
            Fctf[n] *= Fctf[n];
        }
    }

    if (applyTilt) {

        const double my_tilt_x = mdt.containsLabel(EMDL::IMAGE_BEAMTILT_X) ?
            mdt.getValue<double>(EMDL::IMAGE_BEAMTILT_X, particle) : 0.0;

        const double my_tilt_y = mdt.containsLabel(EMDL::IMAGE_BEAMTILT_Y) ?
            mdt.getValue<double>(EMDL::IMAGE_BEAMTILT_Y, particle) : 0.0;

        selfApplyBeamTilt(F2D, my_tilt_x, my_tilt_y, lambda, Cs, angpix, sh);
    }

    bproj.set2DFourierTransform(F2D, A3D, &Fctf);
}

void LegacyObservationModel::setAnisoTilt(double xx, double xy, double yy) {
    beamtilt_xx = xx;
    beamtilt_xy = xy;
    beamtilt_yy = yy;
    anisoTilt = true;
}

double LegacyObservationModel::angToPix(double a, int s) {
    return s * angpix / a;
}

double LegacyObservationModel::pixToAng(double p, int s) {
    return s * angpix / p;
}

bool LegacyObservationModel::containsAllNeededColumns(const MetaDataTable &mdt) {
    return mdt.containsLabel(EMDL::ORIENT_ORIGIN_X) &&
           mdt.containsLabel(EMDL::ORIENT_ORIGIN_Y) &&
           mdt.containsLabel(EMDL::ORIENT_ROT) &&
           mdt.containsLabel(EMDL::ORIENT_TILT) &&
           mdt.containsLabel(EMDL::ORIENT_PSI) &&
           mdt.containsLabel(EMDL::PARTICLE_RANDOM_SUBSET);
}
