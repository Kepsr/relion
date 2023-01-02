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

#include "defocus_helper.h"

#include <src/jaz/slice_helper.h>
#include <src/jaz/img_proc/filter_helper.h>
#include <src/jaz/refinement_helper.h>
#include <src/jaz/optimization/nelder_mead.h>
#include <src/jaz/gravis/t4Matrix.h>
#include "src/jaz/ctf_helper.h"

#include <src/projector.h>

using namespace gravis;

RFLOAT DefocusHelper::findDefocus1D(
        const Image<Complex> &prediction,
        const Image<Complex> &observation,
        const Image<RFLOAT>& weight,
        const CTF &ctf0, ObservationModel *obsModel, int opticsGroup,
        RFLOAT angpix,
        double *destU, double *destV,
        RFLOAT range, int steps,
        int recDepth, RFLOAT recScale)
{
    const RFLOAT delta = ctf0.DeltafV - ctf0.DeltafU;
    CTF ctf(ctf0);

    double minErr = std::numeric_limits<double>::max();
    double bestU = ctf0.DeltafU;
    double bestV = ctf0.DeltafV;

    for (int s = -steps/2; s <= steps/2; s++)
    {
        const RFLOAT u = ctf0.DeltafU + s * range / (steps/2);
        const RFLOAT v = u + delta;

        ctf.DeltafU = u;
        ctf.DeltafV = v;
        ctf.initialise();

        double err = RefinementHelper::squaredDiff(
            prediction, observation, ctf, obsModel, opticsGroup,
            angpix, weight
        );

        if (err < minErr)
        {
            minErr = err;
            bestU = u;
            bestV = v;
        }
    }

    if (recDepth > 0)
    {
        ctf.DeltafU = bestU;
        ctf.DeltafV = bestV;
        ctf.initialise();

        findDefocus1D(prediction, observation, weight,
                ctf, obsModel, opticsGroup, angpix, destU, destV,
                range/recScale, steps,
                recDepth-1, recScale);
    }
    else
    {
        *destU = bestU;
        *destV = bestV;
    }

    return minErr;
}

void DefocusHelper::findAstigmatismNM(
        const Image<Complex> &prediction,
        const Image<Complex> &observation,
        const Image<RFLOAT> &weight,
        const CTF &ctf0, ObservationModel *obsModel, int opticsGroup,
        RFLOAT angpix,
        double *destU, double *destV, double *destPhi)
{
    AstigmatismOptimizationAcc opt(prediction, observation, weight, ctf0, obsModel, opticsGroup, false, false, angpix);

    std::vector<double> initial = opt.getInitialParams();

    std::vector<double> params = NelderMead::optimize(initial, opt, 50.0, 1.0, 1000);

    *destU = opt.getU(params);
    *destV = opt.getV(params);
    *destPhi = opt.getPhi(params);
}

void DefocusHelper::findAstigmatismAndPhaseNM(
        const std::vector<Image<Complex>>& prediction,
        const std::vector<Image<Complex>>& observation,
        const Image<RFLOAT> &weight,
        const CTF &ctf0, ObservationModel *obsModel, int opticsGroup,
        RFLOAT angpix,
        double *destU, double *destV, double *destPhi, double *destPhase)
{
    AstigmatismOptimizationAcc opt(prediction, observation, weight, ctf0, obsModel, true, false, angpix);

    std::vector<double> initial = opt.getInitialParams();

    std::vector<double> params = NelderMead::optimize(initial, opt, 5.0, 0.01, 100000);

    *destU     = opt.getU(params);
    *destV     = opt.getV(params);
    *destPhi   = opt.getPhi(params);
    *destPhase = opt.getPhase(params);
}

void DefocusHelper::findAstigmatismPhaseAndCsNM(
        const std::vector<Image<Complex>>& prediction,
        const std::vector<Image<Complex>>& observation,
        const Image<RFLOAT> &weight,
        const CTF &ctf0, ObservationModel *obsModel, int opticsGroup,
        RFLOAT angpix,
        double *destU, double *destV, double *destPhi, double *destPhase, double* destCs)
{
    AstigmatismOptimizationAcc opt(prediction, observation, weight, ctf0, obsModel, true, true, angpix);

    std::vector<double> initial = opt.getInitialParams();

    std::vector<double> params = NelderMead::optimize(initial, opt, 5.0, 0.01, 100000);

    *destU     = opt.getU(params);
    *destV     = opt.getV(params);
    *destPhi   = opt.getPhi(params);
    *destPhase = opt.getPhase(params);
    *destCs    = opt.getCs(params);
}

void DefocusHelper::findAstigmatismNM(
        const std::vector<Image<Complex>>& prediction,
        const std::vector<Image<Complex>>& observation,
        const Image<RFLOAT> &weight,
        const CTF &ctf0, ObservationModel *obsModel, int opticsGroup,
        RFLOAT angpix,
        double *destU, double *destV, double *destPhi)
{
    AstigmatismOptimizationAcc opt(prediction, observation, weight, ctf0, obsModel, false, false, angpix);

    std::vector<double> initial = opt.getInitialParams();

    std::vector<double> params = NelderMead::optimize(initial, opt, 50.0, 1.0, 1000);

    *destU = opt.getU(params);
    *destV = opt.getV(params);
    *destPhi = opt.getPhi(params);
}

std::vector<d2Vector> DefocusHelper::diagnoseDefocus(
        const Image<Complex> &prediction,
        const Image<Complex> &observation,
        const Image<RFLOAT> &weight,
        const CTF &ctf0, ObservationModel *obsModel, int opticsGroup,
        RFLOAT angpix,
        double range, int steps, int threads)
{
    const RFLOAT delta = ctf0.DeltafV - ctf0.DeltafU;

    std::vector<d2Vector> out(steps);

    #pragma omp parallel for num_threads(threads)
    for (int s = 0; s < steps; s++)
    {
        CTF ctf(ctf0);
        const RFLOAT u = ctf0.DeltafU + (s - steps/2) * range / (double)steps;
        const RFLOAT v = u + delta;

        ctf.DeltafU = u;
        ctf.DeltafV = v;
        ctf.initialise();

        out[s][0] = u;
        out[s][1] = RefinementHelper::squaredDiff(
            prediction, observation, ctf, obsModel, opticsGroup, angpix, weight
        );
    }

    return out;
}

AstigmatismOptimization::AstigmatismOptimization(
    const Image<Complex>& prediction, const Image<Complex>& observation,
    const Image<RFLOAT>& weight, const CTF& ctf0, ObservationModel *obsModel, int opticsGroup, RFLOAT angpix
):
    prediction(prediction), observation(observation),
    weight(weight), ctf0(ctf0), obsModel(obsModel), opticsGroup(opticsGroup), angpix(angpix) {}

double AstigmatismOptimization::f(const std::vector<double>& x) const
{
    CTF ctf(ctf0);

    ctf.DeltafU = x[0];
    ctf.DeltafV = x[1];
    ctf.azimuthal_angle = x[2];
    ctf.initialise();

    return RefinementHelper::squaredDiff(prediction, observation, ctf, obsModel, opticsGroup, angpix, weight);
}

AstigmatismOptimizationAcc::AstigmatismOptimizationAcc(
        const Image<Complex>& prediction,
        const Image<Complex>& observation,
        const Image<RFLOAT>& weight,
        const CTF& ctf0,
        ObservationModel *obsModel,
        int opticsGroup,
        bool phaseShift,
        bool spherAberr,
        RFLOAT angpix,
        RFLOAT phiScale,
        RFLOAT csScale)
:   ctf0(ctf0),
    obsModel(obsModel),
    phaseShift(phaseShift),
    spherAberr(spherAberr),
    angpix(angpix),
    phiScale(phiScale),
    csScale(csScale)
{
    const long w = prediction.data.xdim;
    const long h = prediction.data.ydim;

    data = Image<Complex>(w,h);

    for (long y = 0; y < h; y++)
    for (long x = 0; x < w; x++) {
        Complex vx = direct::elem(prediction.data, x, y);
        const Complex vy = direct::elem(observation.data, x, y);
        const RFLOAT vw = direct::elem(weight.data, x, y);

        const RFLOAT x2 = vx.real*vx.real + vx.imag*vx.imag;

        const RFLOAT yxb = x2 > 0.0? (vy.real*vx.real + vy.imag*vx.imag)/x2 : 0.0;
        const RFLOAT wp = vw * x2;

        direct::elem(data.data, x, y) = Complex(yxb, wp);
    }
}

// Complex dot product
inline RFLOAT dot(Complex a, Complex b) {
    return a.real * b.real + a.imag * b.imag;
}

AstigmatismOptimizationAcc::AstigmatismOptimizationAcc(
    const std::vector<Image<Complex>>& prediction,
    const std::vector<Image<Complex>>& observation,
    const Image<RFLOAT>& weight,
    const CTF& ctf0,
    ObservationModel *obsModel,
    bool phaseShift, bool spherAberr,
    RFLOAT angpix, RFLOAT phiScale, RFLOAT csScale
): 
ctf0(ctf0), obsModel(obsModel), phaseShift(phaseShift), spherAberr(spherAberr), 
angpix(angpix), phiScale(phiScale), csScale(csScale) {
    const long w = prediction[0].data.xdim;
    const long h = prediction[0].data.ydim;

    data = Image<Complex>::zeros(w, h);

    for (int i = 0; i < prediction.size(); i++) {
        for (long y = 0; y < h; y++)
        for (long x = 0; x < w; x++) {
                  Complex vx = direct::elem(prediction[i].data, x, y);
            const Complex vy = direct::elem(observation[i].data, x, y);
            const RFLOAT vw = direct::elem(weight.data, x, y);

            const RFLOAT x2 = vx.real * vx.real + vx.imag * vx.imag;

            const RFLOAT yxb = x2 > 0.0 ? dot(vx, vy) / x2 : 0.0;
            const RFLOAT wp = vw * x2;

            direct::elem(data.data, x, y) += Complex(yxb, wp);
        }
    }
}

double AstigmatismOptimizationAcc::f(
    const std::vector<double> &x, void* tempStorage
) const {
    CTF ctf (ctf0);

    ctf.DeltafU = x[0];
    ctf.DeltafV = x[1];
    ctf.azimuthal_angle = x[2] / phiScale;

    if (phaseShift) ctf.phase_shift = x[3] / phiScale;
    if (spherAberr) ctf.Cs = x[3 + phaseShift] / csScale;

    ctf.initialise();

    const long w = data.data.xdim;
    const long h = data.data.ydim;

    Image<RFLOAT> ctfImg(w, h);
    ctfImg() = CtfHelper::getFftwImage(ctf, w, h, h, h, angpix, obsModel, opticsGroup);

    double out = 0.0;

    for (long y = 0; y < h; y++)
    for (long x = 0; x < w; x++) {
        Complex vd = direct::elem(data.data, x, y);

        RFLOAT vm = ctfImg(y, x);
        RFLOAT dx = vd.real - vm;
        out += vd.imag * dx * dx;
    }
    return out;
}

double AstigmatismOptimizationAcc::getU(const std::vector<double> &x) {
    return x[0];
}

double AstigmatismOptimizationAcc::getV(const std::vector<double> &x) {
    return x[1];
}

double AstigmatismOptimizationAcc::getPhi(const std::vector<double> &x) {
    return x[2] / phiScale;
}

double AstigmatismOptimizationAcc::getPhase(const std::vector<double> &x) {
    return x[3] / phiScale;
}

double AstigmatismOptimizationAcc::getCs(const std::vector<double> &x) {
    return x[phaseShift ? 4 : 3] / csScale;
}

std::vector<double> AstigmatismOptimizationAcc::getInitialParams() {
    std::vector<double> params;
    params.reserve(3 + phaseShift + spherAberr);
    params = {ctf0.DeltafU, ctf0.DeltafV, ctf0.azimuthal_angle * phiScale};
    if (phaseShift) params.push_back(ctf0.phase_shift * phiScale);
    if (spherAberr) params.push_back(ctf0.Cs * csScale);
    return params;
}
