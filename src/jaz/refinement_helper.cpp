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

#include <src/jaz/refinement_helper.h>
#include <src/jaz/slice_helper.h>
#include <src/projector.h>
#include <src/jaz/img_proc/filter_helper.h>
#include <src/jaz/optimization/nelder_mead.h>
#include <src/jaz/gravis/t4Matrix.h>
#include <src/jaz/vtk_helper.h>

namespace constrain {

    void upper(int &x, int ceil) {
        if (x > ceil) { x = ceil; }
        // x = min(x, ceil);
    }

    void lower(int &x, int floor) {
        if (x < floor) { x = floor; }
        // x = max(x, floor);
    }

}

using namespace gravis;

void RefinementHelper::drawFSC(
    const MetaDataTable *mdt, std::vector<double> &dest1D,
    Image<RFLOAT> &dest, double thresh
) {
    const int n = mdt->numberOfObjects();
    const int w = 2 * (n - 1);
    const int h = 2 * (n - 1);

    dest1D = std::vector<double>(n);

    for (int i = 0; i < n; i++) {
        int idx   = mdt->getValue<int>   (EMDL::SPECTRAL_IDX,         i);
        dest1D[i] = mdt->getValue<double>(EMDL::POSTPROCESS_FSC_TRUE, i);
        if (dest1D[i] < thresh) { dest1D[i] = 0.0; }
    }

    dest = Image<RFLOAT>(n, h);

    for (int y = 0; y < h; y++)
    for (int x = 0; x < n; x++) {
        double xx = x;
        double yy = y > h / 2 ? y - h : y;

        double r = sqrt(xx * xx + yy * yy);
        int ri = (int) (r + 0.5);
        if (ri > w / 2) { ri = w / 2; }

        direct::elem(dest.data, x, y) = dest1D[ri];
    }
}

void RefinementHelper::computeSNR(const MetaDataTable *mdt, Image<RFLOAT> &dest, double eps) {
    const int n = mdt->numberOfObjects();
    const int w = 2 * (n - 1);
    const int h = 2 * (n - 1);

    std::vector<double> snr(n);

    for (int i = 0; i < n; i++) {
        int    idx = mdt->getValue<int>   (EMDL::SPECTRAL_IDX,         i);
        double fsc = mdt->getValue<double>(EMDL::POSTPROCESS_FSC_TRUE, i);

        if (fsc > 1.0 - eps) { fsc = 1.0 - eps; }
        //else if (fsc < eps) fsc = 0.0;

        snr[i] = fsc / (1.0 - fsc);
    }

    dest = Image<RFLOAT>(n, h);

    for (int y = 0; y < h; y++)
    for (int x = 0; x < n; x++) {
        double xx = x;
        double yy = y > h / 2 ? y - h : y;

        double r = sqrt(xx * xx + yy * yy);
        int ri = (int) (r + 0.5);
        if (ri > w / 2) { ri = w / 2; }

        direct::elem(dest.data, x, y) = snr[ri];
    }
}

void RefinementHelper::computeSigInvSq(
    const MetaDataTable *mdt, const std::vector<double> &signalPow,
    Image<RFLOAT> &dest, double eps
) {
    const int n = mdt->numberOfObjects();
    const int w = 2 * (n - 1);
    const int h = 2 * (n - 1);

    std::vector<double> sigInvSq(n);

    for (int i = 0; i < n; i++) {
        int    idx = mdt->getValue<int>(EMDL::SPECTRAL_IDX, i);
        double fsc = mdt->getValue<double>(EMDL::POSTPROCESS_FSC_TRUE, i);

        if (fsc > 1.0 - eps) { fsc = 1.0 - eps; }
        // else if (fsc < eps) { fsc = 0.0; }

        double snr = fsc / (1.0 - fsc);
        double sigPow = signalPow[i];
        if (sigPow < eps) { sigPow = eps; }

        sigInvSq[i] = snr / sigPow;
    }

    dest = Image<RFLOAT>(n,h);

    for (int y = 0; y < h; y++)
    for (int x = 0; x < n; x++) {
        double xx = x;
        double yy = y > h / 2 ? y - h : y;

        double r = sqrt(xx * xx + yy * yy);
        int ri = (int) (r + 0.5);
        if (ri > w / 2) { ri = w / 2; }

        direct::elem(dest.data, x, y) = sigInvSq[ri];
    }
}

Image<RFLOAT> RefinementHelper::correlation(
    const Image<Complex> &prediction, const Image<Complex> &observation
) {
    const long w = prediction.data.xdim;
    const long h = prediction.data.ydim;

    Image<RFLOAT> out(w, h);

    for (long y = 0; y < h; y++)
    for (long x = 0; x < w; x++) {
        Complex vx = direct::elem(prediction .data, x, y);
        Complex vy = direct::elem(observation.data, x, y);

        // Dot product
        direct::elem(out.data, x, y) = vy.real * vx.real + vy.imag * vx.imag;
    }

    return out;
}

// Basically, map correlation
Image<RFLOAT> RefinementHelper::correlation(
    const std::vector<Image<Complex>> &predictions,
    const std::vector<Image<Complex>> &observations
) {
    const long w = predictions[0].data.xdim;
    const long h = predictions[0].data.ydim;
    const long c = predictions.size();

    Image<RFLOAT> out = Image<RFLOAT>::zeros(w, h);

    for (long i = 0; i < c; i++) {
        for (long y = 0; y < h; y++)
        for (long x = 0; x < w; x++) {
            Complex vx = direct::elem(predictions [i].data, x, y);
            Complex vy = direct::elem(observations[i].data, x, y);

            // Dot product
            direct::elem(out.data, y, x) += vy.real * vx.real + vy.imag * vx.imag;
        }
    }

    return out;
}

void RefinementHelper::addToQR(
    const Image<Complex> &prediction,
    const Image<Complex> &observation,
    Image<Complex> &q, Image<RFLOAT> &r
) {
    const long w = prediction.data.xdim;
    const long h = prediction.data.ydim;

    for (long y = 0; y < h; y++)
    for (long x = 0; x < w; x++) {
        Complex vx = direct::elem(prediction .data, x, y);
        Complex vy = direct::elem(observation.data, x, y);

        direct::elem(q.data, x, y) += vy.conj() * vx;
        direct::elem(r.data, x, y) += vx.norm();
    }
}

void RefinementHelper::addToPQR(
    const Image<Complex> &prediction,
    const Image<Complex> &observation,
    Image<RFLOAT> &p, Image<Complex> &q, Image<RFLOAT> &r
) {
    const long w = prediction.data.xdim;
    const long h = prediction.data.ydim;

    for (long y = 0; y < h; y++)
    for (long x = 0; x < w; x++) {
        Complex vx = direct::elem(prediction.data,  x, y);
        Complex vy = direct::elem(observation.data, x, y);

        direct::elem(p.data, x, y) += vx.norm();
        direct::elem(q.data, x, y) += vy.conj() * vx;
        direct::elem(r.data, x, y) += vx.norm();
    }
}

double RefinementHelper::squaredDiff(
    const Image<Complex> &prediction,
    const Image<Complex> &observation,
    CTF &ctf, ObservationModel *obsModel,
    RFLOAT angpix, const Image<RFLOAT> &weight
) {
    const long w = prediction.data.xdim;
    const long h = prediction.data.ydim;

    Image<RFLOAT> ctfImg(w, h);
    ctfImg() = ctf.getFftwImage(w, h, h, h, angpix, obsModel);

    double out = 0.0;
    for (long y = 0; y < h; y++)
    for (long x = 0; x < w; x++) {
        Complex vx = direct::elem(prediction.data, x, y);
        const Complex vy = direct::elem(observation.data, x, y);
        const RFLOAT  vw = direct::elem(weight     .data, x, y);

        RFLOAT vm = ctfImg(y, x);
        out += vw * (vy - vm * vx).norm();
    }
    return out;
}

double RefinementHelper::squaredDiff(
    const std::vector<Image<Complex>> &predictions,
    const std::vector<Image<Complex>> &observations,
    CTF &ctf, ObservationModel *obsModel, RFLOAT angpix, const Image<RFLOAT> &weight
) {
    double out = 0.0;
    for (long i = 0; i < predictions.size(); i++) {
        out += squaredDiff(predictions[i], observations[i], ctf, obsModel, angpix, weight);
    }
    return out;
}
