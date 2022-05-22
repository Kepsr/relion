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

#include <src/jaz/fsc_helper.h>
#include <src/jaz/gravis/t2Matrix.h>
#include <src/jaz/optimization/nelder_mead.h>

using namespace gravis;

/**
 * @brief Divide dividend by divisor, unless divisor is zero, in which case return zero.
 * 
 * @tparam T 
 * @param dividend 
 * @param divisor 
 * @return T 
 */
template <typename T>
static T safedivide(T dividend, T divisor) {
    if (divisor == 0.0) return 0.0;
    return dividend / divisor;
}

void FscHelper::computeFscTable(
    const std::vector<std::vector<Image<Complex> > > &frames,
    const std::vector<Image<Complex> > &predictions,
    Image<RFLOAT> &table, Image<RFLOAT> &weight
) {
    const int w = predictions[0].data.xdim;
    const int pc = frames.size();
    const int fc = frames[0].size();

    table  = Image<RFLOAT>::zeros(w, fc);
    weight = Image<RFLOAT>::zeros(w, fc);

    Image<RFLOAT> weight1 = Image<RFLOAT>::zeros(w, fc);
    Image<RFLOAT> weight2 = Image<RFLOAT>::zeros(w, fc);

    for (int p = 0; p < pc; p++)
    for (int f = 0; f < fc; f++) {
        FOR_ALL_ELEMENTS_IN_FFTW_TRANSFORM(frames[p][f]()) {
            int idx = round(sqrt(kp * kp + ip * ip + jp * jp));

            if (idx >= w) {
                continue;
            }

            Complex z1 = direct::elem(frames[p][f](),   i, j, k);
            Complex z2 = direct::elem(predictions[p](), i, j, k);

            direct::elem(table.data,   idx, f) += z1.real * z2.real + z1.imag * z2.imag;
            direct::elem(weight1.data, idx, f) += z1.norm();
            direct::elem(weight2.data, idx, f) += z2.norm();
        }
    }

    for (int f = 0; f < fc; f++)
    for (int x = 0; x < w;  x++) {
        RFLOAT w1 = direct::elem(weight1.data, x, f);
        RFLOAT w2 = direct::elem(weight2.data, x, f);
        RFLOAT ww = sqrt(w1 * w2);

        direct::elem(weight.data, x, f)  = ww;
        direct::elem(table.data,  x, f) /= ww;
    }
}

void FscHelper::computeFscRow(
    const MultidimArray<Complex> &data0, const MultidimArray<Complex> &data1,
    int row, Image<RFLOAT> &table, Image<RFLOAT> &weight
) {
    const int w = data0.xdim;

    std::vector<double> weight1(w, 0.0), weight2(w, 0.0), data(w, 0.0);

    FOR_ALL_ELEMENTS_IN_FFTW_TRANSFORM(data0) {
        int idx = round(sqrt(kp * kp + ip * ip + jp * jp));

        if (idx >= w) {
            continue;
        }

        Complex z1 = direct::elem(data0, i, j, k);
        Complex z2 = direct::elem(data1, i, j, k);

        data[idx] += z1.real * z2.real + z1.imag * z2.imag;
        weight1[idx] += z1.norm();
        weight2[idx] += z2.norm();
    }

    for (int x = 0; x < w; x++) {
        if (x >= table.data.xdim) continue;

        RFLOAT ww = sqrt(weight1[x] * weight2[x]);

        direct::elem(table.data,  x, row) = safedivide(data[x], ww);
        direct::elem(weight.data, x, row) = ww;
    }
}

void FscHelper::initFscTable(
    int kc, int tc, Image<RFLOAT> &table, 
    Image<RFLOAT> &weight0, Image<RFLOAT> &weight1
) {
    table   = Image<RFLOAT>::zeros(kc, tc);
    weight0 = Image<RFLOAT>::zeros(kc, tc);
    weight1 = Image<RFLOAT>::zeros(kc, tc);
}

void FscHelper::updateFscTable(
    const std::vector<Image<Complex> > &frames,
    const Image<Complex> &prediction, double scale,
    Image<RFLOAT> &table,
    Image<RFLOAT> &weight0,
    Image<RFLOAT> &weight1
) {
    const int fc = frames.size();

    for (int f = 0; f < fc; f++) {
        updateFscTable(frames[f], f, prediction, scale, table, weight0, weight1);
    }
}

void FscHelper::updateFscTable(
    const Image<Complex> &frame, int f,
    const Image<Complex> &prediction, double scale,
    Image<RFLOAT> &table,
    Image<RFLOAT> &weight0, Image<RFLOAT> &weight1
) {
    const int w = prediction.data.xdim;

    FOR_ALL_ELEMENTS_IN_FFTW_TRANSFORM(frame()) {
        int idx = round(scale * sqrt(kp * kp + ip * ip + jp * jp));

        if (idx >= w) {
            continue;
        }

        Complex z1 = direct::elem(frame(),      i, j, k);
        Complex z2 = direct::elem(prediction(), i, j, k);

        direct::elem(table.data,   idx, f) += z1.real * z2.real + z1.imag * z2.imag;
        direct::elem(weight0.data, idx, f) += z1.norm();
        direct::elem(weight1.data, idx, f) += z2.norm();
    }
}

void FscHelper::updateFscTableVelWgh(
    const std::vector<Image<Complex> > &frames,
    const std::vector<d2Vector> &velocities,
    const Image<Complex> &prediction,
    Image<RFLOAT> &table,
    Image<RFLOAT> &weight0,
    Image<RFLOAT> &weight1
) {
    const int w = prediction.data.xdim;
    const int fc = frames.size();

    for (int f = 0; f < fc; f++) {
        FOR_ALL_ELEMENTS_IN_FFTW_TRANSFORM(frames[f]()) {
            int idx = round(sqrt(kp * kp + ip * ip + jp * jp));

            if (idx >= w) {
                continue;
            }

            double kv = (ip * velocities[f].y + jp * velocities[f].x)/(double)w;
            double wgh = kv < 1e-20 ? 1.0 : sin(PI * kv) / PI * kv;
            // double wgh = exp(-0.5*kv*kv/0.5);

            Complex z1 = direct::elem(frames[f](), i, j, k);
            Complex z2 = wgh * direct::elem(prediction(), i, j, k);

            direct::elem(table.data,   idx, f) += z1.real * z2.real + z1.imag * z2.imag;
            direct::elem(weight0.data, idx, f) += z1.norm();
            direct::elem(weight1.data, idx, f) += z2.norm();
        }
    }
}

void FscHelper::updateVelFscTable(
    const std::vector<Image<Complex> > &frames,
    const std::vector<d2Vector> &velocities,
    const Image<Complex> &prediction,
    Image<RFLOAT> &table,
    Image<RFLOAT> &weight0,
    Image<RFLOAT> &weight1,
    int kmin, int kmax
) {
    const int w = prediction.data.xdim;
    const int fc = frames.size();

    if (kmax < 0) { kmax = w; }

    for (int f = 0; f < fc; f++) {
        FOR_ALL_ELEMENTS_IN_FFTW_TRANSFORM(frames[f]()) {
            int idx = round(sqrt(kp * kp + ip * ip + jp * jp));

            if (idx >= w || idx < kmin || idx > kmax) {
                continue;
            }

            double kv = ip * velocities[f].y + jp * velocities[f].x;
            int kvidx = round(kv);

            if (kvidx < 0) { kvidx = -kvidx; }
            if (kvidx >= table().xdim) continue;

            Complex z1 = direct::elem(frames[f](),  i, j, k);
            Complex z2 = direct::elem(prediction(), i, j, k);

            direct::elem(table.data,   kvidx, f) += z1.real * z2.real + z1.imag * z2.imag;
            direct::elem(weight0.data, kvidx, f) += z1.norm();
            direct::elem(weight1.data, kvidx, f) += z2.norm();
        }
    }
}

void FscHelper::mergeFscTables(
    const std::vector<Image<RFLOAT> > &tables,
    const std::vector<Image<RFLOAT> > &weights0,
    const std::vector<Image<RFLOAT> > &weights1,
    Image<RFLOAT> &table, Image<RFLOAT> &weight
) {
    const int w = tables[0].data.xdim;
    const int fc = tables[0].data.ydim;
    const int mgc = tables.size();

    Image<RFLOAT> tableSum   = Image<RFLOAT>::zeros(w, fc);
    Image<RFLOAT> weightSum0 = Image<RFLOAT>::zeros(w, fc);
    Image<RFLOAT> weightSum1 = Image<RFLOAT>::zeros(w, fc);

    table  = Image<RFLOAT>(w, fc);
    weight = Image<RFLOAT>(w, fc);

    for (int m = 0; m < mgc; m++)
    for (int f = 0; f < fc;  f++)
    for (int x = 0; x < w;   x++) {
        direct::elem(tableSum.data,   x, f) += direct::elem(tables[m].data,   x, f);
        direct::elem(weightSum0.data, x, f) += direct::elem(weights0[m].data, x, f);
        direct::elem(weightSum1.data, x, f) += direct::elem(weights1[m].data, x, f);
    }

    for (int f = 0; f < fc; f++)
    for (int x = 0; x < w; x++) {
        RFLOAT w1 = direct::elem(weightSum0.data, x, f);
        RFLOAT w2 = direct::elem(weightSum1.data, x, f);
        RFLOAT ww = sqrt(w1 * w2);

        direct::elem(weight.data, x, f) = ww;
        direct::elem(table.data,  x, f) = safedivide(direct::elem(tableSum.data, x, f), ww);
    }
}

double FscHelper::computeTsc(
    const std::vector<Image<RFLOAT>> &tables,
    const std::vector<Image<RFLOAT>> &weights0,
    const std::vector<Image<RFLOAT>> &weights1,
    int k0, int k1
) {
    //const int kc = tables[0].data.xdim;
    const int fc = tables[0].data.ydim;
    const int tc = tables.size();

    double tSum = 0.0, w0Sum = 0.0, w1Sum = 0.0;

    for (int t = 0;  t < tc; t++)
    for (int f = 0;  f < fc; f++)
    for (int k = k0; k < k1; k++) {
        tSum  += tables  [t](f, k);
        w0Sum += weights0[t](f, k);
        w1Sum += weights1[t](f, k);
    }

    double ww = sqrt(w0Sum * w1Sum);

    return safedivide(tSum, ww);
}

void FscHelper::computeNoiseSq(
    std::vector<std::vector<Image<Complex> > > frames,
    std::vector<Image<Complex> > predictions, Image<RFLOAT> &sigma2
) {
    const int w = predictions[0].data.xdim;
    const int pc = frames.size();
    const int fc = frames[0].size();

    sigma2              = Image<RFLOAT>::zeros(w, fc);
    Image<RFLOAT> count = Image<RFLOAT>::zeros(w, fc);

    for (int p = 0; p < pc; p++)
    for (int f = 0; f < fc; f++) {
        FOR_ALL_ELEMENTS_IN_FFTW_TRANSFORM(frames[p][f]()) {
            int idx = round(sqrt(kp * kp + ip * ip + jp * jp));

            if (idx >= w) continue;

            Complex z1 = direct::elem(frames[p][f](), i, j, k);
            Complex z2 = direct::elem(predictions[p](), i, j, k);

            direct::elem(sigma2.data, idx, f) += (z2 - z1).norm();
            direct::elem(count.data,  idx, f) += 1.0;
        }
    }

    for (int f = 0; f < fc; f++)
    for (int x = 0; x < w; x++) {
        if (direct::elem(count.data,  x, f) > 0.0) {
            direct::elem(sigma2.data, x, f) /= direct::elem(count.data, x, f);
        }
    }
}

Image<RFLOAT> FscHelper::computeSignalSq(
    const Image<RFLOAT> &sigma2, const Image<RFLOAT> &frc
) {
    const int kc = sigma2.data.xdim;
    const int fc = sigma2.data.ydim;

    Image<RFLOAT> out(kc, fc);

    const RFLOAT eps = 1e-3;

    for (int f = 0; f < fc; f++)
    for (int k = 0; k < kc; k++) {
        RFLOAT c  = direct::elem(frc.data,    k, f);
        RFLOAT s2 = direct::elem(sigma2.data, k, f);

        if (c < eps) {
            direct::elem(out.data, k, f) = 0.0;
        } else {
            if (c > 1.0 - eps) { c = 1.0 - eps; }
            RFLOAT snr2 = c / (1.0 - c);
            direct::elem(out.data, k, f) = snr2 * s2;
        }
    }

    return out;
}

std::vector<d2Vector> FscHelper::fitBfactorsNM(
    const Image<RFLOAT> &tau2, const Image<RFLOAT> &weight, int cutoff
) {
    const int kc = tau2.data.xdim;
    const int ic = tau2.data.ydim;
    std::vector<d2Vector> out(ic);

    const double Bw = 0.001, Cw = 10.0;

    Image<RFLOAT> tauRel(kc, ic);
    std::vector<RFLOAT> tauAvg(kc, 0.0);

    for (int k = 0; k < kc; k++) {
        double ta = 0.0;
        double ts = 0.0;

        for (int f = 0; f < ic; f++) {
            double t2 = direct::elem(tau2.data, k, f);

            if (t2 >= 0.0) {
                ta += t2;
                ts += 1.0;
            }
        }

        if (ts > 0.0) {
            tauAvg[k] = ta / ts;
        }

        for (int f = 0; f < ic; f++) {
            double t2 = direct::elem(tau2.data, k, f);
            direct::elem(tauRel.data, k, f) = t2 >= 0.0 ? sqrt(t2) / tauAvg[k] : 0.0;
        }
    }

    for (int f = 0; f < ic; f++) {
        std::cout << f << ": ";
        BFactorFit bf(tauRel, weight, f, cutoff, Bw, Cw);

        std::vector<double> initial(2);
       // initial[0] = -0.001;
        //initial[1] = -10.0;
        initial[0] = -0.001 / Bw;
        initial[1] = -10.0 / Cw;

        std::vector<double> opt = NelderMead::optimize(initial, bf, 0.01, 0.000001, 1000000);

        out[f][0] = opt[0] * Bw;
        out[f][1] = opt[1] * Cw;

        std::cout << out[f][0] << ", " << out[f][1]  << "\n";
    }

    return out;
}

std::vector<d2Vector> FscHelper::fitBfactors(
    const Image<RFLOAT> &table, const Image<RFLOAT> &weight
) {
    const int kc = table.data.xdim;
    const int ic = table.data.ydim;

    std::vector<d2Vector> out(ic);
    std::vector<RFLOAT> avgFsc(kc, 0.0);

    for (int k = 0; k < kc; k++) {
        RFLOAT wsum = 0.0;

        for (int i = 0; i < ic; i++) {
            if (direct::elem(table.data, i, k) < 0.0) continue;
            avgFsc[k] += direct::elem(weight.data, i, k) * direct::elem(table.data, i, k);
            wsum += direct::elem(weight.data, i, k);
        }

        avgFsc[k] /= wsum;
    }

    const double eps = 1e-3;

    for (int i = 0; i < ic; i++) {
        d2Matrix A(0.0, 0.0, 0.0, 0.0);
        d2Vector b(0.0, 0.0);

        for (int k = 1; k < kc; k++) {
            RFLOAT fsc = direct::elem(table.data, i, k);
            RFLOAT fscA = avgFsc[k];

            if (fsc < eps || fscA < eps) continue;

            RFLOAT t = sqrt((fsc - fsc * fscA) / (fscA - fsc * fscA));

            if (t < eps) continue;

            RFLOAT w = t * t * direct::elem(weight.data, i, k);
            RFLOAT x = k * k;

            A(0, 0) += w * x * x;
            A(1, 0) += w * x;
            A(0, 1) += w * x;
            A(1, 1) += w;

            b[0] += w * log(t) * x;
            b[1] += w * log(t);
        }

        A.invert();
        d2Vector mq = A * b;

        out[i][0] = -4.0 * mq[0];
        out[i][1] = exp(mq[1]);
    }

    return out;
}

Image<RFLOAT> FscHelper::tauRatio(
    const Image<RFLOAT> &table, const Image<RFLOAT> &weight
) {
    const int kc = table.data.xdim;
    const int ic = table.data.ydim;

    Image<RFLOAT> out(kc,ic);
    std::vector<RFLOAT> avgFsc(kc, 0.0);

    for (int k = 0; k < kc; k++) {
        RFLOAT wsum = 0.0;

        for (int i = 0; i < ic; i++) {
            if (direct::elem(table.data, i, k) < 0.0) continue;

            avgFsc[k] += direct::elem(weight.data, i, k) * direct::elem(table.data, i, k);
            wsum += direct::elem(weight.data, i, k);
        }

        avgFsc[k] /= wsum;
    }

    const double eps = 1e-3;

    for (int i = 0; i < ic; i++) {
        for (int k = 0; k < kc; k++) {
            RFLOAT fsc = direct::elem(table.data, i, k);
            RFLOAT fscA = avgFsc[k];

            if (fsc < eps || fscA < eps) continue;

            direct::elem(out.data, i, k) = (fsc - fsc * fscA) / (fscA - fsc * fscA);
        }
    }

    return out;
}

void FscHelper::computeBfactors(
    const std::vector<d2Vector>& bfacs, Image<RFLOAT> &table
) {
    const int kc = table.data.xdim;
    const int ic = table.data.ydim;

    for (int i = 0; i < ic; i++)
    for (int k = 0; k < kc; k++) {
        const double Bf = bfacs[i][0];
        const double Cf = bfacs[i][1];

        direct::elem(table.data, i, k) = exp(Bf * k * k / 4.0 + Cf);
    }
}

std::vector<double> FscHelper::powerSpectrum3D(const Image<Complex> &img) {
    const int sh = img.data.xdim;
    const int s = img.data.ydim;

    std::vector<double> sum(sh, 0.0), wgh(sh, 0.0);

    for (int z = 0; z < img.data.zdim; z++)
    for (int y = 0; y < img.data.ydim; y++)
    for (int x = 0; x < img.data.xdim; x++) {
        const double xx = x;
        const double yy = y < sh ? y : y - s;
        const double zz = z < sh ? z : z - s;

        const double r = sqrt(xx * xx + yy * yy + zz * zz);
        const int ri = round(r);

        if (ri < sh) {
            sum[ri] += direct::elem(img.data, x, y, z).norm();
            wgh[ri] += 1.0;
        }
    }

    for (int i = 0; i < sh; i++) {
        if (wgh[i] > 0.0) {
            sum[i] /= wgh[i];
        }
    }

    return sum;
}

BFactorFit::BFactorFit(
    const Image<RFLOAT> &tau2, const Image<RFLOAT> &weight,
    int frame, int cutoff, double Bscale, double Cscale
): tau2(tau2), weight(weight), frame(frame), cutoff(cutoff), Bscale(Bscale), Cscale(Cscale) {}

double BFactorFit::f(const std::vector<double> &x, void *tempStorage) const {
    double Bf = x[0] * Bscale;
    double Cf = x[1] * Cscale;

    double sum = 0.0;

    const int kc = tau2.data.xdim;

    for (int k = cutoff; k < kc; k++) {
        double w = direct::elem(weight.data, frame, k);
        RFLOAT t2 = direct::elem(tau2.data, frame, k);
        double pv = exp(Bf * k * k / 4.0 + Cf);

        double e = pv - t2;

        sum += w * e * e;
    }

    return sum;
}

