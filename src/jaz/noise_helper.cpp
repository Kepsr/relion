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

#include "src/jaz/noise_helper.h"
#include "src/jaz/vtk_helper.h"
#include "src/jaz/distribution_helper.h"
#include "src/jaz/img_proc/filter_helper.h"
#include <omp.h>


Image<RFLOAT> NoiseHelper::predictCCNoise(
    Projector &prj, double sigma2,
    double nsamples_ppp, int max_nsamples, int nangles, Image<RFLOAT> &dmgWeight,
    CTF ctf0, ObservationModel *obsModel, int opticsGroup,
    double defocusMu, double defocusSigma, double angpix, int thread_num
) {
    const int s = prj.ori_size;
    const int sh = s / 2 + 1;

    const int vbins = 1024;
    const bool binValues = true;

    int goodAngles = 0;

    Image<RFLOAT> confusion = Image<RFLOAT>::zeros(s, s);

    while (goodAngles < nangles) {
        double dx = 2.0 * rand() / (double) RAND_MAX - 1.0;
        double dy = 2.0 * rand() / (double) RAND_MAX - 1.0;
        double dz = 2.0 * rand() / (double) RAND_MAX - 1.0;

        double dd = dx * dx + dy * dy + dz * dz;

        if (dd > 1.0) continue;

        goodAngles++;

        if (goodAngles % 10 == 0) std::cout << goodAngles << "/" << nangles << "\n";

        Vector<RFLOAT> dm = vectorR3(dx, dy, dz);

        RFLOAT rot, tilt;
        Euler::direction2angles(dm, rot, tilt);

        Matrix<RFLOAT> A3D = Euler::angles2matrix(rot, tilt, 0.0);

        Image<Complex> spec (prj.get2DFourierTransform(sh, s, 1, A3D));

        double defocus = DistributionHelper::sampleGauss(defocusMu, defocusSigma);

        ctf0.DeltafU = defocus;
        ctf0.DeltafV = defocus;
        ctf0.initialise();

        FilterHelper::modulate(spec(), ctf0, obsModel, opticsGroup, angpix);

        Image<Complex> ccspec (s, sh);
        for (long int j = 0; j < s; j++)
        for (long int i = 0; i < sh; i++) {
            direct::elem(  spec.data, i, j) *= sqrt(direct::elem(dmgWeight.data, i, j));
            direct::elem(ccspec.data, i, j) = direct::elem(spec.data, i, j).norm();
        }

        FourierTransformer ft;
        Image<RFLOAT> mu0 (ft.inverseFourierTransform(ccspec.data));
        Image<RFLOAT> img (ft.inverseFourierTransform(spec.data));

        const Image<RFLOAT> mu = mu0;

        double varScale = 0.0;
        for (auto const& x: img.data) {
            varScale += x * x;
        }
        varScale /= s * s;

        const double sigma2CC = sigma2 * varScale;
        const double sigmaCC = sqrt(sigma2CC);

        const auto minmax = std::minmax_element(mu.data.begin(), mu.data.end());
        const RFLOAT vMin = *minmax.first;
        const RFLOAT vMax = *minmax.second;

        std::vector<double> plausibleVals;
        std::vector<std::pair<int, int>> plausiblePixels;
        plausibleVals.reserve(s * sh);
        plausiblePixels.reserve(s * sh);

        for (long int j = 0; j < sh; j++)
        for (long int i = 0; i < s;  i++) {
            const double x = direct::elem(mu.data, i, j);
            if (vMax - x <= 6.0 * sigmaCC) {
                plausibleVals.push_back(x);
                plausiblePixels.emplace_back(i, j);
            }
        }

        const int npp = plausibleVals.size();
        if (npp == 1) {
            std::pair<int, int> origin = plausiblePixels[0];
            direct::elem(confusion.data, origin.second, origin.first) += 1.0;
            continue;
        }

        const int nsamples = std::min((int) (npp * nsamples_ppp), max_nsamples);

        std::vector<int> hitsPerBin (vbins, 0), ppixelsPerBin (vbins, 0);

        const double floorBin = vMax - vMin < 6.0 * sigmaCC ? vMin : vMax - 6.0 * sigmaCC;
        const double binRange = vMax - floorBin;

        if (binValues) {
            for (double const& x: mu.data) {
                if (vMax - x <= 6.0 * sigmaCC) {
                    const int i = (vbins - 1) * (x - floorBin) / binRange;
                    ppixelsPerBin[i]++;
                }
            }
        }

        const int threads = thread_num;

        const double max_double = std::numeric_limits<RFLOAT>::max();

        if (threads > 1) {

            const int stride = 2048;
            std::vector<unsigned int> randThreadStates (threads * stride);
            for (int i = 0; i < threads * stride; i += stride) {
                randThreadStates[i] = rand();
            }

            #pragma omp parallel num_threads(threads)
            {
                int threadnum = omp_get_thread_num();

                for (int i = 0; i < nsamples / threads + 1; i++) {
                    double vmax = -max_double;
                    int jmax = 0;

                    for (long int j = 0; j < npp; j++) {

                        const double m = plausibleVals[j];
                        const double u1 = rand_r(&randThreadStates[stride * threadnum]) / (double) RAND_MAX;
                        const double u2 = rand_r(&randThreadStates[stride * threadnum]) / (double) RAND_MAX;
                        const double z = sqrt(-2.0 * log(u1)) * cos(2.0 * PI * u2);
                        const double v = m + sigmaCC * z;

                        if (v > vmax) {
                            vmax = v;
                            jmax = j;
                        }
                    }

                    if (binValues) {
                        const int b = (vbins - 1) * (plausibleVals[jmax] - floorBin) / binRange;

                        #pragma omp atomic
                        hitsPerBin[b]++;
                    } else {
                        const int i = plausiblePixels[jmax].first;
                        const int j = plausiblePixels[jmax].second;
                        #pragma omp atomic
                        direct::elem(confusion.data, i, j) += 1.0;
                    }
                }
            }
        } else {
            for (int i = 0; i < nsamples; i++) {
                double vmax = -max_double;
                int jmax = 0;

                for (long int j = 0; j < npp; j++) {

                    const double m = plausibleVals[j];
                    const double u1 = rand() / (double) RAND_MAX;
                    const double u2 = rand() / (double) RAND_MAX;
                    const double z = sqrt(-2.0 * log(u1)) * cos(2.0 * PI * u2);
                    const double v =  m + sigmaCC * z;

                    if (v > vmax) {
                        vmax = v;
                        jmax = j;
                    }
                }

                if (binValues) {
                    const int b = (vbins - 1) * (plausibleVals[jmax] - floorBin) / binRange;
                    hitsPerBin[b]++;
                } else {
                    const int i = plausiblePixels[jmax].first;
                    const int j = plausiblePixels[jmax].second;
                    direct::elem(confusion.data, i, j) += 1.0;
                }
            }
        }

        if (binValues) {
            for (long int j = 0; j < sh; j++)
            for (long int i = 0; i < s;  i++) {

                const double x = direct::elem(mu.data, i, j);

                if (x < floorBin) continue;

                const int b = (vbins - 1) * (x - floorBin) / binRange;

                if (b >= 0 && hitsPerBin[b] > 0) {
                    direct::elem(confusion.data, i, j) += hitsPerBin[b] / (double) ppixelsPerBin[b];
                }
            }
        }
    }

    return confusion;
}

std::vector<double> NoiseHelper::radialAverage(Image<RFLOAT> &map, bool half) {
    const int w = map.data.xdim;
    const int h = map.data.ydim;
    const int b = (w + h) / 4;

    std::vector<double> out (b, 0.0);
    std::vector<double> wgh (b, 0.0);

    const int ha = half ? h / 2 : h;

    for (int j = 0; j < ha; j++)
    for (int i = 0; i < w;  i++) {

        const auto v = direct::elem(map.data, i, j);
        const double x = i < w / 2.0 ? i : i - w;
        const double y = j < h / 2.0 ? j : j - h;
        const double r = std::hypot(x, y);
        const int ri = r;
        const double rf = r - ri;

        if (ri < b) {
            out[ri] += (1.0 - rf) * v;
            wgh[ri] += (1.0 - rf);
        }

        if (ri < b - 1) {
            out[ri + 1] += rf * v;
            wgh[ri + 1] += rf;
        }
    }

    for (int i = 0; i < b; i++)
        if (wgh[i] > 0.0)
            out[i] /= wgh[i];

    return out;
}

Image<RFLOAT> NoiseHelper::radialMap(std::vector<double> &radAvg, bool centered) {
    const int b = radAvg.size();
    const int s = 2 * b;

    MultidimArray<RFLOAT> out (s, s);

    for (int j = 0; j < s; j++)
    for (int i = 0; i < s; i++) {

        double x, y;
        if (centered) {
            x = i - s / 2;
            y = j - s / 2;
        } else {
            x = i < s / 2 ? i : i - s;
            y = j < s / 2 ? j : j - s;
        }

        const double r = std::hypot(x, y);
        const int ri = r;          // Integral part
        const double rf = r - ri;  // Fractional part

        direct::elem(out, i, j) = ri < b - 1 ?
            rf * radAvg[ri + 1] + (1 - rf) * radAvg[ri] :
            radAvg[b - 1];
    }

    return {out};
}

std::vector<Complex> NoiseHelper::radialAverage(Image<Complex> &map, bool skipAxes) {
    const int w = map.data.xdim;
    const int h = map.data.ydim;
    const int b = w;

    std::vector<Complex> out (b, 0);
    std::vector<double>  wgh (b, 0);

    for (int j = 0; j < h; j++)
    for (int i = 0; i < w; i++) {

        if (skipAxes && (i == 0 || j == 0)) continue;

        const double x = i, y = j < h / 2.0 ? j : j - h;
        const double r = std::hypot(x, y);
        const int ri = r;
        const double rf = r - ri;

        if (ri < b) {
            out[ri] += (1 - rf) * direct::elem(map.data, i, j);
            wgh[ri] += (1 - rf);
        }

        if (ri < b - 1) {
            out[ri + 1] += rf * direct::elem(map.data, i, j);
            wgh[ri + 1] += rf;
        }
    }

    for (int i = 0; i < b; i++)
        if (wgh[i] > 0)
            out[i] /= wgh[i];

    return out;
}

Image<Complex> NoiseHelper::radialMap(std::vector<Complex> &radAvg) {
    const int b = radAvg.size();
    const int s = 2 * b - 2;

    MultidimArray<Complex> out (b, s);

    for (int j = 0; j < s; j++)
    for (int i = 0; i < b; i++) {

        const double x = i, y = j < s / 2 ? j : j - s;

        const double r = std::hypot(x, y);
        const int ri = r;
        const double rf = r - ri;

        direct::elem(out, i, j) = ri < b - 1 ?
            rf * radAvg[ri + 1] + (1 - rf) * radAvg[ri] :
            radAvg[b - 1];
    }

	return {out};
}

std::vector<std::pair<double,double>> NoiseHelper::radialAverageAndStdDevFFTW(Image<RFLOAT> &map) {
    const int w = map.data.xdim;
    const int h = map.data.ydim;
    const int b = w;

    std::vector<double> avg (b, 0.0), wgh (b, 0.0), var (b, 0.0);

    for (int j = 0; j < h; j++)
    for (int i = 0; i < w; i++) {
        const double x = i, y = j < h / 2.0 ? j : j - h;
        const int r = std::hypot(x, y) + 0.5;

        if (r < b) {
            avg[r] += direct::elem(map.data, i, j);
            wgh[r] += 1.0;
        }
    }

    for (int i = 0; i < b; i++)
        if (wgh[i] > 0.0)
            avg[i] /= wgh[i];

	for (int j = 0; j < h; j++)
    for (int i = 0; i < w; i++) {

        const double x = i, y = j < h / 2.0 ? j : j - h;
		const int r = std::hypot(x, y) + 0.5;
		const double mu = avg[r];
		const double v = direct::elem(map.data, i, j) - mu;

        if (r < b) { var[r] += v * v; }
    }

	for (int i = 0; i < b; i++)
        if (wgh[i] > 1.0)
            var[i] /= wgh[i] - 1;

	std::vector<std::pair<double, double>> out;
    out.reserve(b);
	for (int i = 0; i < b; i++) {
		out.emplace_back(avg[i], sqrt(var[i]));
	}
    return out;
}

std::vector<double> NoiseHelper::radialWeight(int w, int h, bool half) {
    const int b = (w + h) / 4;
    const int ha = half ? h / 2 : h;

    std::vector<double> wgh (b, 0.0);

    for (int j = 0; j < ha; j++)
    for (int i = 0; i < w;  i++) {
        const double x = i < w / 2.0 ? i : i - w;
        const double y = j < h / 2.0 ? j : j - h;
        const int r = std::hypot(x, y) + 0.5;

        if (r < b)
            wgh[r] += 1.0;
    }

    return wgh;
}

std::vector<double> NoiseHelper::fill(Image<RFLOAT>& confusion, double lambda, int iterations) {
    const int s0 = confusion.data.xdim;

    const std::vector<double> ra = radialAverage(confusion, true);
    const std::vector<double> rw = radialWeight(s0, s0, true);

    const int s = ra.size();

    std::vector<double> w;
    w.reserve(s);
    for (int i = 0; i < s; i++)
        w.push_back(rw[i] * rw[i] * ra[i] * ra[i] + 1.0);

    std::vector<double> v = ra, vn = ra;

    for (int it = 0; it < iterations; it++) {
        for (int i = 1; i < s - 1; i++) {
            vn[i] = (lambda * 0.5 * (v[i + 1] + v[i - 1]) + w[i] * ra[i])
                  / (lambda + w[i]);
        }

        vn[s - 1] = (lambda * 0.5 * v[s - 2] + w[s - 1] * ra[s - 1])
                  / (lambda * 0.5 + w[s - 1]);

        std::copy_n(vn.begin(), s, v.begin());
    }

    return v;
}

Image<RFLOAT> NoiseHelper::normalize(const Image<RFLOAT> &confusion) {
    const int w = confusion.data.xdim, h = confusion.data.ydim;
    const double sum = std::accumulate(confusion.data.begin(),
                                       confusion.data.end(), 0.0);

    if (sum <= 0.0)
        return {MultidimArray<RFLOAT>::zeros(w, h)};
    else
        return {confusion.data / sum};
}

void NoiseHelper::testVariance(Image<RFLOAT> img) {
    const int s = img.data.xdim;
    const int sh = s / 2 + 1;

    Image<Complex> spec (MultidimArray<Complex>::zeros(s, sh));

    FourierTransformer ft;
    spec() = ft.FourierTransform(img());

    Image<Complex> ccspec(s, sh);
    for (long int yy = 0; yy < s;  yy++)
    for (long int xx = 0; xx < sh; xx++) {
        direct::elem(ccspec.data, xx, yy) = direct::elem(spec.data, xx, yy).norm();
    }

    Image<RFLOAT> mu (ft.inverseFourierTransform(ccspec.data));

    const double varScale = std::accumulate(img.data.begin(), img.data.end(), 0.0,
        [] (double const& acc, double const& x) { return acc + x * x; }
    ) / (s * s);

    Image<RFLOAT> imgD = img;
    const double sig2 = 2.0;
    const double sig = sqrt(sig2);

    Image<RFLOAT> varImg = img;
    std::fill_n(varImg.data.begin(), varImg.data.end(), 0.0);

    Image<Complex> imgDs(s, sh), ccDs(s, sh);

    const int N = 10000;
    for (int i = 0; i < N; i++) {
        if (i % 10 == 0) std::cout << i << "\n";

        for (long int i = 0; i < imgD.data.size(); i++)
            imgD.data[i] = DistributionHelper::sampleGauss(img.data[i], sig);

        imgDs() = ft.FourierTransform(imgD());

        for (long int i = 0; i < ccDs.data.size(); i++)
            ccDs.data[i] = spec.data[i] * imgDs.data[i].conj();

        MultidimArray<RFLOAT> ccD = ft.inverseFourierTransform(ccDs.data);

        // varImg += (ccD - mu) * (ccD - mu) / N;
        for (long int i = 0; i < s; i++) {
            const double d = ccD[i] - mu.data[i];
            varImg.data[i] += d * d / N;
        }

    }

    const double varSum = std::accumulate(varImg.data.begin(),
                                          varImg.data.end(), 0.0) / (s * s);

    std::cout << varSum << " vs. " << sig2 * varScale << "\n";
}

void NoiseHelper::testColorVariance(Image<RFLOAT> img, std::vector<double> sig2) {
    const int s = img.data.xdim;
    const int sh = s / 2 + 1;

    FourierTransformer ft;
    Image<Complex> spec (ft.FourierTransform(img()));

    VtkHelper::writeVTK(spec, "debug/spec.vtk");

    for (long int j = 0; j < s;  j++)
    for (long int i = 0; i < sh; i++) {
        if (i == 0 && j == 0) continue;

        const double y = j < sh ? j : j - s;
        const double x = i;
        const int r = std::hypot(x, y);

        if (r >= sh) direct::elem(spec.data, i, j) = Complex(0.0, 0.0);
    }


    double varPred = 0.0;

    for (long int j = 0; j < s;  j++)
    for (long int i = 0; i < sh; i++) {
        if (i == 0 && j == 0) continue;

        const double y = j < sh ? j : j - s;
        const double x = i;
        const int r = std::hypot(x, y);

        if (r < sh && r > 0) {
            const double a = direct::elem(spec.data, i, j).norm() / sig2[r];
            varPred += i == 0 ? a : 2.0 * a;
        }
    }

    // varPred *= s * s;

    std::cout << "pred: " << varPred << "\n";

    std::vector<double> sig;
    sig.reserve(sh);
    for (double const& x: sig2)
        sig.push_back(sqrt(x));

    MultidimArray<RFLOAT> varImg (img.data.xdim, img.data.ydim);
    std::fill(varImg.begin(), varImg.end(), 0.0);

    Image<Complex> ccDs (s, sh);

    const int N = 10000;
    const double sqrt_half = sqrt(0.5);

    double varTest = 0.0;

    for (int i = 0; i < N; i++) {
        if (i % 100 == 0) std::cout << i << "\n";

        for (long int j = 0; j < s;  j++)
        for (long int i = 0; i < sh; i++) {

            const double y = j < sh ? j : j - s;
            const double x = i;
            const int r = std::hypot(x, y);

            if (r < sh && r > 0) {

                const double r0 = DistributionHelper::sampleGauss(0, sqrt_half);
                const double r1 = DistributionHelper::sampleGauss(0, sqrt_half);
                Complex z0 = direct::elem(spec.data, i, j);
                Complex z1 = Complex(r0, r1) * z0 / sig[r];

                direct::elem(ccDs.data, i, j) = z1;

                if (i == 0 && j >= sh) {
                    direct::elem(ccDs.data, i, j) = direct::elem(ccDs.data, i, s - j).conj();
                }

                varTest += i == 0 ? direct::elem(ccDs.data, i, j).norm() :
                              2.0 * direct::elem(ccDs.data, i, j).norm();
            } else {
                direct::elem(ccDs.data, i, j) = 0.0;
            }
        }

        MultidimArray<RFLOAT> ccD = ft.inverseFourierTransform(ccDs.data);

        for (long int i = 0; i < s; i++) {
            const double d = ccD[i];
            varImg[i] += d * d / (N * varPred);
        }
    }

    varTest /= N;

    VtkHelper::writeVTK({varImg}, "debug/varImg_cn.vtk");

    const double varSum = std::accumulate(varImg.begin(),
                                          varImg.end(), 0.0) / (s * s);

    std::cout << varSum << " @ " << varPred << " vs. " << varTest << "\n";

}

void NoiseHelper::testParseval() {
    const int s = 512;
    const int sh = s / 2 + 1;

    MultidimArray<RFLOAT>  real (s, s);
    MultidimArray<Complex> freq (s, sh);

    for (int i = 0; i < 10; i++) {

        double varr = 0.0;
        for (RFLOAT& x: real) {
            x = DistributionHelper::sampleGauss(0, 2);
            varr += x * x;
        }
        varr /= s * s;

        FourierTransformer ft;
        freq = ft.FourierTransform(real);

        double var = 0.0;
        for (int j = 0; j < s;  j++)
        for (int i = 0; i < sh; i++) {
            const Complex z = direct::elem(freq, i, j);
            var += i == 0 ? z.norm() : 2 * z.norm();
        }
        var /= s * s;

        std::cout << varr << " vs. " << var << " (" << var * (double) (s * s) << ")\n";
        // varr = varf * A
    }

    std::cout << "\n";

    const double sqrt_2 = sqrt(2.0);
    for (int iter = 0; iter < 10; iter++) {

        double varf = 0.0;
        for (int j = 0; j < s;  j++)
        for (int i = 0; i < sh; i++) {
            direct::elem(freq, i, j).real = DistributionHelper::sampleGauss(0, sqrt_2);
            direct::elem(freq, i, j).imag = i > 0 && i < sh - 1 ?
                DistributionHelper::sampleGauss(0, sqrt_2) : 0.0;

            varf += i == 0 ? direct::elem(freq, i, j).norm() :
                       2.0 * direct::elem(freq, i, j).norm();
        }
        varf /= s * s;

        FourierTransformer ft;
        real = ft.inverseFourierTransform(freq);

        double var = 0.0;
        for (RFLOAT const& x: real)
            var += x * x;
        var /= s * s;

        std::cout << varf << " vs. " << var << " (" << var / (double) (s * s) << ")\n";
        // varf = varr / A
    }
}
