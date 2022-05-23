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
    CTF ctf0, double defocusMu, double defocusSigma, double angpix, int thread_num
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

        Matrix1D<RFLOAT> dm = VECTOR_R3(dx, dy, dz);

        RFLOAT rot, tilt;
        Euler::direction2angles(dm, rot, tilt);

        Matrix2D<RFLOAT> A3D = Euler::angles2matrix(rot, tilt, 0.0);

        Image<Complex> spec = Image<Complex>::zeros(sh, s);

        prj.get2DFourierTransform(spec.data, A3D);

        double defocus = DistributionHelper::sampleGauss(defocusMu, defocusSigma);

        ctf0.DeltafU = defocus;
        ctf0.DeltafV = defocus;
        ctf0.initialise();

        FilterHelper::modulate(spec(), ctf0, angpix);

        Image<Complex> ccspec(sh, s);
        for (long int yy = 0; yy < s; yy++)
        for (long int xx = 0; xx < sh; xx++) {
            direct::elem(  spec.data, xx, yy) *= sqrt(direct::elem(dmgWeight.data, xx, yy));
            direct::elem(ccspec.data, xx, yy) = direct::elem(spec.data, xx, yy).norm();
        }

        Image<RFLOAT> mu0(s,s), img(s,s);

        FourierTransformer ft;
        ft.inverseFourierTransform(ccspec.data, mu0.data);
        ft.inverseFourierTransform(spec.data, img.data);

        const Image<RFLOAT> mu = mu0;

        double varScale = 0.0;

        for (long int yy = 0; yy < s; yy++)
        for (long int xx = 0; xx < s; xx++) {
            double m = direct::elem(img.data, xx, yy) / (s * s);
            varScale += m * m;
        }

        const double sigma2CC = sigma2*varScale;
        const double sigmaCC = sqrt(sigma2CC);

        RFLOAT vMin = +std::numeric_limits<RFLOAT>::max();
        RFLOAT vMax = -std::numeric_limits<RFLOAT>::max();

        for (long int y = 0; y < s; y++)
        for (long int x = 0; x < s; x++) {
            double v = direct::elem(mu.data, x, y);
            if (v > vMax) { vMax = v; }
            if (v < vMin) { vMin = v; }
        }


        std::vector<double> plausibleVals(0);
        std::vector<std::pair<int,int>> plausiblePixels(0);
        plausibleVals.reserve(s*sh);
        plausiblePixels.reserve(s*sh);

        for (long int y = 0; y < sh; y++)
        for (long int x = 0; x < s;  x++) {
            double m = direct::elem(mu.data, x, y);
            double dm = vMax - m;

            if (dm <= 6.0 * sigmaCC) {
                plausibleVals.push_back(m);
                plausiblePixels.push_back(std::make_pair(x,y));
            }
        }

        const int npp = plausibleVals.size();

        if (npp == 1) {
            std::pair<int, int> origin = plausiblePixels[0];
            direct::elem(confusion.data, origin.second, origin.first) += 1.0;

            continue;
        }

        const int nsamples = std::min((int)(npp*nsamples_ppp), max_nsamples);

        std::vector<int> hitsPerBin(vbins, 0), ppixelsPerBin(vbins, 0);

        const double floorBin = vMax - vMin < 6.0 * sigmaCC ? vMin : vMax - 6.0 * sigmaCC;
        const double binRange = vMax - floorBin;

        if (binValues) {
            for (long int y = 0; y < sh; y++)
            for (long int x = 0; x < s;  x++) {
                double m = direct::elem(mu.data, x, y);
                double dm = vMax - m;

                if (dm <= 6.0 * sigmaCC) {
                    const int b = (vbins - 1) * (m - floorBin) / binRange;
                    ppixelsPerBin[b]++;
                }
            }
        }

        const int threads = thread_num;

        const double max_double = std::numeric_limits<RFLOAT>::max();

        if (threads > 1) {
            const int stride = 2048;
            std::vector<unsigned int> randThreadStates(threads * stride);

            for (int i = 0; i < threads; i++) {
                randThreadStates[stride * i] = rand();
            }

            #pragma omp parallel num_threads(threads)
            {
                int threadnum = omp_get_thread_num();

                for (int i = 0; i < nsamples / threads + 1; i++) {
                    double vmax = -max_double;
                    int jmax = 0;

                    for (long int j = 0; j < npp; j++) {
                        double m = plausibleVals[j];

                        double u1 = rand_r(&randThreadStates[stride * threadnum]) / (double) RAND_MAX;
                        double u2 = rand_r(&randThreadStates[stride * threadnum]) / (double) RAND_MAX;

                        double z = sqrt(-2.0 * log(u1)) * cos(2.0 * PI * u2);

                        double v = m + sigmaCC * z;

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
                        int x = plausiblePixels[jmax].first;
                        int y = plausiblePixels[jmax].second;

                        #pragma omp atomic
                        direct::elem(confusion.data, x, y) += 1.0;
                    }
                }
            }
        } else {
            for (int i = 0; i < nsamples; i++) {
                double vmax = -max_double;
                int jmax = 0;

                for (long int j = 0; j < npp; j++) {
                    double m = plausibleVals[j];

                    double u1 = rand() / (double) RAND_MAX;
                    double u2 = rand() / (double) RAND_MAX;

                    double z = sqrt(-2.0 * log(u1)) * cos(2.0 * PI * u2);

                    double v =  m + sigmaCC * z;

                    if (v > vmax) {
                        vmax = v;
                        jmax = j;
                    }
                }

                if (binValues) {
                    const int b = (vbins - 1) * (plausibleVals[jmax] - floorBin) / binRange;
                    hitsPerBin[b]++;
                } else {
                    int x = plausiblePixels[jmax].first;
                    int y = plausiblePixels[jmax].second;

                    direct::elem(confusion.data, x, y) += 1.0;
                }
            }
        }

        if (binValues) {
            for (long int y = 0; y < sh; y++)
            for (long int x = 0; x < s;  x++) {
                const double m = direct::elem(mu.data, x, y);

                if (m < floorBin) continue;

                const int b = (vbins - 1) * (m - floorBin) / binRange;

                if (b >= 0 && hitsPerBin[b] > 0) {
                    direct::elem(confusion.data, x, y) += hitsPerBin[b] / (double) ppixelsPerBin[b];
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

    std::vector<double> out(b, 0.0);
    std::vector<double> wgh(b, 0.0);

    const int ha = half ? h / 2 : h;

    for (int yy = 0; yy < ha; yy++)
    for (int xx = 0; xx < w;  xx++) {
        double x = xx < w / 2.0 ? xx : xx - w;
        double y = yy < h / 2.0 ? yy : yy - h;
        double rd = sqrt(x * x + y * y);

        int r0 = rd;
        double a = rd - r0;

        if (r0 < b) {
            out[r0] += (1.0 - a) * direct::elem(map.data, xx, yy);
            wgh[r0] += 1.0 - a;
        }

        if (r0 < b - 1) {
            out[r0 + 1] += a * direct::elem(map.data, xx, yy);
            wgh[r0 + 1] += a;
        }
    }

    for (int i = 0; i < b; i++) { if (wgh[i] > 0.0) { out[i] /= wgh[i]; } }

    return out;
}

Image<RFLOAT> NoiseHelper::radialMap(std::vector<double> &radAvg, bool centered) {
    const int b = radAvg.size();
    const int s = 2 * b;

    Image<RFLOAT> out(s, s);

    for (int yy = 0; yy < s; yy++)
    for (int xx = 0; xx < s; xx++) {
        double x, y;

        if (centered) {
            x = xx - s / 2;
            y = yy - s / 2;
        } else {
            x = xx < s / 2 ? xx : xx - s;
            y = yy < s / 2 ? yy : yy - s;
        }

        double rd = sqrt(x * x + y * y);

        int r0 = rd;
        double a = rd - r0;

        if (r0 < b - 1) {
            direct::elem(out.data, xx, yy) = a * radAvg[r0 + 1] + (1 - a) * radAvg[r0];
        } else {
            direct::elem(out.data, xx, yy) = radAvg[b - 1];
        }
    }

    return out;
}

std::vector<Complex> NoiseHelper::radialAverage(Image<Complex> &map, bool skipAxes) {
    const int w = map.data.xdim;
    const int h = map.data.ydim;
    const int b = w;

    std::vector<Complex> out(b, 0.0);
    std::vector<double> wgh(b, 0.0);

    for (int yy = 0; yy < h; yy++)
    for (int xx = 0; xx < w; xx++) {
        if (skipAxes && (xx == 0 || yy == 0)) continue;

        double x = xx;
        double y = yy < h / 2.0 ? yy : yy - h;
        double rd = sqrt(x * x + y * y);

        int r0 = rd;
        double a = rd - r0;

        if (r0 < b) {
            out[r0] += (1.0 - a) * direct::elem(map.data, xx, yy);
            wgh[r0] += 1.0 - a;
        }

        if (r0 < b-1) {
            out[r0 + 1] += a * direct::elem(map.data, xx, yy);
            wgh[r0 + 1] += a;
        }
    }

    for (int i = 0; i < b; i++) if (wgh[i] > 0.0) { out[i] /= wgh[i]; }

    return out;
}

Image<Complex> NoiseHelper::radialMap(std::vector<Complex> &radAvg) {
    const int b = radAvg.size();
    const int s = 2*b - 2;
    Image<Complex> out(b,s);

    for (int yy = 0; yy < s; yy++)
    for (int xx = 0; xx < b; xx++) {
        double x = xx;
        double y = yy < s / 2 ? yy : yy - s;

        double rd = sqrt(x * x + y * y);

        int r0 = rd;
        double a = rd - r0;

        if (r0 < b - 1) {
            direct::elem(out.data, xx, yy) = a * radAvg[r0 + 1] + (1 - a) * radAvg[r0];
        } else {
            direct::elem(out.data, xx, yy) = radAvg[b - 1];
        }
    }

	return out;
}

std::vector<std::pair<double,double>> NoiseHelper::radialAverageAndStdDevFFTW(Image<RFLOAT> &map) {
    const int w = map.data.xdim;
    const int h = map.data.ydim;
    const int b = w;

    std::vector<double> avg(b, 0.0);
	std::vector<double> wgh(b, 0.0);
	std::vector<double> var(b, 0.0);

    for (int yy = 0; yy < h; yy++)
    for (int xx = 0; xx < w; xx++) {
        double x = xx;
        double y = yy < h / 2.0 ? yy : yy - h;
        double rd = sqrt(x * x + y * y);

        int r = rd + 0.5;

        if (r < b) {
            avg[r] += direct::elem(map.data, xx, yy);
            wgh[r] += 1.0;
        }
    }

    for (int i = 0; i < b; i++) if (wgh[i] > 0.0) { avg[i] /= wgh[i]; }

	for (int yy = 0; yy < h; yy++)
    for (int xx = 0; xx < w; xx++) {
        double x = xx;
        double y = yy < h / 2.0 ? yy : yy - h;
        double rd = sqrt(x * x + y * y);

		int r = rd + 0.5;

		double mu = avg[r];
		double v = direct::elem(map.data, xx, yy) - mu;

        if (r < b) { var[r] += v * v; }
    }

	for (int i = 0; i < b; i++) if (wgh[i] > 1.0) { var[i] /= wgh[i] - 1; }

	std::vector<std::pair<double,double>> out(b);

	for (int i = 0; i < b; i++) {
		out[i] = std::make_pair(avg[i], sqrt(var[i]));
	}

    return out;
}

std::vector<double> NoiseHelper::radialWeight(int w, int h, bool half) {
    const int b = (w + h) / 4;

    std::vector<double> wgh(b, 0.0);

    const int ha = half ? h / 2 : h;

    for (int yy = 0; yy < ha; yy++)
    for (int xx = 0; xx < w;  xx++) {
        double x = xx < w / 2.0 ? xx : xx - w;
        double y = yy < h / 2.0 ? yy : yy - h;
        double rd = sqrt(x * x + y * y);

        int r = rd + 0.5;

        if (r < b) { wgh[r] += 1.0; }
    }

    return wgh;
}

std::vector<double> NoiseHelper::fill(Image<RFLOAT>& confusion, double lambda, int iterations) {
    const int s0 = confusion.data.xdim;

    std::vector<double> ra = radialAverage(confusion, true);
    std::vector<double> rw = radialWeight(s0, s0, true);

    const int s = ra.size();

    std::vector<double> w(s);

    for (int x = 0; x < s; x++) {
        w[x] = rw[x] * rw[x] * ra[x] * ra[x] + 1.0;
    }

    std::vector<double> v = ra, vn = ra;

    for (int it = 0; it < iterations; it++) {
        for (int x = 1; x < s-1; x++) {
            vn[x] = (lambda * 0.5 * (v[x + 1] + v[x - 1]) + w[x] * ra[x]) 
                  / (lambda + w[x]);
        }

        vn[s - 1] = (lambda * 0.5 * v[s - 2] + w[s - 1] * ra[s - 1]) 
                  / (lambda * 0.5 + w[s - 1]);

        for (int x = 0; x < s; x++) { v[x] = vn[x]; }
    }

    return v;
}

Image<RFLOAT> NoiseHelper::normalize(const Image<RFLOAT> &confusion) {
    const int w = confusion.data.xdim;
    const int h = confusion.data.ydim;

    double sum = 0.0;
    for (int y = 0; y < h; y++)
    for (int x = 0; x < w; x++) {
        sum += direct::elem(confusion.data, x, y);
    }

    if (sum <= 0.0) return Image<RFLOAT>::zeros(w, h);

    Image<RFLOAT> out(w, h);

    for (int y = 0; y < h; y++)
    for (int x = 0; x < w; x++) {
        direct::elem(out.data, x, y) = direct::elem(confusion.data, x, y) / sum;
    }

    return out;
}

void NoiseHelper::testVariance(Image<RFLOAT> img) {
    const int s = img.data.xdim;
    const int sh = s / 2 + 1;

    Image<Complex> spec = Image<Complex>::zeros(sh, s);

    FourierTransformer ft;
    ft.FourierTransform(img(), spec());

    Image<Complex> ccspec(sh, s);
    for (long int yy = 0; yy < s;  yy++)
    for (long int xx = 0; xx < sh; xx++) {
        direct::elem(ccspec.data, xx, yy) = direct::elem(spec.data, xx, yy).norm();
    }

    Image<RFLOAT> mu(s,s);
    ft.inverseFourierTransform(ccspec.data, mu.data);

    double varScale = 0.0;

    for (long int yy = 0; yy < s; yy++)
    for (long int xx = 0; xx < s; xx++) {
        double m = direct::elem(img.data, xx, yy) / (s * s);
        varScale += m * m;
    }

    Image<RFLOAT> imgD = img;
    const double sig2 = 2.0;
    const double sig = sqrt(sig2);

    Image<RFLOAT> varImg = img;
    varImg.data.initZeros();

    Image<RFLOAT> ccD(s,s);
    Image<Complex> imgDs(sh,s), ccDs(sh,s);

    const int N = 10000;
    for (int i = 0; i < N; i++) {
        if (i % 10 == 0) std::cout << i << "\n";

        for (long int yy = 0; yy < s; yy++)
        for (long int xx = 0; xx < s; xx++) {
            double v = direct::elem(img.data, xx, yy);
            direct::elem(imgD.data, xx, yy) = DistributionHelper::sampleGauss(v, sig);
        }

        ft.FourierTransform(imgD(), imgDs());

        for (long int yy = 0; yy < s;  yy++)
        for (long int xx = 0; xx < sh; xx++) {
            direct::elem(ccDs.data, xx, yy) = 
                direct::elem(spec.data, xx, yy)
              * direct::elem(imgDs.data, xx, yy).conj();
        }

        ft.inverseFourierTransform(ccDs.data, ccD.data);

        for (long int yy = 0; yy < s; yy++)
        for (long int xx = 0; xx < s; xx++) {
            double m0 = direct::elem(mu.data, xx, yy);
            double md = direct::elem(ccD.data, xx, yy);
            double d = md - m0;

            direct::elem(varImg.data, xx, yy) += d * d / N;
        }

    }

    double varSum = 0.0;

    for (long int yy = 0; yy < s; yy++)
    for (long int xx = 0; xx < s; xx++) {
        varSum += direct::elem(varImg.data, xx, yy);
    }

    varSum /= s * s;

    std::cout << varSum << " vs. " << sig2 * varScale << "\n";
}

void NoiseHelper::testColorVariance(Image<RFLOAT> img, std::vector<double> sig2) {
    const int s = img.data.xdim;
    const int sh = s / 2 + 1;

    Image<Complex> spec(sh, s);

    FourierTransformer ft;
    ft.FourierTransform(img(), spec());

    VtkHelper::writeVTK(spec, "debug/spec.vtk");

    for (long int y = 0; y < s;  y++)
    for (long int x = 0; x < sh; x++) {
        if (x == 0 && y == 0) continue;

        const double yy = y < sh ? y : y - s;
        const double xx = x;

        const int r = sqrt(xx * xx + yy * yy);

        if (r >= sh) direct::elem(spec.data, x, y) = Complex(0.0, 0.0);
    }


    double varPred = 0.0;

    for (long int y = 0; y < s;  y++)
    for (long int x = 0; x < sh; x++) {
        if (x == 0 && y == 0) continue;

        const double yy = y < sh ? y : y - s;
        const double xx = x;

        const int r = sqrt(xx * xx + yy * yy);

        if (r < sh && r > 0) {
            double a = direct::elem(spec.data, x, y).norm() / sig2[r];
            varPred += x == 0 ? a : 2.0 * a;
        }
    }

    // varPred *= s * s;

    std::cout << "pred: " << varPred << "\n";

    std::vector<double> sig(sh);

    for (int i = 0; i < sh; i++) { sig[i] = sqrt(sig2[i]); }

    Image<RFLOAT> varImg = img;
    varImg.data.initZeros();

    Image<RFLOAT> ccD(s, s);
    Image<Complex> ccDs(sh, s);

    const int N = 10000;
    const double sqrtH = sqrt(0.5);

    double varTest = 0.0;

    for (int i = 0; i < N; i++) {
        if (i % 100 == 0) std::cout << i << "\n";

        for (long int y = 0; y < s;  y++)
        for (long int x = 0; x < sh; x++) {
            const double yy = y < sh ? y : y - s;
            const double xx = x;

            const int r = sqrt(xx * xx + yy * yy);

            if (r < sh && r > 0) {
                Complex z0 = direct::elem(spec.data, x, y);

                double r0 = DistributionHelper::sampleGauss(0, sqrtH);
                double r1 = DistributionHelper::sampleGauss(0, sqrtH);
                Complex z1 = Complex(r0, r1) * z0 / sig[r];

                direct::elem(ccDs.data, x, y) = z1;

                if (x == 0 && y >= sh) {
                    direct::elem(ccDs.data, x, y) = direct::elem(ccDs.data, x, s - y).conj();
                }

                varTest += x == 0 ? direct::elem(ccDs.data, x, y).norm() :
                              2.0 * direct::elem(ccDs.data, x, y).norm();
            } else {
                direct::elem(ccDs.data, x, y) = 0.0;
            }
        }

        ft.inverseFourierTransform(ccDs.data, ccD.data);

        for (long int y = 0; y < s; y++)
        for (long int x = 0; x < s; x++) {
            double d = direct::elem(ccD.data, x, y);

            direct::elem(varImg.data, x, y) += d * d / (N * varPred);
        }
    }

    varTest /= N;

    VtkHelper::writeVTK(varImg, "debug/varImg_cn.vtk");

    double varSum = 0.0;

    for (long int yy = 0; yy < s; yy++)
    for (long int xx = 0; xx < s; xx++) {
        varSum += direct::elem(varImg.data, xx, yy);
    }

    varSum /= s * s;

    std::cout << varSum << " @ " << varPred << " vs. " << varTest << "\n";

}

void NoiseHelper::testParseval() {
    const int s = 512;
    const int sh = s / 2 + 1;

    Image<RFLOAT> real(s,s);
    Image<Complex> freq(sh,s);

    for (int i = 0; i < 10; i++) {
        double varr = 0.0;

        for (int y = 0; y < s; y++)
        for (int x = 0; x < s; x++) {
            direct::elem(real.data, x, y) = DistributionHelper::sampleGauss(0, 2);

            double v = direct::elem(real.data, x, y);

            varr += v * v;
        }

        varr /= s * s;

        FourierTransformer ft;
        ft.FourierTransform(real(), freq());

        double var = 0.0;

        for (int y = 0; y < s;  y++)
        for (int x = 0; x < sh; x++) {
            const Complex z = direct::elem(freq.data, x, y);
            var += x == 0 ? z.norm() : 2 * z.norm();
        }

        var /= s * s;

        std::cout << varr << " vs. " << var << " (" << var * (double) (s * s) << ")\n";
        // varr = varf*A
    }

    std::cout << "\n";

    for (int i = 0; i < 10; i++) {
        double varf = 0.0;

        for (int y = 0; y < s;  y++)
        for (int x = 0; x < sh; x++) {
            direct::elem(freq.data, x, y).real = DistributionHelper::sampleGauss(0, sqrt(2.0));
            if (x > 0 && x < sh - 1) {
                direct::elem(freq.data, x, y).imag = DistributionHelper::sampleGauss(0, sqrt(2.0));
            } else {
                direct::elem(freq.data, x, y).imag = 0.0;
            }

            varf += x == 0 ? direct::elem(freq.data, x, y).norm() :
                             direct::elem(freq.data, x, y).norm();
        }

        varf /= s * s;

        FourierTransformer ft;
        ft.inverseFourierTransform(freq(), real());

        double var = 0.0;

        for (int y = 0; y < s; y++)
        for (int x = 0; x < s; x++) {
            double v = direct::elem(real.data, x, y);
            var += v * v;
        }

        var /= s * s;

        std::cout << varf << " vs. " << var << " (" << var / (double) (s * s) << ")\n";
        // varr/A = varf
    }
}
