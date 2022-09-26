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

#include <src/jaz/slice_helper.h>
#include <src/jaz/img_proc/filter_helper.h>
#include <src/jaz/interpolation.h>

using namespace gravis;

void SliceHelper::affineTransform(const Image<RFLOAT> &img, d4Matrix A, Image<RFLOAT> &dest) {
    d4Matrix Ai = A;
    Ai.invert();

    dest.data.resize(1, 1, img.data.ydim, img.data.xdim);

    for (long int y = 0; y < dest.data.ydim; y++)
    for (long int x = 0; x < dest.data.xdim; x++) {
        d4Vector s0(x,y,0,1);
        d4Vector s1 = Ai * s0;

        direct::elem(dest.data, x, y) = Interpolation::linearXY(img, s1.x, s1.y, 0);
    }
}

Image<RFLOAT> SliceHelper::downsample(Image<RFLOAT> &img, int factor) {
    Image<RFLOAT> dest(img.data.xdim / factor, img.data.ydim / factor);
    double q = dest.data.xdim / (double) img.data.xdim;

    Image<RFLOAT> slice0(img.data.xdim, img.data.ydim, 1);

    FilterHelper::lowPassFilter(img, 0.9*q, q, slice0);
    subsample(slice0, dest);
    return dest;
}

void SliceHelper::downsampleSlices(const Image<RFLOAT> &img, Image<RFLOAT> &dest) {
    double q = dest.data.xdim / (double) img.data.xdim;

    Image<RFLOAT> slice0(img.data.xdim, img.data.ydim, 1);
    Image<RFLOAT> slice1(dest.data.xdim, dest.data.ydim, 1);

    for (long int n = 0; n < img.data.ndim; n++) {
        std::cout << n << "/" << img.data.ndim << "\n";

        auto slice0 = getStackSlice(img, n);
        FilterHelper::lowPassFilter(slice0, 0.9*q, q, slice0);
        subsample(slice0, slice1);
        insertStackSlice(slice1, dest, n);
    }
}

void SliceHelper::downsampleSlicesReal(const Image<RFLOAT> &img, Image<RFLOAT> &dest) {
    double q = dest.data.xdim / (double) img.data.xdim;

    Image<RFLOAT> sliceT(img.data.xdim, img.data.ydim, 1);
    Image<RFLOAT> slice1(dest.data.xdim, dest.data.ydim, 1);

    for (long int n = 0; n < img.data.ndim; n++) {
        auto slice0 = getStackSlice(img, n);
        FilterHelper::separableGaussianXYZ(slice0, sliceT, 1.5/q);
        subsample(sliceT, slice1);
        insertStackSlice(slice1, dest, n);
    }
}

void SliceHelper::lowPassFilterSlicewise(Image<RFLOAT> &img, double maxFreq0, double maxFreq1) {

    for (long int n = 0; n < img.data.ndim; n++) {
        lowPassFilterSlice(img, n, maxFreq0, maxFreq1);
    }
}

void SliceHelper::lowPassFilterSlice(Image<RFLOAT> &img, long int n, double maxFreq0, double maxFreq1) {
    auto slice0 = getStackSlice(img, n);
    FilterHelper::lowPassFilter(slice0, maxFreq0, maxFreq1, slice0);
    insertStackSlice(slice0, img, n);
}

void SliceHelper::subsample(const Image<RFLOAT> &img, Image<RFLOAT> &dest) {
    double q = img.data.xdim / (double) dest.data.xdim;

    for (long int y = 0; y < dest.data.ydim; y++)
    for (long int x = 0; x < dest.data.xdim; x++) {
        direct::elem(dest.data, x, y) = direct::elem(img.data, (long int)(q*x + 0.5), (long int)(q*y + 0.5));
    }
}

template <int D>
double average(const MultidimArray<RFLOAT> &arr);

template <>
double average<2>(const MultidimArray<RFLOAT> &arr) {
    double sum = 0.0;
    for (long int y = 0; y < arr.ydim; y++)
    for (long int x = 0; x < arr.xdim; x++) {
        sum += direct::elem(arr, x, y);
    }
    return sum / (arr.xdim * arr.ydim);
}

template <>
double average<3>(const MultidimArray<RFLOAT> &arr) {
    double sum = 0.0;
    for (long int z = 0; z < arr.zdim; z++)
    for (long int y = 0; y < arr.ydim; y++)
    for (long int x = 0; x < arr.xdim; x++) {
        sum += direct::elem(arr, x, y, z);
    }
    return sum / (arr.xdim * arr.ydim * arr.zdim);
}

namespace SliceHelper {

template <>
void avgPad<2>(const Image<RFLOAT> &src, Image<RFLOAT> &dest, double ratio) {

    const int padX = ratio * src.data.xdim;
    const int padY = ratio * src.data.ydim;

    const double avg = average<2>(src.data);

    dest.data.resize(src.data.xdim + 2 * padX, src.data.ydim + 2 * padY);
    dest.data.fill(avg);

    for (long int y = 0; y < src.data.ydim; y++)
    for (long int x = 0; x < src.data.xdim; x++) {
        direct::elem(dest.data, x + padX, y + padY) = direct::elem(src.data, x, y);
    }
}

template <>
void avgPad<3>(const Image<RFLOAT> &src, Image<RFLOAT> &dest, double ratio) {

    const int padX = ratio * src.data.xdim;
    const int padY = ratio * src.data.ydim;
    const int padZ = ratio * src.data.zdim;

    const double avg = average<3>(src.data);

    dest.data.resize(src.data.xdim + 2 * padX, src.data.ydim + 2 * padY, src.data.zdim + 2 * padZ);
    dest.data.fill(avg);

    for (long int z = 0; z < src.data.zdim; z++)
    for (long int y = 0; y < src.data.ydim; y++)
    for (long int x = 0; x < src.data.xdim; x++) {
        direct::elem(dest.data, x + padX, y + padY, z + padZ) = direct::elem(src.data, x, y, z);
    }
}

};

void SliceHelper::halveSpectrum2D(const Image<Complex> &src, Image<Complex> &dest) {

    const int wd = src.data.xdim / 2 + 1;
    const int hd = src.data.ydim;

    dest.data.resize(wd, hd);

    const int xo = src.data.xdim / 2 + 1;
    const int yo = src.data.ydim / 2 + 1;

    for (long int y = 0; y < hd; y++)
    for (long int x = 0; x < wd; x++) {
        /*if (x == 0) {
            direct::elem(dest.data, 0, y) = direct::elem(src.data, xo, y);
        } else if (xo + x < src.data.xdim) {
            direct::elem(dest.data, x, y) = 0.5 * (direct::elem(src.data, xo + x, y)
                                                      + direct::elem(src.data, xo - x, yo - (y - yo)));
        } else {
            direct::elem(dest.data, x, y) = direct::elem(src.data, xo - x, yo - (y - yo));
        }*/

        const int yr = y;
        const int yw = (yr + yo) % hd;

        direct::elem(dest.data, x, y) = direct::elem(src.data, xo + x, yw);
    }
}

void SliceHelper::extractSpectralSlice(
    Image<Complex> &src, Image<RFLOAT> &dest, d3Matrix proj,
    d2Vector volCentImg, double oversample
) {
    const int wi = dest.data.xdim;
    const int hi = dest.data.ydim;

    const double wv = src.data.xdim;
    const double hv = src.data.ydim;
    const double dv = src.data.zdim;

    const double wios = oversample * wi;
    const double hios = oversample * hi;

    const int wiosI = (int) wios / 2 + 1;
    const int hiosI = (int) hios;

    const int ciosX = (int) wios / 2;
    const int ciosY = (int) hios / 2;

    MultidimArray<Complex> dest2  (wiosI, hiosI);
    MultidimArray<Complex> weight (wiosI, hiosI);

    d2Vector shift (volCentImg.x - ciosX, volCentImg.y - ciosY);

    for (long int y = 0; y < dest2.ydim; y++)
    for (long int x = 0; x < dest2.xdim; x++) {
        d3Vector pi((double)x/(double)(wiosI-1), 2.0*(double)y/(double)hiosI, 0.0);

        if (pi.y >= 1.0) pi.y -= 2.0;

        if (pi.norm2() > 1.0) {
            direct::elem(dest2, x, y) = Complex(0, 0);
            continue;
        }

        d3Vector pv = proj * pi;

        bool conj = false;

        if (pv.x < 0.0) {
            pv = -pv;
            conj = true;
        }

        if (pv.norm2() > 1.0) {
            direct::elem(dest2, x, y) = Complex(0,0);
            continue;
        }

        double xxd = (wv - 1) * pv.x;
        double yyd = hv * pv.y / 2.0;
        double zzd = dv * pv.z / 2.0;

        double ax = std::abs(xxd);
        double ay = std::abs(yyd);
        double az = std::abs(zzd);

        double phi = - PI * (pi.x * shift.x + pi.y * shift.y);
        const auto z0 = Complex::unit(phi);

        if (ax < 1.0 && ay < 1.0 && az < 1.0) {
            direct::elem(weight, x, y) = z0 * Complex((1.0 - ax) * (1.0 - ay) * (1.0 - az), 0.0);
        } else {
            direct::elem(weight, x, y) = Complex(0.0, 0.0);
        }

        if (yyd < 0.0) yyd += hv;
        if (zzd < 0.0) zzd += dv;

        direct::elem(dest2, x, y) = z0 * (
            conj ? Interpolation::linearFFTW3D(src, xxd, yyd, zzd).conj() :
                   Interpolation::linearFFTW3D(src, xxd, yyd, zzd)
        );
    }

    FourierTransformer transformer;
    auto dest2r  = transformer.inverseFourierTransform(dest2);
    auto weightr = transformer.inverseFourierTransform(weight);

    // In principle, we could center dest.data as we assign into it
    for (long int i = 0; i < dest.data.size(); i++)
        dest.data[i] = dest2r[i] / weightr[i];

    CenterFFT(dest.data, +1);

}

struct insert_spectral_slices_stage {

    std::vector<Image<Complex>> srcSpectra;
    std::vector<d2Vector> shifts;
    std::vector<double> thz;
    double wif, hif;

    insert_spectral_slices_stage(
        int ic,
        std::vector<Image<RFLOAT>> &src,
        const std::vector<gravis::d2Vector> &volCentImg,
        const std::vector<gravis::d3Matrix> &proj,
        double wv, double hv, double dv,
        double wir, double hir,
        double imgPad
    ): srcSpectra(0), shifts(0), thz(0) {
        exec(ic, src, volCentImg, proj, wv, hv, dv, wir, hir, imgPad);
    }

    private:

    void exec(
        int ic,
        std::vector<Image<RFLOAT>> &src,
        const std::vector<gravis::d2Vector> &volCentImg,
        const std::vector<gravis::d3Matrix> &proj,
        double wv, double hv, double dv,
        double wir, double hir,
        double imgPad
    ) {
        srcSpectra.reserve(ic);
        shifts.reserve(ic);
        thz.reserve(ic);
        Image<RFLOAT> img;
        FourierTransformer transformer;
        for (int i = 0; i < ic; i++) {
            SliceHelper::avgPad<2>(src[i], img, imgPad);

            CenterFFT(img.data, -1);
            srcSpectra.emplace_back(transformer.FourierTransform(img.data));

            shifts.emplace_back(volCentImg[i].x - wir / 2.0, volCentImg[i].y - hir / 2.0);

            thz.push_back(0.5 * d3Vector(wv * proj[i](2, 0), hv * proj[i](2, 1), dv * proj[i](2, 2)).length());
        }
        wif = srcSpectra[0].data.xdim;
        hif = srcSpectra[0].data.ydim;
    }

};


void SliceHelper::insertSpectralSlices(
    std::vector<Image<RFLOAT>> &src,
    std::vector<gravis::d3Matrix> proj,
    std::vector<gravis::d2Vector> volCentImg,
    Image<Complex> &dest, double thickness, double thicknessSlope, double imgPad
) {
    const double wv = dest.data.xdim;
    const double hv = dest.data.ydim;
    const double dv = dest.data.zdim;

    const double wir = src[0].data.xdim;
    const double hir = src[0].data.ydim;

    const int ic = src.size();

    insert_spectral_slices_stage subroutine
        (ic, src, volCentImg, proj, wv, hv, dv, wir, hir, imgPad);

    const auto &srcSpectra = subroutine.srcSpectra;
    const auto &shifts     = subroutine.shifts;
    const auto &thz        = subroutine.thz;

    const double &wif = subroutine.wif;
    const double &hif = subroutine.hif;

    for (long int z = 0; z < dest.data.zdim; z++)
    for (long int y = 0; y < dest.data.ydim; y++)
    for (long int x = 0; x < dest.data.xdim; x++) {

        d3Vector pv ((double) x / (wv - 1), 2.0 * (double) y / hv, 2.0 * (double) z / dv);

        if (pv.y > 1.0) pv.y -= 2.0;
        if (pv.z > 1.0) pv.z -= 2.0;

        if (pv.norm2() >= 1.0) {
            direct::elem(dest.data, x, y, z) = Complex(0, 0);
            continue;
        }

        const double r = euclid(pv.x, pv.z);

        Complex zs (0.0, 0.0);
        double wgh = 0.0;

        for (int i = 0; i < ic; i++) {

            d3Vector pi3 = proj[i] * pv;

            if (euclidsq(pi3.x, pi3.y) >= 1.0) continue;

            const double za = thz[i] * std::abs(pi3.z);

            const double th_r = thickness + r * thz[i] * thicknessSlope;
            if (za > th_r) continue;

            bool conj = false;

            if (pi3.x < 0.0) {
                pi3 = -pi3;
                conj = true;
            }

            double xi = (wif - 1) * pi3.x;
            double yi = hif * pi3.y / 2.0;

            if (yi < 0.0) yi += hif;

            const double phi = PI * (pi3.x * shifts[i].x + pi3.y * shifts[i].y);
            const Complex z0 (cos(phi), sin(phi));

            // double wgi = (1.0 - za / th_r) * (thickness / th_r);
            const double wgi = 1.0 - za / th_r;
            const Complex zz = z0 * Interpolation::linearFFTW2D(srcSpectra[i], xi, yi);

            zs += wgi * (conj ? zz.conj() : zz);
            wgh += wgi;
        }

        if (wgh > 1.0) zs /= wgh;

        direct::elem(dest.data, x, y, z) = zs;
    }
}


void SliceHelper::insertWeightedSpectralSlices(
    std::vector<Image<RFLOAT>> &src,
    std::vector<gravis::d3Matrix> proj,
    std::vector<gravis::d2Vector> volCentImg,
    std::vector<double> imgWeights,
    Image<Complex> &dest, double thickness, double imgPad
) {
    const double wv = dest.data.xdim;
    const double hv = dest.data.ydim;
    const double dv = dest.data.zdim;

    const double wir = src[0].data.xdim;
    const double hir = src[0].data.ydim;

    const int ic = src.size();

    insert_spectral_slices_stage subroutine
        (ic, src, volCentImg, proj, wv, hv, dv, wir, hir, imgPad);

    const auto &srcSpectra = subroutine.srcSpectra;
    const auto &shifts     = subroutine.shifts;
    const auto &thz        = subroutine.thz;

    const double &wif = subroutine.wif;
    const double &hif = subroutine.hif;

    for (long int z = 0; z < dest.data.zdim; z++)
    for (long int y = 0; y < dest.data.ydim; y++)
    for (long int x = 0; x < dest.data.xdim; x++) {

        d3Vector pv ((double) x / (wv - 1), 2.0 * (double) y / hv, 2.0 * (double) z / dv);

        if (pv.y > 1.0) pv.y -= 2.0;
        if (pv.z > 1.0) pv.z -= 2.0;

        if (pv.norm2() >= 1.0) {
            direct::elem(dest.data, x, y, z) = Complex(0, 0);
            continue;
        }

        Complex zs (0.0, 0.0);
        double wgh = 0.0;

        for (int i = 0; i < ic; i++) {

            d3Vector pi3 = proj[i] * pv;

            if (euclidsq(pi3.x, pi3.y) >= 1.0) continue;

            const double za = thz[i] * std::abs(pi3.z);

            const double th_r = thickness;
            if (za > th_r) continue;

            bool conj = false;

            if (pi3.x < 0.0) {
                pi3 = -pi3;
                conj = true;
            }

            double xi = (wif - 1) * pi3.x;
            double yi = hif * pi3.y / 2.0;

            if (yi < 0.0) yi += hif;

            const double phi = PI * (pi3.x * shifts[i].x + pi3.y * shifts[i].y);
            const Complex z0 (cos(phi), sin(phi));

            const double wgi = imgWeights[i] * (1.0 - za / th_r);
            const Complex zz = z0 * Interpolation::linearFFTW2D(srcSpectra[i], xi, yi);

            zs += wgi * (conj ? zz.conj() : zz);
            wgh += wgi;
        }

        if (wgh > 1.0) zs /= wgh;

        direct::elem(dest.data, x, y, z) = zs;
    }
}

Image<RFLOAT> SliceHelper::getStackSlice(const Image<RFLOAT> &src, long int s) {
    const long int w = src.data.xdim;
    const long int h = src.data.ydim;
    Image<RFLOAT> dest (w, h);
    for (long int y = 0; y < h; y++)
    for (long int x = 0; x < w; x++) {
        direct::elem(dest.data, x, y) = direct::elem(src.data, x, y, 0, s);
    }
    return dest;
}

template <typename T>
Image<RFLOAT> SliceHelper::getStackSlices(const Image<T> &src, long int s, long int ndim) {
    const long int w = src.data.xdim;
    const long int h = src.data.ydim;
    Image<RFLOAT> dest (w, h, 1, ndim);
    for (long int n = 0; n < ndim; n++)
    for (long int y = 0; y < h; y++)
    for (long int x = 0; x < w; x++) {
        direct::elem(dest.data, x, y, 0, n) = direct::elem(src.data, x, y, 0, n + s);
    }
    return dest;
}

template <typename T>
void SliceHelper::insertStackSlice(const Image<T> &src, Image<T> &dest, long int n) {
    if (src.data.xdim != dest.data.xdim || src.data.ydim != dest.data.ydim) {
        REPORT_ERROR("SliceHelper::extractSlice: image size mismatch.\n");
    }

    for (long int y = 0; y < src.data.ydim; y++)
    for (long int x = 0; x < src.data.xdim; x++) {
        direct::elem(dest.data, x, y, 0, n) = direct::elem(src.data, x, y);
    }
}

template <typename T>
void SliceHelper::insertZSlice(const Image<T> &src, Image<T> &dest, long int z) {
    if (src.data.xdim != dest.data.xdim || src.data.ydim != dest.data.ydim) {
        REPORT_ERROR("SliceHelper::extractSlice: image size mismatch.\n");
    }

    for (long int y = 0; y < src.data.ydim; y++)
    for (long int x = 0; x < src.data.xdim; x++) {
        direct::elem(dest.data, x, y, z, 0) = direct::elem(src.data, x, y);
    }
}

template <typename T>
Image<T> SliceHelper::consolidate(const std::vector<Image<T>> &src, bool toN) {
    const int w = src[0].data.xdim;
    const int h = src[0].data.ydim;
    const int ic = src.size();

    const int zc = toN? 1 : ic;
    const int nc = toN? ic : 1;

    Image<double> out (w, h, zc, nc);

    for (int i = 0; i < ic; i++) {
        if (src[i].data.xdim != w || src[i].data.ydim != h) {
            REPORT_ERROR("SliceHelper::consolidate(): images are of unequal size.\n");
        }

        if (toN) insertStackSlice(src[i], out, i);
        else         insertZSlice(src[i], out, i);
    }

    return out;
}

void SliceHelper::stat(const Image<RFLOAT> &img) {
    std::cout << "xdim: " << img.data.xdim << "\n"
              << "ydim: " << img.data.ydim << "\n"
              << "zdim: " << img.data.zdim << "\n"
              << "ndim: " << img.data.ndim << "\n";
}
