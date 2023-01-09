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

#include <src/jaz/img_proc/filter_helper.h>
#include <src/jaz/image_log.h>
#include <src/jaz/index_sort.h>
#include <src/jaz/gravis/t2Vector.h>
#include <src/jaz/tensor2x2.h>
#include <src/jaz/interpolation.h>
#include "src/jaz/ctf_helper.h"

#include <limits>

extern "C" {
    #include <src/jaz/d3x3/dsyev2.h>
}

using namespace gravis;

template<typename T>
std::vector<T> make_kernel(T sigma, int n) {
    std::vector<T> kernel;
    kernel.reserve(2 * n + 1);
    const T s2 = sigma * sigma;
    for (int i = 0; i != 2 * n + 1; i++)
        kernel.push_back(exp(-0.5 * (i - n) * (i - n) / s2));
    return kernel;
}

template <typename valarray_like>
inline void normalise(valarray_like& kernel) {
    const double sum = std::accumulate(kernel.begin(), kernel.end(), 0.0);
    for (double& x: kernel) x /= sum;  // Now kernel sums to 1
}

template <typename T, typename valarray_like>
void convolve_x_nowrap(const Volume<T>& src, const valarray_like& kernel, Volume<T>& dest) {
    const size_t mid = kernel.size() / 2;
    dest.resize(src);
    for (size_t k = 0; k < src.dimz; k++)
    for (size_t j = 0; j < src.dimy; j++)
    for (size_t i = 0; i < src.dimx; i++) {

        T v = 0;
        double m = 0;
        for (long int l = 0; l != kernel.size(); l++) {
            const long int ii = i + l - mid;
            if (ii < 0 || ii >= src.dimx) continue;
            v += kernel[l] * src(ii, j, k);
            m += kernel[l];
        }

        dest(i, j, k) = v / m;
    }
}

template <typename T, typename valarray_like>
void convolve_x_nowrap(const MultidimArray<T>& src, const valarray_like& kernel, MultidimArray<T>& dest) {
    const size_t mid = kernel.size() / 2;
    dest.resize(src);
    for (size_t k = 0; k < src.zdim; k++)
    for (size_t j = 0; j < src.ydim; j++)
    for (size_t i = 0; i < src.xdim; i++) {

        T v = 0;
        double m = 0;
        for (long int l = 0; l != kernel.size(); l++) {
            const long int ii = i + l - mid;
            if (ii < 0 || ii >= src.xdim) continue;
            v += kernel[l] * direct::elem(src, ii, j, k);
            m += kernel[l];
        }

        direct::elem(dest, i, j, k) = v / m;
    }
}

template <typename T, typename valarray_like>
void convolve_x_wrap(const MultidimArray<T>& src, const valarray_like& kernel, MultidimArray<T>& dest) {
    const size_t mid = kernel.size() / 2;
    dest.resize(src);
    for (size_t k = 0; k < src.zdim; k++)
    for (size_t j = 0; j < src.ydim; j++)
    for (size_t i = 0; i < src.xdim; i++) {

        T v = 0;
        for (int l = 0; l < kernel.size(); l++) {
            const int ii = (i + l - mid + src.xdim) % src.xdim;
            v += kernel[l] * direct::elem(src, ii, j, k);
        }

        direct::elem(dest, i, j, k) = v;
    }
}

template <typename T, typename valarray_like>
void convolve_y_nowrap(const MultidimArray<T>& src, const valarray_like& kernel, MultidimArray<T>& dest) {
    const size_t mid = kernel.size() / 2;
    dest.resize(src);
    for (size_t k = 0; k < src.zdim; k++)
    for (size_t j = 0; j < src.ydim; j++)
    for (size_t i = 0; i < src.xdim; i++) {

        T v = 0;
        double m = 0;
        for (long int l = 0; l != kernel.size(); l++) {
            const long int jj = j - mid + l;
            if (jj < 0 || jj >= src.ydim) continue;
            v += kernel[l] * direct::elem(src, i, jj, k);
            m += kernel[l];
        }

        direct::elem(dest, i, j, k) = v / m;
    }
}

template <typename T, typename valarray_like>
void convolve_y_nowrap(const Volume<T>& src, const valarray_like& kernel, Volume<T>& dest) {
    const size_t mid = kernel.size() / 2;
    dest.resize(src);
    for (size_t k = 0; k < src.dimz; k++)
    for (size_t j = 0; j < src.dimy; j++)
    for (size_t i = 0; i < src.dimx; i++) {

        T v = 0;
        double m = 0;
        for (long int l = 0; l != kernel.size(); l++) {
            const long int jj = j - mid + l;
            if (jj < 0 || jj >= src.dimy) continue;
            v += kernel[l] * src(i, jj, k);
            m += kernel[l];
        }

        dest(i, j, k) = v / m;
    }
}

template <typename T, typename valarray_like>
void convolve_y_wrap(const MultidimArray<T>& src, const valarray_like& kernel, MultidimArray<T>& dest) {
    const size_t mid = kernel.size() / 2;
    dest.resize(src);
    for (size_t k = 0; k < src.zdim; k++)
    for (size_t j = 0; j < src.ydim; j++)
    for (size_t i = 0; i < src.xdim; i++) {

        T v = 0;
        for (long int l = 0; l != kernel.size(); l++) {
            const long int jj = (j - mid + l + src.ydim) % src.ydim;
            v += kernel[l] * direct::elem(src, i, jj, k);
        }

        direct::elem(dest, i, j, k) = v;
    }
}

template <typename T, typename valarray_like>
void convolve_z_nowrap(const MultidimArray<T>& src, const valarray_like& kernel, MultidimArray<T>& dest) {
    const size_t mid = kernel.size() / 2;
    dest.resize(src);
    for (size_t k = 0; k < src.zdim; k++)
    for (size_t j = 0; j < src.ydim; j++)
    for (size_t i = 0; i < src.xdim; i++) {

        T v = 0;
        double m = 0;
        for (long int l = 0; l != kernel.size(); l++) {
            const long int kk = k - mid + l;
            if (kk < 0 || kk >= src.zdim) continue;
            v += kernel[l] * direct::elem(src, i, j, kk);
            m += kernel[l];
        }

        direct::elem(dest, i, j, k) = v / m;
    }
}

template <typename T, typename valarray_like>
void convolve_z_nowrap(const Volume<T>& src, const valarray_like& kernel, Volume<T>& dest) {
    const size_t mid = kernel.size() / 2;
    dest.resize(src);
    for (size_t k = 0; k < src.dimz; k++)
    for (size_t j = 0; j < src.dimy; j++)
    for (size_t i = 0; i < src.dimx; i++) {

        T v = 0;
        double m = 0;
        for (long int l = 0; l != kernel.size(); l++) {
            const long int kk = k - mid + l;
            if (kk < 0 || kk >= src.dimz) continue;
            v += kernel[l] * src(i, j, kk);
            m += kernel[l];
        }

        dest(i, j, k) = v / m;
    }
}

template <typename T, typename valarray_like>
void convolve_z_wrap(const MultidimArray<T>& src, const valarray_like& kernel, MultidimArray<T>& dest) {
    const size_t mid = kernel.size() / 2;
    dest.resize(src);
    for (size_t k = 0; k < src.zdim; k++)
    for (size_t j = 0; j < src.ydim; j++)
    for (size_t i = 0; i < src.xdim; i++) {

        T v = 0;
        for (long int l = 0; l != kernel.size(); l++) {
            const long int kk = (k - mid + l + src.zdim) % src.zdim;
            v += kernel[l] * direct::elem(src, i, j, kk);
        }

        direct::elem(dest, i, j, k) = v;
    }
}

Image<RFLOAT> FilterHelper::separableGaussianX(
    const Image<RFLOAT> &src, RFLOAT sigma, int n, bool wrap
) {
    if (n < 0) { n = 2 * sigma + 0.5; }
    auto kernel = make_kernel(sigma, n);
    MultidimArray<RFLOAT> work1;
    if (wrap) {
        normalise(kernel);
        convolve_x_wrap(src.data, kernel, work1);
    } else {
        convolve_x_nowrap(src.data, kernel, work1);
    }
    return {std::move(work1)};
}

Image<RFLOAT> FilterHelper::separableGaussianXY(
    const Image<RFLOAT> &src, RFLOAT sigma, int n, bool wrap
) {
    if (n < 0) { n = 2 * sigma + 0.5; }
    auto kernel = make_kernel(sigma, n);
    MultidimArray<RFLOAT> work1, work2;
    if (wrap) {
        normalise(kernel);
        convolve_x_wrap(src.data, kernel, work1);
        convolve_y_wrap(work1,    kernel, work2);
    } else {
        convolve_x_nowrap(src.data, kernel, work1);
        convolve_y_nowrap(work1,    kernel, work2);
    }
    return {std::move(work2)};
}

Image<RFLOAT> FilterHelper::separableGaussianXYZ(
    const Image<RFLOAT> &src, RFLOAT sigma, int n, bool wrap
) {
    if (n < 0) { n = 2 * sigma + 0.5; }
    auto kernel = make_kernel(sigma, n);
    MultidimArray<RFLOAT> work1, work2;
    if (wrap) {
        normalise(kernel);
        convolve_x_wrap(src.data, kernel, work1);
        convolve_y_wrap(work1,    kernel, work2);
        convolve_z_wrap(work2,    kernel, work1);
    } else {
        convolve_x_nowrap(src.data, kernel, work1);
        convolve_y_nowrap(work1,    kernel, work2);
        convolve_z_nowrap(work2,    kernel, work1);
    }
    return {std::move(work1)};
}

template <typename C, typename valarray_like>
void convolve_x_fftw(const MultidimArray<C>& src, const valarray_like& kernel, MultidimArray<C>& dest) {
    const size_t mid = kernel.size() / 2;
    dest.resize(src);
    for (size_t k = 0; k < src.zdim; k++)
    for (size_t j = 0; j < src.ydim; j++)
    for (size_t i = 0; i < src.xdim; i++) {

        C v = 0;
        for (long int l = 0; l != kernel.size(); l++) {
            long ii = i - mid + l, jj = j, kk = k;
            if (ii < 0) {
                ii = -ii;
                jj = src.ydim - j % src.ydim;
                kk = src.zdim - k % src.zdim;
            } else if (ii >= src.xdim) {
                ii = 2 * src.xdim - 1 - ii;
                jj = src.ydim - j % src.ydim;
                kk = src.zdim - k % src.zdim;
            }
            v += kernel[l] * direct::elem(src, ii, jj, kk);
        }

        direct::elem(dest, i, j, k) = v;
    }
}

template <typename C, typename valarray_like>
void convolve_x_fftw_2(const MultidimArray<C>& src, const valarray_like& kernel, MultidimArray<C>& dest) {
    const size_t mid = kernel.size() / 2;
    dest.resize(src);
    for (size_t k = 0; k < src.zdim; k++)
    for (size_t j = 0; j < src.ydim; j++)
    for (size_t i = 0; i < src.xdim; i++) {

        C v = 0;
        for (long int l = 0; l != kernel.size(); l++) {

            long ii = i - mid + l, jj = j, kk = k;
            if (ii < 0) {
                ii = -ii;
                jj = src.ydim - j % src.ydim;
                kk = src.zdim - k % src.zdim;
            } else if (ii >= src.xdim) {
                ii = 2 * src.xdim - 2 - ii;  // Surely ii = 2 * src.xdim - 1 - ii?
                jj = src.ydim - j % src.ydim;
                kk = src.zdim - k % src.zdim;
            }
            C vv = direct::elem(src, ii, jj, kk);  // cf convolve_x_fftw
            if (i - mid + l < 0 || i - mid + l >= src.xdim) vv = vv.conj();
            v += kernel[l] * vv;
        }

        direct::elem(dest, i, j, k) = v;
    }
}

template <typename C, typename valarray_like>
void convolve_y_fftw(const MultidimArray<C>& src, const valarray_like& kernel, MultidimArray<C>& dest) {
    const size_t mid = kernel.size() / 2;
    dest.resize(src);
    for (size_t k = 0; k < src.zdim; k++)
    for (size_t j = 0; j < src.ydim; j++)
    for (size_t i = 0; i < src.xdim; i++) {

        C v = 0;
        for (long int l = 0; l != kernel.size(); l++) {
            const long int jj = (j - mid + l + src.ydim) % src.ydim;
            v += kernel[l] * direct::elem(src, i, jj, k);
        }

        direct::elem(dest, i, j, k) = v;
    }
}

template <typename C, typename valarray_like>
void convolve_z_fftw(const MultidimArray<C>& src, const valarray_like& kernel, MultidimArray<C>& dest) {
    const size_t mid = kernel.size() / 2;
    dest.resize(src);
    for (size_t k = 0; k < src.zdim; k++)
    for (size_t j = 0; j < src.ydim; j++)
    for (size_t i = 0; i < src.xdim; i++) {

        C v = 0;
        for (long int l = 0; l != kernel.size(); l++) {
            const long int kk = (k - mid + l + src.zdim) % src.zdim;
            v += kernel[l] * direct::elem(src, i, j, kk);
        }

        direct::elem(dest, i, j, k) = v;
    }
}

void FilterHelper::separableGaussianFreq(
    const MultidimArray<Complex> &src,
    MultidimArray<Complex> &dest, double sigma, int k
) {
    if (k < 0) { k = 2 * sigma + 0.5; }
    auto kernel = make_kernel(sigma, k);

    MultidimArray<Complex> temp;
    convolve_x_fftw(src,  kernel, dest);
    convolve_y_fftw(dest, kernel, temp);
    convolve_z_fftw(temp, kernel, dest);
}

void FilterHelper::separableGaussianFreqXY(
    const MultidimArray<Complex> &src,
    MultidimArray<Complex> &dest, double sigma, int k
) {
    if (k < 0) { k = 2 * sigma + 0.5; }
    auto kernel = make_kernel(sigma, k);

    MultidimArray<Complex> temp;
    convolve_x_fftw_2(src, kernel, temp);
    convolve_y_fftw(temp, kernel, dest);
}

void FilterHelper::drawTestPattern(Image<RFLOAT> &img, int squareSize) {
    for (long int k = 0; k < img.data.zdim; k++)
    for (long int j = 0; j < img.data.ydim; j++)
    for (long int i = 0; i < img.data.xdim; i++) {
        const int x = i / squareSize % 2;
        const int y = j / squareSize % 2;
        const int z = k / squareSize % 2;
        direct::elem(img.data, i, j, k) = (x + y + z) % 2;
    }
}

void FilterHelper::drawTestPattern(Volume<RFLOAT>& volume, int squareSize) {
    for (size_t k = 0; k < volume.dimz; k++)
    for (size_t j = 0; j < volume.dimy; j++)
    for (size_t i = 0; i < volume.dimx; i++) {
        const int x = i / squareSize % 2;
        const int y = j / squareSize % 2;
        const int z = k / squareSize % 2;
        volume(i, j, k) = (x + y + z) % 2;
    }
}

Image<RFLOAT> FilterHelper::expImg(const Image<RFLOAT> &img, double scale) {
    MultidimArray<RFLOAT> out (img.data.xdim, img.data.ydim, img.data.zdim, img.data.ndim);
    for (long int i = 0; i != out.size(); i++)
        out[i] = exp(scale * img.data[i]);
    return {std::move(out)};
}

Image<RFLOAT> FilterHelper::logImg(const Image<RFLOAT> &img, double scale, double thresh) {
    MultidimArray<RFLOAT> out (img.data.xdim, img.data.ydim, img.data.zdim, img.data.ndim);
    for (long int i = 0; i != out.size(); i++)
        out[i] = log(scale * std::max(img.data[i], thresh));
    return {std::move(out)};
}

template <typename R>
Image<R> FilterHelper::padCorner2D(const Image<R> &img, int wF, int hF, R fill) {
    const int w0 = img.data.xdim, h0 = img.data.ydim;
    MultidimArray<R> out (wF, hF);
    for (int j = 0; j < hF; j++)
    for (int i = 0; i < wF; i++) {
        int ii = i < wF / 2 ? i : i - wF;
        int jj = j < hF / 2 ? j : j - hF;
        if (ii < w0 / 2 && jj < h0 / 2 && ii >= -w0 / 2 && jj >= -h0 / 2) {
            if (ii < 0) ii += w0;
            if (jj < 0) jj += h0;
            direct::elem(out, i, j) = direct::elem(img.data, ii, jj);
        } else {
            direct::elem(out, i, j) = fill;
        }
    }
    return {std::move(out)};
}

template <typename C>
Image<C> FilterHelper::padCorner2D_fftw(Image<C> &img, int wF, int hF, C fill) {
    const int w0 = img.data.xdim, h0 = img.data.ydim;
    MultidimArray<C> out (wF, wF);
    for (int j = 0; j < wF; j++)
    for (int i = 0; i < wF; i++) {
        int ii = i;
        int jj = j < wF / 2 ? j : j - wF;
        if (ii < w0 / 2 && jj < h0 / 2 && ii >= -w0 / 2 && jj >= -h0 / 2) {
            if (jj < 0) jj += h0;
            direct::elem(out, i, j) = direct::elem(img.data, ii, jj);
        } else {
            direct::elem(out, i, j) = fill;
        }
    }
    return {std::move(out)};
}

template <typename R>
Image<R> FilterHelper::cropCorner2D(const Image<R> &img, int wF, int hF, R fill) {
    const int w0 = img.data.xdim, h0 = img.data.ydim;
    if (wF > w0 || hF > h0) return img;
    MultidimArray<R> out (wF, hF);
    for (int j = 0; j < h0; j++)
    for (int i = 0; i < w0; i++) {
        int ii = i < w0 / 2 ? i : i - w0;
        int jj = j < h0 / 2 ? j : j - h0;
        if (ii < wF / 2 && jj < hF / 2 && ii >= -wF / 2 && jj >= -hF / 2) {
            if (ii < 0) ii += wF;
            if (jj < 0) jj += hF;
            direct::elem(out, ii, jj) = direct::elem(img.data, i, j);
        }
    }
    return {std::move(out)};
}

template <typename C>
Image<C> FilterHelper::cropCorner2D_fftw(const Image<C> &img, int wF, int hF, C fill) {
    const int w0 = img.data.xdim, h0 = img.data.ydim;
    MultidimArray<C> out (wF, hF);
    for (int j = 0; j < h0; j++)
    for (int i = 0; i < w0; i++) {
        int ii = i;
        int jj = j < h0 / 2 ? j : j - h0;
        if (ii < wF && jj < hF / 2 && jj >= -hF / 2) {
            if (jj < 0) jj += hF;
            direct::elem(out, ii, jj) = direct::elem(img.data, i, j);
        }
    }
    return {std::move(out)};
}

Image<RFLOAT> FilterHelper::zeroOutsideCorner2D(const Image<RFLOAT> &img, double radius) {
    const double r2 = radius * radius;
    const int w = img.data.xdim, h = img.data.ydim;
    MultidimArray<RFLOAT> out (w, h);
    for (int j = 0; j < h; j++)
    for (int i = 0; i < w; i++) {
        const int x = i < w / 2 ? i : i - w;
        const int y = j < h / 2 ? j : j - h;
        const int norm = x * x + y * y;
        direct::elem(out, i, j) = norm <= r2 ?
            direct::elem(img.data, i, j) : 0.0;
    }
    return {std::move(out)};
}

void FilterHelper::GaussianEnvelopeCorner2D(Image<RFLOAT> &img, double sigma) {
    const int w = img.data.xdim, h = img.data.ydim;
    const double s2 = sigma * sigma;
    for (int j = 0; j < h; j++)
    for (int i = 0; i < w; i++) {
        const double x = i < w / 2 ? i : i - w;
        const double y = j < h / 2 ? j : j - h;
        const double norm = x * x + y * y;
        direct::elem(img.data, i, j) *= exp(-0.5 * norm / s2);
    }
}

template <typename R>
Image<R> FilterHelper::raisedCosEnvCorner2D(const Image<R> &img, double radIn, double radOut) {
    const int w = img.data.xdim, h = img.data.ydim;
    MultidimArray<R> out (w, h);
    for (int j = 0; j < h; j++)
    for (int i = 0; i < w; i++) {
        const double x = i < w / 2 ? i : i - w;
        const double y = j < h / 2 ? j : j - h;
        const double r = std::hypot(x, y);
        if (r < radIn) {
            direct::elem(out, i, j) = direct::elem(img.data, i, j);
        } else if (r < radOut) {
            const double t = (r - radIn) / (radOut - radIn);
            const double a = raised_cos(PI * t);
            direct::elem(out, i, j) = a * direct::elem(img.data, i, j);
        } else {
            direct::elem(out, i, j) = 0.0;
        }
    }
    return {std::move(out)};
}

Image<RFLOAT> FilterHelper::raisedCosEnvCorner3D(Image<RFLOAT> &img, double radIn, double radOut) {
    const int w = img.data.xdim, h = img.data.ydim, d = img.data.zdim;
    MultidimArray<RFLOAT> out (w, h, d);
    for (int k = 0; k < d; k++)
    for (int j = 0; j < h; j++)
    for (int i = 0; i < w; i++) {
        const double x = i < w / 2 ? i : i - w;
        const double y = j < h / 2 ? j : j - h;
        const double z = k < d / 2 ? k : k - d;
        const double r = hypot(x, y, z);
        if (r < radIn) {
            direct::elem(out, i, j, k) = direct::elem(img.data, i, j, k);
        } else if (r < radOut) {
            const double t = (r - radIn) / (radOut - radIn);
            const double a = raised_cos(PI * t);
            direct::elem(out, i, j, k) = a * direct::elem(img.data, i, j, k);
        } else {
            direct::elem(out, i, j, k) = 0.0;
        }
    }
    return {std::move(out)};
}

Image<RFLOAT> FilterHelper::raisedCosEnvFreq2D(const Image<RFLOAT> &img, double radIn, double radOut) {
    const int w = img.data.xdim, h = img.data.ydim;
    MultidimArray<RFLOAT> out (w, h);
    for (int j = 0; j < h; j++)
    for (int i = 0; i < w; i++) {
        const double x = i;
        const double y = j <= h / 2 ? j : j - h;
        const double r = std::hypot(x, y);
        if (r < radIn) {
            direct::elem(out, i, j) = direct::elem(img.data, i, j);
        } else if (r < radOut) {
            const double t = (r - radIn) / (radOut - radIn);
            const double a = raised_cos(PI * t);
            direct::elem(out, i, j) = a * direct::elem(img.data, i, j);
        } else {
            direct::elem(out, i, j) = 0.0;
        }
    }
    return {std::move(out)};
}

Image<RFLOAT> FilterHelper::raisedCosEnvRingFreq2D(
    const Image<RFLOAT> &img,
    double rad0, double rad1, double stepWidth
) {
    const int w = img.data.xdim, h = img.data.ydim;
    MultidimArray<RFLOAT> out (w, h);
    for (int j = 0; j < h; j++)
    for (int i = 0; i < w; i++) {
        const double x = i;
        const double y = j <= h / 2 ? j : j - h;
        const double r = std::hypot(x, y);
        const double r0 = rad0 > 0.0 ? r - rad0 : stepWidth / 2;
        const double r1 = rad1 - r;
        double re = 2.0 * std::min(r0, r1) / stepWidth;
        if (re > 1.0) {
            direct::elem(out, i, j) = direct::elem(img.data, i, j);
        } else if (re > -1.0) {
            const double t = re * 0.5 + 0.5;
            const double a = raised_cos(PI * t);
            direct::elem(out, i, j) = a * direct::elem(img.data, i, j);
        } else {
            direct::elem(out, i, j) = 0.0;
        }
    }
    return {std::move(out)};
}

void FilterHelper::lowPassFilter(
    Image<RFLOAT> &img,
    double maxFreq0, double maxFreq1,
    Image<RFLOAT> &dest
) {
    FourierTransformer ft;
    MultidimArray<Complex> &imgFreq = ft.FourierTransform(img());
    lowPassFilterSpectrum(imgFreq, maxFreq0, maxFreq1);
    dest.data = ft.inverseFourierTransform(imgFreq);
}

void FilterHelper::lowPassFilterSpectrum(
    MultidimArray<Complex> &spectrum, double maxFreq0, double maxFreq1
) {
    for (long int k = 0; k < Zsize(spectrum); k++)
    for (long int j = 0; j < Ysize(spectrum); j++)
    for (long int i = 0; i < Xsize(spectrum); i++) {
        double x =       i / (double) spectrum.xdim;
        double y = 2.0 * j / (double) spectrum.ydim;
        double z = 2.0 * k / (double) spectrum.zdim;

        if (y > 1.0) y = 2.0 - y;
        if (z > 1.0) z = 2.0 - z;

        double r = hypot(x, y, z);

        if (r > maxFreq1) {
            direct::elem(spectrum, i, j, k) = 0.0;
        } else if (r > maxFreq0) {
            const double t = (r - maxFreq0) / (maxFreq1 - maxFreq0);
            const double q = raised_cos(PI * t);
            direct::elem(spectrum, i, j, k) *= q;
        }
    }
}

RFLOAT FilterHelper::averageValue(Image<RFLOAT> &img) {
    const RFLOAT sum = std::accumulate(img.data.begin(), img.data.end(), 0.0);
    return sum / img.data.size();
}

void FilterHelper::phaseFlip(
    Image<RFLOAT> &img,
    CTF &ctf, ObservationModel *obsModel, int opticsGroup,
    RFLOAT angpix, Image<RFLOAT> &dest
) {
    FourierTransformer ft;
    MultidimArray<Complex> &imgFreq = ft.FourierTransform(img());

    const RFLOAT xs = img.data.xdim * angpix;
    const RFLOAT ys = img.data.ydim * angpix;

    for (long int j = 0; j < Ysize(imgFreq); j++)
    for (long int i = 0; i < Xsize(imgFreq); i++) {
        RFLOAT x = j;
        RFLOAT y = i < imgFreq.ydim / 2 ? i : i - imgFreq.ydim;
        if (obsModel) obsModel->magnify(x, y, obsModel->getMagMatrix(opticsGroup));
        if (ctf(x / xs, y / ys) < 0)
            direct::elem(imgFreq, i, j) = -direct::elem(imgFreq, i, j);
    }

    // dest.data.resize(img.data.xdim, img.data.ydim);
    dest.data = ft.inverseFourierTransform(imgFreq);
}

void FilterHelper::applyBeamTilt(
    Image<RFLOAT> &img, RFLOAT beamtilt_x, RFLOAT beamtilt_y,
    RFLOAT lambda, RFLOAT Cs, RFLOAT angpix, int s, Image<RFLOAT> &dest
) {
    FourierTransformer ft;
    MultidimArray<Complex> &imgFreq = ft.FourierTransform(img.data);
    selfApplyBeamTilt(imgFreq, beamtilt_x, beamtilt_y, lambda, Cs, angpix, s);
    dest.data = ft.inverseFourierTransform(imgFreq);
}

void FilterHelper::modulate(
    Image<RFLOAT> &img,
    const CTF &ctf, ObservationModel *obsModel, int opticsGroup,
    RFLOAT angpix, Image<RFLOAT> &dest
) {
    FourierTransformer ft;
    // Can we pass a reference to ft.fFourier and prolong its lifetime?
    Image<Complex> imgFreq (ft.FourierTransform(img.data));
    modulate(imgFreq, ctf, obsModel, opticsGroup, angpix, dest);
}

void FilterHelper::modulate(
    Image<Complex> &imgFreq,
    const CTF &ctf, ObservationModel *obsModel, int opticsGroup,
    RFLOAT angpix, Image<RFLOAT> &dest
) {
    modulate(imgFreq.data, ctf, obsModel, opticsGroup, angpix);
    FourierTransformer ft;
    // dest.data.resize(2 * (w - 1), h);
    dest.data = ft.inverseFourierTransform(imgFreq.data);
}

void FilterHelper::modulate(
    MultidimArray<Complex> &imgFreq,
    const CTF& ctf, ObservationModel *obsModel, int opticsGroup,
    RFLOAT angpix
) {
    const int w = imgFreq.xdim, h = imgFreq.ydim;
    const MultidimArray<RFLOAT> ctfImg = CtfHelper::getFftwImage(ctf, w, h, h, h, angpix, obsModel, opticsGroup);
    imgFreq *= ctfImg;
}

void FilterHelper::drawCtf(
    CTF &ctf, ObservationModel *obsModel, int opticsGroup,
    RFLOAT angpix,
    Image<Complex> &dest
) {
    const int w = dest.data.xdim, h = dest.data.ydim;
    const MultidimArray<RFLOAT> ctfImg = CtfHelper::getFftwImage(ctf, w, h, h, h, angpix, obsModel, opticsGroup);
    dest.data = ctfImg;
}

void FilterHelper::wienerFilter(
    Image<RFLOAT> &img,
    CTF &ctf, ObservationModel *obsModel, int opticsGroup,
    RFLOAT angpix, RFLOAT eps, RFLOAT Bfac, Image<RFLOAT> &dest
) {
    FourierTransformer ft;
    MultidimArray<Complex> &imgFreq = ft.FourierTransform(img());

    RFLOAT xs = (RFLOAT) img.data.xdim * angpix;
    RFLOAT ys = (RFLOAT) img.data.ydim * angpix;

    for (long int j = 0; j < Ysize(imgFreq); j++)
    for (long int i = 0; i < Xsize(imgFreq); i++) {
        const int x = j;
        const int y = i < imgFreq.ydim / 2 ? i : i - imgFreq.ydim;

        RFLOAT x2 = x, y2 = y;
        if (obsModel) obsModel->magnify(x2, y2, obsModel->getMagMatrix(opticsGroup));
        RFLOAT c = ctf(x2 / xs, y2 / ys);
        if (Bfac > 0.0) { c *= exp(-Bfac * (x * x + y * y) / 4.0); }

        direct::elem(imgFreq, i, j) = c * direct::elem(imgFreq, i, j) / (c * c + eps);
    }

    // dest.data.resize(img.data);
    dest.data = ft.inverseFourierTransform(imgFreq);
}

void FilterHelper::richardsonLucy(
    Image<RFLOAT> &img, CTF &ctf, ObservationModel *obsModel, int opticsGroup, RFLOAT angpix, RFLOAT eps, int iterations, Image<RFLOAT> &dest
) {
    const int w = img.data.xdim, h = img.data.ydim;

    Image<RFLOAT> img1 (w, h, 1, 1), img1M (w, h, 1, 1), imgR (w, h, 1, 1), imgRM (w, h, 1, 1);

    double Bfac = w / 4.0;

    const double min = img.data.size() ?
        *std::min_element(img.data.begin(), img.data.end()) : 0;

    Image<RFLOAT> img0 (img.data + (min - 10));

    wienerFilter(img0, ctf, obsModel, opticsGroup, angpix, eps, Bfac, img1);

    VtkHelper::writeVTK(img1, "rl_it0.vtk");

    for (int it = 0; it < iterations; it++) {
        // img1 = img1 * conv(psf, img / conv(psf, img1) )
        //      = img1 * IFT( ctf * FT(img / IFT( ctf * FT(img1) ) ) )
        //      = img1 * ctf_mod( img / ctf_mod(img1) )

        modulate(img1, ctf, obsModel, opticsGroup, angpix, img1M);
        wienerDivide(img0, img1M, eps, imgR);
        modulate(imgR, ctf, obsModel, opticsGroup, angpix, imgRM);
        img1.data *= imgRM.data;

        VtkHelper::writeVTK(img1, "rl_it" + std::to_string(it + 1) + ".vtk");
    }
}

void FilterHelper::rampFilter(Image<RFLOAT> &img, RFLOAT s0, RFLOAT t1, double ux, double uy, Image<RFLOAT> &dest) {
    FourierTransformer ft;
    MultidimArray<Complex> &imgFreq = ft.FourierTransform(img.data);

    for (long int j = 0; j < Ysize(imgFreq); j++)
    for (long int i = 0; i < Xsize(imgFreq); i++) {
        const int x = j;
        const int y = i < imgFreq.ydim / 2 ? i : i - imgFreq.ydim;

        RFLOAT t = std::abs(x * ux + y * uy);
        RFLOAT s = t < t1 ? s0 + (1 - s0) * t / t1 : 1.0;

        direct::elem(imgFreq, i, j) *= s;
    }

    // dest.data.resize(img.data);
    dest.data = ft.inverseFourierTransform(imgFreq);
}

void FilterHelper::rampFilter3D(Image<Complex> &img, RFLOAT s0, RFLOAT t1, double tx, double ty, double tz) {
    const d3Vector ta (tx, ty, tz);
    for (long int k = 0; k < Zsize(img.data); k++)
    for (long int j = 0; j < Ysize(img.data); j++)
    for (long int i = 0; i < Xsize(img.data); i++) {
        const int x = j;
        const int y = i < img.data.ydim / 2 ? i : i - img.data.ydim;
        const int z = k < img.data.zdim / 2 ? k : k - img.data.zdim;

        d3Vector p (x, y, z);
        d3Vector q = p - p.dot(ta) * ta;
        const double t = q.length();

        RFLOAT s = t < t1 ? s0 + (1 - s0) * t / t1 : 1.0;

        direct::elem(img.data, i, j, k) *= s;
    }
}

void FilterHelper::doubleRampFilter3D(Image<Complex> &img, RFLOAT s0, RFLOAT t1, double tx, double ty, double tz) {
    const d3Vector ta (tx, ty, tz);
    for (long int k = 0; k < Zsize(img.data); k++)
    for (long int j = 0; j < Ysize(img.data); j++)
    for (long int i = 0; i < Xsize(img.data); i++) {
        const int x = j;
        const int y = i < img.data.ydim / 2 ? i : i - img.data.ydim;
        const int z = k < img.data.zdim / 2 ? k : k - img.data.zdim;

        d3Vector p (x, y, z);
        d3Vector q = p - p.dot(ta)*ta;
        const double t = q.length();

        RFLOAT s = t < t1 ? s0 + (1 - s0) * t / t1 : 1.0 + t1 - t;
        if (s < 0) s = 0;

        direct::elem(img.data, i, j, k) = s * direct::elem(img.data, i, j, k);
    }
}

Image<RFLOAT> FilterHelper::getPhase(const Image<Complex> &img) {
    const long w = img.data.xdim, h = img.data.ydim, d = img.data.zdim;
    MultidimArray<RFLOAT> dest (w, h, d);
    for (long int i = 0; i != dest.size(); i++)
        dest[i] = img.data[i].norm() > 0 ? img.data[i].arg() : 0;
    return {std::move(dest)};
}

Image<RFLOAT> FilterHelper::getAbs(const Image<Complex> &img) {
    const long w = img.data.xdim, h = img.data.ydim, d = img.data.zdim;
    MultidimArray<RFLOAT> dest (w, h, d);
    for (long int i = 0; i != dest.size(); i++)
        dest[i] = img.data[i].abs();
    return {std::move(dest)};
}

Image<RFLOAT> FilterHelper::getReal(const Image<Complex> &img) {
    const long w = img.data.xdim, h = img.data.ydim, d = img.data.zdim;
    MultidimArray<RFLOAT> dest (w, h, d);
    for (long int i = 0; i != dest.size(); i++) {
        dest[i] = img.data[i].real;
    }
    return {std::move(dest)};
}

Image<RFLOAT> FilterHelper::getImag(const Image<Complex> &img) {
    const long w = img.data.xdim, h = img.data.ydim, d = img.data.zdim;
    MultidimArray<RFLOAT> dest (w, h, d);
    for (long int i = 0; i != dest.size(); i++) {
        dest[i] = img.data[i].imag;
    }
    return {std::move(dest)};
}

void FilterHelper::powerSpectrum2D(Image<RFLOAT> &img, Volume<RFLOAT>& spectrum) {
    FourierTransformer ft;
    MultidimArray<Complex> &imgFreq = ft.FourierTransform(img.data);
    spectrum.resize(imgFreq.xdim, imgFreq.ydim, 1);
    for (long int j = 0; j < Ysize(imgFreq); j++)
    for (long int i = 0; i < Xsize(imgFreq); i++) {
        Complex z = direct::elem(imgFreq, i, j);
        spectrum(j, i, 0) = z.abs();
    }
}

void FilterHelper::equiphaseAverage2D(const Volume<RFLOAT> &src, Volume<RFLOAT> &dest) {

    const int n = src.dimx;
    std::vector<RFLOAT> val (n, 0.0), wgh (n, 0.0);

    for (long int y = 0; y < src.dimy; y++)
    for (long int x = 0; x < src.dimx; x++) {
        double id = y < src.dimy / 2 ?
            std::hypot(x, y) :
            std::hypot(x, src.dimy - y);

        int i = id;
        double f = id - i;

        if (i >= 0 && i < n) {
            val[i] += (1.0 - f) * src(x, y, 0);
            wgh[i] += (1.0 - f);
        }

        if (i >= -1 && i < n - 1) {
            val[i + 1] += f * src(x, y, 0);
            wgh[i + 1] += f;
        }
    }

    for (long int i = 0; i < n; i++)
        if (wgh[i] > 0.0)
            val[i] /= wgh[i];

    dest.resize(src);

    for (long int y = 0; y < src.dimy; y++)
    for (long int x = 0; x < src.dimx; x++) {

        const double r = y < src.dimy / 2 ?
            std::hypot(x, y) :
            std::hypot(x, src.dimy - y);
        const int ri = r;
        const double rf = r - ri;

        if (ri >= 0 && ri < n - 1) {
            dest(x, y, 0) = (1.0 - rf) * val[ri] + rf * val[ri + 1];
        }
    }
}

void FilterHelper::threshold(Image<RFLOAT> &src, RFLOAT t, Image<RFLOAT> &dest) {
    for (long int i = 0; i != src.data.size(); i++)
        dest.data[i] = (RFLOAT) (src.data[i] > t);
}
void FilterHelper::linearTransform(
    Image<RFLOAT> &src, RFLOAT m, RFLOAT q, Image<RFLOAT> &dest
) {
    for (long int i = 0; i != src.data.size(); i++)
        dest.data[i] = m * src.data[i] + q;
}

void FilterHelper::sumUp(
    const std::vector<Image<RFLOAT>> &srcs, Image<RFLOAT> &dest
) {
    const int w = srcs[0].data.xdim;
    const int h = srcs[0].data.ydim;
    const int d = srcs[0].data.zdim;
    const int m = srcs[0].data.ndim;
    dest.data.resize(w, h, d, m);
    for (const auto& src: srcs) {

        if (
            src.data.xdim != w || src.data.ydim != h ||
            src.data.zdim != d || src.data.ndim != m
        ) REPORT_ERROR("FilterHelper::sumUp(): image dimension mismatch.\n");

        dest.data += src.data;
    }
}

double FilterHelper::L1distance(
    const Image<RFLOAT> &i0, const Image<RFLOAT> &i1,
    int x0, int y0, int w, int h
) {

    if (w < 0) w = i0.data.xdim;
    if (h < 0) h = i0.data.ydim;

    double d = 0.0;
    for (long int n = 0; n < i0.data.ndim; n++)
    for (long int z = 0; z < i0.data.zdim; z++)
    for (long int y = y0; y < y0 + h; y++)
    for (long int x = x0; x < x0 + w; x++) {
        const RFLOAT a = direct::elem(i0.data, x, y, z, n);
        const RFLOAT b = direct::elem(i1.data, x, y, z, n);
        d += std::abs(b - a);
    }
    return d;
}

double FilterHelper::L2distance(
    const Image<RFLOAT> &i0, const Image<RFLOAT> &i1, int x0, int y0, int w, int h
) {

    if (w < 0) w = i0.data.xdim;
    if (h < 0) h = i0.data.ydim;

    double d = 0.0;
    for (long int n = 0; n < i0.data.ndim; n++)
    for (long int z = 0; z < i0.data.zdim; z++)
    for (long int y = y0; y < y0 + h; y++)
    for (long int x = x0; x < x0 + w; x++) {
        const RFLOAT a = direct::elem(i0.data, x, y, z, n);
        const RFLOAT b = direct::elem(i1.data, x, y, z, n);
        d += (b - a) * (b - a);
    }
    return d;
}

double FilterHelper::NCC(
    const Image<RFLOAT> &i0, const Image<RFLOAT> &i1,
    int x0, int y0, int w, int h
) {
    double d = 0.0;

    if (w < 0) { w = i0.data.xdim; }
    if (h < 0) { h = i0.data.ydim; }

    double mu0 = 0.0, mu1 = 0.0;
    const int cnt = i0.data.size();
    for (long int n = 0; n < i0.data.ndim; n++)
    for (long int z = 0; z < i0.data.zdim; z++)
    for (long int y = y0; y < y0 + h; y++)
    for (long int x = x0; x < x0 + w; x++) {
        const RFLOAT v0 = direct::elem(i0.data, x, y, z, n);
        const RFLOAT v1 = direct::elem(i1.data, x, y, z, n);
        mu0 += v0;
        mu1 += v1;
    }
    mu0 /= cnt;
    mu1 /= cnt;

    double sig0 = 0.0, sig1 = 0.0;
    for (long int n = 0; n < i0.data.ndim; n++)
    for (long int z = 0; z < i0.data.zdim; z++)
    for (long int y = y0; y < y0 + h; y++)
    for (long int x = x0; x < x0 + w; x++) {
        const RFLOAT v0 = direct::elem(i0.data, x, y, z, n) - mu0;
        const RFLOAT v1 = direct::elem(i1.data, x, y, z, n) - mu1;
        sig0 += v0 * v0;
        sig1 += v1 * v1;
    }
    sig0 = sqrt(sig0 / (cnt - 1));
    sig1 = sqrt(sig1 / (cnt - 1));

    double ncc = 0.0;
    for (long int n = 0; n < i0.data.ndim; n++)
    for (long int z = 0; z < i0.data.zdim; z++)
    for (long int y = y0; y < y0 + h; y++)
    for (long int x = x0; x < x0 + w; x++) {
        const RFLOAT v0 = (direct::elem(i0.data, x, y, z, n) - mu0);
        const RFLOAT v1 = (direct::elem(i1.data, x, y, z, n) - mu1);
        ncc += v0 * v1;
    }

    return ncc / (sig0 * sig1 * cnt);
}

void FilterHelper::wienerDivide(Image<RFLOAT> &num, Image<RFLOAT> &denom, RFLOAT eps, Image<RFLOAT> &dest) {
    dest.data.resize(num.data);
    for (long int i = 0; i < num.data.size(); i++) {
        const RFLOAT d = denom.data[i];
        dest.data[i] = d * num.data[i] / (d * d + eps);
    }
}

void FilterHelper::divideExcessive(
    Image<Complex> &dividend, Volume<RFLOAT> &divisor, RFLOAT theta,
    Image<Complex> &dest
) {
    for (long int n = 0; n < dividend.data.ndim; n++)
    for (long int z = 0; z < dividend.data.zdim; z++)
    for (long int y = 0; y < dividend.data.ydim; y++)
    for (long int x = 0; x < dividend.data.xdim; x++) {
        const RFLOAT t = divisor(x, y, z) / theta;
        direct::elem(dest.data, x, y, z, n) = t > 1 ?
            direct::elem(dividend.data, x, y, z, n) / t :
            direct::elem(dividend.data, x, y, z, n);
    }
}

void FilterHelper::wienerDeconvolve(
    Image<Complex> &dividend, Image<Complex> &divisor, RFLOAT theta,
    Image<Complex> &dest
) {
    for (long int z = 0; z < dividend.data.zdim; z++)
    for (long int y = 0; y < dividend.data.ydim; y++)
    for (long int x = 0; x < dividend.data.xdim; x++) {
        Complex zz = direct::elem(divisor.data, x, y, z);
        Complex z0 = direct::elem(dividend.data, x, y, z);

        /*
        std::cout << "z0 = " << z0.real << " + " << z0.imag << " * i\n";
        std::cout << "zz = " << zz.real << " + " << zz.imag << " * i\n";
        std::cout << "zzB * z0 = " << (zz.conj() * z0).real << " + " << (zz.conj() * z0).imag << " * i\n";
        std::cout << "((zz.conj() * zz).real + theta) = " << ((zz.conj() * zz).real + theta) << " * i\n";
        */

        // direct::elem(dest.data, x, y, z) = (zz.conj() * z0) / ((zz.conj() * zz).real + theta);

        direct::elem(dest.data, x, y, z) = (zz.real * z0) / (zz.real * zz.real + theta);

        /*
        RFLOAT t = zz.abs() / theta;

        if (t > 1) {
            direct::elem(dest.data, x, y, z) = direct::elem(dividend.data, x, y, z) / t;
        } else {
            direct::elem(dest.data, x, y, z) = direct::elem(dividend.data, x, y, z);
        }
        */
    }
}

void FilterHelper::extract2D(
    const Image<RFLOAT> &src, Image<RFLOAT> &dest,
    long int x0, long int y0,
    long int w, long int h
) {
    dest.data.resize(w, h);
    for (long int y = 0; y < h; y++)
    for (long int x = 0; x < w; x++) {
        const long int xx = x0 + x;
        const long int yy = y0 + y;
        direct::elem(dest.data, x, y) =
            xx >= 0 && xx < src.data.xdim &&
            yy >= 0 && yy < src.data.ydim ?
            direct::elem(src.data, xx, yy) : 0;
    }
}

void FilterHelper::extract(
    const Volume<RFLOAT> &src,
    Volume<RFLOAT> &dest,
    long int x0, long int y0, long int z0,
    long int w, long int h, long int d
) {
    dest.resize(w, h, d);

    for (long int z = 0; z < d; z++)
    for (long int y = 0; y < h; y++)
    for (long int x = 0; x < w; x++) {
        long int xx = x0 + x;
        long int yy = y0 + y;
        long int zz = z0 + z;

        if (
            xx >= 0 && xx < src.dimx &&
            yy >= 0 && yy < src.dimy &&
            zz >= 0 && zz < src.dimz
        ) {
            dest(x, y, z) = src(xx, yy, zz);
        }
    }
}

void FilterHelper::signedDist(const Image<RFLOAT> &src, Image<RFLOAT> &dest) {
    dest = Image<RFLOAT>(src.data.xdim, src.data.ydim, src.data.zdim);

    Image<RFLOAT>
    ggp(src.data.xdim, src.data.ydim, src.data.zdim),
    ggn(src.data.xdim, src.data.ydim, src.data.zdim),
    gp(src.data.xdim,  src.data.ydim, src.data.zdim),
    gn(src.data.xdim,  src.data.ydim, src.data.zdim),
    hp(src.data.xdim,  src.data.ydim, src.data.zdim),
    hn(src.data.xdim,  src.data.ydim, src.data.zdim),
    s(src.data.xdim,   src.data.ydim, src.data.zdim);

    double rmax2 = 4.0 * (src.data.xdim * src.data.xdim + src.data.ydim * src.data.ydim + src.data.zdim * src.data.zdim);

    for (long int z = 0; z < dest.data.zdim; z++)
    for (long int y = 0; y < dest.data.ydim; y++) {
        direct::elem(ggp.data, 0, y, z) = rmax2;
        direct::elem(ggn.data, 0, y, z) = rmax2;

        for (long int x = 1; x < dest.data.xdim; x++) {
            if (direct::elem(src.data, x, y, z) < 0.0) {
                direct::elem(ggp.data, x, y, z) = 0;

                double d = sqrt(direct::elem(ggn.data, x - 1, y, z)) + 1.0;
                direct::elem(ggn.data, x, y, z) = d * d;
            } else {
                direct::elem(ggn.data, x, y, z) = 0;

                double d = sqrt(direct::elem(ggp.data, x - 1, y, z)) + 1.0;
                direct::elem(ggp.data, x, y, z) = d * d;
            }
        }

        direct::elem(gp.data, dest.data.xdim - 1, y, z) = direct::elem(ggp.data, dest.data.xdim - 1, y, z);
        direct::elem(gn.data, dest.data.xdim - 1, y, z) = direct::elem(ggn.data, dest.data.xdim - 1, y, z);

        for (long int x = dest.data.xdim - 2; x >= 0; x--) {
            double dp = sqrt(direct::elem(gp.data, x + 1, y, z)) + 1.0;
            double ddp = dp * dp;

            double dn = sqrt(direct::elem(gn.data, x + 1, y, z)) + 1.0;
            double ddn = dn * dn;

            if (ddp < direct::elem(ggp.data, x, y, z)) {
                direct::elem(gp.data, x, y, z) = ddp;
            } else {
                direct::elem(gp.data, x, y, z) = direct::elem(ggp.data, x, y, z);
            }

            if (ddn < direct::elem(ggn.data, x, y, z)) {
                direct::elem(gn.data, x, y, z) = ddn;
            } else {
                direct::elem(gn.data, x, y, z) = direct::elem(ggn.data, x, y, z);
            }
        }
    }

    for (long int z = 0; z < dest.data.zdim; z++)
    for (long int y = 0; y < dest.data.ydim; y++)
    for (long int x = 0; x < dest.data.xdim; x++) {
        long int rp = sqrt(direct::elem(gp.data, x, y, z));
        long int rn = sqrt(direct::elem(gn.data, x, y, z));

        double minValP = rmax2;
        double minValN = rmax2;

        for (long int yy = y - rp; yy <= y + rp; yy++) {
            if (yy < 0 || yy >= dest.data.ydim) continue;

            double dy = yy - y;
            double vgp = direct::elem(gp.data, x, yy, z) + dy * dy;

            if (vgp < minValP) { minValP = vgp; }
        }

        for (long int yy = y - rn; yy <= y + rn; yy++) {
            if (yy < 0 || yy >= dest.data.ydim) continue;

            double dy = yy - y;
            double vgn = direct::elem(gn.data, x, yy, z) + dy * dy;

            if (vgn < minValN) { minValN = vgn; }
        }

        direct::elem(hp.data, x, y, z) = minValP;
        direct::elem(hn.data, x, y, z) = minValN;
    }

    for (long int z = 0; z < dest.data.zdim; z++)
    for (long int y = 0; y < dest.data.ydim; y++)
    for (long int x = 0; x < dest.data.xdim; x++) {
        if (direct::elem(src.data, x, y, z) < 0.0) {
            direct::elem(dest.data, x, y, z) = -sqrt(direct::elem(hn.data, x, y, z));
        } else {
            direct::elem(dest.data, x, y, z) = +sqrt(direct::elem(hp.data, x, y, z));
        }
    }
}

void FilterHelper::erode3x3(Image<RFLOAT> &src, Image<RFLOAT> &dest) {
    for (long int z = 0; z < src.data.zdim; z++)
    for (long int y = 0; y < src.data.ydim; y++)
    for (long int x = 0; x < src.data.xdim; x++) {
        double v = std::numeric_limits<double>::max();

        for (long int zz = z - 1; zz <= z + 1; zz++)
        for (long int yy = y - 1; yy <= y + 1; yy++)
        for (long int xx = x - 1; xx <= x + 1; xx++) {
            if (xx >= 0 && xx < src.data.xdim &&
                zz >= 0 && zz < src.data.zdim &&
                yy >= 0 && yy < src.data.ydim)
                v = std::min(direct::elem(src.data, xx, yy, zz), v);
        }

        direct::elem(dest.data, x, y, z) = v;
    }
}

void FilterHelper::localMinima(Image<RFLOAT> &src, RFLOAT thresh, Image<RFLOAT> &dest) {
    dest.data.resize(src.data.xdim, src.data.ydim, src.data.zdim);
    for (long int k = 0; k < src.data.zdim; k++)
    for (long int j = 0; j < src.data.ydim; j++)
    for (long int i = 0; i < src.data.xdim; i++) {
        if (direct::elem(src.data, i, j, k) > thresh) {
            direct::elem(dest.data, i, j, k) = 0.0;
            continue;
        }

        double v = std::numeric_limits<double>::max();
        for (long int kk = k - 1; kk <= k + 1; kk++)
        for (long int jj = j - 1; jj <= j + 1; jj++)
        for (long int ii = i - 1; ii <= i + 1; ii++) {
            if (ii >= 0 && ii < src.data.xdim &&
                jj >= 0 && jj < src.data.ydim &&
                kk >= 0 && kk < src.data.zdim)
                v = std::min(direct::elem(src.data, ii, jj, kk), v);
        }

        direct::elem(dest.data, i, j, k) = direct::elem(src.data, i, j, k) == v;
    }
}

std::vector<gravis::d3Vector> FilterHelper::localMinima(
    Image<RFLOAT> &src, RFLOAT thresh
) {
    std::vector<gravis::d3Vector> indices;
    for (long int k = 0; k < src.data.zdim; k++)
    for (long int j = 0; j < src.data.ydim; j++)
    for (long int i = 0; i < src.data.xdim; i++) {

        if (direct::elem(src.data, i, j, k) > thresh) continue;

        double v = std::numeric_limits<double>::max();
        for (long int kk = k - 1; kk <= k + 1; kk++)
        for (long int jj = j - 1; jj <= j + 1; jj++)
        for (long int ii = i - 1; ii <= i + 1; ii++) {
            if (ii >= 0 && ii < src.data.xdim &&
                jj >= 0 && jj < src.data.ydim &&
                kk >= 0 && kk < src.data.zdim)
                v = std::min(direct::elem(src.data, ii, jj, kk), v);
        }

        if (v == direct::elem(src.data, i, j, k))
            indices.emplace_back(i, j, k);
    }
    return indices;
}

void FilterHelper::centralGradient(
    const Volume<RFLOAT> &src, Volume<t3Vector<RFLOAT>> &dest
) {
    const size_t dimx = src.dimx;
    const size_t dimy = src.dimy;
    const size_t dimz = src.dimz;

    dest.resize(dimx, dimy, dimz);

    FOR_ALL_VOXELS(src) {

        // DRY !!!

        if (dimx == 0) {
            dest(x, y, z).x = 0;
        } else if (x == 0) {
            dest(x, y, z).x = src(x + 1, y, z) - src(x, y, z);
        } else if (x < dimx - 1) {
            dest(x, y, z).x = 0.5 * (src(x + 1, y, z) - src(x - 1, y, z));
        } else {
            dest(x, y, z).x = src(x, y, z) - src(x - 1, y, z);
        }

        if (dimy == 0) {
            dest(x, y, z).y = 0;
        } else if (y == 0) {
            dest(x, y, z).y = src(x, y + 1, z) - src(x, y, z);
        } else if (y < dimy - 1) {
            dest(x, y, z).y = 0.5 * (src(x, y + 1, z) - src(x, y - 1, z));
        } else {
            dest(x, y, z).y = src(x, y, z) - src(x, y - 1, z);
        }

        if (dimz == 0) {
            dest(x, y, z).z = 0;
        } else if (z == 0) {
            dest(x, y, z).z = src(x, y, z + 1) - src(x, y, z);
        } else if (z < dimz - 1) {
            dest(x, y, z).z = 0.5 * (src(x, y, z + 1) - src(x, y, z - 1));
        } else {
            dest(x, y, z).z = src(x, y, z) - src(x, y, z - 1);
        }

    }
}

t3Vector<RFLOAT> FilterHelper::centralGradient(const Volume<RFLOAT> &src, size_t x, size_t y, size_t z) {
    t3Vector<RFLOAT> out;

    if (src.dimx == 0) {
        out.x = 0;
    } else if (x == 0) {
        out.x = src(x + 1, y, z) - src(x, y, z);
    } else if (x < src.dimx - 1) {
        out.x = 0.5 * (src(x + 1, y, z) - src(x - 1, y, z));
    } else {
        out.x = src(x, y, z) - src(x - 1, y, z);
    }

    if (src.dimy == 0) {
        out.y = 0;
    } else if (y == 0) {
        out.y = src(x, y + 1, z) - src(x, y, z);
    } else if (y < src.dimy - 1) {
        out.y = 0.5 * (src(x, y + 1, z) - src(x, y - 1, z));
    } else {
        out.y = src(x, y, z) - src(x, y - 1, z);
    }

    if (src.dimz == 0) {
        out.z = 0;
    } else if (z == 0) {
        out.z = src(x, y, z + 1) - src(x, y, z);
    } else if (z < src.dimz - 1) {
        out.z = 0.5 * (src(x, y, z + 1) - src(x, y, z - 1));
    } else {
        out.z = src(x, y, z) - src(x, y, z - 1);
    }

    return out;
}

MultidimArray<Complex> FilterHelper::FriedelExpand(const MultidimArray<Complex> &half) {
    const int wh = half.xdim;
    const int h  = half.ydim;
    const int d  = half.zdim;
    const int c  = half.ndim;
    const int w = 2 * (wh - 1);

    MultidimArray<Complex> out (w, h, d, c);
    for (int l = 0; l < c; l++)
    for (int k = 0; k < d; k++)
    for (int j = 0; j < h; j++) {
        const int kk = (d - k) % d;
        const int jj = (h - j) % h;
        for (int i = 0; i < wh; i++)
            direct::elem(out, i, j, k, l) = direct::elem(half, i, j, k, l);
        for (int i = wh; i < w; i++)
            direct::elem(out, i, j, k, l) = direct::elem(half, w - i, jj, kk, l).conj();
    }
    return out;
}

Image<RFLOAT> FilterHelper::normaliseToUnitInterval(const Image<RFLOAT> &img) {
    RFLOAT minimum, maximum;
    if (img.data.size()) {
        const auto minmax = std::minmax_element(img.data.begin(), img.data.end());
        minimum = *minmax.first;
        maximum = *minmax.second;
    }
    /// TODO: Expression templates
    MultidimArray<RFLOAT> out (img.data.xdim, img.data.ydim, img.data.zdim, img.data.ndim);
    for (long int i = 0; i != img.data.size(); i++)
        out[i] = (img.data[i] - minimum) / (maximum - minimum);
    return {std::move(out)};
}

Image<RFLOAT> FilterHelper::normaliseToUnitIntervalSigned(const Image<RFLOAT> &img) {
    RFLOAT max_abs = 0;
    for (const RFLOAT& x: img.data)
        max_abs = std::max(std::abs(x), max_abs);
    return {img.data / max_abs};
}

void FilterHelper::uniqueInfluenceMask(std::vector<gravis::d3Vector> pts, Image<RFLOAT> &dest, Image<RFLOAT> &indexDest, RFLOAT thresh) {
    const long int w = dest.data.xdim, h = dest.data.ydim;
    const long int pc = pts.size();

    indexDest = Image<RFLOAT> (w, h);

    const double t2 = thresh * thresh;

    for (long int y = 0; y < h; y++)
    for (long int x = 0; x < w; x++) {

        int count_closer = 0;
        int last_i = -1;
        for (long int p = 0; p < pc; p++) {
            const d2Vector d (x - pts[p].x, y - pts[p].y);
            if (d.norm2() < t2) {
                count_closer++;
                last_i = p;
            }
        }

        if (count_closer == 1) {
            direct::elem(dest.data,      x, y) = 1.0;
            direct::elem(indexDest.data, x, y) = last_i;
        } else {
            direct::elem(dest.data,      x, y) = 0.0;
            direct::elem(indexDest.data, x, y) = last_i;
        }
    }
}

void FilterHelper::polarRemap(
    d2Vector pos, const Image<RFLOAT> &src, Image<RFLOAT> &dest,
    const Image<RFLOAT> &mask, Image<RFLOAT> &maskDest, int phiRes, int rRes, double rMax
) {
    const long int w = src.data.xdim, h = src.data.ydim;

    dest    .data.resize(phiRes, rRes);
    maskDest.data.resize(phiRes, rRes);

    for (long int ri = 0; ri < rRes; ri++)
    for (long int p = 0; p < phiRes; p++) {
        const double r   = rMax * ri    / (double) rRes;
        const double phi = 2.0 * PI * p / (double) phiRes;

        // This could be sped up with SINCOS
        d2Vector pp = pos + r * d2Vector(cos(phi), sin(phi));

        int ppnnx = pp.x + 0.5;
        int ppnny = pp.y + 0.5;

        if (
            ppnnx > 0 && ppnnx < w - 1 &&
            ppnny > 0 && ppnny < h - 1 &&
            direct::elem(mask.data, ppnnx, ppnny) > 0.5
        ) {
            direct::elem(dest.data,     p, ri) = Interpolation::linearXY(src, pp.x, pp.y, 0);
            direct::elem(maskDest.data, p, ri) = 1.0;
        } else {
            direct::elem(dest.data,     p, ri) = 0.0;
            direct::elem(maskDest.data, p, ri) = 0.0;
        }
    }
}

void FilterHelper::polarRemap(
    d2Vector pos, const Image<RFLOAT> &distTransf, const Image<RFLOAT> &src, Image<RFLOAT> &dest,
    const Image<RFLOAT> &mask, Image<RFLOAT> &maskDest, int phiRes, int rRes, double rMax
) {

    dest    .data.resize(phiRes, rRes);
    maskDest.data.resize(phiRes, rRes);

    for (long int i = 0; i != phiRes * rRes; i++) {
        dest    .data[i] = 0.0;
        maskDest.data[i] = 0.0;
    }

    const int x0 = pos.x - rMax + 0.5;
    const int x1 = pos.x + rMax + 0.5;
    const int y0 = pos.y - rMax + 0.5;
    const int y1 = pos.y + rMax + 0.5;

    const long int w = src.data.xdim, h = src.data.ydim;
    for (int j = y0; j <= y1; j++)
    for (int i = x0; i <= x1; i++) {
        const double dx = i - pos.x;
        const double dy = j - pos.y;

        if (i < 1 || i >= w - 1 || j < 1 || j >= h - 1 || dx == 0.0 && dy == 0.0)
            continue;

        double phiR = std::atan2(dy, dx);

        if (phiR < 0.0) phiR += 2.0 * PI;

        const double phiD = phiRes * phiR / (2.0 * PI);
        const int phi0 =  (int) phiD      % phiRes;
        const int phi1 = ((int) phiD + 1) % phiRes;
        const double phiF = phiD - (double) phi0;

        const double rD = rRes * direct::elem(distTransf.data, i, j) / rMax;
        const int r0 = rD;
        const int r1 = rD + 1;
        const double rF = rD - r0;

        const double v = direct::elem(src.data, i, j);

        if (r0 >= 0 && r0 < rRes) {
            direct::elem(dest.data,     phi0, r0) += (1.0 - rF) * (1.0 - phiF) * v;
            direct::elem(maskDest.data, phi0, r0) += (1.0 - rF) * (1.0 - phiF);

            direct::elem(dest.data,     phi1, r0) += (1.0 - rF) * phiF * v;
            direct::elem(maskDest.data, phi1, r0) += (1.0 - rF) * phiF;
        }

        if (r1 >= 0 && r1 < rRes) {
            direct::elem(dest.data,     phi0, r1) += rF * (1.0 - phiF) * v;
            direct::elem(maskDest.data, phi0, r1) += rF * (1.0 - phiF);

            direct::elem(dest.data,     phi1, r1) += rF * phiF * v;
            direct::elem(maskDest.data, phi1, r1) += rF * phiF;
        }

    }

    for (long int i = 0; i != rRes * phiRes; i++)
        if (maskDest.data[i] > 0.0)
            dest.data[i] /= maskDest.data[i];
}

Image<RFLOAT> FilterHelper::cartToPolar(const Image<RFLOAT> &img) {
    const int w = img.data.xdim;
    const int h = img.data.ydim;

    const double half_w = w / 2.0;
    const double half_h = h / 2.0;

    const double cx = half_w + 1.0;
    const double cy = half_h + 0.5;

    const double two_pi = 2.0 * PI;

    MultidimArray<RFLOAT> out (two_pi * half_w + 1, half_w);
    for (int j = 0; j < h; j++)
    for (int i = 0; i < w; i++) {
        const double r   = half_w * j / h;
        const double phi = two_pi * i / w;
        const double x = cx + r * cos(phi);
        const double y = cy + r * sin(phi);
        direct::elem(out, i, j) = Interpolation::cubicXY(img, x, y, 0, 0);
    }
    return {std::move(out)};
}

Image<RFLOAT> FilterHelper::polarToCart(const Image<RFLOAT> &img) {
    const int wp = img.data.xdim;
    const int hp = img.data.ydim;

    const double w0h = hp;

    const double cx = w0h + 1;
    const double cy = w0h + 0.5;

    const int w = 2.0 * w0h;
    const int h = 2.0 * w0h;

    MultidimArray<RFLOAT> out (w, h);
    const double two_pi = 2.0 * PI;
    for (int y = 0; y < h; y++)
    for (int x = 0; x < w; x++) {
        const double xd = x - cx;
        const double yd = y - cy;
        const double r = std::hypot(xd, yd);
        double phi = atan2(yd, xd);
        if (phi < 0.0) phi += two_pi;
        direct::elem(out, x, y) = Interpolation::cubicXY(img, wp * phi / two_pi, r, 0, 0);
    }
    return {std::move(out)};
}

Image<RFLOAT> FilterHelper::polarBlur(const Image<RFLOAT> &img, double sigma) {
    const Image<RFLOAT> polar = FilterHelper::cartToPolar(img);
    const Image<RFLOAT> blurred = separableGaussianX(polar, sigma, 2.0 * sigma + 0.5, true);
    return FilterHelper::polarToCart(blurred);
}

Image<RFLOAT> FilterHelper::sectorBlend(const Image<RFLOAT> &img0, const Image<RFLOAT> &img1, int sectors) {

    const int w = img0.data.xdim, h = img0.data.ydim;
    if (img0.data.sameShape(img1.data)) {
        std::cerr << "FilterHelper::sectorBlend: unequal image size: " << w << "x" << h
                  << " vs. " << img1.data.xdim << "x" << img1.data.ydim << "\n";
        REPORT_ERROR("FilterHelper::sectorBlend: unequal image size.");
    }

    MultidimArray<RFLOAT> out (w, h);
    const double cx = w / 2.0;
    const double cy = h / 2.0;
    for (int j = 0; j < h; j++)
    for (int i = 0; i < w; i++) {
        const double xd = i - cx;
        const double yd = j - cy;
        const double phi = atan2(yd, xd) + PI;
        const double a = sectors * phi / (2.0 * PI);
        direct::elem(out, i, j) = a - (int) a < 0.5 ? img0(i, j) : img1(i, j);
    }
    return {std::move(out)};
}

void FilterHelper::diffuseAlongIsocontours2D(
    const Image<RFLOAT> &src, const Image<RFLOAT> &guide,
    Image<RFLOAT> &dest, int iters, RFLOAT sigma, RFLOAT lambda, RFLOAT delta
) {
    const long int w = src.data.xdim, h = src.data.ydim;

    const bool sobel = true;  // ?

    dest.data.resize(w, h);
    std::copy_n(src.data.begin(), w * h, dest.data.begin());

    Volume<d2Vector> flux (w, h, 1);
    flux.fill(d2Vector(0, 0));

    using tensor_t = Tensor2x2<RFLOAT>;
    Volume<tensor_t> D0 (w, h, 1), D (w, h, 1), J (w, h, 1);
    D0.fill({0.0});
    D .fill({0.0});
    J .fill({0.0});

    for (long int y = 1; y < h - 1; y++)
    for (long int x = 1; x < w - 1; x++) {
        d2Vector g;

        if (sobel) {
            const double gxp = 0.25 * direct::elem(guide.data, x + 1, y - 1)
                             + 0.5  * direct::elem(guide.data, x + 1, y    )
                             + 0.25 * direct::elem(guide.data, x + 1, y + 1);
            const double gxn = 0.25 * direct::elem(guide.data, x - 1, y - 1)
                             + 0.5  * direct::elem(guide.data, x - 1, y    )
                             + 0.25 * direct::elem(guide.data, x - 1, y + 1);
            const double gyp = 0.25 * direct::elem(guide.data, x - 1, y + 1)
                             + 0.5  * direct::elem(guide.data, x,     y + 1)
                             + 0.25 * direct::elem(guide.data, x + 1, y + 1);
            const double gyn = 0.25 * direct::elem(guide.data, x - 1, y - 1)
                             + 0.5  * direct::elem(guide.data, x,     y - 1)
                             + 0.25 * direct::elem(guide.data, x + 1, y - 1);
            g.set(0.5 * (gxp - gxn), 0.5 * (gyp - gyn));
        } else {
            const double gxp = direct::elem(guide.data, x + 1, y);
            const double gxn = direct::elem(guide.data, x - 1, y);
            const double gyp = direct::elem(guide.data, x, y + 1);
            const double gyn = direct::elem(guide.data, x, y - 1);
            g.set(0.5 * (gxp - gxn), 0.5 * (gyp - gyn));
        }

        D0(x, y, 0) = tensor_t::dyadicProduct(g, g);
    }

    separableGaussian(D0, D, sigma);

    // Volume<RFLOAT> dbg0 (w, h, 1), dbg1 (w, h, 1);

    for (long int y = 0; y < h; y++)
    for (long int x = 0; x < w; x++) {
        t2Matrix<RFLOAT> DxyR = D(x, y, 0).toMatrix();
        d2Matrix Dxy (DxyR(0, 0), DxyR(0, 1), DxyR(1, 0), DxyR(1, 1));

        double qx, qy, l0, l1;
        dsyev2(Dxy(0, 0), Dxy(0, 1), Dxy(1, 1), &l0, &l1, &qx, &qy);

        const double dl = l0 - l1;
        const RFLOAT ani = 1.0 - exp(-0.5 * dl * dl / (lambda * lambda));

        const d2Vector f (-qy, qx);

        // dbg0(x, y, 0) = f.length();

        J(x, y, 0) = ani * tensor_t::dyadicProduct(f, f);
    }

    // VtkHelper::writeVTK(dbg0, "f_len.vtk");

    for (int it = 0; it < iters; it++) {

        #if JAZ_USE_OPENMP
        #pragma omp parallel for
        #endif
        for (long int y = 1; y < h - 1; y++)
        for (long int x = 1; x < w - 1; x++) {

            d2Vector g;
            if (sobel) {
                const double gxp = 0.25 * direct::elem(dest.data, x + 1, y - 1)
                                 + 0.5  * direct::elem(dest.data, x + 1, y    )
                                 + 0.25 * direct::elem(dest.data, x + 1, y + 1);
                const double gxn = 0.25 * direct::elem(dest.data, x - 1, y - 1)
                                 + 0.5  * direct::elem(dest.data, x - 1, y    )
                                 + 0.25 * direct::elem(dest.data, x - 1, y + 1);
                const double gyp = 0.25 * direct::elem(dest.data, x - 1, y + 1)
                                 + 0.5  * direct::elem(dest.data, x,     y + 1)
                                 + 0.25 * direct::elem(dest.data, x + 1, y + 1);
                const double gyn = 0.25 * direct::elem(dest.data, x - 1, y - 1)
                                 + 0.5  * direct::elem(dest.data, x,     y - 1)
                                 + 0.25 * direct::elem(dest.data, x + 1, y - 1);
                g.set(0.5 * (gxp - gxn), 0.5 * (gyp - gyn));
            } else {
                const double gxp = direct::elem(guide.data, x + 1, y);
                const double gxn = direct::elem(guide.data, x - 1, y);
                const double gyp = direct::elem(guide.data, x, y + 1);
                const double gyn = direct::elem(guide.data, x, y - 1);
                g.set(0.5 * (gxp - gxn), 0.5 * (gyp - gyn));
            }

            flux(x, y, 0) = J(x, y, 0).toMatrix() * g;
        }

        #if JAZ_USE_OPENMP
        #pragma omp parallel for
        #endif
        for (long int y = 1; y < h - 1; y++)
        for (long int x = 1; x < w - 1; x++) {
            const double div = flux(x + 1, y, 0).x - flux(x - 1, y, 0).x;
                             + flux(x, y + 1, 0).y - flux(x, y - 1, 0).y;
            direct::elem(dest.data, x, y) += delta * div;
        }
    }
}

void FilterHelper::EED_2D(const Image<RFLOAT> &src, Image<RFLOAT> &dest, int iters, double sigma, double delta, double tau) {

    const long int w = src.data.xdim, h = src.data.ydim;
    std::copy_n(src.data.begin(), w * h, dest.data.begin());

    const Image<RFLOAT> smooth = separableGaussianXY(dest, sigma);

    Volume<d2Vector> flux (w, h, 1);
    flux.fill(d2Vector(0, 0));

    const double tt = tau * tau;

    for (int it = 0; it < iters; it++) {

        #if JAZ_USE_OPENMP
        #pragma omp parallel for
        #endif
        for (long int y = 1; y < h - 1; y++)
        for (long int x = 1; x < w - 1; x++) {
            d2Vector g(
                direct::elem(dest.data,   x + 1, y) - direct::elem(dest.data,   x, y),
                direct::elem(dest.data,   x, y + 1) - direct::elem(dest.data,   x, y));
            d2Vector gs(
                direct::elem(smooth.data, x + 1, y) - direct::elem(smooth.data, x, y),
                direct::elem(smooth.data, x, y + 1) - direct::elem(smooth.data, x, y));

            const double iso = exp(-0.5 * gs.norm2() / tt);

            const double gsl = gs.length();
            if (gsl > 0.0) gs /= gsl;

            const d2Vector gn = g.dot(gs) * gs;
            const d2Vector gp = g - gn;

            flux(x, y, 0) = iso * g + (1.0 - iso) * gp;
        }


        #if JAZ_USE_OPENMP
        #pragma omp parallel for
        #endif
        for (long int y = 1; y < h - 1; y++)
        for (long int x = 1; x < w - 1; x++) {
            const double div = flux(x, y, 0).x - flux(x - 1, y, 0).x;
                             + flux(x, y, 0).y - flux(x, y - 1, 0).y;
            direct::elem(dest.data, x, y) += delta * div;
        }
    }
}

void FilterHelper::descendTV(const Image<RFLOAT> &src, Image<RFLOAT> &dest, double delta) {

    const long int w = src.data.xdim, h = src.data.ydim, d = src.data.zdim;
    std::vector<RFLOAT> vals;
    vals.reserve(6);

    for (long int z = 0; z < d; z++)
    for (long int y = 0; y < h; y++)
    for (long int x = 0; x < w; x++) {
        const double v0 = direct::elem(src.data, x, y, z);

        vals.clear();

        if (x > 0)     vals.push_back(direct::elem(src.data, x - 1, y, z));
        if (x < w - 1) vals.push_back(direct::elem(src.data, x + 1, y, z));
        if (y > 0)     vals.push_back(direct::elem(src.data, x, y - 1, z));
        if (y < h - 1) vals.push_back(direct::elem(src.data, x, y + 1, z));
        if (z > 0)     vals.push_back(direct::elem(src.data, x, y, z - 1));
        if (z < d - 1) vals.push_back(direct::elem(src.data, x, y, z + 1));

        std::vector<int> order = IndexSort<RFLOAT>::sortIndices(vals);
        const int c = vals.size();

        double vm;
        if (vals.size() % 2 == 0) {
            vm = 0.5 * (vals[order[c / 2]] + vals[order[c / 2 - 1]]);
        } else {
            vm = vals[order[c / 2]];
        }

        if (std::abs(v0 - vm) < delta) {
            direct::elem(dest.data, x, y, z) = vm;
        } else if (v0 < vm) {
            direct::elem(dest.data, x, y, z) = v0 + delta;
        } else {
            direct::elem(dest.data, x, y, z) = v0 - delta;
        }
    }
}

void FilterHelper::descendTV2(
    const Image<RFLOAT> &src, Image<RFLOAT> &dest,
    Volume<gravis::d3Vector>& xi, Volume<RFLOAT>& uBar,
    int iters, double sigma, double tau
) {
    const long int w = src.data.xdim, h = src.data.ydim, d = src.data.zdim;

    fwdGrad(src, xi);

    for (long int z = 0; z < d; z++)
    for (long int y = 0; y < h; y++)
    for (long int x = 0; x < w; x++) {
        uBar(x, y, z)                    = direct::elem(src.data, x, y, z);
        direct::elem(dest.data, x, y, z) = direct::elem(src.data, x, y, z);
    }

    for (int it = 0; it < iters; it++) {
        for (long int z = 0; z < d; z++)
        for (long int y = 0; y < h; y++)
        for (long int x = 0; x < w; x++) {
            d3Vector gradUBar;

            if (w == 1) {
                gradUBar.x = 0;
            } else if (x < w - 1) {
                gradUBar.x = uBar(x + 1, y, z) - uBar(x, y, z);
            } else {
                gradUBar.x = uBar(x, y, z) - uBar(x - 1, y, z);
            }

            if (h == 1) {
                gradUBar.y = 0;
            } else if (y < h - 1) {
                gradUBar.y = uBar(x, y + 1, z) - uBar(x, y, z);
            } else {
                gradUBar.y = uBar(x, y, z) - uBar(x, y - 1, z);
            }

            if (d == 1) {
                gradUBar.z = 0;
            } else if (z < d - 1) {
                gradUBar.z = uBar(x, y, z + 1) - uBar(x, y, z);
            } else {
                gradUBar.z = uBar(x, y, z) - uBar(x, y, z - 1);
            }

            const d3Vector nextXi = xi(x, y, z) + sigma * gradUBar;
            const double nxl = nextXi.length();
            xi(x, y, z) = nxl > 0.0 ? nextXi / nxl : d3Vector(0, 0, 0);
        }

        for (long int z = 0; z < d; z++)
        for (long int y = 0; y < h; y++)
        for (long int x = 0; x < w; x++) {

            double divXi = 0.0;
            if (x > 0)
                divXi += xi(x, y, z).x - xi(x - 1, y, z).x;
            if (y > 0)
                divXi += xi(x, y, z).y - xi(x, y - 1, z).y;
            if (z > 0)
                divXi += xi(x, y, z).z - xi(x, y, z - 1).z;

            const double du = tau * divXi;

            direct::elem(dest.data, x, y, z) += du;
            uBar(x, y, z) = direct::elem(dest.data, x, y, z) + 0.5 * du;
        }
    }
}

void FilterHelper::segmentTV(
    const Image<RFLOAT> &src, Image<RFLOAT> &dest,
    Volume<gravis::d3Vector>& xi, Volume<RFLOAT>& uBar,
    int iters, double sigma, double tau, double nu
) {
    const long int w = src.data.xdim, h = src.data.ydim, d = src.data.zdim;

    xi.fill(d3Vector(0.0, 0.0, 0.0));
    uBar.fill(0.0);

    for (long int z = 0; z < d; z++)
    for (long int y = 0; y < h; y++)
    for (long int x = 0; x < w; x++) {
        uBar(x, y, z) = direct::elem(src.data, x, y, z);
        direct::elem(dest.data, x, y, z) = 0.0;
    }

    for (int it = 0; it < iters; it++) {
        #if JAZ_USE_OPENMP
        #pragma omp parallel for
        #endif
        for (long int z = 0; z < d; z++)
        for (long int y = 0; y < h; y++)
        for (long int x = 0; x < w; x++) {
            d3Vector gradUBar;

            if (w == 1) {
                gradUBar.x = 0;
            } else if (x < w - 1) {
                gradUBar.x = uBar(x + 1, y, z) - uBar(x, y, z);
            } else {
                gradUBar.x = uBar(x, y, z) - uBar(x - 1, y, z);
            }

            if (h == 1) {
                gradUBar.y = 0;
            } else if (y < h - 1) {
                gradUBar.y = uBar(x, y + 1, z) - uBar(x, y, z);
            } else {
                gradUBar.y = uBar(x, y, z) - uBar(x, y - 1, z);
            }

            if (d == 1) {
                gradUBar.z = 0;
            } else if (z < d - 1) {
                gradUBar.z = uBar(x, y, z + 1) - uBar(x, y, z);
            } else {
                gradUBar.z = uBar(x, y, z) - uBar(x, y, z - 1);
            }

            const d3Vector nextXi = xi(x, y, z) + sigma * gradUBar;
            const double nxl = nextXi.length();
            xi(x, y, z) = nxl > 0.0 ? nextXi / nxl : d3Vector(0, 0, 0);
        }

        #if JAZ_USE_OPENMP
        #pragma omp parallel for
        #endif
        for (long int z = 0; z < d; z++)
        for (long int y = 0; y < h; y++)
        for (long int x = 0; x < w; x++) {

            double divXi = 0.0;
            if (x > 0)
                divXi += xi(x, y, z).x - xi(x - 1, y, z).x;
            if (y > 0)
                divXi += xi(x, y, z).y - xi(x, y - 1, z).y;
            if (z > 0)
                divXi += xi(x, y, z).z - xi(x, y, z - 1).z;

            const double u = direct::elem(dest.data, x, y, z);
            const double du = tau * (nu * divXi + direct::elem(src.data, x, y, z));
            const double nextU = clamp(u + du, 0.0, 1.0);

            direct::elem(dest.data, x, y, z) = nextU;
            uBar(x, y, z) = 2.0 * nextU - u;
        }
    }
}

void FilterHelper::segmentTVAniso2D(
    const Image<RFLOAT> &src, Image<RFLOAT> &dest,
    Volume<gravis::d2Vector>& xi, Volume<RFLOAT>& uBar,
    int iters, double sigma, double tau, double nu,
    double rho, double theta, double alpha
) {
    const long int w = src.data.xdim, h = src.data.ydim;

    const Image<RFLOAT> smooth = separableGaussianXY(src, rho);

    Volume<d2Vector> smoothGrad (w, h, 1);
    fwdGrad2D(smooth, smoothGrad);

    xi.fill(d2Vector(0.0, 0.0));
    uBar.fill(0.0);

    Volume<d2Matrix> D (w, h, 1);

    const double tt = theta * theta;

    for (long int y = 0; y < h; y++)
    for (long int x = 0; x < w; x++) {
        uBar(x, y, 0) = direct::elem(src.data, x, y);
        direct::elem(dest.data, x, y) = 0.0;

        d2Vector gs = smoothGrad(x, y, 0);
        const double iso = exp(-0.5 * gs.norm2() / tt);

        const double gsl = gs.length();
        if (gsl > 0.0) gs /= gsl;

        d2Matrix I;

        // d2Matrix G = E - d2Matrix(gs.x*gs.x, gs.y*gs.x, gs.x*gs.y, gs.y*gs.y);
        // G x = x - (x dot gs) gs

        d2Matrix F = d2Matrix(gs.x * gs.x, gs.y * gs.x, gs.x * gs.y, gs.y * gs.y);

        d2Matrix G = sqrt(alpha) * F + sqrt((3.0 - alpha) / 2.0) * (I - F);

        D(x, y, 0) = sqrt(nu) * (iso * I + (1.0 - iso) * G);
    }


    for (int it = 0; it < iters; it++) {
        #if JAZ_USE_OPENMP
        #pragma omp parallel for
        #endif
        for (long int y = 0; y < h; y++)
        for (long int x = 0; x < w; x++) {
            d2Vector gradUBar;

            if (w == 1) {
                gradUBar.x = 0;
            } else if (x < w - 1) {
                gradUBar.x = uBar(x + 1, y, 0) - uBar(x, y, 0);
            } else {
                gradUBar.x = uBar(x, y, 0) - uBar(x - 1, y, 0);
            }

            if (h == 1) {
                gradUBar.y = 0;
            } else if (y < h - 1) {
                gradUBar.y = uBar(x, y + 1, 0) - uBar(x, y, 0);
            } else {
                gradUBar.y = uBar(x, y, 0) - uBar(x, y - 1, 0);
            }

            d2Vector nextXi = D(x, y, 0) * (xi(x, y, 0) + sigma * D(x, y, 0) * gradUBar);

            double nxl = nextXi.length();

            xi(x, y, 0) = nxl > 0.0 ? nextXi / nxl : d2Vector(0, 0);
        }

        #if JAZ_USE_OPENMP
        #pragma omp parallel for
        #endif
        for (long int y = 0; y < h; y++)
        for (long int x = 0; x < w; x++) {

            double divXi = 0.0;
            if (x > 0)
                divXi += xi(x, y, 0).x - xi(x - 1, y, 0).x;
            if (y > 0)
                divXi += xi(x, y, 0).y - xi(x, y - 1, 0).y;

            const double u = direct::elem(dest.data, x, y);
            const double du = tau * (divXi + direct::elem(src.data, x, y));
            const double nextU = clamp(u + du, 0.0, 1.0);

            direct::elem(dest.data, x, y) = nextU;
            uBar(x, y, 0) = 2.0 * nextU - u;
        }
    }
}


void FilterHelper::fwdGrad(const Image<RFLOAT> &u, Volume<gravis::d3Vector> &dest) {
    const long int w = u.data.xdim, h = u.data.ydim, d = u.data.zdim;

    for (long int z = 0; z < d; z++)
    for (long int y = 0; y < h; y++)
    for (long int x = 0; x < w; x++) {
        if (w == 1) {
            dest(x, y, z).x = 0;
        } else if (x < w - 1) {
            dest(x, y, z).x = direct::elem(u.data, x + 1, y, z) - direct::elem(u.data, x, y, z);
        } else {
            dest(x, y, z).x = direct::elem(u.data, x, y, z) - direct::elem(u.data, x - 1, y, z);
        }

        if (h == 1) {
            dest(x, y, z).y = 0;
        } else if (y < h - 1) {
            dest(x, y, z).y = direct::elem(u.data, x, y + 1, z) - direct::elem(u.data, x, y, z);
        } else {
            dest(x, y, z).y = direct::elem(u.data, x, y, z) - direct::elem(u.data, x, y - 1, z);
        }

        if (d == 1) {
            dest(x, y, z).z = 0;
        } else if (z < d - 1) {
            dest(x, y, z).z = direct::elem(u.data, x, y, z + 1) - direct::elem(u.data, x, y, z);
        } else {
            dest(x, y, z).z = direct::elem(u.data, x, y, z) - direct::elem(u.data, x, y, z - 1);
        }
    }
}


void FilterHelper::fwdGrad2D(const Image<RFLOAT> &u, Volume<d2Vector> &dest) {
    const long int w = u.data.xdim, h = u.data.ydim, d = 1;

    for (long int z = 0; z < d; z++)
    for (long int y = 0; y < h; y++)
    for (long int x = 0; x < w; x++) {
        if (w == 1) {
            dest(x, y, z).x = 0;
        } else if (x < w - 1) {
            dest(x, y, z).x = direct::elem(u.data, x + 1, y, z) - direct::elem(u.data, x, y, z);
        } else {
            dest(x, y, z).x = direct::elem(u.data, x, y, z) - direct::elem(u.data, x - 1, y, z);
        }

        if (h == 1) {
            dest(x, y, z).y = 0;
        } else if (y < h - 1) {
            dest(x, y, z).y = direct::elem(u.data, x, y + 1, x, z) - direct::elem(u.data, x, y, z);
        } else {
            dest(x, y, z).y = direct::elem(u.data, x, y, z) - direct::elem(u.data, x, y - 1, z);
        }
    }
}

void FilterHelper::centralGrad2D(const Image<RFLOAT> &u, Volume<d2Vector> &dest) {
    const long int w = u.data.xdim, h = u.data.ydim, d = 1;

    for (long int z = 0; z < d; z++)
    for (long int y = 0; y < h; y++)
    for (long int x = 0; x < w; x++) {
        if (w == 1) {
            dest(x, y, z).x = 0;
        } else if (x < w - 1 && x > 0) {
            dest(x, y, z).x = (direct::elem(u.data, x + 1, y, z, 0) - direct::elem(u.data, x - 1, y, z, 0)) / 2.0;
        } else if (x == 0) {
            dest(x, y, z).x = direct::elem(u.data, x + 1, y, z, 0) - direct::elem(u.data, x, y, z);
        } else if (x == w - 1) {
            dest(x, y, z).x = direct::elem(u.data, x, y, z) - direct::elem(u.data, x - 1, y, z, 0);
        }

        if (h == 1) {
            dest(x, y, z).y = 0;
        } else if (y < h - 1 && y > 0) {
            dest(x, y, z).y = (direct::elem(u.data, x, y + 1, z, 0) - direct::elem(u.data, x, y - 1, z, 0)) / 2;
        } else if (y == 0) {
            dest(x, y, z).y = direct::elem(u.data, x, y + 1, z, 0) - direct::elem(u.data, x, y, z);
        } else if (y == h - 1) {
            dest(x, y, z).y = direct::elem(u.data, x, y, z) - direct::elem(u.data, x, y - 1, z, 0);
        }
    }

}

void FilterHelper::centralGrad2D(
    const Image<Complex> &u, Volume<d2Vector> &destRe, Volume<d2Vector> &destIm
) {
    const long int w = u.data.xdim, h = u.data.ydim, d = 1;

    for (long int z = 0; z < d; z++)
    for (long int y = 0; y < h; y++)
    for (long int x = 0; x < w; x++) {
        if (w == 1) {
            destRe(x, y, z).x = 0;
            destIm(x, y, z).x = 0;
        } else if (x < w - 1 && x > 0) {
            destRe(x, y, z).x = (direct::elem(u.data, x + 1, y, z, 0).real - direct::elem(u.data, x - 1, y, z, 0).real) / 2.0;
            destIm(x, y, z).x = (direct::elem(u.data, x + 1, y, z, 0).imag - direct::elem(u.data, x - 1, y, z, 0).imag) / 2.0;
        } else if (x == 0) {
            destRe(x, y, z).x = direct::elem(u.data, x + 1, y, z, 0).real - direct::elem(u.data, x, y, z).real;
            destIm(x, y, z).x = direct::elem(u.data, x + 1, y, z, 0).imag - direct::elem(u.data, x, y, z).imag;
        } else if (x == w - 1) {
            destRe(x, y, z).x = direct::elem(u.data, x, y, z).real - direct::elem(u.data, x - 1, y, z, 0).real;
            destIm(x, y, z).x = direct::elem(u.data, x, y, z).imag - direct::elem(u.data, x - 1, y, z, 0).imag;
        }

        if (h == 1) {
            destRe(x, y, z).y = 0;
            destIm(x, y, z).y = 0;
        } else if (y < h - 1 && y > 0) {
            destRe(x, y, z).y = (direct::elem(u.data, x, y + 1, z, 0).real - direct::elem(u.data, x, y - 1, z, 0).real) / 2;
            destIm(x, y, z).y = (direct::elem(u.data, x, y + 1, z, 0).imag - direct::elem(u.data, x, y - 1, z, 0).imag) / 2;
        } else if (y == 0) {
            destRe(x, y, z).y = direct::elem(u.data, x, y + 1, z, 0).real - direct::elem(u.data, x, y, z).real;
            destIm(x, y, z).y = direct::elem(u.data, x, y + 1, z, 0).imag - direct::elem(u.data, x, y, z).imag;
        } else if (y == h - 1) {
            destRe(x, y, z).y = direct::elem(u.data, x, y, z).real - direct::elem(u.data, x, y - 1, z, 0).real;
            destIm(x, y, z).y = direct::elem(u.data, x, y, z).imag - direct::elem(u.data, x, y - 1, z, 0).imag;
        }
    }

}

void FilterHelper::blendSoft(
    const Image<Complex> &src0, const Image<Complex> &src1,
    const Volume<RFLOAT>& mask, Image<Complex> &dest, RFLOAT bias1
) {
    const long int w = src0.data.xdim;
    const long int h = src0.data.ydim;
    const long int d = src0.data.zdim;

    for (long int z = 0; z < d; z++)
    for (long int y = 0; y < h; y++)
    for (long int x = 0; x < w; x++) {
        const Complex v0 = direct::elem(src0.data, x, y, z);
        const Complex v1 = direct::elem(src1.data, x, y, z);
        const RFLOAT m = mask(x, y, z);
        direct::elem(dest.data, x, y, z) = (v0 + bias1 * m * v1) / (1.0 + bias1 * m);
    }
}

double FilterHelper::totalVariation(const Image<RFLOAT> &src) {

    double sum = 0.0;
    for (long int j = 0; j < Ysize(src.data); j++)
    for (long int i = 0; i < Xsize(src.data); i++) {
        if (i == src.data.ydim - 1 || j == src.data.xdim - 1) continue;

        const double dx = direct::elem(src.data, i, j + 1) - direct::elem(src.data, i, j);
        const double dy = direct::elem(src.data, i + 1, j) - direct::elem(src.data, i, j);
        const double dtv = std::hypot(dx, dy);
        sum += dtv;
    }
    return sum;
}

double FilterHelper::totalLogVariation(const Image<RFLOAT> &src, double delta) {
    double sum = 0.0;

    for (long int j = 0; j < Ysize(src.data); j++)
    for (long int i = 0; i < Xsize(src.data); i++) {
        if (i == src.data.ydim - 1 || j == src.data.xdim - 1) continue;

        double v0 = direct::elem(src.data, i, j);
        double vx = direct::elem(src.data, i, j + 1);
        double vy = direct::elem(src.data, i + 1, j);

        double dx = vx - v0;
        double dy = vy - v0;

        double dtv = log(delta + std::hypot(dx, dy));

        sum += dtv;
    }

    return sum;
}

Image<RFLOAT> FilterHelper::averageX(Image<RFLOAT> copy, const Image<RFLOAT> &mask) {
    const long int& w = copy.data.xdim, h = copy.data.ydim, d = copy.data.zdim;
    for (size_t k = 0; k < d; k++)
    for (size_t j = 0; j < h; j++) {

        RFLOAT sum = 0.0, avg = 0.0;
        for (size_t i = 0; i < w; i++) {
            avg += direct::elem(mask.data, i, j, k) * direct::elem(copy.data, i, j, k);
            sum += direct::elem(mask.data, i, j, k);
        }
        if (sum > 0.0) avg /= sum;

        for (size_t i = 0; i < w; i++) {
            direct::elem(copy.data, i, j, k) = avg;
        }
    }
    return copy;
}

template<typename T>
void FilterHelper::binomial3x3_2D(const Image<T>& src, Image<T>& dest, bool wrap) {
    const std::array<double, 3> kernel = {0.25, 0.5, 0.25};
    MultidimArray<T> temp;
    if (wrap) {
        convolve_x_wrap(src.data, kernel, temp);
        convolve_y_wrap(temp,     kernel, dest.data);
    } else {
        convolve_x_nowrap(src.data, kernel, temp);
        convolve_y_nowrap(temp,     kernel, dest.data);
    }
}

template<typename T>
void FilterHelper::separableGaussian(const Volume<T>& src, Volume<T>& dest, double sigma, int k) {
    if (k < 0) k = 2 * sigma + 0.5;
    auto kernel = make_kernel(sigma, k);
    Volume<T> temp;
    convolve_x_nowrap(src,  kernel, dest);
    convolve_y_nowrap(dest, kernel, temp);
    convolve_z_nowrap(temp, kernel, dest);
}

// Manual template instantiation
using R = RFLOAT;
using C = Complex;
template Image<R> FilterHelper::cropCorner2D     (const Image<R> &img, int wF, int hF, R fill);
template Image<C> FilterHelper::cropCorner2D_fftw(const Image<C> &img, int wF, int hF, C fill);
