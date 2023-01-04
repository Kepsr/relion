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

#ifndef FILTER_HELPER_H
#define FILTER_HELPER_H

#include <src/image.h>
#include <src/fftw.h>
#include <src/strings.h>
#include <src/ctf.h>
#include <src/jaz/gravis/t3Vector.h>
#include <src/jaz/gravis/t2Vector.h>
#include <src/jaz/gravis/t2Matrix.h>
#include <src/jaz/gravis/t3Matrix.h>
#include <src/jaz/volume.h>
#include <src/jaz/obs_model.h>

namespace FilterHelper {

    template<typename T>
    void binomial3x3_2D(const Image<T>& src, Image<T>& dest, bool wrap = false);

    template<typename T>
    void separableGaussian(const Volume<T>& src, Volume<T>& dest, double sigma, int k = -1);

    template<typename T>
    void separableGaussian(const MultidimArray<T>& src, MultidimArray<T>& dest, double sigma, int k = -1);

    template<typename T>
    void separableGaussianWrap(const MultidimArray<T>& src, MultidimArray<T>& dest, double sigma, int k = -1);

    void separableGaussianFreq(
            const MultidimArray<Complex>& src, MultidimArray<Complex>& dest,
            double sigma, int k = -1);

    void separableGaussianFreqXY(
            const MultidimArray<Complex>& src, MultidimArray<Complex>& dest,
            double sigma, int k = -1);

    void drawTestPattern(Image<RFLOAT>& img, int squareSize);
    void drawTestPattern(Volume<RFLOAT>& volume, int squareSize);

    Image<RFLOAT> expImg(const Image<RFLOAT>& img, double scale = 1.0);
    Image<RFLOAT> logImg(const Image<RFLOAT>& img, double scale = 1.0, double thresh = 1e-20);

    template <typename R>
    Image<R> padCorner2D     (const Image<R> &img, int w, int h, R fill = 0);
    template <typename C>
    Image<C> padCorner2D_fftw(const Image<C> &img, int w, int h, C fill = 0);

    template <typename R>
    Image<R> cropCorner2D     (const Image<R> &img, int w, int h, R fill = 0);
    template <typename C>
    Image<C> cropCorner2D_fftw(const Image<C> &img, int w, int h, C fill = 0);

    Image<RFLOAT> zeroOutsideCorner2D(const Image<RFLOAT>& img, double radius);
    void GaussianEnvelopeCorner2D(Image<RFLOAT>& img, double sigma);
    template <typename R>
    Image<R> raisedCosEnvCorner2D(const Image<R>& img, double radIn, double radOut);
    Image<RFLOAT> raisedCosEnvCorner3D(Image<RFLOAT>& img, double radIn, double radOut);
    Image<RFLOAT> raisedCosEnvFreq2D(const Image<RFLOAT>& img, double radIn, double radOut);
    Image<RFLOAT> raisedCosEnvRingFreq2D(const Image<RFLOAT>& img, double rad0, double rad1, double stepWidth);


    void lowPassFilter(Image<RFLOAT>& img, double maxFreq0, double maxFreq1, Image<RFLOAT>& dest);
    void lowPassFilterSpectrum(MultidimArray<Complex>& spectrum, double maxFreq0, double maxFreq1);

    RFLOAT averageValue(Image<RFLOAT>& img);

    void phaseFlip(Image<RFLOAT>& img, CTF& ctf, ObservationModel *obsModel, int opticsGroup, RFLOAT angpix, Image<RFLOAT>& dest);
    void applyBeamTilt(Image<RFLOAT>& img, RFLOAT beamtilt_x, RFLOAT beamtilt_y,
                                RFLOAT lambda, RFLOAT Cs, RFLOAT angpix, int s, Image<RFLOAT>& dest);
    void modulate      (Image<RFLOAT>& img, const CTF& ctf, ObservationModel *obsModel, int opticsGroup, RFLOAT angpix, Image<RFLOAT>& dest);
    void modulate (Image<Complex>& imgFreq, const CTF& ctf, ObservationModel *obsModel, int opticsGroup, RFLOAT angpix, Image<RFLOAT>& dest);
    void modulate(MultidimArray<Complex>& imgFreq, const CTF& ctf, ObservationModel *obsModel, int opticsGroup, RFLOAT angpix);
    void drawCtf                           (CTF& ctf, ObservationModel *obsModel, int opticsGroup, RFLOAT angpix, Image<Complex>& dest);
    void wienerFilter  (Image<RFLOAT>& img, CTF& ctf, ObservationModel *obsModel, int opticsGroup, RFLOAT angpix, RFLOAT eps, RFLOAT Bfac, Image<RFLOAT>& dest);
    void richardsonLucy(Image<RFLOAT>& img, CTF& ctf, ObservationModel *obsModel, int opticsGroup, RFLOAT angpix, RFLOAT eps, int iterations, Image<RFLOAT>& dest);
    void rampFilter(Image<RFLOAT>& img, RFLOAT s0, RFLOAT t1, double ux, double uy, Image<RFLOAT>& dest);
    void rampFilter3D(Image<Complex>& img, RFLOAT s0, RFLOAT t1, double tx, double ty, double tz);
    void doubleRampFilter3D(Image<Complex>& img, RFLOAT s0, RFLOAT t1, double tx, double ty, double tz);

    Image<RFLOAT> getPhase(const Image<Complex>& img);
    Image<RFLOAT> getAbs  (const Image<Complex>& img);
    Image<RFLOAT> getReal (const Image<Complex>& img);
    Image<RFLOAT> getImag (const Image<Complex>& img);

    void powerSpectrum2D(Image<RFLOAT>& img, Volume<RFLOAT>& spectrum);
    void equiphaseAverage2D(const Volume<RFLOAT>& src, Volume<RFLOAT>& dest);

    void threshold(Image<RFLOAT>& src, RFLOAT t, Image<RFLOAT>& dest);
    void linearTransform(Image<RFLOAT>& src, RFLOAT m, RFLOAT q, Image<RFLOAT>& dest);
    void sumUp(const std::vector<Image<RFLOAT> >& src, Image<RFLOAT>& dest);

    double L1distance(const Image<RFLOAT>& i0, const Image<RFLOAT>& i1, int x0 = 0, int y0 = 0, int w = -1, int h = -1);
    double L2distance(const Image<RFLOAT>& i0, const Image<RFLOAT>& i1, int x0 = 0, int y0 = 0, int w = -1, int h = -1);
    double NCC(const Image<RFLOAT>& i0, const Image<RFLOAT>& i1, int x0 = 0, int y0 = 0, int w = -1, int h = -1);

    void wienerDivide(Image<RFLOAT>& num, Image<RFLOAT>& denom, RFLOAT eps, Image<RFLOAT>& dest);
    void divideExcessive(Image<Complex>& num, Volume<RFLOAT>& denom, RFLOAT theta, Image<Complex>& dest);
    void wienerDeconvolve(Image<Complex>& num, Image<Complex>& denom, RFLOAT theta, Image<Complex>& dest);

    void extract2D(const Image<RFLOAT>& src, Image<RFLOAT>& dest,
                        long int x0, long int y0,
                        long int w, long int h);

    void extract(const Volume<RFLOAT>& src,
                        Volume<RFLOAT>& dest,
                        long int x0, long int y0, long int z0,
                        long int w, long int h, long int d);

    void signedDist(const Image<RFLOAT>& src, Image<RFLOAT>& dest);
    void erode3x3(Image<RFLOAT>& src, Image<RFLOAT>& dest);
    void localMinima(Image<RFLOAT>& src, RFLOAT thresh, Image<RFLOAT>& dest);
    std::vector<gravis::d3Vector> localMinima(Image<RFLOAT>& src, RFLOAT thresh);

    void uniqueInfluenceMask(std::vector<gravis::d3Vector> pts, Image<RFLOAT>& dest, Image<RFLOAT>& indexDest, RFLOAT thresh);
    void polarRemap(gravis::d2Vector pos, const Image<RFLOAT>& src, Image<RFLOAT>& dest,
                            const Image<RFLOAT>& mask, Image<RFLOAT>& maskDest, int phiRes, int rRes, double rMax);
    void polarRemap(gravis::d2Vector pos, const Image<RFLOAT>& distTransf, const Image<RFLOAT>& src, Image<RFLOAT>& dest,
                            const Image<RFLOAT>& mask, Image<RFLOAT>& maskDest, int phiRes, int rRes, double rMax);

    Image<RFLOAT> cartToPolar(const Image<RFLOAT>& img);
    Image<RFLOAT> polarToCart(const Image<RFLOAT>& img);
    Image<RFLOAT> polarBlur(const Image<RFLOAT>& img, double sigma);
    Image<RFLOAT> sectorBlend(const Image<RFLOAT>& img0, const Image<RFLOAT>& img1, int sectors);


    void diffuseAlongIsocontours2D(const Image<RFLOAT>& src, const Image<RFLOAT>& guide,
                                            Image<RFLOAT>& dest, int iters, RFLOAT sigma, RFLOAT lambda, RFLOAT delta);

    void EED_2D(const Image<RFLOAT>& src, Image<RFLOAT>& dest, int iters, double sigma, double delta, double tau);

    void descendTV(const Image<RFLOAT>& src, Image<RFLOAT>& dest, double delta);
    void descendTV2(const Image<RFLOAT>& src, Image<RFLOAT>& dest,
                            Volume<gravis::d3Vector>& xi, Volume<RFLOAT>& uBar,
                            int iters, double sigma, double tau);

    void segmentTV(const Image<RFLOAT>& src, Image<RFLOAT>& dest,
                            Volume<gravis::d3Vector>& xi, Volume<RFLOAT>& uBar,
                            int iters, double sigma, double tau, double nu);

    void segmentTVAniso2D(const Image<RFLOAT>& src, Image<RFLOAT>& dest,
                            Volume<gravis::d2Vector>& xi, Volume<RFLOAT>& uBar,
                            int iters, double sigma, double tau, double nu,
                            double rho, double theta, double alpha);

    void fwdGrad(const Image<RFLOAT>& u, Volume<gravis::d3Vector>& dest);
    void fwdGrad2D(const Image<RFLOAT>& u, Volume<gravis::d2Vector>& dest);
    void centralGrad2D(const Image<RFLOAT>& u, Volume<gravis::d2Vector>& dest);
    void centralGrad2D(const Image<Complex>& u, Volume<gravis::d2Vector>& destRe, Volume<gravis::d2Vector>& destIm);

    void blendSoft(const Image<Complex>& src0, const Image<Complex>& src1,
                            const Volume<RFLOAT>& mask, Image<Complex>& dest, RFLOAT bias1 = 1.0);

    double totalVariation(const Image<RFLOAT>& src);
    double totalLogVariation(const Image<RFLOAT>& src, double delta = 1.0);

    Image<RFLOAT> separableGaussianX  (const Image<RFLOAT>& src, RFLOAT sigma, int k = -1, bool wrap = false);
    Image<RFLOAT> separableGaussianXY (const Image<RFLOAT>& src, RFLOAT sigma, int k = -1, bool wrap = false);
    Image<RFLOAT> separableGaussianXYZ(const Image<RFLOAT>& src, RFLOAT sigma, int k = -1, bool wrap = false);

    Image<RFLOAT> averageX(Image<RFLOAT> copy, const Image<RFLOAT> &mask);

    void centralGradient(const Volume<RFLOAT>& src, Volume<gravis::t3Vector<RFLOAT> >& dest);
    gravis::t3Vector<RFLOAT> centralGradient(const Volume<RFLOAT>& src, size_t x, size_t y, size_t z);

    MultidimArray<Complex> FriedelExpand(const MultidimArray<Complex>& half);

    Image<RFLOAT> normaliseToUnitInterval(const Image<RFLOAT>& img);
    Image<RFLOAT> normaliseToUnitIntervalSigned(const Image<RFLOAT>& img);

}

template<typename T>
void FilterHelper::separableGaussian(const Volume<T>& src, Volume<T>& dest, double sigma, int k)
{
    if (k < 0)
    {
        k = (int)(2*sigma + 0.5);
    }

    dest.resize(src);

    std::vector<double> kernel(2*k+1);
    const double s2 = sigma*sigma;

    for (int i = -k; i <= k; i++)
    {
        kernel[i+k] = exp(-i*i/s2);
    }

    Volume<T> temp(src.dimx, src.dimy, src.dimz);

    for (size_t z = 0; z < src.dimz; z++)
    for (size_t y = 0; y < src.dimy; y++)
    for (size_t x = 0; x < src.dimx; x++)
    {
        T v = 0;
        double m = 0;

        for (long int i = -k; i <= k; i++)
        {
            const long int xx = x + i;
            if (xx < 0 || xx >= src.dimx) continue;

            v += kernel[i+k] * src(xx,y,z);
            m += kernel[i+k];
        }

        dest(x,y,z) = v/m;
    }

    for (size_t z = 0; z < src.dimz; z++)
    for (size_t y = 0; y < src.dimy; y++)
    for (size_t x = 0; x < src.dimx; x++)
    {
        T v = 0;
        double m = 0;

        for (long int i = -k; i <= k; i++)
        {
            const long int yy = y + i;
            if (yy < 0 || yy >= src.dimy) continue;

            v += kernel[i+k] * dest(x,yy,z);
            m += kernel[i+k];
        }

        temp(x,y,z) = v/m;
    }

    for (size_t z = 0; z < src.dimz; z++)
    for (size_t y = 0; y < src.dimy; y++)
    for (size_t x = 0; x < src.dimx; x++)
    {
        T v = 0;
        double m = 0;

        for (long int i = -k; i <= k; i++)
        {
            const long int zz = z + i;
            if (zz < 0 || zz >= src.dimz) continue;

            v += kernel[i+k] * temp(x,y,zz);
            m += kernel[i+k];
        }

        dest(x,y,z) = v/m;
    }
}

template<typename T>
void FilterHelper::separableGaussian(const MultidimArray<T>& src, MultidimArray<T>& dest, double sigma, int k)
{
    if (k < 0)
    {
        k = (int)(2*sigma + 0.5);
    }

    dest.reshape(src);

    std::vector<double> kernel(2*k+1);
    const double s2 = sigma*sigma;

    for (int i = -k; i <= k; i++)
    {
        kernel[i+k] = exp(-0.5*i*i/s2);
    }

    MultidimArray<T> temp(src.zdim, src.ydim, src.xdim);

    for (size_t z = 0; z < src.zdim; z++)
    for (size_t y = 0; y < src.ydim; y++)
    for (size_t x = 0; x < src.xdim; x++)
    {
        T v = 0;
        double m = 0;

        for (long int i = -k; i <= k; i++)
        {
            const long int xx = x + i;
            if (xx < 0 || xx >= src.xdim) continue;

            v += kernel[i+k] * direct::elem(src, y, xx, z, 0);
            m += kernel[i+k];
        }

        direct::elem(dest, y, x, z, 0) = v/m;
    }

    for (size_t z = 0; z < src.zdim; z++)
    for (size_t y = 0; y < src.ydim; y++)
    for (size_t x = 0; x < src.xdim; x++)
    {
        T v = 0;
        double m = 0;

        for (long int i = -k; i <= k; i++)
        {
            const long int yy = y + i;
            if (yy < 0 || yy >= src.ydim) continue;

            v += kernel[i+k] * direct::elem(dest, yy, x, z, 0);
            m += kernel[i+k];
        }

        direct::elem(temp, y, x, z, 0) = v / m;
    }

    for (size_t z = 0; z < src.zdim; z++)
    for (size_t y = 0; y < src.ydim; y++)
    for (size_t x = 0; x < src.xdim; x++)
    {
        T v = 0;
        double m = 0;

        for (long int i = -k; i <= k; i++)
        {
            const long int zz = z + i;
            if (zz < 0 || zz >= src.zdim) continue;

            v += kernel[i+k] * direct::elem(temp, y, x, zz, 0);
            m += kernel[i+k];
        }

        direct::elem(dest, y, x, z, 0) = v/m;
    }
}

template<typename T>
void FilterHelper::separableGaussianWrap(const MultidimArray<T>& src, MultidimArray<T>& dest, double sigma, int k)
{
    if (k < 0)
    {
        k = (int)(2*sigma + 0.5);
    }

    dest.reshape(src);

    std::vector<double> kernel(2*k+1);
    const double s2 = sigma*sigma;

    for (int i = -k; i <= k; i++)
    {
        kernel[i+k] = exp(-0.5*i*i/s2);
    }

    MultidimArray<T> temp(src.zdim, src.ydim, src.xdim);

    for (size_t z = 0; z < src.zdim; z++)
    for (size_t y = 0; y < src.ydim; y++)
    for (size_t x = 0; x < src.xdim; x++)
    {
        T v = 0;
        double m = 0;

        for (long int i = -k; i <= k; i++)
        {
            const long int xx = (src.xdim + x + i) % src.xdim;

            v += kernel[i+k] * direct::elem(src, y, xx, z, 0);
            m += kernel[i+k];
        }

        direct::elem(dest, y, x, z, 0) = v/m;
    }

    for (size_t z = 0; z < src.zdim; z++)
    for (size_t y = 0; y < src.ydim; y++)
    for (size_t x = 0; x < src.xdim; x++)
    {
        T v = 0;
        double m = 0;

        for (long int i = -k; i <= k; i++)
        {
            const long int yy = (src.ydim + y + i) % src.ydim;

            v += kernel[i+k] * direct::elem(dest, yy, x, z, 0);
            m += kernel[i+k];
        }

        direct::elem(temp, y, x, z, 0) = v/m;
    }

    for (size_t z = 0; z < src.zdim; z++)
    for (size_t y = 0; y < src.ydim; y++)
    for (size_t x = 0; x < src.xdim; x++)
    {
        T v = 0;
        double m = 0;

        for (long int i = -k; i <= k; i++)
        {
            const long int zz = (src.zdim + z + i) % src.zdim;

            v += kernel[i+k] * direct::elem(temp, y, x, zz, 0);
            m += kernel[i+k];
        }

        direct::elem(dest, y, x, z, 0) = v/m;
    }
}

#endif
