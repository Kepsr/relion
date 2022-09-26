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

#ifndef SLICE_HELPER_H
#define SLICE_HELPER_H

#include <src/image.h>
#include <src/fftw.h>
#include <src/strings.h>
#include <src/jaz/volume.h>
#include <src/jaz/gravis/t4Matrix.h>
#include <src/jaz/gravis/t2Vector.h>


namespace SliceHelper {

    void affineTransform(const Image<RFLOAT> &img, gravis::d4Matrix A, Image<RFLOAT> &dest);

    Image<RFLOAT> downsample(Image<RFLOAT> &img, int factor);
    void downsampleSlices(const Image<RFLOAT> &img, Image<RFLOAT> &dest);
    void downsampleSlicesReal(const Image<RFLOAT> &img, Image<RFLOAT> &dest);

    void lowPassFilterSlicewise(Image<RFLOAT> &img, double maxFreq0, double maxFreq1);
    void lowPassFilterSlice(Image<RFLOAT> &img, long int n, double maxFreq0, double maxFreq1);

    void subsample(const Image<RFLOAT> &img, Image<RFLOAT> &dest);

    template <int D>
    void avgPad(const Image<RFLOAT> &src, Image<RFLOAT> &dest, double ratio);

    void halveSpectrum2D(const Image<Complex> &src, Image<Complex> &dest);

    void extractSpectralSlice(
        Image<Complex> &src, Image<RFLOAT> &dest,
        gravis::d3Matrix proj, gravis::d2Vector volCentImg, double oversample = 4.0);

    void insertSpectralSlices(
        std::vector<Image<RFLOAT> > &src,
        std::vector<gravis::d3Matrix> proj,
        std::vector<gravis::d2Vector> volCentImg,
        Image<Complex> &dest, double thickness = 1.0, double thicknessSlope = 0.0, double imgPad = 0.5);

    void insertWeightedSpectralSlices(
        std::vector<Image<RFLOAT> > &src,
        std::vector<gravis::d3Matrix> proj,
        std::vector<gravis::d2Vector> volCentImg,
        std::vector<double> imgWeights,
        Image<Complex> &dest, double thickness = 1.0, double imgPad = 0.5);

    Image<RFLOAT> getStackSlice(const Image<RFLOAT> &src, long int s);

    Image<RFLOAT> getStackSlices(const Image<RFLOAT> &src, long int s, long int ndim);

    template <typename T>
    void insertStackSlice(const Image<T> &src, Image<T> &dest, long int n);

    template <typename T>
    void insertZSlice(const Image<T> &src, Image<T> &dest, long int z);

    template <typename T>
    Image<T> consolidate(const std::vector<Image<T> > &src, bool toN = false);

    void stat(const Image<RFLOAT> &img);

};

#endif
