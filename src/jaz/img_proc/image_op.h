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

#ifndef IMAGE_OPS_H
#define IMAGE_OPS_H

#include <src/image.h>

namespace ImageOp {

    template<typename T1, typename T2>
    void linearCombination(const Image<T1> &src0, const Image<T1> &src1, T2 a0, T2 a1, Image<T1> &dest);

    template<typename T1, typename T2, typename T3>
    void multiply(const Image<T1> &i0, const Image<T2>& i1, Image<T3>& dest);

    template<typename T1, typename T2>
    void multiplyBy(Image<T1> &dest, const Image<T2>& i1);

    template<typename T1, typename T2>
    void linearCombination(const Image<T1> &src0, T1 src1, T2 a0, T2 a1, Image<T1> &dest);

    template<typename T1, typename T2>
    void linearCombination(const MultidimArray<T1> &src0, const MultidimArray<T1> &src1, T2 a0, T2 a1, MultidimArray<T1> &dest);

    template<typename T1, typename T2, typename T3>
    void multiply(const MultidimArray<T1> &i0, const MultidimArray<T2>& i1, MultidimArray<T3>& dest);

    template<typename T1, typename T2>
    void linearCombination(const MultidimArray<T1> &src0, T1 src1, T2 a0, T2 a1, MultidimArray<T1> &dest);

    template<typename T>
    void flipX(const MultidimArray<T> &src, MultidimArray<T> &dest);

    template<typename T>
    void flipY(const MultidimArray<T> &src, MultidimArray<T> &dest);

    template<typename T>
    void flipZ(const MultidimArray<T> &src, MultidimArray<T> &dest);

    template <typename T>
    void invert_hand(const MultidimArray<T> &src, MultidimArray<T> &dest);

    template <typename T>
    void flipYAxis(MultidimArray<T> &array);

    template<typename T1>
    void rotate90(const MultidimArray<T1> &src0, MultidimArray<T1> &dest);

    template<typename T1>
    void rotate180(const MultidimArray<T1> &src0, MultidimArray<T1> &dest);

    template<typename T1>
    void rotate270(const MultidimArray<T1> &src0, MultidimArray<T1> &dest);

};

template<typename T1, typename T2>
void ImageOp::linearCombination(const Image<T1> &src0, const Image<T1> &src1, T2 a0, T2 a1, Image<T1> &dest) {
    for (long int i = 0; i < src0.data.size(); i++) {
        dest.data[i] = a0 * src0.data[i] + a1 * src1.data[i];
    }
}

template<typename T1, typename T2, typename T3>
void ImageOp::multiply(const Image<T1> &i0, const Image<T2>& i1, Image<T3>& dest) {
    dest = Image<T3>(i0.data.xdim, i0.data.ydim, i0.data.zdim, i0.data.ndim);
    for (long int i = 0; i < i0.data.size(); i++) {
        dest.data[i] = i0.data[i] * i1.data[i];
    }
}

template<typename T1, typename T2>
void ImageOp::multiplyBy(Image<T1> &dest, const Image<T2>& i1) {
    for (long int i = 0; i < dest.data.size(); i++) {
        dest.data[i] *= i1.data[i];
    }
}

template<typename T1, typename T2>
void ImageOp::linearCombination(const Image<T1> &src0, T1 src1, T2 a0, T2 a1, Image<T1> &dest) {
    for (long int i = 0; i < src0.data.size(); i++) {
        dest.data[i] = a0 * src0.data[i] + a1 * src1;
    }
}

template<typename T1, typename T2>
void ImageOp::linearCombination(const MultidimArray<T1> &src0, const MultidimArray<T1> &src1, T2 a0, T2 a1, MultidimArray<T1> &dest) {
    for (long int i = 0; i < src0.size(); i++) {
        dest[i] = a0 * src0[i] + a1 * src1[i];
    }
}

template<typename T1, typename T2, typename T3>
void ImageOp::multiply(const MultidimArray<T1> &i0, const MultidimArray<T2>& i1, MultidimArray<T3>& dest) {
    dest = MultidimArray<T3>(i0.xdim, i0.ydim, i0.zdim, i0.ndim);
    for (long int i = 0; i < i0.size(); i++) {
        dest[i] = i0[i] * i1[i];
    }
}

template<typename T1, typename T2>
void ImageOp::linearCombination(const MultidimArray<T1> &src0, T1 src1, T2 a0, T2 a1, MultidimArray<T1> &dest) {
    for (long int i = 0; i < src0.size(); i++) {
        dest[i] = a0 * src0[i] + a1 * src1;
    }
}

// flip 'left-to-right'
// MotionCor2's FlipGain 2
template<typename T1>
void ImageOp::flipX(const MultidimArray<T1> &src, MultidimArray<T1> &dest) {
    dest.reshape(src);
    for (long int n = 0; n < src.ndim; n++)
    for (long int z = 0; z < src.zdim; z++)
    for (long int y = 0; y < src.ydim; y++)
    for (long int x = 0; x < src.xdim; x++) {
        direct::elem(dest, x, y, z, n) = direct::elem(src, src.xdim - 1 - x, y, z, n);
    }
}

// flip 'upside down'
// MotionCor2's FlipGain 1
template<typename T1>
void ImageOp::flipY(const MultidimArray<T1> &src, MultidimArray<T1> &dest) {
    dest.reshape(src);
    for (long int n = 0; n < src.ndim; n++)
    for (long int z = 0; z < src.zdim; z++)
    for (long int y = 0; y < src.ydim; y++)
    for (long int x = 0; x < src.xdim; x++) {
        direct::elem(dest, x, y, z, n) = direct::elem(src, x, src.ydim - 1 - y, z, n);
    }
}


template <typename T>
void ImageOp::flipZ(const MultidimArray<T> &src, MultidimArray<T> &dest) {
    dest.reshape(src);
    for (long int n = 0; n < src.ndim; n++)
    for (long int z = 0; z < src.zdim; z++)
    for (long int y = 0; y < src.ydim; y++)
    for (long int x = 0; x < src.xdim; x++) {
        direct::elem(dest, x, y, z, n) = direct::elem(src, x, y, src.zdim - 1 - z, n);
    }
}

template <typename T>
void ImageOp::invert_hand(const MultidimArray<T> &src, MultidimArray<T> &dest) {
    dest.reshape(src);
    for (long int k = 0; k < Zsize(src); k++)
    for (long int j = 0; j < Ysize(src); j++)
    for (long int i1 = 0; i1 < Xsize(src); i1++) {
        long int i2 = Xsize(src) - i1;
        direct::elem(dest, i1, j, k) = direct::elem(src, i2, j, k);
    }
}

template <typename T>
void ImageOp::flipYAxis(MultidimArray<T> &array) {
    const int ylim = array.ydim / 2, z = 0;
    for (int n = 0; n < array.zdim; n++)
    for (int y1 = 0; y1 < ylim; y1++) {
        const int y2 = array.ydim - 1 - y1;
        for (int x = 0; x < array.xdim; x++) {
            /// TODO: memcpy or pointer arithmetic is probably faster
            std::swap(
                direct::elem(array, x, y1, z, n),
                direct::elem(array, x, y2, z, n)
            );
        }
    }
}

// This is equivalent to MotionCor2's -RotGain 1.
// In relion_display, this looks clock-wise.
template<typename T1>
void ImageOp::rotate90(const MultidimArray<T1> &src0, MultidimArray<T1> &dest) {
    dest.reshape(src0.xdim, src0.ydim, src0.zdim, src0.ndim);
    for (long int n = 0; n < src0.ndim; n++)
    for (long int z = 0; z < src0.zdim; z++)
    for (long int y = 0; y < src0.ydim; y++)
    for (long int x = 0; x < src0.xdim; x++) {
        direct::elem(dest, x, src0.ydim - 1 - y, z, n) = direct::elem(src0, x, y, z, n);
    }
}

// MotionCor2's RotGain 2
template<typename T1>
void ImageOp::rotate180(const MultidimArray<T1> &src0, MultidimArray<T1> &dest) {
    dest.reshape(src0);
    for (long int n = 0; n < src0.ndim; n++)
    for (long int z = 0; z < src0.zdim; z++)
    for (long int y = 0; y < src0.ydim; y++)
    for (long int x = 0; x < src0.xdim; x++) {
        direct::elem(dest, src0.xdim - 1 - x, src0.ydim - 1 - y, z, n) = direct::elem(src0, x, y, z, n);
    }
}

// MotionCor2's RotGain 3
template<typename T1>
void ImageOp::rotate270(const MultidimArray<T1> &src0, MultidimArray<T1> &dest) {
    dest.reshape(src0.ydim, src0.xdim, src0.zdim, src0.ndim);
    for (long int n = 0; n < src0.ndim; n++)
    for (long int z = 0; z < src0.zdim; z++)
    for (long int y = 0; y < src0.ydim; y++)
    for (long int x = 0; x < src0.xdim; x++) {
        direct::elem(dest, src0.xdim - 1 - x, y, z, n) = direct::elem(src0, x, y, z, n);
    }
}

#endif
