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

#ifndef TENSOR_2X2_H
#define TENSOR_2X2_H

#include <cmath>
#include <vector>
#include <src/error.h>
#include <src/jaz/index_sort.h>
#include <src/jaz/gravis/t2Vector.h>
#include <src/jaz/gravis/t2Matrix.h>

extern "C" {
    #include <src/jaz/d3x3/dsyev2.h>
}

/* Symmetric 2x2 matrix to be used as e.g. a structure tensor of a 2D image */
template <class T>
class Tensor2x2 {

    public:

    T xx, xy, yy;

    Tensor2x2() {}

    Tensor2x2(T t): xx(t), xy(t), yy(t) {}

    Tensor2x2(T xx, T xy, T yy): xx(xx), xy(xy), yy(yy) {}

    static Tensor2x2<T> dyadicProduct(const gravis::t2Vector<T>& v0, const gravis::t2Vector<T>& v1);

    void diagonalize(gravis::t2Vector<T>& eigenvalues, gravis::t2Matrix<T>& eigenvectors) const;

    gravis::t2Matrix<T> toMatrix() const;

    Tensor2x2<T>& operator += (const Tensor2x2<T>& arg) {
        xx += arg.xx;
        xy += arg.xy;
        yy += arg.yy;
        return *this;
    }

    Tensor2x2& operator -= (const Tensor2x2& arg) {
        xx -= arg.xx;
        xy -= arg.xy;
        yy -= arg.yy;
        return *this;
    }

    Tensor2x2& operator *= (T arg) {
        xx *= arg;
        xy *= arg;
        yy *= arg;
        return *this;
    }

    Tensor2x2& operator /= (T arg) {
        xx /= arg;
        xy /= arg;
        yy /= arg;
        return *this;
    }

};

template <class T>
Tensor2x2<T> Tensor2x2<T>::dyadicProduct(const gravis::t2Vector<T>& v0, const gravis::t2Vector<T>& v1) {
    return {v0.x * v1.x, v0.x * v1.y, v0.y * v1.y};
}

template <class T>
void Tensor2x2<T>::diagonalize(gravis::t2Vector<T>& eigenvalues, gravis::t2Matrix<T>& eigenvectors) const {
    dsyev2(xx, xy, yy, &eigenvalues[0], &eigenvalues[1], &eigenvectors(0, 0), &eigenvectors(0, 1));
    eigenvectors(1, 0) = -eigenvectors(0, 1);
    eigenvectors(1, 1) = +eigenvectors(0, 0);
}

template <class T>
gravis::t2Matrix<T> Tensor2x2<T>::toMatrix() const {
    return {xx, xy, xy, yy};
}

template <class T>
inline Tensor2x2<T> operator - (const Tensor2x2<T>& v1) {
    return {-v1.xx, -v1.xy, -v1.yy};
}

template <class T>
inline Tensor2x2<T> operator + (Tensor2x2<T> v1, const Tensor2x2<T>& v2) {
    return v1 += v2;
}

template <class T>
inline Tensor2x2<T> operator * (Tensor2x2<T> v1, T arg) {
    return v1 *= arg;
}

template <class T>
inline Tensor2x2<T> operator / (Tensor2x2<T> v1, T arg) {
    return v1 /= arg;
}

template <class T>
inline Tensor2x2<T> operator * (const Tensor2x2<T>& v, float f) {
    return v *= f;
}

template <class T>
inline Tensor2x2<T> operator * (const Tensor2x2<T>& v, double f) {
    return v *= f;
}

template <class T>
inline Tensor2x2<T> operator - (Tensor2x2<T> v1, const Tensor2x2<T>& v2) {
    return v1 -= v2;
}

template <class T>
inline Tensor2x2<T> operator * (float f, const Tensor2x2<T>& v) {
    return {f * v.xx, f * v.xy, f * v.yy};
}

template <class T>
inline Tensor2x2<T> operator * (double f, const Tensor2x2<T>& v) {
    return {f * v.xx, f * v.xy, f * v.yy};
}

template <class T>
inline std::ostream& operator<< (std::ostream& os, const Tensor2x2<T>& arg) {
    os << std::setprecision(17)  << "["
        << std::setw(8) << arg.xx << ", "
        << std::setw(8) << arg.xy << ", "
        << std::setw(8) << arg.yy << "]";
    return os;
}

#endif
