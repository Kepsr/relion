/***************************************************************************
 *
 * Author: "Sjors H.W. Scheres"
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
/***************************************************************************
 *
 * Authors:     Carlos Oscar S. Sorzano (coss@cnb.csic.es)
 *
 * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
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
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
 * 02111-1307  USA
 *
 *  All comments concerning this program package may be sent to the
 *  e-mail address 'xmipp@cnb.csic.es'
 ***************************************************************************/

#include "src/multidim_array.h"
#include <unordered_map>

// Show a complex array ---------------------------------------------------
std::ostream& operator << (std::ostream &ostrm, const MultidimArray<Complex> &v) {

    if (v.xdim == 0) {
        ostrm << "NULL MultidimArray\n";
    } else {
        ostrm << std::endl;
    }

    for (int l = 0; l < Nsize(v); l++) {
        if (Nsize(v) > 1) ostrm << "Image No. " << l << std::endl;
        for (int k = Zinit(v); k <= Zlast(v); k++) {
            if (Zsize(v) > 1) { ostrm << "Slice No. " << k << std::endl; }
            for (int j = Yinit(v); j <= Ylast(v); j++) {
                for (int i = Xinit(v); i <= Xlast(v); i++)
                    ostrm << "(" << v.elem(i, j, k).real << "," << v.elem(i, j, k).imag << ")" << ' ';
                ostrm << std::endl;
            }
        }
    }
    return ostrm;
}

template <typename T>
void threshold_abs_above(T *ptr, T a, T b) {
    if (abs(*ptr) > a) { *ptr = b * sgn(*ptr); }
}

template <typename T>
void threshold_abs_below(T *ptr, T a, T b) {
    if (abs(*ptr) < a) { *ptr = b * sgn(*ptr); }
}

template <typename T>
void threshold_above(T *ptr, T a, T b) {
    if (*ptr > a) { *ptr = b; }
}

template <typename T>
void threshold_below(T *ptr, T a, T b) {
    if (*ptr < a) { *ptr = b; }
}

template <typename T>
void threshold_range(T *ptr, T a, T b) {
    if (*ptr < a) { *ptr = a; } else
    if (*ptr > b) { *ptr = b; }
}

template <typename T>
void MultidimArray<T>::threshold(const std::string &type, T a, T b, MultidimArray<int> *mask) {

    static const std::unordered_map<std::string, void (*)(T *ptr, T a, T b)> s2f {
        {"abs_above", &threshold_abs_above},
        {"abs_below", &threshold_abs_below},
        {"above",     &threshold_above},
        {"below",     &threshold_below},
        {"range",     &threshold_range},
    };

    const auto it = s2f.find(type);
    if (it == s2f.end())
        REPORT_ERROR(static_cast<std::string>("Threshold: mode not supported (" + type + ")"));
    const auto f = it->second;

    int *maskptr = mask ? mask->begin() : nullptr;
    for (T *ptr = begin(); ptr != end(); ++ptr, ++maskptr) {
        if (!mask || *maskptr > 0) f(ptr, a, b);
    }
}

template <typename T>
MultidimArray<T> MultidimArray<T>::operator + (const T scalar) const {
    auto copy (*this);
    return copy += scalar;
}

template <typename T>
MultidimArray<T> MultidimArray<T>::operator - (const T scalar) const {
    auto copy (*this);
    return copy -= scalar;
}

template <typename T>
MultidimArray<T> MultidimArray<T>::operator * (const T scalar) const {
    auto copy (*this);
    return copy *= scalar;
}

template <typename T>
MultidimArray<T> MultidimArray<T>::operator / (const T scalar) const {
    auto copy (*this);
    return copy /= scalar;
}

template <typename T>
MultidimArray<T>& MultidimArray<T>::operator += (const T scalar) {
    for (auto &x : *this) { x += scalar; }
    return *this;
}

template <typename T>
MultidimArray<T>& MultidimArray<T>::operator -= (const T scalar) {
    for (auto &x : *this) { x -= scalar; }
    return *this;
}

template <typename T>
MultidimArray<T>& MultidimArray<T>::operator *= (const T scalar) {
    for (auto &x : *this) { x *= scalar; }
    return *this;
}

template <typename T>
MultidimArray<T>& MultidimArray<T>::operator /= (const T scalar) {
    for (auto &x : *this) { x /= scalar; }
    return *this;
}

/** Array by array
 *
 * This function must take two vectors of the same size, and operate element
 * by element according to the operation required. This is the function
 * which really implements the operations. Simple calls to it perform much
 * faster than calls to the corresponding operators. Although it is supposed
 * to be a hidden function not useable by normal programmers.
 *
 */
template <typename T>
inline MultidimArray<T>& pointwise(
    MultidimArray<T> &arg1, const MultidimArray<T> &arg2,
    T (*operation)(T, T)
) {
    if (!arg1.sameShape(arg2)) {
        arg1.printShape();
        arg2.printShape();
        REPORT_ERROR((std::string) "Array_by_array: different shapes");
    }
    for (
        T *ptr1 = arg1.data, *ptr2 = arg2.data,
          *end = arg1.data + arg1.xdim * arg1.ydim * arg1.zdim;
        ptr1 != end;
        ++ptr1, ++ptr2
    ) {
        *ptr1 = operation(*ptr1, *ptr2);
    }
    return arg1;
}

template <typename T>
MultidimArray<T> MultidimArray<T>::operator + (const MultidimArray<T> &arg) const {
    auto copy (*this);
    return copy += arg;
}

template <typename T>
MultidimArray<T> MultidimArray<T>::operator - (const MultidimArray<T> &arg) const {
    auto copy (*this);
    return copy -= arg;
}

template <typename T>
MultidimArray<T> MultidimArray<T>::operator * (const MultidimArray<T> &arg) const {
    auto copy (*this);
    return copy *= arg;
}

template <typename T>
MultidimArray<T> MultidimArray<T>::operator / (const MultidimArray<T> &arg) const {
    auto copy (*this);
    return copy /= arg;
}

template <typename T>
MultidimArray<T>& MultidimArray<T>::operator += (const MultidimArray<T> &arg) {
    return pointwise(*this, arg, +[] (T x, T y) -> T { return x + y; });
}

template <typename T>
MultidimArray<T>& MultidimArray<T>::operator -= (const MultidimArray<T> &arg) {
    return pointwise(*this, arg, +[] (T x, T y) -> T { return x - y; });
}

template <typename T>
MultidimArray<T>& MultidimArray<T>::operator *= (const MultidimArray<T> &arg) {
    return pointwise(*this, arg, +[] (T x, T y) -> T { return x * y; });
}

template <typename T>
MultidimArray<T>& MultidimArray<T>::operator /= (const MultidimArray<T> &arg) {
    return pointwise(*this, arg, +[] (T x, T y) -> T { return x / y; });
}

template <typename T>
MultidimArray<T> operator + (const T scalar, const MultidimArray<T> &input) {
    auto copy (input);
    for (auto &x : copy) { x = scalar + x; }
    return copy;
}

template <typename T>
MultidimArray<T> operator - (const T scalar, const MultidimArray<T> &input) {
    auto copy (input);
    for (auto &x : copy) { x = scalar - x; }
    return copy;
}

template <typename T>
MultidimArray<T> operator * (const T scalar, const MultidimArray<T> &input) {
    auto copy (input);
    for (auto &x : copy) { x = scalar * x; }
    return copy;
}

template <typename T>
MultidimArray<T> operator / (const T scalar, const MultidimArray<T> &input) {
    auto copy (input);
    for (auto &x : copy) { x = scalar / x; }
    return copy;
}

template class MultidimArray<float>;
template class MultidimArray<double>;
template class MultidimArray<unsigned short>;
template class MultidimArray<short>;
template class MultidimArray<unsigned char>;
template class MultidimArray<signed char>;
template class MultidimArray<int>;
template class MultidimArray<Complex>;

// Required for compilation of ml_model.cpp
template MultidimArray<double> operator * (double scalar, const MultidimArray<double> &input);
