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
            for (int i = Yinit(v); i <= Ylast(v); i++) {
                for (int j = Xinit(v); j <= Xlast(v); j++)
                    ostrm << "(" << v.elem(k, i, j).real << "," << v.elem(k, i, j).imag << ")" << ' ';
                ostrm << std::endl;
            }
        }
    }
    return ostrm;
}


/** Array (vector) by scalar.
 *
 * Take a vector and a constant,
 * and apply the appropriate operation element-wise.
 * This is the function which really implements the operations.
 * Simple calls to it perform much faster than calls to the corresponding operators.
 * It is supposed to be hidden from users.
 *
 * This function is not ported to Python.
 */
template <typename T, typename Op>
inline MultidimArray<T> arrayByScalar(
    const MultidimArray<T> &input, T scalar, MultidimArray<T> &output,
    Op operation
) {
    if (!output.data || !output.sameShape(input)) { output.resize(input); }
    T *iptr = input.data, *optr = output.data;
    // These two pointers will move through (respectively) input and output.
    // *iptr will be used to assign *optr.
    for (long int n = 0; n < input.xdim * input.ydim * input.zdim; ++n, ++optr, ++iptr) {
        *optr = operation(*iptr, scalar);
    }
    return output;
}

template <typename T>
MultidimArray<T> MultidimArray<T>::operator + (const T scalar) const {
    MultidimArray<T> output;
    return arrayByScalar(*this, scalar, output, [] (const T x, T y) { return x + y; });
}

template <typename T>
MultidimArray<T> MultidimArray<T>::operator - (const T scalar) const {
    MultidimArray<T> output;
    return arrayByScalar(*this, scalar, output, [] (const T x, T y) { return x - y; });
}

template <typename T>
MultidimArray<T> MultidimArray<T>::operator * (const T scalar) const {
    MultidimArray<T> output;
    return arrayByScalar(*this, scalar, output, [] (const T x, T y) { return x * y; });
}

template <typename T>
MultidimArray<T> MultidimArray<T>::operator / (const T scalar) const {
    MultidimArray<T> output;
    return arrayByScalar(*this, scalar, output, [] (const T x, T y) { return x / y; });
}

template <typename T>
MultidimArray<T> MultidimArray<T>::operator += (const T scalar) {
    return arrayByScalar(*this, scalar, *this, [] (const T x, T y) { return x + y; });
}

template <typename T>
MultidimArray<T> MultidimArray<T>::operator -= (const T scalar) {
    return arrayByScalar(*this, scalar, *this, [] (const T x, T y) { return x - y; });
}

template <typename T>
MultidimArray<T> MultidimArray<T>::operator *= (const T scalar) {
    return arrayByScalar(*this, scalar, *this, [] (const T x, T y) { return x * y; });
}

template <typename T>
MultidimArray<T> MultidimArray<T>::operator /= (const T scalar) {
    return arrayByScalar(*this, scalar, *this, [] (const T x, T y) { return x / y; });
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
template <typename T, typename Op>
inline MultidimArray<T> arrayByArray(
    const MultidimArray<T> &arg1, const MultidimArray<T> &arg2,
    MultidimArray<T> &output,
    Op operation
) {
    if (!arg1.sameShape(arg2)) {
        arg1.printShape();
        arg2.printShape();
        REPORT_ERROR((std::string) "Array_by_array: different shapes");
    }
    if (!output.data || !output.sameShape(arg1)) { output.resize(arg1); }
    T *arg1ptr, *arg2ptr, *optr;
    long int n;
    for (
        n = 0, optr = output.data, arg1ptr = arg1.data, arg2ptr = arg2.data;
        n < arg1.xdim * arg1.ydim * arg1.zdim;
        ++n, ++arg1ptr, ++arg2ptr, ++optr
    ) {
        *optr = operation(*arg1ptr, *arg2ptr);
    }
}

template <typename T>
MultidimArray<T> MultidimArray<T>::operator + (const MultidimArray<T> &arg) const {
    MultidimArray<T> output;
    return arrayByArray(*this, arg, output, [] (T x, T y) { return x + y; });
}

template <typename T>
MultidimArray<T> MultidimArray<T>::operator - (const MultidimArray<T> &arg) const {
    MultidimArray<T> output;
    return arrayByArray(*this, arg, output, [] (T x, T y) { return x - y; });
}

template <typename T>
MultidimArray<T> MultidimArray<T>::operator * (const MultidimArray<T> &arg) const {
    MultidimArray<T> output;
    return arrayByArray(*this, arg, output, [] (T x, T y) { return x * y; });
}

template <typename T>
MultidimArray<T> MultidimArray<T>::operator / (const MultidimArray<T> &arg) const {
    MultidimArray<T> output;
    return arrayByArray(*this, arg, output, [] (T x, T y) { return x / y; });
}

template <typename T>
MultidimArray<T> MultidimArray<T>::operator += (const MultidimArray<T> &arg) {
    return arrayByArray(*this, arg, *this, [] (T x, T y) { return x + y; });
}

template <typename T>
MultidimArray<T> MultidimArray<T>::operator -= (const MultidimArray<T> &arg) {
    return arrayByArray(*this, arg, *this, [] (T x, T y) { return x - y; });
}

template <typename T>
MultidimArray<T> MultidimArray<T>::operator *= (const MultidimArray<T> &arg) {
    return arrayByArray(*this, arg, *this, [] (T x, T y) { return x * y; });
}

template <typename T>
MultidimArray<T> MultidimArray<T>::operator /= (const MultidimArray<T> &arg) {
    return arrayByArray(*this, arg, *this, [] (T x, T y) { return x / y; });
}

/** Scalar by array.
 *
 * This function must take one scalar and a vector, and operate element by
 * element according to the operation required. This is the function which
 * really implements the operations. Simple calls to it perform much faster
 * than calls to the corresponding operators. Although it is supposed to
 * be a hidden function not useable by normal programmers.
 *
 * This function is not ported to Python.
 */
template <typename T, typename Op>
inline MultidimArray<T> scalarByArray(
    const T scalar,
    const MultidimArray<T> &input,
    MultidimArray<T> &output,
    Op operation
) {
    if (!output.data || !output.sameShape(input)) { output.resize(input); }
    T *iptr = output.data, *optr = input.data;
    for (long int n = 0; n < input.xdim * input.ydim * input.zdim; ++n, ++optr, ++iptr)
        *optr = operation(scalar, *iptr);
    return output;
}

template <typename T>
MultidimArray<T> operator + (const T scalar, const MultidimArray<T> &input) {
    MultidimArray<T> output;
    return scalarByArray(scalar, input, output, [] (const T x, T y) { return x + y; });
}

template <typename T>
MultidimArray<T> operator - (const T scalar, const MultidimArray<T> &input) {
    MultidimArray<T> output;
    return scalarByArray(scalar, input, output, [] (const T x, T y) { return x - y; });
}

template <typename T>
MultidimArray<T> operator * (const T scalar, const MultidimArray<T> &input) {
    MultidimArray<T> output;
    return scalarByArray(scalar, input, output, [] (const T x, T y) { return x * y; });
}

template <typename T>
MultidimArray<T> operator / (const T scalar, const MultidimArray<T> &input) {
    MultidimArray<T> output;
    return scalarByArray(scalar, input, output, [] (const T x, T y) { return x / y; });
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
