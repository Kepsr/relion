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
                    ostrm << "(" << A3D_ELEM(v, k, i, j).real << "," << A3D_ELEM(v, k, i, j).imag <<")" << ' ';
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
    return coreArrayByScalar(input, scalar, output, operation);
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

template class MultidimArray<RFLOAT>;
template class MultidimArray<short>;
