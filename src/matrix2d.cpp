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

/* Is diagonal ------------------------------------------------------------- */
#include "src/matrix2d.h"

/* Interface to numerical recipes: svbksb ---------------------------------- */
void svbksb(
    Matrix2D<RFLOAT> &u, Matrix1D<RFLOAT> &w, Matrix2D<RFLOAT> &v,
    Matrix1D<RFLOAT> &b, Matrix1D<RFLOAT> &x
) {
    // Call to the numerical recipes routine. Results will be stored in X
    svbksb(
        u.adaptForNumericalRecipes2(),
        w.adaptForNumericalRecipes(),
        v.adaptForNumericalRecipes2(),
        u.mdimy, u.mdimx,
        b.adaptForNumericalRecipes(),
        x.adaptForNumericalRecipes()
    );
}

template <typename T>
void Matrix2D<T>::setSmallValuesToZero(RFLOAT accuracy) {
    for (auto &x : *this)
    if (abs(x) < accuracy) { x = 0.0; }
}

template <typename T>
T Matrix2D<T>::max() const {
    if (mdim <= 0) return static_cast<T>(0);

    T maxval = mdata[0];
    for (auto &x : *this) if (x > maxval) { maxval = x; }
    return maxval;
}

template <typename T>
T Matrix2D<T>::min() const {
    if (mdim <= 0) return static_cast<T>(0);

    T minval = mdata[0];
    for (auto &x : *this) if (x < minval) { minval = x; }
    return minval;
}

template<typename T>
void Matrix2D<T>::inv(Matrix2D<T> &result) const {

    if (mdimx == 0 || mdimy == 0)
        REPORT_ERROR("Inverse: Matrix is empty");
    // Initialise output
    result.initZeros(mdimx, mdimy);

    if (mdimx == 3 && mdimy == 3) {
        int a, b, c, d;
        for (int i = 0; i <= 2; i++)
        for (int j = 0; j <= 2; j++) {
                a = (j - 1) % 3;
                b = (i - 1) % 3;
                c = (j + 1) % 3;
                d = (i + 1) % 3;
                result.at(i, j) = at(a, b) * at(c, d) - at(a, d) * at(c, b);
        }
        // Multiply first column of `this` with first row of `result`
        RFLOAT divisor = at(0, 0) * result.at(0, 0) 
                       + at(1, 0) * result.at(0, 1) 
                       + at(2, 0) * result.at(0, 2);
        result /= divisor;
    } else if (mdimx == 2 && mdimy == 2) {
        int sign, a, b;
        for (int i = 0; i <= 1; i++)
        for (int j = 0; j <= 1; j++) {
            sign = (i + j) % 2 == 0 ? 1 : -1;
            a = (j + 1) % 2;  // logical negation
            b = (i + 1) % 2;
            result.at(i, j) = sign * at(a, b);
        }
        RFLOAT divisor = at(0, 0) * at(1, 1) - at(0, 1) * at(1, 0);
        result /= divisor;
    } else {
        // Perform SVD
        Matrix2D<RFLOAT> u, v;
        Matrix1D<RFLOAT> w;
        svdcmp(*this, u, w, v); // *this = U * W * V^t

        RFLOAT tol = max() * std::max(mdimx, mdimy) * 1e-14;

        // Compute W^-1
        bool invertible = false;
        for (RFLOAT &x : w) {
            if (abs(x) > tol) {
                x = 1.0 / x;
                invertible = true;
            } else {
                x = 0.0;
            }
        }

        if (!invertible) return;

        // Compute V*W^-1
        for (int i = 0; i < v.mdimy; i++)
        for (int j = 0; j < v.mdimx; j++)
        v.at(i, j) *= w[j];

        // Compute inverse
        for (int i = 0; i < mdimx; i++)
        for (int j = 0; j < mdimy; j++)
        for (int k = 0; k < mdimx; k++)
        result.at(i, j) += (T) v.at(i, k) * u.at(j, k);
    }
}

// Solve a system of linear equations (Ax = b) by SVD
template<typename T>
void solve(
    const Matrix2D<T> &A, const Matrix1D<T> &b,
    Matrix1D<RFLOAT> &result, RFLOAT tolerance
) {
    if (A.mdimx == 0)
        REPORT_ERROR("Solve: Matrix is empty");

    /*if (A.mdimx != A.mdimy)
        REPORT_ERROR("Solve: Matrix is not square");*/

    if (A.mdimy != b.vdim)
        REPORT_ERROR("Solve: Differently sized Matrix and Vector");

    /*if (b.isRow())
        REPORT_ERROR("Solve: Incorrect vector shape");*/

    // First perform SVD
    // Xmipp interface that calls to svdcmp of numerical recipes
    Matrix2D<RFLOAT> u, v;
    Matrix1D<RFLOAT> w;
    svdcmp(A, u, w, v);

    // Check if eigenvalues of the SVD are acceptable.
    // If a value is lower than the tolerance, it is made zero,
    // to improve the routine's precision.
    for (RFLOAT &x : w) if (x < tolerance) { x = 0; }

    // Set size of matrices
    result.resize(b.vdim);

    // Xmipp interface that calls to svdksb of numerical recipes
    Matrix1D<RFLOAT> bd;
    typeCast(b, bd);
    svbksb(u, w, v, bd, result);
}

// https://isocpp.org/wiki/faq/templates#separate-template-fn-defn-from-decl
template void solve(
    const Matrix2D<RFLOAT>&, const Matrix1D<RFLOAT>&,
    Matrix1D<RFLOAT>&, RFLOAT 
);

template class Matrix2D<RFLOAT>;
