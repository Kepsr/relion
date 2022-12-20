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
    Matrix<RFLOAT> &u, Vector<RFLOAT> &w, Matrix<RFLOAT> &v,
    Vector<RFLOAT> &b, Vector<RFLOAT> &x
) {
    // Call the numerical recipes routine
    svbksb(
        nr::adaptForNumericalRecipes2(u),
        nr::adaptForNumericalRecipes(w),
        nr::adaptForNumericalRecipes2(v),
        u.nrows(), u.ncols(),
        nr::adaptForNumericalRecipes(b),
        nr::adaptForNumericalRecipes(x)
    );
    // Results are stored in x
}

template<typename T>
void svdcmp(
    const Matrix<T> &a,
    Matrix<RFLOAT> &u, Vector<RFLOAT> &w, Matrix<RFLOAT> &v
) {
    // svdcmp only works with RFLOAT
    u.resize(a.ncols(), a.nrows());
    std::copy(a.begin(), a.end(), u.begin());
    w.resize(u.ncols());
    std::fill(w.begin(), w.end(), 0);
    v.resize(u.ncols(), u.ncols());
    std::fill(u.begin(), u.end(), 0);
    // Call the numerical recipes routine
    svdcmp(u.data(), a.nrows(), a.ncols(), w.data(), v.data());
}

template<typename T>
Matrix<T> Matrix<T>::inv() const {

    if (size() == 0)
        REPORT_ERROR("Inverse: Matrix is empty");
    // Initialise output
    Matrix<T> inverse (ncols(), nrows());
    std::fill(inverse.begin(), inverse.end(), 0);

    if (ncols() == 3 && nrows() == 3) {
        int a, b, c, d;
        for (int i = 0; i <= 2; i++)
        for (int j = 0; j <= 2; j++) {
            a = (j - 1) % 3;
            b = (i - 1) % 3;
            c = (j + 1) % 3;
            d = (i + 1) % 3;
            inverse.at(i, j) = at(a, b) * at(c, d) - at(a, d) * at(c, b);
        }
        // Multiply first column of `this` with first row of `inverse`
        RFLOAT divisor = at(0, 0) * inverse.at(0, 0) 
                       + at(1, 0) * inverse.at(0, 1) 
                       + at(2, 0) * inverse.at(0, 2);
        return inverse / divisor;
    } else if (ncols() == 2 && nrows() == 2) {
        for (int i = 0; i <= 1; i++)
        for (int j = 0; j <= 1; j++) {
            int a = (j + 1) % 2;
            int b = (i + 1) % 2;
            inverse.at(i, j) = (i + j) % 2 ? -at(a, b) : +at(a, b);
        }
        RFLOAT divisor = at(0, 0) * at(1, 1) - at(0, 1) * at(1, 0);
        return inverse / divisor;
    } else {
        // Perform SVD
        Matrix<RFLOAT> u, v;
        Vector<RFLOAT> w;
        svdcmp(*this, u, w, v); // *this = U * W * V^t

        const T maximum = *std::max_element(begin(), end());
        RFLOAT tol = maximum * std::max(ncols(), nrows()) * 1e-14;

        // Compute W^-1
        bool invertible = false;
        for (RFLOAT &x: w) {
            if (abs(x) > tol) {
                x = 1.0 / x;
                invertible = true;
            } else {
                x = 0.0;
            }
        }

        if (!invertible) return inverse;

        // Compute V*W^-1
        for (int i = 0; i < v.nrows(); i++)
        for (int j = 0; j < v.ncols(); j++)
            v.at(i, j) *= w[j];

        // Compute inverse
        for (int i = 0; i < ncols(); i++)
        for (int j = 0; j < nrows(); j++)
        for (int k = 0; k < ncols(); k++)
            inverse.at(i, j) += (T) v.at(i, k) * u.at(j, k);

        return inverse;
    }
}

template<typename T>
void solve(Matrix<T> A, const Matrix<T> &b, Matrix<T> &x) {
    if (A.size() == 0)
        REPORT_ERROR("Solve: Matrix is empty");

    if (A.ncols() != A.nrows())
        REPORT_ERROR("Solve: Matrix is not square");

    if (A.nrows() != b.nrows())
        REPORT_ERROR("Solve: Different sizes of A and b");

    // Solve
    x = b;
    gaussj(
        nr::adaptForNumericalRecipes2(A), A.nrows(),
        nr::adaptForNumericalRecipes2(x), x.ncols()
    );
}

// Solve a system of linear equations (Ax = b) by SVD
template<typename T>
void solve(const Matrix<T> &A, Vector<T> b, Vector<RFLOAT> &x, RFLOAT epsilon) {
    if (A.ncols() == 0)
        REPORT_ERROR("Solve: Matrix is empty");

    // if (A.ncols() != A.nrows())
    //     REPORT_ERROR("Solve: Matrix is not square");

    if (A.nrows() != b.size())
        REPORT_ERROR("Solve: Differently sized Matrix and Vector");

    // if (b.isRow())
    //     REPORT_ERROR("Solve: Incorrect vector shape");

    // First perform SVD
    // Xmipp interface that calls to svdcmp of numerical recipes
    Matrix<RFLOAT> u, v;
    Vector<RFLOAT> w;
    svdcmp(A, u, w, v);

    // Check if eigenvalues of the SVD are acceptable.
    // To improve the routine's precision,
    // values smaller than epsilon are made zero.
    setSmallValuesToZero(w.begin(), w.end(), epsilon);

    // Set size of matrices
    x.resize(b.size());

    // Xmipp interface that calls to svdksb of numerical recipes
    svbksb(u, w, v, b, x);
}

// https://isocpp.org/wiki/faq/templates#separate-template-fn-defn-from-decl
template class Matrix<RFLOAT>;
template void solve(const Matrix<RFLOAT>&, Vector<RFLOAT>, Vector<RFLOAT>&, RFLOAT);
template void solve(Matrix<RFLOAT>, const Matrix<RFLOAT>&, Matrix<RFLOAT>&);
