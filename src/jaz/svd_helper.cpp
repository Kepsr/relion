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

#include <src/jaz/svd_helper.h>
#include <src/jaz/index_sort.h>

void SvdHelper::decompose(
    const Matrix<RFLOAT> &A,
    Matrix<RFLOAT> &U,
    Vector<RFLOAT> &S,
    Matrix<RFLOAT> &Vt
) {
    Matrix<RFLOAT> U0, Vt0;
    Vector<RFLOAT> S0;
    svdcmp(A, U0, S0, Vt0);

    const int rc = A.nrows(), cc = A.ncols();

    std::vector<RFLOAT> Svec (cc);
    std::copy_n(S0.begin(), cc, Svec.begin());

    std::vector<int> order = IndexSort<RFLOAT>::sortIndices(Svec);

    U = Matrix<RFLOAT>(rc, cc);
    S = Vector<RFLOAT>(cc);
    Vt = Matrix<RFLOAT>(cc, cc);

    for (int i = 0; i < cc; i++) {
        const int j = order[cc - i - 1];

        for (int c = 0; c < cc; c++) {
            Vt(c, i) = Vt0(c, j);
        }

        S[i] = S0[j];

        for (int r = 0; r < rc; r++) {
            U(r, i) = U0(r, j);
        }
    }
}
