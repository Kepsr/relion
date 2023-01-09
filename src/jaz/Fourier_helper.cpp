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

#include <src/jaz/Fourier_helper.h>

void FourierHelper::FourierShift2D(MultidimArray<Complex> &img, RFLOAT xshift, RFLOAT yshift) {
    const long w = img.xdim;
    const long h = img.ydim;

    xshift /= h;
    yshift /= h;

    if (abs(xshift) < Xmipp::epsilon<RFLOAT>() && abs(yshift) < Xmipp::epsilon<RFLOAT>()) return;

    for (long int j = 0; j < h; j++)
    for (long int i = 0; i < w; i++) {
        const RFLOAT x = i;
        const RFLOAT y = j < w ? j : j - h;
        const Complex z = Complex::unit(-2.0 * PI * (x * xshift + y * yshift));
        direct::elem(img, i, j) *= z;
    }
}

void FourierHelper::FourierShift2D(
    MultidimArray<RFLOAT> &img, RFLOAT xshift, RFLOAT yshift
) {
    FourierTransformer ft;
    MultidimArray<Complex> &imgC = ft.FourierTransform(img);
    // FourierShift2D(imgC, xshift, yshift);
    shiftImageInFourierTransform(imgC, img.ydim, xshift, yshift);
    img = ft.inverseFourierTransform(imgC);
}
