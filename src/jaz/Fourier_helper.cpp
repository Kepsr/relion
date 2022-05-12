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

    if (abs(xshift) < Xmipp::epsilon && abs(yshift) < Xmipp::epsilon) {
        return;
    }

    for (long int yy = 0; yy < h; yy++)
    for (long int xx = 0; xx < w; xx++) {
        RFLOAT x = xx;
        RFLOAT y = yy < w ? yy : yy - h;

        Complex c1;
        RFLOAT dotp = -2.0 * PI * (x * xshift + y * yshift);
        #ifdef RELION_SINGLE_PRECISION
        SINCOSF(dotp, &c1.imag, &c1.real);
        #else
        SINCOS(dotp, &c1.imag, &c1.real);
        #endif

        Complex c2 = direct::elem(img, yy, xx);

        direct::elem(img, yy, xx) = Complex(
            c1.real * c2.real - c1.imag * c2.imag,  // c1 dot conj c2
            c1.real * c2.imag + c1.imag * c2.real   // c1 dot (i conj c2)
        );
    }
}

void FourierHelper::FourierShift2D(
    MultidimArray<RFLOAT> &img, RFLOAT xshift, RFLOAT yshift
) {
    FourierTransformer ft;
    MultidimArray<Complex> imgC;

    ft.FourierTransform(img, imgC);
    //FourierShift2D(imgC, xshift, yshift);
    shiftImageInFourierTransform(imgC, imgC, img.ydim, xshift, yshift);

    ft.inverseFourierTransform(imgC, img);
}
