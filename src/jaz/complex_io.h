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

#ifndef COMPLEX_IO_H
#define COMPLEX_IO_H

#include <src/image.h>
#include <src/complex.h>
#include <string>

namespace ComplexIO {

    template <typename T>
    void write(const MultidimArray<tComplex<T> > &img, std::string fnBase, std::string fnSuffix) {
        Image<RFLOAT> realimg(img.xdim, img.ydim, img.zdim, img.ndim);
        Image<RFLOAT> imagimg(img.xdim, img.ydim, img.zdim, img.ndim);

        FOR_ALL_DIRECT_NZYX_ELEMENTS_IN_MULTIDIMARRAY(img) {
            direct::elem(realimg.data, i, j, k, l) = direct::elem(img, i, j, k, l).real;
            direct::elem(imagimg.data, i, j, k, l) = direct::elem(img, i, j, k, l).imag;
        }

        realimg.write(fnBase + "_real" + fnSuffix);
        imagimg.write(fnBase + "_imag" + fnSuffix);
    }

    template <typename T>
    void read(Image<tComplex<T> > &img, std::string fnBase, std::string fnSuffix) {
        Image<RFLOAT> realimg = Image<RFLOAT>::from_filename(fnBase + "_real" + fnSuffix);
        Image<RFLOAT> imagimg = Image<RFLOAT>::from_filename(fnBase + "_imag" + fnSuffix);

        img = Image<Complex>(realimg.data.xdim, realimg.data.ydim, realimg.data.zdim, realimg.data.ndim);
        FOR_ALL_DIRECT_NZYX_ELEMENTS_IN_MULTIDIMARRAY(img.data) {
            direct::elem(img.data, i, j, k, l) = Complex(
                direct::elem(realimg.data, i, j, k, l),
                direct::elem(imagimg.data, i, j, k, l)
            );
        }
    }

};

#endif
