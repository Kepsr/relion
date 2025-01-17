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
/*
 * backprojector.cpp
 *
 *  Created on: 24 Aug 2010
 *      Author: scheres
 */

#include "src/backprojector.h"

void BackProjector::initialiseDataAndWeight(int current_size) {
    initialiseData(current_size);
    weight.resize(data);
}

void BackProjector::initZeros(int current_size) {
    initialiseDataAndWeight(current_size);
    data.initZeros();
    weight.initZeros();
}

inline void fillm(Matrix<RFLOAT> *magMatrix, RFLOAT &m00, RFLOAT &m10, RFLOAT &m01, RFLOAT &m11) {

    if (magMatrix != 0) {
        m00 = (*magMatrix)(0, 0);
        m10 = (*magMatrix)(1, 0);
        m01 = (*magMatrix)(0, 1);
        m11 = (*magMatrix)(1, 1);
    } else {
        m00 = 1.0;
        m10 = 0.0;
        m01 = 0.0;
        m11 = 1.0;
    }

}

void BackProjector::backproject2Dto3D(
    const MultidimArray<Complex> &f2d,
    const Matrix<RFLOAT> &A,
    const MultidimArray<RFLOAT> *Mweight,
    RFLOAT r_ewald_sphere, RFLOAT curvature,
    Matrix<RFLOAT> *magMatrix
) {

    RFLOAT m00, m10, m01, m11;
    fillm(magMatrix, m00, m10, m01, m11);

    // Use the inverse matrix
    Matrix<RFLOAT> Ainv = A.inv();

    // Go from the 2D slice coordinates to the 3D coordinates
    Ainv *= (RFLOAT) padding_factor;  // take scaling into account directly

    // max_r2 and min_r2_nn are defined in 3D-space
    const int max_r2    = round(r_max    * padding_factor) * round(r_max    * padding_factor);
    const int min_r2_nn = round(r_min_nn * padding_factor) * round(r_min_nn * padding_factor);

    // precalculated coefficients for ellipse determination (see further down)

    // Make sure A contains 2D distortion (lowercase 2D, uppercase 3D):
    // Matrix-multiply Ainv(3×2) and m(2×2) to get Am
    const RFLOAT Am_Xx = Ainv(0, 0) * m00 + Ainv(0, 1) * m10;
    const RFLOAT Am_Xy = Ainv(0, 0) * m01 + Ainv(0, 1) * m11;
    const RFLOAT Am_Yx = Ainv(1, 0) * m00 + Ainv(1, 1) * m10;
    const RFLOAT Am_Yy = Ainv(1, 0) * m01 + Ainv(1, 1) * m11;
    const RFLOAT Am_Zx = Ainv(2, 0) * m00 + Ainv(2, 1) * m10;
    const RFLOAT Am_Zy = Ainv(2, 0) * m01 + Ainv(2, 1) * m11;

    // Precompute AmTAm into ATA:
    const RFLOAT AtA_xx = Am_Xx * Am_Xx + Am_Yx * Am_Yx + Am_Zx * Am_Zx;
    const RFLOAT AtA_xy = Am_Xx * Am_Xy + Am_Yx * Am_Yy + Am_Zx * Am_Zy;
    const RFLOAT AtA_yy = Am_Xy * Am_Xy + Am_Yy * Am_Yy + Am_Zy * Am_Zy;
    const RFLOAT AtA_xy2 = AtA_xy * AtA_xy;

    // #define DEBUG_BACKP
    #ifdef DEBUG_BACKP
    std::cerr << " Xsize(f2d)= " << Xsize(f2d) << std::endl;
    std::cerr << " Ysize(f2d)= " << Ysize(f2d) << std::endl;
    std::cerr << " Xsize(data)= " << Xsize(data) << std::endl;
    std::cerr << " Ysize(data)= " << Ysize(data) << std::endl;
    std::cerr << " Xinit(data)= " << Xinit(data) << std::endl;
    std::cerr << " Yinit(data)= " << Yinit(data) << std::endl;
    std::cerr << " Zinit(data)= " << Zinit(data) << std::endl;
    std::cerr << " r_max= "<< r_max << std::endl;
    std::cerr << " Ainv= " << Ainv << std::endl;
    #endif

    // Precalculate inverse of Ewald sphere diameter outside loop
    RFLOAT inv_diam_ewald = r_ewald_sphere > 0.0 ? curvature / (2.0 * r_ewald_sphere) : 0.0;

    const int s  = Ysize(f2d);
    const int sh = Xsize(f2d);

    for (int i = 0; i < s; i++) {
        int y, first_allowed_x;

        if (i < sh) {
            y = i;
            first_allowed_x = 0;
        } else {
            y = i - s;
            // x == 0 plane is stored twice in the FFTW format. Don't set it twice in backprojection!
            first_allowed_x = 1;
        }

        // Only iterate over the ellipse in the 2D-image corresponding to the sphere in 3D.
        // Find the x-range inside that ellipse for every given y:
        // |A*v|^2 <= R^2    (for v = (x,y)^t)
        // = v^t A^t A v =: v^t AtA v
        //   <=>
        // (AtA_xx) x^2 + (2 AtA_xy y) x + (AtA_yy y^2 - R^2) <= 0   (quadratic eq. in x)
        //   <=>
        // x in [q - d, q + d],
        // where: q := -AtA_xy y / AtA_xx,
        //        d := sqrt((AtA_xy y)^2 - AtA_xx (AtA_yy y^2 - R^2)) / AtA_xx

        RFLOAT discr = AtA_xy2 * y * y - AtA_xx * (AtA_yy * y * y - max_r2);

        if (discr < 0.0) continue;  // no points inside ellipse for this y

        RFLOAT d = sqrt(discr) / AtA_xx;
        RFLOAT q = - AtA_xy * y / AtA_xx;

        int first_x = ceil(q - d);
        int last_x = floor(q + d);

        if (first_x < first_allowed_x) { first_x = first_allowed_x; }
        if (last_x > sh - 1) { last_x = sh - 1; }

        for (int x = first_x; x <= last_x; x++) {
            // Get the value from the input image
            Complex my_val = direct::elem(f2d, i, x);

            // Get the weight
            RFLOAT my_weight = Mweight ? direct::elem(*Mweight, i, x) : 1.0;

            if (my_weight <= 0.0) continue;

            /*
            In our implementation, (x, y) are not scaled because:

            x_on_ewald = x * r / sqrt(x * x + y * y + r * r)
                       = x / sqrt(1 + (x * x + y * y) / (r * r))
                       ~ x * (1 - (x * x + y * y) / (2 * r * r) + O(1/r^4)) # binomial expansion
                       = x + O(1/r^2)

            same for y_on_ewald

            z_on_ewald = r - r * r / sqrt(x * x + y * y + r * r)
                       ~ r - r * (1 - (x * x + y * y) / (2 * r * r) + O(1/r^4)) # binomial expansion
                       = (x * x + y * y) / (2 * r) + O(1/r^3)

            The error is < 0.0005 reciprocal voxel even for extreme cases
            like 200kV, 1500 A particle, 1 A / pix.
            */

            // Get logical coordinates in the 3D map.
            // Make sure that the Ewald sphere is spherical even under anisotropic mag
            // by first undistorting (x,y) to obtain the true frequencies (xu,yu)

            RFLOAT xu = m00 * x + m01 * y;
            RFLOAT yu = m10 * x + m11 * y;

            RFLOAT z_on_ewaldp = inv_diam_ewald * (xu * xu + yu * yu);

            RFLOAT xp = Ainv(0, 0) * xu + Ainv(0, 1) * yu + Ainv(0, 2) * z_on_ewaldp;
            RFLOAT yp = Ainv(1, 0) * xu + Ainv(1, 1) * yu + Ainv(1, 2) * z_on_ewaldp;
            RFLOAT zp = Ainv(2, 0) * xu + Ainv(2, 1) * yu + Ainv(2, 2) * z_on_ewaldp;

            double r2_3D = xp * xp + yp * yp + zp * zp;

            // redundant:
            if (r2_3D > max_r2) continue;

            if (interpolator == TRILINEAR || r2_3D < min_r2_nn) {
                bool is_neg_x = xp < 0;

                // Only asymmetric half is stored
                if (is_neg_x) {
                    // Get complex conjugated hermitian symmetry pair
                    xp = -xp;
                    yp = -yp;
                    zp = -zp;
                }

                // Trilinear interpolation (with physical coords)
                // Subtract Yinit and Zinit to accelerate access to data (Xinit=0)
                // In that way use direct::elem, rather than elem
                int x0 = floor(xp);
                RFLOAT fx = xp - x0;
                int x1 = x0 + 1;

                int y0 = floor(yp);
                RFLOAT fy = yp - y0;
                y0 -= Yinit(data);
                int y1 = y0 + 1;

                int z0 = floor(zp);
                RFLOAT fz = zp - z0;
                z0 -= Zinit(data);
                int z1 = z0 + 1;

                if (
                    x0 < 0 || x0 + 1 >= data.xdim ||
                    y0 < 0 || y0 + 1 >= data.ydim ||
                    z0 < 0 || z0 + 1 >= data.zdim
                ) continue;

                RFLOAT mfx = 1.0 - fx;
                RFLOAT mfy = 1.0 - fy;
                RFLOAT mfz = 1.0 - fz;

                RFLOAT dd000 = mfz * mfy * mfx;
                RFLOAT dd001 = mfz * mfy *  fx;
                RFLOAT dd010 = mfz *  fy * mfx;
                RFLOAT dd011 = mfz *  fy *  fx;
                RFLOAT dd100 =  fz * mfy * mfx;
                RFLOAT dd101 =  fz * mfy *  fx;
                RFLOAT dd110 =  fz *  fy * mfx;
                RFLOAT dd111 =  fz *  fy *  fx;

                if (is_neg_x) { my_val = conj(my_val); }

                // Store slice in 3D weighted sum
                direct::elem(data, x0, y0, z0) += dd000 * my_val;
                direct::elem(data, x1, y0, z0) += dd001 * my_val;
                direct::elem(data, x0, y1, z0) += dd010 * my_val;
                direct::elem(data, x1, y1, z0) += dd011 * my_val;
                direct::elem(data, x0, y0, z1) += dd100 * my_val;
                direct::elem(data, x1, y0, z1) += dd101 * my_val;
                direct::elem(data, x0, y1, z1) += dd110 * my_val;
                direct::elem(data, x1, y1, z1) += dd111 * my_val;
                // Store corresponding weights
                direct::elem(weight, x0, y0, z0) += dd000 * my_weight;
                direct::elem(weight, x1, y0, z0) += dd001 * my_weight;
                direct::elem(weight, x0, y1, z0) += dd010 * my_weight;
                direct::elem(weight, x1, y1, z0) += dd011 * my_weight;
                direct::elem(weight, x0, y0, z1) += dd100 * my_weight;
                direct::elem(weight, x1, y0, z1) += dd101 * my_weight;
                direct::elem(weight, x0, y1, z1) += dd110 * my_weight;
                direct::elem(weight, x1, y1, z1) += dd111 * my_weight;

            } else if (interpolator == NEAREST_NEIGHBOUR) {
                int x0 = round(xp);
                int y0 = round(yp);
                int z0 = round(zp);

                bool is_neg_x = x0 < 0;

                if (is_neg_x) {
                    // Get complex conjugated hermitian symmetry pair
                    x0 = -x0;
                    y0 = -y0;
                    z0 = -z0;
                }

                const int xr = x0 - Xinit(data);
                const int yr = y0 - Yinit(data);
                const int zr = z0 - Zinit(data);

                if (
                    xr < 0 || xr >= data.xdim ||
                    yr < 0 || yr >= data.ydim ||
                    zr < 0 || zr >= data.zdim
                ) continue;

                if (is_neg_x) {
                    direct::elem(data,   xr, yr, zr) += conj(my_val);
                    direct::elem(weight, xr, yr, zr) += my_weight;
                } else {
                    direct::elem(data,   xr, yr, zr) += my_val;
                    direct::elem(weight, xr, yr, zr) += my_weight;
                }
            } else {
                REPORT_ERROR("FourierInterpolator::backproject%%ERROR: unrecognized interpolator ");
            }
        }
    }
}

void BackProjector::backproject1Dto2D(
    const MultidimArray<Complex> &f1d,
    const Matrix<RFLOAT> &A,
    const MultidimArray<RFLOAT> *Mweight
) {
    Matrix<RFLOAT> Ainv = A.inv();
    Ainv *= (RFLOAT) padding_factor;  // take scaling into account directly

    const int r_max_src = Xsize(f1d) - 1;

    const int r_max_ref   = r_max * padding_factor;
    const int r_max_ref_2 = r_max_ref * r_max_ref;

    // currently not used for some reason
    //const int r_min_NN_ref_2 = r_min_nn * r_min_nn * padding_factor * padding_factor;

    for (int x = 0; x <= r_max_src; x++) {
        RFLOAT my_weight = Mweight ? direct::elem(*Mweight, x) : 1.0;
        if (my_weight <= 0.0) continue;

        Complex my_val = direct::elem(f1d, x);

        // Get logical coordinates in the 3D map
        RFLOAT xp = Ainv(0, 0) * x;
        RFLOAT yp = Ainv(1, 0) * x;

        const RFLOAT r_ref_2 = xp * xp + yp * yp;

        if (r_ref_2 > r_max_ref_2) continue;

        if (interpolator == TRILINEAR /* && r_ref_2 < r_min_NN_ref_2*/) {
            // Only asymmetric half is stored
            const bool is_neg_x = xp < 0;

            if (is_neg_x) {
                // Get complex conjugated hermitian symmetry pair
                xp = -xp;
                yp = -yp;
            }

            // Trilinear interpolation (with physical coords)
            // Subtract Yinit to accelerate access to data (Xinit=0)
            // In that way use direct::elem, rather than A2D_ELEM
            const int x0 = floor(xp);
            const RFLOAT fx = xp - x0;
            const int x1 = x0 + 1;

            int y0 = floor(yp);
            const RFLOAT fy = yp - y0;
            y0 -= Yinit(data);
            const int y1 = y0 + 1;

            const RFLOAT mfx = 1.0 - fx;
            const RFLOAT mfy = 1.0 - fy;

            const RFLOAT dd00 = mfy * mfx;
            const RFLOAT dd01 = mfy *  fx;
            const RFLOAT dd10 =  fy * mfx;
            const RFLOAT dd11 =  fy *  fx;

            if (is_neg_x) { my_val = conj(my_val); }

            // Store slice in 3D weighted sum
            direct::elem(data, x0, y0) += dd00 * my_val;
            direct::elem(data, x1, y0) += dd01 * my_val;
            direct::elem(data, x0, y1) += dd10 * my_val;
            direct::elem(data, x1, y1) += dd11 * my_val;

            // Store corresponding weights
            direct::elem(weight, x0, y0) += dd00 * my_weight;
            direct::elem(weight, x1, y0) += dd01 * my_weight;
            direct::elem(weight, x0, y1) += dd10 * my_weight;
            direct::elem(weight, x1, y1) += dd11 * my_weight;

        } else if (interpolator == NEAREST_NEIGHBOUR ) {
            const int x0 = round(xp);
            const int y0 = round(yp);

            if (x0 < 0) {
                data.elem(-x0, -y0) += conj(my_val);
                weight.elem(-x0, -y0) += my_weight;
            } else {
                data.elem(x0, y0) += my_val;
                weight.elem(x0, y0) += my_weight;
            }
        } else {
            REPORT_ERROR("FourierInterpolator::backproject1Dto2D%%ERROR: unrecognized interpolator ");
        }
    }
}

void BackProjector::backrotate2D(
    const MultidimArray<Complex> &f2d,
    const Matrix<RFLOAT> &A,
    const MultidimArray<RFLOAT> *Mweight,
    Matrix<RFLOAT> *magMatrix
) {
    Matrix<RFLOAT> Ainv = A.inv();
    Ainv *= (RFLOAT) padding_factor;  // take scaling into account directly

    RFLOAT m00, m10, m01, m11;
    fillm(magMatrix, m00, m10, m01, m11);

    const int r_max_ref   = r_max * padding_factor;
    const int r_max_ref_2 = r_max_ref * r_max_ref;

    int min_r2_nn = r_min_nn * r_min_nn * padding_factor * padding_factor;

    // precalculated coefficients for ellipse determination (see further down)

    // Make sure A contains 2D distortion (lowercase 2D, uppercase 3D):
    // Matrix-multiply Ainv(3×2) and m(2×2) to get Am
    const RFLOAT Am_Xx = Ainv(0, 0) * m00 + Ainv(0, 1) * m10;
    const RFLOAT Am_Xy = Ainv(0, 0) * m01 + Ainv(0, 1) * m11;
    const RFLOAT Am_Yx = Ainv(1, 0) * m00 + Ainv(1, 1) * m10;
    const RFLOAT Am_Yy = Ainv(1, 0) * m01 + Ainv(1, 1) * m11;

    // Precompute AmTAm into ATA:
    const RFLOAT AtA_xx = Am_Xx * Am_Xx + Am_Yx * Am_Yx;
    const RFLOAT AtA_xy = Am_Xx * Am_Xy + Am_Yx * Am_Yy;
    const RFLOAT AtA_yy = Am_Xy * Am_Xy + Am_Yy * Am_Yy;
    const RFLOAT AtA_xy2 = AtA_xy * AtA_xy;

    // #define DEBUG_BACKROTATE
    #ifdef DEBUG_BACKROTATE
    std::cerr << " Xsize(f2d)= "<< Xsize(f2d) << std::endl;
    std::cerr << " Ysize(f2d)= "<< Ysize(f2d) << std::endl;
    std::cerr << " Xsize(data)= "<< Xsize(data) << std::endl;
    std::cerr << " Ysize(data)= "<< Ysize(data) << std::endl;
    std::cerr << " Xinit(data)= "<< Xinit(data) << std::endl;
    std::cerr << " Yinit(data)= "<< Yinit(data) << std::endl;
    std::cerr << " Zinit(data)= "<< Zinit(data) << std::endl;
    std::cerr << " max_r= "<< r_max << std::endl;
    std::cerr << " Ainv= " << Ainv << std::endl;
    #endif
    const int s  = Ysize(f2d);
    const int sh = Xsize(f2d);

    for (int i = 0; i < s; i++) {
        int y, first_allowed_x;

        if (i < sh) {
            y = i;
            first_allowed_x = 0;
        } else {
            y = i - s;
            // x == 0 plane is stored twice in the FFTW format. Don't set it twice in backprojection!
            first_allowed_x = 1;
        }

        // Only iterate over the ellipse in the 2D-image corresponding to the sphere in 3D.
        // Find the x-range inside that ellipse for every given y:
        // |A*v|^2 <= R^2    (for v = (x,y)^t)
        // = v^t A^t A v =: v^t AtA v
        //   <=>
        // (AtA_xx) x^2 + (2 AtA_xy y) x + (AtA_yy y^2 - R^2) <= 0   (quadratic eq. in x)
        //   <=>
        // x in [q - d, q + d],
        // where: q := -AtA_xy y / AtA_xx,
        //        d := sqrt((AtA_xy y)^2 - AtA_xx (AtA_yy y^2 - R^2)) / AtA_xx

        RFLOAT discr = AtA_xy2 * y * y - AtA_xx * (AtA_yy * y * y - r_max_ref_2);

        if (discr < 0.0) continue;  // no points inside ellipse for this y

        RFLOAT d = sqrt(discr) / AtA_xx;
        RFLOAT q = - AtA_xy * y / AtA_xx;

        int first_x = ceil(q - d);
        int last_x = floor(q + d);

        if (first_x < first_allowed_x) { first_x = first_allowed_x; }
        if (last_x > sh - 1) { last_x = sh - 1; }

        for (int x = first_x; x <= last_x; x++) {

            RFLOAT my_weight = Mweight ? direct::elem(*Mweight, x, i) : 1.0;
            if (my_weight <= 0.0) continue;

            // Get the relevant value in the input image
            Complex my_val = direct::elem(f2d, x, i);

            // Get logical coordinates in the 3D map
            RFLOAT xu = m00 * x + m01 * y;
            RFLOAT yu = m10 * x + m11 * y;

            RFLOAT xp = Ainv(0, 0) * xu + Ainv(0, 1) * yu;
            RFLOAT yp = Ainv(1, 0) * xu + Ainv(1, 1) * yu;

            RFLOAT r_ref_2 = xp * xp + yp * yp;

            if (interpolator == TRILINEAR || r_ref_2 < min_r2_nn) {
                const bool is_neg_x = xp < 0;
                // Only asymmetric half is stored
                if (is_neg_x) {
                    // Get complex conjugated hermitian symmetry pair
                    xp = -xp;
                    yp = -yp;
                }

                // Trilinear interpolation (with physical coords)
                // Subtract Yinit to accelerate access to data (Xinit=0)
                // In that way use direct::elem, rather than A2D_ELEM
                const int x0 = floor(xp);
                const RFLOAT fx = xp - x0;
                const int x1 = x0 + 1;

                int y0 = floor(yp);
                const RFLOAT fy = yp - y0;
                y0 -= Yinit(data);
                const int y1 = y0 + 1;

                const RFLOAT mfx = 1.0 - fx;
                const RFLOAT mfy = 1.0 - fy;

                const RFLOAT dd00 = mfy * mfx;
                const RFLOAT dd01 = mfy *  fx;
                const RFLOAT dd10 =  fy * mfx;
                const RFLOAT dd11 =  fy *  fx;

                if (is_neg_x) { my_val = conj(my_val); }

                // Store slice in 3D weighted sum
                direct::elem(data, x0, y0) += dd00 * my_val;
                direct::elem(data, x1, y0) += dd01 * my_val;
                direct::elem(data, x0, y1) += dd10 * my_val;
                direct::elem(data, x1, y1) += dd11 * my_val;

                // Store corresponding weights
                direct::elem(weight, x0, y0) += dd00 * my_weight;
                direct::elem(weight, x1, y0) += dd01 * my_weight;
                direct::elem(weight, x0, y1) += dd10 * my_weight;
                direct::elem(weight, x1, y1) += dd11 * my_weight;

            } else if (interpolator == NEAREST_NEIGHBOUR) {
                const int x0 = round(xp);
                const int y0 = round(yp);

                if (x0 < 0) {
                      data.elem(-x0, -y0) += conj(my_val);
                    weight.elem(-x0, -y0) += my_weight;
                } else {
                      data.elem(x0, y0) += my_val;
                    weight.elem(x0, y0) += my_weight;
                }
            } else {
                REPORT_ERROR("FourierInterpolator::backrotate2D%%ERROR: unrecognized interpolator ");
            }
        }
    }
}

void BackProjector::backrotate3D(
    const MultidimArray<Complex> &f3d,
    const Matrix<RFLOAT> &A,
    const MultidimArray<RFLOAT> *Mweight
) {
    // f3d should already be in the right size (ori_size,orihalfdim)
    // AND the points outside max_r should already be zero.

    Matrix<RFLOAT> Ainv = A.inv();
    Ainv *= (RFLOAT)padding_factor;  // take scaling into account directly

    const int r_max_src = Xsize(f3d) - 1;
    const int r_max_src_2 = r_max_src * r_max_src;

    const int r_max_ref = r_max * padding_factor;
    const int r_max_ref_2 = r_max_ref * r_max_ref;

    const int r_min_NN_ref_2 = r_min_nn * r_min_nn * padding_factor * padding_factor;

    // #define DEBUG_BACKROTATE
    #ifdef DEBUG_BACKROTATE
    std::cerr << " Xsize(f3d)= "<< Xsize(f3d) << std::endl;
    std::cerr << " Ysize(f3d)= "<< Ysize(f3d) << std::endl;
    std::cerr << " Xsize(data)= "<< Xsize(data) << std::endl;
    std::cerr << " Ysize(data)= "<< Ysize(data) << std::endl;
    std::cerr << " Xinit(data)= "<< Xinit(data) << std::endl;
    std::cerr << " Yinit(data)= "<< Yinit(data) << std::endl;
    std::cerr << " Zinit(data)= "<< Zinit(data) << std::endl;
    std::cerr << " max_r= "<< r_max << std::endl;
    std::cerr << " Ainv= " << Ainv << std::endl;
    #endif

    for (int k = 0; k < Zsize(f3d); k++) {
        int z, x_min;

        // Don't search beyond square with side max_r
        if (k <= r_max_src) {
            z = k;
            x_min = 0;
        } else {
            z = k - Zsize(f3d);
            /// TODO: still check this better in the 3D case!!!
            // x==0 (y,z)-plane is stored twice in the FFTW format. Don't set it twice in BACKPROJECTION!
            x_min = 1;
        }

        int z2 = z * z;

        for (int i = 0; i < Ysize(f3d); i++) {
            int y = i <= r_max_src ? i : i - Ysize(f3d);
            int y2 = y * y;

            const RFLOAT yz2 = y2 + z2;

            // avoid negative square root
            if (yz2 > r_max_src_2) continue;

            const int x_max = floor(sqrt(r_max_src_2 - yz2));

            for (int x = x_min; x <= x_max; x++) {
                // Get logical coordinates in the 3D map
                RFLOAT xp = Ainv(0, 0) * x + Ainv(0, 1) * y + Ainv(0, 2) * z;
                RFLOAT yp = Ainv(1, 0) * x + Ainv(1, 1) * y + Ainv(1, 2) * z;
                RFLOAT zp = Ainv(2, 0) * x + Ainv(2, 1) * y + Ainv(2, 2) * z;

                const int r_ref_2 = xp * xp + yp * yp + zp * zp;

                if (r_ref_2 > r_max_ref_2) continue;

                // Get the weight
                RFLOAT my_weight = Mweight ? direct::elem(*Mweight, x, i, k) : 1.0;
                if (my_weight <= 0.0) continue;

                Complex my_val = direct::elem(f3d, x, i, k);

                if (interpolator == TRILINEAR || r_ref_2 < r_min_NN_ref_2) {
                    // Only asymmetric half is stored
                    bool is_neg_x = xp < 0;

                    if (is_neg_x) {
                        // Get complex conjugated hermitian symmetry pair
                        xp = -xp;
                        yp = -yp;
                        zp = -zp;
                    }

                    // Trilinear interpolation (with physical coords)
                    // Subtract Yinit to accelerate access to data (Xinit=0)
                    // In that way use direct::elem, rather than elem
                    const int x0 = floor(xp);
                    const RFLOAT fx = xp - x0;
                    const int x1 = x0 + 1;

                    int y0 = floor(yp);
                    const RFLOAT fy = yp - y0;
                    y0 -= Yinit(data);
                    const int y1 = y0 + 1;

                    int z0 = floor(zp);
                    const RFLOAT fz = zp - z0;
                    z0 -= Zinit(data);
                    const int z1 = z0 + 1;

                    const RFLOAT mfx = 1.0 - fx;
                    const RFLOAT mfy = 1.0 - fy;
                    const RFLOAT mfz = 1.0 - fz;

                    const RFLOAT dd000 = mfz * mfy * mfx;
                    const RFLOAT dd001 = mfz * mfy *  fx;
                    const RFLOAT dd010 = mfz *  fy * mfx;
                    const RFLOAT dd011 = mfz *  fy *  fx;
                    const RFLOAT dd100 =  fz * mfy * mfx;
                    const RFLOAT dd101 =  fz * mfy *  fx;
                    const RFLOAT dd110 =  fz *  fy * mfx;
                    const RFLOAT dd111 =  fz *  fy *  fx;

                    if (is_neg_x) {
                        my_val = conj(my_val);
                    }

                    // Store slice in 3D weighted sum
                    direct::elem(data, x0, y0, z0) += dd000 * my_val;
                    direct::elem(data, x1, y0, z0) += dd001 * my_val;
                    direct::elem(data, x0, y1, z0) += dd010 * my_val;
                    direct::elem(data, x1, y1, z0) += dd011 * my_val;
                    direct::elem(data, x0, y0, z1) += dd100 * my_val;
                    direct::elem(data, x1, y0, z1) += dd101 * my_val;
                    direct::elem(data, x0, y1, z1) += dd110 * my_val;
                    direct::elem(data, x1, y1, z1) += dd111 * my_val;

                    // Store corresponding weights
                    direct::elem(weight, x0, y0, z0) += dd000 * my_weight;
                    direct::elem(weight, x1, y0, z0) += dd001 * my_weight;
                    direct::elem(weight, x0, y1, z0) += dd010 * my_weight;
                    direct::elem(weight, x1, y1, z0) += dd011 * my_weight;
                    direct::elem(weight, x0, y0, z1) += dd100 * my_weight;
                    direct::elem(weight, x1, y0, z1) += dd101 * my_weight;
                    direct::elem(weight, x0, y1, z1) += dd110 * my_weight;
                    direct::elem(weight, x1, y1, z1) += dd111 * my_weight;

                } else if (interpolator == NEAREST_NEIGHBOUR) {
                    const int x0 = round(xp);
                    const int y0 = round(yp);
                    const int z0 = round(zp);

                    if (x0 < 0) {
                        data.elem(-x0, -y0, -z0) += conj(my_val);
                        weight.elem(-x0, -y0, -z0) += my_weight;
                    } else {
                        data.elem(x0, y0, z0) += my_val;
                        weight.elem(x0, y0, z0) += my_weight;
                    }
                } else {
                    REPORT_ERROR("BackProjector::backrotate3D%%ERROR: unrecognized interpolator ");
                }
            }
        }
    }
}

std::pair<MultidimArray<Complex>, MultidimArray<RFLOAT>> BackProjector::getLowResDataAndWeight(int lowres_r_max) {

    const int lowres_r2_max = round(padding_factor * lowres_r_max) * round(padding_factor * lowres_r_max);
    const int lowres_pad_size = 2 * round(padding_factor * lowres_r_max) + 3;

    // Check lowres_r_max is not too big
    if (lowres_r_max > r_max)
        REPORT_ERROR("BackProjector::getLowResDataAndWeight%%ERROR: lowres_r_max is bigger than r_max");

    // Initialize lowres_data and low_res_weight arrays
    MultidimArray<Complex> lowres_data;
    MultidimArray<RFLOAT>  lowres_weight;

    if (ref_dim == 2) {
        lowres_data  .resize(lowres_pad_size, lowres_pad_size / 2 + 1);
        lowres_weight.resize(lowres_pad_size, lowres_pad_size / 2 + 1);
    } else {
        lowres_data  .resize(lowres_pad_size, lowres_pad_size, lowres_pad_size / 2 + 1);
        lowres_weight.resize(lowres_pad_size, lowres_pad_size, lowres_pad_size / 2 + 1);
    }
    lowres_data.setXmippOrigin().xinit = 0;
    lowres_weight.setXmippOrigin().xinit = 0;

    // fill lowres arrays with relevant values
    FOR_ALL_ELEMENTS_IN_ARRAY3D(lowres_data, i, j, k) {
        if (hypot2(i, j, k) <= lowres_r2_max) {
            lowres_data  .elem(i, j, k) = data  .elem(i, j, k);
            lowres_weight.elem(i, j, k) = weight.elem(i, j, k);
        }
    }
    return {lowres_data, lowres_weight};
}

void BackProjector::setLowResDataAndWeight(
    MultidimArray<Complex> &lowres_data, MultidimArray<RFLOAT> &lowres_weight,
    int lowres_r_max
) {

    const int lowres_r2_max = round(padding_factor * lowres_r_max) * round(padding_factor * lowres_r_max);
    const int lowres_pad_size = 2 * round(padding_factor * lowres_r_max) + 3;

    // Check lowres_r_max is not too big
    if (lowres_r_max > r_max)
        REPORT_ERROR("BackProjector::getLowResDataAndWeight%%ERROR: lowres_r_max is bigger than r_max");

    // Check sizes of lowres_data and lowres_weight
    if (
        Ysize(lowres_data) != lowres_pad_size ||
        Xsize(lowres_data) != lowres_pad_size / 2 + 1 ||
        Zsize(lowres_data) != lowres_pad_size && ref_dim == 3
    )
        REPORT_ERROR("BackProjector::setLowResDataAndWeight%%ERROR: lowres_data is not of expected size...");
    if (
        Ysize(lowres_weight) != lowres_pad_size ||
        Xsize(lowres_weight) != lowres_pad_size / 2 + 1 ||
        Zsize(lowres_weight) != lowres_pad_size && ref_dim == 3
    )
        REPORT_ERROR("BackProjector::setLowResDataAndWeight%%ERROR: lowres_weight is not of expected size...");

    // Re-set origin to the expected place
    lowres_data.setXmippOrigin().xinit = 0;
    lowres_weight.setXmippOrigin().xinit = 0;

    // Overwrite data and weight with the lowres arrays
    FOR_ALL_ELEMENTS_IN_ARRAY3D(lowres_data, i, j, k) {
        if (hypot2(i, j, k) <= lowres_r2_max) {
            data  .elem(i, j, k) = lowres_data  .elem(i, j, k);
            weight.elem(i, j, k) = lowres_weight.elem(i, j, k);
        }
    }
}

MultidimArray<Complex> BackProjector::getDownsampledAverage(bool divide) const {

    MultidimArray<Complex> avg;

    // Pre-set down_data and down_weight sizes
    const int down_size = 2 * (r_max + 1) + 1;

    // Short side of data array
    switch (ref_dim) {

        case 2:
        avg.initZeros(down_size, down_size / 2 + 1);
        break;

        case 3:
        avg.initZeros(down_size, down_size, down_size / 2 + 1);
        break;

        default:
        REPORT_ERROR("BackProjector::getDownsampledAverage%%ERROR: Dimension of the data array should be 2 or 3");
    }
    // Set origin in the y.z-center, but on the left side for x.
    avg.setXmippOrigin().xinit = 0;
    auto down_weight = MultidimArray<RFLOAT>::zeros(avg);
    // down_weight same size as down_data

    // Now calculate the down-sized sum
    FOR_ALL_ELEMENTS_IN_ARRAY3D(data, i, j, k) {
        const int ip = round((RFLOAT) i / padding_factor);
        const int jp = round((RFLOAT) j / padding_factor);
        const int kp = round((RFLOAT) k / padding_factor);

        // TMP
        // #define CHECK_SIZE
        #ifdef CHECK_SIZE
        if (
            ip < Yinit(avg) || ip > Ylast(avg)
            jp < Xinit(avg) || jp > Xlast(avg) ||
            kp < Zinit(avg) || kp > Zlast(avg) ||
        ) {
            std::cerr << " kp= " << kp << " ip= " << ip << " jp= " << jp << std::endl;
            avg.printShape();
            REPORT_ERROR("BackProjector::getDownsampledAverage: indices out of range");
        }
        #endif
        avg.elem(ip, jp, kp) += data.elem(i, j, k);
        down_weight.elem(ip, jp, kp) += divide ? weight.elem(i, j, k) : 1.0;
    }

    // Calculate the straightforward average in the downsampled arrays
    for (long int n = 0; n < avg.size(); n++)
        avg[n] = down_weight[n] == 0.0 ? Complex(0.0) :
        avg[n] / down_weight[n];

    return avg;
}

MultidimArray<RFLOAT> BackProjector::calculateDownSampledFourierShellCorrelation(
    const MultidimArray<Complex> &avg1, const MultidimArray<Complex> &avg2
) const {
    if (!avg1.sameShape(avg2))
        REPORT_ERROR("ERROR BackProjector::calculateDownSampledFourierShellCorrelation: two arrays have different sizes");

    const long int n = ori_size / 2 + 1;
    auto num  = MultidimArray<RFLOAT>::zeros(n);
    auto den1 = MultidimArray<RFLOAT>::zeros(n);
    auto den2 = MultidimArray<RFLOAT>::zeros(n);
    auto fsc  = MultidimArray<RFLOAT>::zeros(n);

    FOR_ALL_ELEMENTS_IN_ARRAY3D(avg1, i, j, k) {
        const RFLOAT R = hypot((double) i, j, k);

        if (R > r_max) continue;

        const int idx = round(R);

        const Complex z1 = avg1.elem(i, j, k);
        const Complex z2 = avg2.elem(i, j, k);

        num .elem(idx) += dot(z1, z2);
        den1.elem(idx) += z1.norm();
        den2.elem(idx) += z2.norm();
    }

    for (int i = Xinit(fsc); i <= Xlast(fsc); i++) {
        if (den1.elem(i) * den2.elem(i) > 0.0) {
            fsc.elem(i) = num.elem(i) / sqrt(den1.elem(i) * den2.elem(i));
        }
    }

    // Always set zero-resolution shell to FSC=1
    // Raimond Ravelli reported a problem with FSC=1 at res=0 on 13feb2013...
    // (because of a suboptimal normalisation scheme, but anyway)
    fsc.elem(0) = 1.0;
    return fsc;
}

void BackProjector::updateSSNRarrays(
    RFLOAT tau2_fudge,
    MultidimArray<RFLOAT> &tau2_io,
    MultidimArray<RFLOAT> &sigma2_out,
    MultidimArray<RFLOAT> &data_vs_prior_out,
    MultidimArray<RFLOAT> &fourier_coverage_out,
    const MultidimArray<RFLOAT>& fsc,
    bool update_tau2_with_fsc,
    bool iswhole
) {
    // never rely on references (handed to you from the outside) for computation:
    // they could be the same (i.e. reconstruct(..., dummy, dummy, dummy, dummy, ...); )
    MultidimArray<RFLOAT> tau2 = tau2_io;
    const int max_r2 = round(r_max * padding_factor) * round(r_max * padding_factor);
    RFLOAT oversampling_correction = ref_dim == 3 ? padding_factor * padding_factor * padding_factor : padding_factor * padding_factor;

    // First calculate the radial average of the (inverse of the) power of the noise in the reconstruction
    // This is the left-hand side term in the nominator of the Wiener-filter-like update formula
    // and it is stored inside the weight vector
    // Then, if (do_map) add the inverse of tau2-spectrum values to the weight
    MultidimArray<RFLOAT> sigma2  = MultidimArray<RFLOAT>::zeros(ori_size / 2 + 1);
    MultidimArray<RFLOAT> counter = MultidimArray<RFLOAT>::zeros(ori_size / 2 + 1);
    FOR_ALL_ELEMENTS_IN_ARRAY3D(weight, i, j, k) {
        const RFLOAT r2 = hypot2(i, j, k);
        if (r2 < max_r2) {
            int ires = round(sqrt(r2) / padding_factor);
            RFLOAT invw = oversampling_correction * weight.elem(i, j, k);
            direct::elem(sigma2,  ires) += invw;
            direct::elem(counter, ires) += 1.0;
        }
    }

    // Average (inverse of) sigma2 in reconstruction
    for (long int i = 0; i < Xsize(sigma2); i++) {
        double x = direct::elem(sigma2, i);
        if (x > 1e-10) {
            x = direct::elem(counter, i) / x;
        } else if (x == 0) {
            x = 0.0;
        } else {
            std::cerr << " direct::elem(sigma2, i)= " << x << std::endl;
            REPORT_ERROR("BackProjector::reconstruct: ERROR: unexpectedly small, yet non-zero sigma2 value, this should not happen...a");
        }
    }

    tau2.reshape(ori_size / 2 + 1);
    MultidimArray<RFLOAT> data_vs_prior    = MultidimArray<RFLOAT>::zeros(ori_size / 2 + 1);
    MultidimArray<RFLOAT> fourier_coverage = MultidimArray<RFLOAT>::zeros(ori_size / 2 + 1);
    counter.initZeros(ori_size / 2 + 1);
    if (update_tau2_with_fsc) {
        // Then calculate new tau2 values, based on the FSC
        if (!fsc.sameShape(sigma2) || !fsc.sameShape(tau2)) {
            fsc.printShape(std::cerr);
            tau2.printShape(std::cerr);
            sigma2.printShape(std::cerr);
            REPORT_ERROR("ERROR BackProjector::reconstruct: sigma2, tau2 and fsc have different sizes");
        }
        for (long int i = 0; i < Xsize(sigma2); i++) {
            // FSC cannot be negative or zero for conversion into tau2
            RFLOAT myfsc = std::max(0.001, direct::elem(fsc, i));
            if (iswhole) {
                // Factor two because of twice as many particles
                // Sqrt-term to get 60-degree phase errors....
                myfsc = sqrt(2.0 * myfsc / (myfsc + 1.0));
            }
            myfsc = std::min(0.999, myfsc);
            RFLOAT myssnr = myfsc / (1.0 - myfsc) * tau2_fudge;
            // Sjors 29nov2017 try tau2_fudge for pulling harder on Refine3D runs...
            RFLOAT fsc_based_tau = myssnr * direct::elem(sigma2, i);
            direct::elem(tau2, i) = fsc_based_tau;
            // data_vs_prior is merely for reporting: it is not used for anything in the reconstruction
            direct::elem(data_vs_prior, i) = myssnr;
        }
    }

    // Now accumulate data_vs_prior if (!update_tau2_with_fsc)
    // Also accumulate fourier_coverage
    FOR_ALL_ELEMENTS_IN_ARRAY3D(weight, i, j, k) {
        int r2 = hypot2(i, j, k);
        if (r2 < max_r2) {
            int ires = round(sqrt((RFLOAT) r2) / padding_factor);
            RFLOAT invw = weight.elem(i, j, k);

            RFLOAT invtau2;
            if (direct::elem(tau2, ires) > 0.0) {
                // Calculate inverse of tau2
                invtau2 = 1.0 / (oversampling_correction * tau2_fudge * direct::elem(tau2, ires));
            } else if (direct::elem(tau2, ires) == 0.0) {
                // If tau2 is zero, use small value instead
                invtau2 = 1000.0 / invw;
            } else {
                std::cerr << " sigma2= " << sigma2 << std::endl;
                std::cerr << " fsc= " << fsc << std::endl;
                std::cerr << " tau2= " << tau2 << std::endl;
                REPORT_ERROR("ERROR BackProjector::reconstruct: Negative or zero values encountered for tau2 spectrum!");
            }

            // Keep track of spectral evidence-to-prior ratio and remaining noise in the reconstruction
            if (!update_tau2_with_fsc) {
                direct::elem(data_vs_prior, ires) += invw / invtau2;
            }

            // Keep track of the coverage in Fourier space
            if (invw / invtau2 >= 1.0) {
                direct::elem(fourier_coverage, ires) += 1.0;
            }

            direct::elem(counter, ires) += 1.0;

        }
    }

    // Average data_vs_prior
    if (!update_tau2_with_fsc) {
        for (long int i = 0; i < data_vs_prior.xdim; i++) {
            RFLOAT x = direct::elem(data_vs_prior, i);
            RFLOAT n = direct::elem(counter, i);
            if (i > r_max) {
                x = 0.0;
            } else if (n < 0.001) {
                x = 999.0;
            } else {
                x /= n;
            }
        }
    }

    // Calculate Fourier coverage in each shell
    for (long int i = 0; i < Xsize(fourier_coverage); i++) {
        if (direct::elem(counter, i) > 0.0)
            direct::elem(fourier_coverage, i) /= direct::elem(counter, i);
    }

    // Send back the output
    tau2_io = tau2;
    sigma2_out = sigma2;
    data_vs_prior_out = data_vs_prior;
    fourier_coverage_out = fourier_coverage;
}

MultidimArray<RFLOAT> BackProjector::externalReconstruct(
    const FileName &fn_out,
    MultidimArray<RFLOAT> &fsc_halves_io,
    MultidimArray<RFLOAT> &tau2_io,
    MultidimArray<RFLOAT> &sigma2_ref,
    MultidimArray<RFLOAT> &data_vs_prior,
    bool iswhole,
    RFLOAT tau2_fudge,
    int verb
) {

    MultidimArray<RFLOAT> fsc_halves = fsc_halves_io;
    MultidimArray<RFLOAT> tau2 = tau2_io;

    const int max_r2 = round(r_max * padding_factor) * round(r_max * padding_factor);
    int padoridim = round(padding_factor * ori_size);

    // Write out data array
    auto data_decentered = Projector::decentered(data, max_r2,
        pad_size / 2 + 1, pad_size, ref_dim == 3 ? pad_size : 1);
    windowFourierTransform(data_decentered, padoridim);
    ComplexIO::write(data_decentered, fn_out + "_external_reconstruct_data", ".mrc");

    // Write out weight array
    auto weight_decentered = Projector::decentered(weight, max_r2,
        pad_size / 2 + 1, pad_size, ref_dim == 3 ? pad_size : 1);
    windowFourierTransform(weight_decentered, padoridim);
    Image<RFLOAT>(weight_decentered).write(fn_out + "_external_reconstruct_weight.mrc");

    // Write out STAR file for input to external reconstruction program
    MetaDataTable MDlist;
    MDlist.name = "external_reconstruct_general";
    MDlist.isList = true;
    const long int i = MDlist.addObject();
    const FileName fn_recons   = fn_out + "_external_reconstruct.mrc";
    const FileName fn_out_star = fn_out + "_external_reconstruct_out.star";
    MDlist.setValue(EMDL::OPTIMISER_EXTERNAL_RECONS_DATA_REAL, fn_out + "_external_reconstruct_data_real.mrc", i);
    MDlist.setValue(EMDL::OPTIMISER_EXTERNAL_RECONS_DATA_IMAG, fn_out + "_external_reconstruct_data_imag.mrc", i);
    MDlist.setValue(EMDL::OPTIMISER_EXTERNAL_RECONS_WEIGHT,    fn_out + "_external_reconstruct_weight.mrc",    i);
    MDlist.setValue(EMDL::OPTIMISER_EXTERNAL_RECONS_RESULT, fn_recons, i);
    MDlist.setValue(EMDL::OPTIMISER_EXTERNAL_RECONS_NEWSTAR, fn_out_star, i);
    MDlist.setValue(EMDL::MLMODEL_TAU2_FUDGE_FACTOR, tau2_fudge, i);
    MDlist.setValue(EMDL::MLMODEL_PADDING_FACTOR, padding_factor, i);
    MDlist.setValue(EMDL::MLMODEL_DIMENSIONALITY, ref_dim, i);
    MDlist.setValue(EMDL::MLMODEL_ORIGINAL_SIZE, ori_size, i);
    MDlist.setValue(EMDL::MLMODEL_CURRENT_SIZE, 2 * r_max, i);

    MetaDataTable MDtau;
    MDtau.name = "external_reconstruct_tau2";
    for (int i = 0; i < Xsize(tau2); i++) {
        MDtau.addObject();
        MDtau.setValue(EMDL::SPECTRAL_IDX, i, i);
        MDtau.setValue(EMDL::MLMODEL_TAU2_REF, tau2.elem(i), i);
        MDtau.setValue(EMDL::MLMODEL_FSC_HALVES_REF, fsc_halves.elem(i), i);
    }

    const FileName fn_star = fn_out + "_external_reconstruct.star";
    {
        std::ofstream fh(fn_star.c_str(), std::ios::out);
        if (!fh)
            REPORT_ERROR((std::string) "BackProjector::externalReconstruct: Cannot write file: " + fn_star);
        MDlist.write(fh);
        MDtau.write(fh);
    }

    // Make the system call: program name plus the STAR file for the external reconstruction program as its first argument
    const char *my_exec = getenv ("RELION_EXTERNAL_RECONSTRUCT_EXECUTABLE");
    if (!my_exec) { my_exec = DEFAULT_EXTERNAL_RECONSTRUCT; }
    const std::string command = std::string(my_exec) + " " + fn_star;

    if (verb > 0)
        std::cout << std::endl << " + Making system call for external reconstruction: " << command << std::endl;

    if (system(command.c_str())) {
        REPORT_ERROR(" ERROR: there was something wrong with system call: " + command);
    } else if (verb > 0) {
        std::cout << " + External reconstruction finished successfully, reading result back in ... " << std::endl;
    }

    // Read the resulting map back into memory
    MultidimArray<RFLOAT> vol_out =
        Image<RFLOAT>::from_filename(fn_recons).data.setXmippOrigin();

    if (exists(fn_out_star)) {

        MetaDataTable MDnewtau;
        MDnewtau.read(fn_out_star);

        if (!MDnewtau.containsLabel(EMDL::SPECTRAL_IDX))
            REPORT_ERROR("ERROR: external reconstruct output STAR file does not contain spectral idx!");

        // Directly update tau2 spectrum
        for (long int i : MDnewtau) {
            int idx = MDnewtau.getValue<int>(EMDL::SPECTRAL_IDX, i);

            if (idx >= Xsize(tau2_io)) continue;

            if (MDnewtau.containsLabel(EMDL::MLMODEL_TAU2_REF)) {
                tau2_io.elem(idx)       = MDnewtau.getValue<RFLOAT>(EMDL::MLMODEL_TAU2_REF, i);
                data_vs_prior.elem(idx) = tau2_io.elem(idx) / sigma2_ref.elem(idx);
            } else if (MDnewtau.containsLabel(EMDL::POSTPROCESS_FSC_GENERAL)) {
                idx                     = MDnewtau.getValue<int>(EMDL::SPECTRAL_IDX, i);
                fsc_halves_io.elem(idx) = MDnewtau.getValue<RFLOAT>(EMDL::POSTPROCESS_FSC_GENERAL, i);

                RFLOAT myfsc = std::max(0.001, fsc_halves_io.elem(idx));
                if (iswhole) {
                    // Factor two because of twice as many particles
                    // Sqrt-term to get 60-degree phase errors....
                    myfsc = sqrt(2.0 * myfsc / (myfsc + 1.0));
                }
                myfsc = std::min(0.999, myfsc);
                const RFLOAT myssnr = myfsc / (1.0 - myfsc) * tau2_fudge;
                tau2_io.elem(idx) = myssnr * sigma2_ref.elem(idx);
                data_vs_prior.elem(idx) = myssnr;
            } else {
                REPORT_ERROR("ERROR: output STAR file from external reconstruct does not contain tau2 or FSC array");
            }
        }

        if (verb > 0)
            std::cout << " + External reconstruction successfully updated external tau2 array ... " << std::endl;
    }

    return vol_out;
}

MultidimArray<RFLOAT> BackProjector::reconstruct(
    int max_iter_preweight, bool do_map,
    const MultidimArray<RFLOAT> &tau2,
    RFLOAT tau2_fudge, RFLOAT normalise,
    int minres_map, bool printTimes,
    Image<RFLOAT> *weight_out
) {
    #ifdef TIMING
    Timer ReconTimer;
    int ReconS[] = {
        ReconTimer.setNew(" RcS1_Init "),
        ReconTimer.setNew(" RcS2_Shape&Noise "),
        ReconTimer.setNew(" RcS2.5_Regularize "),
        ReconTimer.setNew(" RcS3_skipGridding "),
        ReconTimer.setNew(" RcS4_doGridding_norm "),
        ReconTimer.setNew(" RcS5_doGridding_init "),
        ReconTimer.setNew(" RcS6_doGridding_iter "),
        ReconTimer.setNew(" RcS7_doGridding_apply "),
        ReconTimer.setNew(" RcS8_blobConvolute "),
        ReconTimer.setNew(" RcS9_blobResize "),
        ReconTimer.setNew(" RcS10_blobSetReal "),
        ReconTimer.setNew(" RcS11_blobSetTemp "),
        ReconTimer.setNew(" RcS12_blobTransform "),
        ReconTimer.setNew(" RcS13_blobCenterFFT "),
        ReconTimer.setNew(" RcS14_blobNorm1 "),
        ReconTimer.setNew(" RcS15_blobSoftMask "),
        ReconTimer.setNew(" RcS16_blobNorm2 "),
        ReconTimer.setNew(" RcS17_WindowReal "),
        ReconTimer.setNew(" RcS18_GriddingCorrect "),
        ReconTimer.setNew(" RcS19_tauInit "),      // Unused
        ReconTimer.setNew(" RcS20_tausetReal "),   // Unused
        ReconTimer.setNew(" RcS21_tauTransform "), // Unused
        ReconTimer.setNew(" RcS22_tautauRest "),   // Unused
        ReconTimer.setNew(" RcS23_tauShrinkToFit "),
        ReconTimer.setNew(" RcS24_extra "),        // Unused
    };
    #endif

    RFLOAT oversampling_correction;
    FourierTransformer transformer;
    const int max_r2 = round(r_max * padding_factor) * round(r_max * padding_factor);
    /// TODO: Turn this block into an init function.
    MultidimArray<RFLOAT> vol_out;
    auto &Fconv = [&] () -> MultidimArray<Complex> & {
    ifdefTIMING(TicToc tt (ReconTimer, ReconS[0]);)
    oversampling_correction = ref_dim == 3 ?
        padding_factor * padding_factor * padding_factor :
        padding_factor * padding_factor;

    // #define DEBUG_RECONSTRUCT
    #ifdef DEBUG_RECONSTRUCT
    Image<RFLOAT>(weight).write("reconstruct_initial_weight.spi");
    std::cerr << " pad_size= " << pad_size << " padding_factor= " << padding_factor << " max_r2= " << max_r2 << std::endl;
    #endif

    // Set Fconv to the right size
    if (ref_dim == 2) {
        vol_out.setDimensions(pad_size, pad_size);
    } else {
        // Too costly to actually allocate the space
        // Trick transformer with the right dimensions
        vol_out.setDimensions(pad_size, pad_size, pad_size);
    }

    transformer.setReal(vol_out);  // Fake set real. 1. Allocate space for Fconv 2. calculate plans.
    vol_out.clear();  // Reset dimensions to 0
    return transformer.getFourier();
    }();

    MultidimArray<RFLOAT> Fweight;
    Fweight.reshape(Fconv);
    {
    ifdefTIMING(TicToc tt (ReconTimer, ReconS[1]);)
    // Go from projector-centered to FFTW-uncentered
    Fweight = Projector::decentered(weight, max_r2, pad_size / 2 + 1, pad_size, pad_size);
    }

    {
    ifdefTIMING(TicToc tt (ReconTimer, ReconS[2]);)
    // Apply MAP-additional term to the Fweight array
    // This will regularise the actual reconstruction
    if (do_map) {
        // Then, add the inverse of tau2-spectrum values to the weight
        FOR_ALL_ELEMENTS_IN_FFTW_TRANSFORM(Fconv) {
            int r2 = hypot2(ip, jp, kp);
            if (r2 < max_r2) {
                int ires = round(sqrt((RFLOAT) r2) / padding_factor);
                RFLOAT invw = direct::elem(Fweight, i, j, k);

                RFLOAT invtau2;
                if (direct::elem(tau2, ires) > 0.0) {
                    // Calculate inverse of tau2
                    invtau2 = 1.0 / (oversampling_correction * tau2_fudge * direct::elem(tau2, ires));
                } else if (direct::elem(tau2, ires) < 1e-20) {
                    // If tau2 is zero, use small value instead
                    invtau2 = invw > 1e-20 ? 1.0 / (0.001 * invw) : 0.0;
                } else {
                    std::cerr << " tau2= " << tau2 << std::endl;
                    REPORT_ERROR("ERROR BackProjector::reconstruct: Negative or zero values encountered for tau2 spectrum!");
                }

                // Only for (ires >= minres_map) add Wiener-filter like term
                if (ires >= minres_map) {
                    // Add the inverse-of-tau2_class term
                    // and store the new weight again in Fweight
                    direct::elem(Fweight, i, j, k) = invw += invtau2;
                }
            }
        }
    }
    }

    if (skip_gridding) {
        ifdefTIMING(ReconTimer.tic(ReconS[3]);)  // BUG: No toc!
        Fconv = Projector::decentered(data, max_r2, pad_size / 2 + 1, pad_size, pad_size);

        // Prevent divisions by zero: set Fweight to at least 1/1000th of the radially averaged weight at that resolution
        // beyond r_max, set Fweight to at least 1/1000th of the radially averaged weight at r_max;
        MultidimArray<RFLOAT> radavg_weight(r_max), counter(r_max);
        const int round_max_r2 = round(r_max * padding_factor * r_max * padding_factor);
        FOR_ALL_ELEMENTS_IN_FFTW_TRANSFORM(Fweight) {
            const int r2 = hypot2(ip, jp, kp);
            // Note that (r < ires) != (r2 < max_r2), because max_r2 = round(r_max * padding_factor)^2.
            // We have to use round_max_r2 = round((r_max * padding_factor)^2).
            // e.g. k = 0, i = 7, j = 28, max_r2 = 841, r_max = 16, padding_factor = 18.
            if (r2 < round_max_r2) {
                const int ires = floor(sqrt((RFLOAT) r2) / padding_factor);
                if (ires >= Xsize(radavg_weight)) {
                    std::cerr << " k= " << k << " i= " << i << " j= " << j << std::endl;
                    std::cerr << " ires= " << ires << " Xsize(radavg_weight)= " << Xsize(radavg_weight) << std::endl;
                    REPORT_ERROR("BUG: ires >=Xsize(radavg_weight) ");
                }
                direct::elem(radavg_weight, ires) += direct::elem(Fweight, i, j, k);
                direct::elem(counter, ires) += 1.0;
            }
        }

        // Calculate 1/1000th of radial averaged weight
        for (long int i = 0; i < Xsize(radavg_weight); i++) {
            if (
                direct::elem(counter,       i) > 0.0 ||
                direct::elem(radavg_weight, i) > 0.0
            ) {
                direct::elem(radavg_weight, i) /= 1000.0 * direct::elem(counter, i);
            } else {
                std::cerr << " counter= " << counter << std::endl;
                std::cerr << " radavg_weight= " << radavg_weight << std::endl;
                REPORT_ERROR("BUG: zeros in counter or radavg_weight!");
            }
        }

        bool have_warned = false;
        // perform std::max on all weight elements, and do division of data/weight
        FOR_ALL_ELEMENTS_IN_FFTW_TRANSFORM(Fweight) {
            const int ires = floor(hypot((double) ip, jp, kp) / padding_factor);
            const RFLOAT weight =  std::max(
                direct::elem(Fweight, i, j, k),
                direct::elem(radavg_weight, ires < r_max ? ires : r_max - 1)
            );
            if (weight == 0 && !have_warned) {
                std::cerr << " WARNING: ignoring divide by zero in skip_gridding: "
                    << "ires = " << ires
                    << " ip = " << ip << " jp = " << jp << " kp = " << kp << std::endl
                    << " max_r2 = " << max_r2 << " r_max = " << r_max << " "
                    << "padding_factor = " << padding_factor
                    << " round(sqrt(max_r2)) = " << round(sqrt(max_r2)) << " "
                    << "round(r_max * padding_factor) = " << round(r_max * padding_factor)
                    << std::endl;
                have_warned = true;
            } else {
                direct::elem(Fconv, i, j, k) /= weight;
            }
        }
    } else {

        {
        ifdefTIMING(TicToc tt (ReconTimer, ReconS[4]);)
        // Divide both data and Fweight by normalisation factor to prevent FFT's with very large values....
        #ifdef DEBUG_RECONSTRUCT
        std::cerr << " normalise= " << normalise << std::endl;
        #endif
        Fweight /= normalise;
        data /= (Complex) normalise;
        }

        // Fnewweight can become too large for a float: always keep this one in double-precision
        MultidimArray<double> Fnewweight = [&] () {
        ifdefTIMING(TicToc tt (ReconTimer, ReconS[5]);)
        // Initialise Fnewweight with 1's and 0's. (also see comments below)
        /// XXX: But this is changing weight!?
        FOR_ALL_ELEMENTS_IN_ARRAY3D(weight, i, j, k) {
            weight.elem(i, j, k) = hypot2(i, j, k) < max_r2;
        }
        return decentered(weight, max_r2, pad_size / 2 + 1, pad_size, pad_size);
        }();

        // Iterative algorithm as in Eq. 14 in Pipe & Menon (1999)
        // or Eq. 4 in Matej (2001)
        for (int iter = 0; iter < max_iter_preweight; iter++) {
            // std::cout << "    iteration " << iter + 1 << "/" << max_iter_preweight << "\n";
            {
            ifdefTIMING(TicToc tt (ReconTimer, ReconS[6]);)
            // Set Fnewweight * Fweight in the transformer
            // In Matej et al (2001), weights w_P^i are convoluted with the kernel,
            // and the initial w_P^0 are 1 at each sampling point
            // Here the initial weights are also 1 (see initialisation Fnewweight above),
            // but each "sampling point" counts "Fweight" times!
            // That is why Fnewweight is multiplied by Fweight prior to the convolution

            Fconv = Fnewweight * Fweight;

            // convolute through Fourier-transform (as both grids are rectangular)
            // Note that convoluteRealSpace acts on the complex array inside the transformer
            convoluteBlobRealSpace(transformer, false);

            RFLOAT w, corr_min = LARGE_NUMBER, corr_max = -LARGE_NUMBER, corr_avg = 0.0, corr_nn = 0.0;

            FOR_ALL_ELEMENTS_IN_FFTW_TRANSFORM(Fconv) {
                if (hypot2(ip, jp, kp) < max_r2) {

                    // Make sure no division by zero can occur....
                    w = std::max(1e-6, abs(direct::elem(Fconv, i, j, k)));
                    // Monitor min, max and avg conv_weight
                    corr_min = std::min(corr_min, w);
                    corr_max = std::max(corr_max, w);
                    corr_avg += w;
                    corr_nn += 1.0;
                    // Apply division of Eq. [14] in Pipe & Menon (1999)
                    direct::elem(Fnewweight, i, j, k) /= w;
                }
            }

            }

            #ifdef DEBUG_RECONSTRUCT
            std::cerr << " PREWEIGHTING ITERATION: " << iter + 1 << " OF " << max_iter_preweight << std::endl;
            // report of maximum and minimum values of current conv_weight
            std::cerr << " corr_avg= " << corr_avg / corr_nn << std::endl;
            std::cerr << " corr_min= " << corr_min << std::endl;
            std::cerr << " corr_max= " << corr_max << std::endl;
            #endif
        }

        {
        ifdefTIMING(TicToc tt (ReconTimer, ReconS[7]);)

        #ifdef DEBUG_RECONSTRUCT
        Image<double> img (Fnewweight);
        img.write("reconstruct_gridding_weight.spi");
        for (long int n = 0; n < Fconv.size(); n++) {
            img()[n] = abs(Fconv[n]);
        }
        img.write("reconstruct_gridding_correction_term.spi");
        #endif

        // Clear memory
        Fweight.clear();

        // Note that Fnewweight now holds the approximation of the inverse of the weights on a regular grid

        // Now do the actual reconstruction with the data array
        // Apply the iteratively determined weight
        Fconv = Projector::decentered(data, max_r2, pad_size / 2 + 1, pad_size, pad_size);
        for (long int n = 0; n < Fconv.size(); n++) {
            #ifdef RELION_SINGLE_PRECISION
            // Prevent numerical instabilities in single-precision reconstruction with very unevenly sampled orientations
            if (Fnewweight[n] > 1e20) { Fnewweight[n] = 1e20; }
            #endif
            Fconv[n] *= Fnewweight[n];
        }

        // Clear memory
        Fnewweight.clear();

        }
    }

    // Gridding theory says one now has to interpolate the fine grid onto the coarse one using a blob kernel
    // and then do the inverse transform and divide by the FT of the blob (i.e. do the gridding correction)
    // In practice, this gives all types of artefacts (perhaps I never found the right implementation?!)
    // Therefore, window the Fourier transform and then do the inverse transform
    // #define RECONSTRUCT_CONVOLUTE_BLOB
    #ifdef RECONSTRUCT_CONVOLUTE_BLOB

    // Apply the same blob-convolution as above to the data array
    // Mask real-space map beyond its original size to prevent aliasing in the downsampling step below
    {
    ifdefTIMING(TicToc tt (ReconTimer, ReconS[8]);)
    convoluteBlobRealSpace(transformer, true);
    }

    {
    ifdefTIMING(TicToc tt (ReconTimer, ReconS[9]);)
    // Now just pick every 3rd pixel in Fourier-space (i.e. down-sample)
    // and do a final inverse FT
    if (ref_dim == 2) {
        vol_out.resize(ori_size, ori_size);
    } else {
        vol_out.resize(ori_size, ori_size, ori_size);
    }
    }

    auto Ftmp = [&] () -> MultidimArray<Complex> {
    ifdefTIMING(TicToc tt (ReconTimer, ReconS[10]);)
    FourierTransformer transformer;  // Different transformer from outer scope, to preserve Fconv.
    transformer.setReal(vol_out);
    return transformer.getFourier();  // std::move?
    }();

    {
    ifdefTIMING(TicToc tt (ReconTimer, ReconS[11]);)
    FOR_ALL_ELEMENTS_IN_FFTW_TRANSFORM(Ftmp) {
        direct::elem(Ftmp, i, j, k) = hypot2(ip, jp, kp) < r_max * r_max ?
            FFTW::elem(Fconv, ip * padding_factor, jp * padding_factor, kp * padding_factor) :
            0.0;
        }
    }

    // Any particular reason why 13 comes before 12?

    {
    ifdefTIMING(TicToc tt (ReconTimer, ReconS[13]);)
    CenterFFTbySign(Ftmp);
    }

    {
    ifdefTIMING(TicToc tt (ReconTimer, ReconS[12]);)
    // inverse FFT leaves result in vol_out
    transformer2.inverseFourierTransform();
    }

    {
    ifdefTIMING(TicToc tt (ReconTimer, ReconS[14]);)
    // Un-normalize FFTW (because original FFTs were done with the size of 2D FFTs)
    if (ref_dim == 3) { vol_out /= ori_size; }
    }

    {
    ifdefTIMING(TicToc tt (ReconTimer, ReconS[15]);)
    // Mask out corners to prevent aliasing artefacts
    softMaskOutsideMap(vol_out);
    }

    {
    ifdefTIMING(TicToc tt (ReconTimer, ReconS[16]);)
    // Gridding correction for the blob
    RFLOAT normftblob = tab_ftblob(0.0);
    FOR_ALL_ELEMENTS_IN_ARRAY3D(vol_out, i, j, k) {

        RFLOAT r = hypot(i, j, k);
        RFLOAT rval = r / (ori_size * padding_factor);
        vol_out.elem(i, j, k) /= tab_ftblob(rval) / normftblob;
        // if (k == 0 && i == 0)
        // 	std::cerr << " j= " << j << " rval= " << rval << " tab_ftblob(rval) / normftblob= " << tab_ftblob(rval) / normftblob << std::endl;
    }
    }

    #else

    // rather than doing the blob-convolution to downsample the data array, do a windowing operation:
    // This is the same as convolution with a SINC. It seems to give better maps.
    // Then just make the blob look as much as a SINC as possible....
    // The "standard" r1.9, m2 and a15 blob looks quite like a sinc until the first zero (perhaps that's why it is standard?)
    //for (RFLOAT r = 0.1; r < 10.0; r += 0.01) {
    //	RFLOAT sinc_theta = sinc(PI * r / padding_factor );
    //	std::cout << " r= " << r << " sinc= " << sinc_theta << " blob= " << blob_val(r, blob) << std::endl;
    //}

    // Now do inverse FFT and window to original size in real-space
    // Pass the transformer to prevent making and clearing a new one before clearing the one declared above....
    // The latter may give memory problems as detected by electric fence....
    {
    ifdefTIMING(TicToc tt (ReconTimer, ReconS[17]);)
    windowToOridimRealSpace(transformer, vol_out, printTimes);
    }

    #endif

    #ifdef DEBUG_RECONSTRUCT
    Image<double>(vol_out).write("reconstruct_before_gridding_correction.spi");
    #endif

    // Correct for the linear/nearest-neighbour interpolation that led to the data array
    {
    ifdefTIMING(TicToc tt (ReconTimer, ReconS[18]);)
    griddingCorrect(vol_out);
    }

    {
    ifdefTIMING(TicToc tt (ReconTimer, ReconS[23]);)
    // Completely empty the transformer object
    transformer.cleanup();
    // Now can use extra mem to move data into smaller array space
    vol_out.shrinkToFit();
    }

    #ifdef TIMING
    if (printTimes) ReconTimer.printTimes(true);
    #endif

    #ifdef DEBUG_RECONSTRUCT
    std::cerr << "done with reconstruct" << std::endl;
    #endif

    if (!weight_out) return vol_out;

    weight_out->data    = MultidimArray<RFLOAT>(ori_size / 2 + 1, ori_size, ori_size);
    Image<RFLOAT> count = Image<RFLOAT>::zeros (ori_size / 2 + 1, ori_size, ori_size);

    // downsample while considering padding:

    for (long int z = 0; z < Fweight.zdim; z++)
    for (long int y = 0; y < Fweight.ydim; y++)
    for (long int x = 0; x < Fweight.xdim; x++) {
        int xl = x;
        int yl = y < Fweight.ydim / 2 ? y : y - Fweight.ydim;
        int zl = z < Fweight.zdim / 2 ? z : z - Fweight.zdim;

        if (
            xl == Fweight.xdim - 1 ||
            yl == Fweight.ydim / 2 || yl == -Fweight.ydim / 2 - 1 ||
            zl == Fweight.zdim / 2 || zl == -Fweight.zdim / 2 - 1
        ) continue;

        int xx =        round(xl / padding_factor);
        int yy = ((int) round(yl / padding_factor) + ori_size) % ori_size;
        int zz = ((int) round(zl / padding_factor) + ori_size) % ori_size;

        if (
               xx >= 0 && xx < ori_size / 2 + 1
            && yy >= 0 && yy < ori_size
            && zz >= 0 && zz < ori_size
        ) {
            direct::elem(weight_out->data, xx, yy, zz) += direct::elem(Fweight, x, y, z);
            direct::elem(count.data,       xx, xx, zz) += 1.0;
        }
    }

    const double pad3 = padding_factor * padding_factor * padding_factor;

    for (long int z = 0; z < ori_size; z++)
    for (long int y = 0; y < ori_size; y++)
    for (long int x = 0; x < ori_size / 2 + 1; x++) {
        const RFLOAT c = direct::elem(count.data, x, y, z);

        if (c > 0.0) {
            direct::elem(weight_out->data, x, y, z) *= pad3 / c;
        }
    }
    return vol_out;
}

void BackProjector::symmetrise(
    int nr_helical_asu, RFLOAT helical_twist, RFLOAT helical_rise, int threads
) {
    // First make sure the input arrays are obeying Hermitian symmetry,
    // which is assumed in the rotation operators of both helical and point group symmetry
    enforceHermitianSymmetry();

    // Then apply helical and point group symmetry (order irrelevant?)
    applyHelicalSymmetry(nr_helical_asu, helical_twist, helical_rise);

    applyPointGroupSymmetry(threads);
}

void BackProjector::enforceHermitianSymmetry() {
    // Make sure all points are only included once.
    for (int iz = Zinit(data); iz <= Zlast(data); iz++)
    for (int iy = iz >= 0;     iy <= Ylast(data); iy++) {
        Complex fsum = data.elem(0, +iy, +iz) + conj(data.elem(0, -iy, -iz));
        data.elem(0, +iy, +iz) =      fsum;
        data.elem(0, -iy, -iz) = conj(fsum);
        RFLOAT sum = weight.elem(0, +iy, +iz) + weight.elem(0, -iy, -iz);
        weight.elem(0, +iy, +iz) = sum;
        weight.elem(0, -iy, -iz) = sum;
    }
}

void BackProjector::applyHelicalSymmetry(
    int nr_helical_asu, RFLOAT helical_twist, RFLOAT helical_rise
) {
    if (nr_helical_asu < 2 || ref_dim != 3) return;

    int rmax2 = round(r_max * padding_factor) * round(r_max * padding_factor);

    MultidimArray<RFLOAT>  sum_weight = weight;
    MultidimArray<Complex> sum_data   = data;
    // First symmetry operator (not stored in SL) is the identity matrix
    int h_min = -nr_helical_asu / 2;
    int h_max = -h_min + nr_helical_asu % 2;
    for (int hh = h_min; hh < h_max; hh++) {
        if (hh == 0) continue;  // h == 0 is done before the loop (where sum_data = data)

        Matrix<RFLOAT> R = rotation3DMatrix(-hh * helical_twist, 'Z');
        setSmallValuesToZero(R.begin(), R.end());
        // TODO: invert rotation matrix?

        // Loop over all points in the output (i.e. rotated, or summed) array
        FOR_ALL_ELEMENTS_IN_ARRAY3D(sum_weight, i, j, k) {
            RFLOAT x = i;  // Xinit(sum_weight) is zero!
            RFLOAT y = j;
            RFLOAT z = k;
            RFLOAT r2 = hypot2(x, y, z);
            if (r2 <= rmax2) {
                // coords_output(x, y) = A * coords_input (xp, yp)
                // {xp, yp zp} = R({x, y, z})
                RFLOAT xp = x * R(0, 0) + y * R(0, 1) + z * R(0, 2);
                RFLOAT yp = x * R(1, 0) + y * R(1, 1) + z * R(1, 2);
                RFLOAT zp = x * R(2, 0) + y * R(2, 1) + z * R(2, 2);

                const bool is_neg_x = xp < 0;
                // Only asymmetric half is stored
                if (is_neg_x) {
                    // Get complex conjugated hermitian symmetry pair
                    xp = -xp;
                    yp = -yp;
                    zp = -zp;
                }

                // Trilinear interpolation (with physical coords)
                // Subtract X/Y/Zinit for faster data access
                // In that way use direct::elem, rather than elem
                int x0 = floor(xp);
                RFLOAT fx = xp - x0;
                x0 -= Xinit(data);  // no-op since Xinit(data) is 0
                int x1 = x0 + 1;

                int y0 = floor(yp);
                RFLOAT fy = yp - y0;
                y0 -= Yinit(data);
                int y1 = y0 + 1;

                int z0 = floor(zp);
                RFLOAT fz = zp - z0;
                z0 -= Zinit(data);
                int z1 = z0 + 1;

                #ifdef CHECK_SIZE
                if (
                    x0 < 0 || y0 < 0 || z0 < 0 ||
                    x1 < 0 || y1 < 0 || z1 < 0 ||
                    x0 >= Xsize(data) || y0 >= Ysize(data) || z0 >= Zsize(data) ||
                    x1 >= Xsize(data) || y1 >= Ysize(data) || z1 >= Zsize(data)
                ) {
                    std::cerr << " x0= " << x0 << " y0= " << y0 << " z0= " << z0 << std::endl;
                    std::cerr << " x1= " << x1 << " y1= " << y1 << " z1= " << z1 << std::endl;
                    data.printShape();
                    REPORT_ERROR("BackProjector::applyPointGroupSymmetry: checksize!!!");
                }
                #endif
                // First interpolate (complex) data
                Complex d000 = direct::elem(data, x0, y0, z0);
                Complex d001 = direct::elem(data, x1, y0, z0);
                Complex d010 = direct::elem(data, x0, y1, z0);
                Complex d011 = direct::elem(data, x1, y1, z0);
                Complex d100 = direct::elem(data, x0, y0, z1);
                Complex d101 = direct::elem(data, x1, y0, z1);
                Complex d110 = direct::elem(data, x0, y1, z1);
                Complex d111 = direct::elem(data, x1, y1, z1);

                Complex dx00 = LIN_INTERP(fx, d000, d001);
                Complex dx01 = LIN_INTERP(fx, d100, d101);
                Complex dx10 = LIN_INTERP(fx, d010, d011);
                Complex dx11 = LIN_INTERP(fx, d110, d111);
                Complex dxy0 = LIN_INTERP(fy, dx00, dx10);
                Complex dxy1 = LIN_INTERP(fy, dx01, dx11);

                // Take complex conjugated for half with negative x
                Complex ddd = LIN_INTERP(fz, dxy0, dxy1);

                if (is_neg_x) { ddd = conj(ddd); }

                // Also apply a phase shift for helical translation along Z
                if (abs(helical_rise) > 0.0) {
                    RFLOAT zshift = -hh * helical_rise / (ori_size * padding_factor);
                    RFLOAT phase = 2 * PI * z * zshift;
                    const Complex X = Complex::unit(phase);
                    const Complex Z = Complex(X.real * ddd.real, X.imag * ddd.imag);
                    ddd = Complex(
                        Z.real - Z.imag,
                        (X.real + X.imag) * (ddd.real + ddd.imag) - Z.real - Z.imag
                    );
                }
                // Accumulated sum of the data term
                sum_data.elem(i, j, k) += ddd;

                // Then interpolate (real) weight
                RFLOAT dd000 = direct::elem(weight, x0, y0, z0);
                RFLOAT dd001 = direct::elem(weight, x1, y0, z0);
                RFLOAT dd010 = direct::elem(weight, x0, y1, z0);
                RFLOAT dd011 = direct::elem(weight, x1, y1, z0);
                RFLOAT dd100 = direct::elem(weight, x0, y0, z1);
                RFLOAT dd101 = direct::elem(weight, x1, y0, z1);
                RFLOAT dd110 = direct::elem(weight, x0, y1, z1);
                RFLOAT dd111 = direct::elem(weight, x1, y1, z1);

                RFLOAT ddx00 = LIN_INTERP(fx, dd000, dd001);
                RFLOAT ddx01 = LIN_INTERP(fx, dd100, dd101);
                RFLOAT ddx10 = LIN_INTERP(fx, dd010, dd011);
                RFLOAT ddx11 = LIN_INTERP(fx, dd110, dd111);
                RFLOAT ddxy0 = LIN_INTERP(fy, ddx00, ddx10);
                RFLOAT ddxy1 = LIN_INTERP(fy, ddx01, ddx11);

                sum_weight.elem(i, j, k) += LIN_INTERP(fz, ddxy0, ddxy1);

            }
        }
    }

    data  = sum_data;
    weight = sum_weight;
}

void BackProjector::applyPointGroupSymmetry(int threads) {

    // #define DEBUG_SYMM
    #ifdef DEBUG_SYMM
    std::cerr << " SL.SymsNo()= " << SL.SymsNo() << std::endl;
    std::cerr << " SL.true_symNo= " << SL.true_symNo << std::endl;
    #endif

    int rmax2 = round(r_max * padding_factor) * round(r_max * padding_factor);
    if (SL.SymsNo() > 0 && ref_dim == 3) {
        Matrix<RFLOAT> L(4, 4), R(4, 4);  // A matrix from the list
        // First symmetry operator (not stored in SL) is the identity matrix
        MultidimArray<RFLOAT>  sum_weight = weight;
        MultidimArray<Complex> sum_data   = data;

        // Loop over all other symmetry operators
        for (int isym = 0; isym < SL.SymsNo(); isym++) {
            SL.get_matrices(isym, L, R);
            #ifdef DEBUG_SYMM
            std::cerr << " isym= " << isym << " R= " << R << std::endl;
            #endif
            // Loop over all points in the output (i.e. rotated, or summed) array

            #pragma omp parallel for num_threads(threads)
            for (long int k = Zinit(sum_weight); k <= Zlast(sum_weight); k++)
            for (long int i = Yinit(sum_weight); i <= Ylast(sum_weight); i++)
            for (long int j = Xinit(sum_weight); j <= Xlast(sum_weight); j++) {
                RFLOAT x = j;  // Xinit(sum_weight) is zero!
                RFLOAT y = i;
                RFLOAT z = k;
                RFLOAT r2 = hypot2(x, y, z);

                if (r2 <= rmax2) {
                    // coords_output(x,y) = A * coords_input (xp,yp)
                    RFLOAT xp = x * R(0, 0) + y * R(0, 1) + z * R(0, 2);
                    RFLOAT yp = x * R(1, 0) + y * R(1, 1) + z * R(1, 2);
                    RFLOAT zp = x * R(2, 0) + y * R(2, 1) + z * R(2, 2);

                    const bool is_neg_x = xp < 0;

                    // Only asymmetric half is stored
                    if (is_neg_x) {
                        // Get complex conjugated hermitian symmetry pair
                        xp = -xp;
                        yp = -yp;
                        zp = -zp;
                    }

                    // Trilinear interpolation (with physical coords)
                    // Subtract Yinit and Zinit to accelerate access to data (Xinit=0)
                    // In that way use direct::elem, rather than elem
                    int x0 = floor(xp);
                    RFLOAT fx = xp - x0;
                    // x0 -= Xinit(data);
                    int x1 = x0 + 1;

                    int y0 = floor(yp);
                    RFLOAT fy = yp - y0;
                    y0 -= Yinit(data);
                    int y1 = y0 + 1;

                    int z0 = floor(zp);
                    RFLOAT fz = zp - z0;
                    z0 -= Zinit(data);
                    int z1 = z0 + 1;

                    #ifdef CHECK_SIZE
                    if (
                        x0 < 0 || y0 < 0 || z0 < 0 ||
                        x1 < 0 || y1 < 0 || z1 < 0 ||
                        x0 >= Xsize(data) || y0  >= Ysize(data) || z0 >= Zsize(data) ||
                        x1 >= Xsize(data) || y1  >= Ysize(data) || z1 >= Zsize(data)
                    ) {
                        std::cerr << " x0= " << x0 << " y0= " << y0 << " z0= " << z0 << std::endl;
                        std::cerr << " x1= " << x1 << " y1= " << y1 << " z1= " << z1 << std::endl;
                        data.printShape();
                        REPORT_ERROR("BackProjector::applyPointGroupSymmetry: checksize!!!");
                    }
                    #endif
                    // First interpolate (complex) data
                    Complex d000 = direct::elem(data, x0, y0, z0);
                    Complex d001 = direct::elem(data, x1, y0, z0);
                    Complex d010 = direct::elem(data, x0, y1, z0);
                    Complex d011 = direct::elem(data, x1, y1, z0);
                    Complex d100 = direct::elem(data, x0, y0, z1);
                    Complex d101 = direct::elem(data, x1, y0, z1);
                    Complex d110 = direct::elem(data, x0, y1, z1);
                    Complex d111 = direct::elem(data, x1, y1, z1);

                    Complex dx00 = LIN_INTERP(fx, d000, d001);
                    Complex dx01 = LIN_INTERP(fx, d100, d101);
                    Complex dx10 = LIN_INTERP(fx, d010, d011);
                    Complex dx11 = LIN_INTERP(fx, d110, d111);

                    Complex dxy0 = LIN_INTERP(fy, dx00, dx10);
                    Complex dxy1 = LIN_INTERP(fy, dx01, dx11);

                    // Take complex conjugated for half with negative x
                    sum_data.elem(i, j, k) += is_neg_x ?
                        conj(LIN_INTERP(fz, dxy0, dxy1)) : LIN_INTERP(fz, dxy0, dxy1);

                    // Then interpolate (real) weight
                    RFLOAT dd000 = direct::elem(weight, x0, y0, z0);
                    RFLOAT dd001 = direct::elem(weight, x1, y0, z0);
                    RFLOAT dd010 = direct::elem(weight, x0, y1, z0);
                    RFLOAT dd011 = direct::elem(weight, x1, y1, z0);
                    RFLOAT dd100 = direct::elem(weight, x0, y0, z1);
                    RFLOAT dd101 = direct::elem(weight, x1, y0, z1);
                    RFLOAT dd110 = direct::elem(weight, x0, y1, z1);
                    RFLOAT dd111 = direct::elem(weight, x1, y1, z1);

                    RFLOAT ddx00 = LIN_INTERP(fx, dd000, dd001);
                    RFLOAT ddx01 = LIN_INTERP(fx, dd100, dd101);
                    RFLOAT ddx10 = LIN_INTERP(fx, dd010, dd011);
                    RFLOAT ddx11 = LIN_INTERP(fx, dd110, dd111);

                    RFLOAT ddxy0 = LIN_INTERP(fy, ddx00, ddx10);
                    RFLOAT ddxy1 = LIN_INTERP(fy, ddx01, ddx11);

                    sum_weight.elem(i, j, k) += LIN_INTERP(fz, ddxy0, ddxy1);

                }
            }
        }

        data   = sum_data;
        weight = sum_weight;
        // Average
        // The division should only be done if we would search all (C1) directions, not if we restrict the angular search!
        /*
        for (long int n = 0; n < data.size(); n++) {
            data[n] = sum_data[n] / (RFLOAT)(SL.SymsNo() + 1);
            weight[n] = sum_weight[n] / (RFLOAT)(SL.SymsNo() + 1);
        }
        */
    }
}

void BackProjector::convoluteBlobRealSpace(FourierTransformer &transformer, bool do_mask) {

    // Set up right dimension of real-space array
    // TODO: resize this according to r_max!!!
    MultidimArray<RFLOAT> Mconv (pad_size, pad_size, ref_dim == 2 ? 1 : pad_size);

    // inverse FFT
    transformer.setReal(Mconv);
    transformer.inverseFourierTransform();

    // Blob normalisation in Fourier space
    RFLOAT normftblob = tab_ftblob(0.0);

    // TMP DEBUGGING
    //struct blobtype blob;
    //blob.order = 0;
    //blob.radius = 1.9 * padding_factor;
    //blob.alpha = 15;

    // Multiply with FT of the blob kernel
    const int padhdim = pad_size / 2;
    for (long int k = 0; k < Zsize(Mconv); k++)
    for (long int j = 0; j < Ysize(Mconv); j++)
    for (long int i = 0; i < Xsize(Mconv); i++) {
        int kp = k < padhdim ? k : k - pad_size;
        int jp = j < padhdim ? j : j - pad_size;
        int ip = i < padhdim ? i : i - pad_size;
        RFLOAT rval = hypot((double) ip, jp, kp) / (ori_size * padding_factor);
        //if (kp==0 && ip==0 && jp > 0)
        //	std::cerr << " jp= " << jp << " rval= " << rval << " tab_ftblob(rval) / normftblob= " << tab_ftblob(rval) / normftblob << " ori_size/2= " << ori_size/2 << std::endl;
        // In the final reconstruction: mask the real-space map beyond its original size to prevent aliasing ghosts
        // Note that rval goes until 1/2 in the oversampled map
        if (do_mask && 2.0 * padding_factor * rval > 1.0) {
            direct::elem(Mconv, i, j, k) = 0.0;
        } else {
            direct::elem(Mconv, i, j, k) *= tab_ftblob(rval) / normftblob;
        }
    }

    // forward FFT to go back to Fourier-space
    transformer.FourierTransform();
}

void BackProjector::windowToOridimRealSpace(
    FourierTransformer &transformer, MultidimArray<RFLOAT> &Mout, bool printTimes
) {

    #ifdef TIMING
    /// FIXME: These are out of order
    Timer OriDimTimer;
    int OrD1_getFourier = OriDimTimer.setNew(" OrD1_getFourier ");
    int OrD2_windowFFT  = OriDimTimer.setNew(" OrD2_windowFFT ");
    int OrD3_reshape    = OriDimTimer.setNew(" OrD3_reshape ");
    int OrD4_setReal    = OriDimTimer.setNew(" OrD4_setReal ");
    int OrD5_invFFT     = OriDimTimer.setNew(" OrD5_invFFT ");
    int OrD6_centerFFT  = OriDimTimer.setNew(" OrD6_centerFFT ");
    int OrD7_window     = OriDimTimer.setNew(" OrD7_window ");
    int OrD8_norm       = OriDimTimer.setNew(" OrD8_norm ");
    int OrD9_softMask   = OriDimTimer.setNew(" OrD9_softMask ");
    #endif

    auto &Fin = [&] () -> MultidimArray<Complex> & {
        ifdefTIMING(TicToc tt (OriDimTimer, OrD1_getFourier);)
        // Why are we timing this? It's just a member access.
        return transformer.getFourier();
    }();

    int padoridim;  // Size of padded real-space volume
    {
    ifdefTIMING(TicToc tt (OriDimTimer, OrD2_windowFFT);)
    padoridim = round(padding_factor * ori_size);
    // Enforce divisibility by 2
    padoridim += padoridim % 2;

    // #define DEBUG_WINDOWORIDIMREALSPACE
    #ifdef DEBUG_WINDOWORIDIMREALSPACE
    MultidimArray<RFLOAT> tmp (Xsize(Fin), Ysize(Fin), Zsize(Fin));
    for (long int n = 0; n < Fin.size(); n++) {
        tmp()[n] = abs(Fin[n]);
    }
    Image<RFLOAT>(tmp).write("windoworidim_Fin.spi");
    #endif

    // Resize incoming complex array to the correct size
    windowFourierTransform(Fin, padoridim);
    }

    const RFLOAT normfft = [&] () {
        ifdefTIMING(TicToc tt (OriDimTimer, OrD3_reshape);)
        if (ref_dim == 2) {
            Mout.reshape(padoridim, padoridim);
            Mout.setXmippOrigin();
            return padding_factor * padding_factor
                * (data_dim == 2 ? 1 : ori_size);
        } else {
            Mout.reshape(padoridim, padoridim, padoridim);
            Mout.setXmippOrigin();
            return padding_factor * padding_factor * padding_factor
                * (data_dim == 3 ? 1 : ori_size);
        }
    }();

    #ifdef DEBUG_WINDOWORIDIMREALSPACE
    tmp.reshape(Xsize(Fin), Ysize(Fin), Zsize(Fin));
    for (long int n = 0; n < Fin.size(); n++) {
        tmp[n] = abs(Fin[n]);
    }
    Image<RFLOAT>(tmp).write("windoworidim_Fresized.spi");
    #endif

    // Shift the map back to its origin
    {
    ifdefTIMING(TicToc tt (OriDimTimer, OrD6_centerFFT);)
    CenterFFTbySign(Fin);
    }

    {
    ifdefTIMING(TicToc tt (OriDimTimer, OrD4_setReal);)
    transformer.setReal(Mout);
    }

    // Do the inverse FFT
    {
    ifdefTIMING(TicToc tt (OriDimTimer, OrD5_invFFT);)
    #ifdef TIMING
    if (printTimes)
        std::cout << std::endl << "FFTrealDims = (" << transformer.fReal->xdim << " , " << transformer.fReal->ydim << " , " << transformer.fReal->zdim << " ) " << std::endl;
    #endif
    transformer.inverseFourierTransform();
    }

    // Mout = transformer.inverseFourierTransform(Fin);

    Fin.clear();
    transformer.fReal = nullptr;  // Make sure to re-calculate fftw plan
    Mout.setXmippOrigin();

    #ifdef DEBUG_WINDOWORIDIMREALSPACE
    Image<RFLOAT>(Mout).write("windoworidim_Munwindowed.spi");
    #endif

    // Window in real-space
    {
    ifdefTIMING(TicToc tt (OriDimTimer, OrD7_window);)
    const long int init = Xmipp::init(ori_size), last = Xmipp::last(ori_size);
    Mout = (ref_dim == 2 ?
        Mout.windowed(init, init, last, last) :
        Mout.windowed(init, init, init, last, last, last)).setXmippOrigin();
    }

    // Normalisation factor of FFTW
    // The Fourier Transforms are all "normalised" for 2D transforms of size = ori_size * ori_size
    {
    ifdefTIMING(TicToc tt (OriDimTimer, OrD8_norm);)
    Mout /= normfft;
    }

    #ifdef DEBUG_WINDOWORIDIMREALSPACE
    Image<RFLOAT>(Mout).write("windoworidim_Mwindowed.spi");
    #endif

    // Mask out corners to prevent aliasing artefacts
    {
    ifdefTIMING(TicToc tt (OriDimTimer, OrD9_softMask);)
    softMaskOutsideMap(Mout);
    }

    #ifdef DEBUG_WINDOWORIDIMREALSPACE
    Image<RFLOAT>(Mout).write("windoworidim_Mwindowed_masked.spi");
    // FourierTransformer ttf;
    Fin = FourierTransformer{}.FourierTransform(Mout);
    imtt().resize(Zsize(Fin), Ysize(Fin), Xsize(Fin));
    for (long int n = 0; n < Fin.size(); n++) {
        imtt()[n] = abs(Fin[n]);
    }
    imtt.write("windoworidim_Fnew.spi");
    #endif

    #ifdef TIMING
    if (printTimes) OriDimTimer.printTimes(true);
    #endif
}
