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

#ifdef TIMING
#define RCTIC(timer, label) (timer.tic(label))
#define RCTOC(timer, label) (timer.toc(label))
#else
#define RCTIC(timer, label)
#define RCTOC(timer, label)
#endif

#define RCTICTOC(timer, label, block) RCTIC(timer, label); { block; } RCTOC(timer, label);


void BackProjector::initialiseDataAndWeight(int current_size) {
    initialiseData(current_size);
    weight.resize(data);
}

void BackProjector::initZeros(int current_size) {
    initialiseDataAndWeight(current_size);
    data.initZeros();
    weight.initZeros();
}

void BackProjector::backproject2Dto3D(
    const MultidimArray<Complex > &f2d,
    const Matrix2D<RFLOAT> &A,
    const MultidimArray<RFLOAT> *Mweight,
    RFLOAT r_ewald_sphere, bool is_positive_curvature,
    Matrix2D<RFLOAT> *magMatrix
) {
    RFLOAT m00, m10, m01, m11;

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

    // Use the inverse matrix
    Matrix2D<RFLOAT> Ainv;
    Ainv = A.inv();

    // Go from the 2D slice coordinates to the 3D coordinates
    Ainv *= (RFLOAT) padding_factor;  // take scaling into account directly

    // max_r2 and min_r2_nn are defined in 3D-space
    const int max_r2    = round(r_max    * padding_factor) * round(r_max    * padding_factor);
    const int min_r2_nn = round(r_min_nn * padding_factor) * round(r_min_nn * padding_factor);

    // precalculated coefficients for ellipse determination (see further down)

    // first, make sure A contains 2D distortion (lowercase 2D, uppercase 3D):
    const RFLOAT Am_Xx = Ainv(0, 0) * m00 + Ainv(0, 1) * m10;
    const RFLOAT Am_Xy = Ainv(0, 0) * m01 + Ainv(0, 1) * m11;
    const RFLOAT Am_Yx = Ainv(1, 0) * m00 + Ainv(1, 1) * m10;
    const RFLOAT Am_Yy = Ainv(1, 0) * m01 + Ainv(1, 1) * m11;
    const RFLOAT Am_Zx = Ainv(2, 0) * m00 + Ainv(2, 1) * m10;
    const RFLOAT Am_Zy = Ainv(2, 0) * m01 + Ainv(2, 1) * m11;

    // next, precompute (Am)^t Am into AtA:
    const RFLOAT AtA_xx = Am_Xx * Am_Xx + Am_Yx * Am_Yx + Am_Zx * Am_Zx;
    const RFLOAT AtA_xy = Am_Xx * Am_Xy + Am_Yx * Am_Yy + Am_Zx * Am_Zy;
    const RFLOAT AtA_yy = Am_Xy * Am_Xy + Am_Yy * Am_Yy + Am_Zy * Am_Zy;
    const RFLOAT AtA_xy2 = AtA_xy * AtA_xy;

    //#define DEBUG_BACKP
    #ifdef DEBUG_BACKP
    std::cerr << " XSIZE(f2d)= "<< XSIZE(f2d) << std::endl;
    std::cerr << " YSIZE(f2d)= "<< YSIZE(f2d) << std::endl;
    std::cerr << " XSIZE(data)= "<< XSIZE(data) << std::endl;
    std::cerr << " YSIZE(data)= "<< YSIZE(data) << std::endl;
    std::cerr << " STARTINGX(data)= "<< STARTINGX(data) << std::endl;
    std::cerr << " STARTINGY(data)= "<< STARTINGY(data) << std::endl;
    std::cerr << " STARTINGZ(data)= "<< STARTINGZ(data) << std::endl;
    std::cerr << " r_max= "<< r_max << std::endl;
    std::cerr << " Ainv= " << Ainv << std::endl;
    #endif

    // precalculate inverse of Ewald sphere diameter
    RFLOAT inv_diam_ewald = r_ewald_sphere > 0.0 ? 1.0 / 2.0 * r_ewald_sphere : 0.0;

    if (!is_positive_curvature) {
        inv_diam_ewald *= -1.0;
    }

    const int s  = YSIZE(f2d);
    const int sh = XSIZE(f2d);

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
            Complex my_val = DIRECT_A2D_ELEM(f2d, i, x);

            // Get the weight
            RFLOAT my_weight = Mweight ? DIRECT_A2D_ELEM(*Mweight, i, x) : 1.0;

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
            if (r2_3D > max_r2) {
                continue;
            }

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
                // Subtract STARTINGY and STARTINGZ to accelerate access to data (STARTINGX=0)
                // In that way use DIRECT_A3D_ELEM, rather than A3D_ELEM
                int x0 = floor(xp);
                RFLOAT fx = xp - x0;
                int x1 = x0 + 1;

                int y0 = floor(yp);
                RFLOAT fy = yp - y0;
                y0 -=  STARTINGY(data);
                int y1 = y0 + 1;

                int z0 = floor(zp);
                RFLOAT fz = zp - z0;
                z0 -= STARTINGZ(data);
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
                DIRECT_A3D_ELEM(data, z0, y0, x0) += dd000 * my_val;
                DIRECT_A3D_ELEM(data, z0, y0, x1) += dd001 * my_val;
                DIRECT_A3D_ELEM(data, z0, y1, x0) += dd010 * my_val;
                DIRECT_A3D_ELEM(data, z0, y1, x1) += dd011 * my_val;
                DIRECT_A3D_ELEM(data, z1, y0, x0) += dd100 * my_val;
                DIRECT_A3D_ELEM(data, z1, y0, x1) += dd101 * my_val;
                DIRECT_A3D_ELEM(data, z1, y1, x0) += dd110 * my_val;
                DIRECT_A3D_ELEM(data, z1, y1, x1) += dd111 * my_val;
                // Store corresponding weights
                DIRECT_A3D_ELEM(weight, z0, y0, x0) += dd000 * my_weight;
                DIRECT_A3D_ELEM(weight, z0, y0, x1) += dd001 * my_weight;
                DIRECT_A3D_ELEM(weight, z0, y1, x0) += dd010 * my_weight;
                DIRECT_A3D_ELEM(weight, z0, y1, x1) += dd011 * my_weight;
                DIRECT_A3D_ELEM(weight, z1, y0, x0) += dd100 * my_weight;
                DIRECT_A3D_ELEM(weight, z1, y0, x1) += dd101 * my_weight;
                DIRECT_A3D_ELEM(weight, z1, y1, x0) += dd110 * my_weight;
                DIRECT_A3D_ELEM(weight, z1, y1, x1) += dd111 * my_weight;

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

                const int xr = x0 - STARTINGX(data);
                const int yr = y0 - STARTINGY(data);
                const int zr = z0 - STARTINGZ(data);

                if (
                    xr < 0 || xr >= data.xdim ||
                    yr < 0 || yr >= data.ydim ||
                    zr < 0 || zr >= data.zdim
                ) continue;

                if (is_neg_x) {
                    DIRECT_A3D_ELEM(data,   zr, yr, xr) += conj(my_val);
                    DIRECT_A3D_ELEM(weight, zr, yr, xr) += my_weight;
                } else {
                    DIRECT_A3D_ELEM(data,   zr, yr, xr) += my_val;
                    DIRECT_A3D_ELEM(weight, zr, yr, xr) += my_weight;
                }
            } else {
                REPORT_ERROR("FourierInterpolator::backproject%%ERROR: unrecognized interpolator ");
            }
        }
    }
}

void BackProjector::backproject1Dto2D(
    const MultidimArray<Complex> &f1d,
    const Matrix2D<RFLOAT> &A,
    const MultidimArray<RFLOAT> *Mweight
) {
    Matrix2D<RFLOAT> Ainv = A.inv();
    Ainv *= (RFLOAT) padding_factor;  // take scaling into account directly

    const int r_max_src = XSIZE(f1d) - 1;

    const int r_max_ref   = r_max * padding_factor;
    const int r_max_ref_2 = r_max_ref * r_max_ref;

    // currently not used for some reason
    //const int r_min_NN_ref_2 = r_min_nn * r_min_nn * padding_factor * padding_factor;

    for (int x = 0; x <= r_max_src; x++) {
        RFLOAT my_weight = Mweight ? DIRECT_A1D_ELEM(*Mweight, x) : 1.0;
        if (my_weight <= 0.0) continue;

        Complex my_val = DIRECT_A1D_ELEM(f1d, x);

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
            // Subtract STARTINGY to accelerate access to data (STARTINGX=0)
            // In that way use DIRECT_A2D_ELEM, rather than A2D_ELEM
            const int x0 = floor(xp);
            const RFLOAT fx = xp - x0;
            const int x1 = x0 + 1;

            int y0 = floor(yp);
            const RFLOAT fy = yp - y0;
            y0 -=  STARTINGY(data);
            const int y1 = y0 + 1;

            const RFLOAT mfx = 1.0 - fx;
            const RFLOAT mfy = 1.0 - fy;

            const RFLOAT dd00 = mfy * mfx;
            const RFLOAT dd01 = mfy *  fx;
            const RFLOAT dd10 =  fy * mfx;
            const RFLOAT dd11 =  fy *  fx;

            if (is_neg_x) {
                my_val = conj(my_val);
            }

            // Store slice in 3D weighted sum
            DIRECT_A2D_ELEM(data, y0, x0) += dd00 * my_val;
            DIRECT_A2D_ELEM(data, y0, x1) += dd01 * my_val;
            DIRECT_A2D_ELEM(data, y1, x0) += dd10 * my_val;
            DIRECT_A2D_ELEM(data, y1, x1) += dd11 * my_val;

            // Store corresponding weights
            DIRECT_A2D_ELEM(weight, y0, x0) += dd00 * my_weight;
            DIRECT_A2D_ELEM(weight, y0, x1) += dd01 * my_weight;
            DIRECT_A2D_ELEM(weight, y1, x0) += dd10 * my_weight;
            DIRECT_A2D_ELEM(weight, y1, x1) += dd11 * my_weight;

        } else if (interpolator == NEAREST_NEIGHBOUR ) {
            const int x0 = round(xp);
            const int y0 = round(yp);

            if (x0 < 0) {
                A2D_ELEM(data,   -y0, -x0) += conj(my_val);
                A2D_ELEM(weight, -y0, -x0) += my_weight;
            } else {
                A2D_ELEM(data,   y0, x0) += my_val;
                A2D_ELEM(weight, y0, x0) += my_weight;
            }
        } else {
            REPORT_ERROR("FourierInterpolator::backproject1Dto2D%%ERROR: unrecognized interpolator ");
        }
    }
}

void BackProjector::backrotate2D(
    const MultidimArray<Complex > &f2d,
    const Matrix2D<RFLOAT> &A,
    const MultidimArray<RFLOAT> *Mweight,
    Matrix2D<RFLOAT>* magMatrix
) {
    Matrix2D<RFLOAT> Ainv = A.inv();
    Ainv *= (RFLOAT) padding_factor;  // take scaling into account directly

    RFLOAT m00, m10, m01, m11;

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

    const int r_max_ref   = r_max * padding_factor;
    const int r_max_ref_2 = r_max_ref * r_max_ref;

    int min_r2_nn = r_min_nn * r_min_nn * padding_factor * padding_factor;

    // precalculated coefficients for ellipse determination (see further down)

    // first, make sure A contains 2D distortion (lowercase 2D, uppercase 3D):
    const RFLOAT Am_Xx = Ainv(0, 0) * m00 + Ainv(0, 1) * m10;
    const RFLOAT Am_Xy = Ainv(0, 0) * m01 + Ainv(0, 1) * m11;
    const RFLOAT Am_Yx = Ainv(1, 0) * m00 + Ainv(1, 1) * m10;
    const RFLOAT Am_Yy = Ainv(1, 0) * m01 + Ainv(1, 1) * m11;

    // next, precompute (Am)^t Am into AtA:
    const RFLOAT AtA_xx = Am_Xx * Am_Xx + Am_Yx * Am_Yx;
    const RFLOAT AtA_xy = Am_Xx * Am_Xy + Am_Yx * Am_Yy;
    const RFLOAT AtA_yy = Am_Xy * Am_Xy + Am_Yy * Am_Yy;
    const RFLOAT AtA_xy2 = AtA_xy * AtA_xy;

    //#define DEBUG_BACKROTATE
    #ifdef DEBUG_BACKROTATE
    std::cerr << " XSIZE(f2d)= "<< XSIZE(f2d) << std::endl;
    std::cerr << " YSIZE(f2d)= "<< YSIZE(f2d) << std::endl;
    std::cerr << " XSIZE(data)= "<< XSIZE(data) << std::endl;
    std::cerr << " YSIZE(data)= "<< YSIZE(data) << std::endl;
    std::cerr << " STARTINGX(data)= "<< STARTINGX(data) << std::endl;
    std::cerr << " STARTINGY(data)= "<< STARTINGY(data) << std::endl;
    std::cerr << " STARTINGZ(data)= "<< STARTINGZ(data) << std::endl;
    std::cerr << " max_r= "<< r_max << std::endl;
    std::cerr << " Ainv= " << Ainv << std::endl;
    #endif
    const int s  = YSIZE(f2d);
    const int sh = XSIZE(f2d);

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

        if (first_x < first_allowed_x) first_x = first_allowed_x;
        if (last_x > sh - 1) last_x = sh - 1;

        for (int x = first_x; x <= last_x; x++) {

            RFLOAT my_weight = Mweight ? DIRECT_A2D_ELEM(*Mweight, i, x) : 1.0;
            if (my_weight <= 0.0) continue;

            // Get the relevant value in the input image
            Complex my_val = DIRECT_A2D_ELEM(f2d, i, x);

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
                // Subtract STARTINGY to accelerate access to data (STARTINGX=0)
                // In that way use DIRECT_A2D_ELEM, rather than A2D_ELEM
                const int x0 = floor(xp);
                const RFLOAT fx = xp - x0;
                const int x1 = x0 + 1;

                int y0 = floor(yp);
                const RFLOAT fy = yp - y0;
                y0 -=  STARTINGY(data);
                const int y1 = y0 + 1;

                const RFLOAT mfx = 1.0 - fx;
                const RFLOAT mfy = 1.0 - fy;

                const RFLOAT dd00 = mfy * mfx;
                const RFLOAT dd01 = mfy *  fx;
                const RFLOAT dd10 =  fy * mfx;
                const RFLOAT dd11 =  fy *  fx;

                if (is_neg_x) {
                    my_val = conj(my_val);
                }

                // Store slice in 3D weighted sum
                DIRECT_A2D_ELEM(data, y0, x0) += dd00 * my_val;
                DIRECT_A2D_ELEM(data, y0, x1) += dd01 * my_val;
                DIRECT_A2D_ELEM(data, y1, x0) += dd10 * my_val;
                DIRECT_A2D_ELEM(data, y1, x1) += dd11 * my_val;

                // Store corresponding weights
                DIRECT_A2D_ELEM(weight, y0, x0) += dd00 * my_weight;
                DIRECT_A2D_ELEM(weight, y0, x1) += dd01 * my_weight;
                DIRECT_A2D_ELEM(weight, y1, x0) += dd10 * my_weight;
                DIRECT_A2D_ELEM(weight, y1, x1) += dd11 * my_weight;

            } else if (interpolator == NEAREST_NEIGHBOUR) {
                const int x0 = round(xp);
                const int y0 = round(yp);

                if (x0 < 0) {
                    A2D_ELEM(data,   -y0, -x0) += conj(my_val);
                    A2D_ELEM(weight, -y0, -x0) += my_weight;
                } else {
                    A2D_ELEM(data,   y0, x0) += my_val;
                    A2D_ELEM(weight, y0, x0) += my_weight;
                }
            } else {
                REPORT_ERROR("FourierInterpolator::backrotate2D%%ERROR: unrecognized interpolator ");
            }
        }
    }
}

void BackProjector::backrotate3D(
    const MultidimArray<Complex> &f3d,
    const Matrix2D<RFLOAT> &A,
    const MultidimArray<RFLOAT> *Mweight
) {
    // f3d should already be in the right size (ori_size,orihalfdim)
    // AND the points outside max_r should already be zero.

    Matrix2D<RFLOAT> Ainv = A.inv();
    Ainv *= (RFLOAT)padding_factor;  // take scaling into account directly

    const int r_max_src = XSIZE(f3d) - 1;
    const int r_max_src_2 = r_max_src * r_max_src;

    const int r_max_ref = r_max * padding_factor;
    const int r_max_ref_2 = r_max_ref * r_max_ref;

    const int r_min_NN_ref_2 = r_min_nn * r_min_nn * padding_factor * padding_factor;

    //#define DEBUG_BACKROTATE
    #ifdef DEBUG_BACKROTATE
    std::cerr << " XSIZE(f3d)= "<< XSIZE(f3d) << std::endl;
    std::cerr << " YSIZE(f3d)= "<< YSIZE(f3d) << std::endl;
    std::cerr << " XSIZE(data)= "<< XSIZE(data) << std::endl;
    std::cerr << " YSIZE(data)= "<< YSIZE(data) << std::endl;
    std::cerr << " STARTINGX(data)= "<< STARTINGX(data) << std::endl;
    std::cerr << " STARTINGY(data)= "<< STARTINGY(data) << std::endl;
    std::cerr << " STARTINGZ(data)= "<< STARTINGZ(data) << std::endl;
    std::cerr << " max_r= "<< r_max << std::endl;
    std::cerr << " Ainv= " << Ainv << std::endl;
    #endif

    for (int k = 0; k < ZSIZE(f3d); k++) {
        int z, x_min;

        // Don't search beyond square with side max_r
        if (k <= r_max_src) {
            z = k;
            x_min = 0;
        } else {
            z = k - ZSIZE(f3d);
            /// TODO: still check this better in the 3D case!!!
            // x==0 (y,z)-plane is stored twice in the FFTW format. Don't set it twice in BACKPROJECTION!
            x_min = 1;
        }

        int z2 = z * z;

        for (int i = 0; i < YSIZE(f3d); i++) {
            int y = i <= r_max_src ? i : i - YSIZE(f3d);
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
                RFLOAT my_weight = Mweight ? DIRECT_A3D_ELEM(*Mweight, k, i, x) : 1.0;
                if (my_weight <= 0.0) continue;

                Complex my_val = DIRECT_A3D_ELEM(f3d, k, i, x);

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
                    // Subtract STARTINGY to accelerate access to data (STARTINGX=0)
                    // In that way use DIRECT_A3D_ELEM, rather than A3D_ELEM
                    const int x0 = floor(xp);
                    const RFLOAT fx = xp - x0;
                    const int x1 = x0 + 1;

                    int y0 = floor(yp);
                    const RFLOAT fy = yp - y0;
                    y0 -=  STARTINGY(data);
                    const int y1 = y0 + 1;

                    int z0 = floor(zp);
                    const RFLOAT fz = zp - z0;
                    z0 -=  STARTINGZ(data);
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
                    DIRECT_A3D_ELEM(data, z0, y0, x0) += dd000 * my_val;
                    DIRECT_A3D_ELEM(data, z0, y0, x1) += dd001 * my_val;
                    DIRECT_A3D_ELEM(data, z0, y1, x0) += dd010 * my_val;
                    DIRECT_A3D_ELEM(data, z0, y1, x1) += dd011 * my_val;
                    DIRECT_A3D_ELEM(data, z1, y0, x0) += dd100 * my_val;
                    DIRECT_A3D_ELEM(data, z1, y0, x1) += dd101 * my_val;
                    DIRECT_A3D_ELEM(data, z1, y1, x0) += dd110 * my_val;
                    DIRECT_A3D_ELEM(data, z1, y1, x1) += dd111 * my_val;

                    // Store corresponding weights
                    DIRECT_A3D_ELEM(weight, z0, y0, x0) += dd000 * my_weight;
                    DIRECT_A3D_ELEM(weight, z0, y0, x1) += dd001 * my_weight;
                    DIRECT_A3D_ELEM(weight, z0, y1, x0) += dd010 * my_weight;
                    DIRECT_A3D_ELEM(weight, z0, y1, x1) += dd011 * my_weight;
                    DIRECT_A3D_ELEM(weight, z1, y0, x0) += dd100 * my_weight;
                    DIRECT_A3D_ELEM(weight, z1, y0, x1) += dd101 * my_weight;
                    DIRECT_A3D_ELEM(weight, z1, y1, x0) += dd110 * my_weight;
                    DIRECT_A3D_ELEM(weight, z1, y1, x1) += dd111 * my_weight;

                } else if (interpolator == NEAREST_NEIGHBOUR) {
                    const int x0 = round(xp);
                    const int y0 = round(yp);
                    const int z0 = round(zp);

                    if (x0 < 0) {
                        A3D_ELEM(data,   -z0, -y0, -x0) += conj(my_val);
                        A3D_ELEM(weight, -z0, -y0, -x0) += my_weight;
                    } else {
                        A3D_ELEM(data,   z0, y0, x0) += my_val;
                        A3D_ELEM(weight, z0, y0, x0) += my_weight;
                    }
                } else {
                    REPORT_ERROR("BackProjector::backrotate3D%%ERROR: unrecognized interpolator ");
                }
            }
        }
    }
}

void BackProjector::getLowResDataAndWeight(
    MultidimArray<Complex> &lowres_data, MultidimArray<RFLOAT> &lowres_weight,
    int lowres_r_max
) {

    const int lowres_r2_max = round(padding_factor * lowres_r_max) * round(padding_factor * lowres_r_max);
    const int lowres_pad_size = 2 * round(padding_factor * lowres_r_max) + 3;

    // Check lowres_r_max is not too big
    if (lowres_r_max > r_max)
        REPORT_ERROR("BackProjector::getLowResDataAndWeight%%ERROR: lowres_r_max is bigger than r_max");

    // Initialize lowres_data and low_res_weight arrays
    lowres_data.clear();
    lowres_weight.clear();

    if (ref_dim == 2) {
        lowres_data  .resize(lowres_pad_size, lowres_pad_size / 2 + 1);
        lowres_weight.resize(lowres_pad_size, lowres_pad_size / 2 + 1);
    } else {
        lowres_data  .resize(lowres_pad_size, lowres_pad_size, lowres_pad_size / 2 + 1);
        lowres_weight.resize(lowres_pad_size, lowres_pad_size, lowres_pad_size / 2 + 1);
    }
    lowres_data.setXmippOrigin();
    lowres_data.xinit = 0;
    lowres_weight.setXmippOrigin();
    lowres_weight.xinit = 0;

    // fill lowres arrays with relevant values
    FOR_ALL_ELEMENTS_IN_ARRAY3D(lowres_data) {
        if (k * k + i * i + j * j <= lowres_r2_max) {
            A3D_ELEM(lowres_data,   k, i, j) = A3D_ELEM(data,   k , i, j);
            A3D_ELEM(lowres_weight, k, i, j) = A3D_ELEM(weight, k , i, j);
        }
    }
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
        YSIZE(lowres_data) != lowres_pad_size || 
        XSIZE(lowres_data) != lowres_pad_size / 2 + 1 ||
        ZSIZE(lowres_data) != lowres_pad_size && ref_dim == 3
    )
        REPORT_ERROR("BackProjector::setLowResDataAndWeight%%ERROR: lowres_data is not of expected size...");
    if (
        YSIZE(lowres_weight) != lowres_pad_size || 
        XSIZE(lowres_weight) != lowres_pad_size / 2 + 1 ||
        ZSIZE(lowres_weight) != lowres_pad_size && ref_dim == 3
    )
        REPORT_ERROR("BackProjector::setLowResDataAndWeight%%ERROR: lowres_weight is not of expected size...");

    // Re-set origin to the expected place
    lowres_data.setXmippOrigin();
    lowres_data.xinit = 0;
    lowres_weight.setXmippOrigin();
    lowres_weight.xinit = 0;

    // Overwrite data and weight with the lowres arrays
    FOR_ALL_ELEMENTS_IN_ARRAY3D(lowres_data) {
        if (k * k + i * i + j * j <= lowres_r2_max) {
            A3D_ELEM(data,   k, i, j) = A3D_ELEM(lowres_data,   k , i, j);
            A3D_ELEM(weight, k, i, j) = A3D_ELEM(lowres_weight, k , i, j);
        }
    }
}

void BackProjector::getDownsampledAverage(
    MultidimArray<Complex>& avg, bool divide
) const {
    MultidimArray<RFLOAT> down_weight;

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
    avg.setXmippOrigin();
    avg.xinit=0;
    // Resize down_weight the same as down_data
    down_weight.initZeros(avg);

    // Now calculate the down-sized sum
    int kp, ip, jp;
    FOR_ALL_ELEMENTS_IN_ARRAY3D(data) {
        kp = round((RFLOAT) k / padding_factor);
        ip = round((RFLOAT) i / padding_factor);
        jp = round((RFLOAT) j / padding_factor);

        // TMP
        //#define CHECK_SIZE
        #ifdef CHECK_SIZE
        if (
            ip < STARTINGY(avg) || ip > FINISHINGY(avg) 
            jp < STARTINGX(avg) || jp > FINISHINGX(avg) || 
            kp < STARTINGZ(avg) || kp > FINISHINGZ(avg) || 
        ) {
            std::cerr << " kp= " << kp << " ip= " << ip << " jp= " << jp << std::endl;
            avg.printShape();
            REPORT_ERROR("BackProjector::getDownsampledAverage: indices out of range");
        }
        #endif
        A3D_ELEM(avg, kp, ip, jp) += A3D_ELEM(data, k , i, j);
        A3D_ELEM(down_weight, kp, ip, jp) += (divide? A3D_ELEM(weight, k , i, j) : 1.0);
    }

    // Calculate the straightforward average in the downsampled arrays
    FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(avg) {
        if (DIRECT_MULTIDIM_ELEM(down_weight, n) > 0.0) {
            DIRECT_MULTIDIM_ELEM(avg, n) /= DIRECT_MULTIDIM_ELEM(down_weight, n);
        } else {
            DIRECT_MULTIDIM_ELEM(avg, n) = 0.0;
        }
    }
}

void BackProjector::calculateDownSampledFourierShellCorrelation(
    const MultidimArray<Complex>& avg1,
    const MultidimArray<Complex>& avg2,
    MultidimArray<RFLOAT>& fsc
) const {
    if (!avg1.sameShape(avg2))
        REPORT_ERROR("ERROR BackProjector::calculateDownSampledFourierShellCorrelation: two arrays have different sizes");

    MultidimArray<RFLOAT> num, den1, den2;
    num.initZeros(ori_size / 2 + 1);
    den1.initZeros(num);
    den2.initZeros(num);
    fsc.initZeros(num);

    FOR_ALL_ELEMENTS_IN_ARRAY3D(avg1) {
        const RFLOAT R = sqrt(k * k + i * i + j * j);

        if (R > r_max) continue;

        int idx = round(R);

        Complex z1 = A3D_ELEM(avg1, k, i, j);
        Complex z2 = A3D_ELEM(avg2, k, i, j);

        RFLOAT nrmz1 = z1.norm();
        RFLOAT nrmz2 = z2.norm();

        num(idx) += z1.real * z2.real + z1.imag * z2.imag;
        den1(idx) += nrmz1;
        den2(idx) += nrmz2;
    }

    FOR_ALL_ELEMENTS_IN_ARRAY1D(fsc) {
        if (den1(i) * den2(i) > 0.0) {
            fsc(i) = num(i) / sqrt(den1(i) * den2(i));
        }
    }

    // Always set zero-resolution shell to FSC=1
    // Raimond Ravelli reported a problem with FSC=1 at res=0 on 13feb2013...
    // (because of a suboptimal normalisation scheme, but anyway)
    fsc(0) = 1.0;
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
    MultidimArray<RFLOAT> sigma2, data_vs_prior, fourier_coverage;
    MultidimArray<RFLOAT> tau2 = tau2_io;
    MultidimArray<RFLOAT> counter;
    const int max_r2 = round(r_max * padding_factor) * round(r_max * padding_factor);
    RFLOAT oversampling_correction = ref_dim == 3 ? padding_factor * padding_factor * padding_factor : padding_factor * padding_factor;

    // First calculate the radial average of the (inverse of the) power of the noise in the reconstruction
    // This is the left-hand side term in the nominator of the Wiener-filter-like update formula
    // and it is stored inside the weight vector
    // Then, if (do_map) add the inverse of tau2-spectrum values to the weight
    sigma2.initZeros(ori_size / 2 + 1);
    counter.initZeros(ori_size / 2 + 1);
    FOR_ALL_ELEMENTS_IN_ARRAY3D(weight) {
        const int r2 = k * k + i * i + j * j;
        if (r2 < max_r2) {
            int ires = round(sqrt((RFLOAT) r2) / padding_factor);
            RFLOAT invw = oversampling_correction * A3D_ELEM(weight, k, i, j);
            DIRECT_A1D_ELEM(sigma2, ires) += invw;
            DIRECT_A1D_ELEM(counter, ires) += 1.0;
        }
    }

    // Average (inverse of) sigma2 in reconstruction
    FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY1D(sigma2) {
        double x = DIRECT_A1D_ELEM(sigma2, i);
        if (x > 1e-10) {
            x = DIRECT_A1D_ELEM(counter, i) / x;
        } else if (x == 0) {
            x = 0.0;
        } else {
            std::cerr << " DIRECT_A1D_ELEM(sigma2, i)= " << x << std::endl;
            REPORT_ERROR("BackProjector::reconstruct: ERROR: unexpectedly small, yet non-zero sigma2 value, this should not happen...a");
        }
    }

    tau2.reshape(ori_size / 2 + 1);
    data_vs_prior.initZeros(ori_size / 2 + 1);
    fourier_coverage.initZeros(ori_size / 2 + 1);
    counter.initZeros(ori_size / 2 + 1);
    if (update_tau2_with_fsc) {
        // Then calculate new tau2 values, based on the FSC
        if (!fsc.sameShape(sigma2) || !fsc.sameShape(tau2)) {
            fsc.printShape(std::cerr);
            tau2.printShape(std::cerr);
            sigma2.printShape(std::cerr);
            REPORT_ERROR("ERROR BackProjector::reconstruct: sigma2, tau2 and fsc have different sizes");
        }
        FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY1D(sigma2) {
            // FSC cannot be negative or zero for conversion into tau2
            RFLOAT myfsc = std::max(0.001, DIRECT_A1D_ELEM(fsc, i));
            if (iswhole) {
                // Factor two because of twice as many particles
                // Sqrt-term to get 60-degree phase errors....
                myfsc = sqrt(2.0 * myfsc / (myfsc + 1.0));
            }
            myfsc = std::min(0.999, myfsc);
            RFLOAT myssnr = myfsc / (1.0 - myfsc);
            // Sjors 29nov2017 try tau2_fudge for pulling harder on Refine3D runs...
            myssnr *= tau2_fudge;
            RFLOAT fsc_based_tau = myssnr * DIRECT_A1D_ELEM(sigma2, i);
            DIRECT_A1D_ELEM(tau2, i) = fsc_based_tau;
            // data_vs_prior is merely for reporting: it is not used for anything in the reconstruction
            DIRECT_A1D_ELEM(data_vs_prior, i) = myssnr;
        }
    }

    // Now accumulate data_vs_prior if (!update_tau2_with_fsc)
    // Also accumulate fourier_coverage
    FOR_ALL_ELEMENTS_IN_ARRAY3D(weight) {
        int r2 = k * k + i * i + j * j;
        if (r2 < max_r2) {
            int ires = round(sqrt((RFLOAT) r2) / padding_factor);
            RFLOAT invw = A3D_ELEM(weight, k, i, j);

            RFLOAT invtau2;
            if (DIRECT_A1D_ELEM(tau2, ires) > 0.0) {
                // Calculate inverse of tau2
                invtau2 = 1.0 / (oversampling_correction * tau2_fudge * DIRECT_A1D_ELEM(tau2, ires));
            } else if (DIRECT_A1D_ELEM(tau2, ires) == 0.0) {
                // If tau2 is zero, use small value instead
                invtau2 = 1.0 / (0.001 * invw);
            } else {
                std::cerr << " sigma2= " << sigma2 << std::endl;
                std::cerr << " fsc= " << fsc << std::endl;
                std::cerr << " tau2= " << tau2 << std::endl;
                REPORT_ERROR("ERROR BackProjector::reconstruct: Negative or zero values encountered for tau2 spectrum!");
            }

            // Keep track of spectral evidence-to-prior ratio and remaining noise in the reconstruction
            if (!update_tau2_with_fsc) {
                DIRECT_A1D_ELEM(data_vs_prior, ires) += invw / invtau2;
            }

            // Keep track of the coverage in Fourier space
            if (invw / invtau2 >= 1.0) {
                DIRECT_A1D_ELEM(fourier_coverage, ires) += 1.0;
            }

            DIRECT_A1D_ELEM(counter, ires) += 1.0;

        }
    }

    // Average data_vs_prior
    if (!update_tau2_with_fsc) {
        for (long int i = 0; i < data_vs_prior.xdim; i++) {
            RFLOAT x = DIRECT_A1D_ELEM(data_vs_prior, i);
            RFLOAT n = DIRECT_A1D_ELEM(counter, i);
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
    FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY1D(fourier_coverage) {
        if (DIRECT_A1D_ELEM(counter, i) > 0.0)
            DIRECT_A1D_ELEM(fourier_coverage, i) /= DIRECT_A1D_ELEM(counter, i);
    }

    // Send back the output
    tau2_io = tau2;
    sigma2_out = sigma2;
    data_vs_prior_out = data_vs_prior;
    fourier_coverage_out = fourier_coverage;
}

void BackProjector::externalReconstruct(
    MultidimArray<RFLOAT> &vol_out,
    FileName &fn_out,
    MultidimArray<RFLOAT> &fsc_halves_io,
    MultidimArray<RFLOAT> &tau2_io,
    MultidimArray<RFLOAT> &sigma2_ref,
    MultidimArray<RFLOAT> &data_vs_prior,
    bool iswhole,
    RFLOAT tau2_fudge,
    int verb
) {

    FileName fn_recons = fn_out + "_external_reconstruct.mrc";
    FileName fn_star = fn_out + "_external_reconstruct.star";
    FileName fn_out_star = fn_out + "_external_reconstruct_out.star";
    MultidimArray<RFLOAT> fsc_halves = fsc_halves_io;
    MultidimArray<RFLOAT> tau2 = tau2_io;

    const int max_r2 = round(r_max * padding_factor) * round(r_max * padding_factor);
    int padoridim = round(padding_factor * ori_size);

    // Write out data array
    Image<Complex> Idata;
    if (ref_dim == 2) {
        Idata().resize(pad_size, pad_size / 2 + 1);
    } else {
        Idata().resize(pad_size, pad_size, pad_size / 2 + 1);
    }
    Projector::decenter(data, Idata(), max_r2);
    windowFourierTransform(Idata(), padoridim);
    ComplexIO::write(Idata(), fn_out + "_external_reconstruct_data", ".mrc");
    Idata.clear();

    // Write out weight array
    Image<RFLOAT> Iweight;
    if (ref_dim == 2) {
        Iweight().resize(pad_size, pad_size / 2 + 1);
    } else {
        Iweight().resize(pad_size, pad_size, pad_size / 2 + 1);
    }
    Projector::decenter(weight, Iweight(), max_r2);
    windowFourierTransform(Iweight(), padoridim);
    Iweight.write(fn_out+"_external_reconstruct_weight.mrc");
    Iweight.clear();

    // Write out STAR file for input to external reconstruction program
    MetaDataTable MDlist, MDtau;

    MDlist.setName("external_reconstruct_general");
    MDlist.isList = true;
    MDlist.addObject();
    MDlist.setValue(EMDL::OPTIMISER_EXTERNAL_RECONS_DATA_REAL, fn_out + "_external_reconstruct_data_real.mrc");
    MDlist.setValue(EMDL::OPTIMISER_EXTERNAL_RECONS_DATA_IMAG, fn_out + "_external_reconstruct_data_imag.mrc");
    MDlist.setValue(EMDL::OPTIMISER_EXTERNAL_RECONS_WEIGHT, fn_out + "_external_reconstruct_weight.mrc");
    MDlist.setValue(EMDL::OPTIMISER_EXTERNAL_RECONS_RESULT, fn_recons);
    MDlist.setValue(EMDL::OPTIMISER_EXTERNAL_RECONS_NEWSTAR, fn_out_star);
    MDlist.setValue(EMDL::MLMODEL_TAU2_FUDGE_FACTOR, tau2_fudge);
    MDlist.setValue(EMDL::MLMODEL_PADDING_FACTOR, padding_factor);
    MDlist.setValue(EMDL::MLMODEL_DIMENSIONALITY, ref_dim);
    MDlist.setValue(EMDL::MLMODEL_ORIGINAL_SIZE, ori_size);
    MDlist.setValue(EMDL::MLMODEL_CURRENT_SIZE, 2 * r_max);

    MDtau.setName("external_reconstruct_tau2");
    for (int ii = 0; ii < XSIZE(tau2); ii++) {
        MDtau.addObject();
        MDtau.setValue(EMDL::SPECTRAL_IDX, ii);
        MDtau.setValue(EMDL::MLMODEL_TAU2_REF, tau2(ii));
        MDtau.setValue(EMDL::MLMODEL_FSC_HALVES_REF, fsc_halves(ii));
    }

    {
        std::ofstream fh(fn_star.c_str(), std::ios::out);
        if (!fh)
            REPORT_ERROR((std::string)"BackProjector::externalReconstruct: Cannot write file: " + fn_star);
        MDlist.write(fh);
        MDtau.write(fh);
    }

    // Make the system call: program name plus the STAR file for the external reconstruction program as its first argument
    const char *my_exec = getenv ("RELION_EXTERNAL_RECONSTRUCT_EXECUTABLE");
    if (!my_exec) { my_exec = DEFAULT_EXTERNAL_RECONSTRUCT; }
    std::string command = std::string(my_exec) + " " + fn_star;

    if (verb > 0)
        std::cout << std::endl << " + Making system call for external reconstruction: " << command << std::endl;

    int res = system(command.c_str());
    if (res) {
        REPORT_ERROR(" ERROR: there was something wrong with system call: " + command);
    } else if (verb > 0) {
        std::cout << " + External reconstruction finished successfully, reading result back in ... " << std::endl;
    }

    // Read the resulting map back into memory
    Iweight.read(fn_recons);
    vol_out = Iweight();
    vol_out.setXmippOrigin();

    if (exists(fn_out_star)) {

        MetaDataTable MDnewtau;
        MDnewtau.read(fn_out_star);

        if (!MDnewtau.containsLabel(EMDL::SPECTRAL_IDX))
            REPORT_ERROR("ERROR: external reconstruct output STAR file does not contain spectral idx!");

        // Directly update tau2 spectrum
        FOR_ALL_OBJECTS_IN_METADATA_TABLE(MDnewtau) {
            int idx = MDnewtau.getValue<int>(EMDL::SPECTRAL_IDX);

            if (idx >= XSIZE(tau2_io)) continue;

            if (MDnewtau.containsLabel(EMDL::MLMODEL_TAU2_REF)) {
                tau2_io(idx)       = MDnewtau.getValue<RFLOAT>(EMDL::MLMODEL_TAU2_REF);
                data_vs_prior(idx) = tau2_io(idx) / sigma2_ref(idx);
            } else if (MDnewtau.containsLabel(EMDL::POSTPROCESS_FSC_GENERAL)) {
                idx                = MDnewtau.getValue<int>(EMDL::SPECTRAL_IDX);
                fsc_halves_io(idx) = MDnewtau.getValue<RFLOAT>(EMDL::POSTPROCESS_FSC_GENERAL);

                RFLOAT myfsc = std::max(0.001, fsc_halves_io(idx));
                if (iswhole) {
                    // Factor two because of twice as many particles
                    // Sqrt-term to get 60-degree phase errors....
                    myfsc = sqrt(2.0 * myfsc / (myfsc + 1.0));
                }
                myfsc = std::min(0.999, myfsc);
                RFLOAT myssnr = myfsc / (1.0 - myfsc);
                myssnr *= tau2_fudge;
                tau2_io(idx) = myssnr * sigma2_ref(idx);
                data_vs_prior(idx) = myssnr;
            } else {
                REPORT_ERROR("ERROR: output STAR file from external reconstruct does not contain tau2 or FSC array");
            }
        }

        if (verb > 0)
            std::cout << " + External reconstruction successfully updated external tau2 array ... " << std::endl;
    }

}

void BackProjector::reconstruct(
    MultidimArray<RFLOAT> &vol_out,
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
    RCTICTOC(ReconTimer, ReconS[0], ({
    oversampling_correction = ref_dim == 3 ?
        padding_factor * padding_factor * padding_factor : 
        padding_factor * padding_factor;

    // #define DEBUG_RECONSTRUCT
    #ifdef DEBUG_RECONSTRUCT
    Image<RFLOAT> ttt;
    FileName fnttt;
    ttt() = weight;
    ttt.write("reconstruct_initial_weight.spi");
    std::cerr << " pad_size= " << pad_size << " padding_factor= " << padding_factor << " max_r2= " << max_r2 << std::endl;
    #endif

    // Set Fconv to the right size
    if (ref_dim == 2) {
        vol_out.setDimensions(pad_size, pad_size, 1, 1);
    } else {
        // Too costly to actually allocate the space
        // Trick transformer with the right dimensions
        vol_out.setDimensions(pad_size, pad_size, pad_size, 1);
    }

    transformer.setReal(vol_out);  // Fake set real. 1. Allocate space for Fconv 2. calculate plans.
    }));

    /// TODO: Find some way of putting these two lines back inside the above block.
    MultidimArray<Complex> &Fconv = transformer.getFourierReference();
    vol_out.clear();  // Reset dimensions to 0

    MultidimArray<RFLOAT> Fweight;
    RCTICTOC(ReconTimer, ReconS[1], ({
    // Go from projector-centered to FFTW-uncentered
    Fweight.reshape(Fconv);
    Projector::decenter(weight, Fweight, max_r2);
    }))

    RCTICTOC(ReconTimer, ReconS[2], ({
    // Apply MAP-additional term to the Fweight array
    // This will regularise the actual reconstruction
    if (do_map) {
        // Then, add the inverse of tau2-spectrum values to the weight
        FOR_ALL_ELEMENTS_IN_FFTW_TRANSFORM(Fconv) {
            int r2 = kp * kp + ip * ip + jp * jp;
            if (r2 < max_r2) {
                int ires = round(sqrt((RFLOAT) r2) / padding_factor);
                RFLOAT invw = DIRECT_A3D_ELEM(Fweight, k, i, j);

                RFLOAT invtau2;
                if (DIRECT_A1D_ELEM(tau2, ires) > 0.0) {
                    // Calculate inverse of tau2
                    invtau2 = 1.0 / (oversampling_correction * tau2_fudge * DIRECT_A1D_ELEM(tau2, ires));
                } else if (DIRECT_A1D_ELEM(tau2, ires) < 1e-20) {
                    // If tau2 is zero, use small value instead
                    invtau2 = invw > 1e-20 ? 1.0 / (0.001 * invw) : 0.0;
                } else {
                    std::cerr << " tau2= " << tau2 << std::endl;
                    REPORT_ERROR("ERROR BackProjector::reconstruct: Negative or zero values encountered for tau2 spectrum!");
                }

                // Only for (ires >= minres_map) add Wiener-filter like term
                if (ires >= minres_map) {
                    // Now add the inverse-of-tau2_class term
                    invw += invtau2;
                    // Store the new weight again in Fweight
                    DIRECT_A3D_ELEM(Fweight, k, i, j) = invw;
                }
            }
        }
    }
    }))

    if (skip_gridding) {
        RCTIC(ReconTimer, ReconS[3]);  // BUG: No RCTOC!
        Fconv.initZeros();  // Clear the input volume
        Projector::decenter(data, Fconv, max_r2);

        // Prevent divisions by zero: set Fweight to at least 1/1000th of the radially averaged weight at that resolution
        // beyond r_max, set Fweight to at least 1/1000th of the radially averaged weight at r_max;
        MultidimArray<RFLOAT> radavg_weight(r_max), counter(r_max);
        const int round_max_r2 = round(r_max * padding_factor * r_max * padding_factor);
        FOR_ALL_ELEMENTS_IN_FFTW_TRANSFORM(Fweight) {
            const int r2 = kp * kp + ip * ip + jp * jp;
            // Note that (r < ires) != (r2 < max_r2), because max_r2 = round(r_max * padding_factor)^2.
            // We have to use round_max_r2 = round((r_max * padding_factor)^2).
            // e.g. k = 0, i = 7, j = 28, max_r2 = 841, r_max = 16, padding_factor = 18.
            if (r2 < round_max_r2) {
                const int ires = floor(sqrt((RFLOAT)r2) / padding_factor);
                if (ires >= XSIZE(radavg_weight)) {
                    std::cerr << " k= " << k << " i= " << i << " j= " << j << std::endl;
                    std::cerr << " ires= " << ires << " XSIZE(radavg_weight)= " << XSIZE(radavg_weight) << std::endl;
                    REPORT_ERROR("BUG: ires >=XSIZE(radavg_weight) ");
                }
                DIRECT_A1D_ELEM(radavg_weight, ires) += DIRECT_A3D_ELEM(Fweight, k, i, j);
                DIRECT_A1D_ELEM(counter, ires) += 1.0;
            }
        }

        // Calculate 1/1000th of radial averaged weight
        FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY1D(radavg_weight) {
            if (
                DIRECT_A1D_ELEM(counter,       i) > 0.0 || 
                DIRECT_A1D_ELEM(radavg_weight, i) > 0.0
            ) {
                DIRECT_A1D_ELEM(radavg_weight, i) /= 1000.0 * DIRECT_A1D_ELEM(counter, i);
            } else {
                std::cerr << " counter= " << counter << std::endl;
                std::cerr << " radavg_weight= " << radavg_weight << std::endl;
                REPORT_ERROR("BUG: zeros in counter or radavg_weight!");
            }
        }

        bool have_warned = false;
        // perform std::max on all weight elements, and do division of data/weight
        FOR_ALL_ELEMENTS_IN_FFTW_TRANSFORM(Fweight) {
            const int r2 = kp * kp + ip * ip + jp * jp;
            const int ires = floor(sqrt((RFLOAT) r2) / padding_factor);
            const RFLOAT weight =  std::max(
                DIRECT_A3D_ELEM(Fweight, k, i, j), 
                DIRECT_A1D_ELEM(radavg_weight, ires < r_max ? ires : r_max - 1)
            );
            if (weight == 0) {
                if (!have_warned) {
                    std::cerr << " WARNING: ignoring divide by zero in skip_gridding: ires = " << ires << " kp = " << kp << " ip = " << ip << " jp = " << jp << std::endl;
                    std::cerr << " max_r2 = " << max_r2 << " r_max = " << r_max << " padding_factor = " << padding_factor
                              << " round(sqrt(max_r2)) = " << round(sqrt(max_r2)) << " round(r_max * padding_factor) = " << round(r_max * padding_factor) << std::endl;
                    have_warned = true;
                }
            } else {
                DIRECT_A3D_ELEM(Fconv, k, i, j) /= weight;
            }
        }
    } else {

        RCTICTOC(ReconTimer, ReconS[4], ({
        // Divide both data and Fweight by normalisation factor to prevent FFT's with very large values....
        #ifdef DEBUG_RECONSTRUCT
        std::cerr << " normalise= " << normalise << std::endl;
        #endif
        FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(Fweight) {
            DIRECT_MULTIDIM_ELEM(Fweight, n) /= normalise;
        }
        FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(data) {
            DIRECT_MULTIDIM_ELEM(data, n) /= normalise;
        }
        }))

        MultidimArray<double> Fnewweight;
        RCTICTOC(ReconTimer, ReconS[5], ({
        // Initialise Fnewweight with 1's and 0's. (also see comments below)
        FOR_ALL_ELEMENTS_IN_ARRAY3D(weight) {
            A3D_ELEM(weight, k, i, j) = k * k + i * i + j * j < max_r2;
        }
        // Fnewweight can become too large for a float: always keep this one in double-precision
        Fnewweight.reshape(Fconv);
        decenter(weight, Fnewweight, max_r2);
        }))

        // Iterative algorithm as in Eq. 14 in Pipe & Menon (1999)
        // or Eq. 4 in Matej (2001)
        for (int iter = 0; iter < max_iter_preweight; iter++) {
            // std::cout << "    iteration " << iter + 1 << "/" << max_iter_preweight << "\n";
            RCTICTOC(ReconTimer, ReconS[6], ({
            // Set Fnewweight * Fweight in the transformer
            // In Matej et al (2001), weights w_P^i are convoluted with the kernel,
            // and the initial w_P^0 are 1 at each sampling point
            // Here the initial weights are also 1 (see initialisation Fnewweight above),
            // but each "sampling point" counts "Fweight" times!
            // That is why Fnewweight is multiplied by Fweight prior to the convolution

            FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(Fconv) {
                DIRECT_MULTIDIM_ELEM(Fconv, n) = DIRECT_MULTIDIM_ELEM(Fnewweight, n) * DIRECT_MULTIDIM_ELEM(Fweight, n);
            }

            // convolute through Fourier-transform (as both grids are rectangular)
            // Note that convoluteRealSpace acts on the complex array inside the transformer
            convoluteBlobRealSpace(transformer, false);

            RFLOAT w, corr_min = LARGE_NUMBER, corr_max = -LARGE_NUMBER, corr_avg = 0.0, corr_nn = 0.0;

            FOR_ALL_ELEMENTS_IN_FFTW_TRANSFORM(Fconv) {
                if (kp * kp + ip * ip + jp * jp < max_r2) {

                    // Make sure no division by zero can occur....
                    w = std::max(1e-6, abs(DIRECT_A3D_ELEM(Fconv, k, i, j)));
                    // Monitor min, max and avg conv_weight
                    corr_min = std::min(corr_min, w);
                    corr_max = std::max(corr_max, w);
                    corr_avg += w;
                    corr_nn += 1.0;
                    // Apply division of Eq. [14] in Pipe & Menon (1999)
                    DIRECT_A3D_ELEM(Fnewweight, k, i, j) /= w;
                }
            }

            }))

            #ifdef DEBUG_RECONSTRUCT
            std::cerr << " PREWEIGHTING ITERATION: " << iter + 1 << " OF " << max_iter_preweight << std::endl;
            // report of maximum and minimum values of current conv_weight
            std::cerr << " corr_avg= " << corr_avg / corr_nn << std::endl;
            std::cerr << " corr_min= " << corr_min << std::endl;
            std::cerr << " corr_max= " << corr_max << std::endl;
            #endif
        }

        RCTICTOC(ReconTimer, ReconS[7], ({

        #ifdef DEBUG_RECONSTRUCT
        Image<double> tttt;
        tttt()=Fnewweight;
        tttt.write("reconstruct_gridding_weight.spi");
        FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(Fconv) {
            DIRECT_MULTIDIM_ELEM(ttt(), n) = abs(DIRECT_MULTIDIM_ELEM(Fconv, n));
        }
        ttt.write("reconstruct_gridding_correction_term.spi");
        #endif

        // Clear memory
        Fweight.clear();

        // Note that Fnewweight now holds the approximation of the inverse of the weights on a regular grid

        // Now do the actual reconstruction with the data array
        // Apply the iteratively determined weight
        Fconv.initZeros();  // Clear the input volume
        Projector::decenter(data, Fconv, max_r2);
        FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(Fconv) {
            #ifdef RELION_SINGLE_PRECISION
            // Prevent numerical instabilities in single-precision reconstruction with very unevenly sampled orientations
            #define BIGNUM 1e20
            if (DIRECT_MULTIDIM_ELEM(Fnewweight, n) > BIGNUM) { 
                DIRECT_MULTIDIM_ELEM(Fnewweight, n) = BIGNUM; 
            }
            #undef BIGNUM
            #endif
            DIRECT_MULTIDIM_ELEM(Fconv, n) *= DIRECT_MULTIDIM_ELEM(Fnewweight, n);
        }

        // Clear memory
        Fnewweight.clear();

        }))
    }

    // Gridding theory says one now has to interpolate the fine grid onto the coarse one using a blob kernel
    // and then do the inverse transform and divide by the FT of the blob (i.e. do the gridding correction)
    // In practice, this gives all types of artefacts (perhaps I never found the right implementation?!)
    // Therefore, window the Fourier transform and then do the inverse transform
    //#define RECONSTRUCT_CONVOLUTE_BLOB
    #ifdef RECONSTRUCT_CONVOLUTE_BLOB

    // Apply the same blob-convolution as above to the data array
    // Mask real-space map beyond its original size to prevent aliasing in the downsampling step below
    RCTICTOC(ReconTimer, ReconS[8], ({
    convoluteBlobRealSpace(transformer, true);
    }))

    RCTICTOC(ReconTimer, ReconS[9], ({
    // Now just pick every 3rd pixel in Fourier-space (i.e. down-sample)
    // and do a final inverse FT
    if (ref_dim == 2) {
        vol_out.resize(ori_size, ori_size);
    } else {
        vol_out.resize(ori_size, ori_size, ori_size);
    }
    }))

    RCTICTOC(ReconTimer, ReconS[10]], ({
    FourierTransformer transformer2;
    transformer2.setReal(vol_out);  // cannot use the first transformer because Fconv is inside there!!
    MultidimArray<Complex> Ftmp;
    transformer2.getFourierAlias(Ftmp);
    }))

    RCTICTOC(ReconTimer, ReconS[11], ({
    FOR_ALL_ELEMENTS_IN_FFTW_TRANSFORM(Ftmp) {
        if (kp * kp + ip * ip + jp * jp < r_max * r_max) {
            DIRECT_A3D_ELEM(Ftmp, k, i, j) = FFTW_ELEM(Fconv, kp * padding_factor, ip * padding_factor, jp * padding_factor);
        } else {
            DIRECT_A3D_ELEM(Ftmp, k, i, j) = 0.0;
        }
    }
    }))

    // Any particular reason why 13 comes before 12?

    RCTICTOC(ReconTimer, ReconS[13], ({
    CenterFFTbySign(Ftmp);
    }))

    RCTICTOC(ReconTimer, ReconS[12], ({
    // inverse FFT leaves result in vol_out
    transformer2.inverseFourierTransform();
    }))

    RCTICTOC(ReconTimer, ReconS[14], ({
    // Un-normalize FFTW (because original FFTs were done with the size of 2D FFTs)
    if (ref_dim == 3) { vol_out /= ori_size; }
    }))

    RCTICTOC(ReconTimer, ReconS[15], ({
    // Mask out corners to prevent aliasing artefacts
    softMaskOutsideMap(vol_out);
    }))

    RCTICTOC(ReconTimer, ReconS[16], ({
    // Gridding correction for the blob
    RFLOAT normftblob = tab_ftblob(0.0);
    FOR_ALL_ELEMENTS_IN_ARRAY3D(vol_out) {

        RFLOAT r = sqrt((RFLOAT) (k * k + i * i + j * j));
        RFLOAT rval = r / (ori_size * padding_factor);
        A3D_ELEM(vol_out, k, i, j) /= tab_ftblob(rval) / normftblob;
        // if (k == 0 && i == 0)
        // 	std::cerr << " j= " << j << " rval= " << rval << " tab_ftblob(rval) / normftblob= " << tab_ftblob(rval) / normftblob << std::endl;
    }
    }))

    #else

    // rather than doing the blob-convolution to downsample the data array, do a windowing operation:
    // This is the same as convolution with a SINC. It seems to give better maps.
    // Then just make the blob look as much as a SINC as possible....
    // The "standard" r1.9, m2 and a15 blob looks quite like a sinc until the first zero (perhaps that's why it is standard?)
    //for (RFLOAT r = 0.1; r < 10.; r+=0.01)
    //{
    //	RFLOAT sinc = sin(PI * r / padding_factor ) / ( PI * r / padding_factor);
    //	std::cout << " r= " << r << " sinc= " << sinc << " blob= " << blob_val(r, blob) << std::endl;
    //}

    // Now do inverse FFT and window to original size in real-space
    // Pass the transformer to prevent making and clearing a new one before clearing the one declared above....
    // The latter may give memory problems as detected by electric fence....
    RCTICTOC(ReconTimer, ReconS[17], ({
    windowToOridimRealSpace(transformer, vol_out, printTimes);
    }))

    #endif

    #ifdef DEBUG_RECONSTRUCT
    ttt() = vol_out;
    ttt.write("reconstruct_before_gridding_correction.spi");
    #endif

    // Correct for the linear/nearest-neighbour interpolation that led to the data array
    RCTICTOC(ReconTimer, ReconS[18], ({ griddingCorrect(vol_out); }))

    RCTICTOC(ReconTimer, ReconS[23], ({
    // Completely empty the transformer object
    transformer.cleanup();
    // Now can use extra mem to move data into smaller array space
    vol_out.shrinkToFit();
    }))

    #ifdef TIMING
    if (printTimes) ReconTimer.printTimes(true);
    #endif


    #ifdef DEBUG_RECONSTRUCT
    std::cerr << "done with reconstruct" << std::endl;
    #endif

    if (weight_out != 0) {
        weight_out->data = MultidimArray<RFLOAT>(1, ori_size, ori_size, ori_size / 2 + 1);

        Image<RFLOAT> count(ori_size / 2 + 1, ori_size, ori_size);
        count.data.initZeros();

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
            ) {
                continue;
            }

            int xx =        round(xl / padding_factor);
            int yy = ((int) round(yl / padding_factor) + ori_size) % ori_size;
            int zz = ((int) round(zl / padding_factor) + ori_size) % ori_size;

            if (
                   xx >= 0 && xx < ori_size / 2 + 1
                && yy >= 0 && yy < ori_size
                && zz >= 0 && zz < ori_size
            ) {
                DIRECT_A3D_ELEM(weight_out->data, zz, yy, xx) += DIRECT_A3D_ELEM(Fweight, z, y, x);
                DIRECT_A3D_ELEM(count.data, zz, yy, xx) += 1.0;
            }
        }

        const double pad3 = padding_factor * padding_factor * padding_factor;

        for (long int z = 0; z < ori_size; z++)
        for (long int y = 0; y < ori_size; y++)
        for (long int x = 0; x < ori_size / 2 + 1; x++) {
            const RFLOAT c = DIRECT_A3D_ELEM(count.data, z, y, x);

            if (c > 0.0) {
                DIRECT_A3D_ELEM(weight_out->data, z, y, x) *= pad3 / c;
            }
        }
    }
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
    for (int iz = STARTINGZ(data); iz <=FINISHINGZ(data); iz++) {
        // Make sure all points are only included once.
        for (int iy = iz >= 0; iy <= FINISHINGY(data); iy++) {
            // I just need to sum the two points, not divide by 2!
            Complex fsum = A3D_ELEM(data, iz, iy, 0) + conj(A3D_ELEM(data, -iz, -iy, 0));
            A3D_ELEM(data,  iz,  iy, 0) = fsum;
            A3D_ELEM(data, -iz, -iy, 0) = conj(fsum);
            RFLOAT sum = A3D_ELEM(weight, iz, iy, 0) + A3D_ELEM(weight, -iz, -iy, 0);
            A3D_ELEM(weight,  iz,  iy, 0) = sum;
            A3D_ELEM(weight, -iz, -iy, 0) = sum;
        }
    }
}

void BackProjector::applyHelicalSymmetry(
    int nr_helical_asu, RFLOAT helical_twist, RFLOAT helical_rise
) {
    if (nr_helical_asu < 2 || ref_dim != 3)
        return;

    int rmax2 = round(r_max * padding_factor) * round(r_max * padding_factor);

    Matrix2D<RFLOAT> R(4, 4);  // A matrix from the list
    MultidimArray<RFLOAT> sum_weight;
    MultidimArray<Complex> sum_data;
    bool is_neg_x;

    // First symmetry operator (not stored in SL) is the identity matrix
    sum_weight = weight;
    sum_data = data;
    int h_min = -nr_helical_asu / 2;
    int h_max = -h_min + nr_helical_asu % 2;
    for (int hh = h_min; hh < h_max; hh++) {
        if (hh != 0) {
            // h == 0 is done before the for loop (where sum_data = data)
            RFLOAT rot_ang = hh * (-helical_twist);
            rotation3DMatrix(rot_ang, 'Z', R);
            R.setSmallValuesToZero();  // TODO: invert rotation matrix?

            // Loop over all points in the output (i.e. rotated, or summed) array
            FOR_ALL_ELEMENTS_IN_ARRAY3D(sum_weight) {
                RFLOAT x = j;  // STARTINGX(sum_weight) is zero!
                RFLOAT y = i;
                RFLOAT z = k;
                RFLOAT r2 = x * x + y * y + z * z;
                if (r2 <= rmax2) {
                    // coords_output(x, y) = A * coords_input (xp, yp)
                    RFLOAT xp = x * R(0, 0) + y * R(0, 1) + z * R(0, 2);
                    RFLOAT yp = x * R(1, 0) + y * R(1, 1) + z * R(1, 2);
                    RFLOAT zp = x * R(2, 0) + y * R(2, 1) + z * R(2, 2);

                    is_neg_x = xp < 0;
                    // Only asymmetric half is stored
                    if (is_neg_x) {
                        // Get complex conjugated hermitian symmetry pair
                        xp = -xp;
                        yp = -yp;
                        zp = -zp;
                    }

                    // Trilinear interpolation (with physical coords)
                    // Subtract STARTINGY and STARTINGZ to accelerate access to data (STARTINGX=0)
                    // In that way use DIRECT_A3D_ELEM, rather than A3D_ELEM
                    int x0 = floor(xp);
                    RFLOAT fx = xp - x0;
                    int x1 = x0 + 1;

                    int y0 = floor(yp);
                    RFLOAT fy = yp - y0;
                    y0 -=  STARTINGY(data);
                    int y1 = y0 + 1;

                    int z0 = floor(zp);
                    RFLOAT fz = zp - z0;
                    z0 -= STARTINGZ(data);
                    int z1 = z0 + 1;

                    #ifdef CHECK_SIZE
                    if (
                        x0 < 0 || y0 < 0 || z0 < 0 ||
                        x1 < 0 || y1 < 0 || z1 < 0 ||
                        x0 >= XSIZE(data) || y0 >= YSIZE(data) || z0 >= ZSIZE(data) ||
                        x1 >= XSIZE(data) || y1 >= YSIZE(data) || z1 >= ZSIZE(data)
                    ) {
                        std::cerr << " x0= " << x0 << " y0= " << y0 << " z0= " << z0 << std::endl;
                        std::cerr << " x1= " << x1 << " y1= " << y1 << " z1= " << z1 << std::endl;
                        data.printShape();
                        REPORT_ERROR("BackProjector::applyPointGroupSymmetry: checksize!!!");
                    }
                    #endif
                    // First interpolate (complex) data
                    Complex d000 = DIRECT_A3D_ELEM(data, z0, y0, x0);
                    Complex d001 = DIRECT_A3D_ELEM(data, z0, y0, x1);
                    Complex d010 = DIRECT_A3D_ELEM(data, z0, y1, x0);
                    Complex d011 = DIRECT_A3D_ELEM(data, z0, y1, x1);
                    Complex d100 = DIRECT_A3D_ELEM(data, z1, y0, x0);
                    Complex d101 = DIRECT_A3D_ELEM(data, z1, y0, x1);
                    Complex d110 = DIRECT_A3D_ELEM(data, z1, y1, x0);
                    Complex d111 = DIRECT_A3D_ELEM(data, z1, y1, x1);

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
                        RFLOAT dotp = 2 * PI * z * zshift;
                        Complex X = Complex(cos(dotp), sin(dotp));  // cis(dotp)
                        Complex Z = Complex(X.real * ddd.real, X.imag * ddd.imag);  // X dot ddd
                        ddd = Complex(
                            Z.real - Z.imag, 
                            (X.real + X.imag) * (ddd.real + ddd.imag) - Z.real - Z.imag
                        );
                    }
                    // Accumulated sum of the data term
                    A3D_ELEM(sum_data, k, i, j) += ddd;

                    // Then interpolate (real) weight
                    RFLOAT dd000 = DIRECT_A3D_ELEM(weight, z0, y0, x0);
                    RFLOAT dd001 = DIRECT_A3D_ELEM(weight, z0, y0, x1);
                    RFLOAT dd010 = DIRECT_A3D_ELEM(weight, z0, y1, x0);
                    RFLOAT dd011 = DIRECT_A3D_ELEM(weight, z0, y1, x1);
                    RFLOAT dd100 = DIRECT_A3D_ELEM(weight, z1, y0, x0);
                    RFLOAT dd101 = DIRECT_A3D_ELEM(weight, z1, y0, x1);
                    RFLOAT dd110 = DIRECT_A3D_ELEM(weight, z1, y1, x0);
                    RFLOAT dd111 = DIRECT_A3D_ELEM(weight, z1, y1, x1);

                    RFLOAT ddx00 = LIN_INTERP(fx, dd000, dd001);
                    RFLOAT ddx01 = LIN_INTERP(fx, dd100, dd101);
                    RFLOAT ddx10 = LIN_INTERP(fx, dd010, dd011);
                    RFLOAT ddx11 = LIN_INTERP(fx, dd110, dd111);
                    RFLOAT ddxy0 = LIN_INTERP(fy, ddx00, ddx10);
                    RFLOAT ddxy1 = LIN_INTERP(fy, ddx01, ddx11);

                    A3D_ELEM(sum_weight, k, i, j) +=  LIN_INTERP(fz, ddxy0, ddxy1);

                }
            }       
        }
    }

    data = sum_data;
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
        Matrix2D<RFLOAT> L(4, 4), R(4, 4);  // A matrix from the list
        MultidimArray<RFLOAT> sum_weight;
        MultidimArray<Complex> sum_data;

        // First symmetry operator (not stored in SL) is the identity matrix
        sum_weight = weight;
        sum_data = data;
        // Loop over all other symmetry operators
        for (int isym = 0; isym < SL.SymsNo(); isym++) {
            SL.get_matrices(isym, L, R);
            #ifdef DEBUG_SYMM
            std::cerr << " isym= " << isym << " R= " << R << std::endl;
            #endif
            // Loop over all points in the output (i.e. rotated, or summed) array

            #pragma omp parallel for num_threads(threads)
            for (long int k = STARTINGZ(sum_weight); k <= FINISHINGZ(sum_weight); k++)
            for (long int i = STARTINGY(sum_weight); i <= FINISHINGY(sum_weight); i++)
            for (long int j = STARTINGX(sum_weight); j <= FINISHINGX(sum_weight); j++) {
                RFLOAT x = j;  // STARTINGX(sum_weight) is zero!
                RFLOAT y = i;
                RFLOAT z = k;
                RFLOAT r2 = x * x + y * y + z * z;

                if (r2 <= rmax2) {
                    // coords_output(x,y) = A * coords_input (xp,yp)
                    RFLOAT xp = x * R(0, 0) + y * R(0, 1) + z * R(0, 2);
                    RFLOAT yp = x * R(1, 0) + y * R(1, 1) + z * R(1, 2);
                    RFLOAT zp = x * R(2, 0) + y * R(2, 1) + z * R(2, 2);

                    bool is_neg_x = xp < 0;

                    // Only asymmetric half is stored
                    if (is_neg_x) {
                        // Get complex conjugated hermitian symmetry pair
                        xp = -xp;
                        yp = -yp;
                        zp = -zp;
                    }

                    // Trilinear interpolation (with physical coords)
                    // Subtract STARTINGY and STARTINGZ to accelerate access to data (STARTINGX=0)
                    // In that way use DIRECT_A3D_ELEM, rather than A3D_ELEM
                    int x0 = floor(xp);
                    RFLOAT fx = xp - x0;
                    int x1 = x0 + 1;

                    int y0 = floor(yp);
                    RFLOAT fy = yp - y0;
                    y0 -=  STARTINGY(data);
                    int y1 = y0 + 1;

                    int z0 = floor(zp);
                    RFLOAT fz = zp - z0;
                    z0 -= STARTINGZ(data);
                    int z1 = z0 + 1;

                    #ifdef CHECK_SIZE
                    if (
                        x0 < 0 || y0 < 0 || z0 < 0 ||
                        x1 < 0 || y1 < 0 || z1 < 0 ||
                        x0 >= XSIZE(data) || y0  >= YSIZE(data) || z0 >= ZSIZE(data) ||
                        x1 >= XSIZE(data) || y1  >= YSIZE(data) || z1 >= ZSIZE(data)
                    ) {
                        std::cerr << " x0= " << x0 << " y0= " << y0 << " z0= " << z0 << std::endl;
                        std::cerr << " x1= " << x1 << " y1= " << y1 << " z1= " << z1 << std::endl;
                        data.printShape();
                        REPORT_ERROR("BackProjector::applyPointGroupSymmetry: checksize!!!");
                    }
                    #endif
                    // First interpolate (complex) data
                    Complex d000 = DIRECT_A3D_ELEM(data, z0, y0, x0);
                    Complex d001 = DIRECT_A3D_ELEM(data, z0, y0, x1);
                    Complex d010 = DIRECT_A3D_ELEM(data, z0, y1, x0);
                    Complex d011 = DIRECT_A3D_ELEM(data, z0, y1, x1);
                    Complex d100 = DIRECT_A3D_ELEM(data, z1, y0, x0);
                    Complex d101 = DIRECT_A3D_ELEM(data, z1, y0, x1);
                    Complex d110 = DIRECT_A3D_ELEM(data, z1, y1, x0);
                    Complex d111 = DIRECT_A3D_ELEM(data, z1, y1, x1);

                    Complex dx00 = LIN_INTERP(fx, d000, d001);
                    Complex dx01 = LIN_INTERP(fx, d100, d101);
                    Complex dx10 = LIN_INTERP(fx, d010, d011);
                    Complex dx11 = LIN_INTERP(fx, d110, d111);

                    Complex dxy0 = LIN_INTERP(fy, dx00, dx10);
                    Complex dxy1 = LIN_INTERP(fy, dx01, dx11);

                    // Take complex conjugated for half with negative x
                    A3D_ELEM(sum_data, k, i, j) += is_neg_x ?
                        conj(LIN_INTERP(fz, dxy0, dxy1)) : LIN_INTERP(fz, dxy0, dxy1);

                    // Then interpolate (real) weight
                    RFLOAT dd000 = DIRECT_A3D_ELEM(weight, z0, y0, x0);
                    RFLOAT dd001 = DIRECT_A3D_ELEM(weight, z0, y0, x1);
                    RFLOAT dd010 = DIRECT_A3D_ELEM(weight, z0, y1, x0);
                    RFLOAT dd011 = DIRECT_A3D_ELEM(weight, z0, y1, x1);
                    RFLOAT dd100 = DIRECT_A3D_ELEM(weight, z1, y0, x0);
                    RFLOAT dd101 = DIRECT_A3D_ELEM(weight, z1, y0, x1);
                    RFLOAT dd110 = DIRECT_A3D_ELEM(weight, z1, y1, x0);
                    RFLOAT dd111 = DIRECT_A3D_ELEM(weight, z1, y1, x1);

                    RFLOAT ddx00 = LIN_INTERP(fx, dd000, dd001);
                    RFLOAT ddx01 = LIN_INTERP(fx, dd100, dd101);
                    RFLOAT ddx10 = LIN_INTERP(fx, dd010, dd011);
                    RFLOAT ddx11 = LIN_INTERP(fx, dd110, dd111);

                    RFLOAT ddxy0 = LIN_INTERP(fy, ddx00, ddx10);
                    RFLOAT ddxy1 = LIN_INTERP(fy, ddx01, ddx11);

                    A3D_ELEM(sum_weight, k, i, j) += LIN_INTERP(fz, ddxy0, ddxy1);

                }
            }
        }

        data = sum_data;
        weight = sum_weight;
        // Average
        // The division should only be done if we would search all (C1) directions, not if we restrict the angular search!
        /*
        FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(data)
        {
            DIRECT_MULTIDIM_ELEM(data, n) = DIRECT_MULTIDIM_ELEM(sum_data, n) / (RFLOAT)(SL.SymsNo() + 1);
            DIRECT_MULTIDIM_ELEM(weight, n) = DIRECT_MULTIDIM_ELEM(sum_weight, n) / (RFLOAT)(SL.SymsNo() + 1);
        }
        */
    }

}

void BackProjector::convoluteBlobRealSpace(FourierTransformer &transformer, bool do_mask) {

    MultidimArray<RFLOAT> Mconv;
    int padhdim = pad_size / 2;

    // Set up right dimension of real-space array
    // TODO: resize this according to r_max!!!
    if (ref_dim == 2) {
        Mconv.reshape(pad_size, pad_size);
    } else {
        Mconv.reshape(pad_size, pad_size, pad_size);
    }

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
    FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY3D(Mconv) {
        int kp = k < padhdim ? k : k - pad_size;
        int ip = i < padhdim ? i : i - pad_size;
        int jp = j < padhdim ? j : j - pad_size;
        RFLOAT rval = sqrt((RFLOAT) (kp * kp + ip * ip + jp * jp)) / (ori_size * padding_factor);
        //if (kp==0 && ip==0 && jp > 0)
        //	std::cerr << " jp= " << jp << " rval= " << rval << " tab_ftblob(rval) / normftblob= " << tab_ftblob(rval) / normftblob << " ori_size/2= " << ori_size/2 << std::endl;
        // In the final reconstruction: mask the real-space map beyond its original size to prevent aliasing ghosts
        // Note that rval goes until 1/2 in the oversampled map
        if (do_mask && rval > 1.0 / (2.0 * padding_factor)) {
            DIRECT_A3D_ELEM(Mconv, k, i, j) = 0.0;
        } else {
            DIRECT_A3D_ELEM(Mconv, k, i, j) *= tab_ftblob(rval) / normftblob;
        }
    }

    // forward FFT to go back to Fourier-space
    transformer.FourierTransform();
}

void BackProjector::windowToOridimRealSpace(
    FourierTransformer &transformer, MultidimArray<RFLOAT> &Mout, bool printTimes
) {

    #ifdef TIMING
    Timer OriDimTimer;
    int OriDims[] = {
        OriDimTimer.setNew(" OrD1_getFourier "),
        OriDimTimer.setNew(" OrD2_windowFFT "),
        OriDimTimer.setNew(" OrD3_reshape "),
        OriDimTimer.setNew(" OrD4_setReal "),
        OriDimTimer.setNew(" OrD5_invFFT "),
        OriDimTimer.setNew(" OrD6_centerFFT "),
        OriDimTimer.setNew(" OrD7_window "),
        OriDimTimer.setNew(" OrD8_norm "),
        OriDimTimer.setNew(" OrD9_softMask "),
    };
    #endif

    // Because Fin is a reference, it must be initialised at declaration time,
    // and so we cannot put it in its own scope (as the RCTICTOC macro would require).
    RCTIC(OriDimTimer, OriDims[0]);
    MultidimArray<Complex> &Fin = transformer.getFourierReference();
    RCTOC(OriDimTimer, OriDims[0]);

    int padoridim;  // Size of padded real-space volume
    RCTICTOC(OriDimTimer, OriDims[1], ({
    padoridim = round(padding_factor * ori_size);
    // Enforce divisibility by 2
    padoridim += padoridim % 2;

    // #define DEBUG_WINDOWORIDIMREALSPACE
    #ifdef DEBUG_WINDOWORIDIMREALSPACE
    Image<RFLOAT> tt;
    tt().reshape(ZSIZE(Fin), YSIZE(Fin), XSIZE(Fin));
    FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(Fin) {
        DIRECT_MULTIDIM_ELEM(tt(), n) = abs(DIRECT_MULTIDIM_ELEM(Fin, n));
    }
    tt.write("windoworidim_Fin.spi");
    #endif

    // Resize incoming complex array to the correct size
    windowFourierTransform(Fin, padoridim);
    }))

    RFLOAT normfft;
    RCTICTOC(OriDimTimer, OriDims[2], ({
    if (ref_dim == 2) {
        Mout.reshape(padoridim, padoridim);
        normfft = (
            data_dim == 2 ? 
            padding_factor * padding_factor : 
            padding_factor * padding_factor * ori_size
        );
    } else {
        Mout.reshape(padoridim, padoridim, padoridim);
        normfft = (
            data_dim == 3 ?
            padding_factor * padding_factor * padding_factor :
            padding_factor * padding_factor * padding_factor * ori_size
        );
    }
    Mout.setXmippOrigin();
    }))

    #ifdef DEBUG_WINDOWORIDIMREALSPACE
    tt().reshape(ZSIZE(Fin), YSIZE(Fin), XSIZE(Fin));
    FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(Fin) {
        DIRECT_MULTIDIM_ELEM(tt(), n) = abs(DIRECT_MULTIDIM_ELEM(Fin, n));
    }
    tt.write("windoworidim_Fresized.spi");
    #endif

    // Shift the map back to its origin
    RCTICTOC(OriDimTimer, OriDims[5], ({ CenterFFTbySign(Fin); }))

    RCTICTOC(OriDimTimer, OriDims[3], ({ transformer.setReal(Mout); }))

    // Do the inverse FFT
    RCTICTOC(OriDimTimer, OriDims[4], ({
    #ifdef TIMING
    if (printTimes)
        std::cout << std::endl << "FFTrealDims = (" << transformer.fReal->xdim << " , " << transformer.fReal->ydim << " , " << transformer.fReal->zdim << " ) " << std::endl;
    #endif
    transformer.inverseFourierTransform();
    }))

    // transformer.inverseFourierTransform(Fin, Mout);

    Fin.clear();
    transformer.fReal = NULL;  // Make sure to re-calculate fftw plan
    Mout.setXmippOrigin();

    #ifdef DEBUG_WINDOWORIDIMREALSPACE
    tt() = Mout;
    tt.write("windoworidim_Munwindowed.spi");
    #endif

    // Window in real-space
    RCTICTOC(OriDimTimer, OriDims[6], ({
    if (ref_dim == 2) {
        Mout.window(
            Xmipp::init(ori_size), Xmipp::init(ori_size),
            Xmipp::last(ori_size), Xmipp::last(ori_size)
        );
    } else {
        Mout.window(
            Xmipp::init(ori_size), Xmipp::init(ori_size), Xmipp::init(ori_size),
            Xmipp::last(ori_size), Xmipp::last(ori_size), Xmipp::last(ori_size)
        );
    }
    Mout.setXmippOrigin();
    }))

    // Normalisation factor of FFTW
    // The Fourier Transforms are all "normalised" for 2D transforms of size = ori_size x ori_size
    RCTICTOC(OriDimTimer, OriDims[7], ({ Mout /= normfft; }))

    #ifdef DEBUG_WINDOWORIDIMREALSPACE
    tt() = Mout;
    tt.write("windoworidim_Mwindowed.spi");
    #endif

    // Mask out corners to prevent aliasing artefacts
    RCTICTOC(OriDimTimer, OriDims[8], ({ softMaskOutsideMap(Mout); }))

    #ifdef DEBUG_WINDOWORIDIMREALSPACE
    tt() = Mout;
    tt.write("windoworidim_Mwindowed_masked.spi");
    FourierTransformer ttf;
    ttf.FourierTransform(Mout, Fin);
    tt().resize(ZSIZE(Fin), YSIZE(Fin), XSIZE(Fin));
    FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(Fin) {
        DIRECT_MULTIDIM_ELEM(tt(), n) = abs(DIRECT_MULTIDIM_ELEM(Fin, n));
    }
    tt.write("windoworidim_Fnew.spi");
    #endif

    #ifdef TIMING
    if (printTimes) OriDimTimer.printTimes(true);
    #endif
}
