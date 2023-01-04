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
 *              Sjors H.W. Scheres
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

#ifndef TRANSFORMATIONS_H
#define TRANSFORMATIONS_H

#include "src/multidim_array.h"
#include "src/euler.h"

const bool IS_INV = true;
const bool IS_NOT_INV = false;
const bool DONT_WRAP = false;
const bool WRAP = true;

#ifndef DBL_EPSILON
#define DBL_EPSILON 1e-50
#endif

/// @defgroup GeometricalTransformations Geometrical transformations
/// @ingroup DataLibrary
//@{
/** Creates a rotational matrix (3x3) for images
 * @ingroup GeometricalTransformations
 *
 * The rotation angle is in degrees.
 * m must have been already resized to 3x3
 *
 * @code
 * m = rotation2DMatrix(60);
 * @endcode
 */
Matrix<RFLOAT> rotation2DMatrix(RFLOAT ang, bool homogeneous=true);

/** Creates a translational matrix (3x3) for images
 * @ingroup GeometricalTransformations
 *
 * The shift is given as a R2 vector (shift_X, shift_Y). An exception is thrown
 * if the displacement is not a R2 vector.
 *
 * @code
 * // Displacement of 1 pixel to the right
 * m = translation2DMatrix(vectorR2(1, 0));
 * @endcode
 */
Matrix<RFLOAT> translation2DMatrix(const Vector<RFLOAT> &v);

/** Creates a rotational matrix (4x4) for volumes around system axis
 * @ingroup GeometricalTransformations
 *
 * The rotation angle is in degrees, and the rotational axis is either 'X', 'Y'
 * or 'Z'. An exception is thrown if the axis given is not one of these.
 *
 * The returned matrices are respectively alpha degrees around Z
 *
 * @code
 * [ cos(A) -sin(A)     0   ]
 * [ sin(A)  cos(A)     0   ]
 * [   0       0        1   ]
 * @endcode
 *
 * alpha degrees around Y
 * @code
 * [ cos(A)    0    -sin(A) ]
 * [   0       1       0    ]
 * [ sin(A)    0     cos(A) ]
 * @endcode
 *
 * alpha degrees around X
 * @code
 * [   1       0       0    ]
 * [   0     cos(A) -sin(A) ]
 * [   0     sin(A)  cos(A) ]
 * @endcode
 *
 * @code
 * m = rotation3DMatrix(60, 'X');
 * @endcode
 */
Matrix<RFLOAT> rotation3DMatrix(RFLOAT ang, char axis, bool homogeneous=true);

/** Creates a rotational matrix (4x4) for volumes around any axis
 * @ingroup GeometricalTransformations
 *
 * The rotation angle is in degrees, and the rotational axis is given as a R3
 * vector. An exception is thrown if the axis is not a R3 vector. The axis needs
 * not to be unitary.
 *
 * @code
 * m = rotation3DMatrix(60, vectorR3(1, 1, 1));
 * @endcode
 */
Matrix<RFLOAT> rotation3DMatrix(RFLOAT ang, const Vector<RFLOAT> &axis, bool homogeneous=true);

/** Matrix which transforms the given axis into Z
 * @ingroup GeometricalTransformations
 *
 * A geometrical transformation matrix (4x4) is returned such that the given
 * axis is rotated until it is aligned with the Z axis. This is very useful in
 * order to produce rotational matrices, for instance, around any axis.
 *
 * @code
 * Matrix<RFLOAT> A = alignWithZ(axis);
 * return A.transpose() * rotation3DMatrix(ang, 'Z') * A;
 * @endcode
 *
 * The returned matrix is such that A*axis=Z, where Z and axis are column
 * vectors.
 */
void alignWithZ(const Vector<RFLOAT> &axis, Matrix<RFLOAT> &m, bool homogeneous=true);

/** Creates a translational matrix (4x4) for volumes
 * @ingroup GeometricalTransformations
 *
 * The shift is given as a R3 vector (shift_X, shift_Y, shift_Z). An exception
 * is thrown if the displacement is not a R3 vector.
 *
 * @code
 * // Displacement of 2 pixels down
 * m = translation3DMatrix(vectorR3(0, 0, 2));
 * @endcode
 */
Matrix<RFLOAT> translation3DMatrix(const Vector<RFLOAT> &v);

/** Creates a scaling matrix (4x4) for volumes
 * @ingroup GeometricalTransformations
 *
 * The scaling factors for the different axis must be given as a vector. So
 * that, XX(sc)=scale for X axis, YY(sc)=...
 */
Matrix<RFLOAT> scale3DMatrix(const Vector<RFLOAT> &sc, bool homogeneous=true);

/** Applies a geometrical transformation.
 * @ingroup GeometricalTransformations
 *
 * Apply a geometrical transformation defined by the matrix A (RFLOAT (4x4)!!
 * ie, in homogeneous R3 coordinates) to the volume V.
 * The result is stored in V2 (it cannot be the same as the input volume).
 * An exception is thrown if the transformation matrix is not 4x4.
 *
 * Structure of the transformation matrix: It should have the following
 * components
 *
 * r11 r12 r13 x
 * r21 r22 r23 y
 * r31 r32 r33 z
 * 0   0   0   1
 *
 * where (x,y,z) is the translation desired, and Rij are the components of
 * the rotation matrix R. If you want to apply a scaling factor to the
 * transformation, then multiply r11, r22 and r33 by it.
 *
 * The result volume (with ndim=1) is resized to the same
 * dimensions as V if V2 is empty (0x0) at the beginning, if it
 * is not, ie, if V2 has got some size then only those values in
 * the volume are filled, this is very useful for resizing the
 * volume, then you manually resize the output volume to the
 * desired size and then call this routine.
 *
 * The relationship between the output coordinates and the input ones are
 *
 * @code
 * out = A * in
 * (x, y, z) = A * (x', y', z')
 * @endcode
 *
 * This function works independently from the logical indexing of each
 * matrix, it sets the logical center and the physical center of the image
 * and work with these 2 coordinate spaces. At the end the original logical
 * indexing of each matrix is kept.
 *
 * The procedure followed goes from coordinates in the output volume
 * to the ones in the input one, so the inverse of the A matrix is
 * needed. There is a flag telling if the given matrix is already
 * the inverse one or the normal one. If it is the normal one internally
 * the matrix is inversed. If you are to do many "rotations" then
 * some time is spent in inverting the matrix. Normally the matrix is the
 * normal one.
 *
 * There is something else to tell about the geometrical tranformation.
 * The value of the voxel in the output volume is computed via
 * bilinear interpolation in the input volume. If any of the voxels
 * participating in the interpolation falls outside the input volume,
 * then automatically the corresponding output voxel is set to 0, unless
 * that the do_wrap flag has been set to 1. In this case if the voxel
 * falls out by the right hand then it is "wrapped" and the corresponding
 * voxel in the left hand is used. The same is appliable to top-bottom.
 * Usually wrap mode is off. Wrap mode is interesting for translations
 * but not for rotations, for example.
 *
 * The inverse mode and wrapping mode should be taken by default by the
 * routine, g++ seems to have problems with template functions outside
 * a class with default parameters. So, I'm sorry, you will have to
 * put them always. The usual combination is
 *
 * applyGeometry(..., IS_NOT_INV, DONT_WRAP).
 *
 * Although you can also use the constants IS_INV, or WRAP.
 *
 * @code
 * Matrix<RFLOAT> A = Matrix<RFLOAT>::identity(4);
 * A = applyGeometry(A, V);
 * @endcode
 */
template<typename T>
MultidimArray<T> applyGeometry(
    const MultidimArray<T> &V,
    const Matrix<RFLOAT> A,
    bool inv, bool do_wrap,
    T outside = 0
) {

    if (V.getDim() == 2 && (A.ncols() != 3 || A.nrows() != 3))
        REPORT_ERROR("ApplyGeometry: 2D transformation matrix is not 3×3");

    if (V.getDim() == 3 && (A.ncols() != 4 || A.nrows() != 4))
        REPORT_ERROR("ApplyGeometry: 3D transformation matrix is not 4×4");

    if (A.isIdentity())
        return V;

    if (Xsize(V) == 0)
        return {};

    MultidimArray<RFLOAT> result;
    result.resize(V);

    Matrix<RFLOAT> Ainv;
    const Matrix<RFLOAT> *Aptr = &A;
    if (!inv) {
        Ainv = A.inv();
        Aptr = &Ainv;
    }
    const Matrix<RFLOAT> &Aref = *Aptr;

    // For scalings the output matrix is resized outside to the final
    // size instead of being resized inside the routine with the
    // same size as the input matrix
    if (Xsize(result) == 0)
        result.resize(V);

    if (V.getDim() == 2) {
        // 2D transformation

        // Find center and limits of image
        int cen_x  = Xsize(result) / 2;
        int cen_y  = Ysize(result) / 2;
        int cen_xp = Xsize(V) / 2;
        int cen_yp = Ysize(V) / 2;
        RFLOAT minxp = -cen_xp;
        RFLOAT minyp = -cen_yp;
        RFLOAT maxxp = Xsize(V) - cen_xp - 1;
        RFLOAT maxyp = Ysize(V) - cen_yp - 1;

        // Now we go from the output image to the input image.
        // ie, for any pixel in the output image
        // we calculate which are the corresponding ones in
        // the original image, make an interpolation with them and put this value
        // at the output pixel

        #ifdef DEBUG_APPLYGEO
        std::cout << "A\n" << Aref << std::endl
        << "(cen_x ,cen_y )=(" << cen_x  << "," << cen_y  << ")\n"
        << "(cen_xp,cen_yp)=(" << cen_xp << "," << cen_yp << ")\n"
        << "(min_xp,min_yp)=(" << minxp  << "," << minyp  << ")\n"
        << "(max_xp,max_yp)=(" << maxxp  << "," << maxyp  << ")\n";
        #endif

        for (int j = 0; j < Ysize(result); j++) {
            // Calculate position of the beginning of the row in the output image
            RFLOAT x =   - cen_x;
            RFLOAT y = j - cen_y;

            // Calculate this position in the input image according to the
            // geometrical transformation
            // they are related by
            // coords_output(=x,y) = A * coords_input (=xp,yp)
            RFLOAT xp = x * Aref(0, 0) + y * Aref(0, 1) + Aref(0, 2);
            RFLOAT yp = x * Aref(1, 0) + y * Aref(1, 1) + Aref(1, 2);

            for (int i = 0; i < Xsize(result); i++) {

                #ifdef DEBUG_APPLYGEO

                std::cout << "Computing "       "(" << i << "," << j << ")\n"
                          << "   (y, x) ="      "(" << y << "," << x << ")\n"
                << "   before wrapping (y',x')=""(" << yp << "," << xp << ")"
                << std::endl;
                #endif
                // If the point is outside the image, apply a periodic extension of the image
                // What exits by one side enters by the other
                bool interp = do_wrap ||
                    !Xmipp::lt(xp, minxp) && !Xmipp::gt(xp, maxxp) &&
                    !Xmipp::lt(yp, minyp) && !Xmipp::gt(yp, maxyp);

                if (do_wrap) {

                    if (Xmipp::lt(xp, minxp) || Xmipp::gt(xp, maxxp))
                        xp = wrap(xp, minxp - 0.5, maxxp + 0.5);

                    if (Xmipp::lt(yp, minyp) || Xmipp::gt(yp, maxyp))
                        yp = wrap(yp, minyp - 0.5, maxyp + 0.5);

                }

                #ifdef DEBUG_APPLYGEO
                std::cout << "   after wrapping (y',x')=(" << yp << "," << xp << ")\n"
                          << "   Interp = " << interp << std::endl;
                // The following line sounds dangerous...
                // x++;
                #endif
                direct::elem(result, i, j) = interp ? interpolate_sub(
                    V, xp, yp, cen_xp, cen_yp, do_wrap
                ) : outside;

                // Compute new point inside input image
                xp += Aref(0, 0);
                yp += Aref(1, 0);
            }
        }
    } else {
        // 3D transformation

        // Find center of MultidimArray
        int cen_z  = result.zdim / 2;
        int cen_y  = result.ydim / 2;
        int cen_x  = result.xdim / 2;
        int cen_zp = V.zdim / 2;
        int cen_yp = V.ydim / 2;
        int cen_xp = V.xdim / 2;
        RFLOAT minxp = -cen_xp;
        RFLOAT minyp = -cen_yp;
        RFLOAT minzp = -cen_zp;
        RFLOAT maxxp = V.xdim - cen_xp - 1;
        RFLOAT maxyp = V.ydim - cen_yp - 1;
        RFLOAT maxzp = V.zdim - cen_zp - 1;

        #ifdef DEBUG
        std::cout
        << "Geometry 2 center=("
        << cen_z  << "," << cen_y  << "," << cen_x  << ")\n"
        << "Geometry 1 center=("
        << cen_zp << "," << cen_yp << "," << cen_xp << ")\n"
        << "           min=("
        << minzp  << "," << minyp  << "," << minxp  << ")\n"
        << "           max=("
        << maxzp  << "," << maxyp  << "," << maxxp  << ")\n"
        ;
        #endif

        // Now we go from the output MultidimArray to the input MultidimArray, ie, for any
        // voxel in the output MultidimArray we calculate which are the corresponding
        // ones in the original MultidimArray, make an interpolation with them and put
        // this value at the output voxel

        // result is not initialised to 0 because all its pixels are rewritten
        for (int k = 0; k < result.zdim; k++)
        for (int j = 0; j < result.ydim; j++) {
            // Calculate position of the beginning of the row in the output MultidimArray
            RFLOAT x =   - cen_x;
            RFLOAT y = j - cen_y;
            RFLOAT z = k - cen_z;

            // Calculate this position in the input image according to the
            // geometrical transformation they are related by
            // coords_output(=x,y) = A * coords_input (=xp,yp)
            RFLOAT xp = x * Aref(0, 0) + y * Aref(0, 1) + z * Aref(0, 2) + Aref(0, 3);
            RFLOAT yp = x * Aref(1, 0) + y * Aref(1, 1) + z * Aref(1, 2) + Aref(1, 3);
            RFLOAT zp = x * Aref(2, 0) + y * Aref(2, 1) + z * Aref(2, 2) + Aref(2, 3);

            for (int i = 0; i < result.xdim; i++) {

                bool show_debug = false;
                #ifdef DEBUG
                show_debug =
                    i == 0 && j == 0 && k == 0 ||
                    i == result.xdim - 1 && j == result.ydim - 1 && k == result.zdim - 1;

                if (show_debug)
                    std::cout << "(x,y,z)-->(xp,yp,zp)= "
                    << "(" << x  << "," << y  << "," << z  << ") "
                    << "(" << xp << "," << yp << "," << zp << ")\n";
                #endif

                // If the point is outside the volume, apply a periodic
                // extension of the volume, what exits by one side enters by
                // the other
                bool interp = do_wrap ||
                    !Xmipp::lt(xp, minxp) && !Xmipp::gt(xp, maxxp) &&
                    !Xmipp::lt(yp, minyp) && !Xmipp::gt(yp, maxyp) &&
                    !Xmipp::lt(zp, minzp) && !Xmipp::gt(zp, maxzp);

                if (do_wrap) {

                    if (
                        Xmipp::lt(xp, minxp) || Xmipp::gt(xp, maxxp)
                    ) { xp = wrap(xp, minxp - 0.5, maxxp + 0.5); }

                    if (
                        Xmipp::lt(yp, minyp) || Xmipp::gt(yp, maxyp)
                    ) { yp = wrap(yp, minyp - 0.5, maxyp + 0.5); }

                    if (
                        Xmipp::lt(zp, minzp) || Xmipp::gt(zp, maxzp)
                    ) { zp = wrap(zp, minzp - 0.5, maxzp + 0.5); }

                }

                direct::elem(result, i, j, k) = interp ? interpolate_sub(
                    V, xp, yp, zp, cen_xp, cen_yp, cen_zp, show_debug
                ) : outside;

                // Compute new point inside input image
                xp += Aref(0, 0);
                yp += Aref(1, 0);
                zp += Aref(2, 0);
            }
        }
    }
    return result;
}

template <typename T>
T interpolate_sub(
    const MultidimArray<T>& V,
    RFLOAT xp, RFLOAT yp,
    int cen_xp, int cen_yp,
    bool do_wrap
);

template <typename T>
T interpolate_sub(
    const MultidimArray<T> &V,
    RFLOAT xp, RFLOAT yp, RFLOAT zp,
    int cen_xp, int cen_yp, int cen_zp,
    bool show_debug
);

/** Rotate an array around a given system axis.
 * @ingroup GeometricalTransformations
 *
 * The rotation angle is in degrees, and the rotational axis is either
 * 'X', 'Y' or 'Z' for 3D arrays and only 'Z' for 2D arrays. An
 * exception is thrown if the axis given is not one of these.
 *
 * @code
 * Vout = rotate(Vin, 60);
 * @endcode
 */

template<typename T>
MultidimArray<T> rotate(
    const MultidimArray<T> &V,
    RFLOAT ang, char axis = 'Z',
    bool do_wrap = DONT_WRAP, T outside = 0
) {
    Matrix<RFLOAT> m;
    switch (V.getDim()) {
        case 2: m = rotation2DMatrix(ang); break;  // axis is irrelevant
        case 3: m = rotation3DMatrix(ang, axis); break;
        default: REPORT_ERROR("rotate ERROR: rotate only valid for 2D or 3D arrays");
    }
    return applyGeometry(V, m, IS_NOT_INV, do_wrap, outside);
}

/** Translate an array.
 * @ingroup GeometricalTransformations
 *
 * The shift is given as a R2 or R3 vector (shift_X, shift_Y, shift_Z) for 2D and 3D arrays, respectively.
 * An exception is thrown if the displacement is not a R3 vector.
 *
 * @code
 * // Displacement of 2 pixels down
 * V2 = translate(V, vectorR3(0, 0, 2));
 * @endcode
 */
template<typename T>
MultidimArray<T> translate(
    const MultidimArray<T> &V,
    const Vector<RFLOAT> &v,
    bool do_wrap = WRAP, T outside = 0
) {
    Matrix<RFLOAT> m;
    switch (V.getDim()) {
        case 2: m = translation2DMatrix(v); break;
        case 3: m = translation3DMatrix(v); break;
        default: REPORT_ERROR("translate ERROR: translate only valid for 2D or 3D arrays");
    }
    return applyGeometry(V, m, IS_NOT_INV, do_wrap, outside);
}

/** Translate center of mass to center
 * @ingroup GeometricalTransformations
 *
 * If the input has very high values, it is better to rescale it to be
 * between 0 and 1.
 */
template<typename T>
MultidimArray<T> translateCenterOfMassToCenter(
    const MultidimArray<T> &V,
    bool do_wrap = WRAP,
    bool verb = false
) {
    MultidimArray<T> V2 = V;
    V2.setXmippOrigin();
    Vector<RFLOAT> center;
    V2.centerOfMass(center);
    if (verb) {
    	std::cout << " Center of mass: x= " << XX(center) << " y= " << YY(center) << " z= " << ZZ(center) << std::endl;
    }
    return translate(V, -center, do_wrap, 0.0);
}

/** Scales to a new size.
 * @ingroup GeometricalTransformations
 *
 * The volume is scaled (resampled) to fill a new size. It is not the
 * same as "window" in this same class. The size can be larger or smaller
 * than the actual one.
 *
 * @code
 * Vout = scaleToSize(Vin, 128, 128, 128);
 * @endcode
 */
template<typename T>
MultidimArray<T> scaleToSize(
    const MultidimArray<T> &V,
    int Xdim, int Ydim, int Zdim = 1
) {

    Matrix<RFLOAT> tmp;
    switch (V.getDim()) {

        case 2:
        tmp.resize(3, 3);
        tmp.setIdentity();
        tmp(0, 0) = (RFLOAT) Xdim / (RFLOAT) Xsize(V);
        tmp(1, 1) = (RFLOAT) Ydim / (RFLOAT) Ysize(V);
        break;

        case 3:
        tmp.resize(4, 4);
        tmp.setIdentity();
        tmp(0, 0) = (RFLOAT) Xdim / (RFLOAT) Xsize(V);
        tmp(1, 1) = (RFLOAT) Ydim / (RFLOAT) Ysize(V);
        tmp(2, 2) = (RFLOAT) Zdim / (RFLOAT) Zsize(V);
        break;

        default:
        REPORT_ERROR("scaleToSize ERROR: scaleToSize only valid for 2D or 3D arrays");

    }

    return applyGeometry(V, tmp, IS_NOT_INV, WRAP, (T) 0);
}

/** Does a radial average of a 2D/3D image, around the voxel where is the origin.
 * @ingroup GeometricalTransformations
 *
 * A vector radial_mean is returned where:
 * - the first element is the mean of the voxels whose
 *   distance to the origin is (0-1),
 * - the second element is the mean of the voxels
 *   whose distance to the origin is (1-2)
 * - and so on.
 *
 * A second vector radial_count is returned containing the number of voxels
 * over which each radial average was calculated.
 *
 * Sjors nov2003: if rounding=true, element=round(distance);
 * - so the first element is the mean of the voxels whose distance to the
 *   origin is (0.5-1.5),
 * - the second element is the mean of the voxels whose distance to the origin
 *   is (1.5-2.5)
 * - and so on.
 */
template<typename T>
void radialAverage(
    const MultidimArray<T> &m,
    Vector<int> &center_of_rot,
    MultidimArray<T> &radial_mean,
    MultidimArray<int> &radial_count,
    bool rounding = false
) {

    // If center_of_rot was written for 2D image
    if (center_of_rot.size() < 3) { center_of_rot.resize(3); }

    // First determine the maximum distance that one should expect,
    // to set the dimension of the radial average vector.
    std::vector<int> distances;
    distances.reserve(8);

    RFLOAT x = Xinit(m) - XX(center_of_rot),
           y = Yinit(m) - YY(center_of_rot),
           z = Zinit(m) - ZZ(center_of_rot);
    distances.push_back(floor(hypot(x, y, z)));

    x = Xlast(m) - XX(center_of_rot);
    distances.push_back(floor(hypot(x, y, z)));

    y = Ylast(m) - YY(center_of_rot);
    distances.push_back(floor(hypot(x, y, z)));

    x = Xinit(m) - XX(center_of_rot);
    distances.push_back(floor(hypot(x, y, z)));

    z = Zlast(m) - ZZ(center_of_rot);
    distances.push_back(floor(hypot(x, y, z)));

    x = Xlast(m) - XX(center_of_rot);
    distances.push_back(floor(hypot(x, y, z)));

    y = Yinit(m) - YY(center_of_rot);
    distances.push_back(floor(hypot(x, y, z)));

    x = Xinit(m) - XX(center_of_rot);
    distances.push_back(floor(hypot(x, y, z)));

    const int dim = ceil(*std::max_element(distances.begin(), distances.end())) + 1 + rounding;

    // Define the vectors
    radial_mean.resize(dim);
    radial_mean.initZeros();
    radial_count.resize(dim);
    radial_count.initZeros();

    Vector<RFLOAT> idx(3);
    // Perform the radial sum and count pixels that contribute to every distance
    FOR_ALL_ELEMENTS_IN_ARRAY3D(m, i, j, k) {
        XX(idx) = i - XX(center_of_rot);
        YY(idx) = j - YY(center_of_rot);
        ZZ(idx) = k - ZZ(center_of_rot);

        // Determine distance to the center
        int distance = rounding ? round(idx.modulus()) : floor(idx.modulus());

        // Sum te value to the pixels with the same distance
        radial_mean.elem(distance) += m.elem(i, j, k);

        // Count the pixel
        radial_count.elem(distance)++;
    }

    // Perform the mean
    for (int i = Xinit(radial_mean); i <= Xlast(radial_mean); i++) {
        radial_mean.elem(i) /= (T) radial_count.elem(i);
    }
}

//@}
#endif
