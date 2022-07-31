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

#include "src/transformations.h"

/* Rotation 2D ------------------------------------------------------------- */
Matrix2D<RFLOAT> rotation2DMatrix(RFLOAT ang, bool homogeneous) {

    int n = homogeneous ? 3 : 2;
    Matrix2D<RFLOAT> result(n, n);

    if (homogeneous) {
        result.at(0, 2) = 0;
        result.at(1, 2) = 0;
        result.at(2, 0) = 0;
        result.at(2, 1) = 0;
        result.at(2, 2) = 1;
    }

    ang = radians(ang);  // Why don't we just get passed angles in radians?
    RFLOAT cosang = cos(ang);
    RFLOAT sinang = sin(ang);

    result.at(0, 0) =  cosang;
    result.at(0, 1) = -sinang;
    result.at(1, 0) =  sinang;
    result.at(1, 1) =  cosang;
}

// Decompose w into an integral part and a fractional part
static int subtract_integral(RFLOAT &w) {
    int i = w;
    w -= i;
    return i;
};

template <typename T>
T interpolate_sub(
    const MultidimArray<T> &V1,
    RFLOAT xp, RFLOAT yp,
    int cen_xp, int cen_yp,
    bool do_wrap
) {
    // Linear interpolation

    // Calculate the integer position in input image.
    // Be careful that it is not the nearest but the one at the top left corner
    // of the interpolation square.
    // ie, (0.7, 0.7) would give (0, 0).
    // Also calculate weights for point (m1 + 1, n1 + 1).
    RFLOAT wx = xp + cen_xp;
    const int m1 = subtract_integral(wx);
          int m2 = m1 + 1;

    RFLOAT wy = yp + cen_yp;
    const int n1 = subtract_integral(wy);
          int n2 = n1 + 1;

    // In case m2 and n2 are out of bounds
    if (do_wrap) {
        if (m2 >= Xsize(V1)) { m2 = 0; }
        if (n2 >= Ysize(V1)) { n2 = 0; }
    }

    #ifdef DEBUG_APPLYGEO
    std::cout << "   From "
        << "(" << m1 << "," << n1 << ") and "
        << "(" << m2 << "," << n2 << ")\n"
        << "   wx= " << wx << " wy= " << wy << std::endl;
    #endif

    // Perform interpolation
    // If wx == 0, then the rightmost point is useless for this interpolation,
    // and might not even be defined if m1 == xdim - 1.
    T tmp = ((1 - wy) * (1 - wx) * direct::elem(V1, m1, n1));

    if (m2 < V1.xdim)
        tmp += (T) ((1 - wy) * wx * direct::elem(V1, m2, n1));

    if (n2 < V1.ydim) {
        tmp += (T) (wy * (1 - wx) * direct::elem(V1, m1, n2));

        if (m2 < V1.xdim)
            tmp += (T) (wy * wx * direct::elem(V1, m2, n2));
    }

    #ifdef DEBUG_APPYGEO
    std::cout << "   val= " << tmp << std::endl;
    #endif
    return tmp;
}

template <typename T>
T interpolate_sub(
    const MultidimArray<T> &V1,
    RFLOAT xp, RFLOAT yp, RFLOAT zp,
    int cen_xp, int cen_yp, int cen_zp,
    bool show_debug
) {
    // Linear interpolation

    // Calculate the integer position in input volume,
    // be careful that it is not the nearest but the one at the
    // top left corner of the interpolation square.
    // ie (0.7, 0.7) would give (0, 0)
    // Also calculate weights for point (m1 + 1, n1 + 1)
    RFLOAT wx = xp + cen_xp;
    const int m1 = subtract_integral(wx);
    const int m2 = m1 + 1;

    RFLOAT wy = yp + cen_yp;
    const int n1 = subtract_integral(wy);
    const int n2 = n1 + 1;

    RFLOAT wz = zp + cen_zp;
    const int o1 = subtract_integral(wz);
    const int o2 = o1 + 1;

    #ifdef DEBUG
    if (show_debug) {
        std::cout << "After wrapping(xp,yp,zp)= "
        << "(" << xp << "," << yp << "," << zp << ")\n";
        std::cout << "(m1,n1,o1)-->(m2,n2,o2)="
        << "(" << m1 << "," << n1 << "," << o1 << ") "
        << "(" << m2 << "," << n2 << "," << o2 << ")\n";
        std::cout << "(wx,wy,wz)="
        << "(" << wx << "," << wy << "," << wz << ")\n";
    }
    #endif

    // Perform interpolation
    // If wx == 0, then the rightmost point is useless for this interpolation,
    // and might not even be defined if m1 == xdim - 1.
    T tmp = ((1 - wz) * (1 - wy) * (1 - wx) * direct::elem(V1, m1, n1, o1));

    if (m2 < V1.xdim)
        tmp += (T) ((1 - wz) * (1 - wy) * wx * direct::elem(V1, m2, n1, o1));

    if (n2 < V1.ydim) {
        tmp += (T) ((1 - wz) * wy * (1 - wx) * direct::elem(V1, m1, n2, o1));
        if (m2 < V1.xdim)
            tmp += (T) ((1 - wz) * wy * wx * direct::elem(V1, m2, n2, o1));
    }

    if (o2 < V1.zdim) {
        tmp += (T) (wz * (1 - wy) * (1 - wx) * direct::elem(V1, m1, n1, o2));
        if (m2 < V1.xdim)
        tmp += (T) (wz * (1 - wy) * wx * direct::elem(V1, m2, n1, o2));
        if (n2 < V1.ydim) {
            tmp += (T) (wz * wy * (1 - wx) * direct::elem(V1, m1, n2, o2));
            if (m2 < V1.xdim)
                tmp += (T) (wz * wy * wx * direct::elem(V1, m2, n2, o2));
        }
    }

    #ifdef DEBUG
    if (show_debug)
        std::cout <<
        "tmp1=" << direct::elem(V1, m1, n1, o1) << " "
        << (T) ((1 - wz) * (1 - wy) * (1 - wx) * direct::elem(V1, m1, n1, o1))
        << "\n" <<
        "tmp2=" << direct::elem(V1, m2, n1, o1) << " "
        << (T) ((1 - wz) * (1 - wy) * wx * direct::elem(V1, m2, n1, o1))
        << "\n" <<
        "tmp3=" << direct::elem(V1, m1, n2, o1) << " "
        << (T) ((1 - wz) * wy * (1 - wx) * direct::elem(V1, m1, n2, o1))
        << "\n" <<
        "tmp4=" << direct::elem(V1, m2, n2, o1) << " "
        << (T) ((1 - wz) * wy * wx * direct::elem(V1, m1, n1, o2))
        << "\n" <<
        "tmp6=" << direct::elem(V1, m2, n1, o2) << " "
        << (T) (wz * (1 - wy) * wx * direct::elem(V1, m2, n1, o2))
        << "\n" <<
        "tmp7=" << direct::elem(V1, m1, n2, o2) << " "
        << (T) (wz * wy * (1 - wx) * direct::elem(V1, m1, n2, o2))
        << "\n" <<
        "tmp8=" << direct::elem(V1, m2, n2, o2) << " "
        << (T) (wz * wy * wx * direct::elem(V1, m2, n2, o2))
        << "\n" <<
        "tmp= " << tmp << std::endl;
    #endif
    return tmp;
}

// Manual template instantiation

template RFLOAT interpolate_sub(
    const MultidimArray<RFLOAT>&, RFLOAT, RFLOAT, int, int, bool
);

template RFLOAT interpolate_sub(
    const MultidimArray<RFLOAT>&, RFLOAT, RFLOAT, RFLOAT, int, int, int, bool
);

/* Translation 2D ---------------------------------------------------------- */
Matrix2D<RFLOAT> translation2DMatrix(const Matrix1D<RFLOAT> &v) {
    // if (v.size() != 2)
    //    REPORT_ERROR("Translation2D_matrix: vector is not in R2");
    Matrix2D<RFLOAT> result;
    result.initIdentity(3);
    result.at(0, 2) = XX(v);
    result.at(1, 2) = YY(v);
    return result;
}

/* Rotation 3D around the system axes -------------------------------------- */
Matrix2D<RFLOAT> rotation3DMatrix(
    RFLOAT ang, char axis, bool homogeneous
) {
    const int n = homogeneous ? 4 : 3;
    auto result = Matrix2D<RFLOAT>::zeros(n, n);
    if (homogeneous) { result.at(3, 3) = 1; }

    RFLOAT cosa = cos(radians(ang));
    RFLOAT sina = sin(radians(ang));

    switch (axis) {

        case 'Z':
        result.at(0, 0) =  cosa;
        result.at(0, 1) = -sina;
        result.at(1, 0) =  sina;
        result.at(1, 1) =  cosa;
        result.at(2, 2) = 1;
        break;

        case 'Y':
        result.at(0, 0) =  cosa;
        result.at(0, 2) = -sina;
        result.at(2, 0) =  sina;
        result.at(2, 2) =  cosa;
        result.at(1, 1) = 1;
        break;

        case 'X':
        result.at(1, 1) =  cosa;
        result.at(1, 2) = -sina;
        result.at(2, 1) =  sina;
        result.at(2, 2) =  cosa;
        result.at(0, 0) = 1;
        break;

        default:
        REPORT_ERROR("rotation3DMatrix: Unknown axis");

    }
    return result;
}

/* Align a vector with Z axis */
void alignWithZ(
    const Matrix1D<RFLOAT> &axis, Matrix2D<RFLOAT> &result, bool homogeneous
) {
    if (axis.size() != 3)
        REPORT_ERROR("alignWithZ: Axis is not in R3");
    if (homogeneous) {
        result.initZeros(4, 4);
        result.at(3, 3) = 1;
    } else {
        result.initZeros(3, 3);
    }
    Matrix1D<RFLOAT> Axis(axis);
    Axis.normalise();

    // Compute length of the projection on YZ plane
    RFLOAT proj_mod = sqrt(YY(Axis) * YY(Axis) + ZZ(Axis) * ZZ(Axis));
    if (proj_mod > Xmipp::epsilon) {   
        // proj_mod != 0
        // Build Matrix result, which makes the turning axis coincident with Z
        result.at(0, 0) = proj_mod;
        result.at(0, 1) = -XX(Axis) * YY(Axis) / proj_mod;
        result.at(0, 2) = -XX(Axis) * ZZ(Axis) / proj_mod;
        result.at(1, 0) = 0;
        result.at(1, 1) = ZZ(Axis) / proj_mod;
        result.at(1, 2) = -YY(Axis) / proj_mod;
        result.at(2, 0) = XX(Axis);
        result.at(2, 1) = YY(Axis);
        result.at(2, 2) = ZZ(Axis);
    } else {
        // I know that the Axis is the X axis, EITHER POSITIVE OR NEGATIVE!!
        result.at(0, 0) = 0;
        result.at(0, 1) = 0;
        result.at(0, 2) = -sgn_nozero(XX(Axis));
        result.at(1, 0) = 0;
        result.at(1, 1) = 1;
        result.at(1, 2) = 0;
        result.at(2, 0) = sgn_nozero(XX(Axis));
        result.at(2, 1) = 0;
        result.at(2, 2) = 0;
    }
}

/* Rotation 3D around any axis -------------------------------------------- */
Matrix2D<RFLOAT> rotation3DMatrix(
    RFLOAT ang, const Matrix1D<RFLOAT> &axis,
    bool homogeneous
) {
    // Compute a matrix which makes the turning axis coincident with Z
    // And turn around this axis
    Matrix2D<RFLOAT> A;
    alignWithZ(axis, A, homogeneous);
    const Matrix2D<RFLOAT> R = rotation3DMatrix(ang, 'Z', homogeneous);
    return A.transpose() * R * A;
}

/* Translation 3D ---------------------------------------------------------- */
Matrix2D<RFLOAT> translation3DMatrix(const Matrix1D<RFLOAT> &v) {
    if (v.size() != 3)
        REPORT_ERROR("Translation3D_matrix: vector is not in R3");
    Matrix2D<RFLOAT> result;
    result.initIdentity(4);
    result.at(0, 3) = XX(v);
    result.at(1, 3) = YY(v);
    result.at(2, 3) = ZZ(v);
    return result;
}

/* Scale 3D ---------------------------------------------------------------- */
Matrix2D<RFLOAT> scale3DMatrix(
    const Matrix1D<RFLOAT> &sc, bool homogeneous
) {
    if (sc.size() != 3)
        REPORT_ERROR("Scale3D_matrix: vector is not in R3");

    const int n = homogeneous ? 4 : 3;
    auto result = Matrix2D<RFLOAT>::zeros(n, n);
    if (homogeneous) { result.at(3, 3) = 1; }

    result.at(0, 0) = XX(sc);
    result.at(1, 1) = YY(sc);
    result.at(2, 2) = ZZ(sc);
    return result;
}
