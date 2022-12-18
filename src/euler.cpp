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

#include <iostream>
#include <math.h>

#include "src/euler.h"
#include "src/funcs.h"

/* Euler angles --> matrix ------------------------------------------------- */
Matrix<RFLOAT> Euler::angles2matrix(
    RFLOAT alpha, RFLOAT beta, RFLOAT gamma,
    bool homogeneous
) {

    Matrix<RFLOAT> A;

    if (homogeneous) {
        A.resize(4, 4);
        std::fill(A.begin(), A.end(), 0);
        A.at(3, 3) = 1;
    } else {
        A.resize(3, 3);
    }

    alpha = radians(alpha);
    beta  = radians(beta);
    gamma = radians(gamma);

    RFLOAT cosa = cos(alpha);
    RFLOAT cosb = cos(beta);
    RFLOAT cosg = cos(gamma);
    RFLOAT sina = sin(alpha);
    RFLOAT sinb = sin(beta);
    RFLOAT sing = sin(gamma);

    // https://en.wikipedia.org/wiki/Euler_angles#Rotation_matrix
    // ZYZ

    A(0, 0) =  cosa * cosb * cosg - sina * sing;
    A(0, 1) =  sina * cosb * cosg + cosa * sing;
    A(0, 2) = -sinb * cosg;
    A(1, 0) = -cosa * cosb * sing - sina * cosg;
    A(1, 1) = -sina * cosb * sing + cosa * cosg;
    A(1, 2) =  sinb * sing;
    A(2, 0) =  cosa * sinb;
    A(2, 1) =  sina * sinb;
    A(2, 2) =  cosb;
    return A;
}

/* Euler direction --------------------------------------------------------- */
Vector<RFLOAT> Euler::angles2direction(RFLOAT alpha, RFLOAT beta) {

    alpha = radians(alpha);
    beta  = radians(beta);

    RFLOAT cosa = cos(alpha);
    RFLOAT cosb = cos(beta);
    RFLOAT sina = sin(alpha);
    RFLOAT sinb = sin(beta);

    return vectorR3(cosa * sinb, sina * sinb, cosb);
}

/* Euler direction2angles ------------------------------- */
//gamma is useless but I keep it for simmetry
//with Euler::direction
void Euler::direction2angles(
    Vector<RFLOAT> &v0,
    RFLOAT &alpha, RFLOAT &beta
) {
	// Aug25,2015 - Shaoda
	// This function can recover tilt (b) as small as 0.0001 degrees
	// It replaces a more complicated version in the code before Aug2015
    Vector<RFLOAT> v;

    // Make sure the vector is normalised
    v.resize(3);
    v = v0;
    v.normalise();

    // Tilt (b) should be [0, +180] degrees. Rot (a) should be [-180, +180] degrees
    alpha = degrees(atan2(v[1], v[0])); // 'atan2' returns an angle within [-pi, +pi] radians for rot
    beta  = degrees(acos(v[2])); // 'acos' returns an angle within [0, +pi] radians for tilt

    // The following is done to keep in line with the results from old codes
    // If tilt (b) = 0 or 180 degrees, sin(b) = 0, rot (a) cannot be calculated from the direction
    if (fabs(beta) < 0.001 || fabs(beta - 180.0) < 0.001) {alpha = 0.0;}

    return;

}

/* Matrix --> Euler angles ------------------------------------------------- */
#define CHECK
// #define DEBUG_EULER
angles_t Euler::matrix2angles(const Matrix<RFLOAT> &A) {

    if (A.ncols() != 3 || A.nrows() != 3)
        REPORT_ERROR("Euler::matrix2angles: The Euler matrix is not 3×3");

    RFLOAT alpha, beta, gamma;

    RFLOAT abs_sb = sqrt(A(0, 2) * A(0, 2) + A(1, 2) * A(1, 2));
    RFLOAT sign_sb;
    if (abs_sb > 16 * FLT_EPSILON) {
        gamma = atan2(A(1, 2), -A(0, 2));
        alpha = atan2(A(2, 1), +A(2, 0));
        if (abs(sin(gamma)) < FLT_EPSILON) {
            // cos(gamma) is very close to either +1 or -1
            sign_sb = sgn_nozero(-A(0, 2)) * sgn(cos(gamma));
        } else {
            // sin(gamma) is not zero
            sign_sb = sgn_nozero(A(1, 2)) * sgn(sin(gamma));
        }
        // if (sin(alpha) < FLT_EPSILON) {
        //     sign_sb = sgn_nozero(-A(0, 2) / cos(gamma));
        // } else {
        //     sign_sb = sgn(sin(alpha)) * sgn_nozero(A(2, 1));
        // }
        beta = atan2(abs_sb * sign_sb * abs_sb, A(2, 2));
    } else {
        if (A(2, 2) >= 0) {
            // Let's consider the matrix as a rotation around Z
            alpha = 0;
            beta  = 0;
            gamma = atan2(-A(1, 0), A(0, 0));
        } else {
            alpha = 0;
            beta  = PI;
            gamma = atan2(A(1, 0), -A(0, 0));
        }
    }

    gamma = degrees(gamma);
    beta  = degrees(beta);
    alpha = degrees(alpha);

    #ifdef DEBUG_EULER
    std::cout << "abs_sb " << abs_sb << std::endl;
    std::cout << "A(1, 2) " << A(1, 2) << " A(0, 2) " << A(0, 2) << " gamma " << gamma << std::endl;
    std::cout << "A(2, 1) " << A(2, 1) << " A(2, 0) " << A(2, 0) << " alpha " << alpha << std::endl;
    std::cout << "sign sb " << sign_sb << " A(2, 2) " << A(2, 2) << " beta "  << beta  << std::endl;
    #endif
    return { alpha, beta, gamma };
}
#undef CHECK

#ifdef NEVERDEFINED
// Michael's method
void Euler::matrix2angles(
    Matrix<RFLOAT> A, 
    RFLOAT *alpha, RFLOAT *beta, RFLOAT *gamma
) {

    RFLOAT abs_sb;
           if (abs(A(1, 1)) > FLT_EPSILON) {
        abs_sb = sqrt((-A(2, 2) * A(1, 2) * A(2, 1) - A(0, 2) * A(2, 0)) / A(1, 1));
    } else if (abs(A(0, 1)) > FLT_EPSILON) {
        abs_sb = sqrt((-A(2, 1) * A(2, 2) * A(0, 2) + A(2, 0) * A(1, 2)) / A(0, 1));
    } else if (abs(A(0, 0)) > FLT_EPSILON) {
        abs_sb = sqrt((-A(2, 0) * A(2, 2) * A(0, 2) - A(2, 1) * A(1, 2)) / A(0, 0));
    } else {
        EXIT_ERROR(1, "Don't know how to extract angles");
    }

    if (abs_sb > FLT_EPSILON) {
        *beta  = atan2(abs_sb, A(2, 2));
        *alpha = atan2(A(2, 1) / abs_sb, A(2, 0) / abs_sb);
        *gamma = atan2(A(1, 2) / abs_sb, -A(0, 2) / abs_sb);
    } else {
        *alpha = 0;
        *beta  = 0;
        *gamma = atan2(A(1, 0), A(0, 0));
    }

    *gamma = degrees(*gamma);
    *beta  = degrees(*beta);
    *alpha = degrees(*alpha);
}
#endif

/* Euler up-down correction ------------------------------------------------ */
angles_t Euler::up_down(RFLOAT rot, RFLOAT tilt, RFLOAT psi)  {
    return { rot, tilt + 180, -psi - 180 };
}

/* Same view, differently expressed ---------------------------------------- */
angles_t Euler::another_set(RFLOAT rot, RFLOAT tilt, RFLOAT psi) {
    return { rot + 180, -tilt, psi - 180 };
}

/* Euler mirror Y ---------------------------------------------------------- */
angles_t Euler::mirrorY(RFLOAT rot, RFLOAT tilt, RFLOAT psi) {
    return { rot, tilt + 180, -psi };
}

/* Euler mirror X ---------------------------------------------------------- */
angles_t Euler::mirrorX(RFLOAT rot, RFLOAT tilt, RFLOAT psi) {
    return { rot, tilt + 180, 180 - psi };
}

/* Euler mirror XY --------------------------------------------------------- */
angles_t Euler::mirrorXY(RFLOAT rot, RFLOAT tilt, RFLOAT psi) {
    return { rot, tilt, 180 + psi };
}

/* Apply a transformation matrix to Euler angles --------------------------- */
angles_t Euler::apply_transf(
    const Matrix<RFLOAT> &L, const Matrix<RFLOAT> &R,
    RFLOAT rot, RFLOAT tilt, RFLOAT psi
) {
    Matrix<RFLOAT> euler = Euler::angles2matrix(rot, tilt, psi);
    return Euler::matrix2angles(L.matmul(euler).matmul(R));
}

/* Rotate (3D) MultidimArray with 3 Euler angles ------------------------------------- */
Matrix<RFLOAT> Euler::rotation3DMatrix(RFLOAT rot, RFLOAT tilt, RFLOAT psi) {
    return Euler::angles2matrix(rot, tilt, psi, true);
}
