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

#include "src/ctf.h"
#include "src/args.h"
#include "src/fftw.h"
#include <src/jaz/gravis/t2Matrix.h>

using namespace gravis;

void CTF::initialise() {
    // Change units
    RFLOAT local_Cs = Cs * 1e7;
    RFLOAT local_kV = kV * 1e3;
    rad_azimuth = radians(azimuthal_angle);

    // Average focus and deviation
    defocus_average   = -(DeltafU + DeltafV) * 0.5;
    defocus_deviation = -(DeltafU - DeltafV) * 0.5;

    // lambda=h/sqrt(2*m*e*kV)
    //    h: Planck constant
    //    m: electron mass
    //    e: electron charge
    // lambda=0.387832/sqrt(kV * (1.+0.000978466*kV)); // Hewz: Angstroms
    // lambda=h/sqrt(2*m*e*kV)
    lambda = 12.2643247 / sqrt(local_kV * (1.0 + local_kV * 0.978466e-6));
    // See http://en.wikipedia.org/wiki/Electron_diffraction

    // Helpful constants
    // ICE: X(u)=-PI/2*deltaf(u)*lambda*u^2+PI/2*Cs*lambda^3*u^4
    //          = K1*deltaf(u)*u^2         +K2*u^4
    K1 = PI / 2 * 2 * lambda;
    K2 = PI / 2 * local_Cs * lambda * lambda * lambda;
    K3 = atan(Q0 / sqrt(1 - Q0 * Q0));
    K4 = -Bfac / 4;

    // Phase shift in radians
    K5 = radians(phase_shift);

    if (Q0 < 0 || Q0 > 1)
        REPORT_ERROR("CTF::initialise ERROR: AmplitudeContrast Q0 cannot be smaller than zero or larger than one!");

    if (abs(DeltafU) < 1e-6 && abs(DeltafV) < 1e-6 && abs(Q0) < 1e-6 && abs(Cs) < 1e-6)
        REPORT_ERROR("CTF::initialise: ERROR: CTF initialises to all-zero values. Was a correct STAR file provided?");

    // express astigmatism as a bilinear form:

    const double sin_az = sin(rad_azimuth);
    const double cos_az = cos(rad_azimuth);

    d2Matrix Q (cos_az, +sin_az, -sin_az, cos_az);
    d2Matrix Qt(cos_az, -sin_az, +sin_az, cos_az);
    d2Matrix D (-DeltafU, 0.0, 0.0, -DeltafV);

    d2Matrix A = Qt * D * Q;

    Axx = A(0, 0);
    Axy = A(0, 1);
    Ayy = A(1, 1);
}

RFLOAT CTF::getGamma(RFLOAT X, RFLOAT Y) const {
    const RFLOAT u2 = X * X + Y * Y;
    const RFLOAT u4 = u2 * u2;
    return K1 * astigDefocus(X, Y) + K2 * u4 - K5 - K3;
}

RFLOAT CTF::getCtfFreq(RFLOAT X, RFLOAT Y) {
    const RFLOAT u2 = X * X + Y * Y;
    const RFLOAT u = sqrt(u2);
    const RFLOAT deltaf = getDeltaF(X, Y);
    return 2 * K1 * deltaf * u + 4 * K2 * u * u * u;
}

t2Vector<RFLOAT> CTF::getGammaGrad(RFLOAT X, RFLOAT Y) const {

    const RFLOAT u2 = X * X + Y * Y;
    // RFLOAT u4 = u2 * u2;

    // u4 = (X² + Y²)²
    // du4/dx = 2 (X² + Y²) 2 X = 4 (X³ + XY²) = 4 u2 X

    return t2Vector<RFLOAT>(
        2 * (K1 * Axx * X + K1 * Axy * Y + 2 * K2 * u2 * X),
        2 * (K1 * Ayy * Y + K1 * Axy * X + 2 * K2 * u2 * Y)
    );
}

std::vector<double> CTF::getK() {
    // offset by one to maintain indices (K[1] = K1)
    return {0, K1, K2, K3, K4, K5};
}

double CTF::getAxx() {
    return Axx;
}

double CTF::getAxy() {
    return Axy;
}

double CTF::getAyy() {
    return Ayy;
}
