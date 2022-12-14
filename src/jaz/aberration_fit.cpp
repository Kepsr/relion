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

#include <src/jaz/aberration_fit.h>
#include <src/jaz/image_log.h>


OriginalBasis AberrationFit::fitBasic(
    Image<RFLOAT> phase, Image<RFLOAT> weight, double angpix
) {
    Matrix2D<RFLOAT> A = Matrix2D<RFLOAT>::zeros(5, 5);
    Matrix1D<RFLOAT> b {0, 0, 0, 0, 0};

    const int sh = phase.data.xdim;
    const int s  = phase.data.ydim;

    const double as = angpix * s;

    OriginalBasis basis;
    std::vector<double> vals (5);

    for (int yi = 0; yi < s; yi++)
    for (int xi = 0; xi < sh; xi++) {
        const double x = xi / as;
        const double y = (yi > sh ? yi - s : yi) / as;

        basis.getBasisValues(x, y, &vals[0]);

        const double v = phase(yi, xi);
        const double w = weight(yi, xi);

        for (int r = 0; r < 5; r++) {
            b[r] += w * w * vals[r] * v;

            for (int c = 0; c < 5; c++) {
                A(r, c) += w * w * vals[r] * vals[c];
            }
        }
    }

    const double tol = 1e-20;
    Matrix1D<RFLOAT> sol (5);
    solve(A, b, sol, tol);

    std::copy(sol.begin(), sol.end(), basis.coefficients.begin());

    return basis;
}

Image<RFLOAT> AberrationFit::draw(AberrationBasis *fit, double angpix, int s) {
    const int sh = s / 2 + 1;
    const double as = angpix * s;

    Image<RFLOAT> vis(s, sh);

    std::vector<double> vals(fit->coefficients.size(), 0.0);

    for (int yi = 0; yi < s;  yi++)
    for (int xi = 0; xi < sh; xi++) {
        const double x = xi / as;
        const double y = (yi > sh ? yi - s : yi) / as;

        fit->getBasisValues(x, y, &vals[0]);

        double v = 0.0;

        for (int i = 0; i < 5; i++) {
            v += fit->coefficients[i] * vals[i];
        }

        vis.data.elem(yi, xi) = v;
    }

    return vis;
}

AberrationBasis::AberrationBasis(int dims): coefficients(dims, 0.0) {}

void AberrationBasis::offsetCtf(MetaDataTable &mdt, int particle) {
    // identical to CTF::read() and CTF::initialize():
    double kV, DeltafU, DeltafV, azimuthal_angle, Cs, scale, Q0, phase_shift;

    #define TRYGETVALUE(variable, label, defaultvalue) \
    try { \
        variable = mdt.getValue<RFLOAT>(label, particle); \
    } catch (const char* errmsg) { \
        variable = defaultvalue; \
    }

    TRYGETVALUE(kV,              EMDL::CTF_VOLTAGE,        200);
    TRYGETVALUE(DeltafU,         EMDL::CTF_DEFOCUSU,       0);
    TRYGETVALUE(DeltafV,         EMDL::CTF_DEFOCUSV,       DeltafU);
    TRYGETVALUE(azimuthal_angle, EMDL::CTF_DEFOCUS_ANGLE,  0);
    TRYGETVALUE(Cs,              EMDL::CTF_CS,             0);
    // TRYGETVALUE(scale,           EMDL::CTF_SCALEFACTOR,    1);
    TRYGETVALUE(Q0,              EMDL::CTF_Q0,             0);
    // TRYGETVALUE(phase_shift,     EMDL::CTF_PHASESHIFT,     0);

    #undef TRYGETVALUE

    // std::cout << DeltafU << ", " << DeltafV << " @ " << azimuthal_angle << "°, " << Cs << ", " << Q0 << "\n";

    double local_Cs = Cs * 1e7;
    double local_kV = kV * 1e3;
    double rad_azimuth = radians(azimuthal_angle);

    double defocus_average   = (DeltafU + DeltafV) * -0.5;
    double defocus_deviation = (DeltafU - DeltafV) * -0.5;
    double lambda = 12.2643247 / sqrt(local_kV * (1.0 + local_kV * 0.978466e-6));

    double K1 = (PI / 2) * 2 * lambda;
    double K2 = (PI / 2) * local_Cs * lambda * lambda * lambda;
    double K3 = atan(Q0 / sqrt(1 - Q0 * Q0));

    _offsetCtf(
        local_Cs, lambda, rad_azimuth, defocus_average, defocus_deviation,
        K1, K2, K3, mdt, particle
    );
}

OriginalBasis::OriginalBasis(): AberrationBasis(5) {}

void OriginalBasis::getBasisValues(double x, double y, double *dest) {
    dest[0] = 1.0;           // phase shift
    dest[1] = x * x + y * y; // defocus
    dest[2] = x * x - y * y; // oblique astigmatism
    dest[3] = x * y;         // vertical astigmatism
    dest[4] = (x * x + y * y) * (x * x + y * y); // primary spherical
}

void OriginalBasis::_offsetCtf(
    double local_Cs, double lambda,
    double rad_azimuth, double defocus_average, double defocus_deviation,
    double K1, double K2, double K3, MetaDataTable &mdt, int particle
) {
    /* from ctf.h:

           RFLOAT u2 = X * X + Y * Y;
           RFLOAT u4 = u2 * u2;
           RFLOAT deltaf = defocus_average + defocus_deviation*cos(2 * (atan2(Y, X) - rad_azimuth))

           argument = K1 * deltaf * u2 + K2 * u4 - K5 - K3

                K1 = PI / 2 * 2 * lambda;
                K2 = PI / 2 * local_Cs * lambda * lambda * lambda;
                K3 = atan(Q0 / sqrt(1 - Q0 * Q0));
                K5 = radians(phase_shift);

                local_Cs = Cs * 1e7;

       astigmatism/defocus:

           K1 * deltaf * u2
           = K1 * defocus_average * u2 + defocus_deviation * K1 * cos(2 * (phi - rad_azimuth)) * u2
           = K1 * defocus_average * u2 + defocus_deviation * K1 * cos(2 * phi - 2 * rad_azimuth) * u2
           = K1 * defocus_average * u2 + defocus_deviation * K1 * [cos(2 * phi) * cos(2 * rad_azimuth) + sin(2 * phi) * sin(2 * rad_azimuth)] * u2
           = K1 * defocus_average * u2 + defocus_deviation * K1 * [(cos(phi) * cos(phi) - sin(phi) * sin(phi)) * cos(2 * rad_azimuth) + 2 sin(phi) * cos(phi) * sin(2 * rad_azimuth)] * u2
           = K1 * defocus_average * u2 + defocus_deviation * K1 * [(X * X - Y * Y) * cos(2 * rad_azimuth) + 2 * Y * X sin(2 * rad_azimuth)]
           = b1 (X * X + Y * Y) + b2 (X * X - Y * Y) + b3 (X * Y)

           where:  b1 =     K1 * defocus_average
                   b2 =     K1 * defocus_deviation * cos(2 * rad_azimuth)
                   b3 = 2 * K1 * defocus_deviation * sin(2 * rad_azimuth)

                   <=>

                   defocus_average = b1 / (PI * lambda)
                   defocus_deviation = sqrt(b2 * b2 + b3 * b3 / 4) / (PI * lambda)
                   rad_azimuth = atan2(b3 / 2, b2) / 2                        */

    double b1 =     K1 * defocus_average                          + coefficients[1];
    double b2 =     K1 * defocus_deviation * cos(2 * rad_azimuth) + coefficients[2];
    double b3 = 2 * K1 * defocus_deviation * sin(2 * rad_azimuth) + coefficients[3];

    double new_defocus_average   = b1 / (PI * lambda);
    double new_defocus_deviation = sqrt(b2 * b2 + b3 * b3 / 4)/(PI * lambda);
    double new_rad_azimuth       = atan2(b3 / 2.0, b2) / 2.0;

    double azimuthal_angle = degrees(new_rad_azimuth);
    double DeltafU = -new_defocus_deviation - new_defocus_average;
    double DeltafV = +new_defocus_deviation - new_defocus_average;

/*     spherical aberration:

           K2 * u4 = b4 * u4
            <=>
           PI / 2 * local_Cs * lambda * lambda * lambda = b4;
            <=>
           local_Cs = 2 * b4 / (PI * lambda * lambda * lambda)
            <=>
           Cs = 1e-7 * 2 * b4 / (PI * lambda * lambda * lambda)            */

    double b4 = PI * lambda * lambda * lambda * local_Cs / 2.0 + coefficients[4];
    double Cs = 1e-7 * 2.0 * b4 / (PI * lambda * lambda * lambda);

    /*     phase shift / amp. contrast:

               K3 = atan(Q0 / sqrt(1 - Q0 * Q0))
                <=>
               Q0 / sqrt(1 - Q0 * Q0) = tan(K3)
                <=>
               Q0 * Q0 = (1 - Q0 * Q0) * tan(K3) * tan(K3)
                <=>
               Q0 * Q0 * (1 + tan(K3) * tan(K3)) = tan(K3) * tan(K3)
                <=>
               Q0 = sqrt(tan(K3) * tan(K3) / (1 + tan(K3) * tan(K3)))
    */

    double b0 = K3 - coefficients[0];

    if (b0 < 0) {
        double phase_shift;
        try {
            phase_shift = mdt.getValue<double>(EMDL::CTF_PHASESHIFT, particle);
        } catch (const char *errmsg) {
            phase_shift = 0;
        }

        phase_shift = phase_shift - degrees(coefficients[0]);

        mdt.setValue(EMDL::CTF_PHASESHIFT, phase_shift, particle);
    } else {
        double t0 = tan(b0);
        double Q0 = sqrt(t0 * t0 / (1 + t0 * t0));

        mdt.setValue(EMDL::CTF_Q0, Q0, particle);
    }

    mdt.setValue(EMDL::CTF_DEFOCUSU,      DeltafU,         particle);
    mdt.setValue(EMDL::CTF_DEFOCUSV,      DeltafV,         particle);
    mdt.setValue(EMDL::CTF_DEFOCUS_ANGLE, azimuthal_angle, particle);
    mdt.setValue(EMDL::CTF_CS,            Cs,              particle);

    //std::cout << DeltafU << ", " << DeltafV << " @ " << azimuthal_angle << "°, " << Cs << ", " << Q0 << "\n\n";
}
