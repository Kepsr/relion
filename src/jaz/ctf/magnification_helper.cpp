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

#include "magnification_helper.h"

#include <src/jaz/slice_helper.h>
#include "src/jaz/ctf_helper.h"
#include <src/projector.h>
#include <src/jaz/img_proc/filter_helper.h>
#include <src/jaz/optimization/nelder_mead.h>
#include <src/jaz/gravis/t4Matrix.h>
#include <src/jaz/d3x3/dsyev2.h>

using namespace gravis;

Matrix<RFLOAT> MagnificationHelper::polarToMatrix(
    double scaleMajor, double scaleMinor, double angleDeg
) {
    // based on definition by T. Nakane

    Matrix<RFLOAT> out(2, 2);

    const double angle = radians(angleDeg);
    const double si = sin(angle), co = cos(angle);
    const double si2 = si * si, co2 = co * co;

    /*
         Out = Rot(angle) * Diag(scale_major, scale_minor) * Rot(-angle), where
        Rot(angle) = [[cos(angle), -sin(angle)], [sin(angle), cos(angle)]]

        [  c s ] [ j 0 ] [ c -s ]
        [ -s c ] [ 0 n ] [ s  c ]
        =
        [  c s ] [ jc -js ]
        [ -s c ] [ ns  nc ]
        =
        [  jcc+nss -jcs+ncs ]
        [ -jcs+ncs  jss+ncc ]
    */

    out(0, 0) =  scaleMajor * co2 + scaleMinor * si2;
    out(1, 1) =  scaleMajor * si2 + scaleMinor * co2;
    out(0, 1) = (-scaleMajor + scaleMinor) * si * co;
    out(1, 0) = out(0, 1);

    return out;
}

void MagnificationHelper::matrixToPolar(
    const Matrix<RFLOAT>& mat,
    RFLOAT& scaleMajor, RFLOAT& scaleMinor, RFLOAT& angleDeg
) {
    matrixToPolar(
        d2Matrix(mat(0, 0), mat(0, 1), mat(1, 0), mat(1, 1)),
        scaleMajor, scaleMinor, angleDeg
    );
}

void MagnificationHelper::matrixToPolar(
    const d2Matrix& mat,
    RFLOAT& scaleMajor, RFLOAT& scaleMinor, RFLOAT& angleDeg
) {
    const double m00 = mat(0, 0);
    const double m11 = mat(1, 1);
    const double m01 = 0.5 * (mat(0, 1) + mat(1, 0));

    double ev0, ev1, cs, sn;

    dsyev2(m00, m01, m11, &ev0, &ev1, &cs, &sn);

    scaleMajor = ev0;
    scaleMinor = ev1;
    angleDeg = degrees(atan2(sn,cs));

    return;
}

void MagnificationHelper::updateScaleFreq(
    const Image<Complex> &prediction,
    const Volume<t2Vector<Complex>>& predGradient,
    const Image<Complex> &observation,
    CTF &ctf,  ObservationModel *obsModel,
    double angpix, Volume<Equation2x2> &eqs,
    bool do_ctf_padding
) {
    const long w = prediction.data.xdim;
    const long h = prediction.data.ydim;

    /*Volume<gravis::d2Vector> gradReal(w, h,1), gradImg(w, h,1);

    FilterHelper::centralGrad2D(prediction, gradReal, gradImg);*/

    Image<RFLOAT> ctfImg(w, h);
    ctfImg() = CtfHelper::getFftwImage(
        ctf,
        w, h, h, h, angpix,
        obsModel,
        false, false, false, true, do_ctf_padding
    );

    for (long y = 0; y < h; y++)
    for (long x = 0; x < w; x++) {
        Complex vx = direct::elem(prediction .data, x, y);
        Complex vy = direct::elem(observation.data, x, y);

        double c = ctfImg(y,x);

        gravis::d2Vector gr(predGradient(x, y, 0).x.real, predGradient(x, y, 0).y.real);
        gravis::d2Vector gi(predGradient(x, y, 0).x.imag, predGradient(x, y, 0).y.imag);

        eqs(x, y, 0).Axx += c * c * (gr.x * gr.x + gi.x * gi.x);
        eqs(x, y, 0).Axy += c * c * (gr.x * gr.y + gi.x * gi.y);
        eqs(x, y, 0).Ayy += c * c * (gr.y * gr.y + gi.y * gi.y);

        eqs(x, y, 0).bx += c * (gr.x * (vy.real - c * vx.real) + gi.x * (vy.imag - c * vx.imag));
        eqs(x, y, 0).by += c * (gr.y * (vy.real - c * vx.real) + gi.y * (vy.imag - c * vx.imag));
    }
}

void MagnificationHelper::updateScaleReal(
    const Image<Complex> &prediction,
    const Image<Complex> &observation,
    const Image<RFLOAT>& snr,
    CTF &ctf, ObservationModel *obsModel,
    double angpix,
    Volume<Equation2x2> &eqs,
    bool do_ctf_padding
) {
    const long ww = 2 * (observation.data.xdim - 1);
    const long w = prediction.data.xdim;
    const long h = prediction.data.ydim;

    Image<Complex> pred2(w, h), obs2(w, h);

    Image<RFLOAT> ctfImg(w, h);
    ctfImg() = CtfHelper::getFftwImage(
        ctf, w, h, h, h, angpix, obsModel,
        false, false, false, true, do_ctf_padding
    );

    for (long y = 0; y < h; y++)
    for (long x = 0; x < w; x++) {
        Complex vx = direct::elem(prediction.data, x, y);
        Complex vy = direct::elem(observation.data, x, y);
        Complex sn = direct::elem(snr.data, x, y);
        double c = ctfImg(y, x);

        direct::elem(pred2.data, x, y) = sn * c * vx;
        direct::elem(obs2 .data, x, y) = sn * vy;
    }

    Image<RFLOAT> realPred(ww, h), realObs(ww, h);
    realPred.data = FourierTransformer{}.inverseFourierTransform(pred2.data);
    realObs.data  = FourierTransformer{}.inverseFourierTransform(obs2.data);

    Volume<gravis::d2Vector> grad(ww, h, 1);
    FilterHelper::centralGrad2D(realPred, grad);

    for (long y = 0; y < h; y++)
    for (long x = 0; x < ww; x++) {
        double vx = direct::elem(realPred.data, x, y);
        double vy = direct::elem(realObs .data, x, y);

        gravis::d2Vector g = grad(x, y, 0);

        eqs(x, y, 0).Axx += g.x * g.x;
        eqs(x, y, 0).Axy += g.x * g.y;
        eqs(x, y, 0).Ayy += g.y * g.y;

        eqs(x, y, 0).bx += g.x * (vy - vx);
        eqs(x, y, 0).by += g.y * (vy - vx);
    }
}

void MagnificationHelper::solvePerPixel(
    const Volume<Equation2x2> &eqs,
    Image<RFLOAT> &vx, Image<RFLOAT> &vy
) {
    const long w = eqs.dimx;
    const long h = eqs.dimy;

    vx = Image<RFLOAT>(w, h);
    vy = Image<RFLOAT>(w, h);

    for (long y = 0; y < h; y++)
    for (long x = 0; x < w; x++) {
        Equation2x2 eq = eqs(x, y, 0);

        gravis::d2Vector b(eq.bx, eq.by);
        gravis::d2Matrix A;
        A(0, 0) = eq.Axx;
        A(0, 1) = eq.Axy;
        A(1, 0) = eq.Axy;
        A(1, 1) = eq.Ayy;

        double det = A(0, 0) * A(1, 1) - A(0, 1) * A(1, 0);

        if (det == 0.0) {
            direct::elem(vx.data, x, y) = 0.0;
            direct::elem(vy.data, x, y) = 0.0;
        } else {
            gravis::d2Matrix Ai = A;
            Ai.invert();

            gravis::d2Vector xx = Ai * b;

            direct::elem(vx.data, x, y) = xx.x;
            direct::elem(vy.data, x, y) = xx.y;
        }
    }
}

Matrix<RFLOAT> MagnificationHelper::solveLinearlyFreq(
    const Volume<Equation2x2> &eqs,
    const Image<RFLOAT>& snr,
    Image<RFLOAT> &vx, Image<RFLOAT> &vy
) {
    Matrix<RFLOAT> mat(2, 2);

    const long w = eqs.dimx;
    const long h = eqs.dimy;

    vx = Image<RFLOAT>(w, h);
    vy = Image<RFLOAT>(w, h);

    d4Vector b(0.0, 0.0, 0.0, 0.0);

    d4Matrix A(
        0.0, 0.0, 0.0, 0.0,
        0.0, 0.0, 0.0, 0.0,
        0.0, 0.0, 0.0, 0.0,
        0.0, 0.0, 0.0, 0.0
    );

    for (long yi = 1; yi < h - 1; yi++)
    for (long xi = 1; xi < w - 1; xi++) {
        Equation2x2 eq = eqs(xi, yi, 0);
        double wgh = direct::elem(snr.data, xi, yi);

        const double x = xi;
        const double y = yi < w ? yi : yi - h;

        A(0, 0) += wgh * x * x * eq.Axx;
        A(0, 1) += wgh * x * y * eq.Axx;
        A(0, 2) += wgh * x * x * eq.Axy;
        A(0, 3) += wgh * x * y * eq.Axy;

        A(1, 0) += wgh * x * y * eq.Axx;
        A(1, 1) += wgh * y * y * eq.Axx;
        A(1, 2) += wgh * x * y * eq.Axy;
        A(1, 3) += wgh * y * y * eq.Axy;

        A(2, 0) += wgh * x * x * eq.Axy;
        A(2, 1) += wgh * x * y * eq.Axy;
        A(2, 2) += wgh * x * x * eq.Ayy;
        A(2, 3) += wgh * x * y * eq.Ayy;

        A(3, 0) += wgh * x * y * eq.Axy;
        A(3, 1) += wgh * y * y * eq.Axy;
        A(3, 2) += wgh * x * y * eq.Ayy;
        A(3, 3) += wgh * y * y * eq.Ayy;

        b[0] += wgh * x * eq.bx;
        b[1] += wgh * y * eq.bx;
        b[2] += wgh * x * eq.by;
        b[3] += wgh * y * eq.by;
    }

    d4Matrix Ai = A;
    Ai.invert();

    d4Vector opt = Ai * b;

    mat(0, 0) = opt[0] + 1.0;
    mat(0, 1) = opt[1];
    mat(1, 0) = opt[2];
    mat(1, 1) = opt[3] + 1.0;

	// std::cout << opt[0] << ", " << opt[1] << "\n"
    //           << opt[2] << ", " << opt[3] << "\n";

    for (long yi = 0; yi < h; yi++)
    for (long xi = 0; xi < w; xi++) {
        const double x = xi;
        const double y = yi < w ? yi : yi - h;

        direct::elem(vx.data, xi, yi) = opt[0] * x + opt[1] * y;
        direct::elem(vy.data, xi, yi) = opt[2] * x + opt[3] * y;
    }

    return mat;
}

void MagnificationHelper::readEQs(std::string path, Volume<Equation2x2> &eqs) {

    Image<RFLOAT> Axx = Image<RFLOAT>::from_filename(path + "_Axx.mrc");
    Image<RFLOAT> Axy = Image<RFLOAT>::from_filename(path + "_Axy.mrc");
    Image<RFLOAT> Ayy = Image<RFLOAT>::from_filename(path + "_Ayy.mrc");

    Image<RFLOAT> bx  = Image<RFLOAT>::from_filename(path + "_bx.mrc");
    Image<RFLOAT> by  = Image<RFLOAT>::from_filename(path + "_by.mrc");

    const long w = Axx.data.xdim;
    const long h = Axx.data.ydim;

    eqs.resize(w, h, 1);

    for (long y = 0; y < h; y++)
    for (long x = 0; x < w; x++) {
        eqs(x, y, 0).Axx = direct::elem(Axx.data, x, y);
        eqs(x, y, 0).Axy = direct::elem(Axy.data, x, y);
        eqs(x, y, 0).Ayy = direct::elem(Ayy.data, x, y);

        eqs(x, y, 0).bx = direct::elem(bx.data, x, y);
        eqs(x, y, 0).by = direct::elem(by.data, x, y);
    }
}

void MagnificationHelper::writeEQs(
    const Volume<Equation2x2> &eqs, std::string path
) {
    const long w = eqs.dimx;
    const long h = eqs.dimy;

    Image<RFLOAT> Axx(w, h);
    Image<RFLOAT> Axy(w, h);
    Image<RFLOAT> Ayy(w, h);
    Image<RFLOAT> bx(w, h);
    Image<RFLOAT> by(w, h);

    for (long y = 0; y < h; y++)
    for (long x = 0; x < w; x++) {
        Equation2x2 eq = eqs(x, y, 0);

        direct::elem(Axx.data, x, y) = eq.Axx;
        direct::elem(Axy.data, x, y) = eq.Axy;
        direct::elem(Ayy.data, x, y) = eq.Ayy;

        direct::elem(bx.data, x, y) = eq.bx;
        direct::elem(by.data, x, y) = eq.by;
    }

    Axx.write(path + "_Axx.mrc");
    Axy.write(path + "_Axy.mrc");
    Ayy.write(path + "_Ayy.mrc");
    bx.write(path + "_bx.mrc");
    by.write(path + "_by.mrc");
}

void MagnificationHelper::updatePowSpec(
    const Image<Complex> &prediction,
    const Image<Complex> &observation,
    CTF &ctf, ObservationModel *obsModel,
    double angpix,
    Image<RFLOAT> &powSpecPred,
    Image<RFLOAT> &powSpecObs,
    bool do_ctf_padding
) {
    const long w = prediction.data.xdim;
    const long h = prediction.data.ydim;

    Image<RFLOAT> ctfImg(w, h);
    ctfImg() = CtfHelper::getFftwImage(ctf, w, h, h, h, angpix, obsModel, false, false, false, true, do_ctf_padding);

    for (long y = 0; y < h; y++)
    for (long x = 0; x < w; x++) {
        const double xf = x;
        const double yf = y < w ? y : y - h;

        Complex vx = direct::elem(prediction.data, x, y);
        Complex vy = direct::elem(observation.data, x, y);
        double c = ctfImg(y,x);

        direct::elem(powSpecPred.data, x, y) += (c * vx).abs();
        direct::elem(powSpecObs .data, x, y) += vy.abs();
    }
}

void MagnificationHelper::adaptAstigmatism(
    const std::vector<Matrix<RFLOAT>>& dMs,
    std::vector<MetaDataTable>& partMdts,
    bool perParticle, ObservationModel* obsModel
) {
    const long int mc = partMdts.size();
    const int ogc = dMs.size();

    std::vector<d2Matrix> M(ogc), Mi(ogc), Mit(ogc);

    for (int og = 0; og < ogc; og++) {
        M[og] = d2Matrix(dMs[og](0, 0), dMs[og](0, 1), dMs[og](1, 0), dMs[og](1, 1));

        Mi[og] = M[og];
        Mi[og].invert();

        Mit[og] = Mi[og];
        Mit[og].transpose();
    }

    for (long int m = 0; m < mc; m++) {
        const int pc = partMdts[m].size();

        std::vector<d2Matrix> A(pc), D(pc), Q(pc);

        for (long int p = 0; p < pc; p++) {

            double deltafU = partMdts[m].getValue<double>(EMDL::CTF_DEFOCUSU, p);
            double deltafV = partMdts[m].getValue<double>(EMDL::CTF_DEFOCUSV, p);
            double phiDeg  = partMdts[m].getValue<double>(EMDL::CTF_DEFOCUS_ANGLE, p);

            const double phi = radians(phiDeg);

            const double si = sin(phi);
            const double co = cos(phi);

            Q[p] = d2Matrix(co, si, -si, co);
            D[p] = d2Matrix(-deltafU, 0.0, 0.0, -deltafV);
            d2Matrix Qt(co, -si, si, co);

            A[p] = Qt * D[p] * Q[p];
        }

        if (perParticle) {
            for (long int p = 0; p < pc; p++) {
                int og = partMdts[m].getValue<int>(EMDL::IMAGE_OPTICS_GROUP, p) - 1;

                d2Matrix A2 = Mit[og] * A[p] * Mi[og];

                RFLOAT deltafU_neg, deltafV_neg, phiDeg;
                matrixToPolar(A2, deltafU_neg, deltafV_neg, phiDeg);

                partMdts[m].setValue(EMDL::CTF_DEFOCUSU, -deltafU_neg, p);
                partMdts[m].setValue(EMDL::CTF_DEFOCUSV, -deltafV_neg, p);
                partMdts[m].setValue(EMDL::CTF_DEFOCUS_ANGLE, phiDeg, p);
            }
        } else {
            // keep difference between deltafU and deltafV, as well as the azimuth angle,
            // constant for all particles in the same micrograph and optics group
            std::vector<int> optGroups = obsModel->getOptGroupsPresent(partMdts[m]);
            const int cc = optGroups.size();

            std::vector<int> groupToIndex(obsModel->numberOfOpticsGroups() + 1, -1);

            for (int i = 0; i < cc; i++) {
                groupToIndex[optGroups[i] + 1] = i;
            }

            for (int c = 0; c < cc; c++) {
                const int og = optGroups[c];

                d2Matrix A_mean(0.0, 0.0, 0.0, 0.0);

                for (long int p = 0; p < pc; p++) {
                    int ogp = partMdts[m].getValue<int>(EMDL::IMAGE_OPTICS_GROUP, p) - 1;

                    if (ogp == og) {
                        A_mean += A[p] * (1.0 / (double) pc);
                    }
                }


                A_mean = Mit[og] * A_mean * Mi[og];

                double deltafU_mean_neg, deltafV_mean_neg, co, si;
                dsyev2(
                    A_mean(0, 0), A_mean(1, 0), A_mean(1, 1),
                    &deltafU_mean_neg, &deltafV_mean_neg, &co, &si
                );

                d2Matrix Q2 (co, si, -si, co);
                d2Matrix Qt2(co, -si, si, co);

                double meanDef_mean = 0.5 * (deltafU_mean_neg + deltafV_mean_neg);

                for (long int p = 0; p < pc; p++) {
                    int ogp = partMdts[m].getValue<int>(EMDL::IMAGE_OPTICS_GROUP, p) - 1;

                    if (ogp == og) {
                        d2Matrix Ap2 = Mit[og] * A[p] * Mi[og];

                        double deltafU_p_neg, deltafV_p_neg, cop, sip;
                        dsyev2(
                            Ap2(0,0), Ap2(1,0), Ap2(1,1),
                            &deltafU_p_neg, &deltafV_p_neg, &cop, &sip
                        );

                        double meanDef_p = 0.5 * (deltafU_p_neg + deltafV_p_neg);

                        d2Matrix Dp2(
                            deltafU_mean_neg - meanDef_mean + meanDef_p, 0.0,
                            0.0, deltafV_mean_neg - meanDef_mean + meanDef_p
                        );

                        d2Matrix Apa2 = Qt2 * Dp2 * Q2;

                        RFLOAT deltafU_pa_neg, deltafV_pa_neg, phiDeg;
                        matrixToPolar(Apa2, deltafU_pa_neg, deltafV_pa_neg, phiDeg);

                        partMdts[m].setValue(EMDL::CTF_DEFOCUSU, -deltafU_pa_neg, p);
                        partMdts[m].setValue(EMDL::CTF_DEFOCUSV, -deltafV_pa_neg, p);
                        partMdts[m].setValue(EMDL::CTF_DEFOCUS_ANGLE, phiDeg, p);
                    }
                }
            }
        }
    }
}

Equation2x2::Equation2x2(): Axx(0.0), Axy(0.0), Ayy(0.0), bx(0.0), by(0.0) {}
