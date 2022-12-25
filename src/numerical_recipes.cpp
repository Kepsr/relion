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
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <limits>
#include <algorithm>
#include "src/complex.h"
#include "src/numerical_recipes.h"

/* NUMERICAL UTILITIES ----------------------------------------------------- */
namespace nr {

    void error(const char error_text[]) {
        fprintf(stderr, "Numerical Recipes run-time error...\n");
        fprintf(stderr, "%s\n", error_text);
        fprintf(stderr, "...now exiting to system...\n");
        exit(1);
    }

}

// BESSEL FUNCTIONS --------------------------------------------------------
/* CO: They may not come in the numerical recipes but it is not a bad
   idea to put them here, in fact they come from Gabor's group in Feb'84     */
//  Bessel function of the first kind for order 0 (J_{0})
RFLOAT bessj0(RFLOAT x) {
    const RFLOAT ax = fabs(x);
    if (ax < 8) {
        const RFLOAT x2 = x * x;
        // Polynomial in x2
        return (
                     5.7568490574e10
            + x2 * (-1.3362590354e10
            + x2 * ( 6.516196407e8
            + x2 * (-1.121442418e7
            + x2 * ( 7.739233017e4
            + x2 *  -1.849052456e2
        ))))) / (
                    5.7568490411e10
            + x2 * (1.029532985e9
            + x2 * (9.494680718e6
            + x2 * (5.927264853e4
            + x2 * (2.678532712e2
            + x2
        )))));
    } else {
        const RFLOAT z = 8.0 / ax;
        const RFLOAT z2 = z * z;
        const RFLOAT xx = ax - 0.785398164;  // ax - 45 degrees
        return sqrt(0.636619772 / ax) * (cos(xx) * (1
            + z2 * (-0.1098628627e-2
            + z2 * ( 0.2734510407e-4
            + z2 * (-0.2073370639e-5
            + z2 *   0.2093887211e-6
        )))) - z * sin(xx) * (
                    -0.1562499995e-1
            + z2 * ( 0.1430488765e-3
            + z2 * (-0.6911147651e-5
            + z2 * ( 0.7621095161e-6
            + z2 *  -0.934935152e-7
        )))));
    }
}

//............................................................................
// Modified Bessel function of the first kind for order 0
RFLOAT bessi0(RFLOAT x) {
    RFLOAT ans;
    const RFLOAT ax = fabs(x);
    if (ax < 3.75) {
        const RFLOAT y = (x * x) / (3.75 * 3.75);
        ans = 1.0 
            + y * (3.5156229 
            + y * (3.0899424 
            + y * (1.2067492
            + y * (2.659732e-1 
            + y * (3.60768e-2
            + y *  4.5813e-3)))));
    } else {
        const RFLOAT y = 3.75 / ax;
        ans = (exp(ax) / sqrt(ax)) * (0.39894228
            + y * (+1.328592e-2
            + y * (+2.25319e-3 
            + y * (-1.57565e-3
            + y * (+9.16281e-3
            + y * (-2.057706e-2
            + y * (+2.635537e-2
            + y * (-1.647633e-2
            + y *  +3.92377e-3))))))));
    }
    return ans;
}

//............................................................................
RFLOAT bessi1(RFLOAT x) {
    RFLOAT ans;
    const RFLOAT ax = fabs(x);
    if (ax < 3.75) {
        const RFLOAT y = (x * x) / (3.75 * 3.75);
        ans = ax * (0.5
        + y * (0.87890594
        + y * (0.51498869
        + y * (0.15084934
        + y * (0.2658733e-1
        + y * (0.301532e-2
        + y *  0.32411e-3))))));
    } else {
        const RFLOAT y = 3.75 / ax;
        ans =       0.2282967e-1
            + y * (-0.2895312e-1
            + y * (+0.1787654e-1
            + y *  -0.420059e-2));
        ans =       0.39894228
            + y * (-0.3988024e-1
            + y * (-0.362018e-2
            + y * ( 0.163801e-2
            + y * (-0.1031555e-1
            + y * ans))));
        ans *= exp(ax) / sqrt(ax);
    }
    return copysign(ans, x);
}

/* General Bessel functions ------------------------------------------------ */
RFLOAT chebev(RFLOAT a, RFLOAT b, const RFLOAT c[], int m, RFLOAT x) {
    if ((x - a) * (x - b) > 0)
        nr::error("x not in range in routine chebev");

    const RFLOAT y = (2 * x - a - b) / (b - a);
    const RFLOAT two_y = 2 * y;
    RFLOAT d = 0, dd = 0, sv;
    for (int j = m - 1; j >= 1; j--) {
        sv = d;
        d  = two_y * d - dd + c[j];
        dd = sv;
    }
    return y * d - dd + 0.5 * c[0];
}

void beschb(RFLOAT x, RFLOAT& gam1, RFLOAT& gam2, RFLOAT& gamma_1add, RFLOAT& gamma_1sub) {
    static const RFLOAT c1[] {
        -1.142022680371172e0, 6.516511267076e-3,
        3.08709017308e-4, -3.470626964e-6, 6.943764e-9,
        3.6780e-11, -1.36e-13
    };
    static const RFLOAT c2[] {
        1.843740587300906e0, -7.6852840844786e-2,
        1.271927136655e-3, -4.971736704e-6, -3.3126120e-8,
        2.42310e-10, -1.70e-13, -1.0e-15
    };
    const int NUSE1 = 5, NUSE2 = 5;
    RFLOAT xx = 8 * x * x - 1;
    gam1 = chebev(-1, 1, c1, NUSE1, xx);
    gam2 = chebev(-1, 1, c2, NUSE2, xx);
    gamma_1add = gam2 - x * gam1;
    gamma_1sub = gam2 + x * gam1;
}

RFLOAT inline at_least(RFLOAT x, RFLOAT min) {
    return fabs(x) < min ? min : x;
}

Complex inline at_least(Complex x, RFLOAT min) {
    return fabs(x.real) + fabs(x.imag) < min ? Complex(min, x.imag) : x;
}

void besseljy(RFLOAT x, RFLOAT nu, RFLOAT& J, RFLOAT& Y, RFLOAT& Jp, RFLOAT& Yp) {
    if (x <= 0 || nu < 0)
        nr::error("bad arguments in besseljy");
    const int MAXIT = 10000;
    const RFLOAT XMIN = 2, EPS = std::numeric_limits<RFLOAT>::epsilon(), FPMIN = 1e-30;
    const int nl = x < XMIN ? nu + 0.5 : std::max<int>(0, nu - x + 1.5);
    // nl is the number of downward recurrences of the J's and upward recurrences of Y's.
    // For x < XMIN, mu lies in the interval (-0.5, +0.5).
    // For x >= XMIN, it is chosen so that x is greater than the turning point.
    const RFLOAT mu = nu - nl, mumu = mu * mu;
    const RFLOAT ix = 1 / x, two_ix = 2 / x, half_x = x / 2;
    const RFLOAT W = two_ix / PI;  // The Wronskian

    /** Evaluate CF1 by modified Lentz's method
     *
     * J'{nu} / J{nu} = nu / x - J{nu + 1} / J{nu}
     *                = nu / x - 1 /
     *                  2 * (nu + 1) / x - 1 /
     *                  2 * (nu + 2) / x - ...
     *
     * The partial numerators (a) are all -1.
     * The partial denominators (b) are given by 2 * (nu + i) / x for i = 1, 2, 3...
     */
    RFLOAT CF1 = at_least(nu * ix, FPMIN);
    RFLOAT b = nu * two_ix, C = CF1, D = 0, delta;
    int nr_sign_changes = 0;
    {
    int i = 0;
    for (; i != MAXIT; i++) {
        // a = -1
        b += two_ix;
        C = at_least(b - 1 / C, FPMIN);  // b + a / C
        D = 1 / at_least(b - D, FPMIN);  // 1 / (b + a * D)
        if (D < 0) ++nr_sign_changes;
        CF1 *= delta = C * D;
        if (fabs(delta - 1) < EPS) break;
    }
    if (i == MAXIT)
        nr::error("x too large in besseljy; try asymptotic expansion");
    }

    // Initialize Jnu and J'nu for downward recurrence.
    RFLOAT Jl = nr_sign_changes % 2 ? -FPMIN : FPMIN, Jpl = Jl * CF1;
    // Store values for later rescaling.
    const RFLOAT Jl1 = Jl, Jp1 = Jpl;
    // Downward recurrence
    RFLOAT nix = nu * ix, Jll;  // From nu/x to mu/x
    for (int i = 0; i < nl; i++) {
        Jll = nix * Jl + Jpl;  // J{nu-1} = nu/x J{nu} + J'{nu}
        Jpl = Jll * (nix -= ix) - Jl;  // J'{nu-1} = (nu/x - 1/x) J{nu-1} - J{nu}
        Jl  = Jll;
    }
    if (Jl == 0) Jl = EPS;
    const RFLOAT fmu = Jpl / Jl;  // Now have unnormalized J{mu} and J'{mu}
    RFLOAT Jmu, Ymu, Y1;
    if (x < XMIN) {  // Use series
        const RFLOAT pimu = PI * mu;
        const RFLOAT isinc_pimu = fabs(pimu) < EPS ? 1 : pimu / sin(pimu);  // 1 / sinc(pimu)
        const RFLOAT half_pimu = 0.5 * pimu;
        const RFLOAT sinc_half_pimu = fabs(half_pimu) < EPS ? 1 : sin(half_pimu) / half_pimu;
        const RFLOAT log_two_over_x = log(two_ix);
        /*
         * sigma = nu ln(2/x)
         * Gamma{1}(nu) = 1/2nu [ 1/Gamma(1 - nu) - 1/Gamma(1 + nu) ]
         * Gamma{2}(nu) = 1/2   [ 1/Gamma(1 - nu) + 1/Gamma(1 + nu) ]
         */
        const RFLOAT sigma = mu * log_two_over_x;
        const RFLOAT exp_sigma = exp(sigma);  // (x/2)**-mu
        const RFLOAT cosh_sigma = cosh(sigma);
        const RFLOAT sinch_sigma = fabs(sigma) < EPS ? 1 : sinh(sigma) / sigma;

        // Evaluate Gamma{1} and Gamma{2} by Chebyshev expansion for |x| <= 1/2
        // Also calculates 1/Gamma(1 + x) and 1/Gamma(1 - x)
        RFLOAT G1, G2, iG1addx, iG1subx;
        beschb(mu, G1, G2, iG1addx, iG1subx);

        const RFLOAT r = PI * half_pimu * sinc_half_pimu * sinc_half_pimu;  // 2/nu sin2(nu pi/2)
        const RFLOAT neg_half_x_sq = half_x * -half_x;
        RFLOAT f = 2 / PI * isinc_pimu * (cosh_sigma * G1 + sinch_sigma * log_two_over_x * G2);
        RFLOAT p = exp_sigma / (iG1addx * PI);
        RFLOAT q = 1 / (exp_sigma * iG1subx * PI);
        RFLOAT g = f + r * q, h = p, c = 1;
        RFLOAT S1 = g, S2 = h, delta1, delta2;
        int k = 1;
        for (; k != MAXIT; k++) {
            f = (k * f + p + q) / (k * k - mumu);  // f{k} = k f{k-1} + p{k-1} + q{k-1}
            p /= k + mu;    // p{k} = p{k-1} / (k - nu)
            q /= k - mu;    // q{k} = q{k-1} / (k + nu)
            g = f + r * q;  // g{k} = f{k} + 2/nu sin2(nu pi/2) q{k}
            h = p - k * g;  // h{k} = p{k} - k g{k}
            c *= neg_half_x_sq / k;
            S1 += delta1 = c * g;
            S2 += delta2 = c * h;
            if (fabs(delta1) < EPS + EPS * fabs(S1)) break;
        }
        if (k == MAXIT)
            nr::error("bessy series failed to converge");

        Ymu = -S1;
        Y1  = -S2 * two_ix;
        const RFLOAT Ymup = mu * ix * Ymu - Y1;
        Jmu = W / (Ymup - fmu * Ymu);
    } else {
        /** Evaluate CF2 by modified Lentz's method
         *
         * p + iq = (Jnu' - iYnu') / (Jnu - iYnu)
         *        = -1/2x + i + i/x * (1/2) ** 2 - nu ** 2 /
         *          2 * (x +  i)    + (3/2) ** 2 - nu ** 2 /
         *          2 * (x + 2i)
         *
         * The partial numerators are 1/4 - nu**2 + j**2 - j for j in 1, 2, 3...
         * The partial denominators are 2 * (x + ij) for j in 1, 2, 3...
         */
        Complex CF2 (-0.5 * ix, 1);
        RFLOAT  a = 0.25 - mumu;
        Complex b (2 * x, 2);
        Complex C = b + a * ix * Complex(0, 1) / CF2, D = 1 / b, delta;
        CF2 *= delta = C * D;
        int i = 1;
        for (; i != MAXIT; i++) {
            a      += 2 * i;
            b.imag += 2;

            C =     at_least(b + a / C, FPMIN);  // C =      b + a / C
            D = 1 / at_least(b + a * D, FPMIN);  // D = 1 / (b + a * D)

            CF2 *= delta = C * D;
            if (fabs(delta.real - 1) + fabs(delta.imag) < EPS) break;
        }
        if (i == MAXIT)
            nr::error("cf2 failed in besseljy");

        const RFLOAT gamma = (CF2.real - fmu) / CF2.imag;
        Jmu = copysign(sqrt(W / (CF2.imag * (gamma * gamma + 1))), Jl);
        Ymu = Jmu * gamma;
        const RFLOAT Ymup = Ymu * (CF2.real + CF2.imag / gamma);
        Y1 = mu * ix * Ymu - Ymup;
    }
    const RFLOAT scale = Jmu / Jl;
    // Scale original J{nu} and J'{nu}
    J  = Jl1 * scale;
    Jp = Jp1 * scale;
    // Upward recurrence of Y{nu}
    RFLOAT Y2;
    for (int i = 1; i <= nl; i++) {
        Y2  = (mu + i) * two_ix * Y1 - Ymu;  // Y{nu+1} = 2 nu/x Y{nu} - Y{nu-1}
        Ymu = Y1;
        Y1  = Y2;
    }
    Y  = Ymu;
    Yp = nu * ix * Ymu - Y1;  // Y'{nu} = nu/x Y{nu} - Y{nu+1}
}

//............................................................................
RFLOAT bessi0_5(RFLOAT x) {
    return x == 0 ? 0 : sqrt(2 / (PI * x)) * sinh(x);
}

RFLOAT bessi1_5(RFLOAT x) {
    return x == 0 ? 0 : sqrt(2 / (PI * x)) * (cosh(x) - sinh(x) / x);
}
RFLOAT bessi2(RFLOAT x) {
    return x == 0 ? 0 : bessi0(x) - 2 / x * bessi1(x);
}

RFLOAT bessi2_5(RFLOAT x) {
    return x == 0 ? 0 : bessi0_5(x) - 3 / x * bessi1_5(x);
}

RFLOAT bessi3(RFLOAT x) {
    return x == 0 ? 0 : bessi1(x) - 4 / x * bessi2(x);
}

RFLOAT bessi3_5(RFLOAT x) {
    return x == 0 ? 0 : bessi1_5(x) - 5 / x * bessi2_5(x);
}

RFLOAT bessi4(RFLOAT x) {
    return x == 0 ? 0 : bessi2(x) - 6 / x * bessi3(x);
}

namespace Bessel {

    RFLOAT J(RFLOAT x, RFLOAT nu) {
        RFLOAT J, Y, Jp, Yp;
        besseljy(x, nu, J, Y, Jp, Yp);
        return J;
    }

}

/* Special functions ------------------------------------------------------- */
RFLOAT gammln(RFLOAT xx) {
    static RFLOAT const cof[6] {
        76.18009173, -86.50532033, 24.01409822,
        -1.231739516, 0.120858003e-2, -0.536382e-5
    };

    RFLOAT x = xx - 1;
    RFLOAT tmp = x + 5.5;
    tmp -= (x + 0.5) * log(tmp);
    RFLOAT ser = 1;
    for (int j = 0; j <= 5; j++) {
        x += 1;
        ser += cof[j] / x;
    }
    return -tmp + log(2.50662827465 * ser);
}


RFLOAT betai(RFLOAT a, RFLOAT b, RFLOAT x) {
    if (x < 0 || x > 1)
        nr::error("Bad x in routine BETAI");
    const RFLOAT bt = x == 0.0 || x == 1 ? 0 :
        exp(gammln(a + b) - gammln(a) - gammln(b) + a * log(x) + b * log(1 - x));
    if (x < (a + 1) / (a + b + 2))
        return bt * betacf(a, b, x) / a;
    else
        return 1 - bt * betacf(b, a, 1 - x) / b;
}

RFLOAT betacf(RFLOAT a, RFLOAT b, RFLOAT x) {
    const int ITMAX = 100;
    const RFLOAT EPS = 3e-7;
    RFLOAT qab = a + b, qap = a + 1, qam = a - 1;
    RFLOAT az = 1, am = 1, bz = 1 - qab * x / qap, bm = 1;
    for (int m = 1; m <= ITMAX; m++) {
        const int two_m = m + m;
        RFLOAT d = m * (b - m) * x / ((qam + two_m) * (a + two_m));
        const RFLOAT ap = az + d * am;
        const RFLOAT bp = bz + d * bm;
        d = -(a + m) * (qab + m) * x / ((qap + two_m) * (a + two_m));
        const RFLOAT app = ap + d * az;
        const RFLOAT bpp = bp + d * bz;
        RFLOAT aold = az;
        am = ap  / bpp;
        bm = bp  / bpp;
        az = app / bpp;
        bz = 1;
        if (fabs(az - aold) < EPS * fabs(az)) return az;
    }
    nr::error("a or b too big, or ITMAX too small in BETACF");
    return 0;
}

/* Powell optimization ------------------------------------------------------------ */

void mnbrak(
    RFLOAT& ax, RFLOAT& bx, RFLOAT& cx,
    RFLOAT& fa, RFLOAT& fb, RFLOAT& fc, RFLOAT(*func)(RFLOAT *, void*),
    void *prm, int ncom, RFLOAT *pcom, RFLOAT *xicom
) {
    const RFLOAT GOLD = 1.618034;
    const RFLOAT GLIMIT = 100;
    const RFLOAT TINY = 1e-20;
    RFLOAT ulim, u, r, q, fu, dum;
    RFLOAT* xt = ask_vector<RFLOAT>(1, ncom);

    for (int j = 1; j <= ncom; j++)
        xt[j] = pcom[j] + ax * xicom[j];
    fa = func(xt, prm);
    for (int j = 1; j <= ncom; j++)
        xt[j] = pcom[j] + bx * xicom[j];
    fb = func(xt, prm);
    if (fb > fa) {
        std::swap(ax, bx);
        std::swap(fa, fb);
    }
    cx = bx + GOLD * (bx - ax);
    for (int j = 1; j <= ncom; j++)
        xt[j] = pcom[j] + cx * xicom[j];
    fc = func(xt, prm);
    while (fb > fc) {
        r = (bx - ax) * (fb - fc);
        q = (bx - cx) * (fb - fa);
        u = bx - ((bx - cx) * q - (bx - ax) * r) /
            (2 * copysign(std::max(fabs(q - r), TINY), q - r));
        ulim = bx + GLIMIT * (cx - bx);
        if ((bx - u) * (u - cx) > 0) {
            for (int j = 1; j <= ncom; j++)
                xt[j] = pcom[j] + u * xicom[j];
            fu = func(xt, prm);
            if (fu < fc) {
                ax = bx;
                bx = u;
                fa = fb;
                fb = fu;
                return;
            } else if (fu > fb) {
                cx = u;
                fc = fu;
                return;
            }
            u = cx + GOLD * (cx - bx);
            for (int j = 1; j <= ncom; j++)
                xt[j] = pcom[j] + u * xicom[j];
            fu = func(xt, prm);
        } else if ((cx - u) * (u - ulim) > 0) {
            for (int j = 1; j <= ncom; j++)
                xt[j] = pcom[j] + u * xicom[j];
            fu = func(xt, prm);
            if (fu < fc) {
                bx = cx; cx = u; u = cx + GOLD * (cx - bx);
                for (int j = 1; j <= ncom; j++)
                    xt[j] = pcom[j] + u * xicom[j];
                RFLOAT aux = func(xt, prm);
                fb = fc; fc = fu; fu = aux;
            }
        } else if ((u - ulim) * (ulim - cx) >= 0) {
            u = ulim;
            for (int j = 1; j <= ncom; j++)
                xt[j] = pcom[j] + u * xicom[j];
            fu = func(xt, prm);
        } else {
            u = cx + GOLD * (cx - bx);
            for (int j = 1; j <= ncom; j++)
                xt[j] = pcom[j] + u * xicom[j];
            fu = func(xt, prm);
        }
        ax = bx; bx = cx; cx = u;
        fa = fb; fb = fc; fc = fu;
    }
    free_vector<RFLOAT>(xt, 1, ncom);
}

RFLOAT brent(
    RFLOAT ax, RFLOAT bx, RFLOAT cx, RFLOAT(*func)(RFLOAT *,void*),
    void *prm, RFLOAT tol, RFLOAT *xmin,
    int ncom, RFLOAT *pcom, RFLOAT *xicom
) {
    const RFLOAT ZEPS = 1e-10;
    const RFLOAT CGOLD = 0.3819660;
    const int ITMAX = 100;
    RFLOAT d, etemp, fu, p, q, r, tol1, tol2, u, xm;
    RFLOAT e = 0;
    RFLOAT* xt = ask_vector<RFLOAT>(1, ncom);

    RFLOAT a = ax < cx ? ax : cx;
    RFLOAT b = ax > cx ? ax : cx;
    RFLOAT v, w, x, fv, fw, fx;
    v = w = x = bx;
    for (int j = 1; j <= ncom; j++)
        xt[j] = pcom[j] + x * xicom[j];
    fx = func(xt, prm);
    fv = fw = fx;
    for (int i = 1; i <= ITMAX; i++) {
        xm = 0.5 * (a + b);
        tol2 = 2 * (tol1 = tol * fabs(x) + ZEPS);
        if (fabs(x - xm) <= tol2 - 0.5 * (b - a)) {
            *xmin = x;
            free_vector<RFLOAT>(xt, 1, ncom);
            return fx;
        }
        if (fabs(e) > tol1) {
            r = (x - w) * (fx - fv);
            q = (x - v) * (fx - fw);
            p = (x - v) * q - (x - w) * r;
            q = 2 * (q - r);
            if (q > 0.0) p = -p;
            q = fabs(q);
            etemp = e;
            e = d;
            if (fabs(p) >= fabs(0.5 * q * etemp) || p <= q * (a - x) || p >= q * (b - x)) {
                d = CGOLD * (e = x >= xm ? a - x : b - x);
            } else {
                d = p / q;
                u = x + d;
                if (u - a < tol2 || b - u < tol2)
                    d = copysign(tol1, xm - x);
            }
        } else {
            d = CGOLD * (e = (x >= xm ? a - x : b - x));
        }
        u = (fabs(d) >= tol1 ? x + d : x + copysign(tol1, d));
        for (int j = 1; j <= ncom; j++)
            xt[j] = pcom[j] + u * xicom[j];
        fu = func(xt, prm);
        if (fu <= fx) {
            if (u >= x) {
                a = x;
            } else {
                b = x;
            }
             v =  w;  w =  x;  x =  u;
            fv = fw; fw = fx; fx = fu;
        } else {
            if (u < x) {
                a = u;
            } else {
                b = u;
            }
            if (fu <= fw || w == x) {
                v = w;
                w = u;
                fv = fw;
                fw = fu;
            } else if (fu <= fv || v == x || v == w) {
                v = u;
                fv = fu;
            }
        }
    }
    nr::error("Too many iterations in brent");
    *xmin = x;
    free_vector<RFLOAT>(xt, 1, ncom);
    return fx;
}

void linmin(
    RFLOAT *p, RFLOAT *xi, int n, RFLOAT &fret,
    RFLOAT(*func)(RFLOAT *, void*), void *prm
) {
    const RFLOAT TOL = 2e-4;
    const int ncom = n;
    RFLOAT* pcom  = ask_vector<RFLOAT>(1, n);
    RFLOAT* xicom = ask_vector<RFLOAT>(1, n);
    for (int j = 1; j <= n; j++) {
        pcom [j] = p [j];
        xicom[j] = xi[j];
    }
    RFLOAT ax = 0, xx = 1, bx = 2, fa, fx, fb;
    mnbrak(ax, xx, bx, fa, fx, fb, func, prm, ncom, pcom, xicom);
    RFLOAT xmin;
    fret = brent(ax, xx, bx, func, prm, TOL, &xmin, ncom, pcom, xicom);
    for (int j = 1; j <= n; j++) {
        p[j] += xi[j] *= xmin;
    }
    free_vector<RFLOAT>(xicom, 1, n);
    free_vector<RFLOAT>(pcom,  1, n);
}

void powell(
    RFLOAT *p, RFLOAT *xi, int n, RFLOAT ftol, int &iter,
    RFLOAT &fret, RFLOAT(*func)(RFLOAT *, void *), void *prm,
    bool show
) {
    const int ITMAX = 200;
    int i, ibig, j;
    RFLOAT t, fptt, fp, del;

    RFLOAT* pt  = ask_vector<RFLOAT>(1, n);
    RFLOAT* ptt = ask_vector<RFLOAT>(1, n);
    RFLOAT* xit = ask_vector<RFLOAT>(1, n);
    fret = func(p, prm);
    std::copy_n(pt + 1, n, p + 1);

    for (iter = 1;; iter++) {
        /* By coss ----- */
        if (show) {
            std::cout << iter << " (" << p[1];
            for (int co = 2; co <= n; co++)
                std::cout << "," << p[co];
            std::cout << ")--->" << fret << std::endl;
        }
        /* ------------- */

        fp = fret;
        ibig = 0;
        del = 0;
        for (i = 1; i <= n; i++) {
            bool any_nonzero = false; // CO
            for (j = 1; j <= n; j++) {
                xit[j] = xi[j * n + i];
                any_nonzero |= xit[j] != 0;
            }
            if (any_nonzero) {
                fptt = fret;
                linmin(p, xit, n, fret, func, prm);
                if (fabs(fptt - fret) > del) {
                    del = fabs(fptt - fret);
                    ibig = i;
                }
                /* By coss ----- */
                if (show) {
                    std::cout << "   (";
                    if (i == 1) { std::cout << "***"; }
                    std::cout << p[1];
                    for (int co = 2; co <= n; co++) {
                        std::cout << ",";
                        if (co == i) { std::cout << "***"; }
                        std::cout << p[co];
                    }
                    std::cout << ")--->" << fret << std::endl;
                }
                /* ------------- */
            }
        }
        if (2 * fabs(fp - fret) <= ftol * (fabs(fp) + fabs(fret))) {
            free_vector<RFLOAT>(xit, 1, n);
            free_vector<RFLOAT>(ptt, 1, n);
            free_vector<RFLOAT>(pt,  1, n);
            return;
        }
        if (iter == ITMAX)
            nr::error("Too many iterations in routine POWELL");
        for (j = 1; j <= n; j++) {
            ptt[j] = 2 * p[j] - pt[j];
            xit[j] =     p[j] - pt[j];
            pt [j] =     p[j];
        }
        fptt = func(ptt, prm);
        if (fptt < fp) {
            t = 2 * (fp - 2.0 * fret + fptt) * (fp - fret - del) * (fp - fret - del) - del * (fp - fptt) * (fp - fptt);
            if (t < 0) {
                linmin(p, xit, n, fret, func, prm);
                std::copy_n(xit + 1, n, xi + n + ibig);
            }
        }
    }
}

// Used in singular value decomposition (SVD) ---------------------------------
// https://en.wikipedia.org/wiki/Singular_value_decomposition
RFLOAT Pythag(RFLOAT a, RFLOAT b) {
    // Return the hypotenuse length of a right triangle with side lengths a and b.
    // Pythagoras' theorem: a^2 + b^2 = c^2 => c = sqrt(a^2 + b^2)
    // The difficulty with the naive implementation sqrt(a * a + b * b) is that x2 or y2 may overflow or underflow.
    // https://en.wikipedia.org/wiki/Pythagorean_addition#Implementation
    // The original code was copied from XmippCore's Bilib library:
    // https://github.com/I2PC/xmippCore/blob/devel/core/bilib/linearalgebra.cc
    // (btw, c++11 has std::hypot)

    // Don't waste time computing sqrt(0 * 0 + 0 * 0)
    if (a == 0 && b == 0) return 0;

    RFLOAT greater = fabs(a), lesser = fabs(b);
    if (greater < lesser) std::swap(greater, lesser);
    return greater * sqrt(1 + (lesser * lesser) / (greater * greater));
}

void svdcmp(RFLOAT *U, int Lines, int Columns, RFLOAT *W, RFLOAT *V) {
    // https://en.wikipedia.org/wiki/Singular_value_decomposition
    // Decompose a matrix into two square unitary matrice U and V,
    // and a diagonal matrix W.
    // #define SVDMAXITER 1000000
    const int MaxIterations = 1000000;

    RFLOAT* rv1 = ask_vector<RFLOAT>(0, Columns * Columns - 1);
    RFLOAT Norm = 0, Scale = 0, g = 0;
    long l = 0;
    for (long i = 0; i < Columns; i++) {
        l = i + 1;
        rv1[i] = Scale * g;
        RFLOAT s = 0;
        g = Scale = 0;
        if (i < Lines) {
            for (long k = i; k < Lines; k++) {
                Scale += fabs(U[k * Columns + i]);
            }
            if (Scale != 0) {
                for (long k = i; k < Lines; k++) {
                    U[k * Columns + i] /= Scale;
                    s += U[k * Columns + i] * U[k * Columns + i];
                }
                RFLOAT f = U[i * Columns + i];
                g = -copysign(sqrt(s), f);
                RFLOAT h = f * g - s;
                U[i * Columns + i] = f - g;
                for (long j = l; j < Columns; j++) {
                    RFLOAT s = 0;
                    for (long k = i; k < Lines; k++) {
                        s += U[k * Columns + i] * U[k * Columns + j];
                    }
                    f = s / h;
                    for (long k = i; k < Lines; k++) {
                        U[k * Columns + j] += f * U[k * Columns + i];
                    }
                }
                for (long k = i; k < Lines; k++) {
                    U[k * Columns + i] *= Scale;
                }
            }
        }
        W[i] = Scale * g;
        g = s = Scale = 0;
        if (i < Lines && i != Columns - 1) {
            for (long k = l; k < Columns; k++) {
                Scale += fabs(U[i * Columns + k]);
            }
            if (Scale != 0) {
                for (long k = l; k < Columns; k++) {
                    U[i * Columns + k] /= Scale;
                    s += U[i * Columns + k] * U[i * Columns + k];
                }
                RFLOAT f = U[i * Columns + l];
                g = -copysign(sqrt(s), f);
                RFLOAT h = f * g - s;
                U[i * Columns + l] = f - g;
                for (long k = l; k < Columns; k++) {
                    rv1[k] = U[i * Columns + k] / h;
                }
                for (long j = l; j < Lines; j++) {
                    RFLOAT s = 0;
                    for (long k = l; k < Columns; k++) {
                        s += U[j * Columns + k] * U[i * Columns + k];
                    }
                    for (long k = l; k < Columns; k++) {
                        U[j * Columns + k] += s * rv1[k];
                    }
                }
                for (long k = l; k < Columns; k++) {
                    U[i * Columns + k] *= Scale;
                }
            }
        }
        Norm = std::max(Norm, fabs(W[i]) + fabs(rv1[i]));
    }
    for (long i = Columns - 1; i >= 0; i--) {
        if (i < Columns - 1) {
            if (g != 0) {
                for (long j = l; j < Columns; j++) {
                    V[j * Columns + i] = U[i * Columns + j] / (U[i * Columns + l] * g);
                }
                for (long j = l; j < Columns; j++) {
                    RFLOAT s = 0;
                    for (long k = l; k < Columns; k++) {
                        s += U[i * Columns + k] * V[k * Columns + j];
                    }
                    if (s != 0)
                    for (long k = l; k < Columns; k++) {
                        V[k * Columns + j] += s * V[k * Columns + i];
                    }
                }
            }
            for (long j = l; j < Columns; j++) {
                V[i * Columns + j] = V[j * Columns + i] = 0;
            }
        }
        V[i * Columns + i] = 1;
        g = rv1[i];
        l = i;
    }
    for (long i = std::min(Lines, Columns) - 1; i >= 0; i--) {
        long l = i + 1;
        g = W[i];
        for (long j = l; j < Columns; j++) {
            U[i * Columns + j] = 0;
        }
        if (g != 0) {
            g = 1 / g;
            for (long j = l; j < Columns; j++) {
                RFLOAT s = 0;
                for (long k = l; k < Lines; k++) {
                    s += U[k * Columns + i] * U[k * Columns + j];
                }
                RFLOAT f = s * g / U[i * Columns + i];
                if (f != 0)
                for (long k = i; k < Lines; k++) {
                    U[k * Columns + j] += f * U[k * Columns + i];
                }
            }
            for (long j = i; j < Lines; j++) {
                U[j * Columns + i] *= g;
            }
        } else {
            for (long j = i; j < Lines; j++) {
                U[j * Columns + i] = 0;
            }
        }
        U[i * Columns + i] += 1;
    }
    for (long k = Columns - 1; k >= 0; k--) {
        for (long its = 1; its <= MaxIterations; its++) {
            bool flag = true;
            long nm = 0;
            long l = 0;
            for (l = k; l >= 0; l--) {
                nm = l - 1;
                if (fabs(rv1[l]) + Norm == Norm) {
                    flag = false;
                    break;
                }
                if (fabs(W[nm]) + Norm == Norm) {
                    break;
                }
            }
            if (flag) {
                RFLOAT c = 0, s = 1;
                for (long i = l; i <= k; i++) {
                    RFLOAT f = s * rv1[i];
                    rv1[i] *= c;
                    if (fabs(f) + Norm == Norm) {
                        break;
                    }
                    g = W[i];
                    RFLOAT h = Pythag(f, g);
                    W[i] = h;
                    h = 1 / h;
                    c = g * h;
                    s = -f * h;
                    for (long j = 0; j < Lines; j++) {
                        RFLOAT y = U[j * Columns + nm];
                        RFLOAT z = U[j * Columns + i];
                        U[j * Columns + nm] = y * c + z * s;
                        U[j * Columns + i]  = z * c - y * s;
                    }
                }
            }
            RFLOAT z = W[k];
            if (l == k) {
                if (z < 0) {
                    W[k] = -z;
                    for (long j = 0; j < Columns; j++) {
                        V[j * Columns + k] = -V[j * Columns + k];
                    }
                }
                break;
            }
            if (its == MaxIterations) {
                free_vector<RFLOAT>(rv1, 0, Columns * Columns - 1);
                return;
            }
            RFLOAT x = W[l];
            nm = k - 1;
            RFLOAT y = W[nm];
            g = rv1[nm];
            RFLOAT h = rv1[k];
            RFLOAT f = ((y - z) * (y + z) + (g - h) * (g + h)) / (2 * h * y);
            g = Pythag(f, 1);
            f = ((x - z) * (x + z) + h * (y / (f + copysign(fabs(g), f)) - h)) / x;
            RFLOAT c = 1, s = 1;
            for (long j = l; j <= nm; j++) {
                long i = j + 1;
                g = rv1[i];
                y = W[i];
                h = s * g;
                g = c * g;
                z = Pythag(f, h);
                rv1[j] = z;
                c = f / z;
                s = h / z;
                f = x * c + g * s;
                g = g * c - x * s;
                h = y * s;
                y *= c;
                for (long jj = 0; jj < Columns; jj++) {
                    RFLOAT x = V[jj * Columns + j];
                    RFLOAT z = V[jj * Columns + i];
                    V[jj * Columns + j] = x * c + z * s;
                    V[jj * Columns + i] = z * c - x * s;
                }
                z = Pythag(f, h);
                W[j] = z;
                if (z != 0) {
                    z = 1 / z;
                    c = f * z;
                    s = h * z;
                }
                f = c * g + s * y;
                x = c * y - s * g;
                for (long jj = 0; jj < Lines; jj++) {
                    RFLOAT y = U[jj * Columns + j];
                    RFLOAT z = U[jj * Columns + i];
                    U[jj * Columns + j] = y * c + z * s;
                    U[jj * Columns + i] = z * c - y * s;
                }
            }
            rv1[l] = 0;
            rv1[k] = f;
            W[k] = x;
        }
    }
    free_vector<RFLOAT>(rv1, 0, Columns * Columns - 1);
}

void svbksb(RFLOAT *u, RFLOAT *w, RFLOAT *v, int m, int n, RFLOAT *b, RFLOAT *x) {
    RFLOAT* tmp = ask_vector<RFLOAT>(1, n);
    for (int j = 1; j <= n; j++) {
        RFLOAT s = 0;
        if (w[j]) {
            for (int i = 1; i <= m; i++)
                s += u[i * n + j] * b[i];
            s /= w[j];
        }
        tmp[j] = s;
    }
    for (int j = 1; j <= n; j++) {
        RFLOAT s = 0;
        for (int jj = 1; jj <= n; jj++)
            s += v[j * n + jj] * tmp[jj];
        x[j] = s;
    }
    free_vector<RFLOAT>(tmp, 1, n);
}


/* Gamma function ---------------------------------------------------------- */

// Series expansion of Gamma function
void gser(RFLOAT& gamser, RFLOAT a, RFLOAT x, RFLOAT& gln) {
    const int ITMAX = 100;
    const RFLOAT EPS = 3.0e-7;

    gln = gammln(a);

    if (x < 0) nr::error("x less than 0 in routine gser");

    if (x == 0) {
        gamser = 0;
        return;
    }

    RFLOAT sum, delta, ap = a;
    sum = delta = 1 / a;
    for (int n = 1; n <= ITMAX; n++) {
        sum += delta *= x / ++ap;
        if (fabs(delta) < fabs(sum) * EPS) {
            gamser = sum * exp(-x + a * log(x) - gln);
            return;
        }
    }
    nr::error("a too large, ITMAX too small in routine gser");
    return;
}

// Continued-fraction evaluation of Gamma function
void gcf(RFLOAT& gammcf, RFLOAT a, RFLOAT x, RFLOAT& gln) {
    const int ITMAX = 100;
    const RFLOAT EPS = 3e-7;
    const RFLOAT FPMIN = 1e-30;

    gln = gammln(a);
    RFLOAT an;
    RFLOAT b = x + 1 - a;
    RFLOAT C = 1 / FPMIN, D = 1 / b, Delta;
    RFLOAT h = D;
    int i = 1;
    for (; i <= ITMAX; i++) {
        an = i * (a - i);
        b += 2;
        C =     at_least(b + an / C, FPMIN);
        D = 1 / at_least(an * D + b, FPMIN);
        h *= Delta = C * D;
        if (fabs(Delta - 1) < EPS) break;
    }
    if (i > ITMAX)
        nr::error("a too large, ITMAX too small in gcf");
    gammcf = exp(-x + a * log(x) - gln) * h;
}

RFLOAT gammp(RFLOAT a, RFLOAT x) {
    if (x < 0 || a <= 0) nr::error("Invalid arguments in routine gammp");
    RFLOAT gln;
    if (x < a + 1) {
        RFLOAT gamser;
        gser(gamser, a, x, gln);
        return gamser;
    } else {
        RFLOAT gammcf;
        gcf(gammcf, a, x, gln);
        return 1 - gammcf;
    }
}
