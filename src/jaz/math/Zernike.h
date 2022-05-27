#ifndef ZERNIKE_H
#define ZERNIKE_H

namespace Zernike {

    // Zernike polynomial
    struct Z {

        const int m, n;  // Where abs(m) <= n

        Z(int m, int n);

        static Z fromOddIndex(int i);
        static Z fromEvenIndex(int i);

        double operator () (double rho, double phi);
        double cart(double x, double y);

    };

    // Radial polynomial
    struct R {

        const int m, n;  // Where 0 <= m <= n

        R(int m, int n);

        double operator () (double rho);

    };

    int numberOfEvenCoeffs(int n_max);

    int numberOfOddCoeffs(int n_max);

};

#endif
