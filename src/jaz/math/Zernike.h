#ifndef ZERNIKE_H
#define ZERNIKE_H

#include <vector>

class Zernike {

    public:

    struct MN { int m, n; };

    static double Z(int m, int n, double rho, double phi);
    static double Z_cart(int m, int n, double x, double y);
    static double R(int m, int n, double rho);

    static MN evenIndexToMN(int i);
    static int numberOfEvenCoeffs(int n_max);

    static MN oddIndexToMN(int i);
    static int numberOfOddCoeffs(int n_max);

    private:

    // Cache for coefficients of the radial polynomials
    static std::vector<std::vector<std::vector<double>>> R_coeffs;

    static void prepCoeffs(int n);

};

#endif
