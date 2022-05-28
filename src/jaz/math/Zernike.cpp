#include "Zernike.h"
#include "src/error.h"
#include <sstream>
#include <cmath>
#include <vector>

using std::vector;

inline double safe_atan2(double y, double x) {
  return x == 0.0 && y == 0.0 ? 0.0 : atan2(y, x);
}

// Cache for coefficients of the radial polynomials
vector<vector<vector<double>>> R_coeffs (0);

long int factorial(int k) {
    // Alternatively: tgamma(k + 1)
    long int out = 1;
    for (int i = 2; i <= k; i++) { out *= i; }
    return out;
}

void resize_coefficients(int N) {
    vector<vector<vector<double>>> newCoeffs (N + 1);

    // Copy all of R_coeffs into the start of newCoeffs
    for (int n = 0; n < R_coeffs.size(); n++) {
        newCoeffs[n] = R_coeffs[n];
    }

    // Fill newCoeffs from where R_coeffs left off
    for (int n = R_coeffs.size(); n <= N; n++) {
        newCoeffs[n] = vector<vector<double>> (n + 1);

        for (int m = 0; m <= n; m++) {

            div_t diff_div2 = std::div(n - m, 2);
            if (diff_div2.rem == 1) continue;
            div_t sum_div2 = std::div(n + m, 2);

            newCoeffs[n][m] = vector<double> (diff_div2.quot + 1);

            for (int k = 0; k <= diff_div2.quot; k++) {
                int sign = 1 - 2 * (k % 2);  // Negative for even k, positive for odd k.
                newCoeffs[n][m][k] =
                      (double) (sign * factorial(n - k))
                    / (double) (factorial(k) * factorial(sum_div2.quot - k) * factorial(diff_div2.quot - k));
            }
        }
    }

    R_coeffs = newCoeffs;
}

Zernike::Z::Z(int m, int n): m{m}, n{n} {
    if (abs(m) > n) {
        REPORT_ERROR_STR("Requirement violated: abs(m) <= n (m = " << m << ", n = " << n << ")\n");
    }
}

// Value of a Zernike polynomial at radial distance rho and azimuthal angle phi
double Zernike::Z::operator () (double rho, double phi) {
    if (m >= 0) {
        return R{+m, n}(rho) * cos(+m * phi);  // Even
    } else {
        return R{-m, n}(rho) * sin(-m * phi);  // Odd
    }
}

double Zernike::Z::cart(double x, double y) {
    return (*this)(sqrt(x * x + y * y), safe_atan2(y, x));
}

Zernike::R::R(int m, int n): m{m}, n{n} {
    if (0 > m || m > n) {
        REPORT_ERROR_STR("Requirement violated: 0 <= m <= n (m = " << m << ", n = " << n << ")\n");
    }
}

// Value of a radial polynomial at radial distance rho
double Zernike::R::operator () (double rho) {

    if ((n - m) % 2 == 1) return 0.0;

    if (R_coeffs.size() <= n) {
        resize_coefficients(n);
    }

    double r = 0.0;
    for (int k = 0; k <= (n - m) / 2; k++) {
        r += R_coeffs[n][m][k] * pow(rho, n - 2 * k);
    }

    return r;
}

Zernike::Z Zernike::Z::fromEvenIndex(int i) {
    const int k = sqrt((double) i);

    int m = 2 * (i - k * k - k);
    int n = 2 * k;
    return { m, n };
}

int Zernike::numberOfEvenCoeffs(int n_max) {
    const int l = n_max / 2;
    return (l + 1) * (l + 1);
}

Zernike::Z Zernike::Z::fromOddIndex(int i) {
    const int k = (sqrt(1 + 4 * i) - 1.0) / 2.0;
    const int i0 = k * k + k;

    int n = 2 * k + 1;
    int m = 2 * (i - i0) - n;
    return { m, n };
}

int Zernike::numberOfOddCoeffs(int n_max) {
    const int l = (n_max - 1) / 2 + 1;
    return l * l + l;
}
