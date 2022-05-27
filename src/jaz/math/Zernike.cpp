#include "Zernike.h"
#include "src/error.h"
#include <sstream>
#include <cmath>

std::vector<std::vector<std::vector<double>>> Zernike::R_coeffs (0);

double Zernike::Z(int m, int n, double rho, double phi)
{
    if (m >= 0)
    {
        return R(m, n, rho) * cos(m * phi);
    }
    else
    {
        return R(-m, n, rho) * sin(-m * phi);
    }
}

double Zernike::Z_cart(int m, int n, double x, double y) {
    return Z(m, n, sqrt(x * x + y * y), x == 0 && y == 0 ? 0.0 : atan2(y,x));
}

/** Radial polynomials
 */
double Zernike::R(int m, int n, double rho)
{
    if (m > n)
    {
        REPORT_ERROR_STR("Zernike::R: illegal argument: m = " << m << ", n = " << n << ".\n");
    }

    if ((n - m) % 2 == 1) return 0.0;

    if (R_coeffs.size() <= n)
    {
        prepCoeffs(n);
    }

    double r = 0.0;
    for (int k = 0; k <= (n - m) / 2; k++) {
        r += R_coeffs[n][m][k] * pow(rho, n - 2*k);
    }

    return r;
}

Zernike::MN Zernike::evenIndexToMN(int i)
{
    const int k = (int)sqrt((double)i);

    int m = 2*(i - k*k - k);
    int n = 2*k;
    return { m, n };
}

int Zernike::numberOfEvenCoeffs(int n_max)
{
    const int l = n_max / 2;
    return (l + 1) * (l + 1);
}

Zernike::MN Zernike::oddIndexToMN(int i)
{
    const int k = (sqrt(1 + 4 * i) - 1.0) / 2.0;
    const int i0 = k*k + k;

    int n = 2 * k + 1;
    int m = 2 * (i - i0) - n;
    return { m, n };
}

int Zernike::numberOfOddCoeffs(int n_max)
{
    const int l = (n_max - 1) / 2 + 1;
    return l * l + l;
}

long int factorial(int k)
{
    // @TODO: replace by tgamma(k+1) once C++11 becomes available
    long int out = 1;
    for (int i = 2; i <= k; i++) { out *= i; }
    return out;
}

void Zernike::prepCoeffs(int N)
{
    std::vector<std::vector<std::vector<double>>> newCoeffs(N + 1);

    // Copy all of R_coeffs into the start of newCoeffs
    for (int n = 0; n < R_coeffs.size(); n++)
    {
        newCoeffs[n] = R_coeffs[n];
    }

    // Fill newCoeffs from where R_coeffs left off
    for (int n = R_coeffs.size(); n <= n; n++)
    {
        newCoeffs[n] = std::vector<std::vector<double>>(n + 1);

        for (int m = 0; m <= n; m++)
        {
            if ((n - m) % 2 == 1) continue;

            newCoeffs[n][m] = std::vector<double>((n - m) / 2 + 1);

            for (int k = 0; k <= (n-m)/2; k++)
            {
                newCoeffs[n][m][k] =
                      (double) ((1 - 2 * (k % 2)) * factorial(n - k))
                    / (double) (factorial(k) * factorial((n + m) / 2 - k) * factorial((n - m) / 2 - k));
            }
        }
    }

    R_coeffs = newCoeffs;
}
