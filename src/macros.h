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

#ifndef MACROS_H
#define MACROS_H

#define RELION_SHORT_VERSION "3.1.3"
extern const char *g_RELION_VERSION;

#include <math.h>
#include <signal.h>
#include "src/pipeline_control.h"
#include "src/error.h"

#ifndef _CYGWIN
    #ifdef __APPLE__
        #include <limits.h>
    #else
        #include <values.h>
    #endif
#endif

#ifndef MINFLOAT
    #define MINFLOAT -1e30
#endif
#ifndef MAXFLOAT
    #define MAXFLOAT  1e30
#endif

#ifdef RELION_SINGLE_PRECISION
    #define RFLOAT float
    #define LARGE_NUMBER 99e36
    #define MY_MPI_DOUBLE  MPI_FLOAT
    #define MY_MPI_COMPLEX MPI_C_COMPLEX
#else
    #define RFLOAT double
    #define LARGE_NUMBER 99e99
    #define MY_MPI_DOUBLE  MPI_DOUBLE
    #define MY_MPI_COMPLEX MPI_C_DOUBLE_COMPLEX
#endif

#if defined CUDA and DEBUG_CUDA
    #define CRITICAL(string) raise(SIGSEGV);
#else
    #define CRITICAL(string) REPORT_ERROR(string);
#endif

//#define DEBUG
//#define DEBUG_CHECKSIZES

/// @defgroup Macros Macros
/// @ingroup DataLibrary
//@{
/// @name Constants
//@{

/** Pi
 * @ingroup MacrosConstants
 */
#ifndef PI
#define PI 3.14159265358979323846
#endif

/** Equal accuracy
 *
 * An epsilon for determining whether two values should be considered equal.
 * Value depends on whether RELION is operating in single or in double precision.
 * Required to correctly find symmetry subgroups.
 */
#ifdef RELION_SINGLE_PRECISION
#define XMIPP_EQUAL_ACCURACY 1e-4
#else
#define XMIPP_EQUAL_ACCURACY 1e-6
#endif
//@}

/// @name Numerical functions
//@{
/** Absolute value
 *
 * Valid for any numeric type (int, short, float, etc).
 *
 * @code
 * magnitude = ABS(x);
 * @endcode
 */
#ifndef ABS
#define ABS(x) ((x) >= 0 ? (x) : -(x))
#endif


/** signum
 *
 * The sign of a value.
 * Valid for any type T which supports 
 * - total ordering via operator <
 * - and construction from 0.
 * If n is positive, sgn(n) returns +1.
 * If n is negative, sgn(n) returns -1.
 * If n is 0,        sgn(n) returns  0.
 * (where 0, +1, -1 are of type T)
 *
 * @code
 * if (sgn(x) == -1) { std::cout << "x is negative" << std::endl; }
 * @endcode
 */
template <typename T>
inline int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}

template <typename T>
inline int sgn_nozero(T val) {
    return val >= T(0) ? +1 : -1;
}

#ifndef XMIPP_MIN
/** Minimum
 * @code
 * lesser = XMIPP_MIN(x, y);
 * @endcode
 */
#define XMIPP_MIN(x, y) ((x) <= (y) ? (x) : (y))
#endif

#ifndef XMIPP_MAX
/** Maximum
 * @code
 * greater = XMIPP_MAX(x, y);
 * @endcode
 */
#define XMIPP_MAX(x, y) ((x) <= (y) ? (y) : (x))
#endif

/** Round to next integer
 *
 * Valid for any kind of numbers (int, short, float, etc). The result is of type
 * integer.
 *
 * @code
 * a = ROUND(-0.8); // a = -1
 * a = ROUND(-0.2); // a = 0
 * a = ROUND(0.2); // a = 0
 * a = ROUND(0.8); // a = 1
 * @endcode
 */
#ifndef ROUND
#define ROUND(x) ((x) > 0 ? (int) ((x) + 0.5) : (int) ((x) - 0.5))
#endif

/** Round up
 *
 * Valid for any numeric type (int, short, float, etc). 
 * Return an int.
 *
 * @code
 * a = CEIL(-0.8); // a = 0
 * a = CEIL(-0.2); // a = 0
 * a = CEIL(0.2); // a = 1
 * a = CEIL(0.8); // a = 1
 * @endcode
 */
#define CEIL(x) ((x) == (int) (x) ? (int) (x) : ((x) > 0 ? (int) ((x) + 1) : (int) (x)))

/** Round down
 *
 * Valid for any numeric type (int, short, float, etc).
 * Return an int.
 *
 * @code
 * a = FLOOR(-0.8); // a = -1
 * a = FLOOR(-0.2); // a = -1
 * a = FLOOR(0.2); // a = 0
 * a = FLOOR(0.8); // a = 0
 * @endcode
 */
#define FLOOR(x) ((x) == (int) (x) ? (int) (x) : ((x) > 0 ? (int) (x) : (int) ((x) - 1)))

/** Return the fractional part of a value
 *
 * The fractional part of 3.7 is 0.7 and of -3.7 is -0.7.
 */
#define FRACTION(x) ((x) - (int) (x))

/** Clip in a saturation fashion
 *
 * CLIP is a macro which acts like a saturation curve, a value x is "clipped" to
 * a range defined by x0 and xF, for example the output values for the following
 * x and CLIP(x,-2,2) would be
 *
 * @code
 * x = ... -8 -7 -6 -5 -4 -3 -2 -1 0 1 2 3 4 5 6 7 8 ...
 * output = ... -2 -2 -2 -2 -2 -2 -2 -1 0 1 2 2 2 2 2 2 2 ...
 * @endcode
 */
#define CLIP(x, x0, xF) ((x) < (x0) ? (x0) : ((x) > (xF) ? (xF) : (x)))

/** Wrapping for integers
 *
 * intWRAP performs a wrapping in the integer set, when the cycle is finsihed it
 * begins again. For example, for intWRAP(x,-2,2) would be
 *
 * @code
 * x = ... -8 -7 -6 -5 -4 -3 -2 -1  0  1  2  3  4  5  6  7  8 ...
 * output = ...  2 -2 -1  0  1  2 -2 -1  0  1  2 -2 -1  0  1  2 -2 ...
 * @endcode
 */
#define intWRAP(x, x0, xF) (((x) >= (x0) && (x) <= (xF)) ? (x) : ((x) < (x0)) ? ((x) - (int)(((x) - (x0) + 1) / ((xF) - (x0) + 1) - 1) *  ((xF) - (x0) + 1)) : ((x) - (int)(((x) - (xF) - 1) / ((xF) - (x0) + 1) + 1) * ((xF) - (x0) + 1)))

/** Wrapping for real numbers
 *
 * realWRAP is used to keep a floating number between a range with a wrapping
 * fashion. For instance, it is used in trigonometry to say that an angle of
 * 5*PI is the same as PI, ie, to keep an angle in the range 0...2*PI
 *
 * @code
 * Corrected_angle = realWRAP(angle, 0, 2*PI);
 * @endcode
 */
#define realWRAP(x, x0, xF) ((x0) <= (x) && (x) <= (xF) ? \
    (x) : ((x) < (x0)) ? \
    ((x) - (int)(((x) - (x0)) / ((xF) - (x0)) - 1) * ((xF) - (x0))) : ((x) - (int)(((x) - (xF)) / ((xF) - (x0)) + 1) * ((xF) - (x0))))

/** Degrees to radians
 *
 * @code
 * angle_in_radians = DEG2RAD(ang_in_degrees);
 * @endcode
 */
#define DEG2RAD(theta) ((theta) * PI / 180.0)

/** Radians to degrees
 *
 * @code
 * angle_in_degrees = RAD2DEG(ang_in_radians);
 * @endcode
 */
#define RAD2DEG(theta) ((theta) * 180.0 / PI)

/** cos(theta), where theta is in degrees
 *
 * @code
 * if (COSD(90) == 0)
 *     std::cout << "This is in degrees!\n";
 * @endcode
 */
#define COSD(theta) cos(DEG2RAD(theta))

/** acos(x) in degrees
 *
 * @code
 * if (ACOSD(0.5) == 60)
 *     std::cout << "This is in degrees!\n";
 * @endcode
 */
#define ACOSD(x) RAD2DEG(acos(x))

/** sin(theta), where theta is in degrees
 *
 * @code
 * if (SIND(90) == 1)
 *     std::cout << "This is in degrees!\n";
 * @endcode
 */
#define SIND(theta) sin(DEG2RAD(theta))

/** acos(x) in degrees
 *
 * @code
 * if (ASIND(0.5) == 30.)
 *     std::cout << "This is in degrees!\n";
 * @endcode
 */
#define ASIND(theta) asin(theta) * 180.0 / PI

/** SINC function
 *
 * The sinc function is defined as sin(PI * x)/(PI * x).
 */
#define SINC(x) (((x) < 0.0001 && (x) > -0.0001) ? 1 : sin(PI * (x)) / (PI * (x)))

#if defined HAVE_SINCOS || defined DOXGEN

    /** Sincos function
     *
     *  Wrappper to make sincos(x, &sinval, &cosval) work for all compilers.
     */
    #define SINCOS(x, s, c) sincos(x, s, c)

    /** Sincosf function
     *
     *  Wrappper to make sincosf(x,&sinval,&cosval) work for all compilers.
     */
    #define SINCOSF(x, s, c) sincosf(x, s, c)

#elif defined HAVE___SINCOS
    // Use __sincos and __sincosf instead (primarily clang)

        #define SINCOS(x, s, c) __sincos(x, s, c)
        #define SINCOSF(x, s, c) __sincosf(x, s, c)

#else
    // Neither sincos or __sincos available, use raw functions.

        static void SINCOS(double x, double *s, double *c) {
            *s = sin(x); 
            *c = cos(x);
        }
        static void SINCOSF(float x, float *s, float *c) {
            *s = sinf(x); 
            *c = cosf(x);
        }

#endif

/** Returns next positive power_class of 2
 *
 * x is supposed to be positive, but it needn't be an integer.
 *
 * @code
 * next_power = NEXT_POWER_OF_2(1000); // next_power = 1024
 * @endcode
 */
#define NEXT_POWER_OF_2(x) (1 << ceil(log((RFLOAT) x) / log(2.0) - XMIPP_EQUAL_ACCURACY))

/** Linear interpolation
 *
 * From low (when a = 0) to high (when a = 1). 
 * The following value is returned (equal to (a * h) + ((1 - a) * l)
 */
#define LIN_INTERP(a, l, h) ((l) + ((h) - (l)) * (a))

/// @name Miscellaneous
//@{

/** Swap two values a and b (of the same size in memory).
 *
 * A new scope is entered, where a temporary variable is declared.
 * The data is exchanged and the scope is exited.
 */
#define SWAP(x, y) { \
    size_t stmp = sizeof(x); \
    unsigned char tmp[stmp == sizeof(y) ? (signed) stmp : -1]; \
    memcpy(tmp,  &(y), stmp); \
    memcpy(&(y), &(x), stmp); \
    memcpy(&(x),  tmp, stmp); \
}


namespace Xmipp {

    // The first index of an Xmipp volume/image/array of size 'size'.
    inline long int init(int size) {
        return -(long int) ((float) size / 2.0);
    }

    // The last index of an Xmipp volume/image/array of size 'size'.
    inline long int last(int size) {
        return size + init(size) - 1;
    }

}


static void PRINT_VERSION_INFO() {

    std::cout << "RELION version: " << g_RELION_VERSION << " "

    #if defined(DEBUG) || defined(DEBUG_CUDA)
    << "(debug-build) "
    #endif

    << std::endl << "Precision: BASE="

    #ifdef RELION_SINGLE_PRECISION
    << "single"
    #else
    << "double"
    #endif

    #if defined(CUDA) || defined(ALTCPU)

        #ifdef CUDA
        << ", CUDA-ACC="
        #endif

        #ifdef ALTCPU
        << ", VECTOR-ACC="
        #endif

        #ifdef ACC_DOUBLE_PRECISION
        << "double "
        #else
        << "single "
        #endif

    #endif

    << std::endl << std::endl;

}

//@}
//@}
#endif
