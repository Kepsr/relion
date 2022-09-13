#pragma once

#include <limits>
#include "src/multidim_array.h"

/// @name Statistical functions
//@{

template <typename T>
T max(const MultidimArray<T> &arr) {

    if (arr.size() <= 0) return static_cast<T>(0);
    return *std::max_element(arr.begin(), arr.end());

}

template <typename T>
T min(const MultidimArray<T> &arr) {

    if (arr.size() <= 0) return static_cast<T>(0);
    return *std::min_element(arr.begin(), arr.end());

}

// offset2indices;

// inline int f(int i, int j = 0, int k = 0, int l = 0) {
//     return g(i - v.xinit, j - v.yinit, k - v.zinit);
// }

template <typename T>
std::array<unsigned long int, 4> offset_to_direct_indices(
    const MultidimArray<T> &arr, long int offset
) {
    ldiv_t div;
    div = std::div(offset, arr.xdim * arr.ydim * arr.zdim);
    const unsigned long int l = div.quot;
    div = std::div(div.rem, arr.xdim * arr.ydim);
    const unsigned long int k = div.quot;
    div = std::div(div.rem, arr.xdim);
    const unsigned long int j = div.quot, i = div.rem;
    return {i, j, k, l};
}

template <typename T>
inline std::array<long int, 4> direct_indices_to_Xmipp_indices(
    const MultidimArray<T> &arr, long int i, long int j = 0, long int k = 0, long int l = 0
) {
    return {i + arr.xinit, j + arr.yinit, k + arr.zinit, l};
}

// Indices of the minimum element
template <typename T>
std::array<long int, 4> minIndex(const MultidimArray<T> &arr) {
    const auto it = std::min_element(arr.begin(), arr.end());
    const auto direct_indices = offset_to_direct_indices(it - arr.begin());
    return direct_indices_to_Xmipp_indices(direct_indices[0], direct_indices[1], direct_indices[2], direct_indices[3]);
}

// Indices of the maximum element
template <typename T>
std::array<long int, 4> maxIndex(const MultidimArray<T> &arr) {
    const auto it = std::max_element(arr.begin(), arr.end());
    const auto direct_indices = offset_to_direct_indices(it - arr.begin());
    return direct_indices_to_Xmipp_indices(direct_indices[0], direct_indices[1], direct_indices[2], direct_indices[3]);
}

// Minimum and maximum of the values in the array
template <typename T>
std::pair<T, T> minmax(const MultidimArray<T> &arr) {

    if (arr.size() <= 0) return {T(), T()};  // Value-initialise min and max

    const auto pair = std::minmax_element(arr.begin(), arr.end());
    return {*pair.first, *pair.second};
}

/** Average of the values in the array.
 *
 * Regardless of the type of the array, return an RFLOAT.
 */
template <typename T>
RFLOAT average(const MultidimArray<T> &arr) {

    if (arr.size() <= 0) return 0;
    // Arithmetic mean
    return static_cast<RFLOAT>(std::accumulate(arr.begin(), arr.end(), RFLOAT(0))) / arr.size();
}

/** Standard deviation of the values in the array.
 *
 * The returned value is always RFLOAT, regardless of the type of the array.
 *
 * stddev(N) = sqrt( sum for (int i = 0; i < N; i++) {(x[i] - mean(N)) ** 2} * 1 / N)
 */
template <typename T>
RFLOAT computeStddev(const MultidimArray<T> &arr) {
    const long int N = arr.size();
    if (N <= 1) return 0;

    RFLOAT avg = 0, stddev = 0;
    #ifdef RELION_SINGLE_PRECISION
        // Two passes through the data, as single-precision is not enough for a single pass
        // Also: averages of large arrays will give trouble: compute median first.
        RFLOAT median = 0.0;
        if (N > 1e6)
            median = arr.median();

        RFLOAT sumofdeviations = 0;
        for (const auto &x : arr) {
            sumofdeviations += static_cast<RFLOAT>(x) - median;
        }
        avg = median + sumofdeviations / N;

        RFLOAT sumofsquareddeviations = 0;
        for (auto y : arr) {
            RFLOAT x = static_cast<RFLOAT>(y);
            RFLOAT dev = x - avg;
            sumofsquareddeviations += dev * dev;
        }

        // Foreseeing numerical instabilities
        stddev = N <= 0 ? 0 : sqrt(static_cast<RFLOAT>(abs(sumofsquareddeviations / N - 1)));
    #else
        RFLOAT total = 0;
        RFLOAT sumofsquares = 0;

        for (auto y : arr) {
            RFLOAT x = static_cast<RFLOAT>(y);
            total += x;
            sumofsquares += x * x;
        }

        avg = total / N;

        if (N > 1) {
            // stddev(X) = sqrt(E[X ** 2] - E[X] ** 2)
            // RFLOAT var = (sumofsquares - avg * avg * N) / (N - 1);
            RFLOAT var = (sumofsquares / N - avg * avg) * N / (N - 1);
            // Unbiased sample variance
            // Foreseeing numerical instabilities
            stddev = sqrt(static_cast<RFLOAT>(abs(var)));
        } else {
            stddev = 0;
        }
    #endif
    return stddev;
}

template<typename T>
struct Stats {

    RFLOAT avg, stddev;
    T min, max;

    /** Print statistics
     *
     * No end of line character is written after this print out.
     *
     * @code
     * const auto stats = computeStats(arr);
     * std::cout << "Statistics: ";
     * stats.print(std::cout);
     * std::cout << std::endl;
     * @endcode
     */
    void print(std::ostream &out = std::cout) const {

        out.setf(std::ios::showpoint);
        const int old_prec = out.precision(7);

        out << " min= " << std::setw(9) << min
            << " max= " << std::setw(9) << max
            << " avg= " << std::setw(9) << avg
            << " dev= " << std::setw(9) << stddev
            << std::setprecision(old_prec);
    }

};

/** Compute statistics.
 *
 * Return the average, standard deviation, minimum and maximum.
 * A single pass through the entire array makes this faster
 * than separately computing average, standard deviation, minimum, and maximum.
 */
template <typename T>
Stats<T> computeStats(const MultidimArray<T> &arr) throw (const char*) {

    if (arr.size() <= 0)
        throw "Statistics cannot be computed for a dimensionless array!";

    double sumx = 0;
    double sumxx = 0;

    Stats<T> stats;
    stats.min = +std::numeric_limits<double>::max();
    stats.max = -std::numeric_limits<double>::max();

    // Make one pass through the array.
    for (auto xx : arr) {
        double x = static_cast<double>(xx);
        sumx  += x;
        sumxx += x * x;

            if (xx > stats.max) { stats.max = xx; } else
            if (xx < stats.min) { stats.min = xx; }
    }

    const long int N = arr.size();

    stats.avg = sumx / N;
    if (N > 1) {
        // Biased sample variance = E[X**2] - E[X]**2
        double var = (sumxx / N - stats.avg * stats.avg) * N / (N - 1);
        // Unbiased sample variance = biased sample variance * N / (N - 1)
        // Foreseeing numerical instabilities
        stats.stddev = sqrt(static_cast<RFLOAT>(abs(var)));
    } else {
        // Sample variance is undefined for N <= 1.
        stats.stddev = 0;
    }

    return stats;
}

/** Median
 *
 * @code
 * med = median(v1);
 * @endcode
 */
template <typename T>
RFLOAT median(const MultidimArray<T> &arr) throw (const char*) {

    if (arr.size <= 0) throw "Cannot take the median of an empty collection!";

    // Copy the array into a std::vector
    std::vector<T> copy (arr.begin(), arr.end());
    std::sort(copy.begin(), copy.end());

    const long int N = copy.size();
    if (N % 2 == 0) {
        return (RFLOAT) (copy[N / 2 - 1] + copy[N / 2]) / 2.0;
    } else {
        return copy[N / 2];
    }
}

/** Scale the array so that it spans the interval [minF,maxF]
 *
 * @code
 * const auto range = v.rangeAdjust(0, 1);
 * // The array now spans the interval [0,1] (provided range != 0)
 * @endcode
 */
template <typename T>
T rangeAdjust(MultidimArray<T> &arr, T minF, T maxF) {

    if (arr.size() <= 0) return;

    const auto interval = minmax(arr);
    const T range = interval.second - interval.first;
    const RFLOAT slope = range == 0 ? 0 :
        static_cast<RFLOAT>(maxF - minF) /
        static_cast<RFLOAT>(range);
    // If range == 0, then the array is flat.
    // No scaling can bring it to span the interval [minF, maxF].
    // We choose to just fill it with minF.
    for (auto &x : arr) {
        // a + b * x
        x = minF + static_cast<T>(slope * static_cast<RFLOAT>(x - interval.min));
    }
    return range;  // Let the caller know the range
}

// For use in rangeAdjust
template <typename T>
std::pair<T, T> minmax_inside_mask(const MultidimArray<T> &arr, const MultidimArray<int> &mask) {

    RFLOAT min = arr[0], max = arr[0];
    int *maskptr = mask.data;

    for (T *ptr = arr.begin(); ptr != arr.end(); ++ptr, ++maskptr) {
        if (*maskptr) {
            min = std::min(min, (RFLOAT) *ptr);
            max = std::max(max, (RFLOAT) *ptr);
        }
    }
    return {min, max};
}

/** Adjust the array so that within the mask it spans the interval [minF,maxF]
 *
 * A linear operation is performed on the values of the array
 * so that the values of the array are comprissed between the two values set.
 * The linear transformation is computed within the mask, but it is applied everywhere.
 *
 * @code
 * @code
 * const auto range = v.rangeAdjust(0, 1);
 * // The array now spans the interval [0,1] (provided range != 0)
 * @endcode
 */
template <typename T>
void rangeAdjust(MultidimArray<T> &arr, const MultidimArray<int> &mask, T minF, T maxF) {

    if (arr.size() <= 0) return;

    const auto interval = minmax_inside_mask(arr, mask);
    const T range = interval.second - interval.first;
    const RFLOAT slope = range == 0 ? 0 :
        static_cast<RFLOAT>(maxF - minF) /
        static_cast<RFLOAT>(range);
    // If range == 0, then the array is flat.
    // No scaling can bring it to span the interval [minF, maxF].
    // We choose to just fill it with minF.
    for (auto &x : arr) {
        // a + b * x
        x = minF + static_cast<T>(slope * static_cast<RFLOAT>(x - interval.min));
    }
    return range;  // Let the caller know the range
}


/** Adjust the range of the array to the range of another array in
    a least squares sense.
*
* A linear operation is performed on the values of the array so
* after it, the values of the self array are as similar as possible
* (L2 sense) to the values of the array shown as sample
*/

// As written this will only work for [T = RFLOAT]
// nevertheless since this is used is better
// to use T than RFLOAT or will create problem for int multidim arrays
template <typename T>
void rangeAdjust(
    MultidimArray<T> &arr,
    const MultidimArray<T> &target, const MultidimArray<int> *mask = nullptr
) {

    if (arr.size() <= 0) return;

    T *targetptr = target.data;
    int *maskptr = mask ? mask->data : nullptr;
    RFLOAT N = 0, sumx = 0, sumy = 0, sumxx = 0, sumxy = 0;
    for (T *ptr = arr.begin(); ptr != arr.end(); ++ptr, ++targetptr, ++maskptr) {
        if (!mask || *maskptr != 0) {
            N++;
            T x = *ptr; T y = *targetptr;
            sumx += x; sumxx += x * x;
            sumy += y; sumxy += x * y;
        }
    }
    RFLOAT slope = (N * sumxy - sumx * sumy) / (N * sumxx - sumx * sumx);
    RFLOAT intercept = sumy / N - slope * sumx / N;
    for (auto &x : arr) {
        // a + b * x
        x = static_cast<RFLOAT>(intercept + slope * static_cast<RFLOAT>(x));
    }
}

/** Adjust the average and stddev of the array to given values.
 *
 * A linear operation is performed on the values of the array,
 * after which the array's average shall be avgF
 * and its standard deviation shall be stddevF.
 * The array itself is modified.
 *
 * @code
 * v.statisticsAdjust(0, 1);
 * // Now the array has mean 0 and stddev 1.
 * @endcode
 */
template <typename T>
void statisticsAdjust(MultidimArray<T> &arr, RFLOAT avgF, RFLOAT stddevF) {

    if (arr.size() == 0) return;

    const auto stats = computeStats(arr);

    const RFLOAT a = stats.stddev == 0 ? 0 :
        static_cast<RFLOAT>(stddevF) / static_cast<RFLOAT>(stats.stddev);
    const RFLOAT b = static_cast<RFLOAT>(avgF) - a * static_cast<RFLOAT>(stats.avg);

    for (auto &x : arr) {
        x = static_cast<T>(a * static_cast<RFLOAT>(x) + b);
    }
}
