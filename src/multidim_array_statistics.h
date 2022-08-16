#pragma once

#include "multidim_array.h"

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

/** 4D Indices for the minimum element.
 *
 * Return the index of the minimum element of an array.
 * array(i,j,k,l). Returns -1 if the array is empty
 */
template <typename T>
T minIndex(const MultidimArray<T> &arr, long int &imin, long int &jmin, long int &kmin, long int &lmin) {
    if (arr.xdim == 0) {
        imin = jmin = kmin = lmin = -1;
        return 0;
    }

    imin = arr.firstX();
    jmin = arr.firstY();
    kmin = arr.firstZ();
    lmin = 0;
    T minval = arr.elem(imin, jmin, kmin, lmin);

    FOR_ALL_NZYX_ELEMENTS_IN_MULTIDIMARRAY(arr) {
        if (arr.elem(i, j, k, l) > minval) {
            minval = arr.elem(i, j, k, l);
            imin = i;
            jmin = j;
            kmin = k;
            lmin = l;
        }
    }

    return minval;
}

/** 3D Indices for the minimum element.
 *
 * Call to the 4D function
 */
template <typename T>
T minIndex(const MultidimArray<T> &arr, long int &imin, long int &jmin, long int &kmin) {
    long int zero = 0;
    return minIndex(arr, imin, jmin, kmin, zero);
}

/** 2D Indices for the minimum element.
 *
 * Call to the 4D function
 */
template <typename T>
T minIndex(const MultidimArray<T> &arr, long int &imin, long int &jmin) {
    long int zero = 0;
    return minIndex(arr, imin, jmin, zero, zero);
}

/** 1D Indices for the minimum element.
 *
 * Call to the 4D function
 */
template <typename T>
T minIndex(const MultidimArray<T> &arr, long int &imin) {
    long int zero = 0;
    return minIndex(arr, imin, zero, zero, zero);
}

/** 4D Indices for the maximum element.
 *
 * Return the index of the maximum element of an array.
 * array(i,j,k,l). Returns -1 if the array is empty
 */
template <typename T>
T maxIndex(const MultidimArray<T> &arr, long int &imax, long int &jmax, long int &kmax, long int &lmax) {
    if (arr.xdim == 0) {
        lmax = kmax = imax = jmax = -1;
        return 0;
    }

    imax = arr.firstX();
    jmax = arr.firstY();
    kmax = arr.firstZ();
    lmax = 0;
    T maxval = arr.elem(imax, jmax, kmax, lmax);

    FOR_ALL_NZYX_ELEMENTS_IN_MULTIDIMARRAY(arr) {
        if (arr.elem(i, j, k, l) > maxval) {
            maxval = arr.elem(i, j, k, l);
            imax = i;
            jmax = j;
            kmax = k;
            lmax = l;
        }
    }

    return maxval;
}

/** 3D Indices for the maximum element.
 *
 * Call to the 4D function
 */
template <typename T>
T maxIndex(const MultidimArray<T> &arr, long int &imax, long int &jmax, long int &kmax) {
    long int _;
    return maxIndex(arr, imax, jmax, kmax, _);
}

/** 2D Indices for the maximum element.
 *
 * Call to the 4D function
 */
template <typename T>
T maxIndex(const MultidimArray<T> &arr, long int &imax, long int &jmax) {
    long int _;
    return maxIndex(arr, imax, jmax, _, _);
}

/** 1D Indices for the maximum element.
 *
 * Call to the 4D function
 */
template <typename T>
T maxIndex(const MultidimArray<T> &arr, long int &imax) {
    long int _;
    return maxIndex(arr, imax, _, _, _);
}

/** Minimum and maximum of the values in the array.
 *
 * As RFLOATs.
 */
template <typename T>
MinMax minmax(const MultidimArray<T> &arr) {

    RFLOAT min, max;

    if (arr.size() <= 0) return MinMax(min, max); // Uninitialised RFLOATs

    min = max = static_cast<RFLOAT>(arr[0]);
    for (const auto &x : arr) {
        if (x < min) { min = static_cast<RFLOAT>(x); } else
        if (x > max) { max = static_cast<RFLOAT>(x); }
    }
    return MinMax(min, max);
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
        int old_prec = out.precision(7);

        out << " min= "; out.width(9); out << min;
        out << " max= "; out.width(9); out << max;
        out << " avg= "; out.width(9); out << avg;
        out << " dev= "; out.width(9); out << stddev;

        out.precision(old_prec);
    }

};

/** Compute statistics.
 *
 * Return the average, standard deviation, minimum and maximum.
 * A single pass through the entire array makes this faster
 * than separately computing average, standard deviation, minimum, and maximum.
 */
template <typename T>
Stats<T> computeStats(const MultidimArray<T> &arr) {

    if (arr.size() <= 0)
        throw "Statistics cannot be computed for a dimensionless array!";

    double sumx = 0;
    double sumxx = 0;

    Stats<T> stats;
    stats.min = +std::numeric_limits<double>::max();
    stats.max = -std::numeric_limits<double>::max();

    // Make one pass through the array.
    for (const auto &xx : arr) {
        double x = static_cast<double>(xx);
        sumx  += x;
        sumxx += x * x;

            if (xx > stats.max) { stats.max = xx; } else
            if (xx < stats.min) { stats.min = xx; }
    }

    long int N = arr.size();

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
RFLOAT median(const MultidimArray<T> &arr) {

    if (arr.xdim == 0) throw "Cannot take the median of an empty collection!";

    // Copy the array
    auto copy = arr;

    // Sort indices
    copy.sort();

    long int N = arr.size();
    if (N % 2 == 0) {
        return (RFLOAT) (copy[N / 2 - 1] + copy[N / 2]) / 2.0;
    } else {
        return copy[N / 2];
    }
}

/** Adjust the range of the array to a given one.
 *
 * Scale the values of the array
 * so that they lie between the two values set.
 * Modify the array itself.
 *
 * @code
 * v.rangeAdjust(0, 1);
 * // The array is now ranging from 0 to 1
 * @endcode
 */
template <typename T>
void rangeAdjust(MultidimArray<T> &arr, T minF, T maxF) {

    if (arr.size() <= 0) return;

    const MinMax range = minmax(arr);

    // If range.min == range.max, the vector is a constant one,
    // so the only possible transformation is to a fixed minF
    RFLOAT slope = range.min == range.max ? 0 :
        static_cast<RFLOAT>(maxF - minF) /
        static_cast<RFLOAT>(range.max - range.min);

    for (auto &x : arr) {
        // a + b * x
        x = minF + static_cast<T>(slope * static_cast<RFLOAT>(x - range.min));
    }
}

// For use in rangeAdjust
template <typename T>
MinMax maskminmax(MultidimArray<T> &arr, MultidimArray<int> &mask) {
    RFLOAT min, max;
    int *maskptr = mask.data;

    min = max = arr[0];

    for (T *ptr = arr.begin(); ptr != arr.end(); ++ptr, ++maskptr) {
        if (*maskptr) {
            min = std::min(min, (RFLOAT) *ptr);
            max = std::max(max, (RFLOAT) *ptr);
        }
    }
    return MinMax(min, max);
}

/** Adjust the range of the array to a given one within a mask.
 *
 * A linear operation is performed on the values of the array
 * so that the values of the array are comprissed between the two values set.
 * The actual array is modified itself.
 * The linear transformation is computed within the mask, but it is applied everywhere.
 *
 * @code
 * v.rangeAdjust(0, 1, mask);
 * // The array is now ranging from 0 to 1
 * @endcode
 */
// This function must be explictly implemented outside
template <typename T>
void rangeAdjust(MultidimArray<T> &arr, MultidimArray<int> &mask, T minF, T maxF) {

    if (arr.size() <= 0) return;

    MinMax range = maskminmax(arr, mask);

    // If range.min == range.max, the vector is a constant one,
    // so the only possible transformation is to a fixed minF
    RFLOAT slope = range.min == range.max ? 0 :
        static_cast<RFLOAT>(maxF - minF) /
        static_cast<RFLOAT>(range.max - range.min);

    for (auto &x : arr) {
        // a + b * x
        x = minF + static_cast<T>(slope * static_cast<RFLOAT>(x - range.min));
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
// This function must be explictly implemented outside.
template <typename T>
void statisticsAdjust(MultidimArray<T> &arr, RFLOAT avgF, RFLOAT stddevF) {

    if (arr.size() == 0) return;

    Stats<T> stats = computeStats(arr);

    const RFLOAT a = stats.stddev == 0 ? 0 :
        static_cast<RFLOAT>(stddevF) / static_cast<RFLOAT>(stats.stddev);
    const RFLOAT b = static_cast<RFLOAT>(avgF) - a * static_cast<RFLOAT>(stats.avg);

    for (auto &x : arr) {
        x = static_cast<T>(a * static_cast<RFLOAT>(x) + b);
    }
}
