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

#ifndef MATRIX1D_H_
#define MATRIX1D_H_

#include <numeric>
#include <cstring>  // memcpy, memset
#include "src/funcs.h"

extern int bestPrecision(float F, int _width);
extern std::string floatToString(float F, int _width, int _prec);

template <typename T> class Matrix2D;

/** @defgroup Vectors Matrix1D Vectors
 * @ingroup DataLibrary
*/
//@{
/** @name Vectors speed up macros
 *
 * This macros are defined to allow high speed in critical parts of your program.
 * They shouldn't be used systematically
 * as usually there is no checking on the correctness of the operation you are performing.
 * Speed comes from three facts:
 * 1. They are macros and no function call is performed
 * (although most critical functions are inline functions).
 * 2. There is no checking on the correctness of the operation
 * (it could be wrong and you are not warned of it).
 * 3. Destination vectors are not returned,
 * saving time in the copy constructor and in the creation/destruction of temporary vectors.
 */
//@{

// Convention: { 0, 1, 2 } <-> { X, Y, Z }

template <typename T>
class Matrix1D;

/** Access to X component
 * @code
 * XX(v) = 1;
 * val = XX(v);
 * @endcode
 */
template <typename T>
inline T& XX(Matrix1D<T> &v) { return v[0]; }

/** Access to Y component
 * @code
 * YY(v) = 1;
 * val = YY(v);
 * @endcode
 */
template <typename T>
inline T& YY(Matrix1D<T> &v) { return v[1]; }

/** Access to Z component
 * @code
 * ZZ(v) = 1;
 * val = ZZ(v);
 * @endcode
 */
template <typename T>
inline T& ZZ(Matrix1D<T> &v) { return v[2]; }

/** Matrix1D class.*/
template<typename T>
class Matrix1D {

    /// The array itself
    T *vdata;

    /// Number of elements
    int vdim;

    /// Whether row vector or column vector
    enum class Mode: bool { row, column } mode;

    public:

    /// @name Constructors
    //@{

    /** Default / dimension constructor
     *
     * Construct a vector of the given size.
     * You can choose between a column (default) or a row vector.
     * The data will be zeroed out.
     * Default construction will produce an empty column vector.
     */
    Matrix1D(int dim = 0, Mode mode = Mode::column): vdim(0), vdata(nullptr), mode(mode) {
        resize(dim);
    }

    /** Copy constructor
     */
    Matrix1D(const Matrix1D<T> &other): Matrix1D(other.size(), other.mode) {
        std::copy(other.begin(), other.end(), begin());
    }

    /** Type-casting copy constructor
     */
    template<typename U>
    Matrix1D(const Matrix1D<U> &other): Matrix1D(other.size(), other.mode) {
        std::copy(other.begin(), other.end(), begin());
    }

    Matrix1D(std::initializer_list<T> list): Matrix1D(list.size()) {
        std::copy(list.begin(), list.end(), begin());
    }

    /** Destructor
     */
    ~Matrix1D() { clear(); }

    /** Assignment
     */
    template <typename U>
    Matrix1D<T>& operator = (const Matrix1D<U> &other) {
        if (this != &other) {
            resize(other.size());
            std::copy(other.begin(), other.end(), begin());
            mode = other.mode;
        }
        return *this;
    }
    //@}

    /** Clear.
     */
    void clear() {
        delete[] vdata;
        vdata = nullptr;
        vdim = 0;
    }

    ///@name Size and shape of Matrix1D
    //@{

    /** Resize to a given size
     *
     * Resize the array to the given size.
     * When the size decreases, the data get truncated.
     * When the size increases, the data get zero-padded.
     * Might throw a std::bad_alloc.
     */
    inline void resize(int new_vdim) {

        if (new_vdim == size()) return;

        if (new_vdim <= 0) {
            clear();
            return;
        }

        T *new_vdata;
        try {
            new_vdata = new T[new_vdim];
        } catch (std::bad_alloc &) {
            REPORT_ERROR("Allocate: No space left");
        }

        // Copy old data
        memcpy(new_vdata, vdata, sizeof(T) * std::min(vdim, new_vdim));
        // Zero-pad if necessary
        memset(new_vdata + vdim, 0, sizeof(T) * std::max(0, new_vdim - vdim));

        // Delete old data
        clear();

        // Assign new data
        vdata = new_vdata;
        vdim  = new_vdim;
    }

    /** The size of the vector
     */
    inline int size() const { return vdim; }

    inline T* data() { return vdata; }

    /** Is this a row vector?
     */
    inline int isRow() const { return mode == Mode::row; }

    /** Is this a column vector?
     */
    inline int isCol() const { return mode == Mode::column; }

    /** Make the vector a row vector
     */
    inline void setRow() { mode = Mode::row; }

    /** Make the vector a column vector
     */
    inline void setCol() { mode = Mode::column; }
    //@}

    /// @name Matrix1D operators
    //@{

    Matrix1D<T> operator + (T rhs) const {
        return Matrix1D<T>(*this) += rhs;
    }

    Matrix1D<T> operator - (T rhs) const {
        return Matrix1D<T>(*this) -= rhs;
    }

    Matrix1D<T> operator * (T rhs) const {
        return Matrix1D<T>(*this) *= rhs;
    }

    Matrix1D<T> operator / (T rhs) const {
        return Matrix1D<T>(*this) /= rhs;
    }

    friend Matrix1D<T> operator + (T lhs, Matrix1D<T> rhs) {
        for (auto& x: rhs) { x = lhs + x; }
        return rhs;
    }

    friend Matrix1D<T> operator - (T lhs, Matrix1D<T> rhs) {
        for (auto& x: rhs) { x = lhs - x; }
        return rhs;
    }

    friend Matrix1D<T> operator * (T lhs, Matrix1D<T> rhs) {
        for (auto& x: rhs) { x = lhs * x; }
        return rhs;
    }

    friend Matrix1D<T> operator / (T lhs, Matrix1D<T> rhs) {
        for (auto& x: rhs) { x = lhs / x; }
        return rhs;
    }

    Matrix1D<T>& operator += (T rhs) { for (auto& x: *this) { x += rhs; } return *this; }

    Matrix1D<T>& operator -= (T rhs) { for (auto& x: *this) { x -= rhs; } return *this; }

    Matrix1D<T>& operator *= (T rhs) { for (auto& x: *this) { x *= rhs; } return *this; }

    Matrix1D<T>& operator /= (T rhs) { for (auto& x: *this) { x /= rhs; } return *this; }

    Matrix1D<T> operator + (const Matrix1D<T> &other) const {
        return Matrix1D<T>(*this) += other;
    }

    Matrix1D<T> operator - (const Matrix1D<T> &other) const {
        return Matrix1D<T>(*this) -= other;
    }

    Matrix1D<T> operator * (const Matrix1D<T> &other) const {
        return Matrix1D<T>(*this) *= other;
    }

    Matrix1D<T> operator / (const Matrix1D<T> &other) const {
        return Matrix1D<T>(*this) /= other;
    }

    Matrix1D<T>& operator += (const Matrix1D<T> &other) {
        if (size() != other.size())
            REPORT_ERROR("Not same sizes in vector summation");
        for (int i = 0; i < size(); i++) { (*this)[i] += other[i]; }
        return *this;
    }

    Matrix1D<T>& operator -= (const Matrix1D<T> &other) {
        if (size() != other.size())
            REPORT_ERROR("Not same sizes in vector subtraction");
        for (int i = 0; i < size(); i++) { (*this)[i] -= other[i]; }
        return *this;
    }

    Matrix1D<T>& operator *= (const Matrix1D<T> &other) {
        if (size() != other.size())
            REPORT_ERROR("Not same sizes in vector multiplication");
        for (int i = 0; i < size(); i++) { (*this)[i] *= other[i]; }
        return *this;
    }

    Matrix1D<T>& operator /= (const Matrix1D<T> &other) {
        if (size() != other.size())
            REPORT_ERROR("Not same sizes in vector division");
        for (int i = 0; i < size(); i++) { (*this)[i] /= other[i]; }
        return *this;
    }

    /** Negation
     */
    Matrix1D<T> operator - () const {
        Matrix1D<T> copy (*this);
        for (auto& x: copy) { x = -x; }
        return copy;
    }

    /** Vector by matrix
     *
     * Algebraic vector by matrix multiplication.
     * This function is actually implemented in xmippMatrices2D.
     */
    Matrix1D<T> operator * (const Matrix2D<T> &M);  // Defined in matrix2D.h

    void operator *= (const Matrix2D<T> &M) {
        *this = *this * M;
    }

    /** Vector element access
     *
     * Returns the value of a vector logical iterator.
     * In our example we could access from v(-2) to v(2).
     * The elements can be used either by value or by reference.
     *
     * @code
     * v(-2) = 1;
     * val = v(-2);
     * @endcode
     */
          T& operator () (int i)       { return vdata[i]; }
    const T& operator () (int i) const { return vdata[i]; }

          T& operator [] (int i)       { return vdata[i]; }
    const T& operator [] (int i) const { return vdata[i]; }
    //@}

    /// @name Utilities for Matrix1D
    //@{

    /** Index for the maximum element.
     *
     * This function returns the index of the maximum element of an matrix1d.
     * Returns -1 if the array is empty
     */
    int maxIndex() const {
        if (size() == 0) return -1;
        return std::max_element(begin(), end()) - begin();
    }

    /** Index for the minimum element.
     *
     * Returns the index of the minimum element of an matrix1d.
     * Returns -1 if the array is empty
     */
    int minIndex() const {
        if (size() == 0) return -1;
        return std::min_element(begin(), end()) - begin();
    }

    /** Algebraic transpose of vector
     */
    Matrix1D<T> transpose() const {
        auto t (*this);
        t.mode = static_cast<Mode>(!static_cast<bool>(mode));
        return t;
    }

    struct iterator {

        /**
         * This struct allows us to concisely loop over vectors
         * with a range-based for.
         * @code
         * for (auto &x : v) {
         *     // do something with or to x
         * }
         * @endcode
         * The struct is meant to mimic a forward iterator.
         * It is essentially a wrapper around a pointer,
         * which can be made to traverse the heap-allocated memory
         * belonging to our Matrix1D<>.
         */

        using iterator_category = std::forward_iterator_tag;
        using difference_type   = std::ptrdiff_t;
        using value_type        = T;
        using pointer           = T*;
        using reference         = T&;

        T *ptr;

        iterator(T *ptr): ptr(ptr) {}

        iterator& operator ++ () { ptr++; return *this; }

        iterator operator ++ (int) { iterator copy (*this); ++ptr; return copy; }

        bool operator == (iterator other) const { return ptr == other.ptr; }

        bool operator != (iterator other) const { return ptr != other.ptr; }

        T& operator * () { return *ptr; }

        T* operator -> () { return ptr; }

    };

          iterator begin()       { return vdata; }
    const iterator begin() const { return vdata; }
          iterator end()       { return vdata + vdim; }
    const iterator end() const { return vdata + vdim; }

    /** Sum of vector values.
     */
    RFLOAT sum() const {
        return std::accumulate(begin(), end(), 0);
    }

    inline RFLOAT mean() const { return sum() / (RFLOAT) size(); }

    /** Sum of squared vector values.
     */
    RFLOAT sum2() const {
        return std::accumulate(begin(), end(), 0, [] (T acc, T x) { return acc + x * x; });
    }

    /** Modulus (magnitude) of the vector
     *
     * This modulus is defined as the square root of the sum of the squared
     * components. Euclidean norm of the vector.
     */
    RFLOAT modulus() const { return sqrt(sum2()); }

    /** Angle of the vector
     *
     * Supposing this vector is in R2 this function returns the angle of this
     * vector with X axis, ie, atan2(YY(v), XX(v))
     */
    RFLOAT angle() {
        return atan2((RFLOAT) (*this)[1], (RFLOAT) (*this)[0]);
    }

    // Normalise vector
    void normalise() {
        const RFLOAT m = modulus();
        if (abs(m) > Xmipp::epsilon) {
            *this *= (T) (1.0 / m);
        } else {
            std::fill(begin(), end(), 0);
        }
    }

    /** Compute numerical derivative
     *
     * https://en.wikipedia.org/wiki/Five-point_stencil
     * The numerical derivative is of the same size as the input vector.
     * However, the first two and the last two samples are set to 0,
     * because this method does not predict the derivative there.
     */
    Matrix1D<RFLOAT> numericalDerivative() {
        Matrix1D<RFLOAT> result = zeros(size());
        for (int i = 2; i <= this->xdim - 1 - 2; i++)
            result[i] = (- (*this)[i + 2] + 8 * (*this)[i + 1]
                         + (*this)[i - 2] - 8 * (*this)[i - 1]) / 12;
        return result;
    }

    /** Output to output stream.*/
    friend std::ostream& operator << (std::ostream &ostrm, const Matrix1D<T> &v) {
        if (v.size() == 0) { ostrm << "NULL Array"; }
        ostrm << '\n';

        const int prec = bestPrecision(v.max(), 10);

        for (const auto &x : v) {
            ostrm << floatToString((RFLOAT) x, 10, prec) << '\n';
        }
        return ostrm;
    }

    //@}
};

/**@name Vector Related functions
 * These functions are not methods of Matrix1D
 */

/** Creates vector in R2
 */
template <typename T>
inline Matrix1D<T> vectorR2(T x, T y) {
    return {x, y};
}

/** Creates vector in R3
 */
template <typename T>
inline Matrix1D<T> vectorR3(T x, T y, T z) {
    return {x, y, z};
}

// This function is only needed for single-precision compilation
#ifdef RELION_SINGLE_PRECISION
Matrix1D<float> vectorR3(double xx, double yy, double zz);
#endif

/** Dot product.
 * Given any two vectors in Rn (n-dimensional vector), this function returns the
 * dot product of both. If the vectors are not of the same size or shape then an
 * exception is thrown. The dot product is defined as the sum of the component
 * by component multiplication.
 *
 * For the R3 vectors (V1x,V1y,V1z), (V2x, V2y, V2z) the result is V1x*V2x +
 * V1y*V2y + V1z*V2z.
 *
 * @code
 * Matrix1D< RFLOAT > v1(1000);
 * v1.init_random(0, 10, "gaussian");
 * std::cout << "The power_class of this vector should be 100 and is " <<
 *     dotProduct(v1, v1) << std::endl;
 * @endcode
 */
template<typename T>
T dotProduct(const Matrix1D<T> &v1, const Matrix1D<T> &v2) {
    if (v1.size() != v2.size())
        REPORT_ERROR("Dot product: vectors of different size or shape");

    T accumulate = 0;
    for (int i = 0; i < v1.size(); i++) {
        accumulate += v1[i] * v2[i];
    }
    return accumulate;
}

/** Vector product in R3.
 * This function takes two R3 vectors and compute their vectorial product. For
 * two vectors (V1x,V1y,V1z), (V2x, V2y, V2z) the result is (V1y*V2z-V1z*v2y,
 * V1z*V2x-V1x*V2z, V1x*V2y-V1y*V2x). Pay attention that this operator is not
 * conmutative. An exception is thrown if the vectors are not of the same shape
 * or they don't belong to R3.
 *
 * @code
 * Matrix1D< T > X = vectorR3(1, 0, 0), Y = vector_R3(0, 1, 0);
 * std::cout << "X*Y=Z=" << crossProduct(X,Y).transpose() << std::endl;
 * @endcode
 */
template<typename T>
Matrix1D<T> crossProduct(const Matrix1D<T> &v1, const Matrix1D<T> &v2) {
    if (v1.size() != 3 || v2.size() != 3)
        REPORT_ERROR("Vector_product: vectors are not in R3");

    if (v1.isRow() != v2.isRow())
        REPORT_ERROR("Vector_product: vectors are of different shape");

    return {
        v1(1) * v2(2) - v1(2) * v2(1),
        v1(2) * v2(0) - v1(0) * v2(2),
        v1(0) * v2(1) - v1(1) * v2(0)
    };
}

/** Vector product in R3.
 * This function computes the vector product of two R3 vectors.
 * No check is performed, it is assumed that the output vector
 * is already resized
 *
 */
template<typename T>
void crossProduct(
    const Matrix1D<T> &v1, const Matrix1D<T> &v2, Matrix1D<T> &result
) {
    result[0] = v1[1] * v2[2] - v1[2] * v2[1];
    result[1] = v1[2] * v2[0] - v1[0] * v2[2];
    result[2] = v1[0] * v2[1] - v1[1] * v2[0];
}

/** Sort two vectors.
  * v1 and v2 must be of the same shape, if not an exception is thrown. After
  * calling this function all components in v1 are the minimum between the
  * corresponding components in v1 and v2, and all components in v2 are the
  * maximum.
  *
  * For instance, XX(v1)=MIN(XX(v1), XX(v2)), XX(v2)=MAX(XX(v1), XX(v2)). Notice
  * that both vectors are modified. This function is very useful for sorting two
  * corners. After calling it you can certainly perform a non-empty for (from
  * corner1 to corner2) loop.
  */
template<typename T>
void sortTwoVectors(Matrix1D<T> &v1, Matrix1D<T> &v2) {
    if (v1.size() != v2.size())
        REPORT_ERROR("sortTwoVectors: vectors are not of the same shape");

    // For two vectors v1 and v2, both of size N,
    // for all i in the interval [0, N),
    // v1[i] <= v2[i]
    for (int i = 0; i < v1.size(); i++) {
        if (v1[i] > v2[i]) std::swap(v1[i], v2[i]);
    }
}

//@}

#endif /* MATRIX1D_H_ */
