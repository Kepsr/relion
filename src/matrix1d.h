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
inline T& XX(const Matrix1D<T> &v) { return v[0]; }

/** Access to Y component
 * @code
 * YY(v) = 1;
 * val = YY(v);
 * @endcode
 */
template <typename T>
inline T& YY(const Matrix1D<T> &v) { return v[1]; }

/** Access to Z component
 * @code
 * ZZ(v) = 1;
 * val = ZZ(v);
 * @endcode
 */
template <typename T>
inline T& ZZ(const Matrix1D<T> &v) { return v[2]; }

/** Matrix1D class.*/
template<typename T>
class Matrix1D {

    public:

    /// The array itself
    T *vdata;

    /// Number of elements
    int vdim;

    /// false: column vector (default)
    /// true: row vector
    bool row;

    /// @name Constructors
    //@{
    /** Empty constructor
     *
     * The empty constructor creates a vector with no memory associated,
     * origin=0, size=0, no statistics, ... You can choose between a column
     * vector (by default), or a row one.
     *
     * @code
     * Matrix1D<RFLOAT> v1;
     * Matrix1D<RFLOAT> v1(true);
     * // both are examples of empty column vectors
     *
     * Matrix1D<int> v1(false);
     * // empty row vector
     * @endcode
     */
    Matrix1D(bool column = true): vdim(0), vdata(nullptr), row(!column) {}

    Matrix1D(std::initializer_list<T> list): Matrix1D((int) list.size()) {
        auto it = list.begin();
        for (int i = 0; i < size(); i++) {
            (*this)[i] = *(it + i);
        }
    }

    /** Dimension constructor
     *
     * The dimension constructor creates a vector with memory associated (but
     * not assigned to anything, could be full of garbage) origin=0, size=the
     * given one. You can choose between a column vector (by default), or a row
     * one.
     *
     * @code
     * Matrix1D<RFLOAT> v1(6);
     * Matrix1D<RFLOAT> v1(6, 'y');
     * // both are examples of column vectors of dimensions 6
     *
     * Matrix1D<int> v1('n');
     * // empty row vector
     * @endcode
     */
    Matrix1D(int dim, bool column = true): Matrix1D(column) {
        resize(dim);
    }

    /** Copy constructor
     */
    Matrix1D(const Matrix1D<T> &other): Matrix1D() {
        *this = other;
    }

    /** Conversion from one type to another
     * e.g. from Matrix1D<int> to Matrix1D<RFLOAT>
     * Each element is static_cast from U to T.
     */
    template<typename U>
    Matrix1D(const Matrix1D<U> &other): Matrix1D(other.size()) {
        for (int i = 0; i < other.size(); i++)
            (*this)[i] = static_cast<T>(other[i]);
    }

    /** Destructor
     */
    ~Matrix1D() { clear(); }

    /** Assignment
     */
    Matrix1D<T>& operator = (const Matrix1D<T> &other) {
        if (this != &other) {
            resize(other.size());
            for (int i = 0; i < size(); i++)
                (*this)[i] = other[i];
            row = other.row;
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
     * An exception is thrown if there is no memory.
     *
     * @code
     * V1.resize(3, 3, 2);
     * @endcode
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
        for (int i = 0; i < std::min(size(), new_vdim); i++) {
            new_vdata[i] = vdata[i];
        }
        // Zero-pad if necessary
        for (int i = size(); i < new_vdim; i++) {
            new_vdata[i] = 0;
        }

        // Delete old data
        clear();

        // Assign new data
        vdata = new_vdata;
        vdim  = new_vdim;
    }

    /** The size of the vector
     */
    inline int size() const { return vdim; }

    /** Is this a row vector?
     *
     * @code
     * if (v.isRow())
     *     std::cout << "v is a row vector\n";
     * @endcode
     */
    inline int isRow() const { return row; }

    /** Is this a column vector?
     *
     * @code
     * if (v.isCol())
     *     std::cout << "v is a column vector\n";
     * @endcode
     */
    inline int isCol() const { return !row; }

    /** Make the vector a row vector
     */
    void setRow() { row = true; }

    /** Make the vector a column vector
     */
    void setCol() { row = false; }
    //@}

    /// @name Initialization of Matrix1D values
    //@{
    /** Same value in all components.
     * @warning Used only once in symmetries.cpp
     *
     * The constant must be of a type compatible with the array type, ie,
     * you cannot assign a RFLOAT to an integer array without a casting.
     * It is not an error if the array is empty, then nothing is done.
     *
     * @code
     * v.initConstant(3.14);
     * @endcode
     */
    void initConstant(T x) {
        for (auto &y : *this) { y = x; }
    }

    /** Initialize to zeros with current size.
     *
     * All values are set to 0. The current size and origin are kept. 
     * If the array is empty, nothing is done.
     *
     * @code
     * v.initZeros();
     * @endcode
     */
    inline void initZeros() { memset(vdata, 0, size() * sizeof(T)); }

    /** Initialize to zeros with a given size.
     */
    void initZeros(int new_vdim) {
        resize(new_vdim);
        initZeros();
    }

    /** Initialize to zeros following a pattern.
     *
     * All values are set to 0, and the origin and size of the pattern are
     * adopted.
     *
     * @code
     * v2.initZeros(v1);
     * @endcode
     */
    template <typename T1>
    void initZeros(const Matrix1D<T1> &op) {
        resize(op);
        initZeros();
    }

    // Constructor
    static Matrix1D<T> zeros(int n) {
        Matrix1D<T> v (n);
        v.initZeros();
        return v;
    }
    //@}

    /// @name Matrix1D operators
    //@{
    Matrix1D<T> operator * (T x) const {
        Matrix1D<T> tmp(*this);
        for (int i = 0; i < size(); i++) { tmp[i] = (*this)[i] * x; }
        return tmp;
    }

    Matrix1D<T> operator / (T x) const {
        Matrix1D<T> tmp(*this);
        for (int i = 0; i < size(); i++) { tmp[i] = (*this)[i] / x; }
        return tmp;
    }

    Matrix1D<T> operator + (T x) const {
        Matrix1D<T> tmp(*this);
        for (int i = 0; i < size(); i++) { tmp[i] = (*this)[i] + x; }
        return tmp;
    }

    Matrix1D<T> operator - (T x) const {
        Matrix1D<T> tmp(*this);
        for (int i = 0; i < size(); i++) { tmp[i] = (*this)[i] - x; }
        return tmp;
    }

    friend Matrix1D<T> operator * (T x, const Matrix1D<T> &m) {
        Matrix1D<T> tmp(m);
        for (int i = 0; i < m.size(); i++) { tmp[i] = x * m[i]; }
        return tmp;
    }

    friend Matrix1D<T> operator / (T x, const Matrix1D<T> &m) {
        Matrix1D<T> tmp(m);
        for (int i = 0; i < m.size(); i++) { tmp[i] = x / m[i]; }
        return tmp;
    }

    friend Matrix1D<T> operator + (T x, const Matrix1D<T> &m) {
        Matrix1D<T> tmp(m);
        for (int i = 0; i < m.size(); i++) { tmp[i] = x + m[i]; }
        return tmp;
    }

    /** Vector summation
     *
     * @code
     * A += B;
     * @endcode
     */
    void operator += (const Matrix1D<T> &op1) const {
        if (size() != op1.size()) 
            REPORT_ERROR("Not same sizes in vector summation");

        for (int i = 0; i < size(); i++) { (*this)[i] += op1[i]; }
    }

    friend Matrix1D<T> operator - (T op1, const Matrix1D<T> &op2) {
        Matrix1D<T> tmp(op2);
        for (int i = 0; i < op2.size(); i++) { tmp[i] = op1 - op2[i]; }
        return tmp;
    }

    void operator *= (T x) { for (int i = 0; i < size(); i++) { (*this)[i] *= x; } }

    void operator /= (T x) { for (int i = 0; i < size(); i++) { (*this)[i] /= x; } }

    void operator += (T x) { for (int i = 0; i < size(); i++) { (*this)[i] += x; } }

    void operator -= (T x) { for (int i = 0; i < size(); i++) { (*this)[i] -= x; } }

    /** v3 = v1 * v2.
     */
    Matrix1D<T> operator * (const Matrix1D<T> &other) const {
        if (size() != other.size())
            REPORT_ERROR("Not same sizes in vector multiplication");

        Matrix1D<T> tmp(other);
        for (int i = 0; i < size(); i++) { tmp[i] = (*this)[i] * other[i]; }
        return tmp;
    }

    Matrix1D<T> operator / (const Matrix1D<T> &other) const {
        if (size() != other.size())
            REPORT_ERROR("Not same sizes in vector division");

        Matrix1D<T> tmp(other);
        for (int i = 0; i < size(); i++) { tmp[i] = (*this)[i] / other[i]; }
        return tmp;
    }

    Matrix1D<T> operator + (const Matrix1D<T> &other) const {
        if (size() != other.size())
            REPORT_ERROR("Not same sizes in vector summation");

        Matrix1D<T> tmp(other);
        for (int i = 0; i < size(); i++) { tmp[i] = (*this)[i] + other[i]; }
        return tmp;
    }

    Matrix1D<T> operator - (const Matrix1D<T> &other) const {
        if (size() != other.size())
            REPORT_ERROR("Not same sizes in vector subtraction");

        Matrix1D<T> tmp(other);
        for (int i = 0; i < size(); i++) { tmp[i] = (*this)[i] - other[i]; }
        return tmp;
    }

    void operator *= (const Matrix1D<T> &other) {
        if (size() != other.size())
            REPORT_ERROR("Not same sizes in vector multiplication");

        for (int i = 0; i < size(); i++) { (*this)[i] *= other[i]; }
    }

    void operator /= (const Matrix1D<T> &other) {
        if (size() != other.size())
            REPORT_ERROR("Not same sizes in vector division");

        for (int i = 0; i < size(); i++) { (*this)[i] /= other[i]; }
    }

    void operator += (const Matrix1D<T> &other) {
        if (size() != other.size())
            REPORT_ERROR("Not same sizes in vector summation");

        for (int i = 0; i < size(); i++) { (*this)[i] += other[i]; }
    }

    void operator -= (const Matrix1D<T> &other) {
        if (size() != other.size())
            REPORT_ERROR("Not same sizes in vector subtraction");

        for (int i = 0; i < size(); i++) { (*this)[i] -= other[i]; }
    }

    /** Negation
     */
    Matrix1D<T> operator - () const {
        Matrix1D<T> tmp(*this);
        for (int i = 0; i < size(); i++) { tmp[i] = -(*this)[i]; }
        return tmp;
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
    T& operator () (int i) const { return vdata[i]; }

    T& operator [] (int i) const { return vdata[i]; }
    //@}

    /// @name Utilities for Matrix1D
    //@{

    /** Produce a vector suitable for working with Numerical Recipes
     *
     * This function must be used only as a preparation for routines which need
     * that the first physical index is 1 and not 0 as it usually is in C. In
     * fact the vector provided for Numerical recipes is exactly this same one
     * but with the indices changed.
     *
     * This function is not ported to Python.
     */
    T* adaptForNumericalRecipes() const { return vdata - 1; }

    /** Index for the maximum element.
     *
     * This function returns the index of the maximum element of an matrix1d.
     * Returns -1 if the array is empty
     */
    int maxIndex() const {
        if (size() == 0) return -1;

        int imax = 0;
        T maxval = *begin();
        for (int i = 0; i < size(); i++) {
            if ((*this)[i] > maxval) { imax = i; }
        }
        return imax;
    }

    /** Index for the minimum element.
     *
     * Returns the index of the minimum element of an matrix1d.
     * Returns -1 if the array is empty
     */
    int minIndex() const {
        if (size() == 0) return -1;

        int imin = 0;
        T minval = *begin();
        for (int i = 0; i < size(); i++) {
            if ((*this)[i] < minval) { imin = i; }
        }
        return imin;
    }

    /** Algebraic transpose of vector
     */
    Matrix1D<T> transpose() const { 
        auto t (*this);
        t.row = !row;
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

        T *ptr;

        iterator(T *ptr): ptr(ptr) {}

        iterator& operator ++ () { ptr++; return *this; }

        iterator operator ++ (int) { iterator ret = *this; ++(*this); return ret; }

        bool operator == (iterator other) const { return ptr == other.ptr; }

        bool operator != (iterator other) const { return ptr != other.ptr; }

        T& operator * () { return *ptr; }

    };

    iterator begin() const { return iterator(&(*this)[0]); }
    iterator end() const { return iterator(&(*this)[size()]); }

    /** Sum of vector values.
     *
     * This function returns the sum of all internal values.
     *
     * @code
     * RFLOAT sum = m.sum();
     * @endcode
     */
    RFLOAT sum() const {
        return std::accumulate(begin(), end(), 0);
    }

    inline RFLOAT mean() const { return sum() / (RFLOAT) size(); }

    /** Sum of squared vector values.
     *
     * This function returns the sum of all internal values to the second power.
     */
    RFLOAT sum2() const {
        return std::accumulate(begin(), end(), 0, [] (T acc, T x) { return acc + x * x; });
    }

    /** Modulus (magnitude) of the vector
     *
     * This modulus is defined as the square root of the sum of the squared
     * components. Euclidean norm of the vector.
     *
     * @code
     * RFLOAT mod = v.modulus();
     * @endcode
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

    RFLOAT max() const {
        RFLOAT x = abs(*begin());  // Why abs?
        for (const auto &y : *this) { x = std::max(x, y); }
        return x;
    }

    // Normalise vector
    void normalise() {
        RFLOAT m = modulus();
        if (abs(m) > Xmipp::epsilon) {
            *this *= (T) (1.0 / m);
        } else {
            initZeros();  // Why?
        }
    }

    // Reverse data
    void reverse() {
        for (int i = 0; 2 * i <= size() - 1; i++) {
            std::swap((*this)[i], (*this)[size() - 1 - i]);
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
