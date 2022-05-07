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

#include "src/funcs.h"
#include "src/filename.h"

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

/** For all elements in the array
 * This macro is used to generate loops for the vector in an easy manner. It
 * defines an internal index 'i' which ranges the vector using its mathematical
 * definition (ie, logical access).
 *
 * @code
 * FOR_ALL_ELEMENTS_IN_MATRIX1D(v)
 *     std::cout << v(i) << " ";
 * @endcode
 */
#define FOR_ALL_ELEMENTS_IN_MATRIX1D(v) for (int i = 0; i < v.vdim; i++)

/** X dimension of the matrix
 */
#define VEC_XSIZE(m) ((m).vdim)

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

/** Creates vector in R2
 * @code
 * v = VECTOR_R2(1.0, 2.0);
 * @endcode
 */
template <typename T>
Matrix1D<T> VECTOR_R2(T x, T y) {
    Matrix1D<T> v(2);
    XX(v) = x; 
    YY(v) = y;
    return v;
}

/** Creates vector in R3
 * @code
 * v = VECTOR_R2(1.0, 2.0, 1.0);
 * @endcode
 */
template <typename T>
Matrix1D<T> VECTOR_R3(T x, T y, T z) {
    Matrix1D<T> v(3);
    XX(v) = x; 
    YY(v) = y; 
    ZZ(v) = z;
    return v;
}

/** Adding two R2 vectors (a=b+c)
 * @code
 * MultidimArray< RFLOAT > a(2), b(2), c(2);
 * ...;
 * V2_PLUS_V2(a, b, c);
 * @endcode
 */
#define V2_PLUS_V2(a, b, c) { \
    XX(a) = XX(b) + XX(c); \
    YY(a) = YY(b) + YY(c); \
}

/** Substracting two R2 vectors (a=b-c)
 * @code
 * MultidimArray< RFLOAT > a(2), b(2), c(2);
 * ...;
 * V2_MINUS_V2(a, b, c);
 * @endcode
 */
#define V2_MINUS_V2(a, b, c) { \
        XX(a) = XX(b) - XX(c); \
        YY(a) = YY(b) - YY(c); }

/** Adding/substracting a constant to a R2 vector (a=b-k).
 * @code
 * MultidimArray< RFLOAT > a(2), b(2);
 * RFLOAT k;
 * ...;
 * V2_PLUS_CT(a, b, k);
 *
 * MultidimArray< RFLOAT > a(2), b(2);
 * RFLOAT k;
 * ...;
 * V2_PLUS_CT(a, b, -k);
 * @endcode
 */
#define V2_PLUS_CT(a, b, k) { \
        XX(a) = XX(b) + (k); \
        YY(a) = YY(b) + (k); }

/** Multiplying/dividing by a constant a R2 vector (a=b*k)
 * @code
 * MultidimArray< RFLOAT > a(2), b(2);
 * RFLOAT k;
 * ...;
 * V2_BY_CT(a, b, k);
 *
 * MultidimArray< RFLOAT > a(2), b(2);
 * RFLOAT k;
 * ...;
 * V2_BY_CT(a, b, 1/k);
 * @endcode
 */
#define V2_BY_CT(a, b, k) { \
        XX(a) = XX(b) * (k); \
        YY(a) = YY(b) * (k); }

/** Adding two R3 vectors (a = b + c)
 * @code
 * MultidimArray<RFLOAT> a(3), b(3), c(3);
 * ...;
 * V3_PLUS_V3(a, b, c);
 * @endcode
 */
#define V3_PLUS_V3(a, b, c) { \
    XX(a) = XX(b) + XX(c); \
    YY(a) = YY(b) + YY(c); \
    ZZ(a) = ZZ(b) + ZZ(c); \
}

/** Substracting two R3 vectors (a=b-c)
 * @code
 * MultidimArray< RFLOAT > a(3), b(3), c(3);
 * ...;
 * V3_MINUS_V3(a, b, c);
 * @endcode
 */
#define V3_MINUS_V3(a, b, c) { \
    XX(a) = XX(b) - XX(c); \
    YY(a) = YY(b) - YY(c); \
    ZZ(a) = ZZ(b) - ZZ(c); \
}

/** Adding/substracting a constant to a R3 vector (a=b-k)
 * @code
 * MultidimArray< RFLOAT > a(3), b(3);
 * RFLOAT k;
 * ...;
 * V3_PLUS_CT(a, b, k);
 *
 * MultidimArray< RFLOAT > a(3), b(3);
 * RFLOAT k;
 * ...;
 * V3_PLUS_CT(a, b, -k);
 * @endcode
 */
#define V3_PLUS_CT(a, b, c) { \
        XX(a) = XX(b) + (c); \
        YY(a) = YY(b) + (c); \
        ZZ(a) = ZZ(b) + (c); }

/** Multiplying/dividing by a constant a R3 vector (a=b*k)
 * @code
 * MultidimArray< RFLOAT > a(3), b(3);
 * RFLOAT k;
 * ...;
 * V3_BY_CT(a, b, k);
 *
 * MultidimArray< RFLOAT > a(3), b(3);
 * RFLOAT k;
 * ...;
 * V3_BY_CT(a, b, 1/k);
 * @endcode
 */
#define V3_BY_CT(a, b, c) { \
        XX(a) = XX(b) * (c); \
        YY(a) = YY(b) * (c); \
        ZZ(a) = ZZ(b) * (c); }

/** Direct access to vector element
 */
#define VEC_ELEM(v,i) ((v).vdata[(i)])
//@}

/** Matrix1D class.*/
template<typename T>
class Matrix1D {

    public:

    /// The array itself
    T *vdata;

    /// Number of elements
    int vdim;

    /// <0=column vector (default), 1=row vector
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
    Matrix1D(bool column = true) {
        coreInit();
        row = !column;
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
    Matrix1D(int dim, bool column = true) {
        coreInit();
        row = !column;
        resize(dim);
    }

    /** Copy constructor
     *
     * The created vector is a perfect copy of the input vector but with a
     * different memory assignment.
     *
     * @code
     * Matrix1D<RFLOAT> v2(v1);
     * @endcode
     */
    Matrix1D(const Matrix1D<T> &v) {
        coreInit();
        *this = v;
    }

    /** Destructor.
     */
    ~Matrix1D() { coreDeallocate(); }

    /** Assignment.
     *
     * You can build as complex assignment expressions as you like. Multiple
     * assignment is allowed.
     *
     * @code
     * v1 = v2 + v3;
     * v1 = v2 = v3;
     * @endcode
     */
    Matrix1D<T>& operator = (const Matrix1D<T> &other) {
        if (this != &other) {
            resize(other);
            for (int i = 0; i < size(); i++) { (*this)[i] = other[i]; }
            row = other.row;
        }
        return *this;
    }
    //@}

    /// @name Core memory operations for Matrix1D
    //@{
    /** Clear.
     */
    void clear() {
        coreDeallocate();
        coreInit();
    }

    /** Core init.
     * Initialize everything to 0
     */
    void coreInit() {
        vdim = 0;
        row = false;
        vdata = NULL;
    }

    /** Core allocate.
     */
    inline void coreAllocate(int _vdim) {
        if (_vdim <= 0) {
            clear();
            return;
        }

        vdim = _vdim;
        vdata = new T[vdim];
        if (!vdata) REPORT_ERROR("Allocate: No space left");
    }

    /** Core deallocate.
     * Free all vdata.
     */
    inline void coreDeallocate() {
        delete[] vdata;
        vdata = NULL;
    }
    //@}

    ///@name Size and shape of Matrix1D
    //@{
    /** Resize to a given size
     *
     * This function resize the actual array to the given size. The origin is
     * not modified. If the actual array is larger than the pattern then the
     * values outside the new size are lost, if it is smaller then 0's are
     * added. An exception is thrown if there is no memory.
     *
     * @code
     * V1.resize(3, 3, 2);
     * @endcode
     */
    inline void resize(int Xdim) {

        if (Xdim == size()) return;

        if (Xdim <= 0) {
            clear();
            return;
        }

        T *new_vdata;
        try {
            new_vdata = new T[Xdim];
        } catch (std::bad_alloc &) {
            REPORT_ERROR("Allocate: No space left");
        }

        // Copy vdata into new_vdata
        for (int j = 0; j < Xdim; j++) {
            new_vdata[j] = j >= size() ? 0 : vdata[j];
            // Fill with 0 if out of bounds
        }

        // deallocate old vector
        coreDeallocate();

        // assign *this vector to the newly created
        vdata = new_vdata;
        vdim = Xdim;
    }

    /** Resize according to a pattern.
     *
     * This function resize the actual array to the same size
     * as the input pattern. If the actual array is larger than the pattern
     * then the trailing values are lost, if it is smaller then 0's are
     * added at the end
     *
     * @code
     * v2.resize(v1);
     * // v2 has got now the same structure as v1
     * @endcode
     */
    template<typename T1>
    inline void resize(const Matrix1D<T1> &other) { resize(other.size()); }

    /** Same shape.
     *
     * Returns true if this object has got the same shape (origin and size)
     * as the argument
     */
    template <typename T1>
    inline bool sameShape(const Matrix1D<T1> &other) const {
        return size() == other.size();
    }

    /** Returns the size of this vector
     *
     * @code
     * int nn = a.size();
     * @endcode
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

    /** Forces the vector to be a row vector
     *
     * @code
     * v.setRow();
     * @endcode
     */
    void setRow() { row = true; }

    /** Forces the vector to be a column vector
     *
     * @code
     * v.setCol();
     * @endcode
     */
    void setCol() { row = false; }
    //@}

    /// @name Initialization of Matrix1D values
    //@{
    /** Same value in all components.
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
        for (int i = 0; i < size(); i++) { (*this)[i] = x; }
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
    void initZeros(int Xdim) {
        resize(Xdim);
        memset(vdata, 0, size() * sizeof(T));
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
        memset(vdata, 0, size() * sizeof(T));
    }

    // Constructor
    static Matrix1D<T> zeros(int n) {
        Matrix1D<T> v(n);
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

    /** Vector subtraction
     *
     * @code
     * A -= B;
     * @endcode
     */
    void operator -= (const Matrix1D<T> &op1) const {
        if (size() != op1.size())
            REPORT_ERROR("Not same sizes in vector subtraction");

        for (int i = 0; i < size(); i++) { (*this)[i] -= op1[i]; }
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
     *
     * It is used to build arithmetic expressions.
     * You can make a minus of anything as long as it is correct semantically.
     *
     * @code
     * v1 = -v2;
     * v1 = -v2.transpose();
     * @endcode
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
     * Returns the value of a vector logical position.
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

    /** Kill an array produced for Numerical Recipes.
     *
     * Nothing needs to be done in fact.
     *
     * This function is not ported to Python.
     */
    void killAdaptationForNumericalRecipes(T* m) const {
        // Do nothing
    }

    /** CEILING
     *
     * Round each array element up.
     */
    void selfCEIL() {
        for (int i = 0; i < size(); i++) { (*this)[i] = ceil((*this)[i]); }
    }

    /** FLOOR
     *
     * Round each array element down.
     */
    void selfFLOOR() {
        for (int i = 0; i < size(); i++) { (*this)[i] = floor((*this)[i]); }
    }

    /** ROUND
     *
     * Round each array element.
     */
    void selfROUND() {
        for (int i = 0; i < size(); i++) { (*this)[i] = round((*this)[i]); }
    }

    /** Index for the maximum element.
     *
     * This function returns the index of the maximum element of an matrix1d.
     * Returns -1 if the array is empty
     */
    int maxIndex() const {
        if (size() == 0) return -1;

        int imax = 0;
        T maxval = (*this)[0];
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
        T minval = (*this)[0];
        for (int i = 0; i < size(); i++) {
            if ((*this)[i] < minval) { imin = i; }
        }
        return imin;
    }

    /** Algebraic transpose of vector
     *
     * You can use the transpose in as complex expressions as you like.
     * The original vector is not changed.
     *
     * @code
     * v2 = v1.transpose();
     * @endcode
     */
    Matrix1D<T> transpose() const {
        Matrix1D<T> tmp(*this);
        tmp.selfTranspose();
        return tmp;
    }

    /** Algebraic transpose of vector
     *
     * The same as before but the result is stored in this same object.
     */
    void selfTranspose() { row = !row; }

    /** Sum of vector values.
     *
     * This function returns the sum of all internal values.
     *
     * @code
     * RFLOAT sum = m.sum();
     * @endcode
     */
    RFLOAT sum() const {
        RFLOAT sum = 0;
        for (int i = 0; i < size(); i++) { sum += (*this)[i]; }
        return sum;
    }

    inline RFLOAT mean() const { return sum() / (RFLOAT) size(); }

    /** Sum of squared vector values.
     *
     * This function returns the sum of all internal values to the second
     * power_class.
     *
     * @code
     * RFLOAT sum2 = m.sum2();
     * @endcode
     */
    RFLOAT sum2() const {
        RFLOAT sum = 0;
        for (int i = 0; i < size(); i++) { sum += (*this)[i] * (*this)[i]; }
        return sum;
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
        RFLOAT x = abs((*this)[0]);  // Why abs?
        for (int i = 0; i < size(); i++) {
            x = std::max(x, (*this)[i]);
        }
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
     * The numerical derivative is of the same size as the input vector.
     * However, the first two and the last two samples are set to 0,
     * because the numerical method is not able to correctly estimate the
     * derivative there.
     */
    Matrix1D<RFLOAT> numericalDerivative() {
        Matrix1D<RFLOAT> result = zeros(size());
        const RFLOAT one_twelfth = 1.0 / 12.0;
        for (int i = (this->xinit) + 2; i <= (this->xinit + this->xdim - 1) - 2; i++)  // Wrong type
            result[i] = one_twelfth * (
                - (*this)[i + 2] + 8 * (*this)[i + 1] 
                + (*this)[i + 2] - 8 * (*this)[i - 1]  // BUG? Should [i + 2] be [i - 2]?
            );
        return result;
    }

    /** Output to output stream.*/
    friend std::ostream& operator << (std::ostream &ostrm, const Matrix1D<T> &v) {
        if (v.size() == 0) { ostrm << "NULL Array"; }
        ostrm << '\n';

        int prec = bestPrecision(v.max(), 10);

        for (int i = 0; i < v.size(); i++) {
            ostrm << floatToString((RFLOAT) v[i], 10, prec) << '\n';
        }
        return ostrm;
    }

    //@}
};

/**@name Vector Related functions
 * These functions are not methods of Matrix1D
 */

/** Creates vector in R2.
 * After this function the vector is (x,y) in R2.
 *
 * @code
 * Matrix1D< RFLOAT > v = vectorR2(1, 2);
 * @endcode
 */
Matrix1D<RFLOAT> vectorR2(RFLOAT x, RFLOAT y);

/** Creates vector in R3.
 * After this function the vector is (x,y,z) in R3.
 *
 * @code
 * Matrix1D< RFLOAT > v = vectorR2(1, 2, 1);
 * @endcode
 */
Matrix1D<RFLOAT> vectorR3(RFLOAT x, RFLOAT y, RFLOAT z);

// This function is only needed for single-precision compilation
#ifdef RELION_SINGLE_PRECISION
Matrix1D<float> vectorR3(double xx, double yy, double zz);
#endif

/** Creates an integer vector in Z3.
 */
Matrix1D<int> vectorR3(int x, int y, int z);

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
    if (!v1.sameShape(v2))
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
 * std::cout << "X*Y=Z=" << vectorProduct(X,Y).transpose() << std::endl;
 * @endcode
 */
template<typename T>
Matrix1D<T> vectorProduct(const Matrix1D<T> &v1, const Matrix1D<T> &v2) {
    if (v1.size() != 3 || v2.size() != 3)
        REPORT_ERROR("Vector_product: vectors are not in R3");

    if (v1.isRow() != v2.isRow())
        REPORT_ERROR("Vector_product: vectors are of different shape");

    Matrix1D<T> result(3);
    result(0) = v1(1) * v2(2) - v1(2) * v2(1);
    result(1) = v1(2) * v2(0) - v1(0) * v2(2);
    result(2) = v1(0) * v2(1) - v1(1) * v2(0);
    return result;
}

/** Vector product in R3.
 * This function computes the vector product of two R3 vectors.
 * No check is performed, it is assumed that the output vector
 * is already resized
 *
 */
template<typename T>
void vectorProduct(
    const Matrix1D<T> &v1, const Matrix1D<T> &v2, Matrix1D<T> &result
) {
    result(0) = v1(1) * v2(2) - v1(2) * v2(1);
    result(1) = v1(2) * v2(0) - v1(0) * v2(2);
    result(2) = v1(0) * v2(1) - v1(1) * v2(0);
 }

template <typename T>
struct LesserGreater {

    T lesser, greater;

    LesserGreater(T x, T y) {
        if (x > y) {
            greater = x; lesser = y;
        } else {
            greater = y; lesser = x;
        }
    }

};

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
    if (!v1.sameShape(v2))
        REPORT_ERROR("sortTwoVectors: vectors are not of the same shape");

    for (int i = 0; i < v1.size(); i++) {
        LesserGreater<T> lg (v1[i], v2[i]);
        v1[i] = lg.lesser;
        v2[i] = lg.greater;
    }
 }

/** Conversion from one type to another.
  * If we have an integer array and we need a RFLOAT one, we can use this
  * function. The conversion is done through a type casting of each element
  * If n >= 0, only the nth volumes will be converted, otherwise all NSIZE volumes
  */
template<typename T1, typename T2>
void typeCast(const Matrix1D<T1> &v1, Matrix1D<T2> &v2) {
    if (v1.size() == 0) {
        v2.clear();
        return;
    }

    v2.resize(v1.size());
    for (int i = 0; i < v1.size(); i++) { v2[i] = static_cast<T2>(v1[i]); }
}
//@}

#endif /* MATRIX1D_H_ */
