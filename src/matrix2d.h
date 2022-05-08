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

#ifndef MATRIX2D_H_
#define MATRIX2D_H_

#include <string.h>
#include <iomanip>
#include "src/matrix1d.h"
#include "src/filename.h"

/** @defgroup Matrices Matrix2D Matrices
 * @ingroup DataLibrary
 */
//@{
/** @name Matrices speed up macros */
//@{

/** For all elements in the array
 *
 * This macro simplifies looping over a matrix's values.
 * It ranges over the matrix's indices and binds them to 'i' and 'j'.
 *
 * @code
 * FOR_ALL_ELEMENTS_IN_MATRIX2D(m) {
 *     std::cout << m(i, j) << " ";
 * }
 * @endcode
 */
#define FOR_ALL_ELEMENTS_IN_MATRIX2D(m) \
    for (int i = 0; i < (m).mdimy; i++) \
        for (int j = 0; j < (m).mdimx; j++)

// Matrix2D class
template<typename T>
class Matrix2D {

    public:

    // The array itself
    T *mdata;

    // Number of elements in X
    int mdimx;

    // Number of elements in Y
    int mdimy;

    // Total number of elements
    int mdim;

    /// @name Constructors
    /// @{
    // Empty constructor
    Matrix2D() {
        coreInit();
    }

    // Dimension constructor
    Matrix2D(int Ydim, int Xdim) {
        coreInit();
        resize(Ydim, Xdim);
    }

    // Copy constructor
    Matrix2D(const Matrix2D<T> &v) {
        coreInit();
        *this = v;
    }

    // Destructor
    ~Matrix2D() { coreDeallocate(); }

    Matrix2D<T>& operator = (const Matrix2D<T> &op1) {
        if (&op1 != this) {
            resize(op1);
            memcpy(mdata, op1.mdata, op1.mdim * sizeof(T));
        }
        return *this;
    }
    //@}

    /// @name Core memory operations for Matrix2D
    //@{
    // Clear
    void clear() {
        coreDeallocate();
        coreInit();
    }

    void coreInit() {
        // Initialise everything to 0
        mdimx = mdimy = mdim = 0;
        mdata = NULL;
    }

    void coreAllocate(int _mdimy, int _mdimx) {
        if (_mdimy <= 0 || _mdimx <= 0) {
            clear();
            return;
        }

        mdimx = _mdimx;
        mdimy = _mdimy;
        mdim = _mdimx * _mdimy;
        mdata = new T[mdim];
        if (!mdata) REPORT_ERROR("coreAllocate: No space left");
    }

    void coreDeallocate() {
        delete[] mdata;
        mdata = NULL;
    }
    //@}

    /// @name Size and shape of Matrix2D
    //@{
    // Resize to a given size
    void resize(int Ydim, int Xdim) {

        if (Xdim == mdimx && Ydim == mdimy) return;

        if (Xdim <= 0 || Ydim <= 0) {
            clear();
            return;
        }

        T *new_mdata;

        try {
            new_mdata = new T[Xdim * Ydim];
        } catch (std::bad_alloc &) {
            REPORT_ERROR("Allocate: No space left");
        }

        for (int i = 0; i < Ydim; i++)
        for (int j = 0; j < Xdim; j++)
        // Copy needed elements
        // Fill with 0 if necessary
        new_mdata[i * Xdim + j] = i >= mdimy || j >= mdimx ? 0 : mdata[i * mdimx + j];

        // Deallocate old vector
        coreDeallocate();

        // Assign *this vector to the newly created
        mdata = new_mdata;
        mdimx = Xdim;
        mdimy = Ydim;
        mdim = Xdim * Ydim;
    }

    /** Resize according to a pattern.
    *
    * This function resize the actual array to the same size and origin
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
    void resize(const Matrix2D<T1> &v) { resize(v.mdimy, v.mdimx); }

    // Extract submatrix and assign to this object
    void submatrix(int i0, int j0, int iF, int jF) {
        if (i0 < 0 || j0 < 0 || iF >= mdimy || jF >= mdimx)
            REPORT_ERROR("Submatrix indices out of bounds");
        Matrix2D<T> result(iF - i0 + 1, jF - j0 + 1);

        FOR_ALL_ELEMENTS_IN_MATRIX2D(result)
        result.at(i, j) = at(i + i0, j + j0);

        *this = result;
    }

    /** Same shape.
    *
    * Returns true if this object has got the same shape (origin and size)
    * than the argument
    */
    template <typename T1>
    bool sameShape(const Matrix2D<T1> &op) const {
        return mdimx == op.mdimx && mdimy == op.mdimy;
    }

    // X dimension
    inline int Xdim() const { return mdimx; }

    // Y dimension
    inline int Ydim() const { return mdimy; }
    //@}

    /// @name Initialise Matrix2D values
    //@{
    /** Same value in all components.
    *
    * The constant must be of a type compatible with the array type, ie,
    * you cannot  assign a RFLOAT to an integer array without a casting.
    * It is not an error if the array is empty, then nothing is done.
    *
    * @code
    * v.initConstant(3.14);
    * @endcode
    */
    void initConstant(T val) {
        for (int j = 0; j < mdim; j++) { mdata[j] = val; }
    }

    /** Initialise to zeros with current size.
    *
    * All values are set to 0. The current size and origin are kept. It is not
    * an error if the array is empty, then nothing is done.
    *
    * @code
    * v.initZeros();
    * @endcode
    */
    void initZeros() {
        memset(mdata, 0, mdimx * mdimy * sizeof(T));
    }

    // Initialise to zeros with a given size
    void initZeros(int Ydim, int Xdim) {
        resize(Ydim, Xdim);
        memset(mdata, 0, mdimx * mdimy * sizeof(T));
    }

    static Matrix2D<T> zeros(int Ydim, int Xdim) {
        Matrix2D<T> m(Ydim, Xdim);
        m.initZeros();
        return m;
    }

    /** Initialise to zeros following a pattern.
    *
    * All values are set to 0, and the origin and size of the pattern are
    * adopted.
    *
    * @code
    * v2.initZeros(v1);
    * @endcode
    */
    template <typename T1>
    void initZeros(const Matrix2D<T1> &op) {
        resize(op);
        memset(mdata, 0, mdimx * mdimy * sizeof(T));
    }

    /** 2D Identity matrix of current size
    *
    * If actually the matrix is not square then an identity matrix is
    * generated of size (Xdim x Xdim).
    *
    * @code
    * m.initIdentity();
    * @endcode
    */
    void initIdentity() {
        initIdentity(mdimx);
    }

    /** 2D Identity matrix of a given size
    *
    * A (dim x dim) identity matrix is generated.
    *
    * @code
    * m.initIdentity(3);
    * @endcode
    */
    void initIdentity(int dim) {
        initZeros(dim, dim);
        for (int i = 0; i < dim; i++) { at(i, i) = 1; }
    }
    //@}

    /// @name Operators for Matrix2D
    //@{

    /** Matrix element access
     *
     * @code
     * m.at(0, 0) = 1;
     * std::cout << m.at(0, 0) << std::endl;
     * @endcode
     */
    inline T& at(int i, int j) const { return mdata[i * mdimx + j]; }

    inline T& operator () (int i, int j) { return at(i, j); }

    inline const T& operator () (int i, int j) const { return at(i, j); }

    // v3 = v1 * k
    Matrix2D<T> operator * (T op1) const {
        Matrix2D<T> tmp(*this);
        for (int i = 0; i < mdim; i++) { tmp.mdata[i] = mdata[i] * op1; }
        return tmp;
    }

    // v3 = v1 / k
    Matrix2D<T> operator / (T op1) const {
        Matrix2D<T> tmp(*this);
        for (int i = 0; i < mdim; i++) { tmp.mdata[i] = mdata[i] / op1; }
        return tmp;
    }

    // v3 = k * v2
    friend Matrix2D<T> operator * (T op1, const Matrix2D<T> &op2) {
        Matrix2D<T> tmp(op2);
        for (int i = 0; i < op2.mdim; i++) { tmp.mdata[i] = op1 * op2.mdata[i]; }
        return tmp;
    }

    // v3 *= k
    void operator *= (T op1) {
        for (int i = 0; i < mdim; i++) { mdata[i] *= op1; }
    }

    // v3 /= k
    void operator /=  (T op1) {
        for (int i = 0; i < mdim; i++) { mdata[i] /= op1; }
    }

    /** Matrix * vector multiplication
    *
    * @code
    * v2 = A*v1;
    * @endcode
    */
    Matrix1D<T> operator * (const Matrix1D<T> &op1) const {
        Matrix1D<T> result;

        if (mdimx != op1.size()) {
            std::cerr << " mdimx= " << mdimx << " opp1.size()= " << op1.size() << std::endl;
            REPORT_ERROR("Incompatible sizes in matrix by vector");
        }

        if (!op1.isCol())
            REPORT_ERROR("Vector is not a column");

        result.initZeros(mdimy);

        for (int i = 0; i < mdimy; i++)
        for (int j = 0; j < op1.size(); j++)
        result[i] += (*this)(i, j) * op1[j];

        result.setCol();
        return result;
    }

    /** Matrix by Matrix multiplication
    *
    * @code
    * C = A*B;
    * @endcode
    */
    Matrix2D<T> operator * (const Matrix2D<T> &op1) const {
        if (mdimx != op1.mdimy)
            REPORT_ERROR("Not compatible sizes in matrix multiplication");

        Matrix2D<T> result = Matrix2D<T>::zeros(mdimy, op1.mdimx);
        for (int i = 0; i < mdimy; i++)
        for (int j = 0; j < op1.mdimx; j++)
        for (int k = 0; k < mdimx; k++)
        result(i, j) += (*this)(i, k) * op1(k, j);
        return result;
    }

    /** Matrix addition
    *
    * @code
    * C = A + B;
    * @endcode
    */
    Matrix2D<T> operator + (const Matrix2D<T> &op1) const {
        Matrix2D<T> result;
        if (mdimx != op1.mdimx || mdimy != op1.mdimy)
            REPORT_ERROR("operator+: Not same sizes in matrix addition");

        result.initZeros(mdimy, mdimx);
        for (int i = 0; i < mdimy; i++)
        for (int j = 0; j < mdimx; j++)
        result(i, j) = (*this)(i, j) + op1(i, j);

        return result;
    }

    /** In-place matrix addition
    *
    * @code
    * A += B;
    * @endcode
    */
    void operator += (const Matrix2D<T> &op1) const {
        if (mdimx != op1.mdimx || mdimy != op1.mdimy)
            REPORT_ERROR("operator+=: Not same sizes in matrix addition");

        for (int i = 0; i < mdimy; i++)
        for (int j = 0; j < mdimx; j++)
        at(i, j) += op1.at(i, j);
    }

    /** Matrix subtraction
    *
    * @code
    * C = A - B;
    * @endcode
    */
    Matrix2D<T> operator - (const Matrix2D<T> &op1) const {
        Matrix2D<T> result;
        if (mdimx != op1.mdimx || mdimy != op1.mdimy)
            REPORT_ERROR("operator-: Not same sizes in matrix subtraction");

        result.initZeros(mdimy, mdimx);
        for (int i = 0; i < mdimy; i++)
        for (int j = 0; j < mdimx; j++)
        result(i, j) = (*this)(i, j) - op1(i, j);

        return result;
    }

    /** In-place matrix subtraction
    *
    * @code
    * A -= B;
    * @endcode
    */
    void operator -= (const Matrix2D<T> &op1) const {
        if (mdimx != op1.mdimx || mdimy != op1.mdimy)
            REPORT_ERROR("operator-=: Not same sizes in matrix subtraction");

        for (int i = 0; i < mdimy; i++)
        for (int j = 0; j < mdimx; j++)
        at(i, j) -= op1.at(i, j);
    }

    /** Equality
    *
    * Returns true if this object has the same shape (origin and size)
    * as the argument and the same values (to within machine epsilon).
    */
    bool equal(
        const Matrix2D<T> &op, RFLOAT accuracy = Xmipp::epsilon
    ) const {
        if (!sameShape(op)) 
            return false;

        for (int i = 0; i < mdimy; i++)
        for (int j = 0; j < mdimx; j++)
        if (abs((*this)(i, j) - op(i, j)) > accuracy) return false;

        return true;
    }
    //@}

    /// @name Utilities for Matrix2D
    //@{
    // Set very small values (abs(val) < accuracy) equal to zero
    void setSmallValuesToZero(RFLOAT accuracy = Xmipp::epsilon) {
        for (int i = 0; i < mdimy; i++)
        for (int j = 0; j < mdimx; j++)
        if (abs((*this)(i, j)) < accuracy) { (*this)(i, j) = 0.0; }
    }

    /// @name Utilities for Matrix2D
    //@{

    // Greatest value in an array
    T max() const {
        if (mdim <= 0) return static_cast<T>(0);

        T maxval = mdata[0];
        for (int n = 0; n < mdim; n++)
        if (mdata[n] > maxval) { maxval = mdata[n]; }

        return maxval;
    }

    // Least value in an array
    T min() const {
        if (mdim <= 0) return static_cast<T>(0);

        T minval = mdata[0];
        for (int n = 0; n < mdim; n++)
        if (mdata[n] < minval) { minval = mdata[n]; }

        return minval;
    }

    /** Produce a 2D array suitable for working with Numerical Recipes
    *
    * This function must be used only as a preparation for routines which need
    * that the first physical index is 1 and not 0 as it usually is in C. New
    * memory is needed to hold the new RFLOAT pointer array.
    */
    T** adaptForNumericalRecipes() const {
        T **m = NULL;
        ask_Tmatrix(m, 1, mdimy, 1, mdimx);

        for (int i = 0; i < mdimy; i++)
        for (int j = 0; j < mdimx; j++)
        m[i + 1][j + 1] = mdata[i * mdimx + j];

        return m;
    }

    /** Produce a 1D pointer suitable for working with Numerical Recipes (2)
    *
    * This function meets the same goal as the one before,
    * however this one works with 2D arrays as a single pointer.
    * result[i * Xdim + j]
    * result[1 * Xdim + 1] points to the first element of the array,
    */
    T* adaptForNumericalRecipes2() const {
        return mdata - 1 - mdimx;
    }

    // Load 2D array from numerical recipes result
    void loadFromNumericalRecipes(T **m, int Ydim, int Xdim) {
        resize(Ydim, Xdim);

        for (int i = 1; i <= Ydim; i++)
        for (int j = 1; j <= Xdim; j++)
        (*this)(i - 1, j - 1) = m[i][j];
    }

    // Kill a 2D array produced for numerical recipes
    void killAdaptationForNumericalRecipes(T **m) const {
        // Free the allocated memory
        free_Tmatrix(m, 1, mdimy, 1, mdimx);
    }

    void killAdaptationForNumericalRecipes2(T **m) const {
        // Do nothing
    }

    // Write this matrix to file
    void write(const FileName &fn) const {
        std::ofstream fhOut(fn.c_str());
        if (!fhOut)
            REPORT_ERROR((std::string) "write: Cannot open " + fn + " for output");
        fhOut << *this;
    }

    // Show matrix
    friend std::ostream& operator << (std::ostream &ostrm, const Matrix2D<T> &v) {
        if (v.Xdim() == 0 || v.Ydim() == 0) {
            ostrm << "NULL matrix\n";
        } else {
            ostrm << std::endl;
            RFLOAT max_val = v.max();
            int epsilon = bestPrecision(max_val, 10);
            for (int i = 0; i < v.Ydim(); i++) {
                for (int j = 0; j < v.Xdim(); j++) {
                    ostrm << std::setw(13) << floatToString((RFLOAT) v(i, j), 10, epsilon) << ' ';
                }
                ostrm << std::endl;
            }
        }
        return ostrm;
    }

    /** Makes a matrix from a vector
    *
    * The origin of the matrix is set such that it has one of the index origins
    * (X or Y) to the same value as the vector, and the other set to 0
    * according to the shape.
    *
    * @code
    * Matrix2D<RFLOAT> m = fromVector(v);
    * @endcode
    */
    void fromVector(const Matrix1D<T> &op1) {
        // Null vector => Null matrix
        if (op1.size() == 0) {
            clear();
            return;
        }

        // Look at shape and copy values
        if (op1.isRow()) {
            resize(1, op1.size());
            for (int j = 0; j < op1.size(); j++)
                at(0, j) = op1[j];
        } else {
            resize(op1.size(), 1);
            for (int i = 0; i < op1.size(); i++)
                at(i, 0) = op1[i];
        }
    }

    /** Make a vector from a matrix
    *
    * An exception is thrown if the matrix is not a single row or a single
    * column. The origin of the vector is set according to the one of the
    * matrix.
    *
    * @code
    * Matrix1D<RFLOAT> v;
    * m.toVector(v);
    * @endcode
    */
    void toVector(Matrix1D<T> &op1) const {
        // Null matrix => Null vector
        if (mdimx == 0 || mdimy == 0) {
            op1.clear();
            return;
        }

        // If matrix is not a vector, produce an error
        if (mdimx != 1 && mdimy != 1)
            REPORT_ERROR("toVector: Matrix cannot be converted to vector");

        // Look at shape and copy values
        if (mdimy == 1) {
            // Row vector
            op1.resize(mdimx);

            for (int j = 0; j < mdimx; j++) { op1[j] = at(0, j); }

            op1.setRow();
        } else {
            // Column vector
            op1.resize(mdimy);

            for (int i = 0; i < mdimy; i++) { op1[i] = at(i, 0); }

            op1.setCol();
        }
    }

    // Copy matrix to stl::vector
    void copyToVector(std::vector<T> &v) {
        v.assign(mdata, mdata + mdim);
    }

    // Copy stl::vector to matrix
    void copyFromVector(std::vector<T> &v, int Xdim, int Ydim) {
        resize(Ydim, Xdim);
        copy(v.begin(), v.begin() + v.size(), mdata);
    }

    /** Get row
    *
    * This function returns a row vector corresponding to the choosen
    * row inside the nth 2D matrix, the numbering of the rows is also
    * logical not physical.
    *
    * @code
    * std::vector<RFLOAT> v;
    * m.getRow(-2, v);
    * @endcode
    */
    void getRow(int i, Matrix1D<T> &v) const {
        if (mdimx == 0 || mdimy == 0) {
            v.clear();
            return;
        }

        if (i < 0 || i >= mdimy)
            REPORT_ERROR("getRow: index out of matrix bounds");
            // std::out_of_range ?

        for (int j = 0; j < mdimx; j++) { v[j] = at(i, j); }

        v.setRow();
    }

    /** Get Column
    *
    * This function returns a column vector corresponding to the
    * choosen column.
    *
    * @code
    * std::vector<RFLOAT> v;
    * m.getCol(-1, v);
    * @endcode
    */
    void getCol(int j, Matrix1D<T> &v) const {
        if (mdimx == 0 || mdimy == 0) {
            v.clear();
            return;
        }

        if (j < 0 || j >= mdimx)
            REPORT_ERROR("getCol: index outside matrix bounds");
            // std::out_of_range() ?

        for (int i = 0; i < mdimy; i++) { v(i) = at(i, j); }

        v.setCol();
    }

    /** Set Row
    *
    * This function sets a row vector corresponding to the choosen row in the 2D Matrix
    *
    * @code
    * m.setRow(-2, m.row(1)); // Copies row 1 in row -2
    * @endcode
    */
    void setRow(int i, const Matrix1D<T> &v) {
        if (mdimx == 0 || mdimy == 0)
            REPORT_ERROR("setRow: Target matrix is empty");

        if (i < 0 || i >= mdimy)
            REPORT_ERROR("setRow: Matrix subscript (i) out of range");

        if (v.size()  != mdimx)
            REPORT_ERROR("setRow: Vector dimension different from matrix one");

        if (!v.isRow())
            REPORT_ERROR("setRow: Not a row vector in assignment");

        for (int j = 0; j < mdimx; j++)
            at(i, j) = v(j);
    }

    /** Set Column
    *
    * This function sets a column vector corresponding to the choosen column
    * inside matrix.
    *
    * @code
    * m.setCol(0, (m.row(1)).transpose()); // Copies row 1 in column 0
    * @endcode
    */
    void setCol(int j, const Matrix1D<T> &v) {
        if (mdimx == 0 || mdimy == 0)
            REPORT_ERROR("setCol: Target matrix is empty");

        if (j < 0 || j >= mdimx)
            REPORT_ERROR("setCol: Matrix subscript (j) out of range");

        if (v.size() != mdimy)
            REPORT_ERROR("setCol: Vector dimension different from matrix one");

        if (!v.isCol())
            REPORT_ERROR("setCol: Not a column vector in assignment");

        for (int i = 0; i < mdimy; i++)
            at(i, j) = v[i];
    }

    /** Matrix determinant
    *
    * An exception is thrown if the matrix is empty or non-square.
    *
    * @code
    * RFLOAT det = m.det();
    * @endcode
    */
    T det() const {
        // (see Numerical Recipes, Chapter 2 Section 5)
        if (mdimx == 0 || mdimy == 0)
            REPORT_ERROR("determinant: Matrix is empty");

        if (mdimx != mdimy)
            REPORT_ERROR("determinant: Matrix is not square");

        for (int i = 0; i < mdimy; i++) {
            bool all_zeros = true;
            for (int j = 0; j < mdimx; j++) {
                if (abs(at(i, j)) > Xmipp::epsilon) {
                    all_zeros = false;
                    break;
                }
            }
            if (all_zeros) {
                return 0;
            }
        }

        // Perform decomposition
        Matrix1D<int> indx;
        T d;
        Matrix2D<T> LU;
        ludcmp(*this, LU, indx, d);

        // Calculate determinant
        for (int i = 0; i < mdimx; i++)
        d *= (T) LU.at(i, i);

        return d;
    }

    /** Algebraic transpose of a matrix
    *
    * You can use the transpose in as complex expressions as you like.
    * The origin of the vector is not changed.
    *
    * @code
    * v2 = v1.transpose();
    * @endcode
    */
    Matrix2D<T> transpose() const {
        Matrix2D<T> result(mdimx, mdimy);
        FOR_ALL_ELEMENTS_IN_MATRIX2D(result)
            result.at(i, j) = at(j, i);
        return result;
    }

    /** Matrix pseudoinverse
    * https://en.wikipedia.org/wiki/Moore–Penrose_inverse
    *
    * Compute the pseudoinverse of a matrix
    * using SVD.
    *
    * @code
    * Matrix2D<RFLOAT> m1_inv;
    * m1.inv(m1_inv);
    * @endcode
    */
    void inv(Matrix2D<T> &result) const;
    // Why modify a pre-existing matrix?
    // Why not just return a new one?

    // Matrix inverse
    Matrix2D<T> inv() const {
        Matrix2D<T> result;
        inv(result);
        return result;
    }

    /** Test for identity matrix
    *
    * @code
    * if (m.isIdentity())
    *     std::cout << "The matrix is identity\n";
    * @endcode
    */
    bool isIdentity() const {
        for (int i = 0; i < mdimy; i++)
        for (int j = 0; j < mdimx; j++) {
            T elem = at(i, j);
            if (abs(i == j ? elem - 1.0 : elem) > Xmipp::epsilon)
                return false;
        }
        return true;
    }
    //@}

};

// vector * matrix
// Declared in matrix1D.h
template<typename T>
Matrix1D<T> Matrix1D<T>::operator * (const Matrix2D<T> &M) {

    if (size() != M.mdimy)
        REPORT_ERROR("Not compatible sizes in matrix by vector");

    if (!isRow())
        REPORT_ERROR("Vector is not a row");

    Matrix1D<T> result = Matrix1D<T>::zeros(M.mdimx);
    for (int j = 0; j < M.mdimx; j++)
    for (int i = 0; i < M.mdimy; i++)
    result[j] += (*this)[i] * M.at(i, j);

    result.setRow();
    return result;
}

/** @name Matrix-related functions
 * These functions are not methods of Matrix2D
 */
//@{
// LU Decomposition
template<typename T>
void ludcmp(const Matrix2D<T> &A, Matrix2D<T> &LU, Matrix1D<int> &indx, T &d) {
    LU = A;
    indx.resize(A.mdimx);
    ludcmp(
        LU.adaptForNumericalRecipes2(), A.mdimx,
        indx.adaptForNumericalRecipes(), &d
    );
}

// LU Backsubstitution
template<typename T>
void lubksb(const Matrix2D<T> &LU, Matrix1D<int> &indx, Matrix1D<T> &b) {
    lubksb(
        LU.adaptForNumericalRecipes2(), indx.size(),
        indx.adaptForNumericalRecipes(), b.adaptForNumericalRecipes()
    );
}

// SVD Backsubstitution
void svbksb(
    Matrix2D<RFLOAT> &u, Matrix1D<RFLOAT> &w, Matrix2D<RFLOAT> &v,
    Matrix1D<RFLOAT> &b, Matrix1D<RFLOAT> &x
);

// Singular Value Decomposition (from numerical_recipes)
template<typename T>
void svdcmp(
    const Matrix2D<T> &a,
    Matrix2D<RFLOAT> &u, Matrix1D<RFLOAT> &w, Matrix2D<RFLOAT> &v
) {
    // svdcmp only works with RFLOAT
    typeCast(a, u);

    // Set size of matrices
    w.initZeros(u.mdimx);
    v.initZeros(u.mdimx, u.mdimx);

    // Call the numerical recipes routine
    svdcmp(u.mdata, u.mdimy, u.mdimx, w.vdata, v.mdata);
}

// Solve a system of linear equations (Ax = b) by SVD
template<typename T>
void solve(
    const Matrix2D<T> &A, const Matrix1D<T> &b,
    Matrix1D<RFLOAT> &result, RFLOAT tolerance
);

// Solve a system of linear equations (Ax=b), where x and b are matrices,
// by SVD Decomposition (through Gauss-Jordan numerical recipes)
template<typename T>
void solve(const Matrix2D<T> &A, const Matrix2D<T> &b, Matrix2D<T> &result) {
    if (A.mdimx == 0)
        REPORT_ERROR("Solve: Matrix is empty");

    if (A.mdimx != A.mdimy)
        REPORT_ERROR("Solve: Matrix is not square");

    if (A.mdimy != b.mdimy)
        REPORT_ERROR("Solve: Different sizes of A and b");

    // Solve
    result = b;
    Matrix2D<T> Aux = A;
    gaussj(
        Aux.adaptForNumericalRecipes2(), Aux.mdimy,
        result.adaptForNumericalRecipes2(), b.mdimx
    );
}


/** Least-squares rigid transformation between two sets of 3D coordinates
 *
RFLOAT lsq_rigid_body_transformation(std::vector<Matrix1D<RFLOAT> > &set1, std::vector<Matrix1D<RFLOAT> > &set2,
        Matrix2D<RFLOAT> &Rot, Matrix1D<RFLOAT> &trans) {
    Matrix2D<RFLOAT> A;
    Matrix1D<RFLOAT> avg1, avg2;

    if (set1.size() != set2.size())
        REPORT_ERROR("lsq_rigid_body_transformation ERROR: unequal set size");

    // Calculate average of set1 and set2
    avg1 = vectorR3(0., 0., 0.);
    avg2 = vectorR3(0., 0., 0.);
    for (int i = 0; i < set1.size(); i++) {
        if (set1[i].vdim != 3)
            REPORT_ERROR("lsq_rigid_body_transformation ERROR: not a 3-point set1");
        if (set2[i].vdim != 3)
            REPORT_ERROR("lsq_rigid_body_transformation ERROR: not a 3-point set2");
        avg1 += set1[i];
        avg2 += set2[i];
    }
    avg1 /= (RFLOAT)set1.size();
    avg2 /= (RFLOAT)set1.size();

    A.initZeros(3, 3);
    Rot.initZeros(4,4);
    for (int i = 0; i < set1.size(); i++) {
        // fill A
        A(0, 0) += (XX(set1[i]) - XX(avg1)) * (XX(set2[i]) - XX(avg2));
        A(0, 1) += (XX(set1[i]) - XX(avg1)) * (YY(set2[i]) - YY(avg2));
        A(0, 2) += (XX(set1[i]) - XX(avg1)) * (ZZ(set2[i]) - ZZ(avg2));
        A(1, 0) += (YY(set1[i]) - YY(avg1)) * (XX(set2[i]) - XX(avg2));
        A(1, 1) += (YY(set1[i]) - YY(avg1)) * (YY(set2[i]) - YY(avg2));
        A(1, 2) += (YY(set1[i]) - YY(avg1)) * (ZZ(set2[i]) - ZZ(avg2));
        A(2, 0) += (ZZ(set1[i]) - ZZ(avg1)) * (XX(set2[i]) - XX(avg2));
        A(2, 1) += (ZZ(set1[i]) - ZZ(avg1)) * (YY(set2[i]) - YY(avg2));
        A(2, 2) += (ZZ(set1[i]) - ZZ(avg1)) * (ZZ(set2[i]) - ZZ(avg2));
    }

    Matrix2D<RFLOAT> U, V;
    Matrix1D<RFLOAT> w;

    // TODO: check inverse, transpose etc etc!!!

    // Optimal rotation
    svdcmp(A, U, w, V);
    Rot = V.transpose() * U;

    // Optimal translation
    trans = avg1 - Rot * avg2;

    // return the squared difference term
    RFLOAT error = 0.;
    for (int i = 0; i < set1.size(); i++)
    {
        error += (Rot * set2[i] + trans - set1[i]).sum2();
    }

    return error;

}
*/

/** Type casting
 *
 * If we have an integer array and we need a RFLOAT one, we can use this function.
 * The conversion is done by type casting each element.
 * If n >= 0, only the nth volumes will be converted, otherwise all NSIZE volumes.
 */
template<typename T1, typename T2>
void typeCast(const Matrix2D<T1> &v1,  Matrix2D<T2> &v2) {
    if (v1.mdim == 0) {
        v2.clear();
        return;
    }

    v2.resize(v1);
    for (unsigned long int n = 0; n < v1.mdim; n++)
        v2.mdata[n] = static_cast<T2>(v1.mdata[n]);
}
//@}
//@}
#endif /* MATRIX2D_H_ */
