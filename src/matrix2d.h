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

//@{

template<typename T>
class Matrix {

    // The array itself
    T* mdata;

    // Number of columns, rows
    int mdimx, mdimy;

    public:

    /// @name Constructors
    /// @{

    // Default / dimension constructor
    Matrix(int m = 0, int n = 0): mdimx(0), mdimy(0), mdata(nullptr) {
        resize(m, n);
    }

    // Copy constructor
    Matrix(const Matrix<T> &other): Matrix(other.ncols(), other.nrows()) {
        std::copy(other.begin(), other.end(), begin());
    }

    // Type-casting copy constructor
    template<typename U>
    Matrix(const Matrix<U> &other): Matrix(other.ncols(), other.nrows()) {
        std::copy(other.begin(), other.end(), begin());
    }

    // Destructor
    ~Matrix() { clear(); }

    friend void swap(Matrix<T> &lhs, Matrix<T> &rhs) {
        std::swap(lhs.mdimx, rhs.mdimx);
        std::swap(lhs.mdimy, rhs.mdimy);
        std::swap(lhs.mdata, rhs.mdata);
    }

    Matrix<T>& operator = (Matrix<T> rhs) {
        swap(*this, rhs);
        return *this;
    }
    //@}

    // Clear
    void clear() {
        delete[] mdata;
        mdata = nullptr;
        mdimx = mdimy = 0;
    }

    /// @name Size and shape of Matrix
    //@{
    // Resize to a given size
    void resize(int new_mdimx, int new_mdimy) {

        if (new_mdimx == mdimx && new_mdimy == mdimy) return;

        if (new_mdimx <= 0 || new_mdimy <= 0) {
            clear();
            return;
        }

        T *new_mdata;
        try {
            new_mdata = new T[new_mdimx * new_mdimy];
        } catch (std::bad_alloc &) {
            REPORT_ERROR("Allocate: No space left");
        }

        for (int i = 0; i < new_mdimy; i++)
        for (int j = 0; j < new_mdimx; j++)
        // Copy old data
        // Zero-pad where necessary
            new_mdata[i * new_mdimx + j] = i >= mdimy || j >= mdimx ? 0 :
                mdata[i * mdimx + j];

        // Deallocate old data
        clear();

        // Assign member variables
        mdata = new_mdata;
        mdimx = new_mdimx;
        mdimy = new_mdimy;
    }

    // Extract submatrix
    Matrix<T> submatrix(int i0, int j0, int iF, int jF) const {
        if (i0 < 0 || j0 < 0 || iF >= nrows() || jF >= ncols())
            REPORT_ERROR("Submatrix indices out of bounds");

        Matrix<T> A (iF - i0 + 1, jF - j0 + 1);
        for (int i = 0; i < A.nrows(); i++) 
        for (int j = 0; j < A.ncols(); j++)
            A.at(i, j) = at(i + i0, j + j0);
        return A;
    }

    int size() const { return mdimx * mdimy; }

    const T* data() const { return mdata; }
          T* data()       { return mdata; }

    std::pair<int, int> shape() const { return {mdimx, mdimy}; }

    // Number of columns
    inline int ncols() const { return mdimx; }

    // Number of rows
    inline int nrows() const { return mdimy; }

    //@}

    /// @name Initialise Matrix values
    //@{

    static Matrix<T> zeros(int m, int n) {
        Matrix<T> A (m, n);
        std::fill(A.begin(), A.end(), 0);
        return A;
    }

    static Matrix<T> identity(int n) {
        Matrix<T> A (n, n);
        A.setIdentity();
        return A;
    }

    /** Identity matrix
    * 1 along diagonal, 0 everywhere else.
    * The matrix will not be resized.
    */
    void setIdentity() {
        const auto n = std::min(ncols(), nrows());
        std::fill(begin(), end(), 0);
        for (int i = 0; i < n; i++) { at(i, i) = 1; }
    }
    //@}

    /// @name Operators for Matrix
    //@{

    /// Subscripting

    inline T  at(int i, int j) const { return mdata[i * mdimx + j]; }
    inline T& at(int i, int j)       { return mdata[i * mdimx + j]; }

    inline T  operator () (int i, int j) const { return at(i, j); }
    inline T& operator () (int i, int j)       { return at(i, j); }

    inline T  operator [] (int i) const { return mdata[i]; }
    inline T& operator [] (int i)       { return mdata[i]; }

    Matrix<T>& operator += (T rhs) {
        for (auto& x: *this) { x += rhs; }
        return *this;
    }

    Matrix<T>& operator -= (T rhs) {
        for (auto& x: *this) { x -= rhs; }
        return *this;
    }

    Matrix<T>& operator *= (T rhs) {
        for (auto& x: *this) { x *= rhs; }
        return *this;
    }

    Matrix<T>& operator /= (T rhs) {
        for (auto& x: *this) { x /= rhs; }
        return *this;
    }

    Matrix<T> matmul(const Matrix<T>& rhs) const {
        const Matrix<T>& lhs = *this;
        if (lhs.ncols() != rhs.nrows())
            REPORT_ERROR("Incompatible shapes in matrix multiplication");

        Matrix<T> product = zeros(lhs.nrows(), rhs.ncols());
        for (int i = 0; i < lhs.nrows(); i++)
        for (int j = 0; j < rhs.ncols(); j++)
        for (int k = 0; k < lhs.ncols(); k++)
            product.at(i, j) += lhs.at(i, k) * rhs.at(k, j);
        return product;
    }

    Matrix<T>& operator += (const Matrix<T> &rhs) {
        if (shape() != rhs.shape())
            REPORT_ERROR("operator+=: Not same sizes in matrix addition");
        for (int i = 0; i < size(); i++) mdata[i] += rhs.mdata[i];
        return *this;
    }

    Matrix<T>& operator -= (const Matrix<T> &rhs) {
        if (shape() != rhs.shape())
            REPORT_ERROR("operator-=: Not same sizes in matrix subtraction");
        for (int i = 0; i < size(); i++) mdata[i] -= rhs.mdata[i];
        return *this;
    }

    /** Equality
    *
    * Returns true if this object has the same shape as the argument
    * and the same values (to within machine epsilon).
    */
    bool equal(const Matrix<T> &other, RFLOAT accuracy = Xmipp::epsilon) const {
        if (shape() != other.shape()) 
            return false;

        for (int i = 0; i < nrows(); i++)
        for (int j = 0; j < ncols(); j++)
        if (abs((*this)(i, j) - other(i, j)) > accuracy) return false;

        return true;
    }
    //@}

    inline const T* begin() const { return mdata; }
    inline       T* begin()       { return mdata; }

    inline const T* end() const { return mdata + mdimx * mdimy; }
    inline       T* end()       { return mdata + mdimx * mdimy; }

    /// @name Utilities for Matrix
    //@{

    /** Construct from a vector
    *
    * The origin of the matrix is set such that it has one of the index origins
    * (X or Y) to the same value as the vector, and the other set to 0
    * according to the shape.
    */
    explicit Matrix(const Vector<T> &v): Matrix() {
        if (v.isRow()) {
            resize(v.size(), 1);
        } else {
            resize(1, v.size());
        }
        std::copy(v.begin(), v.end(), begin());
    }

    /** Construct a vector
    *
    * An exception is thrown if the matrix is not a row or column vector.
    */
    explicit operator Vector<T>() const {
        if (nrows() == 1) {
            Vector<T> v (ncols(), VectorMode::row);  // Row vector
            std::copy(begin(), end(), v.begin());
            return v;
        }
        if (ncols() == 1) {
            Vector<T> v (nrows(), VectorMode::column);  // Column vector
            std::copy(begin(), end(), v.begin());
            return v;
        }
        // Otherwise, throw
        REPORT_ERROR("toVector: Matrix cannot be converted to vector");
    }

    operator std::vector<T>() const { return {begin(), end()}; }

    Matrix(const std::vector<T>& v, int m, int n): Matrix(m, n) {
        // m * n had better be at least v.size()
        copy(v.begin(), v.end(), mdata);
    }

    /** Get row
    *
    * This function returns a row vector corresponding to the choosen
    * row inside the nth 2D matrix, the numbering of the rows is also
    * logical not physical.
    */
    Vector<T> getRow(int i) const {
        if (i < 0 || i >= nrows())
            REPORT_ERROR("getRow: index out of matrix bounds");
            // std::out_of_range ?

        Vector<T> v (ncols(), VectorMode::row);
        for (int j = 0; j < ncols(); j++) { v[j] = at(i, j); }
        return v;
    }

    /** Get column
    *
    * This function returns a column vector corresponding to the
    * choosen column.
    */
    Vector<T> getCol(int j) const {
        if (j < 0 || j >= ncols())
            REPORT_ERROR("getCol: index outside matrix bounds");
            // std::out_of_range() ?

        Vector<T> v (nrows(), VectorMode::column);
        for (int i = 0; i < nrows(); i++) { v[i] = at(i, j); }
        return v;
    }

    /** Set row
    *
    * This function sets a row vector corresponding to the choosen row in the 2D Matrix
    *
    * @code
    * m.setRow(-2, m.row(1));  // Copies row 1 in row -2
    * @endcode
    */
    void setRow(int i, const Vector<T> &v) {
        if (i < 0 || i >= nrows())
            REPORT_ERROR("setRow: Matrix subscript (i) out of range");

        if (v.size() != ncols())
            REPORT_ERROR("setRow: Vector dimension different from matrix one");

        if (!v.isRow())
            REPORT_ERROR("setRow: Not a row vector in assignment");

        for (int j = 0; j < ncols(); j++)
            at(i, j) = v[j];
    }

    /** Set column
    *
    * This function sets a column vector corresponding to the choosen column
    * inside matrix.
    *
    * @code
    * m.setCol(0, m.row(1).transpose());  // Copies row 1 in column 0
    * @endcode
    */
    void setCol(int j, const Vector<T> &v) {
        if (j < 0 || j >= ncols())
            REPORT_ERROR("setCol: Matrix subscript (j) out of range");

        if (v.size() != nrows())
            REPORT_ERROR("setCol: Vector dimension different from matrix one");

        if (!v.isCol())
            REPORT_ERROR("setCol: Not a column vector in assignment");

        for (int i = 0; i < nrows(); i++)
            at(i, j) = v[i];
    }

    /** Matrix determinant
    *
    * An exception is thrown if the matrix is empty or non-square.
    */
    T det() const {
        // (see Numerical Recipes, Chapter 2 Section 5)
        if (ncols() != nrows())
            REPORT_ERROR("determinant: Matrix is not square");

        if (size() == 0) return 1;

        // If any row is all zeros, then the determinant is zero.
        for (int i = 0; i < nrows(); i++) {
            for (int j = 0; j < ncols(); j++) {
                if (abs(at(i, j)) > Xmipp::epsilon) goto next;
            }
            return 0;
            next: continue;
        }

        // Perform decomposition
        Vector<int> indx;
        T d;
        Matrix<T> LU;
        ludcmp(*this, LU, indx, d);

        // Calculate determinant
        for (int i = 0; i < ncols(); i++)
            d *= LU.at(i, i);

        return d;
    }

    /** Algebraic transpose of a matrix
    */
    Matrix<T> transpose() const {
        /// XXX: Could be done in place.
        Matrix<T> t (ncols(), nrows());
        for (int i = 0; i < t.nrows(); i++)
        for (int j = 0; j < t.ncols(); j++)
            t.at(i, j) = at(j, i);
        return t;
    }

    /** Matrix pseudoinverse
    * https://en.wikipedia.org/wiki/Moore–Penrose_inverse
    *
    * Compute the pseudoinverse of a matrix by SVD.
    */
    Matrix<T> inv() const;

    // Test for identity matrix
    bool isIdentity() const {
        for (int i = 0; i < nrows(); i++)
        for (int j = 0; j < ncols(); j++) {
            if (abs(i == j ? at(i, j) - 1.0 : at(i, j)) > Xmipp::epsilon)
                return false;
        }
        return true;
    }
    //@}

};

namespace nr {

    template <typename T>
    T* adaptForNumericalRecipes(Vector<T>& vec) {
        return vec.data() - 1;
    }

    /** Produce a 2D array suitable for working with Numerical Recipes
    *
    * This function must be used only as a preparation for routines which need
    * that the first physical index is 1 and not 0 as it usually is in C. New
    * memory is needed to hold the new pointer array.
    */
   template <typename T>
    T** adaptForNumericalRecipes(const Matrix<T>& mat) {
        T** ptr = ask_matrix<T>(1, mat.nrows(), 1, mat.ncols());
        for (int i = 0; i < mat.nrows(); i++)
        for (int j = 0; j < mat.ncols(); j++)
            ptr[i + 1][j + 1] = mat.at(i, j);
        return ptr;
    }

    // Kill a 2D array produced for numerical recipes
    template <typename T>
    void killAdaptationForNumericalRecipes(const Matrix<T>& mat, T **ptr) {
        // Free the allocated memory
        free_matrix<T>(ptr, 1, mat.nrows(), 1, mat.ncols());
    }

    // Load 2D array from numerical recipes result
    template <typename T>
    void loadFromNumericalRecipes(Matrix<T> &mat, T **ptr, int m, int n) {
        mat.resize(m, n);
        for (int i = 1; i <= n; i++)
        for (int j = 1; j <= m; j++)
            mat.at(i - 1, j - 1) = ptr[i][j];
    }

    /** Produce a 1D pointer suitable for working with Numerical Recipes (2)
    *
    * This function meets the same goal as the one before,
    * however this one works with 2D arrays as a single pointer.
    * ptr[i * mat.ncols() + j] = mat.at(i - 1, j - 1)
    */
   template <typename T>
    T* adaptForNumericalRecipes2(Matrix<T>& mat) {
        return mat.data() - 1 - mat.nrows();
    }

}

// Free functions

template <typename T>
Vector<T> matmul(const Matrix<T> &lhs, const Vector<T> &rhs) {

    if (lhs.ncols() != rhs.size())
        REPORT_ERROR("Incompatible shapes in matrix by vector multiplication");

    if (!rhs.isCol())
        REPORT_ERROR("Right operand is not a column vector");

    Vector<T> result (lhs.nrows(), VectorMode::column);
    for (int i = 0; i < result.size(); i++) {
        result[i] = 0;
        for (int j = 0; j < rhs.size(); j++)
            result[i] += lhs(i, j) * rhs[j];
    }
    return result;
}

template<typename T>
Vector<T> matmul(const Vector<T> &lhs, const Matrix<T> &rhs) {

    if (lhs.size() != rhs.nrows())
        REPORT_ERROR("Incompatible shapes in vector by matrix multiplication");

    if (!lhs.isRow())
        REPORT_ERROR("Left operand is not a row vector");

    Vector<T> result (rhs.ncols(), VectorMode::row);
    for (int j = 0; j < rhs.ncols(); j++) {
        result[j] = 0;
        for (int i = 0; i < rhs.nrows(); i++)
            result[j] += lhs[i] * rhs(i, j);
    }
    return result;
}

/** @name Matrix-related functions
 * These functions are not methods of Matrix
 */
//@{
// LU Decomposition
template<typename T>
void ludcmp(const Matrix<T> &A, Matrix<T> &LU, Vector<int> &indx, T &d) {
    LU = A;
    indx.resize(A.ncols());
    ludcmp(
        nr::adaptForNumericalRecipes2(LU), A.ncols(),
        indx.data() - 1, &d
    );
}

// LU Backsubstitution
template<typename T>
void lubksb(const Matrix<T> &LU, Vector<int> &indx, Vector<T> &b) {
    lubksb(
        nr::adaptForNumericalRecipes2(LU), indx.size(),
        indx.data() - 1, b.data() - 1
    );
}

// SVD Backsubstitution
void svbksb(
    Matrix<RFLOAT> &u, Vector<RFLOAT> &w, Matrix<RFLOAT> &v,
    Vector<RFLOAT> &b, Vector<RFLOAT> &x
);

// Singular Value Decomposition (from numerical_recipes)
template<typename T>
void svdcmp(
    const Matrix<T> &a,
    Matrix<RFLOAT> &u, Vector<RFLOAT> &w, Matrix<RFLOAT> &v
);

// Solve a system of linear equations (Ax = b) by SVD
template<typename T>
void solve(const Matrix<T> &A, Vector<T> b, Vector<RFLOAT> &x, RFLOAT epsilon = Xmipp::epsilon);

// Solve a system of linear equations (Ax=b), where x and b are matrices,
// by SVD Decomposition (through Gauss-Jordan numerical recipes)
template<typename T>
void solve(Matrix<T> A, const Matrix<T> &b, Matrix<T> &x);

/** Least-squares rigid transformation between two sets of 3D coordinates
 *
RFLOAT lsq_rigid_body_transformation(std::vector<Vector<RFLOAT> > &set1, std::vector<Vector<RFLOAT> > &set2,
        Matrix<RFLOAT> &Rot, Vector<RFLOAT> &trans) {
    Matrix<RFLOAT> A;
    Vector<RFLOAT> avg1, avg2;

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

    A.resize(3, 3);
    std::fill(A.begin(), A.end(), 0);
    Rot.resize(4, 4);
    std::fill(Rot.begin(), Rot.end(), 0);
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

    Matrix<RFLOAT> U, V;
    Vector<RFLOAT> w;

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

template <typename T>
Matrix<T> operator + (Matrix<T> lhs, T rhs) {
    return lhs += rhs;
}

template <typename T>
Matrix<T> operator - (Matrix<T> lhs, T rhs) {
    return lhs -= rhs;
}

template <typename T>
Matrix<T> operator * (Matrix<T> lhs, T rhs) {
    return lhs *= rhs;
}

template <typename T>
Matrix<T> operator / (Matrix<T> lhs, T rhs) {
    return lhs /= rhs;
}

template <typename T>
Matrix<T> operator * (T lhs, Matrix<T> rhs) {
    for (auto& x: rhs) x = lhs * x;
    return rhs;
}

template <typename T>
Matrix<T> operator + (Matrix<T> lhs, const Matrix<T> &rhs) {
    return lhs += rhs;
}

template <typename T>
Matrix<T> operator - (Matrix<T> lhs, const Matrix<T> &rhs) {
    return lhs -= rhs;
}

template <typename It>
void setSmallValuesToZero(It here, const It& end, RFLOAT epsilon = Xmipp::epsilon) {
    for (; here != end; ++here)
        if (abs(*here) < epsilon)
            *here = 0;
}

// Show matrix
template <typename T>
std::ostream& operator << (std::ostream &ostrm, const Matrix<T> &m) {
    if (m.size() == 0) {
        ostrm << "Empty matrix\n";
    } else {
        ostrm << '\n';
        const T maximum = *std::max_element(m.begin(), m.end());
        const int epsilon = bestPrecision(maximum, 10);
        for (int i = 0; i < m.nrows(); i++) {
        for (int j = 0; j < m.ncols(); j++) {
            ostrm << std::setw(13) << floatToString((RFLOAT) m.at(i, j), 10, epsilon) << ' ';
            }
            ostrm << std::endl;
        }
    }
    return ostrm;
}

// Write matrix to file
template <typename T>
void write(const Matrix<T>& matrix, const FileName& fn) {
    std::ofstream fhOut (fn.c_str());
    if (!fhOut)
        REPORT_ERROR((std::string) "write: Cannot open " + fn + " for output");
    fhOut << matrix;
}

//@}
//@}
#endif /* MATRIX2D_H_ */
