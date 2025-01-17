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

#ifndef MULTIDIM_ARRAY_H
#define MULTIDIM_ARRAY_H

#include "src/allocators.h"
#include "src/funcs.h"
#include "src/error.h"
#include "src/args.h"
#include "src/matrix1d.h"
#include "src/matrix2d.h"
#include "src/complex.h"

extern int bestPrecision(float F, int _width);
extern std::string floatToString(float F, int _width, int _prec);

/// @defgroup MultidimensionalArrays Multidimensional Arrays
/// @ingroup DataLibrary
//@{
/** @name MultidimArraysSpeedUp Speed up macros
 *
 * This macros are defined to allow high speed in critical parts of your
 * program. They shouldn't be used systematically as usually there is no
 * checking on the correctness of the operation you are performing. Speed comes
 * from three facts: first, they are macros and no function call is performed
 * (although most of the critical functions are inline functions), there is no
 * checking on the correctness of the operation (it could be wrong and you are
 * not warned of it), and destination vectors are not returned saving time in
 * the copy constructor and in the creation/destruction of temporary vectors.
 */
//@{

template <typename T, typename Allocator=relion_aligned_mallocator>
class MultidimArray;

/** Return the first X valid logical index
 */
template <typename T>
inline long int Xinit(const MultidimArray<T> &v) { return v.xinit; }

/** Return the last X valid logical index
 */
template <typename T>
inline long int Xlast(const MultidimArray<T> &v) { return v.xinit + v.xdim - 1; }

/** Return the first Y valid logical index
 */
template <typename T>
inline long int Yinit(const MultidimArray<T> &v) { return v.yinit; }

/** Return the last Y valid logical index
 */
template <typename T>
inline long int Ylast(const MultidimArray<T> &v) { return v.yinit + v.ydim - 1; }

/** Return the first Z valid logical index
 */
template <typename T>
inline long int Zinit(const MultidimArray<T> &v) { return v.zinit; }

/** Return the last Z valid logical index
 */
template <typename T>
inline long int Zlast(const MultidimArray<T> &v) { return v.zinit + v.zdim - 1; }

/** Access to X dimension (size)
 */
template <typename T>
inline long int Xsize(const MultidimArray<T> &v) { return v.xdim; }

/** Access to Y dimension (size)
 */
template <typename T>
inline long int Ysize(const MultidimArray<T> &v) { return v.ydim; }

/** Access to Z dimension (size)
 */
template <typename T>
inline long int Zsize(const MultidimArray<T> &v) { return v.zdim; }

/** Access to N dimension (size)
 */
template <typename T>
inline long int Nsize(const MultidimArray<T> &v) { return v.ndim; }

/** Looping
 *
 * @code
 * for (long int n = 0; n < v.size(); n++) {
 *     std::cout << v[n] << " ";
 * }
 * @endcode
 */

/** For all direct elements in the array
 *
 * This macro is used to generate loops for the array in an easy manner.
 * It binds the names provided to direct indices which iterate over the data.
 *
 * @code
 * FOR_ALL_DIRECT_NZYX_ELEMENTS_IN_MULTIDIMARRAY(v, i, j, k, l) {
 *     std::cout << direct::elem(v, i, j, k, l) << " ";
 * }
 * @endcode
 */
#define FOR_ALL_DIRECT_NZYX_ELEMENTS_IN_MULTIDIMARRAY(V, i, j, k, l) \
    for (long int l = 0; l < Nsize(V); l++) \
    for (long int k = 0; k < Zsize(V); k++) \
    for (long int j = 0; j < Ysize(V); j++) \
    for (long int i = 0; i < Xsize(V); i++)

/** For all elements in the array
 *
 * This macro is used to generate loops for the array in an easy manner.
 * It binds the names provided to logical indices which iterate over the data.
 *
 * This macro is used to easily loop through a matrix.
 *
 * @code
 * FOR_ALL_NZYX_ELEMENTS_IN_MULTIDIMARRAY(v, i, j, k, l) {
 *     std::cout << v.elem(i, j, k, l) << " ";
 * }
 * @endcode
 */
#define FOR_ALL_NZYX_ELEMENTS_IN_MULTIDIMARRAY(V, i, j, k, l) \
    for (long int l = 0; l < Nsize(V); l++) \
    for (long int k = Zinit(V); k <= Zlast(V); k++) \
    for (long int j = Yinit(V); j <= Ylast(V); j++) \
    for (long int i = Xinit(V); i <= Xlast(V); i++)

/** For all elements in the array.
 *
 * This macro is used to generate loops for the volume in an easy way.
 * It binds the names provided to logical indices which iterate over the data.
 *
 * @code
 * FOR_ALL_ELEMENTS_IN_ARRAY3D(V, i, j, k) {
 *     std::cout << V.elem(i, j, k) << " ";
 * }
 * @endcode
 */
#define FOR_ALL_ELEMENTS_IN_ARRAY3D(V, i, j, k) \
    for (long int k = Zinit(V); k <= Zlast(V); k++) \
    for (long int j = Yinit(V); j <= Ylast(V); j++) \
    for (long int i = Xinit(V); i <= Xlast(V); i++)

/** For all elements in the array
 *
 * This macro is used to easily loop through a matrix.
 * It binds the names provided to logical indices which iterate over the data.
 *
 * @code
 * FOR_ALL_ELEMENTS_IN_ARRAY2D(m, i, j) {
 *     std::cout << m.elem(i, j) << " ";
 * }
 * @endcode
 */
#define FOR_ALL_ELEMENTS_IN_ARRAY2D(m, i, j) \
    for (long int j = Yinit(m); j <= Ylast(m); j++) \
    for (long int i = Xinit(m); i <= Xlast(m); i++)

/** Direct access
 *
 * Be careful.
 * These functions go against the array library philosophy,
 * and should therefore be avoided.
 * They give physical/direct access into an array
 * without taking logical position into account.
 * Arrays usually follow the C convention of starting index 0.
 *
 * For instance, you can loop over a vector by physical index:
 * @code
 * for (long int i = 0; i < Xsize(v); i++)
 *     foo(direct::elem(v, i));
 * @endcode
 *
 * You can loop over a matrix by physical index:
 * @code
 * for (long int j = 0; j < Ysize(m); j++)
 * for (long int i = 0; i < Xsize(m); i++)
 *     foo(direct::elem(m, i, j));
 * @endcode
 *
 * You can loop over a volume by physical index:
 * @code
 * for (long int k = 0; k < Zsize(V); k++)
 * for (long int j = 0; j < Ysize(V); j++)
 * for (long int i = 0; i < Xsize(V); i++)
 *     foo(direct::elem(m, i, j, k));
 * @endcode
 */
namespace direct {

    // Direct array access
    template <typename T>
    inline const T& elem(const MultidimArray<T> &v, long int i, long int j = 0, long int k = 0, long int l = 0) {
        return v.data[
            i +
            j * v.xdim +
            k * v.xdim * v.ydim +
            l * v.xdim * v.ydim * v.zdim];
    }

    template <typename T>
    inline T& elem(MultidimArray<T> &v, long int i, long int j = 0, long int k = 0, long int l = 0) {
        return v.data[
            i +
            j * v.xdim +
            k * v.xdim * v.ydim +
            l * v.xdim * v.ydim * v.zdim];
    }

};
//@}

/// Xmipp arrays
template<typename T, typename Allocator>
class MultidimArray {

    using  index_t =          long int;
    using uindex_t = unsigned long int;

    Allocator allocator;  // Either allocate memory or map to a file

    public:

    T *data;  // Pointer to heap-allocated memory

    /** Number of elements in X/Y/Z/N
     * Conventions:
     * - The Z dimension splits our array into slices
     * - The N dimension splits our array into images
     */
    uindex_t xdim, ydim, zdim, ndim;

    /// TODO: Manage access to xdim, ydim, zdim, ndim.

    // Array size (number of elements)
    inline unsigned long int size() const { return xdim * ydim * zdim * ndim; }

    // X/Y/Zinit
    index_t xinit, yinit, zinit;

    /** "Iterator" support
     *
     * begin returns a pointer to the start of the data.
     * end returns a pointer just beyond the data.
     * Together, they can be used to iterate over the array:
     *
     * @code
     * for (const auto *ptr = begin(); ptr != end(); ++ptr) {
     *     foo(*ptr);
     * }
     * @endcode
     */
    T* begin() const { return data; }
    T* end() const { return data + size(); }

    // Logical array access
    inline const T& elem(index_t i, index_t j = 0, index_t k = 0, index_t l = 0) const {
        return direct::elem(*this, i - xinit, j - yinit, k - zinit, l);
    }

    inline T& elem(index_t i, index_t j = 0, index_t k = 0, index_t l = 0) {
        return direct::elem(*this, i - xinit, j - yinit, k - zinit, l);
    }

    /// @name Constructors
    //@{

    // Default ctor
    MultidimArray():
    xdim(0), ydim(0), zdim(0), ndim(0), xinit(0), yinit(0), zinit(0),
    data(nullptr), allocator(Allocator()) {}

    /** Size ctor
     * Construct an array (heap-allocate memory) and fill it with zeros.
     */
    MultidimArray(long int Xdim, long int Ydim = 1, long int Zdim = 1, long int Ndim = 1) {
        coreInit();
        resize(Xdim, Ydim, Zdim, Ndim);
    }

    // Copy ctor
    MultidimArray(const MultidimArray<T> &other):
    xdim(other.xdim), ydim(other.ydim), zdim(other.zdim), ndim(other.ndim),
    xinit(other.xinit), yinit(other.yinit), zinit(other.zinit),
    data(nullptr), allocator(Allocator()) {
        resize(other);
        memcpy(data, other.data, sizeof(T) * size());
    }

    // Copy ctor with type cast
    template <typename U>
    MultidimArray(const MultidimArray<U> &other) {
        coreInit();
        *this = other;  //  This is not defined!
    }

    // Move ctor
    MultidimArray(MultidimArray<T> &&other) noexcept:
    xdim(other.xdim), ydim(other.ydim), zdim(other.zdim), ndim(other.ndim),
    xinit(other.xinit), yinit(other.yinit), zinit(other.zinit),
    data(other.data), allocator(Allocator()) { other.data = nullptr; }

    /** Constructor from a Vector.
     * The Size constructor creates an array with memory associated,
     * and fills it with zeros.
     */
    MultidimArray(const Vector<T> &v) {
        coreInit();
        resize(v.size());
        for (long int i = 0; i < v.size(); i++)
            (*this)[i] = v[i];
    }

    /** Constructor from vector
     * This will create a 1D MultidimArray
     * the size and elements will be copied from
     * the std::vector
     */
    MultidimArray(const std::vector<T> &v) {
        coreInit();
        resize(v.size());
        for (long int i = 0; i < v.size(); i++)
            (*this)[i] = v[i];
    }

    // Dtor
    ~MultidimArray() {
        coreDeallocate();
    }

    void clear() {
        coreDeallocate();
        coreInit();
    }
    //@}

    /// @name Core memory operations
    //@{

    /** Core init.
     * Initialize everything to 0
     */
    void coreInit() {
        xdim = ydim = zdim = ndim = 0;
        setOrigin();
        data = nullptr;
        allocator = relion_aligned_mallocator();
    }

    /** Core allocate with dimensions.
     */
    void coreAllocate(uindex_t xdim, uindex_t ydim, uindex_t zdim, uindex_t ndim) {
        if (xdim == 0 || ydim == 0 || zdim == 0 || ndim == 0) {
            clear();
            return;
        }

        this->xdim = xdim;
        this->ydim = ydim;
        this->zdim = zdim;
        this->ndim = ndim;

        coreAllocate();
    }

    /** Core allocate according to current size.
     * Do nothing if data is not nullptr.
     */
    void coreAllocate() {
        if (data) return;  // Memory already allocated
        data = (T*) allocator.allocate(size() * sizeof(T));
    }

    bool getMmap() { return typeid(allocator) == typeid(mmapper); }

    /** Core deallocate.
     * Free all data.
     * Essential to avoid memory leaks.
     */
    void coreDeallocate() {
        if (!data) return;  // Memory not yet allocated
        allocator.deallocate(data, size() * sizeof(T));
        data = nullptr;
    }

    //@}

    /// @name Size
    //@{

    /** Sets new 4D dimensions.
     *
     *  Note that the dataArray is NOT resized. This should be done separately with coreAllocate()
     *
     */
    void setDimensions(index_t xdim = 1, index_t ydim = 1, index_t zdim = 1, index_t ndim = 1) {
        this->xdim = xdim;
        this->ydim = ydim;
        this->zdim = zdim;
        this->ndim = ndim;
    }

    /** NOTE: When setXdim/setYdim/setZdim/setNdim are used, the array is not resized.
     * This should be done separately with coreAllocate()
     */
    void setXdim(index_t xdim) { this->xdim = xdim; }
    void setYdim(index_t ydim) { this->ydim = ydim; }
    void setZdim(index_t zdim) { this->zdim = zdim; }
    void setNdim(index_t ndim) { this->ndim = ndim; }

    /** Copy the shape parameters
     */
    void copyShape(const MultidimArray<T> &other) {
        xdim = other.xdim;
        ydim = other.ydim;
        zdim = other.zdim;
        ndim = other.ndim;
        setOrigin(other.xinit, other.yinit, other.zinit);
    }

    /** Shrink to fit
     *
     * This function resize the memory occupied by the array
     * As a delayed mechanism for mem release
     * This is important, because, shrinking requires additional
     * memory. That will exceed mem limit at peak. This function
     * thus delayed the peak and helps to keep mem usage within
     * limits.
     */
    void shrinkToFit() {
        if (!data || typeid(Allocator) != typeid(relion_aligned_mallocator) || size() == 0)
            return;
        size_t n = sizeof(T) * size();
        T* old_data = data;
        data = (T*) allocator.allocate(n);
        memcpy(data, old_data, n);
        allocator.deallocate(old_data, n);
    }

    /**
     * The functions reshape and shrinkToFit were added
     * on suggestion by Yunxiao Zhang (5 April 2016)
     */

    /** Adjust array to a given shape
     *
     * This function will resize the array to the given size.
     * No data will be copied/moved to the new space.
     * If shape is unchanged, then so is the data.
     * Otherwise, data is almost always destroyed.
     */
    void reshape(uindex_t xdim = 1, uindex_t ydim = 1, uindex_t zdim = 1, uindex_t ndim = 1) {
        if (data && size() == xdim * ydim * zdim * ndim) {
            setDimensions(xdim, ydim, zdim, ndim);
            return;
        }

        if (xdim == 0 || ydim == 0 || zdim == 0 || ndim == 0) {
            clear();
            return;
        }

        coreDeallocate();
        coreAllocate(xdim, ydim, zdim, ndim);
    }

    /** Adjust shape to match the target array
     *
     * No guarantee about the data stored
     */
    template<typename U>
    void reshape(const MultidimArray<U> &other) {
        reshape(other.xdim, other.ydim, other.zdim, other.ndim);
        setOrigin(other.xinit, other.yinit, other.zinit);
    }

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
    void resizeNoCp(uindex_t xdim = 1, uindex_t ydim = 1, uindex_t zdim = 1, uindex_t ndim = 1) {

        const size_t new_size = xdim * ydim * zdim * ndim;
        if (data && size() == new_size)
            return;

        if (new_size == 0) {
            clear();
            return;
        }

        // data can be nullptr even when xdim etc are not zero
        // (This can happen for reading of images...)
        // In that case, initialize data to zeros.
        if (size() > 0 && !data) {
            coreAllocate();
            return;
        }

        // Ask for memory
        T *new_data;

        try {
            new_data = (T*) allocator.allocate(new_size * sizeof(T));
        } catch (std::bad_alloc &) {
            REPORT_ERROR("Allocate: No space left");
        }

        // deallocate old vector
        coreDeallocate();

        // assign *this vector to the newly created
        setDimensions(xdim, ydim, zdim, ndim);
        data = new_data;
    }

    /** Resize to a given size
     *
     * This function resize the actual array to the given size.
     * The origin is not modified.
     * If the actual array is larger than the pattern
     * then the values outside the new size are lost.
     * If it is smaller then 0's are added.
     * If there is no memory, throw an exception.
     *
     * @code
     * V1.resize(3, 3, 2);
     * @endcode
     */
    MultidimArray<T>& resize(uindex_t xdim = 1, uindex_t ydim = 1, uindex_t zdim = 1, uindex_t ndim = 1) {
        size_t new_size = xdim * ydim * zdim * ndim;
        if (data && new_size == size()) {
            setDimensions(xdim, ydim, zdim, ndim);
            return *this;
        }

        if (new_size == 0) {
            clear();
            return *this;
        }

        // data can be nullptr even when dimensions are not zero
        // (This can happen for reading of images...)
        // In that case, initialize data to zeros.
        if (size() > 0 && !data) {
            coreAllocate();
            return *this;
        }

        // Ask for memory
        T *new_data = nullptr;

        try {
            new_data = (T*) allocator.allocate(new_size * sizeof(T));
        } catch (std::bad_alloc &) {
            REPORT_ERROR("Allocate: No space left");
        }

        // Copy needed elements, fill with 0 if necessary
        for (index_t l = 0; l < ndim; l++)
        for (index_t k = 0; k < zdim; k++)
        for (index_t j = 0; j < ydim; j++)
        for (index_t i = 0; i < xdim; i++) {
            // 0 if out of bounds
            new_data[i + j * xdim + k * xdim * ydim + l * xdim * ydim * zdim] = (
                i >= this->xdim || j >= this->ydim || k >= this->zdim
            ) ? (T) 0.0 : direct::elem(*this, i, j, k);
        }

        // deallocate old vector
        coreDeallocate();

        // assign *this vector to the newly created
        data = new_data;
        setDimensions(xdim, ydim, zdim, ndim);
        return *this;
    }

    /** Resize to match
     *
     * This function resizes the array
     * to match the size and origin of another.
     * If the array starts off too large then it is truncated.
     * If it starts off too small then 0's are added at the end
     *
     * @code
     * v2.resize(v1);
     * // Now v2 has the same structure as v1
     * @endcode
     */
    template<typename U>
    void resize(const MultidimArray<U> &other) {
        resize(other.xdim, other.ydim, other.zdim, other.ndim);
        setOrigin(other.xinit, other.yinit, other.zinit);
    }

    /** Return the array's dimensions in X/Y/Z/N,
     * in other words, its shape.
     * It is wise to have the widest dimension first.
     * Redundant dimensions have width 1:
     * (x, 1, 1, 1) or (x, y, 1, 1).
     */
    std::array<uindex_t, 4> getDimensions() const {
        return { xdim, ydim, zdim, ndim };
    }

    std::array<index_t, 3> getOrigin() const {
        return { xinit, yinit, zinit };
    }

    void setOrigin(index_t xinit = 0, index_t yinit = 0, index_t zinit = 0) {
        this->xinit = xinit;
        this->yinit = yinit;
        this->zinit = zinit;
    }

    /** The dimensionality of an array.
     *
     * The number of indices needed to select an element.
     *
     * @code
     * int dim = V.getDim();
     * @endcode
     */
    inline int getDim() const {
        return zdim > 1 ? 3 : ydim > 1 ? 2 : xdim > 1 ? 1 : 0;
    }

    /** Check dimension.
     *
     * Is this the array's dimension?
     * Print an error message if not.
     */
    #define checkDimension(dim) checkDimensionWithDebug(dim, __FILE__, __LINE__)
    void checkDimensionWithDebug(int dim, const char *file, int line) const {
        if (getDim() != dim) {
            std::cerr << " Check for dimension: " << dim << std::endl;
            std::cerr << "MultidimArray shape: ";
            printShape(std::cerr);
            std::cerr << std::endl;
            std::cerr << "Check called from file " << file << " line " << line << std::endl;
            exit(1);
        }
    }

    /** Put a 3D window to the nth volume
     *
     * The volume is windowed within the two positions given to this function.
     * Indexes always refer to logical indices. If a position is outside the
     * actual matrix range then the matrix is padded init_value until the
     * new position is reached. In the following example suppose that m1
     * is the following and that the origin is (-1,-1,-1).
     *
     * @code
     * slice 0
     * [01 02 03          [
     *  04 05 06           04 05 06 0
     *  07 08 09]          07 08 09 0]
     *
     * ----->
     *
     * slice 1
     * [11 12 13          [
     *  14 15 16           14 15 16 0
     *  17 18 19]          17 18 19 0]
     * @endcode
     *
     * @code
     * V1.window(0, 0, -1, 1, 1, 2);
     * @endcode
     */
    MultidimArray<T> windowed(
        index_t x0, index_t xF, index_t y0, index_t yF, index_t z0, index_t zF,
        T init_value = 0, long n = 0
    ) const {
        MultidimArray<T> result (xF - x0 + 1, yF - y0 + 1, zF - z0 + 1);
        result.setOrigin(x0, y0, z0);
        for (index_t k = z0; k <= zF; k++)
        for (index_t j = y0; j <= yF; j++)
        for (index_t i = x0; i <= xF; i++) {
            result.elem(i, j, k) = inside(i, j, k) ? elem(i, j, k, n) : init_value;
        }
        return result;
    }

    /** Put a 2D window to the nth matrix
     *
     * The matrix is windowed within the two positions given to this function.
     * Indexes always refer to logical indices. If a position is outside the
     * actual matrix range then the matrix is padded with init_value until the
     * new position is reached. In the following examples suppose that m1 is the
     * following and that the origin is (-1,-1).
     *
     * @code
     *      [1 2 3               [1 2 3 0
     * m1 =  4 5 6    --->  m1 =  4 5 6 0
     *       7 8 9]               7 8 9 0]
     *
     * @endcode
     *
     * @code
     * m1.window(-1, -1, 1, 2);
     * @endcode
     */
    MultidimArray<T> windowed(
        index_t x0, index_t xF, index_t yF, index_t y0,
        T init_value = 0, long n = 0
    ) const {
        MultidimArray<T> result (xF - x0 + 1, yF - y0 + 1);
        result.setOrigin(x0, y0);
        FOR_ALL_ELEMENTS_IN_ARRAY2D(result, i, j) {
            result.elem(i, j) = inside(i, j) ? elem(i, j, 0, n) : init_value;
        }
        return result;
    }

    /** Put a 1D window to the nth vector
     *
     * The vector is windowed within the two indices given to this function.
     * Indexes always refer to logical indices. If an index is outside the
     * actual vector range then the vector is padded winit_value. In the
     * following examples suppose that v1=[-2 -1 0 1 2] and that the origin is
     * -2.
     *
     * @code
     * v1.window(-1, 2);  // v1=[-1 0 1 2]; v1.firstX() == -1
     *
     * v1.window(-3, 1);  // v1=[0 -2 -1 0 1]; v1.firstX() == -3
     * @endcode
     */
    MultidimArray<T> windowed(
        index_t x0, index_t xF, T init_value = 0, long n = 0
    ) const {
        MultidimArray<T> result (xF - x0 + 1);
        result.setOrigin(x0);
        for (index_t i = x0; i <= xF; i++) {
            result.elem(i) = inside(i) ? elem(i, 0, 0, n) : init_value;
        }
        return result;
    }

    /** Print shape of multidimensional array.
     *
     * This function shows the size, starting and finishing indices of the
     * given array.
     * No end of line is printed at either the beginning or the end.
     *
     * @code
     * v.printShape();
     *
     * std::ofstream fh;
     * ...;
     * v.printShape(fh);
     * @endcode
     */
    void printShape(std::ostream &out = std::cout) const {
        if (ndim > 1)
            out << " Number of images = " << ndim;

        switch (getDim()) {

            case 3:
            out << " Size(Z,Y,X): " << zdim << "×" << ydim << "×" << xdim
            << " k=[" << firstZ() << ".." << lastZ() << "]"
            << " i=[" << firstY() << ".." << lastY() << "]"
            << " j=[" << firstX() << ".." << lastX() << "]";
            break;

            case 2:
            out << " Size(Y,X): " << ydim << "×" << xdim
            << " i=[" << firstY() << ".." << lastY() << "]"
            << " j=[" << firstX() << ".." << lastX() << "]";
            break;

            case 1:
            out << " Size(X): " << xdim
            << " j=[" << firstX() << ".." << lastX() << "]";
            break;

            default:
            out << " Empty MultidimArray!";

        }
        out << "\n";
    }

    /** sameShape
     *
     * Do these two arrays have the same shape (dimensions and origin)?
     */
    template <typename U>
    inline bool sameShape(const MultidimArray<U> &other) const {
        return getDimensions() == other.getDimensions() && // Same dimensions
               getOrigin()     == other.getOrigin();       // Same origin
    }

    /** inside for 3D matrices */
    inline bool inside(index_t i, index_t j, index_t k) const {
        return inXbounds(i, *this) && inYbounds(j, *this) && inZbounds(k, *this);
    }

    /** inside for 2D matrices */
    inline bool inside(index_t i, index_t j) const {
        return inXbounds(i, *this) && inYbounds(j, *this);
    }

    /** inside for 1D matrices */
    inline bool inside(index_t i) const {
        return inXbounds(i, *this);
    }

    /** outside for 3D matrices */
    bool outside(index_t i, index_t j, index_t k) const {
        return !inside(i, j, k);
    }

    /** outside for 2D matrices */
    bool outside(index_t i, index_t j) const {
        return !inside(i, j);
    }

    /** outside for 1D matrices */
    bool outside(index_t i) const {
        return !inside(i);
    }

    /** Is this logical index vector outside the bounds of this array? */
    bool outside(const Vector<RFLOAT> &r) const {
        int rsize = r.size();
        if (rsize < 1) REPORT_ERROR(std::string(__func__) + ": index vector has too few components");

        switch (rsize) {
            case 1:
            return outside(r[0]);
            case 2:
            return outside(r[0], r[1]);
            case 3:
            return outside(r[0], r[1], r[2]);
            default:
            REPORT_ERROR(std::string(__func__) + ": index vector has too many components");
        }
    }

    /** Return Y dimension. */
    inline index_t rowNumber() const { return ydim; }

    /** Return X dimension. */
    inline index_t colNumber() const { return xdim; }

    // Wouldn't it be nice to have this in the Xmipp namespace? As in:
    // Xmipp::setOrigin(arr);

    /** Set logical origin in Xmipp fashion.
     *
     * Adjust the starting points in the array
     * so that the center of the array is defined in the Xmipp fashion.
     *
     * @code
     * V.setXmippOrigin();
     * @endcode
     */
    MultidimArray<T>& setXmippOrigin() {
        xinit = Xmipp::init(xdim);
        yinit = Xmipp::init(ydim);
        zinit = Xmipp::init(zdim);
        return *this;
    }

    // First logical X index
    inline index_t firstX() const { return Xmipp::init(xdim); }

    // Last logical X index
    inline index_t  lastX() const { return Xmipp::last(xdim); }

    // First logical Y index
    inline index_t firstY() const { return Xmipp::init(ydim); }

    // Last logical Y index
    inline index_t  lastY() const { return Xmipp::last(ydim); }

    // First logical Z index
    inline index_t firstZ() const { return Xmipp::init(zdim); }

    // Last logical Z index
    inline index_t  lastZ() const { return Xmipp::last(zdim); }

    /** IsCorner (in 2D or 3D matrix)
     *
     * Is this logical index a corner of this array?
     */
    bool isCorner(const Vector<RFLOAT> &v) const {

        if (v.size() < 2)
            REPORT_ERROR(std::string(__func__) + ": index vector has too few components");

        switch (xdim) {

            case 2:
            return (v[0] == firstX() || v[0] == lastX()) &&
                   (v[1] == firstY() || v[1] == lastY());

            case 3:
            return (v[0] == firstX() || v[0] == lastX()) &&
                   (v[1] == firstY() || v[1] == lastY()) &&
                   (v[2] == firstZ() || v[2] == lastZ());

            default:
            REPORT_ERROR(std::string(__func__) + ": index vector has too many components");

        }

    }
    //@}

    ///@name Access to the pixel values
    //@{

    T& operator [] (index_t n) const { return data[n]; }

    /** Volume element access by integer vector.
     *
     * Returns the value of a matrix logical position, but this time the
     * element position is determined by a R3 vector. The elements can be used
     * either by value or by reference. An exception is thrown if the index is
     * outside the logical range. Pay attention in the following example that
     * we are accessing the same element as in the previous function but, now
     * we have to give first the X position because we are building first a
     * vector of the form (x,y,z).
     *
     * @code
     * V(vectorR3(1, -2, 0)) = 1;
     * val = V(vectorR3(1, -2, 0));
     * @endcode
     */
    T& operator()(const Vector<index_t> &v) {
        switch (v.size()) {
            case 1:
            return elem(v[0]);
            case 2:
            return elem(v[0], v[1]);
            case 3:
            return elem(v[0], v[1], v[2]);
            default:
            REPORT_ERROR("Matrix dimensions must be 1, 2, or 3");
        }
    }

    /** 3D element access by index (getVoxel).
    *
    * Same function as operator() but with a name. Needed by swig.
    *
    */
    inline T getVoxel(index_t k, index_t i, index_t j) const {
        return elem(i, j, k);
    }

    /** 3D element access by index (setVoxel).
    *
    * Same function as operator() but with a name. Needed by swig.
    *
    */
    inline void setVoxel(index_t k, index_t i, index_t j, T newval) {
        elem(i, j, k) = newval;
    }

    /** Matrix element access by index
     *
     * Returns the value of a matrix logical position. In our example we could
     * access from v(-2,-1) to v(2,1). The elements can be used either by value
     * or by reference. An exception is thrown if the index is outside the
     * logical range. The first argument is the Y position and the second the X
     * position.
     *
     * @code
     * m(-2, 1) = 1;
     * val = m(-2, 1);
     * @endcode
     */
    inline const T& operator()(index_t i, index_t j) const {
        return elem(i, j);
    }

    /** Vector element access
     *
     * Returns the value of a vector logical position.
     * In our example we could access from v(-2) to v(2).
     * The elements can be used either by value or by reference.
     * An exception is thrown if the index is outside the logical range.
     *
     * @code
     * v(-2) = 1;
     * val = v(-2);
     * @endcode
     */
    inline const T& operator()(index_t i) const {
        return elem(i);
    }

    /** Get a single 1,2 or 3D image from a multi-image array
     *
     * This function extracts a single-image array from a multi-image one.
     * @code
     * V.getImage(0, m);
     * @endcode
     */
    void getImage(long n, MultidimArray<T> &M) const {
        if (xdim == 0) {
            M.clear();
            return;
        }

        if (n > ndim)
            REPORT_ERROR(" Multidimarray getImage: n larger than ndim (out of bounds)");

        M.resize(xdim, ydim, zdim);
        for (index_t k = 0; k < Zsize(M); k++)
        for (index_t j = 0; j < Ysize(M); j++)
        for (index_t i = 0; i < Xsize(M); i++) {
            direct::elem(M, i, j) = direct::elem(*this, i, j, k, n);
        }

        M.setOrigin(firstX(), firstY(), firstZ());
    }

    /** Set a single 1,2 or 3D image in a multi-image array
     *
     * Set a single-image array in a multi-image one.
     * @code
     * V.setImage(0, m);
     * @endcode
     */
    void setImage(long n, const MultidimArray<T> &M) {
        if (xdim == 0)
            return;

        if (n < 0 || n > ndim)
            REPORT_ERROR("setImage: MultidimArray subscript (n) out of range");

        if (M.zdim != zdim || M.ydim != ydim || M.xdim != xdim)
            REPORT_ERROR("setImage: MultidimArray dimensions different from the input image ones");

        for (index_t k = 0; k < Zsize(M); k++)
        for (index_t j = 0; j < Ysize(M); j++)
        for (index_t i = 0; i < Xsize(M); i++)
            direct::elem(*this, i, j, k, n) = direct::elem(M, i, j, k);

    }

    /** 2D Slice access for reading.
     *
     * Return a slice (a 2D matrix) corresponding to the choosen
     * slice inside the nth 3D matrix, the numbering of the slices is also logical not
     * physical. This function differs from the previous one in that this one
     * cuts and assign in a single step instead of in two steps, as in
     * the previous example.
     *
     * @code
     * V.slice(0, m);
     * @endcode
     */
    void getSlice(index_t k, MultidimArray<T> &M, char axis = 'Z', long n = 0) const {
        if (xdim == 0) {
            M.clear();
            return;
        }

        switch (axis) {

            case 'Z':
            if (!inZbounds(k, (*this)))
                REPORT_ERROR(std::string(__func__) + ": Multidim index out of range");

            k -= firstZ();
            M.resize(xdim, ydim);
            for (index_t j = 0; j < Ysize(M); j++)
            for (index_t i = 0; i < Xsize(M); i++) {
                direct::elem(M, i, j) = direct::elem(*this, j, i, k, n);
            }
            M.setOrigin(firstX(), firstY());
            break;

            case 'Y':
            if (!inYbounds(k, (*this)))
                REPORT_ERROR(std::string(__func__) + ": Multidim index out of range");

            k -= firstY();
            M.resize(xdim, zdim);
            for (index_t j = 0; j < Ysize(M); j++)
            for (index_t i = 0; i < Xsize(M); i++) {
                direct::elem(M, i, j) = direct::elem(*this, k, j, i, n);
            }
            M.setOrigin(firstX(), firstY());
            break;

            case 'X':
            if (!inXbounds(k, (*this)))
                REPORT_ERROR(std::string(__func__) + ": Multidim index out of range");

            k -= firstX();
            M.resize(ydim, zdim);
            for (index_t j = 0; j < Ysize(M); j++)
            for (index_t i = 0; i < Xsize(M); i++) {
                direct::elem(M, i, j) = direct::elem(*this, j, k, i, n);
            }
            M.setOrigin(firstY(), firstZ());
            break;

            default:
            REPORT_ERROR(std::string(__func__) + ": unsupported axis " + axis);

        }
    }

    /** Slice access for writing.
     *
     * Set a 2D matrix corresponding to the choosen slice inside the nth
     * volume, the numbering of the slices is logical not physical.
     *
     * @code
     * // Copies slice 0 in slice 1
     * V.setSlice(1, (V.slice(0)));
     * @endcode
     */
    void setSlice(index_t k, const MultidimArray<T> &v, long n = 0) {
        if (xdim == 0)
            return;

        if (k < firstZ() || k > lastZ())
            REPORT_ERROR(std::string(__func__) + ": MultidimArray subscript (k) out of range");

        if (v.ydim != ydim || v.xdim != xdim)
            REPORT_ERROR(std::string(__func__) + ": MultidimArray dimensions different from the matrix ones");

        k -= firstZ();

        for (index_t j = 0; j < Ysize(v); j++)
        for (index_t i = 0; i < Xsize(v); i++)
            direct::elem(*this, i, j, k, n) = direct::elem(v, i, j);
    }

    /** Get Column
     *
     * Return a column vector corresponding to the choosen column.
     *
     * @code
     * std::vector<RFLOAT> v;
     * m.getCol(-1, v);
     * @endcode
     */
    void getCol(index_t j, MultidimArray<T> &v) const {
            if (xdim == 0 || ydim == 0) {
                v.clear();
                return;
        }

        if (j < 0 || j >= xdim)
            REPORT_ERROR("getCol: Matrix subscript (j) greater than matrix dimension");

        v.resize(ydim);
        for (index_t i = 0; i < ydim; i++)
            v.elem(i) = elem(i, j);
    }

    /** Set Column
     *
     * Set a column vector corresponding to the choosen column
     * inside matrix.
     *
     * @code
     * m.setCol(0, (m.row(1)).transpose());  // Copies row 1 in column 0
     * @endcode
     */
    void setCol(index_t j, const MultidimArray<T> &v) {
        if (xdim == 0 || ydim == 0)
            REPORT_ERROR("setCol: Target matrix is empty");

        if (j < 0 || j>= xdim)
            REPORT_ERROR("setCol: Matrix subscript (j) out of range");

        if (v.xdim != ydim)
            REPORT_ERROR("setCol: Vector dimension different from matrix one");

        for (index_t i = 0; i < ydim; i++)
            elem(i, j) = v.elem(i);
    }

    /** Get row
     *
     * Return a row vector corresponding to the choosen row inside the nth 2D matrix.
     * The numbering of the rows is logical not physical.
     *
     * @code
     * std::vector<RFLOAT> v;
     * m.getRow(-2, v);
    * @endcode
    */
    void getRow(index_t i, MultidimArray<T> &v) const {
        if (xdim == 0 || ydim == 0) {
            v.clear();
            return;
        }

        if (i < 0 || i >= ydim)
            REPORT_ERROR("getRow: Matrix subscript (i) greater than matrix dimension");

        v.resize(xdim);
        for (index_t j = 0; j < xdim; j++)
            v.elem(j) = elem(i, j);
    }

    /** Set Row
     *
     * Set a row vector corresponding to the choosen row in the 2D Matrix
     *
     * @code
     * m.setRow(-2, m.row(1));  // Copies row 1 in row -2
     * @endcode
     */
    void setRow(index_t i, const MultidimArray<T> &v) {
        if (xdim == 0 || ydim == 0)
            REPORT_ERROR("setRow: Target matrix is empty");

        if (i < 0 || i >= ydim)
            REPORT_ERROR("setRow: Matrix subscript (i) out of range");

        if (v.xdim != xdim)
            REPORT_ERROR("setRow: Vector dimension different from matrix one");

        for (index_t j = 0; j < xdim; j++)
            elem(i, j) = v.elem(j);
    }

    // Logical to physical index translation (3D)
    std::array<long int, 3> toPhysical(index_t i, index_t j, index_t k) const {
        return {i - firstX(), j - firstY(), k - firstZ()};
    }

    // Logical to physical index translation (2D)
    std::array<long int, 2> toPhysical(index_t i, index_t j) const {
        return {i - firstX(), j - firstY()};
    }

    // Logical to physical index translation (1D)
    std::array<long int, 1> toPhysical(index_t i) const {
        return {i - firstX()};
    }

    // Physical to logical index translation (3D)
    std::array<long int, 3> toLogical(index_t i, index_t j, index_t k) const {
        return {i + firstX(), j + firstY(), k + firstZ()};
    }

    // Physical to logical index translation (2D)
    std::array<long int, 2> toLogical(index_t i, index_t j) const {
        return {i + firstX(), j + firstY()};
    }

    // Physical to logical index translation (1D)
    std::array<long int, 1> toLogical(index_t i) const {
        return {i + firstX()};
    }

    //@}

    //@}

    /** @name Operations between two arrays
     *
     * These operations are performed between two arrays of the same type.
     * (typeCast can be used to convert an array to the appropriate type.)
     * The result will be of the same type as the operands.
     *
     * In this kind of operations each component of the lhs array is reassigned
     * according to the operation in question and its corresponding component in the rhs array.
     * It is important that both have got the same size and shape.
     * The result will have the same size and shape as the two operands.
     */
    //@{

    MultidimArray<T>& operator += (const MultidimArray<T> &arg);

    MultidimArray<T>& operator -= (const MultidimArray<T> &arg);

    MultidimArray<T>& operator *= (const MultidimArray<T> &arg);

    MultidimArray<T>& operator /= (const MultidimArray<T> &arg);

    //@}

    /// @name Initialization
    /// @{

    /** Assign from scalar
     *
     * The scalar must be of the same type as the array type.
     * You cannot assign a float to an int array without a casting.
     *
     * @code
     * MultidimArray<float> v (64, 64, 64);
     * v = 3.14;
     * @endcode
     */
    MultidimArray<T>& operator = (T scalar) {
        for (auto &x : *this) { x = scalar; }
        return *this;
    }

    /** Initialize to zeros with current size.
     *
     * All values are set to 0. The current size and origin are kept. It is not
     * an error if the array is empty, then nothing is done.
     *
     * @code
     * v.initZeros();
     * @endcode
     */
    inline void initZeros() {
        memset(data, 0, size() * sizeof(T));
    }

    template <typename U>
    inline void fill(U u) {
        for (auto& t: *this) t = u;
    }

    /** Initialize to zeros with a given size.
     */
    inline void initZeros(long int Xdim, long int Ydim = 1, long int Zdim = 1, long int Ndim = 1) {
        if (xdim != Xdim || ydim != Ydim || zdim != Zdim || ndim != Ndim)
            reshape(Xdim, Ydim, Zdim, Ndim);
        initZeros();
    }

    template <typename U>
    static MultidimArray<T> zeros(const MultidimArray<U> &other) {
        MultidimArray<T> self (other);
        self.initZeros();
        return self;
    }

    static MultidimArray<T> zeros(long int Xdim, long int Ydim = 1, long int Zdim = 1, long int Ndim = 1) {
        MultidimArray<T> arr (Xdim, Ydim, Zdim, Ndim);
        arr.initZeros();
        return arr;
    }

    static MultidimArray<T> ones(long int Xdim, long int Ydim = 1, long int Zdim = 1, long int Ndim = 1) {
        MultidimArray<T> arr (Xdim, Ydim, Zdim, Ndim);
        memset(arr.data, 1, arr.size() * sizeof(T));
        return arr;
    }

    /** Linear initialization (only for 1D)
     *
     * The 1D vector is filled with values increasing/decreasing linearly within a
     * range or at given steps.
     *
     * Increment functionality: The default increment is 1, the initial point is
     * incremented by this value until the upper limit is reached. This is the
     * default working mode for the function.
     *
     * @code
     * v1.initLinear(1, 3);  // v1=[1 2 3]
     * v1.initLinear(1.5, 3.1);  // v1=[1.5 2.5]
     * v1.initLinear(0, 10, 3);  // v1=[0 3 6 9]
     * v1.initLinear(0, 10, 3, "incr");  // v1=[0 3 6 9]
     * @endcode
     *
     * Step functionality: The given range is divided in as many points as
     * indicated (in the example 6 points).
     *
     * @code
     * v1.initLinear(0, 10, 6, "steps");  // v1=[0 2 4 6 8 10]
     * @endcode
     */
    void initLinear(T minF, T maxF, int n = 1, const std::string& mode = "incr") {
        RFLOAT slope;
        int steps;

        if (mode == "incr") {
            steps = 1 + floor((RFLOAT) abs(maxF - minF) / (RFLOAT) n);
            slope = n * sgn(maxF - minF);  // maxF and minF should not be equal
        } else if (mode == "steps") {
            steps = n;
            slope = (maxF - minF) / (steps - 1);
        } else {
            REPORT_ERROR("Init_linear: Mode not supported (" + mode + ")");
        }

        if (steps == 0) {
            clear();
        } else {
            reshape(steps);
            for (int i = 0; i < steps; i++) {
                elem(i) = (T) ((RFLOAT) minF + slope * i);
            }
        }
    }

    /** Initialize with random values.
     *
     * This function allows you to initialize the array with a set of random
     * values picked from a uniform random distribution or a gaussian one. You
     * must choose two parameters for each, for the uniform distribution they
     * mean the range where to generate the random numbers, while in the
     * gaussian case they are the mean and the standard deviation. By default
     * the uniform distribution is selected. The size and origin of the array
     * are not modified.
     *
     * @code
     * v.randomize(0, 1, "uniform");
     * // uniform distribution between 0 and 1
     *
     * v.randomize(0, 1, "gaussian");
     * // gaussian distribution with mean 0 and stddev 1
     * @endcode
     */
    void randomize(RFLOAT op1, RFLOAT op2, const std::string &mode = "uniform") {
        if (mode == "uniform") {
            for (T *ptr = begin(); ptr != end(); ++ptr)
                *ptr = static_cast<T>(rnd_unif(op1, op2));
        } else if (mode == "gaussian") {
            for (T *ptr = begin(); ptr != end(); ++ptr)
                *ptr = static_cast<T>(rnd_gaus(op1, op2));
        } else {
            REPORT_ERROR(static_cast<std::string>("InitRandom: Mode not supported (" + mode + ")"));
        }
    }

    /** Add noise to actual values.
     *
     * This function add some noise to the actual values of the array according
     * to a certain random distribution. You must choose two parameters for
     * each, for the uniform distribution they mean the range where to generate
     * the random numbers, while in the gaussian case they are the mean and the
     * standard deviation. By default the uniform distribution is selected. The
     * size and origin of the array are not modified. The array itself is
     * modified.
     *
     * @code
     * v1.addNoise(0, 1);
     * // uniform distribution between 0 and 1
     *
     * v1.addNoise(0, 1, "uniform");
     * // the same
     *
     * v1.addNoise(0, 1, "gaussian");
     * // gaussian distribution with 0 mean and stddev=1
     *
     * v1.addNoise(0, 1, "student", 3);
     * // t-student distribution with 0 mean and stddev=1, and 3 degrees of freedom
     *

    * @endcode
    */
    void addNoise(
        RFLOAT op1, RFLOAT op2,
        const std::string &mode = "uniform",
        RFLOAT df = 3.0
    ) const {
        if (mode == "uniform") {
            for (T *ptr = begin(); ptr != end(); ++ptr)
                *ptr += static_cast<T>(rnd_unif(op1, op2));
        } else if (mode == "gaussian") {
            for (T *ptr = begin(); ptr != end(); ++ptr)
                *ptr += static_cast<T>(rnd_gaus(op1, op2));
        } else if (mode == "student") {
            for (T *ptr = begin(); ptr != end(); ++ptr)
                *ptr += static_cast<T>(rnd_student_t(df, op1, op2));
        } else {
            REPORT_ERROR(static_cast< std::string >("AddNoise: Mode not supported (" + mode + ")"));
        }
    }
    //@}

    /** @name Utilities
     *
     * Here you have several easy functions to manage the values of
     * the array.
     */
    //@{

    /** Computes the center of mass of the nth slice
     */
    void centerOfMass(Vector<RFLOAT> &CoM, MultidimArray<int> *mask = nullptr, long int n = 0) {

            CoM.resize(3);
            std::fill(CoM.begin(), CoM.end(), 0);
            RFLOAT mass = 0;
            FOR_ALL_ELEMENTS_IN_ARRAY3D(*this, i, j, k) {
                if ((!mask || mask->elem(i, j, k, n)) && elem(i, j, k) > 0) {
                XX(CoM) += i * elem(i, j, k, n);
                YY(CoM) += j * elem(i, j, k, n);
                ZZ(CoM) += k * elem(i, j, k, n);

                mass += elem(i, j, k, n);
            }
        }

        if (mass != 0) { CoM /= mass; }

        // Resize CoM to the correct dimensionality
        CoM.resize(getDim());

    }

    /** Get the indices that sort the 1D vector elements (original array intact)
     *
     * @code
     * MultidimArray<long int> idx(v1.size());
     * v1.sorted_index(idx);
     * @endcode
     */
    // TODO: check this function!
    void sorted_index(MultidimArray<long> &idx) const {
        checkDimension(1);
        std::vector<std::pair<T, long int>> pairs;
        pairs.reserve(xdim);
        for (long int n = 0; n < size(); n++) {
            pairs.emplace_back((*this)[n], n);
        }
        // Sort the pairs (order determined by first element)
        std::sort(pairs.begin(), pairs.end());
        idx.resize(xdim);
        // Fill the output array with the second elements of the sorted pairs
        for (long int n = 0; n < idx.size(); n++) {
            idx[n] = pairs[n].second;
        }
    }

    /** Several thresholding.
     *
     * Apply a threshold to the array, the object is modified itself. There
     * are several kinds of thresholding and you must specify it, the values
     * given in the fuction have different meanings according to the threshold
     * applied.
     *
     * abs_above: if |x| > a => x = b
     * abs_below: if |x| < a => x = b
     * above:     if  x  > a => x = b
     * below:     if  x  < a => x = b
     * range:     if  x  < a => x = a   and    if x > b => x = b
     *
     * @code
     * v.threshold("abs_above", 10, 10);
     * // any value whose absolute value is above 10 will be substituted by
     * // -10 (if it is negative) or 10 (if it is positive)
     *
     * v.threshold("abs_below", 0.1, 0);
     * // any value whose absolute value is below 0.1 will be substituted by
     * // -0 (if it is negative) or 0 (if it is positive)
     *
     * v.threshold("above", 10, 10);
     * // any value above 10 will be substituted by 10
     *
     * v.threshold("below", -10, -10);
     * // any value below -10 will be substituted by -10
     *
     * v.threshold("range", 0, 1);
     * // v is "saturated" by values 0 and 1, any value outside this range
     * // will be substituted by its nearest border
     * @endcode
     */
    void threshold(const std::string &type, T a, T b, MultidimArray<int> *mask = nullptr);

    /** Count with threshold.
     *
     * Return the number of elements meeting the threshold
     * condition.
     */
    long int countThreshold(const std::string& type, T a, T b, MultidimArray<int> *mask = nullptr) {
        int mode =
            type == "abs_above" ? 1 :
            type == "abs_below" ? 2 :
            type == "above"     ? 3 :
            type == "below"     ? 4 :
            type == "range"     ? 5 :
                                  0;
        if (!mode) REPORT_ERROR(static_cast<std::string>("CountThreshold: mode not supported (" + type + ")"));

        long int ret = 0;
        int *maskptr = mask ? mask->begin() : nullptr;
        for (T *ptr = begin(); ptr != end(); ++ptr, ++maskptr) {
            if (!mask || maskptr > 0) {
                switch (mode) {

                    case 1: if (abs(*ptr) > a) { ret++; } break;

                    case 2: if (abs(*ptr) < a) { ret++; } break;

                    case 3: if (*ptr > a) { ret++; } break;

                    case 4: if (*ptr < a) { ret++; } break;

                    case 5: if (*ptr >= a && *ptr <= b) { ret++; } break;

                }
            }
        }
        return ret;
    }

    /** Substitute a value by another.
     *
     * Substitute an old value by a new one. The accuracy is used to say if
     * the value in the array is equal to the old value. Set it to 0 for
     * perfect accuracy.
     */
    void substitute(
        T oldv, T newv, RFLOAT accuracy = Xmipp::epsilon<RFLOAT>(),
        MultidimArray<int> *mask = nullptr
    ) {
        int *maskptr = mask ? mask->begin() : nullptr;
        for (T *ptr = begin(); ptr != end(); ++ptr, ++maskptr) {
            if (
                (!mask || *maskptr > 0)
                && abs(*ptr - oldv) <= accuracy
            ) { *ptr = newv; }
        }
    }

    /** Substitute a given value by a sample from a Gaussian distribution.
     *
     * Substitute  a given value by a sample from a Gaussian distribution.
     * The accuracy is used to say if the value in the array is equal
     * to the old value.  Set it to 0 for perfect accuracy.
     */
    void randomSubstitute(
        T oldv, T avgv, T sigv, RFLOAT accuracy = Xmipp::epsilon<RFLOAT>(),
        MultidimArray<int> *mask = nullptr
    ) {
        int *maskptr = mask ? mask->begin() : nullptr;
        for (T *ptr = begin(); ptr != end(); ++ptr, ++maskptr) {
            if (
                (!mask || *maskptr > 0)
                && abs(*ptr - oldv) <= accuracy
            ) { *ptr = rnd_gaus(avgv, sigv); }
        }
    }

    /** Binarize.
     *
     * Binarize (set to 1 or 0) each value in a volume,
     * according to whether it is greater than val + accuracy.
     * Use threshold to get a very powerful binarization.
     */
    void binarize(
        RFLOAT val = 0, RFLOAT accuracy = Xmipp::epsilon<RFLOAT>(),
        MultidimArray<int> *mask = nullptr
    ) {
        int *maskptr = mask ? mask->begin() : nullptr;
        for (T *ptr = begin(); ptr != end(); ++ptr, ++maskptr) {
            if (!mask || *maskptr > 0) {
                *ptr = *ptr > val + accuracy;
            }
        }
    }

    #if defined(__APPLE__)
    #undef MIN
    #undef MAX
    #endif

    /** MAX
     *
     * Each component of the result is the maximum of the correspoing
     * components of the two input arrays. They must have the same shape, if
     * not an exception is thrown
     */
    friend void MAX(
        const MultidimArray<T> &v1, const MultidimArray<T> &v2,
        MultidimArray<T> &result
    ) {
        if (!v1.sameShape(v2))
            REPORT_ERROR("MAX: arrays of different shape");

        result.resize(v1);
        for (long int n = 0; n < result.size(); n++)
            result[n] = std::max(v1[n], v2[n]);
    }

    /** MIN
     *
     * Each component of the result is the minimum of the correspoing
     * components of the two input arrays. They must have the same shape, if
     * not an exception is thrown
     */
    friend void MIN(
        const MultidimArray<T> &v1, const MultidimArray<T> &v2,
        MultidimArray<T> &result
    ) {
        if (!v1.sameShape(v2))
            REPORT_ERROR("MIN: arrays of different shape");

        result.resize(v1);
        for (long int n = 0; n < result.size(); n++)
            result[n] = std::min(v1[n], v2[n]);
    }

    /** Sqrt.
     *
     * Each component of the result is the square root of the original
     * component.
     */
    void selfSQRT() {
        for (T *ptr = begin(); ptr != end(); ++ptr)
        *ptr = static_cast<T>(sqrt(static_cast<RFLOAT>(*ptr)));
    }

        /** Sum of matrix values.
         *
         * Return the sum of all internal values.
         *
         * @code
         * RFLOAT sum = m.sum();
         * @endcode
     */
    RFLOAT sum() const {
        RFLOAT sum = 0;
        for (T *ptr = begin(); ptr != end(); ++ptr)
        sum += *ptr;
        return sum;
    }

    /** Sum of squared vector values.
     *
     * Return the sum of all internal values to the second
     * power_class.
     *
     * @code
     * RFLOAT sum2 = m.sum2();
     * @endcode
     */
    RFLOAT sum2() const {
        RFLOAT sum = 0;
        for (T *ptr = begin(); ptr != end(); ++ptr)
        sum += *ptr * *ptr;
        return sum;
    }

    /** Log10.
     *
     * Each component of the result is the log10 of the original components.
     */
    void selfLog10() {
        for (T *ptr = begin(); ptr != end(); ++ptr)
        *ptr = static_cast<T>(log10(static_cast<RFLOAT>(*ptr)));
    }

    /** Reverse matrix values over X axis, keep in this object.
     *
     * Maybe better with an example:
     *
     * @code
     * slice 0
     * [01 02 03          [07 08 09
     *  04 05 06           04 05 06
     *  07 08 09]          01 02 03]
     *
     * ----->
     *
     * slice 1
     * [11 12 13          [17 18 19
     *  14 15 16           14 15 16
     *  17 18 19]          11 12 13]
     * @endcode
     *
     */
    void selfReverseX() {
        long int xsize = xdim;
        long int halfSizeX = (long int) (xsize - 1) / 2;
        long int ysize = ydim;
        long int zsize = zdim;
        long int nsize = ndim;
        //0 column should be handled in a different way
        //for even and odd matrices
        long int start_x = !(xsize % 2);

        for (long int l = 0; l < nsize; l++)
        for (long int k = 0; k < zsize; k++)
        for (long int j = 0; j < ysize; j++)
        for (long int i = start_x; i <=  halfSizeX; i++) {
            std::swap(
                direct::elem(*this, i,         j, k, l),
                direct::elem(*this, xsize - i, j, k, l)
            );
        }
        //NOTE: line x=0 should not be modified since gets itself by wrapping
        //NOTE center hyper-plane  (dimx/2,y,z)  should remain unchanged

    }

    /** Reverse matrix values over Y axis, keep in this object.
     *
     * Maybe better with an example:
     *
     * @code
     * slice 0
     * [01 02 03          [03 02 01
     *  04 05 06           06 05 04
     *  07 08 09]          09 08 07]
     *
     * ----->
     *
     * slice 1
     * [11 12 13          [13 12 11
     *  14 15 16           16 15 14
     *  17 18 19]          19 18 17]
     * @endcode
     *
     */
    void selfReverseY() {
        long int xsize = xdim;
        long int ysize = ydim;
        long int halfSizeY = (long int) (ysize - 1) / 2;
        long int zsize = zdim;
        long int nsize = ndim;
        //0 column should be handled in a different way
        //for even and odd matrices
        long int start_y = !(ysize % 2);

        for (long int l = 0; l < nsize; l++)
        for (long int k = 0; k < zsize; k++)
        for (long int j = start_y; j <= halfSizeY; j++)
        for (long int i = 0; i < xsize; i++) {
            std::swap(
                direct::elem(*this, i,         j, k, l),
                direct::elem(*this, i, ysize - j, k, l)
            );
        }

        yinit = -lastY();
    }

    /** Reverse matrix values over Z axis, keep in this object.
     *
     * Maybe better with an example:
     *
     * @code
     * slice 0
     * [01 02 03          [11 12 13
     *  04 05 06           14 15 16
     *  07 08 09]          17 18 19]
     *
     *  ----->
     *
     * slice 1
     * [11 12 13          [01 02 03
     *  14 15 16           04 05 06
     *  17 18 19]          07 08 09]
     * @endcode
     *
     */
    void selfReverseZ() {
        long int xsize = xdim;
        long int ysize = ydim;
        long int zsize = zdim;
        long int halfSizeZ = (long int) (zsize - 1) / 2;
        long int nsize = ndim;
        // 0 column should be handled in a different way
        // for even and odd matrices
        long int start_z = !(zsize % 2);

        for (int l = 0; l < nsize; l++)
        for (int k = start_z; k <= halfSizeZ; k++)
        for (int j = 0; j < ysize; j++)
        for (int i = 0; i < xsize; i++) {
            std::swap(
                direct::elem(*this, i, j, k,         l),
                direct::elem(*this, i, j, zsize - k, l)
            );
        }

        zinit = -lastZ();
    }

    /** Write to an ASCII file.
     */
    void write(const FileName &fn) const {
        std::ofstream ofs(fn.c_str(), std::ios::out);
        if (!ofs)
            REPORT_ERROR(static_cast<std::string>("MultidimArray::write: File " + fn + " cannot be opened for output"));

        ofs << *this;
    }

    /** Read from an ASCII file.
     */
    void read(const FileName &fn) {
        std::ifstream ifs(fn.c_str(), std::ios::in);
        if (!ifs)
            REPORT_ERROR(static_cast<std::string>("MultidimArray::read: Cannot read File " + fn));

        ifs >> *this;
    }
    //@}

    /// @name Operators
    /// @{

    // Copy/move assignment
    MultidimArray<T> &operator = (MultidimArray<T> other) {
        copyShape(other);
        data = other.data;
        other.data = nullptr;
        return *this;
    }

    // Unary minus
    MultidimArray<T> operator - () const {
        auto copy (*this);
        for (auto &x : copy) { x = -x; }
        return copy;
    }

    /** Input from input stream.
     *
     * Actual size of the array is used to know how many values must be read.
     *
     * @code
     * v.<3);
     * std::cin >> v;
     * @endcode
     *
     * This function is not ported to Python.
     */
    friend std::istream& operator >> (std::istream& in, MultidimArray<T> &v) {
        for (auto &x : v) { in >> x; }
        return in;
    }

    /** Equality.
     *
     * Do these two objects have the same shape (origin and size)
     * and the same values (to within accuracy)?
     */
    bool equal(const MultidimArray<T> &other, RFLOAT accuracy = Xmipp::epsilon<RFLOAT>()) const {

        if (!sameShape(other) || !data || !other.data) return false;

        for (long int n = 0; n < size(); n++) {
            if (abs((*this)[n] - other[n]) > accuracy)
                return false;
        }

        return true;
    }

    //@}
};

/// @name Functions for all multidimensional arrays
/// @{

/** Conversion from one type to another.
 *
 * If we have an integer array and we need a RFLOAT one,
 * we can use this function.
 * The conversion is done by static_cast-ing each element.
 */
template<typename T, typename U>
MultidimArray<U> typeCast(const MultidimArray<T>& vt) {
    MultidimArray<U> vu;
    if (vt.size() == 0) {
        vu.clear();
        return vu;
    }

    vu.resize(vt);
    T *ptrt;
    U *ptru;
    for (ptrt = vt.begin(), ptru = vu.begin(); ptrt != vt.end(); ++ptrt, ++ptru) {
        *ptru = static_cast<U>(*ptrt);
    }
    return vu;
}

template<typename T>
bool operator == (const MultidimArray<T> &op1, const MultidimArray<T> &op2) {
    return op1.equal(op2);
}

template<typename T>
bool operator != (const MultidimArray<T> &op1, const MultidimArray<T> &op2) {
    return !(op1 == op2);
}

/** Reduce both volumes to a common size.
 *
 * Search the range of logical indices for which both volumes have valid values,
 * and cut both to that size, the corresponding origin is automatically
 * computed.
 *
 * @code
 * MultidimArray<RFLOAT> V1(4, 5, 3);
 * V1.xinit = -2;
 * V1.yinit = -2;
 * V1.zinit = -2;
 *
 * MultidimArray<RFLOAT> V2(4, 2, 3);
 * V2.xinit = 0;
 * V2.yinit = 0;
 * V2.zinit = 0;
 *
 * // V1 and V2 range from (0,0,0)=(z,y,x) to (1,1,0)
 * cutToCommonSize(V1, V2);
 * @endcode
 */
template<typename T>
void cutToCommonSize(MultidimArray<T> &V1, MultidimArray<T> &V2) {
    long int z0 = std::max(V1.firstZ(), V2.firstZ());
    long int zF = std::min(V1.lastZ(),  V2.lastZ());
    long int y0 = std::max(V1.firstY(), V2.firstY());
    long int yF = std::min(V1.lastY(),  V2.lastY());
    long int x0 = std::max(V1.firstX(), V2.firstX());
    long int xF = std::min(V1.lastX(),  V2.lastX());

    V1.window(z0, y0, x0, zF, yF, xF);
    V2.window(z0, y0, x0, zF, yF, xF);
}

/** Output to output stream.
 * This function is not ported to Python.
 */
template <typename T>
std::ostream& operator << (std::ostream& ostrm, const MultidimArray<T> &v) {

    if (v.xdim == 0) {
        ostrm << "Empty array";
    }
    ostrm << '\n';

    T max_val = abs(direct::elem(v, 0, 0, 0));
    for (T *ptr = v.begin(); ptr != v.end(); ++ptr) {
        max_val = std::max(max_val, (T) std::abs(*ptr));
    }

    int prec = bestPrecision(max_val, 10);

    if (v.ydim == 1 && v.zdim == 1) {
        for (long int j = v.firstX(); j <= v.lastX(); j++) {
            ostrm << floatToString((RFLOAT) v.elem(0, j, 0), 10, prec) << std::endl;
        }
    } else {
        for (long int l = 0; l < v.ndim; l++) {
            if (v.ndim > 1)
                ostrm << "Image No. " << l << std::endl;
            for (long int k = v.firstZ(); k <= v.lastZ(); k++) {
                if (v.zdim > 1)
                    ostrm << "Slice No. " << k << std::endl;
                for (long int i = v.firstY(); i <= v.lastY(); i++) {
                    for (long int j = v.firstX(); j <= v.lastX(); j++) {
                        ostrm << floatToString((RFLOAT) v.elem(i, j, k), 10, prec) << ' ';
                    }
                    ostrm << std::endl;
                }
            }
        }
    }

    return ostrm;
}

template <typename T>
static bool inXbounds (long int i, const MultidimArray<T> &arr) {
    return arr.firstX() <= i && i <= arr.lastX();
}

template <typename T>
static bool inYbounds (long int i, const MultidimArray<T> &arr) {
    return arr.firstY() <= i && i <= arr.lastY();
}

template <typename T>
static bool inZbounds (long int i, const MultidimArray<T> &arr) {
    return arr.firstZ() <= i && i <= arr.lastZ();
}

//@}

// Template specialisation for MultidimArray<Complex>
template<>
std::ostream& operator << (std::ostream &ostrm, const MultidimArray<Complex> &v);
//@}

/** Arithmetic
 *
 */

// If either the lhs or rhs is an rvalue reference,
// we can steal its resources instead of copy constructing our result.
// If neither is an rvalue reference,
// the overload that takes the lhs by value will copy construct the resources we need.

// In order for this to be fast in the correct case,
// we leave it up to the caller to check
// that lhs and rhs are of the right size and shape
// and to live with the consequences if they are not.

// Move if you can, copy if you must

template <typename T>
inline MultidimArray<T> operator + (MultidimArray<T> lhs, const MultidimArray<T> &rhs) {
    const size_t n = std::min(lhs.size(), rhs.size());
    for (long int i = 0; i < n; i++) { lhs[i] += rhs[i]; }
    return std::move(lhs);
}

template <typename T>
inline MultidimArray<T> operator + (const MultidimArray<T> &lhs, MultidimArray<T> &&rhs) {
    const size_t n = std::min(lhs.size(), rhs.size());
    for (long int i = 0; i < n; i++) { rhs[i] = lhs[i] + rhs[i]; }  // Do not assume commutativity
    return std::move(rhs);
}

template <typename T>
inline MultidimArray<T> operator - (MultidimArray<T> lhs, const MultidimArray<T> &rhs) {
    const size_t n = std::min(lhs.size(), rhs.size());
    for (long int i = 0; i < n; i++) { lhs[i] -= rhs[i]; }
    return std::move(lhs);
}

template <typename T>
inline MultidimArray<T> operator - (const MultidimArray<T> &lhs, MultidimArray<T> &&rhs) {
    const size_t n = std::min(lhs.size(), rhs.size());
    for (long int i = 0; i < n; i++) { rhs[i] = lhs[i] - rhs[i]; }
    return std::move(rhs);
}

template <typename T>
inline MultidimArray<T> operator * (MultidimArray<T> lhs, const MultidimArray<T> &rhs) {
    const size_t n = std::min(lhs.size(), rhs.size());
    for (long int i = 0; i < n; i++) { lhs[i] *= rhs[i]; }
    return std::move(lhs);
}

template <typename T>
inline MultidimArray<T> operator * (const MultidimArray<T> &lhs, MultidimArray<T> &&rhs) {
    const size_t n = std::min(lhs.size(), rhs.size());
    for (long int i = 0; i < n; i++) { rhs[i] = lhs[i] * rhs[i]; }
    return std::move(rhs);
}

template <typename T>
inline MultidimArray<T> operator / (MultidimArray<T> lhs, const MultidimArray<T> &rhs) {
    const size_t n = std::min(lhs.size(), rhs.size());
    for (long int i = 0; i < n; i++) { lhs[i] /= rhs[i]; }
    return std::move(lhs);
}

template <typename T>
inline MultidimArray<T> operator / (const MultidimArray<T> &lhs, MultidimArray<T> &&rhs) {
    const size_t n = std::min(lhs.size(), rhs.size());
    for (long int i = 0; i < n; i++) { rhs[i] = lhs[i] / rhs[i]; }
    return std::move(rhs);
}

/** @name Scalar operations
 *
 * These operations are between an array and a scalar (of the array's value type).
 * The result will be of the same type as the operands.
 *
 * In this kind of operation each element of the lhs array is reassigned
 * according to the operation in question and the given scalar.
 * The result will have the same shape as the input array.
 */
//@{

template <typename T>
inline MultidimArray<T>& operator += (MultidimArray<T> &lhs, T scalar) {
    for (auto &x : lhs) x += scalar; return lhs;
}

template <typename T>
inline MultidimArray<T>& operator -= (MultidimArray<T> &lhs, T scalar) {
    for (auto &x : lhs) x -= scalar; return lhs;
}

template <typename T>
inline MultidimArray<T>& operator *= (MultidimArray<T> &lhs, T scalar) {
    for (auto &x : lhs) x *= scalar; return lhs;
}

template <typename T>
inline MultidimArray<T>& operator /= (MultidimArray<T> &lhs, T scalar) {
    for (auto &x : lhs) x /= scalar; return lhs;
}

template <typename T>
inline MultidimArray<T> operator + (MultidimArray<T> lhs, T scalar) {
    return std::move(lhs += scalar);
}

template <typename T>
inline MultidimArray<T> operator - (MultidimArray<T> lhs, T scalar) {
    return std::move(lhs -= scalar);
}

template <typename T>
inline MultidimArray<T> operator * (MultidimArray<T> lhs, T scalar) {
    return std::move(lhs *= scalar);
}

template <typename T>
inline MultidimArray<T> operator / (MultidimArray<T> lhs, T scalar) {
    return std::move(lhs /= scalar);
}

//@}

/** @name Scalar-by-array operations
 *
 * These operations are between an array and a scalar (of the array's value type).
 * The result will be of the same type as the array.
 *
 * The operation in question is performed with the scalar on the lhs and each element of the array on the rhs.
 * The result will have the same shape as the input array.
 *
 */

template <typename T>
inline MultidimArray<T> operator + (T scalar, MultidimArray<T> rhs) {
    for (auto &x : rhs) x = scalar + x; return std::move(rhs);
}

template <typename T>
inline MultidimArray<T> operator - (T scalar, MultidimArray<T> rhs) {
    for (auto &x : rhs) x = scalar - x; return std::move(rhs);
}

template <typename T>
inline MultidimArray<T> operator * (T scalar, MultidimArray<T> rhs) {
    for (auto &x : rhs) x = scalar * x; return std::move(rhs);
}

template <typename T>
inline MultidimArray<T> operator / (T scalar, MultidimArray<T> rhs) {
    for (auto &x : rhs) x = scalar / x; return std::move(rhs);
}

#endif
