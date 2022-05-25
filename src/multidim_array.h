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

#include <typeinfo>
#include <sys/mman.h>
#include <fcntl.h>
#include <unistd.h>
#include "src/funcs.h"
#include "src/error.h"
#include "src/args.h"
#include "src/matrix1d.h"
#include "src/matrix2d.h"
#include "src/complex.h"
#include <limits>
using std::pair;

template<typename T>
struct Stats {

    RFLOAT avg, stddev;
    T min, max;

    /** Print statistics
     *
     * No end of line character is written after this print out.
     *
     * @code
     * Stats<RFLOAT> stats = arr.computeStats();
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

// Intel MKL provides an FFTW-like interface, so this is enough.
#include <fftw3.h>
#define RELION_ALIGNED_MALLOC fftw_malloc
#define RELION_ALIGNED_FREE fftw_free

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

/** Return the first X valid logical index
 */
template <typename T>
class MultidimArray;

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
 * This macro is used to generate loops for the array in an easy
 * manner. It defines internal indices 'l', 'k','i' and 'j' which
 * ranges over the n volume using its physical definition.
 *
 * @code
 * FOR_ALL_DIRECT_NZYX_ELEMENTS_IN_MULTIDIMARRAY(v) {
 *     std::cout << direct::elem(v, i, j, k, l) << " ";
 * }
 * @endcode
 */
#define FOR_ALL_DIRECT_NZYX_ELEMENTS_IN_MULTIDIMARRAY(V) \
    for (long int l = 0; l < Nsize(V); l++) \
    for (long int k = 0; k < Zsize(V); k++) \
    for (long int j = 0; j < Ysize(V); j++) \
    for (long int i = 0; i < Xsize(V); i++)

/** For all elements in the array
 *
 * This macro is used to generate loops for the array in an easy
 * manner. It defines internal indices 'l', 'k','i' and 'j' which
 * ranges over the n volume using its logical definition.
 *
 * @code
 * FOR_ALL_NZYX_ELEMENTS_IN_MULTIDIMARRAY(v) {
 *     std::cout << v.elem(i, j, k, l) << " ";
 * }
 * @endcode
 */
#define FOR_ALL_NZYX_ELEMENTS_IN_MULTIDIMARRAY(V) \
    for (long int l = 0; l < Nsize(V); l++) \
    for (long int k = Zinit(V); k <= Zlast(V); k++) \
    for (long int j = Yinit(V); j <= Ylast(V); j++) \
    for (long int i = Xinit(V); i <= Xlast(V); i++)

/** For all direct elements in the array, pointer version
 *
 * This macro is used to generate loops for the array in an easy manner. It
 * defines an internal index 'k' which goes over the slices and 'n' that
 * goes over the pixels in each slice. Each element can be accessed through
 * an external pointer called ptr.
 *
 * @code
 * T *ptr = NULL;
 * long int n;
 * FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY_ptr(v, n, ptr) {
 *     std::cout << *ptr << " ";
 * }
 * @endcode
 */
#define FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY_ptr(v, n, ptr) \
    for ((n) = 0, (ptr) = (v).data; (n) < (v).size(); ++(n), ++(ptr))

/** For all elements in the array.
 *
 * This macro is used to generate loops for the volume in an easy way. It
 * defines internal indices 'k','i' and 'j' which ranges the volume using its
 * mathematical definition (ie, logical access).
 *
 * @code
 * FOR_ALL_ELEMENTS_IN_ARRAY3D(V) {
 *     std::cout << V(k, i, j) << " ";
 * }
 * @endcode
 */
#define FOR_ALL_ELEMENTS_IN_ARRAY3D(V) \
    for (long int k = Zinit(V); k <= Zlast(V); k++) \
    for (long int j = Yinit(V); j <= Ylast(V); j++) \
    for (long int i = Xinit(V); i <= Xlast(V); i++)

/** For all direct elements in the array.
 *
 * This macro is used to generate loops for the volume in an easy way. It
 * defines internal indices 'k','i' and 'j' which ranges the volume using its
 * physical definition.
 *
 * @code
 * FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY3D(V) {
 *     std::cout << direct::elem(m, i, j, k) << " ";
 * }
 * @endcode
 */
#define FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY3D(V) \
    for (long int k = 0; k < Zsize(V); k++) \
    for (long int j = 0; j < Ysize(V); j++) \
    for (long int i = 0; i < Xsize(V); i++)

/** For all elements in the array
 *
 * This macro is used to easily loop through a matrix.
 * It defines logical indices 'i' and 'j' which range over the matrix.
 *
 * @code
 * FOR_ALL_ELEMENTS_IN_ARRAY2D(m) {
 *     std::cout << m(i, j) << " ";
 * }
 * @endcode
 */
#define FOR_ALL_ELEMENTS_IN_ARRAY2D(m) \
    for (long int j = Yinit(m); j <= Ylast(m); j++) \
    for (long int i = Xinit(m); i <= Xlast(m); i++)

/** For all elements in the array, accessed physically
 *
 * This macro is used to easily loop through a matrix using physical indices.
 * It defines physical indices 'i' and 'j' which range over the matrix.
 *
 * @code
 * FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY2D(m) {
 *     std::cout << direct::elem(m, i, j) << " ";
 * }
 * @endcode
 */
#define FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY2D(m) \
    for (long int j = 0; j < Ysize(m); j++) \
    for (long int i = 0; i < Xsize(m); i++)

/** For all elements in the array, accessed physically
 *
 * This macro is used to generate loops for the vector in an easy way using
 * physical indices. It defines internal the index 'i' which ranges the vector
 * using its physical definition.
 *
 * @code
 * FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY1D(v) {
 *     std::cout << direct::elem(v, i) << " ";
 * }
 * @endcode
 */
#define FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY1D(v) for (long int i = 0; i < v.xdim; i++)

/** Direct access
 *
 * Be careful. 
 * These functions give physical/direct access to an array element
 * without taking its logical position into account.
 * These functions go against the array library philosophy,
 * and should therefore be avoided.
 * Arrays usually follow the C convention of starting index == 0.
 */
namespace direct {

    // Direct access into a 1D array (vector)
    template <typename T>
    inline T& elem(const MultidimArray<T> &v, long int i) {
        return v.data[i];
    }

    // Direct access into a 2D array (matrix)
    template <typename T>
    inline T& elem(const MultidimArray<T> &v, long int i, long int j) {
        return v.data[i + j * v.xdim];
    }

    // Direct access into a 3D array
    template <typename T>
    inline T& elem(const MultidimArray<T> &v, long int i, long int j, long int k) {
        return v.data[i + j * v.xdim + k * v.xdim * v.ydim];
    }

    // Direct access into a 4D array
    template <typename T>
    inline T& elem(MultidimArray<T> v, long int i, long int j, long int k, long int l) {
        return v.data[i + j * v.xdim + k * v.xdim * v.ydim + l * v.xdim * v.ydim * v.zdim];
    }

};
//@}

// Forward declarations ====================================================
template<typename T>
class MultidimArray;

struct Origin {

    long int x, y, z;

    Origin(long int x, long int y, long int z): x(x), y(y), z(z) {}

    bool operator == (Origin other) {
        return x == other.x && y == other.y && z == other.z;
    }

    bool operator != (Origin other) {
        return !((*this) == other);
    }

};

/** Template class for Xmipp arrays.
  * This class provides physical and logical access.
*/
template<typename T>
class MultidimArray {

    public:

    /* The array itself.
    The array is always a 3D array (Z,Y,X) (i.e. a 3rd-order tensor).

    For vectors  (1st-order tensors), the shape of the array is (1, 1, m) (so size = m).
    For matrices (2nd-order tensors), the shape of the array is (1, n, m) (so size = m × n).

    The pixel (i,j) (y,x) is at the position data[i * Xdim + j] or data[y * Xdim + x]
    */
    T *data;

    // Destroy data
    bool destructible;

    struct Dimensions {
        // Essentially a 4D vector.
        // For vectors and matrices, the higher order dimensions will be 1:
        // (x, 1, 1) or (x, y, 1).

        long int x, y, z, n;

        Dimensions(long int x, long int y, long int z, long int n):
        x(x), y(y), z(z), n(n) {}

        template <typename T2>
        bool operator == (T2 other) {
            return x == other.x && y == other.y && z == other.z && n == other.n;
        }

        template <typename T2>
        bool operator != (T2 other) {
            return !((*this) == other);
        }

    };

    /** Number of elements in X/Y/Z/N
     * Conventions:
     * - The Z dimension splits our array into slices
     * - The N dimension splits our array into images
     */
    unsigned long int xdim, ydim, zdim, ndim;

    /// TODO: Manage access to xdim, ydim, zdim, ndim.

    // Array size (number of elements)
    inline unsigned long int size() const { return xdim * ydim * zdim * ndim; }

    // X/Y/Zinit
    long int xinit, yinit, zinit;

    T* begin() { return data; }
    T* end() { return &data[size()]; }

    // Logical array access

    // Vector element
    inline T& elem(long int i) const {
        return direct::elem(*this, i - xinit);
    }

    // Matrix element
    inline T& elem(long int i, long int j) const {
        return direct::elem(*this, i - xinit, j - yinit);
    }

    // Volume element
    inline T& elem(long int i, long int j, long int k) const {
        return direct::elem(*this, i - xinit, j - yinit, k - zinit);
    }

    // Multidim element
    inline T& elem(long int i, long int j, long int k, long int l) const {
        return direct::elem(*this, i - xinit, j - yinit, k - zinit, l);
    }

    private:

    // Allocation-related member variables

    bool mmapOn;  // Whether to allocate memory or map to a file
    FileName mapFile;  // Mapped file name
    int mFd;  // Mapped file handler
    long int allocated_size;  // Number of elements in allocated memory

    T* attempt_mmap(FileName &mapFile, int &mFd, off_t offset) {
        mapFile.initRandom(8);
        mapFile = mapFile.addExtension("tmp");

        if ((mFd = open(mapFile.c_str(), O_RDWR | O_CREAT | O_TRUNC, S_IRUSR | S_IWUSR | S_IRGRP | S_IWGRP)) == -1)
            REPORT_ERROR((std::string) "MultidimArray<T>::" + __func__ + ": Error creating map file.");

        if (lseek(mFd, offset, SEEK_SET) == -1 || ::write(mFd, "", 1) == -1) {
            // Use global ::write (conflict with MultidimArray<T>::write)
            close(mFd);
            REPORT_ERROR((std::string) "MultidimArray<T>::" + __func__ + ": Error 'stretching' the map file.");
        }

        T *ptr;
        if ((ptr = (T*) mmap(0, size() * sizeof(T), PROT_READ | PROT_WRITE, MAP_SHARED, mFd, 0)) == (void*) -1)
            REPORT_ERROR((std::string) "MultidimArray<T>::" + __func__ + ": mmap failed.");
        return ptr;

    }

    public:

    /// @name Constructors
    //@{

    /** Empty constructor.
     * Create an empty array with no memory associated.
     */
    MultidimArray() {
        coreInit();
    }

    /** Size constructor with 4D size.
     * The Size constructor creates an array with memory associated,
     * and fills it with zeros.
     */
    MultidimArray(long int Ndim, long int Zdim, long int Ydim, long int Xdim) {
        coreInit();
        resize(Xdim, Ydim, Zdim, Ndim);
    }

    /** Size constructor with 3D size.
     * The Size constructor creates an array with memory associated,
     * and fills it with zeros.
     */
    MultidimArray(long int Zdim, long int Ydim, long int Xdim) {
        coreInit();
        resize(Xdim, Ydim, Zdim);
    }

    /** Size constructor with 2D size.
     * The Size constructor creates an array with memory associated,
     * and fills it with zeros.
     */
    MultidimArray(long int Ydim, long int Xdim) {
        coreInit();
        resize(Xdim, Ydim);
    }

    /** Size constructor with 1D size.
     * The Size constructor creates an array with memory associated,
     * and fills it with zeros.
     */
    MultidimArray(long int Xdim) {
        coreInit();
        resize(Xdim);
    }

    /** Copy constructor
     *
     * The created volume is a perfect copy of the input array but with a
     * different memory assignment.
     *
     * @code
     * MultidimArray<RFLOAT> V2(V1);
     * @endcode
     */
    MultidimArray(const MultidimArray<T> &V, bool parent=false) {
        coreInit();
        if (parent) {
            copyShape(V);
            coreAllocate();
        } else {
            *this = V;
        }
    }

    template <typename T2>
    MultidimArray<T>(const MultidimArray<T2> &other) {
        coreInit();
        *this = other;
    }

    /** Constructor from a Matrix1D.
     * The Size constructor creates an array with memory associated,
     * and fills it with zeros.
     */
    MultidimArray(const Matrix1D<T> &v) {
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

    /** Destructor.
     */
    ~MultidimArray() {
        coreDeallocate();
    }

    /** Clear.
     */
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
        ndim = 1;
        zdim = 1;
        ydim = 1;
        xdim = 0;
        zinit = 0;
        yinit = 0;
        xinit = 0;
        data = NULL;
        allocated_size = 0;
        destructible = true;
        mmapOn = false;
        mFd = 0;
    }

    /** Core allocate with dimensions.
     */
    void coreAllocate(unsigned long int _xdim, unsigned long int _ydim, unsigned long int _zdim, unsigned long int _ndim) {
        if (_xdim == 0 || _ydim == 0 || _zdim == 0 || _ndim == 0) {
            clear();
            return;
        }

        if (data)
            REPORT_ERROR("Do not allocate space for an image if you have not first deallocated it!");
        // Why? coreAllocate() will do this check too.

        ndim = _ndim;
        zdim = _zdim;
        ydim = _ydim;
        xdim = _xdim;

        coreAllocate();
    }

    /** Core allocate without dimensions.
     *
     * The dimensions should be set beforehand.
     */
    void coreAllocate() {

        if (data)
            REPORT_ERROR("Do not allocate space for an image if you have not first deallocated it!");

        _allocate_memory();
    }

    void _allocate_memory() {
        if (mmapOn) {
            data = attempt_mmap(mapFile, mFd, size() * sizeof(T));
        } else {
            data = (T*) RELION_ALIGNED_MALLOC(size() * sizeof(T));
            if (!data) REPORT_ERROR("Allocate: No space left");
        }
        allocated_size = size();
    }

    /** Core allocate without dimensions.
     *
     * The dimensions should be set beforehand.
     */
    void coreAllocateReuse() {

        if (data && size() <= allocated_size) {
            return;  // Memory already allocated
        } else if (size() > allocated_size) {
            coreDeallocate();
        }

        _allocate_memory();
    }

    void setMmap(bool mmap) { mmapOn = mmap; }

    bool getMmap() { return mmapOn; }

    /** Core deallocate.
     * Free all data.
     * Essential to avoid memory leaks.
     */
    void coreDeallocate() {
        if (data && destructible) {
            if (mmapOn) {
                munmap(data, allocated_size * sizeof(T));
                close(mFd);
                remove(mapFile.c_str());
            } else {
                RELION_ALIGNED_FREE(data);
            }
        }
        data = NULL;
        allocated_size = 0;
    }

    // A job for shared_ptr?

    /** Alias a multidimarray.
     *
     * Treat the multidimarray as if it were a volume. The data is not copied
     * into new memory, but a pointer to the multidimarray is copied.
     * You should not make any operation on this volume so the
     * memory locations are changed
     */
    void alias(const MultidimArray<T> &other) {
        coreDeallocate();  // Otherwise there may be a memory leak!
        copyShape(other);
        data = other.data;
        destructible = false;
    }

    /** Move from a multidimarray.
     *
     * Treat the multidimarray as if it were a volume. The data is not copied
     * into new memory, but a pointer to the multidimarray is copied.
     *
     * After the operation, the operand m will become an alias of this array.
     * Same operation as alias, but reverse the relation between the two arrays
     */
    void moveFrom(MultidimArray<T> &other) {
        coreDeallocate(); // Otherwise there may be a memory leak!
        copyShape(other);
        data = other.data;
        destructible = true;
        allocated_size = other.allocated_size;
        other.destructible = false;
        other.allocated_size = 0;
    }

    //@}

    /// @name Size
    //@{

    /** Sets new 4D dimensions.
     *
     *  Note that the dataArray is NOT resized. This should be done separately with coreAllocate()
     *
     */
    void setDimensions(long int Xdim, long int Ydim, long int Zdim, long int Ndim) {
        ndim = Ndim;
        zdim = Zdim;
        ydim = Ydim;
        xdim = Xdim;
    }

    /** NOTE: When setXdim/setYdim/setZdim/setNdim are used, the array is not resized.
     * This should be done separately with coreAllocate()
     */

    void setNdim(long int Ndim) { ndim = Ndim; }

    void setZdim(long int Zdim) { zdim = Zdim; }

    void setYdim(long int Ydim) { ydim = Ydim; }

    void setXdim(long int Xdim) { xdim = Xdim; }

    /** Copy the shape parameters
     *
     */
    void copyShape(const MultidimArray<T> &m) {
        ndim = m.ndim;
        zdim = m.zdim;
        ydim = m.ydim;
        xdim = m.xdim;
        zinit = m.zinit;
        yinit = m.yinit;
        xinit = m.xinit;
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
        if (!destructible)
            REPORT_ERROR("Array marked as indestructible!");
        if (!data || mmapOn || size() <= 0 || allocated_size <= size())
            return;
        T* old_data = data;
        data = (T*) RELION_ALIGNED_MALLOC(sizeof(T) * size());
        memcpy(data, old_data, sizeof(T) * size());
        RELION_ALIGNED_FREE(old_data);
        allocated_size = size();
    }

    /**
     * The functions reshape, moveFrom and shrinkToFit were added
     * on suggestion by Yunxiao Zhang (5 April 2016)
     */

    /** Adjust array to a given shape
     *
     * This function will resize the array to the given size.
     * No data will be copied/moved to the new space.
     * If shape is unchanged, then so is the data.
     * Otherwise, data is almost always destroyed.
     */
    void reshape(long Xdim = 1, long Ydim = 1, long Zdim = 1, long Ndim = 1) {
        if (data && allocated_size == Xdim * Ydim * Zdim * Ndim) {
            setDimensions(Xdim, Ydim, Zdim, Ndim);
            return;
        }

        if (Xdim <= 0 || Ydim <= 0 || Zdim <= 0 || Ndim <= 0) {
            clear();
            return;
        }

        coreDeallocate();
        coreAllocate(Xdim, Ydim, Zdim, Ndim);
    }

    void dynamic_reshape(long dim, int dimensionality) {
        // Dynamic dispatch
               if (dimensionality == 1) {
            reshape(dim);
        } else if (dimensionality == 2) {
            reshape(dim, dim);
        } else if (dimensionality == 3) {
            reshape(dim, dim, dim);
        } else if (dimensionality == 4) {
            reshape(dim, dim, dim, dim);
        }
    }

    /** Adjust shape to match the target array
     *
     * No guarantee about the data stored
     */
    template<typename T1>
    void reshape(const MultidimArray<T1> &v) {
        if (
            ndim != v.ndim || xdim != v.xdim ||
            ydim != v.ydim || zdim != v.zdim ||
            !data
        ) { reshape(v.xdim, v.ydim, v.zdim, v.ndim); }

        xinit = v.xinit;
        yinit = v.yinit;
        zinit = v.zinit;
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
    void resizeNoCp(unsigned long int Xdim = 1, unsigned long int Ydim = 1, unsigned long int Zdim = 1, unsigned long int Ndim = 1) {

        size_t NZYXdim = Ndim * Zdim * Ydim * Xdim;
        if (NZYXdim == allocated_size && data)
            return;

        if (NZYXdim == 0) {
            clear();
            return;
        }

        // data can be NULL even when xdim etc are not zero
        // (This can happen for reading of images...)
        // In that case, initialize data to zeros.
        if (size() > 0 && !data) {
            coreAllocate();
            return;
        }

        // Ask for memory
        int new_mFd = 0;
        FileName new_mapFile;

        T *new_data;

        try {
            if (mmapOn) {
                new_data = attempt_mmap(new_mapFile, new_mFd, NZYXdim * sizeof(T) - 1);
            } else {
                new_data = (T*) RELION_ALIGNED_MALLOC(NZYXdim * sizeof(T));
            }
        } catch (std::bad_alloc &) {
            REPORT_ERROR("Allocate: No space left");
        }

        // deallocate old vector
        coreDeallocate();

        // assign *this vector to the newly created
        setDimensions(Xdim, Ydim, Zdim, Ndim);
        data = new_data;
        mFd = new_mFd;
        mapFile = new_mapFile;
        allocated_size = size();
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
    void resize(unsigned long int Xdim = 1, unsigned long int Ydim = 1, unsigned long int Zdim = 1, unsigned long int Ndim = 1) {
        size_t NZYXdim = Ndim * Zdim * Ydim * Xdim;
        if (allocated_size == NZYXdim && data) {
            setDimensions(Xdim, Ydim, Zdim, Ndim);
            allocated_size = size();
            return;
        }

        if (NZYXdim == 0) {
            clear();
            return;
        }

        // data can be NULL even when xdim etc are not zero
        // (This can happen for reading of images...)
        // In that case, initialize data to zeros.
        if (size() > 0 && !data) {
            coreAllocate();
            return;
        }

        // Ask for memory
        T *new_data;
        int new_mFd = 0;
        FileName new_mapFile;

        try {
            if (mmapOn) {
                new_data = attempt_mmap(new_mapFile, new_mFd, NZYXdim * sizeof(T) - 1);
            } else {
                new_data = (T*) RELION_ALIGNED_MALLOC(NZYXdim * sizeof(T));
            }
        } catch (std::bad_alloc &) {
            REPORT_ERROR("Allocate: No space left");
        }

        // Copy needed elements, fill with 0 if necessary
        for (long int l = 0; l < Ndim; l++)
        for (long int k = 0; k < Zdim; k++)
        for (long int j = 0; j < Ydim; j++)
        for (long int i = 0; i < Xdim; i++) {
            // 0 if out of bounds
            new_data[i + j * Xdim + k * Xdim * Ydim + l * Xdim * Ydim * Zdim] = (
                k >= zdim || i >= ydim || j >= xdim
            ) ? (T) 0.0 : direct::elem(*this, i, j, k);
        }

        // deallocate old vector
        coreDeallocate();

        // assign *this vector to the newly created
        data = new_data;
        setDimensions(Xdim, Ydim, Zdim, Ndim);
        mFd = new_mFd;
        mapFile = new_mapFile;
        allocated_size = size();
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
    void resize(const MultidimArray<T1> &other) {
        if (
            ndim != other.ndim || xdim != other.xdim ||
            ydim != other.ydim || zdim != other.zdim ||
            !data
        ) { resize(other.xdim, other.ydim, other.zdim, other.ndim); }

        xinit = other.xinit;
        yinit = other.yinit;
        zinit = other.zinit;
    }

    /** Return a struct of the array's X/Y/Z/N dimensions.
     *
     * Could also be considered the "size" of the array.
     *
     * @code
     * dimensions = V.getDimensions();
     * @endcode
     */
    Dimensions getDimensions() const {
        return Dimensions(xdim, ydim, zdim, ndim);
    }

    Origin getOrigin() const {
        return Origin(xinit, yinit, zinit);
    }

    /** The dimension of an array.
     *
     * The number of indices needed to select an element.
     *
     * @code
     * int dim = V.getDim();
     * @endcode
     */
    inline int getDim() const {
        if (size() < 1) return 0; // 0-dimensional (dimensionless) array (scalar)
        if (zdim > 1) return 3;      // 3-dimensional array
        if (ydim > 1) return 2;      // 2-dimensional array (matrix)
        return 1;                    // 1-dimensional array (vector)
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

    /** Generic window routine (dim independent)
     *
     * This function will call to 3D,2D or 1D specific window routines
     */
    void window(
        long int n0, long int z0, long int y0, long int x0,
        long int nF, long int zF, long int yF, long int xF,
        T init_value = 0, long n = 0
    ) {
        if (ndim > 1)
            REPORT_ERROR("stack windowing not implemented");
        if (zdim > 1) {
            //call 3Dwindow
            window(z0, y0, x0, zF, yF, xF, init_value, n);
        } else if (ydim > 1) {
            //call 2Dwindow
            window(y0, x0, yF, xF, init_value, n);
        } else if (xdim > 1) {
            //call 1Dwindow
            window(x0, xF, init_value, n);
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
    void window(
        MultidimArray<T> &result,
        long int z0, long int y0, long int x0,
        long int zF, long int yF, long int xF,
        T init_value = 0, long n = 0
    ) const {
        result.resize(xF - x0 + 1, yF - y0 + 1, zF - z0 + 1);
        result.xinit = x0;
        result.yinit = y0;
        result.zinit = z0;

        for (long int k = z0; k <= zF; k++)
        for (long int j = y0; j <= yF; j++)
        for (long int i = x0; i <= xF; i++) {
            result.elem(i, j, k) = inside(i, j, k) ?
                elem(i, j, k, n) : init_value;
        }
    }

    // As above but acts on itself
    void window(
        long int z0, long int y0, long int x0,
        long int zF, long int yF, long int xF,
        T init_value = 0, long n = 0
    ) {
        MultidimArray<T> result;
        window(result, z0, y0, x0, zF, yF, xF, init_value, n);
        moveFrom(result);
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
    void window(
        MultidimArray<T> &result,
        long int y0, long int x0,
        long int yF, long int xF,
        T init_value = 0, long n = 0
    ) const {
        result.resize(xF - x0 + 1, yF - y0 + 1);
        result.xinit = x0;
        result.yinit = y0;

        FOR_ALL_ELEMENTS_IN_ARRAY2D(result) {
            result.elem(i, j) = inside(i, j) ?
                elem(i, j, 0, n) : init_value;
        }
    }

    // As above but acts on itself
    void window(
        long int y0, long int x0,
        long int yF, long int xF,
        T init_value = 0, long n = 0
    ) {
        MultidimArray<T> result;
        window(result, y0, x0, yF, xF, init_value, n);
        *this = result;

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
     * v1.window(-1, 2); // v1=[-1 0 1 2]; v1.firstX() == -1
     *
     * v1.window(-3, 1); // v1=[0 -2 -1 0 1]; v1.firstX() == -3
     * @endcode
     */
    void window(
        MultidimArray<T> &result,
        long int x0,
        long int xF,
        T init_value = 0, long n = 0
    ) const {
        result.resize(xF - x0 + 1);
        result.xinit = x0;

        for (long int i = x0; i <= xF; i++) {
            result.elem(i) = inside(i) ?
                elem(i, 0, 0, n) : init_value;
            }
    }

    // As above but acts on itself
    void window(long int x0, long int xF, T init_value = 0, long n = 0) {
        MultidimArray<T> result;
        window(result, x0, xF, init_value, n);
        *this = result;
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
    template <typename T1>
    inline bool sameShape(const MultidimArray<T1> &other) const {
        return getDimensions() == other.getDimensions() && // Same dimensions
               getOrigin()     == other.getOrigin();       // Same origin
    }

    /** inside for 3D matrices */
    bool inside(long int i, long int j, long int k) const {
        // Why j <-> X and i <-> Y?
        return inXbounds(j, *this) && inYbounds(i, *this) && inZbounds(k, *this);
    }

    /** inside for 2D matrices */
    bool inside(long int i, long int j) const {
        // Why j <-> X and i <-> Y?
        return inXbounds(j, *this) && inYbounds(i, *this);
    }

    /** inside for 1D matrices */
    bool inside(long int i) const {
        return inXbounds(i, *this);
    }

    /** outside for 3D matrices */
    bool outside(long int i, long int j, long int k) const {
        return !inside(i, j, k);
    }

    /** outside for 2D matrices */
    bool outside(long int i, long int j) const {
        return !inside(i, j);
    }

    /** outside for 1D matrices */
    bool outside(long int i) const {
        return !inside(i);
    }

    /** Is this logical index vector outside the bounds of this array? */
    bool outside(const Matrix1D<RFLOAT> &r) const {
        int rsize = r.size();
        if (rsize < 1) REPORT_ERROR(std::string(__func__) + ": index vector has too few components");

        switch (rsize) {
            case 1:
            return outside(XX(r));
            case 2:
            return outside(YY(r), XX(r));
            case 3:
            return outside(YY(r), XX(r), ZZ(r));
            default:
            REPORT_ERROR(std::string(__func__) + ": index vector has too many components");
        }
    }

    /** Return Y dimension. */
    inline long int rowNumber() const { return ydim; }

    /** Return X dimension. */
    inline long int colNumber() const { return xdim; }

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
    void setXmippOrigin() {
        zinit = Xmipp::init(zdim);
        yinit = Xmipp::init(ydim);
        xinit = Xmipp::init(xdim);
    }

    /** The first logical Z index. */
    inline long int firstZ() const {
        return Xmipp::init(zdim);
    }

    /** The last logical Z index. */
    inline long int lastZ() const {
        return Xmipp::last(zdim);
    }

    /** The first logical Y index. */
    inline long int firstY() const {
        return Xmipp::init(ydim);
    }

    /** The last logical Y index. */
    inline long int lastY() const {
        return Xmipp::last(ydim);
    }

    /** The first logical X index. */
    inline long int firstX() const {
        return Xmipp::init(xdim);
    }

    /** The last logical X index. */
    inline long int lastX() const {
        return Xmipp::last(xdim);
    }

    /** IsCorner (in 2D or 3D matrix)
     *
     * Is this logical index a corner of this array?
     */
    bool isCorner(const Matrix1D<RFLOAT> &v) const {

        if (v.size() < 2)
            REPORT_ERROR(std::string(__func__) + ": index vector has too few components");

        switch (xdim) {

            case 2:
            return (XX(v) == firstX() || XX(v) == lastX()) &&
                   (YY(v) == firstY() || YY(v) == lastY());

            case 3:
            return (XX(v) == firstX() || XX(v) == lastX()) &&
                   (YY(v) == firstY() || YY(v) == lastY()) &&
                   (ZZ(v) == firstZ() || ZZ(v) == lastZ());

            default:
            REPORT_ERROR(std::string(__func__) + ": index vector has too many components");

        }

    }
    //@}

    ///@name Access to the pixel values
    //@{

    T& operator [] (long int n) const { return data[n]; }

    T* begin() const { return data; }
    T* end() const { return &data[size()]; }

    /** Volume element access by RFLOAT vector.
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
    T& operator()(const Matrix1D<RFLOAT> &v) const {
        switch (v.size()) {
            case 1:
            return elem(round(XX(v)));
            case 2:
            return elem(round(XX(v)), round(YY(v)));
            case 3:
            return elem(round(XX(v)), round(YY(v)), round(ZZ(v)));
            default:
            REPORT_ERROR("Matrix dimensions must be 1, 2, or 3");
        }
    }

    /** Volume element access by integer vector. */
    T& operator()(const Matrix1D<long int> &v) const {
        switch (v.size()) {
            case 1:
            return elem(XX(v));
            case 2:
            return elem(XX(v), YY(v));
            case 3:
            return elem(XX(v), YY(v), ZZ(v));
            default:
            REPORT_ERROR("Matrix dimensions must be 1, 2, or 3");
        }
    }

    /** 4D element access by index.
    *
    * Returns the value of a matrix logical position. In our example we could
    * access from v(0, 0,-2,-1) to v(0, 1,2,1). The elements can be used either by
    * value or by reference. An exception is thrown if the index is outside
    * the logical range. Be careful that the argument order is (Z,Y,X).
    *
    * @code
    * V(0, 0, -2, 1) = 1;
    * val = V(0, 0, -2, 1);
    * @endcode
    */
    inline T& operator()(long n, long int k, long int i, long int j) const {
        return elem(i, j, k, n);
    }

    /** 3D element access by index.
    *
    * Returns the value of a matrix logical position. In our example we could
    * access from v(0,-2,-1) to v(1,2,1). The elements can be used either by
    * value or by reference. An exception is thrown if the index is outside
    * the logical range. Be careful that the argument order is (Z,Y,X).
    *
    * @code
    * V(0, -2, 1) = 1;
    * val = V(0, -2, 1);
    * @endcode
    */
    inline T& operator()(long int k, long int i, long int j) const {
        return elem(i, j, k);
    }

    /** 3D element access by index (getVoxel).
    *
    * Same function as operator() but with a name. Needed by swig.
    *
    */
    inline T getVoxel(long int k, long int i, long int j) const {
        return elem(i, j, k);
    }

    /** 3D element access by index (setVoxel).
    *
    * Same function as operator() but with a name. Needed by swig.
    *
    */
    inline void setVoxel(long int k, long int i, long int j, T newval) {
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
    inline T& operator()(long int i, long int j) const {
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
    inline T& operator()(long int i) const {
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
        FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY3D(M) {
            direct::elem(M, i, j) = direct::elem(*this, i, j, k, n);
        }

        M.xinit = firstX();
        M.yinit = firstY();
        M.zinit = firstZ();
    }

    /** Set a single 1,2 or 3D image in a multi-image array
     *
     * Set a single-image array in a multi-image one.
     * @code
     * V.setImage(0, m);
     * @endcode
     */
    void setImage(long n, MultidimArray<T> &M) const {
        if (xdim == 0)
            return;

        if (n < 0 || n > ndim)
            REPORT_ERROR("setImage: MultidimArray subscript (n) out of range");

        if (M.zdim != zdim || M.ydim != ydim || M.xdim != xdim)
            REPORT_ERROR("setImage: MultidimArray dimensions different from the input image ones");

        FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY3D(M)
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
    void getSlice(long int k, MultidimArray<T> &M, char axis = 'Z', long n = 0) const {
        if (xdim == 0) {
            M.clear();
            return;
        }

        switch (axis) {

            case 'Z':
            if (!inZbounds(k, (*this)))
                REPORT_ERROR(std::string(__func__) + ": Multidim subscript (k) out of range");

            k -= firstZ();
            M.resize(xdim, ydim);
            FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY2D(M) {
                direct::elem(M, i, j) = direct::elem(*this, j, i, k, n);
            }
            M.xinit = firstX();
            M.yinit = firstY();
            break;

            case 'Y':
            if (!inYbounds(k, (*this)))
                REPORT_ERROR(std::string(__func__) + ": Multidim subscript (i) out of range");

            k -= firstY();
            M.resize(xdim, zdim);
            FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY2D(M) {
                direct::elem(M, i, j) = direct::elem(*this, k, j, i, n);
            }
            M.xinit = firstX();
            M.yinit = firstZ();
            break;

            case 'X':
            if (!inXbounds(k, (*this)))
                REPORT_ERROR(std::string(__func__) + ": Multidim subscript (j) out of range");

            k -= firstX();
            M.resize(ydim, zdim);
            FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY2D(M) {
                direct::elem(M, i, j) = direct::elem(*this, j, k, i, n);
            }
            M.xinit = firstY();
            M.yinit = firstZ();
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
    void setSlice(long int k, const MultidimArray<T> &v, long n = 0) {
        if (xdim == 0)
            return;

        if (k < firstZ() || k > lastZ())
            REPORT_ERROR(std::string(__func__) + ": MultidimArray subscript (k) out of range");

        if (v.ydim != ydim || v.xdim != xdim)
            REPORT_ERROR(std::string(__func__) + ": MultidimArray dimensions different from the matrix ones");

        k -= firstZ();

        FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY2D(v)
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
    void getCol(long int j, MultidimArray<T> &v) const {
            if (xdim == 0 || ydim == 0) {
                v.clear();
                return;
        }

        if (j < 0 || j >= xdim)
            REPORT_ERROR("getCol: Matrix subscript (j) greater than matrix dimension");

        v.resize(ydim);
        for (long int i = 0; i < ydim; i++)
            v(i) = (*this)(i, j);
    }

    /** Set Column
     *
     * Set a column vector corresponding to the choosen column
     * inside matrix.
     *
     * @code
     * m.setCol(0, (m.row(1)).transpose()); // Copies row 1 in column 0
     * @endcode
     */
    void setCol(long int j, const MultidimArray<T> &v) {
        if (xdim == 0 || ydim == 0)
            REPORT_ERROR("setCol: Target matrix is empty");

        if (j < 0 || j>= xdim)
            REPORT_ERROR("setCol: Matrix subscript (j) out of range");

        if (v.xdim != ydim)
            REPORT_ERROR("setCol: Vector dimension different from matrix one");

        for (long int i = 0; i < ydim; i++)
            (*this)(i, j) = v(i);
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
    void getRow(long int i, MultidimArray<T> &v) const {
        if (xdim == 0 || ydim == 0) {
            v.clear();
            return;
        }

        if (i < 0 || i >= ydim)
            REPORT_ERROR("getRow: Matrix subscript (i) greater than matrix dimension");

        v.resize(xdim);
        for (long int j = 0; j < xdim; j++)
            v(j) = (*this)(i, j);
    }

    /** Set Row
     *
     * Set a row vector corresponding to the choosen row in the 2D Matrix
     *
     * @code
     * m.setRow(-2, m.row(1)); // Copies row 1 in row -2
     * @endcode
     */
    void setRow(long int i, const MultidimArray<T> &v) {
        if (xdim == 0 || ydim == 0)
            REPORT_ERROR("setRow: Target matrix is empty");

        if (i < 0 || i >= ydim)
            REPORT_ERROR("setRow: Matrix subscript (i) out of range");

        if (v.xdim != xdim)
            REPORT_ERROR("setRow: Vector dimension different from matrix one");

        for (long int j = 0; j < xdim; j++)
            (*this)(i, j) = v(j);
    }

    /** 3D Logical to physical index translation.
     *
     * Return the physical position of a logical one.
     *
     * @code
     * m.toPhysical(k_log, i_log, j_log, k_phys, i_phys, j_phys);
     * @endcode
     */
    void toPhysical(
        long int  k_log,  long int  i_log,  long int  j_log,
        long int &k_phys, long int &i_phys, long int &j_phys
    ) const {
        k_phys = k_log - firstZ();
        i_phys = i_log - firstY();
        j_phys = j_log - firstX();
    }

    /** 3D Physical to logical index translation.
     *
     * Return the logical position of a physical one.
     *
     * @code
     * m.toLogical(i_phys, j_phys, i_log, j_log);
     * @endcode
     */
    void toLogical(
        long int  k_phys, long int  i_phys, long int  j_phys,
        long int &k_log,  long int &i_log,  long int &j_log
    ) const {
        k_log = k_phys + firstZ();
        i_log = i_phys + firstY();
        j_log = j_phys + firstX();
    }

    /** 2D Logical to physical index translation
     *
     * Return the physical position of a logical one.
     *
     * @code
     * m.toPhysical(i_log, j_log, i_phys, j_phys);
     * @endcode
     */
    void toPhysical(
        long int  i_log,  long int  j_log,
        long int &i_phys, long int &j_phys
    ) const {
        i_phys = i_log - firstY();
        j_phys = j_log - firstX();
    }

    /** 2D Physical to logical index translation
     *
     * Return the logical position of a physical one.
     *
     * @code
     * m.toLogical(i_phys, j_phys, i_log, j_log);
     * @endcode
     */
    void toLogical(
        long int  i_phys, long int j_phys,
        long int &i_log,  long int &j_log
    ) const {
        i_log = i_phys + firstY();
        j_log = j_phys + firstX();
    }

    /** 1D Logical to physical index translation
     *
     * Return the physical position of a logical one.
     *
     * @code
     * v.toPhysical(i_log, i_phys);
     * @endcode
     */
    void toPhysical(long int i_log, long int &i_phys) const {
        i_phys = i_log - firstX();
    }

    /** 1D Physical to logical index translation.
     *
     * Return the logical position of a physical one.
     *
     * @code
     * v.toLogical(i_phys, i_log);
     * @endcode
     */
    void toLogical(long int i_phys, long int &i_log) const {
        i_log = i_phys + firstX();
    }

    /** Interpolates the value of the nth 3D matrix M at the point (x,y,z).
     *
     * (x,y,z) are in logical coordinates.
     */
    T interpolatedElement3D(RFLOAT x, RFLOAT y, RFLOAT z, T outside_value = (T) 0, long int n = 0) {
        long int x0 = floor(x);
        RFLOAT fx = x - x0;
        long int x1 = x0 + 1;

        long int y0 = floor(y);
        RFLOAT fy = y - y0;
        long int y1 = y0 + 1;

        long int z0 = floor(z);
        RFLOAT fz = z - z0;
        long int z1 = z0 + 1;

        T d000 = outside(y0, x0, z0) ? outside_value : elem(x0, y0, z0, n);
        T d001 = outside(y0, x1, z0) ? outside_value : elem(x1, y0, z0, n);
        T d010 = outside(y1, x0, z0) ? outside_value : elem(x0, y1, z0, n);
        T d011 = outside(y1, x1, z0) ? outside_value : elem(x1, y1, z0, n);
        T d100 = outside(y0, x0, z1) ? outside_value : elem(x0, y0, z1, n);
        T d101 = outside(y0, x1, z1) ? outside_value : elem(x1, y0, z1, n);
        T d110 = outside(y1, x0, z1) ? outside_value : elem(x0, y1, z1, n);
        T d111 = outside(y1, x1, z1) ? outside_value : elem(x1, y1, z1, n);

        RFLOAT dx00 = LIN_INTERP(fx, (RFLOAT) d000, (RFLOAT) d001);
        RFLOAT dx01 = LIN_INTERP(fx, (RFLOAT) d100, (RFLOAT) d101);
        RFLOAT dx10 = LIN_INTERP(fx, (RFLOAT) d010, (RFLOAT) d011);
        RFLOAT dx11 = LIN_INTERP(fx, (RFLOAT) d110, (RFLOAT) d111);
        RFLOAT dxy0 = LIN_INTERP(fy, (RFLOAT) dx00, (RFLOAT) dx10);
        RFLOAT dxy1 = LIN_INTERP(fy, (RFLOAT) dx01, (RFLOAT) dx11);

        return (T) LIN_INTERP(fz, dxy0, dxy1);
    }

    /** Interpolates the value of the nth 2D matrix M at the point (x,y)
     *
     * Bilinear interpolation. (x,y) are in logical coordinates.
     */
    inline T interpolatedElement2D(RFLOAT x, RFLOAT y, T outside_value = (T) 0, long int n = 0) const {
        long int x0 = floor(x);
        RFLOAT fx = x - x0;
        long int x1 = x0 + 1;
        long int y0 = floor(y);
        RFLOAT fy = y - y0;
        long int y1 = y0 + 1;

        T d00 = outside(y0, x0) ? outside_value : elem(x0, y0, 0, n);
        T d10 = outside(y1, x0) ? outside_value : elem(x0, y1, 0, n);
        T d11 = outside(y1, x1) ? outside_value : elem(x1, y1, 0, n);
        T d01 = outside(y0, x1) ? outside_value : elem(x1, y0, 0, n);

        RFLOAT d0 = (T) LIN_INTERP(fx, (RFLOAT) d00, (RFLOAT) d01);
        RFLOAT d1 = (T) LIN_INTERP(fx, (RFLOAT) d10, (RFLOAT) d11);
        return (T) LIN_INTERP(fy, d0, d1);
    }
    //@}

    /// @name Statistics functions
    //@{

    T max() const {

        if (size() <= 0) return static_cast<T>(0);

        T maxval = data[0];

        T *ptr;
        long int n;
        FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY_ptr(*this, n, ptr) {
            if (*ptr > maxval) { maxval = *ptr; }
        }

        return maxval;
    }

    T min() const {

        if (size() <= 0) return static_cast<T>(0);

        T minval = data[0];

        T *ptr;
        long int n;
        FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY_ptr(*this, n, ptr) {
            if (*ptr < minval) {minval = *ptr;}
        }

        return minval;
    }

    /** 4D Indices for the minimum element.
     *
     * Return the index of the minimum element of an array.
     * array(l,k,i,j). Returns -1 if the array is empty
     */
    T minIndex(long int &lmin, long int &kmin, long int &imin, long int &jmin) const {
        if (xdim == 0) {
            lmin = kmin = imin = jmin = -1;
            return 0;
        }

        kmin = firstZ();
        imin = firstY();
        jmin = firstX();
        lmin = 0;
        T minval = elem(imin, jmin, kmin, lmin);


        FOR_ALL_NZYX_ELEMENTS_IN_MULTIDIMARRAY(*this) {
            if (elem(i, j, k, l) > minval) {
                minval = elem(i, j, k, l);
                lmin = l;
                kmin = k;
                imin = i;
                jmin = j;
            }
        }

        return minval;
    }

    /** 3D Indices for the minimum element.
     *
     * Call to the 4D function
     */
    T minIndex(long int &kmin, long int &imin, long int &jmin) const {
        long int zero = 0;
        return minIndex(zero, kmin, imin, jmin);
    }

    /** 2D Indices for the minimum element.
     *
     * Call to the 4D function
     */
    T minIndex(long int &imin, long int &jmin) const {
        long int zero = 0;
        return minIndex(zero, zero, imin, jmin);
    }

    /** 1D Indices for the minimum element.
     *
     * Call to the 4D function
     */
    T minIndex(long int &jmin) const {
        long int zero = 0;
        return minIndex(zero, zero, zero, jmin);
    }

    /** 4D Indices for the maximum element.
     *
     * Return the index of the maximum element of an array.
     * array(l,k,i,j). Returns -1 if the array is empty
     */
    T maxIndex(long int &lmax, long int &kmax, long int &imax, long int &jmax) const {
        if (xdim == 0) {
            lmax = kmax = imax = jmax = -1;
            return 0;
        }

        kmax = firstZ();
        imax = firstY();
        jmax = firstX();
        lmax = 0;
        T maxval = elem(imax, jmax, kmax, lmax);

        FOR_ALL_NZYX_ELEMENTS_IN_MULTIDIMARRAY(*this) {
            if (elem(i, j, k, l) > maxval) {
                maxval = elem(i, j, k, l);
                lmax = l;
                kmax = k;
                imax = i;
                jmax = j;
            }
        }

        return maxval;
    }

    /** 3D Indices for the maximum element.
     *
     * Call to the 4D function
     */
    T maxIndex(long int &kmax, long int &imax, long int &jmax) const {
        long int dum;
        return maxIndex(dum, kmax, imax, jmax);
    }

    /** 2D Indices for the maximum element.
     *
     * Call to the 4D function
     */
    T maxIndex(long int &imax, long int &jmax) const {
        long int dum;
        return maxIndex(dum, dum, imax, jmax);
    }

    /** 1D Indices for the maximum element.
     *
     * Call to the 4D function
     */
    T maxIndex(long int &jmax) const {
        long int dum;
        return maxIndex(dum, dum, dum, jmax);
    }

    struct MinMax {
        RFLOAT min, max;
        MinMax(RFLOAT min, RFLOAT max): min(min), max(max) {}
    };

    /** Minimum and maximum of the values in the array.
     *
     * As RFLOATs.
     */
    MinMax minmax() const {

        RFLOAT min, max;

        if (size() <= 0) return MinMax(min, max); // Uninitialised RFLOATs

        T *ptr;
        long int n;
        FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY_ptr(*this, n, ptr) {
            T val = *ptr;
            if (n == 0) { min = max = static_cast<RFLOAT>(val); }
            else if (val < min) { min = static_cast<RFLOAT>(val); }
            else if (val > max) { max = static_cast<RFLOAT>(val); }
        }
        return MinMax(min, max);
    }

    /** Average of the values in the array.
     *
     * Regardless of the type of the array, return an RFLOAT.
     */
    RFLOAT average() const {
        // Arithmetic mean
        if (size() <= 0) return 0;

        RFLOAT sum = 0;

        T *ptr = NULL;
        long int n;
        // Fold/reduce
        FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY_ptr(*this, n, ptr) {
            sum += static_cast<RFLOAT>(*ptr);
        }
        return sum / size();
    }

    /** Standard deviation of the values in the array.
     *
     * The returned value is always RFLOAT, regardless of the type of the array.
     *
     * stddev(N) = sqrt( sum for (int i = 0; i < N; i++) {(x[i] - mean(N)) ** 2} * 1 / N)
     */
    RFLOAT computeStddev() const {
        if (size() <= 1)
            return 0;

        RFLOAT avg = 0, stddev = 0;
        long int N = size();

        T *ptr = NULL;
        long int n;

        #ifdef RELION_SINGLE_PRECISION
            // Two passes through the data, as single-precision is not enough for a single pass
            // Also: averages of large arrays will give trouble: compute median first.
            RFLOAT median = 0.0;
            if (N > 1e6)
                median = median();

            RFLOAT sumofdeviations = 0;
            FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY_ptr(*this, n, ptr) {
                RFLOAT val = static_cast<RFLOAT>(*ptr);
                sumofdeviations += val - median;
            }
            avg = median + sumofdeviations / N;

            RFLOAT sumofsquareddeviations = 0;
            FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY_ptr(*this, n, ptr) {
                RFLOAT x = static_cast<RFLOAT>(*ptr);
                RFLOAT dev = x - avg;
                sumofsquareddeviations += dev * dev;
            }

            RFLOAT var = 0;

            if (N > 1) {
                var = sumofsquareddeviations / N - 1;
                // Foreseeing numerical instabilities
                stddev = sqrt(static_cast<RFLOAT>(abs(var)));
            } else {
                stddev = 0;
            }
        #else
            RFLOAT total = 0;
            RFLOAT sumofsquares = 0;

            FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY_ptr(*this, n, ptr) {
                RFLOAT x = static_cast<RFLOAT>(*ptr);
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

    /** Compute statistics.
     *
     * Return the average, standard deviation, minimum and maximum.
     * A single pass through the entire array makes this faster
     * than separately computing average, standard deviation, minimum, and maximum.
     */
    Stats<T> computeStats() const {

        if (size() <= 0)
            throw "Statistics cannot be computed for a dimensionless array!";

        double sumx = 0;
        double sumxx = 0;

        Stats<T> stats;
        stats.min =  std::numeric_limits<double>::max();
        stats.max = -std::numeric_limits<double>::max();

        T *ptr;
        long int n;
        // Make one pass through the array.
        FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY_ptr(*this, n, ptr) {
            T val = *ptr;
            double x = static_cast<double>(val);
            sumx  += x;
            sumxx += x * x;

                 if (val > stats.max) { stats.max = val; }
            else if (val < stats.min) { stats.min = val; }
        }

        long int N = size();

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
     * med = v1.median();
     * @endcode
     */
    RFLOAT median() const {

        if (xdim == 0)
            return 0;

        if (xdim == 1)
            return (*this)[0];

        // Copy *this
        long int N = size();
        MultidimArray<RFLOAT> temp(N);
        for (long int n = 0; n < size(); n++) {
            temp[n] = (*this)[n];
        }

        // Sort indices
        temp.sort();

        // Get median
        if (N % 2 == 0) {
            return (temp[N / 2 - 1] + temp[N / 2]) * 0.5;
        } else {
            return temp[N / 2];
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
    void rangeAdjust(T minF, T maxF) {

        if (size() <= 0) return;

        MinMax range = minmax();

        // If range.min == range.max, the vector is a constant one,
        // so the only possible transformation is to a fixed minF
        RFLOAT slope = range.min == range.max ? 0 :
            static_cast<RFLOAT>(maxF - minF) /
            static_cast<RFLOAT>(range.max - range.min);

        T *ptr;
        long int n;
        FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY_ptr(*this, n, ptr) {
            // a + b * x
            *ptr = minF + static_cast<T>(slope * static_cast<RFLOAT>(*ptr - range.min));
        }
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
    void rangeAdjust(T minF, T maxF, MultidimArray<int> &mask) {

        if (size() <= 0) return;

        MinMax range = maskminmax(mask);

        // If range.min == range.max, the vector is a constant one,
        // so the only possible transformation is to a fixed minF
        RFLOAT slope = range.min == range.max ? 0 :
            static_cast<RFLOAT>(maxF - minF) /
            static_cast<RFLOAT>(range.max - range.min);

        T *ptr;
        long int n;
        FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY_ptr(*this, n, ptr) {
            // a + b * x
            *ptr = minF + static_cast<T>(slope * static_cast<RFLOAT>(*ptr - range.min));
        }
    }

    // For use in rangeAdjust
    MinMax maskminmax(MultidimArray<int> &mask) {
        RFLOAT min, max;
        int *maskptr = mask.data;

        T *ptr;
        long int n;
        FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY_ptr(*this, n, ptr) {
            if (*maskptr) {
                T val = *ptr;
                if (n == 0) {
                    min = max = (RFLOAT) val;
                } else {
                    min = std::min(min, (RFLOAT) val);
                    max = std::max(max, (RFLOAT) val);
                }
            }
            maskptr++;
        }
        return MinMax(min, max);
    }

    /** Adjust the range of the array to the range of another array in
        a least squares sense.
    *
    * A linear operation is performed on the values of the array so
    * after it, the values of the self array are as similar as possible
    * (L2 sense) to the values of the array shown as sample
    */

    // As written this will only work for T=RFLOAT
    // nevertheless since this is used is better
    // to use T than RFLOAT or will create problem for int multidim arrays
    void rangeAdjust(
        const MultidimArray<T> &target, const MultidimArray<int> *mask = NULL
    ) {

        if (size() <= 0) return;

        T *targetptr = target.data;
        int *maskptr = mask ? mask->data : NULL;
        T *ptr;
        long int n;
        RFLOAT N = 0, sumx = 0, sumy = 0, sumxx = 0, sumxy = 0;
        FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY_ptr(*this, n, ptr) {
            if (!mask || *maskptr != 0) {
                N++;
                T x = *ptr; T y = *targetptr;
                sumx += x; sumxx += x * x;
                sumy += y; sumxy += x * y;
            }
            targetptr++;
            if (mask) { maskptr++; }
        }
        RFLOAT slope = (N * sumxy - sumx * sumy) / (N * sumxx - sumx * sumx);
        RFLOAT intercept = sumy / N - slope * sumx / N;
        FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY_ptr(*this, n, ptr) {
            // a + b * x
            *ptr = static_cast<RFLOAT>(intercept + slope * static_cast<RFLOAT>(*ptr));
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
    void statisticsAdjust(RFLOAT avgF, RFLOAT stddevF) {

        if (size() == 0) return;

        Stats<T> stats = computeStats();

        RFLOAT a = stats.stddev == 0 ? 0 :
                   static_cast<RFLOAT>(stddevF)  / static_cast<RFLOAT>(stats.stddev);
        RFLOAT b = static_cast<RFLOAT>(avgF) - a * static_cast<RFLOAT>(stats.avg);

        T *ptr;
        long int n;
        FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY_ptr(*this, n, ptr) {
            *ptr = static_cast<T>(a * static_cast<RFLOAT>(*ptr) + b);
        }
    }
    //@}

    /** @name Array "by" array operations.
     *
     * These are operations that are performed between two arrays of the
     * SAME internal type (two integer vectors, two RFLOAT matrices, ...).
     * If they are not of the same type,
     * typeCast should be used to convert one of the arrays to the desired type.
     * The result must have been defined to be of the same type as the operands.
     *
     * In this kind of operations each element of array 1 is operated with its
     * homologous in array 2, it is very important that both have got the
     * same size and starting origins. The result has also got the same
     * shape as the two operated arrays and its former content is lost.
     */
    //@{

    MultidimArray<T> operator + (const MultidimArray<T> &arg) const;

    MultidimArray<T> operator - (const MultidimArray<T> &arg) const;

    MultidimArray<T> operator * (const MultidimArray<T> &arg) const;

    MultidimArray<T> operator / (const MultidimArray<T> &arg) const;

    MultidimArray<T> operator += (const MultidimArray<T> &arg);

    MultidimArray<T> operator -= (const MultidimArray<T> &arg);

    MultidimArray<T> operator *= (const MultidimArray<T> &arg);

    MultidimArray<T> operator /= (const MultidimArray<T> &arg);
    //@}

    /** @name Array "by" scalar operations
     *
     * Operations are between an array and a scalar (of the same type as the array).
     * The result must have been defined to be of the same type as the operands.
     *
     * In this kind of operation each element of array 1 is operated with the given constant.
     * The result has also got the same shape as the input array and its former content is lost
     *
     * Now would be a good time for ad-hoc polymorphism.
     */
    //@{

    MultidimArray<T> operator + (const T scalar) const;

    MultidimArray<T> operator - (const T scalar) const;

    MultidimArray<T> operator * (const T scalar) const;

    MultidimArray<T> operator / (const T scalar) const;

    MultidimArray<T> operator += (const T scalar);

    MultidimArray<T> operator -= (const T scalar);

    MultidimArray<T> operator *= (const T scalar);

    MultidimArray<T> operator /= (const T scalar);
    //@}

    /** @name Scalar "by" array operations
     *
     * These operations are between a scalar (of the same type as the array)
     * and an array. The result must have been defined to be of the same type
     * as the operand. The former content of the result array is lost after
     * the operation.
     *
     * In this kind of operations the constant is operated with each element
     * of array 2. The result has also got the same shape as the input array
     * and its former content is lost
     */
    //@{

    template <typename A>
    friend MultidimArray<A> operator + (const A scalar, const MultidimArray<A> &input);

    template <typename A>
    friend MultidimArray<A> operator - (const A scalar, const MultidimArray<A> &input);

    template <typename A>
    friend MultidimArray<A> operator * (const A scalar, const MultidimArray<A> &input);

    template <typename A>
    friend MultidimArray<A> operator / (const A scalar, const MultidimArray<A> &input);
    //@}

    /// @name Initialization
    /// @{

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
        T *ptr;
        long int n;
        FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY_ptr(*this, n, ptr) { *ptr = val; }
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
    void initZeros(const MultidimArray<T1> &op) {
        if (!data || !sameShape(op)) { reshape(op); }
        memset(data, 0, size() * sizeof(T));
    }

    template <typename T2>
    static MultidimArray<T> zeros(const MultidimArray<T2> &arr1) {
        MultidimArray<T> arr2(arr1);
        arr2.initZeros();
        return arr2;
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

    /** Initialize to zeros with a given size.
     */
    void initZeros(long int Xdim) {
        initZeros(1, 1, 1, Xdim);
    }

    static MultidimArray<T> zeros(long int Xdim) {
        MultidimArray<T> arr(Xdim);
        arr.initZeros();
        return arr;
    }

    /** Initialize to zeros with a given size.
     */
    void initZeros(long int Ydim, long int Xdim) {
        initZeros(1, 1, Ydim, Xdim);
    }

    static MultidimArray<T> zeros(long int Ydim, long int Xdim) {
        MultidimArray<T> arr(Ydim, Xdim);
        arr.initZeros();
        return arr;
    }

    /** Initialize to zeros with a given size.
     */
    void initZeros(long int Zdim, long int Ydim, long int Xdim) {
        initZeros(1, Zdim, Ydim, Xdim);
    }

    static MultidimArray<T> zeros(long int Zdim, long int Ydim, long int Xdim) {
        MultidimArray<T> arr(Zdim, Ydim, Xdim);
        arr.initZeros();
        return arr;
    }

    /** Initialize to zeros with a given size.
     */
    inline void initZeros(long int Ndim, long int Zdim, long int Ydim, long int Xdim) {
        if (xdim != Xdim || ydim != Ydim || zdim != Zdim || ndim != Ndim)
            reshape(Xdim, Ydim, Zdim, Ndim);
        memset(data, 0, size() * sizeof(T));
    }

    static MultidimArray<T> zeros(long int Ndim, long int Zdim, long int Ydim, long int Xdim) {
        MultidimArray<T> arr(Xdim, Ydim, Zdim, Ndim);
        arr.initZeros();
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
     * v1.initLinear(1, 3); // v1=[1 2 3]
     * v1.initLinear(1.5, 3.1); // v1=[1.5 2.5]
     * v1.initLinear(0, 10, 3); // v1=[0 3 6 9]
     * v1.initLinear(0, 10, 3, "incr"); // v1=[0 3 6 9]
     * @endcode
     *
     * Step functionality: The given range is divided in as many points as
     * indicated (in the example 6 points).
     *
     * @code
     * v1.initLinear(0, 10, 6, "steps"); // v1=[0 2 4 6 8 10]
     * @endcode
     */
    void initLinear(T minF, T maxF, int n = 1, const std::string& mode = "incr") {
        RFLOAT slope;
        int steps;

        if (mode == "incr") {
            steps = 1 + floor((RFLOAT) abs(maxF - minF) / (RFLOAT) n);
            slope = n * sgn(maxF - minF); // maxF and minF should not be equal
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
     * v.initRandom(0, 1);
     * // uniform distribution between 0 and 1
     *
     * v.initRandom(0, 1, "uniform");
     * // the same
     *
     * v.initRandom(0, 1, "gaussian");
     * // gaussian distribution with 0 mean and stddev=1
     * @endcode
     */
    void initRandom(RFLOAT op1, RFLOAT op2, const std::string &mode = "uniform") {
        T *ptr = NULL;
        long int n;
               if (mode == "uniform") {
            FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY_ptr(*this, n, ptr)
            { *ptr = static_cast<T>(rnd_unif(op1, op2)); }
        } else if (mode == "gaussian") {
            FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY_ptr(*this, n, ptr)
            { *ptr = static_cast<T>(rnd_gaus(op1, op2)); }
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
        T *ptr = NULL;
        unsigned long int n;
               if (mode == "uniform") {
            FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY_ptr(*this, n, ptr)
            { *ptr += static_cast<T>(rnd_unif(op1, op2)); }
        } else if (mode == "gaussian") {
            FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY_ptr(*this, n, ptr)
            { *ptr += static_cast<T>(rnd_gaus(op1, op2)); }
        } else if (mode == "student") {
            FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY_ptr(*this, n, ptr)
            { *ptr += static_cast<T>(rnd_student_t(df, op1, op2)); }
        } else {
            REPORT_ERROR( static_cast< std::string >("AddNoise: Mode not supported (" + mode + ")"));
        }
    }
    //@}

    /** @name Utilities
     *
     * Here you have several easy functions to manage the values of
     * the array.
     */
    //@{

    /** Produce a 3D array suitable for working with Numerical Recipes.
     *
     * This function must be used only as a preparation for routines which need
     * that the first physical index is 1 and not 0 as it usually is in C. New
     * memory is needed to hold the new RFLOAT pointer array.
     */
    T*** adaptForNumericalRecipes3D(long int n = 0) const {
        T ***m = NULL;
        ask_Tvolume(m, 1, zdim, 1, ydim, 1, xdim);

        FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY3D(*this)
            m[i + 1][j + 1][k + 1] = direct::elem(*this, i, j, k, n);

        return m;
    }

    /** Kill a 3D array produced for numerical recipes.
     */
    void killAdaptationForNumericalRecipes3D(T ***m) const {
        free_Tvolume(m, 1, zdim, 1, ydim, 1, xdim);
    }

    /** Produce a 2D array suitable for working with Numerical Recipes
     *
     * This function must be used only as a preparation for routines which need
     * that the first physical index is 1 and not 0 as it usually is in C. New
     * memory is needed to hold the new RFLOAT pointer array.
     */
    T** adaptForNumericalRecipes2D(long int n = 0) const {
        T **m = NULL;
        ask_Tmatrix(m, 1, ydim, 1, xdim);

        FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY2D(*this)
            m[i + 1][j + 1] = direct::elem(*this, i, j, 0, n);

        return m;
    }

    /** Produce a 1D pointer suitable for working with Numerical Recipes (2)
     *
     * This function meets the same goal as the one before, however this one
     * work with 2D arrays as a single pointer. The first element of the array
     * is pointed by result[1*Xdim+1], and in general result[i*Xdim+j]
     */
    T* adaptForNumericalRecipes22D() const { return data - 1 - xdim; }

    /** Load 2D array from numerical recipes result.
     */
    void loadFromNumericalRecipes2D(T **m, long int Ydim, long int Xdim) {
        resize(Xdim, Ydim);

        for (long int i = 1; i <= Ydim; i++)
        for (long int j = 1; j <= Xdim; j++)
        { (*this)(i - 1, j - 1) = m[i][j]; }
    }

    /** Kill a 2D array produced for numerical recipes
     *
     * The allocated memory is freed.
     */
    void killAdaptationForNumericalRecipes2D(T** m) const {
        free_Tmatrix(m, 1, ydim, 1, xdim);
    }

    /** Kill a 2D array produced for numerical recipes, 2.
     *
     * Nothing needs to be done.
     */
    void killAdaptationForNumericalRecipes22D(T** m) const {}

    /** Produce a 1D array suitable for working with Numerical Recipes
     *
     * This function must be used only as a preparation for routines which need
     * that the first physical index is 1 and not 0 as it usually is in C. In
     * fact the vector provided for Numerical recipes is exactly this same one
     * but with the indices changed.
     *
     * This function is not ported to Python.
     */
    T* adaptForNumericalRecipes1D() const { return data - 1; }

    /** Kill a 1D array produced for Numerical Recipes.
     *
     * Nothing needs to be done in fact.
     *
     * This function is not ported to Python.
     */
    void killAdaptationForNumericalRecipes1D(T* m) const {}

    /** Computes the center of mass of the nth array
     */
    void centerOfMass(Matrix1D<RFLOAT> &center, void *mask = NULL, long int n = 0) {
            center.initZeros(3);
            MultidimArray<int> *imask = (MultidimArray<int>*) mask;

            RFLOAT mass = 0;
            FOR_ALL_ELEMENTS_IN_ARRAY3D(*this) {
                if ((!imask || imask->elem(i, j, k, n)) && elem(i, j, k) > 0) {
                XX(center) += j * elem(i, j, k, n);
                YY(center) += i * elem(i, j, k, n);
                ZZ(center) += k * elem(i, j, k, n);

                mass += elem(i, j, k, n);
            }
        }

        if (mass != 0) { center /= mass; }

        // Resize center to the correct dimensionality
        int dim = getDim();
        if (dim == 1 || dim == 2) { center.resize(dim); }

    }

    /** Sort 1D vector elements
     *
     * Sort in ascending order the vector elements. You can use the "selfReverseX"
     * function to sort in descending order.
     *
     * @code
     * v1.sort();
     * @endcode
     */
    void sort() {
        checkDimension(1);
        std::sort(data, data + xdim);
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
        // Set up a vector of pairs
        std::vector<std::pair<T, long int> > vp;
        vp.reserve(xdim);
        for (long int n = 0; n < size(); n++) {
            vp.push_back(std::make_pair((*this)[n], n));
        }
        // Sort on the first elements of the pairs
        std::sort(vp.begin(), vp.end());
        idx.resize(xdim);
        // Fill the output array with the second elements of the sorted vp
        for (long int n = 0; n < idx.size(); n++) {
            idx[n] = vp[n].second;
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
    void threshold(
        const std::string &type, T a, T b,
        MultidimArray<int> *mask = NULL
    ) {

        int mode =
            type == "abs_above" ?  1 :
            type == "abs_below" ?  2 :
            type == "above"     ?  3 :
            type == "below"     ?  4 :
            type == "range"     ?  5 :
                                   0;

        if (mode == 0)
            REPORT_ERROR(static_cast<std::string>("Threshold: mode not supported (" + type + ")"));

        T *ptr;
        long int n;
        FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY_ptr(*this, n, ptr) {
            // Hopefully the compiler will hoist this loop-invariant switch block
            if (!mask || (*mask)[n] > 0) {
                switch (mode) {

                    case 1: if (abs(*ptr) > a) { *ptr = b * sgn(*ptr); } break;

                    case 2: if (abs(*ptr) < a) { *ptr = b * sgn(*ptr); } break;

                    case 3: if (*ptr > a) { *ptr = b; } break;

                    case 4: if (*ptr < a) { *ptr = b; } break;

                    case 5: if (*ptr < a) { *ptr = a; } else
                            if (*ptr > b) { *ptr = b; } break;

                }
            }
        }
    }

    /** Count with threshold.
     *
     * Return the number of elements meeting the threshold
     * condition.
     */
    long int countThreshold(const std::string& type, T a, T b, MultidimArray<int> *mask = NULL) {
        int mode =
            type == "abs_above" ? 1 :
            type == "abs_below" ? 2 :
            type == "above"     ? 3 :
            type == "below"     ? 4 :
            type == "range"     ? 5 :
                                  0;
        if (!mode) REPORT_ERROR(static_cast<std::string>("CountThreshold: mode not supported (" + type + ")"));

        long int ret = 0;

        T *ptr;
        long int n;
        FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY_ptr(*this, n, ptr) {
            if (!mask || (*mask)[n] > 0) {
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
        T oldv, T newv, RFLOAT accuracy = Xmipp::epsilon,
        MultidimArray<int> *mask = NULL
    ) {
        T *ptr;
        long int n;
        FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY_ptr(*this, n, ptr) {
            if (
                (!mask || (*mask)[n] > 0)
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
        T oldv, T avgv, T sigv, RFLOAT accuracy = Xmipp::epsilon,
        MultidimArray<int> *mask = NULL
    ) {
        T *ptr;
        long int n;
        FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY_ptr(*this, n, ptr) {
            if (
                (!mask || (*mask)[n] > 0)
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
        RFLOAT val = 0, RFLOAT accuracy = Xmipp::epsilon,
        MultidimArray<int> *mask = NULL
    ) {
        T *ptr;
        long int n;
        FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY_ptr(*this, n, ptr) {
            if (!mask || (*mask)[n] > 0) {
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
        T *ptr;
        long int n;
        FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY_ptr(*this, n, ptr)
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
        T *ptr;
        long int n;
        FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY_ptr(*this, n, ptr)
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
        T *ptr;
        long int n;
        FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY_ptr(*this, n, ptr)
        sum += *ptr * *ptr;
        return sum;
    }

    /** Log10.
     *
     * Each component of the result is the log10 of the original components.
     */
    void selfLog10() {
        T *ptr;
        long int n;
        FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY_ptr(*this, n, ptr)
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

    /** Extracts the 1D profile between two points in a 2D array
     *
     * Given two logical indices,
     * use bilinear interpolation to return N samples of the line between them.
     */
    void profile(
        long int x0, long int y0, long int xF, long int yF, long int N,
        MultidimArray<RFLOAT> &profile
    ) const {
        checkDimension(2);
        profile.initZeros(N);
        RFLOAT tx_step = (RFLOAT)(xF - x0) / (N - 1);
        RFLOAT ty_step = (RFLOAT)(yF - y0) / (N - 1);
        RFLOAT tx = x0, ty = y0;

        for (long int i = 0; i < N; i++) {
            profile(i) = interpolatedElement2D(tx, ty);
            tx += tx_step;
            ty += ty_step;
        }
    }

    /** Show using gnuplot
     *
     * This function uses gnuplot to plot this vector. You must supply the
     * xlabel, ylabel, and title.
     */
    void showWithGnuPlot(const std::string& xlabel, const std::string& title) {
        checkDimension(1);

        FileName fn_tmp;
        fn_tmp.initRandom(10);
        MultidimArray<T>::write(static_cast<std::string>("PPP") + fn_tmp + ".txt");

        std::ofstream fh_gplot((static_cast<std::string>("PPP") + fn_tmp + ".gpl").c_str());
        if (!fh_gplot)
            REPORT_ERROR(static_cast<std::string>("vector::showWithGnuPlot: Cannot open PPP") + fn_tmp + ".gpl for output");
        fh_gplot << "set xlabel \"" + xlabel + "\"\n";
        fh_gplot << "plot \"PPP" + fn_tmp + ".txt\" title \"" + title +
        "\" w l\n";
        fh_gplot << "pause 300 \"\"\n";
        system((static_cast<std::string>("(gnuplot PPP") + fn_tmp + ".gpl; rm PPP" + fn_tmp + ".txt PPP" + fn_tmp + ".gpl) &").c_str());
    }

    /** Edit with xmipp_editor.
     *
     * This function generates a random filename starting with PPP and
     * edits it with xmipp_editor. After closing the editor the file is
     * removed.
     */
    void edit() {
        FileName nam;
        nam.initRandom(15);

        nam = static_cast<std::string>("PPP" + nam + ".txt");
        write(nam);

        system((static_cast<std::string>("xmipp_edit -i " + nam + " -remove &").c_str()));
    }

    /* Write to a binary file
    */
    void writeBinary(const FileName &fn) const {
        std::ofstream ofs(fn.c_str(), std::ios::out | std::ios::binary);
        if (!ofs)
            REPORT_ERROR(static_cast<std::string>("MultidimArray::write: File " + fn + " cannot be opened for output"));

        T *ptr;
        unsigned long int n;
        FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY_ptr(*this, n, ptr)
            ofs.write(reinterpret_cast<char*>(ptr), sizeof(T));
    }

    /** Read from a binary file.
     * The array must be previously resized to the correct size.
     */
    void readBinary(const FileName &fn) {
        std::ifstream ifs(fn.c_str(), std::ios::in | std::ios::binary);
        if (!ifs)
            REPORT_ERROR(static_cast<std::string>("MultidimArray::read: File " + fn + " not found"));

        T *ptr;
        unsigned long int n;
        FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY_ptr(*this, n, ptr)
            ifs.read(reinterpret_cast<char*>(ptr), sizeof(T));

    }

    /** Read from a binary file, while summing to the existing array
     * The array must be previously resized to the correct size.
     */
    void readBinaryAndSum(const FileName &fn) {
        std::ifstream ifs(fn.c_str(), std::ios::in | std::ios::binary);
        if (!ifs)
            REPORT_ERROR(static_cast<std::string>("MultidimArray::read: File " + fn + " not found"));

        T *ptr;
        unsigned long int n;
        FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY_ptr(*this, n, ptr) {
            T val;
            ifs.read(reinterpret_cast<char*>(&val), sizeof(T));
            *ptr += val;
        }
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

    /** Assignment.
     *
     * You can build as complex assignment expressions as you like.
     * Multiple assignment is allowed.
     *
     * @code
     * v1 = v2 + v3;
     * v1 = v2 = v3;
     * @endcode
     *
     * This function is ported to Python as assign.
     */
    MultidimArray<T> &operator = (const MultidimArray<T> &op1) {
        if (&op1 != this) {
            if (!data || !sameShape(op1)) resize(op1);
            memcpy(data, op1.data, op1.size() * sizeof(T));
        }
        return *this;
    }

    /** Unary minus.
     *
     * It is used to build arithmetic expressions. You can make a minus
     * of anything as long as it is correct semantically.
     *
     * @code
     * v1 = -v2;
     * v1 = -v2.transpose();
     * @endcode
     */
    MultidimArray<T> operator - () const {
        MultidimArray<T> tmp(*this);
        T *ptr;
        long int n;
        FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY_ptr(tmp, n, ptr) { *ptr = -(*ptr); }
        return tmp;
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
        T *ptr;
        long int n;
        FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY_ptr(v, n, ptr) { in >> *ptr; }
        return in;
    }

    /** Equality.
     *
     * Do these two objects have the same shape (origin and size)
     * and the same values (to within accuracy)?
     */
    bool equal(const MultidimArray<T> &other, RFLOAT accuracy = Xmipp::epsilon) const {

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
 * If we have an integer array and we need a RFLOAT one, we can use this
 * function. The conversion is done through a type casting of each element
 * If n >= 0, only the nth volumes will be converted, otherwise all ndim volumes
 */
template<typename T1, typename T2>
void typeCast(const MultidimArray<T1>& v1,  MultidimArray<T2>& v2, long n = -1) {
    if (v1.size() == 0) {
        v2.clear();
        return;
    }

    if (n < 0) {
        v2.resize(v1);
        T1 *ptr1;
        long int n;
        FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY_ptr(v1, n, ptr1) {
            v2[n] = static_cast<T2>(*ptr1);
        }
    } else {
        v2.resize(v1.xdim, v1.ydim, v1.zdim);
        FOR_ALL_DIRECT_ELEMENTS_IN_ARRAY3D(v2) {
            direct::elem(v2, i, j, k) = static_cast<T2>(direct::elem(v1, i, j, k, n));
        }
    }
}

/** Force positive.
 *  Apply a median filter to negative values. Leave positive values.
 */
void forcePositive(MultidimArray<RFLOAT> &V);

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
    T *ptr;
    long int n;
    FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY_ptr(v, n, ptr) {
        max_val = std::max(max_val, (T) std::abs(*ptr));
    }

    int prec = bestPrecision(max_val, 10);

    if (v.ydim == 1 && v.zdim == 1) {
        for (long int j = v.firstX(); j <= v.lastX(); j++) {
            ostrm << floatToString((RFLOAT) v.elem(0, j, 0), 10, prec)
            << std::endl;
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
#endif
