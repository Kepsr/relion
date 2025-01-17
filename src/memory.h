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

#ifndef _XMIPP_MEMORY
#define _XMIPP_MEMORY

#include "src/error.h"

/* Memory managing --------------------------------------------------------- */
///@defgroup MemoryManaging Memory management for numerical recipes
/// @ingroup DataLibrary
//@{
/** Ask memory for any type vector.
    The valid values range from v[nl] to v[nh]. If no memory is available
    an exception is thrown. nullptr is returned if nh is not greater than nl*/
template <typename T>
T* ask_vector(int nl, int nh) {
    if (nh - nl + 1 <= 1) return nullptr;
    T* v = (T*) malloc((unsigned) (nh - nl + 1) * sizeof(T));
    if (!v) REPORT_ERROR("allocation failure in vector()");
    return v - nl;
}

/** Free memory associated to any type vector. */
template <typename T>
void free_vector(T* v, int nl, int nh) {
    if (v) free((char*) (v + nl));
}

/** Ask memory for any type matrix.
    The valid values range from v[nrl][ncl] to v[nrh][nch].
    If no memory is available an exception is thrown. NULL is returned if any
    nh is not greater than its nl*/
template <typename T>
T** ask_matrix(int nrl, int nrh, int ncl, int nch) {

    const int a = nrh - nrl + 1, b = nch - ncl + 1;
    if (a <= 1 || b <= 1) return nullptr;

    T** m = (T**) malloc((unsigned) a * sizeof(T*));
    if (!m)
        REPORT_ERROR("allocation failure 1 in matrix()");
    m -= nrl;

    for (int i = nrl; i <= nrh; i++) {
        m[i] = (T*) malloc((unsigned) b * sizeof(T));
        if (!m[i])
            REPORT_ERROR("allocation failure 2 in matrix()");
        m[i] -= ncl;
    }
    return m;
}

/** Free memory associated to any type matrix.
    After freeing v=NULL*/
template <class T>
void free_matrix(T** m, int nrl, int nrh, int ncl, int nch) {
    if (m) {
        for (int i = nrh; i >= nrl; i--)
            free((char*) (m[i] + ncl));
        free((char*) (m + nrl));
    }
}

/** Ask memory for any type volume.
    The valid values range from v[nsl][nrl][ncl] to v[nsh][nrh][nch].
    If no memory is available an exception is thrown. NULL is returned if any
    nh is not greater than its nl. */
template <class T> void ask_Tvolume(
    T *** &m, int nsl, int nsh, int nrl, int nrh, int ncl, int nch
) {
    if (nsh - nsl + 1 > 1 && nrh - nrl + 1 > 1 && nch - ncl + 1 > 1) {
        m = (T ***) malloc((unsigned) (nsh - nsl + 1) * sizeof(T**));
        if (!m)
            REPORT_ERROR("allocation failure 1 in matrix()");
        m -= nsl;

        for (int k = nsl; k <= nsh; k++) {
            m[k] = (T **) malloc((unsigned) (nrh - nrl + 1) * sizeof(T*));
            if (!m[k])
                REPORT_ERROR("allocation failure 2 in matrix()");
            m[k] -= nrl;

            for (int i = nrl;i <= nrh;i++) {
                m[k][i] = (T *) malloc((unsigned) (nch - ncl + 1) * sizeof(T));
                if (!m[k][i])
                    REPORT_ERROR("allocation failure 2 in matrix()");
                m[k][i] -= ncl;
            }
        }
    }
    else m = NULL;
}

/** Free memory associated to any type volume.
    After freeing v=NULL*/
template <class T> void free_Tvolume(
    T *** &m, int nsl, int nsh, int nrl, int nrh, int ncl, int nch
) {
    if (m) {
        for (int k = nsh; k >= nsl; k--) {
        for (int i = nrh; i >= nrl; i--)
        free((char*) (m[k][i] + ncl));
        free((char*) (m[k] + nrl));
        }
        free((char*) (m + nsl));
        m = NULL;
    }
}

template <typename T>
struct callocator {

    /** Allocate memory
     *
     * A thin wrapper around calloc.
     * returns a T* pointing to allocated memory.
     * Successfully allocated memory is zero-initialised.
     * If calloc returns nullptr, an exception is thrown.
     *
     */
    static T* allocate(size_t memsize) {
        T *ptr = (T *) calloc(memsize, sizeof(T));
        if (!ptr) {
            std::cerr << "Failed to allocate " <<  memsize << " bytes." << std::endl;
            REPORT_ERROR("Error in callocator::allocate");
        }
        return ptr;
    }

    /** Free allocated memory
     * Adapted from Bsofts bfree
     *
     * A thin wrapper around free.
     * Frees the memory pointed to by ptr.
     * On success, returns 0.
     * If memsize is negative, does nothing and returns -1.
     *
    */
    static int deallocate(T *ptr, size_t memsize) {
        if (memsize < 1) return -1;
        free(ptr);
        return 0;
    }

};

//@}
#endif
