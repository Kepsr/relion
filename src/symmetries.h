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
/* ------------------------------------------------------------------------- */
/* SYMMETRIES                                                                */
/* ------------------------------------------------------------------------- */
#ifndef _SYMMETRIES_HH
#define _SYMMETRIES_HH

#include "src/matrix1d.h"
#include "src/matrix2d.h"
#include "src/euler.h"
#include "src/funcs.h"
#include "src/args.h"

/**
 * @defgroup SymmetryLists Symmetry handling
 * @ingroup DataLibrary
 *
 *  The symmetry lists are, simply, lists of 2D matrices. It's the way of
 *  taking symmetry into account in the reconstruction programs. The
 *  symmetry list must contain matrices which express equivalent views to
 *  the actual one due to the underlying volume symmetry. The identity matrix
 *  is not within the list. You know that symmetry matrices should form
 *  a subgroup, when reading a file the subgroup is automatically computed
 *  and, when you add or remove a new matrix, the subgroup must be
 *  manually computed.
 */

// @{
// point group symmetries
namespace pg { enum {

    // 200..225
    CI = 200,  // Point group Ci (inversion symmetry)
    CS,        // Point group Cs (reflection symmetry)
    CN,        // Point group series Cn (n-fold rotational symmetry)
    CNV,       // Point group series Cnv (pyramidal symmetry)
    CNH,       // Point group series Cnh (prismatic symmetry)
    SN,        // Point group series S2n (rotoreflection symmetry)
    DN,        // Point group series Dn (dihedral symmetry)
    DNV,       // Point group series Dnv
    DNH,       // Point group series Dnh
    T,         // Point group T (chiral tetrahedral symmetry)
    TD,        // Point group Td (full tetrahedral symmetry)
    TH,        // Point group Th (pyritohedral symmetry)
    O,         // Point group O (chiral octahedral symmetry)
    OH,        // Point group Oh (full octahedral symmetry)
    I,         // Point group I (chiral icosahedral symmetry) [default xmipp icosahedaral symmetry]
    IH,        // Point group Ih (full icosahedral symmetry)
    I1,        // no crowther 222
    I2,        // crowther 222-> default in xmipp
    I3,        // 52 as used by spider
    I4,        // another 52
    I5,        // another another 52 (used by EMBL-matfb)
    I1H,       // no crowther 222, + mirror plane
    I2H,       // crowther 222-> default in xmipp+ mirror plane
    I3H,       // 52 as used by spider+ mirror plane
    I4H,       // another 52+ mirror plane
    I5H,       // another another 52 (used by EMBL-matfb)+ mirror plane

}; }

/** Number of an image in the reconstruction list.
    This macro returns the index of a symmetry image (after the symmetry matrix
    number sym_no) within a list where the first images are true images and the
    last ones, the symmetrized copies (all copies of a same image are
    together). The total number of real images is numIMG, and i is the index
    within this first numIMG images of the image we want to symmetrize The
    first image in the list is the number 0 */
#define SYMINDEX(SL, sym_no, i, numIMG) \
    numIMG + SL.__L.nrows() / 4 * i + sym_no

/** Symmetry List class.
    Internally the symmetry list class is implemented as a single 2D matrix,
    where every 4 rows (remember that in 3D the geometrical transformation
    matrices are 4x4) comprise a symmetry matrix. Access, and ways to modify
    the symmetry list are supplied. Remind that any symmetry is expressed
    in terms of two matrices L and R, so that any Euler matrix must be
    transformed by L*Euler*R resulting into a new perspective of the volume
    which is equivalent to the original one.

    The typical use of the symmetry lists is to read the symmetry file, and
    do nothing else but reading matrices from it.

    The symmetry file format is
    @code
    # This is a comment
    # The following line is a 6-fold rotational symmetry axis along Z-axis.
    # The fold is the number of times that the volume can be rotated along
    # the symmetry axis giving the same view from different view points.
    # the structure for the rotational axis is
    # rot_axis      <fold> <X0> <Y0> <Z0>
    # mirror_plane         <X0> <Y0> <Z0>
    rot_axis      6 0 0 1
    mirror_plane    0 0 1
    @endcode
*/
class SymList {

    public:
    // L and R matrices
    Matrix2D<RFLOAT> __L, __R;
    Matrix1D<int>    __chain_length;

    // As the symmetry elements form a subgroup, this is the number of
    // true symmetry elements belonging to the list, the rest of
    // the list are simply the elements to fill the subgroup
    int true_symNo;

    // Number of Axis, mirrors, ...
    int __sym_elements;

    public:

    /** Create an empty list.
     *  The 2D matrices are 0 × 0.
     *  @code
     *  SymList SL;
     *  @endcode
     */
    SymList() {
        __sym_elements = true_symNo = 0;
    }

    /** Translate a string 'symstring' to a symmetry group.
     *  Return value indicates whether the string was recognised as a symmetry group.
     *  See URL http://xmipp.cnb.uam.es/twiki/bin/view/Xmipp/Symmetry for details.
     *  Symmetry information is output to pgGroun and pgOrder.
     */
    bool isSymmetryGroup(const std::string &symstring, int &pgGroup, int &pgOrder);

    /** fill fileContent with symmetry information */
    void fill_symmetry_class(
        const FileName symmetry, int pgGroup, int pgOrder,
        std::vector<std::string> &fileContent
    );

    /** Create a Symmetry List from a Symmetry file.
     *  All the subgroup elements are computed automatically.
     *  @code
     *  SymList SL("sym.txt");
     *  @endcode
     */
    SymList(const FileName& fn_sym) {
        read_sym_file(fn_sym);
    }

    /** Get matrices from the symmetry list.
     *  The number of matrices inside the list is given by SymsNo.
     *  This function return the 4x4 transformation matrices associated to
     *  the one in the list which occupies the position 'i'. The matrix
     *  numbering within the list starts at 0. The output transformation
     *  matrices is given as a pointer to gain speed.
     *  @code
     *  for (i = 0; i < SL.SymsNo(); i++) {
     *      SL.get_matrices(i, L, R);
     *      ...
     *  }
     *  @endcode
     */
    void get_matrices(int i, Matrix2D<RFLOAT> &L, Matrix2D<RFLOAT> &R) const;

    /** Set a couple of matrices in the symmetry list.
        The number of matrices inside the list is given by SymsNo.
        This function sets the 4x4 transformation matrices associated to
        the one in the list which occupies the position 'i'. The matrix
        numbering within the list starts at 0.
        @code
           for (int i = 0; i < SL.SymsNo(); i++) {
               SL.set_matrix(i, L, R);
               ...
           }
        @endcode */
    void set_matrices(int i, const Matrix2D<RFLOAT> &L, const Matrix2D<RFLOAT> &R);

    /** Read a symmetry file into a symmetry list.
        The former symmetry list is overwritten with the new one.
        All the subgroup members are added to the list.
        If the accuracy is negative then the subgroup is not generated.
        return symmetry group
        @code
        SL.read_sym_file("sym.txt");
        @endcode
     */
    int read_sym_file(FileName fn_sym);

    /** Add symmetry matrices to the symmetry list.
        The given matrix must specify a point of view equivalent to the
        actual point of view. The matrices are added to the subgroup generator
        but the subgroup is not updated, you must do it manually using
        compute_subgroup. What is more, the subgroup after the insertion
        is corrupted.

        The chain length is the number of single matrices multiplication of
        which the inserted one is compound.*/
    void add_matrices(const Matrix2D<RFLOAT> &L, const Matrix2D<RFLOAT> &R, int chain_length);

    /** Compute subgroup for this structure.
        After adding or setting a matrix, the subgroup information
        is lost, you must recalculate it using this function. The different
        matrices are multiplied until no more different matrices are produced.
        The accuracy is used in order to compare when two matrix elements are
        the same.

        So far, all the shifts associated to generated matrices are set to 0*/
    void compute_subgroup();

    /** Number of symmetry matrices inside the structure.
        This is the number of all the matrices inside the subgroup.
        @code
           for (i = 0; i < SL.SymsNo(); i++) {
               SL.get_matrix(i, A);
               ...
           }
        @endcode */
    int SymsNo() const {
        return __L.nrows() / 4;
    }

    /** Number of symmetry matrices which generated the structure.
        This is the number of the matrices which generated the structure,
        notice that it should be always less or equal to the total number
        of matrices in the subgroup. */
    int TrueSymsNo() const {
        return true_symNo;
    }

    // Write the symmetry definition (plus all rotation matrices) to an output stream.
    void writeDefinition(std::ostream &outstream, FileName fn_sym);

    // Return the area of the non redundant part of the Ewald sphere
    RFLOAT non_redundant_ewald_sphere(int pgGroup, int pgOrder);

};


// Symmetrise a 3D map according to the specified symmetry
void symmetriseMap(MultidimArray<RFLOAT> &img, FileName &fn_sym, bool do_wrap = false);

//@}
#endif
