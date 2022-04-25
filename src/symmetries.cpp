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

#include <stdio.h>

#include "src/symmetries.h"

// Read Symmetry file ======================================================
// Crystal symmetry matrices from http://cci.lbl.gov/asu_gallery/
int SymList::read_sym_file(FileName fn_sym) {

    // Open file ---------------------------------------------------------
    FILE *fpoii;
    std::vector<std::string> fileContent;
    char line[80];
    int pgGroup = 0, pgOrder = 0;
    if (!(fpoii = fopen(fn_sym.c_str(), "r"))) {
        // Check if reserved word and return group and order
        if (isSymmetryGroup(fn_sym.removeDirectories(), pgGroup, pgOrder)) {
            fill_symmetry_class(fn_sym, pgGroup, pgOrder, fileContent);
        } else {
            REPORT_ERROR(
                (std::string) "SymList::read_sym_file: Can't open file or do not recognize symmetry group " + fn_sym
            );
        }
    } else {
        while (fgets(line, 79, fpoii)) {
            if (line[0] == ';' || line[0] == '#' || line[0] == '\0') continue;
            fileContent.push_back(line);
        }
        fclose(fpoii);
    }

    // Count the number of symmetries ------------------------------------
    true_symNo = 0;
    // count number of axis and mirror planes. It will help to identify
    // the crystallographic symmetry

    int no_axis, no_mirror_planes, no_inversion_points;
    no_axis = no_mirror_planes = no_inversion_points = 0;

    for (int n = 0; n < fileContent.size(); n++) {
        strcpy(line, fileContent[n].c_str());
        char *token = firstToken(line);
        if (!token) {
            std::cout << line;
            std::cout << "Wrong line in symmetry file, the line is skipped\n";
            continue;
        }
        if (strcmp(token, "rot_axis") == 0) {
            token = nextToken();
            true_symNo += textToInteger(token) - 1;
            no_axis++;
        } else if (strcmp(token, "mirror_plane") == 0) {
            true_symNo++;
            no_mirror_planes++;
        } else if (strcmp(token, "inversion") == 0) {
            true_symNo += 1;
            no_inversion_points = 1;
        }
    }
    // Ask for memory
    __L.resize(4 * true_symNo, 4);
    __R.resize(4 * true_symNo, 4);
    __chain_length.resize(true_symNo);
    __chain_length.initConstant(1);

    // Read symmetry parameters
    int i = 0;
    Matrix1D<RFLOAT> axis(3);
    Matrix2D<RFLOAT> L(4, 4), R(4, 4);
    for (int n = 0; n < fileContent.size(); n++) {
        strcpy(line,fileContent[n].c_str());
        char *token = firstToken(line);
        // Rotational axis ---------------------------------------------------
        if (strcmp(token, "rot_axis") == 0) {
            token = nextToken();
            int fold = textToInteger(token);
            token = nextToken();
            XX(axis) = textToDouble(token);
            token = nextToken();
            YY(axis) = textToDouble(token);
            token = nextToken();
            ZZ(axis) = textToDouble(token);
            RFLOAT ang_incr = 360.0 / fold;
            RFLOAT rot_ang;
            L.initIdentity();
            for (int j = 1, rot_ang = ang_incr; j < fold; j++, rot_ang += ang_incr) {
                rotation3DMatrix(rot_ang, axis, R);
                R.setSmallValuesToZero();
                set_matrices(i++, L, R.transpose());
            }
            __sym_elements++;
            // inversion ------------------------------------------------------
        } else if (strcmp(token, "inversion") == 0) {
            L.initIdentity();
            L(2, 2) = -1;
            R.initIdentity();
            R(0, 0) = -1.0;
            R(1, 1) = -1.0;
            R(2, 2) = -1.0;
            set_matrices(i++, L, R);
            __sym_elements++;
            // mirror plane -------------------------------------------------------------
        } else if (strcmp(token, "mirror_plane") == 0) {
            token = nextToken();
            XX(axis) = textToFloat(token);
            token = nextToken();
            YY(axis) = textToFloat(token);
            token = nextToken();
            ZZ(axis) = textToFloat(token);
            L.initIdentity();
            L(2, 2) = -1;
            Matrix2D<RFLOAT> A;
            alignWithZ(axis, A);
            A = A.transpose();
            R = A * L * A.inv();
            L.initIdentity();
            set_matrices(i++, L, R);
            __sym_elements++;
        }
    }

    compute_subgroup();

    return pgGroup;
}

// Get matrix ==============================================================
void SymList::get_matrices(int i, Matrix2D<RFLOAT> &L, Matrix2D<RFLOAT> &R) const {
    L.initZeros(4, 4);
    R.initZeros(4, 4);
    for (int k = 4 * i; k < 4 * i + 4; k++)
    for (int l = 0;     l < 4;         l++) {
        L(k - 4 * i, l) = __L(k, l);
        R(k - 4 * i, l) = __R(k, l);
    }
}

// Set matrix ==============================================================
void SymList::set_matrices(int i, const Matrix2D<RFLOAT> &L, const Matrix2D<RFLOAT> &R) {
    for (int k = 4 * i; k < 4 * i + 4; k++)
    for (int l = 0;     l < 4;         l++) {
        __L(k, l) = L(k - 4 * i, l);
        __R(k, l) = R(k - 4 * i, l);
    }
}

// Matrix addition ============================================================
void SymList::add_matrices(
    const Matrix2D<RFLOAT> &L, const Matrix2D<RFLOAT> &R, int chain_length
) {

    if (
        MAT_XSIZE(L) != 4 || MAT_YSIZE(L) != 4 || 
        MAT_XSIZE(R) != 4 || MAT_YSIZE(R) != 4
    ) REPORT_ERROR( "SymList::add_matrix: Transformation matrix is not 4x4");

    if (TrueSymsNo() == SymsNo()) {
        __L.resize(MAT_YSIZE(__L) + 4, 4);
        __R.resize(MAT_YSIZE(__R) + 4, 4);
        __chain_length.resize(__chain_length.size() + 1);
    }

    set_matrices(true_symNo, L, R);
    __chain_length(__chain_length.size() - 1) = chain_length;
    true_symNo++;
}

// Compute subgroup ========================================================
bool found_not_tried(
    const Matrix2D<int> &tried, int &i, int &j, int true_symNo
) {
    i = j = 0;
    int n = 0;
    while (n != MAT_YSIZE(tried)) {
        if (tried(i, j) == 0 && !(i >= true_symNo && j >= true_symNo))
            return true;
        if (i != n) {
            // Move downwards
            i++;
        } else {
            // Move leftwards
            j--;
            if (j == -1) {
                n++;
                j = n;
                i = 0;
            }
        }
    }
    return false;
}

// #define DEBUG
void SymList::compute_subgroup() {
    Matrix2D<RFLOAT> I(4, 4);
    I.initIdentity();
    Matrix2D<RFLOAT> L1(4, 4), R1(4, 4), L2(4, 4), R2(4, 4), newL(4, 4), newR(4, 4);
    Matrix2D<int>    tried(true_symNo, true_symNo);
    int i, j;
    int new_chain_length;
    while (found_not_tried(tried, i, j, true_symNo)) {
        tried(i, j) = 1;

        get_matrices(i, L1, R1);
        get_matrices(j, L2, R2);
        newL = L1 * L2;
        newR = R1 * R2;
        new_chain_length = __chain_length(i) + __chain_length(j);
        Matrix2D<RFLOAT> newR3 = newR;
        newR3.resize(3, 3);
        if (newL.isIdentity() && newR3.isIdentity()) continue;

        // Try to find it in current ones
        bool found = false;
        for (int l = 0; l < SymsNo(); l++) {
            get_matrices(l, L1, R1);
            if (newL.equal(L1) && newR.equal(R1)) {
                found = true;
                break;
            }
        }

        if (!found) {
            // #define DEBUG
            #ifdef DEBUG
            std::cout << "Matrix size " << tried.Xdim() << " "
            << "trying " << i << " " << j << " "
            << "chain length=" << new_chain_length << std::endl;
            std::cout << "Result R Sh\n" << newR;
            #endif
            // #undef DEBUG
            newR.setSmallValuesToZero();
            newL.setSmallValuesToZero();
            add_matrices(newL, newR, new_chain_length);
            tried.resize(MAT_YSIZE(tried) + 1, MAT_XSIZE(tried) + 1);
        }
    }
}

static inline int ctoi(char c) {
    return int(c) - 48;  // ASCII hack
}

bool SymList::isSymmetryGroup(const std::string &symstring, int &pgGroup, int &pgOrder) {

    int size = symstring.size();

    // The shortest point group names have 1 char:  T/I/O
    // The longest  point group names have 4 chars: Cxxv/Cxxh/Dxxv/Dxxh (where xx is a 2-digit decimal string)
    if (size < 1 || size > 4) {
        pgGroup = -1;
        pgOrder = -1;
        return false;
    }

    // #define DEBUG7
    #ifdef DEBUG7
    #define SUCCESS std::cerr << "pgGroup" << pgGroup << " pgOrder " << pgOrder << std::endl; return true;
    #else
    #define SUCCESS return true;
    #endif
    #undef DEBUG7

    // Copy symstring into a char* and map to upper case
    char Gs[size + 1];
    strcpy(Gs, symstring.c_str()); 
    for (int i = 0; i < size; i++) {
        Gs[i] = toupper(Gs[i]);
    }

    switch (Gs[0]) {

        case 'C':
        if (isdigit(Gs[1]) && size == 2) {
            pgGroup = pg::CN;
            pgOrder = ctoi(Gs[1]);
            SUCCESS
        } else if (isdigit(Gs[1]) && isdigit(Gs[2]) && size == 3) {
            pgGroup = pg::CN;
            pgOrder = atoi(symstring.substr(1, 2).c_str());
            SUCCESS
        } else if (Gs[1] == 'I' && size == 2) {
            pgGroup = pg::CI;
            pgOrder = -1;
            SUCCESS
        } else if (Gs[1] == 'S' && size == 2) {
            pgGroup = pg::CS;
            pgOrder = -1;
            SUCCESS
        } else if (isdigit(Gs[1]) && Gs[2] == 'H' && size == 3) {
            pgGroup = pg::CNH;
            pgOrder = ctoi(Gs[1]);
            SUCCESS
        } else if (isdigit(Gs[1]) && isdigit(Gs[2]) && Gs[3] == 'H' && size == 4) {
            pgGroup = pg::CNH;
            pgOrder = atoi(symstring.substr(1, 2).c_str());
            SUCCESS
        } else if (isdigit(Gs[1]) && Gs[2] == 'V' && size == 3) {
            pgGroup = pg::CNV;
            pgOrder = ctoi(Gs[1]);
            SUCCESS
        } else if (isdigit(Gs[1]) && isdigit(Gs[2]) && Gs[3] == 'V' && size == 4) {
            pgGroup = pg::CNV;
            pgOrder = atoi(symstring.substr(1, 2).c_str());
            SUCCESS
        }

        case 'S':
        if (Gs[0] == 'S' && isdigit(Gs[1]) && size == 2) {
            pgGroup = pg::SN;
            pgOrder = ctoi(Gs[1]);
            SUCCESS
        } else if (Gs[0] == 'S' && isdigit(Gs[1]) && isdigit(Gs[2]) && size == 3) {
            pgGroup = pg::SN;
            pgOrder = atoi(symstring.substr(1, 2).c_str());
            SUCCESS
        }

        case 'D':
        if (isdigit(Gs[1]) && size == 2) {
            pgGroup = pg::DN;
            pgOrder = ctoi(Gs[1]);
            SUCCESS
        } else if (isdigit(Gs[1]) && isdigit(Gs[2]) && size == 3) {
            pgGroup = pg::DN;
            pgOrder = atoi(symstring.substr(1, 2).c_str());
            SUCCESS
        } else if (isdigit(Gs[1]) && Gs[2] == 'V' && size == 3) {
            pgGroup = pg::DNV;
            pgOrder = ctoi(Gs[1]);
            SUCCESS
        } else if (isdigit(Gs[1]) && isdigit(Gs[2]) && Gs[3] == 'V' && size == 4) {
            pgGroup = pg::DNV;
            pgOrder = atoi(symstring.substr(1, 2).c_str());
            SUCCESS
        } else if (isdigit(Gs[1]) && Gs[2] == 'H' && size == 3) {
            pgGroup = pg::DNH;
            pgOrder = ctoi(Gs[1]);
            SUCCESS
        } else if (isdigit(Gs[1]) && isdigit(Gs[2]) && Gs[3] == 'H' && size == 4) {
            pgGroup = pg::DNH;
            pgOrder = atoi(symstring.substr(1, 2).c_str());
            SUCCESS
        }

        case 'T':
        if (size == 1) {
            pgGroup = pg::T;
            pgOrder = -1;
            SUCCESS
        } else if (Gs[1] == 'D' && size == 2) {
            pgGroup = pg::TD;
            pgOrder = -1;
            SUCCESS
        } else if (Gs[1] == 'H' && size == 2) {
            pgGroup = pg::TH;
            pgOrder = -1;
            SUCCESS
        }

        case 'O':
        if (size == 1) {
            pgGroup = pg::O;
            pgOrder = -1;
            SUCCESS
        } else if (Gs[1] == 'H' && size == 2) {
            pgGroup = pg::OH;
            pgOrder = -1;
            SUCCESS
        }
    
        case 'I':
        if (size == 1) {
            pgGroup = pg::I;
            pgOrder = -1;
            SUCCESS
        } else if (Gs[1] == '1' && size == 2) {
            pgGroup = pg::I1;
            pgOrder = -1;
            SUCCESS
        } else if (Gs[1] == '2' && size == 2) {
            pgGroup = pg::I2;
            pgOrder = -1;
            SUCCESS
        } else if (Gs[1] == '3' && size == 2) {
            pgGroup = pg::I3;
            pgOrder = -1;
            SUCCESS
        } else if (Gs[1] == '4' && size == 2) {
            pgGroup = pg::I4;
            pgOrder = -1;
            SUCCESS
        } else if (Gs[1] == '5' && size == 2) {
            pgGroup = pg::I5;
            pgOrder = -1;
            SUCCESS
        } else if (Gs[1] == 'H' && size == 2) {
            pgGroup = pg::IH;
            pgOrder = -1;
            SUCCESS
        } else if (Gs[1] == '1' && Gs[2] == 'H' && size == 3) {
            pgGroup = pg::I1H;
            pgOrder = -1;
            SUCCESS
        } else if (Gs[1] == '2' && Gs[2] == 'H' && size == 3) {
            pgGroup = pg::I2H;
            pgOrder = -1;
            SUCCESS
        } else if (Gs[1] == '3' && Gs[2] == 'H' && size == 3) {
            pgGroup = pg::I3H;
            pgOrder = -1;
            SUCCESS
        } else if (Gs[1] == '4' && Gs[2] == 'H' && size == 3) {
            pgGroup = pg::I4H;
            pgOrder = -1;
            SUCCESS
        } else if (Gs[1] == '5' && Gs[2] == 'H' && size == 3) {
            pgGroup = pg::I5H;
            pgOrder = -1;
            SUCCESS
        }

        default:
        return false;

    }
    #undef SUCCESS
}

void SymList::fill_symmetry_class(
    const FileName symmetry, int pgGroup, int pgOrder, 
    std::vector<std::string> &fileContent
) {
    fileContent.clear();

    switch (pgGroup) {

        case pg::CN:
        fileContent.push_back("rot_axis " + integerToString(pgOrder) + " 0 0 1");
        break;

        case pg::CI:
        fileContent.push_back("inversion ");
        break;

        case pg::CS:
        fileContent.push_back("mirror_plane 0 0 1");
        break;

        case pg::CNV:
        fileContent.push_back("rot_axis " + integerToString(pgOrder) + " 0 0 1");
        fileContent.push_back("mirror_plane 0 1 0");
        break;

        case pg::CNH:
        fileContent.push_back("rot_axis " + integerToString(pgOrder) + " 0 0 1");
        fileContent.push_back("mirror_plane 0 0 1");
        break;

        case pg::SN:
        // i.e. S2n
        if (pgOrder % 2 == 1) {
            std::cerr << "ERROR: order for SN group must be even" << std::endl;
            exit(0);
        }
        fileContent.push_back("rot_axis " + integerToString(pgOrder / 2) + " 0 0 1");
        fileContent.push_back("inversion ");
        break;

        case pg::DN:
        if (pgOrder > 1)
        fileContent.push_back("rot_axis " + integerToString(pgOrder) + " 0 0 1");
        fileContent.push_back("rot_axis 2 1 0 0");
        break;

        case pg::DNV:
        if (pgOrder > 1)
        fileContent.push_back("rot_axis " + integerToString(pgOrder) + " 0 0 1");
        fileContent.push_back("rot_axis 2 1 0 0");
        fileContent.push_back("mirror_plane 1 0 0");
        break;

        case pg::DNH:
        if (pgOrder > 1)
        fileContent.push_back("rot_axis " + integerToString(pgOrder) + " 0 0 1");
        fileContent.push_back("rot_axis 2 1 0 0");
        fileContent.push_back("mirror_plane 0 0 1");
        break;

        case pg::T:
        fileContent.push_back("rot_axis 3  0. 0. 1.");
        fileContent.push_back("rot_axis 2 0. 0.816496 0.577350");
        break;

        case pg::TD:
        fileContent.push_back("rot_axis 3  0. 0. 1.");
        fileContent.push_back("rot_axis 2 0. 0.816496 0.577350");
        fileContent.push_back("mirror_plane 1.4142136 2.4494897 0.0000000");
        break;

        case pg::TH:
        fileContent.push_back("rot_axis 3  0. 0. 1.");
        fileContent.push_back("rot_axis 2 0. -0.816496 -0.577350");
        fileContent.push_back("inversion");
        break;

        case pg::O:
        fileContent.push_back("rot_axis 3  .5773502  .5773502 .5773502");
        fileContent.push_back("rot_axis 4 0 0 1");
        break;

        case pg::OH:
        fileContent.push_back("rot_axis 3  .5773502  .5773502 .5773502");
        fileContent.push_back("rot_axis 4 0 0 1");
        fileContent.push_back("mirror_plane 0 1 1");
        break;

        case pg::I:
        case pg::I2:
        fileContent.push_back("rot_axis 2  0 0 1");
        fileContent.push_back("rot_axis 5  0.525731114  0 0.850650807");
        fileContent.push_back("rot_axis 3  0 0.356822076 0.934172364");
        break;

        case pg::I1:
        fileContent.push_back("rot_axis 2  1  	   0	       0");
        fileContent.push_back("rot_axis 5 0.85065080702670 0 -0.5257311142635");
        fileContent.push_back("rot_axis 3 0.9341723640 0.3568220765 0");
        break;

        case pg::I3:
        fileContent.push_back("rot_axis 2  -0.5257311143 0 0.8506508070");
        fileContent.push_back("rot_axis 5  0. 0. 1.");
        fileContent.push_back("rot_axis 3  -0.4911234778630044, 0.3568220764705179, 0.7946544753759428");
        break;

        case pg::I4:
        fileContent.push_back("rot_axis 2  0.5257311143 0 0.8506508070");
        fileContent.push_back("rot_axis 5  0.8944271932547096 0 0.4472135909903704");
        fileContent.push_back("rot_axis 3  0.4911234778630044 0.3568220764705179 0.7946544753759428");
        break;

        case pg::I5:
        std::cerr << "ERROR: Symmetry pg::I5 not implemented" << std::endl;
        exit(0);

        case pg::IH:
        case pg::I2H:
        fileContent.push_back("rot_axis 2  0 0 1");
        fileContent.push_back("rot_axis 5  0.525731114  0 0.850650807");
        fileContent.push_back("rot_axis 3  0 0.356822076 0.934172364");
        fileContent.push_back("mirror_plane 1 0 0");
        break;

        case pg::I1H:
        fileContent.push_back("rot_axis 2  1  	   0	       0");
        fileContent.push_back("rot_axis 5 0.85065080702670 0 -0.5257311142635");
        fileContent.push_back("rot_axis 3 0.9341723640 0.3568220765 0");
        fileContent.push_back("mirror_plane 0 0 -1");
        break;

        case pg::I3H:
        fileContent.push_back("rot_axis 2  -0.5257311143 0 0.8506508070");
        fileContent.push_back("rot_axis 5  0. 0. 1.");
        fileContent.push_back("rot_axis 3  -0.4911234778630044, 0.3568220764705179, 0.7946544753759428");
        fileContent.push_back("mirror_plane 0.850650807 0  0.525731114");

        break;
        case pg::I4H:
        fileContent.push_back("rot_axis 2  0.5257311143 0 0.8506508070");
        fileContent.push_back("rot_axis 5  0.8944271932547096 0 0.4472135909903704");
        fileContent.push_back("rot_axis 3  0.4911234778630044 0.3568220764705179 0.7946544753759428");
        fileContent.push_back("mirror_plane 0.850650807 0 -0.525731114");
        break;

        case pg::I5H:
        std::cerr << "ERROR: Symmetry pg::I5H not implemented" << std::endl;
        exit(0);

        default:
        std::cerr << "ERROR: Symmetry " << symmetry  << "is not known" << std::endl;
        exit(0);

    }

    // #define DEBUG5
    #ifdef DEBUG5
    for (int n = 0; n < fileContent.size(); n++)
        std::cerr << fileContent[n] << std::endl;
    std::cerr << "fileContent.size()" << fileContent.size() << std::endl;
    #endif
    #undef DEBUG5
}

void SymList::writeDefinition(std::ostream &outstream, FileName fn_sym) {
    read_sym_file(fn_sym);
    Matrix2D<RFLOAT> L(3, 3), R(3, 3);
    outstream << " ++++ Using symmetry group " << fn_sym << ", with the following " << SymsNo() + 1 << " transformation matrices:" << std::endl;
    R.initIdentity();
    outstream << " R(1)= " << R;
    for (int isym = 0; isym < SymsNo(); isym++) {
        get_matrices(isym, L, R);
        R.resize(3, 3);
        L.resize(3, 3);
        if (!L.isIdentity())
            outstream << " L(" << isym + 2 << ")= " << L;
        outstream << " R(" << isym + 2 << ")= " << R;
        RFLOAT alpha, beta, gamma;
        Euler_matrix2angles(R, alpha, beta, gamma);
        outstream << "     Euler angles: " << alpha << " " << beta << " " << gamma << std::endl;
    }

}

RFLOAT SymList::non_redundant_ewald_sphere(int pgGroup, int pgOrder) {

    switch (pgGroup) {

        case pg::CN:
        case pg::SN:
        return 4.0 * PI / pgOrder;

        case pg::CI:
        case pg::CS:
        return 4.0 * PI / 2.0;

        case pg::CNV:
        case pg::CNH:
        case pg::DN:
        return 4.0 * PI / pgOrder / 2;

        case pg::DNV:
        case pg::DNH:
        return 4.0 * PI / pgOrder / 4;

        case pg::T:
        return 4.0 * PI / 12;

        case pg::TD:
        case pg::TH:
        case pg::O:
        return 4.0 * PI / 24;

        case pg::OH:
        return 4.0 * PI / 48;

        case pg::I:
        case pg::I1:
        case pg::I2:
        case pg::I3:
        case pg::I4:
        case pg::I5:
        return 4.0 * PI / 60;

        case pg::IH:
        case pg::I1H:
        case pg::I2H:
        case pg::I3H:
        case pg::I4H:
        case pg::I5H:
        return 4.0 * PI / 120;

        default:
        std::cerr 
        << "ERROR: Symmetry group, order=" << pgGroup 
        << " " <<  pgOrder << "is not known" << std::endl;
        exit(0);

    }
}

void symmetriseMap(MultidimArray<RFLOAT> &img, FileName &fn_sym, bool do_wrap) {

    if (img.getDim() != 3)
        REPORT_ERROR("symmetriseMap ERROR: symmetriseMap can only be run on 3D maps!");

    img.setXmippOrigin();

    SymList SL;
    SL.read_sym_file(fn_sym);

    Matrix2D<RFLOAT> L(4, 4), R(4, 4); // A matrix from the list
    MultidimArray<RFLOAT> sum = img;
    MultidimArray<RFLOAT> aux;
    aux.resize(img);

    for (int isym = 0; isym < SL.SymsNo(); isym++) {
        SL.get_matrices(isym, L, R);
        applyGeometry(img, aux, R, IS_INV, do_wrap);
        sum += aux;
    }

    // Overwrite the input
    img = sum / (SL.SymsNo() + 1);

}
