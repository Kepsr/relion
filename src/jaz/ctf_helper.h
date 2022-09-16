/***************************************************************************
 *
 * Author: "Jasenko Zivanov"
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

#ifndef CTF_HELPER_H
#define CTF_HELPER_H

#include <vector>
#include <src/ctf.h>
#include <src/image.h>
#include <src/metadata_table.h>
#include <src/jaz/optimization/optimization.h>
#include <src/jaz/volume.h>
#include <src/jaz/gravis/t2Matrix.h>
#include "src/jaz/obs_model.h"

namespace CtfHelper {

    std::vector<CTF> loadCtffind4(
        std::string path, int imageCount,
        double voltage = 300.0, double Cs = 2.2,
        double Q0 = 0.1, double Bfac = 0.0,
        double scale = 1.0
    );

    CTF setFromFile(
        std::stringstream& line,
        double voltage, double Cs, double Q0, double Bfac, double scale
    );

    /** Set all values explicitly in 3.1 */
    void setValuesByGroup(
        CTF &ctf, ObservationModel *obs, int opticsGroup,
        RFLOAT defU, RFLOAT defV, RFLOAT defAng,
        RFLOAT Bfac = 0.0, RFLOAT scale = 1.0, RFLOAT phase_shift = 0.0
    );

    // Read CTF parameters from particle table partMdt and optics table opticsMdt.
    void readByGroup(
        CTF &ctf, const MetaDataTable &partMdt, ObservationModel *obs, long int particle = -1
    );

    CTF makeCTF(const MetaDataTable &mdt, long int objectID = -1);
    CTF makeCTF(const MetaDataTable &MD1, const MetaDataTable &MD2, long int objectID = -1);
    CTF makeCTF(const MetaDataTable &partMdt, ObservationModel *obs, long int particle = -1);

    CTF makeCTF(
        ObservationModel *obs, int opticsGroup,
        RFLOAT defU, RFLOAT defV, RFLOAT defAng,
        RFLOAT Bfac = 0.0, RFLOAT scale = 1.0, RFLOAT phase_shift = 0.0
    );

    // Read from a MetaDataTable
    void read(CTF &ctf, const MetaDataTable &MD, long int objectID = -1);

    /** Read CTF parameters from MetaDataTables MD1 and MD2 (deprecated).
     * If a parameter is not found in MD1 it is tried to be read from MD2.
     * If it is also not found in the second then a default value is used.
     * This is useful if micrograph-specific parameters are stored in a separate MD from the image-specific parameters.
     */
    void read(
        CTF &ctf, const MetaDataTable &MD1, const MetaDataTable &MD2,
        long int objectID = -1
    );

    RFLOAT readValue(
        EMDL::EMDLabel label, RFLOAT defaultVal,
        long int particle, int opticsGroup,
        const MetaDataTable &partMdt, const ObservationModel *obs
    );

    // Write to a MetaDataTable
    void write(const CTF &ctf, MetaDataTable &MD);

    // Write to an output stream
    void write(const CTF &ctf, std::ostream &out);

    // Generate (Fourier-space, i.e. FFTW format) image with all CTF values.
    // The dimensions of the result array should have been set correctly already
    MultidimArray<RFLOAT> getFftwImage(
        const CTF &ctf,
        long int Xdim, long int Ydim, int orixdim, int oriydim, RFLOAT angpix,
        ObservationModel *obsModel, int opticsGroup,
        bool do_abs = false, bool do_only_flip_phases = false,
        bool do_intact_until_first_peak = false, bool do_damping = true,
        bool do_ctf_padding = false, bool do_intact_after_first_peak = false
    );

    MultidimArray<RFLOAT> getFftwImage_padded(
        const CTF &ctf,
        long int Xdim, long int Ydim, int orixdim, int oriydim, RFLOAT angpix,
        ObservationModel *obsModel, int opticsGroup,
        bool do_abs = false, bool do_only_flip_phases = false,
        bool do_intact_until_first_peak = false, bool do_damping = true,
        bool do_intact_after_first_peak = false
    );

    // Get a complex image with the CTFP/Q values, where the angle is in degrees between the Y-axis and the CTFP/Q sector line
    MultidimArray<Complex> getCTFPImage(
        const CTF &ctf,
        long int Xdim, long int Ydim, int orixdim, int oriydim, RFLOAT angpix,
        ObservationModel *obsModel, int opticsGroup,
        bool is_positive, float angle
    );


    // Generate a centered image (with Hermitian symmetry)
    // The dimensions of the result array should have been set correctly already
    MultidimArray<RFLOAT> getCenteredImage(
        const CTF &ctf, long int Xdim, long int Ydim, RFLOAT angpix,
        ObservationModel *obsModel, int opticsGroup,
        bool do_abs = false, bool do_only_flip_phases = false,
        bool do_intact_until_first_peak = false, bool do_damping = true,
        bool do_intact_after_first_peak = false
    );


    // Generate a 1D profile along defocusAngle
    // The dimensions of the result array should have been set correctly already, i.e. at the image size!
    MultidimArray<RFLOAT> get1DProfile(
        const CTF &ctf, long int Xdim, long int Ydim, RFLOAT angle, RFLOAT angpix,
        ObservationModel *obsModel, int opticsGroup,
        bool do_abs = false, bool do_only_flip_phases = false,
        bool do_intact_until_first_peak = false, bool do_damping = true,
        bool do_intact_after_first_peak = false
    );

    // Calculate weight W for Ewald-sphere curvature correction: apply this to the result from getFftwImage
    void applyWeightEwaldSphereCurvature(
        CTF &ctf,
        MultidimArray<RFLOAT> &result, int orixdim, int oriydim, RFLOAT angpix,
        ObservationModel *obsModel, int opticsGroup,
        RFLOAT particle_diameter
    );

    void applyWeightEwaldSphereCurvature_new(
        CTF &ctf,
        MultidimArray<RFLOAT> &result, int orixdim, int oriydim, RFLOAT angpix,
        ObservationModel *obsModel, int opticsGroup,
        RFLOAT particle_diameter
    );

    // Calculate weight W for Ewald-sphere curvature correction: apply this to the result from getFftwImage
    void applyWeightEwaldSphereCurvature_noAniso(
        CTF &ctf,
        MultidimArray<RFLOAT> &result, int orixdim, int oriydim, RFLOAT angpix,
        ObservationModel *obsModel, int opticsGroup,
        RFLOAT particle_diameter
    );

};

#endif
