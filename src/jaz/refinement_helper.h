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

#ifndef REFINEMENT_HELPER_H
#define REFINEMENT_HELPER_H

#include <vector>
#include <src/ctf.h>
#include <src/image.h>
#include <src/metadata_table.h>
#include <src/jaz/optimization/optimization.h>
#include <src/jaz/volume.h>
#include <src/jaz/gravis/t2Matrix.h>
#include "src/jaz/obs_model.h"

namespace RefinementHelper {

    Image<RFLOAT> drawFSC(const MetaDataTable* mdt, std::vector<double> &dest1D,
                            double thresh = 0.143);

    Image<RFLOAT> computeSNR(const MetaDataTable* mdt, double eps = 1e-15);
    Image<RFLOAT> computeSigInvSq(const MetaDataTable* mdt, const std::vector<double> &signalPow,
                                    double eps = 1e-15);

    Image<RFLOAT> correlation(const Image<Complex> &prediction,
                                            const Image<Complex> &observation);

    Image<RFLOAT> correlation(const std::vector<Image<Complex> > &prediction,
                                            const std::vector<Image<Complex> > &observation);

    void addToQR(
            const Image<Complex> &prediction, const Image<Complex> &observation,
            Image<Complex> &q, Image<RFLOAT> &r);

    void addToPQR(
            const Image<Complex> &prediction,
            const Image<Complex> &observation,
            Image<RFLOAT> &p, Image<Complex> &q, Image<RFLOAT> &r);

    double squaredDiff(
            const Image<Complex> &prediction, const Image<Complex> &observation,
            CTF &ctf, ObservationModel *obsModel, int opticsGroup,
            RFLOAT angpix, const Image<RFLOAT> &weight);

    double squaredDiff(
            const std::vector<Image<Complex>> &predictions,
            const std::vector<Image<Complex>> &observations,
            CTF &ctf, ObservationModel *obsModel, int opticsGroup,
            RFLOAT angpix, const Image<RFLOAT> &weight);

};

#endif
