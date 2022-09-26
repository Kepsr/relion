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

#ifndef STACK_HELPER_H
#define STACK_HELPER_H

#include <vector>
#include "src/image.h"
#include "src/metadata_table.h"
#include "src/jaz/optimization/optimization.h"
#include "src/jaz/gravis/t2Matrix.h"
#include "src/jaz/obs_model.h"

namespace StackHelper {

    std::vector<MetaDataTable> splitByMicrographName(const MetaDataTable& mdt);

    MetaDataTable merge(const std::vector<MetaDataTable>& mdts);

    std::vector<MetaDataTable> splitByStack(const MetaDataTable& mdt);

    std::vector<Image<RFLOAT>> loadStack(
        const MetaDataTable &mdt, std::string path = "", int threads = 1);

    std::vector<Image<Complex>> loadStackFS(
        const MetaDataTable& mdt,
        std::string path = "",
        int threads = 1,
        bool centerParticle = false,
        ObservationModel* obs = 0);

    void saveStack(std::vector<Image<RFLOAT>>& stack, std::string fn);

    std::vector<std::vector<Image<RFLOAT>>> loadMovieStack(
        const MetaDataTable &mdt, const std::string &moviePath);

    // For movies in file
    std::vector<std::vector<Image<Complex>>> extractMovieStackFS(
        const MetaDataTable &mdt,
        Image<RFLOAT>* gainRef, MultidimArray<bool>* defectMask, std::string movieFn,
        double outPs, double coordsPs, double moviePs, double dataPs,
        int squareSize, int threads,
        bool loadData = true, int firstFrame = 0, int lastFrame = -1,
        RFLOAT hot = -1.0, bool verbose = false, bool saveMemory = false,
        const std::vector<std::vector<gravis::d2Vector>>* offsets_in = 0,
        std::vector<std::vector<gravis::d2Vector>>* offsets_out = 0);

    // For movies in memory
    std::vector<std::vector<Image<Complex>>> extractMovieStackFS(
        const MetaDataTable &mdt, std::vector<MultidimArray<float>> &mgStack,
        double outPs, double coordsPs, double moviePs, double dataPs,
        int squareSize, int threads,
        bool loadData = true,
        bool verbose = false,
        const std::vector<std::vector<gravis::d2Vector>>* offsets_in = 0,
        std::vector<std::vector<gravis::d2Vector>>* offsets_out = 0);

    // Map FT over stack
    std::vector<Image<Complex>> FourierTransform(
        const std::vector<Image<RFLOAT>> &stack);

    // Map IFT over stack
    std::vector<Image<RFLOAT>> inverseFourierTransform(
        const std::vector<Image<Complex>> &stack);

    Image<RFLOAT> toSingleImage(const std::vector<Image<RFLOAT>> &stack);

    void varianceNormalize(
        std::vector<Image<Complex>>& movie,
        bool circleCropped = false);

    std::vector<double> powerSpectrum(
        const std::vector<std::vector<Image<Complex>>>& stack);

    std::vector<double> varSpectrum(
        const std::vector<std::vector<Image<Complex>>>& stack);

    std::vector<double> powerSpectrum(
        const std::vector<std::vector<Image<Complex>>>& obs,
        const std::vector<Image<Complex>>& signal);

};

#endif
