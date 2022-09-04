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

#include <unistd.h>
#include <string.h>

#include <src/image.h>
#include <src/backprojector.h>
#include <src/funcs.h>
#include <src/memory.h>
#include <src/euler.h>
#include <src/time.h>
#include "src/jaz/complex_io.h"



class ext_recons_parameters {

    public:

    FileName fn_star, fn_recons, fn_data_real, fn_data_imag, fn_weight, fn_out;
    MultidimArray<RFLOAT> tau2;
    RFLOAT tau2_fudge;
    float padding_factor;
    int ori_size, current_size, ref_dim;
    int verb;
    bool skip_gridding, do_map;

    void read(int argc, char **argv) {
        if (argc < 2) {
            REPORT_ERROR("  Usage: relion_external_reconstruct input.star");
        }
        const FileName fn_star = argv[1];
        if (fn_star.getExtension() != "star") {
            REPORT_ERROR(" ERROR: input argument does not have a .star extension.");
        }
        skip_gridding = checkParameter(argc, argv, "--skip_gridding");
        do_map = !checkParameter(argc, argv, "--no_map");
        fn_out =  checkParameter(argc, argv, "--o") ? getParameter(argc, argv, "--o") : "";

        MetaDataTable MDlist, MDtau;
        MDlist.read(fn_star, "external_reconstruct_general");
        const long int i = MDlist.size() - 1;
        fn_data_real   = MDlist.getValue<std::string>(EMDL::OPTIMISER_EXTERNAL_RECONS_DATA_REAL, i);
        fn_data_imag   = MDlist.getValue<std::string>(EMDL::OPTIMISER_EXTERNAL_RECONS_DATA_IMAG, i);
        fn_weight      = MDlist.getValue<std::string>(EMDL::OPTIMISER_EXTERNAL_RECONS_WEIGHT, i);
        fn_recons      = MDlist.getValue<std::string>(EMDL::OPTIMISER_EXTERNAL_RECONS_RESULT, i);
        tau2_fudge     = MDlist.getValue<RFLOAT>(EMDL::MLMODEL_TAU2_FUDGE_FACTOR, i);
        padding_factor = MDlist.getValue<float>(EMDL::MLMODEL_PADDING_FACTOR, i);
        ref_dim        = MDlist.getValue<int>(EMDL::MLMODEL_DIMENSIONALITY, i);
        ori_size       = MDlist.getValue<int>(EMDL::MLMODEL_ORIGINAL_SIZE, i);
        current_size   = MDlist.getValue<int>(EMDL::MLMODEL_CURRENT_SIZE, i);

        if (!fn_out.empty()) { fn_recons = fn_out; }

        MDtau.read(fn_star, "external_reconstruct_tau2");
        tau2.resize(MDtau.size());
        for (long int i : MDtau) {
            tau2[i] = MDtau.getValue<RFLOAT>(EMDL::MLMODEL_TAU2_REF, i);
        }

    }

    void reconstruct() {
        BackProjector BP(ori_size, ref_dim, "C1", TRILINEAR, padding_factor);
        BP.initZeros(current_size);

        if (skip_gridding) BP.skip_gridding = skip_gridding;

        Image<Complex> Idata;
        const std::string fn_ext = "." + fn_data_real.getExtension();
        const std::string fn_root = fn_data_real.beforeFirstOf("_real");
        ComplexIO::read(Idata, fn_root, fn_ext);
        const auto Iweight = Image<RFLOAT>::from_filename(fn_weight);

        // Could there be a 1-pixel different in size? use FOR_ALL_ELEMENTS_IN_FFTW_TRANSFORM to be safe
        const int r_max = current_size / 2;
        const int r_max2 = round(r_max * padding_factor) * round(r_max * padding_factor);
        FOR_ALL_ELEMENTS_IN_FFTW_TRANSFORM(Idata()) {
            if (euclidsq(ip, jp, kp) < r_max2) {
                BP.data  .elem(ip, jp, kp) = direct::elem(Idata(),   i, j, k);
                BP.weight.elem(ip, jp, kp) = direct::elem(Iweight(), i, j, k);
            }
        }

        Image<RFLOAT>(BP.reconstruct(10, do_map, tau2, tau2_fudge)).write(fn_recons);
    }
};


int main(int argc, char *argv[]) {
    ext_recons_parameters prm;

    try {
        prm.read(argc, argv);
        prm.reconstruct();
    } catch (RelionError XE) {
        std::cerr << XE;
        return RELION_EXIT_FAILURE;
    }

    return RELION_EXIT_SUCCESS;
}
