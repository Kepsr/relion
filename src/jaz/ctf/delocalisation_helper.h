#ifndef DELOCALISATION_HELPER_H
#define DELOCALISATION_HELPER_H

#include <src/ctf.h>

namespace DelocalisationHelper {

    void maskOutsideBox(
        const CTF& ctf, ObservationModel *obsModel, int opticsGroup,
        double radius, double angpix, int s_orig,
        MultidimArray<RFLOAT>& fftwCtfImg,
        double offsetx, double offsety
    );

    Image<RFLOAT> plotDelocalisation(
        const CTF& ctf, ObservationModel *obsModel, int opticsGroup,
        Image<RFLOAT>& mask, double angpix
    );

};

#endif
