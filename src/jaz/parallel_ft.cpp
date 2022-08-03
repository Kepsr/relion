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

#include "src/macros.h"
#include "src/jaz/parallel_ft.h"
#include "src/fftw.h"
#include "src/args.h"
#include <string.h>
#include <math.h>

static pthread_mutex_t fftw_plan_mutex_par = PTHREAD_MUTEX_INITIALIZER;

// #define DEBUG_PLANS

// Constructors and destructors --------------------------------------------
ParFourierTransformer::ParFourierTransformer() {}

ParFourierTransformer::~ParFourierTransformer() {
    clear();
    #ifdef DEBUG_PLANS
    std::cerr << "CLEARED this= " << this << std::endl;
    #endif
}

ParFourierTransformer::ParFourierTransformer(const ParFourierTransformer& op) {
    plans_are_set = false;
    clear();     // Clear current object
    *this = op;  // Make an extact copy of op
}

// Initialization ----------------------------------------------------------
void ParFourierTransformer::setReal(MultidimArray<RFLOAT> &input) {
    bool recomputePlan =
        !fReal ||
        // dataPtr != input.data ||
        !fReal->sameShape(input);

    fFourier.reshape(Xsize(input) / 2 + 1, Ysize(input), Zsize(input));
    fReal = &input;

    if (recomputePlan) {
        const int rank = get_array_rank(input);
        const int *const n = new_n(input, rank);

        // Destroy any forward or backward plans
        destroyPlans();

        // Make new plans
        plans_are_set = true;

        {
        pthread_lock_guard guard (&fftw_plan_mutex_par);
        fPlanForward = FFTW_PLAN_DFT_R2C(
            rank, n, fReal->data, (FFTW_COMPLEX*) fFourier.data, FFTW_ESTIMATE
        );
        fPlanBackward = FFTW_PLAN_DFT_C2R(
            rank, n, (FFTW_COMPLEX*) fFourier.data, fReal->data, FFTW_ESTIMATE
        );
        }

        if (!fPlanForward || !fPlanBackward)
            REPORT_ERROR("FFTW plans cannot be created");

        #ifdef DEBUG_PLANS
        std::cerr << " SETREAL fPlanForward= " << fPlanForward << " fPlanBackward= " << fPlanBackward  <<" this= " << this << std::endl;
        #endif

        delete[] n;
        dataPtr = fReal->data;
    }
}

void ParFourierTransformer::setReal(MultidimArray<Complex> &input) {
    bool recomputePlan =
        !fComplex ||
        // complexDataPtr != input.data ||
        !fComplex->sameShape(input);

    fFourier.resize(input);
    fComplex = &input;

    if (recomputePlan) {
        const int rank = get_array_rank(input);
        const int *const n = new_n(input, rank);

        // Destroy any forward or backward plans
        destroyPlans();

        plans_are_set = true;

        {
        pthread_lock_guard guard (&fftw_plan_mutex_par);
        fPlanForward = FFTW_PLAN_DFT(
            rank, n, (FFTW_COMPLEX*) fComplex->data,
            (FFTW_COMPLEX*) fFourier.data, FFTW_FORWARD, FFTW_ESTIMATE
        );
        fPlanBackward = FFTW_PLAN_DFT(
            rank, n, (FFTW_COMPLEX*) fFourier.data,
            (FFTW_COMPLEX*) fComplex->data, FFTW_BACKWARD, FFTW_ESTIMATE
        );
        }

        if (!fPlanForward || !fPlanBackward)
            REPORT_ERROR("FFTW plans cannot be created");

        delete[] n;
        complexDataPtr = fComplex->data;
    }
}

void ParFourierTransformer::setFourier(const MultidimArray<Complex> &inputFourier) {
    memcpy(
        fFourier.data, inputFourier.data,
        inputFourier.size() * sizeof(Complex)
    );
}

// Transform ---------------------------------------------------------------
void ParFourierTransformer::Transform(int sign) {
    if (sign == FFTW_FORWARD) {
        FFTW_EXECUTE_DFT_R2C(
            fPlanForward, fReal->data,
            (FFTW_COMPLEX*) fFourier.data
        );
        // Normalisation of the transform
        unsigned long int size = 0;
        if (fReal) {
            size = fReal->size();
        } else if (fComplex) {
            size = fComplex->size();
        } else {
            REPORT_ERROR("No complex nor real data defined");
        }

        for (auto &x : fFourier) { x /= size; }
    } else if (sign == FFTW_BACKWARD) {
        FFTW_EXECUTE_DFT_C2R(
            fPlanBackward, (FFTW_COMPLEX*)
            fFourier.data, fReal->data
        );
    }
}

template <typename T>
void ParFourierTransformer::getCompleteFourier(MultidimArray<T> &V) const {
    V.reshape(*fReal);
    switch (get_array_rank(*fReal)) {

        case 1:
        for (long int i = 0; i < Xsize(V); i++) {
            if (i < Xsize(fFourier)) {
                direct::elem(V, i) = direct::elem(fFourier, i);
            } else {
                direct::elem(V, i) = conj(direct::elem(
                    fFourier,
                    Xsize(*fReal) - i
                ));
            }
        }
        break;

        case 2:
        for (long int j = 0; j < Ysize(V); j++)
        for (long int i = 0; i < Xsize(V); i++) {
            if (j < Xsize(fFourier)) {
                direct::elem(V, i, j) = direct::elem(fFourier, i, j);
            } else {
                direct::elem(V, i, j) = conj(direct::elem(
                    fFourier,
                     Xsize(*fReal) - j,
                    (Ysize(*fReal) - i) % Ysize(*fReal)
                ));
            }
        }
        break;

        case 3:
            for (long int k = 0; k < Zsize(V); k++)
            for (long int j = 0; j < Ysize(V); j++)
            for (long int i = 0; i < Xsize(V); i++) {
                if (j < Xsize(fFourier)) {
                    direct::elem(V, i, j, k) = direct::elem(fFourier, i, j, k);
                } else {
                    direct::elem(V, i, j, k) = conj(direct::elem(
                        fFourier,
                         Xsize(*fReal) - j,
                        (Ysize(*fReal) - i) % Ysize(*fReal),
                        (Zsize(*fReal) - k) % Zsize(*fReal)
                    ));
                }
            }
            break;

    }
}

template <typename T>
void ParFourierTransformer::setFromCompleteFourier(const MultidimArray<T>& V) {
    switch (get_array_rank(*fReal)) {

        case 1:
        for (long int i = 0; i < Xsize(fFourier); i++) {
            direct::elem(fFourier, i) = direct::elem(V, i);
        }
        break;

        case 2:
        for (long int j = 0; j < Ysize(fFourier); j++)
        for (long int i = 0; i < Xsize(fFourier); i++) {
            direct::elem(fFourier, i, j) = direct::elem(V, i, j);
        }
        break;

        case 3:
        for (long int k = 0; k < Zsize(fFourier); k++)
        for (long int j = 0; j < Ysize(fFourier); j++)
        for (long int i = 0; i < Xsize(fFourier); i++)
            direct::elem(fFourier, i, j, k) = direct::elem(V, i, j, k);
        break;

    }
}

#undef FFTW_COMPLEX
#undef FFTW_PLAN_DFT
#undef FFTW_PLAN_DFT_R2C
#undef FFTW_PLAN_DFT_C2R
#undef FFTW_EXECUTE_DFT_R2C
#undef FFTW_EXECUTE_DFT_C2R
#undef FFTW_CLEANUP
#undef FFTW_DESTROY_PLAN
