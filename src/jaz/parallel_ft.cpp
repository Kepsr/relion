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

struct thread_guard {

    thread_guard() {
        pthread_mutex_lock(&fftw_plan_mutex_par);
    }

    ~thread_guard() {
        pthread_mutex_unlock(&fftw_plan_mutex_par);
    }

};

// #define DEBUG_PLANS

static inline int get_n(const MultidimArray<RFLOAT> &arr) {
    return Zsize(arr) == 1 ? Ysize(arr) == 1 ? 1 : 2 : 3;
}

// Constructors and destructors --------------------------------------------

ParFourierTransformer::ParFourierTransformer(): plans_are_set(false) {
    init();
    #ifdef DEBUG_PLANS
    std::cerr << "INIT this= " << this << std::endl;
    #endif
}

ParFourierTransformer::~ParFourierTransformer() {
    clear();
    #ifdef DEBUG_PLANS
    std::cerr << "CLEARED this= " << this << std::endl;
    #endif
}

ParFourierTransformer::ParFourierTransformer(const ParFourierTransformer& op): plans_are_set(false) {
    // Clear current object
    clear();
    // New object is an extact copy of op
    *this = op;
}

void ParFourierTransformer::init() {
    fReal          = nullptr;
    fComplex       = nullptr;
    fPlanForward   = nullptr;
    fPlanBackward  = nullptr;
    dataPtr        = nullptr;
    complexDataPtr = nullptr;
}

void ParFourierTransformer::clear() {
    fFourier.clear();
    // Clean-up all other FFTW-allocated memory
    destroyPlans();
    // Initialise all pointers to nullptr
    init();
}

void ParFourierTransformer::cleanup() {
    // First clear object and destroy plans
    clear();
    // Then clean up all the junk fftw keeps lying around
    // SOMEHOW THE FOLLOWING IS NOT ALLOWED WHEN USING MULTPLE TRANSFORMER OBJECTS....
    FFTW_CLEANUP();

    #ifdef DEBUG_PLANS
    std::cerr << "CLEANED-UP this= " << this << std::endl;
    #endif
}

void ParFourierTransformer::destroyPlans() {
    // Anything to do with plans has to be protected for threads!
    thread_guard g;
    if (plans_are_set) {
        FFTW_DESTROY_PLAN(fPlanForward);
        FFTW_DESTROY_PLAN(fPlanBackward);
        plans_are_set = false;
    }
}

// Initialization ----------------------------------------------------------
const MultidimArray<RFLOAT> &ParFourierTransformer::getReal() const {
    return *fReal;
}

const MultidimArray<Complex > &ParFourierTransformer::getComplex() const {
    return *fComplex;
}


void ParFourierTransformer::setReal(MultidimArray<RFLOAT> &input) {
    bool recomputePlan =
        !fReal ||
        // dataPtr != input.data ||
        !fReal->sameShape(input);

    fFourier.reshape(Xsize(input) / 2 + 1, Ysize(input), Zsize(input));
    fReal = &input;

    if (recomputePlan) {
        int ndim = get_n(input);
        int *N = new int[ndim];
        switch (ndim) {

            case 1:
            N[0] = Xsize(input);
            break;

            case 2:
            N[0] = Ysize(input);
            N[1] = Xsize(input);
            break;

            case 3:
            N[0] = Zsize(input);
            N[1] = Ysize(input);
            N[2] = Xsize(input);
            break;

        }

        // Destroy any forward or backward plans
        destroyPlans();

        // Make new plans
        plans_are_set = true;

        {
        thread_guard g;
        fPlanForward = FFTW_PLAN_DFT_R2C(
            ndim, N, fReal->data,
            (FFTW_COMPLEX*) fFourier.data, FFTW_ESTIMATE
        );
        fPlanBackward = FFTW_PLAN_DFT_C2R(
            ndim, N, (FFTW_COMPLEX*) fFourier.data,
            fReal->data, FFTW_ESTIMATE
        );
        }

        if (!fPlanForward || !fPlanBackward)
            REPORT_ERROR("FFTW plans cannot be created");

        #ifdef DEBUG_PLANS
        std::cerr << " SETREAL fPlanForward= " << fPlanForward << " fPlanBackward= " << fPlanBackward  <<" this= " << this << std::endl;
        #endif

        delete[] N;
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
        int ndim = get_n(input);
        int *N = new int[ndim];
        switch (ndim) {

            case 1:
            N[0] = Xsize(input);
            break;

            case 2:
            N[0] = Ysize(input);
            N[1] = Xsize(input);
            break;

            case 3:
            N[0] = Zsize(input);
            N[1] = Ysize(input);
            N[2] = Xsize(input);
            break;

        }

        // Destroy any forward or backward plans
        destroyPlans();

        plans_are_set = true;

        {
        thread_guard g;
        fPlanForward = FFTW_PLAN_DFT(
            ndim, N, (FFTW_COMPLEX*) fComplex->data,
            (FFTW_COMPLEX*) fFourier.data, FFTW_FORWARD, FFTW_ESTIMATE
        );
        fPlanBackward = FFTW_PLAN_DFT(
            ndim, N, (FFTW_COMPLEX*) fFourier.data,
            (FFTW_COMPLEX*) fComplex->data, FFTW_BACKWARD, FFTW_ESTIMATE
        );
        }

        if (!fPlanForward || !fPlanBackward)
            REPORT_ERROR("FFTW plans cannot be created");

        delete[] N;
        complexDataPtr = fComplex->data;
    }
}

void ParFourierTransformer::setFourier(const MultidimArray<Complex> &inputFourier) {
    memcpy(
        fFourier.data, inputFourier.data,
        inputFourier.size() * 2 * sizeof(RFLOAT)
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

void ParFourierTransformer::FourierTransform() {
    Transform(FFTW_FORWARD);
}

void ParFourierTransformer::inverseFourierTransform() {
    Transform(FFTW_BACKWARD);
}

template <typename T>
void ParFourierTransformer::getCompleteFourier(MultidimArray<T> &V) const {
    V.reshape(*fReal);
    switch (get_n(*fReal)) {

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
    switch (get_n(*fReal)) {

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


// Enforce Hermitian symmetry ---------------------------------------------
void ParFourierTransformer::enforceHermitianSymmetry() {
    long int yHalf = Ysize(*fReal) / 2;
    if (Ysize(*fReal) % 2 == 0) { yHalf--; }
    long int zHalf = Zsize(*fReal) / 2;
    if (Zsize(*fReal) % 2 == 0) { zHalf--; }

    switch (get_n(*fReal)) {

        case 2:
        for (long int j = 1; j <= yHalf; j++) {
            long int jsym = wrap(-j, 0, Ysize(*fReal) - 1);
            Complex mean = 0.5 * (direct::elem(fFourier, 0, j) + conj(direct::elem(fFourier, 0, jsym)));
            direct::elem(fFourier, 0, j) = mean;
            direct::elem(fFourier, 0, jsym) = conj(mean);
        }
        break;

        case 3:
        for (long int k = 0; k < Zsize(*fReal); k++) {
            long int ksym = wrap(-k, 0, Zsize(*fReal) - 1);
            for (long int j = 1; j <= yHalf; j++) {
                long int jsym = wrap(-j, 0, Ysize(*fReal) - 1);
                Complex mean = 0.5 * (direct::elem(fFourier, 0, j, k) + conj(direct::elem(fFourier, 0, jsym, ksym)));
                direct::elem(fFourier, 0, j,    k) = mean;
                direct::elem(fFourier, 0, jsym, ksym) = conj(mean);
            }
        }
        for (long int k = 1; k <= zHalf; k++) {
            long int ksym = wrap(-k, 0, Zsize(*fReal) - 1);
            Complex mean = 0.5 * (direct::elem(fFourier, 0, 0, k) + conj(direct::elem(fFourier, 0, 0, ksym)));
            direct::elem(fFourier, 0, 0, k) = mean;
            direct::elem(fFourier, 0, 0, ksym) = conj(mean);
        }
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
