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

#ifdef RELION_SINGLE_PRECISION
typedef fftwf_complex FFTW_COMPLEX;
#define FFTW_PLAN_DFT        fftwf_plan_dft
#define FFTW_PLAN_DFT_R2C    fftwf_plan_dft_r2c
#define FFTW_PLAN_DFT_C2R    fftwf_plan_dft_c2r
#define FFTW_EXECUTE_DFT_R2C fftwf_execute_dft_r2c
#define FFTW_EXECUTE_DFT_C2R fftwf_execute_dft_c2r
#define FFTW_CLEANUP         fftwf_cleanup
#define FFTW_DESTROY_PLAN    fftwf_destroy_plan
#else
typedef fftw_complex FFTW_COMPLEX;
#define FFTW_PLAN_DFT        fftw_plan_dft
#define FFTW_PLAN_DFT_R2C    fftw_plan_dft_r2c
#define FFTW_PLAN_DFT_C2R    fftw_plan_dft_c2r
#define FFTW_EXECUTE_DFT_R2C fftw_execute_dft_r2c
#define FFTW_EXECUTE_DFT_C2R fftw_execute_dft_c2r
#define FFTW_CLEANUP         fftw_cleanup
#define FFTW_DESTROY_PLAN    fftw_destroy_plan
#endif

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
    fReal          = NULL;
    fComplex       = NULL;
    fPlanForward   = NULL;
    fPlanBackward  = NULL;
    dataPtr        = NULL;
    complexDataPtr = NULL;
}

void ParFourierTransformer::clear() {
    fFourier.clear();
    // Clean-up all other FFTW-allocated memory
    destroyPlans();
    // Initialise all pointers to NULL
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
    pthread_mutex_lock(&fftw_plan_mutex_par);

    if (plans_are_set) {
        FFTW_DESTROY_PLAN(fPlanForward);
        FFTW_DESTROY_PLAN(fPlanBackward);
        plans_are_set = false;
    }

    pthread_mutex_unlock(&fftw_plan_mutex_par);
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
        // dataPtr != MULTIDIM_ARRAY(input) || 
        !fReal->sameShape(input);

    fFourier.reshape(Zsize(input), Ysize(input), Xsize(input) / 2 + 1);
    fReal = &input;

    if (recomputePlan) {
        int ndim = 3;
        if (Zsize(input) == 1) {
            ndim = 2;
            if (Ysize(input) == 1)
                ndim = 1;
        }
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

        pthread_mutex_lock(&fftw_plan_mutex_par);
        fPlanForward = FFTW_PLAN_DFT_R2C(
            ndim, N, MULTIDIM_ARRAY(*fReal),
            (FFTW_COMPLEX*) MULTIDIM_ARRAY(fFourier), FFTW_ESTIMATE
        );
        fPlanBackward = FFTW_PLAN_DFT_C2R(
            ndim, N, (FFTW_COMPLEX*) MULTIDIM_ARRAY(fFourier),
            MULTIDIM_ARRAY(*fReal), FFTW_ESTIMATE
        );
        pthread_mutex_unlock(&fftw_plan_mutex_par);

        if (!fPlanForward || !fPlanBackward)
            REPORT_ERROR("FFTW plans cannot be created");

        #ifdef DEBUG_PLANS
        std::cerr << " SETREAL fPlanForward= " << fPlanForward << " fPlanBackward= " << fPlanBackward  <<" this= " << this << std::endl;
        #endif

        delete [] N;
        dataPtr = MULTIDIM_ARRAY(*fReal);
    }
}

void ParFourierTransformer::setReal(MultidimArray<Complex> &input) {
    bool recomputePlan = 
        !fComplex || 
        // complexDataPtr != MULTIDIM_ARRAY(input) ||
        !fComplex->sameShape(input);

    fFourier.resize(input);
    fComplex = &input;

    if (recomputePlan) {
        int ndim = 3;
        if (Zsize(input) == 1) {
            ndim = 2;
            if (Ysize(input) == 1)
                ndim = 1;
        }
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

        pthread_mutex_lock(&fftw_plan_mutex_par);
        fPlanForward = FFTW_PLAN_DFT(
            ndim, N, (FFTW_COMPLEX*) MULTIDIM_ARRAY(*fComplex),
            (FFTW_COMPLEX*) MULTIDIM_ARRAY(fFourier), FFTW_FORWARD, FFTW_ESTIMATE
        );
        fPlanBackward = FFTW_PLAN_DFT(
            ndim, N, (FFTW_COMPLEX*) MULTIDIM_ARRAY(fFourier),
            (FFTW_COMPLEX*) MULTIDIM_ARRAY(*fComplex), FFTW_BACKWARD, FFTW_ESTIMATE
        );
        pthread_mutex_unlock(&fftw_plan_mutex_par);

        if (!fPlanForward || !fPlanBackward)
            REPORT_ERROR("FFTW plans cannot be created");

        delete [] N;
        complexDataPtr = MULTIDIM_ARRAY(*fComplex);
    }
}

void ParFourierTransformer::setFourier(const MultidimArray<Complex> &inputFourier) {
    memcpy(
        MULTIDIM_ARRAY(fFourier), MULTIDIM_ARRAY(inputFourier),
        MULTIDIM_SIZE(inputFourier) * 2 * sizeof(RFLOAT)
    );
}

// Transform ---------------------------------------------------------------
void ParFourierTransformer::Transform(int sign) {
    if (sign == FFTW_FORWARD) {
        FFTW_EXECUTE_DFT_R2C(
            fPlanForward, MULTIDIM_ARRAY(*fReal),
            (FFTW_COMPLEX*) MULTIDIM_ARRAY(fFourier)
        );
        // Normalisation of the transform
        unsigned long int size = 0;
        if (fReal) {
            size = MULTIDIM_SIZE(*fReal);
        } else if (fComplex) {
            size = MULTIDIM_SIZE(*fComplex);
        } else {
            REPORT_ERROR("No complex nor real data defined");
        }

        FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(fFourier) {
            DIRECT_MULTIDIM_ELEM(fFourier,n) /= size;
        }
    } else if (sign == FFTW_BACKWARD) {
        FFTW_EXECUTE_DFT_C2R(
            fPlanBackward, (FFTW_COMPLEX*) 
            MULTIDIM_ARRAY(fFourier), MULTIDIM_ARRAY(*fReal)
        );
    }
}

void ParFourierTransformer::FourierTransform() {
    Transform(FFTW_FORWARD);
}

void ParFourierTransformer::inverseFourierTransform() {
    Transform(FFTW_BACKWARD);
}

// Enforce Hermitian symmetry ---------------------------------------------
void ParFourierTransformer::enforceHermitianSymmetry() {
    int ndim = 3;
    if (Zsize(*fReal) == 1) {
        ndim = 2;
        if (Ysize(*fReal) == 1)
            ndim = 1;
    }
    long int yHalf = Ysize(*fReal) / 2;
    if (Ysize(*fReal) % 2 == 0) { yHalf--; }
    long int zHalf = Zsize(*fReal) / 2;
    if (Zsize(*fReal) % 2 == 0) { zHalf--; }

    switch (ndim) {

        case 2:
        for (long int i = 1; i <= yHalf; i++) {
            long int isym = wrap(-i, 0, Ysize(*fReal) - 1);
            Complex mean = 0.5 * (DIRECT_A2D_ELEM(fFourier, i, 0) + conj(DIRECT_A2D_ELEM(fFourier, isym, 0)));
            DIRECT_A2D_ELEM(fFourier, i,    0) = mean;
            DIRECT_A2D_ELEM(fFourier, isym, 0) = conj(mean);
        }
        break;

        case 3:
        for (long int k = 0; k < Zsize(*fReal); k++) {
            long int ksym = wrap(-k, 0, Zsize(*fReal) - 1);
            for (long int i = 1; i <= yHalf; i++) {
                long int isym = wrap(-i, 0, Ysize(*fReal) - 1);
                Complex mean = 0.5 * (DIRECT_A3D_ELEM(fFourier,k,i,0) + conj(DIRECT_A3D_ELEM(fFourier, ksym, isym, 0)));
                DIRECT_A3D_ELEM(fFourier, k,    i,    0) = mean;
                DIRECT_A3D_ELEM(fFourier, ksym, isym, 0) = conj(mean);
            }
        }
        for (long int k = 1; k <= zHalf; k++) {
            long int ksym = wrap(-k, 0, Zsize(*fReal) - 1);
            Complex mean = 0.5 * (DIRECT_A3D_ELEM(fFourier, k, 0, 0) + conj(DIRECT_A3D_ELEM(fFourier, ksym, 0, 0)));
            DIRECT_A3D_ELEM(fFourier, k,    0, 0) = mean;
            DIRECT_A3D_ELEM(fFourier, ksym, 0, 0) = conj(mean);
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
