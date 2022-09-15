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
 * Authors:    Roberto Marabini					(roberto@cnb.csic.es)
 *			   Carlos Oscar S. Sorzano			(coss@cnb.csic.es)
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
 *	All comments concerning this program package may be sent to the
 *	e-mail address 'xmipp@cnb.csic.es'
 ***************************************************************************/

#include "src/macros.h"
#include "src/fftw.h"
#include "src/args.h"
#include <string.h>
#include <math.h>
#include "src/multidim_array_statistics.h"

using RealArray = MultidimArray<RFLOAT>;
using ComplexArray = MultidimArray<Complex>;

// Anything to do with plans has to be protected for threads
static pthread_mutex_t fftw_plan_mutex = PTHREAD_MUTEX_INITIALIZER;

// #define TIMING_FFTW
#ifdef TIMING_FFTW
#include "src/time.h"
Timer timer_fftw;
const int TIMING_FFTW_PLAN      = timer_fftw.setNew("fftw - plan");
const int TIMING_FFTW_EXECUTE   = timer_fftw.setNew("fftw - exec");
const int TIMING_FFTW_NORMALISE = timer_fftw.setNew("fftw - normalise");
const int TIMING_FFTW_COPY      = timer_fftw.setNew("fftw - copy");
#define ifdefTIMING_FFTW(statement) statement
#else
#define ifdefTIMING_FFTW(statement)
#endif

// #define DEBUG_PLANS

// Constructors and destructors --------------------------------------------
FourierTransformer::FourierTransformer(): plans_are_set(false) {
    init();
    #ifdef DEBUG_PLANS
    std::cerr << "INIT this= " << this << std::endl;
    #endif
}

FourierTransformer::~FourierTransformer() {
    clear();
    #ifdef DEBUG_PLANS
    std::cerr << "CLEARED this= " << this << std::endl;
    #endif
}

void FourierTransformer::init() {
    fReal          = nullptr;
    fComplex       = nullptr;
    fPlanForward   = nullptr;
    fPlanBackward  = nullptr;
    dataPtr        = nullptr;
    complexDataPtr = nullptr;
}

void FourierTransformer::clear() {
    fFourier.clear();
    // Clean-up all other FFTW-allocated memory
    destroyPlans();
    // Initialise all pointers to nullptr
    init();
}

void FourierTransformer::cleanup() {
    // First clear object and destroy plans
    clear();
    // Then clean up all the junk fftw keeps lying around
    // SOMEHOW THE FOLLOWING IS NOT ALLOWED WHEN USING MULTPLE TRANSFORMER OBJECTS....
    FFTW_CLEANUP();

    #ifdef DEBUG_PLANS
    std::cerr << "CLEANED-UP this= " << this << std::endl;
    #endif

}

void FourierTransformer::destroyPlans() {
    if (plans_are_set) {
        pthread_lock_guard guard (&fftw_plan_mutex);
        FFTW_DESTROY_PLAN(fPlanForward);
        FFTW_DESTROY_PLAN(fPlanBackward);
        plans_are_set = false;
    }
}

const MultidimArray<RFLOAT> &FourierTransformer::getReal() const { return *fReal; }

const MultidimArray<Complex> &FourierTransformer::getComplex() const { return *fComplex; }

// Initialization ----------------------------------------------------------
void FourierTransformer::setReal(const MultidimArray<RFLOAT> &input, bool force_new_plans) {

    const bool recomputePlan =
        force_new_plans ||
        !fReal ||
        dataPtr != input.data ||
        !fReal->sameShape(input) ||
        Xsize(fFourier) != Xsize(input) / 2 + 1 ||
        complexDataPtr != fFourier.data;

    if (!recomputePlan) return;

    fFourier.reshape(Xsize(input) / 2 + 1, Ysize(input), Zsize(input));
    fReal = &input;

    const int rank = get_array_rank(input);
    const int *const n = new_n(input, rank);

    // Destroy both forward and backward plans if they already exist
    destroyPlans();

    // Make new plans
    {
    ifdefTIMING_FFTW(TicToc tt (timer_fftw, TIMING_FFTW_PLAN);)
    pthread_lock_guard guard (&fftw_plan_mutex);
    fPlanForward  = FFTW_PLAN_DFT_R2C(rank, n, fReal->data,
        (FFTW_COMPLEX*) fFourier.data, FFTW_ESTIMATE);
    fPlanBackward = FFTW_PLAN_DFT_C2R(rank, n, (FFTW_COMPLEX*) fFourier.data,
                          fReal->data, FFTW_ESTIMATE);
    }

    delete[] n;

    if (!fPlanForward || !fPlanBackward)
        REPORT_ERROR("FFTW plans cannot be created");

    #ifdef DEBUG_PLANS
    std::cerr << " SETREAL fPlanForward= " << fPlanForward << " fPlanBackward= " << fPlanBackward  << " this= " << this << std::endl;
    #endif

    plans_are_set = true;
    dataPtr = fReal->data;
    complexDataPtr = fFourier.data;

}

void FourierTransformer::setReal(const MultidimArray<Complex> &input, bool force_new_plans) {

    const bool recomputePlan =
        force_new_plans ||
        !fComplex ||
        complexDataPtr != input.data ||
        !fComplex->sameShape(input);

    fFourier.resize(input);
    fComplex = &input;

    if (!recomputePlan) return;

    const int rank = get_array_rank(input);
    const int *const n = new_n(input, rank);

    // Destroy both forward and backward plans if they already exist
    destroyPlans();

    {
    ifdefTIMING_FFTW(TicToc tt (timer_fftw, TIMING_FFTW_PLAN);)
    pthread_lock_guard guard (&fftw_plan_mutex);
    fPlanForward  = FFTW_PLAN_DFT(rank, n, (FFTW_COMPLEX*) fComplex->data,
        (FFTW_COMPLEX*) fFourier.data,  FFTW_FORWARD,  FFTW_ESTIMATE);
    fPlanBackward = FFTW_PLAN_DFT(rank, n, (FFTW_COMPLEX*) fFourier.data,
        (FFTW_COMPLEX*) fComplex->data, FFTW_BACKWARD, FFTW_ESTIMATE);
    }

    delete[] n;

    if (!fPlanForward || !fPlanBackward)
        REPORT_ERROR("FFTW plans cannot be created");

    plans_are_set = true;
    complexDataPtr = fComplex->data;

}

void FourierTransformer::setFourier(const MultidimArray<Complex> &inputFourier) {
    ifdefTIMING_FFTW(TicToc tt (timer_fftw, TIMING_FFTW_COPY);)

    if (!fFourier.sameShape(inputFourier)) {
        std::cerr << " fFourier= "; fFourier.printShape(std::cerr);
        std::cerr << " inputFourier= "; inputFourier.printShape(std::cerr);
        REPORT_ERROR("BUG: incompatible shaped in setFourier part of FFTW transformer");
    }

    memcpy(fFourier.data, inputFourier.data, inputFourier.size() * sizeof(Complex));

}

template <typename T>
void FourierTransformer::setFromCompleteFourier(const MultidimArray<T> &V) {
    switch (get_array_rank(V)) {
        case 1:
        for (long int i = 0; i < Xsize(fFourier); i++)
            direct::elem(fFourier, i) = direct::elem(V, i);
        break;
        case 2:
        for (long int j = 0; j < Ysize(fFourier); j++)
        for (long int i = 0; i < Xsize(fFourier); i++)
            direct::elem(fFourier, i, j) = direct::elem(V, i, j);
        break;
        case 3:
        for (long int k = 0; k < Zsize(fFourier); k++)
        for (long int j = 0; j < Ysize(fFourier); j++)
        for (long int i = 0; i < Xsize(fFourier); i++)
            direct::elem(fFourier, i, j, k) = direct::elem(V, i, j, k);
        break;
    }
}

inline unsigned long int getsize(const FourierTransformer &t) {
    if (t.fReal)    return t.fReal->size();
    if (t.fComplex) return t.fComplex->size();
    REPORT_ERROR("No data defined");
}

// Transform ---------------------------------------------------------------
void FourierTransformer::Transform(int sign) {
    switch (sign) {

        case FFTW_FORWARD:
        {
        ifdefTIMING_FFTW(TicToc tt (timer_fftw, TIMING_FFTW_EXECUTE);)
        FFTW_EXECUTE_DFT_R2C(fPlanForward, fReal->data, (FFTW_COMPLEX*) fFourier.data);
        }

        // Normalise the transform
        {
        ifdefTIMING_FFTW(TicToc tt (timer_fftw, TIMING_FFTW_NORMALISE);)
        const RFLOAT n = getsize(*this);
        for (auto &x : fFourier) { x /= n; }
        /// TODO: fFourier /= (RFLOAT) getsize(*this);
        }
        return;

        case FFTW_BACKWARD:
        {
        ifdefTIMING_FFTW(TicToc tt (timer_fftw, TIMING_FFTW_EXECUTE);)
        FFTW_EXECUTE_DFT_C2R(
            fPlanBackward, (FFTW_COMPLEX*) fFourier.data, fReal->data
        );
        }
        return;

    }
}

void FourierTransformer::FourierTransform() {
    Transform(FFTW_FORWARD);
}

void FourierTransformer::inverseFourierTransform() {
    Transform(FFTW_BACKWARD);
}

// Inforce Hermitian symmetry ---------------------------------------------
void FourierTransformer::enforceHermitianSymmetry() {
    const long int yHalf = Ysize(*fReal) / 2 + Ysize(*fReal) % 2 - 1;
    const long int zHalf = Zsize(*fReal) / 2 + Zsize(*fReal) % 2 - 1;
    switch (get_array_rank(*fReal)) {
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

MultidimArray<RFLOAT> randomizePhasesBeyond(MultidimArray<RFLOAT> v, int index) {
    FourierTransformer transformer;
    auto &FT = transformer.FourierTransform(v);
    randomizePhasesBeyond(FT, index);
    return transformer.inverseFourierTransform(FT);
}

void randomizePhasesBeyond(MultidimArray<Complex> &FT, int index) {
    const int index2 = index * index;
    FOR_ALL_ELEMENTS_IN_FFTW_TRANSFORM(FT) {
        if (euclidsq(ip, jp, kp) >= index2) {
            Complex &x = direct::elem(FT, i, j, k);
            const RFLOAT mag = abs(x);
            const RFLOAT phase = rnd_unif(0.0, 2.0 * PI);
            x = Complex(mag * cos(phase), mag * sin(phase));
        }
    }
}

// Fourier ring correlation -----------------------------------------------
// from precalculated Fourier Transforms, and without sampling rate etc.
MultidimArray<RFLOAT> getFSC(
    const MultidimArray<Complex> &FT1,
    const MultidimArray<Complex> &FT2
) {
    if (!FT1.sameShape(FT2))
        REPORT_ERROR("fourierShellCorrelation ERROR: MultidimArrays have different shapes!");

    auto num  = MultidimArray<RFLOAT>::zeros(Xsize(FT1));
    auto den1 = MultidimArray<RFLOAT>::zeros(Xsize(FT1));
    auto den2 = MultidimArray<RFLOAT>::zeros(Xsize(FT1));
    auto fsc  = MultidimArray<RFLOAT>::zeros(Xsize(FT1));
    FOR_ALL_ELEMENTS_IN_FFTW_TRANSFORM(FT1) {
        int idx = round(euclid(ip, jp, kp));
        if (idx >= Xsize(FT1)) continue;
        Complex z1 = direct::elem(FT1, i, j, k);
        Complex z2 = direct::elem(FT2, i, j, k);
        RFLOAT absz1 = abs(z1);
        RFLOAT absz2 = abs(z2);
        num.elem(idx) += (conj(z1) * z2).real;
        den1.elem(idx) += absz1 * absz1;
        den2.elem(idx) += absz2 * absz2;
    }

    for (int i = Xinit(fsc); i <= Xlast(fsc); i++) {
        fsc.elem(i) = num.elem(i) / sqrt(den1.elem(i) * den2.elem(i));
    }
    return fsc;

}

MultidimArray<RFLOAT> getFSC(
    MultidimArray<RFLOAT> &m1,
    MultidimArray<RFLOAT> &m2
) {
    FourierTransformer transformer;
    MultidimArray<Complex> FT1 = transformer.FourierTransform(m1);
    MultidimArray<Complex> FT2 = transformer.FourierTransform(m2);
    return getFSC(FT1, FT2);
}

std::pair<MultidimArray<RFLOAT>, MultidimArray<RFLOAT>> getAmplitudeCorrelationAndDifferentialPhaseResidual(
    const MultidimArray<Complex> &FT1, const MultidimArray<Complex> &FT2
) {

    MultidimArray<int> radial_count(Xsize(FT1));
    auto mu1   = MultidimArray<RFLOAT>::zeros(radial_count);
    auto mu2   = MultidimArray<RFLOAT>::zeros(radial_count);
    auto sig1  = MultidimArray<RFLOAT>::zeros(radial_count);
    auto sig2  = MultidimArray<RFLOAT>::zeros(radial_count);
    auto num   = MultidimArray<RFLOAT>::zeros(radial_count);
    auto acorr = MultidimArray<RFLOAT>::zeros(radial_count);  // Amplitude correlation
    auto dpr   = MultidimArray<RFLOAT>::zeros(radial_count);  // Differential phase residual
    FOR_ALL_ELEMENTS_IN_FFTW_TRANSFORM(FT1) {
        // Amplitudes
        int idx = round(euclid(ip, jp, kp));
        if (idx >= Xsize(FT1)) continue;
        RFLOAT abs1 = abs(direct::elem(FT1, i, j, k));
        RFLOAT abs2 = abs(direct::elem(FT2, i, j, k));
        mu1.elem(idx) += abs1;
        mu2.elem(idx) += abs2;
        radial_count.elem(idx)++;

        // Phases
        RFLOAT phas1 = degrees(direct::elem(FT1, i, j, k).arg());
        RFLOAT phas2 = degrees(direct::elem(FT2, i, j, k).arg());
        RFLOAT delta_phas = phas1 - phas2;
        if (delta_phas > +180.0) { delta_phas -= 360.0; }
        if (delta_phas < -180.0) { delta_phas += 360.0; }
        dpr.elem(idx) += delta_phas * delta_phas * (abs1 + abs2);
        num.elem(idx) += abs1 + abs2;
    }

    // Get average amplitudes in each shell for both maps
    for (int i = Xinit(mu1); i <= Xlast(mu1); i++) {
        if (radial_count.elem(i) > 0) {
            mu1.elem(i) /= radial_count.elem(i);
            mu2.elem(i) /= radial_count.elem(i);
            dpr.elem(i) = sqrt(dpr.elem(i) / num.elem(i));
        }
    }

    // Now calculate Pearson's correlation coefficient
    FOR_ALL_ELEMENTS_IN_FFTW_TRANSFORM(FT1) {
        const int idx = round(euclid(ip, jp, kp));
        if (idx >= Xsize(FT1)) continue;
        const RFLOAT z1 = abs(direct::elem(FT1, i, j, k)) - mu1.elem(idx);
        const RFLOAT z2 = abs(direct::elem(FT2, i, j, k)) - mu2.elem(idx);
        acorr.elem(idx) += z1 * z2;
        sig1.elem(idx)  += z1 * z1;
        sig2.elem(idx)  += z2 * z2;
    }

    for (int i = Xinit(acorr); i <= Xlast(acorr); i++) {
        const RFLOAT divisor = sqrt(sig1.elem(i) * sig2.elem(i));
        if (divisor > 0.0) { acorr.elem(i) /= divisor; } else { acorr.elem(i) = 1.0; }
    }

    return {acorr, dpr};
}

std::vector<RFLOAT> cosDeltaPhase(
    const MultidimArray<Complex> &FT1,
    const MultidimArray<Complex> &FT2
) {
    std::vector<RFLOAT> radial_count (Xsize(FT1), 0);
    std::vector<RFLOAT> cos_phi      (Xsize(FT1), 0);

    FOR_ALL_ELEMENTS_IN_FFTW_TRANSFORM(FT1) {
        int idx = round(euclid(ip, jp, kp));
        if (idx >= Xsize(FT1)) continue;

        RFLOAT delta_phase = direct::elem(FT1, i, j, k).arg() 
                           - direct::elem(FT2, i, j, k).arg();
        cos_phi[idx] += cos(delta_phase);
        radial_count[idx]++;
    }

    for (int i = 0; i < cos_phi.size(); i++) {
        if (radial_count[i] > 0) { cos_phi[i] /= radial_count[i]; }
    }

    return cos_phi;
}

std::pair<MultidimArray<RFLOAT>, MultidimArray<RFLOAT>> getAmplitudeCorrelationAndDifferentialPhaseResidual(
    MultidimArray<RFLOAT> &m1, MultidimArray<RFLOAT> &m2
) {
    FourierTransformer transformer;
    MultidimArray<Complex> FT1 = transformer.FourierTransform(m1);
    MultidimArray<Complex> FT2 = transformer.FourierTransform(m2);
    return getAmplitudeCorrelationAndDifferentialPhaseResidual(FT1, FT2);
}

/*
void selfScaleToSizeFourier(long int Ydim, long int Xdim, MultidimArray<RFLOAT>& Mpmem, int nThreads) {

    // Mmem = *this
    // memory for fourier transform output
    // Perform the Fourier transform
    FourierTransformer transformerM;
    transformerM.setThreadsNumber(nThreads);
    MultidimArray<Complex> MmemFourier = transformerM.FourierTransform(Mpmem);

    // Create space for the downsampled image and its Fourier transform
    Mpmem.resize(Xdim, Ydim);
    FourierTransformer transformerMp;
    transformerMp.setReal(Mpmem);
    MultidimArray<Complex> &MpmemFourier = transformerMp.getFourier();
    long int ihalf = std::min((Ysize(MpmemFourier)/2+1),(Ysize(MmemFourier)/2+1));
    long int xsize = std::min((Xsize(MmemFourier)),(Xsize(MpmemFourier)));
    // Init with zero
    MpmemFourier.initZeros();
    for (long int i = 0; i < ihalf; i++)
    for (long int j = 0; j < xsize; j++) {
        MpmemFourier(i, j) = MmemFourier(i, j);
    }
    for (long int i = Ysize(MpmemFourier) - 1, n = 1; n < ihalf - 1; i--, n++) {
        long int ip = Ysize(MmemFourier) - n;
        for (long int j = 0; j<xsize; j++)
            MpmemFourier(i, j) = MmemFourier(ip, j);
    }

    // Transform data
    transformerMp.inverseFourierTransform();
}
*/

void getAbMatricesForShiftImageInFourierTransform(
    MultidimArray<Complex> &in, MultidimArray<Complex> &out,
    RFLOAT oridim, RFLOAT xshift, RFLOAT yshift, RFLOAT zshift
) {
    out.resize(in);
    RFLOAT x, y, z;
    switch (in.getDim()) {

        case 1:
        xshift /= -oridim;
        for (long int i = 0; i < Xsize(in); i++) {
            x = i;
            direct::elem(out, i) = Complex::unit(2 * PI * (x * xshift));
        }
        break;

        case 2:
        xshift /= -oridim;
        yshift /= -oridim;
        for (long int j = 0; j < Xsize(in); j++)
        for (long int i = 0; i < Xsize(in); i++) {
            x = i;
            y = j;
            direct::elem(out, i, j) = Complex::unit(2 * PI * (x * xshift + y * yshift));
        }
        for (long int j = Ysize(in) - 1; j >= Xsize(in); j--) {
        y = j - Ysize(in);
        for (long int i = 0; i < Xsize(in); i++) {
        x = i;
        direct::elem(out, i, j) = Complex::unit(2 * PI * (x * xshift + y * yshift));
        }
        }
        break;

        case 3:
        xshift /= -oridim;
        yshift /= -oridim;
        zshift /= -oridim;
        for (long int k = 0; k < Zsize(in); k++) {
        z = k < Xsize(in) ? k : k - Zsize(in);
        for (long int j = 0; j < Ysize(in); j++) {
        y = j < Xsize(in) ? j : j - Ysize(in);
        for (long int i = 0; i < Xsize(in); i++) {
        x = i;
        direct::elem(out, i, j, k) = Complex::unit(2 * PI * (x * xshift + y * yshift + z * zshift));
        }
        }
        }
        break;

        default:
        REPORT_ERROR("getAbMatricesForShiftImageInFourierTransform ERROR: dimension should be 1, 2 or 3!");

    }
}

void shiftImageInFourierTransformWithTabSincos(
    MultidimArray<Complex> &in,
    MultidimArray<Complex> &out,
    RFLOAT oridim, long int newdim,
    TabSine &tabsin, TabCosine &tabcos,
    RFLOAT xshift, RFLOAT yshift, RFLOAT zshift
) {

    if (&in == &out)
        REPORT_ERROR("shiftImageInFourierTransformWithTabSincos ERROR: Input and output images should be different!");
    // Check size of the input array
    if (Ysize(in) > 1 && Ysize(in) / 2 + 1 != Xsize(in))
        REPORT_ERROR("shiftImageInFourierTransformWithTabSincos ERROR: the Fourier transform should be of an image with equal sizes in all dimensions!");

    long int newhdim = newdim / 2 + 1;
    if (newhdim > Xsize(in))
        REPORT_ERROR("shiftImageInFourierTransformWithTabSincos ERROR: 'newdim' should be equal or smaller than the size of the original array!");

    // Initialise output array
    out.clear();
    switch (in.getDim()) {
        case 2: out.initZeros(newdim,         newhdim); break;
        case 3: out.initZeros(newdim, newdim, newhdim); break;
        default: REPORT_ERROR("shiftImageInFourierTransformWithTabSincos ERROR: dimension should be 2 or 3!");
    }

    if (in.getDim() == 2) {
        xshift /= -oridim;
        yshift /= -oridim;
        if (abs(xshift) < Xmipp::epsilon && abs(yshift) < Xmipp::epsilon) {
            out = windowFourierTransform(in, newdim);
            return;
        }

        FOR_ALL_ELEMENTS_IN_FFTW_TRANSFORM2D(out) {
            RFLOAT dotp = 2.0 * PI * (ip * xshift + jp * yshift);

            Complex X = direct::elem(in, i, j);
            RFLOAT a = tabcos(dotp);
            RFLOAT b = tabsin(dotp);
            Complex Y = Complex(a, b);
            direct::elem(out, i, j) = Complex(
                Y.real * X.real - Y.imag * X.imag, 
                Y.real * X.imag + Y.imag * X.real
            );
        }
    } else if (in.getDim() == 3) {
        xshift /= -oridim;
        yshift /= -oridim;
        zshift /= -oridim;
        if (abs(xshift) < Xmipp::epsilon && abs(yshift) < Xmipp::epsilon && abs(zshift) < Xmipp::epsilon) {
            out = windowFourierTransform(in, newdim);
            return;
        }

        FOR_ALL_ELEMENTS_IN_FFTW_TRANSFORM(out) {
            RFLOAT dotp = 2.0 * PI * (ip * xshift + jp * yshift + kp * zshift);

            Complex X = direct::elem(in, i, j, k);
            RFLOAT a = tabcos(dotp);
            RFLOAT b = tabsin(dotp);
            Complex Y = Complex(a, b);
            direct::elem(out, i, j, k) = Complex(
                Y.real * X.real - Y.imag * X.imag, 
                Y.real * X.imag + Y.imag * X.real
            );
        }
    }
}

// Shift an image through phase-shifts in its Fourier Transform (without pretabulated sine and cosine)
void shiftImageInFourierTransform(
    const MultidimArray<Complex> &in, MultidimArray<Complex> &out,
    RFLOAT oridim, RFLOAT xshift, RFLOAT yshift, RFLOAT zshift
) {
    out.resize(in);
    switch (in.getDim()) {

        case 1:
        xshift /= -oridim;
        if (abs(xshift) < Xmipp::epsilon) {
            out = in;
            return;
        }
        for (long int i = 0; i < Xsize(in); i++) {
            RFLOAT x = i;
            Complex X = direct::elem(in, i);
            Complex Y = Complex::unit(2 * PI * (x * xshift));
            direct::elem(out, i) = Complex(
                X.real * Y.real - X.imag * Y.imag,  // X dot conj Y
                X.imag * Y.real + X.real * Y.imag   // (i conj X) dot Y
            );
        }
        return;

        case 2:
        xshift /= -oridim;
        yshift /= -oridim;
        if (abs(xshift) < Xmipp::epsilon && abs(yshift) < Xmipp::epsilon) {
            out = in;
            return;
        }
        for (long int j = 0; j < Xsize(in); j++)
        for (long int i = 0; i < Xsize(in); i++) {
            RFLOAT x = i, y = j;
            Complex X = direct::elem(in, i, j);
            Complex Y = Complex::unit(2 * PI * (x * xshift + y * yshift));
            direct::elem(out, i, j) = Complex(
                X.real * Y.real - X.imag * Y.imag,  // X dot conj Y
                X.imag * Y.real + X.real * Y.imag   // (i conj X) dot Y
            );
        }
        for (long int j = Ysize(in) - 1; j >= Xsize(in); j--) {
        RFLOAT y = j - Ysize(in);
        for (long int i = 0; i < Xsize(in); i++) {
        RFLOAT x = i;
        Complex X = direct::elem(in, i, j);
        Complex Y = Complex::unit(2 * PI * (x * xshift + y * yshift));
        direct::elem(out, i, j) = Complex(
            X.real * Y.real - X.imag * Y.imag,  // X dot conj Y
            X.real * Y.imag + X.imag * Y.real   // X dot (i conj Y)
        );
        }
        }
        return;

        case 3:
        xshift /= -oridim;
        yshift /= -oridim;
        zshift /= -oridim;
        if (abs(xshift) < Xmipp::epsilon && abs(yshift) < Xmipp::epsilon && abs(zshift) < Xmipp::epsilon) {
            out = in;
            return;
        }
        for (long int k = 0; k < Zsize(in); k++) {
        RFLOAT z = k < Xsize(in) ? k : k - Zsize(in);
        for (long int j = 0; j < Ysize(in); j++) {
        RFLOAT y = j < Xsize(in) ? j : j - Ysize(in);
        for (long int i = 0; i < Xsize(in); i++) {
            RFLOAT x = i;
            Complex X = direct::elem(in, i, j, k);
            Complex Y = Complex::unit(2 * PI * (x * xshift + y * yshift + z * zshift));
            direct::elem(out, i, j, k) = Complex(
                X.real * Y.real - X.imag * Y.imag,  // X dot conj Y
                X.real * Y.imag + X.imag * Y.real   // X dot (i conj Y)
            );
        }
        }
        }
        return;

        default:
        REPORT_ERROR("shiftImageInFourierTransform ERROR: dimension should be 1, 2 or 3!");
    }
}

void shiftImageInFourierTransform(
    MultidimArray<Complex> &in_out,
    RFLOAT oridim, RFLOAT xshift, RFLOAT yshift, RFLOAT zshift
) {
    shiftImageInFourierTransform(in_out, in_out, oridim, xshift, yshift, zshift);
}

MultidimArray<RFLOAT> getSpectrum(const MultidimArray<RFLOAT> &Min, const int spectrum_type) {

    const int xsize = Xsize(Min);
    // Takanori: The above line should be Xsize(Min) / 2 + 1 but for compatibility reasons, I keep this as it is.
    auto spectrum = MultidimArray<RFLOAT>::zeros(xsize);

    FourierTransformer transformer;
    auto &FT = transformer.FourierTransform(Min);
    auto f = spectrum_type == AMPLITUDE_SPECTRUM ? [] (Complex x) { return abs(x); } :
                                                   [] (Complex x) { return norm(x); };
    std::vector<RFLOAT> count (xsize, 0);
    FOR_ALL_ELEMENTS_IN_FFTW_TRANSFORM(FT) {
        long int idx = round(euclid(ip, jp, kp));
        spectrum[idx] += f(direct::elem(FT, i, j, k));
        count[idx]++;
    }

    for (long int i = 0; i < xsize; i++) {
        if (count[i] > 0) { spectrum[i] /= (RFLOAT) count[i]; }
    }
    return spectrum;
}

inline RFLOAT safelydivide(RFLOAT dividend, RFLOAT divisor) {
    if (divisor == 0.0) return 1.0;
    return dividend / divisor;
}

void divideBySpectrum(
    MultidimArray<RFLOAT> &Min,
    const MultidimArray<RFLOAT> &spectrum
) {

    MultidimArray<RFLOAT> inverse_spectrum (spectrum);
    // Is this really better than just doing the division normally?
    for (RFLOAT &x : inverse_spectrum) { x = safelydivide(1.0, x); }
    multiplyBySpectrum(Min, inverse_spectrum);
}

void multiplyBySpectrum(
    MultidimArray<RFLOAT> &Min,
    const MultidimArray<RFLOAT> &spectrum
) {

    FourierTransformer transformer;
    // RFLOAT dim3 = Xsize(Min) * Ysize(Min) * Zsize(Min);
    MultidimArray<Complex> &FT = transformer.FourierTransform(Min);
    FOR_ALL_ELEMENTS_IN_FFTW_TRANSFORM(FT) {
        const long int idx = round(euclid(ip, jp, kp));
        direct::elem(FT, i, j, k) *= spectrum[idx];  // * dim3;
    }
    transformer.inverseFourierTransform();

}

MultidimArray<RFLOAT> whitenSpectrum(
    const MultidimArray<RFLOAT> &Min, int spectrum_type, bool leave_origin_intact
) {

    auto spectrum = getSpectrum(Min, spectrum_type);
    if (!leave_origin_intact) { spectrum[0] = 1.0; }
    MultidimArray<RFLOAT> Mout = Min;
    divideBySpectrum(Mout, spectrum);
    return Mout;
}

MultidimArray<RFLOAT> adaptSpectrum(
    const MultidimArray<RFLOAT> &Min,
    const MultidimArray<RFLOAT> &spectrum_ref,
    int spectrum_type,
    bool leave_origin_intact
) {

    auto spectrum = spectrum_ref / getSpectrum(Min, spectrum_type);
    if (!leave_origin_intact) { spectrum[0] = 1.0; }
    MultidimArray<RFLOAT> Mout = Min;
    multiplyBySpectrum(Mout, spectrum);
    return Mout;
}

/** Kullback-Leibler divergence */
RFLOAT getKullbackLeiblerDivergence(
    MultidimArray<Complex> &Fimg,
    MultidimArray<Complex> &Fref, MultidimArray<RFLOAT> &sigma2,
    MultidimArray<RFLOAT> &p_i, MultidimArray<RFLOAT> &q_i, int highshell, int lowshell
) {
    // First check dimensions are OK
    if (!Fimg.sameShape(Fref))
        REPORT_ERROR("getKullbackLeiblerDivergence ERROR: Fimg and Fref are not of the same shape.");

    if (highshell < 0) { highshell = Xsize(Fimg) - 1; }
    if (lowshell < 0) { lowshell = 0; }

    if (highshell > Xsize(sigma2))
        REPORT_ERROR("getKullbackLeiblerDivergence ERROR: highshell is larger than size of sigma2 array.");

    if (highshell < lowshell)
        REPORT_ERROR("getKullbackLeiblerDivergence ERROR: highshell is smaller than lowshell.");

    // Initialize the histogram
    int histogram_size = 101;
    int histogram_origin = histogram_size / 2;
    RFLOAT sigma_max = 10.0;
    RFLOAT histogram_factor = histogram_origin / sigma_max;
    MultidimArray<int> histogram = MultidimArray<int>::zeros(histogram_size);

    // This way this will work in both 2D and 3D
    FOR_ALL_ELEMENTS_IN_FFTW_TRANSFORM(Fimg) {
        int ires = round(euclid(ip, jp, kp));
        if (ires >= lowshell && ires <= highshell) {
            // Use FT of masked image for noise estimation!
            Complex diff = direct::elem(Fref, i, j, k) - direct::elem(Fimg, i, j, k);
            RFLOAT sigma = sqrt(direct::elem(sigma2, ires));

            // Divide by standard deviation to normalise all the difference
            diff /= sigma;

            // Histogram runs from -10 sigma to +10 sigma
            diff += Complex(sigma_max, sigma_max);

            // Make histogram on the fly
            // Real part
            int ihis = round(diff.real * histogram_factor);
            if (ihis < 0) {
                ihis = 0;
            } else if (ihis >= histogram_size) {
                ihis = histogram_size - 1;
            }
            histogram.elem(ihis)++;
            // Imaginary part
            ihis = round(diff.imag * histogram_factor);
            if (ihis < 0) {
                ihis = 0;
            } else if (ihis > histogram_size) {
                ihis = histogram_size;
            }
            histogram.elem(ihis)++;

        }
    }

    // Normalise the histogram and the discretised analytical Gaussian
    RFLOAT norm = (RFLOAT) histogram.sum();
    RFLOAT gaussnorm = 0.0;
    for (int i = 0; i < histogram_size; i++) {
        RFLOAT x = (RFLOAT)i / histogram_factor;
        gaussnorm += gaussian1D(x - sigma_max, 1.0 , 0.0);
    }

    // Now calculate the actual Kullback-Leibler divergence
    RFLOAT kl_divergence = 0.0;
    p_i.resize(histogram_size);
    q_i.resize(histogram_size);
    for (int i = 0; i < histogram_size; i++) {
        // Data distribution
        p_i.elem(i) = (RFLOAT) histogram.elem(i) / norm;
        // Theoretical distribution
        RFLOAT x = (RFLOAT) i / histogram_factor;
        q_i.elem(i) = gaussian1D(x - sigma_max, 1.0 , 0.0) / gaussnorm;

        if (p_i.elem(i) > 0.0)
            kl_divergence += p_i.elem(i) * log (p_i.elem(i) / q_i.elem(i));
    }
    return kl_divergence / (RFLOAT) histogram_size;
}

void resizeMap(MultidimArray<RFLOAT> &img, int newsize) {

    FourierTransformer transformer;
    const auto &FT = transformer.FourierTransform(img);
    const auto FT2 = windowFourierTransform(FT, newsize);
    img = transformer.inverseFourierTransform(FT2);

}

void applyBFactorToMap(
    MultidimArray<Complex> &FT, int ori_size, RFLOAT bfactor, RFLOAT angpix
) {
    const RFLOAT Nyquist = 0.5 / angpix;
    FOR_ALL_ELEMENTS_IN_FFTW_TRANSFORM(FT) {
        int r2 = euclidsq(ip, jp, kp);
        RFLOAT res = sqrt((RFLOAT) r2) / (ori_size * angpix); // get resolution in 1/Angstrom
        if (res <= Nyquist) {
            // Apply B-factor sharpening until Nyquist, then low-pass filter later on (with a soft edge)
            direct::elem(FT, i, j, k) *= exp(res * res * -bfactor / 4.0);
        } else {
            direct::elem(FT, i, j, k) = 0.0;
        }
    }
}

void applyBFactorToMap(
    MultidimArray<RFLOAT> &img, RFLOAT bfactor, RFLOAT angpix
) {
    FourierTransformer transformer;
    MultidimArray<Complex> &FT = transformer.FourierTransform(img);
    applyBFactorToMap(FT, Xsize(img), bfactor, angpix);
    transformer.inverseFourierTransform();
}

void LoGFilterMap(MultidimArray<Complex> &FT, int ori_size, RFLOAT sigma, RFLOAT angpix) {

    // Calculate sigma in reciprocal pixels (input is in Angstroms) and pre-calculate its square
    // Factor of 1/2 because input is diameter, and filter uses radius
    RFLOAT isigma2 = 0.5 * ori_size * angpix / sigma;
    isigma2 *= isigma2;

    // Gunn Pattern Recognition 32 (1999) 1463-1472
    // The Laplacian filter is: 1/(PI*sigma2)*(r^2/2*sigma2 - 1) * exp(-r^2/(2*sigma2))
    // and its Fourier Transform is: r^2 * exp(-0.5*r2/isigma2);
    // Then to normalise for different scales: divide by isigma2;
    FOR_ALL_ELEMENTS_IN_FFTW_TRANSFORM(FT) {
        RFLOAT r2 = euclidsq((RFLOAT) ip, (RFLOAT) jp, (RFLOAT) kp);
        direct::elem(FT, i, j, k) *= exp(-0.5 * r2 / isigma2) * r2 / isigma2;
    }

}

void LoGFilterMap(MultidimArray<RFLOAT> &img, RFLOAT sigma, RFLOAT angpix) {
    // Make this work for maps (or more likely 2D images) that have unequal X and Y dimensions
    img.setXmippOrigin();
    int my_xsize = Xsize(img);
    int my_ysize = Ysize(img);
    int my_size = std::max(my_xsize, my_ysize);
    if (my_xsize != my_ysize) {
        if (img.getDim() == 2) {
            int my_small_size = std::min(my_xsize, my_ysize);
            const auto stats = computeStats(img);

            img = img.windowed(
                Xmipp::init(my_size), Xmipp::last(my_size),
                Xmipp::init(my_size), Xmipp::last(my_size));
            if (my_small_size == my_xsize) {
                FOR_ALL_ELEMENTS_IN_ARRAY2D(img, i, j) {
                    if (j < Xmipp::init(my_small_size) || j > Xmipp::last(my_small_size))
                        img.elem(i, j) = rnd_gaus(stats.avg, stats.stddev);
                }
            } else {
                FOR_ALL_ELEMENTS_IN_ARRAY2D(img, i, j) {
                    if (i < Xmipp::init(my_small_size) || i > Xmipp::last(my_small_size))
                        img.elem(i, j) = rnd_gaus(stats.avg, stats.stddev);
                }
            }
        } else {
            REPORT_ERROR("lowPassFilterMap: filtering of non-cube maps is not implemented...");
        }
    }
    FourierTransformer transformer;
    MultidimArray<Complex> &FT = transformer.FourierTransform(img);
    LoGFilterMap(FT, Xsize(img), sigma, angpix);
    transformer.inverseFourierTransform();
    img.setXmippOrigin();
    if (my_xsize != my_ysize) {
        if (img.getDim() == 2) {
            img = img.windowed(
                Xmipp::init(my_ysize), Xmipp::last(my_ysize),
                Xmipp::init(my_xsize), Xmipp::last(my_xsize));
        } else {
            REPORT_ERROR("lowPassFilterMap: filtering of non-cube maps is not implemented...");
        }
    }
}

void lowPassFilterMap(
    MultidimArray<Complex> &FT, int ori_size,
    RFLOAT low_pass, RFLOAT angpix, int filter_edge_width, bool do_highpass_instead
) {

    // Which resolution shell is the filter?
    const int ires_filter = round((ori_size * angpix) / low_pass);
    const int filter_edge_halfwidth = filter_edge_width / 2;

    // Soft-edge: from 1 shell less to one shell more:
    const RFLOAT edge_low  = std::max(0.0,                (ires_filter - filter_edge_halfwidth) / (RFLOAT) ori_size); // in 1/pix
    const RFLOAT edge_high = std::min((double) Xsize(FT), (ires_filter + filter_edge_halfwidth) / (RFLOAT) ori_size); // in 1/pix
    const RFLOAT edge_width = edge_high - edge_low;

    // Put a raised cosine from edge_low to edge_high
    FOR_ALL_ELEMENTS_IN_FFTW_TRANSFORM(FT) {
        const RFLOAT r2 = euclidsq(ip, jp, kp);
        const RFLOAT res = sqrt(r2) / ori_size; // get resolution in 1/pixel

        if (do_highpass_instead) {
            if (res < edge_low) {
                direct::elem(FT, i, j, k) = 0.0;
            } else if (res > edge_high) {
                continue;
            } else {
                direct::elem(FT, i, j, k) *= 0.5 * (1.0 - cos(PI * (res - edge_low) / edge_width));
            }
        } else {
            if (res < edge_low) {
                continue;
            } else if (res > edge_high) {
                direct::elem(FT, i, j, k) = 0.0;
            } else {
                direct::elem(FT, i, j, k) *= 0.5 * (1.0 + cos(PI * (res - edge_low) / edge_width));
            }
        }
    }
}

void lowPassFilterMap(
    MultidimArray<RFLOAT> &img, RFLOAT low_pass, RFLOAT angpix, int filter_edge_width
) {
    // Make this work for maps (or more likely 2D images) that have unequal X and Y dimensions
    img.setXmippOrigin();
    int my_xsize = Xsize(img);
    int my_ysize = Ysize(img);
    int my_size = std::max(my_xsize, my_ysize);
    if (my_xsize != my_ysize) {
        if (img.getDim() == 2) {
            int my_small_size = std::min(my_xsize, my_ysize);
            const auto stats = computeStats(img);
            img = img.windowed(
                Xmipp::init(my_size), Xmipp::last(my_size), 
                Xmipp::init(my_size), Xmipp::last(my_size));
            if (my_small_size == my_xsize) {
                FOR_ALL_ELEMENTS_IN_ARRAY2D(img, i, j) {
                    if (i < Xmipp::init(my_small_size) || i > Xmipp::last(my_small_size))
                        img.elem(i, j) = rnd_gaus(stats.avg, stats.stddev);
                }
            } else {
                FOR_ALL_ELEMENTS_IN_ARRAY2D(img, i, j) {
                    if (j < Xmipp::init(my_small_size) || j > Xmipp::last(my_small_size))
                        img.elem(i, j) = rnd_gaus(stats.avg, stats.stddev);
                }
            }
        } else {
            REPORT_ERROR("lowPassFilterMap: filtering of non-cube maps is not implemented...");
        }
    }
    FourierTransformer transformer;
    MultidimArray<Complex> &FT = transformer.FourierTransform(img);
    lowPassFilterMap(FT, Xsize(img), low_pass, angpix, filter_edge_width, false);
    transformer.inverseFourierTransform();
    img.setXmippOrigin();
    if (my_xsize != my_ysize) {
        if (img.getDim() == 2) {
            img = img.windowed(
                Xmipp::init(my_ysize), Xmipp::last(my_xsize),
                Xmipp::init(my_xsize), Xmipp::last(my_ysize));
        } else {
            REPORT_ERROR("lowPassFilterMap: filtering of non-cube maps is not implemented...");
        }
    }
}

void highPassFilterMap(
    MultidimArray<RFLOAT> &img, RFLOAT low_pass, RFLOAT angpix, int filter_edge_width
) {
    FourierTransformer transformer;
    auto &FT = transformer.FourierTransform(img);
    lowPassFilterMap(FT, Xsize(img), low_pass, angpix, filter_edge_width, true);
    transformer.inverseFourierTransform();
}

void directionalFilterMap(
    MultidimArray<Complex> &FT, int ori_size,
    RFLOAT low_pass, RFLOAT angpix, std::string axis, int filter_edge_width
) {

    // Which resolution shell is the filter?
    int ires_filter = round(ori_size * angpix / low_pass);
    int filter_edge_halfwidth = filter_edge_width / 2;

    // Soft-edge: from 1 shell less to one shell more:
    RFLOAT edge_low  = std::max(0.0,                (ires_filter - filter_edge_halfwidth) / (RFLOAT) ori_size); // in 1/pix
    RFLOAT edge_high = std::min((double) Xsize(FT), (ires_filter + filter_edge_halfwidth) / (RFLOAT) ori_size); // in 1/pix
    RFLOAT edge_width = edge_high - edge_low;

    if (axis == "x" || axis == "X") {
        FOR_ALL_ELEMENTS_IN_FFTW_TRANSFORM(FT) {

            RFLOAT r2 = ip * ip;
            /// BUG: sqrt(r2 / ori_size) or sqrt(r2) / ori_size ?
            RFLOAT res = sqrt(r2 / (RFLOAT) ori_size); // get resolution in 1/pixel

            if (res < edge_low) {
                continue;
            } else if (res > edge_high) {
                direct::elem(FT, i, j, k) = 0.0;
            } else {
                direct::elem(FT, i, j, k) *= raised_cos(PI * (res - edge_low) / edge_width);
            }
        }
    } else if (axis == "y" || axis == "Y") {
        FOR_ALL_ELEMENTS_IN_FFTW_TRANSFORM(FT) {

            RFLOAT r2 = jp * jp;
            RFLOAT res = sqrt(r2) / ori_size; // get resolution in 1/pixel

            if (res < edge_low) {
                continue;
            } else if (res > edge_high) {
                direct::elem(FT, i, j, k) = 0.0;
            } else {
                direct::elem(FT, i, j, k) *= raised_cos(PI * (res - edge_low) / edge_width);
            }
        }
    } else if  (axis == "z" || axis == "Z") {
        FOR_ALL_ELEMENTS_IN_FFTW_TRANSFORM(FT) {

            RFLOAT r2 = kp * kp;
            RFLOAT res = sqrt(r2) / ori_size; // get resolution in 1/pixel

            if (res < edge_low) {
                continue;
            } else if (res > edge_high) {
                direct::elem(FT, i, j, k) = 0.;
            } else {
                direct::elem(FT, i, j, k) *= raised_cos(PI * (res - edge_low) / edge_width);
            }
        }
    }
}

void directionalFilterMap(
    MultidimArray<RFLOAT> &img,
    RFLOAT low_pass, RFLOAT angpix,
    std::string axis, int filter_edge_width
) {
    // Make this work for maps (or more likely 2D images) that have unequal X and Y dimensions
    img.setXmippOrigin();
    int my_xsize = Xsize(img);
    int my_ysize = Ysize(img);
    int my_size = std::max(my_xsize, my_ysize);
    if (my_xsize != my_ysize) {
        if (img.getDim() == 2) {
            int my_small_size = std::min(my_xsize, my_ysize);
            const auto stats = computeStats(img);
            img = img.windowed(
                Xmipp::init(my_size), Xmipp::last(my_size),
                Xmipp::init(my_size), Xmipp::last(my_size));
            if (my_small_size == my_xsize) {
                FOR_ALL_ELEMENTS_IN_ARRAY2D(img, i, j) {
                    if (i < Xmipp::init(my_small_size) || i > Xmipp::last(my_small_size))
                        img.elem(i, j) = rnd_gaus(stats.avg, stats.stddev);
                }
            } else {
                FOR_ALL_ELEMENTS_IN_ARRAY2D(img, i, j) {
                    if (j < Xmipp::init(my_small_size) || j > Xmipp::last(my_small_size))
                        img.elem(i, j) = rnd_gaus(stats.avg, stats.stddev);
                }
            }
        } else {
            REPORT_ERROR("lowPassFilterMap: filtering of non-cube maps is not implemented...");
        }
    }
    FourierTransformer transformer;
    auto &FT = transformer.FourierTransform(img);
    directionalFilterMap(FT, Xsize(img), low_pass, angpix, axis, filter_edge_width);
    transformer.inverseFourierTransform();
    img.setXmippOrigin();
    if (my_xsize != my_ysize) {
        if (img.getDim() == 2) {
            img = img.windowed(
                Xmipp::init(my_ysize), Xmipp::last(my_ysize),
                Xmipp::init(my_xsize), Xmipp::last(my_xsize));
        } else {
            REPORT_ERROR("lowPassFilterMap: filtering of non-cube maps is not implemented...");
        }
    }
}

void applyBeamTilt(
    const MultidimArray<Complex> &Fin, MultidimArray<Complex> &Fout,
    RFLOAT beamtilt_x, RFLOAT beamtilt_y,
    RFLOAT wavelength, RFLOAT Cs, RFLOAT angpix, int ori_size
) {
    Fout = Fin;
    selfApplyBeamTilt(Fout, beamtilt_x, beamtilt_y, wavelength, Cs, angpix, ori_size);
}

void selfApplyBeamTilt(
    MultidimArray<Complex> &Fimg, RFLOAT beamtilt_x, RFLOAT beamtilt_y,
    RFLOAT wavelength, RFLOAT Cs, RFLOAT angpix, int ori_size
) {
    if (Fimg.getDim() != 2)
        REPORT_ERROR("applyBeamTilt can only be done on 2D Fourier Transforms!");

    RFLOAT boxsize = angpix * ori_size;
    RFLOAT factor = 0.360 * Cs * 10000000 * wavelength * wavelength / (boxsize * boxsize * boxsize);
    FOR_ALL_ELEMENTS_IN_FFTW_TRANSFORM2D(Fimg) {
        RFLOAT delta_phase = factor * (ip * ip + jp * jp) * (ip * beamtilt_x + jp * beamtilt_y);
        Complex A = direct::elem(Fimg, i, j);
        RFLOAT mag = sqrt(A.real * A.real + A.imag * A.imag);
        RFLOAT phas = atan2(A.imag, A.real) + radians(delta_phase); // apply phase shift!
        direct::elem(Fimg, i, j) = Complex(mag * cos(phas), mag * sin(phas));
    }
}

void selfApplyBeamTilt(
    MultidimArray<Complex> &Fimg,
    RFLOAT beamtilt_x, RFLOAT beamtilt_y,
    RFLOAT beamtilt_xx, RFLOAT beamtilt_xy, RFLOAT beamtilt_yy,
    RFLOAT wavelength, RFLOAT Cs, RFLOAT angpix, int ori_size
) {
    if (Fimg.getDim() != 2)
        REPORT_ERROR("applyBeamTilt can only be done on 2D Fourier Transforms!");

    RFLOAT boxsize = angpix * ori_size;
    RFLOAT factor = 0.360 * Cs * 10000000 * wavelength * wavelength / (boxsize * boxsize * boxsize);

    FOR_ALL_ELEMENTS_IN_FFTW_TRANSFORM2D(Fimg) {
        // (ip beamtilt_x + jp beamtilt_y) ** 2
        RFLOAT q = beamtilt_xx * ip * ip + 2.0 * beamtilt_xy * ip * jp + beamtilt_yy * jp * jp;

        RFLOAT delta_phase = factor * q * (ip * beamtilt_x + jp * beamtilt_y);
        Complex X = direct::elem(Fimg, i, j);
        RFLOAT mag = X.abs();
        RFLOAT phas = X.arg() + radians(delta_phase); // apply phase shift!
        direct::elem(Fimg, i, j) = Complex(mag * cos(phas), mag * sin(phas));
    }
}

MultidimArray<RFLOAT> padAndFloat2DMap(const MultidimArray<RFLOAT> &v, int factor) {

    // Check dimensions
    const auto dimensions = v.getDimensions();
    if (dimensions[2] > 1 || dimensions[3] > 1)
        REPORT_ERROR("fftw.cpp::padAndFloat2DMap(): ERROR MultidimArray should be 2D.");
    if (dimensions[0] * dimensions[1] <= 16)
        REPORT_ERROR("fftw.cpp::padAndFloat2DMap(): ERROR MultidimArray is too small.");
    if (factor <= 1)
        REPORT_ERROR("fftw.cpp::padAndFloat2DMap(): ERROR Padding factor should be larger than 1.");

    // Calculate background and border values
    RFLOAT bg_val, bg_pix, bd_val, bd_pix;
    bg_val = bg_pix = bd_val = bd_pix = 0.0;
    for (long int j = 0; j < Ysize(v); j++)
    for (long int i = 0; i < Xsize(v); i++) {
        bg_val += direct::elem(v, i, j);
        bg_pix += 1.0;
        if (i == 0 || j == 0 || i == Xsize(v) - 1 || j == Ysize(v) - 1) {
            bd_val += direct::elem(v, i, j);
            bd_pix += 1.0;
        }
    }
    if (bg_pix < 1.0 || bd_pix < 1.0) {
        REPORT_ERROR("fftw.cpp::padAndFloat2DMap(): ERROR MultidimArray is too small.");
    }
    bg_val /= bg_pix;
    bd_val /= bd_pix;
    // DEBUG
    // std::cout << "bg_val = " << bg_val << ", bg_pix = " << bg_pix << std::endl;
    // std::cout << "bd_val = " << bd_val << ", bd_pix = " << bd_pix << std::endl;

    // Pad and float output MultidimArray (2× original size by default)
    long int box_len = std::max(dimensions[0], dimensions[1]) * factor;
    MultidimArray<RFLOAT> out (box_len, box_len);
    out = bd_val - bg_val;
    out.setXmippOrigin();
    for (long int j = 0; j < Ysize(v); j++)
    for (long int i = 0; i < Xsize(v); i++) {
        direct::elem(out, i, j) = direct::elem(v, i, j) - bg_val;
    }
    return out;
}

MultidimArray<RFLOAT> amplitudeOrPhaseMap(
    const MultidimArray<RFLOAT> &v, int output_map_type
) {
    // Pad and float
    MultidimArray<RFLOAT> amp = padAndFloat2DMap(v);
    if (Xsize(amp) != Ysize(amp) || Zsize(amp) > 1 || Nsize(amp) > 1)
        REPORT_ERROR("fftw.cpp::amplitudeOrPhaseMap(): ERROR MultidimArray should be 2D square.");
    long int XYdim = Xsize(amp);

    // Fourier Transform
    FourierTransformer transformer;
    auto FT = transformer.FourierTransform(amp);
    CenterFFTbySign(FT);

    const auto f = output_map_type == AMPLITUDE_MAP ? +[] (Complex x) { return         x.abs() ; } :
                   output_map_type == PHASE_MAP     ? +[] (Complex x) { return degrees(x.arg()); } :
                   nullptr;
    if (!f) REPORT_ERROR("fftw.cpp::amplitudeOrPhaseMap(): ERROR Unknown type of output map.");

    // Write to output files
    amp.setXmippOrigin();
    amp.initZeros(XYdim, XYdim);
    long int maxr2 = (XYdim - 1) * (XYdim - 1) / 4;
    FOR_ALL_ELEMENTS_IN_FFTW_TRANSFORM2D(FT) {
        if (
            ip > Xinit(amp) && ip < Xlast(amp) &&
            jp > Yinit(amp) && jp < Ylast(amp) &&
            euclidsq(ip, jp) < maxr2
        ) {
            RFLOAT val = f(FFTW::elem(FT, ip, jp));
            amp.elem(ip, jp) = amp.elem(-ip, -jp) = val;
        }
    }
    amp.elem(0, 0) = 0.0;
    return amp;
}

void helicalLayerLineProfile(
    const MultidimArray<RFLOAT> &v, std::string title, std::string fn_eps
) {
    std::vector<RFLOAT> ampl_list, ampr_list, nr_pix_list;

    // TODO: DO I NEED TO ROTATE THE ORIGINAL MULTIDIMARRAY BY 90 DEGREES ?

    // Pad and float
    MultidimArray<RFLOAT> out = padAndFloat2DMap(v);
    if (
        Xsize(out) != Ysize(out) || Zsize(out) > 1 || Nsize(out) > 1
    ) REPORT_ERROR("fftw.cpp::helicalLayerLineProfile(): ERROR MultidimArray should be 2D square.");
    long int XYdim = Xsize(out);

    // Fourier Transform
    FourierTransformer transformer;
    auto FT = transformer.FourierTransform(out);
    CenterFFTbySign(FT);

    // Statistics
    out.setXmippOrigin();
    long int maxr2 = (XYdim - 1) * (XYdim - 1) / 4;
    ampl_list.resize(Xsize(FT) + 2);
    ampr_list.resize(Xsize(FT) + 2);
    nr_pix_list.resize(Xsize(FT) + 2);
    for (int ii = 0; ii < ampl_list.size(); ii++)
        ampl_list[ii] = ampr_list[ii] = nr_pix_list[ii] = 0.0;

    FOR_ALL_ELEMENTS_IN_FFTW_TRANSFORM2D(FT) {
        if (euclidsq(ip, jp) < maxr2 && ip > 0) {
            nr_pix_list[jp] += 1.0;
            ampl_list[jp] += FFTW::elem(FT,  ip, jp).abs();
            ampr_list[jp] += FFTW::elem(FT, -ip, jp).abs();
        }
    }
    CDataSet dataSetAmpl, dataSetAmpr;
    RFLOAT linewidth = 1.0;
    std::string figTitle = "Helical Layer Line Profile - " + title;
    std::string yTitle = "Reciprocal pixels (padded box size = " + integerToString(XYdim) + ")";
    for (int ii = 0; ii < 3 * ampl_list.size() / 4 + 1; ii++) {
        if (nr_pix_list[ii] < 1.0) break;  // CHECK: IS THIS CORRECT?
        dataSetAmpl.AddDataPoint(CDataPoint(ii, log(ampl_list[ii] / nr_pix_list[ii])));
        dataSetAmpr.AddDataPoint(CDataPoint(ii, log(ampr_list[ii] / nr_pix_list[ii])));
    }
    dataSetAmpl.SetDrawMarker(false);
    dataSetAmpl.SetLineWidth(linewidth);
    dataSetAmpl.SetDatasetColor(1.0, 0.0, 0.0);
    dataSetAmpl.SetDatasetTitle("ln(amplitudes) (left)");
    dataSetAmpr.SetDrawMarker(false);
    dataSetAmpr.SetLineWidth(linewidth);
    dataSetAmpr.SetDatasetColor(0.0, 1.0, 0.0);
    dataSetAmpr.SetDatasetTitle("ln(amplitudes) (right)");
    CPlot2D *plot2D = new CPlot2D(figTitle);
    plot2D->SetXAxisSize(600);
    plot2D->SetYAxisSize(400);
    plot2D->SetXAxisTitle(yTitle);
    plot2D->SetYAxisTitle("ln(amplitudes)");
    plot2D->AddDataSet(dataSetAmpl);
    plot2D->AddDataSet(dataSetAmpr);
    plot2D->OutputPostScriptPlot(fn_eps);
    delete plot2D;
}

MultidimArray<RFLOAT> generateBinaryHelicalFourierMask(
    long int xdim, long int ydim, long int zdim,
    std::vector<RFLOAT> exclude_begin,
    std::vector<RFLOAT> exclude_end,
    RFLOAT angpix
) {
    if (exclude_begin.size() != exclude_end.size())
        REPORT_ERROR("BUG: generateHelicalFourierMask: provide start-end resolutions for each shell.");

    auto mask = MultidimArray<RFLOAT>::ones(xdim, ydim, zdim);

    const bool is_2d = mask.getDim() == 2;
    FOR_ALL_ELEMENTS_IN_FFTW_TRANSFORM(mask) {
        RFLOAT res = is_2d ?
            jp == 0 ? 999.0 : Ysize(mask) * angpix / fabs(jp) : // helical axis along X-axis, so only jp matters!
            kp == 0 ? 999.0 : Zsize(mask) * angpix / fabs(kp);  // helical axis along Z-axis, so only kp matters!

        for (int ishell = 0; ishell < exclude_begin.size(); ishell++) {
            if (res <= exclude_begin[ishell] && res >= exclude_end[ishell]) {
                direct::elem(mask, i, j, k) = 0.0;
            }
        }
    }
    return mask;
}
