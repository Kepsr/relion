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
 * Authors:    Roberto Marabini	(roberto@cnb.csic.es)
 *             Carlos Oscar S. Sorzano	(coss@cnb.csic.es)
 *
 * Unidad de Bioinformatica of Centro Nacional de Biotecnologia, CSIC
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
 * 02111-1307 USA
 *
 * All comments concerning this program package may be sent to the
 * e-mail address 'xmipp@cnb.csic.es'
 ***************************************************************************/

#ifndef __RELIONFFTW_H
#define __RELIONFFTW_H

#include <fftw3.h>
#include "src/multidim_array.h"
#include "src/funcs.h"
#include "src/tabfuncs.h"
#include "src/complex.h"
#include "src/CPlot2D.h"

#ifdef RELION_SINGLE_PRECISION

typedef fftwf_complex FFTW_COMPLEX;
#define FFTW_PLAN            fftwf_plan
#define FFTW_PLAN_DFT        fftwf_plan_dft
#define FFTW_PLAN_DFT_R2C    fftwf_plan_dft_r2c
#define FFTW_PLAN_DFT_C2R    fftwf_plan_dft_c2r
#define FFTW_EXECUTE_DFT_R2C fftwf_execute_dft_r2c
#define FFTW_EXECUTE_DFT_C2R fftwf_execute_dft_c2r
#define FFTW_CLEANUP         fftwf_cleanup
#define FFTW_DESTROY_PLAN    fftwf_destroy_plan

#else

typedef fftw_complex FFTW_COMPLEX;
#define FFTW_PLAN            fftw_plan
#define FFTW_PLAN_DFT        fftw_plan_dft
#define FFTW_PLAN_DFT_R2C    fftw_plan_dft_r2c
#define FFTW_PLAN_DFT_C2R    fftw_plan_dft_c2r
#define FFTW_EXECUTE_DFT_R2C fftw_execute_dft_r2c
#define FFTW_EXECUTE_DFT_C2R fftw_execute_dft_c2r
#define FFTW_CLEANUP         fftw_cleanup
#define FFTW_DESTROY_PLAN    fftw_destroy_plan

#endif

//#define TIMING_FFTW
#ifdef TIMING_FFTW
    #include "src/time.h"
    extern Timer timer_fftw;
#endif

#ifdef FAST_CENTERFFT	// defined if ALTCPU=on *AND* Intel Compiler used
#include "src/acc/cpu/cuda_stubs.h"
#include "src/acc/settings.h"
#include "src/acc/cpu/cpu_settings.h"
#include "src/acc/cpu/cpu_kernels/helper.h"
#include <tbb/parallel_for.h>
#endif

/** @defgroup FourierW FFTW Fourier transforms
  * @ingroup DataLibrary
  */

/** For all direct elements in the complex array in FFTW format.
 *
 * This macro can be used to write loops over the indices of a Fourier-space volume.
 * [ Xsize(FT) = Xsize(V) / 2 + 1 ] [ Uncentered DC component ]
 * The indices i, j, k range the data "physically".
 * The indices ip, jp, kp range the data "logically".
 * It will also work for 1D and 2D FFTW transforms.
 *
 * @code
 * FOR_ALL_ELEMENTS_IN_FFTW_TRANSFORM(V) {
 *     int r2 = ip * ip + jp * jp + kp * kp;
 *     std::cout << "element at physical coords: " << i << " " << j << " " << k << " has value: " << direct::elem(m, i, j, k) << std::endl;
 *     std::cout << "its logical coords are: " << ip << " " << jp << " " << kp << std::endl;
 *     std::cout << "its distance from the origin = " << sqrt(r2) << std::endl;
 * }
 * @endcode
 */
#define FOR_ALL_ELEMENTS_IN_FFTW_TRANSFORM(V) \
    for (long int k = 0, kp = 0; k < Zsize(V); k++, kp = Fourier::K(V, k)) \
    for (long int j = 0, jp = 0; j < Ysize(V); j++, jp = Fourier::J(V, j)) \
    for (long int i = 0, ip = 0; i < Xsize(V); i++, ip = Fourier::I(V, i))

namespace Fourier {

    template <typename T>
    inline int I(const MultidimArray<T> &arr, int i) { return i; }
    template <typename T>
    inline int J(const MultidimArray<T> &arr, int j) { return j < Xsize(arr) ? j : j - Ysize(arr); }
    template <typename T>
    inline int K(const MultidimArray<T> &arr, int k) { return k < Xsize(arr) ? k : k - Zsize(arr); }

};

/** For all direct elements in the complex array in FFTW format.
 * The same as above, but now only for 2D images
 */
// FOR_i_j_ip_jp_IN_FFTW_TRANSFORM2D
#define FOR_ALL_ELEMENTS_IN_FFTW_TRANSFORM2D(V) \
    for (long int j = 0, jp = 0; j < Ysize(V); j++, jp = Fourier::J(V, j)) \
    for (long int i = 0, ip = 0; i < Xsize(V); i++, ip = Fourier::I(V, i))

namespace FFTW {

    /** Logical access to FFTW element in 3D
     *
     * @code
     * FFTW::elem(V, -1, -2, 1) = 1;
     * val = FFTW::elem(V, -1, -2, 1);
     * @endcode
     */
    template <typename T>
    T& elem(MultidimArray<T> &V, int ip, int jp, int kp) {
        return direct::elem(V, ip,
                            jp < 0 ? jp + Ysize(V) : jp,
                            kp < 0 ? kp + Zsize(V) : kp);
    }

    template <typename T>
    const T& elem(const MultidimArray<T> &V, int ip, int jp, int kp) {
        return direct::elem(V, ip,
                            jp < 0 ? jp + Ysize(V) : jp,
                            kp < 0 ? kp + Zsize(V) : kp);
    }

    /** Logical access to FFTW element in 2D
     *
     * @code
     * FFTW::elem(V, -2, 1) = 1;
     * val = FFTW::elem(V, -2, 1);
     * @endcode
     */
    template <typename T>
    T& elem(MultidimArray<T> &V, int ip, int jp) {
        return direct::elem(V, ip,
                            jp < 0 ? jp + Ysize(V) : jp);
    }

    template <typename T>
    const T& elem(const MultidimArray<T> &V, int ip, int jp) {
        return direct::elem(V, ip,
                            jp < 0 ? jp + Ysize(V) : jp);
    }

}

/** Fourier Transformer class.
 * @ingroup FourierW
 *
 * The memory for the Fourier transform is handled by this object.
 * However, the memory for the real space image is handled externally
 * and this object only has a pointer to it.
 *
 * Here you have an example of use
 * @code
 * FourierTransformer transformer;
 * MultidimArray<Complex> &Vfft = transformer.FourierTransform(V());
 * MultidimArray<RFLOAT> Vmag;
 * Vmag.resize(Vfft);
 * FOR_ALL_ELEMENTS_IN_ARRAY3D(Vmag, i, j, k)
    *	Vmag(k, i, j) = 20 * log10(abs(Vfft(k, i, j)));
 * @endcode
 */
class FourierTransformer {

    public:

    template <typename T>
    static inline int get_array_rank(const MultidimArray<T> &arr) {
        return Zsize(arr) == 1 ? Ysize(arr) == 1 ? 1 : 2 : 3;
    }

    static int *new_n(const MultidimArray<RFLOAT> &arr, int rank) {
        int *n = new int[rank];
        const std::array<long int, 3> v { Xsize(arr), Ysize(arr), Zsize(arr) };
        for (int i = 0; i < rank; i++) { n[rank - i - 1] = v[i]; }
        return n;
    }

    // Pointer to a real array
    const MultidimArray<RFLOAT> *fReal;

    // Pointer to a complex array
    const MultidimArray<Complex> *fComplex;

    // Fourier array
    MultidimArray<Complex> fFourier;

    // fftw Forward plan
    FFTW_PLAN fPlanForward;

    // fftw Backward plan
    FFTW_PLAN fPlanBackward;

    bool plans_are_set;

    public:

    /** Default constructor */
    FourierTransformer();

    /** Destructor */
    ~FourierTransformer();

    /** Compute the Fourier transform of a MultidimArray, 2D and 3D.
     */
    template <typename T>
    MultidimArray<tComplex<T>> &FourierTransform(const MultidimArray<T> &v) {
        setReal(v);
        Transform(FFTW_FORWARD);
        return fFourier;
    }

    /** Compute the Fourier transform.
        The data is taken from the matrix with which the object was
        created. */
    void FourierTransform();

    /** Enforce Hermitian symmetry.
        If the Fourier transform risks losing Hermitian symmetry,
        use this function to restore it. */
    void enforceHermitianSymmetry(MultidimArray<Complex> &array);
    /// WARNING: Unused

    /** Compute the inverse Fourier transform.
        The result is stored in the same real data that was passed for
        the forward transform. The Fourier coefficients are taken from
        the internal Fourier coefficients */
    void inverseFourierTransform();

    /** Compute the inverse Fourier transform.
        New data is provided for the Fourier coefficients and the output
        can be any MultiDimarray (1D, 2D or 3D). */
    template <typename T>
    MultidimArray<T> inverseFourierTransform(const MultidimArray<tComplex<T> > &V) {
        MultidimArray<T> IFT (2 * (Xsize(V) - 1), Ysize(V), Zsize(V));
        setReal(IFT);
        setFourier(V);
        Transform(FFTW_BACKWARD);
        return IFT;
    }

    /** Get Fourier coefficients. */
    MultidimArray<Complex>& getFourier() { return fFourier; }

    /** Return a complete Fourier transform (two halves).
    */
    MultidimArray<Complex> getCompleteFourier() {
        MultidimArray<Complex> V;
        V.reshape(*fReal);
        switch (get_array_rank(V)) {
            case 1:
            for (long int i = 0; i < Xsize(V); i++) {
                direct::elem(V, i) = i < Xsize(fFourier) ? direct::elem(fFourier, i) :
                    conj(direct::elem(fFourier, Xsize(*fReal) - i));
            } break;
            case 2:
            for (long int j = 0; j < Ysize(V); j++)
            for (long int i = 0; i < Xsize(V); i++) {
                direct::elem(V, i, j) = i < Xsize(fFourier) ? direct::elem(fFourier, i, j) :
                    conj(direct::elem(fFourier, Xsize(*fReal) - i, (Ysize(*fReal) - j) % Ysize(*fReal)));
            } break;
            case 3:
            for (long int k = 0; k < Zsize(V); k++)
            for (long int j = 0; j < Ysize(V); j++)
            for (long int i = 0; i < Xsize(V); i++) {
                direct::elem(V, i, j, k) = i < Xsize(fFourier) ? direct::elem(fFourier, i, j, k) :
                    conj(direct::elem(fFourier, Xsize(*fReal) - i, (Ysize(*fReal) - j) % Ysize(*fReal), (Zsize(*fReal) - k) % Zsize(*fReal)));
            } break;
        }
        return V;
    }

    /** Set one half of the FT in fFourier from the input complete Fourier transform (two halves).
        The fReal and fFourier already should have the right sizes
    */
    template <typename T>
    void setFromCompleteFourier(const MultidimArray<T> &V);

    public:

    /* Pointer to the real array with which the plan was computed */
    RFLOAT *dataPtr;

    /* Pointer to the complex array with which the plan was computed */
    Complex *complexDataPtr;

    /* Initialise all pointers to NULL */
    void init();

    /** Clear object */
    void clear();

    /** This calls fftw_cleanup.
    */
    void cleanup();

    /** Destroy both forward and backward fftw planes (mutex locked */
    void destroyPlans();

    /** Computes the transform, specified in Init() function
        If normalization=true the forward transform is normalized
        (no normalization is made in the inverse transform)
        If normalize=false no normalization is performed and therefore
        the image is scaled by the number of pixels.
    */
    void Transform(int sign);

    /** Get the Multidimarray that is being used as input. */
    const MultidimArray<RFLOAT> &getReal() const { return *fReal; }
    const MultidimArray<Complex> &getComplex() const { return *fComplex; }

    /** Set a Multidimarray for input.
        The data of img will be the one of fReal. In forward
        transforms it is not modified, but in backward transforms,
        the result will be stored in img. This means that the size
        of img cannot change between calls. */
    void setReal(const MultidimArray<RFLOAT> &img);

    /** Set a Multidimarray for input.
        The data of img will be the one of fComplex. In forward
        transforms it is not modified, but in backward transforms,
        the result will be stored in img. This means that the size
        of img cannot change between calls. */
    void setReal(const MultidimArray<Complex> &img);

    // Call FFTW_PLAN_DFT_R2C and FFTW_PLAN_DFT_C2R
    void computePlans(const MultidimArray<RFLOAT> &input);

    // Call FFTW_PLAN_DFT
    void computePlans(const MultidimArray<Complex> &input);

    /** Set a Multidimarray for the Fourier transform.
        The values of the input array are copied in the internal array.
        It is assumed that the container for the real image as well as
        the one for the Fourier array are already resized.
        No plan is updated. */
    void setFourier(const MultidimArray<Complex> &imgFourier);

};

// Randomize phases beyond the given F-space shell (index) of R-space input image
MultidimArray<RFLOAT> randomizePhasesBeyond(MultidimArray<RFLOAT> I, int index);

// Randomize phases beyond the given F-space shell (index) of F-space input image
void randomizePhasesBeyond(MultidimArray<Complex> &FT, int index);

template <typename T>
MultidimArray<T>& CenterFFTbySign(MultidimArray<T> &v) {
    // This technique does not work when the dimensions of iFFT(v) are of odd size.
    // Unfortunately, this cannot be checked within this function...
    // Forward and backward shifts are equivalent.
    FOR_ALL_ELEMENTS_IN_ARRAY3D(v, i, j, k) {
        if ((i ^ j ^ k) & 1)  // If the least significant bit of i ^ j ^ k is 1 (i.e. if i + j + k is odd)
            direct::elem(v, i, j, k) *= -1;
    }
    return v;
}

/** Center an array,
 * so that its zero-frequency component lies at the origin.
 * Like numpy.fft.fftshift.
 */
template <typename T>
void CenterFFT(MultidimArray<T> &v, int sign = +1);

// Window an FFTW-centered Fourier-transform to a given size
template<class T>
MultidimArray<T> windowFourierTransform(const MultidimArray<T> &in, long int newdim) {
    // Check size of the input array
    if (Ysize(in) > 1 && Ysize(in) / 2 + 1 != Xsize(in))
        REPORT_ERROR("windowFourierTransform ERROR: the Fourier transform should be of an image with equal sizes in all dimensions!");
    long int newhdim = newdim / 2 + 1;

    // If same size, just return input
    // Sjors 5dec2017: only check for xdim is not enough, even/off ydim leaves ambiguity for dim>1
    if (newhdim == Xsize(in) && newdim == Ysize(in))
        return in;

    const int rank = in.getDim();
    if (rank < 1 || rank > 3)
        REPORT_ERROR("windowFourierTransform ERROR: dimension should be 1, 2 or 3!");

    // Apply a windowing operation
    // Initialise output array
    auto out = MultidimArray<T>::zeros(
        newhdim, rank >= 2 ? newdim : 1, rank >= 3 ? newdim : 1
    );

    if (newhdim > Xsize(in)) {
        long int max_r2 = (Xsize(in) - 1) * (Xsize(in) - 1);
        FOR_ALL_ELEMENTS_IN_FFTW_TRANSFORM(in) {
            // Make sure windowed FT has nothing in the corners, otherwise we end up with an asymmetric FT!
            if (euclidsq(ip, jp, kp) <= max_r2)
                FFTW::elem(out, ip, jp, kp) = FFTW::elem(in, ip, jp, kp);
        }
    } else {
        FOR_ALL_ELEMENTS_IN_FFTW_TRANSFORM(out) {
            FFTW::elem(out, ip, jp, kp) = FFTW::elem(in, ip, jp, kp);
        }
    }
    return out;
}

// A resize operation in Fourier-space (i.e. changing the sampling of the Fourier Transform) by windowing in real-space
// If recenter=true, the real-space array will be recentered to have its origin at the origin of the FT
template<class T>
void resizeFourierTransform(const MultidimArray<T> &in, MultidimArray<T> &out, long int newdim, bool do_recenter=true) {
    // Check size of the input array
    if (Ysize(in) > 1 && Ysize(in) / 2 + 1 != Xsize(in))
        REPORT_ERROR("windowFourierTransform ERROR: the Fourier transform should be of an image with equal sizes in all dimensions!");
    long int newhdim = newdim / 2 + 1;
    long int olddim = 2 * (Xsize(in) - 1);

    // If same size, just return input
    if (newhdim == Xsize(in)) {
        out = in;
        return;
    }

    // Otherwise apply a windowing operation
    FourierTransformer transformer;
    long int x0, y0, z0, xF, yF, zF;
    x0 = y0 = z0 = Xmipp::init(newdim);
    xF = yF = zF = Xmipp::last(newdim);

    // in might be a MultidimArray<RFLOAT>,
    // so we make a copy Fin of in, of type MultidimArray<Complex>
    MultidimArray<Complex> Fin;
    Fin.reshape(Xsize(in), Ysize(in), Zsize(in));
    for (long int i = 0; i < in.size(); i++) {
        Fin.data[i] = in.data[i];
    }

    // Initialise output array
    MultidimArray<RFLOAT> Min;
    switch (in.getDim()) {

        case 1:
        Min.reshape(olddim);
        y0 = yF = z0 = zF = 0;
        break;

        case 2:
        Min.reshape(olddim, olddim);
        z0 = zF = 0;
        break;

        case 3:
        Min.reshape(olddim, olddim, olddim);
        break;

        default:
        REPORT_ERROR("resizeFourierTransform ERROR: dimension should be 1, 2 or 3!");

    }

    Min = transformer.inverseFourierTransform(Fin).setXmippOrigin();
    if (do_recenter) {
        CenterFFT(Min, -1);
    }

    // Now do the actual windowing in real-space
    Min = Min.windowed(x0, xF, y0, yF, z0, zF).setXmippOrigin();

    // If upsizing: mask the corners to prevent aliasing artefacts
    if (newdim > olddim) {
        FOR_ALL_ELEMENTS_IN_ARRAY3D(Min, i, j, k) {
            if (euclidsq(i, j, k) > olddim * olddim / 4) {
                Min.elem(i, j, k) = 0.0;
            }
        }
    }

    // Recenter FFT back again
    if (do_recenter)
        CenterFFT(Min, +1);

    // And do the inverse Fourier transform
    transformer.clear();
    out = transformer.FourierTransform(Min);
}

/** Fourier-Ring-Correlation between two multidimArrays using FFT
 * From precalculated Fourier Transforms
 * Simpler I/O than above.
 */
MultidimArray<RFLOAT> getFSC(
    const MultidimArray<Complex> &FT1,
    const MultidimArray<Complex> &FT2
);

/** Fourier-Ring-Correlation between two multidimArrays using FFT
 * @ingroup FourierOperations
 * Simpler I/O than above.
 */
MultidimArray<RFLOAT> getFSC(
    MultidimArray<RFLOAT> &m1,
    MultidimArray<RFLOAT> &m2
);

std::pair<MultidimArray<RFLOAT>, MultidimArray<RFLOAT>> getAmplitudeCorrelationAndDifferentialPhaseResidual(
    const MultidimArray<Complex> &FT1, const MultidimArray<Complex> &FT2
);

std::pair<MultidimArray<RFLOAT>, MultidimArray<RFLOAT>> getAmplitudeCorrelationAndDifferentialPhaseResidual(
    MultidimArray<RFLOAT> &m1, MultidimArray<RFLOAT> &m2
);

std::vector<RFLOAT> cosDeltaPhase(const MultidimArray<Complex> &FT1, const MultidimArray<Complex> &FT2);

// Get precalculated AB-matrices for on-the-fly shift calculations (without tabulated sine and cosine)
void getAbMatricesForShiftImageInFourierTransform(MultidimArray<Complex> &in, MultidimArray<Complex> &out,
                                                  RFLOAT oridim, RFLOAT shift_x, RFLOAT shift_y, RFLOAT shift_z = 0.);

MultidimArray<Complex> shiftImageInFourierTransform(
    const MultidimArray<Complex> &in,
    RFLOAT oridim, RFLOAT xshift, RFLOAT yshift, RFLOAT zshift = 0.0
);

void shiftImageInFourierTransform(
    MultidimArray<Complex> &in_out,
    RFLOAT oridim, RFLOAT xshift, RFLOAT yshift, RFLOAT zshift = 0.0
);

// Shift an image through phase-shifts in its Fourier Transform (without tabulated sine and cosine)
// Note that in and out may be the same array, in that case in is overwritten with the result
// if oridim is in pixels, xshift, yshift and zshift should be in pixels as well!
// or both can be in Angstroms
void shiftImageInFourierTransform(
    const MultidimArray<Complex> &in, MultidimArray<Complex> &out,
    RFLOAT oridim, RFLOAT shift_x, RFLOAT shift_y, RFLOAT shift_z = 0.0
);

void shiftImageInFourierTransformWithTabSincos(MultidimArray<Complex> &in, MultidimArray<Complex> &out,
                                               RFLOAT oridim, long int newdim,
                                               TabSine& tabsin, TabCosine& tabcos,
                                               RFLOAT xshift, RFLOAT yshift, RFLOAT zshift = 0.);

static RFLOAT amplitude(Complex x) { return abs(x); }

static RFLOAT power(Complex x) { return norm(x); }

#define AMPLITUDE_MAP 0
#define PHASE_MAP 1

/** Get the amplitude or power spectrum of the map in Fourier space.
 * @ingroup FourierOperations
    i.e. the radial average of the (squared) amplitudes of all Fourier components
*/
MultidimArray<RFLOAT> getSpectrum(const MultidimArray<RFLOAT> &Min,
                                  RFLOAT(*spectrum_type)(Complex) = power);

/** Divide the input map in Fourier-space by the spectrum provided.
 *  @ingroup FourierOperations
*/
void divideBySpectrum(MultidimArray<RFLOAT> &Min, const MultidimArray<RFLOAT> &spectrum);

/** Multiply the input map in Fourier-space by the spectrum provided.
 *  @ingroup FourierOperations
*/
void multiplyBySpectrum(MultidimArray<RFLOAT> &Min, const MultidimArray<RFLOAT> &spectrum);

/** Perform a whitening of the amplitude/power_class spectrum of a 3D map
 *  @ingroup FourierOperations
    If leave_origin_intact==true, the origin pixel will remain untouched
*/
MultidimArray<RFLOAT> whitenSpectrum(
    const MultidimArray<RFLOAT> &Min,
    RFLOAT(*spectrum_type)(Complex) = amplitude, bool leave_origin_intact=false
);

/** Adapts Min to have the same spectrum as spectrum_ref
 *  @ingroup FourierOperations
    If only_amplitudes==true, the amplitude rather than the power_class spectrum will be equalized
*/
MultidimArray<RFLOAT> adaptSpectrum(
    const MultidimArray<RFLOAT> &Min, const MultidimArray<RFLOAT> &spectrum_ref,
    RFLOAT(*spectrum_type)(Complex) = amplitude, bool leave_origin_intact=false
);

/** Kullback-Leibler divergence */
RFLOAT getKullbackLeiblerDivergence(MultidimArray<Complex> &Fimg,
                                    MultidimArray<Complex> &Fref, MultidimArray<RFLOAT> &sigma2,
                                    MultidimArray<RFLOAT> &p_i, MultidimArray<RFLOAT> &q_i,
                                    int highshell = -1, int lowshell = -1);


// Resize a map by windowing its Fourier Transform
void resizeMap(MultidimArray<RFLOAT> &img, int newsize);

// Apply a B-factor to a map in Fourier space
void applyBFactorToMap(MultidimArray<Complex> &FT, int ori_size, RFLOAT bfactor, RFLOAT angpix);

// Apply a B-factor to a map in real space
void applyBFactorToMap(MultidimArray<RFLOAT> &img, RFLOAT bfactor, RFLOAT angpix);

// Apply a Laplacian-of-Gaussian filter to a map in Fourier space
void LoGFilterMap(MultidimArray<Complex> &FT, int ori_size, RFLOAT sigma, RFLOAT angpix);

// Apply a Laplacian-of-Gaussian filter to a map in real space
void LoGFilterMap(MultidimArray<RFLOAT> &img, RFLOAT sigma, RFLOAT angpix);

// Low-pass filter a map in Fourier space
void lowPassFilterMap(MultidimArray<Complex> &FT, int size,
                      RFLOAT low_pass, RFLOAT angpix, int filter_edge_width = 2);

void highPassFilterMap(MultidimArray<Complex> &FT, int size,
                      RFLOAT low_pass, RFLOAT angpix, int filter_edge_width = 2);

// Low-pass and high-pass filter a map in real space
void lowPassFilterMap(MultidimArray<RFLOAT> &img, RFLOAT low_pass, RFLOAT angpix, int filter_edge_width = 2);
void highPassFilterMap(MultidimArray<RFLOAT> &img, RFLOAT low_pass, RFLOAT angpix, int filter_edge_width = 2);

// Directional filter a map in Fourier space
void directionalFilterMap(MultidimArray<Complex> &FT, int ori_size,
                          RFLOAT low_pass, RFLOAT angpix, int axis = 0, int filter_edge_width = 2);
void directionalFilterMap(MultidimArray<RFLOAT> &img, RFLOAT low_pass, RFLOAT angpix, int axis = 0, int filter_edge_width = 2);

/*
 *	Beamtilt x and y are given in mradians
 *	Wavelength in Angstrom, Cs in mm
 *	Phase shifts caused by the beamtilt will be calculated and applied to Fimg
 */
void selfApplyBeamTilt(MultidimArray<Complex> &Fimg, RFLOAT beamtilt_x, RFLOAT beamtilt_y,
                       RFLOAT wavelength, RFLOAT Cs, RFLOAT angpix, int ori_size);

/* same as above, but for the anisotropic coma model*/
void selfApplyBeamTilt(MultidimArray<Complex> &Fimg,
                       RFLOAT beamtilt_x, RFLOAT beamtilt_y,
                       RFLOAT beamtilt_xx, RFLOAT beamtilt_xy, RFLOAT beamtilt_yy,
                       RFLOAT wavelength, RFLOAT Cs, RFLOAT angpix, int ori_size);

void applyBeamTilt(const MultidimArray<Complex> &Fin, MultidimArray<Complex> &Fout, RFLOAT beamtilt_x, RFLOAT beamtilt_y,
                   RFLOAT wavelength, RFLOAT Cs, RFLOAT angpix, int ori_size);

MultidimArray<RFLOAT> padAndFloat2DMap(const MultidimArray<RFLOAT> &v, int factor = 2);

MultidimArray<RFLOAT> amplitudeOrPhaseMap(const MultidimArray<RFLOAT> &v, int output_map_type);

void helicalLayerLineProfile(const MultidimArray<RFLOAT> &v, std::string title, std::string fn_eps);

MultidimArray<RFLOAT> generateBinaryHelicalFourierMask(
    long int xdim, long int zdim, long int ydim,
    std::vector<RFLOAT> exclude_begin, std::vector<RFLOAT> exclude_end, RFLOAT angpix
);

template <class T>
void cropInFourierSpace(MultidimArray<T> &Fref, MultidimArray<T> &Fbinned) {
    const int nfx = Xsize(Fref), nfy = Ysize(Fref);
    const int new_nfx = Xsize(Fbinned), new_nfy = Ysize(Fbinned);
    const int half_new_nfy = new_nfy / 2;

    if (new_nfx > nfx || new_nfy > nfy) REPORT_ERROR("Invalid size given to cropInFourierSpace");

    for (int y = 0; y < half_new_nfy; y++) {
        for (int x = 0; x < new_nfx; x++) {
            direct::elem(Fbinned, x, y) =  direct::elem(Fref, x, y);
        }
    }
    for (int y = half_new_nfy; y < new_nfy; y++) {
        for (int x = 0; x < new_nfx; x++) {
            direct::elem(Fbinned, x, y) =  direct::elem(Fref, x, nfy - new_nfy + y);
        }
    }
}

static MultidimArray<RFLOAT> real(const MultidimArray<Complex>& arr) {
    MultidimArray<RFLOAT> out (arr.xdim, arr.ydim, arr.zdim, arr.ndim);
    for (long int n = 0; n < out.size(); n++) {
        out[n] = arr[n].real;
    }
    return out;
}

static MultidimArray<RFLOAT> imag(const MultidimArray<Complex>& arr) {
    MultidimArray<RFLOAT> out (arr.xdim, arr.ydim, arr.zdim, arr.ndim);
    for (long int n = 0; n < out.size(); n++) {
        out[n] = arr[n].imag;
    }
    return out;
}

struct pthread_lock_guard {

    pthread_mutex_t *mutex;

    pthread_lock_guard(pthread_mutex_t *mutex): mutex(mutex) {
        pthread_mutex_lock(mutex);
    }

    ~pthread_lock_guard() {
        pthread_mutex_unlock(mutex);
    }

};

#endif // __RELIONFFTW_H
