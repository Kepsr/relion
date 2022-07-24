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
 * This macro is used to generate loops for the volume in an easy way. It
 * defines internal indexes 'k','i' and 'j' which ranges the volume using its
 * physical definition. It also defines 'kp', 'ip' and 'jp', which are the logical coordinates
 * It also works for 1D or 2D FFTW transforms
 *
 * @code
 * FOR_ALL_ELEMENTS_IN_FFTW_TRANSFORM(V) {
 *	int r2 = jp*jp + ip*ip + kp*kp;
 *
 *	std::cout << "element at physical coords: "<< i<<" "<<j<<" "<<k<<" has value: " << direct::elem(m, i, j, k) << std::endl;
 *	std::cout << "its logical coords are: "<< ip<<" "<<jp<<" "<<kp<<std::endl;
 *	std::cout << "its distance from the origin = "<<sqrt(r2)<<std::endl;
 *
 * }
 * @endcode
 */
#define FOR_ALL_ELEMENTS_IN_FFTW_TRANSFORM(V) \
    for (long int k = 0, kp = 0; k < Zsize(V); k++, kp = k < Xsize(V) ? k : k - Zsize(V)) \
    for (long int j = 0, jp = 0; j < Ysize(V); j++, jp = j < Xsize(V) ? j : j - Ysize(V)) \
    for (long int i = 0, ip = 0; i < Xsize(V); i++, ip = i)

/** For all direct elements in the complex array in FFTW format.
 * The same as above, but now only for 2D images (this saves some time as k is not sampled
 */
// FOR_i_j_ip_jp_IN_FFTW_TRANSFORM2D
#define FOR_ALL_ELEMENTS_IN_FFTW_TRANSFORM2D(V) \
    for (long int j = 0, jp = 0; j < Ysize(V); j++, jp = j < Xsize(V) ? j : j - Ysize(V)) \
    for (long int i = 0, ip = 0; i < Xsize(V); i++, ip = i)

namespace FFTW {

    /** Logical access to FFTW element in 3D
     *
     * @code
     * FFTW::elem(V, -1, -2, 1) = 1;
     * val = FFTW::elem(V, -1, -2, 1);
     * @endcode
     */
    template <typename T>
    T& elem(const MultidimArray<T> &V, int ip, int jp, int kp) {
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
    T& elem(const MultidimArray<T> &V, int ip, int jp) {
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
 * MultidimArray<Complex> Vfft;
 * transformer.FourierTransform(V(),Vfft,false);
 * MultidimArray<RFLOAT> Vmag;
 * Vmag.resize(Vfft);
 * FOR_ALL_ELEMENTS_IN_ARRAY3D(Vmag)
 *	Vmag(k,i,j)=20*log10(abs(Vfft(k,i,j)));
 * @endcode
 */
class FourierTransformer {

    public:

    // Pointer to a real array
    MultidimArray<RFLOAT> *fReal;

    // Pointer to a complex array
    MultidimArray<Complex> *fComplex;

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
        If getCopy is false, an alias to the transformed data is returned.
        This is a faster option since a copy of all the data is avoided,
        but you need to be careful that an inverse Fourier transform may
        change the data.
    */
    template <typename T>
    void FourierTransform(MultidimArray<T> &v, MultidimArray<tComplex<T> > &V, bool getCopy=true, bool force_new_plans = false) {
        setReal(v, force_new_plans);
        Transform(FFTW_FORWARD);
        if (getCopy) {
            getFourierCopy(V);
        } else {
            getFourierAlias(V);
        }
    }

    template <typename T>
    MultidimArray<tComplex<T> > FourierTransform(MultidimArray<T> &v, bool force_new_plans = false) {
        MultidimArray<tComplex<T> > V;
        setReal(v, force_new_plans);
        Transform(FFTW_FORWARD);
        getFourierCopy(V);
        return V;
    }

    /** Compute the Fourier transform.
        The data is taken from the matrix with which the object was
        created. */
    void FourierTransform();

    /** Inforce Hermitian symmetry.
        If the Fourier transform risks of losing Hermitian symmetry,
        use this function to renforce it. */
    void enforceHermitianSymmetry();

    /** Compute the inverse Fourier transform.
        The result is stored in the same real data that was passed for
        the forward transform. The Fourier coefficients are taken from
        the internal Fourier coefficients */
    void inverseFourierTransform();

    /** Compute the inverse Fourier transform.
        New data is provided for the Fourier coefficients and the output
        can be any matrix1D, 2D or 3D. It is important that the output
        matrix is already resized to the right size before entering
        in this function. */
    template <typename T>
    void inverseFourierTransform(const MultidimArray<tComplex<T> > &V, MultidimArray<T> &v) {
        setReal(v);
        setFourier(V);
        Transform(FFTW_BACKWARD);
    }

    template <typename T>
    MultidimArray<tComplex<T> > inverseFourierTransform(const MultidimArray<tComplex<T> > &V) {
        MultidimArray<T> v;
        setReal(v);
        setFourier(V);
        Transform(FFTW_BACKWARD);
        return v;
    }

    /** Get Fourier coefficients. */
    template <typename T>
    void getFourierAlias(T& V) { V.alias(fFourier); }

    /** Get Fourier coefficients. */
    MultidimArray<Complex>& getFourierReference() { return fFourier; }

    /** Get Fourier coefficients. */
    template <typename T>
    void getFourierCopy(T& V) {
        V.reshape(fFourier);
        memcpy(
            V.data, fFourier.data,
            fFourier.size() * 2 * sizeof(RFLOAT)
        );
    }

    /** Return a complete Fourier transform (two halves).
    */
    template <typename T>
    void getCompleteFourier(T& V) {
        V.reshape(*fReal);
        int ndim = 3;
        if (Zsize(*fReal) == 1) {
            ndim = 2;
            if (Ysize(*fReal) == 1)
                ndim = 1;
        }
        switch (ndim) {
            case 1:
            for (long int i = 0; i < Xsize(V); i++) {
                direct::elem(V, i) = i < Xsize(fFourier) ? direct::elem(fFourier, i) :
                    conj(direct::elem(fFourier, Xsize(*fReal) - i));
            } break;
            case 2:
                for (long int j = 0; j < Ysize(V); j++)
                for (long int i = 0; i < Xsize(V); i++) {
                    direct::elem(V, i, j) = j < Xsize(fFourier) ? direct::elem(fFourier, i, j) :
                        conj(direct::elem(fFourier, Xsize(*fReal) - j, (Ysize(*fReal) - i) % Ysize(*fReal)));
                } break;
            case 3:
                for (long int k = 0; k < Zsize(V); k++)
                for (long int j = 0; j < Ysize(V); j++)
                for (long int i = 0; i < Xsize(V); i++) {
                    direct::elem(V, i, j, k) = j < Xsize(fFourier) ? direct::elem(fFourier, i, j, k) :
                        conj(direct::elem(fFourier, Xsize(*fReal) - j, (Ysize(*fReal) - i) % Ysize(*fReal), (Zsize(*fReal) - k) % Zsize(*fReal)));
                } break;
        }
    }

    /** Set one half of the FT in fFourier from the input complete Fourier transform (two halves).
        The fReal and fFourier already should have the right sizes
    */
    template <typename T>
    void setFromCompleteFourier(T& V) {
        int ndim = 3;
        if (Zsize(*fReal) == 1) {
            ndim = 2;
            if (Ysize(*fReal) == 1)
                ndim = 1;
        }
        switch (ndim) {
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

    public:

    /* Pointer to the array of RFLOATs with which the plan was computed */
    RFLOAT *dataPtr;

    /* Pointer to the array of complex<RFLOAT> with which the plan was computed */
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
    const MultidimArray<RFLOAT> &getReal() const;
    const MultidimArray<Complex> &getComplex() const;

    /** Set a Multidimarray for input.
        The data of img will be the one of fReal. In forward
        transforms it is not modified, but in backward transforms,
        the result will be stored in img. This means that the size
        of img cannot change between calls. */
    void setReal(MultidimArray<RFLOAT> &img, bool force_new_plans = false);

    /** Set a Multidimarray for input.
        The data of img will be the one of fComplex. In forward
        transforms it is not modified, but in backward transforms,
        the result will be stored in img. This means that the size
        of img cannot change between calls. */
    void setReal(MultidimArray<Complex> &img, bool force_new_plans = false);

    /** Set a Multidimarray for the Fourier transform.
        The values of the input array are copied in the internal array.
        It is assumed that the container for the real image as well as
        the one for the Fourier array are already resized.
        No plan is updated. */
    void setFourier(const MultidimArray<Complex> &imgFourier);
};

// Randomize phases beyond the given F-space shell (index) of R-space input image
void randomizePhasesBeyond(MultidimArray<RFLOAT> &I, int index);

// Randomize phases beyond the given F-space shell (index) of F-space input image
// void randomizePhasesBeyond(MultidimArray<Complex> &v, int index);

template <typename T>
void CenterFFTbySign(MultidimArray <T> &v) {
    // This technique does not work when the sizes of dimensions of iFFT(v) are odd.
    // Unfortunately, this cannot be checked within this function...
    // Forward and backward shifts are equivalent.

    FOR_ALL_ELEMENTS_IN_ARRAY3D(v) {
    // NOTE: != has higher precedence than & in C as pointed out in GitHub issue #637.
    // So (k ^ i ^ j) & 1 != 0 is not good (fortunately in this case the behaviour happened to be the same)
        if (((i ^ j ^ k) & 1) != 0) // if ODD
            direct::elem(v, i, j, k) *= -1;
    }
}


/** Center an array, to have its origin at the origin of the FFTW
 *
 */
template <typename T>
void CenterFFT(MultidimArray<T>& v, bool forward) {
    #ifndef FAST_CENTERFFT
    if (v.getDim() == 1) {
        // 1D
        MultidimArray<T> aux;
        int l, shift;

        l = Xsize(v);
        aux.reshape(l);
        shift = l / 2;

        if (!forward) { shift = -shift; }

        // Shift the input in an auxiliary vector
        for (int i = 0; i < l; i++) {
            int ip = i + shift;

                 if (ip <  0) { ip += l; }
            else if (ip >= l) { ip -= l; }

            aux(ip) = direct::elem(v, i);
        }

        // Copy the vector
        for (int i = 0; i < l; i++)
            direct::elem(v, i) = direct::elem(aux, i);
    } else if (v.getDim() == 2) {
        // 2D

        // Shift in the X direction
        int l = Xsize(v);
        MultidimArray<T> aux;
        aux.reshape(l);
        int shift = l / 2;

        if (!forward) { shift = -shift; }

        for (int j = 0; j < Ysize(v); j++) {
            // Shift the input in an auxiliar vector
            for (int i = 0; i < l; i++) {
                int ip = i + shift;

                     if (ip <  0) { ip += l; }
                else if (ip >= l) { ip -= l; }

                aux(ip) = direct::elem(v, i, j);
            }

            // Copy the vector
            for (int i = 0; i < l; i++)
                direct::elem(v, i, j) = direct::elem(aux, i);
        }

        // Shift in the Y direction
        l = Ysize(v);
        aux.reshape(l);
        shift = l / 2;

        if (!forward) { shift = -shift; }

        for (int i = 0; i < Xsize(v); i++) {
            // Shift the input in an auxiliar vector
            for (int j = 0; j < j; i++) {
                int jp = j + shift;

                     if (jp <  0) { jp += l; }
                else if (jp >= l) { jp -= l; }

                aux(jp) = direct::elem(v, i, j);
            }

            // Copy the vector
            for (int j = 0; j < l; j++)
                direct::elem(v, i, j) = direct::elem(aux, j);
        }
    } else if (v.getDim() == 3) {
        // 3D
        MultidimArray<T> aux;
        int l, shift;

        // Shift in the X direction
        l = Xsize(v);
        aux.reshape(l);
        shift = l / 2;

        if (!forward) { shift = -shift; }

        for (int k = 0; k < Zsize(v); k++)
        for (int j = 0; j < Ysize(v); j++) {
            // Shift the input in an auxiliar vector
            for (int i = 0; i < l; i++) {
                int ip = i + shift;

                     if (ip <  0) { ip += l; }
                else if (ip >= l) { ip -= l; }

                aux(ip) = direct::elem(v, i, j, k);
            }

            // Copy the vector
            for (int i = 0; i < l; i++)
                direct::elem(v, i, j, k) = direct::elem(aux, i);
        }

        // Shift in the Y direction
        l = Ysize(v);
        aux.reshape(l);
        shift = l / 2;

        if (!forward) { shift = -shift; }

        for (int k = 0; k < Zsize(v); k++)
        for (int i = 0; i < Xsize(v); i++) {
            // Shift the input in an auxiliar vector
            for (int j = 0; j < l; j++) {
                int jp = j + shift;

                     if (jp <  0) { jp += l; }
                else if (jp >= l) { jp -= l; }

                aux(jp) = direct::elem(v, i, j, k);
            }

            // Copy the vector
            for (int j = 0; j < l; j++)
                direct::elem(v, i, j, k) = direct::elem(aux, j);
        }

        // Shift in the Z direction
        l = Zsize(v);
        aux.reshape(l);
        shift = l / 2;

        if (!forward) { shift = -shift; }

        for (int j = 0; j < Ysize(v); j++)
        for (int i = 0; i < Xsize(v); i++) {
            // Shift the input in an auxiliar vector
            for (int k = 0; k < l; k++) {
                int kp = k + shift;

                     if (kp <  0) { kp += l; }
                else if (kp >= l) { kp -= l; }

                aux(kp) = direct::elem(v, i, j, k);
            }

            // Copy the vector
            for (int k = 0; k < l; k++)
                direct::elem(v, i, j, k) = direct::elem(aux, k);
        }
    } else {
        v.printShape();
        REPORT_ERROR("CenterFFT ERROR: Dimension should be 1, 2 or 3");
    }
    #else // FAST_CENTERFFT
    if (v.getDim() == 1) {
        // 1D

        int l = Xsize(v);
        MultidimArray<T> aux;
        aux.reshape(l);
        int shift = l / 2;

        if (!forward) { shift = -shift; }

        // Shift the input in an auxiliary vector
        for (int i = 0; i < l; i++) {
            int ip = i + shift;

                 if (ip <  0) { ip += l; }
            else if (ip >= l) { ip -= l; }

            aux(ip) = direct::elem(v, i);
        }

        // Copy the vector
        for (int i = 0; i < l; i++)
            direct::elem(v, i) = direct::elem(aux, i);
    } else if (v.getDim() == 2 ) {
        int batchSize = 1;
        int xSize = Xsize(v);
        int ySize = Ysize(v);

        int xshift = xSize / 2;
        int yshift = ySize / 2;

        if (!forward) {
            xshift = -xshift;
            yshift = -yshift;
        }

        size_t image_size = xSize * ySize;
        size_t isize2 = image_size / 2;
        int blocks = ceilf((float) (image_size / (float) (2 * CFTT_BLOCK_SIZE)));

		// for (int i = 0; i < blocks; i++) {
        tbb::parallel_for(0, blocks, [&](int i) {
            size_t pixel_start = i * CFTT_BLOCK_SIZE;
            size_t pixel_end = (i + 1) * CFTT_BLOCK_SIZE;
            if (pixel_end > isize2) { pixel_end = isize2; }

            CpuKernels::centerFFT_2D<T>(
                batchSize, pixel_start, pixel_end, v.data,
                (size_t) xSize * ySize, xSize, ySize, xshift, yshift
            );
        }
        );
    } else if (v.getDim() == 3) {
        int  batchSize = 1;
        int xSize = Xsize(v);
        int ySize = Ysize(v);
        int zSize = Zsize(v);

        if(zSize > 1) {
            int xshift = xSize / 2;
            int yshift = ySize / 2;
            int zshift = zSize / 2;

            if (!forward) {
                xshift = -xshift;
                yshift = -yshift;
                zshift = -zshift;
            }

            size_t image_size = xSize * ySize * zSize;
            size_t isize2 = image_size / 2;
            int block =ceilf((float) (image_size / (float) (2 * CFTT_BLOCK_SIZE)));
			// for (int i = 0; i < block; i++){
            tbb::parallel_for(0, block, [&](int i) {
                size_t pixel_start = i * CFTT_BLOCK_SIZE;
                size_t pixel_end = (i + 1) * CFTT_BLOCK_SIZE;
                if (pixel_end > isize2) { pixel_end = isize2; }

                CpuKernels::centerFFT_3D<T>(
                    batchSize, pixel_start, pixel_end, v.data,
                    (size_t) xSize * ySize * zSize, xSize, ySize, zSize, xshift, yshift, zshift
                );
            }
            );
        } else {
            int xshift = xSize / 2;
            int yshift = ySize / 2;

            if (!forward) {
                xshift = -xshift;
                yshift = -yshift;
            }

            size_t image_size = xSize * ySize;
            size_t isize2 = image_size / 2;
            int blocks = ceilf((float) (image_size / (float) (2 * CFTT_BLOCK_SIZE)));
			// for (int i = 0; i < blocks; i++) {
            tbb::parallel_for(0, blocks, [&](int i) {
                size_t pixel_start = i * CFTT_BLOCK_SIZE;
                size_t pixel_end = (i + 1) * CFTT_BLOCK_SIZE;
                if (pixel_end > isize2) { pixel_end = isize2; }

                CpuKernels::centerFFT_2D<T>(
                    batchSize, pixel_start, pixel_end, v.data,
                    (size_t) xSize * ySize, xSize, ySize, xshift, yshift
                );
            }
            );
        }
    } else {
        v.printShape();
        REPORT_ERROR("CenterFFT ERROR: Dimension should be 1, 2 or 3");
    }
    #endif	// FAST_CENTERFFT
}



// Window an FFTW-centered Fourier-transform to a given size
template<class T>
void windowFourierTransform(MultidimArray<T> &in, MultidimArray<T> &out, long int newdim) {
    // Check size of the input array
    if (Ysize(in) > 1 && Ysize(in) / 2 + 1 != Xsize(in))
        REPORT_ERROR("windowFourierTransform ERROR: the Fourier transform should be of an image with equal sizes in all dimensions!");
    long int newhdim = newdim / 2 + 1;

    // If same size, just return input
    // Sjors 5dec2017: only check for xdim is not enough, even/off ydim leaves ambiguity for dim>1
    if (newdim == Ysize(in) && newhdim == Xsize(in)) {
        out = in;
        return;
    }

    // Otherwise apply a windowing operation
    // Initialise output array
    switch (in.getDim()) {

        case 1:
        out.initZeros(newhdim);
        break;

        case 2:
        out.initZeros(newdim, newhdim);
        break;

        case 3:
        out.initZeros(newdim, newdim, newhdim);
        break;

        default:
        REPORT_ERROR("windowFourierTransform ERROR: dimension should be 1, 2 or 3!");

    }

    if (newhdim > Xsize(in)) {
        long int max_r2 = (Xsize(in) - 1) * (Xsize(in) - 1);
        FOR_ALL_ELEMENTS_IN_FFTW_TRANSFORM(in) {
            // Make sure windowed FT has nothing in the corners, otherwise we end up with an asymmetric FT!
            if (kp * kp + ip * ip + jp * jp <= max_r2)
                FFTW::elem(out, ip, jp, kp) = FFTW::elem(in, ip, jp, kp);
        }
    } else {
        FOR_ALL_ELEMENTS_IN_FFTW_TRANSFORM(out) {
            FFTW::elem(out, ip, jp, kp) = FFTW::elem(in, ip, jp, kp);
        }
    }
}

// Same as above, acts on the input array directly
template<class T>
void windowFourierTransform(MultidimArray<T> &V, long int newdim) {
    // Check size of the input array
    if (Ysize(V) > 1 && Ysize(V)/2 + 1 != Xsize(V))
        REPORT_ERROR("windowFourierTransform ERROR: the Fourier transform should be of an image with equal sizes in all dimensions!");
    long int newhdim = newdim / 2 + 1;

    // If same size, just return input
    // Sjors 5dec2017: only check for xdim is not enough, even/off ydim leaves ambiguity for dim>1
    if (newdim == Ysize(V) && newhdim == Xsize(V)) return;

    MultidimArray<T> tmp;
    windowFourierTransform<T>(V, tmp, newdim);
    V.moveFrom(tmp);
}

// A resize operation in Fourier-space (i.e. changing the sampling of the Fourier Transform) by windowing in real-space
// If recenter=true, the real-space array will be recentered to have its origin at the origin of the FT
template<class T>
void resizeFourierTransform(MultidimArray<T> &in, MultidimArray<T> &out, long int newdim, bool do_recenter=true) {
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
    MultidimArray<Complex> Fin;
    MultidimArray<RFLOAT> Min;
    FourierTransformer transformer;
    long int x0, y0, z0, xF, yF, zF;
    x0 = y0 = z0 = Xmipp::init(newdim);
    xF = yF = zF = Xmipp::last(newdim);

    // Initialise output array
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

    // This is to handle RFLOAT-valued input arrays
    Fin.reshape(Zsize(in), Ysize(in), Xsize(in));
    for (long int i = 0, end = in.size(); i < end; i++) {
        Fin.data[i] = in.data[i];
    }
    transformer.inverseFourierTransform(Fin, Min);
    Min.setXmippOrigin();
    if (do_recenter) {
        CenterFFT(Min, false);
    }

    // Now do the actual windowing in real-space
    Min.window(z0, y0, x0, zF, yF, xF);
    Min.setXmippOrigin();

    // If upsizing: mask the corners to prevent aliasing artefacts
    if (newdim > olddim) {
        FOR_ALL_ELEMENTS_IN_ARRAY3D(Min) {
            if (k * k + i * i + j * j > olddim * olddim / 4) {
                Min.elem(i, j, k) = 0.0;
            }
        }
    }

    // Recenter FFT back again
    if (do_recenter)
        CenterFFT(Min, true);

    // And do the inverse Fourier transform
    transformer.clear();
    transformer.FourierTransform(Min, out);
}

/** Fourier-Ring-Correlation between two multidimArrays using FFT
 * From precalculated Fourier Transforms
 * Simpler I/O than above.
 */
void getFSC(MultidimArray<Complex> &FT1,
            MultidimArray<Complex> &FT2,
            MultidimArray<RFLOAT> &fsc);

/** Fourier-Ring-Correlation between two multidimArrays using FFT
 * @ingroup FourierOperations
 * Simpler I/O than above.
 */
void getFSC(MultidimArray<RFLOAT> & m1,
            MultidimArray<RFLOAT> & m2,
            MultidimArray<RFLOAT> &fsc);

void getAmplitudeCorrelationAndDifferentialPhaseResidual(MultidimArray<Complex> &FT1, MultidimArray<Complex> &FT2,
                                                         MultidimArray<RFLOAT> &acorr, MultidimArray<RFLOAT> &dpr);

void getAmplitudeCorrelationAndDifferentialPhaseResidual(MultidimArray<RFLOAT> &m1, MultidimArray<RFLOAT> &m2,
                                                         MultidimArray<RFLOAT> &acorr, MultidimArray<RFLOAT> &dpr);

MultidimArray<RFLOAT> cosDeltaPhase(MultidimArray<Complex> &FT1, MultidimArray<Complex> &FT2);

// Get precalculated AB-matrices for on-the-fly shift calculations (without tabulated sine and cosine)
void getAbMatricesForShiftImageInFourierTransform(MultidimArray<Complex> &in, MultidimArray<Complex> &out,
                                                  RFLOAT oridim, RFLOAT shift_x, RFLOAT shift_y, RFLOAT shift_z = 0.);

void shiftImageInFourierTransformWithTabSincos(MultidimArray<Complex> &in, MultidimArray<Complex> &out,
                                               RFLOAT oridim, long int newdim,
                                               TabSine& tabsin, TabCosine& tabcos,
                                               RFLOAT xshift, RFLOAT yshift, RFLOAT zshift = 0.);

// Shift an image through phase-shifts in its Fourier Transform (without tabulated sine and cosine)
// Note that in and out may be the same array, in that case in is overwritten with the result
// if oridim is in pixels, xshift, yshift and zshift should be in pixels as well!
// or both can be in Angstroms
void shiftImageInFourierTransform(MultidimArray<Complex> &in, MultidimArray<Complex> &out,
                                  RFLOAT oridim, RFLOAT shift_x, RFLOAT shift_y, RFLOAT shift_z = 0.);

#define POWER_SPECTRUM 0
#define AMPLITUDE_SPECTRUM 1
#define AMPLITUDE_MAP 0
#define PHASE_MAP 1

/** Get the amplitude or power_class spectrum of the map in Fourier space.
 * @ingroup FourierOperations
    i.e. the radial average of the (squared) amplitudes of all Fourier components
*/
void getSpectrum(MultidimArray<RFLOAT> &Min,
                 MultidimArray<RFLOAT> &spectrum,
                 int spectrum_type=POWER_SPECTRUM);

/** Divide the input map in Fourier-space by the spectrum provided.
 *  @ingroup FourierOperations
    If leave_origin_intact==true, the origin pixel will remain untouched
*/
void divideBySpectrum(MultidimArray<RFLOAT> &Min, MultidimArray<RFLOAT> &spectrum, bool leave_origin_intact=false);

/** Multiply the input map in Fourier-space by the spectrum provided.
 *  @ingroup FourierOperations
    If leave_origin_intact==true, the origin pixel will remain untouched
*/
void multiplyBySpectrum(MultidimArray<RFLOAT> &Min, MultidimArray<RFLOAT> &spectrum, bool leave_origin_intact=false);

/** Perform a whitening of the amplitude/power_class spectrum of a 3D map
 *  @ingroup FourierOperations
    If leave_origin_intact==true, the origin pixel will remain untouched
*/
void whitenSpectrum(MultidimArray<RFLOAT> &Min, MultidimArray<RFLOAT> &Mout, int spectrum_type=AMPLITUDE_SPECTRUM, bool leave_origin_intact=false);

/** Adapts Min to have the same spectrum as spectrum_ref
 *  @ingroup FourierOperations
    If only_amplitudes==true, the amplitude rather than the power_class spectrum will be equalized
*/
void adaptSpectrum(MultidimArray<RFLOAT> &Min, MultidimArray<RFLOAT> &Mout,
                   const MultidimArray<RFLOAT> &spectrum_ref, int spectrum_type=AMPLITUDE_SPECTRUM, bool leave_origin_intact=false);

/** Kullback-Leibler divergence */
RFLOAT getKullbackLeiblerDivergence(MultidimArray<Complex> &Fimg,
                                    MultidimArray<Complex> &Fref, MultidimArray<RFLOAT> &sigma2,
                                    MultidimArray<RFLOAT> &p_i, MultidimArray<RFLOAT> &q_i,
                                    int highshell = -1, int lowshell = -1);


// Resize a map by windowing it's Fourier Transform
void resizeMap(MultidimArray<RFLOAT> &img, int newsize);

// Apply a B-factor to a map (given it's Fourier transform)
void applyBFactorToMap(MultidimArray<Complex> &FT, int ori_size, RFLOAT bfactor, RFLOAT angpix);

// Apply a B-factor to a map (given it's real-space array)
void applyBFactorToMap(MultidimArray<RFLOAT> &img, RFLOAT bfactor, RFLOAT angpix);

// Apply a Laplacian-of-Gaussian filter to a map (given it's Fourier transform)
void LoGFilterMap(MultidimArray<Complex> &FT, int ori_size, RFLOAT sigma, RFLOAT angpix);

// Apply a Laplacian-of-Gaussian filter to a map (given it's real-space array)
void LoGFilterMap(MultidimArray<RFLOAT> &img, RFLOAT sigma, RFLOAT angpix);

// Low-pass filter a map (given it's Fourier transform)
void lowPassFilterMap(MultidimArray<Complex> &FT, int ori_size,
                      RFLOAT low_pass, RFLOAT angpix, int filter_edge_width = 2, bool do_highpass_instead = false);

// Low-pass and high-pass filter a map (given it's real-space array)
void lowPassFilterMap(MultidimArray<RFLOAT> &img, RFLOAT low_pass, RFLOAT angpix, int filter_edge_width = 2);
void highPassFilterMap(MultidimArray<RFLOAT> &img, RFLOAT low_pass, RFLOAT angpix, int filter_edge_width = 2);

// Directional filter a map (given it's Fourier transform)
void directionalFilterMap(MultidimArray<Complex> &FT, int ori_size,
                          RFLOAT low_pass, RFLOAT angpix, std::string axis = "x", int filter_edge_width = 2);
void directionalFilterMap(MultidimArray<RFLOAT> &img, RFLOAT low_pass, RFLOAT angpix, std::string axis = "x", int filter_edge_width = 2);

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

void padAndFloat2DMap(const MultidimArray<RFLOAT> &v, MultidimArray<RFLOAT> &out, int factor = 2);

void amplitudeOrPhaseMap(const MultidimArray<RFLOAT> &v, MultidimArray<RFLOAT> &amp, int output_map_type);

void helicalLayerLineProfile(const MultidimArray<RFLOAT> &v, std::string title, std::string fn_eps);

void generateBinaryHelicalFourierMask(MultidimArray<RFLOAT> &mask, std::vector<RFLOAT> exclude_begin, std::vector<RFLOAT> exclude_end, RFLOAT angpix);

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

#endif // __RELIONFFTW_H
