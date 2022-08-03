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

#ifndef PAR_FFTW_H
#define PAR_FFTW_H

#include <fftw3.h>
#include "src/multidim_array.h"
#include "src/funcs.h"
#include "src/tabfuncs.h"
#include "src/complex.h"
#include "src/fftw.h"
#include "src/CPlot2D.h"

/*
  Parallelizable version of FourierTransformer from src/fftw.h:

  FFTW plans are managed globally and only one can be computed at a time.
  The original class recomputes its plan each time a new Image is passed
  to its FourierTransform() or inverseFourierTransform().
  As a consequence, when working on sets of multiple images, only one
  FT can run at a time.

  This class only recomputes the plans if the size of the image changes.
  The same plan is reused for different images of the same size.

  Otherwise, the two classes are identical.

          -- J. Zivanov, Feb. 9th 2018
*/
class ParFourierTransformer: public FourierTransformer {

    public:

    ParFourierTransformer();

    ~ParFourierTransformer();

    // Copy ctor
    ParFourierTransformer(const ParFourierTransformer& op);

    /** Compute the Fourier transform of a MultidimArray, 2D and 3D.
     */
    template <typename T>
    MultidimArray<tComplex<T> > &FourierTransform(MultidimArray<T> &v) {
        setReal(v);
        Transform(FFTW_FORWARD);
        return fFourier;
    }

    /** Compute the Fourier transform.
        The data is taken from the matrix with which the object was
        created. */
    void FourierTransform();

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
    template <typename T, typename T1>
    void inverseFourierTransform(const T& V, T1& v) {
        setReal(v);
        setFourier(V);
        Transform(FFTW_BACKWARD);
    }

    /** Get Fourier coefficients. */
    MultidimArray<Complex> &getFourier() { return fFourier; }

    /** Return a complete Fourier transform (two halves). */
    template <typename T>
    void getCompleteFourier(MultidimArray<T>& V) const;

    /** Set one half of the FT in fFourier from the input complete Fourier transform (two halves).
        fReal and fFourier should already be of the right size.
    */
    template <typename T>
    void setFromCompleteFourier(const MultidimArray<T>& V);

    public:

    /* Initialise all pointers to NULL */
    void init();

    void Transform(int sign);

    /** Set a Multidimarray for input.
        The data of img will be the one of fReal. In forward
        transforms it is not modified, but in backward transforms,
        the result will be stored in img. This means that the size
        of img cannot change between calls. */
    void setReal(MultidimArray<RFLOAT> &img);

    /** Set a Multidimarray for input.
        The data of img will be the one of fComplex. In forward
        transforms it is not modified, but in backward transforms,
        the result will be stored in img. This means that the size
        of img cannot change between calls. */
    void setReal(MultidimArray<Complex> &img);

    /** Set a Multidimarray for the Fourier transform.
        The values of the input array are copied in the internal array.
        It is assumed that the container for the real image as well as
        the one for the Fourier array are already resized.
        No plan is updated. */
    void setFourier(const MultidimArray<Complex> &imgFourier);

};

#endif
