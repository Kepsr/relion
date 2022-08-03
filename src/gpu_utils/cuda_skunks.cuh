#ifndef CUDA_SKUNKS_CUH_
#define CUDA_SKUNKS_CUH_

#include "src/projector.h"
#include "src/multidim_array.h"
#include "src/fftw.h"

void computeFourierTransformMap(Projector *P, MultidimArray<RFLOAT> &vol_in, MultidimArray<RFLOAT> &power_spectrum, int current_size = -1, int nr_threads = 1, bool do_gridding = true) {
    MultidimArray<RFLOAT> Mpad;
    FourierTransformer transformer;
    RFLOAT normfft;

    // Size of padded real-space volume
    int padoridim = P->padding_factor * P->ori_size;

    // Initialize data array of the oversampled transform
    P->ref_dim = vol_in.getDim();

    // Make Mpad
    switch (P->ref_dim) {

        case 2:
        Mpad.initZeros(padoridim, padoridim);
        normfft = P->padding_factor * P->padding_factor;
        break;

        case 3:
        Mpad.initZeros(padoridim, padoridim, padoridim);
        normfft = P->data_dim == 3 ?
            P->padding_factor * P->padding_factor * P->padding_factor :
            P->padding_factor * P->padding_factor * P->padding_factor * P->ori_size;
        break;

        default:
        REPORT_ERROR("Projector::computeFourierTransformMap%%ERROR: Dimension of the data array should be 2 or 3");

    }

    // First do a gridding pre-correction on the real-space map:
    // Divide by the inverse Fourier transform of the interpolator in Fourier-space
    // 10feb11: at least in 2D case, this seems to be the wrong thing to do!!!
    // TODO: check what is best for subtomo!
    if (do_gridding)// && data_dim != 3)
        P->griddingCorrect(vol_in);

    // Pad translated map with zeros
    vol_in.setXmippOrigin();
    Mpad.setXmippOrigin();
    FOR_ALL_ELEMENTS_IN_ARRAY3D(vol_in) // This will also work for 2D
        Mpad.elem(i, j, k) = vol_in.elem(i, j, k);

    // Translate padded map to put origin of FT in the center
    CenterFFT(Mpad, true);

    // Calculate the oversampled Fourier transform
    MultidimArray<Complex> &Faux = transformer.FourierTransform(Mpad);

    // Free memory: Mpad no longer needed
    Mpad.clear();

    // Resize data array to the right size and initialise to zero
    P->initZeros(current_size);

    // Fill data only for those points with distance to origin less than max_r
    // (other points will be zero because of initZeros() call above
    // Also calculate radial power spectrum
    power_spectrum.initZeros(P->ori_size / 2 + 1);
    MultidimArray<RFLOAT> counter = MultidimArray<RFLOAT>::zeros(power_spectrum);

    int max_r2 = P->r_max * P->r_max * P->padding_factor * P->padding_factor;
    FOR_ALL_ELEMENTS_IN_FFTW_TRANSFORM(Faux) {
        // This will also work for 2D
        int r2 = kp * kp + ip * ip + jp * jp;
        // The Fourier Transforms are all "normalised" for 2D transforms of size = ori_size x ori_size
        if (r2 <= max_r2) {
            // Set data array
            P->data.elem(ip, jp, kp) = direct::elem(Faux, i, j, k) * normfft;

            // Calculate power spectrum
            int ires = round(sqrt((RFLOAT) r2) / P->padding_factor);
            // Factor two because of two-dimensionality of the complex plane
            direct::elem(power_spectrum, ires) += norm(P->data.elem(ip, jp, kp)) / 2.0;
            direct::elem(counter, ires) += 1.0;
        }
    }

    // Calculate radial average of power spectrum
    for (long int i = 0; i < Xsize(power_spectrum); i++) {
        if (direct::elem(counter, i) < 1.0) {
            direct::elem(power_spectrum, i) = 0.0;
        } else {
            direct::elem(power_spectrum, i) /= direct::elem(counter, i);
        }
    }
}

#endif //CUDA_SKUNKS_CUH_
