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
#include "src/projector.h"
#include "src/jaz/gravis/t3Vector.h"
#include <src/time.h>
#ifdef CUDA
#include <cufft.h>
#include <sys/types.h>
#include <sys/ipc.h>
#include <sys/sem.h>
#define SEMKEY 826976737978L /* key value for semget() */
#define PERMS 0666
#endif
// #define DEBUG

// #define PROJ_TIMING
#ifdef PROJ_TIMING
Timer proj_timer;
int TIMING_TOP    = proj_timer.setNew("PROJECTOR - computeFourierTransformMap");
int TIMING_GRID   = proj_timer.setNew("PROJECTOR - gridCorr");
int TIMING_PAD    = proj_timer.setNew("PROJECTOR - padTransMap");
int TIMING_CENTER = proj_timer.setNew("PROJECTOR - centerFFT");
int TIMING_TRANS  = proj_timer.setNew("PROJECTOR - transform");
int TIMING_FAUX   = proj_timer.setNew("PROJECTOR - Faux");
int TIMING_POW    = proj_timer.setNew("PROJECTOR - power_spectrum");
int TIMING_INIT1  = proj_timer.setNew("PROJECTOR - init1");
int TIMING_INIT2  = proj_timer.setNew("PROJECTOR - init2");
#endif

using namespace gravis;


void trilinear(
    RFLOAT &xp, RFLOAT &yp, RFLOAT &zp, const MultidimArray<Complex> &data,
    MultidimArray<Complex> &f2d, int x, int i
) {
    // Only asymmetric half is stored

    bool is_neg_x = xp < 0;

    if (is_neg_x) {
        // Get complex conjugate hermitian symmetry pair
        xp = -xp;
        yp = -yp;
        zp = -zp;
    }

    // Trilinear interpolation (with physical coords)
    // Subtract Yinit and Zinit to accelerate access to data (Xinit=0)
    // In that way use direct::elem, rather than A3D_ELEM
    const int x0 = floor(xp);
    const RFLOAT fx = xp - x0;
    const int x1 = x0 + 1;

    int y0 = floor(yp);
    const RFLOAT fy = yp - y0;
    y0 -= Yinit(data);
    const int y1 = y0 + 1;

    int z0 = floor(zp);
    const RFLOAT fz = zp - z0;
    z0 -= Zinit(data);
    const int z1 = z0 + 1;

    // Avoid reading outside the box
    if (
        x0 < 0 || x0 + 1 >= data.xdim ||
        y0 < 0 || y0 + 1 >= data.ydim ||
        z0 < 0 || z0 + 1 >= data.zdim
    ) return;

    // Matrix access can be accelerated through pre-calculation of z0*xydim etc.
    const Complex d000 = direct::elem(data, x0, y0, z0);
    const Complex d001 = direct::elem(data, x1, y0, z0);
    const Complex d010 = direct::elem(data, x0, y1, z0);
    const Complex d011 = direct::elem(data, x1, y1, z0);
    const Complex d100 = direct::elem(data, x0, y0, z1);
    const Complex d101 = direct::elem(data, x1, y0, z1);
    const Complex d110 = direct::elem(data, x0, y1, z1);
    const Complex d111 = direct::elem(data, x1, y1, z1);

    // Set the interpolated value in the 2D output array
    const Complex dx00 = LIN_INTERP(fx, d000, d001);
    const Complex dx01 = LIN_INTERP(fx, d100, d101);
    const Complex dx10 = LIN_INTERP(fx, d010, d011);
    const Complex dx11 = LIN_INTERP(fx, d110, d111);
    const Complex dxy0 = LIN_INTERP(fy, dx00, dx10);
    const Complex dxy1 = LIN_INTERP(fy, dx01, dx11);

    direct::elem(f2d, i, x) = LIN_INTERP(fz, dxy0, dxy1);

    // Take complex conjugate for half with negative x
    if (is_neg_x) direct::elem(f2d, i, x) = conj(direct::elem(f2d, i, x));
}


void trilinear(
    RFLOAT &xp, RFLOAT &yp, const MultidimArray<Complex> &data,
    MultidimArray<Complex> &f1d, int x
) {
    // Only asymmetric half is stored
    const bool is_neg_x = xp < 0;

    if (is_neg_x) {
        // Get complex conjugate hermitian symmetry pair
        xp = -xp;
        yp = -yp;
    }

    // Trilinear interpolation (with physical coords)
    // Subtract Yinit to accelerate access to data (Xinit=0)
    // In that way use direct::elem, rather than A3D_ELEM
    const int x0 = floor(xp);
    const RFLOAT fx = xp - x0;
    const int x1 = x0 + 1;

    int y0 = floor(yp);
    const RFLOAT fy = yp - y0;
    y0 -= Yinit(data);
    const int y1 = y0 + 1;

    // Matrix access can be accelerated through pre-calculation of z0*xydim etc.
    const Complex d00 = direct::elem(data, x0, y0);
    const Complex d01 = direct::elem(data, x1, y0);
    const Complex d10 = direct::elem(data, x0, y1);
    const Complex d11 = direct::elem(data, x1, y1);

    // Set the interpolated value in the 2D output array
    const Complex dx0 = LIN_INTERP(fx, d00, d01);
    const Complex dx1 = LIN_INTERP(fx, d10, d11);

    direct::elem(f1d, x) = LIN_INTERP(fy, dx0, dx1);

    // Take complex conjugate for half with negative x
    if (is_neg_x) direct::elem(f1d, x) = conj(direct::elem(f1d, x));
}

void nearest_neighbour(
    RFLOAT &xp, RFLOAT &yp, RFLOAT &zp, const MultidimArray<Complex> &data,
    MultidimArray<Complex> &f2d,
    int x, int i
) {
    int x0 = round(xp);
    int y0 = round(yp);
    int z0 = round(zp);

    const bool is_neg_x = x0 < 0;

    if (is_neg_x) {
        // Get complex conjugate hermitian symmetry pair
        x0 = -x0;
        y0 = -y0;
        z0 = -z0;
    }

    const int xr = x0 - Xinit(data);
    const int yr = y0 - Yinit(data);
    const int zr = z0 - Zinit(data);

    if (
        xr < 0 || xr >= data.xdim ||
        yr < 0 || yr >= data.ydim ||
        zr < 0 || zr >= data.zdim
    ) return;

    direct::elem(f2d, i, x) = is_neg_x ?
        conj(direct::elem(data, xr, yr, zr)) : data.elem(xr, yr, zr);
        /// XXX: Should we be mixing direct::elem and MultidimArray::elem?
}

void nearest_neighbour(
    RFLOAT &xp, RFLOAT &yp, const MultidimArray<Complex> &data, 
    MultidimArray<Complex> &f1d, int x
) {
    const int x0 = round(xp);
    const int y0 = round(yp);

    direct::elem(f1d, x) = x0 < 0 ?
        conj(data.elem(-x0, -y0)) : data.elem(x0, y0);
}


void Projector::initialiseData(int current_size) {
    // By default r_max is half ori_size
    r_max = (current_size < 0 ? ori_size : current_size) / 2;

    // Never allow r_max beyond Nyquist...
    r_max = std::min(r_max, ori_size / 2);

    // Set pad_size
    pad_size = 2 * round(padding_factor * r_max) + 3;

    // Short side of data array
    switch (ref_dim) {

        case 2:
        data.resize(pad_size, pad_size / 2 + 1);
        break;

        case 3:
        data.resize(pad_size, pad_size, pad_size / 2 + 1);
        break;

        default:
        REPORT_ERROR("Projector::resizeData%%ERROR: Dimension of the data array should be 2 or 3");

    }

    // Set origin in the y.z-center, but on the left side for x.
    data.setXmippOrigin();
    data.xinit = 0;

}

void Projector::initZeros(int current_size) {
    initialiseData(current_size);
    data.initZeros();
}

long int Projector::getSize() {
    // Short side of data array
    switch (ref_dim) {

        case 2:
        return pad_size * (pad_size / 2 + 1);

        case 3:
        return pad_size * pad_size * (pad_size / 2 + 1);

        default:
        REPORT_ERROR("Projector::resizeData%%ERROR: Dimension of the data array should be 2 or 3");

    }
}


// Fill data array with oversampled Fourier transform, and calculate its power spectrum
void Projector::computeFourierTransformMap(
    MultidimArray<RFLOAT> &vol_in, MultidimArray<RFLOAT> &power_spectrum,
    int current_size, int nr_threads, bool do_gridding, bool do_heavy, int min_ires,
    const MultidimArray<RFLOAT> *fourier_mask, bool do_gpu
) {
    #ifdef CUDA
    static int semid = -1;
    CudaCustomAllocator *allocator;
    AccPtrFactory ptrFactory;
    AccPtr<RFLOAT>  dMpad;
    AccPtr<Complex> dFaux;
    AccPtr<RFLOAT>  dvol;
    AccPtr<Complex> ddata;
    AccPtr<RFLOAT> dfourier_mask;
    AccPtr<RFLOAT> dpower_spectrum;
    AccPtr<RFLOAT> dcounter;
    const cufftType cufft_type = sizeof(RFLOAT) == sizeof(double) ? CUFFT_D2Z : CUFFT_R2C;
    struct sembuf op_lock[2] = {
        0, 0, 0, /* wait for sem #0 to become 0 */
        0, 1, SEM_UNDO /* then increment sem #0 by 1 */
    };
    struct sembuf op_unlock[1] = { 0, -1, (IPC_NOWAIT | SEM_UNDO) /* decrement sem #0 by 1 (sets it to 0) */ };
    int fmXsz, fmYsz, fmZsz;
    size_t mem_req, ws_sz;
    #endif
    RFLOAT normfft;
    int padoridim;
    int n[3];
    bool do_fourier_mask = !!fourier_mask;
    {
    ifdefTIMING(TicToc tt (proj_timer, TIMING_TOP);)

    MultidimArray<Complex> Faux;
    MultidimArray<RFLOAT> Mpad;
    FourierTransformer transformer;
    {
    ifdefTIMING(TicToc tt (proj_timer, TIMING_INIT1);)

    // Size of padded real-space volume
    padoridim = round(padding_factor * ori_size);
    // make sure padoridim is even
    padoridim += padoridim % 2;
    // Re-calculate padding factor
    padding_factor = (float) padoridim / (float) ori_size;

    // Initialize data array of the oversampled transform
    ref_dim = vol_in.getDim();

    // Make Mpad
    switch (ref_dim) {

        case 2:
        if (do_heavy) {
            Mpad.initZeros(padoridim, padoridim);
        } else {
            Mpad.reshape(padoridim, padoridim);
        }
        normfft = (RFLOAT) (
            data_dim == 2 ?
            padding_factor * padding_factor :
            padding_factor * padding_factor * ori_size
        );
        break;

        case 3:
        if (do_heavy) {
            Mpad.initZeros(padoridim, padoridim, padoridim);
        } else {
            Mpad.reshape(padoridim, padoridim, padoridim);
        }
        normfft = (RFLOAT) (
            data_dim == 3 ?
            padding_factor * padding_factor * padding_factor :
            padding_factor * padding_factor * padding_factor * ori_size
        );
        break;

        default:
        REPORT_ERROR((std::string) "Projector::" + __func__ + "%%ERROR: Dimension of the data array should be 2 or 3");
    }
    #ifdef CUDA

    const int Faux_sz = (padoridim / 2 + 1) * padoridim * (ref_dim == 3 ? padoridim : 1);
    n[0] = padoridim; n[1] = padoridim; n[2] = padoridim;

    mem_req = (size_t) 1024;
    if (do_heavy && do_gpu) {
        cufftEstimateMany(ref_dim, n, nullptr, 0, 0, nullptr, 0, 0, cufft_type, 1, &ws_sz);

        mem_req = (size_t) sizeof(RFLOAT)  * vol_in.size()  // dvol
                + (size_t) sizeof(Complex) * Faux_sz        // dFaux
                + (size_t) sizeof(RFLOAT)  * Mpad.size()    // dMpad
                + ws_sz + 4096;                             // workspace for cuFFT + extra space for alignment
    }

    CudaCustomAllocator *allocator = nullptr;
    if (do_gpu && do_heavy) {
        cudaDeviceProp devProp;
        if (semid < 0) {
            int devid;
            HANDLE_ERROR(cudaGetDevice(&devid));
            if ((semid = semget(SEMKEY + devid, 1, IPC_CREAT | PERMS)) < 0)
                REPORT_ERROR("semget error");
        }
        if (semop(semid, &op_lock[0], 2) < 0)
            REPORT_ERROR("semop lock error");

        size_t mem_free, mem_tot;
        HANDLE_ERROR(cudaMemGetInfo(&mem_free, &mem_tot));
        if (mem_free > mem_req) {
            allocator = new CudaCustomAllocator(mem_req, (size_t) 16);
        } else {
            do_gpu = false; // change local copy of do_gpu variable
            if (semop(semid, &op_unlock[0], 1) < 0)
                REPORT_ERROR("semop unlock error");
        }
    }
    ptrFactory = allocator;
    dMpad = ptrFactory.make<RFLOAT>(Mpad.size());
    dFaux = ptrFactory.make<Complex>(Faux_sz);
    dvol  = ptrFactory.make<RFLOAT>(vol_in.size());
    if (do_heavy && do_gpu) {
        dvol.setHostPtr(vol_in.data);
        dvol.accAlloc();
        dvol.cpToDevice();
    }
    #endif
    }

    {
    ifdefTIMING(TicToc tt (proj_timer, TIMING_GRID);)
    // First do a gridding pre-correction on the real-space map:
    // Divide by the inverse Fourier transform of the interpolator in Fourier-space
    // 10feb11: at least in 2D case, this seems to be the wrong thing to do!!!
    // TODO: check what is best for subtomo!
    if (do_gridding) {
        // && data_dim != 3)
        if (do_heavy) {
            #ifdef CUDA
            if (do_gpu) {
                vol_in.setXmippOrigin();
                run_griddingCorrect(
                    ~dvol, interpolator, (RFLOAT) (ori_size * padding_factor), r_min_nn,
                    Xsize(vol_in), Ysize(vol_in), Zsize(vol_in)
                );
            } else
            #endif
            griddingCorrect(vol_in);
        } else {
            vol_in.setXmippOrigin();
        }
    }
    }

    {
    ifdefTIMING(TicToc tt (proj_timer, TIMING_PAD);)
    // Pad translated map with zeros
    vol_in.setXmippOrigin();
    Mpad.setXmippOrigin();
    if (do_heavy) {
        #ifdef CUDA
        if (do_gpu) {
            dMpad.accAlloc();
            run_padTranslatedMap(
                ~dvol, ~dMpad,
                Xinit(vol_in), Xlast(vol_in), Yinit(vol_in), Ylast(vol_in), Zinit(vol_in), Zlast(vol_in),   // Input  dimensions
                Xinit(Mpad),   Xlast(Mpad),   Yinit(Mpad),   Ylast(Mpad),   Zinit(Mpad),   Zlast(Mpad)      // Output dimensions
            );
            dvol.freeDevice();
        }
        else
        #endif
        FOR_ALL_ELEMENTS_IN_ARRAY3D(vol_in) // This will also work for 2D
            Mpad.elem(i, j, k) = vol_in.elem(i, j, k);
    }
    }

    const size_t cudanormfft = (size_t) padoridim * (size_t) padoridim * (ref_dim == 3 ? (size_t) padoridim : 1);
    {
    ifdefTIMING(TicToc tt (proj_timer, TIMING_TRANS);)
    // Calculate the oversampled Fourier transform
    if (do_heavy)
        #ifdef CUDA
        if (do_gpu) {
            dFaux.accAlloc();
            cufftResult err;
            cufftHandle plan;

            err = cufftCreate(&plan);
            if (err != CUFFT_SUCCESS)
                REPORT_ERROR("failed to create cufft plan");
            cufftSetAutoAllocation(plan, 0); // do not allocate work area
            //Allocate space with smart allocator
            AccPtr<char> fft_ws = ptrFactory.make<char>(ws_sz);
            fft_ws.accAlloc();
            cufftSetWorkArea(plan, ~fft_ws);
            err = cufftMakePlanMany(plan, ref_dim, n, nullptr, 0, 0, nullptr, 0, 0, cufft_type, 1, &ws_sz);
            if (err != CUFFT_SUCCESS)
                REPORT_ERROR("failed to create cufft plan");

            // do inverse FFT (dMpad->dFaux)
            if (sizeof(RFLOAT) == sizeof(double))
                err = cufftExecD2Z(plan, (cufftDoubleReal*) ~dMpad, (cufftDoubleComplex*) ~dFaux);
            else
                err = cufftExecR2C(plan, (cufftReal*) ~dMpad, (cufftComplex*) ~dFaux);
            if (err != CUFFT_SUCCESS)
                REPORT_ERROR("failed to exec fft");
            // deallocate plan, free mem
            cufftDestroy(plan);
            fft_ws.free();

            if (ref_dim == 2) Faux.reshape(padoridim / 2 + 1, padoridim);
            if (ref_dim == 3) Faux.reshape(padoridim / 2 + 1, padoridim, padoridim);

            scale((RFLOAT*) ~dFaux, 2 * dFaux.getSize(), 1.0 / (RFLOAT) cudanormfft);
        } else
        #endif
        Faux = transformer.FourierTransform(Mpad);
    }

    {
    ifdefTIMING(TicToc tt (proj_timer, TIMING_CENTER);)
    // Translate padded map to put origin of FT in the center
    if (do_heavy)
        #ifdef CUDA
        if (do_gpu)
        run_CenterFFTbySign(~dFaux, Xsize(Faux), Ysize(Faux), Zsize(Faux));
        else
        #endif
        CenterFFTbySign(Faux);
    }

    {
    ifdefTIMING(TicToc tt (proj_timer, TIMING_INIT2);)
    // Free memory: Mpad no longer needed
    #ifdef CUDA
    dMpad.free();
    #endif
    Mpad.clear();

    // Resize data array to the right size and initialise to zero
    initZeros(current_size);

    // Fill data only for those points with distance to origin less than max_r
    // (other points will be zero because of initZeros() call above
    // Also calculate radial power spectrum
    #ifdef CUDA
    int fourier_mask_sz = do_fourier_mask ? fourier_mask->size() : 16;
    ddata = ptrFactory.make<Complex>(data.size());
    dfourier_mask = ptrFactory.make<RFLOAT>(fourier_mask_sz);
    dpower_spectrum = ptrFactory.make<RFLOAT>(ori_size / 2 + 1);
    dcounter = ptrFactory.make<RFLOAT>(ori_size / 2 + 1);

    if (do_heavy && do_gpu) {
        ddata.accAlloc();
        dpower_spectrum.accAlloc();
        dcounter.accAlloc();
        ddata.deviceInit(0);
        dpower_spectrum.deviceInit(0);
        dcounter.deviceInit(0);
        dfourier_mask.accAlloc();
        fmXsz = fmYsz = fmZsz = 0;
        if (do_fourier_mask) {
            dfourier_mask.setHostPtr(fourier_mask->data);
            dfourier_mask.cpToDevice();
            fmXsz = Xsize(*fourier_mask);
            fmYsz = Ysize(*fourier_mask);
            fmZsz = Zsize(*fourier_mask);
        }
    }
    #endif
    }
    // Maybe power_spectrum.fill(0)?
    power_spectrum = MultidimArray<RFLOAT>::zeros(ori_size / 2 + 1);
    auto counter   = MultidimArray<RFLOAT>::zeros(ori_size / 2 + 1);

    {
    ifdefTIMING(TicToc tt (proj_timer, TIMING_FAUX);)
    int max_r2 = round(r_max * padding_factor) * round(r_max * padding_factor);
    int min_r2 = min_ires > 0 ? round(min_ires * padding_factor) * round(min_ires * padding_factor) : -1;

    if (do_heavy) {
        RFLOAT weight = 1.0;
        #ifdef CUDA
        if (do_gpu) {
            run_calcPowerSpectrum(
                ~dFaux, padoridim, ~ddata, Ysize(data), ~dpower_spectrum, ~dcounter,
                max_r2, min_r2, cudanormfft, padding_factor, weight,
                ~dfourier_mask, fmXsz, fmYsz, fmZsz, do_fourier_mask, ref_dim == 3
            );
            ddata.setHostPtr(data.data);
            ddata.cpToHost();
            dfourier_mask.free();
        } else
        #endif
        FOR_ALL_ELEMENTS_IN_FFTW_TRANSFORM(Faux) {
            // This will also work for 2D
            int r2 = euclidsq(ip, jp, kp);
            // The Fourier Transforms are all "normalised" for 2D transforms of size = ori_size x ori_size
            // Set data array
            if (r2 <= max_r2) {
                if (do_fourier_mask) {
                    weight = FFTW::elem(
                        *fourier_mask,
                        round(ip / padding_factor),
                        round(jp / padding_factor),
                        round(kp / padding_factor)
                    );
                }
                // Set data array
                data.elem(ip, jp, kp) = weight * direct::elem(Faux, i, j, k) * normfft;

                // Calculate power spectrum
                int ires = round(sqrt((RFLOAT) r2) / padding_factor);
                direct::elem(power_spectrum, ires) += norm(data.elem(ip, jp, kp)) / 2.0;
                // Divide by two because of the complex plane's two-dimensionality
                direct::elem(counter,        ires) += weight;

                // High-pass filter the reference
                if (r2 <= min_r2) { data.elem(ip, jp, kp) = 0; }
            }
        }
    }
    }

    /*
    FourierTransformer ft2;
    MultidimArray<Complex> Faux2(padding_factor * ori_size, (padding_factor * ori_size)/2+1);
    Image<RFLOAT> tt2(padding_factor * ori_size, padding_factor * ori_size);
    decenter(data, Faux2, max_r2);
    CenterFFTbySign(Faux2);
    windowFourierTransform(Faux2, padding_factor * ori_size);
    ft2.inverseFourierTransform(Faux2, tt2());
    tt2().setXmippOrigin();
    tt2().window(Xmipp::init(ori_size), Xmipp::init(ori_size), Xmipp::last(ori_size), Xmipp::last(ori_size));
    tt2.write("Fdata_proj.spi");
    std::cerr << "written Fdata_proj.spi" << std::endl;
    REPORT_ERROR("STOP");
    */

    {
    ifdefTIMING(TicToc tt (proj_timer, TIMING_POW);)
    // Calculate radial average of power spectrum
    if (do_heavy) {
        #ifdef CUDA
        if (do_gpu) {
            run_updatePowerSpectrum(~dcounter, dcounter.getSize(), ~dpower_spectrum);
            dpower_spectrum.setHostPtr(power_spectrum.data);
            dpower_spectrum.cpToHost();
        } else
        #endif
        for (long int i = 0; i < Xsize(power_spectrum); i++) {
            if (counter[i] < 1.0) {
                power_spectrum[i] = 0.0;
            } else {
                power_spectrum[i] /= counter[i];
            }
        }
    }
    }

    }
    #ifdef CUDA
    ddata.free();
    dpower_spectrum.free();
    dcounter.free();
    dvol.free();
    dMpad.free();
    dFaux.free();

    if (allocator) {
        delete allocator;
        if (semop(semid, &op_unlock[0], 1) < 0) REPORT_ERROR("semop unlock error");
    }
    #endif

    #ifdef PROJ_TIMING
    proj_timer.printTimes(false);
    #endif

}

void Projector::griddingCorrect(MultidimArray<RFLOAT> &vol_in) {
    // Correct real-space map by dividing it by the Fourier transform of the interpolator(s)
    vol_in.setXmippOrigin();
    FOR_ALL_ELEMENTS_IN_ARRAY3D(vol_in) {
        const RFLOAT r = euclid(i, j, k);
        if (r > 0.0) {
            const RFLOAT rval = r / (ori_size * padding_factor);
            const RFLOAT sinc_theta = sinc(PI * rval);
            // RFLOAT ftblob = blob_Fourier_val(rval, blob) / blob_Fourier_val(0.0, blob);
            // Interpolation (goes with "interpolator") to go from arbitrary to fine grid
            if (interpolator == NEAREST_NEIGHBOUR && r_min_nn == 0) {
                // NN interpolation is convolution with a rectangular pulse, which FT is a sinc function
                vol_in.elem(i, j, k) /= sinc_theta;
            } else if (interpolator == TRILINEAR || interpolator == NEAREST_NEIGHBOUR && r_min_nn > 0) {
                // trilinear interpolation is convolution with a triangular pulse, which FT is a sinc^2 function
                vol_in.elem(i, j, k) /= sinc_theta * sinc_theta;
            } else {
                REPORT_ERROR((std::string) "BUG Projector::" + __func__ + ": unrecognised interpolator scheme.");
            }
        // #define DEBUG_GRIDDING_CORRECT
        #ifdef DEBUG_GRIDDING_CORRECT
        if (i == 0 && j == 0 && k == 0)
            std::cerr << " j= " << j << " sinc= " << sinc_theta << std::endl;
        #endif
        }
    }
}

void Projector::project(MultidimArray<Complex> &f2d, const Matrix2D<RFLOAT> &A) const {
    // f2d should already be in the right size (ori_size,orihalfdim)
    // AND the points outside r_max should already be zero...
    // f2d.initZeros();

    // Use the inverse matrix

    const Matrix2D<RFLOAT> Ainv = A.inv() * (RFLOAT) padding_factor;  // Take scaling directly into account

    // The f2d image may be smaller than r_max, in that case also make sure not to fill the corners!
    const int r_max_out = Xsize(f2d) - 1;
    const int r_max_ref = r_max * padding_factor;

    const int r_min_NN_ref_2 = r_min_nn * r_min_nn * padding_factor * padding_factor;

    // #define DEBUG
    #ifdef DEBUG
    std::cerr << " Xsize(f2d)= "  << Xsize(f2d)  << '\n'
              << " Ysize(f2d)= "  << Ysize(f2d)  << '\n'
              << " Xsize(data)= " << Xsize(data) << '\n'
              << " Ysize(data)= " << Ysize(data) << '\n'
              << " Xinit(data)= " << Xinit(data) << '\n'
              << " Yinit(data)= " << Yinit(data) << '\n'
              << " Zinit(data)= " << Zinit(data) << '\n'
              << " max_r= "       << r_max       << '\n'
              << " Ainv= "        << Ainv        << std::endl;
    #endif

    const int r_max_out_2 = r_max_out * r_max_out;
    const int r_max_ref_2 = r_max_ref * r_max_ref;
    for (int i = 0; i < Ysize(f2d); i++) {
        const int y = i <= r_max_out ? i : i - Ysize(f2d);

        const int x_max = floor(sqrt(r_max_out_2 - y * y));

        for (int x = 0; x <= x_max; x++) {
            // Guaranteed that: Pythag(x, y) < r_max_out

            // Get logical coordinates in the 3D map
            RFLOAT xp = Ainv(0, 0) * x + Ainv(0, 1) * y;
            RFLOAT yp = Ainv(1, 0) * x + Ainv(1, 1) * y;
            RFLOAT zp = Ainv(2, 0) * x + Ainv(2, 1) * y;

            const RFLOAT r_ref_2 = euclidsq(xp, yp, zp);

            if (r_ref_2 > r_max_ref_2) continue;

            if (interpolator == TRILINEAR || r_ref_2 < r_min_NN_ref_2) {
                trilinear(xp, yp, zp, data, f2d, x, i);
            } else if (interpolator == NEAREST_NEIGHBOUR) {
                nearest_neighbour(xp, yp, zp, data, f2d, x, i);
            } else {
                REPORT_ERROR((std::string) "Unrecognized interpolator in Projector::" + __func__);
            }

        }
    }

    #ifdef DEBUG
    std::cerr << "done with project..." << std::endl;
    #endif
}


void Projector::projectGradient(Volume<t2Vector<Complex>>& img_out, Matrix2D<RFLOAT>& At) {
    const int sh = img_out.dimx;
    const int s  = img_out.dimy;

    const Matrix2D<RFLOAT> Ainv = At.inv() * (RFLOAT) padding_factor;  // Take scaling directly into account

    // Go from the 2D slice coordinates to the 3D coordinates

    for (int yy = 0; yy < s; yy++) {
        const double y = yy < sh ? yy : yy - s;

        for (int xx = 0; xx < sh; xx++) {
            const double x = xx;

            if (x * x + y * y > sh * sh) continue;

            // Get logical coordinates in the 3D map
            double xp = Ainv(0, 0) * x + Ainv(0, 1) * y;
            double yp = Ainv(1, 0) * x + Ainv(1, 1) * y;
            double zp = Ainv(2, 0) * x + Ainv(2, 1) * y;

            const bool is_neg_x = xp < 0;

            // Only asymmetric half is stored
            if (is_neg_x) {
                // Get complex conjugate hermitian symmetry pair
                xp = -xp;
                yp = -yp;
                zp = -zp;
            }

            // Trilinear interpolation (with physical coords)
            // Subtract Yinit and Zinit to accelerate access to data (Xinit = 0)
            // In that way use direct::elem, rather than A3D_ELEM

            int x0 = floor(xp);
            double fx = xp - x0;
            // x0 -= Xinit(data);
            int x1 = x0 + 1;

            int y0 = floor(yp);
            double fy = yp - y0;
            y0 -= Yinit(data);
            int y1 = y0 + 1;

            int z0 = floor(zp);
            double fz = zp - z0;
            z0 -= Zinit(data);
            int z1 = z0 + 1;

            if (
                x0 < 0 || x0 + 1 >= data.xdim ||
                y0 < 0 || y0 + 1 >= data.ydim ||
                z0 < 0 || z0 + 1 >= data.zdim
            ) {
                img_out(xx, yy, 0) = t2Vector<Complex>(Complex(0.0, 0.0), Complex(0.0, 0.0));
                continue;
            }

            Complex v000 = direct::elem(data, x0, y0, z0);
            Complex v001 = direct::elem(data, x1, y0, z0);
            Complex v010 = direct::elem(data, x0, y1, z0);
            Complex v011 = direct::elem(data, x1, y1, z0);
            Complex v100 = direct::elem(data, x0, y0, z1);
            Complex v101 = direct::elem(data, x1, y0, z1);
            Complex v110 = direct::elem(data, x0, y1, z1);
            Complex v111 = direct::elem(data, x1, y1, z1);

            Complex v00 = LIN_INTERP(fx, v000, v001);
            Complex v10 = LIN_INTERP(fx, v100, v101);
            Complex v01 = LIN_INTERP(fx, v010, v011);
            Complex v11 = LIN_INTERP(fx, v110, v111);

            Complex v0 = LIN_INTERP(fy, v00, v01);
            Complex v1 = LIN_INTERP(fy, v10, v11);

            // Complex v = LIN_INTERP(fz, v0, v1);

            Complex v00_dx = v001 - v000;
            Complex v10_dx = v101 - v100;
            Complex v01_dx = v011 - v010;
            Complex v11_dx = v111 - v110;
            Complex v0_dx = LIN_INTERP(fy, v00_dx, v01_dx);
            Complex v1_dx = LIN_INTERP(fy, v10_dx, v11_dx);
            Complex v_dx  = LIN_INTERP(fz, v0_dx,  v1_dx);

            Complex v0_dy = v01 - v00;
            Complex v1_dy = v11 - v10;
            Complex v_dy = LIN_INTERP(fz, v0_dy, v1_dy);

            Complex v_dz = v1 - v0;

            t3Vector<Complex> grad3D(v_dx, v_dy, v_dz);

            // Take complex conjugate for half with negative x
            if (is_neg_x) {
                grad3D.x = -(grad3D.x).conj();
                grad3D.y = -(grad3D.y).conj();
                grad3D.z = -(grad3D.z).conj();
            }

            img_out(xx, yy, 0).x = Ainv(0, 0) * grad3D.x + Ainv(1, 0) * grad3D.y + Ainv(2, 0) * grad3D.z;
            img_out(xx, yy, 0).y = Ainv(0, 1) * grad3D.x + Ainv(1, 1) * grad3D.y + Ainv(2, 1) * grad3D.z;
        }
    }
}

void Projector::project2Dto1D(MultidimArray<Complex> &f1d, const Matrix2D<RFLOAT> &A) const {
    // f1d should already be in the right size (ori_size,orihalfdim)
    // AND the points outside r_max should already be zero...
    // f1d.initZeros();

    const Matrix2D<RFLOAT> Ainv = A.inv() * (RFLOAT) padding_factor;  // Take scaling directly into account

    // The f2d image may be smaller than r_max, in that case also make sure not to fill the corners!
    const int r_max_out = Xsize(f1d) - 1;

    const int r_max_ref = r_max * padding_factor;
    const int r_max_ref_2 = r_max_ref * r_max_ref;

    const int r_min_NN_ref_2 = r_min_nn * r_min_nn * padding_factor * padding_factor;

    for (int x = 0; x <= r_max_out; x++) {
        // Get logical coordinates in the 2D map
        RFLOAT xp = Ainv(0, 0) * x;
        RFLOAT yp = Ainv(1, 0) * x;

        const RFLOAT r_ref_2 = xp * xp + yp * yp;

        if (r_ref_2 > r_max_ref_2) continue;

        if (interpolator == TRILINEAR || r_ref_2 < r_min_NN_ref_2) {
            trilinear(xp, yp, data, f1d, x);
        } else if (interpolator == NEAREST_NEIGHBOUR ) {
            nearest_neighbour(xp, yp, data, f1d, x);
        } else {
            REPORT_ERROR((std::string) "Unrecognized interpolator in Projector::" + __func__);
        }
    }
}


void Projector::rotate2D(MultidimArray<Complex> &f2d, const Matrix2D<RFLOAT> &A) const {
    // f2d should already be in the right size (ori_size,orihalfdim)
    // AND the points outside max_r should already be zero...
    // f2d.initZeros();

    // Use the inverse matrix
    const Matrix2D<RFLOAT> Ainv = A.inv() * (RFLOAT) padding_factor;  // Take scaling directly into account

    // The f2d image may be smaller than r_max, in that case also make sure not to fill the corners!
    const int r_max_out = Xsize(f2d) - 1;
    const int r_max_ref = r_max * padding_factor;

    const int r_min_NN_ref_2 = r_min_nn * r_min_nn * padding_factor * padding_factor;

    #ifdef DEBUG
    std::cerr << " Xsize(f2d)= "  << Xsize(f2d)  << '\n'
              << " Ysize(f2d)= "  << Ysize(f2d)  << '\n'
              << " Xsize(data)= " << Xsize(data) << '\n'
              << " Ysize(data)= " << Ysize(data) << '\n'
              << " Xinit(data)= " << Xinit(data) << '\n'
              << " Yinit(data)= " << Yinit(data) << '\n'
              << " Zinit(data)= " << Zinit(data) << '\n'
              << " max_r= " << r_max << '\n'
              << " Ainv= "  << Ainv  << std::endl;
    #endif

    const int r_max_out_2 = r_max_out * r_max_out;
    const int r_max_ref_2 = r_max_ref * r_max_ref;
    for (int i = 0; i < Ysize(f2d); i++) {
        const int y = i <= r_max_out ? i : i - Ysize(f2d);

        const int x_max = floor(sqrt(r_max_out_2 - y * y));

        for (int x = 0; x <= x_max; x++) {
            // Pythag(x, y) guaranteed to be < r_max_out

            RFLOAT xp = Ainv(0,0) * x + Ainv(0,1) * y;
            RFLOAT yp = Ainv(1,0) * x + Ainv(1,1) * y;

            const int r_ref_2 = xp * xp + yp * yp;

            if (r_ref_2 > r_max_ref_2) continue;

            if (interpolator == TRILINEAR || r_ref_2 < r_min_NN_ref_2) {
                const bool is_neg_x = xp < 0;

                // Only asymmetric half is stored
                if (is_neg_x) {
                    // Get complex conjugate hermitian symmetry pair
                    xp = -xp;
                    yp = -yp;
                }

                // Trilinear interpolation (with physical coords)
                // Subtract Yinit to accelerate access to data (Xinit=0)
                // In that way use direct::elem, rather than A3D_ELEM
                const int x0 = floor(xp);
                const RFLOAT fx = xp - x0;
                const int x1 = x0 + 1;

                int y0 = floor(yp);
                const RFLOAT fy = yp - y0;
                y0 -= Yinit(data);
                const int y1 = y0 + 1;

                // Matrix access can be accelerated through pre-calculation of z0*xydim etc.
                const Complex d00 = direct::elem(data, x0, y0);
                const Complex d01 = direct::elem(data, x1, y0);
                const Complex d10 = direct::elem(data, x0, y1);
                const Complex d11 = direct::elem(data, x1, y1);

                // Set the interpolated value in the 2D output array
                const Complex dx0 = LIN_INTERP(fx, d00, d01);
                const Complex dx1 = LIN_INTERP(fx, d10, d11);

                direct::elem(f2d, x, i) = LIN_INTERP(fy, dx0, dx1);

                // Take complex conjugate for half with negative x
                if (is_neg_x) direct::elem(f2d, x, i) = conj(direct::elem(f2d, x, i));

            } else if (interpolator == NEAREST_NEIGHBOUR) {
                /// NOTE: Unused
                const int x0 = round(xp);
                const int y0 = round(yp);

                direct::elem(f2d, x, i) = x0 < 0 ?
                    conj(data.elem(-x0, -y0)) : data.elem(x0, y0);

            } else {
                REPORT_ERROR("Unrecognized interpolator in Projector::project");
            }
        }
    }
}


void Projector::rotate3D(MultidimArray<Complex> &f3d, const Matrix2D<RFLOAT> &A) const {
    // f3d should already be in the right size (ori_size,orihalfdim)
    // AND the points outside max_r should already be zero
    // f3d.initZeros();

    // Use the inverse matrix
    const Matrix2D<RFLOAT> Ainv = A.inv() * (RFLOAT) padding_factor;  // Take scaling directly into account

    const int r_max_out = Xsize(f3d) - 1;
    const int r_max_ref = r_max * padding_factor;

    const int r_min_NN_ref_2 = r_min_nn * r_min_nn * padding_factor * padding_factor;

    #ifdef DEBUG
    std::cerr << " Xsize(f3d)= "  << Xsize(f3d)  << '\n'
              << " Ysize(f3d)= "  << Ysize(f3d)  << '\n'
              << " Xsize(data)= " << Xsize(data) << '\n'
              << " Ysize(data)= " << Ysize(data) << '\n'
              << " Xinit(data)= " << Xinit(data) << '\n'
              << " Yinit(data)= " << Yinit(data) << '\n'
              << " Zinit(data)= " << Zinit(data) << '\n'
              << " max_r= " << r_max << '\n'
              << " Ainv= "  << Ainv  << std::endl;
    #endif

    const int r_max_out_2 = r_max_out * r_max_out;
    const int r_max_ref_2 = r_max_ref * r_max_ref;
    for (int k = 0; k < Zsize(f3d); k++) {
    const int z = k <= r_max_out ? k : k - Zsize(f3d);
    for (int i = 0; i < Ysize(f3d); i++) {
    const int y = i <= r_max_out ? i : i - Ysize(f3d);

        const RFLOAT yyzz = y * y + z * z;

        // avoid negative square root
        if (yyzz > r_max_out_2) continue;

        const int x_max = floor(sqrt(r_max_out_2 - yyzz));

        for (int x = 0; x <= x_max; x++) {
            // Get logical coordinates in the 3D map
            RFLOAT xp = Ainv(0, 0) * x + Ainv(0, 1) * y + Ainv(0, 2) * z;
            RFLOAT yp = Ainv(1, 0) * x + Ainv(1, 1) * y + Ainv(1, 2) * z;
            RFLOAT zp = Ainv(2, 0) * x + Ainv(2, 1) * y + Ainv(2, 2) * z;

            const int r_ref_2 = euclidsq(xp, yp, zp);

            if (r_ref_2 > r_max_ref_2) continue;

            if (interpolator == TRILINEAR || r_ref_2 < r_min_NN_ref_2) {
                // Only asymmetric half is stored
                const bool is_neg_x = xp < 0;

                if (is_neg_x) {
                    // Get complex conjugate hermitian symmetry pair
                    xp = -xp;
                    yp = -yp;
                    zp = -zp;
                }

                // Trilinear interpolation (with physical coords)
                // Subtract Yinit to accelerate access to data (Xinit=0)
                // In that way use direct::elem, rather than A3D_ELEM
                const int x0 = floor(xp);
                const RFLOAT fx = xp - x0;
                // x0 -= Xinit(data);
                const int x1 = x0 + 1;

                int y0 = floor(yp);
                const RFLOAT fy = yp - y0;
                y0 -= Yinit(data);
                const int y1 = y0 + 1;

                int z0 = floor(zp);
                const RFLOAT fz = zp - z0;
                z0 -= Zinit(data);
                const int z1 = z0 + 1;

                // Matrix access can be accelerated through pre-calculation of z0*xydim etc.
                const Complex d000 = direct::elem(data, x0, y0, z0);
                const Complex d001 = direct::elem(data, x1, y0, z0);
                const Complex d010 = direct::elem(data, x0, y1, z0);
                const Complex d011 = direct::elem(data, x1, y1, z0);
                const Complex d100 = direct::elem(data, x0, y0, z1);
                const Complex d101 = direct::elem(data, x1, y0, z1);
                const Complex d110 = direct::elem(data, x0, y1, z1);
                const Complex d111 = direct::elem(data, x1, y1, z1);

                // Set the interpolated value in the 2D output array
                // interpolate in x
                const Complex dx00 = LIN_INTERP(fx, d000, d001);
                const Complex dx01 = LIN_INTERP(fx, d100, d101);
                const Complex dx10 = LIN_INTERP(fx, d010, d011);
                const Complex dx11 = LIN_INTERP(fx, d110, d111);
                // interpolate in y
                const Complex dxy0 = LIN_INTERP(fy, dx00, dx10);
                const Complex dxy1 = LIN_INTERP(fy, dx01, dx11);
                // interpolate in z
                direct::elem(f3d, x, i, k) = LIN_INTERP(fz, dxy0, dxy1);

                // Take complex conjugate for half with negative x
                if (is_neg_x) direct::elem(f3d, x, i, k) = conj(direct::elem(f3d, x, i, k));

            } else if (interpolator == NEAREST_NEIGHBOUR) {
                const int x0 = round(xp);
                const int y0 = round(yp);
                const int z0 = round(zp);

                direct::elem(f3d, x, i, k) = x0 < 0 ?
                    conj(data.elem(-x0, -y0, -z0)) : data.elem(x0, y0, z0);

            } else {
                REPORT_ERROR((std::string) "Unrecognized interpolator in Projector::" + __func__);
            }
        }
    }
    }
}
