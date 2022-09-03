/***************************************************************************
 *
 * Author: "Takanori Nakane"
 * MRC Laboratory of Molecular Biology
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * This complete copyright notice must be included in any revised version of the
 * source code. Additional authorship citations may be added, but existing
 * author citations must be preserved.
 ***************************************************************************/

#ifndef RWTIFF_H
#define RWTIFF_H

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "src/image.h"

inline DataType determine_datatype(uint16 bitsPerSample, uint16 sampleFormat) throw (int) {
    switch (bitsPerSample) {

        case 8:
        if (sampleFormat == SAMPLEFORMAT_UINT) return UChar;
        if (sampleFormat == SAMPLEFORMAT_INT) return SChar;
        break;

        case 16:
        if (sampleFormat == SAMPLEFORMAT_UINT) return UShort;
        if (sampleFormat == SAMPLEFORMAT_INT) return Short;
        break;

        case 32:
        if (sampleFormat == SAMPLEFORMAT_IEEEFP) return Float;

        default:
        throw 1;
    }
}

// I/O prototypes
/** TIFF Reader
  * @ingroup TIFF
*/
template <typename T>
int Image<T>::readTIFF(
    TIFF *ftiff, long int img_select,
    bool readdata, bool isStack, const FileName &name
) {
    // #define DEBUG_TIFF
    #ifdef DEBUG_TIFF
    printf("DEBUG readTIFF: Reading TIFF file. img_select %d\n", img_select);
    #endif

    // libtiff uses uint16 and uin32

    uint32 width, length;  // Apparent file dimensions
    if (
        TIFFGetField(ftiff, TIFFTAG_IMAGEWIDTH,  &width)  != 1 ||
        TIFFGetField(ftiff, TIFFTAG_IMAGELENGTH, &length) != 1
    ) REPORT_ERROR("The input TIFF file does not have the width or height field.");

    // True image dimensions
    std::array<unsigned long int, 4> dims { width, length, 1, 1 };

    uint16 bitsPerSample, sampleFormat;
    TIFFGetFieldDefaulted(ftiff, TIFFTAG_BITSPERSAMPLE, &bitsPerSample);
    TIFFGetFieldDefaulted(ftiff, TIFFTAG_SAMPLEFORMAT,  &sampleFormat);

    // Find the number of frames
    while (TIFFSetDirectory(ftiff, dims[3]) != 0) { dims[3]++; }
    // and go back to the start
    TIFFSetDirectory(ftiff, 0);

    #ifdef DEBUG_TIFF
    printf(
        "TIFF width %d, length %d, nDim %d, sample format %d, bits per sample %d\n",
        width, length, dims[3], sampleFormat, bitsPerSample
    );
    #endif

    // Detect 4-bit packed TIFFs. This is IMOD's own extension.
    // It is not easy to detect this format. Here we check only the image size.
    // See IMOD's iiTIFFCheck() in libiimod/iitif.c and sizeCanBe4BitK2SuperRes() in libiimod/mrcfiles.c.
    const bool packed_4bit = (bitsPerSample == 8 && (
        width == 5760 && length == 8184  || width == 8184  && length == 5760 || // K3 SR: 11520 x 8184
        width == 4092 && length == 11520 || width == 11520 && length == 4092 ||
        width == 3710 && length == 7676  || width == 7676  && length == 3710 || // K2 SR: 7676 x 7420
        width == 3838 && length == 7420  || width == 7420  && length == 3838
    ));

    if (packed_4bit) dims[0] *= 2;

    DataType datatype;
    try {
        datatype = packed_4bit ? UHalf : determine_datatype(bitsPerSample, sampleFormat);
    } catch (int) {
        std::cerr << "Unsupported TIFF format in " << name << ": sample format = " << sampleFormat << ", bits per sample = " << bitsPerSample << std::endl;
        REPORT_ERROR("Unsupported TIFF format.\n");
    }

    const long int i = MDMainHeader.index();
    MDMainHeader.setValue(EMDL::IMAGE_DATATYPE, (int) datatype, i);

    uint16 resolutionUnit;
    float xResolution;
    if (
        TIFFGetField(ftiff, TIFFTAG_RESOLUTIONUNIT, &resolutionUnit) == 1 &&
        TIFFGetField(ftiff, TIFFTAG_XRESOLUTION,    &xResolution)    == 1
    ) {
        // We don't support anistropic pixel size
        if (resolutionUnit == RESUNIT_INCH) {
            MDMainHeader.setValue(EMDL::IMAGE_SAMPLINGRATE_X, RFLOAT(2.54E8 / xResolution), i);  // 1 inch = 2.54 cm
            MDMainHeader.setValue(EMDL::IMAGE_SAMPLINGRATE_Y, RFLOAT(2.54E8 / xResolution), i);
        } else if (resolutionUnit == RESUNIT_CENTIMETER) {
            MDMainHeader.setValue(EMDL::IMAGE_SAMPLINGRATE_X, RFLOAT(1.00E8 / xResolution), i);
            MDMainHeader.setValue(EMDL::IMAGE_SAMPLINGRATE_Y, RFLOAT(1.00E8 / xResolution), i);
        }
        #ifdef DEBUG_TIFF
        std::cout << "resolutionUnit = " << resolutionUnit << " xResolution = " << xResolution << std::endl;
        std::cout << "pixel size = " << (RFLOAT) MDMainHeader.getValue(EMDL::IMAGE_SAMPLINGRATE_X) << std::endl;
        #endif
    }

    // TODO: TIFF is always a stack, isn't it?
    if (isStack) {
        dims[2] = 1;
        replaceNsize = dims[3];
        if (img_select >= (int) dims[3]) {
            std::string Num1str = std::to_string(img_select + 1);
            std::string Num2str = std::to_string(dims[3]);
            REPORT_ERROR((std::string) "readTIFF: Image number " + Num1str + " exceeds stack size " + Num2str + " of image " + name);
        }
    } else {
        replaceNsize = 0;
    }

    // Map the parameters
    if (isStack) {
        dims[2] = 1;
        if (img_select != -1) {
            dims[3] = 1;
        }
    }
    data.setDimensions(dims[0], dims[1], dims[2], dims[3]);
    data.coreAllocateReuse();

    /*
    if (header->mx && header->a != 0)  // ux
        MDMainHeader.setValue(EMDL::IMAGE_SAMPLINGRATE_X, (RFLOAT) header->a / header->mx);
    if (header->my && header->b != 0)  // yx
        MDMainHeader.setValue(EMDL::IMAGE_SAMPLINGRATE_Y, (RFLOAT) header->b / header->my);
    if (header->mz && header->c != 0)  // zx
        MDMainHeader.setValue(EMDL::IMAGE_SAMPLINGRATE_Z, (RFLOAT) header->c / header->mz);
    */

    if (readdata) {
        if (img_select == -1) { img_select = 0; }  // img_select starts from 0

        size_t haveread_n = 0;
        for (int i = 0; i < dims[3]; i++) {
            TIFFSetDirectory(ftiff, img_select);

            // Make sure image property is consistent for all frames

            uint32 cur_width, cur_length;
            if (
                TIFFGetField(ftiff, TIFFTAG_IMAGEWIDTH,  &cur_width)  != 1 ||
                TIFFGetField(ftiff, TIFFTAG_IMAGELENGTH, &cur_length) != 1
            ) REPORT_ERROR(name + ": The input TIFF file does not have the width or height field.");

            uint16 cur_bitsPerSample, cur_sampleFormat;
            TIFFGetFieldDefaulted(ftiff, TIFFTAG_BITSPERSAMPLE, &cur_bitsPerSample);
            TIFFGetFieldDefaulted(ftiff, TIFFTAG_SAMPLEFORMAT,  &cur_sampleFormat);

            if (
                cur_width         != width         ||
                cur_length        != length        ||
                cur_bitsPerSample != bitsPerSample ||
                cur_sampleFormat  != sampleFormat
            ) REPORT_ERROR(name + ": All frames in a TIFF should have same width, height and pixel format.\n");

            tsize_t stripSize = TIFFStripSize(ftiff);
            tstrip_t numberOfStrips = TIFFNumberOfStrips(ftiff);
            tdata_t buf = _TIFFmalloc(stripSize);
            #ifdef DEBUG_TIFF
            size_t readsize_n = stripSize * 8 / bitsPerSample;
            std::cout << "TIFF stripSize=" << stripSize << " numberOfStrips=" << numberOfStrips << " readsize_n=" << readsize_n << std::endl;
            #endif
            for (tstrip_t strip = 0; strip < numberOfStrips; strip++) {
                tsize_t actually_read = TIFFReadEncodedStrip(ftiff, strip, buf, stripSize);
                if (actually_read == -1)
                    REPORT_ERROR((std::string) "Failed to read an image data from " + name);
                tsize_t actually_read_n = actually_read * 8 / bitsPerSample;
                #ifdef DEBUG_TIFF
                std::cout << "Reading strip: " << strip << "actually read byte:" << actually_read << std::endl;
                #endif
                if (packed_4bit) {
                    actually_read_n *= 2; // convert physical size to logical size
                }
                castPage2T((char*) buf, data.data + haveread_n, datatype, actually_read_n);
                haveread_n += actually_read_n;
            }

            _TIFFfree(buf);
            img_select++;
        }

        /* Flip the Y axis.

           In an MRC file, the origin is bottom-left, +X to the right, +Y to the top.
           (c.f. Fig. 2 of Heymann et al, JSB 2005 https://doi.org/10.1016/j.jsb.2005.06.001
           IMOD's interpretation http://bio3d.colorado.edu/imod/doc/mrc_format.txt)
           3dmod (from IMOD) and e2display.py (from EMAN2) display like this.

           relion_display has the origin at top-left, +X to the right, +Y to the bottom.
           GIMP and ImageJ display in this way as well.
           A TIFF file, with TIFFTAG_ORIENTATION = 1 (default), shares this convention.

           So, the origin and the direction of the Y axis are the opposite between MRC and TIFF.
           IMOD, EMAN2, SerialEM and MotionCor2 flip the Y axis whenever they read or write a TIFF file.
           We follow this.
        */

        const int ylim = dims[1] / 2, z = 0;
        for (int n = 0; n < dims[3]; n++)
        for (int y1 = 0; y1 < ylim; y1++) {
            const int y2 = dims[1] - 1 - y1;
            for (int x = 0; x < dims[0]; x++) {
                /// TODO: memcpy or pointer arithmetic is probably faster
                std::swap(
                    direct::elem(data, x, y1, z, n),
                    direct::elem(data, x, y2, z, n)
                );
            }
        }
    }

    return 0;
}

#endif
