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
#include <memory>
#include "src/image.h"
#include "src/rwSPIDER.h"
#include "src/rwIMAGIC.h"
#include "src/rwMRC.h"
#include "src/rwTIFF.h"

// #define DEBUG_REGULARISE_HELICAL_SEGMENTS

// Z-score of x given mean mu and standard deviation sigma
inline RFLOAT Z(RFLOAT x, RFLOAT mu, RFLOAT sigma) {
    return (x - mu) / sigma;
}

// Get size of datatype
size_t gettypesize(DataType type) throw (RelionError) {

    switch (type) {

        case UChar:
        case SChar:
        return sizeof(char);

        case UShort:
        case Short:
        return sizeof(short);

        case UInt:
        case Int:
        return sizeof(int);

        case Float:
        return sizeof(float);

        case Double:
        return sizeof(RFLOAT);

        case Boolean:
        return sizeof(bool);

        case UHalf:
        unsigned long size;
        REPORT_ERROR("Logic error: UHalf (4-bit) needs special consideration. Don't use this function.");
        return size;

        default:
        return 0;
    }
}

// int DataType::String2Int(std::string s) {
//     toLower(s);
//     if (!strcmp(s.c_str(), "uchar"))  return UChar;
//     if (!strcmp(s.c_str(), "ushort")) return UShort;
//     if (!strcmp(s.c_str(), "short"))  return Short;
//     if (!strcmp(s.c_str(), "uint"))   return UInt;
//     if (!strcmp(s.c_str(), "int"))    return Int;
//     if (!strcmp(s.c_str(), "float"))  return Float;
//     REPORT_ERROR("DataType::String2int; unknown datatype");
// }

template <typename T>
int Image<T>::allocatePage(
    size_t pagesize, size_t off,
    DataType datatype, size_t datatypesize
) {
    // pagesize: size of object
    static const size_t pagemax = 0x40000000;  // 1 GB (1 << 30)
    const size_t size = std::max(pagesize, pagemax) * sizeof(char);
    const auto deleter = [size] (char *ptr) { callocator<char>::deallocate(ptr, size); };
    const auto page = std::unique_ptr<char, decltype(deleter)>(callocator<char>::allocate(size), deleter);
    // Because we requested XYsize to be even for UHalf, this is always safe.
    if (fseek(fimg, off, SEEK_SET) != 0) return -1;

    size_t pixel_progress = 0;  // Number of pixels processed so far
    for (size_t n = 0; n < Nsize(data); n++) {
        for (size_t j = 0; j < pagesize; j += pagemax) {
            // Read next page
            // Divide pages larger than pagemax
            size_t readsize = std::min(pagesize - j, pagemax);

            size_t readsize_n = datatype == UHalf ? readsize * 2 : readsize / datatypesize;

            #ifdef DEBUG
            std::cout << "NX = " << Xsize(data) << " NY = " << Ysize(data) << " NZ = " << Zsize(data) << std::endl;
            std::cout << "pagemax = " << pagemax << " pagesize = " << pagesize  << " readsize = " << readsize << " readsize_n = " << readsize_n << std::endl;
            #endif

            // Read page from disk
            if (fread(page.get(), readsize, 1, fimg) != 1) return -2;

            // Swap bytes if required
            if (swap) swapPage(page.get(), readsize, gettypesize(datatype), swap);
            // Cast to T
            castFromPage(data.data + pixel_progress, page.get(), datatype, readsize_n);
            pixel_progress += readsize_n;
        }
        if (pad > 0) {
            // fread(padpage, pad, 1, fimg);
            if (fseek(fimg, pad, SEEK_CUR) != 0) return -1;
        }
    }
    return 0;
}

// Some image-specific operations
void normalise(
    Image<RFLOAT> &I,
    int bg_radius,
    RFLOAT white_dust_stddev, RFLOAT black_dust_stddev,
    bool do_ramp, bool is_helical_segment,
    RFLOAT helical_mask_tube_outer_radius_pix,
    RFLOAT tilt_deg, RFLOAT psi_deg
) {

    if (2 * bg_radius > Xsize(I()))
        REPORT_ERROR("normalise ERROR: 2*bg_radius is larger than image size!");

    if (is_helical_segment && 2 * (helical_mask_tube_outer_radius_pix + 1) > Xsize(I()))
        REPORT_ERROR("normalise ERROR: Diameter of helical tube is larger than image size!");

    if (is_helical_segment && I().getDim() == 2) { tilt_deg = 0.0; }

    if (white_dust_stddev > 0.0 || black_dust_stddev > 0.0) {
        // Calculate initial avg and stddev values
        const auto stats = calculateBackgroundAvgStddev(
            I, bg_radius,
            is_helical_segment, helical_mask_tube_outer_radius_pix,
            tilt_deg, psi_deg
        );

        // Remove white and black noise
        if (white_dust_stddev > 0.0)
        removeDust(I, true, white_dust_stddev, stats.avg, stats.stddev);
        if (black_dust_stddev > 0.0)
        removeDust(I, false, black_dust_stddev, stats.avg, stats.stddev);

    }

    if (do_ramp) {
        subtractBackgroundRamp(
            I, bg_radius,
            is_helical_segment, helical_mask_tube_outer_radius_pix,
            tilt_deg, psi_deg
        );
    }

    // Calculate avg and stddev (also redo if dust was removed!)
    const auto stats = calculateBackgroundAvgStddev(
        I, bg_radius,
        is_helical_segment, helical_mask_tube_outer_radius_pix,
        tilt_deg, psi_deg
    );

    if (stats.stddev < 1e-10) {
        std::cerr << " WARNING! Stddev of image " << I.name() << " is zero! Skipping normalisation..." << std::endl;
    } else {
        // Normalise
        for (auto &x : I()) { x = Z(x, stats.avg, stats.stddev); }
    }
}

Stats<RFLOAT> calculateBackgroundAvgStddev(
    Image<RFLOAT> &I,
    int bg_radius,
    bool is_helical_segment,
    RFLOAT helical_mask_tube_outer_radius_pix,
    RFLOAT tilt_deg, RFLOAT psi_deg
) {
    int bg_radius2 = bg_radius * bg_radius;
    RFLOAT sum = 0.0, sum_of_squares = 0.0, avg = 0.0, stddev = 0.0;
    long int n = 0;

    if (is_helical_segment) {
        int dim = I().getDim();
        if (dim != 2 && dim != 3)
            REPORT_ERROR("image.cpp::calculateBackgroundAvgStddev(): 2D or 3D image is required!");
        if (dim == 2) { tilt_deg = 0.0; }

        Matrix1D<RFLOAT> coords;
        // Init coords
        coords.clear();
        coords.resize(3);
        coords.initZeros();

        // Rotate the particle (helical axes are X and Z for 2D and 3D segments respectively)
        Matrix2D<RFLOAT> A = Euler::angles2matrix(0.0, tilt_deg, psi_deg).transpose();

        // Refer to the code in calculateBackgroundAvgStddev() for 3D implementation

        #ifdef DEBUG_REGULARISE_HELICAL_SEGMENTS
        int angle = round(fabs(psi_deg));
        FileName fn_test = integerToString(angle);
        if (psi_deg < 0.0) { fn_test = fn_test.addExtension("neg"); }
        fn_test = fn_test.addExtension("mrc");
        Image<RFLOAT> img_test;
        img_test.clear();
        img_test().resize(I());
        img_test().initZeros();
        std::cout << "FileName = " << fn_test.c_str() << std::endl;
        #endif

        // Calculate avg in the background pixels
        FOR_ALL_ELEMENTS_IN_ARRAY3D(I(), i, j, k) {
            // X, Y, Z coordinates
            ZZ(coords) = dim == 3 ? (RFLOAT) k : 0.0;
            YY(coords) =            (RFLOAT) i;
            XX(coords) =            (RFLOAT) j;
            // Rotate
            coords = A * coords;

            // Distance from the point to helical axis (perpendicular to X axis)
            RFLOAT d = dim == 3 ? sqrt(YY(coords) * YY(coords) + XX(coords) * XX(coords)) : abs(YY(coords));

            if (d > helical_mask_tube_outer_radius_pix) {
                RFLOAT x = I().elem(i, j, k);
                sum += x;
                sum_of_squares += x * x;
                n += 1;

                #ifdef DEBUG_REGULARISE_HELICAL_SEGMENTS
                img_test().elem(i, j, k) = 1.0;  // Mark bg pixels as 1, others as 0
                #endif
            }
        }
        if (n < 1) {
            REPORT_ERROR("image.cpp::calculateBackgroundAvgStddev(): No pixels in background are found. Radius of helical mask is too large.");
        }

        avg    = sum / (RFLOAT) n;
        stddev = sqrt(sum_of_squares / (RFLOAT) n - avg * avg);

        #ifdef DEBUG_REGULARISE_HELICAL_SEGMENTS
        img_test.write(fn_test);
        #endif
    } else {
        // Calculate avg in the background pixels
        FOR_ALL_ELEMENTS_IN_ARRAY3D(I(), i, j, k) {
            if (i * i + j * j + k * k > bg_radius2) {
                RFLOAT x = I().elem(i, j, k);
                sum += x;
                sum_of_squares += x * x;
                n += 1;
            }
        }
        if (n < 1) {
            REPORT_ERROR("image.cpp::calculateBackgroundAvgStddev(): No pixels in background are found. Radius of circular mask is too large.");
        }

        avg = sum / (RFLOAT) n;
        stddev = sqrt(sum_of_squares / (RFLOAT) n - avg * avg);
    }

    return { avg, stddev };
}

void subtractBackgroundRamp(
    Image<RFLOAT> &I,
    int bg_radius,
    bool is_helical_segment,
    RFLOAT helical_mask_tube_outer_radius_pix,
    RFLOAT tilt_deg, RFLOAT psi_deg
) {

    int bg_radius2 = bg_radius * bg_radius;
    fit_point3D point;
    std::vector<fit_point3D> allpoints;
    RFLOAT avgbg, stddevbg, minbg, maxbg;

    if (I().getDim() == 3)
        REPORT_ERROR("ERROR %% calculateBackgroundRamp is not implemented for 3D data!");

    if (is_helical_segment) {
        // not implemented for 3D data
        if (I().getDim() == 2) { tilt_deg = 0.0; }

        auto coords = Matrix1D<RFLOAT>::zeros(3);

        // Rotate the particle (helical axes are X and Z for 2D and 3D segments respectively)
        // Since Z = 0, tilt_deg does not matter
        Matrix2D<RFLOAT> A = Euler::angles2matrix(0.0, tilt_deg, psi_deg).transpose();

        FOR_ALL_ELEMENTS_IN_ARRAY2D(I(), i, j) {
            // not implemented for 3D data
            XX(coords) = (RFLOAT) i;
            YY(coords) = (RFLOAT) j;
            ZZ(coords) = 0.0;
            // Rotate
            coords = A * coords;
            if (abs(YY(coords)) > helical_mask_tube_outer_radius_pix) {
                // Not implemented for 3D data
                point.x = i;
                point.y = j;
                point.z = I().elem(i, j);
                point.w = 1.0;
                allpoints.push_back(point);
            }
        }
        if (allpoints.size() < 5) {
            REPORT_ERROR("image.cpp::subtractBackgroundRamp(): Less than 5 pixels in background are found. Radius of helical mask is too large.");
        }
    } else {
        FOR_ALL_ELEMENTS_IN_ARRAY2D(I(), i, j) {
            if (i * i + j * j > bg_radius2) {
                point.x = i;
                point.y = j;
                point.z = I().elem(i, j);
                point.w = 1.0;
                allpoints.push_back(point);
            }
        }
    }

    RFLOAT pA, pB, pC;
    fitLeastSquaresPlane(allpoints, pA, pB, pC);

    // Substract the plane from the image
    FOR_ALL_ELEMENTS_IN_ARRAY2D(I(), i, j) {
        I().elem(i, j) -= pA * i + pB * j + pC;
    }

}

void removeDust(
    Image<RFLOAT> &I, bool is_white, RFLOAT thresh, RFLOAT avg, RFLOAT stddev
) {
    FOR_ALL_ELEMENTS_IN_ARRAY3D(I(), i, j, k) {
        RFLOAT aux = I().elem(i, j, k);
        if (
             is_white && aux - avg >  thresh * stddev ||
            !is_white && aux - avg < -thresh * stddev
        ) {
            I().elem(i, j, k) = rnd_gaus(avg, stddev);
        }
    }
}


void invert_contrast(Image<RFLOAT> &I) {
    for (auto &x : I()) { x = -x; }
}

void rescale(Image<RFLOAT> &I, int mysize) {
    int olddim = Xsize(I());

    resizeMap(I(), mysize);

    // Try to rescale entries in I.MDmainheader
    const long int i = I.header.size() - 1;
    try {
        const auto sx = I.header.getValue<RFLOAT>(EMDL::IMAGE_SAMPLINGRATE_X, i) * (RFLOAT) olddim / (RFLOAT) mysize;
        I.header.setValue(EMDL::IMAGE_SAMPLINGRATE_X, sx, i);
    } catch (const char *errmsg) {}
    try {
        const auto sy = I.header.getValue<RFLOAT>(EMDL::IMAGE_SAMPLINGRATE_Y, i) * (RFLOAT) olddim / (RFLOAT) mysize;
        I.header.setValue(EMDL::IMAGE_SAMPLINGRATE_Y, sy, i);
    } catch (const char *errmsg) {}
    if (I().getDim() == 3)
    try {
        const auto sz = I.header.getValue<RFLOAT>(EMDL::IMAGE_SAMPLINGRATE_Z, i) * (RFLOAT) olddim / (RFLOAT) mysize;
        I.header.setValue(EMDL::IMAGE_SAMPLINGRATE_Z, sz, i);
    } catch (const char *errmsg) {}
}

void rewindow(Image<RFLOAT> &I, int size) {
    // Check 2D or 3D dimensionality
    switch (I().getDim()) {

        case 2: I() = I().windowed(
            Xmipp::init(size), Xmipp::last(size),
            Xmipp::init(size), Xmipp::last(size)
        ); break;

        case 3: I() = I().windowed(
            Xmipp::init(size), Xmipp::last(size),
            Xmipp::init(size), Xmipp::last(size),
            Xmipp::init(size), Xmipp::last(size)
        ); break;

    }
}

std::pair<RFLOAT, RFLOAT> getImageContrast(
    MultidimArray<RFLOAT> &image, RFLOAT minval, RFLOAT maxval, RFLOAT &sigma_contrast
) {
    // First check whether to apply sigma-contrast, 
    // i.e. set minval and maxval to the mean +/- sigma_contrast times the stddev
    bool redo_minmax = sigma_contrast > 0.0 || minval != maxval;

    if (sigma_contrast > 0.0 || minval == maxval) {
        const auto stats = computeStats(image);
        if (sigma_contrast > 0.0) {
            minval = stats.avg - stats.stddev * sigma_contrast;
            maxval = stats.avg + stats.stddev * sigma_contrast;
        } else {
            minval = stats.min;
            maxval = stats.max;
        }
    }

    if (redo_minmax) {
        // Constrain image data to the interval [minval, maxval]
        for (auto &x : image) {
            if (x > maxval) x = maxval;
            if (x < minval) x = minval;
        }
    }

    return { minval, maxval };
}

template <typename T>
int Image<T>::read(
    const FileName &name, bool readdata, long int select_img,
    bool mapData, bool is_2D
) {
    if (name.empty())
        REPORT_ERROR("ERROR: trying to read image with empty file name!");
    fImageHandler hFile;
    hFile.openFile(name);
    return _read(name, hFile, readdata, select_img, mapData, is_2D);
    // fImageHandler's destructor will close the file
}

template <typename T>
void Image<T>::write(
    FileName name, long int select_img, bool isStack,
    int mode
) {
    const FileName &fname = name.empty() ? filename : name;
    fImageHandler hFile;
    hFile.openFile(name, mode);
    _write(fname, hFile, select_img, isStack, mode);
    // fImageHandler's destructor will close the file
}

template <typename T>
int Image<T>::_read(
    const FileName &name, fImageHandler &hFile,
    bool readdata, long int select_img, bool mapData, bool is_2D
) {
    int err = 0;  // Exit code (negative is bad)

    // readdata indicates whether to read the data or only the header
    // We won't read the individual header and the data if not necessary

    // Check whether to map the data or not
    mmapper = mapData ? new image_mmapper : nullptr;

    FileName ext_name = hFile.ext_name;
    fimg = hFile.fimg;
    fhed = hFile.fhed;

    long int dump;
    name.decompose(dump, filename);
    // Subtract 1 to have numbering 0...N-1 instead of 1...N
    if (dump > 0) { dump--; }
    filename = name;

    if (select_img == -1) { select_img = dump; }

    #undef DEBUG
    // #define DEBUG
    #ifdef DEBUG
        std::cerr << "READ\n" <<
        "name=" << name << std::endl;
        std::cerr << "ext= " << ext_name << std::endl;
        std::cerr << " now reading: " << filename << " readdata= " << readdata
        << " select_img "  << select_img << std::endl;
    #endif
    #undef DEBUG

    // Clear the header before reading
    header.clear();
    header.addObject();

    if (
        ext_name.contains("spi") || ext_name.contains("xmp") ||
        ext_name.contains("stk") || ext_name.contains("vol")
    ) {
        const DataType datatype = readSPIDER(select_img);
        if (!isStack || readdata)
        err = readData(select_img, datatype);
    } else if (ext_name.contains("mrcs") ||
               ext_name.contains("mrc") && is_2D) {
        // MRC stack MUST go BEFORE plain MRC
        const DataType datatype = readMRC(select_img, true, name);
        if (!isStack || readdata)
        err = readData(select_img, datatype);
    } else if (ext_name.contains("tif")) {
        err = readTIFF(hFile.ftiff, select_img, readdata, true, name);
    } else if (select_img >= 0 && ext_name.contains("mrc")) {
        REPORT_ERROR("Image::read ERROR: stacks of images in MRC-format should have extension .mrcs; .mrc extensions are reserved for 3D maps.");
    } else if (ext_name.contains("mrc")) {
        // MRC 3D map
        const DataType datatype = readMRC(select_img, false, name);
        if (!isStack || readdata)
        err = readData(select_img, datatype);
    } else if (ext_name.contains("img") || ext_name.contains("hed")) {
        const DataType datatype = readIMAGIC(select_img);
        if (!isStack || readdata)
        err = readData(select_img, datatype);
    } else if (ext_name.contains("dm")) {
        REPORT_ERROR("The Digital Micrograph format (DM3, DM4) is not supported. You can convert it to MRC by other programs, for example, dm2mrc in IMOD.");
    } else if (ext_name.contains("eer") || ext_name.contains("ecc")) {
        REPORT_ERROR("BUG: EER movies should be handled by EERRenderer, not by Image.");
    } else {
        const DataType datatype = readSPIDER(select_img);
        if (!isStack || readdata)
        err = readData(select_img, datatype);
    }
    return err;
}

template <typename T>
void Image<T>::_write(
    const FileName &name, fImageHandler &hFile,
    long int select_img, bool isStack, int mode
) {
    int err = 0;

    FileName ext_name = hFile.ext_name;
    fimg = hFile.fimg;
    fhed = hFile.fhed;
    _exists = hFile.exist;

    filename = name;

    long int aux;
    FileName filNamePlusExt(name);
    name.decompose(aux, filNamePlusExt);
    // Subtract 1 to have numbering 0...N-1 instead of 1...N
    if (aux > 0)
        aux--;

    if (select_img == -1)
        select_img = aux;

    size_t found = filNamePlusExt.find_first_of("%");

    std::string imParam = "";

    if (found != std::string::npos) {
        imParam = filNamePlusExt.substr(found + 1).c_str();
        filNamePlusExt = filNamePlusExt.substr(0, found) ;
    }

    found = filNamePlusExt.find_first_of(":");
    if (found != std::string::npos)
        filNamePlusExt = filNamePlusExt.substr(0, found);

    // #define DEBUG
    #ifdef DEBUG

    std::cerr << "write" << std::endl;
    std::cerr << "extension for write= " << ext_name << std::endl;
    std::cerr << "filename= " << filename << std::endl;
    std::cerr << "mode= " << mode << std::endl;
    std::cerr << "isStack= " << isStack << std::endl;
    std::cerr << "select_img= " << select_img << std::endl;
    #endif
    #undef DEBUG
    // Check that image is not empty
    if (getSize() < 1)
        REPORT_ERROR("write Image ERROR: image is empty!");

    replaceNsize = 0; // reset replaceNsize in case image is reused
    if (select_img == -1 && mode == WRITE_REPLACE) {
        REPORT_ERROR("write: Please specify object to be replaced");
    } else if (!_exists && mode == WRITE_REPLACE) {
        std::stringstream replace_number;
        replace_number << select_img;
        REPORT_ERROR((std::string)
            "Cannot replace object number: " + replace_number.str()
            + " in file " + filename + ". It does not exist"
        );
    } else if (_exists && (mode == WRITE_REPLACE || mode == WRITE_APPEND)) {
        // CHECK FOR INCONSISTENCIES BETWEEN data.xdim and x, etc???
        Dimensions dimensions = this->getDimensions();

        Image<T> auxI;
        auxI.read(filNamePlusExt, false);
        Dimensions auxIdimensions = auxI.getDimensions();
        replaceNsize = auxIdimensions.n;
        if (
            dimensions.x != auxIdimensions.x ||
            dimensions.y != auxIdimensions.y ||
            dimensions.z != auxIdimensions.z
        )
            REPORT_ERROR("write: target and source objects have different size");
        if (mode == WRITE_REPLACE && select_img > auxIdimensions.n)
            REPORT_ERROR("write: cannot replace image stack is not large enough");
        if (auxI.replaceNsize < 1)
            REPORT_ERROR("write: output file is not an stack");
    } else if (!_exists && mode == WRITE_APPEND) {} else if (mode == WRITE_READONLY) {
        // If new file we are in the WRITE_OVERWRITE mode
        REPORT_ERROR((std::string) "File " + name + " opened in read-only mode. Cannot write.");
    }

    // Select format
    if (
        ext_name.contains("spi") || ext_name.contains("xmp") ||
        ext_name.contains("stk") || ext_name.contains("vol")
    ) {
        err = writeSPIDER(select_img, isStack, mode);
    } else if (ext_name.contains("mrcs")) {
        writeMRC(select_img, true, mode);
    } else if (ext_name.contains("mrc")) {
        writeMRC(select_img, false, mode);
    } else if (ext_name.contains("img") || ext_name.contains("hed")) {
        writeIMAGIC(select_img, mode);
    } else {
        err = writeSPIDER(select_img, isStack, mode);
    }
    if (err < 0) {
        std::cerr << " Filename = " << filename << " Extension= " << ext_name << std::endl;
        REPORT_ERROR((std::string) "Error writing file " + filename + " Extension= " + ext_name);
    }

    /* If initially the file did not exist, once the first image is written, then the file exists
        */
    if (!_exists)
        hFile.exist = _exists = true;
}

template <typename T>
int Image<T>::readTiffInMemory(
    void *buf, size_t size, bool readdata, long int select_img,
    bool mapData, bool is_2D
) {
    int err = 0;

    TiffInMemory handle;
    handle.buf = (unsigned char *) buf;
    handle.size = size;
    handle.pos = 0;

    // Check whether to map the data or not
    mmapper = mapData ? new image_mmapper : nullptr;

    // Just clear the header before reading
    header.clear();
    header.addObject();

    TIFF* ftiff = TIFFClientOpen(
        "in-memory-tiff", "r", (thandle_t) &handle,
        TiffInMemoryReadProc,  TiffInMemoryWriteProc, TiffInMemorySeekProc,
        TiffInMemoryCloseProc, TiffInMemorySizeProc,  TiffInMemoryMapFileProc,
        TiffInMemoryUnmapFileProc
    );
    err = readTIFF(ftiff, select_img, readdata, true, "in-memory-tiff");
    TIFFClose(ftiff);

    return err;
}

// Manually instantiate classes derivable from the Image class template.
// Required to avoid linker errors.
// https://www.cs.technion.ac.il/users/yechiel/c++-faq/separate-template-class-defn-from-decl.html
template class Image<int>;
template class Image<float>;
template class Image<double>;
template class Image<signed short>;
template class Image<unsigned short>;
template class Image<signed char>;
template class Image<unsigned char>;
