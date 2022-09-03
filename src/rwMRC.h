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
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * This complete copyright notice must be included in any revised version of the
 * source code. Additional authorship citations may be added, but existing
 * author citations must be preserved.
 ***************************************************************************/
/*
        Based on rwMRC.h
        Header file for reading and writing MRC files
        Format: 3D crystallographic image file format for the MRC package
        Author: Bernard Heymann
        Created: 19990321 Modified: 20030723
*/

#ifndef RWMRC_H
#define RWMRC_H

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "src/image.h"
const int MRCSIZE = 1024; // Minimum size of the MRC header (when nsymbt = 0)

///@defgroup MRC MRC File format
///@ingroup ImageFormats

/** MRC Old Header
  * @ingroup MRC
*/
struct MRCheadold {
    // file header for MRC data
    int nx;              //  0   0       image size
    int ny;              //  1   4
    int nz;              //  2   8
    int mode;            //  3           0=uchar,1=short,2=float
    int nxStart;         //  4           unit cell offset
    int nyStart;         //  5
    int nzStart;         //  6
    int mx;              //  7           unit cell size in voxels
    int my;              //  8
    int mz;              //  9
    float a;             // 10   40      cell dimensions in A
    float b;             // 11
    float c;             // 12
    float alpha;         // 13           cell angles in degrees
    float beta;          // 14
    float gamma;         // 15
    int mapc;            // 16           column axis
    int mapr;            // 17           row axis
    int maps;            // 18           section axis
    float amin;          // 19           minimum density value
    float amax;          // 20   80      maximum density value
    float amean;         // 21           average density value
    int ispg;            // 22           space group number
    int nsymbt;          // 23           bytes used for sym. ops. table
    float extra[29];     // 24           user-defined info
    float xOrigin;       // 53           phase origin in pixels
    float yOrigin;       // 54
    int nlabl;           // 55           number of labels used
    char labels[10][80]; // 56-255       10 80-character labels
};

/** MRC Header
  * @ingroup MRC
*/
struct MRChead {
    // file header for MRC data
    int nx;              //  0   0       image size
    int ny;              //  1   4
    int nz;              //  2   8
    int mode;            //  3           0=char,1=short,2=float
    int nxStart;         //  4           unit cell offset
    int nyStart;         //  5
    int nzStart;         //  6
    int mx;              //  7           unit cell size in voxels
    int my;              //  8
    int mz;              //  9
    float a;             // 10   40      cell dimensions in A
    float b;             // 11
    float c;             // 12
    float alpha;         // 13           cell angles in degrees
    float beta;          // 14
    float gamma;         // 15
    int mapc;            // 16           column axis
    int mapr;            // 17           row axis
    int maps;            // 18           section axis
    float amin;          // 19           minimum density value
    float amax;          // 20   80      maximum density value
    float amean;         // 21           average density value
    int ispg;            // 22           space group number
    int nsymbt;          // 23           bytes used for sym. ops. table
    float extra[25];     // 24           user-defined info
    float xOrigin;       // 49           phase origin in pixels
    float yOrigin;       // 50
    float zOrigin;       // 51
    char map[4];         // 52       identifier for map file ("MAP ")
    char machst[4];      // 53           machine stamp
    float arms;          // 54       RMS deviation
    int nlabl;           // 55           number of labels used
    char labels[800];    // 56-255       10 80-character labels
};

// For determination of machine stamp
enum {
    UNKNOWN_SYSTEM,
    BIGIEEE,
    LITTLEIEEE,
    LITTLEVAX,
};

static int systype() {
    char *test = (char*) askMemory(12);
    int *itest = (int*) test;
    float *ftest = (float*) test;
    memcpy(test, "jbh     ", 8);

    int type = UNKNOWN_SYSTEM;
    if (itest[0] == 1784834080 && fabs(ftest[0] - 6.84272e+25) < 1e+21)
        type = BIGIEEE;
    if (itest[0] == 543711850 && fabs(ftest[0] - 1.96837e-19) < 1e-23)
        type = LITTLEIEEE;

    freeMemory(test, 12);

    return type;
}

static DataType determine_datatype(int mode, int nx, int ny) throw (RelionError) {

    switch (mode) {

        case 12:
        REPORT_ERROR("RELION 3.1 does not support half-precision floating point numbers (MRC mode 12). Please use later versions.");

        case 101:
        // This is SerialEM's non-standard extension.
        // https://bio3d.colorado.edu/imod/doc/mrc_format.txt
        // http://bio3d.colorado.edu/SerialEM/hlp/html/hidd_k2_save_options.htm
        if (nx % 2 == 1 && ny % 2 == 1)
        REPORT_ERROR("Currently we support 4-bit MRC (mode 101) only when nx * ny is even.");
        return UHalf;

    }

    switch (mode % 5) {

        case 0:
        return SChar; // Changed to SIGNED in RELION 3.1 to be compatible with the official specification and SerialEM

        case 1:
        return Short;

        case 2:
        return Float;

        case 3:
        case 4:
        REPORT_ERROR("readMRC: RELION can only read real-space images.");

        default:
        return SChar;

    }
}

// I/O prototypes
/** MRC Reader
  * @ingroup MRC
*/
template <typename T>
int Image<T>::readMRC(long int img_select, bool isStack, const FileName &name) throw (RelionError) {
    #undef DEBUG
    // #define DEBUG
    #ifdef DEBUG
    printf("DEBUG readMRC: Reading MRC file\n");
    #endif

    MRChead *header = (MRChead*) askMemory(sizeof(MRChead));
    if (fread(header, MRCSIZE, 1, fimg) < 1)
        REPORT_ERROR("rwMRC: error in reading header of image " + name);

    // Determine byte order and swap bytes if from little-endian machine
    swap = 0;
    if (abs(header->mode) > SWAPTRIG || abs(header->nx) > SWAPTRIG) {
        #ifdef DEBUG
        fprintf(stderr, "Warning: Swapping header byte order for 4-byte types\n");
        #endif
        swap = 1;
        const int extent = MRCSIZE - 800;  // exclude labels from swapping
        for (int i = 0; i < extent; i += 4) swapbytes((char*) header + i, 4);
    }

    std::array<int, 4> dims {header->nx, header->ny, header->nz, 1};

    if (isStack) {
        std::swap(dims[2], dims[3]);
        replaceNsize = dims[3];
        if (img_select >= dims[3]) {
            // img_select starts from 0, while dims[3] from 1
            REPORT_ERROR((std::string) "readMRC: Image number " + std::to_string(img_select + 1) + " exceeds stack size " + std::to_string(dims[3]) + " of image " + name);
        }
    } else {
        replaceNsize = 0;
    }

    // Map the parameters
    if (isStack) dims[2] = 1;
    if (!isStack || img_select != -1) dims[3] = 1;

    data.setDimensions(dims[0], dims[1], dims[2], dims[3]);

    const DataType datatype = determine_datatype(header->mode, dims[0], dims[1]);

    offset = MRCSIZE + header->nsymbt;

    const long int i = MDMainHeader.index();
    MDMainHeader.setValue(EMDL::IMAGE_STATS_MIN,    (RFLOAT) header->amin, i);
    MDMainHeader.setValue(EMDL::IMAGE_STATS_MAX,    (RFLOAT) header->amax, i);
    MDMainHeader.setValue(EMDL::IMAGE_STATS_AVG,    (RFLOAT) header->amean, i);
    MDMainHeader.setValue(EMDL::IMAGE_STATS_STDDEV, (RFLOAT) header->arms, i);
    MDMainHeader.setValue(EMDL::IMAGE_DATATYPE, (int) datatype, i);

    if (header->mx && header->a != 0)  // ux
        MDMainHeader.setValue(EMDL::IMAGE_SAMPLINGRATE_X, (RFLOAT)header->a / header->mx, i);
    if (header->my && header->b != 0)  // yx
        MDMainHeader.setValue(EMDL::IMAGE_SAMPLINGRATE_Y, (RFLOAT)header->b / header->my, i);
    if (header->mz && header->c != 0)  // zx
        MDMainHeader.setValue(EMDL::IMAGE_SAMPLINGRATE_Z, (RFLOAT)header->c / header->mz, i);

    if (isStack && !dataflag) {
        // Don't read the individual header and the data if not necessary
        freeMemory(header, sizeof(MRChead));
        return 0;
    }

    // #define DEBUG
    #ifdef DEBUG
    MDMainHeader.write(std::cerr);
    MD[0].write(std::cerr);
    #endif

    freeMemory(header, sizeof(MRChead));

    return readData(fimg, img_select, datatype, 0);
}

/** MRC Writer
  * @ingroup MRC
*/
template <typename T>
int Image<T>::writeMRC(long int img_select, bool isStack, int mode) {
    MRChead *header = (MRChead *) askMemory(sizeof(MRChead));

    // Map the parameters
    strncpy(header->map, "MAP ", 4);
    // Set CCP4 machine stamp
    switch (systype()) {
        case BIGIEEE:
        header->machst[0] = header->machst[1] = 17;
        break;
        case LITTLEIEEE:
        header->machst[0] = 68;
        header->machst[1] = 65;
        break;
        case LITTLEVAX:
        header->machst[0] = 34;
        header->machst[1] = 65;
        break;
        case UNKNOWN_SYSTEM:
        REPORT_ERROR("Unkown system type in writeMRC machine stamp determination.");
        default:
        break;
    }

    /// FIXME: TO BE DONE WITH rwCCP4!!
    // set_CCP4_machine_stamp(header->machst);
    const long int Xdim = Xsize(data), Ydim = Ysize(data), 
                   Zdim = Zsize(data), Ndim = Nsize(data);
    long int imgStart = 0, imgEnd = Ndim;
    if (img_select != -1) {
        imgStart = img_select;
        imgEnd = img_select + 1;
    }
    if (mode == WRITE_APPEND || mode == WRITE_REPLACE) {
        imgStart = 0;
        imgEnd = 1;
    }
    header->nx = Xdim;
    header->ny = Ydim;
    header->nz = isStack ? Ndim : Zdim;

    // Convert T to datatype
    DataType output_type;
    if (
        typeid(T) == typeid(unsigned char) ||
        typeid(T) == typeid(signed char)
    ) {
        header->mode = 0;
        output_type = UChar;
    } else if (typeid(T) == typeid(short)) {
        header->mode = 1;
        output_type = Short;
    } else if (
        typeid(T) == typeid(RFLOAT) ||
        typeid(T) == typeid(float)  ||
        typeid(T) == typeid(int)
    ) {
        header->mode = 2;
        output_type = Float;
    } else {
        REPORT_ERROR("ERROR write MRC image: invalid typeid(T)");
    }

    // Set this to zero till we decide if we want to update it
    header->mx = header->nx;  // (int) (ua/ux + 0.5);
    header->my = header->ny;  // (int) (ub/uy + 0.5);
    header->mz = header->nz;  // (int) (uc/uz + 0.5);
    header->mapc = 1;
    header->mapr = 2;
    header->maps = 3;

    // TODO: fix this!
    header->a = (float) 0.0;  // ua;
    header->b = (float) 0.0;  // ub;
    header->c = (float) 0.0;  // uc;
    header->alpha = (float) 90.0;
    header->beta  = (float) 90.0;
    header->gamma = (float) 90.0;
    header->xOrigin = (float) 0.0;
    header->yOrigin = (float) 0.0;
    header->zOrigin = (float) 0.0;
    header->nxStart = (int) 0;
    header->nyStart = (int) 0;
    header->nzStart = (int) 0;

    if (!MDMainHeader.empty()) {

        const long int i = MDMainHeader.index();

        // If MDMainHeader contains none of these,
        // we will be looping through `data` more times than necessary.
        header->amin  = MDMainHeader.containsLabel(EMDL::IMAGE_STATS_MIN)    ? MDMainHeader.getValue<float>(EMDL::IMAGE_STATS_MIN,    i) : (float) min(data);
        header->amax  = MDMainHeader.containsLabel(EMDL::IMAGE_STATS_MAX)    ? MDMainHeader.getValue<float>(EMDL::IMAGE_STATS_MAX,    i) : (float) max(data);
        header->amean = MDMainHeader.containsLabel(EMDL::IMAGE_STATS_AVG)    ? MDMainHeader.getValue<float>(EMDL::IMAGE_STATS_AVG,    i) : (float) average(data);
        header->arms  = MDMainHeader.containsLabel(EMDL::IMAGE_STATS_STDDEV) ? MDMainHeader.getValue<float>(EMDL::IMAGE_STATS_STDDEV, i) : (float) computeStddev(data);

        // int nxStart, nyStart, nzStart;

        // RFLOAT aux;
        // aux = MDMainHeader.getValue<float>(EMDL::ORIENT_ORIGIN_X, i);
        // if (std::isfinite(nxStart = aux - 0.5)) { header->nxStart = nxStart; }

        if (MDMainHeader.containsLabel(EMDL::IMAGE_SAMPLINGRATE_X)) {
        const float srx = MDMainHeader.getValue<float>(EMDL::IMAGE_SAMPLINGRATE_X, i),
                    xOrigin = header->nxStart * srx, a = header->nx * srx;
        if (std::isfinite(xOrigin)) { header->xOrigin = xOrigin; }
        if (std::isfinite(a))       { header->a       = a; }
        }

        // aux = MDMainHeader.getValue<float>(EMDL::ORIENT_ORIGIN_Y, i, aux);
        // if (std::isfinite(nyStart = aux - 0.5)) { header->nyStart = nyStart; }

        if (MDMainHeader.containsLabel(EMDL::IMAGE_SAMPLINGRATE_Y))  {
        const float sry = MDMainHeader.getValue<float>(EMDL::IMAGE_SAMPLINGRATE_Y, i),
                    yOrigin = header->nyStart * sry, b = header->ny * sry;
        if (std::isfinite(yOrigin)) { header->yOrigin = yOrigin; }
        if (std::isfinite(b))       { header->b       = b; }
        }

        // aux = MDMainHeader.getValue<float>(EMDL::ORIENT_ORIGIN_Z, i);
        // if (std::isfinite(nzStart = aux - 0.5)) { header->nzStart = nzStart; }

        if (MDMainHeader.containsLabel(EMDL::IMAGE_SAMPLINGRATE_Z))  {
        const float srz = MDMainHeader.getValue<float>(EMDL::IMAGE_SAMPLINGRATE_Z, i),
                    zOrigin = header->nzStart * srz, c = header->nz * srz;
        if (std::isfinite(zOrigin)) { header->zOrigin = zOrigin; }
        if (std::isfinite(c))       { header->c       = c; }
        }

    }

    header->nsymbt = 0;

    // Create label "Relion version    date time"
    const int MRC_LABEL_LEN = 80;
    header->nlabl = 1;

    char label[MRC_LABEL_LEN] = "Relion ";
    time_t rawtime;
    time(&rawtime);
    struct tm *timeinfo = localtime(&rawtime);

    #ifdef PACKAGE_VERSION
    strcat(label, PACKAGE_VERSION);
    #endif
    strcat(label, "   ");
    strftime(label + strlen(label), MRC_LABEL_LEN - strlen(label), "%d-%b-%y  %R:%S", timeinfo);
    strncpy(header->labels, label, MRC_LABEL_LEN);

    // strncpy(header->labels, p->label.c_str(), 799);

    offset = MRCSIZE + header->nsymbt;
    const size_t datasize_n = Xdim * Ydim * Zdim;
    const size_t datasize = datasize_n * gettypesize(output_type);

    // #define DEBUG
    #ifdef DEBUG
    printf("DEBUG rwMRC: Offset = %ld,  Datasize_n = %ld\n", offset, datasize_n);
    #endif

    // For multi-image files
    if (mode == WRITE_APPEND && isStack) {
        header->nz = replaceNsize + 1;
    }
    //else header-> is correct

    //locking
    struct flock fl;

    fl.l_type   = F_WRLCK;  /* F_RDLCK, F_WRLCK, F_UNLCK    */
    fl.l_whence = SEEK_SET; /* SEEK_SET, SEEK_CUR, SEEK_END */
    fl.l_start  = 0;        /* Offset from l_whence         */
    fl.l_len    = 0;        /* length, 0 = to EOF           */
    fl.l_pid    = getpid(); /* our PID                      */

    //BLOCK
    fl.l_type   = F_WRLCK;
    fcntl(fileno(fimg), F_SETLKW, &fl); /* locked */

    // Write header
    if (mode == WRITE_OVERWRITE || mode == WRITE_APPEND)
        fwrite(header, MRCSIZE, 1, fimg);
    freeMemory(header, sizeof(MRChead));

    // write only once, ignore select_img
    char* fdata = (char*) askMemory(datasize);
    // think about writing in several chunks

    if (Nsize(data) == 1 && mode == WRITE_OVERWRITE) {
        castPage2Datatype(fdata, data.data, output_type, datasize_n);
        fwrite(fdata, datasize, 1, fimg);
    } else {
        if (mode == WRITE_APPEND) {
            fseek(fimg, 0, SEEK_END);
        } else if (mode == WRITE_REPLACE) {
            fseek(fimg, offset + datasize * img_select, SEEK_SET);
        }

        for (size_t i = imgStart; i < imgEnd; i++) {
            castPage2Datatype(fdata, data.data + i * datasize_n, output_type, datasize_n);
            fwrite(fdata, datasize, 1, fimg);
        }
    }

    // Unlock the file
    fl.l_type = F_UNLCK;
    fcntl(fileno(fimg), F_SETLK, &fl); /* unlocked */

    freeMemory(fdata, datasize);

    return 0;
}
#endif
