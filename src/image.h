/***************************************************************************
 *
 * Author: "Sjors H.W. Scheres", "Takanori Nakane"
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
 * Authors: Sjors H.W. Scheres (scheres@cnb.csic.es)
 *
 * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
 *
 * Part of this module has been developed by Lorenzo Zampighi and Nelson Tang
 * Dept. Physiology of the David Geffen School of Medicine
 * Univ. of California, Los Angeles.
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
 * All comments concerning this program package may be sent to the
 * e-mail address 'xmipp@cnb.csic.es'
 ***************************************************************************/

#ifndef IMAGE_H
#define IMAGE_H

#include <typeinfo>
#include <fcntl.h>
#include <sys/mman.h>
#include <fcntl.h>
#include <unistd.h>
#include <tiffio.h>
#include "src/funcs.h"
#include "src/memory.h"
#include "src/filename.h"
#include "src/multidim_array.h"
#include "src/transformations.h"
#include "src/metadata_table.h"
#include "src/fftw.h"

/// @defgroup Images Images
//@{

// An ad-hoc representation of the type system
// Allows for image data type inspection/manipulation at run time
enum DataType {
    Unknown_Type,   // Undefined data type
    UChar,          // Unsigned character or byte type
    SChar,          // Signed character (for CCP4)
    UShort,         // Unsigned integer (2-byte)
    Short,          // Signed integer (2-byte)
    UInt,           // Unsigned integer (4-byte)
    Int,            // Signed integer (4-byte)
    Long,           // Signed integer (4 or 8 byte, depending on system)
    Float,          // Floating point (4-byte)
    Double,         // Double precision floating point (8-byte)
    Boolean,        // Boolean (1-byte?)
    UHalf,          // Signed 4-bit integer (SerialEM extension)
};

// Convert string to int corresponding to value in enum
// int DataType::String2Int(std::string s);

/** Returns memory size of datatype
 */
unsigned long gettypesize(DataType type);

/** WriteMode
 * To indicate the writing behavior.
 */
enum WriteMode {
    WRITE_OVERWRITE, // Forget about the old file and overwrite it
    WRITE_APPEND,	 // Append and object at the end of a stack, so far can not append stacks
    WRITE_REPLACE,	 // Replace a particular object with another
    WRITE_READONLY	 // Read-only
};


static std::string writemode2string(WriteMode mode, bool exist) {

    switch (mode) {

        case WRITE_READONLY:
        return "r";

        case WRITE_OVERWRITE:
        return "w";

        case WRITE_APPEND:
        return exist ? "r+" : "w+";  // w+ will destroy file contents. We don't want that.

        case WRITE_REPLACE:
        return "r+";

        default:
        throw "Invalid write mode!";

    }

}

extern "C" {

    typedef struct TiffInMemory {
        unsigned char *buf;
        tsize_t size;
        toff_t pos;
    } TiffInMemory;

    static tsize_t TiffInMemoryReadProc(thandle_t handle, tdata_t buf, tsize_t read_size) {
        TiffInMemory *tiff_handle = (TiffInMemory *) handle;
        #ifdef TIFF_DEBUG
        std::cout << "TiffInMemoryReadProc: read_size = " << read_size << " cur_pos = " << tiff_handle->pos << " buf_size = " << tiff_handle->size << std::endl;
        #endif
        if (tiff_handle->pos + read_size >= tiff_handle->size)
            REPORT_ERROR("TiffInMemoryReadProc: seeking beyond the end of the buffer.");

        memcpy(buf, tiff_handle->buf + tiff_handle->pos, read_size);
        tiff_handle->pos += read_size;

        return read_size;
    }

    static tsize_t TiffInMemoryWriteProc(thandle_t handle, tdata_t buf, tsize_t write_size) {
        #ifdef TIFF_DEBUG
        REPORT_ERROR("TiffInMemoryWriteProc: Not implemented.");
        #endif
        return -1;
    }

    static toff_t TiffInMemorySeekProc(thandle_t handle, toff_t offset, int whence) {
        TiffInMemory *tiff_handle = (TiffInMemory*) handle;
        #ifdef TIFF_DEBUG
        std::cout << "TiffInMemorySeekProc: offset = " << offset << " cur_pos = " << tiff_handle->pos << " buf_size = " << tiff_handle->size << std::endl;
        #endif
        switch (whence) {

            case SEEK_SET:
            tiff_handle->pos = 0;
            break;

            case SEEK_CUR:
            tiff_handle->pos += offset;
            break;

            case SEEK_END:
            REPORT_ERROR("TIFFInMemorySeekProc: SEEK_END is not supported.");
            // break; // intentional to suppress compiler warnings.

        }

        if (tiff_handle->pos >= tiff_handle->size)
            REPORT_ERROR("TIFFInMemorySeekProc: seeking beyond the end of the buffer.");

        return 0;
    }

    static int TiffInMemoryCloseProc(thandle_t handle) {
        #ifdef TIFF_DEBUG
        std::cout << __func__ << std::endl;
        #endif
        return 0;
    }

    static toff_t TiffInMemorySizeProc(thandle_t handle) {
        #ifdef TIFF_DEBUG
        std::cout << __func__ << std::endl;
        #endif
        return ((TiffInMemory *) handle)->size;
    }

    static int TiffInMemoryMapFileProc(thandle_t handle, tdata_t *base, toff_t *size) {
        TiffInMemory *tiff_handle = (TiffInMemory *) handle;
        #ifdef TIFF_DEBUG
        std::cout << __func__ << std::endl;
        #endif

        *base = tiff_handle->buf;
        *size = tiff_handle->size;

        return 1;
    }

    static void TiffInMemoryUnmapFileProc(thandle_t handle, tdata_t base, toff_t size) {
        #ifdef TIFF_DEBUG
            std::cout << __func__ << std::endl;
        #endif

        return;
    }
}

/** File handler class
 * This struct is used to share the File handlers with Image Collection class
 */
class fImageHandler {

    public:

    FILE     *fimg;	// Image File handler
    FILE     *fhed;	// Image File header handler
    TIFF     *ftiff;
    FileName  ext_name; // Filename extension
    bool	  exist;    // Does the file exist?

    // Empty constructor
    fImageHandler() {
        fimg = NULL;
        fhed = NULL;
        ftiff = NULL;
        ext_name = "";
        exist = false;
    }

    ~fImageHandler() { closeFile(); }

    void openFile(const FileName &name, int mode = WRITE_READONLY) {

        // Close any file that was left open in this handler
        if (fimg || fhed) closeFile();

        FileName fileName, headName = "";
        // get the format, checking for possible format specifier before suffix
        // getFileFormat("file.spi") will return "spi"
        // getFileFormat("file.spi:mrc") will return "mrc"
        // getFileFormat("file") will return ""
        ext_name = name.getFileFormat();

        long int dump;
        name.decompose(dump, fileName);
        // Subtract 1 to have numbering 0...N-1 instead of 1...N
        if (dump > 0) { dump--; }

        // create the filename from a possible input format specifier (file.spi:mrc means "it's called .spi, but it's really a .mrc")
        // file.spi:mrc -> file.spi
        fileName = fileName.removeFileFormat();

        size_t found = fileName.find_first_of("%");
        if (found != std::string::npos) { fileName = fileName.substr(0, found); }

        exist = exists(fileName);

        if (mode == WRITE_READONLY and !exist)
        REPORT_ERROR((std::string) "Can't read file " + fileName + ". It doesn't exist!");

        std::string wmstr = writemode2string((WriteMode) mode, exist); // Write mode string

        if (ext_name.contains("img") || ext_name.contains("hed")) {
            fileName = fileName.withoutExtension();
            headName = fileName.addExtension("hed");
            fileName = fileName.addExtension("img");
        } else if (ext_name == "") {
            ext_name = "spi"; // SPIDER is default format
            fileName = fileName.addExtension(ext_name);
        }

        bool isTiff = ext_name.contains("tif");
        if (isTiff && mode != WRITE_READONLY)
            REPORT_ERROR((std::string) "TIFF is supported only for reading");

        // Open image file
        if (
             isTiff && !(ftiff = TIFFOpen(fileName.c_str(), "r")) ||
            !isTiff && !(fimg  = fopen   (fileName.c_str(), wmstr.c_str()))
        ) {
            REPORT_ERROR((std::string) "Image::" + __func__ + " cannot open: " + name);
        }

        if (headName != "") {
            if (!(fhed = fopen(headName.c_str(), wmstr.c_str())))
                REPORT_ERROR((std::string) "Image::" + __func__ + " cannot open: " + headName);
        } else {
            fhed = NULL;
        }

    }

    // Close file (if open)
    void closeFile() {
        ext_name = "";
        exist = false;

        // Check whether the file was closed already
        if (!fimg && !fhed && !ftiff) return;

        bool isTiff = ext_name.contains("tif");
        if (isTiff && ftiff) {
            TIFFClose(ftiff);
            ftiff = NULL;
        }

        if (!isTiff && fclose(fimg) != 0) {
            REPORT_ERROR((std::string) "Cannot close image file ");
        } else {
            fimg = NULL;
        }

        if (fhed && fclose(fhed) != 0) {
            REPORT_ERROR((std::string) "Cannot close header file ");
        } else {
            fhed = NULL;
        }
    }

};

/** Swapping trigger.
 * Threshold file z size above which bytes are swapped.
 */
#define SWAPTRIG 65535

/** Template class for images.
 * The image class is the general image handling class.
 */
template<typename T>
class Image {

    public:

    MultidimArray<T> data;      // Image data
    MetaDataTable MDMainHeader; // File metadata

    // Why int x, y, z but long int n?
    struct Dimensions { int x, y, z; long int n; };

    private:

    FileName filename; // File name
    FILE *fimg; // Image File handler
    FILE *fhed; // Image File header handler
    bool stayOpen; // To maintain the image file open after read/write
    bool dataflag;  // Control reading of the data
    unsigned long i; // Current image number (may be > Nsize)
    unsigned long offset; // Data offset
    int swap; // Perform byte swapping upon reading
    long int replaceNsize; // Stack size in the replace case
    bool _exists;  // Does the target file exist? 0 if file does not exist or is not a stack.

    // Allocation-related member variables

    bool mmapOn; // Whether to alloc memory or map to a file
    int mFd; // Mapped file handle
    size_t mappedSize; // Size of the mapped file

    void attempt_mmap(T *data, FileName &mapFile, int &mFd, size_t pagesize) {

        // if ((mFd = open(mapFile.c_str(), O_RDWR, S_IREAD | S_IWRITE)) == -1)
        if ((mFd = open(mapFile.c_str(), O_RDWR, S_IRUSR | S_IWUSR)) == -1)
            REPORT_ERROR((std::string) "Image<T>::" + __func__ + ": Error opening the image file.");

        mappedSize = pagesize + offset;

        char *map;
        if ((map = (char*) mmap(0, mappedSize, PROT_READ | PROT_WRITE, MAP_SHARED, mFd, 0)) == (void*) -1)
            REPORT_ERROR((std::string) "Image<T>::" + __func__ + ": mmap of image file failed.");
        data = reinterpret_cast<T*> (map + offset);

    }

    public:

    /** Empty constructor
     *
     * An empty image is created.
     *
     * @code
     * Image<RFLOAT> I;
     * @endcode
     */
    Image() {
        mmapOn = false;
        clear();
        MDMainHeader.addObject();
    }

    /** Constructor with size
     *
     * A blank image (0.0 filled) is created with the given size. Pay attention
     * to the dimension order: Y and then X.
     *
     * @code
     * Image I(64,64);
     * @endcode
     */
    Image(long int Xdim, long int Ydim, long int Zdim = 1, long int Ndim = 1) {
        mmapOn = false;
        clear();
        data.resize(Ndim, Zdim, Ydim, Xdim);
        MDMainHeader.addObject();
    }

    inline static Image<T> zeros(long int Xdim, long int Ydim, long int Zdim = 1, long int Ndim = 1) {
        Image<T> img(Xdim, Ydim, Zdim, Ndim);
        img.data.initZeros();
    }

    /** Clear.
     * Zero-initialize.
     */
    void clear() {
        if (mmapOn) {
            munmap(data.data - offset, mappedSize);
            close(mFd);
            data.data = NULL;
        } else {
            data.clear();
        }

        dataflag = false;
        i = 0;
        filename = "";
        offset = 0;
        swap = 0;
        clearHeader();
        replaceNsize = 0;
        mmapOn = false;
    }

    // Clear the image header
    void clearHeader() { MDMainHeader.clear(); }

    ~Image() { clear(); }

    // Read/write functions for different file formats

    int readSPIDER(long int img_select);

    int writeSPIDER(long int select_img=-1, bool isStack = false, int mode = WRITE_OVERWRITE);

    int readMRC(long int img_select, bool isStack = false, const FileName &name = "");

    int writeMRC(long int img_select, bool isStack = false, int mode = WRITE_OVERWRITE);

    int readIMAGIC(long int img_select);

    void writeIMAGIC(long int img_select = -1, int mode = WRITE_OVERWRITE);

    int readTIFF(
        TIFF *ftiff, long int img_select, 
        bool readdata = false, bool isStack = false, const FileName &name = ""
    );

    /** Is this file an image?
     *
     *	Check whether a real-space image can be read.
     */
    bool isImage(const FileName &name) { return !read(name, false); }

    // Rename the image
    void rename(const FileName &name) { filename = name; }

    /** General read function
     * you can read a single image from a single image file
     * or a single image file from an stack, in the second case
     * the select slide may come in the image name or in the select_img parameter
     * file name takes precedence over select_img
     * If -1 is given the whole object is read
     * The number before @ in the filename is 1-indexed, while select_img is 0-indexed.
     */
    int read(
        const FileName &name, bool readdata = true, long int select_img = -1, 
        bool mapData = false, bool is_2D = false
    );

    /** Read from an open file
     */
    int readFromOpenFile(
        const FileName &name, fImageHandler &hFile, long int select_img, 
        bool is_2D = false
    ) {
        int err = _read(name, hFile, true, select_img, false, is_2D);
        // Reposition file pointer for a next read
        rewind(fimg);
        return err;
    }

    /** General write function
     * select_img = which slice should I replace
     * overwrite = 0, append slice
     * overwrite = 1 overwrite slice
     *
     * NOTE:
     *	select_img has higher priority than the number before "@" in the name.
     *	select_img counts from 0, while the number before "@" in the name from 1!
     */
    void write(
        FileName name = "", long int select_img = -1, bool isStack = false,
        int mode = WRITE_OVERWRITE
    );

    /** Cast a page of data from type dataType to type Tdest
     *	  input pointer char *
     */
    void castPage2T(char *page, T *ptrDest, DataType datatype, size_t pageSize) {
        switch (datatype) {

            case Unknown_Type:
            REPORT_ERROR("ERROR: datatype is Unknown_Type");

            case UChar:
            cast_page_from<unsigned char>(page, ptrDest, pageSize);
            break;

            case SChar:
            cast_page_from<signed char>(page, ptrDest, pageSize);
            break;

            case UShort:
            cast_page_from<unsigned short>(page, ptrDest, pageSize);
            break;

            case Short:
            cast_page_from<short>(page, ptrDest, pageSize);
            break;

            case UInt:
            cast_page_from<unsigned int>(page, ptrDest, pageSize);
            break;

            case Int:
            cast_page_from<int>(page, ptrDest, pageSize);
            break;

            case Long:
            cast_page_from<long>(page, ptrDest, pageSize);
            break;

            case Float:
            cast_page_from<float>(page, ptrDest, pageSize);
            break;

            case Double:
            cast_page_from<RFLOAT>(page, ptrDest, pageSize);
            break;

            case UHalf:
            if (pageSize % 2 != 0) {
                REPORT_ERROR((std::string) "Logic error in " + __func__ + "; for UHalf, pageSize must be even.");
            }

            for (size_t i = 0, ilim = pageSize / 2; i < ilim; i++) {
                // Here we are assuming the fill-order is LSB2MSB according to IMOD's
                // iiProcessReadLine() in libiimod/mrcsec.c.
                // The default fill-order in the TIFF specification is MSB2LSB
                // but IMOD assumes LSB2MSB even for TIFF.
                // See IMOD's iiTIFFCheck() in libiimod/iitif.c.
                ptrDest[i * 2    ] = (T) ( page[i]       & 15); // 1111 = 1+2+4+8 = 15
                ptrDest[i * 2 + 1] = (T) ((page[i] >> 4) & 15);
            }
            break;

            default:
            std::cerr << "Datatype= " << datatype << std::endl;
            REPORT_ERROR(" ERROR: cannot cast datatype to T");

        }
    }

    template <typename T2>
    void cast_page_from(char *page, T *dest, size_t pageSize) {
        if (typeid(T) == typeid(T2)) {
            memcpy(dest, page, pageSize * sizeof(T));
        } else {
            T2 *src = (T2 *) page;
            for (size_t i = 0; i < pageSize; i++) {
                dest[i] = (T) src[i];
            }
        }
    }

    /** Cast page from T to datatype
     *  input pointer char *
     */
    void castPage2Datatype(char *page, T *srcPtr, DataType datatype, size_t pageSize) {
        switch (datatype) {

            case Float:
            cast_page_to<float>(page, srcPtr, pageSize);
            break;

            case Double:
            cast_page_to<RFLOAT>(page, srcPtr, pageSize);
            break;

            case Short: 
            cast_page_to<short>(page, srcPtr, pageSize);
            break;

            case UShort:
            cast_page_to<unsigned short>(page, srcPtr, pageSize);
            break;

            case UChar:
            cast_page_to<unsigned char>(page, srcPtr, pageSize);
            break;

            default:
            std::cerr << "outputDatatype= " << datatype << std::endl;
            REPORT_ERROR(" ERROR: cannot cast T to outputDatatype");

            }
    }

    template <typename T2>
    void cast_page_to(char *page, T *src, size_t pageSize) {
        if (typeid(T) == typeid(T2)) {
            memcpy(page, src, pageSize * sizeof(T));
        } else {
            T2 *dest = (T2*) page;
            for (size_t i = 0; i < pageSize; i++) {
                dest[i] = (T2) src[i];
            }
        }
    }

    // Check file Datatype is same as T type to use mmap.
    bool checkMmapT(DataType datatype) {

        switch (datatype) {
            case Unknown_Type:
            REPORT_ERROR("ERROR: datatype is Unknown_Type");
            case UChar:
            return typeid(T) == typeid(unsigned char);
            case SChar:
            return typeid(T) == typeid(signed char);
            case UShort:
            return typeid(T) == typeid(unsigned short);
            case Short:
            return typeid(T) == typeid(short);
            case UInt:
            return typeid(T) == typeid(unsigned int);
            case Int:
            return typeid(T) == typeid(int);
            case Long:
            return typeid(T) == typeid(long);
            case Float:
            return typeid(T) == typeid(float);
            case Double:
            return typeid(T) == typeid(RFLOAT);
            default:
            std::cerr << "Datatype= " << datatype << std::endl;
            REPORT_ERROR(" ERROR: cannot cast datatype to T");
        }
        //  int *iTemp = (int*) map;
        //  ptrDest = reinterpret_cast<T*> (iTemp);
    }

    /** Write an entire page as datatype
     *
     * A page of datasize_n elements T is cast to datatype and written to fimg
     * The memory for the casted page is allocated and freed internally.
     */
    void writePageAsDatatype(FILE *fimg, DataType datatype, size_t datasize_n ) {
        size_t datasize = datasize_n * gettypesize(datatype);
        char * fdata = (char *) askMemory(datasize);
        castPage2Datatype(fdata, data.data, datatype, datasize_n);
        fwrite(fdata, datasize, 1, fimg);
        freeMemory(fdata, datasize);
    }

    /** Swap an entire page
      * input pointer char *
      */
    void swapPage(char *page, size_t pageNrElements, DataType datatype) {
        unsigned long datatypesize = gettypesize(datatype);
        #ifdef DEBUG
            std::cerr << "DEBUG " << __func__ << ": Swapping image data with swap= "
            << swap << " datatypesize= " << datatypesize
            << " pageNrElements " << pageNrElements
            << " datatype " << datatype
            << std::endl;
        #endif

        // Swap bytes if required
        if (swap >= 1) {
            unsigned long increment = swap == 1 ? datatypesize : swap;
            for (unsigned long i = 0; i < pageNrElements; i += increment)
                swapbytes(page + i, increment);
        }
    }

    // Read the raw data
    int readData(FILE *fimg, long int select_img, DataType datatype, unsigned long pad) {
        // #define DEBUG
        #ifdef DEBUG
        std::cerr << "entering " << __func__ << std::endl;
        std::cerr << " " << __func__ << " flag= " << dataflag << std::endl;
        #endif

        if (!dataflag) return 0;

        size_t datatypesize; // bytes
        size_t pagesize; // bytes
        if (datatype == UHalf) {
            if (Xsize(data) * Ysize(data) % 2 != 0) {
                REPORT_ERROR("For UHalf, Xsize(data) * Ysize(data) must be even.");
            }
            /// BUG: datatypesize not set
            pagesize = Xsize(data) * Ysize(data) * Zsize(data) / 2;
        } else {
            datatypesize = gettypesize(datatype);
            pagesize = Xsize(data) * Ysize(data) * Zsize(data) * datatypesize;
        }
        size_t haveread_n = 0; // number of pixels (not necessarily bytes!) processed so far

        // data.mmapOn takes priority over this->mmapOn
        if (data.getMmap() || datatype == UHalf) { mmapOn = false; }

        // Flag to know that data is not going to be mapped although mmapOn is true
        if (mmapOn && !checkMmapT(datatype)) {
            std::cout << "WARNING: Image Class. File datatype and image declaration not compatible with mmap. Loading into memory." << std::endl;
            mmapOn = false;
            mFd = -1;
        }

        if (mmapOn) {
            if (Nsize(data) > 1) {
                REPORT_ERROR(
                    (std::string) "Image<T>::" + __func__ + ": mmap with multiple \
                    images file not compatible. Try selecting a unique image."
                );
            }
            fclose(fimg);
            attempt_mmap(data.data, filename, mFd, pagesize);
        } else {
            // Reset select to get the correct offset
            if (select_img < 0) { select_img = 0; }

            // Allocate memory for image data 
            // (Assume xdim, ydim, zdim and ndim are already set
            // if memory already allocated use it (no resize allowed)
            data.coreAllocateReuse();
            size_t myoffset = offset + select_img * (pagesize + pad);
            // #define DEBUG

            #ifdef DEBUG
            data.printShape();
            printf("DEBUG: Page size: %ld offset= %d \n", pagesize, offset);
            printf("DEBUG: Swap = %d  Pad = %ld  Offset = %ld\n", swap, pad, offset);
            printf("DEBUG: myoffset = %d select_img= %d \n", myoffset, select_img);
            #endif

            const size_t pagemax = 1073741824; // 1 GB
            char *page = askMemory(std::max(pagesize, pagemax) * sizeof(char));

            // Because we requested XYsize to be even for UHalf, this is always safe.
            if (fseek(fimg, myoffset, SEEK_SET) != 0) return -1;

            for (size_t n = 0; n < Nsize(data); n++) {
                for (size_t j = 0; j < pagesize; j += pagemax) {
                    // pagesize size of object
                    // Read next page. Divide pages larger than pagemax
                    size_t readsize = pagesize - j;
                    if (readsize > pagemax) { readsize = pagemax; }

                    size_t readsize_n = datatype == UHalf ? readsize * 2 : readsize / datatypesize;

                    #ifdef DEBUG
                    std::cout << "NX = " << Xsize(data) << " NY = " << Ysize(data) << " NZ = " << Zsize(data) << std::endl;
                    std::cout << "pagemax = " << pagemax << " pagesize = " << pagesize  << " readsize = " << readsize << " readsize_n = " << readsize_n << std::endl;
                    #endif

                    // Read page from disc
                    if (fread(page, readsize, 1, fimg) != 1) return -2;

                    // swap per page
                    if (swap) swapPage(page, readsize, datatype);
                    // cast to T per page
                    castPage2T(page, data.data + haveread_n, datatype, readsize_n);
                    haveread_n += readsize_n;
                }
                if (pad > 0) {
                    // fread(padpage, pad, 1, fimg);
                    if (fseek(fimg, pad, SEEK_CUR) != 0) return -1;
                }
            }
            // if (pad > 0) { freeMemory(padpage, pad * sizeof(char)); }
            if (page) { freeMemory(page, pagesize * sizeof(char)); }

        #ifdef DEBUG
        printf("DEBUG img_read_data: Finished reading and converting data\n");
        #endif

        }
        return 0;
    }

    /** Data access
     *
     * This operator can be used to access the data multidimarray.
     * In this way 
     * we could resize an image just by resizing its associated matrix:
     * @code
     * image().resize(128, 128);
     * @endcode
     * or we could add two images by adding their matrices.
     * @code
     * image1() = image2() + image3();
     * @endcode
     */
    MultidimArray<T>& operator () () {
        /// NOTE: [jhooker] Is this the most intuitive way to acces data?
        // Why not just make data public?
        return data;
    }

    const MultidimArray<T>& operator () () const {
        return data;
    }

    /** Pixel access
    *
    * This operator is used to access a pixel within a 2D image. This is a
    * logical access, so you could access to negative positions if the image
    * has been defined so (see the general explanation for the class).
    *
    * @code
    * std::cout << "Grey level of pixel (-3,-3) of the image = " << I(-3, -3)
    * << std::endl;
    *
    * I(-3, -3) = I(-3, -2);
    * @endcode
    */
    T& operator () (int i, int j) const {
        // return data.elem(i, j);
        return data.data[(j - data.xinit) + data.xdim * (i - data.yinit)];
    }

    #ifdef IMGPIXEL
    /** Set pixel
     * (direct access) needed by swig
     */
    void setPixel(int i, int j, T v) {
        IMGPIXEL(*this, i, j) = v;
    }

    /** Get pixel
     * (direct acces) needed by swig
     */
    T getPixel(int i, int j) const {
        return IMGPIXEL(*this, i, j);
    }
    #endif

    /** Voxel access
     *
     * This operator is used to access a voxel within a 3D image. This is a
     * logical access, so you could access to negative positions if the image
     * has been defined so (see the general explanation for the class).
     *
     * @code
     * std::cout << "Grey level of pixel (-3,-3, 1) of the volume = " << I(-3, -3, 1)
     * << std::endl;
     *
     * I(-3, -3, 1) = I(-3, -2, 0);
     * @endcode
     */
    inline T& operator () (int k, int i, int j) const {
        return data.elem(i, j, k);
    }

    /** const reference to filename
     *
     * @code
     * std::cout << "Image name = " << I.name() << std::endl;
     * @endcode
     */
    const FileName& name() const {
        return filename;
    }

    /** Get Image dimensions
     */
    Dimensions getDimensions() const {
        Dimensions dimensions;
        dimensions.x = Xsize(data);
        dimensions.y = Ysize(data);
        dimensions.z = Zsize(data);
        dimensions.n = Nsize(data);
        return dimensions;
    }

    long unsigned int getSize() const {
        return data.size();
    }

    /* Is there label in the main header */
    bool mainContainsLabel(EMDL::EMDLabel label) const {
        return MDMainHeader.containsLabel(label);
    }

    /** Data type
     *
     * @code
     * std::cout << "datatype= " << dataType() << std::endl;
     * @endcode
     */
    int dataType() const {
        return MDMainHeader.getValue<int>(EMDL::IMAGE_DATATYPE);
    }

    /** Sampling rate in X
    *
    * @code
    * std::cout << "sampling= " << samplingRateX() << std::endl;
    * @endcode
    */
    RFLOAT samplingRateX(const long int n = 0) const {
        try {
            return MDMainHeader.getValue<RFLOAT>(EMDL::IMAGE_SAMPLINGRATE_X);
        } catch (const char *errmsg) {
            return 1.0;
        }
    }

    /** Sampling rate in Y
    *
    * @code
    * std::cout << "sampling= " << samplingRateY() << std::endl;
    * @endcode
    */
    RFLOAT samplingRateY(const long int n = 0) const {
        try {
            return MDMainHeader.getValue<RFLOAT>(EMDL::IMAGE_SAMPLINGRATE_Y);
        } catch (const char *errmsg) {
            return 1.0;
        }
    }

    // Set file name
    void setName(const FileName &_filename) {
        filename = _filename;
    }

    // Set image statistics in the main header
    void setStatisticsInHeader() {
        Stats<T> statistics = data.computeStats();
        MDMainHeader.setValue(EMDL::IMAGE_STATS_AVG,    statistics.avg);
        MDMainHeader.setValue(EMDL::IMAGE_STATS_STDDEV, statistics.stddev);
        MDMainHeader.setValue(EMDL::IMAGE_STATS_MIN,    statistics.min);
        MDMainHeader.setValue(EMDL::IMAGE_STATS_MAX,    statistics.max);
    }

    void setSamplingRateInHeader(RFLOAT rate_x, RFLOAT rate_y = -1.0, RFLOAT rate_z = -1.0) {
        MDMainHeader.setValue(EMDL::IMAGE_SAMPLINGRATE_X, rate_x);
        if (rate_y < 0.0) rate_y = rate_x;
        MDMainHeader.setValue(EMDL::IMAGE_SAMPLINGRATE_Y, rate_y);
        if (Zsize(data) > 1) {
            if (rate_z < 0.0) rate_z = rate_x;
            MDMainHeader.setValue(EMDL::IMAGE_SAMPLINGRATE_Z, rate_z);
        }
    }

    // Show image properties
    friend std::ostream& operator << (std::ostream& o, const Image<T>& I) {
        o << "Image type   : ";
        o << "Real-space image" << std::endl;

        o << "Reversed	   : ";
        o << (I.swap ? "TRUE" : "FALSE") << std::endl;

        o << "Data type    : ";
        switch (I.dataType()) {
            case Unknown_Type:
            o << "Undefined data type";
            break;

            case UChar:
            o << "Unsigned character or byte type";
            break;

            case SChar:
            o << "Signed character (for CCP4)";
            break;

            case UShort:
            o << "Unsigned integer (2-byte)";
            break;

            case Short:
            o << "Signed integer (2-byte)";
            break;

            case UInt:
            o << "Unsigned integer (4-byte)";
            break;

            case Int:
            o << "Signed integer (4-byte)";
            break;

            case Long:
            o << "Signed integer (4 or 8 byte, depending on system)";
            break;

            case Float:
            o << "Floating point (4-byte)";
            break;

            case Double:
            o << "Double precision floating point (8-byte)";
            break;

            case Boolean:
            o << "Boolean (1-byte?)";
            break;

            case UHalf:
            o << "4-bit integer";
            break;

            default:
            break;

        }
        o << std::endl;

        o << "dimensions   : " << Nsize(I()) << " x " << Zsize(I()) << " x " << Ysize(I()) << " x " << Xsize(I());
        o << "	(noObjects x slices x rows x columns)" << std::endl;
        return o;
    }

    /** Sum this object with other file and keep in this object
      */
    void sumWithFile(const FileName &fn) {
        Image<T> aux;
        aux.read(fn);
        (*this)() += aux();
    }

    int readTiffInMemory(
        void* buf, size_t size, bool readdata = true, long int select_img = -1,
        bool mapData = false, bool is_2D = false
    );

    private:

    int _read(
        const FileName &name, fImageHandler &hFile, bool readdata = true, long int select_img = -1,
        bool mapData = false, bool is_2D = false
    );

    void _write(
        const FileName &name, fImageHandler &hFile, long int select_img=-1, 
        bool isStack = false, int mode = WRITE_OVERWRITE
    );

};

// Some image-specific operations

// For image normalisation
void normalise(
    Image<RFLOAT> &I,
    int bg_radius,
    RFLOAT white_dust_stddev, RFLOAT black_dust_stddev,
    bool do_ramp, bool is_helical_segment = false,
    RFLOAT helical_mask_tube_outer_radius_pix = -1.0,
    RFLOAT tilt_deg = 0.0, RFLOAT psi_deg = 0.0
);

Stats<RFLOAT> calculateBackgroundAvgStddev(
    Image<RFLOAT> &I,
    int bg_radius,
    bool is_helical_segment = false,
    RFLOAT helical_mask_tube_outer_radius_pix = -1.0,
    RFLOAT tilt_deg = 0.0, RFLOAT psi_deg = 0.0
);

void subtractBackgroundRamp(
    Image<RFLOAT> &I,
    int bg_radius,
    bool is_helical_segment = false,
    RFLOAT helical_mask_tube_outer_radius_pix = -1.0,
    RFLOAT tilt_deg = 0.0, RFLOAT psi_deg = 0.0
);

// For dust removal
void removeDust(
    Image<RFLOAT> &I, bool is_white, RFLOAT thresh, 
    RFLOAT avg, RFLOAT stddev
);

// for contrast inversion
void invert_contrast(Image<RFLOAT> &I);

// for image re-scaling
void rescale(Image<RFLOAT> &I, int mysize);

// for image re-windowing
void rewindow(Image<RFLOAT> &I, int mysize);

/// @defgroup ImageFormats Image Formats
/// @ingroup Images
// Functions belonging to this topic are commented in rw*.h
//@}

MultidimArray<RFLOAT>::MinMax getImageContrast(MultidimArray<RFLOAT> &image, RFLOAT minval, RFLOAT maxval, RFLOAT &sigma_contrast);

namespace ColourScheme {

    // Different paths through colourspace
    enum {
        greyscale, black_grey_red, blue_grey_white, 
        blue_grey_red, rainbow, cyan_black_yellow
    };

}

inline void greyToRGB(const int colour_scheme, const unsigned char grey, unsigned char &red, unsigned char &green, unsigned char &blue) {
    switch (colour_scheme) {

        case ColourScheme::greyscale:
        red = green = blue = grey;
        return;

        case ColourScheme::black_grey_red:
        if (grey >= 128) {
            red = 255; 
            blue = green = floor((RFLOAT) (255 - grey) * 2);
        } else { 
            red = green = blue = floor((RFLOAT) (grey * 2.0)); 
        }
        return;

        case ColourScheme::blue_grey_white:
        if (grey >= 128) {
            red = green = blue = floor((RFLOAT)((grey - 128) * 2));
        } else { 
            red = 0; 
            blue = green = floor((RFLOAT) (255 - 2 * grey));
        }
        return;

        case ColourScheme::blue_grey_red: {
            const RFLOAT a = grey / 85.0; // group
            const int X = floor(a);	//this is the integer part
            const unsigned char Y = floor(255 * (a - X)); //fractional part from 0 to 255
            switch (X) {
                case 0: red =   0; green = 255 - Y; blue = 255 - Y; break;
                case 1: red =   Y; green =       Y; blue =       Y; break;
                case 2: red = 255; green = 255 - Y; blue = 255 - Y; break;
                case 3: red = 255; green =       0; blue =       0; break;
            }
            return;
        }

        case ColourScheme::rainbow: {
            const RFLOAT a = (255 - grey) / 64.0; //invert and group
            const int X = floor(a);
            const unsigned char Y = floor(255 * (a - X)); //fractional part from 0 to 255
            switch(X) {
                case 0: red =     255; green =       Y; blue =   0; break;
                case 1: red = 255 - Y; green =     255; blue =   0; break;
                case 2: red =       0; green =     255; blue =   Y; break;
                case 3: red =       0; green = 255 - Y; blue = 255; break;
                case 4: red =       0; green =       0; blue = 255; break;
            }
            return;
        }
        case ColourScheme::cyan_black_yellow: {
            const RFLOAT d_rb = 3 * (grey - 128);
            const RFLOAT d_g = 3 * (std::abs(grey - 128) - 42);
            red   = (unsigned char) floor(std::min(255.0, std::max(0.0,  d_rb)));
            green = (unsigned char) floor(std::min(255.0, std::max(0.0,  d_g )));
            blue  = (unsigned char) floor(std::min(255.0, std::max(0.0, -d_rb)));
            return;
        }
    }
}
#endif
