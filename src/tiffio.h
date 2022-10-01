#include <tiffio.h>
#include "src/error.h"

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
