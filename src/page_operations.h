#pragma once

#include <functional>

struct uhalf_t { unsigned char bits: 4; };

namespace transcription {
    // Given a page of page_size bytes,
    // reverse the bytes of each word,
    // each of size word_size.
    static void swapPage(char *page, size_t page_size, size_t word_size) {
        #ifdef DEBUG
            std::cerr << "DEBUG " << __func__ << ": "
            << "Swapping image data with swap= " << swap
            << " datatypesize= " << word_size
            << " pageNrElements " << n
            << std::endl;
        #endif
        if (word_size < 2) return;
        for (size_t i = 0; i < page_size; i += word_size)
            swapbytes(page + i, word_size);
    }

    template <typename T, typename U>
    static void page_cast_copy(T *dest, U *src, unsigned long int n) {
        if (typeid(T) == typeid(U)) {
            memcpy(dest, src, n * sizeof(T));
        } else {
            for (size_t i = 0; i < n; i++) {
                dest[i] = (T) src[i];
            }
        }
    }

    // Unfortunately, we cannot partially specialise template functions
    template <typename T, typename U=unsigned char>
    static void page_cast_copy_half(T *dest, U *src, unsigned long int n) throw (RelionError) {

        if (n % 2 != 0) {
            REPORT_ERROR((std::string) "Logic error in " + __func__ + "; for UHalf, pageSize must be even.");
        }

        for (size_t i = 0; 2 * i < n; i++) {
            // Here we are assuming the fill-order is LSB2MSB according to IMOD's
            // iiProcessReadLine() in libiimod/mrcsec.c.
            // The default fill-order in the TIFF specification is MSB2LSB
            // but IMOD assumes LSB2MSB even for TIFF.
            // See IMOD's iiTIFFCheck() in libiimod/iitif.c.
            dest[i * 2    ] = (T) ( (src)[i]       & 0b1111);
            dest[i * 2 + 1] = (T) (((src)[i] >> 4) & 0b1111);
        }
    }

    // Cast a page of n data from type U (index u) to type T
    template <typename T>
    static void castFromPage(T *dest, char *page, std::type_index u, unsigned long int n) {

        if (u == typeid(void)) {
            REPORT_ERROR("ERROR: datatype is Unknown_Type");
        }

        if (u == typeid(unsigned char)) {
            using U = unsigned char;
            page_cast_copy<T, U>(dest, (U*) page, n);
            return;
        }

        if (u == typeid(signed char)) {
            using U = signed char;
            page_cast_copy<T, U>(dest, (U*) page, n);
            return;
        }

        if (u == typeid(unsigned short)) {
            using U = unsigned short;
            page_cast_copy<T, U>(dest, (U*) page, n);
            return;
        }

        if (u == typeid(short)) {
            using U = short;
            page_cast_copy<T, U>(dest, (U*) page, n);
            return;
        }

        if (u == typeid(unsigned int)) {
            using U = unsigned int;
            page_cast_copy<T, U>(dest, (U*) page, n);
            return;
        }

        if (u == typeid(int)) {
            using U = int;
            page_cast_copy<T, U>(dest, (U*) page, n);
            return;
        }

        if (u == typeid(long)) {
            using U = long;
            page_cast_copy<T, U>(dest, (U*) page, n);
            return;
        }

        if (u == typeid(float)) {
            using U = float;
            page_cast_copy<T, U>(dest, (U*) page, n);
            return;
        }

        if (u == typeid(double)) {
            using U = RFLOAT;
            page_cast_copy<T, U>(dest, (U*) page, n);
            return;
        }

        if (u == typeid(uhalf_t)) {
            using U = unsigned char;
            page_cast_copy_half<T, U>(dest, (U*) page, n);
            return;
        }

        std::cerr << "Datatype= " << u.name() << std::endl;
        REPORT_ERROR(" ERROR: cannot cast datatype to T");

    }

    // Cast a page of n data from type T to type U (index u)
    template <typename T>
    static void castToPage(char *page, T *src, std::type_index u, unsigned long int n) {

        if (u == typeid(float)) {
            using U = float;
            page_cast_copy<U, T>((U*) page, src, n);
            return;
        }

        if (u == typeid(RFLOAT)) {
            using U = RFLOAT;
            page_cast_copy<U, T>((U*) page, src, n);
            return;
        }

        if (u == typeid(short)) {
            using U = short;
            page_cast_copy<U, T>((U*) page, src, n);
            return;
        }

        if (u == typeid(unsigned short)) {
            using U = unsigned short;
            page_cast_copy<U, T>((U*) page, src, n);
            return;
        }

        if (u == typeid(unsigned char)) {
            using U = unsigned char;
            page_cast_copy<U, T>((U*) page, src, n);
            return;
        }

        std::cerr << "outputDatatype= " << u.name() << std::endl;
        REPORT_ERROR(" ERROR: cannot cast T to outputDatatype");

    }

    template <typename T>
    static int copyViaPage(MultidimArray<T> &data, FILE *fimg, size_t bytes_per_slice, std::type_index index_u, size_t size_u, long off, long pad, bool swap) {
        static const size_t pagemax = 0x40000000;  // 1 GB (1 << 30)
        const size_t pagesize = std::min(bytes_per_slice, pagemax);
        const auto deleter = [pagesize] (char *ptr) { callocator<char>::deallocate(ptr, pagesize); };
        const auto page = std::unique_ptr<char, decltype(deleter)>(callocator<char>::allocate(pagesize), deleter);
        // Move the file pointer off bytes from the start of the file
        if (fseek(fimg, off, SEEK_SET) != 0) return -1;

        std::function<size_t (size_t)> how_many_pixels;
        if (index_u == typeid(uhalf_t))
            how_many_pixels = [      ] (size_t bytes) -> size_t { return bytes * 2; };
        else
            how_many_pixels = [size_u] (size_t bytes) -> size_t { return bytes / size_u; };

        unsigned long int pixel_progress = 0;  // Number of pixels processed so far
        for (size_t n = 0; n < Nsize(data); n++) {

            for (size_t j = 0; j < bytes_per_slice; j += pagemax) {
                // Read no more than than pagemax bytes in one go
                const size_t readsize = std::min(bytes_per_slice - j, pagemax);
                const unsigned long int num_pixels = how_many_pixels(readsize);

                #ifdef DEBUG
                std::cout << "NX = " << Xsize(data) << " NY = " << Ysize(data) << " NZ = " << Zsize(data) << std::endl;
                std::cout << "pagemax = " << pagemax << " bytes_per_slice = " << bytes_per_slice  << " readsize = " << readsize << " num_pixels = " << num_pixels << std::endl;
                #endif

                // Read into page
                if (fread(page.get(), readsize, 1, fimg) != 1) return -2;
                // Maybe swap bytes (won't work for UHalf)
                if (swap) swapPage(page.get(), readsize, size_u);
                // Copy into data, while casting to T
                castFromPage(data.data + pixel_progress, page.get(), index_u, num_pixels);
                pixel_progress += num_pixels;
            }
            // Advance the file pointer by pad bytes to get to the next slice
            if (fseek(fimg, pad, SEEK_CUR) != 0) return -1;

        }
        return 0;
    }

};
