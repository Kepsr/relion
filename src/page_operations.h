#pragma once

struct uhalf_t { unsigned char bits: 4; };

namespace pages {
        
    // Swap a page of n elements, each of size size
    static void swapPage(char *page, size_t n, size_t size, size_t swap) {
        #ifdef DEBUG
            std::cerr << "DEBUG " << __func__ << ": Swapping image data with swap= "
            << swap << " datatypesize= " << size
            << " pageNrElements " << n
            << std::endl;
        #endif

        const size_t di = swap == 1 ? size : swap;
        for (size_t i = 0; i < n; i += di) swapbytes(page + i, di);
    }

    template <typename T, typename U>
    static void page_cast_copy(T *dest, U *src, size_t size) {
        if (typeid(T) == typeid(U)) {
            memcpy(dest, src, size * sizeof(T));
        } else {
            for (size_t i = 0; i < size; i++) {
                dest[i] = (T) src[i];
            }
        }
    }

    // Unfortunately, we cannot partially specialise template functions
    template <typename T, typename U=unsigned char>
    static void page_cast_copy_half(T *dest, U *src, size_t size) throw (RelionError) {

        if (size % 2 != 0) {
            REPORT_ERROR((std::string) "Logic error in " + __func__ + "; for UHalf, pageSize must be even.");
        }

        for (size_t i = 0; 2 * i < size; i++) {
            // Here we are assuming the fill-order is LSB2MSB according to IMOD's
            // iiProcessReadLine() in libiimod/mrcsec.c.
            // The default fill-order in the TIFF specification is MSB2LSB
            // but IMOD assumes LSB2MSB even for TIFF.
            // See IMOD's iiTIFFCheck() in libiimod/iitif.c.
            dest[i * 2    ] = (T) ( (src)[i]       & 0b1111);
            dest[i * 2 + 1] = (T) (((src)[i] >> 4) & 0b1111);
        }
    }

    // Cast a page of data from type U (index u) to type T
    template <typename T>
    void castFromPage(T *dest, char *page, std::type_index u, size_t pageSize) {

        if (u == typeid(void)) {
            REPORT_ERROR("ERROR: datatype is Unknown_Type");
        }

        if (u == typeid(unsigned char)) {
            using U = unsigned char;
            page_cast_copy<T, U>(dest, (U*) page, pageSize);
            return;
        }

        if (u == typeid(signed char)) {
            using U = signed char;
            page_cast_copy<T, U>(dest, (U*) page, pageSize);
            return;
        }

        if (u == typeid(unsigned short)) {
            using U = unsigned short;
            page_cast_copy<T, U>(dest, (U*) page, pageSize);
            return;
        }

        if (u == typeid(short)) {
            using U = short;
            page_cast_copy<T, U>(dest, (U*) page, pageSize);
            return;
        }

        if (u == typeid(unsigned int)) {
            using U = unsigned int;
            page_cast_copy<T, U>(dest, (U*) page, pageSize);
            return;
        }

        if (u == typeid(int)) {
            using U = int;
            page_cast_copy<T, U>(dest, (U*) page, pageSize);
            return;
        }

        if (u == typeid(long)) {
            using U = long;
            page_cast_copy<T, U>(dest, (U*) page, pageSize);
            return;
        }

        if (u == typeid(float)) {
            using U = float;
            page_cast_copy<T, U>(dest, (U*) page, pageSize);
            return;
        }

        if (u == typeid(double)) {
            using U = RFLOAT;
            page_cast_copy<T, U>(dest, (U*) page, pageSize);
            return;
        }

        if (u == typeid(uhalf_t)) {
            using U = unsigned char;
            page_cast_copy_half<T, U>(dest, (U*) page, pageSize);
            return;
        }

        std::cerr << "Datatype= " << u.name() << std::endl;
        REPORT_ERROR(" ERROR: cannot cast datatype to T");

    }

    // Cast a page of data from type T to type U (index u)
    template <typename T>
    void castToPage(char *page, T *src, std::type_index u, size_t pageSize) {

        if (u == typeid(float)) {
            using U = float;
            page_cast_copy<U, T>((U*) page, src, pageSize);
            return;
        }

        if (u == typeid(RFLOAT)) {
            using U = RFLOAT;
            page_cast_copy<U, T>((U*) page, src, pageSize);
            return;
        }

        if (u == typeid(short)) {
            using U = short;
            page_cast_copy<U, T>((U*) page, src, pageSize);
            return;
        }

        if (u == typeid(unsigned short)) {
            using U = unsigned short;
            page_cast_copy<U, T>((U*) page, src, pageSize);
            return;
        }

        if (u == typeid(unsigned char)) {
            using U = unsigned char;
            page_cast_copy<U, T>((U*) page, src, pageSize);
            return;
        }

        std::cerr << "outputDatatype= " << u.name() << std::endl;
        REPORT_ERROR(" ERROR: cannot cast T to outputDatatype");

    }

    template <typename T>
    int allocatePage(FILE *fimg, size_t pagesize, size_t off, std::type_index index_u, size_t size_u, MultidimArray<T> &data, size_t swap, size_t pad) {
        // pagesize: size of object
        static const size_t pagemax = 0x40000000;  // 1 GB (1 << 30)
        const size_t memsize = std::max(pagesize, pagemax) * sizeof(char);
        const auto deleter = [memsize] (char *ptr) { callocator<char>::deallocate(ptr, memsize); };
        const auto page = std::unique_ptr<char, decltype(deleter)>(callocator<char>::allocate(memsize), deleter);
        // Because we requested XYsize to be even for UHalf, this is always safe.
        if (fseek(fimg, off, SEEK_SET) != 0) return -1;

        std::function<size_t (size_t)> pixels_per_page;
        if (index_u == typeid(uhalf_t))
            pixels_per_page = [      ] (size_t readsize) -> size_t { return readsize * 2; };
        else
            pixels_per_page = [size_u] (size_t readsize) -> size_t { return readsize / size_u; };

        size_t pixel_progress = 0;  // Number of pixels processed so far
        for (size_t n = 0; n < Nsize(data); n++) {
            for (size_t j = 0; j < pagesize; j += pagemax) {
                // Read next page
                // Divide pages larger than pagemax
                size_t readsize = std::min(pagesize - j, pagemax);
                size_t num_pixels = pixels_per_page(readsize);

                #ifdef DEBUG
                std::cout << "NX = " << Xsize(data) << " NY = " << Ysize(data) << " NZ = " << Zsize(data) << std::endl;
                std::cout << "pagemax = " << pagemax << " pagesize = " << pagesize  << " readsize = " << readsize << " num_pixels = " << num_pixels << std::endl;
                #endif

                // Read page from disk
                if (fread(page.get(), readsize, 1, fimg) != 1) return -2;
                // Swap bytes if required
                if (swap) swapPage(page.get(), readsize, size_u, swap);
                // Cast to T
                castFromPage(data.data + pixel_progress, page.get(), index_u, num_pixels);
                pixel_progress += num_pixels;
            }
            if (pad > 0) {
                // fread(padpage, pad, 1, fimg);
                if (fseek(fimg, pad, SEEK_CUR) != 0) return -1;
            }
        }
        return 0;
    }

};
