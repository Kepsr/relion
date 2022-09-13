#pragma once
// Intel MKL provides an FFTW-like interface, so this is enough.
#include <fftw3.h>
#include <fcntl.h>
#include <unistd.h>
#include <sys/mman.h>
#include "src/filename.h"
#include "src/error.h"

// Allocator types for MultidimArray

struct relion_aligned_mallocator {

    static void* allocate(size_t n) {
        void *ptr = fftw_malloc(n);
        if (!ptr) REPORT_ERROR("Allocate: No space left");
        return ptr;
    }

    static void deallocate(void* ptr, size_t size) {
        fftw_free(ptr);
    }

};

struct mmapper {

    FileName mapFile;  // Mapped file name
    int mFd;           // Mapped file handler

    mmapper(): mapFile(""), mFd(0) {}

    void* allocate(size_t size, off_t offset = 0) {
        mapFile.initRandom(8);
        mapFile = mapFile.addExtension("tmp");
        mFd = open(mapFile.c_str(), O_RDWR | O_CREAT | O_TRUNC, S_IRUSR | S_IWUSR | S_IRGRP | S_IWGRP);
        if (mFd == -1)
            REPORT_ERROR((std::string) "MultidimArray<T>::" + __func__ + ": Error creating map file.");

        if (lseek(mFd, offset, SEEK_SET) == -1 || write(mFd, "", 1) == -1) {
            close(mFd);
            REPORT_ERROR((std::string) "MultidimArray<T>::" + __func__ + ": Error 'stretching' the map file.");
        }

        void *ptr = mmap(0, size, PROT_READ | PROT_WRITE, MAP_SHARED, mFd, 0);
        if (ptr == (void*) -1)
            REPORT_ERROR((std::string) "MultidimArray<T>::" + __func__ + ": mmap failed.");
        return ptr;

    }

    void deallocate(void* ptr, size_t size) {
        munmap(ptr, size);
        close(mFd);
        remove(mapFile.c_str());
    }

};
