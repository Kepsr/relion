#ifndef ACC_PTR_H_
#define ACC_PTR_H_

#include "src/acc/settings.h"
#ifdef CUDA
#include "src/acc/cuda/cuda_settings.h"
#include <cuda_runtime.h>
#include "src/acc/cuda/custom_allocator.cuh"
#include "src/acc/cuda/cuda_mem_utils.h"
#include "src/acc/cuda/shortcuts.cuh"
#else
#include "src/acc/cpu/cpu_settings.h"
#endif

#include <signal.h>
#include <fstream>
#include <iostream>
#include <vector>
#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <cstring>

#include "src/macros.h"
#include "src/error.h"
#include "src/parallel.h"

#ifndef MEM_ALIGN
#define MEM_ALIGN 64
#endif

#ifdef DEBUG_CUDA
const bool __debug_cuda__ = true;
#else
const bool __debug_cuda__ = false;
#endif

namespace acc {

    enum Type {unset, cuda, cpu};

    #ifdef CUDA
    const Type type = cuda;
    #else
    const Type type = cpu;
    #endif

};

#ifdef CUDA
using StreamType     = cudaStream_t;
using AllocatorType  = CudaCustomAllocator;
using AllocationType = CudaCustomAllocator::Alloc;
#else
// Dummy types
using StreamType     = float;
using AllocatorType  = double;
using AllocationType = double;
#endif

template <typename T, acc::Type accType = acc::type>
class AccPtr {

    protected:

    AllocatorType *allocator;
    AllocationType *alloc;
    StreamType stream;

    size_t size;  // Size used when copying data from and to device
    T *hPtr, *dPtr;  // Host and device pointers
    bool mustFreeHost, mustFreeDevice;  // Does the host/device need freeing?

    public:

    inline void setDoFreeHost(bool p) {
        mustFreeHost = p;
    }

    static void HandleDebugFatal(const char *err, const char *file, int line) {
        fprintf(stderr, "DEBUG ERROR: %s in %s:%d\n", err, file, line);
        fflush(stdout);
        #ifdef DEBUG_CUDA
        raise(SIGSEGV);
        #else
        CRITICAL(ERRGPUKERN);
        #endif
    }

    static void HandleDebugInformational(const char *err, const char *file, int line) {
            fprintf(stderr, "POSSIBLE ISSUE: %s in %s:%d\n", err, file, line);
            fflush(stdout);
    }

    /*======================================================
                          CONSTRUCTORS
    ======================================================*/

    AccPtr(
        AllocatorType *allocator = nullptr, StreamType stream = cudaStreamPerThread,
        size_t size = 0, T *hptr = nullptr, T *dptr = nullptr
    ):
        allocator(allocator), alloc(nullptr), stream(stream),
        size(size), hPtr(hptr), dPtr(dptr), mustFreeHost(false), mustFreeDevice(false)
    {
        if (size && hPtr && posix_memalign((void**) &hPtr, MEM_ALIGN, sizeof(T) * size))
            CRITICAL(RAMERR);
    }

    AccPtr(
        size_t size, AllocatorType *allocator = nullptr,
        StreamType stream = cudaStreamPerThread,
        T *hptr = nullptr, T *dptr = nullptr
    ):
        allocator(allocator), alloc(nullptr), stream(stream),
        size(size), hPtr(hptr), dPtr(dptr), mustFreeHost(false), mustFreeDevice(false)
    {
        if (size && hPtr && posix_memalign((void**) &hPtr, MEM_ALIGN, sizeof(T) * size))
            CRITICAL(RAMERR);
    }

    AccPtr(size_t size, StreamType stream): AccPtr(nullptr, stream, size) {
        hostAlloc();
    }

    // Copy ctor
    AccPtr(const AccPtr &ptr):
        size(ptr.size), hPtr(ptr.hPtr), dPtr(ptr.dPtr), mustFreeHost(false), mustFreeDevice(false),
        allocator(ptr.allocator), alloc(nullptr), stream(ptr.stream)
    {}

    AccPtr(const AccPtr<T> &ptr, size_t start_idx, size_t size):
        size(size), hPtr(ptr.hPtr + start_idx), dPtr(ptr.dPtr + start_idx), mustFreeHost(false), mustFreeDevice(false),
        allocator(ptr.allocator), alloc(nullptr), stream(ptr.stream)
    {}

    /*======================================================
                         METHOD BODY
    ======================================================*/

    void markReadyEvent() {
        #ifdef CUDA
        if (accType == acc::cuda) {
            if (__debug_cuda__ && !alloc)
                HandleDebugFatal("markReadyEvent called on null allocation.\n", __FILE__, __LINE__);
            alloc->markReadyEvent(stream);
        }
        #endif
    }

    /**
     * Allocate memory on device
     */
    void deviceAlloc() {
        #ifdef CUDA
        if (accType == acc::cuda) {
            if (__debug_cuda__) {
                if (size == 0)
                    HandleDebugFatal("deviceAlloc called with size == 0", __FILE__, __LINE__);
                if (mustFreeDevice)
                    HandleDebugFatal("Device double allocation.\n", __FILE__, __LINE__);
            }
            alloc = allocator->alloc(size * sizeof(T));
            dPtr = (T*) alloc->getPtr();
            mustFreeDevice = true;
        }
        #endif
    }

    /**
     * Allocate memory on host
     */
    void hostAlloc() {
        if (__debug_cuda__) {
            if (size == 0)
                HandleDebugFatal("hostAlloc called with size == 0", __FILE__, __LINE__);
            if (mustFreeHost)
                HandleDebugFatal("Host double allocation.\n", __FILE__, __LINE__);
        }
        // TODO - alternatively, this could be aligned std::vector
        if (posix_memalign((void**) &hPtr, MEM_ALIGN, sizeof(T) * size))
            CRITICAL(RAMERR);
        memset(hPtr, 0, sizeof(T) * size);
        mustFreeHost = true;
    }

    void allAlloc() {
        deviceAlloc();
        hostAlloc();
    }

    void accAlloc() {
        accType == acc::cuda ? deviceAlloc() : hostAlloc();
    }

    // Allocate storage of a new size for the array
    // Resize while retaining as much of the original contents as possible
    void resizeHost(size_t newSize) {

        /// TODO: alternatively, this could be aligned std::vector
        T* new_hptr;
        if (posix_memalign((void**) &new_hptr, MEM_ALIGN, sizeof(T) * newSize))
            CRITICAL(RAMERR);

        if (size > 0 && hPtr) {
            // Copy what we can
            memcpy(new_hptr, hPtr, sizeof(T) * std::min(size, newSize));
            // Initialize any remaining memory
            if (newSize > size)
                memset(new_hptr, 0, sizeof(T) * (newSize - size));
        } else {
            // Nothing to copy
            memset(new_hptr, 0, sizeof(T) * newSize);
        }

        if (__debug_cuda__) {
            if (dPtr)
                HandleDebugFatal("resizeHost: Resizing host with present device allocation.\n", __FILE__, __LINE__);
            if (newSize == 0)
                HandleDebugInformational("resizeHost: Array resized to size zero (permitted with fear).  Something may break downstream\n", __FILE__, __LINE__);
        }
        freeHost();
        setSize(newSize);
        setHostPtr(new_hptr);
        mustFreeHost = true;
    }

    /**
     * Initialise device memory with provided value
     */
    void deviceInit(int value = 0) {
        #ifdef CUDA
        if (accType == acc::cuda) {
            if (__debug_cuda__) {
                if (!dPtr)
                    HandleDebugFatal("cudaMemInit requested before allocation in deviceInit().\n", __FILE__, __LINE__);
            }
            cudaMemInit<T>(dPtr, value, size, stream);
        }
        #endif
    }

    /**
     * Initialise host memory with provided value
     */
    void hostInit(int value = 0) {
        if (__debug_cuda__) {
            if (!hPtr)
                HandleDebugFatal("Memset requested before allocation in hostInit().\n", __FILE__, __LINE__);
        }
        memset(hPtr, value, size * sizeof(T));
    }

    /**
     * Initialise memory with provided value
     */
    void accInit(int value = 0) {
        accType == acc::cuda ? deviceInit(value) : hostInit(value);
    }

    /**
     * Initialise all used memory with provided value
     */
    void allInit(int value = 0) {
        hostInit(value);
        deviceInit(value);
    }

    /**
     * Copy bytes from host to device
     */
    void cpToDevice() {
        #ifdef CUDA
        if (accType == acc::cuda) {
            if (__debug_cuda__) {
                if (!dPtr)
                    HandleDebugFatal("cpToDevice() called before allocation.\n", __FILE__, __LINE__);
                if (!hPtr)
                    HandleDebugFatal("null host pointer in cpToDevice().\n", __FILE__, __LINE__);
            }
            CudaShortcuts::cpyHostToDevice<T>(hPtr, dPtr, size, stream);
        }
        #endif
    }

    /**
     * alloc and copy
     */
    void putOnDevice() {
        deviceAlloc();
        cpToDevice();
    }

    /**
     * Copy bytes from device to host
     */
    void cpToHost() {
        #ifdef CUDA
        if (accType == acc::cuda) {
            if (__debug_cuda__) {
                if (!dPtr)
                    HandleDebugFatal("cp_to_host() called before device allocation.\n", __FILE__, __LINE__);
                if (!hPtr)
                    HandleDebugFatal("null host pointer in cp_to_host().\n", __FILE__, __LINE__);
            }
            cudaCpyDeviceToHost<T>(dPtr, hPtr, size, stream);
        }
        #endif
    }

    /**
     * Copy a number (size) of bytes from device pointer to the provided new device pointer
     */
    void cpOnDevice(T* const dest) {
        #ifdef CUDA
        if (accType == acc::cuda) {
            if (__debug_cuda__) {
                if (!dest)
                    HandleDebugFatal("nullptr given in cpOnDevice(dest).\n", __FILE__, __LINE__);
            }
            CudaShortcuts::cpyDeviceToDevice(dPtr, dest, size, stream);
        }
        #endif
    }

    /**
     * Copy bytes into another host pointer
     */
    void cpOnHost(T* const dest) {
        if (__debug_cuda__) {
            if (!dest)
                HandleDebugFatal("nullptr given in cp_on_host(dest).\n", __FILE__, __LINE__);
            if (!hPtr)
                HandleDebugFatal("Null input pointer given in cp_on_host(hPtr).\n", __FILE__, __LINE__);
        }
        memcpy(dest, hPtr, size * sizeof(T));
    }

    void cpOnAcc(T* const dest) {
        accType == acc::cuda ? cpOnDevice(dest) : cpOnHost(dest);
    }

    void streamSync() {
        #ifdef CUDA
        if (accType == acc::cuda)
            DEBUG_HANDLE_ERROR(cudaStreamSynchronize(stream));
        #endif
    }

    T getAccValueAt(size_t idx) {
        return accType == acc::cuda ? getDeviceAt(idx) : hPtr[idx];
    }

    T getDeviceAt(size_t idx) {
        #ifdef CUDA
        if (accType == acc::cuda) {
            T value;
            cudaCpyDeviceToHost<T>(&dPtr[idx], &value, 1, stream);
            streamSync();
            return value;
        } else
        #endif
        return {};
    }

    void dumpDeviceToFile(const std::string &fileName) {
        std::ofstream f (fileName.c_str());
        #ifdef CUDA
        if (accType == acc::cuda) {
            T tmp[size];
            cudaCpyDeviceToHost<T>(dPtr, tmp, size, stream);
            streamSync();
            for (unsigned i = 0; i < size; i++)
                f << tmp[i] << std::endl;
        } else
        #endif
        {
            f << "Pointer has no device support." << std::endl;
        }
    }

    void dumpHostToFile(const std::string &fileName) {
        std::ofstream f (fileName.c_str());
        for (unsigned i = 0; i < size; i ++)
            f << hPtr[i] << std::endl;
    }

    void dumpAccToFile(const std::string &fileName) {
        accType == acc::cuda ?  dumpDeviceToFile(fileName) : dumpHostToFile(fileName);
    }

    /**
     * Delete device data
     */
    void freeDevice() {
        if (!mustFreeDevice) return;
        #ifdef CUDA
        if (accType == acc::cuda) {
            #ifdef DEBUG_CUDA
            if (!dPtr)
                HandleDebugFatal("Free device memory was called on nullptr in free_device().\n", __FILE__, __LINE__);
            #endif

            if (alloc->getReadyEvent() == 0)
                alloc->markReadyEvent(stream);
            alloc->doFreeWhenReady();
            alloc = nullptr;

			// DEBUG_HANDLE_ERROR(cudaFree(dPtr));

            dPtr = nullptr;
            mustFreeDevice = false;
        }
        #endif
    }

    /**
     * Delete host data
     */
    void freeHost() {
        if (!mustFreeHost) return;
        if (__debug_cuda__ && !hPtr)
            HandleDebugFatal("free_host() called on nullptr.\n", __FILE__, __LINE__);
        mustFreeHost = false;
        ::free(hPtr);
        hPtr = nullptr;
    }

    void free() {
        freeHost();
        freeDevice();
    }

    ~AccPtr() {
        free();
    }

    /*======================================================
                       GETTERS AND SETTERS
    ======================================================*/

    inline bool willFreeHost() const {
        return mustFreeHost;
    }

    inline bool willFreeDevice() const {
        return mustFreeDevice;
    }

    inline void setStream(StreamType s) {
        stream = s;
    }

    inline StreamType getStream() const {
        return stream;
    }

    inline void setSize(size_t s) {
        size = s;
    }

    inline size_t getSize() const {
        return size;
    }

    inline T *getDevicePtr() const {
        return dPtr;
    }

    inline T *getHostPtr() const {
        return hPtr;
    }

    inline T *getAccPtr() const {
        return accType == acc::cuda ? dPtr : hPtr;
    }

    void setAllocator(AllocatorType *a) {
        freeDevice();
        allocator = a;
    };

    inline AllocatorType *getAllocator() const {
        return allocator;
    }

    void setDevicePtr(T *ptr) {
        if (__debug_cuda__)
            if (mustFreeDevice)
                HandleDebugFatal("Device pointer set without freeing the old one.\n", __FILE__, __LINE__);
        dPtr = ptr;
    }

    void setHostPtr(T *ptr) {
        if (__debug_cuda__)
            if (mustFreeHost)
                HandleDebugFatal("Host pointer set without freeing the old one.\n", __FILE__, __LINE__);
        hPtr = ptr;
    }

    void setAccPtr(T *ptr) {
        accType == acc::cuda ? setDevicePtr(ptr) : setHostPtr(ptr);
    }

    template <typename U>
    AccPtr<U, accType> make(size_t s = 0) {
        return AccPtr<U, accType>(allocator, stream, s);
    }

};


class AccPtrBundle: public AccPtr<unsigned char> {

    using byte_t = unsigned char;

    size_t current_packed_pos;

    public:

    template <acc::Type accType = acc::type>
    AccPtrBundle(StreamType stream, AllocatorType *allocator, size_t size = 0):
        AccPtr<byte_t, accType>(allocator, stream, size),
        current_packed_pos(0) {}

    template <typename T>
    void pack(AccPtr<T> &ptr) {
        #ifdef CUDA
        #ifdef DEBUG_CUDA
        if (current_packed_pos + ptr.getSize() > size)
            HandleDebugFatal("Packing exceeds bundle total size.\n", __FILE__, __LINE__);
        if (!hPtr)
            HandleDebugFatal("Pack called on null host pointer.\n", __FILE__, __LINE__);
        #endif
        if (ptr.getHostPtr())
            memcpy(&hPtr[current_packed_pos], ptr.getHostPtr(), ptr.getSize() * sizeof(T));
        ptr.freeHost();
        ptr.setHostPtr  ((T*) &hPtr[current_packed_pos]);
        ptr.setDevicePtr((T*) &dPtr[current_packed_pos]);

        current_packed_pos += ptr.getSize() * sizeof(T);
        #else
        if (!ptr.getHostPtr()) ptr.hostAlloc();
        #endif
    }

    // Override allocation methods and block for no device

    void allAlloc() {
        #ifdef CUDA
        AccPtr<byte_t>::allAlloc();
        #endif
    }

    void hostAlloc() {
        #ifdef CUDA
        AccPtr<byte_t>::hostAlloc();
        #endif
    }

};

template <acc::Type accType = acc::cuda>
class AccPtrFactory {

    StreamType stream;

    public:

    AllocatorType *allocator;

    AccPtrFactory(): allocator(nullptr), stream(0) {}

    AccPtrFactory(AllocatorType *alloc, StreamType stream = 0):
        allocator(alloc), stream(stream)
    {}

    template <typename T>
    AccPtr<T> make(size_t size = 0, StreamType stream = cudaStreamPerThread) {
        return AccPtr<T, accType>(size, allocator, stream);
    }

    AccPtrBundle makeBundle(size_t size = 0) {
        return AccPtrBundle(stream, allocator, size);
    }

};

#endif
