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

#define ACC_PTR_DEBUG_FATAL( err ) (HandleAccPtrDebugFatal( err, __FILE__, __LINE__ ))
static void HandleAccPtrDebugFatal( const char *err, const char *file, int line )
{
        fprintf(stderr, "DEBUG ERROR: %s in %s:%d\n", err, file, line );
        fflush(stdout);
#ifdef DEBUG_CUDA
        raise(SIGSEGV);
#else
        CRITICAL(ERRGPUKERN);
#endif
}

#define ACC_PTR_DEBUG_INFO( err ) (HandleAccPtrDebugInformational( err, __FILE__, __LINE__ ))
static void HandleAccPtrDebugInformational( const char *err, const char *file, int line )
{
        fprintf(stderr, "POSSIBLE ISSUE: %s in %s:%d\n", err, file, line );
        fflush(stdout);
}

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

template <typename T>
class AccPtr
{

    protected:

    AllocatorType *allocator;
    AllocationType *alloc;
    StreamType stream;

    acc::Type accType;

    size_t size;  // Size used when copying data from and to device
    T *hPtr, *dPtr;  // Host and device pointers
    bool doFreeHost, doFreeDevice;  // Does the host/device need to be freed?

    public:

    inline void setDoFreeHost(bool p) {
        doFreeHost = p;
    }

    /*======================================================
                          CONSTRUCTORS
    ======================================================*/

    AccPtr(
        AllocatorType *allocator = nullptr, StreamType stream = cudaStreamPerThread,
        acc::Type accType = acc::type, size_t size = 0,
        T *h_start = nullptr, T *d_start = nullptr
    ):
        size(size), hPtr(h_start), dPtr(d_start), doFreeHost(false),
        doFreeDevice(false), allocator(allocator), alloc(nullptr), stream(stream),
        accType(accType)
    {
        if (size && hPtr && posix_memalign((void**) &hPtr, MEM_ALIGN, sizeof(T) * size))
            CRITICAL(RAMERR);
    }

    AccPtr(
        size_t size, AllocatorType *allocator = nullptr,
        StreamType stream = cudaStreamPerThread,
        acc::Type accType = acc::type,
        T *h_start = nullptr, T *d_start = nullptr
    ):
        size(size), hPtr(h_start), dPtr(d_start), doFreeHost(false),
        doFreeDevice(false), allocator(allocator), alloc(nullptr), stream(stream),
        accType(accType)
    {
        if (size && hPtr && posix_memalign((void**) &hPtr, MEM_ALIGN, sizeof(T) * size))
            CRITICAL(RAMERR);
    }

    AccPtr(size_t size, StreamType stream):
        size(size), dPtr(nullptr), doFreeHost(true),
        doFreeDevice(false), allocator(nullptr), alloc(nullptr), stream(stream),
        accType(acc::type)
    {
        if (posix_memalign((void**) &hPtr, MEM_ALIGN, sizeof(T) * size))
            CRITICAL(RAMERR);
    }

    // Copy ctor
    AccPtr(const AccPtr &ptr):
        size(ptr.size), hPtr(ptr.hPtr), dPtr(ptr.dPtr), doFreeHost(false),
        doFreeDevice(false), allocator(ptr.allocator), alloc(nullptr), stream(ptr.stream),
        accType(ptr.accType)
    {}

    AccPtr(const AccPtr<T> &ptr, size_t start_idx, size_t size):
        size(size), hPtr(&ptr.hPtr[start_idx]), dPtr(&ptr.dPtr[start_idx]), doFreeHost(false),
        doFreeDevice(false), allocator(ptr.allocator), alloc(nullptr), stream(ptr.stream),
        accType(ptr.accType)
    {}

    /*======================================================
                         METHOD BODY
    ======================================================*/

    void markReadyEvent()
    {
        #ifdef CUDA
        if (accType == acc::cuda)
        {
            #ifdef DEBUG_CUDA
            if (!alloc)
                ACC_PTR_DEBUG_FATAL("markReadyEvent called on null allocation.\n");
            #endif
            alloc->markReadyEvent(stream);
        }
        #endif
    }

    /**
     * Allocate memory on device
     */
    void deviceAlloc()
    {
        #ifdef CUDA
        if (accType == acc::cuda)
        {
            #ifdef DEBUG_CUDA
            if (size == 0)
                ACC_PTR_DEBUG_FATAL("deviceAlloc called with size == 0");
            if (doFreeDevice)
                ACC_PTR_DEBUG_FATAL("Device double allocation.\n");
            #endif
            doFreeDevice = true;

            alloc = allocator->alloc(size * sizeof(T));
            dPtr = (T*) alloc->getPtr();
        }
        #endif
    }

    /**
     * Allocate memory on host
     */
    void hostAlloc()
    {
        #ifdef DEBUG_CUDA
        if (size == 0)
            ACC_PTR_DEBUG_FATAL("deviceAlloc called with size == 0");
        if (doFreeHost)
            ACC_PTR_DEBUG_FATAL("Host double allocation.\n");
        #endif
        doFreeHost = true;
        // TODO - alternatively, this could be aligned std::vector
        if (posix_memalign((void**) &hPtr, MEM_ALIGN, sizeof(T) * size))
            CRITICAL(RAMERR);
    }

    void allAlloc()
    {
        deviceAlloc();
        hostAlloc();
    }

    void accAlloc()
    {
        if (accType == acc::cuda)
            deviceAlloc();
        else
            hostAlloc();
    }

    // Allocate storage of a new size for the array
    void resizeHost(size_t newSize)
    {
        #ifdef DEBUG_CUDA
        if (size == 0)
            ACC_PTR_DEBUG_INFO("Resizing from size zero (permitted).\n");
        #endif
        // TODO - alternatively, this could be aligned std::vector
        T* newArr;
        if (posix_memalign((void**) &newArr, MEM_ALIGN, sizeof(T) * newSize))
            CRITICAL(RAMERR);
        memset(newArr, 0x0, sizeof(T) * newSize);

        #ifdef DEBUG_CUDA
        if (dPtr)
            ACC_PTR_DEBUG_FATAL("resizeHost: Resizing host with present device allocation.\n");
        if (newSize == 0)
            ACC_PTR_DEBUG_INFO("resizeHost: Array resized to size zero (permitted with fear).  Something may break downstream\n");
        #endif
        freeHost();
        setSize(newSize);
        setHostPtr(newArr);
        doFreeHost = true;
    }

    // Resize retaining as much of the original contents as possible
    void resizeHostCopy(size_t newSize)
    {
        #ifdef DEBUG_CUDA
        // if (size == 0)
        // 	ACC_PTR_DEBUG_INFO("Resizing from size zero (permitted).\n");
        #endif
        // TODO - alternatively, this could be aligned std::vector
        T* newArr;
        if (posix_memalign((void**) &newArr, MEM_ALIGN, sizeof(T) * newSize))
            CRITICAL(RAMERR);

        // Copy in what we can from the original matrix
        if (size > 0 && hPtr)
        {
            memcpy(newArr, hPtr, sizeof(T) * std::min(size, newSize));

            // Initialize remaining memory if any
            if (newSize > size)
            {
                memset(newArr, 0x0, sizeof(T) * (newSize - size));
            }
        }

        // There was nothing from before to copy - clear new memory
        if (!hPtr)
        {
            memset(newArr, 0x0, sizeof(T) * newSize);
        }

        #ifdef DEBUG_CUDA
        if (dPtr)
            ACC_PTR_DEBUG_FATAL("resizeHostCopy: Resizing host with present device allocation.\n");
        if (newSize == 0)
            ACC_PTR_DEBUG_INFO("resizeHostCopy: Array resized to size zero (permitted with fear).  Something may break downstream\n");
        #endif
        freeHost();
        setSize(newSize);
        setHostPtr(newArr);
        doFreeHost = true;
    }

    /**
     * Initiate device memory with provided value
     */
    void deviceInit(int value)
    {
        #ifdef CUDA
        if (accType == acc::cuda)
        {
            #ifdef DEBUG_CUDA
            if (!dPtr)
                ACC_PTR_DEBUG_FATAL("Memset requested before allocation in deviceInit().\n");
            #endif
            cudaMemInit<T>(dPtr, value, size, stream);
        }
        #endif
    }

    /**
     * Initiate host memory with provided value
     */
    void hostInit(int value)
    {
        #ifdef DEBUG_CUDA
        if (!hPtr)
            ACC_PTR_DEBUG_FATAL("Memset requested before allocation in hostInit().\n");
        #endif
        memset(hPtr, value, size * sizeof(T));
    }

    /**
     * Initiate memory with provided value
     */
    void accInit(int value)
    {
        if (accType == acc::cuda)
            deviceInit(value);
        else
            hostInit(value);
    }

    /**
     * Initiate all used memory with provided value
     */
    void allInit(int value)
    {
        hostInit(value);
        if (accType == acc::cuda)
            deviceInit(value);
    }

    /**
     * Copy a number (size) of bytes to device stored in the host pointer
     */
    void cpToDevice()
    {
        #ifdef CUDA
        if (accType == acc::cuda)
        {
            #ifdef DEBUG_CUDA
            if (!dPtr)
                ACC_PTR_DEBUG_FATAL("cpToDevice() called before allocation.\n");
            if (!hPtr)
                ACC_PTR_DEBUG_FATAL("null host pointer in cpToDevice().\n");
            #endif
            CudaShortcuts::cpyHostToDevice<T>(hPtr, dPtr, size, stream);
        }
        #endif
    }

    /**
     * Copy a number (size) of bytes to device stored in the provided host pointer
     */
    void cpToDevice(T * hostPtr)
    {
        if (accType == acc::cuda)
        {
            #ifdef DEBUG_CUDA
            if (!hostPtr)
                ACC_PTR_DEBUG_FATAL("nullptr given in cpToDevice(hostPtr).\n");
            #endif
            hPtr = hostPtr;
            cpToDevice();
        }
    }

    /**
     * alloc and copy
     */
    void putOnDevice()
    {
        deviceAlloc();
        cpToDevice();
    }

    /**
     * Copy a number (size) of bytes from device to the host pointer
     */
    void cpToHost()
    {
        #ifdef CUDA
        if (accType == acc::cuda)
        {
            #ifdef DEBUG_CUDA
            if (!dPtr)
                ACC_PTR_DEBUG_FATAL("cp_to_host() called before device allocation.\n");
            if (!hPtr)
                ACC_PTR_DEBUG_FATAL("null host pointer in cp_to_host().\n");
            #endif
            cudaCpyDeviceToHost<T>(dPtr, hPtr, size, stream);
        }
        #endif
    }

    /**
     * Copy a number (thisSize) of bytes from device to the host pointer
     */
    void cpToHost(size_t thisSize)
    {
        #ifdef CUDA
        if (accType == acc::cuda)
        {
            #ifdef DEBUG_CUDA
            if (!dPtr)
                ACC_PTR_DEBUG_FATAL("cp_to_host(thisSize) called before device allocation.\n");
            if (!hPtr)
                ACC_PTR_DEBUG_FATAL("null host pointer in cp_to_host(thisSize).\n");
            #endif
            cudaCpyDeviceToHost<T>(dPtr, hPtr, thisSize, stream);
        }
        #endif
    }

    /**
     * Copy a number (thisSize) of bytes from device to a specific host pointer
     */
    void cpToHost(T* hstPtr, size_t thisSize)
    {
        #ifdef CUDA
        if (accType == acc::cuda)
        {
            #ifdef DEBUG_CUDA
            if (!dPtr)
                ACC_PTR_DEBUG_FATAL("cp_to_host(hstPtr, thisSize) called before device allocation.\n");
            if (!hstPtr)
                ACC_PTR_DEBUG_FATAL("null host pointer in cp_to_host(hstPtr, thisSize).\n");
            #endif
            cudaCpyDeviceToHost<T>(dPtr, hstPtr, thisSize, stream);
        }
        #endif
    }

    /**
     * Copy a number (size) of bytes from device to the host pointer
     */
    void cpToHostOnStream(StreamType s)
    {
        #ifdef CUDA
        if (accType == acc::cuda)
        {
            #ifdef DEBUG_CUDA
            if (!dPtr)
                ACC_PTR_DEBUG_FATAL("cp_to_host_on_stream(s) called before device allocation.\n");
            if (!hPtr)
                ACC_PTR_DEBUG_FATAL("null host pointer in cp_to_host_on_stream(s).\n");
            #endif
            cudaCpyDeviceToHost<T>(dPtr, hPtr, size, s);
        }
        #endif
    }

    /**
     * Copy a number (size) of bytes from device pointer to the provided new device pointer
     */
    void cpOnDevice(T * const dstDevPtr)
    {
        #ifdef CUDA
        if (accType == acc::cuda)
        {
            #ifdef DEBUG_CUDA
            if (!dstDevPtr)
                ACC_PTR_DEBUG_FATAL("nullptr given in cpOnDevice(dstDevPtr).\n");
            #endif
            CudaShortcuts::cpyDeviceToDevice(dPtr, dstDevPtr, size, stream);
        }
        #endif
    }

    /**
     * Copy bytes into another host pointer
     */
    void cpOnHost(T * const dstDevPtr)
    {
        #ifdef DEBUG_CUDA
        if (!dstDevPtr)
            ACC_PTR_DEBUG_FATAL("nullptr given in cp_on_host(dstDevPtr).\n");
        if (!hPtr)
            ACC_PTR_DEBUG_FATAL("Null input pointer given in cp_on_host(hPtr).\n");
        #endif
        memcpy(dstDevPtr, hPtr, size * sizeof(T));
    }

    void cpOnAcc(T * const dstDevPtr)
    {
        if (accType == acc::cuda)
            cpOnDevice(dstDevPtr);
        else
            cpOnHost(dstDevPtr);
    }

    inline void cpOnAcc(const AccPtr<T> &devPtr)
    {
        cpOnAcc(devPtr.dPtr);
    }

    /**
     * Host data quick access
     */
    const T& operator[](size_t idx) const
    {
        #ifdef DEBUG_CUDA
        if (!hPtr)
            ACC_PTR_DEBUG_FATAL("const operator[] called with null host pointer.\n");
        #endif
        return hPtr[idx];
    };

    /**
     * Host data quick access
     */
    T& operator[](size_t idx)
    {
        #ifdef DEBUG_CUDA
        if (!hPtr)
            ACC_PTR_DEBUG_FATAL("operator[] called with null host pointer.\n");
        #endif
        return hPtr[idx];
    };

    /**
     * Device data quick access
     */
    T& operator()(size_t idx)
    {
        #ifdef DEBUG_CUDA
        if (!dPtr)
            ACC_PTR_DEBUG_FATAL("operator(idx) called with null acc pointer.\n");
        #endif
        return dPtr[idx];
    };


    /**
     * Device data quick access
     */
    const T& operator()(size_t idx) const
    {
        #ifdef DEBUG_CUDA
        if (!dPtr)
            ACC_PTR_DEBUG_FATAL("operator(idx) called with null acc pointer.\n");
        #endif
        return dPtr[idx];
    };

    /**
     * Raw data pointer quick access
     */
    T* operator()()
    {
        // TODO - this could cause considerable confusion given the above operators.  But it
        // also simplifies code that uses it.   What to do...
        if (accType == acc::cuda)
        {
            #ifdef DEBUG_CUDA
            if (!dPtr)
                ACC_PTR_DEBUG_FATAL("operator() called with null device pointer.\n");
            #endif
            return dPtr;
        }
        else
        {
            #ifdef DEBUG_CUDA
            if (!hPtr)
                ACC_PTR_DEBUG_FATAL("operator() called with null host pointer.\n");
            #endif
            return hPtr;
        }
    };

    T* operator~()
    {
        // TODO - this could cause considerable confusion given the above operators.  But it
        // also simplifies code that uses it.   What to do...
        if (accType == acc::cuda)
        {
            #ifdef DEBUG_CUDA
            if (!dPtr)
                ACC_PTR_DEBUG_FATAL("DEBUG_WARNING: \"kernel cast\" on null device pointer.\n");
            #endif
            return dPtr;
        }
        else
        {
            #ifdef DEBUG_CUDA
            if (!hPtr)
                ACC_PTR_DEBUG_FATAL("DEBUG_WARNING: \"kernel cast\" on null host pointer.\n");
            #endif
        return hPtr;
        }
    }

    void streamSync()
    {
        #ifdef CUDA
        if (accType == acc::cuda)
            DEBUG_HANDLE_ERROR(cudaStreamSynchronize(stream));
        #endif
    }

    T getAccValueAt(size_t idx)
    {
        #ifdef CUDA
        if (accType == acc::cuda)
        {
            T value;
            cudaCpyDeviceToHost<T>(&dPtr[idx], &value, 1, stream);
            streamSync();
            return value;
        }
        else
        #endif
            return hPtr[idx];
    }

    T getDeviceAt(size_t idx)
    {
        #ifdef CUDA
        if (accType == acc::cuda)
        {
            T value;
            cudaCpyDeviceToHost<T>(&dPtr[idx], &value, 1, stream);
            streamSync();
            return value;
        }
        #else
        return nullptr;
        #endif
    }

    void dumpDeviceToFile(const std::string &fileName)
    {

        #ifdef CUDA
        if (accType == acc::cuda)
        {
            T *tmp = new T[size];
            cudaCpyDeviceToHost<T>(dPtr, tmp, size, stream);

            std::ofstream f (fileName.c_str());
            streamSync();
            for (unsigned i = 0; i < size; i ++)
                f << tmp[i] << std::endl;
            delete [] tmp;
        }
        else
        #endif
        {
            std::ofstream f (fileName.c_str());
            f << "Pointer has no device support." << std::endl;
        }
    }

    void dumpHostToFile(const std::string &fileName)
    {
        std::ofstream f (fileName.c_str());
        for (unsigned i = 0; i < size; i ++)
            f << hPtr[i] << std::endl;
    }

    void dumpAccToFile(const std::string &fileName)
    {
        if (accType == acc::cuda)
            dumpDeviceToFile(fileName);
        else
            dumpHostToFile(fileName);
    }

    /**
     * Delete device data
     */
    void freeDevice()
    {
        if (!doFreeDevice) return;
        #ifdef CUDA
        if (accType == acc::cuda)
        {
            #ifdef DEBUG_CUDA
            if (!dPtr)
                ACC_PTR_DEBUG_FATAL("Free device memory was called on nullptr in free_device().\n");
            #endif
            doFreeDevice = false;

            if (alloc->getReadyEvent() == 0)
                alloc->markReadyEvent(stream);
            alloc->doFreeWhenReady();
            alloc = nullptr;

			// DEBUG_HANDLE_ERROR(cudaFree(dPtr));

            dPtr = nullptr;
        }
        #endif
    }

    /**
     * Delete host data
     */
    void freeHost()
    {
        if (!doFreeHost) return;
        #ifdef DEBUG_CUDA
        if (!hPtr)
            ACC_PTR_DEBUG_FATAL("free_host() called on nullptr.\n");
        #endif
        doFreeHost = false;
        ::free(hPtr);  // free(nullptr) is a no-op
        hPtr = nullptr;
    }

    void free()
    {
        freeHost();
        freeDevice();
    }

    ~AccPtr()
    {
        free();
    }

    /*======================================================
                       GETTERS AND SETTERS
    ======================================================*/

    inline bool willFreeHost() const
    {
        return doFreeHost;
    }

    inline bool willFreeDevice() const
    {
        return doFreeDevice;
    }

    inline void setStream(StreamType s)
    {
        stream = s;
    }

    inline StreamType getStream() const
    {
        return stream;
    }

    inline void setSize(size_t s)
    {
        size = s;
    }

    inline size_t getSize() const
    {
        return size;
    }

    inline T *getDevicePtr() const
    {
        return dPtr;
    }

    inline T *getHostPtr() const
    {
        return hPtr;
    }

    inline T *getAccPtr() const
    {
        if (accType == acc::cuda)
            return dPtr;
        else
            return hPtr;
    }

    void setAllocator(AllocatorType *a)
    {
        freeDevice();
        allocator = a;
    };

    inline AllocatorType *getAllocator() const
    {
        return allocator;
    }

    void setDevicePtr(T *ptr)
    {
        #ifdef DEBUG_CUDA
        if (doFreeDevice)
            ACC_PTR_DEBUG_FATAL("Device pointer set without freeing the old one.\n");
        #endif
        dPtr = ptr;
    }

    void setDevicePtr(const AccPtr<T> &ptr)
    {
        #ifdef DEBUG_CUDA
        if (!ptr.dPtr)
            ACC_PTR_DEBUG_FATAL("Device pointer is not set.\n");
        #endif
        setDevicePtr(ptr.dPtr);
    }

    void setHostPtr(T *ptr)
    {
        #ifdef DEBUG_CUDA
        if (doFreeHost)
            ACC_PTR_DEBUG_FATAL("Host pointer set without freeing the old one.\n");
        #endif
        hPtr = ptr;
    }

    void setHostPtr(const AccPtr<T> &ptr)
    {
        #ifdef DEBUG_CUDA
        if (!ptr.hPtr)
            ACC_PTR_DEBUG_FATAL("Host pointer is not set.\n");
        #endif
        setHostPtr(ptr.hPtr);
    }

    inline void setAccPtr(const AccPtr<T> &ptr)
    {
        setAccPtr(ptr.hPtr);
    }

    void setAccPtr(T *ptr)
    {
        if (accType == acc::cuda)
            setDevicePtr(ptr);
        else
            setHostPtr(ptr);
    }

    template <typename U>
    AccPtr<U> make(size_t s = 0)
    {
        return AccPtr<U>(allocator, stream, accType, s);
    }

};

using AccPtrBundleByte = unsigned char;

class AccPtrBundle: public AccPtr<AccPtrBundleByte>
{
    private:

    size_t current_packed_pos;

    public:

    AccPtrBundle(StreamType stream, AllocatorType *allocator, size_t size = 0, acc::Type accType = acc::type):
        AccPtr<AccPtrBundleByte>(allocator, stream, accType, size),
        current_packed_pos(0)
    {}

    template <typename T>
    void pack(AccPtr<T> &ptr)
    {
        #ifdef CUDA
        #ifdef DEBUG_CUDA
        if (current_packed_pos + ptr.getSize() > size)
            ACC_PTR_DEBUG_FATAL("Packing exceeds bundle total size.\n");
        if (!hPtr)
            ACC_PTR_DEBUG_FATAL("Pack called on null host pointer.\n");
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

    void allAlloc()
    {
        #ifdef CUDA
        AccPtr<AccPtrBundleByte>::allAlloc();
        #endif
    }

    void hostAlloc()
    {
        #ifdef CUDA
        AccPtr<AccPtrBundleByte>::hostAlloc();
        #endif
    }

};

class AccPtrFactory
{

    public:

    AccPtrFactory(acc::Type accT = acc::unset):  // Why not acc::cuda?
        allocator(nullptr), stream(0), accType(accT)
    {}

    AccPtrFactory(AllocatorType *alloc, StreamType stream = 0):
        allocator(alloc), stream(stream), accType(acc::cuda)
    {}

    template <typename T>
    AccPtr<T> make(size_t size = 0, StreamType stream = cudaStreamPerThread)
    {
        return AccPtr<T>(size, allocator, stream, accType);
    }

    AccPtrBundle makeBundle(size_t size = 0)
    {
        return AccPtrBundle(stream, allocator, size, accType);
    }

    private:

    AllocatorType *allocator;
    StreamType stream;
    acc::Type accType;

};

#endif
