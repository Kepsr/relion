#ifndef CUDA_CUSTOM_ALLOCATOR_CUH_
#define CUDA_CUSTOM_ALLOCATOR_CUH_

#ifdef CUDA
#include <cuda_runtime.h>
#include "src/acc/cuda/cuda_settings.h"
#include "src/acc/cuda/cuda_mem_utils.h"
#endif

#include <signal.h>
#include <fstream>
#include <iostream>
#include <vector>
#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#include "src/macros.h"
#include "src/error.h"
#include "src/parallel.h"

#ifdef CUSTOM_ALLOCATOR_MEMGUARD
#include <execinfo.h>
#include <cxxabi.h>
#endif

#ifdef DUMP_CUSTOM_ALLOCATOR_ACTIVITY
#define CUSTOM_ALLOCATOR_REGION_NAME( name ) (fprintf(stderr, "\n%s", name))
#else
#define CUSTOM_ALLOCATOR_REGION_NAME( name ) //Do nothing
#endif

class CudaCustomAllocator {

    using byte_t = unsigned char;

    static const unsigned GUARD_SIZE =   4;
    static const byte_t GUARD_VALUE  = 145;
    static const int ALLOC_RETRY     = 500;

    public:

    class Alloc {

        friend class CudaCustomAllocator;

        private:

        byte_t *ptr;
        size_t size;
        Alloc *prev, *next;
        bool free, freeWhenReady;
        cudaEvent_t readyEvent;  // Event record used for auto free

        #ifdef CUSTOM_ALLOCATOR_MEMGUARD
        byte_t *guardPtr;
        void *backtrace[20];
        size_t backtraceSize;
        #endif

        Alloc(size_t size = 0, bool free = false):
            ptr(nullptr), prev(nullptr), next(nullptr),
            size(size), free(free), readyEvent(nullptr), freeWhenReady(false)
        {}

        ~Alloc() {
            prev = nullptr;
            next = nullptr;
            ptr  = nullptr;
            if (readyEvent)
                DEBUG_HANDLE_ERROR(cudaEventDestroy(readyEvent));
        }

        public:

        inline byte_t *getPtr() { return ptr; }

        inline size_t getSize() const { return size; }

        inline bool isFree() const { return free; }

        inline void doFreeWhenReady() { freeWhenReady = true; }

        inline cudaEvent_t getReadyEvent() const { return readyEvent; }

        inline void markReadyEvent(cudaStream_t stream = 0) {
            //TODO add a debug warning if event already set
            DEBUG_HANDLE_ERROR(cudaEventCreate(&readyEvent));
            DEBUG_HANDLE_ERROR(cudaEventRecord(readyEvent, stream));
        }

    };

    private:

    Alloc *first;
    size_t totalSize;
    size_t alignmentSize;

    bool cache;

    pthread_mutex_t mutex;

    // Look for the first suited space
    inline Alloc *_getFirstSuitedFree(size_t size) {
        for (Alloc *a = first; a; a = a->next)
            if (a->free && a->size >= size)
                return a;
        return nullptr;
    }

    // Free allocs with recorded ready events
    inline bool _syncReadyEvents() {
        bool anything_ready = false;
        for (Alloc *a = first; a; a = a->next) {
            if (!a->free && a->freeWhenReady && a->readyEvent) {
                DEBUG_HANDLE_ERROR(cudaEventSynchronize(a->readyEvent));
                anything_ready = true;
            }
        }
        return anything_ready;
    }

    // Free allocs with recorded ready events
    inline bool _freeReadyAllocs() {
        bool anything_freed = false;
        for (Alloc *curr = first, *next = curr->next; next; curr = next, next = curr->next) {
            if (!curr->free && curr->freeWhenReady && curr->readyEvent) {
                const cudaError_t e = cudaEventQuery(curr->readyEvent);
                if (e == cudaSuccess) {
                    _free(curr);
                    next = first;  // List modified, restart
                    anything_freed = true;
                } else if (e != cudaErrorNotReady) {
                    _printState();
                    HandleError(e, __FILE__, __LINE__ );
                }
            }
        }
        return anything_freed;
    }

    inline size_t _getTotalFreeSpace() {
        if (cache) {
            size_t total = 0;
            for (Alloc *a = first; a; a = a->next) {
                if (a->free) total += a->size;
            }
            return total;
        } else {
            size_t free, total;
            DEBUG_HANDLE_ERROR(cudaMemGetInfo(&free, &total));
            return free;
        }
    }

    inline size_t _getTotalUsedSpace() {
        size_t total = 0;
        for (Alloc *a = first; a; a = a->next) {
            if (!a->free) total += a->size;
        }
        return total;
    }

    size_t _getNumberOfAllocs() {
        size_t total = 0;
        for (Alloc *a = first; a; a = a->next) {
            if (!a->free) total++;
        }
        return total;
    }

    inline size_t _getLargestContinuousFreeSpace() {
        if (cache) {
            size_t largest = 0;
            for (Alloc *a = first; a; a = a->next) {
                if (a->free && a->size > largest)
                    largest = a->size;
            }
            return largest;
        } else {
            return _getTotalFreeSpace();
        }
    }

    inline void _printState() {
        size_t total = 0;
        for (Alloc *a = first; a; a = a->next) {
            if (a->free) {
                printf("[%luB] ", a->size);
            } else if (a->freeWhenReady) {
                printf("<%luB> ", a->size);
            } else {
                printf("(%luB) ", a->size);
            }
            total += a->size;
        }
        printf("= %luB\n", total);
        fflush(stdout);
    }

    inline void _free(Alloc* a) {
        // printf("free: %u ", a->size);
        // _printState();

        #ifdef CUSTOM_ALLOCATOR_MEMGUARD
        size_t guardCount = a->size - (a->guardPtr - a->ptr);
        byte_t *guards = new byte_t[guardCount];
        cudaStream_t stream = 0;
        cudaCpyDeviceToHost<byte_t>(a->guardPtr, guards, guardCount, stream);
        DEBUG_HANDLE_ERROR(cudaStreamSynchronize(stream));
        for (int i = 0; i < guardCount; i ++)
            if (guards[i] != GUARD_VALUE) {
                fprintf (stderr, "ERROR: CORRUPTED byte_t GUARDS DETECTED\n");

                char ** messages = backtrace_symbols(a->backtrace, a->backtraceSize);

                // skip first stack frame (points here)
                for (int i = 1; i < a->backtraceSize && messages; ++i) {
                    char *mangled_name = 0, *offset_begin = 0, *offset_end = 0;

                    // find parantheses and +address offset surrounding mangled name
                    for (char *p = messages[i]; *p; ++p) {
                        if (*p == '(') {
                            mangled_name = p;
                        } else if (*p == '+') {
                            offset_begin = p;
                        } else if (*p == ')') {
                            offset_end = p;
                            break;
                        }
                    }

                    // if the line could be processed, attempt to demangle the symbol
                    if (mangled_name && offset_begin && offset_end && mangled_name < offset_begin) {
                        *mangled_name++ = '\0';
                        *offset_begin++ = '\0';
                        *offset_end++   = '\0';

                        int status;
                        char *real_name = abi::__cxa_demangle(mangled_name, 0, 0, &status);

                        // if demangling is successful, output the demangled function name
                        if (status == 0) {
                            std::cerr << "[bt]: (" << i << ") " << messages[i] << " : "
                                      << real_name << "+" << offset_begin << offset_end
                                      << std::endl;

                        } else {
                            // otherwise, output the mangled function name
                            std::cerr << "[bt]: (" << i << ") " << messages[i] << " : "
                                      << mangled_name << "+" << offset_begin << offset_end
                                      << std::endl;
                        }
						// free(real_name);
                    } else {
                        // otherwise, print the whole line
                        std::cerr << "[bt]: (" << i << ") " << messages[i] << std::endl;
                    }
                }
                std::cerr << std::endl;

				// free(messages);

                exit(EXIT_FAILURE);
            }
        delete[] guards;
        #endif

        a->free = true;

        if (cache) {
            //Previous neighbor is free, concatenate
            if (a->prev && a->prev->free) {
                //Resize and set pointer
                a->size += a->prev->size;
                a->ptr = a->prev->ptr;

                // Fetch secondary neighbor
                Alloc *ppL = a->prev->prev;

                // Remove primary neighbor
                if (!ppL) { //If the previous is first in chain
                    first = a;
                } else {
                    ppL->next = a;
                }

                delete a->prev;

                // Attach secondary neighbor
                a->prev = ppL;
            }

            // Next neighbor is free, concatenate
            if (a->next && a->next->free) {
                //Resize and set pointer
                a->size += a->next->size;

                // Fetch secondary neighbor
                Alloc *nnL = a->next->next;

                // Remove primary neighbor
                if (nnL) { nnL->prev = a; }
                delete a->next;

                // Attach secondary neighbor
                a->next = nnL;
            }
        } else {
            DEBUG_HANDLE_ERROR(cudaFree(a->ptr ));
            a->ptr = nullptr;

            if (a->prev) {
                a->prev->next = a->next;
            } else {
                first = a->next;  // This is the first link
            }

            if (a->next) { a->next->prev = a->prev; }

            delete a;
        }
    };

    void _clear() {
        DEBUG_HANDLE_ERROR(cudaFree(first->ptr));
        first->ptr = nullptr;
        for (Alloc *a = first, *nL = a->next; a; a = nL, nL = a->next) {
            delete a;
        }
    }

    public:

    CudaCustomAllocator(size_t size, size_t alignment):
        totalSize(size), alignmentSize(alignment),
        first(new Alloc(size, true)), cache(size > 0), mutex()
    {
        if (cache)
            HANDLE_ERROR(cudaMalloc((void**) &first->ptr, totalSize));
        const int mutex_error = pthread_mutex_init(&mutex, nullptr);
        if (mutex_error) {
            printf("ERROR: Mutex could not be created for alloactor. CODE: %d.\n", mutex_error);
            fflush(stdout);
            CRITICAL(ERR_CAMUX);
        }
    }

    void resize(size_t size) {
        Lock ml (&mutex);
        _clear();
        totalSize = size;
        first = new Alloc(totalSize, true);
        if (cache = totalSize > 0) {
            HANDLE_ERROR(cudaMalloc((void**) &first->ptr, totalSize));
        }
    }

    inline Alloc* alloc(size_t requestedSize) {
        Lock ml (&mutex);

        _freeReadyAllocs();

		// printf("alloc: %u ", size);
		// _printState();

        size_t size = requestedSize;

        #ifdef CUSTOM_ALLOCATOR_MEMGUARD
        //Ad byte-guards
        size += alignmentSize * GUARD_SIZE; //Ad an integer multiple of alignment size as byte guard size
        #endif

        #ifdef DUMP_CUSTOM_ALLOCATOR_ACTIVITY
        fprintf(stderr, " %.4f", 100.0 * (float) size / (float) totalSize);
        #endif

        Alloc *newAlloc = nullptr;

        if (cache) {
            // To prevent misaligned memory
            size = alignmentSize * ceilf((float) size / (float) alignmentSize);

            Alloc *curAlloc = _getFirstSuitedFree(size);

            // If out of memory
            if (!curAlloc) {
                #ifdef DEBUG_CUDA
                size_t spaceDiff = _getTotalFreeSpace();
                #endif
                // Try to recover before throwing error
                for (int i = 0; i <= ALLOC_RETRY; i++) {
                    if (_syncReadyEvents() && _freeReadyAllocs()) {
                        curAlloc = _getFirstSuitedFree(size); //Is there space now?
                        if (curAlloc) break; //Success
                    } else {
                        usleep(10000); // 10 ms, Order of magnitude of largest kernels
                    }
                }
                #ifdef DEBUG_CUDA
                spaceDiff =  _getTotalFreeSpace() - spaceDiff;
                printf("DEBUG_INFO: Out of memory handled by waiting for unfinished tasks, which freed %lu B.\n", spaceDiff);
                #endif

                // Did we manage to recover?
                if (!curAlloc) {
                    printf(
                        "ERROR: CudaCustomAllocator out of memory\n [requestedSpace:             %lu B]\n [largestContinuousFreeSpace: %lu B]\n [totalFreeSpace:             %lu B]\n",
                        (unsigned long) size, (unsigned long) _getLargestContinuousFreeSpace(), (unsigned long) _getTotalFreeSpace()
                    );

                    _printState();

                    fflush(stdout);
                    CRITICAL(ERRCUDACAOOM);
                }
            }

            if (curAlloc->size == size) {
                curAlloc->free = false;
                newAlloc = curAlloc;
            } else {  // Or curAlloc->size is smaller than size
                // Setup new pointer
                newAlloc = new Alloc(size, false);
                newAlloc->next = curAlloc;
                newAlloc->ptr = curAlloc->ptr;

                // Modify old pointer
                curAlloc->ptr = curAlloc->ptr + size;
                curAlloc->size -= size;

                // Insert new allocation region into chain
                if (!curAlloc->prev) {
                    // If the first allocation region
                    first = newAlloc;
                } else {
                    curAlloc->prev->next = newAlloc;
                }
                newAlloc->prev = curAlloc->prev;
                newAlloc->next = curAlloc;
                curAlloc->prev = newAlloc;
            }
        } else {
            newAlloc = new Alloc(size);
            DEBUG_HANDLE_ERROR(cudaMalloc((void**) &newAlloc->ptr, size));

            //Just add to start by replacing first
            newAlloc->next = first;
            first->prev = newAlloc;
            first = newAlloc;
        }

        #ifdef CUSTOM_ALLOCATOR_MEMGUARD
        newAlloc->backtraceSize = backtrace(newAlloc->backtrace, 20);
        newAlloc->guardPtr = newAlloc->ptr + requestedSize;
        cudaStream_t stream = 0;
        cudaMemInit<byte_t>(newAlloc->guardPtr, GUARD_VALUE, size - requestedSize, stream); //TODO switch to specialized stream
        DEBUG_HANDLE_ERROR(cudaStreamSynchronize(stream));
        #endif

        return newAlloc;
    };

    ~CudaCustomAllocator() {
        {
            Lock ml (&mutex);
            _clear();
        }
        pthread_mutex_destroy(&mutex);
    }

    // Thread-safe wrapper functions

    inline void free(Alloc* a) {
        Lock ml (&mutex);
        _free(a);
    }

    inline void syncReadyEvents() {
        Lock ml (&mutex);
        _syncReadyEvents();
    }

    inline void freeReadyAllocs() {
        Lock ml (&mutex);
        _freeReadyAllocs();
    }

    size_t getTotalFreeSpace() {
        Lock ml (&mutex);
        return _getTotalFreeSpace();
    }

    size_t getTotalUsedSpace() {
        Lock ml (&mutex);
        return _getTotalUsedSpace();
    }

    size_t getNumberOfAllocs() {
        Lock ml (&mutex);
        return _getNumberOfAllocs();
    }

    size_t getLargestContinuousFreeSpace() {
        Lock ml (&mutex);
        return _getLargestContinuousFreeSpace();
    }

    void printState() {
        Lock ml (&mutex);
        _printState();
    }

};

template <typename T, bool CustomAlloc=true>
class CudaGlobalPtr {

    CudaCustomAllocator *allocator;
    CudaCustomAllocator::Alloc *alloc;
    cudaStream_t stream;

    public:
    size_t size;  // Size used when copying data from and to device
    T *h_ptr, *d_ptr;  // Host and device pointers
    bool h_do_free, d_do_free;  // True if host or device needs to be freed

    /*======================================================
                CONSTRUCTORS WITH ALLOCATORS
    ======================================================*/

    inline CudaGlobalPtr(CudaCustomAllocator *allocator):
        size(0), h_ptr(0), d_ptr(0), h_do_free(false),
        d_do_free(false), allocator(allocator), alloc(0), stream(cudaStreamPerThread) {};

    inline CudaGlobalPtr(cudaStream_t stream, CudaCustomAllocator *allocator):
        size(0), h_ptr(0), d_ptr(0), h_do_free(false),
        d_do_free(false), allocator(allocator), alloc(0), stream(stream) {};

    inline CudaGlobalPtr(size_t size, CudaCustomAllocator *allocator):
        size(size), h_ptr(new T[size]), d_ptr(0), h_do_free(true),
        d_do_free(false), allocator(allocator), alloc(0), stream(cudaStreamPerThread) {};

    inline CudaGlobalPtr(size_t size, cudaStream_t stream, CudaCustomAllocator *allocator):
        size(size), h_ptr(new T[size]), d_ptr(0), h_do_free(true),
        d_do_free(false), allocator(allocator), alloc(0), stream(stream) {};

    inline CudaGlobalPtr(T * h_start, size_t size, CudaCustomAllocator *allocator):
        size(size), h_ptr(h_start), d_ptr(0), h_do_free(false),
        d_do_free(false), allocator(allocator), alloc(0), stream(cudaStreamPerThread) {};

    inline CudaGlobalPtr(T * h_start, size_t size, cudaStream_t stream, CudaCustomAllocator *allocator):
        size(size), h_ptr(h_start), d_ptr(0), h_do_free(false),
        d_do_free(false), allocator(allocator), alloc(0), stream(cudaStreamPerThread) {};

    inline CudaGlobalPtr(T * h_start, T * d_start, size_t size, CudaCustomAllocator *allocator):
        size(size), h_ptr(h_start), d_ptr(d_start), h_do_free(false),
        d_do_free(false), allocator(allocator), alloc(0), stream(cudaStreamPerThread) {};

    inline CudaGlobalPtr(T * h_start, T * d_start, size_t size, cudaStream_t stream, CudaCustomAllocator *allocator):
        size(size), h_ptr(h_start), d_ptr(d_start), h_do_free(false),
        d_do_free(false), allocator(allocator), alloc(0), stream(stream) {};

    /*======================================================
               CONSTRUCTORS WITHOUT ALLOCATORS
    ======================================================*/

    inline CudaGlobalPtr():
        size(0), h_ptr(0), d_ptr(0), h_do_free(false),
        d_do_free(false), allocator(0), alloc(0), stream(cudaStreamPerThread) {};

    inline CudaGlobalPtr(cudaStream_t stream):
        size(0), h_ptr(0), d_ptr(0), h_do_free(false),
        d_do_free(false), allocator(0), alloc(0), stream(stream) {};

    inline CudaGlobalPtr(size_t size):
        size(size), h_ptr(new T[size]), d_ptr(0), h_do_free(true),
        d_do_free(false), allocator(0), alloc(0), stream(cudaStreamPerThread) {};

    inline CudaGlobalPtr(size_t size, cudaStream_t stream):
        size(size), h_ptr(new T[size]), d_ptr(0), h_do_free(true),
        d_do_free(false), allocator(0), alloc(0), stream(stream) {};

    inline CudaGlobalPtr(T * h_start, size_t size):
        size(size), h_ptr(h_start), d_ptr(0), h_do_free(false),
        d_do_free(false), allocator(0), alloc(0), stream(0) {};

    inline CudaGlobalPtr(T * h_start, size_t size, cudaStream_t stream):
        size(size), h_ptr(h_start), d_ptr(0), h_do_free(false),
        d_do_free(false), allocator(0), alloc(0), stream(cudaStreamPerThread) {};

    inline CudaGlobalPtr(T * h_start, T * d_start, size_t size):
        size(size), h_ptr(h_start), d_ptr(d_start), h_do_free(false),
        d_do_free(false), allocator(0), alloc(0), stream(cudaStreamPerThread) {};

    inline CudaGlobalPtr(T * h_start, T * d_start, size_t size, cudaStream_t stream):
        size(size), h_ptr(h_start), d_ptr(d_start), h_do_free(false),
        d_do_free(false), allocator(0), alloc(0), stream(stream) {};

    /*======================================================
           CONSTRUCTORS WITH OTHER GLOBAL POINTERS
    ======================================================*/

    inline CudaGlobalPtr(const CudaGlobalPtr<T> &ptr):
        size(ptr.size), h_ptr(ptr.h_ptr), d_ptr(ptr.d_ptr), h_do_free(false),
        d_do_free(false), allocator(ptr.allocator), alloc(0), stream(ptr.stream) {};

    inline CudaGlobalPtr(const CudaGlobalPtr<T> &ptr, size_t start_idx, size_t size):
        size(size), h_ptr(&ptr.h_ptr[start_idx]), d_ptr(&ptr.d_ptr[start_idx]), h_do_free(false),
        d_do_free(false), allocator(ptr.allocator), alloc(0), stream(ptr.stream) {};


    /*======================================================
                        OTHER STUFF
    ======================================================*/

    CudaCustomAllocator *getAllocator() {return allocator; };
    cudaStream_t &getStream() {return stream; };
    void setStream(cudaStream_t s) { stream = s; };

    void setSize(size_t s) { size = s; };
    size_t getSize() { return size; };


    void setAllocator(CudaCustomAllocator *a) {
        free_device_if_set();
        allocator = a;
    };

    void markReadyEvent() {
        #ifdef DEBUG_CUDA
        if (!alloc) printf("DEBUG_WARNING: markReadyEvent called on null allocation.\n");
        #endif
        alloc->markReadyEvent(stream);
    }

    void setDevPtr(T *ptr) {
        #ifdef DEBUG_CUDA
        if (d_do_free) printf("DEBUG_WARNING: Device pointer set without freeing the old one.\n");
        #endif
        d_ptr = ptr;
    };

    void setDevPtr(const CudaGlobalPtr<T> &ptr) {
        #ifdef DEBUG_CUDA
        if (!ptr.d_ptr) printf("DEBUG_WARNING: Device pointer is not set.\n");
        #endif
        setHstPtr(ptr.d_ptr);
    };

    void setHstPtr(T *ptr) {
        #ifdef DEBUG_CUDA
        if (h_do_free) printf("DEBUG_WARNING: Host pointer set without freeing the old one.\n");
        #endif
        h_ptr = ptr;
    };

    void setHstPtr(const CudaGlobalPtr<T> &ptr) {
        #ifdef DEBUG_CUDA
        if (!ptr.h_ptr) printf("DEBUG_WARNING: Host pointer is not set.\n");
        #endif
        setHstPtr(ptr.h_ptr);
    };

    /**
     * Allocate memory on device
     */
    inline void device_alloc() {
        #ifdef DEBUG_CUDA
        if (size == 0) printf("DEBUG_WARNING: device_alloc called with size == 0");
        if (d_do_free) printf("DEBUG_WARNING: Device double allocation.\n");
        #endif
        d_do_free = true;
        if (CustomAlloc) {
            alloc = allocator->alloc(size * sizeof(T));
            d_ptr = (T*) alloc->getPtr();
        } else {
            DEBUG_HANDLE_ERROR(cudaMalloc((void**) &d_ptr, size * sizeof(T)));
        }
    }

    /**
     * Allocate memory on device with given size
     */
    inline void device_alloc(size_t newSize) {
        size = newSize;
        device_alloc();
    }

    /**
     * Allocate memory on host
     */
    inline void host_alloc() {
        #ifdef DEBUG_CUDA
        if (size == 0) printf("DEBUG_WARNING: device_alloc called with size == 0");
        if (h_do_free) printf("DEBUG_WARNING: Host double allocation.\n");
        #endif
        h_do_free = true;
        h_ptr = new T[size];
    }

    /**
     * Allocate memory on host with given size
     */
    inline void host_alloc(size_t newSize) {
        size = newSize;
        host_alloc();
    }

    void resize_host(size_t newSize) {
        #ifdef DEBUG_CUDA
        if (size == 0) printf("DEBUG_WARNING: Resizing from size zero (permitted).\n");
        #endif
        T* newArr = new T[newSize];
        memcpy(newArr, h_ptr, newSize * sizeof(T) );

        size = newSize;
        #ifdef DEBUG_CUDA
        if (d_ptr) printf("DEBUG_WARNING: Resizing host with present device allocation.\n");
        #endif
        free_host();
        setHstPtr(newArr);
        h_do_free = true;
    }

    /**
     * Initiate device memory with provided value
     */
    inline void device_init(int value) {
        #ifdef DEBUG_CUDA
        if (d_ptr == 0) printf("DEBUG_WARNING: Memset requested before allocation in device_init().\n");
        #endif
        cudaMemInit<T>(d_ptr, value, size, stream);
    }

    /**
     * Copy a number (size) of bytes to device stored in the host pointer
     */
    inline void cp_to_device() {
        #ifdef DEBUG_CUDA
        if (d_ptr == 0) printf("DEBUG_WARNING: cp_to_device() called before allocation.\n");
        if (h_ptr == 0) printf("DEBUG_WARNING: nullptr host pointer in cp_to_device().\n");
        #endif
        cudaCpyHostToDevice<T>(h_ptr, d_ptr, size, stream);
    }

    /**
     * Copy a number (size) of bytes to device stored in the provided host pointer
     */
    inline void cp_to_device(T * hostPtr) {
        #ifdef DEBUG_CUDA
        if (!hostPtr) printf("DEBUG_WARNING: Null-pointer given in cp_to_device(hostPtr).\n");
        #endif
        h_ptr = hostPtr;
        cp_to_device();
    }

    /**
     * Copy a number (size) of bytes from device pointer to the provided new device pointer
     */
    inline void cp_on_device(T * dstDevPtr) {
        #ifdef DEBUG_CUDA
        if (!dstDevPtr)
        printf("DEBUG_WARNING: Null-pointer given in cp_on_device(dstDevPtr).\n");
        #endif
        cudaCpyDeviceToDevice(d_ptr, dstDevPtr, size, stream);
    }

    /**
     * Copy a number (size) of bytes from device pointer to the provided new device pointer
     */
    inline void cp_on_device(CudaGlobalPtr<T> &devPtr) {
        #ifdef DEBUG_CUDA
        if (devPtr.size == 0) printf("DEBUG_WARNING: Zero size on provided pointer in cp_on_device.\n");
        #endif
        cp_on_device(devPtr.d_ptr);
    }

    /**
     * alloc and copy
     */
    inline void put_on_device() {
        device_alloc();
        cp_to_device();
    }

    /**
     * alloc size and copy
     */
    inline void put_on_device(size_t newSize) {
        size = newSize;
        device_alloc();
        cp_to_device();
    }


    /**
     * Copy a number (size) of bytes from device to the host pointer
     */
    inline void cp_to_host() {
        #ifdef DEBUG_CUDA
        if (!d_ptr) printf("DEBUG_WARNING: cp_to_host() called before device allocation.\n");
        if (!h_ptr) printf("DEBUG_WARNING: nullptr host pointer in cp_to_host().\n");
        #endif
        cudaCpyDeviceToHost<T>(d_ptr, h_ptr, size, stream);
    }

    /**
     * Copy a number (thisSize) of bytes from device to the host pointer
     */
    inline void cp_to_host(size_t thisSize) {
        #ifdef DEBUG_CUDA
        if (!d_ptr) printf("DEBUG_WARNING: cp_to_host(thisSize) called before device allocation.\n");
        if (!h_ptr) printf("DEBUG_WARNING: nullptr host pointer in cp_to_host(thisSize).\n");
        #endif
        cudaCpyDeviceToHost<T>(d_ptr, h_ptr, thisSize, stream);
    }

    /**
     * Copy a number (thisSize) of bytes from device to a specific host pointer
     */
    inline void cp_to_host(T* hstPtr, size_t thisSize) {
        #ifdef DEBUG_CUDA
        if (!d_ptr) printf("DEBUG_WARNING: cp_to_host(hstPtr, thisSize) called before device allocation.\n");
        if (!hstPtr) printf("DEBUG_WARNING: nullptr host pointer in cp_to_host(hstPtr, thisSize).\n");
        #endif
        cudaCpyDeviceToHost<T>(d_ptr, hstPtr, thisSize, stream);
    }

    /**
     * Copy a number (size) of bytes from device to the host pointer
     */
    inline void cp_to_host_on_stream(cudaStream_t s) {
        #ifdef DEBUG_CUDA
        if (!d_ptr) printf("DEBUG_WARNING: cp_to_host_on_stream(s) called before device allocation.\n");
        if (!h_ptr) printf("DEBUG_WARNING: nullptr host pointer in cp_to_host_on_stream(s).\n");
        #endif
        cudaCpyDeviceToHost<T>(d_ptr, h_ptr, size, s);
    }

    /**
     * Host data quick access
     */
    inline T& operator[](size_t idx) { return h_ptr[idx]; };


    /**
     * Host data quick access
     */
    inline const T& operator[](size_t idx) const { return h_ptr[idx]; };

    /**
     * Device data quick access
     */
    inline T& operator()(size_t idx) { return d_ptr[idx]; };


    /**
     * Device data quick access
     */
    inline const T& operator()(size_t idx) const { return d_ptr[idx]; };

    /**
     * Device pointer quick access
     */
    inline T* operator~() {
        #ifdef DEBUG_CUDA
        if (d_ptr == 0) printf("DEBUG_WARNING: \"kernel cast\" on null pointer.\n");
        #endif
        return d_ptr;
    };

    inline void streamSync() {
        DEBUG_HANDLE_ERROR(cudaStreamSynchronize(stream));
    }

    inline T getDeviceAt(size_t idx) {
        T value;
        cudaCpyDeviceToHost<T>(&d_ptr[idx], &value, 1, stream);
        streamSync();
        return value;
    }


    void dump_device_to_file(std::string fileName) {
        T *tmp = new T[size];
        cudaCpyDeviceToHost<T>(d_ptr, tmp, size, stream);

        std::ofstream f;
        f.open(fileName.c_str());
        streamSync();
        for (unsigned i = 0; i < size; i++) {
            f << tmp[i] << std::endl;
        }
        f.close();
        delete [] tmp;
    }

    void dump_host_to_file(std::string fileName) {
        std::ofstream f;
        f.open(fileName.c_str());
        for (unsigned i = 0; i < size; i ++)
            f << h_ptr[i] << std::endl;
        f.close();
    }

    /**
     * Delete device data
     */
    inline void free_device() {
        #ifdef DEBUG_CUDA
        if (d_ptr == 0)
            printf("DEBUG_WARNING: Free device memory was called on nullptr pointer in free_device().\n");
        #endif
        d_do_free = false;
        if (CustomAlloc) {
            if (alloc->getReadyEvent() == 0) 
                alloc->markReadyEvent(stream);
            alloc->doFreeWhenReady();

            alloc = nullptr;
        } else {
            DEBUG_HANDLE_ERROR(cudaFree(d_ptr));
        }
        d_ptr = 0;
    }

    /**
     * Delete host data
     */
    inline void free_host() {
        #ifdef DEBUG_CUDA
        if (h_ptr == 0) {
            printf("DEBUG_ERROR: free_host() called on nullptr pointer.\n");
            exit(EXIT_FAILURE );
        }
        #endif
        h_do_free = false;
        delete [] h_ptr;
        h_ptr = 0;
    }

    inline void free_host_if_set() {
        if (h_do_free)
            free_host();
    }

    inline void free_device_if_set() {
        if (d_do_free)
            free_device();
    }

    /**
     * Delete both device and host data
     */
    inline void free() {
        free_device();
        free_host();
    }

    inline void free_if_set() {
        free_host_if_set();
        free_device_if_set();
    }

    inline ~CudaGlobalPtr() {
        free_if_set();
    }

};

template <typename T>
class cudaStager {

    public:

    CudaGlobalPtr<T> AllData;
    size_t size;  // size of allocated host-space (AllData.size dictates the amount of memory copied to/from the device)

    /*======================================================
                CONSTRUCTORS WITH ALLOCATORS
    ======================================================*/

    inline cudaStager(CudaCustomAllocator *allocator):
        AllData(allocator), size(0) {};

    inline cudaStager(CudaCustomAllocator *allocator, size_t newSize):
        AllData(newSize,allocator), size(newSize) {
        AllData.size = 0;
    };

    /*======================================================
                CONSTRUCTORS WITHOUT ALLOCATORS
    ======================================================*/

    inline cudaStager(): AllData(), size(0) {};

    inline cudaStager(size_t newSize):
        AllData(newSize), size(newSize) {
        AllData.size = 0;
    };

    public:

    void prepare_host() {

        if (size == 0) {
            printf("trying to host-alloc a stager with size=0");
            CRITICAL(ERR_STAGEMEM);
        }
        size_t temp_size = AllData.size;
        AllData.size = size;
        if (!AllData.h_ptr) {
            AllData.host_alloc();
        } else {
            printf("WARNING : host_alloc when host-ptr is non-null");
        }
        AllData.size = temp_size;
    }

    void prepare_host(size_t alloc_size) {
        if (size == 0) {
            printf("trying to device-alloc a stager with size=0");
            CRITICAL(ERR_STAGEMEM);
        }
        size_t temp_size = AllData.size;
        AllData.size = alloc_size;
        if (!AllData.h_ptr) {
            AllData.host_alloc();
        } else {
            printf("WARNING : host_alloc when host-ptr is non-null");
        }
        AllData.size = temp_size;
    }

    void prepare_device() {
        if (size == 0) {
            printf("trying to host-alloc a stager with size=0");
            CRITICAL(ERR_STAGEMEM);
        }
        size_t temp_size = AllData.size;
        AllData.size = size;
        if (!AllData.d_ptr) {
            AllData.device_alloc();
        } else {
            printf("WARNING : device_alloc when dev-ptr is non-null");
        }
        AllData.size = temp_size;
    }

    void prepare_device(size_t alloc_size) {
        if (size == 0) {
            printf("trying to device-alloc a stager with size=0");
            CRITICAL(ERR_STAGEMEM);
        }
        size_t temp_size = AllData.size;
        AllData.size = alloc_size;
        if (!AllData.d_ptr) {
            AllData.device_alloc();
        } else {
            printf("WARNING : device_alloc when dev-ptr is non-null");
        }
        AllData.size = temp_size;
    }

    void prepare() {
         prepare_host();
         prepare_device();
    }

    void prepare(size_t alloc_size) {
         prepare_host(alloc_size);
         prepare_device(alloc_size);
    }

    void stage(CudaGlobalPtr<T> &input) {
        if (AllData.size + input.size > size) {
            printf("trying to stage more than stager can fit");
            printf(" (attempted to stage %lu addtionally to the allready staged %lu, and total host-allocated capacity is %lu ",input.size,AllData.size,size);
            exit(EXIT_FAILURE );
        }

        for (size_t i = 0 ; i < input.size; i++)
            AllData.h_ptr[AllData.size + i] = input.h_ptr[i];

        // reset the staged object to this new position (TODO: disable for pinned mem)
        if (input.h_ptr && input.h_do_free) input.free_host_if_set();
        input.h_ptr = &AllData.h_ptr[AllData.size];
        input.d_ptr = &AllData.d_ptr[AllData.size];

        AllData.size += input.size;
    }

    void cp_to_device() {
        AllData.cp_to_device();
    }

    void cp_to_host() {
        AllData.cp_to_host();
    }

};

#endif
