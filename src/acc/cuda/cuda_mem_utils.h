#ifndef CUDA_DEVICE_MEM_UTILS_H_
#define CUDA_DEVICE_MEM_UTILS_H_

#ifdef CUDA
#include <cuda_runtime.h>
#include <curand.h>
#include "src/acc/cuda/cuda_settings.h"
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
#include "src/complex.h"

/**
 * Print cuda device memory info
 */
static void cudaPrintMemInfo() {
    size_t free, total;
    DEBUG_HANDLE_ERROR(cudaMemGetInfo(&free, &total));
	const float MiB = 1024 * 1024;
    const float free_hr  = free  / MiB;
    const float total_hr = total / MiB;
    printf(
        "free %.2fMiB, total %.2fMiB, used %.2fMiB\n", 
        free_hr, total_hr, total_hr - free_hr
    );
}

template <typename T>
static inline
void cudaCpyHostToDevice(T *src, T *dest, size_t size) {
    DEBUG_HANDLE_ERROR(cudaMemcpy(dest, src, size * sizeof(T), cudaMemcpyHostToDevice));
};

template <typename T>
static inline
void cudaCpyHostToDevice(T *src, T *dest, size_t size, cudaStream_t stream) {
    DEBUG_HANDLE_ERROR(cudaMemcpyAsync(dest, src, size * sizeof(T), cudaMemcpyHostToDevice, stream));
};

template <typename T>
static inline
void cudaCpyDeviceToHost(T *src, T *dest, size_t size) {
    DEBUG_HANDLE_ERROR(cudaMemcpy(dest, src, size * sizeof(T), cudaMemcpyDeviceToHost));
};

template <typename T>
static inline
void cudaCpyDeviceToHost(T *src, T *dest, size_t size, cudaStream_t stream) {
    DEBUG_HANDLE_ERROR(cudaMemcpyAsync(dest, src, size * sizeof(T), cudaMemcpyDeviceToHost, stream));
};

template <typename T>
static inline
void cudaCpyDeviceToDevice(T *src, T *dest, size_t size, cudaStream_t stream) {
    DEBUG_HANDLE_ERROR(cudaMemcpyAsync(dest, src, size * sizeof(T), cudaMemcpyDeviceToDevice, stream));
};

template <typename T>
static inline
void cudaMemInit(T *ptr, T value, size_t size) {
    DEBUG_HANDLE_ERROR(cudaMemset(ptr, value, size * sizeof(T)));
};

template <typename T>
static inline
void cudaMemInit(T *ptr, T value, size_t size, cudaStream_t &stream) {
    DEBUG_HANDLE_ERROR(cudaMemsetAsync(ptr, value, size * sizeof(T), stream));
};

#endif
