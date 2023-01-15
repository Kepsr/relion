#ifndef ACC_SETTINGS_H_
#define ACC_SETTINGS_H_

#include "src/macros.h"

#ifdef ACC_DOUBLE_PRECISION

using XFLOAT = double;
#ifndef CUDA
typedef struct {
    XFLOAT x;
    XFLOAT y;
} double2;
#endif
namespace acc {
    using xfloat = double;
    // using xfloat2 = double2;
    using Complex = double2;
}

#else

using XFLOAT = float;
#ifndef CUDA
typedef struct {
    XFLOAT x;
    XFLOAT y;
} float2;
#endif
namespace acc {
    using xfloat = float;
    // using xfloat2 = float2;
    using Complex = float2;
}

#endif

#ifdef ALTCPU
#ifndef CUDA
typedef double CudaCustomAllocator;
#define cudaStreamPerThread 0
#endif
#endif

#endif /* ACC_SETTINGS_H_ */
