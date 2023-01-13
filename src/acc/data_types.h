#ifndef ACC_DATA_TYPES_H_
#define ACC_DATA_TYPES_H_

#include "src/acc/acc_ptr.h"
#include "src/multidim_array.h"

namespace AccDataTypes {

template<typename T, acc::Type accType = acc::type>
class Image {

    using ptr_t = AccPtr<T, accType>;
    using factory_t = AccPtrFactory<accType>;

    private:

    int x, y, z;
    bool fourier; // Is this a Fourier space data array?

    public:

    ptr_t ptr;

    /*======================================================
                         CONSTRUCTORS
    ======================================================*/

    Image(factory_t &f):
        x(0), y(0), z(0), fourier(false),
        ptr(f.template make<T>()) {}

    Image(factory_t &f, int xdim, int ydim = 1, int zdim = 1):
        x(xdim), y(ydim), z(zdim), fourier(false),
        ptr(f.template make<T>(xdim * ydim * zdim)) {}

    Image(factory_t &f, int box_dim, bool is_fourier, bool is3D) {
        setSize(box_dim, is_fourier, is3D);
        ptr(f.template make<T>(x * y * z));
    }

    /*======================================================
                           METHODS
    ======================================================*/

    int getx() const { return x; }
    int gety() const { return y; }
    int getz() const { return z; }
    int getxy() const { return x * y; }
    int getxyz() const { return ptr.getSize(); }
    bool is3D() const { return z > 1; }

    void setSize(int box_dim, bool is_fourier, bool is3D) {
        if (fourier = is_fourier) {
            x = box_dim / 2 + 1;
            y = box_dim;
            z = is3D ? box_dim : 1;
        }
        ptr.setSize(x * y * z);
    }

    void setSize(int xdim, int ydim = 1, int zdim = 1) {
        ptr.setSize(xdim * ydim * zdim);
    }

    template <typename T1>
    void setSize(MultidimArray<T1> img) {
        ptr.setSize(img.xdim * img.xdim * img.xdim);
    }

    template <typename T1>
    void setHost(MultidimArray<T1> &img) {
        if (img.xdim != x || img.ydim != y || img.zdim != z) {
            if (img.size() > this->ptr.getSize()) {
                this->ptr.free();
                setSize(img);
                this->ptr.hostAlloc();
            } else {
                setSize(img);
            }
        }

        if (!this->ptr.getHostPtr()) this->ptr.hostAlloc();

        T *ptr = this->ptr.getHostPtr();

        if (sizeof(T) == sizeof(T1)) {
            memcpy(ptr, img.data, sizeof(T) * img.size());
        } else {
            std::copy_n(img.data, img.size(), ptr);
        }
    }

    template <typename U>
    void getHost(MultidimArray<U> &img) {

        if (img.size() != this->ptr.getSize()) {
            if (img.size() == 0) {
                img.resize(x, y, z);
            } else {
                CRITICAL("Trying to fill host-array with data from an array with different size!")
            }
        }
        T *ptr = this->ptr.getHostPtr();

        if (sizeof(T) == sizeof(U)) {
            memcpy(img.data, ptr, sizeof(T) * img.size());
        } else {
            std::copy_n(ptr, img.size(), img.data);
        }
    }

};

}

#endif
