#ifndef ACC_DATA_TYPES_H_
#define ACC_DATA_TYPES_H_

#include "src/acc/acc_ptr.h"
#include "src/multidim_array.h"

namespace AccDataTypes {

template<typename T>
class Image: public AccPtr<T> {

    private:

    int x, y, z;
    bool fourier; // Is this a Fourier space data array?

    public:

    /*======================================================
                         CONSTRUCTORS
    ======================================================*/

    template <acc::Type accType>
    Image(AccPtrFactory<accType> &f):
        AccPtr<T>(f.template make<T>()),
        x(0), y(0), z(0), fourier(false) {}

    template <acc::Type accType>
    Image(int xdim, AccPtrFactory<accType> &f):
        AccPtr<T>(f.template make<T>(xdim)),
        x(xdim), y(1), z(1), fourier(false) {}

    template <acc::Type accType>
    Image(int xdim, int ydim, AccPtrFactory<accType> &f):
        AccPtr<T>(f.template make<T>(xdim * ydim)),
        x(xdim), y(ydim), z(1), fourier(false) {}

    template <acc::Type accType>
    Image(int xdim, int ydim, int zdim, AccPtrFactory<accType> &f):
        AccPtr<T>(f.template make<T>(xdim * ydim * zdim)),
        x(xdim), y(ydim), z(zdim), fourier(false) {}

    template<typename T1, acc::Type accType>
    Image(MultidimArray<T1> img, AccPtrFactory<accType> &f):
        AccPtr<T>(f.template make<T>(img.size())),
        x(img.xdim), y(img.ydim), z(img.zdim),
        fourier(false) {}

    template <acc::Type accType>
    Image(int box_dim, bool is_fourier, bool is3D, AccPtrFactory<accType> &f) {
        setSize(box_dim, is_fourier, is3D);
        AccPtr<T>(f.template make<T>(x * y * z));
    }

    /*======================================================
                           METHODS
    ======================================================*/

    int getx() { return x; }
    int gety() { return y; }
    int getz() { return z; }
    int getxy()  { return x * y; }
    int getxyz() { return AccPtr<T>::getSize(); }

    bool is3D() { return z > 1; }

    void setSize(int box_dim, bool is_fourier, bool is3D) {
        if (fourier = is_fourier) {
            x = box_dim / 2 + 1;
            y = box_dim;
            z = is3D ? box_dim : 1;
        }
        AccPtr<T>::setSize(x * y * z);
    }

    void setSize(int xdim, int ydim = 1, int zdim = 1) {
        AccPtr<T>::setSize(xdim * ydim * zdim);
    }

    template <typename T1>
    void setSize(MultidimArray<T1> img) {
        AccPtr<T>::setSize(img.xdim * img.xdim * img.xdim);
    }

    template <typename T1>
    void setHost(MultidimArray<T1> &img) {
        if (img.xdim != x || img.ydim != y || img.zdim != z) {
            if (img.size() > AccPtr<T>::getSize()) {
                AccPtr<T>::free();
                setSize(img);
                AccPtr<T>::hostAlloc();
            } else {
                setSize(img);
            }
        }

        if (!AccPtr<T>::getHostPtr()) AccPtr<T>::hostAlloc();

        T *ptr = AccPtr<T>::getHostPtr();

        if (sizeof(T) == sizeof(T1)) {
            memcpy(ptr, img.data, sizeof(T) * img.size());
        } else {
            for (unsigned long i = 0; i < img.size(); i++)
                ptr[i] = (T) img.data[i];
        }
    }

    template <typename T1>
    void getHost(MultidimArray<T1> &img) {

        if (img.size() != AccPtr<T>::getSize()) {
            if (img.size() == 0) {
                img.resize(x, y, z);
            } else {
                CRITICAL("Trying to fill host-array with data from an array with different size!")
            }
        }
        T *ptr = AccPtr<T>::getHostPtr();

        if (sizeof(T) == sizeof(T1)) {
            memcpy(img.data, ptr, sizeof(T) * img.size());
        } else {
            for (unsigned long i = 0; i < img.size(); i++)
                img.data[i] = (T1) ptr[i];
        }
    }

};

}

#endif
