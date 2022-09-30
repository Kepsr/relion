#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <vector>
#include <algorithm>
#include <omp.h>

#include <src/image.h>

#include <tiffio.h>

class EERRenderer {

    private:

    FileName fn_movie;
    TIFF *ftiff;

    bool ready;
    bool is_legacy;
    bool is_7bit;
    bool read_data;

    std::vector<long long> frame_starts, frame_sizes;
    unsigned char *buf;

    static const char EER_FOOTER_OK[];
    static const char EER_FOOTER_ERR[];
    static const int EER_IMAGE_WIDTH, EER_IMAGE_HEIGHT, EER_IMAGE_PIXELS;
    static const unsigned int EER_LEN_FOOTER;
    static const uint16_t TIFF_COMPRESSION_EER8bit, TIFF_COMPRESSION_EER7bit;

    int eer_upsampling;
    int nframes;
    int preread_start, preread_end;
    long long file_size;
    void readLegacy(FILE *fh);
    void lazyReadFrames();

    template <typename T>
    void render16K(MultidimArray<T> &image, std::vector<unsigned int> &positions, std::vector<unsigned char> &symbols, int n_electrons);

    template <typename T>
    void render8K(MultidimArray<T> &image, std::vector<unsigned int> &positions, std::vector<unsigned char> &symbols, int n_electrons);

    template <typename T>
    void render4K(MultidimArray<T> &image, std::vector<unsigned int> &positions, std::vector<unsigned char> &symbols, int n_electrons);

    static TIFFErrorHandler prevTIFFWarningHandler;

    public:

    EERRenderer();
    ~EERRenderer();

    // TODO: Implement proper copy constructors. Currently, they are disabled to prevent memory corruption.
    EERRenderer(const EERRenderer&) {
        REPORT_ERROR("Copy constructor for EERRenderer not implemented yet.");
    }

    EERRenderer& operator=(const EERRenderer&) {
        REPORT_ERROR("Copy assignment operator for EERRenderer not implemented yet.");
    }

    // Wrapper to the default TIFF warning handler to suppress EER private tag warnings
    static void TIFFWarningHandler(const char* module, const char* fmt, va_list ap);
    static void silenceTIFFWarnings();

    // 1-indexed
    void setFramesOfInterest(int start, int end) {
        if (is_legacy)
            return;

        if (read_data)
            REPORT_ERROR("Logic error in EERRenderer::setFramesOfInterest(). This must be set before rendering.");
        preread_start = start - 1;
        preread_end = end - 1;
    }

    void read(FileName _fn_movie, int eer_upsampling = 2);

    int getNFrames();
    int getWidth();
    int getHeight();

    // Frame indices are 1-indexed.
    // image is cleared.
    // This function is thread-safe (except for timing).
    // It is caller's responsibility to make sure type T does not overflow.
    template <typename T>
    long long renderFrames(int frame_start, int frame_end, MultidimArray<T> &image);

    // The gain reference for EER is not multiplicative! So the inverse is taken here.
    // 0 means defect.
    template <typename T>
    static void loadEERGain(FileName fn_gain, MultidimArray<T> &gain, int eer_upsampling) {
        const bool is_multiplicative = (fn_gain.getExtension() == "gain");
        if (is_multiplicative) {
            silenceTIFFWarnings();
            fn_gain += ":tif";
        }

        Image<T> original;
        original.read(fn_gain, true, 0, nullptr, true); // explicitly use the first page
        const int nx_in = Xsize(original());
        const int ny_in = Ysize(original());
        const long long size_out = EER_IMAGE_WIDTH * eer_upsampling;

        // Revert Y flip in TIFF reader
        if (is_multiplicative) {
            const int ylim = ny_in / 2;
            for (int y1 = 0; y1 < ylim; y1++) {
                const int y2 = ny_in - 1 - y1;
                for (int x = 0; x < nx_in; x++) {
                    std::swap(
                        direct::elem(original(), x, y1), 
                        direct::elem(original(), x, y2)
                    );
                }
            }
        } 

        if (
            eer_upsampling == 2 && nx_in == EER_IMAGE_WIDTH && ny_in == EER_IMAGE_HEIGHT
            // gain = 4K and grid = 8K
        ) {
            gain.initZeros(size_out, size_out);
            for (long int j = 0; j < Ysize(gain); j++)
            for (long int i = 0; i < Xsize(gain); i++)
                direct::elem(gain, i, j) = direct::elem(original(), i / 2, j / 2);
        } else if (
            (eer_upsampling == 1 && nx_in == EER_IMAGE_WIDTH && ny_in == EER_IMAGE_HEIGHT) || // gain = 4K and grid = 4K
            (eer_upsampling == 2 && nx_in == EER_IMAGE_WIDTH * 2 && ny_in == EER_IMAGE_HEIGHT * 2)  // gain = 8K and grid = 8K
        ) {
            gain = original();
        } else if (
            eer_upsampling == 1 && nx_in == EER_IMAGE_WIDTH * 2 && ny_in == EER_IMAGE_HEIGHT * 2
            // gain = 8K and grid = 4K
        ) {
            gain.initZeros(size_out, size_out);
            for (long int j = 0; j < Ysize(original()); j++)
            for (long int i = 0; i < Xsize(original()); i++)
                direct::elem(gain, i / 2, j / 2) += direct::elem(original(), i, j);
        } else {
            std::cerr << "Size of input gain: X = " << nx_in << " Y = " << ny_in << " Expected: X = " << size_out << " Y = " << size_out << std::endl;
            REPORT_ERROR("Invalid gain size in EERRenderer::upsampleEERGain()");
        }
        
        if (!is_multiplicative) {
            double avg = gain.sum() / (size_out * size_out);

            for (auto &x : gain) {
                if (x != 0) { x = avg / x; }
            }
        }
    }

    static bool isEER(FileName fn_movie) {
        FileName ext = fn_movie.getExtension();
        return ext == "eer" || ext == "ecc";
    }
};
