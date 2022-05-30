#pragma once

struct rgb_t {
    unsigned char r, g, b;
};

struct ColourScheme {

    typedef unsigned char (*grey_ft)(const rgb_t);

    typedef rgb_t (*rgb_ft)(const unsigned char);

    ColourScheme(rgb_ft greyToRGB = NULL, grey_ft rgbToGrey = NULL): 
    _greyToRGB(greyToRGB), _rgbToGrey(rgbToGrey) {}

    unsigned char rgbToGrey(const rgb_t rgb) {
        if (_rgbToGrey)
        return _rgbToGrey(rgb);
        throw std::exception();
    }

    rgb_t greyToRGB(const unsigned char grey) {
        if (_greyToRGB)
        return _greyToRGB(grey);
        throw std::exception();
    }

    private:

    grey_ft _rgbToGrey;
    rgb_ft _greyToRGB;

};

extern ColourScheme greyscale, rainbow, cyan_black_yellow, blue_grey_red, blue_grey_white, black_grey_red;
