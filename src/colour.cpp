#include <cmath>
#include <algorithm>
#include "src/colour.h"

ColourScheme black_grey_red = ColourScheme([] (const ColourScheme::grey_t grey) -> ColourScheme::rgb_t {
    ColourScheme::grey_t red, green, blue;
    if (grey >= 128) {
        red = 255;
        blue = green = floor((float) (255 - grey) * 2);
    } else {
        red = green = blue = floor((float) (grey * 2.0));
    }
    return {red, green, blue};
}, [] (const ColourScheme::rgb_t rgb) -> ColourScheme::grey_t {
    return floor(rgb.r == 255 ? 255.f - rgb.b / 2.f : rgb.r / 2.f);
});

ColourScheme blue_grey_white = ColourScheme([] (const ColourScheme::grey_t grey) -> ColourScheme::rgb_t {
    ColourScheme::grey_t red, green, blue;
    if (grey >= 128) {
        red = green = blue = floor((float) (2 * (grey - 128)));
    } else {
        red = 0;
        blue = green = floor((float) (255 - 2 * grey));
    }
    return {red, green, blue};
}, [] (const ColourScheme::rgb_t rgb) -> ColourScheme::grey_t {
    return floor(rgb.r == 0 ? 255.f - rgb.b / 2.f : rgb.r / 2.f + 128.f);
});

ColourScheme blue_grey_red = ColourScheme([] (const ColourScheme::grey_t grey) -> ColourScheme::rgb_t {
    const float a = grey / 85.0; // group
    const int X = floor(a);	//this is the integer part
    const ColourScheme::grey_t Y = floor(255 * (a - X)); //fractional part from 0 to 255
    const ColourScheme::grey_t Z = 255 - Y;
    switch (X) {
        case 0:  return {   0, Z, Z };
        case 1:  return {   Y, Y, Y };
        case 2:  return { 255, Z, Z };
        case 3:  return { 255, 0, 0 };
        default: return {   0, 0, 0 };
    }
}, [] (const ColourScheme::rgb_t rgb) -> ColourScheme::grey_t {
    int X;
    ColourScheme::grey_t Y;
    if (rgb.r == 0) {
        X = 0;
        Y = 255 - rgb.b;
    } else if (rgb.r == 255) {
        X = 2;
        Y = 255 - rgb.b;
    } else {
        X = 1;
        Y = rgb.b;
    }
    return ceil(85.f * (X + Y / 256.f));
});

ColourScheme rainbow = ColourScheme([] (const ColourScheme::grey_t grey) -> ColourScheme::rgb_t {
    const float a = (255 - grey) / 64.0; //invert and group
    const int X = floor(a);
    const ColourScheme::grey_t Y = floor(255 * (a - X)); //fractional part from 0 to 255
    const ColourScheme::grey_t Z = 255 - Y;
    switch (X) {
        case 0:  return { 255,   Y,   0 };
        case 1:  return {   Z, 255,   0 };
        case 2:  return {   0, 255,   Y };
        case 3:  return {   0,   Z, 255 };
        case 4:  return {   0,   0, 255 };
        default: return {   0,   0,   0 };
    }
}, [] (const ColourScheme::rgb_t rgb) -> ColourScheme::grey_t {
    int X;
    ColourScheme::grey_t Y;
    if (rgb.r > 0) {
        if (rgb.r == 255) {
            X = 0;
            Y = rgb.g;
        } else {
            X = 1;
            Y = 255 - rgb.r;
        }
    } else if (rgb.g > 0) {
        if (rgb.g == 255) {
            X = 2;
            Y = rgb.b;
        } else {
            X = 3;
            Y = 255 - rgb.g;
        }
    } else {
        Y = 255; X = 4;
    }
    return 255 - ceil(64.f * (X + Y / 255.f));
});

ColourScheme cyan_black_yellow = ColourScheme([] (const ColourScheme::grey_t grey) -> ColourScheme::rgb_t {
    const float d_rb = 3 * (grey - 128);
    const float d_g  = 3 * (std::abs(grey - 128) - 42);
    ColourScheme::grey_t red   = floor(std::min(255.0f, std::max(0.0f,  d_rb)));
    ColourScheme::grey_t green = floor(std::min(255.0f, std::max(0.0f,  d_g )));
    ColourScheme::grey_t blue  = floor(std::min(255.0f, std::max(0.0f, -d_rb)));
    return {red, green, blue};
}, [] (const ColourScheme::rgb_t rgb) -> ColourScheme::grey_t {
    return floor(128 + (rgb.r > 0 ?
        (rgb.r < 255 ? +(rgb.r / 3.f) : +(rgb.g / 3.f) + 42) :
        (rgb.b < 255 ? -(rgb.b / 3.f) : -(rgb.g / 3.f) - 42))
    );
});

ColourScheme greyscale = ColourScheme([] (const ColourScheme::grey_t grey) -> ColourScheme::rgb_t {
    return {grey, grey, grey};
}, nullptr);
