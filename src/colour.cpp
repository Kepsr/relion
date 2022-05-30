#include <cmath>
#include <algorithm>
#include "src/colour.h"

ColourScheme black_grey_red = ColourScheme([] (
    const unsigned char grey
) -> rgb_t {
    unsigned char red, green, blue;
    if (grey >= 128) {
        red = 255;
        blue = green = floor((float) (255 - grey) * 2);
    } else {
        red = green = blue = floor((float) (grey * 2.0));
    }
    return {red, green, blue};
}, [] (
    const unsigned char red,
    const unsigned char green,
    const unsigned char blue
) -> unsigned char {
    return floor(red == 255 ? (float) 255.0 - blue / 2.0 : (float) red / 2.0);
});

ColourScheme blue_grey_white = ColourScheme([] (
    const unsigned char grey
) -> rgb_t {
    unsigned char red, green, blue;
    if (grey >= 128) {
        red = green = blue = floor((float) (2 * (grey - 128)));
    } else {
        red = 0;
        blue = green = floor((float) (255 - 2 * grey));
    }
    return {red, green, blue};
}, [] (
    const unsigned char red, 
    const unsigned char green, 
    const unsigned char blue
) -> unsigned char {
    return floor(red == 0 ? (float) (255.0 - blue) / 2.0 : (float) red / 2.0 + 128.0);
});

ColourScheme blue_grey_red = ColourScheme([] (
    const unsigned char grey
) -> rgb_t {
    const float a = grey / 85.0; // group
    const int X = floor(a);	//this is the integer part
    const unsigned char Y = floor(255 * (a - X)); //fractional part from 0 to 255
    const unsigned char Z = 255 - Y;
    switch (X) {
        case 0: return {   0, Z, Z };
        case 1: return {   Y, Y, Y };
        case 2: return { 255, Z, Z };
        case 3: return { 255, 0, 0 };
    }
}, [] (
    const unsigned char red, 
    const unsigned char green, 
    const unsigned char blue
) -> unsigned char {
    int X;
    unsigned char Y;
    if (red == 0) {
        X = 0;
        Y = 255 - blue;
    } else if (red == 255) {
        X = 2;
        Y = 255 - blue;
    } else {
        X = 1;
        Y = blue;
    }
    return ceil(85 * (X + (float) Y / 256.0));
});

ColourScheme rainbow = ColourScheme([] (
    const unsigned char grey
) -> rgb_t {
    const float a = (255 - grey) / 64.0; //invert and group
    const int X = floor(a);
    const unsigned char Y = floor(255 * (a - X)); //fractional part from 0 to 255
    const unsigned char Z = 255 - Y;
    switch (X) {
        case 0: return { 255,   Y,   0 };
        case 1: return {   Z, 255,   0 };
        case 2: return {   0, 255,   Y };
        case 3: return {   0,   Z, 255 };
        case 4: return {   0,   0, 255 };
    }
}, [] (
    const unsigned char red, 
    const unsigned char green, 
    const unsigned char blue
) -> unsigned char {
    int X;
    unsigned char Y;
    if (red > 0) {
        if (red == 255) {
            X = 0;
            Y = green;
        } else {
            X = 1;
            Y = 255 - red;
        }
    } else if (green > 0) {
        if (green == 255) {
            X = 2;
            Y = blue;
        } else {
            X = 3;
            Y = 255 - green;
        }
    } else {
        Y = 255; X = 4;
    }
    return 255 - ceil(64 * (X + (float) Y / 255.0));
});

ColourScheme cyan_black_yellow = ColourScheme([] (
    const unsigned char grey
) -> rgb_t {
    const float d_rb = 3 * (grey - 128);
    const float d_g  = 3 * (std::abs(grey - 128) - 42);
    unsigned char red   = floor(std::min(255.0f, std::max(0.0f,  d_rb)));
    unsigned char green = floor(std::min(255.0f, std::max(0.0f,  d_g )));
    unsigned char blue  = floor(std::min(255.0f, std::max(0.0f, -d_rb)));
    return {red, green, blue};
}, [] (
    const unsigned char red, 
    const unsigned char green, 
    const unsigned char blue
) -> unsigned char {
    return floor(
        red > 0 ?
        (red  < 255 ?  (float) red  / 3.0 + 128 :  (float) green / 3.0 + 42 + 128) :
        (blue < 255 ? -(float) blue / 3.0 + 128 : -(float) green / 3.0 - 42 + 128)
    );
});

ColourScheme greyscale = ColourScheme([] (
    const unsigned char grey
) -> rgb_t {
    return {grey, grey, grey};
}, nullptr);
