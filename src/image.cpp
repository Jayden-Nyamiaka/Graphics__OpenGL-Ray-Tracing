#include "image.h"

#include <png.h>

const int MAX_INTENSITY = 255;

Image::Image() {
    Image(100, 100);
}

Image::Image(int width, int height) {
    xres = width;
    yres = height;
    pixels = std::vector<Eigen::Vector3f>(width * height);
    for (int i = 0; i < width * height; i++) {
        pixels[i] = Eigen::Vector3f::Ones();
    }
}

void Image::SetPixel(int x, int y, const Eigen::Vector3f &color) {
    pixels[x + xres * y] = color;
}

bool Image::SaveImage(const std::string &path) {
    FILE *fp;
    png_structp png_ptr = NULL;
    png_infop info_ptr = NULL;
    png_byte **row_pointers = NULL;

    fp = fopen(path.c_str(), "wb");
    if (!fp) {
        return false;
    }

    // Create the data and info structures used by libpng
    png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
    info_ptr = png_create_info_struct(png_ptr);
    if (!png_ptr || !info_ptr)
        return false;

    // Set up libpng to write a XRES x YRES image with 8 bits of RGB color depth
    png_set_IHDR(png_ptr, info_ptr, xres, yres, 8,
        PNG_COLOR_TYPE_RGB, PNG_INTERLACE_NONE, PNG_COMPRESSION_TYPE_DEFAULT,
        PNG_FILTER_TYPE_DEFAULT);
    
    // Allocate an array of pointers to each row of pixels in the PNG file
    row_pointers = (png_byte **) png_malloc(png_ptr, yres * sizeof(png_byte *));
    for (int y = 0; y < yres; y++) {
        // Allocate a row of pixels in the PNG file
        png_byte *row = (png_byte *) png_malloc(png_ptr, 3 * xres * sizeof(uint8_t));
        // The screen matrix has (x, y) = (0, 0) in the lower left so that it's
        // easy to translate between it and the "camera sensor" grid, so we have
        // to go through the rows in reverse because of convention
        row_pointers[yres - 1 - y] = row;
        for (int x = 0; x < xres; x++) {
            // Get the pixel's location and fill in its red, green, and blue
            // values, bounding them between 0 and MAX_INTENSITY = 255
            Eigen::Vector3f &pixel = pixels[x + xres * y];
            *row++ = (uint8_t) (fmin(fmax(0.0, pixel.x()), 1.0) * MAX_INTENSITY);
            *row++ = (uint8_t) (fmin(fmax(0.0, pixel.y()), 1.0) * MAX_INTENSITY);
            *row++ = (uint8_t) (fmin(fmax(0.0, pixel.z()), 1.0) * MAX_INTENSITY);
        }
    }

    // Write the PNG to file
    png_init_io(png_ptr, fp);
    png_set_rows(png_ptr, info_ptr, row_pointers);
    png_write_png(png_ptr, info_ptr, PNG_TRANSFORM_IDENTITY, NULL);
    
    // Free each row of the PNG that we allocated
    for (int y = 0; y < yres; y++) {
        png_free(png_ptr, row_pointers[y]);
    }
    
    // Free the array of row pointers
    png_free(png_ptr, row_pointers);
    fclose(fp);

    return true;
}