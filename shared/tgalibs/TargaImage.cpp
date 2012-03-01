// simple C++ class wrapper for the data in LibTarga - designed to
// make things a little bit easier.
//
// class designed by Mike Gleicher 9/5/2007
// Implementation by Chi Man Liu 9/2007
//
// be sure to give proper attribution for this (and any other code)
// that you use in your projects!
//

#include <malloc.h>
#include <string.h>

#include "TargaImage.h"
#include "libtarga.h"

#define BYTES_PER_PIXEL 4

// creates an empty (black, alpha=0) image
// remember that the constructor will copy the pixels, so we need to free them
TargaImage* TargaImage::blankImage(const unsigned int width, const unsigned int height) {
	unsigned char* pixels = (unsigned char*)tga_create(width, height, TGA_TRUECOLOR_32);
	memset(pixels, 0, width * height * BYTES_PER_PIXEL);
	TargaImage *img = new TargaImage(width, height, pixels);
	free(pixels);
	return img;
}

// reads a targa image from a file- returns a pointer to an image if
// (and only if) the read was successful. returns a NULL pointer if
// there is an error. use the targa error check to get the error
// message 
TargaImage* TargaImage::readImage(const char* filename) {
	int w, h;
	unsigned char* pixels = (unsigned char*)tga_load(filename, &w, &h, TGA_TRUECOLOR_32);
	if (pixels == NULL)
		return NULL;
	TargaImage *img = new TargaImage(w, h, pixels);
	img->flip();
	free(pixels);
	return img;
}

// destuctor - remember to call this on anything you create
TargaImage::~TargaImage() {
	free(_pixels);
}


// flip the image in place (re-orders the pixels)
void TargaImage::flip() {
	unsigned int half_height = _height / 2;
	unsigned int bytes_per_row = _width * BYTES_PER_PIXEL;
	for (unsigned int i = 0; i < half_height; ++i) {
		unsigned int another = _height - i - 1;
		for (unsigned int j = 0; j < bytes_per_row; ++j) {
			unsigned char tmp = _pixels[i * bytes_per_row + j];
			_pixels[i * bytes_per_row + j] = _pixels[another * bytes_per_row + j];
			_pixels[another * bytes_per_row + j] = tmp;
		}
	}
}

// write the image to a file - returns 1 if there is success,
// otherwise check the LibTarga error code
int TargaImage::write(const char* filename) const {
	TargaImage *temp = clone();
	temp->flip();
	int retval = tga_write_raw(filename, _width, _height, temp->pixels(), TGA_TRUECOLOR_32); //, 1);
	return retval;
}

// set all the alphas to a specific value
void TargaImage::fillAlpha(const unsigned char value) {
	for (unsigned int i = 0; i < _height; ++i)
		for (unsigned int j = 0; j < _width; ++j)
			_pixels[(i * _width + j) * BYTES_PER_PIXEL + (BYTES_PER_PIXEL - 1)] = value;
}

// this takes the place of a copy constructor. it creates a "deep"
// copy of the image. the caller is responsible for deleting what
// gets created
TargaImage* TargaImage::clone() const {
	return new TargaImage(_width, _height, _pixels);
}

// private constructor
TargaImage::TargaImage(const unsigned int width, const int height, const unsigned char* pixels) {
	_width = width;
	_height = height;
	_pixels = (unsigned char*)tga_create(width, height, TGA_TRUECOLOR_32);
	memcpy(_pixels, pixels, width * height * BYTES_PER_PIXEL);
}

// this null constructor should never be called - it is here to make
// sure it is private and cannot be called accidentally
TargaImage::TargaImage() {
}