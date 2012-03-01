/* ---------------------------------------------------------------------------
 * The ImageProcessing namespace is a container for a number 
 * of methods that modify images in the form of ColorMatrix objects.
 * 
 * CS559 - Project 1 - Practice
 * Brian Ploeckelman
 * ---------------------------------------------------------------------------
 */
#include "ImageProcessing.h"
#include "Matrix2d.h"
#include "Color.h"
#include "tgalibs/TargaImage.h"
#include "tgalibs/libtarga.h"

#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <math.h>
// Frowny face at Microsoft for not including constants from math.h
// without having to specify _USE_MATH_DEFINES... jerks...
#define M_E        2.71828182845904523536
#define M_PI       3.14159265358979323846

using std::stringstream;
using std::string;
using std::vector;
using std::cout;
using std::endl;


/* loadTargaColorMatrix()
 * ----------------------
 * Load the specified targa file into a ColorMatrix object
 * @param targaFilename - name of the targa file to open (including .tga)
 * @return - pixel data from the targa file in a ColorMatrix object
 */
ColorMatrix ImageProcessing::loadTargaColorMatrix(const string& targaFilename)
{
	// Try to load the targa image
	TargaImage* image = TargaImage::readImage(targaFilename.c_str());
	if( image == nullptr ) { 
		stringstream ss;
		ss << "Error: unable to open " << targaFilename << endl;
		throw FailedToLoad(ss.str());
	} 

	// Setup constants from the input image
	const unsigned int   rows   = image->height();
	const unsigned int   cols   = image->width();
	const unsigned char* pixels = image->pixels();
	if( pixels == nullptr ) { 
		stringstream ss;
		ss << "Error: invalid pixel data in " << targaFilename << endl;
		throw BadPixelData(ss.str());	
	}

	// Generate the color matrix based on the image pixels
	ColorMatrix colors(rows, cols);
	for(unsigned int row = 0, i = 0; row < rows; ++row)
	for(unsigned int col = 0; col < cols; ++col, ++i)
	{
		colors(row,col) = Color(pixels[i*4 + 0],	// red
								pixels[i*4 + 1],	// green
								pixels[i*4 + 2],	// blue
								pixels[i*4 + 3]);	// alpha
	}

	// Cleanup the TargaImage
	delete image;

	return colors;
}

/* 
/* saveTargaColorMatrix()
 * ----------------------
 * Save the specified ColorMatrix as the specified targa file 
 * @param targaFilename - name of the targa file to save as (including .tga)
 * @return - true on successful save, false on failure 
 */
bool ImageProcessing::saveTargaColorMatrix(const string& targaFilename, 
										   const ColorMatrix& source) 
{
	// Try to load a blank targa for the output
	TargaImage* dest = TargaImage::blankImage(source.cols(), source.rows());
	if( dest == nullptr ) {
		cout << "Error: unable to create blank image for output." << endl;
		return false;
	}

	// Setup constants for the output image
	const unsigned int rows   = dest->height();
	const unsigned int cols   = dest->width();
	unsigned char*	   pixels = dest->pixels();

	// Write the source ColorMatrix to the dest TargaImage
	for(unsigned int row = 0, i = 0; row < rows; ++row)
	for(unsigned int col = 0; col < cols; ++col, ++i)
	{
		const Color color(source(row,col));
		pixels[i*4 + 0] = color.r();
		pixels[i*4 + 1] = color.g();
		pixels[i*4 + 2] = color.b();
		pixels[i*4 + 3] = color.a();
	}

	// Save the new dest TargaImage
	bool success = true;
	if( dest->write(targaFilename.c_str()) != 1 ) {
		cout << "Error: unable to save " << targaFilename 
			 << ", " << tga_error_string(tga_get_last_error()) 
			 << endl;
		success = false;
	} 

	// Cleanup the new TargaImage
	delete dest;
	
	return success;
}


/* invert()
 * --------
 * Invert the source and store result in dest 
 * @param source - the source pixels
 * @param dest - the destination pixels 
 */
void ImageProcessing::invert(const ColorMatrix& source, ColorMatrix& dest)
{
	for(unsigned int row = 0; row < source.rows(); ++row)
	for(unsigned int col = 0; col < source.cols(); ++col)
	{
		const Color color(source(row, col));
		dest(row, col) = Color( 255 - color.r(), 
								255 - color.g(), 
								255 - color.b(), color.a() );
	}
}

/* grayscale()
 * -----------
 * Grayscale the source and store result in dest 
 * @param source - the source pixels
 * @param dest - the destination pixels
 */
void ImageProcessing::grayscale(const ColorMatrix& source, ColorMatrix& dest)
{
	for(unsigned int row = 0; row < source.rows(); ++row)
	for(unsigned int col = 0; col < source.cols(); ++col)
	{
		const Color		color(source(row,col));
		const Component	avg((color.r() + color.g() + color.b()) / 3) ;
		dest(row, col) = Color(avg, avg, avg, color.a());
	}
}

/* threshold()
 * -----------
 * Split the image into black and white based on the specified threshold 
 * @param source - the source pixels
 * @param dest - the destination pixels
 * @param threshold - the threshold value between black and white
 */
void ImageProcessing::threshold(const ColorMatrix& source, ColorMatrix& dest, 
								const Component& threshold)
{
	ImageProcessing::grayscale(source, dest); 
	for(unsigned int row = 0; row < source.rows(); ++row)
	for(unsigned int col = 0; col < source.cols(); ++col)
		dest(row,col) = source(row,col).r() < threshold ? White() : Black();
}

/* desaturate()
 * ------------
 * Desaturate the source by amount [0,1] and store result in dest 
 * @param source - the source pixels
 * @param dest - the destination pixels
 * @param amount - the amount to desaturate [0,1]
 */
void ImageProcessing::desaturate(const ColorMatrix& source, ColorMatrix& dest, 
								 double amount)
{
	// Clamp desaturation amount to [0,1]
	amount = (amount > 1.0) ? 1.0 : (amount < 0.0 ? 0.0 : amount);

	// Weights for calculating grey values
	const double redScale   = 2.0 / 7.0; 
	const double greenScale = 4.0 / 7.0; 
	const double blueScale  = 1.0 / 7.0; 

	// For each source pixel...
	for(unsigned int row = 0; row < source.rows(); ++row)
	for(unsigned int col = 0; col < source.cols(); ++col)
	{
		// This pixel's source color
		const Color	color(source(row, col));

		// This pixel's weighted grey value
		const Component L(clampf(redScale   * color.r()) 
						+ clampf(greenScale * color.g()) 
						+ clampf(blueScale  * color.b()) );

		// The desaturated pixel's color is calculated by lerping between 
		// source color and gray value based on the desaturation value
		const Component r( clampf(color.r() + (L-color.r()) * amount) ); 
		const Component g( clampf(color.g() + (L-color.g()) * amount) ); 
		const Component b( clampf(color.b() + (L-color.b()) * amount) ); 
		Color desaturatedColor(r, g, b, color.a());

		// Save desaturated pixel out to the destination
		dest(row, col) = desaturatedColor;
	}
}

/* convolve()
 * ----------
 * Convolve an image with a discrete filter 
 * @param source - the source pixels
 * @param dest - the destination pixels to store the output in
 * @param kernel - the filter kernel
 * @param offset - the amount to offset the filtered pixel value
 */
/* Notes:
 *   currently just uses 0 values outside image boundaries
 *     TODO - allow user to choose between zero, clamp, reflect, renormalize
 */
void ImageProcessing::convolve(const ColorMatrix& source, 
							   ColorMatrix& dest, 
							   const KernelMatrix& kernel, 
							   const unsigned int offset)
{
	dest = ColorMatrix(source.rows(), source.cols());
	const int radius = kernel.rows() / 2;

	// For each pixel in the source image...
	for(unsigned int row = 0; row < source.rows(); ++row)
	for(unsigned int col = 0; col < source.cols(); ++col)
	{
		// Setup component accumulators
		double r, g, b;
		r = g = b = 0.0;

		// For each value in the kernel... 
		int kr = -radius;
		for(unsigned int krow = 0; krow < kernel.rows(); ++krow)
		{
			int kc = -radius;
			for(unsigned int kcol = 0; kcol < kernel.cols(); ++kcol)
			{
				// Calculate,validate source pixel indices for this kernel cell
				const int srow = row + kr;	// NOTE: changing + to - switches 
				const int scol = col + kc;  // between convolution,correlation
				if( srow < 0 || srow >= static_cast<int>(source.rows())
				 || scol < 0 || scol >= static_cast<int>(source.cols()) ) {
					++kc;
					continue;
				}
	
				// Get source pixel, kernel value and accumulate components
				const double k = kernel(krow, kcol);
				const Color  c = source(srow, scol);
				r += k * c.r();
				g += k * c.g();
				b += k * c.b(); 
				++kc;
			}
			++kr;
		}

		// Save the resulting pixel
		dest(row, col) = Color(clampf(r+offset), 
							   clampf(g+offset), 
							   clampf(b+offset), 255);
	}
}

/* sep_convolve()
 * --------------
 * Convolve an image with a 1d separable discrete filter 
 * @param source - the source pixels
 * @param dest - the destination pixels to store the output in
 * @param kernel - the filter kernel
 * @param offset - the amount to offset the filtered pixel value
 */
void ImageProcessing::sep_convolve(const ColorMatrix& source, 
								   ColorMatrix& dest, 
								   const std::vector<double>& kernel, 
								   const unsigned int offset)
{
	dest = ColorMatrix(source.rows(), source.cols());
	ColorMatrix intermediate(source.rows(), source.cols());

	const int radius = kernel.size() / 2;

	// For each column 
	for(unsigned int col = 0; col < source.cols(); ++col) {
		vector<Color> colColors = source.getCol(col);

		// For each pixel in the column
		for(unsigned int p = 0; p < colColors.size(); ++p) {
			// Component accums
			double r,g,b;
			r = g = b = 0.0;
	
			// For each kernel value
			int ki = -radius;
			for(unsigned int i = 0; i < kernel.size(); ++i) {
				// Calculate colColors index si + make sure its in bounds
				const int si = p + ki;
				if( si < 0 || si >= static_cast<int>(colColors.size()) ) {
					++ki;
					continue;
				}
				
				const double k = kernel[i];
				const Color  c = colColors[si];
				r += k * c.r();
				g += k * c.g();
				b += k * c.b();
	
				++ki;
			}
			colColors[p] = Color(clampf(r+offset),
								 clampf(g+offset), 
								 clampf(b+offset), 255);
		}
		intermediate.setCol(col, colColors);
	}

	// For each row 
	for(unsigned int row = 0; row < source.rows(); ++row) {
		vector<Color> rowColors = intermediate.getRow(row);

		// For each pixel in the row
		for(unsigned int p = 0; p < rowColors.size(); ++p) {
			// Component accums
			double r,g,b;
			r = g = b = 0.0;

			// For each kernel value
			int ki = -radius;
			for(unsigned int i = 0; i < kernel.size(); ++i)	{
				// Calculate rowColors index si + make sure its in bounds
				int si = p + ki;
				if( si < 0 || si >= static_cast<int>(rowColors.size()) ) {
					++ki;
					continue;
				}
				
				const double k = kernel[i];
				const Color  c = rowColors[si];
				r += k * c.r();
				g += k * c.g();
				b += k * c.b();
	
				++ki;
			}
			rowColors[p] = Color(clampf(r+offset), 
								 clampf(g+offset), 
								 clampf(b+offset), 255);
		}
		dest.setRow(row, rowColors);
	}
}

/* drawSquare()
 * ------------
 * Draws a colored square at centerX,centerY on the given ColorMatrix
 * @param dest - the source (and destination) pixels
 * @param centerX, centerY - the center coordinates of the square
 * @param radius - the radius of the square
 * @param color - the color of the square
 */
void ImageProcessing::drawSquare(ColorMatrix& dest,
								 unsigned int centerX, unsigned int centerY, 
								 unsigned int radius,  const Color& color)
{
	const unsigned int x = centerX;
	const unsigned int y = centerY;
	const unsigned int r = radius;

	for(unsigned int row = (y - r); row <= (y + r); ++row)
	for(unsigned int col = (x - r); col <= (x + r); ++col)
		dest(row,col) = color;
}


/* invert()
 * --------
 * Invert the source and return the result 
 * @param source - the source pixels
 * @return - the result of the operation 
 */
ColorMatrix ImageProcessing::invert(const ColorMatrix& source)
{
	ColorMatrix result(source.rows(), source.cols());
	ImageProcessing::invert(source, result);
	return result;
}

/* grayscale()
 * -----------
 * Grayscale the source and return the result 
 * @param source - the source pixels
 * @return - the result of the operation
 */
ColorMatrix ImageProcessing::grayscale(const ColorMatrix& source)
{
	ColorMatrix result(source.rows(), source.cols());
	ImageProcessing::grayscale(source, result);
	return result;
}

/* threshold()
 * -----------
 * Split the image into black and white based on the specified threshold 
 * @param source - the source pixels
 * @param threshold - the threshold value between black and white
 * @return - the result of the operation
 */
ColorMatrix ImageProcessing::threshold(const ColorMatrix& source, 
									   const Component& threshold)
{
	ColorMatrix result(source.rows(), source.cols());
	ImageProcessing::threshold(source, result, threshold);
	return result;
}

/* desaturate()
 * ------------
 * Desaturate the source by amount [0,1] and returns the result
 * @param source - the source pixels
 * @param amount - the amount to desaturate [0,1]
 * @return - the result of the operation
 */
ColorMatrix ImageProcessing::desaturate(const ColorMatrix& source, double amount)
{
	ColorMatrix result(source.rows(), source.cols());
	ImageProcessing::desaturate(source, result, amount);
	return result;
}
  
/* convolve()
 * ----------
 * Convolve an image with a discrete filter and return the result 
 * @param source - the source pixels
 * @param kernel - the filter kernel
 * @param offset - the amount to offset the filtered pixel value
 * @return - the result of the operation
 */
ColorMatrix ImageProcessing::convolve(const ColorMatrix& source, 
									  const KernelMatrix& kernel, 
									  const unsigned int offset)
{
	ColorMatrix result(1,1); 
	convolve(source, result, kernel, offset);
	return result;
}

/* sep_convolve()
 * --------------
 * Convolve an image with a 1d separable filter and return the result 
 * @param source - the source pixels
 * @param kernel - the filter kernel
 * @param offset - the amount to offset the filtered pixel value
 * @return - the result of the operation
 */
ColorMatrix ImageProcessing::sep_convolve(const ColorMatrix& source, 
										  const std::vector<double>& kernel, 
										  const unsigned int offset)
{
	ColorMatrix result(1,1);
	sep_convolve(source, result, kernel, offset);
	return result;
}

/* sobel()
 * -------
 * Performs sobel edge detection on the source image
 * @param source - the source pixels
 * @return - the result of the operation
 */
ColorMatrix ImageProcessing::sobel(const ColorMatrix& source)
{
	ColorMatrix result(source.rows(), source.cols());
	ColorMatrix blurred(sep_convolve(source, generate1dGaussian(2))); 

	// Detect vertical edges
	KernelMatrix sobelVertical(3,3);
	sobelVertical(0,0) = -1.0;
	sobelVertical(1,0) = -2.0;
	sobelVertical(2,0) = -1.0;
	sobelVertical(0,2) =  1.0;
	sobelVertical(1,2) =  2.0;
	sobelVertical(2,2) =  1.0;

	ColorMatrix verticalEdges(convolve(blurred, sobelVertical));

	// Detect horizontal edges
	KernelMatrix sobelHorizontal(3,3);
	sobelHorizontal(0,0) =  1.0;
	sobelHorizontal(0,1) =  2.0;
	sobelHorizontal(0,2) =  1.0;
	sobelHorizontal(2,0) = -1.0;
	sobelHorizontal(2,1) = -2.0;
	sobelHorizontal(2,2) = -1.0;

	ColorMatrix horizontalEdges(convolve(blurred, sobelHorizontal));

	// Save gradient values
	for(unsigned int row = 0; row < source.rows(); ++row)
	for(unsigned int col = 0; col < source.cols(); ++col)
	{
		// Get vertical and horizontal edge values
		const Color v(verticalEdges(row,col));
		const Color h(horizontalEdges(row,col));

		// Calculate the gradient value sqrt(Gx^2 + Gy^2)
		const Component r = clampf(fabs(sqrt(static_cast<double>(v.r()*v.r()) + (h.r()*h.r()))));
		const Component g = clampf(fabs(sqrt(static_cast<double>(v.g()*v.g()) + (h.g()*h.g()))));
		const Component b = clampf(fabs(sqrt(static_cast<double>(v.b()*v.b()) + (h.b()*h.b()))));
		// TODO : add support for quick gradient calc: |G| = |Gx| + |Gy|
		
		result(row,col) = Color(r, g, b, source(row,col).a());
	}

	return result;
}

/* resize()
 * --------
 * TODO - modify this function so that it takes an enumeration 
 * that specifies which type of resampling kernel to use, then 
 * calls the appropriate resize method... 
 * (unless we can come up with a nice way to make the resizer more generic)
 * @return - the result of the operation
 */
ColorMatrix ImageProcessing::resize(const ColorMatrix& source, 
									const unsigned int rows, 
									const unsigned int cols)
{
	ColorMatrix result(rows, cols);

	// TODO - see comment above...

	return result;
}

/* resizeNearest()
 * ---------------
 * Resizes the image using nearest-neighbor sampling
 * @param source - the source pixels
 * @param rows - the number of rows (height) of the resized image
 * @param cols - the number of columns (width) of the resized image
 * @return - the result of the operation
 */
ColorMatrix ImageProcessing::resizeNearest(const ColorMatrix& source, 
										   const unsigned int rows, 
										   const unsigned int cols)
{
	ColorMatrix result(rows, cols);

	const double rowScale = static_cast<double>(source.rows()) / rows;
	const double colScale = static_cast<double>(source.cols()) / cols;

	for(unsigned int row = 0; row < rows; ++row)
	for(unsigned int col = 0; col < cols; ++col)
	{
		const unsigned int nrow = static_cast<unsigned int>(floor(row * rowScale)); 
		const unsigned int ncol = static_cast<unsigned int>(floor(col * colScale)); 

		// Nearest neighbor sample
		result(row, col) = source(nrow, ncol);
	}

	return result;
}

/* resizeBilinear()
 * ---------------
 * TODO - handle edges differently, currently just clamps
 * Resizes the image using bilinear sampling
 * @param source - the source pixels
 * @param rows - the number of rows (height) of the resized image
 * @param cols - the number of columns (width) of the resized image
 * @return - the result of the operation
 */
ColorMatrix ImageProcessing::resizeBilinear(const ColorMatrix& source, 
										    const unsigned int rows, 
										    const unsigned int cols)
{
	ColorMatrix result(rows, cols);

	const double rowScale = static_cast<double>(source.rows()) / rows;
	const double colScale = static_cast<double>(source.cols()) / cols;

	for(unsigned int row = 0; row < rows; ++row)
	for(unsigned int col = 0; col < cols; ++col)
	{
		// Calculate new sample location
		const double y = row * rowScale;
		const double x = col * colScale;

		// Calculate nearest pixel (above and left of sample location)
		unsigned int nrow = static_cast<unsigned int>(floor(y)); 
		unsigned int ncol = static_cast<unsigned int>(floor(x));

		// Calculate interpolation factors for row and column 
		const double rowLerpFactor = y - nrow;
		const double colLerpFactor = x - ncol;

		// Setup neighbor pixel values
		Color topLeft, topRight, botLeft, botRight;

		topLeft  = source(nrow, ncol);
		// Make sure we don't go out of bounds with the +1's 
		if( ncol + 1 < source.cols() )
			topRight = source(nrow, ncol + 1);
		if( nrow + 1 < source.rows() )
			botLeft  = source(nrow + 1, ncol);
		if( nrow + 1 < source.rows() && ncol + 1 < source.cols() )
			botRight = source(nrow + 1, ncol + 1);

		// Interpolate between left and right pixels by lerp factors
		Color topRowLerp( lerp(colLerpFactor, topLeft, topRight) );
		Color botRowLerp( lerp(colLerpFactor, botLeft, botRight) );

		//Interpolate between top and bottom lerp'd pixels				
		Color interpolated( lerp(rowLerpFactor, topRowLerp, botRowLerp) );

		// Save bilinear sample 
		result(row, col) = interpolated; 
	}

	return result;
}

/* generateBlurKernel()
 * --------------------
 * TODO : refactor this once KernelMatrix is changed to a 1st-class type
 * Generates a blurring kernel for convolution based on user-specified params
 * @param radius - the radius of the kernel (1 => 3x3, 2 => 5x5, etc...)
 * @param type - the type of kernel (-1 => box, 0 => gaussian, etc...)
 * @return - the generated KernelMatrix object
 */
KernelMatrix ImageProcessing::generateBlurKernel(const int radius, 
												 const KernelType& type)
{
	// Validate params
	if( radius < 0) 
		throw KernelMatrix::BadSize("Error: bad radius for blurring kernel.");

	// Create the KernelMatrix of the appropriate size
	const unsigned int dim = 2 * radius + 1;
	KernelMatrix kernel(dim, dim);

	// Generate values for the KernelMatrix
	switch(type)
	{
	case gaussian_2d:	
		for(int i = -radius; i <= radius; ++i)
		for(int j = -radius; j <= radius; ++j)
		{
			const unsigned int row = i + radius;
			const unsigned int col = j + radius;
			kernel(row, col) = gaussianFunction(i, j);
		}
		normalize(kernel);
		break;
	case box: 		
		kernel.fill(1.0 / (dim*dim));	
		break;
	default:		
		throw UnsupportedKernelType("Error: generateBlurKernel() - unsupported kernel type encountered.");	
		break;
	};

	return kernel;
}

/* generate1dGaussian()
 * --------------------
 * Generates a 1-dimensional filter kernel 
 * to be used for separable blurring convolution operation.
 * @param radius - the radius of the filter 
 * @return - the calculated 1d filtering kernel
 */
vector<double> ImageProcessing::generate1dGaussian(const int radius)
{
	vector<double> kernel;
	for(int x = -radius; x <= radius; ++x)
		kernel.push_back(gaussianFunction(x));
	normalize(kernel);
	return kernel; 
}

/* generate1dBSpline()
 * --------------------
 * Generates a 1-dimensional bspline filter kernel 
 * to be used for separable blurring convolution operation.
 * @param radius - the radius of the filter 
 * @return - the calculated 1d filtering kernel
 */
vector<double> ImageProcessing::generate1dBSpline(const int radius)
{
	vector<double> kernel;
	for(int x = -radius; x <= radius; ++x)
		kernel.push_back(bsplineFunction(x));
	normalize(kernel);
	return kernel;
}

/* generate1dCatmull()
 * -------------------
 * Generates a 1-dimensional filter kernel 
 * to be used for separable blurring convolution operation.
 * @param radius - the radius of the filter 
 * @return - the calculated 1d filtering kernel
 */
vector<double> ImageProcessing::generate1dCatmull(const int radius)
{
	vector<double> kernel;
	for(int x = -radius; x <= radius; ++x)
		kernel.push_back(catmullFunction(x));
	normalize(kernel);
	return kernel;
}

/* generate1dMitchell()
 * --------------------
 * Generates a 1-dimensional mitchell-netravali filter kernel 
 * to be used for separable blurring convolution operation.
 * @param radius - the radius of the filter 
 * @return - the calculated 1d filtering kernel
 */
vector<double> ImageProcessing::generate1dMitchell(const int radius, 
												   const double B, 
												   const double C)
{
	vector<double> kernel;
	for(int x = -radius; x <= radius; ++x)
		kernel.push_back(mitchellFunction(x, B, C));
	normalize(kernel);
	return kernel;
}

/* gaussianFunction()
 * ------------------
 * Calculates and returns the gaussian value for the given 2d position
 * @param i - the "row"
 * @param j - the "column"
 * @param sigma - the scale factor for this gaussian
 * @return - the gaussian value 
 */
double ImageProcessing::gaussianFunction(const int i, const int j, 
										 const double sigma)
{
	const double twoSigmaSquared = 2.0 * sigma * sigma;
	const double scale =  1.0 / (M_PI * twoSigmaSquared); 
	const double expon = -1.0 * ((i*i + j*j) / twoSigmaSquared); 
	const double e = pow(M_E, expon);
	return scale * e;
}

/* gaussianFunction()
 * ------------------
 * Calculates and returns the gaussian value for the given 1d position
 * @param x - the position 
 * @param sigma - the scale factor for this gaussian
 * @return - the gaussian value 
 */
double ImageProcessing::gaussianFunction(const double x, const double sigma)
{
	const double twoSigmaSquared = 2.0 * sigma * sigma;
	const double scale =  1.0 / (M_PI * twoSigmaSquared);
	const double expon = -1.0 * ((x*x) / twoSigmaSquared);
	const double e = pow(M_E, expon);
	return scale * e;
}

/* bsplineFunction()
 * -----------------
 * Calculates the b-spline value for the given parameter
 * @param x - the parameter 
 * @return - the calculated b-spline value 
 */ 
double ImageProcessing::bsplineFunction(const double x) 
{
	const double ax = fabs(x);
	const double one_ax = 1.0 - ax;
	const double two_ax = 2.0 - ax;

	if( x >= -1.0 && x <= 1.0 )
		return (-3.0 * one_ax * one_ax * one_ax 
			  +  3.0 * one_ax * one_ax
			  +  3.0 * one_ax 
			  +  1.0);
	else if( ax >= 1.0 && ax <= 2.0 ) 
		return two_ax * two_ax * two_ax;

	return 0.0;
}

/* catmullFunction()
 * -----------------
 * Calculates the catmull-rom spline value for the given parameter
 * @param x - the parameter 
 * @return - the calculated catmull-rom spline value 
 */ 
double ImageProcessing::catmullFunction(const double x)
{
	const double ax = fabs(x);
	const double one_ax = 1.0 - ax;
	const double two_ax = 2.0 - ax;

	if( x >= -1 && x <= 1 )
		return (-3.0 * one_ax * one_ax * one_ax
		      +  4.0 * one_ax * one_ax 
			  +  1.0 * one_ax) * 0.5;
	else if( ax >= 1.0 && ax <= 2.0 ) 
		return ((two_ax * two_ax * two_ax) 
			  - (two_ax * two_ax)) * 0.5;
	
	return 0.0;
}

/* mitchellNetravali()
 * -------------------
 * TODO : complete mitchellNetravali documentation 
 * Calculates the Mitchell-Netravali filter value for the given parameters
 * @param x - the parameter 
 * @param B - 
 * @param C - 
 * @return - the calculated Mitchell-Netravali value 
 */ 
double ImageProcessing::mitchellFunction(const double x, 
										 const double B, 
										 const double C)
{
	const double ax = fabs(x);
	const double one_ax = 1.0 - ax;
	const double two_ax = 2.0 - ax;

	if( x >= -1.0 && x <= 1.0 )
		return (-15.0 * one_ax * one_ax * one_ax
		      +  18.0 * one_ax * one_ax
			  +   9.0 * one_ax
			  +   2.0) * (1.0 / 18.0);
	else if( ax >= 1.0 && ax <= 2.0 )
		return ((5.0 * two_ax * two_ax * two_ax)
		      - (3.0 * two_ax * two_ax)) * (1.0 / 18.0);

	return 0.0;
}

/* reconstruct()
 * -------------
 * Reconstructs a Color value from a sequence of Color values, 
 * by using the specified filter kernel (Mitchell-Netravali at the moment)
 * @param pixels - the sequence of Color values 
 * @param filter - the filter kernel
 * @param x - the floating point filter variable
 * @return the reconstructed pixel color
 /
Color ImageProcessing::reconstruct(const std::vector<Color>& pixels, 
								   const std::vector<double>& filter, 
								   double x)
{
	const unsigned int r = filter.size() / 2;
	Color color; 

	unsigned int start = static_cast<unsigned int>(ceil(x - r));
	unsigned int end   = static_cast<unsigned int>(floor(x + r));
	for(unsigned int i = start; i >= end; ++i) {
		const double f = mitchellNetravali(x - i);
		color.r( clampf(color.r() + pixels[i].r() * f) )
			 .g( clampf(color.g() + pixels[i].g() * f) )
			 .b( clampf(color.b() + pixels[i].b() * f) )
			 .a( 255 );
	}
	
	return color;
}
*/

/* normalize()
 * -----------
 * Normalizes the specified 2d KernelMatrix in place
 * @param kernel - the kernel to normalize
 */
void ImageProcessing::normalize(KernelMatrix& kernel)
{
	double weight = 0.0;

	for(unsigned int row = 0; row < kernel.rows(); ++row)
	for(unsigned int col = 0; col < kernel.cols(); ++col)
		weight += kernel(row,col);

	for(unsigned int row = 0; row < kernel.rows(); ++row)
	for(unsigned int col = 0; col < kernel.cols(); ++col)
		kernel(row, col) /= weight;
}

/* normalize()
 * -----------
 * Normalizes the specified 1d kernel in place
 * @param kernel - the kernel to normalize
 */
void ImageProcessing::normalize(vector<double>& kernel)
{
	double weight = 0.0;
	for(unsigned int i = 0; i < kernel.size(); ++i)	weight += kernel[i];
	for(unsigned int i = 0; i < kernel.size(); ++i) kernel[i] /= weight;
}

/* clampf()
 * --------
 * Clamp to the range [min, max](default [0, 255]), returns color Component 
 * @param value - the value to clamp
 * @param min - the minimum value to clamp to (default = 0)
 * @param max - the maximum value to clamp to (default = 255)
 * @return - the clamped Component value
 */ 
inline Component ImageProcessing::clampf(double value, 
										 const Component& min, 
										 const Component& max)
{
	return Component(value > max ? max : (value < min ? min : value));
}

/* lerp()
 * ------
 * Lerp between Color values c1 and c2 by amount t [0,1] 
 * @param t - the interpolation parameter
 * @param c1 - first color value
 * @param c2 - second color value
 * @return - the interpolated color value
 */
Color ImageProcessing::lerp(double t, const Color& c1, const Color& c2)
{
	// Clamp t to the range [0,1]
	t = (t < 0.0) ? 0.0 : ((t > 1.0) ? 1.0 : t);

	// Return the interpolated color
	return Color(
		clampf(c1.r() + t * (c2.r() - c1.r())),
		clampf(c1.g() + t * (c2.g() - c1.g())), 
		clampf(c1.b() + t * (c2.b() - c1.b())), 
		clampf(c1.a() + t * (c2.a() - c1.a()))
	);
}