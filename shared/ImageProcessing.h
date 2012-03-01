#pragma once
/* ---------------------------------------------------------------------------
 * The ImageProcessing namespace is a container for a number 
 * of methods that modify images in the form of ColorMatrix objects.
 * 
 * CS559 - Project 1 - Practice
 * Brian Ploeckelman
 * ---------------------------------------------------------------------------
 */
#include "Matrix2d.h"
#include "Color.h"

#include <stdexcept>
#include <string>

typedef Matrix2d<Color>  ColorMatrix;
typedef Matrix2d<double> KernelMatrix;
// TODO : make KernelMatrix a real class with a Matrix2d<double> member
// this will allow for nicer construction of kernels, 
// from double[] for example...
// double[] k = { 1.0, 2.0, 1.0, 
//				  2.0, 4.0, 2.0, 
//				  1.0, 2.0, 1.0 };
// KernelMatrix kernel(k);
// also, then the offset/bias and scaling factor can be built in to the data type
// and the number of rows/cols can be validated as odd and equal
// and there can be a getRadius(), and 1d and 2d can be abstracted away


namespace ImageProcessing
{
	/* _____________________________________________________________________ */
	/* ------------------------ Enumerations ------------------------------- */

	enum KernelType { box=0, gaussian, bspline, catmull, mitchell, gaussian_2d };

	/* _____________________________________________________________________ */
	/* ------------------------ Exceptions --------------------------------- */

	// loadTargaColorMatrix() throws FailedToLoad, or BadPixelData exceptions
	class FailedToLoad : public std::runtime_error { public: FailedToLoad(const std::string& what_arg) : std::runtime_error(what_arg) { } };
	class BadPixelData : public std::runtime_error { public: BadPixelData(const std::string& what_arg) : std::runtime_error(what_arg) { } };

	// generateBlurKernel() throws UnsupportedKernelType exception 
	class UnsupportedKernelType : public std::runtime_error { public: UnsupportedKernelType(const std::string& what_arg) : std::runtime_error(what_arg) { } };

	/* _____________________________________________________________________ */
	/* ------------------- Targa loading/saving ---------------------------- */ 

	ColorMatrix loadTargaColorMatrix(const std::string& targaFilename);
	bool		saveTargaColorMatrix(const std::string& targaFilename, const ColorMatrix& source);

	/* _____________________________________________________________________ */
	/* ------------ Image manipulation functions --------------------------- */

	// In-place destination ColorMatrix 
	void invert		(const ColorMatrix& source, ColorMatrix& dest);
	void grayscale	(const ColorMatrix& source, ColorMatrix& dest);
	void threshold	(const ColorMatrix& source, ColorMatrix& dest, const Component& threshold=128);
	void desaturate	(const ColorMatrix& source, ColorMatrix& dest, double amount);
	void convolve	(const ColorMatrix& source, ColorMatrix& dest, const KernelMatrix& kernel, const unsigned int offset=0);
	void sep_convolve(const ColorMatrix& source, ColorMatrix& dest, const std::vector<double>& kernel, const unsigned int offset=0);
	void drawSquare	(ColorMatrix& dest, unsigned int centerX, unsigned int centerY, unsigned int radius,  const Color& color);

	// Returning destination ColorMatrix 
	ColorMatrix invert		(const ColorMatrix& source);
	ColorMatrix grayscale	(const ColorMatrix& source);
	ColorMatrix threshold	(const ColorMatrix& source, const Component& threshold=128);
	ColorMatrix desaturate	(const ColorMatrix& source, double amount);
	ColorMatrix convolve	(const ColorMatrix& source, const KernelMatrix& kernel, const unsigned int offset=0);
	ColorMatrix sep_convolve(const ColorMatrix& source, const std::vector<double>& kernel, const unsigned int offset=0);
	ColorMatrix sobel		(const ColorMatrix& source);

	// TODO : work in progress, needs to be a generic resizer that calls out to the more specific resizing functions...
	ColorMatrix resize        (const ColorMatrix& source, const unsigned int rows, const unsigned int cols);
	ColorMatrix resizeNearest (const ColorMatrix& source, const unsigned int rows, const unsigned int cols);
	ColorMatrix resizeBilinear(const ColorMatrix& source, const unsigned int rows, const unsigned int cols);

	/* _____________________________________________________________________ */
	/* ------------------------ Kernel related ----------------------------- */
	KernelMatrix generateBlurKernel(const int radius, const KernelType& type=box);

	std::vector<double> generate1dGaussian(const int radius);
	std::vector<double> generate1dBSpline (const int radius);
	std::vector<double> generate1dCatmull (const int radius);
	std::vector<double> generate1dMitchell(const int radius, const double B=(1.0/3.0), const double C=(1.0/3.0));

	double gaussianFunction(const int i, const int j, const double sigma=1.0);
	double gaussianFunction(const double x, const double sigma=1.0);
	double bsplineFunction (const double x);
	double catmullFunction (const double x);
	double mitchellFunction(const double x, const double B, const double C);

	void normalize(KernelMatrix& kernel);	
	void normalize(std::vector<double>& kernel);

	/* _____________________________________________________________________ */
	/* -------------------- Utility functions ------------------------------ */

	// Clamp a double value into the unsigned char range needed for a Component value
	Component clampf(double value, const Component& min=0, const Component& max=255);
	
	// Lerp between two color values
	Color lerp(double t, const Color& c1, const Color& c2);
};

