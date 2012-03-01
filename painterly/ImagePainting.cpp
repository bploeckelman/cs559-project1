/* ---------------------------------------------------------------------------
 * The ImagePainting namespace is a container for a number 
 * of methods that perform a painterly-rendering algorithm on an image
 * 
 * CS559 - Project 1 - Painterly
 * Brian Ploeckelman
 * ---------------------------------------------------------------------------
 */
#include "ImagePainting.h"
#include "ImageProcessing.h"

#include <algorithm> // for std::random_shuffle()
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <list>
#include <set>

using namespace ImageProcessing;
using std::stringstream;
using std::string;
using std::vector;
using std::list;
using std::cout;
using std::endl;
using std::set;


/* _______________________________________________ */
/* -------------------- Brush -------------------- */
Brush::Brush(const unsigned int diameter, const Type& type)
	: diameter(diameter), type(type), weights(1, 1)
{ 
	switch(type)
	{
	case Abstract:  weights = loadTargaColorMatrix("brushes/brush-abstract.tga");  break;
	case square:    weights = loadTargaColorMatrix("brushes/brush-square.tga");    break;
	case circle:    weights = loadTargaColorMatrix("brushes/brush-circle.tga");    break;
	case feather:   weights = loadTargaColorMatrix("brushes/brush-feather.tga");   break;
	case backslash: weights = loadTargaColorMatrix("brushes/brush-backslash.tga"); break;
	case foreslash: weights = loadTargaColorMatrix("brushes/brush-foreslash.tga"); break;
	case x:			weights = loadTargaColorMatrix("brushes/brush-x.tga");		   break;
	case squiggle:  
		{ // If this type is specified, pick a random brush-squiggleX.tga image
			int squiggleNum = ((rand() % 7) + 1);
			stringstream ss;
			ss << "brushes/brush-squiggle" << squiggleNum << ".tga";
			weights = loadTargaColorMatrix(ss.str());
			break;
		}
	default: throw UnsupportedBrushType("Error: tried to create a brush of an unsupported type.");
	}
	weights = resizeBilinear(weights, diameter, diameter);
}


/* ________________________________________________ */
/* -------------------- Stroke -------------------- */
Stroke::Stroke(const unsigned int x, 
			   const unsigned int y, 
			   const Brush& brush, 
			   const Color& color)
	: x(x), y(y), brush(brush.size(), brush.getType()), color(color)
{ }


/* apply()
 * -------
 * Applies this stroke to the specified canvas
 * @param canvas - the canvas to draw to
 */
void Stroke::apply(ColorMatrix& canvas)
{
	// The stroke radius
	const int r = brush.size() / 2;

	// Start/end positions of canvas region corresponding to this stroke
	const int rowStart = y - r;
	const int rowEnd   = y + r;
	const int colStart = x - r;
	const int colEnd   = x + r;

	// Brush indices
	unsigned int brow = 0;
	unsigned int bcol = 0;

	// For each pixel in the canvas region corresponding to this stroke
	for(int row = rowStart; row < rowEnd; ++row)
	for(int col = colStart; col < colEnd; ++col)
	{
		// Don't try to write out of the canvas bounds 
		if( row >= 0 && row < static_cast<int>(canvas.rows()) 
		 && col >= 0 && col < static_cast<int>(canvas.cols()) ) { 
			// Scale the stroke color by the brush weight at that pixel
			const double w = brush.weight(brow, bcol);
			Color strokeColor = color;
			strokeColor.r(clampf(w * strokeColor.r()))
					   .g(clampf(w * strokeColor.g()))
					   .b(clampf(w * strokeColor.b()))
					   .a(clampf(w * strokeColor.a()));

			// Calculate the destination color by blending between 
			// the canvas color and the stroke color based on the stroke's alpha
			const double srcA = strokeColor.a() / 255.0;
			const Color canvasColor = canvas(row, col);
			Color destColor;
			destColor.r(clampf(srcA * strokeColor.r() + (1 - srcA) * canvasColor.r()))
					 .g(clampf(srcA * strokeColor.g() + (1 - srcA) * canvasColor.g()))
					 .b(clampf(srcA * strokeColor.b() + (1 - srcA) * canvasColor.b()))
					 .a(255);

			// Apply the color to the canvas
			canvas(row, col) = destColor; 

			// Increment the brush indices 
			++bcol;
			if( bcol >= brush.size() ) { 
				bcol = 0;
				++brow;
			}
		}
	}
}

/* __________________________________________________ */
/* -------------------- Painting -------------------- */
const bool Painting::debug = false;

Painting::Painting(const ColorMatrix& source)
	: source(source), 
	canvas(source.rows(), source.cols()), 
	brushes(), 
	gridSizeFactor(0.25), 
	errorThreshold(100.0)
{ }

/* paint()
 * -------
 * Generate a painterly rendering of this Painting object's source image
 * using the current set of brushes. 
 */ 
void Painting::paint()
{
	cout << "Painting image..." << endl;

	// Create a new canvas with a constant color (opaque black in this case)
	canvas = ColorMatrix(source.rows(), source.cols());
	canvas.fill(Color(255,255,255,255)); 
	// NOTE : the algorithm calls for the canvas to be painted a special color 
	// C such that: "the difference between C and any color is MAXINT"
	// I don't understand how the difference between a constant color 
	// and any other arbitrary color can be a constant value... 

	// For each brush, from largest to smallest
	list<Brush>::iterator it  = brushes.begin();
	list<Brush>::iterator end = brushes.end();
	for(; it != end; ++it) {
		const Brush& brush = *it;

		// Blur the source image with a gaussian 
		// kernel based on the current brush size
		cout << "\tGenerating blurred reference image..." << endl; 
		ColorMatrix reference(1,1);
		sep_convolve(source, reference, generate1dGaussian(brush.size() / 2));

		// Save the reference image if in debug mode 
		if( debug ) {
			static int refNum = 0;
			stringstream ss; 
			ss << "output-reference" << refNum++ << ".tga";
			saveTargaColorMatrix(ss.str(), reference);
		}

		// Paint a new layer based on the reference with the current brush
		paintLayer(reference, brush);		
	}

	cout << "Painting completed." << endl;
}

/* paintLayer()
 * ------------
 * Paints a layer onto the canvas based on the reference image and current brush.
 * This is the workhorse of the painterly algorithm.
 * @param reference - the reference image to use
 * @param brush - the brush to use
 */
void Painting::paintLayer(const ColorMatrix& reference, const Brush& brush) 
{
	cout << "\tPainting new layer: brush size = " << brush.size() << endl;

	// Create an empty set of strokes for this layer
	vector<Stroke> strokes;

	cout << "\t\tCalculating difference image..." << endl;

	// Generate an image of the difference between the canvas and reference image
	ColorMatrix diff(canvas.rows(), canvas.cols());
	difference(canvas, reference, diff);	

	// Save difference images if in debug mode
	if( debug ) {
		static int diffNum = 0; 
		stringstream ss;
		ss << "output-difference" << diffNum++ << ".tga";
		saveTargaColorMatrix(ss.str(), diff);
	}

	// Set the grid step size based on the brush size and a constant factor
	const unsigned int grid = static_cast<unsigned int>(ceil(gridSizeFactor * brush.size()));

	// Generate a set of strokes for this layer
	for(unsigned int y = 0; y < canvas.rows(); y += grid)
	for(unsigned int x = 0; x < canvas.cols(); x += grid)
	{
		// Get a copy of the local region from the difference image
		ColorMatrix D(grid, grid);
		copyRegion(diff, x, y, grid, D);

/* Note: this is only for serious debugging... it creates LOTS of extra images
		if( debug ) {
			static unsigned int dNum = 0;
			stringstream ss;
			ss << "output-D" << dNum++ << ".tga";
			saveTargaColorMatrix(ss.str(), D);
		}
*/
		// Calculate the overall error for this region
		const double areaError = sumDifferences(D) / (grid * grid);

		// If the error for this region is above the threshold: 
		// Find the indices for the largest error value in the image, 
		// and make a stroke with the current brush at that location 
		// from the reference image.
		if( areaError > errorThreshold ) { 
			// Find indices of largest local error in difference region
			const unsigned int dx = maxErrorColIndex(D);
			const unsigned int dy = maxErrorRowIndex(D);
			// Find corresponding indices for reference image
			const unsigned int x1 = x - (grid / 2) + dx;
			const unsigned int y1 = y - (grid / 2) + dy;
			// Create a new stroke at the specified location 
			strokes.push_back( makeStroke(brush, x1, y1, reference) );
		}
	}

	// Randomize the Stroke vector
	std::random_shuffle(strokes.begin(), strokes.end());
	
	cout << "\t\tRendering " << strokes.size() << " strokes..." << endl;

	// Paint all the strokes in a random order
	vector<Stroke>::iterator it  = strokes.begin();
	vector<Stroke>::iterator end = strokes.end();
	for(; it != end; ++it) 
		(*it).apply(canvas);

	// Save intermediate canvas states if in debug mode
	if( debug )	{
		static int intermedNum = 0;
		stringstream ss;
		ss << "canvas-intermediate" << intermedNum++ << ".tga";
		saveTargaColorMatrix(ss.str(), canvas);
	}

	cout << "\tLayer painted." << endl;
}

/* difference()
 * ------------
 * Calculates the difference image between two source images
 * uses |rgb1 - rgb2| = sqrt((r1-r2)^2 + (g1-g2)^2 + (b1-b2)^2))
 * @param c1 - the first source image
 * @param c2 - the second source image
 * @param diff - the image to store the difference values to
 */ 
void Painting::difference(const ColorMatrix& c1, const ColorMatrix& c2, ColorMatrix& diff)
{
	// Validate that the sizes are the same
	if( c1.rows() != c2.rows() || c1.cols() != c2.cols() ) 
		throw ColorMatrix::BadSize("Error: Painting::difference() - images are not the same size.");

	// For each pixel, calculate the difference
	for(unsigned int row = 0; row < diff.rows(); ++row)
	for(unsigned int col = 0; col < diff.cols(); ++col) {
		const Color cc1 = c1(row, col);
		const Color cc2 = c2(row, col);

		// Difference = sqrt((r1-r2)^2 + (g1-g2)^2 + (b1-b2)^2)
		const Component c = clampf( abs(sqrt(static_cast<double>(
			(cc1.r() - cc2.r()) * (cc1.r() - cc2.r()) 
		  + (cc1.g() - cc2.g()) * (cc1.g() - cc2.g())
		  + (cc1.b() - cc2.b()) * (cc1.b() - cc2.b()) ))) );
		diff(row, col) = Color(c, c, c, 255);
	}
}

/* sumDifferences()
 * ----------------
 * Calculates the sum of intensity values in a difference image
 * using only the red component of a particular pixel 
 * (because difference images are grayscale)
 * @param diff - the difference image
 * @return - the sum of all the difference intensities in the image
 */
double Painting::sumDifferences(const ColorMatrix& diff)
{
	double result = 0.0;

	for(unsigned int row = 0; row < diff.rows(); ++row)
	for(unsigned int col = 0; col < diff.cols(); ++col) {
		const Color pixel = diff(row, col);
		
		// Ignore pixels marked as transparent, 
		// as they are outside the image bounds
		if( pixel.a() != 0 )
			result += pixel.r();
	}

	return result; 
}

// TODO - copyRegion() should probably be moved to ImageProcessing namespace
/* copyRegion()
 * ------------
 * Copies values from the source image to the dest image 
 * from the region source[x-(grid/2)..x+(grid/2), y-(grid/2)..y+(grid/2)].
 * Note: for values that would fall outside the source image, 
 * dest pixels are marked as transparent black.
 * @param source - the source of the color data to copy
 * @param x,y - the center indices of the region to copy
 * @param grid - the diameter of the region to copy
 * @param dest - the destination image to copy color data to 
 *				 (requires size of [grid,grid])
 */
void Painting::copyRegion(const ColorMatrix& source, 
					const unsigned int x, 
					const unsigned int y, 
					const unsigned int grid, 
					ColorMatrix& dest)
{
	// Validate dest size
	if( dest.rows() != grid || dest.cols() != grid )
		dest = ColorMatrix(grid, grid);

	for(unsigned int row = 0; row < grid; ++row)
	for(unsigned int col = 0; col < grid; ++col)
	{
		// Calculate source pixel index corresponding to region index
		const int srow = y - (grid / 2) + row;
		const int scol = x - (grid / 2) + col;

		// If pixel index is out of bounds of the source image, 
		// Set the corresponding region color value to transparent 
		if( srow < 0 || srow >= static_cast<int>(source.rows()) 
		 || scol < 0 || scol >= static_cast<int>(source.cols()) )
			dest(row, col) = Color(0,0,0,0);
		else { // Copy the pixel from source to dest
			dest(row, col) = source(srow, scol);
			dest(row, col).a(255);
		}
	}
}

/* maxErrorRowIndex()
 * ------------------
 * Calculates the row index of the maximum 
 * error value from the specified difference region
 * @param diff - the difference region to examine 
 * @return - the row index to the max error value from the difference region
 */
unsigned int Painting::maxErrorRowIndex(const ColorMatrix& diff)
{
	unsigned int rowIndex = 0;
	Component    maxError = 0;
	for(unsigned int row = 0; row < diff.rows(); ++row)
	for(unsigned int col = 0; col < diff.cols(); ++col)
	{
		const Component error = diff(row, col).r();
		if( error > maxError ) {
			maxError = error;
			rowIndex = row;
		}
	}
	return rowIndex;
}

/* maxErrorColIndex()
 * ------------------
 * Calculates the column index of the maximum 
 * error value from the specified difference region
 * @param diff - the difference region to examine 
 * @return - the column index to the max error value from the difference region
 */
unsigned int Painting::maxErrorColIndex(const ColorMatrix& diff)
{
	unsigned int colIndex = 0;
	Component    maxError = 0;
	for(unsigned int row = 0; row < diff.rows(); ++row)
	for(unsigned int col = 0; col < diff.cols(); ++col)
	{
		const Component error = diff(row, col).r();
		if( error > maxError ) {
			maxError = error;
			colIndex = col;
		}
	}
	return colIndex;
}

/* makeStroke()
 * ------------
 * @param brush - the brush to use for this stroke
 * @param x,y - the pixel indices to center the stroke on
 * @param reference - the reference image to get a color value from 
 * @return - the newly created Stroke object
 */
Stroke Painting::makeStroke(const Brush& brush, 
							const unsigned int x, 
							const unsigned int y, 
							const ColorMatrix& reference)
{
	return Stroke(x, y, brush, reference(y, x));
}

/* generate()
 * ----------
 * A public method to generate a new painterly rendering of 
 * this Painting object's current source image using 
 * the current set of Brushs and specified constant values
 * @param gridFactor - is the grid scale factor for each layer
 * @param errorT - is the error threshold value for each layer
 */
void Painting::generate(const vector<unsigned int>& brushesAndSizes, 
						const double gridFactor, const double errorT)
{
	gridSizeFactor = gridFactor;
	errorThreshold = errorT;
	
	// Generate a brush set based on interleaved type/size data
	brushes.clear();
	if( !brushesAndSizes.empty() ) {  // Use brushesAndSizes data
		vector<unsigned int>::const_iterator it  = brushesAndSizes.begin();
		vector<unsigned int>::const_iterator end = brushesAndSizes.end();
		for(; it != end; ++it) 
		{
			const unsigned int b = *it;
			const unsigned int s = (++it != end) ? *it : 8;	
			const Brush::Type  type = static_cast<Brush::Type>(b);
			brushes.push_back(Brush(s, type));
		}
	} else { // Use default brush set
		brushes.push_back(Brush(8,Brush::circle));
		brushes.push_back(Brush(4,Brush::circle));
		brushes.push_back(Brush(2,Brush::circle));
	}

	paint();
}