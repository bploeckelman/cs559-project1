#pragma once
/* ---------------------------------------------------------------------------
 * The ImagePainting namespace is a container for a number 
 * of methods that perform a painterly-rendering algorithm on an image
 * 
 * CS559 - Project 1 - Painterly
 * Brian Ploeckelman
 * ---------------------------------------------------------------------------
 */
#include "ImageProcessing.h"

#include <stdexcept>
#include <string>
#include <list>


/* _______________________________________________ */
/* -------------------- Brush -------------------- */
class Brush
{
public:
	enum Type { 
		square = 0, circle, feather, 
		backslash, foreslash, x, 
		Abstract, squiggle
	};

private:
	Type		type;
	ColorMatrix weights;
	unsigned int diameter;

public:
	Brush(const unsigned int diameter=8, const Type& type=squiggle); 

	// Accessors
	unsigned int size()    const { return diameter; } 
	Type		 getType() const { return type;     }
	Color		 getWeight(unsigned int row, unsigned int col) const 
	{ return weights(row, col); }

	/* Converts the grayscale color value at (row,col) 
	 * to an intensity value [0.0..1.0] 
	 */
	double weight(unsigned int row, unsigned int col) const {
		if( row > diameter || col > diameter ) 
			return 0.0; 
		else
			return weights(row, col).r() / 255.0;
	}
};
// Brush constructor throws UnsupportedBrushType exception
class UnsupportedBrushType : public std::runtime_error 
{
public:
	UnsupportedBrushType(const std::string& what_arg) 
		: std::runtime_error(what_arg) { }
};


/* ________________________________________________ */
/* -------------------- Stroke -------------------- */
class Stroke
{
private:
	unsigned int x,y;
	Brush brush;
	Color color;

public:
	Stroke(const unsigned int x, 
		   const unsigned int y, 
		   const Brush& brush, 
		   const Color& color);

	void Stroke::apply(ColorMatrix& canvas);

	unsigned int getX() const { return x; }
	unsigned int getY() const { return y; }
	const Brush& getBrush() const { return brush; }
	const Color& getColor() const { return color; }
};


/* __________________________________________________ */
/* -------------------- Painting -------------------- */
class Painting
{
private:
	static const bool debug;	// If true: saves intermediate images

	ColorMatrix source;
	ColorMatrix canvas;
	std::list<Brush> brushes;

	double gridSizeFactor;
	double errorThreshold;

	void paint();
	void paintLayer(const ColorMatrix& reference, const Brush& brush);
	
	void difference(const ColorMatrix& c1, 
					const ColorMatrix& c2, 
					ColorMatrix& diff);
	double sumDifferences(const ColorMatrix& diff);

	// TODO - consolidate maxErrorRow/ColIndex() methods into 1 method
	unsigned int maxErrorRowIndex(const ColorMatrix& diff);
	unsigned int maxErrorColIndex(const ColorMatrix& diff);

	Stroke makeStroke(const Brush& brush, 
					  const unsigned int x, 
					  const unsigned int y, 
					  const ColorMatrix& reference);
	
	void copyRegion(const ColorMatrix& source, 
					const unsigned int x, 
					const unsigned int y, 
					const unsigned int grid, 
					ColorMatrix& dest);

public:
	Painting(const ColorMatrix& source);

	void generate(const std::vector<unsigned int>& brushesAndSizes, 
				  const double gridFactor=0.25, const double errorT=75);

	const ColorMatrix& sourceCopy() const { return source; }
	const ColorMatrix& canvasCopy() const { return canvas; } 
};