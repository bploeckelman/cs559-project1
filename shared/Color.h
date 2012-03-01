#pragma once
/* ----------------------------------------------
	Color class, stores rgba components
	Author: Brian Ploeckelman
   ---------------------------------------------- */
typedef unsigned char Component;	// Each color component is a byte

class Color
{
private: 
	Component r_, g_, b_, a_;

public:
	Color() : r_(0), g_(0), b_(0), a_(255) { }
	Color(Component r, Component g, Component b, Component a)
		: r_(r), g_(g), b_(b), a_(a) { }

	// Get component values
	Component r() const { return r_; }	
	Component g() const { return g_; }	
	Component b() const { return b_; }	
	Component a() const { return a_; }	

	// Set component values (can be chained)
	Color& r(Component r) { r_ = r; return *this; }
	Color& g(Component g) { g_ = g; return *this; }
	Color& b(Component b) { b_ = b; return *this; }
	Color& a(Component a) { a_ = a; return *this; }
}; 

// Pre-defined Colors
static const Color Red    () { static Color c(255,   0,   0, 255); return c; }
static const Color Orange () { static Color c(255, 165,   0, 255); return c; }
static const Color Yellow () { static Color c(255, 255,   0, 255); return c; }
static const Color Green  () { static Color c(  0, 255,   0, 255); return c; }
static const Color Blue   () { static Color c(  0,   0, 255, 255); return c; }
static const Color Indigo () { static Color c( 75,   0, 130, 255); return c; }
static const Color Violet () { static Color c(238, 130, 238, 255); return c; }
static const Color Magenta() { static Color c(255,   0, 255, 255); return c; }
static const Color Black  () { static Color c(  0,   0,   0, 255); return c; }
static const Color White  () { static Color c(255, 255, 255, 255); return c; }