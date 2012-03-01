/* ---------------------------------------------------------------------------
 * CS559 - Project 1 - Image resizing 
 * Brian Ploeckelman
 * 
 * This program takes a targa image as input and resizes it, 
 * saving the resulting color data back out as a targa file.
 */
#include "ImageProcessing.h"

#include <iostream>
#include <string>
#include <conio.h> // for _getch()


int main(int argc, char* argv[])
{
	using namespace ImageProcessing;
	using std::string;
	using std::cout;
	using std::endl;

	// Validate command line args
	if( argc < 5 || argc > 6 ) { 
		cout << "Usage: " 
			 << "resize INFILE OUTFILE X Y [N - optional]"           << endl 	
			 << "[X Y] - the new dimensions to resize the image to." << endl 
			 << "[N]   - the type of resizing to perform: "          << endl 
			 << "     0: nearest-neighbor"                           << endl 
		     << "     1: bilinear"                                   << endl 
			 << "  none: bilinear (default)"                 << endl << endl
			 << "Press any key to quit..."                           << endl;
		_getch();
		return 1;
	}

	// Get command line args
	const string infile  = argv[1];
	const string outfile = argv[2];
	const int    cols    = atoi(argv[3]);
	const int    rows    = atoi(argv[4]);
	const int n = (argc == 6) ? atoi(argv[5]) : 1; // 1 - bilinear default

	// Validate new dimensions
	if( rows < 1 || cols < 1 ) {
		cout << "Error: new dimensions must be non-zero positive values." << endl
			 << "Press any key to quit..." << endl;
		_getch(); 
		return 1;
	}

	// Validate filter type
	if( n != 0 && n != 1 ) { 
		cout << "Error: filter type can only be 0 (nearest) or 1 (bilinear)." << endl
			 << "Press any key to quit..." << endl;
		_getch();
		return 1;
	}

	try { 
		cout << "Loading input image: " << infile << "..." << endl;
		ColorMatrix source( loadTargaColorMatrix(infile) );

		cout << "Resizing the input image..." << endl;
		ColorMatrix output(1,1);
		if( n == 0 )
			output = resizeNearest(source, rows, cols);
		else
			output = resizeBilinear(source, rows, cols);

		cout << "Saving resized image..." << endl;
		bool success = saveTargaColorMatrix(outfile, output);
		cout << (success ? "Targa image saved: " : "Unable to save targa image: ") 
			 << outfile << endl;
	}
	catch(const ColorMatrix::BadSize& e)			{ cout << e.what() << endl; } 
	catch(const KernelMatrix::BadSize& e)			{ cout << e.what() << endl; } 
	catch(const ColorMatrix::BoundsViolation& e)	{ cout << e.what() << endl; }
	catch(const FailedToLoad& e)					{ cout << e.what() << endl; }
	catch(const BadPixelData& e)					{ cout << e.what() << endl; }
	catch(const UnsupportedKernelType& e)			{ cout << e.what() << endl; } 

	cout << "Press any key to quit..." << endl;
	_getch();

	return 0;
}