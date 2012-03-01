/* ---------------------------------------------------------------------------
 * CS559 - Project 1 - Image painterly rendering
 * Brian Ploeckelman
 * 
 * This program takes a targa image as input and renders it as a painting, 
 * saving the resulting color data back out as a targa file.
 */
#include "ImageProcessing.h"
#include "ImagePainting.h"

#include <iostream>
#include <vector>
#include <string>
#include <conio.h> // for _getch()


int main(int argc, char* argv[])
{
	using namespace ImageProcessing;
	using std::vector;
	using std::string;
	using std::cout;
	using std::endl;

	// Validate command line args
	if( argc < 3 || (argc > 3 && argc < 5) ) {  
		cout << "Usage: " 
			 << "painterly INFILE OUTFILE [G T [N S N S..]]"<< endl 
			 << "[G] - the grid scale factor"		<< endl
			 << "none: 0.25 (default)"		<< endl << endl 
			 << "[T] - the error threshold"			<< endl 
			 << "none: 75 (default)"		<< endl << endl 
			 << "[N S] - the brush type to use, and its size:" << endl
			 << "   0: square"						<< endl
			 << "   1: circle"						<< endl
			 << "   2: soft circle"					<< endl
			 << "   3: backslash"					<< endl
			 << "   4: forward slash"				<< endl
			 << "   5: x"							<< endl 
			 << "   6: abstract"					<< endl
			 << "   7: random stroke"				<< endl 
			 << "none: 3 circles, sized 8 4 2"      << endl 
			 << "  ie: 7 8 7 4 7 2(default)"<< endl << endl
			 << "Press any key to quit..."		    << endl;
		_getch();
		return 1;
	}

	// Get command line args
	const string infile  = argv[1];
	const string outfile = argv[2];
	const double gridFactor		= (argc > 3) ? atof(argv[3]) : 0.25;
	const double errorThreshold = (argc > 4) ? atoi(argv[4]) : 100.0;

	// Handle brush-size pair command line args
	vector<unsigned int> brushesAndSizes;
	for(int i = 5; i < argc; i += 2) { 
		const unsigned int n = (argc >= i)   ? atoi(argv[i])   : 1; // default brush = circle
		const unsigned int s = (argc >= i+1) ? atoi(argv[i+1]) : 8; // default size = 8
		brushesAndSizes.push_back(n);
		brushesAndSizes.push_back(s);
	}

	try { 
		cout << "Loading input image: " << infile << "..." << endl;
		ColorMatrix source( loadTargaColorMatrix(infile) );

		// Generate the painterly rendering of the source image
		Painting painting(source);
		painting.generate(brushesAndSizes, gridFactor, errorThreshold);

		cout << "Saving painted image..." << endl;
		bool success = saveTargaColorMatrix(outfile, painting.canvasCopy());
		cout << (success ? "Targa image saved: " : "Unable to save targa image: ") 
			 << outfile << endl;
	} 
	catch(const ColorMatrix::BadSize& e)			{ cout << e.what() << endl; } 
	catch(const KernelMatrix::BadSize& e)			{ cout << e.what() << endl; } 
	catch(const ColorMatrix::BoundsViolation& e)	{ cout << e.what() << endl; }
	catch(const FailedToLoad& e)					{ cout << e.what() << endl; }
	catch(const BadPixelData& e)					{ cout << e.what() << endl; }
	catch(const UnsupportedKernelType& e)			{ cout << e.what() << endl; } 
	catch(const UnsupportedBrushType& e)			{ cout << e.what() << endl; } 

	cout << "Press any key to quit..." << endl;
	_getch();

	return 0;
}