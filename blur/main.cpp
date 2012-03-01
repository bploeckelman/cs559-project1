/* ---------------------------------------------------------------------------
 * CS559 - Project 1 - Image blurring 
 * Brian Ploeckelman
 * 
 * This program takes a targa image as input, and blurs it 
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
	if( argc < 4 || argc > 5 ) { 
		cout << "Usage: " 
			 << "blur INFILE OUTFILE R [N - optional]" << endl 
			 << "[R] - the blur kernel radius "        << endl 
			 << "[N] - Kernel type specifier: "        << endl 
			 << "   0 = box (slow)"                    << endl 
			 << "   1 = gaussian"                      << endl  
			 << "   2 = b-spline"                      << endl 
			 << "   3 = catmull-rom spline"            << endl 
			 << "   4 = mitchell-netravali"            << endl 
			 << "   5 = gaussian (2-d slow)"           << endl 
			 << "none = mitchell-netravali (default)"  << endl;
		cout << endl << "Press any key to quit..."     << endl;
		_getch();
		return 1;
	}

	// Get command line args
	const string infile  = argv[1];
	const string outfile = argv[2];
	const int    r       = atoi(argv[3]);
	
	// Handle the optional command-line arg, use default kernel type if no arg
	int	n = (argc == 5) ? atoi(argv[4]) : 4;  // 4 == mitchell-netravali
	KernelType kernelType = static_cast<KernelType>(n);

	try { 
		cout << "Loading input image: " << infile << "..." << endl;
		ColorMatrix source( loadTargaColorMatrix(infile) );

		cout << "Blurring the input image..." << endl;
		ColorMatrix output(1,1);
		switch(kernelType)
		{
		case gaussian: sep_convolve(source, output, generate1dGaussian(r)); 	break;
		case bspline:  sep_convolve(source, output, generate1dBSpline(r));		break;
		case catmull:  sep_convolve(source, output, generate1dCatmull(r));		break;
		case mitchell: sep_convolve(source, output, generate1dMitchell(r));		break;
		default:   	   convolve(source, output, generateBlurKernel(r, kernelType));
		}

		cout << "Saving blurred image..." << endl;
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