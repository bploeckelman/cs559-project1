CS559 - Project 1 - Resize program
Author: Brian Ploeckelman
Date:   2-27-2012
--------------------------------

Usage:
------
resize INFILE OUTFILE X Y [N - optional]
[X Y] - the new dimensions for the image (must be > 0)
[N]   - the type of resizing to perform
     0: nearest-neighbor
	 1: bilinear
  none: Bilinear (default)


This program resizes an input targa file to the new dimensions specified 
at the command-line, and saves the result as a targa file with the filename 
specified by OUTFILE. As the program is running, text is displayed informing 
the user about what processing the program is currently doing.

The program supports nearest-neighbor and bilinear filtering.  
Many other kernels are available in the code (same set as is available for 
the blurring program), but resampling using these kernels has not yet 
been implemented.

Image boundaries are handled by ignoring values outside the image, 
although the program is extensible enough that with a bit more time 
other boundary-handling techniques could be implemented.

The result images tagged as "-bad" use nearest-neighbor sampling, 
whereas the images without that tag use bilinear sampling.