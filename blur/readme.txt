CS559 - Project 1 - Blur program
Author: Brian Ploeckelman
Date:   2-27-2012
--------------------------------

Usage:
------
blur INFILE OUTFILE R [N - optional]
[R] - the blur kernel radius
[N] - the kernel type specifier:
   0: box (2-d slow)
   1: gaussian
   2: b-spline
   3: catmull-rom spline
   4: mitchell-netravali cubic (with default B=C=1/3)
   5: gaussian (2-d slow)
none: mitchell-netravali cubic (default)


This program blurs an input targa file with a kernel of a specified radius 
and saves the result as a targa file with filename specified by OUTFILE. 
As the program is running, text is displayed informing the user about what 
processing the program is currently doing.

The program allows users to select between 6 different blurring kernels, 
and if no kernel type is specified at the command-line, then the default 
mitchell-netravali cubic kernel with the specified radius is used.

Supported kernel types:
box, gaussian (1d + 2d), b-spline, catmull-rom spline, mitchell-netravali cubic

Image boundaries are handled by ignoring values outside the image, 
although the program is extensible enough that with a bit more time 
other boundary-handling techniques could be implemented