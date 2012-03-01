CS559 - Project 1 - Painterly program
Author: Brian Ploeckelman
Date:   2-27-2012
--------------------------------

Usage: 
------
painterly INFILE OUTFILE [G T [N S N S ..]]
[G] - the grid scale factor
none: 0.25 (default)
[T] - the error threshold
none: 75 (default)
[N S] - the brush type to use, and its size
     0: square
	 1: circle
	 2: soft circle
	 3: back slash
	 4: forward slash
	 5: x
	 6: abstract
	 7: random stroke
  none: 3 circles, sized 8 4 2 (ie. 7 8 7 4 7 2) (default)


This program takes an input targa image and uses the Hertzmann painterly 
rendering algorithm to make it look like it was painted.  The result is 
saved as a targa file with filename specified by OUTFILE.  As the program 
is running, text is displayed informing the user about what processing the 
program is currently doing.

The program allows the user to specify a variety of options that define 
how the "painting" is accomplished, all of these extra options are optional.

G - the grid scale factor, default value 0.25
T - the error threshold, default value 75

These two options, if specified, can be followed by a series of pairs of 
integral numbers.  The first number in a pair specifies a type of brush, 
the second number in a pair specifies the size of its corresponding brush.

For example: 
a 16 pixel random stroke, followed by 
a 10 pixel forward slash, followed by 
an 8 pixel soft circle

painterly input.tga output.tga 0.25 75 7 16 4 10 2 8

As many brushes as desired can be specified this way, and they will be 
applied to the canvas in the order specified, however using many brushes 
will result in the program taking longer to complete (but can make for 
very interesting results).

1) The painting algorithm is a straighforward implementation of the Hertzmann 
algorithm as specified in the paper "Painterly Rendering with Curved Brush 
Strokes of Multiple Sizes" (Hertzmann, 1998). The curved brush stroke 
technique described in the paper is not implemented.   

The way the program decides where to place brush strokes is by calculating 
a difference image between the current state of the canvas and a reference 
image that is blurred with a filter whose strength corresponds to the size 
of the brush used for the current layer.  Then overall error is calculated 
for many regions by summing difference intensities in a region whose size 
corresponds to the size of the current brush multiplied by the grid scaling 
factor for the painting. If this overall region error (divided by the 
square of the grid size) is above the threshold for the painting, then a 
stroke is created that is centered on the corresponding region.  All strokes 
are then randomly applied to the canvas, completing the current layer.

The way the program creates brush strokes is based on the user's selection of 
brush types and sizes (or the default values if the user provides none).  
Note: the "random stroke" option selects randomly between the 7 different 
variants of brush-squiggleX.tga when creating new strokes, the other brush 
types behave as expected for new strokes.

The look of the brush stroke depends on the specified style; but in general 
brushes are created from the pixel intensities of a brush targa image 
that corresponds to the type specified.

2/3) Given the variety of brush types which can be used, and the essentially 
unlimited order in which they can be applied, the program is capable of 
producing interesting results for nearly any type of input.  However it may 
take some experimentation to see which combination of brush types and sizes 
produce the nicest results.

The program supports arbitrary brushes defined as targa files where white 
indicates a place where color should be applied, other pixels are transparent. 
See the targa files in the Build/brushes subdirectory for examples of how 
to create these files.  While arbitrary brushes can be created easily, the 
program does not yet support adding new types of brushes without making 
changes in the code.  Future versions of this program could build a list of 
brushes in the appropriate subdirectory and allow the user to pick from them.

Note: the slowest part of the rendering technique is the calculation of the 
difference image for each layer, so if it seems like it is sitting at that step 
for a long time, be patient :)

