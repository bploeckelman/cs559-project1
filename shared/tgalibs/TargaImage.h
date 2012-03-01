#pragma once
// simple C++ class wrapper for the data in LibTarga - designed to
// make things a little bit easier.
//
// class designed by Mike Gleicher 9/5/2007
// be sure to give proper attribution for this (and any other code)
// that you use in your projects!
//

// design ideas:
// - we only provide the barest minimum of functionality. students
//   should extend this class (if they want to) by subclassing it
//   (although that has issues with the construction style), or adding
//   to it.
// - we don't use any fancy C++ stuff (exceptions, templates, ...), so
//   we use C style error checking
// - we don't want to have "invalid" images, so there is no default
//   constructor (in fact, we'll protect it from accidental use)
// - we don't want to copy the images (since they're big) except when
//   we know what we're doing, so we don't have a copy constructor
// - there are static functions that dynamically create instances,
//   this way they can fail (for example when loading from a file),
//   rather than creating an invalid image. the constructors are
//   private, any construction must be done by the creator functions
// - because we generally will allocate them dynamically, it will be
//   the callers responsibility to delete things when finished
// - the above few things are stylistically very "C" like, but this is
//   a wrapper on a C library
// - to make life easier, we'll always assume that we want 4 bits per
//   pixel. when a 3 channel image is loaded in, the alpha value will
//   be set to 100% (e.g. 255) for each pixel
// - 0,0 refers to the top left pixel. this is NOT the way that the
//   targa file stores things. things will be flipped when they get
//   read in, and they will be written out correctly
// - the pixels are stored densely - there is no row padding
// - we copy the data that is passed into the constructor. this isn't
//   necessary (we could have the class take ownership) for the current
//   functionality, but could simplify things in the future.

class TargaImage {
  // constructors are private - use the construction functions!
 public:
  //// creating and destroying
  // creates an empty (black, alpha=0) image
  static TargaImage* blankImage(const unsigned int width,
								const unsigned int height);
  // reads a targa image from a file- returns a pointer to an image if
  // (and only if) the read was successful. returns a NULL pointer if
  // there is an error. use the targa error check to get the error
  // message 
  static TargaImage* readImage(const char* filename);

  // destuctor - remember to call this on anything you create
  ~TargaImage();

 public:
  //// accessors (to the private data) - you should never change these
  //   things after it is created!
  inline unsigned int width() const { return _width; };
  inline unsigned int height() const { return _height; };
  inline unsigned char* pixels() const { return _pixels;};

  //// utility routines - we only provide a very basic set
  // flip the image in place (re-orders the pixels)
  void flip();
  // write the image to a file - returns 1 if there is success,
  // otherwise check the LibTarga error code
  int write(const char* filename) const;
  // set all the alphas to a specific value
  void fillAlpha(const unsigned char value);

  // this takes the place of a copy constructor. it creates a "deep"
  // copy of the image. the caller is responsible for deleting what
  // gets created
  TargaImage* clone() const;

  // constructors are private - use the construction functions!
 private:
  // this constructor is useful if you used the targa reading function
  // (or the targa blank creator). these are really the only ways to
  // get new images
  // note: this allocates its own memory, so make sure you free the thing
  //       that you pass into this
  TargaImage(const unsigned int width, const int height, 
			 const unsigned char* pixels);
  // this null constructor should never be called - it is here to make
  // sure it is private and cannot be called accidentally
  TargaImage();

  // data members are private
 private:
  unsigned int _width;
  unsigned int _height;
  unsigned char* _pixels;		/* pointer to width*height*4 bytes */
};