#pragma once
/* ----------------------------------------------
	Matrix2d is a template class that implements 
	a two-dimensional matrix of a particular 
	type of object.

	Based on example from c++ faq lite section [16.19]
	@ http://www.parashift.com/c++-faq-lite/

	Adapted by: Brian Ploeckelman
	Added: 
	- more descriptive exception handling 
	- fill() method
	- getRow/Col() methods
	- setRow/Col() methods
    --------------------------------------------- */
#include <stdexcept>
#include <sstream>
#include <string>
#include <vector>


template<typename T>
class Matrix2d
{
public:
	// Standard constructor
	Matrix2d(unsigned int numRows, unsigned int numCols);
	// Throws a BadSize object if either size is zero
	class BadSize : public std::runtime_error 
	{
	public:
		BadSize(const std::string& what_arg)
			: std::runtime_error(what_arg) { }
	};

	// Access methods to get the (row, col)th entry
	T&		 operator()(unsigned int row, unsigned int col);
	T const& operator()(unsigned int row, unsigned int col) const;
	// Throws a BoundsViolation object if i or j are out of bounds
	class BoundsViolation : public std::runtime_error 
	{ 
	public:
		BoundsViolation(const std::string& what_arg)
			: std::runtime_error(what_arg) { }
	};

	// Utility methods
	void fill(const T& data);

	// Fill a specified row or column with the given values
	void setRow(const unsigned int row, const std::vector<T>& rowValues);
	void setCol(const unsigned int col, const std::vector<T>& colValues);

	// Get a copy of a specified row or column
	const std::vector<T> getRow(const unsigned int row) const;
	const std::vector<T> getCol(const unsigned int col) const;

	// Get the number of rows or columns
	const unsigned int rows() const; 
	const unsigned int cols() const;

private:
	std::vector< std::vector<T> > data_;
};

// Include inline function definitions 

/* Standard constructor */
template <typename T>
inline Matrix2d<T>::Matrix2d(unsigned int numRows, unsigned int numCols)
	: data_(numRows)
{
	if( numRows == 0 || numCols == 0 )
		throw BadSize("Error: Matrix2d - 0 is an invalid dimension.");
	for(unsigned int i = 0; i < numRows; ++i)
		data_[i].resize(numCols);
}

/* Element access, non-const */
template<typename T>
inline T& Matrix2d<T>::operator()(unsigned int row, unsigned int col) 
{
	if( row >= rows() || col >= cols() ) 
		throw BoundsViolation("Error: Matrix2d element access - index out of bounds.");
	return data_[row][col];
}

/* Element access, const */
template<typename T>
inline T const& Matrix2d<T>::operator()(unsigned int row, unsigned int col) const 
{
	if( row >= rows() || col >= cols() ) 
		throw BoundsViolation("Error: Matrix2d element access - index out of bounds.");
	return data_[row][col];
}

/* Fill the Matrix2d with the specified data */
template <typename T>
inline void Matrix2d<T>::fill(const T& data) 
{
	for(unsigned int row = 0; row < rows(); ++row)
	for(unsigned int col = 0; col < cols(); ++col)
		(*this)(row,col) = data;
}

/* Sets values for the specified row to the values in rowValues */
template <typename T>
inline void Matrix2d<T>::setRow(const unsigned int row, const std::vector<T>& rowValues)
{
	// Validate row index and number of values
	if( row >= rows() || rowValues.size() > cols() )
		throw BoundsViolation("Error: Matrix2d::setRow() - index out of bounds.");

	// Copy values into specified row
	for(unsigned int i = 0; i < cols(); ++i)
		(*this)(row,i) = rowValues[i];
}

/* Sets values for the specified col to the values in colValues */
template <typename T>
inline void Matrix2d<T>::setCol(const unsigned int col, const std::vector<T>& colValues)
{
	// Validate column index and number of values
	if( col >= cols() || colValues.size() > rows() )
		throw BoundsViolation("Error: Matrix2d::setCol() - index out of bounds.");

	// Copy values into specified column
	for(unsigned int i = 0; i < rows(); ++i)
		(*this)(i,col) = colValues[i];
}

/* Returns the data from a specific row */
template <typename T>
inline const std::vector<T> Matrix2d<T>::getRow(const unsigned int row) const
{
	if( row >= rows() ) 
		throw BoundsViolation("Error: Matrix2d::getRow() - index out of bounds.");
	return data_[row];
}

/* Returns the data from a specific column */
template <typename T>
inline const std::vector<T> Matrix2d<T>::getCol(const unsigned int col) const 
{
	if( col >= cols() )
		throw BoundsViolation("Error: Matrix2d::getCol() - index out of bounds.");

	// Create a new vector for column data
	std::vector<T> column;	

	// Fill the vector with data from the specified column
	for(unsigned int i = 0; i < rows(); ++i)
		column.push_back((*this)(i,col));

	return column;
}

/* Number of rows */
template <typename T>
inline const unsigned int Matrix2d<T>::rows() const 
{
	return data_.size();
}

/* Number of columns */
template<typename T>
inline const unsigned int Matrix2d<T>::cols() const
{
	return data_[0].size();
}