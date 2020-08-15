
#ifndef VECTOR_H
#define VECTOR_H

#include <stddef.h>
#include <stdbool.h>

#define ROW true
#define COL false
#define EPSILON 0.0000000000001 // value used for double comparison

// A Vector is a wrapped one dimensional array that is either a 
// row or a column, depending on the value of type (using macros defined above).
// So, a Vector is a member of either R^n or R^(1Xn).
typedef struct Vector {
    size_t dim;
    bool type;
    double *elements;
} Vector;

// Create a Vector that is initialized to all zeroes
Vector make_vector(size_t dim, bool type);

// Provide a deep copy of given Vector
Vector copy_vector(Vector vec);

// Create a Vector using the given array of vals
Vector make_vector_vals(double *vals, size_t dim, bool type);

// Create a Vector filled with random values in the range [-max_abs_val, max_abs_val]
Vector make_vector_rand(size_t dim, bool type, unsigned int max_abs_val);

// Create a new Vector that deletes all repeat values from given Vector
Vector make_vector_unique(Vector vec);

// Delete memory associated with a Vector. Must be done with every created Vector
void delete_vector(Vector vec);

// Append vec2 to vec1, given that they are the same type
Vector combine(Vector vec1, Vector vec2);

// Print the contents of given Vector to the console. It will be printed either
// horizontally or vertically depending on if it is a row or a column
void print_vector(Vector vec);

// Fill the given Vector with zeroes
void zeroes_vec(Vector vec);

// Fill the given Vector with ones
void ones_vec(Vector vec);

// Fill the given Vector with the value
void fill_vec(Vector vec, double value);

// Set the element of the Vector at index to the given value
void set_element_vec(double value, Vector vec, size_t index);

// Set a subsection of vec to the given sub Vector
void set_subvector(Vector sub_vec, Vector vec, size_t index);

// Get the element in the Vector at index
double get_element_vec(Vector vec, size_t index);

// Get an array containing the values of the given Vector
double * get_elements_vec(Vector vec);

// Get the sub Vector starting at index
Vector get_subvector(Vector vec, size_t index);

// Get the index of the pivot of the given Vector. If the Vector
// has no pivot, the dimension of the Vector is returned
size_t get_pivot_index(Vector vec);

// Return a singleton NULL_VECTOR that is used akin to NULL ptr
Vector get_null_vec(void);

// Get the n dimensional standard basis. For dim == 3:
//      1       0       0
//      0       1       0
//      0       0       1
Vector * get_ndim_standard_basis(size_t dim);

// Get the norm of a Vector, which is the square root of the
// Vector dotted with itself
double get_norm_vec(Vector vec);

// Get the angle between two Vectors, using the fact that
// cos(theta) == dot(vec1, vec2) / ||vec1||*||vec2||
// Angle is returned in degress
double get_angle_vec(Vector vec1, Vector vec2);

// Get the angle between two Vectors, using the fact that
// cos(theta) == dot(vec1, vec2) / ||vec1||*||vec2||
// Angle is returned in radians
double get_angle_rads_vec(Vector vec1, Vector vec2);

// Get the distance between two vectors. Distance is defined
// as ||vec2 - vec1||
double get_distance_vec(Vector vec1, Vector vec2);

// Get the unit vector pointed in the direction of the given Vector.
// This is defined as vec / ||vec||
Vector get_unit_vec(Vector vec);

// Get the scalar projection of vec1 onto vec2
double get_scalar_proj_vec(Vector vec1, Vector vec2);

// Get the vector projection of vec1 onto vec2
Vector get_vector_proj_vec(Vector vec1, Vector vec2);

// Get the projection of vec onto the span of basis, which is 
// the Vector in basis that is closest to vec (min distance)
Vector get_projection_vec(Vector *basis, size_t dim, Vector vec);

// Check if index is in range of Vector's dimension
bool in_range_vec(Vector vec, size_t index);

// Check if two Vectors have the same dimension
bool same_size_vec(Vector vec1, Vector vec2);

// Check if two Vectors have the same type (ie ROW and COL)
bool same_type(Vector vec1, Vector vec2);

// Check to see if two Vectors are of the same type and have
// the same dimension
bool same_size_and_type(Vector lhs, Vector rhs);

// Check if two Vectors are equal, which implies that they have the same type
// and dimension, and that Ai == Bi for all i
bool equal_vec(Vector lhs, Vector rhs);

// Check if two Vectors are capable of being dotted
bool can_dot(Vector lhs, Vector rhs);

// Determine if a set of Vectors is linearly independent
bool are_lin_ind(Vector *vecs, size_t dim);

// Determine if a set of Vectors is linearly dependent
bool are_lin_dep(Vector *vecs, size_t dim);

// Determine if a set of Vectors is a basis of R^n
bool is_ndim_basis(Vector *vecs, size_t dim);

// Determine if two Vectors are orthogonal, meaning dot(vec1, vec2) == 0
bool are_orthogonal_vec(Vector vec1, Vector vec2);

// Determine if every Vector in vecs1 is orthogonal to every vector in vecs2
bool are_orthogonal_sets_vec(Vector *vecs1, Vector *vecs2, size_t dim1, size_t dim2);

// Determine if every Vector in vecs is orthogonal to every other Vector in vecs
bool is_orthogonal_set_vec(Vector *vecs, size_t dim);

// Determine if every Vector in vecs is orthogonal to every other Vector in vecs
// and whether every Vector in vecs has norm == 1
bool is_orthonormal_set_vec(Vector *vecs, size_t dim);

// Determine if a Vector is a unit Vector, meaning that norm == 1
bool is_unit_vec(Vector vec);

// Check to see if Vector is a row (ie type == ROW)
bool is_row(Vector vec);

// Check to see if Vector is a column (ie type == COL)
bool is_col(Vector vec);

// Check to see if every element in vec is zero
bool is_zero_vec(Vector vec);

// Check to see if vec is the singleton NULL_VECTOR
bool is_null_vec(Vector vec);

// Return transpose of the given Vector, which just entails changing
// its type from ROW->COL or COL->ROW
Vector transpose_vec(Vector vec);

// Perform element wise addition of two Vectors
Vector add_vec(Vector lhs, Vector rhs);

// Perform element wise subtraction of two Vectors
Vector subtract_vec(Vector lhs, Vector rhs);

// Multiply every element in Vector by given scalar, return result in new Vector
Vector smultiply_vec(Vector vec, double scalar);

// Multiply every element in Vector by given scalar, inplace
void smultiply_inplace_vec(Vector vec, double scalar);

// Perform the dot product of two Vectors
double dot(Vector lhs, Vector rhs);

// Perform the dot product of two Vectors, given that the lhs
// Vector has already been transposed by calling function
double dot_t(Vector lhs, Vector rhs);

// Given a Vector in basis1, return that vector described by basis2
Vector change_ndim_basis(Vector *basis1, Vector *basis2, Vector vec, size_t dim);

// Given the coordinated of a Vector in basis1, return the coordinates of that
// Vector in basis2
Vector change_ndim_basis_coords(Vector *basis1, Vector *basis2, Vector coords1, size_t dim);

// Get the coordinates of the Vector in the given basis, returned as a Vector
Vector get_coords(Vector *basis, size_t dim, Vector vec);

// Get the coordinates of the Vector in the given basis, returned as a double array
double * get_coords_dbls(Vector *basis, size_t dim, Vector vec);

// Return the result of the linear combination of the vectors in basis
// using the coefficients in Vector coords
Vector lin_combination(Vector *basis, size_t dim, Vector coords);

// Return the result of the linear combination of the vectors in basis
// using the coefficients in double array coords
Vector lin_combination_dbls(Vector *basis, size_t dim, double *coords);

// Orthonormalize the set of Vectors using the Gram-Schmidt process
Vector * orthonormalize_set_vec(Vector *vecs, size_t dim);

#endif // VECTOR_H
