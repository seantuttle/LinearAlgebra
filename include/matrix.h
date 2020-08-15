// File containing functions to create and manipulate Matrix
// structures. Every Matrix that is created must be accompanied by
// a subsequent call to delete_matrix() in order to ensure proper
// memory management

#ifndef MATRIX_H
#define MATRIX_H

#include "vector.h"
#include <stddef.h>
#include <stdbool.h>

// In this implementation, a Matrix is defined as a group of rows.
// The columns are then determined implicitly from these rows.
// Both columns and rows employ zero based indexing, so the top left
// entry in a Matrix is at row 0 and column 0.
typedef struct Matrix {
    size_t num_rows;
    size_t num_cols;
    Vector *rows;
} Matrix;

// Creates a Matrix with all values initialized to 0
Matrix make_matrix(size_t num_rows, size_t num_cols);

// Provide a deep copy of parameter Matrix
Matrix copy_matrix(Matrix mtx);

// Create a matrix with the given set of rows. All rows must have same dimension
Matrix make_matrix_rows(Vector *rows, size_t num_rows);

// Create a matrix with the given set of columns. All columns must have same dimension
Matrix make_matrix_cols(Vector *cols, size_t num_cols);

// Create a matrix with the given values. A caller should create vals as a double
// array of size[num_rows * num_size], and this function will do the work of
// partitioning the rows and columns as needed
Matrix make_matrix_vals(double *vals, size_t num_rows, size_t num_cols);

// Create a matrix with random values in the range [-max_abs_val, max_abs_val]
Matrix make_matrix_rand(size_t num_rows, size_t num_cols, unsigned int max_abs_val);

// Create an identity Matrix with the given dimension
Matrix make_identity(size_t dim);

// Delete memory associated with a Matrix. Must be done with every created Matrix
void delete_matrix(Matrix mtx);

// Append mtx2 to mtx1 horizontally. They must have the same number of rows
Matrix hcombine(Matrix mtx1, Matrix mtx2);

// Append or prepend vec to mtx horizontally, depending on the mtx_on_left param. 
// They must have the same number of rows.
Matrix hcombine_v(Matrix mtx, Vector vec, bool mtx_on_left);

// Append mtx2 to mtx1 vertically. They must have the same number of rows
Matrix vcombine(Matrix mtx1, Matrix mtx2);

// Append or prepend vec to mtx vertically, depending on the mtx_on_left param. 
// They must have the same number of rows.
Matrix vcombine_v(Matrix mtx, Vector vec, bool mtx_on_top);

// Output Matrix in row-column format to the console screen
void print_matrix(Matrix mtx);

// Output the specified row to the console screen
void print_row(Matrix mtx, size_t row);

// Output the specified column to the console screen
void print_col(Matrix mtx, size_t col);

// Fill Matrix with zeroes
void zeroes_mtx(Matrix mtx);

// Fill Matrix with ones
void ones_mtx(Matrix mtx);

// Fill Matrix with given value
void fill_mtx(Matrix mtx, double value);

// Turn Matrix into identity Matrix, inplace. The Matrix must be square
void identity(Matrix mtx);

// Set element at the specified index of the Matrix to the given value
void set_element_mtx(double value, Matrix mtx, size_t row, size_t col);

// Set the submatrix starting at row_index and col_index to the provided submatrix.
// The provided submatrix must be able to fit in the remaining space of mtx
void set_submatrix_true(Matrix sub_mtx, Matrix mtx, size_t row_index, size_t col_index);

// Set the specified row of the matrix to the parameter row
void set_row(Vector row, Matrix mtx, size_t row_index);

// Set the specified column of the matrix to the parameter column
void set_col(Vector col, Matrix mtx, size_t col_index);

// Set a subsection of the specified row to the parameter subrow
void set_subrow(Vector subrow, Matrix mtx, size_t row_index, size_t replace_index);

// Set a subsection of the specified column to the parameter subcolumn
void set_subcol(Vector subcol, Matrix mtx, size_t col_index, size_t replace_index);

// Swap the two rows in the given Matrix, provided they are both in range
void swap_rows(Matrix mtx, size_t row1_index, size_t row2_index);

// Swap the two columns in the given Matrix, provided they are both in range
void swap_cols(Matrix mtx, size_t col1_index, size_t col2_index);

// Create a Matrix equal to the parameter mtx except for that the row at row_index is deleted
Matrix remove_row(Matrix mtx, size_t row_index);

// Create a Matrix equal to the parameter mtx except for that the column at col_index is deleted
Matrix remove_col(Matrix mtx, size_t col_index);

// Get the element at the specified index
double get_element_mtx(Matrix mtx, size_t row_index, size_t col_index);

// Get a pointer sized (num_rows * num_cols) that contains the values of the Matrix.
// This pointer will need to be freed
double * get_elements_mtx(Matrix mtx);

// Get submatrix created by deleting specified row and column
Matrix get_submatrix_det(Matrix mtx, size_t row_del_index, size_t col_del_index);

// Get submatrix starting at the specified index and ending at the bottom right corner
Matrix get_submatrix_true(Matrix mtx, size_t row_start_index, size_t col_start_index);

// Get the specified row
Vector get_row(Matrix mtx, size_t row_index);

// Get the specified column
Vector get_col(Matrix mtx, size_t col_index);

// Get a subsection of the specified row
Vector get_subrow(Matrix mtx, size_t row_index, size_t start_index);

// Get a subsection of the specified column
Vector get_subcol(Matrix mtx, size_t col_index, size_t start_index);

// Get a pointer containing copies of the rows of the matrix
Vector * get_rows(Matrix mtx);

// Get a pointer containing copies of the columns of the matrix
Vector * get_cols(Matrix mtx);

// Get the companion Matrix of a given monic polynomial that has the form
// x^n+Cn-1*x^n-1+...+C1*x+C0=0. The parameter coefs is a Vector structure
// containing the n coefficients C0...Cn-1, in that order. A companion Matrix
// of a given monic polynomial has the following form:
//      0 0 0 ... -C0
//      1 0 0 ... -C1
//      0 1 0 ... -C2
//          ...
//      0 0 ... 1 -Cn-1
// where the ones are the subdiagonal. The eigenvalues of this Matrix are
// the roots of the given monic polynomial
Matrix get_companion_matrix(Vector coefs);

// Return a singleton NULL_MATRIX that is used akin to NULL ptr
Matrix get_null_mtx(void);

// Check if num_cols == num_rows
bool is_square(Matrix mtx);

// Check if given index is contained within Matrix
bool in_range_mtx(Matrix mtx, size_t row_index, size_t col_index);

// Check if two matrices are the same size
bool same_size_mtx(Matrix mtx1, Matrix mtx2);

// Check if two matrices are equal. Equality implies that they
// are the same size and Aij == Bij for all i, j
bool equal_mtx(Matrix lhs, Matrix rhs);

// Check if lhs.num_cols == rhs.num_rows
bool can_multiply(Matrix lhs, Matrix rhs);

// Check if a matrix is invertible by comparing its determinant to 0
bool is_invertible(Matrix  mtx);

// Check if a Matrix is in row echelon form
bool is_ref(Matrix mtx);

// Check if a Matrix is in reduced row echelon form
bool is_rref(Matrix mtx);

// Check to see if A == A^T
bool is_symmetric(Matrix mtx);

// Check if all entries below main diagonal are zero
bool is_upper_triangular(Matrix mtx);

// Check if all entries above main diagonal are zero
bool is_lower_triangular(Matrix mtx);

// Check if all entries below and above main diagonal are zero
bool is_diagonal(Matrix mtx);

// Check if Matrix is diagonal and all diagonal entries == 1
bool is_identity(Matrix mtx);

// Check if inv_mtx is the inverse of orig_mtx.
// That is to say, check if inv_mtx * orig_mtx == Identity
bool is_inverse(Matrix inv_mtx, Matrix orig_mtx);

// Check if two Matrices are orthogonal (ie <mtx1, mtx2> == 0)
bool are_orthogonal_mtx(Matrix mtx1, Matrix mtx2);

// Check if each Matrix in mtxs1 is orthogonal to each Matrix in mtxs2
bool are_orthogonal_sets_mtx(Matrix *mtxs1, Matrix *mtxs2, size_t dim1, size_t dim2);

// Check if each Matrix in mtxs is orthogonal to every other Matrix in mtxs
bool is_orthogonal_set_mtx(Matrix *mtxs, size_t dim);

// Check if each Matrix in mtxs is orthogonal to every other Matrix in mtxs
// and if each Matrix has norm 1 (ie ||mtx||==sqrt(<mtx, mtx>) == 1)
bool is_orthonormal_set_mtx(Matrix *mtxs, size_t dim);

// Check if Matrix has norm 1 (ie ||mtx||==sqrt(<mtx, mtx>) == 1)
bool is_unit_mtx(Matrix mtx);

// Check if each and every entry in a Matrix is zero
bool is_zero_mtx(Matrix mtx);

// Check if A^2 == A, the definition of idempotence
bool is_idempotent(Matrix mtx);

// Check if A^2 == Identity, the definition of involution
bool is_involution(Matrix mtx);

// Check if Matrix is equal to the singleton NULL_MATRIX
bool is_null_mtx(Matrix mtx);

// Return the transpose of the given Matrix
Matrix transpose_mtx(Matrix mtx);

// Perform element-wise addition of two Matrices and return the result
Matrix add_mtx(Matrix lhs, Matrix rhs);

// Perform element-wise subtraction of two Matrices and return the result
Matrix subtract_mtx(Matrix lhs, Matrix rhs);

// Multiply each element in a Matrix by a scalar and return result in new Matrix
Matrix smultiply_mtx(Matrix mtx, double scalar);

// Multiply each element in a Matrix by a scalar, inplace
void smultiply_inplace_mtx(Matrix mtx, double scalar);

// Perform standard Matrix multiplication
Matrix mmultiply(Matrix lhs, Matrix rhs);

// Multiply a Matrix by a Vector, with the order determined by mtx_on_left param
Vector vmultiply(Matrix mtx, Vector vec, bool mtx_on_left);

// Raise the given Matrix the power defined by exponent
Matrix expon(Matrix mtx, unsigned int exponent);

// Find the determinant of the given Matrix (must be square)
double determinant(Matrix mtx);

// Get the row echelon form of given Matrix returned in a new Matrix
Matrix get_ref(Matrix mtx);

// Get the row echelon form of given Matrix inplace
void get_ref_inplace(Matrix mtx);

// Reduce mtx1 to its row echelon form while applying the same row operations
// to mtx2. The result of performing these operations on mtx2 is returned
Matrix get_ref_mirror(Matrix mtx1, Matrix mtx2);

// Get the reduced row echelon form of given Matrix returned in a new Matrix
Matrix get_rref(Matrix mtx);

// Get the reduced row echelon form of given Matrix inplace
void get_rref_inplace(Matrix mtx);

// Reduce mtx1 to its reduced row echelon form while applying the same row operations
// to mtx2. The result of performing these operations on mtx2 is returned
Matrix get_rref_mirror(Matrix mtx1, Matrix mtx2);

// Reduce the given Matrix to a similar Matrix of upper Hessenberg form.
// Upper Hessenberg form means that all entries below the subdiagonal are zero.
Matrix get_hessenberg(Matrix mtx);

// Reduce the given Matrix to a similar Matrix of upper Hessenberg form, inplace.
// Upper Hessenberg form means that all entries below the subdiagonal are zero.
void get_hessenberg_inplace(Matrix mtx);

// Get the standard basis for R^(dimXdim). For dim == 2, this would be:
//      1 0     0 1     0 0     0 0
//      0 0     0 0     1 0     0 1
Matrix * get_nndim_basis(size_t dim);

// Get the basis of the null space of the given Matrix
Vector * get_nullspace_basis(const Matrix mtx);

// Get the basis of the row space of the given Matrix
Vector * get_rowspace_basis(const Matrix mtx);

// Get the basis of the column space of the given Matrix
Vector * get_colspace_basis(const Matrix mtx);

// Get the basis of the column complament space of the given Matrix.
// By Fundamental Subspaces, Range(A)^C == Null(A^T)
Vector * get_colcomp_basis(const Matrix mtx);

// Get the basis of the null complament space of the given Matrix.
// By Fundamental Subspaces, Null(A)^C == Range(A^T)
Vector * get_nullcomp_basis(const Matrix mtx);

// Get the trace of the given Matrix. The trace is the sum of
// all of the diagonal entries
double get_trace(Matrix mtx);

// Get the rank of the Matrix. Rank is the dimension of a matrix's row
// and column spaces
size_t get_rank(Matrix mtx);

// Get the nullity of the Matrix. Nullity is the dimension of a matrix's
// null space. Rank-Nullity theorem states that Rank+Nullity == num_cols
size_t get_nullity(Matrix mtx);

// Determine the inverse of the Matrix, if one exists
Matrix get_inverse(Matrix mtx);

// Get the standard inner product for R^(mXn), written <mtx1, mtx2>
double get_inner_product(Matrix mtx1, Matrix mtx2);

// Get the outer product of two vectors, which produces a Matrix
Matrix get_outer_product(Vector vec1, Vector vec2);

// Get the Euclidean norm of a Matrix. This norm is defined as sqrt(<mtx, mtx>)
double get_norm_mtx(Matrix mtx);

// Get the angle between two matrices. This is done using cosine, which in an inner
// product space is defined as: cos(theta) = <mtx1, mtx2> / (||mtx1|| * ||mtx2||)
// Angle is returned in degress.
double get_angle_mtx(Matrix mtx1, Matrix mtx2);

// Get the angle between two matrices. This is done using cosine, which in an inner
// product space is defined as: cos(theta) = <mtx1, mtx2> / (||mtx1|| * ||mtx2||)
// Angle is returned in radians.
double get_angle_rads_mtx(Matrix mtx1, Matrix mtx2);

// Get the distance between two Matrices. In an inner product space, distance is
// defined as ||mtx2 - mtx1||
double get_distance_mtx(Matrix mtx1, Matrix mtx2);

// Get the unit Matrix pointing the direction of the given Matrix. This is done
// by dividing the Matrix by its norm
Matrix get_unit_mtx(Matrix mtx);

// Get the scalar projection of mtx1 onto mtx2
double get_scalar_proj_mtx(Matrix mtx1, Matrix mtx2);

// Get the vector projection of mtx1 onto mtx2
Matrix get_vector_proj_mtx(Matrix mtx1, Matrix mtx2);

// Get the projection of mtx onto the span of basis, which is 
// the Matrix in basis that is closest to mtx (min distance)
Matrix get_projection_mtx(Matrix *basis, size_t dim, Matrix mtx);

// Orthonormalize a set of Matrices using the Gram-Schmidt process
Matrix * orthonormalize_set_mtx(Matrix *mtxs, size_t dim);

// Decompose the given Matrix into the multiplication of two Matrices
// Q and R, where Q has orthonormal columns and R is upper triangular.
// This is donoe using Householder transformations
void qr_decomposition(Matrix mtx, Matrix *qr);

// Decompose the given Matrix into the multiplication of two Matrices
// R and Q, where Q has orthonormal columns and R is upper triangular.
// This is done using the Gram-Schmidt process
void rq_decomposition(Matrix mtx, Matrix *rq);

// Get the eigenvalues of the given Matrix. This is done by first reducing the
// matrix to its upper Hessenberg form, which is similar (meaning they have the same
// eigenvalues). Then a single (Wilkinson) shifted QR algorithm with deflation is
// applied to this Hessenberg matrix, yielding the eigenvalues of the original Matrix.
// Note that this library was written without support for complex arithmetic, so that
// limits the usefullness of eigenvalues and their related vectors and spaces
Vector get_eigenvalues(Matrix mtx);

// Get the basis of the eigenspace corresponding to the given Matrix and eigenvalue.
// That is to say, Null(mtx - eigenvalue * Identity)
Vector * get_eigenspace_basis(Matrix mtx, double eigenvalue);

// Get the dimension of the corresponding eigenspace
size_t get_eigenspace_basis_dim(Matrix mtx, double eigenvalue);

// Get an eigenvectors corresponding to the given eigenvalues. Note that there
// are infinitely many eigenvectors to choose from, so an arbitrary one is returned
Vector get_eigenvector(Matrix mtx, double eigenvalue);

// Get the set of eigenvectors that correspond to the eigenvalues
Vector * get_eigenvectors(Matrix mtx, Vector eigenvalues);

// Decompose a Matrix into the multiplication of three Matrixes X, D, and X^-1,
// where D is diagonal and the values along the diagonal are the eigenvalues of X
// and X is made up of the corresponding eigenvectors as its columns
void diagonalize(Matrix mtx, Matrix *xdx_1);

#endif // MATRIX_H
