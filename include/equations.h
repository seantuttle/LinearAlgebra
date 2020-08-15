#ifndef EQUATIONS_H
#define EQUATIONS_H

#include "vector.h"
#include "matrix.h"

// Solve a linear system, first trying to do so by inversion and then trying
// least squares, before settling for either choosing an arbitrary solution
// if there are infinitely many or returing NULL_VECTOR if 
// the system is incompatible
Vector solve_system(const Matrix mtx);
Vector solve_system_v(const Matrix mtx, const Vector vec);

// Attempt to solve equation Ax=b by trying to invert A and returning
// A^(-1)b, returning NULL_VECTOR if Matrix is not invertible
Vector solve_invertible_system(const Matrix mtx);
Vector solve_invertible_system_v(const Matrix mtx, const Vector vec);

// Attempt to solve equation A^TAx=A^Tb by back subsitution using the
// QR decomposition of the given Matrix
Vector solve_least_squares(const Matrix mtx);
Vector solve_least_squares_v(const Matrix mtx, const Vector vec);

// Solve the given monic polynomial by finding the eigenvalues
// of its companion Matrix (see get_companion_matrix in matrix.h)
Vector solve_polynomial(Vector coefs);

#endif // EQUATIONS_H
