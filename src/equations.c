
#include "equations.h"
#include <math.h>
#include <stdio.h>

static Vector make_solution(size_t dim)
{
    if(dim > 0) {
        return make_vector(dim, COL);
    } else {
        return get_null_vec();
    }
}

Vector solve_system(const Matrix mtx)
{
    if(!is_null_mtx(mtx)) {
        Vector sln = solve_invertible_system(mtx);
        if(!is_null_vec(sln)) {
            return sln;
        }
        sln = solve_least_squares(mtx);
        if(!is_null_vec(sln)) {
            return solve_least_squares(mtx);
        }

        Matrix rref_mtx = get_rref(mtx);

        if(is_null_mtx(rref_mtx)) {
            return get_null_vec();
        }
        sln = make_solution(mtx.num_cols - 1);

        for(size_t i = 0; i < rref_mtx.num_rows; ++i) {
            size_t pivot_index = get_pivot_index(rref_mtx.rows[i]);
            if(pivot_index < rref_mtx.num_cols - 1) {
                double value = get_element_vec(rref_mtx.rows[i], rref_mtx.num_cols - 1);
                set_element_vec(value, sln, pivot_index);
                for(size_t j = pivot_index + 1; j < rref_mtx.num_cols - 1; ++j) {
                    value = get_element_vec(rref_mtx.rows[i], j);
                    if(fabs(value) > EPSILON) {
                        set_element_vec(0, sln, j);
                    }
                }
            } else {
                delete_matrix(rref_mtx);
                if(pivot_index == rref_mtx.num_cols - 1) {
                    delete_vector(sln);
                    return get_null_vec();
                } else {
                    return sln;
                }
            }
        }
        delete_matrix(rref_mtx);

        return sln;
    }  else {
        return get_null_vec();
    }
}

Vector solve_system_v(const Matrix mtx, const Vector vec)
{
    if(!is_null_mtx(mtx) && !is_null_vec(vec) && mtx.num_rows == vec.dim && is_col(vec)) {
        Matrix sys_of_equations = hcombine_v(mtx, vec, true);
        Vector sln = solve_system(sys_of_equations);
        delete_matrix(sys_of_equations);

        return sln;
    } else {
        return get_null_vec();
    }
}

Vector solve_invertible_system(const Matrix mtx)
{
    if(!is_null_mtx(mtx) && mtx.num_cols == mtx.num_rows + 1) {
        Matrix sub_mtx = remove_col(mtx, mtx.num_cols - 1);
        if(is_invertible(sub_mtx)) {
            Matrix inverse = get_inverse(sub_mtx);
            Vector last_col = get_col(mtx, mtx.num_cols - 1);

            Vector sln = vmultiply(inverse, last_col, true);

            delete_matrix(inverse);
            delete_matrix(sub_mtx);
            delete_vector(last_col);

            return sln;
        } else {
            delete_matrix(sub_mtx);
            return get_null_vec();
        }
    } else {
        return get_null_vec();
    }
}

Vector solve_invertible_system_v(const Matrix mtx, const Vector vec)
{
    if(!is_null_mtx(mtx) && !is_null_vec(vec) && is_square(mtx) && mtx.num_rows == vec.dim && is_col(vec)) {
        Matrix sys_of_equations = hcombine_v(mtx, vec, true);
        Vector sln = solve_system(sys_of_equations);
        delete_matrix(sys_of_equations);

        return sln;
    } else {
        return get_null_vec();
    }
}

Vector solve_least_squares(const Matrix mtx)
{
    if(!is_null_mtx(mtx)) {
        Matrix sub_mtx = remove_col(mtx, mtx.num_cols - 1);
        Vector sub_vec = get_col(mtx, mtx.num_cols - 1);
        Vector sln = solve_least_squares_v(sub_mtx, sub_vec);

        delete_matrix(sub_mtx);
        delete_vector(sub_vec);

        return sln;
    } else {
        return get_null_vec();
    }
}

Vector solve_least_squares_v(const Matrix mtx, const Vector vec)
{
    if(!is_null_mtx(mtx) && !is_null_vec(vec) && mtx.num_rows == vec.dim && is_col(vec)) {
        Vector sln;
        if(mtx.num_rows >= mtx.num_cols) {
            Matrix qr[2];
            qr_decomposition(mtx, qr);
            print_matrix(mtx);
            printf("%f\n", determinant(mtx));
            print_matrix(qr[0]);
            puts("");
            print_matrix(qr[1]);
            Matrix q_transpose = transpose_mtx(qr[0]);
            Matrix r_inverse = get_inverse(qr[1]);

            Matrix intermediate_mtx = mmultiply(r_inverse, q_transpose);
            sln = vmultiply(intermediate_mtx, vec, true);

            delete_matrix(qr[0]);
            delete_matrix(qr[1]);
            delete_matrix(q_transpose);
            delete_matrix(r_inverse);
            delete_matrix(intermediate_mtx);
        } else {
            sln = get_null_vec();
        }
        return sln;
    } else {
        return get_null_vec();
    }
}

Vector solve_polynomial(Vector coefs)
{
    if(!is_null_vec(coefs)) {
        Matrix comp_mtx = get_companion_matrix(coefs);
        Vector roots = get_eigenvalues(comp_mtx);
        delete_matrix(comp_mtx);

        return roots;
    } else {
        return get_null_vec();
    }
}
