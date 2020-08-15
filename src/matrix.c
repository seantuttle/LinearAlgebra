
#include "matrix.h"
#include "equations.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

const Matrix NULL_MATRIX = {0, 0, NULL}; // singleton value for Matrix NULL

Matrix make_matrix(size_t num_rows, size_t num_cols)
{
    if(num_rows > 0 && num_cols > 0) {
        Vector size_check = make_vector(num_cols, ROW);
        Matrix new_mtx = {num_rows, num_cols, malloc(num_rows * sizeof(size_check))};
        if(!new_mtx.rows) {
            return get_null_mtx();
        }

        for(size_t i = 0; i < num_rows; ++i) {
            new_mtx.rows[i] = make_vector(num_cols, ROW);
        }
        zeroes_mtx(new_mtx);
        delete_vector(size_check);
        
        return new_mtx;
    } else {
        return get_null_mtx();
    }
}

Matrix copy_matrix(Matrix mtx)
{
    Matrix new_mtx = make_matrix(mtx.num_rows, mtx.num_cols);

    for(size_t i = 0; i < new_mtx.num_rows; ++i) {
        Vector row = copy_vector(mtx.rows[i]);
        set_row(row, new_mtx, i);
    }

    return new_mtx;
}

Matrix make_matrix_rows(Vector *rows, size_t num_rows)
{
    if(num_rows > 0 && rows) {
        size_t init_num_cols = rows[0].dim;
        Matrix mtx = make_matrix(num_rows, init_num_cols);
        for(size_t i = 0; i < num_rows; ++i) {
            if(rows[i].dim == init_num_cols) {
                set_row(rows[i], mtx, i);
            } else {
                delete_matrix(mtx);
                return get_null_mtx();
            }
        }
        return mtx;
    } else {
        return get_null_mtx();
    }
}

Matrix make_matrix_cols(Vector *cols, size_t num_cols)
{
    if(num_cols > 0 && cols) {
        size_t init_num_rows = cols[0].dim;
        Matrix mtx = make_matrix(init_num_rows, num_cols);
        for(size_t i = 0; i < num_cols; ++i) {
            if(cols[i].dim == init_num_rows) {
                set_col(cols[i], mtx, i);
            } else {
                delete_matrix(mtx);
                return get_null_mtx();
            }
        }
        return mtx;
    } else {
        return get_null_mtx();
    }
}

Matrix make_matrix_vals(double *vals, size_t num_rows, size_t num_cols)
{
    if(num_rows > 0 && num_cols > 0 && vals) {
        size_t vals_ctr = 0;
        Matrix mtx = make_matrix(num_rows, num_cols);
        for(size_t i = 0; i < mtx.num_rows; ++i) {
            for(size_t j = 0; j < mtx.num_cols; ++j) {
                set_element_mtx(vals[vals_ctr], mtx, i, j);
                ++vals_ctr;
            }
        }

        return mtx;
    } else {
        return get_null_mtx();
    }
}

Matrix make_matrix_rand(size_t num_rows, size_t num_cols, unsigned int max_abs_val)
{
    Matrix mtx = make_matrix(num_rows, num_cols);
    if(!is_null_mtx(mtx)) {
        srand((unsigned) time(NULL));
        int max_val_standin = (2 * max_abs_val) + 1;
        for(size_t i = 0; i < mtx.num_rows; ++i) {
            for(size_t j = 0; j < mtx.num_cols; ++j) {
                double value = rand() % max_val_standin;
                set_element_mtx(value - max_abs_val, mtx, i, j);
            }
        }
    }
    return mtx;
}

Matrix make_identity(size_t dim)
{
    if(dim > 0) {
        Matrix ident = make_matrix(dim, dim);
        identity(ident);

        return ident;
    } else {
        return get_null_mtx();
    }
}

void delete_matrix(Matrix mtx)
{
    for(size_t i = 0; i < mtx.num_rows; ++i) {
        delete_vector(mtx.rows[i]);
    }
    free(mtx.rows);
}

Matrix hcombine(Matrix mtx1, Matrix mtx2)
{
    if(!is_null_mtx(mtx1) && !is_null_mtx(mtx2) && mtx1.num_rows == mtx2.num_rows) {
        Matrix combo = make_matrix(mtx1.num_rows, mtx1.num_cols + mtx2.num_cols);
        set_submatrix_true(mtx1, combo, 0, 0);
        set_submatrix_true(mtx2, combo, 0, mtx1.num_cols);

        return combo;
    } else {
        return get_null_mtx();
    }
}

Matrix hcombine_v(Matrix mtx, Vector vec, bool mtx_on_left)
{
    if(is_col(vec) && vec.dim == mtx.num_rows) {
        Matrix combo = make_matrix(mtx.num_rows, mtx.num_cols + 1);
        Matrix vec_mtx = make_matrix_cols(&vec, 1);
        if(mtx_on_left) {
            set_submatrix_true(mtx, combo, 0, 0);
            set_submatrix_true(vec_mtx, combo, 0, mtx.num_cols);
        } else {
            set_submatrix_true(vec_mtx, combo, 0, 0);
            set_submatrix_true(mtx, combo, 0, 1);
        }
        delete_matrix(vec_mtx);
        return combo;
    } else {
        return get_null_mtx();
    }
}

Matrix vcombine(Matrix mtx1, Matrix mtx2)
{
    if(!is_null_mtx(mtx1) && !is_null_mtx(mtx2) && mtx1.num_cols == mtx2.num_cols) {
        Matrix combo = make_matrix(mtx1.num_rows + mtx2.num_rows, mtx1.num_cols);
        set_submatrix_true(mtx1, combo, 0, 0);
        set_submatrix_true(mtx2, combo, mtx1.num_rows, 0);

        return combo;
    } else {
        return get_null_mtx();
    }
}

Matrix vcombine_v(Matrix mtx, Vector vec, bool mtx_on_top)
{
   if(is_row(vec) && vec.dim == mtx.num_cols) {
        Matrix combo = make_matrix(mtx.num_rows, mtx.num_cols + 1);
        Matrix vec_mtx = make_matrix_rows(&vec, 1);
        if(mtx_on_top) {
            set_submatrix_true(mtx, combo, 0, 0);
            set_submatrix_true(vec_mtx, combo, mtx.num_rows, 1);
        } else {
            set_submatrix_true(vec_mtx, combo, 0, 0);
            set_submatrix_true(mtx, combo, 1, 0);
        }
        delete_matrix(vec_mtx);
        return combo;
    } else {
        return get_null_mtx();
    }
}

void print_matrix(Matrix mtx) 
{
    for(size_t i = 0; i < mtx.num_rows; ++i) {
        print_vector(mtx.rows[i]);
    }
}

void print_row(Matrix mtx, size_t row_index)
{
    print_vector(mtx.rows[row_index]);
}

void print_col(Matrix mtx, size_t col_index)
{
    Vector col = get_col(mtx, col_index);
    print_vector(col);
    delete_vector(col);
}

void zeroes_mtx(Matrix mtx)
{
    fill_mtx(mtx, 0);
}

 void ones_mtx(Matrix mtx)
 {
     fill_mtx(mtx, 1);
 }

void fill_mtx(Matrix mtx, double value)
{
    for(size_t i = 0; i < mtx.num_rows; ++i){
        fill_vec(mtx.rows[i], value);
    }  
}

void identity(Matrix mtx)
{
    if(is_square(mtx)) {
        for(size_t i = 0; i < mtx.num_rows; ++i) {
            for(size_t j = 0; j < mtx.num_rows; ++j) {
                if(i == j) {
                    set_element_mtx(1, mtx, i, j);
                } else {
                    set_element_mtx(0, mtx, i, j);
                }
            }
        }
    }
}

void set_element_mtx(double value, Matrix mtx, size_t row_index, size_t col_index)
{
    if(in_range_mtx(mtx, row_index, col_index)) {
        set_element_vec(value, mtx.rows[row_index], col_index);
    }
}

void set_submatrix_true(Matrix sub_mtx, Matrix mtx, size_t row_index, size_t col_index)
{
    if(row_index + sub_mtx.num_rows <= mtx.num_rows && col_index + sub_mtx.num_cols <= mtx.num_cols) {
        for(size_t i = 0; i < sub_mtx.num_rows; ++i) {
            for(size_t j = 0; j < sub_mtx.num_cols; ++j) {
                double value = get_element_mtx(sub_mtx, i, j);
                set_element_mtx(value, mtx, i + row_index, j + col_index);
            }
        }
    }
}

void set_row(Vector row, Matrix mtx, size_t row_index)
{
    if(row_index < mtx.num_rows) {
        delete_vector(mtx.rows[row_index]);
        mtx.rows[row_index] = copy_vector(row);
    }
}

void set_col(Vector col, Matrix mtx, size_t col_index)
{
    if(col_index < mtx.num_cols) {
        for(size_t i = 0; i < mtx.num_rows; ++i) {
            double value = get_element_vec(col, i);
            set_element_mtx(value, mtx, i, col_index);
        }
    }
}

void set_subrow(Vector subrow, Matrix mtx, size_t row_index, size_t replace_index)
{
    if(row_index < mtx.num_rows && subrow.dim + replace_index <= mtx.num_cols) {
        for(size_t i = 0; i < subrow.dim; ++i) {
            double value = get_element_vec(subrow, i);
            set_element_mtx(value, mtx, row_index, replace_index + i);
        }
    } 
}

void set_subcol(Vector subcol, Matrix mtx, size_t col_index, size_t replace_index)
{
    if(col_index < mtx.num_cols && subcol.dim + replace_index <= mtx.num_rows) {
        for(size_t i = 0; i < subcol.dim; ++i) {
            double value = get_element_vec(subcol, i);
            set_element_mtx(value, mtx, col_index, replace_index + i);
        }
    } 
}

void swap_rows(Matrix mtx, size_t row1_index, size_t row2_index)
{
    if(row1_index != row2_index && row1_index < mtx.num_rows && row2_index < mtx.num_rows) {
        Vector temp_row = get_row(mtx, row1_index);
        set_row(mtx.rows[row2_index], mtx, row1_index);
        set_row(temp_row, mtx, row2_index);
    }
}

void swap_cols(Matrix mtx, size_t col1_index, size_t col2_index)
{
    if(col1_index != col2_index && col1_index < mtx.num_cols && col2_index < mtx.num_cols) {
        Vector col1 = get_col(mtx, col1_index);
        Vector col2 = get_col(mtx, col2_index);

        set_col(col2, mtx, col1_index);
        set_col(col1, mtx, col2_index);

        delete_vector(col1);
        delete_vector(col2);
    } 
}

Matrix remove_row(Matrix mtx, size_t row_index)
{
    if(!is_null_mtx(mtx) && row_index < mtx.num_rows) {
        Matrix new_mtx = make_matrix(mtx.num_rows - 1, mtx.num_cols);
        size_t new_mtx_ctr = 0;
        for(size_t i = 0; i < mtx.num_rows; ++i) {
            if(i != row_index) {
                set_row(mtx.rows[i], new_mtx, new_mtx_ctr);
                ++new_mtx_ctr;
            }
        }
        return new_mtx;
    } else {
        return get_null_mtx();
    }
}

Matrix remove_col(Matrix mtx, size_t col_index)
{
    if(!is_null_mtx(mtx) && col_index < mtx.num_cols) {
        Matrix new_mtx = make_matrix(mtx.num_rows, mtx.num_cols - 1);
        size_t new_mtx_ctr = 0;
        for(size_t i = 0; i < mtx.num_cols; ++i) {
            if(i != col_index) {
                Vector col = get_col(mtx, i);
                set_col(col, new_mtx, new_mtx_ctr);
                delete_vector(col);
                ++new_mtx_ctr;
            }
        }
        return new_mtx;
    } else {
        return get_null_mtx();
    }
}

double get_element_mtx(Matrix mtx, size_t row_index, size_t col_index)
{
    return get_element_vec(mtx.rows[row_index], col_index);
}

double * get_elements_mtx(Matrix mtx)
{
    double *elements = malloc(mtx.num_rows * mtx.num_cols * sizeof *elements);
    size_t elements_ctr = 0;
    if(elements) {
        for(size_t i = 0; i < mtx.num_rows; ++i) {
            for(size_t j = 0; j < mtx.num_cols; ++j) {
                elements[elements_ctr] = get_element_mtx(mtx, i, j);
                ++elements_ctr;
            }
        }
    }
    return elements;
}

Matrix get_submatrix_det(Matrix mtx, size_t row_del_index, size_t col_del_index)
{
    if(mtx.num_rows > 1 && mtx.num_rows > 1) {
        Matrix sub_mtx = make_matrix(mtx.num_rows - 1, mtx.num_cols - 1);
        size_t sub_mtx_row_i = 0;
        size_t sub_mtx_col_i = 0;

        for(size_t i = 0; i < mtx.num_rows; ++i) {
            if(i != row_del_index) {
                for(size_t j = 0; j < mtx.num_cols; ++j) {
                    if(j != col_del_index) {
                        double value = get_element_mtx(mtx, i, j);
                        set_element_mtx(value, sub_mtx, sub_mtx_row_i, sub_mtx_col_i);
                        ++sub_mtx_col_i;
                    }
                }
                sub_mtx_col_i = 0;
                ++sub_mtx_row_i;
            }
        }
        return sub_mtx;
    } else {
        return get_null_mtx();
    }
}

Matrix get_submatrix_true(Matrix mtx, size_t row_start_index, size_t col_start_index)
{
    if(!is_null_mtx(mtx) && in_range_mtx(mtx, row_start_index, col_start_index)) {
        size_t sub_num_rows = mtx.num_rows - row_start_index;
        size_t sub_num_cols = mtx.num_cols - col_start_index;
        Matrix sub_mtx = make_matrix(sub_num_rows, sub_num_cols);

        for(size_t i = 0; i < sub_mtx.num_rows; ++i) {
            for(size_t j = 0; j < sub_mtx.num_cols; ++j) {
                size_t val_row_index = i + row_start_index;
                size_t val_col_index = j + col_start_index;
                double value = get_element_mtx(mtx, val_row_index, val_col_index);
                set_element_mtx(value, sub_mtx, i, j);
            }
        }
        return sub_mtx;
    } else {
        return get_null_mtx();
    }
}

Vector get_row(Matrix mtx, size_t row_index)
{
    if(row_index < mtx.num_rows) {
        Vector row = copy_vector(mtx.rows[row_index]);
        return row;
    } else {
        return get_null_vec();
    }
}

Vector get_col(Matrix mtx, size_t col_index)
{
    if(col_index < mtx.num_cols) {
        Vector col = make_vector(mtx.num_rows, COL);
        for(size_t i = 0; i < mtx.num_rows; ++i) {
            double value = get_element_mtx(mtx, i, col_index);
            set_element_vec(value, col, i);
        }

        return col;
    } else {
        return get_null_vec();
    }
}

Vector get_subrow(Matrix mtx, size_t row_index, size_t start_index)
{
    if(!is_null_mtx(mtx) && row_index < mtx.num_rows && start_index < mtx.num_cols) {
        return get_subvector(mtx.rows[row_index], start_index);
    } else {
        return get_null_vec();
    }
}

Vector get_subcol(Matrix mtx, size_t col_index, size_t start_index)
{
    if(!is_null_mtx(mtx) && col_index < mtx.num_cols && start_index < mtx.num_rows) {
        Vector col = get_col(mtx, col_index);
        Vector sub_col = get_subvector(col, start_index);
        delete_vector(col);

        return sub_col;
    } else {
        return get_null_vec();
    }
}

Vector * get_rows(Matrix mtx)
{
    if(!is_null_mtx(mtx)) {
        Vector *row_vecs = malloc(mtx.num_rows * sizeof(mtx.rows[0]));
        if(row_vecs) {
            for(size_t i = 0; i < mtx.num_rows; ++i) {
                row_vecs[i] = get_row(mtx, i);
            }
        }
        return row_vecs;
    } else {
        return NULL;
    }
}

Vector * get_cols(Matrix mtx)
{
    if(!is_null_mtx(mtx)) {
        Vector reference = make_vector(mtx.num_rows, COL);
        Vector *col_vecs = malloc(mtx.num_cols * sizeof(reference));
        delete_vector(reference);
        if(col_vecs) {
            for(size_t i = 0; i < mtx.num_cols; ++i) {
                col_vecs[i] = get_col(mtx, i);
            }
        }
        return col_vecs;
    } else {
        return NULL;
    }
}

Matrix get_companion_matrix(Vector coefs)
{
    if(!is_null_vec(coefs)) {
        Matrix comp_mtx = make_matrix(coefs.dim, coefs.dim);
        for(size_t i = 0; i < comp_mtx.num_rows; ++i) {
            for(size_t j = 0; j < comp_mtx.num_cols; ++j) {
                if(j == comp_mtx.num_cols - 1) {
                    double value = get_element_vec(coefs, i);
                    set_element_mtx(-1 * value, comp_mtx, i, j);
                } else if(j == i - 1) {
                    set_element_mtx(1, comp_mtx, i, j);
                }
            }
        }
        return comp_mtx;
    } else {
        return get_null_mtx();
    }
}

Matrix get_null_mtx(void)
{
    return NULL_MATRIX;
}

bool is_square(Matrix mtx)
{
    return mtx.num_cols == mtx.num_rows && !is_null_mtx(mtx);
}

bool in_range_mtx(Matrix mtx, size_t row, size_t col)
{
    return (row < mtx.num_rows) && (col < mtx.num_cols);
}

bool same_size_mtx(Matrix mtx1, Matrix mtx2)
{
    return (mtx1.num_rows == mtx2.num_rows) && (mtx1.num_cols == mtx2.num_cols);
}

bool equal_mtx(Matrix lhs, Matrix rhs)
{
    if(same_size_mtx(lhs, rhs)) {
        for(size_t i = 0; i < lhs.num_rows; ++i) {
            if(!equal_vec(lhs.rows[i], rhs.rows[i])) {
                return false;
            }
        }
        return true;
    } else {
        return false;
    }
}

bool can_multiply(Matrix lhs, Matrix rhs)
{
    return lhs.num_cols == rhs.num_rows;
}

bool is_invertible(Matrix  mtx)
{
    if(is_square(mtx) && !is_null_mtx(mtx)) {
        double det = determinant(mtx);

        return fabs(det) > EPSILON;
    } else {
        return false;
    }
}

bool is_ref(Matrix mtx)
{
    if(is_zero_mtx(mtx)) {
        return true;
    }
    for(size_t i = 0; i < mtx.num_rows; ++i) {
        size_t row_piv_index = get_pivot_index(mtx.rows[i]);
        if(row_piv_index < mtx.num_cols) {
            if(fabs(get_element_mtx(mtx, i, row_piv_index) - 1) < EPSILON) {
                for(size_t j = i + 1; j < mtx.num_rows; ++j) {
                    for(size_t k = 0; k <= row_piv_index; ++k) {
                        if(fabs(get_element_mtx(mtx, j, k)) > EPSILON) {
                            return false;
                        }
                    }
                }
            } else {
                return false;
            }
        } else {
            for(size_t j = i + 1; j < mtx.num_rows; ++j) {
                if(!is_zero_vec(mtx.rows[j])) {
                    return false;
                }
            }
            return true;
        }
    }
    return true;
}

bool is_rref(Matrix mtx)
{
    if(is_zero_mtx(mtx)) {
        return true;
    }
    if(is_ref(mtx)) {
        for(size_t i = 0; i < mtx.num_rows; ++i) {
            size_t row_piv_index = get_pivot_index(mtx.rows[i]);
            if(row_piv_index < mtx.num_cols) {
                Vector col = get_col(mtx, row_piv_index);
                for(size_t j = 0; j < i; ++j) {
                    if(fabs(get_element_vec(col, j)) > EPSILON) {
                        return false;
                    }
                }
                delete_vector(col);
            } else {
                return true;
            }
        }
        return true;
    } else {
        return false;
    }
}

bool is_symmetric(Matrix mtx)
{
    if(!is_null_mtx(mtx)) {
        Matrix t_mtx = transpose_mtx(mtx);
        bool is_symm = equal_mtx(mtx, t_mtx);
        delete_matrix(t_mtx);

        return is_symm;
    } else {
        return false;
    }
}

static bool is_triangular(Matrix mtx, bool (*requirement_func)(size_t, size_t, double))
{
    if(is_square(mtx)) {
        for(size_t i = 0; i < mtx.num_rows; ++i) {
            for(size_t j = 0; j < mtx.num_cols; ++j) {
                double value = get_element_mtx(mtx, i, j);
                if(requirement_func(i, j, fabs(value))) {
                    return false;
                }
            }
        }
        return true;
    } else {
        return false;
    }
}

static bool upper_triangular_requirement(size_t row_i, size_t col_i, double abs_val)
{
    return row_i > col_i && abs_val > EPSILON;
}

bool is_upper_triangular(Matrix mtx)
{
    return is_triangular(mtx, upper_triangular_requirement);
}

static bool lower_triangular_requirement(size_t row_i, size_t col_i, double abs_val)
{
    return row_i < col_i && abs_val > EPSILON;
}

bool is_lower_triangular(Matrix mtx)
{
    return is_triangular(mtx, lower_triangular_requirement);
}

static bool diagonal_check(Matrix mtx, bool (*requirement_func)(double))
{
    if(is_square(mtx)) {
        for(size_t i = 0; i < mtx.num_rows; ++i) {
            for(size_t j = 0; j < mtx.num_cols; ++j) {
                double value = get_element_mtx(mtx, i, j);
                if((i == j && requirement_func(value)) || (i != j && fabs(value) > EPSILON)) {
                    return false;
                }
            }
        }
        return true;
    } else {
        return false;
    }
}

static bool is_diagonal_requirement(double value)
{
    return fabs(value) < EPSILON;
}

bool is_diagonal(Matrix mtx)
{
    return diagonal_check(mtx, is_diagonal_requirement);
}

static bool is_identity_requirement(double value)
{
    return fabs(value - 1) > EPSILON;
}

bool is_identity(Matrix mtx)
{
    return diagonal_check(mtx, is_identity_requirement);
}

bool is_inverse(Matrix inv_mtx, Matrix orig_mtx)
{
    if(is_square(orig_mtx) && same_size_mtx(orig_mtx, inv_mtx)) {
        Matrix product = mmultiply(orig_mtx, inv_mtx);
        bool is_inv = is_identity(product);
        delete_matrix(product);

        return is_inv;
    } else {
        return false;
    }
}

bool are_orthogonal_mtx(Matrix mtx1, Matrix mtx2)
{
    if(same_size_mtx(mtx1, mtx2) && !is_null_mtx(mtx1)) {
        return fabs(get_inner_product(mtx1, mtx2)) < EPSILON;
    } else {
        return false;
    }
}

bool are_orthogonal_sets_mtx(Matrix *mtxs1, Matrix *mtxs2, size_t dim1, size_t dim2)
{
    if(mtxs1 && mtxs2 && dim1 > 0 && dim2 > 0) {
        for(size_t i = 0; i < dim1; ++i) {
            for(size_t j = 0; j < dim2; ++j) {
                if(!are_orthogonal_mtx(mtxs1[i], mtxs2[j])) {
                    return false;
                }
            }
        }
        return true;        
    } else {
        return false;
    }
}

static bool check_orthonormal_or_orthogonal(Matrix *mtxs, size_t dim, bool (*func)(Matrix, Matrix))
{
    if(mtxs && dim > 0) {
        for(size_t i = 0; i < dim; ++i) {
            for(size_t j = 0; j < dim; ++j) {
                if(i != j && !func(mtxs[i], mtxs[j])) {
                    return false;
                }
            }
        }
        return true;
    } else {
        return false;
    }
}

bool is_orthogonal_set_mtx(Matrix *mtxs, size_t dim)
{
    return check_orthonormal_or_orthogonal(mtxs, dim, are_orthogonal_mtx);
}

static bool orthonormal_requirement(Matrix mtx1, Matrix mtx2)
{
    return are_orthogonal_mtx(mtx1, mtx2) && is_unit_mtx(mtx1) && is_unit_mtx(mtx2);
}

bool is_orthonormal_set_mtx(Matrix *mtxs, size_t dim)
{
    return check_orthonormal_or_orthogonal(mtxs, dim, orthonormal_requirement);
}

bool is_unit_mtx(Matrix mtx)
{
    if(!is_null_mtx(mtx)) {
        return fabs(get_norm_mtx(mtx) - 1) < EPSILON;
    } else {
        return false;
    }
}

bool is_zero_mtx(Matrix mtx)
{
    for(size_t i = 0; i < mtx.num_rows; ++i) {
        if(!is_zero_vec(mtx.rows[i])) {
            return false;
        }
    }
    return true;
}

bool is_idempotent(Matrix mtx)
{
    if(!is_null_mtx(mtx) && is_square(mtx)) {
        Matrix product = expon(mtx, 2);
        bool is_idem = equal_mtx(mtx, product);
        delete_matrix(product);

        return is_idem;
    } else {
        return false;
    }
}

bool is_involution(Matrix mtx)
{
    if(!is_null_mtx(mtx) && is_square(mtx)) {
        Matrix product = expon(mtx, 2);
        bool is_inv = is_identity(mtx);
        delete_matrix(product);

        return is_inv;
    } else {
        return false;
    }   
}

bool is_null_mtx(Matrix mtx)
{
    return equal_mtx(mtx, get_null_mtx());
}

Matrix transpose_mtx(Matrix mtx)
{
    Matrix t_mtx = make_matrix(mtx.num_cols, mtx.num_rows);

    for(size_t i = 0; i < mtx.num_rows; ++i) {
        for(size_t j = 0; j < mtx.num_cols; ++j) {
            double value = get_element_mtx(mtx, i, j);
            set_element_mtx(value, t_mtx, j, i);
        }
    }

    return t_mtx;
}

Matrix add_mtx(Matrix lhs, Matrix rhs)
{
    if(same_size_mtx(lhs, rhs)) {
        Matrix sum = make_matrix(lhs.num_rows, lhs.num_cols);

        for(size_t i = 0; i < lhs.num_rows; ++i) {
            Vector sum_row = add_vec(lhs.rows[i], rhs.rows[i]);
            set_row(sum_row, sum, i);
        }
        return sum;
    } else {
        return get_null_mtx();
    }
}

Matrix subtract_mtx(Matrix lhs, Matrix rhs)
{
    Matrix neg_rhs = smultiply_mtx(rhs, -1);
    Matrix difference = add_mtx(lhs, neg_rhs);

    delete_matrix(neg_rhs);

    return difference;
}

Matrix smultiply_mtx(Matrix mtx, double scalar)
{
    Matrix product = copy_matrix(mtx);

    smultiply_inplace_mtx(product, scalar);

    return product;
}

void smultiply_inplace_mtx(Matrix mtx, double scalar)
{
    for(size_t i = 0; i < mtx.num_rows; ++i) {
        smultiply_inplace_vec(mtx.rows[i], scalar);
    }
}

Matrix mmultiply(Matrix lhs, Matrix rhs)
{
    if(can_multiply(lhs, rhs)) {
        Matrix product = make_matrix(lhs.num_rows, rhs.num_cols);

        for(size_t i = 0; i < product.num_rows; ++i) {
            for(size_t j = 0; j < product.num_cols; ++j) {
                Vector lhs_row = get_row(lhs, i);
                Vector rhs_col = get_col(rhs, j);
                double value = dot_t(lhs_row, rhs_col);

                delete_vector(lhs_row);
                delete_vector(rhs_col);

                if(!isnan(value)) {
                    set_element_mtx(value, product, i, j);
                } else {
                    delete_matrix(product);
                    return get_null_mtx();
                }
            }
        }
        return product;
    } else {
        return get_null_mtx();
    }
}

Vector vmultiply(Matrix mtx, Vector vec, bool mtx_on_left)
{
    Matrix vec_mtx;
    if(vec.type == ROW) {
        vec_mtx = make_matrix_rows(&vec, 1);
    } else {
        vec_mtx = make_matrix_cols(&vec, 1);
    }
    
    Vector result_vec;
    Matrix result_mtx;
    if(mtx_on_left) {
        if(can_multiply(mtx, vec_mtx)) {
            result_mtx = mmultiply(mtx, vec_mtx);
        } else {
            delete_matrix(vec_mtx);
            return get_null_vec();
        }
    } else {
        if(can_multiply(vec_mtx, mtx)) {
            result_mtx = mmultiply(vec_mtx, mtx);
        } else {
            delete_matrix(vec_mtx);
            return get_null_vec();
        }
    }
    if(result_mtx.num_rows < result_mtx.num_cols) {
        result_vec = get_row(result_mtx, 0);
    } else {
        result_vec = get_col(result_mtx, 0);
    }
    delete_matrix(result_mtx);
    delete_matrix(vec_mtx);

    return result_vec;
}

Matrix expon(Matrix mtx, unsigned int exponent)
{
    if(is_square(mtx) && !is_null_mtx(mtx)) {
        Matrix power = copy_matrix(mtx);
        
        if(!is_diagonal(power)) {
            for(size_t i = 1; i < exponent; ++i) {
                power = mmultiply(power, mtx);
            }
        } else {
            for(size_t i = 0; i < power.num_rows; ++i) {
                double old_value = get_element_mtx(power, i, i);
                double new_value = pow(old_value, exponent);
                set_element_mtx(new_value, power, i, i);
            }
        }

        return power;
    } else {
        return get_null_mtx();
    }
}

double determinant(Matrix mtx)
{   
    if(!is_square(mtx) || is_null_mtx(mtx)) {
        return 0;
    } else if(mtx.num_cols == 1 && mtx.num_rows == 1) {
        return get_element_mtx(mtx, 0, 0);
    } else {
        double det = 0;

        for(size_t i = 0; i < mtx.num_cols; ++i) {
            double mult = get_element_mtx(mtx, 0, i);
            if(fabs(mult) > EPSILON) {
                if(i % 2 == 1) {
                    mult *= -1;
                }
                Matrix minor_mtx = get_submatrix_det(mtx, 0, i);
                double cofactor = mult * determinant(minor_mtx);

                det += cofactor;
                delete_matrix(minor_mtx);
            }
        }

        return det;
    }
}

static size_t get_row_index_with_closest_pivot(Matrix mtx)
{
    size_t min_piv_index = get_pivot_index(mtx.rows[0]);
    size_t row_index = 0;
    for(size_t i = 1; i < mtx.num_rows; ++i) {
        size_t curr_piv_index = get_pivot_index(mtx.rows[i]);
        if(curr_piv_index < min_piv_index) {
            min_piv_index = curr_piv_index;
            row_index = i;
        }
    }
    return row_index;
}

static void swap_and_set_rows(Matrix mtx, size_t row_index, double set_value) {
    if(row_index < mtx.num_rows) {
        swap_rows(mtx, 0, row_index);
        if(fabs(set_value) > EPSILON) {
            smultiply_inplace_vec(mtx.rows[0], (1 / set_value));
        }
    }
}

static void elim_row(Matrix mtx, Vector subtrahend, size_t row_index) {
    if(row_index < mtx.num_rows) {
        Vector new_row = subtract_vec(mtx.rows[row_index], subtrahend);
        set_row(new_row, mtx, row_index);
        delete_vector(new_row);
    }
}

Matrix get_ref(Matrix mtx)
{
    if(!is_null_mtx(mtx)) {
        Matrix ref_mtx = copy_matrix(mtx);
        get_ref_inplace(ref_mtx);

        return ref_mtx;
    } else {
        return get_null_mtx();
    }
}

void get_ref_inplace(Matrix mtx)
{
    if(is_null_mtx(mtx)) {
        return;
    } else if(mtx.num_rows == 1 || mtx.num_cols == 1) {
        if(mtx.num_rows == 1) {
            size_t pivot_index = get_pivot_index(mtx.rows[0]);
            if(pivot_index < mtx.num_cols) {
                double value = get_element_mtx(mtx, 0, pivot_index);
                smultiply_inplace_vec(mtx.rows[0], (1 / value));
            }
        } else {
            size_t row_min_piv_index = get_row_index_with_closest_pivot(mtx);
            if(row_min_piv_index < mtx.num_cols) {
                double value = get_element_mtx(mtx, row_min_piv_index, 0);
                swap_and_set_rows(mtx, row_min_piv_index, value);
                for(size_t i = 1; i < mtx.num_rows; ++i) {
                    double val_below_pivot = get_element_mtx(mtx, i, 0);
                    Vector subtrahend = smultiply_vec(mtx.rows[0], val_below_pivot);
                    elim_row(mtx, subtrahend, i);
                    delete_vector(subtrahend);
                }
            }
        }
    } else {
        size_t row_min_piv_index = get_row_index_with_closest_pivot(mtx);
        size_t  pivot_index = get_pivot_index(mtx.rows[row_min_piv_index]);
        if(pivot_index < mtx.num_cols) {        
            double pivot_val = get_element_mtx(mtx, row_min_piv_index, pivot_index);
            swap_and_set_rows(mtx, row_min_piv_index, pivot_val);

            for(size_t i = 1; i < mtx.num_rows; ++i) {
                double val_below_pivot = get_element_mtx(mtx, i, 0);
                Vector subtrahend = smultiply_vec(mtx.rows[0], val_below_pivot);
                elim_row(mtx, subtrahend, i);
                delete_vector(subtrahend);
            }
            Matrix sub_mtx = get_submatrix_true(mtx, 1, pivot_index + 1);
            get_ref_inplace(sub_mtx);
            set_submatrix_true(sub_mtx, mtx, 1, pivot_index + 1);
            delete_matrix(sub_mtx);
        }
    }
}

Matrix get_ref_mirror(Matrix mtx1, Matrix mtx2)
{
    if(!is_null_mtx(mtx1) && !is_null_mtx(mtx2) && mtx1.num_rows == mtx2.num_rows) {
        Matrix partitioned_mtx = hcombine(mtx1, mtx2);
        get_ref_inplace(partitioned_mtx);

        Matrix ref_mirror = get_submatrix_true(partitioned_mtx, 0, mtx1.num_cols);
        delete_matrix(partitioned_mtx);

        return ref_mirror;
    } else {
        return get_null_mtx();
    }
}

Matrix get_rref(Matrix mtx)
{
    if(!is_null_mtx(mtx)) {
        Matrix rref_mtx = copy_matrix(mtx);
        get_rref_inplace(rref_mtx);

        return rref_mtx;
    } else {
        return get_null_mtx();
    }
}

void get_rref_inplace(Matrix mtx)
{
    get_ref_inplace(mtx);
    if(!is_null_mtx(mtx) && !is_zero_mtx(mtx)) {    
        for(size_t i = 0; i < mtx.num_rows; ++i) {
            size_t pivot_index = get_pivot_index(mtx.rows[i]);
            if(pivot_index < mtx.num_cols) {
                double element = get_element_mtx(mtx, i, pivot_index);
                if(fabs(element - 1) < EPSILON) {
                    for(size_t j = 0; j < i; ++j) {
                        double val_above_pivot = get_element_mtx(mtx, j, pivot_index);
                        Vector subtrahend = smultiply_vec(mtx.rows[i], val_above_pivot);
                        elim_row(mtx, subtrahend, j);
                        delete_vector(subtrahend);
                    }
                }
            }
        }
    }
}

Matrix get_rref_mirror(Matrix mtx1, Matrix mtx2)
{
    if(!is_null_mtx(mtx1) && !is_null_mtx(mtx2) && mtx1.num_rows == mtx2.num_rows) {
        Matrix partitioned_mtx = hcombine(mtx1, mtx2);
        get_rref_inplace(partitioned_mtx);

        Matrix rref_mirror = get_submatrix_true(partitioned_mtx, 0, mtx1.num_cols);
        delete_matrix(partitioned_mtx);

        return rref_mirror;
    } else {
        return get_null_mtx();
    }
}

static Matrix make_householder(Matrix mtx, size_t minor_index, bool for_hess)
{
    Matrix householder_mtx = make_identity(mtx.num_rows);
    Matrix minor_mtx = get_submatrix_true(mtx, minor_index, minor_index);
    Matrix ident = make_identity(minor_mtx.num_rows);

    Vector minor_col;
    if(for_hess) {
        minor_col = get_subcol(mtx, minor_index - 1, minor_index);
    } else {
        minor_col = get_col(minor_mtx, 0);
    }
    double value = get_element_vec(minor_col, 0) - get_norm_vec(minor_col);
    set_element_vec(value, minor_col, 0);

    Vector unit_vec = get_unit_vec(minor_col);
    Matrix unit_vec_o_product = get_outer_product(unit_vec, unit_vec);
    smultiply_inplace_mtx(unit_vec_o_product, 2);

    Matrix hh_sub_mtx = subtract_mtx(ident, unit_vec_o_product);
    set_submatrix_true(hh_sub_mtx, householder_mtx, minor_index, minor_index);

    delete_matrix(minor_mtx);
    delete_matrix(ident);
    delete_matrix(unit_vec_o_product);
    delete_matrix(hh_sub_mtx);
    delete_vector(minor_col);
    delete_vector(unit_vec);

    return householder_mtx;
}

Matrix get_hessenberg(Matrix mtx)
{
    Matrix hessenberg_mtx = copy_matrix(mtx);
    get_hessenberg_inplace(hessenberg_mtx);

    return hessenberg_mtx;
}

void get_hessenberg_inplace(Matrix mtx)
{
    if(!is_null_mtx(mtx) && is_square(mtx)) {
        for(size_t i = 0; i < mtx.num_rows - 1; ++i) {
            Matrix householder_mtx = make_householder(mtx, i + 1, true);

            Matrix intermediate_mtx1 = mmultiply(householder_mtx, mtx);
            Matrix intermediate_mtx2 = mmultiply(intermediate_mtx1, householder_mtx);

            set_submatrix_true(intermediate_mtx1, mtx, 0, 0);

            delete_matrix(householder_mtx);
            delete_matrix(intermediate_mtx1);
            delete_matrix(intermediate_mtx2);
        }
    } 
}

Matrix * get_nndim_basis(size_t dim)
{
    if(dim > 0) {
        Matrix reference = make_matrix(dim, dim);
        Matrix *basis = malloc(dim * dim * sizeof(reference));
        delete_matrix(reference);
        if(!basis) {
            return NULL;
        }


        for(size_t i = 0; i < dim * dim; ++i) {
            basis[i] = make_matrix(dim, dim);
            size_t row_index = i / dim;
            size_t col_index = i % dim;
            set_element_mtx(1, basis[i], row_index, col_index);
        }
        return basis;
    } else {
        return NULL;
    }
}

Vector * get_nullspace_basis(const Matrix mtx)
{
    Matrix rref_mtx = get_rref(mtx);
    if(!is_null_mtx(mtx) && !is_null_mtx(rref_mtx)) {
        size_t dim_basis = 0;
        size_t cols_with_pivot[mtx.num_cols];
        for(size_t i = 0; i < mtx.num_cols; ++i) {
            cols_with_pivot[i] = mtx.num_rows + mtx.num_cols;
        }
        for(size_t i = 0; i < mtx.num_rows; ++i) {
            size_t pivot_index = get_pivot_index(rref_mtx.rows[i]);
            if(pivot_index < mtx.num_cols) {
                cols_with_pivot[pivot_index] = i;
                for(size_t j = pivot_index + 1; j < mtx.num_cols; ++j) {
                    double value = get_element_mtx(rref_mtx, i, j);
                    if(fabs(value) > EPSILON && cols_with_pivot[j] == mtx.num_cols + mtx.num_rows) {
                        cols_with_pivot[j] = mtx.num_rows + dim_basis;
                        ++dim_basis;
                    } else if(cols_with_pivot[j] == mtx.num_cols + mtx.num_rows) {
                        bool is_free_zeroes = true;
                        for(size_t k = i; k < mtx.num_rows; ++k) {
                            if(fabs(get_element_mtx(rref_mtx, k, j)) > EPSILON) {
                                is_free_zeroes = false;
                                break;
                            }
                        }
                        if(is_free_zeroes) {
                            cols_with_pivot[j] = mtx.num_rows + dim_basis;
                            ++dim_basis;                           
                        }
                    }
                }
            } else {
                size_t num_pre_zero_cols = 0;
                for(size_t j = 0; j < mtx.num_cols; ++j) {
                    if(cols_with_pivot[j] == (mtx.num_rows + mtx.num_cols)) {
                        cols_with_pivot[j] = mtx.num_rows + num_pre_zero_cols;
                        ++num_pre_zero_cols;
                        ++dim_basis;
                    } else if(cols_with_pivot[j] >= mtx.num_rows) {
                        cols_with_pivot[j] += num_pre_zero_cols;
                    }
                }
            }
        }
        Vector *basis = malloc(dim_basis * sizeof(rref_mtx.rows[0]));
        if(!basis || dim_basis == 0) {
            delete_matrix(rref_mtx);
            return NULL;
        }
        for(size_t i = 0; i < dim_basis; ++i) {
            basis[i] = make_vector(mtx.num_cols, COL);
        }
        
        for(size_t i = 0; i < mtx.num_cols; ++i) {
            size_t row_piv_index = cols_with_pivot[i];
            if(row_piv_index < mtx.num_rows) {
                for(size_t j = i + 1; j < mtx.num_cols; ++j) {
                    double value = get_element_mtx(rref_mtx, row_piv_index, j);
                    if(fabs(value) > EPSILON) {
                        size_t basis_index = cols_with_pivot[j] - mtx.num_rows;
                        set_element_vec(-1 * value, basis[basis_index], i);
                    }
                }
            } else if(row_piv_index < mtx.num_cols + mtx.num_rows) {
                size_t basis_index = row_piv_index - mtx.num_rows;
                set_element_vec(1, basis[basis_index], i);
            }
        }
        delete_matrix(rref_mtx);
        
        return basis;
    } else {
        delete_matrix(rref_mtx);
        return NULL;
    }
}

Vector * get_rowspace_basis(const Matrix mtx)
{
    Matrix rref_mtx = get_rref(mtx);
    if(!is_null_mtx(mtx) && !is_null_mtx(rref_mtx)) {
        size_t basis_dim_ctr = 0;
        Vector *basis = malloc(mtx.num_rows * sizeof(mtx.rows[0]));
        if(!basis) {
            return NULL;
        }

        for(size_t i = 0; i < rref_mtx.num_rows; ++i) {
            size_t pivot_index = get_pivot_index(rref_mtx.rows[i]);
            if(pivot_index < rref_mtx.num_cols) {
                basis[basis_dim_ctr] = get_row(rref_mtx, i);
                ++basis_dim_ctr;
            }
        }
        delete_matrix(rref_mtx);

        Vector *resized_basis = realloc(basis, basis_dim_ctr * sizeof(mtx.rows[0]));
        if(resized_basis) {
            return resized_basis;
        } else {
            return basis;
        }
    } else {
        return NULL;
    }
}

Vector * get_colspace_basis(const Matrix mtx)
{
    Matrix rref_mtx = get_rref(mtx);
    if(!is_null_mtx(mtx) && !is_null_mtx(rref_mtx)) {
        size_t basis_dim_ctr = 0;
        Vector reference = get_col(rref_mtx, 0);
        Vector *basis = malloc(mtx.num_cols * sizeof(reference));
        if(!basis) {
            delete_vector(reference);
            return NULL;
        }

        for(size_t i = 0; i < rref_mtx.num_rows; ++i) {
            size_t pivot_index = get_pivot_index(rref_mtx.rows[i]);
            if(pivot_index < rref_mtx.num_cols) {
                basis[basis_dim_ctr] = get_col(mtx, pivot_index);
                ++basis_dim_ctr;
            }
        }
        Vector *resized_basis = realloc(basis, basis_dim_ctr * sizeof(reference));

        delete_matrix(rref_mtx);
        delete_vector(reference);
        
        if(resized_basis) {
            return resized_basis;
        } else {
            return basis;
        }        
    } else {
        return NULL;
    }
}

Vector * get_colcomp_basis(const Matrix mtx)
{
    Matrix t_mtx = transpose_mtx(mtx);
    Vector *basis = get_nullspace_basis(t_mtx);
    delete_matrix(t_mtx);

    return basis;
}

Vector * get_nullcomp_basis(const Matrix mtx)
{
    Matrix t_mtx = transpose_mtx(mtx);
    Vector *basis = get_colspace_basis(t_mtx);
    delete_matrix(t_mtx);

    return basis;
}

double get_trace(Matrix mtx)
{
    if(!is_null_mtx(mtx) && is_square(mtx)) {
        double trace = 0;
        for(size_t i = 0; i < mtx.num_rows; ++i) {
            double diagonal_value = get_element_mtx(mtx, i, i);
            trace += diagonal_value;
        }
        return trace;
    } else {
        return nan("");
    }
}

size_t get_rank(Matrix mtx)
{
    if(!is_null_mtx(mtx)) {
        Matrix rref_mtx = get_rref(mtx);

        size_t rnk = 0;
        for(size_t i = 0; i < rref_mtx.num_rows; ++i) {
            if(!is_zero_vec(rref_mtx.rows[i])) {
                ++rnk;
            }
        }
        delete_matrix(rref_mtx);

        return rnk;
    } else {
        return mtx.num_cols + mtx.num_rows;
    }
}

size_t get_nullity(Matrix mtx)
{
    if(!is_null_mtx(mtx)) {
        return mtx.num_cols - get_rank(mtx);
    } else {
        return mtx.num_cols;
    }
}

Matrix get_inverse(Matrix mtx)
{
    if(is_invertible(mtx)) {
        Matrix ident = make_identity(mtx.num_cols);

        Matrix inverse = get_rref_mirror(mtx, ident);

        delete_matrix(ident);

        return inverse;
    } else {
        return get_null_mtx();
    }
}

double get_inner_product(Matrix mtx1, Matrix mtx2)
{
    if(same_size_mtx(mtx1, mtx2) && !is_null_mtx(mtx1)) {
        double product = 0;
        for(size_t i = 0; i < mtx1.num_rows; ++i) {
            for(size_t j = 0; j < mtx1.num_cols; ++j) {
                double val1 = get_element_mtx(mtx1, i, j);
                double val2 = get_element_mtx(mtx2, i, j);

                product += (val1 * val2);
            }
        }
        return product;
    } else {
        return nan("");
    }
}

Matrix get_outer_product(Vector vec1, Vector vec2)
{
    if(!is_null_vec(vec1) && !is_null_vec(vec2)) {
        Matrix product = make_matrix(vec1.dim, vec2.dim);
        for(size_t i = 0; i < product.num_rows; ++i) {
            for(size_t j = 0; j < product.num_cols; ++j) {
                double value = get_element_vec(vec1, i) * get_element_vec(vec2, j);
                set_element_mtx(value, product, i, j);
            }
        }
        return product;
    } else {
        return get_null_mtx();
    }
}

double get_norm_mtx(Matrix mtx)
{
    if(!is_null_mtx(mtx)) {
        double norm_square = get_inner_product(mtx, mtx);
        return sqrt(norm_square);
    } else {
        return -1;
    }
}

double get_angle_mtx(Matrix mtx1, Matrix mtx2)
{
    double angle_rads = get_angle_rads_mtx(mtx1, mtx2);
    if(angle_rads >= 0) {
        return angle_rads * (180 / acos(-1));
    } else {
        return -1;
    }
}

double get_angle_rads_mtx(Matrix mtx1, Matrix mtx2)
{
    if(same_size_mtx(mtx1, mtx2) && !is_null_mtx(mtx1)) {
        double inner_product = get_inner_product(mtx1, mtx2);
        double mag1 = get_norm_mtx(mtx1);
        double mag2 = get_norm_mtx(mtx2);
        if(!isnan(inner_product) && mag1 > 0 && mag2 > 0) {
            return acos(inner_product / (mag1 * mag2));
        } else {
            return -1;
        }
    } else {
        return -1;
    }
}

double get_distance_mtx(Matrix mtx1, Matrix mtx2)
{
    if(same_size_mtx(mtx1, mtx2) && !is_null_mtx(mtx1) && !is_zero_mtx(mtx1) && !is_zero_mtx(mtx2)) {
        Matrix distance_mtx = subtract_mtx(mtx1, mtx2);
        double distance = get_norm_mtx(distance_mtx);
        delete_matrix(distance_mtx);

        return distance;
    } else {
        return -1;
    }
}

Matrix get_unit_mtx(Matrix mtx)
{
    if(!is_null_mtx(mtx) && !is_zero_mtx(mtx)) {
        return smultiply_mtx(mtx, (1 / get_norm_mtx(mtx)));
    } else {
        return get_null_mtx();
    }
}

double get_scalar_proj_mtx(Matrix mtx1, Matrix mtx2)
{
    if(same_size_mtx(mtx1, mtx2) && !is_null_mtx(mtx1)) {
        double inner_product = get_inner_product(mtx1, mtx2);
        double mag2 = get_norm_mtx(mtx2);
        if(!isnan(inner_product) && mag2 > 0) {
            return inner_product / mag2;
        } else {
            return -1;
        }
    } else {
        return -1;
    }
}

Matrix get_vector_proj_mtx(Matrix mtx1, Matrix mtx2)
{
    if(same_size_mtx(mtx1, mtx2) && !is_null_mtx(mtx1)) {
        double scalar_proj =  get_scalar_proj_mtx(mtx1, mtx2);
        Matrix unit_mtx = get_unit_mtx(mtx2);
        if(scalar_proj >= 0 && !is_null_mtx(unit_mtx)) {
            Matrix vector_proj = smultiply_mtx(unit_mtx, scalar_proj);
            delete_matrix(unit_mtx);

            return vector_proj;
        } else {
            delete_matrix(unit_mtx);
            return get_null_mtx();
        }
    } else {
        return get_null_mtx();
    }
}

Matrix get_projection_mtx(Matrix *basis, size_t dim, Matrix mtx)
{
    if(basis && dim > 0 && !is_null_mtx(mtx)) {
        Matrix projection = make_matrix(mtx.num_rows, mtx.num_cols);
        for(size_t i = 0; i < dim; ++i) {
            double inner_product = get_inner_product(mtx, basis[i]);
            if(isnan(inner_product)) {
                delete_matrix(projection);
                return get_null_mtx();
            }
            Matrix addend = smultiply_mtx(basis[i], inner_product);
            Matrix temp_mtx = copy_matrix(projection);
            delete_matrix(projection);

            projection = add_mtx(temp_mtx, addend);

            delete_matrix(temp_mtx);
            delete_matrix(addend);
        }
        return projection;
    } else {
        return get_null_mtx();
    }
}

Matrix * orthonormalize_set_mtx(Matrix *mtxs, size_t dim)
{
    if(mtxs && dim > 0) {
        Matrix *orth_mtxs = malloc(dim * sizeof(mtxs[0]));
        if(!orth_mtxs) {
            return NULL;
        }
        orth_mtxs[0] = get_unit_mtx(mtxs[0]);

        for(size_t i = 1; i < dim; ++i) {
            Matrix projection = get_projection_mtx(orth_mtxs, i, mtxs[i]);
            Matrix orth_mtx = subtract_mtx(mtxs[i], projection);

            orth_mtxs[i] = get_unit_mtx(orth_mtx);

            delete_matrix(projection);
            delete_matrix(orth_mtx);
        }
        return orth_mtxs;
    } else {
        return NULL;
    }
}

void qr_decomposition(Matrix mtx, Matrix *qr)
{
    if(!is_null_mtx(mtx) && mtx.num_rows >= mtx.num_cols) {
        Matrix q_tot = make_identity(mtx.num_rows);
        Matrix r_tot = copy_matrix(mtx);
        Matrix intermediate_mtx = copy_matrix(mtx);
        size_t num_iterations = mtx.num_cols;
        if(mtx.num_cols > mtx.num_rows - 1) {
            num_iterations = mtx.num_rows - 1;
        }
        for(size_t i = 0; i < num_iterations; ++i) {
            Matrix q_i = make_householder(intermediate_mtx, i, false);
            Matrix intermediate_mtx_temp = copy_matrix(intermediate_mtx);
            delete_matrix(intermediate_mtx);

            intermediate_mtx = mmultiply(q_i, intermediate_mtx_temp);

            Matrix q_tot_temp = copy_matrix(q_tot);
            Matrix q_i_transpose = transpose_mtx(q_i);
            delete_matrix(q_tot);
            q_tot = mmultiply(q_tot_temp, q_i_transpose);

            Matrix r_tot_temp = copy_matrix(r_tot);
            delete_matrix(r_tot);
            r_tot = mmultiply(q_i, r_tot_temp);

            delete_matrix(q_i);
            delete_matrix(q_tot_temp);
            delete_matrix(q_i_transpose);
            delete_matrix(r_tot_temp);
            delete_matrix(intermediate_mtx_temp);
        }
        if(is_square(mtx)) {
            qr[0] = copy_matrix(q_tot);
            qr[1] = copy_matrix(r_tot);
        } else {
            qr[0] = make_matrix(mtx.num_rows, mtx.num_cols);
            qr[1] = make_matrix(mtx.num_cols, mtx.num_cols);
            for(size_t i = 0; i < mtx.num_rows; ++i) {
                for(size_t j = 0; j < mtx.num_cols; ++j) {
                    if(i < mtx.num_cols) {
                        double r_val = get_element_mtx(r_tot, i, j);
                        set_element_mtx(r_val, qr[1], i, j);
                    }
                    double q_val = get_element_mtx(q_tot, i, j);
                    set_element_mtx(q_val, qr[0], i, j);
                }
            }
        }

        delete_matrix(q_tot);
        delete_matrix(r_tot);
        delete_matrix(intermediate_mtx);
    } else {
        qr[0] = get_null_mtx();
        qr[1] = get_null_mtx();
    }
}

void rq_decomposition(Matrix mtx, Matrix *rq)
{
    if(!is_null_mtx(mtx) && mtx.num_rows >= mtx.num_cols && get_rank(mtx) == mtx.num_cols) {
        Vector *rows = malloc(mtx.num_rows * sizeof(mtx.rows[0]));
        if(rows) {
            for(size_t i = 0; i < mtx.num_rows; ++i) {
                rows[i] = get_row(mtx, mtx.num_rows - i - 1);
            }
            
            Vector *orth_rows = orthonormalize_set_vec(rows, mtx.num_rows);
            for(size_t i = 0; i < mtx.num_rows; ++i) {
                delete_vector(rows[i]);
            }
            free(rows);
            if(orth_rows) {
                Matrix q = make_matrix_rows(orth_rows, mtx.num_rows);
                Matrix q_t = transpose_mtx(q);
                Matrix r = mmultiply(mtx, q_t);

                rq[0] = r;
                rq[1] = q;

                delete_matrix(q_t);
                for(size_t i = 0; i < mtx.num_rows; ++i) {
                    delete_vector(orth_rows[i]);
                }
                free(orth_rows);

                return;
            }
        }
    }
    rq[0] = get_null_mtx();
    rq[1] = get_null_mtx();  
}

static double calculate_wilkinson_shift(double a, double b, double c)
{
    if(fabs(b) < EPSILON) {
        return c;
    }
    double delta = (a - c) / 2;
    int sign;
    if(delta > 0) {
        sign = 1;
    } else {
        sign = -1;
    }

    return c - sign * pow(b, 2)/(fabs(delta) + sqrt(pow(delta, 2) + pow(b, 2)));
}

Vector get_eigenvalues(Matrix mtx)
{
    if(!is_null_mtx(mtx) && is_square(mtx)) {
        Vector eigenvalues = make_vector(1, ROW);
        if(mtx.num_rows == 1) {
            double eigenvalue = get_element_mtx(mtx, 0, 0);
            set_element_vec(eigenvalue, eigenvalues, 0);
        } else {
            Matrix iter_mtx = get_hessenberg(mtx);
            double cutoff_val = 1;
            while(fabs(cutoff_val) > pow(10, -10)) {
                double a = get_element_mtx(iter_mtx, mtx.num_rows - 2, mtx.num_rows - 2);
                double b = get_element_mtx(iter_mtx, mtx.num_rows - 1, mtx.num_rows - 2);
                if(fabs(b) < EPSILON) {
                    break;
                }
                double c = get_element_mtx(iter_mtx, mtx.num_rows - 1, mtx.num_rows - 1);
                double shift = calculate_wilkinson_shift(a, b, c);
                Matrix shift_amount_mtx = make_identity(mtx.num_rows);
                smultiply_inplace_mtx(shift_amount_mtx, shift);

                Matrix shifted_mtx = subtract_mtx(mtx, shift_amount_mtx);

                Matrix qr[2];
                qr_decomposition(shifted_mtx, qr);

                delete_matrix(iter_mtx);
                
                Matrix r_times_q = mmultiply(qr[1], qr[0]);
                iter_mtx = add_mtx(r_times_q, shift_amount_mtx);

                cutoff_val = get_element_mtx(iter_mtx, mtx.num_rows - 1, mtx.num_cols - 2);

                delete_matrix(shift_amount_mtx);
                delete_matrix(shifted_mtx);
                delete_matrix(r_times_q);
                delete_matrix(qr[0]);
                delete_matrix(qr[1]);
            }
            double eigenvalue = get_element_mtx(iter_mtx, mtx.num_rows - 1, mtx.num_cols - 1);
            set_element_vec(eigenvalue, eigenvalues, 0);

            Matrix sub_mtx = get_submatrix_det(iter_mtx, mtx.num_rows - 1, mtx.num_cols - 1);
            Vector other_eigenvalues = get_eigenvalues(sub_mtx);

            Vector temp1 = combine(eigenvalues, other_eigenvalues);
            Vector temp2 = make_vector_unique(temp1);
            delete_vector(eigenvalues);
            eigenvalues = copy_vector(temp2);

            delete_matrix(iter_mtx);
            delete_matrix(sub_mtx);
            delete_vector(other_eigenvalues);
            delete_vector(temp1);
            delete_vector(temp2);
        }

        return eigenvalues;
    } else {
        return get_null_vec();
    }
}

static Matrix get_eigenspace_mtx(Matrix mtx, double eigenvalue)
{
    Matrix subtrahend = make_identity(mtx.num_cols);
    smultiply_inplace_mtx(subtrahend, eigenvalue);
    Matrix eigenspace_mtx = subtract_mtx(mtx, subtrahend);

    delete_matrix(subtrahend);

    return eigenspace_mtx;
}

Vector * get_eigenspace_basis(Matrix mtx, double eigenvalue)
{
    if(!is_null_mtx(mtx) && is_square(mtx)) {
        Matrix eigenspace_mtx = get_eigenspace_mtx(mtx, eigenvalue);

        Vector *basis = get_nullspace_basis(eigenspace_mtx);

        delete_matrix(eigenspace_mtx);

        return basis;
    } else {
        return NULL;
    }
}

size_t get_eigenspace_basis_dim(Matrix mtx, double eigenvalue)
{
    if(!is_null_mtx(mtx) && is_square(mtx)) {
        Matrix eigenspace_mtx = get_eigenspace_mtx(mtx, eigenvalue);

        size_t dim_basis = get_nullity(eigenspace_mtx);

        delete_matrix(eigenspace_mtx);

        return dim_basis;
    } else {
        return 0;
    }
}



Vector get_eigenvector(Matrix mtx, double eigenvalue)
{
    if(!is_null_mtx(mtx) && is_square(mtx)) {
        Vector *basis = get_eigenspace_basis(mtx, eigenvalue);
        size_t dim_basis = get_eigenspace_basis_dim(mtx, eigenvalue);

        double coords[dim_basis];
        for(size_t i = 0; i < dim_basis; ++i) {
            coords[i] = 1;
        }

        Vector eigenvector = lin_combination_dbls(basis, dim_basis, coords);

        for(size_t i = 0; i < dim_basis; ++i) {
            delete_vector(basis[i]);
        }
        free(basis);

        return eigenvector;
    } else {
        return get_null_vec();
    }
}

Vector * get_eigenvectors(Matrix mtx, Vector eigenvalues)
{
    if(!is_null_mtx(mtx) && !is_null_vec(eigenvalues) && is_square(mtx)) {
        Vector *eigenvectors = malloc(eigenvalues.dim * sizeof(mtx.rows[0]));
        if(eigenvectors) {        
            for(size_t i = 0; i < eigenvalues.dim; ++i) {
                double eigenvalue = get_element_vec(eigenvalues, i);
                eigenvectors[i] = get_eigenvector(mtx, eigenvalue);
            }
        }
        return eigenvectors;
    } else {
        return NULL;
    }
}

void diagonalize(Matrix mtx, Matrix *xdx_1)
{
    if(!is_null_mtx(mtx) && is_square(mtx)) {
        Vector eigenvalues = get_eigenvalues(mtx);
        if(eigenvalues.dim == mtx.num_cols) {
            Vector *eigenvectors = get_eigenvectors(mtx, eigenvalues);
            if(eigenvectors) {
                Matrix x = make_matrix_cols(eigenvectors, eigenvalues.dim);
                Matrix x_1 = get_inverse(x);

                Matrix d = make_matrix(mtx.num_rows, mtx.num_cols);
                for(size_t i = 0; i < eigenvalues.dim; ++i) {
                    double eigenvalue = get_element_vec(eigenvalues, i);
                    set_element_mtx(eigenvalue, d, i, i);
                }
                for(size_t i = 0; i < eigenvalues.dim; ++i) {
                    delete_vector(eigenvectors[i]);
                }
                free(eigenvectors);

                xdx_1[0] =copy_matrix(x);
                xdx_1[1] =copy_matrix(d);
                xdx_1[2] =copy_matrix(x_1);

                delete_matrix(x);
                delete_matrix(d);
                delete_matrix(x_1);
                delete_vector(eigenvalues);

                return;
            }
        }
        delete_vector(eigenvalues);
    }
}
