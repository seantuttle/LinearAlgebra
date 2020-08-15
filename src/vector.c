
#include "vector.h"
#include "matrix.h"
#include "equations.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

const Vector NULL_VECTOR = {0, ROW, NULL}; // singleton value for Vector NULL

// Number of decimal places is abs val of this constant - 1
const int PRINT_EPSILON_POW = -3;

Vector make_vector(size_t dim, bool type)
{
    if(dim > 0) {
        Vector new_vec = {dim, type, malloc(dim * sizeof(double))};
        if(!new_vec.elements) {
            return get_null_vec();
        }

        zeroes_vec(new_vec); // fill new Vector with zeroes

        return new_vec;
    } else {
        return get_null_vec();
    }
}

Vector copy_vector(Vector vec)
{
    if(!is_null_vec(vec)) {
        Vector new_vec = make_vector(vec.dim, vec.type);
        for(size_t i = 0; i < vec.dim; ++i) {
            double value = get_element_vec(vec, i);
            set_element_vec(value, new_vec, i);
        }

        return new_vec;
    } else {
        return get_null_vec();
    }
}

Vector make_vector_vals(double *vals, size_t dim, bool type)
{
    if(dim > 0 && vals) {
        Vector vec = make_vector(dim, type);
        for(size_t i = 0; i < vec.dim; ++i) {
            set_element_vec(vals[i], vec, i);
        }

        return vec;
    } else {
        return get_null_vec();
    }
}

Vector make_vector_rand(size_t dim, bool type, unsigned int max_abs_val)
{
    Vector vec = make_vector(dim, type);
    if(!is_null_vec(vec)) {
        int max_val_standin = (2 * max_abs_val) + 1;
        srand((unsigned) time(NULL)); // set the seed to the current time
        for(size_t i = 0; i < vec.dim; ++i) {
            double value = rand() % max_val_standin; // in range [0, 2*max_abs_val]
            set_element_vec(value - max_abs_val, vec, i);
        }
    }
    return vec;
}

Vector make_vector_unique(Vector vec)
{
    if(!is_null_vec(vec)) {
        size_t num_unique_vals = 0;
        double unique_vals[vec.dim];
        for(size_t i = 0; i < vec.dim; ++i) {
            unique_vals[i] = nan("");
        }

        // go through and make new Vector that has all
        // the values of original Vector but no repeats
        for(size_t i = 0; i < vec.dim; ++i) {
            bool already_added = false;
            double value = get_element_vec(vec, i);
            for(size_t i = 0; i < num_unique_vals; ++i) {
                if(fabs(unique_vals[i] - value) < EPSILON) {
                    already_added = true;
                    break;
                }
            }
            if(!already_added) {
                unique_vals[num_unique_vals] = value;
                ++num_unique_vals;
            }
        }
        return make_vector_vals(unique_vals, num_unique_vals, vec.type);
    } else {
        return get_null_vec();
    }
}

void delete_vector(Vector vec)
{
    free(vec.elements);
}

Vector combine(Vector vec1, Vector vec2)
{
    if(same_type(vec1, vec2) && !is_null_vec(vec1) && !is_null_vec(vec2)) {
        Vector combo = make_vector(vec1.dim + vec2.dim, vec1.type);
        set_subvector(vec1, combo, 0);
        set_subvector(vec2, combo, vec1.dim);

        return combo;
    } else {
        return get_null_vec();
    }
}

void print_vector(Vector vec)
{
    for(size_t i = 0; i < vec.dim; ++i) {
        double value = get_element_vec(vec, i);
        
        // Make sure that it prints 0 and not -0
        if(fabs(value) < pow(10, PRINT_EPSILON_POW)) {
            value = 0;
        }
        printf("%.2f", value);
        if(is_row(vec)) {
            printf("\t"); // print in a row
        } else {
            printf("\n"); // print in a column
        }
    }
    if(is_row(vec)) {
        printf("\n");
    }
}

void zeroes_vec(Vector vec)
{
    fill_vec(vec, 0);
}

void ones_vec(Vector vec)
{
    fill_vec(vec, 1);
}

void fill_vec(Vector vec, double value)
{
    for(size_t i = 0; i < vec.dim; ++i) {
        set_element_vec(value, vec, i);
    }
}

void set_element_vec(double value, Vector vec, size_t index)
{
    if(in_range_vec(vec, index)) {
        vec.elements[index] = value;
    }
}

void set_subvector(Vector sub_vec, Vector vec, size_t index)
{
    // make sure that new subvector fits, then copy over values
    if(index + sub_vec.dim <= vec.dim) {
        for(size_t i = 0; i < sub_vec.dim; ++i) {
            double value = get_element_vec(sub_vec, i);
            set_element_vec(value, vec, i + index);
        }
    }
}

double get_element_vec(Vector vec, size_t index)
{
    return vec.elements[index];
}

double * get_elements_vec(Vector vec)
{
    double *elements = malloc(vec.dim * sizeof *elements);
    if(elements) {
        for(size_t i = 0; i < vec.dim; ++i) {
            elements[i] = get_element_vec(vec, i);
        }
    }
    return elements;
}

Vector get_subvector(Vector vec, size_t index)
{
    if(!is_null_vec(vec) && in_range_vec(vec, index)) {
        // make new Vector the size of subvector starting at index
        Vector sub_vec = make_vector(vec.dim - index, vec.type);
        for(size_t i = 0; i < sub_vec.dim; ++i) {
            double value = get_element_vec(vec, i + index);
            set_element_vec(value, sub_vec, i);
        }
        return sub_vec;
    } else {
        return get_null_vec();
    }
}

size_t get_pivot_index(Vector vec)
{
    // pivot is first nonzero element in a Vector
    for(size_t i = 0; i < vec.dim; ++i) {
        if(fabs(get_element_vec(vec, i)) > EPSILON) {
            return i;
        }
    }
    return vec.dim; // pivot does not exist
}

Vector get_null_vec(void)
{
    return NULL_VECTOR;
}

Vector * get_ndim_standard_basis(size_t dim)
{
    Vector reference = make_vector(dim, COL); // size reference
    Vector *basis = malloc(dim * sizeof(reference));
    delete_vector(reference);

    for(size_t i = 0; i < dim; ++i) {
        Vector basis_vector = make_vector(dim, COL);
        set_element_vec(1, basis_vector, i); // make standard basis vector
        basis[i] = basis_vector;
    }
    return basis;
}

double get_norm_vec(Vector vec)
{
    if(!is_null_vec(vec)) {
        double norm_square = dot(vec, vec);
        return sqrt(norm_square);
    } else {
        return -1; // illegal value, norm can't be < 0
    }
}

double get_angle_vec(Vector vec1, Vector vec2)
{
    double angle_rads = get_angle_rads_vec(vec1, vec2);
    if(angle_rads >= 0) {
        return angle_rads * (180 / acos(-1)); // acos(-1) == PI
    } else {
        return -1; // illegal value, angle can't be < 0
    }
}

double get_angle_rads_vec(Vector vec1, Vector vec2)
{
    if(same_size_and_type(vec1, vec2) && !is_null_vec(vec1)) {
        double dot_product = dot(vec1, vec2);
        double mag1 = get_norm_vec(vec1);
        double mag2 = get_norm_vec(vec2);
        if(!isnan(dot_product) && mag1 > 0 && mag2 > 0) {
            return acos(dot_product / (mag1 * mag2));
        } else {
            return -1; // illegal value, angle can't be < 0
        }
    } else {
        return -1; // illegal value, angle can't be < 0
    }
}

double get_distance_vec(Vector vec1, Vector vec2)
{
    if(same_size_and_type(vec1, vec2) && !is_null_vec(vec1) && !is_zero_vec(vec1) && !is_zero_vec(vec2)) {
        Vector distance_vec = subtract_vec(vec1, vec2);
        double distance = get_norm_vec(distance_vec);
        delete_vector(distance_vec);

        return distance;
    } else {
        return -1;
    }
}

Vector get_unit_vec(Vector vec)
{
    if(!is_null_vec(vec) && !is_zero_vec(vec)) {
        return smultiply_vec(vec, (1 / get_norm_vec(vec)));
    } else {
        return get_null_vec();
    }
}

double get_scalar_proj_vec(Vector vec1, Vector vec2)
{
    if(same_size_and_type(vec1, vec2) && !is_null_vec(vec1)) {
        double dot_product = dot(vec1, vec2);
        double mag2 = get_norm_vec(vec2);
        if(!isnan(dot_product) && mag2 > 0) {
            return dot_product / mag2;
        } else {
            return -1;
        }
    } else {
        return -1;
    }
}

Vector get_vector_proj_vec(Vector vec1, Vector vec2)
{
    if(same_size_and_type(vec1, vec2) && !is_null_vec(vec1)) {
        double scalar_proj =  get_scalar_proj_vec(vec1, vec2);
        Vector unit_vec = get_unit_vec(vec2);
        if(scalar_proj >= 0 && !is_null_vec(unit_vec)) {
            Vector vector_proj = smultiply_vec(unit_vec, scalar_proj);
            delete_vector(unit_vec);

            return vector_proj;
        } else {
            delete_vector(unit_vec);
            return get_null_vec();
        }
    } else {
        return get_null_vec();
    }
}

Vector get_projection_vec(Vector *basis, size_t dim, Vector vec)
{
    if(basis && dim > 0 && !is_null_vec(vec)) {
        Vector projection = make_vector(vec.dim, vec.type);
        for(size_t i = 0; i < dim; ++i) {
            double dot_product = dot(vec, basis[i]);
            if(isnan(dot_product)) {
                delete_vector(projection);
                return get_null_vec();
            }
            Vector addend = smultiply_vec(basis[i], dot_product);
            Vector temp_vec = copy_vector(projection);
            delete_vector(projection);

            projection = add_vec(temp_vec, addend);

            delete_vector(temp_vec);
            delete_vector(addend);
        }
        return projection;
    } else {
        return get_null_vec();
    }
}

bool in_range_vec(Vector vec, size_t index)
{
    return index < vec.dim;
}

bool same_size_vec(Vector vec1, Vector vec2)
{
    return vec1.dim == vec2.dim;
}

bool same_type(Vector vec1, Vector vec2)
{
    return vec1.type == vec2.type;
}

bool same_size_and_type(Vector lhs, Vector rhs)
{
    return same_size_vec(lhs, rhs) && same_type(lhs, rhs);
}

bool equal_vec(Vector lhs, Vector rhs)
{
    if(same_size_and_type(lhs, rhs)) {
        for(size_t i = 0; i < lhs.dim; ++i) {
            double lhs_val = get_element_vec(lhs, i);
            double rhs_val = get_element_vec(rhs, i);
            if(fabs(lhs_val - rhs_val) > EPSILON) {
                return false;
            }
        }
        return true;
    } else {
        return false;
    }
}

bool can_dot(Vector lhs, Vector rhs)
{
    if((lhs.type == ROW) && (rhs.type == COL) && !is_null_vec(lhs)) {
        return same_size_vec(lhs, rhs);
    } else {
        return false;
    }
}

bool are_lin_ind(Vector *vecs, size_t dim)
{
    if(vecs && dim > 0) {
        Matrix mtx;
        if(is_col(vecs[0])) {
            mtx = make_matrix_cols(vecs, dim);
        } else {
            mtx = make_matrix_rows(vecs, dim);
        }
        size_t rank = get_rank(mtx);

        delete_matrix(mtx);

        return rank == dim;
    } else {
        return false;
    }
}

bool are_lin_dep(Vector *vecs, size_t dim)
{
    return !are_lin_ind(vecs, dim);
}

bool is_ndim_basis(Vector *vecs, size_t dim)
{
    if(vecs && dim > 0) {
        if(vecs[0].dim == dim) {
            return are_lin_ind(vecs, dim);
        } else {
            return false;
        }
    } else {
        return false;
    }   
}

bool are_orthogonal_vec(Vector vec1, Vector vec2)
{
    if(same_size_and_type(vec1, vec2) && !is_null_vec(vec1)) {
        double dot_product = dot(vec1, vec2);
        if(!isnan(dot_product)) {
            return fabs(dot_product) <EPSILON;
        } else {
            return false;
        }
    } else {
        return false;
    }
}

bool are_orthogonal_sets_vec(Vector *vecs1, Vector *vecs2, size_t dim1, size_t dim2)
{
    if(vecs1 && vecs2 && dim1 > 0 && dim2 > 0) {
        for(size_t i = 0; i < dim1; ++i) {
            for(size_t j = 0; j < dim2; ++j) {
                if(!are_orthogonal_vec(vecs1[i], vecs2[j])) {
                    return false;
                }
            }
        }
        return true;
    } else {
        return false;
    }
}

static bool check_orthonormal_or_orthogonal(Vector *vecs, size_t dim, bool (*func)(Vector, Vector))
{
    if(vecs && dim > 0) {
        for(size_t i = 0; i < dim; ++i) {
            for(size_t j = 0; j < dim; ++j) {
                if(i != j && !func(vecs[i], vecs[j])) {
                    return false;
                }
            }
        }
        return true;
    } else {
        return false;
    }
}

bool is_orthogonal_set_vec(Vector *vecs, size_t dim)
{
    return check_orthonormal_or_orthogonal(vecs, dim, are_orthogonal_vec);
}

static bool orthonormal_requirement(Vector vec1, Vector vec2)
{
    return are_orthogonal_vec(vec1, vec2) && is_unit_vec(vec1) && is_unit_vec(vec2);
}

bool is_orthonormal_set_vec(Vector *vecs, size_t dim)
{
    return check_orthonormal_or_orthogonal(vecs, dim, orthonormal_requirement);
}

bool is_unit_vec(Vector vec)
{
    if(!is_null_vec(vec)) {
        return fabs(get_norm_vec(vec) - 1) < EPSILON;
    } else {
        return false;
    }
}

bool is_row(Vector vec)
{
    return vec.type == ROW;
}

bool is_col(Vector vec)
{
    return !is_row(vec);
}

bool is_zero_vec(Vector vec)
{
    if(!is_null_vec(vec)) {
        for(size_t i = 0; i < vec.dim; ++i) {
            if(fabs(get_element_vec(vec, i)) > EPSILON) {
                return false;
            }
        }
        return true;
    } else {
        return false;
    }
}

bool is_null_vec(Vector vec)
{
    return equal_vec(vec, get_null_vec());
}

Vector transpose_vec(Vector vec)
{
    bool type;
    Vector t_vec = copy_vector(vec);

    if(vec.type == ROW) {
        type = COL;
    } else {
        type = ROW;
    }
    t_vec.type = type;

    return t_vec;
}

Vector add_vec(Vector lhs, Vector rhs)
{
    if(same_size_and_type(lhs, rhs) && !is_null_vec(lhs)) {
        Vector sum = make_vector(lhs.dim, lhs.type);

        for(size_t i = 0; i < lhs.dim; ++i) {
            double value = get_element_vec(lhs, i) + get_element_vec(rhs, i);
            set_element_vec(value, sum, i);
        }

        return sum;
    } else {
        return get_null_vec();
    }
}

Vector subtract_vec(Vector lhs, Vector rhs)
{
    Vector neg_rhs = smultiply_vec(rhs, -1);
    Vector difference = add_vec(lhs, neg_rhs);
    delete_vector(neg_rhs);

    return difference;
}

Vector smultiply_vec(Vector vec, double scalar)
{
    Vector product = copy_vector(vec);
    smultiply_inplace_vec(product, scalar);

    return product;
}

void smultiply_inplace_vec(Vector vec, double scalar)
{
    for(size_t i = 0; i < vec.dim; ++i) {
        double value = scalar * get_element_vec(vec, i);
        set_element_vec(value, vec, i);
    }
}

double dot(Vector lhs, Vector rhs)
{
    if(same_size_and_type(lhs, rhs)) {
        double product = 0;
        for(size_t i = 0; i < lhs.dim; ++i) {
            product += (get_element_vec(lhs, i) * get_element_vec(rhs, i));
        }
        return product;
    } else {
        return nan("");
    }
}

double dot_t(Vector lhs, Vector rhs)
{
    Vector t_lhs = transpose_vec(lhs);
    double product = dot(t_lhs, rhs);
    delete_vector(t_lhs);

    return product;
}

Vector change_ndim_basis(Vector *basis1, Vector *basis2, Vector vec, size_t dim)
{
    if(is_ndim_basis(basis1, dim) && is_ndim_basis(basis2, dim) && is_col(vec)) {
        Vector old_coords = get_coords(basis1, dim, vec);
        Vector new_coords = change_ndim_basis_coords(basis1, basis2, old_coords, dim);
        delete_vector(old_coords);

        if(!is_null_vec(new_coords)) {
            Vector new_vec = lin_combination(basis2, dim, new_coords);
            delete_vector(new_coords);

            return new_vec;
        } else {
            delete_vector(new_coords);
            return get_null_vec();
        }
    } else {
        return get_null_vec();
    }
}

Vector change_ndim_basis_coords(Vector *basis1, Vector *basis2, Vector coords, size_t dim)
{
    if(is_ndim_basis(basis1, dim) && is_ndim_basis(basis2, dim) && is_col(coords)) {
        Matrix t_mtx1 = make_matrix_cols(basis1, dim);

        Matrix t_mtx2_intermediate = make_matrix_cols(basis2, dim);
        Matrix t_mtx2 = get_inverse(t_mtx2_intermediate);
        delete_matrix(t_mtx2_intermediate);

        if(!is_null_mtx(t_mtx1) && !is_null_mtx(t_mtx2)) {
            Matrix transition_mtx = mmultiply(t_mtx2, t_mtx1);
            Vector new_coords = vmultiply(transition_mtx, coords, true);
            
            delete_matrix(t_mtx1);
            delete_matrix(t_mtx2);
            delete_matrix(transition_mtx);

            return new_coords;

        } else {
            delete_matrix(t_mtx1);
            delete_matrix(t_mtx2);
            return get_null_vec();
        }
    } else {
        return get_null_vec();
    }
}

Vector get_coords(Vector *basis, size_t dim, Vector vec)
{
    double *coords = get_coords_dbls(basis, dim, vec);

    if(coords) {
        return make_vector_vals(coords, dim, COL);
    } else {
        return get_null_vec();
    }
}

double * get_coords_dbls(Vector *basis, size_t dim, Vector vec)
{
    if(basis && !is_null_vec(vec) && dim > 0) {
        Matrix basis_mtx = make_matrix_cols(basis, dim);
        Vector sln = solve_system_v(basis_mtx, vec);
        double *sln_vals = get_elements_vec(sln);
        
        delete_matrix(basis_mtx);
        delete_vector(sln);

        return sln_vals;
    } else {
        return NULL;
    }
    
}

Vector lin_combination(Vector *basis, size_t dim, Vector coords)
{
    if(coords.dim == dim && dim > 0) {
        return lin_combination_dbls(basis, dim, coords.elements);
    } else {
        return get_null_vec();
    }
}

Vector lin_combination_dbls(Vector *basis, size_t dim, double *coords)
{
    if(basis && coords && dim > 0) {
        Vector lin_combo = make_vector(basis[0].dim, COL);
        for(size_t i = 0; i < dim; ++i) {
            Vector lin_combo_intermediate = copy_vector(lin_combo);
            Vector multiplier = smultiply_vec(basis[i], coords[i]);
            delete_vector(lin_combo);
            lin_combo = add_vec(lin_combo_intermediate, multiplier);

            delete_vector(lin_combo_intermediate);
            delete_vector(multiplier);
        }
        return lin_combo;
    } else {
        return get_null_vec();
    }
}

Vector * orthonormalize_set_vec(Vector *vecs, size_t dim)
{
    if(vecs && dim > 0) {
        Vector *orth_vecs = malloc(dim * sizeof(vecs[0]));
        if(orth_vecs) {        
            orth_vecs[0] = get_unit_vec(vecs[0]);

            for(size_t i = 1; i < dim; ++i) {
                Vector projection = get_projection_vec(orth_vecs, i, vecs[i]);
                Vector orth_vec = subtract_vec(vecs[i], projection);

                orth_vecs[i] = get_unit_vec(orth_vec);

                delete_vector(projection);
                delete_vector(orth_vec);
            }
        }
        return orth_vecs;
    } else {
        return NULL;
    }
}
