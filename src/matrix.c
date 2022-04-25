#include "matrix.h"
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <omp.h>

// Include SSE intrinsics
#if defined(_MSC_VER)
#include <intrin.h>
#elif defined(__GNUC__) && (defined(__x86_64__) || defined(__i386__))
#include <immintrin.h>
#include <x86intrin.h>
#endif

/* Below are some intel intrinsics that might be useful
 * void _mm256_storeu_pd (double * mem_addr, __m256d a)
 * __m256d _mm256_set1_pd (double a)
 * __m256d _mm256_set_pd (double e3, double e2, double e1, double e0)
 * __m256d _mm256_loadu_pd (double const * mem_addr)
 * __m256d _mm256_add_pd (__m256d a, __m256d b)
 * __m256d _mm256_sub_pd (__m256d a, __m256d b)
 * __m256d _mm256_fmadd_pd (__m256d a, __m256d b, __m256d c)
 * __m256d _mm256_mul_pd (__m256d a, __m256d b)
 * __m256d _mm256_cmp_pd (__m256d a, __m256d b, const int imm8)
 * __m256d _mm256_and_pd (__m256d a, __m256d b)
 * __m256d _mm256_max_pd (__m256d a, __m256d b)
*/

/* Generates a random double between low and high */
double rand_double(double low, double high) {
    double range = (high - low);
    double div = RAND_MAX / range;
    return low + (rand() / div);
}

/* Generates a random matrix */
void rand_matrix(matrix *result, unsigned int seed, double low, double high) {
    srand(seed);
    for (int i = 0; i < result->rows; i++) {
        for (int j = 0; j < result->cols; j++) {
            set(result, i, j, rand_double(low, high));
        }
    }
}

/*
 * Returns the double value of the matrix at the given row and column.
 * You may assume `row` and `col` are valid. Note that the matrix is in row-major order.
 */
double get(matrix *mat, int row, int col) {
    // Task 1.1 TODO
    return mat->data[row*col + col]; // need to check whether I need to subtract 1 from row and col
}

/*
 * Sets the value at the given row and column to val. You may assume `row` and
 * `col` are valid. Note that the matrix is in row-major order.
 */
void set(matrix *mat, int row, int col, double val) {
    // Task 1.1 TODO
    mat->data[row*col + col] = val;
}

/*
 * Allocates space for a matrix struct pointed to by the double pointer mat with
 * `rows` rows and `cols` columns. You should also allocate memory for the data array
 * and initialize all entries to be zeros. `parent` should be set to NULL to indicate that
 * this matrix is not a slice. You should also set `ref_cnt` to 1.
 * You should return -1 if either `rows` or `cols` or both have invalid values. Return -2 if any
 * call to allocate memory in this function fails.
 * Return 0 upon success.
 */
int allocate_matrix(matrix **mat, int rows, int cols) {
    // Task 1.2 TODO
    // HINTS: Follow these steps.
    // 1. Check if the dimensions are valid. Return -1 if either dimension is not positive. (check)
    // 2. Allocate space for the new matrix struct. Return -2 if allocating memory failed. (check)
    // 3. Allocate space for the matrix data, initializing all entries to be 0. Return -2 if allocating memory failed.  (check)
    // 4. Set the number of rows and columns in the matrix struct according to the arguments provided. (check)
    // 5. Set the `parent` field to NULL, since this matrix was not created from a slice. (check)
    // 6. Set the `ref_cnt` field to 1. (check)
    // 7. Store the address of the allocated matrix struct at the location `mat` is pointing at.
    // 8. Return 0 upon success.

    if (rows <= 0 || cols <= 0) {
        return -1;
    }
    matrix* matrix_struct = (matrix*)malloc(sizeof(matrix));
    if (matrix_struct == NULL) {
        return -2;
    }

    // stored in row major fashion
    // matrix_struct->data = (double*)calloc(rows * cols, sizeof(double)); // does it need +1 for the terminator?
    matrix_struct->data = (double*)malloc(rows * cols * sizeof(double));

    if (matrix_struct->data == NULL) {
        return -2;
    }

    matrix_struct->parent = NULL;
    matrix_struct->ref_cnt = 1;
    matrix_struct->rows = rows;
    matrix_struct->cols = cols;
    mat[0] = matrix_struct; // use & ?

    for (int curr_row = 0; curr_row < rows; curr_row += 1) {
        for (int curr_col = 0; curr_col < cols; curr_col += 1) {
            set(matrix_struct, curr_row, curr_col, 0);
        }
    }

    return 0;

}

/*
 * You need to make sure that you only free `mat->data` if `mat` is not a slice and has no existing slices,
 * or that you free `mat->parent->data` if `mat` is the last existing slice of its parent matrix and its parent
 * matrix has no other references (including itself).
 */
void deallocate_matrix(matrix *mat) {
    // Task 1.3 TODO
    // HINTS: Follow these steps.
    // 1. If the matrix pointer `mat` is NULL, return.
    // 2. If `mat` has no parent: decrement its `ref_cnt` field by 1. If the `ref_cnt` field becomes 0, then free `mat` and its `data` field.
    // 3. Otherwise, recursively call `deallocate_matrix` on `mat`'s parent, then free `mat`.
    if (mat == NULL) {
        return;
    }
    if (mat->parent == NULL) {
        mat->ref_cnt -= 1;
    }
    if (mat->ref_cnt == 0) {
        free(mat->data);
        free(mat);
    } else {
        deallocate_matrix(mat->parent);
        // may have an issue if all data points to same thing --> don't free mat->data
        deallocate_matrix(mat); // or just free mat?
    }
}

/*
 * Allocates space for a matrix struct pointed to by `mat` with `rows` rows and `cols` columns.
 * Its data should point to the `offset`th entry of `from`'s data (you do not need to allocate memory)
 * for the data field. `parent` should be set to `from` to indicate this matrix is a slice of `from`
 * and the reference counter for `from` should be incremented. Lastly, do not forget to set the
 * matrix's row and column values as well.
 * You should return -1 if either `rows` or `cols` or both have invalid values. Return -2 if any
 * call to allocate memory in this function fails.
 * Return 0 upon success.
 * NOTE: Here we're allocating a matrix struct that refers to already allocated data, so
 * there is no need to allocate space for matrix data.
 */
int allocate_matrix_ref(matrix **mat, matrix *from, int offset, int rows, int cols) {
    // Task 1.4 TODO
    // HINTS: Follow these steps.
    // 1. Check if the dimensions are valid. Return -1 if either dimension is not positive.
    if (rows <= 0 || cols <= 0) {
        return -1;
    }
    // 2. Allocate space for the new matrix struct. Return -2 if allocating memory failed.
    matrix* new_mat = (matrix*)malloc(sizeof(matrix)); // might not want to calloc
    if (new_mat == NULL) {
        return -2;
    }

    new_mat->data = (double*)malloc(sizeof(double) * rows * cols);
    if (new_mat->data == NULL) {
        return -2;
    }
    // new_mat->data = from->data[offset..from->rows * from->cols];

    // 3. Set the `data` field of the new struct to be the `data` field of the `from` struct plus `offset`.
    // double get(matrix *mat, int row, int col)
    //int curr_col = offset % cols;
    //int curr_row = (int)(offset / cols);

    for (int counter = 0; counter < rows * cols; counter += 1) {
        double old = from->data[counter + offset];
        new_mat->data[counter] = old;
    }

    // 4. Set the number of rows and columns in the new struct according to the arguments provided.
    new_mat->cols = cols;
    new_mat->rows = rows;

    // 5. Set the `parent` field of the new struct to the `from` struct pointer.
    new_mat->parent = from;

    // 6. Increment the `ref_cnt` field of the `from` struct by 1.
    from->ref_cnt += 1;
    // 7. Store the address of the allocated matrix struct at the location `mat` is pointing at.
    mat[0] = new_mat;
    // 8. Return 0 upon success.
    return 0;
}

/*
 * Sets all entries in mat to val. Note that the matrix is in row-major order.
 */
void fill_matrix(matrix *mat, double val) {
    // Task 1.5 TODO
    for (int curr_row = 0; curr_row < mat->rows; curr_row += 1) {
        for (int curr_col = 0; curr_col < mat->cols; curr_col += 1) {
            set(mat, curr_row, curr_col, val);
        }
    }
}

/*
 * Store the result of taking the absolute value element-wise to `result`.
 * Return 0 upon success.
 * Note that the matrix is in row-major order.
 */
int abs_matrix(matrix *result, matrix *mat) {
    // Task 1.5 TODO
    double curr_val = -999;
    for (int curr_row = 0; curr_row < mat->rows; curr_row += 1) {
        for (int curr_col = 0; curr_col < mat->cols; curr_col += 1) {
            curr_val = fabs(get(mat, curr_row, curr_col));
            set(result, curr_row, curr_col, curr_val);
        }
    }
    return 0;
}

/*
 * (OPTIONAL)
 * Store the result of element-wise negating mat's entries to `result`.
 * Return 0 upon success.
 * Note that the matrix is in row-major order.
 */
int neg_matrix(matrix *result, matrix *mat) {
    // Task 1.5 TODO
    double curr_val = -999;
    for (int curr_row = 0; curr_row < mat->rows; curr_row += 1) {
        for (int curr_col = 0; curr_col < mat->cols; curr_col += 1) {
            curr_val = get(mat, curr_row, curr_col) * -1;
            set(result, curr_row, curr_col, curr_val);
        }
    }
    return 0;
}

/*
 * Store the result of adding mat1 and mat2 to `result`.
 * Return 0 upon success.
 * You may assume `mat1` and `mat2` have the same dimensions.
 * Note that the matrix is in row-major order.
 */
int add_matrix(matrix *result, matrix *mat1, matrix *mat2) {
    // Task 1.5 TODO
    double curr_mat1 = -999; 
    double curr_mat2 = -999;
    for (int curr_row = 0; curr_row < result->rows; curr_row += 1) {
        for (int curr_col = 0; curr_col < result->cols; curr_col += 1) {
            curr_mat1 = get(mat1, curr_row, curr_col);
            curr_mat2 = get(mat2, curr_row, curr_col);
            set(result, curr_row, curr_col, curr_mat1 + curr_mat2);
        }
    }
    return 0;
}

/*
 * (OPTIONAL)
 * Store the result of subtracting mat2 from mat1 to `result`.
 * Return 0 upon success.
 * You may assume `mat1` and `mat2` have the same dimensions.
 * Note that the matrix is in row-major order.
 */
int sub_matrix(matrix *result, matrix *mat1, matrix *mat2) {
    // Task 1.5 TODO
    double curr_mat1 = -999; 
    double curr_mat2 = -999;
    for (int curr_row = 0; curr_row < result->rows; curr_row += 1) {
        for (int curr_col = 0; curr_col < result->cols; curr_col += 1) {
            curr_mat1 = get(mat1, curr_row, curr_col);
            curr_mat2 = get(mat2, curr_row, curr_col);
            set(result, curr_row, curr_col, curr_mat1 - curr_mat2);
        }
    }
    return 0;
}

/*
 * Store the result of multiplying mat1 and mat2 to `result`.
 * Return 0 upon success.
 * Remember that matrix multiplication is not the same as multiplying individual elements.
 * You may assume `mat1`'s number of columns is equal to `mat2`'s number of rows.
 * Note that the matrix is in row-major order.
 */
int mul_matrix(matrix *result, matrix *mat1, matrix *mat2) {
    // Task 1.6 TODO
    // matrix multiplication --> mat1 col == mat2 row
    // new output = mat1 row x mat2 col
    // each entry in result = row of mat 1 * col of mat 2
    double total = 0;
    int curr_col = 0;
    int curr_row = 0;

    while(curr_row < result->rows) {
        if (curr_col == result->cols) {
            curr_col = 0;
            curr_row += 1;
        }

        // get total
        for (int counter = 0; counter < result->rows; counter += 1) {
            total += get(mat1, curr_row, counter) * get(mat2, counter, curr_col);
        }

        set(result, curr_row, curr_col, total);
        curr_col += 1;
    }
    return 0;
}

/*
 * Store the result of raising mat to the (pow)th power to `result`.
 * Return 0 upon success.
 * Remember that pow is defined with matrix multiplication, not element-wise multiplication.
 * You may assume `mat` is a square matrix and `pow` is a non-negative integer.
 * Note that the matrix is in row-major order.
 */
int pow_matrix(matrix *result, matrix *mat, int pow) {
    // Task 1.6 TODO
    // if assigning, may need to malloc, realloc space for each 
    // need to determine if dimensions of result are good

    if (pow == 0) {
        for (int curr_row = 0; curr_row < result->rows; curr_row += 1) {
            for (int curr_col = 0; curr_col < result->cols; curr_col += 1) {
                if (curr_row == curr_col) {
                    set(result, curr_row, curr_col, 1);
                } else {
                    set(result, curr_row, curr_col, 0);
                }
            }
        }
    } else if (pow == 1) {
        for (int curr_row = 0; curr_row < result->rows; curr_row += 1) {
            for (int curr_col = 0; curr_col < result->cols; curr_col += 1) {
                set(result, curr_row, curr_col, get(mat, curr_row, curr_col));
            }
        }
    } else {
        mul_matrix(result, mat, mat);
        pow_matrix(result, result, pow - 1);
    }
    return 0;
}

