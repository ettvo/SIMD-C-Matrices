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

// TODO: start task 2
// program uses 256-bit instructions = 4 doubles
// work on 4 doubles at a time BUT need tail cases for 
// situations w/o 4 values at a time
// can go with the approach in store_simd4 ptr or use pointer arithmetic


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
    int index = mat->cols * row + col;
    return mat->data[index]; // need to check whether I need to subtract 1 from row and col
}

/*
 * Sets the value at the given row and column to val. You may assume `row` and
 * `col` are valid. Note that the matrix is in row-major order.
 */
void set(matrix *mat, int row, int col, double val) {
    // Task 1.1 TODO
    int index = mat->cols * row + col;
    mat->data[index] = val;
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
    //matrix_struct->data = (double*)malloc(rows * cols * sizeof(double));
    matrix_struct->data = (double*)calloc(rows * cols, sizeof(double));

    if (matrix_struct->data == NULL) {
        return -2;
    }

    matrix_struct->parent = NULL;
    matrix_struct->ref_cnt = 1;
    matrix_struct->rows = rows;
    matrix_struct->cols = cols;
    // mat[0] = matrix_struct; // use & ?
    *mat = matrix_struct;

    //for (int curr_row = 0; curr_row < rows; curr_row += 1) {
    //    for (int curr_col = 0; curr_col < cols; curr_col += 1) {
    //        set(matrix_struct, curr_row, curr_col, 0);
    //    }
    //}

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
        if (mat->ref_cnt == 0) {
            free(mat->data);
            free(mat);
        }
    } else {
        deallocate_matrix(mat->parent);
        // may have an issue if all data points to same thing --> don't free mat->data
        free(mat); // or just free mat?
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

    new_mat->data = (double*)(from->data + offset);

    // 3. Set the `data` field of the new struct to be the `data` field of the `from` struct plus `offset`.
    // double get(matrix *mat, int row, int col)
    //int curr_col = offset % cols;
    //int curr_row = (int)(offset / cols);

    // 4. Set the number of rows and columns in the new struct according to the arguments provided.
    new_mat->cols = cols;
    new_mat->rows = rows;
    new_mat->ref_cnt = 1; // what to put here?

    // 5. Set the `parent` field of the new struct to the `from` struct pointer.
    new_mat->parent = from;

    // 6. Increment the `ref_cnt` field of the `from` struct by 1.
    from->ref_cnt += 1;
    // 7. Store the address of the allocated matrix struct at the location `mat` is pointing at.
    // mat[0] = new_mat;
    *mat = new_mat;
    // 8. Return 0 upon success.
    return 0;
}

/*
 * Sets all entries in mat to val. Note that the matrix is in row-major order.
 */
void fill_matrix(matrix *mat, double val) {
    // Task 1.5 TODO
    /* non-optimized 
    for (int curr_row = 0; curr_row < mat->rows; curr_row += 1) {
        for (int curr_col = 0; curr_col < mat->cols; curr_col += 1) {
            set(mat, curr_row, curr_col, val);
        }
    }
    */
    //int counter_4 = (mat->rows * mat->cols)/4;
    //int counter_tail = (mat->rows * mat->cols) % 4;
    //double* ptr = (double*)malloc(sizeof(double*));

    /* Took time 0.1 sec
    #pragma omp parallel 
    {   
        double fill[4] = {val, val, val, val};
        __m256d temp = _mm256_loadu_pd(fill);
        for (int counter = 0; counter < ((mat->rows * mat->cols)/4); counter += 1) {
            // temp = _mm256_loadu_pd((__m256d*) (mat->data + 4 * counter));
            _mm256_storeu_pd((mat->data + 4 * counter), temp);
        }
    }

    #pragma omp parallel 
    {
        double fill[4] = {val, val, val, val};
        __m256d temp = _mm256_loadu_pd(fill);
        double tail_arr[4];
        _mm256_storeu_pd(tail_arr, temp);
        for (int counter = 0; counter < ((mat->rows * mat->cols) % 4); counter += 1) {
            mat->data[((mat->rows * mat->cols)/4) * 4 + counter] = tail_arr[counter];
        }
    }
    */
   /*
   double fill[4] = {val, val, val, val};
    __m256d temp = _mm256_loadu_pd(fill);

   // Took time 0.06 sec
   #pragma omp parallel for 
        for (int counter = 0; counter < ((mat->rows * mat->cols)/4); counter += 1) {
            // temp = _mm256_loadu_pd((__m256d*) (mat->data + 4 * counter));
            _mm256_storeu_pd((mat->data + 4 * counter), temp);
        }
    #pragma omp parallel for 
        for (int counter = 0; counter < ((mat->rows * mat->cols) % 4); counter += 1) {
            mat->data[((mat->rows * mat->cols)/4) * 4 + counter] = val;
        }
    */
   double fill[4] = {val, val, val, val};
    __m256d temp = _mm256_loadu_pd(fill);
    for (int counter = 0; counter < ((mat->rows * mat->cols)/4); counter += 1) {
        // temp = _mm256_loadu_pd((__m256d*) (mat->data + 4 * counter));
        _mm256_storeu_pd((mat->data + 4 * counter), temp);
    }
    for (int counter = 0; counter < ((mat->rows * mat->cols) % 4); counter += 1) {
        mat->data[((mat->rows * mat->cols)/4) * 4 + counter] = val;
    }
    //free(ptr);
}

/*
 * Store the result of taking the absolute value element-wise to `result`.
 * Return 0 upon success.
 * Note that the matrix is in row-major order.
 */
int abs_matrix(matrix *result, matrix *mat) {
    // Task 1.5 TODO
    /* non-optimized 
    double curr_val = -999;
    for (int curr_row = 0; curr_row < mat->rows; curr_row += 1) {
        for (int curr_col = 0; curr_col < mat->cols; curr_col += 1) {
            curr_val = fabs(get(mat, curr_row, curr_col));
            set(result, curr_row, curr_col, curr_val);
        }
    }
    */

    int counter_4 = (mat->rows * mat->cols)/4;
    int counter_tail = (mat->rows * mat->cols) % 4;
    // took 0.08 sec
        //double* ptr = (double*)malloc(sizeof(double*));
    __m256d temp;
    double neg[4] = {-1, -1, -1, -1};
    __m256d negator = _mm256_loadu_pd(neg);
    __m256d negated_arr;
    // floating point representation = has a sign bit
    #pragma omp parallel for
        for (int counter = 0; counter < counter_4; counter += 1) {
            //temp = _mm256_loadu_pd((__m256d*) (mat->data + 4 * counter));
            temp = _mm256_loadu_pd((mat->data + 4 * counter));
            negated_arr = _mm256_mul_pd(negator, temp);
            temp = _mm256_max_pd(negated_arr, temp);
            _mm256_storeu_pd((result->data + 4 * counter), temp);
        }

    temp = _mm256_loadu_pd((mat->data + 4 * counter_4));
    negated_arr = _mm256_mul_pd(negator, temp);
    temp = _mm256_max_pd(negated_arr, temp);
    double tail_arr[4];
    _mm256_storeu_pd(tail_arr, temp);
    //for (int counter = 0; counter < counter_tail; counter += 1) {
    //    result->data[counter_4 * 4 + counter] = tail_arr[counter];
    //}
    
    switch (counter_tail) {
        case 3:
            result->data[counter_4 * 4 + 2] = tail_arr[2];
        case 2:
            result->data[counter_4 * 4 + 1] = tail_arr[1];
        case 1:
            result->data[counter_4 * 4 + 0] = tail_arr[0];
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
    /* non-optimized 
    double curr_val = -999;
    for (int curr_row = 0; curr_row < mat->rows; curr_row += 1) {
        for (int curr_col = 0; curr_col < mat->cols; curr_col += 1) {
            curr_val = get(mat, curr_row, curr_col) * -1;
            set(result, curr_row, curr_col, curr_val);
        }
    }
    */
    int counter_4 = (mat->rows * mat->cols)/4;
    int counter_tail = (mat->rows * mat->cols) % 4;
    //double* ptr = (double*)malloc(sizeof(double*));
    double neg[4] = {-1, -1, -1, -1};
    __m256d temp;
    __m256d negator = _mm256_loadu_pd(neg);
    // floating point representation = has a sign bit
    for (int counter = 0; counter < counter_4; counter += 1) {
        temp = _mm256_loadu_pd((mat->data + 4 * counter));
        temp = _mm256_mul_pd(negator, temp);
        _mm256_storeu_pd((result->data + 4 * counter), temp);
    }
    
    temp = _mm256_loadu_pd((mat->data + 4 * counter_4));
    temp = _mm256_mul_pd(negator, temp);
    double tail_arr[4];
    _mm256_storeu_pd(tail_arr, temp);
    for (int counter = 0; counter < counter_tail; counter += 1) {
        result->data[counter_4 * 4 + counter] = tail_arr[counter];
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
    /* non-optimized
    double curr_mat1 = -999; 
    double curr_mat2 = -999; 
    for (int curr_row = 0; curr_row < result->rows; curr_row += 1) {
        for (int curr_col = 0; curr_col < result->cols; curr_col += 1) {
            curr_mat1 = get(mat1, curr_row, curr_col);
            curr_mat2 = get(mat2, curr_row, curr_col);
            set(result, curr_row, curr_col, curr_mat1 + curr_mat2);
        }
    }
    */
    int counter_4 = (result->rows * result->cols)/4;
    int counter_tail = (result->rows * result->cols) % 4;
    __m256d temp;
    __m256d mat1_load;
        //double* ptr = (double*)malloc(sizeof(double*));
        // floating point representation = has a sign bit
    #pragma omp parallel for
        for (int counter = 0; counter < counter_4; counter += 1) {
            mat1_load = _mm256_loadu_pd((mat1->data + 4 * counter));
            temp = _mm256_loadu_pd((mat2->data + 4 * counter));
            temp = _mm256_add_pd(temp, mat1_load);
            _mm256_storeu_pd((result->data + 4 * counter), temp);
        }
        //for (int counter = 0; counter < counter_tail; counter += 1) {
        //    result->data[counter_4 * 4 + counter] = mat1->data[counter_4 * 4 + counter] + mat2->data[counter_4 * 4 + counter];
        //}
    mat1_load = _mm256_loadu_pd((mat1->data + 4 * counter_4));
    temp = _mm256_loadu_pd((mat2->data + 4 * counter_4));
    temp = _mm256_add_pd(temp, mat1_load);
    double tail_arr[4];
    _mm256_storeu_pd(tail_arr, temp);
    //for (int counter = 0; counter < counter_tail; counter += 1) {
    //    result->data[counter_4 * 4 + counter] = tail_arr[counter];
    //}
    switch (counter_tail) {
        case 3:
            result->data[counter_4 * 4 + 2] = tail_arr[2];
        case 2:
            result->data[counter_4 * 4 + 1] = tail_arr[1];
        case 1:
            result->data[counter_4 * 4 + 0] = tail_arr[0];
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
    /* non-optimized 
    double curr_mat1 = -999; 
    double curr_mat2 = -999;
    for (int curr_row = 0; curr_row < result->rows; curr_row += 1) {
        for (int curr_col = 0; curr_col < result->cols; curr_col += 1) {
            curr_mat1 = get(mat1, curr_row, curr_col);
            curr_mat2 = get(mat2, curr_row, curr_col);
            set(result, curr_row, curr_col, curr_mat1 - curr_mat2);
        }
    }
    */
    int counter_4 = (result->rows * result->cols)/4;
    int counter_tail = (result->rows * result->cols) % 4;
    //double* ptr = (double*)malloc(sizeof(double*));
    __m256d temp;
    __m256d mat1_load;
    // floating point representation = has a sign bit
    for (int counter = 0; counter < counter_4; counter += 1) {
        mat1_load = _mm256_loadu_pd((mat1->data + 4 * counter));
        temp = _mm256_loadu_pd((mat2->data + 4 * counter));
        temp = _mm256_sub_pd(mat1_load, temp);
        _mm256_storeu_pd((result->data + 4 * counter), temp);
    }

    mat1_load = _mm256_loadu_pd((mat1->data + 4 * counter_4));
    temp = _mm256_loadu_pd((mat2->data + 4 * counter_4));
    temp = _mm256_sub_pd(mat1_load, temp);
    double tail_arr[4];
    _mm256_storeu_pd(tail_arr, temp);
    //for (int counter = 0; counter < counter_tail; counter += 1) {
    //    result->data[counter_4 * 4 + counter] = mat1->data[counter_4 * 4 + counter] + mat2->data[counter_4 * 4 + counter];
    //}
    for (int counter = 0; counter < counter_tail; counter += 1) {
        result->data[counter_4 * 4 + counter] = tail_arr[counter];
    }

    //for (int counter = 0; counter < counter_tail; counter += 1) {
    //    result->data[counter_4 * 4 + counter] = mat1->data[counter_4 * 4 + counter] - mat2->data[counter_4 * 4 + counter];
    //}
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
    if (mat1->cols != mat2->rows) {
        printf("Mul dimension error");
        return -1;
    }

    __m256d mat1_load;
    __m256d mat2_load;
    __m256d temp_total;
    double total;
    double multiplied[4];
    matrix *temp_result = NULL;
    allocate_matrix(&temp_result, result->rows, result->cols);

    matrix *transpose2 = NULL;
    allocate_matrix(&transpose2, mat2->cols, mat2->rows);
    transpose_matrix(transpose2, mat2);

    int counter_4 = (mat1->cols)/4;
    int counter_tail = (mat1->cols) % 4;
    // or (int counter = 0; counter < mat->cols * mat->rows; counter += 1) {
    //    curr_row = counter/(mat->cols);
    //    curr_col = counter - (curr_row * mat->cols);
    int curr_row;
    int curr_col;

    //#pragma omp parallel for

    #pragma omp parallel for
    for (int counter = 0; counter < result->rows * result->cols; counter += 1) {
        curr_row = counter/(result->cols);
        curr_col = counter - (curr_row * result->cols);
        total = 0;
        //#pragma omp parallel for reduction (+ : total)
        for (int counter = 0; counter < counter_4; counter += 1) {
            mat1_load = _mm256_loadu_pd((mat1->data + 4 * counter + curr_row * mat1->cols));
            mat2_load = _mm256_loadu_pd((transpose2->data + 4 * counter + curr_col * transpose2->cols)); 
            temp_total = _mm256_mul_pd(mat1_load, mat2_load);
            _mm256_storeu_pd(multiplied, temp_total);
            total += multiplied[0] + multiplied[1] + multiplied[2] + multiplied[3];
        }
        mat1_load = _mm256_loadu_pd((mat1->data + 4 * counter_4 + curr_row * mat1->cols));
        mat2_load = _mm256_loadu_pd((transpose2->data + 4 * counter_4 + curr_col * transpose2->cols));
            
        temp_total = _mm256_mul_pd(mat1_load, mat2_load);
        _mm256_storeu_pd(multiplied, temp_total);
             
        switch (counter_tail) {
        case 3:
            total += multiplied[2];
        case 2:
            total += multiplied[1];
        case 1:
            total += multiplied[0];
        }

        set(result, curr_row, curr_col, total);

    }


    /* no OpenMP
    for(int curr_row = 0; curr_row < result->rows; curr_row += 1) {
        for (int curr_col = 0; curr_col < result->cols; curr_col += 1) {
            total = 0;
            for (int counter = 0; counter < counter_4; counter += 1) {
                mat1_load = _mm256_loadu_pd((mat1->data + 4 * counter + curr_row * mat1->cols));
                mat2_load = _mm256_loadu_pd((transpose2->data + 4 * counter + curr_col * transpose2->cols)); 
                temp_total = _mm256_mul_pd(mat1_load, mat2_load);
                _mm256_storeu_pd(multiplied, temp_total);
                total += multiplied[0] + multiplied[1] + multiplied[2] + multiplied[3];
            }
            mat1_load = _mm256_loadu_pd((mat1->data + 4 * counter_4 + curr_row * mat1->cols));
            mat2_load = _mm256_loadu_pd((transpose2->data + 4 * counter_4 + curr_col * transpose2->cols));
            
            temp_total = _mm256_mul_pd(mat1_load, mat2_load);
            _mm256_storeu_pd(multiplied, temp_total);
            for (int counter = 0; counter < counter_tail; counter += 1) {
                total += multiplied[counter];
            }
            set(result, curr_row, curr_col, total);
        }
    }
    */

    deallocate_matrix(temp_result);
    deallocate_matrix(transpose2);

    /* non-optimized 
    double total = 0;

    matrix *temp = NULL;
    allocate_matrix(&temp, result->rows, result->cols);
    //fill_matrix(temp, 0);

    for(int curr_row = 0; curr_row < result->rows; curr_row += 1) { // TODO: fix the segmentation error for non-square matrices
        for (int curr_col = 0; curr_col < result->cols; curr_col += 1) { // segmentation error likely due to changes in this part
        // get total
            for (int counter = 0; counter < mat1->cols; counter += 1) {
                //printf("mat1 val: %f for row %d and col %d, mat2 val: %f for row %d and col %d\n", get(mat1, curr_row, counter), curr_row, counter, get(mat2, counter, curr_col), counter, curr_col);
                total += get(mat1, curr_row, counter) * get(mat2, counter, curr_col); // issue is with non-square dimensions
                // need to check that the given item has the given row or counter
            }
            // mat1 rows, mat2 cols, and mat1 cols == mat2 rows number sums per entry
            //printf("total %f at row %d, col %d\n", total, curr_row, curr_col);
            set(temp, curr_row, curr_col, total);
            total = 0;
        }
    }

    double curr_val;
    for (int curr_row = 0; curr_row < result->rows; curr_row += 1) {
        for (int curr_col = 0; curr_col < result->cols; curr_col += 1) {
            curr_val = get(temp, curr_row, curr_col);
            set(result, curr_row, curr_col, curr_val);
        }
    }
    deallocate_matrix(temp);
    */

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
        int counter = pow - 1;

        matrix *temp = NULL;
        allocate_matrix(&temp, mat->rows, mat->cols);
        mul_matrix(temp, mat, mat);

        while (counter > 2) {
            mul_matrix(result, temp, result);
            counter -= 2;
        }
        
        while (counter > 1) {
            mul_matrix(result, mat, result);
            counter -= 1;
        }
        deallocate_matrix(temp);
    }
    return 0;
}

/*** 
 * HELPER FUNCTIONS 
 ***/

/* 
Transposes a given MAT and stores it in RESULT.
*/ 
// can probably SIMD optimize this one
int transpose_matrix(matrix *result, matrix *mat) {
    if (result->cols != mat->rows || result->rows != mat->cols) {
        return -1;
    }
    double val;
    int curr_row = 0;
    int curr_col = 0;
    ///*
    //#pragma omp parallel for
    for (int counter = 0; counter < mat->cols * mat->rows; counter += 1) {
        curr_row = counter/(mat->cols);
        curr_col = counter - (curr_row * mat->cols);
        val = get(mat, curr_row, curr_col);
        set(result, curr_col, curr_row, val); 
    }
    //*/

    /*
    for (int mat_row = 0; mat_row < mat->rows; mat_row += 1) {
        for (int mat_col = 0; mat_col < mat->cols; mat_col += 1) {
            val = get(mat, mat_row, mat_col);
            set(result, mat_col, mat_row, val); 
        }
    }
    */
    return 0; 
}

/* Returns an array of 4 entries from the given matrix 
   starting at the given coordinates. If there are less than
   4 entries left, get_simd4 returns a matrix with the remaining
   entries left and zeroes in the remaining spaces. 
   Assumes row and col are valid coordinates for an entry in mat. */
double* get_simd4_ptr(matrix* mat, int row, int col) {
    int total_entries = (mat->rows + 1) * (mat->cols + 1);
    // if it doesn't work due to incorrect if-else, just add 1 
    if (total_entries - 4 > (row + 1) * (col+ 1)) { // at least 4 entries left
        int start_index = mat->cols * row + col;
        double result[4] = {mat->data[start_index], mat->data[start_index + 1], mat->data[start_index + 2], mat->data[start_index + 3]};
        double* ptr = (double*)malloc(sizeof(double*));
        ptr = result;
        return ptr;
    } else { // less than 4 entrries left
        double result[4] = {0, 0, 0, 0};
        int start_index = mat->cols * row + col;
        for (int counter = 0; counter < total_entries - row * col; counter += 1) {
            result[counter] = mat->data[start_index + counter];
        }
        double* ptr = (double*)malloc(sizeof(double*));
        ptr = result;
        return ptr;
    }
}

/* Stores a pointer to an array of 4 entries from the given matrix 
   starting at the given coordinates to DEST. If there are less than
   4 entries left, get_simd4 stores a pointer to a matrix with the remaining
   entries left and zeroes in the remaining spaces to DEST. 
   Assumes row and col are valid coordinates for an entry in mat. */
void store_simd4_ptr(double* dest, matrix* mat, int row, int col) {

    int total_entries = (mat->rows + 1) * (mat->cols + 1);
    // if it doesn't work due to incorrect if-else, just add 1 
    if (total_entries - 4 > (row + 1) * (col + 1)) { // at least 4 entries left
        int start_index = mat->cols * row + col;
        double result[4] = {mat->data[start_index], mat->data[start_index + 1], mat->data[start_index + 2], mat->data[start_index + 3]};
        dest = result;
    } else { // less than 4 entrries left
        double result[4] = {0, 0, 0, 0};
        int start_index = mat->cols * row + col;
        for (int counter = 0; counter < total_entries - row * col; counter += 1) {
            result[counter] = mat->data[start_index + counter];
        }
        dest = result;
    }
}


void check_mul_index(matrix *mat1, matrix *mat2) {
    __m256d mat1_load;
    __m256d mat2_load;
    
    double curr_vals[4];

    matrix *transpose2 = NULL;
    allocate_matrix(&transpose2, mat2->cols, mat2->rows);
    transpose_matrix(transpose2, mat2);

    int counter_4 = (mat1->cols)/4;
    //int counter_tail = (mat1->cols) % 4;

    //printf("\nApproach is (data + 4 * counter + 4 * curr_row * mat1->cols)\n");
    //printf("counter_4 = %d\n", counter_4);
    
    for(int curr_row = 0; curr_row < mat1->rows; curr_row += 1) {
        for (int curr_col = 0; curr_col < mat2->cols; curr_col += 1) {
            for (int counter = 0; counter < counter_4; counter += 1) {
                // counter should focus only on the # per row, not anything else 
                mat1_load = _mm256_loadu_pd((mat1->data + 4 * counter + 4 * curr_row * mat1->cols));
                _mm256_storeu_pd(curr_vals, mat1_load); 
                //printf("at row %d, col %d, counter %d, mat1 == [%f, %f, %f, %f]\n", curr_row, curr_col, counter, curr_vals[0], curr_vals[1], curr_vals[2], curr_vals[3]);
                mat2_load = _mm256_loadu_pd((transpose2->data + 4 * counter + 4 * curr_col * transpose2->cols));
                _mm256_storeu_pd(curr_vals, mat2_load); 
                //printf("at row %d, col %d, counter %d, mat2 == [%f, %f, %f, %f]\n", curr_row, curr_col, counter, curr_vals[0], curr_vals[1], curr_vals[2], curr_vals[3]);
            }
            
            mat1_load = _mm256_loadu_pd((mat1->data + curr_row * mat1->cols));
            mat2_load = _mm256_loadu_pd((transpose2->data + curr_col * transpose2->cols));

            _mm256_storeu_pd(curr_vals, mat1_load); 
            //printf("at row %d, col %d, counter %d, mat1 == [%f, %f, %f, %f]\n", curr_row, curr_col, counter_4, curr_vals[0], curr_vals[1], curr_vals[2], curr_vals[3]);
            _mm256_storeu_pd(curr_vals, mat2_load); 
            //printf("at row %d, col %d, counter %d, mat2 == [%f, %f, %f, %f]\n\n", curr_row, curr_col, counter_4, curr_vals[0], curr_vals[1], curr_vals[2], curr_vals[3]);
            
        }
    }

    deallocate_matrix(transpose2);
}