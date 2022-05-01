#include "CUnit/CUnit.h"
#include "CUnit/Basic.h"
#include "../src/matrix.h"
#include <stdio.h>

/* Test Suite setup and cleanup functions: */
int init_suite(void) { return 0; }
int clean_suite(void) { return 0; }

/************* Test case functions ****************/

void get_test(void) {
  matrix *mat = NULL;
  allocate_matrix(&mat, 2, 2);
  for (int i = 0; i < 2; i++) {
    for (int j = 0; j < 2; j++) {
      set(mat, i, j, i * 2 + j);
    }
  }
  CU_ASSERT_EQUAL(get(mat, 0, 0), 0);
  CU_ASSERT_EQUAL(get(mat, 0, 1), 1);
  CU_ASSERT_EQUAL(get(mat, 1, 0), 2);
  CU_ASSERT_EQUAL(get(mat, 1, 1), 3);
  deallocate_matrix(mat);
}

void set_test(void) {
  matrix *mat = NULL;
  allocate_matrix(&mat, 2, 2);
  for (int i = 0; i < 2; i++) {
    for (int j = 0; j < 2; j++) {
      set(mat, i, j, i * 2 + j);
    }
  }
  CU_ASSERT_EQUAL(get(mat, 0, 0), 0);
  CU_ASSERT_EQUAL(get(mat, 0, 1), 1);
  CU_ASSERT_EQUAL(get(mat, 1, 0), 2);
  CU_ASSERT_EQUAL(get(mat, 1, 1), 3);
  deallocate_matrix(mat);
}

void alloc_fail_test(void) {
  matrix *mat = NULL;
  CU_ASSERT_EQUAL(allocate_matrix(&mat, 0, 0), -1);
  CU_ASSERT_EQUAL(allocate_matrix(&mat, 0, 1), -1);
  CU_ASSERT_EQUAL(allocate_matrix(&mat, 1, 0), -1);
}

void alloc_success_test(void) {
  matrix *mat = NULL;
  CU_ASSERT_EQUAL(allocate_matrix(&mat, 3, 2), 0);
  CU_ASSERT_EQUAL(mat->parent, NULL);
  CU_ASSERT_EQUAL(mat->ref_cnt, 1);
  CU_ASSERT_EQUAL(mat->rows, 3);
  CU_ASSERT_EQUAL(mat->cols, 2);
  CU_ASSERT_NOT_EQUAL(mat->data, NULL);
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 2; j++) {
      CU_ASSERT_EQUAL(get(mat, i, j), 0);
    }
  }
  deallocate_matrix(mat);
}

void alloc_ref_fail_test(void) {
  matrix *mat = NULL;
  matrix *from = NULL;
  CU_ASSERT_EQUAL(allocate_matrix_ref(&mat, from, 0, 0, 0), -1);
  CU_ASSERT_EQUAL(allocate_matrix_ref(&mat, from, 0, 0, 1), -1);
  CU_ASSERT_EQUAL(allocate_matrix_ref(&mat, from, 0, 1, 0), -1);
}

void alloc_ref_success_test(void) {
  matrix *mat = NULL;
  matrix *from = NULL;
  allocate_matrix(&from, 3, 2);
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 2; j++) {
      set(from, i, j, i * 2 + j);
    }
  }
  CU_ASSERT_EQUAL(allocate_matrix_ref(&mat, from, 2, 2, 2), 0);
  CU_ASSERT_PTR_EQUAL(mat->data, from->data + 2);
  CU_ASSERT_PTR_EQUAL(mat->parent, from);
  CU_ASSERT_EQUAL(mat->parent->ref_cnt, 2);
  CU_ASSERT_EQUAL(mat->rows, 2);
  CU_ASSERT_EQUAL(mat->cols, 2);
  deallocate_matrix(from);
  deallocate_matrix(mat);
}

void dealloc_null_test(void) {
  matrix *mat = NULL;
  deallocate_matrix(mat); // Test the null case doesn't crash
}

void add_test(void) {
  matrix *result = NULL;
  matrix *mat1 = NULL;
  matrix *mat2 = NULL;
  CU_ASSERT_EQUAL(allocate_matrix(&result, 2, 2), 0);
  CU_ASSERT_EQUAL(allocate_matrix(&mat1, 2, 2), 0);
  CU_ASSERT_EQUAL(allocate_matrix(&mat2, 2, 2), 0);
  for (int i = 0; i < 2; i++) {
    for (int j = 0; j < 2; j++) {
      set(mat1, i, j, i * 2 + j);
      set(mat2, i, j, i * 2 + j);
    }
  }
  add_matrix(result, mat1, mat2);
  for (int i = 0; i < 2; i++) {
    for (int j = 0; j < 2; j++) {
      CU_ASSERT_EQUAL(get(result, i, j), 2 * (i * 2 + j));
    }
  }
  deallocate_matrix(result);
  deallocate_matrix(mat1);
  deallocate_matrix(mat2);
}

/* (OPTIONAL) Uncomment the following sub_test if you have decided to implement it in matrix.c. */
void sub_test(void) {
  matrix *result = NULL;
  matrix *mat1 = NULL;
  matrix *mat2 = NULL;
  CU_ASSERT_EQUAL(allocate_matrix(&result, 2, 2), 0);
  CU_ASSERT_EQUAL(allocate_matrix(&mat1, 2, 2), 0);
  CU_ASSERT_EQUAL(allocate_matrix(&mat2, 2, 2), 0);
  for (int i = 0; i < 2; i++) {
    for (int j = 0; j < 2; j++) {
      set(mat1, i, j, i * 2 + j);
      set(mat2, i, j, (i * 2 + j) * 3);
    }
  }
  sub_matrix(result, mat1, mat2);
  for (int i = 0; i < 2; i++) {
    for (int j = 0; j < 2; j++) {
      CU_ASSERT_EQUAL(get(result, i, j), (-2) * (i * 2 + j));
    }
  }
  deallocate_matrix(result);
  deallocate_matrix(mat1);
  deallocate_matrix(mat2);
}


/* (OPTIONAL) Uncomment the following neg_test if you have decided to implement it in matrix.c. */
void neg_test(void) {
  matrix *result = NULL;
  matrix *mat = NULL;
  CU_ASSERT_EQUAL(allocate_matrix(&result, 2, 2), 0);
  CU_ASSERT_EQUAL(allocate_matrix(&mat, 2, 2), 0);
  for (int i = 0; i < 2; i++) {
    for (int j = 0; j < 2; j++) {
      set(mat, i, j, i * 2 + j);
    }
  }
  neg_matrix(result, mat);
  for (int i = 0; i < 2; i++) {
    for (int j = 0; j < 2; j++) {
      CU_ASSERT_EQUAL(get(result, i, j), -(i * 2 + j));
    }
  }
  deallocate_matrix(result);
  deallocate_matrix(mat);
} 

void abs_test(void) {
  matrix *result = NULL;
  matrix *mat = NULL;
  CU_ASSERT_EQUAL(allocate_matrix(&result, 2, 2), 0);
  CU_ASSERT_EQUAL(allocate_matrix(&mat, 2, 2), 0);
  for (int i = 0; i < 2; i++) {
    for (int j = 0; j < 2; j++) {
      if (j % 2 == 0)
        set(mat, i, j, i * 2 + j);
      else
        set(mat, i, j, -(i * 2 + j));
    }
  }
  abs_matrix(result, mat);
  for (int i = 0; i < 2; i++) {
    for (int j = 0; j < 2; j++) {
      CU_ASSERT_EQUAL(get(result, i, j), i * 2 + j);
    }
  }
  deallocate_matrix(result);
  deallocate_matrix(mat);
}

void mul_square_test(void) {
  matrix *result = NULL;
  matrix *mat1 = NULL;
  matrix *mat2 = NULL;
  CU_ASSERT_EQUAL(allocate_matrix(&result, 3, 3), 0);
  CU_ASSERT_EQUAL(allocate_matrix(&mat1, 3, 3), 0);
  CU_ASSERT_EQUAL(allocate_matrix(&mat2, 3, 3), 0);
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      set(mat1, i, j, i * 3 + j + 1);
      set(mat2, i, j, i * 3 + j + 1);
    }
  }
  mul_matrix(result, mat1, mat2);
  CU_ASSERT_EQUAL(get(result, 0, 0), 30);
  printf("Wanted %d at row %d, col %d, but got %f\n", 36, 0, 1, get(result, 0, 1));
  CU_ASSERT_EQUAL(get(result, 0, 1), 36);

  printf("Wanted %d at row %d, col %d, but got %f\n", 42, 0, 2, get(result, 0, 2));
  CU_ASSERT_EQUAL(get(result, 0, 2), 42);
  CU_ASSERT_EQUAL(get(result, 1, 0), 66);
  CU_ASSERT_EQUAL(get(result, 1, 1), 81);
  CU_ASSERT_EQUAL(get(result, 1, 2), 96);
  CU_ASSERT_EQUAL(get(result, 2, 0), 102);
  CU_ASSERT_EQUAL(get(result, 2, 1), 126);
  CU_ASSERT_EQUAL(get(result, 2, 2), 150);
  deallocate_matrix(result);
  deallocate_matrix(mat1);
  deallocate_matrix(mat2);
}

void mul_non_square_test(void) {
  matrix *result = NULL;
  matrix *mat1 = NULL;
  matrix *mat2 = NULL;
  CU_ASSERT_EQUAL(allocate_matrix(&result, 3, 3), 0);
  CU_ASSERT_EQUAL(allocate_matrix(&mat1, 3, 2), 0);
  CU_ASSERT_EQUAL(allocate_matrix(&mat2, 2, 3), 0);
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 2; j++) {
      set(mat1, i, j, i * 2 + j + 1);
    }
  }
  for (int i = 0; i < 2; i++) {
    for (int j = 0; j < 3; j++) {
      set(mat2, i, j, i * 3 + j + 1);
    }
  }
  mul_matrix(result, mat1, mat2);
  CU_ASSERT_EQUAL(get(result, 0, 0), 9);
  CU_ASSERT_EQUAL(get(result, 0, 1), 12);
  CU_ASSERT_EQUAL(get(result, 0, 2), 15);
  CU_ASSERT_EQUAL(get(result, 1, 0), 19);
  CU_ASSERT_EQUAL(get(result, 1, 1), 26);
  CU_ASSERT_EQUAL(get(result, 1, 2), 33);
  CU_ASSERT_EQUAL(get(result, 2, 0), 29);
  CU_ASSERT_EQUAL(get(result, 2, 1), 40);
  CU_ASSERT_EQUAL(get(result, 2, 2), 51);
  deallocate_matrix(result);
  deallocate_matrix(mat1);
  deallocate_matrix(mat2);
}

void pow_test(void) {
  matrix *result = NULL;
  matrix *mat = NULL;
  CU_ASSERT_EQUAL(allocate_matrix(&result, 2, 2), 0);
  CU_ASSERT_EQUAL(allocate_matrix(&mat, 2, 2), 0);
  set(mat, 0, 0, 1);
  set(mat, 0, 1, 1);
  set(mat, 1, 0, 1);
  set(mat, 1, 1, 0);
  pow_matrix(result, mat, 3);
  //printf("Expected %d but got %f at row %d, col %d\n", 3, get(result, 0, 0), 0, 0);
  CU_ASSERT_EQUAL(get(result, 0, 0), 3);
  //printf("Expected %d but got %f at row %d, col %d\n", 2, get(result, 0, 1), 0, 1);
  CU_ASSERT_EQUAL(get(result, 0, 1), 2);
  //printf("Expected %d but got %f at row %d, col %d\n", 2, get(result, 1, 0), 1, 0);
  CU_ASSERT_EQUAL(get(result, 1, 0), 2);
  CU_ASSERT_EQUAL(get(result, 1, 1), 1);
  pow_matrix(result, mat, 10);
  CU_ASSERT_EQUAL(get(result, 0, 0), 89);
  CU_ASSERT_EQUAL(get(result, 0, 1), 55);
  CU_ASSERT_EQUAL(get(result, 1, 0), 55);
  CU_ASSERT_EQUAL(get(result, 1, 1), 34);
  deallocate_matrix(result);
  deallocate_matrix(mat);
}


void transpose_test(void) {
  matrix *mat1 = NULL; // small, different dimensions
  
  CU_ASSERT_EQUAL(allocate_matrix(&mat1, 1, 3), 0);

  set(mat1, 0, 0, 1);
  set(mat1, 0, 1, 2);
  set(mat1, 0, 2, 3);
  
  matrix *temp_result = NULL;
  allocate_matrix(&temp_result, mat1->cols, mat1->rows);
  fill_matrix(temp_result, 0);
  transpose_matrix(temp_result, mat1);

  CU_ASSERT_EQUAL(get(temp_result, 0, 0), 1);
  CU_ASSERT_EQUAL(get(temp_result, 1, 0), 2);
  CU_ASSERT_EQUAL(get(temp_result, 2, 0), 3);

  deallocate_matrix(temp_result);
  deallocate_matrix(mat1);


  matrix *mat2 = NULL; // medium, different dimensions

  CU_ASSERT_EQUAL(allocate_matrix(&mat2, 3, 5), 0);
  
  set(mat2, 0, 0, 1);
  set(mat2, 0, 1, 2);
  set(mat2, 0, 2, 3);
  set(mat2, 0, 3, 4);
  set(mat2, 0, 4, 5);

  set(mat2, 1, 0, 6);
  set(mat2, 1, 1, 7);
  set(mat2, 1, 2, 8);
  set(mat2, 1, 3, 9);
  set(mat2, 1, 4, 10);

  set(mat2, 2, 0, 11);
  set(mat2, 2, 1, 12);
  set(mat2, 2, 2, 13);
  set(mat2, 2, 3, 14);
  set(mat2, 2, 4, 15);

  temp_result = NULL;

  allocate_matrix(&temp_result, mat2->cols, mat2->rows);
  fill_matrix(temp_result, 0);
  transpose_matrix(temp_result, mat2);

  CU_ASSERT_EQUAL(get(temp_result, 0, 0), 1);
  CU_ASSERT_EQUAL(get(temp_result, 1, 0), 2);
  CU_ASSERT_EQUAL(get(temp_result, 2, 0), 3);
  CU_ASSERT_EQUAL(get(temp_result, 3, 0), 4);
  CU_ASSERT_EQUAL(get(temp_result, 4, 0), 5);

  CU_ASSERT_EQUAL(get(temp_result, 0, 1), 6);
  CU_ASSERT_EQUAL(get(temp_result, 1, 1), 7);
  CU_ASSERT_EQUAL(get(temp_result, 2, 1), 8);
  CU_ASSERT_EQUAL(get(temp_result, 3, 1), 9);
  CU_ASSERT_EQUAL(get(temp_result, 4, 1), 10);

  CU_ASSERT_EQUAL(get(temp_result, 0, 2), 11);
  CU_ASSERT_EQUAL(get(temp_result, 1, 2), 12);
  CU_ASSERT_EQUAL(get(temp_result, 2, 2), 13);
  CU_ASSERT_EQUAL(get(temp_result, 3, 2), 14);
  CU_ASSERT_EQUAL(get(temp_result, 4, 2), 15);
  
  deallocate_matrix(temp_result);
  deallocate_matrix(mat2);

  
  matrix *mat3 = NULL; // same dimensions
  CU_ASSERT_EQUAL(allocate_matrix(&mat3, 4, 4), 0);

  set(mat3, 0, 0, 1);
  set(mat3, 0, 1, 2);
  set(mat3, 0, 2, 3);
  set(mat3, 0, 3, 4);

  set(mat3, 1, 0, 5);
  set(mat3, 1, 1, 6);
  set(mat3, 1, 2, 7);
  set(mat3, 1, 3, 8);

  set(mat3, 2, 0, 9);
  set(mat3, 2, 1, 10);
  set(mat3, 2, 2, 11);
  set(mat3, 2, 3, 12);
  
  set(mat3, 3, 0, 13);
  set(mat3, 3, 1, 14);
  set(mat3, 3, 2, 15);
  set(mat3, 3, 3, 16);

  temp_result = NULL;

  allocate_matrix(&temp_result, mat3->cols, mat3->rows);
  fill_matrix(temp_result, 0);
  transpose_matrix(temp_result, mat3);

  CU_ASSERT_EQUAL(get(temp_result, 0, 0), 1);
  CU_ASSERT_EQUAL(get(temp_result, 1, 0), 2);
  CU_ASSERT_EQUAL(get(temp_result, 2, 0), 3);
  CU_ASSERT_EQUAL(get(temp_result, 3, 0), 4);

  CU_ASSERT_EQUAL(get(temp_result, 0, 1), 5);
  CU_ASSERT_EQUAL(get(temp_result, 1, 1), 6);
  CU_ASSERT_EQUAL(get(temp_result, 2, 1), 7);
  CU_ASSERT_EQUAL(get(temp_result, 3, 1), 8);

  CU_ASSERT_EQUAL(get(temp_result, 0, 2), 9);
  CU_ASSERT_EQUAL(get(temp_result, 1, 2), 10);
  CU_ASSERT_EQUAL(get(temp_result, 2, 2), 11);
  CU_ASSERT_EQUAL(get(temp_result, 3, 2), 12);

  CU_ASSERT_EQUAL(get(temp_result, 0, 3), 13);
  CU_ASSERT_EQUAL(get(temp_result, 1, 3), 14);
  CU_ASSERT_EQUAL(get(temp_result, 2, 3), 15);
  CU_ASSERT_EQUAL(get(temp_result, 3, 3), 16);
  
  deallocate_matrix(temp_result);
  deallocate_matrix(mat3);

}

void mul_comp_test1(void) {
  matrix *result = NULL;
  matrix *mat1 = NULL;
  matrix *mat2 = NULL;

  CU_ASSERT_EQUAL(allocate_matrix(&result, 2, 2), 0);
  CU_ASSERT_EQUAL(allocate_matrix(&mat1, 2, 2), 0);
  CU_ASSERT_EQUAL(allocate_matrix(&mat2, 2, 2), 0);

  set(mat1, 0, 0, 3);
  set(mat1, 0, 1, 7);
  set(mat1, 1, 0, 4);
  set(mat1, 1, 1, 9);

  set(mat2, 0, 0, 6);
  set(mat2, 0, 1, 2);
  set(mat2, 1, 0, 5);
  set(mat2, 1, 1, 8);

  mul_matrix(result, mat1, mat2);
  printf("expected %d at row %d, col %d, but got %f\n", 53, 0, 0, get(result, 0, 0));
  CU_ASSERT_EQUAL(get(result, 0, 0), 53);
  CU_ASSERT_EQUAL(get(result, 0, 1), 62);
  CU_ASSERT_EQUAL(get(result, 1, 0), 69);
  CU_ASSERT_EQUAL(get(result, 1, 1), 80);

  deallocate_matrix(result);
  deallocate_matrix(mat1);
  deallocate_matrix(mat2);
  
  result = NULL;
  mat1 = NULL;
  mat2 = NULL;

  CU_ASSERT_EQUAL(allocate_matrix(&result, 3, 3), 0);
  CU_ASSERT_EQUAL(allocate_matrix(&mat1, 3, 3), 0);
  CU_ASSERT_EQUAL(allocate_matrix(&mat2, 3, 3), 0);

  set(mat1, 0, 0, 12);
  set(mat1, 0, 1, 8);
  set(mat1, 0, 2, 4);

  set(mat1, 1, 0, 3);
  set(mat1, 1, 1, 17);
  set(mat1, 1, 2, 14);

  set(mat1, 2, 0, 9);
  set(mat1, 2, 1, 8);
  set(mat1, 2, 2, 10);


  set(mat2, 0, 0, 5);
  set(mat2, 0, 1, 19);
  set(mat2, 0, 2, 3);

  set(mat2, 1, 0, 6);
  set(mat2, 1, 1, 15);
  set(mat2, 1, 2, 9);

  set(mat2, 2, 0, 7);
  set(mat2, 2, 1, 8);
  set(mat2, 2, 2, 16);

  mul_matrix(result, mat1, mat2);
  CU_ASSERT_EQUAL(get(result, 0, 0), 136);
  CU_ASSERT_EQUAL(get(result, 0, 1), 380);
  CU_ASSERT_EQUAL(get(result, 0, 2), 172);
  CU_ASSERT_EQUAL(get(result, 1, 0), 215);
  CU_ASSERT_EQUAL(get(result, 1, 1), 424);
  CU_ASSERT_EQUAL(get(result, 1, 2), 386);
  CU_ASSERT_EQUAL(get(result, 2, 0), 163);
  CU_ASSERT_EQUAL(get(result, 2, 1), 371);
  CU_ASSERT_EQUAL(get(result, 2, 2), 259);

  deallocate_matrix(result);
  deallocate_matrix(mat1);
  deallocate_matrix(mat2);

  // subtest 3
  result = NULL;
  mat1 = NULL;
  mat2 = NULL;

  CU_ASSERT_EQUAL(allocate_matrix(&result, 3, 3), 0);
  CU_ASSERT_EQUAL(allocate_matrix(&mat1, 3, 3), 0);
  CU_ASSERT_EQUAL(allocate_matrix(&mat2, 3, 3), 0);

  set(mat1, 0, 0, 12);
  set(mat1, 0, 1, 8);
  set(mat1, 0, 2, 4);

  set(mat1, 1, 0, 3);
  set(mat1, 1, 1, 17);
  set(mat1, 1, 2, 14);

  set(mat1, 2, 0, 9);
  set(mat1, 2, 1, 8);
  set(mat1, 2, 2, 10);


  set(mat2, 0, 0, -5);
  set(mat2, 0, 1, -19);
  set(mat2, 0, 2, -3);

  set(mat2, 1, 0, -6);
  set(mat2, 1, 1, -15);
  set(mat2, 1, 2, -9);

  set(mat2, 2, 0, -7);
  set(mat2, 2, 1, -8);
  set(mat2, 2, 2, -16);

  mul_matrix(result, mat1, mat2);
  CU_ASSERT_EQUAL(get(result, 0, 0), -136);
  CU_ASSERT_EQUAL(get(result, 0, 1), -380);
  CU_ASSERT_EQUAL(get(result, 0, 2), -172);
  CU_ASSERT_EQUAL(get(result, 1, 0), -215);
  CU_ASSERT_EQUAL(get(result, 1, 1), -424);
  CU_ASSERT_EQUAL(get(result, 1, 2), -386);
  CU_ASSERT_EQUAL(get(result, 2, 0), -163);
  CU_ASSERT_EQUAL(get(result, 2, 1), -371);
  CU_ASSERT_EQUAL(get(result, 2, 2), -259);

  deallocate_matrix(result);
  deallocate_matrix(mat1);
  deallocate_matrix(mat2);

}


/************* Test Runner Code goes here **************/

int main (void)
{
  Py_Initialize(); // Need to call this so that Python.h functions won't segfault
  CU_pSuite pSuite = NULL;

  /* initialize the CUnit test registry */
  if (CU_initialize_registry() != CUE_SUCCESS)
    return CU_get_error();

  /* add a suite to the registry */
  pSuite = CU_add_suite("mat_test_suite", init_suite, clean_suite);
  if (pSuite == NULL) {
    CU_cleanup_registry();
    return CU_get_error();
  }

   /* add the tests to the suite */
   if ((CU_add_test(pSuite, "add_test", add_test) == NULL) ||
        /* (OPTIONAL) Uncomment the following lines if you have implemented sub_matrix and neg_matrix. */
        (CU_add_test(pSuite, "sub_test", sub_test) == NULL) ||
        (CU_add_test(pSuite, "neg_test", neg_test) == NULL) ||
        (CU_add_test(pSuite, "mul_square_test", mul_square_test) == NULL) ||
        (CU_add_test(pSuite, "mul_non_square_test", mul_non_square_test) == NULL) ||
        (CU_add_test(pSuite, "abs_test", abs_test) == NULL) ||
        (CU_add_test(pSuite, "pow_test", pow_test) == NULL) ||
        (CU_add_test(pSuite, "alloc_fail_test", alloc_fail_test) == NULL) ||
        (CU_add_test(pSuite, "alloc_success_test", alloc_success_test) == NULL) ||
        (CU_add_test(pSuite, "alloc_ref_fail_test", alloc_ref_fail_test) == NULL) ||
        (CU_add_test(pSuite, "alloc_ref_success_test", alloc_ref_success_test) == NULL) ||
        (CU_add_test(pSuite, "dealloc_null_test", dealloc_null_test) == NULL) ||
        (CU_add_test(pSuite, "get_test", get_test) == NULL) ||
        (CU_add_test(pSuite, "set_test", set_test) == NULL) ||
        (CU_add_test(pSuite, "transpose_test", transpose_test) == NULL) ||
        (CU_add_test(pSuite, "mul_comp_test1", mul_comp_test1) == NULL)
     )
   {
      CU_cleanup_registry();
      return CU_get_error();
   }

  // Run all tests using the basic interface
  //CU_basic_set_mode(CU_BRM_NORMAL);
  CU_basic_set_mode(CU_BRM_VERBOSE);
  CU_basic_run_tests();
  printf("\n");
  CU_basic_show_failures(CU_get_failure_list());
  printf("\n\n");

  /* Clean up registry and return */
  CU_cleanup_registry();
  return CU_get_error();
}
