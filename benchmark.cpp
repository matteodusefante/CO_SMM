// CO_SMM

#include <tgmath.h>

#include "utils.h"
// DEPENDECIES cascading.h utils.h
#include "cascading.h"
#include "coalesced.h"
#include "sparsegen.h"

#include "algorithm.h"

#include <chrono>
#include <ctime>
#include <iostream>
#include <math.h>
#include <string>

#define probability 10

int main(__attribute__((unused)) int argc,
         __attribute__((unused)) char *argv[]) {

   // std::srand(std::time(nullptr));
   // std::srand(1);

   // typedef int32_t e_type; // matrix entry type is unsigned 32-bit integers
   // typedef utils::entry<e_type> ae_type;

   // size_t A_c = 5;
   // size_t C_c = 5;
   // size_t n = 100;
   // // size_t k;
   // size_t A_h = A_c * n;
   // size_t C_h = C_c * n;

   // utils::entry<utils::e_type> *A_matrix = new
   // utils::entry<utils::e_type>[A_h]; utils::entry<utils::e_type> *C_matrix =
   // new utils::entry<utils::e_type>[C_h];

   // utils::e_type *buff = new utils::e_type[n];

   // for (size_t index = 0; index < n; ++index)
   //    buff[index] = index;

   // std::srand(1);

   // for (size_t column = 0; column < n; ++column) {
   //    std::random_shuffle(buff, buff + n);
   //    for (size_t row = 0; row < n; ++row) {
   //       A_matrix[column * n + row] = utils::entry<utils::e_type>(
   //           rand() % RAND_MAX - RAND_MAX, buff[row], column);
   //    }
   // }
   //
   // for (size_t row = 0; row < n; ++row) {
   //    std::random_shuffle(buff, buff + n);
   //    for (size_t column = 0; column < n; ++column) {
   //       C_matrix[row * n + column] = utils::entry<utils::e_type>(
   //           rand() % RAND_MAX - RAND_MAX, row, buff[column]);
   //    }
   // }

   // typedef cascading::augmented_entry<utils::entry<utils::e_type>> ae_type;
   // ae_type *index_vector = new ae_type[100];
   //
   // for (size_t index = 0; index < 100; ++index)
   //    index_vector[index] = cascading::augmented_entry<utils::e_type>();

   // cascading::augmented_entry<ae_type> *A_index_vector =
   //     new cascading::augmented_entry<ae_type>[n];
   //
   // cascading::augment_matrix(A_matrix, A_index_vector, A_h, n);

   std::srand(std::time(nullptr));

   size_t d = pow(2, 2);
   size_t n = pow(2, 4);
   size_t A_nnz = 0;
   size_t C_nnz = 0;
   size_t AC_nnz = 0;
   size_t k = pow(d, 2) * n / log2(n);

   utils::e_type **A_matrix = new utils::e_type *[n];
   utils::e_type **C_matrix = new utils::e_type *[n];
   utils::e_type **AC_matrix = new utils::e_type *[n];

   for (size_t i = 0; i < n; ++i) {
      A_matrix[i] = new utils::e_type[n]();
      C_matrix[i] = new utils::e_type[n]();
      AC_matrix[i] = new utils::e_type[n]();
   }

   // for (size_t i = 0; i < n; ++i) {
   //    for (size_t j = 0; j < n; ++j) {
   //       A_matrix[i][j] = 1;
   //       C_matrix[i][j] = 1;
   //    }
   // }
   sparsegen::generate_matrix(A_matrix, C_matrix, n, d, k);

   utils::matrix_multiplication(AC_matrix, A_matrix, C_matrix, n);

   for (size_t i = 0; i < n; ++i)
      for (size_t j = 0; j < n; ++j) {
         if (A_matrix[i][j] != 0)
            A_nnz++;
         if (C_matrix[i][j] != 0)
            C_nnz++;
      }

   auto *A_sparse = utils::sparsify_matrix(A_matrix, A_nnz, n);
   auto *C_sparse = utils::sparsify_matrix(C_matrix, C_nnz, n);

   utils::sparse_transpose(A_sparse, A_nnz, n, utils::COLUMN_MAJOR);
   utils::sparse_transpose(C_sparse, C_nnz, n, utils::ROW_MAJOR);

   utils::print_sparse_matrix(A_sparse, A_nnz);
   utils::print_sparse_matrix(C_sparse, A_nnz);

   utils::sparse_prefix_sum(A_sparse, A_nnz, utils::COLUMN);
   utils::sparse_prefix_sum(C_sparse, C_nnz, utils::ROW);

   utils::print_sparse_matrix(A_sparse, A_nnz);
   utils::print_sparse_matrix(C_sparse, A_nnz);

   utils::print_sparse_as_dense(A_sparse, A_nnz, n);
   utils::print_sparse_as_dense(C_sparse, C_nnz, n);

   // auto A_index = coalesced::range_coalesced(A_sparse, A_nnz, n, utils::ROW);
   // auto C_index = coalesced::range_coalesced(C_sparse, C_nnz, n,
   // utils::COLUMN);
   auto A_index = cascading::augment_matrix(A_sparse, A_nnz, n, utils::ROW);
   auto C_index = cascading::augment_matrix(C_sparse, C_nnz, n, utils::COLUMN);

   algorithm::matrix_product(A_index, C_index, n);

   printf("A_nnz+C_nnz=%ld+%ld=%ld, AC_nnz=%ld, h/log n=%lf\n", A_nnz, C_nnz,
          A_nnz + C_nnz, AC_nnz, ((A_nnz + C_nnz) / log2(d)));

#ifdef DEBUG
   // utils::print_matrix(matrix, d);
   // utils::print_matrix(A_matrix, n);
   // utils::print_matrix(C_matrix, n);

#endif

   return 0;

   /***********/

   // for (size_t index1 = 0; index1 < n; ++index1) {
   //    for (size_t index2 = 0; index2 < n; ++index2)
   //       printf("%d\t", A_matrix[index1][index2]);
   //    printf("\n");
   // }
   // printf("\n");
   //
   // for (size_t index1 = 0; index1 < n; ++index1) {
   //    for (size_t index2 = 0; index2 < n; ++index2)
   //       printf("%d\t", C_matrix[index1][index2]);
   //    printf("\n");
   // }
   // printf("\n");

   // for (size_t index1 = 0; index1 < n; ++index1) {
   //    for (size_t index2 = 0; index2 < n; ++index2)
   //       printf("%d\t", AC_matrix[index1][index2]);
   //    printf("\n");
   // }
   // printf("\n");

   //    uint32_t length = 100;
   //    uint32_t **out = new uint32_t *[length];
   //    uint32_t **test = new uint32_t *[length];
   //    uint32_t **matrix = new uint32_t *[length];
   //    uint32_t **count_out = new uint32_t *[length];
   //
   //    std::chrono::duration<double_t> ttime;
   //    std::chrono::high_resolution_clock::time_point start_time,
   //    stop_time;
   //
   //    for (size_t i = 0; i < length; ++i) {
   //       out[i] = new uint32_t[length]();
   //       test[i] = new uint32_t[length]();
   //       matrix[i] = new uint32_t[length]();
   //       count_out[i] = new uint32_t[length]();
   //    }
   //
   //    // std::srand(std::time(nullptr));
   //
   //    std::srand(1);
   //
   //    for (size_t i = 0; i < length; ++i) {
   //       for (size_t j = 0; j < length; ++j) {
   //          if (i == j)
   //             matrix[i][j] = 1;
   //          else
   //             matrix[i][j] = 1; // std::rand() % 10 == 1 ? 1 : 0;
   //       }
   //    }
   //
   //    start_time = std::chrono::high_resolution_clock::now();
   //
   //    for (size_t i = 0; i < length; ++i)
   //       for (size_t j = 0; j < length; ++j) {
   //          for (size_t k = 0; k < length; ++k) {
   //             if (matrix[i][k] != 0 && matrix[k][j] != 0) {
   //                test[i][j] = 1;
   //                break;
   //             }
   //          }
   //          // test[i][j] += matrix[i][k] * matrix[k][j];
   //          // test[i][j] = test[i][j] > 0 ? 1 : 0;
   //       }
   //
   //    stop_time = std::chrono::high_resolution_clock::now();
   //    ttime = std::chrono::duration_cast<std::chrono::duration<double_t>
   //    /**/>(
   //        stop_time - start_time);
   //
   //    std::cout << ttime.std::chrono::duration<double_t>::count() <<
   //    std::endl;
   //
   //    std::cout << "Completed" << std::endl;
   //
   //    start_time = std::chrono::high_resolution_clock::now();
   //
   //    for (size_t i = 0; i < length; ++i)
   //       for (size_t j = 0; j < length; ++j) {
   //          if (i == j) {
   //             out[i][j] = 1;
   //             // count_out[i][j]++;
   //          } else {
   //             if (matrix[i][j] != 0) {
   //                for (size_t k = 0; k < length; ++k) {
   //                   // count_out[i][k]++;
   //                   // if (out[j][k] == 0)
   //                   if ((matrix[j][k] != 0))
   //                      out[i][k] = 1;
   //                }
   //             }
   //          }
   //       }
   //
   //    stop_time = std::chrono::high_resolution_clock::now();
   //    ttime = std::chrono::duration_cast<std::chrono::duration<double_t>
   //    /**/>(
   //        stop_time - start_time);
   //
   //    std::cout << ttime.std::chrono::duration<double_t>::count() <<
   //    std::endl;
   //
   //    std::cout << "Completed" << std::endl;
   //
   // #ifdef DEBUG
   //    std::cout << "Matrix" << std::endl;
   //    for (size_t i = 0; i < length; ++i) {
   //       for (size_t j = 0; j < length; ++j) {
   //          std::cout << matrix[i][j] << "  ";
   //       }
   //       std::cout << std::endl;
   //    }
   //
   //    std::cout << std::endl << "Out" << std::endl;
   //    for (size_t i = 0; i < length; ++i) {
   //       for (size_t j = 0; j < length; ++j) {
   //          std::cout << out[i][j] << "  ";
   //       }
   //       std::cout << std::endl;
   //    }
   //
   //    std::cout << std::endl << "Test" << std::endl;
   //
   //    for (size_t i = 0; i < length; ++i) {
   //       for (size_t j = 0; j < length; ++j) {
   //          std::cout << test[i][j] << "  ";
   //       }
   //       std::cout << std::endl;
   //    }
   //    std::cout << std::endl;
   //
   //    std::cout << std::endl << "Count" << std::endl;
   //    for (size_t i = 0; i < length; ++i) {
   //       for (size_t j = 0; j < length; ++j) {
   //          std::cout << count_out[i][j] << "  ";
   //       }
   //       std::cout << std::endl;
   //    }
   // #endif
   //
   //    for (size_t i = 0; i < length; ++i)
   //       for (size_t j = 0; j < length; ++j)
   //          if (test[i][j] != out[i][j]) {
   //             std::cout << "Error" << std::endl;
   //             return 0;
   //          }
   //
   //    std::cout << "Exit" << std::endl;
   //
   //    for (size_t i = 0; i < length; ++i) {
   //       delete[] out[i];
   //       delete[] test[i];
   //       delete[] matrix[i];
   //    }
   //
   //    delete[] out;
   //    delete[] test;
   //    delete[] matrix;
   //
   //    return 0;
}
