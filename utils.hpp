
#include <algorithm>
#include <iostream>
#include <stdlib.h>
#include <string>
#include <utility>
#include <vector>

#include <cassert>
// uncomment to disable assert()
// #define NDEBUG

#define DEBUG
#define DEBUG_VERBOSE

// #ifdef DEBUG
// #define DEBUG_STDERR(x) (std::cerr << (x))
// #define DEBUG_STDOUT(x) (std::cout << (x))
// #else
// #define DEBUG_STDERR(x)
// #define DEBUG_STDOUT(x)
// #endif

namespace utils {

static const size_t COLUMN_MAJOR = 0;
static const size_t ROW_MAJOR = 1;
static const size_t COLUMN = COLUMN_MAJOR;
static const size_t ROW = ROW_MAJOR;

static const size_t INHERITED = 0;
static const size_t NATIVE = 1;

template <typename T> inline T sparse_transpose(const T &matrix, matrix::index_t n);

// template <typename First, typename... Rest> void debug(First &&first, Rest &&... rest)
// {
//    std::cout << std::forward<First>(first) << " ";
//    debug(std::forward<Rest>(rest)...);
// }

/****************************************************************************************/

template <typename T> struct return_index {

   matrix::index_t operator()(const T &i, const T &j, size_t layout) const {
      if (layout == COLUMN_MAJOR)
         return j;
      else if (layout == ROW_MAJOR)
         return i;
      else
         exit(EXIT_FAILURE);
   }

   matrix::index_t operator()(const T &entry, size_t layout) const {
      if (layout == COLUMN_MAJOR)
         return entry.j;
      else if (layout == ROW_MAJOR)
         return entry.i;
      else
         exit(EXIT_FAILURE);
   }
};

/****************************************************************************************/

// multiply a mxn matrix with a nxp matrix
template <typename T>
inline T matrix_multiplication(const T &A_matrix, const T &C_matrix, matrix::index_t m,
                               matrix::index_t n, matrix::index_t p) {

   matrix::dense<matrix::entry_t> AC_matrix(n);

   for (matrix::index_t i = 0; i < m; ++i)
      for (matrix::index_t j = 0; j < p; ++j)
         for (matrix::index_t k = 0; k < n; ++k)
            AC_matrix[i][j] += A_matrix[i][k] * C_matrix[k][j];

   return AC_matrix;
}

template <typename T>
inline T matrix_multiplication(const T &A_matrix, const T &C_matrix, matrix::index_t n) {

   return matrix_multiplication(A_matrix, C_matrix, n, n, n);
}
/****************************************************************************************/

template <typename T> inline void print_matrix(const T &matrix) {

   for (const auto &vec : matrix) {
      for (const auto &entry : vec)
         std::cout << entry << "\t";
      std::cout << std::endl;
   }
}

/****************************************************************************************/

template <typename T> inline void print_sparse_matrix(const T &matrix) {

   matrix.layout ? std::cout << "ROW " : std::cout << "COLUMN ";
   std::cout << "MAJOR" << std::endl;

   for (const auto &vec : matrix) {
      for (const auto &entry : vec)
         std::cout << "(" << entry.i << "," << entry.j << "," << entry.a << ") ";
      std::cout << std::endl;
   }
}

/****************************************************************************************/

template <typename T>
inline void print_sparse_as_dense(const T &matrix, matrix::index_t n) {

   return_index<matrix::entry<matrix::entry_t>> matrix_index;
   T temp;
   if (matrix.layout == utils::COLUMN_MAJOR)
      temp = sparse_transpose(matrix, n); // transpose only if in CML
   else
      temp = matrix;

   for (const auto &vec : temp) {
      auto it = vec.begin();
      for (matrix::index_t index = 0; index < n; ++index)
         if ((it != vec.end()) && (matrix_index(*it, 1 - temp.layout) == index))
            std::cout << it++->a << "\t";
         else
            std::cout << "0\t";
      std::cout << std::endl;
   }
}

/****************************************************************************************/

template <typename T> inline void dense_transpose(T &matrix, matrix::index_t n) {

   for (matrix::index_t i = 0; i < n; ++i)
      for (matrix::index_t j = 0; j < n; ++j)
         std::swap(matrix[i][j], matrix[j][i]);
}

/****************************************************************************************/

template <typename T> inline void sparse_prefix_sum(T &matrix) {

   for (auto &vector : matrix) {
      matrix::entry_t sum = 0;
      for (auto &entry : vector)
         entry.a = sum += entry.a;
   }
}

/****************************************************************************************/

template <typename T> inline T sparse_transpose(const T &matrix, matrix::index_t n) {

   T temp_matrix(n);
   return_index<matrix::entry<matrix::entry_t>> matrix_index;
   temp_matrix.layout = 1 - matrix.layout;

   for (auto &vec : matrix)
      for (auto &entry : vec)
         temp_matrix[matrix_index(entry, temp_matrix.layout)].push_back(entry);
   return temp_matrix;
}

/****************************************************************************************/

template <typename T, typename S>
inline void vector_matrix_multiplication(const T &matrix, const T &out_matrix,
                                         const S &vector, matrix::index_t h) {

   for (matrix::index_t index = 0; index < h; ++index)
      out_matrix[index]->a = matrix[index] * vector[matrix[index]->i];
}

/****************************************************************************************/

template <typename T, typename S>
inline void matrix_vector_multiplication(const T &matrix, const S &out_vector,
                                         const S &vector, matrix::index_t h) {

   for (matrix::index_t index = 0; index < h; ++index)
      out_vector[matrix[index]->i] += matrix[index] * vector[matrix[index]->i];
}

/****************************************************************************************/

template <typename T>
inline void swap_rows(T &matrix, matrix::index_t i, matrix::index_t j,
                      matrix::index_t n) {

   // matrix[i].vector::swap(matrix[j]);

   for (matrix::index_t index = 0; index < n; ++index)
      std::swap(matrix[i][index], matrix[j][index]);
}

/****************************************************************************************/

template <typename T>
inline void swap_columns(T &matrix, matrix::index_t i, matrix::index_t j,
                         matrix::index_t n) {

   for (matrix::index_t index = 0; index < n; ++index)
      std::swap(matrix[index][i], matrix[index][j]);
}

/****************************************************************************************/

template <typename T> inline void randomize_matrix(T &matrix, matrix::index_t n) {

   for (matrix::index_t i = 0; i < n; ++i) {
      matrix::entry_t r = std::rand() % 100;
      for (matrix::index_t j = 0; j < n; ++j) {
         matrix[i][j] *= r;
      }
   }

   std::vector<unsigned int> perm(n);
   std::iota(std::begin(perm), std::end(perm), 1);

   std::random_shuffle(std::begin(perm), std::end(perm));

   for (matrix::index_t i = 0; i < n; ++i) {
      swap_rows(matrix, i, perm[i], n);
      swap_columns(matrix, i, perm[i], n);
   }
}

/****************************************************************************************/

template <typename T>
inline matrix::sparse<matrix::entry<matrix::entry_t>>
sparsify_matrix(const T &matrix, matrix::index_t n, ssize_t layout) {

   matrix::sparse<matrix::entry<matrix::entry_t>> sparse_matrix(n, layout);

   for (matrix::index_t i = 0; i < n; ++i)
      for (matrix::index_t j = 0; j < n; ++j) {
         if (matrix[i][j] != 0) {
            if (layout == ROW_MAJOR) {
               sparse_matrix[i].emplace_back(
                   matrix::entry<matrix::entry_t>(i, j, matrix[i][j]));
            } else {
               sparse_matrix[j].emplace_back(
                   matrix::entry<matrix::entry_t>(i, j, matrix[i][j]));
            }
         }
      }
   return sparse_matrix;
}

} // namespace utils
