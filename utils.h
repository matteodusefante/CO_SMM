
#include <algorithm>
#include <stdlib.h>
#include <string>

// uncomment to disable assert()
// #define NDEBUG
#include <cassert>

#define DEBUG
#define DEBUG_VERBOSE

namespace utils {

typedef int32_t e_type;

static const ssize_t COLUMN_MAJOR = 0;
static const ssize_t ROW_MAJOR = 1;
static const ssize_t COLUMN = COLUMN_MAJOR;
static const ssize_t ROW = ROW_MAJOR;

static const size_t INHERITED = 0;
static const size_t NATIVE = 1;

/****************************************************************************************/

template <typename T> class entry {

 public:
   entry(){};
   ~entry(){};
   entry(size_t i, size_t j, const T &a);

   entry<T> &operator=(const entry<T> &rhs) {
      this->i = rhs.i;
      this->j = rhs.j;
      this->a = rhs.a;
      return *this; // return the result by reference
   }

   size_t i, j;
   T a;

   void set_i(size_t i) { this->i = i; }
   void set_j(size_t j) { this->j = j; }
   void set_a(const T &a) { this->a = a; }
};

template <typename T>
entry<T>::entry(size_t i, size_t j, const T &a) : i(i), j(j), a(a) {
   ;
}

typedef entry<e_type> ee;

/****************************************************************************************/

template <typename T> struct sort_by {
   size_t operator()(const T &entry, ssize_t type) const {
      if (type == COLUMN_MAJOR)
         return entry.j;
      else if (type == ROW_MAJOR)
         return entry.i;
      else
         exit(EXIT_FAILURE);
   }
};
/****************************************************************************************/

template <typename T>
inline void matrix_multiplication(const T &AC_matrix, const T &A_matrix,
                                  const T &C_matrix, size_t n) {

   for (size_t i = 0; i < n; ++i)
      for (size_t j = 0; j < n; ++j)
         for (size_t k = 0; k < n; ++k)
            AC_matrix[i][j] += A_matrix[i][k] * C_matrix[k][j];
}

/****************************************************************************************/

template <typename T> inline void print_matrix(const T &matrix, size_t n) {

   for (size_t i = 0; i < n; ++i) {
      for (size_t j = 0; j < n; ++j)
         printf("%d\t", matrix[i][j]);
      printf("\n");
   }
   printf("\n");
}

/****************************************************************************************/

template <typename T>
inline void print_sparse_matrix(const T &matrix, size_t h) {

   for (size_t index = 0; index < h; ++index)
      printf("(%ld,%ld,%d) ", matrix[index].i, matrix[index].j,
             matrix[index].a);
   printf("\n\n");
}

/****************************************************************************************/

template <typename T>
inline void print_sparse_as_dense(const T &matrix, size_t h, size_t n) {

   bool found = false;
   for (size_t i = 0; i < n; ++i) {
      for (size_t j = 0; j < n; ++j) {
         for (size_t k = 0; k < h; ++k) {
            if (matrix[k].entry == nullptr)
               continue;
            if ((matrix[k].entry->i == i) && (matrix[k].entry->j == j) &&
                matrix[k].entry_type == INHERITED) {
               found = true;
               printf("%dH\t", matrix[k].entry->a);
            }
            if ((matrix[k].entry->i == i) && (matrix[k].entry->j == j) &&
                matrix[k].entry_type == NATIVE) {
               found = true;
               printf("%d\t", matrix[k].entry->a);
            }
         }
         if (!found)
            printf("0\t");
         found = false;
      }
      printf("\n");
   }
   printf("\n");
}

/****************************************************************************************/
// Template specialization for utils::ee *

inline void print_sparse_as_dense(utils::ee *const &matrix, size_t h,
                                  size_t n) {

   bool found = false;
   for (size_t i = 0; i < n; ++i) {
      for (size_t j = 0; j < n; ++j) {
         for (size_t k = 0; k < h; ++k) {
            if ((matrix[k].i == i) && (matrix[k].j == j)) {
               found = true;
               printf("%d\t", matrix[k].a);
            }
         }
         if (!found)
            printf("0\t");
         found = false;
      }
      printf("\n");
   }
   printf("\n");
}

/****************************************************************************************/

template <typename T> inline void dense_transpose(const T &matrix, size_t n) {

   for (size_t i = 0; i < n; ++i)
      for (size_t j = 0; j < n; ++j)
         std::swap(matrix[i][j], matrix[i][j]);
}

/****************************************************************************************/

template <typename T>
inline void dense_column_prefix_sum(const T &matrix, size_t n) {

   for (size_t j = 0; j < n; ++j)
      for (size_t i = 1; i < n; ++i)
         matrix[i][j] += matrix[i - 1][j];
}

/****************************************************************************************/

template <typename T>
inline void sparse_prefix_sum(const T &matrix, size_t h, size_t order) {

   // sort_by returns either row or column index depending on order
   utils::sort_by<utils::ee> index_type;

   size_t index = 0;
   size_t sum = 0;

   for (auto *entry = &matrix[0]; entry != &matrix[h]; ++entry) {
      if (index_type(*entry, order) == index){
         sum = entry->a += sum;
      } else {
         index = index_type(*entry, order);
         sum = entry->a ;
      }
   }
}

/****************************************************************************************/

template <typename T>
inline void sparse_transpose(const T &matrix, size_t h, size_t n, size_t type) {

   size_t *counters = new size_t[n + 1]();
   utils::sort_by<utils::entry<utils::e_type>> sb;

   for (size_t index = 0; index < h; ++index)
      counters[sb(matrix[index], type) + 1]++;

   for (size_t index = 1; index < n + 1; ++index)
      counters[index] += counters[index - 1];

   utils::entry<utils::e_type> *transpose = new utils::entry<utils::e_type>[h];

   for (size_t index = 0; index < h; ++index)
      transpose[counters[sb(matrix[index], type)]++] = matrix[index];

   std::copy(transpose, transpose + h, matrix);
   delete[] transpose;
   delete[] counters;
}

/****************************************************************************************/

template <typename T, typename S>
inline void vector_matrix_multiplication(const T &matrix, const T &out_matrix,
                                         const S &vector, size_t h) {

   for (size_t index = 0; index < h; ++index)
      out_matrix[index]->a = matrix[index] * vector[matrix[index]->i];
}

/****************************************************************************************/

template <typename T, typename S>
inline void matrix_vector_multiplication(const T &matrix, const S &out_vector,
                                         const S &vector, size_t h) {

   for (size_t index = 0; index < h; ++index)
      out_vector[matrix[index]->i] += matrix[index] * vector[matrix[index]->i];
}

/****************************************************************************************/

template <typename T>
inline size_t inner_product(const T &vec1, const T &vec2, size_t n) {

   size_t temp = 0;

   for (size_t index = 0; index < n; ++index)
      temp += vec1[index] * vec2[index];

   printf("= %ld \n", temp);

   return temp;
}

/****************************************************************************************/

template <typename T>
inline void swap_rows(const T &matrix, size_t i, size_t j, size_t n) {

   for (size_t index = 0; index < n; ++index)
      std::swap(matrix[i][index], matrix[j][index]);
}

/****************************************************************************************/

template <typename T>
inline void swap_columns(const T &matrix, size_t i, size_t j, size_t n) {

   for (size_t index = 0; index < n; ++index)
      std::swap(matrix[index][i], matrix[index][j]);
}

/****************************************************************************************/

template <typename T>
inline void randomize_matrix(const T &A_matrix, const T &C_matrix, size_t n) {

   size_t r, s;
   for (size_t i = 0; i < n; ++i) {
      r = std::rand() % 100;
      s = std::rand() % 100;
      for (size_t j = 0; j < n; ++j) {
         A_matrix[i][j] *= r;
         C_matrix[j][i] *= s;
      }
   }

   T *perm = new T[n];

   for (size_t i = 0; i < n; ++i)
      perm[i] = i;

   for (size_t i = 0; i < n; ++i) {
      size_t index = (std::rand() % (n - i)) + i;
      std::swap(perm[index], perm[i]);
   }

   for (size_t i = 0; i < n; ++i) {
      swap_rows(A_matrix, i, perm[i], n);
      swap_columns(C_matrix, i, perm[i], n);
   }

   for (size_t i = 0; i < n; ++i) {
      size_t index = (std::rand() % (n - i)) + i;
      std::swap(perm[index], perm[i]);
   }

   delete[] perm;
}

/****************************************************************************************/

template <typename T>
inline utils::entry<T> *sparsify_matrix(T **matrix, size_t n) {

   utils::entry<T> *temp = new utils::entry<T>[(size_t)pow(n, 2)];
   size_t h = 0;
   for (size_t i = 0; i < n; ++i)
      for (size_t j = 0; j < n; ++j)
         if (matrix[i][j] != 0)
            temp[h++] = utils::entry<T>(i, j, matrix[i][j]);

   utils::entry<T> *sparse_matrix = new utils::entry<T>[h];
   std::copy(temp, temp + h, sparse_matrix);
   delete[] temp;
   return sparse_matrix;
}

/****************************************************************************************/

template <typename T>
inline utils::entry<T> *sparsify_matrix(T **matrix, size_t h, size_t n) {

   utils::entry<T> *sparse_matrix = new utils::entry<T>[h];

   size_t index = 0;
   for (size_t i = 0; i < n; ++i)
      for (size_t j = 0; j < n; ++j)
         if (matrix[i][j] != 0)
            sparse_matrix[index++] = utils::entry<T>(i, j, matrix[i][j]);

   return sparse_matrix;
}

} // namespace utils
