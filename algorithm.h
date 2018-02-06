
namespace algorithm {

/****************************************************************************************/

inline utils::e_type *
extract_sketch(coalesced::bucket<utils::ee> **const &index_vector, size_t i1,
               size_t i2, size_t n, size_t order) {

   return coalesced::extract_sketch(index_vector, i1, i2, n, order);
}

/****************************************************************************************/

inline utils::e_type *
extract_sketch(cascading::augmented_entry<utils::ee> *const &index_vector,
               size_t i1, size_t i2, size_t n, size_t order) {

   return cascading::extract_sketch(index_vector, i1, i2, n, order);
}

/****************************************************************************************/

template <class T>
inline void matrix_product(const T &A_index, const T &C_index, size_t i1,
                           size_t i2, size_t j1, size_t j2, size_t n) {

   printf("[%ld,%ld]x[%ld,%ld] \n", i1, i2, j1, j2);

   size_t im, jm;

   if (i1 == i2 && j1 == j2) {
      printf("[%ld,%ld]x[%ld,%ld] \n", i1, i2, j1, j2);
      utils::inner_product(
          algorithm::extract_sketch(A_index, i1, i2, n, utils::ROW),
          algorithm::extract_sketch(C_index, j1, j2, n, utils::COLUMN), n);
      return;
   } else if (i1 == i2 && j1 < j2) {
      jm = j1 + (j2 - j1) / 2;
      if (utils::inner_product(
              algorithm::extract_sketch(A_index, i1, i2, n, utils::ROW),
              algorithm::extract_sketch(C_index, j1, jm, n, utils::COLUMN),
              n) != 0)
         algorithm::matrix_product(A_index, C_index, i1, i2, j1, jm, n);
      printf("[%ld,%ld]x[%ld,%ld] \n", i1, i2, jm + 1, j2);
      if (utils::inner_product(
              algorithm::extract_sketch(A_index, i1, i2, n, utils::ROW),
              algorithm::extract_sketch(C_index, jm + 1, j2, n, utils::COLUMN),
              n) != 0)
         algorithm::matrix_product(A_index, C_index, i1, i2, jm + 1, j2, n);
   } else if (i1 < i2 && j1 == j2) {
      im = i1 + (i2 - i1) / 2;
      printf("[%ld,%ld]x[%ld,%ld] \n", i1, im, j1, j2);
      if (utils::inner_product(
              algorithm::extract_sketch(A_index, i1, im, n, utils::ROW),
              algorithm::extract_sketch(C_index, j1, j2, n, utils::COLUMN),
              n) != 0)
         algorithm::matrix_product(A_index, C_index, i1, im, j1, j2, n);
      printf("[%ld,%ld]x[%ld,%ld] \n", im + 1, i2, j1, j2);
      if (utils::inner_product(
              algorithm::extract_sketch(A_index, im + 1, i2, n, utils::ROW),
              algorithm::extract_sketch(C_index, j1, j2, n, utils::COLUMN),
              n) != 0)
         algorithm::matrix_product(A_index, C_index, im + 1, i2, j1, j2, n);
   } else {
      im = i1 + (i2 - i1) / 2;
      jm = j1 + (j2 - j1) / 2;
      printf("[%ld,%ld]x[%ld,%ld] \n", i1, im, j1, jm);
      if (utils::inner_product(
              algorithm::extract_sketch(A_index, i1, im, n, utils::ROW),
              algorithm::extract_sketch(C_index, j1, jm, n, utils::COLUMN),
              n) != 0)
         algorithm::matrix_product(A_index, C_index, i1, im, j1, jm, n);
      printf("[%ld,%ld]x[%ld,%ld] \n", i1, im, jm + 1, j2);
      if (utils::inner_product(
              algorithm::extract_sketch(A_index, i1, im, n, utils::ROW),
              algorithm::extract_sketch(C_index, jm + 1, j2, n, utils::COLUMN),
              n) != 0)
         algorithm::matrix_product(A_index, C_index, i1, im, jm + 1, j2, n);
      printf("[%ld,%ld]x[%ld,%ld] \n", im + 1, i2, j1, jm);
      if (utils::inner_product(
              algorithm::extract_sketch(A_index, im + 1, i2, n, utils::ROW),
              algorithm::extract_sketch(C_index, j1, jm, n, utils::COLUMN),
              n) != 0)
         algorithm::matrix_product(A_index, C_index, im + 1, i2, j1, jm, n);
      printf("[%ld,%ld]x[%ld,%ld] \n", im + 1, i2, jm + 1, j2);
      if (utils::inner_product(
              algorithm::extract_sketch(A_index, im + 1, i2, n, utils::ROW),
              algorithm::extract_sketch(C_index, jm + 1, j2, n, utils::COLUMN),
              n) != 0)
         algorithm::matrix_product(A_index, C_index, im + 1, i2, jm + 1, j2, n);
   }
}

/****************************************************************************************/

template <typename T>
inline void matrix_product(const T &A_index, const T &C_index, size_t n) {

   algorithm::matrix_product(A_index, C_index, 0, n - 1, 0, n - 1, n);
}

} // namespace algorithm