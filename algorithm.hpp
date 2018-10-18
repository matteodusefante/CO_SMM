
namespace algorithm {

/****************************************************************************************/

// inline matrix::entry_t &
// extract_sketch(coalesced::bucket<utils::ee> **const &index_vector, size_t i1, size_t
// i2,
//                size_t n, size_t layout) {
//
//    return coalesced::extract_sketch(index_vector, i1, i2, n, layout);
// }

inline std::vector<matrix::entry_t>
extract_sketch(const coalesced::coalesced_buckets<matrix::coo_entry> &ds,
               matrix::index_t i1, matrix::index_t i2, matrix::index_t n,
               ssize_t layout) {

   return coalesced::extract_sketch(ds, i1, i2, n, layout);
}

/****************************************************************************************/

inline std::vector<matrix::entry_t>
extract_sketch(cascading::augmented_matrix<matrix::coo_entry> &ds, matrix::index_t i1,
               matrix::index_t i2, matrix::index_t n, ssize_t layout) {

   return cascading::extract_sketch(ds, i1, i2, n, layout);
}

/****************************************************************************************/

template <typename T>
void matrix_product(T &A_ds, T &C_ds, matrix::index_t i1, matrix::index_t i2,
                    matrix::index_t j1, matrix::index_t j2, matrix::index_t n) {

   std::cout << "[" << i1 << "," << i2 << "]x[" << j1 << "," << j2 << "]" << std::endl;

   matrix::index_t im, jm;

   if (i1 == i2 && j1 == j2) {
      std::cout << "[" << i1 << "," << i2 << "]x[" << j1 << "," << j2 << "]" << std::endl;
      matrix::inner_product(extract_sketch(A_ds, i1, i2, n, utils::ROW),
                           extract_sketch(C_ds, j1, j2, n, utils::COLUMN));
      return;
   } else if (i1 == i2 && j1 < j2) {
      jm = j1 + (j2 - j1) / 2;
      std::cout << "[" << i1 << "," << i2 << "]x[" << j1 << "," << jm << "]" << std::endl;
      if (matrix::inner_product(extract_sketch(A_ds, i1, i2, n, utils::ROW),
                               extract_sketch(C_ds, j1, jm, n, utils::COLUMN)) != 0)
         matrix_product(A_ds, C_ds, i1, i2, j1, jm, n);
      std::cout << "[" << i1 << "," << i2 << "]x[" << jm + 1 << "," << j2 << "]"
                << std::endl;
      if (matrix::inner_product(extract_sketch(A_ds, i1, i2, n, utils::ROW),
                               extract_sketch(C_ds, jm + 1, j2, n, utils::COLUMN)) != 0)
         matrix_product(A_ds, C_ds, i1, i2, jm + 1, j2, n);
   } else if (i1 < i2 && j1 == j2) {
      im = i1 + (i2 - i1) / 2;
      std::cout << "[" << i1 << "," << im << "]x[" << j1 << "," << j2 << "]" << std::endl;
      if (matrix::inner_product(extract_sketch(A_ds, i1, im, n, utils::ROW),
                               extract_sketch(C_ds, j1, j2, n, utils::COLUMN)) != 0)
         matrix_product(A_ds, C_ds, i1, im, j1, j2, n);
      std::cout << "[" << i1 << "," << im + 1 << "]x[" << j1 << "," << j2 << "]"
                << std::endl;
      if (matrix::inner_product(extract_sketch(A_ds, im + 1, i2, n, utils::ROW),
                               extract_sketch(C_ds, j1, j2, n, utils::COLUMN)) != 0)
         matrix_product(A_ds, C_ds, im + 1, i2, j1, j2, n);
   } else {
      im = i1 + (i2 - i1) / 2;
      jm = j1 + (j2 - j1) / 2;
      // DEBUG_STDOUT("[" + i1 + "," + i2 + "]x[" + j1 + "," + j2 + "]\n");
      std::cout << "[" << i1 << "," << im << "]x[" << j1 << "," << jm << "]" << std::endl;
      if (matrix::inner_product(extract_sketch(A_ds, i1, im, n, utils::ROW),
                               extract_sketch(C_ds, j1, jm, n, utils::COLUMN)) != 0)
         matrix_product(A_ds, C_ds, i1, im, j1, jm, n);
      std::cout << "[" << i1 << "," << im << "]x[" << jm + 1 << "," << j2 << "]"
                << std::endl;
      if (matrix::inner_product(extract_sketch(A_ds, i1, im, n, utils::ROW),
                               extract_sketch(C_ds, jm + 1, j2, n, utils::COLUMN)) != 0)
         matrix_product(A_ds, C_ds, i1, im, jm + 1, j2, n);
      std::cout << "[" << i1 << "," << im + 1 << "]x[" << j1 << "," << jm << "]"
                << std::endl;
      if (matrix::inner_product(extract_sketch(A_ds, im + 1, i2, n, utils::ROW),
                               extract_sketch(C_ds, j1, jm, n, utils::COLUMN)) != 0)
         matrix_product(A_ds, C_ds, im + 1, i2, j1, jm, n);
      std::cout << "[" << i1 << "," << im + 1 << "]x[" << j1 << "," << jm + 1 << "]"
                << std::endl;
      if (matrix::inner_product(extract_sketch(A_ds, im + 1, i2, n, utils::ROW),
                               extract_sketch(C_ds, jm + 1, j2, n, utils::COLUMN)) != 0)
         matrix_product(A_ds, C_ds, im + 1, i2, jm + 1, j2, n);
   }
}

/****************************************************************************************/

template <typename T> void matrix_product(T &A_ds, T &C_ds, size_t n) {

   matrix_product(A_ds, C_ds, 0, n - 1, 0, n - 1, n);
}

} // namespace algorithm