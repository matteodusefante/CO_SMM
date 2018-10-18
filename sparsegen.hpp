

namespace sparsegen {

template <typename T>
inline T kronecker_product(const T &A_matrix, const T &B_matrix, matrix::index_t A_m,
                           matrix::index_t A_n, matrix::index_t B_m,
                           matrix::index_t B_n) {

   T AB_matrix(A_m * B_m, A_n * B_n); // A_m * B_m x A_n * B_n matrix

   typename T::iterator iter_type;

   for (matrix::index_t A_i = 0; A_i < A_m; ++A_i)
      for (matrix::index_t A_j = 0; A_j < A_n; ++A_j)
         for (matrix::index_t B_i = 0; B_i < B_m; ++B_i)
            for (matrix::index_t B_j = 0; B_j < B_n; ++B_j)
               AB_matrix[A_i * B_m + B_i][A_j * B_n + B_j] =
                   A_matrix[A_i][A_j] * B_matrix[B_i][B_j];
   return AB_matrix;
}

/****************************************************************************************/

matrix::dense<matrix::entry_t> identity(matrix::index_t n) {

   matrix::dense<matrix::entry_t> I_matrix(n);

   for (matrix::index_t i = 0; i < n; ++i)
      I_matrix[i][i] = 1;

   return I_matrix;
}

/****************************************************************************************/

template <typename T> T merge(const T &source1, const T &source2, matrix::index_t n) {

   T destination(n);

   for (matrix::index_t i = 0; i < n / 2; ++i)
      for (matrix::index_t j = 0; j < n; ++j)
         destination[i][j] = source1[i][j];

   for (matrix::index_t i = 0; i < n / 2; ++i)
      for (matrix::index_t j = 0; j < n; ++j)
         destination[i + n / 2][j] = source2[i][j];

   return destination;
}

/****************************************************************************************/

// n should be a power of 2
matrix::dense<matrix::entry_t> haar_matrix(matrix::index_t n) {

   if (n == 2) {
      matrix::dense<matrix::entry_t> m(2);
      m[0][0] = m[0][1] = m[1][0] = 1;
      m[1][1] = -1;
      return m;
   } else {
      matrix::dense<matrix::entry_t> mm(1, 2); // 1x2 matrix
      matrix::dense<matrix::entry_t> mn(1, 2); // 1x2 matrix
      mm[0][0] = mn[0][0] = mm[0][1] = 1;
      mn[0][1] = -1;

      return merge(kronecker_product(haar_matrix(n / 2), mm, n / 2, n / 2, 1, 2),
                   kronecker_product(identity(n / 2), mn, n / 2, n / 2, 1, 2), n);
   }
}

/****************************************************************************************/

template <typename T>
void plant_entries(T &A_matrix, T &C_matrix, matrix::index_t n, matrix::index_t k) {

   std::vector<unsigned int> buff(n);
   std::iota(std::begin(buff), std::end(buff), 1);

   std::random_shuffle(std::begin(buff), std::end(buff));

   for (matrix::entry_t ent = 0; ent < k; ++ent) {
      matrix::entry_t r = buff[ent];
      for (matrix::index_t i = 0; i < n; ++i) {
         if (A_matrix[r][i] != 0 && C_matrix[i][r] != 0) {
            A_matrix[r][i] = 0;
            break;
         }
         if (A_matrix[r][i] != 0 && C_matrix[i][r] == 0) {
            C_matrix[r][i] = 1;
            break;
         }
         if (A_matrix[r][i] == 0 && C_matrix[i][r] != 0) {
            A_matrix[r][i] = 1;
            break;
         }
      }
   }
}

/****************************************************************************************/

std::pair<matrix::dense<matrix::entry_t>, matrix::dense<matrix::entry_t>>
generate_matrix(matrix::index_t n, matrix::index_t d,
                __attribute__((unused)) matrix::index_t k) {

   matrix::dense<matrix::entry_t> matrix = sparsegen::haar_matrix(d);
   matrix::dense<matrix::entry_t> A_matrix(n);
   matrix::dense<matrix::entry_t> C_matrix(n);

   for (matrix::index_t i = 0; i < d / 2; ++i) {
      for (matrix::index_t j = 0; j < d; ++j) {
         A_matrix[i][j] = matrix[i][j];
         A_matrix[i + d / 2][j] = matrix[i][j];
         C_matrix[j][i] = matrix[i + d / 2][j];
         C_matrix[j][i + d / 2] = matrix[i + d / 2][j];
      }
   }

   for (matrix::index_t rep = 0; rep < n / d - 1; ++rep) {
      for (matrix::index_t i = rep * d; i < rep * d + d; ++i) {
         for (matrix::index_t j = rep * d; j < rep * d + d; ++j) {
            A_matrix[i + d][j + d] = A_matrix[i][j];
            C_matrix[i + d][j + d] = C_matrix[i][j];
         }
      }
   }
   utils::randomize_matrix(A_matrix, n);
   utils::randomize_matrix(C_matrix, n);
   plant_entries(A_matrix, C_matrix, n, 2 * k / d);
   return std::pair<matrix::dense<matrix::entry_t>, matrix::dense<matrix::entry_t>>(
       A_matrix, C_matrix);
}

/****************************************************************************************/
// TO DO: RASMUS GENERATOR
// template <typename T>

} // namespace sparsegen
