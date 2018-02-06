

namespace sparsegen {

template <typename T>
inline T kronecker_product(const T &A_matrix, const T &B_matrix, size_t A_m,
                           size_t A_n, size_t B_m, size_t B_n) {

   size_t m = A_m * B_m;
   size_t n = A_n * B_n;

   utils::e_type **AB_matrix = new utils::e_type *[m];

   for (size_t index = 0; index < m; ++index)
      AB_matrix[index] = new utils::e_type[n]();

   for (size_t A_i = 0; A_i < A_m; ++A_i)
      for (size_t A_j = 0; A_j < A_n; ++A_j)
         for (size_t B_i = 0; B_i < B_m; ++B_i)
            for (size_t B_j = 0; B_j < B_n; ++B_j)
               AB_matrix[A_i * B_m + B_i][A_j * B_n + B_j] =
                   A_matrix[A_i][A_j] * B_matrix[B_i][B_j];
   return AB_matrix;
}

/****************************************************************************************/

utils::e_type **identity(size_t n) {

   utils::e_type **I_matrix = new utils::e_type *[n];

   for (size_t i = 0; i < n; ++i) {
      I_matrix[i] = new utils::e_type[n]();
      I_matrix[i][i] = 1;
   }
   return I_matrix;
}

/****************************************************************************************/

template <typename T> T merge(const T &source1, const T &source2, size_t n) {

   utils::e_type **destination = new utils::e_type *[n];

   for (size_t index = 0; index < n; ++index)
      destination[index] = new utils::e_type[n]();

   for (size_t i = 0; i < n / 2; ++i)
      for (size_t j = 0; j < n; ++j)
         destination[i][j] = source1[i][j];

   for (size_t i = 0; i < n / 2; ++i)
      for (size_t j = 0; j < n; ++j)
         destination[i + n / 2][j] = source2[i][j];

   return destination;
}

/****************************************************************************************/

utils::e_type **haar_matrix(size_t n) {

   if (n == 2) {
      utils::e_type **m = new utils::e_type *[2];
      m[0] = new utils::e_type[2];
      m[1] = new utils::e_type[2];
      m[0][0] = m[0][1] = m[1][0] = 1;
      m[1][1] = -1;
      return m;
   } else {

      utils::e_type **mm = new utils::e_type *[1];
      utils::e_type **mn = new utils::e_type *[1];
      mm[0] = new utils::e_type[2];
      mn[0] = new utils::e_type[2];
      mm[0][0] = mn[0][0] = mm[0][1] = 1;
      mn[0][1] = -1;

      return sparsegen::merge(
          sparsegen::kronecker_product(sparsegen::haar_matrix(n / 2), mm, n / 2,
                                       n / 2, 1, 2),
          sparsegen::kronecker_product(sparsegen::identity(n / 2), mn, n / 2,
                                       n / 2, 1, 2),
          n);
   }
}

/****************************************************************************************/

template <typename T>
void plant_entries(const T &A_matrix, const T &C_matrix, size_t n, size_t k) {

   size_t *buff = new size_t[n];

   for (size_t index = 0; index < n; ++index)
      buff[index] = index;

   std::random_shuffle(buff, buff + n);

   for (size_t ent = 0; ent < k; ++ent) {
      size_t r = buff[ent];
      for (size_t i = 0; i < n; ++i) {
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
   delete[] buff;
}

/****************************************************************************************/

template <typename T>
void generate_matrix(const T &A_matrix, const T &C_matrix, size_t n, size_t d,
                     size_t k) {

   T matrix = sparsegen::haar_matrix(d);

   for (size_t i = 0; i < d / 2; ++i) {
      for (size_t j = 0; j < d; ++j) {
         A_matrix[i][j] = matrix[i][j];
         A_matrix[i + d / 2][j] = matrix[i][j];
         C_matrix[j][i] = matrix[i + d / 2][j];
         C_matrix[j][i + d / 2] = matrix[i + d / 2][j];
      }
   }

   for (size_t rep = 0; rep < n / d - 1; ++rep) {
      for (size_t i = rep * d; i < rep * d + d; ++i) {
         for (size_t j = rep * d; j < rep * d + d; ++j) {
            A_matrix[i + d][j + d] = A_matrix[i][j];
            C_matrix[i + d][j + d] = C_matrix[i][j];
         }
      }
   }
   for (size_t i = 0; i < d / 2; ++i)
      delete[] matrix[i];
   delete[] matrix;

   // utils::randomize_matrix(A_matrix, C_matrix, n);
   // sparsegen::plant_entries(A_matrix, C_matrix, n, 2 * k / d);
}

/****************************************************************************************/
//TO DO: RASMUS GENERATOR
// template <typename T>
// void rasmus_generator(const T &A_matrix, const T &C_matrix, size_t n, size_t
// d,
//                       size_t k) {

// p = int(sys.argv[1])
// q = int(sys.argv[2])
//
// def indicator(i,p):
//     for j in xrange(p):
//         if j==i: print q,
//         else: print 0,
//
// for (size_t i = 0; i < p; ++i1) {
//    for (size_t i = 0; i < p; ++i2) {
//       A_matrix[i][j] =
//    }
// }

// for i1 in xrange(p):
//     for i2 in xrange(p):
//         print (p-1)*"1 " + "1",
//         indicator(i1,p)
//         for j in xrange(p):
//             indicator((j*i1+i2)%p,p)
//         print ""
// }

} // namespace sparsegen
