
#include <algorithm>
#include <cstdlib>
// #include <cmath>

namespace cascading {

// static const size_t INHERITED = 0
// static const size_t NATIVE = 1

/****************************************************************************************/

template <typename T> class augmented_entry {

 public:
   typedef augmented_entry<T> S;
   augmented_entry() {
      this->entry = nullptr;
      this->predecessor = nullptr;
      this->successor = nullptr;
      this->bridge = nullptr;
      this->entry_type = utils::INHERITED;
   };
   ~augmented_entry(){};

   // copy overload, augmented_entry = augmented_entry
   S &operator=(const S &rhs) {
      if (this->entry == nullptr)
         this->entry = new T();
      memcpy(this->entry, rhs.entry, sizeof(utils::entry<T>));
      return *this; // return the result by reference
   }

   // copy overload, augmented_entry = entry
   S &operator=(const T &rhs) {
      if (this->entry == nullptr)
         this->entry = new T();
      memcpy(this->entry, &rhs, sizeof(utils::entry<T>));
      this->entry_type = utils::NATIVE;
      return *this; // return the result by reference
   }

   utils::entry<utils::e_type> *entry, *predecessor, *successor;
   const S *bridge;
   bool entry_type; // 0 inherited, 1 native

   // void set_entry(T *entry) { this->entry = entry; }
   void set_predecessor(const T &predecessor) {
      this->predecessor = &predecessor;
   }
   void set_successor(const T &successor) { this->successor = &successor; }
   void set_bridge(const S &bridge) { this->bridge = &bridge; }
   // void set_entry_type(bool entry_type) { this->entry_type = entry_type; }

   // void set_entry(T &entry, bool entry_type) {
   //    this->entry = &entry;
   //    this->entry_type = entry_type;
   // }
};

/****************************************************************************************/

template <typename T, typename S>
inline cascading::augmented_entry<utils::entry<utils::e_type>> *
augment_matrix(const T &matrix, const S &pointers, size_t h, size_t n,
               size_t order) {

   typedef cascading::augmented_entry<utils::entry<utils::e_type>> ae_type;

   // sort_by returns either row or column index depending on order
   utils::sort_by<utils::entry<utils::e_type>> index_type;

   ae_type *index_vector = new ae_type[n]();
   ae_type *augmented_matrix = new ae_type[3 * h]();
   size_t *augmented_pointers = new size_t[n]();

   size_t even;
   size_t matr_index = pointers[1]; // n-1 column
   size_t curr_index(0), prev_index, curr_position, bridge_pt;

   // copy the last column in the augmented_matrix
   // last column is from pointers[1] to h=nnz(A)
   while (matr_index < h)
      augmented_matrix[curr_index++] = matrix[matr_index++];

   // column n-1 of augmented_matrix starts at curr_index
   augmented_pointers[1] = curr_index; // n-1 is 0

   for (size_t index = 1; index < n; ++index) {   // row or column index
      prev_index = augmented_pointers[index - 1]; // pt for inherited entries
      matr_index = pointers[index + 1];           // pt for native entries
                                        // nb pointers are in descending order
      // every column starts with even = false
      even = false;

      // the predecessor in row/column j is the first element in row/column j+1
      bridge_pt = prev_index;

      // current index from pointers[index] to pointers[index + 1]
      // previous index from augmented_pointers[index - 1] to
      // augmented_pointers[index]
      while ((prev_index < augmented_pointers[index]) ||
             (matr_index < pointers[index])) {

         assert(augmented_matrix[prev_index].entry != nullptr);

         // curr_position is the minimum index between
         // augmented_matrix and matrix in current row/column
         curr_position = std::min(
             prev_index < augmented_pointers[index] &&
                     augmented_matrix[prev_index].entry != nullptr
                 ? index_type(*augmented_matrix[prev_index].entry, order)
                 : n + 1,
             matr_index < pointers[index]
                 ? index_type(matrix[matr_index], order)
                 : n + 1);

         assert(0 <= curr_position && curr_position <= n);

         // the entry is native
         if (index_type(matrix[matr_index], order) == curr_position &&
             matr_index < pointers[index]) {
            augmented_matrix[curr_index] = matrix[matr_index];

#ifdef DEBUG_VERBOSE
            printf("(%ld,%ld,%d) native added (%ld,%ld,%d) at %ld\n",
                   matrix[matr_index].i, matrix[matr_index].j,
                   matrix[matr_index].a, augmented_matrix[curr_index].entry->i,
                   augmented_matrix[curr_index].entry->j,
                   augmented_matrix[curr_index].entry->a, curr_index);
#endif
            // if matrix and augmented_matrix have an entry at i
            // increase counter for augmented_matrix
            //  obj.entry->i can segfault if obj is nullptr
            if (augmented_matrix[prev_index].entry != nullptr &&
                index_type(*augmented_matrix[prev_index].entry, order) ==
                    index_type(matrix[matr_index], order)) {
               even = !even;
               prev_index++;
            }
            // curr_index++;
            matr_index++;
            // otherwise the entry is inherited
         } else if (index_type(*augmented_matrix[prev_index].entry, order) ==
                        curr_position &&
                    prev_index < augmented_pointers[index]) {
            if (even) {
               augmented_matrix[curr_index] = augmented_matrix[prev_index];
               // decreasing index means copy to previous row/column
               order ? augmented_matrix[curr_index].entry->set_j(
                           index_type(*augmented_matrix[prev_index].entry,
                                      utils::COLUMN) -
                           1)
                     : augmented_matrix[curr_index].entry->set_i(
                           index_type(*augmented_matrix[prev_index].entry,
                                      utils::ROW) -
                           1);
               // add bridge to next row column 
                   if (index_type(*augmented_matrix[curr_index].entry,
                                         order) ==
                              index_type(*augmented_matrix[bridge_pt].entry,
                                         order)) augmented_matrix[curr_index]
                       .set_bridge(augmented_matrix[bridge_pt]);
               if (index_type(*augmented_matrix[curr_index].entry, order) ==
                   index_type(*augmented_matrix[bridge_pt + 1].entry, order))
                  augmented_matrix[curr_index].set_bridge(
                      augmented_matrix[++bridge_pt]);
               else
                  bridge_pt++;
#ifdef DEBUG_VERBOSE
               printf("(%ld,%ld,%d) inherited added (%ld,%ld,%d) at %ld\n",
                      augmented_matrix[prev_index].entry->i,
                      augmented_matrix[prev_index].entry->j,
                      augmented_matrix[prev_index].entry->a,
                      augmented_matrix[curr_index].entry->i,
                      augmented_matrix[curr_index].entry->j,
                      augmented_matrix[curr_index].entry->a, curr_index);
#endif
            }
            curr_index++;
            // prev_index and even are updated regardless of assignment
            prev_index++;
            even = !even;
         }
#ifdef DEBUG_VERBOSE
         utils::print_sparse_as_dense(augmented_matrix, h, n);
#endif
      }
      // pointer to beginning of current row/column in augmented_matrix
      augmented_pointers[index + 1] = curr_index;
   }

   // compute the bridges between index_vector and
   // first row/column of augmented_matrix
   // for (ssize_t index = n - 1; index >= 0; --index) {
   size_t index = n;
   bridge_pt = augmented_pointers[n] - 1;
   while (index > 0 && bridge_pt >= augmented_pointers[n - 1]) {
      if (index_type(*augmented_matrix[bridge_pt].entry, order) > index - 1)
         index_vector[index - 1].set_bridge(augmented_matrix[--bridge_pt]);
      else
         index_vector[index - 1].set_bridge(augmented_matrix[bridge_pt]);
      index--;
      // (index_type(*augmented_matrix[bridge_pt].entry, order) <= index)

      // bridge_pt = augmented_pointers[n - 1];
      // for (size_t index = 0; index < n; ++index) {
      //    if (index_type(*augmented_matrix[bridge_pt].entry, order) ==
      //    index)
      //       index_vector[index].set_bridge(augmented_matrix[bridge_pt]);
      //    if (bridge_pt + 1 < augmented_pointers[n] &&
      //        index_type(*augmented_matrix[bridge_pt].entry, order) <
      //        index
      //        && index_type(*augmented_matrix[bridge_pt + 1].entry,
      //        order)
      //        == index)
      //       index_vector[index].set_bridge(augmented_matrix[++bridge_pt]);
      //    if (bridge_pt + 1 == augmented_pointers[n] &&
      //        index_type(*augmented_matrix[bridge_pt].entry, order) <
      //        index)
      //       index_vector[index].set_bridge(augmented_matrix[bridge_pt]);
#ifdef DEBUG_VERBOSE
      if (index_vector[index].bridge == nullptr)
         printf("%ld -> nullptr\n", index);
      else
         printf("%ld -> (%ld,%ld,%d)\n", index,
                index_vector[index].bridge->entry->i,
                index_vector[index].bridge->entry->j,
                index_vector[index].bridge->entry->a);
#endif
   }

   return index_vector;
}

/****************************************************************************************/

template <typename T>
inline cascading::augmented_entry<utils::entry<utils::e_type>> *
augment_matrix(const T &matrix, size_t h, size_t n, size_t order) {

   // sort_by returns either row or column index depending on order
   utils::sort_by<utils::entry<utils::e_type>> index_type;

   // pointers store a pointer to the beginning of each row/column
   // in decreasing order
   size_t *pointers = new size_t[n + 1]();

   // 1 - ROW = COLUMN
   // 1 - COLUMN = ROW
   size_t pt = 0;
   for (size_t index = 0; index < h; ++index) {
      if (index_type(matrix[index], 1 - order) != pt) {
         pointers[n - pt - 1] = index;
         pt = index_type(matrix[index], 1 - order);
      }
   }
   pointers[0] = h;

   return cascading::augment_matrix(matrix, pointers, h, n, order);
}

/****************************************************************************************/

template <typename T>
inline utils::e_type *extract_sketch(const T &index_vector, size_t i1,
                                     size_t i2, size_t n) {

   utils::e_type *sketch = new utils::e_type[n]();
   T ei1 = index_vector[std::min(i1 - 1, (size_t)0)]; // easier with i1-1
   T ei2 = index_vector[i2];

   for (size_t index = 0; index < n + 1; ++index) {
      ei1 = ei1->bridge;
      ei2 = ei2->bridge;
      sketch[index] = ei2->a - ei1->a;
   }
   return sketch;
}

/****************************************************************************************/

template <typename T>
inline void matrix_product(const T &A_index_vector, const T &C_index_vector,
                           size_t i1, size_t i2, size_t j1, size_t j2,
                           size_t n) {

   size_t im, jm;

   if (i1 == i2 && j1 == j2) {
      utils::inner_product(
          cascading::extract_sketch(A_index_vector, i1, i2, n),
          cascading::extract_sketch(C_index_vector, j1, j2, n));
   } else if (i1 == i2 && j1 < j2) {
      jm = j1 + j2 / 2;
      if (utils::inner_product(
              cascading::extract_sketch(A_index_vector, i1, i2, n),
              cascading::extract_sketch(C_index_vector, j1, jm, n)) != 0)
         cascading::matrix_product(A_index_vector, C_index_vector, i1, i2, j1,
                                   jm, n);
      if (utils::inner_product(
              cascading::extract_sketch(A_index_vector, i1, i2, n),
              cascading::extract_sketch(C_index_vector, jm, j2, n)) != 0)
         cascading::matrix_product(A_index_vector, C_index_vector, i1, i2, j1,
                                   jm, n);
   } else if (i1 < i2 && j1 == j2) {
      im = i1 + i2 / 2;
      if (utils::inner_product(
              cascading::extract_sketch(A_index_vector, i1, im, n),
              cascading::extract_sketch(C_index_vector, j1, j2, n)) != 0)
         cascading::matrix_product(A_index_vector, C_index_vector, i1, im, j1,
                                   j2, n);
      if (utils::inner_product(
              cascading::extract_sketch(A_index_vector, im, i2, n),
              cascading::extract_sketch(C_index_vector, j1, j2, n)) != 0)
         cascading::matrix_product(A_index_vector, C_index_vector, im, i2, j1,
                                   j2, n);
   } else {
      im = i1 + i2 / 2;
      jm = j1 + j2 / 2;
      if (utils::inner_product(
              cascading::extract_sketch(A_index_vector, i1, im, n),
              cascading::extract_sketch(C_index_vector, j1, jm, n)) != 0)
         cascading::matrix_product(A_index_vector, C_index_vector, i1, im, j1,
                                   jm, n);
      if (utils::inner_product(
              cascading::extract_sketch(A_index_vector, i1, im, n),
              cascading::extract_sketch(C_index_vector, jm, j2, n)) != 0)
         cascading::matrix_product(A_index_vector, C_index_vector, i1, im, jm,
                                   j2, n);
      if (utils::inner_product(
              cascading::extract_sketch(A_index_vector, im, i2, n),
              cascading::extract_sketch(C_index_vector, j1, jm, n)) != 0)
         cascading::matrix_product(A_index_vector, C_index_vector, im, i2, j1,
                                   jm, n);
      if (utils::inner_product(
              cascading::extract_sketch(A_index_vector, im, i2, n),
              cascading::extract_sketch(C_index_vector, jm, j2, n)) != 0)
         cascading::matrix_product(A_index_vector, C_index_vector, im, i2, jm,
                                   j2, n);
   }
}

/****************************************************************************************/

// template <typename T>
// inline void matrix_product(const T &A_matrix, const T &C_matrix, size_t n)
// {}

} // namespace cascading