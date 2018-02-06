
#include <algorithm>
#include <cstdlib>

namespace cascading {

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

   utils::ee *entry, *predecessor, *successor;
   const S *bridge;
   bool entry_type; // 0 utils::INHERITED, 1 utils::NATIVE

   void set_predecessor(T &predecessor) { this->predecessor = &predecessor; }
   void set_successor(T &successor) { this->successor = &successor; }
   void set_bridge(const S &bridge) { this->bridge = &bridge; }
};

/****************************************************************************************/

template <typename T, typename S>
inline cascading::augmented_entry<utils::ee> *
augment_matrix(const T &matrix, const S &pointers, size_t h, size_t n,
               size_t order) {

   // sort_by returns either row or column index depending on order
   utils::sort_by<utils::ee> index_type;

   auto *index_vector = new cascading::augmented_entry<utils::ee>[n]();
   auto *augmented_matrix = new cascading::augmented_entry<utils::ee>[3 * h]();
   size_t *augmented_pointers = new size_t[n]();

   size_t prev_index, curr_position, bridge_pt, native_pt, inherited_pt, even,
       matr_index(pointers[1]), curr_index(0);

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
      bridge_pt = augmented_pointers[index - 1];
      native_pt = -1;
      inherited_pt = -1;

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

            // between two nonconsecutive native entries (one or more inherited
            // entries may be inbetween) we add pointers to
            // predecessor/successor
            if (native_pt > 0 && curr_index != native_pt + 1) {
               augmented_matrix[curr_index].set_predecessor(
                   *augmented_matrix[native_pt].entry);
               augmented_matrix[native_pt].set_successor(
                   *augmented_matrix[curr_index].entry);
            }
            native_pt = curr_index;

            // add pointers from the previous inherited entries to
            // the new native entry
            if (inherited_pt > 0)
               for (size_t ii = curr_index - 1; ii >= pointers[index + 1];
                    --ii) {
                  if (augmented_matrix[ii].entry_type == utils::NATIVE)
                     break;
                  else
                     augmented_matrix[ii].set_successor(
                         *augmented_matrix[curr_index].entry);
               }

            // add bridge to next row/column
            // bridge is from augmented_matrix[curr_index] to
            // augmented_matrix[bridge_pt]
            if (bridge_pt + 1 < augmented_pointers[index] &&
                index_type(*augmented_matrix[curr_index].entry, order) ==
                    index_type(*augmented_matrix[bridge_pt + 1].entry, order))
               augmented_matrix[curr_index].set_bridge(
                   augmented_matrix[++bridge_pt]);
            else if (index_type(*augmented_matrix[curr_index].entry, order) >=
                     index_type(*augmented_matrix[bridge_pt].entry, order))
               augmented_matrix[curr_index].set_bridge(
                   augmented_matrix[bridge_pt]);

#ifdef DEBUG_VERBOSE
            printf("(%ld,%ld,%d) native added (%ld,%ld,%d) at %ld\n",
                   matrix[matr_index].i, matrix[matr_index].j,
                   matrix[matr_index].a, augmented_matrix[curr_index].entry->i,
                   augmented_matrix[curr_index].entry->j,
                   augmented_matrix[curr_index].entry->a, curr_index);
            if (augmented_matrix[curr_index].bridge != nullptr)
               printf("(%ld,%ld,%d) bridge (%ld,%ld,%d)\n",
                      augmented_matrix[curr_index].entry->i,
                      augmented_matrix[curr_index].entry->j,
                      augmented_matrix[curr_index].entry->a,
                      augmented_matrix[curr_index].bridge->entry->i,
                      augmented_matrix[curr_index].bridge->entry->j,
                      augmented_matrix[curr_index].bridge->entry->a);
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
            curr_index++;
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

               if (native_pt > 0)
                  augmented_matrix[curr_index].set_predecessor(
                      *augmented_matrix[native_pt].entry);

               inherited_pt = curr_index;

               // add bridge to next row/column
               // bridge is from augmented_matrix[curr_index] to
               // augmented_matrix[bridge_pt]
               if (bridge_pt + 1 < augmented_pointers[index] &&
                   index_type(*augmented_matrix[curr_index].entry, order) ==
                       index_type(*augmented_matrix[bridge_pt + 1].entry,
                                  order))
                  augmented_matrix[curr_index].set_bridge(
                      augmented_matrix[++bridge_pt]);
               else if (index_type(*augmented_matrix[curr_index].entry,
                                   order) >=
                        index_type(*augmented_matrix[bridge_pt].entry, order))
                  augmented_matrix[curr_index].set_bridge(
                      augmented_matrix[bridge_pt]);
#ifdef DEBUG_VERBOSE
               printf("(%ld,%ld,%d) inherited added (%ld,%ld,%d) at %ld\n",
                      augmented_matrix[prev_index].entry->i,
                      augmented_matrix[prev_index].entry->j,
                      augmented_matrix[prev_index].entry->a,
                      augmented_matrix[curr_index].entry->i,
                      augmented_matrix[curr_index].entry->j,
                      augmented_matrix[curr_index].entry->a, curr_index);
               if (augmented_matrix[curr_index].bridge != nullptr)
                  printf("(%ld,%ld,%d) bridge (%ld,%ld,%d)\n",
                         augmented_matrix[curr_index].entry->i,
                         augmented_matrix[curr_index].entry->j,
                         augmented_matrix[curr_index].entry->a,
                         augmented_matrix[curr_index].bridge->entry->i,
                         augmented_matrix[curr_index].bridge->entry->j,
                         augmented_matrix[curr_index].bridge->entry->a);
#endif
               curr_index++;
            }
            // prev_index and even are updated regardless of assignment
            prev_index++;
            even = !even;
         }
#ifdef DEBUG_VERBOSE
         utils::print_sparse_as_dense(augmented_matrix, 3 * h, n);
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
inline cascading::augmented_entry<utils::ee> *
augment_matrix(const T &matrix, size_t h, size_t n, size_t order) {

   // sort_by returns either row or column index depending on order
   utils::sort_by<utils::ee> index_type;

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
                                     size_t i2, size_t n, size_t order) {

   utils::e_type *sketch = new utils::e_type[n]();
   auto ei2 = index_vector[i2];

   for (size_t index = 0; index < n; ++index) {
      ei2 = *ei2.bridge;
      sketch[index] = (&ei2)->entry->a;
   }

   i1 = std::max((ssize_t)(i1 - 1), (ssize_t)0); // easier with i1-1
   if (i1 == 0)
      return sketch;
   auto ei1 = index_vector[i1];

   for (size_t index = 0; index < n; ++index) {
      ei1 = *ei1.bridge;
      sketch[index] -= (&ei1)->entry->a;
   }
   for (size_t index = 0; index < n; ++index)
      printf("%d ", sketch[index]);
   printf("\n");
   // int c = getchar();
   return sketch;
}

} // namespace cascading