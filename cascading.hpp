
#include <algorithm>
#include <cstdlib>

namespace cascading {

/****************************************************************************************/

template <class T> class augmented_entry {

 public:
   typedef augmented_entry<T> S;

   const T *entry, *predecessor;
   S *bridge;
   bool entry_type; // 0 utils::INHERITED, 1 utils::NATIVE

   augmented_entry()
       : entry(nullptr), predecessor(nullptr), bridge(nullptr),
         entry_type(utils::NATIVE) {}
   augmented_entry(const T &rhs_entry)
       : entry(&rhs_entry), predecessor(nullptr), bridge(nullptr),
         entry_type(utils::NATIVE) {}
   ~augmented_entry(){
       // delete entry;
       // delete predecessor;
       // delete successor;
       // delete bridge;
   };

   // copy overload, augmented_entry = augmented_entry
   S &operator=(const S &rhs) {
      // std::cout << "augmented_entry = augmented_entry" << std::endl;  // for debug
      memcpy(this, &rhs, sizeof(S));
      return *this; // return the result by reference
   }
   // copy overload, augmented_entry = entry
   S &operator=(const T &rhs) {
      // std::cout << "augmented_entry = entry" << std::endl;  // for debug
      memcpy(this->entry, &rhs, sizeof(T));
      this->entry_t = utils::NATIVE;
      return *this; // return the result by reference
   }

   // augmented_entry(const S &&rhs) noexcept
   //     : std::exchange(rhs, nullptr) {} // move constructor
   //
   // S &operator=(const S &&rhs) noexcept { // move assignment
   //    std::swap(bridge, rhs);
   //    return *this;
   // }
}; // namespace cascading

template <typename T> class augmented_matrix {

 public:
   std::vector<augmented_entry<T> *> index_vector;
   std::vector<std::vector<augmented_entry<T>>> matrix;

   augmented_matrix(matrix::index_t n)
       : index_vector(n), matrix(n, std::vector<augmented_entry<T>>(0)) {} // constructor
   ~augmented_matrix() {}                                                  // destructor
   // augmented_matrix(const augmented_matrix &other) {}                 // copy
   // constructor augmented_matrix(augmented_matrix &&other) noexcept {}             //
   // move constructor augmented_matrix &operator=(const coalesced_buckets &other) {}
   // // copy assignment augmented_matrix &operator=(coalesced_buckets &&other) noexcept
   // {} // move assignment
};

/****************************************************************************************/

template <typename T>
augmented_matrix<matrix::coo_entry> augment_matrix(const T &matrix, matrix::index_t n) {

   // return_index returns either row or column index depending on matrix.layout
   utils::return_index<matrix::coo_entry> matrix_index;

   // data structure representing the augmented matrix
   augmented_matrix<matrix::coo_entry> ds(n);

   // the last vector is the same as the original matrix
   for (auto &e : matrix[n - 1])
      ds.matrix[n - 1].push_back(e);

   // we start from row/column n-2 and we inherit from row/column n-1 down to row/column 0
   for (matrix::index_t index = n - 2; index >= 0; --index) {

      auto curr_pt = std::begin(matrix[index]); // pointer to current row/col vector
      auto next_pt = std::begin(ds.matrix[index + 1]); // pointer to next row/col vector

      bool even = false; // we only inherit entries with even rank
      while ((curr_pt != std::end(matrix[index])) &&
             (next_pt != std::end(ds.matrix[index + 1]))) {

         // if we are ahead with inherited entries
         if (matrix_index(*curr_pt, 1 - matrix.layout) <
             matrix_index(*next_pt->entry, 1 - matrix.layout)) {
            augmented_entry<matrix::coo_entry> et(*curr_pt);
            if (next_pt != std::begin(ds.matrix[index + 1]))
               et.bridge = &(*(next_pt - 1));
            ds.matrix[index].push_back(et);
            ++curr_pt;

            // if we are at the same level
         } else if (matrix_index(*curr_pt, 1 - matrix.layout) ==
                    matrix_index(*next_pt->entry, 1 - matrix.layout)) {
            augmented_entry<matrix::coo_entry> et(*curr_pt);
            et.bridge = &(*next_pt);
            ds.matrix[index].push_back(et);
            even = !even;
            ++next_pt;
            ++curr_pt;

            // otherwise we are behind with inherited entries
         } else {
            if (even) {
               augmented_entry<matrix::coo_entry> et(*next_pt);
               et.bridge = &(*next_pt);
               et.predecessor = &(*curr_pt);
               et.entry_type = utils::INHERITED;
               ds.matrix[index].push_back(et);
            }
            ++next_pt;
            even = !even;
         }
      }
      // after the while loop curr_pt == std::end(matrix[index]) or
      // next_pt == std::end(ds.matrix[index + 1]) (or both). Therefore,
      // only once of the following for-loops will be executed

      // flush the remaining native entries into the current row/column
      for (; curr_pt != std::end(matrix[index]); ++curr_pt) {
         augmented_entry<matrix::coo_entry> et(*curr_pt);
         et.bridge = &(*(next_pt - 1)); // next_pt is std::end(matrix[index])
         ds.matrix[index].push_back(et);
      }

      // flush the remaining inherited entries of even rank into the current row/column
      for (; next_pt != std::end(ds.matrix[index + 1]); ++next_pt) {
         if (even) {
            augmented_entry<matrix::coo_entry> et(*next_pt->entry);
            et.bridge = &(*next_pt);
            et.predecessor = &(*(curr_pt - 1)); // curr_pt is std::end(matrix[index])
            et.entry_type = utils::INHERITED;
            ds.matrix[index].push_back(et);
         }
         even = !even;
      }
      // no more entries (neither native nor inherited) left
      assert((curr_pt == std::end(matrix[index])) &&
             (next_pt == std::end(ds.matrix[index + 1])));
   }

   // index vector indices the first row/colum of the data structure
   auto next_pt = std::begin(ds.matrix[0]);
   for (matrix::index_t index = 0; index < n; ++index)
      if (index <= matrix_index(*next_pt->entry, 1 - matrix.layout)) {
         ds.index_vector[index] = &(*next_pt);
         if (next_pt + 1 != std::end(ds.matrix[0]))
            ++next_pt;
      }

   return ds;
}

/****************************************************************************************/

template <typename T>
inline std::vector<matrix::entry_t> extract_sketch(T &ds, matrix::index_t i1,
                                                   matrix::index_t i2, matrix::index_t n,
                                                   ssize_t layout) {

   utils::return_index<matrix::coo_entry> matrix_index;

   assert(i2 >= i1);
   std::vector<matrix::entry_t> sketch(n, 0);

   i1 = std::max((i1 - 1), 0); // easier with i1-1

   auto i2_pt = ds.index_vector[i2]; // index matrix in i2
   auto i1_pt = ds.index_vector[i1]; // index matrix in i1

   for (matrix::index_t index = 0; index < n; ++index) {

      if (i2_pt == nullptr) {          // if there is no bridge to i2 (i2_pt == null)
         i2_pt = &ds.matrix[index][0]; // we consider the beginning of the vector
      } else {
         if (matrix_index(*i2_pt->entry, layout) <= i2) {
            if (i2_pt->entry_type == utils::NATIVE) { // if it's native
               sketch[index] += i2_pt->entry->a;      // just take the value
            } else {                                  // otherwise it's an inherited entry
               if (i2_pt->predecessor != nullptr)     // take the native predecessor
                  sketch[index] += i2_pt->predecessor->a;
            }
         }
      }

      if (i1_pt == nullptr) {
         i1_pt = &ds.matrix[index][0]; // we consider the beginning of the vector
      } else {
         if (matrix_index(*i1_pt->entry, layout) <= i1) {
            if (i1_pt->entry_type == utils::NATIVE) { // if it's native
               sketch[index] -= i1_pt->entry->a;      // just take the value
            } else {                                  // otherwise it's an inherited entry
               if (i1_pt->predecessor != nullptr)     // take the native predecessor
                  sketch[index] -= i1_pt->predecessor->a;
            }
         }
      }

      assert(i2_pt != nullptr && i1_pt != nullptr);

      i2_pt = i2_pt->bridge; // jump to next row/column
      i1_pt = i1_pt->bridge; // jump to next row/column
   }

#ifndef DEBUG_VERBOSE
   for (auto &x : sketch)
      std::cout << x << " ";
   std::cout << std::endl;
#endif
   return sketch;
}

/****************************************************************************************/

// std::vector<cascading::augmented_entry<matrix::entry_t> *> index_vector(n);
// auto *index_vector = new cascading::augmented_entry<utils::ee>[n]();
// auto *augmented_matrix = new cascading::augmented_entry<utils::ee>[3 * h]();
// size_t *augmented_pointers = new size_t[n + 1]();
// std::vector<matrix::index_t> augmented_pointers(n + 1, 0);
//
// size_t prev_index, curr_position, bridge_pt, even, matr_index(pointers[1]),
//     curr_index(0);
// ssize_t native_pt, inherited_pt;

// copy the last column in the augmented_matrix
// last column is from pointers[1] to h=nnz(A)

//    while (matr_index < h)
//       ds.matrix.push_back(matrix[matr_index++]);
//
//    // column n-1 of augmented_matrix starts at curr_index
//    augmented_pointers[1] = curr_index;
//
//    for (size_t index = 1; index < n; ++index) {   // row or column index
//       prev_index = augmented_pointers[index - 1]; // pt for inherited entries
//       matr_index = pointers[index + 1];           // pt for native entries
//                                                   // nb pointers are in descending
//                                                   order
//
//       // every column starts with even = false
//       even = false;
//
//       // the predecessor in row/column index
//       // is the first element in row/column index+1
//       bridge_pt = augmented_pointers[index - 1];
//       native_pt = -1;
//       inherited_pt = -1;
//
//       // current index from pointers[index] to pointers[index + 1]
//       // previous index from augmented_pointers[index - 1] to
//       // augmented_pointers[index]
//       while ((prev_index < augmented_pointers[index]) || (matr_index <
//       pointers[index])) {
//
//          assert(ds.matrix[index][prev_index].entry != nullptr);
//
//          // curr_position is the minimum index between
//          // augmented_matrix and matrix in current row/column
//          curr_position = std::min(
//              prev_index < augmented_pointers[index]
//                  ? matrix_index(*ds.matrix[index][prev_index].entry, matrix.layout)
//                  : n + 1,
//              matr_index < pointers[index]
//                  ? matrix_index(matrix[matr_index], matrix.layout)
//                  : n + 1);
//
//          assert(0 <= curr_position && curr_position <= n);
//
//          // the entry is native
//          if (matrix_index(matrix[matr_index], matrix.layout) == curr_position &&
//              matr_index < pointers[index]) {
//             ds.matrix[curr_index] = matrix[matr_index];
//
//             // between two nonconsecutive native entries (one or more inherited
//             // entries may be inbetween) we add pointers to
//             // predecessor/successor
//             // if (native_pt != -1 && curr_index != native_pt + 1) {
//             //    augmented_matrix[curr_index].set_predecessor(
//             //        *augmented_matrix[native_pt].entry);
//             //    augmented_matrix[native_pt].set_successor(
//             //        *augmented_matrix[curr_index].entry);
//             // }
//             native_pt = curr_index;
//
//             // add pointers from the previous inherited entries to
//             // the new native entry
//             // if (inherited_pt != -1){
//             //    for (size_t ii = curr_index - 1; ii >= pointers[index + 1];
//             //         --ii) {
//             //       if (augmented_matrix[ii].entry_t == utils::NATIVE)
//             //          break;
//             //       else
//             //          augmented_matrix[ii].set_successor(
//             //              *augmented_matrix[curr_index].entry);
//             //    }
//             // }
//
//             // add bridge to next row/column
//             // bridge is from augmented_matrix[curr_index] to
//             // augmented_matrix[bridge_pt]
//             while (bridge_pt < augmented_pointers[index] &&
//                    matrix_index(*ds.matrix[bridge_pt].entry, matrix.layout) <
//                        matrix_index(*ds.matrix[curr_index].entry, matrix.layout))
//                bridge_pt++;
//
//             if (matrix_index(*ds.matrix[bridge_pt].entry, matrix.layout) <=
//                 matrix_index(*ds.matrix[curr_index].entry, matrix.layout))
//                ds.matrix[curr_index].bridge = &ds.matrix[bridge_pt];
//
// #ifdef DEBUG_VERBOSE
//             printf("(%ld,%ld,%d) native added (%ld,%ld,%d) at %ld\n",
//                    matrix[matr_index].i, matrix[matr_index].j, matrix[matr_index].a,
//                    ds.matrix[curr_index].entry->i, ds.matrix[curr_index].entry->j,
//                    ds.matrix[curr_index].entry->a, curr_index);
//             if (ds.matrix[curr_index].bridge != nullptr)
//                printf("(%ld,%ld,%d) bridge (%ld,%ld,%d)\n",
//                       ds.matrix[curr_index].entry->i,
//                       ds.matrix[curr_index].entry->j,
//                       ds.matrix[curr_index].entry->a,
//                       ds.matrix[curr_index].bridge->entry->i,
//                       ds.matrix[curr_index].bridge->entry->j,
//                       ds.matrix[curr_index].bridge->entry->a);
// #endif
//             // if matrix and augmented_matrix have an entry at i
//             // increase counter for augmented_matrix
//             // obj.entry->i can segfault if obj is nullptr
//             if (ds.matrix[prev_index].entry != nullptr &&
//                 matrix_index(*ds.matrix[prev_index].entry, matrix.layout) ==
//                     matrix_index(matrix[matr_index], matrix.layout)) {
//                even = !even;
//                prev_index++;
//             }
//             curr_index++;
//             matr_index++;
//
//             // otherwise the entry is inherited
//          } else if (matrix_index(*ds.matrix[prev_index].entry, matrix.layout) ==
//                         curr_position &&
//                     prev_index < augmented_pointers[index]) {
//             if (even) {
//                ds.matrix[curr_index] = ds.matrix[prev_index];
//                // decreasing index means copy to previous row/column
//                matrix.layout
//                    ? ds.matrix[curr_index].entry->j =
//                          matrix_index(*ds.matrix[prev_index].entry, utils::COLUMN) -
//                          1
//                    : ds.matrix[curr_index].entry->i =
//                          matrix_index(*ds.matrix[prev_index].entry, utils::ROW) - 1;
//
//                // if (native_pt != -1)
//                //    augmented_matrix[curr_index].set_predecessor(
//                //        *augmented_matrix[native_pt].entry);
//
//                inherited_pt = curr_index;
//
//                // add bridge to next row/column
//                // bridge is from augmented_matrix[curr_index] to
//                // augmented_matrix[bridge_pt]
//                ds.matrix[curr_index].bridge = &ds.matrix[prev_index];
//
// #ifdef DEBUG_VERBOSE
//                printf("(%ld,%ld,%d) inherited added (%ld,%ld,%d) at %ld\n",
//                       ds.matrix[prev_index].entry->i,
//                       ds.matrix[prev_index].entry->j,
//                       ds.matrix[prev_index].entry->a,
//                       ds.matrix[curr_index].entry->i,
//                       ds.matrix[curr_index].entry->j,
//                       ds.matrix[curr_index].entry->a, curr_index);
//                if (ds.matrix[curr_index].bridge != nullptr)
//                   printf("(%ld,%ld,%d) bridge (%ld,%ld,%d)\n",
//                          ds.matrix[curr_index].entry->i,
//                          ds.matrix[curr_index].entry->j,
//                          ds.matrix[curr_index].entry->a,
//                          ds.matrix[curr_index].bridge->entry->i,
//                          ds.matrix[curr_index].bridge->entry->j,
//                          ds.matrix[curr_index].bridge->entry->a);
// #endif
//                curr_index++;
//             }
//             // prev_index and even are updated regardless of assignment
//             prev_index++;
//             even = !even;
//          }
// #ifdef DEBUG_VERBOSE
//          utils::print_sparse_as_dense(ds.matrix, n);
// #endif
//       }
//       // pointer to beginning of current row/column in augmented_matrix
//       augmented_pointers[index + 1] = curr_index;
//    }
//
//    // compute the bridges between index_vector and
//    // first row/column of augmented_matrix
//    // for (ssize_t index = n - 1; index >= 0; --index) {
//    size_t index = n;
//    bridge_pt = augmented_pointers[n] - 1;
//    while (index > 0 && bridge_pt >= augmented_pointers[n - 1]) {
//       if (matrix_index(*ds.vector[bridge_pt].entry, matrix.layout) > index - 1)
//          ds.index_vector[index - 1]->bridge = &ds.vector[--bridge_pt];
//       else
//          ds.index_vector[index - 1]->bridge = &ds.vector[bridge_pt];
//       index--;
//
// #ifdef DEBUG_VERBOSE
//       if (ds.index_vector[index]->bridge == nullptr)
//          printf("%ld -> nullptr\n", index);
//       else
//          printf("%ld -> (%ld,%ld,%d)\n", index,
//          ds.index_vector[index]->bridge->entry->i,
//                 ds.index_vector[index]->bridge->entry->j,
//                 ds.index_vector[index]->bridge->entry->a);
// #endif
//    }
//
// #ifdef DEBUG_VERBOSE
//    // ei2 of type cascading::augmented_entry<utils::ee>
//    auto ei2 = ds.index_vector[4];
//    for (size_t index = 0; index < n; ++index) {
//       if (ei2->entry != nullptr)
//          printf("jumping from (%ld,%ld,%d) ", ei2->entry->i, ei2->entry->j,
//                 ei2->entry->a);
//       else
//          printf("jumping from null    ");
//       if (ei2->bridge == nullptr)
//          printf("bridge is null\n");
//       ei2 = ei2->bridge;
//       printf("assignment done\n");
//       if (ei2->entry == nullptr)
//          printf("yes\n");
//       if (ei2->entry != nullptr)
//          printf("to (%ld,%ld,%d) ", ei2->entry->i, ei2->entry->j, ei2->entry->a);
//       else
//          printf("to null    ");
//       printf("\n");
//       // if (ei2.predecessor == nullptr) {
//       //    printf("(np) \n");
//       // } else
//       //    sketch[index] = ei2.predecessor->a;
//    }
// #endif
//
//    exit(0);
//
//    return ds;
// }

/****************************************************************************************/

template <typename T>
inline matrix::entry_t extract_sketchh(const T &index_vector, matrix::index_t i1,
                                       matrix::index_t i2, matrix::index_t n) {

   matrix::entry_t *sketch = new matrix::entry_t[n]();
   // ei2 of type cascading::augmented_entry<utils::ee>
   auto ei2 = index_vector[i2];

   printf("Extracting sketch...\n");

   printf("i2 = %ld\n", i2);

   for (matrix::index_t index = 0; index < n; ++index) {
      if (ei2.entry != nullptr)
         printf("jumping from (%ld,%ld,%d) ", ei2.entry->i, ei2.entry->j, ei2.entry->a);
      else
         printf("jumping from null    ");
      ei2 = *ei2.bridge;
      printf("to (%ld,%ld,%d) ", ei2.entry->i, ei2.entry->j, ei2.entry->a);
      if (ei2.entry_type == utils::INHERITED) {
         printf("inherited ");
         if (ei2.predecessor == nullptr) {
            printf("(np) \n");
            sketch[index] = 0;
         } else
            sketch[index] = ei2.predecessor->a;
      } else
         sketch[index] = ei2.entry->a;
   }
   printf("****\n");

   printf("i1 = %ld\n", i1);

   i1 = std::max((i1 - 1), 0); // easier with i1-1
   if (i1 == 0)
      return *sketch;
   auto ei1 = index_vector[i1];

   for (matrix::index_t index = 0; index < n; ++index) {
      ei1 = *ei1.bridge;
      if (ei1.entry_t == utils::INHERITED)
         sketch[index] -= ei1.predecessor->a;
      else
         sketch[index] = ei1.entry->a;
   }
   for (matrix::index_t index = 0; index < n; ++index)
      printf("%d ", sketch[index]);
   printf("\n");
   // int c = getchar();
   return *sketch;
}

} // namespace cascading