// #include "coalesced.hpp"

/****************************************************************************************/

template <typename T>
std::vector<coalesced::bucket<utils::coo_entry> *> coalesced::range_coalesced(const T &matrix,
                                                                   size_t h, size_t n) {

   size_t ds_size = std::floor(h / n);
   std::vector<size_t> counters(n + 1, 0);
   std::vector<size_t> splitters(n, 0);
   std::vector<coalesced::bucket<utils::coo_entry>> coalesced_buckets(ds_size);
   std::vector<coalesced::bucket<utils::coo_entry> *> index_vector(n);

   // return_index returns either row or column index depending on layout
   utils::return_index<utils::coo_entry> matrix_index;

   // counting vector where counters[i] is the number
   // of nnz in row/column i
   for (auto &vec : matrix)
      for (auto &entry : vec)
         counters[matrix_index(entry, matrix.layout) + 1]++;

   // prefix sum over the counting vector
   for (auto it = std::begin(counters); it != std::end(counters); ++it)
      *(it + 1) += *it;

   // compute ds_size splitters where between splitters[i]
   // and splitters[i+1] there are at most 2n entries
   size_t bucket_pt = 1;
   for (utils::index_type index = 0; index < n; ++index) {
      if (counters[index] >= n * bucket_pt)
         ++bucket_pt;
      splitters[index] = bucket_pt - 1;
      index_vector[index] = &coalesced_buckets[bucket_pt - 1]; // pointer to bucket
   }

   // partition the entries in buckets based on splitters
   // TO DO: the following is NOT cache efficient!
   for (auto &vec : matrix)
      for (auto &entry : vec)
         coalesced_buckets[splitters[matrix_index(entry, matrix.layout)]].v.push_back(
             entry);

   // updating dummy entries so they point to the max
   // prefix sum less than the current row/column
   auto it_zero = std::begin(coalesced_buckets[0].v);
   for (utils::index_type index = 0; index < n; ++index) {
      // looking for the entry with index "index"
      while ((it_zero + 1) != std::end(coalesced_buckets[0].v) &&
             matrix_index(*(it_zero + 1), 1 - matrix.layout) <= index)
         ++it_zero;
      if (matrix_index(*it_zero, 1 - matrix.layout) != index) {
         if (matrix.layout == utils::COLUMN_MAJOR)
            coalesced_buckets[0].w.emplace_back(
                matrix::entry<utils::entry_type>(0, index, 0));
         else
            coalesced_buckets[0].w.emplace_back(
                matrix::entry<utils::entry_type>(index, 0, 0));
      }
   }

   // update inherited entry vector of current bucket
   // with max prefix sum of previous bucket
   for (utils::index_type bk = 0; bk < ds_size - 1; ++bk) {
      auto it_prev = std::begin(coalesced_buckets[bk].v);
      auto it_curr = std::begin(coalesced_buckets[bk + 1].v);
      auto it_curr_aux = std::begin(coalesced_buckets[bk].w);
      for (utils::index_type index = 0; index < n; ++index) {
         while ((it_curr + 1) != std::end(coalesced_buckets[bk + 1].v) &&
                matrix_index(*(it_curr + 1), 1 - matrix.layout) <= index)
            ++it_curr;
         while ((it_prev + 1) != std::end(coalesced_buckets[bk].v) &&
                matrix_index(*(it_prev + 1), 1 - matrix.layout) <= index)
            ++it_prev;

         // if (matrix_index(*it_curr, 1 - matrix.layout) == index);
         if (matrix_index(*it_prev, 1 - matrix.layout) == index)
            coalesced_buckets[bk + 1].w.emplace_back(
                matrix::entry<utils::entry_type>(it_prev->j, it_prev->j, it_prev->a));
         else if (matrix_index(*it_prev, 1 - matrix.layout) == index)
            coalesced_buckets[bk + 1].w.emplace_back(matrix::entry<utils::entry_type>(
                it_curr_aux->j, it_curr_aux->j, it_curr_aux->a));
         else
            coalesced_buckets[bk + 1].w.emplace_back(
                matrix::entry<utils::entry_type>(index, index, 0));

         if (it_curr_aux + 1 != std::end(coalesced_buckets[bk].w))
            ++it_curr_aux;
      }
   }

#ifdef DEBUG_VERBOSE
   for (auto it = std::begin(coalesced_buckets[0].v);
        it != std::end(coalesced_buckets[0].v); ++it)
      std::cout << "(" << it->i << "," << it->j << "," << it->a << ")";

   for (size_t bk = 0; bk < ds_size; ++bk) {
      std::cout << std::endl << "v : ";
      for (auto it = std::begin(coalesced_buckets[bk].v);
           it != std::end(coalesced_buckets[bk].v); ++it)
         std::cout << "(" << it->i << "," << it->j << "," << it->a << ")";
      std::cout << std::endl << "w : ";
      for (auto it = std::begin(coalesced_buckets[bk].w);
           it != std::end(coalesced_buckets[bk].w); ++it)
         std::cout << "(" << it->i << "," << it->j << "," << it->a << ")";
   }
#endif

   return index_vector;
}

/****************************************************************************************/
// template <typename T>
// inline bucket<utils::ee> **aug(const T &matrix, size_t h, size_t n,
//                                size_t layout) {}
/****************************************************************************************/

template <typename T>
std::vector<utils::entry_type>
coalesced::extract_sketch(const T &index_vector, utils::index_type i1,
                          utils::index_type i2, size_t n, size_t layout) {

   // return_index returns either row or column index depending on layout
   utils::return_index<utils::coo_entry> matrix_index;

   assert(i2 >= i1);
   std::vector<utils::entry_type> sketch(n, 0);
   auto it_v_i2 = std::begin(index_vector[i2]->v);
   auto it_w_i2 = std::begin(index_vector[i2]->w);

   for (auto &it : index_vector[i2]->v)
      std::cout << "(" << it.i << "," << it.j << "," << it.a << ")";
   printf("\n");

   for (auto &it : index_vector[i2]->w)
      std::cout << "(" << it.i << "," << it.j << "," << it.a << ")";
   printf("\n");

   // printf("here\n");
   //
   // for (auto it = (index_vector[i2])->begin(); it != (index_vector[i2])->end(); ++it)
   // printf("(%ld,%ld,%d) ", it->i, it->j, it->a);
   // printf("%d\n", *it_2[0]);
   //
   for (utils::index_type index = 0; index < n; ++index) {
      while (matrix_index(*it_v_i2, 1 - layout) == index &&
             matrix_index(*it_v_i2, layout) <= i2)
         sketch[index] += it_v_i2++->a;

      while (matrix_index(*it_w_i2, 1 - layout) == index)
         sketch[index] += it_w_i2++->a;

      for (auto &i : sketch)
         printf("%d ", i);
      printf("\n");
      // while (matrix_index(*it_v_i2, 1 - layout) < index &&
      //        matrix_index(*it_v_i2, layout) <= i2 &&
      //        it_v_i2 != std::end(index_vector[i2]->v))
      //    ++it_v_i2;
      //
      // if (matrix_index(*it_v_i2, 1 - layout) == index &&
      //     matrix_index(*it_v_i2, layout) <= i2)
      //    sketch[index] += it_v_i2->a;
      // else
      //    sketch[index] += it_w_i2->a;
      // ++it_w_i2;
   }
   printf("here\n");
   //
   // printf("here\n");
   //
   // // for (size_t index = 0; index < n; ++index)
   // //    printf("%d ", sketch[index]);
   // // printf("\n");
   //
   i1 = std::max((ssize_t)(i1 - 1), (ssize_t)0); // easier with i1-1
   if (i1 == 0)
      return sketch;
   auto it_v_i1 = std::begin(index_vector[i1]->v);
   auto it_w_i1 = std::begin(index_vector[i1]->w);

   for (utils::index_type index = 0; index < n; ++index) {
      while (matrix_index(*it_v_i1, 1 - layout) < index &&
             matrix_index(*it_v_i1, layout) <= i1)
         ++it_v_i1;

      if (matrix_index(*it_v_i1, 1 - layout) == index &&
          matrix_index(*it_v_i1, layout) <= i1)
         sketch[index] -= it_v_i1->a;
      else
         sketch[index] -= it_w_i1->a;
      ++it_w_i1;
   }

   printf("here\n");

   return sketch;
}
