#include <vector>

namespace coalesced {

template <typename T> class bucket {

 public:
   bucket(){};
   ~bucket(){};

   // bucket main entries
   std::vector<T> v;
   // auxiliary bucket for inherited entries
   std::vector<T> w;

   bucket(size_t n) { v.resize(n); };
};

/****************************************************************************************/

template <typename T>
inline bucket<utils::ee> **range_coalesced(const T &matrix, size_t h, size_t n,
                                           size_t order) {

   size_t ds_size = std::floor(h / n);
   size_t *counters = new size_t[n + 1]();
   size_t *splitters = new size_t[n]();
   auto *coalesced_buckets = new coalesced::bucket<utils::ee>[ds_size];
   auto **index_vector = new coalesced::bucket<utils::ee> *[n];

   // sort_by returns either row or column index depending on order
   utils::sort_by<utils::ee> index_type;

   // counting vector where counters[i] is the number
   // of nnz in row/column i
   for (size_t index = 0; index < h; ++index)
      counters[index_type(matrix[index], order) + 1]++;

   // prefix sum over the counting vector
   for (size_t index = 0; index < n; ++index)
      counters[index + 1] += counters[index];

   // compute ds_size splitters where between splitters[i]
   // and splitters[i+1] there are at most 2n entries
   size_t bucket_pt = 1;
   for (size_t index = 0; index < n; ++index) {
      if (counters[index] >= n * bucket_pt)
         ++bucket_pt;
      splitters[index] = bucket_pt - 1;
      index_vector[index] = &coalesced_buckets[bucket_pt - 1];
   }

   // partition the entries in buckets based on splitters
   // TO DO: the following is NOT cache efficient!
   for (size_t index = 0; index < h; ++index)
      coalesced_buckets[splitters[index_type(matrix[index], order)]]
          .v.push_back(matrix[index]);

   // updating dummy entries so they point to the max
   // prefix sum less than the current row/column
   auto it_zero = std::begin(coalesced_buckets[0].v);
   for (size_t index = 0; index < n; ++index) {
      while ((it_zero + 1) != coalesced_buckets[0].v.end() &&
             index_type(*(it_zero + 1), 1 - order) <= index)
         ++it_zero;
      if (index_type(*it_zero, 1 - order) != index)
         coalesced_buckets[0].w.push_back(*new utils::ee(index, index, 0));
   }

   // update inherited entry vector of current bucket
   // with max prefix sum of previous bucket
   for (size_t bk = 0; bk < ds_size - 1; ++bk) {
      auto it_next = coalesced_buckets[bk + 1].v.begin();
      auto it_curr_v = coalesced_buckets[bk].v.begin();
      auto it_curr_w = coalesced_buckets[bk].w.begin();
      for (size_t index = 0; index < n; ++index) {
         while ((it_next + 1) != coalesced_buckets[bk + 1].v.end() &&
                index_type(*(it_next + 1), 1 - order) <= index)
            ++it_next;
         while ((it_curr_v + 1) != coalesced_buckets[bk].v.end() &&
                index_type(*(it_curr_v + 1), 1 - order) <= index)
            ++it_curr_v;

         if (index_type(*it_next, 1 - order) == index)
            ;
         else if (index_type(*it_curr_v, 1 - order) == index)
            coalesced_buckets[bk + 1].w.push_back(
                *new utils::ee(it_curr_v->j, it_curr_v->j, it_curr_v->a));
         else if (index_type(*it_curr_w, 1 - order) == index)
            coalesced_buckets[bk + 1].w.push_back(
                *new utils::ee(it_curr_w->j, it_curr_w->j, it_curr_w->a));
         else
            coalesced_buckets[bk + 1].w.push_back(
                *new utils::ee(index, index, 0));

         if (it_curr_w + 1 != coalesced_buckets[bk].w.end())
            ++it_curr_w;
      }
   }

#ifdef DEBUG_VERBOSE
   for (auto it = std::begin(coalesced_buckets[0].v);
        it != std::end(coalesced_buckets[0].v); ++it)
      printf("(%ld,%ld,%d) ", it->i, it->j, it->a);

   for (size_t bk = 0; bk < ds_size; ++bk) {
      printf("\n v : ");
      for (std::vector<utils::ee>::iterator it =
               std::begin(coalesced_buckets[bk].v);
           it != std::end(coalesced_buckets[bk].v); ++it)
         printf("(%ld,%ld,%d) ", it->i, it->j, it->a);
      printf("\n w : ");
      for (std::vector<utils::ee>::iterator it =
               std::begin(coalesced_buckets[bk].w);
           it != std::end(coalesced_buckets[bk].w); ++it)
         printf("(%ld,%ld,%d) ", it->i, it->j, it->a);
      printf("\n***\n");
   }
#endif

   return index_vector;
}

/****************************************************************************************/
// template <typename T>
// inline bucket<utils::ee> **aug(const T &matrix, size_t h, size_t n,
//                                size_t order) {}
/****************************************************************************************/

template <typename T>
inline utils::e_type *extract_sketch(const T &index_vector, size_t i1,
                                     size_t i2, size_t n, size_t order) {

   // sort_by returns either row or column index depending on order
   utils::sort_by<utils::ee> index_type;

   assert(i2 >= i1);
   utils::e_type *sketch = new utils::e_type[n]();
   auto it_v2 = index_vector[i2]->v.begin();
   auto it_w2 = index_vector[i2]->w.begin();

   // for (auto it = std::begin(index_vector[i2]->v);
   //      it != std::end(index_vector[i2]->v); ++it)
   //    printf("(%ld,%ld,%d) ", it->i, it->j, it->a);

   for (size_t index = 0; index < n; ++index) {
      while (index_type(*it_v2, 1 - order) < index &&
             index_type(*it_v2, order) <= i2)
         ++it_v2;

      if (index_type(*it_v2, 1 - order) == index &&
          index_type(*it_v2, order) <= i2) {
         sketch[index] = it_v2->a;
      } else {
         sketch[index] = it_w2->a;
         ++it_w2;
      }
   }
   // for (size_t index = 0; index < n; ++index)
   //    printf("%d ", sketch[index]);
   // printf("\n");

   i1 = std::max((ssize_t)(i1 - 1), (ssize_t)0); // easier with i1-1
   if (i1 == 0)
      return sketch;
   auto it_v1 = index_vector[i1]->v.begin();
   auto it_w1 = index_vector[i1]->w.begin();

   for (size_t index = 0; index < n; ++index) {
      while (index_type(*it_v1, 1 - order) < index &&
             index_type(*it_v1, order) <= i1)
         ++it_v1;

      if (index_type(*it_v1, 1 - order) == index &&
          index_type(*it_v1, order) <= i1) {
         sketch[index] -= it_v1->a;
      } else {
         sketch[index] -= it_w1->a;
         ++it_w1;
      }
   }

   return sketch;
}

} // namespace coalesced