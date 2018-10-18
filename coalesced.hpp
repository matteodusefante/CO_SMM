#include <vector>

namespace coalesced {

template <typename T> struct bucket {
   std::vector<T> native; // for native entries
   std::vector<T> inhert; // for inherited entries
};

template <typename T> class coalesced_buckets {

 public:
   std::vector<bucket<T> *> index_vector;
   std::vector<bucket<T>> buckets;

   coalesced_buckets(matrix::index_t n, ssize_t ds_size)
       : index_vector(n), buckets(ds_size) {}
   ~coalesced_buckets() {
      for (auto &el : index_vector)
         delete el;
   }
   // coalesced_buckets(const coalesced_buckets &other) // copy constructor
   //     : coalesced_buckets(other.cstring) {}
   // coalesced_buckets(coalesced_buckets &&other) noexcept // move constructor
   //     : cstring(std::exchange(other.cstring, nullptr)) {}
   // coalesced_buckets &operator=(const coalesced_buckets &other) { // copy assignment
   //    return *this = coalesced_buckets(other);
   // }
   // coalesced_buckets &operator=(coalesced_buckets &&other) noexcept { // move
   // assignment
   //    std::swap(cstring, other.cstring);
   //    return *this;
   // }
};

/****************************************************************************************/

template <typename T>
inline coalesced_buckets<matrix::coo_entry> range_coalesced(const T &matrix, ssize_t h,
                                                            matrix::index_t n) {

   ssize_t ds_size = std::floor(h / n);
   std::vector<ssize_t> counters(n + 1, 0);
   std::vector<ssize_t> splitters(n, 0);
   coalesced_buckets<matrix::coo_entry> ds(n, ds_size);

   // return_index returns either row or column index depending on layout
   utils::return_index<matrix::coo_entry> matrix_index;

   // counting vector where counters[i] is the number
   // of nnz in row/column i
   for (auto &vec : matrix)
      for (auto &entry : vec)
         counters[matrix_index(entry, matrix.layout) + 1]++;

   // prefix sum over the counting vector
   std::partial_sum(std::begin(counters), std::end(counters), std::begin(counters));

   // compute ds_size splitters where between splitters[i]
   // and splitters[i+1] there are at most 2n entries
   ssize_t bucket_pt = 1;
   for (matrix::index_t index = 0; index < n; ++index) {
      if (counters[index] >= n * bucket_pt)
         ++bucket_pt;
      splitters[index] = bucket_pt - 1;
      ds.index_vector[index] = &ds.buckets[bucket_pt - 1]; // pointer to bucket
   }

   // partition the entries in buckets based on splitters
   // TO DO: the following is NOT cache efficient!
   for (auto &vec : matrix)
      for (auto &entry : vec)
         ds.buckets[splitters[matrix_index(entry, matrix.layout)]].native.push_back(
             entry);

   // updating dummy entries so they point to the max
   // prefix sum less than the current row/column
   auto it_zero = std::begin(ds.buckets[0].native);
   for (matrix::index_t index = 0; index < n; ++index) {
      // looking for the entry with index "index"
      while ((it_zero + 1) != std::end(ds.buckets[0].native) &&
             matrix_index(*(it_zero + 1), 1 - matrix.layout) <= index)
         ++it_zero;
      if (matrix_index(*it_zero, 1 - matrix.layout) != index) {
         if (matrix.layout == utils::COLUMN_MAJOR)
            ds.buckets[0].inhert.emplace_back(
                matrix::entry<matrix::entry_t>(0, index, 0));
         else
            ds.buckets[0].inhert.emplace_back(
                matrix::entry<matrix::entry_t>(index, 0, 0));
      }
   }

   // update inherited entry vector of current bucket
   // with max prefix sum of previous bucket
   for (matrix::index_t bk = 0; bk < ds_size - 1; ++bk) {
      auto it_prev = std::begin(ds.buckets[bk].native);
      auto it_curr = std::begin(ds.buckets[bk + 1].native);
      auto it_curr_aux = std::begin(ds.buckets[bk].inhert);

      for (matrix::index_t index = 0; index < n; ++index) {
         while ((it_curr + 1) != std::end(ds.buckets[bk + 1].native) &&
                matrix_index(*(it_curr + 1), 1 - matrix.layout) <= index)
            ++it_curr;
         while ((it_prev + 1) != std::end(ds.buckets[bk].native) &&
                matrix_index(*(it_prev + 1), 1 - matrix.layout) <= index)
            ++it_prev;

         // if (matrix_index(*it_curr, 1 - matrix.layout) == index);
         if (matrix_index(*it_prev, 1 - matrix.layout) == index)
            ds.buckets[bk + 1].inhert.emplace_back(
                matrix::entry<matrix::entry_t>(it_prev->j, it_prev->j, it_prev->a));
         else if (matrix_index(*it_prev, 1 - matrix.layout) == index)
            ds.buckets[bk + 1].inhert.emplace_back(matrix::entry<matrix::entry_t>(
                it_curr_aux->j, it_curr_aux->j, it_curr_aux->a));
         else
            ds.buckets[bk + 1].inhert.emplace_back(
                matrix::entry<matrix::entry_t>(index, index, 0));

         if (it_curr_aux + 1 != std::end(ds.buckets[bk].inhert))
            ++it_curr_aux;
      }
   }

#ifdef DEBUG_VERBOSE
   for (auto it = std::begin(ds.buckets[0].native); it != std::end(ds.buckets[0].native);
        ++it)
      std::cout << "(" << it->i << "," << it->j << "," << it->a << ")";

   for (ssize_t bk = 0; bk < ds_size; ++bk) {
      std::cout << std::endl << "native : ";
      for (auto it = std::begin(ds.buckets[bk].native);
           it != std::end(ds.buckets[bk].native); ++it)
         std::cout << "(" << it->i << "," << it->j << "," << it->a << ")";
      std::cout << std::endl << "inhert : ";
      for (auto it = std::begin(ds.buckets[bk].inhert);
           it != std::end(ds.buckets[bk].inhert); ++it)
         std::cout << "(" << it->i << "," << it->j << "," << it->a << ")";
   }
#endif

   return ds;
}

/****************************************************************************************/

template <typename T>
inline std::vector<matrix::entry_t>
extract_sketch(const T &ds, matrix::index_t i1, matrix::index_t i2,
               matrix::index_t n, ssize_t layout) {

   // return_index returns either row or column index depending on layout
   utils::return_index<matrix::coo_entry> matrix_index;

   assert(i2 >= i1);
   std::vector<matrix::entry_t> sketch(n, 0);

   for (auto &el : ds.index_vector[i2]->native)
      if (matrix_index(el, 1 - layout) <= i2)
         sketch[matrix_index(el, layout)] += el.a;

   for (auto &el : ds.index_vector[i2]->inhert)
      if (matrix_index(el, 1 - layout) <= i2)
         sketch[matrix_index(el, layout)] += el.a;

   // for (auto &i : sketch)
   //    printf("%d ", i);
   // printf("\n");

   i1 = std::max((i1 - 1), 0); // easier with i1-1
   assert(0 <= i1 && i1 <= n);

   if (i1 == 0)
      return sketch;

   for (auto &el : ds.index_vector[i1]->native)
      if (matrix_index(el, 1 - layout) <= i1)
         sketch[matrix_index(el, layout)] += el.a;

   for (auto &el : ds.index_vector[i1]->inhert)
      if (matrix_index(el, 1 - layout) <= i1)
         sketch[matrix_index(el, layout)] += el.a;

   return sketch;
}

} // namespace coalesced