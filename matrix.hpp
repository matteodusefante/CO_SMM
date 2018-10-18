#include <numeric>

namespace matrix {

typedef int32_t entry_t;
typedef int32_t index_t;

template <typename T> class entry {

 public:
   index_t i, j;
   T a;

   entry() : i(-1), j(-1), a(0) {} // dummy entry as default
   ~entry() {}
   entry(index_t i, index_t j, const T &a) : i(i), j(j), a(a) {}

   entry<T> &operator=(const entry<T> &rhs) {
      this->i = rhs.i;
      this->j = rhs.j;
      this->a = rhs.a;
      return *this; // return the result by reference
   }

   // void set_i(size_t i) { this->i = i; }
   // void set_j(size_t j) { this->j = j; }
   // void set_a(const T &a) { this->a = a; }
};

typedef entry<entry_t> coo_entry;

/****************************************************************************************/

// mxn dense matrix
template <typename T> class dense : public std::vector<std::vector<T>> {
   typedef std::vector<std::vector<T>> matrix_type;
   index_t m, n;

 public:
   dense() {}
   dense(index_t n) : matrix_type(n, std::vector<T>(n, 0)) {}
   dense(index_t m, index_t n) : matrix_type(m, std::vector<T>(n, 0)) {}
   ~dense() {}
};

/****************************************************************************************/

// m sparse column/row vectors
template <typename T> class sparse : public std::vector<std::vector<T>> {
   typedef std::vector<std::vector<T>> matrix_type;
   index_t n;

 public:
   ssize_t layout;

   sparse() {}
   sparse(index_t m) : matrix_type(m, std::vector<T>(0)) {}
   sparse(index_t m, ssize_t layout)
       : matrix_type(m, std::vector<T>(0)), layout(layout) {}
   ~sparse() {}
};

/****************************************************************************************/

template <typename T> inline entry_t inner_product(const T &vec1, const T &vec2) {

   index_t temp = 0;
   auto it_vec1 = std::begin(vec1);
   auto it_vec2 = std::begin(vec2);

   for (; it_vec1 != std::end(vec1) || it_vec2 != std::end(vec2); ++it_vec1, ++it_vec2)
      temp += (*it_vec1) * (*it_vec2);

#ifdef DEBUG_VERBOSE
   std::cout << "= " << temp << std::endl;
#endif

   return temp;
}

} // namespace matrix