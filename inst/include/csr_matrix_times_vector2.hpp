// Sparse (CSR, zero-based v and u) matrix times dense vector. Since Stan 2.26, data vectors arrive as Eigen::Map
// and parameter vectors as var matrices, so both arguments are templated on any Eigen vector type.
template <typename TW, typename TB>
inline Eigen::Matrix<stan::return_type_t<TW, TB>, Eigen::Dynamic, 1>
csr_matrix_times_vector2(const int& m,
                         const int& n,
                         const TW& w,
                         const std::vector<int>& v,
                         const std::vector<int>& u,
                         const TB& b,
                         std::ostream* pstream__) {
  using T2 = stan::value_type_t<TW>;
  const Eigen::Matrix<T2, Eigen::Dynamic, 1> w_ = w;
  const Eigen::Matrix<stan::value_type_t<TB>, Eigen::Dynamic, 1> b_ = b;
  Eigen::Map<const Eigen::SparseMatrix<T2, Eigen::RowMajor> >
    sm(m, n, w_.size(), &u[0], &v[0], w_.data());
  return sm * b_;
}

