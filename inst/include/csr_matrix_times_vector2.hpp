template <typename T2__, typename T5__>
inline
Eigen::Matrix<typename boost::math::tools::promote_args<T2__, T5__>::type,
              Eigen::Dynamic, 1>
csr_matrix_times_vector2(const int& m,
                         const int& n,
                         const Eigen::Matrix<T2__, Eigen::Dynamic, 1>& w,
                         const std::vector<int>& v,
                         const std::vector<int>& u,
                         const Eigen::Matrix<T5__, Eigen::Dynamic, 1>& b,
                         std::ostream* pstream__) {
  Eigen::Map<const Eigen::SparseMatrix<T2__,Eigen::RowMajor> >
    sm(m, n, w.size(), &u[0], &v[0], &w[0]);
  return sm * b;
}

