#pragma once

#include <Eigen/Core>
#include <Eigen/SparseCore>

namespace pemu::linalg {

using Scalar = double;
using Index = Eigen::Index;

using SparseMatrix = Eigen::SparseMatrix<Scalar, Eigen::ColMajor, Index>;

using Vector = Eigen::Vector<Scalar, Eigen::Dynamic>;

using VectorRef = Eigen::Ref<Vector>;

using ConstVectorRef = Eigen::Ref<const Vector>;

struct LinearSystem {
  SparseMatrix A;
  Vector b;
};

}  // namespace pemu::linalg