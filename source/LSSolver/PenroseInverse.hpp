#pragma once
#include <Eigen/dense>

using Eigen::Matrix;
template<class _Scalar>
class PenroseInverseSolver{
	using _tDM = Matrix<_Scalar, Eigen::Dynamic, Eigen::Dynamic>;
public:
	_tDM operator()(const _tDM& m) {
		if (m.rows() <= m.cols())
			return m.transpose() * ((m * m.transpose()).inverse());
		return (m.transpose() * m).inverse() * m.transpose();
	}
};