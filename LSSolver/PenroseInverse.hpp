#pragma once
#include <Eigen/dense>

using Eigen::Matrix;
class PenroseInverseSolver{
	using _tDM = Matrix<double, Eigen::Dynamic, Eigen::Dynamic>;
public:
	_tDM operator()(const _tDM& m) {
		if (m.rows() <= m.cols())
			return m.transpose() * ((m * m.transpose()).inverse());
		return (m.transpose() * m).inverse() * m.transpose();
	}
};