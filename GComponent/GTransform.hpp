#pragma once

#include <eigen3/Eigen/Dense>
#include <tuple>

namespace GComponent {
using Eigen::Matrix;
using Eigen::Vector3d;
using Eigen::Vector4d;
using Eigen::Matrix3d;
using Eigen::Matrix4d;
using std::pair;
using vec3d = Vector3d;
using vec4d = Vector4d;
using twistd = Eigen::Matrix<double, 6, 1>;
using SE3d = Eigen::Matrix4d;
using AdMatrixd = Eigen::Matrix<double, 6, 6>;
using DynMatrixd = Matrix<double, Eigen::Dynamic, Eigen::Dynamic>;

template<typename _Scaler>
class so3 :public Matrix<_Scaler, 3, 1>
{

};

template<typename _Scaler>
class SO3 : public Matrix<_Scaler, 3, 3>
{

};

template<typename _Scaler>
class se3 : public Matrix<_Scaler, 6, 1>
{

};

template<typename _Scaler>
class SE3 : public Matrix<_Scaler, 4, 4>
{

};

template<typename _Scaler>
class AdjointSE3 : public Matrix<_Scaler, 6, 6>
{

};

// Rt + p
inline vec3d affineProduct(const Matrix4d & lmat, const vec3d & point)
{
    return lmat.block(0, 0, 3, 3)* point + lmat.block(0, 3, 3, 1);
}

inline pair<vec3d, double> GetRotateAxisAngle(const vec3d& _w)
{
	return std::make_pair(_w.normalized(), _w.norm());
}

inline Matrix3d CrossMatrix(const vec3d& _x)
{
		
	Matrix3d mat;
	mat.setZero();
	mat(0, 1) = -_x[2];
	mat(0, 2) = _x[1];
	mat(1, 0) = _x[2];
	mat(1, 2) = -_x[0];
	mat(2, 0) = -_x[1];
	mat(2, 1) = _x[0];
	return mat;
}

inline vec3d CrossVec3(const Matrix3d& _skew_sym_mat)
{
    return vec3d(
		_skew_sym_mat(2, 1), 
		_skew_sym_mat(0, 2), 
		_skew_sym_mat(1, 0));
}

inline Matrix3d Roderigues(const Matrix3d& _cross_m, double _theta)
{
	return Matrix3d(
		Matrix3d::Identity() + 
		sin(_theta) * _cross_m + 
		(1 - cos(_theta)) * (_cross_m * _cross_m));
}

inline Matrix3d Roderigues(const vec3d& _axis, double _theta)
{
	Matrix3d _cross_m = CrossMatrix(_axis.normalized());
	return Roderigues(_cross_m, _theta);
}

inline Matrix3d Roderigues(const vec3d& _v)
{
	double _theta = _v.norm();
	Matrix3d _cross_m = CrossMatrix(_v.normalized());
	return Roderigues(_cross_m, _theta);
}

inline Matrix3d DiffRoderigues(const Matrix3d& _cross_m, double _theta)
{
    return Matrix3d(
        Matrix3d::Identity() +
        cos(_theta) * _cross_m +
        (1 + sin(_theta)) * (_cross_m * _cross_m));
}

inline vec3d LogMapSO3Toso3(const Matrix3d& mat)
{
	double theta = acos((mat.trace() - 1.) * .5);
		
	if (theta < 1e-5)
        return vec3d::Zero();
		
	if ((EIGEN_PI - theta) < 1e-5)
		return 
		EIGEN_PI 
		* (1 / sqrt(2 * (1 + mat(2, 2)))) 
        * vec3d(mat(0, 2), mat(1, 2), 1 + mat(2, 2));
		
    vec3d vec = CrossVec3((mat - mat.transpose()) / (2. * sin(theta)));
	return (theta * vec);
}


inline SE3d InverseSE3(const SE3d& mat)
{
	SE3d _inv_mat = SE3d::Identity();
	
	_inv_mat.block(0, 0, 3, 3) = mat.block(0, 0, 3, 3).transpose();
	_inv_mat.block(0, 3, 3, 1) = 
		-_inv_mat.block(0, 0, 3, 3) * mat.block(0, 3, 3, 1);

	return _inv_mat;
}

inline pair<twistd, double> GetTwistAxisAngle(const twistd& _t)
{
    const vec3d& w = _t.block(0, 0, 3, 1);
    const vec3d& v = _t.block(3, 0, 3, 1);

	double _theta = w.norm();
	if (_theta < 1e-5)
		_theta = v.norm();

	return std::make_pair(_t / _theta, _theta);
}

inline SE3d ExpMapping(const twistd& _norm_t, double _theta)
{
    const vec3d& w = _norm_t.block(0, 0, 3, 1);
    const vec3d& v = _norm_t.block(3, 0, 3, 1);

	SE3d mat = SE3d::Identity();
	Matrix3d _cross_axis = CrossMatrix(w.normalized());
	mat.block(0, 0, 3, 3) = Roderigues(_cross_axis, _theta);

	Matrix3d G = _theta * Matrix3d::Identity()
		+ (1 - cos(_theta)) * _cross_axis
		+ (_theta - sin(_theta)) * _cross_axis * _cross_axis;
	mat.block(0, 3, 3, 1) = G * v;

	return mat;
}

inline SE3d ExpMapping(const twistd & _t)
{
	auto && [S, _theta] = GetTwistAxisAngle(_t);
	return ExpMapping(S, _theta);
}

inline SE3d DiffExpMapping(const twistd& _norm_t, double _theta)
{
    const vec3d& w = _norm_t.block(0, 0, 3, 1);
    const vec3d& v = _norm_t.block(3, 0, 3, 1);

    SE3d mat = SE3d::Identity();
    Matrix3d _cross_axis = CrossMatrix(w.normalized());
    mat.block(0, 0, 3, 3) = DiffRoderigues(_cross_axis, _theta);

    Matrix3d G =
        Matrix3d::Identity()
        + (1 + sin(_theta)) * _cross_axis
        + (1 - cos(_theta)) * _cross_axis * _cross_axis;
    mat.block(0, 3, 3, 1) = G * v;

    return mat;
}

inline twistd ScrewToTwist(const vec3d & _q, const vec3d & _s, double h)
{
	twistd _normal_twist = twistd::Zero();
	_normal_twist.block(0, 0, 3, 1) = _s.normalized();
	_normal_twist.block(3, 0, 3, 1) = -_s.cross(_q) + h * _s.normalized();
	return _normal_twist;
}

inline twistd LogMapSE3Tose3(const SE3d& _T)
{
	const Matrix3d& R = _T.block(0, 0, 3, 3);
    const vec3d& p = _T.block(0, 3, 3, 1);
	
    vec3d w = LogMapSO3Toso3(R);
	auto&& [axis, theta] = GetRotateAxisAngle(w);
	Matrix3d _cross_axis = CrossMatrix(axis);
	
	Matrix3d G_Inv = 1 / theta * Matrix3d::Identity()
		- 0.5 * _cross_axis
		+ (1 / theta - 0.5 / tan(0.5 * theta)) * _cross_axis * _cross_axis;
	if (abs(theta) < 1e-5)
		G_Inv = Matrix3d::Identity();
	else
		G_Inv *= theta;
    vec3d v = G_Inv * p;

	twistd t = twistd::Zero();
	t.block(0, 0, 3, 1) = w;
	t.block(3, 0, 3, 1) = v;

	return t;
}

inline AdMatrixd Adjoint(const SE3d & T)
{
	const Matrix3d& R = T.block(0, 0, 3, 3);
	const Matrix3d _cross_p = CrossMatrix(T.block(0, 3, 3, 1));

	AdMatrixd adMat = AdMatrixd::Zero();
	adMat.block(0, 0, 3, 3) =
	adMat.block(3, 3, 3, 3) = R;
	adMat.block(3, 0, 3, 3) = _cross_p * R;
		
	return adMat;
}

}
