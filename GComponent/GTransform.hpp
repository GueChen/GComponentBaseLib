#pragma once

#include <eigen3/Eigen/Dense>
#include <tuple>

namespace GComponent {
using Eigen::Matrix;
using Eigen::Vector;
using Eigen::Vector3d;
using Eigen::Vector4d;
using Eigen::Matrix3d;
using Eigen::Matrix4d;
using std::pair;
using DynMatrixd = Matrix<double, Eigen::Dynamic, Eigen::Dynamic>;
template<class _Scaler>
using SE3      = Matrix<_Scaler, 4, 4>;
template<class _Scaler>
using twist    = Vector<_Scaler, 6>;
template<class _Scaler>
using AdMatrix = Matrix<_Scaler, 6, 6>;
using SE3d   = SE3<double>;
using twistd = twist<double>;

//template<typename _Scaler>
//class so3 :public Matrix<_Scaler, 3, 1>
//{
//
//};
//
//template<typename _Scaler>
//class SO3 : public Matrix<_Scaler, 3, 3>
//{
//
//};
//
//template<typename _Scaler>
//class se3 : public Matrix<_Scaler, 6, 1>
//{
//
//};
//
//template<typename _Scaler>
//class SE3 : public Matrix<_Scaler, 4, 4>
//{
//
//};
//
//template<typename _Scaler>
//class AdjointSE3 : public Matrix<_Scaler, 6, 6>
//{
//
//};

// Rt + p
template<class Scaler>
inline Vector<Scaler, 3> affineProduct(const Matrix<Scaler, 4, 4> & lmat, const Vector<Scaler, 3> & point)
{
    return static_cast<Vector<Scaler, 3>>(lmat.block(0, 0, 3, 3)* point + lmat.block(0, 3, 3, 1));
}
template<class Scaler>
inline pair<Vector<Scaler, 3>, Scaler> GetRotateAxisAngle(const Vector<Scaler, 3>& _w)
{
	return std::make_pair(_w.normalized(), _w.norm());
}

template<class Scaler>
Matrix<Scaler, 3, 3> CrossMatrix(const Vector<Scaler,3>& _x)
{
		
	Matrix<Scaler, 3, 3> mat;
	mat.setZero();
	mat(0, 1) = -_x[2];
	mat(0, 2) = _x[1];
	mat(1, 0) = _x[2];
	mat(1, 2) = -_x[0];
	mat(2, 0) = -_x[1];
	mat(2, 1) = _x[0];
	return mat;
}

template<class Scaler>
inline Vector<Scaler, 3> CrossVec3(const Matrix<Scaler, 3, 3>& _skew_sym_mat)
{
    return Vector<Scaler, 3>(
		_skew_sym_mat(2, 1), 
		_skew_sym_mat(0, 2), 
		_skew_sym_mat(1, 0));
}

template<class Scaler>
inline Matrix<Scaler, 3, 3> Roderigues(const Matrix<Scaler, 3, 3>& _cross_m, Scaler _theta)
{
	return Matrix<Scaler, 3, 3>(
		Matrix<Scaler, 3, 3>::Identity() +
		sin(_theta) * _cross_m + 
		(1 - cos(_theta)) * (_cross_m * _cross_m));
}

template<class Scaler>
inline Matrix<Scaler, 3, 3> Roderigues(const Vector<Scaler, 3>& _axis, Scaler _theta)
{
	Matrix<Scaler, 3, 3> 
		_cross_m = CrossMatrix(_axis.normalized());
	return Roderigues(_cross_m, _theta);
}

template<class Scaler>
inline Matrix<Scaler, 3, 3> Roderigues(const Vector<Scaler, 3>& _v)
{
	Scaler _theta = _v.norm();
	Matrix<Scaler, 3, 3> _cross_m = CrossMatrix(_v.normalized());
	return Roderigues(_cross_m, _theta);
}

template<class _Scaler>
inline Matrix<_Scaler, 3, 3> DiffRoderigues(const Matrix<_Scaler, 3, 3>& _cross_m, _Scaler _theta)
{
    return  Matrix<_Scaler, 3, 3>(
		Matrix<_Scaler, 3, 3>::Identity() +
        cos(_theta) * _cross_m +
        (1 + sin(_theta)) * (_cross_m * _cross_m));
}

template<class Scaler>
Vector<Scaler, 3> LogMapSO3Toso3(const Matrix<Scaler, 3, 3>& mat)
{
	double theta = acos((mat.trace() - 1.) * .5);
		
	if (theta < 1e-5)
        return Vector<Scaler, 3>::Zero();
		
	if ((EIGEN_PI - theta) < 1e-5)
		return 
		EIGEN_PI 
		* (1 / sqrt(2 * (1 + mat(2, 2)))) 
        * Vector<Scaler, 3>(mat(0, 2), mat(1, 2), 1 + mat(2, 2));
		
	Vector<Scaler, 3> vec = CrossVec3(
		static_cast<Matrix<Scaler,3, 3>>((mat - mat.transpose()) / (2. * sin(theta))));
	return (theta * vec);
}

template<class _Scaler>
inline SE3<_Scaler> InverseSE3(const SE3<_Scaler>& mat)
{
	SE3<_Scaler> _inv_mat = SE3<_Scaler>::Identity();
	
	_inv_mat.block(0, 0, 3, 3) = mat.block(0, 0, 3, 3).transpose();
	_inv_mat.block(0, 3, 3, 1) = -_inv_mat.block(0, 0, 3, 3) * mat.block(0, 3, 3, 1);
	return _inv_mat;
}

template<class _Scaler>
pair<twist<_Scaler>, _Scaler> GetTwistAxisAngle(const twist<_Scaler>& _t)
{
    const Vector<_Scaler, 3>& w = _t.block(0, 0, 3, 1);
    const Vector<_Scaler, 3>& v = _t.block(3, 0, 3, 1);

	_Scaler _theta = w.norm();
	if (_theta < 1e-5) _theta = v.norm();

	return std::make_pair(_theta < 1e-5?_t:_t / _theta, _theta);
}

template<class _Scaler>
SE3<_Scaler> ExpMapping(const twist<_Scaler>& _norm_t, _Scaler _theta)
{
    const Vector<_Scaler, 3>& w = _norm_t.block(0, 0, 3, 1);
    const Vector<_Scaler, 3>& v = _norm_t.block(3, 0, 3, 1);

	SE3<_Scaler> mat		= SE3<_Scaler>::Identity();
	Matrix<_Scaler, 3, 3> _cross_axis	
							= CrossMatrix(w.normalized());
	mat.block(0, 0, 3, 3)   = Roderigues(_cross_axis, _theta);

	Matrix<_Scaler, 3, 3> 
		G = 
		_theta * Matrix<_Scaler, 3, 3>::Identity()
		+ (1 - cos(_theta)) * _cross_axis
		+ (_theta - sin(_theta)) * _cross_axis * _cross_axis;
	mat.block(0, 3, 3, 1) = G * v;

	return mat;
}

template<class _Scaler>
inline SE3<_Scaler> ExpMapping(const twist<_Scaler> & _t)
{
	auto && [S, _theta] = GetTwistAxisAngle(_t);
	return ExpMapping(S, _theta);
}

template<class _Scaler>
SE3<_Scaler> DiffExpMapping(const twist<_Scaler>& _norm_t, _Scaler _theta)
{
    const Vector<_Scaler, 3>& w = _norm_t.block(0, 0, 3, 1);
    const Vector<_Scaler, 3>& v = _norm_t.block(3, 0, 3, 1);

    SE3<_Scaler> mat = SE3<_Scaler>::Identity();
    Matrix<_Scaler, 3, 3> _cross_axis = CrossMatrix(w.normalized());
    mat.block(0, 0, 3, 3) = DiffRoderigues(_cross_axis, _theta);

	Matrix<_Scaler, 3, 3> 
		G =
		Matrix<_Scaler, 3, 3>::Identity()
        + (1 + sin(_theta)) * _cross_axis
        + (1 - cos(_theta)) * _cross_axis * _cross_axis;
    mat.block(0, 3, 3, 1) = G * v;

    return mat;
}

template<class _Scaler>
inline twist<_Scaler> ScrewToTwist(const Vector<_Scaler, 3> & _q, const Vector<_Scaler, 3>& _s, _Scaler h)
{
	twist<_Scaler> _normal_twist	= twist<_Scaler>::Zero();
	_normal_twist.block(0, 0, 3, 1) = _s.normalized();
	_normal_twist.block(3, 0, 3, 1) = -_s.cross(_q) + h * _s.normalized();
	return _normal_twist;
}

template<class _Scaler>
twist<_Scaler> LogMapSE3Tose3(const SE3<_Scaler>& _T)
{
	const Matrix<_Scaler, 3, 3>& R = _T.block(0, 0, 3, 3);
    const Vector<_Scaler, 3>   & p = _T.block(0, 3, 3, 1);
	
	Vector<_Scaler, 3> w = LogMapSO3Toso3(R);
	auto&& [axis, theta] = GetRotateAxisAngle(w);
	Matrix<_Scaler, 3, 3> _cross_axis = CrossMatrix(axis);
	
	Matrix<_Scaler, 3, 3>
		G_Inv = 
		1. / theta * Matrix<_Scaler, 3, 3>::Identity()
		- 0.5 * _cross_axis
		+ (1. / theta - 0.5 / tan(0.5 * theta)) * _cross_axis * _cross_axis;
	if (abs(theta) < 1e-5)
		G_Inv = Matrix<_Scaler, 3, 3> ::Identity();
	else
		G_Inv *= theta;
	Vector<_Scaler, 3> v = G_Inv * p;

	twist<_Scaler> t = twist<_Scaler>::Zero();
	t.block(0, 0, 3, 1) = w;
	t.block(3, 0, 3, 1) = v;

	return t;
}

template<class _Scaler>
AdMatrix<_Scaler> Adjoint(const SE3<_Scaler> & T)
{
	const Matrix<_Scaler, 3, 3>& R = T.block(0, 0, 3, 3);
	const Matrix<_Scaler, 3, 3> _cross_p = CrossMatrix(
		static_cast<Vector<_Scaler, 3>>(T.block(0, 3, 3, 1)));

	AdMatrix<_Scaler> adMat = AdMatrix<_Scaler>::Zero();
	adMat.block(0, 0, 3, 3) =
	adMat.block(3, 3, 3, 3) = R;
	adMat.block(3, 0, 3, 3) = _cross_p * R;
		
	return adMat;
}

}
