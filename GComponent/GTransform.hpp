#ifndef _GTRANSFORM_HPP
#define _GTRANSFORM_HPP

#include <eigen3/Eigen/Dense>
#include <tuple>

namespace GComponent {
using Eigen::Matrix;
using Eigen::Vector;
using std::pair;
using DynMatrixd = Matrix<double, Eigen::Dynamic, Eigen::Dynamic>;
template<class _Scaler>
using SE3	 = Eigen::Matrix4<_Scaler>;
template<class _Scaler>
using SO3	 = Eigen::Matrix3<_Scaler>;
template<class _Scaler>
using Twist  = Vector<_Scaler, 6>;
template<class _Scaler>
using AdMatrix = Matrix<_Scaler, 6, 6>;
using SE3d   = SE3<double>;
using Twistd = Twist<double>;
using std::tuple;


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

template<class _Scaler>
Matrix<_Scaler, 3, 3> Shear(const Vector<_Scaler, 3>& shear) 
{
	Matrix<_Scaler, 3, 3> shear_mat = Matrix<_Scaler, 3, 3>::Identity();
	shear_mat(0, 1) = shear.x();
	shear_mat(0, 2) = shear.y();
	shear_mat(1, 2) = shear.z();
	return shear_mat;
}

template<class _Scaler>
Matrix<_Scaler, 3, 3> ShearInverse(const Vector<_Scaler, 3>& shear) 
{
	Matrix<_Scaler, 3, 3> shear_mat = Matrix<_Scaler, 3, 3>::Identity();
	shear_mat(0, 1) = -shear.x();
	shear_mat(0, 2) = shear.x() * shear.z() - shear.y();
	shear_mat(1, 2) = -shear.z();
	return shear_mat;
}

template<class _Scaler>
Matrix<_Scaler, 3, 3> Scale(const Vector<_Scaler, 3>& scale) 
{
	Matrix<_Scaler, 3, 3> scale_mat = Matrix<_Scaler, 3, 3>::Zero();
	scale_mat(0, 0) = scale.x();
	scale_mat(1, 1) = scale.y();
	scale_mat(2, 2) = scale.z();
	return scale_mat;
}

template<class _Scaler>
Matrix<_Scaler, 3, 3> ScaleInverse(const Vector<_Scaler, 3>& scale) 
{
	Matrix<_Scaler, 3, 3> scale_mat = Matrix<_Scaler, 3, 3>::Zero();
	if (abs(scale.x()) > 1e-8) scale_mat(0, 0) = 1.0 / scale.x();	
	if (abs(scale.y()) > 1e-8) scale_mat(1, 1) = 1.0 / scale.y();
	if (abs(scale.z()) > 1e-8) scale_mat(2, 2) = 1.0 / scale.z();
	return scale_mat;
}

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
Matrix<Scaler, 3, 3> CrossMatrix(const Vector<Scaler,3>& x)
{
		
	Matrix<Scaler, 3, 3> mat;
	mat.setZero();
	mat(0, 1) = -x[2];
	mat(0, 2) = x[1];
	mat(1, 0) = x[2];
	mat(1, 2) = -x[0];
	mat(2, 0) = -x[1];
	mat(2, 1) = x[0];
	return mat;
}

template<class Scaler>
inline Vector<Scaler, 3> CrossVec3(const Matrix<Scaler, 3, 3>& skew_sym_mat)
{
    return Vector<Scaler, 3>(
		skew_sym_mat(2, 1), 
		skew_sym_mat(0, 2), 
		skew_sym_mat(1, 0));
}

template<class Scaler>
inline Matrix<Scaler, 3, 3> Roderigues(const Matrix<Scaler, 3, 3>& cross_m, Scaler theta)
{
	return Matrix<Scaler, 3, 3>(
		Matrix<Scaler, 3, 3>::Identity() +
		sin(theta) * cross_m + 
		(1 - cos(theta)) * (cross_m * cross_m));
}

template<class Scaler>
inline Matrix<Scaler, 3, 3> Roderigues(const Vector<Scaler, 3>& axis_, Scaler _theta)
{
	return Roderigues(CrossMatrix(axis_.normalized()), _theta);
}

template<class Scaler>
inline Matrix<Scaler, 3, 3> Roderigues(const Vector<Scaler, 3>& v)
{
	return Roderigues(CrossMatrix(v.normalized()), v.norm());
}

template<class _Scaler>
inline Matrix<_Scaler, 3, 3> DiffRoderigues(const Matrix<_Scaler, 3, 3>& _cross_m, _Scaler _theta)
{
    return  Matrix<_Scaler, 3, 3>(
		Matrix<_Scaler, 3, 3>::Identity() +
        cos(_theta) * _cross_m +
        (1 + sin(_theta)) * (_cross_m * _cross_m));
}

template<class _Scaler>
Vector<_Scaler, 3> LogMapSO3Toso3(const Matrix<_Scaler, 3, 3>& mat)
{	
	_Scaler cos_val = (mat.trace() - 1.) * 0.5;
	if (cos_val > 1.0)  cos_val = 1.0;
	if (cos_val < -1.0) cos_val = -1.0;

	_Scaler theta = acos(cos_val);
		
	if (theta < 1e-5)
        return Vector<_Scaler, 3>::Zero();
		
	if ((EIGEN_PI - theta) < 1e-5)
		return 
		EIGEN_PI 
		* (1 / sqrt(2 * (1 + mat(2, 2)))) 
        * Vector<_Scaler, 3>(mat(0, 2), mat(1, 2), 1 + mat(2, 2));
		
	Vector<_Scaler, 3> vec = CrossVec3(
		static_cast<Matrix<_Scaler,3, 3>>((mat - mat.transpose()) / (2. * sin(theta))));
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
pair<Twist<_Scaler>, _Scaler> GetTwistAxisAngle(const Twist<_Scaler>& _t)
{
    const Vector<_Scaler, 3>& w = _t.block(0, 0, 3, 1);
    const Vector<_Scaler, 3>& v = _t.block(3, 0, 3, 1);

	_Scaler _theta = w.norm();
	if (_theta < 1e-5) _theta = v.norm();

	return std::make_pair(_theta < 1e-5?_t:_t / _theta, _theta);
}

template<class _Scaler>
SE3<_Scaler> ExpMapping(const Twist<_Scaler>& _norm_t, _Scaler _theta)
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
inline SE3<_Scaler> ExpMapping(const Twist<_Scaler> & _t)
{
	auto && [S, _theta] = GetTwistAxisAngle(_t);
	return ExpMapping(S, _theta);
}

template<class _Scaler>
SE3<_Scaler> DiffExpMapping(const Twist<_Scaler>& _norm_t, _Scaler _theta)
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
inline Twist<_Scaler> ScrewToTwist(const Vector<_Scaler, 3> & q, const Vector<_Scaler, 3>& w, _Scaler h = 0.0)
{
	Twist<_Scaler> _normal_twist	= Twist<_Scaler>::Zero();
	_normal_twist.block(0, 0, 3, 1) = w.normalized();
	_normal_twist.block(3, 0, 3, 1) = q.cross(w) + h * w.normalized();
	return _normal_twist;
}

template<class _Scaler>
Twist<_Scaler> LogMapSE3Tose3(const SE3<_Scaler>& _T)
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

	Twist<_Scaler> t = Twist<_Scaler>::Zero();
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

/// <summary>
/// Decomposing a 3 x 3 matrix into two 3 x 3 matrices of orthogonal matrix Q and upper triangular matrix R
/// <para>将一个可逆 3 x 3 实矩阵分解为 Q R 两个矩阵，其中 Q 为单位正交阵， R 为上三角矩阵 </para>
/// </summary>
/// <typeparam name="_Scaler"></typeparam>
/// <param name="r_mat"></param>
/// <returns></returns>
template<class _Scaler>
pair<SO3<_Scaler>, Matrix<_Scaler, 3, 3>> QRDecompositionMat3(const Matrix<_Scaler, 3, 3>& mat)
{
	auto GetInv = [](_Scaler val)->_Scaler {
		return abs(val) > std::numeric_limits<_Scaler>::epsilon() ? 1.0 / val : 0.0;
	};
	Vector<_Scaler, 3>
		m_cow_1 = mat.block(0, 0, 3, 1),
		m_cow_2 = mat.block(0, 1, 3, 1),
		m_cow_3 = mat.block(0, 2, 3, 1);
	
	// Step1: Schmidt Orthogonalization get Q Matrix
	Vector<_Scaler, 3> Q_cow_1, Q_cow_2, Q_cow_3;
	_Scaler inv_length = 0.0;

	inv_length = GetInv(m_cow_1.norm());
	Q_cow_1 = m_cow_1 * inv_length;

	Q_cow_2 = m_cow_2 - Q_cow_1.dot(m_cow_2) * Q_cow_1;
	inv_length = GetInv(Q_cow_2.norm());
	Q_cow_2 = Q_cow_2 * inv_length;

	Q_cow_3 = m_cow_3 - Q_cow_1.dot(m_cow_3) * Q_cow_1 - Q_cow_2.dot(m_cow_3) * Q_cow_2;
	inv_length = GetInv(Q_cow_3.norm());
	Q_cow_3 = Q_cow_3 * inv_length;

	SO3<_Scaler>		Q;
	Q.block(0, 0, 3, 1) = Q_cow_1;
	Q.block(0, 1, 3, 1) = Q_cow_2;
	Q.block(0, 2, 3, 1) = Q_cow_3;

	if (Q.determinant() < 0.0) Q = -Q;

	// Step2: Get R use R = Q^t * M
	Matrix<_Scaler, 3, 3> R = Q.transpose() * mat;

	return { Q, R };
}

/// <summary>
/// Decoposing a 3 x 3 upper triangular Matrix into two 3 x 1 vectors of scale and shear 
/// <para>分解 3 x 3 上三角实矩阵为 缩放 和 剪切 两个 3 x 1 向量 </para>
/// </summary>
/// <typeparam name="_Scaler">{float/double} precision</typeparam>
/// <param name="ut_mat"> 　
/// <para>　　{const ref Matrix3} </para>
/// <para>　　upper triangular Matrix 上三角实矩阵</para>
/// </param>
/// <returns></returns>
template<class _Scaler>
inline pair<Vector<_Scaler, 3>, Vector<_Scaler, 3>> SSDecompositionMat3(const Matrix<_Scaler, 3, 3>& ut_mat) 
{	
	return  {Vector<_Scaler, 3>(ut_mat(0, 0),					// scale
								ut_mat(1, 1),
								ut_mat(2, 2)),
			 Vector<_Scaler, 3>(ut_mat(0, 1) / ut_mat(0, 0),	// shear
								ut_mat(0, 2) / ut_mat(0, 0),
								ut_mat(1, 2) / ut_mat(1, 1))};	
}

template<class _Scaler>
tuple<SO3<_Scaler>, Vector<_Scaler, 3>, Vector<_Scaler, 3>> RSSDecompositionMat3(const Matrix<_Scaler, 3, 3>& mat)
{
	// Step1: get The QR Decomposition
	auto [Q, R] = QRDecompositionMat3(mat);
	
	// Step2: Decompose R into scale and shear
	auto [scale, shear] = SSDecompositionMat3(R);

	//     Q ∈ SO(3)
	return { Q, scale, shear};
}

template<class _Scaler>
tuple<Vector<_Scaler, 3>, Vector<_Scaler, 3>, Vector<_Scaler, 3>, Vector<_Scaler, 3>> TRSSDecompositionMat4(const Matrix<_Scaler, 4, 4>& mat) {
	auto [rot_mat, scale, shear] = RSSDecompositionMat3(static_cast<Matrix<_Scaler, 3, 3>>(mat.block(0, 0, 3, 3)));
	return { static_cast<Vector<_Scaler, 3>>(mat.block(0, 3, 3, 1)),
			 LogMapSO3Toso3(rot_mat),
			 scale,
			 shear};			
}

template<class _Scaler>
pair<Matrix<_Scaler, 3, 3>, Vector<_Scaler, 3>> RtDecompositionMat4(const Matrix<_Scaler, 4, 4>& mat)
{
	Matrix<_Scaler, 3, 3> rot_mat	= mat.block(0, 0, 3, 3);
	Vector<_Scaler, 3>	  trans_vec = mat.block(0, 3, 3, 1);
	return { rot_mat, trans_vec };
}

template<class _Scaler>
pair<Vector<_Scaler, 3>, Vector<_Scaler, 3>> rtDecompositionMat4(const Matrix<_Scaler, 4, 4>& mat)
{
	Vector<_Scaler, 3>	  
		trans_vec = mat.block(0, 3, 3, 1),
		rot_vec   = LogMapSO3Toso3(static_cast<Matrix<_Scaler, 3, 3>>(mat.block(0, 0, 3, 3)));
	return {rot_vec, trans_vec};
}

}

#endif
