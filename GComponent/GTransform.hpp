/**
 *  @file  	GTransform.hpp
 *  @brief 	some transform functions use Lie group and decompose theory.
 *  @author Gue Chen<guechen@buaa.edu.cn>
 *  @date 	Nov  xxx, 2021 create this file
*			July 4th, 2022 using concept to reconstruct functions
 **/
#ifndef _GTRANSFORM_HPP
#define _GTRANSFORM_HPP

#include <Concept/gconcept.hpp>

#include <eigen3/Eigen/Dense>
#include <tuple>

namespace GComponent {
using Eigen::Matrix;
using Eigen::Vector;
using std::pair;
using std::tuple;
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

/// <summary>
/// make 3 x 3 shear matrix form 3 x 1 shear vector
/// <para>
/// 使用 3 x 1 剪切向量获取 3 x 3 剪切矩阵
/// </para>
/// </summary>
/// <param name="shear">	cref	{vec3}	3 x 1 shear vctor	3 x 1 剪切向量		</param>
/// <returns>				val		{mat3}	3 x 3 shear matrix	3 x 3 剪切矩阵		</returns>
template<Vec3Convertible Derived>
auto Shear(const Eigen::MatrixBase<Derived>& shear)
{
	using Mat3 = Eigen::Matrix3<typename Derived::Scalar>;
	Mat3 shear_mat = Mat3::Identity();
	shear_mat(0, 1) = shear.x();
	shear_mat(0, 2) = shear.y();
	shear_mat(1, 2) = shear.z();
	return shear_mat;
}

/// <summary>
/// make 3 x 3 inverse matrix of shear from 3 x 1 shear vector
/// <para>
/// 使用 3 x 1 剪切向量获取 3 x 3 剪切逆矩阵
/// </para>
/// </summary>
/// <param name="shear">	cref	{vec3}	3 x 1 shear vctor			3 x 1 剪切向量			</param>
/// <returns>				val		{mat3}	3 x 3 inverse shear matrix	3 x 3 剪切逆矩阵		</returns>
template<Vec3Convertible Derived>
auto ShearInverse(const Eigen::MatrixBase<Derived>& shear) 
{
	using Mat3 = Eigen::Matrix3<typename Derived::Scalar>;
	Mat3 shear_mat  = Mat3::Identity();
	shear_mat(0, 1) = -shear.x();
	shear_mat(0, 2) = shear.x() * shear.z() - shear.y();
	shear_mat(1, 2) = -shear.z();
	return shear_mat;
}

/// <summary>
/// make 3 x 3 scale matrix form 3 x 1 scale vector
/// <para>
/// 使用 3 x 1 缩放向量获取 3 x 3 缩放矩阵
/// </para>
/// </summary>
/// <param name="scale">	cref	{vec3}	3 x 1 scale vctor	3 x 1 缩放向量		</param>
/// <returns>				val		{mat3}	3 x 3 scale matrix	3 x 3 缩放矩阵		</returns>
template<Vec3Convertible Derived>
auto Scale(const Eigen::MatrixBase<Derived>& scale)
{
	using Mat3 = Eigen::Matrix3<typename Derived::Scalar>;
	Mat3 scale_mat = Mat3::Zero();
	scale_mat(0, 0) = scale.x();
	scale_mat(1, 1) = scale.y();
	scale_mat(2, 2) = scale.z();
	return scale_mat;
}

/// <summary>
/// make 3 x 3 inverse scale matrix form 3 x 1 scale vector
/// <para>
/// 使用 3 x 1 缩放向量获取 3 x 3 缩放逆矩阵
/// </para>
/// </summary>
/// <param name="scale">	cref	{vec3}	3 x 1 scale vctor			3 x 1 缩放向量			</param>
/// <returns>				val		{mat3}	3 x 3 inverse scale matrix	3 x 3 缩放逆矩阵		</returns>
template<Vec3Convertible Derived>
auto ScaleInverse(const Eigen::MatrixBase<Derived>& scale)
{
	using Mat3 = Eigen::Matrix3<typename Derived::Scalar>;
	Mat3 scale_mat = Mat3::Zero();
	if (abs(scale.x()) > 1e-8) scale_mat(0, 0) = 1.0 / scale.x();	
	if (abs(scale.y()) > 1e-8) scale_mat(1, 1) = 1.0 / scale.y();
	if (abs(scale.z()) > 1e-8) scale_mat(2, 2) = 1.0 / scale.z();
	return scale_mat;
}

/// <summary>
/// Caculator homogeneous Mat with a R^3 vector return as a vector, represent translate and rotate on an R^3 vector as Rt + p
/// <para>
/// 对一个 3 x 1 向量施加齐次矩阵变换，几何直观上代表着对一个点施加旋转与平移变换
/// </para>
/// </summary>
/// <param name="mat">	cref	{SE3}	homogeneous matrix in SE3	齐次变换矩阵					</param>
/// <param name="vec">	cref	{vec3}	vector in R^3				3 x 1 向量						</param>
/// <returns>			val		{vec3}	new vector after transform	变换后的新向量					</returns>
template<Mat4Convertible DerivedMat, Vec3Convertible DerivedVec> 
requires MatScalarEquivalence<DerivedMat, DerivedVec>
inline Eigen::Vector3<typename DerivedMat::Scalar> 
AffineProduct(const Eigen::MatrixBase<DerivedMat>& mat, const Eigen::MatrixBase<DerivedVec>& vec)
{
    return mat.block(0, 0, 3, 3) * vec + mat.block(0, 3, 3, 1);
}

/// <summary>
/// Decompose angle axis vector to normlized axis direction with a real number angle in [rad.]
/// <para>
/// 将轴角表示向量解耦为正则化的轴向量与以实数表示的弧度角
/// </para>
/// </summary>
/// <param name="vec">	cref	{vec3}		[in]	an 3 x 1 Vector							3 x 1 向量			</param>
/// <returns>			val		{vec3, R}	[out]	norm angle dir and real number angle	正则化轴向与实数角	</returns>
template<Vec3Convertible Derived>
inline 
pair<Eigen::Vector3<typename Derived::Scalar>, typename Derived::Scalar> 
GetRotateAxisAngle(const Eigen::MatrixBase<Derived>& vec)
{
	return std::make_pair(vec.normalized(), vec.norm());
}

/// <summary>
/// Hat operator on so3 to convert element to 3 x 3 skew-symmatrix from 3 x 1 vector
/// <para>
/// Hat 算子将 so3 从 3 x 1 向量表述转换到 3 x 3 反对称矩阵表述
/// </para>
/// </summary>
/// <param name="vec">	cref	{so3}	[in]	an so3 vector form			3 x 1 轴角表述			</param>
/// <returns>			val		{[so3]}	[out]	an skew-symmatrix matrix	3 x 3 反对称矩阵		</returns>
template<Vec3Convertible Derived>
auto Hat(const Eigen::MatrixBase<Derived>& vec)
{
	using Mat3 = Eigen::Matrix3<typename Derived::Scalar>;
	Mat3 mat  = Mat3::Zero();
	mat(0, 1) = -vec[2]; 
	mat(0, 2) = vec[1];
	mat(1, 0) = vec[2];
	mat(1, 2) = -vec[0];
	mat(2, 0) = -vec[1];
	mat(2, 1) = vec[0];
	return mat;
}

/// <summary>
/// Vee operator on so3 to convert element to 3 x 1 vector from 3 x 3 skew-symmatrix 
/// <para>
/// Vee 算子将 so3 从 3 x 3 反对称矩阵表述转换到 3 x 1 向量表述
/// </para>
/// </summary>
/// <param name="skew_sym_mat">	cref	{[so3]}	[in]	an skew-symmatrix matrix	3 x 3 反对称矩阵	</param>
/// <returns>					val		{so3}	[out]	so3 element vector form		3 x 1 轴角表述		</returns>
template<Mat3Convertible Derived>
inline Eigen::Vector3<typename Derived::Scalar> 
Vee(const Eigen::MatrixBase<Derived>& skew_sym_mat)
{
    return Eigen::Vector3<typename Derived::Scalar>(skew_sym_mat(2, 1), skew_sym_mat(0, 2), skew_sym_mat(1, 0));
}

/// <summary>
/// Using Roderigues formula to get SO3 element from [so3] axis and real number angle in [rad.]
/// <para>
/// 使用 Roderigues 公式从反对称矩阵转轴和实数弧度转角获取对应旋转矩阵
/// </para>
/// </summary>
/// <param name="skew_sym_m">	cref	{[so3]}	[in]	an skew-symmatrix matrix	3x3 反对称矩阵	</param>
/// <param name="theta">		val		{R}		[in]	a real number angle [rad.]	实数弧度角		</param>
/// <returns>					val		{SO3}	[out]	rotation matrix				旋转矩阵		</returns>
template<Mat3Convertible Derived, ScalarEquivalence<Derived> U>
inline Eigen::Matrix3<typename Derived::Scalar> 
Roderigues(const Eigen::MatrixBase<Derived>& skew_sym_m, U theta)
{
	using Mat3 = Eigen::Matrix3<typename Derived::Scalar>;
	return Mat3::Identity() + sin(theta) * skew_sym_m + (1 - cos(theta)) * (skew_sym_m * skew_sym_m);
}

/// <summary>
/// Using Roderigues formula to get SO3 element from so3 axis and real number angle in [rad.]
/// <para>
/// 使用 Roderigues 公式从转轴向量和实数弧度转角获取对应旋转矩阵
/// </para>
/// </summary>
/// <param name="axis">			cref	{so3}	[in]	an 3 x 1 vector				3 x 1 转轴向量	</param>
/// <param name="theta">		val		{R}		[in]	a real number angle [rad.]	实数弧度角		</param>
/// <returns>					val		{SO3}	[out]	rotation matrix				旋转矩阵		</returns>
template<Vec3Convertible Derived, ScalarEquivalence<Derived> U>
inline auto Roderigues(const Eigen::MatrixBase<Derived>& axis, U theta)
{
	return Roderigues(Hat(axis.normalized()), theta);
}

/// <summary>
/// Using Roderigues formula to get SO3 matrix from so3 axis angle vector with unit [rad.]
/// <para>
/// 使用 Roderigues 公式从轴角向量获取对应旋转矩阵
/// </para>
/// </summary>
/// <param name="v">	cref	{so3}	[in]	3 x 1 axis angle vecotr in [rad.]	3 x 1 轴角向量	</param>
/// <returns>			val		{SO3}	[out]	3 x 3 rotate matrix					3 x 3 旋转矩阵	</returns>
template<Vec3Convertible Derived>
inline auto Roderigues(const Eigen::MatrixBase<Derived>& v)
{
	return Roderigues(Hat(v.normalized()), v.norm());
}

template<Mat3Convertible Derived, ScalarEquivalence<Derived> U>
inline Eigen::Matrix3<typename Derived::Scalar> 
DiffRoderigues(const Eigen::MatrixBase<Derived>& skew_sym_m, U theta)
{
	using Mat3 = Eigen::Matrix3<typename Derived::Scalar>;
    return  Mat3::Identity() + cos(theta) * skew_sym_m + (1 + sin(theta)) * (skew_sym_m * skew_sym_m);
}

template<Mat3Convertible Derived>
Eigen::Vector3<typename Derived::Scalar> 
LogMapSO3Toso3(const Eigen::MatrixBase<Derived>& mat)
{	
	assert(mat.rows() == 3 && mat.cols() == 3 && "logrithimic: SO3 -> so3 , size not matching");
	using Scalar = typename Derived::Scalar;
	using Vec3	 = Eigen::Vector3<Scalar>;

	// Clamp cos_val back to [-1.0, 1.0] in avoiding numerical error
	Scalar cos_val = std::max(
					 std::min(0.5 * (mat.trace() - 1.), 1.0), -1.0);
	Scalar theta  = acos(cos_val);
	
	// check whther theta near zero or pi
	if (abs(theta) < 1e-5) return Vec3::Zero();
	if (abs(EIGEN_PI - theta) < 1e-5) {
		// TODO: consider mat(2, 2) == -1 and deal with it
		return EIGEN_PI * (1 / sqrt(2 * (1 + mat(2, 2)))) * Vec3(mat(0, 2), mat(1, 2), 1 + mat(2, 2));
	}
	
	// normal situation
	return theta * Vee(mat - mat.transpose()) / (2. * sin(theta));
}

template<Mat4Convertible Derived>
Eigen::Matrix4<typename Derived::Scalar> 
InverseSE3(const Eigen::MatrixBase<Derived>& mat)
{
	using SE3	 = Eigen::Matrix4<typename Derived::Scalar>;
	SE3 _inv_mat = SE3::Identity();
	
	_inv_mat.block(0, 0, 3, 3) = mat.block(0, 0, 3, 3).transpose();
	_inv_mat.block(0, 3, 3, 1) = -_inv_mat.block(0, 0, 3, 3) * mat.block(0, 3, 3, 1);

	return _inv_mat;
}

template<Vec6Convertible Derived>
pair<Twist<typename Derived::Scalar>, typename Derived::Scalar> 
GetTwistAxisAngle(const Eigen::MatrixBase<Derived>& t)
{
	using Scalar = typename Derived::Scalar;
	using Vec3   = Eigen::Vector3<Scalar>;
	using Twist  = Eigen::Vector<Scalar, 6>;
	
	Twist	xi		= t;
    Vec3	w		= t.block(0, 0, 3, 1), 
			v		= t.block(3, 0, 3, 1);

	// caculator theta and check whether not contain rotate or close to zero
	Scalar	theta	= w.norm();
	if (abs(theta) < 1e-5) theta	= v.norm();
	if (abs(theta) > 1e-5) xi		/= theta;

	return std::make_pair(xi, theta);
}

template<Vec6Convertible Derived, ScalarEquivalence<Derived> U>
SE3<typename Derived::Scalar> 
ExpMapping(const Eigen::MatrixBase<Derived>& norm_t, U theta)
{
	using Scalar = typename Derived::Scalar;
	using Vec3	 = Eigen::Vector3<Scalar>;
	using Mat3   = Eigen::Matrix3<Scalar>;
	using SE3	 = Eigen::Matrix4<Scalar>;

    Vec3	w				= norm_t.block(0, 0, 3, 1), 
			v				= norm_t.block(3, 0, 3, 1);
	SE3		mat				= SE3::Identity();

	// Get the rotation SO3 part
	Mat3	skew_sym_m		= Hat(w.normalized());
	mat.block(0, 0, 3, 3)	= Roderigues(skew_sym_m, theta);

	// Get the translation T3 part
	Mat3	G = theta * Mat3::Identity() + (1 - cos(theta)) * skew_sym_m + 
				(theta - sin(theta)) * skew_sym_m * skew_sym_m;
	mat.block(0, 3, 3, 1) = G * v;
	return mat;
}

template<Vec6Convertible Derived>
inline auto ExpMapping(const Eigen::MatrixBase<Derived> & t)
{
	auto && [xi, theta] = GetTwistAxisAngle(t);
	return ExpMapping(xi, theta);
}

template<Vec6Convertible Derived, ScalarEquivalence<Derived> U>
SE3<typename Derived::Scalar> 
DiffExpMapping(const Eigen::MatrixBase<Derived>& norm_t, U theta)
{
	using Scalar = typename Derived::Scalar;
	using Vec3	 = Eigen::Vector3<Scalar>;
	using Mat3   = Eigen::Matrix3<Scalar>;
	using SE3    = Eigen::Matrix4<Scalar>;

    Vec3	w				= norm_t.block(0, 0, 3, 1), 
			v				= norm_t.block(3, 0, 3, 1);
    SE3		mat				= SE3::Identity();

	// Get rotation part in differential form
    Mat3	skew_sym_m		= Hat(w.normalized());
    mat.block(0, 0, 3, 3)	= DiffRoderigues(skew_sym_m, theta);

	// Get translation part in differential form
	Mat3	G				= Mat3::Identity() + (1 + sin(theta)) * skew_sym_m + 
							  (1 - cos(theta)) * skew_sym_m * skew_sym_m;
    mat.block(0, 3, 3, 1) = G * v;

    return mat;
}

template<Vec3Convertible Derived, ScalarEquivalence<Derived> U>
Twist<typename Derived::Scalar> 
ScrewToTwist(const Eigen::MatrixBase<Derived> & q, const Eigen::MatrixBase<Derived>& w, U h = 0.0)
{
	using Twist = Eigen::Vector<typename Derived::Scalar, 6>;
	Twist _normal_twist				= Twist::Zero();
	_normal_twist.block(0, 0, 3, 1) = w.normalized();
	_normal_twist.block(3, 0, 3, 1) = q.cross(w) + h * w.normalized();
	return _normal_twist;
}

template<Mat4Convertible Derived>
Twist<typename Derived::Scalar> 
LogMapSE3Tose3(const Eigen::MatrixBase<Derived>& T)
{
	assert(T.rows() == 4 && T.cols() == 4 && "logrithmic: SE3 -> se3, matrix size not matching");
	using Scalar = typename Derived::Scalar;
	using Vec3	 = Eigen::Vector3<Scalar>;
	using Mat3   = Eigen::Matrix3<Scalar>;
	using Twist  = Eigen::Vector<Scalar, 6>;

    Vec3	t				= T.block(0, 3, 3, 1),	 
			w				= LogMapSO3Toso3(T.block(0, 0, 3, 3));
	auto&& [axis, theta]	= GetRotateAxisAngle(w);
	Mat3	skew_sym_m		= Hat(axis);
	
	Scalar	cot_half_theta	= abs(theta) < 1e-4 ? 1.0 : 0.5 * theta / tan(0.5 * theta);
	Mat3	G_Inv			= Mat3::Identity() - 0.5 * theta * skew_sym_m +
							  (1 - cot_half_theta) * skew_sym_m * skew_sym_m;
	Vec3	v				= G_Inv * t;

	Twist	xi				= Twist::Zero();
	xi.block(0, 0, 3, 1)	= w;
	xi.block(3, 0, 3, 1)	= v;

	return xi;
}

template<Mat4Convertible Derived>
AdMatrix<typename Derived::Scalar> Adjoint(const Eigen::MatrixBase<Derived> & T)
{
	assert(T.rows() == 4 && T.cols() == 4 && "adjoint: SE3 -> adjSE3, matrix size not matching");
	using Scalar = typename Derived::Scalar;
	using Vec3	 = Eigen::Vector3<Scalar>;
	using Mat3	 = Eigen::Matrix3<Scalar>;
	using AdjMat = Eigen::Matrix<Scalar, 6, 6>;

	Mat3	R				= T.block(0, 0, 3, 3);
	Mat3	skew_sym_p		= Hat(T.block(0, 3, 3, 1));
	AdjMat	adMat			= AdjMat::Zero();

	adMat.block(0, 0, 3, 3) = adMat.block(3, 3, 3, 3) 
							= R;
	adMat.block(3, 0, 3, 3) = skew_sym_p * R;
		
	return adMat;
}

/// <summary>
/// Decomposing 3 x 3 matrix into two 3 x 3 matrices of orthogonal matrix Q and upper triangular matrix R
/// <para>
/// 将一个可逆 3 x 3 实矩阵分解为 Q R 两个矩阵，其中 Q 为单位正交阵， R 为上三角矩阵 
/// </para>
/// </summary>
/// <param name="mat"></param>
/// <returns></returns>
template<Mat3Convertible Derived>
pair<SO3<typename Derived::Scalar>, Eigen::Matrix3<typename Derived::Scalar>> 
QRDecompositionMat3(const Eigen::MatrixBase<Derived>& mat)
{
	using Scalar = typename Derived::Scalar;
	using Vec3 = Eigen::Vector3<Scalar>;
	using Mat3 = Eigen::Matrix3<Scalar>;

	auto GetInv = [](auto val)->Scalar {
		return abs(val) > std::numeric_limits<Scalar>::epsilon() ? 1.0 / val : 0.0;
	};
	Vec3		m_cow_1		= mat.block(0, 0, 3, 1),
				m_cow_2		= mat.block(0, 1, 3, 1),
				m_cow_3		= mat.block(0, 2, 3, 1);
	
	// Step1: Schmidt Orthogonalization get Q Matrix
	Vec3		Q_cow_1, Q_cow_2, Q_cow_3;
	Scalar		inv_length	= 0.0;

	inv_length			= GetInv(m_cow_1.norm());
	Q_cow_1				= m_cow_1 * inv_length;

	Q_cow_2				= m_cow_2 - Q_cow_1.dot(m_cow_2) * Q_cow_1;
	inv_length			= GetInv(Q_cow_2.norm());
	Q_cow_2				= Q_cow_2 * inv_length;

	Q_cow_3				= m_cow_3 - Q_cow_1.dot(m_cow_3) * Q_cow_1 - Q_cow_2.dot(m_cow_3) * Q_cow_2;
	inv_length			= GetInv(Q_cow_3.norm());
	Q_cow_3				= Q_cow_3 * inv_length;

	SO3<Scalar>	Q;
	Q.block(0, 0, 3, 1) = Q_cow_1;
	Q.block(0, 1, 3, 1) = Q_cow_2;
	Q.block(0, 2, 3, 1) = Q_cow_3;
	if (Q.determinant() < 0.0) Q = -Q;

	// Step2: Get R use R = Q^t * M
	Mat3		R		= Q.transpose() * mat;

	return { Q, R };
}

/// <summary>
/// Decoposing a 3 x 3 upper triangular Matrix into two 3 x 1 vectors of scale and shear 
/// <para>分解 3 x 3 上三角实矩阵为 缩放 和 剪切 两个 3 x 1 向量 </para>
/// </summary>
/// <typeparam name="_Scaler">	{float/double} precision	</typeparam>
/// <param name="ut_mat"> 　
/// <para>　　{const ref Matrix3} </para>
/// <para>　　upper triangular Matrix 上三角实矩阵</para>
/// </param>
/// <returns></returns>
template<Mat3Convertible Derived>
pair<Eigen::Vector3<typename Derived::Scalar>, Eigen::Vector3<typename Derived::Scalar>> 
SSDecompositionMat3(const Eigen::MatrixBase<Derived>& ut_mat)
{	
	using Vec3 = Eigen::Vector3<typename Derived::Scalar>;
	/*		 x,								   y,							z								*/	
	return  {Vec3(ut_mat(0, 0),				   ut_mat(1, 1),				ut_mat(2, 2)),					// scale vector
			 Vec3(ut_mat(0, 1) / ut_mat(0, 0), ut_mat(0, 2) / ut_mat(0, 0), ut_mat(1, 2) / ut_mat(1, 1))};	// shear vector
}

template<Mat3Convertible Derived>
tuple<SO3<typename Derived::Scalar>, Eigen::Vector3<typename Derived::Scalar>, Eigen::Vector3<typename Derived::Scalar>> 
RSSDecompositionMat3(const Eigen::MatrixBase<Derived>& mat)
{
	// Step1: Decompose mat into Q and R matrices 
	auto [Q, R]			= QRDecompositionMat3(mat);
	// Step2: Decompose R into scale and shear vectors
	auto [scale, shear] = SSDecompositionMat3(R);

	return { Q, scale, shear};
}

template<Mat4Convertible Derived>
pair<Eigen::Matrix3<typename Derived::Scalar>, Eigen::Vector3<typename Derived::Scalar>> 
RtDecompositionMat4(const Eigen::MatrixBase<Derived>& mat)
{
	using Vec3 = Eigen::Vector3<typename Derived::Scalar>;
	using Mat3 = Eigen::Matrix3<typename Derived::Scalar>;

	Mat3	R	= mat.block(0, 0, 3, 3);
	Vec3	t	= mat.block(0, 3, 3, 1);

	return { R, t };
}

template<Mat4Convertible Derived>
tuple<Eigen::Vector3<typename Derived::Scalar>, Eigen::Vector3<typename Derived::Scalar>, 
	  Eigen::Vector3<typename Derived::Scalar>, Eigen::Vector3<typename Derived::Scalar>>
TRSSDecompositionMat4(const Eigen::MatrixBase<Derived>& mat) 
{
	using Vec3 = Eigen::Vector3<typename Derived::Scalar>;
	using Mat3 = Eigen::Matrix3<typename Derived::Scalar>;

	auto [R, t]					 = RtDecompositionMat4(mat);
	auto [rot_mat, scale, shear] = RSSDecompositionMat3(R);
	return { t, LogMapSO3Toso3(rot_mat), scale, shear};			
}

template<Mat4Convertible Derived>
pair<Eigen::Vector3<typename Derived::Scalar>, Eigen::Vector3<typename Derived::Scalar>> 
rtDecompositionMat4(const Eigen::MatrixBase<Derived>& mat)
{
	using Vec3 = Eigen::Vector3<typename Derived::Scalar>;
	using Mat3 = Eigen::Matrix3<typename Derived::Scalar>;

	Vec3	t	= mat.block(0, 3, 3, 1),
			r	= LogMapSO3Toso3(mat.block(0, 0, 3, 3));
	return {r, t};
}

}

#endif	// _GTRANSFORM_HPP
