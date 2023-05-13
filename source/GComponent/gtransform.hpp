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
#include <GComponent/types.h>

#include <Eigen/Core>
#include <Eigen/LU>
#include <Eigen/Geometry>

#include <tuple>

namespace GComponent {

using Eigen::Matrix;
using Eigen::Vector;
using std::pair;
using std::tuple;


/// <summary>
/// make 3 x 3 shear matrix form 3 x 1 shear vector
/// <para>
/// ʹ�� 3 x 1 ����������ȡ 3 x 3 ���о���
/// </para>
/// </summary>
/// <param name="shear">	cref	{vec3}	3 x 1 shear vctor	3 x 1 ��������		</param>
/// <returns>				val		{mat3}	3 x 3 shear matrix	3 x 3 ���о���		</returns>
template<Vec3Convertible Derived>
auto Shear(const Eigen::MatrixBase<Derived>& shear)
{
	assert(shear.rows() == 3 && shear.cols() == 1 && "Shear operator: R^3 -> R^3 x 3, only support 3 x 1 vector, the size not matching");
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
/// ʹ�� 3 x 1 ����������ȡ 3 x 3 ���������
/// </para>
/// </summary>
/// <param name="shear">	cref	{vec3}	3 x 1 shear vctor			3 x 1 ��������			</param>
/// <returns>				val		{mat3}	3 x 3 inverse shear matrix	3 x 3 ���������		</returns>
template<Vec3Convertible Derived>
auto ShearInverse(const Eigen::MatrixBase<Derived>& shear) 
{
	assert(shear.rows() == 3 && shear.cols() == 1 && 
		   "Inv Shear operator: R^3 -> R^3 x 3, only support 3 x 1 vector, the size not matching");
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
/// ʹ�� 3 x 1 ����������ȡ 3 x 3 ���ž���
/// </para>
/// </summary>
/// <param name="scale">	cref	{vec3}	3 x 1 scale vctor	3 x 1 ��������		</param>
/// <returns>				val		{mat3}	3 x 3 scale matrix	3 x 3 ���ž���		</returns>
template<Vec3Convertible Derived>
auto Scale(const Eigen::MatrixBase<Derived>& scale)
{
	assert(scale.rows() == 3 && scale.cols() == 1 && 
		   "Scale operator: R^3 -> R^3 x 3, only support 3 x 1 vector, the size not matching");
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
/// ʹ�� 3 x 1 ����������ȡ 3 x 3 ���������
/// </para>
/// </summary>
/// <param name="scale">	cref	{vec3}	3 x 1 scale vctor			3 x 1 ��������			</param>
/// <returns>				val		{mat3}	3 x 3 inverse scale matrix	3 x 3 ���������		</returns>
template<Vec3Convertible Derived>
auto ScaleInverse(const Eigen::MatrixBase<Derived>& scale)
{
	assert(scale.rows() == 3 && scale.cols() == 1 && 
		   "Inv Scale operator: R^3 -> R^3 x 3, only support 3 x 1 vector, the size not matching");
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
/// ��һ�� 3 x 1 ����ʩ����ξ���任������ֱ���ϴ����Ŷ�һ����ʩ����ת��ƽ�Ʊ任
/// </para>
/// </summary>
/// <param name="mat">	cref	{SE3}	homogeneous matrix in SE3	��α任����					</param>
/// <param name="vec">	cref	{vec3}	vector in R^3				3 x 1 ����						</param>
/// <returns>			val		{vec3}	new vector after transform	�任���������					</returns>
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
/// ����Ǳ�ʾ��������Ϊ���򻯵�����������ʵ����ʾ�Ļ��Ƚ�
/// </para>
/// </summary>
/// <param name="vec">	cref	{vec3}		[in]	an 3 x 1 Vector							3 x 1 ����			</param>
/// <returns>			val		{vec3, R}	[out]	norm angle dir and real number angle	����������ʵ����	</returns>
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
/// Hat ���ӽ� so3 �� 3 x 1 ��������ת���� 3 x 3 ���Գƾ������
/// </para>
/// </summary>
/// <param name="vec">	cref	{so3}	[in]	an so3 vector form			3 x 1 ��Ǳ���			</param>
/// <returns>			val		{[so3]}	[out]	an skew-symmatrix matrix	3 x 3 ���Գƾ���		</returns>
template<Vec3Convertible Derived>
auto Hat(const Eigen::MatrixBase<Derived>& vec)
{
	assert(vec.cols() == 1 && vec.rows() == 3 && 
		   "Hat operator: R^3 -> [so3], only support 3 x 1 vector, the size not matching");
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
/// Vee ���ӽ� so3 �� 3 x 3 ���Գƾ������ת���� 3 x 1 ��������
/// </para>
/// </summary>
/// <param name="skew_sym_mat">	cref	{[so3]}	[in]	an skew-symmatrix matrix	3 x 3 ���Գƾ���	</param>
/// <returns>					val		{so3}	[out]	so3 element vector form		3 x 1 ��Ǳ���		</returns>
template<Mat3Convertible Derived>
inline Eigen::Vector3<typename Derived::Scalar> 
Vee(const Eigen::MatrixBase<Derived>& skew_sym_mat)
{
    return Eigen::Vector3<typename Derived::Scalar>(skew_sym_mat(2, 1), skew_sym_mat(0, 2), skew_sym_mat(1, 0));
}

/// <summary>
/// Using Roderigues formula to get SO3 element from [so3] axis and real number angle in [rad.]
/// <para>
/// ʹ�� Roderigues ��ʽ�ӷ��Գƾ���ת���ʵ������ת�ǻ�ȡ��Ӧ��ת����
/// </para>
/// </summary>
/// <param name="skew_sym_m">	cref	{[so3]}	[in]	an skew-symmatrix matrix	3x3 ���Գƾ���	</param>
/// <param name="theta">		val		{R}		[in]	a real number angle [rad.]	ʵ�����Ƚ�		</param>
/// <returns>					val		{SO3}	[out]	rotation matrix				��ת����		</returns>
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
/// ʹ�� Roderigues ��ʽ��ת��������ʵ������ת�ǻ�ȡ��Ӧ��ת����
/// </para>
/// </summary>
/// <param name="axis">			cref	{so3}	[in]	an 3 x 1 vector				3 x 1 ת������	</param>
/// <param name="theta">		val		{R}		[in]	a real number angle [rad.]	ʵ�����Ƚ�		</param>
/// <returns>					val		{SO3}	[out]	rotation matrix				��ת����		</returns>
template<Vec3Convertible Derived, ScalarEquivalence<Derived> U>
inline auto Roderigues(const Eigen::MatrixBase<Derived>& axis, U theta)
{
	return Roderigues(Hat(axis.normalized()), theta);
}

/// <summary>
/// Using Roderigues formula to get SO3 matrix from so3 axis angle vector with unit [rad.]
/// <para>
/// ʹ�� Roderigues ��ʽ�����������ȡ��Ӧ��ת����
/// </para>
/// </summary>
/// <param name="v">	cref	{so3}	[in]	3 x 1 axis angle vecotr in [rad.]	3 x 1 �������	</param>
/// <returns>			val		{SO3}	[out]	3 x 3 rotate matrix					3 x 3 ��ת����	</returns>
template<Vec3Convertible Derived>
inline auto Roderigues(const Eigen::MatrixBase<Derived>& v)
{
	return Roderigues(Hat(v.normalized()), v.norm());
}

/// <summary>
/// Using Roderigues formula's diiferential formula to get matrix from so3 axis angle vector with unit [rad.]
/// <para>
/// ���� Roderigues ��ʽ�����������ȡ��Ӧ��΢����ת����
/// </para>
/// </summary>
/// <param name="skew_sym_m">	cref	{[so3]}		[in]	3 x 3 skew symmetric matrix		3 x 3 ��Ǿ���		</param>
/// <param name="theta">		val		{R}			[in]	a real number angle in [rad.]	ʵ�����Ƚ�			</param>
/// <returns>					val		{R^3 x 3}	[out]	3 x 3 differential matrix		3 x 3 ��ת����΢��	</returns>
template<Mat3Convertible Derived, ScalarEquivalence<Derived> U>
inline Eigen::Matrix3<typename Derived::Scalar> 
DiffRoderigues(const Eigen::MatrixBase<Derived>& skew_sym_m, U theta)
{
	using Mat3 = Eigen::Matrix3<typename Derived::Scalar>;
    return  Mat3::Identity() + cos(theta) * skew_sym_m + sin(theta) * (skew_sym_m * skew_sym_m);
}

/// <summary>
/// Logrithimic map : get a so3 vector from a SO3 rotation matrix
/// <para>
/// ʹ�ö���ӳ�����ת�����ȡ��Ӧ���������
/// </para>
/// </summary>
/// <param name="mat">	cref	{SO3}	[in]	3 x 3 rotation matrix			3 x 3 ��ת����		</param>
/// <returns>			val		{so3}	[out]	3 x 1 angle axis vector	[rad.]	3 x 1 �����������	</returns>
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

/// <summary>
/// get inverse of SE3 homogeneous transform matrix
/// <para>
/// ��ȡ��α任�������
/// </para>
/// </summary>
/// <param name="mat">	cref	{SE3}	[in]	4 x 4 homogeneous transform matrix		4 x 4 ��α任����	</param>
/// <returns>			val		{SE3}	[out]	inverse of 4 x 4 input matrix			4 x 4 ����������	</returns>
template<Mat4Convertible Derived>
auto InverseSE3(const Eigen::MatrixBase<Derived>& mat)
{
	using SE3	 = Eigen::Matrix4<typename Derived::Scalar>;
	SE3 _inv_mat = SE3::Identity();
	
	_inv_mat.block(0, 0, 3, 3) = mat.block(0, 0, 3, 3).transpose();
	_inv_mat.block(0, 3, 3, 1) = -_inv_mat.block(0, 0, 3, 3) * mat.block(0, 3, 3, 1);

	return _inv_mat;
}

/// <summary>
/// Decompose a twist into axis and angle components
/// <para>
/// ��һ����������Ϊ��Ƿ�������ʽ
/// </para>
/// </summary>
/// <param name="t">	cref	{twist}		[in]	6 x 1 twist vector								6 x 1 ����					</param>
/// <returns>			val		{twist, R}	[out]	normalized twist direction with angle in [rad.] ��������������ʵ�����Ƚ�	</returns>
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

/// <summary>
/// Exponential map : get a SE3 homogeneous transform matrix form a se3 vector with angle
/// <para>
/// ָ��ӳ�� : ��һ�� se3 �����뻡�Ƚǻ�ȡ SE3 ��α任����
/// </para>
/// </summary>
/// <param name="norm_t">	cref	{se3}	[in]	normalized 6 x 1 twist vector		���򻯵� 6 x 1 ����	</param>
/// <param name="theta">	val		{R}		[in]	the real number angle in [rad.]		����ʵ����			</param>
/// <returns>				val		{SE3}	[out]	4 x 4 homogeneous transform matrix	4 x 4 ��α任����	</returns>
template<Vec6Convertible Derived, ScalarEquivalence<Derived> U>
auto ExpMapping(const Eigen::MatrixBase<Derived>& norm_t, U theta)
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

/// <summary>
/// Exponential map : get a SE3 homogeneous transform matrix form a se3 vector
/// <para>
/// ָ��ӳ�� : ��һ�� se3 ������ȡ SE3 ��α任����
/// </para>
/// </summary>
/// <param name="t">	cref	{se3}	[in]	6 x 1 twist vector in [rad.]		������ 6 x 1 ����	</param>
/// <returns>			val		{SE3}	[out]	4 x 4 homogeneous transform matrix	4 x 4 ��α任����	</returns>
template<Vec6Convertible Derived>
inline auto ExpMapping(const Eigen::MatrixBase<Derived> & t)
{
	auto && [xi, theta] = GetTwistAxisAngle(t);
	return ExpMapping(xi, theta);
}

/// <summary>
/// Differential Exponential map: get a differential theta transform matrix from se3 vector and real number anlge [rad.]
/// <para>
/// ΢��ָ��ӳ�� : ��ȡ��α任������ڽǶȵ�΢��
/// </para>
/// </summary>
/// <param name="norm_t">	cref	{se3}		[in]	normalized 6 x 1 twist vector		���򻯵� 6 x 1 ����	</param>
/// <param name="theta">	val		{R}			[in]	the real number angle in [rad.]		����ʵ����			</param>
/// <returns>				val		{R^3 x 3}	[out]	the differential matrix				΢�־���			</returns>
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
	Mat3	G				= Mat3::Identity() + sin(theta) * skew_sym_m + 
							  (1 - cos(theta)) * skew_sym_m * skew_sym_m;
    mat.block(0, 3, 3, 1) = G * v;

    return mat;
}

/// <summary>
/// get the twist vector from screw q, w and h components
/// <para>
/// ��������λ�á�������ھ�����л�ȡ����
/// </para>
/// </summary>
/// <param name="q">	cref	{R^3}	[in]	screw axis position		��������һ��λ��	</param>
/// <param name="w">	cref	{R^3}	[in]	screw axis direction	�����᷽��			</param>
/// <param name="h">	val		{R}		[in]	screw pitch				�ھ�				</param>
/// <returns>			val		{se3}	[out]	6 x 1 twist vector		6 x 1 ��������		</returns>
template<Vec3Convertible Derived, ScalarEquivalence<Derived> U>
auto ScrewToTwist(const Eigen::MatrixBase<Derived> & q, const Eigen::MatrixBase<Derived>& w, U h)
{
	using Twist = Eigen::Vector<typename Derived::Scalar, 6>;
	Twist xi			 = Twist::Zero();
	xi.block(0, 0, 3, 1) = w.normalized();
	xi.block(3, 0, 3, 1) = q.cross(w) + h * w.normalized();
	return xi;
}

/// <summary>
/// get the twist vector from screw q and w components without h(.equal h = 0)
/// <para>
/// ��������λ�á������л�ȡ�������ȼ��ڽھ�Ϊ 0 �����
/// </para>
/// </summary>
/// <param name="q">	cref	{R^3}	[in]	screw axis position		��������һ��λ��	</param>
/// <param name="w">	cref	{R^3}	[in]	screw axis direction	�����᷽��			</param>
/// <returns>			val		{se3}	[out]	6 x 1 twist vector		6 x 1 ��������		</returns>
template<Vec3Convertible Derived>
auto ScrewToTwist(const Eigen::MatrixBase<Derived>& q, const Eigen::MatrixBase<Derived>& w)
{
	using Twist = Eigen::Vector<typename Derived::Scalar, 6>;
	Twist xi			 = Twist::Zero();
	xi.block(0, 0, 3, 1) = w.normalized();
	xi.block(3, 0, 3, 1) = q.cross(w);
	return xi;
}

/// <summary>
/// Logrithimic map : get a se3 twist from a SE3 homogeneous transform matrix
/// <para>
/// ʹ�ö���ӳ�����α任�����ȡ��Ӧ������
/// </para>
/// </summary>
/// <param name="T">	cref	{SE3}	[in]	4 x 4 homogeneous transform matrix	��α任����	</param>
/// <returns>			val		{se3}	[out]	6 x 1 twist vector					6 x 1 ����		</returns>
template<Mat4Convertible Derived>
auto LogMapSE3Tose3(const Eigen::MatrixBase<Derived>& T)
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

/// <summary>
/// Adjoint map : get R^6 x 6 adjoint matrix from SE3 4 x 4 homogeneous transform matrix
/// <para>
/// ����任���� 4 x 4 ��α任�����ȡ��Ӧ�� 6 x 6 �������
/// </para>
/// </summary>
/// <param name="T">	cref	{SE3}		[in]	4 x 4 homogeneous transform matrix	��α任����	</param>
/// <returns>			val		{R^6 x 6}	[out]	6 x 6 adjoint matrix				����任����	</returns>
template<Mat4Convertible Derived>
AdMatrix<typename Derived::Scalar> Adjoint(const Eigen::MatrixBase<Derived> & T)
{
	assert(T.rows() == 4 && T.cols() == 4 && "adjoint: SE3 -> adjSE3, matrix size not matching");
	using Scalar = typename Derived::Scalar;
	using Vec3	 = Eigen::Vector3<Scalar>;
	using Mat3	 = Eigen::Matrix3<Scalar>;
	using AdjMat = Eigen::Matrix<Scalar, 6, 6>;

	Mat3	R				= T.block(0, 0, 3, 3);
	Mat3	skew_sym_p		= Hat(static_cast<Vec3>(T.block(0, 3, 3, 1)));
	AdjMat	adMat			= AdjMat::Zero();

	adMat.block(0, 0, 3, 3) = adMat.block(3, 3, 3, 3) 
							= R;
	adMat.block(3, 0, 3, 3) = skew_sym_p * R;
		
	return adMat;
}

/// <summary>
/// Decomposing 3 x 3 matrix into two 3 x 3 matrices of orthogonal matrix Q and upper triangular matrix R using Schmidt Orthogonalization
/// <para>
/// ʹ�� Schmidt ����������һ������ 3 x 3 ʵ����Ϊ Q R ������������ Q Ϊ��λ������ R Ϊ�����Ǿ��� 
/// </para>
/// </summary>
/// <param name="mat">	cref	{R^3x3}			[in]	3 x 3 real number matrix				3 x 3 ʵ����				</param>
/// <returns>			val		{SO3, R^3x3}	[out]	rotate matrix and upper triangle matrix ��ת���������������Ǿ���	</returns>
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
/// Decompose a 3 x 3 upper triangular Matrix into two 3 x 1 vectors of scale and shear 
/// <para>
/// ���� 3 x 3 ������ʵ����Ϊ���� 3 x 1 ���������źͼ��� 
/// </para>
/// </summary>
/// <param name="ut_mat">	cref	{R^3x3}		[in]	3 x 3 upper triangle matrix				3 x 3 �����Ǿ���	</param>
/// <returns>				val		{R^3, R^3}	[out]	pair of 3 x 1 vectors: scale and shear	�������������		</returns>
template<Mat3Convertible Derived>
pair<Eigen::Vector3<typename Derived::Scalar>, Eigen::Vector3<typename Derived::Scalar>> 
SSDecompositionMat3(const Eigen::MatrixBase<Derived>& ut_mat)
{	
	using Vec3 = Eigen::Vector3<typename Derived::Scalar>;
	/*		 x,								   y,							z								*/	
	return  {Vec3(ut_mat(0, 0),				   ut_mat(1, 1),				ut_mat(2, 2)),					// scale vector
			 Vec3(ut_mat(0, 1) / ut_mat(0, 0), ut_mat(0, 2) / ut_mat(0, 0), ut_mat(1, 2) / ut_mat(1, 1))};	// shear vector
}

/// <summary>
/// Decompose a real number 3 x 3 matrix into Rotate matrix, scale vector and shear vector
/// <para>
/// ���� 3 x 3 ʵ����Ϊһ����ת���������������������
/// </para>
/// </summary>
/// <param name="mat">		cref	{R^3x3}			[in]	3 x 3 real number matrix			3 x 3 ʵ������			</param>
/// <returns>				val		{SO3, R^3, R^3}	[out]	tuple of rot, scale and shear		��ת�����������Ԫ��	</returns>
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

/// <summary>
/// Decompose 4 x 4(R^3x3 x t^3) affine transform matrix into R^3x3 matrix and t^3 vector
/// <para> 
/// �� 4 x 4 ����任�������Ϊ 3 x 3 ������ 3 x 1 ����
/// </para>
/// </summary>
/// <param name="mat">	cref	{R^4x4}			[in]	4 x 4 affine transform matrix	4 x 4 ����任����	</param>
/// <returns>			val		{R^3x3, R^3}	[out]	pair of R matrix and t vecotr	R ������ t ������	</returns>
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

/// <summary>
/// Decompose 4 x 4 affine transform matrix into trans, rot, scale, shear components tuple
/// <para>
/// ���� 4 x 4 �������Ϊ����ƽ�ơ���ת����������е� 3 x 1 ����Ԫ��
/// </para>
/// </summary>
/// <param name="mat">	cref	{R^4x4}			[in]	4 x 4 afiine transform matrix				4 x 4 ����任����	</param>
/// <returns>			val		{R^3,...}_4		[out]	3 x 1 vector tuple[t, rot, scale, shear]	3 x 1 ����Ԫ��		</returns>
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

/// <summary>
/// Decompose 4 x 4 homogeneous transform matrix into 3 x 1 vectors of trans and rot
/// <para>
/// ���� 4 x 4 ��α任����Ϊƽ������ת���� 3 x 1 ����
/// </para>
/// </summary>
/// <param name="mat">	cref	{SE3}		[in]	4 x 4 homogeneous transform matrix		4 x 4 ��α任����	</param>
/// <returns>			val		{R^3, R^3}	[out]	3 x 1 vectors pair of trans and rot		3 x 1 ����Ԫ�� {ƽ��, ��ת}</returns>
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


/*____________________________________Rotation Process ��ת����________________________________*/
template<class Scalar>
Eigen::Matrix3<Scalar> RotZ(Scalar q) {
	Eigen::Matrix3<Scalar> rot;
	rot << cos(q), -sin(q), 0,
		   sin(q), cos(q),  0,
			    0,      0,  1;
	return rot;
}

template<class Scalar>
Eigen::Matrix3<Scalar> RotY(Scalar q) {
	Eigen::Matrix3<Scalar> rot;
	rot << cos(q), 0, sin(q),
			    0, 1,      0,
		  -sin(q), 0, cos(q);
	return rot;
}

template<class Scalar>
Eigen::Matrix3<Scalar> RotX(Scalar q) {
	Eigen::Matrix3<Scalar> rot;
	rot <<1,      0,       0,
		  0, cos(q), -sin(q),
		  0, sin(q),  cos(q);		
	return rot;
}

template<Vec3Convertible Derived>
Eigen::Vector3<typename Derived::Scalar> 
ToZYXEuler(const Eigen::MatrixBase<Derived>& vec) {
	assert(vec.rows() == 3 && vec.cols() == 1 && "so3 only support 3 x 1 vector, the size not matching");
	using Scalar = typename Derived::Scalar;
	using Vec3   = Eigen::Vector3<Scalar>;
		
	Scalar angle  = vec.norm();	
	Vec3   axis   = vec.normalized();
				
	Scalar sin_v = std::sin(angle / 2.0);
	Scalar cos_v = std::cos(angle / 2.0);

	Scalar x = axis.x() * sin_v;
	Scalar y = axis.y() * sin_v;
	Scalar z = axis.z() * sin_v;
	Scalar w = cos_v;

	Scalar eulerX = std::atan2(2 * (x * w + y * z), w * w - x * x - y * y + z * z);
	Scalar eulerY = std::asin(-2 * (x * z - y * w));
	Scalar eulerZ = std::atan2(2 * (x * y + z * w), w * w + x * x - y * y - z * z);

	Vec3 euler(eulerX, eulerY, eulerZ);
	return euler;
}

template<Vec3Convertible Derived>
Eigen::Vector3<typename Derived::Scalar>
FromZYXEuler(const Eigen::MatrixBase<Derived>& vec) {
	assert(vec.rows() == 3 && vec.cols() == 1 && "ZYXEuler only support 3 x 1 vector, the size not matching");
	using Scalar = typename Derived::Scalar;
	
	Eigen::Matrix3<Scalar> mat = RotZ(vec.z()) * RotY(vec.y()) * RotX(vec.x());
	return LogMapSO3Toso3(mat);
}

template<Vec3Convertible Derived>
Eigen::Vector3<typename Derived::Scalar>
ToXYZEuler(const Eigen::MatrixBase<Derived>& vec) {
	assert(vec.rows() == 3 && vec.cols() == 1 && "so3 only support 3 x 1 vector, the size not matching");

	using Scalar = typename Derived::Scalar;
	using Vec3   = Eigen::Vector3<Scalar>;
	using Mat3   = Eigen::Matrix3<Scalar>;

	Scalar angle  = vec.norm();	
	Vec3   axis   = vec.normalized();
				
	Eigen::AngleAxis<Scalar> angle_axis(angle, axis);
	Mat3 mat = angle_axis.toRotationMatrix();

	Scalar euler_x = atan2(-mat(1, 2), mat(2, 2));
	Scalar euler_y = asin(mat(0, 2));
	Scalar euler_z = atan2(-mat(0, 1), mat(0, 0));

	Vec3 euler(euler_x, euler_y, euler_z);
	return euler;
}

template<Vec3Convertible Derived>
Eigen::Vector3<typename Derived::Scalar>
FromXYZEuler(const Eigen::MatrixBase<Derived>& vec) {
	assert(vec.rows() == 3 && vec.cols() == 1 && "XYZEuler only support 3 x 1 vector, the size not matching");
	using Scalar = typename Derived::Scalar;
	
	Eigen::Matrix3<Scalar> mat = RotX(vec.x()) * RotY(vec.y()) * RotZ(vec.z());
	return LogMapSO3Toso3(mat);
}

}

#endif	// _GTRANSFORM_HPP
