#ifndef _GGEOMETRY_HPP
#define _GGEOMETRY_HPP

#include <Concept/gconcept.hpp>
#include <GComponent/types.h>
#include <GComponent/gtransform.hpp>
#include <LSSolver/LinearSystemSolver.hpp>

#include <Eigen/Core>
#include <Eigen/LU>

#include <functional>
#include <algorithm>
#include <numeric>
#include <vector>

// TODO: add comments
namespace GComponent
{

using std::tuple;
using std::vector;
using std::function;
using Eigen::Vector3d;
using Eigen::Matrix3d;
using Eigen::Matrix4d;



/// <summary>
/// 插值函数，目前仅支持列向量插值，待完善为行列均支持
/// </summary>
/// <typeparam name="_AnyVec"></typeparam>
/// <param name="v1"></param>
/// <param name="v2"></param>
/// <param name="t"></param>
/// <returns></returns>
template<class _AnyVec>
_AnyVec Lerp(const _AnyVec & v1, const _AnyVec& v2, double t)
{
	static_assert(_AnyVec::ColsAtCompileTime == 1, "The Function not Support Matrix");
	const unsigned ROW = v1.rows();
	/* This Step is to Deal with Dynamic Vec */
	_AnyVec tmp = _AnyVec::Zero();
    for (unsigned i = 0; i < ROW; ++i)
	{
		tmp[i] = (1 - t) * v1[i] + t * v2[i];
	}
	return tmp;
}

inline double Lerp(double d1, double d2, double t)
{
    return (1 - t) * d1 + t * d2;
}

/// <summary>
/// (1 - t) * a + t * b
/// <para>
/// 线性插值
/// </para>
/// </summary>
/// <param name="a"></param>
/// <param name="b"></param>
/// <param name="t"></param>
/// <returns></returns>
inline float Lerp(float a, float b, float t) {
	return (1 - t) * a + t * b;
}

/// <summary>
/// 将线段空间进行 n 等分函数
/// </summary>
/// <param name="lowerBound"></param>
/// <param name="upperBound"></param>
/// <param name="n"></param>
/// <returns></returns>
template<class Scalar>
vector<Scalar> Linspace(const Scalar lowerBound, const Scalar upperBound,const int kN)
{
	vector<Scalar> rval(kN);
	const Scalar   kStep	= (upperBound - lowerBound) / (kN - 1);
	double		   val_cur  = lowerBound;
	for (auto& val : rval) {
		val		=  val_cur;
		val_cur += kStep;
	}
	
	return rval;
}

template<Vec3Convertible Derived>
Eigen::Vector3<typename Derived::Scalar> 
GetCenterOfCircle(const Eigen::MatrixBase<Derived>& p1, 
				  const Eigen::MatrixBase<Derived>& p2, 
				  const Eigen::MatrixBase<Derived>& p3)
{
	using Scalar  = typename Derived::Scalar;
	using Mat3    = Eigen::Matrix3<Scalar>;
	using Vec3    = Eigen::Vector3<Scalar>;
	using RowVec3 = Eigen::RowVector3<Scalar>;
				
	RowVec3 vec1			  = (p2 - p1).transpose(),
		    vec2			  = (p3 - p1).transpose(),
		    vec_cross		  = vec1.cross(vec2);

	Mat3	aux_mat;
	aux_mat.block(0, 0, 1, 3) = vec1;
	aux_mat.block(1, 0, 1, 3) = vec2;
	aux_mat.block(2, 0, 1, 3) = vec_cross;
		
	Vec3 b = Vec3(vec1.dot((p1 + p2) * static_cast<Scalar>(0.5)),
				  vec2.dot((p3 + p1) * static_cast<Scalar>(0.5)),
				  vec_cross.dot(p3));	
	Vec3 x = aux_mat.inverse() * b;
	return x;
}

/// <summary>
/// 获取从 v1 到 v2 的旋转向量轴角表示
/// v1 -> v2
/// v2 = [rot]^ * v1
/// </summary>
/// <typeparam name="_Scaler"></typeparam>
/// <param name="v1"></param>
/// <param name="v2"></param>
/// <returns></returns>
template<class _Scaler>
Vector<_Scaler, 3> GetRotateAxisAngleFrom2Vec(const Vector<_Scaler, 3>& v1, const Vector<_Scaler, 3>& v2)
{
	Vector<_Scaler, 3> v1_norm = v1.normalized(),
					   v2_norm = v2.normalized();
	if ((v1_norm - v2_norm).norm() < std::numeric_limits<_Scaler>::epsilon()) return Vector<_Scaler, 3>::Zero();
	Vector<_Scaler, 3> n = v1_norm.cross(v2_norm).normalized();
	_Scaler angle = acos(
		std::max(std::min(v1_norm.dot(v2_norm), static_cast<_Scaler>(1.0)), static_cast<_Scaler>(-1.0)));

	return (angle * n).eval();
}

template<Vec3Convertible Derived>
auto GetRadiusOfCircle(const Eigen::MatrixBase<Derived>& p1, const Eigen::MatrixBase<Derived>& p2, const Eigen::MatrixBase<Derived>& p3)
{
	using Vec3 = Eigen::Vector3<typename Derived::Scalar>;
    Vec3 coc   = GetCenterOfCircle(p1, p2, p3);
    return ( p1 - coc).norm();
}

template<Vec3Convertible Derived>
typename Derived::Scalar 
GetArcDeltaOfCircle(const Eigen::MatrixBase<Derived>& p_ini, 
					const Eigen::MatrixBase<Derived>& p_mid,
					const Eigen::MatrixBase<Derived>& p_end)
{
	using Scalar = typename Derived::Scalar;
	using Vec3   = Eigen::Vector3<Scalar>;

    Vec3   axis     = ((p_mid - p_ini).cross(p_end - p_ini)).normalized(),
		   center   = GetCenterOfCircle(p_ini, p_mid, p_end);
    Vec3   vec_ini  = (p_ini - center).normalized(),
           vec_end  = (p_end - center).normalized();

    Scalar angle	= acos(vec_ini.dot(vec_end));
	bool   same_dir = axis.dot(vec_ini.cross(vec_end)) > 0.;
    return same_dir? angle: 2 * EIGEN_PI - angle;
}

template<Mat4Convertible Derived>
function<Twist<typename Derived::Scalar>(typename Derived::Scalar)> 
GetScrewLineFunction(const Eigen::MatrixBase<Derived>& T_ini, const Eigen::MatrixBase<Derived>& T_end)
{
	using Scalar = typename Derived::Scalar;
	using Vec3   = Eigen::Vector3<Scalar>;
	using Twist  = Eigen::Vector<Scalar, 6>;
	using SE3    = Eigen::Matrix4<Scalar>;

    const Vec3& pos_ini = T_ini.block(0, 3, 3, 1), &pos_end = T_end.block(0, 3, 3, 1);
    Twist		t_ini   = LogMapSE3Tose3(T_ini),    t_end   = LogMapSE3Tose3(T_end);
    const Vec3& w_ini   = t_ini.block(0, 0, 3, 1), &w_end   = t_end.block(0, 0, 3, 1);

    auto lin_func = [pos_ini = pos_ini, pos_end = pos_end, w_ini = w_ini, w_end = w_end](Scalar t)->Twist{
        SE3   T_cur   = SE3::Identity();
        Vec3  pos_cur = Lerp(pos_ini, pos_end, t);
        Vec3  w_cur   = Lerp(w_ini, w_end, t);
        T_cur.block(0, 0, 3, 3) = Roderigues(w_cur);
        T_cur.block(0, 3, 3, 1) = pos_cur;
        return LogMapSE3Tose3(T_cur);
    };

    return lin_func;
}

template<Vec6Convertible Derived>
auto 
GetScrewLineFunction(const Eigen::MatrixBase<Derived> & t_ini, const Eigen::MatrixBase<Derived>& t_end)
{
    return GetScrewLineFunction(ExpMapping(t_ini), ExpMapping(t_end));
}

template<Vec3Convertible Derived>
auto
GetCircleFunction(const Eigen::MatrixBase<Derived>& p_ini, 
				  const Eigen::MatrixBase<Derived>& p_end, 
				  const Eigen::MatrixBase<Derived>& p_mid)
{
	using Scalar = typename Derived::Scalar;
	using Vec3   = Eigen::Vector3<Scalar>;

	Vec3   axis	  = ((p_mid - p_ini).cross(p_end - p_ini)).normalized(),
		   center = GetCenterOfCircle(p_ini, p_mid, p_end);
	Scalar radius = (p_ini - center).norm();

	Vec3  vec_ini = (p_ini - center).normalized(),
		  vec_end = (p_end - center).normalized();

	bool   same_dir = axis.dot(vec_ini.cross(vec_end)) > 0.;
	Scalar angle    = acos(vec_ini.dot(vec_end));
	angle = same_dir ? angle : 2 * EIGEN_PI - angle;

	auto circle_func = [total  = angle,  axis = axis,		center = center, 
						radius = radius, vec_ini = vec_ini]
	(Scalar t)->Vec3 {
		return center + radius * (Roderigues(axis, t * total) * vec_ini);
	};
	return circle_func;
}

template<Mat4Convertible MatDerived, Vec3Convertible VecDerived>
requires ScalarSame<MatDerived, VecDerived>
auto
GetCircleFunction(const Eigen::MatrixBase<MatDerived>&   T_ini,
				  const Eigen::MatrixBase<MatDerived>&   T_end,
				  const Eigen::MatrixBase<VecDerived>&   p_mid){
	using Scalar = typename MatDerived::Scalar;
	using SE3    = Eigen::Matrix4<Scalar>;
	using Vec3   = Eigen::Vector3<Scalar>;
	using Twist  = Eigen::Vector<Scalar, 6>;

	Vec3  p_ini = T_ini.block(0, 3, 3, 1), p_end = T_end.block(0, 3, 3, 1);
	Vec3  w_ini = LogMapSO3Toso3(T_ini.block(0, 0, 3, 3)),
		  w_end = LogMapSO3Toso3(T_end.block(0, 0, 3, 3));
	auto pos_func	 = GetCircleFunction(p_ini, p_end, p_mid);
	auto circle_func = [PosFunc = pos_func, w_ini = w_ini, w_end = w_end]
	(Scalar t)->Twist {
        SE3 T_cur = SE3::Identity();
        T_cur.block(0, 3, 3, 1) = PosFunc(t);
        T_cur.block(0, 0, 3, 3) = Roderigues(Lerp(w_ini, w_end, t));
        return LogMapSE3Tose3(T_cur);
	};
	return circle_func;
}

template<Vec6Convertible TwistDerived, Vec3Convertible VecDerived>
requires ScalarSame<TwistDerived, VecDerived>
auto
GetCircleFunction(const Eigen::MatrixBase<TwistDerived>& t_ini, 
				  const Eigen::MatrixBase<TwistDerived>& t_end, 
				  const Eigen::MatrixBase<VecDerived>&   p_mid)
{
	using Scalar = typename TwistDerived::Scalar;
	using SE3    = Eigen::Matrix4<Scalar>;
	using Vec3   = Eigen::Vector3<Scalar>;
	using Twist  = Eigen::Vector<Scalar, 6>;

    SE3   T_ini = ExpMapping(t_ini),	   T_end = ExpMapping(t_end);
    Vec3  p_ini = T_ini.block(0, 3, 3, 1), p_end = T_end.block(0, 3, 3, 1);
	Vec3  w_ini = t_ini.block(0, 0, 3, 1), w_end = t_end.block(0, 0, 3, 1);

	auto pos_func	 = GetCircleFunction(p_ini, p_end, p_mid);
	auto circle_func = [PosFunc = pos_func, w_ini = w_ini, w_end = w_end]
	(Scalar t)->Twist {
        SE3 T_cur = SE3::Identity();
        T_cur.block(0, 3, 3, 1) = PosFunc(t);
        T_cur.block(0, 0, 3, 3) = Roderigues(Lerp(w_ini, w_end, t));
        return LogMapSE3Tose3(T_cur);
	};
	return circle_func;
}

template<Vec3Convertible Derived>
auto
GetCubicSplineFunction(const vector<Derived>& poses, 
					   double M0 = 0, 
					   double Mn = 0)
{
	using  Scalar  = typename Derived::Scalar;
	using  MatX    = Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>;
	using  Vec3	   = Eigen::Vector3<Scalar>;
	using  RowVec3 = Eigen::RowVector3<Scalar>;
	int    n = poses.size();
	MatX   A_pad, b;
	vector<Scalar> t = Linspace(static_cast<Scalar>(0.0), static_cast<Scalar>(1.0), n);
	
	MatX Ms = MatX::Zero(n, 3);
		
	GaussSeidelSolver<Scalar> solver;

	for (int i = 0; i < 3; ++i)
	{
		/* Expandding a Padding to Fit Fomulation */
		A_pad = MatX::Zero(n - 1, n - 1);
		b     = MatX::Zero(n - 2, 1);
		Scalar h_last = t[1] - t[0];
		Scalar b_last = 6 / h_last * (poses[1][i] - poses[0][i]);

		b(0, 0) -= M0 * h_last;
		for (int j = 0; j < n - 2; ++j)
		{
			Scalar h_next = t[j + 2] - t[j + 1];
			Scalar b_next = 6 / h_next * (poses[j + 2][i] - poses[j + 1][i]);

			A_pad(j, j) = 2 * (h_last + h_next);
			A_pad(j, j + 1) = A_pad(j + 1, j) = h_next;
			b(j, 0) = b_next - b_last;

			b_last = b_next;
			h_last = h_next;
		}	
		b(n - 3, 0) -= Mn * h_last;

		// resize to normal size 
		const MatX & A = A_pad.block(0, 0, n - 2, n - 2);
		Ms.block(1, i, n - 2, 1) = solver(A, b);
	}
	
	Ms.block(0, 0, 1, 3)     = RowVec3(M0, M0, M0);
	Ms.block(n - 1, 0, 1, 3) = RowVec3(Mn, Mn, Mn);

	auto si_func = []
	(const pair<Scalar, Scalar>& M, const pair<Scalar, Scalar>& x_delta,
	 const pair<Scalar, Scalar>& y, Scalar hi)->Scalar {
		return (
			M.first  / (6. * hi) * pow(x_delta.second, 3) +
			M.second / (6. * hi) * pow(x_delta.first, 3)  +
			(y.second / hi - M.second * hi / 6.) * x_delta.first +
			(y.first / hi - M.first * hi / 6.) * x_delta.second);
	};

	auto cubic_spline_func = [x = t, y= poses, Ms = Ms, func = si_func]
	(Scalar t)->Vec3{
		static Scalar step = 1. / (x.size() - 1);	
		int  idx = t >= 1 ?x.size() - 2 : t / step;
		Vec3 point;
		const Scalar hi			= x[idx + 1] - x[idx],
					 delta_next = x[idx + 1] - t,
					 delta_cur  = t - x[idx];
		for (int i = 0; i < 3; ++i)
		{
			point[i] = func(
				{ Ms(idx, i), Ms(idx + 1, i) }, 
				{ delta_cur, delta_next }, 
				{ y[idx][i], y[idx + 1][i] }, 
				hi);
		}
		return point;
	};

	return cubic_spline_func;
}

template<Vec3Convertible Derived>
auto
GetCubicSplineFunction(const Eigen::MatrixBase<Derived>& p_ini, 
					   const Eigen::MatrixBase<Derived>& p_end, 
					   const vector<Eigen::Vector3<typename Derived::Scalar>>& pList, 
					   double M0 = 0, double Mn = 0)
{
	using Vec3 = Eigen::Vector3<typename Derived::Scalar>;
	vector<Vec3> temp{p_ini, p_end};
	temp.insert(temp.end() - 1, pList.begin(), pList.end());	
	return GetCubicSplineFunction(temp, M0, Mn);
}

template<Vec6Convertible Derived>
function<Twist<typename Derived::Scalar>(typename Derived::Scalar)>
GetCubicSplineFunction(const Eigen::MatrixBase<Derived>& t_ini,
					   const Eigen::MatrixBase<Derived>& t_end,
					   const vector<Eigen::Vector3<typename Derived::Scalar>>& pList, 
					   double M0 = 0, double Mn = 0)
{
	using Scalar = typename Derived::Scalar;
	using SE3    = Eigen::Matrix4<Scalar>;
	using Vec3   = Eigen::Vector3<Scalar>;
	using Twist  = Eigen::Vector<Scalar, 6>;

    SE3  T_ini = ExpMapping(t_ini),		  T_end = ExpMapping(t_end);
    Vec3 p_ini = T_ini.block(0, 3, 3, 1), p_end = T_end.block(0, 3, 3, 1);
    Vec3 w_ini = t_ini.block(0, 0, 3, 1), w_end = t_end.block(0, 0, 3, 1);

	auto pos_func		   = GetCubicSplineFunction(p_ini, p_end, pList, M0, Mn);
	auto cubic_spline_func = [pos_func= pos_func, w_ini = w_ini, w_end = w_end]
	(Scalar t)->Twist {
        SE3 T_cur = SE3::Identity();
        T_cur.block(0, 3, 3, 1) = pos_func(t);
        T_cur.block(0, 0, 3, 3) = Roderigues(Lerp(w_ini, w_end, t));
        return LogMapSE3Tose3(T_cur);
	};

	return cubic_spline_func;
}

template<Vec3Convertible Derived>
Eigen::Vector3<typename Derived::Scalar>
DeCasteljau(const vector<Derived>& pRest, double t)
{
	using Scalar = typename Derived::Scalar;
	using Vec3   = Eigen::Vector3<Scalar>;

	if (pRest.size() == 1)
	{
		return pRest.back();
	}
	vector<Vec3> pRest_new(pRest.size() - 1);
	std::transform(pRest.begin(), pRest.end() - 1, pRest.begin() + 1, pRest_new.begin(),
		[t = t](auto& num1, auto& num2)->Vec3{
		return Lerp(num1, num2, t);
	});

	return DeCasteljau(pRest_new, t);
}

template<class Scalar>
function<Eigen::Vector3<Scalar>(Scalar)>
GetBezierSplineFunction(const vector<Eigen::Vector3<Scalar>> & pList)
{
	return [pRest = pList](Scalar t)->Vec3d {
		return DeCasteljau(pRest, t);
	};
}

template<class Scalar>
function<Eigen::Vector3<Scalar>(Scalar)>
GetBezierInterSplineFunction(const vector<Eigen::Vector3<Scalar>>& pList, double insertLength = 0.18)
{
	using Vec3 = Eigen::Vector3<Scalar>;

	const int n = pList.size();
	vector<Vec3> control_points(3 * n - 2);

	std::for_each(control_points.begin(), control_points.begin() + 2, [val = pList.front()](auto & num) {num = val; });
	std::for_each(control_points.end() - 2, control_points.end(), [val = pList.back()](auto & num) {num = val; });
	for (int i = 1; i < n - 1; ++i)
	{
		Vec3 delta = pList[i + 1] - pList[i - 1];
		control_points[3 * i]     = pList[i];
		control_points[3 * i - 1] = pList[i] - insertLength * delta;
		control_points[3 * i + 1] = pList[i] + insertLength * delta;
	}

	vector<function<Vec3(Scalar)>> spline_funcs(n - 1);
	{
		auto itf = spline_funcs.begin();
		auto itP = control_points.begin();
		for (; itf != spline_funcs.end(); ++itf, itP += 3)
		{
			*itf = GetBezierSplineFunction(
				vector<Vec3>(itP, itP + 4)
			);
		}
	}
    
	auto bezier_splin_func = [splin_funcs = spline_funcs]
	(Scalar t)->Vec3 {
		static Scalar step = 1. / splin_funcs.size();
		int idx = t > 1 ? splin_funcs.size() - 1 : t / step;
		Scalar t_input = t / step - idx;
		return splin_funcs[idx](t_input);
	};

	return bezier_splin_func;
}

template<class Scalar>
function<Eigen::Vector3<Scalar>(Scalar)>
GetBezierInterSplineFunction(const Eigen::Vector3<Scalar>& p_ini, 
							 const Eigen::Vector3<Scalar>& p_end, 
							 const Eigen::Vector3<Scalar>& pList, 
							 double insertLength = 0.18)
{
	using Vec3 = Eigen::Vector3<Scalar>;
	vector<Vec3> temp{ p_ini, p_end };
	temp.insert(temp.begin() + 1, pList.begin(), pList.end());
	return GetBezierInterSplineFunction(temp, insertLength);
}

template<Vec6Convertible T, Vec3Convertible U>
requires ScalarSame<T, U>
function<Twist<typename T::Scalar>(typename T::Scalar)>
GetBezierInterSplineFunction(const Eigen::MatrixBase<T>& t_ini, 
							 const Eigen::MatrixBase<T>& t_end, 
							 const vector<U>& pList, 
							 double insertLength = 0.18)
{
	using Scalar = typename T::Scalar;
	using Vec3   = Eigen::Vector3<Scalar>;
	using Twist  = Eigen::Vector<Scalar, 6>;
	using SE3    = Eigen::Matrix4<Scalar>;

	Vec3 p_ini = t_ini.block(3, 0, 3, 1), p_end = t_end.block(3, 0, 3, 1);
	Vec3 w_ini = t_ini.block(0, 0, 3, 1), w_end = t_end.block(0, 0, 3, 1);

	auto pos_func = GetBezierInterSplineFunction(p_ini, p_end, pList, insertLength);
	auto bezier_splin_func = [pos_func = pos_func, w_ini = w_ini, w_end = w_end]
	(Scalar t)->Twist {
		SE3 T_cur = SE3::Identity();
		T_cur.block(0, 3, 3, 1) = pos_func(t);
		T_cur.block(0, 0, 3, 3) = Roderigues(Lerp(w_ini, w_end, t));
		return LogMapSE3Tose3(T_cur);
	};
	return bezier_splin_func;
}


enum class BSplineNodeDefinition { Uniform, QuasiUniform, NonUniform };
template<class _AnyVec>
class BSpline :function<_AnyVec(double)>
{
public:
	BSpline(const vector<_AnyVec>& pList, uint32_t order = 3, bool IsInterpolation = false, BSplineNodeDefinition mode = BSplineNodeDefinition::Uniform);
	_AnyVec operator()(double t);

private:
	void Elevation();
	void CalculateControlPoints(const vector<_AnyVec>& pList);

	int curOrder = 0;
	vector<vector<function<double(double)>>> _BaseFunc;
	vector<double>  _uList;
	vector<_AnyVec> _PList;

};

template<class _AnyVec>
BSpline<_AnyVec>::BSpline(const vector<_AnyVec>& pList, uint32_t order, bool isInterpolation, BSplineNodeDefinition mode) :
	_PList(pList), _BaseFunc(order + 1)
{

	switch (const unsigned
		N = pList.size() + order + 1,
		PadNum = order,
		RegionNum = N - 2 * PadNum - 1
		; mode)
	{
    case BSplineNodeDefinition::Uniform: {
		const double BorderPadding = 1.0 / RegionNum * PadNum;
		_uList = Linspace(-BorderPadding, 1.0 + BorderPadding, N);
		break;
	}
    case BSplineNodeDefinition::QuasiUniform: {
		auto temp = Linspace(0.0, 1.0, RegionNum + 1);
		_uList.insert(_uList.cbegin(), PadNum, 0.0);
		_uList.insert(_uList.cbegin() + PadNum, PadNum, 1.0);
		_uList.insert(_uList.cbegin() + PadNum, temp.begin(), temp.end());
		break;
	}
    case BSplineNodeDefinition::NonUniform:

		break;
	}


	while (curOrder <= order)
	{
		Elevation();
	}


	if (isInterpolation)
	{
		if (mode != BSplineNodeDefinition::QuasiUniform)
		{
			assert(mode == BSplineNodeDefinition::QuasiUniform);
			throw("Interpolation Only using in QuasiUniform");
		}
		CalculateControlPoints(pList);
	}
}


template<class _AnyVec>
_AnyVec BSpline<_AnyVec>::operator()(double t) {
	_AnyVec val = _AnyVec::Zero();
	auto PIter = _PList.begin();
	auto FIter = _BaseFunc[curOrder - 1].begin();

	for (; PIter != _PList.end(); ++PIter, ++FIter)
	{
		val += (*FIter)(t) * (*PIter);
	}
	return val;

	/* STL Version */
	/*_AnyVec ZeroVal = _AnyVec::Zero();
	return std::inner_product(_PList.begin(), _PList.end(), _BaseFunc[curOrder - 1].begin(), ZeroVal,
		[](auto Val1, auto addVal)->_AnyVec {
			return Val1 + addVal;
		},
		[t](auto & point, auto& Fun)->_AnyVec{
			return Fun(t) * point;
		});*/
		/* No Test But May get lower efficiency Cause sum function */
}

template<class _AnyVec>
void
BSpline<_AnyVec>::Elevation()
{
	auto& curFuncs = _BaseFunc[curOrder];
	const int N = _uList.size() - curOrder - 1;
	curFuncs.resize(N);
	if (curOrder == 0)
	{
		for (int i = 0; i < N; ++i)
		{
			curFuncs[i] =
				[this, idx = i](double t) {
				return _uList[idx] <= t && t < _uList[idx + 1] ? 1 : 0;
			};
		}
	}
	else
	{
		for (int i = 0; i < N; ++i)
		{
			curFuncs[i] =
				[this, k = curOrder, idx = i](double t) {
				double prev = _uList[idx + k] - _uList[idx];
				double next = _uList[idx + k + 1] - _uList[idx + 1];
				prev += prev == 0 ? 1 : 0;
				next += next == 0 ? 1 : 0;
				return (t - _uList[idx]) / prev * _BaseFunc[k - 1][idx](t) + (_uList[idx + k + 1] - t) / next * _BaseFunc[k - 1][idx + 1](t);
			};
		}
	}
	++curOrder;
}

template<class _AnyVec>
void
BSpline<_AnyVec>::CalculateControlPoints(const vector<_AnyVec>& pList_int)
{
	const int k = curOrder - 1;
	auto& BaseFuncs = _BaseFunc[k];
	const int N = pList_int.size();
	const int DIM = pList_int.front().rows();
	const int M = BaseFuncs.size();
	auto InterList = Linspace(0, 1, pList_int.size());
	DynMatrixd GetThis;
	GetThis = DynMatrixd::Zero(N - 2, M - 2);
	DynMatrixd Residual, catBorder;
	Residual = DynMatrixd::Zero(N - 2, 2);
	catBorder = DynMatrixd::Zero(2, 1);

	for (int i = 0; i < GetThis.rows(); ++i)
	{
		double t_temp = InterList[i + 1];
		for (int j = 0; j < GetThis.cols(); ++j)
		{
			GetThis(i, j) = BaseFuncs[j + 1](t_temp);
		}
		Residual(i, 0) = BaseFuncs.front()(t_temp);
		Residual(i, 1) = BaseFuncs.back()(t_temp);
	}

	GaussSeidelSolver solve;
	vector<_AnyVec> controlPoints(N - 2);
	for (int dim = 0; dim < DIM; ++dim)
	{
		catBorder(0, 0) = pList_int.front()(dim, 0);
		catBorder(1, 0) = pList_int.back()(dim, 0);

		DynMatrixd b = DynMatrixd::Zero(N - 2, 1);
		std::transform(pList_int.begin() + 1, pList_int.end() - 1, b.data(), [dim](auto& a) {return a(dim, 0); });
		b -= Residual * catBorder;

		auto x = solve(GetThis, b);
		for (int i = 0; i < x.rows(); ++i)
		{
			controlPoints[i](dim, 0) = x(i, 0);
		}
	}
	controlPoints.insert(controlPoints.begin(), pList_int.front());
	controlPoints.insert(controlPoints.end(), pList_int.back());
	_PList = controlPoints;
}

inline function<Twistd(double t)>
GetDecompositionFunction(const Twistd & t_ini, const Twistd & t_end)
{
    const SE3d  T_ini = ExpMapping(t_ini),
                T_end = ExpMapping(t_end);
    const Vec3d & w_ini = t_ini.block(0, 0, 3, 1),
                & v_ini = T_ini.block(0, 3, 3, 1),
                & w_end = t_end.block(0, 0, 3, 1),
                & v_end = T_end.block(0, 3, 3, 1);
	Matrix3d R_ini = Roderigues(w_ini), R_end = Roderigues(w_end);
	auto LineFun = [R_ini = R_ini, v_ini = v_ini,
			        R_end = R_end, v_end = v_end]
	(double t)->Twistd {
		Twistd twist;
        twist.block(0, 0, 3, 1) = LogMapSO3Toso3((R_ini * Roderigues(
														 (LogMapSO3Toso3(
														 (R_ini.transpose() * R_end).eval()) * t).eval())).eval());
		twist.block(3, 0, 3, 1) = Lerp(v_ini, v_end, t);
		return twist;
	};
	return LineFun;
}

}

#endif // !_GGEOMETRY_HPP
