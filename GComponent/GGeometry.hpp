#ifndef _GGEOMETRY_HPP
#define _GGEOMETRY_HPP

#include <eigen3/Eigen/Dense>
#include <functional>
#include <algorithm>
#include <numeric>
#include <vector>

#include <GComponent/GTransform.hpp>

#include <LSSolver/LinearSystemSolver.hpp>

namespace GComponent
{
using std::function;
using std::tuple;
using std::vector;
using DynMatrixd = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>;
using Twistd = Eigen::Matrix<double, 6, 1>;
using Vec3d = Eigen::Vector3d;
using Vec3dT = Eigen::Matrix<double, 1, 3>;
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
inline _AnyVec Lerp(const _AnyVec & v1, const _AnyVec& v2, double t)
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

inline double Lerp(const double & d1, const double& d2, double t)
{
    return (1 - t) * d1 + t * d2;
}

/// <summary>
/// 将线段空间进行 n 等分函数
/// </summary>
/// <param name="lowerBound"></param>
/// <param name="upperBound"></param>
/// <param name="n"></param>
/// <returns></returns>
inline vector<double> Linspace(const double lowerBound, const double upperBound,const int n)
{
	vector<double> rval(n);
	const double step = (upperBound - lowerBound) / (n - 1);
	double val_cur = lowerBound;
	for (auto it = rval.begin(); it != rval.end(); ++it)
	{
		*it = val_cur;
		val_cur += step;
	}
	return rval;
}

/// <summary>
/// 空间三点获取圆心函数
/// </summary>
/// <param name="p_1"></param>
/// <param name="p_2"></param>
/// <param name="p_3"></param>
/// <returns></returns>
inline Vec3d GetCenterOfCircle(const Vec3d & p_1, const Vec3d & p_2, const Vec3d & p_3)
{
	Matrix3d GetThis;	Vec3d x,b;
	
	Vec3dT vec_1 = (p_2 - p_1).transpose(),
		   vec_2 = (p_3 - p_1).transpose(),
		   vec_cross = vec_1.cross(vec_2);
	GetThis.block(0, 0, 1, 3) = vec_1;
	GetThis.block(1, 0, 1, 3) = vec_2;
	GetThis.block(2, 0, 1, 3) = vec_cross;
	b = Vec3d(vec_1.dot((p_1 + p_2) / 2.0f),
		      vec_2.dot((p_3 + p_1) / 2.0f),
			  vec_cross.dot(p_3));
	
	x = GetThis.inverse() * b;
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


// FIXME: 与GetCircleFunction冗余
/// <summary>
/// 空间三点获取圆弧函数
/// </summary>
/// <param name="p_1"></param>
/// <param name="p_2"></param>
/// <param name="p_3"></param>
/// <returns></returns>
inline double GetRadiusOfCircle(const Vec3d & p_1, const Vec3d& p_2, const Vec3d &p_3)
{
    Vec3d coc = GetCenterOfCircle(p_1, p_2, p_3);
    return ( p_1 - coc).norm();
}

// FIXME: 与GetCircleFunction冗余
/// <summary>
/// 空间三点获取圆弧差角函数
/// </summary>
/// <param name="p_1"></param>
/// <param name="p_2"></param>
/// <param name="p_3"></param>
/// <returns></returns>
inline double GetArcDeltaOfCircle(const Vec3d & p_1, const Vec3d& p_2, const Vec3d &p_3)
{
    Vec3d rotateAxis = ((p_1 - p_2).cross(p_3 - p_1)).normalized();
    Vec3d coc = GetCenterOfCircle(p_1, p_2, p_3);
    Vec3d vec_ini = (p_1 - coc).normalized(),
          vec_end = (p_3 - coc).normalized();

    double angle =  acos(vec_ini.dot(vec_end));

    return rotateAxis.dot(vec_ini.cross(vec_end)) > 0.?
                angle:
                2* EIGEN_PI - angle;
}

/// <summary>
/// 螺旋轴直线插值函数
/// </summary>
/// <param name="t_ini"></param>
/// <param name="t_end"></param>
/// <returns></returns>
inline function<Twistd(double)> GetScrewLineFunction(const SE3d& T_ini, const SE3d& T_end)
{
    const Vec3d
        & pos_ini = T_ini.block(0, 3, 3, 1),
        & pos_end = T_end.block(0, 3, 3, 1);
    Twistd
          t_ini   = LogMapSE3Tose3(T_ini),
          t_end   = LogMapSE3Tose3(T_end);
    const Vec3d
        & w_ini   = t_ini.block(0, 0, 3, 1),
        & w_end   = t_end.block(0, 0, 3, 1);

    auto LineFunction = [pos_ini = pos_ini, pos_end = pos_end, w_ini = w_ini, w_end = w_end](double t){
        SE3d  T_cur   = SE3d::Identity();
        Vec3d pos_cur = Lerp(pos_ini, pos_end, t);
        Vec3d w_cur   = Lerp(w_ini, w_end, t);
        T_cur.block(0, 0, 3, 3) = Roderigues(w_cur);
        T_cur.block(0, 3, 3, 1) = pos_cur;
        return LogMapSE3Tose3(T_cur);
    };

    return LineFunction;
}

inline function<Twistd(double)> GetScrewLineFunction(const Twistd & t_ini, const Twistd & t_end)
{
    return GetScrewLineFunction(ExpMapping(t_ini), ExpMapping(t_end));
}


/// <summary>
/// 三点画圆弧函数
/// </summary>
/// <param name="p_ini"></param>
/// <param name="p_end"></param>
/// <param name="p_mid"></param>
/// <returns></returns>
inline function<Vec3d(double)> GetCircleFunction(const Vec3d& p_ini, const Vec3d& p_end, const Vec3d& p_mid)
{
	Vec3d rotateAxis = ((p_mid - p_ini).cross(p_end - p_ini)).normalized(),
		  coC = GetCenterOfCircle(p_ini, p_mid, p_end);
	double radius = (p_ini - coC).norm();

	Vec3d vec_ini = (p_ini - coC).normalized(),
		  vec_end = (p_end - coC).normalized();

	double theta_delta =
		rotateAxis.dot(vec_ini.cross(vec_end)) > 0.?
		acos(vec_ini.dot(vec_end)) :
		-acos(vec_ini.dot(vec_end)) + 2 * EIGEN_PI;

	auto CircleFun = [
		total = theta_delta,
		rotateAxis = rotateAxis,
		coCircle = coC,
		radius = radius,
		vec_ini = vec_ini](double t)->Vec3d {
		return coCircle + radius * (Roderigues(rotateAxis, t * total) * vec_ini);
	};
	return CircleFun;
}

inline
function<Twistd(double)>
GetCircleFunction(const Twistd& t_ini, const Twistd& t_end, const Vec3d& p_mid)
{
    SE3d T_ini = ExpMapping(t_ini),
         T_end = ExpMapping(t_end);
    const Vec3d& p_ini = T_ini.block(0, 3, 3, 1),
                 p_end = T_end.block(0, 3, 3, 1);
	const Vec3d& w_ini = t_ini.block(0, 0, 3, 1),
	             w_end = t_end.block(0, 0, 3, 1);
	auto CirclePosFunc = GetCircleFunction(p_ini, p_end, p_mid);
	auto CircleFun = [PosFunc = CirclePosFunc,
		              w_ini = w_ini, w_end = w_end](double t)->Twistd {
        SE3d T_cur = SE3d::Identity();
        T_cur.block(0, 3, 3, 1) = PosFunc(t);
        T_cur.block(0, 0, 3, 3) = Roderigues(Lerp(w_ini, w_end, t));
        return LogMapSE3Tose3(T_cur);
	};
	return CircleFun;
}

inline
function<Vec3d(double)>
GetCubicSplineFunction(const vector<Vec3d>& pList, double M0 = 0, double Mn = 0)
{
	int n = pList.size();
	DynMatrixd A_pad, b;
	vector<double> t = Linspace(0, 1, n);
	
	DynMatrixd Ms = DynMatrixd::Zero(n, 3);
		
	GaussSeidelSolver solver;

	for (int i = 0; i < 3; ++i)
	{
		/* Expandding a Padding to Fit Fomulation */
		A_pad = DynMatrixd::Zero(n - 1, n - 1);
		b = DynMatrixd::Zero(n - 2, 1);
		double h_last = t[1] - t[0];
		double b_last = 6 / h_last * (pList[1][i] - pList[0][i]);

		b(0, 0) -= M0 * h_last;
		for (int j = 0; j < n - 2; ++j)
		{
			double h_next = t[j + 2] - t[j + 1];
			double b_next = 6 / h_next * (pList[j + 2][i] - pList[j + 1][i]);

			A_pad(j, j) = 2 * (h_last + h_next);
			A_pad(j, j + 1) = A_pad(j + 1, j) = h_next;
			b(j, 0) = b_next - b_last;

			b_last = b_next;
			h_last = h_next;
		}	
		b(n - 3, 0) -= Mn * h_last;

		/* resize to The Normal Size */
		const DynMatrixd & GetThis = A_pad.block(0, 0, n - 2, n - 2);
		Ms.block(1, i, n - 2, 1) = solver(GetThis, b);
	}
	
	Ms.block(0, 0, 1, 3) = Vec3dT(M0, M0, M0);
	Ms.block(n - 1, 0, 1, 3) = Vec3dT(Mn, Mn, Mn);
    //std::cout << "Ms:=\n" << Ms << std::endl << std::endl;

	auto SiFun = []
	(const pair<double, double>& M,
	 const pair<double, double>& x_delta,
	 const pair<double, double>& y,
	 double hi)->double {
		return (
			M.first  / (6. * hi) * pow(x_delta.second, 3) +
			M.second / (6. * hi) * pow(x_delta.first, 3)  +
			(y.second / hi - M.second * hi / 6.) * x_delta.first +
			(y.first / hi - M.first * hi / 6.) * x_delta.second);
	};

	auto SlineFun = [x = t, y= pList, Ms = Ms, Fun = SiFun]
	(double t)->Vec3d{
		static double step = 1. / (x.size() - 1);	
		int idx = t >= 1 ?x.size() - 2 : t / step;
		Vec3d point;
		const double hi = x[idx + 1] - x[idx],
					 delta_next = x[idx + 1] - t,
					delta_cur  = t - x[idx];
		for (int i = 0; i < 3; ++i)
		{
			point[i] = Fun(
				{ Ms(idx, i), Ms(idx + 1, i) }, 
				{ delta_cur, delta_next }, 
				{ y[idx][i], y[idx + 1][i] }, 
				hi);
		}
		return point;
	};

	return SlineFun;
}

inline function<Vec3d(double)>
GetCubicSplineFunction(const Vec3d& p_ini, const Vec3d& p_end, const vector<Vec3d>& pList, double M0 = 0, double Mn = 0)
{
	vector<Vec3d> temp{p_ini, p_end};
	temp.insert(temp.end() - 1, pList.begin(), pList.end());
	return GetCubicSplineFunction(temp, M0, Mn);
}

inline function<Twistd(double)>
GetCubicSplineFunction(const Twistd& t_ini, const Twistd& t_end, const vector<Vec3d>& pList, double M0 = 0, double Mn = 0)
{
    SE3d T_ini = ExpMapping(t_ini),
         T_end = ExpMapping(t_end);
    const Vec3d& p_ini = T_ini.block(0, 3, 3, 1),
                 p_end = T_end.block(0, 3, 3, 1);
    const Vec3d& w_ini = t_ini.block(0, 0, 3, 1),
                 w_end = t_end.block(0, 0, 3, 1);
	auto PosFunc = GetCubicSplineFunction(p_ini, p_end, pList, M0, Mn);
	auto SplineFunc = [PosFunc= PosFunc, 
					   w_ini = w_ini, w_end = w_end](double t)->Twistd {
        SE3d T_cur = SE3d::Identity();
        T_cur.block(0, 3, 3, 1) = PosFunc(t);
        T_cur.block(0, 0, 3, 3) = Roderigues(Lerp(w_ini, w_end, t));
        return LogMapSE3Tose3(T_cur);
	};
	return SplineFunc;
}


template<class _AnyVec>
inline _AnyVec
DeCasteljau(const vector<_AnyVec>& pRest, double t)
{
	if (pRest.size() == 1)
	{
		return pRest.back();
	}
	vector<_AnyVec> pRest_new(pRest.size() - 1);
	std::transform(pRest.begin(), pRest.end() - 1, pRest.begin() + 1, pRest_new.begin(),
		[t = t](auto& num1, auto& num2) {
		return Lerp(num1, num2, t);
	});
	return DeCasteljau(pRest_new, t);
}


inline function<Vec3d(double)>
GetBezierSplineFunction(const vector<Vec3d> & pList)
{
	return [pRest = pList](double t)->Vec3d {
		return DeCasteljau(pRest, t);
	};
}

inline function<Vec3d(double)>
GetBezierInterSplineFunction(const vector<Vec3d>& pList, double insertLength = 0.18)
{
	const int n = pList.size();
	vector<Vec3d> contrlPoints(3 * n - 2);

	std::for_each(contrlPoints.begin(), contrlPoints.begin() + 2, [val = pList.front()](auto & num) {num = val; });
	std::for_each(contrlPoints.end() - 2, contrlPoints.end(), [val = pList.back()](auto & num) {num = val; });
	for (int i = 1; i < n - 1; ++i)
	{
		Vec3d delta = pList[i + 1] - pList[i - 1];
		contrlPoints[3 * i] = pList[i];
		contrlPoints[3 * i - 1] = pList[i] - insertLength * delta;
		contrlPoints[3 * i + 1] = pList[i] + insertLength * delta;
	}

	vector<function<Vec3d(double)>> SplineFuncs(n - 1);
	{
		auto itf = SplineFuncs.begin();
		auto itP = contrlPoints.begin();
		for (; itf != SplineFuncs.end(); ++itf, itP += 3)
		{
			*itf = GetBezierSplineFunction(
				vector<Vec3d>(itP, itP + 4)
			);
		}
	}
    
	auto BezierSplineFunc = [SplineFuncs = SplineFuncs]
	(double t)->Vec3d {
		static double step = 1. / SplineFuncs.size();

		int idx = t > 1 ? SplineFuncs.size() - 1 : t / step;

		double t_input = t / step - idx;
		return SplineFuncs[idx](t_input);
	};

	return BezierSplineFunc;
}

inline function<Vec3d(double)>
GetBezierInterSplineFunction(const Vec3d& p_ini, const Vec3d& p_end, const vector<Vec3d>& pList, double insertLength = 0.18)
{
	vector<Vec3d> temp{ p_ini, p_end };
	temp.insert(temp.begin() + 1, pList.begin(), pList.end());
	return GetBezierInterSplineFunction(temp, insertLength);
}

inline function<Twistd(double)>
GetBezierInterSplineFunction(const Twistd& t_ini, const Twistd& t_end, const vector<Vec3d>& pList, double insertLength = 0.18)
{
	const Vec3d& p_ini = t_ini.block(3, 0, 3, 1),
				 p_end = t_end.block(3, 0, 3, 1);
	const Vec3d& w_ini = t_ini.block(0, 0, 3, 1),
				 w_end = t_end.block(0, 0, 3, 1);
	auto PosFunc = GetBezierInterSplineFunction(p_ini, p_end, pList, insertLength);
	auto SplineFunc = [PosFunc = PosFunc,
		w_ini = w_ini, w_end = w_end](double t)->Twistd {
		Twistd twist;
		twist.block(0, 0, 3, 1) = Lerp(w_ini, w_end, t);
		twist.block(3, 0, 3, 1) = PosFunc(t);
		return twist;
	};
	return SplineFunc;
}


enum class BSplineNodeDefinition { Uniform, QuasiUniform, NonUniform };
template<class _AnyVec>
class BSpline :function<_AnyVec(double)>
{
public:
	BSpline(const vector<_AnyVec>& pList, unsigned order = 3, bool IsInterpolation = false, BSplineNodeDefinition mode = BSplineNodeDefinition::Uniform);
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
BSpline<_AnyVec>::BSpline(const vector<_AnyVec>& pList, unsigned order, bool isInterpolation, BSplineNodeDefinition mode) :
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
