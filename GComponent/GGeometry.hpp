//#define _CONSOLE_DEBUG
#include <eigen3/Eigen/Dense>
#include <functional>
#include <algorithm>
#include <numeric>
#include <vector>

#include "GTransform.hpp"
#include "../LSSolver/LinearSystemSolver.hpp"

namespace GComponent
{
using std::function;
using std::tuple;
using std::vector;
using DynMatrixd = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>;
using twistd = Eigen::Matrix<double, 6, 1>;
using vec3d = Eigen::Vector3d;
using vec3Td = Eigen::Matrix<double, 1, 3>;

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

vec3d GetCenterOfCircle(const vec3d & p_1, const vec3d & p_2, const vec3d & p_3)
{
	Matrix3d A;	vec3d x,b;
	
	vec3Td vec_1 = (p_2 - p_1).transpose(),
		   vec_2 = (p_3 - p_1).transpose(),
		   vec_cross = vec_1.cross(vec_2);
	A.block(0, 0, 1, 3) = vec_1;
	A.block(1, 0, 1, 3) = vec_2;
	A.block(2, 0, 1, 3) = vec_cross;
	b = vec3d(vec_1.dot((p_1 + p_2) / 2.0f),
		      vec_2.dot((p_3 + p_1) / 2.0f),
			  vec_cross.dot(p_3));
	
	x = A.inverse() * b;
	return x;
}

function<twistd(double)> GetScrewLineFunction(const twistd & t_ini, const twistd & t_end)
{
	return [t_ini = t_ini, t_end = t_end](double t)->twistd {
		return Lerp(t_ini, t_end, t); };
}

function<twistd(double)> GetScrewLineFunction(const SE3d& T_ini, const SE3d& T_end)
{
	return GetScrewLineFunction(LogMapSE3Tose3(T_ini), LogMapSE3Tose3(T_end));
}

function<vec3d(double)> GetCircleFunction(const vec3d& p_ini, const vec3d& p_end, const vec3d& p_mid)
{
	vec3d rotateAxis = ((p_mid - p_ini).cross(p_end - p_ini)).normalized(),
		  coC = GetCenterOfCircle(p_ini, p_mid, p_end);
	double radius = (p_ini - coC).norm();

	vec3d vec_ini = (p_ini - coC).normalized(),
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
		vec_ini = vec_ini](double t)->vec3d {
		return coCircle + radius * (Roderigues(rotateAxis, t * total) * vec_ini);
	};
	return CircleFun;
}

function<twistd(double)> GetCircleFunction(const twistd& t_ini, const twistd& t_end, const vec3d& p_mid)
{
	const vec3d& p_ini = t_ini.block(3, 0, 3, 1), 
		         p_end = t_end.block(3, 0, 3, 1);
	const vec3d& w_ini = t_ini.block(0, 0, 3, 1),
	             w_end = t_end.block(0, 0, 3, 1);
	auto CirclePosFunc = GetCircleFunction(p_ini, p_end, p_mid);
	auto CircleFun = [PosFunc = CirclePosFunc,
		              w_ini = w_ini, w_end = w_end](double t)->twistd {
		twistd twist;
		twist.block(0, 0, 3, 1) = Lerp(w_ini, w_end, t);
		twist.block(3, 0, 3, 1) = PosFunc(t);
		return twist;
	};
	return CircleFun;
}


function<vec3d(double)> GetCubicSplineFunction(const vector<vec3d>& pList, double M0 = 0, double Mn = 0)
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
		const DynMatrixd & A = A_pad.block(0, 0, n - 2, n - 2);
		Ms.block(1, i, n - 2, 1) = solver(A, b);
	}
	
	Ms.block(0, 0, 1, 3) = vec3Td(M0, M0, M0);
	Ms.block(n - 1, 0, 1, 3) = vec3Td(Mn, Mn, Mn);
	std::cout << "Ms:=\n" << Ms << std::endl << std::endl;


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
	(double t)->vec3d{
		static double step = 1. / (x.size() - 1);	
		int idx = t > 1 ?x.size() - 2 : t / step;
		vec3d point;
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

function<vec3d(double)> GetCubicSplineFunction(const vec3d& p_ini, const vec3d& p_end, const vector<vec3d>& pList, double M0 = 0, double Mn = 0)
{
	vector<vec3d> temp{p_ini, p_end};
	temp.insert(temp.end() - 1, pList.begin(), pList.end());
	return GetCubicSplineFunction(temp, M0, Mn);
}

function<twistd(double)> GetCubicSplineFunction(const twistd& t_ini, const twistd& t_end, const vector<vec3d>& pList, double M0 = 0, double Mn = 0)
{
	const vec3d& p_ini = t_ini.block(3, 0, 3, 1),
				 p_end = t_end.block(3, 0, 3, 1);
	const vec3d& w_ini = t_ini.block(0, 0, 3, 1),
				 w_end = t_end.block(0, 0, 3, 1);
	auto PosFunc = GetCubicSplineFunction(p_ini, p_end, pList, M0, Mn);
	auto SplineFunc = [PosFunc= PosFunc, 
					   w_ini = w_ini, w_end = w_end](double t)->twistd {
		twistd twist;
		twist.block(0, 0, 3, 1) = Lerp(w_ini, w_end, t);
		twist.block(3, 0, 3, 1) = PosFunc(t);
		return twist;
	};
	return SplineFunc;
}


template<class _AnyVec>
inline _AnyVec DeCasteljau(const vector<_AnyVec>& pRest, double t)
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


function<vec3d(double)> GetBezierSplineFunction(const vector<vec3d> & pList)
{
	return [pRest = pList](double t)->vec3d {
		return DeCasteljau(pRest, t);
	};
}

function<vec3d(double)> GetBezierInterSplineFunction(const vector<vec3d>& pList, double insertLength = 0.18)
{
	const int n = pList.size();
	vector<vec3d> contrlPoints(3 * n - 2);

	std::for_each(contrlPoints.begin(), contrlPoints.begin() + 2, [val = pList.front()](auto & num) {num = val; });
	std::for_each(contrlPoints.end() - 2, contrlPoints.end(), [val = pList.back()](auto & num) {num = val; });
	for (int i = 1; i < n - 1; ++i)
	{
		vec3d delta = pList[i + 1] - pList[i - 1];
		contrlPoints[3 * i] = pList[i];
		contrlPoints[3 * i - 1] = pList[i] - insertLength * delta;
		contrlPoints[3 * i + 1] = pList[i] + insertLength * delta;
	}

	vector<function<vec3d(double)>> SplineFuncs(n - 1);
	{
		auto itf = SplineFuncs.begin();
		auto itP = contrlPoints.begin();
		for (; itf != SplineFuncs.end(); ++itf, itP += 3)
		{
			*itf = GetBezierSplineFunction(
				vector<vec3d>(itP, itP + 4)
			);
		}
	}
    
	auto BezierSplineFunc = [SplineFuncs = SplineFuncs]
	(double t)->vec3d {
		static double step = 1. / SplineFuncs.size();

		int idx = t > 1 ? SplineFuncs.size() - 1 : t / step;

		double t_input = t / step - idx;
		return SplineFuncs[idx](t_input);
	};

	return BezierSplineFunc;
}

function<vec3d(double)> GetBezierInterSplineFunction(const vec3d& p_ini, const vec3d& p_end, const vector<vec3d>& pList, double insertLength = 0.18)
{
	vector<vec3d> temp{ p_ini, p_end };
	temp.insert(temp.begin() + 1, pList.begin(), pList.end());
	return GetBezierInterSplineFunction(temp, insertLength);
}

function<twistd(double)> GetBezierInterSplineFunction(const twistd& t_ini, const twistd& t_end, const vector<vec3d>& pList, double insertLength = 0.18)
{
	const vec3d& p_ini = t_ini.block(3, 0, 3, 1),
				 p_end = t_end.block(3, 0, 3, 1);
	const vec3d& w_ini = t_ini.block(0, 0, 3, 1),
				 w_end = t_end.block(0, 0, 3, 1);
	auto PosFunc = GetBezierInterSplineFunction(p_ini, p_end, pList, insertLength);
	auto SplineFunc = [PosFunc = PosFunc,
		w_ini = w_ini, w_end = w_end](double t)->twistd {
		twistd twist;
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
	vector<double> _uList;
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
		using enum BSplineNodeDefinition;
	case Uniform: {
		const double BorderPadding = 1.0 / RegionNum * PadNum;
		_uList = Linspace(-BorderPadding, 1.0 + BorderPadding, N);
		break;
	}
	case QuasiUniform: {
		auto temp = Linspace(0.0, 1.0, RegionNum + 1);
		_uList.insert(_uList.cbegin(), PadNum, 0.0);
		_uList.insert(_uList.cbegin() + PadNum, PadNum, 1.0);
		_uList.insert(_uList.cbegin() + PadNum, temp.begin(), temp.end());
		break;
	}
	case NonUniform:

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
void BSpline<_AnyVec>::Elevation()
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
void BSpline<_AnyVec>::CalculateControlPoints(const vector<_AnyVec>& pList_int)
{
	const int k = curOrder - 1;
	auto& BaseFuncs = _BaseFunc[k];
	const int N = pList_int.size();
	const int DIM = pList_int.front().rows();
	const int M = BaseFuncs.size();
	auto InterList = Linspace(0, 1, pList_int.size());
	DynMatrixd A;
	A = DynMatrixd::Zero(N - 2, M - 2);
	DynMatrixd Residual, catBorder;
	Residual = DynMatrixd::Zero(N - 2, 2);
	catBorder = DynMatrixd::Zero(2, 1);

	for (int i = 0; i < A.rows(); ++i)
	{
		double t_temp = InterList[i + 1];
		for (int j = 0; j < A.cols(); ++j)
		{
			A(i, j) = BaseFuncs[j + 1](t_temp);
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

		auto x = solve(A, b);
		for (int i = 0; i < x.rows(); ++i)
		{
			controlPoints[i](dim, 0) = x(i, 0);
		}
	}
	controlPoints.insert(controlPoints.begin(), pList_int.front());
	controlPoints.insert(controlPoints.end(), pList_int.back());
	_PList = controlPoints;
}

function<twistd(double t)> GetDecompositionFunction(const twistd & t_ini, const twistd & t_end)
{
	const vec3d& w_ini = t_ini.block(0, 0, 3, 1),
	     		& v_ini = t_ini.block(3, 0, 3, 1),
			& w_end = t_end.block(0, 0, 3, 1),
			& v_end = t_end.block(3, 0, 3, 1);
	Matrix3d R_ini = Roderigues(w_ini), R_end = Roderigues(w_end);
	auto LineFun = [R_ini = R_ini, v_ini = v_ini,
			        R_end = R_end, v_end = v_end]
	(double t)->twistd {
		twistd twist;
		twist.block(0, 0, 3, 1) = LogMapSO3Toso3(R_ini * Roderigues(LogMapSO3Toso3(R_ini.transpose() * R_end) * t));
		twist.block(3, 0, 3, 1) = Lerp(v_ini, v_end, t);
		return twist;
	};
	return LineFun;
}


}
