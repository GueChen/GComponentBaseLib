#ifndef _LINEARSYSTEM_SOLVER_HPP
#define _LINEARSYSTEM_SOLVER_HPP

#include <Eigen/Dense>
#include <Eigen/SVD>

#include <LSSolver/PenroseInverse.hpp>
#include <GComponent/GNumerical.hpp>

#include <ranges>
#include <algorithm>
#include <execution>

using Eigen::Matrix;
using Eigen::Vector;

template<class T, class U>
concept Equivalence = std::is_convertible_v<T, U> | std::is_same_v<T, U>;

template<class _Scaler>
inline Eigen::Vector<_Scaler, Eigen::Dynamic> ClampMaxAbs(const Eigen::Vector<_Scaler, Eigen::Dynamic>& vec, _Scaler factor)
{
	return vec.norm() > factor ? vec.normalized() * factor : vec;
}
template<class _Scaler, class _OtherScaler> requires std::is_convertible_v<_OtherScaler, _Scaler>
inline Eigen::Vector<_Scaler, Eigen::Dynamic> ClampMaxAbs(const Eigen::Vector<_Scaler, Eigen::Dynamic>& vec, _OtherScaler factor) {
	return vec.norm() > factor ? vec.normalized() * factor : vec;
}

template<class _Scaler>
class LinearSystemSolver {
	using _DynMat = Matrix<_Scaler, Eigen::Dynamic, Eigen::Dynamic>;
	using _Vec = Matrix<_Scaler, Eigen::Dynamic, 1>;
public:
	LinearSystemSolver()		  = default;
	virtual ~LinearSystemSolver() = default;
	virtual _DynMat operator()(const _DynMat& A, const _DynMat& b) const = 0;

};

template<class _Scaler>
class DynamicLeastNormSolver : public LinearSystemSolver<_Scaler>{
	using _Tdyn = Matrix<_Scaler, Eigen::Dynamic, Eigen::Dynamic>;

	using SVD   = Eigen::JacobiSVD<_Tdyn>;
public:
	DynamicLeastNormSolver()	= default;
	~DynamicLeastNormSolver()	= default;

	_Tdyn operator()(const _Tdyn& A, const _Tdyn& b) const override{
		SVD svd(A, Eigen::ComputeFullU | Eigen::ComputeFullV);
		return svd.solve(b);
	}
};
template<class _Scaler>
using LNSolver = DynamicLeastNormSolver<_Scaler>;


template<class _Scaler>
class DynamicWeightedLeastNormSolver: public LinearSystemSolver<_Scaler> {
	using _Tdyn = Matrix<_Scaler, Eigen::Dynamic, Eigen::Dynamic>;
	using SVD	= Eigen::JacobiSVD<_Tdyn>;
public:
    explicit DynamicWeightedLeastNormSolver(const _Tdyn& weight_mat) : _weight_mat_inverse(weight_mat.inverse()) {}
	~DynamicWeightedLeastNormSolver() = default;

	_Tdyn operator()(const _Tdyn& A, const _Tdyn& b) const override {
		SVD svd(A * _weight_mat_inverse, Eigen::ComputeFullU | Eigen::ComputeFullV);
		return _weight_mat_inverse * svd.solve(b);
	}
private:
	_Tdyn _weight_mat_inverse;
};
template<class _Scaler>
using WLNSolver = DynamicWeightedLeastNormSolver<_Scaler>;

template<class _Scaler>
class DynamicDampedLeastSquareSolver : public LinearSystemSolver<_Scaler> {
	using _Tdyn = Matrix<_Scaler, Eigen::Dynamic, Eigen::Dynamic>;
public:
	template<class _OtherScaler> requires std::is_convertible_v<_OtherScaler, _Scaler>
	explicit DynamicDampedLeastSquareSolver(_OtherScaler factor) :damped_factor_(std::max(factor, static_cast<_OtherScaler>(MIN_FACTOR))) {}
	explicit DynamicDampedLeastSquareSolver(_Scaler		 factor) :damped_factor_(std::max(factor, MIN_FACTOR)) {}

	~DynamicDampedLeastSquareSolver() = default;

	template<class _OtherScaler> requires std::is_convertible_v<_OtherScaler, _Scaler>
	inline void SetDampedFactor(_OtherScaler factor) { damped_factor_ = std::max(factor, static_cast<_OtherScaler>(MIN_FACTOR)); }
	inline void SetDampedFactor(_Scaler		 factor) { damped_factor_ = std::max(factor, MIN_FACTOR); }

	_Tdyn operator()(const _Tdyn& A, const _Tdyn& b) const override {
		return A.transpose() * (A * A.transpose() + std::pow(damped_factor_, 2) * _Tdyn::Identity(A.rows(), A.rows())).inverse() * b;
	}

private:
	_Scaler damped_factor_ = 1e-5f;

	static constexpr _Scaler MIN_FACTOR = 1e-10;
};
template<class _Scaler>
using DLSSolver = DynamicDampedLeastSquareSolver<_Scaler>;

template<class _Scaler>
class DynamicAdaptiveDampedLeastSquareSolver : public LinearSystemSolver<_Scaler> {
	using _Tdyn		= Matrix<_Scaler, Eigen::Dynamic, Eigen::Dynamic>;
	using _Vecdyn	= Matrix<_Scaler, Eigen::Dynamic, 1>;
	using SVD		= Eigen::JacobiSVD<_Tdyn>;
public:
	DynamicAdaptiveDampedLeastSquareSolver()  = default;	
	~DynamicAdaptiveDampedLeastSquareSolver() = default;

	template<class _OtherScaler> requires std::is_convertible_v<_OtherScaler, _Scaler>
	inline void SetFactor(_OtherScaler	factor) { scaler_ = std::max(factor, static_cast<_OtherScaler>(MIN_FACTOR)); }
	inline void SetFactor(_Scaler		factor) { scaler_ = std::max(factor, MIN_FACTOR); }

	_Tdyn operator()(const _Tdyn& A, const _Tdyn& b) const override {
		assert(b.cols() == 1);
		int n = A.cols();
		SVD svd(A.transpose() * A + std::pow(scaler_ * b.norm(), 2) * _Tdyn::Identity(n, n)
			, Eigen::ComputeThinU | Eigen::ComputeThinV);

		return svd.solve(A.transpose() * b);
	}

private:
	_Scaler scaler_ = 0.001;
	static constexpr _Scaler MIN_FACTOR = 1e-15;
};
template<class _Scaler>
using ADLSSolver = DynamicAdaptiveDampedLeastSquareSolver<_Scaler>;

template<class _Scaler>
class DynamicJacobianTransposeSolver : public LinearSystemSolver<_Scaler> {
	using _Tdyn   = Matrix<_Scaler, Eigen::Dynamic, Eigen::Dynamic>;
	using _Vecdyn = Matrix<_Scaler, Eigen::Dynamic, 1>;

public:
	DynamicJacobianTransposeSolver()  = default;
	~DynamicJacobianTransposeSolver() = default;

	template<class _OtherScaler> requires std::is_convertible_v<_OtherScaler, _Scaler>
	inline void SetFactor(_OtherScaler	factor)	{ scaler_ = std::max(factor, static_cast<_OtherScaler>(MIN_FACTOR)); }
	inline void SetFactor(_Scaler		factor) { scaler_ = std::max(factor, MIN_FACTOR); }

	_Tdyn operator()(const _Tdyn& A, const _Tdyn& b) const override {
		_Vecdyn grad		  = A.transpose() * b;
		_Vecdyn auxiliary_vec = A * grad;
		return scaler_ * static_cast<_Vecdyn>(b).dot(auxiliary_vec) / auxiliary_vec.norm() * grad;
	}

private:
	_Scaler scaler_ = 0.001;
	static constexpr _Scaler MIN_FACTOR = 1e-15;
};
template<class _Scaler>
using JacobiTSolver = DynamicJacobianTransposeSolver<_Scaler>;

template<class _Scaler>
class DynamicSelectivelyDampedLeastSquareSolver :public LinearSystemSolver<_Scaler> {
	using _Tdyn		= Matrix<_Scaler, Eigen::Dynamic, Eigen::Dynamic>;
	using _Vecdyn	= Matrix<_Scaler, Eigen::Dynamic, 1>;
	using SVD		= Eigen::JacobiSVD<_Tdyn>;
public:
	DynamicSelectivelyDampedLeastSquareSolver()  = default;
	~DynamicSelectivelyDampedLeastSquareSolver() = default;

	_Tdyn operator()(const _Tdyn& A, const _Tdyn& b) const override {
		
		assert(b.cols() == 1);
		size_t m = A.rows(), n = A.cols();
		SVD svd(A, Eigen::ComputeThinU | Eigen::ComputeThinV);
		
		_Vecdyn A_norm;
		A_norm.resize(n);
		std::transform(std::execution::par_unseq, A.colwise().begin(), A.colwise().end(), A_norm.begin(), [](auto col) {
			return 1.0 / col.norm();
			});
		size_t  rank = svd.rank();
		_Vecdyn phi  = _Vecdyn::Zero(n);		
		
		for (int i = 0; i < rank; ++i) {
			_Vecdyn v		  = svd.matrixV().col(i),
					u		  = svd.matrixU().col(i);
			_Scaler inv_sigma = 1.0 / svd.singularValues()[i];

			_Scaler alpha	  = u.dot(static_cast<_Vecdyn>(b));
			_Vecdyn phi_tmp	  = inv_sigma * alpha * v;
			
			_Scaler N		  = u.norm();	
			_Scaler M		  = inv_sigma * std::inner_product(
					v.begin(), v.end(), A_norm.begin(), static_cast<_Scaler>(0.0), 
					std::plus<>{}, [](auto& v_val, auto& norm_scaler)->_Scaler {
						return abs(v_val) * norm_scaler;
					}
			);

			_Scaler gamma	  = std::min(_Scaler(1.), N / M);
			phi_tmp			  = ClampMaxAbs(phi_tmp, gamma);
			phi				 += phi_tmp;
		}

		return ClampMaxAbs(phi, GammaMax);
		//return phi;
	}

private:
	static constexpr _Scaler GammaMax = GComponent::MyPI / 4.0;
};
template<class _Scaler>
using SDLSSolver = DynamicSelectivelyDampedLeastSquareSolver<_Scaler>;


template<unsigned Cow, unsigned Row>
class StaticInverseSolver {
	using _TA = Matrix<double, Cow, Row>;
	using _Tb = Matrix<double, Cow, 1>;
	using _Tx = Matrix<double, Row, 1>;
public:
	inline _Tx operator()(const _TA& A, const _Tb& b)
	{
		assert(A.rows() <= A.cols() && "The Matrix Rows can't more then Cols, Equation has No Solution.");
		if(A.rows() == A.cols())
			return A.inverse() * b;					
		PenroseInverseSolver<double> pinv;
		return pinv(A) * b;
	}
};

template<class _Scalar>
class DynamicInverseSolver {
	using _Tdyn = Matrix<_Scalar, Eigen::Dynamic, Eigen::Dynamic>;
	using _Self = DynamicInverseSolver;
public:
	inline _Tdyn operator()(const _Tdyn& A, const _Tdyn& b) 
	{
		assert(A.rows() <= A.cols() && "The Matrix Rows can't more then Cols, Equation has No Solution.");
		if(A.rows() == A.cols())
			return A.inverse() * b;		
		PenroseInverseSolver<_Scalar> pinv;
		return pinv(A) * b;
	}

	DynamicInverseSolver(_Self& other)  = delete;
	DynamicInverseSolver(_Self&& other) = delete;
	_Self& operator=(_Self& other)		= delete;
	_Self& operator=(_Self&& other)		= delete;
};

template<class _Scalar>
class JacobiSolver {
	using _Tdyn = Matrix<_Scalar, Eigen::Dynamic, Eigen::Dynamic>;
public:
	inline _Tdyn operator()(const _Tdyn& A, const _Tdyn& b, int _MaxIteration = 100, double _MaxResidual = 1e-8) const
	{
		_Tdyn LU = A, D = _Tdyn::Zero(A.cols(), A.rows());
		for (int i = 0; i < A.rows(); ++i)
		{
			std::swap(D(i, i),LU(i, i));
			D(i, i) = 1.0 / D(i, i);
		}
		_Tdyn x = _Tdyn::Zero(A.cols(), 1);
		for (int iter = 0; iter < _MaxIteration; ++iter)
		{
			x = D * (b - LU * x);
			_Scalar residual = (A * x - b).lpNorm<2>();
			if (residual < _MaxResidual)
			{
				break;
			}
		}
		return x;
	}
};

template<class _Scalar>
class GaussSeidelSolver {
	using _Tdyn = Matrix<_Scalar, Eigen::Dynamic, Eigen::Dynamic>;
public:
	_Tdyn operator()(const _Tdyn& A, const _Tdyn& b, int _MaxIteration = 100, double _MaxResidual = 1e-8) const
	{
		const unsigned ROW = A.rows(), COL = A.cols();
		_Tdyn LD = _Tdyn::Zero(ROW, COL), U = LD;
        for (unsigned i = 0; i < ROW; ++i)
		{
            for (unsigned j = 0; j < COL; ++j)
			{
				if (j > i)
				{
					U(i, j) = A(i, j);
				}
				else
				{
					LD(i, j) = A(i, j);
				}
			}
		}
		PenroseInverseSolver<_Scalar> pinv;
		auto invLD = pinv(LD);
		_Tdyn x = _Tdyn::Zero(A.cols(), 1);
		for (int iter = 0; iter < _MaxIteration; ++iter)
		{
			x = invLD * (b - U * x);
			_Scalar residual = (A * x - b).lpNorm<2>();
			if (residual < _MaxResidual)
			{
				break;
			}
		}
		return x;
	}
#ifdef _NO_OPT
	inline _Tdyn operator()(const _Tdyn& A, const _Tdyn& b, int _MaxIteration = 100, double _MaxResidual = 1e-8)
	{
		const unsigned ROW = A.rows(), COL = A.cols();
		_Tdyn LU = A, D = Eigen::MatrixXd::Zero(COL, ROW);
		for (int i = 0; i < ROW; ++i)
		{
			std::swap(D(i, i), LU(i, i));
			D(i, i) = 1.0 / D(i, i);
		}
		_Tdyn _x = Eigen::MatrixXd::Zero(A.cols(), 1);
		for (int iter = 0; iter < _MaxIteration; ++iter)
		{
			for (int i = 0; i < ROW; ++i)
			{
				_x(i, 0) = D(i, i) * (b(i, 0) - (LU.block(i, 0, i + 1, ROW) * _x).value());
			}
			double residual = (A * _x - b).lpNorm<2>();
			if (residual < _MaxResidual)
			{
				std::cout << "Gauss-Seidel optimizal iterations:" << iter << std::endl;
				break;
			}
		}
		return _x;
	}
#endif 
};

template<class _Scalar>
class SOR_Solver {
	using _Tdyn = Matrix<_Scalar, Eigen::Dynamic, Eigen::Dynamic>;
public:
	explicit SOR_Solver(_Scalar _o = 1.25) :_omega(_o) {}
	
	_Tdyn operator()(const _Tdyn& A, const _Tdyn& b, int _MaxIteration = 100, double _MaxResidual = 1e-8)
	{
		const int ROW = A.rows(), COL = A.cols();
		_Tdyn LD = _Tdyn::Zero(ROW, COL), U = LD;
		for (int i = 0; i < ROW; ++i)
		{
			for (int j = 0; j < COL; ++j)
			{
				if (j > i)
				{
					U(i, j) = A(i, j);
				}
				else
				{
					LD(i, j) = A(i, j);
				}
			}
		}
		PenroseInverseSolver<_Scalar> pinv;
		auto invLD = pinv(LD);
		_Tdyn x = Eigen::MatrixX<_Scalar>::Zero(A.cols(), 1);
		for (int iter = 0; iter < _MaxIteration; ++iter)
		{
			x = (1 - _omega) * x + _omega * invLD * (b - U * x);
			_Scalar residual = (A * x - b).lpNorm<2>();
			if (residual < _MaxResidual)
			{                
				break;
			}
		}
		return x;
	}

	_Scalar _omega;
};

template<class _Scalar>
class ConjugateGradientMethodSolver {
	using _matType = Matrix<_Scalar, Eigen::Dynamic, Eigen::Dynamic>;
	using _vecType = Matrix<_Scalar, Eigen::Dynamic, 1>;
public:
	_vecType operator()(const _matType& A, const _vecType& b, int _MaxIteration = 100, double _MaxResidual = 1e-8)
	{
		const int COL = A.cols();
	
		_vecType _d, _r, _x = _vecType::Zero(COL, 1);
		_d = _r = b - A * _x;

		for (int iter = 0; iter < _MaxIteration; ++iter)
		{
			if (_r.lpNorm<2>() < _MaxResidual)
			{

				break;
			}
			auto tmp    = A * _d;
			auto _rIner = _r.squaredNorm();
			auto alpha  = _rIner / _d.dot(tmp);
			_x += alpha * _d;
			_r -= alpha * tmp;

			auto beta = _r.squaredNorm() / _rIner;
			_d = _r + beta * _d;
		}
		return _x;
	}
};

#endif // !_LINEARSYSTEM_SOLVER_HPP
