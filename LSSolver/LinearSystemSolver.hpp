#pragma once
#include <Eigen/LU>
#include <Eigen/dense>
#include <Eigen/core>

#include "PenroseInverse.hpp"
using Eigen::Matrix;

template<unsigned Cow, unsigned Row>
class StaticInverseSolver {
	using _TA = Matrix<double, Cow, Row>;
	using _Tb = Matrix<double, Cow, 1>;
	using _Tx = Matrix<double, Row, 1>;
public:
	inline _Tx operator()(const _TA& A, const _Tb& b)
	{
		if(A.rows() == A.cols())
			return A.inverse() * b;
		if (A.rows() > A.cols())
			assert("The Matrix Rows can't more then Cols, Equation has No Solution.");
		PenroseInverseSolver pinv;
		return pinv(A) * b;
	}
};

class DynamicInverseSolver {
	using _Tdyn = Matrix<double, Eigen::Dynamic, Eigen::Dynamic>;
public:
	inline _Tdyn operator()(const _Tdyn& A, const _Tdyn& b) 
	{
		if(A.rows() == A.cols())
			return A.inverse() * b;
		if (A.rows() > A.cols())
			assert("The Matrix Rows can't more then Cols, Equation has No Solution.");
		PenroseInverseSolver pinv;
		return pinv(A) * b;
	}
};

class JacobiSolver {
	using _Tdyn = Matrix<double, Eigen::Dynamic, Eigen::Dynamic>;
public:
	inline _Tdyn operator()(const _Tdyn& A, const _Tdyn& b, int _MaxIteration = 100, double _MaxResidual = 1e-8) 
	{
		_Tdyn LU = A, D = Eigen::MatrixXd::Zero(A.cols(), A.rows());
		for (int i = 0; i < A.rows(); ++i)
		{
			std::swap(D(i, i),LU(i, i));
			D(i, i) = 1.0 / D(i, i);
		}
		_Tdyn _x = Eigen::MatrixXd::Zero(A.cols(), 1);
		for (int iter = 0; iter < _MaxIteration; ++iter)
		{
			_x = D * (b - LU * _x);
			double residual = (A * _x - b).lpNorm<2>();
			if (residual < _MaxResidual)
			{
				//std::cout << "Jacobi iterations:" << iter << std::endl;
				break;
			}
		}
		return _x;
	}
};

class GaussSeidelSolver {
	using _Tdyn = Matrix<double, Eigen::Dynamic, Eigen::Dynamic>;
public:
	inline _Tdyn operator()(const _Tdyn& A, const _Tdyn& b, int _MaxIteration = 100, double _MaxResidual = 1e-8)
	{
		const unsigned ROW = A.rows(), COL = A.cols();
		_Tdyn LD = Eigen::MatrixXd::Zero(ROW, COL), U = LD;
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
		PenroseInverseSolver pinv;
		auto invLD = pinv(LD);
		_Tdyn _x = Eigen::MatrixXd::Zero(A.cols(), 1);
		for (int iter = 0; iter < _MaxIteration; ++iter)
		{
			_x = invLD * (b - U * _x);
			double residual = (A * _x - b).lpNorm<2>();
			if (residual < _MaxResidual)
			{
                //std::cout << "Gauss-Seidel optimizal iterations:" << iter << std::endl;
				break;
			}
		}
		return _x;
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

class SOR_Solver {
	using _Tdyn = Matrix<double, Eigen::Dynamic, Eigen::Dynamic>;
public:
	explicit SOR_Solver(double _o = 1.25) :_omega(_o) {}
	
	inline _Tdyn operator()(const _Tdyn& A, const _Tdyn& b, int _MaxIteration = 100, double _MaxResidual = 1e-8)
	{
		const unsigned ROW = A.rows(), COL = A.cols();
		_Tdyn LD = Eigen::MatrixXd::Zero(ROW, COL), U = LD;
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
		PenroseInverseSolver pinv;
		auto invLD = pinv(LD);
		_Tdyn _x = Eigen::MatrixXd::Zero(A.cols(), 1);
		for (int iter = 0; iter < _MaxIteration; ++iter)
		{
			_x = (1 - _omega) * _x + _omega * invLD * (b - U * _x);
			double residual = (A * _x - b).lpNorm<2>();
			if (residual < _MaxResidual)
			{
                //std::cout << "Gauss-Seidel optimizal iterations:" << iter << std::endl;
				break;
			}
		}
		return _x;
	}

	double _omega;
};

class ConjugateGradientMethodSolver {
	using _matType = Matrix<double, Eigen::Dynamic, Eigen::Dynamic>;
	using _vecType = Matrix<double, Eigen::Dynamic, 1>;
public:
	_vecType operator()(const _matType& A, const _vecType& b, int _MaxIteration = 100, double _MaxResidual = 1e-8)
	{
		const unsigned ROW = A.rows(), COL = A.cols();
	
		_vecType _d, _r, _x = _vecType::Zero(COL, 1);
		_d = _r = b - A * _x;

		for (int iter = 0; iter < _MaxIteration; ++iter)
		{
			if (_r.lpNorm<2>() < _MaxResidual)
			{

				break;
			}
			auto tmp = A * _d;
			auto _rIner = _r.squaredNorm();
			auto alpha = _rIner / _d.dot(tmp);
			_x += alpha * _d;
			_r -= alpha * tmp;

			auto beta = _r.squaredNorm() / _rIner;
			_d = _r + beta * _d;
		}
		return _x;
	}
};
