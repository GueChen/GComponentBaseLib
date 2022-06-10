#ifndef _GROBOTKINEMATIC_INL_HPP
#define _GROBOTKINEMATIC_INL_HPP

#include <GComponent/grobotkinematic.h>

namespace GComponent {
namespace RobotKinematic {

template<class _Scaler>
bool Transforms(vector<SE3<_Scaler>>& out_transforms, const vector<Twist<_Scaler>>& expcoords, const DynVec<_Scaler>& thetas)
{
	int n = expcoords.size();
	if (thetas.rows() != n) return false;

	out_transforms.resize(n);
	std::transform(std::execution::par_unseq, expcoords.begin(), expcoords.end(), thetas.begin(),
		out_transforms.begin(),
		[](const Twist<_Scaler>& epi, const _Scaler& theta) {
			return ExpMapping(epi, theta);
		});

	return true;
}

template<class _Scaler>
bool ForwardKinematic(SE3<_Scaler>& out_mat, const SE3<_Scaler>& ini_mat, const vector<Twist<_Scaler>>& expcoords, const DynVec<_Scaler>& thetas)
{
	if (thetas.rows() != expcoords.size()) return false;

	vector<SE3<_Scaler>> adj_matrices;
	Transforms(adj_matrices, expcoords, thetas);
	ForwardKinematic(out_mat, ini_mat, adj_matrices);
	return true;
}
template<class _Scaler>
bool ForwardKinematic(SE3<_Scaler>& out_matrix, const SE3<_Scaler>& ini_mat, const vector<SE3<_Scaler>>& adj_matrices)
{
	out_matrix.setIdentity();
	for (auto& mat : adj_matrices) {
		out_matrix *= mat;
	}
	out_matrix *= ini_mat;
	return true;
}

template<class _Scaler>
bool Jacobian(DynMat<_Scaler>& out_jacobian, const vector<Twist<_Scaler>>& expcoords, const DynVec<_Scaler>& thetas)
{
	if (expcoords.size() != thetas.rows()) return false;
	vector<SE3<_Scaler>> adj_matrices;
	Transforms(adj_matrices, expcoords, thetas);
	Jacobian(out_jacobian, expcoords, adj_matrices);
	return true;
}
template<class _Scaler>
bool Jacobian(DynMat<_Scaler>& out_jacobian, const vector<Twist<_Scaler>>& expcoords, const vector<SE3<_Scaler>>& adj_matrices)
{
	if (expcoords.size() != adj_matrices.size()) return false;
	SE3<_Scaler> local_mat = SE3<_Scaler>::Identity();
	JacobianWithSE3(out_jacobian, local_mat, expcoords, adj_matrices);
	return true;
}

template<class _Scaler>
bool JacobianWithSE3(DynMat<_Scaler>& out_jacobian, SE3<_Scaler>& out_matrix, const vector<Twist<_Scaler>>& expcoords, const vector<SE3<_Scaler>>& adj_matrices)
{
	if (expcoords.size() != adj_matrices.size()) return false;
	const int N = expcoords.size();
	out_matrix.setIdentity();
	out_jacobian.resize(6, N);
	for (int i = 0; i < N; ++i) {
		out_jacobian.block(0, i, 6, 1) = Adjoint(out_matrix) * expcoords[i];
		out_matrix *= adj_matrices[i];
	}
	return true;
}

template<class _Scaler>
bool NullSpaceProjection(DynMat<_Scaler>& out_projection_matrix, const vector<Twist<_Scaler>>& expcoords, const DynVec<_Scaler>& thetas)
{
	if (expcoords.size() != thetas.rows()) return false;
	const int N = expcoords.size();
	DynMat<_Scaler> J;
	Jacobian(J, expcoords, thetas);
	Eigen::JacobiSVD<DynMat<_Scaler>> svd(J, Eigen::ComputeFullU | Eigen::ComputeFullV);
	DynMat<_Scaler> pinv_J  = svd.solve(DynMat<_Scaler>::Identity(6, 6));
	out_projection_matrix = DynMat<_Scaler>::Identity(N, N) - pinv_J * J;
	return true;
}

template<class _Scaler>
bool InverseKinematic(DynVec<_Scaler>& out_thetas, const SE3<_Scaler>& zero_mat, const vector<Twist<_Scaler>>& expcoords, const SE3<_Scaler>& goal_mat, const DynVec<_Scaler>& init_guess, const IKSolver<_Scaler>& solver, const double Precision, const int MaxIteration, const double Scaler)
{
	if (init_guess.size() != expcoords.size()) return false;
	const double		InvScaler = 1.0f / Scaler;

	SE3<_Scaler>		mat_cur,						// current SE3
						mat_new;						// trial error SE3
	vector<SE3<_Scaler>>adj_matrices;					// adjoint matrices
	DynVec<_Scaler>		thetas		= init_guess,		// update val
						thetas_new;
	DynMat<_Scaler>		A,								// Jacobian
						x,								// val update direction
						b,								// residual vec
						b_new;							
	_Scaler				residual		= 0.0f,
						residual_new	= 0.0f,
						decay			= 1.0f,			// step length													
						tolerance		= 0.0f;																
	int					iter_count		= 0;																
	bool				precision_flag  = false;
	/* Measure residual between Desire and Cur 
	*	b   = log(desire - cur) 
	*	res = norm(b) 
	* */	   	 
	Transforms(adj_matrices, expcoords, init_guess);
	ForwardKinematic(mat_cur, zero_mat, adj_matrices);
	thetas	 = init_guess;
	b		 = LogMapSE3Tose3((goal_mat * InverseSE3(mat_cur)).eval());
	residual = b.norm();
	
	/* Use this to trial and error until the residual under tolerance */
	auto trial_error = [&]() {
		thetas_new = thetas + decay * x;
		
		Transforms(adj_matrices, expcoords, thetas_new);
		ForwardKinematic(mat_new, zero_mat, adj_matrices);
		b_new		 = LogMapSE3Tose3((goal_mat * InverseSE3(mat_new)).eval());
		residual_new = b_new.norm();
	};

	/* Accept Step and Update result */
	auto accept = [&]() {
		thetas	  = thetas_new;		
		mat_cur   = mat_new;
		residual  = residual_new;
		b		  = b_new;
	};

	/* Main Loop Solving the SQP Problem */
	while (iter_count < MaxIteration && residual > Precision)
	{
		Jacobian(A, expcoords, adj_matrices);
		
		x = solver(A, b);
		trial_error();

		while(residual < tolerance + residual_new) {			
			decay *= Scaler;
			if (decay < Precision) {
				break;
			}
			trial_error();			
		}

		if (decay < Precision || abs(residual_new - residual) < 1e-5) break;
		if (decay < 1.0f) {
			decay = 1.0f;
		}

		accept();
		
		tolerance = -residual / 20.0;
	
		++iter_count;
	}

	/* Accept the Result or Refuse */
	if (iter_count > MaxIteration && residual > Precision) {
		return false;
	}
	else {
		out_thetas = thetas;
		return true;
	}
}

}
}

#endif
