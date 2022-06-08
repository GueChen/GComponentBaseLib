#include "GComponent/grobotkinematic.h"

namespace GComponent{

namespace RobotKinematic{
bool Transforms(vector<SE3<float>>& out_transforms, const vector<Twist<float>>& expcoords, const DynVec<float>& thetas)
{
	int n = expcoords.size();
	if (thetas.rows() != n) return false;

	out_transforms.resize(n);
	std::transform(std::execution::par_unseq, expcoords.begin(), expcoords.end(), thetas.begin(),
		out_transforms.begin(),
		[](const Twist<float>& epi, const float& theta) {
			return ExpMapping(epi, theta);
		});

	return true;
}

bool ForwardKinematic(SE3<float>& out_mat, const SE3<float>& ini_mat, const vector<Twist<float>>& expcoords, const DynVec<float>& thetas)
{
	if (thetas.rows() != expcoords.size()) return false;

	vector<SE3<float>> adj_matrices; 
	Transforms(adj_matrices, expcoords, thetas);
	ForwardKinematic(out_mat, ini_mat, adj_matrices);
	return true;
}
bool ForwardKinematic(SE3<float>& out_matrix, const SE3<float>& ini_mat, const vector<SE3<float>>& adj_matrices)
{
	out_matrix.setIdentity();
	for (auto& mat : adj_matrices) {
		out_matrix *= mat;
	}
	out_matrix *= ini_mat;
	return true;
}

bool Jacobian(DynMat<float>& out_jacobian, const vector<Twist<float>>& expcoords, const DynVec<float>& thetas)
{
	if (expcoords.size() != thetas.rows()) return false;
	vector<SE3<float>> adj_matrices;
	Transforms(adj_matrices, expcoords, thetas);
	Jacobian(out_jacobian, expcoords, adj_matrices);
	return true;
}
bool Jacobian(DynMat<float>& out_jacobian, const vector<Twist<float>>& expcoords, const vector<SE3<float>>& adj_matrices)
{
	if (expcoords.size() != adj_matrices.size()) return false;
	SE3<float> local_mat = SE3<float>::Identity();
	JacobianWithSE3(out_jacobian, local_mat, expcoords, adj_matrices);
	return true;
}

bool JacobianWithSE3(DynMat<float>& out_jacobian, SE3<float>& out_matrix, const vector<Twist<float>>& expcoords, const vector<SE3<float>>& adj_matrices)
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

bool NullSpaceProjection(DynMat<float>& out_projection_matrix, const vector<Twist<float>>& expcoords, const DynVec<float>& thetas)
{
	if (expcoords.size() != thetas.rows()) return false;
	const int N = expcoords.size();
	DynMat<float> J;
	Jacobian(J, expcoords, thetas);
	Eigen::JacobiSVD<DynMat<float>> svd(J, Eigen::ComputeFullU | Eigen::ComputeFullV);
	DynMat<float> pinv_J  = svd.solve(DynMat<float>::Identity(6, 6));
	out_projection_matrix = DynMat<float>::Identity(N, N) - pinv_J * J;
	return true;
}

bool InverseKinematic(DynVec<float>& out_thetas, const SE3<float>& zero_mat, const vector<Twist<float>>& expcoords, const SE3<float>& goal_mat, const DynVec<float>& init_guess, const IKSolver<float>& solver, const float Precision, const int MaxIteration, const float Scaler)
{
	if (init_guess.size() != expcoords.size()) return false;
	const float			InvScaler = 1.0f / Scaler;

	SE3<float>			mat_cur,						// current SE3
						mat_new;						// trial error SE3
	vector<SE3<float>>  adj_matrices;					// adjoint matrices
	DynVec<float>		thetas		= init_guess,		// update val
						thetas_new;
	DynMat<float>		A,								// Jacobian
						x,								// val update direction
						b,								// residual vec
						b_new;							
	float				residual		= 0.0f,			
						residual_new	= 0.0f,
						decay			= 1.0f,			// step length													
						tolerance		= -1e2f;																
	int					iter_count		= 0;																
	
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
			trial_error();
			if (decay * x.norm() < Precision) break;
		}

		if (residual < tolerance + residual_new) break;
		if (decay < 1.0f) {
			decay *= InvScaler;
		}

		accept();
		
		tolerance = residual / 10.0f;
	
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

SE3d StandardDH(double alpha, double a, double theta0, double d, double theta)
{
	Affine3d transform; 
	transform.setIdentity();

	transform.rotate(AngleAxisd(theta + theta0, Vec3d::UnitZ()));
	transform.translate(d * Vec3d::UnitZ());
	transform.rotate(AngleAxisd(alpha, Vec3d::UnitX()));
	transform.translate(a * Vec3d::UnitX());

	return transform.matrix();
}

SE3f StandardDH(float alpha, float a, float theta0, float d, float theta)
{
	Affine3f transform;
	transform.setIdentity();

	transform.rotate(AngleAxisf(theta + theta0, Vec3f::UnitZ()));
	transform.translate(d * Vec3f::UnitZ());
	transform.rotate(AngleAxisf(alpha, Vec3f::UnitX()));
	transform.translate(a * Vec3f::UnitX());

	return transform.matrix();
}

SE3d ModifiedDH(double alpha, double a, double theta0, double d, double theta)
{
	Affine3d transform; 
	transform.setIdentity();

	transform.rotate(AngleAxisd(alpha, Vec3d::UnitX()));
	transform.translate(a * Vec3d::UnitX());
	transform.rotate(AngleAxisd(theta + theta0, Vec3d::UnitZ()));
	transform.translate(d * Vec3d::UnitZ());

	return transform.matrix();
}

SE3f ModifiedDH(float alpha, float a, float theta0, float d, float theta)
{
	Affine3f transform;
	transform.setIdentity();

	transform.rotate(AngleAxisf(alpha, Vec3f::UnitX()));
	transform.translate(a * Vec3f::UnitX());
	transform.rotate(AngleAxisf(theta + theta0, Vec3f::UnitZ()));
	transform.translate(d * Vec3f::UnitZ());
	
	return transform.matrix();
}

}
}
