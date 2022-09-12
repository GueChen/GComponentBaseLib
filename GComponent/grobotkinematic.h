/**
 *  @file  	gkinematic.h
 *  @brief 	This file contains Robot Kinematic Caculation Methods.
 *  @author Gue Chen<guechen@buaa.edu.cn>
 *  @date 	Nov 27th, 2021
 *  @update Jun  2nd, 2022	reform the file name add some overwrites
 **/
#ifndef _GROBOTKINE_H
#define _GROBOTKINE_H

#include <LSSolver/LinearSystemSolver.hpp>
#include <GComponent/GTransform.hpp>
#include <GComponent/GNumerical.hpp>
#include <GComponent/GGeometry.hpp>
#include <GComponent/heuristic_ik_solver.hpp>

#include <Eigen/dense>
#include <Eigen/Geometry>
#include <Eigen/SVD>

#include <algorithm>
#include <execution>
#include <vector>

namespace GComponent {
///	  Type Alias 类型别名
using std::vector;
using SE3d  = Eigen::Matrix4d;
using SE3f  = Eigen::Matrix4f;
using Vec3d = Eigen::Vector3d;
using Vec3f = Eigen::Vector3f;
using Eigen::Affine3f;
using Eigen::Affine3d;
using Eigen::AngleAxisd;
using Eigen::AngleAxisf;
template<class _Scaler>
using DynVec = Eigen::Vector<_Scaler, Eigen::Dynamic>;
template<class _Scaler>
using DynMat = Eigen::Matrix<_Scaler, Eigen::Dynamic, Eigen::Dynamic>;

namespace RobotKinematic{
template<class _Scaler>
using IKSolver = LinearSystemSolver<_Scaler>;

template<class _Scaler>
class NoIKSolver : public IKSolver<_Scaler> {
	using _DynMat = Matrix<_Scaler, Eigen::Dynamic, Eigen::Dynamic>;
public:
	_DynMat operator()(const _DynMat& A, const _DynMat& b) override
	{
		return _DynMat::Zero(A.cols(), 1);
	}
};

enum class SolveResult {
	Sucess				= 0,
	SizeNotMatch		= 1,
	OutOfLimitation		= 2,
	PrecisionNotMatch	= 3
};

/// <summary>
/// Get the ajoint transform matrices under exponential coordinates with specific thetas
/// <para>
/// 使用指定角度值获取特定旋量轴下的相邻变换矩阵
/// </para>
/// </summary>
/// <param name="out_transforms">	ref		{DynMats}	[out]	result							计算结果，邻接变换矩阵	</param>
/// <param name="expcoords">		cref	{Twists}	[in]	exponential coordinates			指数坐标				</param>
/// <param name="thetas">			cref	{DynVec}	[in]	the value of thetas				指定角度值向量			</param>
/// <returns>								{bool}		[out]	flag of solving result			求解结果标志			</returns>
template<class _Scaler>
bool	  Transforms  (vector<SE3<_Scaler>>&			 out_transforms,
					   const vector<Twist<_Scaler>>&	 expcoords,
					   const DynVec<_Scaler>&			 thetas);

// TODO: add comments
template<class _Scalar>
bool	  Differential(vector<SE3<_Scalar>>&			 out_transforms,
					   const vector<Twist<_Scalar>>&	 expcoords,
					   const DynVec<_Scalar>&			 thetas);

/// <summary>
/// Solving the end effect Transform using Exponential Product under exponential coordinates with specific thetas
/// <para>
/// 使用指数积公式求解指定角度值下特定末端点的齐次变换矩阵
/// </para>
/// </summary>
/// <param name="out_mat">			ref		{SE3}		[out]	result							计算结果，末端位姿矩阵	</param>
/// <param name="zero_mat">			cref	{SE3}		[in]	zero-point end transform matrix 零点末端位姿矩阵		</param>
/// <param name="expcoords">		cref	{Twists}	[in]	exponential coordinates			指数坐标				</param>
/// <param name="thetas">			cref	{DynVec}	[in]	the value of thetas				指定角度值向量			</param>
/// <returns>								{bool}		[out]	flag of solving result			求解结果标志			</returns>
template<class _Scaler>
bool	  ForwardKinematic(SE3<_Scaler>&					out_mat, 
						   const SE3<_Scaler>&				zero_mat,
						   const vector<Twist<_Scaler>>&	expcoords, 
						   const DynVec<_Scaler>&			thetas);
template<class _Scaler>
bool	  ForwardKinematic(SE3<_Scaler>&					out_mat,
						   const SE3<_Scaler>&				ini_mat,
						   const vector<SE3<_Scaler>>&		adj_matrices);

/// <summary>
/// Get Analytic Jacobian with Specific thetas under exponential coordinates
/// <para>
/// 使用指定角度计算解析雅可比
/// </para>
/// </summary>
/// <param name="out_jacobian">		ref		{DynMat}	[out]	result							计算结果，雅可比矩阵	</param>
/// <param name="expcoords">		cref	{Twists}	[in]	exponential coordinates			指数坐标				</param>
/// <param name="thetas">			cref	{DynVec}	[in]	the value of thetas				指定角度值向量			</param>
/// <returns>								{bool}		[out]	flag of solving result			求解结果标志			</returns>
template<class _Scaler>
bool	  Jacobian(DynMat<_Scaler>&							out_jacobian,
							const vector<Twist<_Scaler>>&	expcoords,
							const DynVec<_Scaler>&			thetas);
template<class _Scaler>
bool	  Jacobian(DynMat<_Scaler>&							out_jacobian,
							const vector<Twist<_Scaler>>&	expcoords,
							const vector<SE3<_Scaler>>&		adj_matrices);

/// <summary>
/// Get Null Space Projection Matrix with Specific thetas under exponential coordinates
/// <para>
/// 计算指定角度下的零空间投影矩阵
/// </para>
/// </summary>
/// <param name="out_matrix">		ref		{DynMat}	[out]	result							计算结果，零空间投影矩阵</param>
/// <param name="expcoords">		cref	{Twists}	[in]	exponential coordinates			指数坐标				</param>
/// <param name="thetas">			cref	{DynVec}	[in]	the value of thetas				指定角度值向量			</param>
/// <returns>								{bool}		[out]	flag of solving result			求解结果标志			</returns>
template<class _Scaler>
bool	  NullSpaceProjection(DynMat<_Scaler>&				out_projection_matrix,
						const vector<Twist<_Scaler>>&		expcoords,
						const DynVec<_Scaler>&				thetas);

/// <summary>
/// Solving Inverse Kinematic Problem using numerical method called Least-Squared Method.
/// <para>
/// 使用最小二乘方法求解逆运动学数值解
/// </para>
/// </summary>
/// <param name="out_thetas">	ref		{DynVec}		[out]	result							计算结果		</param>
/// <param name="zero_mat">		cref	{SE3}			[in]	zero-point end transform matrix 零点末端位姿矩阵</param>
/// <param name="expcoords">	cref    {Twists}		[in]	exponential coordinates			指数坐标		</param>
/// <param name="goal_mat">		cref	{SE3}			[in]	goal transform matrix			目标末端位姿矩阵</param>
/// <param name="init_guess">	cref	{DynVec}		[in]	initial guess of result theta   初值点猜测值	</param>
/// <param name="solver">		cref	{IKSolver}		[in]	linear system solver			线性系统求解器	</param>
/// <param name="Precision">	const	{Scaler}		[in]	precision						求解精度		</param>
/// <param name="MaxIteration"> const	{int}			[in]	max iteration limit				最大迭代次数	</param>
/// <param name="Scaler">		const	{Scaler}		[in]	decay scaler					步长衰减值		</param>
/// <returns>							{bool}			[out]	flag of solving result			求解结果标志	</returns>
template<class _Scaler>
bool	  InverseKinematic(DynVec<_Scaler>&				out_thetas,
						   const SE3<_Scaler>&			zero_mat,
						   const vector<Twist<_Scaler>>&expcoords,
						   const SE3<_Scaler>&			goal_mat, 
						   const DynVec<_Scaler>&		init_guess,
						   const IKSolver<_Scaler>&		solver, 
						   const double					Precision	 = 1e-5f, 
						   const int					MaxIteration = 50,
						   const double					Scaler		 = 0.3f);

/// <summary>
/// 
/// </summary>
/// <param name="out_thetas">	ref		{DynVec}		[out]	result							计算结果		</param>
/// <param name="zero_mat">		cref	{SE3}			[in]	zero-point end transform matrix 零点末端位姿矩阵</param>
/// <param name="expcoords">	cref    {Twists}		[in]	exponential coordinates			指数坐标		</param>
/// <param name="goal_mat">		cref	{SE3}			[in]	goal transform matrix			目标末端位姿矩阵</param>
/// <param name="init_guess">	cref	{DynVec}		[in]	initial guess of result theta   初值点猜测值	</param>
/// <param name="solver">		cref	{HerusticSolver}[in]	heuristic solver				启发式求解器	</param>
/// <param name="Precision">	const	{Scaler}		[in]	precision						求解精度		</param>
/// <param name="MaxIteration"> const	{int}			[in]	max iteration limit				最大迭代次数	</param>
/// <param name="Scaler">		const	{Scaler}		[in]	decay scaler					步长衰减值		</param>
/// <returns>							{bool}			[out]	flag of solving result			求解结果标志	</returns>
template<class _Scaler>
bool	  InverseKinematicHeuristic(
						   DynVec<_Scaler>&				out_thetas,
						   const SE3<_Scaler>&			zero_mat,
						   const vector<Twist<_Scaler>>&expcoords,
						   const SE3<_Scaler>&			goal_mat,
						   const DynVec<_Scaler>&		init_guess,
						   const HeuristicInverseKinematicSolver<_Scaler>&
														solver,
						   const double					Precision	= 1e-5f,
						   const int					MaxIteration= 50);


template<class _Scaler>
bool	  JacobianWithSE3(DynMat<_Scaler>&				out_jacobian,
						   SE3<_Scaler> &				out_matrix,
						   const vector<Twist<_Scaler>>&expcoords,
						   const vector<SE3<_Scaler>>&  adj_matrices);

/// <summary>
/// Get transform matrix using standard Denavit-Hartenberg method
/// <para>
/// 使用标准 DH 方法获取位姿变换矩阵
/// </para>
/// </summary>
/// <param name="alpha">	{double}	 [in]  joint twist angle(rad.)			关节扭角(rad.)		</param>
/// <param name="a">		{double}	 [in]  link length(m.)					连杆长度(m.)		</param>
/// <param name="theta0">	{double}	 [in]  initial joint angle bais(rad.)	关节转角偏转(rad.)	</param>
/// <param name="d">		{double}	 [in]  link bias(m.)					连杆偏置(m.)		</param>
/// <param name="theta">	{double}	 [in]  joint angle(rad.)				关节转角(rad.)		</param>
/// <returns>				{SE3 double} [out] transform matrix					位姿变换矩阵		</returns>
SE3d StandardDH(double alpha, double a, double theta0, double d, double theta = 0.0);
SE3f StandardDH(float alpha, float a, float theta0, float d, float theta = 0.0f);

/// <summary>
/// Get transform matrix using modified Denavit-Hartenberg method
/// <para>
/// 使用改进 DH 方法获取位姿变换矩阵
/// </para>
/// <param name="alpha">	{double}	 [in]  joint twist angle(rad.)			关节扭角(rad.)		</param>
/// <param name="a">		{double}	 [in]  link length(m.)					连杆长度(m.)		</param>
/// <param name="theta0">	{double}	 [in]  initial joint angle bais(rad.)	关节转角偏转(rad.)	</param>
/// <param name="d">		{double}	 [in]  link bias(m.)					连杆偏置(m.)		</param>
/// <param name="theta">	{double}	 [in]  joint angle(rad.)				关节转角(rad.)		</param>
/// <returns>				{SE3 double} [out] transform matrix					位姿变换矩阵		</returns>
SE3d ModifiedDH(double alpha, double a, double theta0, double d, double theta = 0.0);
SE3f ModifiedDH(float alpha, float a, float theta0, float d, float theta = 0.0f);

}
}

#include <GComponent/grobotkinematic-inl.hpp>

#endif
