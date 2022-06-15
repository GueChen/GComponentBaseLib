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

#include <Eigen/dense>
#include <Eigen/Geometry>
#include <Eigen/SVD>

#include <algorithm>
#include <execution>
#include <vector>

namespace GComponent {
///	  Type Alias ���ͱ���
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

enum class SolveResult {
	Sucess				= 0,
	SizeNotMatch		= 1,
	OutOfLimitation		= 2,
	PrecisionNotMatch	= 3
};

/// <summary>
/// Get the ajoint transform matrices under exponential coordinates with specific thetas
/// <para>
/// ʹ��ָ���Ƕ�ֵ��ȡ�ض��������µ����ڱ任����
/// </para>
/// </summary>
/// <param name="out_transforms">	ref		{DynMats}	[out]	result							���������ڽӱ任����	</param>
/// <param name="expcoords">		cref	{Twists}	[in]	exponential coordinates			ָ������				</param>
/// <param name="thetas">			cref	{DynVec}	[in]	the value of thetas				ָ���Ƕ�ֵ����			</param>
/// <returns>								{bool}		[out]	flag of solving result			�������־			</returns>
template<class _Scaler>
bool	  Transforms(vector<SE3<_Scaler>>&					out_transforms,
						   const vector<Twist<_Scaler>>&	expcoords,
						   const DynVec<_Scaler>&			thetas);

/// <summary>
/// Solving the end effect Transform using Exponential Product under exponential coordinates with specific thetas
/// <para>
/// ʹ��ָ������ʽ���ָ���Ƕ�ֵ���ض�ĩ�˵����α任����
/// </para>
/// </summary>
/// <param name="out_mat">			ref		{SE3}		[out]	result							��������ĩ��λ�˾���	</param>
/// <param name="zero_mat">			cref	{SE3}		[in]	zero-point end transform matrix ���ĩ��λ�˾���		</param>
/// <param name="expcoords">		cref	{Twists}	[in]	exponential coordinates			ָ������				</param>
/// <param name="thetas">			cref	{DynVec}	[in]	the value of thetas				ָ���Ƕ�ֵ����			</param>
/// <returns>								{bool}		[out]	flag of solving result			�������־			</returns>
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
/// ʹ��ָ���Ƕȼ�������ſɱ�
/// </para>
/// </summary>
/// <param name="out_jacobian">		ref		{DynMat}	[out]	result							���������ſɱȾ���	</param>
/// <param name="expcoords">		cref	{Twists}	[in]	exponential coordinates			ָ������				</param>
/// <param name="thetas">			cref	{DynVec}	[in]	the value of thetas				ָ���Ƕ�ֵ����			</param>
/// <returns>								{bool}		[out]	flag of solving result			�������־			</returns>
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
/// ����ָ���Ƕ��µ���ռ�ͶӰ����
/// </para>
/// </summary>
/// <param name="out_matrix">		ref		{DynMat}	[out]	result							����������ռ�ͶӰ����</param>
/// <param name="expcoords">		cref	{Twists}	[in]	exponential coordinates			ָ������				</param>
/// <param name="thetas">			cref	{DynVec}	[in]	the value of thetas				ָ���Ƕ�ֵ����			</param>
/// <returns>								{bool}		[out]	flag of solving result			�������־			</returns>
template<class _Scaler>
bool	  NullSpaceProjection(DynMat<_Scaler>&				out_projection_matrix,
						const vector<Twist<_Scaler>>&		expcoords,
						const DynVec<_Scaler>&				thetas);

/// <summary>
/// Solving Inverse Kinematic Problem using numerical method called Least-Squared Method.
/// <para>
/// ʹ����С���˷���������˶�ѧ��ֵ��
/// </para>
/// </summary>
/// <param name="out_thetas">	ref		{DynVec}		[out]	result							������		</param>
/// <param name="zero_mat">		cref	{SE3}			[in]	zero-point end transform matrix ���ĩ��λ�˾���</param>
/// <param name="expcoords">	cref    {Twists}		[in]	exponential coordinates			ָ������		</param>
/// <param name="goal_mat">		cref	{SE3}			[in]	goal transform matrix			Ŀ��ĩ��λ�˾���</param>
/// <param name="init_guess">	cref	{DynVec}		[in]	initial guess of result theta   ��ֵ��²�ֵ	</param>
/// <param name="solver">		cref	{IKSolver}		[in]	linear system solver			����ϵͳ�����	</param>
/// <param name="Precision">	const	{Scaler}		[in]	precision						��⾫��		</param>
/// <param name="MaxIteration"> const	{int}			[in]	max iteration limit				����������	</param>
/// <param name="Scaler">		const	{Scaler}		[in]	decay scaler					����˥��ֵ		</param>
/// <returns>							{bool}			[out]	flag of solving result			�������־	</returns>
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

template<class _Scaler>
bool	  JacobianWithSE3(DynMat<_Scaler>&				out_jacobian,
						   SE3<_Scaler> &					out_matrix,
						   const vector<Twist<_Scaler>>&	expcoords,
						   const vector<SE3<_Scaler>>&    adj_matrices);

SE3d StandardDH(double alpha, double a, double theta0, double d, double theta = 0.0);
SE3f StandardDH(float alpha, float a, float theta0, float d, float theta = 0.0f);

SE3d ModifiedDH(double alpha, double a, double theta0, double d, double theta = 0.0);
SE3f ModifiedDH(float alpha, float a, float theta0, float d, float theta = 0.0f);

}
}

#include <GComponent/grobotkinematic-inl.hpp>

#endif
