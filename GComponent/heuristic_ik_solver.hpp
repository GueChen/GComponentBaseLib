#ifndef _HEURISTIC_IK_SOLVER_HPP
#define _HEURISTIC_IK_SOLVER_HPP

#include <GComponent/GTransform.hpp>
#include <GComponent/GNumerical.hpp>
#include <GComponent/GGeometry.hpp>

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

namespace RobotKinematic {
/// <summary>
/// inverse kinematic heuristic solver abstract base class
/// <para>逆运动学启发式求解器抽象基类，需实现 operator()(Args...)</para>
/// <para>派生子类</para>
/// <para>CCD - Cyclic coordinate descent</para>
/// <para> 1). ForwardCyclicCoordinateDescentSolver</para>
/// <para> 2). BackwardCyclicCoordinateDescentSolver</para>
/// </summary>
/// <typeparam name="_Scaler"></typeparam>
template<class _Scaler>
class HeuristicInverseKinematicSolver {
public:
	HeuristicInverseKinematicSolver()			= default;
	virtual ~HeuristicInverseKinematicSolver()  = default;
	virtual void operator()(SE3<_Scaler>&					mat_cur, 
							DynVec<_Scaler>&				thetas, 
							vector<SE3<_Scaler>>&			adj_matrices,
							const Vector<_Scaler, 3>&		goal_pos,
							const SE3<_Scaler>&				zero_mat,
							const vector<Twist<_Scaler>>&	exp_coords) const = 0;

protected:
	void IterateOnce(
					 SE3<_Scaler>&				mat_cur, 						
					 _Scaler&					relate_theta,
					 SE3<_Scaler>&				relate_adj_matrix,
					 const SE3<_Scaler>&		transform_matrix,
					 const Vector<_Scaler, 3>&	goal_pos,
					 const Twist<_Scaler>&		expcoord) const
	{
		// Get The Current Axis Message
		Twist<_Scaler>	   cur_axis		   = Adjoint(transform_matrix) * expcoord;
		Vector<_Scaler, 3> end_effect_pose = mat_cur.block(0, 3, 3, 1);
		Vector<_Scaler, 3> axis_dir		   = cur_axis.block(0, 0, 3, 1), 
						   axis_v		   = cur_axis.block(3, 0, 3, 1);
		Vector<_Scaler, 3> axis_pose	   = axis_dir.cross(axis_v);
		
		// 
		Vector<_Scaler, 3> u = end_effect_pose - axis_pose,
						   v = goal_pos - axis_pose;
		u -= axis_dir * axis_dir.dot(u);
		v -= axis_dir * axis_dir.dot(v);
		_Scaler	delta = /*Clamp(axis_dir.dot(GetRotateAxisAngleFrom2Vec(u, v)), _Scaler(-MyPI / 4.0), _Scaler(MyPI / 4.0));*/
			axis_dir.dot(GetRotateAxisAngleFrom2Vec(u, v));

		mat_cur				 = ExpMapping(cur_axis, delta) * mat_cur;
		relate_theta		+= delta;
		relate_adj_matrix	 = ExpMapping(expcoord, relate_theta);
	}
};

/// <summary>
/// Cyclic coordinate descent solver using forward iteration form, form base to top.
/// <para>前向 CCD 求解器，自基座至末端</para>
/// </summary>
/// <typeparam name="_Scaler"></typeparam>
template<class _Scaler>
class ForwardCyclicCoordinateDescentSolver : public HeuristicInverseKinematicSolver<_Scaler> {
public:
	ForwardCyclicCoordinateDescentSolver()  = default;
	~ForwardCyclicCoordinateDescentSolver() = default;

	void operator()(SE3<_Scaler>&					mat_cur, 
					DynVec<_Scaler>&				thetas, 
					vector<SE3<_Scaler>>&			adj_matrices,
					const Vector<_Scaler, 3>&		goal_pos,
					const SE3<_Scaler>&				zero_mat,
					const vector<Twist<_Scaler>>&	exp_coords) const override
	{
		SE3<_Scaler> transform_matrix = SE3<_Scaler>::Identity();
		for (int i = 0; i < thetas.size(); ++i) {	
			HeuristicInverseKinematicSolver<_Scaler>::IterateOnce(mat_cur, thetas[i], adj_matrices[i], transform_matrix, goal_pos, exp_coords[i]);
			transform_matrix *= adj_matrices[i];
		}
	}
};

/// <summary>
/// Cyclic coordinate descent solver using backward iteration form, from top to base.
/// <para>后向 CCD 求解器，自末端至基座</para>
/// </summary>
/// <typeparam name="_Scaler"></typeparam>
template<class _Scaler>
class BackwardCyclicCoordinateDescentSolver : public HeuristicInverseKinematicSolver<_Scaler> {
public:
	BackwardCyclicCoordinateDescentSolver()  = default;
	~BackwardCyclicCoordinateDescentSolver() = default;

	void operator()(SE3<_Scaler>&					mat_cur, 
					DynVec<_Scaler>&				thetas, 
					vector<SE3<_Scaler>>&			adj_matrices,
					const Vector<_Scaler, 3>&		goal_pos,
					const SE3<_Scaler>&				zero_mat,
					const vector<Twist<_Scaler>>&	exp_coords) const override
	{
		SE3<_Scaler> transform_matrix = mat_cur * InverseSE3(zero_mat);
		for (int i = thetas.size() - 1; i >= 0; --i) {
			transform_matrix *= InverseSE3(adj_matrices[i]);
			HeuristicInverseKinematicSolver<_Scaler>::IterateOnce(mat_cur, thetas[i], adj_matrices[i], transform_matrix, goal_pos, exp_coords[i]);
		}
	}
};
}
}
#endif // !_HEURISTIC_IK_SOLVER_HPP

