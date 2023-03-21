#ifndef __TYPES_GCOM_H
#define __TYPES_GCOM_H
#include <Eigen/Core>


template<class _Scaler>
using SE3 = Eigen::Matrix4<_Scaler>;
template<class _Scaler>
using SO3 = Eigen::Matrix3<_Scaler>;
template<class _Scaler>
using Twist = Eigen::Vector<_Scaler, 6>;
template<class _Scaler>
using AdMatrix = Eigen::Matrix<_Scaler, 6, 6>;

using SE3d	     = SE3<double>;
using SE3f       = SE3<float>;
using Twistd     = Twist<double>;
using Twistf     = Twist<float>;
using SO3f       = SO3<float>;
using SO3d		 = SO3<double>;
using Vec3d      = Eigen::Vector3d;
using Vec3f      = Eigen::Vector3f;
using Vec3dT     = Eigen::Matrix<double, 1, 3>;
using Twistd     = Eigen::Matrix<double, 6, 1>;
using DynMatrixd = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>;

#endif // !__TYPES_GCOM_H
