#ifndef _GMATRIX_OPERATION_H
#define _GMATRIX_OPERATION_H

#include <Eigen/Dense>

#include <tuple>

namespace GComponent{

using std::tuple;
template<class _Scaler>
using Mat3 = Eigen::Matrix3<_Scaler>;
template<class _Scaler>
using Vec3 = Eigen::Vector3<_Scaler>;
template<class _Scaler>
using Mat4 = Eigen::Matrix4<_Scaler>;
template<class _Scaler>
using Vec4 = Eigen::Vector4<_Scaler>;
template<class _Scaler>
using SO3  = Mat3<_Scaler>;
template<class _Scaler>
using SE3  = Mat4<_Scaler>;



}
#endif // !_GMATRIX_OPERATION_H
