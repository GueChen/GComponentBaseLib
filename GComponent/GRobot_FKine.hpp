#pragma once
#include <Eigen/dense>
#include <Eigen/Geometry>

namespace GComponent {
using SE3d = Eigen::Matrix4d;
using Eigen::Vector3d;

inline SE3d StandardDH(double _alpha, double _a, double _theta, double _d, double _theta0 = 0.0f)
{
	Eigen::Transform<double, 3, Eigen::Affine> transform; transform.setIdentity();
	
	transform.rotate(Eigen::AngleAxisd(_theta + _theta0, Vector3d(0.0f, 0.0f, 1.0f)));
	transform.translate(Vector3d(0.0f, 0.0f, _d));
	transform.rotate(Eigen::AngleAxisd(_alpha, Vector3d(1.0f, 0.0f, 0.0f)));
	transform.translate(Vector3d(_a, 0.0f, 0.0f));

	return transform.matrix();
}

inline SE3d ModifiledDH(double _alpha, double _a, double _theta, double _d, double _theta0 = 0.0f)
{
	Eigen::Transform<double, 3, Eigen::Affine> transform; transform.setIdentity();

	transform.rotate(Eigen::AngleAxisd(_alpha, Vector3d(1.0f, 0.0f, 0.0f)));
	transform.translate(Vector3d(_a, 0.0f, 0.0f));
	transform.rotate(Eigen::AngleAxisd(_theta + _theta0, Vector3d(0.0f, 0.0f, 1.0f)));
	transform.translate(Vector3d(0.0f, 0.0f, _d));
	
	return transform.matrix();
}

template<class ... Args>
void ForwardKinematic(Args... parms)
{

}
}