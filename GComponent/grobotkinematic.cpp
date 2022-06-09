#include "GComponent/grobotkinematic.h"

namespace GComponent{

namespace RobotKinematic{



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
