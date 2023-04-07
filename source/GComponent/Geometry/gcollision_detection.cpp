#include "GComponent/Geometry/gcollision_detection.h"
#include "GComponent/Geometry/gdistance.h"

#include <GComponent/gtransform.hpp>
#include <GComponent/GNumerical.hpp>

#include <Eigen/Geometry>

constexpr static const float kEpsilon = 1e-5;
/*
*  
*				Methods Implementations
*
**/
bool GComponent::IntersectSphereSphere(float radius_A, const Vec3f& trans_A, float radius_B, const Vec3f& trans_B)
{
	return (trans_A - trans_B).norm() + kEpsilon < (radius_A + radius_B);
}

bool GComponent::IntersectSphereCapsule(float radius_sphere, const Vec3f& trans_sphere, float radius_capsule, float half_height_capsule, const Vec3f& trans_cap, const Vec3f rot_cap)
{	
	return IntersectSphereCapsule(radius_sphere, trans_sphere,
								  radius_capsule, half_height_capsule, trans_cap, Roderigues(rot_cap));
}

bool GComponent::IntersectSphereCapsule(float radius_sphere, const Vec3f& trans_sphere, float radius_capsule, float half_height_capsule, const Vec3f& trans_cap, const SO3f& rot_cap)
{
	Vec3f half_vector = half_height_capsule * rot_cap * Vec3f::UnitZ();
	float safe_distance = radius_sphere + radius_capsule;
	return SqrDistPointSeg(trans_sphere, trans_cap - half_vector, trans_cap + half_vector) <= safe_distance * safe_distance;;
}

uint32_t GComponent::PointOutsideOfPlane4(const Vec3f& a, const Vec3f& b, const Vec3f& c, const Vec3f& d)
{
	uint32_t bit_map = 0;	// effective component is first 4 bit |-----| w | z | y | x |
							//						bit count:	  32   5  4   3   2   1
	const Vec3f ab = b - a,
				ac = c - a,
				ad = d - a,
				bd = d - b,
				bc = c - b;

	const Vec3f v0 = ab.cross(ac),	// abc plane normal vector
				v1 = ac.cross(ad),
				v2 = ad.cross(ab),
				v3 = bd.cross(bc);

	bit_map |= ((v0.dot(a) * v0.dot(d) >= -std::numeric_limits<float>::epsilon()) << 0);		// this part indicate whether oa, od in the same dir for norm of abc plane
	bit_map |= ((v1.dot(a) * v1.dot(b) >= -std::numeric_limits<float>::epsilon()) << 1);
	bit_map |= ((v2.dot(a) * v2.dot(c) >= -std::numeric_limits<float>::epsilon()) << 2);
	bit_map |= ((v3.dot(a) * v3.dot(b) >= -std::numeric_limits<float>::epsilon()) << 3);

	return bit_map;
}

// reference from Nvidia PhysX GuIntersectionBoxBox
bool GComponent::IntersectOBBOBB(const Vec3f& half_a, const Vec3f& trans_a, const Vec3f& rot_a, const Vec3f& half_b, const Vec3f& trans_b, const Vec3f& rot_b)
{	
	return IntersectOBBOBB(half_a, trans_a, Roderigues(rot_a), half_b, trans_b, Roderigues(rot_b));
}

bool GComponent::IntersectOBBOBB(const Vec3f& half_a, const Vec3f& trans_a, const SO3f& rot_a, const Vec3f& half_b, const Vec3f& trans_b, const SO3f& rot_b)
{
	const SO3f  mat_a = rot_a, mat_b = rot_b;							// rotation matrix
	const Vec3f t_w = trans_b - trans_a;								// center bias vector
	const Vec3f t_a = mat_a.transpose() * t_w;							// A frame center vector expression
	const SO3f  mat = mat_a.transpose() * mat_b;						// A frame B rot expression => R_B^A
		
	SO3f		mat_abs;
	float		ra, rb, t;
	
	// precaculate absolute value of mat
	for (int i = 0; i < 3; ++i) 
	for (int j = 0; j < 3; ++j) {
		mat_abs(i, j) = fabs(mat(i, j)) + 1e-6f;
	}

	// check A axis in A frame
	// 在以 A 为参考系的坐标下检测 A 的三个轴向
	for (int i = 0; i < 3; ++i) {
		ra = half_a(i);			
		rb = half_b.dot(mat_abs.row(i));
		t  = fabs(t_a(i));
		if (t > ra + rb) return false;
	}

	// check B axis in B frame
	// 在以 B 为参考系的坐标下检测 B 的三个轴向
	for (int i = 0; i < 3; ++i) {
		ra = half_a.dot(mat_abs.col(i));
		rb = half_b(i);
		t = fabs(mat.col(i).dot(t_a));
		if (t > ra + rb) return false;
	};

	// 9 cross product 
	// t = Ax x Bx
	// 以 Ax 轴与 Bx 轴为纵深展开的方向向量，其它的同理
	ra = half_a(1) * mat_abs(2, 0) + half_a(2) * mat_abs(1, 0);
	rb = half_b(1) * mat_abs(0, 2) + half_b(2) * mat_abs(0, 1);
	t  = fabs(t_a(1) * mat(2, 0) - t_a(2) * mat(1, 0));
	if (t > ra + rb) return false;

	// t = Ax x By
	ra = half_a(1) * mat_abs(2, 1) + half_a(2) * mat_abs(1, 1);
	rb = half_b(0) * mat_abs(0, 2) + half_b(2) * mat_abs(0, 0);
	t  = fabs(t_a(1) * mat(2, 1) - t_a(2) * mat(1, 1));
	if (t > ra + rb) return false;

	// t = Ax x Bz
	ra = half_a(1) * mat_abs(2, 2) + half_a(2) * mat_abs(1, 2);
	rb = half_b(0) * mat_abs(0, 1) + half_b(1) * mat_abs(0, 0);
	t  = fabs(t_a(1) * mat(2, 2) - t_a(2) * mat(1, 2));
	if (t > ra + rb) return false;

	// t = Ay x Bx
	ra = half_a(0) * mat_abs(2, 0) + half_a(2) * mat_abs(0, 0);
	rb = half_b(1) * mat_abs(1, 2) + half_b(2) * mat_abs(1, 1);
	t  = fabs(t_a(0) * mat(2, 0) - t_a(2) * mat(0, 0));
	if (t > ra + rb) return false;

	// t = Ay x By
	ra = half_a(0) * mat_abs(2, 1) + half_a(2) * mat_abs(0, 1);
	rb = half_b(0) * mat_abs(1, 2) + half_b(2) * mat_abs(1, 0);
	t  = fabs(t_a(0) * mat(2, 1) - t_a(2) * mat(0, 1));
	if (t > ra + rb) return false;

	// t = Ay x Bz
	ra = half_a(0) * mat_abs(2, 2) + half_a(2) * mat_abs(0, 2);
	rb = half_b(0) * mat_abs(1, 1) + half_b(1) * mat_abs(1, 0);
	t  = fabs(t_a(0) * mat(2, 2) - t_a(2) * mat(0, 2));
	if (t > ra + rb) return false;

	// t = Az x Bx
	ra = half_a(0) * mat_abs(1, 0) + half_a(1) * mat_abs(0, 0);
	rb = half_b(1) * mat_abs(2, 2) + half_b(2) * mat_abs(2, 1);
	t  = fabs(t_a(0) * mat(1, 0) - t_a(1) * mat(0, 0));
	if (t > ra + rb) return false;

	// t = Az x By
	ra = half_a(0) * mat_abs(1, 1) + half_a(1) * mat_abs(0, 1);
	rb = half_b(0) * mat_abs(2, 2) + half_b(2) * mat_abs(2, 0);
	t  = fabs(t_a(0) * mat(1, 1) - t_a(1) * mat(0, 1));
	if (t > ra + rb) return false;

	// t = Az x Bz
	ra = half_a(0) * mat_abs(1, 2) + half_a(1) * mat_abs(0, 2);
	rb = half_b(0) * mat_abs(2, 1) + half_b(1) * mat_abs(2, 0);
	t  = fabs(t_a(0) * mat(1, 2) - t_a(1) * mat(0, 2));
	if (t > ra + rb) return false;

	return true;
}

bool GComponent::IntersectOBBSphere(const Vec3f& half_box, const Vec3f& trans_box, const Vec3f& rot_box, float radius_sphere, const Vec3f& trans_sphere)
{			
	return IntersectOBBSphere(half_box, trans_box, Roderigues(rot_box), radius_sphere, trans_sphere);
}

bool GComponent::IntersectOBBSphere(const Vec3f& half_box, const Vec3f& trans_box, const SO3f& rot_box, float radius_sphere, const Vec3f& trans_sphere)
{
	const Vec3f t_w	   = trans_sphere - trans_box;
	const Vec3f t_box  = rot_box.transpose() * t_w;
	Vec3f t_align      = t_box;
	bool outside = false;

	// check if sphere origin is in the box
	// 检查球心是否在盒子内部并找出至圆心最近点
	// +------+
	// |	  |			 p1 align
	// |   *o |			 p2 ori
	// |   |--*1-----*2	
	// +------+			 
	// set the outside axis align to border
	if (t_box.x() < -half_box.x()) {
		t_align.x() = -half_box.x();
		outside = true;
	}
	else if (t_box.x() > half_box.x()) {
		t_align.x() = half_box.x();
		outside = true;
	}

	if (t_box.y() < -half_box.y()) {
		t_align.y() = -half_box.y();
		outside = true;
	}
	else if (t_box.y() > half_box.y()) {
		t_align.y() = half_box.y();
		outside = true;
	}

	if (t_box.z() < -half_box.z()) {
		t_align.z() = -half_box.z();
		outside = true;
	}
	else if (t_box.z() > half_box.z()) {
		t_align.z() = half_box.z();
		outside = true;
	}

	if (outside) {		
		// caculate the distance in sphere radius
		if (radius_sphere < (t_box - t_align).norm()) {
			return false;
		}
	}
	return true;
}

bool GComponent::IntersectOBBCapsule(const Vec3f& half_box, const Vec3f& trans_box, const Vec3f& rot_box, 
									  float radius, float half_height, const Vec3f& trans_cap, const Vec3f& rot_cap)
{
	Vec3f half_vector = half_height * Roderigues(rot_cap) * Vec3f::UnitZ();
	return GComponent::SqrDistBoxSeg(half_box, trans_box, rot_box, trans_cap - half_vector, trans_cap + half_vector) <= radius * radius;
}

bool GComponent::IntersectOBBCapsule(const Vec3f& half_box, const Vec3f& trans_box, const SO3f& rot_box, 
									 float radius, float half_height, const Vec3f& trans_cap, const SO3f& rot_cap)
{
	Vec3f half_vector = half_height * rot_cap * Vec3f::UnitZ();
	return GComponent::SqrDistBoxSeg(half_box, trans_box, LogMapSO3Toso3(rot_box), 
									 trans_cap - half_vector, trans_cap + half_vector) <= radius * radius;
}

bool GComponent::IntersectCapsuleCapsule(float radius_a, float half_height_a, const Vec3f& trans_a, const SO3f& rot_a, float radius_b, float half_height_b, const Vec3f& trans_b, const SO3f& rot_b)
{
	Vec3f half_vector_a = half_height_a * rot_a * Vec3f::UnitZ();
	Vec3f half_vector_b = half_height_b * rot_b * Vec3f::UnitZ();
	float safe_dist = radius_a + radius_b;
	return SqrDistSegSeg(trans_a - half_vector_a, trans_a + half_vector_a, 
						 trans_b - half_vector_b, trans_b + half_vector_b) <= safe_dist * safe_dist;
}

