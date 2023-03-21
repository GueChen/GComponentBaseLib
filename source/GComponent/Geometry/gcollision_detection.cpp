#include "GComponent/Geometry/gcollision_detection.h"
#include "GComponent/Geometry/gdistance.h"

#include <GComponent/gtransform.hpp>
#include <GComponent/GNumerical.hpp>

#include <Eigen/Geometry>

static Vec3f NearestPointTriangle(Vec3f* Q, const uint32_t& out_plane_4, uint32_t* indices, uint32_t& size)
{
	float	 sqrt_dist = std::numeric_limits<float>::max();
	uint32_t calc_indices[3] = { 0, 1, 2 };
	Vec3f    closest = Vec3f::Zero();

	auto compare_nearest = [&](uint32_t i1, uint32_t i2, uint32_t i3) {
		calc_indices[0] = i1;
		calc_indices[1] = i2;
		calc_indices[2] = i3;

		uint32_t	calc_size = 3;
		Vec3f		calc_closest;

		const float calc_sqrt_dist = GComponent::NearestTriangleBaryCentric(Q[i1], Q[i2], Q[i3], calc_indices, calc_size, calc_closest);

		if (sqrt_dist > calc_sqrt_dist) {
			closest = calc_closest;
			sqrt_dist = calc_sqrt_dist;

			indices[0] = calc_indices[0];
			indices[1] = calc_indices[1];
			indices[2] = calc_indices[2];
			size = calc_size;
		}
	};

	if ((out_plane_4 >> 0) & 1) {	// a point out of plane case
		sqrt_dist = GComponent::NearestTriangleBaryCentric(Q[0], Q[1], Q[2], indices, size, closest);
	}

	if ((out_plane_4 >> 1) & 1) {	// b point out of plane case
		compare_nearest(0, 2, 3);
	}

	if ((out_plane_4 >> 2) & 1) {
		compare_nearest(0, 3, 1);
	}

	if ((out_plane_4 >> 3) & 1) {
		compare_nearest(1, 3, 2);
	}

	return closest;
}


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

	const Vec3f v0 = ab.cross(ac),
				v1 = ac.cross(ad),
				v2 = ad.cross(ab),
				v3 = bd.cross(bc);

	bit_map |= ((v0.dot(a) * v0.dot(d) > 0) << 0);
	bit_map |= ((v1.dot(a) * v1.dot(b) > 0) << 1);
	bit_map |= ((v2.dot(a) * v2.dot(c) > 0) << 2);
	bit_map |= ((v3.dot(a) * v3.dot(b) > 0) << 3);

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

Vec3f GComponent::NearestSimplex(Vec3f* Q, Vec3f* A, Vec3f* B, uint32_t& size, Vec3f& support)
{
	enum SimplexShape : uint32_t {
		ePoint = 1, eSegment, eTriangle, eTetrahedron
	};
	switch (size) {
	case ePoint:	   return support;
	case eSegment:	   return NearestSegment	(Q, size);
	case eTriangle:	   return NearestTriangle	(Q, A, B, size);
	case eTetrahedron: return NearestTetrahedron(Q, A, B, size);
	default:
		assert(false && "Nearest Simplex Size Error");
	}
	return support;
}

Vec3f GComponent::NearestSegment(Vec3f* Q, uint32_t& size)
{
	Vec3f &a = Q[0], &b = Q[1];
	
	Vec3f ab = b - a,
		  ao =   - a;
	
	float sq_ab_length = ab.squaredNorm();								// 长度过短退化为点
	if (sq_ab_length < std::numeric_limits<float>::epsilon()) {	// degenerate to single point
		size = 1;
		return Q[0];
	}

	float scale = Clamp(ao.dot(ab) / sq_ab_length, 0.0f, 1.0f);

	return a + scale * ab;
}

Vec3f GComponent::NearestTriangle(Vec3f* Q, Vec3f* A, Vec3f* B, uint32_t& size)
{
	Vec3f& a = Q[0], & b = Q[1], & c = Q[2];
	Vec3f ab	= b - a,
		  ac	= c - a,
		  ab_ac = ab.cross(ac);

	float area  = ab_ac.norm();
	if (area < std::numeric_limits<float>::epsilon()) {	// three point on one same line 
		size = 2;										// degenerate to segment
		return NearestSegment(Q, size);
	}

	uint32_t calc_size;
	uint32_t indices[3] = { 0, 1, 2 };
	Vec3f	 closest;

	NearestTriangleBaryCentric(a, b, c, indices, calc_size, closest);
	
	if (calc_size != 3) {
		const Vec3f q0 = Q[indices[0]], q1 = Q[indices[1]],
					a0 = A[indices[0]], a1 = A[indices[1]],
					b0 = B[indices[0]], b1 = B[indices[1]];
		
		Q[0] = q0; Q[1] = q1;
		A[0] = a0; A[1] = a1;
		B[0] = b0; B[1] = b1;

		size = calc_size;
	}

	return closest;
}

Vec3f GComponent::NearestTetrahedron(Vec3f* Q, Vec3f* A, Vec3f* B, uint32_t& size)
{
	const Vec3f a = Q[0],
				b = Q[1],
				c = Q[2],
				d = Q[3];

	const Vec3f ab = b - a,
				ac = c - a,
				ad = d - a,
				n  = ab.cross(ac).normalized();
	if (abs(n.dot(ad)) < 1e-4) {	// degenerated to triangle case
		size = 3;
		return NearestTriangle(Q, A, B, size);
	}

	bool out_plane = PointOutsideOfPlane4(a, b, c, d);
	if (out_plane) {
		return Vec3f::Zero();
	}

	uint32_t indices[3] = { 0, 1, 2 };
	const Vec3f closest = NearestPointTriangle(Q, out_plane, indices, size);

	const Vec3f q0 = Q[indices[0]], q1 = Q[indices[1]], q2 = Q[indices[2]], 
				a0 = A[indices[0]], a1 = A[indices[1]], a2 = A[indices[2]],
				b0 = B[indices[0]], b1 = B[indices[1]], b2 = B[indices[2]];

	Q[0] = q0; Q[1] = q1; Q[2] = q2;
	A[0] = a0; A[1] = a1; A[2] = a2;
	B[0] = b0; B[1] = b1; B[2] = b2;

	return closest;
}

float GComponent::NearestTriangleBaryCentric(Vec3f& a, Vec3f& b, Vec3f& c, uint32_t* indices, uint32_t& size, Vec3f& closest)
{
	size = 3;
	const Vec3f ab = b - a,
				ac = c - a,
				bc = c - b,
				n  = ab.cross(ac);
	
	const float n_sqrt = n.dot(n);
	if (n_sqrt < std::numeric_limits<float>::epsilon()) {
		return std::numeric_limits<float>::max();
	}

	const Vec3f cross_bc = b.cross(c),
				cross_ca = c.cross(a),
				cross_ab = a.cross(b);

	const float va = n.dot(cross_bc),
				vb = n.dot(cross_ca),
				vc = n.dot(cross_ab);
	
	// size = 3
	if (va > 0 && vb > 0 && vc > 0) {
		closest = n.dot(a) / n_sqrt * n;
		return closest.dot(closest);
	}

	const Vec3f ao = -a,
				bo = -b,
				co = -c;

	const float d1 = ab.dot(ao),
				d2 = ac.dot(ao),
				d3 = ab.dot(bo),
				d4 = ac.dot(bo),
				d5 = ab.dot(co),
				d6 = ac.dot(co);

	const float unom   = d4 - d3,
				udenom = d5 - d6;

	size = 2;
	if (vc <= 0 && d1 >= 0 && d3 <= 0) {		// the closest point in ab edge segment
		float recip = d1 - d3;
		if (abs(recip) > std::numeric_limits<float>::epsilon()) {
			closest = d1 * recip * ab + a; 
		}
		else {
			closest = a;
		}
		return closest.dot(closest);
	}
	else if (va <= 0 && d4 >= d3 && 0 <= d6) {	// the closest point in bc edge segment
		float recip = unom + udenom;
		if (abs(recip) > std::numeric_limits<float>::epsilon()) {
			closest = d2 / recip * bc + b;
		}
		else {
			closest = b;
		}
		indices[0] = indices[1];
		indices[1] = indices[2];
		return closest.dot(closest);
	}
	else if (vb <= 0 && d2 >= 0 && d6 <= 0) {	// the closest point in ac edge segment
		float recip = d2 - d6;
		if (recip > std::numeric_limits<float>::epsilon()) {
			closest = d2 / recip * ac + a;
		}
		else {
			closest = a;
		}
		indices[1] = indices[2];
		return closest.dot(closest);
	}

	size = 1;
	if (d1 <= 0 && d2 <= 0) {				// the closest point outside of a
		closest = a;
		
	}
	else if (d3 >= 0 && d3 >= d4) {
		closest = b;
		indices[0] = indices[1];
	}
	else {
		closest = c;
		indices[0] = indices[2];
	}

	return closest.dot(closest);
}

bool GComponent::IntersectEPA()
{
	return false;
}
