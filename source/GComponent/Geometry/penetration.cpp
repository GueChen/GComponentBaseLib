#include "GComponent/Geometry/gcollision_detection.h"
#include "GComponent/Geometry/gdistance.h"

#include "GComponent/gtransform.hpp"
#include "GComponent/GGeometry.hpp"

namespace GComponent{

bool GComponent::PenetrationSphereSphere(PenetrationOutput& output, float radius_A, const Vec3f& trans_A, float radius_B, const Vec3f& trans_B)
{
	Vec3f diff	   = trans_B - trans_A,
		  diff_dir = diff.normalized();

	output.normal	 = diff_dir;
	output.depth	 = radius_A + radius_B - diff.norm();
	output.closest_a = trans_A + radius_A * diff_dir;
	output.closest_b = trans_B - radius_B * diff_dir;

	return output.depth >= 0;
}

bool PenetrationSphereCapsule(PenetrationOutput& output, float radius_sphere, const Vec3f& trans_sphere, float radius_capsule, float half_height_capsule, const Vec3f& trans_cap, const SO3f& rot_cap)
{
	Vec3f half_vec    = half_height_capsule * rot_cap * Vec3f::UnitZ();	
	Vec3f seg_a	      = trans_cap - half_vec, seg_b = trans_cap + half_vec;
	Vec3f closest_seg = seg_a;

	// get the closest point on segment and sqr distance
	float sqr = SqrDistPointSeg(trans_sphere, seg_a, seg_b, &closest_seg);

	float marg_sum = radius_sphere + radius_capsule;
	
	Vec3f diff     = closest_seg - trans_sphere,
		  diff_dir = diff.normalized();

	output.normal	 = diff_dir;
	output.depth	 = marg_sum - diff.norm();
	output.closest_a = trans_sphere + radius_sphere * diff_dir;
	output.closest_b = closest_seg - radius_capsule * diff_dir;

	return output.depth >= 0;
}

static bool BoxBoxAxisTest(PenetrationOutput& output, 
						   const Vec3f& half_A, const Vec3f& trans_A, const SO3f& rot_A,
						   const Vec3f& half_B, const Vec3f& trans_B, const SO3f& rot_B,
						   const Vec3f& axis)
{
	const float center_A = trans_A.dot(axis),
				extend_A = (rot_A * axis).cwiseAbs().dot(half_A);
	const float center_B = trans_B.dot(axis),
				extend_B = (rot_B * axis).cwiseAbs().dot(half_B);
	const float s_A = center_A - extend_A, e_A = center_A + extend_A;
	const float s_B = center_B - extend_B, e_B = center_B + extend_B;
	
	if (e_A < s_B || e_B < s_A) {
		output.depth = 0.0f; // replace with a negative value
		return false;
	}
	
	const float dist = std::min(e_B - s_A, e_A - s_B);

	if (dist < output.depth) {
		output.normal = axis;
		output.depth  = dist;
	}
	return true;
}

bool PenetrationOBBOBB(PenetrationOutput& output, 
					   const Vec3f& half_A, const Vec3f& trans_A, const SO3f& rot_A, 
					   const Vec3f& half_B, const Vec3f& trans_B, const SO3f& rot_B)
{	
	output.depth = std::numeric_limits<float>::max();
	// test all A axis 
	for (auto & axis_A : rot_A.colwise())
	if (not BoxBoxAxisTest(output,
	    half_A, trans_A, rot_A,
	    half_B, trans_B, rot_B,
	    axis_A)){
		return false;
	}

	// test all B axis
	for (auto& axis_B : rot_B.colwise()) 
	if(not BoxBoxAxisTest(output,
	   half_A, trans_A, rot_A,
	   half_B, trans_B,  rot_B,
	   axis_B)){
		return false;
	}

	// test cross axis
	for (auto& axis_A : rot_A.colwise()) for (auto& axis_B : rot_B.colwise()) {
		Vec3f cross = axis_A.cross(axis_B);		
		if (cross.squaredNorm() > 1e-5) {			
			cross = cross.normalized();
			if (not BoxBoxAxisTest(output,
				half_A, trans_A, rot_A,
				half_B, trans_B, rot_B,
				cross)) {
				return false;
			}
		}
	}

	if (output.normal.dot(trans_A - trans_B) >= 0.0f) {
		output.normal = -output.normal;
	}

	return true;
}

bool PenetrationOBBSphere(PenetrationOutput& output, const Vec3f& half_box, const Vec3f& trans_box, const SO3f& rot_box, float radius_sphere, const Vec3f& trans_sphere)
{	
	Vec3f closest_on_box;
	float sqr = SqrDistBoxPoint(half_box, trans_box, rot_box, trans_sphere, &closest_on_box);
	// outside box
	if (sqr > 0) {
		Vec3f diff		 = closest_on_box - trans_sphere,
			  diff_dir	 = diff.normalized();

		output.normal	 = diff_dir;
		output.depth	 = radius_sphere - diff.norm();
		output.closest_a = trans_sphere - radius_sphere * diff_dir;
		output.closest_b = closest_on_box;

		return output.depth >= 0;
	}
	// inside box find the closest face
	else {		
		const Vec3f t_w = trans_sphere - trans_box;		// relative box to sphere trans world
		const Vec3f t_l = rot_box.transpose() * t_w;	// rel trans in box local coordinate	

		// find closest surface
		const Vec3f dist_to_suf = half_box - t_l.cwiseAbs();		
		int   min_idx = 0;
		float min_dist = dist_to_suf.x();
		for (int i = 0; i < 3; ++i) if(dist_to_suf(i) < min_dist){
			min_idx  = i;
			min_dist = dist_to_suf(i);
		}

		// caculate box closest point in box local coordinate
		const Vec3f sign_mask = t_l.cwiseSign();
		Vec3f diff_dir  = rot_box * Vec3f::Unit(min_idx).cwiseProduct(sign_mask);
		Vec3f box_closest	  = t_l;
		box_closest(min_idx)  = sign_mask(min_idx) * half_box(min_idx);
		if (diff_dir.norm() < 1e-5) {
			diff_dir = Vec3f::Unit(min_idx);
		}

		output.normal = diff_dir;
		output.depth  = -min_dist - radius_sphere;
		output.closest_a = trans_sphere + radius_sphere * diff_dir;
		output.closest_b = trans_box + rot_box * box_closest;

		return true;
	}
	
}

static bool BoxCapsuleAxisTest(PenetrationOutput& output,
							   const Vec3f& half_box, const Vec3f& trans_box, const SO3f& rot_box, 
							   const Vec3f& p_a, const Vec3f& p_b, float radius,
							   const Vec3f& axis) {
	// make capsule project on axis
	float s_seg = p_a.dot(axis), e_seg = p_b.dot(axis);
	if (s_seg > e_seg) std::swap(s_seg, e_seg);
	s_seg -= radius, e_seg += radius;

	// make box project on axis
	const float center_box = trans_box.dot(axis),
				extend_box = (rot_box * axis).cwiseAbs().dot(half_box);
	const float s_box = center_box - extend_box, e_box = center_box + extend_box;

	if (e_seg < s_box || e_box < s_seg) return false;

	float depth = std::min(e_box - s_seg, e_seg - s_box);
	if (depth < output.depth) {
		output.depth  = depth;
		output.normal = axis;
	}
	return true;
}

static bool ComputeMTDOBBCapsule(PenetrationOutput& output, const Vec3f& half_box, const Vec3f& trans_box, const SO3f& rot_box, const Vec3f& p_a, const Vec3f& p_b, float radius) {	
	output.depth = std::numeric_limits<float>::max();

	// check the box axis
	for (auto& axis_box : rot_box.colwise()) {
		if (not BoxCapsuleAxisTest(output, 
			half_box, trans_box, rot_box, 
			p_a, p_b, radius, 
			axis_box)) {
			output.depth = 0.0f;
			return false;
		}		
	}

	// check the capsule box edge combines
	const Vec3f cap_seg_dir = (p_b - p_a).normalized();
	for (auto& axis_box : rot_box.colwise()) {
		Vec3f cross = cap_seg_dir.cross(axis_box);
		if (cross.squaredNorm() > 1e-5f) {
			cross = cross.normalized();
			if (not BoxCapsuleAxisTest(output,
				half_box, trans_box, rot_box,
				p_a, p_b, radius,
				axis_box)) {
				output.depth = 0.0f;
				return false;
			}
		}
	}

	if (output.normal.dot(0.5f * (p_a + p_b) - trans_box) >= 0.0f) {
		output.normal = -output.normal;
	}
	return true;
}

bool PenetrationOBBCapsule(PenetrationOutput& output, const Vec3f& half_box, const Vec3f& trans_box, const SO3f& rot_box, float radius, float half_height, const Vec3f& trans_cap, const SO3f& rot_cap)
{
	const Vec3f half_vec = half_height * rot_cap * Vec3f::UnitZ();
	const Vec3f p_a = trans_cap - half_vec,
				p_b = trans_cap + half_vec;
	Vec3f on_seg, on_box;
	const float sq_d = SqrDistBoxSeg(half_box, trans_box, rot_box, 
									 p_a, p_b, 
									 &on_box, &on_seg);		
	// capsule seg not in box	
	//	 +-----+    
	//  /box  /|   |
	// +-----+*+   * on_seg
	// |     |/    |
	// +-----+
	if (sq_d != 0) {
		Vec3f diff = on_box - on_seg,
			diff_dir = diff.normalized();
		output.normal = diff_dir;
		output.depth  = radius - sqrtf(sq_d);
		output.closest_a = on_seg + radius * diff_dir;
		output.closest_b = on_box;

		return output.depth >= 0.0f;
	}
	// capsule seg in box
	//     +-----+	
	//	  /     /|
	//	 +-----+ |		seg in box 
	// --|     |--------  dist = 0
	//   |     | +
	//	 |     |/ 
	//	 +-----+
	return  ComputeMTDOBBCapsule(output, half_box, trans_box, rot_box,
								 p_a, p_b, radius);
}

bool PenetrationCapsuleCapsule(PenetrationOutput& output, float radius_a, float half_height_a, const Vec3f& trans_a, const SO3f& rot_a, float radius_b, float half_height_b, const Vec3f& trans_b, const SO3f& rot_b)
{	
	const Vec3f half_a = half_height_a * rot_a * Vec3f::UnitZ(),
				half_b = half_height_b * rot_b * Vec3f::UnitZ();
	
	const Vec3f p_a0 = trans_a - half_a, p_a1 = trans_a + half_a,
				p_b0 = trans_b - half_b, p_b1 = trans_b + half_b;

	Vec3f on_seg_a, on_seg_b;
	const float sq_d = SqrDistSegSeg(p_a0, p_a1, p_b0, p_b1, &on_seg_a, &on_seg_b);

	Vec3f diff = on_seg_b - on_seg_a,
		  diff_dir = diff.normalized();
			
	output.normal = diff_dir;
	output.depth  = radius_a + radius_b - sqrtf(sq_d);
	output.closest_a = on_seg_a + radius_a * diff_dir;
	output.closest_b = on_seg_b - radius_b * diff_dir;

	return output.depth >= 0.0f;
}


/*___________________________________Box Contact Related Method____________________________________________*/
/**
* 
* now just simple copy from NVIDIA/Omniverse/PhysX5/GuContactBoxBox.cpp
* TODO:
*	+ use better abstract to encapsule the invoke interface
*   + now for caculation exposing too many detatils
* 
* */

struct ContactPoint {
	Vec3f point;
	Vec3f normal;
	float seperation;

	ContactPoint(const Vec3f& p, const Vec3f& n, float sep) :
		point(p), normal(n), seperation(sep)
	{}
};
using PenetrationPoints = std::vector<ContactPoint>;


static float InYZ(float y, float z, const std::array<Vec3f*, 4>& face) {
	
	float pre_y = face[3]->y(), pre_z = face[3]->z();
	for (auto& vert : face) {	// judge the (y, z) is the inner point
		const float cur_y = vert->y(), cur_z = vert->z();
		if ((cur_y - pre_y) * (z - pre_z) - (cur_z - pre_z) * (y - pre_y) >= 0) return -1.0;
		pre_y = cur_y, pre_z = cur_z;
	}
	float x  = face[0]->x();
	
	float ay = y - face[0]->y(), az = z - face[0]->z();	
	Vec3f b = *face[1] - *face[0];
	x += b.x() * (ay * b.y() + az * b.z()) / b.squaredNorm();

	b = *face[3] - *face[0];
	x += b.x() * (ay * b.y() + az * b.z()) / b.squaredNorm();

	return x;
}

static bool ContactBoxPlane(PenetrationPoints& contacts, const Vec3f& normal, float y_border, float z_border, const Vec3f& half, const SE3f& trans_0, const SE3f& trans_1) {
	const SE3f trans_1to0	  = InverseSE3(trans_0) * trans_1;
	auto [rot_1to0, pos_1to0] = RtDecompositionMat4(trans_1to0);

	struct VertexInfo {
		bool in_area   = false;
		bool penetrate = false;
	};
	/*	cube motioned reference from GuContactBoxBox.cpp
	       6 +------+ 7
	        /|     /|
	       / |    / |
	      / 4+---/--+5
	    2+------+3 /    y  z
	     | /    | /     | /
	     |/     |/      |/
	    0+------+1      *---x
	*/

	std::array<Vec3f, 8>	  pos;
	std::array<VertexInfo, 8> infos;
	{
		std::array<Vec3f, 3> es;
		for (int i = 0; i < 3; ++i) es[i] = rot_1to0.col(i) * half(i);		
		pos[0] = pos_1to0 - es[0] - es[1] - es[2];
		pos[1] = pos_1to0 + es[0] - es[1] - es[2];
		pos[2] = pos_1to0 - es[0] + es[1] - es[2];
		pos[3] = pos_1to0 + es[0] + es[1] - es[2];
		pos[4] = pos_1to0 - es[0] - es[1] + es[2];
		pos[5] = pos_1to0 + es[0] - es[1] + es[2];
		pos[6] = pos_1to0 - es[0] + es[1] + es[2];
		pos[7] = pos_1to0 + es[0] + es[1] + es[2];		
	}

	// check foreach whether point in region 
	for (uint32_t i = 0; i < 8; i++) {
		if (pos[i].x() < 0) continue;				
		infos[i].penetrate = true;
		// the vertex which in [0,+) [-y,y] [-z,z]
		if (fabs(pos[i].y()) <= y_border && 
			fabs(pos[i].z()) <= z_border) {
			infos[i].in_area = true;
			contacts.emplace_back(pos[i], normal, -pos[i].x());			
		}
	}

	constexpr static const uint32_t edges[] = {
		0,1, 1,3, 3,2, 2,0, 4,5, 5,7, 7,6, 6,4, 0,4, 1,5, 2,6, 3,7
	};
	for (auto ptr = std::begin(edges); ptr != std::end(edges); ptr += 2) {
		uint32_t p1_idx = *ptr, p2_idx = *(ptr + 1);
		VertexInfo&i1 = infos[p1_idx], &i2 = infos[p2_idx];
		Vec3f& p1 = pos[p1_idx], &p2 = pos[p2_idx];
		if (i1.penetrate || i2.penetrate) {
			// edge cross from the quad
			//     quad
			//	+--------+ 
			//  |    p1+-|--+ p2	edge
			//  |     /  |   \
			//  |in area |    not in area
			//  +--------+
			if (!i1.in_area || !i2.in_area) {
				
				// Test whether intesect on axis y
				if (p1.y() > p2.y()) std::swap(p1, p2);
				
				// Test +Y border
				if (p1.y() < y_border && p2.y() >= y_border) { 
					float factor  = (y_border - p1.y()) / (p2.y() - p1.y());
					float z = Lerp(p1.z(), p2.z(), factor);					
					float x = Lerp(p1.x(), p2.x(), factor);
					if (fabs(z) <= z_border && x >= 0.0f) {	// judge the intesect point whether in [-z, z]						
						contacts.emplace_back(Vec3f(x, y_border, z), normal, -x);						
					}
				}
				
				// Test -Y border
				if (p1.y() < -y_border && p2.y() >= -y_border) {
					float factor = (-y_border - p1.y()) / (p2.y() - p1.y());
					float z = Lerp(p1.z(), p2.z(), factor);
					float x = Lerp(p1.x(), p2.x(), factor);
					if (fabs(z) <= z_border && x >= 0.0f) {
						contacts.emplace_back(Vec3f(x, -y_border, z), normal, -x);						
					}
				}

				if (p1.z() > p2.z()) std::swap(p1, p2);
				if (p1.z() < z_border && p2.z() >= z_border) {	// TEST +Z				
					float factor = (z_border - p1.z()) / (p2.z() - p1.z());
					float y = Lerp(p1.y(), p2.y(), factor);
					float x = Lerp(p1.x(), p2.x(), factor);
					if (fabs(y) <= y_border && x >= 0.0f) {						
						contacts.emplace_back(Vec3f(x, y, z_border), normal, -x);						
					}
				}

				if (p1.z() < -z_border && p2.z() >= -z_border) {	// TEST -Z				
					float factor = (-z_border - p1.z()) / (p2.z() - p1.z());
					float y = Lerp(p1.y(), p2.y(), factor);
					float x = Lerp(p1.x(), p2.x(), factor);
					if (fabs(y) <= y_border && x >= 0.0f) {							
						contacts.emplace_back(Vec3f(x, y, -z_border), normal, -x);
					}
				}
			}

			if ((!i1.penetrate && !i2.in_area) || (!i1.in_area && !i2.penetrate)) {
				float factor = (-p1.x()) / (p2.x() - p1.x());
				float y = Lerp(p1.y(), p2.y(), factor);
				float z = Lerp(p1.z(), p2.z(), factor);
				if (fabs(y) <= y_border && fabs(z) <= z_border) {
					contacts.emplace_back(Vec3f(0, y, z), normal, 0);
				}
			}
		}
	}

	constexpr static const uint32_t faces[][4] = {
		{0,1,3,2}, {1,5,7,3}, {5,4,6,7}, {4,0,2,6}, {2,3,7,6}, {0,4,5,1}
	};
	int add = 0;
	for (auto ptr = std::begin(faces); ptr != std::end(faces) && add != 0x0f; ++ptr) {
		const uint32_t* face_p = *ptr;
		std::array<VertexInfo, 4> is = {
			infos[face_p[0]], infos[face_p[1]], infos[face_p[2]], infos[face_p[3]]
		};
		std::array<Vec3f*, 4> vs = {
			&pos[face_p[0]], &pos[face_p[1]], &pos[face_p[2]], &pos[face_p[3]]
		};
		if (is[0].penetrate && is[1].penetrate && is[2].penetrate && is[3].penetrate) {
			if (!is[0].in_area || !is[1].in_area || !is[2].in_area || !is[3].in_area) {
				if (!(add & 1)) if (float x = InYZ(-y_border, -z_border, vs); x >= 0) {
					add |= 1; contacts.emplace_back(Vec3f(x, -y_border, -z_border), normal, -x);						
				}				
				if (!(add & 2)) if (float x = InYZ(y_border, -z_border, vs); x >= 0) {
					add |= 2; contacts.emplace_back(Vec3f(x, y_border, -z_border), normal, -x);
				}				
				if (!(add & 4)) if (float x = InYZ(-y_border, z_border, vs); x >= 0) {
					add |= 4; contacts.emplace_back(Vec3f(x, -y_border, z_border), normal, -x);
				}				
				if (!(add & 8)) if (float x = InYZ(y_border, z_border, vs); x >= 0) {
					add |= 8; contacts.emplace_back(Vec3f(x, y_border, z_border), normal, -x);					
				}
			}
		}				
	}
	
	// transfer from local to world
	for (auto& p : contacts) {
		p.point = AffineProduct(trans_0, p.point);
	}

	return not contacts.empty();

}

static bool ContactOBBOBB(PenetrationOutput& output, const Vec3f& half_A, const Vec3f& trans_A, const SO3f& rot_A, const Vec3f& half_B, const Vec3f& trans_B, const SO3f& rot_B)
{
	Vec3f t		 = trans_B - trans_A;				// local trans centriod B from A
	SO3f  B_to_A = rot_A.transpose() * rot_B;		// trans box point from local to A
	SO3f  abs_B_to_A = B_to_A.cwiseAbs();

	Vec3f t_A	  = rot_A.transpose() * t;
	Vec3f abs_t_A = t_A.cwiseAbs();
	
	std::array<float, 6> test_axis_d;
	auto axis_d_cur = test_axis_d.begin();

	// test A standard axis
	{
		Vec3f proj_to_A = abs_B_to_A * half_B;
		for (int i = 0; i < 3; ++i) {			
			*axis_d_cur = half_A(i) + proj_to_A(i) - abs_t_A(i);
			if (*axis_d_cur++ < 0) {
				// seperate in A axis{i}
				return false;
			}
		}
	}

	Vec3f t_B = rot_B.transpose() * t;
	Vec3f abs_t_B = t_B.cwiseAbs();

	// test B standard axis
	{
		Vec3f proj_to_B = abs_B_to_A.transpose() * half_A;
		for (int i = 0; i < 3; ++i) {
			*axis_d_cur = half_B(i) + proj_to_B(i) - abs_t_B(i);
			if (*axis_d_cur++ < 0) {
				// seperate in B axis{i}
				return false;
			}
		}
	}

	// test cross axis see also 
	{
		float t, ra, rb;
		// t = Ax x Bx
		ra = half_A(1) * abs_B_to_A(2, 0) + half_A(2) * abs_B_to_A(1, 0);
		rb = half_B(1) * abs_B_to_A(0, 2) + half_B(2) * abs_B_to_A(0, 1);
		t = fabs(t_A(2) * B_to_A(1, 0) - t_A(1) * B_to_A(2, 0));
		if (t > ra + rb) return false;

		// t = Ax x By
		ra = half_A(1) * abs_B_to_A(2, 1) + half_A(2) * abs_B_to_A(1, 1);
		rb = half_B(0) * abs_B_to_A(0, 2) + half_B(2) * abs_B_to_A(0, 0);
		t = fabs(t_A(1) * B_to_A(2, 1) - t_A(2) * B_to_A(1, 1));
		if (t > ra + rb) return false;

		// t = Ax x Bz
		ra = half_A(1) * abs_B_to_A(2, 2) + half_A(2) * abs_B_to_A(1, 2);
		rb = half_B(0) * abs_B_to_A(0, 1) + half_B(1) * abs_B_to_A(0, 0);
		t = fabs(t_A(1) * B_to_A(2, 2) - t_A(2) * B_to_A(1, 2));
		if (t > ra + rb) return false;

		// t = Ay x Bx
		ra = half_A(0) * abs_B_to_A(2, 0) + half_A(2) * abs_B_to_A(0, 0);
		rb = half_B(1) * abs_B_to_A(1, 2) + half_B(2) * abs_B_to_A(1, 1);
		t = fabs(t_A(0) * B_to_A(2, 0) - t_A(2) * B_to_A(0, 0));
		if (t > ra + rb) return false;

		// t = Ay x By
		ra = half_A(0) * abs_B_to_A(2, 1) + half_A(2) * abs_B_to_A(0, 1);
		rb = half_B(0) * abs_B_to_A(1, 2) + half_B(2) * abs_B_to_A(1, 0);
		t = fabs(t_A(0) * B_to_A(2, 1) - t_A(2) * B_to_A(0, 1));
		if (t > ra + rb) return false;

		// t = Ay x Bz
		ra = half_A(0) * abs_B_to_A(2, 2) + half_A(2) * abs_B_to_A(0, 2);
		rb = half_B(0) * abs_B_to_A(1, 1) + half_B(1) * abs_B_to_A(1, 0);
		t = fabs(t_A(0) * B_to_A(2, 2) - t_A(2) * B_to_A(0, 2));
		if (t > ra + rb) return false;

		// t = Az x Bx
		ra = half_A(0) * abs_B_to_A(1, 0) + half_A(1) * abs_B_to_A(0, 0);
		rb = half_B(1) * abs_B_to_A(2, 2) + half_B(2) * abs_B_to_A(2, 1);
		t = fabs(t_A(0) * B_to_A(1, 0) - t_A(1) * B_to_A(0, 0));
		if (t > ra + rb) return false;

		// t = Az x By
		ra = half_A(0) * abs_B_to_A(1, 1) + half_A(1) * abs_B_to_A(0, 1);
		rb = half_B(0) * abs_B_to_A(2, 2) + half_B(2) * abs_B_to_A(2, 0);
		t = fabs(t_A(0) * B_to_A(1, 1) - t_A(1) * B_to_A(0, 1));
		if (t > ra + rb) return false;

		// t = Az x Bz
		ra = half_A(0) * abs_B_to_A(1, 2) + half_A(1) * abs_B_to_A(0, 2);
		rb = half_B(0) * abs_B_to_A(2, 1) + half_B(1) * abs_B_to_A(2, 0);
		t = fabs(t_A(0) * B_to_A(1, 2) - t_A(1) * B_to_A(0, 2));
		if (t > ra + rb) return false;
	}

	// Find which axis is the minmum sperate axis
	float min_dist = std::numeric_limits<float>::max();
	int   min_idx  = -1;
	for (int i = 0; i < 6; ++i) {
		float d = test_axis_d[i];
		if (d >= 0.0f && d < min_dist) {
			min_idx = i, min_dist = d;
		}
	}

	enum MinAxis { A_x = 0, A_y, A_z, B_x, B_y, B_z }
	which = static_cast<MinAxis>(min_idx);

	SE3f  trans_0 = SE3f::Identity(), 
		  trans_1 = SE3f::Identity();
	Vec3f axis;	
	std::vector<ContactPoint> points;

	if (which <= A_z) {
		trans_1.block(0, 0, 3, 3) = rot_B;
		trans_1.block(0, 3, 3, 1) = trans_B;
	}
	else {
		trans_1.block(0, 0, 3, 3) = rot_A;
		trans_1.block(0, 3, 3, 1) = trans_A;
	}

	switch (which) {
	case A_x:
		axis = rot_A.col(0);
		trans_0.block(0, 0, 3, 3) = rot_A;
		if (t_A(which) < 0) {					
			axis = -axis;
			trans_0.block(0, 0, 3, 2) = -trans_0.block(0, 0, 3, 2);
		}
		trans_0.block(0, 3, 3, 1) = trans_A - half_A(which) * axis;				
		ContactBoxPlane(points, axis, half_A.y(), half_A.z(), half_B, trans_0, trans_1);
		break;
	case A_y:
		axis = rot_A.col(1);
		trans_0.block(0, 0, 3, 1) = rot_A.col(1);
		trans_0.block(0, 1, 3, 1) = rot_A.col(2);
		trans_0.block(0, 2, 3, 1) = rot_A.col(0);
		if (t_A.y() < 0) {
			axis = -axis;
			trans_0.block(0, 0, 3, 2) = -trans_0.block(0, 0, 3, 2);
		}
		trans_0.block(0, 3, 3, 1) = trans_A - half_A.y() * axis;
		ContactBoxPlane(points, axis, half_A.z(), half_A.x(), half_B, trans_0, trans_1);
		break;
	case A_z:
		axis = rot_A.col(2);
		trans_0.block(0, 0, 3, 1) = rot_A.col(2);
		trans_0.block(0, 1, 3, 1) = rot_A.col(0);
		trans_0.block(0, 2, 3, 1) = rot_A.col(1);
		if (t_A.z() < 0) {
			axis = -axis;
			trans_0.block(0, 0, 3, 2) = -trans_0.block(0, 0, 3, 2);
		}
		trans_0.block(0, 3, 3, 1) = trans_A - half_A.z() * axis;
		ContactBoxPlane(points, axis, half_A.x(), half_A.y(), half_B, trans_0, trans_1);
		break;
	case B_x:
		axis = -rot_B.col(0);
		trans_0.block(0, 0, 3, 3) = rot_B;
		if (t_B.x() >= 0) {
			axis = -axis;
			trans_0.block(0, 0, 3, 2) = -trans_0.block(0, 0, 3, 2);
		}
		trans_0.block(0, 3, 3, 1) = trans_B + half_B.x() * axis;
		ContactBoxPlane(points, axis, half_B.z(), half_B.x(), half_A, trans_0, trans_1);
		break;
	case B_y:
		axis = -rot_B.col(1);
		trans_0.block(0, 0, 3, 1) = rot_B.col(1);
		trans_0.block(0, 1, 3, 1) = rot_B.col(2);
		trans_0.block(0, 2, 3, 1) = rot_B.col(0);
		if (t_B.y() >= 0) {
			axis = -axis;
			trans_0.block(0, 0, 3, 2) = -trans_0.block(0, 0, 3, 2);
		}
		trans_0.block(0, 3, 3, 1) = trans_B + half_B.y() * axis;
		ContactBoxPlane(points, axis, half_B.z(), half_B.x(), half_A, trans_0, trans_1);
		break;
	case B_z:
		axis = -rot_B.col(2);
		trans_0.block(0, 0, 3, 1) = rot_B.col(2);
		trans_0.block(0, 1, 3, 1) = rot_B.col(0);
		trans_0.block(0, 2, 3, 1) = rot_B.col(1);
		if (t_B.z() >= 0) {
			axis = -axis;
			trans_0.block(0, 0, 3, 2) = -trans_0.block(0, 0, 3, 2);
		}
		trans_0.block(0, 3, 3, 1) = trans_B + half_B.z() * axis;
		ContactBoxPlane(points, axis, half_B.x(), half_B.y(), half_A, trans_0, trans_1);
		break;
	default:
		assert(false && "Never happens");
		return false;
	}

	output.normal = axis;
	output.depth = min_dist;
	// closest
	return !points.empty();
}

}