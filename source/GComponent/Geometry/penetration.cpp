#include "GComponent/Geometry/gcollision_detection.h"
#include "GComponent/Geometry/gdistance.h"

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

bool PenetrationOBBSphere(PenetrationOutput& output, const Vec3f& half_box, const Vec3f& trans_box, const SO3f& rot_box, float radius_sphere, const Vec3f& trans_sphere)
{	
	Vec3f closest_on_box;
	float sqr = SqrDistBoxPoint(half_box, trans_box, rot_box, trans_sphere, &closest_on_box);
	// outside box
	if (sqr > 0) {
		Vec3f diff		 = trans_sphere - closest_on_box,
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

bool PenetrationOBBOBB(PenetrationOutput& output, const Vec3f& half_A, const Vec3f& trans_A, const SO3f& rot_A, const Vec3f& half_B, const Vec3f& trans_B, const SO3f& rot_B)
{
	return false;
}


bool PenetrationOBBCapsule(PenetrationOutput& output, const Vec3f& half_box, const Vec3f& trans_box, const SO3f& rot_box, float radius, float half_height, const Vec3f& trans_cap, const SO3f& rot_cap)
{
	return false;
}

bool PenetrationCapsuleCapsule(PenetrationOutput& output, float radius_a, float half_height_a, const Vec3f& trans_a, const SO3f& rot_a, float radius_b, float half_height_b, const Vec3f& trans_b, const SO3f& rot_b)
{
	return false;
}

}