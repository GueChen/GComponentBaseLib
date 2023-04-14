#include "GComponent/Geometry/gcollision_detection.h"
#include "GComponent/Geometry/barycentric_coordinates.h"
#include "GComponent/GNumerical.hpp"

namespace GComponent{

static Vec3f NearestPointTriangle(std::vector<Vec3f>& simiplex, const uint32_t& out_plane_4, uint32_t* indices, uint32_t& size)
{
	float	 sqrt_dist		 = std::numeric_limits<float>::max();
	Vec3f    closest		 = Vec3f::Zero();
	uint32_t calc_indices[3] = { 0, 1, 2 };

	auto compare_nearest = [&](uint32_t i1, uint32_t i2, uint32_t i3) {
		calc_indices[0] = i1;
		calc_indices[1] = i2;
		calc_indices[2] = i3;

		uint32_t	calc_size = 3;
		Vec3f		calc_closest;

		const float calc_sqrt_dist = NearestTriangleBaryCentric(simiplex[i1], simiplex[i2], simiplex[i3], calc_indices, calc_size, calc_closest);

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
		sqrt_dist = NearestTriangleBaryCentric(simiplex[0], simiplex[1], simiplex[2], indices, size, closest);
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

static inline void PopBack(std::vector<Vec3f>& simplex, std::vector<Vec3f>& simplex_a, std::vector<Vec3f>& simplex_b) {
	simplex  .pop_back();
	simplex_a.pop_back();
	simplex_b.pop_back();
}

Vec3f NearestSimplex(std::vector<Vec3f>& simplex, std::vector<Vec3f>& simplex_a, std::vector<Vec3f>& simplex_b, Vec3f& support) {
	enum SimplexShape : uint32_t {
		ePoint = 1, eSegment, eTriangle, eTetrahedron
	} size = static_cast<SimplexShape>(simplex.size());

	switch (size) {
	case ePoint:	   return support;
	case eSegment:	   return NearestSegment	(simplex, simplex_a, simplex_b);
	case eTriangle:	   return NearestTriangle	(simplex, simplex_a, simplex_b);
	case eTetrahedron: return NearestTetrahedron(simplex, simplex_a, simplex_b);
	default:
		assert(false && "Nearest Simplex Size Error");
	}
	return support;
}

Vec3f NearestSegment(std::vector<Vec3f>& simplex, std::vector<Vec3f>& simplex_a, std::vector<Vec3f>& simplex_b) {
	
	Vec3f &a = simplex[0], &b = simplex[1];
	
	Vec3f ab = b - a,
		  ao =   - a;
	
	float sq_ab_length = ab.squaredNorm();						// 长度过短退化为点
	if (sq_ab_length < std::numeric_limits<float>::epsilon()) {	// degenerate to single point		
		PopBack(simplex, simplex_a, simplex_b);
		return simplex.front();
	}

	float scale = Clamp(ao.dot(ab) / sq_ab_length, 0.0f, 1.0f);

	return a + scale * ab;
}

Vec3f NearestTriangle(std::vector<Vec3f>& simplex, std::vector<Vec3f>& simplex_a, std::vector<Vec3f>& simplex_b) {
	Vec3f& a = simplex[0], & b = simplex[1], & c = simplex[2];
	Vec3f ab	= b - a,
		  ac	= c - a,
		  ab_ac = ab.cross(ac);

	if (float norm_length = ab_ac.squaredNorm();		// three point on one same line 
		norm_length < std::numeric_limits<float>::epsilon()) {															
		PopBack(simplex, simplex_a, simplex_b);			// degenerate to segment
		return NearestSegment(simplex, simplex_a, simplex_b);
	}

	uint32_t calc_size;
	uint32_t indices[3] = { 0, 1, 2 };
	Vec3f	 closest;

	NearestTriangleBaryCentric(a, b, c, indices, calc_size, closest);
	
	if (calc_size != 3) {
		std::vector<Vec3f> new_simplex	(calc_size),
						   new_simplex_a(calc_size),
						   new_simplex_b(calc_size);
		for (int i = 0; i < calc_size; ++i) {
			new_simplex[i]	 = simplex[indices[i]];
			new_simplex_a[i] = simplex_a[indices[i]];
			new_simplex_b[i] = simplex_b[indices[i]];
		}
		simplex  .swap(new_simplex);
		simplex_a.swap(new_simplex_a);
		simplex_b.swap(new_simplex_b);
	}

	return closest;
}

Vec3f NearestTetrahedron(std::vector<Vec3f>& simplex, std::vector<Vec3f>& simplex_a, std::vector<Vec3f>& simplex_b) {
	const Vec3f a = simplex[0],
				b = simplex[1],
				c = simplex[2],
				d = simplex[3];

	const Vec3f ab = b - a,
				ac = c - a,
				ad = d - a,
				n  = ab.cross(ac).normalized();
	if (abs(n.dot(ad)) < std::numeric_limits<float>::epsilon()) {	// degenerated to triangle case
		PopBack(simplex, simplex_a, simplex_b);
		return NearestTriangle(simplex, simplex_a, simplex_b);
	}

	uint32_t out_plane = PointOutsideOfPlane4(a, b, c, d);
	if (!out_plane) {											    // checking whether all false
		return Vec3f::Zero();
	}

	uint32_t	indices[3] = { 0, 1, 2 };
	uint32_t	calc_size;
	const Vec3f closest = NearestPointTriangle(simplex, out_plane, indices, calc_size);

	std::vector<Vec3f> new_simplex	(calc_size),
					   new_simplex_a(calc_size),
					   new_simplex_b(calc_size);
	for (int i = 0; i < calc_size; ++i) {
		new_simplex[i]	 = simplex[indices[i]];
		new_simplex_a[i] = simplex_a[indices[i]];
		new_simplex_b[i] = simplex_b[indices[i]];
	}
	simplex  .swap(new_simplex);
	simplex_a.swap(new_simplex_a);
	simplex_b.swap(new_simplex_b);
	
	return closest;
}


void GComponent::GetClosestPoint(Vec3f& closest_a, Vec3f&  closest_b, 
								 const std::vector<Vec3f>& simplex, 
								 const std::vector<Vec3f>& simplex_a, 
								 const std::vector<Vec3f>& simplex_b, 
								 const Vec3f& closest){
	enum SimplexShape : uint32_t {
		ePoint = 1, eSegment, eTriangle
	} size = static_cast<SimplexShape>(simplex.size());
	switch (size) {
	case ePoint: {
		closest_a = simplex_a.front();
		closest_b = simplex_b.front();
		return;
	}
	case eSegment: {
		float alpha;
		BarycentricCoordSegment(alpha, closest, simplex[0], simplex[1]);
		closest_a = (1.0 - alpha) * simplex_a[0] + alpha * simplex_a[1];
		closest_b = (1.0 - alpha) * simplex_b[0] + alpha * simplex_b[1];
		return;
	}
	case eTriangle: {
		float alpha, beta, gamma;
		BarycentricCoordTriangle(beta, gamma, closest, simplex[0], simplex[1], simplex[2]);
		alpha = 1.0 - beta - gamma;
		closest_a = alpha * simplex_a[0] + beta * simplex_a[1] + gamma * simplex_a[2];
		closest_b = alpha * simplex_b[0] + beta * simplex_b[1] + gamma * simplex_b[2];
		return;
	}
	default:
		assert(false && "This case should never happen");
	}
}
}