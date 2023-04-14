#include <GComponent/Geometry/gcollision_detection.h>
#include <GComponent/GNumerical.hpp>


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

/*____________________________________GJK Related Methods____________________________________*/
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

	// based on voronoi region distance algorithm
	float scale = Clamp(ao.dot(ab) / sq_ab_length, 0.0f, 1.0f);

	return a + scale * ab;
}

Vec3f GComponent::NearestTriangle(Vec3f* Q, Vec3f* A, Vec3f* B, uint32_t& size)
{
	Vec3f& a = Q[0], & b = Q[1], & c = Q[2];
	Vec3f ab	= b - a,
		  ac	= c - a,
		  ab_ac = ab.cross(ac);

	if (float norm_length = ab_ac.squaredNorm();		// three point on one same line 
		norm_length < std::numeric_limits<float>::epsilon()) {	
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
	if (abs(n.dot(ad)) < std::numeric_limits<float>::epsilon()) {	// degenerated to triangle case
		size = 3;
		return NearestTriangle(Q, A, B, size);
	}

	uint32_t out_plane = PointOutsideOfPlane4(a, b, c, d);
	if (!out_plane) {											    // checking whether all false		
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
	if (va >= 0 && vb >= 0 && vc >= 0) {
		//		  dist * dir
		closest = (n.dot(a) / n_sqrt) * n;
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
			closest = d1 / recip * ab + a; 
		}
		else {
			closest = a;
		}
		return closest.dot(closest);
	}
	else if (va <= 0 && d4 >= d3 && d6 <= d5) {	// the closest point in bc edge segment
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
		if (abs(recip) > std::numeric_limits<float>::epsilon()) {
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



