#include "GComponent/Geometry/barycentric_coordinates.h"

namespace GComponent{

void BarycentricCoordSegment(float& alpha, const Vec3f& p, const Vec3f& a, const Vec3f& b){
	const Vec3f pa = a - p,
				pb = b - p;
	const Vec3f d  = pa - pb;
	if (float denom = d.squaredNorm();
		denom >= 0) {
		alpha = -pa.dot(d) / denom;
	}
	else {
		alpha = 0;
	}

}

void BarycentricCoordTriangle(float& beta, float& gamma, 
							  const Vec3f& p, 
							  const Vec3f& a, const Vec3f& b, const Vec3f& c){
	const Vec3f ab = b - a, ac = c - a,
				pa = a - p, pb = b - p, pc = c - p;

	const Vec3f area_abc = ab.cross(ac),
				area_pbc = pb.cross(pc), area_pca = pc.cross(pa), area_pab = pa.cross(pb);
	
	const float va = area_abc.dot(area_pbc), vb = area_abc.dot(area_pca), vc = area_abc.dot(area_pab);

	if (float denom = va + vb + vc;
		denom != 0) {
		beta  = vb / denom;
		gamma = vc / denom;
	}
	else {	// special case, degenerating to seg/point
		beta = gamma = 0;
	}
}

}