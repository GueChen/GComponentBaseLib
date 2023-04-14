#ifndef __GJK_CONVEX_H
#define __GJK_CONVEX_H

#include <GComponent/types.h>

namespace GComponent {

enum GJKStatus {
	NON_INTERSECT,	
	CONTACT,	

	EPA_CONTACT,
	EPA_DEGENERATE,
	EPA_FAIL
};

struct GJKOutput {	
	Vec3f closest_a	 = Vec3f::Zero();
	Vec3f closest_b	 = Vec3f::Zero();
	Vec3f normal	 = Vec3f::Zero();
	Vec3f search_dir = Vec3f::Zero();
	float depth		 = 0.0f;
};

/// <summary>
/// Base class for GJK Convexhull collision check
/// </summary>
class GJKConvex {
public:
	GJKConvex() = default;
	virtual Vec3f Support(const Vec3f& dir) const = 0;
};

}


#endif // !__GJK_CONVEX_H
