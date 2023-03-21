#ifndef __GDISTANCE_H
#define __GDISTANCE_H

#include <GComponent/types.h>

namespace GComponent{

float SqrDistPointSeg(const Vec3f& point,
					  const Vec3f& seg_p_0, const Vec3f& seg_p_1,
					  Vec3f* closet_on_seg = nullptr);

//__declspec(dllimport)
float SqrDistBoxSeg(const Vec3f& half_box, const Vec3f& trans_box, const Vec3f& rot_box,
					const Vec3f& seg_p_0,  const Vec3f& seg_p_1,
					Vec3f* closet_on_box = nullptr,
					Vec3f* closet_on_seg = nullptr);

//__declspec(dllimport)
float SqrDistBoxPoint(const Vec3f& half_box, const Vec3f& trans_box, const Vec3f& rot_box, 
					  const Vec3f& point, 
					  Vec3f* closest_p_on_box = nullptr);

//__declspec(dllimport)
float SqrDistSegSeg(const Vec3f& a_p0, const Vec3f& a_p1, 
					const Vec3f& b_p0, const Vec3f& b_p1,
					Vec3f* closest_on_a = nullptr,
					Vec3f* closest_on_b = nullptr);

}
#endif // __GDISTANCE_H