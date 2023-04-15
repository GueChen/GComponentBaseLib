#ifndef __GDISTANCE_H
#define __GDISTANCE_H

#include <GComponent/types.h>

namespace GComponent{

/// <summary>
/// 计算点 {point} 到线段 {seg} 的最小距离平方，若 {closest_on_seg} 非 nullptr，输出最近点到该参数
/// </summary>
/// <param name="point">空间点</param>
/// <param name="seg_p_0">线段端点 0</param>
/// <param name="seg_p_1">线段端点 1</param>
/// <param name="closet_on_seg">线段最近点指针，可为空</param>
/// <returns>空间点到线段的平方距离</returns>
float SqrDistPointSeg(const Vec3f& point,
					  const Vec3f& seg_p_0, const Vec3f& seg_p_1,
					  Vec3f* closet_on_seg = nullptr);

/// <summary>
/// 计算线段 {seg} 到 {box} 的最小平方距离，若 {closest_*} 非空，输出最近点参数至该参数中
/// </summary>
/// <param name="half_box"></param>
/// <param name="trans_box"></param>
/// <param name="rot_box"></param>
/// <param name="seg_p_0"></param>
/// <param name="seg_p_1"></param>
/// <param name="closet_on_box"></param>
/// <param name="closet_on_seg"></param>
/// <returns></returns>
float SqrDistBoxSeg(const Vec3f& half_box, const Vec3f& trans_box, const Vec3f& rot_box,
					const Vec3f& seg_p_0,  const Vec3f& seg_p_1,
					Vec3f* closet_on_box = nullptr,
					Vec3f* closet_on_seg = nullptr);
float SqrDistBoxSeg(const Vec3f& half_box, const Vec3f& trans_box, const SO3f& rot_box,
					const Vec3f& seg_p_0,  const Vec3f& seg_p_1,
					Vec3f* closet_on_box = nullptr,
					Vec3f* closet_on_seg = nullptr);
/// <summary>
/// 计算包围盒 {box} 到点 {point} 的最小平方距离，若 {closest_*} 非空，输出最近点至该参数
/// </summary>
/// <param name="half_box"></param>
/// <param name="trans_box"></param>
/// <param name="rot_box"></param>
/// <param name="point"></param>
/// <param name="closest_p_on_box"></param>
/// <returns></returns>
float SqrDistBoxPoint(const Vec3f& half_box, const Vec3f& trans_box, const Vec3f& rot_box, 
					  const Vec3f& point, 
					  Vec3f* closest_p_on_box = nullptr);
float SqrDistBoxPoint(const Vec3f& half_box, const Vec3f& trans_box, const SO3f& rot_mat,
					  const Vec3f& point,
					  Vec3f* closest_p_on_box = nullptr);

/// <summary>
/// 计算线段 {seg_a} {seg_b} 间的最小距离平方，若 {closest_*} 非空，输出对应参数至该位置
/// </summary>
/// <param name="a_p0"></param>
/// <param name="a_p1"></param>
/// <param name="b_p0"></param>
/// <param name="b_p1"></param>
/// <param name="closest_on_a"></param>
/// <param name="closest_on_b"></param>
/// <returns></returns>
float SqrDistSegSeg(const Vec3f& a_p0, const Vec3f& a_p1, 
					const Vec3f& b_p0, const Vec3f& b_p1,
					Vec3f* closest_on_a = nullptr,
					Vec3f* closest_on_b = nullptr);


}
#endif // __GDISTANCE_H