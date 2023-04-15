#ifndef __GDISTANCE_H
#define __GDISTANCE_H

#include <GComponent/types.h>

namespace GComponent{

/// <summary>
/// ����� {point} ���߶� {seg} ����С����ƽ������ {closest_on_seg} �� nullptr���������㵽�ò���
/// </summary>
/// <param name="point">�ռ��</param>
/// <param name="seg_p_0">�߶ζ˵� 0</param>
/// <param name="seg_p_1">�߶ζ˵� 1</param>
/// <param name="closet_on_seg">�߶������ָ�룬��Ϊ��</param>
/// <returns>�ռ�㵽�߶ε�ƽ������</returns>
float SqrDistPointSeg(const Vec3f& point,
					  const Vec3f& seg_p_0, const Vec3f& seg_p_1,
					  Vec3f* closet_on_seg = nullptr);

/// <summary>
/// �����߶� {seg} �� {box} ����Сƽ�����룬�� {closest_*} �ǿգ���������������ò�����
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
/// �����Χ�� {box} ���� {point} ����Сƽ�����룬�� {closest_*} �ǿգ������������ò���
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
/// �����߶� {seg_a} {seg_b} �����С����ƽ������ {closest_*} �ǿգ������Ӧ��������λ��
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