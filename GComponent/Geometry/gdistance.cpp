
#include "gdistance.h"

#include <GComponent/gtransform.hpp>

namespace GComponent {

static float SqrDistFaceToPoint(int i0, int i1, int i2, const Vec3f& half_box, const Vec3f& pos_box_cor, Vec3f& pos, const Vec3f& dir, float* seg_param) {
	Vec3f p_to_neg_half = pos + half_box;
	float sqrt_dist = 0.0f;
	if (dir(i0) * p_to_neg_half(i1) >= dir(i1) * pos_box_cor(i0)) {
		if (dir(i0) * p_to_neg_half(i2) >= dir(i2) * pos_box_cor(i0)) {
			// 
			if (seg_param) {
				pos(i0) = half_box(i0);
				*seg_param = -pos_box_cor(i0) / dir(i0);
				pos(i1) += *seg_param * dir(i1);
				pos(i2) += *seg_param * dir(i2);
			}
		}
		else {
			//
			float line_sqr = dir(i0) * dir(i0) + dir(i2) * dir(i2);
			float temp = line_sqr * p_to_neg_half(i1) - dir(i1) * (dir(i0) * pos_box_cor(i0) + dir(i2) * p_to_neg_half(i2));
			if (temp <= 2.0f * line_sqr * half_box(i1)) {
				float t = temp / line_sqr;
				line_sqr += dir(i1) * dir(i1);
				temp = p_to_neg_half(i1) - t;
				float delta = dir(i0) * pos_box_cor(i0) + dir(i1) * temp + dir(i2) * p_to_neg_half(i2);
				float param = -delta / line_sqr;
				sqrt_dist += pos_box_cor(i0) * pos_box_cor(i0) + temp * temp + p_to_neg_half(i2) * p_to_neg_half(i2) + delta * param;
				if (seg_param) {
					*seg_param = param;
					pos(i0) = half_box(i0);
					pos(i1) = t - half_box(i1);
					pos(i2) = -half_box(i2);
				}
			}
			else {
				line_sqr += dir(i1) * dir(i1);
				float delta = dir(i0) * pos_box_cor(i0) + dir(i1) * pos_box_cor(i1) + dir(i2) * p_to_neg_half(i2);
				float param = -delta / line_sqr;
				sqrt_dist += pos_box_cor(i0) * pos_box_cor(i0) + pos_box_cor(i1) * pos_box_cor(i1) + p_to_neg_half(i2) * p_to_neg_half(i2) + delta * param;
				if (seg_param) {
					*seg_param = param;
					pos(i0) = half_box(i0);
					pos(i1) = half_box(i1);
					pos(i2) = -half_box(i2);
				}
			}
		}
	}
	else{
		if (dir(i0) * p_to_neg_half(i2) >= dir(i2) * pos_box_cor(i0)) {
			float line_sqr = dir(i0) * dir(i0) + dir(i1) * dir(i1);
			float temp = line_sqr * p_to_neg_half(i2) - dir(i2) * (dir(i0) * pos_box_cor(i0) + dir(i1) * p_to_neg_half(i1));
			if (temp <= 2.0f * line_sqr * half_box(i2)) {
				float t = temp / line_sqr;
				line_sqr += dir(i2) * dir(i2);
				temp = p_to_neg_half(i2) - t;
				float delta = dir(i0) * pos_box_cor(i0) + dir(i1) * p_to_neg_half(i1) + dir(i2) * temp;
				float param = -delta / line_sqr;
				sqrt_dist += pos_box_cor(i0) * pos_box_cor(i0) + p_to_neg_half(i1) * p_to_neg_half(i1) + temp * temp + delta * param;
				if (seg_param) {
					*seg_param = param;
					pos(i0) = half_box(i0);
					pos(i1) = -half_box(i1);
					pos(i2) = t - half_box(i2);
				}
			}
			else {
				line_sqr += dir(i2) * dir(i2);
				float delta = dir(i0) * pos_box_cor(i0) + dir(i1) * p_to_neg_half(i1) + dir(i2) * pos_box_cor(i2);
				float param = -delta / line_sqr;
				sqrt_dist += pos_box_cor(i0) * pos_box_cor(i0) + p_to_neg_half(i1) * p_to_neg_half(i1) + pos_box_cor(i2) * pos_box_cor(i2) + delta * param;
				if (seg_param) {
					*seg_param = param;
					pos(i0) = half_box(i0);
					pos(i1) = -half_box(i1);
					pos(i2) = half_box(i2);
				}
			}
		}
		else {
			float line_sqr = dir(i0) * dir(i0) + dir(i2) * dir(i2);
			float temp = line_sqr * p_to_neg_half(i1) - dir(i1) * (dir(i0) * pos_box_cor(i0) + dir(i2) * p_to_neg_half(i2));
			if (temp >= 0.0f) {
				if (temp <= 2.0f * line_sqr * half_box(i1)) {
					float t = temp / line_sqr;
					line_sqr += dir(i1) * dir(i1);
					temp = p_to_neg_half(i1) - t;
					float delta = dir(i0) * pos_box_cor(i0) + dir(i1) * temp + dir(i2) * p_to_neg_half(i2);
					float param = -delta / line_sqr;
					sqrt_dist += pos_box_cor(i0) * pos_box_cor(i0) + temp * temp + p_to_neg_half(i2) * p_to_neg_half(i2) + delta * param;
					if (seg_param) {
						*seg_param = param;
						pos(i0) = half_box(i0);
						pos(i1) = t - half_box(i1);
						pos(i2) = -half_box(i2);
					}
				}
				else {
					line_sqr += dir(i1) * dir(i1);
					float delta = dir(i0) * pos_box_cor(i0) + dir(i1) * pos_box_cor(i1) + dir(i2) * p_to_neg_half(i2);
					float param = -delta / line_sqr;
					sqrt_dist += pos_box_cor(i0) * pos_box_cor(i0) + pos_box_cor(i1) * pos_box_cor(i1) + p_to_neg_half(i2) * p_to_neg_half(i2) + delta * param;
					if (seg_param) {
						*seg_param = param;
						pos(i0) = half_box(i0);
						pos(i1) = half_box(i1);
						pos(i2) = -half_box(i2);
					}
				}
				return sqrt_dist;
			}

			line_sqr = dir(i0) * dir(i0) + dir(i1) * dir(i1);
			temp = line_sqr * p_to_neg_half(i2) - dir(i2) * (dir(i0) * pos_box_cor(i0) + dir(i1) * p_to_neg_half(i1));
			if (temp >= 0.0f)
			{
				// v(i2)-edge is closest
				if (temp <= 2.0f * line_sqr * half_box(i2))
				{
					float t = temp / line_sqr;
					line_sqr += dir(i2) * dir(i2);
					temp = p_to_neg_half(i2) - t;
					float delta = dir(i0) * pos_box_cor(i0) + dir(i1) * p_to_neg_half(i1) + dir(i2) * temp;
					float param = -delta / line_sqr;
					sqrt_dist += pos_box_cor(i0) * pos_box_cor(i0) + p_to_neg_half(i1) * p_to_neg_half(i1) + temp * temp + delta * param;

					if (seg_param)
					{
						*seg_param = param;
						pos(i0) = half_box(i0);
						pos(i1) = -half_box(i1);
						pos(i2) = t - half_box(i2);
					}
				}
				else
				{
					line_sqr += dir(i2) * dir(i2);
					float delta = dir(i0) * pos_box_cor(i0) + dir(i1) * p_to_neg_half(i1) + dir(i2) * pos_box_cor(i2);
					float param = -delta / line_sqr;
					sqrt_dist += pos_box_cor(i0) * pos_box_cor(i0) + p_to_neg_half(i1) * p_to_neg_half(i1) + pos_box_cor(i2) * pos_box_cor(i2) + delta * param;

					if (seg_param)
					{
						*seg_param = param;
						pos(i0) = half_box(i0);
						pos(i1) = -half_box(i1);
						pos(i2) = half_box(i2);
					}
				}
				return sqrt_dist;
			}

			// (v(i1),v(i2))-corner is closest
			line_sqr += dir(i2) * dir(i2);
			float delta = dir(i0) * pos_box_cor(i0) + dir(i1) * p_to_neg_half(i1) + dir(i2) * p_to_neg_half(i2);
			float param = -delta / line_sqr;
			sqrt_dist += pos_box_cor(i0) * pos_box_cor(i0) + p_to_neg_half(i1) * p_to_neg_half(i1) + p_to_neg_half(i2) * p_to_neg_half(i2) + delta * param;

			if (seg_param)
			{
				*seg_param = param;
				pos(i0) = half_box(i0);
				pos(i1) = -half_box(i1);
				pos(i2) = -half_box(i2);
			}			
		}
	}
	return sqrt_dist;
}

static float SqrDistBoxThreeDirLine(const Vec3f& half_box, Vec3f& line_pos, Vec3f& line_dir, float* seg_param) {
	// 
	// Not make sense this part 
	// ----------------------------------------------
	Vec3f pos_box_cor = line_pos - half_box;
	float sqrt_dist;

	float prod_px_dy = pos_box_cor.x() * line_dir.y();
	float prod_py_dx = pos_box_cor.y() * line_dir.x();
	if (prod_px_dy >= prod_py_dx) {									// x <====> y		
		float prod_px_dz = pos_box_cor.x() * line_dir.z();
		float prod_pz_dx = pos_box_cor.z() * line_dir.x();
		if (prod_px_dz >= prod_pz_dx) {								// x <====> z
			// first intersects x = half_x
			sqrt_dist = SqrDistFaceToPoint(0, 1, 2, half_box, pos_box_cor, line_pos, line_dir, seg_param);
		}
		else {
			// first intersects z = half_z
			sqrt_dist = SqrDistFaceToPoint(2, 0, 1, half_box, pos_box_cor, line_pos, line_dir, seg_param);
		}
	}
	else {
		float prod_py_dz = pos_box_cor.y() * line_dir.z();
		float prod_pz_dy = pos_box_cor.z() * line_dir.y();
		if (prod_py_dz >= prod_pz_dy) {
			// intersects y = half_y
			sqrt_dist = SqrDistFaceToPoint(1, 2, 0, half_box, pos_box_cor, line_pos, line_dir, seg_param);
		}
		else {
			// intersects z = half_z
			sqrt_dist = SqrDistFaceToPoint(2, 0, 1, half_box, pos_box_cor, line_pos, line_dir, seg_param);
		}

	}
	return sqrt_dist;
	
}

static float SqrDistBoxTwoDirLine(int i0, int i1, int i2, const Vec3f& half_box, Vec3f& line_pos, Vec3f& line_dir, float* seg_param) {
	// 退化为一个平面的例子
	// 直线存在两个可能的发射方向，则使用除法对比两个方向谁率先抵达目标平面
	// dist_0 / dir_0 >? dist_1 / dir_1
	// trick : 使用乘法化除为乘
	// -----------------------------------------------------------------------
	float sqrt_dist = 0.0f;
	float delta_i0 = line_pos(i0) - half_box(i0);	// 注意方向与通常定义的正向向反
	float delta_i1 = line_pos(i1) - half_box(i1);

	float prod_i0 = delta_i0 * line_dir(i1);
	float prod_i1 = delta_i1 * line_dir(i0);

	if (prod_i0 >= prod_i1) {						// i0 方向直线率先抵达，大的负值意味着更短时间
		line_pos[i0] = half_box[i0];
		float delta_pos_to_neg_half_i1 = line_pos(i1) + half_box(i1); // -(-half_box(i1))
		float delta = prod_i0 - line_dir(i0) * delta_pos_to_neg_half_i1;
		if (delta >= 0.0f) {
			float sqr = line_dir(i0) * line_dir(i0) + line_dir(i1) * line_dir(i1);
			sqrt_dist += delta * delta / sqr;
			if (seg_param) {
				line_pos(i1) = -half_box(i1);
				*seg_param   = -(line_dir(i0) * delta_i0 + line_dir(i1) * delta_pos_to_neg_half_i1) / sqr;
			}
		}
		else {
			if (seg_param) {
				line_pos(i1) -= prod_i0  / line_dir(i0);
				*seg_param   = -delta_i0 / line_dir(i0);
			}
		}
	}
	else {
		line_pos(i1) = half_box(i1);
		float delta_pos_to_neg_half_i0 = line_pos(i0) + half_box(i0); 
		float delta = prod_i1 - line_dir(i1) * delta_pos_to_neg_half_i0;
		if (delta >= 0.0f) {
			float sqr = line_dir(i0) * line_dir(i0) + line_dir(i1) * line_dir(i1);
			sqrt_dist += delta * delta / sqr;
			if (seg_param) {
				line_pos(i0) = -half_box(i0);
				*seg_param = -(line_dir(i0) * delta_pos_to_neg_half_i0 + line_dir(i1) * delta_i1) / sqr;
			}
		}
		else {
			if (seg_param) {
				line_pos(i0) -= prod_i1 / line_dir(i1); 
				*seg_param = -delta_i1 / line_dir(i1);
			}
		}
	}
	
	if (line_pos(i2) < -half_box(i2)) {
		float delta = -half_box(i2) - line_pos(i2);
		sqrt_dist += delta * delta;
		line_pos(i2) = -half_box(i2);
	}
	else if (line_pos(i2) > half_box(i2)) {
		float delta = -half_box(i2) + line_pos(i2);
		sqrt_dist += delta * delta;
		line_pos(i2) = half_box(i2);
	}

	return sqrt_dist;
}

static float SqrDistBoxOneDirLine(int i0, int i1, int i2, const Vec3f& half_box, Vec3f& line_pos, Vec3f& line_dir, float* seg_param) {
	// 考虑直线，仅在一个自由度可控时,距离由其余两个为约束方向控制
	// +----+   
	// |    |			line_param < 0.0f	
	// +----+====		返回后可确定原点最近
	//      | d				
	//      *===*p-->   非负方向同理，总能在直线上找到一点穿过平面，线段上最近点由返回参数判别
	
	if (seg_param) {		// 计算线段参数
		*seg_param = (half_box(i0) - line_pos(i0)) / line_dir(i0);
	}

	float sqrt_dist = 0.0f;	
	if (line_pos(i1) < -half_box(i1)) {
		float delta = -half_box(i1) - line_pos(i1);
		sqrt_dist += delta * delta;
		line_pos(i1) = -half_box(i1);
	}
	else if(line_pos(i1) > half_box(i1)) {
		float delta = -half_box(i1) + line_pos(i1);
		sqrt_dist += delta * delta;
		line_pos(i1) = half_box(i1);
	}
	if (line_pos(i2) < -half_box(i2)) {
		float delta = -half_box(i2) - line_pos(i2);
		sqrt_dist += delta * delta;
		line_pos(i2) = -half_box(i2);
	}
	else if (line_pos(i2) > half_box(i2)) {
		float delta = -half_box(i2) + line_pos(i2);
		sqrt_dist += delta * delta;
		line_pos(i2) = half_box(i2);
	}
	return sqrt_dist;
}

// https://github.com/NVIDIA-Omniverse/PhysX/blob/release/104.0/physx/source/geomutils/src/distance/GuDistanceSegmentBox.cpp 
// 
// 计算直线与盒子的最近点，通过将方向统一后考虑直线方向约束的条件来作优化
// 
// ---------------------------------------------------------------------------------------------------------------------------------
static float SqrDistBoxLine(const Vec3f& half_box, const Vec3f& trans_box, const Vec3f& rot_box,
							const Vec3f& seg_ori, const Vec3f& seg_dir,
							Vec3f* closet_on_box, float* seg_param) {
	//assert(false && "No Implementation");
	const SO3f  rot_mat_inv = Roderigues(rot_box).transpose();
	Vec3f ori_local = rot_mat_inv * (seg_ori - trans_box);
	Vec3f dir_local = rot_mat_inv * seg_dir;

	bool ref[3] = {false};
	for (unsigned int i = 0; i < 3; ++i) {
		if (dir_local(i) < 0.0f) {
			ori_local(i) = -ori_local(i);
			dir_local(i) = -dir_local(i);
			ref[i] = true;
		}		
	}
	
	float sqrt_dist = 0.0f;

	if (dir_local(0) == 0.0f) {
		if (dir_local(1) == 0.0f) {
			if (dir_local(2) == 0.0f) {
				if (seg_param) *seg_param = 0.0f;				
				sqrt_dist = SqrDistBoxPoint(half_box, trans_box, rot_box, seg_ori, closet_on_box);		// 直线退化为一个点
			}
			else {
				sqrt_dist = SqrDistBoxOneDirLine(2, 0, 1, half_box, ori_local, dir_local, seg_param);	// 直线被约束在仅存 x 方向
			}
		}
		else {
			if (dir_local(2) == 0.0f) {
				sqrt_dist = SqrDistBoxOneDirLine(1, 0, 2, half_box, ori_local, dir_local, seg_param);	// 直线被约束在仅存 y 方向
			}
			else {
				sqrt_dist = SqrDistBoxTwoDirLine(1, 2, 0, half_box, ori_local, dir_local, seg_param);	// 直线被约束在仅存 y z 方向
			}
		}
	}
	else {
		if (dir_local(1) == 0.0f) {
			if (dir_local(2) == 0.0f) {
				sqrt_dist = SqrDistBoxOneDirLine(0, 1, 2, half_box, ori_local, dir_local, seg_param);	// 直线被约束在仅存 z 方向
			}
			else {
				sqrt_dist = SqrDistBoxTwoDirLine(0, 2, 1, half_box, ori_local, dir_local, seg_param);	// 直线被约束在仅存 x z 方向
			}
		}
		else {
			if (dir_local(2) == 0.0f) {
				sqrt_dist = SqrDistBoxTwoDirLine(0, 1, 2, half_box, ori_local, dir_local, seg_param);	// 直线被约束在仅存 x y 方向
			}
			else {
				sqrt_dist = SqrDistBoxThreeDirLine(half_box, ori_local, dir_local, seg_param);			// 直线存在 x y z 方向
			}
		}
	}

	if (closet_on_box) {
		for (unsigned int i = 0; i < 3; ++i) {
			if (ref[i]) {
				ori_local(i) = -ori_local(i);
			}
		}
		*closet_on_box = ori_local;
	}

	return sqrt_dist;
}

float SqrDistBoxSeg(const Vec3f& half_box, const Vec3f& trans_box, const Vec3f& rot_box, 
					const Vec3f& seg_p_0, const Vec3f& seg_p_1, 
					Vec3f* closet_on_box, Vec3f* closet_on_seg)
{
	Vec3f closest_box;
	float seg_params;
	// 计算直线条件下的最近点

	float sqr_dist_line = SqrDistBoxLine(half_box, trans_box, rot_box, seg_p_0, seg_p_1 - seg_p_0, &closest_box, &seg_params);
	if (seg_params > 1.0f) {
		// negative dir case, p0 is closest point
		// 最近点在直线正方向且插值权重大于 1.0 ，因此线段最近点在p1
		//			     *p0
		// +--------+	/
		// |---*o---|  *p1
		// +--------+
		if (closet_on_seg) {
			*closet_on_seg = seg_p_1;
		}
		return SqrDistBoxPoint(half_box, trans_box, rot_box, seg_p_1, closet_on_box);
	}
	else if (seg_params >= 0.0f) {
		// closest point on segment
		// 最近点在两点之间，可直接返回结果
		//			     *p0
		// +--------+	 |
		// |---*o---|    *p1
		// +--------+
		if (closet_on_seg) {
			*closet_on_seg = (1.0 - seg_params) * seg_p_0 + seg_params * seg_p_1;
		}
		if (closet_on_box) {
			*closet_on_box = closest_box;
		}
		return sqr_dist_line;
	}
	else {
		// negative dir case, p0 is closest point
		// 最近点在直线负方向，因此线段最近点在p0
		//			     *p1
		// +--------+	/
		// |---*o---|  *p0
		// +--------+
		if (closet_on_seg) {
			*closet_on_seg = seg_p_0;
		}
		return SqrDistBoxPoint(half_box, trans_box, rot_box, seg_p_0, closet_on_box);
	}

	return 0.0f;
}

float SqrDistBoxPoint(const Vec3f& half_box, const Vec3f& trans_box, const Vec3f& rot_box,	// box params
					  const Vec3f& point,													// point params
					  Vec3f* closest_p_on_box)
{
	const Vec3f t_w     = point - trans_box;
	const SO3f  rot_mat = Roderigues(rot_box);
	const Vec3f t_box   = rot_mat.transpose() * t_w;
	
	Vec3f closest = Vec3f::Zero();
	float sqr_dist = 0.0f;

	for (int i = 0; i < 3; ++i) {
		//  check which region is point belonged
		//  判别检测点属于哪个区域
		//  -*p-|---*o--|---  
		//  
		//  *p 检测点   *o 盒子中心  | 盒子边界		
		// --------------------------------------------------------
		// 检测点为负值且小于盒子边界
		if (t_box(i) < -half_box(i)) {
			float delta = half_box(i) + t_box(i);			// half_i > 0 , t_box_i < 0
			sqr_dist  += delta * delta;
			closest(i) = -half_box(i);
		}
		else if (t_box(i) > half_box(i)) {
			float delta = t_box(i) - half_box(i);			// half_i > 0 , t_box_i > 0
			sqr_dist += delta * delta;
			closest(i) = half_box(i);
		}
	}
	if (closest_p_on_box) {
		*closest_p_on_box = rot_mat * closest;
	}
	return sqr_dist;
}

}


