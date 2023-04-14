/**
 *  @file  	geom_utils.h
 *  @brief 	some common type as input and output definitions was setted here.
 *  @author Gue Chen<guechen@buaa.edu.cn>
 *  @date 	Apr 9th, 2023
 **/
#ifndef __G_GEOM_UTILS_H
#define __G_GEOM_UTILS_H

#include <GComponent/types.h>

namespace GComponent {
enum   ContactStatus {

};

struct PenetrationOutput {
	Vec3f closest_a = Vec3f::Zero();
	Vec3f closest_b = Vec3f::Zero();
	Vec3f normal	= Vec3f::Zero();	
	float depth		= 0.0f;
};

} // !namespace GComponent

#endif // !__G_GEOM_UTILS_H
