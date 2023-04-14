/**
 *  @file  	barycentric_coordinates.h
 *  @brief 	this file define two methods to caculate barycentric coordinates on segment/triangle.
 *  @author Gue Chen<guechen@buaa.edu.cn>
 *  @date 	Apr 8th, 2023
 **/
#ifndef __GBARYCENTRIC_COORDINATES_H
#define __GBARYCENTRIC_COORDINATES_H

#include <GComponent/types.h>
#include <Eigen/Dense>

namespace GComponent{
void BarycentricCoordSegment(float& alpha,							
							 const Vec3f& p, 
							 const Vec3f& s0, const Vec3f& s1);

void BarycentricCoordTriangle(float& beta, float& gamma, 
							  const Vec3f& p,
		  					  const Vec3f& t0, const Vec3f& t1, const Vec3f& t2);

}

#endif // !__GBARYCENTRIC_COORDINATES_H