#ifndef __GCOLLISION_DETECTION_H
#define __GCOLLISION_DETECTION_H

#include <GComponent/types.h>
#include <Concept/gconcept.hpp>

#include <limits>

#ifdef __GJK_CHECKING_PRINTING
#include <iostream>
#include <format>
#endif

namespace GComponent { 



bool IntersectSphereSphere(float radius_A, const Vec3f& trans_A,
						   float radius_B, const Vec3f& trans_B);

bool IntersectSphereCapsule(float radius_sphere,  const Vec3f& trans_sphere,
							float radius_capsule, float half_height_capsule, const Vec3f& trans_cap, const Vec3f rot_cap);
bool IntersectSphereCapsule(float radius_sphere,  const Vec3f& trans_sphere,
							float radius_capsule, float half_height_capsule, const Vec3f& trans_cap, const SO3f& rot_cap);

/// <summary>
/// Detects whther a pair of bounding boxes intersects with each other
/// <para>
/// 检测一对包围盒是否相交
/// </para>
/// </summary>
/// <param name="half_A">	cref	{vec3}	[in]	  A 盒的 3 x 1 半长向量		</param>
/// <param name="trans_A">	cref	{vec3}	[in]	  A 盒的 3 x 1 偏移向量	    </param>
/// <param name="rot_A">	cref	{vec3}	[in]	  A 盒的 3 x 1 旋转向量	    </param>
/// <param name="half_B">	cref	{vec3}	[in]	  B 盒的 3 x 1 半长向量		</param>
/// <param name="trans_B">  cref	{vec3}	[in]	  B 盒的 3 x 1 偏移向量		</param>
/// <param name="rot_B">    cref	{vec3}	[in]	  B 盒的 3 x 1 旋转向量		</param>
/// <returns>				val		{bool}	[out]	  相交检测的结果			</returns>
bool IntersectOBBOBB	 (const Vec3f& half_A, const Vec3f& trans_A, const Vec3f& rot_A, 
						  const Vec3f& half_B, const Vec3f& trans_B, const Vec3f& rot_B);
bool IntersectOBBOBB	 (const Vec3f& half_A, const Vec3f& trans_A, const SO3f& rot_A, 
						  const Vec3f& half_B, const Vec3f& trans_B, const SO3f& rot_B);

/// <summary>
/// Detects if a pair of obb and sphere is intersected with each other
/// <para>
/// 检测一对包围盒和球体是否相交
/// </para>
/// </summary>
/// <param name="half_box">			cref	{vec3}	[in]	  包围盒的 3 x 1 半长向量	</param>
/// <param name="trans_box">		cref	{vec3}	[in]	  包围盒的 3 x 1 偏移向量	</param>
/// <param name="rot_box">			cref	{vec3}	[in]	  包围盒的 3 x 1 旋转向量	</param>
/// <param name="radius_sphere">	val		{float} [in]	  球体的半径长度			</param>
/// <param name="trans_sphere">		cref	{vec3}	[in]	  球体的 3 x 1 偏移向量		</param>
/// <returns>						val		{bool}	[out]	  相交检测的结果			</returns>
bool IntersectOBBSphere  (const Vec3f& half_box, const Vec3f& trans_box, const Vec3f& rot_box, 
						  float radius_sphere,   const Vec3f& trans_sphere);
bool IntersectOBBSphere  (const Vec3f& half_box, const Vec3f& trans_box, const SO3f& rot_box, 
						  float radius_sphere,   const Vec3f& trans_sphere);

/// <summary>
/// Detects if obb and capsule pair is intersected with each other
/// <para>
/// 检测包围盒和胶囊体对是否相交
/// </para>
/// </summary>
/// <param name="half_box">   cref	{vec3}	[in]	包围盒的 3 x 1 半长向量	</param>
/// <param name="trans_box">  cref	{vec3}	[in]	包围盒的 3 x 1 偏移向量	</param>
/// <param name="rot_box">    cref	{vec3}	[in]	包围盒的 3 x 1 旋转向量	</param>
/// <param name="radius">     val   {float} [in]	胶囊体的半径长度        </param>
/// <param name="half_height">val   {float} [in]	胶囊体的半长轴长度		</param>
/// <param name="trans_cap">  cref	{vec3}  [in]	胶囊体的 3 x 1 偏移向量 </param>
/// <param name="rot_cap">    cref  {vec3}  [in]	胶囊体的 3 x 1 旋转向量 </param>
/// <returns>				  val   {bool}  [out]   相交检测的结果		    </returns>
bool IntersectOBBCapsule(const Vec3f& half_box, const Vec3f& trans_box, const Vec3f& rot_box, 
						  float radius, float half_height, const Vec3f& trans_cap, const Vec3f& rot_cap);
bool IntersectOBBCapsule(const Vec3f& half_box, const Vec3f& trans_box, const SO3f& rot_box, 
						  float radius, float half_height, const Vec3f& trans_cap, const SO3f& rot_cap);

bool IntersectCapsuleCapsule(float radius_a,	   float half_height_a, 
							 const Vec3f& trans_a, const SO3f& rot_a,
							 float radius_b,	   float half_height_b,
							 const Vec3f& trans_b, const SO3f& rot_b);

/*__________________________________________ GJK Related Methods ________________________________________*/
template<IsConvex ConvexA, IsConvex ConvexB>
bool IntersectGJK		(const ConvexA& a, const ConvexB& b, const Vec3f& search_dir);

Vec3f NearestSimplex	(Vec3f* Q, Vec3f* A, Vec3f* B, uint32_t& size, Vec3f& support);

Vec3f NearestSegment	(Vec3f* Q, uint32_t& size);

Vec3f NearestTriangle	(Vec3f* Q, Vec3f* A, Vec3f* B, uint32_t& size);

Vec3f NearestTetrahedron(Vec3f* Q, Vec3f* A, Vec3f* B, uint32_t& size);

float NearestTriangleBaryCentric(Vec3f& a, Vec3f& b, Vec3f& c, uint32_t* indices, uint32_t& size, Vec3f& closest);

uint32_t PointOutsideOfPlane4(const Vec3f& a, const Vec3f& b, const Vec3f& c, const Vec3f& d);

bool IntersectEPA();

template<IsConvex ConvexA, IsConvex ConvexB>
bool IntersectGJK(const ConvexA& a, const ConvexB& b, const Vec3f& search_dir)
{
	Vec3f Q[4], A[4], B[4];
	uint32_t simplex_size = 0;

	float dist        = std::numeric_limits<float>::max();
	float prev_dist;

	Vec3f closest     = search_dir;
	Vec3f closest_dir = closest.normalized();

	float eps = 1e-4f;
#ifdef __GJK_CHECKING_PRINTING
	int   iter = 1;
#endif
	do {
		prev_dist = dist;

		Vec3f support_a = a.Support(-closest_dir);
		Vec3f support_b = b.Support(closest_dir);
		Vec3f support   = support_a - support_b;
		
		assert(simplex_size < 4, "Simplex dimension not greater than 4");

		A[simplex_size] = support_a;
		B[simplex_size] = support_b;
		Q[simplex_size] = support;
		++simplex_size;
		
		closest = NearestSimplex(Q, A, B, simplex_size, support);
		dist = closest.norm();
		closest_dir = closest / dist;

#ifdef __GJK_CHECKING_PRINTING
		std::cout << std::format(
			"iter {:<3}:\n", iter++
		);
		std::cout << "closest   : " << closest.transpose()   << std::endl;
		std::cout << "support_a : " << support_a.transpose() << std::endl;
		std::cout << "support_b : " << support_b.transpose() << std::endl;
		std::cout << "support   : " << support.transpose()   << std::endl;
#endif
	} while (prev_dist > dist && dist > eps);

	return dist <= eps;
}

}
#endif // !__GCOLLISION_DETECTION_H


