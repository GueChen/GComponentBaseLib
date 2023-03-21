#include "GComponent/Geometry/gcollision_detection.h"
#include "GComponent/GNumerical.hpp"

#include <iostream>
#include <format>
#include <random>

#include <gtest/gtest.h>

using namespace std;
using namespace GComponent;

constexpr const float kEPS = 1e-4f;

TEST(TestNearestSegment, SegmentOrthogonalTest) {
	// the nearest point must be orthogonal with seg or on the sides
	Vec3f	 seg[2];
	uint32_t size = 2;

	seg[0] = Vec3f(1.0f, 1.0f,  1.0f);
	seg[1] = Vec3f(5.0f, -1.0f, -3.0f);
	
	auto closest = NearestSegment(&seg[0], size);
	ASSERT_FALSE((closest - seg[0]).squaredNorm() < kEPS) << "closest should not on point A\n";
	ASSERT_FALSE((closest - seg[1]).squaredNorm() < kEPS) << "closest should not on point B\n";
	EXPECT_TRUE(closest.dot(seg[1] - seg[0]) < kEPS)	  << "closest should be orthogonal with AB";
	
}

TEST(TestNearestSegment, SegmentZeroTest) {
	Vec3f	 seg[2];
	uint32_t size = 2;
	seg[0] = Vec3f(1.0f, 1.0f, 1.0f);
	seg[1] = Vec3f(-1.0f, -1.0f, -1.0f);
	auto closest = NearestSegment(&seg[0], size);
	EXPECT_TRUE(closest.squaredNorm() < kEPS)			  << "closest should be close to zero\n";
}

TEST(TestNearestSegment, SegmentOneSideTest) {
	Vec3f	 seg[2];
	uint32_t size = 2;
	seg[0] = Vec3f(1.0f, 2.0f, 3.0f);
	seg[1] = Vec3f(5.0f, 5.0f, 5.0f);
	
	auto closest = NearestSegment(&seg[0], size);
	EXPECT_TRUE((closest - seg[0]).squaredNorm() < kEPS)  << "closest should close to point A\n";

	seg[1] = Vec3f(-1.0f, 2.0f, -3.0f);
	seg[0] = Vec3f(-5.0f, 5.0f, -5.0f);

	closest = NearestSegment(&seg[0], size);
	EXPECT_TRUE((closest - seg[1]).squaredNorm() < kEPS) << "closest should close to point B\n";
}

TEST(TestNearestSegment, RandomTest) {
	
	Vec3f	 seg[2];
	uint32_t size = 2;
	
	random_device					 r;
	default_random_engine			 el(r());
	uniform_real_distribution<float> dis(-99.99f, 99.99f);

#define Gen dis(el)

	seg[0] = Vec3f(Gen, Gen, Gen);
	seg[1] = Vec3f(Gen, Gen, Gen);

	auto closest = NearestSegment(&seg[0], size);

	Vec3f ab = seg[0] - seg[1];
	float sq_ab = (seg[0] - seg[1]).squaredNorm();
	float dot = closest.dot(ab);
	float eps = sq_ab < 1e3f ? kEPS : (sq_ab < 1e4f ? 1e-3f : 1e-2f);
	if (bool dot_result = (dot < eps);
		dot_result) {
		EXPECT_TRUE(dot_result);
	}
	else {		
		EXPECT_TRUE((closest - seg[1]).squaredNorm() < eps || 
					(closest - seg[0]).squaredNorm() < eps) 
			<< std::format("{:<6.4f} val greater than {:<6.5f}, should on one side\n"
						   "l_ab^2 = {:8.4f}",
						   dot, eps, sq_ab);
	}

#undef Gen
}

int main(int argc, char* argv[]) {
	cout << format("Test From {:}\n", __FILE__);
	testing::InitGoogleTest(&argc, argv);
	return RUN_ALL_TESTS();
}