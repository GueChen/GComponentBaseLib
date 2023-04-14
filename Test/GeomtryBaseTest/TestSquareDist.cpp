#include "GComponent/Geometry/gdistance.h"

#include <gtest/gtest.h>

#include <Eigen/Dense>

#include <iostream>
#include <format>



using namespace std;
using namespace GComponent;

TEST(SquareDistanceTest, TestBoxPoint) {
	Vec3f box_half  = 1.0f * Vec3f::Ones();
	Vec3f box_trans = Vec3f(2.0f, 0.0f, 0.0f);
	Vec3f box_rot   = Vec3f::Zero();
	Vec3f closest, p, expect;
	float dist;
	// inner case
	p = Vec3f(1.5f, 0.6f, -0.35f);
	dist = SqrDistBoxPoint(box_half, box_trans, box_rot, p, &closest);
	expect = p;
	EXPECT_EQ(closest, expect);

	// one outter case
	p = Vec3f::Zero();		
	dist = SqrDistBoxPoint(box_half, box_trans, box_rot, p, &closest);		
	expect = Vec3f::UnitX();
	EXPECT_EQ(closest, expect);

	// one outter with other paramter
	p = Vec3f(0.0f, 0.8f, -0.4f);	
	dist = SqrDistBoxPoint(box_half, box_trans, box_rot, p, &closest);
	expect = Vec3f(1.0f, 0.8f, -0.4f);
	EXPECT_EQ(closest, expect);

	// two outter
	p = Vec3f(0.0f, 1.2f, -0.55f);
	dist = SqrDistBoxPoint(box_half, box_trans, box_rot, p, &closest);
	expect = Vec3f(1.0f, 1.0f, -0.55f);
	EXPECT_EQ(closest, expect);

	// three outter
	p = Vec3f(-50.0f, -6.4f, 5.2f);
	dist = SqrDistBoxPoint(box_half, box_trans, box_rot, p, &closest);
	expect = Vec3f(1.0f, -1.0f, 1.0f);
	EXPECT_EQ(closest, expect);
}

int main(int argc, char* argv[]) {
	cout << format("Test From {:}\n", __FILE__);
	testing::InitGoogleTest(&argc, argv);
	return RUN_ALL_TESTS();
}