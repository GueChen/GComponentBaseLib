#include "GComponent/Geometry/barycentric_coordinates.h"

#include <gtest/gtest.h>

#include <iostream>
#include <format>

using namespace std;
using namespace GComponent;

TEST(BarycentrcCoordTest, TestTriangle) {
	float beta, gamma;
	Vec3f a = Vec3f::Random(),
		  b = Vec3f::Random(),
		  c = Vec3f::Random();
	Vec3f p = a;
	BarycentricCoordTriangle(beta, gamma, p, a, b, c);
	EXPECT_EQ(beta, 0);
	EXPECT_EQ(gamma, 0);
	EXPECT_EQ(p, (1 - beta - gamma) * a + beta * b + gamma * c);
}

int main(int argc, char* argv[]) {
	cout << format("Test From {:}\n", __FILE__);
	testing::InitGoogleTest(&argc, argv);
	return RUN_ALL_TESTS();
}