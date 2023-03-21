#include "GComponent/Geometry/gcollision_detection.h"

#include <iostream>
#include <format>
#include <random>

#include <gtest/gtest.h>

using namespace std;
using namespace GComponent;

constexpr const float kEPS = 1e-4f;

TEST(TestNearestTriangle, TriangleOrthogonalTest) {
	Vec3f    tri[3];
	Vec3f    A[3], B[3];	// used just for input
	uint32_t size = 3;

	tri[0] = Vec3f(1.0f, 1.0f,  1.0f);
	tri[1] = Vec3f(1.0f, -5.0f, -1.0f);
	tri[2] = Vec3f(1.0f, 5.0f,  1.0f);

	auto closest = NearestTriangle(&tri[0], &A[0], &B[0], size);

	for (int i = 0; i < 3; ++i) {
		Vec3f edge = tri[(i + 1) % 3] - tri[i];
		EXPECT_LE(edge.dot(closest), kEPS) << format("edge {:} is not orthogonal with closest\n", i + 1)<<
			setiosflags(ios::fixed) << setprecision(4) 
			<< "closest: " << closest.transpose() << "\n"
			<< "edge   : " << edge.transpose()    << "\n";
	}
}


int main(int argc, char* argv[]) {
	cout << format("Test From {:}\n", __FILE__);
	testing::InitGoogleTest(&argc, argv);
	return RUN_ALL_TESTS();
}