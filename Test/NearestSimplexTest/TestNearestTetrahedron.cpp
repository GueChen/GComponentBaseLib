#include "GComponent/Geometry/gcollision_detection.h"
#include "GComponent/gtransform.hpp"

#include <iostream>
#include <format>
#include <random>

#include <gtest/gtest.h>
using namespace std;
using namespace GComponent;

constexpr const float kEPS = 1e-4f;

TEST(TestNearestTetrahedron, Degen2TriangleTest) {
	Vec3f	 te[4], A[4], B[4];
	uint32_t size = 4;
	Vec3f	 closest;
	uint32_t index[3];

	random_device					 r;
	default_random_engine			 el(r());
	uniform_real_distribution<float> dis(-3.14f, 3.14f); //[-pi, pi]

#define IndexSet(idx1, idx2, i3)\
	do {					\
		index[0] = idx1;		\
		index[1] = idx2;		\
		index[2] = i3;		\
	} while(0)

	auto make_test = [&](int i1, int i2, int i3, std::string plane_name) {
		
		uint32_t size_copy = size;
		Vec3f	 te_copy[4]; std::ranges::transform(te, te_copy, [](auto&& val) {return val; });
		
		closest = NearestTetrahedron(te_copy, A, B, size_copy);

		EXPECT_GE(closest.squaredNorm(), 0.0f) << "closest should not close to zero\n";
		EXPECT_EQ(size_copy, 3) << "size should be 3, one face close\n";
		IndexSet(i1, i2, i3);
		for (int i = 0; i < 3; ++i) {
			uint32_t idx1 = index[(i + 1) % 3], idx2 = index[i];
			Vec3f edge = (te[idx1] - te[idx2]).normalized();
			EXPECT_LE(edge.dot(closest), kEPS) << 
				format("closest should be orthogonal with plane {:}\n"
					   "cur not orthorgnal with edge constructed by {:}-{:}\n",
						"abc", idx1, idx2)
				<< setiosflags(ios::fixed) << setprecision(4)
				<< "closest: " << closest.transpose() << "\n"
				<< "edge   : " << edge.transpose()	  << "\n";
		}
	};
	auto add_random_rotation = [&]() {
		Vec3f random_rotation(dis(el), dis(el), dis(el));
		SO3f  R = Roderigues(random_rotation);
		for (auto& vert : te) vert = R * vert;
	};

	te[0] = Vec3f(1.0f,  0.0f,  1.0f);
	te[1] = Vec3f(1.0f,  -5.0f, -1.0f);
	te[2] = Vec3f(1.0f,  5.0f,  -1.0f);
	te[3] = Vec3f(10.0f, 0.0f,  0.0f);
	make_test(0, 1, 2, "abc");

	add_random_rotation();
	swap(te[3], te[0]);
	make_test(1, 2, 3, "bcd");
	
	add_random_rotation();
	swap(te[0], te[1]);
	make_test(0, 2, 3, "acd");

	add_random_rotation();
	swap(te[1], te[2]);
	make_test(0, 1, 3, "abd");
}

TEST(TestNearestTetrahedron, Degen2SegmentTest) {
	Vec3f	 te[4], A[4], B[4];
	uint32_t size = 4;
	Vec3f	 closest;
	uint32_t index[3];
	float dot_result = 0.0f;

	te[0] = Vec3f(1.0f, 0.0f,  8.0f);
	te[1] = Vec3f(1.0f, -5.0f, 5.0f);
	te[2] = Vec3f(1.0f, 5.0f,  5.0f);
	te[3] = Vec3f(15.0f, 0.0f, 7.0f);

	closest = NearestTetrahedron(te, A, B, size);

	
	dot_result = closest.dot((te[2] - te[1]).normalized());
	EXPECT_EQ(size, 2) << "size should be 2 degenerating to a seg\n";
	EXPECT_LE(dot_result, kEPS) <<
		format("should be orthorgnal with {:} edge\n"
			"real dot val = {:6.4f} > {:6.4f}\n",
			"bc", dot_result, kEPS);
}

TEST(TestNearestTetrahedron, Degen2PointTest) {
	Vec3f	 te[4], A[4], B[4];
	uint32_t size = 4;
	Vec3f	 closest;

	float result = 0.0f;

	te[0] = Vec3f(5.0f, 0.0f,   5.0f);
	te[1] = Vec3f(5.0f, -5.0f, -5.0f);
	te[2] = Vec3f(5.0f, 5.0f,  -5.0f);
	te[3] = Vec3f(3.0f, 0.0f,  0.0f);

	closest = NearestTetrahedron(te, A, B, size);

	result = (closest - te[3]).squaredNorm();
	EXPECT_EQ(size, 1) << "size should be 2 degenerating to a seg\n";
	EXPECT_LE(result, kEPS) <<
		format("should be close to point {:}\n"
			"real dot val = {:6.4f} > {:6.4f}\n",
			"D", result, kEPS);
}

TEST(TestNearestTetrahedron, ContainZeroTest) {
	Vec3f	 te[4], A[4], B[4];
	uint32_t size = 4;
	Vec3f	 closest;

	float result = 0.0f;

	te[0] = Vec3f(5.0f, 0.0f, 5.0f);
	te[1] = Vec3f(5.0f, -5.0f, -5.0f);
	te[2] = Vec3f(5.0f, 5.0f, -5.0f);
	te[3] = Vec3f(-5.0f, 0.0f, 0.0f);

	closest = NearestTetrahedron(te, A, B, size);

	result = (closest).squaredNorm();	
	EXPECT_LE(result, kEPS) <<
		format("should be close to zero point\n"
			"real dot val = {:6.4f} > {:6.4f}\n",
			 result, kEPS);
}

int main(int argc, char*argv[]) {
	cout << format("{:}-Test From : {:}\n","Nearest Tetrahedron", __FILE__);
	testing::InitGoogleTest(&argc, argv);
	return RUN_ALL_TESTS();
}