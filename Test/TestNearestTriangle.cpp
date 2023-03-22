#include "GComponent/Geometry/gcollision_detection.h"
#include "GComponent/gtransform.hpp"

#include <iostream>
#include <format>
#include <random>

#include <gtest/gtest.h>

using namespace std;
using namespace GComponent;

constexpr const float kEPS = 1e-4f;

TEST(TestNearestTriangle, TriangleOrthogonalTest) {
	Vec3f    tri[3], A[3], B[3];	// used just for input
	uint32_t size	= 3;
	Vec3f    closest;
	auto do_orthorgnal_test = [&]() {
		closest = NearestTriangle(&tri[0], &A[0], &B[0], size);
		for (int i = 0; i < 3; ++i) {
			Vec3f edge = tri[(i + 1) % 3] - tri[i];
			EXPECT_LE(edge.dot(closest), kEPS) << format("edge {:} is not orthogonal with closest\n", i + 1) <<
				setiosflags(ios::fixed) << setprecision(4)
				<< "closest: " << closest.transpose() << "\n"
				<< "edge   : " << edge.transpose() << "\n";
		}
	};

	auto make_swap_test = [&](int idx1, int idx2) {
		do_orthorgnal_test();
		for (int i = 0; i < 3; ++i) {
			tri[i](idx1) *= -1.0f;
			tri[i](idx2) *= -1.0f;
		}
		do_orthorgnal_test();
		for (int i = 0; i < 3; ++i) {
			swap(tri[i](idx1), tri[i](idx2));
		}
		do_orthorgnal_test();
	};

	// x direction Test
	tri[0] = Vec3f(1.0f, 0.0f,  1.0f);
	tri[1] = Vec3f(1.0f, -5.0f, -1.0f);
	tri[2] = Vec3f(1.0f, 5.0f,  -1.0f);	
	make_swap_test(1, 2);

	for (int i = 0; i < 3; ++i) {
		tri[i] -= Vec3f::UnitX() * 2.0f;
	}	
	make_swap_test(1, 2);

	// test Y direction
	tri[0] = Vec3f(0.0f,  1.0f, 1.0f);
	tri[1] = Vec3f(-5.0f, 1.0f, -1.0f);
	tri[2] = Vec3f(5.0f,  1.0f, -1.0f);
	make_swap_test(0, 2);

	for (int i = 0; i < 3; ++i) {
		tri[i] -= Vec3f::UnitY() * 2.0f;
	}
	make_swap_test(0, 2);
	
	// test Z direction
	tri[0] = Vec3f(0.0f,  1.0f, 1.0f);
	tri[1] = Vec3f(-5.0f,-1.0f, 1.0f);
	tri[2] = Vec3f(5.0f, -1.0f, 1.0f);
	make_swap_test(1, 0);
	for (int i = 0; i < 3; ++i) {
		tri[i] -= Vec3f::UnitZ() * 2.0f;
	}
	make_swap_test(1, 0);
}

TEST(TestNearestTriangle, TestEdgeCase) {
	Vec3f    tri[3], A[3], B[3];	// used just for input
	uint32_t size = 3;
	Vec3f    closest;
	float    dot_result = 0.0f;

	random_device					 r;
	default_random_engine			 el(r());
	uniform_real_distribution<float> dis(-3.14f, 3.14f); //[-pi, pi]

	auto test_edge = [&](int idx1, int idx2, std::string edge_name) {
		assert(idx1 < 3 && idx2 < 3);
		// In NearestTriangle the ori triangle may be changed,
		// therefore using a copy form ori tri not directly make compare
		Vec3f tri_copy[3]; std::ranges::transform(tri, tri_copy, [](auto&& val) {return val; });

		closest = NearestTriangle(&tri_copy[0], &A[0], &B[0], size);
		dot_result = closest.dot((tri[idx2] - tri[idx1]).normalized());
		EXPECT_EQ(size, 2) << format("should degenerate to {:} edge\n", edge_name);
		EXPECT_LE(dot_result, kEPS) << 
				format("should be orthorgnal with {:} edge\n"
					   "real dot val = {:6.4f} > {:6.4f}\n", 
						edge_name, dot_result, kEPS);
	};
	auto add_random_rotation = [&]() {
		Vec3f random_rotation(dis(el), dis(el), dis(el));
		SO3f  R = Roderigues(random_rotation);
		for (auto& vert : tri) vert = R * vert;
	};


	tri[0] = Vec3f(1.0f, 0.0f, 15.0f);
	tri[1] = Vec3f(1.0f, -5.0f, 9.0f);
	tri[2] = Vec3f(1.0f, 5.0f,  9.0f);
	//      a
	//		/\
	//     /  \
	//    /__*_\ closest
	//	b	 |	c
	// 		 |
    //	     * ori
	//	
	test_edge(2, 1, "bc");	
	
	// add a random rotation should not change the closest point on bc edge
	add_random_rotation();
	test_edge(2, 1, "bc");

	// clockwise swaping 
	swap(tri[0], tri[1]);
	swap(tri[1], tri[2]);
	test_edge(1, 0, "ab");
	add_random_rotation();
	test_edge(1, 0, "ab");

	// clockwise swaping 
	swap(tri[0], tri[1]);
	swap(tri[1], tri[2]);
	test_edge(0, 2, "ac");
	add_random_rotation();
	test_edge(0, 2, "ac");
}

TEST(TestNearestTriangle, TestPointCase) {
	Vec3f    tri[3], A[3], B[3];	// used just for input
	uint32_t size = 3;
	Vec3f    closest;

	auto test_close_to = [&](int idx, std::string p_name) {
		// In NearestTriangle the ori triangle may be changed,
		// therefore using a copy form ori tri not directly make compare
		Vec3f tri_copy[3]; std::ranges::transform(tri, tri_copy, [](auto&& val) {return val; });

		closest = NearestTriangle(&tri_copy[0], &A[0], &B[0], size);
		EXPECT_EQ(size, 1) << "should degenerating to a point\n";
		EXPECT_LE((closest - tri[idx]).squaredNorm(), kEPS) <<
			format("should close on {:} point\n", p_name) << setiosflags(ios::fixed) << setprecision(4)
			<< format("{:<7}: ", p_name) << tri[0].transpose() << "\n"
			<< "closest: " << closest.transpose() << "\n";
	};

	// case A
	tri[0] = Vec3f(1.0f, 0.0f, -3.0f);
	tri[1] = Vec3f(1.0f, -5.0f, -8.0f);
	tri[2] = Vec3f(1.0f, 5.0f, -8.0f);
	test_close_to(0, "A");

	tri[0] = Vec3f(1.0f, 0.0f,  3.0f);
	tri[1] = Vec3f(1.0f, -5.0f, 8.0f);
	tri[2] = Vec3f(1.0f, 5.0f,  8.0f);
	test_close_to(0, "A");

	// case B	
	swap(tri[0], tri[1]);	
	test_close_to(1, "B");

	// case C	
	swap(tri[1], tri[2]);
	test_close_to(2, "C");
}

int main(int argc, char* argv[]) {
	cout << format("Test From {:}\n", __FILE__);
	testing::InitGoogleTest(&argc, argv);
	return RUN_ALL_TESTS();
}