#include "GComponent/gtransform.hpp"

#include <gtest/gtest.h>

#include <iostream>
#include <format>

using namespace std;
using namespace Eigen;
using namespace GComponent;

constexpr const float kEPS = 1e-4f;

TEST(TestEulerAngle, EulerToso3Test) {
	Vector3f rot   = Vector3f::Zero();
	Vector3f euler = ToZYXEuler(rot);
	Vector3f ret   = FromZYXEuler(euler);

	EXPECT_EQ(rot, ret);

	rot   = Vector3f::Random();
	//std::cout << "angle :" << rot.norm() << std::endl;
	//std::cout << "axis  :" << (rot / rot.norm()).transpose() << std::endl;
	euler = ToZYXEuler(rot);
	//std::cout << "euler :" << euler.transpose() << std::endl;
	ret   = FromZYXEuler(euler);		
	EXPECT_LE((rot - ret).norm(), kEPS);

	rot = Vector3f::Random();

	euler = ToXYZEuler(rot);

	ret = FromXYZEuler(euler);
	EXPECT_LE((rot - ret).norm(), kEPS) << 
		"euler:" << euler.transpose() << "\n" << 
		"rotat:" << rot.transpose() << "\n" << 
		"compu:" << ret.transpose();
}

int main(int argc, char* argv[]) {
	cout << format("Test From {:}\n", __FILE__);
	testing::InitGoogleTest(&argc, argv);
	return RUN_ALL_TESTS();
}