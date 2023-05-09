#include <iostream>
#include <fstream>
#include <filesystem>
#include <format>

#include "GComponent/GGeometry.hpp"

using namespace std;
using namespace Eigen;
using namespace GComponent;

int main() {
	vector<Vector3f> poses(11);

	ofstream file("scatter.txt");
	cout << "inserted points:\n";
	for (int x = 0;  auto & pos : poses) {
		pos = Vector3f::Random() * 2.0f;
		pos.x() = x++;
		file <<std::format("{:} {:} {:}\n", pos.x(), pos.y(), pos.z());
	}
	
	auto run_and_write = [](auto&& func, std::string name){
		filesystem::path path = name;
		ofstream file(path);	
		if (!file.is_open()) {
			return;
		}
		for (float t = 0.0f; t <= 1.0f; t += 0.01f) {
			auto intepolated = func(t);
			file << format("{:} {:} {:}\n", intepolated.x(), intepolated.y(), intepolated.z());
		}
		file.close();
	};


	cout << "cubic intepolated:\n";
	auto cubic_spline = GetCubicSplineFunction(poses);
	run_and_write(cubic_spline, "cubic_curves.txt");


	cout << "bezier intepolated:\n";
	auto bezier_spline = GetBezierInterSplineFunction(poses);
	run_and_write(bezier_spline, "bezier_curves.txt");


	cout << "bspline intepolated:\n";
	BSpline bspline(poses, 2, true, BSplineNodeDefinition::QuasiUniform);
	run_and_write(bspline, "bspline_curves.txt");

}
