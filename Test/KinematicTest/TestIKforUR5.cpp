#include "GComponent/grobotkinematic.h"

#include "GComponent/GNumerical.hpp"
#include "GComponent/gtransform.hpp"
#include "LSSolver/LinearSystemSolver.hpp"

#include <iostream>
#include <chrono>

using namespace GComponent;
using namespace Eigen;
using namespace std;
using namespace RobotKinematic;
constexpr double kPi = 3.1415926535;

int main() {
	int N = 1000;
	vector<Vector<double, Dynamic>> thetas(N);

	for (int i = 0; i < thetas.size(); ++i) {
		auto vals = GenUniformRandom<double, 6>(-kPi, kPi);
		thetas[i] = Map<Vector<double, 6>>(vals.data());
	}

	vector<Twist<double>> exp_coords(6);

	exp_coords[0] << 0, 0,  1, 0, 0, 0;
	exp_coords[1] << 0, -1, 0, 89, 0, 0;
	exp_coords[2] << 0, -1, 0, 481, 0, 0;
	exp_coords[3] << 0, -1, 0, 906, 0, 0;
	exp_coords[4] << 0, 0,  1, -109, 0, 0;
	exp_coords[5] << 0, -1, 0, 1001, 0, 0;

	SE3d M = SE3d::Identity();
	M.block(0, 3, 3, 1) = Vector3d(0, -191, 1001);
	
	auto test = [&](auto&& solver) mutable{
		uint64_t count = 0;
				
		for (auto& theta : thetas) {
			SE3d goal;
			ForwardKinematic(goal, M, exp_coords, theta);

			DynVec  <double> out_theta;
			DynVec<double>   init_guess = DynVec<double>::Zero(6);
			InverseKinematic(out_theta,
				M, exp_coords,
				goal, init_guess,
				solver);

			SE3d cur = SE3d::Identity();
			ForwardKinematic(cur, M, exp_coords, out_theta);

			if (double res = LogMapSE3Tose3(InverseSE3(cur) * goal).squaredNorm(); res < 1e-5) {
				++count;
			}

		}

		cout << count / (double)N << endl;
	};
	
	
	LNSolver<double>   solver;
	DLSSolver<double>  dls_solver(1e-5);
	ADLSSolver<double> adls_solver;
	SDLSSolver<double> sdls_solver;
	adls_solver.SetFactor(1e-5);
	
	test(solver);
	test(dls_solver);
	test(adls_solver);
	test(sdls_solver);
	

}