#pragma once
#include <array>
namespace GComponent {
	using std::array;
	
	template<unsigned GENNUM, unsigned NUM>
	array<array<double, NUM>, GENNUM> GenUniformRandoms(double, double);

	template<unsigned NUM>
	array<double, NUM>GenUniformRandom(double lowerBound, double upperBound)
	{
		return GenUniformRandoms<1, NUM>(lowerBound, upperBound)[0];
	}

	template<unsigned GENNUM, unsigned NUM>
	array<array<double, NUM>, GENNUM> GenUniformRandoms(double lowerBound, double upperBound)
	{
		array<array<double, NUM>, GENNUM> NUMSS;
		static std::random_device rd;
		static std::mt19937_64 gen(rd());
		std::uniform_real_distribution<> dis(lowerBound, upperBound);
		for (auto& nums : NUMSS)
		{
			for (auto& num : nums)
			{
				num = dis(gen);
			}
		}
		return NUMSS;
	}

	
};