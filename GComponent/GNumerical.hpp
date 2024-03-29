#ifndef _GNUMERICAL_HPP
#define _GNUMERICAL_HPP

#include <array>
#include <random>

namespace GComponent {

using std::array;

constexpr double MyDegree = 180.0;
constexpr double MyPI     = 3.1415926535897932384626433832795028841971693993751058209749445923078164062862089986280;

//template<unsigned GENNUM, unsigned NUM>
//array<array<double, NUM>, GENNUM> GenUniformRandoms(double, double);

template<class _Scaler, unsigned GENNUM, unsigned NUM>
array<array<_Scaler, NUM>, GENNUM> GenUniformRandoms(_Scaler lowerBound, _Scaler upperBound)
{
    array<array<_Scaler, NUM>, GENNUM> NUMSS;
    static std::random_device rd;
    static std::mt19937_64 gen(rd());
    std::uniform_real_distribution<> dis(lowerBound, upperBound);
    for (auto& nums : NUMSS) for (auto& num : nums)
    {
        num = dis(gen);
    }
    return NUMSS;
}

template<class _Scaler, unsigned NUM>
array<_Scaler, NUM>GenUniformRandom(_Scaler lowerBound, _Scaler upperBound)
{
    return GenUniformRandoms<_Scaler, 1, NUM>(lowerBound, upperBound)[0];
}

template<class Scaler>
Scaler Clamp(Scaler val, Scaler lower_bound, Scaler upper_bound)
{
    if(val <= lower_bound)
    {
        return lower_bound;
    }
    else if(val <upper_bound)
    {
        return val;
    }
    else
    {
        return upper_bound;
    }
}

inline double ToStandarAngle(double val)
{
    while(val < -MyPI) val += 2 * MyPI;
    while(val > MyPI)  val -= 2 * MyPI;
    return val;
}

constexpr inline double RadiusToDegree(double val)
{
    return val * MyDegree / MyPI;
}

constexpr inline double  DegreeToRadius(double val)
{
    return val * MyPI / MyDegree;
}

};

#endif // !_GNUMERICAL
