/**
 *  @file  	gconcept.hpp
 *  @brief 	an eigen extending template constrain with cpp20.
 *  @author Gue Chen<guechen@buaa.edu.cn>
 *  @date 	July 3rd, 2022
 **/
#ifndef __G_CONCEPT_HPP
#define __G_CONCEPT_HPP

#include <Eigen/Core>

#include <concepts>

namespace GComponent{
template<typename T, typename U>
concept EigenSameRows = T::RowsAtCompileTime == Eigen::Dynamic || U::RowsAtCompileTime == Eigen::Dynamic || T::RowsAtCompileTime == U::RowsAtCompileTime;

template<typename T, typename U>
concept EigenSameCols = T::ColsAtCompileTime == Eigen::Dynamic || U::ColsAtCompileTime == Eigen::Dynamic || T::ColsAtCompileTime == U::ColsAtCompileTime;

template<typename T, typename U>
concept EigenSameSize = EigenSameRows<T, U> && EigenSameCols<T, U>;

template<typename Derived, int RowAtCompile, int ColAtComile>
concept MatConvertible = EigenSameSize<Eigen::Matrix<typename Derived::Scalar, RowAtCompile, ColAtComile>, Derived>;
// convertible_to<Eigen::Matrix<typename Derived::Scalar, ColAtCompile, RowAtComile>, Derived> // not effect cause : Matrix(MatrixBase<Derived>)
template<typename Derived, int RowAtCompile>
concept VecConvertible = MatConvertible<Derived, RowAtCompile, 1>;

template<typename T, typename U>
concept ScalarSame = std::same_as<typename T::Scalar, typename U::Scalar>;

template<typename T, typename Derived>
concept ScalarEquivalence = std::convertible_to<T, typename Derived::Scalar> || std::same_as<T, typename Derived::Scalar>;

template<typename Derived, typename OtherDerived>
concept MatScalarEquivalence = std::same_as<typename Derived::Scalar, typename OtherDerived::Scalar>;

// Special concepts for most usually types
// matrix concepts
template<typename Derived>
concept Mat3Convertible = MatConvertible<Derived, 3, 3>;
template<typename Derived>
concept Mat4Convertible = MatConvertible<Derived, 4, 4>;
template<typename Derived>
concept Mat6Convertible = MatConvertible<Derived, 6, 6>;
// vector concepts
template<typename Derived>
concept Vec3Convertible = VecConvertible<Derived, 3>;
template<typename Derived>
concept Vec4Convertible = VecConvertible<Derived, 4>;
template<typename Derived>
concept Vec6Convertible = VecConvertible<Derived, 6>;

template<typename ConvexType>
concept IsConvex = requires(ConvexType c) { {c.Support(Eigen::Vector3f{})}->std::convertible_to<Eigen::Vector3f>; };
} // namespace GComponent

#endif // !__G_CONCEPT_HPP
