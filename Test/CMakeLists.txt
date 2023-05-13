
find_package(GTest REQUIRED)
set(CommonInclude ${GTest_INCLUDE_DIR} ${common_lib_INCLUDE_DIR})
set(CommonLib     common_lib GTest::gtest GTest::gtest_main GTest::gmock GTest::gmock_main)

enable_testing()

add_subdirectory("NearestSimplexTest")
add_subdirectory("GeomtryBaseTest")
add_subdirectory("CurvesTest")
add_subdirectory("KinematicTest")
add_subdirectory("RotationTest")