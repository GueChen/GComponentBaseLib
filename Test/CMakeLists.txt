cmake_minimum_required (VERSION 3.8)

enable_testing()

find_package(GTest REQUIRED)
set(CommonInclude ${GTest_INCLUDE_DIR} ${common_lib_INCLUDE_DIR})
set(CommonLib     common_lib GTest::gtest GTest::gtest_main GTest::gmock GTest::gmock_main)

add_executable(SegmentTest     TestNearestSegment.cpp)
add_executable(TriangleTest    TestNearestTriangle.cpp)
add_executable(TetrahedronTest TestNearestTetrahedron.cpp)

set_property(TARGET		   SegmentTest PROPERTY CXX_STANDARD 20)
set_property(TARGET		   SegmentTest PROPERTY LINKER_LANGUAGE CXX)
target_include_directories(SegmentTest PRIVATE ${CommonInclude})
target_link_libraries(	   SegmentTest PRIVATE ${CommonLib})

set_property(TARGET		   TriangleTest PROPERTY CXX_STANDARD 20)
set_property(TARGET		   TriangleTest PROPERTY LINKER_LANGUAGE CXX)
target_include_directories(TriangleTest PRIVATE ${CommonInclude})
target_link_libraries(	   TriangleTest PRIVATE ${CommonLib})

set_property(TARGET		   TetrahedronTest PROPERTY CXX_STANDARD 20)
set_property(TARGET		   TetrahedronTest PROPERTY LINKER_LANGUAGE CXX)
target_include_directories(TetrahedronTest PRIVATE ${CommonInclude})
target_link_libraries(	   TetrahedronTest PRIVATE ${CommonLib})

add_test(AllTestsGeometry SegmentTest TriangleTest TetrahedronTest)

