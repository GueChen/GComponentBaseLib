cmake_minimum_required (VERSION 3.8)

add_executable(SegmentTest     "TestNearestSegment.cpp")
add_executable(TriangleTest    "TestNearestTriangle.cpp")
add_executable(TetrahedronTest "TestNearestTetrahedron.cpp")

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

add_test(NearestSimplexTest SegmentTest TriangleTest TetrahedronTest)