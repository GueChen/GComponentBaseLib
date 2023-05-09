cmake_minimum_required (VERSION 3.8)

# help Add Testing function
# ARGV0 -- executable name
# ARGV1 -- compile files
function(AddKinematicTest arg)

	add_executable(${ARGV0} ${ARGV1})

	set_property(TARGET		   ${ARGV0} PROPERTY CXX_STANDARD 20)
	set_property(TARGET		   ${ARGV0} PROPERTY LINKER_LANGUAGE CXX)
	target_include_directories(${ARGV0} PRIVATE ${CommonInclude})
	target_link_libraries(	   ${ARGV0} PRIVATE ${CommonLib})

	add_test(KinematicTest ${ARGV0})
endfunction()

AddKinematicTest(TestIK4UR5 TestIKforUR5.cpp)