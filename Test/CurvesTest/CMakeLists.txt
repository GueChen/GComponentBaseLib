# mex function is too dificult
# add_subdirectory("MATLAB_Test")

function(AddCurvesTest arg)
	add_executable(${ARGV0} ${ARGV1})

	set_property(TARGET		   ${ARGV0} PROPERTY CXX_STANDARD 20)
	set_property(TARGET		   ${ARGV0} PROPERTY LINKER_LANGUAGE CXX)
	target_include_directories(${ARGV0} PRIVATE ${CommonInclude})
	target_link_libraries(	   ${ARGV0} PRIVATE ${CommonLib})

	#add_test(AddBaseCurvesTest ${ARGV0})
endfunction()

AddCurvesTest(CubicCurvesFunction "TestCubicCurvesFunction.cpp")