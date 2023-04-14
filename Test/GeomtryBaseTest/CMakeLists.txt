cmake_minimum_required (VERSION 3.8)

function(AddGeometryTest arg)
	add_executable(${ARGV0} ${ARGV1})

	set_property(TARGET		   ${ARGV0} PROPERTY CXX_STANDARD 20)
	set_property(TARGET		   ${ARGV0} PROPERTY LINKER_LANGUAGE CXX)
	target_include_directories(${ARGV0} PRIVATE ${CommonInclude})
	target_link_libraries(	   ${ARGV0} PRIVATE ${CommonLib})

	add_test(GeomtryBaseTest ${ARGV0})
endfunction()

#______________________ Bary Centric Test ______________________________________________________
AddGeometryTest(BarycentricTest "TestBarycentricCoords.cpp")

#_______________________ Square Distance Test Box with Point ___________________________________
AddGeometryTest(SqrDistTest "TestSquareDist.cpp")


