cmake_minimum_required (VERSION 3.8)

set(TARGET_NAME common_lib)

file(GLOB_RECURSE HEADER_FILES "source/*.h" "source/*.hpp")
file(GLOB_RECURSE SOURCE_FILES "source/*.c" "source/*.cpp")
#____________________ Keep File Path Structrue ________________________________________
source_group(TREE "${CMAKE_CURRENT_SOURCE_DIR}" FILES ${HEADER_FILES} ${SOURCE_FILES})

find_package(glad REQUIRED)
find_package(GLM REQUIRED)
find_package(glfw3 REQUIRED)
find_package(Eigen3 REQUIRED)

add_library(${TARGET_NAME} ${SOURCE_FILES} ${HEADER_FILES})

set_property(TARGET ${TARGET_NAME} PROPERTY CXX_STANDARD z)
set_property(TARGET ${TARGET_NAME} PROPERTY LINKER_LANGUAGE CXX)

target_include_directories(${TARGET_NAME} PUBLIC
	source
	${GLM_INCLUDE_DIR}
	${GLFW3_INCLUDE_DIR}
	${EIGEN_INCLUDE_DIR}
)


target_link_libraries(${TARGET_NAME} PUBLIC
	glad::glad
	glm::glm
	glfw
	Eigen3::Eigen
)

target_include_directories(
	${TARGET_NAME} 
	PUBLIC
	"$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>"
)

add_subdirectory("Test")