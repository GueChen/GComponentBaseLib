cmake_minimum_required (VERSION 3.8)

set(TARGET_NAME common_lib)

file(GLOB_RECURSE HEADER_FILES "*.h" "*.hpp")
file(GLOB_RECURSE SOURCE_FILES "*.c" "*.cpp")

find_package(glad REQUIRED)
find_package(GLM REQUIRED)
find_package(glfw3 REQUIRED)
find_package(Eigen3 REQUIRED)

add_library(${TARGET_NAME} ${SOURCE_FILES} ${HEADER_FILES})

set_property(TARGET ${TARGET_NAME} PROPERTY CXX_STANDARD 20)
set_property(TARGET ${TARGET_NAME} PROPERTY LINKER_LANGUAGE CXX)

target_include_directories(${TARGET_NAME} PUBLIC
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