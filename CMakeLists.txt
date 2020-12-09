cmake_minimum_required(VERSION 3.1)
project(monte-carlo-geometry-processing)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
set(CMAKE_CXX_STANDARD 17)

option(LIBIGL_WITH_OPENGL "Use OpenGL" ON)
option(LIBIGL_WITH_OPENGL_GLFW "Use GLFW" ON)
option(LIBIGL_WITH_OPENGL_GLFW_IMGUI "Use IMGUI" ON)

find_package(LIBIGL REQUIRED QUIET)

file(GLOB SRCFILES *.cpp)
add_executable(${PROJECT_NAME} main.cpp)
target_link_libraries(${PROJECT_NAME} igl::core igl::opengl igl::opengl_glfw igl::opengl_glfw_imgui)
