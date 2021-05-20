cmake_minimum_required(VERSION 3.10)
project(SquareMindCpp LANGUAGES C CXX)

if(POLICY CMP0072)
  cmake_policy(SET CMP0072 NEW)
endif()

# Default Build Type to be Release
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release")
endif()

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_LIST_DIR}/build)

# C++14 for dependencies and ABI compatibility
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# System dependencies
set(CMAKE_FOLDER "third_party")
find_package(OpenGL REQUIRED QUIET)
find_package(Threads REQUIRED QUIET)

# Third party dependencies
add_subdirectory(third_party/eigen)
add_subdirectory(third_party/glad)
add_subdirectory(third_party/glfw)
add_subdirectory(third_party/imgui)
add_subdirectory(third_party/libigl)

# build our app with C++17
set(CMAKE_CXX_STANDARD 17)
add_subdirectory(src)
