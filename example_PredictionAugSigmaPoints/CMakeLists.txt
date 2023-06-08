cmake_minimum_required(VERSION 2.8 FATAL_ERROR)

add_definitions(-std=c++11)

#set(CXX_FLAGS "-Wall")
#"-Wextra" "-ggdb" "-O0" "-Wunreachable-code"
set(CXX_FLAGS "-Wall" "-pedantic")
set(CMAKE_CXX_FLAGS, "${CXX_FLAGS}")

project(example_PredictionAugSigmaPoints)

find_package(OpenCV 4.1 REQUIRED)
find_package(Iconv REQUIRED)
find_package(Eigen3 3.4.0 REQUIRED NO_MODULE)
if (TARGET Eigen3::Eigen)
  # Use the imported target
endif (TARGET Eigen3::Eigen)
find_package(Qt5Core REQUIRED)
find_package(Qt5Test REQUIRED)
find_package(Qt5Widgets REQUIRED)
find_package(Qt5 REQUIRED COMPONENTS Concurrent)
include_directories(${Eigen_INCLUDE_DIRS})
#include_directories(include)    do not have include in this folder
include_directories("${CMAKE_SOURCE_DIR}/jsoncpp"
    "${CMAKE_SOURCE_DIR}/jsoncpp/json"
    "${CMAKE_SOURCE_DIR}/jsonreader")
include_directories(${OpenCV_INCLUDE_DIRS})
link_directories(${OpenCV_LIBRARY_DIRS})
add_definitions(${OpenCV_DEFINITIONS})
find_package(Ceres REQUIRED)

# Executable for create matrix exercise
add_executable (example_PredictionAugSigmaPoints main.cpp ukf.cpp)
target_link_libraries(example_PredictionAugSigmaPoints Qt5::Test)
target_link_libraries(example_PredictionAugSigmaPoints Qt5::Concurrent)
target_link_libraries (example_PredictionAugSigmaPoints ${OpenCV_LIBS})
