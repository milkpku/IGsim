cmake_minimum_required(VERSION 2.8)

PROJECT(unit_test)

find_package(Eigen)
include_directories(${EIGEN_INCLUDE_DIRS})

find_package(libigl)
include_directories(${LIBIGL_INCLUDE_DIR})

enable_testing()
set(GTEST_ROOT ${VC_ENV_DIR}/lib64/)
find_package(GTest)
if (GTEST_FOUND)
else()
  set(GTEST_BOTH_LIBRARIES ${VC_ENV_DIR}/lib64/gtest.lib
                           ${VC_ENV_DIR}/gtest_main.lib
                           ${VC_ENV_DIR}/gtestd.lib
                           ${VC_ENV_DIR}/gtest_maind.lib)
endif (GTEST_FOUND)
include_directories(${GTEST_INCLUDE_DIRS})

file(GLOB test_list *.cpp)
add_executable(runUnitTests ${test_list})

set(SHARED_PATH ${CMAKE_SOURCE_DIR}/shared)
add_definitions("-DSHARED_PATH=\"${SHARED_PATH}\"")

target_link_libraries(runUnitTests ${GTEST_BOTH_LIBRARIES})
add_test(runUnitTests runUnitTests)
