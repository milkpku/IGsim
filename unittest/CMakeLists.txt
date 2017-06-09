cmake_minimum_required(VERSION 2.8)

PROJECT(unit_test)

if (WIN32)
  include_directories(C:/VC_env/include/)
  include_directories(/usr/include/ /usr/local/include/)
endif(WIN32)

find_package(Eigen)
include_directories(${EIGEN_INCLUDE_DIRS})

find_package(libigl)
include_directories(${LIBIGL_INCLUDE_DIR})

enable_testing()
find_package(GTest)
if (GTEST_FOUND)
else()
  set(GTEST_BOTH_LIBRARIES C:/VC_env/lib64/gtest.lib
                           C:/VC_env/lib64/gtest_main.lib)
endif (GTEST_FOUND)
include_directories(${GTEST_INCLUDE_DIRS})


file(GLOB test_list *.cpp)
add_executable(runUnitTests ${test_list})


target_link_libraries(runUnitTests ${GTEST_BOTH_LIBRARIES})
add_test(runUnitTests runUnitTests)
