cmake_minimum_required(VERSION 2.8)
project(pybind)

### Compilation flags ###
if (MSVC)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /std:c++11 /MP /bigobj /w") ## Enable parallel compilation
    set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})
else()
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-deprecated-declarations -Wno-unused-parameter -Wno-deprecated-register -Wno-return-type-c-linkage")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wno-deprecated-declarations -Wno-unused-parameter -Wno-deprecated-register -Wno-return-type-c-linkage")
    set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "../")
endif()

set(Python_ADDITIONAL_VERSIONS 3.4 3.5 3.6)
find_package(PythonLibs REQUIRED)
find_package(PythonInterp REQUIRED)

string(TOUPPER "${CMAKE_BUILD_TYPE}" U_CMAKE_BUILD_TYPE)
if (UNIX)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -fPIC")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fPIC")
    if (NOT ${U_CMAKE_BUILD_TYPE} MATCHES DEBUG)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fvisibility=hidden -flto")
    endif()
endif()

include_directories(${PYTHON_INCLUDE_DIR} include)

#include pybind
include_directories(C:/VC_env/include)
include_directories(D:/v-likma/playground/IGsim/include/)


### Prepare building environment
add_library(pysim SHARED
        python_shared.cpp
        py_sim.cpp
        py_doc.cpp)

set_target_properties(pysim PROPERTIES PREFIX "")
set_target_properties(pysim PROPERTIES LIBRARY_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR})

target_link_libraries(pysim ${IGSIM_LIBRARIES})

if (WIN32)
    if (MSVC)
    # Enforce size-based optimization and link time code generation on MSVC (~30% smaller binaries in experiments)
    set_target_properties(pysim PROPERTIES COMPILE_FLAGS "/Os /GL")
    set_target_properties(pysim PROPERTIES LINK_FLAGS "/LTCG")
    endif()

    # .PYD file extension on Windows
    set_target_properties(pysim PROPERTIES SUFFIX ".pyd")

    target_link_libraries(pysim ${PYTHON_LIBRARIES})
endif()