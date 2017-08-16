cmake_minimum_required(VERSION 2.8)

project(IGsim)

set(VC_ENV_DIR C:/VC_env)

if (WIN32)
    include_directories(${VC_ENV_DIR}/include/)
    include_directories(/usr/include/ /usr/local/include/)
endif(WIN32)

add_subdirectory(unittest)
