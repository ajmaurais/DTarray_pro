
cmake_minimum_required(VERSION 3.15)
project(DTArray_pro)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/cmake")

#set(${DTARRAY_PRO_TARGET} DTarray)

set(CMAKE_CXX_STANDARD 11)

include_directories(include)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DSHARE_DIR=\\\"${SHARE_DIR}\\\"")

#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DGIT_VERSION="\"4bd5ba1\"" -DGIT_DATE="\"Tue Mar 10 18:32:25 2020\""")

SET(CMAKE_BUILD_TYPE Debug)

add_executable(DTarray
        src/calcMW.cpp
        src/dtafilter.cpp
        src/locReport.cpp
        src/main.cpp
        src/params.cpp
        src/saintOutput.cpp )

#add utils library
include(FindPkgConfig)
find_package(utils REQUIRED)
target_include_directories(DTarray PUBLIC ${UTILS_INCLUDE_DIR})
target_link_libraries(DTarray ${UTILS_LIB})

