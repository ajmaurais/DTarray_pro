
cmake_minimum_required(VERSION 3.15)
project(DTArray_pro)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/cmake")

set(CMAKE_CXX_STANDARD 11)

set(GIT_PRE_CONFIGURE_FILE ${CMAKE_CURRENT_SOURCE_DIR}/cmake/git_info.h.in)
set(GIT_POST_CONFIGURE_FILE ${CMAKE_CURRENT_BINARY_DIR}/include/git_info.h)

include_directories(include)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DSHARE_DIR=\\\"${SHARE_DIR}\\\"")

SET(CMAKE_BUILD_TYPE Debug)

add_executable(DTarray
        src/calcMW.cpp
        src/dtafilter.cpp
        src/locReport.cpp
        src/main.cpp
        src/params.cpp
        src/saintOutput.cpp)

#add utils library
include(FindPkgConfig)
find_package(utils REQUIRED)
target_include_directories(DTarray PUBLIC ${UTILS_INCLUDE_DIR}
        PUBLIC ${CMAKE_CURRENT_BINARY_DIR}/include)
target_link_libraries(DTarray ${UTILS_LIB})

#git tracking
option(TRACK_GIT "Should git tracking be enabeled?" ON)
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/include)
if(TRACK_GIT)
    include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/git_watcher.cmake)
    add_dependencies(DTarray check_git_repository)
    set_directory_properties(DIRECTORY APPEND PROPERTY ADDITIONAL_MAKE_CLEAN_FILES "git-state")
else()
    function(WriteBlankGitHeader)
        set(GIT_RETRIEVED_STATE "false")
        set(GIT_HEAD_SHA1 "GIT_NOT_FOUND")
        set(GIT_IS_DIRTY "false")
        set(GIT_LAST_COMMIT_DATE "GIT_NOT_FOUND")
        configure_file("${GIT_PRE_CONFIGURE_FILE}" "${GIT_POST_CONFIGURE_FILE}" @ONLY)
    endfunction()
    WriteBlankGitHeader()
endif()

