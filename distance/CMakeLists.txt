cmake_minimum_required(VERSION 3.3)

project(dist)

option(CXX "enable C++ compilation" ON)
enable_language(CXX)

set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR})
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release")
endif()
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")

file(GLOB MISC misc/*.hpp misc/*.cpp)
file(GLOB DIST distance_algorithms/*.hpp distance_algorithms/*.cpp)
file(GLOB EXTER external/*.hpp external/*.cpp)
file(GLOB CYCLE cycle/*.hpp cycle/*.cpp)
file(GLOB HEUR heur/*.hpp heur/*.cpp)

set(CMAKE_CXX_STANDARD 11)

if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
        set(CXX_FILESYSTEM_LIBRARIES "stdc++fs")
endif()

#################################################################################################
# dist
#################################################################################################

find_package(OpenMP REQUIRED)
add_executable(${CMAKE_PROJECT_NAME} main_dist.cpp ${MISC} ${DIST} ${EXTER} ${CYCLE})
target_compile_options(dist PRIVATE -Wall PRIVATE "${OpenMP_CXX_FLAGS}")
target_link_libraries(dist PUBLIC ${CXX_FILESYSTEM_LIBRARIES} PRIVATE "${OpenMP_CXX_FLAGS}")

#################################################################################################
# dist_test
#################################################################################################

enable_testing()
find_package(OpenMP REQUIRED)
file(GLOB TEST quickcheck/quickcheck/*.hh)
add_executable(dist_test test.cpp ${MISC} ${DIST} ${TEST} ${EXTER} ${CYCLE} ${HEUR})
target_compile_options(dist_test PRIVATE -Wall PRIVATE "${OpenMP_CXX_FLAGS}")
target_link_libraries(dist_test PUBLIC ${CXX_FILESYSTEM_LIBRARIES} PRIVATE "${OpenMP_CXX_FLAGS}")

add_test(NAME dist_test COMMAND dist_test)

#################################################################################################
# decomposition
#################################################################################################

find_package(OpenMP REQUIRED)
add_executable(dec main_dec.cpp ${MISC} ${CYCLE} ${HEUR})
target_compile_options(dec PRIVATE -Wall PRIVATE "${OpenMP_CXX_FLAGS}")
target_link_libraries(dec PUBLIC ${CXX_FILESYSTEM_LIBRARIES} PRIVATE "${OpenMP_CXX_FLAGS}")

#################################################################################################


if(${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_CURRENT_SOURCE_DIR})
    include(FeatureSummary)
    feature_summary(WHAT ALL)
endif()
