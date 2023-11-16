cmake_minimum_required(VERSION 3.0)
project(voronoiOpt)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS "-fopenmp -DBOOST_LOG_DYN_LINK")
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

set(sys_include_directories
        /usr/local/include/
        /usr/include/
        ${Python_INCLUDE_DIRS}
)

set(3rd_include_directories
        ${CMAKE_SOURCE_DIR}/3rd/eigen-3.4.0
        ${CMAKE_SOURCE_DIR}/3rd/libigl-2.4.0/include
        ${CMAKE_SOURCE_DIR}/3rd/spdlog-1.12.0/include
)

message(STATUS "3rd_include_directories: " ${3rd_include_directories})

include_directories(
        ${3rd_include_directories}
        ${sys_include_directories}
        ${CMAKE_SOURCE_DIR}/CVT/include
        ${CMAKE_SOURCE_DIR}/mma
        ${CMAKE_SOURCE_DIR}/src
)

link_directories(
        /usr/local/lib
        /usr/lib
        ${CMAKE_SOURCE_DIR}/lib/Release
)

find_package(TBB REQUIRED)
find_package(SuiteSparse REQUIRED)
find_package(Boost COMPONENTS filesystem log log_setup program_options REQUIRED)
find_package(OpenMesh REQUIRED)
find_package(CGAL REQUIRED)
find_package(OpenMP REQUIRED)
find_package(Threads REQUIRED)


add_executable(work ${CMAKE_SOURCE_DIR}/src/main.cpp)
target_link_libraries(work PUBLIC CBN_dev CVT mma_mma TBB::tbb ${SUITESPARSE_LIBRARIES}
        CGAL::CGAL OpenMeshCore geogram geogram_gfx ${Boost_LIBRARIES} ${Boost_LOG_LIBRARY}
        ${Boost_LOG_SETUP_LIBRARY} fmt gmp mpfr ${CMAKE_THREAD_LIBS_INIT})
