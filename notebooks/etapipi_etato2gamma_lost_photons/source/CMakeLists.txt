cmake_minimum_required(VERSION 3.10 FATAL_ERROR)
project(PiPiX VERSION 1.0)
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()
set(CMAKE_CXX_FLAGS "-Wall")
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")
find_package(Boost COMPONENTS program_options REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})
find_package(kfcmd REQUIRED)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/env.sh.in env.sh)
file(GLOB SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp)
add_executable(kf-2pix-run ${SOURCES})
target_include_directories(kf-2pix-run
  PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:include>)
target_link_libraries(kf-2pix-run
  PUBLIC
  kfcmd::kfcmd_hypos
  ${Boost_LIBRARIES})
install(TARGETS kf-2pix-run DESTINATION bin)
install(FILES ${PROJECT_BINARY_DIR}/env.sh DESTINATION bin)
