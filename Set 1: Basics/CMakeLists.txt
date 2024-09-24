cmake_minimum_required(VERSION 3.8)

project(set1)

set(SOURCE main.cpp)

set(LIBRARIES convertions.cpp)

add_library(LIBS STATIC ${LIBRARIES})

add_executable(set1 ${SOURCE})

target_link_libraries(set1 LIBS)