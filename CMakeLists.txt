cmake_minimum_required(VERSION 3.6)
project(IDO4PoissonEq)

set(CMAKE_CXX_STANDARD 11)

set(SOURCE_FILES
        inc/grid.h
        src/grid.c
        src/main.c)
#
#add_executable(IDO4PoissonEq ${SOURCE_FILES})

#execute_process(COMMAND make WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
add_custom_target(IDO4PoissonEq COMMAND make -C ${IDO4PoissonEq_SOURCE_DIR}
        CLION_EXE_DIR=${PROJECT_BINARY_DIR})