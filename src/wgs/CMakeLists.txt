cmake_minimum_required(VERSION 3.3.0)

project(CAPG VERSION 1.0)

set(CMAKE_C_COMPILER gcc)
set(CMAKE_C_STANDARD 17)

# -DCMAKE_BUILD_TYPE=Debug for debug version
#set(CMAKE_BUILD_TYPE Release)	# this appears to be Release
set(CMAKE_DEBUG_POSTFIX d)	# append d to executable names to signal it is debug version

set(RMATH_PATH ".")

include(FindPackageHandleStandardArgs)

aux_source_directory(. DIR_SRCS)

# below modified from RStudio's FindLibR.cmake to find R include and lib path 
# https://github.com/rstudio/rstudio/blob/master/cmake/modules/FindLibR.cmake 
function(find_rlib)
    if(APPLE)
    
       find_library(LIBR_LIBRARIES R)
       
       if(LIBR_LIBRARIES MATCHES ".*\\.framework")
       set(LIBR_HOME "${LIBR_LIBRARIES}/Resources" CACHE PATH "R home directory")
          set(LIBR_INCLUDE_DIRS "${LIBR_HOME}/include" CACHE PATH "R include directory")
          set(LIBR_LIB_DIR "${LIBR_HOME}/lib" CACHE PATH "R lib directory")
          set(LIBR_EXECUTABLE "${LIBR_HOME}/bin/R" CACHE PATH "R executable")
       else()
          get_filename_component(_LIBR_LIBRARIES "${LIBR_LIBRARIES}" REALPATH)
          get_filename_component(_LIBR_LIBRARIES_DIR "${_LIBR_LIBRARIES}" PATH)
          set(LIBR_EXECUTABLE "${_LIBR_LIBRARIES_DIR}/../bin/R")
          execute_process(
             COMMAND ${LIBR_EXECUTABLE} "--slave" "--vanilla" "-e" "cat(R.home())"
                       OUTPUT_VARIABLE LIBR_HOME
          )
          set(LIBR_HOME ${LIBR_HOME} CACHE PATH "R home directory")
          set(LIBR_INCLUDE_DIRS "${LIBR_HOME}/include" CACHE PATH "R include directory")
          set(LIBR_LIB_DIR "${LIBR_HOME}/lib" CACHE PATH "R lib directory")
       endif()
    
    else()  # UNIX
    
    # find executable
          find_program(LIBR_EXECUTABLE R)
          if(LIBR_EXECUTABLE-NOTFOUND)
             message(STATUS "Unable to locate R executable")
          endif()
    
    # ask R for the home path
          if(NOT LIBR_HOME)
             execute_process(
                COMMAND ${LIBR_EXECUTABLE} "--slave" "--vanilla" "-e" "cat(R.home())"
                          OUTPUT_VARIABLE LIBR_HOME
             )
             if(LIBR_HOME)
               set(LIBR_HOME ${LIBR_HOME} CACHE PATH "R home directory")
             endif()
          endif()
    
     # ask R for the include dir
          if(NOT LIBR_INCLUDE_DIRS)
             execute_process(
                COMMAND ${LIBR_EXECUTABLE} "--slave" "--no-save" "-e" "cat(R.home('include'))"
                OUTPUT_VARIABLE LIBR_INCLUDE_DIRS
             )
             if(LIBR_INCLUDE_DIRS)
               set(LIBR_INCLUDE_DIRS ${LIBR_INCLUDE_DIRS} CACHE PATH "R include directory")
             endif()
          endif()
    
      # ask R for the lib dir
          if(NOT LIBR_LIB_DIR)
             execute_process(
                COMMAND ${LIBR_EXECUTABLE} "--slave" "--no-save" "-e" "cat(R.home('lib'))"
                OUTPUT_VARIABLE LIBR_LIB_DIR
             )
              if(LIBR_LIB_DIR)
               set(LIBR_LIB_DIR ${LIBR_LIB_DIR} CACHE PATH "R lib directory")           
             endif()
          endif()
    
    
       # sometimes lib and include are not under R.home()
       get_filename_component(_LIBR_EXECUTABLE "${LIBR_EXECUTABLE}" PATH)
       get_filename_component(_LIBR_EXECUTABLE_DIR "${_LIBR_EXECUTABLE}" PATH)
       #message(STATUS "Found R upper directory: ${_LIBR_EXECUTABLE_DIR}")
    
    endif()
endfunction()



set(C_OBJ ./align.c ./capg.c ./cmdline.c ./dnorm.c ./error.c ./fmax2.c ./io.c
    ./lgammacor.c ./nuc.c ./pchisq.c ./pick_reads.c ./qual.c ./sequence.c ./vcf.c
    ./bd0.c ./chebyshev.c ./cospi.c ./dpois.c ./fastq.c ./gamma.c ./lgamma.c
    ./mlutils.c ./order.c ./pgamma.c ./pnorm.c ./sam.c ./stirlerr.c)
set(C_OBJ_SIM ./simulator.c ./sequence.c ./nuc.c ./qual.c ./error.c ./fastq.c ./align.c ./io.c ./cmdline.c)

# compile option 
#add_compile_options(-std=c11 -D_POSIX_C_SOURCE=200112L -D_XOPEN_SOURCE=500 -O3)
#add_compile_options(-Wall -Wextra -pedantic)

add_executable(capg_wgs ${C_OBJ})
set_target_properties(capg_wgs PROPERTIES DEBUG_POSTFIX ${CMAKE_DEBUG_POSTFIX})
target_link_libraries(capg_wgs -lm) # math

find_rlib()
message(STATUS "Found R: ${LIBR_EXECUTABLE}")
message(STATUS "Found R home: ${LIBR_HOME}")
message(STATUS "Found R include directory: ${LIBR_INCLUDE_DIRS}")
message(STATUS "Found R lib directory: ${LIBR_LIB_DIR}")
if(LIBR_EXECUTABLE-NOTFOUND)
else()
	add_executable(capg_sim ${C_OBJ_SIM})
	set_target_properties(capg_sim PROPERTIES DEBUG_POSTFIX ${CMAKE_DEBUG_POSTFIX})

	target_include_directories(capg_sim PUBLIC ${LIBR_INCLUDE_DIRS})
	target_link_directories(capg_sim PUBLIC ${LIBR_LIB_DIR} ${RMATH_PATH})

	# Rmath
	find_library(RMATH_LIBRARY NAMES Rmath HINTS ${RMATH_PATH} ${LIBR_LIB_DIR} ${_LIBR_EXECUTABLE_DIR}/lib ${_LIBR_EXECUTABLE_DIR}/lib64)
	find_package_handle_standard_args(LibRmath DEFAULT_MSG RMATH_LIBRARY)
	if(RMATH_LIBRARY)
		target_link_libraries(capg_sim PUBLIC ${RMATH_LIBRARY})
	else()
		message(FATAL_ERROR "Failed to find libRmath !")
	endif()
	target_link_libraries(capg_sim PUBLIC -lm)
endif()
