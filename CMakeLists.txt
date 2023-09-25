# -----------------------------------------
# CMAKE file to compile 'crowndens' and
# install it in CrownDensity directory
# -----------------------------------------
# For Unix Makefile build system for Debug:
# mkdir debug
# cd  debug
# cmake .. -DCMAKE_BUILD_TYPE=Debug
# Optionally explicitely set c++ compiler:
# cmake ..  -DCMAKE_BUILD_TYPE=Debug  -DCMAKE_C_COMPILER=clang -DCMAKE_CXX_COMPILER=clang++ -DCMAKE_CXX_FLAGS=-stdlib=libc++
# make install
# ------------------------------------------
# For Unix Makefile build system for Release:
# mkdir release
# cd release
# cmake .. -DCMAKE_BUILD_TYPE=Release
# make install
# ---------------------------------------------
# 1. For Xcode IDE first create Xcode project file:
# mkdir xcode
# cd xcode
# cmake .. -G Xcode
# Open Xcode IDE:
# open crowndens.xcodeproj
# 2. Build the 'crowndens' Product in Xcode for debugging.
# See also Xcode -> Product (in menu bar) -> Scheme is set to 'crowndens' to allow Run to debug.
# Copy necessary *.fun files and *.txt parameter files to xcode/Debug where 'lignum-forest´ is located.
# Otherwise hard coded files are not found in the program.
# Set command line parameters for 'lignum-forest' in Product (menu bar) -> Scheme -> Edit Schmeme -> Arguments.
# Divide the command line into practical parts for debugging from Arguments -> '+'.
# -------------------------------------------------------------------------------
# Graphviz
# Generate Graphviz files showing the dependency graph between the targets and libraries in a project:
# mkdir release
# cd release
# mkdir graphviz
# cmake .. -DCMAKE_BUILD_TYPE=Release -graphviz release/CrownDensity.dot
# Create pdf file of project dependencies:
# dot -Kneato -Goverlap=scale -Tpdf -o CrownDensity.pdf graphviz/CrownDensity.dot
# --------------------------------------------------------------------------------
cmake_minimum_required(VERSION 3.12)
# cmake 3.23 for Xcode 13 recommends the policy CMP0114 (warning output)
cmake_policy(SET CMP0114 NEW)
project(crowndens DESCRIPTION "Single tree in homogeneous forest")
message("Project " ${CMAKE_PROJECT_NAME})

# CMake can generate Graphviz files showing the dependencies between the targets in a project
include(CMakeGraphVizOptions)
set(GRAPHVIZ_CUSTOM_TARGETS TRUE)
# Lignum core projects as an ExternalProject. 'lignum-core' target is compiled first
# Note: lignum-core  has to be configured with cmake first. See lignum-core/CMakeLists.txt
include(ExternalProject)

# Uncomment L-system to be used
# Comment out other L-systems
#set(LSYSTEMFILE pine-em98.L)
#set(LSYSTEMSRC  pine-em98.cpp)
#set(LSYSTEMFILE pine-em98-branch-A.L)
#set(LSYSTEMSRC  pine-em98-branch-A.cpp)
set(LSYSTEMFILE pine-em98-branch-B.L)
set(LSYSTEMSRC  pine-em98-branch-B.cpp)

ExternalProject_add(lignum-core
SOURCE_DIR ${CMAKE_CURRENT_LIST_DIR}/..
DOWNLOAD_COMMAND "")
# The '*.cmake' files are produced with 'cmake' configuration for 'lignum-core' (see lignum-core/CMakeLists.txt).  
include(../c++adt/lib/cxxadt.cmake)
include(../Firmament/lib/sky.cmake)
include(../stl-lignum/lib/LGM.cmake)
include(../stl-voxelspace/lib/voxel.cmake)
include(../LEngine/lib/L.cmake)
include(../LEngine/lib/l2c.cmake)
if(CMAKE_HOST_SYSTEM_NAME MATCHES "Darwin")
  message("Darwin (macOS): Compiling Visual")
  include(../Graphics/lib/Visual.cmake)
endif()

find_package(HDF5 COMPONENTS C CXX HL REQUIRED)
if(CMAKE_HOST_SYSTEM_NAME MATCHES "Darwin")
  message("SYSTEM: Apple platform")
  find_package(Qt4 COMPONENTS QtOpenGL QtXml QtGui QtCore REQUIRED)
  find_package(OpenGL REQUIRED) 
  # cmake 3.23 output erroneusly prints 'glut not found'. It is found. 
  find_package(GLUT REQUIRED)
elseif(CMAKE_HOST_SYSTEM_NAME MATCHES "Linux")
  message("SYSTEM: Linux platform")
  find_package(HDF5 COMPONENTS C CXX HL REQUIRED)
  find_package(Qt5 COMPONENTS Xml Core  REQUIRED)
else()
  message(${CMAKE_HOST_SYSTEM_NAME} " not supported")
  message(FATAL_ERROR "Exiting cmake")
endif()

#find_package(OpenGL REQUIRED)
# cmake 3.23 output erroneusly prints 'glut not found'. It is found. 
#find_package(GLUT REQUIRED)

set(SOURCE_FILES main.cc
  crowndensity_globals.cc
  globalvariables.cc
  ${LSYSTEMSRC}
  ../LignumForest/branchfunctor.cc
  ../LignumForest/src/metabolism.cc
  src/radiation.cc
  ../LignumForest/src/space.cc
)

# Rebuild L system compiler with custom command for 'l2c'
add_custom_command(
OUTPUT ${CMAKE_CURRENT_LIST_DIR}/${LSYSTEMSRC}
DEPENDS ${CMAKE_CURRENT_LIST_DIR}/${LSYSTEMFILE} l2c
COMMAND ${CMAKE_CURRENT_LIST_DIR}/../LEngine/bin/l2c ${CMAKE_CURRENT_LIST_DIR}/${LSYSTEMFILE} ${CMAKE_CURRENT_LIST_DIR}/${LSYSTEMSRC}
)

# Main target, the executable 'crowndens'.
add_executable(crowndens ${SOURCE_FILES})
target_include_directories(crowndens PUBLIC . include
${CMAKE_CURRENT_LIST_DIR}/../LignumForest
${CMAKE_CURRENT_LIST_DIR}/../LignumForest/include
${CMAKE_CURRENT_LIST_DIR}/../stl-voxelspace/include
${CMAKE_CURRENT_LIST_DIR}/../Pine
${CMAKE_CURRENT_LIST_DIR}/../XMLTree
${CMAKE_CURRENT_LIST_DIR}/../LEngine/include
)

# This add_dependencies for crowndens seem to generate time stamp based necessary automatic rebuilds
# of 'crowndens' and 'lignum-core' projects for Unix Makefile system, i.e. if a file is changed then
# that file and files depending on it will be recompiled with 'make'.
# Note: 'make clean' is required first before 'make'.
#
# Xcode project created for 'crowndens' can compile and rebuild 'crowndens' (Product in Xcode) when needed.
# However 'lignum-core' dependency (Targets in Xcode) is compiled once ignoring after
# that time stamp based changes in Lignum core project files when rebuilding crowndens.
add_dependencies(crowndens lignum-core)

if(CMAKE_HOST_SYSTEM_NAME MATCHES "Darwin")
  target_link_libraries(crowndens
    HDF5::HDF5 Qt4::QtOpenGL Qt4::QtXml Qt4::QtGui Qt4::QtCore
    GLUT::GLUT OpenGL::GLU OpenGL::GL
    Visual L voxel sky LGM cxxadt
  )
elseif(CMAKE_HOST_SYSTEM_NAME MATCHES "Linux")
  target_link_libraries(crowndens
    HDF5::HDF5 Qt5::Xml Qt5::Core
    L voxel sky LGM cxxadt
  )
  message("REMEMBER to type \"make install\" to install lignum-forest") 
else()
  message(${CMAKE_HOST_SYSTEM_NAME} " not supported")
  message(FATAL_ERROR "Exiting cmake")
endif()

install(TARGETS crowndens DESTINATION "${CMAKE_CURRENT_LIST_DIR}")
