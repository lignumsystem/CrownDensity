# cmake file to compile CrownDensity
# For Unix Makefile build system:
# mkdir build
# cd  build
# cmake ..
# make
# For Xcode IDE:
# mkdir xcode
# cd xcode
# cmake .. -G Xcode
# open crowndens.xcodeproj   #Opens Xcode IDE
# Build the 'crowndens' Product in Xcode. See also Product (menu bar) -> Scheme is set to 'crowndens' to allow Run to debug.
# Copy necessary *.fun files and *.txt parameter files to xcode/Debug where 'crowndensÂ´ is located.
# Otherwise hard coded files are not found in the program.
# Set command line parameters for 'crowndens' in Product (menu bar) -> Scheme -> Edit Schmeme -> Arguments.
# Divide the command line into practical parts for debugging from Arguments -> '+'.
cmake_minimum_required(VERSION 3.12)
# cmake 3.23 for Xcode 13 recommends the policy CMP0114 (warning output)
cmake_policy(SET CMP0114 NEW)
project(crowndens DESCRIPTION "Single tree in homogeneous forest")

# Lignum core projects as an ExternalProject. 'lignum-core' target is compiled first
# Note: lignum-core  has to be configured with cmake first. See lignum-core/CMakeLists.txt
include(ExternalProject)
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
include(../Graphics/lib/Visual.cmake)

find_package(HDF5 COMPONENTS C CXX HL REQUIRED)
find_package(Qt4 COMPONENTS QtOpenGL QtXml QtGui QtCore REQUIRED)
find_package(OpenGL REQUIRED)
# cmake 3.23 output erroneusly prints 'glut not found'. It is found. 
find_package(GLUT REQUIRED)

set(SOURCE_FILES main.cc
globalvariables.cc
pine-em98.cpp
../LignumForest/branchfunctor.cc
../LignumForest/src/metabolism.cc
src/radiation.cc
../LignumForest/src/space.cc
)

# Rebuild L system compiler with custom command for 'l2c'
add_custom_command(
OUTPUT ${CMAKE_CURRENT_LIST_DIR}/pine-em98.cpp
DEPENDS ${CMAKE_CURRENT_LIST_DIR}/pine-em98.L l2c
COMMAND ${CMAKE_CURRENT_LIST_DIR}/../LEngine/bin/l2c ${CMAKE_CURRENT_LIST_DIR}/pine-em98.L  ${CMAKE_CURRENT_LIST_DIR}/pine-em98.cpp
)

# Main target, the executable 'crowndens'.
add_executable(crowndens ${SOURCE_FILES})
target_include_directories(crowndens PUBLIC include
${CMAKE_CURRENT_LIST_DIR}/../LignumForest
${CMAKE_CURRENT_LIST_DIR}/../LignumForest/include
${CMAKE_CURRENT_LIST_DIR}/../stl-voxelspace/include
${CMAKE_CURRENT_LIST_DIR}/../Pine
${CMAKE_CURRENT_LIST_DIR}/../XMLTree
${CMAKE_CURRENT_LIST_DIR}/../LEngine/include
)

# This add_dependencies for crowndens seem to generate time stamp based necessary automatic rebuilds
# of 'crowndens' and 'lignum-core' projects for Unix Makefile system, i.e. if a file is changed then
# that file and files depending on it will be recompiled.
# Note: 'make clean' is required first before 'make'.
#
# Xcode project created for 'crowndens' can compile and rebuild 'crowndens' (Product in Xcode) when needed.
# However 'lignum-core' dependency (Targets in Xcode) is compiled once ignoring after
# that time stamp based changes in Lignum core project files when rebuilding crowndens.
add_dependencies(crowndens lignum-core)

target_link_libraries(crowndens
HDF5::HDF5 Qt4::QtOpenGL Qt4::QtXml Qt4::QtGui Qt4::QtCore
GLUT::GLUT OpenGL::GLU OpenGL::GL
Visual L voxel sky LGM cxxadt
)

install(TARGETS crowndens DESTINATION "${CMAKE_CURRENT_LIST_DIR}")
