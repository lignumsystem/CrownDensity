
# /Library/Developer/CommandLineTools/SDKs/MacOSX11.3.sdk/
cmake_minimum_required(VERSION 3.12)
project(crowndens DESCRIPTION "Single tree in homogeneous forest")

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
find_package(GLUT REQUIRED)

set(SOURCE_FILES main.cc
globalvariables.cc
pine-em98.cpp
../LignumForest/branchfunctor.cc
../LignumForest/src/metabolism.cc
src/radiation.cc
../LignumForest/src/space.cc
)

# add_library(sky STATIC IMPORTED GLOBAL)
# add_library(L STATIC IMPORTED GLOBAL)
# add_library(voxel STATIC IMPORTED GLOBAL)
# add_library(LGM STATIC IMPORTED GLOBAL)
# add_library(cxxadt STATIC IMPORTED GLOBAL)

# add_directory("${CMAKE_CURRENT_LIST_DIR}/../Firmament")
# add_directory("${CMAKE_CURRENT_LIST_DIR}/../stl-lignum")
# add_directory("${CMAKE_CURRENT_LIST_DIR}/../stl-voxelspace")
# add_directory("${CMAKE_CURRENT_LIST_DIR}/../LEngine")

# include_directories(${CMAKE_CURRENT_LIST_DIR}
# ${CMAKE_CURRENT_LIST_DIR}/include
# ${CMAKE_CURRENT_LIST_DIR}/../LignumForest
# ${CMAKE_CURRENT_LIST_DIR}/../LignumForest/include
# ${CMAKE_CURRENT_LIST_DIR}/../c++adt/include
# ${CMAKE_CURRENT_LIST_DIR}/../Firmament/include
# ${CMAKE_CURRENT_LIST_DIR}/../stl-lignum/include
# ${CMAKE_CURRENT_LIST_DIR}/../stl-voxelspace/include
# ${CMAKE_CURRENT_LIST_DIR}/../LEngine/include
# ${CMAKE_CURRENT_LIST_DIR}/../Graphics
# ${CMAKE_CURRENT_LIST_DIR}/../Pine
# ${CMAKE_CURRENT_LIST_DIR}/../XMLTree
# /opt/local/include
# )

add_custom_command(
OUTPUT ${CMAKE_CURRENT_LIST_DIR}/pine-em98.cpp
DEPENDS ${CMAKE_CURRENT_LIST_DIR}/pine-em98.L l2c
COMMAND ${CMAKE_CURRENT_LIST_DIR}/../LEngine/bin/l2c ${CMAKE_CURRENT_LIST_DIR}/pine-em98.L  ${CMAKE_CURRENT_LIST_DIR}/pine-em98.cpp
)

add_executable(crowndens ${SOURCE_FILES})
target_include_directories(crowndens PUBLIC include
${CMAKE_CURRENT_LIST_DIR}/../LignumForest
${CMAKE_CURRENT_LIST_DIR}/../LignumForest/include
${CMAKE_CURRENT_LIST_DIR}/../stl-voxelspace/include
${CMAKE_CURRENT_LIST_DIR}/../Pine
${CMAKE_CURRENT_LIST_DIR}/../XMLTree
${CMAKE_CURRENT_LIST_DIR}/../LEngine/include
)
target_link_libraries(crowndens HDF5::HDF5 Qt4::QtOpenGL Qt4::QtXml Qt4::QtGui Qt4::QtCore GLUT::GLUT OpenGL::GLU OpenGL::GL  Visual L voxel sky LGM cxxadt)

