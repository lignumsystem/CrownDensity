######################################################################
# Automatically generated by qmake (2.00a) Thu Aug 17 15:33:51 2006
######################################################################
CONFIG -= app_bundle
CONFIG += qt
QT += xml
TEMPLATE = app
TARGET = crowndens 
INCLUDEPATH += . /opt/local/include include ../c++adt/include ../stl-lignum/include ../Firmament/include ../stl-voxelspace/include ../LEngine/include ../Pine ../XMLTree  ../Graphics
INCLUDEPATH += ../LignumForest ../LignumForest/include 
DEPENDPATH += $$INCLUDEPATH
LIBS += -L/opt/local/lib -L../c++adt/lib -L../stl-lignum/lib -L../Firmament/lib -L../LEngine/lib -L../stl-voxelspace/lib
LIBS += -lsky -lL -lvoxel -lLGM  -lcxxadt
LIBS += -lhdf5_cpp -lhdf5 -lz -ldl -lm
# Input
unix{
   system(../LEngine/bin/l2c pine-em98.L pine-em98.cpp){
     SOURCES += pine-em98.cpp
   }
   else{
     message("unix pine-em98.L failed")
   }
}


win32{
    system(..\LEngine\bin\l2c pine-em98.L pine-em98.cpp){
      SOURCES  += pine-em98.cpp
    }
    else{
      message("win32 pine-em98.L failed")
    }
}
     
macx:LIBS +=  -L../Graphics -lVisual  -framework GLUT -framework OpenGL
win32:CONFIG += console
HEADERS += 
#include/DiameterGrowth.h \
#           include/CrownDensityGlobals.h\
#           include/ScotsPine.h \
#           include/SomeFunctors.h include/RadiationCrownDens.h \
#           include/Space.h include/Palubicki_functors.h include/ByBranches.h
           
SOURCES +=  main.cc globalvariables.cc crowndensity_globals.cc\
           ../LignumForest/branchfunctor.cc \
           ../LignumForest/src/metabolism.cc src/radiation.cc \
           ../LignumForest/src/space.cc
#src/bybranches.cc
