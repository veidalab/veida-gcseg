cmake_minimum_required(VERSION 2.8)
 
PROJECT(GC2D)
 
find_package(VTK REQUIRED)
include(${VTK_USE_FILE})
 
add_executable(GC2D MACOSX_BUNDLE GC2D)
 
if(VTK_LIBRARIES)
  target_link_libraries(GC2D ${VTK_LIBRARIES})
else()
  target_link_libraries(GC2D vtkHybrid vtkWidgets)
endif()