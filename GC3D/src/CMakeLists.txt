cmake_minimum_required(VERSION 2.8)
 
PROJECT(GC3D)
 
find_package(VTK REQUIRED)
include(${VTK_USE_FILE})
 
add_executable(GC3D MACOSX_BUNDLE GC3D)
 
if(VTK_LIBRARIES)
  target_link_libraries(GC3D ${VTK_LIBRARIES})
else()
  target_link_libraries(GC3D vtkHybrid vtkWidgets)
endif()