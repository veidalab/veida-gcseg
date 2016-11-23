cmake_minimum_required(VERSION 2.8)
 
PROJECT(GC3Dslices)
 
find_package(VTK REQUIRED)
include(${VTK_USE_FILE})
 
add_executable(GC3Dslices MACOSX_BUNDLE GC3Dslices)
 
if(VTK_LIBRARIES)
  target_link_libraries(GC3Dslices ${VTK_LIBRARIES})
else()
  target_link_libraries(GC3Dslices vtkHybrid vtkWidgets)
endif()