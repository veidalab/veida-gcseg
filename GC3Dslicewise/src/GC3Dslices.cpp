#include <math.h>
#include<stdio.h>
#include <iostream>

#include <vtkSmartPointer.h>
#include <vtkExtractVOI.h>
#include <vtkImageData.h>
#include <vtkImageReader2.h>
#include <vtkInteractorStyleImage.h>
#include <vtkImageActor.h>
#include <vtkImageMapper3D.h>
#include <vtkImageThreshold.h>

#include <vtkRenderWindow.h>
#include <vtkRenderer.h>
#include <vtkCamera.h>
#include <vtkRenderWindowInteractor.h>

#include "graphcutAlgorithm.h"

segmentation * seg1;
vtkImageData* img;
char* string_arr0;
char* string_arr;
char* string_arr2;

int main(int argc, char *argv[])
{
   int marked = 0;
   string_arr0 = "C:\\Users\\user\\Desktop\\MKyan\\Tests\\Test3D\\bin\\Foot0.jpg";
   string_arr =  "C:\\Users\\user\\Desktop\\MKyan\\Tests\\Test3D\\bin\\Foot.jpg";
   string_arr2 = "C:\\Users\\user\\Desktop\\MKyan\\Tests\\Test3D\\bin\\Foot2.jpg";
  
     seg1 = new segmentation(string_arr0,string_arr,string_arr2,marked);
						seg1->segment();

  return EXIT_SUCCESS;
}
