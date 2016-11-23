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
char* string_arr;

int main(int argc, char *argv[])
{

/* cout << "Getting File Information..." << endl;

	 vtkSmartPointer<vtkImageReader2> reader =
  vtkSmartPointer<vtkImageReader2>::New();
reader->SetFileName("Foot.raw");
reader->SetFileDimensionality(3);
reader->SetDataExtent(0, 255, 0, 255, 0, 255);
reader->SetDataSpacing(1.0,1.0,1.0);
reader->SetDataOrigin(0.0, 0.0, 0.0);
reader->SetDataScalarTypeToUnsignedChar();
reader->SetDataByteOrderToBigEndian();
reader->UpdateWholeExtent();
 


cout << "Extracting a slice..." << endl;
int dims=255;
  //extract a slice from the volume
  vtkSmartPointer<vtkExtractVOI> extractSlice =
      vtkSmartPointer<vtkExtractVOI>::New();
  extractSlice->SetInputConnection(reader->GetOutputPort());
  extractSlice->SetVOI(0,dims,0,dims,50,50);
  extractSlice->Update(); 

 // string_arr = "C:\\Users\\user\\Desktop\\MKyan\\readraw\\bin\\Foot.raw";
  img = extractSlice->GetOutput();*/
  seg1 = new segmentation("Foot.raw");
					seg1->segment();
/*
    // Create actors
  vtkSmartPointer<vtkImageActor> inputActor =
    vtkSmartPointer<vtkImageActor>::New();
  inputActor->GetMapper()->SetInputConnection(extractSlice->GetOutputPort());
 
  //  one render window
  vtkSmartPointer<vtkRenderWindow> renderWindow =
    vtkSmartPointer<vtkRenderWindow>::New();
  renderWindow->SetSize(600, 500);
 
  // And one interactor
  vtkSmartPointer<vtkRenderWindowInteractor> interactor =
    vtkSmartPointer<vtkRenderWindowInteractor>::New();
  interactor->SetRenderWindow(renderWindow);
 
  // Setup renderer
  vtkSmartPointer<vtkRenderer> Renderer =
    vtkSmartPointer<vtkRenderer>::New();
  renderWindow->AddRenderer(Renderer);
  Renderer->SetBackground(.6, .5, .4);
  Renderer->AddActor(inputActor);
  Renderer->ResetCamera();
  renderWindow->Render();
  interactor->Start(); */


  return EXIT_SUCCESS;
}
