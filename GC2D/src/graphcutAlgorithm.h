//Test2D

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <ctype.h>
#include <iostream>
#include<math.h>
#include <fstream>
#include <time.h>
#include "graphcut.h"

#include <vtkImageData.h>
#include <vtkSmartPointer.h>
#include <vtkImageReader2.h>
#include <vtkImageProperty.h>
#include <vtkExtractVOI.h>
#include <vtkDataArray.h>
#include <vtkDoubleArray.h>
#include <vtkJPEGWriter.h>
#include <vtkInteractorStyleImage.h>
#include <vtkImageActor.h>
#include <vtkImageMapper3D.h>
#include <vtkRenderWindow.h>
#include <vtkRenderer.h>
#include <vtkCamera.h>
#include <vtkRenderWindowInteractor.h>

using std::ofstream;
using std::ifstream;
using std::cout;


#define K 1000
#define Ki 1
//#define sigma2 1000
#define const1 100
#define const2 1000
#define distance 50

/* A simple bit-vector. Everyone's written their own-- here's mine.
   It's endian neutral so life is much easier to deal with on multiple
   platforms.  It probably would be more efficient to do the
   operations on words, but I don't want to get into the endian
   issue. Allows you to set / clear and query a value. All one needs
   for sets.

   Walter Bell (wbell@cs.cornell.edu) 
*/
class segmentation 
{
public:
  /* Create a new bit vector with n indices. All bit values will be
     initialized to clear.
  */
//vtkImageData* finalimg;
// vtkImageData* finalimgColor;
// vtkDataArray* img_data;
// vtkDataArray* img_dataColor;
 vtkImageData* img;
	ofstream info;
	int imgheight, imgwidth;
	int npixels;
	int i, j,idx, intensity,counter;
	int new_idx;
	ifstream inFile;
	ofstream outFile;
	int max_thresh, min_thresh;
	int max_threshred, min_threshred;
	double value;
	double int_p, int_q;
	int height, width;
	int row, column;
	bool error_p;
	//double sigma;
	bool * redArray;
	int redArrIdx;
	bool redFlag;
	int * inImgData;
	double *pixelValue;
	char * file;
	clock_t t1,t2,t3;

	
  segmentation(char * filename)
  {
	max_thresh = min_thresh = 0;
	max_threshred = min_threshred = 0;
	value = 0;
	file = filename;
	redArrIdx = 0;
    
  }

  int segment()
  {
	outFile.open("C:\\Users\\user\\Desktop\\MKyan\\Marie\\GC2D\\bin\\Results\\final.txt");
	info.open("C:\\Users\\user\\Desktop\\MKyan\\Marie\\GC2D\\bin\\Results\\info.txt");

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

  
    img = extractSlice->GetOutput();
	imgheight = 256;
	imgwidth = 256;
	height = imgheight;
	width = imgwidth;
	npixels = imgheight* imgwidth;	
	counter = 0;
	inImgData = new int[npixels];

	//int* dim = img->GetDimensions();
 for (int z = 50; z < 51; z++)
    {for (int y = 0; y < 256; y++)
      {for (int x = 0; x < 256; x++) 
        {
			i=x;
			j=255-y;
			unsigned char* pixel = static_cast<unsigned char*>(img->GetScalarPointer(i,j,z));
        inImgData[counter] = int (pixel[0]);
           counter++;}
  //std::cout << inImgData[counter] << " ";
       }
    }

 t1=clock(); 

	info << "Number of Pixels: " << npixels << "\n";

	BKCutGraph* graph = new BKCutGraph(npixels);

	//inFile.open("C:\\Users\\user\\Desktop\\MKyan\\Test\\bin\\Results\\bseeds_idx.txt");
	inFile.open("C:\\Users\\user\\Desktop\\MKyan\\Marie\\GC2D\\bin\\bseeds_idx.txt");
	 if(!inFile) {
    cout << "Cannot open file.\n";
    //return 1;
  }
	 else{
		 cout << "Retrieving blue seeds.\n";
    
  }

	while (inFile >> idx) 
	{
		graph->add_weight(idx,graph->SOURCE,K);  
		//graph->add_weight(idx,graph->SINK,-K); 
		//for the given pixel, try to find intensity, to set up some thresholds
		min_thresh = 255;
		//i = idx / width ;//row
		//j = idx % width;//column
		//img_data = cvGet2D(img,i,j); // get the (i,j) pixel value
		//intensity = inImgData[idx];
		//if (max_thresh < intensity)
			//max_thresh = intensity;
		//if (min_thresh > intensity)
			//min_thresh = intensity;
		/*now min_thresh contains the lowest intensity value of the pixel.
		max_thresh contains the highest intensity value of the pixel*/ 
		
		//info << "Index:" << idx << " => Row:" << i << ", Column:" << j << "\n";
		

		/*outFile << "We are analyzing the following Image Pixel right now:";		
		outFile << "Index:" << idx << " => Row:" << i << ", Column:" << j << ", Intensity:" << intensity ;
		outFile << "\n";*/
    }
	info << "min:" << min_thresh << " max:" << max_thresh << "\n";
	inFile.close();

	////info.close();
	//info.open("C:\\results\\info.txt");

	//inFile.open("C:\\Users\\user\\Desktop\\MKyan\\Test\\bin\\Results\\rseeds_idx.txt");
	inFile.open("C:\\Users\\user\\Desktop\\MKyan\\Marie\\GC2D\\bin\\rseeds_idx.txt");
	 if(!inFile) {
    cout << "Cannot open file.\n";
  }
	  else{
		 cout << "Retrieving red seeds.\n";
    
  }
	redArray = new bool[npixels];
	for(counter=0;counter<npixels;counter++)
		redArray[counter] = false;
	idx = 0;
	while (inFile >> idx) 
	{		
		redArray[idx] = true;
		//i = idx / width ;//row
		//j = idx % width;//column
		//img_data = cvGet2D(img,i,j); // get the (i,j) pixel value
		//intensity = inImgData[idx];
		min_threshred = 255;
		//graph->add_weight(idx,graph->SOURCE,-K); 
		graph->add_weight(idx,graph->SINK,K); 
		if (max_threshred < intensity)
			max_threshred = intensity;
		if (min_threshred > intensity)
			min_threshred = intensity;
    }       
	info << "RED min:" << min_threshred << " max:" << max_threshred << "\n";
    inFile.close();
	redArrIdx = 0;

	i = 0;
	redFlag = false;

	while(i<npixels)
	{
		if(((int)((i/width)%2)) == 1)
		{
			//we are looking for only odd numbered rows
			for(j=0;j<width;j++)
			{
				//loop through all columns of this row
				//UP LINK
				if ((i-width) >= 0) //check if this is the first row (not usually possible)
				{
					//row = i/width;
					//idx = (row*width)+j;
					//row indicates the current row and j indicates the current column
					//img_data = cvGet2D(img,row,j); // get the (i,j) pixel value
					//int_p = img_data.val[0];
					int_p= inImgData[i];
					//row = (i-width)/width;
					//idx = (row*width)+j;
					//img_data = cvGet2D(img,row,j); // get the (i,j) pixel value
					int_q = inImgData[i-width];
					//info << "Row:" << row << " Column: " << j ;
					//info << " intensity1:" << int_p;
					//info << " intensity2: " << int_q;
					/*code Naaz 29th Jan 2009*/
					if(!redArray[i])
					{
						
						if(int_p <= max_thresh && int_p >= min_thresh)
						{
							info << "taken0";
							if(int_p <= max_threshred && int_p >= min_threshred)
							{
								info << "taken1";
								
								//graph->add_weight(i,graph->SOURCE,1);
								//graph->add_weight(i,graph->SINK,1);
							}
							else
							{
								info << "taken2";
								//good candidate for foreground
								graph->add_weight(i,graph->SOURCE,Ki);
								//graph->add_weight(i,graph->SINK,-Ki);
							}
						}
						else
						{
							info << "taken3";
							if(int_p <= max_threshred && int_p >= min_threshred)
							{
								info << "taken4";
								//graph->add_weight(i,graph->SOURCE,-Ki); 
								graph->add_weight(i,graph->SINK,Ki);
							}
						}
					}
					else
					{
						info << "taken5";
						graph->add_weight(i,graph->SOURCE,-K); 
						graph->add_weight(i,graph->SINK,K);
					}
					if(!redArray[i-width])
					{
						if(int_q <= max_thresh && int_q >= min_thresh)
						{
							//good candidate for foreground
							graph->add_weight(i-width,graph->SOURCE,Ki);
							//graph->add_weight(i-width,graph->SINK,-Ki);
						}
						else
						{
							//graph->add_weight(i-width,graph->SOURCE,-Ki); 
							graph->add_weight(i-width,graph->SINK,Ki);
						}
					}
					else
					{
						//graph->add_weight(i-width,graph->SOURCE,-K); 
						graph->add_weight(i-width,graph->SINK,K);
					}
				
					
					value = calculate_weight(int_p, int_q);
					graph->add_edge(i, (i-width),value, value);	
					
				}
				
				if (i < (width*(height-1))) //down link
				{
					//row = i / width;
					//idx = (row*width)+j;
					//int_p= inImgData[idx];
					int_p= inImgData[i];
					int_q = inImgData[i+width];
					//row = i / width;
					//img_data = cvGet2D(img,row,j); // get the (i,j) pixel value
					//int_p = img_data.val[0];
					//row = (i+width) / width;
					//img_data = cvGet2D(img,row,j); // get the (i,j) pixel value
					//int_q = img_data.val[0];
					if(!redArray[i])
					{
						if(int_p <= max_thresh && int_p >= min_thresh)
						{
							if(int_p <= max_threshred && int_p >= min_threshred)
							{						
								//good candidate for foreground
								graph->add_weight(i,graph->SOURCE,1);
								graph->add_weight(i,graph->SINK,1);
							}
							else
							{
								graph->add_weight(i,graph->SOURCE,Ki);
								//graph->add_weight(i,graph->SINK,-Ki);
							}
						}
						else
						{
							if(int_p <= max_threshred && int_p >= min_threshred)
							{
								//graph->add_weight(i,graph->SOURCE,-Ki); 
								graph->add_weight(i,graph->SINK,Ki);
							}
						}
					}
					else
					{
						//graph->add_weight(i,graph->SOURCE,-K); 
						graph->add_weight(i,graph->SINK,K);
					}
					if(!redArray[i+width])
					{
						if(int_q <= max_thresh && int_q >= min_thresh)
						{
							//good candidate for foreground
							graph->add_weight(i+width,graph->SOURCE,Ki);
							//graph->add_weight(i+width,graph->SINK,-Ki);
						}
						else
						{
							//graph->add_weight(i+width,graph->SOURCE,-Ki); 
							graph->add_weight(i+width,graph->SINK,Ki);
						}
					}
					else
					{
						//graph->add_weight(i+width,graph->SOURCE,-K); 
						graph->add_weight(i+width,graph->SINK,K);
					}
					//value = -((int_p - int_q)*(int_p - int_q))/sigma2;
					//int_p = exp(value)* sigma;
					//int_q = exp(value)* sigma;
					//graph->add_edge(i, (i+width),int_p, int_q);
					value = calculate_weight(int_p, int_q);
					graph->add_edge(i, (i+width),value, value);				
				}
				i++;
			}
		}
		else
		{
			i = i+width;
		}
	}

	for(i=0; i<npixels; i++)
	{
		//we are constructing only right and left links, alternate nodes only
		
		if(((int)i%2) == 1)
		{		
			if (((int)(i%width)) != 0) //left link
			{
				/*row = i / width;
				column = i % width;
				img_data = cvGet2D(img,row,column); // get the (i,j) pixel value
				int_p = img_data.val[0];
				row = (i-1) / width;
				img_data = cvGet2D(img,row,column); // get the (i,j) pixel value
				int_q = img_data.val[0];*/
				int_p= inImgData[i];
				int_q = inImgData[i-1];
				if(!redArray[i])
				{
					if(int_p <= max_thresh && int_p >= min_thresh)
					{
						if(int_p <= max_threshred && int_p >= min_threshred)
						{						
							//good candidate for foreground
							graph->add_weight(i,graph->SOURCE,1);
							graph->add_weight(i,graph->SINK,1);
						}
						else
						{
							graph->add_weight(i,graph->SOURCE,Ki);
							//graph->add_weight(i,graph->SINK,-Ki);
						}
					}
					else
					{
						if(int_p <= max_threshred && int_p >= min_threshred)
						{
							//graph->add_weight(i,graph->SOURCE,-Ki); 
							graph->add_weight(i,graph->SINK,Ki);
						}
					}
				}
				else
				{
					//graph->add_weight(i,graph->SOURCE,-K); 
					graph->add_weight(i,graph->SINK,K);
				}
				if(!redArray[i-1])
				{
					if(int_q <= max_thresh && int_q >= min_thresh)
					{
						//good candidate for foreground
						graph->add_weight(i-1,graph->SOURCE,Ki);
						//graph->add_weight(i-1,graph->SINK,-Ki);
					}
					else
					{
						//graph->add_weight(i-1,graph->SOURCE,-Ki); 
						graph->add_weight(i-1,graph->SINK,Ki);
					}
				}
				else
				{
					//graph->add_weight(i-1,graph->SOURCE,-K); 
					graph->add_weight(i-1,graph->SINK,K);
				}
				
				value = calculate_weight(int_p, int_q);
				graph->add_edge(i, (i-1),value, value);
			}
			if (((i+1)%width) != 0) //right link
			{
				/*row = i / width;
				column = i % width;
				img_data = cvGet2D(img,row,column); // get the (i,j) pixel value
				int_p = img_data.val[0];
				row = (i+1) / width;
				img_data = cvGet2D(img,row,column); // get the (i,j) pixel value
				int_q = img_data.val[0];*/
				int_p= inImgData[i];
				int_q = inImgData[i+1];
				if(!redArray[i])
				{
					if(int_p <= max_thresh && int_p >= min_thresh)
					{
						if(int_p <= max_threshred && int_p >= min_threshred)
						{						
							//good candidate for foreground
							graph->add_weight(i,graph->SOURCE,1);
							graph->add_weight(i,graph->SINK,1);
						}
						else
						{
							graph->add_weight(i,graph->SOURCE,Ki);
							//graph->add_weight(i,graph->SINK,-Ki);
						}
					}
					else
					{
						if(int_p <= max_threshred && int_p >= min_threshred)
						{
							//graph->add_weight(i,graph->SOURCE,-Ki); 
							graph->add_weight(i,graph->SINK,Ki);
						}
					}
				}
				else
				{
					//graph->add_weight(i,graph->SOURCE,-K); 
					graph->add_weight(i,graph->SINK,K);
				}
				if(!redArray[i+1])
				{
					if(int_q <= max_thresh && int_q >= min_thresh)
					{
						//good candidate for foreground
						graph->add_weight(i+1,graph->SOURCE,Ki);
						//graph->add_weight(i+1,graph->SINK,-Ki);
					}
					else
					{
						//graph->add_weight(i+1,graph->SOURCE,-Ki); 
						graph->add_weight(i+1,graph->SINK,Ki);
					}
				}
				else
				{
					//graph->add_weight(i+1,graph->SOURCE,-K); 
					graph->add_weight(i+1,graph->SINK,K);
				}
					
					value = calculate_weight(int_p, int_q);
					graph->add_edge(i, (i+1),value, value);
			}
		}
	}
	delete [] redArray;

  t2=clock();    

	graph->maxflow();

	t3=clock();

	float diff ((float)t2-(float)t1);  //time taken for graphcut prep
	float diff1 ((float)t3-(float)t2);  //time taken for graphcut
	float seconds = diff / CLOCKS_PER_SEC;
	 cout<<seconds<<endl;
	float seconds1 = diff1 / CLOCKS_PER_SEC;
    cout<<seconds1<<endl;


	error_p = graph->error();

// int* dimn = img->GetDimensions();

 

	/*vtkSmartPointer<vtkImageData> finalimg = 
    vtkSmartPointer<vtkImageData>::New();
  finalimg->SetDimensions(256,256,1);
  finalimg->SetSpacing(1.0,1.0,1.0);
  finalimg->SetOrigin(0.0, 0.0, 0.0);
  finalimg->AllocateScalars(VTK_DOUBLE ,1);*/
  
 vtkSmartPointer<vtkImageData> finalimgColor= 
    vtkSmartPointer<vtkImageData>::New();
  finalimgColor->SetDimensions(256,256,1);
  finalimgColor->SetSpacing(1.0,1.0,1.0);
  finalimgColor->SetOrigin(0.0, 0.0, 0.0);
  finalimgColor->AllocateScalars(VTK_DOUBLE ,3);

  vtkSmartPointer<vtkDoubleArray> img_data = 
    vtkSmartPointer<vtkDoubleArray>::New();
 img_data -> SetNumberOfComponents(1);

 vtkSmartPointer<vtkDoubleArray> img_dataColor = 
    vtkSmartPointer<vtkDoubleArray>::New();
 img_dataColor -> SetNumberOfComponents(3);

  /*vtkDataArray* img_data;
   img_data -> SetNumberOfComponents(1);

  vtkDataArray* img_dataColor;
  img_dataColor -> SetNumberOfComponents(3);*/

 pixelValue = new double[npixels];  //for original image
 
	for(idx=0;idx<npixels;idx++)
	{
		outFile << graph->what_segment(idx) << "\n";

	 if (graph->what_segment(idx) == 1)//SOURCE
        {		
			//img_data ->InsertNextTuple1(inImgData[idx]);
			img_dataColor ->InsertNextTuple3(inImgData[idx],0.0,0.0);  //blue for foreground
			//pixelValue[idx]= inImgData[idx];
	       }
     else
        {   
			//img_data ->InsertNextTuple1(0.0);
			img_dataColor ->InsertNextTuple3(0.0,0.0,inImgData[idx]);  //red for background
			//pixelValue[idx] =0.0;
			
          }
   // std::cout << pixelValue[idx] << " ";
 }   
       
	counter = 0;
	for (int z = 0; z < 1; z++)
    {for (int y = 0; y < 256; y++)
      {for (int x = 0; x < 256; x++)  
	{//double* pixelA = static_cast<double*>(finalimg->GetScalarPointer(x,y,z));
        // pixel[0]= inImgData[count]; count++;
      // pixelA[0]= pixelValue[counter]; 

		j = 255 - y;

	   double* pixelB = static_cast<double*>(finalimgColor->GetScalarPointer(x,j,z));
	   double *Tups = img_dataColor -> GetTuple3(counter);
	   pixelB[0]= Tups[0];
	   pixelB[1]= Tups[1];
	   pixelB[2]= Tups[2];

	   counter++; 
	}}}


	/*//Test 
	for (int z = 1; z < 2; z++)
    {for (int y = 50; y < 55; y++)
      {for (int x = 60; x < 65; x++)
	{double* pixel = static_cast<double*>(finalimg->GetScalarPointer(x,y,z));
	std::cout << pixel[0] << " ";
	double* pixel1 = static_cast<double*>(finalimgColor->GetScalarPointer(x,y,z));
	std::cout << pixel1[0] << " ";std::cout << pixel1[1] << " ";std::cout << pixel1[2] << " ";
	}}}*/

	delete [] inImgData;  
		 
		/* vtkSmartPointer<vtkJPEGWriter> writer =
			vtkSmartPointer<vtkJPEGWriter>::New(); 
		writer->SetFileName( "Final.jpg" );  
		writer->SetInputData( finalimg); 
		writer->Write(); */

	vtkSmartPointer<vtkImageProperty> Property =
    vtkSmartPointer<vtkImageProperty>::New();
	Property->SetAmbient(0.4);
	Property->SetDiffuse(0.6);
	Property->SetOpacity(1.0);

	 // Create actors
  vtkSmartPointer<vtkImageActor> inputActor1 = 
	  vtkSmartPointer<vtkImageActor>::New();
  inputActor1->GetMapper()->SetInputData(img);
  inputActor1->SetProperty(Property);

  vtkSmartPointer<vtkImageActor> inputActor2 = 
	  vtkSmartPointer<vtkImageActor>::New();
  inputActor2->GetMapper()->SetInputData(finalimgColor);
  inputActor2->SetProperty(Property);
 
  // There will be one render window
  vtkSmartPointer<vtkRenderWindow> renderWindow =
    vtkSmartPointer<vtkRenderWindow>::New();
  renderWindow->SetSize(700, 700);
 
  // And one interactor
  vtkSmartPointer<vtkRenderWindowInteractor> interactor =
    vtkSmartPointer<vtkRenderWindowInteractor>::New();
  interactor->SetRenderWindow(renderWindow);

   // Define viewport ranges
  // (xmin, ymin, xmax, ymax)
  double leftViewport[4] = {0.0, 0.0, 0.5, 1.0};
  double rightViewport[4] = {0.5, 0.0, 1.0, 1.0};
 
  // Setup both renderers
  vtkSmartPointer<vtkRenderer> leftRenderer = 
    vtkSmartPointer<vtkRenderer>::New();
  renderWindow->AddRenderer(leftRenderer);
  leftRenderer->SetViewport(leftViewport);
  leftRenderer->SetBackground(.6, .5, .4);  
  leftRenderer->AddActor(inputActor1);

  vtkSmartPointer<vtkRenderer> rightRenderer = 
    vtkSmartPointer<vtkRenderer>::New();
  renderWindow->AddRenderer(rightRenderer);
  rightRenderer->SetViewport(rightViewport);
  rightRenderer->SetBackground(.4, .5, .6);  
  rightRenderer->AddActor(inputActor2);
 
  leftRenderer->ResetCamera();
 
  rightRenderer->ResetCamera();
  rightRenderer->GetActiveCamera()->Azimuth(30);
  rightRenderer->GetActiveCamera()->Elevation(30);
 
  renderWindow->Render();
  interactor->Start();
 
  return EXIT_SUCCESS;
 
	}

   double calculate_weight(int p_int, int q_int)
{
    double power, final_value, temp;    
    temp = p_int - q_int;
    power = -(temp*temp/const1);
    final_value = exp(power)* const2;
    return (final_value);
}   
};