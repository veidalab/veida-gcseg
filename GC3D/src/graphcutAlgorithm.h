#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <ctype.h>
#include <iostream>
#include <fstream>
#include<math.h>
#include <time.h>

using std::cout;
using std::ofstream;
using std::ifstream;

#include "graphcut.h"

#include <vtkImageData.h>
#include <vtkSmartPointer.h>
#include <vtkImageReader2.h>
#include <vtkExtractVOI.h>
#include <vtkDataArray.h>
#include <vtkDoubleArray.h>
#include <vtkUnsignedCharArray.h>
#include <vtkJPEGWriter.h>
#include <vtkInteractorStyleImage.h>
#include <vtkLookupTable.h>
#include <vtkVolumeMapper.h>
#include <vtkSmartVolumeMapper.h>
#include <vtkGPUVolumeRayCastMapper.h>
#include <vtkFixedPointVolumeRayCastMapper.h>
#include <vtkColorTransferFunction.h>
#include <vtkPiecewiseFunction.h>
#include <vtkImageMapToColors.h>
#include <vtkVolumeProperty.h>
#include <vtkRenderWindow.h>
#include <vtkRenderer.h>
#include <vtkCamera.h>
#include <vtkRenderWindowInteractor.h>

#include <vtkImageProperty.h>
#include <vtkImageActor.h>
#include <vtkImageMapper3D.h>

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
	vtkImageData* img0;
	vtkImageData* img1;
	vtkImageData* img2;
	vtkImageData* imgOrigin;
	ofstream info;
	int imgheight, imgwidth;
	int npixels;
	int i, j,idx,k,intensity,counter;
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
	bool * redArray;
	int redArrIdx;
	bool redFlag;
	char * file0;
	char * file;
	char * file2;
	int * inImgData0;
	int * inImgData1;
	int * inImgData2;
	int * mask0;
	int * mask1;
	int * mask2;
	int marked;

	clock_t t1,t2;
	
  segmentation(char * filename0, char * filename, char* filename2, int markedImage)
  {
	max_thresh = min_thresh = 0;
	max_threshred = min_threshred = 0;
	value = 0;
	file0 = filename0;
	file = filename;
	file2 = filename2;
	marked = markedImage;
	redArrIdx = 0;
    
  }

  int segment()
  {
	outFile.open("C:\\Users\\user\\Desktop\\MKyan\\Marie\\GC3D\\bin\\final.txt");
	info.open("C:\\Users\\user\\Desktop\\MKyan\\Marie\\GC3D\\bin\\info.txt");
	
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
 
cout << "Extracting slices..." << endl;
int dims=255;
  //extract a slice from the volume
  vtkSmartPointer<vtkExtractVOI> extractSlice0 =
      vtkSmartPointer<vtkExtractVOI>::New();
  extractSlice0->SetInputConnection(reader->GetOutputPort());
  extractSlice0->SetVOI(0,dims,0,dims,49,49);
  extractSlice0->Update();
  img0 = extractSlice0->GetOutput();

  vtkSmartPointer<vtkExtractVOI> extractSlice1 =
      vtkSmartPointer<vtkExtractVOI>::New();
  extractSlice1->SetInputConnection(reader->GetOutputPort());
  extractSlice1->SetVOI(0,dims,0,dims,50,50);
  extractSlice1->Update();
  img1 = extractSlice1->GetOutput();

  vtkSmartPointer<vtkExtractVOI> extractSlice2 =
      vtkSmartPointer<vtkExtractVOI>::New();
  extractSlice2->SetInputConnection(reader->GetOutputPort());
  extractSlice2->SetVOI(0,dims,0,dims,51,51);
  extractSlice2->Update();
  img2 = extractSlice2->GetOutput();

  vtkSmartPointer<vtkExtractVOI> extractSlice =
      vtkSmartPointer<vtkExtractVOI>::New();
  extractSlice->SetInputConnection(reader->GetOutputPort());
  extractSlice->SetVOI(0,dims,0,dims,49,51);
  extractSlice->Update();
  imgOrigin = extractSlice->GetOutput();

	imgheight = dims+1;
	imgwidth = dims+1;
	height = imgheight;
	width = imgwidth;
	npixels = imgheight* imgwidth;	
	k=0;

	inImgData0 = new int[npixels];	
	inImgData1 = new int[npixels];
	inImgData2 = new int[npixels];

	int counter=0;

    {for (int y = 0; y < 256; y++)
      {for (int x = 0; x < 256; x++) 
        { j=255-y;
			unsigned char* pixel0 = static_cast<unsigned char*>(img0->GetScalarPointer(x,j,49));
        inImgData0[counter] = int (pixel0[0]);
		unsigned char* pixel1 = static_cast<unsigned char*>(img1->GetScalarPointer(x,j,50));
        inImgData1[counter] = int (pixel1[0]);
		unsigned char* pixel2 = static_cast<unsigned char*>(img2->GetScalarPointer(x,j,51));
        inImgData2[counter] = int (pixel2[0]);
		counter++;}}}	
	
	info << "Number of Pixels: " << npixels << "\n";

	 t1=clock();

	BKCutGraph* graph = new BKCutGraph(npixels*3);	

	//inFile.open("C:\\Users\\user\\Desktop\\MKyan\\seg2D\\segmentation1\\results\\bseeds_idx.txt");
	inFile.open("C:\\Users\\user\\Desktop\\MKyan\\Marie\\GC3D\\bin\\bseeds_idx.txt");
	if(inFile.is_open())
	{cout << "Reading blue seeds \n";}
	else
	{cout << "Error Reading File... \n";}

	//Now the seeds are for img (i.e. the middle slice, so we need to add an offset
	while (inFile >> idx) 
	{
		graph->add_weight((idx+(npixels*marked)),graph->SOURCE,K);  
		//for the given pixel, try to find intensity, to set up some thresholds
		min_thresh = 255;
    }
	info << "min:" << min_thresh << " max:" << max_thresh << "\n";
	inFile.close();

	//inFile.open("C:\\Users\\user\\Desktop\\MKyan\\seg2D\\segmentation1\\results\\rseeds_idx.txt");
	inFile.open("C:\\Users\\user\\Desktop\\MKyan\\Marie\\GC3D\\bin\\rseeds_idx.txt");
	if(inFile.is_open())
	{cout << "Reading red seeds \n";}
	else
	{cout << "Error Reading File... \n";}

	redArray = new bool[npixels];
	for(counter=0;counter<npixels;counter++)
		redArray[counter] = false;
	idx = 0;
	while (inFile >> idx) 
	{		
		redArray[idx] = true;		
		min_threshred = 255;
		
		graph->add_weight((idx+(npixels*marked)),graph->SINK,K); 
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

	while(i<(npixels*3))
	{
		k = (int)i/npixels;
		if(((int)(((i-(k*npixels))/width)%2)) == 1)
		{
			k = (int)i/npixels;
			//we are looking for only odd numbered rows
			for(j=0;j<width;j++)
			{
				//loop through all columns of this row
				//UP LINK
				if ((i-(k*npixels)-width) >= 0) //check if this is the first row (not usually possible)
				{
					
					if(k == 0)
					{
						int_p= inImgData0[i];
						int_q = inImgData0[i-width];
					}
					else if(k==1) 
					{
						int_p= inImgData1[i-(npixels*k)];
						int_q = inImgData1[i-(npixels*k)-width];
					}
					else if(k==2) 
					{
						int_p= inImgData2[i-(npixels*k)];
						int_q = inImgData2[i-(npixels*k)-width];
					}
					
					if(!redArray[i-(npixels*k)])
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
					if(!redArray[i-(npixels*k)])
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
				
				if (i-(npixels*k) < (width*(height-1))) //down link
				{
					if(k == 0)
					{
						int_p= inImgData0[i];
						int_q = inImgData0[i+width];
					}
					else if(k==1) 
					{
						int_p= inImgData1[i-(npixels*k)];
						int_q = inImgData1[i-(npixels*k)+width];
					}
					else if(k==2) 
					{
						int_p= inImgData2[i-(npixels*k)];
						int_q = inImgData2[i-(npixels*k)+width];
					}
					if(!redArray[i-(npixels*k)])
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
					if(!redArray[i-(npixels*k)+ width])
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
						graph->add_weight(i+width,graph->SINK,K);
					}
					
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

	for(i=0; i<(npixels*3); i++)
	{
		if(((int)(i-(npixels*k))%2) == 1)
		{		
			if (((int)((i-(npixels*k))%width)) != 0) //left link
			{
				if(k == 0)
					{
						int_p= inImgData0[i];
						int_q = inImgData0[i-1];
					}
					else if(k==1) 
					{
						int_p= inImgData1[i-(npixels*k)];
						int_q = inImgData1[i-(npixels*k)-1];
					}
					else if(k==2) 
					{
						int_p= inImgData2[i-(npixels*k)];
						int_q = inImgData2[i-(npixels*k)-1];
					}
				if(!redArray[i-(npixels*k)])
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
				if(!redArray[i-(npixels*k)-1])
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
			if ((((i-(npixels*k))+1)%width) != 0) //right link
			{
				if(k == 0)
					{
						int_p= inImgData0[i];
						int_q = inImgData0[i+1];
					}
					else if(k==1) 
					{
						int_p= inImgData1[i-(npixels*k)];
						int_q = inImgData1[i-(npixels*k)+1];
					}
					else if(k==2) 
					{
						int_p= inImgData2[i-(npixels*k)];
						int_q = inImgData2[i-(npixels*k)+1];
					}
				if(!redArray[(i-(npixels*k))])
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
				if(!redArray[(i-(npixels*k)+1)])
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

	for(idx=0;idx<npixels;idx++)
	{
		int_p = inImgData0[idx];
		int_q = inImgData1[idx];
		value = calculate_weight(int_p, int_q);
		graph->add_edge(idx, (idx+npixels),value, value);
		int_p = inImgData1[idx];
		int_q = inImgData2[idx];
		value = calculate_weight(int_p, int_q);
		graph->add_edge(idx+npixels, (idx+(npixels*2)),value, value);
	}
	delete [] redArray;
    
	graph->maxflow();

	error_p = graph->error();

	t2=clock();
    float diff ((float)t2-(float)t1);
	float seconds = diff / CLOCKS_PER_SEC;
    cout<<seconds<<endl;
 
     mask0 = new int[npixels];	
	 mask1 = new int[npixels];
	 mask2 = new int[npixels];

	  for(idx=0;idx<npixels;idx++)
	{
		outFile << graph->what_segment(idx) << "\n";

	 if (graph->what_segment(idx) == 1)//SOURCE
        { mask0[idx] = 1;
	    }
     else
        {mask0[idx] = 0;
		}

	 if (graph->what_segment(idx+npixels) == 1)//SOURCE
	 { mask1[idx] = 1;
	 }
	 else
	 {mask1[idx] =0;
	 }

	 if (graph->what_segment(idx+(npixels*2)) == 1)//SOURCE
	 {mask2[idx] = 1;
	 }
	 else
	 {mask2[idx] =0;
	 }
 } 

	vtkSmartPointer<vtkImageData> masking= 
    vtkSmartPointer<vtkImageData>::New();
  masking->SetDimensions(256,256,3);
  masking->SetSpacing(1.0,1.0,1.0);
  masking->SetOrigin(0.0, 0.0, 0.0);
  masking->AllocateScalars(VTK_UNSIGNED_CHAR ,1); 

  
  counter = 0;

    {for (int y = 0; y < 256; y++)
      {for (int x = 0; x < 256; x++)  
	{j=255-y;
	   unsigned char* pixel0 = static_cast<unsigned char*>(masking->GetScalarPointer(x,j,0));
	   int Tups0 = mask0[counter];
	   pixel0[0]= Tups0;

	  unsigned char* pixel1 = static_cast<unsigned char*>(masking->GetScalarPointer(x,j,1));
	   int Tups1 = mask1[counter];
	   pixel1[0]= Tups1;

	   unsigned char* pixel2 = static_cast<unsigned char*>(masking->GetScalarPointer(x,j,2));
	   int Tups2 = mask2[counter];
	   pixel2[0]= Tups2;

	   counter++; 
	}}}

	delete [] inImgData0;
	delete [] inImgData1;
	delete [] inImgData2;
	delete []  mask0;
	delete []  mask1;
	delete []  mask2;

	outFile.close();
	info.close();

	//3D rendering

	vtkSmartPointer<vtkRenderer> ren =
		vtkSmartPointer<vtkRenderer>::New();
	vtkSmartPointer<vtkRenderWindow> renWin =
		vtkSmartPointer<vtkRenderWindow>::New();
	renWin->AddRenderer(ren);
	vtkSmartPointer<vtkRenderWindowInteractor> iren =
		vtkSmartPointer<vtkRenderWindowInteractor>::New();
	iren->SetRenderWindow(renWin);

	vtkSmartPointer<vtkGPUVolumeRayCastMapper> volumeMapper =
		vtkSmartPointer<vtkGPUVolumeRayCastMapper>::New();
	volumeMapper->SetBlendModeToComposite();
	volumeMapper->SetInputConnection(extractSlice->GetOutputPort());
	volumeMapper->SetMaskTypeToLabelMap(); //SetMaskTypeToLabelMap(); SetMaskTypeToBinary();
	volumeMapper->SetMaskInput(masking);
	volumeMapper->SetAutoAdjustSampleDistances(0);
	volumeMapper->SetMaskBlendFactor(1.0);  //1.0 for mask color transfer fn; 0.0 for standard color transfer fn
	

	//double distance = 0.3477 / 2.0;
    //volumeMapper->SetSampleDistance(static_cast<float>(distance));

	vtkSmartPointer<vtkVolume> volume =
		vtkSmartPointer<vtkVolume>::New();

	vtkSmartPointer<vtkColorTransferFunction> colorFun =
		vtkSmartPointer<vtkColorTransferFunction>::New();
	vtkSmartPointer<vtkPiecewiseFunction>opacityFun =
		vtkSmartPointer<vtkPiecewiseFunction>::New();

    colorFun->AddRGBPoint(0, 1.0, 0.0, 0.0);
	//colorFun->AddRGBPoint(1, 1.0, 0.0, 0.0);
	colorFun->AddRGBPoint(1, 0.0, 0.0, 1.0);
	//colorFun->AddRGBPoint(155.0, 1.0, 0.0, 0.0);

	opacityFun->AddPoint(0, 0.15);
	opacityFun->AddPoint(1, 0.65);
	//opacityFun->AddPoint(255, 0.95);
	

	vtkSmartPointer<vtkVolumeProperty> property =
		vtkSmartPointer<vtkVolumeProperty>::New();
	property->SetColor(colorFun);
	property->SetScalarOpacity(opacityFun);
	property->SetInterpolationTypeToLinear();

	volume->SetProperty(property);
	volume->SetMapper(volumeMapper);	

	/*property->ShadeOn();
	property->SetAmbient(0.45);
	property->SetDiffuse(0.55);
	property->SetSpecular(0.6);
	property->SetSpecularPower(8.5); */

	renWin->SetSize(600, 600);
	renWin->Render();

	// Add the volume to the scene
	ren->AddVolume(volume);
	ren->SetBackground(.1, .2, .3);
	ren->ResetCamera();

	// interact with data
	renWin->Render();

	iren->Start();

	opacityFun->Delete();
	colorFun->Delete();
	property->Delete();

	volume->Delete();
	volumeMapper->Delete();

	ren->Delete();
	renWin->Delete();
	iren->Delete();
	masking->Delete();
 
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
