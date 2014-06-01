//
//  ImgToMesh.h
//  ReebGraph
//
//  Created by Stephen John Russell on 18/03/2014.
//  Copyright (c) 2014 Stephen John Russell. All rights reserved.
//

#ifndef __ReebGraph__ImgToMesh__
#define __ReebGraph__ImgToMesh__

#include "vtkActor.h"
#include "vtkCellPicker.h"
#include "vtkCleanPolyData.h"
#include "vtkColor.h"
#include "vtkDataSetAttributes.h"
#include "vtkDataSetMapper.h"
#include "vtkDoubleArray.h"
#include "vtkEdgeListIterator.h"
#include "vtkExtractSelection.h"
#include "vtkIdList.h"
#include "vtkIdTypeArray.h"
#include "vtkImageMagnitude.h"
#include "vtkImageToPolyDataFilter.h"
#include "vtkImageQuantizeRGBToIndex.h"
#include "vtkInteractorStyleTrackballCamera.h"
#include "vtkLookupTable.h"
#include "vtkObjectFactory.h"
#include "vtkPNGReader.h"
#include "vtkPointData.h"
#include "vtkPolyData.h"
#include "vtkPolyDataConnectivityFilter.h"
#include "vtkPolyDataMapper.h"
#include "vtkPolyDataReader.h"
#include "vtkPolyDataWriter.h"
#include "vtkProperty.h"
#include "vtkRenderer.h"
#include "vtkRendererCollection.h"
#include "vtkRenderWindow.h"
#include "vtkRenderWindowInteractor.h"
#include "vtkSelection.h"
#include "vtkSelectionNode.h"
#include "vtkSmartPointer.h"
#include "vtkThresholdPoints.h"
#include "vtkTriangle.h"
#include "vtkTriangleFilter.h"
#include "vtkUnstructuredGrid.h"
#include "vtkUnstructuredGridToReebGraphFilter.h"
#include "vtkWindowToImageFilter.h"


using namespace std;

class ImgToMesh {
    
public:
    
    string fileName;
    string meshName;
    
    vtkSmartPointer<vtkPNGReader> reader;
    vtkSmartPointer<vtkImageQuantizeRGBToIndex> quantizer;
    vtkSmartPointer<vtkImageToPolyDataFilter> img2data;
    vtkSmartPointer<vtkTriangleFilter> triangleFilter;
    vtkSmartPointer<vtkPolyDataConnectivityFilter> connector;
    vtkSmartPointer<vtkCleanPolyData> cleaner;
    vtkSmartPointer<vtkPolyData> mesh;
    
    void init();
    
    bool loadFile();
    bool cullMesh();
    void displayMesh();
    
};

#endif /* defined(__ReebGraph__ImgToMesh__) */
