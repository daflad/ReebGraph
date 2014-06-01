//
//  ReebComparator.h
//  ReebGraph
//
//  Created by Stephen John Russell on 25/05/2014.
//  Copyright (c) 2014 Stephen John Russell. All rights reserved.
//

#ifndef __ReebGraph__ReebComparator__
#define __ReebGraph__ReebComparator__

#include <iostream>
#include "vtkReebGraph.h"
#include "vtkPolyData.h"
#include "vtkTable.h"
#include "vtkActor.h"
#include "vtkAreaContourSpectrumFilter.h"
#include "vtkCamera.h"
#include "vtkDataSetAttributes.h"
#include "vtkDoubleArray.h"
#include "vtkEdgeListIterator.h"
#include "vtkIdList.h"
#include "vtkLight.h"
#include "vtkObjectFactory.h"
#include "vtkPointData.h"
#include "vtkPolyData.h"
#include "vtkPolyDataMapper.h"
#include "vtkPolyDataToReebGraphFilter.h"
#include "vtkProperty.h"
#include "vtkPNGWriter.h"
#include "vtkReebGraph.h"
#include "vtkReebGraphSurfaceSkeletonFilter.h"
#include "vtkReebGraphSimplificationFilter.h"
#include "vtkReebGraphSimplificationMetric.h"
#include "vtkReebGraphToJoinSplitTreeFilter.h"
#include "vtkReebGraphVolumeSkeletonFilter.h"
#include "vtkRenderer.h"
#include "vtkRenderWindow.h"
#include "vtkRenderWindowInteractor.h"
#include "vtkSmartPointer.h"
#include "vtkSphereSource.h"
#include "vtkTable.h"
#include "vtkTriangle.h"
#include "vtkUnstructuredGrid.h"
#include "vtkUnstructuredGridToReebGraphFilter.h"
#include "vtkVariantArray.h"
#include "vtkVolumeContourSpectrumFilter.h"
#include "vtkWindowToImageFilter.h"
#include "vtkDijkstraGraphGeodesicPath.h"

using namespace std;

// Compare two reeb graphs
class ReebComparator {
    
public:
    
    int minDeg;
    double minDist;
    
    void init(int, double);
    
    float compareReebGraphs(vtkReebGraph*, vtkReebGraph*, vtkPolyData*, vtkPolyData*);
    
private:
    
    
    
};

#endif /* defined(__ReebGraph__ReebComparator__) */
