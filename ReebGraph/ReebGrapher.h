//
//  ReebGrapher.h
//  ReebGraph
//
//  Created by Stephen John Russell on 19/03/2014.
//  Copyright (c) 2014 Stephen John Russell. All rights reserved.
//

#ifndef __ReebGraph__ReebGrapher__
#define __ReebGraph__ReebGrapher__

#include <iostream>
#include "GeodesicCalc.h"
#include "AreaSimplificationMetric.h"
#include <ctime>
#include "vtkPolyData.h"
#include "vtkTable.h"
#include "vtkActor.h"
#include "vtkDoubleArray.h"
#include "vtkEdgeListIterator.h"
#include "vtkIdList.h"
#include "vtkPointData.h"
#include "vtkPolyData.h"
#include "vtkPolyDataMapper.h"
#include "vtkPolyDataToReebGraphFilter.h"
#include "vtkProperty.h"
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
#include "vtkVariantArray.h"

// Build and display reeb graph
class ReebGrapher {
    
public:
    
    double ma;
    double mi;
    
    int DisplayReebGraph(vtkReebGraph *g);
    
    int DisplaySurfaceSkeleton();
    
    void init(vtkPolyData *surfaceMesh);
    
    void attachScalarField(int dimention);
    
    void buildReebGraph();
    
    vtkPolyData *mesh;
    vtkReebGraphSimplificationFilter *surfaceSimplification;
    vtkReebGraph *surfaceReebGraph;
    vtkReebGraph *simplifiedSurfaceReebGraph;
    vtkTable *skeleton;
};

#endif /* defined(__ReebGraph__ReebGrapher__) */
