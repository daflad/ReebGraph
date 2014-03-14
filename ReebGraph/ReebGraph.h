//
//  ReebGraph.h
//  ReebGraph
//
//  Created by Stephen John Russell on 11/03/2014.
//  Copyright (c) 2014 Stephen John Russell. All rights reserved.
//

#ifndef __ReebGraph__ReebGraph__
#define __ReebGraph__ReebGraph__

#include <iostream>
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
#include "vtkSphereSource.h"
#include "vtkTable.h"
#include "vtkTriangle.h"
#include "vtkUnstructuredGrid.h"
#include "vtkUnstructuredGridToReebGraphFilter.h"
#include "vtkVariantArray.h"
#include "vtkVolumeContourSpectrumFilter.h"
#include "vtkWindowToImageFilter.h"

#include <map>

class ReebGraph {
    
public:
    
    void computeReebGraph();
};

#endif /* defined(__ReebGraph__ReebGraph__) */
