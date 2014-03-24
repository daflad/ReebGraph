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
#include "vtkReebGraph.h"
#include "vtkPolyData.h"
#include "vtkTable.h"

class ReebGrapher {
    
public:
    
    int DisplayReebGraph(vtkReebGraph *g);
    
    int DisplaySurfaceSkeleton(vtkPolyData *surfaceMesh, vtkTable *skeleton);
    
    
    
};

#endif /* defined(__ReebGraph__ReebGrapher__) */
