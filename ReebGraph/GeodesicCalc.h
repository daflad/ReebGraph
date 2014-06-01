//
//  GeodesicCalc.h
//  ReebGraph
//
//  Created by Stephen John Russell on 31/03/2014.
//  Copyright (c) 2014 Stephen John Russell. All rights reserved.
//

#ifndef __ReebGraph__GeodesicCalc__
#define __ReebGraph__GeodesicCalc__

#include <iostream>
#include "vtkDoubleArray.h"
#include "vtkIdList.h"
#include "vtkPolyData.h"
#include "vtkSmartPointer.h"
#include "vtkDijkstraGraphGeodesicPath.h"
#include "vtkMath.h"
#include "vtkPointData.h"


//Helper class to get distances for scalar function of reeb graph
class GeodesicCalc {
    
public:
    
    void calculateDistances(vtkPolyData *surfaceMesh);
    void calculateAverageDistances(vtkPolyData *surfaceMesh);
    
private:
};

#endif /* defined(__ReebGraph__GeodesicCalc__) */
