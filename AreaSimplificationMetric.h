//
//  AreaSimplificationMetric.h
//  ReebGraph
//
//  Created by Stephen John Russell on 19/03/2014.
//  Copyright (c) 2014 Stephen John Russell. All rights reserved.
//

#ifndef __ReebGraph__AreaSimplificationMetric__
#define __ReebGraph__AreaSimplificationMetric__

#include <iostream>
#include <map>
#include "vtkReebGraphSimplificationMetric.h"
#include "vtkDataSet.h"
#include "vtkDataArray.h"
#include "vtkAbstractArray.h"
#include "vtkIdList.h"
#include "vtkTriangle.h"

class AreaSimplificationMetric : public vtkReebGraphSimplificationMetric{
public:
    vtkTypeMacro(AreaSimplificationMetric, vtkReebGraphSimplificationMetric);
    static AreaSimplificationMetric* New();
    double ComputeMetric(vtkDataSet *mesh, vtkDataArray *scalarField,
                         vtkIdType startCriticalPoint, vtkAbstractArray *vertexList,
                         vtkIdType endCriticalPoint);
};


#endif /* defined(__ReebGraph__AreaSimplificationMetric__) */
