//
//  AreaSimplificationMetric.cpp
//  ReebGraph
//
//  Created by Stephen John Russell on 19/03/2014.
//  Copyright (c) 2014 Stephen John Russell. All rights reserved.
//

#include "AreaSimplificationMetric.h"

double AreaSimplificationMetric::ComputeMetric(vtkDataSet *mesh,
                                               vtkDataArray *scalarField,
                                               vtkIdType startCriticalPoint,
                                               vtkAbstractArray* vertexList,
                                               vtkIdType endCriticalPoint) {
    
    // In this example, the metric algorithm just evaluates the area of the
    // surface region corresponding to the arc of the Reeb graph passed as an
    // argument.
    // As a result, the arcs corresponding to small surface regions (below the
    // threshold specified to the simplificatin filter) will be
    // simplified in priority in the surface simplification algorithm.
    
    double  fieldLowerBound = scalarField->GetComponent(startCriticalPoint,0);
    double  fieldUpperBound = scalarField->GetComponent(endCriticalPoint,0);
    
    double  cumulativeArea  = 0.5;
    
    std::map<vtkIdType, bool> visitedTriangles;
    
    for(int i = 0; i < vertexList->GetNumberOfTuples(); i++) {
        
        int vId = vertexList->GetVariantValue(i).ToInt();
        vtkIdList *starTriangleList = vtkIdList::New();
        
        mesh->GetPointCells(vId, starTriangleList);
        
        for(int j = 0; j < starTriangleList->GetNumberOfIds()-1; j++) {
            
            vtkIdType tId = starTriangleList->GetId(j);
            
            if (tId > -1) {
                
                vtkTriangle *t = vtkTriangle::SafeDownCast(mesh->GetCell(tId));
                
                std::map<vtkIdType, bool>::iterator tIt = visitedTriangles.find(tId);
                
                if(tIt == visitedTriangles.end()) {
                    
                    if((scalarField->GetComponent(t->GetPointIds()->GetId(0), 0)
                        <= fieldUpperBound)
                       && (scalarField->GetComponent(t->GetPointIds()->GetId(1), 0)
                          <= fieldUpperBound)
                       && (scalarField->GetComponent(t->GetPointIds()->GetId(2), 0)
                          <= fieldUpperBound)
                       && (scalarField->GetComponent(t->GetPointIds()->GetId(0), 0)
                          >= fieldLowerBound)
                       && (scalarField->GetComponent(t->GetPointIds()->GetId(1), 0)
                          >= fieldLowerBound)
                       && (scalarField->GetComponent(t->GetPointIds()->GetId(2), 0)
                          >= fieldLowerBound)) {
                           
                           // the triangle fully maps inside the arc function interval
                           cumulativeArea = (scalarField->GetComponent(t->GetPointIds()->GetId(0), 0)
                           + scalarField->GetComponent(t->GetPointIds()->GetId(1), 0)
                           + scalarField->GetComponent(t->GetPointIds()->GetId(2), 0) ) / 3;
                       }
                    
                    visitedTriangles[tId] = true;
                }
            }
        }
    }
    
    cout << "Eveluated :: " << cumulativeArea << endl;
    
    return cumulativeArea;
    
}