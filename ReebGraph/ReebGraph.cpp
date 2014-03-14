//
//  ReebGraph.cpp
//  ReebGraph
//
//  Created by Stephen John Russell on 11/03/2014.
//  Copyright (c) 2014 Stephen John Russell. All rights reserved.
//

#include "ReebGraph.h"


void ReebGraph::computeReebGraph() {
    
    int errorCode;
    
    cout << endl
    << "Reeb Graph Tests ========================== Surface Mesh Tests" << endl;
    
    // Loading the mesh
    vtkPolyData *surfaceMesh = vtkPolyData::New();
    LoadSurfaceMesh(surfaceMesh);
    
    // Attaching a height scalar field to it
    vtkDoubleArray *surfaceScalarField = vtkDoubleArray::New();
    surfaceScalarField->SetNumberOfTuples(surfaceMesh->GetNumberOfPoints());
    for(vtkIdType vId = 0; vId < surfaceMesh->GetNumberOfPoints(); vId++)
    {
        double *p = (double *) malloc(sizeof(double)*3);
        surfaceMesh->GetPoint(vId, p);
        double scalarValue = p[1];
        // add a bit of noise for the split tree test
        if(vId == 2) scalarValue -= 10*scalarValue;
        
        surfaceScalarField->SetTuple1(vId, scalarValue);
        free(p);
    }
    surfaceMesh->GetPointData()->SetScalars(surfaceScalarField);
    
    cout << "   Test 2D.1 Reeb graph computation... " << endl;
    vtkPolyDataToReebGraphFilter *surfaceReebGraphFilter =
    vtkPolyDataToReebGraphFilter::New();
    surfaceReebGraphFilter->SetInputData(surfaceMesh);
    surfaceReebGraphFilter->Update();
    vtkReebGraph *surfaceReebGraph = surfaceReebGraphFilter->GetOutput();
    cout << "      Test 2D.1 ";
    
    if(surfaceReebGraph->GetNumberOfEdges() == 12)
        cout << "OK!" << endl;
    else
    {
        cout << "Failed!" << endl;
        return EXIT_FAILURE;
    }
    
    cout << "   Test 2D.2 Customized Reeb graph simplification... " << endl;
    vtkReebGraphSimplificationFilter *surfaceSimplification =
    vtkReebGraphSimplificationFilter::New();
    
    AreaSimplificationMetric *metric = AreaSimplificationMetric::New();
    metric->SetLowerBound(0);
    // determining the maximum area
    double globalArea = 0;
    for(int i = 0; i < surfaceMesh->GetNumberOfCells(); i++)
    {
        vtkTriangle *t = vtkTriangle::SafeDownCast(surfaceMesh->GetCell(i));
        globalArea += t->ComputeArea();
    }
    metric->SetUpperBound(globalArea);
    surfaceSimplification->SetSimplificationMetric(metric);
    
    surfaceSimplification->SetInputData(surfaceReebGraph);
    surfaceSimplification->SetSimplificationThreshold(0.01);
    surfaceSimplification->Update();
    vtkReebGraph *simplifiedSurfaceReebGraph = surfaceSimplification->GetOutput();
    metric->Delete();
    
    cout << "      Test 2D.2 ";
    if(simplifiedSurfaceReebGraph->GetNumberOfEdges() == 12)
        cout << "OK!" << endl;
    else
    {
        cout << "Failed!" << endl;
        return EXIT_FAILURE;
    }
    
    cout << "   Test 2D.3 Reeb graph traversal..." << endl;
    errorCode = DisplayReebGraph(simplifiedSurfaceReebGraph);
    cout << "      Test 2D.3 ";
    if(!errorCode)
        cout << "OK!" << endl;
    else
    {
        cout << "Failed! (code " << errorCode << ")" << endl;
        return EXIT_FAILURE;
    }
    
    cout << "   Test 2D.4 Reeb graph based surface skeleton... " << endl;
    vtkReebGraphSurfaceSkeletonFilter *surfaceSkeletonFilter =
    vtkReebGraphSurfaceSkeletonFilter::New();
    surfaceSkeletonFilter->SetInputData(0, surfaceMesh);
    surfaceSkeletonFilter->SetInputConnection(1, surfaceSimplification->GetOutputPort());
    surfaceSkeletonFilter->SetNumberOfSamples(5);
    surfaceSkeletonFilter->Update();
    vtkTable *surfaceSkeleton = surfaceSkeletonFilter->GetOutput();
    errorCode = DisplaySurfaceSkeleton(surfaceMesh, surfaceSkeleton);
    cout << "      Test 2D.4 ";
    if(surfaceSkeleton->GetNumberOfColumns() == 12)
        cout << "OK!" << endl;
    else
    {
        cout << "Failed!" << endl;
        return EXIT_FAILURE;
    }
    
    cout << "   Test 2D.5 Area contour spectrum..." << endl;
    vtkAreaContourSpectrumFilter *areaSpectrumFilter =
    vtkAreaContourSpectrumFilter::New();
    areaSpectrumFilter->SetInputData(0, surfaceMesh);
    areaSpectrumFilter->SetInputConnection(1, surfaceSimplification->GetOutputPort());
    areaSpectrumFilter->SetArcId(0);
    areaSpectrumFilter->Update();
    vtkTable *areaSpectrum = areaSpectrumFilter->GetOutput();
    cout << "      Test 2D.5 ";
    if(areaSpectrum->GetNumberOfRows() == 100)
        cout << "OK!" << endl;
    else
    {
        cout << "Failed!" << endl;
        return EXIT_FAILURE;
    }
    
    
    cout << "   Test 2D.6 Reeb graph to split tree filter..." << endl;
    cout << "      Not currently tested..." << endl;
    //   vtkReebGraphToJoinSplitTreeFilter *splitTreeFilter =
    //     vtkReebGraphToJoinSplitTreeFilter::New();
    //   splitTreeFilter->SetInput(0, surfaceMesh);
    //   splitTreeFilter->SetInput(1, simplifiedSurfaceReebGraph);
    //   splitTreeFilter->SetIsSplitTree(true);
    //   splitTreeFilter->Update();
    //   vtkReebGraph *splitTree = splitTreeFilter->GetOutput();
    //   DisplayReebGraph(splitTree);
    //   cout << "      Test 2D.6 ";
    //   if(splitTree->GetNumberOfEdges() == 3)
    //     cout << "OK!" << endl;
    //   else
    //     {
    //     cout << "Failed!" << endl;
    //     return EXIT_FAILURE;
    //     }
    
    //   splitTreeFilter->Delete();
    areaSpectrumFilter->Delete();
    surfaceSkeletonFilter->Delete();
    surfaceSimplification->Delete();
    surfaceReebGraphFilter->Delete();
    surfaceScalarField->Delete();
    surfaceMesh->Delete();

    
    
}