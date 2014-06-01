//
//  ReebGrapher.cpp
//  ReebGraph
//
//  Created by Stephen John Russell on 19/03/2014.
//  Copyright (c) 2014 Stephen John Russell. All rights reserved.
//

#include "ReebGrapher.h"

// init class vars
void ReebGrapher::init(vtkPolyData *surfaceMesh) {
    mesh = surfaceMesh;
    ma = 0;
    mi = 100000;
}

// Build Reeb Graph
void ReebGrapher::buildReebGraph() {

    // Filter data
    vtkPolyDataToReebGraphFilter *surfaceReebGraphFilter =
    vtkPolyDataToReebGraphFilter::New();
    surfaceReebGraphFilter->SetInputData(mesh);
    surfaceReebGraphFilter->Update();
    surfaceReebGraph = surfaceReebGraphFilter->GetOutput();
   // Simplify
    vtkReebGraphSimplificationFilter *surfaceSimplification =
    vtkReebGraphSimplificationFilter::New();
    // Attach metric
    AreaSimplificationMetric *metric = AreaSimplificationMetric::New();
    metric->SetLowerBound(0);
    // determining the maximum area
    double globalArea = 0;
    for(int i = 0; i < mesh->GetNumberOfCells(); i++)
    {
        vtkTriangle *t = vtkTriangle::SafeDownCast(mesh->GetCell(i));
        double a = t->ComputeArea();
        globalArea += a;
    }
    metric->SetUpperBound(globalArea);
    surfaceSimplification->SetSimplificationMetric(metric);
    
    surfaceSimplification->SetInputData(surfaceReebGraph);

    // Set threshold for area simplification
    surfaceSimplification->SetSimplificationThreshold(0.001);
    surfaceSimplification->Update();
    simplifiedSurfaceReebGraph = surfaceSimplification->GetOutput();
    
    
    metric->Delete();
    surfaceReebGraphFilter->Delete();
    
    // Make skeleton
    vtkReebGraphSurfaceSkeletonFilter *surfaceSkeletonFilter =
    vtkReebGraphSurfaceSkeletonFilter::New();
    surfaceSkeletonFilter->SetInputData(0, mesh);
    surfaceSkeletonFilter->SetInputConnection(1, surfaceSimplification->GetOutputPort());
    surfaceSkeletonFilter->SetNumberOfSamples(7);
    surfaceSkeletonFilter->Update();
    skeleton = surfaceSkeletonFilter->GetOutput();
}

// Attach either scalar function values
void ReebGrapher::attachScalarField(int dimention) {
    // Let see how long it takes!!
    std::clock_t    start;
    start = std::clock();
    // chose x, y or z from mesh values
    if (dimention < 3) {
        // list fot normalisation
        std::vector<double> temp;
        vtkDoubleArray *scalarField = vtkDoubleArray::New();
        scalarField->SetNumberOfTuples(mesh->GetNumberOfPoints());
        // Retieve values
        for(vtkIdType vId = 0; vId < mesh->GetNumberOfPoints(); vId++)
        {
            double *p = (double *) malloc(sizeof(double)*3);
            mesh->GetPoint(vId, p);
            double scalarValue = p[dimention];
            
            if (scalarValue > ma) {
                ma = scalarValue;
            }
            if (scalarValue < mi) {
                mi = scalarValue;
            }
            temp.push_back(scalarValue);
            free(p);
        }
        // Normalise & add to mesh
        for (int i = 0 ; i < temp.size(); i++) {
            scalarField->SetTuple1(i, temp[i]/ma);
        }
        mesh->GetPointData()->SetScalars(scalarField);
    } else {
        // Go to GeodesicCalc and get distance
        GeodesicCalc gd;
        gd.calculateDistances(mesh);
    }
    // Stop time
    std::cout << "Time: " << (std::clock() - start) / (double)CLOCKS_PER_SEC << " s" << std::endl;
}


// Display textual reeb graph.
int ReebGrapher::DisplayReebGraph(vtkReebGraph *g) {
    
    // Vetrex or critical point info
    vtkDataArray *vertexInfo = vtkDataArray::SafeDownCast(g->GetVertexData()->GetAbstractArray("Vertex Ids"));
    
    if(!vertexInfo) {
        return 1;
    }
    
    // Edge or arc info
    
    vtkVariantArray *edgeInfo = vtkVariantArray::SafeDownCast(g->GetEdgeData()->GetAbstractArray("Vertex Ids"));
    if(!edgeInfo) {
        return 2;
    }
    
    // Whats there??
    // List vertexes
    cout << "   Reeb graph nodes:" << endl;
    for(int i = 0; i < vertexInfo->GetNumberOfTuples(); i++)
    cout << "      Node #" << i << ") VertexMeshId: "
    << ((int) *(vertexInfo->GetTuple(i))) << endl;
    // list arcs
    cout << "   Reeb graph arcs:" << endl;
    vtkEdgeListIterator *eIt = vtkEdgeListIterator::New();
    g->GetEdges(eIt);
    
    do {
        vtkEdgeType e = eIt->Next();
        vtkAbstractArray *deg2NodeList = edgeInfo->GetPointer(e.Id)->ToArray();
        cout << "     Arc #" << e.Id << ": "
        << *(vertexInfo->GetTuple(e.Source)) << " -> "
        << *(vertexInfo->GetTuple(e.Target)) << " ("
        << deg2NodeList->GetNumberOfTuples() << " degree-2 nodes)" << endl;
    } while(eIt->HasNext());
    
    eIt->Delete();
    
    return 0;
}

// Visual Reeb Graph display
int ReebGrapher::DisplaySurfaceSkeleton() {
    
    // Rendering setting
    vtkRenderer *renderer = vtkRenderer::New();
    // Rendering window
    vtkRenderWindow *renderWindow1 = vtkRenderWindow::New();
    renderWindow1->AddRenderer(renderer);
    renderWindow1->SetSize(800, 800);
    renderWindow1->SetWindowName("Reeb Graph");
    // Interact with data
    vtkRenderWindowInteractor *windowInteractor =
    vtkRenderWindowInteractor::New();
    windowInteractor->SetRenderWindow(renderWindow1);
    // Map data
    vtkPolyDataMapper *surfaceMapper = vtkPolyDataMapper::New();
    surfaceMapper->SetInputData(mesh);
    
    vtkActor *surfaceActor = vtkActor::New();
    surfaceActor->SetMapper(surfaceMapper);
    surfaceActor->GetProperty()->SetOpacity(0.4);
    
    vtkSmartPointer<vtkPolyDataMapper> mapper =
    vtkSmartPointer<vtkPolyDataMapper>::New();
    mapper->SetInputData(mesh);
    
    vtkSmartPointer<vtkActor> actor =
    vtkSmartPointer<vtkActor>::New();
    actor->SetMapper(mapper);
    actor->GetProperty()->SetRepresentationToWireframe();
    actor->GetProperty()->SetOpacity(0.04);
    
//    renderer->AddActor(actor);
    renderer->AddActor(surfaceActor);
    
    // Actual display of the skeleton
    vtkSphereSource *nodeSphere = vtkSphereSource::New();
    nodeSphere->SetThetaResolution(50);
    nodeSphere->SetPhiResolution(2);
    nodeSphere->SetRadius(.2);
    
    vtkPolyDataMapper *sphereMapper = vtkPolyDataMapper::New();
    sphereMapper->SetInputConnection(nodeSphere->GetOutputPort());
    
    // 2 nodes per arc of the skeleton
    vtkActor **nodeActors = (vtkActor **) malloc(sizeof(vtkActor *)*
                                                 2*skeleton->GetNumberOfColumns());
    
    int     sampleId = 0;
    double  *point = (double *) malloc(sizeof(double)*3);
    double  *pointE = (double *) malloc(sizeof(double)*3);
    
    vtkPolyData *embeddedSkeleton = vtkPolyData::New();
    embeddedSkeleton->Allocate();
    
    vtkPoints   *skeletonSamples = vtkPoints::New();
    skeletonSamples->SetNumberOfPoints(
                                       skeleton->GetNumberOfColumns()*skeleton->GetNumberOfRows());
    
    double c = 0;
    double dd = (double)1.0 / skeleton->GetNumberOfColumns();
    
    for(int i = 0; i < skeleton->GetNumberOfColumns(); i++) {
        vtkDoubleArray *arc = vtkDoubleArray::SafeDownCast(skeleton->GetColumn(i));
        
        // critical point at the origin of the arc
        
        arc->GetTupleValue(0, point);
        nodeActors[2*i] = vtkActor::New();
        nodeActors[2*i]->SetMapper(sphereMapper);
        nodeActors[2*i]->GetProperty()->SetColor(0, c, 1 - c);
        nodeActors[2*i]->SetPosition(point);
        renderer->AddActor(nodeActors[2*i]);
        
        arc->GetTupleValue(arc->GetNumberOfTuples() - 1, pointE);
        nodeActors[2*i + 1] = vtkActor::New();
        nodeActors[2*i + 1]->SetMapper(sphereMapper);
        nodeActors[2*i + 1]->GetProperty()->SetColor(0, c, 1 - c);
        nodeActors[2*i + 1]->SetPosition(point);
        renderer->AddActor(nodeActors[2*i + 1]);
        
        c += dd;
        
        // now add the samples to the skeleton polyData
        int initialSampleId = sampleId;
        
        for(int j = 0; j < arc->GetNumberOfTuples(); j++) {
            arc->GetTupleValue(j, point);
            skeletonSamples->SetPoint(sampleId, point);
            sampleId++;
        }
        for(int j = 1; j < arc->GetNumberOfTuples(); j++) {
            vtkIdType samplePair[2];
            samplePair[0] = j - 1 + initialSampleId;
            samplePair[1] = j + initialSampleId;
            embeddedSkeleton->InsertNextCell(VTK_LINE, 2, samplePair);
        }
    }
    // Update points
    embeddedSkeleton->SetPoints(skeletonSamples);
    free(point);
    skeletonSamples->Delete();
    
    vtkPolyDataMapper *lineMapper = vtkPolyDataMapper::New();
    lineMapper->SetInputData(embeddedSkeleton);
    
    vtkActor          *skeletonActor = vtkActor::New();
    
    skeletonActor->SetMapper(lineMapper);
    skeletonActor->GetProperty()->SetColor(0, 1, 0);
    skeletonActor->GetProperty()->SetLineWidth(5);
    renderer->AddActor(skeletonActor);
    
    // Interactive mode
    windowInteractor->Initialize();
    windowInteractor->Start();
    
    skeletonActor->Delete();
    lineMapper->Delete();
    
    // Remove refs
    for(int i = 0; i < 2*skeleton->GetNumberOfColumns(); i++) {
        nodeActors[i]->Delete();
    }
    
    embeddedSkeleton->Delete();
    free(nodeActors);
    sphereMapper->Delete();
    nodeSphere->Delete();
    surfaceActor->Delete();
    surfaceMapper->Delete();
    windowInteractor->Delete();
    renderWindow1->Delete();
    renderer->Delete();
//    surfaceSkeletonFilter->Delete();
    return 0;
}



