//
//  ReebGrapher.cpp
//  ReebGraph
//
//  Created by Stephen John Russell on 19/03/2014.
//  Copyright (c) 2014 Stephen John Russell. All rights reserved.
//

#include "ReebGrapher.h"

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

int ReebGrapher::DisplayReebGraph(vtkReebGraph *g) {
    
    vtkDataArray *vertexInfo = vtkDataArray::SafeDownCast(g->GetVertexData()->GetAbstractArray("Vertex Ids"));
    
    if(!vertexInfo) return 1;
    
    vtkVariantArray *edgeInfo = vtkVariantArray::SafeDownCast(
                                                              g->GetEdgeData()->GetAbstractArray("Vertex Ids"));
    if(!edgeInfo) return 2;
    
    cout << "   Reeb graph nodes:" << endl;
    for(int i = 0; i < vertexInfo->GetNumberOfTuples(); i++)
    cout << "      Node #" << i << ") VertexMeshId: "
    << ((int) *(vertexInfo->GetTuple(i))) << endl;
    
    cout << "   Reeb graph arcs:" << endl;
    vtkEdgeListIterator *eIt = vtkEdgeListIterator::New();
    g->GetEdges(eIt);
    do{
        vtkEdgeType e = eIt->Next();
        vtkAbstractArray *deg2NodeList = edgeInfo->GetPointer(e.Id)->ToArray();
        cout << "     Arc #" << e.Id << ": "
        << *(vertexInfo->GetTuple(e.Source)) << " -> "
        << *(vertexInfo->GetTuple(e.Target)) << " ("
        << deg2NodeList->GetNumberOfTuples() << " degree-2 nodes)" << endl;
    }while(eIt->HasNext());
    eIt->Delete();
    
    return 0;
}

int ReebGrapher::DisplaySurfaceSkeleton(vtkPolyData *surfaceMesh, vtkTable *skeleton)
{
    
    // Rendering setting
    vtkRenderer *renderer = vtkRenderer::New();
    
    vtkRenderWindow *renderWindow1 = vtkRenderWindow::New();
    renderWindow1->AddRenderer(renderer);
    renderWindow1->SetSize(800, 800);
    
    vtkRenderWindowInteractor *windowInteractor =
    vtkRenderWindowInteractor::New();
    windowInteractor->SetRenderWindow(renderWindow1);
    
    vtkPolyDataMapper *surfaceMapper = vtkPolyDataMapper::New();
    surfaceMapper->SetInputData(surfaceMesh);
    
    vtkActor *surfaceActor = vtkActor::New();
    surfaceActor->SetMapper(surfaceMapper);
    surfaceActor->GetProperty()->SetOpacity(0.4);
    
    vtkSmartPointer<vtkPolyDataMapper> mapper =
    vtkSmartPointer<vtkPolyDataMapper>::New();
    mapper->SetInputData(surfaceMesh);
    
    vtkSmartPointer<vtkActor> actor =
    vtkSmartPointer<vtkActor>::New();
    actor->SetMapper(mapper);
    actor->GetProperty()->SetRepresentationToWireframe();
    
    renderer->AddActor(actor);
    renderer->AddActor(surfaceActor);
    
    // Actual display of the skeleton
    vtkSphereSource *nodeSphere = vtkSphereSource::New();
    nodeSphere->SetThetaResolution(50);
    nodeSphere->SetPhiResolution(20);
    nodeSphere->SetRadius(0.03);
    
    vtkPolyDataMapper *sphereMapper = vtkPolyDataMapper::New();
    sphereMapper->SetInputConnection(nodeSphere->GetOutputPort());
    
    // 2 nodes per arc of the skeleton
    vtkActor **nodeActors = (vtkActor **) malloc(sizeof(vtkActor *)*
                                                 2*skeleton->GetNumberOfColumns());
    
    int     sampleId = 0;
    double  *point = (double *) malloc(sizeof(double)*3);
    
    vtkPolyData *embeddedSkeleton = vtkPolyData::New();
    embeddedSkeleton->Allocate();
    
    vtkPoints   *skeletonSamples = vtkPoints::New();
    skeletonSamples->SetNumberOfPoints(
                                       skeleton->GetNumberOfColumns()*skeleton->GetNumberOfRows());
    
    for(int i = 0; i < skeleton->GetNumberOfColumns(); i++)
    {
        vtkDoubleArray *arc = vtkDoubleArray::SafeDownCast(skeleton->GetColumn(i));
        
        // critical point at the origin of the arc
        arc->GetTupleValue(0, point);
        nodeActors[2*i] = vtkActor::New();
        nodeActors[2*i]->SetMapper(sphereMapper);
        nodeActors[2*i]->GetProperty()->SetColor(0, 0, 1);
        nodeActors[2*i]->SetPosition(point);
        renderer->AddActor(nodeActors[2*i]);
        
        arc->GetTupleValue(arc->GetNumberOfTuples() - 1, point);
        nodeActors[2*i + 1] = vtkActor::New();
        nodeActors[2*i + 1]->SetMapper(sphereMapper);
        nodeActors[2*i + 1]->GetProperty()->SetColor(0, 0, 1);
        nodeActors[2*i + 1]->SetPosition(point);
        renderer->AddActor(nodeActors[2*i + 1]);
        
        // now add the samples to the skeleton polyData
        int initialSampleId = sampleId;
        for(int j = 0; j < arc->GetNumberOfTuples(); j++)
        {
            arc->GetTupleValue(j, point);
            skeletonSamples->SetPoint(sampleId, point);
            sampleId++;
        }
        for(int j = 1; j < arc->GetNumberOfTuples(); j++)
        {
            vtkIdType samplePair[2];
            samplePair[0] = j - 1 + initialSampleId;
            samplePair[1] = j + initialSampleId;
            embeddedSkeleton->InsertNextCell(VTK_LINE, 2, samplePair);
        }
    }
    embeddedSkeleton->SetPoints(skeletonSamples);
    free(point);
    skeletonSamples->Delete();
    
    vtkPolyDataMapper *lineMapper = vtkPolyDataMapper::New();
    lineMapper->SetInputData(embeddedSkeleton);
    
    vtkActor          *skeletonActor = vtkActor::New();
    
    skeletonActor->SetMapper(lineMapper);
    skeletonActor->GetProperty()->SetColor(0, 1, 0);
    skeletonActor->GetProperty()->SetLineWidth(2);
    renderer->AddActor(skeletonActor);
    
    // Interactive mode
    windowInteractor->Initialize();
    windowInteractor->Start();
    
    skeletonActor->Delete();
    lineMapper->Delete();
    for(int i = 0; i < 2*skeleton->GetNumberOfColumns(); i++)
        nodeActors[i]->Delete();
    embeddedSkeleton->Delete();
    free(nodeActors);
    sphereMapper->Delete();
    nodeSphere->Delete();
    surfaceActor->Delete();
    surfaceMapper->Delete();
    windowInteractor->Delete();
    renderWindow1->Delete();
    renderer->Delete();
    
    cout << "HERE HERE HERE HERE HERE HERE HERE HERE HERE HERE HERE HERE " << endl;
    
    return 0;
}


