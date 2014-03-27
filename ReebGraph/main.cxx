#define vtkRenderingCore_AUTOINIT 4(vtkInteractionStyle,vtkRenderingFreeType,vtkRenderingFreeTypeOpenGL,vtkRenderingOpenGL)
#define vtkRenderingVolume_AUTOINIT 1(vtkRenderingVolumeOpenGL)

#include <map>
#include <vector>
#include "ImgToMesh.hpp"
#include "AreaSimplificationMetric.h"
#include "ReebGrapher.h"
#include "vtkDijkstraGraphGeodesicPath.h"

vtkStandardNewMacro(AreaSimplificationMetric);

int main(int vtkNotUsed(argc), char* vtkNotUsed(argv)[] ) {
    ImgToMesh im;
    ReebGrapher rg;
    
    im.init();
    
    int errorCode;
    
    cout << endl
    << "Reeb Graph Tests ========================== Surface Mesh Tests" << endl;
    
    // Loading the mesh
    vtkPolyData *surfaceMesh = vtkPolyData::New();
    if (im.loadFile()) {
        surfaceMesh->DeepCopy(im.mesh);
    } else {
        return EXIT_FAILURE;
    }
    
    // Attaching a height scalar field to it
    vtkDoubleArray *surfaceScalarField = vtkDoubleArray::New();
    surfaceScalarField->SetNumberOfTuples(surfaceMesh->GetNumberOfPoints());
    
    double max = 0;
    double min = 1000000;
    
    for(vtkIdType i = 0; i < surfaceMesh->GetNumberOfPoints(); i++) {
        cout << "Val :: " << surfaceScalarField->GetTuple(i)[0] << endl;
//        if (surfaceScalarField->GetTuple(i)[0] == 0) {
        
            vector<double> dist;
            
            for(vtkIdType j = 0; j < surfaceMesh->GetNumberOfPoints(); j++) {
                
                vtkSmartPointer<vtkDijkstraGraphGeodesicPath> dijkstra =
                vtkSmartPointer<vtkDijkstraGraphGeodesicPath>::New();
                dijkstra->SetInputData(surfaceMesh);
                dijkstra->SetStartVertex(i);
                dijkstra->SetEndVertex(j);
                dijkstra->Update();
                vtkIdList *l = dijkstra->GetIdList();
                double d = 0;
                for (int k = 1; k < l->GetNumberOfIds(); k++) {
                    l->GetId(k);
                    double *pA = (double *) malloc(sizeof(double)*3);
                    double *pB = (double *) malloc(sizeof(double)*3);
                    surfaceMesh->GetPoint(k - 1, pA);
                    surfaceMesh->GetPoint(k, pB);
                    // Find the squared distance between the points.
                    double squaredDistance = vtkMath::Distance2BetweenPoints(pA, pB);
                    
                    // Take the square root to get the Euclidean distance between the points.
                    d += sqrt(squaredDistance);
                    free(pA);
                    free(pB);
                }
                dist.push_back(d);
                
            }
            double avgDist = 0;
            for (int j = 0; j < dist.size(); j++) {
                avgDist += dist[j];
            }
            avgDist /= dist.size();
            surfaceScalarField->SetTuple1(i, avgDist);
            cout << "Done :: " << i << " Val :: " << surfaceScalarField->GetTuple(i)[0] << endl;;
//        }
    }
    surfaceMesh->GetPointData()->SetScalars(surfaceScalarField);
    
    vtkSmartPointer<vtkPolyDataWriter> writer =
    vtkSmartPointer<vtkPolyDataWriter>::New();
    writer->SetInputData(surfaceMesh);
    writer->SetFileName("gun.vtk");
    writer->Update();
    writer->Write();
    
    cout << "   Test 2D.1 Reeb graph computation... " << endl;
    vtkPolyDataToReebGraphFilter *surfaceReebGraphFilter =
    vtkPolyDataToReebGraphFilter::New();
    surfaceReebGraphFilter->SetInputData(surfaceMesh);
    surfaceReebGraphFilter->Update();
    vtkReebGraph *surfaceReebGraph = surfaceReebGraphFilter->GetOutput();
    cout << "      Test 2D.1 ";
    
    if(surfaceReebGraph->GetNumberOfEdges() > 0)
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
    metric->SetLowerBound(min);
    metric->SetUpperBound(max);
    surfaceSimplification->SetSimplificationMetric(metric);
    
    surfaceSimplification->SetInputData(surfaceReebGraph);
    //surfaceSimplification->SetSimplificationThreshold(0.4);
    surfaceSimplification->Update();
    vtkReebGraph *simplifiedSurfaceReebGraph = surfaceSimplification->GetOutput();
    metric->Delete();
    
    cout << "      Test 2D.2 ";
    if(simplifiedSurfaceReebGraph->GetNumberOfEdges() > 0)
        cout << "OK!" << endl;
    else
    {
        cout << "Failed!" << endl;
        return EXIT_FAILURE;
    }
    
    cout << "   Test 2D.3 Reeb graph traversal..." << endl;
    errorCode = rg.DisplayReebGraph(simplifiedSurfaceReebGraph);
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
    errorCode = rg.DisplaySurfaceSkeleton(surfaceMesh, surfaceSkeleton);
    cout << "      Test 2D.4 ";
    if(surfaceSkeleton->GetNumberOfColumns() > 0)
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
    if(areaSpectrum->GetNumberOfRows() > 0)
        cout << "OK!" << endl;
    else
    {
        cout << "Failed!" << endl;
        return EXIT_FAILURE;
    }
    
    areaSpectrumFilter->Delete();
    surfaceSkeletonFilter->Delete();
    surfaceSimplification->Delete();
    surfaceReebGraphFilter->Delete();
    surfaceMesh->Delete();
    
    return 0;
}

