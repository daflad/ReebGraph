#define vtkRenderingCore_AUTOINIT 4(vtkInteractionStyle,vtkRenderingFreeType,vtkRenderingFreeTypeOpenGL,vtkRenderingOpenGL)
#define vtkRenderingVolume_AUTOINIT 1(vtkRenderingVolumeOpenGL)

#include <map>
#include "ImgToMesh.hpp"
#include "AreaSimplificationMetric.h"
#include "ReebGrapher.h"

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
    for(vtkIdType vId = 0; vId < surfaceMesh->GetNumberOfPoints(); vId++) {
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
//    surfaceSimplification->SetSimplificationThreshold(1 - 0.0000001);
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

