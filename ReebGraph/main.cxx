#define vtkRenderingCore_AUTOINIT 4(vtkInteractionStyle,vtkRenderingFreeType,vtkRenderingFreeTypeOpenGL,vtkRenderingOpenGL)
#define vtkRenderingVolume_AUTOINIT 1(vtkRenderingVolumeOpenGL)

#include <map>
#include <vector>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include "ImgToMesh.hpp"
#include "AreaSimplificationMetric.h"
#include "ReebGrapher.h"
#include "vtkDijkstraGraphGeodesicPath.h"

vtkStandardNewMacro(AreaSimplificationMetric);

double ma = 300;
double mi = 0;

double  loadDistancesFromFile(string fn, vtkDoubleArray *surfaceScalarField);
void    saveDistancesToFile(string fn, vtkDoubleArray *surfaceScalarField);

vector<double> avg;

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
    
    double aMax = loadDistancesFromFile("gun_dist.txt", surfaceScalarField);
    double aMin = 10000;
    
    long np = surfaceMesh->GetNumberOfPoints();
    
    if (aMax == -1) {
        
//        vector<vector<double>> distLUT;
//        vtkSmartPointer<vtkDijkstraGraphGeodesicPath> dijkstra =
//        vtkSmartPointer<vtkDijkstraGraphGeodesicPath>::New();
//        dijkstra->SetInputData(surfaceMesh);
//        
//        for(vtkIdType i = 0; i < np; i++) {
//
//            vector<double> dist;
//            double avgDist = 0;
//            for(vtkIdType j = 0; j < np; j++) {
//                double d = 0;
//                if (j < distLUT.size() && i < distLUT[j].size()) {
//                    d = distLUT[j][i];
//                } else {
//
//                    dijkstra->SetStartVertex(i);
//                    dijkstra->SetEndVertex(j);
//                    dijkstra->Update();
//                    vtkIdList *l = dijkstra->GetIdList();
//
//                    for (int k = 1; k < l->GetNumberOfIds(); k++) {
//                        vtkIdType idA = l->GetId(k - 1);
//                        vtkIdType idB = l->GetId(k);
//                        double *pA = (double *) malloc(sizeof(double)*3);
//                        double *pB = (double *) malloc(sizeof(double)*3);
//                        surfaceMesh->GetPoint(idA, pA);
//                        surfaceMesh->GetPoint(idB, pB);
//                        // Find the squared distance between the points.
//                        double squaredDistance = vtkMath::Distance2BetweenPoints(pA, pB);
//
//                        // Take the square root to get the Euclidean distance between the points.
//                        d += sqrt(squaredDistance);
//                        free(pA);
//                        free(pB);
//                    }
//
//                }
//                if (d > ma) {
//                    ma = d;
//                }
//                if (d < mi) {
//                    mi = d;
//                }
//                avgDist += d;
//                dist.push_back(d);
//            }
//            distLUT.push_back(dist);
//            avgDist /= (double)np;
//            if (avgDist > aMax) {
//                aMax = avgDist;
//            }
//            if (avgDist < aMin) {
//                aMin = avgDist;
//            }
//            avg.push_back(avgDist);
//            cout << "Done :: " << i << "/" << np << " ~ " << avgDist << endl;
//        }
//        for (int i = 0; i < np; i++) {
//            surfaceScalarField->SetTuple1(i, (avg[i] - aMin)/(aMax-aMin));
//        }
        double mMax = 0;
        vector<double> sv;
        vtkDoubleArray *surfaceScalarField = vtkDoubleArray::New();
        surfaceScalarField->SetNumberOfTuples(surfaceMesh->GetNumberOfPoints());
        for(vtkIdType vId = 0; vId < surfaceMesh->GetNumberOfPoints(); vId++) {
            double *p = (double *) malloc(sizeof(double)*3);
            surfaceMesh->GetPoint(vId, p);
            double scalarValue = p[1];
            
            if (scalarValue > mMax) {
                mMax = scalarValue;
            }
            
            sv.push_back(scalarValue);
            free(p);
        }
        for (int i = 0; i < sv.size(); i++) {
            surfaceScalarField->SetTuple1(i, sv[i]/mMax);
        }
        surfaceMesh->GetPointData()->SetScalars(surfaceScalarField);
        
        saveDistancesToFile("gun_dist.txt", surfaceScalarField);
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
    cout << "Max :: " << ma << " Min :: " << mi << endl;
    metric->SetLowerBound(mi);
    metric->SetUpperBound(ma);
    surfaceSimplification->SetSimplificationMetric(metric);
    
    surfaceSimplification->SetInputData(surfaceReebGraph);
    surfaceSimplification->SetSimplificationThreshold(0.005);
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


double loadDistancesFromFile(string fn, vtkDoubleArray *surfaceScalarField) {
    double max = 0;
    string line;
    ifstream myfile (fn);
    if (myfile.is_open())
    {
        int i = 0;
        int c = 0;
        while ( getline (myfile,line) )
        {
            if (c == 0) {
                ma = atof(line.c_str());
                
            }
            if (c == 1) {
                mi = atof(line.c_str());
                
            }
            if (c > 1) {
                double v = atof(line.c_str());
                surfaceScalarField->SetTuple1(i++, v);
                if (v > max) {
                    max = v;
                }
                cout << line << " :: " << v << '\n';
            }
            c++;
        }
        
        myfile.close();
        return max;
    } else {
        cout << "Unable to open file" << endl;
        return -1;
    }
    
}

void saveDistancesToFile(string fn, vtkDoubleArray *surfaceScalarField) {
    ofstream myfile;
    myfile.open (fn);
    myfile << ma << "\n" << mi << "\n";
    for (int i = 0; i < surfaceScalarField->GetNumberOfTuples(); i++) {
        myfile << surfaceScalarField->GetTuple(i)[0] << "\n";
    }
    myfile.close();
}

