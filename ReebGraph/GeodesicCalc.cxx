//
//  GeodesicCalc.cxx
//  ReebGraph
//
//  Created by Stephen John Russell on 18/03/2014.
//  Copyright (c) 2014 Stephen John Russell. All rights reserved.
//

#include "GeodesicCalc.h"

using namespace std;

// Calculate the distance from a single point to all other points along the mesh.
void GeodesicCalc::calculateDistances(vtkPolyData *surfaceMesh) {
    
    // Values to be attache dto the mesh
    vtkDoubleArray *scalarField = vtkDoubleArray::New();
    scalarField->SetNumberOfTuples(surfaceMesh->GetNumberOfPoints());
    
    // Algorithm to find shortst path between points
    vtkSmartPointer<vtkDijkstraGraphGeodesicPath> dijkstra =
    vtkSmartPointer<vtkDijkstraGraphGeodesicPath>::New();
    dijkstra->SetInputData(surfaceMesh);

    // Max X & Y locations
    double maxX = 0;
    double maxY = 0;
    vtkIdType index = 0;
    
    // Look for a unique ppint location.
    for(vtkIdType i = 0; i < surfaceMesh->GetNumberOfPoints(); i++) {
        double *p = (double *) malloc(sizeof(double)*3);
        surfaceMesh->GetPoint(i, p);
        if (p[0] > maxX && p[1] > maxY) {
            maxX = p[1];
            maxX = p[0];
            index = i;
        }
    }
    // list of distances to the point
    vector<double> dist;
    double maxDist = 0;
    // Loop to all other points
    for(vtkIdType i = 0; i < surfaceMesh->GetNumberOfPoints(); i++) {
        double d = 0;
        // Update Algorithm start and end points
        dijkstra->SetStartVertex(index);
        dijkstra->SetEndVertex(i);
        dijkstra->Update();
        
        //Recover list of points
        vtkIdList *l = dijkstra->GetIdList();
        d = l->GetNumberOfIds();
        for (int k = 1; k < l->GetNumberOfIds(); k++) {
            
            vtkIdType idA = l->GetId(k - 1);
            vtkIdType idB = l->GetId(k);
            
            double *pA = (double *) malloc(sizeof(double)*3);
            double *pB = (double *) malloc(sizeof(double)*3);
            
            surfaceMesh->GetPoint(idA, pA);
            surfaceMesh->GetPoint(idB, pB);
            
            // Find the squared distance between the points.
            double squaredDistance = vtkMath::Distance2BetweenPoints(pA, pB);
            
            // Take the square root to get the Euclidean distance between the points.
            d += sqrt(squaredDistance);
            free(pA);
            free(pB);
        }
        // Update Max
        if (d > maxDist) {
            maxDist = d;
        }
        // Add to list
        dist.push_back(d);
    }
    
    // normalisation
    for (int i = 0; i < dist.size(); i++) {
        scalarField->SetTuple1(i, dist[i] / maxDist);
    }
    //Add tp mesh
    surfaceMesh->GetPointData()->SetScalars(scalarField);
    
}


//Calc all geo distances.
void GeodesicCalc::calculateAverageDistances(vtkPolyData *surfaceMesh) {
    // list of average values
    vector<double> avg;
    // list of individual pont values
    vector<vector<double>> distLUT;
    long np = surfaceMesh->GetNumberOfPoints();
    // list of scal values
    vtkDoubleArray *scalarField = vtkDoubleArray::New();
    scalarField->SetNumberOfTuples(surfaceMesh->GetNumberOfPoints());
    // Shortest path algorithm
    vtkSmartPointer<vtkDijkstraGraphGeodesicPath> dijkstra =
    vtkSmartPointer<vtkDijkstraGraphGeodesicPath>::New();
    dijkstra->SetInputData(surfaceMesh);
    // max and min
    double aMax = 0;
    double aMin = 1000000;
    
    // Loop over each point to find average disance
    for(vtkIdType i = 0; i < np; i++) {
        
        // list of current distances for lookup
        vector<double> dist;
        double avgDist = 0;
        // Loop again to vidit all nodes once.
        for(vtkIdType j = 0; j < np; j++) {
            double d = 0;
            // Update location of start & end verts
            dijkstra->SetStartVertex(i);
            dijkstra->SetEndVertex(j);
            dijkstra->Update();
            vtkIdList *l = dijkstra->GetIdList();
            d = l->GetNumberOfIds();
            //loop length of path
            for (int k = 1; k < l->GetNumberOfIds(); k++) {
                vtkIdType idA = l->GetId(k - 1);
                // if distance is known, fetch and break out of loop.
                if (idA < distLUT.size()) {
                    if (j < distLUT[idA].size()) {
                        d += distLUT[idA][j ];
                        break;
                    }
                }
                // otherwise best calculate
                vtkIdType idB = l->GetId(k);
                double *pA = (double *) malloc(sizeof(double)*3);
                double *pB = (double *) malloc(sizeof(double)*3);
                surfaceMesh->GetPoint(idA, pA);
                surfaceMesh->GetPoint(idB, pB);
                // Find the squared distance between the points.
                double squaredDistance = vtkMath::Distance2BetweenPoints(pA, pB);

                // Take the square root to get the Euclidean distance between the points.
                d += sqrt(squaredDistance); //sqrt();
                free(pA);
                free(pB);
            }
            // Update avg dist & store
            avgDist += d;
            dist.push_back(d);
        }
        // Update list for node
        distLUT.push_back(dist);
        // Normalise
        avgDist /= (double)np;
        // Updae max & min
        if (avgDist > aMax) {
            aMax = avgDist;
        }
        if (avgDist < aMin) {
            aMin = avgDist;
        }
        // Store averages
        avg.push_back(avgDist);
        // Update user this is a long process!!!!!!
        cout << "Done :: " << i << "/" << np << " ~ " << avgDist << endl;
    }
    // Normalise & add to mesh
    for (int i = 0; i < np; i++) {
        scalarField->SetTuple1(i, (avg[i] - aMin)/(aMax - aMin));
    }
    surfaceMesh->GetPointData()->SetScalars(scalarField);
}