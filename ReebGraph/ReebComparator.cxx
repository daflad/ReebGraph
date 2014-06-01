//
//  ReebComparator.cpp
//  ReebGraph
//
//  Created by Stephen John Russell on 25/05/2014.
//  Copyright (c) 2014 Stephen John Russell. All rights reserved.
//

#include "ReebComparator.hpp"
#include "vtkDataSet.h"

float ReebComparator::compareReebGraphs(vtkReebGraph *first, vtkReebGraph *second, vtkPolyData *fir, vtkPolyData *sec) {
    
    // Check which is bigger
    if (fir->GetNumberOfVerts() < sec->GetNumberOfVerts()) {
        vtkReebGraph *temp = first;
        first = second;
        second = temp;
    }
    
    // Get vert & arc info
    vtkVariantArray *edgeInfo = vtkVariantArray::SafeDownCast(first->GetEdgeData()->GetAbstractArray("Vertex Ids"));
    vtkVariantArray *edgeInfoTwo = vtkVariantArray::SafeDownCast(second->GetEdgeData()->GetAbstractArray("Vertex Ids"));
    
    vtkEdgeListIterator *eIt = vtkEdgeListIterator::New();
    first->GetEdges(eIt);
    vtkEdgeListIterator *eItS = vtkEdgeListIterator::New();
    second->GetEdges(eItS);
    
    // list of edges from both graphs
    std::vector<vtkEdgeType> edgeList;
    std::vector<vtkEdgeType> edgeListS;
    
    // Get values
    do {
        vtkEdgeType e = eIt->Next();
        edgeList.push_back(e);
    } while(eIt->HasNext());
    do {
        vtkEdgeType e = eItS->Next();
        edgeListS.push_back(e);
    } while(eItS->HasNext());
    
    // Current max sub-graph
    int subGraphSize = 0;
    
    // Go to all nodes
    for (int i = 0; i < edgeList.size(); i++) {
        
        // Get Distance and number of degrees of arc for first graph
        vtkEdgeType ed = edgeList[i];
        vtkAbstractArray *deg2NodeListE = edgeInfo->GetPointer(ed.Id)->ToArray();
        int deg = (int)deg2NodeListE->GetNumberOfTuples();
        double * p1 = (double *) malloc(sizeof(double)*3);
        double * p2 = (double *) malloc(sizeof(double)*3);
        fir->vtkDataSet::GetPoint(ed.Source, p1);
        fir->vtkDataSet::GetPoint(ed.Target, p2);
        double squaredDistance = vtkMath::Distance2BetweenPoints(p1, p2);
        double dist = sqrt(squaredDistance);
        
        // Visit all point in second graph
        for (int j = 0; j < edgeListS.size(); j++) {
            // Get Distance and number of degrees of arc for second graph
            int tempSubGraphSize = 0;
            vtkEdgeType eS = edgeListS[j];
            vtkAbstractArray *deg2NodeListES = edgeInfoTwo->GetPointer(eS.Id)->ToArray();
            int degS = (int)deg2NodeListES->GetNumberOfTuples();
            double * p1a = (double *) malloc(sizeof(double)*3);
            double * p2b = (double *) malloc(sizeof(double)*3);
            sec->vtkDataSet::GetPoint(eS.Source, p1a);
            sec->vtkDataSet::GetPoint(eS.Target, p2b);
            double squaredDistanceE = vtkMath::Distance2BetweenPoints(p1a, p2b);
            double sdist = sqrt(squaredDistanceE);

            // Check within bounds
            if (deg - degS < minDeg && deg - degS > -minDeg) {
                if (dist - sdist < minDist && dist - sdist > -minDist) {
                    //printf("Dist : %f\n", sdist - dist);
                    tempSubGraphSize++;
                    bool matchLost = false;
                    // check not at end of graph
                    if (i + 1 < edgeList.size() && j + 1 < edgeListS.size()) {
                        // new list locations
                        int t = i + 1;
                        int s = j + 1;
                        // loop until sub graph runs out
                        while (!matchLost) {
                            // Get Distance and number of degrees of arc for first graph
                            vtkEdgeType et = edgeList[t++];
                            vtkAbstractArray *deg2NodeListET = edgeInfo->GetPointer(et.Id)->ToArray();
                            int degE = (int)deg2NodeListET->GetNumberOfTuples();
                            double * p1t = (double *) malloc(sizeof(double)*3);
                            double * p2t = (double *) malloc(sizeof(double)*3);
                            fir->vtkDataSet::GetPoint(et.Source, p1t);
                            fir->vtkDataSet::GetPoint(et.Target, p2t);
                            double squaredDistanceT = vtkMath::Distance2BetweenPoints(p1t, p2t);
                            double distt = sqrt(squaredDistanceT);

                            // Get Distance and number of degrees of arc for second graph
                            vtkEdgeType eSS = edgeListS[s++];
                            vtkAbstractArray *deg2NodeListESS = edgeInfoTwo->GetPointer(eSS.Id)->ToArray();
                            int degSS = (int)deg2NodeListESS->GetNumberOfTuples();
                            double * p1at = (double *) malloc(sizeof(double)*3);
                            double * p2bt = (double *) malloc(sizeof(double)*3);
                            sec->vtkDataSet::GetPoint(eSS.Source, p1at);
                            sec->vtkDataSet::GetPoint(eSS.Target, p2bt);
                            double squaredDistanceEt = vtkMath::Distance2BetweenPoints(p1at, p2bt);
                            double sdistt = sqrt(squaredDistanceEt);
                            
                            // Check within bounds
                            if (degE - degSS  < minDeg && degE - degSS > -minDeg) {
                                if (distt - sdistt  < minDist  && distt - sdistt > -minDist) {
                                    //printf(" Dist : %f\n", sdistt - distt);
                                    tempSubGraphSize++;
                                }
                            } else {
                                // Lets get out of here
                                matchLost = true;
                            }
                            // Dont got too far!!
                            if (t > edgeList.size() - 1 || s > edgeListS.size() - 1) {
                                matchLost = true;
                            }
                        }
                    }
                    // Update max graph size
                    if (tempSubGraphSize > subGraphSize) {
                        subGraphSize = tempSubGraphSize;
                    }
                }
            }
        }
    }
    
    // calulate max and min vals
    int maxG = (int)edgeListS.size();
    int minG = (int)edgeList.size();
    if (edgeList.size() > edgeListS.size()) {
        maxG = (int) edgeList.size();
        minG = (int) edgeListS.size();
    }
    float sgs = (float)subGraphSize/maxG;
    cout << "Sub Graph :: " << subGraphSize << "\nMax Point :: " << maxG << endl;
    cout << "Min Point :: " << minG << endl;
    return sgs;
}

// init class vars
void ReebComparator::init(int d, double di) {
    minDeg = d;
    minDist = di;
}