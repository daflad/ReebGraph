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
#include "ReebComparator.hpp"

vtkStandardNewMacro(AreaSimplificationMetric);

// Main function comparing two reeb graphs
int main(int argc, char* argv[] ) {
    
    // First Image
    ImgToMesh im;
    im.init();
    im.meshName = argv[1];
    im.loadFile();
    // First Graph
    ReebGrapher rg;
    rg.init(im.mesh);
    rg.attachScalarField(atoi(argv[3]));
    rg.buildReebGraph();
    rg.DisplaySurfaceSkeleton();
    // Second Image
    ImgToMesh im2;
    im2.init();
    im2.meshName = argv[2];
    im2.loadFile();
    // Second Graph
    ReebGrapher rgB;
    rgB.init(im2.mesh);
    rgB.attachScalarField(atoi(argv[3]));
    rgB.buildReebGraph();
    rgB.DisplaySurfaceSkeleton();
    // Compare Graphs
    ReebComparator rc;
    rc.init(atoi(argv[4]), atoi(argv[5]));
    float sgs = rc.compareReebGraphs(rg.simplifiedSurfaceReebGraph, rgB.simplifiedSurfaceReebGraph, rg.mesh, rgB.mesh);
    // Show result
    printf("Similarity :: %f", sgs);
    
    return 0;
}

