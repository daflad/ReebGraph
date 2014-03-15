#define vtkRenderingCore_AUTOINIT 4(vtkInteractionStyle,vtkRenderingFreeType,vtkRenderingFreeTypeOpenGL,vtkRenderingOpenGL)
#define vtkRenderingVolume_AUTOINIT 1(vtkRenderingVolumeOpenGL)

#include <map>

#include "vtkActor.h"
#include "vtkAreaContourSpectrumFilter.h"
#include "vtkCamera.h"
#include "vtkCleanPolyData.h"
#include "vtkClipPolyData.h"
#include "vtkColor.h"
#include "vtkDataSetAttributes.h"
#include "vtkDoubleArray.h"
#include "vtkEdgeListIterator.h"
#include "vtkIdList.h"
#include "vtkImageToPolyDataFilter.h"
#include "vtkImageQuantizeRGBToIndex.h"
#include "vtkLight.h"
#include "vtkLookupTable.h"
#include "vtkObjectFactory.h"
#include "vtkPNGReader.h"
#include "vtkPointData.h"
#include "vtkPolyData.h"
#include "vtkPolyDataConnectivityFilter.h"
#include "vtkPolyDataMapper.h"
#include "vtkPolyDataToReebGraphFilter.h"
#include "vtkProperty.h"
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
#include "vtkTriangleFilter.h"
#include "vtkUnstructuredGrid.h"
#include "vtkUnstructuredGridToReebGraphFilter.h"
#include "vtkVariantArray.h"
#include "vtkVolumeContourSpectrumFilter.h"
#include "vtkWindowToImageFilter.h"



// loading code generated automatically.
// put whatever mesh you like here.
int LoadSurfaceMesh(vtkPolyData *sMesh)
{
    sMesh->Allocate();
    
    vtkSmartPointer<vtkPNGReader> reader =
    vtkSmartPointer<vtkPNGReader>::New();
    reader->SetFileName("/Users/sjr/Pictures/gun/gun_1s.png");
    
    vtkSmartPointer<vtkImageQuantizeRGBToIndex> quant =
    vtkSmartPointer<vtkImageQuantizeRGBToIndex>::New();
    quant->SetInputConnection(reader->GetOutputPort());
    quant->SetNumberOfColors(16);
    
    vtkSmartPointer<vtkImageToPolyDataFilter> i2pd =
    vtkSmartPointer<vtkImageToPolyDataFilter>::New();
    i2pd->SetInputConnection(quant->GetOutputPort());
    i2pd->SetLookupTable(quant->GetLookupTable());
    i2pd->SetColorModeToLUT();
    i2pd->SetOutputStyleToPolygonalize();
    i2pd->SetError(0);
    i2pd->DecimationOn();
    i2pd->SetDecimationError(0.0);
    i2pd->SetSubImageSize(5);
    
    vtkSmartPointer<vtkClipPolyData> clip =
    vtkSmartPointer<vtkClipPolyData>::New();
    clip->SetGenerateClipScalars(1);
    clip->SetInputConnection(i2pd->GetOutputPort());
    clip->GenerateClippedOutputOn();
    clip->InsideOutOn();
    clip->Update();
    
    vtkSmartPointer<vtkTriangleFilter> tfa =
    vtkSmartPointer<vtkTriangleFilter>::New();
    tfa->SetInputConnection(clip->GetOutputPort());
    tfa->PassLinesOff();
    tfa->PassVertsOff();
    tfa->Update();
    
    
    vtkSmartPointer<vtkPolyDataConnectivityFilter> oCon = vtkSmartPointer<vtkPolyDataConnectivityFilter>::New();
    oCon->SetInputData( tfa->GetOutput() );
    oCon->SetExtractionModeToAllRegions();
    oCon->Update();
    
    vtkSmartPointer<vtkCleanPolyData> oC2 = vtkSmartPointer<vtkCleanPolyData>::New();
    oC2->SetInputData( oCon->GetOutput() );
    oC2->Update();
    
    sMesh->DeepCopy(oC2->GetOutput());
    
    cout << "Cells :: " << sMesh->GetNumberOfCells() << endl;
    cout << "Points :: " << sMesh->GetNumberOfPoints() << endl;
    
    return 1;
}


class AreaSimplificationMetric : public vtkReebGraphSimplificationMetric{
public:
    vtkTypeMacro(AreaSimplificationMetric, vtkReebGraphSimplificationMetric);
    static AreaSimplificationMetric* New();
    double ComputeMetric(vtkDataSet *mesh, vtkDataArray *scalarField,
                         vtkIdType startCriticalPoint, vtkAbstractArray *vertexList,
                         vtkIdType endCriticalPoint);
};

vtkStandardNewMacro(AreaSimplificationMetric);

double AreaSimplificationMetric::ComputeMetric(vtkDataSet *mesh,
                                               vtkDataArray *scalarField, vtkIdType startCriticalPoint,
                                               vtkAbstractArray* vertexList, vtkIdType endCriticalPoint)
{
    // In this example, the metric algorithm just evaluates the area of the
    // surface region corresponding to the arc of the Reeb graph passed as an
    // argument.
    // As a result, the arcs corresponding to small surface regions (below the
    // threshold specified to the simplificatin filter) will be
    // simplified in priority in the surface simplification algorithm.
    
    double  fieldLowerBound = scalarField->GetComponent(startCriticalPoint,0),
    fieldUpperBound = scalarField->GetComponent(endCriticalPoint,0);
    
    double  cumulativeArea = 0;
    
    std::map<vtkIdType, bool> visitedTriangles;
    
    for(int i = 0; i < vertexList->GetNumberOfTuples(); i++)
    {
        int vId = vertexList->GetVariantValue(i).ToInt();
        vtkIdList *starTriangleList = vtkIdList::New();
        
        mesh->GetPointCells(vId, starTriangleList);
        
        for(int j = 0; j < starTriangleList->GetNumberOfIds(); j++)
        {
            vtkIdType tId = starTriangleList->GetId(j);
            vtkTriangle *t = vtkTriangle::SafeDownCast(mesh->GetCell(tId));
            std::map<vtkIdType, bool>::iterator tIt = visitedTriangles.find(tId);
            if(tIt == visitedTriangles.end())
            {
                if((scalarField->GetComponent(t->GetPointIds()->GetId(0), 0)
                    <= fieldUpperBound)
                   &&(scalarField->GetComponent(t->GetPointIds()->GetId(1), 0)
                      <= fieldUpperBound)
                   &&(scalarField->GetComponent(t->GetPointIds()->GetId(2), 0)
                      <= fieldUpperBound)
                   &&(scalarField->GetComponent(t->GetPointIds()->GetId(0), 0)
                      >= fieldLowerBound)
                   &&(scalarField->GetComponent(t->GetPointIds()->GetId(1), 0)
                      >= fieldLowerBound)
                   &&(scalarField->GetComponent(t->GetPointIds()->GetId(2), 0)
                      >= fieldLowerBound))
                {
                    // the triangle fully maps inside the arc function interval
                    cumulativeArea += t->ComputeArea();
                }
                visitedTriangles[tId] = true;
            }
        }
        
        starTriangleList->Delete();
    }
    double c = 1 - ( cumulativeArea/(this->UpperBound - this->LowerBound));
    cout << "Eveluated :: " << c << endl;
    return c;
}

int DisplayReebGraph(vtkReebGraph *g)
{
    vtkDataArray *vertexInfo = vtkDataArray::SafeDownCast(
                                                          g->GetVertexData()->GetAbstractArray("Vertex Ids"));
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

int DisplaySurfaceSkeleton(vtkPolyData *surfaceMesh, vtkTable *skeleton)
{
    
    // Rendering setting
    vtkRenderer *renderer = vtkRenderer::New();
    
    vtkRenderWindow *renderWindow = vtkRenderWindow::New();
    renderWindow->AddRenderer(renderer);
    renderWindow->SetSize(800, 800);
    
    vtkRenderWindowInteractor *windowInteractor =
    vtkRenderWindowInteractor::New();
    windowInteractor->SetRenderWindow(renderWindow);
    
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
    
    windowInteractor->Initialize();
    
    // Interactive mode
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
    renderWindow->Delete();
    renderer->Delete();
    
    return 0;
}

int main(int vtkNotUsed(argc), char* vtkNotUsed(argv)[] )
{
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
    surfaceSimplification->SetSimplificationThreshold(0.01);
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

