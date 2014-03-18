//
//  ImgToMesh.cpp
//  ReebGraph
//
//  Created by Stephen John Russell on 18/03/2014.
//  Copyright (c) 2014 Stephen John Russell. All rights reserved.
//

#include "ImgToMesh.hpp"

std::vector<int> cellIds;

// Catch mouse events
class MouseInteractorStyle : public vtkInteractorStyleTrackballCamera
{
public:
    static MouseInteractorStyle* New();
    
    MouseInteractorStyle()
    {
        selectedMapper = vtkSmartPointer<vtkDataSetMapper>::New();
        selectedActor = vtkSmartPointer<vtkActor>::New();
    }
    
    virtual void OnLeftButtonDown()
    {
        // Get the location of the click (in window coordinates)
        int* pos = this->GetInteractor()->GetEventPosition();
        
        vtkSmartPointer<vtkCellPicker> picker =
        vtkSmartPointer<vtkCellPicker>::New();
        picker->SetTolerance(0.0005);
        
        // Pick from this location.
        picker->Pick(pos[0], pos[1], 0, this->GetDefaultRenderer());
        
        double* worldPosition = picker->GetPickPosition();
        std::cout << "Cell id is: " << picker->GetCellId() << std::endl;
        
        if(picker->GetCellId() != -1)
        {
            cellIds.push_back((int) picker->GetCellId());
            std::cout << "Pick position is: " << worldPosition[0] << " " << worldPosition[1]
            << " " << worldPosition[2] << endl;
            
            vtkSmartPointer<vtkIdTypeArray> ids =
            vtkSmartPointer<vtkIdTypeArray>::New();
            ids->SetNumberOfComponents(1);
            ids->InsertNextValue(picker->GetCellId());
            
            vtkSmartPointer<vtkSelectionNode> selectionNode =
            vtkSmartPointer<vtkSelectionNode>::New();
            selectionNode->SetFieldType(vtkSelectionNode::CELL);
            selectionNode->SetContentType(vtkSelectionNode::INDICES);
            selectionNode->SetSelectionList(ids);
            
            vtkSmartPointer<vtkSelection> selection =
            vtkSmartPointer<vtkSelection>::New();
            selection->AddNode(selectionNode);
            
            vtkSmartPointer<vtkExtractSelection> extractSelection =
            vtkSmartPointer<vtkExtractSelection>::New();
            extractSelection->SetInputData(0, this->Data);
            extractSelection->SetInputData(1, selection);
            extractSelection->Update();
            
            // In selection
            vtkSmartPointer<vtkUnstructuredGrid> selected =
            vtkSmartPointer<vtkUnstructuredGrid>::New();
            selected->ShallowCopy(extractSelection->GetOutput());
            
            std::cout << "There are " << selected->GetNumberOfPoints()
            << " points in the selection." << std::endl;
            std::cout << "There are " << selected->GetNumberOfCells()
            << " cells in the selection." << std::endl;
            
            selectedMapper->SetInputData(selected);
            
            selectedActor->SetMapper(selectedMapper);
            selectedActor->GetProperty()->EdgeVisibilityOn();
            selectedActor->GetProperty()->SetEdgeColor(1,0,0);
            selectedActor->GetProperty()->SetLineWidth(3);
            
            this->Interactor->GetRenderWindow()->GetRenderers()->GetFirstRenderer()->AddActor(selectedActor);
            
        }
        // Forward events
        vtkInteractorStyleTrackballCamera::OnLeftButtonDown();
    }
    
    vtkSmartPointer<vtkPolyData> Data;
    vtkSmartPointer<vtkDataSetMapper> selectedMapper;
    vtkSmartPointer<vtkActor> selectedActor;
    
};

void ImgToMesh::init() {
    
    fileName = "/Users/sjr/Pictures/gun/gun_0s.png";
    
    reader = vtkSmartPointer<vtkPNGReader>::New();
    quantizer = vtkSmartPointer<vtkImageQuantizeRGBToIndex>::New();
    
    img2data = vtkSmartPointer<vtkImageToPolyDataFilter>::New();
    
    
    triangleFilter = vtkSmartPointer<vtkTriangleFilter>::New();
    
    connector = vtkSmartPointer<vtkPolyDataConnectivityFilter>::New();
    
    
    cleaner = vtkSmartPointer<vtkCleanPolyData>::New();
    
    mesh = vtkSmartPointer<vtkPolyData>::New();
    mesh->Allocate();
    
}

vtkStandardNewMacro(MouseInteractorStyle);

bool ImgToMesh::loadFile() {
    
    reader->SetFileName("/Users/sjr/Pictures/gun/gun_0s.png");
    reader->Update();
    quantizer->SetInputConnection(reader->GetOutputPort());
    quantizer->SetNumberOfColors(2);
    quantizer->Update();
    img2data->SetInputConnection(quantizer->GetOutputPort());
    img2data->SetLookupTable(quantizer->GetLookupTable());
    img2data->SetColorModeToLUT();
    img2data->SetSmoothing(1);
    img2data->SetOutputStyleToPolygonalize();
    img2data->SetError(0);
    img2data->DecimationOn();
    img2data->SetDecimationError(0.0);
    img2data->SetSubImageSize(5);
    img2data->Update();
    mesh->DeepCopy(img2data->GetOutput());
    
    bool culled = false;
    
    while (!culled) {
        displayMesh();
        
        mesh->BuildLinks();
        for (int i = 0; i < cellIds.size(); i++) {
            mesh->DeleteCell(cellIds[i]);
        }
        mesh->RemoveDeletedCells();
        if (cellIds.size() > 0) {
            cellIds.clear();
        } else {
            culled = true;
        }
    }
    
    cout << "Cells :: " << mesh->GetNumberOfCells() << endl;
    cout << "Points :: " << mesh->GetNumberOfPoints() << endl;
    
    triangleFilter->SetInputData( mesh );
    triangleFilter->PassLinesOff();
    triangleFilter->PassVertsOff();
    triangleFilter->Update();
    
    connector->SetInputConnection(triangleFilter->GetOutputPort());
    connector->SetExtractionModeToAllRegions();
    connector->Update();
    
    cleaner->SetInputConnection( connector->GetOutputPort() );
    cleaner->Update();
    
    mesh->DeepCopy(cleaner->GetOutput());
    
    displayMesh();
    
    return true;
}

void ImgToMesh::displayMesh() {
    vtkRenderer *renderer = vtkRenderer::New();
    
    renderer->SetBackground(0.3,0.3,1);
    
    vtkRenderWindow *renderWindow = vtkRenderWindow::New();
    renderWindow->AddRenderer(renderer);
    renderWindow->SetSize(800, 800);
    
    vtkRenderWindowInteractor *windowInteractor =
    vtkRenderWindowInteractor::New();
    windowInteractor->SetRenderWindow(renderWindow);
    
    vtkSmartPointer<vtkPolyDataMapper> mapper =
    vtkSmartPointer<vtkPolyDataMapper>::New();
    mapper->SetInputData(mesh);
    
    vtkSmartPointer<vtkActor> actor =
    vtkSmartPointer<vtkActor>::New();
    actor->SetMapper(mapper);
    actor->GetProperty()->SetRepresentationToWireframe();
    
    windowInteractor->Initialize();
    vtkSmartPointer<MouseInteractorStyle> style =
    vtkSmartPointer<MouseInteractorStyle>::New();
    style->SetDefaultRenderer(renderer);
    style->Data = mesh;
    
    windowInteractor->SetInteractorStyle(style);
    
    renderer->AddActor(actor);
    renderer->ResetCamera();
    
    windowInteractor->Start();
}