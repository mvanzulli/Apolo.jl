using Apolo

using WriteVTK

dcm_dir = "./test/DICOMImages/"

img = load_dicom(dcm_dir)


filename = "foo"
# function vtkStrGridPlot( grid, nodalMagnitude, filename )
intensity = getintensity(img)
x,y,z = build_coordinates(img)

  # vtk_grid( filename, x, y, z,  append = false, ascii=true ) do vtk
  vtk_grid( filename, x,y,z) do vtk
    #for i in (1:length(nodalMagnitudes))
    #  vtk[ string("intensity",i), VTKPointData() ] = nodalMagnitudes[i]
    #end
    vtk[ "intensity", VTKPointData() ] = intensity
  end



intensity  = getintensity(imgs)