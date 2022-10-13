########################################
# Functions to read and plot vtk files #
########################################

# Import dependencies to overlead

# Add libraries to use
# Internal
using Apolo.Images
# External
using WriteVTK: VTKPointData, vtk_grid

# Export file functions
export vtk_write

##########
# Images #
##########
"Function to write a medical image with structured grid into a .vtk file"
function vtk_write(med_img::MedicalImage, filename = get_patient_name(med_img))
    # get image coordinates
    x, y, z = build_coordinates(med_img)
    # get image intensity
    intensity = getintensity(med_img)
    # plot vtk 
    vtk_grid( filename, x, y, z) do vtk
        vtk[ "Intensity", VTKPointData() ] = intensity
    end
end