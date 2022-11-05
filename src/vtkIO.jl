########################################
# Functions to read and plot vtk files #
########################################

using Apolo.Images

using WriteVTK: VTKPointData, vtk_grid
using ReadVTK: VTKFile, get_whole_extent, get_origin, get_spacing, get_point_data, get_data

export vtk_structured_write, load_vtk_img

##########
# Images #
##########

"Write a .VTK structured given a field function."
function vtk_structured_write(
    coords::Vector{<:AbstractRange},
    fieldfunction::Function,
    fieldname::Symbol=:generic,
    filename::String="generic"
)
    fieldarray = [fieldfunction(c...) for c in Iterators.product(coords...)]
    vtk_structured_write(coords, fieldarray, fieldname, filename)
end

"Write a vtk structured grid given a 3D scalar array"
function vtk_structured_write(
    coords::Vector{<:AbstractRange},
    fieldarray::Array,
    fieldname::Symbol=:generic,
    filename::String="generic"
)
    # Check dimensions.
    mapreduce(length, isequal, coords) || throw(ArgumentError("Check dimensions"))

    # Generate a VTK file.
    vtk_grid(filename, coords...) do vtk
        vtk[String(fieldname), VTKPointData()] = fieldarray
    end
end

"""
Reads a .VTK image with an structured grid.
"""
function load_vtk_img(path_img::String)
    # Check if it has the .vti extension and if not add it.
    path_with_extension = if occursin(".vti", path_img)
        path_img
    else
        path_img * ".vti"
    end
    vtk = VTKFile(path_with_extension)

    # Extract image dimensions.
    start_img_grid = Tuple(get_origin(vtk))
    spacing_img = Tuple(get_spacing(vtk))
    start_img = @. start_img_grid - spacing_img / 2

    # Compute the number of voxels and if it is zero then delete it.
    num_voxels = get_whole_extent(vtk)[2:2:end]
    zero_voxels_axis = findall(iszero, num_voxels)
    [popat!(num_voxels, i) for i in zero_voxels_axis]

    if length(num_voxels) == 2
        num_pixels_img = Tuple(num_voxels .+ (1, 1))
    elseif length(num_voxels) == 3
        num_pixels_img = Tuple(num_voxels .+ (1, 1, 1))
    else
        throw(ArgumentError("Dimension error"))
    end

    # Build intensity array.
    intensity_vec = get_data(get_point_data(vtk)["intensity"])
    intensity_array = reshape(intensity_vec, num_pixels_img) |> collect

    vtk_img = VTKImage(
        intensity_array, spacing_img, start_img, path_img, ferrite_grid=true
    )

    return vtk_img
end

"Function to write a medical image with structured grid into a .vtk file"
function vtk_write(med_img::MedicalImage, filename=get_patient_name(med_img))
    # Get image coordinates.
    x, y, z = build_coordinates(med_img)
    # Get image intensity.
    intensity = getintensity(med_img)
    # Plot VTK.
    vtk_grid(filename, x, y, z) do vtk
        vtk["Intensity", VTKPointData()] = intensity
    end
end
