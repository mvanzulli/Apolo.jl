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

"Write a .VTK structured given a scalar_field function."
function vtk_structured_write(
    coords::Vector{<:AbstractRange},
    fieldfunction::Function,
    fieldname::Symbol = :generic,
    filename::String = "generic"
    )

    # dimension
    dimension = length(coords)

    if dimension == 3
        # Extract ranges
        x, y, z = coords

        scalar_array = [fieldfunction(xᵢ,yᵢ,zᵢ) for xᵢ in x for yᵢ in y for zᵢ in z]

        # generate a vtk
        vtk_grid(filename, x, y, z) do vtk
            vtk[ String(fieldname), VTKPointData() ] = scalar_array
        end

    elseif dimension == 2
        # Extract ranges
        x, y = coords

        scalar_array = [fieldfunction(xᵢ,yᵢ) for xᵢ in x for yᵢ in y ]

        # generate a vtk
        vtk_grid(filename, x, y) do vtk
            vtk[ String(fieldname), VTKPointData() ] = scalar_array
        end

    end

end

"Write a vtk structured grid given a 3D scalar array"
function vtk_structured_write(
    coords::Vector{<:AbstractRange},
    fieldarray::Array,
    fieldname::Symbol = :generic,
    filename::String = "generic"
    )

    dimension = length(coords)


    if dimension == 3

        # extract coordinates in each direction
        x, y, z = coords

        (length(x), length(y), length(z)) != size(fieldarray) && throw(ArgumentError("Check dimensions"))
        # Extract ranges
        x, y, z = coords

        # generate a vtk
        vtk_grid(filename, x, y, z) do vtk
            vtk[ String(fieldname), VTKPointData() ] = fieldarray
        end

    elseif dimension == 2
        # Extract ranges
        x, y = coords

        (length(x), length(y)) != size(fieldarray) && throw(ArgumentError("Check dimensions"))

        # generate a vtk
        vtk_grid(filename, x, y) do vtk
            vtk[ String(fieldname), VTKPointData() ] = fieldarray
        end

    end

end


"""
Reads a .VTK image with an structured grid.
"""
function load_vtk_img(path_img::String)

    # Checks if it has the .vti extension and if not added
    extension = ".vti"
    if occursin(extension, path_img)
        path_with_extension = path_img
        vtk = VTKFile(path_with_extension)
    else
        path_with_extension = path_img * ".vti"
        vtk = VTKFile(path_with_extension)
    end

    # Extract image dimensions
    start_img_grid = Tuple(get_origin(vtk))
    spacing_img = Tuple(get_spacing(vtk))
    start_img = start_img_grid .- spacing_img./2

    # Computes the number of voxels and if it is zero then delete it
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

    # Builds intensity array
    intensity_vec = get_data( get_point_data(vtk)["intensity"] )
    intensity_array = reshape(intensity_vec, num_pixels_img) |> collect

    vtk_img = VTKImage(
        intensity_array, spacing_img, start_img, path_img, ferrite_grid = true
        )

    return vtk_img

end




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
