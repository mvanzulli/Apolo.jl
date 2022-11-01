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
export vtk_structured_write

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
