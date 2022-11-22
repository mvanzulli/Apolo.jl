########################################
# Functions to read and plot vtk files #
########################################

using Apolo.ForwardProblem.FerriteSolver: dofhandler
using Apolo.Images
using NaturalSort: natural
using WriteVTK: VTKPointData, vtk_grid, vtk_point_data
using ReadVTK: VTKFile, get_whole_extent, get_origin, get_spacing, get_point_data, get_data

export write_vtk_fsol, vtk_structured_write, vtk_structured_write_sequence, load_vtk_img, load_vtk_sequence_imgs

###############
# .VTI format #
###############

"Write a vtk structured grid given a 2D or 3D scalar array."
function vtk_structured_write(
    coords::AbstractVector,
    fieldarray::FA,
    fieldname::Symbol=:generic,
    filename::String="generic",
    filedir::String="./",
) where {FA<:Union{AbstractArray{<:Number,3},AbstractArray{<:Number,2}}}

    # Check dimensions.
    prod(length.(coords)) == length(fieldarray) || throw(ArgumentError("Check dimensions"))

    # Generate a VTK file.
    vtk_grid(joinpath(filedir, filename), coords...) do vtk
        vtk[String(fieldname), VTKPointData()] = fieldarray
    end
end

"Write a .VTK structured given a field function."
function vtk_structured_write(
    coords::AbstractVector,
    fieldfunction::Function,
    fieldname::Symbol=:generic,
    filename::String="generic",
    filedir::String="./",
)
    fieldarray = [fieldfunction(c...) for c in Iterators.product(coords...)]
    vtk_structured_write(coords, fieldarray, fieldname, filename, filedir)
end

"Write a vtk structured grid given a 4D scalar array."
function vtk_structured_write_sequence(
    coords::AbstractVector,
    fieldarray::Array{<:Number,D},
    fieldname::Symbol=:generic,
    filename::String="generic",
    filedir::String="./",
) where {D}

    # Plot a sequence of .vti
    for iseq in 1:size(fieldarray, D)
        # add the current time to the filename
        filename_slice = filename * "_" * string(iseq)
        #create a variable colomn argumets
        colums_aux = Vector{Function}(undef, D - 1)
        fill!(colums_aux, :)
        #create the array
        slice_int_array = view(fieldarray, colums_aux..., iseq)
        vtk_structured_write(coords, slice_int_array, fieldname, filename_slice, filedir)

    end
end

"Write a vtk structured grid given a 3D/4D field function."
function vtk_structured_write_sequence(
    vars::AbstractVector,
    fieldfunction::Function,
    fieldname::Symbol=:generic,
    filename::String="generic",
    filedir::String="./",
)

    fieldarray = [fieldfunction(v...) for v in Iterators.product(vars...)]
    coords = vars[1:length(vars)-1]
    vtk_structured_write_sequence(coords, fieldarray, fieldname, filename, filedir)

    return nothing
end

""" Reads a .VTK image with an and create a grid. """
function load_vtk_img(filename::String)
    # Check if it has the .vti extension and if not add it.
    path_with_extension = if occursin(".vti", filename)
        filename
    else
        filename * ".vti"
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
        intensity_array, spacing_img, start_img, filename, ferrite_grid=true
    )

    return vtk_img
end

""" Reads a .VTK image with an and create a grid. """
function load_vtk_img(filename::String, grid::AbstractStructuredGrid)
    # Check if it has the .vti extension and if not add it.
    path_with_extension = if occursin(".vti", filename)
        filename
    else
        filename * ".vti"
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
    intensity_array = reshape(intensity_vec, num_pixels_img)

    vtk_img = VTKImage(
        collect(intensity_array), grid, spacing_img, start_img, filename)

    return vtk_img
end

"Loads a VTK sequence of images"
function load_vtk_sequence_imgs(folder_path::String)

    # Extract.vti or .vtk list
    files = readdir(folder_path)
    vtkfiles = sort(files[endswith.(files, r".vti|vtk")], lt=natural)
    vtkfiles = joinpath.(folder_path, vtkfiles)

    # Load reference image
    ref_img = load_vtk_img(vtkfiles[1])
    images_grid = grid(ref_img)

    # Build the vector of images and fill it
    imgs = [ref_img]
    for def_img_fname in vtkfiles[2:end]
        def_img = load_vtk_img(def_img_fname, images_grid)
        push!(imgs, def_img)
    end


    return imgs
end

"Generates the vtk of a forward problem solution."
function write_vtk_fsol(
    sol::ForwardProblemSolution{S},
    dir::String,
    filename::String,
) where {S<:FerriteForwardSolver,D,T}

    # Create and eval the dof handler
    dh = dofhandler(sol)

    vtk_grid(dir * filename, dh) do vtkfile
        vtk_point_data(vtkfile, dh, dofsvals(sol))
    end

    return "VTK generated at :$dir"
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
