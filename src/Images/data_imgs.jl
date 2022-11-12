########################################################################
# Main types and functions to handle with a squence of measured images #
########################################################################

using ..Geometry: AbstractStructuredGrid, dimension, node_type
using ..Images: AbstractDataMeasured, AbstractImage
using Ferrite: addnodeset!, getnodeset, getcoordinates

export reference_img, deformed_imgs, roi_nodes, roi_nodes_coords

""" Image data struct.
### Fields:

- `vec_img`      -- vector of images
- `roi`          -- region of interest
- `grid`         -- region of interest
"""
struct ImageData{G<:AbstractStructuredGrid, Img<:AbstractImage, R, T<:AbstractVector} <:AbstractDataMeasured{Img}
    vec_img::Vector{Img}
    roi::R
    grid::G
    mtime::T
end

""
function ImageData(vec_img::Vector{I}, roi::R, t) where {R, I<:AbstractImage}

    fgrid = grid(vec_img[1])
    addnodeset!(grid(fgrid), "roi", roi)

    return ImageData(vec_img, roi, fgrid, t)
end

"Returns the reference image."
reference_img(img_data::ImageData) = img_data.vec_img[begin]

"Returns deformed images."
deformed_imgs(img_data::ImageData) = img_data.vec_img[2:end]

"Returns in a vector the node coordinates for each node in the roi."
function roi_nodes(img_data::ImageData{FG}) where {FG <: FerriteStructuredGrid}

    fgrid = grid(grid(img_data))

    return getnodeset(fgrid, "roi")

end

"Returns in a vector the node coordinates for each node in the roi."
function roi_nodes_coords(img_data::ImageData{FG}) where {FG <: FerriteStructuredGrid}

    # roi indexes
    roi_nodes_idx =  roi_nodes(img_data)

    # roi coordinates
    fgrid = grid(grid(img_data))
    roi_coords = Vector{NTuple{dimension(grid(img_data)), node_type(grid(img_data))}}(undef,0)
    for roi_index in roi_nodes_idx
        push!(roi_coords, Tuple(getcoordinates(fgrid.nodes[roi_index])))
    end

    return roi_coords
end
