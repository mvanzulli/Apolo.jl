"""
Module defining different type of data measured. The model parameters will be computed
to match the measured data.
"""
module DataMeasured

using Apolo.Geometry: AbstractStructuredGrid
using Apolo.Geometry: dimension, grid, node_type
using Apolo.Geometry.FerriteGrids: FerriteStructuredGrid
using Apolo.Images: AbstractImage

using Ferrite: addnodeset!, getnodeset, getcoordinates

import Apolo.Geometry: grid
import Apolo.InverseProblem: roi

export AbstractDataMeasured, ImageData
export reference_img, deformed_imgs, roi_nodes, roi_nodes_coords


""" Abstract supertype for data measured.

The following methods are provided by the interface:

- `grid(datam)`       -- returns the grid where the data is measured.
- elapsed_time(datam) -- returns the elapsed time where the data was measured.
"""
abstract type AbstractDataMeasured{T} end

"Gets the grid where the data is measured"
grid(datam::AbstractDataMeasured) = datam.grid

"Gets the region where the data is measured"
roi(datam::AbstractDataMeasured) = datam.roi

"Gets the time where the data is measured"
elapsed_time(datam::AbstractDataMeasured) = datam.mtime

""" Image data struct.
### Fields:

- `vec_img`      -- vector of images
- `roi`          -- region of interest
- `grid`         -- region of interest
"""
struct ImageData{G<:AbstractStructuredGrid,Img<:AbstractImage,R,T<:AbstractVector} <: AbstractDataMeasured{Img}
    vec_img::Vector{Img}
    roi::R
    grid::G
    mtime::T
end

"Image data constructor wth a vector of images and a region of interest `roi`."
function ImageData(vec_img::Vector{I}, roi::R, t) where {R,I<:AbstractImage}

    fgrid = grid(vec_img[1])
    addnodeset!(grid(fgrid), "roi", roi)

    return ImageData(vec_img, roi, fgrid, t)
end

"Returns the reference image."
reference_img(img_data::ImageData) = img_data.vec_img[begin]

"Returns deformed images."
deformed_imgs(img_data::ImageData) = img_data.vec_img[2:end]

"Returns in a vector the node coordinates for each node in the roi."
function roi_nodes(img_data::ImageData{FG}) where {FG<:FerriteStructuredGrid}

    fgrid = grid(grid(img_data))

    return getnodeset(fgrid, "roi")

end

"Returns in a vector the node coordinates for each node in the roi."
function roi_nodes_coords(img_data::ImageData{FG}) where {FG<:FerriteStructuredGrid}

    # roi indexes
    roi_nodes_idx = roi_nodes(img_data)

    # roi coordinates
    fgrid = grid(grid(img_data))
    roi_coords = Vector{NTuple{dimension(grid(img_data)),node_type(grid(img_data))}}(undef, 0)
    for roi_index in roi_nodes_idx
        push!(roi_coords, Tuple(getcoordinates(fgrid.nodes[roi_index])))
    end

    return roi_coords
end

end #endmodule
