"""
Module defining analytic images (an image with analytical intensity function).
"""
module AnalyticImages

using ..Images: AbstractImage

import ..Images: finish, _eval_intensity, intensity, start

export AnalyticImage

""" Analytic image struct.

### Fields:
- `intensity`  -- Intensity function depending on x, y, z
- `start`      -- Start coordinates (considering the start located at [1,1,1] )
- `finish`    -- Space in each direction between each pixel
"""
struct AnalyticImage{D,T,G<:Nothing} <: AbstractImage{D,T,G}
    intensity::Function
    start::NTuple{D,T}
    finish::NTuple{D,T}
    grid::G
    function AnalyticImage(
        intensity::Function,
        start::NTuple{D,T},
        finish::NTuple{D,T},
    ) where {D,T}

        try
            intensity(Tuple(rand(D))...)
        catch
            throw(ArgumentError("The intensity function has not D nargs"))
        end

        new{D,T,Nothing}(intensity, start, finish, nothing)

    end

end

"Returns the final point of an analytic `a_img`."
finish(a_img::AnalyticImage) = a_img.finish

"Returns the length of an analytic image `a_img`."
Base.length(a_img::AnalyticImage) = finish(a_img) .- start(a_img)

"Internal function to evaluate the image intensity at a generic point `p` with an `offset`."
function _eval_intensity(
    p::NTuple{D,T},
    a_img::AnalyticImage{D};
    offset::NTuple{D,T}=Tuple(zeros(T, D))
) where {D,T}

    p = p .+ offset
    intensity_p = if p ∈ a_img
        intensity_function = intensity(a_img)
        intensity_function(p...)
    else
        # @warn "p + offset = $p is not inside the img frame"
    end

    return intensity_p
end

"Internal function to evaluate the image intensity at a vector of point `vec_points` with an `offset`."
function _eval_intensity(
    vec_points::Vector{NTuple{D,T}},
    a_img::AnalyticImage{D};
    offset::NTuple{D,T}=Tuple(zeros(T, D))
) where {D,T}

    intensity_vec = Vector{T}(undef, length(vec_points))

    intensity_function = intensity(a_img)

    # Check if p is inside the grid image
    for (num_point, p) in enumerate(vec_points)
        p = p .+ offset
        if p ∈ a_img
            intensity_vec[num_point] = intensity_function(p...)
        else
            # @warn "p + offset = $p is not inside the img frame"
            return missing
        end
    end

    return intensity_vec

end

end #endmodule
