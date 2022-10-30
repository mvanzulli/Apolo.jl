###########################################################
# Main types and functions to handle with analytic images #
###########################################################

export AnalyticImage

""" Analytic image struct.

### Fields:
- `intensity`  -- Intensity function depending on x, y, z
- `start`      -- Start coordinates (considering the start located at [1,1,1] )
- `finish`    -- Space in each direction between each pixel
"""
struct AnalyticImage{D,T} <: AbstractImage{D,T}
    intensity::Function
    start::NTuple{D,T}
    finish::NTuple{D,T}
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

        new{D,T}(intensity, start, finish)

    end

end

Base.length(a_img::AnalyticImage) = finish(a_img) .- start(a_img)

"Internal function to evaluate the image intensity at a generic point"
function _eval_intensity(
    p::NTuple{D,T},
    a_img::AnalyticImage{D};
    offset::NTuple{D,T} = Tuple(zeros(T,D))
) where {D,T}

    p = p .+ offset
    intensity_p = if p ⊂ a_img
        intensity_function = intensity(a_img)
        intensity_function(p...)
    else
        throw(ArgumentError("p + offset = $p is not inside the img frame"))
    end

    return intensity_p
end

"Internal function to evaluate the image intensity at a vector of points"
function _eval_intensity(
    vec_points::Vector{NTuple{D,T}},
    a_img::AnalyticImage{D};
    offset::NTuple{D,T} = Tuple(zeros(T,D))
) where {D,T}

    intensity_vec = Vector{T}(undef, length(vec_points))

    intensity_function = intensity(a_img)

    # Check if p is inside the grid image
    for (num_point, p) in enumerate(vec_points)
        p = p .+ offset
        if p ⊂ a_img
            intensity_vec[num_point] = intensity_function(p...)
        else
            throw(Warning("p + offset = $p is not inside the img frame"))
            return  missing
        end
    end

    return intensity_vec

end
