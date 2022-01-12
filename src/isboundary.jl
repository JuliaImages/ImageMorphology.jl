"""
    isboundary!(img::AbstractArray; background = 0, dims = coords_spatial(A), kwargs...)

Finds the boundaries that are just within each object, replacing the original image.
`background` is the scalar value of the background pixels which will not be marked as boundaries.
Keyword arguments are passed to `extremefilt!` which include `dims` indicating the dimension(s)
over which to discover boundaries.

See `isboundary` for examples.
"""
function isboundary!(img::AbstractArray{Bool,N};
                         background=false,
                         kwargs...
                        ) where N
    # Find regions where there is at least one background pixel
    # centered on a foreground pixel
    if background == true
        # background_img = img
        # foreground_img = .!img
        img .= .!img .& extremefilt!(copy(img), |; kwargs...)
        # Alternatively, with ImageFiltering.jl
        # img .= .!img .& mapwindow(any, img, (3,3))
    elseif background == false
        # background_img = .!img
        # foreground_img = img
        img .&= extremefilt!(.!img, |; kwargs...)
        # Alternatively, with ImageFiltering.jl
        # img .&= mapwindow(any, .!img, (3,3))
    else
        # This should be the same as isboundary_thick
        # background is neither true or false, use more generic algorithm
        background = Int8(-1) # Any value other than true or false will do
        #allequal(x,y) = ifelse(x == y, x, background)
        allequal(x,y) = x == y ? x : background
        # the entire image is the foreground
        img .= (extremefilt!(Int8.(img), allequal; kwargs...) .== background)
    end
    return img
end
function isboundary!(img::AbstractArray{T};
                          background = zero(T),
                          kwargs...) where T
    # Find regions where that are not homogeneous (all equal)
    # centered on a foreground pixel
    # background_img = img .== background
    # foreground_img = img .!= background
    #allequal(x,y) = ifelse(x == y, x, background)
    allequal(x,y) = x == y ? x : background
    img .= (extremefilt!(copy(img), allequal; kwargs...) .== background) .& (img .!= background)
    # Alternatively, with ImageFiltering.jl
    # allequal(x) = all(==(first(x)), @view(x[2:end]))
    # .!mapwindow(allequal, A, (3,3)) .& (A .!= background)
    return img
end

"""
    isboundary(img::AbstractArray; background = 0, dims = coords_spatial(A), kwargs...)

Finds the boundaries that are just within each object.
`background` is the scalar value of the background pixels which will not be marked as boundaries.
Keyword arguments are passed to `extremefilt!` which include `dims` indicating the dimension(s)
over which to discover boundaries.

See also `isboundary_thick`.

# Examples

```@meta
DocTestSetup = quote
    import ImageMorphology: isboundary
end
```

```jldoctest
julia> A = zeros($Int, 16, 16); A[4:8, 4:8] .= 5; A[4:8, 9:12] .= 6; A[10:12,13:15] .= 3; A[10:12,3:6] .= 9; A
16×16 Matrix{$Int}:
 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
 0  0  0  5  5  5  5  5  6  6  6  6  0  0  0  0
 0  0  0  5  5  5  5  5  6  6  6  6  0  0  0  0
 0  0  0  5  5  5  5  5  6  6  6  6  0  0  0  0
 0  0  0  5  5  5  5  5  6  6  6  6  0  0  0  0
 0  0  0  5  5  5  5  5  6  6  6  6  0  0  0  0
 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
 0  0  9  9  9  9  0  0  0  0  0  0  3  3  3  0
 0  0  9  9  9  9  0  0  0  0  0  0  3  3  3  0
 0  0  9  9  9  9  0  0  0  0  0  0  3  3  3  0
 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0

julia> isboundary(A)
16×16 Matrix{$Int}:
 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
 0  0  0  1  1  1  1  1  1  1  1  1  0  0  0  0
 0  0  0  1  0  0  0  1  1  0  0  1  0  0  0  0
 0  0  0  1  0  0  0  1  1  0  0  1  0  0  0  0
 0  0  0  1  0  0  0  1  1  0  0  1  0  0  0  0
 0  0  0  1  1  1  1  1  1  1  1  1  0  0  0  0
 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
 0  0  1  1  1  1  0  0  0  0  0  0  1  1  1  0
 0  0  1  0  0  1  0  0  0  0  0  0  1  0  1  0
 0  0  1  1  1  1  0  0  0  0  0  0  1  1  1  0
 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0

julia> isboundary(A .!= 0)
16×16 BitMatrix:
 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
 0  0  0  1  1  1  1  1  1  1  1  1  0  0  0  0
 0  0  0  1  0  0  0  0  0  0  0  1  0  0  0  0
 0  0  0  1  0  0  0  0  0  0  0  1  0  0  0  0
 0  0  0  1  0  0  0  0  0  0  0  1  0  0  0  0
 0  0  0  1  1  1  1  1  1  1  1  1  0  0  0  0
 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
 0  0  1  1  1  1  0  0  0  0  0  0  1  1  1  0
 0  0  1  0  0  1  0  0  0  0  0  0  1  0  1  0
 0  0  1  1  1  1  0  0  0  0  0  0  1  1  1  0
 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0

julia> isboundary(A .!= 0; dims = 1)
16×16 BitMatrix:
 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
 0  0  0  1  1  1  1  1  1  1  1  1  0  0  0  0
 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
 0  0  0  1  1  1  1  1  1  1  1  1  0  0  0  0
 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
 0  0  1  1  1  1  0  0  0  0  0  0  1  1  1  0
 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
 0  0  1  1  1  1  0  0  0  0  0  0  1  1  1  0
 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0

julia> isboundary(A .!= 0; dims = 2)
16×16 BitMatrix:
 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
 0  0  0  1  0  0  0  0  0  0  0  1  0  0  0  0
 0  0  0  1  0  0  0  0  0  0  0  1  0  0  0  0
 0  0  0  1  0  0  0  0  0  0  0  1  0  0  0  0
 0  0  0  1  0  0  0  0  0  0  0  1  0  0  0  0
 0  0  0  1  0  0  0  0  0  0  0  1  0  0  0  0
 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
 0  0  1  0  0  1  0  0  0  0  0  0  1  0  1  0
 0  0  1  0  0  1  0  0  0  0  0  0  1  0  1  0
 0  0  1  0  0  1  0  0  0  0  0  0  1  0  1  0
 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
```
"""
function isboundary(img::AbstractArray{T};
                         background = zero(T),
                         kwargs...
                        ) where T
    isboundary!(copy(img); background = background, kwargs...)
end

## Alternate implementations using dilate and erode
# The implementations below are based on scikit-image

#=
Copyright (C) 2021, Mark Kittisopikul and JuliaImages Github Organization
Copyright (C) 2019, the scikit-image team
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are
met:

 1. Redistributions of source code must retain the above copyright
    notice, this list of conditions and the following disclaimer.
 2. Redistributions in binary form must reproduce the above copyright
    notice, this list of conditions and the following disclaimer in
    the documentation and/or other materials provided with the
    distribution.
 3. Neither the name of skimage nor the names of its contributors may be
    used to endorse or promote products derived from this software without
    specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE AUTHOR ``AS IS'' AND ANY EXPRESS OR
IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY DIRECT,
INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,
STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING
IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
POSSIBILITY OF SUCH DAMAGE.
=#

# Specific implementation of `isboundary_thick` using dilate and erode similar to scikit-image, uses extremefilt! twice
function isboundary_thick_dilate_erode(img::AbstractArray; kwargs...)
    # https://github.com/scikit-image/scikit-image/blob/d44ceda6241cb23a22dc8abf09a05090ed14da7f/skimage/segmentation/boundaries.py#L165-L167
    return dilate(img; kwargs...) .!= erode(img; kwargs...)
end
function isboundary_dilate_erode(img::AbstractArray{T};
                                      background = zero(T),
                                      kwargs...
                                     ) where T
    # https://github.com/scikit-image/scikit-image/blob/d44ceda6241cb23a22dc8abf09a05090ed14da7f/skimage/segmentation/boundaries.py#L165-L170
    thick_boundaries = isboundary_thick_dilate_erode(img; kwargs...)
    foreground_img = img .!= background
    return thick_boundaries .& foreground_img
end
"""
    isboundary_thick(img::AbstractArray; dims = coords_spatial(img), kwargs...)

Find thick boundaries that are just outside and just inside the objects.
This is a union of the inner and outer boundaries.
Keyword `dims` indicates over which dimensions to look for boundaries. This `dims` and additional keywords `kwargs` are passed to `extremefilt!`.

# Examples

```@meta
DocTestSetup = quote
    import ImageMorphology: isboundary_thick
end
```

```jldoctest
julia> A = zeros($Int, 16, 16); A[4:8, 4:8] .= 5; A[4:8, 9:12] .= 6; A[10:12,13:15] .= 3; A[10:12,3:6] .= 9; A
16×16 Matrix{$Int}:
 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
 0  0  0  5  5  5  5  5  6  6  6  6  0  0  0  0
 0  0  0  5  5  5  5  5  6  6  6  6  0  0  0  0
 0  0  0  5  5  5  5  5  6  6  6  6  0  0  0  0
 0  0  0  5  5  5  5  5  6  6  6  6  0  0  0  0
 0  0  0  5  5  5  5  5  6  6  6  6  0  0  0  0
 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
 0  0  9  9  9  9  0  0  0  0  0  0  3  3  3  0
 0  0  9  9  9  9  0  0  0  0  0  0  3  3  3  0
 0  0  9  9  9  9  0  0  0  0  0  0  3  3  3  0
 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0

julia> isboundary_thick(A)
16×16 BitMatrix:
 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
 0  0  1  1  1  1  1  1  1  1  1  1  1  0  0  0
 0  0  1  1  1  1  1  1  1  1  1  1  1  0  0  0
 0  0  1  1  0  0  0  1  1  0  0  1  1  0  0  0
 0  0  1  1  0  0  0  1  1  0  0  1  1  0  0  0
 0  0  1  1  0  0  0  1  1  0  0  1  1  0  0  0
 0  0  1  1  1  1  1  1  1  1  1  1  1  0  0  0
 0  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1
 0  1  1  1  1  1  1  0  0  0  0  1  1  1  1  1
 0  1  1  0  0  1  1  0  0  0  0  1  1  0  1  1
 0  1  1  1  1  1  1  0  0  0  0  1  1  1  1  1
 0  1  1  1  1  1  1  0  0  0  0  1  1  1  1  1
 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0

 julia> isboundary_thick(A) .& (A .!= 0)
16×16 BitMatrix:
 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
 0  0  0  1  1  1  1  1  1  1  1  1  0  0  0  0
 0  0  0  1  0  0  0  1  1  0  0  1  0  0  0  0
 0  0  0  1  0  0  0  1  1  0  0  1  0  0  0  0
 0  0  0  1  0  0  0  1  1  0  0  1  0  0  0  0
 0  0  0  1  1  1  1  1  1  1  1  1  0  0  0  0
 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
 0  0  1  1  1  1  0  0  0  0  0  0  1  1  1  0
 0  0  1  0  0  1  0  0  0  0  0  0  1  0  1  0
 0  0  1  1  1  1  0  0  0  0  0  0  1  1  1  0
 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0

 julia> isboundary_thick(A) == isboundary(A; background = -1)
true

julia> isboundary_thick(A) .& (A .!= 0) == isboundary(A) # inner boundaries
true

julia> isboundary_thick(A .!= 0) .& (A .== 0)  == isboundary(A .== 0) # outer boundaries
true
 ```
"""
const isboundary_thick = isboundary_thick_dilate_erode