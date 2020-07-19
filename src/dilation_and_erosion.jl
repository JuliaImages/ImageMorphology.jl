# Erode and dilate support 3x3 regions only (and higher-dimensional generalizations).

"""
```
imgd = dilate(img, [region])
```

perform a max-filter over nearest-neighbors. The
default is 8-connectivity in 2d, 27-connectivity in 3d, etc. You can specify the
list of dimensions that you want to include in the connectivity, e.g., `region =
[1,2]` would exclude the third dimension from filtering.

# Examples

```jldoctest; setup = :(using ImageMorphology), filter = r"Array{Float64,2}|Matrix{Float64}"
julia> img = zeros(5, 5); img[3, 3] = 1.; img
5×5 Array{Float64,2}:
 0.0  0.0  0.0  0.0  0.0
 0.0  0.0  0.0  0.0  0.0
 0.0  0.0  1.0  0.0  0.0
 0.0  0.0  0.0  0.0  0.0
 0.0  0.0  0.0  0.0  0.0

julia> dilate(img)
5×5 Array{Float64,2}:
 0.0  0.0  0.0  0.0  0.0
 0.0  1.0  1.0  1.0  0.0
 0.0  1.0  1.0  1.0  0.0
 0.0  1.0  1.0  1.0  0.0
 0.0  0.0  0.0  0.0  0.0

julia> dilate(img, 1)
5×5 Array{Float64,2}:
 0.0  0.0  0.0  0.0  0.0
 0.0  0.0  1.0  0.0  0.0
 0.0  0.0  1.0  0.0  0.0
 0.0  0.0  1.0  0.0  0.0
 0.0  0.0  0.0  0.0  0.0
```
"""
dilate(img::AbstractArray, region=coords_spatial(img)) = dilate!(copy(img), region)

"""
```
imge = erode(img, [region])
```

perform a min-filter over nearest-neighbors. The
default is 8-connectivity in 2d, 27-connectivity in 3d, etc. You can specify the
list of dimensions that you want to include in the connectivity, e.g., `region =
[1,2]` would exclude the third dimension from filtering.

# Examples

```jldoctest; setup = :(using ImageMorphology), filter = r"Array{Float64,2}|Matrix{Float64}"
julia> img = zeros(5, 5); img[2:4, 2:4] .= 1.; img
5×5 Array{Float64,2}:
 0.0  0.0  0.0  0.0  0.0
 0.0  1.0  1.0  1.0  0.0
 0.0  1.0  1.0  1.0  0.0
 0.0  1.0  1.0  1.0  0.0
 0.0  0.0  0.0  0.0  0.0

julia> erode(img)
5×5 Array{Float64,2}:
 0.0  0.0  0.0  0.0  0.0
 0.0  0.0  0.0  0.0  0.0
 0.0  0.0  1.0  0.0  0.0
 0.0  0.0  0.0  0.0  0.0
 0.0  0.0  0.0  0.0  0.0

julia> erode(img, 1)
5×5 Array{Float64,2}:
 0.0  0.0  0.0  0.0  0.0
 0.0  0.0  0.0  0.0  0.0
 0.0  1.0  1.0  1.0  0.0
 0.0  0.0  0.0  0.0  0.0
 0.0  0.0  0.0  0.0  0.0
```
"""
erode(img::AbstractArray, region=coords_spatial(img)) = erode!(copy(img), region)

dilate!(maxfilt, region=coords_spatial(maxfilt)) = extremefilt!(maxfilt, max, region)
erode!(minfilt, region=coords_spatial(minfilt)) = extremefilt!(minfilt, min, region)
function extremefilt!(A::AbstractArray, select::Function, region=coords_spatial(A))
    inds = axes(A)
    for d = 1:ndims(A)
        if size(A, d) == 1 || !in(d, region)
            continue
        end
        Rpre = CartesianIndices(inds[1:d-1])
        Rpost = CartesianIndices(inds[d+1:end])
        _extremefilt!(A, select, Rpre, inds[d], Rpost)
    end
    A
end

@noinline function _extremefilt!(A, select, Rpre, inds, Rpost)
    # TODO: improve the cache efficiency
    for Ipost in Rpost, Ipre in Rpre
        # first element along dim
        i1 = first(inds)
        a2, a3 = A[Ipre, i1, Ipost], A[Ipre, i1+1, Ipost]
        A[Ipre, i1, Ipost] = select(a2, a3)
        # interior along dim
        for i = i1+2:last(inds)
            a1, a2 = a2, a3
            a3 = A[Ipre, i, Ipost]
            A[Ipre, i-1, Ipost] = select(select(a1, a2), a3)
        end
        # last element along dim
        A[Ipre, last(inds), Ipost] = select(a2, a3)
    end
    A
end

"""
`imgo = opening(img, [region])` performs the `opening` morphology operation, equivalent to `dilate(erode(img))`.
`region` allows you to control the dimensions over which this operation is performed.

# Examples

```jldoctest; setup = :(using ImageMorphology), filter = r"Array{Float64,2}|Matrix{Float64}"
julia> img = zeros(5, 5); img[1, 1] = 1.; img[3:5, 3:5] .= 1.; img
5×5 Array{Float64,2}:
 1.0  0.0  0.0  0.0  0.0
 0.0  0.0  0.0  0.0  0.0
 0.0  0.0  1.0  1.0  1.0
 0.0  0.0  1.0  1.0  1.0
 0.0  0.0  1.0  1.0  1.0

julia> opening(img)
5×5 Array{Float64,2}:
 0.0  0.0  0.0  0.0  0.0
 0.0  0.0  0.0  0.0  0.0
 0.0  0.0  1.0  1.0  1.0
 0.0  0.0  1.0  1.0  1.0
 0.0  0.0  1.0  1.0  1.0
```
"""
opening(img::AbstractArray, region=coords_spatial(img)) = opening!(copy(img), region)
opening!(img::AbstractArray, region=coords_spatial(img)) = dilate!(erode!(img, region),region)

"""
`imgc = closing(img, [region])` performs the `closing` morphology operation, equivalent to `erode(dilate(img))`.
`region` allows you to control the dimensions over which this operation is performed.

# Examples

```jldoctest; setup = :(using ImageMorphology), filter = r"Array{Float64,2}|Matrix{Float64}"
julia> img = zeros(7, 7); img[3:5, 3:5] .= 1.; img[4, 4] = 0.; img
7×7 Array{Float64,2}:
 0.0  0.0  0.0  0.0  0.0  0.0  0.0
 0.0  0.0  0.0  0.0  0.0  0.0  0.0
 0.0  0.0  1.0  1.0  1.0  0.0  0.0
 0.0  0.0  1.0  0.0  1.0  0.0  0.0
 0.0  0.0  1.0  1.0  1.0  0.0  0.0
 0.0  0.0  0.0  0.0  0.0  0.0  0.0
 0.0  0.0  0.0  0.0  0.0  0.0  0.0

julia> closing(img)
7×7 Array{Float64,2}:
 0.0  0.0  0.0  0.0  0.0  0.0  0.0
 0.0  0.0  0.0  0.0  0.0  0.0  0.0
 0.0  0.0  1.0  1.0  1.0  0.0  0.0
 0.0  0.0  1.0  1.0  1.0  0.0  0.0
 0.0  0.0  1.0  1.0  1.0  0.0  0.0
 0.0  0.0  0.0  0.0  0.0  0.0  0.0
 0.0  0.0  0.0  0.0  0.0  0.0  0.0
```
"""
closing(img::AbstractArray, region=coords_spatial(img)) = closing!(copy(img), region)
closing!(img::AbstractArray, region=coords_spatial(img)) = erode!(dilate!(img, region),region)

"""
`imgth = tophat(img, [region])` performs `top hat` of an image,
which is defined as the image minus its morphological opening.
`region` allows you to control the dimensions over which this operation is performed.

# Examples

```jldoctest; setup = :(using ImageMorphology), filter = r"Array{Float64,2}|Matrix{Float64}"
julia> img = zeros(5, 5); img[1, 1] = 1.; img[3:5, 3:5] .= 1.; img
5×5 Array{Float64,2}:
 1.0  0.0  0.0  0.0  0.0
 0.0  0.0  0.0  0.0  0.0
 0.0  0.0  1.0  1.0  1.0
 0.0  0.0  1.0  1.0  1.0
 0.0  0.0  1.0  1.0  1.0

julia> tophat(img)
5×5 Array{Float64,2}:
 1.0  0.0  0.0  0.0  0.0
 0.0  0.0  0.0  0.0  0.0
 0.0  0.0  0.0  0.0  0.0
 0.0  0.0  0.0  0.0  0.0
 0.0  0.0  0.0  0.0  0.0
```
"""
tophat(img::AbstractArray, region=coords_spatial(img)) = img - opening(img, region)

"""
`imgbh = bothat(img, [region])` performs `bottom hat` of an image,
which is defined as its morphological closing minus the original image.
`region` allows you to control the dimensions over which this operation is performed.

# Examples

```jldoctest; setup = :(using ImageMorphology), filter = r"Array{Float64,2}|Matrix{Float64}"
julia> img = zeros(7, 7); img[3:5, 3:5] .= 1.; img[4, 4] = 0.; img
7×7 Array{Float64,2}:
 0.0  0.0  0.0  0.0  0.0  0.0  0.0
 0.0  0.0  0.0  0.0  0.0  0.0  0.0
 0.0  0.0  1.0  1.0  1.0  0.0  0.0
 0.0  0.0  1.0  0.0  1.0  0.0  0.0
 0.0  0.0  1.0  1.0  1.0  0.0  0.0
 0.0  0.0  0.0  0.0  0.0  0.0  0.0
 0.0  0.0  0.0  0.0  0.0  0.0  0.0

julia> bothat(img)
7×7 Array{Float64,2}:
 0.0  0.0  0.0  0.0  0.0  0.0  0.0
 0.0  0.0  0.0  0.0  0.0  0.0  0.0
 0.0  0.0  0.0  0.0  0.0  0.0  0.0
 0.0  0.0  0.0  1.0  0.0  0.0  0.0
 0.0  0.0  0.0  0.0  0.0  0.0  0.0
 0.0  0.0  0.0  0.0  0.0  0.0  0.0
 0.0  0.0  0.0  0.0  0.0  0.0  0.0
```
"""
bothat(img::AbstractArray, region=coords_spatial(img)) = closing(img, region) - img

"""
`imgmg = morphogradient(img, [region])` returns morphological gradient of the image,
which is the difference between the dilation and the erosion of a given image.
`region` allows you to control the dimensions over which this operation is performed.

# Examples

```jldoctest; setup = :(using ImageMorphology), filter = r"Array{Float64,2}|Matrix{Float64}"
julia> img = zeros(7, 7); img[3:5, 3:5] .= 1.; img
7×7 Array{Float64,2}:
 0.0  0.0  0.0  0.0  0.0  0.0  0.0
 0.0  0.0  0.0  0.0  0.0  0.0  0.0
 0.0  0.0  1.0  1.0  1.0  0.0  0.0
 0.0  0.0  1.0  1.0  1.0  0.0  0.0
 0.0  0.0  1.0  1.0  1.0  0.0  0.0
 0.0  0.0  0.0  0.0  0.0  0.0  0.0
 0.0  0.0  0.0  0.0  0.0  0.0  0.0

julia> morphogradient(img)
7×7 Array{Float64,2}:
 0.0  0.0  0.0  0.0  0.0  0.0  0.0
 0.0  1.0  1.0  1.0  1.0  1.0  0.0
 0.0  1.0  1.0  1.0  1.0  1.0  0.0
 0.0  1.0  1.0  0.0  1.0  1.0  0.0
 0.0  1.0  1.0  1.0  1.0  1.0  0.0
 0.0  1.0  1.0  1.0  1.0  1.0  0.0
 0.0  0.0  0.0  0.0  0.0  0.0  0.0
```
"""
morphogradient(img::AbstractArray, region=coords_spatial(img)) = dilate(img, region) - erode(img, region)

"""
`imgml = morpholaplace(img, [region])` performs `Morphological Laplacian` of an image,
which is defined as the arithmetic difference between the internal and the external gradient.
`region` allows you to control the dimensions over which this operation is performed.

# Examples

```jldoctest; setup = :(using ImageMorphology), filter = r"Array{Float64,2}|Matrix{Float64}"
julia> img = zeros(7, 7); img[3:5, 3:5] .= 1.; img[4, 4] = 0.; img
7×7 Array{Float64,2}:
 0.0  0.0  0.0  0.0  0.0  0.0  0.0
 0.0  0.0  0.0  0.0  0.0  0.0  0.0
 0.0  0.0  1.0  1.0  1.0  0.0  0.0
 0.0  0.0  1.0  0.0  1.0  0.0  0.0
 0.0  0.0  1.0  1.0  1.0  0.0  0.0
 0.0  0.0  0.0  0.0  0.0  0.0  0.0
 0.0  0.0  0.0  0.0  0.0  0.0  0.0

julia> morpholaplace(img)
7×7 Array{Float64,2}:
 0.0  0.0   0.0   0.0   0.0  0.0  0.0
 0.0  1.0   1.0   1.0   1.0  1.0  0.0
 0.0  1.0  -1.0  -1.0  -1.0  1.0  0.0
 0.0  1.0  -1.0   1.0  -1.0  1.0  0.0
 0.0  1.0  -1.0  -1.0  -1.0  1.0  0.0
 0.0  1.0   1.0   1.0   1.0  1.0  0.0
 0.0  0.0   0.0   0.0   0.0  0.0  0.0
```
"""
morpholaplace(img::AbstractArray, region=coords_spatial(img)) = dilate(img, region) + erode(img, region) - 2img
