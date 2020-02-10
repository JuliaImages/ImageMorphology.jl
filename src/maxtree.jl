"""
Max-tree morphological representation of an image.

# Details
Let's consider a *thresholding* operation,
```
    mask = [val ≥ threshold for val in image]
```
One can identify the connected components (the sets of neighboring true values)
in `mask`. When *image thresholding* is sequentially applied for all possible
thresholds, it generates a collection of connected components that could be
organized into a hierarchical structure called *component tree*. Consider 1D
"image" with values 1, 2 and 3:
```
       2233233312223322
```
The connected components would be
```
    1: AAAAAAAAAAAAAAAA
    2: BBBBBBBB.CCCCCCC
    3: ..DD.EEE....FF..
```
Here, the letters are the labels of the resulting connected components,
and `.` specifies that the pixel value is below the threshold.
In this example, the corresponding *component tree* is:
```
      A
     ⭩ ⭨
    B   C
   ⭩ ⭨   ⭨
  D   E   F
```

A *max-tree* is an efficient representation of the *component tree*.
A connected component ``C`` at threshold level ``t`` is represented by the
single *reference pixel* ``r`` from this level (`image[r] == t`), which is the
parent to all other pixels of ``C`` and also to the *reference pixels* of the
connected components at higher thresholds, which are the children of ``C``.
In our example, the reference pixels (denoted by the letter of the corresponding
component) would be:
```
    1: ........A.......
    2: B........C......
    3: ..D..E......F...
```
I.e.

| Comp | Ref.Pixel |
| ---- | ---------:|
| *A* |  9 |
| *B* |  1 |
| *C* | 10 |
| *D* |  3 |
| *E* |  6 |
| *F* | 13 |

So the whole max-tree could be encoded as a vector of indices of parent pixels:
```
9  1  1  3  1  1  6  6  9  9 10 10 10 13 10 10
```

The *max-tree* is the basis for many morphological operators,
namely connected operators. Unlike morphological openings and closings, these
operators do not require a fixed structuring element, but rather act with a
flexible structuring element that meets a certain criterion.

# See also
[`area_opening`](@ref), [`area_closing`](@ref),
[`diameter_opening`](@ref), [`diameter_closing`](@ref).

# References
1. Salembier, P., Oliveras, A., & Garrido, L. (1998). *Antiextensive Connected Operators for Image and Sequence Processing*. IEEE Transactions on Image Processing, 7(4), 555-570.
   > https://doi.org/10.1109/83.663500
2. Berger, C., Geraud, T., Levillain, R., Widynski, N., Baillard, A., Bertin, E. (2007). *Effective Component Tree Computation with Application to Pattern Recognition in Astronomical Imaging*. In International Conference on Image Processing (ICIP), 41-44.
   > https://doi.org/10.1109/ICIP.2007.4379949
3. Najman, L., & Couprie, M. (2006). *Building the component tree in quasi-linear time*. IEEE Transactions on Image Processing, 15(11), 3531-3539.
   > https://doi.org/10.1109/TIP.2006.877518
4. Carlinet, E., & Geraud, T. (2014). *A Comparative Review of Component Tree Computation Algorithms*. IEEE Transactions on Image Processing, 23(9), 3885-3895.
   > https://doi.org/10.1109/TIP.2014.2336551
"""
struct MaxTree{N}
    """
    If `false`, the tree is constructed from the smallest to the largest values,
    the tree root being the smallest (darkest) pixel.
    Otherwise, if `true`, the direction is from large to small values.
    """
    rev::Bool

    """
    Array of the same shape as the source image.
    Each value is the linear index of the parent pixel in the max-tree.
    """
    parentindices::Array{Int,N}

    """
    The order of the elements in the tree, from root to leaves.
    Each element corresponds to the linear index of the pixel in the image.
    A parent of the pixel always comes before the pixel itself, i.e.
    `parentindices[traverse[1]] == traverse[1]` (the root pixel) and
    `parentindices[traverse[k]] ∈ traverse[1:k-1]` for all `k >= 2`.
    """
    traverse::Vector{Int}

    # create uninitialized MaxTree for image
    MaxTree{N}(image::GenericGrayImage{<:Any, N}, rev::Bool) where N =
        new{N}(rev, Array{Int,N}(undef, size(image)),
               Vector{Int}(undef, length(image)))
end

Base.ndims(::Type{<:MaxTree{N}}) where N = N
Base.ndims(maxtree::MaxTree) = ndims(typeof(maxtree))
Base.length(maxtree::MaxTree) = length(maxtree.parentindices)
Base.size(maxtree::MaxTree) = size(maxtree.parentindices)
Base.:(==)(a::MaxTree, b::MaxTree) =
    (ndims(a) == ndims(b)) && (a.rev == b.rev) &&
    (a.parentindices == b.parentindices) &&
    (a.traverse == b.traverse)

"""
    root_index(maxtree::MaxTree) -> Int

Linear index of the root pixel.
"""
root_index(maxtree::MaxTree) = maxtree.traverse[1]

# Gets the path to the root of the max-tree (defined by `ancestors`) starting from
# `node`. The indices of the nodes along the path are stored into `path`.
# Returns the index of the root
@inline function rootpath!(path::Vector{Int},
                           ancestors::AbstractArray{<:Integer},
                           node::Integer)
    empty!(path)
    p = node
    @label next
    push!(path, p)
    @inbounds q = ancestors[p]
    if p != q
        p = q
        @goto next
    end
    return p
end

# Corrects `maxtree.parentindices`, so that the parent of every pixel
# is a canonical node. I.e. for each connected component at a level l,
# all pixels point to the same representative, which in turn points to the
# representative pixel at the next level.
#
# A canonical node ``p`` is either a root of the max-tree, or
# ``image(parent(p)) != image(p)``
function canonize!(maxtree::MaxTree, image::GenericGrayImage)
    @assert size(maxtree) == size(image)
    parentindices = maxtree.parentindices
    @inbounds for p in maxtree.traverse
        q = parentindices[p]
        if image[q] == image[parentindices[q]]
            parentindices[p] = parentindices[q]
        end
    end
    return maxtree
end

# Checks whether the `pixel` + `offset` is still a valid pixel
# contained inside the image area defined by `axes`.
isvalid_offset(offset::CartesianIndex{N}, pixel::CartesianIndex{N},
               axes::NTuple{N}) where N =
    Base.checkbounds_indices(Bool, axes, (pixel + offset,))

# generates a vector of offsets to the CartesianIndex that defines
# the neighborhood of the pixel in an N-dimensional image.
# `connectivity` is the maximal number of orthogonal steps
# (+1/-1 to a single dimension) that are required to reach a neighbor.
# Note that the corresponding offsets of linear indices would be sorted in ascending order
# ensuring efficient memory access pattern.
function neighbor_cartesian_offsets(::Type{CartesianIndex{N}}, connectivity::Integer) where N
    (connectivity >= 1) || throw(ArgumentError("connectivity should be positive integer"))
    ci1 = oneunit(CartesianIndex{N})
    return [ci for ci in -ci1:ci1 if 1 <= sum(abs, Tuple(ci)) <= connectivity]
end

# convert offsets of cartesian indices into linear offsets of the given array
linear_offsets(offsets::AbstractVector{<:CartesianIndex},
               arr::AbstractArray) =
    dot.(Tuple.(offsets), Ref(strides(arr)))

"""
    rebuild!(maxtree::MaxTree, image::GenericGrayImage,
             neighbors::AbstractVector{CartesianIndex}) -> maxtree

Rebuilds the `maxtree` for the `image` using `neighbors` as the pixel
connectivity specification.

# Details
The pixels in the connected components generated by the method should be
connected to each other by a path through neighboring pixels.
The pixels ``p_1`` and ``p_2`` are neighbors, if `neighbors` array contains
``d``, such that ``p_2 = p_1 + d``.

# See also
[`MaxTree`](@ref)
"""
function rebuild!(maxtree::MaxTree{N},
                  image::GenericGrayImage{<:Any, N},
                  #=mask::AbstractArray{Bool, N},=#
                  neighbors::AbstractVector{CartesianIndex{N}}) where N
    check_maxtree(maxtree, image, rev=maxtree.rev)
    # pixels need to be sorted according to their gray level.
    sortperm!(maxtree.traverse, vec(image), rev=maxtree.rev)

    neighbor_offsets = linear_offsets(neighbors, maxtree.parentindices)
    p_indices = fill!(maxtree.parentindices, 0) |> vec # 0=uninitialized parent
    p_cis = CartesianIndices(maxtree.parentindices)
    p_axes = axes(maxtree.parentindices)
    r_indices = fill!(similar(p_indices), 0) # temporary array containing current roots
    rootpath = Vector{Int}() # temp array to hold the path to the current root within ancestors

    # traverse the array in reversed order (from highest value to lowest value)
    @inbounds for p in Iterators.Reverse(maxtree.traverse)
        p_ci = p_cis[p]
        p_indices[p] = r_indices[p] = p # it's a new root

        for i in eachindex(neighbors)
            # check that p is in the mask and that it's neighbor is a valid image pixel
            (#=mask[p] && =#isvalid_offset(neighbors[i], p_ci, p_axes)) || continue
            index = p + neighbor_offsets[i] # linear index of the neighbor
            (p_indices[index] == 0) && continue # ignore neighbor without parent (this means image[index] < image[p])
            root = rootpath!(rootpath, r_indices, index) # neighbor's root
            if (root != p) || (length(rootpath) > 1) # attach neighbor's root to the p
                p_indices[root] = p
                r_indices[rootpath] .= p # also attach the ancestral nodes of neighbor to p
            end
        end
    end
    return canonize!(maxtree, image)
end

"""
    MaxTree(image::GenericGrayImage; connectivity=1, rev=false) -> MaxTree

Constructs the *max-tree* of the `image`.

# Arguments
- `connectivity::Integer=1`: defines the pixel neighborhood used to construct
  the connected components. The value is the maximum number of orthogonal steps
  to reach a neighbor. In 2D, it is 1 for a 4-neighborhood and 2 for a
  8-neighborhood. See [`rebuild!](@ref).
- `rev::Bool=false`: if `false`, the max-tree is traversed from the darkest
  (the root node) to the brightest, otherwise it's traversed from the brightest
  (the root) to the darkest.

# Examples
We create a small sample image (Figure 1 from [4]) and build the max-tree.

```jldoctest
julia> image = [15 13 16; 12 12 10; 16 12 14]
3×3 Array{Int64,2}:
 15  13  16
 12  12  10
 16  12  14

julia> mtree = MaxTree(image, connectivity=2)
MaxTree{2}(false, [4 2 4; 8 2 8; 2 2 2], [8, 2, 5, 6, 4, 9, 1, 3, 7])
```
"""
function MaxTree(image::GenericGrayImage;
                 connectivity::Integer=1, rev::Bool=false)
    N = ndims(image)
    # User defined masks are not allowed, as there might be more than one
    # connected component in the mask (and therefore not a single tree that
    # represents the image). Mask here is an image that is 0 on the border
    # and 1 everywhere else
    #=
    mask = np.ones(image.shape)
    for k in range(len(image.shape)):
        np.moveaxis(mask, k, 0)[0] = 0
        np.moveaxis(mask, k, 0)[-1] = 0
    =#
    neighbors = neighbor_cartesian_offsets(CartesianIndex{N}, connectivity)
    return rebuild!(MaxTree{N}(image, rev), #=mask,=# image, neighbors)
end

"""
    areas(maxtree::MaxTree) -> Array{Int}

Computes the areas of all `maxtree` components.

# Returns
The array of the same shape as the original image. The `i`-th element is
the area (in pixels) of the component that is represented by the reference
pixel with index `i`.

# See also
[`diameters`](@ref), [`area_opening`](@ref), [`area_closing`](@ref).
"""
function areas(maxtree::MaxTree)
    areas = fill(1, size(maxtree)) # start with 1-pixel areas
    parentindices = maxtree.parentindices
    @inbounds for p in Iterators.Reverse(maxtree.traverse)
        q = parentindices[p]
        (q != p) && (areas[q] += areas[p])
    end
    return areas
end

"""
    boundingboxes(maxtree::MaxTree) -> Array{NTuple{2, CartesianIndex}}

Computes the minimal bounding boxes of all `maxtree` components.

# Returns
The array of the same shape as the original image. The `i`-th element is
the tuple of the minimal and maximal cartesian indices for the bounding box
of the component that is represented by the reference pixel with index `i`.

# See also
[`diameters`](@ref).
"""
function boundingboxes(maxtree::MaxTree{N}) where N
    # initialize bboxes
    bboxes = [(ci, ci) for ci in CartesianIndices(maxtree.parentindices)]
    @inbounds for p in Iterators.Reverse(maxtree.traverse)
        q = maxtree.parentindices[p]
        bbp, bbq = bboxes[p], bboxes[q]
        bboxes[q] = (min(bbq[1], bbp[1]), max(bbq[2], bbp[2]))
    end
    return bboxes
end

"""
    diameters(maxtree::MaxTree) -> Array{Int}

Computes the "diameters" of all `maxtree` components.

"Diameter" of the *max-tree* connected component is the length of
the widest side of the component's bounding box.

# Returns
The array of the same shape as the original image. The `i`-th element is the
"diameter" of the component that is represented by the reference pixel with
index `i`.

# See also
[`boundingboxes`](@ref), [`areas`](@ref),
[`diameter_opening`](@ref), [`diameter_closing`](@ref).
"""
diameters(maxtree::MaxTree) =
    [maximum(Tuple(bbox[2] - bbox[1])) + 1 for bbox in boundingboxes(maxtree)]

"""
    filter_components!(output::GenericGrayImage, image::GenericGrayImage,
                       maxtree::MaxTree, attrs::AbstractVector,
                       min_attr, all_below_min) -> output

Filters the connected components of the `image` and stores the result in `output`.

The ``output`` is the copy of the ``image`` exluding the connected components,
whose attribute value is below `min_attr`. That is, the pixels of the exluded
component are reset to the value of the reference pixel of its first valid
ancestor (the connected component with the attribute value greater or equal
to `min_attr`).

# Arguments
- `maxtree::MaxTree`: pre-built max-tree representation of the `image`
- `attrs::AbstractVector`: `attrs[i]` is the attribute value for the ``i``-th
   component of the tree (``i`` being the linear index of its *reference pixel*)
- `all_below_min`: the value to fill the `output` if all attributes of all
  components (including the root one) are below `min_attr`

# Details
This function is the basis for [`area_opening`](@ref), [`diameter_opening`](@ref)
and similar transformations.
E.g. for [`area_opening`](@ref) the attribute is the area of the components.
In this case, the max-tree components of the `output` have area no smaller
than `min_attr` pixels.

The method assumes that the attribute values are monotone with respect to the
components hieararchy, i.e. ``attrs[i] <= attrs[maxtree.parentindices[i]]`` for
each `i`.
"""
function filter_components!(output::GenericGrayImage, image::GenericGrayImage,
                            maxtree::MaxTree,
                            attrs::AbstractArray, min_attr,
                            all_below_min = zero(eltype(output)))
    # should have been already checked by higher-level functions
    @assert size(output) == size(image)
    @assert size(image) == size(maxtree)
    @assert length(attrs) == length(maxtree)

    p_root = root_index(maxtree)
    # if p_root is below min_attr, then all others are as well
    (attrs[p_root] < min_attr) && return fill!(output, all_below_min)
    output[p_root] = image[p_root]
    parentindices = maxtree.parentindices

    for p in maxtree.traverse
        q = parentindices[p]
        # this means p is not canonical
        # in other words, it has a parent that has the
        # same image value.
        if image[p] == image[q]
            output[p] = output[q]
            continue
        end
        # if attribute is below the threshold, stop at the lower parent,
        # otherwise continue to reconstruct the image
        output[p] = attrs[p] < min_attr ? output[q] : image[p]
    end

    return output
end

function check_output_image(output::GenericGrayImage, image::GenericGrayImage)
    (size(output) == size(image)) ||
        throw(DimensionMismatch("The sizes of the output and the input image do not match"))
end

# checks if the provided maxtree is compatible with the given options
# or builds the new maxtree if none was given
function check_maxtree(maxtree::Union{MaxTree, Nothing},
                       image::GenericGrayImage;
                       connectivity::Integer=1,
                       rev::Bool=false)
    (maxtree === nothing) && return MaxTree(image, connectivity=connectivity, rev=rev)
    (size(maxtree) == size(image)) ||
        throw(DimensionMismatch("The sizes of the max-tree and the input image do not match"))
    (maxtree.rev == rev) ||
        throw(ArgumentError("The traversal order of the given max-tree is different from the requested one"))
    return maxtree
end

"""
    area_opening!(output, image;
                  min_area=64, connectivity=1, maxtree=nothing) -> output

Performs in-place *area opening* of the `image` and stores the result in `output`.
See [`area_opening`](@ref) for the detailed description of the method.
"""
function area_opening!(output::GenericGrayImage, image::GenericGrayImage;
                       min_area::Number=64, connectivity::Integer=1,
                       maxtree::Union{MaxTree, Nothing}=nothing)
    check_output_image(output, image)
    _maxtree = check_maxtree(maxtree, image, connectivity=connectivity, rev=false)
    return filter_components!(output, image, _maxtree, areas(_maxtree), min_area)
end

"""
    area_opening(image; min_area=64, connectivity=1, maxtree=nothing) -> Array

Performs an *area opening* of the `image`.

*Area opening* replaces all bright components of an image that
have a surface smaller than `min_area` with the darker value taken from
their first ancestral component (in *max-tree* representation of `image`)
that has the area no smaller than `min_area`.

# Details
Area opening is similar to morphological opening (see [`opening`](@ref)),
but instead of using a fixed structuring element (e.g. disk) it employs
small (less than `min_area`) components of the *max-tree*. Consequently,
the `area_opening` with `min_area = 1` is the identity transformation.

In the binary case, area opening is equivalent to `remove_small_objects`;
this operator is thus extended to gray-level images.

# Arguments
- `image::GenericGrayImage`: the ``N``-dimensional input image
- `min_area::Number=64`: the smallest size (in pixels) of the image component
  to keep intact
- `connectivity::Integer=1`: the neighborhood connectivity. The maximum number
  of orthogonal steps to reach a neighbor of the pixel.
  In 2D, it is 1 for a 4-neighborhood and 2 for a 8-neighborhood.
- `maxtree::MaxTree=nothing`: optional pre-built *max-tree*.
  Note that if `maxtree` is specified, the `connectivity` option is ignored.

# Returns
An array of the same type and shape as the `image`.

# See also
[`area_opening!`](@ref), [`area_closing`](@ref),
[`diameter_opening`](@ref),
[`MaxTree`](@ref), [`opening`](@ref)

# References
1. Vincent, L. (1993). *Grayscale area openings and closings, their efficient implementation and applications*, Proc. of EURASIP Workshop on Mathematical Morphology and its Applications to Signal Processing, Barcelona, Spain, 22-27
2. Soille, P. (2003). Chapter 6 *Geodesic Metrics* of *Morphological Image Analysis: Principles and Applications*, 2nd edition, Springer.
   > https://doi.org/10.1007/978-3-662-05088-0
3. Salembier, P., Oliveras, A., & Garrido, L. (1998). *Antiextensive Connected Operators for Image and Sequence Processing*. IEEE Transactions on Image Processing, 7(4), 555-570.
   > https://doi.org/10.1109/83.663500
4. Najman, L., & Couprie, M. (2006). *Building the component tree in quasi-linear time*. IEEE Transactions on Image Processing, 15(11), 3531-3539.
   > https://doi.org/10.1109/TIP.2006.877518
5. Carlinet, E., & Geraud, T. (2014). *A Comparative Review of Component Tree Computation Algorithms*. IEEE Transactions on Image Processing, 23(9), 3885-3895.
   > https://doi.org/10.1109/TIP.2014.2336551

# Examples
Creating a test image `f` (quadratic function with a maximum in the center and
4 additional local maxima):
```jldoctest
julia> w = 12;
julia> f = [20 - 0.2*((x - w/2)^2 + (y-w/2)^2) for x in 0:w, y in 0:w];
julia> f[3:4, 2:6] .= 40; f[3:5, 10:12] .= 60; f[10:12, 3:5] .= 80;
julia> f[10:11, 10:12] .= 100; f[11, 11] = 100;
```
Area opening of `f`:
```jldoctest
julia> f_aopen = area_opening(f, min_area=8, connectivity=1)
```
The peaks with a surface smaller than 8 are removed.
"""
area_opening(image::GenericGrayImage; kwargs...) =
    area_opening!(similar(image), image; kwargs...)

"""
    diameter_opening!(output, image; min_diameter=8,
                      connectivity=1, maxtree=nothing) -> output

Performs in-place *diameter opening* of the `image` and stores the result in `output`.
See [`diameter_opening`](@ref) for the detailed description of the method.
"""
function diameter_opening!(output::GenericGrayImage, image::GenericGrayImage;
                           maxtree::Union{MaxTree, Nothing} = nothing,
                           min_diameter=8, connectivity=1)
    check_output_image(output, image)
    _maxtree = check_maxtree(maxtree, image, connectivity=connectivity, rev=false)
    return filter_components!(output, image, _maxtree, diameters(_maxtree), min_diameter)
end

"""
    diameter_opening(image; min_diameter=8, connectivity=1,
                     maxtree=nothing) -> Array

Performs a *diameter opening* of the `image`.

*Diameter opening* replaces all bright structures of an image that have
the diameter (the widest dimension of their bounding box) smaller than
`min_diameter` with the darker value taken from their first ancestral component
(in *max-tree* representation of `image`) that has the diameter no smaller than
`min_diameter`.

The operator is also called *Bounding Box Opening*. In practice,
the result is similar to a *morphological opening*, but long and thin
structures are not removed.

# Arguments
- `image::GenericGrayImage`: the ``N``-dimensional input image
- `min_diameter::Number=8`: the minimal length (in pixels) of the widest
  dimension of the bounding box of the image component to keep intact
- `connectivity::Integer=1`: the neighborhood connectivity. The maximum number
  of orthogonal steps to reach a neighbor of the pixel. In 2D, it is 1 for
  a 4-neighborhood and 2 for a 8-neighborhood.
- `maxtree::MaxTree=nothing`: optional pre-built *max-tree*. Note that,
  when `maxtree` is specified, the `connectivity` option is ignored.

# Returns
An array of the same type and shape as the `image`.

# See also
[`diameter_opening!`](@ref), [`diameter_closing`](@ref),
[`area_opening`](@ref),
[`MaxTree`](@ref), [`opening`](@ref)

# References
1. Walter, T., & Klein, J.-C. (2002). *Automatic Detection of Microaneurysms in Color Fundus Images of the Human Retina by Means of the Bounding Box Closing*. In A. Colosimo, P. Sirabella, A. Giuliani (Eds.), *Medical Data Analysis. Lecture Notes in Computer Science*, vol 2526, 210-220. Springer Berlin Heidelberg.
   > https://doi.org/10.1007/3-540-36104-9_23
2. Carlinet, E., & Geraud, T. (2014). *A Comparative Review of Component Tree Computation Algorithms*. IEEE Transactions on Image Processing, 23(9), 3885-3895.
   > https://doi.org/10.1109/TIP.2014.2336551

# Examples
Creating a test image `f` (quadratic function with a maximum in the center and
4 additional local maxima):
```jldoctest
julia> w = 12;
julia> f = [20 - 0.2*((x - w/2)^2 + (y-w/2)^2) for x in 0:w, y in 0:w];
julia> f[3:4, 2:6] .= 40; f[3:5, 10:12] .= 60; f[10:12, 3:5] .= 80;
julia> f[10:11, 10:12] .= 100; f[11, 11] = 100;
```
Diameter opening of `f`:
```jldoctest
julia> f_dopen = diameter_opening(f, min_diameter=3, connectivity=1)
```
The peaks with a maximal diameter of 2 or less are removed.
For the remaining peaks the widest side of the bounding box is at least 3.
"""
diameter_opening(image::GenericGrayImage; kwargs...) =
    diameter_opening!(similar(image), image; kwargs...)

"""
    area_closing!(output, image; min_area=64, connectivity=1,
                  maxtree=nothing) -> output

Performs in-place *area closing* of the `image` and stores the result in `output`.
See [`area_closing`](@ref) for the detailed description of the method.
"""
function area_closing!(output::GenericGrayImage, image::GenericGrayImage;
                       min_area::Number=64, connectivity::Integer=1,
                       maxtree::Union{MaxTree, Nothing}=nothing)
    check_output_image(output, image)
    _maxtree = check_maxtree(maxtree, image, connectivity=connectivity, rev=true)
    return filter_components!(output, image, _maxtree, areas(_maxtree), min_area)
end

"""
    area_closing(image; min_area=64, connectivity=1, maxtree=nothing) -> Array

Performs an *area closing* of the `image`.

*Area closing* replaces all dark components of an image that
have a surface smaller than `min_area` with the brighter value taken from
their first ancestral component (in *max-tree* representation of `image`)
that has the area no smaller than `min_area`.

# Details
*Area closing* is the dual operation to *area opening* (see [`area_opening`](@ref)).
It is similar to morphological closings (see [`closing`](@ref)),
but instead of using a fixed structuring element (e.g. disk) it employs
small (less than `min_area`) components of the *max-tree*. Consequently,
the `area_closing` with `min_area = 1` is the identity transformation.

In the binary case, area closing is equivalent to
`remove_small_holes`; this operator is thus extended to gray-level images.

# Arguments
- `image::GenericGrayImage`: the ``N``-dimensional input image
- `min_area::Number=64`: the smallest size (in pixels) of the image component
  to keep intact
- `connectivity::Integer=1`: the neighborhood connectivity. The maximum number
  of orthogonal steps to reach a neighbor of the pixel. In 2D, it is 1 for
  a 4-neighborhood and 2 for a 8-neighborhood.
- `maxtree::MaxTree=nothing`: optional pre-built *max-tree*. Note that,
  when `maxtree` is specified, the `connectivity` option is ignored.

# Returns
An array of the same type and shape as the `image`.

# See also
[`area_closing!`](@ref), [`area_opening`](@ref),
[`diameter_closing`](@ref),
[`MaxTree`](@ref), [`closing`](@ref)

# References
1. Vincent, L. (1993). *Grayscale area openings and closings, their efficient implementation and applications*, Proc. of EURASIP Workshop on Mathematical Morphology and its Applications to Signal Processing, Barcelona, Spain, 22-27
2. Soille, P. (2003). Chapter 6 *Geodesic Metrics* of *Morphological Image Analysis: Principles and Applications*, 2nd edition, Springer.
   > https://doi.org/10.1007/978-3-662-05088-0
3. Salembier, P., Oliveras, A., & Garrido, L. (1998). *Antiextensive Connected Operators for Image and Sequence Processing*. IEEE Transactions on Image Processing, 7(4), 555-570.
   > https://doi.org/10.1109/83.663500
4. Najman, L., & Couprie, M. (2006). *Building the component tree in quasi-linear time*. IEEE Transactions on Image Processing, 15(11), 3531-3539.
   > https://doi.org/10.1109/TIP.2006.877518
5. Carlinet, E., & Geraud, T. (2014). *A Comparative Review of Component Tree Computation Algorithms*. IEEE Transactions on Image Processing, 23(9), 3885-3895.
   > https://doi.org/10.1109/TIP.2014.2336551

# Examples
Creating a test image `f` (quadratic function with a minimum in the center and
4 additional local minima):
```jldoctest
julia> w = 12;
julia> f = [180 + 0.2*((x - w/2)^2 + (y-w/2)^2) for x in 0:w, y in 0:w];
julia> f[3:4, 2:6] .= 40; f[3:5, 10:12] .= 60; f[10:12, 3:5] .= 80;
julia> f[10:11, 10:12] .= 100; f[11, 11] = 100;
```
Area closing of `f`:
```jldoctest
julia> f_aclose = area_closing(f, min_area=8, connectivity=1)
```
All small minima are removed, and the remaining minima have at least
a size of 8.
"""
area_closing(image::GenericGrayImage; kwargs...) =
    area_closing!(similar(image), image; kwargs...)

"""
    diameter_closing!(output, image; min_diameter=8, connectivity=1,
                      maxtree=nothing) -> output

Performs in-place *diameter closing* of the `image` and stores the result in `output`.
See [`diameter_closing`](@ref) for the detailed description of the method.
"""
function diameter_closing!(output::GenericGrayImage, image::GenericGrayImage;
                           min_diameter::Number=8, connectivity::Integer=1,
                           maxtree::Union{MaxTree, Nothing} = nothing)
    check_output_image(output, image)
    _maxtree = check_maxtree(maxtree, image, connectivity=connectivity, rev=true)
    return filter_components!(output, image, _maxtree, diameters(_maxtree), min_diameter)
end

"""
    diameter_closing(image; min_diameter=8, connectivity=1,
                     maxtree=nothing) -> Array

Performs a *diameter closing* of the `image`.

*Diameter closing* replaces all dark structures of an image that have
the diameter (the widest dimension of their bounding box) smaller than
`min_diameter` with the brighter value taken from their first ancestral component
(in *max-tree* representation of `image`) that has the diameter no smaller than
`min_diameter`.

# Arguments
- `image::GenericGrayImage`: the ``N``-dimensional input image
- `min_diameter::Number=8`: the minimal length (in pixels) of the widest
  dimension of the bounding box of the image component to keep intact
- `connectivity::Integer=1`: the neighborhood connectivity. The maximum number
  of orthogonal steps to reach a neighbor of the pixel. In 2D, it is 1 for
  a 4-neighborhood and 2 for a 8-neighborhood.
- `maxtree::MaxTree=nothing`: optional pre-built *max-tree*. Note that,
  when `maxtree` is specified, the `connectivity` option is ignored.

# Returns
An array of the same type and shape as the `image`.

# See also
[`diameter_closing!`](@ref), [`diameter_opening`](@ref),
[`area_closing`](@ref),
[`MaxTree`](@ref), [`closing`](@ref)

# References
1. Walter, T., & Klein, J.-C. (2002). *Automatic Detection of Microaneurysms in Color Fundus Images of the Human Retina by Means of the Bounding Box Closing*. In A. Colosimo, P. Sirabella, A. Giuliani (Eds.), *Medical Data Analysis. Lecture Notes in Computer Science*, vol 2526, 210-220. Springer Berlin Heidelberg.
   > https://doi.org/10.1007/3-540-36104-9_23
2. Carlinet, E., & Geraud, T. (2014). *A Comparative Review of Component Tree Computation Algorithms*. IEEE Transactions on Image Processing, 23(9), 3885-3895.
   > https://doi.org/10.1109/TIP.2014.2336551

# Examples
Creating a test image `f` (quadratic function with a minimum in the center and
4 additional local minima):
```jldoctest
julia> w = 12;
julia> f = [180 + 0.2*((x - w/2)^2 + (y-w/2)^2) for x in 0:w, y in 0:w];
julia> f[3:4, 2:6] .= 40; f[3:5, 10:12] .= 60; f[10:12, 3:5] .= 80;
julia> f[10:11, 10:12] .= 100; f[11, 11] = 100;
```
Area closing of `f`:
```jldoctest
julia> f_dclose = diameter_closing(f, min_diameter=3, connectivity=1)
```
All small minima with a diameter of 2 or less are removed.
For the remaining minima the widest bounding box side is at least 3.
"""
diameter_closing(image::GenericGrayImage; kwargs...) =
    diameter_closing!(similar(image), image; kwargs...)

# Calculates the local maxima or minima of the image using the max-tree
# representation (maxima if maxtree.rev=false, minima otherwise).
# It's not the most effecient method for extrema calculation, but could be
# useful if the max-tree is already calculated.
# Each minima/maxima region is labeled with the unique id (1 is the global
# maximum/minimum).
function local_extrema!(output::GenericGrayImage,
                        image::GenericGrayImage,
                        maxtree::MaxTree)
    (size(image) == size(maxtree)) || throw(DimensionMismatch())
    (size(output) == size(image)) || throw(DimensionMismatch())

    next_label = oneunit(eltype(output))
    fill!(output, next_label) # initialize output (just has to be non-zero)
    @inbounds for p in Iterators.reverse(maxtree.traverse)
        q = maxtree.parentindices[p]
        # if p is canonical (parent has a different value)
        if image[p] != image[q]
            output[q] = 0
            # if output[p] was the parent of some other canonical
            # pixel, it has been set to zero. Only the leaves
            # (local maxima) are thus > 0.
            if output[p] > 0
                output[p] = next_label
                next_label += 1
            end
        end
    end

    @inbounds for p in Iterators.reverse(maxtree.traverse)
        q = maxtree.parentindices[p]
        # if p is not canonical (parent has the same value)
        if image[p] == image[q]
            # in this case we propagate the value
            output[p] = output[q]
        end
    end

    return output
end

"""
    local_maxima!(output, image; connectivity=1, maxtree=nothing) -> output

Detects the local maxima of `image` and stores the result in `output`.
See [`local_maxima`](@ref) for the detailed description of the method.
"""
function local_maxima!(output::GenericGrayImage, image::GenericGrayImage;
                       connectivity::Integer=1,
                       maxtree::Union{MaxTree, Nothing} = nothing)
    check_output_image(output, image)
    _maxtree = check_maxtree(maxtree, image, connectivity=connectivity, rev=false)
    return local_extrema!(output, image, _maxtree)
end

"""
    local_maxima(image::GenericGrayImage; connectivity=1, maxtree=nothing) -> Array

Determines and labels all *local maxima* of the `image`.

# Details
The *local maximum* is defined as the connected set of pixels that have the
same value, which is greater than the values of all pixels in direct
neighborhood of the set.

Technically, the implementation is based on the *max-tree* representation
of an image. It's beneficial if the max-tree is already computed,
otherwise [`Images.findlocalmaxima`](@ref) would be more efficient.

# Arguments
- `image::GenericGrayImage`: the ``N``-dimensional input image
- `connectivity::Integer=1`: the neighborhood connectivity.
  The maximum number of orthogonal steps to reach a neighbor of the pixel.
  In 2D, it is 1 for a 4-neighborhood and 2 for a 8-neighborhood.
- `maxtree::MaxTree=nothing`: optional pre-built *max-tree*. Note that,
  when `maxtree` is specified, the `connectivity` option is ignored.

# Returns
An integer array of the same shape as the `image`. Pixels that are not local
maxima have 0 value. Pixels of the same local maximum share the same positive
value (the local maximum id).

# See also
[`MaxTree`](@ref), [`local_maxima!`](@ref), [`local_minima`](@ref),
[`Images.findlocalmaxima`](@ref)

# Examples
Create `f` (quadratic function with a maximum in the center and
4 additional constant maxima):
```jldoctest
julia> w = 10;
julia> f = [20 - 0.2*((x - w/2)^2 + (y-w/2)^2) for x in 0:w, y in 0:w];
julia> f[3:5, 3:5] .= 40; f[3:5, 8:10] .= 60; f[8:10, 3:5] .= 80; f[8:10, 8:10] .= 100;
```
Get all local maxima of ``f``:
```jldoctest
julia> f_maxima = local_maxima(f)
```
The resulting image contains the 4 labeled local maxima.
"""
local_maxima(image::GenericGrayImage; connectivity::Integer=1,
             maxtree::Union{MaxTree, Nothing} = nothing) =
    local_maxima!(similar(image, Int), image, connectivity=connectivity, maxtree=maxtree)

"""
    local_minima!(output, image; connectivity=1, maxtree=nothing) -> output

Detects the local minima of `image` and stores the result in `output`.
See [`local_minima`](@ref) for the detailed description of the method.
"""
function local_minima!(output::GenericGrayImage, image::GenericGrayImage;
                       connectivity::Integer=1,
                       maxtree::Union{MaxTree, Nothing} = nothing)
    check_output_image(output, image)
    _maxtree = check_maxtree(maxtree, image, connectivity=connectivity, rev=true)
    return local_extrema!(output, image, _maxtree)
end

"""
    local_minima(image::GenericGrayImage; connectivity=1, maxtree=nothing) -> Array

Determines and labels all *local minima* of the `image`.

# Details
The *local minimum* is defined as the connected set of pixels that have the
same value, which is less than the values of all pixels in direct
neighborhood of the set.

Technically, the implementation is based on the *max-tree* representation
of an image. It's beneficial if the max-tree is already computed,
otherwise [`Images.findlocalminima`](@ref) would be more efficient.

# Arguments
- `image::GenericGrayImage`: the ``N``-dimensional input image
- `connectivity::Integer=1`: the neighborhood connectivity. The maximum number
  of orthogonal steps to reach a neighbor of the pixel. In 2D, it is 1 for
  a 4-neighborhood and 2 for a 8-neighborhood.
- `maxtree::MaxTree=nothing`: optional pre-built *max-tree*. Note that, when
  `maxtree` is specified, the `connectivity` option is ignored.

# Returns
An integer array of the same shape as the `image`. Pixels that are not local
minima have 0 value. Pixels of the same local minimum share the same positive
value (the local minimum id).

# See also
[`MaxTree`](@ref), [`local_minima!`](@ref), [`local_maxima`](@ref),
[`Images.findlocalminima`](@ref)

# Examples
Create `f` (quadratic function with a minimum in the center and
4 additional constant minimum):
```jldoctest
julia> w = 10;
julia> f = [180 + 0.2*((x - w/2)^2 + (y-w/2)^2) for x in 0:w, y in 0:w];
julia> f[3:5, 3:5] .= 40; f[3:5, 8:10] .= 60; f[8:10, 3:5] .= 80; f[8:10, 8:10] .= 100;
```
Calculate all local minima of `f`:
```jldoctest
julia> f_minima = local_minima(f)
```
The resulting image contains the labeled local minima.
"""
local_minima(image::GenericGrayImage; connectivity::Integer=1,
             maxtree::Union{MaxTree, Nothing} = nothing) =
    local_minima!(similar(image, Int), image, connectivity=connectivity, maxtree=maxtree)
