"""
Max-tree morphological representation of an image.

# Type parameters
- `N`: the dimensions of the source image
- `A`: the type of the `axes(image)`

# Details
Let's consider image *thresholding* operation. The result is the image mask s.t.
``image[p] ≥ threshold`` for all ``p`` in the mask. This mask could also be
viewed as a set of connected components. When *image thresholding* is
sequentially applied for all possible thresholds, it generates a collection of
connected components that could be organized into a hierarchical structure
called *component tree*. A connected component at one threshold is a parent to
a component at the higher threshold if the latter is a subset of the first.

A *max-tree* is an efficient representation of the *component tree*.
A connected component at one level is represented by the single *reference pixel*
from this level, which is the parent to all other pixels of that level and
to the *reference pixel* of the level above.

The *max-tree* is the basis for many morphological operators,
namely connected operators. Unlike morphological openings and closings, these
operators do not require a fixed structuring element, but rather act with a
flexible structuring element that meets a certain criterion.

# See also
[`area_opening`](@ref), [`area_closing`](@ref),
[`diameter_opening`](@ref), [`diameter_closing`](@ref).

# References
 - Salembier, P., Oliveras, A., & Garrido, L. (1998). Antiextensive
   Connected Operators for Image and Sequence Processing.
   IEEE Transactions on Image Processing, 7(4), 555-570.
   :DOI:10.1109/83.663500
 - Berger, C., Geraud, T., Levillain, R., Widynski, N., Baillard, A.,
   Bertin, E. (2007). Effective Component Tree Computation with
   Application to Pattern Recognition in Astronomical Imaging.
   In International Conference on Image Processing (ICIP) (pp. 41-44).
   :DOI:10.1109/ICIP.2007.4379949
 - Najman, L., & Couprie, M. (2006). Building the component tree in
   quasi-linear time. IEEE Transactions on Image Processing, 15(11),
   3531-3539.
   :DOI:10.1109/TIP.2006.877518
 - Carlinet, E., & Geraud, T. (2014). A Comparative Review of
   Component Tree Computation Algorithms. IEEE Transactions on Image
   Processing, 23(9), 3885-3895.
   :DOI:10.1109/TIP.2014.2336551
"""
struct MaxTree{N,A}
    """
    Axes of the source image.
    """
    axes::A

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
    parents::Array{Int,N}

    """
    The order of the elements in the tree, from root to leaves.
    Each element corresponds to the linear index of the pixel in the image.
    A parent of a pixel always comes before the pixel itself, i.e.
    if index `i` comes before `j`, then `j`-th pixel cannot be
    the parent of `i`-th pixel.
    """
    traverse::Vector{Int}

    # create uninitialized MaxTree for image
    MaxTree{N}(image::AbstractArray{<:Any, N}, rev::Bool) where N =
        new{N, typeof(axes(image))}(axes(image), rev, similar(image, Int),
                                    Vector{Int}(undef, length(image)))
end

Base.ndims(::Type{<:MaxTree{N}}) where N = N
Base.ndims(maxtree::MaxTree) = ndims(typeof(maxtree))
Base.axes(maxtree::MaxTree) = maxtree.axes
Base.length(maxtree::MaxTree) = length(maxtree.parents)
Base.size(maxtree::MaxTree) = size(maxtree.parents)
Base.:(==)(a::MaxTree{N, A}, b::MaxTree{N, A}) where {N, A} =
    (a.rev == b.rev) && (a.axes == b.axes) &&
    (a.parents == b.parents) && (a.traverse == b.traverse)
Base.isequal(a::MaxTree{N, A}, b::MaxTree{N, A}) where {N, A} =
    isequal(a.rev, b.rev) && isequal(a.axes, b.axes) &&
    isequal(a.parents, b.parents) && isequal(a.traverse, b.traverse)

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

# Corrects `maxtree.parents`, so that the parent of every pixel
# is a canonical node. I.e. for each connected component at a level l,
# all pixels point to the same representative, which in turn points to the
# representative pixel at the next level.
#
# A canonical node ``p`` is either a root of the max-tree, or
# ``image(parent(p)) != image(p)``
function canonize!(maxtree::MaxTree{N}, image::AbstractArray{<:Any, N}) where N
    @assert size(maxtree) == size(image)
    parents = maxtree.parents
    @inbounds for p in maxtree.traverse
        q = parents[p]
        if image[q] == image[parents[q]]
            parents[p] = parents[q]
        end
    end
    return maxtree
end

# Checks whether the `pixel` + `offset` is still a valid pixel
# contained inside the image area defined by `axes`.
@generated isvalid_offset(offset::NTuple{N, Int}, pixel::CartesianIndex{N},
                          axes::NTuple{N}) where N =
    quote
        Base.Cartesian.@nall $N i -> in(pixel[i] + offset[i], axes[i])
    end

"""
    rebuild!(maxtree::MaxTree, image::AbstractArray,
             neighbors::AbstractVector{NTuple}) -> maxtree

Rebuilds the `maxtree` for the `image` using `neighbors` as the pixel
connectivity specification.

# Arguments
- `neighbors`: the vector of the offsets (``o_i``) to the cartesian index of
   the pixel (``p``), so that ``p + o_i`` is the cartesian index of the ``i``-th
   neighbor of ``p``.

# Details
The pixels in the connected components generated by the method should be
connected to each other by a path through neighboring pixels (as defined by
`neighbors`).

# See also
[`MaxTree`](@ref)
"""
function rebuild!(maxtree::MaxTree{N},
                  image::AbstractArray{<:Number, N},
                  #=mask::AbstractArray{Bool, N},=#
                  neighbors::AbstractVector{NTuple{N, Int}}) where N
    # pixels need to be sorted according to their gray level.
    sortperm!(maxtree.traverse, vec(image), rev=maxtree.rev)
    # some sufficiently inner point
    center = CartesianIndex{N}(ntuple(i -> 1 + size(image, i) ÷ 2, Val{N}()))
    # offsets for linear indices (note that they are sorted in ascending order)
    neighbor_offsets = getindex.(Ref(LinearIndices(image)),
        [CartesianIndex{N}(ntuple(i -> center[i] + n[i], Val{N}())) for n in neighbors]) .-
        getindex(LinearIndices(image), center)

    # initialization of the image parent
    parents = fill!(maxtree.parents, 0) |> vec
    roots = fill!(similar(parents), 0) # temporary array containing current roots
    rootpath = Vector{Int}() # temp array to hold the path to the current root within ancestors
    cindexes = CartesianIndices(maxtree.parents)

    # traverse the array in reversed order (from highest value to lowest value)
    @inbounds for p in Iterators.Reverse(maxtree.traverse)
        p_ci = cindexes[p]
        parents[p] = roots[p] = p # it's a new root

        for i in eachindex(neighbors)
            # check that p is in the mask and that it's neighbor is a valid image pixel
            (#=mask[p] && =#isvalid_offset(neighbors[i], p_ci, maxtree.axes)) || continue
            index = p + neighbor_offsets[i] # linear index of the neighbor
            (parents[index] == 0) && continue # ignore neighbor without parent (= its value is lower)
            root = rootpath!(rootpath, roots, index) # neighbor's root
            if (root != p) || (length(rootpath) > 1) # attach neighbor's root to the p
                parents[root] = p
                roots[rootpath] .= p # also attach the ancestral nodes of neighbor to p
            end
        end
    end
    return canonize!(maxtree, image)
end

# generates a vector of offsets to the CartesianIndex that defines
# the neighborhood of the pixel in an N-dimensional image.
# `connectivity` is the maximal number of orthogonal steps
# (+1/-1 to a single dimension) that are required to reach a neighbor.
function neighbor_cartesian_offsets(::Type{CartesianIndex{N}}, connectivity::Integer) where N
    cis = CartesianIndices(ntuple(_ -> 3, Val{N}()))
    offsets = [Tuple(ci) .- 2 for ci in cis]
    return filter!(offset -> 0 < sum(abs, offset) <= connectivity, vec(offsets))
end

"""
    MaxTree(image::AbstractArray; [connectivity=1], [rev=false]) -> MaxTree

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
MaxTree{2,Tuple{Base.OneTo{Int64},Base.OneTo{Int64}}}((Base.OneTo(3), Base.OneTo(3)), false, [4 2 4; 8 2 8; 2 2 2], [8, 2, 5, 6, 4, 9, 1, 3, 7])
```
"""
function MaxTree(image::AbstractArray{<:Number, N};
                 connectivity::Integer=1, rev::Bool=false) where N
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
    areas(maxtree::MaxTree) -> Vector{Int}

Computes the areas of all `maxtree` components.

# Returns
The vector of component areas. The `i`-th element is the area (in pixels) of
the component that is represented by the reference pixel with linear index `i`.

# See also
[`diameters`](@ref), [`area_opening`](@ref), [`area_closing`](@ref).
"""
function areas(maxtree::MaxTree)
    areas = fill(1, length(maxtree)) # start with 1-pixel areas
    parents = maxtree.parents
    @inbounds for p in Iterators.Reverse(maxtree.traverse)
        q = parents[p]
        (q != p) && (areas[q] += areas[p])
    end
    return areas
end

"""
    boundingboxes(maxtree::MaxTree) -> Matrix{Int}

Computes the minimal bounding boxes of all `maxtree` components.

# Returns
The matrix, where the `i`-th column encodes the bounding box of the subtree
with the root at the `i`-th pixel. The first ``N`` elements of the column
define the minimal cartesian index of the bounding box, the last ``N``
elements are its maximal cartesian index.

# See also
[`diameters`](@ref).
"""
function boundingboxes(maxtree::MaxTree{N}) where N
    # initialize bboxes
    bboxes = Matrix{Int}(undef, 2N, length(maxtree))
    @inbounds for (i, ci) in enumerate(CartesianIndices(maxtree.axes))
        offset = 2N * (i-1)
        bboxes[offset .+ (1:N)] .= Tuple(ci)
        bboxes[offset .+ ((N+1):2N)] .= Tuple(ci)
    end

    @inbounds for p in Iterators.Reverse(maxtree.traverse)
        q = maxtree.parents[p]
        for i in 1:N
            bboxes[i, q] = min(bboxes[i, q], bboxes[i, p])
        end
        for i in (N+1):2N
            bboxes[i, q] = max(bboxes[i, q], bboxes[i, p])
        end
    end
    return bboxes
end

"""
    diameters(maxtree::MaxTree) -> Vector{Int}

Computes the "diameters" of all `maxtree` components.

"Diameter" of the *max-tree* connected component is the length of
the widest side of the component's bounding box.

# Returns
The `i`-th element of the result is the "diameter" of the `i`-th component.

# See also
[`boundingboxes`](@ref), [`areas`](@ref),
[`diameter_opening`](@ref), [`diameter_closing`](@ref).
"""
diameters(maxtree::MaxTree{N}) where N =
    dropdims(mapslices(boundingboxes(maxtree), dims=1) do bbox
        maximum(ntuple(i -> bbox[i+N] - bbox[i] + 1, Val{N}()))
    end, dims=1)

"""
    direct_filter!(output::AbstractArray, image::AbstractArray,
                   maxtree::MaxTree, attrs::AbstractVector, min_attr) -> output

Applies a direct filtering of the `image` and stores the result in `output`.

Produces an image in which all components of its max-tree representation
have the specified attribute value no less than `min_attr`. It's done by
replacing the pixels of the component that has smaller attribute value with
the value of the reference pixel of its first valid ancestor.

# Arguments
- `maxtree::MaxTree{N}`: pre-built max-tree of the `image`
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
"""
function direct_filter!(output::AbstractArray, image::AbstractArray,
                        maxtree::MaxTree,
                        attrs::AbstractVector, min_attr,
                        all_below_min = zero(eltype(output)))
    # should have been already checked by higher-level functions
    @assert axes(output) == axes(image)
    @assert axes(image) == axes(maxtree)
    @assert length(attrs) == length(maxtree)

    p_root = root_index(maxtree)
    # if p_root is below min_attr, then all others are as well
    (attrs[p_root] < min_attr) && return fill!(output, all_below_min)
    output[p_root] = image[p_root]
    parents = maxtree.parents

    for p in maxtree.traverse
        q = parents[p]
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

function check_output_image(output::AbstractArray, image::AbstractArray)
    (size(output) == size(image)) ||
        throw(DimensionMismatch("The sizes of the output and the input image do not match"))
end

# checks if the provided maxtree is compatible with the given options
# or builds the new maxtree if none was given
function check_maxtree(maxtree::Union{MaxTree, Nothing},
                       image::AbstractArray;
                       connectivity::Integer = 0,
                       rev::Bool = false)
    (maxtree === nothing) && return MaxTree(image, connectivity=connectivity, rev=rev)
    (axes(maxtree) == axes(image)) ||
        throw(DimensionMismatch("The axes of the max-tree and the input image do not match"))
    (maxtree.rev == rev) ||
        throw(ArgumentError("The traversal order of the given max-tree is different from the requested one"))
    return maxtree
end

"""
    area_opening!(output, image;
                  [min_area=64], [connectivity=1], [maxtree=nothing]) -> output

Performs in-place *area opening* of the `image` and stores the result in `output`.
See [`area_opening`](@ref) for the detailed description of the method.
"""
function area_opening!(output::AbstractArray, image::AbstractArray;
                       min_area::Number=64, connectivity::Integer=1,
                       maxtree::Union{MaxTree, Nothing}=nothing)
    check_output_image(output, image)
    _maxtree = check_maxtree(maxtree, image, connectivity=connectivity, rev=false)
    return direct_filter!(output, image, _maxtree, areas(_maxtree), min_area)
end

"""
    area_opening(image; [min_area=64], [connectivity=1],
                 [maxtree=nothing]) -> Array

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
- `image::AbstractArray`: the ``N``-dimensional input image
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
- Vincent L., Proc. "Grayscale area openings and closings,
  their efficient implementation and applications",
  EURASIP Workshop on Mathematical Morphology and its
  Applications to Signal Processing, Barcelona, Spain, pp.22-27,
  May 1993.
- Soille, P., "Morphological Image Analysis: Principles and
  Applications" (Chapter 6), 2nd edition (2003), ISBN 3540429883.
  DOI:10.1007/978-3-662-05088-0
- Salembier, P., Oliveras, A., & Garrido, L. (1998). Antiextensive
  Connected Operators for Image and Sequence Processing.
  IEEE Transactions on Image Processing, 7(4), 555-570.
  DOI:10.1109/83.663500
- Najman, L., & Couprie, M. (2006). Building the component tree in
  quasi-linear time. IEEE Transactions on Image Processing, 15(11),
  3531-3539.
  DOI:10.1109/TIP.2006.877518
- Carlinet, E., & Geraud, T. (2014). A Comparative Review of
  Component Tree Computation Algorithms. IEEE Transactions on Image
  Processing, 23(9), 3885-3895.
  DOI:10.1109/TIP.2014.2336551

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
area_opening(image::AbstractArray; kwargs...) =
    area_opening!(similar(image), image; kwargs...)

"""
    diameter_opening!(output, image; [min_diameter=8],
                      [connectivity=1], [maxtree=nothing]) -> output

Performs in-place *diameter opening* of the `image` and stores the result in `output`.
See [`diameter_opening`](@ref) for the detailed description of the method.
"""
function diameter_opening!(output::AbstractArray, image::AbstractArray;
                           maxtree::Union{MaxTree, Nothing} = nothing,
                           min_diameter=8, connectivity=1)
    check_output_image(output, image)
    _maxtree = check_maxtree(maxtree, image, connectivity=connectivity, rev=false)
    return direct_filter!(output, image, _maxtree, diameters(_maxtree), min_diameter)
end

"""
    diameter_opening(image; [min_diameter=8], [connectivity=1],
                     [maxtree=nothing]) -> Array

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
- `image::AbstractArray`: the ``N``-dimensional input image
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
- Walter, T., & Klein, J.-C. (2002). Automatic Detection of
  Microaneurysms in Color Fundus Images of the Human Retina by Means
  of the Bounding Box Closing. In A. Colosimo, P. Sirabella,
  A. Giuliani (Eds.), Medical Data Analysis. Lecture Notes in Computer
  Science, vol 2526, pp. 210-220. Springer Berlin Heidelberg.
  :DOI:`10.1007/3-540-36104-9_23`
- Carlinet, E., & Geraud, T. (2014). A Comparative Review of
  Component Tree Computation Algorithms. IEEE Transactions on Image
  Processing, 23(9), 3885-3895.
  :DOI:`10.1109/TIP.2014.2336551`

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
diameter_opening(image::AbstractArray; kwargs...) =
    diameter_opening!(similar(image), image; kwargs...)

"""
    area_closing!(output, image; [min_area=64], [connectivity=1],
                  [maxtree=nothing]) -> output

Performs in-place *area closing* of the `image` and stores the result in `output`.
See [`area_closing`](@ref) for the detailed description of the method.
"""
function area_closing!(output::AbstractArray, image::AbstractArray;
                       min_area::Number=64, connectivity::Integer=1,
                       maxtree::Union{MaxTree, Nothing}=nothing)
    check_output_image(output, image)
    _maxtree = check_maxtree(maxtree, image, connectivity=connectivity, rev=true)
    return direct_filter!(output, image, _maxtree, areas(_maxtree), min_area)
end

"""
    area_closing(image; [min_area=64], [connectivity=1],
                 [maxtree=nothing]) -> Array

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
- `image::AbstractArray`: the ``N``-dimensional input image
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
- Vincent L., Proc. "Grayscale area openings and closings,
  their efficient implementation and applications",
  EURASIP Workshop on Mathematical Morphology and its
  Applications to Signal Processing, Barcelona, Spain, pp.22-27,
  May 1993.
- Soille, P., "Morphological Image Analysis: Principles and
  Applications" (Chapter 6), 2nd edition (2003), ISBN 3540429883.
  DOI:10.1007/978-3-662-05088-0
- Salembier, P., Oliveras, A., & Garrido, L. (1998). Antiextensive
  Connected Operators for Image and Sequence Processing.
  IEEE Transactions on Image Processing, 7(4), 555-570.
  DOI:10.1109/83.663500
- Najman, L., & Couprie, M. (2006). Building the component tree in
  quasi-linear time. IEEE Transactions on Image Processing, 15(11),
  3531-3539.
  DOI:10.1109/TIP.2006.877518
- Carlinet, E., & Geraud, T. (2014). A Comparative Review of
  Component Tree Computation Algorithms. IEEE Transactions on Image
  Processing, 23(9), 3885-3895.
  DOI:10.1109/TIP.2014.2336551

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
area_closing(image::AbstractArray; kwargs...) =
    area_closing!(similar(image), image; kwargs...)

"""
    diameter_closing!(output, image; [min_diameter=8], [connectivity=1],
                      [maxtree=nothing]) -> output

Performs in-place *diameter closing* of the `image` and stores the result in `output`.
See [`diameter_closing`](@ref) for the detailed description of the method.
"""
function diameter_closing!(output::AbstractArray, image::AbstractArray;
                           min_diameter::Number=8, connectivity::Integer=1,
                           maxtree::Union{MaxTree, Nothing} = nothing)
    check_output_image(output, image)
    _maxtree = check_maxtree(maxtree, image, connectivity=connectivity, rev=true)
    return direct_filter!(output, image, _maxtree, diameters(_maxtree), min_diameter)
end

"""
    diameter_closing(image; [min_diameter=8], [connectivity=1],
                     [maxtree=nothing]) -> Array

Performs a *diameter closing* of the `image`.

*Diameter closing* replaces all dark structures of an image that have
the diameter (the widest dimension of their bounding box) smaller than
`min_diameter` with the brighter value taken from their first ancestral component
(in *max-tree* representation of `image`) that has the diameter no smaller than
`min_diameter`.

# Arguments
- `image::AbstractArray`: the ``N``-dimensional input image
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
- Walter, T., & Klein, J.-C. (2002). Automatic Detection of
  Microaneurysms in Color Fundus Images of the Human Retina by Means
  of the Bounding Box Closing. In A. Colosimo, P. Sirabella,
  A. Giuliani (Eds.), Medical Data Analysis. Lecture Notes in Computer
  Science, vol 2526, pp. 210-220. Springer Berlin Heidelberg.
  :DOI:`10.1007/3-540-36104-9_23`
- Carlinet, E., & Geraud, T. (2014). A Comparative Review of
  Component Tree Computation Algorithms. IEEE Transactions on Image
  Processing, 23(9), 3885-3895.
  :DOI:`10.1109/TIP.2014.2336551`

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
diameter_closing(image::AbstractArray; kwargs...) =
    diameter_closing!(similar(image), image; kwargs...)

# Calculates the local maxima or minima of the image using the max-tree
# representation (maxima if maxtree.rev=false, minima otherwise).
# It's not the most effecient method for extrema calculation, but could be
# useful if the max-tree is already calculated.
# Each minima/maxima region is labeled with the unique id (1 is the global
# maximum/minimum).
function local_extrema!(output::AbstractArray{<:Any, N},
                        image::AbstractArray{<:Any, N},
                        maxtree::MaxTree{N}) where N
    (axes(image) == axes(maxtree)) || throw(DimensionMismatch())
    (size(output) == size(image)) || throw(DimensionMismatch())

    next_label = one(eltype(output))
    fill!(output, next_label) # initialize output (just has to be non-zero)
    @inbounds for p in Iterators.reverse(maxtree.traverse)
        q = maxtree.parents[p]
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
        q = maxtree.parents[p]
        # if p is not canonical (parent has the same value)
        if image[p] == image[q]
            # in this case we propagate the value
            output[p] = output[q]
        end
    end

    return output
end

"""
    local_maxima!(output, image; [connectivity=1], [maxtree=nothing]) -> output

Detects the local maxima of `image` and stores the result in `output`.
See [`local_maxima`](@ref) for the detailed description of the method.
"""
function local_maxima!(output::AbstractArray, image::AbstractArray;
                       connectivity::Integer=1,
                       maxtree::Union{MaxTree, Nothing} = nothing)
    check_output_image(output, image)
    _maxtree = check_maxtree(maxtree, image, connectivity=connectivity, rev=false)
    return local_extrema!(output, image, _maxtree)
end

"""
    local_maxima(image::AbstractArray; [connectivity=1], [maxtree=nothing]) -> Array

Determines and labels all *local maxima* of the `image`.

# Details
The *local maximum* is defined as the connected set of pixels that have the
same value, which is grater then the values of all pixels in direct
neighborhood of the set.

Technically, the implementation is based on the *max-tree* representation
of an image. It's efficient if the max-tree representation
has already been computed. Otherwise, it is preferable to use
the function `local_maxima`.

# Arguments
- `image::AbstractArray`: the ``N``-dimensional input image
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
[`MaxTree`](@ref), [`local_maxima!`](@ref), [`local_minima`](@ref)

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
local_maxima(image::AbstractArray; connectivity::Integer=1,
             maxtree::Union{MaxTree, Nothing} = nothing) =
    local_maxima!(similar(image, Int), image, connectivity=connectivity, maxtree=maxtree)

"""
    local_minima!(output, image; [connectivity=1], [maxtree=nothing]) -> output

Detects the local minima of `image` and stores the result in `output`.
See [`local_minima`](@ref) for the detailed description of the method.
"""
function local_minima!(output::AbstractArray, image::AbstractArray;
                       connectivity::Integer=1,
                       maxtree::Union{MaxTree, Nothing} = nothing)
    check_output_image(output, image)
    _maxtree = check_maxtree(maxtree, image, connectivity=connectivity, rev=true)
    return local_extrema!(output, image, _maxtree)
end

"""
    local_minima(image::AbstractArray; [connectivity=1], [maxtree=nothing]) -> Array

Determines and labels all *local minima* of the `image`.

# Details
The *local minimum* is defined as the connected set of pixels that have the
same value, which is less then the values of all pixels in direct
neighborhood of the set.

Technically, the implementation is based on the *max-tree* representation
of an image. It's efficient if the max-tree representation
has already been computed. Otherwise, it is preferable to use
the function `local_minima`.

# Arguments
- `image::AbstractArray`: the ``N``-dimensional input image
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
[`MaxTree`](@ref), [`local_minima!`](@ref), [`local_maxima`](@ref)

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
local_minima(image::AbstractArray; connectivity::Integer=1,
             maxtree::Union{MaxTree, Nothing} = nothing) =
    local_minima!(similar(image, Int), image, connectivity=connectivity, maxtree=maxtree)
