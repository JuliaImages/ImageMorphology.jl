"""max_tree.py - max_tree representation of images.

This module provides operators based on the max-tree representation of images.
A grayscale image can be seen as a pile of nested sets, each of which is the
result of a threshold operation. These sets can be efficiently represented by
max-trees, where the inclusion relation between connected components at
different levels are represented by parent-child relationships.

These representations allow efficient implementations of many algorithms, such
as attribute operators. Unlike morphological openings and closings, these
operators do not require a fixed structuring element, but rather act with a
flexible structuring element that meets a certain criterion.

This implementation provides functions for:
1. max-tree generation
2. area openings / closings
3. diameter openings / closings
4. local maxima

References:
    .. [1] Salembier, P., Oliveras, A., & Garrido, L. (1998). Antiextensive
           Connected Operators for Image and Sequence Processing.
           IEEE Transactions on Image Processing, 7(4), 555-570.
           :DOI:10.1109/83.663500
    .. [2] Berger, C., Geraud, T., Levillain, R., Widynski, N., Baillard, A.,
           Bertin, E. (2007). Effective Component Tree Computation with
           Application to Pattern Recognition in Astronomical Imaging.
           In International Conference on Image Processing (ICIP) (pp. 41-44).
           :DOI:10.1109/ICIP.2007.4379949
    .. [3] Najman, L., & Couprie, M. (2006). Building the component tree in
           quasi-linear time. IEEE Transactions on Image Processing, 15(11),
           3531-3539.
           :DOI:10.1109/TIP.2006.877518
    .. [4] Carlinet, E., & Geraud, T. (2014). A Comparative Review of
           Component Tree Computation Algorithms. IEEE Transactions on Image
           Processing, 23(9), 3885-3895.
           :DOI:10.1109/TIP.2014.2336551
"""

"""_max_tree.pyx - building a max-tree from an image.

This is an implementation of the max-tree, which is a morphological
representation of the image. Many morphological operators can be built
from this representation, namely attribute openings and closings.

This file also contains implementations of max-tree based filters and
functions to characterize the tree components.
"""

"""
Max-tree representation of an image.
"""
struct MaxTree{N,A}
    """
    Axes of the source image.
    """
    axes::A

    """
    Array of the same shape as the source image.
    Each value is the linear index of the parent pixel in the max-tree.
    """
    parents::Array{Int,N}

    """
    Vector of the same length as the source image.
    Each element corresponds to linear index of the pixel in the image.
    The vector encodes the order of elements in the tree:
    a parent of a pixel always comes before the element itself.
    More formally: if `i` comes before `j`, then `j` cannot be the parent of `i`.
    """
    traverse::Vector{Int}

    # create uninitialized MaxTree for image
    MaxTree{N}(image::AbstractArray{<:Any, N}) where N =
        new{N, typeof(axes(image))}(axes(image), similar(image, Int),
                                    Vector{Int}(undef, length(image)))
end

Base.axes(maxtree::MaxTree) = maxtree.axes
Base.length(maxtree::MaxTree) = length(maxtree.parents)
Base.size(maxtree::MaxTree) = size(maxtree.parents)
root_index(maxtree::MaxTree) = maxtree.traverse[1]

"""
Get the path to root of the current tree starting from `index`.

Here, we do without path compression and accept the higher complexity, but
the function is inline and avoids some overhead induced by its recursive
version.

Parameters
----------
parent : The array containing parent relationships.
index :  The index of which we want to find the root.

Returns
-------
root : int
    The root found from ``index``.
"""
@inline function rootpath!(path::Vector{Int}, parents::AbstractArray{<:Integer}, index::Integer)
    empty!(path)
    p = index
    @label next
    push!(path, p)
    @inbounds q = parents[p]
    if p != q
        p = q
        @goto next
    end
    return p
end

"""
Corrects the given `maxtree`, so that every node's parent is a canonical node.

The parent of a non-canonical pixel is a canonical pixel.
The parent of a canonical pixel is also a canonical pixel with a different
value. There is exactly one canonical pixel for each component in the
component tree.

Parameters
----------
maxtree : The max-tree of the same shape as the image
image : The source image of `maxtree`
"""
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

"""
Checks whether a neighbor of a given pixel is inside the image.

Parameters
----------
offset: the offsets of each cartesian index component of the pixel
index : index of the image pixel
axes : image axes

Returns
-------
is_neighbor : uint8
    0 if the neighbor falls outside the image, 1 otherwise.
"""
@generated function isvalid_offset(offset::NTuple{N, Int},
                                   index::CartesianIndex{N},
                                   axes::NTuple{N}) where N
    quote
        Base.Cartesian.@nall $N i -> in(index[i] + offset[i], axes[i])
    end
end

isvalid_offset(offset::NTuple{N, Int}, index::Integer, axes::NTuple{N}) where N =
    isvalid_offset(offset, LinearIndices(axes)[index], axes)

"""
Compute the area of all max-tree components.

This attribute is used for area opening and closing
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
Compute the bounding box extension of all max-tree components.

This attribute is used for diameter opening and closing.
"""
function boundingboxes(maxtree::MaxTree{N}) where N
    # initialize bboxes
    bboxes = Matrix{Int}(undef, 2N, length(maxtree))
    for (bbox, ci) in zip(eachcol(bboxes), CartesianIndices(maxtree.axes))
        @inbounds bbox[1:N] .= Tuple(ci)
        @inbounds bbox[(N+1):2N] .= Tuple(ci)
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

diameters(maxtree::MaxTree{N}) where N =
    [maximum(ntuple(i -> bbox[i+N] - bbox[i] + 1, Val{N}()))
     for bbox in eachcol(boundingboxes(maxtree))]

"""
Find the local maxima in image from the max-tree representation.

Parameters
----------

image : array of arbitrary type
    The flattened image pixels.
output : array of the same shape and type as image.
    The output image must contain only ones.
parent : array of int
    Image of the same shape as the input image. The value
    at each pixel is the parent index of this pixel in the max-tree
    reprentation.
sorted_indices : array of int
    List of length = number of pixels. Each element
    corresponds to one pixel index in the image. It encodes the order
    of elements in the tree: a parent of a pixel always comes before
    the element itself. More formally: i < j implies that j cannot be
    the parent of i.
"""
# _max_tree_local_maxima cacluates the local maxima from the max-tree
# representation this is interesting if the max-tree representation has
# already been calculated for other reasons. Otherwise, it is not the most
# efficient method. If the parameter label is True, the minima are labeled.
function local_maxima!(output::AbstractArray{<:Any, N},
                       image::AbstractArray{<:Any, N},
                       maxtree::MaxTree{N}) where N
    (size(image) == size(maxtree)) || throw(DimensionMismatch())
    (size(output) == size(image)) || throw(DimensionMismatch())

    next_label = 1
    for p in Iterators.reverse(maxtree.traverse)
        q = parent[p]
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

    for p in Iterators.reverse(maxtree.traverse)
        q = parent[p]
        # if p is not canonical (parent has the same value)
        if image[p] == image[q]
            # in this case we propagate the value
            output[p] = output[q]
        end
    end

    return output
end

"""
Apply a direct filtering.

This produces an image in which for all possible thresholds, each connected
component has the specified attribute value greater than that threshold.
This is the basic function called by :func:`area_opening`,
:func:`diameter_opening`, and similar.

For :func:`area_opening`, for instance, the attribute is the area.  In this
case, an image is produced for which all connected components for all
thresholds have at least an area (pixel count) of the threshold given by
the user.

Parameters
----------

image : array
    The flattened image pixels.
output : array, same size and type as `image`
    The array into which to write the output values. **This array will be
    modified in-place.**
maxtree : max-tree of `image`
attrs : array of float
    Contains the attributes computed for the max-tree.
min_attr : float
    The threshold to be applied to the attribute.
"""
function direct_filter!(output::AbstractArray{T1, N},
                        image::AbstractArray{T2, N},
                        maxtree::MaxTree{N},
                        attrs::AbstractVector{<:Number},
                        min_attr::Number) where {T1, T2, N}
    @assert axes(output) == axes(image)
    @assert axes(image) == axes(maxtree)
    @assert length(attrs) == length(maxtree)
    p_root = root_index(maxtree)
    output[p_root] = attrs[p_root] < min_attr ? 0 : image[p_root]
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

# _max_tree is the main function. It allows to construct a max
# tree representation of the image.
"""
Rebuilds `maxtree` to be the max-tree for `image`.

Parameters
----------
image : array
    The flattened image pixels.
mask : array of int
    An array of the same shape as `image` where each pixel contains a
    nonzero value if it is to be considered for the filtering.  NOTE: it is
    *essential* that the border pixels (those with neighbors falling
    outside the volume) are all set to zero, or segfaults could occur.
structure : array of int
    A list of coordinate offsets to compute the raveled coordinates of each
    neighbor from the raveled coordinates of the current pixel.
"""
function rebuild!(maxtree::MaxTree{N},
                  image::AbstractArray{<:Number, N},
                  #=mask::AbstractArray{Bool, N},=#
                  neighbors::AbstractVector{NTuple{N, Int}};
                  rev::Bool = false) where N
    # pixels need to be sorted according to their gray level.
    sortperm!(maxtree.traverse, vec(image), rev=rev)
    # some sufficiently inner point
    center = CartesianIndex{N}(ntuple(i -> size(image, i) ÷ 2, Val{N}()))
    # offsets for linear indices (note that they are sorted in ascending order)
    neighbor_offsets = getindex.(Ref(LinearIndices(image)),
        [CartesianIndex{N}(ntuple(i -> center[i] + n[i], Val{N}())) for n in neighbors]) .-
        getindex(LinearIndices(image), center)

    # initialization of the image parent
    parents = fill!(maxtree.parents, 0) |> vec
    roots = fill!(similar(parents), 0) # temporary array containing current roots
    rootpath = Vector{Int}() # temp array to hold the path to the current root
    cindexes = CartesianIndices(maxtree.parents)

    # traverse the array in reversed order (from highest value to lowest value)
    @inbounds for p in Iterators.Reverse(maxtree.traverse)
        p_ci = cindexes[p]
        parents[p] = roots[p] = p # it's a new root

        for i in eachindex(neighbors)
            # check that p is in the mask and that it's neighbor is a valid image pixel
            (#=mask[p] && =#isvalid_offset(neighbors[i], p_ci, maxtree.axes)) || continue
            index = p + neighbor_offsets[i] # linear index of the neighbor
            # ignore unset parent (it's value is higher)
            (parents[index] > 0) || continue
            root = rootpath!(rootpath, roots, index) # neighbor's root
            if (root != p) || (length(rootpath) > 1) # attach neighbor's root to the p
                parents[root] = p
                roots[rootpath] .= p
            end
        end
    end
    # In a canonized max-tree, each parent is a canonical pixel,
    # i.e. for each connected component at a level l, all pixels point
    # to the same representative which in turn points to the representative
    # pixel at the next level.
    return canonize!(maxtree, image)
end

"""
Convert any valid connectivity to a structuring element and offset.

Parameters
----------
image_dim : int
    The number of dimensions of the input image.
connectivity : int, array, or None
    The neighborhood connectivity. An integer is interpreted as in
    ``scipy.ndimage.generate_binary_structure``, as the maximum number
    of orthogonal steps to reach a neighbor. An array is directly
    interpreted as a structuring element and its shape is validated against
    the input image shape. ``None`` is interpreted as a connectivity of 1.
offset : tuple of int, or None
    The coordinates of the center of the structuring element.
Returns
-------
c_connectivity : array of bool
    The structuring element corresponding to the input `connectivity`.
offset : array of int
    The offset corresponding to the center of the structuring element.
Raises
------
ValueError:
    If the image dimension and the connectivity or offset dimensions don't
    match.
"""
function neighbor_cartesian_offsets(::Type{CartesianIndex{N}}, connectivity::Integer) where N
    cis = CartesianIndices(ntuple(_ -> 3, Val{N}()))
    offsets = [Tuple(ci) .- 2 for ci in cis]
    return filter!(offset -> 0 < sum(abs, offset) <= connectivity, vec(offsets))
end

"""
Build the max tree from an image.

Component trees represent the hierarchical structure of the connected
components resulting from sequential thresholding operations applied to an
image. A connected component at one level is parent of a component at a
higher level if the latter is included in the first. A max-tree is an
efficient representation of a component tree. A connected component at
one level is represented by one reference pixel at this level, which is
parent to all other pixels at that level and to the reference pixel at the
level above. The max-tree is the basis for many morphological operators,
namely connected operators.

Parameters
----------
image: ndarray
    The input image for which the max-tree is to be calculated.
    This image can be of any type.
connectivity: unsigned int, optional
    The neighborhood connectivity. The integer represents the maximum
    number of orthogonal steps to reach a neighbor. In 2D, it is 1 for
    a 4-neighborhood and 2 for a 8-neighborhood. Default value is 1.

Returns
-------
parent: ndarray, int64
    Array of same shape as image. The value of each pixel is the index of
    its parent in the ravelled array.

tree_traverser: 1D array, int64
    The ordered pixel indices (referring to the ravelled array). The pixels
    are ordered such that every pixel is preceded by its parent (except for
    the root which has no parent).

References
----------
.. [1] Salembier, P., Oliveras, A., & Garrido, L. (1998). Antiextensive
       Connected Operators for Image and Sequence Processing.
       IEEE Transactions on Image Processing, 7(4), 555-570.
       :DOI:`10.1109/83.663500`
.. [2] Berger, C., Geraud, T., Levillain, R., Widynski, N., Baillard, A.,
       Bertin, E. (2007). Effective Component Tree Computation with
       Application to Pattern Recognition in Astronomical Imaging.
       In International Conference on Image Processing (ICIP) (pp. 41-44).
       :DOI:`10.1109/ICIP.2007.4379949`
.. [3] Najman, L., & Couprie, M. (2006). Building the component tree in
       quasi-linear time. IEEE Transactions on Image Processing, 15(11),
       3531-3539.
       :DOI:`10.1109/TIP.2006.877518`
.. [4] Carlinet, E., & Geraud, T. (2014). A Comparative Review of
       Component Tree Computation Algorithms. IEEE Transactions on Image
       Processing, 23(9), 3885-3895.
       :DOI:`10.1109/TIP.2014.2336551`

Examples
--------
We create a small sample image (Figure 1 from [4]) and build the max-tree.

>>> image = np.array([[15, 13, 16], [12, 12, 10], [16, 12, 14]])
>>> P, S = max_tree(image, connectivity=2)
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
    return rebuild!(MaxTree{N}(image), #=mask,=# image, neighbors, rev=rev)
end

"""
Perform an area opening of the image.

Area opening removes all bright structures of an image with
a surface smaller than area_threshold.
The output image is thus the largest image smaller than the input
for which all local maxima have at least a surface of
area_threshold pixels.

Area openings are similar to morphological openings, but
they do not use a fixed structuring element, but rather a deformable
one, with surface = area_threshold. Consequently, the area_opening
with area_threshold=1 is the identity.

In the binary case, area openings are equivalent to
remove_small_objects; this operator is thus extended to gray-level images.

Technically, this operator is based on the max-tree representation of
the image.

Parameters
----------
- image: ndarray
    The input image for which the area_opening is to be calculated.
    This image can be of any type.
- min_area: unsigned int
    The size parameter (number of pixels). The default value is arbitrarily
    chosen to be 64.
- connectivity: unsigned int, optional
    The neighborhood connectivity. The integer represents the maximum
    number of orthogonal steps to reach a neighbor. In 2D, it is 1 for
    a 4-neighborhood and 2 for a 8-neighborhood. Default value is 1.
- parent: ndarray, int64, optional
    Parent image representing the max tree of the image. The
    value of each pixel is the index of its parent in the ravelled array.
tree_traverser: 1D array, int64, optional
    The ordered pixel indices (referring to the ravelled array). The pixels
    are ordered such that every pixel is preceded by its parent (except for
    the root which has no parent).

# Returns
output: ndarray
    Output image of the same shape and type as the input image.

# See also
[`area_closing`](@ref), [`diameter_opening`](@ref), [`diameter_closing`](@ref),
[`MaxTree`](@ref)

# References
> Vincent L., Proc. "Grayscale area openings and closings,
  their efficient implementation and applications",
  EURASIP Workshop on Mathematical Morphology and its
  Applications to Signal Processing, Barcelona, Spain, pp.22-27,
  May 1993.
> Soille, P., "Morphological Image Analysis: Principles and
  Applications" (Chapter 6), 2nd edition (2003), ISBN 3540429883.
  DOI:10.1007/978-3-662-05088-0
> Salembier, P., Oliveras, A., & Garrido, L. (1998). Antiextensive
  Connected Operators for Image and Sequence Processing.
  IEEE Transactions on Image Processing, 7(4), 555-570.
  DOI:10.1109/83.663500
> Najman, L., & Couprie, M. (2006). Building the component tree in
  quasi-linear time. IEEE Transactions on Image Processing, 15(11),
  3531-3539.
  DOI:10.1109/TIP.2006.877518
> Carlinet, E., & Geraud, T. (2014). A Comparative Review of
  Component Tree Computation Algorithms. IEEE Transactions on Image
  Processing, 23(9), 3885-3895.
  DOI:10.1109/TIP.2014.2336551

# Examples

We create an image (quadratic function with a maximum in the center and
4 additional local maxima.
```juliacmd
w = 12
x, y = np.mgrid[0:w,0:w]
f = 20 - 0.2*((x - w/2)**2 + (y-w/2)**2)
f[2:3,1:5] = 40; f[2:4,9:11] = 60; f[9:11,2:4] = 80
f[9:10,9:11] = 100; f[10,10] = 100
f = f.astype(np.int)
```
We can calculate the area opening:
```juliacmd
open = area_opening(f, 8, connectivity=1)
```
The peaks with a surface smaller than 8 are removed.
"""
function area_opening!(output::AbstractArray{<:Any, N}, image::AbstractArray{<:Any, N};
                       min_area::Number=64, connectivity::Integer=1,
                       maxtree::Union{MaxTree{N}, Nothing}=nothing) where N
    _maxtree = isnothing(maxtree) ? MaxTree(image, connectivity=connectivity) : maxtree
    return direct_filter!(output, image, _maxtree, areas(_maxtree), min_area)
end

area_opening(image::AbstractArray; kwargs...) =
    area_opening!(similar(image), image; kwargs...)

"""
Perform a diameter opening of the image.

Diameter opening removes all bright structures of an image with
maximal extension smaller than diameter_threshold. The maximal
extension is defined as the maximal extension of the bounding box.
The operator is also called Bounding Box Opening. In practice,
the result is similar to a morphological opening, but long and thin
structures are not removed.

Technically, this operator is based on the max-tree representation of
the image.

Parameters
----------
image: ndarray
    The input image for which the area_opening is to be calculated.
    This image can be of any type.
diameter_threshold: unsigned int
    The maximal extension parameter (number of pixels). The default value
    is 8.
connectivity: unsigned int, optional
    The neighborhood connectivity. The integer represents the maximum
    number of orthogonal steps to reach a neighbor. In 2D, it is 1 for
    a 4-neighborhood and 2 for a 8-neighborhood. Default value is 1.
parent: ndarray, int64, optional
    Parent image representing the max tree of the image. The
    value of each pixel is the index of its parent in the ravelled array.
tree_traverser: 1D array, int64, optional
    The ordered pixel indices (referring to the ravelled array). The pixels
    are ordered such that every pixel is preceded by its parent (except for
    the root which has no parent).

Returns
-------
output: ndarray
    Output image of the same shape and type as the input image.

See also
--------
skimage.morphology.area_opening
skimage.morphology.area_closing
skimage.morphology.diameter_closing
skimage.morphology.max_tree

References
----------
.. [1] Walter, T., & Klein, J.-C. (2002). Automatic Detection of
       Microaneurysms in Color Fundus Images of the Human Retina by Means
       of the Bounding Box Closing. In A. Colosimo, P. Sirabella,
       A. Giuliani (Eds.), Medical Data Analysis. Lecture Notes in Computer
       Science, vol 2526, pp. 210-220. Springer Berlin Heidelberg.
       :DOI:`10.1007/3-540-36104-9_23`
.. [2] Carlinet, E., & Geraud, T. (2014). A Comparative Review of
       Component Tree Computation Algorithms. IEEE Transactions on Image
       Processing, 23(9), 3885-3895.
       :DOI:`10.1109/TIP.2014.2336551`

Examples
--------
We create an image (quadratic function with a maximum in the center and
4 additional local maxima.

>>> w = 12
>>> x, y = np.mgrid[0:w,0:w]
>>> f = 20 - 0.2*((x - w/2)**2 + (y-w/2)**2)
>>> f[2:3,1:5] = 40; f[2:4,9:11] = 60; f[9:11,2:4] = 80
>>> f[9:10,9:11] = 100; f[10,10] = 100
>>> f = f.astype(np.int)

We can calculate the diameter opening:

>>> open = diameter_opening(f, 3, connectivity=1)

The peaks with a maximal extension of 2 or less are removed.
The remaining peaks have all a maximal extension of at least 3.
"""
function diameter_opening!(output::AbstractArray{<:Any, N}, image::AbstractArray{<:Any, N};
                           maxtree::Union{MaxTree{N}, Nothing} = nothing,
                           min_diameter=8, connectivity=1) where N
    _maxtree = isnothing(maxtree) ? MaxTree(image, connectivity=connectivity) : maxtree
    return direct_filter!(output, image, _maxtree, diameters(_maxtree), min_diameter)
end

diameter_opening(image::AbstractArray; kwargs...) =
    diameter_opening!(similar(image), image; kwargs...)

"""
Perform an area closing of the image.

Area closing removes all dark structures of an image with
a surface smaller than area_threshold.
The output image is larger than or equal to the input image
for every pixel and all local minima have at least a surface of
area_threshold pixels.

Area closings are similar to morphological closings, but
they do not use a fixed structuring element, but rather a deformable
one, with surface = area_threshold.

In the binary case, area closings are equivalent to
remove_small_holes; this operator is thus extended to gray-level images.

Technically, this operator is based on the max-tree representation of
the image.

Parameters
----------
image: ndarray
    The input image for which the area_closing is to be calculated.
    This image can be of any type.
area_threshold: unsigned int
    The size parameter (number of pixels). The default value is arbitrarily
    chosen to be 64.
connectivity: unsigned int, optional
    The neighborhood connectivity. The integer represents the maximum
    number of orthogonal steps to reach a neighbor. In 2D, it is 1 for
    a 4-neighborhood and 2 for a 8-neighborhood. Default value is 1.
parent: ndarray, int64, optional
    Parent image representing the max tree of the inverted image. The
    value of each pixel is the index of its parent in the ravelled array.
    See Note for further details.
tree_traverser: 1D array, int64, optional
    The ordered pixel indices (referring to the ravelled array). The pixels
    are ordered such that every pixel is preceded by its parent (except for
    the root which has no parent).

Returns
-------
output: ndarray
    Output image of the same shape and type as input image.

See also
--------
skimage.morphology.area_opening
skimage.morphology.diameter_opening
skimage.morphology.diameter_closing
skimage.morphology.max_tree
skimage.morphology.remove_small_objects
skimage.morphology.remove_small_holes

References
----------
.. [1] Vincent L., Proc. "Grayscale area openings and closings,
       their efficient implementation and applications",
       EURASIP Workshop on Mathematical Morphology and its
       Applications to Signal Processing, Barcelona, Spain, pp.22-27,
       May 1993.
.. [2] Soille, P., "Morphological Image Analysis: Principles and
       Applications" (Chapter 6), 2nd edition (2003), ISBN 3540429883.
       :DOI:`10.1007/978-3-662-05088-0`
.. [3] Salembier, P., Oliveras, A., & Garrido, L. (1998). Antiextensive
       Connected Operators for Image and Sequence Processing.
       IEEE Transactions on Image Processing, 7(4), 555-570.
       :DOI:`10.1109/83.663500`
.. [4] Najman, L., & Couprie, M. (2006). Building the component tree in
       quasi-linear time. IEEE Transactions on Image Processing, 15(11),
       3531-3539.
       :DOI:`10.1109/TIP.2006.877518`
.. [5] Carlinet, E., & Geraud, T. (2014). A Comparative Review of
       Component Tree Computation Algorithms. IEEE Transactions on Image
       Processing, 23(9), 3885-3895.
       :DOI:`10.1109/TIP.2014.2336551`


Examples
--------
We create an image (quadratic function with a minimum in the center and
4 additional local minima.

>>> w = 12
>>> x, y = np.mgrid[0:w,0:w]
>>> f = 180 + 0.2*((x - w/2)**2 + (y-w/2)**2)
>>> f[2:3,1:5] = 160; f[2:4,9:11] = 140; f[9:11,2:4] = 120
>>> f[9:10,9:11] = 100; f[10,10] = 100
>>> f = f.astype(np.int)

We can calculate the area closing:

>>> closed = area_closing(f, 8, connectivity=1)

All small minima are removed, and the remaining minima have at least
a size of 8.


Notes
-----
If a max-tree representation (parent and tree_traverser) are given to the
function, they must be calculated from the inverted image for this
function, i.e.:
>>> P, S = max_tree(invert(f))
>>> closed = diameter_closing(f, 3, parent=P, tree_traverser=S)
"""
function area_closing!(output::AbstractArray{<:Any, N},
                       image::AbstractArray{<:Any, N};
                       min_area::Number=64, connectivity::Integer=1,
                       maxtree::Union{MaxTree{N}, Nothing}=nothing) where N
    _maxtree = isnothing(maxtree) ? MaxTree(image, connectivity=connectivity, rev=true) : maxtree
    return direct_filter!(output, image, _maxtree, areas(_maxtree), min_area)
end

area_closing(image::AbstractArray; kwargs...) =
    area_closing!(similar(image), image; kwargs...)

"""
Perform a diameter closing of the image.

Diameter closing removes all dark structures of an image with
maximal extension smaller than diameter_threshold. The maximal
extension is defined as the maximal extension of the bounding box.
The operator is also called Bounding Box Closing. In practice,
the result is similar to a morphological closing, but long and thin
structures are not removed.

Technically, this operator is based on the max-tree representation of
the image.

Parameters
----------
image: ndarray
    The input image for which the diameter_closing is to be calculated.
    This image can be of any type.
diameter_threshold: unsigned int
    The maximal extension parameter (number of pixels). The default value
    is 8.
connectivity: unsigned int, optional
    The neighborhood connectivity. The integer represents the maximum
    number of orthogonal steps to reach a neighbor. In 2D, it is 1 for
    a 4-neighborhood and 2 for a 8-neighborhood. Default value is 1.
parent: ndarray, int64, optional
    Precomputed parent image representing the max tree of the inverted
    image. This function is fast, if precomputed parent and tree_traverser
    are provided. See Note for further details.
tree_traverser: 1D array, int64, optional
    Precomputed traverser, where the pixels are ordered such that every
    pixel is preceded by its parent (except for the root which has no
    parent). This function is fast, if precomputed parent and
    tree_traverser are provided. See Note for further details.

Returns
-------
output: ndarray
    Output image of the same shape and type as input image.

See also
--------
skimage.morphology.area_opening
skimage.morphology.area_closing
skimage.morphology.diameter_opening
skimage.morphology.max_tree

References
----------
.. [1] Walter, T., & Klein, J.-C. (2002). Automatic Detection of
       Microaneurysms in Color Fundus Images of the Human Retina by Means
       of the Bounding Box Closing. In A. Colosimo, P. Sirabella,
       A. Giuliani (Eds.), Medical Data Analysis. Lecture Notes in Computer
       Science, vol 2526, pp. 210-220. Springer Berlin Heidelberg.
       :DOI:`10.1007/3-540-36104-9_23`
.. [2] Carlinet, E., & Geraud, T. (2014). A Comparative Review of
       Component Tree Computation Algorithms. IEEE Transactions on Image
       Processing, 23(9), 3885-3895.
       :DOI:`10.1109/TIP.2014.2336551`

Examples
--------
We create an image (quadratic function with a minimum in the center and
4 additional local minima.

>>> w = 12
>>> x, y = np.mgrid[0:w,0:w]
>>> f = 180 + 0.2*((x - w/2)**2 + (y-w/2)**2)
>>> f[2:3,1:5] = 160; f[2:4,9:11] = 140; f[9:11,2:4] = 120
>>> f[9:10,9:11] = 100; f[10,10] = 100
>>> f = f.astype(np.int)

We can calculate the diameter closing:

>>> closed = diameter_closing(f, 3, connectivity=1)

All small minima with a maximal extension of 2 or less are removed.
The remaining minima have all a maximal extension of at least 3.


Notes
-----
If a max-tree representation (parent and tree_traverser) are given to the
function, they must be calculated from the inverted image for this
function, i.e.:
>>> P, S = max_tree(invert(f))
>>> closed = diameter_closing(f, 3, parent=P, tree_traverser=S)
"""
function diameter_closing!(output::AbstractArray{<:Any, N},
                           image::AbstractArray{<:Any, N};
                           min_diameter::Number=8, connectivity::Integer=1,
                           maxtree::Union{MaxTree{N}, Nothing} = nothing) where N
    _maxtree = isnothing(maxtree) ? MaxTree(image, connectivity=connectivity, rev=true) : maxtree
    return direct_filter!(output, image, _maxtree, diameters(_maxtree), min_diameter)
end

diameter_closing(image::AbstractArray; kwargs...) =
    diameter_closing!(similar(image), image; kwargs...)

"""
Determine all local maxima of the image.

The local maxima are defined as connected sets of pixels with equal
gray level strictly greater than the gray levels of all pixels in direct
neighborhood of the set. The function labels the local maxima.

Technically, the implementation is based on the max-tree representation
of an image. The function is very efficient if the max-tree representation
has already been computed. Otherwise, it is preferable to use
the function local_maxima.

Parameters
----------
image : ndarray
    The input image for which the maxima are to be calculated.
connectivity: unsigned int, optional
    The neighborhood connectivity. The integer represents the maximum
    number of orthogonal steps to reach a neighbor. In 2D, it is 1 for
    a 4-neighborhood and 2 for a 8-neighborhood. Default value is 1.
parent: ndarray, int64, optional
    The value of each pixel is the index of its parent in the ravelled
    array.
tree_traverser: 1D array, int64, optional
    The ordered pixel indices (referring to the ravelled array). The pixels
    are ordered such that every pixel is preceded by its parent (except for
    the root which has no parent).

Returns
-------
local_max : ndarray, uint64
    Labeled local maxima of the image.

See also
--------
skimage.morphology.local_maxima
skimage.morphology.max_tree

References
----------
.. [1] Vincent L., Proc. "Grayscale area openings and closings,
       their efficient implementation and applications",
       EURASIP Workshop on Mathematical Morphology and its
       Applications to Signal Processing, Barcelona, Spain, pp.22-27,
       May 1993.
.. [2] Soille, P., "Morphological Image Analysis: Principles and
       Applications" (Chapter 6), 2nd edition (2003), ISBN 3540429883.
       :DOI:`10.1007/978-3-662-05088-0`
.. [3] Salembier, P., Oliveras, A., & Garrido, L. (1998). Antiextensive
       Connected Operators for Image and Sequence Processing.
       IEEE Transactions on Image Processing, 7(4), 555-570.
       :DOI:`10.1109/83.663500`
.. [4] Najman, L., & Couprie, M. (2006). Building the component tree in
       quasi-linear time. IEEE Transactions on Image Processing, 15(11),
       3531-3539.
       :DOI:`10.1109/TIP.2006.877518`
.. [5] Carlinet, E., & Geraud, T. (2014). A Comparative Review of
       Component Tree Computation Algorithms. IEEE Transactions on Image
       Processing, 23(9), 3885-3895.
       :DOI:`10.1109/TIP.2014.2336551`

Examples
--------
We create an image (quadratic function with a maximum in the center and
4 additional constant maxima.

>>> w = 10
>>> x, y = np.mgrid[0:w,0:w]
>>> f = 20 - 0.2*((x - w/2)**2 + (y-w/2)**2)
>>> f[2:4,2:4] = 40; f[2:4,7:9] = 60; f[7:9,2:4] = 80; f[7:9,7:9] = 100
>>> f = f.astype(np.int)

We can calculate all local maxima:

>>> maxima = max_tree_local_maxima(f)

The resulting image contains the labeled local maxima.
"""
function local_maxima!(output::AbstractArray{<:Any, N},
                       image::AbstractArray{<:Any, N};
                       connectivity::Integer=1,
                       maxtree::Union{MaxTree{N}, Nothing} = nothing) where N
    _maxtree = isnothing(maxtree) ? MaxTree(image_inv, connectivity=connectivity) : maxtree
    return local_maxima!(output, image, maxtree)
end
