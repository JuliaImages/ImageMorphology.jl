
#utility
function check_image_consistency(im1::GenericGrayImage, im2::GenericGrayImage)
    (size(im1) == size(im2)) ||
        throw(DimensionMismatch("The sizes of the two images do not match"))
    (typeof(im1) == typeof(im2)) ||
        throw("The types of the two images do not match")
end

import Base.isless
import Base.isgreater
import DataStructures.Queue, DataStructures.Deque, DataStructures.enqueue!, DataStructures.dequeue!, DataStructures.dequeue_pair!


"""
    hmaxima(image, connectivity=[], constant) -> Array
The h-maxima transformation suppresses all maxima whose "relative" depth is below a given threshold (constant)
# Arguments
- `image::AbstractArray{T, N}`: where {T<:ImageCore.NumberLike , N} the ``N``-dimensional input image
- `connectivity::AbstractArray{Bool}`: the neighborhood connectivity.
- `constant::T`: the soustracted constant.
# Returns
An array of the same type and shape as the `image`.
# See also
[`hminima`](@ref)
# References
    Morphological image analysis by Soille pg 170-172
"""
function hmaxima(image::AbstractArray{T, N}, connectivity::AbstractArray{Bool}, constant::T) where {T<:ImageCore.NumberLike , N}
    tmp = image.-constant
    return underbuild(tmp,image,connectivity) 
end

"""
    hminima(image, connectivity=[], constant) -> Array
The h-minima transformation suppresses all minima whose "relative" depth is below a given threshold (constant)
# Arguments
- `image::AbstractArray{T, N}`: where {T<:ImageCore.NumberLike , N} the ``N``-dimensional input image
- `connectivity::AbstractArray{Bool}`: the neighborhood connectivity.
- `constant::T`: the added constant.
# Returns
An array of the same type and shape as the `image`.
# See also
[`hmaxima`](@ref)
# References
    Morphological image analysis by Soille pg 170-172
"""
function hminima(image::AbstractArray{T, N}, connectivity::AbstractArray{Bool}, constant::T) where {T<:ImageCore.NumberLike , N}
    tmp = image.+constant
    return overbuild(tmp,image,connectivity)    
end

"""
    underbuild(marker, mask, connectivity=[]) -> Array

Morphological reconstruction by dilatation according to connection type

# Details
This operation performs dilations of marker restricted by the provided mask until idempotence This means that the output image
will be included between mask and min(marker,mask).
# Arguments
- `marker::AbstractArray{T, N}`: where {T<:ImageCore.NumberLike , N} the ``N``-dimensional input image
- `mask::AbstractArray{T, N}`: where {T<:ImageCore.NumberLike , N} the ``N``-dimensional input image
- `connectivity::AbstractArray{Bool}`: the neighborhood connectivity.
# Returns
An array of the same type and shape as the `marker`.
# See also
[`overbuild`](@ref)
# References
  Morphological Grayscale Reconstruction in Image Analysis: Applications and Efficient Algorithms IEEE Trans Image Process. 1993;2(2)
"""
function underbuild(marker::AbstractArray{T, N}, mask::AbstractArray{T, N}, connectivity::AbstractArray{Bool}) where {T<:ImageCore.NumberLike , N}
    check_image_consistency(mask, marker)   
   
    all(in((1,3)), size(connectivity)) || throw(ArgumentError("connectivity must have size 1 or 3 in each dimension"))
    for d = 1:ndims(connectivity)
        size(connectivity, d) == 1 || reverse(connectivity; dims=d) == connectivity || throw(ArgumentError("connectivity must be symmetric"))
    end

    # WARNING
    #  The original algorithm assume that the marker < mask, this assertion is not always true
    #  so instead we extract the start from min(marker,mask)
    # first operation take the infimum between marker and mask
    output = map(min, marker, mask)

    queue = Queue{CartesianIndex{N}}(); 
    deltaoffsets = neighborhood(connectivity)
    uper_deltaoffest = upper_neighborhood(connectivity)
    lower_deltaoffest= lower_neighborhood(connectivity)

    #NOTE we could pad the array with -Inf to speedup forward/backward scan 
    # forward scan
    R = CartesianIndices(axes(marker))
    for i in R
        @inbounds current_max = output[i]
        for Δi in uper_deltaoffest # examine neighborhoods
            ii = i + Δi
            if checkbounds(Bool, R, ii) #check that we are in the image
                @inbounds current_max = max(current_max,output[ii])
            end
        end
        @inbounds output[i] = min(current_max,mask[i])
    end
    # backward scan
    for i in reverse(R)
        @inbounds current_max = output[i]
        for Δi in lower_deltaoffest # examine neighborhoods
            ii = i + Δi
            if checkbounds(Bool, R, ii) #check that we are in the image
                @inbounds current_max = max(current_max,output[ii])
            end
        end
        @inbounds output[i] = min(current_max,mask[i])
        for Δi in lower_deltaoffest # examine neighborhoods
            ii = i + Δi
            if checkbounds(Bool, R, ii) #check that we are in the image
                @inbounds if output[ii]<output[i] && output[ii]<mask[ii]
                    enqueue!(queue, i)
                end
            end
        end
    end
    # Loop until all pixel have been examined 
    while !isempty(queue)
        curr_idx = dequeue!(queue)
        for Δi in deltaoffsets # examine neighborhoods
            ii = curr_idx + Δi
            if checkbounds(Bool, R, ii) #check that we are in the image
                @inbounds if output[ii] < output[curr_idx] && mask[ii] != output[ii] 
                    output[ii] = min(output[curr_idx],mask[ii])
                    enqueue!(queue, ii)
                end
            end
        end
    end
    return output
end

#specialization for binary case
function underbuild(marker::AbstractArray{<:Union{Bool,AbstractGray{Bool}},N}, mask::AbstractArray{<:Union{Bool,AbstractGray{Bool}},N}, connectivity::AbstractArray{Bool})  where {N}
    check_image_consistency(mask, marker)
    all(in((1,3)), size(connectivity)) || throw(ArgumentError("connectivity must have size 1 or 3 in each dimension"))
    for d = 1:ndims(connectivity)
        size(connectivity, d) == 1 || reverse(connectivity; dims=d) == connectivity || throw(ArgumentError("connectivity must be symmetric"))
    end

    output = copy(marker)

    #Use queue to store propagation front
    propagationfront = Queue{CartesianIndex{N}}(); 
    deltaoffsets = neighborhood(connectivity)

    # Initialisation of the queue with contour pixels of marker image
    R = CartesianIndices(axes(output))
    for i in R
        @inbounds if output[i] == 0 continue
        end
        for Δi in deltaoffsets # examine neighborhoods
            ii = i + Δi
            if checkbounds(Bool, R, ii) #check that we are in the image #note here we could add border or simulate then by splitting iteration
                #Ok we walk on the perimeter
                @inbounds if output[ii] == 0 && mask[ii] == 1
                    #push pixel in the queue
                    enqueue!(propagationfront, i)
                    break
                end
            end
        end
    end
    # Loop until all pixel have been examined 
    while !isempty(propagationfront)
        curr_idx = dequeue!(propagationfront)
        for Δi in deltaoffsets # examine neighborhoods
            ii = curr_idx + Δi
            if checkbounds(Bool, R, ii) #check that we are in the image
                #propagate
                @inbounds if output[ii] == 0 && mask[ii] == 1
                    @inbounds output[ii] = 1
                    enqueue!(propagationfront, ii)
                end
            end
        end
    end
    return output
end
"""
    overbuild(marker, mask, connectivity=[]) -> Array

Morphological reconstruction by erosion according to connection type

# Details
This operation performs erosions of marker restricted by the provided mask until idempotence This means that the output image
will be included between mask and max(marker,mask).
# Arguments
- `marker::AbstractArray{T, N}`: where {T<:ImageCore.NumberLike , N} the ``N``-dimensional input image
- `mask::AbstractArray{T, N}`: where {T<:ImageCore.NumberLike , N} the ``N``-dimensional input image
- `connectivity::AbstractArray{Bool}`: the neighborhood connectivity.
# Returns
An array of the same type and shape as the `marker`.
# See also
[`underbuild`](@ref)
# References
  Morphological Grayscale Reconstruction in Image Analysis: Applications and Efficient Algorithms IEEE Trans Image Process. 1993;2(2)
"""
function overbuild(marker::AbstractArray{T, N}, mask::AbstractArray{T, N}, connectivity::AbstractArray{Bool}) where {T<:ImageCore.NumberLike , N}
    check_image_consistency(mask, marker)
    
    all(in((1,3)), size(connectivity)) || throw(ArgumentError("connectivity must have size 1 or 3 in each dimension"))
    for d = 1:ndims(connectivity)
        size(connectivity, d) == 1 || reverse(connectivity; dims=d) == connectivity || throw(ArgumentError("connectivity must be symmetric"))
    end
   
    # WARNING
    #  The original algorithm assume that the marker > mask, this assertion is not always true
    #  so instead we extract the start from max(marker,mask)
    # first operation take the supremum between marker and mask
    output = map(max, marker, mask)

    queue = Queue{CartesianIndex{N}}(); 
    deltaoffsets = neighborhood(connectivity)
    uper_deltaoffest = upper_neighborhood(connectivity)
    lower_deltaoffest= lower_neighborhood(connectivity)

    #NOTE we could pad the array with +Inf to speedup forward/backward scan 
    # forward scan
    R = CartesianIndices(axes(marker))
    for i in R
        @inbounds current_min = output[i]
        for Δi in uper_deltaoffest # examine neighborhoods
            ii = i + Δi
            if checkbounds(Bool, R, ii) #check that we are in the image
                @inbounds current_min = min(current_min,output[ii])
            end
        end
        @inbounds output[i] = max(current_min,mask[i])
    end
    # backward scan
    for i in reverse(R)
        @inbounds current_min = output[i]
        for Δi in lower_deltaoffest # examine neighborhoods
            ii = i + Δi
            if checkbounds(Bool, R, ii) #check that we are in the image
                @inbounds current_min = min(current_min,output[ii])
            end
        end
        @inbounds output[i] = max(current_min,mask[i])
        for Δi in lower_deltaoffest # examine neighborhoods
            ii = i + Δi
            if checkbounds(Bool, R, ii) #check that we are in the image
                @inbounds if output[ii]>output[i] && output[ii]>mask[ii]
                    enqueue!(queue, i)
                end
            end
        end
    end
    # Loop until all pixel have been examined 
    while !isempty(queue)
        curr_idx = dequeue!(queue)
        for Δi in deltaoffsets # examine neighborhoods
            ii = curr_idx + Δi
            if checkbounds(Bool, R, ii) #check that we are in the image
                @inbounds if output[ii] > output[curr_idx] && mask[ii] != output[ii] 
                    output[ii] = max(output[curr_idx],mask[ii])
                    enqueue!(queue, ii)
                end
            end
        end
    end
    return output
end

#specialization for binary case
function overbuild(marker::AbstractArray{<:Union{Bool,AbstractGray{Bool}},N}, mask::AbstractArray{<:Union{Bool,AbstractGray{Bool}},N}, connectivity::AbstractArray{Bool})  where {N}
    marker_inv =  map(!, marker)
    mask_inv = map(!, mask)
    output_inv= underbuild(marker_inv,mask_inv,connectivity)
    return map(!, output_inv)
end

"""
    regional_maxima(image, connectivity=[]) -> Array

Determines all *regional maxima* of the `image`.

# Details
The *regional maximum* is defined as the connected set of pixels that have the
same value, which is greater than the values of all pixels in direct
neighborhood of the set.
This implementation is faster than maxtree approach if maxtree not precomputed
# Arguments
- `image::AbstractArray{T, N}`: where {T<:ImageCore.NumberLike , N} the ``N``-dimensional input image
- `connectivity::AbstractArray{Bool}`: the neighborhood connectivity.
"""
regional_maxima(image::AbstractArray{T, N}, connectivity::AbstractArray{Bool}) where {T<:ImageCore.NumberLike , N} = extract_regional_extrema(image, connectivity, isgreater ) 

"""
    regional_minima(image, connectivity=[]) -> Array

Determines all *regional minima* of the `image`.

# Details
The *regional minima* is defined as the connected set of pixels that have the
same value, which is lower than the values of all pixels in direct
neighborhood of the set.

This implementation is faster than maxtree approach if maxtree not precomputed

# Arguments
- `image::AbstractArray{T, N}`: where {T<:ImageCore.NumberLike , N} the ``N``-dimensional input image
- `connectivity::AbstractArray{Bool}`: the neighborhood connectivity.
"""
regional_minima(image::AbstractArray{T, N}, connectivity::AbstractArray{Bool}) where {T<:ImageCore.NumberLike , N} = extract_regional_extrema(image, connectivity, isless ) 


neighborhood(connectivity) = append!(upper_neighborhood(connectivity), lower_neighborhood(connectivity))

upper_neighborhood(x) = _neighborhood(i->x[i], CartesianIndices(x))
lower_neighborhood(x) = _neighborhood(i->x[i], reverse(CartesianIndices(x)))

function _neighborhood(f, R)
    center = CartesianIndex(map(axes(R)) do ax
        (first(ax) + last(ax)) ÷ 2
    end)
    out = CartesianIndex[]
    for i in R
        i == center && break
        if f(i)
            push!(out,i - center)
        end
    end
    return out
end


#faster than local_minima/maxima from maxtree
function extract_regional_extrema(image::AbstractArray{T, N}, connectivity::AbstractArray{Bool}, comp) where {T<:ImageCore.NumberLike , N}
    #restrict connectivity in order to operator make sens
    all(in((1,3)), size(connectivity)) || throw(ArgumentError("connectivity must have size 1 or 3 in each dimension"))
    for d = 1:ndims(connectivity)
        size(connectivity, d) == 1 || reverse(connectivity; dims=d) == connectivity || throw(ArgumentError("connectivity must be symmetric"))
    end
    center = CartesianIndex(map(axes(connectivity)) do ax
        (first(ax) + last(ax)) ÷ 2
    end)
    output = similar(Array{Bool}, axes(image))
    fill!(output, false)
    visited = similar(Array{Bool}, axes(image))
    fill!(visited, false)
    current_region = Deque{CartesianIndex{N}}() # use deque to store point belonging to current region
    propagationfront = Queue{CartesianIndex{N}}(); #Use queue to store propagation front
    is_region_regional_extremum = 0
    deltaoffsets = neighborhood(connectivity)
    R = CartesianIndices(axes(image))
    for i in R
        # compute the flat zones of the image
        @inbounds if !visited[i] #candidate point aka not visited yet, start flat zone
            # mark point as visited
            @inbounds visited[i] = true
            push!(current_region, i)
            enqueue!(propagationfront, i)
            #default behavior our current region (flat zone) is a extremum until we could proove its not
            is_region_regional_extremum = true
            while !isempty(propagationfront)
                #get indice from propagationfront
                idx_front = dequeue!(propagationfront)
                # get value at current point
                @inbounds centervalue = image[idx_front]
                #examine neighborhood points
                for Δi in deltaoffsets
                    ii = idx_front + Δi
                    if checkbounds(Bool, R, ii) #check that we are in the image
                        @inbounds nl_value = image[ii]
                        if (is_region_regional_extremum && comp(nl_value,centervalue))
                            # no, we are not an extremum anymore, mark is as such and
                            # don't propagate the flat zone on this pixel
                            is_region_regional_extremum = false
                        else
                            # is this pixel in the flat zone ?
                            if nl_value == centervalue
                                # yes, propagate the flat zone to this pixel
                                # (Note: we must always propagate on the *whole* flat zone, so no quick exit even if
                                #  we are not a regional extremum )
                                @inbounds if !visited[ii]
                                    enqueue!(propagationfront, ii)
                                    # mark it as visited
                                    visited[ii] = true
                                    # if we are still in a regional extremum, add this pixel to region
                                    if is_region_regional_extremum
                                        push!(current_region, ii)
                                    end
                                end
                            end
                        end
                    end
                end #for all neighborValue
            end # propagation front empty
            #ok a region is alive, its a regional extremum ?
            if is_region_regional_extremum
                #process the region 
                # note that we could label instead of simply mark regions
                for p in current_region
                    @inbounds output[p] = true
                end
                empty!(current_region)  
            else
                empty!(current_region)  
            end
        end # test we can visit current point
    end # loop along idx image
    return output
end

