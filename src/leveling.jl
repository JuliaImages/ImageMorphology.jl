#utility #TODO move to some core.jl
function check_image_consistency(im1::GenericGrayImage, im2::GenericGrayImage)
    (size(im1) == size(im2)) ||
        throw(DimensionMismatch("The sizes of the two images do not match"))
    (typeof(im1) == typeof(im2)) ||
        throw("The types of the two images do not match")
end

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

import Base.isless
import Base.isgreater
import DataStructures.PriorityQueue, DataStructures.enqueue!, DataStructures.dequeue_pair!

"""
low_leveling(imageref::AbstractArray{T, N}, marker::AbstractArray{T, N}, connectivity::AbstractArray{Bool}) -> L

 A function g is an lower leveling of a function f if and only if for any couple of neighbouring pixels (p, q) gp>gq->gp>=f
 This algorithm modify g to become a low_leveling of f
    
 # Returns
 An array of the same type and shape as the `imageref`.
 
 # See also
 [`low_leveling!`](@ref),[`leveling`](@ref),
 [`underbuild`](@ref),[`overbuild`](@ref)
 
 # References
 1.'From connected operators to levelings' [Meyer,1998] In ISMM '98: Proceedings of the fourth international symposium on Mathematical morphology and its applications to image and signal processingJune 1998 Pages 191–198
 2.'Mathematical Morphology: From Theory to Applications' chap6 [Serra et al.] https://doi.org/10.1002/9781118600788.ch8 
"""
function low_leveling(imageref::AbstractArray{T, N}, marker::AbstractArray{T, N}, connectivity::AbstractArray{Bool}) where {T<:ImageCore.NumberLike , N}
    check_image_consistency(imageref, marker)
    deltaoffsets = neighborhood(connectivity)
    output = copy(marker)
    # declare priority queue, store position of pixel associated with priority 
    pq = PriorityQueue{CartesianIndex{N}, T}()
    R = CartesianIndices(axes(output))
    for i in R
        @inbounds if output[i] < imageref[i]
            @inbounds current_max = output[i]
            for Δi in deltaoffsets # examine neighborhoods
                ii = i + Δi
                if checkbounds(Bool, R, ii) #check that we are in the image
                    @inbounds current_max= max(current_max,output[ii])
                end
            end
            if @inbounds current_max != output[i]
                @inbounds output[i] = min(imageref[i],current_max)
                @inbounds enqueue!(pq, i, output[i])
            end
        # else posponed to the end
        end
    end
    # Loop until all pixel have been examined 
    while !isempty(pq)
        curr_idx, _ = dequeue_pair!(pq)
        for Δi in deltaoffsets # examine neighborhoods
            ii = curr_idx + Δi
            if checkbounds(Bool, R, ii) #check that we are in the image
                @inbounds if output[ii] < output[curr_idx] && output[ii] < imageref[ii]
                    @inbounds new_value= min( imageref[ii] , output[curr_idx])
                    @inbounds if output[ii] != new_value
                        @inbounds output[ii] = new_value
                        pq[ii]=new_value
                    end
                end 
            end
        end
    end
    #Post-processing, so that we have nicer levelings
	#(L-= inf( ref, dilate(mark)) , and L+=sup(ref,erode(mark)),
	#and thus L(ref,mark)=L-( L+(ref,mark),mark) and they commute)
    for i in R
        if output[i] > imageref[i]
            output[i] = imageref[i]
        end
    end
    return output
end

"""
high_leveling(imageref::AbstractArray{T, N}, marker::AbstractArray{T, N}, connectivity::AbstractArray{Bool}) -> L

A function g is an upper leveling of a function f if and only if for any couple of neighbouring pixels (p, q) gp>gq->gp<=f
This algorithm modify g to become a high_leveling of f

# Returns
 An array of the same type and shape as the `imageref`.
 
 # See also
 [`low_leveling!`](@ref),[`leveling`](@ref),
 [`underbuild`](@ref),[`overbuild`](@ref)
 
 # References
 1.'From connected operators to levelings' [Meyer,1998] In ISMM '98: Proceedings of the fourth international symposium on Mathematical morphology and its applications to image and signal processingJune 1998 Pages 191–198
 2.'Mathematical Morphology: From Theory to Applications' chap6 [Serra et al.] https://doi.org/10.1002/9781118600788.ch8 
"""
function high_leveling(imageref::AbstractArray{T, N}, marker::AbstractArray{T, N}, connectivity::AbstractArray{Bool}) where {T<:ImageCore.NumberLike , N}
    check_image_consistency(imageref, marker)
    deltaoffsets = neighborhood(connectivity)
    output = copy(marker)
    # declare priority queue, store position of pixel associated with priority 
    pq = PriorityQueue{CartesianIndex{N}, T}()
    R = CartesianIndices(axes(output))
    for i in R
        @inbounds if output[i] > imageref[i]
            @inbounds current_min = output[i]
            for Δi in deltaoffsets # examine neighborhoods
                ii = i + Δi
                if checkbounds(Bool, R, ii) #check that we are in the image
                    @inbounds current_min= min(current_min,output[ii])
                end
            end
            if @inbounds current_min != output[i]
                @inbounds output[i] = max(imageref[i],current_min)
                @inbounds enqueue!(pq, i, output[i])
            end
        # else posponed to the end
        end
    end
    # Loop until all pixel have been examined 
    while !isempty(pq)
        curr_idx, _ = dequeue_pair!(pq)
        for Δi in deltaoffsets # examine neighborhoods
            ii = curr_idx + Δi
            if checkbounds(Bool, R, ii) #check that we are in the image
                @inbounds if output[ii] > output[curr_idx] && output[ii] > imageref[ii]
                    @inbounds new_value = max( imageref[ii] , output[curr_idx])
                    @inbounds if output[ii] != new_value
                        @inbounds output[ii] = new_value
                        pq[ii] = new_value
                    end
                end 
            end
        end
    end
    #Post-processing, so that we have nicer levelings
	#(L-= inf( ref, dilate(mark)) , and L+=sup(ref,erode(mark)),
	#and thus L(ref,mark)=L-( L+(ref,mark),mark) and they commute)
    for i in R
        if output[i] < imageref[i]
            output[i] = imageref[i]
        end
    end
    return output
end

"""
leveling(imageref::AbstractArray{T, N}, marker::AbstractArray{T, N}, connectivity::AbstractArray{Bool}) -> L

A function g is a leveling of a function f if and only if it is both an upper and a lower leveling of the function f 
This algorithm modify g to become a leveling of f
 
# Returns
 An array of the same type and shape as the `imageref`.
 
 # See also
 [`low_leveling!`](@ref),[`high_leveling`](@ref),
 [`underbuild`](@ref),[`overbuild`](@ref)
 
 # References
 1.'From connected operators to levelings' [Meyer,1998] In ISMM '98: Proceedings of the fourth international symposium on Mathematical morphology and its applications to image and signal processingJune 1998 Pages 191–198
 2.'Mathematical Morphology: From Theory to Applications' chap6 [Serra et al.] https://doi.org/10.1002/9781118600788.ch8 
"""
function leveling(imageref::AbstractArray{T, N}, marker::AbstractArray{T, N}, connectivity::AbstractArray{Bool}) where {T<:ImageCore.NumberLike , N}
    check_image_consistency(imageref, marker)
    tmp = high_leveling(imageref,marker,connectivity)
    return low_leveling(tmp,marker,connectivity)
end