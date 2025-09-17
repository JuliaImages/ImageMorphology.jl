import Base.isless
import Base.isgreater
import DataStructures.PriorityQueue, DataStructures.dequeue_pair!

"""
    low_leveling(op, marker, mask; [dims])
    low_leveling(op, marker, mask, se)

A function g is an lower leveling of a function f if and only if 
for any couple of neighbouring pixels (p, q) g(p)>g(q)->g(p)>=f
This algorithm modify marker to become a high_leveling of ref

The `dims` keyword is used to specify the dimension to process by constructing the box shape
structuring element [`strel_box(marker; dims)`](@ref strel_box). For generic structuring
element, the half-size is expected to be either `0` or `1` along each dimension.

# See also
 [`high_leveling`](@ref),[`leveling`](@ref),
 [`underbuild`](@ref),[`overbuild`](@ref)
 
# References
- [1] 'From connected operators to levelings' [Meyer,1998] In ISMM '98: Proceedings of the fourth international symposium on Mathematical morphology and its applications to image and signal processingJune 1998 Pages 191–198
- [2] 'Mathematical Morphology: From Theory to Applications' chap6 [Serra et al.] https://doi.org/10.1002/9781118600788.ch8 
"""
function low_leveling(ref, marker; dims=coords_spatial(marker))
    return low_leveling(ref, marker, strel_box(marker, dims))
end
function low_leveling(ref, marker, se)
    return low_leveling!(similar(ref), ref, marker, se)
end

function low_leveling!(out, ref, marker; dims=coords_spatial(marker))
    return low_leveling!(out, ref, marker, strel_box(marker, dims))
end
function low_leveling!(out, ref, marker, se)
    return _low_leveling!(out, ref, marker, se)
end
function _low_leveling!(out::AbstractArray{T,N}, ref::AbstractArray{T,N}, marker::AbstractArray{T,N}, se) where {T<:ImageCore.NumberLike,N}
    axes(out) == axes(ref) == axes(marker) || throw(DimensionMismatch("images should have the same axes"))

    # For generic structuring element, the half-size should to be either `0` or `1` along
    # each dimension. See also section 6.2.3 of [2].
    se_size = strel_size(se)
    if length(se_size) != N
        msg = "the input structuring element is not for $N dimensional array, instead it is for $(length(se_size)) dimensional array"
        throw(DimensionMismatch(msg))
    end
    if !all(x -> in(x, (1, 3)), strel_size(se))
        msg = "structuring element with half-size larger than 1 is invalid"
        throw(DimensionMismatch(msg))
    end
    se = strel(CartesianIndex, se)
    out = copy(marker)

    # declare priority queue, store position of pixel associated with priority 
    pq = PriorityQueue{CartesianIndex{N},T}()
    R = CartesianIndices(axes(out))
    for i in R
        if @inbounds out[i] < ref[i]
            current_max = out[i]
            for Δi in se # examine neighborhoods
                ii = i + Δi
                if checkbounds(Bool, R, ii) #check that we are in the image
                    @inbounds current_max = max(current_max, out[ii])
                end
            end
            if @inbounds current_max != out[i]
                out[i] = min(ref[i], current_max)
                push!(pq, i, out[i])
            end
            # else posponed to the end
        end
    end
    # Loop until all pixel have been examined 
    while !isempty(pq)
        curr_idx, _ = dequeue_pair!(pq)
        for Δi in se # examine neighborhoods
            ii = curr_idx + Δi
            if checkbounds(Bool, R, ii) #check that we are in the image
                @inbounds if out[ii] < out[curr_idx] && out[ii] < ref[ii]
                    new_value = min(ref[ii], out[curr_idx])
                    if out[ii] != new_value
                        out[ii] = new_value
                        pq[ii] = new_value
                    end
                end
            end
        end
    end
    #Post-processing, so that we have nicer levelings
    #(L-= inf( ref, dilate(mark)) , and L+=sup(ref,erode(mark)),
    #and thus L(ref,mark)=L-( L+(ref,mark),mark) and they commute)
    out .= min.(out, ref)
    return out
end

"""
    high_leveling(op, marker, mask; [dims])
    high_leveling(op, marker, mask, se)

A function g is an upper leveling of a function f if and only if for any couple of 
neighbouring pixels (p, q) g(p)>g(q)->g(p)<=f
This algorithm modify marker to become a high_leveling of ref

The `dims` keyword is used to specify the dimension to process by constructing the box shape
structuring element [`strel_box(marker; dims)`](@ref strel_box). For generic structuring
element, the half-size is expected to be either `0` or `1` along each dimension.

# See also
 [`low_leveling`](@ref),[`leveling`](@ref),
 [`underbuild`](@ref),[`overbuild`](@ref)
 
# References
- [1] 'From connected operators to levelings' [Meyer,1998] In ISMM '98: Proceedings of the fourth international symposium on Mathematical morphology and its applications to image and signal processingJune 1998 Pages 191–198
- [2] 'Mathematical Morphology: From Theory to Applications' chap6 [Serra et al.] https://doi.org/10.1002/9781118600788.ch8 
"""
function high_leveling(ref, marker; dims=coords_spatial(marker))
    return high_leveling(ref, marker, strel_box(marker, dims))
end
function high_leveling(ref, marker, se)
    return high_leveling!(similar(ref), ref, marker, se)
end

function high_leveling!(out, ref, marker; dims=coords_spatial(marker))
    return high_leveling!(out, ref, marker, strel_box(marker, dims))
end
function high_leveling!(out, ref, marker, se)
    return _high_leveling!(out, ref, marker, se)
end
function _high_leveling!(out::AbstractArray{T,N}, ref::AbstractArray{T,N}, marker::AbstractArray{T,N}, se) where {T<:ImageCore.NumberLike,N}
    axes(out) == axes(ref) == axes(marker) || throw(DimensionMismatch("images should have the same axes"))

    # For generic structuring element, the half-size should to be either `0` or `1` along
    # each dimension. See also section 6.2.3 of [2].
    se_size = strel_size(se)
    if length(se_size) != N
        msg = "the input structuring element is not for $N dimensional array, instead it is for $(length(se_size)) dimensional array"
        throw(DimensionMismatch(msg))
    end
    if !all(x -> in(x, (1, 3)), strel_size(se))
        msg = "structuring element with half-size larger than 1 is invalid"
        throw(DimensionMismatch(msg))
    end

    se = strel(CartesianIndex, se)
    out = copy(marker)

    # declare priority queue, store position of pixel associated with priority 
    pq = PriorityQueue{CartesianIndex{N},T}()
    R = CartesianIndices(axes(out))
    for i in R
        @inbounds if out[i] > ref[i]
            current_min = out[i]
            for Δi in se # examine neighborhoods
                ii = i + Δi
                if checkbounds(Bool, R, ii) #check that we are in the image
                    @inbounds current_min = min(current_min, out[ii])
                end
            end
            if @inbounds current_min != out[i]
                out[i] = max(ref[i], current_min)
                push!(pq, i, out[i])
            end
            # else posponed to the end
        end
    end
    # Loop until all pixel have been examined 
    while !isempty(pq)
        curr_idx, _ = dequeue_pair!(pq)
        for Δi in se # examine neighborhoods
            ii = curr_idx + Δi
            if checkbounds(Bool, R, ii) #check that we are in the image
                @inbounds if out[ii] > out[curr_idx] && out[ii] > ref[ii]
                    new_value = max(ref[ii], out[curr_idx])
                    if out[ii] != new_value
                        out[ii] = new_value
                        pq[ii] = new_value
                    end
                end
            end
        end
    end
    #Post-processing, so that we have nicer levelings
    #(L-= inf( ref, dilate(mark)) , and L+=sup(ref,erode(mark)),
    #and thus L(ref,mark)=L-( L+(ref,mark),mark) and they commute)
    out .= max.(out, ref)
    return out
end

"""
leveling(ref, marker)

A function g is a leveling of a function f if and only if it is both an upper and a lower leveling of the function f 
This algorithm modify g (a copy of the marker) to become a leveling of f (the ref)

 # See also
 [`low_leveling](@ref),[`high_leveling`](@ref),
 [`underbuild`](@ref),[`overbuild`](@ref)
 
 # References
 1.'From connected operators to levelings' [Meyer,1998] In ISMM '98: Proceedings of the fourth international symposium on Mathematical morphology and its applications to image and signal processingJune 1998 Pages 191–198
 2.'Mathematical Morphology: From Theory to Applications' chap6 [Serra et al.] https://doi.org/10.1002/9781118600788.ch8 
"""
function leveling(ref, marker; dims=coords_spatial(marker))
    return leveling(ref, marker, strel_box(marker, dims))
end
function leveling(ref, marker, se)
    return leveling!(similar(ref), ref, marker, se)
end

function leveling!(out, ref, marker; dims=coords_spatial(marker))
    return leveling!(out, ref, marker, strel_box(marker, dims))
end
function leveling!(out, ref, marker, se)
    return _leveling!(out, ref, marker, se)
end

function _leveling!(out, ref, marker, se)
    tmp = _high_leveling!(out, ref, marker, se)
    return _low_leveling!(out, tmp, marker, se)
end
