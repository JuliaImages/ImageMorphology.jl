"""
    ThinAlgo

A thinning algorithm.
"""
abstract type ThinAlgo end

"""
    struct GuoAlgo <: ThinAlgo end

The Guo algorithm evaluates three conditions in order to determine which pixels of
the image should be removed.

The three conditions are explained in the page 361 of **Guo, Z., & Hall, R. W. (1989).
Parallel thinning with two-subiteration algorithms. Communications of the ACM, 32(3), 359-373.**
"""
struct GuoAlgo <: ThinAlgo end

"""
    thinning(img::AbstractArray{Bool}; algo::ThinAlgo=GuoAlgo())

Applies a binary blob thinning operation to achieve a skeletization of the input image.

See also: [`GuoAlgo`](@ref)
"""
function thinning(img::AbstractArray{Bool}; algo::ThinAlgo=GuoAlgo())
    # dispatch appropriate implementation
    return thinning_impl(img, algo)
end

function thinning_impl(img::AbstractArray{Bool,2}, algo::GuoAlgo)
    # pad input image
    prev = falses(size(img) .+ 2)
    prev[2:(end - 1), 2:(end - 1)] = img

    # preallocate memory for the mask
    mask = falses(size(prev))

    # perform single iteration
    it = 1
    curr = copy(prev)
    guo_iteration!(mask, curr, isodd(it))
    curr[mask] .= false

    # loop if necessary
    while prev != curr
        it += 1
        prev .= curr
        guo_iteration!(mask, curr, isodd(it))
        curr[mask] .= false
    end

    # unpad result
    return curr[2:(end - 1), 2:(end - 1)]
end

# update mask with Guo iteration on padded image
function guo_iteration!(mask::AbstractArray{Bool,2}, img::AbstractArray{Bool,2}, odd_iteration::Bool)
    # start clean
    mask .= false

    # loop over pixels and update mask
    @inbounds for j in 2:(size(img, 2) - 1), i in 2:(size(img, 1) - 1)
        img[i, j] || continue
        p1 = img[i - 1, j - 1]
        p2 = img[i - 1, j]
        p3 = img[i - 1, j + 1]
        p4 = img[i, j + 1]
        p5 = img[i + 1, j + 1]
        p6 = img[i + 1, j]
        p7 = img[i + 1, j - 1]
        p8 = img[i, j - 1]
        C = (!p2 && (p3 || p4)) + (!p4 && (p5 || p6)) + (!p6 && (p7 || p8)) + (!p8 && (p1 || p2))
        N1 = (p1 || p2) + (p3 || p4) + (p5 || p6) + (p7 || p8)
        N2 = (p2 || p3) + (p4 || p5) + (p6 || p7) + (p8 || p1)
        N = min(N1, N2)
        if odd_iteration
            O = (p2 || p3 || !p5) && p4
        else
            O = (p6 || p7 || !p1) && p8
        end
        if C == 1 && (2 ≤ N ≤ 3) && !O
            mask[i, j] = true
        end
    end
end
