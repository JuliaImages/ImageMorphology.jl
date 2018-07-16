"""
    ThinAlgo

A thinning operation algorithm.
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
    thinning_impl(img, algo)
end

function thinning_impl(img::AbstractArray{Bool}, algo::GuoAlgo)
    # pad input image
    prev = copy(img)
    curr = copy(img)

    it = 1
    thinning_iteration!(curr, isodd(it))
    while prev != curr
        it += 1
        prev .= curr
        thinning_iteration!(curr, isodd(it))
    end

    curr
end

function thinning_iteration!(img::AbstractArray{Bool,2}, odd_iteration::Bool)
    # pad input image with zeros
    pad = falses(size(img).+2)
    pad[2:end-1,2:end-1] = img
    h, w = size(pad)
    for j=2:w-1, i=2:h-1
        !pad[i,j] && continue
        p1 = pad[i-1,j-1]
        p2 = pad[i-1,j]
        p3 = pad[i-1,j+1]
        p4 = pad[i,j+1]
        p5 = pad[i+1,j+1]
        p6 = pad[i+1,j]
        p7 = pad[i+1,j-1]
        p8 = pad[i,j-1]
        C = (!p2 && (p3 || p4)) + (!p4 && (p5 || p6)) + (!p6 && (p7 || p8)) + (!p8 && (p1 || p2))
        N1 = (p1 || p2) + (p3 || p4) + (p5 || p6) + (p7 || p8)
        N2 = (p2 || p3) + (p4 || p5) + (p6 || p7) + (p8 || p1)
        N = min(N1, N2)
        if odd_iteration
            O = (p2 || p3 || (!p5)) && p4
        else
            O = (p6 || p7 || (!p1)) && p8
        end
        if C == 1 && (2 ≤ N ≤ 3) && !O
            img[i-1,j-1] = false
        end
    end
end
