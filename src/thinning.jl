abstract type ThinAlgo end

doc"""
    struct GuoAlgo <: ThinAlgo end

A thinning operation algorithm.

The procedure carries on until there are no changes between two consecutives iterations.

The algorithm is described in the the page 361 (Algorithm A1) of:
* Guo, Z., & Hall, R. W. (1989). Parallel thinning with two-subiteration algorithms. Communications of the ACM, 32(3), 359-373.
"""
struct GuoAlgo <: ThinAlgo end

"""
```
function thinning(img::AbstractArray{Bool}; algo::ThinAlgo=GuoAlgo())
```
Applies a binary blob thinning operation to achieve a skeletization of the input image.

See also:
* [`GuoAlgo`](@ref)
"""
function thinning(img::AbstractArray{Bool}; algo::ThinAlgo=GuoAlgo())
    # dispatch appropriate implementation
    thinning_impl(img, algo)
end

function thinning_impl(img::AbstractArray{Bool}, algo::GuoAlgo) 
    prev = falses(size(img))
    processed = copy(img)
    it = 0
    while prev != processed
        prev = copy(processed)
        it += 1
        thinning_iteration!(processed, isodd(it))
    end
    return processed
end

"""
```
function thinning_iteration(img::AbstractArray{Bool,2}; odd_iteration::Bool)
```
The thining iteration evaluates three conditions in order to determine which pixels of the image should be removed from img.
The three conditions are explained in the page 361 of:
* Guo, Z., & Hall, R. W. (1989). Parallel thinning with two-subiteration algorithms. Communications of the ACM, 32(3), 359-373.
"""
function thinning_iteration!(img_ori::AbstractArray{Bool,2}, odd_iteration::Bool)
    img = falses(size(img_ori).+2)
    img[2:end-1,2:end-1] = img_ori
    h, w = size(img)
    for j=2:w-1, i=2:h-1
        !img[i,j] && continue
        p1 = img[i-1,j-1]
        p2 = img[i-1,j]
        p3 = img[i-1,j+1]
        p4 = img[i,j+1]
        p5 = img[i+1,j+1]
        p6 = img[i+1,j]
        p7 = img[i+1,j-1]
        p8 = img[i,j-1]
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
            img_ori[i-1,j-1] = false
        end
    end
end
