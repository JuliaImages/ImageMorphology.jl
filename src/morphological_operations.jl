# Erode and dilate support 3x3 regions only (and higher-dimensional generalizations).

"""
```
imgd = dilate(img, [region])
```

perform a max-filter over nearest-neighbors. The
default is 8-connectivity in 2d, 27-connectivity in 3d, etc. You can specify the
list of dimensions that you want to include in the connectivity, e.g., `region =
[1,2]` would exclude the third dimension from filtering.
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
"""
erode(img::AbstractArray, region=coords_spatial(img)) = erode!(copy(img), region)

dilate!(maxfilt, region=coords_spatial(maxfilt)) = extremefilt!(maxfilt, max, region)
erode!(minfilt, region=coords_spatial(minfilt)) = extremefilt!(minfilt, min, region)
function extremefilt!(A::AbstractArray, select::Function, region=coords_spatial(A))
    inds = indices(A)
    for d = 1:ndims(A)
        if size(A, d) == 1 || !in(d, region)
            continue
        end
        Rpre = CartesianRange(inds[1:d-1])
        Rpost = CartesianRange(inds[d+1:end])
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
"""
opening(img::AbstractArray, region=coords_spatial(img)) = opening!(copy(img), region)
opening!(img::AbstractArray, region=coords_spatial(img)) = dilate!(erode!(img, region),region)

"""
`imgc = closing(img, [region])` performs the `closing` morphology operation, equivalent to `erode(dilate(img))`.
`region` allows you to control the dimensions over which this operation is performed.
"""
closing(img::AbstractArray, region=coords_spatial(img)) = closing!(copy(img), region)
closing!(img::AbstractArray, region=coords_spatial(img)) = erode!(dilate!(img, region),region)

"""
`imgth = tophat(img, [region])` performs `top hat` of an image,
which is defined as the image minus its morphological opening.
`region` allows you to control the dimensions over which this operation is performed.
"""
tophat(img::AbstractArray, region=coords_spatial(img)) = img - opening(img, region)

"""
`imgbh = bothat(img, [region])` performs `bottom hat` of an image,
which is defined as its morphological closing minus the original image.
`region` allows you to control the dimensions over which this operation is performed.
"""
bothat(img::AbstractArray, region=coords_spatial(img)) = closing(img, region) - img

"""
`imgmg = morphogradient(img, [region])` returns morphological gradient of the image,
which is the difference between the dilation and the erosion of a given image.
`region` allows you to control the dimensions over which this operation is performed.
"""
morphogradient(img::AbstractArray, region=coords_spatial(img)) = dilate(img, region) - erode(img, region)

"""
`imgml = morpholaplace(img, [region])` performs `Morphological Laplacian` of an image,
which is defined as the arithmetic difference between the internal and the external gradient.
`region` allows you to control the dimensions over which this operation is performed.
"""
morpholaplace(img::AbstractArray, region=coords_spatial(img)) = dilate(img, region) + erode(img, region) - 2img


"""
```
function thinning_iteration(img::AbstractArray{Bool,2}; odd_iteration::Bool)
```
The thining iteration evaluates three conditions in order to determine which pixels of the image should be removed from img.
The three conditions are explained in:
Guo, Z., & Hall, R. W. (1989). Parallel thinning with two-subiteration algorithms. Communications of the ACM, 32(3), 359-373.
"""
function thinning_iteration!(img_ori::AbstractArray{Bool,2}, odd_iteration::Bool)
    marker = trues(size(img_ori))
    img = falses(size(img_ori).+2)
    img[2:end-1,2:end-1]=img_ori
    h, w = size(img)
    for i=2:h-1
        for j=2:w-1
            if !img[i,j]
                continue
            end
            p1 = img[i-1,j-1]
            p2 = img[i-1,j]
            p3 = img[i-1,j+1]
            p4 = img[i,j+1]
            p5 = img[i+1,j+1]
            p6 = img[i+1,j]
            p7 = img[i+1,j-1]
            p8 = img[i,j-1]
            A = (!p2 && (p3 || p4)) + (!p4 && (p5 || p6)) + (!p6 && (p7 || p8)) + (!p8 && (p1 || p2))
            B1 = (p1 || p2) + (p3 || p4) + (p5 || p6) + (p7 || p8)
            B2 = (p2 || p3) + (p4 || p5) + (p6 || p7) + (p8 || p1)
            B = min(B1,B2)
            G3 = false
            if (odd_iteration)
                G3 = (p2 || p3 || (!p5)) && p4
            else
                G3 = (p6 || p7 || (!p1)) && p8
            end
            if (A == 1) && ((2 <= B) && (B <= 3)) && (!G3) 
                marker[i-1,j-1] =false
            end
        end
        
    end
    img_ori[:] = img_ori[:].&(marker[:])
end

doc"""
```julia
function thinning(img::AbstractArray{Bool}; maxiter::Integer=0) 
```

Applies a binary blob thinning operation, to achieve a skeletization of the input image.
The procedure carries on until max_iterations is achieved or when there are no changes between two consecutives iterations.

The algorithm is described in:
1. Guo, Z., & Hall, R. W. (1989). Parallel thinning with two-subiteration algorithms. Communications of the ACM, 32(3), 359-373.
2. Lam, L., Lee, S. W., & Suen, C. Y. (1992). Thinning methodologies-a comprehensive survey. IEEE Transactions on pattern analysis and machine intelligence, 14(9), 869-885.
"""
function thinning(img::AbstractArray{Bool}; maxiter::Integer=0) 
    prev = falses(size(img))
    processed = copy(img)
    it = 0
    while ((maxiter==0) || (it<maxiter)) && (sum(prev.!=processed)>0)
        prev = copy(processed)
        it = it +1
        thinning_iteration!(processed, isodd(it))
    end
    return processed
end
