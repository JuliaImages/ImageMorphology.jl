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
Applies a binary blob thinning operation, to achieve a skeletization of the input image.

The function transforms a binary blob image into a skeletized form using the technique of Zhang-Suen.
"""

function thinning_iteration(img::AbstractArray, iter::Int64)
    marker = falses(size(img))
    h, w = size(img)
    for i in range(2,h-1)
        for j in range(2,w-1)
            p2 = img[i-1][j]
            p3 = img[i-1][j+1]
            p4 = img[i][j+1]
            p5 = img[i+1][j+1]
            p6 = img[i+1][j]
            p7 = img[i+1][j-1]
            p8 = img[i][j-1]
            p9 = img[i-1][j-1]

            A  = (p2 == 0 && p3 == 1) + (p3 == 0 && p4 == 1) +
                 (p4 == 0 && p5 == 1) + (p5 == 0 && p6 == 1) +
                 (p6 == 0 && p7 == 1) + (p7 == 0 && p8 == 1) +
                 (p8 == 0 && p9 == 1) + (p9 == 0 && p2 == 1)

            B = p2 + p3 + p4 + p5 + p6 + p7 + p8 + p9
            m1 = iter == 0 ? (p2 * p4 * p6) : (p2 * p4 * p8)
            m2 = iter == 0 ? (p4 * p6 * p8) : (p2 * p6 * p8)
            if (A == 1 && (B >= 2 && B <= 6) && m1 == 0 && m2 == 0)
                marker[i,j] = 1
            end
        end
    end
    return img & ~marker            
end    

"""
'''
output = thinning(img)
'''
The function transforms a binary blob image into a skeletized form using the technique of Zhang-Suen.

-  img          = Binary input image
"""

function thinning_iteration(img::AbstractArray{Bool,2}, iter::Int64)
    marker = falses(size(img))
    h, w = size(img)
    for i in range(2,h-2)
        for j in range(2,w-2)
            p2 = img[i-1,j]
            p3 = img[i-1,j+1]
            p4 = img[i,j+1]
            p5 = img[i+1,j+1]
            p6 = img[i+1,j]
            p7 = img[i+1,j-1]
            p8 = img[i,j-1]
            p9 = img[i-1,j-1]

            A  = (p2 == 0 && p3 == 1) + (p3 == 0 && p4 == 1) +
                 (p4 == 0 && p5 == 1) + (p5 == 0 && p6 == 1) +
                 (p6 == 0 && p7 == 1) + (p7 == 0 && p8 == 1) +
                 (p8 == 0 && p9 == 1) + (p9 == 0 && p2 == 1)

            B = p2 + p3 + p4 + p5 + p6 + p7 + p8 + p9
            m1 = iter == 0 ? (p2 * p4 * p6) : (p2 * p4 * p8)
            m2 = iter == 0 ? (p4 * p6 * p8) : (p2 * p6 * p8)
            if (A == 1 && (B >= 2 && B <= 6) && m1 == 0 && m2 == 0)
                marker[i,j] = 1
            end
        end
    end
    return img.&(.~marker)            
end    

function thinning(img::AbstractArray)
    processed = img.>0.5
    prev = falses(size(img))

    while(count(prev.!=processed)>0)
        prev = processed
        processed = thinning_iteration(processed, 0)
        processed = thinning_iteration(processed, 1)
    end
    return processed    
end