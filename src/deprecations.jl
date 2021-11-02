@deprecate label_components(A::AbstractArray, region::Union{Dims, AbstractVector{Int}}, bkg = 0) label_components(A; bkg=bkg, dims=region)
@deprecate label_components(A::AbstractArray{T,N}, connectivity::Array{Bool,N}, bkg) where {T,N}   label_components(A, connectivity; bkg=bkg)
@deprecate label_components!(out::AbstractArray{Int}, A::AbstractArray, region::Union{Dims, AbstractVector{Int}}, bkg = 0)  label_components!(out, A; bkg=bkg, dims=region)
@deprecate label_components!(out::AbstractArray{Int,N}, A::AbstractArray{T,N}, connectivity::Array{Bool,N}, bkg) where {T,N}  label_components!(out, A, connectivity; bkg=bkg)

@deprecate imfill(img::AbstractArray{Bool}, interval::Tuple{Real,Real}, dims::Union{Dims, AbstractVector{Int}}) imfill(img, interval; dims=dims)

@deprecate dilate(img::AbstractArray, region)   dilate(img; dims=region)
@deprecate dilate!(img::AbstractArray, region)  dilate!(img; dims=region)
@deprecate erode(img::AbstractArray, region)    erode(img; dims=region)
@deprecate erode!(img::AbstractArray, region)   erode!(img; dims=region)

@deprecate opening(img::AbstractArray, region)  opening(img; dims=region)
@deprecate opening!(img::AbstractArray, region) opening!(img; dims=region)
@deprecate closing(img::AbstractArray, region)  closing(img; dims=region)
@deprecate closing!(img::AbstractArray, region) closing!(img; dims=region)

@deprecate tophat(img::AbstractArray, region)   tophat(img; dims=region)
@deprecate bothat(img::AbstractArray, region)   bothat(img; dims=region)

@deprecate morphogradient(img::AbstractArray, region)   morphogradient(img; dims=region)
@deprecate morpholaplace(img::AbstractArray, region)    morpholaplace(img; dims=region)

import .FeatureTransform: feature_transform
@deprecate feature_transform(img, weights; kwargs...) feature_transform(img; weights=weights, kwargs...)
