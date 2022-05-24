#! format: off
@deprecate label_components(A::AbstractArray, region::Union{Dims, AbstractVector{Int}}, bkg = 0) label_components(A; bkg=bkg, dims=region)
@deprecate label_components(A::AbstractArray{T,N}, connectivity::Array{Bool,N}, bkg) where {T,N}   label_components(A, connectivity; bkg=bkg)
@deprecate label_components!(out::AbstractArray{Int}, A::AbstractArray, region::Union{Dims, AbstractVector{Int}}, bkg = 0)  label_components!(out, A; bkg=bkg, dims=region)
@deprecate label_components!(out::AbstractArray{Int,N}, A::AbstractArray{T,N}, connectivity::Array{Bool,N}, bkg) where {T,N}  label_components!(out, A, connectivity; bkg=bkg)

@deprecate imfill(img::AbstractArray{Bool}, interval::Tuple{Real,Real}, dims::Union{Dims, AbstractVector{Int}}) imfill(img, interval; dims=dims)

@deprecate dilate!(img; kwargs...) dilate!(img, copy(img); kwargs...)
@deprecate erode!(img; kwargs...) erode!(img, copy(img); kwargs...)

@deprecate opening!(img; kwargs...) opening!(img, copy(img), similar(img); kwargs...)
@deprecate closing!(img; kwargs...) closing!(img, copy(img), similar(img); kwargs...)

Base.@deprecate_binding morphogradient mgradient
Base.@deprecate_binding morpholaplace mlaplacian

@deprecate extremefilt!(A, select; kwargs...) extreme_filter!(select, A, copy(A); kwargs...)

import .FeatureTransform: feature_transform
@deprecate feature_transform(img, weights; kwargs...) feature_transform(img; weights=weights, kwargs...)
