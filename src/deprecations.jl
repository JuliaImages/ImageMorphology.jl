@deprecate label_components(A::AbstractArray, region::Union{Dims, AbstractVector{Int64}}, bkg = 0) label_components(A; bkg=bkg, dims=region)
@deprecate label_components(A::AbstractArray{T,N}, connectivity::Array{Bool,N}, bkg) where {T,N}   label_components(A, connectivity; bkg=bkg)
@deprecate label_components!(out::AbstractArray{Int}, A::AbstractArray, region::Union{Dims, AbstractVector{Int64}}, bkg = 0)  label_components!(out, A; bkg=bkg, dims=region)
@deprecate label_components!(out::AbstractArray{Int,N}, A::AbstractArray{T,N}, connectivity::Array{Bool,N}, bkg) where {T,N}  label_components!(out, A, connectivity; bkg=bkg)

@deprecate imfill(img::AbstractArray{Bool}, interval::Tuple{Real,Real}, dims::Union{Dims, AbstractVector{Int}}) imfill(img, interval; dims=dims)
