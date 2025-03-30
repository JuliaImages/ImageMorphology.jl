module DiskArraysExt

using ImageMorphology, DiskArrays

get_extreme_value(::typeof(max), ::Type{T}) where {T} = typemin(T)
get_extreme_value(::typeof(min), ::Type{T}) where {T} = typemax(T)

function ImageMorphology.extreme_filter!(f, out, A::AbstractDiskArray, Ω)
    axes(out) == axes(A) || throw(DimensionMismatch("axes(out) must match axes(A)"))
    ImageMorphology.require_select_function(f, eltype(A))
    val0 = get_extreme_value(f, eltype(A))

    pad = extrema.(axes(Ω))
    PAD1 = CartesianIndex(ntuple(i->pad[i][1], length(pad)))
    PAD2 = CartesianIndex(ntuple(i->pad[i][2], length(pad)))
    ZERO = CartesianIndex(ntuple(_->0, ndims(A)))
    ONE = CartesianIndex(ntuple(_->1, ndims(A)))
    SIZEA = CartesianIndex(size(A))
    chunk_pad = [c+length(range(p...))-1
                 for (c,p) in zip(DiskArrays.max_chunksize(eachchunk(A)),pad)]
    A_chunk = Array{eltype(A)}(undef, chunk_pad...)
    out_chunk = similar(A_chunk)
    for i in eachchunk(A)
        I = minimum(CartesianIndices(i))
        R = maximum(CartesianIndices(i))
        A_chunk .= val0
        A1 = max(ONE, I+PAD1)
        A2 = min(SIZEA, R+PAD2)
        A_chunk[ONE:A2-A1+ONE] .= A[A1:A2]
        ImageMorphology._extreme_filter!(strel_type(Ω), f, out_chunk, A_chunk, Ω)
        OUT1 = ONE-PAD1 + min(ZERO, I+PAD1-ONE)
        OUT2 = OUT1+R-I
        out[I:R] .= out_chunk[OUT1:OUT2]
    end
    return out
end

end
