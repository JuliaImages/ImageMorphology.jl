module ImageMorphologyLoopVectorizationExt

using ImageMorphology, LoopVectorization

using ImageCore
using OffsetArrays


function __init__()
    # let the system know that we can opt-in to use LoopVectorization accelerations by default
    ImageMorphology.USE_SIMD[] = true
end

@static if Sys.WORD_SIZE == 64
    # optimized `extreme_filter` implementation for SEDiamond SE and 2D images
    # Similar approach could be found in
    # Žlaus, Danijel & Mongus, Domen. (2018). In-place SIMD Accelerated Mathematical Morphology. 76-79. 10.1145/3206157.3206176.
    # see https://www.researchgate.net/publication/325480366_In-place_SIMD_Accelerated_Mathematical_Morphology
    function _extreme_filter_diamond_2D!(f::ImageMorphology.MAX_OR_MIN, out::AbstractArray{T}, A, iter) where {T}
        @debug "call the AVX-enabled `extreme_filter` implementation for 2D SEDiamond" fname = _extreme_filter_diamond_2D!
        out_actual = out

        out = OffsetArrays.no_offset_view(out)
        A = OffsetArrays.no_offset_view(A)

        src = if (out === A) || (iter > 1)
            # To avoid the result affected by loop order, we need two arrays
            T.(A)
        elseif eltype(A) != T
            # To avoid incorrect result, we need to ensure the eltypes are the same
            # https://github.com/JuliaImages/ImageMorphology.jl/issues/104
            T.(A)
        else
            A
        end
        # NOTE(johnnychen94): we don't need to do `out .= src` here because it is write-only;
        # all values are generated from the read-only `src`.

        # LoopVectorization currently doesn't understand Gray and N0f8 types, thus we
        # reinterpret to its raw data type UInt8
        src = rawview(channelview(src))
        out = rawview(channelview(out))

        # creating temporaries
        ySize, xSize = size(A)
        tmp = similar(out, eltype(out), ySize)
        tmp2 = similar(tmp)
        tmp3 = similar(tmp)

        # applying radius=r filter is equivalent to applying radius=1 filter r times
        for i in 1:iter
            # NOTE(johnnychen94): this implementation is essentially equivalent to the
            # following loop. Here we buffer the getindex to further improve the performance
            # by explicitly introducing SIMD buffers.

            # for y in 2:(size(src, 2) - 1), x in 2:(size(src, 1) - 1)
            #     tmp = f(f(src[x, y], src[x - 1, y]), src[x + 1, y])
            #     out[x, y] = f(f(tmp, src[x, y - 1]), src[x, y + 1])
            # end

            #compute first edge column
            #translate to clipped connection, x neighborhood of interest, ? neighborhood we don't care, . center
            # x ?
            # . x
            # x ?
            viewprevious = view(src, :, 1)
            # dilate/erode col 1
            _unsafe_shift_arith!(f, tmp2, tmp, tmp3, viewprevious)
            viewnext = view(src, :, 2)
            viewout = view(out, :, 1)
            # inf/sup between dilate/erode col 1 0 and col 2
            LoopVectorization.vmap!(f, viewout, viewnext, tmp2)
            #next->current
            viewcurrent = view(src, :, 2)
            for c in 3:xSize
                # viewprevious c-2
                # viewcurrent  c-1
                # viewnext     c
                # ? x ?
                # x . x
                # ? x ?
                viewout = view(out, :, c - 1)
                viewnext = view(src, :, c)
                # dilate(x-1)/erode(x-1)
                _unsafe_shift_arith!(f, tmp2, tmp, tmp3, viewcurrent)
                @turbo warn_check_args = false for i in eachindex(viewout, viewprevious, tmp2)
                    #sup(x-2,dilate(x-1)),inf(x-2,erode(x-1))
                    #sup(sup(x-2,dilate(x-1),x) || inf(inf(x-2,erode(x-1),x)
                    viewout[i] = f(f(viewprevious[i], tmp2[i]), viewnext[i])
                end
                #current->previous
                viewprevious = view(src, :, c - 1)
                #next->current
                viewcurrent = view(src, :, c)
            end
            #end last column
            #translate to clipped connection
            # ? x
            # x .
            # ? x
            # dilate/erode col x
            viewout = view(out, :, xSize)
            _unsafe_shift_arith!(f, tmp2, tmp, tmp3, viewcurrent)
            LoopVectorization.vmap!(f, viewout, tmp2, viewprevious)
            if iter > 1 && i < iter
                src .= out
            end
        end

        return out_actual
    end
else
    # FIXME(johnnychen94): unknown segfault happens to 32bit machine
    # https://github.com/JuliaImages/ImageMorphology.jl/pull/106
    function _extreme_filter_diamond_2D!(f::ImageMorphology.MAX_OR_MIN, out::AbstractArray, A, iter)
        return ImageMorphology._extreme_filter_diamond_generic!(f, out, A, strel_diamond(A; r=iter))
    end
end

# optimized implementation for SEDiamond -- a typical case of separable filter
function ImageMorphology._extreme_filter_diamond!(f, out, A, Ω::ImageMorphology.SEDiamondArray{N}, ::Val{true}) where {N}
    rΩ = strel_size(Ω) .÷ 2
    if N == 2 && f isa ImageMorphology.MAX_OR_MIN && all(rΩ[1] .== rΩ)
        iter = rΩ[1]
        return _extreme_filter_diamond_2D!(f, out, A, iter)
    else
        return ImageMorphology._extreme_filter_diamond_generic!(f, out, A, Ω)
    end
end

function ImageMorphology._extreme_filter_box!(f, out, A, Ω::ImageMorphology.SEBoxArray{N}, ::Val{true}) where {N}
    rΩ = strel_size(Ω) .÷ 2
    if N == 2 && f isa ImageMorphology.MAX_OR_MIN && all(rΩ[1] .== rΩ)
        iter = rΩ[1]
        return _extreme_filter_box_2D!(f, out, A, iter)
    else
        return ImageMorphology._extreme_filter_generic!(f, out, A, Ω)
    end
end

@static if Sys.WORD_SIZE == 64
    # optimized `extreme_filter` implementation for SEBox SE and 2D images
    # Similar approach could be found in
    # Žlaus, Danijel & Mongus, Domen. (2018). In-place SIMD Accelerated Mathematical Morphology. 76-79. 10.1145/3206157.3206176.
    # see https://www.researchgate.net/publication/325480366_In-place_SIMD_Accelerated_Mathematical_Morphology
    function _extreme_filter_box_2D!(f::ImageMorphology.MAX_OR_MIN, out::AbstractArray{T}, A, iter) where {T}
        @debug "call the AVX-enabled `extreme_filter` implementation for 2D SEBox" fname = _extreme_filter_box_2D!
        out_actual = out

        out = OffsetArrays.no_offset_view(out)
        A = OffsetArrays.no_offset_view(A)

        src = if (out === A) || (iter > 1)
            # To avoid the result affected by loop order, we need two arrays
            T.(A)
        elseif eltype(A) != T
            # To avoid incorrect result, we need to ensure the eltypes are the same
            # https://github.com/JuliaImages/ImageMorphology.jl/issues/104
            T.(A)
        else
            A
        end
        # NOTE(johnnychen94): we don't need to do `out .= src` here because it is write-only;
        # all values are generated from the read-only `src`.

        # LoopVectorization currently doesn't understand Gray and N0f8 types, thus we
        # reinterpret to its raw data type UInt8
        src = rawview(channelview(src))
        out = rawview(channelview(out))

        #creating temporaries
        ySize, xSize = size(A)
        tmp = similar(out, eltype(out), ySize)
        tmp2 = similar(tmp)
        tmp3 = similar(tmp)
        tmp4 = similar(tmp)
        tmp5 = similar(tmp)

        # applying radius=r filter is equivalent to applying radius=1 filter r times
        for i in 1:iter
            #compute first edge column
            #translate to clipped connection, x neighborhood of interest, ? neighborhood we don't care, . center
            # x x
            # . x
            # x x
            viewprevious = view(src, :, 1)
            # dilate/erode col 1
            _unsafe_shift_arith!(f, tmp2, tmp, tmp5, viewprevious)
            viewnext = view(src, :, 2)
            viewout = view(out, :, 1)
            # dilate/erode col 2
            _unsafe_shift_arith!(f, tmp3, tmp, tmp5, viewnext)
            #sup(dilate col 1,dilate col 0) or inf(erode col 1,erode col 0)
            LoopVectorization.vmap!(f, viewout, tmp3, tmp2)
            #next->current
            viewcurrent = view(src, :, 2)
            for c in 3:xSize
                # Invariant of the loop
                # temp2 contains dilation/erosion of previous col x-2
                # temp3 contains dilation/erosion of current col x-1
                # temp4 contains dilation/erosion of next col ie x
                # viewprevious c-2
                # viewcurrent  c-1
                # viewnext     c
                # x x x
                # x . x
                # x x x
                viewout = view(out, :, c - 1)
                viewnext = view(src, :, c)
                # dilate(x)/erode(x)
                _unsafe_shift_arith!(f, tmp4, tmp, tmp5, viewnext)
                @turbo warn_check_args = false for i in eachindex(viewout, tmp4, tmp3, tmp2)
                    #sup(x-2,dilate(x-1)),inf(x-2,erode(x-1))
                    #sup(dilate(x-2),dilate(x-1),dilate(y)) or inf(erode(x-2),erode(x-1),erode(x))
                    viewout[i] = f(f(tmp2[i], tmp3[i]), tmp4[i])
                end
                #swap
                tmp4, tmp3 = tmp3, tmp4
                tmp4, tmp2 = tmp2, tmp4
            end
            #end last column
            #translate to clipped connection
            # x x
            # x .
            # x x
            # dilate/erode col x
            viewout = view(out, :, xSize)
            _unsafe_shift_arith!(f, tmp4, tmp, tmp5, viewcurrent)
            LoopVectorization.vmap!(f, viewout, tmp2, tmp3)
            if iter > 1 && i < iter
                src .= out
            end
        end

        return out_actual
    end
else
    # FIXME(johnnychen94): unknown segfault happens to 32bit machine
    # https://github.com/JuliaImages/ImageMorphology.jl/pull/106
    function _extreme_filter_box_2D!(f::ImageMorphology.MAX_OR_MIN, out::AbstractArray, A, iter)
        return ImageMorphology._extreme_filter_generic!(f, out, A, strel_box(A; r=iter))
    end
end

# ptr level optimized implementation for Real types
# short-circuit all check
# call LoopVectorization directly
function _unsafe_shift_arith!(
    f::ImageMorphology.MAX_OR_MIN,
    out::AbstractVector,
    tmp::AbstractVector,
    tmp2::AbstractVector,
    A::AbstractVector,
) #tmp external to reuse external allocation
    if f === min
        padd = typemax(eltype(out))
    else
        padd = typemin(eltype(out))
    end
    N = length(out)
    #in  src  = [1,2,3,4], padd, N=4
    #out tmp  = [padd,1,2,3]
    #out tmp2 = [1,2,3,padd]
    _unsafe_padded_copyto!(tmp, A, true, N, padd)
    _unsafe_padded_copyto!(tmp2, A, false, N, padd)
    @turbo warn_check_args = false for i in eachindex(out, A, tmp, tmp2)
        out[i] = f(f(A[i], tmp[i]), tmp2[i])
    end
    return out
end

# Shift vector of size N up or down by 1 accorded to dir, pad with v
const SafeArrayTypes{T,N,NA} = Union{Array{T,N},Base.SubArray{T,N,Array{T,NA}}}
function _unsafe_padded_copyto!(dest::SafeArrayTypes{T,1}, src::SafeArrayTypes{T,1}, dir, N, v) where {T}
    # ptr level optimized implementation for Real types
    # short-circuit all check so unsafe
    if dir
        unsafe_store!(pointer(dest), v)
        unsafe_copyto!(pointer(dest, 2), pointer(src, 1), N - 1)
    else
        unsafe_copyto!(pointer(dest, 1), pointer(src, 2), N - 1)
        unsafe_store!(pointer(dest), v, N)
    end
    return dest
end
function _unsafe_padded_copyto!(dest::AbstractVector, src::AbstractVector, dir, N, v)
    if dir
        dest[begin] = v
        copyto!(dest, 2, src, 1, N - 1)
    else
        dest[end] = v
        copyto!(dest, 1, src, 2, N - 1)
    end
    return dest
end


end