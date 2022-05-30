#=
Optimized version of `@. out = f.(A, B)` for a few things:

- there are some overhead when applying broadcasting on small subarray
  (https://github.com/JuliaLang/julia/issues/28126)
- Apply LoopVectorization when possible
=#
function _mapf!(f, out, A, B)
    use_broadcast = !all(ArrayInterface.fast_scalar_indexing, (out, A, B))
    if use_broadcast
        @. out = f(A, B)
    else
        __mapf_loop!(f, out, A, B)
    end
    return out
end

function __mapf_loop!(f, out, A, B)
    # LoopVectorization currently does not support Colorant
    if ArrayInterface.can_avx(f) && _is_gray(out, A, B)
        __mapf_loop_avx!(f, _numeric_array(out), _numeric_array(A), _numeric_array(B))
    else
        __mapf_loop_simd!(f, out, A, B)
    end
    return out
end

function __mapf_loop_avx!(f, out, A, B)
    @tturbo for i in eachindex(out, A, B)
        out[i] = f(A[i], B[i])
    end
    return out
end

# TODO(johnnychen94): remove these when VectorizationBase lower bound >= v0.21.35?
# https://github.com/JuliaSIMD/VectorizationBase.jl/pull/85
__mapf_loop_avx!(::typeof(max), out::AbstractArray{Bool}, A, B) = __mapf_loop_simd!(|, out, A, B)
__mapf_loop_avx!(::typeof(min), out::AbstractArray{Bool}, A, B) = __mapf_loop_simd!(&, out, A, B)

function __mapf_loop_simd!(f, out, A, B)
    @inbounds @simd for i in eachindex(out, A, B)
        out[i] = f(A[i], B[i])
    end
    return out
end

_numeric_array(A::AbstractArray{T}) where {T<:Real} = A
_numeric_array(A::AbstractArray{AT}) where {AT<:Gray} = reinterpret(eltype(AT), A)
_is_gray(As::AbstractArray...) = all(A -> eltype(A) <: Union{Real,Gray}, As)
