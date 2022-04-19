"""
    connectivity_ball(radius; square=true)
    connectivity_ball(radius, Val(N); square=true)

Build the N-dimensional connectivity array by filling the ellipse/ball region with `true`s,
where `true` indicates a connected position to its center. The default dimension `N` is
\$2\$.

The output shape is specified by `radius`. If `radius` is an `Int` number, then the `true`s
region forms a ball, and if `radius` is a tuple of `Int`, then it is typically an ellipse.
The `square` keyword is used to specify if the output array should be square.

!!! info "special case: radius == 1"
    If `radius == 1` then output array will be filled with `true`. For 2-dimensional case,
    this is also known as C8 connectivity, i.e., `connectivity_region(:C8)`.

```jldoctest; setup = :(using ImageMorphology)
julia> using ImageMorphology: connectivity_ball

julia> connectivity_ball(1, Val(2))
3×3 BitMatrix:
 1  1  1
 1  1  1
 1  1  1

julia> connectivity_ball((1, 3)) # ellipse
7×7 BitMatrix:
 0  0  0  0  0  0  0
 0  0  0  0  0  0  0
 0  0  0  1  0  0  0
 1  1  1  1  1  1  1
 0  0  0  1  0  0  0
 0  0  0  0  0  0  0
 0  0  0  0  0  0  0

julia> connectivity_ball((1, 3); square=false) # ellipse
3×7 BitMatrix:
 0  0  0  1  0  0  0
 1  1  1  1  1  1  1
 0  0  0  1  0  0  0
```

See also [`connectivity_region`](@ref ImageMorphology.connectivity_region).
"""
connectivity_ball(radius; kwargs...) = connectivity_ball(radius, Val(2); kwargs...)
function connectivity_ball(r::Int, ::Val{N}; kwargs...) where {N}
    return connectivity_ball(ntuple(_ -> r, N), Val(N); kwargs...)
end
function connectivity_ball(r::NTuple{N,Int}, ::Val{N}; square=true) where {N}
    sz = square ? ntuple(_ -> 2maximum(r) + 1, N) : @. 2r + 1

    # for special case when r==1, it's common to just use `trues(3,3)`
    all(r .== 1) && return trues(sz)

    # compute the ellips region
    Ω = falses(sz)
    c = OffsetArrays.center(Ω)
    @inbounds for i in CartesianIndices(Ω)
        v = mapreduce(+, i.I, c, r) do x, x0, a
            (x - x0)^2 / (a^2)
        end
        if v <= 1 # within the ellipse
            Ω[i] = true
        end
    end
    return Ω
end
# this is not type-inferable, yet we still provide it for convenience
connectivity_ball(radius, n::Int; kwargs...) = connectivity_ball(radius, Val(n); kwargs...)

"""
    connectivity_region(alias::Symbol)
    connectivity_region(alias::Symbol, Val(N))

Return commonly used N-dimensional connectivity region by alias.

```jldoctest; setup = :(using ImageMorphology)
julia> using ImageMorphology: connectivity_region

julia> connectivity_region(:C8) # equivalent to `connectivity_region(8)`
3×3 BitMatrix:
 1  1  1
 1  1  1
 1  1  1
```

Currently supported alias are

| dimension | alias     |
| --------  | --------- |
| 2         | :C4       |
| 2         | :C8       |

See also [`connectivity_ball`](@ref ImageMorphology.connectivity_ball).
"""
connectivity_region(alias::Symbol) = connectivity_region(alias, Val(2))

# for type-stability, always return BitArray{N}
function connectivity_region(alias::Symbol, ::Val{2})
    if alias == :C4
        Ω = [false true false
            true true true
            false true false]
        return BitMatrix(Ω)
    elseif alias == :C8
        return trues(3, 3)
    else
        error("unknown connectivity alias for dimension 2: $alias")
    end
end

connectivity_region(alias::Symbol, ::Val{N}) where {N} = error("unsupported dimension $N")
# this is not type-inferable, yet we still provide it for convenience
connectivity_region(alias, n::Int) = connectivity_region(alias, Val(n))
