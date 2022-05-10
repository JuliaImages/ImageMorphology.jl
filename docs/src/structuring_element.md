```@setup concept_se
using ImageMorphology
using ImageBase
using TestImages
```

# [Structuring element](@id concept_se)

Structuring Element (SE) is the key concept in morphology to indicate the connectivity and the
neighborhood. This page explains this structuring element concept, and how ImageMorphology supports
the general SEs without compromising the performance on the most commonly used special SE cases.

## The erosion example

The erosion `erode` function in its simplest 1-dimensional case can be defined as

$$\varepsilon_A[p] = min(A[p-1], A[p], A[p+1])$$

Because the output value at position $p$ not only depends on its own pixel `A[p]` but also on
its neighborhood values `A[p-1]` and `A[p+1]`, we call this type of operation a _neighborhood
image transformation_.

Now comes the question: **if we try to generalize the `erode` function, what should we do?** --
we would like to generalize the concept of "neighborhood".

## Two neighborhood representations

By saying "$\Omega_p$ is the neighborhood of $p$", we are expressing `p in Ωₚ` in plain Julia. For
performance consideration, this `Ωₚ` is usually generated from the `(p, Ω)` pair. `p` is the center
point that changes during the iteration, and `Ω` is usually a pre-defined and unchanged data which
contains the neighborhood and shape information. We call this `Ω` a _structuring element_. There are
usually two ways to express `Ω`:

- **displacement offset**: a list of `CartesianIndex` to inidcate the offset to the center point `p`
- **connectivity mask**: a bool array mask to indicate the connectivity to the center point `p`

For instance, in the following code block we build a commonly named C4 connectivity in the 2-dimensional case:

```@example concept_se
# displacement offset
Ω_offsets = [
    CartesianIndex(-1, 0),
    CartesianIndex(0, -1),
    CartesianIndex(0, 1),
    CartesianIndex(1, 0),
]

# connectivity mask
Ω_bool = Bool[
    0 1 0
    1 1 1
    0 1 0
]
nothing #hide
```

If `p=CartesianIndex(3, 3)`, then we know `p=CartesianIndex(3, 4)` is in `Ωₚ`.

Now, back to the erosion example. Based on the displacement offset representation, the simplest
generic version of `erode` can be implemented quite simply:

```@example concept_se
# For illustration only, performance can be greatly improved using iteration to eliminate allocations
function my_erode(A, Ω)
    out = similar(A)
    R = CartesianIndices(A)
    for p in R
        Ωₚ = filter!(q->in(q, R), Ref(p) .+ Ω)
        # here we don't assume p in Ωₚ
        out[p] = min(A[p], minimum(A[Ωₚ]))
    end
    return out
end
nothing #hide
```

```@example concept_se
using ImageMorphology
using ImageBase
using TestImages

img = Gray.(testimage("morphology_test_512"))
img = Gray.(img .< 0.8)
img_e = my_erode(img, Ω_offsets)
mosaic(img, img_e; nrow=1)
```

As you may realize, the displacement offset representation is convenient to use when implementing
algorithms, but it is hard to visualize. In contrast, the connectivity mask is not so convenient to
use when implementing algorithms, but it is easy to visualize. For instance, one can very easily
understand the following SE at the first glance:

```@example concept_se
Ω = Bool[1 1 1; 1 1 0; 0 0 0] # hide
```

but not

```@example concept_se
strel(CartesianIndex, Ω) # hide
```

## The `strel` function

This package supports the conversion between different SE representations via the [`strel`](@ref)
helper function. `strel` is the short name for "STRucturing ELement".

To convert a connectivity mask representation to displacement offset representation:

```@example concept_se
Ω_mask = Bool[1 1 1; 1 1 0; 0 0 0] |> centered
Ω_offsets = strel(CartesianIndex, Ω_mask)
```

!!! note "zero-centered mask"
    The mask array is expected to be zero-centered. That means, the axes of a 3×3 mask `axes(se)`
    should be `(-1:1, -1:1)`. The [`centered`](@ref OffsetArrays.centered) function is used to shift
    the center point of the array to `(0, 0, ..., 0)`.
    ```julia
    julia> A = centered([1 2 3; 4 5 6; 7 8 9])
    3×3 OffsetArray(::Matrix{Int64}, -1:1, -1:1) with eltype Int64 with indices -1:1×-1:1:
     1  2  3
     4  5  6
     7  8  9

    julia> A[-1, -1], A[0, 0], A[1, 1] # top-left, center, bottom-right
    (1, 5, 9)
    ```
    This `centered` function comes from [OffsetArrays.jl](https://github.com/JuliaArrays/OffsetArrays.jl)
    and is also exported by ImageMorphology.


And to convert back from a displacement offset representation to connectivity mask representation:

```@example concept_se
strel(Bool, Ω_offsets)
```

Quite simple, right? Thus to make our `my_erode` function more generic, we only need to add one
single line:

```diff
 function my_erode(A, Ω)
     out = similar(A)
+    Ω = strel(CartesianIndex, Ω)
     R = CartesianIndices(A)
```

## Convenient constructors

Among all the SE possibilities, this package provides constructors for two commonly used cases:

- diamond-like constructor: [`strel_diamond`](@ref)
- box-like constructor: [`strel_box`](@ref)

```@repl concept_se
strel_diamond((3, 3)) # immediate neighborhood: C4 connectivity
strel_diamond((3, 3), (1, )) # along the first dimension
strel_box((3, 3)) # all adjacent neighborhood: C8 connectivity
strel_box((3, 3), (1, ))
```

Utilizing these constructors, we can provide an easier-to-use `my_erode(A, [dims])` interface by
adding one more method:

```julia
my_erode(A, dims::Dims=ntuple(identity, ndims(A))) = my_erode(A, strel_diamond(A, dims))
```

!!! tip "Performance tip: keep the array type"
    For the structuring element `Ω` generated from `strel_diamond` and `strel_box`, it is likely
    to hit a fast path if you keep its array type. For instance, `erode(A, strel_diamond(A))` is
    usually faster than `erode(A, Array(strel_diamond(A)))` because more information of the `Ω`
    shape is passed to Julia during coding and compilation.

## Performance optimizations and the `strel_type` function

Thanks to Julia's multiple dispatch mechanism, we can provide all the optimization tricks without
compromising the simple user interface. This can be programmatically done with the help of the
`strel_type` function. For example, if you know a very efficient `erode` implementation for the C4
connectivity SE, then you can add it incrementally:

```julia
using ImageMorphology: MorphologySE, SEDiamond

my_erode(A, dims::Dims) = my_erode(A, strel_diamond(A, dims))
my_erode(A, Ω) = _my_erode(strel_type(Ω), A, Ω)

# the generic implementation we've written above
function _my_erode(::MorphologySE, A, Ω)
   ...
end

# the optimized implementation for SE generated from `strel_diamond` function
function _my_erode(::SEDiamond, A, Ω)
   ...
end

# ... and other optimized versions, if there are
```

In essence, `strel_type` is a trait function to assist the dispatch and code design:

```@repl concept_se
strel_type(Ω_mask)
```

It returns an internal object `SEMask{2}()`. This might look scary at first glance, but it's quite a
simple lookup table that reflects our previous reasoning:

| representation          | element type     | `strel_type` |
| ----------------------- | ---------------- | ------------ |
| displacement offset     | `CartesianIndex` | `SEOffset`   |
| connectivity mask       | `Bool`           | `SEMask`     |
| [`strel_diamond`](@ref) | `Bool`           | `SEDiamond`  |
| [`strel_box`](@ref)     | `Bool`           | `SEBox`      |
