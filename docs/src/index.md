# ImageMorphology.jl

This package provides morphology-related functionalities for structure analysis and image processing.

!!! info "setup"
    ImageMorphology is a sub-package of the umbrella package . You can choose to use either `using
    Images` or `using ImageMorphology` for the functionality.

## Installation

Just like all normal Julia packages, you can use Pkg to install it:

```julia
pkg> add ImageMorphology # hit ] to enter Pkg mode
```

## Quickstart

The perhaps quickest way to use this library is to find a set of useful operators from the
[gallery](@ref op_index) and build your own pipeline from it.

For any advanced usage, we recommend you to read the "Concept" part first, for instance, many
morphology operation supports generic ["structuring element"](@ref concept_se).
