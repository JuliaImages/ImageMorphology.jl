# [Morphology Operators](@id op_index)

This page describes the interfaces for the basic morphology operators that you can use to build your
own pipeline. To easily visualize the effect of different operators, the "blobs" image is used to
build the cover image.

```@example
using TestImages, ImageBase
restrict(testimage("blobs"))
```

!!! note "color image is not supported"
    Unless explicitly used, most operations in this package don't support color images. This is
    because many morphological operation is based on pixel value comparisons `min`, `max`, `>` and
    `<` -- they are not well-defined for `RGB` type.

!!! note "🚧 work in progress"
    This page contains only a subset of exported functions. The missing ones will be added in the
    future. You might need to check the heavy [reference](@ref reference_index) page to find out
    what you need. Contributions are welcome!

{{{democards}}}
