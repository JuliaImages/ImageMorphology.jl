using Documenter
using ImageMorphology
using OffsetArrays
using ImageBase, ImageShow, TestImages

prettyurls = get(ENV, "CI", nothing) == "true"
format = Documenter.HTML(; prettyurls)

pages = Any[
    "index.md",
    "Concepts" => Any["structuring_element.md"],
    "reference.md"
]

makedocs(;
    modules=[ImageMorphology, OffsetArrays],
    format=format,
    sitename="ImageMorphology",
    pages,
    doctest=false
)

deploydocs(; repo="github.com/JuliaImages/ImageMorphology.jl.git")
