using Documenter
using ImageMorphology
using ImageBase, ImageShow, TestImages

prettyurls = get(ENV, "CI", nothing) == "true"
format = Documenter.HTML(; prettyurls)

#! format: off
pages = Any[
    "index.md",
    "Concepts" => Any["structuring_element.md"],
    "reference.md"
]
#! format: on

makedocs(; modules=[ImageMorphology], format=format, sitename="ImageMorphology", pages)

deploydocs(; repo="github.com/JuliaImages/ImageMorphology.jl.git")
