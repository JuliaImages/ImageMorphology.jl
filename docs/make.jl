using Documenter
using DemoCards
using ImageMorphology
using OffsetArrays
using ImageBase, ImageShow, TestImages

operators, operators_cb, demo_assets = makedemos("operators")

prettyurls = get(ENV, "CI", nothing) == "true"
assets = Any[demo_assets]

pages = Any[
    "index.md",
    "Concepts" => Any["structuring_element.md"],
    operators,
    "reference.md"
]

makedocs(;
    modules=[ImageMorphology, OffsetArrays],
    format=Documenter.HTML(; prettyurls, assets, size_threshold=300_000_000),  # mreconstruct has a GIF
    sitename="ImageMorphology",
    pages,
	checkdocs=:exports,
    doctest=false
)
operators_cb()

deploydocs(;
    forcepush = true,
    push_preview = true,
    repo="github.com/JuliaImages/ImageMorphology.jl.git"
)
