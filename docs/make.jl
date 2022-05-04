using Documenter
using ImageMorphology
using ImageCore, ImageShow, TestImages

prettyurls = get(ENV, "CI", nothing) == "true"
format = Documenter.HTML(; prettyurls)

pages = ["index.md", "reference.md"]
makedocs(; modules=[ImageMorphology], format=format, sitename="ImageMorphology", pages=[])

deploydocs(; repo="github.com/JuliaImages/ImageMorphology.jl.git")
