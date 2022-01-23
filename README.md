# ImageMorphology

[![][action-img]][action-url]
[![][pkgeval-img]][pkgeval-url]
[![][codecov-img]][codecov-url]

This package provides morphology-related functionality to the [Images.jl][images-url] project.

## Installation

Get the latest stable release with Julia's package manager:

```julia
Pkg.add("ImageMorphology")
```

## Exported functions

```julia
area_closing
area_closing!
area_opening
area_opening!
areas
bothat
boundingboxes
clearborder
closing
closing!
component_boxes
component_centroids
component_indices
component_lengths
component_subscripts
convexhull
diameter_closing
diameter_closing!
diameter_opening
diameter_opening!
diameters
dilate
dilate!
distance_transform
erode
erode!
feature_transform
imfill
isboundary
isboundary!
label_components
label_components!
local_maxima
local_maxima!
local_minima
local_minima!
morphogradient
morpholaplace
opening
opening!
thinning
tophat
```

## Documentation

Please check the top-level documentation at [Images.jl][images-url].

<!-- URLS -->

[pkgeval-img]: https://juliaci.github.io/NanosoldierReports/pkgeval_badges/I/ImageMorphology.svg
[pkgeval-url]: https://juliaci.github.io/NanosoldierReports/pkgeval_badges/report.html
[action-img]: https://github.com/JuliaImages/ImageMorphology.jl/workflows/Unit%20test/badge.svg
[action-url]: https://action-ci.org/JuliaImages/ImageMorphology.jl
[codecov-img]: https://codecov.io/github/JuliaImages/ImageMorphology.jl/coverage.svg?branch=master
[codecov-url]: https://codecov.io/github/JuliaImages/ImageMorphology.jl?branch=master

[images-url]: https://github.com/JuliaImages/Images.jl
