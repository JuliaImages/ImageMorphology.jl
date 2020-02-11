# ImageMorphology

[![][travis-img]][travis-url]
[![][appveyor-img]][appveyor-url]
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
dilate
erode
opening
closing
tophat
bothat
morphogradient
morpholaplace
thinning
imfill
```

## Documentation

Please check the top-level documentation at [Images.jl][images-url].

<!-- URLS -->

[pkgeval-img]: https://juliaci.github.io/NanosoldierReports/pkgeval_badges/I/ImageMorphology.svg
[pkgeval-url]: https://juliaci.github.io/NanosoldierReports/pkgeval_badges/report.html
[travis-img]: https://travis-ci.org/JuliaImages/ImageMorphology.jl.svg?branch=master
[travis-url]: https://travis-ci.org/JuliaImages/ImageMorphology.jl
[appveyor-img]: https://ci.appveyor.com/api/projects/status/github/JuliaImages/ImageMorphology.jl?svg=true
[appveyor-url]: https://ci.appveyor.com/project/JuliaImages/ImageMorphology-jl
[codecov-img]: https://codecov.io/github/JuliaImages/ImageMorphology.jl/coverage.svg?branch=master
[codecov-url]: https://codecov.io/github/JuliaImages/ImageMorphology.jl?branch=master

[images-url]: https://github.com/JuliaImages/Images.jl
