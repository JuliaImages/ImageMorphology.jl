# ImageMorphology

[![Build Status](https://travis-ci.org/JuliaImages/ImageMorphology.jl.svg?branch=master)](https://travis-ci.org/JuliaImages/ImageMorphology.jl)
[![Build status](https://ci.appveyor.com/api/projects/status/github/JuliaImages/ImageMorphology.jl?branch=master&svg=true)](https://ci.appveyor.com/project/kmsquire/imagemorphology-jl/branch/master)
[![codecov.io](http://codecov.io/github/JuliaImages/ImageMorphology.jl/coverage.svg?branch=master)](http://codecov.io/github/JuliaImages/ImageMorphology.jl?branch=master)

This package provides morphology-related functionality to the [Images.jl](https://github.com/JuliaImages/Images.jl) project.

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

Please check the top-level documentation at [Images.jl](https://github.com/JuliaImages/Images.jl).
