# ImageMorphology

[![Build Status](https://travis-ci.org/juliohm/ImageMorphology.jl.svg?branch=master)](https://travis-ci.org/juliohm/ImageMorphology.jl)
[![codecov.io](http://codecov.io/github/juliohm/ImageMorphology.jl/coverage.svg?branch=master)](http://codecov.io/github/juliohm/ImageMorphology.jl?branch=master)

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
```
