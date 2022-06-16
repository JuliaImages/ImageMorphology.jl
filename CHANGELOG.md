# ImageMorphology

## Version `v0.4.0`

This release introduces a heavy redesign of the codebase to appropriately support structuring
element concept. A lot of APIs are revisited, and documentation are added for those revisited APIs.
Please checkout [docs](https://juliaimages.org/ImageMorphology.jl/stable/) for more information of
those APIs.

The following briefly summarizes the noteworthy changes:

- ![BREAKING][badge-breaking] The return value of `component_boxes`, `component_indices`,
  `component_lengths` and `component_centeroids` is now 0-based OffsetArray ([#96][github-96])
- ![BREAKING][badge-breaking] The return value of `component_boxes` becomes a vector of
  `CartesianIndices` ([#96][github-96])
- ![Deprecation][badge-deprecation] `component_subscripts` is deprecated in favor of
  `component_indices` ([#96][github-96])
- ![Deprecation][badge-deprecation] `extrem_filt!` is deprecated in favor of `extreme_filter`
- ![Deprecation][badge-deprecation] `morphogradient`/`morpholaplace` is deprecated in favor of
  `mgradient`/`mlaplacian` ([#88][github-88])
- ![Feature][badge-feature] A lot of helper functions for structuring element are added: `strel`,
  `strel_chain`, `strel_product`, `strel_type`, `strel_size`, `strel_diamond`, `strel_box`
- ![Feature][badge-feature] Basic morphological operations `erode`, `dilate`, `bothat`, `tophat`,
  `opening`, `closing`, `mgradient`, `mlaplacian` now support generic structuring element inputs.
  The associated in-place versions such as `erode!` are provided as well.
- ![Enhancement][badge-enhancement] The performance for diamond-shape structuring element is
  improved by 10x -- it simply beats MATLAB by 3x. ([#97][github-97])

[github-88]: https://github.com/JuliaImages/ImageMorphology.jl/pull/88
[github-96]: https://github.com/JuliaImages/ImageMorphology.jl/pull/96
[github-97]: https://github.com/JuliaImages/ImageMorphology.jl/pull/97


[badge-breaking]: https://img.shields.io/badge/BREAKING-red.svg
[badge-deprecation]: https://img.shields.io/badge/deprecation-orange.svg
[badge-feature]: https://img.shields.io/badge/feature-green.svg
[badge-enhancement]: https://img.shields.io/badge/enhancement-blue.svg
[badge-bugfix]: https://img.shields.io/badge/bugfix-purple.svg
[badge-security]: https://img.shields.io/badge/security-black.svg
[badge-experimental]: https://img.shields.io/badge/experimental-lightgrey.svg
[badge-maintenance]: https://img.shields.io/badge/maintenance-gray.svg

<!--
# Badges

![BREAKING][badge-breaking]
![Deprecation][badge-deprecation]
![Feature][badge-feature]
![Enhancement][badge-enhancement]
![Bugfix][badge-bugfix]
![Security][badge-security]
![Experimental][badge-experimental]
![Maintenance][badge-maintenance]
-->
