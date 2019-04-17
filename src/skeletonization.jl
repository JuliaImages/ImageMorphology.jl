# checks if pixel is an articulation/critical node
function compute_critical_indices(num_configurations::Integer)
	critical_index_array = Array{Bool, 1}()
	for index = 0:num_configurations - 1
		append!(critical_index_array, maximum(find_components(index)) != maximum(find_components(index & ~2^4)))
	end
	critical_index_array
end

function find_components(index::Integer)
	label_components(reshape(digits(index, base=2, pad=9), (3, 3)), trues(3, 3))
end

# computes degree of cornerness of pixels
function corner_table_lookup(img::AbstractArray{Bool, 2}, corner_table::AbstractArray{T, 1}) where T <: Integer
	corner_score = zeros(Integer, size(img) .+ 2)
	score_table = [256 128 64; 32 16 8; 4 2 1]

	for i = 2:size(img)[1] + 1
		for j = 2:size(img)[2] + 1
			if img[i - 1, j - 1] 
				corner_score[i - 1:i + 1, j - 1:j + 1] = corner_score[i - 1:i + 1, j - 1:j + 1] + score_table
			end
		end
	end
	corner_table[corner_score[2: size(img)[1] + 1, 2: size(img)[2] + 1] .+ 1]
end

# traverse foreground and reduce image to skeleton
function inner_skeleton_loop(img::AbstractArray{Bool, 2}, order::AbstractArray{T, 1}, table::AbstractArray{Bool, 1}) where T <: Integer
	index = findall(img)[order]
	score_table = [256 128 64; 32 16 8; 4 2 1]
	result = falses(size(img))
	padded_image = padarray(img, Fill(0, (1, 1), (1, 1)))
	for i in index
		result[i] = table[sum(padded_image[i[1] - 1:i[1] + 1, i[2] - 1:i[2] + 1] .* score_table) + 1]
		padded_image[i[1], i[2]] = result[i]
	end
	result
end

"""
```
skeletonize(img, MedialAxisTransform())
```

Returns Medial axis of binary image that matches the dimensions of the input binary image.

# Details
Skeletonization is a process for reducing foreground regions in a binary image to a
skeletal remnant that largely preserves the extent and connectivity of the original
region while throwing away most of the original foreground pixels.

Steps of the algorithm :
- A look up table is built to check whether a pixel is to be removed or not. The pixel
is removed if it has more than one neighbour and if removing it doesn't change the connectedness.
- Distance transform is computed as well as degree of cornerness.
- Foreground pixels are traversed as per distance transform and cornerness.

# Arguments
The img parameter needs to be a two-dimensional array of bool type where true represents
foreground and false represents background.

# Example

Compute the skeleton of rectangular object.
```julia

using Images, ImageMorphology

img = Bool.([0 0 0 0 0
			 0 1 1 1 0
			 0 1 1 1 0
			 0 1 1 1 0
			 0 1 1 1 0
			 0 0 0 0 0])

skeleton = skeletonize(img, MedialAxisTransform())

skeleton = [0 0 0 0 0
			0 1 0 1 0
			0 0 1 0 0
			0 0 1 0 0
			0 1 0 1 0
			0 0 0 0 0]
```

# References
[1] http://homepages.inf.ed.ac.uk/rbf/HIPR2/skeleton.htm
[2] H. Blum in 1967, and in 1968 by L. Calabi and W. E. Hartnett.
"""
function skeletonize(img::AbstractArray{Bool, 2}, algo::MedialAxisTransform)

	# check if array has only 0s
	if all(.~img)
		return img
	end

	# check if array has only 1s
	if sum(img) == length(img)
		img[1, 1:end] .= false
		img[1:end - 1, end] .= false
		img[end, 1:end - 1] .= false
		return .~img
	end

	# possible number of configurations of 3x3 kernel
	num_configurations = 512

	# build look-up table
	patterns = [0:1:num_configurations - 1;]
	num_pixels_in_patterns_kernel = sum.(digits.(patterns, base=2))
	center_is_foreground = patterns .& 2^4 .> 0
	table = (center_is_foreground .& compute_critical_indices(num_configurations) .| (num_pixels_in_patterns_kernel .< 3))

	# store distance transform
	distance = distance_transform(feature_transform(.~img))

	# corners handling
	corner_table = 9 .- num_pixels_in_patterns_kernel
	corner_score = corner_table_lookup(img, corner_table)

	# generate order for traversing the shape
	dist_corner_pair = CartesianIndex.(round.(Integer, distance[img] .* 10), corner_score[img])
	first_sort = sortperm(dist_corner_pair, by=x->x[2])
	second_sort = sortperm(dist_corner_pair[first_sort], by=x->x[1])
	traversal_order = first_sort[second_sort]

	inner_skeleton_loop(img, traversal_order, table)
end
