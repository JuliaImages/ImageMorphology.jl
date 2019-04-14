# checks if element is essential for connectivity
function compute_essential_indices()
	essential_index_array = Array{Bool, 1}()
	for index = 0:511
		append!(essential_index_array, maximum(find_components(index)) != maximum(find_components(index & ~2^4)))
	end
	essential_index_array
end

function find_components(index::Int64)
	label_components(reshape(digits(index, base=2, pad=9), (3, 3)), trues(3, 3))
end

function corner_table_lookup(img::AbstractArray{Bool, 2}, corner_table::AbstractArray{Int, 1})
	corner_score = zeros(Int, size(img) .+ 2)
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

function inner_skeleton_loop(img::AbstractArray{Bool, 2}, order::Array{Int, 1}, table::AbstractArray{Bool, 1})
    index = findall(img)[order]
    score_table = [256 128 64; 32 16 8; 4 2 1]
    result = zeros(Int, size(img))
    padded_image = zeros(Int, size(img) .+ 2)
	padded_image[2: size(img)[1] + 1, 2: size(img)[2] + 1] = img
    for i in index
        result[i] = table[sum(padded_image[i[1]:i[1] + 2, i[2]:i[2] + 2] .* score_table) + 1]
        padded_image[i[1] + 1, i[2] + 1] = result[i]
    end
    result
end 

function skeletonize(img::AbstractArray{Bool}, algo::SkeletonizationAlgo)
    skeletonize_impl(img, algo)
end

function skeletonize_impl(img::AbstractArray{Bool, 2}, algo::MedialAxis)

	# check if array has only 0s
	if sum(img) == 0
		return img
	end

	# check if array has only 1s
	if sum(img) == length(img)
		img[1, 1:end] = false
		img[1:end - 1, end] = false
		img[end, 1:end - 1] = false 
		return convert(Array{Bool}, .~img)
	end

	# build look-up table
	patterns = [0:1:511;]
	num_pixels_in_patterns_kernel = sum.(digits.(patterns, base=2))
    center_is_foreground = convert(Array{Bool, 1}, patterns .& 2^4 .> 0)
    table = (center_is_foreground .& compute_essential_indices() .| (num_pixels_in_patterns_kernel .< 3))
    table = convert(Array{Bool, 1}, table)
    
    # store distance transform
    distance = distance_transform(feature_transform(.~img))
    
    # corners handling
    corner_table = 9 .- num_pixels_in_patterns_kernel
    corner_score = corner_table_lookup(img, corner_table)

    # generate order for traversing the shape
    dist_corner_pair = CartesianIndex.(round.(Int, distance[img] .* 10), corner_score[img])
    first_sort = sortperm(dist_corner_pair, by=x->x[2])
    second_sort = sortperm(dist_corner_pair[first_sort], by=x->x[1])
    order = first_sort[second_sort]
  
    convert(Array{Bool, 2}, inner_skeleton_loop(img, order, table))
end
