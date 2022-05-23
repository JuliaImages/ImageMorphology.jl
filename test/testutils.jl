# build a random mask se with extra symmetricity constraint
function rand_se_mask(r, N; symmetric=false)
    se = centered(rand(Bool, ntuple(_ -> r, N)))
    symmetric || return se
    dual_se = copy(se)
    for i in CartesianIndices(dual_se)
        dual_se[i] = se[-i]
    end

    out = @. (Int(dual_se) + Int(se)) ÷ 2
    out[OffsetArrays.center(out)...] = 1
    return Bool.(out)
end
