


function jacobian_is_stable(jac::Matrix)
    ThreadsX.all(real.(eigvals(jac)) .< 0)
end
function fixedpoint_is_stable(fp, nullcline_params)
    if !all(abs.(derive_vector_fn(nullcline_params)(fp)) .< sqrt(eps()))
        @error "Mismatch between loaded data and assumed parameters: $nullcline_params"
    end
    jac_fn = derive_jacobian_fn(nullcline_params)
    jac = jac_fn(fp)
    jacobian_is_stable(jac)
end

function jacobian_is_oscillatory(jac::Matrix)
    any(imag.(eigvals(jac)) .!= 0)
end
function fixedpoint_is_oscillatory(nullcline_params::AbstractNullclineParams, fp)
    @assert manhattan_norm(derive_vector_fn(nullcline_params)(fp)) .< sqrt(eps()) derive_vector_fn(nullcline_params)(fp)
    jac_fn = derive_jacobian_fn(nullcline_params)
    jac = jac_fn(fp)
    jacobian_is_oscillatory(jac)
end

const STABILITY_MARKERS = OffsetArray([:circle, :star, :xcross], -1:1)
function fixedpoint_stability(nullcline_params, fp)
    jac_fn = derive_jacobian_fn(nullcline_params)
    jac = jac_fn(fp)
    if all(real.(eigvals(jac)) .< 0)
        return -1
    elseif all(real.(eigvals(jac)) .== 0)
        return 0
    else
        return 1
    end
end

function count_stable_fps(fp_arr::NamedDimsArray, fp_axes, prototype_name::String, fixed_nt::NamedTuple)
    prototype = get_prototype(prototype_name)
    map(enumerate_nda_dims(fp_arr, fp_axes)) do (nt, fps)
        params = get_nullcline_params(prototype(; fixed_nt..., nt...))
        count(fixedpoint_is_stable.(fps, Ref(params)))
    end
end

function filter_stable_fps(fp_arr::NamedDimsArray{NAMES}, fp_axes::NamedTuple{NAMES}, prototype_name::String, fixed_nt::NamedTuple) where NAMES
    fp_arr = copy(fp_arr)
    filter_stable_fps!(fp_arr, fp_axes, prototype_name, fixed_nt)
    return fp_arr
end
function filter_stable_fps!(fp_arr::NamedDimsArray{NAMES}, fp_axes::NamedTuple{NAMES}, prototype_name::String, fixed_nt::NamedTuple) where NAMES
    prototype = get_prototype(prototype_name)
    for (nt, fps) ∈ enumerate_nda_dims(fp_arr, fp_axes)
        params = get_nullcline_params(prototype(; fixed_nt..., nt...))
        filter_stable_fps!(fps, params)
    end
end
function filter_stable_fps(fps::AbstractVector, params::AbstractNullclineParams)
    filter(fp -> fixedpoint_is_stable(params, fp), fps)
end
function filter_stable_fps!(fps::AbstractVector, params::AbstractNullclineParams)
    keepat!(fps, fixedpoint_is_stable.(fp, Ref(params)))
end
