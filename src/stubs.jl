
function field_functions(params)
    error("undefined")
end

function phase_space_bounds(params)
    error("undefined")
end

# FIXME these two probably don't need to be stubs
# could autothread?
function derive_vector_fn(params)
    field_fns = field_functions(params)
    function du_dv(uv)
        [fn(uv..., params) for fn in field_fns]
    end
end
function derive_vector_fn!(params)
    field_fns = field_functions(params)
    function du_dv!(dudv, uv)
        for i ∈ 1:length(dudv)
            dudv[i] = field_fns[i](uv..., params)
        end
    end
end

# FIXME these two aren't stubs
function derive_jacobian_fn(
    nullcline_params::AbstractNullclineParams
)
    du_dv =  derive_vector_fn(nullcline_params)
    uv -> ForwardDiff.jacobian(du_dv, uv)
end

function derive_jacobian_fn!(
    nullcline_params::AbstractNullclineParams
)
    du_dv! =  derive_vector_fn!(nullcline_params)
    (result, dudv, uv) -> jacobian!(result, du_dv!, dudv, uv)
end

