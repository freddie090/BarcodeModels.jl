function _kw(kwargs, name, default)
    return haskey(kwargs, name) ? kwargs[name] : default
end

_nseed_last(x::Int64) = x
_nseed_last(x::Vector{Int64}) = last(x)

function _validate_tmax_vector_constraints(tmax, t_Pass)
    if (tmax isa AbstractVector) && !(t_Pass isa AbstractVector && isempty(t_Pass))
        error("Vector tmax is only supported when t_Pass is an empty vector (no passage events).")
    end
    return nothing
end

function _validate_tmax_length(tmax, n_rep::Int64)
    if tmax isa AbstractVector
        length(tmax) == n_rep || error("Length of tmax vector must match n_rep.")
    end
    return nothing
end

function _replicate_tmax(tmax, n_rep::Int64, i::Int64)
    _validate_tmax_length(tmax, n_rep)
    return tmax isa AbstractVector ? Float64(tmax[i]) : Float64(tmax)
end

function _copy_respop_params(params::ResPopParams; rho = params.rho, al = params.al, drug_effect = params.drug_effect)
    return ResPopParams(
        b = params.b,
        d = params.d,
        rho = rho,
        mu = params.mu,
        sig = params.sig,
        del = params.del,
        al = al,
        Dc = params.Dc,
        k = params.k,
        psi = params.psi,
        drug_effect = drug_effect
    )
end

function _copy_resdmg_params(params::ResDmgParams; rho = params.rho, drug_effect = params.drug_effect)
    return ResDmgParams(
        b = params.b,
        d = params.d,
        rho = rho,
        mu = params.mu,
        sig = params.sig,
        del = params.del,
        ome = params.ome,
        zet_S = params.zet_S,
        zet_R = params.zet_R,
        Dc = params.Dc,
        k = params.k,
        psi = params.psi,
        drug_effect = drug_effect
    )
end

function _with_drug_effect(model::ResPop, de::Symbol)
    de == model.params.drug_effect && return model
    return ResPop(_copy_respop_params(model.params; drug_effect = de))
end

function _with_drug_effect(model::ResPop_ABM, de::Symbol)
    de == model.params.drug_effect && return model
    return ResPop_ABM(_copy_respop_params(model.params; drug_effect = de); abm = model.abm)
end

function _with_drug_effect(model::ResDmg_ABM, de::Symbol)
    de == model.params.drug_effect && return model
    return ResDmg_ABM(_copy_resdmg_params(model.params; drug_effect = de); abm = model.abm)
end

function _with_drug_effect(model::ResDmg, de::Symbol)
    de == model.params.drug_effect && return model
    params_eff = _copy_resdmg_params(model.params; drug_effect = de)
    return ResDmg(params_eff)
end

