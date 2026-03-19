function _kw(kwargs, name, default)
    return haskey(kwargs, name) ? kwargs[name] : default
end

_nseed_last(x::Int64) = x
_nseed_last(x::Vector{Int64}) = last(x)

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

function _copy_resdmg_params(params::ResDmgParams; rho = params.rho, al = params.al, drug_effect = params.drug_effect)
    return ResDmgParams(
        b = params.b,
        d = params.d,
        rho = rho,
        mu = params.mu,
        sig = params.sig,
        del = params.del,
        al = al,
        ome = params.ome,
        zet = params.zet,
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

