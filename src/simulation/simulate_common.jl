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

"""Normalize a passage schedule into sorted positive passage times."""
function _passage_times(t_Pass::Union{Float64, Vector{Float64}}, tmax::Float64)
    if t_Pass isa AbstractVector
        times = sort(unique(Float64.(t_Pass)))
        any(t -> t <= 0.0, times) && error("All t_Pass values must be > 0.0. Use Float64[] for no passage events.")
        return times
    else
        t_val = Float64(t_Pass)
        t_val > 0.0 || error("t_Pass must be > 0.0. Use Float64[] for no passage events.")
        return [t_val]
    end
end

function _experiment_condition_design(drug_treatment::Bool, inc_control::Bool, n_rep::Int64)
    design = Vector{NamedTuple{(:cond, :treat, :rep), Tuple{String, Bool, Int64}}}()
    if inc_control || !drug_treatment
        for rep in 1:n_rep
            push!(design, (cond = "CO", treat = false, rep = Int64(rep)))
        end
    end
    if drug_treatment
        for rep in 1:n_rep
            push!(design, (cond = "DT", treat = true, rep = Int64(rep)))
        end
    end
    return design
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

function _copy_respop_invivo_params(params::ResPopInVivoParams;
    rho = params.rho, al = params.al, drug_effect = params.drug_effect,
    fEG1 = params.fEG1, pEG = params.pEG, sEG = params.sEG)

    return ResPopInVivoParams(
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
        drug_effect = drug_effect,
        fEG1 = fEG1,
        pEG = pEG,
        sEG = sEG
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

function _with_drug_effect(model::ResPop_ABM_EvBC, de::Symbol)
    de == model.params.drug_effect && return model
    return ResPop_ABM_EvBC(_copy_respop_params(model.params; drug_effect = de); abm = model.abm)
end

function _with_drug_effect(model::ResDmg_ABM_EvBC, de::Symbol)
    de == model.params.drug_effect && return model
    return ResDmg_ABM_EvBC(_copy_resdmg_params(model.params; drug_effect = de); abm = model.abm)
end

function _with_drug_effect(model::ResDmg, de::Symbol)
    de == model.params.drug_effect && return model
    params_eff = _copy_resdmg_params(model.params; drug_effect = de)
    return ResDmg(params_eff)
end

function _with_drug_effect(model::ResPopInVivo, de::Symbol)
    de == model.params.drug_effect && return model
    return ResPopInVivo(_copy_respop_invivo_params(model.params; drug_effect = de))
end

function _with_drug_effect(model::ResPopInVivo_ABM, de::Symbol)
    de == model.params.drug_effect && return model
    return ResPopInVivo_ABM(_copy_respop_invivo_params(model.params; drug_effect = de); abm = model.abm)
end

