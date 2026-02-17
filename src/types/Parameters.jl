# Parameter types for models and simulations

const DRUG_EFFECTS = (:d, :b, :c)

function normalize_drug_effect(drug_effect)
    de = drug_effect isa AbstractString ? Symbol(drug_effect) : drug_effect
    de in DRUG_EFFECTS || error("drug_effect must be :d, :b, or :c")
    return de
end

struct ModelParams
    b::Float64
    d::Float64
    rho::Float64
    mu::Float64
    sig::Float64
    del::Float64
    al::Float64
    Dc::Float64
    k::Float64
    psi::Float64
    drug_effect::Symbol
end

function ModelParams(; b, d, rho=0.0, mu, sig, del, al, Dc, k, psi, drug_effect="d")
    de = normalize_drug_effect(drug_effect)
    return ModelParams(
        Float64(b), Float64(d), Float64(rho), Float64(mu), Float64(sig), Float64(del),
        Float64(al), Float64(Dc), Float64(k), Float64(psi), de
    )
end

function validate_model_params(params::ModelParams)
    0.0 <= params.rho <= 1.0 || error("rho must be between 0 and 1.")
    0.0 <= params.mu <= 1.0 || error("mu must be between 0 and 1.")
    0.0 <= params.sig <= 1.0 || error("sig must be between 0 and 1.")
    0.0 <= params.del <= 1.0 || error("del must be between 0 and 1.")
    0.0 <= params.al <= 1.0 || error("al must be between 0 and 1.")
    0.0 <= (params.al + params.sig) <= 1.0 || error("al and sig must not sum to > 1.0")
    params.drug_effect in DRUG_EFFECTS || error("drug_effect must be :d, :b, or :c")
    if params.drug_effect === :b
        params.Dc <= params.b || error("When drug_effect == :b, Dc must be <= b.")
    end
    return params
end

struct SimParams
    n0::Int64
    t0::Float64
    tmax::Float64
    t_Pass::Union{Float64, Vector{Float64}}
    Nmax::Int64
    Cc::Int64
    treat_ons::Vector{Float64}
    treat_offs::Vector{Float64}
    Nswitch::Int64
    save_at::Float64
    treat::Bool
    n_Pass::Int64
    epsi::Float64
end

function SimParams(; n0, tmax, Nmax, Cc, Nswitch, treat_ons, treat_offs,
                    t0=0.0, t_Pass=-1.0, save_at=0.5, treat=false,
                    n_Pass=1, epsi=100.0)
    return SimParams(
        Int64(n0),
        Float64(t0),
        Float64(tmax),
        t_Pass,
        Int64(Nmax),
        Int64(Cc),
        Vector{Float64}(treat_ons),
        Vector{Float64}(treat_offs),
        Int64(Nswitch),
        Float64(save_at),
        Bool(treat),
        Int64(n_Pass),
        Float64(epsi)
    )
end

struct ExperimentParams
    n0::Int64
    t_exp::Union{Float64, Vector{Float64}}
    tmax::Float64
    t_Pass::Union{Float64, Vector{Float64}}
    Nseed::Union{Int64, Vector{Int64}}
    Nmax::Int64
    Cc::Int64
    treat_ons::Vector{Float64}
    treat_offs::Vector{Float64}
    t_keep::Vector{Float64}
    Nswitch::Int64
    save_at::Float64
    n_Pass::Int64
    n_rep::Int64
    drug_treatment::Bool
    full_sol::Bool
    run_IC::Bool
    IC_n0::Int64
    IC_tmax::Float64
    IC_treat_on::Float64
    run_colony::Bool
    nCol::Int64
    tCol::Float64
    ColNmax::Int64
end

function ExperimentParams(; n0, t_exp, tmax, t_Pass, Nseed, Nmax, Cc,
                          treat_ons, treat_offs, t_keep, Nswitch,
                          save_at=0.5, n_Pass=1, n_rep=4, drug_treatment=true,
                          full_sol=false, run_IC=false, IC_n0=1000,
                          IC_tmax=4.0, IC_treat_on=1.0, run_colony=false,
                          nCol=1000, tCol=12.0, ColNmax=50)
    return ExperimentParams(
        Int64(n0),
        t_exp,
        Float64(tmax),
        t_Pass,
        Nseed,
        Int64(Nmax),
        Int64(Cc),
        Vector{Float64}(treat_ons),
        Vector{Float64}(treat_offs),
        Vector{Float64}(t_keep),
        Int64(Nswitch),
        Float64(save_at),
        Int64(n_Pass),
        Int64(n_rep),
        Bool(drug_treatment),
        Bool(full_sol),
        Bool(run_IC),
        Int64(IC_n0),
        Float64(IC_tmax),
        Float64(IC_treat_on),
        Bool(run_colony),
        Int64(nCol),
        Float64(tCol),
        Int64(ColNmax)
    )
end