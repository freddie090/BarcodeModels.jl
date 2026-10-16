# Parameter types for models and simulations

const RESPOP_DRUG_EFFECTS = (:d, :b, :c)
const RESDMG_DRUG_EFFECTS = (:d, :b, :c)
const TRANSITION_COUPLINGS = (:birth_coupled, :independent_rate)

function normalize_transition_coupling(transition_coupling)
    coupling = transition_coupling isa AbstractString ? Symbol(transition_coupling) : transition_coupling
    coupling in TRANSITION_COUPLINGS || error("transition_coupling must be :birth_coupled or :independent_rate.")
    return coupling
end

function normalize_respop_drug_effect(drug_effect)
    de = drug_effect isa AbstractString ? Symbol(drug_effect) : drug_effect
    de in RESPOP_DRUG_EFFECTS || error("ResPop drug_effect must be :d, :b, or :c")
    return de
end

function normalize_resdmg_drug_effect(drug_effect)
    de = drug_effect isa AbstractString ? Symbol(drug_effect) : drug_effect
    de in RESDMG_DRUG_EFFECTS || error("ResDmg drug_effect must be :d, :b, or :c")
    return de
end

struct ResPopParams
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
    transition_coupling::Symbol
end

function ResPopParams(; b, d, rho=0.0, mu, sig, del, al, Dc, k, psi, drug_effect="d",
                      transition_coupling=:birth_coupled)
    de = normalize_respop_drug_effect(drug_effect)
    return ResPopParams(
        Float64(b), Float64(d), Float64(rho), Float64(mu), Float64(sig), Float64(del),
        Float64(al), Float64(Dc), Float64(k), Float64(psi), de,
        normalize_transition_coupling(transition_coupling)
    )
end

ResPopParams(b, d, rho, mu, sig, del, al, Dc, k, psi, drug_effect) =
    ResPopParams(; b, d, rho, mu, sig, del, al, Dc, k, psi, drug_effect)

struct ResDmgParams
    b::Float64
    d::Float64
    rho::Float64
    mu::Float64
    sig::Float64
    del::Float64
    ome::Float64
    zet_S::Float64
    zet_R::Float64
    Dc::Float64
    k::Float64
    psi::Float64
    drug_effect::Symbol
    transition_coupling::Symbol
end

function ResDmgParams(; b, d, rho=0.0, mu, sig, del, ome, zet_S, zet_R, Dc, k, psi, drug_effect="d",
                      transition_coupling=:birth_coupled)
    de = normalize_resdmg_drug_effect(drug_effect)
    return ResDmgParams(
        Float64(b), Float64(d), Float64(rho), Float64(mu), Float64(sig), Float64(del),
        Float64(ome), Float64(zet_S), Float64(zet_R), Float64(Dc), Float64(k), Float64(psi), de,
        normalize_transition_coupling(transition_coupling)
    )
end

ResDmgParams(b, d, rho, mu, sig, del, ome, zet_S, zet_R, Dc, k, psi, drug_effect) =
    ResDmgParams(; b, d, rho, mu, sig, del, ome, zet_S, zet_R, Dc, k, psi, drug_effect)

function validate_model_params(params::ResPopParams)
    for name in (:b, :d, :Dc, :k)
        rate = getproperty(params, name)
        isfinite(rate) && rate >= 0.0 || error("$(name) must be finite and >= 0.")
    end
    params.b > params.d || error("Baseline population growth rate will lead to contraction or no growth; require b > d.")
    0.0 <= params.rho <= 1.0 || error("rho must be between 0 and 1.")
    normalize_transition_coupling(params.transition_coupling)
    if params.transition_coupling === :birth_coupled
        0.0 <= params.mu <= 1.0 || error("mu must be between 0 and 1.")
        0.0 <= params.sig <= 1.0 || error("sig must be between 0 and 1.")
        0.0 <= params.al <= 1.0 || error("al must be between 0 and 1.")
        0.0 <= (params.al + params.sig) <= 1.0 || error("al and sig must not sum to > 1.0")
    else
        for name in (:mu, :sig, :al)
            rate = getproperty(params, name)
            isfinite(rate) && rate >= 0.0 || error("$(name) must be finite and >= 0 in independent_rate mode.")
        end
    end
    0.0 <= params.del <= 1.0 || error("del must be between 0 and 1.")
    isfinite(params.psi) && params.psi <= 1.0 || error("psi must be finite and <= 1.0.")
    params.drug_effect in RESPOP_DRUG_EFFECTS || error("ResPop drug_effect must be :d, :b, or :c")
    if params.drug_effect === :b
        params.Dc <= params.b || error("When drug_effect == :b, Dc must be <= b.")

        bR = params.b * (1 - params.del)
        psi_scale = 1 - params.psi
        if psi_scale > 0.0
            params.Dc <= (bR / psi_scale) || error("When drug_effect == :b, Dc*(1-psi) must be <= b*(1-del) to keep resistant birth non-negative. Use drug_effect == :c if stronger drug effects are intended.")
        end
    end
    return params
end

function validate_model_params_strict(params::ResPopParams)
    validate_model_params(params)
    params.psi >= 0.0 || error("psi must be between 0 and 1.")
    return params
end

function validate_model_params(params::ResDmgParams)
    for name in (:b, :d, :Dc, :k, :ome, :zet_S, :zet_R)
        rate = getproperty(params, name)
        isfinite(rate) && rate >= 0.0 || error("$(name) must be finite and >= 0.")
    end
    params.b > params.d || error("Baseline population growth rate will lead to contraction or no growth; require b > d.")
    0.0 <= params.rho <= 1.0 || error("rho must be between 0 and 1.")
    normalize_transition_coupling(params.transition_coupling)
    if params.transition_coupling === :birth_coupled
        0.0 <= params.mu <= 1.0 || error("mu must be between 0 and 1.")
        0.0 <= params.sig <= 1.0 || error("sig must be between 0 and 1.")
    else
        for name in (:mu, :sig)
            rate = getproperty(params, name)
            isfinite(rate) && rate >= 0.0 || error("$(name) must be finite and >= 0 in independent_rate mode.")
        end
    end
    0.0 <= params.del <= 1.0 || error("del must be between 0 and 1.")
    isfinite(params.psi) && params.psi <= 1.0 || error("psi must be finite and <= 1.0.")
    params.drug_effect in RESDMG_DRUG_EFFECTS || error("ResDmg drug_effect must be :d, :b, or :c")
    if params.drug_effect === :b
        params.Dc <= params.b || error("When drug_effect == :b, Dc must be <= b.")

        bR = params.b * (1 - params.del)
        psi_scale = 1 - params.psi
        if psi_scale > 0.0
            params.Dc <= (bR / psi_scale) || error("When drug_effect == :b, Dc*(1-psi) must be <= b*(1-del) to keep resistant birth non-negative. Use drug_effect == :c if stronger drug effects are intended.")
        end
    end
    return params
end

function validate_model_params_strict(params::ResDmgParams)
    validate_model_params(params)
    params.psi >= 0.0 || error("psi must be between 0 and 1.")
    if params.zet_R > 0.0
        params.zet_R >= params.zet_S || error("zet_R must be >= zet_S for resistant cells to have at least as high repair propensity as sensitive cells.")
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
    Aswitch::Float64
end

"""Validate the source-population and full event-rate hybrid thresholds."""
function validate_hybrid_switching_thresholds(Nswitch, Aswitch)
    Nswitch isa Integer || error("Nswitch must be an integer.")
    Nswitch > 0 || error("Nswitch must be > 0.")
    Aswitch isa Real || error("Aswitch must be a number.")
    activity_threshold = Float64(Aswitch)
    isfinite(activity_threshold) && activity_threshold > 0 || error("Aswitch must be finite and > 0.")
    return nothing
end

function SimParams(; n0, tmax, Nmax, Cc, Nswitch, treat_ons, treat_offs,
                    t0=0.0, t_Pass=Float64[], save_at=0.5, treat=false,
                    Aswitch=1000.0)
    validate_hybrid_switching_thresholds(Nswitch, Aswitch)
    Int64(Cc) > 0 || error("Cc must be > 0.")

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
        Float64(Aswitch)
    )
end

function validate_abm_params(; skew_lib::Bool, use_lib_probs::Bool, Nbc::Int64, bc_probs::Vector{Float64})
    skew_lib && use_lib_probs && error("use_lib_probs and skew_lib cannot both be true.")
    if use_lib_probs
        Nbc > 0 || error("Nbc must be > 0 when use_lib_probs is true.")
        length(bc_probs) == Nbc || error("bc_probs must have length Nbc when use_lib_probs is true.")
        all(x -> x >= 0.0, bc_probs) || error("bc_probs must contain nonnegative values.")
        isapprox(sum(bc_probs), 1.0; atol = 1e-8) || error("bc_probs must sum to 1.0.")
    end
    return nothing
end

struct ABMParams
    Nbuff::Int64
    t_frac::Float64
    dt_save_at::Float64
    skew_lib::Bool
    use_lib_probs::Bool
    bc_unif::Float64
    Nbc::Int64
    bc_probs::Vector{Float64}
    sub_sample_cells::Bool
    K::Int64
end

function ABMParams(; Nbuff=100000, t_frac=0.005, dt_save_at=0.1,
                   skew_lib=false, use_lib_probs=false, bc_unif=0.0, Nbc=0,
                   bc_probs=Float64[],
                   sub_sample_cells=false, K=0)
    bc_probs_vec = Vector{Float64}(bc_probs)
    validate_abm_params(skew_lib = Bool(skew_lib), use_lib_probs = Bool(use_lib_probs),
                        Nbc = Int64(Nbc), bc_probs = bc_probs_vec)

    return ABMParams(
        Int64(Nbuff),
        Float64(t_frac),
        Float64(dt_save_at),
        Bool(skew_lib),
        Bool(use_lib_probs),
        Float64(bc_unif),
        Int64(Nbc),
        bc_probs_vec,
        Bool(sub_sample_cells),
        Int64(K)
    )
end

"""Shared runtime settings for standard and lineage-aware ABM simulations."""
struct ABMSimParams
    t0::Float64
    tmax::Float64
    Nmax::Int64
    Cc::Int64
    treat_ons::Vector{Float64}
    treat_offs::Vector{Float64}
    dt_save_at::Float64
    R_real::String
    t_frac::Float64
    Passage::Int64
    drug_effect::Symbol
end

function ABMSimParams(; t0 = 0.0, tmax, Nmax, Cc, treat_ons, treat_offs,
    dt_save_at = 0.1, R_real = "b", t_frac = 0.05, Passage = 1, drug_effect = :d)

    return ABMSimParams(
        Float64(t0),
        Float64(tmax),
        Int64(Nmax),
        Int64(Cc),
        Vector{Float64}(treat_ons),
        Vector{Float64}(treat_offs),
        Float64(dt_save_at),
        String(R_real),
        Float64(t_frac),
        Int64(Passage),
        normalize_respop_drug_effect(drug_effect)
    )
end

struct ExperimentParams
    n0::Int64
    t_exp::Union{Float64, Vector{Float64}}
    tmax::Union{Float64, Vector{Float64}}
    t_Pass::Union{Float64, Vector{Float64}}
    Nseed::Union{Int64, Vector{Int64}}
    Nmax::Int64
    Cc::Int64
    treat_ons::Vector{Float64}
    treat_offs::Vector{Float64}
    t_keep::Vector{Float64}
    Nswitch::Int64
    Aswitch::Float64
    save_at::Float64
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

struct SimpleSimParams
    n0::Int64
    tmax::Float64
    Nmax::Int64
    Cc::Int64
    treat_ons::Vector{Float64}
    treat_offs::Vector{Float64}
    Nswitch::Int64
    Aswitch::Float64
    save_at::Float64
    drug_treatment::Bool
end

function validate_simple_params(; n0, tmax, Nmax, Cc,
                                treat_ons, treat_offs,
                                Nswitch, Aswitch,
                                save_at, drug_treatment)
    n0 isa Integer || error("n0 must be an integer.")
    n0 > 0 || error("n0 must be > 0.")

    tmax isa Real || error("tmax must be a number.")
    tmax > 0 || error("tmax must be > 0.")

    Nmax isa Integer || error("Nmax must be an integer.")
    Nmax > 0 || error("Nmax must be > 0.")

    Cc isa Integer || error("Cc must be an integer.")
    Cc > 0 || error("Cc must be > 0.")
    Nmax <= Cc || error("Nmax must be <= Cc.")

    treat_ons isa AbstractVector{<:Real} || error("treat_ons must be a vector of numbers.")
    treat_offs isa AbstractVector{<:Real} || error("treat_offs must be a vector of numbers.")
    length(treat_ons) == length(treat_offs) || error("treat_ons and treat_offs must have matching lengths.")

    validate_hybrid_switching_thresholds(Nswitch, Aswitch)

    save_at isa Real || error("save_at must be a number.")
    save_at > 0 || error("save_at must be > 0.")

    drug_treatment isa Bool || error("drug_treatment must be Bool.")

    return nothing
end

function SimpleSimParams(; n0, tmax, Nmax, Cc,
                         treat_ons=Float64[], treat_offs=Float64[],
                         Nswitch=100, Aswitch=1000.0,
                         save_at=0.5, drug_treatment=true)
    validate_simple_params(
        n0 = n0,
        tmax = tmax,
        Nmax = Nmax,
        Cc = Cc,
        treat_ons = treat_ons,
        treat_offs = treat_offs,
        Nswitch = Nswitch,
        Aswitch = Aswitch,
        save_at = save_at,
        drug_treatment = drug_treatment
    )

    return SimpleSimParams(
        Int64(n0),
        Float64(tmax),
        Int64(Nmax),
        Int64(Cc),
        Vector{Float64}(treat_ons),
        Vector{Float64}(treat_offs),
        Int64(Nswitch),
        Float64(Aswitch),
        Float64(save_at),
        Bool(drug_treatment)
    )
end

function validate_experiment_params(; n0, t_exp, tmax, t_Pass, Nseed, Nmax, Cc,
                                    treat_ons, treat_offs, t_keep, Nswitch,
                                    Aswitch, save_at, n_rep, drug_treatment,
                                    full_sol, run_IC, IC_n0, IC_tmax,
                                    IC_treat_on, run_colony, nCol, tCol, ColNmax)
    n0 isa Integer || error("n0 must be an integer.")
    n0 > 0 || error("n0 must be > 0.")

    (t_exp isa Real || t_exp isa AbstractVector{<:Real}) || error("t_exp must be a number or a vector of numbers.")

    (tmax isa Real || tmax isa AbstractVector{<:Real}) || error("tmax must be a number or a vector of numbers.")
    if tmax isa AbstractVector
        all(x -> x > 0.0, tmax) || error("All tmax values must be > 0.0.")
    else
        tmax > 0 || error("tmax must be > 0.")
    end

    (t_Pass isa Real || t_Pass isa AbstractVector{<:Real}) || error("t_Pass must be a number or a vector of numbers.")
    if t_Pass isa AbstractVector
        all(x -> x > 0.0, t_Pass) || error("All t_Pass values must be > 0.0. Use Float64[] for no passage events.")
    else
        t_Pass > 0.0 || error("t_Pass must be > 0.0. Use Float64[] for no passage events.")
    end

    (Nseed isa Integer || Nseed isa AbstractVector{<:Integer}) || error("Nseed must be an integer or a vector of integers.")
    if Nseed isa Integer
        Nseed > 0 || error("Nseed must be > 0.")
    else
        all(x -> x > 0, Nseed) || error("All Nseed values must be > 0.")
    end

    Nmax isa Integer || error("Nmax must be an integer.")
    Nmax > 0 || error("Nmax must be > 0.")

    Cc isa Integer || error("Cc must be an integer.")
    Cc > 0 || error("Cc must be > 0.")
    Nmax <= Cc || error("Nmax must be <= Cc.")

    treat_ons isa AbstractVector{<:Real} || error("treat_ons must be a vector of numbers.")
    treat_offs isa AbstractVector{<:Real} || error("treat_offs must be a vector of numbers.")
    t_keep isa AbstractVector{<:Real} || error("t_keep must be a vector of numbers.")

    validate_hybrid_switching_thresholds(Nswitch, Aswitch)

    save_at isa Real || error("save_at must be a number.")
    save_at > 0 || error("save_at must be > 0.")

    n_rep isa Integer || error("n_rep must be an integer.")
    n_rep > 0 || error("n_rep must be > 0.")

    drug_treatment isa Bool || error("drug_treatment must be Bool.")
    full_sol isa Bool || error("full_sol must be Bool.")
    run_IC isa Bool || error("run_IC must be Bool.")

    IC_n0 isa Integer || error("IC_n0 must be an integer.")
    IC_n0 > 0 || error("IC_n0 must be > 0.")

    IC_tmax isa Real || error("IC_tmax must be a number.")
    IC_tmax > 0 || error("IC_tmax must be > 0.")

    IC_treat_on isa Real || error("IC_treat_on must be a number.")

    run_colony isa Bool || error("run_colony must be Bool.")

    nCol isa Integer || error("nCol must be an integer.")
    nCol > 0 || error("nCol must be > 0.")

    tCol isa Real || error("tCol must be a number.")
    tCol > 0 || error("tCol must be > 0.")

    ColNmax isa Integer || error("ColNmax must be an integer.")
    ColNmax > 0 || error("ColNmax must be > 0.")

    return nothing
end

function ExperimentParams(; n0, t_exp, tmax, t_Pass, Nseed, Nmax, Cc,
                          treat_ons, treat_offs, t_keep, Nswitch,
                          Aswitch=1000.0,
                          save_at=0.5, n_rep=4, drug_treatment=true,
                          full_sol=false, run_IC=false, IC_n0=1000,
                          IC_tmax=4.0, IC_treat_on=1.0, run_colony=false,
                          nCol=1000, tCol=12.0, ColNmax=50)
    validate_experiment_params(
        n0 = n0,
        t_exp = t_exp,
        tmax = tmax,
        t_Pass = t_Pass,
        Nseed = Nseed,
        Nmax = Nmax,
        Cc = Cc,
        treat_ons = treat_ons,
        treat_offs = treat_offs,
        t_keep = t_keep,
        Nswitch = Nswitch,
        Aswitch = Aswitch,
        save_at = save_at,
        n_rep = n_rep,
        drug_treatment = drug_treatment,
        full_sol = full_sol,
        run_IC = run_IC,
        IC_n0 = IC_n0,
        IC_tmax = IC_tmax,
        IC_treat_on = IC_treat_on,
        run_colony = run_colony,
        nCol = nCol,
        tCol = tCol,
        ColNmax = ColNmax
    )

    if (tmax isa AbstractVector) && !(t_Pass isa AbstractVector && isempty(t_Pass))
        error("Vector tmax is only supported when t_Pass is an empty vector (no passage events).")
    end

    tmax_out = tmax isa AbstractVector ? Vector{Float64}(tmax) : Float64(tmax)

    return ExperimentParams(
        Int64(n0),
        t_exp,
        tmax_out,
        t_Pass,
        Nseed,
        Int64(Nmax),
        Int64(Cc),
        Vector{Float64}(treat_ons),
        Vector{Float64}(treat_offs),
        Vector{Float64}(t_keep),
        Int64(Nswitch),
        Float64(Aswitch),
        Float64(save_at),
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

