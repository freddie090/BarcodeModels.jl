# Parameter types for models and simulations

const RESPOP_DRUG_EFFECTS = (:d, :b, :c)
const RESDMG_DRUG_EFFECTS = (:d, :b, :c)

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
end

function ResPopParams(; b, d, rho=0.0, mu, sig, del, al, Dc, k, psi, drug_effect="d")
    de = normalize_respop_drug_effect(drug_effect)
    return ResPopParams(
        Float64(b), Float64(d), Float64(rho), Float64(mu), Float64(sig), Float64(del),
        Float64(al), Float64(Dc), Float64(k), Float64(psi), de
    )
end

struct ResDmgParams
    b::Float64
    d::Float64
    rho::Float64
    mu::Float64
    sig::Float64
    del::Float64
    al::Float64
    ome::Float64
    zet::Float64
    Dc::Float64
    k::Float64
    psi::Float64
    drug_effect::Symbol
end

function ResDmgParams(; b, d, rho=0.0, mu, sig, del, al, ome, zet, Dc, k, psi, drug_effect="d")
    de = normalize_resdmg_drug_effect(drug_effect)
    return ResDmgParams(
        Float64(b), Float64(d), Float64(rho), Float64(mu), Float64(sig), Float64(del),
        Float64(al), Float64(ome), Float64(zet), Float64(Dc), Float64(k), Float64(psi), de
    )
end

function validate_model_params(params::ResPopParams)
    0.0 <= params.rho <= 1.0 || error("rho must be between 0 and 1.")
    0.0 <= params.mu <= 1.0 || error("mu must be between 0 and 1.")
    0.0 <= params.sig <= 1.0 || error("sig must be between 0 and 1.")
    0.0 <= params.del <= 1.0 || error("del must be between 0 and 1.")
    0.0 <= params.al <= 1.0 || error("al must be between 0 and 1.")
    0.0 <= params.psi <= 1.0 || error("psi must be between 0 and 1.")
    0.0 <= (params.al + params.sig) <= 1.0 || error("al and sig must not sum to > 1.0")
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

function validate_model_params(params::ResDmgParams)
    0.0 <= params.rho <= 1.0 || error("rho must be between 0 and 1.")
    0.0 <= params.mu <= 1.0 || error("mu must be between 0 and 1.")
    0.0 <= params.sig <= 1.0 || error("sig must be between 0 and 1.")
    0.0 <= params.del <= 1.0 || error("del must be between 0 and 1.")
    0.0 <= params.al <= 1.0 || error("al must be between 0 and 1.")
    0.0 <= params.psi <= 1.0 || error("psi must be between 0 and 1.")
    params.ome >= 0.0 || error("ome must be >= 0.")
    params.zet >= 0.0 || error("zet must be >= 0.")
    0.0 <= (params.al + params.sig) <= 1.0 || error("al and sig must not sum to > 1.0")
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
    N_trans_switch::Float64
end

function SimParams(; n0, tmax, Nmax, Cc, Nswitch, treat_ons, treat_offs,
                    t0=0.0, t_Pass=Float64[], save_at=0.5, treat=false,
                    N_trans_switch=1000.0)
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
        Float64(N_trans_switch)
    )
end

struct ABMParams
    Nbuff::Int64
    t_frac::Float64
    dt_save_at::Float64
    skew_lib::Bool
    bc_unif::Float64
    Nbc::Int64
    sub_sample_cells::Bool
    K::Int64
end

function ABMParams(; Nbuff=100000, t_frac=0.005, dt_save_at=0.1,
                   skew_lib=false, bc_unif=0.0, Nbc=0,
                   sub_sample_cells=false, K=0)
    return ABMParams(
        Int64(Nbuff),
        Float64(t_frac),
        Float64(dt_save_at),
        Bool(skew_lib),
        Float64(bc_unif),
        Int64(Nbc),
        Bool(sub_sample_cells),
        Int64(K)
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
    N_trans_switch::Float64
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

function validate_experiment_params(; n0, t_exp, tmax, t_Pass, Nseed, Nmax, Cc,
                                    treat_ons, treat_offs, t_keep, Nswitch,
                                    N_trans_switch, save_at, n_rep, drug_treatment,
                                    full_sol, run_IC, IC_n0, IC_tmax,
                                    IC_treat_on, run_colony, nCol, tCol, ColNmax)
    n0 isa Integer || error("n0 must be an integer.")
    n0 > 0 || error("n0 must be > 0.")

    (t_exp isa Real || t_exp isa AbstractVector{<:Real}) || error("t_exp must be a number or a vector of numbers.")

    tmax isa Real || error("tmax must be a number.")
    tmax > 0 || error("tmax must be > 0.")

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

    Nswitch isa Integer || error("Nswitch must be an integer.")
    Nswitch > 0 || error("Nswitch must be > 0.")

    N_trans_switch isa Real || error("N_trans_switch must be a number.")
    N_trans_switch > 0 || error("N_trans_switch must be > 0.")

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
                          N_trans_switch=1000.0,
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
        N_trans_switch = N_trans_switch,
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
        Float64(N_trans_switch),
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
