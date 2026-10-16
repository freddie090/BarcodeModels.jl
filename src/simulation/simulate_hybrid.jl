function _simulate_experiment_hybrid(model::ResPop, exp::ExperimentParams; kwargs...)
    params = model.params
    R_real = normalize_R_real(_kw(kwargs, :R_real, "l"))
    save_at = _kw(kwargs, :save_at, exp.save_at)
    n_rep = _kw(kwargs, :n_rep, exp.n_rep)
    drug_treatment = _kw(kwargs, :drug_treatment, exp.drug_treatment)
    full_sol = _kw(kwargs, :full_sol, exp.full_sol)
    de = normalize_respop_drug_effect(_kw(kwargs, :drug_effect, params.drug_effect))

    _validate_tmax_vector_constraints(exp.tmax, exp.t_Pass)
    _validate_tmax_length(exp.tmax, n_rep)

    normalize_t_pass(times) = begin
        if times isa AbstractVector
            Vector{Float64}(times)
        else
            t_val = Float64(times)
            t_val < 0.0 ? Float64[] : [t_val]
        end
    end
    n_pass_eff = length(normalize_t_pass(exp.t_Pass)) + 1

    n_pop_obsv = n_rep * (length(exp.t_keep) + n_pass_eff)
    validate_resistant_birth_drug_effect(params, R_real; drug_effect = de)
    model_eff = _with_drug_effect(model, de)

    if params.psi < 0.0
        if full_sol
            error("psi must be between 0 and 1.")
        end
        return Dict("t" => repeat([-1.0], n_pop_obsv), "u" => repeat([-1.0], n_pop_obsv))
    end

    if (exp.t_exp isa AbstractVector) && (exp.Nseed isa AbstractVector)
        @assert length(exp.t_exp) == length(exp.Nseed) "Length of t_exp and Nseed vectors must match."

        nR = Float64(round(params.rho * exp.n0))
        nS = exp.n0 - nR
        nE = 0.0

        for i in 1:(length(exp.t_exp) - 1)
            N_total = Int64(round(nS + nR + nE))
            exp_sim = SimParams(
                n0 = N_total,
                t0 = 0.0,
                tmax = exp.t_exp[i],
                t_Pass = -1.0,
                Nmax = exp.Nmax,
                Cc = exp.Cc,
                Nswitch = exp.Nswitch,
                Aswitch = exp.Aswitch,
                R_real = R_real,
                treat_ons = exp.treat_ons,
                treat_offs = exp.treat_offs,
                save_at = save_at,
                treat = false
            )
            exp_state = ResPopState(nS, nR, nE; gam = 0.0, pass_num = 1)
            exp_sol = run_model_core_hybrid(model_eff, exp_state, exp_sim; treat = false)

            nS = Int64(round(respop_pop_fun(last(exp_sol.u))[1]))
            nR = Int64(round(respop_pop_fun(last(exp_sol.u))[2]))
            nE = Int64(round(respop_pop_fun(last(exp_sol.u))[3]))

            total_cells = nS + nR + nE
            if total_cells < exp.Nseed[i]
                return Dict("t" => repeat([-1.0], n_pop_obsv), "u" => repeat([-1.0], n_pop_obsv))
            end

            counts = multivariate_hypergeometric_draw([nS, nR, nE], exp.Nseed[i])
            nS = Float64(counts[1])
            nR = Float64(counts[2])
            nE = Float64(counts[3])
        end

        N_total = Int64(round(nS + nR + nE))
        exp_sim = SimParams(
            n0 = N_total,
            t0 = 0.0,
            tmax = exp.t_exp[end],
            t_Pass = -1.0,
            Nmax = exp.Nmax,
            Cc = exp.Cc,
            Nswitch = exp.Nswitch,
            Aswitch = exp.Aswitch,
            R_real = R_real,
            treat_ons = exp.treat_ons,
            treat_offs = exp.treat_offs,
            save_at = save_at,
            treat = false
        )
        exp_state = ResPopState(nS, nR, nE; gam = 0.0, pass_num = 1)
        exp_sol = run_model_core_hybrid(model_eff, exp_state, exp_sim; treat = false)

        exp_nS = Int64(round(respop_pop_fun(last(exp_sol.u))[1]))
        exp_nR = Int64(round(respop_pop_fun(last(exp_sol.u))[2]))
        exp_nE = Int64(round(respop_pop_fun(last(exp_sol.u))[3]))
    elseif (exp.t_exp isa Real) && (exp.Nseed isa Integer)
        nR = Float64(round(params.rho * exp.n0))
        nS = exp.n0 - nR
        nE = 0.0

        exp_sim = SimParams(
            n0 = exp.n0,
            t0 = 0.0,
            tmax = exp.t_exp,
            t_Pass = -1.0,
            Nmax = exp.Nmax,
            Cc = exp.Cc,
            Nswitch = exp.Nswitch,
            Aswitch = exp.Aswitch,
            R_real = R_real,
            treat_ons = exp.treat_ons,
            treat_offs = exp.treat_offs,
            save_at = save_at,
            treat = false
        )
        exp_state = ResPopState(nS, nR, nE; gam = 0.0, pass_num = 1)
        exp_sol = run_model_core_hybrid(model_eff, exp_state, exp_sim; treat = false)

        exp_nS = Int64(round(respop_pop_fun(last(exp_sol.u))[1]))
        exp_nR = Int64(round(respop_pop_fun(last(exp_sol.u))[2]))
        exp_nE = Int64(round(respop_pop_fun(last(exp_sol.u))[3]))
    else
        error("t_exp and Nseed must either both be scalars or both be vectors of the same length.")
    end

    nseed_last = _nseed_last(exp.Nseed)
    if (exp_nS + exp_nR + exp_nE) < (nseed_last * n_rep)
        return Dict("t" => repeat([-1.0], n_pop_obsv), "u" => repeat([-1.0], n_pop_obsv))
    end

    rep_phenos = Array{Float64}(undef, 3, n_rep)
    nS = exp_nS
    nR = exp_nR
    nE = exp_nE
    for i in 1:n_rep
        sampled_counts = multivariate_hypergeometric_draw([nS, nR, nE], nseed_last)
        rep_phenos[1, i] = Float64(sampled_counts[1])
        rep_phenos[2, i] = Float64(sampled_counts[2])
        rep_phenos[3, i] = Float64(sampled_counts[3])
        nS -= sampled_counts[1]
        nR -= sampled_counts[2]
        nE -= sampled_counts[3]
    end

    fin_t_outs = Float64[]
    fin_u_outs = Float64[]
    sol_dfs = DataFrame[]

    for i in 1:n_rep
        rep_tmax = _replicate_tmax(exp.tmax, n_rep, i)

        sim = SimParams(
            n0 = nseed_last,
            t0 = 0.0,
            tmax = rep_tmax,
            t_Pass = exp.t_Pass,
            Nmax = exp.Nmax,
            Cc = exp.Cc,
            Nswitch = exp.Nswitch,
            Aswitch = exp.Aswitch,
            R_real = R_real,
            treat_ons = exp.treat_ons,
            treat_offs = exp.treat_offs,
            save_at = save_at,
            treat = drug_treatment,
        )
        state = ResPopState(rep_phenos[1, i], rep_phenos[2, i], rep_phenos[3, i];
                           gam = 0.0, pass_num = 1)
        sol = run_model_core_hybrid(model_eff, state, sim; treat = drug_treatment)

        if full_sol
            sol_df = DataFrame(
                t = sol.t,
                nS = map(x -> x[RESPOP_NS_INDEX], sol.u),
                nR = map(x -> x[RESPOP_NR_INDEX], sol.u),
                nE = map(x -> x[RESPOP_NE_INDEX], sol.u),
                N = map(x -> x[RESPOP_NS_INDEX], sol.u) .+
                    map(x -> x[RESPOP_NR_INDEX], sol.u) .+
                    map(x -> x[RESPOP_NE_INDEX], sol.u),
                rep = i
            )
            push!(sol_dfs, sol_df)
        end

        total_pop = map(u -> sum(respop_pop_fun(u)), sol.u)
        t_outs = Vector{Float64}(undef, length(exp.t_keep) + n_pass_eff)
        u_outs = Vector{Float64}(undef, length(exp.t_keep) + n_pass_eff)

        for j in 1:length(exp.t_keep)
            if !(exp.t_keep[j] in sol.t)
                t_closest_pos = findmin(abs.(sol.t .- exp.t_keep[j]))[2]
                t_outs[j] = sol.t[t_closest_pos]
                u_outs[j] = Float64(round(total_pop[t_closest_pos]))
            else
                t_realised_pos = findlast(sol.t .== exp.t_keep[j])
                t_outs[j] = sol.t[t_realised_pos]
                u_outs[j] = Float64(round(total_pop[t_realised_pos]))
            end
        end

        for j in 1:n_pass_eff
            if sum(respop_pass_fun.(sol.u) .== j) > 0
                t_realised_pos = findlast(respop_pass_fun.(sol.u) .== j)
                idx = length(exp.t_keep) + j
                t_outs[idx] = sol.t[t_realised_pos]
                u_outs[idx] = Float64(round(total_pop[t_realised_pos]))
            else
                idx = length(exp.t_keep) + j
                t_outs[idx] = t_outs[idx - 1]
                u_outs[idx] = u_outs[idx - 1]
            end
        end

        append!(fin_t_outs, t_outs)
        append!(fin_u_outs, u_outs)
    end

    fin_t_outs = Float64.(round.(fin_t_outs, digits = 0))
    if full_sol
        sol_df = isempty(sol_dfs) ? DataFrame() : vcat(sol_dfs...)
        return Dict("t" => fin_t_outs, "u" => fin_u_outs, "sol_df" => sol_df)
    end
    return Dict("t" => fin_t_outs, "u" => fin_u_outs)
end

function _simulate_experiment_hybrid(model::ResDmg, exp::ExperimentParams; kwargs...)
    params = model.params
    R_real = normalize_R_real(_kw(kwargs, :R_real, "l"))
    save_at = _kw(kwargs, :save_at, exp.save_at)
    n_rep = _kw(kwargs, :n_rep, exp.n_rep)
    drug_treatment = _kw(kwargs, :drug_treatment, exp.drug_treatment)
    full_sol = _kw(kwargs, :full_sol, exp.full_sol)
    de = normalize_resdmg_drug_effect(_kw(kwargs, :drug_effect, params.drug_effect))

    _validate_tmax_vector_constraints(exp.tmax, exp.t_Pass)
    _validate_tmax_length(exp.tmax, n_rep)

    normalize_t_pass(times) = begin
        if times isa AbstractVector
            Vector{Float64}(times)
        else
            t_val = Float64(times)
            t_val < 0.0 ? Float64[] : [t_val]
        end
    end
    n_pass_eff = length(normalize_t_pass(exp.t_Pass)) + 1

    n_pop_obsv = n_rep * (length(exp.t_keep) + n_pass_eff)
    validate_resistant_birth_drug_effect(params, R_real; drug_effect = de)
    model_eff = _with_drug_effect(model, de)

    if params.psi < 0.0
        if full_sol
            error("psi must be between 0 and 1.")
        end
        return Dict("t" => repeat([-1.0], n_pop_obsv), "u" => repeat([-1.0], n_pop_obsv))
    end

    if (params.zet_S > params.zet_R) && (params.zet_R != 0.0)
        if full_sol
            error("zet_R must be >= zet_S for resistant cells to have at least as high repair propensity as sensitive cells.")
        end
        return Dict("t" => repeat([-1.0], n_pop_obsv), "u" => repeat([-1.0], n_pop_obsv))
    end

    if (exp.t_exp isa AbstractVector) && (exp.Nseed isa AbstractVector)
        @assert length(exp.t_exp) == length(exp.Nseed) "Length of t_exp and Nseed vectors must match."

        nR = Float64(round(params.rho * exp.n0))
        nS = exp.n0 - nR
        nDS = 0.0
        nDR = 0.0

        for i in 1:(length(exp.t_exp) - 1)
            N_total = Int64(round(nS + nDS + nDR + nR))
            exp_sim = SimParams(
                n0 = N_total,
                t0 = 0.0,
                tmax = exp.t_exp[i],
                t_Pass = -1.0,
                Nmax = exp.Nmax,
                Cc = exp.Cc,
                Nswitch = exp.Nswitch,
                Aswitch = exp.Aswitch,
                R_real = R_real,
                treat_ons = exp.treat_ons,
                treat_offs = exp.treat_offs,
                save_at = save_at,
                treat = false
            )
            exp_state = ResDmgState(nS, nDS, nDR, nR; gam = 0.0, pass_num = 1)
            exp_sol = run_model_core_hybrid(model_eff, exp_state, exp_sim; treat = false)

            nS = Int64(round(resdmg_pop_fun(last(exp_sol.u))[1]))
            nDS = Int64(round(resdmg_pop_fun(last(exp_sol.u))[2]))
            nDR = Int64(round(resdmg_pop_fun(last(exp_sol.u))[3]))
            nR = Int64(round(resdmg_pop_fun(last(exp_sol.u))[4]))

            total_cells = nS + nDS + nDR + nR
            if total_cells < exp.Nseed[i]
                return Dict("t" => repeat([-1.0], n_pop_obsv), "u" => repeat([-1.0], n_pop_obsv))
            end

            counts = multivariate_hypergeometric_draw([nS, nDS, nDR, nR], exp.Nseed[i])
            nS = Float64(counts[1])
            nDS = Float64(counts[2])
            nDR = Float64(counts[3])
            nR = Float64(counts[4])
        end

        N_total = Int64(round(nS + nDS + nDR + nR))
        exp_sim = SimParams(
            n0 = N_total,
            t0 = 0.0,
            tmax = exp.t_exp[end],
            t_Pass = -1.0,
            Nmax = exp.Nmax,
            Cc = exp.Cc,
            Nswitch = exp.Nswitch,
            Aswitch = exp.Aswitch,
            R_real = R_real,
            treat_ons = exp.treat_ons,
            treat_offs = exp.treat_offs,
            save_at = save_at,
            treat = false
        )
        exp_state = ResDmgState(nS, nDS, nDR, nR; gam = 0.0, pass_num = 1)
        exp_sol = run_model_core_hybrid(model_eff, exp_state, exp_sim; treat = false)

        exp_nS = Int64(round(resdmg_pop_fun(last(exp_sol.u))[1]))
        exp_nDS = Int64(round(resdmg_pop_fun(last(exp_sol.u))[2]))
        exp_nDR = Int64(round(resdmg_pop_fun(last(exp_sol.u))[3]))
        exp_nR = Int64(round(resdmg_pop_fun(last(exp_sol.u))[4]))
    elseif (exp.t_exp isa Real) && (exp.Nseed isa Integer)
        nR = Float64(round(params.rho * exp.n0))
        nS = exp.n0 - nR
        nDS = 0.0
        nDR = 0.0

        exp_sim = SimParams(
            n0 = exp.n0,
            t0 = 0.0,
            tmax = exp.t_exp,
            t_Pass = -1.0,
            Nmax = exp.Nmax,
            Cc = exp.Cc,
            Nswitch = exp.Nswitch,
            Aswitch = exp.Aswitch,
            R_real = R_real,
            treat_ons = exp.treat_ons,
            treat_offs = exp.treat_offs,
            save_at = save_at,
            treat = false
        )
        exp_state = ResDmgState(nS, nDS, nDR, nR; gam = 0.0, pass_num = 1)
        exp_sol = run_model_core_hybrid(model_eff, exp_state, exp_sim; treat = false)

        exp_nS = Int64(round(resdmg_pop_fun(last(exp_sol.u))[1]))
        exp_nDS = Int64(round(resdmg_pop_fun(last(exp_sol.u))[2]))
        exp_nDR = Int64(round(resdmg_pop_fun(last(exp_sol.u))[3]))
        exp_nR = Int64(round(resdmg_pop_fun(last(exp_sol.u))[4]))
    else
        error("t_exp and Nseed must either both be scalars or both be vectors of the same length.")
    end

    nseed_last = _nseed_last(exp.Nseed)
    if (exp_nS + exp_nDS + exp_nDR + exp_nR) < (nseed_last * n_rep)
        return Dict("t" => repeat([-1.0], n_pop_obsv), "u" => repeat([-1.0], n_pop_obsv))
    end

    rep_phenos = Array{Float64}(undef, 4, n_rep)
    nS = exp_nS
    nDS = exp_nDS
    nDR = exp_nDR
    nR = exp_nR
    for i in 1:n_rep
        sampled_counts = multivariate_hypergeometric_draw([nS, nDS, nDR, nR], nseed_last)
        rep_phenos[1, i] = Float64(sampled_counts[1])
        rep_phenos[2, i] = Float64(sampled_counts[2])
        rep_phenos[3, i] = Float64(sampled_counts[3])
        rep_phenos[4, i] = Float64(sampled_counts[4])
        nS -= sampled_counts[1]
        nDS -= sampled_counts[2]
        nDR -= sampled_counts[3]
        nR -= sampled_counts[4]
    end

    fin_t_outs = Float64[]
    fin_u_outs = Float64[]
    sol_dfs = DataFrame[]

    for i in 1:n_rep
        rep_tmax = _replicate_tmax(exp.tmax, n_rep, i)

        sim = SimParams(
            n0 = nseed_last,
            t0 = 0.0,
            tmax = rep_tmax,
            t_Pass = exp.t_Pass,
            Nmax = exp.Nmax,
            Cc = exp.Cc,
            Nswitch = exp.Nswitch,
            Aswitch = exp.Aswitch,
            R_real = R_real,
            treat_ons = exp.treat_ons,
            treat_offs = exp.treat_offs,
            save_at = save_at,
            treat = drug_treatment,
        )
        state = ResDmgState(rep_phenos[1, i], rep_phenos[2, i], rep_phenos[3, i], rep_phenos[4, i];
                            gam = 0.0, pass_num = 1)
        sol = run_model_core_hybrid(model_eff, state, sim; treat = drug_treatment)

        if full_sol
            sol_df = DataFrame(
                t = sol.t,
                nS = map(x -> x[RESDMG_NS_INDEX], sol.u),
                nDS = map(x -> x[RESDMG_NDS_INDEX], sol.u),
                nDR = map(x -> x[RESDMG_NDR_INDEX], sol.u),
                nR = map(x -> x[RESDMG_NR_INDEX], sol.u),
                N = map(x -> x[RESDMG_NS_INDEX], sol.u) .+
                    map(x -> x[RESDMG_NDS_INDEX], sol.u) .+
                    map(x -> x[RESDMG_NDR_INDEX], sol.u) .+
                    map(x -> x[RESDMG_NR_INDEX], sol.u),
                rep = i
            )
            push!(sol_dfs, sol_df)
        end

        total_pop = map(u -> sum(resdmg_pop_fun(u)), sol.u)
        t_outs = Vector{Float64}(undef, length(exp.t_keep) + n_pass_eff)
        u_outs = Vector{Float64}(undef, length(exp.t_keep) + n_pass_eff)

        for j in 1:length(exp.t_keep)
            if !(exp.t_keep[j] in sol.t)
                t_closest_pos = findmin(abs.(sol.t .- exp.t_keep[j]))[2]
                t_outs[j] = sol.t[t_closest_pos]
                u_outs[j] = Float64(round(total_pop[t_closest_pos]))
            else
                t_realised_pos = findlast(sol.t .== exp.t_keep[j])
                t_outs[j] = sol.t[t_realised_pos]
                u_outs[j] = Float64(round(total_pop[t_realised_pos]))
            end
        end

        for j in 1:n_pass_eff
            if sum(resdmg_pass_fun.(sol.u) .== j) > 0
                t_realised_pos = findlast(resdmg_pass_fun.(sol.u) .== j)
                idx = length(exp.t_keep) + j
                t_outs[idx] = sol.t[t_realised_pos]
                u_outs[idx] = Float64(round(total_pop[t_realised_pos]))
            else
                idx = length(exp.t_keep) + j
                t_outs[idx] = t_outs[idx - 1]
                u_outs[idx] = u_outs[idx - 1]
            end
        end

        append!(fin_t_outs, t_outs)
        append!(fin_u_outs, u_outs)
    end

    fin_t_outs = Float64.(round.(fin_t_outs, digits = 0))
    if full_sol
        sol_df = isempty(sol_dfs) ? DataFrame() : vcat(sol_dfs...)
        return Dict("t" => fin_t_outs, "u" => fin_u_outs, "sol_df" => sol_df)
    end
    return Dict("t" => fin_t_outs, "u" => fin_u_outs)
end

function _simulate_simple_hybrid(model::ResPop, sim::SimpleSimParams; kwargs...)
    params = model.params
    R_real = normalize_R_real(_kw(kwargs, :R_real, "l"))
    save_at = _kw(kwargs, :save_at, sim.save_at)
    drug_treatment = _kw(kwargs, :drug_treatment, sim.drug_treatment)
    de = normalize_respop_drug_effect(_kw(kwargs, :drug_effect, params.drug_effect))

    validate_resistant_birth_drug_effect(params, R_real; drug_effect = de)
    model_eff = _with_drug_effect(model, de)

    nR = Float64(round(params.rho * sim.n0))
    nS = sim.n0 - nR
    nE = 0.0

    sim_eff = SimParams(
        n0 = sim.n0,
        t0 = 0.0,
        tmax = sim.tmax,
        t_Pass = -1.0,
        Nmax = sim.Nmax,
        Cc = sim.Cc,
        Nswitch = sim.Nswitch,
        Aswitch = sim.Aswitch,
        R_real = R_real,
        treat_ons = sim.treat_ons,
        treat_offs = sim.treat_offs,
        save_at = save_at,
        treat = drug_treatment,
    )
    state = ResPopState(nS, nR, nE; gam = 0.0, pass_num = 1)
    sol = run_model_core_hybrid(model_eff, state, sim_eff; treat = drug_treatment)

    sol_df = DataFrame(
        t = sol.t,
        nS = map(x -> x[RESPOP_NS_INDEX], sol.u),
        nR = map(x -> x[RESPOP_NR_INDEX], sol.u),
        nE = map(x -> x[RESPOP_NE_INDEX], sol.u),
        N = map(x -> x[RESPOP_NS_INDEX], sol.u) .+
            map(x -> x[RESPOP_NR_INDEX], sol.u) .+
            map(x -> x[RESPOP_NE_INDEX], sol.u)
    )

    return Dict("sol_df" => sol_df)
end

function _simulate_simple_hybrid(model::ResDmg, sim::SimpleSimParams; kwargs...)
    params = model.params
    R_real = normalize_R_real(_kw(kwargs, :R_real, "l"))
    save_at = _kw(kwargs, :save_at, sim.save_at)
    drug_treatment = _kw(kwargs, :drug_treatment, sim.drug_treatment)
    de = normalize_resdmg_drug_effect(_kw(kwargs, :drug_effect, params.drug_effect))

    validate_resistant_birth_drug_effect(params, R_real; drug_effect = de)
    model_eff = _with_drug_effect(model, de)

    nR = Float64(round(params.rho * sim.n0))
    nS = sim.n0 - nR
    nDS = 0.0
    nDR = 0.0

    sim_eff = SimParams(
        n0 = sim.n0,
        t0 = 0.0,
        tmax = sim.tmax,
        t_Pass = -1.0,
        Nmax = sim.Nmax,
        Cc = sim.Cc,
        Nswitch = sim.Nswitch,
        Aswitch = sim.Aswitch,
        R_real = R_real,
        treat_ons = sim.treat_ons,
        treat_offs = sim.treat_offs,
        save_at = save_at,
        treat = drug_treatment,
    )
    state = ResDmgState(nS, nDS, nDR, nR; gam = 0.0, pass_num = 1)
    sol = run_model_core_hybrid(model_eff, state, sim_eff; treat = drug_treatment)

    sol_df = DataFrame(
        t = sol.t,
        nS = map(x -> x[RESDMG_NS_INDEX], sol.u),
        nDS = map(x -> x[RESDMG_NDS_INDEX], sol.u),
        nDR = map(x -> x[RESDMG_NDR_INDEX], sol.u),
        nR = map(x -> x[RESDMG_NR_INDEX], sol.u),
        N = map(x -> x[RESDMG_NS_INDEX], sol.u) .+
            map(x -> x[RESDMG_NDS_INDEX], sol.u) .+
            map(x -> x[RESDMG_NDR_INDEX], sol.u) .+
            map(x -> x[RESDMG_NR_INDEX], sol.u)
    )

    return Dict("sol_df" => sol_df)
end

"""Initialize exact resistance/EG marginal counts with a random independent overlap."""
function _seed_invivo_hybrid_counts(n0::Int64, rho::Float64, fEG1::Float64)
    nR = round(Int, n0*rho)
    counts = [n0-nR, nR, 0]
    eg1 = multivariate_hypergeometric_draw(counts, round(Int, n0*fEG1))
    return counts-eg1, eg1
end

function _invivo_state(eg0, eg1; gam = 0.0, pass_num = 1)
    return ResPopInVivoState(eg0[1], eg1[1], eg0[2], eg1[2], eg0[3], eg1[3]; gam, pass_num)
end

function _invivo_counts_from_u(u)
    return round.(Int, u[[2, 4, 6]]), round.(Int, u[[3, 5, 7]])
end

function _invivo_hybrid_sol_df(sol, t_offset::Float64; cond::String, rep::Int64, passage = nothing)
    df = DataFrame(t = sol.t .+ t_offset,
        nS_EG0 = [u[2] for u in sol.u], nS_EG1 = [u[3] for u in sol.u],
        nR_EG0 = [u[4] for u in sol.u], nR_EG1 = [u[5] for u in sol.u],
        nE_EG0 = [u[6] for u in sol.u], nE_EG1 = [u[7] for u in sol.u])
    df.nS = df.nS_EG0 + df.nS_EG1
    df.nR = df.nR_EG0 + df.nR_EG1
    df.nE = df.nE_EG0 + df.nE_EG1
    df.n_EG0 = df.nS_EG0 + df.nR_EG0 + df.nE_EG0
    df.n_EG1 = df.nS_EG1 + df.nR_EG1 + df.nE_EG1
    df.N = df.n_EG0 + df.n_EG1
    df[!, :cond] = fill(cond, nrow(df))
    df[!, :rep] = fill(rep, nrow(df))
    df.passage = passage === nothing ? [round(Int, u[8]) for u in sol.u] : fill(passage, length(sol.t))
    return df
end

function _run_invivo_hybrid_segment(model::ResPopInVivo, eg0, eg1, duration, exp::ExperimentParams;
    R_real = "l", save_at = exp.save_at)
    state = _invivo_state(eg0, eg1; pass_num = 0)
    duration == 0 && return (t = [0.0], u = [to_componentarray(state)])
    sim = SimParams(n0 = sum(eg0)+sum(eg1), tmax = duration, t_Pass = Float64[],
        Nmax = exp.Nmax, Cc = exp.Cc, Nswitch = exp.Nswitch, Aswitch = exp.Aswitch,
        R_real = R_real, save_at = save_at, treat_ons = Float64[], treat_offs = Float64[])
    # Passage zero is an output label, not an available extra passage.
    state = _invivo_state(eg0, eg1)
    return run_model_core_hybrid(model, state, sim; treat = false)
end

function _simulate_experiment_hybrid(model::ResPopInVivo, exp::ExperimentParams; kwargs...)
    R_real = normalize_R_real(_kw(kwargs, :R_real, "l"))
    n_rep = Int64(_kw(kwargs, :n_rep, exp.n_rep))
    n_rep > 0 || error("n_rep must be > 0.")
    treat = Bool(_kw(kwargs, :drug_treatment, exp.drug_treatment))
    inc_control = Bool(_kw(kwargs, :inc_control, exp.inc_control))
    inc_pot = Bool(_kw(kwargs, :inc_pot, exp.inc_pot))
    full_sol = Bool(_kw(kwargs, :full_sol, exp.full_sol))
    save_at = _kw(kwargs, :save_at, exp.save_at)
    de = normalize_respop_drug_effect(_kw(kwargs, :drug_effect, model.params.drug_effect))
    _validate_tmax_vector_constraints(exp.tmax, exp.t_Pass)
    _validate_tmax_length(exp.tmax, n_rep)
    validate_resistant_birth_drug_effect(model.params, R_real; drug_effect = de)
    model_eff = _with_drug_effect(model, de)
    if model.params.psi < 0
        full_sol && error("psi must be between 0 and 1.")
        return Dict("t" => [-1.0], "u" => [-1.0])
    end
    design = _experiment_condition_design(treat, inc_control, n_rep)
    expansion = ResPopInVivo(_copy_respop_invivo_params(model_eff.params; al = 0.0))
    eg0, eg1 = _seed_invivo_hybrid_counts(exp.n0, model.params.rho, model.params.fEG1)
    if exp.t_exp isa AbstractVector && exp.Nseed isa AbstractVector
        length(exp.t_exp) == length(exp.Nseed) && !isempty(exp.t_exp) ||
            error("t_exp and Nseed vectors must be nonempty and have matching lengths.")
        durations, seeds = exp.t_exp, exp.Nseed
    elseif exp.t_exp isa Real && exp.Nseed isa Integer
        durations, seeds = [exp.t_exp], [exp.Nseed]
    else
        error("t_exp and Nseed must both be scalars or both be vectors of matching lengths.")
    end
    dfs = DataFrame[]
    offset = 0.0
    for i in eachindex(durations)
        durations[i] >= 0 || error("Expansion durations must be nonnegative.")
        seg = _run_invivo_hybrid_segment(expansion, eg0, eg1, durations[i], exp; R_real, save_at)
        inc_pot && push!(dfs, _invivo_hybrid_sol_df(seg, offset; cond = "POT", rep = Int64(0), passage = 0))
        offset += last(seg.t)
        eg0, eg1 = _invivo_counts_from_u(last(seg.u))
        if i != lastindex(durations)
            sum(eg0)+sum(eg1) >= seeds[i] || return Dict("t" => [-1.0], "u" => [-1.0])
            draw = multivariate_hypergeometric_draw(vcat(eg0, eg1), seeds[i])
            eg0, eg1 = draw[1:3], draw[4:6]
        end
    end
    nseed = Int64(last(seeds))
    pool = vcat(eg0, eg1)
    sum(pool) >= nseed*length(design) || return Dict("t" => [-1.0], "u" => [-1.0])
    grafts = NamedTuple[]
    reps = NamedTuple[]
    for rec in design
        draw = multivariate_hypergeometric_draw(pool, nseed)
        pool -= draw
        a, b, stats = engraftment_selection(draw[1:3], draw[4:6], model.params.pEG, model.params.sEG)
        _push_invivo_engraftment!(grafts, rec.cond, rec.rep, 1, stats)
        push!(reps, (design = rec, eg0 = a, eg1 = b))
    end
    times, totals, conditions, replicates = Float64[], Float64[], String[], Int64[]
    for batch in reps
        rec = batch.design
        state = _invivo_state(batch.eg0, batch.eg1)
        sim = SimParams(n0 = nseed, tmax = _replicate_tmax(exp.tmax, n_rep, rec.rep),
            t_Pass = exp.t_Pass, Nmax = exp.Nmax, Cc = exp.Cc, Nswitch = exp.Nswitch,
            Aswitch = exp.Aswitch, R_real = R_real, save_at = save_at, treat = rec.treat,
            treat_ons = exp.treat_ons, treat_offs = exp.treat_offs)
        sol = if sum(batch.eg0)+sum(batch.eg1) == 0
            (t = [0.0], u = [to_componentarray(state)])
        else
            on_engraft = (passage, stats) -> _push_invivo_engraftment!(grafts, rec.cond, rec.rep, passage, stats)
            run_model_core_hybrid(model_eff, state, sim; on_engraft)
        end
        full_sol && push!(dfs, _invivo_hybrid_sol_df(sol, 0.0; cond = rec.cond, rep = rec.rep))
        push!(times, last(sol.t)); push!(totals, respop_invivo_total_population(last(sol.u)))
        push!(conditions, rec.cond); push!(replicates, rec.rep)
    end
    out = Dict{String, Any}("t" => times, "u" => totals, "cond" => conditions, "rep" => replicates,
        "engraft_df" => _invivo_engraftment_df(grafts))
    (full_sol || inc_pot) && (out["sol_df"] = isempty(dfs) ? DataFrame() : vcat(dfs...))
    return out
end

function _simulate_simple_hybrid(model::ResPopInVivo, settings::SimpleSimParams; kwargs...)
    R_real = normalize_R_real(_kw(kwargs, :R_real, "l"))
    de = normalize_respop_drug_effect(_kw(kwargs, :drug_effect, model.params.drug_effect))
    validate_model_params_strict(model.params)
    validate_resistant_birth_drug_effect(model.params, R_real; drug_effect = de)
    model_eff = _with_drug_effect(model, de)
    eg0, eg1 = _seed_invivo_hybrid_counts(settings.n0, model.params.rho, model.params.fEG1)
    sim = SimParams(n0 = settings.n0, tmax = settings.tmax, Nmax = settings.Nmax, Cc = settings.Cc,
        Nswitch = settings.Nswitch, Aswitch = settings.Aswitch, R_real = R_real,
        save_at = _kw(kwargs, :save_at, settings.save_at),
        treat = Bool(_kw(kwargs, :drug_treatment, settings.drug_treatment)),
        treat_ons = settings.treat_ons, treat_offs = settings.treat_offs)
    sol = run_model_core_hybrid(model_eff, _invivo_state(eg0, eg1), sim)
    df = _invivo_hybrid_sol_df(sol, 0.0; cond = "DT", rep = Int64(1))
    select!(df, Not([:cond, :rep, :passage]))
    return Dict("sol_df" => df)
end
