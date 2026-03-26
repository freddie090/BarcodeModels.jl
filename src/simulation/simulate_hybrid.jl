function _simulate_experiment_hybrid(model::ResPop, exp::ExperimentParams; kwargs...)
    params = model.params
    save_at = _kw(kwargs, :save_at, exp.save_at)
    n_rep = _kw(kwargs, :n_rep, exp.n_rep)
    drug_treatment = _kw(kwargs, :drug_treatment, exp.drug_treatment)
    full_sol = _kw(kwargs, :full_sol, exp.full_sol)
    run_IC = _kw(kwargs, :run_IC, exp.run_IC)
    IC_n0 = _kw(kwargs, :IC_n0, exp.IC_n0)
    IC_tmax = _kw(kwargs, :IC_tmax, exp.IC_tmax)
    IC_treat_on = _kw(kwargs, :IC_treat_on, exp.IC_treat_on)
    run_colony = _kw(kwargs, :run_colony, exp.run_colony)
    nCol = _kw(kwargs, :nCol, exp.nCol)
    tCol = _kw(kwargs, :tCol, exp.tCol)
    ColNmax = _kw(kwargs, :ColNmax, exp.ColNmax)
    de = normalize_respop_drug_effect(_kw(kwargs, :drug_effect, params.drug_effect))

    @assert !(run_IC && run_colony) "Cannot run IC and colony assays at the same time."
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

    n_pop_obsv = n_rep * ((run_IC * 2) + length(exp.t_keep) + n_pass_eff) + (run_colony * 1)
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
                N_trans_switch = exp.N_trans_switch,
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
            N_trans_switch = exp.N_trans_switch,
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
            N_trans_switch = exp.N_trans_switch,
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

    if run_colony
        tx1_cols = 0
        tx0_cols = 0

        col_nR = max(Float64(round(params.rho * exp.n0)), 0.0)
        col_nS = max((exp.n0 - col_nR), 0.0)
        col_nE = 0.0

        for j in 1:nCol
            total_cells = max(col_nS + col_nR + col_nE, 1e-10)
            cell_probs = [col_nS, col_nR, col_nE] ./ total_cells
            if any(cell_probs .< 0) || sum(cell_probs) == 0
                error("Invalid probability vector: $(cell_probs)")
            end
            sampled_phenos = zeros(Float64, 3)
            sampled_phenos[rand(Categorical(cell_probs))] = 1.0

            col_sim = SimParams(
                n0 = 1,
                t0 = 0.0,
                tmax = tCol,
                t_Pass = -1.0,
                Nmax = ColNmax + 10,
                Cc = exp.Cc,
                Nswitch = exp.Nswitch,
                N_trans_switch = exp.N_trans_switch,
                treat_ons = [1.0],
                treat_offs = [100.0],
                save_at = save_at,
                treat = true
            )
            col_state = ResPopState(sampled_phenos[1], sampled_phenos[2], sampled_phenos[3];
                                   gam = 0.0, pass_num = 1)
            col_sol_tx1 = run_model_core_hybrid(model_eff, col_state, col_sim; treat = true)
            col_sol_tx0 = run_model_core_hybrid(model_eff, col_state, col_sim; treat = false)

            if round(sum(last(respop_pop_fun.(col_sol_tx1.u)))) >= ColNmax
                tx1_cols += 1
            end
            if round(sum(last(respop_pop_fun.(col_sol_tx0.u)))) >= ColNmax
                tx0_cols += 1
            end
        end

        col_prop = tx0_cols == 0 ? 0.0 : (tx1_cols / tx0_cols)
        append!(fin_t_outs, tCol)
        append!(fin_u_outs, col_prop)
    end

    for i in 1:n_rep
        rep_tmax = _replicate_tmax(exp.tmax, n_rep, i)

        if run_IC
            ic_nR = Float64(round(params.rho * IC_n0))
            ic_nS = IC_n0 - ic_nR
            ic_nE = 0.0
            ic_state = ResPopState(ic_nS, ic_nR, ic_nE; gam = 0.0, pass_num = 1)
            ic_sim = SimParams(
                n0 = IC_n0,
                t0 = 0.0,
                tmax = IC_tmax,
                t_Pass = -1.0,
                Nmax = exp.Nmax,
                Cc = exp.Cc,
                Nswitch = exp.Nswitch,
                N_trans_switch = exp.N_trans_switch,
                treat_ons = [IC_treat_on],
                treat_offs = [100.0],
                save_at = save_at,
                treat = true
            )
            IC_sol_tx1 = run_model_core_hybrid(model_eff, ic_state, ic_sim; treat = true)
            IC_sol_tx0 = run_model_core_hybrid(model_eff, ic_state, ic_sim; treat = false)
        end

        sim = SimParams(
            n0 = nseed_last,
            t0 = 0.0,
            tmax = rep_tmax,
            t_Pass = exp.t_Pass,
            Nmax = exp.Nmax,
            Cc = exp.Cc,
            Nswitch = exp.Nswitch,
            N_trans_switch = exp.N_trans_switch,
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
        t_outs = Vector{Float64}(undef, (run_IC * 2) + length(exp.t_keep) + n_pass_eff)
        u_outs = Vector{Float64}(undef, (run_IC * 2) + length(exp.t_keep) + n_pass_eff)

        if run_IC
            t_outs[1] = last(IC_sol_tx1.t)
            u_outs[1] = Float64(round(sum(last(respop_pop_fun.(IC_sol_tx1.u)))))
            t_outs[2] = last(IC_sol_tx0.t)
            u_outs[2] = Float64(round(sum(last(respop_pop_fun.(IC_sol_tx0.u)))))
        end

        for j in 1:length(exp.t_keep)
            if !(exp.t_keep[j] in sol.t)
                t_closest_pos = findmin(abs.(sol.t .- exp.t_keep[j]))[2]
                t_outs[j + (run_IC * 2)] = sol.t[t_closest_pos]
                u_outs[j + (run_IC * 2)] = Float64(round(total_pop[t_closest_pos]))
            else
                t_realised_pos = findlast(sol.t .== exp.t_keep[j])
                t_outs[j + (run_IC * 2)] = sol.t[t_realised_pos]
                u_outs[j + (run_IC * 2)] = Float64(round(total_pop[t_realised_pos]))
            end
        end

        for j in 1:n_pass_eff
            if sum(respop_pass_fun.(sol.u) .== j) > 0
                t_realised_pos = findlast(respop_pass_fun.(sol.u) .== j)
                idx = length(exp.t_keep) + j + (run_IC * 2)
                t_outs[idx] = sol.t[t_realised_pos]
                u_outs[idx] = Float64(round(total_pop[t_realised_pos]))
            else
                idx = length(exp.t_keep) + j + (run_IC * 2)
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
    save_at = _kw(kwargs, :save_at, exp.save_at)
    n_rep = _kw(kwargs, :n_rep, exp.n_rep)
    drug_treatment = _kw(kwargs, :drug_treatment, exp.drug_treatment)
    full_sol = _kw(kwargs, :full_sol, exp.full_sol)
    run_IC = _kw(kwargs, :run_IC, exp.run_IC)
    IC_n0 = _kw(kwargs, :IC_n0, exp.IC_n0)
    IC_tmax = _kw(kwargs, :IC_tmax, exp.IC_tmax)
    IC_treat_on = _kw(kwargs, :IC_treat_on, exp.IC_treat_on)
    run_colony = _kw(kwargs, :run_colony, exp.run_colony)
    nCol = _kw(kwargs, :nCol, exp.nCol)
    tCol = _kw(kwargs, :tCol, exp.tCol)
    ColNmax = _kw(kwargs, :ColNmax, exp.ColNmax)
    de = normalize_resdmg_drug_effect(_kw(kwargs, :drug_effect, params.drug_effect))

    @assert !(run_IC && run_colony) "Cannot run IC and colony assays at the same time."
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

    n_pop_obsv = n_rep * ((run_IC * 2) + length(exp.t_keep) + n_pass_eff) + (run_colony * 1)
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
                N_trans_switch = exp.N_trans_switch,
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
            N_trans_switch = exp.N_trans_switch,
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
            N_trans_switch = exp.N_trans_switch,
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

    if run_colony
        tx1_cols = 0
        tx0_cols = 0

        col_nR = max(Float64(round(params.rho * exp.n0)), 0.0)
        col_nS = max((exp.n0 - col_nR), 0.0)
        col_nDS = 0.0
        col_nDR = 0.0

        for j in 1:nCol
            total_cells = max(col_nS + col_nDS + col_nDR + col_nR, 1e-10)
            cell_probs = [col_nS, col_nDS, col_nDR, col_nR] ./ total_cells
            if any(cell_probs .< 0) || sum(cell_probs) == 0
                error("Invalid probability vector: $(cell_probs)")
            end
            sampled_phenos = zeros(Float64, 4)
            sampled_phenos[rand(Categorical(cell_probs))] = 1.0

            col_sim = SimParams(
                n0 = 1,
                t0 = 0.0,
                tmax = tCol,
                t_Pass = -1.0,
                Nmax = ColNmax + 10,
                Cc = exp.Cc,
                Nswitch = exp.Nswitch,
                N_trans_switch = exp.N_trans_switch,
                treat_ons = [1.0],
                treat_offs = [100.0],
                save_at = save_at,
                treat = true
            )
            col_state = ResDmgState(sampled_phenos[1], sampled_phenos[2], sampled_phenos[3], sampled_phenos[4];
                                    gam = 0.0, pass_num = 1)
            col_sol_tx1 = run_model_core_hybrid(model_eff, col_state, col_sim; treat = true)
            col_sol_tx0 = run_model_core_hybrid(model_eff, col_state, col_sim; treat = false)

            if round(sum(last(resdmg_pop_fun.(col_sol_tx1.u)))) >= ColNmax
                tx1_cols += 1
            end
            if round(sum(last(resdmg_pop_fun.(col_sol_tx0.u)))) >= ColNmax
                tx0_cols += 1
            end
        end

        col_prop = tx0_cols == 0 ? 0.0 : (tx1_cols / tx0_cols)
        append!(fin_t_outs, tCol)
        append!(fin_u_outs, col_prop)
    end

    for i in 1:n_rep
        rep_tmax = _replicate_tmax(exp.tmax, n_rep, i)

        if run_IC
            ic_nR = Float64(round(params.rho * IC_n0))
            ic_nS = IC_n0 - ic_nR
            ic_nDS = 0.0
            ic_nDR = 0.0
            ic_state = ResDmgState(ic_nS, ic_nDS, ic_nDR, ic_nR; gam = 0.0, pass_num = 1)
            ic_sim = SimParams(
                n0 = IC_n0,
                t0 = 0.0,
                tmax = IC_tmax,
                t_Pass = -1.0,
                Nmax = exp.Nmax,
                Cc = exp.Cc,
                Nswitch = exp.Nswitch,
                N_trans_switch = exp.N_trans_switch,
                treat_ons = [IC_treat_on],
                treat_offs = [100.0],
                save_at = save_at,
                treat = true
            )
            IC_sol_tx1 = run_model_core_hybrid(model_eff, ic_state, ic_sim; treat = true)
            IC_sol_tx0 = run_model_core_hybrid(model_eff, ic_state, ic_sim; treat = false)
        end

        sim = SimParams(
            n0 = nseed_last,
            t0 = 0.0,
            tmax = rep_tmax,
            t_Pass = exp.t_Pass,
            Nmax = exp.Nmax,
            Cc = exp.Cc,
            Nswitch = exp.Nswitch,
            N_trans_switch = exp.N_trans_switch,
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
        t_outs = Vector{Float64}(undef, (run_IC * 2) + length(exp.t_keep) + n_pass_eff)
        u_outs = Vector{Float64}(undef, (run_IC * 2) + length(exp.t_keep) + n_pass_eff)

        if run_IC
            t_outs[1] = last(IC_sol_tx1.t)
            u_outs[1] = Float64(round(sum(last(resdmg_pop_fun.(IC_sol_tx1.u)))))
            t_outs[2] = last(IC_sol_tx0.t)
            u_outs[2] = Float64(round(sum(last(resdmg_pop_fun.(IC_sol_tx0.u)))))
        end

        for j in 1:length(exp.t_keep)
            if !(exp.t_keep[j] in sol.t)
                t_closest_pos = findmin(abs.(sol.t .- exp.t_keep[j]))[2]
                t_outs[j + (run_IC * 2)] = sol.t[t_closest_pos]
                u_outs[j + (run_IC * 2)] = Float64(round(total_pop[t_closest_pos]))
            else
                t_realised_pos = findlast(sol.t .== exp.t_keep[j])
                t_outs[j + (run_IC * 2)] = sol.t[t_realised_pos]
                u_outs[j + (run_IC * 2)] = Float64(round(total_pop[t_realised_pos]))
            end
        end

        for j in 1:n_pass_eff
            if sum(resdmg_pass_fun.(sol.u) .== j) > 0
                t_realised_pos = findlast(resdmg_pass_fun.(sol.u) .== j)
                idx = length(exp.t_keep) + j + (run_IC * 2)
                t_outs[idx] = sol.t[t_realised_pos]
                u_outs[idx] = Float64(round(total_pop[t_realised_pos]))
            else
                idx = length(exp.t_keep) + j + (run_IC * 2)
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




