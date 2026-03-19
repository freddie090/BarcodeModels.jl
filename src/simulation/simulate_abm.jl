function _record_abm_outputs!(kmc_out, cells, rep::Int64, curr_P::Int64,
    cell_lin_df_vec::Vector{DataFrame},
    Nvec::Vector{Int64}, nS_vec::Vector{Int64}, nR_vec::Vector{Int64}, nE_vec::Vector{Int64},
    tvec::Vector{Float64}, Pvec::Vector{Int64};
    sub_sample_cells::Bool = false, K::Int64 = 0,
    samp_cell_lin_df_vec::Vector{DataFrame} = DataFrame[])

    live_cells = alive_cells(cells)
    bc_df = get_counts(live_cells, string("DT", rep, "_P", curr_P))
    push!(cell_lin_df_vec, bc_df)

    if sub_sample_cells
        if K <= length(live_cells)
            samp_cells = sample(live_cells, K, replace = false)
            push!(samp_cell_lin_df_vec, get_counts(samp_cells, string("DT", rep, "_P", curr_P)))
        else
            push!(samp_cell_lin_df_vec, get_counts(live_cells, string("DT", rep, "_P", curr_P)))
        end
    end

    update_track_vec!(kmc_out, Nvec, nS_vec, nR_vec, nE_vec, tvec, Pvec)
end

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

function _run_abm_passage_experiment!(
    model::ResPop_ABM,
    cells::Vector{CancerCell};
    t0::Float64,
    tmax::Float64,
    t_Pass::Union{Float64, Vector{Float64}},
    Nseed::Int64,
    Nmax::Int64,
    Cc::Int64,
    treat_ons::Vector{Float64},
    treat_offs::Vector{Float64},
    dt_save_at::Float64,
    Nbuff::Int64,
    R_real::String,
    t_frac::Float64,
    rep::Int64,
    treat::Bool,
    drug_effect::Symbol,
    sub_sample_cells::Bool = false,
    K::Int64 = 0
)
    progress_tol = 0.1

    curr_t = t0
    curr_P = 1
    tP_count = 1
    t_pass_vec = _passage_times(t_Pass, tmax)
    n_pass_eff = length(t_pass_vec) + 1

    function advance_passage_index!(time_now)
        while tP_count <= length(t_pass_vec) && t_pass_vec[tP_count] <= (time_now + progress_tol)
            tP_count += 1
        end
    end

    function compute_next_t()
        if tP_count <= length(t_pass_vec)
            return min(tmax, t_pass_vec[tP_count])
        else
            return tmax
        end
    end

    advance_passage_index!(curr_t)

    next_t = compute_next_t()
    cell_lin_df_vec = DataFrame[]
    samp_cell_lin_df_vec = DataFrame[]
    Nvec = Int64[]
    nS_vec = Int64[]
    nR_vec = Int64[]
    nE_vec = Int64[]
    tvec = Float64[]
    Pvec = Int64[]

    while curr_t <= (tmax + progress_tol)
        curr_t <= (next_t + progress_tol) || error("ABM non-progress invariant failed: curr_t ($(curr_t)) exceeds next_t ($(next_t)).")
        prev_curr_t = curr_t

        sim = ABMSimParams(
            t0 = curr_t,
            tmax = next_t,
            Nmax = Nmax,
            Cc = Cc,
            treat_ons = treat_ons,
            treat_offs = treat_offs,
            dt_save_at = dt_save_at,
            R_real = R_real,
            t_frac = t_frac,
            Passage = curr_P,
            drug_effect = drug_effect
        )
        state = ResPopABMState(cells)
        kmc_out = run_model_core_abm(model, state, sim; treat = treat)
        kmc_last_t = last(kmc_out.tvec)
        isfinite(kmc_last_t) || error("ABM produced a non-finite time value.")

        curr_t_candidate = min(max(kmc_last_t, curr_t), next_t)

        live_count = length(alive_cells(cells))
        if last(kmc_out.Nvec) >= Nmax
            if last(kmc_out.Pvec) >= n_pass_eff
                _record_abm_outputs!(kmc_out, cells, rep, curr_P, cell_lin_df_vec,
                                     Nvec, nS_vec, nR_vec, nE_vec, tvec, Pvec,
                                     sub_sample_cells = sub_sample_cells, K = K,
                                     samp_cell_lin_df_vec = samp_cell_lin_df_vec)
                curr_t = min(max(kmc_last_t, curr_t), tmax)
                break
            else
                _record_abm_outputs!(kmc_out, cells, rep, curr_P, cell_lin_df_vec,
                                     Nvec, nS_vec, nR_vec, nE_vec, tvec, Pvec,
                                     sub_sample_cells = sub_sample_cells, K = K,
                                     samp_cell_lin_df_vec = samp_cell_lin_df_vec)
                if live_count < Nseed
                    break
                else
                    live_cells = alive_cells(cells)
                    cells = sample(live_cells, Nseed, replace = false)
                    extend_with_dead_cells!(cells, Nbuff, make_dead_cell)
                    curr_P += 1
                    curr_t = curr_t_candidate
                    advance_passage_index!(curr_t)
                    next_t = compute_next_t()
                end
            end
        elseif live_count == 0
            _record_abm_outputs!(kmc_out, cells, rep, curr_P, cell_lin_df_vec,
                                 Nvec, nS_vec, nR_vec, nE_vec, tvec, Pvec,
                                 sub_sample_cells = sub_sample_cells, K = K,
                                 samp_cell_lin_df_vec = samp_cell_lin_df_vec)
            curr_t = min(max(kmc_last_t, curr_t), tmax)
            break
        elseif kmc_last_t >= (tmax - progress_tol)
            _record_abm_outputs!(kmc_out, cells, rep, curr_P, cell_lin_df_vec,
                                 Nvec, nS_vec, nR_vec, nE_vec, tvec, Pvec,
                                 sub_sample_cells = sub_sample_cells, K = K,
                                 samp_cell_lin_df_vec = samp_cell_lin_df_vec)
            curr_t = min(max(kmc_last_t, curr_t), tmax)
            break
        elseif tP_count <= length(t_pass_vec) && kmc_last_t >= (t_pass_vec[tP_count] - progress_tol)
            if last(kmc_out.Pvec) < n_pass_eff
                if live_count < Nseed
                    update_track_vec!(kmc_out, Nvec, nS_vec, nR_vec, nE_vec, tvec, Pvec)
                    curr_t = curr_t_candidate
                    advance_passage_index!(curr_t)
                    next_t = compute_next_t()
                else
                    _record_abm_outputs!(kmc_out, cells, rep, curr_P, cell_lin_df_vec,
                                         Nvec, nS_vec, nR_vec, nE_vec, tvec, Pvec,
                                         sub_sample_cells = sub_sample_cells, K = K,
                                         samp_cell_lin_df_vec = samp_cell_lin_df_vec)
                    live_cells = alive_cells(cells)
                    cells = sample(live_cells, Nseed, replace = false)
                    extend_with_dead_cells!(cells, Nbuff, make_dead_cell)
                    curr_P += 1
                    curr_t = curr_t_candidate
                    advance_passage_index!(curr_t)
                    next_t = compute_next_t()
                end
            end
        else
            error("ABM passage loop made no progress. prev_curr_t=$(prev_curr_t), curr_t=$(curr_t), next_t=$(next_t), kmc_last_t=$(kmc_last_t), tP_count=$(tP_count)")
        end

        if curr_t <= (prev_curr_t + progress_tol) && next_t <= (prev_curr_t + progress_tol)
            error("ABM passage loop stalled. prev_curr_t=$(prev_curr_t), curr_t=$(curr_t), next_t=$(next_t), kmc_last_t=$(kmc_last_t), tP_count=$(tP_count)")
        end
    end

    for i in 1:n_pass_eff
        if i != 1 && length(cell_lin_df_vec) < i
            temp_df = deepcopy(cell_lin_df_vec[i - 1])
            rename!(temp_df, [:bc, Symbol("DT", rep, "_P", i)])
            push!(cell_lin_df_vec, temp_df)
        end
    end

    out = Dict(
        "cell_lin_df_vec" => cell_lin_df_vec,
        "Nvec" => Nvec,
        "tvec" => tvec,
        "Pvec" => Pvec,
        "nS_vec" => nS_vec,
        "nR_vec" => nR_vec,
        "nE_vec" => nE_vec
    )
    if sub_sample_cells
        out["sub_samp_cell_lin_df_vec"] = samp_cell_lin_df_vec
    end
    return out
end

function _expand_split_cells_abm(model::ResPop_ABM, exp::ExperimentParams, n_rep::Int64;
    R_real::String = "b",
    drug_effect::Symbol = model.params.drug_effect,
    skew_lib::Bool = model.abm.skew_lib,
    bc_unif::Float64 = model.abm.bc_unif,
    Nbc::Int64 = model.abm.Nbc,
    dt_save_at::Float64 = model.abm.dt_save_at,
    t_frac::Float64 = model.abm.t_frac)

    Nbuff = model.abm.Nbuff
    exp_cells = seed_cells(exp.n0, model.params.rho, Nbuff;
                           skew_lib = skew_lib, bc_unif = bc_unif, Nbc = Nbc)

    expansion_model = ResPop_ABM(_copy_respop_params(model.params; al = 0.0, drug_effect = drug_effect);
                                 abm = model.abm)

    if (exp.t_exp isa Vector{Float64}) && (exp.Nseed isa Vector{Int64})
        @assert length(exp.t_exp) == length(exp.Nseed) "t_exp and Nseed vectors must be of same length"

        for i in 1:(length(exp.t_exp) - 1)
            stage_sim = ABMSimParams(
                t0 = 0.0,
                tmax = exp.t_exp[i],
                Nmax = exp.Nmax,
                Cc = exp.Cc,
                treat_ons = [0.0],
                treat_offs = [0.0],
                dt_save_at = dt_save_at,
                R_real = R_real,
                t_frac = t_frac,
                Passage = 1,
                drug_effect = drug_effect
            )
            run_model_core_abm(expansion_model, ResPopABMState(exp_cells), stage_sim; treat = false)

            exp_cells = alive_cells(exp_cells)
            if length(exp_cells) < exp.Nseed[i]
                error("Not enough cells after expansion at stage $i for bottlenecking.")
            end
            exp_cells = sample(exp_cells, exp.Nseed[i], replace = false)
            extend_with_dead_cells!(exp_cells, Nbuff, make_dead_cell)
        end

        final_sim = ABMSimParams(
            t0 = 0.0,
            tmax = exp.t_exp[end],
            Nmax = exp.Nmax,
            Cc = exp.Cc,
            treat_ons = [0.0],
            treat_offs = [0.0],
            dt_save_at = dt_save_at,
            R_real = R_real,
            t_frac = t_frac,
            Passage = 1,
            drug_effect = drug_effect
        )
        run_model_core_abm(expansion_model, ResPopABMState(exp_cells), final_sim; treat = false)
        exp_cells = alive_cells(exp_cells)
        final_seed = exp.Nseed[end]
    elseif (exp.t_exp isa Float64) && (exp.Nseed isa Int64)
        final_sim = ABMSimParams(
            t0 = 0.0,
            tmax = exp.t_exp,
            Nmax = exp.Nmax,
            Cc = exp.Cc,
            treat_ons = [0.0],
            treat_offs = [0.0],
            dt_save_at = dt_save_at,
            R_real = R_real,
            t_frac = t_frac,
            Passage = 1,
            drug_effect = drug_effect
        )
        run_model_core_abm(expansion_model, ResPopABMState(exp_cells), final_sim; treat = false)
        exp_cells = alive_cells(exp_cells)
        final_seed = exp.Nseed
    else
        error("t_exp and Nseed must both be scalars or both be vectors of equal length.")
    end

    n_rep * final_seed <= length(exp_cells) || error("Not enough cells for $n_rep replicates of size $final_seed.")
    rep_cells = sample(exp_cells, n_rep * final_seed, replace = false)
    rep_cells = reshape(rep_cells, (final_seed, n_rep))

    fin_rep_cells = Vector{Vector{CancerCell}}(undef, n_rep)
    for i in 1:n_rep
        fin_rep_cells[i] = collect(rep_cells[:, i])
        extend_with_dead_cells!(fin_rep_cells[i], Nbuff, make_dead_cell)
    end
    return fin_rep_cells
end

function _simulate_experiment_abm(model::ResPop_ABM, exp::ExperimentParams; kwargs...)
    n_rep = _kw(kwargs, :n_rep, exp.n_rep)
    R_real = _kw(kwargs, :R_real, "b")
    t_frac = _kw(kwargs, :t_frac, model.abm.t_frac)
    just_lin = _kw(kwargs, :just_lin, false)
    de = normalize_respop_drug_effect(_kw(kwargs, :drug_effect, model.params.drug_effect))
    drug_treatment = _kw(kwargs, :drug_treatment, exp.drug_treatment)
    sub_sample_cells = _kw(kwargs, :sub_sample_cells, model.abm.sub_sample_cells)
    K = _kw(kwargs, :K, model.abm.K)
    skew_lib = _kw(kwargs, :skew_lib, model.abm.skew_lib)
    bc_unif = _kw(kwargs, :bc_unif, model.abm.bc_unif)
    Nbc = _kw(kwargs, :Nbc, model.abm.Nbc)
    run_IC = _kw(kwargs, :run_IC, exp.run_IC)
    IC_n0 = _kw(kwargs, :IC_n0, exp.IC_n0)
    IC_tmax = _kw(kwargs, :IC_tmax, exp.IC_tmax)
    IC_treat_on = _kw(kwargs, :IC_treat_on, exp.IC_treat_on)
    run_colony = _kw(kwargs, :run_colony, exp.run_colony)
    nCol = _kw(kwargs, :nCol, exp.nCol)
    tCol = _kw(kwargs, :tCol, exp.tCol)
    ColNmax = _kw(kwargs, :ColNmax, exp.ColNmax)
    dt_save_at = _kw(kwargs, :dt_save_at, model.abm.dt_save_at)

    @assert !(run_IC && run_colony) "Cannot run IC and colony assays at the same time."

    n_pass_eff = length(_passage_times(exp.t_Pass, exp.tmax)) + 1

    model_eff = _with_drug_effect(model, de)
    rep_cells = _expand_split_cells_abm(model_eff, exp, n_rep;
                                        R_real = R_real,
                                        drug_effect = de,
                                        skew_lib = skew_lib,
                                        bc_unif = bc_unif,
                                        Nbc = Nbc,
                                        dt_save_at = dt_save_at,
                                        t_frac = t_frac)

    fin_t_outs = Float64[]
    fin_u_outs = Float64[]
    lin_df_outs = DataFrame[]
    sub_lin_df_outs = DataFrame[]
    sim_dfs = DataFrame[]

    if !just_lin && run_colony
        col_cells_tx1 = seed_cells(nCol, model.params.rho, Int64(1e6))
        col_cells_tx0 = seed_cells(nCol, model.params.rho, Int64(1e6))
        col_sim = ABMSimParams(
            t0 = 0.0,
            tmax = tCol,
            Nmax = exp.Nmax,
            Cc = exp.Cc,
            treat_ons = [1.0],
            treat_offs = [1000.0],
            dt_save_at = dt_save_at,
            R_real = R_real,
            t_frac = t_frac,
            Passage = 1,
            drug_effect = de
        )

        run_model_core_abm(model_eff, ResPopABMState(col_cells_tx1), col_sim; treat = true)
        run_model_core_abm(model_eff, ResPopABMState(col_cells_tx0), col_sim; treat = false)

        col_tx1_bcs = get_counts(alive_cells(col_cells_tx1), "col_tx1")
        col_tx0_bcs = get_counts(alive_cells(col_cells_tx0), "col_tx0")
        tx1_cols = sum(col_tx1_bcs[!, :col_tx1] .> ColNmax)
        tx0_cols = sum(col_tx0_bcs[!, :col_tx0] .> ColNmax)
        col_prop = tx0_cols == 0 ? 0.0 : tx1_cols / tx0_cols

        append!(fin_t_outs, tCol)
        append!(fin_u_outs, col_prop)
    end

    nseed_last = _nseed_last(exp.Nseed)
    for i in 1:n_rep
        if !just_lin && run_IC
            IC_cells_1 = seed_cells(IC_n0, model.params.rho, model.abm.Nbuff)
            IC_cells_0 = seed_cells(IC_n0, model.params.rho, model.abm.Nbuff)

            IC_sim_tx1 = _run_abm_passage_experiment!(
                model_eff, IC_cells_1;
                t0 = 0.0, tmax = IC_tmax, t_Pass = 1000.0,
                Nseed = IC_n0, Nmax = exp.Nmax, Cc = exp.Cc,
                treat_ons = [IC_treat_on], treat_offs = [100.0],
                dt_save_at = dt_save_at, Nbuff = model.abm.Nbuff,
                R_real = R_real, t_frac = t_frac, rep = i,
                treat = true, drug_effect = de
            )

            IC_sim_tx0 = _run_abm_passage_experiment!(
                model_eff, IC_cells_0;
                t0 = 0.0, tmax = IC_tmax, t_Pass = 1000.0,
                Nseed = IC_n0, Nmax = exp.Nmax, Cc = exp.Cc,
                treat_ons = [IC_treat_on], treat_offs = [100.0],
                dt_save_at = dt_save_at, Nbuff = model.abm.Nbuff,
                R_real = R_real, t_frac = t_frac, rep = i,
                treat = false, drug_effect = de
            )
        end

        extend_with_dead_cells!(rep_cells[i], model.abm.Nbuff, make_dead_cell)

        sim = _run_abm_passage_experiment!(
            model_eff, rep_cells[i];
            t0 = 0.0, tmax = exp.tmax, t_Pass = exp.t_Pass,
            Nseed = nseed_last, Nmax = exp.Nmax, Cc = exp.Cc,
            treat_ons = exp.treat_ons, treat_offs = exp.treat_offs,
            dt_save_at = dt_save_at, Nbuff = model.abm.Nbuff,
            R_real = R_real, t_frac = t_frac, rep = i,
            treat = drug_treatment, drug_effect = de,
            sub_sample_cells = sub_sample_cells, K = K
        )

        sim_df = DataFrame(
            t = sim["tvec"],
            N = sim["Nvec"],
            nS = sim["nS_vec"],
            nR = sim["nR_vec"],
            nE = sim["nE_vec"],
            rep = i
        )
        push!(sim_dfs, sim_df)

        push!(lin_df_outs, join_dfs(sim["cell_lin_df_vec"], "bc"))
        if sub_sample_cells
            push!(sub_lin_df_outs, join_dfs(sim["sub_samp_cell_lin_df_vec"], "bc"))
        end

        if !just_lin
            t_outs = Vector{Float64}(undef, (run_IC * 2) + length(exp.t_keep) + n_pass_eff)
            u_outs = Vector{Float64}(undef, (run_IC * 2) + length(exp.t_keep) + n_pass_eff)

            if run_IC
                t_outs[1] = last(IC_sim_tx1["tvec"])
                u_outs[1] = Float64(round(last(IC_sim_tx1["Nvec"])))
                t_outs[2] = last(IC_sim_tx0["tvec"])
                u_outs[2] = Float64(round(last(IC_sim_tx0["Nvec"])))
            end

            for j in 1:length(exp.t_keep)
                idx = j + (run_IC * 2)
                if !(exp.t_keep[j] in sim["tvec"])
                    t_closest_pos = findmin(abs.(sim["tvec"] .- exp.t_keep[j]))[2]
                    t_outs[idx] = sim["tvec"][t_closest_pos]
                    u_outs[idx] = sim["Nvec"][t_closest_pos]
                else
                    t_realised_pos = findlast(sim["tvec"] .== exp.t_keep[j])
                    t_outs[idx] = sim["tvec"][t_realised_pos]
                    u_outs[idx] = sim["Nvec"][t_realised_pos]
                end
            end

            for j in 1:n_pass_eff
                idx = length(exp.t_keep) + j + (run_IC * 2)
                if sum(sim["Pvec"] .== j) > 0
                    t_realised_pos = findlast(sim["Pvec"] .== j)
                    t_outs[idx] = sim["tvec"][t_realised_pos]
                    u_outs[idx] = sim["Nvec"][t_realised_pos]
                else
                    t_outs[idx] = t_outs[idx - 1]
                    u_outs[idx] = u_outs[idx - 1]
                end
            end

            append!(fin_t_outs, t_outs)
            append!(fin_u_outs, u_outs)
        end
    end

    fin_t_outs = round.(fin_t_outs; digits = 0)
    fin_lin_df = join_dfs(lin_df_outs, "bc")
    fin_sol_df = isempty(sim_dfs) ? DataFrame() : vcat(sim_dfs...)

    out = Dict{String, Any}(
        "lin_df" => fin_lin_df,
        "sol_df" => fin_sol_df
    )
    if !just_lin
        out["t"] = fin_t_outs
        out["u"] = fin_u_outs
    end
    if sub_sample_cells
        out["sub_lin_df"] = join_dfs(sub_lin_df_outs, "bc")
    end
    return out
end

function _record_resdmg_abm_outputs!(kmc_out, cells, rep::Int64, curr_P::Int64,
    cell_lin_df_vec::Vector{DataFrame},
    Nvec::Vector{Int64}, nS_vec::Vector{Int64}, nD_vec::Vector{Int64}, nR_vec::Vector{Int64}, nE_vec::Vector{Int64},
    tvec::Vector{Float64}, Pvec::Vector{Int64};
    sub_sample_cells::Bool = false, K::Int64 = 0,
    samp_cell_lin_df_vec::Vector{DataFrame} = DataFrame[])

    live_cells = alive_cells(cells)
    bc_df = get_counts(live_cells, string("DT", rep, "_P", curr_P))
    push!(cell_lin_df_vec, bc_df)

    if sub_sample_cells
        if K <= length(live_cells)
            samp_cells = sample(live_cells, K, replace = false)
            push!(samp_cell_lin_df_vec, get_counts(samp_cells, string("DT", rep, "_P", curr_P)))
        else
            push!(samp_cell_lin_df_vec, get_counts(live_cells, string("DT", rep, "_P", curr_P)))
        end
    end

    update_track_vec_resdmg!(kmc_out, Nvec, nS_vec, nD_vec, nR_vec, nE_vec, tvec, Pvec)
end

function _run_abm_passage_experiment!(
    model::ResDmg_ABM,
    cells::Vector{ResDmgCell};
    t0::Float64,
    tmax::Float64,
    t_Pass::Union{Float64, Vector{Float64}},
    Nseed::Int64,
    Nmax::Int64,
    Cc::Int64,
    treat_ons::Vector{Float64},
    treat_offs::Vector{Float64},
    dt_save_at::Float64,
    Nbuff::Int64,
    R_real::String,
    t_frac::Float64,
    rep::Int64,
    treat::Bool,
    drug_effect::Symbol,
    sub_sample_cells::Bool = false,
    K::Int64 = 0
)
    progress_tol = 0.1

    curr_t = t0
    curr_P = 1
    tP_count = 1
    t_pass_vec = _passage_times(t_Pass, tmax)
    n_pass_eff = length(t_pass_vec) + 1

    function advance_passage_index!(time_now)
        while tP_count <= length(t_pass_vec) && t_pass_vec[tP_count] <= (time_now + progress_tol)
            tP_count += 1
        end
    end

    function compute_next_t()
        if tP_count <= length(t_pass_vec)
            return min(tmax, t_pass_vec[tP_count])
        else
            return tmax
        end
    end

    advance_passage_index!(curr_t)

    next_t = compute_next_t()
    cell_lin_df_vec = DataFrame[]
    samp_cell_lin_df_vec = DataFrame[]
    Nvec = Int64[]
    nS_vec = Int64[]
    nD_vec = Int64[]
    nR_vec = Int64[]
    nE_vec = Int64[]
    tvec = Float64[]
    Pvec = Int64[]

    while curr_t <= (tmax + progress_tol)
        curr_t <= (next_t + progress_tol) || error("ABM non-progress invariant failed: curr_t ($(curr_t)) exceeds next_t ($(next_t)).")
        prev_curr_t = curr_t

        sim = ABMSimParams(
            t0 = curr_t,
            tmax = next_t,
            Nmax = Nmax,
            Cc = Cc,
            treat_ons = treat_ons,
            treat_offs = treat_offs,
            dt_save_at = dt_save_at,
            R_real = R_real,
            t_frac = t_frac,
            Passage = curr_P,
            drug_effect = drug_effect
        )
        state = ResDmgABMState(cells)
        kmc_out = run_model_core_abm(model, state, sim; treat = treat)
        kmc_last_t = last(kmc_out.tvec)
        isfinite(kmc_last_t) || error("ABM produced a non-finite time value.")

        curr_t_candidate = min(max(kmc_last_t, curr_t), next_t)

        live_count = length(alive_cells(cells))
        if last(kmc_out.Nvec) >= Nmax
            if last(kmc_out.Pvec) >= n_pass_eff
                _record_resdmg_abm_outputs!(kmc_out, cells, rep, curr_P, cell_lin_df_vec,
                                            Nvec, nS_vec, nD_vec, nR_vec, nE_vec, tvec, Pvec,
                                            sub_sample_cells = sub_sample_cells, K = K,
                                            samp_cell_lin_df_vec = samp_cell_lin_df_vec)
                curr_t = min(max(kmc_last_t, curr_t), tmax)
                break
            else
                _record_resdmg_abm_outputs!(kmc_out, cells, rep, curr_P, cell_lin_df_vec,
                                            Nvec, nS_vec, nD_vec, nR_vec, nE_vec, tvec, Pvec,
                                            sub_sample_cells = sub_sample_cells, K = K,
                                            samp_cell_lin_df_vec = samp_cell_lin_df_vec)
                if live_count < Nseed
                    break
                else
                    live_cells = alive_cells(cells)
                    cells = sample(live_cells, Nseed, replace = false)
                    extend_with_dead_cells!(cells, Nbuff, make_dead_resdmg_cell)
                    curr_P += 1
                    curr_t = curr_t_candidate
                    advance_passage_index!(curr_t)
                    next_t = compute_next_t()
                end
            end
        elseif live_count == 0
            _record_resdmg_abm_outputs!(kmc_out, cells, rep, curr_P, cell_lin_df_vec,
                                        Nvec, nS_vec, nD_vec, nR_vec, nE_vec, tvec, Pvec,
                                        sub_sample_cells = sub_sample_cells, K = K,
                                        samp_cell_lin_df_vec = samp_cell_lin_df_vec)
            curr_t = min(max(kmc_last_t, curr_t), tmax)
            break
        elseif kmc_last_t >= (tmax - progress_tol)
            _record_resdmg_abm_outputs!(kmc_out, cells, rep, curr_P, cell_lin_df_vec,
                                        Nvec, nS_vec, nD_vec, nR_vec, nE_vec, tvec, Pvec,
                                        sub_sample_cells = sub_sample_cells, K = K,
                                        samp_cell_lin_df_vec = samp_cell_lin_df_vec)
            curr_t = min(max(kmc_last_t, curr_t), tmax)
            break
        elseif tP_count <= length(t_pass_vec) && kmc_last_t >= (t_pass_vec[tP_count] - progress_tol)
            if last(kmc_out.Pvec) < n_pass_eff
                if live_count < Nseed
                    update_track_vec_resdmg!(kmc_out, Nvec, nS_vec, nD_vec, nR_vec, nE_vec, tvec, Pvec)
                    curr_t = curr_t_candidate
                    advance_passage_index!(curr_t)
                    next_t = compute_next_t()
                else
                    _record_resdmg_abm_outputs!(kmc_out, cells, rep, curr_P, cell_lin_df_vec,
                                                Nvec, nS_vec, nD_vec, nR_vec, nE_vec, tvec, Pvec,
                                                sub_sample_cells = sub_sample_cells, K = K,
                                                samp_cell_lin_df_vec = samp_cell_lin_df_vec)
                    live_cells = alive_cells(cells)
                    cells = sample(live_cells, Nseed, replace = false)
                    extend_with_dead_cells!(cells, Nbuff, make_dead_resdmg_cell)
                    curr_P += 1
                    curr_t = curr_t_candidate
                    advance_passage_index!(curr_t)
                    next_t = compute_next_t()
                end
            end
        else
            error("ABM passage loop made no progress. prev_curr_t=$(prev_curr_t), curr_t=$(curr_t), next_t=$(next_t), kmc_last_t=$(kmc_last_t), tP_count=$(tP_count)")
        end

        if curr_t <= (prev_curr_t + progress_tol) && next_t <= (prev_curr_t + progress_tol)
            error("ABM passage loop stalled. prev_curr_t=$(prev_curr_t), curr_t=$(curr_t), next_t=$(next_t), kmc_last_t=$(kmc_last_t), tP_count=$(tP_count)")
        end
    end

    for i in 1:n_pass_eff
        if i != 1 && length(cell_lin_df_vec) < i
            temp_df = deepcopy(cell_lin_df_vec[i - 1])
            rename!(temp_df, [:bc, Symbol("DT", rep, "_P", i)])
            push!(cell_lin_df_vec, temp_df)
        end
    end

    out = Dict(
        "cell_lin_df_vec" => cell_lin_df_vec,
        "Nvec" => Nvec,
        "tvec" => tvec,
        "Pvec" => Pvec,
        "nS_vec" => nS_vec,
        "nD_vec" => nD_vec,
        "nR_vec" => nR_vec,
        "nE_vec" => nE_vec
    )
    if sub_sample_cells
        out["sub_samp_cell_lin_df_vec"] = samp_cell_lin_df_vec
    end
    return out
end

function _expand_split_cells_abm(model::ResDmg_ABM, exp::ExperimentParams, n_rep::Int64;
    R_real::String = "b",
    drug_effect::Symbol = model.params.drug_effect,
    skew_lib::Bool = model.abm.skew_lib,
    bc_unif::Float64 = model.abm.bc_unif,
    Nbc::Int64 = model.abm.Nbc,
    dt_save_at::Float64 = model.abm.dt_save_at,
    t_frac::Float64 = model.abm.t_frac)

    Nbuff = model.abm.Nbuff
    exp_cells = seed_resdmg_cells(exp.n0, model.params.rho, Nbuff;
                                  skew_lib = skew_lib, bc_unif = bc_unif, Nbc = Nbc)

    expansion_model = ResDmg_ABM(_copy_resdmg_params(model.params; al = 0.0, drug_effect = drug_effect);
                                 abm = model.abm)

    if (exp.t_exp isa Vector{Float64}) && (exp.Nseed isa Vector{Int64})
        @assert length(exp.t_exp) == length(exp.Nseed) "t_exp and Nseed vectors must be of same length"

        for i in 1:(length(exp.t_exp) - 1)
            stage_sim = ABMSimParams(
                t0 = 0.0,
                tmax = exp.t_exp[i],
                Nmax = exp.Nmax,
                Cc = exp.Cc,
                treat_ons = [0.0],
                treat_offs = [0.0],
                dt_save_at = dt_save_at,
                R_real = R_real,
                t_frac = t_frac,
                Passage = 1,
                drug_effect = drug_effect
            )
            run_model_core_abm(expansion_model, ResDmgABMState(exp_cells), stage_sim; treat = false)

            exp_cells = alive_cells(exp_cells)
            if length(exp_cells) < exp.Nseed[i]
                error("Not enough cells after expansion at stage $i for bottlenecking.")
            end
            exp_cells = sample(exp_cells, exp.Nseed[i], replace = false)
            extend_with_dead_cells!(exp_cells, Nbuff, make_dead_resdmg_cell)
        end

        final_sim = ABMSimParams(
            t0 = 0.0,
            tmax = exp.t_exp[end],
            Nmax = exp.Nmax,
            Cc = exp.Cc,
            treat_ons = [0.0],
            treat_offs = [0.0],
            dt_save_at = dt_save_at,
            R_real = R_real,
            t_frac = t_frac,
            Passage = 1,
            drug_effect = drug_effect
        )
        run_model_core_abm(expansion_model, ResDmgABMState(exp_cells), final_sim; treat = false)
        exp_cells = alive_cells(exp_cells)
        final_seed = exp.Nseed[end]
    elseif (exp.t_exp isa Float64) && (exp.Nseed isa Int64)
        final_sim = ABMSimParams(
            t0 = 0.0,
            tmax = exp.t_exp,
            Nmax = exp.Nmax,
            Cc = exp.Cc,
            treat_ons = [0.0],
            treat_offs = [0.0],
            dt_save_at = dt_save_at,
            R_real = R_real,
            t_frac = t_frac,
            Passage = 1,
            drug_effect = drug_effect
        )
        run_model_core_abm(expansion_model, ResDmgABMState(exp_cells), final_sim; treat = false)
        exp_cells = alive_cells(exp_cells)
        final_seed = exp.Nseed
    else
        error("t_exp and Nseed must both be scalars or both be vectors of equal length.")
    end

    n_rep * final_seed <= length(exp_cells) || error("Not enough cells for $n_rep replicates of size $final_seed.")
    rep_cells = sample(exp_cells, n_rep * final_seed, replace = false)
    rep_cells = reshape(rep_cells, (final_seed, n_rep))

    fin_rep_cells = Vector{Vector{ResDmgCell}}(undef, n_rep)
    for i in 1:n_rep
        fin_rep_cells[i] = collect(rep_cells[:, i])
        extend_with_dead_cells!(fin_rep_cells[i], Nbuff, make_dead_resdmg_cell)
    end
    return fin_rep_cells
end

function _simulate_experiment_abm(model::ResDmg_ABM, exp::ExperimentParams; kwargs...)
    n_rep = _kw(kwargs, :n_rep, exp.n_rep)
    R_real = _kw(kwargs, :R_real, "b")
    t_frac = _kw(kwargs, :t_frac, model.abm.t_frac)
    just_lin = _kw(kwargs, :just_lin, false)
    de = normalize_resdmg_drug_effect(_kw(kwargs, :drug_effect, model.params.drug_effect))
    drug_treatment = _kw(kwargs, :drug_treatment, exp.drug_treatment)
    sub_sample_cells = _kw(kwargs, :sub_sample_cells, model.abm.sub_sample_cells)
    K = _kw(kwargs, :K, model.abm.K)
    skew_lib = _kw(kwargs, :skew_lib, model.abm.skew_lib)
    bc_unif = _kw(kwargs, :bc_unif, model.abm.bc_unif)
    Nbc = _kw(kwargs, :Nbc, model.abm.Nbc)
    run_IC = _kw(kwargs, :run_IC, exp.run_IC)
    IC_n0 = _kw(kwargs, :IC_n0, exp.IC_n0)
    IC_tmax = _kw(kwargs, :IC_tmax, exp.IC_tmax)
    IC_treat_on = _kw(kwargs, :IC_treat_on, exp.IC_treat_on)
    run_colony = _kw(kwargs, :run_colony, exp.run_colony)
    nCol = _kw(kwargs, :nCol, exp.nCol)
    tCol = _kw(kwargs, :tCol, exp.tCol)
    ColNmax = _kw(kwargs, :ColNmax, exp.ColNmax)
    dt_save_at = _kw(kwargs, :dt_save_at, model.abm.dt_save_at)

    @assert !(run_IC && run_colony) "Cannot run IC and colony assays at the same time."

    n_pass_eff = length(_passage_times(exp.t_Pass, exp.tmax)) + 1

    model_eff = _with_drug_effect(model, de)
    rep_cells = _expand_split_cells_abm(model_eff, exp, n_rep;
                                        R_real = R_real,
                                        drug_effect = de,
                                        skew_lib = skew_lib,
                                        bc_unif = bc_unif,
                                        Nbc = Nbc,
                                        dt_save_at = dt_save_at,
                                        t_frac = t_frac)

    fin_t_outs = Float64[]
    fin_u_outs = Float64[]
    lin_df_outs = DataFrame[]
    sub_lin_df_outs = DataFrame[]
    sim_dfs = DataFrame[]

    if !just_lin && run_colony
        col_cells_tx1 = seed_resdmg_cells(nCol, model.params.rho, Int64(1e6))
        col_cells_tx0 = seed_resdmg_cells(nCol, model.params.rho, Int64(1e6))
        col_sim = ABMSimParams(
            t0 = 0.0,
            tmax = tCol,
            Nmax = exp.Nmax,
            Cc = exp.Cc,
            treat_ons = [1.0],
            treat_offs = [1000.0],
            dt_save_at = dt_save_at,
            R_real = R_real,
            t_frac = t_frac,
            Passage = 1,
            drug_effect = de
        )

        run_model_core_abm(model_eff, ResDmgABMState(col_cells_tx1), col_sim; treat = true)
        run_model_core_abm(model_eff, ResDmgABMState(col_cells_tx0), col_sim; treat = false)

        col_tx1_bcs = get_counts(alive_cells(col_cells_tx1), "col_tx1")
        col_tx0_bcs = get_counts(alive_cells(col_cells_tx0), "col_tx0")
        tx1_cols = sum(col_tx1_bcs[!, :col_tx1] .> ColNmax)
        tx0_cols = sum(col_tx0_bcs[!, :col_tx0] .> ColNmax)
        col_prop = tx0_cols == 0 ? 0.0 : tx1_cols / tx0_cols

        append!(fin_t_outs, tCol)
        append!(fin_u_outs, col_prop)
    end

    nseed_last = _nseed_last(exp.Nseed)
    for i in 1:n_rep
        if !just_lin && run_IC
            IC_cells_1 = seed_resdmg_cells(IC_n0, model.params.rho, model.abm.Nbuff)
            IC_cells_0 = seed_resdmg_cells(IC_n0, model.params.rho, model.abm.Nbuff)

            IC_sim_tx1 = _run_abm_passage_experiment!(
                model_eff, IC_cells_1;
                t0 = 0.0, tmax = IC_tmax, t_Pass = 1000.0,
                Nseed = IC_n0, Nmax = exp.Nmax, Cc = exp.Cc,
                treat_ons = [IC_treat_on], treat_offs = [100.0],
                dt_save_at = dt_save_at, Nbuff = model.abm.Nbuff,
                R_real = R_real, t_frac = t_frac, rep = i,
                treat = true, drug_effect = de
            )

            IC_sim_tx0 = _run_abm_passage_experiment!(
                model_eff, IC_cells_0;
                t0 = 0.0, tmax = IC_tmax, t_Pass = 1000.0,
                Nseed = IC_n0, Nmax = exp.Nmax, Cc = exp.Cc,
                treat_ons = [IC_treat_on], treat_offs = [100.0],
                dt_save_at = dt_save_at, Nbuff = model.abm.Nbuff,
                R_real = R_real, t_frac = t_frac, rep = i,
                treat = false, drug_effect = de
            )
        end

        extend_with_dead_cells!(rep_cells[i], model.abm.Nbuff, make_dead_resdmg_cell)

        sim = _run_abm_passage_experiment!(
            model_eff, rep_cells[i];
            t0 = 0.0, tmax = exp.tmax, t_Pass = exp.t_Pass,
            Nseed = nseed_last, Nmax = exp.Nmax, Cc = exp.Cc,
            treat_ons = exp.treat_ons, treat_offs = exp.treat_offs,
            dt_save_at = dt_save_at, Nbuff = model.abm.Nbuff,
            R_real = R_real, t_frac = t_frac, rep = i,
            treat = drug_treatment, drug_effect = de,
            sub_sample_cells = sub_sample_cells, K = K
        )

        sim_df = DataFrame(
            t = sim["tvec"],
            N = sim["Nvec"],
            nS = sim["nS_vec"],
            nD = sim["nD_vec"],
            nR = sim["nR_vec"],
            nE = sim["nE_vec"],
            rep = i
        )
        push!(sim_dfs, sim_df)

        push!(lin_df_outs, join_dfs(sim["cell_lin_df_vec"], "bc"))
        if sub_sample_cells
            push!(sub_lin_df_outs, join_dfs(sim["sub_samp_cell_lin_df_vec"], "bc"))
        end

        if !just_lin
            t_outs = Vector{Float64}(undef, (run_IC * 2) + length(exp.t_keep) + n_pass_eff)
            u_outs = Vector{Float64}(undef, (run_IC * 2) + length(exp.t_keep) + n_pass_eff)

            if run_IC
                t_outs[1] = last(IC_sim_tx1["tvec"])
                u_outs[1] = Float64(round(last(IC_sim_tx1["Nvec"])))
                t_outs[2] = last(IC_sim_tx0["tvec"])
                u_outs[2] = Float64(round(last(IC_sim_tx0["Nvec"])))
            end

            for j in 1:length(exp.t_keep)
                idx = j + (run_IC * 2)
                if !(exp.t_keep[j] in sim["tvec"])
                    t_closest_pos = findmin(abs.(sim["tvec"] .- exp.t_keep[j]))[2]
                    t_outs[idx] = sim["tvec"][t_closest_pos]
                    u_outs[idx] = sim["Nvec"][t_closest_pos]
                else
                    t_realised_pos = findlast(sim["tvec"] .== exp.t_keep[j])
                    t_outs[idx] = sim["tvec"][t_realised_pos]
                    u_outs[idx] = sim["Nvec"][t_realised_pos]
                end
            end

            for j in 1:n_pass_eff
                idx = length(exp.t_keep) + j + (run_IC * 2)
                if sum(sim["Pvec"] .== j) > 0
                    t_realised_pos = findlast(sim["Pvec"] .== j)
                    t_outs[idx] = sim["tvec"][t_realised_pos]
                    u_outs[idx] = sim["Nvec"][t_realised_pos]
                else
                    t_outs[idx] = t_outs[idx - 1]
                    u_outs[idx] = u_outs[idx - 1]
                end
            end

            append!(fin_t_outs, t_outs)
            append!(fin_u_outs, u_outs)
        end
    end

    fin_t_outs = round.(fin_t_outs; digits = 0)
    fin_lin_df = join_dfs(lin_df_outs, "bc")
    fin_sol_df = isempty(sim_dfs) ? DataFrame() : vcat(sim_dfs...)

    out = Dict{String, Any}(
        "lin_df" => fin_lin_df,
        "sol_df" => fin_sol_df
    )
    if !just_lin
        out["t"] = fin_t_outs
        out["u"] = fin_u_outs
    end
    if sub_sample_cells
        out["sub_lin_df"] = join_dfs(sub_lin_df_outs, "bc")
    end
    return out
end

