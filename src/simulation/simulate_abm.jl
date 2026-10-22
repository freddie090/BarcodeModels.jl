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
    K::Int64 = 0,
    neutral_tracker::Union{Nothing, NeutralMutationTracker} = nothing,
    instant_kill::Union{Nothing, InstantKillState} = nothing
)
    progress_tol = instant_kill === nothing ? 0.1 : 0.0

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

    function sync_passage_schedule!()
        tP_count = max(tP_count, min(curr_P, length(t_pass_vec) + 1))
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
        kmc_out = run_model_core_abm(model, state, sim; treat = treat, neutral_tracker = neutral_tracker, instant_kill = instant_kill)
        kmc_last_t = last(kmc_out.tvec)
        isfinite(kmc_last_t) || error("ABM produced a non-finite time value.")

        curr_t_candidate = min(max(kmc_last_t, curr_t), next_t)

        live_count = length(alive_cells(cells))
        if last(kmc_out.Nvec) >= Nmax
            if last(kmc_out.Pvec) >= n_pass_eff
                _record_abm_outputs!(kmc_out, cells, rep, curr_P, cell_lin_df_vec,
                                     Nvec, nS_vec, nR_vec, nE_vec, tvec, Pvec,
                                     sub_sample_cells = sub_sample_cells, K = K,
                                     samp_cell_lin_df_vec = samp_cell_lin_df_vec,
                                     neutral_tracker = neutral_tracker)
                curr_t = min(max(kmc_last_t, curr_t), tmax)
                break
            else
                _record_abm_outputs!(kmc_out, cells, rep, curr_P, cell_lin_df_vec,
                                     Nvec, nS_vec, nR_vec, nE_vec, tvec, Pvec,
                                     sub_sample_cells = sub_sample_cells, K = K,
                                     samp_cell_lin_df_vec = samp_cell_lin_df_vec,
                                     neutral_tracker = neutral_tracker)
                if live_count < Nseed
                    break
                else
                    live_cells = alive_cells(cells)
                    cells = sample(live_cells, Nseed, replace = false)
                    extend_with_dead_cells!(cells, Nbuff, make_dead_cell)
                    curr_P += 1
                    curr_t = curr_t_candidate
                    sync_passage_schedule!()
                    advance_passage_index!(curr_t)
                    next_t = compute_next_t()
                end
            end
        elseif live_count == 0
            _record_abm_outputs!(kmc_out, cells, rep, curr_P, cell_lin_df_vec,
                                 Nvec, nS_vec, nR_vec, nE_vec, tvec, Pvec,
                                 sub_sample_cells = sub_sample_cells, K = K,
                                 samp_cell_lin_df_vec = samp_cell_lin_df_vec,
                                 neutral_tracker = neutral_tracker)
            curr_t = min(max(kmc_last_t, curr_t), tmax)
            break
        elseif kmc_last_t >= (tmax - progress_tol)
            _record_abm_outputs!(kmc_out, cells, rep, curr_P, cell_lin_df_vec,
                                 Nvec, nS_vec, nR_vec, nE_vec, tvec, Pvec,
                                 sub_sample_cells = sub_sample_cells, K = K,
                                 samp_cell_lin_df_vec = samp_cell_lin_df_vec,
                                 neutral_tracker = neutral_tracker)
            curr_t = min(max(kmc_last_t, curr_t), tmax)
            break
        elseif tP_count < curr_P
            curr_t = curr_t_candidate
            sync_passage_schedule!()
            advance_passage_index!(curr_t)
            next_t = compute_next_t()
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
                                         samp_cell_lin_df_vec = samp_cell_lin_df_vec,
                                         neutral_tracker = neutral_tracker)
                    live_cells = alive_cells(cells)
                    cells = sample(live_cells, Nseed, replace = false)
                    extend_with_dead_cells!(cells, Nbuff, make_dead_cell)
                    curr_P += 1
                    curr_t = curr_t_candidate
                    sync_passage_schedule!()
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
    R_real::String = "l",
    drug_effect::Symbol = model.params.drug_effect,
    skew_lib::Bool = model.abm.skew_lib,
    use_lib_probs::Bool = model.abm.use_lib_probs,
    bc_unif::Float64 = model.abm.bc_unif,
    Nbc::Int64 = model.abm.Nbc,
    bc_probs::Vector{Float64} = model.abm.bc_probs,
    dt_save_at::Float64 = model.abm.dt_save_at,
    t_frac::Float64 = model.abm.t_frac,
    neutral_tracker = _neutral_tracker(model.abm))

    Nbuff = model.abm.Nbuff
    barcode_kwargs = (; skew_lib = skew_lib, use_lib_probs = use_lib_probs,
                       bc_unif = bc_unif, Nbc = Nbc, bc_probs = bc_probs)
    exp_cells = seed_cells(exp.n0, model.params.rho, Nbuff;
                           barcode_kwargs...)
    _seed_neutral_mutations!(exp_cells, model.abm)

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
            run_model_core_abm(expansion_model, ResPopABMState(exp_cells), stage_sim; treat = false, neutral_tracker = neutral_tracker)

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
        run_model_core_abm(expansion_model, ResPopABMState(exp_cells), final_sim; treat = false, neutral_tracker = neutral_tracker)
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
        run_model_core_abm(expansion_model, ResPopABMState(exp_cells), final_sim; treat = false, neutral_tracker = neutral_tracker)
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
    R_real = normalize_R_real(_kw(kwargs, :R_real, "l"))
    t_frac = _kw(kwargs, :t_frac, model.abm.t_frac)
    just_lin = _kw(kwargs, :just_lin, false)
    de = normalize_respop_drug_effect(_kw(kwargs, :drug_effect, model.params.drug_effect))
    drug_treatment = _kw(kwargs, :drug_treatment, exp.drug_treatment)
    sub_sample_cells = _kw(kwargs, :sub_sample_cells, model.abm.sub_sample_cells)
    K = _kw(kwargs, :K, model.abm.K)
    skew_lib = _kw(kwargs, :skew_lib, model.abm.skew_lib)
    bc_unif = _kw(kwargs, :bc_unif, model.abm.bc_unif)
    Nbc = _kw(kwargs, :Nbc, model.abm.Nbc)
    dt_save_at = _kw(kwargs, :dt_save_at, model.abm.dt_save_at)

    _validate_tmax_vector_constraints(exp.tmax, exp.t_Pass)
    _validate_tmax_length(exp.tmax, n_rep)
    _validate_instant_kill_run(model.abm, drug_treatment, exp.tmax isa AbstractVector ? exp.tmax : (exp.tmax,))

    n_pass_eff = exp.tmax isa AbstractVector ? 1 : (length(_passage_times(exp.t_Pass, Float64(exp.tmax))) + 1)

    validate_resistant_birth_drug_effect(model.params, R_real; drug_effect = de)
    model_eff = _with_drug_effect(model, de)
    neutral_tracker = _neutral_tracker(model.abm)
    rep_cells = _expand_split_cells_abm(model_eff, exp, n_rep;
                                        R_real = R_real,
                                        drug_effect = de,
                                        skew_lib = skew_lib,
                                        bc_unif = bc_unif,
                                        Nbc = Nbc,
                                        dt_save_at = dt_save_at,
                                        t_frac = t_frac,
                                        neutral_tracker = neutral_tracker)

    fin_t_outs = Float64[]
    fin_u_outs = Float64[]
    lin_df_outs = DataFrame[]
    sub_lin_df_outs = DataFrame[]
    sim_dfs = DataFrame[]
    kill_rows = model.abm.insta_kill ? NamedTuple[] : nothing

    nseed_last = _nseed_last(exp.Nseed)
    for i in 1:n_rep
        rep_tmax = _replicate_tmax(exp.tmax, n_rep, i)
        instant_kill = _instant_kill_state(model.abm)

        extend_with_dead_cells!(rep_cells[i], model.abm.Nbuff, make_dead_cell)

        sim = _run_abm_passage_experiment!(
            model_eff, rep_cells[i];
            t0 = 0.0, tmax = rep_tmax, t_Pass = exp.t_Pass,
            Nseed = nseed_last, Nmax = exp.Nmax, Cc = exp.Cc,
            treat_ons = exp.treat_ons, treat_offs = exp.treat_offs,
            dt_save_at = dt_save_at, Nbuff = model.abm.Nbuff,
            R_real = R_real, t_frac = t_frac, rep = i,
            treat = drug_treatment, drug_effect = de,
            sub_sample_cells = sub_sample_cells, K = K,
            neutral_tracker = neutral_tracker,
            instant_kill = instant_kill
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
        instant_kill === nothing || push!(kill_rows, _instant_kill_row(instant_kill, i))

        push!(lin_df_outs, join_dfs(sim["cell_lin_df_vec"], "bc"))
        if sub_sample_cells
            push!(sub_lin_df_outs, join_dfs(sim["sub_samp_cell_lin_df_vec"], "bc"))
        end

        if !just_lin
            t_outs = Vector{Float64}(undef, length(exp.t_keep) + n_pass_eff)
            u_outs = Vector{Float64}(undef, length(exp.t_keep) + n_pass_eff)

            for j in 1:length(exp.t_keep)
                idx = j
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
                idx = length(exp.t_keep) + j
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
    kill_rows === nothing || (out["inst_kill_df"] = DataFrame(kill_rows))
    return _add_neutral_outputs!(out, neutral_tracker)
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
    K::Int64 = 0,
    neutral_tracker::Union{Nothing, NeutralMutationTracker} = nothing
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

    function sync_passage_schedule!()
        tP_count = max(tP_count, min(curr_P, length(t_pass_vec) + 1))
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
    nDS_vec = Int64[]
    nDR_vec = Int64[]
    nR_vec = Int64[]
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
        kmc_out = run_model_core_abm(model, state, sim; treat = treat, neutral_tracker = neutral_tracker)
        kmc_last_t = last(kmc_out.tvec)
        isfinite(kmc_last_t) || error("ABM produced a non-finite time value.")

        curr_t_candidate = min(max(kmc_last_t, curr_t), next_t)

        live_count = length(alive_cells(cells))
        if last(kmc_out.Nvec) >= Nmax
            if last(kmc_out.Pvec) >= n_pass_eff
                _record_resdmg_abm_outputs!(kmc_out, cells, rep, curr_P, cell_lin_df_vec,
                                            Nvec, nS_vec, nDS_vec, nDR_vec, nR_vec, tvec, Pvec,
                                            sub_sample_cells = sub_sample_cells, K = K,
                                            samp_cell_lin_df_vec = samp_cell_lin_df_vec,
                                            neutral_tracker = neutral_tracker)
                curr_t = min(max(kmc_last_t, curr_t), tmax)
                break
            else
                _record_resdmg_abm_outputs!(kmc_out, cells, rep, curr_P, cell_lin_df_vec,
                                            Nvec, nS_vec, nDS_vec, nDR_vec, nR_vec, tvec, Pvec,
                                            sub_sample_cells = sub_sample_cells, K = K,
                                            samp_cell_lin_df_vec = samp_cell_lin_df_vec,
                                            neutral_tracker = neutral_tracker)
                if live_count < Nseed
                    break
                else
                    live_cells = alive_cells(cells)
                    cells = sample(live_cells, Nseed, replace = false)
                    extend_with_dead_cells!(cells, Nbuff, make_dead_resdmg_cell)
                    curr_P += 1
                    curr_t = curr_t_candidate
                    sync_passage_schedule!()
                    advance_passage_index!(curr_t)
                    next_t = compute_next_t()
                end
            end
        elseif live_count == 0
            _record_resdmg_abm_outputs!(kmc_out, cells, rep, curr_P, cell_lin_df_vec,
                                        Nvec, nS_vec, nDS_vec, nDR_vec, nR_vec, tvec, Pvec,
                                        sub_sample_cells = sub_sample_cells, K = K,
                                        samp_cell_lin_df_vec = samp_cell_lin_df_vec,
                                        neutral_tracker = neutral_tracker)
            curr_t = min(max(kmc_last_t, curr_t), tmax)
            break
        elseif kmc_last_t >= (tmax - progress_tol)
            _record_resdmg_abm_outputs!(kmc_out, cells, rep, curr_P, cell_lin_df_vec,
                                        Nvec, nS_vec, nDS_vec, nDR_vec, nR_vec, tvec, Pvec,
                                        sub_sample_cells = sub_sample_cells, K = K,
                                        samp_cell_lin_df_vec = samp_cell_lin_df_vec,
                                        neutral_tracker = neutral_tracker)
            curr_t = min(max(kmc_last_t, curr_t), tmax)
            break
        elseif tP_count < curr_P
            curr_t = curr_t_candidate
            sync_passage_schedule!()
            advance_passage_index!(curr_t)
            next_t = compute_next_t()
        elseif tP_count <= length(t_pass_vec) && kmc_last_t >= (t_pass_vec[tP_count] - progress_tol)
            if last(kmc_out.Pvec) < n_pass_eff
                if live_count < Nseed
                    update_track_vec_resdmg!(kmc_out, Nvec, nS_vec, nDS_vec, nDR_vec, nR_vec, tvec, Pvec)
                    curr_t = curr_t_candidate
                    advance_passage_index!(curr_t)
                    next_t = compute_next_t()
                else
                    _record_resdmg_abm_outputs!(kmc_out, cells, rep, curr_P, cell_lin_df_vec,
                                                Nvec, nS_vec, nDS_vec, nDR_vec, nR_vec, tvec, Pvec,
                                                sub_sample_cells = sub_sample_cells, K = K,
                                                samp_cell_lin_df_vec = samp_cell_lin_df_vec,
                                                neutral_tracker = neutral_tracker)
                    live_cells = alive_cells(cells)
                    cells = sample(live_cells, Nseed, replace = false)
                    extend_with_dead_cells!(cells, Nbuff, make_dead_resdmg_cell)
                    curr_P += 1
                    curr_t = curr_t_candidate
                    sync_passage_schedule!()
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
        "nDS_vec" => nDS_vec,
        "nDR_vec" => nDR_vec,
        "nR_vec" => nR_vec,
    )
    if sub_sample_cells
        out["sub_samp_cell_lin_df_vec"] = samp_cell_lin_df_vec
    end
    return out
end

function _expand_split_cells_abm(model::ResDmg_ABM, exp::ExperimentParams, n_rep::Int64;
    R_real::String = "l",
    drug_effect::Symbol = model.params.drug_effect,
    skew_lib::Bool = model.abm.skew_lib,
    use_lib_probs::Bool = model.abm.use_lib_probs,
    bc_unif::Float64 = model.abm.bc_unif,
    Nbc::Int64 = model.abm.Nbc,
    bc_probs::Vector{Float64} = model.abm.bc_probs,
    dt_save_at::Float64 = model.abm.dt_save_at,
    t_frac::Float64 = model.abm.t_frac,
    neutral_tracker = _neutral_tracker(model.abm))

    Nbuff = model.abm.Nbuff
    barcode_kwargs = (; skew_lib = skew_lib, use_lib_probs = use_lib_probs,
                       bc_unif = bc_unif, Nbc = Nbc, bc_probs = bc_probs)
    exp_cells = seed_resdmg_cells(exp.n0, model.params.rho, Nbuff;
                                  barcode_kwargs...)
    _seed_neutral_mutations!(exp_cells, model.abm)

    expansion_model = ResDmg_ABM(_copy_resdmg_params(model.params; drug_effect = drug_effect);
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
            run_model_core_abm(expansion_model, ResDmgABMState(exp_cells), stage_sim; treat = false, neutral_tracker = neutral_tracker)

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
        run_model_core_abm(expansion_model, ResDmgABMState(exp_cells), final_sim; treat = false, neutral_tracker = neutral_tracker)
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
        run_model_core_abm(expansion_model, ResDmgABMState(exp_cells), final_sim; treat = false, neutral_tracker = neutral_tracker)
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
    R_real = normalize_R_real(_kw(kwargs, :R_real, "l"))
    t_frac = _kw(kwargs, :t_frac, model.abm.t_frac)
    just_lin = _kw(kwargs, :just_lin, false)
    de = normalize_resdmg_drug_effect(_kw(kwargs, :drug_effect, model.params.drug_effect))
    drug_treatment = _kw(kwargs, :drug_treatment, exp.drug_treatment)
    sub_sample_cells = _kw(kwargs, :sub_sample_cells, model.abm.sub_sample_cells)
    K = _kw(kwargs, :K, model.abm.K)
    skew_lib = _kw(kwargs, :skew_lib, model.abm.skew_lib)
    bc_unif = _kw(kwargs, :bc_unif, model.abm.bc_unif)
    Nbc = _kw(kwargs, :Nbc, model.abm.Nbc)
    dt_save_at = _kw(kwargs, :dt_save_at, model.abm.dt_save_at)

    _validate_tmax_vector_constraints(exp.tmax, exp.t_Pass)
    _validate_tmax_length(exp.tmax, n_rep)

    n_pass_eff = exp.tmax isa AbstractVector ? 1 : (length(_passage_times(exp.t_Pass, Float64(exp.tmax))) + 1)

    validate_resistant_birth_drug_effect(model.params, R_real; drug_effect = de)
    model_eff = _with_drug_effect(model, de)
    neutral_tracker = _neutral_tracker(model.abm)
    rep_cells = _expand_split_cells_abm(model_eff, exp, n_rep;
                                        R_real = R_real,
                                        drug_effect = de,
                                        skew_lib = skew_lib,
                                        bc_unif = bc_unif,
                                        Nbc = Nbc,
                                        dt_save_at = dt_save_at,
                                        t_frac = t_frac,
                                        neutral_tracker = neutral_tracker)

    fin_t_outs = Float64[]
    fin_u_outs = Float64[]
    lin_df_outs = DataFrame[]
    sub_lin_df_outs = DataFrame[]
    sim_dfs = DataFrame[]

    nseed_last = _nseed_last(exp.Nseed)
    for i in 1:n_rep
        rep_tmax = _replicate_tmax(exp.tmax, n_rep, i)

        extend_with_dead_cells!(rep_cells[i], model.abm.Nbuff, make_dead_resdmg_cell)

        sim = _run_abm_passage_experiment!(
            model_eff, rep_cells[i];
            t0 = 0.0, tmax = rep_tmax, t_Pass = exp.t_Pass,
            Nseed = nseed_last, Nmax = exp.Nmax, Cc = exp.Cc,
            treat_ons = exp.treat_ons, treat_offs = exp.treat_offs,
            dt_save_at = dt_save_at, Nbuff = model.abm.Nbuff,
            R_real = R_real, t_frac = t_frac, rep = i,
            treat = drug_treatment, drug_effect = de,
            sub_sample_cells = sub_sample_cells, K = K,
            neutral_tracker = neutral_tracker
        )

        sim_df = DataFrame(
            t = sim["tvec"],
            N = sim["Nvec"],
            nS = sim["nS_vec"],
            nDS = sim["nDS_vec"],
            nDR = sim["nDR_vec"],
            nR = sim["nR_vec"],
            rep = i
        )
        push!(sim_dfs, sim_df)

        push!(lin_df_outs, join_dfs(sim["cell_lin_df_vec"], "bc"))
        if sub_sample_cells
            push!(sub_lin_df_outs, join_dfs(sim["sub_samp_cell_lin_df_vec"], "bc"))
        end

        if !just_lin
            t_outs = Vector{Float64}(undef, length(exp.t_keep) + n_pass_eff)
            u_outs = Vector{Float64}(undef, length(exp.t_keep) + n_pass_eff)

            for j in 1:length(exp.t_keep)
                idx = j
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
                idx = length(exp.t_keep) + j
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
    return _add_neutral_outputs!(out, neutral_tracker)
end

function _run_abm_simple!(
    model::ResPop_ABM,
    cells::Vector{CancerCell};
    t0::Float64,
    tmax::Float64,
    Nmax::Int64,
    Cc::Int64,
    treat_ons::Vector{Float64},
    treat_offs::Vector{Float64},
    dt_save_at::Float64,
    R_real::String,
    t_frac::Float64,
    rep::Int64,
    treat::Bool,
    drug_effect::Symbol,
    sub_sample_cells::Bool = false,
    K::Int64 = 0,
    neutral_tracker::Union{Nothing, NeutralMutationTracker} = nothing,
    instant_kill::Union{Nothing, InstantKillState} = nothing
)
    sim = ABMSimParams(
        t0 = t0,
        tmax = tmax,
        Nmax = Nmax,
        Cc = Cc,
        treat_ons = treat_ons,
        treat_offs = treat_offs,
        dt_save_at = dt_save_at,
        R_real = R_real,
        t_frac = t_frac,
        Passage = 1,
        drug_effect = drug_effect
    )
    state = ResPopABMState(cells)
    kmc_out = run_model_core_abm(model, state, sim; treat = treat, neutral_tracker = neutral_tracker, instant_kill = instant_kill)

    cell_lin_df_vec = DataFrame[]
    samp_cell_lin_df_vec = DataFrame[]
    Nvec = Int64[]
    nS_vec = Int64[]
    nR_vec = Int64[]
    nE_vec = Int64[]
    tvec = Float64[]
    Pvec = Int64[]

    _record_abm_outputs!(kmc_out, cells, rep, 1, cell_lin_df_vec,
                         Nvec, nS_vec, nR_vec, nE_vec, tvec, Pvec,
                         sub_sample_cells = sub_sample_cells, K = K,
                         samp_cell_lin_df_vec = samp_cell_lin_df_vec,
                         neutral_tracker = neutral_tracker)

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

function _run_abm_simple!(
    model::ResDmg_ABM,
    cells::Vector{ResDmgCell};
    t0::Float64,
    tmax::Float64,
    Nmax::Int64,
    Cc::Int64,
    treat_ons::Vector{Float64},
    treat_offs::Vector{Float64},
    dt_save_at::Float64,
    R_real::String,
    t_frac::Float64,
    rep::Int64,
    treat::Bool,
    drug_effect::Symbol,
    sub_sample_cells::Bool = false,
    K::Int64 = 0,
    neutral_tracker::Union{Nothing, NeutralMutationTracker} = nothing
)
    sim = ABMSimParams(
        t0 = t0,
        tmax = tmax,
        Nmax = Nmax,
        Cc = Cc,
        treat_ons = treat_ons,
        treat_offs = treat_offs,
        dt_save_at = dt_save_at,
        R_real = R_real,
        t_frac = t_frac,
        Passage = 1,
        drug_effect = drug_effect
    )
    state = ResDmgABMState(cells)
    kmc_out = run_model_core_abm(model, state, sim; treat = treat, neutral_tracker = neutral_tracker)

    cell_lin_df_vec = DataFrame[]
    samp_cell_lin_df_vec = DataFrame[]
    Nvec = Int64[]
    nS_vec = Int64[]
    nDS_vec = Int64[]
    nDR_vec = Int64[]
    nR_vec = Int64[]
    tvec = Float64[]
    Pvec = Int64[]

    _record_resdmg_abm_outputs!(kmc_out, cells, rep, 1, cell_lin_df_vec,
                                Nvec, nS_vec, nDS_vec, nDR_vec, nR_vec, tvec, Pvec,
                                sub_sample_cells = sub_sample_cells, K = K,
                                samp_cell_lin_df_vec = samp_cell_lin_df_vec,
                                neutral_tracker = neutral_tracker)

    out = Dict(
        "cell_lin_df_vec" => cell_lin_df_vec,
        "Nvec" => Nvec,
        "tvec" => tvec,
        "Pvec" => Pvec,
        "nS_vec" => nS_vec,
        "nDS_vec" => nDS_vec,
        "nDR_vec" => nDR_vec,
        "nR_vec" => nR_vec,
    )
    if sub_sample_cells
        out["sub_samp_cell_lin_df_vec"] = samp_cell_lin_df_vec
    end
    return out
end

function _simulate_simple_abm(model::ResPop_ABM, sim::SimpleSimParams; kwargs...)
    R_real = normalize_R_real(_kw(kwargs, :R_real, "l"))
    t_frac = _kw(kwargs, :t_frac, model.abm.t_frac)
    de = normalize_respop_drug_effect(_kw(kwargs, :drug_effect, model.params.drug_effect))
    drug_treatment = _kw(kwargs, :drug_treatment, sim.drug_treatment)
    skew_lib = _kw(kwargs, :skew_lib, model.abm.skew_lib)
    use_lib_probs = _kw(kwargs, :use_lib_probs, model.abm.use_lib_probs)
    bc_unif = _kw(kwargs, :bc_unif, model.abm.bc_unif)
    Nbc = _kw(kwargs, :Nbc, model.abm.Nbc)
    bc_probs = _kw(kwargs, :bc_probs, model.abm.bc_probs)
    dt_save_at = _kw(kwargs, :dt_save_at, model.abm.dt_save_at)

    validate_resistant_birth_drug_effect(model.params, R_real; drug_effect = de)
    _validate_instant_kill_run(model.abm, drug_treatment, (sim.tmax,))
    model_eff = _with_drug_effect(model, de)
    instant_kill = _instant_kill_state(model.abm)
    neutral_tracker = _neutral_tracker(model.abm)
    barcode_kwargs = (; skew_lib = skew_lib, use_lib_probs = use_lib_probs,
                       bc_unif = bc_unif, Nbc = Nbc, bc_probs = bc_probs)
    cells = seed_cells(sim.n0, model.params.rho, model.abm.Nbuff;
                       barcode_kwargs...)
    _seed_neutral_mutations!(cells, model.abm)

    sim_out = _run_abm_simple!(
        model_eff, cells;
        t0 = 0.0, tmax = sim.tmax,
        Nmax = sim.Nmax, Cc = sim.Cc,
        treat_ons = sim.treat_ons, treat_offs = sim.treat_offs,
        dt_save_at = dt_save_at,
        R_real = R_real, t_frac = t_frac, rep = 1,
        treat = drug_treatment, drug_effect = de,
        sub_sample_cells = false, K = 0,
        neutral_tracker = neutral_tracker,
        instant_kill = instant_kill
    )

    sol_df = DataFrame(
        t = sim_out["tvec"],
        N = sim_out["Nvec"],
        nS = sim_out["nS_vec"],
        nR = sim_out["nR_vec"],
        nE = sim_out["nE_vec"]
    )

    out = Dict(
        "lin_df" => join_dfs(sim_out["cell_lin_df_vec"], "bc"),
        "sol_df" => sol_df
    )
    instant_kill === nothing || (out["inst_kill_df"] = DataFrame([_instant_kill_row(instant_kill, 1)]))
    return _add_neutral_outputs!(out, neutral_tracker)
end

function _simulate_simple_abm(model::ResDmg_ABM, sim::SimpleSimParams; kwargs...)
    R_real = normalize_R_real(_kw(kwargs, :R_real, "l"))
    t_frac = _kw(kwargs, :t_frac, model.abm.t_frac)
    de = normalize_resdmg_drug_effect(_kw(kwargs, :drug_effect, model.params.drug_effect))
    drug_treatment = _kw(kwargs, :drug_treatment, sim.drug_treatment)
    skew_lib = _kw(kwargs, :skew_lib, model.abm.skew_lib)
    use_lib_probs = _kw(kwargs, :use_lib_probs, model.abm.use_lib_probs)
    bc_unif = _kw(kwargs, :bc_unif, model.abm.bc_unif)
    Nbc = _kw(kwargs, :Nbc, model.abm.Nbc)
    bc_probs = _kw(kwargs, :bc_probs, model.abm.bc_probs)
    dt_save_at = _kw(kwargs, :dt_save_at, model.abm.dt_save_at)

    validate_resistant_birth_drug_effect(model.params, R_real; drug_effect = de)
    model_eff = _with_drug_effect(model, de)
    neutral_tracker = _neutral_tracker(model.abm)
    barcode_kwargs = (; skew_lib = skew_lib, use_lib_probs = use_lib_probs,
                       bc_unif = bc_unif, Nbc = Nbc, bc_probs = bc_probs)
    cells = seed_resdmg_cells(sim.n0, model.params.rho, model.abm.Nbuff;
                              barcode_kwargs...)
    _seed_neutral_mutations!(cells, model.abm)

    sim_out = _run_abm_simple!(
        model_eff, cells;
        t0 = 0.0, tmax = sim.tmax,
        Nmax = sim.Nmax, Cc = sim.Cc,
        treat_ons = sim.treat_ons, treat_offs = sim.treat_offs,
        dt_save_at = dt_save_at,
        R_real = R_real, t_frac = t_frac, rep = 1,
        treat = drug_treatment, drug_effect = de,
        sub_sample_cells = false, K = 0,
        neutral_tracker = neutral_tracker
    )

    sol_df = DataFrame(
        t = sim_out["tvec"],
        N = sim_out["Nvec"],
        nS = sim_out["nS_vec"],
        nDS = sim_out["nDS_vec"],
        nDR = sim_out["nDR_vec"],
        nR = sim_out["nR_vec"]
    )

    out = Dict(
        "lin_df" => join_dfs(sim_out["cell_lin_df_vec"], "bc"),
        "sol_df" => sol_df
    )
    return _add_neutral_outputs!(out, neutral_tracker)
end

function _expand_split_cells_abm(model::ResPopInVivo_ABM, exp::ExperimentParams, n_rep::Int64;
    R_real = "l", drug_effect = model.params.drug_effect, rep_design = _experiment_condition_design(true, false, n_rep),
    skew_lib = model.abm.skew_lib, use_lib_probs = model.abm.use_lib_probs,
    split_after_barcoding = model.abm.split_after_barcoding, bc_unif = model.abm.bc_unif,
    Nbc = model.abm.Nbc, bc_probs = model.abm.bc_probs, dt_save_at = model.abm.dt_save_at,
    t_frac = model.abm.t_frac, inc_pot = false, full_pheno_bc = false, pot_outputs = Dict{String, Any}())
    validate_abm_params(; skew_lib = Bool(skew_lib), use_lib_probs = Bool(use_lib_probs),
        Nbc = Int64(Nbc), bc_probs = Vector{Float64}(bc_probs))
    split_after_barcoding && !use_lib_probs && error("split_after_barcoding requires use_lib_probs=true.")
    nseed = _nseed_last(exp.Nseed)
    barcode_kwargs = (; skew_lib = Bool(skew_lib), use_lib_probs = Bool(use_lib_probs),
        bc_unif = Float64(bc_unif), Nbc = Int64(Nbc), bc_probs = bc_probs)
    n_buff = model.abm.Nbuff
    pot_dfs = DataFrame[]
    if split_after_barcoding
        cells = seed_invivo_cells(exp.n0, 0.0, 0.0, n_buff; barcode_kwargs...)
        _assign_invivo_phenotype_by_barcode!(cells, :R, model.params.rho)
        _assign_invivo_phenotype_by_barcode!(cells, :EG, model.params.fEG1)
        if inc_pot
            sim = ABMSimParams(tmax = 0.0, Nmax = exp.Nmax, Cc = exp.Cc,
                treat_ons = Float64[], treat_offs = Float64[], Passage = 0, R_real = R_real)
            out = run_model_core_abm(model, ResPopInVivoABMState(cells), sim)
            push!(pot_dfs, _invivo_abm_sol_df(out; cond = "POT", rep = Int64(0)))
        end
    else
        cells = seed_invivo_cells(exp.n0, model.params.rho, model.params.fEG1, n_buff; barcode_kwargs...)
        expansion = ResPopInVivo_ABM(_copy_respop_invivo_params(model.params; al = 0.0, drug_effect); abm = model.abm)
        if exp.t_exp isa AbstractVector && exp.Nseed isa AbstractVector
            length(exp.t_exp) == length(exp.Nseed) && !isempty(exp.t_exp) ||
                error("t_exp and Nseed vectors must be nonempty and have matching lengths.")
            durations, seeds = exp.t_exp, exp.Nseed
        elseif exp.t_exp isa Real && exp.Nseed isa Integer
            durations, seeds = [exp.t_exp], [exp.Nseed]
        else
            error("t_exp and Nseed must both be scalars or both be vectors of matching lengths.")
        end
        offset = 0.0
        for i in eachindex(durations)
            durations[i] >= 0 || error("Expansion durations must be nonnegative.")
            sim = ABMSimParams(tmax = Float64(durations[i]), Nmax = exp.Nmax, Cc = exp.Cc,
                treat_ons = Float64[], treat_offs = Float64[], R_real = R_real,
                dt_save_at = dt_save_at, t_frac = t_frac, Passage = 0, drug_effect = drug_effect)
            out = run_model_core_abm(expansion, ResPopInVivoABMState(cells), sim)
            if inc_pot
                df = _invivo_abm_sol_df(out; cond = "POT", rep = Int64(0))
                df.t .+= offset
                push!(pot_dfs, df)
            end
            offset+= out.fin_t
            if i != lastindex(durations)
                live = alive_cells(cells)
                length(live) >= seeds[i] || error("Not enough cells after expansion at stage $i for bottlenecking.")
                cells = sample(live, Int(seeds[i]), replace = false)
                extend_with_dead_cells!(cells, n_buff, make_dead_cell_invivo)
            end
        end
    end
    if inc_pot
        pot_outputs["sol_df"] = vcat(pot_dfs...)
        pot_outputs["lin_df"] = get_counts(alive_cells(cells), "POT_P0")
        full_pheno_bc && (pot_outputs["pheno_bc_df"] = get_pheno_counts(cells, "POT_P0"))
    end
    live = alive_cells(cells)
    n_batches = length(rep_design)
    length(live) >= n_batches*nseed || error("Not enough cells for $n_batches replicate batches of size $nseed.")
    sampled = sample(live, n_batches*nseed, replace = false)
    reps = Vector{Vector{InVivoCancerCell}}()
    grafts = NamedTuple[]
    for (i, rec) in enumerate(rep_design)
        batch = sampled[((i-1)*nseed+1):(i*nseed)]
        batch, stats = engraftment_selection(batch, model.params.pEG, model.params.sEG)
        _push_invivo_engraftment!(grafts, rec.cond, rec.rep, 1, stats)
        extend_with_dead_cells!(batch, n_buff, make_dead_cell_invivo)
        push!(reps, batch)
    end
    return reps, grafts
end

function _run_abm_passage_experiment!(model::ResPopInVivo_ABM, cells::Vector{InVivoCancerCell};
    tmax::Float64, t_Pass, Nseed::Int64, Nmax::Int64, Cc::Int64,
    treat_ons, treat_offs, dt_save_at, Nbuff::Int64, R_real::String, t_frac,
    rep::Int64, treat::Bool, drug_effect::Symbol, cond::String = "DT",
    sub_sample_cells::Bool = false, K::Int64 = 0, full_pheno_bc::Bool = false)
    schedule = _passage_times(t_Pass, tmax)
    n_pass = length(schedule)+1
    schedule_idx = 1
    passage = 1
    time = 0.0
    dfs, lin_dfs, pheno_dfs, sub_dfs = DataFrame[], DataFrame[], DataFrame[], DataFrame[]
    grafts = NamedTuple[]
    while true
        while schedule_idx <= length(schedule) && schedule[schedule_idx] <= time
            schedule_idx+= 1
        end
        boundary = schedule_idx <= length(schedule) ? min(tmax, schedule[schedule_idx]) : tmax
        sim = ABMSimParams(t0 = time, tmax = boundary, Nmax = Nmax, Cc = Cc, treat_ons = treat_ons,
            treat_offs = treat_offs, dt_save_at = dt_save_at, R_real = R_real, t_frac = t_frac,
            Passage = passage, drug_effect = drug_effect)
        out = run_model_core_abm(model, ResPopInVivoABMState(cells), sim; treat = treat)
        push!(dfs, _invivo_abm_sol_df(out; cond = cond, rep = rep))
        time = last(out.tvec)
        live = alive_cells(cells)
        at_cap = length(live) >= Nmax
        finished = isempty(live) || time >= tmax || (at_cap && passage >= n_pass)
        scheduled = !at_cap && time >= boundary && boundary<tmax
        passage_now = !finished && passage<n_pass && (at_cap || (scheduled && length(live) >= Nseed))
        if finished || passage_now
            _record_invivo_abm_readout!(cells, string(cond, rep, "_P", passage), lin_dfs, pheno_dfs, sub_dfs;
                full_pheno_bc, sub_sample_cells, K)
        end
        finished && break
        if passage_now
            cells = sample(live, Nseed, replace = false)
            cells, stats = engraftment_selection(cells, model.params.pEG, model.params.sEG)
            passage+= 1
            _push_invivo_engraftment!(grafts, cond, rep, passage, stats)
            extend_with_dead_cells!(cells, Nbuff, make_dead_cell_invivo)
            schedule_idx = max(schedule_idx, passage)
        elseif scheduled
            schedule_idx+= 1
        else
            error("InVivo ABM segment ended without a stopping or passage event.")
        end
    end
    result = Dict{String, Any}("sol_df" => vcat(dfs...), "lin_df" => join_dfs(lin_dfs, "bc"),
        "t" => time, "u" => Float64(n_alive(cells)), "engraft_rows" => grafts)
    full_pheno_bc && (result["pheno_bc_df"] = join_dfs(pheno_dfs, "bc"))
    sub_sample_cells && (result["sub_lin_df"] = join_dfs(sub_dfs, "bc"))
    return result
end

function _simulate_experiment_abm(model::ResPopInVivo_ABM, exp::ExperimentParams; kwargs...)
    n_rep = Int64(_kw(kwargs, :n_rep, exp.n_rep))
    n_rep>0 || error("n_rep must be > 0.")
    R_real = normalize_R_real(_kw(kwargs, :R_real, "l"))
    de = normalize_respop_drug_effect(_kw(kwargs, :drug_effect, model.params.drug_effect))
    validate_resistant_birth_drug_effect(model.params, R_real; drug_effect = de)
    model_eff = _with_drug_effect(model, de)
    _validate_tmax_vector_constraints(exp.tmax, exp.t_Pass)
    _validate_tmax_length(exp.tmax, n_rep)
    treat = Bool(_kw(kwargs, :drug_treatment, exp.drug_treatment))
    inc_control = Bool(_kw(kwargs, :inc_control, exp.inc_control))
    inc_pot = Bool(_kw(kwargs, :inc_pot, exp.inc_pot))
    just_lin = Bool(_kw(kwargs, :just_lin, false))
    full_pheno_bc = Bool(_kw(kwargs, :full_pheno_bc, model.abm.full_pheno_bc))
    sub_sample_cells = Bool(_kw(kwargs, :sub_sample_cells, model.abm.sub_sample_cells))
    K = Int64(_kw(kwargs, :K, model.abm.K))
    dt_save_at = _kw(kwargs, :dt_save_at, model.abm.dt_save_at)
    t_frac = _kw(kwargs, :t_frac, model.abm.t_frac)
    design = _experiment_condition_design(treat, inc_control, n_rep)
    pot = Dict{String, Any}()
    reps, grafts = _expand_split_cells_abm(model_eff, exp, n_rep; R_real = R_real, drug_effect = de,
        rep_design = design, inc_pot = inc_pot, full_pheno_bc = full_pheno_bc, pot_outputs = pot,
        skew_lib = _kw(kwargs, :skew_lib, model.abm.skew_lib),
        use_lib_probs = _kw(kwargs, :use_lib_probs, model.abm.use_lib_probs),
        split_after_barcoding = _kw(kwargs, :split_after_barcoding, model.abm.split_after_barcoding),
        bc_unif = _kw(kwargs, :bc_unif, model.abm.bc_unif), Nbc = _kw(kwargs, :Nbc, model.abm.Nbc),
        bc_probs = _kw(kwargs, :bc_probs, model.abm.bc_probs), dt_save_at = dt_save_at, t_frac = t_frac)
    dfs, lin_dfs, pheno_dfs, sub_dfs = DataFrame[], DataFrame[], DataFrame[], DataFrame[]
    if inc_pot
        !just_lin && push!(dfs, pot["sol_df"])
        push!(lin_dfs, pot["lin_df"])
        full_pheno_bc && push!(pheno_dfs, pot["pheno_bc_df"])
    end
    times, totals, conditions, replicates = Float64[], Float64[], String[], Int64[]
    for (i, rec) in enumerate(design)
        result = _run_abm_passage_experiment!(model_eff, reps[i];
            tmax = _replicate_tmax(exp.tmax, n_rep, rec.rep), t_Pass = exp.t_Pass,
            Nseed = _nseed_last(exp.Nseed), Nmax = exp.Nmax, Cc = exp.Cc,
            treat_ons = exp.treat_ons, treat_offs = exp.treat_offs, dt_save_at = dt_save_at,
            Nbuff = model.abm.Nbuff, R_real = R_real, t_frac = t_frac, rep = rec.rep, treat = rec.treat,
            drug_effect = de, cond = rec.cond, sub_sample_cells = sub_sample_cells, K = K, full_pheno_bc = full_pheno_bc)
        push!(lin_dfs, result["lin_df"])
        full_pheno_bc && push!(pheno_dfs, result["pheno_bc_df"])
        sub_sample_cells && push!(sub_dfs, result["sub_lin_df"])
        append!(grafts, result["engraft_rows"])
        if !just_lin
            push!(dfs, result["sol_df"])
            push!(times, result["t"]); push!(totals, result["u"])
            push!(conditions, rec.cond); push!(replicates, rec.rep)
        end
    end
    result = Dict{String, Any}("lin_df" => join_dfs(lin_dfs, "bc"), "engraft_df" => _invivo_engraftment_df(grafts))
    full_pheno_bc && (result["pheno_bc_df"] = join_dfs(pheno_dfs, "bc"))
    sub_sample_cells && (result["sub_lin_df"] = join_dfs(sub_dfs, "bc"))
    if !just_lin
        merge!(result, Dict("sol_df" => vcat(dfs...), "t" => times, "u" => totals, "cond" => conditions, "rep" => replicates))
    end
    return result
end

function _simulate_simple_abm(model::ResPopInVivo_ABM, settings::SimpleSimParams; kwargs...)
    R_real = normalize_R_real(_kw(kwargs, :R_real, "l"))
    de = normalize_respop_drug_effect(_kw(kwargs, :drug_effect, model.params.drug_effect))
    model_eff = _with_drug_effect(model, de)
    validate_resistant_birth_drug_effect(model.params, R_real; drug_effect = de)
    cells = seed_invivo_cells(settings.n0, model.params.rho, model.params.fEG1, model.abm.Nbuff;
        skew_lib = Bool(_kw(kwargs, :skew_lib, model.abm.skew_lib)),
        use_lib_probs = Bool(_kw(kwargs, :use_lib_probs, model.abm.use_lib_probs)),
        bc_unif = Float64(_kw(kwargs, :bc_unif, model.abm.bc_unif)),
        Nbc = Int64(_kw(kwargs, :Nbc, model.abm.Nbc)), bc_probs = _kw(kwargs, :bc_probs, model.abm.bc_probs))
    sim = ABMSimParams(tmax = settings.tmax, Nmax = settings.Nmax, Cc = settings.Cc,
        treat_ons = settings.treat_ons, treat_offs = settings.treat_offs, R_real = R_real,
        dt_save_at = _kw(kwargs, :dt_save_at, model.abm.dt_save_at),
        t_frac = _kw(kwargs, :t_frac, model.abm.t_frac), drug_effect = de)
    out = run_model_core_abm(model_eff, ResPopInVivoABMState(cells), sim;
        treat = Bool(_kw(kwargs, :drug_treatment, settings.drug_treatment)))
    result = Dict{String, Any}("lin_df" => get_counts(alive_cells(cells), "DT1_P1"))
    Bool(_kw(kwargs, :full_pheno_bc, model.abm.full_pheno_bc)) &&
        (result["pheno_bc_df"] = get_pheno_counts(cells, "DT1_P1"))
    if Bool(_kw(kwargs, :sub_sample_cells, model.abm.sub_sample_cells))
        live = alive_cells(cells)
        K = Int64(_kw(kwargs, :K, model.abm.K))
        sampled = K <= length(live) ? sample(live, K, replace = false) : live
        result["sub_lin_df"] = get_counts(sampled, "DT1_P1")
    end
    if !Bool(_kw(kwargs, :just_lin, false))
        df = _invivo_abm_sol_df(out; cond = "DT", rep = Int64(1))
        select!(df, Not([:cond, :rep]))
        result["sol_df"] = df
    end
    return result
end
